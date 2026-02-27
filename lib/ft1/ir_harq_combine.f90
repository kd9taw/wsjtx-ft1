! FT1 - 4-CPM turbo equalization mode for WSJT-X
! Copyright (C) 2026 Seth McCall, KD9TAW
!
! This file is part of WSJT-X.
!
! WSJT-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! WSJT-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with WSJT-X. If not, see <https://www.gnu.org/licenses/>.
!
module ir_harq_combine_mod

! IR-HARQ soft LLR combining for FT1.
!
! Manages per-signal LLR buffers across T/R periods.
! When a decode fails on RV0, the turbo-extracted LLRs are stored.
! When RV1 or RV2 arrives at the same frequency, LLRs are combined
! and decoded with the extended LDPC code.
!
! Buffer key: (frequency_bin, T/R parity) with ±10 Hz matching tolerance.
! Expiry: 30 seconds (7-8 T/R periods at 4.0s).
!
! Reference: specs/3d_protocol.md Section 8.7
!            ir_harq_protocol_design.md Section 3.3-3.6

  use ldpc348_91_mod
  implicit none

  integer, parameter :: MAX_HARQ_SLOTS = 100    !Max simultaneous in-progress decodes
  integer, parameter :: FREQ_TOL_HZ = 10        !Frequency matching tolerance (Hz)
  integer, parameter :: EXPIRY_MS = 30000        !Buffer expiry (30 seconds)

  type :: harq_slot
     logical :: active = .false.
     real :: freq = 0.0                          !Signal frequency (Hz)
     integer :: rv_count = -1                    !Highest RV stored (0, 1, or 2)
     integer :: timestamp_ms = 0                 !Timestamp of last update
     real :: llr_rv0(174)                        !LLRs from RV0 turbo decode
     real :: llr_rv1_new(87)                     !New parity LLRs from RV1
     real :: llr_rv1_repeat(87)                  !Repeated systematic LLRs from RV1
     real :: llr_rv2_new(87)                     !New parity LLRs from RV2
     real :: llr_rv2_repeat(87)                  !Repeated systematic LLRs from RV2
  end type harq_slot

  type(harq_slot) :: slots(MAX_HARQ_SLOTS)
  logical :: harq_initialized = .false.

contains

  subroutine harq_init()
    implicit none
    integer :: i
    do i = 1, MAX_HARQ_SLOTS
       slots(i)%active = .false.
       slots(i)%rv_count = -1
    enddo
    harq_initialized = .true.
  end subroutine harq_init


  subroutine harq_store_rv0(freq, llr174, timestamp_ms)
  ! Store LLRs from a failed RV0 decode attempt.
  ! Called after turbo_decode_ft1 returns ntype < 0.

    implicit none
    real, intent(in) :: freq
    real, intent(in) :: llr174(174)
    integer, intent(in) :: timestamp_ms
    integer :: islot

    if(.not.harq_initialized) call harq_init()

    ! Find existing slot for this frequency or allocate new one
    islot = find_slot(freq, timestamp_ms)
    if(islot .le. 0) then
       islot = allocate_slot(timestamp_ms)
       if(islot .le. 0) return          !Buffer full
    endif

    slots(islot)%active = .true.
    slots(islot)%freq = freq
    slots(islot)%rv_count = 0
    slots(islot)%timestamp_ms = timestamp_ms
    slots(islot)%llr_rv0 = llr174

  end subroutine harq_store_rv0


  subroutine harq_combine_rv1(freq, llr174_rv1, timestamp_ms, &
       message77, nharderror, decode_ok)
  ! Combine RV0 + RV1 LLRs and attempt decode with LDPC(261,91).
  !
  ! RV1 transmitted bits: 87 new parity + 87 repeated systematic
  ! Combined vector (261 elements):
  !   Bits 1-174:   llr_rv0 + chase component from RV1 repeated bits
  !   Bits 175-261: new parity LLRs from RV1

    implicit none
    real, intent(in) :: freq
    real, intent(in) :: llr174_rv1(174)
    integer, intent(in) :: timestamp_ms
    integer*1, intent(out) :: message77(77)
    integer, intent(out) :: nharderror
    logical, intent(out) :: decode_ok

    real :: llr_combined(N_EXT1)         !261 combined LLRs
    integer :: islot, iter, i

    decode_ok = .false.
    nharderror = -1
    message77 = 0

    if(.not.harq_initialized) call harq_init()

    islot = find_slot(freq, timestamp_ms)
    if(islot .le. 0) return             !No stored RV0 for this frequency
    if(slots(islot)%rv_count .lt. 0) return

    ! Combine LLRs
    ! First 174 bits: add RV0 LLRs + chase combining from repeated systematic
    llr_combined(1:174) = slots(islot)%llr_rv0(1:174)

    ! RV1 transmitted: bits 1-87 = new parity, bits 88-174 = systematic repeat
    ! Add repeated systematic LLRs to positions 1-87 of combined vector
    do i = 1, 87
       llr_combined(i) = llr_combined(i) + llr174_rv1(87 + i)
    enddo

    ! New parity LLRs go to positions 175-261
    llr_combined(175:261) = llr174_rv1(1:87)

    ! Store RV1 components for potential RV2 combining
    slots(islot)%llr_rv1_new = llr174_rv1(1:87)
    slots(islot)%llr_rv1_repeat = llr174_rv1(88:174)
    slots(islot)%rv_count = 1
    slots(islot)%timestamp_ms = timestamp_ms

    ! Decode with extended LDPC(261,91)
    call bpdecode_ext(llr_combined, N_EXT1, 50, message77, nharderror, iter)

    if(nharderror .ge. 0) then
       decode_ok = .true.
       call harq_clear_slot(islot)
    endif

  end subroutine harq_combine_rv1


  subroutine harq_combine_rv2(freq, llr174_rv2, timestamp_ms, &
       message77, nharderror, decode_ok)
  ! Combine RV0 + RV1 + RV2 LLRs and attempt decode with LDPC(348,91).
  !
  ! Full 348-bit combined vector:
  !   Bits 1-174:   llr_rv0 + RV1 repeat + RV2 repeat
  !   Bits 175-261: RV1 new parity
  !   Bits 262-348: RV2 new parity

    implicit none
    real, intent(in) :: freq
    real, intent(in) :: llr174_rv2(174)
    integer, intent(in) :: timestamp_ms
    integer*1, intent(out) :: message77(77)
    integer, intent(out) :: nharderror
    logical, intent(out) :: decode_ok

    real :: llr_combined(N_EXT2)         !348 combined LLRs
    integer :: islot, iter, i

    decode_ok = .false.
    nharderror = -1
    message77 = 0

    if(.not.harq_initialized) call harq_init()

    islot = find_slot(freq, timestamp_ms)
    if(islot .le. 0) return
    if(slots(islot)%rv_count .lt. 1) return   !Need at least RV0+RV1 stored

    ! Build 348-element combined vector
    ! Start with RV0 + RV1 chase contributions on first 174 bits
    llr_combined(1:174) = slots(islot)%llr_rv0(1:174)

    ! Add RV1 repeated systematic (positions 1-87)
    do i = 1, 87
       llr_combined(i) = llr_combined(i) + slots(islot)%llr_rv1_repeat(i)
    enddo

    ! Add RV2 repeated systematic (positions 1-87)
    do i = 1, 87
       llr_combined(i) = llr_combined(i) + llr174_rv2(87 + i)
    enddo

    ! RV1 new parity -> positions 175-261
    llr_combined(175:261) = slots(islot)%llr_rv1_new(1:87)

    ! RV2 new parity -> positions 262-348
    llr_combined(262:348) = llr174_rv2(1:87)

    ! Store RV2 and update slot
    slots(islot)%llr_rv2_new = llr174_rv2(1:87)
    slots(islot)%llr_rv2_repeat = llr174_rv2(88:174)
    slots(islot)%rv_count = 2
    slots(islot)%timestamp_ms = timestamp_ms

    ! Decode with full mother code LDPC(348,91)
    call bpdecode_ext(llr_combined, N_EXT2, 50, message77, nharderror, iter)

    if(nharderror .ge. 0) then
       decode_ok = .true.
    endif

    ! Always clear slot after RV2 (final attempt)
    call harq_clear_slot(islot)

  end subroutine harq_combine_rv2


  subroutine harq_expire(current_time_ms)
  ! Expire stale buffers older than EXPIRY_MS.

    implicit none
    integer, intent(in) :: current_time_ms
    integer :: i

    if(.not.harq_initialized) return

    do i = 1, MAX_HARQ_SLOTS
       if(slots(i)%active) then
          if(current_time_ms - slots(i)%timestamp_ms .gt. EXPIRY_MS) then
             call harq_clear_slot(i)
          endif
       endif
    enddo

  end subroutine harq_expire


  ! === Internal helper functions ===

  integer function find_slot(freq, timestamp_ms)
  ! Find an active slot matching the given frequency within tolerance.

    implicit none
    real, intent(in) :: freq
    integer, intent(in) :: timestamp_ms
    integer :: i

    find_slot = 0
    do i = 1, MAX_HARQ_SLOTS
       if(slots(i)%active) then
          if(abs(slots(i)%freq - freq) .le. real(FREQ_TOL_HZ)) then
             ! Check not expired
             if(timestamp_ms - slots(i)%timestamp_ms .le. EXPIRY_MS) then
                find_slot = i
                return
             else
                call harq_clear_slot(i)   !Expired, clear it
             endif
          endif
       endif
    enddo

  end function find_slot


  integer function allocate_slot(timestamp_ms)
  ! Find a free slot, or recycle the oldest expired one.

    implicit none
    integer, intent(in) :: timestamp_ms
    integer :: i, oldest_slot, oldest_time

    allocate_slot = 0

    ! First pass: find inactive slot
    do i = 1, MAX_HARQ_SLOTS
       if(.not.slots(i)%active) then
          allocate_slot = i
          return
       endif
    enddo

    ! Second pass: recycle oldest slot
    oldest_slot = 1
    oldest_time = slots(1)%timestamp_ms
    do i = 2, MAX_HARQ_SLOTS
       if(slots(i)%timestamp_ms .lt. oldest_time) then
          oldest_time = slots(i)%timestamp_ms
          oldest_slot = i
       endif
    enddo
    call harq_clear_slot(oldest_slot)
    allocate_slot = oldest_slot

  end function allocate_slot


  subroutine harq_clear_slot(islot)
    implicit none
    integer, intent(in) :: islot
    slots(islot)%active = .false.
    slots(islot)%rv_count = -1
    slots(islot)%freq = 0.0
    slots(islot)%timestamp_ms = 0
    slots(islot)%llr_rv0 = 0.0
    slots(islot)%llr_rv1_new = 0.0
    slots(islot)%llr_rv1_repeat = 0.0
    slots(islot)%llr_rv2_new = 0.0
    slots(islot)%llr_rv2_repeat = 0.0
  end subroutine harq_clear_slot

end module ir_harq_combine_mod
