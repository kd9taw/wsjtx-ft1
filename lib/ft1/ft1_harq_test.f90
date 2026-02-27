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
program ft1_harq_test

! FT1 IR-HARQ AWGN integration test -- validates LDPC(348,91) code quality
! in isolation using simplified BPSK channel (no CPM demodulation).
!
! Tests three redundancy versions:
!   1TX (RV0):       LDPC(174,91) rate 0.523
!   2TX (RV0+RV1):   LDPC(261,91) rate 0.349
!   3TX (RV0+RV1+2): LDPC(348,91) rate 0.261
!
! Usage: ft1_harq_test [ebn0_start ebn0_end ntrials]
!   Defaults: 6 0 500
!
! Expected thresholds (50% decode, Eb/N0 in dB):
!   1TX (R=0.52): ~1.5 dB  (Shannon limit ~0.37 dB)
!   2TX (R=0.35): ~0.0 dB  (Shannon limit ~-0.79 dB)
!   3TX (R=0.26): ~-0.5 dB (Shannon limit ~-1.27 dB)

  use ldpc348_91_mod
  use packjt77
  implicit none

  character*37 msg37, msgsent37
  character arg*12
  integer*1 msgbits(77)
  integer*1 cw348(N_MOTHER), cw174(N_BASE)
  integer*1 decoded77(77)
  integer*1 apmask(N_BASE)
  integer*1 cw_base(N_BASE)
  integer itone(99)

  real :: snr_start, snr_end, snrdb
  real :: sigma, sigma2
  real :: tx_bpsk(N_MOTHER)
  real :: llr_rv0(N_BASE), llr_rv1_new(N_NEW_PER_RV)
  real :: llr_rv1_repeat(N_NEW_PER_RV)
  real :: llr_rv2_new(N_NEW_PER_RV), llr_rv2_repeat(N_NEW_PER_RV)
  real :: llr_combined(N_MOTHER)

  integer :: ntrials, nargs, nsnr, isnr, itrial
  integer :: ndec_1tx, ndec_2tx, ndec_3tx
  integer :: nharderror, iter, ncheck
  integer :: i
  real :: gran  ! External function

  ! Test message
  msg37 = 'CQ W9XYZ EN37'

  ! Default arguments (Eb/N0 in dB, sweep high to low)
  snr_start = 6.0
  snr_end = -2.0
  ntrials = 500

  ! Parse command-line arguments
  nargs = iargc()
  if(nargs .ge. 1) then
     call getarg(1, arg)
     read(arg, *) snr_start
  endif
  if(nargs .ge. 2) then
     call getarg(2, arg)
     read(arg, *) snr_end
  endif
  if(nargs .ge. 3) then
     call getarg(3, arg)
     read(arg, *) ntrials
  endif

  ! Initialize
  call init_ldpc348_91()
  apmask = 0

  ! Encode test message: genft1 gives us 77-bit payload and tones
  call genft1(msg37, 0, msgsent37, msgbits, itone)

  ! encode174_91 adds CRC-14 internally, producing 174-bit codeword
  call encode174_91(msgbits, cw174)

  ! encode348_91 does the same base encoding + extension parity
  ! It takes 91-bit message (77 msg + 14 CRC), but constructs CRC internally
  ! Actually, encode348_91 takes message91 and calls encode174_91(message91(1:77))
  ! So we need the 91-bit message. But encode174_91 constructs it from 77 bits.
  ! For our test, we just need cw174 from encode174_91, then manually compute
  ! extension parity. Let's call encode348_91 with proper input.
  block
    integer*1 :: msg91(K_MOTHER)
    ! Build the 91-bit message (77 msg + CRC14) by reading from cw174
    msg91(1:K_MOTHER) = cw174(1:K_MOTHER)  ! First 91 bits = systematic part
    call encode348_91(msg91, cw348)
  end block

  ! Verify base encoding matches
  if(any(cw348(1:N_BASE) .ne. cw174)) then
     write(*,*) 'ERROR: Base encoding mismatch!'
     stop 1
  endif

  write(*,'(a,a37)') 'Test message: ', msgsent37
  write(*,'(a,i5)') 'Trials per SNR: ', ntrials
  write(*,*)

  call sgran()  ! Seed random number generator

  ! SNR sweep
  nsnr = nint(snr_start - snr_end) + 1
  if(nsnr .gt. 50) nsnr = 50

  write(*,'(a)') ' Es/N0    1TX     2TX     3TX   (% decoded)'
  write(*,'(a)') '  (dB)   R=.52   R=.35   R=.26'
  write(*,'(a)') ' -----   -----   -----   -----'

  do isnr = 1, nsnr
     snrdb = snr_start - real(isnr - 1)

     ! Direct Eb/N0 (dB) to sigma for BPSK AWGN
     ! tx = +/-1, rx = tx + N(0, sigma^2)
     ! Eb/N0 = 1/(2*sigma^2)  => sigma^2 = 1/(2 * 10^(Eb_N0_dB/10))
     sigma2 = 1.0 / (2.0 * 10.0**(snrdb / 10.0))
     sigma = sqrt(sigma2)

     ndec_1tx = 0
     ndec_2tx = 0
     ndec_3tx = 0

     do itrial = 1, ntrials

        ! BPSK modulate full 348-bit codeword: 0 -> +1, 1 -> -1
        do i = 1, N_MOTHER
           tx_bpsk(i) = 1.0 - 2.0 * real(cw348(i))
        enddo

        ! === 1TX test (RV0 only: bits 1-174) ===
        do i = 1, N_BASE
           llr_rv0(i) = -2.0 * (tx_bpsk(i) + sigma * gran()) / sigma2
        enddo

        call bpdecode174_91(llr_rv0, apmask, 50, decoded77, cw_base, &
             nharderror, iter, ncheck)
        if(nharderror .ge. 0) ndec_1tx = ndec_1tx + 1

        ! === 2TX test (RV0 + RV1) ===
        ! RV1 transmits: 87 new parity (175-261) + 87 systematic repeat (1-87)

        ! New parity bits (175-261): fresh noise realization
        do i = 1, N_NEW_PER_RV
           llr_rv1_new(i) = -2.0 * (tx_bpsk(N_BASE + i) + sigma * gran()) &
                / sigma2
        enddo

        ! Systematic repeat (bits 1-87): fresh noise, chase combine with RV0
        do i = 1, N_NEW_PER_RV
           llr_rv1_repeat(i) = -2.0 * (tx_bpsk(i) + sigma * gran()) / sigma2
        enddo

        ! Build combined LLR for LDPC(261,91) decode
        llr_combined(1:N_BASE) = llr_rv0(1:N_BASE)
        do i = 1, N_NEW_PER_RV
           llr_combined(i) = llr_combined(i) + llr_rv1_repeat(i)
        enddo
        llr_combined(N_BASE+1:N_EXT1) = llr_rv1_new(1:N_NEW_PER_RV)

        call bpdecode_ext(llr_combined, N_EXT1, 50, decoded77, &
             nharderror, iter)
        if(nharderror .ge. 0) ndec_2tx = ndec_2tx + 1

        ! === 3TX test (RV0 + RV1 + RV2) ===
        ! RV2 transmits: 87 new parity (262-348) + 87 systematic repeat (1-87)

        ! New parity bits (262-348): fresh noise realization
        do i = 1, N_NEW_PER_RV
           llr_rv2_new(i) = -2.0 * (tx_bpsk(N_EXT1 + i) + sigma * gran()) &
                / sigma2
        enddo

        ! Systematic repeat: more chase combining on bits 1-87
        do i = 1, N_NEW_PER_RV
           llr_rv2_repeat(i) = -2.0 * (tx_bpsk(i) + sigma * gran()) / sigma2
        enddo

        ! Extend combined LLR to LDPC(348,91)
        do i = 1, N_NEW_PER_RV
           llr_combined(i) = llr_combined(i) + llr_rv2_repeat(i)
        enddo
        llr_combined(N_EXT1+1:N_EXT2) = llr_rv2_new(1:N_NEW_PER_RV)

        call bpdecode_ext(llr_combined, N_EXT2, 50, decoded77, &
             nharderror, iter)
        if(nharderror .ge. 0) ndec_3tx = ndec_3tx + 1

     enddo  ! itrial

     write(*,'(f7.1, 3f8.1)') snrdb, &
          100.0 * real(ndec_1tx) / real(ntrials), &
          100.0 * real(ndec_2tx) / real(ntrials), &
          100.0 * real(ndec_3tx) / real(ntrials)

  enddo  ! isnr

  write(*,*)
  write(*,'(a)') 'Test complete.'

end program ft1_harq_test
