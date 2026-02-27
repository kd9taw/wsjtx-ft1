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
module timing_results
! Module-level storage for decode results (accessible from callback)
  implicit none
  integer, parameter :: MAX_SIGNALS = 20
  integer :: ndecoded = 0
  character(len=37) :: dec_msgs(MAX_SIGNALS)
  real :: dec_freqs(MAX_SIGNALS)
  integer :: dec_snrs(MAX_SIGNALS)
  real :: dec_dts(MAX_SIGNALS)
end module timing_results


subroutine timing_callback(this,sync,snr,dt,freq,decoded,nap,qual)
  use ft1_decode
  use timing_results
  implicit none
  class(ft1_decoder), intent(inout) :: this
  real, intent(in) :: sync
  integer, intent(in) :: snr
  real, intent(in) :: dt
  real, intent(in) :: freq
  character(len=37), intent(in) :: decoded
  integer, intent(in) :: nap
  real, intent(in) :: qual

  if(ndecoded .lt. MAX_SIGNALS) then
     ndecoded = ndecoded + 1
     dec_msgs(ndecoded) = decoded
     dec_freqs(ndecoded) = freq
     dec_snrs(ndecoded) = snr
     dec_dts(ndecoded) = dt
  endif

  return
end subroutine timing_callback


program ft1_timing_test

! FT1 Full-Pipeline Timing Budget Validation
!
! Generates N FT1 signals at different frequencies, combines them with
! AWGN, then runs the full decode pipeline (sync -> candidates -> turbo
! -> OSD -> AP -> subtraction, 3 passes) and measures wall-clock time.
!
! Usage: ft1_timing_test [nsignals snrdb ntrials]
!   Defaults: 8 -12.0 1
!
! Budget: decode must complete within 500 ms

  use ft1_decode
  use timer_module, only: timer
  use timer_impl, only: init_timer
  use timing_results
  use packjt77
  include 'ft1_params.f90'

  integer, parameter :: BUDGET_MS = 500

  type(ft1_decoder) :: decoder

  character*37 msgs(MAX_SIGNALS), msgsent
  character*12 arg
  integer itone(NN)
  integer*1 msgbits(77)
  integer*2 iwave(NMAX)
  real wave(NMAX), dd(NMAX)
  real freqs(MAX_SIGNALS), dt_offsets(MAX_SIGNALS)

  ! Timing
  integer :: count_start, count_end, count_rate
  real :: elapsed_ms, margin
  real :: t_min, t_max, t_sum
  integer :: ntotal_decoded

  ! Callback interface
  external timing_callback

  ! Test messages (diverse callsigns)
  character(len=37) :: test_msgs(MAX_SIGNALS)
  data test_msgs / &
       'CQ W1ABC FN42       ', &
       'CQ K2DEF EN74       ', &
       'CQ VE3GHI FN03      ', &
       'CQ JA1KLM PM95      ', &
       'CQ G4NOP IO91       ', &
       'CQ DL5QRS JN48      ', &
       'CQ VK6TUV OF86      ', &
       'CQ ZL2WXY RF72      ', &
       'CQ UA3ABC KO85      ', &
       'CQ PY2DEF GG87      ', &
       'CQ 9A1GHI JN75      ', &
       'CQ OH3JKL KP20      ', &
       'CQ HL5MNO PM37      ', &
       'CQ LU8PQR FF99      ', &
       'CQ YB0STU OI33      ', &
       'CQ A61VWX LL75      ', &
       'CQ HS0YZA NK99      ', &
       'CQ ZS6BCD KG44      ', &
       'CQ SV9EFG KM25      ', &
       'CQ EA8HIJ IL28      ' /

  ! Parse command-line arguments
  nsignals = 8
  snrdb = -12.0
  ntrials = 1
  nargs = iargc()
  if(nargs.ge.1) then
     call getarg(1,arg)
     read(arg,*) nsignals
  endif
  if(nargs.ge.2) then
     call getarg(2,arg)
     read(arg,*) snrdb
  endif
  if(nargs.ge.3) then
     call getarg(3,arg)
     read(arg,*) ntrials
  endif
  nsignals = min(nsignals, MAX_SIGNALS)

  write(*,'(a)') '=== FT1 Full-Pipeline Timing Test ==='
  write(*,'(a,i2,a,f6.1,a,i2,a)') 'Signals: ',nsignals, &
       '  SNR: ',snrdb,' dB  Depth: 3  Trials: ',ntrials
  write(*,'(a)') 'Band: 200-4900 Hz  Budget: 500 ms'
  write(*,*)

  ! Initialize timer infrastructure for per-stage breakdown
  call init_timer('timer.out')
  call sgran()

  fs = 12000.0
  bandwidth_ratio = 2500.0 / (fs/2.0)
  sig = sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)

  ! Set up signal frequencies (400 Hz spacing, starting at 600 Hz)
  ! For 1 signal, use 1500 Hz (band center)
  write(*,'(a)') 'Generating signals...'
  if(nsignals .eq. 1) then
     freqs(1) = 1500.0
  else
     do i = 1, nsignals
        freqs(i) = 600.0 + (i-1) * 400.0
        if(freqs(i) .gt. 4500.0) freqs(i) = 4500.0 - (i-9)*200.0
     enddo
  endif

  t_min = 1.0e30
  t_max = 0.0
  t_sum = 0.0
  ntotal_decoded = 0

  do itrial = 1, ntrials

     ! Build combined audio buffer
     dd = 0.0

     do i = 1, nsignals
        ! Nominal TX start at 0.25s into the 4.0s buffer (3000 samples)
        ! ft1_sync searches j0=20..36 (~0.18-0.32s), fine sync covers wider
        ! Plus random jitter within +/- 20 ms
        call random_number(rr)
        dt_offsets(i) = 0.25 + (rr - 0.5) * 0.040   ! 0.23 to 0.27 s
        noff = nint(dt_offsets(i) * fs)

        ! Encode message
        msgs(i) = test_msgs(i)
        call genft1(msgs(i), 0, msgsent, msgbits, itone)
        msgs(i) = msgsent

        ! Generate CPM waveform
        nwave = NMAX
        wave = 0.0
        call gen_ft1wave(itone, NN, NSPS_NUM, NSPS_DEN, fs, freqs(i), &
             wave, nwave)

        ! Add to buffer at TX offset with scaling
        do j = 1, nwave
           k = j + noff
           if(k .ge. 1 .and. k .le. NMAX) then
              dd(k) = dd(k) + sig * wave(j)
           endif
        enddo

        if(itrial .eq. 1) then
           write(*,'(a,i2,a,a24,f8.1,a,f7.3)') '  ',i,': ', &
                msgs(i)(1:24), freqs(i),' Hz  dt=',dt_offsets(i)
        endif
        if(nsignals .eq. 1) then
           write(*,'(a,i3,a,f6.3,a,i4)') &
                'T',itrial,' dt=',dt_offsets(i),' ib_exp=', &
                nint(dt_offsets(i)*fs/real(NDOWN))
        endif
     enddo

     ! Add AWGN (unit variance)
     do j = 1, NMAX
        dd(j) = dd(j) + gran()
     enddo

     ! Convert to int16 (scale to use dynamic range reasonably)
     do j = 1, NMAX
        iwave(j) = nint(min(32767.0, max(-32768.0, dd(j) * 1000.0)))
     enddo

     ! Run full decode pipeline
     ndecoded = 0
     call system_clock(count_start, count_rate)

     call timer('ft1_time',0)
     call decoder%decode(timing_callback, iwave, 0, 1500, &
          200, 4900, 3, .false., 0, 'W1ABC       ', '            ')
     call timer('ft1_time',1)

     call system_clock(count_end)
     elapsed_ms = 1000.0 * real(count_end - count_start) / real(count_rate)

     t_sum = t_sum + elapsed_ms
     if(elapsed_ms .lt. t_min) t_min = elapsed_ms
     if(elapsed_ms .gt. t_max) t_max = elapsed_ms
     ntotal_decoded = ntotal_decoded + ndecoded

     if(itrial .eq. 1) then
        write(*,*)
        write(*,'(a)') 'Decode results:'
        do i = 1, ndecoded
           write(*,'(a,i2,a,f7.1,a,i4,a,f5.2,a,a)') &
                '  #',i,' f=',dec_freqs(i), &
                '  snr=',dec_snrs(i), &
                '  dt=',dec_dts(i), &
                '  ',trim(dec_msgs(i))
        enddo
     endif
  enddo

  ! Print timing summary
  write(*,*)
  write(*,'(a)') '=== RESULT ==='

  if(ntrials .eq. 1) then
     elapsed_ms = t_sum
     margin = 100.0 * (1.0 - elapsed_ms / real(BUDGET_MS))
     write(*,'(a,i2,a,i2,a,f7.1,a,f7.1,a)',advance='no') &
          'Decoded: ',ndecoded,'/',nsignals, &
          '  Time: ',elapsed_ms,' ms  Budget: ',real(BUDGET_MS),' ms  '
     if(elapsed_ms .le. real(BUDGET_MS)) then
        write(*,'(a,f5.1,a)') 'PASS (margin: ',margin,'%)'
     else
        write(*,'(a,f5.1,a)') 'FAIL (over by ',abs(margin),'%)'
     endif
  else
     t_avg = t_sum / real(ntrials)
     margin = 100.0 * (1.0 - t_max / real(BUDGET_MS))
     write(*,'(a,f7.1,a,f7.1,a,f7.1,a)') &
          'Time (ms): avg=',t_avg,' min=',t_min,' max=',t_max
     write(*,'(a,f7.1,a)',advance='no') &
          'Budget: ',real(BUDGET_MS),' ms  '
     if(t_max .le. real(BUDGET_MS)) then
        write(*,'(a,f5.1,a)') 'PASS (worst-case margin: ',margin,'%)'
     else
        write(*,'(a,f5.1,a)') 'FAIL (worst-case over by ',abs(margin),'%)'
     endif
  endif

  if(ntrials .gt. 1) then
     write(*,'(a,i5,a,i5,a,f6.1,a)') &
          'Decoded: ',ntotal_decoded,'/',nsignals*ntrials, &
          '  Rate: ',100.0*real(ntotal_decoded)/real(nsignals*ntrials),'%'
  endif

  ! Trigger timer summary (k>100 -> print and reset)
  call timer('ft1_time',101)

  write(*,*)
  write(*,'(a)') 'Timer breakdown written to: timer.out'

end program ft1_timing_test
