program ft1_harq_fading_test

! FT1 IR-HARQ fading channel test -- exercises turbo decoder + HARQ combining
! under Watterson HF fading with independent fading per retransmission.
!
! Each RV (0,1,2) is transmitted in a separate T/R period with independent
! fading realization. The receiver pipeline:
!   RV0: Full turbo decode. If success -> done.
!   RV1: BCJR-only demod -> combine with RV0 LLRs -> bpdecode_ext(261)
!   RV2: BCJR-only demod -> combine with RV0+RV1 -> bpdecode_ext(348)
!
! Usage: ft1_harq_fading_test [snr_start snr_end ntrials fspread delay]
!   Defaults: 0 -20 100 1.0 1.0
!   Set fspread=0, delay=0 for AWGN baseline.

  use packjt77
  use ldpc348_91_mod
  include 'ft1_params.f90'

  parameter (NDMAX=NMAX/NDOWN)

  character*37 msg37, msgsent37
  character arg*12
  integer itone_rv0(NN), itone_rv1(NN), itone_rv2(NN)
  integer itmp(ND)
  integer*1 msgbits(77)
  integer*1 message91(91)
  integer*1 cw174(174), cw348(N_MOTHER)
  integer*1 tx174_rv1(174), tx174_rv2(174)
  integer*1 interleaved(174)
  integer*1 msg91(K_MOTHER)
  integer*1 decoded77(77)

  ! OSD fallback variables for HARQ combined decode
  integer*1 :: apmask_osd(174), cw_osd(174), msg91_osd(91)
  integer :: nhd_osd, i_osd
  real :: dmin_osd

  ! Timing variables
  real :: t_start, t_end, t_rv0, t_rv1, t_rv2
  real :: t_rv0_sum, t_rv1_sum, t_rv2_sum
  integer :: n_rv1_calls, n_rv2_calls

  real wave_rv0(NMAX), wave_rv1(NMAX), wave_rv2(NMAX)
  real llr_out(174), llr_rv0(174), llr_rv1(174), llr_rv2(174)
  real llr_combined(N_MOTHER)

  complex cd_rv0(0:NDMAX-1), cd_rv1(0:NDMAX-1), cd_rv2(0:NDMAX-1)
  complex c_analytic_rv0(0:NMAX-1)
  complex c_analytic_rv1(0:NMAX-1)
  complex c_analytic_rv2(0:NMAX-1)

  integer icos_rv1(4), icos_rv2(4)
  data icos_rv1/1,3,2,0/
  data icos_rv2/3,0,2,1/

  msg37 = 'CQ W9XYZ EN37'
  f0 = 1500.0

  ! Parse command-line arguments
  snr_start = 0.0
  snr_end = -20.0
  ntrials = 100
  fspread = 1.0
  delay = 1.0
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
  if(nargs .ge. 4) then
     call getarg(4, arg)
     read(arg, *) fspread
  endif
  if(nargs .ge. 5) then
     call getarg(5, arg)
     read(arg, *) delay
  endif

  fs = 12000.0
  bandwidth_ratio = 2500.0 / (fs / 2.0)

  ! Initialize HARQ code tables
  call init_ldpc348_91()

  ! Encode test message
  call genft1(msg37, 0, msgsent37, msgbits, itone_rv0)

  ! Build LDPC(348,91) mother codeword
  call encode174_91(msgbits, cw174)
  msg91(1:K_MOTHER) = cw174(1:K_MOTHER)
  call encode348_91(msg91, cw348)

  ! Puncture for RV1 and RV2
  call puncture_rv(cw348, 1, tx174_rv1)
  call puncture_rv(cw348, 2, tx174_rv2)

  ! Generate RV1 tones: interleave -> Gray map -> insert Costas
  call ft1_interleave(tx174_rv1, interleaved, 1)
  do i = 1, ND
     is = interleaved(2*i) + 2*interleaved(2*i-1)
     if(is .le. 1) itmp(i) = is
     if(is .eq. 2) itmp(i) = 3
     if(is .eq. 3) itmp(i) = 2
  enddo
  itone_rv1(1:4)   = icos_rv1
  itone_rv1(5:47)  = itmp(1:43)
  itone_rv1(48:51) = icos_rv1
  itone_rv1(52:95) = itmp(44:87)
  itone_rv1(96:99) = icos_rv1

  ! Generate RV2 tones
  call ft1_interleave(tx174_rv2, interleaved, 1)
  do i = 1, ND
     is = interleaved(2*i) + 2*interleaved(2*i-1)
     if(is .le. 1) itmp(i) = is
     if(is .eq. 2) itmp(i) = 3
     if(is .eq. 3) itmp(i) = 2
  enddo
  itone_rv2(1:4)   = icos_rv2
  itone_rv2(5:47)  = itmp(1:43)
  itone_rv2(48:51) = icos_rv2
  itone_rv2(52:95) = itmp(44:87)
  itone_rv2(96:99) = icos_rv2

  ! Generate 3 CPM waveforms at 12 kHz
  nwave = NMAX
  wave_rv0 = 0.0
  wave_rv1 = 0.0
  wave_rv2 = 0.0
  call gen_ft1wave(itone_rv0, NN, NSPS_NUM, NSPS_DEN, fs, f0, wave_rv0, nwave)
  call gen_ft1wave(itone_rv1, NN, NSPS_NUM, NSPS_DEN, fs, f0, wave_rv1, nwave)
  call gen_ft1wave(itone_rv2, NN, NSPS_NUM, NSPS_DEN, fs, f0, wave_rv2, nwave)

  ! Compute original RMS for SNR calibration
  rms_orig = sqrt(sum(wave_rv0(1:NZ)**2) / NZ)

  ! Pre-compute analytic signals (Hilbert transform) for all 3 RVs
  call make_analytic(wave_rv0, c_analytic_rv0, NMAX)
  call make_analytic(wave_rv1, c_analytic_rv1, NMAX)
  call make_analytic(wave_rv2, c_analytic_rv2, NMAX)

  write(*,'(a,a37)') 'Test message: ', msgsent37
  write(*,'(a,f7.1,a)') 'Carrier freq: ', f0, ' Hz'
  write(*,'(a,f6.2,a,f5.2,a)') 'Fading:       fspread=', fspread, &
       ' Hz, delay=', delay, ' ms'
  write(*,'(a,i4,a)') 'Trials:       ', ntrials, ' per SNR point'
  write(*,*)

  call sgran()

  ! SNR sweep
  nsnr = nint(snr_start - snr_end) + 1
  if(nsnr .gt. 50) nsnr = 50

  write(*,'(a)') '   SNR    1TX%    2TX%    3TX%'
  write(*,'(a)') '   ---    ----    ----    ----'

  do isnr = 1, nsnr
     snrdb = snr_start - real(isnr - 1)
     sig = sqrt(2 * bandwidth_ratio) * 10.0**(0.05 * snrdb)

     ndec_1tx = 0
     ndec_2tx = 0
     ndec_3tx = 0
     t_rv0_sum = 0.0
     t_rv1_sum = 0.0
     t_rv2_sum = 0.0
     n_rv1_calls = 0
     n_rv2_calls = 0

     do itrial = 1, ntrials

        ! === RV0: Full turbo decode ===
        call apply_channel(c_analytic_rv0, wave_rv0, NMAX, NZ, fs, f0, &
             fspread, delay, sig, rms_orig, cd_rv0, NDMAX)

        npts = NDMAX
        dt0 = 0.0
        snr_est = snrdb
        ntype = -1
        nharderror = -1
        dmin = 0.0
        message91 = 0
        llr_out = 0.0

        call cpu_time(t_start)
        call turbo_decode_ft1(cd_rv0, npts, f0, dt0, snr_est, llr_out, &
             message91, ntype, nharderror, dmin, 0, ncheck_out)
        call cpu_time(t_end)
        t_rv0_sum = t_rv0_sum + (t_end - t_start) * 1000.0

        if(ntype .ge. 0) then
           ndec_1tx = ndec_1tx + 1
           ndec_2tx = ndec_2tx + 1
           ndec_3tx = ndec_3tx + 1
           cycle  ! Next trial
        endif

        ! RV0 failed -- store turbo-refined LLRs
        llr_rv0 = llr_out

        ! === RV1: BCJR-only demod + HARQ combine ===
        call apply_channel(c_analytic_rv1, wave_rv1, NMAX, NZ, fs, f0, &
             fspread, delay, sig, rms_orig, cd_rv1, NDMAX)

        call cpu_time(t_start)
        call ft1_demod_bcjr(cd_rv1, NDMAX, snrdb, 1, llr_rv1)

        ! Combine: LDPC(261,91)
        ! llr_rv1(1:87) = new parity LLRs, llr_rv1(88:174) = systematic repeat
        llr_combined(1:174) = llr_rv0(1:174)
        do i = 1, 87
           llr_combined(i) = llr_combined(i) + llr_rv1(87 + i)
        enddo
        llr_combined(175:261) = llr_rv1(1:87)

        call bpdecode_ext(llr_combined, N_EXT1, 50, decoded77, &
             nharderror, iter)
        ! OSD fallback on combined base-code LLRs
        if(nharderror .lt. 0) then
           apmask_osd = 0
           do i_osd = 1, 4
              call osd174_91(llr_combined, 91, apmask_osd, i_osd, &
                   msg91_osd, cw_osd, nhd_osd, dmin_osd)
              if(nhd_osd .ge. 0) then
                 nharderror = nhd_osd
                 exit
              endif
           enddo
        endif
        call cpu_time(t_end)
        t_rv1_sum = t_rv1_sum + (t_end - t_start) * 1000.0
        n_rv1_calls = n_rv1_calls + 1
        if(nharderror .ge. 0) then
           ndec_2tx = ndec_2tx + 1
           ndec_3tx = ndec_3tx + 1
           cycle  ! Next trial
        endif

        ! === RV2: BCJR-only demod + HARQ combine ===
        call apply_channel(c_analytic_rv2, wave_rv2, NMAX, NZ, fs, f0, &
             fspread, delay, sig, rms_orig, cd_rv2, NDMAX)

        call cpu_time(t_start)
        call ft1_demod_bcjr(cd_rv2, NDMAX, snrdb, 2, llr_rv2)

        ! Extend to LDPC(348,91)
        do i = 1, 87
           llr_combined(i) = llr_combined(i) + llr_rv2(87 + i)
        enddo
        llr_combined(262:348) = llr_rv2(1:87)

        call bpdecode_ext(llr_combined, N_EXT2, 50, decoded77, &
             nharderror, iter)
        ! OSD fallback on combined base-code LLRs
        if(nharderror .lt. 0) then
           apmask_osd = 0
           do i_osd = 1, 4
              call osd174_91(llr_combined, 91, apmask_osd, i_osd, &
                   msg91_osd, cw_osd, nhd_osd, dmin_osd)
              if(nhd_osd .ge. 0) then
                 nharderror = nhd_osd
                 exit
              endif
           enddo
        endif
        call cpu_time(t_end)
        t_rv2_sum = t_rv2_sum + (t_end - t_start) * 1000.0
        n_rv2_calls = n_rv2_calls + 1
        if(nharderror .ge. 0) then
           ndec_3tx = ndec_3tx + 1
        endif

     enddo  ! itrial

     write(*,'(f6.1, 3f8.1, a, f6.1, a)', advance='no') snrdb, &
          100.0 * real(ndec_1tx) / real(ntrials), &
          100.0 * real(ndec_2tx) / real(ntrials), &
          100.0 * real(ndec_3tx) / real(ntrials), &
          '  RV0=', t_rv0_sum / real(ntrials), 'ms'
     if(n_rv1_calls .gt. 0) then
        write(*,'(a, f6.1, a)', advance='no') &
             ' RV1=', t_rv1_sum / real(n_rv1_calls), 'ms'
     endif
     if(n_rv2_calls .gt. 0) then
        write(*,'(a, f6.1, a)', advance='no') &
             ' RV2=', t_rv2_sum / real(n_rv2_calls), 'ms'
     endif
     write(*,*)

  enddo  ! isnr

  write(*,*)
  write(*,'(a)') 'Test complete.'

end program ft1_harq_fading_test


subroutine make_analytic(wave, c_analytic, nmax)
! Convert real waveform to analytic signal via Hilbert transform.
! Zeros negative frequencies, doubles positive frequencies.
  implicit none
  integer, intent(in)  :: nmax
  real, intent(in)     :: wave(nmax)
  complex, intent(out) :: c_analytic(0:nmax-1)
  integer :: i

  do i = 0, nmax - 1
     c_analytic(i) = cmplx(wave(i+1), 0.0)
  enddo
  call four2a(c_analytic, nmax, 1, -1, 1)       ! Forward FFT
  c_analytic(nmax/2+1:nmax-1) = cmplx(0.0, 0.0) ! Zero negative freqs
  do i = 1, nmax/2 - 1
     c_analytic(i) = 2.0 * c_analytic(i)         ! Double positive freqs
  enddo
  call four2a(c_analytic, nmax, 1, 1, 1)         ! Inverse FFT
  c_analytic = c_analytic / real(nmax)

  return
end subroutine make_analytic


subroutine apply_channel(c_analytic, wave_orig, nmax, nz, fs, f0, &
     fspread, delay, sig, rms_orig, cd_out, ndmax)
! Apply independent fading + AWGN + downsample for one RV.
  implicit none
  integer, intent(in)  :: nmax, nz, ndmax
  complex, intent(in)  :: c_analytic(0:nmax-1)
  real, intent(in)     :: wave_orig(nmax)
  real, intent(in)     :: fs, f0, fspread, delay, sig, rms_orig
  complex, intent(out) :: cd_out(0:ndmax-1)

  real    :: dd_faded(nmax), dd_noisy(nmax)
  complex :: c0(0:nmax-1)
  real    :: rms_faded, sum2, gran
  logical :: newdata
  integer :: i

  if(fspread .ne. 0.0 .or. delay .ne. 0.0) then
     c0 = c_analytic
     call watterson(c0, nmax, nz, fs, delay, fspread)
     do i = 1, nmax
        dd_faded(i) = real(c0(i-1))
     enddo
     rms_faded = sqrt(sum(dd_faded(1:nz)**2) / nz)
     if(rms_faded .gt. 0.0) dd_faded = dd_faded * (rms_orig / rms_faded)
  else
     dd_faded = wave_orig
  endif

  ! Add AWGN
  do i = 1, nmax
     dd_noisy(i) = sig * dd_faded(i) + gran()
  enddo

  ! Downsample to baseband
  newdata = .true.
  call ft1_downsample(dd_noisy, newdata, f0, cd_out)

  ! Normalize to unit power
  sum2 = sum(real(cd_out * conjg(cd_out))) / real(ndmax)
  if(sum2 .gt. 0.0) cd_out = cd_out / sqrt(sum2)

  return
end subroutine apply_channel
