! FT1 - 4-CPM turbo equalization mode for WSJT-X
! Copyright (C) 2026 Seth McCall, KD9TAW
!
! This file is part of WSJT-X.  GPLv3 (see genft1.f90 header).
!
program ft1_rv_detect_test

! Measure FT1 redundancy-version (RV) DETECTION accuracy.
!
! For each SNR and each true RV in {0,1,2}, generate the RV frame with the
! production encoder (genft1_rv), pass it through AWGN (+optional Watterson
! fading), downconvert to baseband, and detect the RV with a COHERENT CPM
! Costas-correlation discriminator (argmax over RV). Reports a confusion-derived
! accuracy + the false "RV0 tagged as RVk" rate (the rate that regresses RV0
! decodes when the tag gates fine-sync). Variant 0 also measures the legacy
! spectrogram detector (ft1_sync candidate(3)) for comparison.
!
! Discriminator metric per RV at timing i0 (rv_metric, imode):
!   G1 uses the known frame-start state (clean, no search bias).
!   G2/G3 entering states are data-dependent (unknown), handled by:
!     imode=0  max  over 16 theta=0 states of |corr|^2  (GLRT, sync1d-style)
!     imode=1  G1 only (no G2/G3 term)
!     imode=2  sum  over 16 states of |corr|^2 (incoherent combine; cancels the
!              common max-bias, leaving the clean signal term -> better at low SNR)
!
! Usage: ft1_rv_detect_test [snr_start snr_end ntrials fspread delay imode marg twin]
!   imode 0/1/2 (default 0); marg confidence margin best>marg*2nd else RV0 (0=off);
!   twin timing half-window in symbols (default 1.0).

  use packjt77
  use cpm_trellis_mod, only: init_cpm_trellis, next_state
  use matched_filter_bank_mod, only: init_matched_filters, mf_bank
  include 'ft1_params.f90'
  parameter (NDMAX=NMAX/NDOWN)

  character*37 msg37, msgsent37
  character arg*12
  integer itone(NN)
  integer*1 msgbits(77)
  integer icos_rv(0:3,0:2)

  real    wave(NMAX), dd(NMAX)
  complex c_analytic(0:NMAX-1)
  complex cd(0:NDMAX-1)

  real    savg(NH1), sbase(NH1)
  real    candidate(3,MAXCAND)

  integer conf_coh(0:2,0:2), conf_spec(0:2,0:2)
  integer nspec_miss
  integer true_rv, idet, isnr, itrial, nsnr, i, ncand, nwave
  real    snr_start, snr_end, snrdb, fspread, delay
  real    fs, f0, bandwidth_ratio, rms_orig, sig
  integer ntrials, imode, det_lo, det_hi
  real    det_marg, twin, nsps_dn0
  real    acc_coh, acc_spec, false0_coh, false0_spec
  integer ncorrect, nfalse0, ntot_rv0
  logical do_spec
  real    toff_sym, rr
  integer toff_ds, ioff_ds, ioff12, coarse_hw
  real    wave_t(NMAX)

  msg37 = 'CQ W9XYZ EN37'
  fs = 12000.0
  f0 = 1500.0
  icos_rv(0:3,0)=(/0,2,3,1/)
  icos_rv(0:3,1)=(/1,3,2,0/)
  icos_rv(0:3,2)=(/3,0,2,1/)

  snr_start=-5.0; snr_end=-16.0; ntrials=200; fspread=0.0; delay=0.0
  imode=0; det_marg=0.0; twin=1.0; toff_sym=0.0
  nargs=iargc()
  if(nargs.ge.1) then; call getarg(1,arg); read(arg,*) snr_start; endif
  if(nargs.ge.2) then; call getarg(2,arg); read(arg,*) snr_end; endif
  if(nargs.ge.3) then; call getarg(3,arg); read(arg,*) ntrials; endif
  if(nargs.ge.4) then; call getarg(4,arg); read(arg,*) fspread; endif
  if(nargs.ge.5) then; call getarg(5,arg); read(arg,*) delay; endif
  if(nargs.ge.6) then; call getarg(6,arg); read(arg,*) imode; endif
  if(nargs.ge.7) then; call getarg(7,arg); read(arg,*) det_marg; endif
  if(nargs.ge.8) then; call getarg(8,arg); read(arg,*) twin; endif
  if(nargs.ge.9) then; call getarg(9,arg); read(arg,*) toff_sym; endif

  nsps_dn0 = real(NSPS_NUM)/(real(NSPS_DEN)*real(NDOWN))
  det_lo = -nint(twin*nsps_dn0)
  det_hi =  nint(twin*nsps_dn0)
  toff_ds = nint(toff_sym*nsps_dn0)                          !random arrival-timing offset range
  coarse_hw = nint((toff_sym+0.6)*nsps_dn0)                  !two-stage coarse timing half-window
  do_spec = (imode.eq.0 .and. det_marg.eq.0.0 .and. abs(twin-2.0).lt.0.01 .and. toff_sym.eq.0.0)

  bandwidth_ratio = 2500.0 / (fs/2.0)
  call sgran()
  call init_cpm_trellis()
  call init_matched_filters(NSS)

  write(*,'(a,a)')     'Test message: ', trim(msg37)
  write(*,'(a,f6.2,a,f5.2,a)') 'Fading:       fspread=',fspread,' Hz, delay=',delay,' ms'
  write(*,'(a,i5,a)')  'Trials:       ', ntrials, ' per (SNR, true RV)'
  write(*,'(a,i2,a,f4.1,a,i4,a,i4,a,f5.2,a,f4.1)') 'Detector:     imode=',imode, &
       '  twin=',twin,' [',det_lo,',',det_hi,']  margin=',det_marg, &
       '  rand-toff(sym)=',toff_sym
  if(imode.eq.3) write(*,'(a,i4,a)') '              two-stage coarse half-window=',coarse_hw,' samp'
  write(*,*)
  write(*,'(a)') '   SNR   coh_acc%  coh_false0%   spec_acc%  spec_false0%  spec_miss%'
  write(*,'(a)') '   ---   --------  -----------   ---------  ------------  ----------'

  nsnr = nint(snr_start - snr_end) + 1
  if(nsnr.gt.60) nsnr=60

  do isnr=1,nsnr
     snrdb = snr_start - real(isnr-1)
     sig = sqrt(2.0*bandwidth_ratio) * 10.0**(0.05*snrdb)
     conf_coh=0; conf_spec=0; nspec_miss=0

     do true_rv=0,2
        call genft1_rv(msg37, 0, true_rv, msgsent37, msgbits, itone)
        wave=0.0; nwave=NMAX
        call gen_ft1wave(itone, NN, NSPS_NUM, NSPS_DEN, fs, f0, wave, nwave)
        rms_orig = sqrt(sum(wave(1:NZ)**2)/NZ)
        if(fspread.ne.0.0 .or. delay.ne.0.0) call make_analytic(wave, c_analytic, NMAX)

        do itrial=1,ntrials
           ! Optional random arrival-timing offset: shift the frame in the buffer
           ! by a random amount (downsampled samples) to test timing robustness.
           ! (Offset runs are AWGN; offset+fading is not combined here.)
           if(toff_ds.gt.0) then
              call random_number(rr)
              ioff_ds = nint((rr-0.5)*2.0*real(toff_ds))
              ioff12 = ioff_ds*NDOWN
              wave_t = 0.0
              do i=1,nwave
                 if(i+ioff12.ge.1 .and. i+ioff12.le.NMAX) wave_t(i+ioff12)=wave(i)
              enddo
           else
              ioff_ds = 0
              wave_t = wave
           endif

           call rv_channel(c_analytic, wave_t, NMAX, NZ, fs, f0, fspread, delay, &
                sig, rms_orig, cd, dd, NDMAX)

           if(imode.eq.-1) then
              call ft1_rv_detect(cd, NDMAX, ioff_ds, idet)   ! PRODUCTION routine
           else if(imode.eq.3) then
              call detect_two_stage(cd, coarse_hw, det_marg, idet)
           else
              call detect_coherent(cd, det_lo, det_hi, imode, det_marg, idet)
           endif
           conf_coh(true_rv, idet) = conf_coh(true_rv, idet) + 1

           if(do_spec) then
              candidate=0.0; ncand=0
              call ft1_sync(dd, 200.0, 2900.0, 1.0, nint(f0), MAXCAND, savg, &
                   candidate, ncand, sbase)
              idet = spec_rv_near(candidate, ncand, f0)
              if(idet.lt.0) then
                 nspec_miss = nspec_miss + 1
              else
                 conf_spec(true_rv, idet) = conf_spec(true_rv, idet) + 1
              endif
           endif
        enddo
     enddo

     ncorrect = conf_coh(0,0)+conf_coh(1,1)+conf_coh(2,2)
     acc_coh = 100.0*real(ncorrect)/real(3*ntrials)
     ntot_rv0 = sum(conf_coh(0,0:2)); nfalse0 = conf_coh(0,1)+conf_coh(0,2)
     false0_coh = 100.0*real(nfalse0)/real(max(1,ntot_rv0))

     ncorrect = conf_spec(0,0)+conf_spec(1,1)+conf_spec(2,2)
     acc_spec = 100.0*real(ncorrect)/real(3*ntrials)
     ntot_rv0 = sum(conf_spec(0,0:2)); nfalse0 = conf_spec(0,1)+conf_spec(0,2)
     false0_spec = 100.0*real(nfalse0)/real(max(1,ntot_rv0))

     write(*,'(f6.1,2x,f7.1,3x,f8.1,4x,f8.1,5x,f8.1,4x,f8.1)') snrdb, &
          acc_coh, false0_coh, acc_spec, false0_spec, &
          100.0*real(nspec_miss)/real(3*ntrials)
  enddo

  write(*,*)
  write(*,'(a)') 'Test complete.'

contains

  subroutine rv_channel(c_analytic, wave_orig, nmx, nz, fs, f0, fspread, delay, &
       sig, rms_orig, cd_out, dd_out, ndmx)
    integer, intent(in)  :: nmx, nz, ndmx
    complex, intent(in)  :: c_analytic(0:nmx-1)
    real,    intent(in)  :: wave_orig(nmx), fs, f0, fspread, delay, sig, rms_orig
    complex, intent(out) :: cd_out(0:ndmx-1)
    real,    intent(out) :: dd_out(nmx)
    real    :: dd_faded(nmx)
    complex :: c0(0:nmx-1)
    real    :: rms_faded, s2, gran
    logical :: newdata
    integer :: i
    if(fspread.ne.0.0 .or. delay.ne.0.0) then
       c0 = c_analytic
       call watterson(c0, nmx, nz, fs, delay, fspread)
       do i=1,nmx
          dd_faded(i)=real(c0(i-1))
       enddo
       rms_faded=sqrt(sum(dd_faded(1:nz)**2)/nz)
       if(rms_faded.gt.0.0) dd_faded=dd_faded*(rms_orig/rms_faded)
    else
       dd_faded=wave_orig
    endif
    do i=1,nmx
       dd_out(i) = sig*dd_faded(i) + gran()
    enddo
    newdata=.true.
    call ft1_downsample(dd_out, newdata, f0, cd_out)
    s2 = sum(real(cd_out*conjg(cd_out)))/real(ndmx)
    if(s2.gt.0.0) cd_out = cd_out/sqrt(s2)
  end subroutine rv_channel

  complex function corr_costas(cd, istart, icos, s0)
  ! Coherent correlation of 4 Costas sync symbols starting at trellis state s0.
    complex, intent(in) :: cd(0:NDMAX-1)
    integer, intent(in) :: istart, icos(0:3), s0
    complex :: csync(4*NSS), z
    integer :: s_state, k, i, j
    real    :: scale_amp
    scale_amp = sqrt(real(NSS))
    s_state = s0; k = 1
    do i=0,3
       do j=1,NSS
          csync(k) = scale_amp*mf_bank(j, s_state, icos(i)); k = k+1
       enddo
       s_state = next_state(s_state, icos(i))
    enddo
    z = cmplx(0.0,0.0)
    if(istart.ge.0 .and. istart+4*NSS-1.le.NDMAX-1) then
       z = sum(cd(istart:istart+4*NSS-1)*conjg(csync(1:4*NSS)))
    endif
    corr_costas = z
  end function corr_costas

  real function rv_metric(cd, i0, icos, imode)
  ! Coherent CPM Costas detection metric at timing i0 for Costas `icos`.
    complex, intent(in) :: cd(0:NDMAX-1)
    integer, intent(in) :: i0, icos(0:3), imode
    complex :: z
    real    :: m, best, ssum, nsps_dn
    integer :: i1, i2, i3, s, ig, grp
    nsps_dn = real(NSPS_NUM)/(real(NSPS_DEN)*real(NDOWN))
    i1 = i0
    i2 = i0 + nint(47.0*nsps_dn)
    i3 = i0 + nint(95.0*nsps_dn)
    z = corr_costas(cd, i1, icos, 1)             ! G1: known start state
    m = real(z)**2 + aimag(z)**2
    if(imode.ne.1) then
       do ig=1,2
          if(ig.eq.1) grp=i2
          if(ig.eq.2) grp=i3
          if(imode.eq.0) then
             best=0.0
             do s=1,16
                z=corr_costas(cd, grp, icos, s)
                best=max(best, real(z)**2+aimag(z)**2)
             enddo
             m=m+best
          else
             ssum=0.0
             do s=1,16
                z=corr_costas(cd, grp, icos, s)
                ssum=ssum + real(z)**2+aimag(z)**2
             enddo
             m=m+ssum
          endif
       enddo
    endif
    rv_metric = m
  end function rv_metric

  subroutine detect_coherent(cd, i0lo, i0hi, imode, marg, irv_det)
  ! argmax over RV of the best (over a timing window) rv_metric. If marg>1,
  ! accept the winner only if best>marg*second, else fall back to RV0.
    complex, intent(in)  :: cd(0:NDMAX-1)
    integer, intent(in)  :: i0lo, i0hi, imode
    real,    intent(in)  :: marg
    integer, intent(out) :: irv_det
    real    :: smet(0:2), v, best, second
    integer :: rv, i0, iwin
    do rv=0,2
       best=0.0
       do i0=i0lo,i0hi
          v = rv_metric(cd, i0, icos_rv(0:3,rv), imode)
          if(v.gt.best) best=v
       enddo
       smet(rv)=best
    enddo
    iwin=0
    if(smet(1).gt.smet(iwin)) iwin=1
    if(smet(2).gt.smet(iwin)) iwin=2
    if(marg.gt.1.0) then
       second=-1.0
       do rv=0,2
          if(rv.ne.iwin .and. smet(rv).gt.second) second=smet(rv)
       enddo
       if(smet(iwin) .lt. marg*second) iwin=0
    endif
    irv_det=iwin
  end subroutine detect_coherent

  subroutine detect_two_stage(cd, hw, marg, irv_det)
  ! Self-contained robust detector. Stage 1: find frame timing RV-agnostically
  ! as the peak of the sum-over-RV sum-states metric over a wide window. Stage 2:
  ! discriminate RV in a narrow window at that timing (sum-states + margin).
    complex, intent(in)  :: cd(0:NDMAX-1)
    integer, intent(in)  :: hw
    real,    intent(in)  :: marg
    integer, intent(out) :: irv_det
    real    :: smet(0:2), tsum, best_t, v, best, second
    integer :: i0, i0star, rv, iwin
    ! Stage 1: RV-agnostic coarse timing
    best_t=-1.0; i0star=0
    do i0=-hw,hw
       tsum = rv_metric(cd,i0,icos_rv(0:3,0),2) + rv_metric(cd,i0,icos_rv(0:3,1),2) &
            + rv_metric(cd,i0,icos_rv(0:3,2),2)
       if(tsum.gt.best_t) then; best_t=tsum; i0star=i0; endif
    enddo
    ! Stage 2: narrow RV discrimination at the timing peak (+-0.5 symbol)
    do rv=0,2
       best=0.0
       do i0=i0star-4,i0star+4
          v=rv_metric(cd,i0,icos_rv(0:3,rv),2)
          if(v.gt.best) best=v
       enddo
       smet(rv)=best
    enddo
    iwin=0
    if(smet(1).gt.smet(iwin)) iwin=1
    if(smet(2).gt.smet(iwin)) iwin=2
    if(marg.gt.1.0) then
       second=-1.0
       do rv=0,2
          if(rv.ne.iwin .and. smet(rv).gt.second) second=smet(rv)
       enddo
       if(smet(iwin).lt.marg*second) iwin=0
    endif
    irv_det=iwin
  end subroutine detect_two_stage

  subroutine make_analytic(wave, c_analytic, nmx)
    integer, intent(in)  :: nmx
    real,    intent(in)  :: wave(nmx)
    complex, intent(out) :: c_analytic(0:nmx-1)
    integer :: i
    do i=0,nmx-1
       c_analytic(i)=cmplx(wave(i+1),0.0)
    enddo
    call four2a(c_analytic, nmx, 1, -1, 1)
    c_analytic(nmx/2+1:nmx-1)=cmplx(0.0,0.0)
    do i=1,nmx/2-1
       c_analytic(i)=2.0*c_analytic(i)
    enddo
    call four2a(c_analytic, nmx, 1, 1, 1)
    c_analytic=c_analytic/real(nmx)
  end subroutine make_analytic

  integer function spec_rv_near(cand, nc, f0)
    real, intent(in)    :: cand(3,MAXCAND), f0
    integer, intent(in) :: nc
    integer :: k, kbest
    real    :: sbest
    kbest=-1; sbest=-1.0e30
    do k=1,nc
       if(abs(cand(1,k)-f0).le.20.0 .and. cand(2,k).gt.sbest) then
          sbest=cand(2,k); kbest=k
       endif
    enddo
    if(kbest.lt.0) then
       spec_rv_near=-1
    else
       spec_rv_near=nint(cand(3,kbest))
    endif
  end function spec_rv_near

end program ft1_rv_detect_test
