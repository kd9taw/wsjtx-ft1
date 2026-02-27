program ft1_fading_test

! FT1 fading channel test -- exercises turbo decoder under Watterson HF fading.
!
! Pipeline: genft1 -> gen_ft1wave -> watterson -> AWGN -> ft1_downsample
!           -> turbo_decode_ft1
!
! Usage: ft1_fading_test [snr_start snr_end ntrials fspread delay]
!   Defaults: 0 -20 100 1.0 1.0
!
! fspread: Doppler spread (Hz), delay: multipath delay (ms)
! Set both to 0 for AWGN baseline.

  use packjt77
  include 'ft1_params.f90'

  parameter (NDMAX=NMAX/NDOWN)

  character*37 msg37,msgsent37
  character arg*12
  integer itone(NN)
  integer*1 msgbits(77)
  integer*1 message91(91)
  real wave(NMAX)
  real dd(NMAX)
  real dd_faded(NMAX)
  real dd_noisy(NMAX)
  real llr_out(174)
  complex cd(0:NDMAX-1)
  complex c0(0:NMAX-1)
  complex c_analytic(0:NMAX-1)
  logical newdata

  msg37='CQ W9XYZ EN37'
  f0=1500.0

  ! Parse command-line arguments
  snr_start=0.0
  snr_end=-20.0
  ntrials=100
  fspread=1.0
  delay=1.0
  nargs=iargc()
  if(nargs.ge.1) then
     call getarg(1,arg)
     read(arg,*) snr_start
  endif
  if(nargs.ge.2) then
     call getarg(2,arg)
     read(arg,*) snr_end
  endif
  if(nargs.ge.3) then
     call getarg(3,arg)
     read(arg,*) ntrials
  endif
  if(nargs.ge.4) then
     call getarg(4,arg)
     read(arg,*) fspread
  endif
  if(nargs.ge.5) then
     call getarg(5,arg)
     read(arg,*) delay
  endif

  fs=12000.0
  bandwidth_ratio=2500.0/(fs/2.0)

  ! Encode the test message
  call genft1(msg37,0,msgsent37,msgbits,itone)
  write(*,'(a,a37)') 'Test message: ',msgsent37
  write(*,'(a,f7.1,a)') 'Carrier freq: ',f0,' Hz'
  write(*,'(a,f6.2,a,f5.2,a)') 'Fading:       fspread=',fspread, &
       ' Hz, delay=',delay,' ms'
  write(*,'(a,i4,a)') 'Trials:       ',ntrials,' per SNR point'
  write(*,*)

  ! Generate clean 4-CPM waveform at 12 kHz
  nwave=NMAX
  wave=0.
  call gen_ft1wave(itone,NN,NSPS_NUM,NSPS_DEN,fs,f0,wave,nwave)

  dd=0.
  dd(1:NMAX)=wave(1:NMAX)

  ! Compute original RMS over signal portion for SNR calibration
  rms_orig=sqrt(sum(dd(1:NZ)**2)/NZ)

  ! Pre-compute analytic signal (Hilbert transform) for correct fading
  ! Real passband wrapped as cmplx(real,0) has energy at +f0 and -f0.
  ! Zeroing negative frequencies gives the analytic signal, so watterson
  ! applies Rayleigh fading (not Gaussian) to the envelope.
  do i=0,NMAX-1
     c_analytic(i)=cmplx(dd(i+1), 0.0)
  enddo
  call four2a(c_analytic, NMAX, 1, -1, 1)    !Forward FFT
  c_analytic(NMAX/2+1:NMAX-1) = cmplx(0.0, 0.0)  !Zero negative freqs
  do i=1,NMAX/2-1
     c_analytic(i) = 2.0 * c_analytic(i)      !Double positive freqs
  enddo
  ! DC (i=0) and Nyquist (i=NMAX/2) unchanged
  call four2a(c_analytic, NMAX, 1, 1, 1)     !Inverse FFT
  c_analytic = c_analytic / real(NMAX)         !Normalize FFT

  call sgran()

  ! SNR sweep
  nsnr=nint(snr_start-snr_end)+1
  if(nsnr.gt.50) nsnr=50

  write(*,'(a)') '   SNR   Decoded  Trials   Rate    Avg_nerr'
  write(*,'(a)') '   ---   -------  ------   ----    --------'

  do isnr=1,nsnr
     snrdb=snr_start - real(isnr-1)
     sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)

     ndec=0
     nerr_total=0

     do itrial=1,ntrials
        ! Apply Watterson fading (new realization each trial)
        if(fspread.ne.0.0 .or. delay.ne.0.0) then
           c0 = c_analytic
           call watterson(c0, NMAX, NZ, fs, delay, fspread)
           do i=1,NMAX
              dd_faded(i)=real(c0(i-1))
           enddo
           ! Renormalize to match original avg power for correct SNR calibration
           rms_faded=sqrt(sum(dd_faded(1:NZ)**2)/NZ)
           if(rms_faded.gt.0.0) dd_faded=dd_faded*(rms_orig/rms_faded)
        else
           dd_faded=dd
        endif

        ! Add AWGN
        do i=1,NMAX
           dd_noisy(i)=sig*dd_faded(i) + gran()
        enddo

        ! Downsample to baseband complex signal
        newdata=.true.
        call ft1_downsample(dd_noisy,newdata,f0,cd)

        ! Normalize signal to unit power per sample
        sum2=sum(real(cd*conjg(cd)))/real(NDMAX)
        if(sum2.gt.0.0) cd=cd/sqrt(sum2)

        ! Run turbo decoder
        npts=NDMAX
        dt0=0.0
        snr_est=snrdb
        ntype=-1
        nharderror=-1
        dmin=0.0
        message91=0
        llr_out=0.0

        call turbo_decode_ft1(cd,npts,f0,dt0,snr_est,llr_out, &
             message91,ntype,nharderror,dmin,0,ncheck_out)

        if(ntype.ge.0) then
           ndec=ndec+1
           if(nharderror.ge.0) nerr_total=nerr_total+nharderror
        endif
     enddo

     rate=real(ndec)/real(ntrials)
     if(ndec.gt.0) then
        avg_nerr=real(nerr_total)/real(ndec)
     else
        avg_nerr=-1.0
     endif
     write(*,'(f6.1,i9,i8,f8.3,f10.1)') snrdb,ndec,ntrials,rate,avg_nerr
  enddo

  write(*,*)
  write(*,'(a)') 'Test complete.'

end program ft1_fading_test
