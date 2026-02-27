program ft1_freq_est_test

! Test accuracy of ft1_freq_est_acorr frequency estimator.
! Generates a signal at f0, adds AWGN, downsamples correctly, then applies
! a known frequency offset to the baseband signal before calling the
! estimator. Sweeps offset from 0 to df_max, reports estimation statistics.
!
! Usage: ft1_freq_est_test [snrdb ntrials df_max df_step]
!   Defaults: -5.0 500 1.5 0.1

  use packjt77
  use cpm_trellis_mod
  use matched_filter_bank_mod
  include 'ft1_params.f90'

  parameter (NDMAX=NMAX/NDOWN)

  character*37 msg37,msgsent37
  character arg*12
  integer itone(NN)
  integer*1 msgbits(77)
  real wave(NMAX)
  real dd(NMAX)
  real dd_noisy(NMAX)
  complex cd(0:NDMAX-1)
  complex cd_shifted(0:NDMAX-1)
  logical newdata
  real*8 twopi
  real df_est, df_error

  twopi=8.d0*atan(1.d0)
  msg37='CQ W9XYZ EN37'
  f0=1500.0

  ! Parse command-line arguments
  snrdb=-5.0
  ntrials=500
  df_max=1.5
  df_step_arg=0.1
  nargs=iargc()
  if(nargs.ge.1) then
     call getarg(1,arg)
     read(arg,*) snrdb
  endif
  if(nargs.ge.2) then
     call getarg(2,arg)
     read(arg,*) ntrials
  endif
  if(nargs.ge.3) then
     call getarg(3,arg)
     read(arg,*) df_max
  endif
  if(nargs.ge.4) then
     call getarg(4,arg)
     read(arg,*) df_step_arg
  endif

  fs=12000.0
  bandwidth_ratio=2500.0/(fs/2.0)
  dt_samp=real(NDOWN)/fs

  ! Encode the test message
  call genft1(msg37,0,msgsent37,msgbits,itone)
  write(*,'(a,a37)') 'Test message: ',msgsent37
  write(*,'(a,f6.1,a,i5)') 'SNR: ',snrdb,' dB   Trials: ',ntrials
  write(*,'(a,f6.3,a,f6.3,a)') 'Offset sweep: 0 to ',df_max, &
       ' Hz, step ',df_step_arg,' Hz'
  write(*,*)

  ! Generate clean 4-CPM waveform
  nwave=NMAX
  wave=0.
  call gen_ft1wave(itone,NN,NSPS_NUM,NSPS_DEN,fs,f0,wave,nwave)
  dd=0.
  dd(1:NMAX)=wave(1:NMAX)

  sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)
  call sgran()

! === Direct diagnostic: compare cd with MF at symbol 0 ===
  call init_cpm_trellis()
  call init_matched_filters(NSS)

  ! Generate noiseless signal, downsample
  dd=0.
  dd(1:NMAX)=wave(1:NMAX)
  newdata=.true.
  call ft1_downsample(dd,newdata,f0,cd)
  sum2=sum(real(cd*conjg(cd)))/real(NDMAX)
  if(sum2.gt.0.0) cd=cd/sqrt(sum2)

  write(*,'(a)') 'RAW cd vs MF at symbol 0 (state 1, a=0):'
  write(*,'(a)') '  k  cd_re     cd_im     mf_re     mf_im     cd_phase  mf_phase'
  do i=0,NSS-1
     write(*,'(i3,6f10.5)') i, real(cd(i)), aimag(cd(i)), &
          real(mf_bank(i+1,1,0)), aimag(mf_bank(i+1,1,0)), &
          atan2(aimag(cd(i)),real(cd(i))), &
          atan2(aimag(mf_bank(i+1,1,0)),real(mf_bank(i+1,1,0)))
  enddo
  write(*,'(a)') 'RAW cd vs MF at symbol 1 (state 1, a=2):'
  write(*,'(a)') '  k  cd_re     cd_im     mf_re     mf_im     cd_phase  mf_phase'
  do i=0,NSS-1
     write(*,'(i3,6f10.5)') i, real(cd(8+i)), aimag(cd(8+i)), &
          real(mf_bank(i+1,1,2)), aimag(mf_bank(i+1,1,2)), &
          atan2(aimag(cd(8+i)),real(cd(8+i))), &
          atan2(aimag(mf_bank(i+1,1,2)),real(mf_bank(i+1,1,2)))
  enddo
  write(*,*)

  write(*,'(a)') '  df_true(Hz)   mean_err    rms_err  max_abs_err'
  write(*,'(a)') '  -----------   --------    -------  -----------'

  ndf=nint(df_max/df_step_arg)+1
  do idf=0,ndf-1
     df=real(idf)*df_step_arg

     err_sum = 0.0
     err_sq_sum = 0.0
     err_abs_max = 0.0

     do itrial=1,ntrials
        ! Add AWGN
        do i=1,NMAX
           dd_noisy(i)=sig*dd(i) + gran()
        enddo

        ! Downsample at correct frequency
        newdata=.true.
        call ft1_downsample(dd_noisy,newdata,f0,cd)

        ! Normalize
        sum2=sum(real(cd*conjg(cd)))/real(NDMAX)
        if(sum2.gt.0.0) cd=cd/sqrt(sum2)

        ! Apply frequency offset (simulating coarse sync residual error)
        do i=0,NDMAX-1
           t=real(i)*dt_samp
           cd_shifted(i)=cd(i)*cmplx(cos(twopi*df*t),sin(twopi*df*t))
        enddo

        ! Call frequency estimator
        call ft1_freq_est_acorr(cd_shifted, NDMAX, 0, df_est)

        df_error = df_est - df
        err_sum = err_sum + df_error
        err_sq_sum = err_sq_sum + df_error**2
        if(abs(df_error) .gt. err_abs_max) err_abs_max = abs(df_error)
     enddo

     err_mean = err_sum / real(ntrials)
     err_rms = sqrt(err_sq_sum / real(ntrials))
     write(*,'(f10.3,3f12.4)') df, err_mean, err_rms, err_abs_max
  enddo

  write(*,*)
  write(*,'(a)') 'Test complete.'

end program ft1_freq_est_test
