program ft1_freq_tolerance_test

! Test turbo decoder frequency tolerance.
! Generates a signal at f0, downsamples correctly, then applies a frequency
! offset df to the baseband signal before calling turbo_decode_ft1.
! Sweeps df from 0 to df_max in steps of df_step.
! Reports decode success rate at each offset.
!
! Usage: ft1_freq_tolerance_test [snrdb ntrials df_max df_step niter_max]
!   Defaults: -5.0 100 1.0 0.05 0

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
  real dd_noisy(NMAX)
  real llr_out(174)
  complex cd(0:NDMAX-1)
  complex cd_shifted(0:NDMAX-1)
  logical newdata
  real*8 twopi

  twopi=8.d0*atan(1.d0)
  msg37='CQ W9XYZ EN37'
  f0=1500.0

  ! Parse command-line arguments
  snrdb=-5.0
  ntrials=100
  df_max=1.0
  df_step=0.05
  niter_max_arg=0
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
     read(arg,*) df_step
  endif
  if(nargs.ge.5) then
     call getarg(5,arg)
     read(arg,*) niter_max_arg
  endif

  fs=12000.0
  bandwidth_ratio=2500.0/(fs/2.0)
  dt_samp=real(NDOWN)/fs

  ! Encode the test message
  call genft1(msg37,0,msgsent37,msgbits,itone)
  write(*,'(a,a37)') 'Test message: ',msgsent37
  write(*,'(a,f6.1,a,i5,a)') 'SNR: ',snrdb,' dB   Trials: ',ntrials,' per offset'
  write(*,'(a,f6.3,a,f6.3,a)') 'Offset sweep: 0 to ',df_max,' Hz, step ',df_step,' Hz'
  write(*,'(a,i3,a)') 'niter_max: ',niter_max_arg, &
       ' (0=full turbo, N=probe mode returning ncheck after N iter)'
  write(*,*)

  ! Generate clean 4-CPM waveform
  nwave=NMAX
  wave=0.
  call gen_ft1wave(itone,NN,NSPS_NUM,NSPS_DEN,fs,f0,wave,nwave)
  dd=0.
  dd(1:NMAX)=wave(1:NMAX)

  sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)
  call sgran()

  if(niter_max_arg.eq.0) then
     write(*,'(a)') '  df(Hz)   Decoded  Trials   Rate'
     write(*,'(a)') '  ------   -------  ------   ----'
  else
     write(*,'(a)') '  df(Hz)   Decoded  Trials   Rate    Avg_ncheck'
     write(*,'(a)') '  ------   -------  ------   ----    ----------'
  endif

  ndf=nint(df_max/df_step)+1
  do idf=0,ndf-1
     df=real(idf)*df_step

     ndec=0
     ncheck_total=0
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

        ! Apply frequency offset (simulating frequency estimation error)
        do i=0,NDMAX-1
           t=real(i)*dt_samp
           cd_shifted(i)=cd(i)*cmplx(cos(twopi*df*t),sin(twopi*df*t))
        enddo

        ! Run turbo decoder
        ntype=-1
        nharderror=-1
        dmin=0.0
        message91=0
        llr_out=0.0
        ncheck_out=83
        call turbo_decode_ft1(cd_shifted,NDMAX,0.0,0.0,snrdb, &
             llr_out,message91,ntype,nharderror,dmin, &
             niter_max_arg,ncheck_out)

        if(ntype.ge.0) ndec=ndec+1
        ncheck_total=ncheck_total+ncheck_out
     enddo

     rate=real(ndec)/real(ntrials)
     avg_ncheck=real(ncheck_total)/real(ntrials)
     if(niter_max_arg.eq.0) then
        write(*,'(f8.3,i9,i8,f8.3)') df,ndec,ntrials,rate
     else
        write(*,'(f8.3,i9,i8,f8.3,f12.1)') df,ndec,ntrials,rate,avg_ncheck
     endif
     if(rate.lt.0.01 .and. df.gt.0.2 .and. niter_max_arg.eq.0) exit
  enddo

  write(*,*)
  write(*,'(a)') 'Test complete.'

end program ft1_freq_tolerance_test
