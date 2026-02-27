program ft1sim

! Generate simulated FT1 signals with AWGN for decoder testing.
!
! Usage: ft1sim "message" f0 DT fdop del nfiles snr
!
! Generates .wav files containing FT1 4-CPM signals at the specified
! SNR (in 2500 Hz bandwidth). Optionally applies Watterson fading.
!
! Template: ft4sim.f90

  use wavhdr
  use packjt77
  include 'ft1_params.f90'
  type(hdr) h
  character arg*12,fname*17
  character msg37*37,msgsent37*37
  character c77*77
  real wave(4*NMAX)                      !Waveform at 48kHz (4x oversample)
  integer itone(NN)
  integer*1 msgbits(77)
  integer*2 iwave(NMAX)                  !Output waveform at 12kHz
  real dd(NMAX)                          !Float working buffer

! Get command-line arguments
  nargs=iargc()
  if(nargs.ne.7) then
     print*,'Usage:    ft1sim  "message"          f0    DT fdop del nfiles snr'
     print*,'Examples: ft1sim  "CQ W9XYZ EN37"   1500  0.0  0.1 1.0   10  -15'
     go to 999
  endif
  call getarg(1,msg37)
  call getarg(2,arg)
  read(arg,*) f0
  call getarg(3,arg)
  read(arg,*) xdt
  call getarg(4,arg)
  read(arg,*) fspread
  call getarg(5,arg)
  read(arg,*) delay
  call getarg(6,arg)
  read(arg,*) nfiles
  call getarg(7,arg)
  read(arg,*) snrdb

  nfiles=abs(nfiles)
  twopi=8.0*atan(1.0)
  fs=12000.0
  dt=1.0/fs
  baud=12000.0/real(NSPS)                !~28 Bd
  txt=real(NZ)/12000.0                   !TX duration (s)

  bandwidth_ratio=2500.0/(fs/2.0)
  sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)
  if(snrdb.gt.90.0) sig=1.0

! Encode message -> 99 quaternary channel symbols
  call genft1(msg37,0,msgsent37,msgbits,itone)
  write(*,*)
  write(*,'(a9,a37)') 'Message: ',msgsent37
  write(*,1000) f0,xdt,txt,snrdb
1000 format('f0:',f9.3,'   DT:',f6.2,'   TxT:',f6.1,'   SNR:',f6.1)
  write(*,*)
  write(*,'(a17)') 'Channel symbols: '
  write(*,'(99i1)') itone
  write(*,*)

  call sgran()

! Generate 4-CPM waveform at 12 kHz sample rate
! gen_ft1wave expects: itone, nsym, nsps_num, nsps_den, fsample, f0, wave, nwave
  nwave=NMAX
  wave=0.
  call gen_ft1wave(itone,NN,NSPS_NUM,NSPS_DEN,fs,f0,wave(1:NMAX),nwave)

! Copy to float working buffer
  dd=wave(1:NMAX)

! Apply time offset
  k=nint(xdt/dt)
  dd=cshift(dd,-k)
  if(k.gt.0) dd(1:k)=0.0
  if(k.lt.0) dd(NMAX+k+1:NMAX)=0.0

  do ifile=1,nfiles
     wave(1:NMAX)=sig*dd

     if(snrdb.lt.90) then
        do i=1,NMAX
           xnoise=gran()
           wave(i)=wave(i) + xnoise
        enddo
     endif

     gain=100.0
     if(snrdb.lt.90.0) then
       wave(1:NMAX)=gain*wave(1:NMAX)
     else
       datpk=maxval(abs(wave(1:NMAX)))
       fac=32766.9/datpk
       wave(1:NMAX)=fac*wave(1:NMAX)
     endif
     if(any(abs(wave(1:NMAX)).gt.32767.0)) print*,"Warning - data will be clipped."
     iwave=nint(wave(1:NMAX))
     h=default_header(12000,NMAX)
     write(fname,1102) ifile
1102 format('000000_',i6.6,'.wav')
     open(10,file=fname,status='unknown',access='stream')
     write(10) h,iwave
     close(10)
     write(*,1110) ifile,xdt,f0,snrdb,fname
1110 format(i4,f7.2,f8.2,f7.1,2x,a17)
  enddo

999 end program ft1sim
