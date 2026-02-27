! FT1
! 4-CPM h=1/2 BT=0.3, LDPC(174,91) code, three 4x4 Costas arrays for sync
! Turbo equalization: iterative BCJR + LDPC BP decoding

parameter (KK=91)                     !Information bits (77 + CRC14)
parameter (ND=87)                     !Data symbols (174 coded bits / 2 bits per symbol)
parameter (NS=12)                     !Sync symbols (3 groups of 4, Costas 4x4)
parameter (NN=NS+ND)                  !Total channel symbols (99)
parameter (NSPS_NUM=3000,NSPS_DEN=7)  !Samples per symbol = 3000/7 = 428.571...
parameter (NSPS=429)                  !Nominal samples per symbol (rounded up)
parameter (NZ=42429)                  !Total TX samples (99 * 3000/7, rounded)
parameter (NMAX=4*12000)              !Samples in iwave (4.0s * 12000 Hz = 48000)
parameter (NFFT1=1024, NH1=NFFT1/2)  !Length of FFTs for spectrogram
parameter (NSTEP=107)                 !Spectrogram step size (~quarter symbol)
parameter (NHSYM=(NMAX-NFFT1)/NSTEP) !Number of spectral columns
parameter (NDOWN=54)                  !Downsample factor (12000/54 = 222.2 Hz, ~8 samp/sym)
parameter (NSS=8)                     !Downsampled samples/symbol: ceil(NSPS/NDOWN)=ceil(7.94)
parameter (MAXCAND=200)               !Maximum sync candidates
