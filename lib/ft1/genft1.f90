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
subroutine genft1(msg0,ichk,msgsent,msgbits,i4tone)

! Encode an FT1 message, producing 99 quaternary channel symbols.
!
! Input:
!   msg0     - requested message to be transmitted
!   ichk     - if ichk=1, return only msgsent (no encoding)
!
! Output:
!   msgsent  - message as it will be decoded
!   msgbits  - 77-bit message array
!   i4tone   - array of 99 channel symbols, values in {0,1,2,3}
!
! Frame structure (99 symbols):
!   S4 + D43 + S4 + D44 + S4
!   G1(0-3) + data(4-46) + G2(47-50) + data(51-94) + G3(95-98)
!
! Transmit chain:
!   77-bit message -> CRC-14 -> LDPC(174,91) encode ->
!   S-random interleave -> Gray map -> insert Costas sync

  use packjt77
  include 'ft1_params.f90'

  character*37 msg0
  character*37 message                    !Message to be generated
  character*37 msgsent                    !Message as it will be received
  character*77 c77
  integer*4 i4tone(NN),itmp(ND)
  integer*1 codeword(2*ND)               !174 coded bits
  integer*1 interleaved(2*ND)            !174 interleaved bits
  integer*1 msgbits(77)
  integer icos_rv0(4),icos_rv1(4),icos_rv2(4)
  integer icos_ft1(4)
  logical unpk77_success

! Costas arrays for the 3 redundancy versions
  data icos_rv0/0,2,3,1/                 !RV0 (initial transmission)
  data icos_rv1/1,3,2,0/                 !RV1 (first retransmission)
  data icos_rv2/3,0,2,1/                 !RV2 (second retransmission)

  message=msg0

  do i=1,37
     if(ichar(message(i:i)).eq.0) then
        message(i:37)=' '
        exit
     endif
  enddo
  do i=1,37                               !Strip leading blanks
     if(message(1:1).ne.' ') exit
     message=message(i+1:)
  enddo

  i3=-1
  n3=-1
  call pack77(message,i3,n3,c77)
  call unpack77(c77,0,msgsent,unpk77_success) !Unpack to get msgsent

  if(ichk.eq.1) go to 999
  read(c77,'(77i1)',err=1) msgbits
  if(unpk77_success) go to 2
1 msgbits=0
  i4tone=0
  msgsent='*** bad message ***                  '
  go to 999

entry get_ft1_tones_from_77bits(msgbits,i4tone)

! LDPC encode: 77 msg bits -> 91 msg+CRC bits -> 174 coded bits
! (encode174_91 adds CRC-14 internally)
2 call encode174_91(msgbits,codeword)

! S-random interleave (length 174, S=9)
! FT1 does NOT use XOR scrambling (no rvec), unlike FT4
  call ft1_interleave(codeword,interleaved,1)

! Gray code mapping: 2 bits -> 1 quaternary symbol
! bits   tone
!  00     0
!  01     1
!  11     2
!  10     3
  do i=1,ND
     is=interleaved(2*i)+2*interleaved(2*i-1)
     if(is.le.1) itmp(i)=is
     if(is.eq.2) itmp(i)=3
     if(is.eq.3) itmp(i)=2
  enddo

! Insert Costas sync arrays and data symbols
! Default: use RV0 Costas array for initial transmission
  icos_ft1=icos_rv0

! Sync group G1: positions 1-4 (Fortran 1-indexed)
  i4tone(1:4)=icos_ft1

! Data block 1: positions 5-47 (43 data symbols)
  i4tone(5:47)=itmp(1:43)

! Sync group G2: positions 48-51
  i4tone(48:51)=icos_ft1

! Data block 2: positions 52-95 (44 data symbols)
  i4tone(52:95)=itmp(44:87)

! Sync group G3: positions 96-99
  i4tone(96:99)=icos_ft1

999 return
end subroutine genft1
