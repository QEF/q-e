!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine davcio (vect, nword, unit, nrec, io)
  !----------------------------------------------------------------------
  !
  !   direct-access vector input/output
  !   read/write nword words starting from the address specified by vect
  !
  USE io_global, ONLY : stdout
  use parameters
implicit none
  integer :: nword, unit, nrec, io
  ! input: the dimension of vect
  ! input: the unit where to read/write
  ! input: the record where to read/write
  ! input: flag if < 0 reading if > 0 writing
  real(kind=DP) :: vect (nword)
  ! inp/out: the vector to read/write
  !
  !   the local variables
  !


  integer :: ios
  ! integer variable for I/O control
  call start_clock ('davcio')
  if (unit.le.0) call errore ('davcio', 'wrong unit', 1)
  if (nrec.le.0) call errore ('davcio', 'wrong record number', 2)
  if (nword.le.0) call errore ('davcio', 'wrong record length', 3)
  ios = 0
  if (io.lt.0) then
     read (unit, rec = nrec, iostat = ios) vect
  elseif (io.gt.0) then
     write (unit, rec = nrec, iostat = ios) vect
  else
     call errore ('davcio', 'nothing to do?', - 1)
  endif
  if (ios.ne.0) then
     WRITE( stdout, * ) ' IOS = ', ios
     call errore ('davcio', 'i/o error in davcio', unit)
  endif
  call stop_clock ('davcio')
  return
end subroutine davcio
