!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine read_config_from_file
  !-----------------------------------------------------------------------

  USE io_global,      ONLY : stdout
  use pwcom
  use io_files,       only: prefix, iunres
  use restart_module, only : readfile_config

  implicit none

  ! parameter indicating from where to restart
  integer :: nat_, ibrav_, ierr
  real(kind=DP) :: alat_, at_(3,3)
  real(kind=DP) :: tau_(3,nat)
  logical exst

  if (trim(startingconfig).ne.'file') return

  WRITE( stdout, '(/5x,"Starting configuration read from file ", a14 )') &
        trim(prefix)//".save"
  !
  !     check if restart file is present, if yes read config parameters
  !
  call readfile_config( iunres, ibrav_, nat_, alat_, at_, tau_, ierr )
  if ( ierr == 1 ) then
     WRITE( stdout, '(/5x,"Failed to open file", a14 )') trim(prefix)//".save"
     WRITE( stdout, '(/5x,"Use input configuration")')
     return
  else if( ierr > 1 ) then
     call errore ('read_config_from_file', 'problems in reading file', 1)
  endif
  !
  !  check if atomic positions from restart file if present
  !
  if (nat_.ne.nat.or.ibrav_.ne.ibrav) then
     WRITE( stdout,*) 'wrong nat ', nat, nat_, ' or ibrav ', ibrav, ibrav_
     call errore('read_config_from_file','wrong nat or ibrav',1)
  endif
  alat = alat_
  at   = at_
  tau(:,1:nat) = tau_(:,1:nat)

  call volume (alat, at(1,1), at(1,2), at(1,3), omega)
  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  if (lmovecell) then
     !
     ! input value of at and omega (currently stored in xxx_old variables)
     ! must be used to initialize G vectors and other things
     ! swap xxx and xxx_old variables and scale the atomic position to the
     ! input cell shape in order to check the symmetry.
     !
     call cryst_to_cart (nat, tau, bg, - 1)
     call swap (9, at, at_old)
     call swap (1, omega, omega_old)
     call cryst_to_cart (nat, tau, at, + 1)
  endif
  !
  return
end subroutine read_config_from_file

!-----------------------------------------------------------------------
subroutine read_config_from_file_old
  !-----------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  use pwcom
  use io_files,   only : prefix

  implicit none

  ! parameter indicating from where to restart
  integer :: nat_, ibrav_, iunit
  logical exst

  if (trim(startingconfig).ne.'file') return

  WRITE( stdout, '(/5x,"Starting configuration read from file ", a14 )') &
        trim(prefix)//".config"
  !
  !     check if restart file is present
  !
  iunit = 1
  call seqopn (iunit, trim(prefix)//".config", 'unformatted', exst)
  if (.not.exst) then
     close (unit = iunit, status = 'delete')
     WRITE( stdout, '(/5x,"Failed to open file", a14 )') trim(prefix)//".config"
     WRITE( stdout, '(/5x,"Use input configuration")')
     return
  endif
  !
  !  read atomic positions from restart file if present
  !
  read (iunit, err = 10, end = 10) ibrav_, nat_

  if (nat_.ne.nat.or.ibrav_.ne.ibrav) then
     WRITE( stdout,*) 'wrong nat ', nat, nat_, ' or ibrav ', ibrav, ibrav_
     call errore('read_config_from_file','wrong nat or ibrav',1)
  endif

  read (iunit, err = 10, end = 10) alat, at, tau
  call volume (alat, at(1,1), at(1,2), at(1,3), omega)
  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  if (lmovecell) then
     !
     ! input value of at and omega (currently stored in xxx_old variables)
     ! must be used to initialize G vectors and other things
     ! swap xxx and xxx_old variables and scale the atomic position to the
     ! input cell shape in order to check the symmetry.
     !
     call cryst_to_cart (nat, tau, bg, - 1)
     call swap (9, at, at_old)
     call swap (1, omega, omega_old)
     call cryst_to_cart (nat, tau, at, + 1)
  endif
  !
  !  close the file for later use
  !
  close (unit = iunit, status = 'keep')
  return

10 call errore ('read_config_from_file', 'problems in reading file', 1)

end subroutine read_config_from_file_old
