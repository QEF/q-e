!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE parameters
  
  IMPLICIT NONE
  SAVE

  !
  !       First all the parameter declaration
  !

  INTEGER, PARAMETER :: &
       ntypx  = 10,     &! max number of different types of atom
       npsx   = ntypx,  &! max number of different PPs (obsolete)
       npk    = 40000,  &! max number of k-points               
       lmaxx  = 3,      &! max non local angular momentum (l=0 to lmaxx)      
       cp_lmax  = lmaxx+1,&! max number of channels !TMP FOR PAW: REMOVE
       nchix  = 6,      & ! max number of atomic wavefunctions per atom
       nwfsx  = 14       ! max number of beta functions per atom
  !
  INTEGER, PARAMETER :: &
       nbrx = 14,          &! max number of beta functions
       lqmax= 2*lmaxx+1,   &! max number of angular momenta of Q
       nqfx = 8             ! max number of coefficients in Q smoothing

  !
  ! ...   More parameter for the CP codes
  !
  INTEGER, PARAMETER :: nacx    = 10         ! max number of averaged 
                                             ! quantities saved to the restart
  INTEGER, PARAMETER :: nsx     = ntypx      ! maximum number of species
  INTEGER, PARAMETER :: natx    = 5000       ! maximum number of atoms
  INTEGER, PARAMETER :: nbndxx  = 10000      ! maximum number of electronic states
  INTEGER, PARAMETER :: npkx    = npk        ! maximum number of K points
  INTEGER, PARAMETER :: ncnsx   = 101        ! maximum number of constraints
  INTEGER, PARAMETER :: nspinx  = 2          ! maximum number of spinors

  INTEGER, PARAMETER :: nhclm   = 4          ! maximum number NH chain length, 
                                             ! nhclm can be easily increased since the restart 
                                             ! file should be able to handle it, perhaps better 
                                             ! to align nhclm by 4
  INTEGER, PARAMETER :: max_nconstr = 100    ! max number of constrains
END MODULE parameters
