!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
module parameters
  !
  !
     IMPLICIT NONE
     SAVE

  !
  !       First all the parameter declaration
  !
  integer , parameter ::                                                 &
       ntypx  = 6,     &! max number of different types of atom
       npsx   = ntypx, &! max number of different PPs (obsolete)
       npk    = 40000, &! max number of k-points               
       lmaxx  = 3,     &! max non local angular momentum       
       nchix  = 6,     &! max number of atomic wavefunctions per atom
       ndm    = 2000    ! max number of points in the atomic radial mesh
  integer , parameter :: DP = kind(0.0d0)

  integer, parameter  :: &
       nbrx = 8,           &! max number of beta functions
       lqmax= 2*lmaxx+1,   &! max number of angular momenta of Q
       nqfm = 8             ! max number of coefficients in Q smoothing

  !
  ! ...   More parameter for the CP codes
  !
  INTEGER, PARAMETER :: mmaxx = 1301    ! maximum mesh size for pseudo
  INTEGER, PARAMETER :: cp_lmax = 4     ! maximum number of channels
                                        ! (s,p,d,f)
  INTEGER, PARAMETER :: nacx  = 10      ! maximum number of averaged 
                                        ! quantities saved to the restart
  INTEGER, PARAMETER :: nsx  = 13       ! maximum number of species
  INTEGER, PARAMETER :: natx  = 599     ! maximum number of atoms
  INTEGER, PARAMETER :: nbndxx = 1000   ! maximum number of electronic states
  INTEGER, PARAMETER :: npkx = 300      ! maximum number of K points
  INTEGER, PARAMETER :: ncnsx = 101     ! maximum number of constraints
  INTEGER, PARAMETER :: nspinx = 2      ! maximum number of spinors
  !
end module parameters

