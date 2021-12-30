!
! Copyright (C) 2004-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE uspp_param
  !
  ! ... Ultrasoft and Norm-Conserving pseudopotential parameters
  !  
  USE pseudo_types, ONLY : pseudo_upf
  IMPLICIT NONE
  SAVE
  !
  INTEGER :: nsp = 0 
  TYPE (pseudo_upf),  ALLOCATABLE, TARGET :: upf(:)
  !! the upf structure contains all info on atomic pseudopotential parameters
  INTEGER, ALLOCATABLE :: nh(:)
  !! number of beta functions, with angular parts, per atomic type 
  INTEGER :: nhm
  !! max number of beta functions, including angular parts, across atoms
  INTEGER :: nbetam
  !! max number of radial beta functions
  INTEGER :: nwfcm
  !! max number of radial atomic wavefunctions across atoms
  INTEGER :: lmaxkb
  !! max angular momentum of beta functions
  INTEGER :: lmaxq
  !! max angular momentum + 1 for Q functions
  !
CONTAINS
  !
  SUBROUTINE init_uspp_dims ()
    !
    !!     calculates the number of beta functions for each atomic type
    !
    IMPLICIT NONE
    !
    INTEGER :: nt, nb
    !
    ! Check is needed, may be called more than once (but it shouldn't be!)
    ! Maybe nh should be allocated when upf is, when upf is read ?
    !
    IF ( .NOT. ALLOCATED(nh) ) ALLOCATE ( nh(nsp) )
    !
    lmaxkb = - 1
    DO nt = 1, nsp
       !
       nh (nt) = 0
       !
       ! do not add any beta projector if pseudo in 1/r fmt (AF)
       !
       IF ( upf(nt)%tcoulombp ) CYCLE 
       !
       DO nb = 1, upf(nt)%nbeta
          nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
          lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
       ENDDO
       !
    ENDDO
    lmaxq = 2*lmaxkb+1
    !
    ! calculate max numbers of beta functions and of atomic wavefunctions
    !
    nhm    = MAXVAL (nh (1:nsp))
    nbetam = MAXVAL (upf(1:nsp)%nbeta)
    nwfcm  = MAXVAL (upf(1:nsp)%nwfc)
    !
  END SUBROUTINE init_uspp_dims
  !
END MODULE uspp_param
