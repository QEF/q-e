!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine tabd (nt, occ_loc)  
  !-----------------------------------------------------------------------
  !
  ! This routine is a table (far from being complete) for the total number
  ! of localized electrons in transition metals or rare earths
  ! (PPs usually are built on non physical configurations)
  !
  USE kinds, ONLY: DP
  USE uspp_param, ONLY: upf
  implicit none
  real(DP) :: occ_loc
  ! output: the total number of d electrons

  integer :: nt
  
  SELECT CASE( TRIM(ADJUSTL(upf(nt)%psd)) )
     !
     ! TRANSITION METALS
     !
     CASE( 'Mn' )
        occ_loc = 5.d0
     !
     CASE( 'Fe' )
        occ_loc = 6.d0
     !
     CASE( 'Co' )
        occ_loc = 7.d0
     !
     CASE( 'Ni' )
        occ_loc = 8.d0
     !
     CASE( 'Cu' )
        occ_loc = 10.d0
     !
     ! RARE EARTHS
     !
     CASE( 'Ce' )
        occ_loc = 2.d0
     !
     ! OTHER ELEMENTS
     !
     CASE( 'C'  )
        occ_loc = 2.d0
     !
     CASE( 'O'  )
        occ_loc = 4.d0
     !
     CASE( 'H'  )
        occ_loc = 1.d0
     !
     ! NOT INSERTED
     !
     CASE DEFAULT
        occ_loc = 0.d0
        call errore ('tabd', 'pseudopotential not yet inserted', 1)
     !
  END SELECT
  
  return
end subroutine tabd

