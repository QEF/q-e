!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION hubbard_occ ( psd )
  !-----------------------------------------------------------------------
  !
  ! This routine is a table (far from being complete) for the total number
  ! of localized electrons in transition metals or rare earths
  ! (PPs usually are built on non physical configurations)
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=2), INTENT(IN) :: psd
  REAL(DP)                     :: hubbard_occ
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
     !
     ! TRANSITION METALS
     !
     CASE( 'Ti', 'Zr', 'Hf' )
        hubbard_occ = 2.d0
        !
     CASE( 'V', 'Nb', 'Ta' )
        hubbard_occ = 3.d0
     !
     CASE( 'Cr', 'Mo', 'W' )
        hubbard_occ = 5.d0
     !
     CASE( 'Mn', 'Tc', 'Re' )
        hubbard_occ = 5.d0
     !
     CASE( 'Fe', 'Ru', 'Os' )
        hubbard_occ = 6.d0
     !
     CASE( 'Co', 'Rh', 'Ir' )
        hubbard_occ = 7.d0
     !
     CASE( 'Ni', 'Pd', 'Pt' )
        hubbard_occ = 8.d0
     !
     CASE( 'Cu', 'Ag', 'Au' )
        hubbard_occ = 10.d0
     !
     CASE( 'Zn', 'Cd', 'Hg' )
        hubbard_occ = 10.d0
     !
     ! RARE EARTHS
     !
     CASE( 'Ce', 'Th' )
        hubbard_occ = 2.d0
     !
     CASE( 'Pr', 'Pa' )
        hubbard_occ = 3.d0
     !
     CASE( 'Nd', 'U'  )
        hubbard_occ = 4.d0
     !
     CASE( 'Pm', 'Np' )
        hubbard_occ = 5.d0
     !
     CASE( 'Sm', 'Pu' )
        hubbard_occ = 6.d0
     !
     CASE( 'Eu', 'Am' )
        hubbard_occ = 6.d0
     !
     CASE( 'Gd', 'Cm' )
        hubbard_occ = 7.d0
     !
     CASE( 'Tb', 'Bk' )
        hubbard_occ = 8.d0
     !
     CASE( 'Dy', 'Cf' )
        hubbard_occ = 9.d0
     !
     CASE( 'Ho', 'Es' )
        hubbard_occ =10.d0
     !
     CASE( 'Er', 'Fm' )
        hubbard_occ =11.d0
     !
     CASE( 'Tm', 'Md' )
        hubbard_occ =12.d0
     !
     CASE( 'Yb', 'No' )
        hubbard_occ =13.d0
     !
     CASE( 'Lu', 'Lr' )
        hubbard_occ =14.d0
     !
     ! OTHER ELEMENTS
     !
     CASE( 'C'  )
        hubbard_occ = 2.d0
     !
     CASE( 'N'  )
        hubbard_occ = 3.d0
     !
     CASE( 'O'  )
        hubbard_occ = 4.d0
     !
     CASE( 'H'  )
        hubbard_occ = 1.d0
     !
     CASE( 'Ga', 'In'  )
        hubbard_occ = 10.d0
     !
     !
     ! NOT INSERTED
     !
     CASE DEFAULT
        hubbard_occ = 0.d0
        call errore ('hubbard_occ', 'pseudopotential not yet inserted', 1)
     !
  END SELECT
  
  RETURN

END FUNCTION hubbard_occ

