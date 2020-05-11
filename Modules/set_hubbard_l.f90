!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
FUNCTION set_hubbard_l( psd ) RESULT( hubbard_l )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_l
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
     !
     ! ... transition metals
     !
     CASE( 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
           'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
           'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'  )
        !
        hubbard_l = 2  
        !
     !
     ! ... rare earths
     !
     CASE('Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', & 
          'Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' )
        !
        hubbard_l = 3
        !
     !
     ! ... other elements
     !
     CASE( 'H' )
        !
        hubbard_l =  0
        !
     CASE( 'C', 'N', 'O', 'As' )
        !
        hubbard_l =  1
        !
     CASE( 'Ga', 'In' )
        !
        hubbard_l =  2
        !
     CASE DEFAULT
        !
        hubbard_l = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_l', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_hubbard_l
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION set_hubbard_l_back( psd ) RESULT( hubbard_l_back )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_l_back
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
     !
     CASE( 'H', 'He', 'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba',&
           'Fr', 'Ra' )
        !
        hubbard_l_back =  -1 ! no background states
        !
     CASE( 'Se' )
        !
        hubbard_l_back =  2
        !
     CASE( 'Zn' ) 
        !
        hubbard_l_back =  1
        !
     CASE DEFAULT
        !
        ! For almost all the elements the background states are s states
        !
        hubbard_l_back = 0
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_hubbard_l_back
