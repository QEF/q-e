!
! Copyright (C) 2001-2004 PWSCF group
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
  SELECT CASE( TRIM( psd ) )
     !
     ! ... transition metals
     !
     CASE( 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu' )
        !
        hubbard_l = 2  
        !
     !
     ! ... rare earths
     !
     CASE( 'Ce' )
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
     CASE( 'C', 'O' )
        !
        hubbard_l =  1
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
END FUNCTION set_Hubbard_l
