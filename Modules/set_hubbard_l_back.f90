!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
     ! ... transition metals
     !

     CASE( 'H', 'He', 'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba',&
           'Fr', 'Ra' )
        !
        hubbard_l_back =  -1 ! no background states
        !
     CASE( 'Se' )
!        !
        hubbard_l_back =  2
        !
     CASE( 'Zn' ) 
        !
        hubbard_l_back =  1
        !
     CASE DEFAULT
        !
        hubbard_l_back = 0 ! for almost all the elements the background states are s states
        CALL infomsg( 'set_hubbard_l_back', 'The default background state for DFT+U is s orbitals' )
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_hubbard_l_back
