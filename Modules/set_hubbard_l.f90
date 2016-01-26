!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
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
     CASE( 'C', 'N', 'O' )
        !
        hubbard_l =  1
        !
     CASE( 'As', 'Ga', 'In' )
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
END FUNCTION set_Hubbard_l
