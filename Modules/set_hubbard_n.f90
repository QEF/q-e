!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
FUNCTION set_hubbard_n( psd ) RESULT( hubbard_n )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_n
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
     !
     ! ... transition metals, 4-th row 
     !
     CASE( 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')  
         hubbard_n=3
     !
     !  ... transition metals, 5-th  row
     !
     CASE( 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd') 
         hubbard_n=4
     ! 
     ! ... transition metals, 6-th  row 
     !
     CASE( 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'  )
     !
        hubbard_n = 5  
     !
     !
     ! ... rare earths (lanthanoid) 
     !
     CASE('Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu' ) 
     !  
        hubbard_n = 4 
     !  ... rare earths (actinoids )
     CASE ('Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' )
        !
        hubbard_n = 5
        !
     !
     ! ... other elements
     !
     CASE( 'H' )
        !
        hubbard_n =  1
        !
     CASE( 'C', 'N', 'O' )
        !
        hubbard_n =  2
        !
     CASE( 'Ga' ) 
        !
        hubbard_n =  3
        !
     CASE ( 'In', 'As' )
        ! 
        hubbard_n = 4  
        ! 
     CASE DEFAULT
        !
        hubbard_n = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_n', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_hubbard_n
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION set_hubbard_n_back( psd ) RESULT( hubbard_n_back )
  !---------------------------------------------------------------------------
  !
  ! IT: Note, currently this routine is not used anywhere. The data reported 
  !     here is not complete.
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_n_back
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
        hubbard_n_back =  -1 ! no background states
        !
     CASE( 'Se' )
        !
        hubbard_n_back =  3
        !
     CASE( 'Zn' ) 
        !
        hubbard_n_back =  3
        !
     CASE DEFAULT
        !
        hubbard_n_back = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_n_back', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_hubbard_n_back
