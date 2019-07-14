!
! Copyright (C) 2019 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE set_occupations( occupations, smearing, degauss, &
     lfixed, ltetra, tetra_type,  smearing_, lgauss, ngauss ) 
  !------------------------------------------------------------------------
  USE kinds, ONLY: dp
  ! 
  IMPLICIT NONE 
  !
  CHARACTER(LEN=*), INTENT(IN)    :: occupations
  CHARACTER(LEN=*), INTENT(IN) :: smearing
  CHARACTER(LEN=*), INTENT(OUT):: smearing_
  REAL(dp), INTENT(INOUT) :: degauss
  LOGICAL, INTENT(OUT) :: lfixed, lgauss, ltetra
  INTEGER, INTENT(OUT) :: tetra_type, ngauss
  !
  lfixed = .FALSE.
  ltetra = .FALSE.
  tetra_type = 0
  smearing_ = 'none'
  lgauss = .FALSE. 
  ngauss = 0
  
  SELECT CASE( trim( occupations ) )
  CASE( 'fixed' )
     !
     IF ( degauss /= 0.D0 ) THEN
        CALL errore( ' set_occupations ', &
             & ' fixed occupations, gauss. broadening ignored', -1 )
        degauss = 0.D0
     ENDIF
     !
  CASE( 'smearing' )
     !
     lgauss = ( degauss > 0.0_dp ) 
     IF ( .NOT. lgauss ) CALL errore( ' set_occupations ',  &
          ' smearing requires gaussian broadening', 1 )
     !
     SELECT CASE ( trim( smearing ) )
     CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
        ngauss = 0
        smearing_ = 'Gaussian'
     CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
        ngauss = 1
        smearing_ = 'Methfessel-Paxton'
     CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
        ngauss = -1
        smearing_ = 'Marzari-Vanderbilt'
     CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
        ngauss = -99
        smearing_ = 'Fermi-Dirac'
     CASE DEFAULT
        CALL errore( ' set_occupations ', &
             ' smearing '//trim(smearing)//' unknown', 1 )
     END SELECT
     !
  CASE( 'tetrahedra' )
     !
     ltetra = .true.
     tetra_type = 0
     !
  CASE( 'tetrahedra_lin', 'tetrahedra-lin')
     !
     ltetra = .true.
     tetra_type = 1
     !
  CASE('tetrahedra_opt', 'tetrahedra-opt')
     !
     ltetra = .true.
     tetra_type = 2
     !
  CASE( 'from_input' )
     !
     ngauss = 0
     lfixed = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'set_occupations', &
          'occupations ' // trim( occupations ) // ' not implemented', 1 )
     !
  END SELECT
  !
END SUBROUTINE set_occupations
