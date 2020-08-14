!
! Copyright (C) 2019 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE set_occupations( occupations, smearing, degauss, &
     lfixed, ltetra, tetra_type, lgauss, ngauss ) 
  !------------------------------------------------------------------------
  ! Set variables describing occupations
  !
  USE kinds, ONLY: dp
  ! 
  IMPLICIT NONE 
  !
  CHARACTER(LEN=*), INTENT(IN)    :: occupations
  ! On output smearing is converted to a "standard" value for printout
  CHARACTER(LEN=*), INTENT(INOUT) :: smearing
  ! On output degauss is set to 0 if smearing is not used
  REAL(dp), INTENT(INOUT) :: degauss
  LOGICAL, INTENT(OUT) :: lfixed, lgauss, ltetra
  INTEGER, INTENT(OUT) :: tetra_type, ngauss
  !
  lfixed = .FALSE.
  ltetra = .FALSE.
  tetra_type = 0
  lgauss = .FALSE. 
  ngauss = 0
  !
  SELECT CASE( trim( occupations ) )
  CASE( 'fixed' )
     !
     IF ( degauss /= 0.D0 ) THEN
        CALL errore( ' set_occupations ', &
             & ' fixed occupations, gaussian broadening (degauss) ignored', -1 )
        degauss = 0.D0
     ENDIF
     smearing = 'none'
     !
  CASE( 'smearing' )
     !
     lgauss = ( degauss > 0.0_dp ) 
     IF ( .NOT. lgauss ) CALL errore( ' set_occupations ',  &
          ' smearing requires a vaklue for gaussian broadening (degauss)', 1 )
     !
     SELECT CASE ( trim( smearing ) )
     CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
        ngauss = 0
        smearing = 'Gaussian'
     CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
        ngauss = 1
        smearing = 'Methfessel-Paxton'
     CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
        ngauss = -1
        smearing = 'Marzari-Vanderbilt'
     CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
        ngauss = -99
        smearing = 'Fermi-Dirac'
     CASE DEFAULT
        CALL errore( ' set_occupations ', &
             ' smearing '//trim(smearing)//' unknown', 1 )
     END SELECT
     !
  CASE( 'tetrahedra' )
     !
     ltetra = .true.
     tetra_type = 0
     smearing = 'none'
     degauss = 0.0_dp
     !
  CASE( 'tetrahedra_lin', 'tetrahedra-lin')
     !
     ltetra = .true.
     tetra_type = 1
     smearing = 'none'
     degauss = 0.0_dp
     !
  CASE('tetrahedra_opt', 'tetrahedra-opt')
     !
     ltetra = .true.
     tetra_type = 2
     smearing = 'none'
     degauss = 0.0_dp
     !
  CASE( 'from_input' )
     !
     lfixed = .true.
     smearing = 'none'
     degauss = 0.0_dp
     !
  CASE DEFAULT
     !
     CALL errore( 'set_occupations', &
          'occupations ' // trim( occupations ) // ' not implemented', 1 )
     !
  END SELECT
  !
END SUBROUTINE set_occupations
!---------------------------------------------------------------------------
FUNCTION schema_smearing( smearing )
  !------------------------------------------------------------------------
  ! Converts smearing to the standard value needed in xml file
  !
  USE kinds, ONLY: dp
  ! 
  IMPLICIT NONE 
  !
  CHARACTER(LEN=*), INTENT(IN) :: smearing
  CHARACTER(LEN=8) :: schema_smearing
  !
  SELECT CASE ( trim( smearing ) )
  CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
     schema_smearing = 'gaussian'
  CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
     schema_smearing = 'mp'
  CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
     schema_smearing = 'mv'
  CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
     schema_smearing = 'fd'
  CASE DEFAULT
     schema_smearing = TRIM(smearing(1:8))
  END SELECT
  !
END FUNCTION schema_smearing
