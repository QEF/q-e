
! Copyright (C) 2019 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw, lxdm )
  !
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: vdw_corr
  LOGICAL, INTENT(out) :: llondon, ldftd3, ts_vdw, lxdm
  !
  SELECT CASE( TRIM( vdw_corr ) )
     !
  CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
     !
     llondon= .TRUE.
     ldftd3 = .FALSE.
     ts_vdw = .FALSE.
     lxdm   = .FALSE.
     !
  CASE( 'grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3' )
     !
     ldftd3 = .TRUE.
     llondon= .FALSE.
     ts_vdw = .FALSE.
     lxdm   = .FALSE.
     !
  CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
     !
     llondon= .FALSE.
     ldftd3 = .FALSE.
     ts_vdw = .TRUE.
     lxdm   = .FALSE.
     !
  CASE( 'XDM', 'xdm' )
     !
     llondon= .FALSE.
     ldftd3 = .FALSE.
     ts_vdw = .FALSE.
     lxdm   = .TRUE.
     !
  CASE DEFAULT
     !
     llondon= .FALSE.
     ldftd3 = .FALSE.
     ts_vdw = .FALSE.
     lxdm   = .FALSE.
     !
  END SELECT
    
END SUBROUTINE set_vdw_corr
