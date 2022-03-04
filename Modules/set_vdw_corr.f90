
! Copyright (C) 2019 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw, mbd_vdw, lxdm )
  USE io_global, ONLY: stdout
  !
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: vdw_corr
  LOGICAL, INTENT(out) :: llondon, ldftd3, ts_vdw, mbd_vdw, lxdm
  !

  llondon= .FALSE.
  ldftd3 = .FALSE.
  ts_vdw = .FALSE.
  mbd_vdw= .FALSE.
  lxdm   = .FALSE.

  SELECT CASE( TRIM( vdw_corr ) )
  CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
     llondon= .TRUE.

  CASE( 'grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3' )
     ldftd3 = .TRUE.

  CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
     ts_vdw = .TRUE.

  CASE( 'MBD', 'mbd', 'many-body-dispersion', 'mbd_vdw' )
     ts_vdw = .TRUE.
     mbd_vdw = .TRUE.

  CASE( 'XDM', 'xdm' )
     lxdm   = .TRUE.

  CASE('none','')

  CASE DEFAULT
     WRITE (stdout,*)
     CALL infomsg('set_vdw_corr','WARNING: unknown vdw correction (vdw_corr): '//TRIM(vdw_corr)//'. No vdw correction used.')
     WRITE (stdout,*)

  END SELECT
    

END SUBROUTINE set_vdw_corr
