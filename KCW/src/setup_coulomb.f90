!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
SUBROUTINE setup_coulomb ( )
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE cell_base,            ONLY : at, alat
  USE io_global,            ONLY : stdout, ionode
  USE exx_base,             ONLY : use_coulomb_vcut_ws
  !
  IMPLICIT NONE
  !
  REAL(DP)     :: atws(3,3)
  REAL(DP)     :: ecutvcut
  TYPE(vcut_type)   :: vcut
  !
  ! build the superperiodicity direct lattice
        !
  use_coulomb_vcut_ws = .true.
  atws = alat * at
  !WRITE(*,*) "NICOLA", alat, at
  !
  atws(:,1) = atws(:,1) * 2
  atws(:,2) = atws(:,2) * 2
  atws(:,3) = atws(:,3) * 2
  !
  ecutvcut = 0.D0 ! For now FIXME
  !
  WRITE(*,*) "NICOLA", ecutvcut
  !
  CALL vcut_init( vcut, atws, ecutvcut )
  !
  IF ( ionode ) CALL vcut_info( stdout, vcut )
  !
END SUBROUTINE
