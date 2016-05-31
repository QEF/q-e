!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_int_forces()
!----------------------------------------------------------------------------
!
!
USE mp,               ONLY : mp_bcast
USE mp_images,        ONLY : intra_image_comm
USE io_global,        ONLY : stdout, ionode, ionode_id
USE kinds,            ONLY : DP
USE cell_base,        ONLY : at, bg, alat, omega
USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, amass
USE fft_base,         ONLY : dfftp
USE fft_interfaces,   ONLY : fwfft
USE gvect,            ONLY : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm
USE lsda_mod,         ONLY : nspin
USE force_mod,        ONLY : force
USE scf,              ONLY : rho
USE vlocal,           ONLY : strf, vloc
USE control_flags,    ONLY : iverbosity, gamma_only
USE martyna_tuckerman, ONLY: do_comp_mt, wg_corr_force
!
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
! aux is used to store a possible additional density
! now defined in real space
!
COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
!
INTEGER  :: ipol, na
  ! counter on polarization
  ! counter on atoms
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_int_forces
