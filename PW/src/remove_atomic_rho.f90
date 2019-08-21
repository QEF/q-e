!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE remove_atomic_rho
  !-----------------------------------------------------------------------
  USE kinds,        ONLY: DP
  USE io_global,    ONLY: stdout
  USE io_files,     ONLY: output_drho, restart_dir
  USE control_flags,ONLY: gamma_only
  USE gvect,        ONLY: ngm, ig_l2g, mill
  USE lsda_mod,     ONLY: nspin
  USE scf,          ONLY: rho
  USE io_base,      ONLY: write_rhog
  USE mp_pools,     ONLY: my_pool_id
  USE mp_bands,     ONLY: my_bgrp_id, root_bgrp_id, &
       root_bgrp, intra_bgrp_comm
  USE cell_base,    ONLY: bg, tpiba
  !
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: drhog(:,:)
  !
  IF ( nspin > 1 ) CALL errore &
       ( 'remove_atomic_rho', 'spin polarization not allowed in drho', 1 )

  WRITE( stdout, '(/5x,"remove atomic charge density from scf rho")')
  !
  !     subtract the old atomic charge density
  !
  ALLOCATE ( drhog( ngm, nspin) )
  CALL atomic_rho_g ( drhog, nspin)
  drhog = rho%of_g - drhog
  !
  IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
       CALL write_rhog( TRIM( restart_dir( ) ) // output_drho, &
       root_bgrp, intra_bgrp_comm, &
       bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
       gamma_only, mill, ig_l2g, drhog )
  !
  DEALLOCATE(drhog)
  !
END SUBROUTINE remove_atomic_rho

