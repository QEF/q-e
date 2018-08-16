!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine remove_atomic_rho
  !-----------------------------------------------------------------------
  USE io_global, ONLY: stdout
  USE io_files, ONLY: output_drho, tmp_dir, prefix, postfix
  USE kinds, ONLY: DP
  USE control_flags,    ONLY : gamma_only
  USE gvect,            ONLY : ig_l2g, mill
  USE fft_base, ONLY: dfftp
  USE fft_rho, ONLY: rho_r2g
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho
  USE io_base, ONLY : write_rhog
  USE mp_pools,         ONLY : my_pool_id
  USE mp_bands,         ONLY : my_bgrp_id, root_bgrp_id, &
       root_bgrp, intra_bgrp_comm
  USE cell_base,        ONLY : bg, tpiba
  !
  implicit none
  CHARACTER(LEN=256) :: filename
  real(DP), allocatable :: work (:,:)
  ! workspace, is the difference between the charge density
  ! and the superposition of atomic charges
  complex(DP), allocatable :: workc(:,:)
  !
  IF ( nspin > 1 ) CALL errore &
       ( 'remove_atomic_rho', 'spin polarization not allowed in drho', 1 )

  WRITE( stdout, '(/5x,"remove atomic charge density from scf rho")')
  !
  !     subtract the old atomic charge density (FIXME: in real space)
  !
  allocate ( work( dfftp%nnr, 1 ) )
  work = 0.d0
  call atomic_rho (work, nspin)
  work = rho%of_r - work
  !
  !      FIXME: move to G-space
  !
  allocate ( workc( dfftp%ngm, 1) )
  CALL rho_r2g ( dfftp, work, workc )
  deallocate(work)
  !
  filename = TRIM(tmp_dir) // TRIM(prefix) // postfix // output_drho
  IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
       CALL write_rhog( filename, root_bgrp, intra_bgrp_comm, &
       bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
       gamma_only, mill, ig_l2g, workc )
  !
  deallocate(work)
  return

end subroutine remove_atomic_rho

