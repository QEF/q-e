!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_corr_gpu( forcescc )
  !-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and read from pseudopotential files
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : msh, rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, g_d, ngl, gl_d, igtongl_d
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE device_fbuff_m,             ONLY : dev_buf
  !
  USE simpsn_gpum,          ONLY : simpsn_gpu_dev
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  !
  implicit none
  !
  real(DP) :: forcescc(3,nat)
  !
#if defined(__CUDA)
  real(DP), pointer, contiguous :: rhocgnt_d(:),  aux_d(:,:), r_d(:), rab_d(:), rhoat_d(:), tau_d(:,:)
  real(DP), allocatable :: vauxr(:)
  complex(DP), allocatable :: vauxg(:,:)
  ! work space
  real(DP) ::  gx, arg, fact, forcesccx, forcesccy, forcesccz
  ! temp factors
  integer :: ir, isup, isdw, ig, nt, na, ndm, ierr, msh_nt, igg, glblock,igl, ierrs(6)
  ATTRIBUTES(DEVICE)  :: rhocgnt_d, aux_d, r_d, rab_d, rhoat_d, tau_d
  !
  !
  ndm = MAXVAL( msh(1:ntyp) )
  CALL dev_buf%lock_buffer( rhocgnt_d, ngl, ierrs(1) )
  CALL dev_buf%lock_buffer( r_d, ndm, ierrs(2) )
  CALL dev_buf%lock_buffer( rab_d, ndm, ierrs(3) )
  CALL dev_buf%lock_buffer( rhoat_d, ndm, ierrs(4) )
  CALL dev_buf%lock_buffer( aux_d, [ndm, ngl], ierrs(5) )
  CALL dev_buf%lock_buffer( tau_d, [3,nat], ierrs(6) )
  IF (ANY(ierrs /= 0)) CALL errore('force_corr_gpu', 'cannot allocate buffers', &
                                                             ABS(MAXVAL(ierrs)) )
  !
  tau_d(1:3,1:nat) = tau(1:3,1:nat)
  !
  allocate( vauxr(dfftp%nnr), vauxg(dfftp%nnr,1) )
  !
  if (nspin == 1 .or. nspin == 4) then
     vauxr(:) = vnew%of_r(:,1)
  else
     isup = 1
     isdw = 2
     vauxr(:) = (vnew%of_r(:,isup) + vnew%of_r(:,isdw)) * 0.5d0
  endif
  !
  !$acc data copyin(vauxr) create(vauxg)
  !
  call rho_r2g( dfftp, vauxr, vauxg )
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  endif
  !
  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
     msh_nt = msh(nt)
     r_d(1:msh_nt) = rgrid(nt)%r(1:msh_nt)
     rab_d(1:msh_nt) = rgrid(nt)%rab(1:msh_nt)
     rhoat_d(1:msh_nt) = upf(nt)%rho_at(1:msh_nt)
     !
     !$cuf kernel do(1)
     do ig = 1, ngl
        gx = sqrt(gl_d(ig)) * tpiba
        do ir = 1, msh_nt
           if (r_d(ir) < 1.0d-8) then
              aux_d(ir,ig) = rhoat_d(ir)
           else
              aux_d(ir,ig) = rhoat_d(ir) * sin(gx*r_d(ir)) / (r_d(ir)*gx)
           endif
        enddo
        !
        call simpsn_gpu_dev( msh_nt, aux_d(:,ig), rab_d, rhocgnt_d(ig) )
        !
     enddo
     !
     do na = 1, nat
        if ( nt == ityp(na) ) then
           forcescc(1:3,na) = 0.0_DP
           forcesccx = 0.0d0
           forcesccy = 0.0d0
           forcesccz = 0.0d0
           !$acc parallel loop
           do ig = gstart, ngm
              arg = (g_d(1,ig) * tau_d(1,na) + g_d(2,ig) * tau_d(2,na) &
                   + g_d(3,ig) * tau_d(3,na) ) * tpi
              forcesccx = forcesccx + REAL(fact * &
                      rhocgnt_d(igtongl_d(ig)) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(1,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccy = forcesccy + REAL(fact * &
                      rhocgnt_d(igtongl_d(ig)) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(2,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccz = forcesccz + REAL(fact * &
                      rhocgnt_d(igtongl_d(ig)) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(3,ig) * tpiba * CONJG(vauxg(ig,1)))
           enddo
           forcescc(1, na) = forcesccx
           forcescc(2, na) = forcesccy
           forcescc(3, na) = forcesccz
        endif
     enddo
  enddo
  !
  !$acc end data
  deallocate( vauxr, vauxg )
  !
  call dev_buf%release_buffer( r_d, ierr )
  call dev_buf%release_buffer( rab_d, ierr )
  call dev_buf%release_buffer( rhoat_d, ierr )
  call dev_buf%release_buffer( aux_d, ierr )
  !
  call dev_buf%release_buffer( tau_d, ierr )
  call dev_buf%release_buffer( rhocgnt_d, ierr )
  !
  call mp_sum( forcescc, intra_bgrp_comm )
#endif
  !
  return
end subroutine force_corr_gpu

