!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_corr_gpu (forcescc)
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
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl_d, igtongl_d
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE device_fbuff_m,             ONLY : dev_buf
  USE gvect_gpum,           ONLY : g_d
  !
  USE simpsn_gpum,          ONLY : simpsn_gpu_dev
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  !
  implicit none
  !
  real(DP) :: forcescc (3, nat)
  !
  real(DP), pointer, contiguous ::   rhocgnt_d(:),  aux_d (:,:), r_d(:), rab_d(:), rhoat_d(:), tau_d(:,:)
  complex(DP),pointer      ::   psic_d(:)
  integer, pointer      :: nl_d (:)
  real(DP), allocatable :: tmp(:)

  ! work space
  real(DP) ::  gx, arg, fact, forcesccx, forcesccy, forcesccz
  ! temp factors
  integer :: ir, isup, isdw, ig, nt, na, ndm, ierr, msh_nt, igg, glblock,igl, ierrs(7)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE)  :: rhocgnt_d, aux_d, psic_d,r_d, rab_d, rhoat_d, nl_d, tau_d
#endif

  ! counters
  !
  ! vnew is V_out - V_in, psic is the temp space
  !
  nl_d => dfftp%nl_d
  CALL dev_buf%lock_buffer(psic_d, size(vnew%of_r,1), ierrs(1) )
  if (nspin == 1 .or. nspin == 4) then
     psic_d(:) = vnew%of_r (:, 1)
  else
     isup = 1
     isdw = 2
     psic_d(:) = (vnew%of_r (:, isup) + vnew%of_r (:, isdw)) * 0.5d0
  end if
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  CALL dev_buf%lock_buffer ( rhocgnt_d, ngl, ierrs(2) )
  CALL dev_buf%lock_buffer ( r_d, ndm, ierrs(3))
  CALL dev_buf%lock_buffer ( rab_d, ndm, ierrs(4) )
  CALL dev_buf%lock_buffer ( rhoat_d, ndm, ierrs(5))
  CALL dev_buf%lock_buffer ( aux_d, [ndm, ngl], ierrs(6) )
  !
  CALL dev_buf%lock_buffer ( tau_d, [3,nat], ierrs(7))
  !
  IF (ANY(ierrs /= 0)) CALL errore('force_corr_gpu', 'cannot allocate buffers', -1)
  !
  tau_d(1:3,1:nat)=tau(1:3,1:nat)
  !
  CALL fwfft ('Rho', psic_d, dfftp)
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
     msh_nt = msh(nt)
     r_d (1:msh_nt) =  rgrid(nt)%r(1:msh_nt)
     rab_d(1:msh_nt) = rgrid(nt)%rab(1:msh_nt)
     rhoat_d (1:msh_nt) =  upf(nt)%rho_at (1:msh_nt)
     !
     !$cuf kernel do(1)
     do ig = 1, ngl
        gx = sqrt (gl_d (ig) ) * tpiba
        do ir = 1, msh_nt
           if (r_d(ir) .lt.1.0d-8) then
              aux_d (ir,ig) = rhoat_d (ir)
           else
              aux_d (ir,ig) = rhoat_d (ir) * sin(gx*r_d(ir)) / (r_d(ir)*gx)
           endif
        enddo

        call simpsn_gpu_dev(msh_nt, aux_d(:,ig), rab_d, rhocgnt_d(ig) )

     enddo
     !
     do na = 1, nat
        if (nt.eq.ityp (na) ) then
           forcescc (1:3, na) = 0.0_DP
           forcesccx = 0.0d0
           forcesccy = 0.0d0
           forcesccz = 0.0d0
           !$cuf kernel do(1)
           do ig = gstart, ngm
              arg = (g_d(1, ig) * tau_d (1, na) + g_d (2, ig) * tau_d (2, na) &
                   + g_d (3, ig) * tau_d (3, na) ) * tpi
              forcesccx = forcesccx+ REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(1,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
              forcesccy = forcesccy + REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(2,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
              forcesccz = forcesccz + REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(3,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
           enddo
           forcescc(1, na) = forcesccx
           forcescc(2, na) = forcesccy
           forcescc(3, na) = forcesccz
        endif
     enddo
  enddo
  !
  call dev_buf%release_buffer ( r_d, ierr )
  call dev_buf%release_buffer ( rab_d, ierr )
  call dev_buf%release_buffer ( rhoat_d, ierr )
  call dev_buf%release_buffer ( aux_d, ierr )
  !
  call dev_buf%release_buffer ( tau_d, ierr )
  call dev_buf%release_buffer ( psic_d, ierr )
  call dev_buf%release_buffer ( rhocgnt_d, ierr )
  !
  call mp_sum(  forcescc, intra_bgrp_comm )
  !
  return
end subroutine force_corr_gpu

