!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE test_f_sum_rule
  !-----------------------------------------------------------------------
  !
  ! ... Test the f-sum rule
  ! ... 2/N_k sum_k sum_n^{occ} 
  ! ...   <u_{nk}| p_{k,\alpha} G_k v_{k,\beta} | u_{nk}> = 
  ! ...   = -N_{el} \delta_{\alpha,\beta}
  ! ...
  ! ... where: p_k = -i\nabla + k and v_k = -i [r, H_k]  
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                          current_k
  USE lsda_mod,                    ONLY : current_spin, lsda, isk
  USE becmod,                      ONLY : becp
  USE pwcom
  USE nmr_module

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
  complex(dp), allocatable, dimension(:,:,:) :: p_evc, vel_evc, g_vel_evc
  real(dp) :: f_sum(3,3), f_sum_k(3,3), q(3)
  integer :: ik, ipol, jpol, ibnd, ig
  complex(dp), external :: ZDOTC

  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3),   &
             vel_evc(npwx,nbnd,3), &
             g_vel_evc(npwx,nbnd,3) )

  ! zero the f-sum
  f_sum(:,:) = 0.d0
  q(:) = 0.d0

  write(stdout, '(5X,''Computing the f-sum rule'')')

  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)

    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw,igk,xk(1,ik),vkb)

    ! read wfcs from file and compute becp
    CALL davcio (evc, nwordwfc, iunwfc, ik, -1)
    call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)

    !q = 0.d0; q(1) = 0.0d0
    !write(stdout, '(5X,''computing wfcs at k + q'')')
    !call compute_u_kq(ik, q)
    !evc = evq
 
    ! compute p_k|evc>, v_k|evc> and G_k v_k|evc>
    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      call greenfunction(ik, vel_evc(1,1,ipol), g_vel_evc(1,1,ipol), q)
    enddo

    ! k-point contribution to the f-sum rule
    f_sum_k = 0.0

    ! loop over cartesian directions
    do jpol = 1, 3
      do ipol = 1, 3
        do ibnd = 1, nbnd_occ (ik)
          f_sum_k(ipol,jpol) = f_sum_k(ipol,jpol) + wg(ibnd,ik) * &
            2.d0 * real(ZDOTC(npw, p_evc(1,ibnd,ipol), 1, &
                                   g_vel_evc(1,ibnd,jpol), 1))
        enddo
      enddo   ! ipol
    enddo   ! jpol

    write(stdout, '(5X,''f-sum rule (ik='',I5,''):'')') ik
    write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum_k

    f_sum(:,:) = f_sum(:,:) + f_sum_k(:,:)
  enddo   ! ik
#ifdef __PARA
  call reduce(9, f_sum)
  call poolreduce(9, f_sum)
#endif
  write(stdout, '(5X,''f-sum rule:'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum(:,:)

  call sym_cart_tensor(f_sum)
  write(stdout, '(5X,''f-sum rule (symmetrized):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum

  deallocate(p_evc, vel_evc, g_vel_evc)

END SUBROUTINE test_f_sum_rule


