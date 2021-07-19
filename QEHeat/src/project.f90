!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
module project_mod

contains
   subroutine project(ipol, dvpsi_save, save_dvpsi)
      !----------------------------------------------------------------------
      !
      ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis
      !            (projected on at(*,ipol) )
      !
      ! dvpsi is READ from file if this_pcxpsi_is_on_file(ik,ipol)=.true.
      ! otherwise dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set)
      !
      USE io_files, ONLY: nwordwfc, tmp_dir
      USE kinds, ONLY: DP
      USE cell_base, ONLY: tpiba2, at
      USE io_global, ONLY: stdout, ionode
      USE klist, ONLY: xk
      USE gvect, ONLY: g, gstart
      USE wvfct, ONLY: npw, npwx, nbnd, g2kin, et
      USE klist, ONLY: igk_k
      USE wavefunctions, ONLY: evc
      USE noncollin_module, ONLY: noncolin, npol
      USE becmod, ONLY: bec_type, becp, calbec, &
                        allocate_bec_type, deallocate_bec_type
      USE uspp, ONLY: okvan, nkb, vkb
      USE uspp_param, ONLY: nh, nhm
      USE ramanm, ONLY: eth_rps
      USE eqv, ONLY: dpsi, dvpsi
      USE lrus, ONLY: becp1
      USE qpoint, ONLY: nksq
      USE units_ph, ONLY: this_pcxpsi_is_on_file, lrcom, iucom, &
                          lrebar, iuebar
      USE control_lr, ONLY: nbnd_occ, alpha_pv
      use mp, ONLY: mp_sum, mp_min, mp_max
      USE mp_global, ONLY: inter_pool_comm, intra_pool_comm
!      USE eqv, ONLY: evq
      implicit none
      integer, intent(IN) :: ipol
      logical, intent(in) :: save_dvpsi
      complex(dp), intent(inout) :: dvpsi_save(:, :, :)
      real(DP):: emin, emax
      integer :: ik
      !
      ! Local variables
      !
      integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, &
                 nrec, is, js, ijs, kbnd, ipw
      ! counters
      type(bec_type) ::becp0
      real(DP), allocatable  :: h_diag(:, :)
      ! the diagonal part of h_scf
      type(bec_type) :: becp2 ! the scalar products
      real(DP) ::   anorm, thresh
      ! preconditioning cut-off
      ! the desired convergence of linter
      logical :: conv_root
      ! true if convergence has been achieved
      real(DP) ::emme(nbnd, nbnd)
!     logical ::l_test, exst

      external ch_psi_all, cg_psi
      !debug
!      allocate (evq(npwx, nbnd))
!      evq = evc
      ik = 1
      dpsi = (0.d0, 0.d0)
      dvpsi = (0.d0, 0.d0)
      call allocate_bec_type(nkb, nbnd, becp2)
      ! calculate the commutator [H,x_ipol]  psi > and store it in dpsi
      ! dvpsi used as workspace
      call allocate_bec_type(nkb, nbnd, becp0)
      CALL calbec(npw, vkb, evc, becp0)
      call commutator_Hx_psi(ik, nbnd, at(:, ipol), becp0, becp2, dpsi)
      call deallocate_bec_type(becp0)
      !    orthogonalize dpsi to the valence subspace: ps = <evc|dpsi>
      !    Apply -P^+_c
      !    NB it uses dvpsi as workspace
      ! manual orthogonalization
      emme = 0.d0
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evc, 2*npwx, dpsi, 2*npwx, 0.d0, emme, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(evc(1, ibnd))*dpsi(1, jbnd))
            end do
         end do
      end if
      call mp_sum(emme, intra_pool_comm)
      call dgemm('N', 'N', 2*npw, nbnd, nbnd, -1.d0, evc, 2*npwx, emme, nbnd, 1.d0, dpsi, 2*npwx)
      !dpsi=-dpsi
      !
      !   dpsi contains P^+_c [H-eS,x] psi_v for the three crystal polarizations
      !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
      !
      do ig = 1, npw
         g2kin(ig) = SUM((g(1:3, igk_k(ig, 1)))**2)*tpiba2
      enddo
      allocate (h_diag(npwx*npol, nbnd))
      h_diag = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! preconditioning
      do ibnd = 1, nbnd
         do ig = 1, npw
            h_diag(ig, ibnd) = 1.d0/max(1.0d0, g2kin(ig))
         enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! computation of alpha_pv
      emin = et(1, 1)
      do ibnd = 1, nbnd
         emin = min(emin, et(ibnd, 1))
      enddo
#ifdef __MPI
      ! find the minimum across pools
      call mp_min(emin, inter_pool_comm)
#endif
      emax = et(1, 1)
      do ibnd = 1, nbnd
         emax = max(emax, et(ibnd, 1))
      enddo
#ifdef __MPI
      ! find the maximum across pools
      call mp_max(emax, inter_pool_comm)
#endif
      alpha_pv = 2.d0*(emax - emin)
      ! avoid zero value for alpha_pv
      alpha_pv = max(alpha_pv, 1.0d-2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! end computation of alpha_pv
      if (save_dvpsi) then
         dvpsi = dvpsi_save(:, :, ipol)
      else
         dvpsi(:, :) = (0.d0, 0.d0)
      end if
      eth_rps = 1.D-10
      thresh = eth_rps
      !dvpsi is the initial estimate of the solution
      call cgsolve_all(ch_psi_all, cg_psi, et(1, 1), dpsi, dvpsi, &
                       h_diag, npwx, npw, thresh, 1, lter, conv_root, anorm, &
                       nbnd, npol)
      if (.not. conv_root) WRITE (stdout, '(5x,"ik",i4," ibnd",i4, &
           & " linter: root not converged ",e10.3)') &
           ik, ibnd, anorm
      if (save_dvpsi) then
         dvpsi_save(:, :, ipol) = dvpsi
      end if
      deallocate (h_diag)
      !
      ! we have now obtained P_c x |psi>.
      !
      IF (nkb > 0) call deallocate_bec_type(becp2)
!      deallocate (evq)
      return
   end subroutine project

end module
