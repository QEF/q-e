MODULE kohn_sham_mod

   USE kinds, ONLY: DP
   use scf_result_mod, only: scf_result, multiple_scf_result, scf_result_set_global_variables

contains

   subroutine init_kohn_sham()
      use qpoint, only: ikqs
      use wvfct, only: nbnd, npwx
      implicit none

      if (.not. allocated(ikqs)) allocate (ikqs(1))
      ikqs(1) = 1

   end subroutine

   subroutine current_kohn_sham(J, J_a, J_b, J_el, dt, scf_all, &
                                dvpsi_save, save_dvpsi, &
                                nbnd, npw, npwx, dffts, evc, g, ngm, gstart, &
                                tpiba2, at, vkb, nkb, xk, igk_k, g2kin, et)
      use kinds, only: DP
      !use wvfct, only: nbnd, npw, npwx
      !use wavefunctions, only: psic, evc
      !use gvect, only: g, ngm, gstart
      !USE cell_base, ONLY: tpiba, omega, tpiba2, alat, at
      !use ions_base, only: tau
      !use uspp, ONLY: vkb, nkb
      !use klist, only: xk, igk_k
      !use wvfct, ONLY: g2kin, et
      use fft_base, only: fft_type_descriptor !,dffts
      use mp, only: mp_sum
      use io_global, only: ionode
      USE constants, ONLY: e2, fpi, pi
      USE fft_interfaces, ONLY: fwfft, invfft
      use becmod
      USE eqv, ONLY: dpsi, dvpsi !TODO: those variables looks like they can be local variables
      USE mp_pools, ONLY: intra_pool_comm
      USE funct, ONLY: get_igcx, get_igcc
      use compute_charge_mod, only: compute_charge
      use project_mod, only: project


      use test_h_psi, only : test

      implicit none

      type(multiple_scf_result), intent(in) :: scf_all
      logical, intent(in) :: save_dvpsi
      complex(dp), intent(inout) :: dvpsi_save(:, :, :)

      INTEGER, intent(in) :: nbnd, npwx
      INTEGER, intent(inout) :: npw, igk_k(:, :)
      TYPE(fft_type_descriptor), intent(inout) :: dffts
      COMPLEX(DP), intent(inout) ::  evc(:, :)
      INTEGER, intent(in) :: ngm, gstart, nkb
      REAL(DP), intent(inout) :: g2kin(:)
      REAL(DP), intent(out) :: J(3), J_a(3), J_b(3), J_el(3)
      REAL(DP), intent(in) ::  tpiba2, g(:, :), at(:, :), &
                              xk(:, :), et(:, :), dt
      COMPLEX(DP), intent(in) :: vkb(:, :)

      !character(LEN=20) :: dft_name
      complex(kind=DP), allocatable ::  evp(:, :), tmp(:, :)
      integer ::  ibnd
      real(kind=DP) :: amodulus
      real(DP) :: emme(nbnd, nbnd), kcurrent(3), &
                  kcurrent_a(3), kcurrent_b(3), ecurrent(3), &
                  sa(nbnd, nbnd), ssa(nbnd, nbnd), sb(nbnd, nbnd), ssb(nbnd, nbnd)

      integer  :: jbnd, ig, ipol
      allocate (tmp(npwx, nbnd))
      allocate (evp(npwx, nbnd))

      !---------------------------------KOHN------------------------------------------------
      call start_clock('kohn_current')
      call scf_result_set_global_variables(scf_all%t_zero) ! this set evc, becp, vkb, tau and vel from saved values
      allocate (dpsi(npwx, nbnd))
      allocate (dvpsi(npwx, nbnd))
! For preconditioning:
      do ig = 1, npw
         igk_k(ig, 1) = ig
         g2kin(ig) = ((xk(1, 1) + g(1, igk_k(ig, 1)))**2 + &
                      (xk(2, 1) + g(2, igk_k(ig, 1)))**2 + &
                      (xk(3, 1) + g(3, igk_k(ig, 1)))**2)*tpiba2
      end do
! init potentials needed to evaluate  H|psi>
      call init_us_1()
      call init_us_2(npw, igk_k(1, 1), xk(1, 1), vkb)
      call allocate_bec_type(nkb, nbnd, becp)
      call calbec(npw, vkb, evc, becp)
      call test

      sa = 0.d0
      sb = 0.d0

! computed sb = < evc_due, evc_uno >, sa = <evc_uno, evc_uno> remove contribution at G=0
! sb is a multiple of the identity if 2 points are used (t_minus==t_zero)

      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, scf_all%t_plus%evc, 2*npwx, scf_all%t_zero%evc, 2*npwx, 0.d0, sa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, scf_all%t_minus%evc, 2*npwx, scf_all%t_zero%evc, 2*npwx, 0.d0, sb, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               sa(ibnd, jbnd) = sa(ibnd, jbnd) - dble(conjg(scf_all%t_plus%evc(1, ibnd))*scf_all%t_zero%evc(1, jbnd))
               sb(ibnd, jbnd) = sb(ibnd, jbnd) - dble(conjg(scf_all%t_minus%evc(1, ibnd))*scf_all%t_zero%evc(1, jbnd))
            end do
         end do
      end if
      call mp_sum(sa, intra_pool_comm)
      call mp_sum(sb, intra_pool_comm)
! computes phi_v^c_punto using <s,s>

      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sa, nbnd, sa, nbnd, 0.d0, ssa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sb, nbnd, sb, nbnd, 0.d0, ssb, nbnd)

      evp = 0.d0
      do ibnd = 1, nbnd
         do jbnd = 1, nbnd
            do ig = 1, npw

               evp(ig, ibnd) = evp(ig, ibnd) + scf_all%t_plus%evc(ig, jbnd)*sa(jbnd, ibnd)
               evp(ig, ibnd) = evp(ig, ibnd) - scf_all%t_zero%evc(ig, jbnd)*ssa(ibnd, jbnd)

               evp(ig, ibnd) = evp(ig, ibnd) - scf_all%t_minus%evc(ig, jbnd)*sb(jbnd, ibnd)
               evp(ig, ibnd) = evp(ig, ibnd) + scf_all%t_zero%evc(ig, jbnd)*ssb(ibnd, jbnd)

            end do
         end do
      end do
      evp(:, :) = evp(:, :)/dt

!   end if

      tmp(:, :) = (0.d0, 0.d0)
      !does H|evp>  i.e. H|phi_v^c_punto> and saves it into tmp:
      call h_psi(npwx, npw, nbnd, evp, tmp)
      kcurrent = 0.d0
      kcurrent_a = 0.d0
      kcurrent_b = 0.d0
      ecurrent = 0.d0
      polariz: do ipol = 1, 3
         call start_clock('project')
         ! computes projection
         call project(ipol, dvpsi_save, save_dvpsi)
         call stop_clock('project')
         call print_clock('project')

         call start_clock('kohn-first')
         emme = 0.d0

         ! scalar product
         call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, dvpsi, 2*npwx, tmp, 2*npwx, 0.d0, emme, nbnd)
         if (gstart == 2) then
            do ibnd = 1, nbnd
               do jbnd = 1, nbnd
                  emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(dvpsi(1, ibnd))*tmp(1, jbnd))
               end do
            end do
         end if
         call mp_sum(emme, intra_pool_comm)

         do ibnd = 1, nbnd
            kcurrent(ipol) = kcurrent(ipol) + emme(ibnd, ibnd)
            kcurrent_a(ipol) = kcurrent_a(ipol) + emme(ibnd, ibnd)
         end do
         call stop_clock('kohn-first')
         call print_clock('kohn-first')

         call start_clock('kohn-second')
         emme = 0.d0
         call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evp, 2*npwx, dvpsi, 2*npwx, 0.d0, emme, nbnd)
         if (gstart == 2) then
            do ibnd = 1, nbnd
               do jbnd = 1, nbnd
                  emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(evp(1, ibnd))*dvpsi(1, jbnd))
               end do
            end do
         end if
         call mp_sum(emme, intra_pool_comm)

         ! et containes the eigvenvalues
         do ibnd = 1, nbnd
            kcurrent(ipol) = kcurrent(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
            kcurrent_b(ipol) = kcurrent_b(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
            ecurrent(ipol) = ecurrent(ipol) + emme(ibnd, ibnd)
         end do

         call stop_clock('kohn-second')
         call print_clock('kohn-second')

         if (ionode) write (*, "('  KOHN POLARIZATION ',I3,' COMPLETED')") ipol

      end do polariz
      J = 0.d0
      J_a = 0.d0
      J_b = 0.d0
      J_el = 0.d0
      do ipol = 1, 3
!     at(:, ipol) / amodulus is the versor along direction ipol
!      amodulus = sqrt(at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2)
         amodulus = at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2
         J(:) = J(:) + 2.d0*at(:, ipol)*real(kcurrent(ipol))/amodulus
         J_a(:) = J_a(:) + 2.d0*at(:, ipol)*real(kcurrent_a(ipol))/amodulus
         J_b(:) = J_b(:) + 2.d0*at(:, ipol)*real(kcurrent_b(ipol))/amodulus
         J_el(:) = J_el(:) + 2.d0*2.d0*at(:, ipol)*real(ecurrent(ipol))/amodulus
!    J(ipol)=J(ipol)+2.d0*real(kcurrent(ipol))
      end do
      call stop_clock('kohn_current')
      call print_clock('kohn_current')
      !npw=npwold
      if (ionode) print *, 'KOHN CURRENT CALCULATED'
!---------------------------------------------------------------------------
      call deallocate_bec_type(becp)
      deallocate (dpsi)
      deallocate (dvpsi)
      deallocate (tmp)
      deallocate (evp)

   end subroutine

END MODULE

