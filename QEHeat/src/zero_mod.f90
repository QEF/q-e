!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE zero_mod

   USE kinds, ONLY: DP
   !use splines
   use fft_base
   use becmod

contains

   subroutine init_zero(tabr, H_g, &
                        nsp, zv, tpiba2, tpiba, omega, at, alat, &
                        ngm, gg, gstart, g, igtongl, gl, ngl, dq, &
                        upf, rgrid, nqxq, intra_bgrp_comm, &
                        nat, ityp)
!called once to init stuff that does not depend on the atomic positions
!tabr and h_g, on output, will be used by all the current_zero calls
      use kinds, only: DP
      USE radial_grids, ONLY: radial_grid_type
      USE uspp_param, ONLY: pseudo_upf
      !use gvect, only: ngm, gg, g, gstart
      !use hartree_mod
      !use ions_base, only: nsp
      implicit none

      real(dp), intent(inout) :: tabr(:, :, :, :), H_g(:, :, :, :)
      integer, intent(in) :: nsp, igtongl(:), ngm, gstart, ngl, nqxq, &
                             intra_bgrp_comm, nat, ityp(:)
      real(dp), intent(in) :: zv(:), tpiba2, tpiba, omega, at(3, 3), alat, &
                              gg(:), g(:, :), gl(:), dq
      type(radial_grid_type), intent(in) :: rgrid(:)
      type(pseudo_upf), intent(in) :: upf(:)

      real(DP), allocatable :: tablocal_hg(:, :, :)

      integer :: isp, iun, a, b
      logical :: exst
      call start_clock('init_zero')
      allocate (tablocal_hg(nqxq, nsp, 2))
      call init_us_1a(tabr, tablocal_hg, rgrid, nsp, zv, omega, nqxq, dq, upf)
      call init_reciprocal_parts_tab(H_g, tablocal_hg, nsp, zv, tpiba2, tpiba, omega, at, alat, &
                                     ngm, gg, gstart, g, igtongl, gl, ngl, dq)
      do a = 1, 3
         do b = 1, 3
            if (a > b) then
               do isp = 1, nsp
                  H_g(:, b, a, isp) = H_g(:, a, b, isp)
               end do
            end if
         end do
      end do
      deallocate (tablocal_hg)
      call stop_clock('init_zero')

   end subroutine

   subroutine current_zero(j_zero, tabr, H_g, &
                           nbnd, npwx, npw, dffts, nsp, zv, nat, ityp, amass, tau, &
                           vel, tpiba, tpiba2, at, alat, omega, psic, evc, ngm, gg, g, gstart, &
                           nkb, vkb, deeq, upf, nh, xk, igk_k, bg, ec_test_)
      !needs tabr and h_g calculated by init_zero
      !output in j_zero
      use kinds, only: DP
      use fft_base, only: fft_type_descriptor
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      use io_global, only: stdout, ionode
      use constants, only: e2, AMU_RY
      !use splines
      !use hartree_mod
      use compute_charge_mod, only: compute_charge
      USE uspp_param, ONLY: pseudo_upf

      implicit none

      real(dp), intent(out) :: j_zero(3)
      real(dp), intent(in) :: tabr(:, :, :, :), H_g(:, :, :, :)
      INTEGER, intent(in) :: nbnd, npwx, nsp, nat, ityp(:), &
                             ngm, gstart, nkb, igk_k(:, :), nh(:)
      INTEGER, intent(inout) :: npw
      TYPE(fft_type_descriptor), intent(inout) :: dffts
      REAL(DP), intent(in) :: zv(:), amass(:), tau(:, :), vel(:, :), &
                              tpiba, omega, tpiba2, alat, at(:, :), &
                              g(:, :), gg(:), bg(:, :)
      COMPLEX(DP), intent(inout) :: psic(:), evc(:, :)
      COMPLEX(DP), intent(inout):: vkb(:, :)
      real(dp), intent(in) :: deeq(:, :, :, :), xk(:, :)
      type(pseudo_upf), intent(in) :: upf(:)
      logical, intent(in) :: ec_test_

!three dimensional auxiliary real vectors
      real(dp), parameter  :: amconv = AMU_RY
      real(DP) ::n(3), u(3), u_pbc(3), u_x(3), u_y(3), u_z(3), value, x(1:3), ics(3)
      !real(DP)  :: dtauij(3)
      real(DP) ::modul, erf_value
      real(DP) ::fac_uno, fac_due, fac_tre, fac_quattro
      real(DP), external :: qe_erf

!counters and indexes
      integer      ::a, b !cartesian components
      integer      ::igm, isp, jgm !reciprocal lattice,types, ions
      integer      ::nr1a, nr2a, nr3a
      integer      ::im, iqq !pseudo mesh, real_mesh_one_coordinate
      integer      ::ii, iv !processors,bands
      !integer      ::ix, iy, iz !real_mesh_three_coordinates
      integer      ::iun !unit
      integer      ::iatom, jatom !atoms
      integer      ::n_x, n_y, n_z, n_maxl
      integer      ::istep

!variables for calling functions
      logical :: exst, l_test, l_scambio, l_scambio_alt

!auxiliary variables
      integer ::err, ir, ieta
      real(DP) ::R
      real(DP), allocatable   :: charge(:)
      complex(DP), allocatable   :: u_g(:, :) !ngm,a
      complex(DP), allocatable ::charge_g(:)
      logical  ::l_non_loc

      if (ionode) write (stdout, *) 'ROUTINE_ZERO BEGINNING'
      call start_clock('routine_zero')
      allocate (charge(dffts%nnr))
      call start_clock('zero_current')
      if (nkb > 0) then
         l_non_loc = .true.
      else
         l_non_loc = .false.
      end if
      npw = npwx
      allocate (charge_g(ngm))
      call compute_charge(psic, evc, npw, nbnd, ngm, dffts, charge, charge_g)
!
!initialization of  u_g, Eq. 27
      allocate (u_g(ngm, 3))
      u_g = 0.d0
      do a = 1, 3
         do b = 1, 3
            do igm = 1, ngm
               do iatom = 1, nat
                  u_g(igm, a) = u_g(igm, a) - vel(b, iatom)*H_g(igm, a, b, ityp(iatom))*&
   &exp(-(0.d0, 1.d0)*DOT_PRODUCT(tpiba*g(1:3, igm), alat*tau(1:3, iatom)))
               end do
            end do
         end do
      end do
!
!computation of the current, , Eq. 26
      j_zero = 0.d0
      do a = 1, 3
         do igm = gstart, ngm
            j_zero(a) = j_zero(a) + 2.d0*dble(charge_g(igm)*conjg(u_g(igm, a)))
         end do
         if (gstart == 2) then
            j_zero(a) = j_zero(a) + dble(charge_g(1)*conjg(u_g(1, a)))
         end if
      end do
      deallocate (charge_g)
      deallocate (u_g)

      call mp_sum(j_zero, intra_pool_comm)
      if (l_non_loc) then
         ! Add non local part, Eq. 34
         call add_nc_curr(j_zero, tabr, nkb, vkb, deeq, upf, nh, vel, nbnd, npw, npwx, evc, &
                          g, tpiba, nat, ityp, nsp, xk, igk_k, ec_test_)
      end if
      j_zero = j_zero*alat
      call stop_clock('zero_current')
      call print_clock('zero_current')
      deallocate (charge)
      if (ionode) print *, 'ZERO CURRENT CALCULATED'
      call stop_clock('routine_zero')
      call print_clock('routine_zero')

   end subroutine current_zero

   subroutine add_nc_curr(current, tabr, nkb, vkb, deeq, upf, nh, vel, nbnd, npw, npwx, evc, &
                          g, tpiba, nat, ityp, nsp, xk, igk_k, ec_test_)
      use kinds, only: DP
      use becmod
      USE uspp_param, ONLY: pseudo_upf
      use init_us_3_mod, only: init_us_3
      USE uspp_init,            ONLY : init_us_2

      implicit none

      real(DP), intent(inout) :: current(3)
      integer, intent(in) :: nkb, nh(:), nbnd, npw, npwx, nat, ityp(:), &
                             nsp, igk_k(:, :)
      complex(DP), intent(in) :: evc(:, :)
      complex(DP), intent(inout) :: vkb(:, :)
      real(dp), intent(in) :: deeq(:, :, :, :), vel(:, :), g(:, :), tpiba, xk(:, :), tabr(:, :, :, :)
      type(pseudo_upf), intent(in) :: upf(:)
      logical, intent(in) :: ec_test_

      !integer ::iun, a, b
      !logical ::l_test
      integer :: ipol, jpol, ibnd
      integer ::ijkb, ikb, ih, na, nt, ipw
      real(DP) ::J_nl(3), J_1(3), J_2(3)
      complex(DP), allocatable ::vkb1(:, :)
      TYPE(bec_type) :: becpr(3)
      TYPE(bec_type) :: becpd(3)
      TYPE(bec_type) :: becprd(3, 3)
      complex(DP), allocatable :: xvkb(:, :, :)
      complex(DP), allocatable :: dvkb(:, :, :)
      complex(DP), allocatable :: xdvkb(:, :, :, :)

      allocate (vkb1(npwx, nkb))
!do nothing if there is no non local potential
      if (nkb <= 0) return
      do nt = 1, nsp
         if ((upf(nt)%typ .eq. "US") .or. (upf(nt)%typ .eq. "PAW")) then
            CALL errore('add_nc_curr', 'US and PAW not implemented', 1)
         end if
      end do
      CALL allocate_bec_type(nkb, nbnd, becp)
      do ipol = 1, 3
         CALL allocate_bec_type(nkb, nbnd, becpr(ipol))
         CALL allocate_bec_type(nkb, nbnd, becpd(ipol))
         do jpol = 1, 3
            CALL allocate_bec_type(nkb, nbnd, becprd(ipol, jpol))
         end do
      end do
      allocate (xvkb(npwx, nkb, 3))
      allocate (dvkb(npwx, nkb, 3))
      allocate (xdvkb(npwx, nkb, 3, 3))
!
!  initialization of vkb (just to be sure) and xvkb
!  - vkb is the Fourier transform of the beta functions,
!  - xvkb is the Fourier tranform of xbeta.
!  Both are needed in Eq. 36
      CALL init_us_2(npw, igk_k(:, 1), xk(:, 1), vkb)
      call init_us_3(npw, xvkb, tabr, ec_test_)
!
! initialization of dvkb and xdvkb
! - dvkb is ( - partial_a beta ) (note the minus sign in Eq. 36, 37)
! - xdvkb is ( - partial_a x_b beta ) (note the minus sign in Eq. 36, 37)
! , via identities in Eq. 37, 38
      dvkb = 0.d0
      do ipol = 1, 3
         do ikb = 1, nkb
            do ipw = 1, npw
               dvkb(ipw, ikb, ipol) = (0.d0, -1.d0)*tpiba*g(ipol, ipw)*vkb(ipw, ikb)
            end do
         end do
      end do
      do jpol = 1, 3
         do ipol = 1, 3
            do ikb = 1, nkb
               do ipw = 1, npw

                  xdvkb(ipw, ikb, ipol, jpol) = (0.d0, -1.d0)*tpiba*g(jpol, ipw)*xvkb(ipw, ikb, ipol)
                  if (ipol == jpol) then
                     xdvkb(ipw, ikb, ipol, jpol) = xdvkb(ipw, ikb, ipol, jpol) + vkb(ipw, ikb)

                  end if
               end do
            end do
         end do
      end do
! Scalar products. Note that here we have MPI communication between the nodes.
      call calbec(npw, vkb, evc, becp)
      do ipol = 1, 3
         CALL calbec(npw, xvkb(1:npwx, 1:nkb, ipol), evc, becpr(ipol))
         do ipw = 1, npw
            do ikb = 1, nkb
               vkb1(ipw, ikb) = dvkb(ipw, ikb, ipol)
            end do
         end do
         CALL calbec(npw, vkb1, evc, becpd(ipol))
         do jpol = 1, 3
            call calbec(npw, xdvkb(1:npwx, 1:nkb, ipol, jpol), evc, becprd(ipol, jpol))
         end do
      end do
! Here we loaded :
! becp -> <phi|beta>
! becpr -> <phi|xbeta>
! becpd -> <phi|-partial_a beta>
! becprd -> <phi|- partial_a x_b beta>
! Now we use these scalar products to evaluate the current, Eq. 35,36
      J_nl = 0.d0
      ! J_1 & J_2 for debugging purposes
      J_1 = 0.d0
      J_2 = 0.d0
      ijkb = 0
      do nt = 1, nsp
         do na = 1, nat
            if (ityp(na) .eq. nt) then
               do ih = 1, nh(nt)
                  ikb = ijkb + ih
                  do ipol = 1, 3
                     do jpol = 1, 3
                        do ibnd = 1, nbnd
                           J_nl(ipol) = J_nl(ipol) + vel(jpol, na)*&
   &becprd(ipol, jpol)%r(ikb, ibnd)*becp%r(ikb, ibnd)*deeq(ih, ih, na, 1)

                           J_1(ipol) = J_1(ipol) + vel(jpol, na)*&
   &becprd(ipol, jpol)%r(ikb, ibnd)*becp%r(ikb, ibnd)*deeq(ih, ih, na, 1)

!                         print*,'corrente non locale: ', J_nl(:)
!                          if (ionode) then
!                              print*,'becpr-ikb-ipol-ibnd-ityp',becpr(ipol)%r(ikb,ibnd),ikb,ipol,ibnd,ityp(na)
!                          end if
!                         print*,'VEL',vel(:,na)
                           J_nl(ipol) = J_nl(ipol) + vel(jpol, na)*&
   &becpr(ipol)%r(ikb, ibnd)*becpd(jpol)%r(ikb, ibnd)*deeq(ih, ih, na, 1)

                           J_2(ipol) = J_2(ipol) + vel(jpol, na)*&
   &becpr(ipol)%r(ikb, ibnd)*becpd(jpol)%r(ikb, ibnd)*deeq(ih, ih, na, 1)
!                       print*,'corrente',J_nl(ipol)
                        end do
                     end do
                  end do
               end do
               ijkb = ijkb + nh(nt)
            end if

         end do
      end do
!
! the factor 2 is there for spin degeneracy
      current(:) = current(:) + 2.d0*J_nl(:)
!
!free memory
      CALL deallocate_bec_type(becp)
      do ipol = 1, 3
         CALL deallocate_bec_type(becpr(ipol))
         CALL deallocate_bec_type(becpd(ipol))
         do jpol = 1, 3
            CALL deallocate_bec_type(becprd(ipol, jpol))
         end do
      end do
      deallocate (xvkb)
      deallocate (dvkb)
      deallocate (xdvkb)
!
!
   end subroutine add_nc_curr

   subroutine allocate_zero
      use ions_base, only: nsp
      use gvect, only: ngm
      !use uspp_data, ONLY: nqxq
      !USE uspp_param, ONLY: nbetam
!
      implicit none
!
      !integer ::isp

      !allocate (H_g(ngm, 3, 3, nsp))
      !allocate (tabr(nqxq, nbetam, nsp, 3))
!
   end subroutine allocate_zero

   subroutine deallocate_zero
!
      implicit none
!
!
   end subroutine deallocate_zero

   SUBROUTINE init_reciprocal_parts_tab(H_g, tablocal_hg, nsp, zv, tpiba2, tpiba, omega, at, alat, &
                                        ngm, gg, gstart, g, igtongl, gl, ngl, dq)
      use kinds, only: DP
      !use ions_base, only: nsp, zv
      !use cell_base, only: tpiba2, tpiba, omega, at, alat
      use constants, only: pi, fpi
      !use mp_world, only: mpime
      use mp, only: mp_sum
      !use gvect, only: ngm, gg, gstart, g, igtongl, gl, ngl
      !use uspp_data, only: dq
      use mp_pools, only: intra_pool_comm
!  use splinelib
      implicit none

      real(dp), intent(inout) :: H_g(:, :, :, :)
      integer, intent(in) :: nsp, igtongl(:), ngm, gstart, ngl
      real(dp), intent(in) :: zv(:), tpiba2, tpiba, omega, at(3, 3), alat, &
                              gg(:), g(:, :), gl(:), dq, tablocal_hg(:, :, :)

      real(DP) :: px, ux, vx, wx, xg
      integer  ::  a, b, igm, igl, it
      integer  :: i0, i1, i2, i3
      !real(DP) :: ra, rb
      real(DP) :: H_g_rad(ngl, 0:1)
      !logical  ::l_plot

!
!Initialization H_g
!
      H_g = 0.d0
      do it = 1, nsp
!initialization  radial part
         do igl = 1, ngl
            xg = sqrt(gl(igl))*tpiba
               px = xg/dq - int(xg/dq)
               ux = 1.d0 - px
               vx = 2.d0 - px
               wx = 3.d0 - px
               i0 = INT(xg/dq) + 1
               i1 = i0 + 1
               i2 = i0 + 2
               i3 = i0 + 3
! tablocal_hg is computed on a 1D in reciprocal space
! and is initialized in  init_us_1a.f90. It is independent of the cell and depends only on the local pseudo.
! H_g_rad is the the radial counterpart of H_g and is computed at the desired G values defined by the reciprocal
! lattice interpolating tablocal_hg
!
               H_g_rad(igl, 0) = tablocal_hg(i0, it, 1)*ux*vx*wx/6.d0 + &
                           &tablocal_hg(i1, it, 1)*px*vx*wx/2.d0 - &
                           &tablocal_hg(i2, it, 1)*px*ux*wx/2.d0 + &
                           &tablocal_hg(i3, it, 1)*px*ux*vx/6.d0
               H_g_rad(igl, 1) = tablocal_hg(i0, it, 2)*ux*vx*wx/6.d0 + &
                           &tablocal_hg(i1, it, 2)*px*vx*wx/2.d0 - &
                           &tablocal_hg(i2, it, 2)*px*ux*wx/2.d0 + &
                           &tablocal_hg(i3, it, 2)*px*ux*vx/6.d0
         end do
! Now we use H_g_rad to evaluate H_g, adding some "structure factors"
         do a = 1, 3
            do b = 1, 3
               if (a >= b) then
! Eq. 31. Add h1 to H_g
                  do igm = gstart, ngm
                     H_g(igm, a, b, it) = g(a, igm)*g(b, igm)/gg(igm)* &
                                        (sqrt(gg(igm))*tpiba*H_g_rad(igtongl(igm), 1) - 2.d0*zv(it)*fpi/omega*2.d0/(gg(igm)*tpiba2))
!
                  end do
!
               end if
! Add Eq. 32. Add h2 to H_g
               if (a == b) then
                  do igm = gstart, ngm
!
                     H_g(igm, a, b, it) = H_g(igm, a, b, it) - &
                                          (H_g_rad(igtongl(igm), 0) - fpi/omega*2.d0*zv(it)/(tpiba2*gg(igm)))
!
                  end do
! G=0 components from Eq. 32
                  if (gstart == 2) then
                     H_g(1, a, b, it) = -H_g_rad(1, 0)
                  end if
!
               end if
            end do
         end do
      end do

   END SUBROUTINE init_reciprocal_parts_tab

   subroutine init_us_1a(tabr, tablocal_hg, &
                         rgrid, nsp, zv, omega, nqxq, dq, upf)
      !----------------------------------------------------------------------
      !
      USE radial_grids, ONLY: radial_grid_type
      USE uspp_param, ONLY: pseudo_upf
      USE kinds, ONLY: DP
      USE constants, ONLY: fpi, sqrt2
      !USE atom, ONLY: rgrid
      !USE ions_base, ONLY: nsp, zv
      !USE cell_base, ONLY: omega
      !use uspp_data, ONLY: nqxq, dq
      USE splinelib
      !USE uspp_param, ONLY: upf
      USE mp_global, ONLY: intra_bgrp_comm
      USE mp, ONLY: mp_sum

      implicit none

      type(radial_grid_type), intent(in) :: rgrid(:)
      integer, intent(in) :: nsp, nqxq
      type(pseudo_upf), intent(in) :: upf(:)
      real(dp), intent(in) :: zv(:), omega, dq
      real(DP), intent(inout) :: tablocal_hg(:, :, :), tabr(:, :, :, :)
      !
      !
      !     here a few local variables
      !
      integer :: nt, nb, l, ir, iq, startq, &
                 lastq, ndm, rm
      ! various counters
      real(DP), allocatable :: aux(:), besr(:)
      ! various work space
      real(DP) ::  pref, qi
      ! the prefactor of the q functions
      ! the prefactor of the beta functions
      ! the modulus of g for each shell
      ! q-point grid for interpolation
      real(DP) ::  vqint
      ! interpolated value
      !integer ::  &
      !             vk
      !
      !real(DP), allocatable :: xdata(:)

      !
      call start_clock('init_us_1a')
      !    Initialization of the variables
      !

!!!!!!! init_us_1a initializes tabr & tabocal
! 1. tabr(nqxq,nbetam,nsp,1:3)
! Contains the integrals f_(la,lb)(q)=\int _0 ^\infty dr r^3 f_la(r) j_lb(q.r) needed in Eq. 41
! - la is the angolar momentum of the beta function (more than one function
! with the same l are allowed )
! - lb is the angular moment of bessel function
! - only lb=la+1, lb=la-1 or lb=la are allowed
!   these 3 cases are represented by the last index
!
! 2. tablocal(nqxq,nsp,0:2)
! The radial integrals needed in Eq. 31 , 32. They depend on the modulus of G (nqxq), on the atom (nsp)
! Third index refers to the bessel function with l=0 (0) or l=1 (1)

      ndm = MAXVAL(upf(:)%kkbeta)
      allocate (aux(ndm))
      allocate (besr(ndm))
!
      pref = fpi/sqrt(omega)
      call divide(intra_bgrp_comm, nqxq, startq, lastq)
      tabr(:, :, :, :) = 0.d0
      do nt = 1, nsp
         do nb = 1, upf(nt)%nbeta
            l = upf(nt)%lll(nb)
            do iq = startq, lastq
               qi = (iq - 1)*dq
!initializzation of table with l_bessel=l
!kkbeta tells where bessel function goes to zero
               call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
               do ir = 1, upf(nt)%kkbeta
                  aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
                  ! upf files contains x times the projector
                  ! that is why  we have only two rgrids r^2)
               enddo
               call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
               tabr(iq, nb, nt, 2) = vqint*pref
!l_bessel=l+1
               call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l + 1, besr)
               do ir = 1, upf(nt)%kkbeta
                  aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
               enddo
               call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
               tabr(iq, nb, nt, 3) = vqint*pref
!l_bessel=l-1, only if l>0
               if (l > 0) then
                  call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l - 1, besr)
                  do ir = 1, upf(nt)%kkbeta
                     aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
                  enddo
                  call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
                  tabr(iq, nb, nt, 1) = vqint*pref
               end if
            enddo
         enddo
      enddo
#ifdef __MPI
      call mp_sum(tabr, intra_bgrp_comm)
#endif
!!!!!!!!!!!!!!
!Initialization of tablocal_hg
!Warning the new grid must be compatible with the one of the local pseudo
      rm = MAXVAL(rgrid(:)%mesh)
      deallocate (besr)
      allocate (besr(rm))
      deallocate (aux)
      allocate (aux(rm))
!
      pref = fpi/omega
      call divide(intra_bgrp_comm, nqxq, startq, lastq)
      tablocal_hg(:, :, :) = 0.d0
      do nt = 1, nsp
         do iq = startq, lastq
            qi = (iq - 1)*dq
!initialization of l=0
            call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 0, besr)
            do ir = 1, rgrid(nt)%mesh
               aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                         besr(ir)*rgrid(nt)%r(ir)
            end do
            call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
            tablocal_hg(iq, nt, 1) = vqint*pref
!initialization of l=1
            call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 1, besr)
            do ir = 1, rgrid(nt)%mesh
               aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                         besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
            end do
            call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
            tablocal_hg(iq, nt, 2) = vqint*pref
         end do
      end do
#ifdef __MPI
      call mp_sum(tablocal_hg, intra_bgrp_comm)
#endif
      ! initialize spline interpolation
      deallocate (besr)
      deallocate (aux)

      call stop_clock('init_us_1a')
      return
   end subroutine init_us_1a

END MODULE zero_mod

