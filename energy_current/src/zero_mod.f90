MODULE zero_mod

   USE kinds, ONLY: DP
   use splines
   use fft_base
   use becmod

   SAVE
!non local potential variables
   logical  ::l_non_loc, l_tab
   real(DP), allocatable :: tabr(:, :, :, :)
   real(DP), allocatable :: tabr_d2y(:, :, :, :)
   real(DP), allocatable :: tablocal_hg(:, :, :)
   real(DP), allocatable :: tablocal_d2y_hg(:, :, :)
   TYPE(bec_type) :: becpr(3)
   TYPE(bec_type) :: becpd(3)
   TYPE(bec_type) :: becprd(3, 3)
   complex(DP), allocatable :: xvkb(:, :, :)
   complex(DP), allocatable :: dvkb(:, :, :)
   complex(DP), allocatable :: xdvkb(:, :, :, :)

   !the component of the current here computed
   real(dp) ::z_current(3), i_current(3), i_current_a(3), i_current_b(3), i_current_c(3), i_current_d(3), i_current_e(3)

   !sum or not i_current_b in i_current
   logical :: add_i_current_b = .false.

   character(len=256) :: vel_input_units = 'PW'

   !input from stdout
   integer        :: natoms !cutoff per somme in griglia reale
   integer        :: n_max !cutoff per somme in griglia reale
   real(kind=DP)  :: eta !ewald factor for convergence

   real(DP) ::I_primo_rec
   real(DP), allocatable   :: H_g(:, :, :, :)
   complex(DP), allocatable   :: u_g(:, :) !ngm,a  !this looks like can be transformed in a local variable
   complex(DP), allocatable ::charge_g(:)
   real(DP), allocatable   :: I_uno_g(:, :, :) !griglia,a,b
   real(DP), allocatable   :: I_due_g(:)
   real(DP)  ::I_primo

   !variable for S_due
   real(DP), allocatable   ::alpha_0_all(:, :), alpha_2_all(:, :)

   !counters and indexes
   integer      ::nr3s_start, nr3s_end !z-planes in processor
   integer      ::nr3s_start_par
contains

   real(kind=DP) function h_(x)
      real(kind=DP) :: x
      real(kind=DP), external :: qe_erfc
      h_ = qe_erfc(x)/x
   end function h_

   real(kind=DP) function hp_(x)
      use constants, only: pi
      real(kind=DP) :: x
      real(kind=DP), external :: qe_erfc
      hp_ = -(2.d0/sqrt(pi))*(1.d0/x)*exp(-x*x) - 1.d0/(x*x)*qe_erfc(x)
   end function hp_

   real(kind=DP) function modulus(vector)
      real(kind=DP) ::vector(3)
      modulus = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)
   end function modulus

   subroutine pbc(vin, vout, alat, at, bg)
      implicit none
      real(DP), intent(in) :: vin(3), alat, at(3,3), bg(3,3)
      real(DP), intent(out) :: vout(3)
      !local
      real(DP):: alatdir ! box side in a given direction
      integer :: i, n
      !vout used also as help variable
      vout(:) = matmul(vin(:)/alat, bg(:, :)) ! vout = vin in S space
      vout(:) = vout(:) - ANINT(vout(:))
      vout(:) = MATMUL(at(:, :), vout(:))*alat ! back from S space to R space
   end subroutine pbc

   subroutine I_uno_value(y, x, a, b, flag, tpiba, at, alat, gstart, g, gg, npw)
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::flag, npw, gstart
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3), tpiba, g(:,:), gg(:), at(3,3),alat
      integer, intent(in)   ::a, b
      integer ::igm
      real(DP) :: comp_iso

      call start_clock('rec')

      y = 0.d0
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0) exit
         y = y + 2.d0*(I_uno_g(igm, a, b)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
      end do
      if (gstart == 2) y = y + I_uno_g(1, a, b)
      call mp_sum(y, intra_pool_comm)
      call stop_clock('rec')
      if (flag == 1) then
         call add_local_uno(y, x, a, b, at, alat)
      end if
!change for opt
!  call I_due_value(comp_iso,x,1)
!  if (a==b) then
!     y=y-comp_iso
!  end if
      !if (a == b) then
      !   call I_due_value(comp_iso, x, 1)
      !   y = y - comp_iso
      !end if

   end subroutine I_uno_value

   subroutine add_local_uno(value, pos, a, b, at, alat)
      use mp_world, only: nproc, mpime
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::a, b
      real(DP), intent(in) ::pos(3), at(3,3), alat
      real(DP), intent(inout) ::value
      real(DP) ::u(3), u_mod, n(3)
      integer  ::n_x, n_y, n_z
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1

      call start_clock('real')
      l_blk = (2*n_max + 1)/nproc
      if (l_blk*nproc < (2*n_max + 1)) l_blk = l_blk + 1
      nbegin = mpime*l_blk - n_max
      nend = nbegin + l_blk - 1
      if (nend > n_max) nend = n_max
      value1 = 0.d0
      do n_x = nbegin, nend
         do n_y = -n_max, n_max
            do n_z = -n_max, n_max
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               u(1:3) = pos(1:3) - n(1:3)
               u_mod = modulus(u)
               value1 = value1 + eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
               if (a == b) then
                  value1 = value1 + sqrt(eta)*h_(sqrt(eta)*u_mod)
               end if
            end do
         end do
      end do
      call mp_sum(value1, intra_pool_comm)
      value = value + value1
      call stop_clock('real')
!!!!!!!!!!!!!!!!!!!
      !call start_clock( 'real' )
      !do n_x=-n_max,n_max
      !   do n_y=-n_max,n_max
      !      do n_z=-n_max,n_max
      !         n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
      !         u(1:3) = pos(1:3)-n(1:3)
      !         u_mod=modulus(u)
      !         value=value+eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
      !         if (a==b) then
      !            value=value+sqrt(eta)*h_(sqrt(eta)*u_mod)
      !         end if
      !      end do
      !   end do
      !end do
      !call stop_clock( 'real' )
   end subroutine add_local_uno

   subroutine I_due_value(y, x, flag, tpiba, at, alat, gstart, g, gg, npw)
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::flag, npw, gstart
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3), at(3,3), alat,g(:,:), gg(:), tpiba
      integer ::igm!,ng_max
      real(DP) ::scalar

      call start_clock('rec')
      y = 0.d0
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0) exit
         y = y + 2.d0*(I_due_g(igm)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
      end do
      if (gstart == 2) y = y + I_due_g(1)
      call mp_sum(y, intra_pool_comm)
      if (flag == 1) then
         call add_local_due(y, x, at, alat)
      end if
      call stop_clock('rec')
   end subroutine i_due_value

   subroutine add_local_due(value, pos, at, alat)
      !This routines computes sum_{L}erfc(sqrt(eta)abs(d-L))/abs(d-L)
      ! the sum is in real space over a small number of shells
      use mp, only: mp_sum
      use mp_world, ONLY: nproc, mpime
      use mp_pools, ONLY: intra_pool_comm
      implicit none
      real(DP), intent(in) ::pos(3), at(3,3), alat
      real(DP), intent(inout) ::value
      real(DP) ::modul, n(3), erf_value
      integer  ::n_x, n_y, n_z
      real(kind=DP), external :: qe_erfc
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1

      call start_clock('real')
      l_blk = (2*n_max + 1)/nproc
      if (l_blk*nproc < (2*n_max + 1)) l_blk = l_blk + 1
      nbegin = mpime*l_blk - n_max
      nend = nbegin + l_blk - 1
      if (nend > n_max) nend = n_max
      value1 = 0.d0
      do n_x = nbegin, nend
         do n_y = -n_max, n_max
            do n_z = -n_max, n_max
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               modul = modulus(pos(1:3) - n(1:3))
               erf_value = qe_erfc(sqrt(eta)*modul)
               value1 = value1 + erf_value/modul
            end do
         end do
      end do
      call mp_sum(value1, intra_pool_comm)
      value = value + value1
      call stop_clock('real')
!!!!!!!!!!
!  call start_clock( 'real' )
!  do n_x=-n_max,n_max         !primo ciclo   su n
!     do n_y=-n_max,n_max      !secondo ciclo su n
!        do n_z=-n_max,n_max   !terzo ciclo   su n
!           n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
!           modul=modulus(pos(1:3)-n(1:3))
!           erf_value=qe_erfc(sqrt(eta)*modul)
!           value=value+ erf_value/modul
!        end do
!     end do
!  end do
!  call stop_clock( 'real' )

   end subroutine add_local_due

!!!!!!!!!!!!!!!!!!!



subroutine init_zero(nsp, zv, tpiba2, tpiba, omega, at, alat, &
                ngm, gg, gstart, g, igtongl, gl, ngl, spline_ps, dq, &
                upf, rgrid, nqxq)
!called once to init stuff that does not depend on the atomic positions
   use kinds, only : DP
   USE radial_grids, ONLY : radial_grid_type
   USE uspp_param, ONLY: pseudo_upf
   !use gvect, only: ngm, gg, g, gstart
   !use hartree_mod
   !use ions_base, only: nsp
   implicit none

   integer, intent(in) :: nsp, igtongl(:), ngm, gstart, ngl, nqxq
   real(dp), intent(in) :: zv(:), tpiba2, tpiba, omega, at(3,3), alat, &
                       gg(:), g(:,:), gl(:), dq
   logical, intent(in) :: spline_ps
   type(radial_grid_type), intent(in) :: rgrid(:)
   type(pseudo_upf), intent(in) :: upf(:)

   integer :: isp, iun, a, b
   character(256) :: pref_box
   integer, external :: find_free_unit
   logical :: exst
   call start_clock('init_zero')
   !call init_us_1all(rgrid, nsp, zv, omega, nqxq, dq, spline_ps, upf)
   call init_us_1()
   call init_us_1a(rgrid, nsp, zv, omega, nqxq, dq, spline_ps, upf)
   call init_reciprocal_parts_tab(nsp, zv, tpiba2, tpiba, omega, at, alat, &
           ngm, gg, gstart, g, igtongl, gl, ngl, spline_ps, dq)
   do a = 1, 3
      do b = 1, 3
         if (a > b) then
            do isp = 1, nsp
               H_g(:, b, a, isp) = H_g(:, a, b, isp)
            end do
            I_uno_g(:, b, a) = I_uno_g(:, a, b)
         end if
      end do
   end do
   call stop_clock('init_zero')

end subroutine


subroutine current_zero(nbnd, npwx, npw, dffts, nsp, zv, nat, ityp, amass, tau, &
                        vel, tpiba, tpiba2, at, alat, omega, psic, evc, ngm, gg, g, gstart, &
                        nkb, vkb, deeq, upf, nh, xk, igk_k, bg )
   use kinds, only: DP
   use fft_base, only: fft_type_descriptor
   use mp, only: mp_sum
   use mp_pools, only: intra_pool_comm
   use io_global, only: stdout, ionode
   use constants, only: e2, AMU_RY
   use splines
   use hartree_mod
   use compute_charge_mod, only : compute_charge
   USE uspp_param, ONLY: pseudo_upf

   implicit none

   INTEGER, intent(in) :: nbnd, npwx, nsp, nat, ityp(:), &
                           ngm, gstart, nkb, igk_k(:,:), nh(:)
   INTEGER, intent(inout) :: npw
   TYPE ( fft_type_descriptor ), intent(inout) :: dffts
   REAL(DP), intent(in) :: zv(:), amass(:), tau(:,:), vel(:,:), &
                           tpiba, omega, tpiba2, alat, at(:,:), &
                           g(:,:), gg(:), bg(:,:)
   COMPLEX(DP), intent(inout) :: psic(:), evc(:,:)
   COMPLEX(DP), intent(in):: vkb(:,:)
   real(dp), intent(in) :: deeq(:,:,:,:), xk(:,:)
   type(pseudo_upf), intent(in) :: upf(:)

!three dimensional auxiliary real vectors
   real(dp), parameter  :: amconv = AMU_RY
   real(DP) ::n(3), u(3), u_pbc(3), u_x(3), u_y(3), u_z(3), value, x(1:3), ics(3)
   real(DP), external :: qe_erfc
   !real(DP)  :: dtauij(3)
   real(DP) ::modul, erf_value
   real(DP) ::fac_uno, fac_due, fac_tre, fac_quattro
   real(DP), external :: qe_erf
   integer, external  :: find_free_unit

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
   character(256) ::filename, pref_box
   real(DP), allocatable :: values(:)
   real(DP), allocatable   :: charge(:)

   if (ionode) write (stdout, *) 'ROUTINE_ZERO BEGINNING'
   call start_clock('routine_zero')
   if (ionode) print *, 'eta', eta
   allocate (charge(dffts%nnr))
   call start_clock('zero_current')
   if (nkb > 0) then
      l_non_loc = .true.
   else
      l_non_loc = .false.
   end if
   npw=npwx
   call compute_charge(psic, evc, npw, nbnd, ngm, dffts, charge, charge_g)
!
!initialization of  u_g
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
!computation of the current
   z_current = 0.d0
   do a = 1, 3
      do igm = gstart, ngm
         z_current(a) = z_current(a) + 2.d0*dble(charge_g(igm)*conjg(u_g(igm, a)))
      end do
      if (gstart == 2) then
         z_current(a) = z_current(a) + dble(charge_g(1)*conjg(u_g(1, a)))
      end if
   end do
   call mp_sum(z_current, intra_pool_comm)
   if (l_non_loc) then
      call add_nc_curr(z_current, nkb, vkb, deeq, upf, nh, vel, nbnd, npw, npwx, evc, &
                       g, tpiba, nat, ityp, nsp, xk, igk_k)
   end if
   z_current = z_current*alat
   call stop_clock('zero_current')
   call print_clock('zero_current')
   deallocate (charge)
   if (ionode) print *, 'ZERO CURRENT CALCULATED'
   call stop_clock('routine_zero')
   call print_clock('routine_zero')

end subroutine current_zero


subroutine current_ionic(current, current_a, current_b, current_c, current_d, current_e, add_current_b, &
                         nat, tau, vel, zv, ityp, alat, at, bg, tpiba, gstart, g, gg, npw, amass)
   use kinds, only: DP
   use constants, only: AMU_RY, e2
   use io_global, only: ionode

   implicit none

   real(dp), intent(out) :: current(3), current_a(3), current_b(3), current_c(3),&
                            current_d(3), current_e(3)
   logical, intent(in) :: add_current_b
   INTEGER, intent(in) :: npw, nat, ityp(:), gstart
   REAL(DP), intent(in) :: zv(:), tau(:,:), vel(:,:), &
                           tpiba, alat, at(:,:), &
                           g(:,:), gg(:), bg(:,:), amass(:)

   real(dp) :: u(3), u_pbc(3), val
   integer :: iatom, jatom, a, b

   current = 0.d0
   current_a = 0.d0
   current_b = 0.d0
   current_c = 0.d0
   current_d = 0.d0
   current_e = 0.d0

   call start_clock('calcolo_i')
   do iatom = 1, nat
      current_a(:) = current_a(:) + vel(:, iatom)*(1./2.*amu_ry*amass(ityp(iatom))*(vel(1, iatom)**2 +&
                       vel(2, iatom)**2 + vel(3, iatom)**2))
      current_b(:) = current_b(:) + 2./3.*e2*zv(ityp(iatom))**2*vel(:, iatom)*I_primo
   end do

   do iatom = 1, nat
      do jatom = 1, nat
         if (iatom > jatom) then
            u(1:3) = (tau(:, iatom) - tau(:, jatom))*alat
            call pbc(u(1:3), u_pbc(1:3), alat, at, bg)
            call I_due_value(val, u_pbc, 1, tpiba, at, alat, gstart, g, gg, npw)
            !current_c(:) = current_c(:) + 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
            current_c(:) = current_c(:) + 1.d0*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
                             (vel(:, iatom) + vel(:, jatom))*val

            do a = 1, 3
               do b = 1, 3
                  if (a > b) then
                     call I_uno_value(val, u_pbc, a, b, 1, tpiba, at, alat, gstart, g, gg, npw)
                     current_e(a) = current_e(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*val
                     current_e(b) = current_e(b) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(a, jatom) + vel(a, iatom))*val
                  end if
                  if (a == b) then
                     call I_uno_value(val, u_pbc, a, b, 1, tpiba, at, alat, gstart, g, gg, npw)
                     current_d(a) = current_d(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*val
                  end if
               end do
            end do
         end if
      end do
   end do
   current_a = current_a * alat**3
   current_b = current_b * alat
   current_c = current_c * alat
   current_d = current_d * alat
   current_e = current_e * alat
   current   = current_a + current_c + current_d + &
                 current_e
   if (add_current_b) &
      current = current + current_b


   call stop_clock('calcolo_i')
   call print_clock('calcolo_i')
   if (ionode) print *, 'IONIC CURRENT CALCULATED'

   end subroutine

subroutine add_nc_curr(current, nkb, vkb, deeq, upf, nh, vel, nbnd, npw, npwx, evc, &
                       g, tpiba, nat, ityp, nsp, xk, igk_k)
   use kinds, only: DP
   use becmod
   USE uspp_param, ONLY: pseudo_upf

   implicit none

   real(DP), intent(inout) :: current(3)
   integer, intent(in) :: nkb, nh(:), nbnd, npw, npwx, nat, ityp(:), &
                          nsp, igk_k(:,:)
   complex(DP), intent(in) :: vkb(:,:), evc(:,:)
   real(dp), intent(in) :: deeq(:,:,:,:), vel(:,:), g(:,:), tpiba, xk(:,:)
   type(pseudo_upf), intent(in) :: upf(:)

   !integer ::iun, a, b
   !logical ::l_test
   integer :: ipol, jpol, ibnd
   integer ::ijkb, ikb, ih, na, nt, ipw
   real(DP) ::J_nl(3), J_1(3), J_2(3)
   complex(DP), allocatable ::vkb1(:, :)
   integer, external       :: find_free_unit

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

!inizializzazione di tab(serve?),tabr ed indici
!   call init_us_1a()
!inizializzazione di vkb (per essere sicuri che lo sia) e xvkb
   CALL init_us_2(npw, igk_k(1, 1), xk(1, 1), vkb)
   call init_us_3(npw, xvkb)
!!
!inizializzazione di dvkb e xdvkb (nb il ciclo su ipw va dentro per essere
!ottimizzato)
!servono altre inizializzazioni?
   dvkb = 0.d0
   do ipol = 1, 3
      do ikb = 1, nkb
         do ipw = 1, npw
            dvkb(ipw, ikb, ipol) = (0.d0, -1.d0)*tpiba*g(ipol, ipw)*vkb(ipw, ikb)
!         dvkb(ipw,ikb,ipol)=(0.d0,-1.d0)*g(ipol,ipw)*vkb(ipw,ikb)
         end do
      end do
   end do
   do jpol = 1, 3
      do ipol = 1, 3
         do ikb = 1, nkb
            do ipw = 1, npw

               xdvkb(ipw, ikb, ipol, jpol) = (0.d0, -1.d0)*tpiba*g(jpol, ipw)*xvkb(ipw, ikb, ipol)
!            xdvkb(ipw,ikb,ipol,jpol)=(0.d0,-1.d0)*g(jpol,ipw)*xvkb(ipw,ikb,ipol)
               if (ipol == jpol) then
                  xdvkb(ipw, ikb, ipol, jpol) = xdvkb(ipw, ikb, ipol, jpol) + vkb(ipw, ikb)

               end if
            end do
         end do
      end do
   end do
!
!
!prodotti scalari (si possono evitare di fare tutti?), qui si esegue comunicazione MPI.
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

!
!a questo punto possiamo usare le quantità calcolate per calcolare la corrente (OpenMP do?)
   J_nl = 0.d0
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
!il fattore due è fatto per la degenrazione di spin
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
   USE us, ONLY: nqxq, spline_ps
   USE uspp_param, ONLY: nbetam
!
   implicit none
!
   !integer ::isp

   allocate (H_g(ngm, 3, 3, nsp))
   allocate (I_uno_g(ngm, 3, 3))
   allocate (tablocal_hg(nqxq, nsp, 0:1))
   if (spline_ps) then
      allocate (tablocal_d2y_hg(nqxq, nsp, 0:1))
   end if
   allocate (I_due_g(ngm))
   allocate (tabr(nqxq, nbetam, nsp, -1:1))
   if (spline_ps) then
      allocate (tabr_d2y(nqxq, nbetam, nsp, -1:1))
   end if
   allocate (charge_g(ngm))
   allocate (u_g(ngm, 3))
!
end subroutine allocate_zero


subroutine deallocate_zero
!
   implicit none
!
   if (allocated(charge_g)) &
      deallocate (charge_g)
   if (allocated(u_g)) &
      deallocate (u_g)
!
end subroutine deallocate_zero

   real(DP) function alpha_0_lr(eta, G_sq, flag)
      real(DP)  ::eta, G_sq
      integer ::flag
      if (flag == 1) then
         alpha_0_lr = exp(-G_sq/(4.d0*eta))/G_sq
      else
         alpha_0_lr = 1.d0/G_sq
      end if
   end function alpha_0_lr


SUBROUTINE init_reciprocal_parts_tab(nsp, zv, tpiba2, tpiba, omega, at, alat, &
                ngm, gg, gstart, g, igtongl, gl, ngl, spline_ps, dq)
   use kinds, only: DP
   !use ions_base, only: nsp, zv
   !use cell_base, only: tpiba2, tpiba, omega, at, alat
   use constants, only: pi, fpi
   !use mp_world, only: mpime
   use mp, only: mp_sum
   !use gvect, only: ngm, gg, gstart, g, igtongl, gl, ngl
   !use us, only: spline_ps, dq
   use mp_pools, only: intra_pool_comm
!  use splinelib
   implicit none

   integer, intent(in) :: nsp, igtongl(:), ngm, gstart, ngl
   real(dp), intent(in) :: zv(:), tpiba2, tpiba, omega, at(3,3), alat, &
                       gg(:), g(:,:), gl(:), dq
   logical, intent(in) :: spline_ps

   real(DP) :: px, ux, vx, wx, xg
   real(DP) :: n(3), modul, erf_value
   integer  ::  a, b, igm, igl, it, n_x, n_y, n_z
   integer  :: i0, i1, i2, i3
   !real(DP) :: ra, rb
   real(DP) :: H_g_rad(ngl, 0:1)
   !logical  ::l_plot
   real(kind=DP), external :: qe_erfc

! I_primo == S_{2,B} (Aris note)
   I_primo = 0.d0
! these sum over n_x n_y n_z are sum_{L!=0}erfc(sqrt(eta*L)*1/L)
   do n_x = -n_max, n_max
      do n_y = -n_max, n_max
         do n_z = -n_max, n_max
            if ((n_x /= 0) .or. (n_y /= 0) .or. (n_z /= 0)) then
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               modul = modulus(n(:))
               erf_value = qe_erfc(sqrt(eta)*modul)
               I_primo = I_primo + erf_value/modul
            end if
         end do
      end do
   end do
!now compute the other addendum of S_{2,B}
   I_primo = I_primo - 2.d0*sqrt(eta/pi)
!
   I_primo_rec = 0.d0
   do igm = gstart, ngm
      I_primo_rec = I_primo_rec + 2.d0*(4.d0*pi)/omega*alpha_0_lr(eta, gg(igm)*tpiba2, 1)
   end do
   call mp_sum(I_primo_rec, intra_pool_comm)
   I_primo_rec = I_primo_rec - pi/(eta*omega)
!
   I_primo = I_primo + I_primo_rec
!
!calcolo di I_due_g
!
   do igm = gstart, ngm
!
      I_due_g(igm) = (4.d0*pi)/omega*alpha_0_lr(eta, gg(igm)*tpiba2, 1)
!
   end do
!
   if (gstart == 2) then
      I_due_g(1) = -pi/(eta*omega)
   end if
!
!analogamente a prima inizializzo I_uno_g
   do a = 1, 3
      do b = 1, 3
         if (a >= b) then
            do igm = gstart, ngm
               I_uno_g(igm, a, b) = (4.d0*pi)/omega*exp(-tpiba2*gg(igm)/(4.d0*eta))/(tpiba2*gg(igm))*&
 & (2 + gg(igm)*tpiba2/(2.d0*eta))*g(a, igm)*g(b, igm)/gg(igm)
            end do
!
            if (gstart == 2) then
               I_uno_g(1, a, b) = 0.d0
            end if
         end if
      end do
   end do
!
!Inizializzo H_g
!
   H_g = 0.d0
   do it = 1, nsp
!inizializzazione parti radiali
      do igl = 1, ngl
         xg = sqrt(gl(igl))*tpiba
         if (spline_ps) then
            CALL errore('init_reciprocal_parts', 'splines not implemented', 1)
         else
            px = xg/dq - int(xg/dq)
            ux = 1.d0 - px
            vx = 2.d0 - px
            wx = 3.d0 - px
            i0 = INT(xg/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
!tablocal_hg e' probabilmente fatto su una griglia 1D ed e' inizializzato in init_us_1a.f90, e' indipendente dalla cell
!H_g_rad serve per fare un integrale radiale, interpolando tablocal_hg
            H_g_rad(igl, 0) = tablocal_hg(i0, it, 0)*ux*vx*wx/6.d0 + &
                        &tablocal_hg(i1, it, 0)*px*vx*wx/2.d0 - &
                        &tablocal_hg(i2, it, 0)*px*ux*wx/2.d0 + &
                        &tablocal_hg(i3, it, 0)*px*ux*vx/6.d0
            H_g_rad(igl, 1) = tablocal_hg(i0, it, 1)*ux*vx*wx/6.d0 + &
                        &tablocal_hg(i1, it, 1)*px*vx*wx/2.d0 - &
                        &tablocal_hg(i2, it, 1)*px*ux*wx/2.d0 + &
                        &tablocal_hg(i3, it, 1)*px*ux*vx/6.d0
         end if
      end do
      do a = 1, 3
         do b = 1, 3
            if (a >= b) then
!
               do igm = gstart, ngm
                  H_g(igm, a, b, it) = g(a, igm)*g(b, igm)/gg(igm)* &
                                       (sqrt(gg(igm))*tpiba*H_g_rad(igtongl(igm), 1) - 2.d0*zv(it)*fpi/omega*2.d0/(gg(igm)*tpiba2))
!
               end do
!
            end if
            if (a == b) then
!
               do igm = gstart, ngm
!
                  H_g(igm, a, b, it) = H_g(igm, a, b, it) - &
                                       (H_g_rad(igtongl(igm), 0) - fpi/omega*2.d0*zv(it)/(tpiba2*gg(igm)))
!
               end do
               if (gstart == 2) then
                  H_g(1, a, b, it) = -H_g_rad(1, 0)
               end if
!
            end if
         end do
      end do
   end do

END SUBROUTINE init_reciprocal_parts_tab


subroutine init_us_1a(rgrid, nsp, zv, omega, nqxq, dq, spline_ps, upf)
   !----------------------------------------------------------------------
   !
   USE radial_grids, ONLY : radial_grid_type
   USE uspp_param, ONLY: pseudo_upf
   USE kinds, ONLY: DP
   USE constants, ONLY: fpi, sqrt2
   !USE atom, ONLY: rgrid
   !USE ions_base, ONLY: nsp, zv
   !USE cell_base, ONLY: omega
   !USE us, ONLY: nqxq, dq, spline_ps
   USE splinelib
   !USE uspp_param, ONLY: upf
   USE mp_global, ONLY: intra_bgrp_comm
   USE mp, ONLY: mp_sum

   implicit none

   type(radial_grid_type), intent(in) :: rgrid(:)
   integer, intent(in) :: nsp, nqxq
   type(pseudo_upf), intent(in) :: upf(:)
   real(dp), intent(in) :: zv(:), omega, dq
   logical, intent(in) :: spline_ps

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

!!!!!!!aggiunta di codice per la zero_current
!tabr(nqxq,nbetam,nsp,-1:1) e tablocal(nqxq,nsp,0:2)
!contiene f_(la,lb)(q)=\int _0 ^\infty dr r^3 f_la(r) j_lb(q.r)
!dove la è il momento angolare della beta function (ce ne possono
!essere più di una con il medesimo l) e
!lb è il momento angolare della funzione di bessel
!si può avere solo lb=la+1 lb=la-1 o lb=la
!e questi tre casi sono rappresentati dall'ultimo indice
!3.71
   ndm = MAXVAL(upf(:)%kkbeta)
   allocate (aux(ndm))
   allocate (besr(ndm))

!
!inizializzazione di tabr
! questo può essere estratto. Il resto della routine è init_us_1.f90 di PW
!
   pref = fpi/sqrt(omega)
   call divide(intra_bgrp_comm, nqxq, startq, lastq)
   tabr(:, :, :, :) = 0.d0
   do nt = 1, nsp
      do nb = 1, upf(nt)%nbeta
         l = upf(nt)%lll(nb)
         do iq = startq, lastq
            qi = (iq - 1)*dq
!inizializzazione tabella con l_bessel=l
            !                  kkbeta dice dove bessel va a zero
            call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
            do ir = 1, upf(nt)%kkbeta
               aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir) ! nel file upf c'è x per il proiettorie
               !quindi qui c'è solo r^2)
            enddo
            call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
            tabr(iq, nb, nt, 0) = vqint*pref
!l_bessel=l+1
            call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l + 1, besr)
            do ir = 1, upf(nt)%kkbeta
               aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
            enddo
            call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
            tabr(iq, nb, nt, 1) = vqint*pref
!l_bessel=l-1, solo se l>0
            if (l > 0) then
               call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l - 1, besr)
               do ir = 1, upf(nt)%kkbeta
                  aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
               enddo
               call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
               tabr(iq, nb, nt, -1) = vqint*pref
            end if
         enddo
      enddo
   enddo
#ifdef __MPI
   call mp_sum(tabr, intra_bgrp_comm)
#endif
   ! initialize spline interpolation
   if (spline_ps) then
      CALL errore('init_us_1a', 'splines for tabr not implemented', 1)
   endif
!!!!!!!!!!!!!!
!inizializzazione di tablocal_hg 3.17 3.20
   !attenzione alla nuova griglia che deve essere compatibile con quella
!dello pseudo locale
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
!inizializzazione con l=0
         call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 0, besr)
         do ir = 1, rgrid(nt)%mesh
            aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                      besr(ir)*rgrid(nt)%r(ir)
         end do
         call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
         tablocal_hg(iq, nt, 0) = vqint*pref
!inizializzazione con l=1
         call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 1, besr)
         do ir = 1, rgrid(nt)%mesh
            aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                      besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
         end do
         call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
         tablocal_hg(iq, nt, 1) = vqint*pref
      end do
   end do
#ifdef __MPI
   call mp_sum(tablocal_hg, intra_bgrp_comm)
#endif
   ! initialize spline interpolation
   if (spline_ps) then
      CALL errore('init_us_1a', 'splines for tablocal not implemented', 1)
   endif
   deallocate (besr)
   deallocate (aux)

   call stop_clock('init_us_1a')
   return
end subroutine init_us_1a


END MODULE zero_mod

