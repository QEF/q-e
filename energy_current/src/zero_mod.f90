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

   !variables depending on the step
   !wavefunction
   !complex(DP),allocatable :: evc_uno(:,:)

   !ion positions and velocities
   !real(DP), allocatable ::ion_pos(:,:)
   !real(DP), allocatable ::ion_vel(:, :)
   !second ion positions and velocities read from input
   !real(DP), allocatable ::ion_pos2(:,:) ! must call convert_tau from ../PW/src/input.f90 to obtain correct units for positions
   character(len=256) :: vel_input_units = 'PW'

   !input from stdout
   integer        :: natoms !cutoff per somme in griglia reale
   integer        :: n_max !cutoff per somme in griglia reale
!  character(len=256) :: status !what to do with the program. "initialize" or "compute"
   real(kind=DP)  :: eta !ewald factor for convergence

   real(DP) ::I_primo_rec
   real(DP), allocatable   :: H_g(:, :, :, :)
   complex(DP), allocatable   :: u_g(:, :) !ngm,a
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!

   real(kind=DP) function modulus(vector)
      real(kind=DP) ::vector(3)
      modulus = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)
   end function modulus

!!!!!!!!!!!!!!!!!!!!!!!!

   real(DP) function alpha_0_lr(eta, G_sq, flag)
      real(DP)  ::eta, G_sq
      integer ::flag
      if (flag == 1) then
         alpha_0_lr = exp(-G_sq/(4.d0*eta))/G_sq
      else
         alpha_0_lr = 1.d0/G_sq
      end if
   end function alpha_0_lr

!!!!!!!!!!!!!!!!!!!!!!!

   real(DP) function alpha_2_lr(eta, G_sq, flag)
      use constants, only: pi
      real(kind=DP), external :: qe_erf
      real(DP)  ::eta, G_sq
      integer ::flag
      if (flag == 1) then
         alpha_2_lr = 1.d0/G_sq*((3.d0*sqrt(pi)*sqrt(eta/G_sq))*qe_erf(sqrt(G_sq/(eta))/2.d0) - exp(-(G_sq)/(4.d0*eta)))
      else
         alpha_2_lr = 2.d0/G_sq
      end if
   end function alpha_2_lr

!!!!!!!!!!!!

   subroutine mid(vin, vout)
      use kinds, only: DP
      use cell_base, only: alat
! apply the minimum image distance in a cubic cell
      implicit none
      real(DP), intent(in) :: vin(3)
      real(DP), intent(out) :: vout(3)
      !local
      integer :: i
      do i = 1, 3
         vout(i) = vin(i) - nint(vin(i)/alat)*alat
      end do
   end subroutine mid

!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine pbc_ortho(vin, vout)
      use kinds, only: DP
      use cell_base, only: alat, at, bg
! apply the minimum image distance in a orthogonal cell
      implicit none
      real(DP), intent(in) :: vin(3)
      real(DP), intent(out) :: vout(3)
      !local
      real(DP):: alatdir ! box side in a given direction
      integer :: i, n
      !vout used also as help variable
      vout(:) = matmul(vin(:)/alat, bg(:, :)) ! vout = vin in S space
      vout(:) = vout(:) - ANINT(vout(:))
      vout(:) = MATMUL(at(:, :), vout(:))*alat ! back from S space to R space
   end subroutine pbc_ortho

   subroutine pbc_ortho_old(vin, vout)
      use kinds, only: DP
      use cell_base, only: alat, at
! apply the minimum image distance in a orthogonal cell
      implicit none
      real(DP), intent(in) :: vin(3)
      real(DP), intent(out) :: vout(3)
      !local
      real(DP):: alatdir ! box side in a given direction
      integer :: i, n
      vout(:) = 0.d0
      do i = 1, 3
         alatdir = alat*at(i, i)
         n = nint(vin(i)/alatdir)
         ! if (vin(i)>=0) then
         !     n=int(vin(i)/alatdir)
         ! else
         !     n=int(vin(i)/alatdir)-1
         ! end if
         vout(i) = vin(i) - dble(n)*alatdir
      end do
   end subroutine pbc_ortho_old

   subroutine pbc(vin, vout)
      use kinds, only: DP
      use cell_base, only: alat
! apply the minimum image distance in a cubic cell
      implicit none
      real(DP), intent(in) :: vin(3)
      real(DP), intent(out) :: vout(3)
      !local
      integer :: i, n
      vout(:) = 0.d0
      do i = 1, 3
         if (vin(i) >= 0) then
            n = int(vin(i)/alat)
         else
            n = int(vin(i)/alat) - 1
         end if
         vout(i) = vin(i) - dble(n)*alat
      end do
   end subroutine pbc

   subroutine check_positions(ion_pos)
      use ions_base, only: tau, nat
      use cell_base, only: alat, at
      real(DP), intent(in) ::ion_pos(3, nat)
      integer ::coord, iatom
      !local variables
      real(DP):: tau_at(3, nat) !

! tau_at(i,iatom) = at(1,i)*tau(1,iatom) +at(2,i)*tau(2,iatom) + at(3,i)*tau(3,iatom)
      tau_at = matmul(transpose(at), tau)
      do iatom = 1, nat
         do coord = 1, 3
            if (abs(tau_at(coord, iatom)*alat - ion_pos(coord, iatom)) > 1.E-4) then
               call errore('check_positions', 'positions from MD and from PW not matching', 1)
            end if
         end do
      end do
   end subroutine check_positions

!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!

   subroutine I_uno_value(y, x, a, b, flag)
      use cell_base, only: tpiba
      use gvect, only: ngm, gstart, g, gl, gg
      use mp, only: mp_sum
      use wvfct, only: npw
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::flag
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3)
      integer, intent(in)   ::a, b
      integer ::igm

!  integer ::ng_max

      real(DP) :: comp_iso

      call start_clock('rec')
!  ng_max=10

      y = 0.d0
! do igm=gstart,ngm
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0) exit
         y = y + 2.d0*(I_uno_g(igm, a, b)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
!     if (gl(igm)>ng_max) exit
      end do
      if (gstart == 2) y = y + I_uno_g(1, a, b)
      call mp_sum(y, intra_pool_comm)
      call stop_clock('rec')
      if (flag == 1) then
         call add_local_uno(y, x, a, b)
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

   subroutine add_local_uno(value, pos, a, b)
      use cell_base, only: at, alat
      use mp_world, only: nproc, mpime
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::a, b
      real(DP), intent(in) ::pos(3)
      real(DP), intent(inout) ::value
      real(DP) ::u(3), u_mod, n(3)
      integer  ::n_x, n_y, n_z
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1
!!!!!!!!!!!!!!!!!!!!!
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

   subroutine I_due_value(y, x, flag)
      use cell_base, only: tpiba
      use gvect, only: ngm, gstart, g, gl, gg
      use mp, only: mp_sum
      use wvfct, only: npw
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::flag
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3)
      integer ::igm!,ng_max
      real(DP) ::scalar
!
      call start_clock('rec')
!  ng_max=10
      y = 0.d0
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0) exit
         y = y + 2.d0*(I_due_g(igm)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
!     if        (gl(igm)>ng_max) exit
      end do
      if (gstart == 2) y = y + I_due_g(1)
      call mp_sum(y, intra_pool_comm)
      if (flag == 1) then
         call add_local_due(y, x)
      end if
      call stop_clock('rec')
   end subroutine i_due_value

   subroutine add_local_due(value, pos)
      use cell_base, only: at, alat
      use mp, only: mp_sum
      use mp_world, ONLY: nproc, mpime
      use mp_pools, ONLY: intra_pool_comm
      implicit none
      real(DP), intent(in) ::pos(3)
      real(DP), intent(inout) ::value
      real(DP) ::modul, n(3), erf_value
      integer  ::n_x, n_y, n_z
      real(kind=DP), external :: qe_erfc
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1

!
!!!!!!!!!!!
      !This routines computes sum_{L}erfc(sqrt(eta)abs(d-L))/abs(d-L)
      ! the sum is in real space over a small number of shells
      call start_clock('real')
      l_blk = (2*n_max + 1)/nproc
      if (l_blk*nproc < (2*n_max + 1)) l_blk = l_blk + 1
      nbegin = mpime*l_blk - n_max
      nend = nbegin + l_blk - 1
      if (nend > n_max) nend = n_max
      value1 = 0.d0
      do n_x = nbegin, nend
         do n_y = -n_max, n_max      !secondo ciclo su n
            do n_z = -n_max, n_max   !terzo ciclo   su n
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



subroutine init_zero()
!called once to init stuff that does not depend on the atomic positions
   use gvect, only: ngm, gg, g, gstart
   use hartree_mod
   use ions_base, only: nsp
   implicit none
   integer :: isp, iun, a, b
   character(256) :: pref_box
   integer, external :: find_free_unit
   logical :: exst
   call start_clock('init_zero')
   call init_us_1all()
   call init_reciprocal_parts_tab()
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


subroutine routine_zero(nbnd, npwx, npw, dffts, nsp, zv, nat, ityp, amass, tau, &
                        vel, tpiba, tpiba2, at, alat, omega, psic, evc, ngm, gg, g, gstart, &
                        nkb, vkb, deeq, upf, nh, xk, igk_k )
   use kinds, only: DP
   !use wvfct, only: nbnd, npwx, npw
   use fft_base, only: fft_type_descriptor !, dffts
   use mp, only: mp_sum
   use mp_pools, only: intra_pool_comm
   use io_global, only: stdout, ionode
   !
   !use ions_base, only: nsp, zv, nat, ityp, amass, tau
   !use dynamics_module, only: vel
   !use cell_base, only: tpiba, tpiba2, at, alat, omega
   !use wavefunctions, only: psic, evc
   !use gvect, only: ngm, gg, g, gstart
   use constants, only: e2, AMU_RY
   !use uspp, only: nkb
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
                           g(:,:), gg(:)
   COMPLEX(DP), intent(inout) :: psic(:), evc(:,:)
   COMPLEX(DP), intent(in):: vkb(:,:)
   real(dp), intent(in) :: deeq(:,:,:,:), xk(:,:)
   type(pseudo_upf), intent(in) :: upf(:)

!three dimensional auxiliary real vectors
   real(dp), parameter  :: amconv = AMU_RY
   real(DP) ::n(3), u(3), u_pbc(3), u_x(3), u_y(3), u_z(3), value, x(1:3), ics(3)
   real(DP), external :: qe_erfc
   real(DP)  :: dtauij(3)
   real(DP) ::modul, erf_value, charge_atom
   real(DP) ::fac_uno, fac_due, fac_tre, fac_quattro
   real(DP), external :: qe_erf
   integer, external  :: find_free_unit

!counters and indexes
   integer      ::a, b !cartesian components
   integer      ::igm, isp, jgm !reciprocal lattice,types, ions
   integer      ::nr1a, nr2a, nr3a
   integer      ::im, iqq !pseudo mesh, real_mesh_one_coordinate
   integer      ::ii, iv !processors,bands
   integer      ::ix, iy, iz !real_mesh_three_coordinates
   integer      ::iun !unit
   integer      ::iatom, jatom !atoms
   integer      ::n_x, n_y, n_z, n_maxl
   integer      ::istep

!variables for calling functions
   logical :: exst, l_test, l_scambio, l_scambio_alt

!auxiliary variables
   integer ::err, ir, ieta, enne
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

!!
!!------ IONIC CURRENT ----------------
!!
   i_current = 0.d0
   i_current_a = 0.d0
   i_current_b = 0.d0
   i_current_c = 0.d0
   i_current_d = 0.d0
   i_current_e = 0.d0

   call start_clock('calcolo_i')
   do iatom = 1, nat
      i_current_a(:) = i_current_a(:) + vel(:, iatom)*(1./2.*amconv*amass(ityp(iatom))*(vel(1, iatom)**2 +&
                       vel(2, iatom)**2 + vel(3, iatom)**2))
      i_current_b(:) = i_current_b(:) + 2./3.*e2*zv(ityp(iatom))**2*vel(:, iatom)*I_primo
   end do

   do iatom = 1, nat
      do jatom = 1, nat
         if (iatom > jatom) then
            u(1:3) = (tau(:, iatom) - tau(:, jatom))*alat
            call pbc_ortho(u(1:3), u_pbc(1:3))
            call I_due_value(value, u_pbc, 1)
            !i_current_c(:) = i_current_c(:) + 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
            i_current_c(:) = i_current_c(:) + 1.d0*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
                             (vel(:, iatom) + vel(:, jatom))*value

            do a = 1, 3
               do b = 1, 3
                  if (a > b) then
                     call I_uno_value(value, u_pbc, a, b, 1)
                     i_current_e(a) = i_current_e(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*value
                     i_current_e(b) = i_current_e(b) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(a, jatom) + vel(a, iatom))*value
                  end if
                  if (a == b) then
                     call I_uno_value(value, u_pbc, a, b, 1)
                     i_current_d(a) = i_current_d(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*value
                  end if
               end do
            end do
         end if
      end do
   end do
   i_current_a = i_current_a * alat**3
   i_current_b = i_current_b * alat 
   i_current_c = i_current_c * alat 
   i_current_d = i_current_d * alat 
   i_current_e = i_current_e * alat 
   i_current   = i_current_a + i_current_c + i_current_d + &
                 i_current_e
   if (add_i_current_b) &
      i_current = i_current + i_current_b
   

   call stop_clock('calcolo_i')
   call print_clock('calcolo_i')
   if (ionode) print *, 'IONIC CURRENT CALCULATED'

300 call stop_clock('routine_zero')
   call print_clock('routine_zero')

end subroutine routine_zero


subroutine add_nc_curr(current, nkb, vkb, deeq, upf, nh, vel, nbnd, npw, npwx, evc, &
                       g, tpiba, nat, ityp, nsp, xk, igk_k)
!
   use kinds, only: DP
   use becmod
   !use uspp, ONLY: nkb, vkb, deeq
   !USE uspp_param, ONLY: upf, nh
   !use hartree_mod, only: evc
   !use dynamics_module, only: vel
   !use wvfct, ONLY: nbnd, npw, npwx
   !use wavefunctions, only: evc
   !use gvect, ONLY: g
   !use cell_base, ONLY: tpiba
   !use ions_base, ONLY: nat, ityp, nsp
   !use klist, only: xk, igk_k
   !use io_global, only: ionode
   USE uspp_param, ONLY: pseudo_upf
!
   implicit none
   real(DP), intent(inout) :: current(3)
   integer, intent(in) :: nkb, nh(:), nbnd, npw, npwx, nat, ityp(:), &
                          nsp, igk_k(:,:)
   complex(DP), intent(in) :: vkb(:,:), evc(:,:)
   real(dp), intent(in) :: deeq(:,:,:,:), vel(:,:), g(:,:), tpiba, xk(:,:)
   type(pseudo_upf), intent(in) :: upf(:)


!
   integer ::iun, a, b
   logical ::l_test, exst
   integer :: ipol, jpol, ipwd, ijk, ibnd
   integer ::ijkb, ikb, ih, na, nt, ipw
   real(DP) ::J_nl(3), J_1(3), J_2(3)
   complex(DP), allocatable ::write (:, :)
   complex(DP), allocatable ::vkb1(:, :)
   integer, external       :: find_free_unit
!
   allocate (vkb1(npwx, nkb))
   allocate (write (npwx, nbnd))
!
!non fare nulla se non vi è potenziale non locale.
   if (nkb <= 0) return
!check se vi sono pseudo ultrasoffici o PAW ed in tal caso blocca il calcolo.
   do nt = 1, nsp
      if ((upf(nt)%typ .eq. "US") .or. (upf(nt)%typ .eq. "PAW")) then
         CALL errore('add_nc_curr', 'US and PAW not implemented', 1)
      end if
   end do
!
!memory initialization
   CALL allocate_bec_type(nkb, nbnd, becp)
   do ipol = 1, 3
      CALL allocate_bec_type(nkb, nbnd, becpr(ipol))
      CALL allocate_bec_type(nkb, nbnd, becpd(ipol))
      do jpol = 1, 3
         CALL allocate_bec_type(nkb, nbnd, becprd(ipol, jpol))
      end do
   end do
!
!allocazione delle vkb (vkb già allocato?)
   allocate (xvkb(npwx, nkb, 3))
   allocate (dvkb(npwx, nkb, 3))
   allocate (xdvkb(npwx, nkb, 3, 3))

!
!inizializzazione di tab(serve?),tabr ed indici
   call init_us_1a()
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
   use ions_base, only: nsp, nat
   use gvect, only: ngm
   use atom, only: rgrid
   USE us, ONLY: nqxq, spline_ps
   USE uspp_param, ONLY: upf, nbetam
   use fft_base, only: dffts
   use wvfct, ONLY: npwx, nbnd
!
   implicit none
!
   integer ::isp

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

END MODULE zero_mod

