subroutine init_zero()
!called once to init stuff that does not depend on the atomic positions
   use io_files, only: nwordwfc, diropn, iunwfc, prefix, tmp_dir
   use gvect, only: ngm, gg, g, gstart
   use zero_mod
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


subroutine routine_zero()
   use kinds, only: DP
   use wvfct, only: nbnd, npwx, npw
   use io_files, only: nwordwfc, diropn, iunwfc, prefix, tmp_dir
   use fft_base, only: dffts
   use mp, only: mp_sum, mp_bcast, mp_get
   use mp_pools, only: intra_pool_comm
   use io_global, only: stdout, ionode_id, ionode
   use ions_base, only: nsp, zv, nat, ityp, amass, tau
   use dynamics_module, only: vel
   use cell_base, only: tpiba, tpiba2
   use uspp_param, only: upf
   use atom, only: rgrid
   use mp_world, only: mpime
   use cell_base, only: at, alat, omega
   use wavefunctions, only: psic, evc
   use fft_interfaces, only: invfft, fwfft
   use gvect, only: ngm, gg, g, gstart
   use constants, only: e2, AMU_RY
   use uspp, only: nkb
   use splines
   use zero_mod
   use hartree_mod
   use compute_charge_mod, only : compute_charge

   implicit none

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

   write (stdout, *) 'ROUTINE_ZERO BEGINNING'
   call start_clock('routine_zero')
   if (ionode) print *, 'eta', eta

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
      call add_nc_curr(z_current)
   end if
   z_current = z_current*alat
   call stop_clock('zero_current')
   call print_clock('zero_current')
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
