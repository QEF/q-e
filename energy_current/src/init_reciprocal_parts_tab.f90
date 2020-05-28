SUBROUTINE init_reciprocal_parts_tab()
   use kinds, only: DP
   use fft_base, only: dffts
   use ions_base, only: nsp, zv
   use cell_base, only: tpiba2, tpiba, omega, at, alat
   use io_global, only: stdout
   use atom, only: rgrid
   use constants, only: pi, fpi
   use mp_world, only: mpime
   use zero_mod
   use fft_interfaces, only: invfft
   use mp, only: mp_sum, mp_min
   use gvect, only: ngm, gg, gstart, g, igtongl, gl, ngl
   use us, only: spline_ps, dq
   use mp_pools, only: intra_pool_comm
!  use splinelib
   implicit none

   real(DP) :: px, ux, vx, wx, xg
   real(DP) :: u(3), u_x(3), u_y(3), u_z(3), n(3), modul, erf_value
   integer  :: isp, a, b, err, i, igm, ii, igl, it, n_x, n_y, n_z
   integer  :: i0, i1, i2, i3
   real(DP) :: integral, ra, rb
   real(DP) :: H_g_rad(ngl, 0:1)
   logical  ::l_plot
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
