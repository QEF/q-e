subroutine allocate_zero
   use zero_mod, only: H_g, I_uno_g, I_due_g, &
                       tabr, tabr_d2y, charge_g, u_g, charge, ion_vel, &
                       tablocal_hg, tablocal_d2y_hg !,ion_pos
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
   allocate (charge(dffts%nnr))
   if (.not. allocated(ion_vel)) &
      allocate (ion_vel(3, nat))
!
end subroutine allocate_zero
