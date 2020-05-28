subroutine deallocate_zero
   use zero_mod, only: H_g, I_uno_g, I_due_g, &
                       tabr, tabr_d2y, charge_g, u_g, charge, ion_vel
!
   implicit none
!
   if (allocated(charge)) &
      deallocate (charge)
   if (allocated(ion_vel)) &
      deallocate (ion_vel)
   if (allocated(charge_g)) &
      deallocate (charge_g)
   if (allocated(u_g)) &
      deallocate (u_g)
!
end subroutine deallocate_zero
