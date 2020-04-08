subroutine deallocate_zero
use zero_mod,   only:H_g,I_uno_g,I_due_g,status,&
tabr,tabr_d2y,charge_g,u_g,charge,ion_vel !,ion_pos
!
implicit none
!
!da completare con deallocazioni per altri status, se necessario
if (status=='compute') then
    deallocate(charge)
    !deallocate(ion_pos)
    deallocate(ion_vel)
    deallocate(charge_g)
    deallocate(u_g) 
end if
!
end subroutine deallocate_zero
