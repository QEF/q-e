subroutine deallocate_zero
use zero_mod,   only:H_g,I_uno_g,I_due_g,status,&
tabr,tabr_d2y,charge_g,evc_uno,u_g,charge,ion_pos,ion_vel
!
implicit none
!
!da completare con deallocazioni per altri status, se necessario
if (status=='compute') then
    deallocate(evc_uno)
    deallocate(charge)
    deallocate(ion_pos)
    deallocate(ion_vel)
    deallocate(charge_g)
    deallocate(u_g) 
end if
!
end subroutine deallocate_zero
