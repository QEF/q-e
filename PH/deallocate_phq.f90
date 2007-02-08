!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------
subroutine deallocate_phq
!----------========-----------------------
!
!  deallocates the variables allocated by allocate_phq
!
  use phcom
  USE noncollin_module, ONLY : m_loc
  USE becmod, ONLY: becp, becp_nc
  USE wavefunctions_module,  ONLY: evc
  USE ramanm, ONLY: ramtns
  use el_phon

  if(allocated(ramtns)) deallocate (ramtns)  
  if (lgamma) then
     if(associated(evq)) nullify(evq)
     if(associated(igkq)) nullify(igkq)
  else
     if(associated(evq)) deallocate(evq)
     if(associated(igkq)) deallocate(igkq)
  end if
  
  if(allocated(dvpsi)) deallocate (dvpsi)    
  if(allocated(dpsi)) deallocate ( dpsi)    
  !
  if(allocated(vlocq)) deallocate (vlocq)
  if(allocated(dmuxc)) deallocate (dmuxc)
  !
  if(allocated(eigqts)) deallocate (eigqts)
  if(allocated(rtau)) deallocate (rtau)
  if(associated(u)) deallocate (u)
  if(associated(ubar)) deallocate (ubar)
  if(allocated(dyn)) deallocate (dyn)
  if(allocated(dyn00)) deallocate (dyn00)
  if(allocated(w2)) deallocate (w2)
  if(associated(t)) deallocate (t)
  if(associated(tmq)) deallocate (tmq)

  if(allocated(npert)) deallocate (npert)    
  if(allocated(zstareu)) deallocate (zstareu)    
  if(allocated(zstareu0)) deallocate (zstareu0)    
  if(allocated(zstarue)) deallocate (zstarue)    
  if(allocated(zstarue0)) deallocate (zstarue0)    

  if(allocated(int1)) deallocate (int1)    
  if(allocated(int2)) deallocate (int2)
  if(allocated(int3)) deallocate (int3)
  if(allocated(int4)) deallocate (int4)
  if(allocated(int5)) deallocate (int5)
  if(allocated(dpqq)) deallocate (dpqq)
  if(allocated(int1_nc)) deallocate(int1_nc)
  if(allocated(int3_nc)) deallocate(int3_nc)
  if(allocated(int4_nc)) deallocate(int4_nc)
  if(allocated(becsum_nc)) deallocate(becsum_nc)
  if(allocated(alphasum_nc)) deallocate(alphasum_nc)
  if(allocated(int2_so)) deallocate(int2_so)
  if(allocated(int5_so)) deallocate(int5_so)
!  if (allocated(dpqq_so)) deallocate(dpqq_so)

  if(allocated(alphasum)) deallocate (alphasum)
  if(allocated(this_dvkb3_is_on_file)) deallocate (this_dvkb3_is_on_file)    


  if(allocated(this_pcxpsi_is_on_file)) deallocate (this_pcxpsi_is_on_file)
  if(allocated(alphap_nc))    deallocate (alphap_nc)
  if(allocated(becp1_nc))     deallocate(becp1_nc)
  if(allocated(becp_nc))      deallocate(becp_nc)

  if(allocated(alphap))    deallocate (alphap)    
  if(allocated(becp1))     deallocate(becp1) 
  if(allocated(becp))      deallocate(becp)
  if(allocated(el_ph_mat)) deallocate (el_ph_mat)    
  if(allocated(m_loc))     deallocate(m_loc)

  if(allocated(drc)) deallocate(drc)

  if(allocated(dvxc_rr)) deallocate (dvxc_rr)    
  if(allocated(dvxc_sr)) deallocate (dvxc_sr)    
  if(allocated(dvxc_ss)) deallocate (dvxc_ss)    
  if(allocated(dvxc_s)) deallocate (dvxc_s)    
  if(allocated(grho)) deallocate (grho)  

  return
end subroutine deallocate_phq
