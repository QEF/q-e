!--------------------------------------------------
subroutine deallocate_part()
!----------===============-------------------------

  use phcom

  if (allocated(comp_irr)) deallocate (comp_irr)    
  if (allocated(ifat)) deallocate (ifat)    
  if (allocated(done_irr)) deallocate (done_irr)    
  if (allocated(list)) deallocate (list)    
  if (allocated(atomo)) deallocate (atomo)    


  return
end subroutine deallocate_part
