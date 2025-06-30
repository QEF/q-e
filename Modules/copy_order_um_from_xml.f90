
SUBROUTINE copy_order_um_from_xml(dftU_obj, nsp, nat, Hubbard_lmax, order_um)
  !! routine that copies order_um values from xml file to allocatable order_um, used 
  !! for the initialization of the order_um variable in ldaU module. 
  !FIXME to be collected with other related copying from dftU object or better read from the scf file
  USE qes_types_module, only: dftU_type
  IMPLICIT NONE
  TYPE(dftU_type),INTENT(IN)  :: dftU_obj
  !! object containing the dftU infortmation read from the XML file
  INTEGER,INTENT(IN)          :: nsp, nat, Hubbard_lmax
  !! number of species defined for this calculation 
  !! number of atoms for this calculation
  !! max l value for species in this calculatio 
  INTEGER, ALLOCATABLE,INTENT(INOUT) :: order_um(:,:,:)  
  !! allocatable 3-rank array to be written dimensions will be 2*Hubbarl_lmax+1, nspin, nat  
  INTEGER :: nspin, iobj, ll, ia, ispin 
  nspin = 1 
  IF (dftU_obj%Hub_m_order(1)%spin_ispresent) nspin = 2  
  IF (ALLOCATED(order_um)) DEALLOCATE(order_um)
  ALLOCATE(order_um(2*Hubbard_lmax+1, nspin, nat)) 
  order_um = 0 
  DO iobj = 1, dftU_obj%ndim_Hub_m_order
     ia = dftU_obj%Hub_m_order(iobj)%atomidx 
     ispin = 1 
     IF (nspin == 2) ispin = dftU_obj%Hub_m_order(iobj)%spin 
     ll = dftU_obj%Hub_m_order(iobj)%size     
     order_um(1:ll,ispin,ia) = dftU_obj%Hub_m_order(iobj)%orderUm(1:ll)
  END DO 
END SUBROUTINE copy_order_um_from_xml 
