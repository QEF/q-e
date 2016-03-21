subroutine make_v_state(numb_v,v)
  USE gvect,                 ONLY : gstart
  USE lsda_mod,              ONLY : nspin
  use wavefunctions_module,  ONLY : evc
  use io_files,  ONLY : prefix, iunwfc
  USE wvfct,    ONLY : nbnd, npwx,npw
  implicit none

  type(v_state) :: v
  integer :: numb_v
  
  integer :: is,ivmax

  v%nspin=nspin
  v%numb_v(:)=numb_v(:)
  v%npw=npw
  v%gstart=gstart

  allocate( evc( npwx, nbnd ) )
  
  if (nspin==1) then
     ivmax= v%numb_v(1)
  else 
     ivmax=max(v%numb_v(1),v%numb_v(2))
  endif

  allocate( v%wfn(v%npw,ivmax,v%nspin)
  allocate( v%esp(ivmax,v%nspin)

  do is=1,nspin
     call davcio(evc,2*nwordwfc,iunwfc,is,-1)
     do iv=1,v%numb_v(is)
        v%wfn(1:v%npw,1:v%numb_v(is),is)=evc(1:v%npw,1:v%numb_v(is))
     enddo  
        v%esp(1:v%numb_v(is),is)=et(1:v%numb_v(is),is)
  enddo

  deallocate(evc)

  


  return
end subroutine

