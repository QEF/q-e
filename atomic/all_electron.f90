!
!---------------------------------------------------------------
subroutine all_electron(ild)
  !---------------------------------------------------------------
  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  use  ld1inc
  use kinds, only : DP
  implicit none

  logical :: ild    ! if true compute log der
  !
  !    compute an initial estimate of the potential
  !
  call starting_potential(ndm,mesh,zval,zed,nwf,oc,nn,ll,r,enl,vxt,vpot,&
       enne,nspin)
  !
  !     solve the eigenvalue self-consistent equation
  !
  call scf
  !
  !  compute total energy
  !
  call elsd (mesh,zed,r,r2,dx,rho,zeta,vxt,vh,nlcc,  &
       nwf,enl,ll,lsd,nspin,oc,ndm,vnl,    &
       etot,ekin,encl,epseu,ehrt,ecxc,evxt)
  !
  !   add sic correction if needed
  !
  if(isic /= 0) call esic  
  !
  !   print results
  !
  call write_results 
  !
  !  compute logarithmic derivative
  !
  if (deld > 0.0_DP .and. ild) call lderiv
  return
end subroutine all_electron
