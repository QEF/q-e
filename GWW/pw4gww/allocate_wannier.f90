! FOR GWW
! Author: P. Umari
! Modified by G. Stenuit
!
subroutine  allocate_wannier
! these subroutines allocate and deallocate  what's needed for the wannier functions

!#ifdef __GWW

  USE wannier_gw,   ONLY :  wannier_centers, wannier_radii, u_trans, w_centers, w_radii, becp_gw, num_nbndc_set,&
       & becp_gw_c, w_centers_c, w_radii_c, nbnd_normal
  USE wvfct,     ONLY :  nbnd
  USE uspp, ONLY : okvan,nkb

  implicit none

!set up nbnd_normal
  if(nbnd_normal<=0 .or. nbnd_normal > nbnd) nbnd_normal = nbnd


  allocate(wannier_centers(3,nbnd_normal))
  allocate(wannier_radii(nbnd_normal))
  allocate(u_trans(nbnd_normal,nbnd_normal))
  allocate(w_centers(3,nbnd_normal))
  allocate(w_radii(nbnd_normal))
  if(okvan) then
     allocate(becp_gw(nkb,nbnd))
     allocate(becp_gw_c(nkb,nbnd))
  endif


  if(num_nbndc_set > 0) then
     allocate(w_centers_c(3,num_nbndc_set))
     allocate(w_radii_c(num_nbndc_set))
  endif


!#endif __GWW
  return
end subroutine allocate_wannier

subroutine  deallocate_wannier
!#ifdef __GWW
  USE wannier_gw,   ONLY :  wannier_centers, wannier_radii, u_trans, w_centers, w_radii, becp_gw, &
       & w_radii_c, w_centers_c, becp_gw_c

  implicit none

  if(allocated(wannier_centers)) deallocate(wannier_centers)
  if(allocated(wannier_radii))   deallocate(wannier_radii)
  if(allocated(u_trans))         deallocate(u_trans)
  if(allocated(w_centers))       deallocate(w_centers)
  if(allocated(w_radii))         deallocate(w_radii)
  if(allocated(becp_gw))         deallocate(becp_gw)
  if(allocated(w_centers_c))       deallocate(w_centers_c)
  if(allocated(w_radii_c))         deallocate(w_radii_c)
  if(allocated(becp_gw_c))         deallocate(becp_gw_c)


!#endif
  return
end subroutine deallocate_wannier
