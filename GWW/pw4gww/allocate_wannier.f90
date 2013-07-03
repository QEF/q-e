!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!these subroutines allocate and deallocate  what's needed for the wannier functions


subroutine  allocate_wannier


  USE wannier_gw,   ONLY :  wannier_centers, wannier_radii, u_trans, w_centers, w_radii, becp_gw,   becp_gw_c, vg_q
  USE wvfct,     ONLY :  nbnd, npw,npwx
  USE uspp, ONLY : okvan,nkb
  USE lsda_mod,             ONLY : nspin


  implicit none


  allocate(wannier_centers(3,nbnd,nspin))
  allocate(wannier_radii(nbnd,nspin))
  allocate(u_trans(nbnd,nbnd,nspin))
  allocate(w_centers(3,nbnd,nspin))
  allocate(w_radii(nbnd,nspin))
  if(okvan) then
     allocate(becp_gw(nkb,nbnd,nspin))
     allocate(becp_gw_c(nkb,nbnd,nspin))
  endif
  allocate(vg_q(npwx))


  return
end subroutine allocate_wannier

subroutine  deallocate_wannier

  USE wannier_gw,   ONLY :  wannier_centers, wannier_radii, u_trans, w_centers, w_radii, becp_gw, &
       & becp_gw_c, vg_q
                                                                                                                             
  implicit none
                                                                                                                             
  if(allocated(wannier_centers)) deallocate(wannier_centers)
  if(allocated(wannier_radii))   deallocate(wannier_radii)
  if(allocated(u_trans))         deallocate(u_trans)
  if(allocated(w_centers))       deallocate(w_centers)
  if(allocated(w_radii))         deallocate(w_radii)
  if(allocated(becp_gw))         deallocate(becp_gw)
  if(allocated(becp_gw_c))       deallocate(becp_gw_c)
  if(allocated(vg_q))            deallocate(vg_q)

  
  return
                                                                                                                             
end subroutine deallocate_wannier

