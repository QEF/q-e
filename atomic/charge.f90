!---------------------------------------------------------------
      subroutine charge(ndm,mesh,nwf,npol,oc,chi,rho,isw)
!---------------------------------------------------------------
!
!   calculate the (spherical) charge density 
!
implicit none
integer, parameter :: dp=kind(1.d0)
integer:: ndm, mesh, nwf, isw(nwf), i, n, npol
real(kind=dp):: oc(nwf),chi(ndm,npol,nwf),rho(ndm,2)
!
rho=0.d0
!
if (npol.eq.1) then
   do n=1,nwf
      do i=1,mesh
         rho(i,isw(n))=rho(i,isw(n))+oc(n)*chi(i,1,n)**2
      enddo
   enddo
else
   do n=1,nwf
      do i=1,mesh
         rho(i,isw(n))=rho(i,isw(n))+oc(n)*(chi(i,1,n)**2+chi(i,2,n)**2)
      enddo
   enddo
endif
!
return
end
