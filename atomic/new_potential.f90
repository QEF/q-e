!
!---------------------------------------------------------------
subroutine new_potential &
     (ndm,mesh,r,r2,sqr,dx,zed,vxt,lsd, &
     nlcc,latt,enne,rhoc,rho,vh,vnew)
  !---------------------------------------------------------------
  !   set up the selfconsistent atomic potential
  !
  use kinds, only : DP
  use funct
  implicit none
  logical :: nlcc, gga
  integer :: ndm,mesh,lsd,latt,i,is,nspin, ierr
  real(kind=dp):: rho(ndm,2),r(ndm),r2(ndm),vxcp(2), &
       sqr(ndm),vnew(ndm,2),vxt(ndm),vh(ndm), rhoc(ndm)
  real(kind=dp):: zed,enne,dx,rh(2),rhc
  real(kind=dp),allocatable:: vgc(:,:), egc(:), rhotot(:)
  real(kind=dp),parameter ::fourpi=4.0_DP*3.141592653589793_DP, e2=2.0_DP

  gga=igcx.ne.0.or.igcc.ne.0
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the total charge
  !
  allocate(rhotot(ndm),stat=ierr)
  do i=1,ndm
     rhotot(i)=rho(i,1)
  enddo
  if (lsd.eq.1) then
     do i=1,ndm
        rhotot(i)=rhotot(i)+rho(i,2)
     enddo
  endif
  call hartree(0,2,mesh,r,r2,sqr,dx,rhotot,vh)
  !
  !    add exchange and correlation potential: LDA or LSDA only
  !
  rhc=0.0_DP
  do i=1,mesh
     vh(i) = e2*vh(i)
     do is=1,nspin
        rh(is) = rho(i,is)/r2(i)/fourpi
     enddo
     if (nlcc) rhc = rhoc(i)/r2(i)/fourpi
     call vxc_t(rh,rhc,lsd,vxcp)
     do is=1,nspin
        vnew(i,is)= - zed*e2/r(i)+vxt(i)+vh(i)+vxcp(is)
     enddo
  end do

  deallocate(rhotot)

  if (.not.gga) goto 1010
  !
  !   add exchange and correlation potential: GGA only
  !
  allocate(vgc(ndm,2),stat=ierr)
  allocate(egc(ndm),stat=ierr)
  call errore('new_potential','allocating vgc and egc',ierr)

  call vxcgc(ndm,mesh,nspin,r,r2,rho,rhoc,vgc,egc)
  do is=1,nspin
     do i=1,mesh
        vnew(i,is)=vnew(i,is)+vgc(i,is)
     enddo
  enddo

  deallocate(egc)
  deallocate(vgc)

1010 if (latt.ne.0) then
     do is=1,nspin
        do i=1,mesh
           vnew(i,is)= min(vnew(i,is),-e2*(zed-enne+1.0_DP)/r(i))
        enddo
     enddo
  end if

  return
end subroutine new_potential
