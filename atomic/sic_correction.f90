!
!---------------------------------------------------------------
subroutine sic_correction(n,vhn1,vhn2,egc) 
  !---------------------------------------------------------------
  !   set up the selfconsistent atomic potential for the charge of
  !   the wavefunction n
  !
  use constants, only: e2
  use ld1inc
  use funct
  implicit none
  integer :: i, is, n
  real(kind=dp):: vhn1(ndm),vhn2(ndm), rh(2), rhc, vxcpl(2), rhotot(ndm,2),&
       fpi 
  real(kind=dp):: vgc(ndm,2), egc(ndm), egc0(ndm)
  real(kind=dp):: exc_t, vxcp(2)
  logical :: gga

  vhn1=0.d0
  vhn2=0.d0
  fpi=16.d0*atan(1.d0)
  gga=igcx.ne.0.or.igcc.ne.0
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the total charge
  !
  rhotot=0.d0
  if (rel.eq.2) then
     do i=1,mesh
        rhotot(i,1)=psi_dir(i,1,n)**2+psi_dir(i,2,n)**2
     enddo
  else
     do i=1,mesh
        rhotot(i,1)=psi(i,n)**2
     enddo
  endif
  !call hartree(0,2*(ll(n)+1),mesh,r,r2,sqr,dx,rhotot,vhn1)
  call hartree(0,2,mesh,r,r2,sqr,dx,rhotot,vhn1)
  !
  !    add exchange and correlation potential: LDA or LSDA only
  !
  rhc=0.d0
  rh=0.d0
  do i=1,mesh
     vhn1(i) = e2*vhn1(i)
     rh(1) = rhotot(i,1)/r2(i)/fpi
     if (nlcc) rhc = rhoc(i)/r2(i)/fpi
     call vxc_t(rh,rhc,lsd,vxcp)
     vhn2(i)= vhn1(i)+vxcp(1)
     egc(i)= exc_t(rh,rhc,lsd)*rhotot(i,1)
  end do

  if (.not.gga) return
  !
  !   add exchange and correlation potential: GGA only
  !
  egc0=egc
  call vxcgc(ndm,mesh,nspin,r,r2,rhotot,rhoc,vgc,egc)
  do i=1,mesh
     vhn2(i)=vhn2(i)+vgc(i,1)
     egc(i)=egc(i)*r2(i)*fpi+egc0(i)
  enddo
  return
end subroutine sic_correction
