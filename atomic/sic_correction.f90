!
!---------------------------------------------------------------
subroutine sic_correction(n,vhn1,vhn2,egc) 
  !---------------------------------------------------------------
  !   set up the orbital-dependent selfconsistent potential generated
  !   by the n-th wavefunction - for self-interaction correction
  !
  use constants, only: e2, fpi
  use ld1inc
  use funct
  implicit none
  integer :: n
  real(kind=dp):: vhn1(ndm),vhn2(ndm), egc(ndm)
  !
  integer :: i, is
  real(kind=dp):: rh(2), rhc, exc_t, vxcp(2)
  real(kind=dp):: vgc(ndm,2),  egc0(ndm), rhotot(ndm,2)
  logical :: gga

  vhn1=0.0_dp
  vhn2=0.0_dp
  gga=igcx.ne.0.or.igcc.ne.0
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the charge of orbital n
  !
  rhotot=0.0_dp
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
  !    add exchange and correlation potential: LDA or LSDA terms
  !
  rhc=0.0_dp
  rh=0.0_dp
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
  !   add gradient-correction terms to exchange-correlation potential
  !
  egc0=egc
  call vxcgc(ndm,mesh,nspin,r,r2,rhotot,rhoc,vgc,egc)
  do i=1,mesh
     vhn2(i)=vhn2(i)+vgc(i,1)
     egc(i)=egc(i)*r2(i)*fpi+egc0(i)
  enddo
  return
end subroutine sic_correction
