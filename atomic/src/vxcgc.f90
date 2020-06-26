!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine vxc_t(lsd,rho,rhoc,exc,vxc)
  !---------------------------------------------------------------
  !
  !  this function returns the XC potential and energy in LDA or
  !  LSDA approximation
  !
  use kinds, only : DP
  use xc_lda_lsda, only: xc
  implicit none
  integer, intent(in)  :: lsd ! 1 in the LSDA case, 0 otherwise
  real(DP), intent(in) :: rho(2), rhoc ! the system density
  real(DP), intent(out):: exc(1), vxc(2)
  !
  integer,  parameter :: length=1
  real(DP), dimension(length) :: ex, ec , arho
  REAL(DP), dimension(length,2) :: rhoaux, vx, vc
  !
  real(DP), parameter :: e2=2.0_dp, eps=1.e-30_dp
  !
  vxc(1) = 0.0_dp
  exc    = 0.0_dp
  !
  if (lsd == 0) then
     !
     !     LDA case
     !
     rhoaux(1,1) = abs(rho(1) + rhoc)
     if (rhoaux(1,1) > eps) then
        !
        CALL xc( length, 1, 1, rhoaux, ex, ec, vx(:,1), vc(:,1) )
        !
        vxc(1) = e2 * ( vx(1,1) + vc(1,1) )
        exc    = e2 * ( ex(1)   + ec(1)   )
     endif
  else
     !
     !     LSDA case
     !
     vxc(2)=0.0_dp
     !
     rhoaux(1,1) = rho(1) + rho(2) + rhoc
     rhoaux(1,2) = rho(1) - rho(2)
     !
     CALL xc( length, 2, 2, rhoaux, ex, ec, vx, vc )
     !
     vxc(1) = e2 * ( vx(1,1) + vc(1,1) )
     vxc(2) = e2 * ( vx(1,2) + vc(1,2) )
     exc    = e2 * ( ex(1)   + ec(1)   )
     !
  endif
  !
  return
  !
end subroutine vxc_t
!
!
!---------------------------------------------------------------
subroutine vxcgc( ndm, mesh, nspin, r, r2, rho, rhoc, vgc, egc, &
                  tau, vtau, iflag )
  !---------------------------------------------------------------
  !
  !
  !     This routine computes the exchange and correlation potential and
  !     energy to be added to the local density, to have the first
  !     gradient correction.
  !     In input the density is rho(r) (multiplied by 4*pi*r2).
  !
  !     The units of the potential are Ry.
  !
  use kinds,     only : DP
  use constants, only : fpi, e2
  use funct,     only : dft_is_meta
  use xc_gga,    only : xc_gcx, change_threshold_gga
  use metagga,   only : tpsscxc
  implicit none
  integer,  intent(in) :: ndm,mesh,nspin,iflag
  real(DP), intent(in) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm)
  real(DP), intent(out):: vgc(ndm,2), egc(ndm)
  real(DP), intent(in) :: tau(ndm,2)
  real(DP), intent(out):: vtau(mesh)

  integer :: i, is, ierr
  real(DP) :: sx, sc, v2c, v1x, v2x, v1c
  !
  REAL(DP) :: grho_v(3,mesh,nspin)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sx_v, sc_v, v2c_ud
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v1x_v, v2x_v, v1c_v, v2c_v
  !
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  real(DP) :: v3x, v3c, de_cc, dv1_cc,dv2_cc
  real(DP) :: segno, arho
  real(DP) :: rh(1), zeta(1), grh2(1), grho2(2)
  real(DP),parameter :: eps=1.e-12_dp, small = 1.E-10_DP

  real(DP), allocatable :: grho(:,:), h(:,:), dh(:), rhoaux(:,:)
  !
  !      First compute the charge and the charge gradient, assumed  
  !      to have spherical symmetry. The gradient is the derivative of
  !      the charge with respect to the modulus of r. 
  !
  allocate(rhoaux(mesh,nspin),stat=ierr)
  allocate(grho(mesh,nspin),stat=ierr)
  allocate(h(mesh,2),stat=ierr)
  allocate(dh(mesh),stat=ierr)
  
  egc=0.0_dp
  vgc=0.0_dp

  do is=1,nspin
     do i=1, mesh
        rhoaux(i,is)=(rho(i,is)+rhoc(i)/nspin)/fpi/r2(i)
     enddo
     call radial_gradient(rhoaux(1,is),grho(1,is),r,mesh,iflag)
  enddo
  !
  do is=1,nspin
     do i=1, mesh
        grho_v(:,i,is) = grho(i,is)/SQRT(3.d0)
     enddo
  enddo
  !
  allocate( sx_v(mesh) , sc_v(mesh)  )
  allocate( v1x_v(mesh,nspin), v2x_v(mesh,nspin) )
  allocate( v1c_v(mesh,nspin), v2c_v(mesh,nspin) )
  IF (nspin==2) allocate( v2c_ud(mesh) )
  !
  if (nspin.eq.1) then
     !
     IF ( dft_is_meta ()  ) THEN
        !
        !  meta-GGA case
        !
        ! for core correction - not implemented
        de_cc = 0.0_dp
        dv1_cc= 0.0_dp
        dv2_cc= 0.0_dp
        !
        vtau(:) = 0.0_dp
        ! 
       do i=1,mesh
           arho=abs(rhoaux(i,1)) 
           segno=sign(1.0_dp,rhoaux(i,1))
           if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
!
! currently there is a single meta-GGA implemented (tpss)
! that calculates all needed terms (LDA, GGA, metaGGA)
!
             call tpsscxc ( arho, grho(i,1)**2, tau(i,1)+tau(i,2), &
                   sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )
              !
              egc(i)=sx+sc+de_cc
              vgc(i,1)= v1x+v1c + dv1_cc
              h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
              vtau(i) = v3x+v3c
          else
              vgc(i,1)=0.0_dp
              egc(i)=0.0_dp
              h(i,1)=0.0_dp
              vtau(i)=0.0_dp
           endif
        end do

     ELSE
        !
        !     GGA case
        !
        CALL change_threshold_gga( small, eps**2 )
        !
        CALL xc_gcx( mesh, nspin, rhoaux, grho_v, sx_v, sc_v, v1x_v, v2x_v, v1c_v, v2c_v )
        !
        egc(1:mesh) = sx_v + sc_v
        vgc(1:mesh,1) = v1x_v(1:mesh,1) + v1c_v(1:mesh,1) 
        h(1:mesh,1) = ( v2x_v(1:mesh,1) + v2c_v(1:mesh,1) ) * grho(1:mesh,1)*r2(1:mesh)
        !
     END IF
     !
  ELSE
     !
     !   this is the \sigma-GGA case
     !
     CALL change_threshold_gga( small, small )
     !
     CALL xc_gcx( mesh, 2, rhoaux, grho_v, sx_v, sc_v, v1x_v, v2x_v, v1c_v, v2c_v, v2c_ud )
     !
     do i = 1, mesh
        egc(i)=sx+sc
        vgc(i,1) = v1x_v(i,1)+v1c_v(i,1)
        vgc(i,2) = v1x_v(i,2)+v1c_v(i,2)
        h(i,1) =((v2x_v(i,1)+v2c_v(i,1))*grho(i,1)+v2c_v(i,1)*grho(i,2))*r2(i)
        h(i,2) =((v2x_v(i,2)+v2c_v(i,1))*grho(i,2)+v2c_v(i,1)*grho(i,1))*r2(i)
        !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
        !                          rho(i,1)*2.0_dp, grho(i,1)**2*4.0_dp, &
        !                          vgc(i,1),  h(i,2)
     enddo
!
  endif
  !
  deallocate( sx_v  , sc_v  )
  deallocate( v1x_v , v2x_v )
  deallocate( v1c_v , v2c_v )
  IF (nspin==2) deallocate( v2c_ud )
  !     
  !     We need the gradient of h to calculate the last part of the exchange
  !     and correlation potential.
  !     
  do is=1,nspin
     call radial_gradient(h(1,is),dh,r,mesh,iflag)
     !
     !     Finally we compute the total exchange and correlation energy and
     !     potential. We put the original values on the charge and multiply
     !     by e^2 = two to have as output Ry units.

     do i=1, mesh
        vgc(i,is)=vgc(i,is)-dh(i)/r2(i)
        vgc(i,is)=e2*vgc(i,is)
        if (is.eq.1) egc(i)=e2*egc(i)
        !            if (is.eq.1.and.i.lt.4) write(6,'(3f20.12)') &
        !                                      vgc(i,1)
     enddo
  enddo
  IF ( dft_is_meta() ) vtau(:) = e2*vtau(:)

  deallocate(dh)
  deallocate(h)
  deallocate(grho)
  deallocate(rhoaux)

  return
end subroutine vxcgc
