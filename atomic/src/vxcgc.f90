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
  use funct, only : xc, xc_spin
  implicit none
  integer, intent(in)  :: lsd
  real(DP), intent(in) :: rho(2), rhoc
  real(DP), intent(out):: exc, vxc(2)
  real(DP):: arho, zeta, vx(2), vc(2), ex, ec
  !
  real(DP), parameter :: e2=2.0_dp, eps=1.e-30_dp

  vxc(1)=0.0_dp
  exc=0.0_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     arho=abs(rho(1)+rhoc)
     if (arho.gt.eps) then      
        call xc(arho,ex,ec,vx(1),vc(1))
        vxc(1)=e2*(vx(1)+vc(1))
        exc   =e2*(ex+ec)
     endif
  else
     !
     !     LSDA case
     !
     vxc(2)=0.0_dp
     arho = abs(rho(1)+rho(2)+rhoc)
     if (arho.gt.eps) then      
        zeta = (rho(1)-rho(2)) / arho
        ! zeta has to stay between -1 and 1, but can get a little
        ! out the bound during the first iterations.
        if (abs(zeta).gt.1.0_dp) zeta = sign(1._dp, zeta)
        call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
        vxc(1) = e2*(vx(1)+vc(1))
        vxc(2) = e2*(vx(2)+vc(2))
        exc    = e2*(ex+ec)
     endif
  endif

  return
end subroutine vxc_t
!
!---------------------------------------------------------------
subroutine vxcgc ( ndm, mesh, nspin, r, r2, rho, rhoc, vgc, egc, &
     tau, vtau, iflag)
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
  use kinds, only : DP
  use constants, only : fpi, e2
  use funct, only : gcxc, gcx_spin, gcc_spin, dft_is_meta, xc
  implicit none
  integer,  intent(in) :: ndm,mesh,nspin,iflag
  real(DP), intent(in) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm)
  real(DP), intent(out):: vgc(ndm,2), egc(ndm)
  real(DP), intent(in) :: tau(ndm,2)
  real(DP), intent(out):: vtau(mesh)

  integer :: i, is, ierr
  real(DP) :: sx,sc,v1x,v2x,v1c,v2c
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  real(DP) :: v3x, v3c, de_cc, dv1_cc,dv2_cc
  real(DP) :: segno, arho
  real(DP) :: rh, zeta, grh2, grho2(2)
  real(DP),parameter :: eps=1.e-12_dp

  real(DP), allocatable :: grho(:,:), h(:,:), dh(:), rhoaux(:,:)
  !
  !      First compute the charge and the charge gradient, assumed  
  !      to have spherical symmetry. The gradient is the derivative of
  !      the charge with respect to the modulus of r. 
  !
  allocate(rhoaux(mesh,2),stat=ierr)
  allocate(grho(mesh,2),stat=ierr)
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
        do i=1,mesh
           arho=abs(rhoaux(i,1)) 
           segno=sign(1.0_dp,rhoaux(i,1))
           if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
              call gcxc(arho,grho(i,1)**2,sx,sc,v1x,v2x,v1c,v2c)
              egc(i)=(sx+sc)*segno
              vgc(i,1)= v1x+v1c
              h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
              !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
              !                          rho(i,1), grho(i,1)**2,  &
              !                          vgc(i,1),h(i,1)
           else
              vgc(i,1)=0.0_dp
              egc(i)=0.0_dp
              h(i,1)=0.0_dp
           endif
        end do
     END IF
  else
     !
     !   this is the \sigma-GGA case
     !       
     do i=1,mesh
        !
        !  NB: the special or wrong cases where one or two charges 
        !      or gradients are zero or negative must
        !      be detected within the gcxc_spin routine
        !
        !    spin-polarised case
        !
        do is = 1, nspin
           grho2(is)=grho(i,is)**2
        enddo

        call gcx_spin (rhoaux(i, 1), rhoaux(i, 2), grho2(1), grho2(2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rhoaux(i, 1) + rhoaux(i, 2)
        if (rh.gt.eps) then
           zeta = (rhoaux (i, 1) - rhoaux (i, 2) ) / rh
           grh2 = (grho (i, 1) + grho (i, 2) ) **2 
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.0_dp
           v1cup = 0.0_dp
           v1cdw = 0.0_dp
           v2c = 0.0_dp
        endif

        egc(i)=sx+sc
        vgc(i,1)= v1xup+v1cup
        vgc(i,2)= v1xdw+v1cdw
        h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
        h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
        !                          rho(i,1)*2.0_dp, grho(i,1)**2*4.0_dp, &
        !                          vgc(i,1),  h(i,2)
     enddo
  endif
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
  IF ( dft_is_meta () ) vtau(:) = e2*vtau(:)

  deallocate(dh)
  deallocate(h)
  deallocate(grho)
  deallocate(rhoaux)

  return
end subroutine vxcgc
