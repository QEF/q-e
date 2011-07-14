!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Original version by Minoru Otani (AIST), Yoshio Miura (Tohoku U.),
! Nicephore Bonet (MIT), Nicola Marzari (MIT), Brandon Wood (LLNL), 
! Tadashi Ogitsu (LLNL)
!
! Contains subroutines for implementation of the ESM (Effective Screening
! Medium Method) developed by M. Otani and O. Sugino (see PRB 73, 115407 
! [2006]).
!
! ESM enables description of a surface slab sandwiched between two 
! semi-infinite media, making it possible to deal with polarized surfaces 
! without using dipole corrections. It is useful for simulating interfaces 
! with vacuum, one or more electrodes, or an electrolyte.
!
! Modified subroutines for calculating the Hartree potential, the local 
! potential, and the Ewald sum are contained here, along with subroutines for
! calculating force contributions based on the modified local potential and 
! Ewald term.
!
!----------------------------------------------------------------------------
MODULE esm
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the variables and subroutines needed for the 
  ! ... EFFECTIVE SCREENING MEDIUM (ESM) METHOD 
  !
  USE kinds, ONLY :  DP
  USE constants, ONLY : pi, eps8, tpi, e2
  SAVE
  !
  LOGICAL :: do_comp_esm=.FALSE.
  INTEGER :: esm_nfit
  REAL(KIND=DP) :: esm_efield, esm_w                 
  CHARACTER (LEN=3) :: esm_bc           
  !
  PUBLIC :: esm_hartree, esm_local, esm_ewald, esm_force_lc, esm_force_ew, &
            esm_printpot, esm_summary

CONTAINS
!
!-----------------------------------------------------------------------
!--------------ESM HARTREE SUBROUTINE-----------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_hartree (rhog, ehart, aux)

  USE gvect,     ONLY : g, nl, nlm, ngm
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE fft_scalar,       ONLY : cft_1z
  USE mp_global,        ONLY : intra_pool_comm, me_pool
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP)  :: rhog(ngm,nspin)   !  n(G)      
  REAL(DP)     :: ehart             !  Hartree energy
  COMPLEX(DP)  :: aux(dfftp%nnr)    !  v_h(G)   
  !
  !    here the local variables
  !
  real(DP)                 :: tt, t(2), zz, gz, z0, gp, gp2, z1, kn, cc, ss, z, L, &
                              z_l, z_r
  integer                  :: ipol, ijk, i, j, k, k1, k2, k3, iz, ng, n1, n2, n3, &
                              nz_r, nz_l
  complex(DP),allocatable  :: work(:,:), rhog3(:,:,:), vg2(:), vg2_in(:), vg3(:,:,:)
  complex(DP)              :: xc, ci, tmp, tmp1, tmp2, tmp3, tmp4, f1, f2, f3, f4, &
                              a0, a1, a2, a3, c_r, c_l, s_r, s_l


  allocate(vg2(dfftp%nr3x),vg2_in(dfftp%nr3x),rhog3(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x))
!
! Map to FFT mesh (nrxx)
  allocate(work(dfftp%nnr,2))
  work(nl(1:ngm),1) = rhog(1:ngm,1)
  if (nspin == 2) work(nl(1:ngm),2) = rhog(1:ngm,2)
   IF ( gamma_only ) THEN
    work(nlm(1:ngm),1) = CONJG(work(nl(1:ngm),1))
    if (nspin == 2) work(nlm(1:ngm),2) = CONJG(work(nl(1:ngm),2))
  ENDIF
! Map to FULL FFT mesh (dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  rhog3(:,:,:)=(0.d0,0.d0)
  do ng=1,ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      IF (n1<1) n1 = n1 + dfftp%nr1
      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      IF (n2<1) n2 = n2 + dfftp%nr2
      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1 
      IF (n3<1) n3 = n3 + dfftp%nr3    
#ifdef __PARA
     ijk=n3+(dfftp%isind(n1+(n2-1)*dfftp%nr1x)-1)*dfftp%nr3x
#else
     ijk=n1+(n2-1)*dfftp%nr1x+(n3-1)*dfftp%nr1x*dfftp%nr2x
#endif
     if( nspin == 2 ) then
        rhog3(n1,n2,n3)=work(ijk,1)+work(ijk,2)
     else
        rhog3(n1,n2,n3)=work(ijk,1)
     endif
  enddo
#ifdef __PARA
  call mp_sum( rhog3, intra_pool_comm )
#endif

  deallocate(work)
! End mapping
!
  allocate(vg3(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x))

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************

  do i = 1, dfftp%nr1x
     k1=i-1
     if (i.gt.dfftp%nr1x/2) k1=i-dfftp%nr1x-1
     do j = 1, dfftp%nr2x
        k2=j-1
        if (j.gt.dfftp%nr2x/2) k2=j-dfftp%nr2x-1
        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = k1 * bg (ipol, 1) + k2 * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)
        if (gp.ge.eps8) then
           tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0)
           vg2(:)=(0.d0,0.d0)
           do iz=1, dfftp%nr3x
              if(iz<=dfftp%nr3x/2) kn=dble(iz-1)     * 2.d0*pi/L
              if(iz> dfftp%nr3x/2) kn=dble(iz-1-dfftp%nr3x) * 2.d0*pi/L
              cc=cos(kn*z0)
              ss=sin(kn*z0)
              vg2(iz)=4.d0*pi*rhog3(i,j,iz)/(gp**2+kn**2)
              if (esm_bc.eq.'bc1') then
                 tmp1=tmp1+rhog3(i,j,iz)*(cc+ci*ss)/(gp-ci*kn)
                 tmp2=tmp2+rhog3(i,j,iz)*(cc-ci*ss)/(gp+ci*kn)
              else if (esm_bc.eq.'bc2') then
                 tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
                 tmp1=tmp1+rhog3(i,j,iz)*(cc+ci*ss)/(gp**2+kn**2)*tmp
                 tmp=((gp-ci*kn)*exp(gp*(z1-z0))+(gp+ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
                 tmp2=tmp2+rhog3(i,j,iz)*(cc-ci*ss)/(gp**2+kn**2)*tmp
              else if (esm_bc.eq.'bc3') then
                 tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
                 tmp1=tmp1+rhog3(i,j,iz)*(cc+ci*ss)/(gp**2+kn**2)*tmp
                 tmp=(gp-ci*kn)/gp
                 tmp2=tmp2+rhog3(i,j,iz)*(cc-ci*ss)/(gp**2+kn**2)*tmp
              endif
           enddo

           vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)  ! Since cft_1z is not in-place
           call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,1,vg2)
  
           do iz=1,dfftp%nr3x
              k3=iz-1
              if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
              z=dble(k3)/dble(dfftp%nr3x)*L
              if (esm_bc.eq.'bc1') then
                 vg2(iz)=vg2(iz)-2.d0*pi/gp*(exp(gp*(z-z0))*tmp1+exp(-gp*(z+z0))*tmp2)
              else if (esm_bc.eq.'bc2') then
                 vg2(iz)=vg2(iz)-4.d0*pi*(exp(gp*(z-z1))-exp(-gp*(z+3.d0*z1)))*tmp1 &
                                 /(1.d0-exp(-4.d0*gp*z1)) &
                                +4.d0*pi*(exp(gp*(z-3.d0*z1))-exp(-gp*(z+z1)))*tmp2 &
                                 /(1.d0-exp(-4.d0*gp*z1))
              else if (esm_bc.eq.'bc3') then
                 vg2(iz)=vg2(iz)-4.d0*pi*exp(gp*(z-z1))*tmp1 &
                                +2.d0*pi*(exp(gp*(z-z0-2.d0*z1))-exp(-gp*(z+z0)))*tmp2
              endif
           enddo

        vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)  ! Since cft_1z is not in-place
        call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,-1,vg2)

        vg3(i,j,1:dfftp%nr3x)=vg2(1:dfftp%nr3x)*2.d0
        endif
     enddo
  enddo


!****For gp=0 case ********************
  tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
  !for smoothing
  f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
  nz_l=dfftp%nr3x/2+1+esm_nfit
  nz_r=dfftp%nr3x/2+1-esm_nfit
  z_l=dble(nz_l-1)*L/dble(dfftp%nr3x)-L
  z_r=dble(nz_r-1)*L/dble(dfftp%nr3x)
  !
  if (esm_bc.eq.'bc1') then
     vg2(1)=-2.d0*pi*z0**2*rhog3(1,1,1)
  else if (esm_bc.eq.'bc2') then
     vg2(1)= 2.d0*pi*(2.d0*z1-z0)*z0*rhog3(1,1,1)
  else if (esm_bc.eq.'bc3') then
     vg2(1)= 2.d0*pi*(4.d0*z1-z0)*z0*rhog3(1,1,1)
  endif
  do iz=2,dfftp%nr3x
     if(iz<=dfftp%nr3x/2) kn=dble(iz-1)     *2.d0*pi/L
     if(iz> dfftp%nr3x/2) kn=dble(iz-1-dfftp%nr3x) *2.d0*pi/L
     cc=cos(kn*z0)
     ss=sin(kn*z0)
     if (esm_bc.eq.'bc1') then
        tmp1=tmp1+rhog3(1,1,iz)*ci*(cc+ci*ss)/kn
        tmp2=tmp2+rhog3(1,1,iz)*ci*(cc-ci*ss)/kn
        tmp3=tmp3+rhog3(1,1,iz)*cc/kn**2
        tmp4=tmp4+(0.d0,0.d0)
     else if (esm_bc.eq.'bc2') then
        tmp1=tmp1+rhog3(1,1,iz)*(cc+ci*ss)/kn**2
        tmp2=tmp2+rhog3(1,1,iz)*(cc-ci*ss)/kn**2
        tmp3=tmp3+rhog3(1,1,iz)*ci*cc/kn
        tmp4=tmp4+rhog3(1,1,iz)*ss/kn
     else if (esm_bc.eq.'bc3') then
        tmp1=tmp1+rhog3(1,1,iz)*(cc+ci*ss)/kn**2
        tmp2=tmp2+rhog3(1,1,iz)*(cc-ci*ss)/kn
        tmp3=tmp3+rhog3(1,1,iz)*(cc+ci*ss)/kn
        tmp4=tmp4+(0.d0,0.d0)
     endif
     vg2(iz)=4.d0*pi*rhog3(1,1,iz)/(kn**2)
     !for smoothing
     c_r=cos(kn*z_r)
     s_r=sin(kn*z_r)
     c_l=cos(kn*z_l)
     s_l=sin(kn*z_l)
     f1=f1+4.d0*pi*   rhog3(1,1,iz)*(c_r+ci*s_r)/kn**2
     f2=f2+4.d0*pi*   rhog3(1,1,iz)*(c_l+ci*s_l)/kn**2
     f3=f3+4.d0*pi*ci*rhog3(1,1,iz)*(c_r+ci*s_r)/kn
     f4=f4+4.d0*pi*ci*rhog3(1,1,iz)*(c_l+ci*s_l)/kn
     !
  enddo

  vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)  ! Since cft_1z is not in-place
  call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,1,vg2)
  
  do iz=1,dfftp%nr3x
     k3=iz-1
     if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
     z=dble(k3)/dble(dfftp%nr3x)*L
     if (esm_bc.eq.'bc1') then
        vg2(iz)=vg2(iz)-2.d0*pi*z**2*rhog3(1,1,1) &
                       -2.d0*pi*(z-z0)*tmp1       &
                       -2.d0*pi*(z+z0)*tmp2       &
                       -4.d0*pi*tmp3              
     else if (esm_bc.eq.'bc2') then
        vg2(iz)=vg2(iz)-2.d0*pi*z**2*rhog3(1,1,1)    &
                       -2.d0*pi*(z+z1)*tmp1/z1       &
                       +2.d0*pi*(z-z1)*tmp2/z1       &
                       -4.d0*pi*z*(z1-z0)/z1*tmp3    &
                       +4.d0*pi*(z1-z0)*tmp4               
     else if (esm_bc.eq.'bc3') then
        vg2(iz)=vg2(iz)-2.d0*pi*(z**2+2.d0*z*z0)*rhog3(1,1,1) &
                       -4.d0*pi*tmp1                          &
                       -4.d0*pi*ci*(z-z0)*tmp2                &
                       -4.d0*pi*ci*(z1-z0)*tmp3              
     endif
  enddo
  !for smoothing
  if (esm_bc.eq.'bc1') then
     f1=f1-2.d0*pi*z_r**2*rhog3(1,1,1) &
          -2.d0*pi*(z_r-z0)*tmp1 &
          -2.d0*pi*(z_r+z0)*tmp2 &
          -4.d0*pi*tmp3
     f1=f1-2.d0*pi*z0**2*rhog3(1,1,1)
     f2=f2-2.d0*pi*z_l**2*rhog3(1,1,1) &
          -2.d0*pi*(z_l-z0)*tmp1 &
          -2.d0*pi*(z_l+z0)*tmp2 &
          -4.d0*pi*tmp3
     f2=f2-2.d0*pi*z0**2*rhog3(1,1,1)
     f3=f3-2.d0*pi*tmp1-2.d0*pi*tmp2-4.d0*pi*z_r*rhog3(1,1,1)
     f4=f4-2.d0*pi*tmp1-2.d0*pi*tmp2-4.d0*pi*z_l*rhog3(1,1,1)
  else if (esm_bc.eq.'bc2') then
     f1=f1-2.d0*pi*z_r**2*rhog3(1,1,1) &
          -2.d0*pi*(z_r+z1)*tmp1/z1 &
          +2.d0*pi*(z_r-z1)*tmp2/z1 &
          -4.d0*pi*z*(z1-z0)/z1*tmp3 &
          +4.d0*pi  *(z1-z0)   *tmp4
     f1=f1+2.d0*pi*(2.d0*z1-z0)*z0*rhog3(1,1,1)
     f2=f2-2.d0*pi*z_l**2*rhog3(1,1,1) &
          -2.d0*pi*(z_l+z1)*tmp1/z1 &
          +2.d0*pi*(z_l-z1)*tmp2/z1 &
          -4.d0*pi*z*(z1-z0)/z1*tmp3 &
          +4.d0*pi  *(z1-z0)   *tmp4
     f2=f2+2.d0*pi*(2.d0*z1-z0)*z0*rhog3(1,1,1)
     f3=f3-4.d0*pi*z_r*rhog3(1,1,1)-2.d0*pi*tmp1/z1+2.d0*pi*tmp2/z1-4.d0*pi*(z1-z0)/z1*tmp3
     f4=f4-4.d0*pi*z_l*rhog3(1,1,1)-2.d0*pi*tmp1/z1+2.d0*pi*tmp2/z1-4.d0*pi*(z1-z0)/z1*tmp3
  else if (esm_bc.eq.'bc3') then
     f1=f1-2.d0*pi*(z_r**2+2.d0*z_r*z0)*rhog3(1,1,1) &
          -4.d0*pi*tmp1 &
          -4.d0*pi*ci*(z_r-z1)*tmp2 &
          -4.d0*pi*ci*(z1 -z0)*tmp3
     f1=f1+2.d0*pi*(4.d0*z1-z0)*z0*rhog3(1,1,1)
     f2=f2-2.d0*pi*(z_l**2+2.d0*z_l*z0)*rhog3(1,1,1) &
          -4.d0*pi*tmp1 &
          -4.d0*pi*ci*(z_l-z1)*tmp2 &
          -4.d0*pi*ci*(z1 -z0)*tmp3
     f2=f2+2.d0*pi*(4.d0*z1-z0)*z0*rhog3(1,1,1)
     f3=f3-2.d0*pi*(2.d0*z_r+2.d0*z0)*rhog3(1,1,1)-4.d0*pi*ci*tmp2
     f4=f4-2.d0*pi*(2.d0*z_l+2.d0*z0)*rhog3(1,1,1)-4.d0*pi*ci*tmp2
  endif
  ! for smoothing
  !factor 2 will be multiplied later (at vg3 <= vg2)
  !f1=f1*2.d0; f2=f2*2.d0; f3=f3*2.d0; f4=f4*2.d0
  z_r=z_r
  z_l=z_l+L
  a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
     +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
  a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
    -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
  a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
     +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
  a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
  do iz=nz_r,nz_l
    z=dble(iz-1)/dble(dfftp%nr3x)*L
    vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
  enddo

  vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)  ! Since cft_1z is not in-place
  call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,-1,vg2)

  vg3(1,1,1:dfftp%nr3x)=vg2(1:dfftp%nr3x)*2.d0

  ehart=0.d0
  do k=1,dfftp%nr3x
     do j=1,dfftp%nr2x
        do i=1,dfftp%nr1x
           ehart=ehart+vg3(i,j,k)*conjg(rhog3(i,j,k))*omega*0.5d0
        enddo
     enddo
  enddo
 
! Map to FFT mesh (nrxx)
  aux=0.0d0
  do ng=1,ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      IF (n1<1) n1 = n1 + dfftp%nr1
      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      IF (n2<1) n2 = n2 + dfftp%nr2
      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1 
      IF (n3<1) n3 = n3 + dfftp%nr3    
#if defined (__PARA) && !defined (__USE_3D_FFT)
     ijk=n3+(dfftp%isind(n1+(n2-1)*dfftp%nr1x)-1)*dfftp%nr3x
#else
     ijk=n1+(n2-1)*dfftp%nr1x+(n3-1)*dfftp%nr1x*dfftp%nr2x
#endif
     aux(ijk)= aux(ijk) + vg3(n1,n2,n3)
  enddo

  deallocate (vg3)
  deallocate (vg2,vg2_in,rhog3)

  RETURN
  END SUBROUTINE esm_hartree

!-----------------------------------------------------------------------
!--------------ESM EWALD SUBROUTINE-------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_ewald ( charge, alpha, ewg)

  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE ions_base,        ONLY : nat, tau, ityp, ntyp=>nsp
  USE uspp_param,       ONLY : upf
  USE fft_base,         ONLY : dfftp

  implicit none
  REAL(DP)                :: charge, alpha, ewg
  !
  !    here the local variables
  !
  real(DP), external      :: qe_erfc, qe_erf
  real(DP)                :: gp2, t(2), gp, sa, z1, z0, L
  integer                 :: i, j, k, k1, k2, k3, ipol, it1, it2, ik, kg(2,dfftp%nr1x*dfftp%nr2x), nkg
  real(DP) :: tt, z, zp, kk1, kk2, g, cc1, cc2, arg1, arg2, t1, t2, ff, argmax

  argmax=0.9*log(huge(1.d0))
  ewg=0.d0
  ik=0
  do i = 1, dfftp%nr1x
     k1=i-1
     if (i.gt.dfftp%nr1x/2) k1=i-dfftp%nr1x-1
     do j = 1, dfftp%nr2x
        k2=j-1
        if (j.gt.dfftp%nr2x/2) k2=j-dfftp%nr2x-1
        ik=ik+1
        kg(1,ik)=k1
        kg(2,ik)=k2
     enddo
  enddo
  nkg=ik

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  g=sqrt(alpha)
  sa=omega/L

  do it1=1,nat
  do 1 it2=1,nat

     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat

     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat

     tt=upf(ityp(it1))%zp*upf(ityp(it2))%zp*2.0*pi/sa

     kk1=0.5d0*(-(z-zp)*qe_erf(g*(z-zp))-exp(-g**2*(z-zp)**2)/g/sqrt(pi))

     if (esm_bc.eq.'bc1') then
        kk2=0.d0
     else if (esm_bc.eq.'bc2') then
        kk2=0.5d0*(z1-z*zp/z1)
     else if (esm_bc.eq.'bc3') then
        kk2=0.5d0*(2.d0*z1-z-zp)
     endif

     if (it1.eq.it2) then

        cc1=0.d0
        cc2=0.d0

        do ik = 2, nkg 

           gp2 = 0.d0
           do ipol = 1, 2
              t (ipol) = kg(1,ik) * bg (ipol, 1) + kg(2,ik) * bg (ipol, 2)
              gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
           enddo
           gp=sqrt(gp2)

           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           cc1=cc1+(t1+t2)/4.d0/gp

           if (esm_bc.eq.'bc1') then
              cc2=0.d0
           else if (esm_bc.eq.'bc2') then
              cc2=cc2+(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                      -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                      /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp
           else if (esm_bc.eq.'bc3') then
              cc2=cc2-exp(gp*(z+zp-2.d0*z1))/2.d0/gp
           endif            
        enddo
        ewg=ewg+tt*(kk1+kk2+cc1+cc2)
        goto  1

     else if (it1.gt.it2) then

        cc1=0.d0
        cc2=0.d0

        do ik = 2, nkg

           gp2 = 0.d0
           do ipol = 1, 2
              t (ipol) = kg(1,ik) * bg (ipol, 1) + kg(2,ik) * bg (ipol, 2)
              gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
           enddo
           gp=sqrt(gp2)

           ff = ( ( kg(1,ik)*bg(1,1)+kg(2,ik)*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( kg(1,ik)*bg(2,1)+kg(2,ik)*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * 2.d0 * pi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp

           if (esm_bc.eq.'bc1') then
              cc2=0.d0
           else if (esm_bc.eq.'bc2') then
              cc2=cc2+cos(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                              -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                              /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp
           else if (esm_bc.eq.'bc3') then
              cc2=cc2+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
           endif
        enddo
        ewg=ewg+tt*(kk1+kk2+cc1+cc2)*2.0

        goto 1

     endif
1 continue
  enddo

  ewg=2.0*ewg

     do it1=1,nat
        ewg=ewg- upf(ityp(it1))%zp **2 * sqrt (8.d0 / tpi * alpha)
     enddo

  return
end subroutine esm_ewald


!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL SUBROUTINE---------------------------
!-----------------------------------------------------------------------
subroutine esm_local (aux)

  USE kinds,            ONLY : DP
  USE constants,        ONLY : pi, eps8
  USE gvect,            ONLY : g, ngm
  USE cell_base,        ONLY : at, bg, alat, tpiba2, tpiba, omega
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE scf,              ONLY : rho
  USE lsda_mod,         ONLY : nspin
  USE fft_scalar,       ONLY : cft_1z
  USE fft_base,         ONLY : dfftp
  USE mp_global,        ONLY : me_pool
  !
  implicit none
  COMPLEX(DP)             :: aux( dfftp%nnr )     ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  !
  !    here the local variables
  !
  complex(DP),allocatable :: vloc3(:,:,:),vg2(:),vg2_in(:)
  real(DP),allocatable    :: rhog(:,:),bgauss(:,:)
  real(DP), external      :: qe_erf, qe_erfc
  real(DP)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2, &
                             z,zp,arg11,arg12,arg21,arg22,v0,tmp,L,argmax, &
                             z_l,z_r
  integer                 :: iz,ig,it,ijk,i,j,k,ipol,k1,k2,k3,ng,n1,n2,n3, &
                             nz_l,nz_r
  complex(DP)             :: cs,cc1,cc2,ci,a0,a1,a2,a3,f1,f2,f3,f4

  argmax=0.9*log(huge(1.d0))
  L =at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)

  allocate(vloc3(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x),vg2(dfftp%nr3x),vg2_in(dfftp%nr3x),bgauss(nat,1))
  do it=1,nat
     bgauss(it,1)=1.d0
  enddo
  sa=omega/L
  v0=esm_efield*z1*2.d0/2.d0 ! factor 1/2: unit Ry. -> hartree
  ci=(0.d0,1.d0)

! for gp!=0
  do i = 1, dfftp%nr1x
     k1=i-1
     if (i.gt.dfftp%nr1x/2) k1=i-dfftp%nr1x-1
     do j = 1, dfftp%nr2x
        k2=j-1
        if (j.gt.dfftp%nr2x/2) k2=j-dfftp%nr2x-1
        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = k1 * bg (ipol, 1) + k2 * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)
        vg2(1:dfftp%nr3x)=(0.d0,0.d0)
        if (gp.ge.eps8) then
           do it=1,nat
              tt=-4.d0*pi*upf(ityp(it))%zp/sa
              pp=-2.d0*pi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
              cc=cos(pp)
              ss=sin(pp)
              cs=CMPLX ( cc, ss, kind=DP )
              zp=tau(3,it)
              if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
              zp=zp*alat
              do iz=1,dfftp%nr3x
                 k3=iz-1
                 if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
                 z=dble(k3)/dble(dfftp%nr3x)*L
                 cc1=(0.d0,0.d0)
                 do ig=1,1
                    tmp=1.d0
                    arg11=-gp*(z-zp)
                    arg11=min(arg11,argmax)
                    arg12= gp/2.d0/tmp-tmp*(z-zp)
                    arg21= gp*(z-zp)
                    arg21=min(arg21,argmax)
                    arg22= gp/2.d0/tmp+tmp*(z-zp)
                    t1=exp(arg11)*qe_erfc(arg12)
                    t2=exp(arg21)*qe_erfc(arg22)
                    cc1=cc1+bgauss(it,ig)*cs*(t1+t2)/4.d0/gp
                 enddo
                 if (esm_bc.eq.'bc1') then
                    cc2=(0.d0,0.d0)
                 else if (esm_bc.eq.'bc2') then
                    cc2=cs*( exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                            -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                            /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp 
                 else if (esm_bc.eq.'bc3') then
                    cc2=cs*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
                 endif
                 vg2(iz) = vg2(iz) + tt*(cc1+cc2)*2.d0 ! factor 2: hartree -> Ry.
              enddo       
           enddo
           vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)
           call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,-1,vg2)
           do iz=1,dfftp%nr3x
              vloc3(i,j,iz)=vg2(iz)
           enddo
        endif
     enddo
  enddo
  vg2(1:dfftp%nr3x)=(0.d0,0.d0)
! for smoothing
  f1=0.d0; f2=0.d0; f3=0.d0; f4=0.d0
  nz_l=dfftp%nr3x/2+1+esm_nfit
  nz_r=dfftp%nr3x/2+1-esm_nfit
  z_l=dble(nz_l-1)*L/dble(dfftp%nr3x)-L
  z_r=dble(nz_r-1)*L/dble(dfftp%nr3x)
! add constant potential (capacitor term)
  do iz=1,dfftp%nr3x
     k3=iz-1
     if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
     z=dble(k3)/dble(dfftp%nr3x)*L
     vg2(iz)=-0.5d0*v0*(z-z1)/z1*2.d0 ! factor 2: hartree -> Ry.
  enddo
  f1=-0.5d0*v0*(z_r-z1)/z1 ! unit: hartree
  f2=-0.5d0*v0*(z_l-z1)/z1 ! unit: hartree
  f3=-0.5d0*v0/z1 ! unit: hartree/a.u.
  f4=-0.5d0*v0/z1 ! unit: harteee/a.u.
! for gp=0
  do it=1,nat
     tt=-4.d0*pi*upf(ityp(it))%zp/sa
     zp=tau(3,it)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     do iz=1,dfftp%nr3x
        k3=iz-1
        if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
        z=dble(k3)/dble(dfftp%nr3x)*L
        cc1=(0.d0,0.d0) 
        do ig=1,1
           tmp=1.d0
           cc1=cc1+bgauss(it,ig)*0.5d0*(-(z-zp)*qe_erf(tmp*(z-zp)) &
                                        -exp(-tmp**2*(z-zp)**2)/tmp/sqrt(pi))
        enddo
        if (esm_bc.eq.'bc1') then
           cc2=(0.d0,0.d0)
        else if (esm_bc.eq.'bc2') then
           cc2=0.5d0*(z1-z*zp/z1)
        else if (esm_bc.eq.'bc3') then
           cc2=0.5d0*(2.d0*z1-z-zp)
        endif
        vg2(iz) = vg2(iz) + tt*(cc1+cc2)*2.d0 ! factor 2: hartree -> Ry.
     enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
     do ig=1,1
        tmp=1.d0
        f1=f1+tt*bgauss(it,ig)*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
                                -exp(-tmp**2*(z_r-zp)**2)/tmp/sqrt(pi))
        f2=f2+tt*bgauss(it,ig)*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
                                -exp(-tmp**2*(z_l-zp)**2)/tmp/sqrt(pi))
        f3=f3-tt*bgauss(it,ig)*0.5d0*qe_erf(tmp*(z_r-zp))
        f4=f4-tt*bgauss(it,ig)*0.5d0*qe_erf(tmp*(z_l-zp))
     enddo
     if(esm_bc.eq.'bc1')then
        f1=f1+tt*0.d0
        f2=f2+tt*0.d0
        f3=f3+tt*0.d0
        f4=f4+tt*0.d0
     elseif(esm_bc.eq.'bc2')then
        f1=f1+tt*0.5d0*(z1-z_r*zp/z1)
        f2=f2+tt*0.5d0*(z1-z_l*zp/z1)
        f3=f3+tt*(-0.5d0*(zp/z1))
        f4=f4+tt*(-0.5d0*(zp/z1))
     elseif(esm_bc.eq.'bc3')then
        f1=f1+tt*0.5d0*(2.d0*z1-z_r-zp)
        f2=f2+tt*0.5d0*(2.d0*z1-z_l-zp)
        f3=f3-tt*0.5d0
        f4=f4-tt*0.5d0
     endif
  enddo
  ! for smoothing
  f1=f1*2.d0; f2=f2*2.d0; f3=f3*2.d0; f4=f4*2.d0 ! factor 2: hartree -> Ry.
  z_r=z_r
  z_l=z_l+L
  a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
     +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
  a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
    -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
  a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
     +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
  a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
  do iz=nz_r,nz_l
    z=dble(iz-1)/dble(dfftp%nr3x)*L
    vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
  enddo
  vg2_in(1:dfftp%nr3x)=vg2(1:dfftp%nr3x)
  call cft_1z(vg2_in,1,dfftp%nr3x,dfftp%nr3x,-1,vg2)
  do iz=1,dfftp%nr3x
     vloc3(1,1,iz)=vg2(iz)
  enddo
  deallocate(vg2,vg2_in,bgauss)

! Map to FFT mesh (nrxx)
  do ng=1,ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      IF (n1<1) n1 = n1 + dfftp%nr1
      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      IF (n2<1) n2 = n2 + dfftp%nr2
      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1 
      IF (n3<1) n3 = n3 + dfftp%nr3    

#if defined (__PARA) && !defined (__USE_3D_FFT)
     ijk=n3+(dfftp%isind(n1+(n2-1)*dfftp%nr1x)-1)*dfftp%nr3x
#else
     ijk=n1+(n2-1)*dfftp%nr1x+(n3-1)*dfftp%nr1x*dfftp%nr2x
#endif
     aux(ijk)= aux(ijk)+ vloc3(n1,n2,n3)
  enddo

  deallocate(vloc3)

  return
  end subroutine esm_local



!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE SUBROUTINE-----------------------
!-----------------------------------------------------------------------
subroutine esm_force_ew ( alpha, forceion ) 

  USE kinds
  USE constants,        ONLY : pi, tpi, e2, eps8
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE fft_base,         ONLY : dfftp

  implicit none
  REAL(DP)                :: alpha
  REAL(DP)                :: forceion(3,nat) 
  !
  !    here the local variables
  !
  real(DP), external      :: qe_erfc, qe_erf
  integer  :: it1, it2, ik, kg(2,dfftp%nr1x*dfftp%nr2x), nkg, ipol, i, j, k, k1, k2, k3
  real(DP) :: tt_for, z, zp, kk1_for, kk2_for, g, for_g(3, nat), gp2, gp, z1, t(2), L
  real(DP) :: cx1_for, cy1_for, cz1_for, cx2_for, cy2_for, cz2_for, arg1, arg2, t1, t2, ff, sa, z0
  real(DP) :: g_b,tauz1,tauz2,gt,tt,gz,argmax

  argmax=0.9*log(huge(1.d0))
  for_g(:,:)=0.d0
  ik=0
  do i = 1, dfftp%nr1x
     k1=i-1
     if (i.gt.dfftp%nr1x/2) k1=i-dfftp%nr1x-1
     do j = 1, dfftp%nr2x
        k2=j-1
        if (j.gt.dfftp%nr2x/2) k2=j-dfftp%nr2x-1
        ik=ik+1
        kg(1,ik)=k1
        kg(2,ik)=k2
     enddo
  enddo
  nkg=ik

  forceion(:,:)=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w) 
  sa=omega/L
  g=sqrt(alpha)

  do  it1=1,nat
  do 1 it2=1,nat

   z=tau(3,it1)
   if (z.gt.at(3,3)*0.5) z=z-at(3,3)
   z=z*alat
   zp=tau(3,it2)
   if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
   zp=zp*alat

   tt_for=upf(ityp(it1))%zp*upf(ityp(it2))%zp*2.d0*2.d0*pi/sa
   kk1_for=0.5d0*qe_erf(g*(z-zp))
   if (esm_bc.eq.'bc1') then
     kk2_for=0.d0
   else if (esm_bc.eq.'bc2') then
     kk2_for=-0.5d0*(z/z1)
   else if (esm_bc.eq.'bc3') then
     kk2_for=-0.5d0
   endif

  
   if (it1.eq.it2) then
     cz1_for=0.d0
     cz2_for=0.d0
     do ik = 2, nkg 
        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = kg(1,ik) * bg (ipol, 1) + kg(2,ik) * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)

        arg1=-gp*(z-zp)
        arg2= gp*(z-zp)
        arg1=min(arg1,argmax)
        arg2=min(arg2,argmax)
        t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
        t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
        cz1_for=0.d0
        if (esm_bc.eq.'bc1') then      
           cz2_for=0.d0
        else if (esm_bc.eq.'bc2') then      
           cz2_for=cz2_for - (exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                             +exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                             /(1.d0-exp(-4.d0*gp*z1))/2.d0
        else if (esm_bc.eq.'bc3') then      
           cz2_for=cz2_for - exp(gp*(z+zp-2.d0*z1))/2.d0
        endif
     enddo
     for_g(3,it2) = for_g(3,it2) &
                  + tt_for*(cz1_for+cz2_for+kk1_for+kk2_for)

     goto 1

   else if (it1.gt.it2) then

     cx1_for=0.d0
     cy1_for=0.d0
     cz1_for=0.d0
     cx2_for=0.d0
     cy2_for=0.d0
     cz2_for=0.d0
     do ik = 2, nkg

        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = kg(1,ik) * bg (ipol, 1) + kg(2,ik) * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)
        ff = ( ( kg(1,ik)*bg(1,1)+kg(2,ik)*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
           +   ( kg(1,ik)*bg(2,1)+kg(2,ik)*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * 2.d0 * pi
        arg1=-gp*(z-zp)
        arg2= gp*(z-zp)
        arg1=min(arg1,argmax)
        arg2=min(arg2,argmax)
        t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
        t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))

        cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*kg(1,ik)
        cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*kg(2,ik)
        cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
        if (esm_bc.eq.'bc1') then
           cx2_for=0.d0
           cy2_for=0.d0
           cz2_for=0.d0
        else if (esm_bc.eq.'bc2') then
           cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                    - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*kg(1,ik)
           cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                    - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*kg(2,ik)
           cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                    + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0
        else if (esm_bc.eq.'bc3') then
           cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*kg(1,ik)
           cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*kg(2,ik)
           cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
        endif
     enddo
     for_g(1,it2)=for_g(1,it2)+tt_for*(cx1_for+cx2_for)
     for_g(2,it2)=for_g(2,it2)+tt_for*(cy1_for+cy2_for)
     for_g(3,it2)=for_g(3,it2) &
                 +tt_for*(cz1_for+cz2_for+kk1_for+kk2_for)
    
     goto 1
   else if (it1.lt.it2) then

     cx1_for=0.d0
     cy1_for=0.d0
     cz1_for=0.d0
     cx2_for=0.d0
     cy2_for=0.d0
     cz2_for=0.d0
     do ik = 2, nkg

        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = kg(1,ik) * bg (ipol, 1) + kg(2,ik) * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)
        
        ff = ( ( kg(1,ik)*bg(1,1)+kg(2,ik)*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
           +   ( kg(1,ik)*bg(2,1)+kg(2,ik)*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * 2.d0 * pi
        arg1=-gp*(z-zp)
        arg2= gp*(z-zp)
        arg1=min(arg1,argmax)
        arg2=min(arg2,argmax)
        t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
        t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))

        cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*kg(1,ik)
        cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*kg(2,ik)
        cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
        if (esm_bc.eq.'bc1') then
           cx2_for=0.d0
           cy2_for=0.d0
           cz2_for=0.d0
        else if (esm_bc.eq.'bc2') then
           cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                    - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*kg(1,ik)
           cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                    - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*kg(2,ik)
           cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                    + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                    /(1.d0-exp(-4.d0*gp*z1))/2.d0
        else if (esm_bc.eq.'bc3') then
           cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*kg(1,ik)
           cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*kg(2,ik)
           cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
        endif
     enddo
     for_g(1,it2)=for_g(1,it2)+tt_for*(cx1_for+cx2_for)
     for_g(2,it2)=for_g(2,it2)+tt_for*(cy1_for+cy2_for)
     for_g(3,it2)=for_g(3,it2) &
                 +tt_for*(cz1_for+cz2_for+kk1_for+kk2_for)

     goto 1

   endif

 1 continue
   enddo

   for_g(:,:)=2.0*for_g(:,:)

   do it1=1,nat
      do ipol=1,2
         forceion(ipol,it1)=(for_g(1,it1)*bg(ipol,1)+for_g(2,it1)*bg(ipol,2))*sqrt(tpiba2)
      enddo
      forceion(3,it1)=for_g(3,it1)
   enddo
   forceion(:,:)=-forceion(:,:)

  return
end subroutine esm_force_ew

!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL-DERIVED FORCE SUBROUTINE-------------
!-----------------------------------------------------------------------
subroutine esm_force_lc ( aux, forcelc )

  USE kinds
  USE constants,        ONLY : pi, eps8
  USE gvect,            ONLY : g, ngm
  USE cell_base,        ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_pool_comm
  USE fft_scalar,       ONLY : cft_1z
  USE fft_base,         ONLY : dfftp

  implicit none
  COMPLEX(DP)             :: aux(dfftp%nnr)       ! aux contains n(G) (input)   
  REAL(DP)                :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  complex(DP),allocatable :: vlocx(:), vlocy(:), vlocdz(:)
  real(DP),allocatable    :: bgauss(:,:),for_tmp(:),for(:,:)
  real(DP), external      :: qe_erf, qe_erfc
  real(DP)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2,z,zp,L,forcelc2(3,nat)
  real(DP)                :: arg11,arg12,arg21,arg22,tmp,r1,r2,fx1,fy1,fz1,fx2,fy2,fz2,argmax
  integer                 :: iz,ig,it,ijk,i,j,k,ipol,k1,k2,k3,ng,n1,n2,n3
  complex(DP),allocatable :: vg2(:),vg2_fx(:),vg2_fy(:),vg2_fz(:),rhog3(:,:,:)
  complex(DP)             :: cx1,cy1,cz1,cx2,cy2,cz2,cc1,cc2

  argmax=0.9*log(huge(1.d0))

! Mat to FULL FFT mesh (dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  allocate(rhog3(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x))
  rhog3(:,:,:)=(0.d0,0.d0)
  do ng=1,ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      IF (n1<1) n1 = n1 + dfftp%nr1
      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      IF (n2<1) n2 = n2 + dfftp%nr2
      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1 
      IF (n3<1) n3 = n3 + dfftp%nr3    
#ifdef __PARA
     ijk=n3+(dfftp%isind(n1+(n2-1)*dfftp%nr1x)-1)*dfftp%nr3x
#else
     ijk=n1+(n2-1)*dfftp%nr1x+(n3-1)*dfftp%nr1x*dfftp%nr2x
#endif
     rhog3(n1,n2,n3)=aux(ijk)
  enddo
#ifdef __PARA
  call mp_sum( rhog3, intra_pool_comm )
#endif

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)

  allocate(vg2(dfftp%nr3x),vg2_fx(dfftp%nr3x),vg2_fy(dfftp%nr3x),vg2_fz(dfftp%nr3x),bgauss(nat,1))
  allocate(for(3,nat))
  do it=1,nat
     bgauss(it,1)=1.d0
  enddo
  sa=omega/L
  for(:,:)=0.d0
  vg2_fx(:)=(0.d0,0.d0)
  vg2_fy(:)=(0.d0,0.d0)
  vg2_fz(:)=(0.d0,0.d0)

!**** for gp!=0 *********
  do i = 1, dfftp%nr1x
     k1=i-1
     if (i.gt.dfftp%nr1x/2) k1=i-dfftp%nr1x-1
     do j = 1, dfftp%nr2x
        k2=j-1
        if (j.gt.dfftp%nr2x/2) k2=j-dfftp%nr2x-1
        gp2 = 0.d0
        do ipol = 1, 2
           t (ipol) = k1 * bg (ipol, 1) + k2 * bg (ipol, 2)
           gp2  = gp2 + t (ipol) * t (ipol) * tpiba2
        enddo
        gp=sqrt(gp2)
        if (gp.ge.eps8) then
           do it=1,nat
              tt=-4.d0*pi*upf(ityp(it))%zp/sa
              pp=-2.d0*pi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
              cc=cos(pp)
              ss=sin(pp)
              zp=tau(3,it)
              if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
              zp=zp*alat
              do iz=1,dfftp%nr3x
                 k3=iz-1
                 if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
                 z=dble(k3)/dble(dfftp%nr3x)*L
                 cx1=(0.d0,0.d0); cy1=(0.d0,0.d0); cz1=(0.d0,0.d0)
                 do ig=1,1
                    tmp=1.d0
                    arg11=-gp*(z-zp)
                    arg11=min(arg11,argmax) 
                    arg12= gp/2.d0/tmp-tmp*(z-zp)
                    arg21= gp*(z-zp)
                    arg21=min(arg21,argmax)
                    arg22= gp/2.d0/tmp+tmp*(z-zp)
                    t1=exp(arg11)*qe_erfc(arg12)
                    t2=exp(arg21)*qe_erfc(arg22)
                    cx1=cx1+bgauss(it,ig)*CMPLX(ss, -cc, kind=DP) &
                        *(t1+t2)/4.d0/gp*k1
                    cy1=cy1+bgauss(it,ig)*CMPLX(ss, -cc, kind=DP) &
                        *(t1+t2)/4.d0/gp*k2
                    cz1=cz1+bgauss(it,ig)*CMPLX(cc, ss, kind=DP)  &
                        *(t1-t2)/4.d0
                 enddo
                 if (esm_bc.eq.'bc1') then
                    cx2=(0.d0,0.d0)
                    cy2=(0.d0,0.d0)
                    cz2=(0.d0,0.d0)
                 else if (esm_bc.eq.'bc2') then
                    cx2=CMPLX(ss, -cc, kind=DP)* &
                       (exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                       -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
                    cy2=CMPLX(ss, -cc, kind=DP)* &
                       (exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                       -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
                    cz2=CMPLX(cc, ss, kind=DP)* &
                       (-exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                        -exp(gp*(z+zp-2.d0*z1))+exp(-gp*(z+zp+2.d0*z1))) &
                       /(1.d0-exp(-4.d0*gp*z1))/2.d0
                 else if (esm_bc.eq.'bc3') then
                    cx2=CMPLX(ss, -cc, kind=DP)* &
                       (-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
                    cy2=CMPLX(ss, -cc, kind=DP)* &
                       (-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
                    cz2=CMPLX(cc, ss, kind=DP)* &
                       (-exp(gp*(z+zp-2.d0*z1)))/2.d0
                 endif
                 vg2_fx(iz) = tt*(cx1+cx2)
                 vg2_fy(iz) = tt*(cy1+cy2)
                 vg2_fz(iz) = tt*(cz1+cz2)
              enddo
              vg2(1:dfftp%nr3x)=vg2_fx(1:dfftp%nr3x)  ! Since cft_1z is not in-place
              call cft_1z(vg2,1,dfftp%nr3x,dfftp%nr3x,-1,vg2_fx)
              vg2(1:dfftp%nr3x)=vg2_fy(1:dfftp%nr3x)  ! Since cft_1z is not in-place
              call cft_1z(vg2,1,dfftp%nr3x,dfftp%nr3x,-1,vg2_fy)
              vg2(1:dfftp%nr3x)=vg2_fz(1:dfftp%nr3x)  ! Since cft_1z is not in-place
              call cft_1z(vg2,1,dfftp%nr3x,dfftp%nr3x,-1,vg2_fz)
              do iz=1,dfftp%nr3x
                 r1= dble(rhog3(i,j,iz))
                 r2=aimag(rhog3(i,j,iz))
                 fx1=dble(  vg2_fx(iz))
                 fy1=dble(  vg2_fy(iz))
                 fz1=dble(  vg2_fz(iz))
                 fx2=aimag( vg2_fx(iz))
                 fy2=aimag( vg2_fy(iz))
                 fz2=aimag( vg2_fz(iz))
                 for(1,it)=for(1,it)-r1*fx1-r2*fx2
                 for(2,it)=for(2,it)-r1*fy1-r2*fy2
                 for(3,it)=for(3,it)-r1*fz1-r2*fz2
              enddo

           enddo
        endif
     enddo
  enddo

!***** for gp==0********
  i=1
  j=1
  allocate(for_tmp(nat))
  for_tmp(:)=0.d0
  vg2_fz(:)=(0.d0,0.d0)
  
  do it=1,nat
     tt=-4.d0*pi*upf(ityp(it))%zp/sa
     zp=tau(3,it)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     do iz=1,dfftp%nr3x
        k3=iz-1
        if (k3.gt.dfftp%nr3x/2) k3=iz-dfftp%nr3x-1
        z=dble(k3)/dble(dfftp%nr3x)*L
        cc1=(0.d0,0.d0)
        do ig=1,1
           tmp=1.d0
           cc1=cc1+bgauss(it,ig)*(0.5d0*qe_erf(tmp*(z-zp)))
        enddo
        if (esm_bc.eq.'bc1') then
           cc2=(0.d0,0.d0)
        else if (esm_bc.eq.'bc2') then
           cc2=-0.5d0*(z/z1)
        else if (esm_bc.eq.'bc3') then
           cc2=-0.5d0
        endif
        vg2_fz(iz) =  tt*(cc1+cc2)

     enddo
     vg2(1:dfftp%nr3x)=vg2_fz(1:dfftp%nr3x)  ! Since cft_1z is not in-place
     call cft_1z(vg2,1,dfftp%nr3x,dfftp%nr3x,-1,vg2_fz)
     do iz=1,dfftp%nr3x
        r1=dble( rhog3(1,1,iz))
        r2=aimag(rhog3(1,1,iz))
        fz1=dble( vg2_fz(iz))
        fz2=aimag(vg2_fz(iz))
        for(3,it)=for(3,it)-r1*fz1-r2*fz2
     enddo
  enddo
  deallocate(vg2,vg2_fx,vg2_fy,vg2_fz,bgauss)

!***** sum short_range part and long_range part in local potential force at cartecian coordinate

  do it=1,nat
     do ipol=1,2
        forcelc(ipol,it)=forcelc(ipol,it)+(for(1,it)*bg(ipol,1)+for(2,it)*bg(ipol,2))*sqrt(tpiba2)*omega*2.d0
     enddo
     forcelc(3,it)=forcelc(3,it)+for(3,it)*omega*2.d0
  enddo

  deallocate(for)

  call setlocal()

  deallocate(rhog3)
  return
end subroutine esm_force_lc

!-----------------------------------------------------------------------
!--------------ESM FINAL PRINTOUT SUBROUTINE----------------------------
!-----------------------------------------------------------------------
!
! Prints out vlocal and vhartree to stdout once electrons are converged
! Format: z, rho(r), v_hartree, v_local, (v_hartree + v_local) 
!
SUBROUTINE esm_printpot ()
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, alat, omega
  USE gvect,                ONLY : ngm
  USE scf,                  ONLY : rho, vltot
  USE lsda_mod,             ONLY : nspin
  USE constants,            ONLY : eps4
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : intra_pool_comm, me_pool
  USE fft_base,             ONLY : grid_gather, dfftp
  USE io_global,            ONLY : ionode, stdout
  !
  IMPLICIT NONE
  !
  REAL(DP)                :: z1,z2,z3,z4,charge,ehart,bohr,rydv,L,area
  REAL(DP), ALLOCATABLE   :: vh(:),work1(:),work2(:),work3(:)
  INTEGER                 :: ix,iy,iz,i,k3

        allocate(vh(dfftp%nnr))
        allocate(work1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
        allocate(work2(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
        allocate(work3(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
        vh(:)=0.d0; work1(:)=0.d0; work2(:)=0.d0; work3(:)=0.d0
        L=alat*at(3,3)
        area=(at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2
        bohr=0.52917720859d0
        rydv=13.6058d0
        CALL v_h (rho%of_g, ehart, charge, vh)
#ifdef __PARA
        call grid_gather( vh, work2 )
        call mp_sum(work2, intra_pool_comm)
        call grid_gather( vltot, work3 )
        call mp_sum(work3, intra_pool_comm)
#else
        work2(1:dfftp%nnr)=vh(1:dfftp%nnr)
        work3(1:dfftp%nnr)=vltot(1:dfftp%nnr)
#endif
        if( nspin == 2 ) then
          vh(:)=rho%of_r(:,1)+rho%of_r(:,2)
        else
          vh(:)=rho%of_r(:,1)
        endif
#ifdef __PARA
        call grid_gather( vh, work1 )
        call mp_sum(work1, intra_pool_comm)
#else
        work1(1:dfftp%nnr)=vh(1:dfftp%nnr)
#endif
        deallocate(vh)
        IF ( ionode ) then
        write(stdout,                                                 &
              FMT = '(/,5x, "ESM Charge and Potential",&
                     &/,5x, "========================",/)' )
        write(stdout, 9051)
        write(stdout, 9052)
! z = position along slab (A)
! rho = planar-summed charge density of slab section (e)
! v_hartree = planar-averaged hartree potential term (eV/A)
! v_local = planar-averaged local potential term (eV/A)
        do iz=dfftp%nr3/2+2,dfftp%nr3
           k3=iz-1-dfftp%nr3
           z1=0.d0;z2=0.d0;z3=0.d0;z4=0.d0
           do iy=1,dfftp%nr2
           do ix=1,dfftp%nr1
              i=ix+(iy-1)*dfftp%nr1+(iz-1)*dfftp%nr1*dfftp%nr2
              z1=z1+work1(i)*area/dble(dfftp%nr1*dfftp%nr2)
              z2=z2+(work2(i)+work3(i))/dble(dfftp%nr1*dfftp%nr2)
              z3=z3+work2(i)/dble(dfftp%nr1*dfftp%nr2)
              z4=z4+work3(i)/dble(dfftp%nr1*dfftp%nr2)
           enddo
           enddo
           write(stdout,'(f9.3,f13.5,2f19.7,f18.7)') &
           dble(k3)/dble(dfftp%nr3)*L*bohr, z1, z3*rydv/bohr, &
             z4*rydv/bohr, z2*rydv/bohr
        enddo
        !do iz=1,dfftp%nr3x/2+1-esm_nfit
        do iz=1,dfftp%nr3/2+1
           k3=iz-1
           z1=0.d0;z2=0.d0;z3=0.d0;z4=0.d0
           do iy=1,dfftp%nr2
           do ix=1,dfftp%nr1
              i=ix+(iy-1)*dfftp%nr1+(iz-1)*dfftp%nr1*dfftp%nr2
              z1=z1+work1(i)*area/dble(dfftp%nr1*dfftp%nr2)
              z2=z2+(work2(i)+work3(i))/dble(dfftp%nr1*dfftp%nr2)
              z3=z3+work2(i)/dble(dfftp%nr1*dfftp%nr2)
              z4=z4+work3(i)/dble(dfftp%nr1*dfftp%nr2)
           enddo
           enddo
           write(stdout,'(f9.3,f13.5,2f19.7,f18.7)') &
           dble(k3)/dble(dfftp%nr3)*L*bohr, z1, z3*rydv/bohr, &
             z4*rydv/bohr, z2*rydv/bohr
        enddo
        deallocate(work1,work2,work3)
        write(stdout,*) 
        ENDIF
9051    FORMAT( 4x,'z (A)',6x,'rho (e)',6x,'Avg v_hartree',8x,&
        &'Avg v_local',2x,'Avg v_hart+v_loc' )
9052    FORMAT(35x,'(eV/A)',13x,'(eV/A)',12x,'(eV/A)',/,4x,& 
 &'==========================================================================' )
END SUBROUTINE esm_printpot
!
!-----------------------------------------------------------------------
!--------------ESM SUMMARY PRINTOUT SUBROUTINE--------------------------
!-----------------------------------------------------------------------
!
! Prints summary of ESM parameters to stdout
!
SUBROUTINE esm_summary ()
      !
      USE io_global,      ONLY : stdout,                               &
                                 ionode
      !
      IMPLICIT NONE
      !
       WRITE( UNIT = stdout,                                          &
              FMT  = '(/,5x, "Effective Screening Medium Method",     &
                      &/,5x, "=================================")' )
      !
      WRITE( UNIT = stdout, FMT = 9051 ) esm_efield
      !
      WRITE( UNIT = stdout, FMT = 9052 ) esm_w
      !
      WRITE( UNIT = stdout, FMT = 9053 ) esm_nfit
      !
      IF( ionode ) THEN
        !
        SELECT CASE( TRIM( esm_bc ) )
        !
        CASE( 'pbc' )
           WRITE( UNIT = stdout,                                     &
             FMT  = '(5x, "Ordinary Periodic Boundary Conditions")' )
        CASE( 'bc1' )
           WRITE( UNIT = stdout,                                     &
             FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-Vacuum")' )
        CASE( 'bc2' )
           WRITE( UNIT = stdout,                                     &
             FMT  = '(5x, "Boundary Conditions: Metal-Slab-Metal")' )
        CASE( 'bc3' )
           WRITE( UNIT = stdout,                                     &
             FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-Metal")' )
        END SELECT
      END IF
      !
      WRITE( stdout, * )
      !
9051 FORMAT( '     field strength (Ry/a.u.)         = ',  F10.2,' ')
9052 FORMAT( '     ESM offset from cell edge (a.u.) = ',  F10.2,' ' )
9053 FORMAT( '     grid points for fit at edges     = ',  I10,' ')

END SUBROUTINE esm_summary

END MODULE esm
