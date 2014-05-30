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
  USE constants, ONLY : pi, tpi, fpi, eps4, eps8, e2
  SAVE
  !
  LOGICAL :: do_comp_esm=.FALSE.
  INTEGER :: esm_nfit
  REAL(KIND=DP) :: esm_efield, esm_w                 
  CHARACTER (LEN=3) :: esm_bc           
  INTEGER, ALLOCATABLE, TARGET :: mill_2d(:,:), imill_2d(:,:)
  INTEGER :: ngm_2d = 0
  !
  PUBLIC :: esm_hartree, esm_local, esm_ewald, esm_force_lc, esm_force_ew, &
            esm_printpot, esm_summary, esm_ggen_2d, esm_deallocate_gvect_2d

CONTAINS

SUBROUTINE esm_deallocate_gvect_2d
  IF( ALLOCATED( mill_2d ) ) DEALLOCATE( mill_2d )
  RETURN
END SUBROUTINE esm_deallocate_gvect_2d

SUBROUTINE esm_ggen_2d()
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, mill
  USE control_flags,    ONLY : gamma_only
  USE fft_scalar,       ONLY : cft_1z
  !
  IMPLICIT NONE
  !
  INTEGER :: n1xh, n2xh, ng, n1, n2, ng_2d
  Logical, ALLOCATABLE :: do_mill_2d(:,:)
  COMPLEX(DP), ALLOCATABLE :: vg2_in(:), vg2(:)
  !
  !     Make g parallel array
  !
  n1xh = dfftp%nr1x/2
  n2xh = dfftp%nr2x/2
  ALLOCATE( do_mill_2d(-n1xh:n1xh,-n2xh:n2xh) )
  do_mill_2d(:,:) = .false.
  
  DO ng = 1, ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     do_mill_2d(n1,n2) = .true.
  ENDDO
  ngm_2d = COUNT( do_mill_2d )
!*** do_mill_2d(h,k) = .true. means there is an h,k vector on this proc
!*** ngm_2d = total number of vectors (h,k) on this proc, excluding duplicates
!*** with different l values
  
  ALLOCATE( mill_2d(2,ngm_2d), imill_2d(-n1xh:n1xh,-n2xh:n2xh) )
  mill_2d(:,:) = 0
  imill_2d(:,:) = 0
  ng_2d = 1
  DO n1 = -n1xh, n1xh
  DO n2 = -n2xh, n2xh
     IF( do_mill_2d(n1,n2) ) THEN
        mill_2d(1,ng_2d) = n1
        mill_2d(2,ng_2d) = n2
        imill_2d(n1,n2) = ng_2d
        ng_2d = ng_2d + 1
     ENDIF
  ENDDO
  ENDDO
  DEALLOCATE(do_mill_2d)  
!**** mill_2d(:,ig) = h,k indices of vector ig
!**** imill_2d(h,k) = 2d index of vector with h,k indices
!**** ng_2d = total number of 2d g vectors on this proc

  RETURN
END SUBROUTINE esm_ggen_2d

!
!-----------------------------------------------------------------------
!--------------ESM HARTREE SUBROUTINE-----------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_hartree (rhog, ehart, aux)

  USE gvect,     ONLY : g, nl, nlm, ngm, mill
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE fft_scalar,       ONLY : cft_1z
  USE mp_bands,         ONLY : intra_bgrp_comm
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
                              z_l, z_r, eh
  integer                  :: ipol, k, k1, k2, k3, iz, ng, n1, n2, n3, &
                              nz_r, nz_l, ng_2d
  complex(DP),allocatable  :: rhog3(:,:), vg2(:), vg2_in(:), vg3(:,:)
  complex(DP)              :: xc, ci, tmp, tmp1, tmp2, tmp3, tmp4, f1, f2, f3, f4, &
                              a0, a1, a2, a3, c_r, c_l, s_r, s_l, rg3

  allocate(vg2(dfftp%nr3),vg2_in(dfftp%nr3),rhog3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + dfftp%nr3    
     if (nspin == 2) then
        rg3 = rhog(ng,1)+rhog(ng,2)
     else
        rg3 = rhog(ng,1)
     endif
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo
! End mapping
!
  allocate(vg3(dfftp%nr3,ngm_2d))
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel do private( k1, k2, gp2, ipol, t, gp, tmp1, tmp2, vg2, iz, kn, &
!$omp                      cc, ss, tmp, vg2_in, k3, z, rg3 )
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0)
     vg2(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)     * tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3) * tpi/L
        cc=cos(kn*z0)
        ss=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        vg2(iz)=fpi*rg3/(gp**2+kn**2)
        if (esm_bc.eq.'bc1') then
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp-ci*kn)
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp+ci*kn)
        else if (esm_bc.eq.'bc2') then
           tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp**2+kn**2)*tmp
           tmp=((gp-ci*kn)*exp(gp*(z1-z0))+(gp+ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp**2+kn**2)*tmp
        else if (esm_bc.eq.'bc3') then
           tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp**2+kn**2)*tmp
           tmp=(gp-ci*kn)/gp
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp**2+kn**2)*tmp
        endif
     enddo
     
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,1,vg2)
     
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        if (esm_bc.eq.'bc1') then
           vg2(iz)=vg2(iz)-tpi/gp*(exp(gp*(z-z0))*tmp1+exp(-gp*(z+z0))*tmp2)
        else if (esm_bc.eq.'bc2') then
           vg2(iz)=vg2(iz)-fpi*(exp(gp*(z-z1))-exp(-gp*(z+3.d0*z1)))*tmp1 &
                         /(1.d0-exp(-4.d0*gp*z1)) &
                          +fpi*(exp(gp*(z-3.d0*z1))-exp(-gp*(z+z1)))*tmp2 &
                         /(1.d0-exp(-4.d0*gp*z1))
        else if (esm_bc.eq.'bc3') then
           vg2(iz)=vg2(iz)-fpi*exp(gp*(z-z1))*tmp1 &
                +tpi*(exp(gp*(z-z0-2.d0*z1))-exp(-gp*(z+z0)))*tmp2
        endif
     enddo
     
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,-1,vg2)
     
     vg3(1:dfftp%nr3,ng_2d)=vg2(1:dfftp%nr3)*2.d0
  enddo
  
!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     !for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     !
     rg3=rhog3(1,ng_2d)
     if (esm_bc.eq.'bc1') then
        vg2(1)=-tpi*z0**2*rg3
     else if (esm_bc.eq.'bc2') then
        vg2(1)= tpi*(2.d0*z1-z0)*z0*rg3
     else if (esm_bc.eq.'bc3') then
        vg2(1)= tpi*(4.d0*z1-z0)*z0*rg3
     endif
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)     *tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3) *tpi/L
        cc=cos(kn*z0)
        ss=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        if (esm_bc.eq.'bc1') then
           tmp1=tmp1+rg3*ci*(cc+ci*ss)/kn
           tmp2=tmp2+rg3*ci*(cc-ci*ss)/kn
           tmp3=tmp3+rg3*cc/kn**2
           tmp4=tmp4+(0.d0,0.d0)
        else if (esm_bc.eq.'bc2') then
           tmp1=tmp1+rg3*(cc+ci*ss)/kn**2
           tmp2=tmp2+rg3*(cc-ci*ss)/kn**2
           tmp3=tmp3+rg3*ci*cc/kn
           tmp4=tmp4+rg3*ss/kn
        else if (esm_bc.eq.'bc3') then
           tmp1=tmp1+rg3*(cc+ci*ss)/kn**2
           tmp2=tmp2+rg3*(cc-ci*ss)/kn
           tmp3=tmp3+rg3*(cc+ci*ss)/kn
           tmp4=tmp4+(0.d0,0.d0)
        endif
        vg2(iz)=fpi*rg3/(kn**2)
        !for smoothing
        c_r=cos(kn*z_r)
        s_r=sin(kn*z_r)
        c_l=cos(kn*z_l)
        s_l=sin(kn*z_l)
        f1=f1+fpi*   rg3*(c_r+ci*s_r)/kn**2
        f2=f2+fpi*   rg3*(c_l+ci*s_l)/kn**2
        f3=f3+fpi*ci*rg3*(c_r+ci*s_r)/kn
        f4=f4+fpi*ci*rg3*(c_l+ci*s_l)/kn
        !
     enddo
     
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,1,vg2)
     
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        if (esm_bc.eq.'bc1') then
           vg2(iz)=vg2(iz)-tpi*z**2*rg3    &
                          -tpi*(z-z0)*tmp1 &
                          -tpi*(z+z0)*tmp2 &
                          -fpi*tmp3              
        else if (esm_bc.eq.'bc2') then
           vg2(iz)=vg2(iz)-tpi*z**2*rg3          &
                          -tpi*(z+z1)*tmp1/z1    &
                          +tpi*(z-z1)*tmp2/z1    &
                          -fpi*z*(z1-z0)/z1*tmp3 &
                          +fpi*(z1-z0)*tmp4               
        else if (esm_bc.eq.'bc3') then
           vg2(iz)=vg2(iz)-tpi*(z**2+2.d0*z*z0)*rg3 &
                          -fpi*tmp1                 &
                          -fpi*ci*(z-z1)*tmp2       &
                          -fpi*ci*(z1-z0)*tmp3
        endif
     enddo
     !for smoothing
     if (esm_bc.eq.'bc1') then
        f1=f1-tpi*z_r**2*rg3 &
             -tpi*(z_r-z0)*tmp1 &
             -tpi*(z_r+z0)*tmp2 &
             -fpi*tmp3
        f1=f1-tpi*z0**2*rg3
        f2=f2-tpi*z_l**2*rg3 &
             -tpi*(z_l-z0)*tmp1 &
             -tpi*(z_l+z0)*tmp2 &
             -fpi*tmp3
        f2=f2-tpi*z0**2*rg3
        f3=f3-tpi*tmp1-tpi*tmp2-fpi*z_r*rg3
        f4=f4-tpi*tmp1-tpi*tmp2-fpi*z_l*rg3
     else if (esm_bc.eq.'bc2') then
        f1=f1-tpi*z_r**2*rg3 &
             -tpi*(z_r+z1)*tmp1/z1 &
             +tpi*(z_r-z1)*tmp2/z1 &
             -fpi*z_r*(z1-z0)/z1*tmp3 &
             +fpi    *(z1-z0)   *tmp4
        f1=f1+tpi*(2.d0*z1-z0)*z0*rg3
        f2=f2-tpi*z_l**2*rg3 &
             -tpi*(z_l+z1)*tmp1/z1 &
             +tpi*(z_l-z1)*tmp2/z1 &
             -fpi*z_l*(z1-z0)/z1*tmp3 &
             +fpi    *(z1-z0)   *tmp4
        f2=f2+tpi*(2.d0*z1-z0)*z0*rg3
        f3=f3-fpi*z_r*rg3-tpi*tmp1/z1+tpi*tmp2/z1-fpi*(z1-z0)/z1*tmp3
        f4=f4-fpi*z_l*rg3-tpi*tmp1/z1+tpi*tmp2/z1-fpi*(z1-z0)/z1*tmp3
     else if (esm_bc.eq.'bc3') then
        f1=f1-tpi*(z_r**2+2.d0*z_r*z0)*rg3 &
             -fpi*tmp1 &
             -fpi*ci*(z_r-z1)*tmp2 &
             -fpi*ci*(z1 -z0)*tmp3
        f1=f1+tpi*(4.d0*z1-z0)*z0*rg3
        f2=f2-tpi*(z_l**2+2.d0*z_l*z0)*rg3 &
             -fpi*tmp1 &
             -fpi*ci*(z_l-z1)*tmp2 &
             -fpi*ci*(z1 -z0)*tmp3
        f2=f2+tpi*(4.d0*z1-z0)*z0*rg3
        f3=f3-tpi*(2.d0*z_r+2.d0*z0)*rg3-fpi*ci*tmp2
        f4=f4-tpi*(2.d0*z_l+2.d0*z0)*rg3-fpi*ci*tmp2
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
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,-1,vg2)
     
     vg3(1:dfftp%nr3,ng_2d)=vg2(1:dfftp%nr3)*2.d0
     
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, k, eh )
  eh = 0d0
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
!$omp atomic
  ehart=ehart+eh
!$omp end parallel
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart *omega*0.5d0
  !
  call mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
  enddo
  if (gamma_only) then
     do ng=1,ngm
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     enddo
  endif 

  deallocate (vg3)
  deallocate (vg2,vg2_in,rhog3)

  RETURN
END SUBROUTINE esm_hartree

!-----------------------------------------------------------------------
!--------------ESM EWALD SUBROUTINE-------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_ewald ( charge, alpha, ewg )

  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE ions_base,        ONLY : nat, tau, ityp, ntyp=>nsp
  USE uspp_param,       ONLY : upf
  USE fft_base,         ONLY : dfftp
  USE control_flags,    ONLY : gamma_only

  implicit none
  REAL(DP)                :: charge, alpha, ewg
  !
  !    here the local variables
  !
  real(DP), external      :: qe_erfc, qe_erf
  real(DP)                :: gp2, t(2), gp, sa, z1, z0, L
  integer                 :: k1, k2, k3, ipol, it1, it2, ng_2d
  real(DP) :: tt, z, zp, kk1, kk2, g, cc1, cc2, arg1, arg2, t1, t2, ff, argmax, ew
#ifdef __OPENMP
  INTEGER :: nth, ith, omp_get_thread_num, omp_get_num_threads
#endif

  argmax=0.9*log(huge(1.d0))
  ewg=0.d0

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  g=sqrt(alpha)
  sa=omega/L
#ifdef __OPENMP
!$omp parallel private( nth, ith, ew, it1, it2, z, zp, tt, kk1, kk2, cc1, cc2, &
!$omp                   ng_2d, k1, k2, gp2, ipol, t, gp, ff, arg1, arg2, t1, t2 )
#endif
#ifdef __OPENMP
  nth=omp_get_num_threads()
  ith=omp_get_thread_num()
#endif
  ew=0d0
  do it1=1,nat
  do it2=1,it1
#ifdef __OPENMP
     if( mod( (it1-1)*it1/2+it2-1, nth) /= ith ) cycle
#endif

     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat

     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     
     tt=upf(ityp(it1))%zp*upf(ityp(it2))%zp*tpi/sa

     kk1=0.5d0*(-(z-zp)*qe_erf(g*(z-zp))-exp(-g**2*(z-zp)**2)/g/sqrt(pi))

     if (esm_bc.eq.'bc1') then
        kk2=0.d0
     else if (esm_bc.eq.'bc2') then
        kk2=0.5d0*(z1-z*zp/z1)
     else if (esm_bc.eq.'bc3') then
        kk2=0.5d0*(2.d0*z1-z-zp)
     endif

     cc1=0.d0
     cc2=0.d0
     
     if (it1.eq.it2) then

        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if( k1==0 .and. k2==0 ) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
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
              cc2=cc2+(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
           endif
        
        enddo

        if( gamma_only ) then
           cc1 = cc1 * 2d0
           cc2 = cc2 * 2d0
        endif
        ew=ew+tt*(cc1+cc2)
        if(gstart==2) ew=ew+tt*(kk1+kk2)

     else

        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if( k1==0 .and. k2==0 ) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
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

        if( gamma_only ) then
           cc1 = cc1 * 2d0
           cc2 = cc2 * 2d0
        endif
        ew=ew+tt*(cc1+cc2)*2d0
        if(gstart==2) ew=ew+tt*(kk1+kk2)*2d0

     endif
  enddo
  enddo
!$omp atomic
  ewg=ewg+ew
#ifdef __OPENMP
!$omp end parallel
#endif

  ewg=2.0*ewg

  if( gstart == 2 ) then
     do it1=1,nat
        ewg=ewg- upf(ityp(it1))%zp **2 * sqrt (8.d0 / tpi * alpha)
     enddo
  endif

  return
end subroutine esm_ewald


!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL SUBROUTINE---------------------------
!-----------------------------------------------------------------------
subroutine esm_local (aux)

  USE kinds,            ONLY : DP
  USE gvect,            ONLY : g, ngm, nl, nlm, mill
  USE control_flags,    ONLY : gamma_only
  USE cell_base,        ONLY : at, bg, alat, tpiba2, tpiba, omega
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE scf,              ONLY : rho
  USE lsda_mod,         ONLY : nspin
  USE fft_scalar,       ONLY : cft_1z
  USE fft_base,         ONLY : dfftp
  !
  implicit none
  COMPLEX(DP)             :: aux( dfftp%nnr )     ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  !
  !    here the local variables
  !
  complex(DP),allocatable :: vloc3(:,:),vg2(:),vg2_in(:)
  real(DP),allocatable    :: rhog(:,:),bgauss(:,:)
  real(DP), external      :: qe_erf, qe_erfc
  real(DP)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2, &
                             z,zp,arg11,arg12,arg21,arg22,v0,tmp,L,argmax, &
                             z_l,z_r
  integer                 :: iz,ig,it,ipol,k1,k2,k3,ng,n1,n2,n3, &
                             nz_l,nz_r, ng_2d
  complex(DP)             :: cs,cc1,cc2,ci,a0,a1,a2,a3,f1,f2,f3,f4

  argmax=0.9*log(huge(1.d0))
  L =at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)

  allocate(vloc3(dfftp%nr3,ngm_2d),vg2(dfftp%nr3),vg2_in(dfftp%nr3),bgauss(nat,1))
  do it=1,nat
     bgauss(it,1)=1.d0
  enddo
  sa=omega/L
  v0=esm_efield*z1*2.d0/2.d0 ! factor 1/2: unit Ry. -> hartree
  ci=(0.d0,1.d0)

! for gp!=0
!$omp parallel do private( k1, k2, gp2, gp, vg2, it, tt, pp, cc, ss, cs, zp, iz, &
!$omp                      k3, z, cc1, ig, tmp, arg11, arg12, arg21, arg22, t1, t2, &
!$omp                      cc2, vg2_in )
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
        
     vg2(1:dfftp%nr3)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*upf(ityp(it))%zp/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX ( cc, ss, kind=DP )
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
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
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,-1,vg2)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg2(iz)
     enddo
  enddo
  
  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then

     vg2(1:dfftp%nr3)=(0.d0,0.d0)
! for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
! add constant potential (capacitor term)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg2(iz)=-0.5d0*v0*(z-z1)/z1*2.d0 ! factor 2: hartree -> Ry.
     enddo
     f1=-0.5d0*v0*(z_r-z1)/z1 ! unit: hartree
     f2=-0.5d0*v0*(z_l-z1)/z1 ! unit: hartree
     f3=-0.5d0*v0/z1 ! unit: hartree/a.u.
     f4=-0.5d0*v0/z1 ! unit: harteee/a.u.
! for gp=0
     do it=1,nat
        tt=-fpi*upf(ityp(it))%zp/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
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
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     vg2_in(1:dfftp%nr3)=vg2(1:dfftp%nr3)
     call cft_1z(vg2_in,1,dfftp%nr3,dfftp%nr3,-1,vg2)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg2(iz)
     enddo
     
  endif ! if( ng_2d > 0 )
  deallocate(vg2,vg2_in,bgauss)
  
! Map to FFT mesh (dfftp%nrx)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
  enddo
  if (gamma_only) then
     do ng=1,ngm
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     enddo
  endif

  deallocate(vloc3)

  return
  end subroutine esm_local



!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE SUBROUTINE-----------------------
!-----------------------------------------------------------------------
subroutine esm_force_ew ( alpha, forceion ) 

  USE kinds
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : gstart

  implicit none
  REAL(DP)                :: alpha
  REAL(DP)                :: forceion(3,nat) 
  !
  !    here the local variables
  !
  real(DP), external      :: qe_erfc, qe_erf
  integer  :: it1, it2, ipol, k1, k2, k3, ng_2d
  integer  :: nth, ith, omp_get_num_threads, omp_get_thread_num
  real(DP) :: t1_for, t2_for, z, zp, kk1_for, kk2_for, g, for_g(3, nat), gp2, gp, z1, t(2), L
  real(DP) :: cx1_for, cy1_for, cz1_for, cx2_for, cy2_for, cz2_for, arg1, arg2, t1, t2, ff
  real(DP) :: sa, z0, g_b,tauz1,tauz2,gt,tt,gz,argmax,for(3, nat)

  argmax=0.9*log(huge(1.d0))
  for_g(:,:)=0.d0
  
  forceion(:,:)=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w) 
  sa=omega/L
  g=sqrt(alpha)

!$omp parallel private( nth, ith, for, z, zp, t1_for, t2_for, kk1_for, kk2_for, &
!$omp                   cz1_for, cz2_for, ng_2d, k1, k2, gp2, gp, arg1, arg2, t1, t2, &
!$omp                   cx1_for, cy1_for, cx2_for, cy2_for, ff )
#ifdef __OPENMP
  nth=omp_get_num_threads()
  ith=omp_get_thread_num()
#endif
  for=0d0
  do it1=1,nat
  do it2=1,nat
#ifdef __OPENMP
     if( mod( (it1-1)*nat+it2-1, nth) /= ith ) cycle
#endif

     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=upf(ityp(it1))%zp*upf(ityp(it2))%zp*fpi/sa*2.d0
     else
        t1_for=upf(ityp(it1))%zp*upf(ityp(it2))%zp*fpi/sa
     endif
     t2_for=upf(ityp(it1))%zp*upf(ityp(it2))%zp*fpi/sa

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
        
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if(k1==0.and.k2==0) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
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
        for(3,it2) = for(3,it2) + t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2) = for(3,it2) + t2_for*(kk1_for+kk2_for)
        endif
        
     else if (it1.gt.it2) then
        
        cx1_for=0.d0
        cy1_for=0.d0
        cz1_for=0.d0
        cx2_for=0.d0
        cy2_for=0.d0
        cz2_for=0.d0
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if(k1==0.and.k2==0) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           
           cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*k1
           cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*k2
           cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
           if (esm_bc.eq.'bc1') then
              cx2_for=0.d0
              cy2_for=0.d0
              cz2_for=0.d0
           else if (esm_bc.eq.'bc2') then
              cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
              cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
              cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                       + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then
              cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
              cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
              cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
           endif
        enddo
        for(1,it2)=for(1,it2)+t1_for*(cx1_for+cx2_for)
        for(2,it2)=for(2,it2)+t1_for*(cy1_for+cy2_for)
        for(3,it2)=for(3,it2)+t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
        endif
        
     else if (it1.lt.it2) then
        
        cx1_for=0.d0
        cy1_for=0.d0
        cz1_for=0.d0
        cx2_for=0.d0
        cy2_for=0.d0
        cz2_for=0.d0
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if(k1==0.and.k2==0) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           
           cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*k1
           cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*k2
           cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
           if (esm_bc.eq.'bc1') then
              cx2_for=0.d0
              cy2_for=0.d0
              cz2_for=0.d0
           else if (esm_bc.eq.'bc2') then
              cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
              cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
              cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                       + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then
              cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
              cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
              cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
           endif
        enddo
        for(1,it2)=for(1,it2)+t1_for*(cx1_for+cx2_for)
        for(2,it2)=for(2,it2)+t1_for*(cy1_for+cy2_for)
        for(3,it2)=for(3,it2)+t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
        endif
     endif
     
  enddo
  enddo
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
!$omp end parallel

  for_g(:,:)=2.0*for_g(:,:)

  do it1=1,nat
     forceion(1,it1)=sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
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
  USE gvect,            ONLY : g, ngm, nl, nlm, mill
  USE cell_base,        ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : nat, tau, ityp
  USE uspp_param,       ONLY : upf
  USE fft_scalar,       ONLY : cft_1z
  USE fft_base,         ONLY : dfftp

  implicit none
  COMPLEX(DP)             :: aux(dfftp%nnr)       ! aux contains n(G) (input)   
  REAL(DP)                :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  real(DP),allocatable    :: bgauss(:,:),for(:,:),for_g(:,:)
  real(DP), external      :: qe_erf, qe_erfc
  real(DP)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2,z,zp,L
  real(DP)                :: arg11,arg12,arg21,arg22,tmp,r1,r2,fx1,fy1,fz1,fx2,fy2,fz2,argmax
  integer                 :: iz,ig,it,ipol,k1,k2,k3,ng,n1,n2,n3,ng_2d
  complex(DP),allocatable :: vg2(:),vg2_fx(:),vg2_fy(:),vg2_fz(:),rhog3(:,:)
  complex(DP)             :: cx1,cy1,cz1,cx2,cy2,cz2,cc1,cc2

  argmax=0.9*log(huge(1.d0))

! Map to FULL FFT mesh (dfftp%nr1x,dfftp%nr2x,dfftp%nr3)
  allocate(rhog3(dfftp%nr3,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
      n1 = mill(1,ng)
      n2 = mill(2,ng)
      ng_2d = imill_2d(n1,n2)
      n3 = mill(3,ng) + 1
      IF (n3<1) n3 = n3 + dfftp%nr3
      rhog3(n3,ng_2d)=aux(nl(ng))
      if (gamma_only .and. n1==0 .and. n2==0) then
         n3 = -mill(3,ng)+1
         IF(n3<1)n3=n3+dfftp%nr3
         rhog3(n3,ng_2d)=aux(nlm(ng))
      endif  
  enddo

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)

  allocate(vg2(dfftp%nr3),vg2_fx(dfftp%nr3),vg2_fy(dfftp%nr3),vg2_fz(dfftp%nr3),bgauss(nat,1))
  allocate(for_g(3,nat))
  do it=1,nat
     bgauss(it,1)=1.d0
  enddo
  sa=omega/L
  for_g(:,:)=0.d0
  vg2_fx(:)=(0.d0,0.d0)
  vg2_fy(:)=(0.d0,0.d0)
  vg2_fz(:)=(0.d0,0.d0)

!**** for gp!=0 *********
!$omp parallel private( k1, k2, gp2, gp, it, tt, pp, cc, ss, zp, iz, &
!$omp                   k3, z, cx1, cy1, cz1, tmp, arg11, arg12, arg21, arg22, &
!$omp                   t1, t2, cx2, cy2, cz2, vg2_fx, vg2_fy, vg2_fz, vg2, &
!$omp                   r1, r2, fx1, fy1, fz1, fx2, fy2, fz2, for )
  allocate(for(3,nat))
  for(:,:)=0.d0
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        IF (gamma_only) THEN
           tt=-fpi*upf(ityp(it))%zp/sa*2.d0
        ELSE 
           tt=-fpi*upf(ityp(it))%zp/sa
        ENDIF 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
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
        vg2(1:dfftp%nr3)=vg2_fx(1:dfftp%nr3)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,dfftp%nr3,dfftp%nr3,-1,vg2_fx)
        vg2(1:dfftp%nr3)=vg2_fy(1:dfftp%nr3)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,dfftp%nr3,dfftp%nr3,-1,vg2_fy)
        vg2(1:dfftp%nr3)=vg2_fz(1:dfftp%nr3)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,dfftp%nr3,dfftp%nr3,-1,vg2_fz)
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
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
  enddo
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
  deallocate(for)
!$omp end critical
!$omp end parallel


!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then

     vg2_fz(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*upf(ityp(it))%zp/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
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
        vg2(1:dfftp%nr3)=vg2_fz(1:dfftp%nr3)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,dfftp%nr3,dfftp%nr3,-1,vg2_fz)
        do iz=1,dfftp%nr3
           r1=dble( rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           fz1=dble( vg2_fz(iz))
           fz2=aimag(vg2_fz(iz))
           for_g(3,it)=for_g(3,it)-r1*fz1-r2*fz2
        enddo
     enddo
     
  endif ! if( ng_2d > 0 )

  deallocate(vg2,vg2_fx,vg2_fy,vg2_fz,bgauss)

!***** sum short_range part and long_range part in local potential force at cartecian coordinate

  do it=1,nat
     forcelc(1,it)=forcelc(1,it)+sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*2.d0
     forcelc(2,it)=forcelc(2,it)+sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*2.d0
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*2.d0
  enddo

  deallocate(for_g)

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
  USE cell_base,            ONLY : at, alat
  USE scf,                  ONLY : rho, vltot
  USE lsda_mod,             ONLY : nspin
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode, stdout
  USE constants,            ONLY : rytoev, bohr_radius_angs
  !
  IMPLICIT NONE
  !
  REAL(DP)                :: z1,z2,z3,z4,charge,ehart,L,area
  REAL(DP), ALLOCATABLE   :: work1(:),work2(:,:),work3(:), work4(:,:)
  INTEGER                 :: ix,iy,iz,izz,i,k3

        allocate(work1(dfftp%nnr))
        allocate(work2(dfftp%nnr,nspin))
        allocate(work3(dfftp%nnr))
        allocate(work4(5,dfftp%nr3))
        work1(:)=0.d0; work2(:,:)=0.d0; work3(:)=0.d0; work4(:,:)=0.d0
        L=alat*at(3,3)
        area=(at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2
        CALL v_h (rho%of_g, ehart, charge, work2)
        work3(1:dfftp%nnr)=vltot(1:dfftp%nnr)
        if( nspin == 2 ) then
           work1(:)=rho%of_r(:,1)+rho%of_r(:,2)
        else
           work1(:)=rho%of_r(:,1)
        endif

! z = position along slab (A)
! rho = planar-summed charge density of slab section (e)
! v_hartree = planar-averaged hartree potential term (eV)
! v_local = planar-averaged local potential term (eV)

!$omp parallel do private( iz, izz, k3, z1, z2, z3, z4, iy, ix, i )
        do iz = 1, dfftp%npp(dfftp%mype+1)
           izz = iz + dfftp%ipp(dfftp%mype+1)
           k3 = izz - 1
           if( k3 > dfftp%nr3/2 ) k3 = k3 - dfftp%nr3
           z1=0.d0;z2=0.d0;z3=0.d0;z4=0.d0
           do iy=1,dfftp%nr2
           do ix=1,dfftp%nr1
              i=ix+(iy-1)*dfftp%nr1+(iz-1)*dfftp%nr1*dfftp%nr2
              z1=z1+work1(i)*area/dble(dfftp%nr1*dfftp%nr2)
              z2=z2+(work2(i,1)+work3(i))/dble(dfftp%nr1*dfftp%nr2)
              z3=z3+work2(i,1)/dble(dfftp%nr1*dfftp%nr2)
              z4=z4+work3(i)/dble(dfftp%nr1*dfftp%nr2)
           enddo
           enddo
           work4(1:5,izz) = (/dble(k3)/dble(dfftp%nr3)*L*bohr_radius_angs, &
                              z1/bohr_radius_angs, z3*rytoev,z4*rytoev, &
                              z2*rytoev/)
        enddo
        !
        call mp_sum(work4, intra_bgrp_comm)
        !
        IF ( ionode ) then
           write(stdout,                                             &
                 FMT = '(/,5x, "ESM Charge and Potential",&
                        &/,5x, "========================",/)' )
           write(stdout, 9051)
           write(stdout, 9052)
           do k3 = dfftp%nr3/2-dfftp%nr3+1, dfftp%nr3/2
              iz = k3 + dfftp%nr3 + 1
              if( iz > dfftp%nr3 ) iz = iz - dfftp%nr3
              write(stdout,'(f9.2,f12.4,2f19.7,f18.7)') work4(1:5,iz)
           enddo
           write(stdout,*) 
        ENDIF
        deallocate(work1,work2,work3,work4)
9051    FORMAT( 4x,'z (A)',3x,'Tot chg (e/A)',3x,'Avg v_hartree',8x,&
        &'Avg v_local',2x,'Avg v_hart+v_loc' )
9052    FORMAT(37x,'(eV)',15x,'(eV)',14x,'(eV)',/,4x,& 
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
