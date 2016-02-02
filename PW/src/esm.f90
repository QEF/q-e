!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Original version by Minoru Otani (AIST), Yoshio Miura (Tohoku U.),
! Nicephore Bonet (MIT), Nicola Marzari (MIT), Brandon Wood (LLNL), 
! Tadashi Ogitsu (LLNL).
! Constant bias potential (constant-mu) method by Minoru Otani (AIST) and
! Nicephore Bonnet (AIST).
!
! Contains SUBROUTINEs for implementation of
! 1) ESM (Effective Screening Medium Method) developed by M. Otani and 
!    O. Sugino (see PRB 73, 115407 [2006])
! 2) Constant-mu method developed by N. Bonnet, T. Morishita, O. Sugino, 
!    and M. Otani (see PRL 109, 266101 [2012]).
!
! ESM enables description of a surface slab sandwiched between two 
! semi-infinite media, making it possible to deal with polarized surfaces 
! without using dipole corrections. It is useful for simulating interfaces 
! with vacuum, one or more electrodes, or an electrolyte.
!
! Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
! description of the system is connected to a potentiostat which preserves
! the Fermi energy of the system as the target Fermi energy (mu).
!
! Modified SUBROUTINEs for calculating the Hartree potential, the local 
! potential, and the Ewald sum are contained here, along with SUBROUTINEs for
! calculating force contributions based on the modified local potential and 
! Ewald term. Constant-mu parts are contained in the fcp.f90.
!
!----------------------------------------------------------------------------
MODULE esm
  !--------------------------------------------------------------------------
  !
  ! ... this module contains the variables and SUBROUTINEs needed for the 
  ! ... EFFECTIVE SCREENING MEDIUM (ESM) METHOD 
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: do_comp_esm, esm_nfit, esm_efield, esm_w, esm_a, esm_bc, &
            mill_2d, imill_2d, ngm_2d, &
            esm_init, esm_hartree, esm_local, esm_ewald, esm_force_lc, &
            esm_force_ew, esm_printpot, esm_summary

  !
  LOGICAL              :: do_comp_esm=.FALSE.
  INTEGER              :: esm_nfit
  REAL(KIND=DP)        :: esm_efield, esm_w, esm_a
  CHARACTER (LEN=3)    :: esm_bc

  INTEGER, ALLOCATABLE :: mill_2d(:,:), imill_2d(:,:)
  INTEGER              :: ngm_2d = 0
  real(DP), external   :: qe_erf, qe_erfc
  !
  CONTAINS
     subroutine esm_hartree( rhog, ehart, aux )
        USE kinds,    ONLY : DP
        USE gvect,    ONLY : ngm
        USE lsda_mod, ONLY : nspin
        USE fft_base, ONLY : dfftp
        IMPLICIT NONE
        real(DP)    :: ehart             !  Hartree energy
        complex(DP) :: rhog(ngm,nspin)   !  n(G)
        complex(DP) :: aux(dfftp%nnr)    !  v_h(G)

        if( esm_bc == 'pbc' ) then
           call esm_hartree_pbc ( rhog, ehart, aux )
        else if ( esm_bc == 'bc1' ) then
           call esm_hartree_bc1 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc2' ) then
           call esm_hartree_bc2 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc3' ) then
           call esm_hartree_bc3 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc4' ) then
           call esm_hartree_bc4 ( rhog, ehart, aux )
        end if

 
     end subroutine esm_hartree

     SUBROUTINE esm_ewaldr( alpha_g, ewr )
        USE kinds,    ONLY : DP
        implicit none
        real(DP), intent(in)  :: alpha_g
        real(DP), intent(out) :: ewr

        if( esm_bc == 'pbc' ) then
           call esm_ewaldr_pbc ( alpha_g, ewr )
        else if ( esm_bc == 'bc1' ) then
           call esm_ewaldr_pbc ( alpha_g, ewr )
        else if ( esm_bc == 'bc2' ) then
           call esm_ewaldr_pbc ( alpha_g, ewr )
        else if ( esm_bc == 'bc3' ) then
           call esm_ewaldr_pbc ( alpha_g, ewr )
        else if ( esm_bc == 'bc4' ) then
           call esm_ewaldr_bc4 ( alpha_g, ewr )
        end if

     END SUBROUTINE esm_ewaldr

     SUBROUTINE esm_ewaldg( alpha_g, ewg )
        USE kinds,    ONLY : DP
        implicit none
        real(DP), intent(in)  :: alpha_g
        real(DP), intent(out) :: ewg

        if( esm_bc == 'pbc' ) then
           call esm_ewaldg_pbc ( alpha_g, ewg )
        else if ( esm_bc == 'bc1' ) then
           call esm_ewaldg_bc1 ( alpha_g, ewg )
        else if ( esm_bc == 'bc2' ) then
           call esm_ewaldg_bc2 ( alpha_g, ewg )
        else if ( esm_bc == 'bc3' ) then
           call esm_ewaldg_bc3 ( alpha_g, ewg )
        else if ( esm_bc == 'bc4' ) then
           call esm_ewaldg_bc4 ( alpha_g, ewg )
        end if

     END SUBROUTINE esm_ewaldg

     SUBROUTINE esm_local( aux )
        USE kinds,    ONLY : DP
        USE fft_base, ONLY : dfftp
        implicit none
        complex(DP), intent(inout) :: aux( dfftp%nnr )

        if( esm_bc == 'pbc' ) then
           call esm_local_pbc ( aux )
        else if ( esm_bc == 'bc1' ) then
           call esm_local_bc1 ( aux )
        else if ( esm_bc == 'bc2' ) then
           call esm_local_bc2 ( aux )
        else if ( esm_bc == 'bc3' ) then
           call esm_local_bc3 ( aux )
        else if ( esm_bc == 'bc4' ) then
           call esm_local_bc4 ( aux )
        end if

     END SUBROUTINE esm_local
     
     SUBROUTINE esm_force_ewr( alpha_g, forceion ) 
        USE kinds,     ONLY : DP
        USE ions_base, ONLY : nat
        implicit none
        real(DP),intent(in)    :: alpha_g
        real(DP),intent(inout) :: forceion(3,nat) 

        if( esm_bc == 'pbc' ) then
           call esm_force_ewr_pbc ( alpha_g, forceion )
        else if ( esm_bc == 'bc1' ) then
           call esm_force_ewr_pbc ( alpha_g, forceion )
        else if ( esm_bc == 'bc2' ) then
           call esm_force_ewr_pbc ( alpha_g, forceion )
        else if ( esm_bc == 'bc3' ) then
           call esm_force_ewr_pbc ( alpha_g, forceion )
        else if ( esm_bc == 'bc4' ) then
           call esm_force_ewr_bc4 ( alpha_g, forceion )
        end if

     END SUBROUTINE esm_force_ewr

     SUBROUTINE esm_force_ewg( alpha_g, forceion )
        USE kinds,     ONLY : DP
        USE ions_base, ONLY : nat
        implicit none
        real(DP), intent(in)    :: alpha_g
        real(DP), intent(out)   :: forceion(3,nat) 

        if( esm_bc == 'pbc' ) then
           call esm_force_ewg_pbc ( alpha_g, forceion )
        else if ( esm_bc == 'bc1' ) then
           call esm_force_ewg_bc1 ( alpha_g, forceion )
        else if ( esm_bc == 'bc2' ) then
           call esm_force_ewg_bc2 ( alpha_g, forceion )
        else if ( esm_bc == 'bc3' ) then
           call esm_force_ewg_bc3 ( alpha_g, forceion )
        else if ( esm_bc == 'bc4' ) then
           call esm_force_ewg_bc4 ( alpha_g, forceion )
        end if

     END SUBROUTINE esm_force_ewg

     SUBROUTINE esm_force_lc( aux, forcelc )
        USE kinds,     ONLY : DP
        USE ions_base, ONLY : nat
        USE fft_base,  ONLY : dfftp
        implicit none
        complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
        real(DP)   , intent(inout) :: forcelc(3,nat)

        if( esm_bc == 'pbc' ) then
           call esm_force_lc_pbc ( aux, forcelc )
        else if ( esm_bc == 'bc1' ) then
           call esm_force_lc_bc1 ( aux, forcelc )
        else if ( esm_bc == 'bc2' ) then
           call esm_force_lc_bc2 ( aux, forcelc )
        else if ( esm_bc == 'bc3' ) then
           call esm_force_lc_bc3 ( aux, forcelc )
        else if ( esm_bc == 'bc4' ) then
           call esm_force_lc_bc4 ( aux, forcelc )
        end if

     END SUBROUTINE esm_force_lc

SUBROUTINE esm_rgen_2d ( dtau, rmax, mxr, at, bg, r, r2, nrm)
  !-----------------------------------------------------------------------
  !
  !   generates neighbours shells (cartesian, in units of lattice parameter)
  !   with length < rmax,and returns them in order of increasing length:
  !      r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
  !   where a1, a2, a3 are primitive lattice vectors. Other input variables:
  !     mxr = maximum number of vectors
  !     at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
  !     bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
  !   Other output variables:
  !     nrm = the number of vectors with r^2 < rmax^2
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(in)   :: mxr
  INTEGER, INTENT(out)  :: nrm
  REAL(DP), INTENT(in)  :: at(3,3), bg(3,3), dtau(3), rmax
  REAL(DP), INTENT(out) :: r(3,mxr), r2(mxr)
  !
  !    and here the local variables
  !
  INTEGER, ALLOCATABLE :: irr (:)
  INTEGER  ::  nm1, nm2, i, j, ipol, ir, indsw, iswap
  real(DP) :: ds(3), dtau0(3)
  real(DP) :: t (3), tt, swap
  real(DP), EXTERNAL :: dnrm2
  !
  !
  nrm = 0
  IF (rmax==0.d0) RETURN

  ! bring dtau into the unit cell centered on the origin - prevents trouble
  ! if atomic positions are not centered around the origin but displaced
  ! far away (remember that translational invariance allows this!)
  !
  ds(:) = matmul( dtau(:), bg(:,:) )
  ds(:) = ds(:) - anint(ds(:))
  dtau0(:) = matmul( at(:,:), ds(:) )
  !
  ALLOCATE (irr( mxr))
  !
  ! these are estimates of the maximum values of needed integer indices
  !
  nm1 = int (dnrm2 (3, bg (1, 1), 1) * rmax) + 2
  nm2 = int (dnrm2 (3, bg (1, 2), 1) * rmax) + 2
  !
  DO i = -nm1, nm1
     DO j = -nm2, nm2
        tt = 0.d0
        DO ipol = 1, 3
          t (ipol) = i*at (ipol, 1) + j*at (ipol, 2) - dtau0(ipol)
           tt = tt + t (ipol) * t (ipol)
        ENDDO
        IF (tt<=rmax**2.and.abs (tt) >1.d-10) THEN
           nrm = nrm + 1
           IF (nrm>mxr) CALL errore ('esm_rgen_2d', 'too many r-vectors', nrm)
           DO ipol = 1, 3
              r (ipol, nrm) = t (ipol)
           ENDDO
           r2 (nrm) = tt
        ENDIF
     ENDDO
  ENDDO
  !
  !   reorder the vectors in order of increasing magnitude
  !
  !   initialize the index inside sorting routine
  !
  irr (1) = 0
  IF (nrm>1) CALL hpsort (nrm, r2, irr)
  DO ir = 1, nrm - 1
20   indsw = irr (ir)
     IF (indsw/=ir) THEN
        DO ipol = 1, 3
           swap = r (ipol, indsw)
           r (ipol, indsw) = r (ipol, irr (indsw) )
           r (ipol, irr (indsw) ) = swap
        ENDDO
        iswap = irr (ir)
        irr (ir) = irr (indsw)
        irr (indsw) = iswap
        GOTO 20
     ENDIF
  ENDDO
  DEALLOCATE(irr)
  !
  RETURN
END SUBROUTINE esm_rgen_2d

SUBROUTINE esm_init()
   USE fft_base, ONLY : dfftp
   USE esm_cft,  ONLY : esm_cft_1z_init
   IMPLICIT NONE
   
   call esm_ggen_2d()
   call esm_cft_1z_init(1,dfftp%nr3,dfftp%nr3)

END SUBROUTINE esm_init

SUBROUTINE esm_ggen_2d()
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, mill
  !
  IMPLICIT NONE
  !
  integer              :: n1xh, n2xh, ng, n1, n2, ng_2d
  Logical,allocatable  :: do_mill_2d(:,:)
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
!*** ngm_2d = total number of vectors (h,k) on this proc, excluding 
!*** duplicates with different l values
  
  IF( .not. ALLOCATED(mill_2d ) ) ALLOCATE( mill_2d(2,ngm_2d) )
  IF( .not. ALLOCATED(imill_2d) ) ALLOCATE( imill_2d(-n1xh:n1xh,-n2xh:n2xh) )
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
SUBROUTINE esm_hartree_pbc (rhog, ehart, aux)
  USE gvect,    ONLY : ngm
  USE lsda_mod, ONLY : nspin
  USE fft_base, ONLY : dfftp
  IMPLICIT NONE
  real(DP)    :: ehart             !  Hartree energy
  complex(DP) :: rhog(ngm,nspin)   !  n(G)
  complex(DP) :: aux(dfftp%nnr)    !  v_h(G)
  
  stop 'esm_hartree must not be called for esm_bc = pbc'

END SUBROUTINE esm_hartree_pbc

SUBROUTINE esm_hartree_bc1(rhog, ehart, aux)

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : nl, nlm, ngm, mill
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE mp_global,        ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm,nspin)   !  n(G)      
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)   
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp1, tmp2, tmp3
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( ng, n1, n2, ng_2d, n3, rg3 )
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
!$omp end parallel do
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel private( k1, k2, t, gp2, gp, tmp1, tmp2, vg, iz, kn, cc0, ss0, &
!$omp                   rg3, vg_r, k3, z, arg1, arg2, ng_2d )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc1
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp-ci*kn)
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp+ci*kn)
     enddo
     vg3(:,ng_2d)=vg(:)

     ! real part
     vg_r(:) = (0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc1
        arg1= gp*(z-z0)
        arg2=-gp*(z+z0)
        vg_r(iz)=-tpi/gp*(exp(arg1)*tmp1+exp(arg2)*tmp2)
     enddo
     
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)
!$omp end parallel

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)=-tpi*z0**2*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*ci*(cc0+ci*ss0)/kn
        tmp2=tmp2+rg3*ci*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*cc0/kn**2
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0)
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc1
        vg_r(iz)=-tpi*z**2*rg3 &
                 -tpi*(z-z0)*tmp1 &
                 -tpi*(z+z0)*tmp2 &
                 -fpi*tmp3
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*z_r**2*rg3 &
        -tpi*(z_r-z0)*tmp1 &
        -tpi*(z_r+z0)*tmp2 &
        -fpi*tmp3
     f2=-tpi*z_l**2*rg3 &
        -tpi*(z_l-z0)*tmp1 &
        -tpi*(z_l+z0)*tmp2 &
        -fpi*tmp3
     f3=-fpi*z_r*rg3 &
        -tpi*tmp1 &
        -tpi*tmp2
     f4=-fpi*z_l*rg3 &
        -tpi*tmp1 &
        -tpi*tmp2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing

     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, eh )
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
  ehart = ehart*omega*0.5d0
  !
  call mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE esm_hartree_bc1

SUBROUTINE esm_hartree_bc2 (rhog, ehart, aux)

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : nl, nlm, ngm, mill
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE mp_global,        ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm,nspin)   !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2, arg3, arg4, arg5
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( ng, n1, n2, ng_2d, n3, rg3 )
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
!$omp end parallel do
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel private( k1, k2, t, gp2, gp, tmp1, tmp2, vg, iz, kn, cc0, ss0, &
!$omp                   rg3, tmp, vg_r, k3, z, arg1, arg2, arg3, arg4, arg5, &
!$omp                   ng_2d )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc2
        arg1= gp*(z1-z0)
        arg2=-gp*(z1-z0)
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp=((gp+ci*kn)*exp(arg1)+(gp-ci*kn)*exp(arg2))/(2.d0*gp)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
        tmp=((gp-ci*kn)*exp(arg1)+(gp+ci*kn)*exp(arg2))/(2.d0*gp)
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:)=(0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc2
        arg1= gp*(z-z1)
        arg2=-gp*(z+z1)
        arg3= gp*(z-3.d0*z1)
        arg4=-gp*(z+3.d0*z1)
        arg5=-4.d0*gp*z1
        vg_r(iz)=-fpi*(exp(arg1)-exp(arg4))*tmp1/(1.d0-exp(arg5)) &
           +fpi*(exp(arg3)-exp(arg2))*tmp2/(1.d0-exp(arg5))
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)
!$omp end parallel

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)= tpi*(2.d0*z1-z0)*z0*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn**2
        tmp3=tmp3+rg3*ci*cc0/kn
        tmp4=tmp4+rg3*ss0/kn
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0)
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg_r(iz)=-tpi*z**2*rg3 &
                 -tpi*(z+z1)*tmp1/z1 &
                 +tpi*(z-z1)*tmp2/z1 &
                 -fpi*z*(z1-z0)/z1*tmp3 &
                 +fpi*(z1-z0)*tmp4
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*z_r**2*rg3 &
        -tpi*(z_r+z1)*tmp1/z1 &
        +tpi*(z_r-z1)*tmp2/z1 &
        -fpi*z_r*(z1-z0)/z1*tmp3 &
        +fpi*(z1-z0)*tmp4
     f2=-tpi*z_l**2*rg3 &
        -tpi*(z_l+z1)*tmp1/z1 &
        +tpi*(z_l-z1)*tmp2/z1 &
        -fpi*z_l*(z1-z0)/z1*tmp3 &
        +fpi*(z1-z0)*tmp4
     f3=-fpi*z_r*rg3 &
        -tpi*tmp1/z1 &
        +tpi*tmp2/z1 &
        -fpi*(z1-z0)/z1*tmp3
     f4=-fpi*z_l*rg3 &
        -tpi*tmp1/z1 &
        +tpi*tmp2/z1 &
        -fpi*(z1-z0)/z1*tmp3
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing
     
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, eh )
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
  ehart = ehart*omega*0.5d0
  !
  call mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE esm_hartree_bc2

SUBROUTINE esm_hartree_bc3 (rhog, ehart, aux)

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : nl, nlm, ngm, mill
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE mp_global,        ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm,nspin)   !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2, arg3
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( ng, n1, n2, ng_2d, n3, rg3 )
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
!$omp end parallel do
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel private( k1, k2, t, gp2, gp, tmp1, tmp2, vg, iz, kn, cc0, ss0, &
!$omp                   rg3, tmp, vg_r, k3, z, arg1, arg2, arg3, ng_2d )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc3
        arg1= gp*(z1-z0)
        arg2=-gp*(z1-z0)
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp=((gp+ci*kn)*exp(arg1)+(gp-ci*kn)*exp(arg2))/(2.d0*gp)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
        tmp=(gp-ci*kn)/gp
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
     enddo
     vg3(:,ng_2d)=vg(:)

     ! real part
     vg_r(:)=(0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc3
        arg1= gp*(z-z1)
        arg2=-gp*(z+z0)
        arg3= gp*(z-z0-2.d0*z1)
        vg_r(iz)=-fpi*exp(arg1)*tmp1+tpi*(exp(arg3)-exp(arg2))*tmp2
     enddo
     
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)
!$omp end parallel

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)= tpi*(4.d0*z1-z0)*z0*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*(cc0+ci*ss0)/kn
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0) 
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg_r(iz)=-tpi*(z**2+2.d0*z*z0)*rg3 &
                 -fpi*tmp1 &
                 -fpi*ci*(z-z1)*tmp2 &
                 -fpi*ci*(z1-z0)*tmp3
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*(z_r**2+2.d0*z_r*z0)*rg3 &
        -fpi*tmp1 &
        -fpi*ci*(z_r-z1)*tmp2 &
        -fpi*ci*(z1 -z0)*tmp3
     f2=-tpi*(z_l**2+2.d0*z_l*z0)*rg3 &
        -fpi*tmp1 &
        -fpi*ci*(z_l-z1)*tmp2 &
        -fpi*ci*(z1 -z0)*tmp3
     f3=-tpi*(2.d0*z_r+2.d0*z0)*rg3 &
        -fpi*ci*tmp2
     f4=-tpi*(2.d0*z_l+2.d0*z0)*rg3 &
        -fpi*ci*tmp2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing

     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, eh )
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
  ehart = ehart*omega*0.5d0
  !
  call mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE esm_hartree_bc3

SUBROUTINE esm_hartree_bc4 (rhog, ehart, aux)

  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : nl, nlm, ngm, mill
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE mp_global,        ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm,nspin)   !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, aaa, cc1, ss1, alpha, beta, &
                              chi, xi, kappa, lambda, arg1, arg2, arg3, &
                              arg4, argr1, argr2, argr3, argr4, argr5
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4, &
                              tmpr1, tmpr2, tmpr3, tmpr4
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:), vr(:), &
                              vr_r(:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( ng, n1, n2, ng_2d, n3, rg3 )
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
!$omp end parallel do
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  aaa=esm_a
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel private( k1, k2, t, gp2, gp, tmp1, tmp2, tmp3, tmp4, tmpr1, &
!$omp                   tmpr2, tmpr3, tmpr4, vg, vr, iz, kn, cc0, ss0, rg3, &
!$omp                   tmp, cc1, ss1, alpha, beta, kappa, xi, chi, lambda, &
!$omp                   vg_r, vr_r, k3, z, arg1, arg2, arg3, arg4, argr1, &
!$omp                   argr2, argr3, argr4, argr5, ng_2d )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3),vr(dfftp%nr3),vr_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0)
     tmp4=(0.d0,0.d0); tmpr1=(0.d0,0.d0); tmpr2=(0.d0,0.d0)
     tmpr3=(0.d0,0.d0); tmpr4=(0.d0,0.d0)
     vr(:)=(0.d0,0.d0); vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc4
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        vr(iz)=fpi*rg3/(gp**2+kn**2+ci*aaa*kn)
        cc1=cos(kn*z1)
        ss1=sin(kn*z1)
        alpha=aaa+gp+sqrt(aaa**2+gp**2)
        beta =aaa+gp-sqrt(aaa**2+gp**2)
        kappa=aaa-gp+sqrt(aaa**2+gp**2)
        xi   =aaa   +sqrt(aaa**2+gp**2)
        chi  =aaa   -sqrt(aaa**2+gp**2)
        lambda=      sqrt(aaa**2+gp**2)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(xi-ci*kn)/alpha
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp+ci*kn)/gp
        tmp3=tmp3+rg3*kappa/alpha*(cc0-ci*ss0)/(gp+ci*kn)/gp
        tmp4=tmp4+rg3*kappa*(cc1+ci*ss1)/(xi-ci*kn)/(gp**2+kn**2)
        tmpr1=tmpr1+rg3*(cc0-ci*ss0)/(gp+ci*kn)/alpha
        tmpr2=tmpr2+rg3*(cc0+ci*ss0)/(xi-ci*kn)/lambda
        tmpr3=tmpr3+rg3*beta/alpha*(cc0+ci*ss0)/(xi-ci*kn)/lambda
        tmpr4=tmpr4+rg3*beta*(cc1+ci*ss1)/(gp+ci*kn) &
           /(gp**2+kn**2+ci*2.d0*aaa*kn)
     enddo
     
     call esm_cft_1z(vg,1,dfftp%nr3,dfftp%nr3,1,vg_r)
     ! bc4
     CALL esm_cft_1z(vr,1,dfftp%nr3,dfftp%nr3,1,vr_r)
     
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc4
        arg1= gp*(z-z1)-xi*(z0-z1)
        arg2=-gp*(z+z0)
        arg3=-gp*(z0+z1)+gp*(z-z1)
        arg4= gp*(z-z1)
        argr1=-gp*(z0+z1)-xi*(z-z1)
        argr2=-xi*(z0-z1)-chi*(z-z1)
        argr3=-xi*(z-z1)-xi*(z0-z1)
        argr4=-xi*(z-z1)
        argr5=-2.d0*aaa*(z-z1)
        if (z < z1) then
           vg_r(iz) = vg_r(iz)-fpi*exp(arg1)*tmp1-tpi*exp(arg2)*tmp2 &
              +tpi*exp(arg3)*tmp3-fpi*exp(arg4)*tmp4
        else
           vg_r(iz) = vr_r(iz)*exp(argr5) &
              -fpi*exp(argr1)*tmpr1-tpi*exp(argr2)*tmpr2 &
              +tpi*exp(argr3)*tmpr3-fpi*exp(argr4)*tmpr4
        endif
     enddo
     
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=vg(:)*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r,vr,vr_r)
!$omp end parallel

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3),vr(dfftp%nr3),vr_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     vg(:)=(0.d0,0.d0); vr(:)=(0.d0,0.d0)
     !for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     !
     rg3=rhog3(1,ng_2d)
     ! bc4
     arg1=-2.d0*aaa*(z0-z1)
     vg(1)= tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
        -pi *(exp(arg1)-1.d0)/aaa**2*rg3
     vr(1)= tpi*(z0+0.5d0/aaa)/aaa*rg3

     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        ! bc4
        cc0= cos(kn*z0)
        ss0= sin(kn*z0)
        cc1=cos(kn*z1)
        ss1=sin(kn*z1)
        tmp1=tmp1+rg3*(cc1+ci*ss1)/(2.d0*aaa-ci*kn)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*(cc0+ci*ss0)/(2.d0*aaa-ci*kn)
        tmp4=tmp4+(0.d0,0.d0)

        vg(iz)=fpi*rg3/(kn**2)
        ! bc4
        vr(iz)=fpi*rg3/(kn**2+ci*2.d0*aaa*kn)

        !for smoothing
        c_r=cos(kn*z_r)
        s_r=sin(kn*z_r)
        c_l=cos(kn*z_l)
        s_l=sin(kn*z_l)
        ! bc4
        f1=f1+fpi*   rg3*(c_r+ci*s_r)/(kn**2+ci*2.d0*aaa*kn)
        f2=f2+fpi*   rg3*(c_l+ci*s_l)/kn**2
        f3=f3+fpi*ci*rg3*(c_r+ci*s_r)/kn
        f4=f4+fpi*ci*rg3*(c_l+ci*s_l)/kn
        !
     enddo
     
     call esm_cft_1z(vg,1,dfftp%nr3,dfftp%nr3,1,vg_r)
     ! bc4
     call esm_cft_1z(vr,1,dfftp%nr3,dfftp%nr3,1,vr_r)
     
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc4
        arg1=-2.d0*aaa*(z0-z1)
        arg2=-2.d0*aaa*(z-z1)
        if (z < z1) then
           vg_r(iz)=vg_r(iz) &
              -fpi*2.d0*aaa*tmp1 &
              -tpi*ci*(2.d0*(z -z1)-1.d0/aaa)*tmp2 &
              -tpi*exp(arg1)/aaa*tmp3 &
              -tpi*z*(z+2.d0*z0)*rg3
        else
           vg_r(iz)=vr_r(iz)*exp(arg2) &
              +tpi*ci*exp(arg2)/aaa*tmp2 &
              -tpi*exp(arg1)/aaa*tmp3 &
              +tpi*exp(arg2)*z/aaa*rg3 &
                    -pi *exp(arg1)/aaa**2*rg3
        endif
     enddo

     !for smoothing
     ! bc4
     arg1=-2.d0*aaa*(z0-z1)
     arg2=-2.d0*aaa*(z_r-z1)
     f1=f1+tpi*(z0+0.5d0/aaa)/aaa*rg3
     f1=f1*exp(arg2) &
        +tpi*ci*exp(arg2)/aaa*tmp2 &
        -tpi*exp(arg1)/aaa*tmp3 &
        +tpi*exp(arg2)*z_r/aaa*rg3 &
        -pi *exp(arg1)/aaa**2*rg3
     f2=f2+tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
        -pi *(exp(arg1)-1.d0)/aaa**2*rg3
     f2=f2 &
        -fpi*2.d0*aaa*tmp1 &
        -tpi*ci*(2.d0*(z_l-z1)-1.d0/aaa)*tmp2 &
        -tpi*exp(arg1)/aaa*tmp3 &
        -tpi*z_l*(z_l+2.d0*z0)*rg3
     f3=f3*exp(arg2) &
        -fpi*ci*exp(arg2)*tmp2 &
        -fpi*(z_r+z0)*exp(arg2)*rg3 
     f4=f4-fpi*ci*tmp2-fpi*(z_l+z0)*rg3   
     ! for smoothing
     !factor e2 will be multiplied later (at vg3 <= vg)
     !f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     
     vg3(:,ng_2d)=vg(:)*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r,vr,vr_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, eh )
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
  ehart = ehart*omega*0.5d0
  !
  call mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE esm_hartree_bc4


FUNCTION esm_ewald()
  !-----------------------------------------------------------------------
  !
  ! Calculates Ewald energy with both G- and R-space terms.
  ! Determines optimal alpha. Should hopefully work for any structure.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi, e2
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE cell_base, ONLY : tpiba2
  USE ions_base, ONLY : zv, nat, ityp
  USE gvect,     ONLY : gcutm
  implicit none

  real(DP) :: esm_ewald
  ! output: the ewald energy
  !
  !    here the local variables
  !
  integer :: na
  ! counter on atoms

  real(DP) :: charge, ewaldg, ewaldr, alpha, upperbound
  ! total ionic charge in the cell
  ! ewald energy computed in reciprocal space
  ! ewald energy computed in real space
  ! alpha term in ewald sum
  ! the maximum radius to consider real space sum

  charge = 0.d0
  do na = 1, nat
     charge = charge+zv( ityp(na) )
  enddo

  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  alpha = 2.9d0
  do
     alpha = alpha - 0.1d0
     if (alpha.le.0.d0) call errore ('esm_ewald', 'optimal alpha not found', 1)
     upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * &
        qe_erfc ( sqrt (tpiba2 * gcutm / 4.d0 / alpha) )
     if ( upperbound < 1.0d-7 ) exit
  end do

  ! G-space sum here.
  ! Determine if this processor contains G=0 and set the constant term
  CALL esm_ewaldg( alpha, ewaldg )

  ! R-space sum here (only for the processor that contains G=0)
  CALL esm_ewaldr( alpha, ewaldr )

  esm_ewald = 0.5d0 * e2 * ( ewaldg + ewaldr )

  call mp_sum(  esm_ewald, intra_bgrp_comm )
  !write( *,'(5x,"alpha used in ewald term: ",f5.2 )')alpha

  return
END FUNCTION esm_ewald


!-----------------------------------------------------------------------
!--------------ESM EWALD RSUM SUBROUTINE--------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_ewaldr_pbc ( alpha_g, ewr )

  USE io_global,        ONLY : stdout
  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE control_flags,    ONLY : iverbosity
  USE mp,               ONLY : mp_rank, mp_size
  USE mp_global,        ONLY : intra_bgrp_comm

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewr
  !
  !    here the local variables
  !
  integer            :: na, nb, nr, nrm, np, ip, ith
  ! counter on atoms
  ! counter on atoms
  ! counter over direct vectors
  ! number of R vectors included in r sum
  integer, parameter :: mxr = 500
  ! the maximum number of R vectors included in r
  real(DP)           :: dtau(3), r(3,mxr), r2(mxr)
  ! the difference tau_s - tau_s'
  ! neighbering shell vector
  ! the square modulus of R_j-tau_s-tau_s'
  ! buffer variable
  ! buffer variable
  !
  ! ESM variables
  !
  real(DP)            :: tmp, fac, ss, ew, rmax0, rr
  !
#if defined __OPENMP
  integer, external   :: OMP_GET_THREAD_NUM
#endif
  !
  ewr = 0.d0

  tmp=sqrt(alpha_g)
  rmax0=4.d0/tmp/alat

  ip = mp_rank( intra_bgrp_comm )
  np = mp_size( intra_bgrp_comm )

!$omp parallel private( na, nb, dtau, fac, r, r2, nrm, nr, &
!$omp                   ss, ew, ith )
  ith = 0
#if defined __OPENMP
  ith = OMP_GET_THREAD_NUM()
#endif

  ew = 0.d0
  do na = ip+1, nat, np
!$omp do
     do nb = 1, nat
        dtau(:)=tau(:,na)-tau(:,nb)
        fac=zv(ityp(nb))*zv(ityp(na))
        !
        ! generates nearest-neighbors shells
        !
        call rgen(dtau,rmax0,mxr,at,bg,r,r2,nrm)
        !
        ! and sum to the real space part
        !
        do nr=1,nrm
           rr=sqrt(r2(nr))*alat
           ew=ew+fac*qe_erfc(tmp*rr)/rr
        enddo
     enddo
!$omp end do
     if( ith /= 0 ) cycle
     ! Here add the other constant term
     ew=ew-zv(ityp(na))**2*tmp/sqrt(pi)*2.d0 ! 2.d0: fit to original code
  enddo
!$omp atomic
  ewr = ewr + ew
!$omp end parallel

END SUBROUTINE esm_ewaldr_pbc

SUBROUTINE esm_ewaldr_bc4 ( alpha_g, ewr )

  USE io_global,        ONLY : stdout
  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE control_flags,    ONLY : iverbosity
  USE mp,               ONLY : mp_rank, mp_size
  USE mp_global,        ONLY : intra_bgrp_comm

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewr
  !
  !    here the local variables
  !
  integer            :: na, nb, nr, nrm, np, ip, ith
  ! counter on atoms
  ! counter on atoms
  ! counter over direct vectors
  ! number of R vectors included in r sum
  integer, parameter :: mxr = 500
  ! the maximum number of R vectors included in r
  real(DP)           :: dtau(3), r(3,mxr), r2(mxr), rxy, rxyz
  ! the difference tau_s - tau_s'
  ! neighbering shell vector
  ! the square modulus of R_j-tau_s-tau_s'
  ! buffer variable
  ! buffer variable
  !
  ! ESM variables
  !
  real(DP)            :: L, z, zp, z0, z1, aaa, tmp, fac, ss, ew, err, ss0, &
                         gpmax, rmax0, rmax, zbuff, znrm, rr
  ! gpmax: upper bound of g_parallel integral
  ! rmax: the maximum radius to consider real space sum
  ! zbuff: smearing width to avoid the singularity of the Force
  ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
  real(DP), parameter :: eps=1.d-11, epsneib=1.d-6
  !
#if defined __OPENMP
  integer, external   :: OMP_GET_THREAD_NUM
#endif
  !
  ewr = 0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  aaa=esm_a
  tmp=sqrt(alpha_g)
  zbuff=1.d0
  !
  ! Define upperbound for g_parallel integral
  err=1.d0; ss0=0.d0; gpmax=1.d0
  do
     gpmax=gpmax+1.d0
     if (gpmax.gt.1000.d0) &
        call errore ('esm_ewaldr', 'optimal gpmax not found', 1)
     call qromb(vl11,aaa,tmp,z1,z1-zbuff,z1-zbuff,0.0_DP,gpmax,ss)
     err=abs(ss-ss0); ss0=ss
     if ( err.lt.eps ) exit
  enddo
  ! Define znrm using the deviation from the constant term in RSUM
  znrm=z1
  do
     znrm=znrm-0.01d0
     if ( znrm.le.-z0 ) &
        call errore ('esm_ewaldr', 'optimal znrm not found', 1)
     call qromb(vl11,aaa,tmp,z1,znrm,znrm,0.0_DP,gpmax,ss)
     err=-2.d0*tmp/sqrt(pi)-ss*2.d0
     if(abs(err).lt.eps) exit
  enddo
  ! Define rmax for real space sum
  rmax=1.d0
  do
     rmax=rmax+1.d0
     if ( rmax.gt.200.d0 ) &
        call errore ('esm_ewaldr', 'optimal rmax not found', 1)
     call qromb(vl11j0,aaa,tmp,z1,z1-zbuff,z1-zbuff,rmax,gpmax,ss)
     err=1.d0/rmax+ss*2.d0
     if(abs(err).lt.epsneib) exit
  enddo
  rmax=rmax/alat
  if (iverbosity > 0) then
     write( stdout, '(5x,"=== Smooth-ESM RSUM parameters (Energy) ===")')
     write( stdout, '(5x,A,F10.2,A)') &
        'Upper bound of g_parallel integral:      ',gpmax,' (1/a.u.)'
     write( stdout, '(5x,A,F10.2,A)') &
        'Boundary for normal RSUM|Smooth-ESM RSUM:',z1-znrm,' (a.u.)'
     write( stdout, '(5x,A,F10.2,A)') &
        'Upper bound of real-space summation:     ',rmax*alat,' (a.u.)'
     write( stdout, '(5x,"===========================================")')
  endif

  ip = mp_rank( intra_bgrp_comm )
  np = mp_size( intra_bgrp_comm )

!$omp parallel private( na, z, nb, zp, dtau, fac, r, r2, nrm, nr, rxy, &
!$omp                   rxyz, ss, ew, ith )
  ith = 0
#if defined __OPENMP
  ith = OMP_GET_THREAD_NUM()
#endif

  ew = 0.d0
  do na = ip+1, nat, np
     z=tau(3,na)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
!$omp do
     do nb = 1, nat
        zp=tau(3,nb)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        dtau(1:2)=tau(1:2,na)-tau(1:2,nb)
        dtau(3)=(z-zp)/alat
        fac=zv(ityp(nb))*zv(ityp(na))
        if ( z < znrm ) then
           if ( zp < znrm ) then ! z in I, zp in I (normal RSUM)
              rmax0=4.d0/tmp/alat
              !
              ! generates nearest-neighbors shells
              !
              call rgen(dtau,rmax0,mxr,at,bg,r,r2,nrm)
              !
              ! and sum to the real space part
              !
              do nr=1,nrm
                 rr=sqrt(r2(nr))*alat
                 ew=ew+fac*qe_erfc(tmp*rr)/rr
              enddo
           elseif ( zp < z1 ) then ! z in I, zp in I
              call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl11j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*(1.d0/rxyz+ss*2.d0)
              enddo
           else ! z in I, zp in II
              call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl12j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*ss*2.d0
              enddo
           endif ! if for zp
        elseif ( z < z1 ) then ! znrm < z < z1
           call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
           if ( zp < z1 ) then ! z in I, zp in I
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl11j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*(1.d0/rxyz+ss*2.d0)
              enddo
           else ! z in I, zp in II
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl12j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*ss*2.d0
              enddo
           endif ! if for zp
        else ! z1 < z
           call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
           if ( zp < z1 ) then ! z in II, zp in I
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl21j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*ss*2.d0
              enddo
           else ! z in II, zp in II
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl22j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 ew=ew+fac*(exp(-aaa*(rxyz+z+zp-2.d0*z1))/rxyz+ss*2.d0)
              enddo
           endif
        endif ! if for z
     enddo
!$omp end do
     if( ith /= 0 ) cycle
     if (z < znrm ) then
        ss=-tmp/sqrt(pi)
     elseif (z < z1) then
        call qromb(vl11,aaa,tmp,z1,z,z,0.0_DP,gpmax,ss)
     else
        call qromb(vl22,aaa,tmp,z1,z,z,0.0_DP,gpmax,ss)
     endif
     ew=ew+zv(ityp(na))**2*ss*2.d0 ! 2.0: fit to original code
  enddo
!$omp atomic
  ewr = ewr + ew
!$omp end parallel

END SUBROUTINE esm_ewaldr_bc4

!-----------------------------------------------------------------------
!--------------ESM EWALD GSUM SUBROUTINE--------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_ewaldg_pbc ( alpha_g, ewg )

  USE constants,        ONLY : tpi
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, tpiba2
  USE ions_base,        ONLY : zv, nat, nsp, ityp
  USE control_flags,    ONLY : gamma_only
  USE gvect,            ONLY : ngm, gg
  USE vlocal,           ONLY : strf

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewg

  integer     :: ng
  real(DP)    :: charge, fact
  complex(DP) :: rhon

  charge = sum( zv( ityp(1:nat) ) )

  ! same of the GSUM part in ewald.f90
  if( gstart == 2 ) then
     ewg = - charge**2 / alpha_g / 4.0d0
  else
     ewg = 0.0d0
  endif
  if( gamma_only ) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do ng = gstart, ngm
     rhon = sum( zv(1:nsp) * CONJG( strf (ng, 1:nsp) ) )
     ewg = ewg + fact * abs (rhon) **2 * &
        exp( - gg(ng) * tpiba2 / alpha_g / 4.d0) / gg(ng) / tpiba2
  enddo
  ewg = 2.d0 * tpi / omega * ewg

END SUBROUTINE esm_ewaldg_pbc

SUBROUTINE esm_ewaldg_bc1 ( alpha_g, ewg )

  USE constants,        ONLY : pi, tpi, fpi
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
  USE control_flags,    ONLY : gamma_only

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewg
  !
  !    here the local variables
  !
  integer  :: k1, k2, it1, it2, ng_2d
  real(DP) :: gp2, t(2), gp, sa, z, zp, z0, L, t1, t2, tt, &
              tmp, cc1, cc2, kk1, kk2, ff, ew,arg001, arg002, &
              arg101, arg102
  
  ewg=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  tmp=sqrt(alpha_g)
  sa=omega/L
!$omp parallel private( ew, it1, it2, z, zp, tt, arg001, arg002, &
!$omp                   arg101, arg102, &
!$omp                   kk1, kk2, t1, t2, cc1, cc2, ng_2d, &
!$omp                   k1, k2, t, gp2, gp, ff )
  ew=0d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     tt=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc1
     arg001=-tmp**2*(z-zp)**2
     arg101= tmp*(z-zp)
     kk1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
     kk2=0.d0

     cc1=0.d0
     cc2=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if( k1==0 .and. k2==0 ) cycle
        t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
        gp2 = sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +  (k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc1
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg101=gp/2.d0/tmp-tmp*(z-zp)
        arg102=gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp
        cc2=0.d0
     enddo

     if( gamma_only ) then
        cc1=cc1*2d0
        cc2=cc2*2d0
     endif
     ew=ew+tt*(cc1+cc2)
     if(gstart==2) ew=ew+tt*(kk1+kk2)
  enddo
  enddo
!$omp end do
!$omp atomic
  ewg=ewg+ew
!$omp end parallel

  return
END SUBROUTINE esm_ewaldg_bc1

SUBROUTINE esm_ewaldg_bc2 ( alpha_g, ewg )

  USE constants,        ONLY : pi, tpi, fpi
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
  USE control_flags,    ONLY : gamma_only

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewg
  !
  !    here the local variables
  !
  integer  :: k1, k2, it1, it2, ng_2d
  real(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
              tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
              arg003, arg004, arg005, arg006, arg007, arg101,  &
              arg102
  
  ewg=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  tmp=sqrt(alpha_g)
  sa=omega/L
!$omp parallel private( ew, it1, it2, z, zp, tt, arg001, arg002, &
!$omp                   arg003, arg005, arg006, arg007, arg101, &
!$omp                   arg102, kk1, kk2, t1, t2, cc1, cc2, ng_2d, &
!$omp                   k1, k2, t, gp2, gp, ff )
  ew=0d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     tt=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc2
     arg001=-tmp**2*(z-zp)**2
     arg101= tmp*(z-zp)
     kk1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
     kk2=0.5d0*(z1-z*zp/z1)

     cc1=0.d0
     cc2=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if( k1==0 .and. k2==0 ) cycle
        t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
        gp2 = sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +  (k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc2
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg003=-gp*(z+zp+2.d0*z1)
        arg004= gp*(z+zp-2.d0*z1)
        arg005=-gp*(z-zp+4.d0*z1)
        arg006= gp*(z-zp-4.d0*z1)
        arg007=-4.d0*gp*z1
        arg101=gp/2.d0/tmp-tmp*(z-zp)
        arg102=gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp
        cc2=cc2+cos(ff)*(exp(arg006)+exp(arg005) &
           -exp(arg004)-exp(arg003) ) &
           /(1.d0-exp(arg007))/2.d0/gp
     enddo

     if( gamma_only ) then
        cc1=cc1*2d0
        cc2=cc2*2d0
     endif
     ew=ew+tt*(cc1+cc2)
     if(gstart==2) ew=ew+tt*(kk1+kk2)
  enddo
  enddo
!$omp end do
!$omp atomic
  ewg=ewg+ew
!$omp end parallel

  return
END SUBROUTINE esm_ewaldg_bc2

SUBROUTINE esm_ewaldg_bc3 ( alpha_g, ewg )

  USE constants,        ONLY : pi, tpi, fpi
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
  USE control_flags,    ONLY : gamma_only

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewg
  !
  !    here the local variables
  !
  integer  :: k1, k2, it1, it2, ng_2d
  real(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
              tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
              arg003, arg101, arg102
  
  ewg=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  tmp=sqrt(alpha_g)
  sa=omega/L
!$omp parallel private( ew, it1, it2, z, zp, tt, arg001, arg002, &
!$omp                   arg003, arg101, arg102, kk1, kk2, t1, t2, &
!$omp                   cc1, cc2, ng_2d, k1, k2, t, gp2, gp, ff )
  ew=0d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     tt=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc3
     arg001=-tmp**2*(z-zp)**2
     arg101= tmp*(z-zp)
     kk1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
     kk2=0.5d0*(2.d0*z1-z-zp)

     cc1=0.d0
     cc2=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if( k1==0 .and. k2==0 ) cycle
        t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
        gp2 = sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +  (k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc3
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg003= gp*(z+zp-2.d0*z1)
        arg101=gp/2.d0/tmp-tmp*(z-zp)
        arg102=gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp
        cc2=cc2+cos(ff)*(-exp(arg003))/2.d0/gp
     enddo

     if( gamma_only ) then
        cc1=cc1*2d0
        cc2=cc2*2d0
     endif
     ew=ew+tt*(cc1+cc2)
     if(gstart==2) ew=ew+tt*(kk1+kk2)
  enddo
  enddo
!$omp end do
!$omp atomic
  ewg=ewg+ew
!$omp end parallel

  return
END SUBROUTINE esm_ewaldg_bc3

SUBROUTINE esm_ewaldg_bc4 ( alpha_g, ewg )

  USE constants,        ONLY : pi, tpi, fpi
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
  USE control_flags,    ONLY : gamma_only

  implicit none
  real(DP), intent(in)  :: alpha_g
  real(DP), intent(out) :: ewg
  !
  !    here the local variables
  !
  integer  :: k1, k2, it1, it2, ng_2d
  real(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
              tmp, cc1, cc2, kk1, kk2, ff, ew,arg001, arg002, &
              arg003, arg005, arg006, arg007, arg008, &
              arg009, arg011, arg101, arg102, arg103, arg104, &
              arg106, arg107, arg109, arg111, arg113, aaa, t3, &
              alpha, beta, kappa, lambda, xi, chi
  
  ewg=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  aaa=esm_a
  tmp=sqrt(alpha_g)
  sa=omega/L
!$omp parallel private( ew, it1, it2, z, zp, tt, arg001, arg002, &
!$omp                   arg003, arg005, arg006, arg007, arg008, arg009, &
!$omp                   arg011, arg101, arg102, arg103, arg104, arg107, &
!$omp                   arg109, arg111, arg113, alpha, beta, kappa, xi, &
!$omp                   chi, lambda, kk1, kk2, t1, t2, t3, cc1, cc2, ng_2d, &
!$omp                   k1, k2, t, gp2, gp, ff )
  ew=0d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     tt=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc4
     arg001=-tmp**2*(z-zp)**2
     arg002=-tmp**2*(z1-zp)**2
     arg005=-2.d0*aaa*(z-z1)
     arg006= aaa**2/tmp**2+2.d0*aaa*(z1-zp)
     arg101= tmp*(z-zp)
     arg102= tmp*(z1-zp)
     arg104= aaa/tmp+tmp*(z-zp)
     arg106= aaa/tmp+tmp*(z1-zp)
     if (z < z1) then
        t1=-(z-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
        t2=0.5d0/aaa*exp_erfc(arg006,arg106)
        t3=0.5d0/aaa-(z-z1)+exp(arg002)/tmp/sqrt(pi) &
           -exp(arg001)/tmp/sqrt(pi)
        kk1=(t1+t2)/2.d0
        kk2=t3/2.d0
     else
        t1=-exp_erfc(arg005,arg101)/aaa 
        t2= exp_erfc(arg006,arg104)/aaa
        t3= exp(arg005)/aaa
        kk1=(t1+t2)/4.d0
        kk2=t3/2.d0
     endif

     cc1=0.d0
     cc2=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if( k1==0 .and. k2==0 ) cycle
        t(1:2) = k1*bg(1:2, 1)+k2*bg(1:2, 2)
        gp2 = sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff = ((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +  (k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc4
        alpha=aaa+gp+sqrt(aaa**2+gp**2)
        beta =aaa+gp-sqrt(aaa**2+gp**2)
        kappa=aaa-gp+sqrt(aaa**2+gp**2)
        xi   =aaa   +sqrt(aaa**2+gp**2)
        chi  =aaa   -sqrt(aaa**2+gp**2)
        lambda=      sqrt(aaa**2+gp**2)
        arg001= gp*(z-zp)
        arg002=-gp*(z-zp)
        arg003= gp*(z+zp-2.d0*z1)
        arg005=-gp*(z1-zp)-xi*(z-z1)
        arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
        arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
        arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
        arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
        arg101= gp/2.d0/tmp+tmp*(z-zp)
        arg102= gp/2.d0/tmp-tmp*(z-zp)
        arg103= gp/2.d0/tmp+tmp*(z1-zp)
        arg104= gp/2.d0/tmp-tmp*(z1-zp)
        arg107= xi/2.d0/tmp+tmp*(z-zp)
        arg109= xi/2.d0/tmp+tmp*(z1-zp)
        arg111=chi/2.d0/tmp+tmp*(z-zp)
        arg113=chi/2.d0/tmp+tmp*(z1-zp)
        if (z < z1) then
           t1= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
           t2= exp_erfc(arg002,arg102) &
              -kappa/alpha*exp_erfc(arg003,arg104)
           t3= exp_erfc(arg006,arg109)/alpha
           cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp
           cc2=cc2+cos(ff)*t3/2.d0
        else
           t1= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
           t2= exp_erfc(arg008,arg107) &
              -beta/alpha*exp_erfc(arg009,arg109)
           t3= exp_erfc(arg005,arg104)/alpha
           cc1=cc1+cos(ff)*(t1+t2)/4.d0/lambda
           cc2=cc2+cos(ff)*t3/2.d0
        endif
     enddo

     if( gamma_only ) then
        cc1=cc1*2d0
        cc2=cc2*2d0
     endif
     ew=ew+tt*(cc1+cc2)
     if(gstart==2) ew=ew+tt*(kk1+kk2)
  enddo
  enddo
!$omp end do
!$omp atomic
  ewg=ewg+ew
!$omp end parallel

  return
END SUBROUTINE esm_ewaldg_bc4



!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL SUBROUTINE---------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_local_pbc (aux)
  USE fft_base, ONLY : dfftp
  implicit none
  complex(DP), intent(inout) :: aux( dfftp%nnr )

  stop 'esm_local must not be called for esm_bc = pbc'

END SUBROUTINE esm_local_pbc

SUBROUTINE esm_local_bc1 (aux)

  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE control_flags,    ONLY : gamma_only
  USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  implicit none
  ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  complex(DP), intent(inout) :: aux( dfftp%nnr )
  !
  !    here the local variables
  !
  real(DP)                :: t(2), tt, gp, gp2, sa, z0, pp, cc, ss, &
                             t1, t2, z, zp, tmp, L, z_l, z_r, &
                             arg001, arg002, arg101, arg102
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                             nz_r, ng_2d
  complex(DP)             :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
  complex(DP),allocatable :: vloc3(:,:), vg(:), vg_r(:)

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0   ! Gaussian width
  allocate(vloc3(dfftp%nr3,ngm_2d))

! for gp!=0
!$omp parallel private( ng_2d, k1, k2, t, gp2, gp, vg, vg_r, it, tt, pp, &
!$omp                   cc, ss, cs, zp, iz, k3, z, cc1, cc2, t1, t2, &
!$omp                   arg001, arg002, arg101, arg102 )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2)=k1*bg(1:2, 1)+k2*bg (1:2, 2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp =sqrt(gp2)
     vg_r(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX (cc, ss,kind=DP)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc1
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           cc1=cs*(t1+t2)/4.d0/gp
           cc2=(0.d0,0.d0)
           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
  enddo
  deallocate(vg,vg_r)
!$omp end parallel
  
  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     vg_r(:)=(0.d0,0.d0)
! for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
! for gp=0
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc1
           arg001=-tmp**2*(z-zp)**2
           arg101= tmp*(z-zp)
           cc1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
           cc2=(0.d0,0.d0)

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc1
        f1=f1+tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
           -exp(-tmp**2*(z_r-zp)**2)/tmp/sqrt(pi))
        f2=f2+tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
           -exp(-tmp**2*(z_l-zp)**2)/tmp/sqrt(pi))
        f3=f3-tt*0.5d0*qe_erf(tmp*(z_r-zp))
        f4=f4-tt*0.5d0*qe_erf(tmp*(z_l-zp))
     enddo
     ! for smoothing
     ! factor e2: hartree -> Ry.
     f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
     
     deallocate(vg,vg_r)
  endif ! if( ng_2d > 0 )
  
! Map to FFT mesh (dfftp%nrx)
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
     if (gamma_only) then
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate(vloc3)

  return
END SUBROUTINE esm_local_bc1

SUBROUTINE esm_local_bc2 (aux)

  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE control_flags,    ONLY : gamma_only
  USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  implicit none
  ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  complex(DP), intent(inout) :: aux( dfftp%nnr )
  !
  !    here the local variables
  !
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                             t1, t2, z, zp, v0, tmp, L, z_l, z_r, &
                             arg001, arg002, arg003, arg004, arg005, &
                             arg006, arg007, arg101, arg102
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                             nz_r, ng_2d
  complex(DP)             :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
  complex(DP),allocatable :: vloc3(:,:), vg(:), vg_r(:)

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0   ! Gaussian width
  z1=z0+esm_w
  v0=esm_efield*z1*2.d0/e2 ! factor 1/e2: unit Ry. -> hartree
  allocate(vloc3(dfftp%nr3,ngm_2d))

! for gp!=0
!$omp parallel private( ng_2d, k1, k2, t, gp2, gp, vg, vg_r, it, tt, pp, &
!$omp                   cc, ss, cs, zp, iz, k3, z, cc1, cc2, t1, t2, &
!$omp                   arg001, arg002, arg003, arg005, arg006, arg007, &
!$omp                   arg101, arg102 )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2)=k1*bg(1:2, 1)+k2*bg (1:2, 2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp =sqrt(gp2)
     vg_r(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX (cc, ss,kind=DP)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc2
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003=-gp*(z+zp+2.d0*z1)
           arg004= gp*(z+zp-2.d0*z1)
           arg005=-gp*(z-zp+4.d0*z1)
           arg006= gp*(z-zp-4.d0*z1)
           arg007=-4.d0*gp*z1
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           cc1=cs*(t1+t2)/4.d0/gp
           cc2=cs*(exp(arg006)+exp(arg005)-exp(arg004)-exp(arg003)) &
              /(1.d0-exp(arg007))/2.d0/gp 

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
  enddo
  deallocate(vg,vg_r)
!$omp end parallel
  
  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     vg_r(:)=(0.d0,0.d0)
! for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
! add constant potential (capacitor term)
     ! bc2
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg_r(iz)=-0.5d0*v0*(z-z1)/z1*e2 ! factor e2: hartree -> Ry.
     enddo
     f1=-0.5d0*v0*(z_r-z1)/z1 ! unit: hartree
     f2=-0.5d0*v0*(z_l-z1)/z1 ! unit: hartree
     f3=-0.5d0*v0/z1 ! unit: hartree/a.u.
     f4=-0.5d0*v0/z1 ! unit: harteee/a.u.

! for gp=0
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc2
           arg001=-tmp**2*(z-zp)**2
           arg101= tmp*(z-zp)
           cc1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
           cc2=0.5d0*(z1-z*zp/z1)

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc2
        f1=f1+tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
           -exp(-tmp**2*(z_r-zp)**2)/tmp/sqrt(pi))
        f2=f2+tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
           -exp(-tmp**2*(z_l-zp)**2)/tmp/sqrt(pi))
        f3=f3-tt*0.5d0*qe_erf(tmp*(z_r-zp))
        f4=f4-tt*0.5d0*qe_erf(tmp*(z_l-zp))
        f1=f1+tt*0.5d0*(z1-z_r*zp/z1)
        f2=f2+tt*0.5d0*(z1-z_l*zp/z1)
        f3=f3+tt*(-0.5d0*(zp/z1))
        f4=f4+tt*(-0.5d0*(zp/z1))
     enddo
     ! for smoothing
     ! factor e2: hartree -> Ry.
     f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
     
     deallocate(vg,vg_r)
  endif ! if( ng_2d > 0 )
  
! Map to FFT mesh (dfftp%nrx)
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
     if (gamma_only) then
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate(vloc3)

  return
END SUBROUTINE esm_local_bc2

SUBROUTINE esm_local_bc3 (aux)

  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE control_flags,    ONLY : gamma_only
  USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  implicit none
  ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  complex(DP), intent(inout) :: aux( dfftp%nnr )
  !
  !    here the local variables
  !
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                             t1, t2, z, zp, tmp, L, z_l, z_r, &
                             arg001, arg002, arg003, arg101, arg102
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                             nz_r, ng_2d
  complex(DP)             :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
  complex(DP),allocatable :: vloc3(:,:), vg(:), vg_r(:)

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0   ! Gaussian width
  z1=z0+esm_w
  allocate(vloc3(dfftp%nr3,ngm_2d))

! for gp!=0
!$omp parallel private( ng_2d, k1, k2, t, gp2, gp, vg, vg_r, it, tt, pp, &
!$omp                   cc, ss, cs, zp, iz, k3, z, cc1, cc2, t1, t2, &
!$omp                   arg001, arg002, arg003, arg101, arg102 )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2)=k1*bg(1:2, 1)+k2*bg (1:2, 2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp =sqrt(gp2)
     vg_r(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX (cc, ss,kind=DP)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc3
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003= gp*(z+zp-2.d0*z1)
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           cc1=cs*(t1+t2)/4.d0/gp
           cc2=cs*(-exp(arg003))/2.d0/gp

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
  enddo
  deallocate(vg,vg_r)
!$omp end parallel
  
  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     vg_r(:)=(0.d0,0.d0)
! for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)

! for gp=0
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc3
           arg001=-tmp**2*(z-zp)**2
           arg101= tmp*(z-zp)
           cc1=0.5d0*(-(z-zp)*qe_erf(arg101)-exp(arg001)/tmp/sqrt(pi))
           cc2=0.5d0*(2.d0*z1-z-zp)

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc3
        f1=f1+tt*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
           -exp(-tmp**2*(z_r-zp)**2)/tmp/sqrt(pi))
        f2=f2+tt*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
           -exp(-tmp**2*(z_l-zp)**2)/tmp/sqrt(pi))
        f3=f3-tt*0.5d0*qe_erf(tmp*(z_r-zp))
        f4=f4-tt*0.5d0*qe_erf(tmp*(z_l-zp))
        f1=f1+tt*0.5d0*(2.d0*z1-z_r-zp)
        f2=f2+tt*0.5d0*(2.d0*z1-z_l-zp)
        f3=f3-tt*0.5d0
        f4=f4-tt*0.5d0
     enddo
     ! for smoothing
     ! factor e2: hartree -> Ry.
     f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo

     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
     
     deallocate(vg,vg_r)
  endif ! if( ng_2d > 0 )
  
! Map to FFT mesh (dfftp%nrx)
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
     if (gamma_only) then
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate(vloc3)

  return
END SUBROUTINE esm_local_bc3

SUBROUTINE esm_local_bc4 (aux)

  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE control_flags,    ONLY : gamma_only
  USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z
  !
  implicit none
  ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  complex(DP), intent(inout) :: aux( dfftp%nnr )
  !
  !    here the local variables
  !
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                             t1, t2, z, zp, tmp, L, z_l, z_r, &
                             arg001, arg002, arg003, arg005, &
                             arg006, arg008, arg009, arg011, &
                             arg101, arg102, arg103, arg104, arg106, &
                             arg107, arg109, arg111, arg113, aaa, t3, &
                             alpha, beta, kappa, lambda, xi, chi
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                             nz_r, ng_2d
  complex(DP)             :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
  complex(DP),allocatable :: vloc3(:,:), vg(:), vg_r(:)

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0   ! Gaussian width
  z1=z0+esm_w
  aaa=esm_a
  allocate(vloc3(dfftp%nr3,ngm_2d))

! for gp!=0
!$omp parallel private( ng_2d, k1, k2, t, gp2, gp, vg, vg_r, it, tt, pp, &
!$omp                   cc, ss, cs, zp, iz, k3, z, cc1, cc2, t1, t2, t3, &
!$omp                   arg001, arg002, arg003, arg005, arg006, arg008, &
!$omp                   arg009, arg011, arg101, arg102, arg103, arg104, &
!$omp                   arg107, arg109, arg111, arg113, alpha, beta, kappa, &
!$omp                   xi, chi, lambda )
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2)=k1*bg(1:2, 1)+k2*bg (1:2, 2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp =sqrt(gp2)
     vg_r(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX (cc, ss,kind=DP)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc4
           alpha=aaa+gp+sqrt(aaa**2+gp**2)
           beta =aaa+gp-sqrt(aaa**2+gp**2)
           kappa=aaa-gp+sqrt(aaa**2+gp**2)
           xi   =aaa   +sqrt(aaa**2+gp**2)
           chi  =aaa   -sqrt(aaa**2+gp**2)
           lambda=      sqrt(aaa**2+gp**2)
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003= gp*(z+zp-2.d0*z1)
           arg005=-gp*(z1-zp)-xi*(z-z1)
           arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
           arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
           arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
           arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
           arg101=  gp/2.d0/tmp+tmp*(z-zp)
           arg102=  gp/2.d0/tmp-tmp*(z-zp)
           arg103=  gp/2.d0/tmp+tmp*(z1-zp)
           arg104=  gp/2.d0/tmp-tmp*(z1-zp)
           arg107=  xi/2.d0/tmp+tmp*(z-zp)
           arg109=  xi/2.d0/tmp+tmp*(z1-zp)
           arg111= chi/2.d0/tmp+tmp*(z-zp)
           arg113= chi/2.d0/tmp+tmp*(z1-zp)
           if (z < z1) then
              t1= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3= exp_erfc(arg006,arg109)/alpha
              cc1=cs*(t1+t2)/4.d0/gp
              cc2=cs*t3/2.d0
           else
              t1= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
              t2= exp_erfc(arg008,arg107) &
                 -beta/alpha*exp_erfc(arg009,arg109)
              t3= exp_erfc(arg005,arg104)/alpha
              cc1=cs*(t1+t2)/4.d0/lambda
              cc2=cs*t3/2.d0
           endif

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     enddo
     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
  enddo
  deallocate(vg,vg_r)
!$omp end parallel
  
  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     vg_r(:)=(0.d0,0.d0)
! for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
! for gp=0
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc4
           arg001=-tmp**2*(z-zp)**2
           arg002=-tmp**2*(z1-zp)**2
           arg005=-2.d0*aaa*(z-z1)
           arg006= aaa**2/tmp**2+2.d0*aaa*(z1-zp)
           arg101= tmp*(z-zp)
           arg102= tmp*(z1-zp)
           arg104= aaa/tmp+tmp*(z-zp)
           arg106= aaa/tmp+tmp*(z1-zp)
           if (z < z1) then
              t1=-(z-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
              t2=0.5d0/aaa*exp_erfc(arg006,arg106)
              t3=0.5d0/aaa-(z-z1)+exp(arg002)/tmp/sqrt(pi) &
                 -exp(arg001)/tmp/sqrt(pi)
              cc1=(t1+t2)/2.d0
              cc2=t3/2.d0
           else
              t1=-exp_erfc(arg005,arg101)/aaa
              t2= exp_erfc(arg006,arg104)/aaa
              t3= exp(arg005)/aaa
              cc1=(t1+t2)/4.d0
              cc2=t3/2.d0
           endif

           vg_r(iz) = vg_r(iz)+tt*(cc1+cc2)*e2 ! factor e2: hartree -> Ry.
        enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc4
        arg002=-tmp**2*(z1-zp)**2
        arg006= aaa**2/tmp**2+2.d0*aaa*(z1-zp)
        arg102= tmp*(z1-zp)
        arg106= aaa/tmp+tmp*(z1-zp)
        !-right only
        arg005=-2.d0*aaa*(z_r-z1)
        arg101= tmp*(z_r-zp)
        arg104= aaa/tmp+tmp*(z_r-zp)
        !--
        t1=-exp_erfc(arg005,arg101)/aaa
        t2= exp_erfc(arg006,arg104)/aaa
        t3= exp(arg005)/aaa
        f1=f1+tt*((t1+t2)/2.d0+t3)/2.d0
        f3=f3-tt*0.5d0*exp(arg005)*(1.d0+qe_erf(arg101))
        !-left only
        arg001=-tmp**2*(z_l-zp)**2
        arg101= tmp*(z_l-zp)
        !--
        t1=-(z_l-zp)*qe_erf(arg101)+(0.5d0/aaa+z1-zp)*qe_erf(arg102)
        t2=0.5d0/aaa*exp_erfc(arg006,arg106)
        t3=0.5d0/aaa-(z_l-z1)+exp(arg002)/tmp/sqrt(pi) &
           -exp(arg001)/tmp/sqrt(pi)
        f2=f2+tt*(t1+t2+t3)/2.d0
        f4=f4-tt*0.5d0*(1.d0+qe_erf(arg101))
     enddo
     ! for smoothing
     ! factor e2: hartree -> Ry.
     f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo

     call esm_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     do iz=1,dfftp%nr3
        vloc3(iz,ng_2d)=vg(iz)
     enddo
     
     deallocate(vg,vg_r)
  endif ! if( ng_2d > 0 )
  
! Map to FFT mesh (dfftp%nrx)
!$omp parallel do private( ng, n1, n2, ng_2d, n3 )
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + dfftp%nr3
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
     if (gamma_only) then
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     endif
  enddo
!$omp end parallel do

  deallocate(vloc3)

  return
END SUBROUTINE esm_local_bc4


SUBROUTINE esm_force_ew( forceion )
  !-----------------------------------------------------------------------
  !
  !  This routine computes the Ewald contribution to the forces,
  !  both the real- and reciprocal-space terms are present
  !
  USE kinds
  USE constants, ONLY : tpi, e2
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE ions_base, ONLY : zv, nat, ityp
  USE gvect,     ONLY : gcutm
  USE cell_base, ONLY : tpiba2
  implicit none

  real(DP), intent(out) :: forceion( 3, nat )
  ! output: the ewald part of the forces
  !
  real(DP) :: alpha, charge, upperbound
  ! the alpha parameter
  ! the total charge
  ! used to determine alpha

  forceion(:,:) = 0.d0
  charge = sum( zv( ityp(1:nat) ) )
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error ON THE ENERGY
  !
  alpha = 2.9d0
  do
     alpha = alpha - 0.1d0
     if( alpha == 0.d0) then
        call errore ('esm_force_ew', 'optimal alpha not found', 1)
     end if
     upperbound = e2 * charge**2 * sqrt (2.d0 * alpha / tpi) * &
        qe_erfc ( sqrt (tpiba2 * gcutm / 4.d0 / alpha) )
     if( upperbound < 1.0d-7) exit
  end do
  !write(*,'(5X,A,F5.2)')'alpha used in esm ewald force :',alpha

  CALL esm_force_ewg ( alpha, forceion )

  CALL esm_force_ewr ( alpha, forceion )

  CALL mp_sum( forceion, intra_bgrp_comm )

  return
END SUBROUTINE esm_force_ew


!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE (RSUM) SUBROUTINE ---------------
!-----------------------------------------------------------------------
SUBROUTINE esm_force_ewr_pbc ( alpha_g, forceion ) 
  USE constants,        ONLY : pi, e2
  USE cell_base,        ONLY : alat, at, bg
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE mp,               ONLY : mp_rank, mp_size
  USE mp_global,        ONLY : intra_bgrp_comm

  implicit none
  integer               :: na, nb, nr, nrm, ip, np, ith
  ! counter on atoms
  ! counter on atoms
  ! counter over direct vectors
  ! number of R vectors included in r sum
  integer, parameter    :: mxr = 1000
  ! the maximum number of R vectors included in r
  real(DP)              :: dtau(3), r(3,mxr), r2(mxr)
  ! the difference tau_s - tau_s'
  ! neighbering shell vector
  ! the square modulus of R_j-tau_s-tau_s'
  real(DP),intent(in)   :: alpha_g
  real(DP),intent(inout):: forceion(3,nat) 
  !
  ! ESM variables
  !
  real(DP)              :: tmp, fac, rmax0, rr
  ! rmax0: the maximum radius to consider real space sum
  real(DP), allocatable :: force(:,:)
#if defined __OPENMP
  integer, external   :: OMP_GET_THREAD_NUM
#endif

  tmp=sqrt(alpha_g)
  rmax0 = 5.d0 / tmp / alat

  ip = mp_rank( intra_bgrp_comm )
  np = mp_size( intra_bgrp_comm )

!$omp parallel private( force, na, nb, dtau, fac, r, r2, nrm, ith )
  ith = 0
#if defined __OPENMP
  ith = OMP_GET_THREAD_NUM()
#endif

  allocate( force(3,nat) )
  force(:,:)=0.d0
  do na = ip+1, nat, np
!$omp do
     do nb = 1, nat
        if (nb.eq.na)cycle
        dtau(:)=tau(:,na)-tau(:,nb)
        fac=zv(ityp(na))*zv(ityp(nb))*e2
        !
        ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
        !
        call rgen (dtau, rmax0, mxr, at, bg, r, r2, nrm)
        !
        ! and sum to the real space part
        !
        do nr=1,nrm
           rr=sqrt(r2(nr))*alat
           force(:,na)=force(:,na) &
              -fac/rr**2*(qe_erfc(tmp*rr)/rr+2.d0*tmp/sqrt(pi) &
              *exp(-tmp**2*rr**2))*r(:,nr)*alat
        enddo
     enddo
!$omp end do
  enddo
!$omp critical
  forceion(:,:)=forceion(:,:)+force(:,:)
!$omp end critical
  deallocate( force )
!$omp end parallel

END SUBROUTINE esm_force_ewr_pbc

SUBROUTINE esm_force_ewr_bc4 ( alpha_g, forceion ) 
  USE io_global,        ONLY : stdout
  USE constants,        ONLY : pi, tpi, fpi, e2
  USE gvect,            ONLY : gstart
  USE cell_base,        ONLY : alat, tpiba2, at, bg
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE control_flags,    ONLY : iverbosity
  USE mp,               ONLY : mp_rank, mp_size
  USE mp_global,        ONLY : intra_bgrp_comm

  implicit none
  integer               :: na, nb, nr, nrm, ipol, ip, np, ith
  ! counter on atoms
  ! counter on atoms
  ! counter over direct vectors
  ! number of R vectors included in r sum
  integer, parameter    :: mxr = 1000
  ! the maximum number of R vectors included in r
  real(DP)              :: dtau(3), r(3,mxr), r2(mxr), rxy, rxyz
  ! the difference tau_s - tau_s'
  ! neighbering shell vector
  ! the square modulus of R_j-tau_s-tau_s'
  ! buffer variable
  ! buffer variable
  real(DP),intent(in)   :: alpha_g
  real(DP),intent(inout):: forceion(3,nat) 
  !
  ! ESM variables
  !
  real(DP)              :: L, z, zp, z0, z1, aaa, tmp, ss, fac, err, ss0, &
                           gpmax, rmax0, rmax, zbuff, znrm, rr
  ! gpmax: upper bound of g_parallel integral
  ! rmax: the maximum radius to consider real space sum
  ! zbuff: smearing width to avoid the singularity of the Force
  ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
  real(DP), parameter :: eps=1.d-11, epsneib=1.d-6
  real(DP), allocatable :: force(:,:)
#if defined __OPENMP
  integer, external   :: OMP_GET_THREAD_NUM
#endif

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  aaa=esm_a
  tmp=sqrt(alpha_g)
  zbuff=1.d0
  !
  ! Define upperbound for g_parallel integral
  err=1.d0; ss0=0.d0; gpmax=1.d0
  do
     gpmax=gpmax+1.d0
     if (gpmax.gt.1000.d0) &
        call errore ('esm_force_ewr', 'optimal gpmax not found', 1)
     call qromb(vl11,aaa,tmp,z1,z1-zbuff,z1-zbuff,0.0_DP,gpmax,ss)
     err=abs(ss-ss0); ss0=ss
     if(err.lt.eps) exit
  enddo
  ! Define znrm using the deviation from the constant term in RSUM
  znrm=z1
  do
     znrm=znrm-0.01d0
     if (znrm.le.-z0) &
        call errore ('esm_force_ewr', 'optimal znrm not found', 1)
     call qromb(vl11,aaa,tmp,z1,znrm,znrm,0.0_DP,gpmax,ss)
     err=-2.d0*tmp/sqrt(pi)-ss*2.d0
     if(abs(err).lt.eps) exit
  enddo
  ! Define rmax for real space sum
  rmax=1.d0
  do
     rmax=rmax+1.d0
     if (rmax.gt.200.d0) &
        call errore ('esm_force_ewr', 'optimal rmax not found', 1)
     call qromb(dvl11j0,aaa,tmp,z1,z1-zbuff,z1-zbuff,rmax,gpmax,ss)
     err=ss
     if(abs(err).lt.epsneib) exit
  enddo
  rmax=rmax/alat
  if (iverbosity > 0) then
     write( stdout, '(5x,"=== Smooth-ESM RSUM parameters (Force) ===")')
     write( stdout, '(5x,A,F10.2,A)') &
        'Upper bound of g_parallel integral:      ',gpmax,' (1/a.u.)'
     write( stdout, '(5x,A,F10.2,A)') &
        'Boundary for normal RSUM|Smooth-ESM RSUM:',z1-znrm,' (a.u.)'
     write( stdout, '(5x,A,F10.2,A)') &
        'Upper bound of real-space summation:     ',rmax*alat,' (a.u.)'
     write( stdout, '(5x,"==========================================")')
  endif
  !
  ip = mp_rank( intra_bgrp_comm )
  np = mp_size( intra_bgrp_comm )

!$omp parallel private( force, na, z, nb, zp, dtau, fac, r, r2, nrm, rxy, &
!$omp                   rxyz, ss, ith )
  ith = 0
#if defined __OPENMP
  ith = OMP_GET_THREAD_NUM()
#endif

  allocate( force(3,nat) )
  force(:,:)=0.d0
  do na = ip+1, nat, np
     z=tau(3,na)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
!$omp do
     do nb = 1, nat
        if (nb.eq.na)cycle
        zp=tau(3,nb)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        dtau(1:2)=tau(1:2,na)-tau(1:2,nb)
        dtau(3)=(z-zp)/alat
        fac=zv(ityp(na))*zv(ityp(nb))*e2
        if ( z < znrm ) then
           if ( zp < znrm ) then ! z in I, zp in I (normal RSUM)
              rmax0 = 5.d0 / tmp / alat
              !
              ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
              !
              call rgen (dtau, rmax0, mxr, at, bg, r, r2, nrm)
              !
              ! and sum to the real space part
              !
              do nr=1,nrm
                 rr=sqrt(r2(nr))*alat
                 do ipol=1,3
                    force(ipol,na)=force(ipol,na) &
                       -fac/rr**2*(qe_erfc(tmp*rr)/rr+2.d0*tmp/sqrt(pi) &
                       *exp(-tmp**2*rr**2))*r(ipol,nr)*alat
                 enddo
              enddo
           elseif ( zp < z1 ) then ! z in I, zp in I
              call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl11j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -fac*(1.d0/rxyz**3+1.d0/rxy*ss)*r(1:2,nr)*alat
                 call qromb(dvl11j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb)-fac*((z-zp)/rxyz**3+ss)
              enddo
           else ! z in I, zp in II
              call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl12j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -fac*ss/rxy*r(1:2,nr)*alat
                 call qromb(dvl12j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb)-fac*ss
              enddo
           endif ! if for zp
        elseif ( z < z1 ) then ! znrm < z < z1
           call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
           if (zp < z1 ) then ! z in I, zp in I
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl11j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -fac*(1.d0/rxyz**3+1.d0/rxy*ss)*r(1:2,nr)*alat
                 call qromb(dvl11j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb)-fac*((z-zp)/rxyz**3+ss)
              enddo
           else ! z in I, zp in II
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl12j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -fac*ss/rxy*r(1:2,nr)*alat
                 call qromb(dvl12j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb)-fac*ss
              enddo
           endif ! if for zp
        else ! z1 < z
           call esm_rgen_2d(dtau,rmax,mxr,at,bg,r,r2,nrm)
           if (zp < z1 ) then ! z in II, zp in I
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 call qromb(vl21j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -fac*ss/rxy*r(1:2,nr)*alat
                 call qromb(dvl21j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb)-fac*ss
              enddo
           else ! z in II, zp in II
              do nr = 1, nrm
                 rxy = sqrt(r2(nr))*alat
                 rxyz = sqrt(r2(nr)+dtau(3)**2)*alat
                 call qromb(vl22j1,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(1:2,nb)=force(1:2,nb) &
                    -(exp(-aaa*(rxyz+z+zp-2.d0*z1))*(aaa+1.d0/rxyz)/rxyz**2 &
                    +ss/rxy)*fac*r(1:2,nr)*alat
                 call qromb(dvl22j0,aaa,tmp,z1,z,zp,rxy,gpmax,ss)
                 force(3,nb)=force(3,nb) &
                    -(exp(-aaa*(rxyz+z+zp-2.d0*z1))*(aaa+1.d0/rxyz)/rxyz**2 &
                    *(z-zp)-aaa*exp(-aaa*(rxyz+z+zp-2.d0*z1))/rxyz+ss)*fac
              enddo
           endif ! if for zp
        endif
     enddo
!$omp end do
     if( ith /= 0 ) cycle
     if (z < znrm) then
        ss=0.d0
     elseif (z < z1) then
        call qromb(dvl11,aaa,tmp,z1,z,z,0.0_DP,gpmax,ss)
     else
        call qromb(dvl22,aaa,tmp,z1,z,z,0.0_DP,gpmax,ss)
     endif
     ! factor e2: hartree -> Ry.
     force(3,na)=force(3,na)-zv(ityp(na))**2*e2*ss
  enddo
!$omp critical
  forceion(:,:)=forceion(:,:)+force(:,:)
!$omp end critical
  deallocate( force )
!$omp end parallel

END SUBROUTINE esm_force_ewr_bc4

!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE (GSUM) SUBROUTINE ---------------
!-----------------------------------------------------------------------
SUBROUTINE esm_force_ewg_pbc ( alpha_g, forceion )

  USE constants,        ONLY : tpi, e2
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
  USE gvect,            ONLY : gstart, ngm, gg, g
  USE vlocal,           ONLY : strf

  implicit none
  real(DP), intent(in)     :: alpha_g
  real(DP), intent(out)    :: forceion(3,nat) 
  integer                  :: nt, ig, na, ipol
  real(DP)                 :: fact, arg, sumnb
  complex(DP), allocatable :: aux (:)

  forceion(:,:)=0.d0

  ! same of the GSUM part in force_ew.f90
  allocate(aux(ngm))
  aux(:) = (0.d0, 0.d0)
  
  do nt = 1, nsp
     do ig = gstart, ngm
        aux (ig) = aux (ig) + zv (nt) * CONJG(strf (ig, nt) )
     enddo
  enddo
  do ig = gstart, ngm
     aux (ig) = aux (ig) * exp ( - gg (ig) * tpiba2 / alpha_g / 4.d0 ) &
        / (gg (ig) * tpiba2)
  enddo
  
  if (gamma_only) then
     fact = 4.d0
  else
     fact = 2.d0
  end if
  do na = 1, nat
     do ig = gstart, ngm
        arg = tpi * (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
           + g (3, ig) * tau (3, na) )
        sumnb = cos (arg) * AIMAG (aux(ig)) - sin (arg) *  DBLE (aux(ig) )
        forceion (1, na) = forceion (1, na) + g (1, ig) * sumnb
        forceion (2, na) = forceion (2, na) + g (2, ig) * sumnb
        forceion (3, na) = forceion (3, na) + g (3, ig) * sumnb 
     enddo
     do ipol = 1, 3
        forceion (ipol, na) = - zv (ityp (na) ) * fact * e2 * tpi**2 / &
           omega / alat * forceion (ipol, na)
     enddo
  enddo
  deallocate (aux)

  return
END SUBROUTINE esm_force_ewg_pbc

SUBROUTINE esm_force_ewg_bc1 ( alpha_g, forceion )

  USE constants,        ONLY : tpi, fpi, e2
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
  USE gvect,            ONLY : gstart, ngm, g

  implicit none
  real(DP), intent(in)    :: alpha_g
  real(DP), intent(out)   :: forceion(3,nat) 
  !
  !    here the local variables
  !
  integer                 :: it1, it2, k1, k2, ng_2d
  real(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                             c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                             ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                             arg001, arg002, arg101, arg102

  forceion(:,:)=0.d0
  for_g(:,:)=0.d0
  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  z1=z0+esm_w
  tmp=sqrt(alpha_g)

!$omp parallel private( for, it1, it2, z, zp, t1_for, t2_for, &
!$omp                   kk1_for, kk2_for, c1_for, c2_for, &
!$omp                   ng_2d, k1, k2, t, gp2, gp, ff, t1, t2, &
!$omp                   arg001, arg002, arg101, arg102 )
  for=0.d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
     else
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     endif
     t2_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc1
     kk1_for=0.5d0*qe_erf(tmp*(z-zp))
     kk2_for=0.d0

     c1_for(:)=0.d0; c2_for(:)=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if(k1==0.and.k2==0) cycle
        t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
        gp2=sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff=((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc1
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg101=gp/2.d0/tmp-tmp*(z-zp)
        arg102=gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        c1_for(1)=c1_for(1)+sin(ff)*(t1+t2)/4.d0/gp*k1
        c1_for(2)=c1_for(2)+sin(ff)*(t1+t2)/4.d0/gp*k2
        c1_for(3)=c1_for(3)+cos(ff)*(t1-t2)/4.d0
     enddo
     for(:,it2)=for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
     if(gstart==2) then
        for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
     endif
     
  enddo
  enddo
!$omp end do
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
!$omp end parallel

  for_g(:,:)=for_g(:,:)*e2 ! factor e2: hartree -> Ry.

  do it1=1,nat
     forceion(1,it1)=-sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=-sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
     forceion(3,it1)=-for_g(3,it1)
  enddo
  
  return
END SUBROUTINE esm_force_ewg_bc1

SUBROUTINE esm_force_ewg_bc2 ( alpha_g, forceion )

  USE constants,        ONLY : tpi, fpi, e2
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
  USE gvect,            ONLY : gstart, ngm, g

  implicit none
  real(DP), intent(in)    :: alpha_g
  real(DP), intent(out)   :: forceion(3,nat) 
  !
  !    here the local variables
  !
  integer                 :: it1, it2, k1, k2, ng_2d
  real(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                             c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                             ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                             arg001, arg002, arg003, arg004, arg005, &
                             arg006, arg007, arg101, arg102

  forceion(:,:)=0.d0
  for_g(:,:)=0.d0
  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  z1=z0+esm_w
  tmp=sqrt(alpha_g)

!$omp parallel private( for, it1, it2, z, zp, t1_for, t2_for, &
!$omp                   kk1_for, kk2_for, c1_for, c2_for, ng_2d, k1, k2, t, &
!$omp                   gp2, gp, ff, t1, t2, arg001, arg002, arg003, arg004, &
!$omp                   arg005, arg006, arg007, arg101, arg102 )
  for=0.d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
     else
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     endif
     t2_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc2
     kk1_for=0.5d0*qe_erf(tmp*(z-zp))
     kk2_for=-0.5d0*(z/z1)

     c1_for(:)=0.d0; c2_for(:)=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if(k1==0.and.k2==0) cycle
        t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
        gp2=sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff=((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc2
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg003=-gp*(z+zp+2.d0*z1)
        arg004= gp*(z+zp-2.d0*z1)
        arg005=-gp*(z-zp+4.d0*z1)
        arg006= gp*(z-zp-4.d0*z1)
        arg007=-4.d0*gp*z1
        arg101=gp/2.d0/tmp-tmp*(z-zp)
        arg102=gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        c1_for(1)=c1_for(1)+sin(ff)*(t1+t2)/4.d0/gp*k1
        c1_for(2)=c1_for(2)+sin(ff)*(t1+t2)/4.d0/gp*k2
        c1_for(3)=c1_for(3)+cos(ff)*(t1-t2)/4.d0
        c2_for(1)=c2_for(1)+sin(ff)*(exp(arg006)+exp(arg005) &
           -exp(arg004)-exp(arg003))/(1.d0-exp(arg007))/2.d0/gp*k1
        c2_for(2)=c2_for(2)+sin(ff)*(exp(arg006)+exp(arg005) &
           -exp(arg004)-exp(arg003))/(1.d0-exp(arg007))/2.d0/gp*k2
        c2_for(3)=c2_for(3)-cos(ff)*(exp(arg006)-exp(arg005) &
           +exp(arg004)-exp(arg003))/(1.d0-exp(arg007))/2.d0
     enddo
     for(:,it2)=for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
     if(gstart==2) then
        for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
     endif
     
  enddo
  enddo
!$omp end do
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
!$omp end parallel

  for_g(:,:)=for_g(:,:)*e2 ! factor e2: hartree -> Ry.

  do it1=1,nat
     forceion(1,it1)=-sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=-sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
     forceion(3,it1)=-for_g(3,it1)
  enddo
  
  return
END SUBROUTINE esm_force_ewg_bc2

SUBROUTINE esm_force_ewg_bc3 ( alpha_g, forceion )

  USE constants,        ONLY : tpi, fpi, e2
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
  USE gvect,            ONLY : gstart, ngm, g

  implicit none
  real(DP), intent(in)    :: alpha_g
  real(DP), intent(out)   :: forceion(3,nat) 
  !
  !    here the local variables
  !
  integer                 :: it1, it2, k1, k2, ng_2d
  real(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                             c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                             ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                             arg001, arg002, arg003, arg101, arg102

  forceion(:,:)=0.d0
  for_g(:,:)=0.d0
  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  z1=z0+esm_w
  tmp=sqrt(alpha_g)

!$omp parallel private( for, it1, it2, z, zp, t1_for, t2_for, kk1_for, &
!$omp                   kk2_for, c1_for, c2_for, ng_2d, k1, k2, t, gp2, &
!$omp                   gp, ff, t1, t2, arg001, arg002, arg003, arg101, arg102 )
  for=0.d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
     else
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     endif
     t2_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc3
     kk1_for=0.5d0*qe_erf(tmp*(z-zp))
     kk2_for=-0.5d0

     c1_for(:)=0.d0; c2_for(:)=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if(k1==0.and.k2==0) cycle
        t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
        gp2=sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff=((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc3
        arg001=-gp*(z-zp)
        arg002= gp*(z-zp)
        arg003= gp*(z+zp-2.d0*z1)
        arg101= gp/2.d0/tmp-tmp*(z-zp)
        arg102= gp/2.d0/tmp+tmp*(z-zp)
        t1=exp_erfc(arg001,arg101)
        t2=exp_erfc(arg002,arg102)
        c1_for(1)=c1_for(1)+sin(ff)*(t1+t2)/4.d0/gp*k1
        c1_for(2)=c1_for(2)+sin(ff)*(t1+t2)/4.d0/gp*k2
        c1_for(3)=c1_for(3)+cos(ff)*(t1-t2)/4.d0
        c2_for(1)=c2_for(1)+sin(ff)*(-exp(arg003))/2.d0/gp*k1
        c2_for(2)=c2_for(2)+sin(ff)*(-exp(arg003))/2.d0/gp*k2
        c2_for(3)=c2_for(3)+cos(ff)*(-exp(arg003))/2.d0
     enddo
     for(:,it2)=for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
     if(gstart==2) then
        for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
     endif
     
  enddo
  enddo
!$omp end do
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
!$omp end parallel

  for_g(:,:)=for_g(:,:)*e2 ! factor e2: hartree -> Ry.

  do it1=1,nat
     forceion(1,it1)=-sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=-sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
     forceion(3,it1)=-for_g(3,it1)
  enddo
  
  return
END SUBROUTINE esm_force_ewg_bc3

SUBROUTINE esm_force_ewg_bc4 ( alpha_g, forceion )

  USE constants,        ONLY : tpi, fpi, e2
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
  USE gvect,            ONLY : gstart, ngm, g

  implicit none
  real(DP), intent(in)    :: alpha_g
  real(DP), intent(out)   :: forceion(3,nat) 
  !
  !    here the local variables
  !
  integer                 :: it1, it2, k1, k2, ng_2d
  real(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                             c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                             ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                             arg001, arg002, arg003, arg004, arg005, &
                             arg006, arg007, arg008, arg009, arg010, &
                             arg011, arg012, arg101, arg102, arg103, &
                             arg104, arg105, arg106, arg107, arg108, &
                             arg109, arg110, arg111, arg112, arg113, &
                             arg114, aaa, t3, alpha, beta, kappa, lambda, &
                             xi, chi
  ! auxiliary space

  forceion(:,:)=0.d0
  for_g(:,:)=0.d0
  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  aaa=esm_a
  z1=z0+esm_w
  tmp=sqrt(alpha_g)

!$omp parallel private( for, it1, it2, z, zp, t1_for, t2_for, &
!$omp                   kk1_for, kk2_for, c1_for, c2_for, &
!$omp                   ng_2d, k1, k2, t, gp2, &
!$omp                   gp, ff, t1, t2, t3, arg001, arg002, arg003, arg004, &
!$omp                   arg005, arg006, arg007, arg008, arg009, arg010, &
!$omp                   arg011, arg012, arg101, arg102, arg103, arg104, &
!$omp                   arg105, arg106, arg107, arg108, arg109, arg110, &
!$omp                   arg111, arg112, arg113, arg114, alpha, beta, kappa, &
!$omp                   xi, chi, lambda )
  for=0.d0
!$omp do
  do it1=1,nat
  do it2=1,nat
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
     else
        t1_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     endif
     t2_for=zv(ityp(it1))*zv(ityp(it2))*fpi/sa
     ! bc4
     arg004=-2.d0*aaa*(zp-z1)
     arg006= aaa**2/tmp**2+2.d0*aaa*(z1-zp)
     arg101= tmp*(z-zp)
     arg102= tmp*(z1-zp)
     arg104= aaa/tmp+tmp*(z-zp)
     arg106= aaa/tmp+tmp*(z1-zp)
     if (z < z1)then  ! factor 1/2 <- non-reciprocality
        if (zp < z1)then
           kk1_for= 0.5d0*(qe_erf(arg101)-qe_erf(arg102))/2.d0 &
              -0.5d0*exp_erfc(arg006,arg106)/2.d0
           kk2_for=-0.5d0*qe_erfc(arg101)/2.d0
        else
           kk1_for= 0.5d0*(qe_erf(arg101)-qe_erf(arg102))/2.d0 &
              -0.5d0*exp_erfc(arg006,arg106)/2.d0
           kk2_for=-0.5d0*exp_erfc(arg004,arg101)/2.d0
        endif
     else
        if ( zp < z1 )then
           kk1_for=-0.5d0*exp_erfc(arg006,arg104)/2.d0
           kk2_for=-0.5d0*qe_erfc(arg101)/2.d0
        else
           kk1_for=-0.5d0*exp_erfc(arg006,arg104)/2.d0
           kk2_for=-0.5d0*exp_erfc(arg004,arg101)/2.d0
        endif
     endif
     c1_for(:)=0.d0; c2_for(:)=0.d0
     do ng_2d = 1, ngm_2d
        k1 = mill_2d(1,ng_2d)
        k2 = mill_2d(2,ng_2d)
        if(k1==0.and.k2==0) cycle
        t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
        gp2=sum(t(:)*t(:))*tpiba2
        gp=sqrt(gp2)
        ff=((k1*bg(1,1)+k2*bg(1,2))*(tau(1,it1)-tau(1,it2))  &
           +(k1*bg(2,1)+k2*bg(2,2))*(tau(2,it1)-tau(2,it2)))*tpi
        ! bc4
        alpha=aaa+gp+sqrt(aaa**2+gp**2)
        beta =aaa+gp-sqrt(aaa**2+gp**2)
        kappa=aaa-gp+sqrt(aaa**2+gp**2)
        xi   =aaa   +sqrt(aaa**2+gp**2)
        chi  =aaa   -sqrt(aaa**2+gp**2)
        lambda=      sqrt(aaa**2+gp**2)
        arg001= gp*(z-zp)
        arg002=-gp*(z-zp)
        arg003= gp*(z+zp-2.d0*z1)
        arg004= gp*(z-z1)+xi*(z1-zp)
        arg005=-gp*(z1-zp)-xi*(z-z1)
        arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
        arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
        arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
        arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
        arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
        arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
        arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
        arg101=  gp/2.d0/tmp+tmp*(z-zp)
        arg102=  gp/2.d0/tmp-tmp*(z-zp)
        arg103=  gp/2.d0/tmp+tmp*(z1-zp)
        arg104=  gp/2.d0/tmp-tmp*(z1-zp)
        arg105=  gp/2.d0/tmp+tmp*(z-z1)
        arg106=  gp/2.d0/tmp-tmp*(z-z1)
        arg107=  xi/2.d0/tmp+tmp*(z-zp)
        arg108=  xi/2.d0/tmp-tmp*(z-zp)
        arg109=  xi/2.d0/tmp+tmp*(z1-zp)
        arg110=  xi/2.d0/tmp-tmp*(z-z1)
        arg111= chi/2.d0/tmp+tmp*(z-zp)
        arg112= chi/2.d0/tmp-tmp*(z-zp)
        arg113= chi/2.d0/tmp+tmp*(z1-zp)
        arg114= chi/2.d0/tmp-tmp*(z-z1)
        if (z < z1) then ! factor 1/2 <- non-reciprocality
           if (zp < z1) then
              t1= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3= exp_erfc(arg006,arg109)/alpha
              c1_for(1)=c1_for(1)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
              c1_for(2)=c1_for(2)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
              t2= exp_erfc(arg001,arg101) &
                 -kappa/alpha*exp_erfc(arg003,arg105)
              t3= exp_erfc(arg007,arg110)/alpha
              c2_for(1)=c2_for(1)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
              c2_for(2)=c2_for(2)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3=-xi/alpha*exp_erfc(arg006,arg109)
              c1_for(3)=c1_for(3)+cos(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0 
              t1= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
              t2=-exp_erfc(arg001,arg101) &
                 -kappa/alpha*exp_erfc(arg003,arg105)
              t3= gp/alpha*exp_erfc(arg007,arg110)
              c2_for(3)=c2_for(3)+cos(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0
           else
              t1= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3= exp_erfc(arg006,arg109)/alpha
              c1_for(1)=c1_for(1)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
              c1_for(2)=c1_for(2)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112)
              t2= exp_erfc(arg010,arg108) &
                 -beta/alpha*exp_erfc(arg009,arg110)
              t3= exp_erfc(arg004,arg105)/alpha
              c2_for(1)=c2_for(1)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k1/2.d0
              c2_for(2)=c2_for(2)+sin(ff)*((t1+t2)/4.d0/gp+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3=-xi/alpha*exp_erfc(arg006,arg109)
              c1_for(3)=c1_for(3)+cos(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0 
              t1= xi*(exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))
              t2=-chi*exp_erfc(arg010,arg108) &
                 +xi*beta/alpha*exp_erfc(arg009,arg110)
              t3=-xi/alpha*exp_erfc(arg004,arg105)
              c2_for(3)=c2_for(3)+cos(ff)*((t1+t2)/4.d0+t3/2.d0)/2.d0
           endif
        else
           if(zp < z1)then
              t1= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
              t2= exp_erfc(arg008,arg107) &
                 -beta/alpha*exp_erfc(arg009,arg109)
              t3= exp_erfc(arg005,arg104)/alpha
              c1_for(1)=c1_for(1)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
              c1_for(2)=c1_for(2)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
              t2= exp_erfc(arg001,arg101) &
                 -kappa/alpha*exp_erfc(arg003,arg105)
              t3= exp_erfc(arg007,arg110)/alpha
              c2_for(1)=c2_for(1)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
              c2_for(2)=c2_for(2)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
              t1= chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
              t2=-xi*exp_erfc(arg008,arg107) &
                 +xi*beta/alpha*exp_erfc(arg009,arg109)
              t3= gp/alpha*exp_erfc(arg005,arg104)
              c1_for(3)=c1_for(3)+cos(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
              t1= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106)
              t2=-exp_erfc(arg001,arg101) &
                 -kappa/alpha*exp_erfc(arg003,arg105)
              t3= gp/alpha*exp_erfc(arg007,arg110)
              c2_for(3)=c2_for(3)+cos(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
            else
              t1= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
              t2= exp_erfc(arg008,arg107) &
                 -beta/alpha*exp_erfc(arg009,arg109)
              t3= exp_erfc(arg005,arg104)/alpha
              c1_for(1)=c1_for(1)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
              c1_for(2)=c1_for(2)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
              t1= exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112)
              t2= exp_erfc(arg010,arg108) &
                 -beta/alpha*exp_erfc(arg009,arg110)
              t3= exp_erfc(arg004,arg105)/alpha
              c2_for(1)=c2_for(1)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k1/2.d0
              c2_for(2)=c2_for(2)+sin(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)*k2/2.d0
              t1= chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
              t2=-xi*exp_erfc(arg008,arg107) &
                 +xi*beta/alpha*exp_erfc(arg009,arg109)
              t3= gp/alpha*exp_erfc(arg005,arg104)
              c1_for(3)=c1_for(3)+cos(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
              t1= xi*(exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))
              t2=-chi*exp_erfc(arg010,arg108) &
                 +xi*beta/alpha*exp_erfc(arg009,arg110)
              t3=-xi/alpha*exp_erfc(arg004,arg105)
              c2_for(3)=c2_for(3)+cos(ff)*((t1+t2)/4.d0/lambda+t3/2.d0)/2.d0
           endif
        endif
     enddo
     for(:,it2)=for(:,it2)+t1_for*(c1_for(:)+c2_for(:))
     if(gstart==2) then
        for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
     endif
     
  enddo
  enddo
!$omp end do
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
!$omp end parallel

  for_g(:,:)=for_g(:,:)*e2 ! factor e2: hartree -> Ry.

  do it1=1,nat
     forceion(1,it1)=-sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=-sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
     forceion(3,it1)=-for_g(3,it1)
  enddo
  
  return
END SUBROUTINE esm_force_ewg_bc4


!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL-DERIVED FORCE SUBROUTINE-------------
!-----------------------------------------------------------------------
SUBROUTINE esm_force_lc_pbc( aux, forcelc )
  USE ions_base, ONLY : nat
  USE fft_base,  ONLY : dfftp
  implicit none
  complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
  real(DP),    intent(inout) :: forcelc(3,nat)

  stop 'esm_force_lc must not be called for esm_bc = pbc'

END SUBROUTINE esm_force_lc_pbc

SUBROUTINE esm_force_lc_bc1 ( aux, forcelc )

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z

  implicit none
  complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)   
  real(DP),    intent(inout) :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
  real(DP),allocatable    :: for(:,:), for_g(:,:)
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                             t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                             arg001, arg002, arg101, arg102
  complex(DP),allocatable :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
  complex(DP)             :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
  allocate(rhog3(dfftp%nr3,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( n1, n2, ng_2d, n3 )
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
!$omp end parallel do

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0
  z1=z0+esm_w

  allocate(for_g(3,nat))
  for_g(:,:)=0.d0

!**** for gp!=0 *********
!$omp parallel private( for, vg_f, vg_f_r, ng_2d, k1, k2, k3, &
!$omp                   gp2, gp, it, tt, pp, cc, ss, zp, iz, z, t1, t2, &
!$omp                   c1, c2, r1, r2, f1, &
!$omp                   f2, arg001, arg002, arg101, arg102 )
  allocate(for(3,nat),vg_f(dfftp%nr3x,3), vg_f_r(dfftp%nr3x,3))
  for(:,:)=0.d0
  vg_f_r(:,:)=(0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        if (gamma_only) then
           tt=-fpi*zv(ityp(it))/sa*2.d0
        else 
           tt=-fpi*zv(ityp(it))/sa
        endif 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc1
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           c1(1)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k1
           c1(2)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k2
           c1(3)=CMPLX(cc, ss,kind=DP)*(t1-t2)/4.d0
           c2(:)=(0.d0,0.d0)
           vg_f_r(iz,:) = tt*(c1(:)+c2(:))
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        call esm_cft_1z(vg_f_r(:,2),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,2))
        call esm_cft_1z(vg_f_r(:,3),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,3))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(:)= dble(vg_f(iz,:))
           f2(:)=aimag(vg_f(iz,:))
           for(:,it)=for(:,it)-r1*f1(:)-r2*f2(:)
        enddo
     enddo
  enddo
!$omp critical
  for_g(:,:)=for_g(:,:)+for(:,:)
!$omp end critical
  deallocate(for,vg_f,vg_f_r)
!$omp end parallel

!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg_f(dfftp%nr3x,1),vg_f_r(dfftp%nr3x,1))
     vg_f_r(:,1)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc1
           cc1=0.5d0*qe_erf(tmp*(z-zp))
           cc2=(0.d0,0.d0)

           vg_f_r(iz,1) = tt*(cc1+cc2)
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(3)= dble(vg_f(iz,1))
           f2(3)=aimag(vg_f(iz,1))
           for_g(3,it)=for_g(3,it)-r1*f1(3)-r2*f2(3)
        enddo
     enddo

     deallocate(vg_f,vg_f_r)
  endif ! if( ng_2d > 0 )


!***** sum short_range part and long_range part in local potential force 
!***** at cartecian coordinate

  do it=1,nat
                                             ! factor e2: hartree -> Ry.
     forcelc(1,it)=forcelc(1,it) &
                  +sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(2,it)=forcelc(2,it) &
                  +sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*e2
  enddo

  deallocate(for_g)

  call setlocal()

  deallocate(rhog3)
  return
END SUBROUTINE esm_force_lc_bc1

SUBROUTINE esm_force_lc_bc2 ( aux, forcelc )

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z

  implicit none
  complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)   
  real(DP),    intent(inout) :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
  real(DP),allocatable    :: for(:,:), for_g(:,:)
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                             t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                             arg001, arg002, arg003, arg005, &
                             arg006, arg008, arg009, arg101, arg102
  complex(DP),allocatable :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
  complex(DP)             :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
  allocate(rhog3(dfftp%nr3,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( n1, n2, ng_2d, n3 )
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
!$omp end parallel do

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0
  z1=z0+esm_w

  allocate(for_g(3,nat))
  for_g(:,:)=0.d0

!**** for gp!=0 *********
!$omp parallel private( for, vg_f, vg_f_r, ng_2d, k1, k2, k3, &
!$omp                   gp2, gp, it, tt, pp, cc, ss, zp, iz, z, t1, t2, &
!$omp                   c1, c2, r1, r2, f1, f2, arg001, arg002, arg003, &
!$omp                   arg005, arg006, arg008, arg009, arg101, arg102 )
  allocate(for(3,nat),vg_f(dfftp%nr3x,3),vg_f_r(dfftp%nr3x,3))
  for(:,:)=0.d0
  vg_f_r(:,:)=(0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        if (gamma_only) then
           tt=-fpi*zv(ityp(it))/sa*2.d0
        else 
           tt=-fpi*zv(ityp(it))/sa
        endif 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc2
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003=-gp*(z+zp+2.d0*z1)
           arg005= gp*(z+zp-2.d0*z1)
           arg006=-gp*(z-zp+4.d0*z1)
           arg008= gp*(z-zp-4.d0*z1)
           arg009=-4.d0*gp*z1
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           c1(1)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k1
           c1(2)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k2
           c1(3)=CMPLX(cc, ss,kind=DP)*(t1-t2)/4.d0
           c2(1)=CMPLX(ss,-cc,kind=DP)*(exp(arg008)+exp(arg006) &
              -exp(arg005)-exp(arg003))/(1.d0-exp(arg009))/2.d0/gp*k1
           c2(2)=CMPLX(ss,-cc,kind=DP)*(exp(arg008)+exp(arg006) &
              -exp(arg005)-exp(arg003))/(1.d0-exp(arg009))/2.d0/gp*k2
           c2(3)=CMPLX(cc, ss,kind=DP)*(-exp(arg008)+exp(arg006) &
              -exp(arg005)+exp(arg003))/(1.d0-exp(arg009))/2.d0
           vg_f_r(iz,:) = tt*(c1(:)+c2(:))
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        call esm_cft_1z(vg_f_r(:,2),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,2))
        call esm_cft_1z(vg_f_r(:,3),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,3))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(:)= dble(vg_f(iz,:))
           f2(:)=aimag(vg_f(iz,:))
           for(:,it)=for(:,it)-r1*f1(:)-r2*f2(:)
        enddo
     enddo
  enddo
!$omp critical
  for_g(:,:)=for_g(:,:)+for(:,:)
!$omp end critical
  deallocate(for,vg_f,vg_f_r)
!$omp end parallel

!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg_f(dfftp%nr3x,1),vg_f_r(dfftp%nr3x,1))
     vg_f_r(:,1)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc2
           cc1=0.5d0*qe_erf(tmp*(z-zp))
           cc2=-0.5d0*(z/z1)
           vg_f_r(iz,1) = tt*(cc1+cc2)
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(3)= dble(vg_f(iz,1))
           f2(3)=aimag(vg_f(iz,1))
           for_g(3,it)=for_g(3,it)-r1*f1(3)-r2*f2(3)
        enddo
     enddo
     deallocate(vg_f,vg_f_r)
  endif ! if( ng_2d > 0 )


!***** sum short_range part and long_range part in local potential force 
!***** at cartecian coordinate

  do it=1,nat
                                             ! factor e2: hartree -> Ry.
     forcelc(1,it)=forcelc(1,it) &
                  +sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(2,it)=forcelc(2,it) &
                  +sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*e2
  enddo

  deallocate(for_g)

  call setlocal()

  deallocate(rhog3)
  return
END SUBROUTINE esm_force_lc_bc2

SUBROUTINE esm_force_lc_bc3 ( aux, forcelc )

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z

  implicit none
  complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)   
  real(DP),    intent(inout) :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
  real(DP),allocatable    :: for(:,:), for_g(:,:)
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                             t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                             arg001, arg002, arg003, arg101, arg102
  complex(DP),allocatable :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
  complex(DP)             :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
  allocate(rhog3(dfftp%nr3,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( n1, n2, ng_2d, n3 )
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
!$omp end parallel do

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0
  z1=z0+esm_w

  allocate(for_g(3,nat))
  for_g(:,:)=0.d0

!**** for gp!=0 *********
!$omp parallel private( for, vg_f, vg_f_r,  ng_2d, k1, k2, k3, &
!$omp                   gp2, gp, it, tt, pp, cc, ss, zp, iz, z, t1, t2, &
!$omp                   c1, c2, r1, r2, f1, f2, &
!$omp                   arg001, arg002, arg003, arg101, arg102 )
  allocate(for(3,nat),vg_f(dfftp%nr3x,3),vg_f_r(dfftp%nr3x,3))
  for(:,:)=0.d0
  vg_f_r(:,:)=(0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        if (gamma_only) then
           tt=-fpi*zv(ityp(it))/sa*2.d0
        else 
           tt=-fpi*zv(ityp(it))/sa
        endif 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc3
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003= gp*(z+zp-2.d0*z1)
           arg101= gp/2.d0/tmp+tmp*(z-zp)
           arg102= gp/2.d0/tmp-tmp*(z-zp)
           t1=exp_erfc(arg002,arg102)
           t2=exp_erfc(arg001,arg101)
           c1(1)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k1
           c1(2)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k2
           c1(3)=CMPLX(cc, ss,kind=DP)*(t1-t2)/4.d0
           c2(1)=CMPLX(ss,-cc,kind=DP)*(-exp(arg003))/2.d0/gp*k1
           c2(2)=CMPLX(ss,-cc,kind=DP)*(-exp(arg003))/2.d0/gp*k2
           c2(3)=CMPLX(cc, ss,kind=DP)*(-exp(arg003))/2.d0

           vg_f_r(iz,:) = tt*(c1(:)+c2(:))
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        call esm_cft_1z(vg_f_r(:,2),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,2))
        call esm_cft_1z(vg_f_r(:,3),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,3))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(:)= dble(vg_f(iz,:))
           f2(:)=aimag(vg_f(iz,:))
           for(:,it)=for(:,it)-r1*f1(:)-r2*f2(:)
        enddo
     enddo
  enddo
!$omp critical
  for_g(:,:)=for_g(:,:)+for(:,:)
!$omp end critical
  deallocate(for,vg_f,vg_f_r)
!$omp end parallel

!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg_f(dfftp%nr3x,1),vg_f_r(dfftp%nr3x,1))
     vg_f_r(:,1)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc3
           cc1=0.5d0*qe_erf(tmp*(z-zp))
           cc2=-0.5d0
           vg_f_r(iz,1) = tt*(cc1+cc2)
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(3)= dble(vg_f(iz,1))
           f2(3)=aimag(vg_f(iz,1))
           for_g(3,it)=for_g(3,it)-r1*f1(3)-r2*f2(3)
        enddo
     enddo

     deallocate(vg_f,vg_f_r)
  endif ! if( ng_2d > 0 )


!***** sum short_range part and long_range part in local potential force 
!***** at cartecian coordinate

  do it=1,nat
                                             ! factor e2: hartree -> Ry.
     forcelc(1,it)=forcelc(1,it) &
                  +sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(2,it)=forcelc(2,it) &
                  +sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*e2
  enddo

  deallocate(for_g)

  call setlocal()

  deallocate(rhog3)
  return
END SUBROUTINE esm_force_lc_bc3

SUBROUTINE esm_force_lc_bc4 ( aux, forcelc )

  USE constants,        ONLY : tpi, fpi, e2
  USE gvect,            ONLY : ngm, nl, nlm, mill
  USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE control_flags,    ONLY : gamma_only
  USE ions_base,        ONLY : zv, nat, tau, ityp
  USE fft_base,         ONLY : dfftp
  USE esm_cft,          ONLY : esm_cft_1z

  implicit none
  complex(DP), intent(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)   
  real(DP),    intent(inout) :: forcelc(3,nat)
  !
  !    here are the local variables
  !
  integer                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
  real(DP),allocatable    :: for(:,:), for_g(:,:)
  real(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                             t2, z, zp, L, tmp, r1, r2, f1(3), &
                             f2(3), arg001, arg002, arg003, arg005, &
                             arg006, arg008, arg009, arg011, arg101, &
                             arg102, arg103, arg104, arg106, arg107, &
                             arg109, arg111, arg113, aaa, t3, alpha, beta, &
                             kappa, lambda, xi, chi
  complex(DP),allocatable :: vg_f(:,:), vg_f_r(:,:), rhog3(:,:)
  complex(DP)             :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
  allocate(rhog3(dfftp%nr3,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
!$omp parallel do private( n1, n2, ng_2d, n3 )
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
!$omp end parallel do

  L=at(3,3)*alat
  sa=omega/L
  z0=L/2.d0
  tmp=1.d0
  z1=z0+esm_w
  aaa=esm_a

  allocate(for_g(3,nat))
  for_g(:,:)=0.d0

!**** for gp!=0 *********
!$omp parallel private( for, vg_f, vg_f_r, ng_2d, k1, k2, k3, &
!$omp                   gp2, gp, it, tt, pp, cc, ss, zp, iz, z, t1, t2, t3, &
!$omp                   c1, c2, r1, r2, f1, f2, &
!$omp                   arg001, arg002, arg003, arg005, &
!$omp                   arg006, arg008, arg009, arg011, arg101, arg102, &
!$omp                   arg103, arg104, arg107, arg109, arg111, arg113, &
!$omp                   alpha, beta, kappa, xi, chi, lambda )
  allocate(for(3,nat),vg_f(dfftp%nr3x,3),vg_f_r(dfftp%nr3x,3))
  for(:,:)=0.d0
  vg_f_r(:,:)=(0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2)=k1*bg(1:2,1)+k2*bg(1:2,2)
     gp2=sum(t(:)*t(:))*tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        if (gamma_only) then
           tt=-fpi*zv(ityp(it))/sa*2.d0
        else 
           tt=-fpi*zv(ityp(it))/sa
        endif 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2)) &
                +tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc4
           alpha = aaa+gp+sqrt(aaa**2+gp**2)
           beta  = aaa+gp-sqrt(aaa**2+gp**2)
           kappa = aaa-gp+sqrt(aaa**2+gp**2)
           xi    = aaa   +sqrt(aaa**2+gp**2)
           chi   = aaa   -sqrt(aaa**2+gp**2)
           lambda=        sqrt(aaa**2+gp**2)
           arg001= gp*(z-zp)
           arg002=-gp*(z-zp)
           arg003= gp*(z+zp-2.d0*z1)
           arg005=-gp*(z1-zp)-xi*(z-z1)
           arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
           arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
           arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
           arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
           arg101=  gp/2.d0/tmp+tmp*(z-zp)
           arg102=  gp/2.d0/tmp-tmp*(z-zp)
           arg103=  gp/2.d0/tmp+tmp*(z1-zp)
           arg104=  gp/2.d0/tmp-tmp*(z1-zp)
           arg107=  xi/2.d0/tmp+tmp*(z-zp)
           arg109=  xi/2.d0/tmp+tmp*(z1-zp)
           arg111= chi/2.d0/tmp+tmp*(z-zp)
           arg113= chi/2.d0/tmp+tmp*(z1-zp)
           if (z < z1) then
              t1= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3= exp_erfc(arg006,arg109)/alpha
              c1(1)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k1
              c1(2)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/gp*k2
              c2(1)=CMPLX(ss,-cc,kind=DP)*t3/2.d0*k1
              c2(2)=CMPLX(ss,-cc,kind=DP)*t3/2.d0*k2
              t1= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101)
              t2= exp_erfc(arg002,arg102) &
                 -kappa/alpha*exp_erfc(arg003,arg104)
              t3=-xi/alpha*exp_erfc(arg006,arg109)
              c1(3)=CMPLX(cc, ss,kind=DP)*(t1+t2)/4.d0
              c2(3)=CMPLX(cc, ss,kind=DP)*t3/2.d0
           else
              t1= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111)
              t2= exp_erfc(arg008,arg107) &
                 -beta/alpha*exp_erfc(arg009,arg109)
              t3= exp_erfc(arg005,arg104)/alpha
              c1(1)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/lambda*k1
              c1(2)=CMPLX(ss,-cc,kind=DP)*(t1+t2)/4.d0/lambda*k2
              c2(1)=CMPLX(ss,-cc,kind=DP)*t3/2.d0*k1
              c2(2)=CMPLX(ss,-cc,kind=DP)*t3/2.d0*k2
              t1= chi*(exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))
              t2=-xi*(exp_erfc(arg008,arg107) &
                 +beta/alpha*exp_erfc(arg009,arg109))
              t3= gp/alpha*exp_erfc(arg005,arg104)
              c1(3)=CMPLX(cc, ss,kind=DP)*(t1+t2)/4.d0/lambda
              c2(3)=CMPLX(cc, ss,kind=DP)*t3/2.d0
           endif
           vg_f_r(iz,:) = tt*(c1(:)+c2(:))
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        call esm_cft_1z(vg_f_r(:,2),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,2))
        call esm_cft_1z(vg_f_r(:,3),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,3))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(:)= dble(vg_f(iz,:))
           f2(:)=aimag(vg_f(iz,:))
           for(:,it)=for(:,it)-r1*f1(:)-r2*f2(:)
        enddo
     enddo
  enddo
!$omp critical
  for_g(:,:)=for_g(:,:)+for(:,:)
!$omp end critical
  deallocate(for,vg_f,vg_f_r)
!$omp end parallel

!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg_f(dfftp%nr3x,1),vg_f_r(dfftp%nr3x,1))
     vg_f_r(:,1)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*zv(ityp(it))/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,dfftp%nr3
           k3=iz-1
           if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
           z=dble(k3)/dble(dfftp%nr3)*L
           ! bc4
           arg006= aaa**2/tmp**2+2.d0*aaa*(z1-zp)
           arg101= tmp*(z-zp)
           arg102= tmp*(z1-zp)
           arg104= aaa/tmp+tmp*(z-zp)
           arg106= aaa/tmp+tmp*(z1-zp)
           if (z < z1)then
              cc1= 0.5d0*(qe_erf(arg101)-qe_erf(arg102))
              cc2=-0.5d0*exp_erfc(arg006,arg106)
           else
              cc1= 0.d0
              cc2=-0.5d0*exp_erfc(arg006,arg104)
           endif
           vg_f_r(iz,1) = tt*(cc1+cc2)
        enddo
        call esm_cft_1z(vg_f_r(:,1),1,dfftp%nr3,dfftp%nr3,-1,vg_f(:,1))
        do iz=1,dfftp%nr3
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           f1(3)= dble(vg_f(iz,1))
           f2(3)=aimag(vg_f(iz,1))
           for_g(3,it)=for_g(3,it)-r1*f1(3)-r2*f2(3)
        enddo
     enddo
     deallocate(vg_f,vg_f_r)
  endif ! if( ng_2d > 0 )


!***** sum short_range part and long_range part in local potential force 
!***** at cartecian coordinate

  do it=1,nat
                                             ! factor e2: hartree -> Ry.
     forcelc(1,it)=forcelc(1,it) &
                  +sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(2,it)=forcelc(2,it) &
                  +sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*e2
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*e2
  enddo

  deallocate(for_g)

  call setlocal()

  deallocate(rhog3)
  return
END SUBROUTINE esm_force_lc_bc4

!-----------------------------------------------------------------------
!--------------ESM FINAL PRINTOUT SUBROUTINE----------------------------
!-----------------------------------------------------------------------
!
! Prints out vlocal and vhartree to stdout once electrons are converged
! Format: z, rho(r), v_hartree, v_local, (v_hartree + v_local) 
!
SUBROUTINE esm_printpot ()
  USE constants,            ONLY : pi, tpi, fpi, eps4, eps8, e2
  USE cell_base,            ONLY : at, alat
  USE scf,                  ONLY : rho, vltot
  USE lsda_mod,             ONLY : nspin
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : intra_bgrp_comm
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : prefix, tmp_dir
  USE constants,            ONLY : rytoev, bohr_radius_angs
  !
  IMPLICIT NONE
  !
  real(DP)                :: z1,z2,z3,z4,charge,ehart,L,area
  real(DP),allocatable    :: work1(:),work2(:,:),work3(:), work4(:,:)
  integer                 :: ix,iy,iz,izz,i,k3
  character (len=256)     :: esm1_file = 'os.esm1'

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
  IF ( ionode ) THEN
     esm1_file = TRIM( tmp_dir ) // TRIM( prefix ) // ".esm1"
     open( UNIT = 4, FILE = esm1_file, STATUS = "UNKNOWN", &
           ACTION = "WRITE" )
     !
     write( UNIT = 4, FMT = 9050 )
     do k3 = dfftp%nr3/2-dfftp%nr3+1, dfftp%nr3/2
        iz = k3 + dfftp%nr3 + 1
        if( iz > dfftp%nr3 ) iz = iz - dfftp%nr3
        write( UNIT = 4, FMT = 9051 ) work4(1:5,iz)
     enddo
     close( UNIT = 4 )
  ENDIF
  deallocate(work1,work2,work3,work4)
9050 FORMAT( '#z (A)',2X,'Tot chg (e/A)',2X,'Avg v_hartree (eV)',2X,&
        &'Avg v_local (eV)',2x,'Avg v_hart+v_loc (eV)' )
9051 FORMAT( F6.2,F14.4,F20.7,F18.7,F18.7 )
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
      USE io_global,        ONLY : stdout, ionode
      USE constants,        ONLY : rytoev, BOHR_RADIUS_ANGS
      !
      IMPLICIT NONE
      !
      IF( .NOT. ionode ) RETURN
      !
      WRITE( UNIT = stdout,                                          &
              FMT  = '(/,5x, "Effective Screening Medium Method",     &
                      &/,5x, "=================================")' )
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
      CASE( 'bc4' )
         WRITE( UNIT = stdout,                                     &
              FMT  = '(5x, "Boundary Conditions: Vacuum-Slab-smooth ESM)")' )
      END SELECT
      !
      IF( esm_efield /= 0.0_DP ) THEN
         WRITE( UNIT = stdout, FMT = 9051 ) esm_efield
      END IF
      !
      IF( esm_w /= 0.0_DP ) THEN
         WRITE( UNIT = stdout, FMT = 9052 ) esm_w*BOHR_RADIUS_ANGS, esm_w
      END IF
      !
      IF (esm_bc .EQ. 'bc4') THEN
        WRITE( UNIT = stdout, FMT = 9054 ) esm_a
      ENDIF
      !
      WRITE( UNIT = stdout, FMT = 9053 ) esm_nfit
      !
      WRITE( stdout, * )
      !
9051  FORMAT( '     field strength                   = ', F8.2,' Ry/a.u.')
9052  FORMAT( '     ESM offset from cell edge        = ', F8.2,' A' &
             /'                                      = ', F8.2,' a.u.')
9053  FORMAT( '     grid points for fit at edges     = ', I8,' ')
9054  FORMAT( '     smoothness parameter             = ', F8.2,' 1/a.u.' )

END SUBROUTINE esm_summary

function vl11j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl11j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg006, arg101, arg102, arg103, arg104, &
                         arg109
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   t1=-exp(arg003)*kappa/alpha
   t2= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
      +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
   t3= exp_erfc(arg006,arg109)/alpha
   vl11j0=(t1/2.d0-(t2/4.d0+gp*t3/2.d0))*dbesj0(gp*rxy)
   !
   return
end function vl11j0

function vl11j1(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl11j1
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg006, arg007, arg101, arg102, arg103, arg104, &
                         arg105, arg106, arg109, arg110
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg106=  gp/2.d0/tmp-tmp*(z-z1)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   t1=-exp(arg003)*kappa/alpha
   t2= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
      +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104) &
      +exp_erfc(arg006,arg109)*2.d0*gp/alpha
   t3= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
      +exp_erfc(arg001,arg101)-kappa/alpha*exp_erfc(arg003,arg105) &
      +exp_erfc(arg007,arg110)*2.d0*gp/alpha
   vl11j1=gp*(t1-(t2+t3)/4.d0)*dbesj1(gp*rxy)
   !
   return
end function vl11j1

function vl12j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl12j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg004, arg006, arg101, arg102, arg103, &
                         arg104, arg109
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   t1= exp(arg004)/alpha
   t2= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
      +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
   t3= exp_erfc(arg006,arg109)/alpha
   vl12j0=(gp*t1-(t2/4.d0+gp*t3/2.d0))*dbesj0(gp*rxy)   
end function vl12j0

function vl12j1(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl12j1
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                         arg001, arg002, arg003, arg004, arg006, arg009, &
                         arg010, arg012, arg101, arg102, arg103, arg104, &
                         arg105, arg108, arg109, arg110, arg112, arg114
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
   arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg108=  xi/2.d0/tmp-tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg112= chi/2.d0/tmp-tmp*(z-zp)
   arg114= chi/2.d0/tmp-tmp*(z-z1)
   t1= exp(arg004)/alpha
   t2= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
      +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104) &
      +exp_erfc(arg006,arg109)*2.d0*gp/alpha
   t3= exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112) &
      +exp_erfc(arg010,arg108)-beta/alpha*exp_erfc(arg009,arg110) &
      +exp_erfc(arg004,arg105)*2.d0*lambda/alpha
   vl12j1=(2.d0*t1*gp**2-(gp*t2+gp**2*t3/lambda)/4.d0)*dbesj1(gp*rxy)
end function vl12j1

function vl21j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl21j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg005, &
                         arg008, arg009, arg011, arg104, arg107, arg109, &
                         arg111, arg113
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   t1= exp(arg005)/alpha
   t2= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
      +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
   t3= exp_erfc(arg005,arg104)/alpha
   vl21j0=gp*(t1-(t2/4.d0/lambda+t3/2.d0))*dbesj0(gp*rxy)
end function vl21j0

function vl21j1(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl21j1
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                         arg001, arg002, arg003, arg005, arg007,arg008, &
                         arg009, arg011, arg101, arg102, arg104, arg105, &
                         arg106, arg107, arg109, arg110, arg111, arg113
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg106=  gp/2.d0/tmp-tmp*(z-z1)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   t1= exp(arg005)/alpha
   t2= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
      +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109) &
      +exp_erfc(arg005,arg104)*2.d0*lambda/alpha
   t3= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
      +exp_erfc(arg001,arg101)-kappa/alpha*exp_erfc(arg003,arg105) &
      +exp_erfc(arg007,arg110)*2.d0*gp/alpha
   vl21j1=(2.d0*t1*gp**2-(gp**2*t2/lambda+gp*t3)/4.d0)*dbesj1(gp*rxy)
end function vl21j1

function vl22j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl22j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                         arg005, arg008, arg009, arg011, arg104, arg107, &
                         arg109, arg111, arg113
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg000=-xi*(z+zp-2.d0*z1)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   t1=-exp(arg000)*beta/alpha
   t2= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
      +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
   t3= exp_erfc(arg005,arg104)/alpha
   vl22j0=gp*(t1/2.d0/lambda-(t2/4.d0/lambda+t3/2.d0))*dbesj0(gp*rxy)
   !
   return
end function vl22j0

function vl22j1(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl22j1
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                         arg004, arg005, arg008, arg009, arg010, arg011, &
                         arg012, arg104, arg105, arg107, arg108, arg109, &
                         arg110, arg111, arg112, arg113, arg114
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg000=-xi*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg108=  xi/2.d0/tmp-tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg112= chi/2.d0/tmp-tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   arg114= chi/2.d0/tmp-tmp*(z-z1)
   t1=-exp(arg000)*beta/alpha
   t2= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
      +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109) &
      +exp_erfc(arg005,arg104)*2.d0*lambda/alpha
   t3= exp_erfc(arg012,arg114)-exp_erfc(arg012,arg112) &
      +exp_erfc(arg010,arg108)-beta/alpha*exp_erfc(arg009,arg110) &
      +exp_erfc(arg004,arg105)*2.d0*lambda/alpha
   vl22j1=gp**2*(t1-(t2+t3)/4.d0)*dbesj1(gp*rxy)/lambda
   !
   return
end function vl22j1

function vl11(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl11
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg006, arg101, arg102, arg103, arg104, &
                         arg109
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   t1=-exp(arg003)*kappa/alpha
   t2= exp_erfc(arg001,arg101)-exp_erfc(arg001,arg103) &
      +exp_erfc(arg002,arg102)-kappa/alpha*exp_erfc(arg003,arg104)
   t3= exp_erfc(arg006,arg109)/alpha
   vl11=t1/2.d0-(t2/4.d0+gp*t3/2.d0)
   !
   return
end function vl11

function vl22(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: vl22
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                         arg005, arg008, arg009, arg011, arg104, arg107, &
                         arg109, arg111, arg113
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg000=-xi*(z+zp-2.d0*z1)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   t1=-exp(arg000)*beta/alpha
   t2= exp_erfc(arg011,arg113)-exp_erfc(arg011,arg111) &
      +exp_erfc(arg008,arg107)-beta/alpha*exp_erfc(arg009,arg109)
   t3= exp_erfc(arg005,arg104)/alpha
   vl22=gp*t1/2.d0/lambda-gp*(t2/4.d0/lambda+t3/2.d0)
   !
   return
end function vl22

function dvl11(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl11
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg006, arg007, arg101, arg102, arg103, &
                         arg104, arg105, arg106, arg109, arg110
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg106=  gp/2.d0/tmp-tmp*(z-z1)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   t1=-exp(arg003)*kappa/alpha
   t2= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
      +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
      -exp_erfc(arg006,arg109)*xi/alpha*2.d0
   t3= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
      -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
      +exp_erfc(arg007,arg110)*gp/alpha*2.d0
   dvl11=gp*(t1-(t2+t3)/4.d0)
   !
   return
end function dvl11

function dvl22(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl22
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, arg000, &
                         arg004, arg005, arg008, arg009, arg010, arg011, &
                         arg012, arg104, arg105, arg107, arg108, arg109, &
                         arg110, arg111, arg112, arg113, arg114, t1, t2, t3
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg000=-xi*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg108=  xi/2.d0/tmp-tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg112= chi/2.d0/tmp-tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   arg114= chi/2.d0/tmp-tmp*(z-z1)
   t1= exp(arg000)*beta*xi/alpha/lambda
   t2= (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
      -exp_erfc(arg008,arg107)*xi/lambda &
      +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda &
      +exp_erfc(arg005,arg104)*gp/alpha*2.d0
   t3= (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda &
      -exp_erfc(arg010,arg108)*chi/lambda &
      +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
      -exp_erfc(arg004,arg105)*xi/alpha*2.d0
   dvl22=gp*(t1-(t2+t3)/4.d0)
   !
   return
end function dvl22

function dvl11j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl11j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                         arg003, arg006, arg007, arg101, arg102, arg103, &
                         arg104, arg105, arg106, arg109, arg110
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg106=  gp/2.d0/tmp-tmp*(z-z1)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   t1=-exp(arg003)*kappa/alpha
   t2= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
      +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
      -exp_erfc(arg006,arg109)*xi/alpha*2.d0
   t3= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
      -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
      +exp_erfc(arg007,arg110)*gp/alpha*2.d0
   dvl11j0=gp*(t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
   !
   return
end function dvl11j0

function dvl12j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl12j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, arg001, &
                         arg002, arg003, arg004, arg006, arg009, arg010, &
                         arg012, arg101, arg102, arg103, arg104, arg105, &
                         arg108, arg109, arg110, arg112, arg114, t1, t2, t3
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg006= aaa/2.d0/tmp**2*xi+gp*(z-z1)+xi*(z1-zp)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
   arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg103=  gp/2.d0/tmp+tmp*(z1-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg108=  xi/2.d0/tmp-tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg112= chi/2.d0/tmp-tmp*(z-zp)
   arg114= chi/2.d0/tmp-tmp*(z-z1)
   t1=-exp(arg004)*xi/alpha
   t2= exp_erfc(arg001,arg103)-exp_erfc(arg001,arg101) &
      +exp_erfc(arg002,arg102)-exp_erfc(arg003,arg104)*kappa/alpha &
      -exp_erfc(arg006,arg109)*xi/alpha*2.d0
   t3= (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda &
      -exp_erfc(arg010,arg108)*chi/lambda &
      +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
      -exp_erfc(arg004,arg105)*xi/alpha*2.d0
   dvl12j0=gp*(2.d0*t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
   !
   return
end function dvl12j0

function dvl21j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl21j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, arg001, &
                         arg002, arg003, arg005, arg007, arg008, arg009, &
                         arg011, arg101, arg102, arg104, arg105, arg106, &
                         arg107, arg109, arg110, arg111, arg113, t1, t2, t3
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg001= gp*(z-zp)
   arg002=-gp*(z-zp)
   arg003= gp*(z+zp-2.d0*z1)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg007= aaa/2.d0/tmp**2*xi-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg101=  gp/2.d0/tmp+tmp*(z-zp)
   arg102=  gp/2.d0/tmp-tmp*(z-zp)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg106=  gp/2.d0/tmp-tmp*(z-z1)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   t1= exp(arg005)*gp/alpha
   t2= (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
      -exp_erfc(arg008,arg107)*xi/lambda &
      +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda &
      +exp_erfc(arg005,arg104)*gp/alpha*2.d0
   t3= exp_erfc(arg002,arg102)-exp_erfc(arg002,arg106) &
      -exp_erfc(arg001,arg101)-exp_erfc(arg003,arg105)*kappa/alpha &
      +exp_erfc(arg007,arg110)*gp/alpha*2.d0
   dvl21j0=gp*(2.d0*t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
   !
   return
end function dvl21j0

function dvl22j0(gp,aaa,tmp,z1,z,zp,rxy)
   USE kinds,            ONLY : DP
   implicit none
   real(DP)            :: dvl22j0
   real(DP),intent(in) :: gp, aaa, tmp, z1, z, zp, rxy
   !local
   real(DP)            ::alpha, beta, kappa, xi, chi, lambda, arg000, &
                         arg004, arg005, arg008, arg009, arg010, arg011, &
                         arg012, arg104, arg105, arg107, arg108, arg109, &
                         arg110, arg111, arg112, arg113, arg114, t1, t2, t3
   !
   alpha=aaa+gp+sqrt(aaa**2+gp**2)
   beta =aaa+gp-sqrt(aaa**2+gp**2)
   kappa=aaa-gp+sqrt(aaa**2+gp**2)
   xi   =aaa   +sqrt(aaa**2+gp**2)
   chi  =aaa   -sqrt(aaa**2+gp**2)
   lambda=      sqrt(aaa**2+gp**2)
   arg000=-xi*(z+zp-2.d0*z1)
   arg004= gp*(z-z1)+xi*(z1-zp)
   arg005=-gp*(z1-zp)-xi*(z-z1)
   arg008= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-chi*(z-z1)
   arg009= aaa/2.d0/tmp**2*xi+xi*(z1-zp)-xi*(z-z1)
   arg010= aaa/2.d0/tmp**2*xi+chi*(z1-zp)-xi*(z-z1)
   arg011= aaa/2.d0/tmp**2*chi+chi*(z1-zp)-xi*(z-z1)
   arg012= aaa/2.d0/tmp**2*chi+xi*(z1-zp)-chi*(z-z1)
   arg104=  gp/2.d0/tmp-tmp*(z1-zp)
   arg105=  gp/2.d0/tmp+tmp*(z-z1)
   arg107=  xi/2.d0/tmp+tmp*(z-zp)
   arg108=  xi/2.d0/tmp-tmp*(z-zp)
   arg109=  xi/2.d0/tmp+tmp*(z1-zp)
   arg110=  xi/2.d0/tmp-tmp*(z-z1)
   arg111= chi/2.d0/tmp+tmp*(z-zp)
   arg112= chi/2.d0/tmp-tmp*(z-zp)
   arg113= chi/2.d0/tmp+tmp*(z1-zp)
   arg114= chi/2.d0/tmp-tmp*(z-z1)
   t1= exp(arg000)*beta*xi/alpha/lambda
   t2= (exp_erfc(arg011,arg111)-exp_erfc(arg011,arg113))*chi/lambda &
      -exp_erfc(arg008,arg107)*xi/lambda &
      +exp_erfc(arg009,arg109)*xi*beta/alpha/lambda &
      +exp_erfc(arg005,arg104)*gp/alpha*2.d0
   t3= (exp_erfc(arg012,arg112)-exp_erfc(arg012,arg114))*xi/lambda &
      -exp_erfc(arg010,arg108)*chi/lambda &
      +exp_erfc(arg009,arg110)*xi*beta/alpha/lambda &
      -exp_erfc(arg004,arg105)*xi/alpha*2.d0
   dvl22j0=gp*(t1-(t2+t3)/4.d0)*dbesj0(gp*rxy)
   !
   return
end function dvl22j0

SUBROUTINE qromb(func,aaa,tmp,z1,z,zp,rxy,b,ss)
   USE kinds,            ONLY : DP
   integer, parameter   ::jmax=20, jmaxp=jmax+1, k=5, km=k-1
   integer              ::j
   real(DP),intent(in)  ::aaa, tmp, z1, z, zp, rxy, b
   real(DP),intent(out) ::ss
   real(DP), parameter  ::a=0.0_DP,eps=1.e-12
   real(DP), external   ::func
   real(DP)             ::dss, h(jmaxp), s(jmaxp)
   !
   ! ss=int_a^b func(gp,aaa,tmp,z1,z,zp,rxy) dgp
   !
   h(1)=1.0_DP
   do j=1,jmax
     call trapzd(func,aaa,tmp,z1,z,zp,rxy,a,b,s(j),j)
     if (j.ge.k) then
       call polint(h(j-km),s(j-km),k,0.0_DP,ss,dss)
       if (abs(ss).le.1.e-8 ) return
       if (abs(dss).le.eps*abs(ss)) return
     endif
     s(j+1)=s(j)
     h(j+1)=0.25*h(j)
   enddo
   stop 'too many steps in qromb'
END SUBROUTINE qromb

SUBROUTINE trapzd(func,aaa,tmp,z1,z,zp,rxy,a,b,s,n)
   USE kinds,            ONLY : DP
   integer,intent(in)::n
   real(DP),intent(in)::aaa, tmp, z1, z, zp, rxy, a, b
   real(DP),intent(inout)::s
   real(DP),external ::func
   integer::it,j
   real(DP)::del, sum, tnm, x

   if (n.eq.1) then
     s=0.5*(b-a)*(func(a,aaa,tmp,z1,z,zp,rxy)+func(b,aaa,tmp,z1,z,zp,rxy))
   else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do j=1,it
       sum=sum+func(x,aaa,tmp,z1,z,zp,rxy)
       x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
   endif
   return
END SUBROUTINE trapzd

SUBROUTINE polint(xa,ya,n,x,y,dy)
   USE kinds,            ONLY : DP
   integer, intent(in)  ::n
   integer, parameter   ::nmax=10
   integer              ::i, m, ns
   real(DP), intent(in) ::x, xa(n), ya(n)
   real(DP), intent(out)::y, dy
   real(DP)             ::den, dif, dift, ho, hp, w, c(nmax), d(nmax)

   ns=1
   dif=abs(x-xa(1))
   do i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
       ns=i
       dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
   enddo
   y=ya(ns)
   ns=ns-1
   do m=1,n-1
     do i=1,n-m
       ho=xa(i)-x
       hp=xa(i+m)-x
       w=c(i+1)-d(i)
       den=ho-hp
       if(den.eq.0.)stop 'failure in polint'
       den=w/den
       d(i)=hp*den
       c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
       dy=c(ns+1)
     else
       dy=d(ns)
       ns=ns-1
     endif
     y=y+dy
   enddo
   return
END SUBROUTINE polint
!
! Bessel J_0(x) function in double precision
!
real(8) FUNCTION dbesj0(x)
   implicit none
   real(8), intent(in) :: x
   real(8), parameter  :: pi4 = 0.78539816339744830962d0
   real(8), parameter  :: a(0:7) = (/ &
      -0.0000000000023655394d0, 0.0000000004708898680d0, &
      -0.0000000678167892231d0, 0.0000067816840038636d0, &
      -0.0004340277777716935d0, 0.0156249999999992397d0, &
      -0.2499999999999999638d0, 0.9999999999999999997d0 /)
   real(8), parameter  :: b(0:12, 0:4) = reshape( (/ &
       0.0000000000626681117d0, -0.0000000022270614428d0, & 
       0.0000000662981656302d0, -0.0000016268486502196d0, & 
       0.0000321978384111685d0, -0.0005005237733315830d0, & 
       0.0059060313537449816d0, -0.0505265323740109701d0, & 
       0.2936432097610503985d0, -1.0482565081091638637d0, & 
       1.9181123286040428113d0, -1.1319199475221700100d0, & 
      -0.1965480952704682000d0, &
       0.0000000000457457332d0, -0.0000000015814772025d0, & 
       0.0000000455487446311d0, -0.0000010735201286233d0, & 
       0.0000202015179970014d0, -0.0002942392368203808d0, & 
       0.0031801987726150648d0, -0.0239875209742846362d0, & 
       0.1141447698973777641d0, -0.2766726722823530233d0, & 
       0.1088620480970941648d0,  0.5136514645381999197d0, & 
      -0.2100594022073706033d0, &
       0.0000000000331366618d0, -0.0000000011119090229d0, & 
       0.0000000308823040363d0, -0.0000006956602653104d0, & 
       0.0000123499947481762d0, -0.0001662951945396180d0, & 
       0.0016048663165678412d0, -0.0100785479932760966d0, & 
       0.0328996815223415274d0, -0.0056168761733860688d0, & 
      -0.2341096400274429386d0,  0.2551729256776404262d0, & 
       0.2288438186148935667d0, & 
       0.0000000000238007203d0, -0.0000000007731046439d0, & 
       0.0000000206237001152d0, -0.0000004412291442285d0, & 
       0.0000073107766249655d0, -0.0000891749801028666d0, & 
       0.0007341654513841350d0, -0.0033303085445352071d0, & 
       0.0015425853045205717d0,  0.0521100583113136379d0, & 
      -0.1334447768979217815d0, -0.1401330292364750968d0, & 
       0.2685616168804818919d0, &
       0.0000000000169355950d0, -0.0000000005308092192d0, & 
       0.0000000135323005576d0, -0.0000002726650587978d0, & 
       0.0000041513240141760d0, -0.0000443353052220157d0, & 
       0.0002815740758993879d0, -0.0004393235121629007d0, & 
      -0.0067573531105799347d0,  0.0369141914660130814d0, & 
       0.0081673361942996237d0, -0.2573381285898881860d0, & 
       0.0459580257102978932d0 /), (/13, 5/) )
   real(8), parameter  :: c(0:13, 0:4) = reshape( (/ &
      -0.00000000003009451757d0, -0.00000000014958003844d0, & 
       0.00000000506854544776d0,  0.00000001863564222012d0, & 
      -0.00000060304249068078d0, -0.00000147686259937403d0, & 
       0.00004714331342682714d0,  0.00006286305481740818d0, & 
      -0.00214137170594124344d0, -0.00089157336676889788d0, & 
       0.04508258728666024989d0, -0.00490362805828762224d0, & 
      -0.27312196367405374426d0,  0.04193925184293450356d0,  &
      -0.00000000000712453560d0, -0.00000000041170814825d0, & 
       0.00000000138012624364d0,  0.00000005704447670683d0, & 
      -0.00000019026363528842d0, -0.00000533925032409729d0, & 
       0.00001736064885538091d0,  0.00030692619152608375d0, & 
      -0.00092598938200644367d0, -0.00917934265960017663d0, & 
       0.02287952522866389076d0,  0.10545197546252853195d0, & 
      -0.16126443075752985095d0, -0.19392874768742235538d0,  &
       0.00000000002128344556d0, -0.00000000031053910272d0, & 
      -0.00000000334979293158d0,  0.00000004507232895050d0, & 
       0.00000036437959146427d0, -0.00000446421436266678d0, & 
      -0.00002523429344576552d0,  0.00027519882931758163d0, & 
       0.00097185076358599358d0, -0.00898326746345390692d0, & 
      -0.01665959196063987584d0,  0.11456933464891967814d0, & 
       0.07885001422733148815d0, -0.23664819446234712621d0,  & 
       0.00000000003035295055d0,  0.00000000005486066835d0, & 
      -0.00000000501026824811d0, -0.00000000501246847860d0, & 
       0.00000058012340163034d0,  0.00000016788922416169d0, & 
      -0.00004373270270147275d0,  0.00001183898532719802d0, & 
       0.00189863342862291449d0, -0.00113759249561636130d0, & 
      -0.03846797195329871681d0,  0.02389746880951420335d0, & 
       0.22837862066532347461d0, -0.06765394811166522844d0,  &
       0.00000000001279875977d0,  0.00000000035925958103d0, & 
      -0.00000000228037105967d0, -0.00000004852770517176d0, & 
       0.00000028696428000189d0,  0.00000440131125178642d0, & 
      -0.00002366617753349105d0, -0.00024412456252884129d0, & 
       0.00113028178539430542d0,  0.00708470513919789080d0, & 
      -0.02526914792327618386d0, -0.08006137953480093426d0, & 
       0.16548380461475971846d0,  0.14688405470042110229d0/), (/14, 5/) )
   real(8), parameter  :: d(0:12, 0:3) = reshape( (/ &
       1.059601355592185731d-14, -2.71150591218550377d-13, & 
       8.6514809056201638d-12,   -4.6264028554286627d-10, & 
       5.0815403835647104d-8,    -1.76722552048141208d-5, & 
       0.16286750396763997378d0,  2.949651820598278873d-13, & 
      -8.818215611676125741d-12,  3.571119876162253451d-10, & 
      -2.631924120993717060d-8,   4.709502795656698909d-6, & 
      -5.208333333333283282d-3, & 
       7.18344107717531977d-15,  -2.51623725588410308d-13, & 
       8.6017784918920604d-12,   -4.6256876614290359d-10, & 
       5.0815343220437937d-8,    -1.76722551764941970d-5, & 
       0.16286750396763433767d0,  2.2327570859680094777d-13, & 
      -8.464594853517051292d-12,  3.563766464349055183d-10, & 
      -2.631843986737892965d-8,   4.709502342288659410d-6, & 
      -5.2083333332278466225d-3, & 
       5.15413392842889366d-15,  -2.27740238380640162d-13, & 
       8.4827767197609014d-12,   -4.6224753682737618d-10, & 
       5.0814848128929134d-8,    -1.76722547638767480d-5, & 
       0.16286750396748926663d0,  1.7316195320192170887d-13, & 
      -7.971122772293919646d-12,  3.544039469911895749d-10, & 
      -2.631443902081701081d-8,   4.709498228695400603d-6, & 
      -5.2083333315143653610d-3, & 
       3.84653681453798517d-15,  -2.04464520778789011d-13, & 
       8.3089298605177838d-12,   -4.6155016158412096d-10, & 
       5.0813263696466650d-8,    -1.76722528311426167d-5, & 
       0.16286750396650065930d0,  1.3797879972460878797d-13, & 
      -7.448089381011684812d-12,  3.512733797106959780d-10, & 
      -2.630500895563592722d-8,   4.709483934775839193d-6, & 
      -5.2083333227940760113d-3 /), (/13, 4/) )
   real(8) :: w, t, y, v, theta
   integer :: k, i
   
   w = abs(x)
   if( w < 1.0d0 ) then
      t = w * w
      y = a(0)
      do i = 1, 7
         y = y * t + a(i)
      end do
   else if( w < 8.5d0 ) then
      t = w * w * 0.0625d0
      k = int(t)
      t = t - (k + 0.5d0)
      y = b(0,k)
      do i = 1, 12
         y = y * t + b(i,k)
      end do
   else if( w < 12.5d0 ) then
      k = int(w)
      t = w - ( k + 0.5d0 )
      k = k - 8
      y = c(0,k)
      do i = 1, 13
         y = y * t + c(i,k)
      end do
   else
      v = 24.0d0 / w
      t = v * v
      k = int(t)
      y = d(0,k)
      do i = 1, 6
         y = y * t + d(i,k)
      end do
      y = y * sqrt(v)
      theta = d(7,k)
      do i = 8, 12
         theta = theta * t + d(i,k)
      end do
      theta = theta * v - pi4
      y = y * cos( w + theta )
   end if
   dbesj0 = y
END FUNCTION dbesj0


real(8) FUNCTION dbesj1(x)
   implicit none
   real(8), intent(in) :: x
   real(8), parameter :: pi4 = 0.78539816339744830962d0
   real(8), parameter :: a(0:7) = (/ &
      -0.00000000000014810349d0,  0.00000000003363594618d0, &
      -0.00000000565140051697d0,  0.00000067816840144764d0, &
      -0.00005425347222188379d0,  0.00260416666666662438d0, &
      -0.06249999999999999799d0,  0.49999999999999999998d0 /)
   real(8), parameter :: b(0:12, 0:4) = reshape( (/ &
       0.00000000000243721316d0, -0.00000000009400554763d0, &
       0.00000000306053389980d0, -0.00000008287270492518d0, &
       0.00000183020515991344d0, -0.00003219783841164382d0, &
       0.00043795830161515318d0, -0.00442952351530868999d0, &
       0.03157908273375945955d0, -0.14682160488052520107d0, &
       0.39309619054093640008d0, -0.47952808215101070280d0, &
       0.14148999344027125140d0, &
       0.00000000000182119257d0, -0.00000000006862117678d0, &
       0.00000000217327908360d0, -0.00000005693592917820d0, &
       0.00000120771046483277d0, -0.00002020151799736374d0, &
       0.00025745933218048448d0, -0.00238514907946126334d0, &
       0.01499220060892984289d0, -0.05707238494868888345d0, &
       0.10375225210588234727d0, -0.02721551202427354117d0, &
      -0.06420643306727498985d0, &
       0.000000000001352611196d0, -0.000000000049706947875d0, &
       0.000000001527944986332d0, -0.000000038602878823401d0, &
       0.000000782618036237845d0, -0.000012349994748451100d0, &
       0.000145508295194426686d0, -0.001203649737425854162d0, &
       0.006299092495799005109d0, -0.016449840761170764763d0, &
       0.002106328565019748701d0,  0.058527410006860734650d0, &
      -0.031896615709705053191d0, &
       0.000000000000997982124d0, -0.000000000035702556073d0, &
       0.000000001062332772617d0, -0.000000025779624221725d0, &
       0.000000496382962683556d0, -0.000007310776625173004d0, &
       0.000078028107569541842d0, -0.000550624088538081113d0, &
       0.002081442840335570371d0, -0.000771292652260286633d0, &
      -0.019541271866742634199d0,  0.033361194224480445382d0, &
       0.017516628654559387164d0, &
       0.000000000000731050661d0, -0.000000000025404499912d0, &
       0.000000000729360079088d0, -0.000000016915375004937d0, &
       0.000000306748319652546d0, -0.000004151324014331739d0, &
       0.000038793392054271497d0, -0.000211180556924525773d0, &
       0.000274577195102593786d0,  0.003378676555289966782d0, &
      -0.013842821799754920148d0, -0.002041834048574905921d0, &
       0.032167266073736023299d0 /), (/13, 5/) )
   real(8), parameter :: c(0:13, 0:4) = reshape( (/ &
       -0.00000000001185964494d0,  0.00000000039110295657d0, &
        0.00000000180385519493d0, -0.00000005575391345723d0, &
       -0.00000018635897017174d0,  0.00000542738239401869d0, &
        0.00001181490114244279d0, -0.00033000319398521070d0, &
       -0.00037717832892725053d0,  0.01070685852970608288d0, &
        0.00356629346707622489d0, -0.13524776185998074716d0, &
        0.00980725611657523952d0,  0.27312196367405374425d0,  &
       -0.00000000003029591097d0,  0.00000000009259293559d0, &
        0.00000000496321971223d0, -0.00000001518137078639d0, &
       -0.00000057045127595547d0,  0.00000171237271302072d0, &
        0.00004271400348035384d0, -0.00012152454198713258d0, &
       -0.00184155714921474963d0,  0.00462994691003219055d0, &
        0.03671737063840232452d0, -0.06863857568599167175d0, &
       -0.21090395092505707655d0,  0.16126443075752985095d0,  & 
       -0.00000000002197602080d0, -0.00000000027659100729d0, &
        0.00000000374295124827d0,  0.00000003684765777023d0, &
       -0.00000045072801091574d0, -0.00000327941630669276d0, &
        0.00003571371554516300d0,  0.00017664005411843533d0, &
       -0.00165119297594774104d0, -0.00485925381792986774d0, &
        0.03593306985381680131d0,  0.04997877588191962563d0, &
       -0.22913866929783936544d0, -0.07885001422733148814d0,  &
        0.00000000000516292316d0, -0.00000000039445956763d0, &
       -0.00000000066220021263d0,  0.00000005511286218639d0, &
        0.00000005012579400780d0, -0.00000522111059203425d0, &
       -0.00000134311394455105d0,  0.00030612891890766805d0, &
       -0.00007103391195326182d0, -0.00949316714311443491d0, &
        0.00455036998246516948d0,  0.11540391585989614784d0, &
       -0.04779493761902840455d0, -0.22837862066532347460d0,  &
        0.00000000002697817493d0, -0.00000000016633326949d0, &
       -0.00000000433134860350d0,  0.00000002508404686362d0, &
        0.00000048528284780984d0, -0.00000258267851112118d0, &
       -0.00003521049080466759d0,  0.00016566324273339952d0, &
        0.00146474737522491617d0, -0.00565140892697147306d0, &
       -0.02833882055679300400d0,  0.07580744376982855057d0, &
        0.16012275906960187978d0, -0.16548380461475971845d0 /), (/14, 5/) )
   real(8), parameter :: d(0:12, 0:3) = reshape( (/ &
       -1.272346002224188092d-14,  3.370464692346669075d-13, &
       -1.144940314335484869d-11,  6.863141561083429745d-10, &
       -9.491933932960924159d-8,   5.301676561445687562d-5,  &
        0.1628675039676399740d0,  -3.652982212914147794d-13, &
        1.151126750560028914d-11, -5.165585095674343486d-10, &
        4.657991250060549892d-8,  -1.186794704692706504d-5,  &
        1.562499999999994026d-2, &
       -8.713069680903981555d-15,  3.140780373478474935d-13, &
       -1.139089186076256597d-11,  6.862299023338785566d-10, &
       -9.491926788274594674d-8,   5.301676558106268323d-5,  &
        0.1628675039676466220d0,  -2.792555727162752006d-13, &
        1.108650207651756807d-11, -5.156745588549830981d-10, &
        4.657894859077370979d-8,  -1.186794650130550256d-5,  &
        1.562499999987299901d-2, & 
       -6.304859171204770696d-15,  2.857249044208791652d-13, &
       -1.124956921556753188d-11,  6.858482894906716661d-10, &
       -9.491867953516898460d-8,   5.301676509057781574d-5,  &
        0.1628675039678191167d0,  -2.185193490132496053d-13, &
        1.048820673697426074d-11, -5.132819367467680132d-10, &
        4.657409437372994220d-8,  -1.186794150862988921d-5,  &
        1.562499999779270706d-2, &
       -4.740417209792009850d-15,  2.578715253644144182d-13, &
       -1.104148898414138857d-11,  6.850134201626289183d-10, &
       -9.491678234174919640d-8,   5.301676277588728159d-5,  &
        0.1628675039690033136d0,  -1.755122057493842290d-13, &
        9.848723331445182397d-12, -5.094535425482245697d-10, &
        4.656255982268609304d-8,  -1.186792402114394891d-5,  &
        1.562499998712198636d-2 /), (/13, 4/) )
   real(8) :: w, t, y, v, theta
   integer :: k, i

   w = abs(x)
   if( w < 1.0d0 ) then
      t = w * w
      y = a(0)
      do i = 1, 7
         y = y * t + a(i)
      end do
      y = y * w
   else if( w < 8.5d0 ) then
      t = w * w * 0.0625d0
      k = int(t)
      t = t - (k + 0.5d0)
      y = b(0,k)
      do i = 1, 12
         y = y * t + b(i,k)
      end do
      y = y * w
   else if( w < 12.5d0 ) then
      k = int(w)
      t = w - (k + 0.5d0)
      k = k - 8
      y = c(0,k)
      do i = 1, 13
         y = y * t + c(i,k)
      end do
   else
      v = 24.0d0 / w
      t = v * v
      k = int(t)
      y = d(0,k)
      do i = 1, 6
         y = y * t + d(i,k)
      end do
      y = y * sqrt(v)
      theta = d(7,k)
      do i = 8, 12
         theta = theta * t + d(i,k)
      end do
      theta = theta * v - pi4
      y = y * sin(w + theta)
   end if
   if( x < 0.0d0 ) y = -y
   dbesj1 = y
END FUNCTION dbesj1

! exp(x) * erfc(y)
! This function is to avoid INFINITY * ZERO for large positive x and y.
real(8) function exp_erfc (x, y)
   implicit none
   real(8), intent(in) :: x, y
   real(8)             :: ym, ym2, nume, deno
   real(8), parameter  :: rtpim = 0.564189583547756279d0 ! 1 / sqrt(PI)
   real(8), parameter  :: r(0:4) = (/ &
      -2.99610707703542174d-3, -4.94730910623250734d-2, &
      -2.26956593539686930d-1, -2.78661308609647788d-1, &
      -2.23192459734184686d-2 /)
   real(8), parameter  :: s(0:4) = (/ &
      1.06209230528467918d-2, 1.91308926107829841d-1, &
      1.05167510706793207d0,  1.98733201817135256d0,  &
      1.00000000000000000d0 /)
   
   if( x < 709.0d0 .or. y < 4.0d0 ) then
      exp_erfc = exp(x) * qe_erfc(y)
   else
      ym  = 1d0 / y
      ym2 = ym ** 2
      nume =( ( ( r(4) * ym2 + r(3) ) * ym2 + r(2) ) * ym2 + r(1) ) * ym2 + r(0)
      deno =( ( ( s(4) * ym2 + s(3) ) * ym2 + s(2) ) * ym2 + s(1) ) * ym2 + s(0)
      exp_erfc = exp( - y**2 + x) * ym * ( rtpim + ym2 * nume / deno )
   end if
   
   return
end function exp_erfc

END MODULE esm

