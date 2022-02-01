!
! Copyright (C) 2002-2008 Quantm-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!----------------------------------------------------------------
#include<cpv_device_macros.h> 
!----------------------------------------------------------------



    SUBROUTINE potential_print_info( iunit )

        USE control_flags, ONLY: iesr

        INTEGER, INTENT(IN) :: iunit

        WRITE(iunit,50)
        WRITE(iunit,115) (2*iesr+1),(2*iesr+1),(2*iesr+1)

   50   FORMAT(//,3X,'Potentials Parameters',/,3X,'---------------------')
  115   FORMAT(   3X,'Ewald sum over ',I1,'*',I1,'*',I1,' cells')

        RETURN
   END SUBROUTINE potential_print_info

!=----------------------------------------------------------------------------=!

  SUBROUTINE cluster_bc( screen_coul, hg, omega, hmat )

      USE kinds,           ONLY: DP
      USE mp_global,       ONLY: me_bgrp
      USE fft_base,        ONLY: dfftp
      USE fft_interfaces,  ONLY: fwfft
      USE constants,       ONLY: gsmall, pi
      USE cell_base,       ONLY: tpiba2, s_to_r, alat
      USE fft_rho

      IMPLICIT NONE
      
      REAL(DP), INTENT(IN) :: hg( dfftp%ngm )
      REAL(DP), INTENT(IN) :: omega, hmat( 3, 3 )
      COMPLEX(DP) :: screen_coul( dfftp%ngm )

      ! ... Locals
      !
      REAL(DP), ALLOCATABLE :: grr(:,:)
      COMPLEX(DP), ALLOCATABLE :: grg(:,:)
      REAL(DP) :: rc, r(3), s(3), rmod, g2, rc2, arg, fact
      INTEGER   :: ig, i, j, k, ir
      INTEGER   :: ir1, ir2, ir3, nr3l

      ir1 = 1
      ir2 = 1
      ir3 = 1
      DO k = 1, me_bgrp
        ir3 = ir3 + dfftp%nr3p( k )
      END DO
      nr3l = dfftp%my_nr3p

      ALLOCATE( grr( dfftp%nnr, 1 ) )
      ALLOCATE( grg( dfftp%nnr, 1 ) )

      grr = 0.0d0

      ! ... Martyna and Tuckerman convergence criterium
      !
      rc  = 7.0d0 / alat
      rc2 = rc**2
      fact  = omega / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
      IF( MOD(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, 2) /= 0 ) fact = -fact

      DO k = 1, nr3l
        s(3) = DBLE ( (k-1) + (ir3 - 1) ) / dfftp%nr3 - 0.5d0
        DO j = 1, dfftp%nr2
          s(2) = DBLE ( (j-1) + (ir2 - 1) ) / dfftp%nr2 - 0.5d0
          DO i = 1, dfftp%nr1
            s(1) = DBLE ( (i-1) + (ir1 - 1) ) / dfftp%nr1 - 0.5d0
            CALL S_TO_R( S, R, hmat )
            rmod = SQRT( r(1)**2 + r(2)**2 + r(3)**2 )
            ir =  i + (j-1)*dfftp%nr1x + (k-1)*dfftp%nr1x*dfftp%nr2x
            IF( rmod < gsmall ) THEN
              grr( ir, 1 ) = fact * 2.0d0 * rc / SQRT( pi )
            ELSE
              grr( ir, 1 ) = fact * erf( rc * rmod ) / rmod
            END IF
          END DO
        END DO
      END DO

      ! grg = FFT( grr )

      CALL rho_r2g( dfftp, grr, grg )

      DO ig = 1, SIZE( screen_coul )
        IF( hg(ig) < gsmall ) THEN
          screen_coul(ig) = grg(1,1) - ( - pi / rc2 )
        ELSE
          g2  = tpiba2 * hg(ig)
          arg = - g2 / ( 4.0d0 * rc2 )
          screen_coul(ig) = grg(ig,1) - ( 4.0d0 * pi * EXP( arg ) / g2 ) 
        END IF
      END DO

      DEALLOCATE( grr, grg )

    RETURN
  END SUBROUTINE cluster_bc


!=----------------------------------------------------------------------------=!


   SUBROUTINE vofps_x( eps, vloc, rhoeg, vps, sfac, omega )

      !  this routine computes:
      !  omega = ht%deth
      !  vloc_ps(ig)  =  (sum over is) sfac(is,ig) * vps(ig,is)
      !
      !  Eps = Fact * omega * (sum over ig) cmplx ( rho_e(ig) ) * vloc_ps(ig)
      !  if Gamma symmetry Fact = 2 else Fact = 1
      !

      USE kinds,              ONLY: DP
      USE io_global,          ONLY: stdout
      USE ions_base,          ONLY: nsp
      USE gvect,              ONLY: gstart
      USE mp_global,          ONLY: intra_bgrp_comm
      USE mp,                 ONLY: mp_sum
      USE fft_base,           ONLY: dfftp

      IMPLICIT NONE

      ! ... Arguments

      REAL(DP),    INTENT(IN)  :: vps(:,:)
      REAL(DP),    INTENT(IN)  :: omega
      COMPLEX(DP), INTENT(OUT) :: vloc(:)
      COMPLEX(DP), INTENT(IN)  :: rhoeg(:)
      COMPLEX(DP), INTENT(IN)  :: sfac(:,:)
      COMPLEX(DP), INTENT(OUT) :: eps

      ! ... Locals

      INTEGER     :: is, ig
      COMPLEX(DP) :: vp

      ! ... Subroutine body ...
      !
      eps   = (0.D0,0.D0)
      !
      DO ig = gstart, dfftp%ngm 

        vp   = (0.D0,0.D0)
        DO is = 1, nsp
          vp = vp + sfac( ig, is ) * vps( ig, is )
        END DO

        vloc(ig) = vp
        eps      = eps  +     vp * CONJG( rhoeg( ig ) )

      END DO
      ! ... 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        vp = (0.D0,0.D0)
        DO is = 1, nsp
          vp = vp + sfac( 1, is) * vps(1, is)
        END DO
        vloc(1) = VP
        eps     = eps + vp * CONJG( rhoeg(1) ) * 0.5d0
      END IF
      !
      eps = 2.D0 * eps  * omega
      !
      CALL mp_sum( eps, intra_bgrp_comm )

      RETURN
   END SUBROUTINE vofps_x


!=----------------------------------------------------------------------------=!

  SUBROUTINE vofloc_x( tscreen, ehte, ehti, eh, vloc, rhoeg, &
                     rhops, vps, sfac, omega, screen_coul )

      !  this routine computes:
      !  omega = ht%deth
      !  rho_e(ig)    =  (sum over iss) rhoeg(ig,iss) 
      !  rho_I(ig)    =  (sum over is) sfac(is,ig) * rhops(ig,is) 
      !  vloc_h(ig)   =  fpi / ( g(ig) * tpiba2 ) * { rho_e(ig) + rho_I(ig) }
      !
      !  Eh  = Fact * omega * (sum over ig) * fpi / ( g(ig) * tpiba2 ) *
      !        { rho_e(ig) + rho_I(ig) } * conjugate { rho_e(ig) + rho_I(ig) }
      !  if Gamma symmetry Fact = 1 else Fact = 1/2
      !
      !  Hatree potential and local pseudopotential
      !  vloc(ig)     =  vloc_h(ig) + vloc_ps(ig) 
      !

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2, tpiba
      USE io_global,          ONLY: stdout
      USE gvect,              ONLY: gstart, gg
      USE ions_base,          ONLY: nsp
      USE mp_global,          ONLY: intra_bgrp_comm
      USE mp,                 ONLY: mp_sum
      USE fft_base,           ONLY: dfftp

      IMPLICIT NONE

      ! ... Arguments

      LOGICAL,     INTENT(IN)    :: tscreen
      REAL(DP),    INTENT(IN)    :: rhops(:,:), vps(:,:)
      COMPLEX(DP), INTENT(INOUT) :: vloc(:)
      COMPLEX(DP), INTENT(IN)    :: rhoeg(:)
      COMPLEX(DP), INTENT(IN)    :: sfac(:,:)
      REAL(DP),    INTENT(OUT)   :: ehte, ehti
      REAL(DP),    INTENT(IN)    :: omega
      COMPLEX(DP), INTENT(OUT)   :: eh
      COMPLEX(DP), INTENT(IN)    :: screen_coul(:)

      ! ... Locals

      INTEGER     :: is, ig
      REAL(DP)    :: fpibg, cost
      COMPLEX(DP) :: rhet, rhog, rp, vscreen

      ! ... Subroutine body ...

      eh    = 0.0d0
      ehte  = 0.0d0
      ehti  = 0.0d0

!$omp parallel do default(shared), private(rp,is,rhet,rhog,fpibg), reduction(+:eh,ehte,ehti)
      DO ig = gstart, dfftp%ngm 

        rp   = (0.D0,0.D0)
        DO is = 1, nsp
          rp = rp + sfac( ig, is ) * rhops( ig, is )
        END DO

        rhet  = rhoeg( ig )
        rhog  = rhet + rp

        IF( tscreen ) THEN
          fpibg     = fpi / ( gg(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          fpibg     = fpi / ( gg(ig) * tpiba2 )
        END IF

        vloc(ig) = vloc(ig)  +  fpibg *        rhog 
        eh       = eh        +  fpibg *        rhog * CONJG(rhog)
        ehte     = ehte      +  fpibg *   DBLE(rhet * CONJG(rhet))
        ehti     = ehti      +  fpibg *   DBLE(  rp * CONJG(rp))

      END DO
      ! ... 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        rp = (0.D0,0.D0)
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        DO IS = 1, nsp
          rp = rp + sfac( 1, is) * rhops(1, is)
        END DO
        rhet    = rhoeg(1)
        rhog    = rhet + rp
        vloc(1) = vloc(1)   +  vscreen *   rhog
        eh      = eh        +  vscreen *        rhog * CONJG(rhog)
        ehte    = ehte      +  vscreen *   DBLE(rhet * CONJG(rhet))
        ehti    = ehti      +  vscreen *   DBLE(  rp * CONJG(rp))
      END IF
      ! ...
      eh   =        eh   * omega
      ehte =        ehte * omega
      ehti =        ehti * omega
      ! ...
      CALL mp_sum(eh  , intra_bgrp_comm)
      CALL mp_sum(ehte, intra_bgrp_comm)
      CALL mp_sum(ehti, intra_bgrp_comm)
      !
      RETURN
  END SUBROUTINE vofloc_x


  SUBROUTINE force_loc_x( tscreen, rhoeg, fion, rhops, vps, ei1, ei2, ei3, &
                        sfac, omega, screen_coul )

      !  this routine computes:
      !
      !  Local contribution to the forces on the ions
      !  eigrx(ig,isa)   = ei1( mill(1,ig), isa)
      !  eigry(ig,isa)   = ei2( mill(2,ig), isa)
      !  eigrz(ig,isa)   = ei3( mill(3,ig), isa)
      !  fpibg           = fpi / ( g(ig) * tpiba2 )
      !  tx_h(ig,is)     = fpibg * rhops(ig, is) * CONJG( rho_e(ig) + rho_I(ig) )
      !  tx_ps(ig,is)    = vps(ig,is) * CONJG( rho_e(ig) )
      !  gx(ig)          = cmplx(0.D0, gx(1,ig),kind=DP) * tpiba
      !  fion(x,isa)     = fion(x,isa) + 
      !      Fact * omega * ( sum over ig, iss) (tx_h(ig,is) + tx_ps(ig,is)) * 
      !      gx(ig) * eigrx(ig,isa) * eigry(ig,isa) * eigrz(ig,isa) 
      !  if Gamma symmetry Fact = 2.0 else Fact = 1
      !

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2, tpiba
      USE io_global,          ONLY: stdout
      USE gvect,              ONLY: mill, gstart, g, gg
      USE ions_base,          ONLY: nat, nsp, ityp
      USE fft_base,           ONLY: dfftp, dffts

      IMPLICIT NONE

      ! ... Arguments

      LOGICAL     :: tscreen
      REAL(DP)    :: fion(:,:)
      REAL(DP)    :: rhops(:,:), vps(:,:)
      COMPLEX(DP) :: rhoeg(:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      COMPLEX(DP) :: ei1(-dfftp%nr1:dfftp%nr1,nat)
      COMPLEX(DP) :: ei2(-dfftp%nr2:dfftp%nr2,nat)
      COMPLEX(DP) :: ei3(-dfftp%nr3:dfftp%nr3,nat)
      REAL(DP)    :: omega
      COMPLEX(DP) :: screen_coul(:)

      ! ... Locals

      INTEGER     :: is, ia, ig, ig1, ig2, ig3
      REAL(DP)    :: fpibg
      COMPLEX(DP) :: cxc, rhet, rhog, vp, gxc, gyc, gzc
      COMPLEX(DP),ALLOCATABLE :: rp(:) 
      INTEGER                 :: ngm_ 
      ! 
      COMPLEX(DP) :: teigr, cnvg, cvn, tx, ty, tz, fx, fy, fz
      COMPLEX(DP), ALLOCATABLE :: ftmp(:,:)
      INTEGER :: s_ngm_ 

      ! ... Subroutine body ... 

      s_ngm_ = dffts%ngm
      ALLOCATE (rp(s_ngm_)) 
  
!

DEV_ACC data present(rhoeg, rhops, mill,g ) copy(fion)  create(rp(1:s_ngm_)) copyin(sfac, screen_coul, gg, vps, ityp,ei1, ei2, ei3) 
!
DEV_OMP parallel default(none) &
DEV_OMP shared(gstart, dffts,sfac, rhops, screen_coul, rhoeg, nsp, gg, tpiba2, tpiba, mill, g, &
DEV_OMP         nat, ityp, vps, ei1, ei2, ei3, tscreen, rp, fion, omega, s_ngm_ ) &
DEV_OMP private(ig, is, rhet, rhog, fpibg, ig1, ig2, ig3, gxc, gyc, gzc, ia, cnvg, cvn, tx, &
DEV_OMP          ty, tz, teigr,fx, fy, fz )
 
 
DEV_ACC parallel loop       
DEV_OMP do  
   DO ig = gstart, s_ngm_
      rp( ig) = (0.d0,0.d0) 
      DO is = 1, nsp 
         rp( ig)   = rp( ig)  + sfac(ig, is) * rhops( ig, is)
      END DO 
   END DO 

DEV_ACC parallel vector_length(128) 
DEV_ACC loop gang private(is, fx,fy,fz) 
DEV_OMP do  
   DO ia = 1, nat
      is = ityp(ia) 
      fx = (0.d0, 0.d0) 
      fy = (0.d0, 0.d0) 
      fz = (0.d0, 0.d0) 
DEV_ACC loop vector private(rhet, rhog, fpibg, ig1, ig2, ig3, gxc,gyc,gzc, cnvg, cvn, &
DEV_ACC&                             tx, ty, tz, teigr) reduction(+:fx,fy,fz)   
      DO ig = gstart, s_ngm_ 
         rhet = rhoeg ( ig ) 
         rhog = rhet + rp ( ig)
         IF ( tscreen ) THEN 
            fpibg = fpi / ( gg(ig) * tpiba2  ) + screen_coul (ig) 
         ELSE  
            fpibg = fpi / ( gg (ig) * tpiba2 ) 
         END IF 
         ig1 = mill (1,ig)
         ig2 = mill (2,ig) 
         ig3 = mill (3,ig) 
         gxc = CMPLX(0.d0,g(1,ig),KIND=DP) 
         gyc = CMPLX(0.d0,g(2,ig),KIND=DP) 
         gzc = CMPLX(0.d0,g(3,ig),KIND=DP)
         cnvg = rhops ( ig, is) * fpibg * CONJG ( rhog ) 
         cvn  = vps ( ig, is ) * CONJG( rhet) 
         tx   = (cnvg + cvn) * gxc 
         ty   = (cnvg + cvn) * gyc 
         tz   = (cnvg + cvn) * gzc 
         teigr = ei1( ig1, ia) * ei2 ( ig2, ia) * ei3 (ig3, ia) 
         fx    = fx + teigr * tx 
         fy    = fy + teigr * ty 
         fz    = fz + teigr * tz
      END DO         
      fion (:,ia) =  fion(:,ia) + [DBLE(fx),DBLE(fy),DBLE(fz)] * 2.d0 * omega * tpiba
   END DO
DEV_ACC end parallel
DEV_ACC end data 
DEV_OMP end parallel 
   DEALLOCATE (rp) 
      RETURN
      END SUBROUTINE force_loc_x


!
!=----------------------------------------------------------------------------=!
   SUBROUTINE vofesr( iesr, esr, desr, fion, taus, tstress, hmat )
!=----------------------------------------------------------------------------=!

      USE kinds,       ONLY : DP
      USE constants,   ONLY : sqrtpm1
      USE cell_base,   ONLY : s_to_r, pbcs
      USE mp_global,   ONLY : nproc_bgrp, me_bgrp, intra_bgrp_comm
      USE mp,          ONLY : mp_sum
      USE ions_base,   ONLY : rcmax, zv, nsp, na, nat, ityp
 
      IMPLICIT NONE

! ... ARGUMENTS 
      
      INTEGER,  INTENT(IN) :: iesr
      REAL(DP), INTENT(IN) :: taus(3,nat)
      REAL(DP) :: ESR
      REAL(DP) :: DESR(6)
      REAL(DP) :: FION(3,nat)
      LOGICAL,  INTENT(IN) :: TSTRESS
      REAL(DP), INTENT(in) :: hmat( 3, 3 )

      INTEGER, EXTERNAL :: ldim_block, gind_block

      
! ... LOCALS 

      INTEGER :: na_loc, ia_s, ia_e, igis
      INTEGER :: k, i, j, is, ia, ib, ix, iy, iz
      LOGICAL :: split, tzero, tshift
      REAL(DP), ALLOCATABLE :: zv2(:,:)
      REAL(DP), ALLOCATABLE :: rc(:,:)  
      REAL(DP), ALLOCATABLE :: fionloc(:,:) 
      REAL(DP)  :: rxlm(3), sxlm(3)
      REAL(DP)  :: xlm,ylm,zlm, xlm0,ylm0,zlm0, erre2, rlm, arg, esrtzero
      REAL(DP)  :: addesr, addpre, repand, fxx
      REAL(DP)  :: rckj_m1
      REAL(DP)  :: zvk, zvj, zv2_kj
      REAL(DP)  :: fact_pre

      INTEGER, DIMENSION(6), PARAMETER :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: BETA  = (/ 1,1,1,2,2,3 /)

      INTEGER :: omp_get_num_threads

! ... SUBROUTINE BODY 

      ALLOCATE( rc( nsp, nsp ) )
      ALLOCATE( zv2( nsp, nsp ) )
      ALLOCATE( fionloc( 3, nat ) )
      rc      = 0.0_DP
      zv2     = 0.0_DP
      fionloc = 0.0_DP

      !  Here pre-compute some factors

      DO j = 1, nsp
        DO k = 1, nsp
          zv2( k, j ) = zv( k ) * zv( j )
          rc ( k, j ) = SQRT( rcmax(k)**2 + rcmax(j)**2 )
        END DO
      END DO

      xlm     = 1.0_DP
      ylm     = 1.0_DP
      zlm     = 1.0_DP
      ESR     = 0.0_DP
      DESR    = 0.0_DP

      !  Distribute the atoms pairs to processors

      NA_LOC = ldim_block( nat, nproc_bgrp, me_bgrp)
      IA_S   = gind_block( 1, nat, nproc_bgrp, me_bgrp )
      IA_E   = IA_S + NA_LOC - 1

!$omp parallel do reduction(+:esr,desr) num_threads(min(max(1,na_loc),omp_get_num_threads())) default(none) &
!$omp private(ia,ib,k,j,zv2_kj,rckj_m1,fact_pre,xlm,ylm,zlm,tzero,xlm0,ylm0,zlm0,ix,iy,iz,sxlm,tshift, &
!$omp         rxlm,erre2,rlm,arg,esrtzero,addesr,addpre,repand,i,fxx ) &
!$omp shared(ia_s,ia_e,nat,ityp,zv2,rc,taus,iesr,hmat,fionloc,tstress,na_loc)
      DO ia = ia_s, ia_e
        DO ib = ia, nat
          k = ityp(ia)
          j = ityp(ib)

          zv2_kj   = zv2(k,j)
          rckj_m1  = 1.0_DP / rc(k,j)
          fact_pre = (2.0_DP * zv2_kj * sqrtpm1) * rckj_m1

          IF( ia.EQ.ib ) THEN      
            ! ...     same atoms
            xlm=0.0_DP; ylm=0.0_DP; zlm=0.0_DP; 
            tzero=.TRUE.
          ELSE
            ! ...     different atoms
            xlm0= taus(1,ia) - taus(1,ib)
            ylm0= taus(2,ia) - taus(2,ib)
            zlm0= taus(3,ia) - taus(3,ib)
            CALL pbcs(xlm0,ylm0,zlm0,xlm,ylm,zlm,1)
            TZERO=.FALSE.
          END IF

          DO IX=-IESR,IESR
            sxlm(1) = XLM + DBLE(IX)
            DO IY=-IESR,IESR
              sxlm(2) = YLM + DBLE(IY)
              DO IZ=-IESR,IESR
                TSHIFT= IX.EQ.0 .AND. IY.EQ.0 .AND. IZ.EQ.0
                IF( .NOT. ( TZERO .AND. TSHIFT ) ) THEN
                  sxlm(3) = ZLM + DBLE(IZ)
                  CALL S_TO_R( sxlm, rxlm, hmat )
                  ERRE2 = rxlm(1)**2 + rxlm(2)**2 + rxlm(3)**2
                  RLM   = SQRT(ERRE2)
                  ARG   = RLM * rckj_m1
                  IF (TZERO) THEN
                     ESRTZERO=0.5D0
                  ELSE
                     ESRTZERO=1.D0
                  END IF
                  ADDESR = ZV2_KJ * erfc(ARG) / RLM
                  ESR    = ESR + ESRTZERO*ADDESR
                  ADDPRE = FACT_PRE * EXP(-ARG*ARG)
                  REPAND = ESRTZERO*(ADDESR + ADDPRE)/ERRE2
                  !
                  DO i = 1, 3
                     fxx = repand * rxlm( i )
                     fionloc( i, ia ) = fionloc( i, ia ) + fxx
                     fionloc( i, ib ) = fionloc( i, ib ) - fxx
                  END DO
                  !
                  IF( tstress ) THEN
                     DO i = 1, 6
                        fxx = repand * rxlm( alpha( i ) ) * rxlm( beta( i ) )
                        desr( i ) = desr( i ) - fxx
                     END DO
                  END IF
                END IF
              END DO    ! IZ
            END DO      ! IY
          END DO        ! IX
        END DO
      END DO
!$omp end parallel do

!
!     each processor add its own contribution to the array FION
!
      !  FION = FION+FIONLOC
      !
      CALL daxpy( 3*nat, 1.0d0, fionloc, 1, fion, 1 )

      CALL mp_sum(esr, intra_bgrp_comm)
     
      DEALLOCATE(rc)
      DEALLOCATE(zv2)
      DEALLOCATE(fionloc)
      
      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE vofesr
!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!
   SUBROUTINE self_vofhar_x( tscreen, self_ehte, vloc, rhoeg, omega, hmat )
!=----------------------------------------------------------------------------=!

      !  adds the hartree part of the self interaction

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE control_flags,      ONLY: gamma_only
      USE cell_base,          ONLY: tpiba2
      USE gvect,              ONLY: gstart, gg
      USE sic_module,         ONLY: sic_epsilon, sic_alpha
      USE mp_global,          ONLY: intra_bgrp_comm
      USE mp,                 ONLY: mp_sum
      USE fft_base,           ONLY: dfftp

      IMPLICIT NONE

      ! ... Arguments
      LOGICAL     :: tscreen
      COMPLEX(DP) :: vloc(:)
      COMPLEX(DP) :: rhoeg(:,:)
      REAL(DP)    :: self_ehte
      REAL(DP), INTENT(IN) :: omega
      REAL(DP), INTENT(IN) :: hmat( 3, 3 )

      ! ... Locals

      INTEGER      :: ig
      REAL(DP)    :: fpibg
      COMPLEX(DP) :: rhog
      COMPLEX(DP) :: ehte
      COMPLEX(DP) :: vscreen
      COMPLEX(DP), ALLOCATABLE :: screen_coul(:)

      ! ... Subroutine body ...

      IF( tscreen ) THEN
        ALLOCATE( screen_coul( dfftp%ngm ) )
        CALL cluster_bc( screen_coul, gg, omega, hmat )
      END IF

      !==  HARTREE ==

      ehte = 0.D0

      DO IG = gstart, dfftp%ngm

        rhog  = rhoeg(ig,1) - rhoeg(ig,2)

        IF( tscreen ) THEN
          FPIBG     = fpi / ( gg(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          FPIBG     = fpi / ( gg(ig) * tpiba2 )
        END IF

        vloc(ig) = fpibg * rhog
        ehte     = ehte   +  fpibg *   rhog * CONJG(rhog)

      END DO
 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        rhog    = rhoeg(1,1) - rhoeg(1,2)
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        vloc(1) = vscreen * rhog
        ehte    = ehte   +  vscreen *  rhog * CONJG(rhog)
      END IF

      ! ...

      IF( .NOT. gamma_only ) THEN
        ehte  = ehte  * 0.5d0
      END IF
      !
      self_ehte = DBLE(ehte) * omega * sic_epsilon
      vloc = vloc * sic_epsilon

      CALL mp_sum( self_ehte, intra_bgrp_comm )

      IF( ALLOCATED( screen_coul ) ) DEALLOCATE( screen_coul )

      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE self_vofhar_x
!=----------------------------------------------------------------------------=!
