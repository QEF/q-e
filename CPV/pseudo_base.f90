!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sat Nov  6 12:44:29 MET 1999
!  ----------------------------------------------
!  BEGIN manual

#undef __BESSEL_TEST

      MODULE pseudo_base

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE nlin_base(ap,hg,wnl)
!  SUBROUTINE nlin_stress_base(ap,hg,wnla)
!  SUBROUTINE nlset_base(ap,wsgset)
!  SUBROUTINE formfn_base(ap,ac,hg,vps,dvps,rhoc1,rhocp,omega)
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE kinds
      USE bessel_functions, ONLY: bessel1, bessel2, bessel3
      USE constants, ONLY: gsmall, fpi, pi 
      USE cell_base, ONLY: tpiba

      IMPLICIT NONE

      SAVE

      PRIVATE

      PUBLIC :: nlin_base, nlin_stress_base, nlset_base, formfn_base
      PUBLIC :: corecor_base, rhops_base


!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE nlin_base(ap, hg, wnl)

        USE pseudo_types, ONLY: pseudo_ncpp
        USE io_global, ONLY: stdout

!  compute Kleinman-Bylander factors
!
!    wnl(ig,l) = Y_l(ig) (integral) j_l(ig,r) V(l,r) r**2 dr
!
!    Y_l(ig) = spherical harmonics at point (k+G)/|k+G|
!    j_l(ig,r) = spherical Bessel function of order l at point (k+G)r
!    V(l,r) = nonlocal pseudopotential for angular momentum l
!
!    l = index of angular momentum
!    ig = index of G vector
!
!    Remember: 
!    i)   j_0( 0 ) = 1 


!  (describe briefly what this routine does...)
!  This subroutine is executed only for non-local potentials
!  ----------------------------------------------

! ...   declare subroutine arguments
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl), INTENT(IN) :: hg(:)
        REAL(dbl), INTENT(OUT) :: wnl(:,:)

! ...   declare other variables
        REAL(dbl), ALLOCATABLE :: fint(:,:)
        REAL(dbl)  :: xg, dx
        INTEGER :: ig, mmax, lnl, l, ind
        REAL(dbl), ALLOCATABLE :: ftest(:,:)

! ...   end of declarations
!  ----------------------------------------------

        lnl  = ap%lnl
        mmax = ap%mesh
        dx = ap%dx

        ALLOCATE(fint(mmax,lnl))

        ! WRITE( stdout,*) 'DEBUG nlin_base, tpiba = ', tpiba

        DO ig = 1, SIZE( wnl, 1 )
          IF( hg(ig) < gsmall ) THEN
            DO l = 1,lnl
! ...         G=0 (Only if l=1, since otherwise the radial Bessel function jl=0)
              IF( ap%indl(l) == 1 ) THEN
                fint(1:mmax,l) = ap%rw(1:mmax)**2 * ap%vrps(1:mmax,l)
                call simpson_fpmd(mmax, fint(:,l), dx, wnl(ig,l))
                ! call simpson(mmax, fint(:,l), ap%rab, wnl(ig,l))
              ELSE
                wnl(ig,l) = 0.d0
              END IF
            END DO
          ELSE
! ...       Bessel functions: j_0(0)=1, j_n(0)=0 for n>0
            xg = SQRT( hg(ig) ) * tpiba
            CALL bessel2(xg, ap%rw, fint, lnl, ap%indl, mmax)
            DO l = 1, lnl
              fint(1:mmax,l) = fint(1:mmax,l) * ap%rw(1:mmax)**2 * ap%vrps(1:mmax,l)
              call simpson_fpmd(mmax, fint(:,l), dx, wnl(ig,l))
              ! call simpson(mmax, fint(:,l), ap%rab, wnl(ig,l))
            END DO

! ...       Bessel Test

#if defined __BESSEL_TEST

            ALLOCATE( ftest( mmax, 3 ) )
            WRITE( 11, &
              fmt="(' ir  ','l ',' q         ',' r        ','FPMD                ','FPMD/CP             ','FPMD/PW')")
            DO l = 1, 4
              CALL bessjl  ( xg, ap%rw, ftest(:,1), l, mmax )
              CALL bess    ( xg, l, mmax, ap%rw, ftest(:,2) )
              CALL sph_bes ( mmax, ap%rw, xg, l-1, ftest(:,3) )
              ind = 2
              WRITE( 11, fmt="(I4,I2,2F10.5,3D19.12)" ) &
                 ind, l, xg, ap%rw(ind), ftest(ind,1), ftest(ind,2)/ftest(ind,1), ftest(ind,3)/ftest(ind,1)
              ind = mmax/2
              WRITE( 11, fmt="(I4,I2,2F10.5,3D19.12)" ) &
                 ind, l, xg, ap%rw(ind), ftest(ind,1), ftest(ind,2)/ftest(ind,1), ftest(ind,3)/ftest(ind,1)
              ind = mmax
              WRITE( 11, fmt="(I4,I2,2F10.5,3D19.12)" ) &
                 ind, l, xg, ap%rw(ind), ftest(ind,1), ftest(ind,2)/ftest(ind,1), ftest(ind,3)/ftest(ind,1)
            END DO 
            DEALLOCATE( ftest )

#endif

          END IF
        END DO

        DEALLOCATE(fint)

        RETURN
      END SUBROUTINE nlin_base

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE nlin_stress_base(ap, hg, wnla)

!  (describe briefly what this routine does...)
!  This subroutine is executed only for non-local potentials
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp

! ...   declare subroutine arguments
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl), INTENT(IN)  :: hg(:)
        REAL(dbl), INTENT(OUT) :: wnla(:,:)

! ...   declare other variables
        REAL(dbl), ALLOCATABLE :: fint(:,:)
        REAL(dbl)  xg, dx
        INTEGER ig,mmax,gstart,lnl
        INTEGER l,ll
        INTEGER ir

! ...   end of declarations
!  ----------------------------------------------

        lnl   = ap%lnl
        mmax  = ap%mesh
        dx  = ap%dx

        ALLOCATE( fint(mmax, lnl) )

        DO ig = 1, SIZE( wnla, 1 )
          IF( hg(ig) < gsmall ) THEN
! ...       G=0 (Only if L=1, since otherwise the radial Bessel function JL=0)
            DO l = 1, lnl
              IF(ap%indl(l).EQ.1) THEN
                fint(1:mmax,l) = ap%rw(1:mmax)**2 * ap%vrps(1:mmax,l) 
                call simpson_fpmd(mmax, fint(:,l), dx, wnla(ig, l))
              ELSE
                wnla(ig, l) = 0.d0
              END IF
            END DO
          ELSE
            xg = SQRT(hg(ig)) * tpiba
            CALL bessel3(xg, ap%rw, fint, lnl, ap%indl, mmax)
            DO l = 1, lnl
              fint(1:mmax,l) = fint(1:mmax,l) * ap%rw(1:mmax)**2 * ap%vrps(1:mmax,l)
              call simpson_fpmd(mmax, fint(:,l), dx, wnla(ig,l))
            END DO
          END IF
        END DO

        DEALLOCATE(fint)

        RETURN
      END SUBROUTINE nlin_stress_base

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE nlset_base(ap, wsgset)

!  This subroutine computes the pseudopotential 
!  constants
!  ap%rw**2 * ap%vrps   = [ ( Vpsnl(r) - Vpsloc(r) )* Rps(r) * r^2 ]  
!                       = [ DVpsnl(r) * Rps(r) * r^2 ]  
!  wsgset = 4 PI (2l+1) / < Rps(r) | DVpsnl(r) | Rps(r) >
!
!  for all l and m
!  

!  This subroutine is executed only for non-local
!  pseudopotentials.
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp

! ...   declare subroutine arguments
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl), INTENT(OUT) :: wsgset(:)

! ...   declare other variables
        INTEGER :: ir,i,igh2,l,ll,igh1,igh,mesh,igau,lnl,ltru
        REAL(dbl)  :: alt,blt,one_by_rcl
        REAL(dbl)  :: wsgl, dx
        REAL(dbl), ALLOCATABLE :: gloc(:), fint(:)

! ...   end of declarations
!  ----------------------------------------------

        lnl  = ap%lnl
        igau = ap%igau
        dx   = ap%dx
        mesh = ap%mesh

        wsgset = 0.0d0
        ALLOCATE(fint(mesh))

! ...   Set the normalizing factor "wsg" 

        igh2 = 0

        DO l = 1, lnl

          !  find out the angular momentum (ll-1) of the component stored in position l

          ll = ap%indl( l )  

          fint( 1:mesh ) = ap%rps( 1:mesh, ll ) * ap%vrps( 1:mesh, l ) * ap%rw( 1:mesh )
          call simpson_fpmd( mesh, fint, dx, wsgl )
          !call simpson(mesh, fint, ap%rab, wsgl)

          !  ltru is the true angular momentum quantum number

          ltru = ll - 1  

          igh1 = igh2 + 1
          igh2 = igh1 + ( 2*ltru + 1 ) - 1   

          DO igh = igh1, igh2
            wsgset( igh ) = 4.0d0 * pi * ( 2*ltru + 1 ) / wsgl
          END DO

        END DO

        DEALLOCATE(fint)

        RETURN
      END SUBROUTINE nlset_base     

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE formfn_base(ap, hg, vps, dvps, omega)

!  this routine computes:
!    the local pseudopotentials form factors (vps)
!    their derivatives with respect to cell degrees of freedom (dvps)
!  pseudopotentials are given as numerical tables
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp
        USE io_global, ONLY: stdout


        IMPLICIT NONE

! ...   declare subroutine arguments
        REAL(dbl), INTENT(OUT) :: vps(:),dvps(:)
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl), INTENT(IN) :: hg(:)
        REAL(dbl), INTENT(IN) :: omega

! ... declare external function

        REAL(dbl) :: erf, erfc
        EXTERNAL erf, erfc

        REAL(dbl) :: cclock
        EXTERNAL  :: cclock

! ...   declare other variables
        REAL(dbl), ALLOCATABLE :: jl(:), djl(:), fint(:), func(:)

        REAL(dbl)  :: cost1, cost2, xg, check, fpibo, s1, s2, stot
        INTEGER :: ig, mmax, ir, gstart


! ...   end of declarations
!  ----------------------------------------------

          stot = 0.0d0

          mmax    = ap%mesh
          ALLOCATE( jl(mmax), djl(mmax), fint(mmax), func(mmax) )

          fpibo = fpi / omega
          cost1 = 2.0d0 * ap%zv / (sqrt(pi)*ap%raggio)

          ! WRITE( stdout,*) ' DEBUG formfn_base ',mmax,clfpibo,cost1  ! DEBUG

          DO ir=1,mmax
            cost2    = ap%zv * erf(ap%rw(ir)/ap%raggio)
            func(ir) = ap%rw(ir)**2 * (ap%rw(ir) * ap%vloc(ir) + cost2)
            IF( ABS( ap%rw(ir) ) > 1.d-8 ) THEN
              check = ap%vloc(ir) + cost2 / ap%rw(ir)
            ELSE
              check = ap%vloc(ir) + cost1
            END IF
            IF( ABS ( check ) < 1.d-8) func(ir)=0.d0
            ! WRITE( stdout,*) ir,func(ir),ap%rw(ir),ap%vloc(ir),cost2,ap%rw(ir) * ap%vloc(ir) ! DEBUG
          END DO
          DO ig = 1, SIZE( hg )
            IF( hg(ig) < gsmall ) THEN
              call simpson_fpmd(mmax, func, ap%dx, vps(ig))
              ! call simpson(mmax, func, ap%rab, vps(ig))
              vps(ig) = fpibo * vps(ig)
              fint(1:mmax) = ap%rw(1:mmax)**2 * func(1:mmax)
              call simpson_fpmd(mmax, fint, ap%dx, dvps(ig))
              ! call simpson(mmax, fint, ap%rab, dvps(ig))
              dvps(ig) = fpibo * dvps(ig) / 6.0d0
            ELSE
              xg = SQRT( hg( ig ) ) * tpiba
              CALL bessel1(xg, ap%rw, jl, djl, mmax)
              fint(1:mmax) = func(1:mmax) * jl(1:mmax)
              call simpson_fpmd(mmax, fint, ap%dx, vps(ig))
              ! call simpson(mmax, fint, ap%rab, vps(ig))
              vps(ig) = fpibo * vps(ig)

              fint(1:mmax) = ap%rw(1:mmax) * func(1:mmax) * djl(1:mmax)
              call simpson_fpmd(mmax, fint, ap%dx, dvps(ig))
              ! call simpson(mmax, fint, ap%rab, dvps(ig))
              dvps(ig) = 0.5d0 * fpibo * dvps(ig) / xg
            END IF
            ! WRITE( stdout,*) ig,vps(ig)  ! DEBUG
          END DO

          ! WRITE( stdout,fmt="(/,'* FORMF_BASE TIMING: ',F18.8)" ) stot

          DEALLOCATE( jl, djl, fint, func )

        RETURN
      END SUBROUTINE formfn_base

!  ----------------------------------------------
      SUBROUTINE corecor_base(ap, hg, rhoc1, rhocp, omega)

!  this routine computes:
!  rhoc1(G) = (integral) rho_cc(r) j_0(r,G) r**2 dr
!           = (integral) rho_cc(r) j_0(r,G) r**2 dr/dx dx
!  rhocp(G) = (integral) rho_cc(r) dj_0(r,G)/dG r**2 dr
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE

! ...   declare subroutine arguments
        REAL(dbl), INTENT(OUT) :: rhoc1(:), rhocp(:)
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl), INTENT(IN) :: hg(:)
        REAL(dbl), INTENT(IN) :: omega

! ...   declare other variables
        REAL(dbl), ALLOCATABLE :: jl(:)
        REAL(dbl), ALLOCATABLE :: fint(:)
        REAL(dbl), ALLOCATABLE :: djl(:)
        REAL(dbl), ALLOCATABLE :: funcc(:)

        REAL(dbl)  :: xg, fpibo
        INTEGER :: ig, mmax

! ...   end of declarations
!  ----------------------------------------------

          IF( .NOT. ap%tnlcc ) THEN
            RETURN
          END IF

          mmax    = ap%mesh
          fpibo = fpi / omega
          ALLOCATE( jl(mmax), fint(mmax), djl(mmax), funcc(mmax) )
          funcc(1:mmax) = ap%rw(1:mmax)**3 * ap%rhoc(1:mmax)

          DO ig = 1, SIZE( rhoc1 )
            IF( hg(ig) < gsmall ) THEN
              call simpson_fpmd(mmax, funcc, ap%dx, rhoc1(ig))
              rhoc1(ig) = fpibo * rhoc1(ig)
              rhocp(1)  = 0.0d0
            ELSE
              xg = SQRT( hg( ig ) ) * tpiba
              CALL bessel1(xg, ap%rw, jl, djl, mmax)
              fint (1:mmax) = funcc(1:mmax) * jl(1:mmax)
              call simpson_fpmd(mmax, fint, ap%dx, rhoc1(ig))
              rhoc1(ig) = fpibo * rhoc1(ig)
              fint (1:mmax) = ap%rw(1:mmax) * funcc(1:mmax) * djl(1:mmax)
              call simpson_fpmd(mmax, fint, ap%dx, rhocp(ig))
              rhocp(ig) = fpibo * rhocp(ig)
            END IF
          END DO
          DEALLOCATE( jl, fint, djl, funcc )
        RETURN
      END SUBROUTINE corecor_base
!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE rhops_base(ap, hg, rhops, omega)

!  this routine computes:
!    the ionic Ewald pseudocharges form factors (rhops)
!  pseudopotentials are given as numerical tables
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE
! ...   declare subroutine arguments
        REAL(dbl),            INTENT(OUT) :: rhops(:)
        REAL(dbl),            INTENT(IN)  :: hg(:)
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        REAL(dbl),            INTENT(IN)  :: omega
! ...   declare other variables
        REAL(dbl)   :: r2max, iondens
        INTEGER     :: ig
! ...   end of declarations
!  ----------------------------------------------
          iondens = - ap%zv / omega
          r2max = ( 0.5d0 * ap%raggio * tpiba )**2
          DO ig = 1, SIZE( rhops )
            rhops( ig ) = iondens * EXP (- r2max * hg(ig) ) 
          END DO

        RETURN
      END SUBROUTINE rhops_base
!  ----------------------------------------------
!  ----------------------------------------------

    END MODULE pseudo_base

