!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

#undef __BESSEL_TEST

!=----------------------------------------------------------------------------=!
   MODULE pseudo_base
!=----------------------------------------------------------------------------=!


      USE kinds
      USE bessel_functions, ONLY: bessel1, bessel2, bessel3
      USE constants, ONLY: gsmall, fpi, pi 
      USE cell_base, ONLY: tpiba

      IMPLICIT NONE

      SAVE

      PRIVATE

      PUBLIC :: nlin_base, nlin_stress_base, nlset_base
      PUBLIC :: corecor_base, compute_rhops, formfn, formfa
      PUBLIC :: compute_eself



!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!



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



!-----------------------------------------------------------------------
      subroutine compute_rhocg( rhocb, r, rab, rho_atc, gb, omegab, &
                         tpibab2, mesh, ngb )

!-----------------------------------------------------------------------

        !  rhoc1(G) = (integral) rho_cc(r) j_0(r,G) r**2 dr
        !           = (integral) rho_cc(r) j_0(r,G) r**2 dr/dx dx

        use kinds,         only: dbl
        use constants,     only: fpi
        use control_flags, only: iprsta
        use io_global,     only: stdout

        implicit none
        
        integer,   intent(in)  :: mesh
        integer,   intent(in)  :: ngb
        real(dbl), intent(out) :: rhocb( ngb )
        real(dbl), intent(in)  :: rho_atc( mesh )
        real(dbl), intent(in)  :: r( mesh )
        real(dbl), intent(in)  :: rab( mesh )
        real(dbl), intent(in)  :: gb( ngb )
        real(dbl), intent(in)  :: omegab
        real(dbl), intent(in)  :: tpibab2
        
        integer :: ig, ir
        real(dbl), allocatable :: fint(:), jl(:)
        real(dbl) :: c, xg
      
        allocate(fint(mesh))
        allocate(jl(mesh))

        c = fpi / omegab
        do ig = 1, ngb
          xg = sqrt( gb(ig) * tpibab2 )
          call sph_bes ( mesh, r(1), xg, 0, jl )
          do ir=1,mesh
            fint(ir)=r(ir)**2*rho_atc(ir)*jl(ir)
          end do
          call simpson_cp90( mesh,fint,rab(1),rhocb(ig))
        end do
        do ig=1,ngb
          rhocb(ig)=c*rhocb(ig)
        end do
        if(iprsta >= 4) &
          WRITE( stdout,'(a,f12.8)') ' integrated core charge= ',omegab*rhocb(1)

        deallocate( jl, fint )

        return
      end subroutine



!-----------------------------------------------------------------------
      subroutine compute_rhops( rhops, drhops, zv, rcmax, g, omega, tpiba2, ngs, tpre )
!-----------------------------------------------------------------------
!
        use kinds, only: dbl
        !
        implicit none
        integer,   intent(in)  :: ngs
        logical,   intent(in)  :: tpre
        real(dbl), intent(in)  ::  g( ngs )
        real(dbl), intent(out) ::  rhops( ngs )
        real(dbl), intent(out) :: drhops( ngs )
        real(dbl), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(dbl) :: r2new
        integer   :: ig
        !
        r2new = 0.25 * tpiba2 * rcmax**2
        do ig = 1, ngs
          rhops(ig) = - zv * exp( -r2new * g(ig) ) / omega
        end do
        if(tpre) then
          drhops( 1:ngs ) = - rhops( 1:ngs ) * r2new / tpiba2
        endif
        !
        return
      end subroutine



!-----------------------------------------------------------------------
      FUNCTION compute_eself( na, zv, rcmax, nsp )
!-----------------------------------------------------------------------
        !
        !     calculation of gaussian selfinteraction
        !
        USE constants, ONLY: pi
        !
        IMPLICIT NONE
        REAL (dbl) :: compute_eself
        !
        INTEGER,    INTENT(IN) :: nsp
        INTEGER,    INTENT(IN) :: na( nsp )
        REAL (dbl), INTENT(IN) :: zv( nsp )
        REAL (dbl), INTENT(IN) :: rcmax( nsp )
        !
        REAL (dbl) :: eself
        INTEGER :: is
        !
        eself = 0.0d0
        DO is = 1, nsp
          eself = eself + DBLE( na( is ) ) * zv( is )**2 / rcmax( is )
        END DO
        eself = eself / SQRT( 2.0d0 * pi )
        !
        compute_eself = eself
        RETURN
      END FUNCTION compute_eself




!-----------------------------------------------------------------------
      subroutine formfn( vps, dvps, r, rab, vloc_at, zv, rcmax, g, omega, &
                         tpiba2, cmesh, mesh, ngs, oldvan, tpre )
!-----------------------------------------------------------------------
!
        !computes the form factors of pseudopotential (vps),
        !         also calculated the derivative of vps with respect to
        !         g^2 (dvps)
        !
        use kinds, only: dbl
        use constants, only: pi, fpi, gsmall
        !
        implicit none
        integer,   intent(in)  :: ngs
        integer,   intent(in)  :: mesh
        logical,   intent(in)  :: oldvan
        logical,   intent(in)  :: tpre
        real(dbl), intent(in)  ::  g( ngs )
        real(dbl), intent(in)  ::  r( mesh )
        real(dbl), intent(in)  ::  rab( mesh )
        real(dbl), intent(in)  ::  vloc_at( mesh )
        real(dbl), intent(out) ::  vps( ngs )
        real(dbl), intent(out) :: dvps( ngs )
        real(dbl), intent(in)  :: zv, rcmax, omega, tpiba2, cmesh
        !
        real(dbl) :: xg
        integer   :: ig, ir, irmax
        real(kind=8), allocatable:: f(:),vscr(:), figl(:)
        real(kind=8), allocatable:: df(:), dfigl(:)
        real(kind=8), external :: erf
!
        allocate( figl(ngs), f(mesh), vscr(mesh) )
        if (tpre) then
           allocate( dfigl(ngs), df(mesh) )
        end if
        !
        !     definition of irmax: gridpoint beyond which potential is zero
        !
        irmax = 0
        do ir = 1, mesh
          if( r( ir ) < 10.0d0 ) irmax = ir
        end do
        !
        do ir = 1, irmax
          vscr(ir) = 0.5d0 * r(ir) * vloc_at(ir) + zv * erf( r(ir) / rcmax )
        end do
        do ir = irmax + 1, mesh
          vscr(ir)=0.0
        end do

        do ig = 1, ngs
          xg = sqrt( g(ig) * tpiba2 )
          if( xg < gsmall ) then
            !
            !     g=0
            !
            do ir = 1, irmax
              f(ir) = vscr(ir) * r(ir)
              if( tpre ) then
                df(ir) = vscr(ir) * r(ir) ** 3
              endif
            end do
            do ir = irmax + 1, mesh
              f(ir)  = 0.0
              if( tpre ) then
                df(ir) = 0.0d0
              end if
            end do
            !
            if ( oldvan ) then
              call herman_skillman_int( mesh, cmesh, f,  figl(ig) )
              if(tpre) call herman_skillman_int( mesh, cmesh, df, dfigl(ig) )
            else
              call simpson_cp90( mesh, f, rab,  figl(ig) )
              if(tpre) call simpson_cp90( mesh, df, rab, dfigl(ig) )
            end if
            !
          else
            !
            !     g>0
            !
            do ir = 1, mesh
              f(ir) = vscr(ir) * sin( r(ir) * xg )
              if( tpre ) then
                df(ir) = vscr(ir) * cos( r(ir) * xg ) * 0.5d0 * r(ir) / xg
              endif
            end do
            !
            if ( oldvan ) then
              call herman_skillman_int( mesh, cmesh, f, figl(ig) )
              if(tpre) call herman_skillman_int( mesh, cmesh, df, dfigl(ig) )
            else
              call simpson_cp90(mesh,f,rab(1),figl(ig))
              if(tpre) call simpson_cp90(mesh,df,rab(1),dfigl(ig))
            end if
            !
          end if
        end do
        !
        do ig = 1, ngs
          xg = sqrt( g(ig) * tpiba2 )
          if( xg < gsmall ) then
            !
            !     g=0
            !
            vps(ig)   = fpi *  figl(ig) / omega
            if(tpre)then
              dvps(ig) = - fpi * dfigl(ig) / omega / 6.0d0  !  limit ( xg -> 0 ) dvps( xgi )
            end if
            !
          else
            !
            !     g>0
            !
            vps(ig)  = fpi *  figl(ig) / ( omega * xg )
            if(tpre)then
              dvps(ig) = fpi * dfigl(ig) / ( omega * xg ) - 0.5 * vps(ig) / (xg*xg)
            endif
          end if
        end do
        !
        deallocate( figl, f, vscr )
        if (tpre) then
           deallocate( dfigl, df )
        end if
        !
      return
      end subroutine




!-----------------------------------------------------------------------
      subroutine formfa( vps, dvps, rc1, rc2, wrc1, wrc2, rcl, al, bl, &
                         zv, rcmax, g, omega, tpiba2, ngs, gstart, tpre )
!-----------------------------------------------------------------------
!
        !computes the form factors of pseudopotential (vps),
        !         also calculated the derivative of vps with respect to
        !         g^2 (dvps)
        !
        ! bhs pseudopotentials (fourier transformed analytically)

        use kinds, only: dbl
        use constants, only: pi, fpi, gsmall
        !
        implicit none
        integer,   intent(in)  :: ngs, gstart
        logical,   intent(in)  :: tpre
        real(dbl), intent(in)  ::  g( ngs )
        real(dbl), intent(in)  ::  rc1, rc2
        real(dbl), intent(in)  ::  wrc1, wrc2
        real(dbl), intent(in)  ::  rcl( 3 ), al( 3 ), bl( 3 )
        real(dbl), intent(out) ::  vps( ngs )
        real(dbl), intent(out) :: dvps( ngs )
        real(dbl), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(dbl) :: r2max, r21, r22, gps, sfp, r2l, ql, el, par, sp
        real(dbl) :: emax, e1, e2, fpibg, dgps, dsfp
        integer   :: ib, ig

        r2max = rcmax**2
        r21   = rc1**2
        r22   = rc2**2

        !
        !     g = 0
        !
        if (gstart == 2) then
          gps = - zv * pi * ( - wrc2 * r22 - wrc1 * r21 + r2max ) / omega
          sfp = 0.0d0
          do ib = 1, 3
            r2l = rcl( ib )**2
            ql  = 0.25d0 * r2l * g(1) * tpiba2
            el  = exp( -ql )
            par = al( ib ) + bl( ib ) * r2l * ( 1.5d0 - ql )
            sp  = ( pi * r2l )**1.5 * el / omega
            sfp = sp * par + sfp
          end do
          vps(1) = gps + sfp
        end if
        !
        !     g > 0
        !
        do ig = gstart, ngs
          !
          emax  = exp ( -0.25d0 * r2max * g(ig) * tpiba2 )
          e1    = exp ( -0.25d0 * r21   * g(ig) * tpiba2 )
          e2    = exp ( -0.25d0 * r22   * g(ig) * tpiba2 )
          fpibg = fpi / ( g(ig) * tpiba2 )
          gps   = - zv * ( wrc1 * e1 - emax + wrc2 * e2 ) / omega
          gps   = gps * fpibg
          !
          if(tpre) then
            dgps = - gps / ( tpiba2 * g(ig) ) +  fpibg * zv *                  &
                 &   ( wrc1 * r21 * e1 - r2max * emax + wrc2 * r22 * e2 ) *    &
                 &   0.25d0 / omega
          end if
          !
          sfp   = 0.0d0
          dsfp  = 0.0d0
          !
          do ib = 1, 3
            r2l   = rcl( ib )**2
            ql    = 0.25d0 * r2l * g(ig) * tpiba2
            par   = al( ib ) + bl( ib ) * r2l * ( 1.5d0 - ql )
            sp    = ( pi * r2l )**1.5d0 * exp( -ql ) / omega
            sfp   = sp * par + sfp
            if(tpre) then
              dsfp  = dsfp - sp * ( par + bl( ib ) * r2l ) * ql / ( tpiba2 * g(ig) )
            end if
          end do
          !
          vps(ig) = sfp + gps
          if(tpre) dvps(ig) = dsfp + dgps
          !
        end do
!
      return
      end subroutine



!=----------------------------------------------------------------------------=!
   END MODULE pseudo_base
!=----------------------------------------------------------------------------=!

