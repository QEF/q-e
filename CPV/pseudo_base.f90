!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
   MODULE pseudo_base
!=----------------------------------------------------------------------------=!


      USE kinds
      USE cp_interfaces, ONLY: bessel2, bessel3
      USE constants, ONLY: gsmall, fpi, pi 
      USE cell_base, ONLY: tpiba

      IMPLICIT NONE

      SAVE

      PRIVATE

      PUBLIC :: nlin_base, nlin_stress_base
      PUBLIC :: compute_rhops, formfn, formfa
      PUBLIC :: compute_eself, compute_rhocg



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
        REAL(DP), INTENT(IN) :: hg(:)
        REAL(DP), INTENT(OUT) :: wnl(:,:)

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: fint(:,:)
        REAL(DP)  :: xg, dx
        INTEGER :: ig, mmax, nbeta, l, ind
        REAL(DP), ALLOCATABLE :: ftest(:,:)

! ...   end of declarations
!  ----------------------------------------------

        nbeta  = ap%nbeta
        mmax = ap%mesh
        dx = ap%dx

        ALLOCATE(fint(mmax,nbeta))

        ! WRITE( stdout,*) 'DEBUG nlin_base, tpiba = ', tpiba

        DO ig = 1, SIZE( wnl, 1 )
          IF( hg(ig) < gsmall ) THEN
            DO l = 1,nbeta
! ...         G=0 (Only if l=1, since otherwise the radial Bessel function jl=0)
              IF( ap%lll(l) == 0 ) THEN
                fint(1:mmax,l) = ap%rw(1:mmax) * ap%vrps(1:mmax,l)
                call simpson_cp90( mmax, fint(1,l), ap%rab(1), wnl(ig,l) )
              ELSE
                wnl(ig,l) = 0.d0
              END IF
            END DO
          ELSE
! ...       Bessel functions: j_0(0)=1, j_n(0)=0 for n>0
            xg = SQRT( hg(ig) ) * tpiba
            CALL bessel2(xg, ap%rw, fint, nbeta, ap%lll, mmax)
            DO l = 1, nbeta
              fint(1:mmax,l) = fint(1:mmax,l) * ap%rw(1:mmax) * ap%vrps(1:mmax,l)
              call simpson_cp90( mmax, fint(1,l), ap%rab(1), wnl(ig,l) )
            END DO

! ...       Bessel Test

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
        REAL(DP), INTENT(IN)  :: hg(:)
        REAL(DP), INTENT(OUT) :: wnla(:,:)

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: fint(:,:)
        REAL(DP)  xg, dx
        INTEGER ig,mmax,gstart,nbeta
        INTEGER l,ll
        INTEGER ir

! ...   end of declarations
!  ----------------------------------------------

        nbeta   = ap%nbeta
        mmax  = ap%mesh
        dx  = ap%dx

        ALLOCATE( fint(mmax, nbeta) )

        DO ig = 1, SIZE( wnla, 1 )
          IF( hg(ig) < gsmall ) THEN
! ...       G=0 (Only if L = 0, since otherwise the radial Bessel function JL=0)
            DO l = 1, nbeta
              IF( ap%lll(l) == 0 ) THEN
                fint(1:mmax,l) = ap%rw(1:mmax) * ap%vrps(1:mmax,l) 
                call simpson_cp90( mmax, fint(1,l), ap%rab(1), wnla(ig,l) )
              ELSE
                wnla(ig, l) = 0.d0
              END IF
            END DO
          ELSE
            xg = SQRT(hg(ig)) * tpiba
            CALL bessel3(xg, ap%rw, fint, nbeta, ap%lll, mmax)
            DO l = 1, nbeta
              fint(1:mmax,l) = fint(1:mmax,l) * ap%rw(1:mmax) * ap%vrps(1:mmax,l)
              call simpson_cp90( mmax, fint(1,l), ap%rab(1), wnla(ig,l) )
            END DO
          END IF
        END DO

        DEALLOCATE(fint)

        RETURN
      END SUBROUTINE nlin_stress_base


!-----------------------------------------------------------------------
      subroutine compute_rhocg( rhocb, drhocb, r, rab, rho_atc, gb, omegab, &
                         tpibab2, mesh, ngb, what )

!-----------------------------------------------------------------------

        !  if what == 0 compute rhocb(G)
        !  if what == 1 compute rhocb(G) and drhocb(G)
        !
        !  rhocb(G) = (integral) rho_cc(r) j_0(r,G) r**2 dr
        !           = (integral) rho_cc(r) j_0(r,G) r**2 dr/dx dx
        ! drhocb(G) = (integral) rho_cc(r) dj_0(r,G)/dG r**2 dr

        use kinds,         only: DP
        use constants,     only: fpi
        use control_flags, only: iprsta
        use io_global,     only: stdout

        implicit none
        
        integer,   intent(in)  :: mesh
        integer,   intent(in)  :: ngb
        integer,   intent(in)  :: what
        real(DP), intent(out) :: rhocb( ngb )
        real(DP), intent(out) :: drhocb( ngb )
        real(DP), intent(in)  :: rho_atc( mesh )
        real(DP), intent(in)  :: r( mesh )
        real(DP), intent(in)  :: rab( mesh )
        real(DP), intent(in)  :: gb( ngb )
        real(DP), intent(in)  :: omegab
        real(DP), intent(in)  :: tpibab2
        
        integer :: ig, ir
        real(DP), allocatable :: fint(:), jl(:), djl(:)
        real(DP) :: c, xg
      
        allocate(fint(mesh))
        allocate(jl(mesh))
        if( what == 1 ) then
          allocate(djl(mesh))
        end if

        if( what < 0 .and. what > 1 ) &
          call errore(" compute_rhocg ", " parameter what is out of range ", 1 )

        c = fpi / omegab
        do ig = 1, ngb
           xg = sqrt( gb(ig) * tpibab2 )
           call sph_bes ( mesh, r(1), xg, 0, jl )
           do ir=1,mesh
              fint(ir)=r(ir)**2*rho_atc(ir)*jl(ir)
           end do
           call simpson_cp90( mesh,fint,rab(1),rhocb(ig))
           if( what == 1 ) then
              ! djl = - d j_0(x) /dx = + j_1(x) 
              call sph_bes ( mesh, r(1), xg, +1, djl )
              do ir=1,mesh
                 fint(ir)=r(ir)**3*rho_atc(ir)*djl(ir)
              end do
              call simpson_cp90( mesh, fint, rab(1), drhocb(ig) )
           end if
        end do
        do ig=1,ngb
           rhocb(ig) = c * rhocb(ig)
        end do
        if( what == 1 ) then
           do ig=1,ngb
              drhocb(ig) = c * drhocb(ig)
           end do
        end if

        if(iprsta >= 4) &
             WRITE( stdout,'(a,f12.8)') ' integrated core charge= ',omegab*rhocb(1)
        
        deallocate( jl, fint )
        if( what == 1 ) then
           deallocate(djl)
        end if

        return
      end subroutine compute_rhocg

!-----------------------------------------------------------------------
      subroutine compute_rhops( rhops, drhops, zv, rcmax, g, omega, tpiba2, ngs, tpre )
!-----------------------------------------------------------------------
!
        use kinds, only: DP
        !
        implicit none
        integer,   intent(in)  :: ngs
        logical,   intent(in)  :: tpre
        real(DP), intent(in)  ::  g( ngs )
        real(DP), intent(out) ::  rhops( ngs )
        real(DP), intent(out) :: drhops( ngs )
        real(DP), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(DP) :: r2new
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
      end subroutine compute_rhops



!-----------------------------------------------------------------------
      FUNCTION compute_eself( na, zv, rcmax, nsp )
!-----------------------------------------------------------------------
        !
        !     calculation of gaussian selfinteraction
        !
        USE constants, ONLY: pi
        !
        IMPLICIT NONE
        REAL (DP) :: compute_eself
        !
        INTEGER,    INTENT(IN) :: nsp
        INTEGER,    INTENT(IN) :: na( nsp )
        REAL (DP), INTENT(IN) :: zv( nsp )
        REAL (DP), INTENT(IN) :: rcmax( nsp )
        !
        REAL (DP) :: eself
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
                         tpiba2, mesh, ngs, oldvan, tpre )
!-----------------------------------------------------------------------
!
        !computes the form factors of pseudopotential (vps),
        !         also calculated the derivative of vps with respect to
        !         g^2 (dvps)
        !
        use kinds, only: DP
        use constants, only: pi, fpi, gsmall
        !
        implicit none
        integer,   intent(in)  :: ngs
        integer,   intent(in)  :: mesh
        logical,   intent(in)  :: oldvan
        logical,   intent(in)  :: tpre
        real(DP), intent(in)  ::  g( ngs )
        real(DP), intent(in)  ::  r( mesh )
        real(DP), intent(in)  ::  rab( mesh )
        real(DP), intent(in)  ::  vloc_at( mesh )
        real(DP), intent(out) ::  vps( ngs )
        real(DP), intent(out) :: dvps( ngs )
        real(DP), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(DP) :: xg
        integer   :: ig, ir, irmax
        real(8), allocatable:: f(:),vscr(:), figl(:)
        real(8), allocatable:: df(:), dfigl(:)
        real(8), external :: erf
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
              call herman_skillman_int( mesh, f, rab, figl(ig) )
              if(tpre) call herman_skillman_int( mesh, df, rab, dfigl(ig) )
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
              call herman_skillman_int( mesh, f, rab(1), figl(ig) )
              if(tpre) call herman_skillman_int( mesh, df, rab(1), dfigl(ig) )
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
      end subroutine formfn




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

        use kinds, only: DP
        use constants, only: pi, fpi, gsmall
        !
        implicit none
        integer,   intent(in)  :: ngs, gstart
        logical,   intent(in)  :: tpre
        real(DP), intent(in)  ::  g( ngs )
        real(DP), intent(in)  ::  rc1, rc2
        real(DP), intent(in)  ::  wrc1, wrc2
        real(DP), intent(in)  ::  rcl( 3 ), al( 3 ), bl( 3 )
        real(DP), intent(out) ::  vps( ngs )
        real(DP), intent(out) :: dvps( ngs )
        real(DP), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(DP) :: r2max, r21, r22, gps, sfp, r2l, ql, el, par, sp
        real(DP) :: emax, e1, e2, fpibg, dgps, dsfp
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
      end subroutine formfa



!=----------------------------------------------------------------------------=!
   END MODULE pseudo_base
!=----------------------------------------------------------------------------=!

