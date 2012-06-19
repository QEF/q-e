!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE pseudo_base
!=----------------------------------------------------------------------------=!


      USE kinds
      USE constants, ONLY: gsmall, fpi, pi 
      USE cell_base, ONLY: tpiba

      IMPLICIT NONE

      SAVE

      PRIVATE

      PUBLIC :: compute_rhops, formfn, formfa
      PUBLIC :: compute_eself, compute_rhocg



!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


    subroutine compute_rhocg( rhocb, drhocb, r, rab, rho_atc, gb, omegab, &
                         tpibab2, mesh, ngb, what )


        !  if what == 0 compute rhocb(G)
        !  if what == 1 compute rhocb(G) and drhocb(G)
        !
        !  rhocb(G) = (integral) rho_cc(r) j_0(r,G) r**2 dr
        !           = (integral) rho_cc(r) j_0(r,G) r**2 dr/dx dx
        ! drhocb(G) = (integral) rho_cc(r) dj_0(r,G)/dG r**2 dr

        use kinds,         only: DP
        use constants,     only: fpi
        use control_flags, only: iverbosity
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
      
!$omp parallel default(none) private(ig,c,xg,fint,jl,djl,ir) &
!$omp          shared(mesh,what,omegab,ngb,tpibab2,gb,r,rho_atc,rhocb,rab,drhocb)

        allocate(fint(mesh))
        allocate(jl(mesh))
        if( what == 1 ) then
          allocate(djl(mesh))
        end if

        if( what < 0 .and. what > 1 ) &
          call errore(" compute_rhocg ", " parameter what is out of range ", 1 )

        c = fpi / omegab
!$omp do
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
!$omp do
        do ig=1,ngb
           rhocb(ig) = c * rhocb(ig)
        end do
        if( what == 1 ) then
!$omp do
           do ig=1,ngb
              drhocb(ig) = c * drhocb(ig)
           end do
        end if
        deallocate( jl, fint )
        if( what == 1 ) then
           deallocate(djl)
        end if
!$omp end parallel

        if(iverbosity > 2) WRITE( stdout,'(a,f12.8)') &
                           ' integrated core charge= ',omegab*rhocb(1)

        return
      end subroutine compute_rhocg



!-----------------------------------------------------------------------



      subroutine compute_rhops( rhops, drhops, zv, rcmax, g, omega, tpiba2, ngs, tpre )
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
        r2new = 0.25d0 * tpiba2 * rcmax**2
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
      subroutine formfn( r, rab, vloc_at, zv, rcmax, g, omega, &
                         tpiba2, mesh, ngs, oldvan, tpre, vps, dv0, dvps )
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
        real(DP), intent(out) :: dv0
        real(DP), intent(in)  :: zv, rcmax, omega, tpiba2
        !
        real(DP) :: xg
        integer   :: ig, ir, irmax
        real(DP), allocatable:: f(:),vscr(:), figl(:)
        real(DP), allocatable:: df(:), dfigl(:)
        real(DP), external :: qe_erf, qe_erfc
!
        allocate( vscr(mesh), figl(ngs) )
        if (tpre) then
           allocate( dfigl(ngs) )
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
          vscr(ir) = 0.5d0 * r(ir) * vloc_at(ir) + zv * qe_erf( r(ir) / rcmax )
        end do
        do ir = irmax + 1, mesh
          vscr(ir)=0.0d0
        end do
        !
        ! ... In CP the G=0 value of the Hartree+local pseudopotential
        ! ... is not set to its correct value, the "alpha Z" term, but
        ! ... to a different value. This has no effect on the energy
        ! ... of a neutral system as long as all terms are consistent
        ! ... but it yields a different alignment of levels and, only 
        ! ... in charged system, a different energy. 
        ! ... dv0 is the correction to the G=0 term in CP needed to
        ! ...  reproduce the results from other PW codes
        !
!$omp parallel default(none) private( ig, xg, ir, f, df ) &
!$omp          shared( irmax, r, rcmax, mesh, oldvan, rab, dv0, tpiba2, g, ngs, vscr, tpre, zv, figl, vps, dvps, omega, dfigl )

        allocate( f(mesh) )
        if (tpre) then
           allocate( df(mesh) )
        end if
        DO ir = 1, irmax
           f(ir) = fpi * ( zv * qe_erfc( r(ir)/rcmax ) ) * r(ir)
        END DO
        DO ir = irmax + 1, mesh
          f(ir)=0.0d0
        END DO

!$omp master
        IF ( oldvan ) THEN
           CALL herman_skillman_int( mesh, f, rab, dv0 )
        ELSE
           CALL simpson_cp90( mesh, f, rab, dv0 )
        END IF
!$omp end master
        !

!$omp do
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
              f(ir)  = 0.0d0
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
!$omp do
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
              dvps(ig) = fpi * dfigl(ig) / ( omega * xg ) - 0.5d0 * vps(ig) / (xg*xg)
            endif
          end if
        end do

        deallocate( f )
        if (tpre) then
           deallocate( df )
        end if

!$omp end parallel
        !
        deallocate( figl, vscr )
        if (tpre) then
           deallocate( dfigl )
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
        ! BHS pseudopotentials (fourier transformed analytically)

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
            sp  = ( pi * r2l )**1.5d0 * el / omega
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

