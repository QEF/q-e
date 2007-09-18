!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE paw_v_xc( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                 nrxx, nl, ngm, g, nspin, alat, omega, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  ! ... Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,         ONLY : domag
  USE kinds,            ONLY : DP
  USE funct

#ifdef EXX
  USE exx,              ONLY: exxalfa
#endif

  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                         nrxx, ngm, nl(ngm)
    !  nspin=1 :unpolarized, =2 :spin-polarized
    ! the FFT indices
    ! the true FFT array dimensions
    ! the total dimension
    ! the number of G vectors
    ! correspondence G <-> FFT
  REAL (DP), INTENT(IN) :: rho(nrxx,nspin), rho_core(nrxx), &
                                g(3,ngm), alat, omega
    ! the valence charge
    ! the core charge
    ! the G vectors
    ! the length of the cell
    ! the volume of the cell
  !
  REAL (DP), INTENT(OUT) :: v(nrxx,nspin), vtxc, etxc
    ! V_xc potential
    ! integral V_xc * rho
    ! E_xc energy
  !
  ! ... local variables
  !
  REAL (DP) :: rhox, arhox, zeta, amag, vs, ex, ec, vx(2), vc(2), &
                    rhoneg(2)
    ! the total charge in each point
    ! the absolute value of the charge
    ! the absolute value of the charge
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, is, ig, ipol
    ! counter on mesh points
    ! counter on spin polarizations
    ! counter on G vectors
    ! counter on nspin
    ! number of points with wrong zeta/charge
  !
  REAL (DP), PARAMETER :: vanishing_charge = 1.D-10, &
                               vanishing_mag    = 1.D-20
  !
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     !
     ! ... spin-unpolarized case
     !
     DO ir = 1, nrxx
        !
        rhox = rho(ir,1) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           CALL xc( arhox, ex, ec, vx(1), vc(1) )
           !
#if defined EXX
           v(ir,1) = e2 * ( (1.d0-exxalfa)*vx(1) + vc(1) )
           !
           etxc = etxc + e2 * ( (1.d0-exxalfa)*ex + ec ) * rhox
#else
           v(ir,1) = e2 * ( vx(1) + vc(1) )
           !
           etxc = etxc + e2 * ( ex + ec ) * rhox
#endif
           !
           vtxc = vtxc + v(ir,1) * rho(ir,1)
           !
        ENDIF
        !
        IF ( rho(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho(ir,1)
        !
     END DO
     !
  ELSE IF (nspin == 2) THEN
     !
     ! ... spin-polarized case
     !
     DO ir = 1, nrxx
        !
        rhox = rho(ir,1) + rho(ir,2) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = ( rho(ir,1) - rho(ir,2) ) / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
           !
           IF ( rho(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho(ir,1)
           IF ( rho(ir,2) < 0.D0 ) rhoneg(2) = rhoneg(2) - rho(ir,2)
           !
               !xc_spin (rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw)
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
#ifdef EXX
           DO is = 1, nspin
              !
              v(ir,is) = e2 * ( (1.d0-exxalfa)*vx(is) + vc(is) )
              !
           END DO
           !
           etxc = etxc + e2 * ( (1.d0-exxalfa)*ex + ec ) * rhox
#else
           DO is = 1, nspin
              !
              v(ir,is) = e2 * ( vx(is) + vc(is) )
              !
           END DO
           !
           etxc = etxc + e2 * ( ex + ec ) * rhox
#endif
           !
           vtxc = vtxc + v(ir,1) * rho(ir,1) + v(ir,2) * rho(ir,2)
           !
        END IF
        !
     END DO
     !
  ELSE IF (nspin==4) THEN
     !
     ! ... noncolinear case
     !
     DO ir = 1,nrxx
        !
        amag = SQRT( rho(ir,2)**2 + rho(ir,3)**2 + rho(ir,4)**2 )
        !
        rhox = rho(ir,1) + rho_core(ir)
        !
        IF ( rho(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho(ir,1)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = amag / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) THEN
              !
              rhoneg(2) = rhoneg(2) + 1.D0 / omega
              !
              zeta = SIGN( 1.D0, zeta )
              !
           END IF
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           vs = 0.5D0 * ( vx(1) + vc(1) - vx(2) - vc(2) )
           !
           v(ir,1) = e2 * ( 0.5D0 * ( vx(1) + vc(1) + vx(2) + vc(2 ) ) )
           !
           IF ( amag > vanishing_mag ) THEN
             !
              DO ipol = 2, 4
                 !
                 v(ir,ipol) = e2 * vs * rho(ir,ipol) / amag
                 !
              END DO
              !
           END IF
           !
           etxc= etxc + e2 * ( ex + ec ) * rhox
           vtxc= vtxc + v(ir,1) * rho(ir,1)
           !
        END IF
        !
     END DO
     !
  END IF
  !
  CALL reduce( 2, rhoneg )
  !
  rhoneg(:) = rhoneg(:) * omega / ( nr1 * nr2 * nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,4X," negative rho (up, down): ",2E10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( nr1 * nr2 * nr3 )
  etxc = omega * etxc / ( nr1 * nr2 * nr3 )
  !
  ! ... add gradient corrections (if any)
  !
  !
  CALL paw_gradcorr( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                 nrxx, nl, ngm, g, alat, omega, nspin, etxc, vtxc, v )
  !
  CALL reduce( 1, vtxc )
  CALL reduce( 1, etxc )
  !
  RETURN
  !
END SUBROUTINE paw_v_xc


!----------------------------------------------------------------------------
SUBROUTINE paw_gradcorr( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                     nrxx, nl, ngm, g, alat, omega, nspin, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  USE constants, ONLY : e2
  USE kinds,     ONLY : DP
  USE funct !,     ONLY : dft_is_gradient, get_igcc  !igcx, igcc
  USE spin_orb, ONLY : domag

#ifdef EXX
  USE exx,       ONLY: lexx, exxalfa
#endif

  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN)    :: nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                                   nrxx, ngm, nl(ngm), nspin
  REAL (DP), INTENT(IN)    :: rho_core(nrxx), g(3,ngm), alat, omega
  REAL (DP), INTENT(OUT)   :: v(nrxx,nspin), vtxc, etxc
  REAL (DP), INTENT(INOUT) :: rho(nrxx,nspin)
  !
  INTEGER :: k, ipol, is, nspin0
  !
  REAL (DP), ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL (DP), ALLOCATABLE :: rhoout(:,:), segni(:), vgg(:,:), vsave(:,:)
  !
  REAL (DP) :: grho2(2), sx, sc, v1x, v2x, v1c, v2c, &
                    v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw , &
                    etxcgc, vtxcgc, segno, arho, fac, zeta, rh, grh2, amag 
  !
  REAL (DP) :: v2cup, v2cdw,  v2cud, rup, rdw, &
                    grhoup, grhodw, grhoud, grup, grdw

  REAL (DP), PARAMETER :: epsr = 1.D-6, &
                               epsg = 1.D-10
  !
  !
  IF ( .not. dft_is_gradient() ) RETURN
  !
  etxcgc = 0.D0
  vtxcgc = 0.D0
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2

  ALLOCATE(    h( 3, nrxx, nspin0) )
  ALLOCATE( grho( 3, nrxx, nspin0) )
  ALLOCATE( rhoout( nrxx, nspin0) )
  IF (nspin==4.AND.domag) THEN
     ALLOCATE( segni( nrxx ) )
     ALLOCATE( vgg( nrxx, nspin0 ) )
     ALLOCATE( vsave( nrxx, nspin ) )
     vsave=v
     v=0.d0
  ENDIF 

  IF (nspin==4.AND.domag) THEN
     CALL compute_rho(rho,rhoout,segni,nrxx) 
  ELSE
     DO is=1,nspin0
        rhoout(:,is) = rho(:,is)
     END DO
  ENDIF
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  fac = 1.D0 / DBLE( nspin0 )
  !
  DO is = 1, nspin0
     !
     rhoout(:,is) = fac * rho_core(:) + rhoout(:,is)
     !
     CALL paw_gradient( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                    rhoout(1,is), ngm, g, nl, alat, grho(1,1,is) )
     !
  END DO
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     DO k = 1, nrxx
        !
        arho = ABS( rhoout(k,1) )
        !
        IF ( arho > epsr ) THEN
           !
           grho2(1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
           !
           IF ( grho2(1) > epsg ) THEN
              !
              segno = SIGN( 1.D0, rhoout(k,1) )
              !
                 ! gcxc (rho, grho, sx, sc, v1x, v2x, v1c, v2c)
              CALL gcxc( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
#if defined (EXX)
              if (lexx) then
                 sx  = (1.d0-exxalfa)*sx
                 v1x = (1.d0-exxalfa)*v1x
                 v2x = (1.d0-exxalfa)*v2x
!                 sc  = (1.d0-exxalfa)*sc
!                 v1c = (1.d0-exxalfa)*v1c
!                 v2c = (1.d0-exxalfa)*v2c
              end if
#endif
              !
              ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
              !
              v(k,1) = v(k,1) + e2 * ( v1x + v1c )
              !
              ! ... h contains :
              !
              ! ...    D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
              !
              h(:,k,1) = e2 * ( v2x + v2c ) * grho(:,k,1)
              !
              vtxcgc = vtxcgc+e2*( v1x + v1c ) * ( rhoout(k,1) - rho_core(k) )
              etxcgc = etxcgc+e2*( sx + sc ) * segno
              !
           END IF
           !
        ELSE
           !
           h(:,k,1) = 0.D0
           !
        END IF
        !
     END DO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
     DO k = 1, nrxx
        !
        rh = rhoout(k,1) + rhoout(k,2)
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
                       grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
#if defined (EXX)
        if (lexx) then
           sx    = (1.d0-exxalfa)*sx
           v1xup = (1.d0-exxalfa)*v1xup
           v1xdw = (1.d0-exxalfa)*v1xdw
           v2xup = (1.d0-exxalfa)*v2xup
           v2xdw = (1.d0-exxalfa)*v2xdw
        end if
#endif
        !
        IF ( rh > epsr ) THEN
           !
           IF ( get_igcc() == 3 ) THEN
              !
              rup = rhoout(k,1)
              rdw = rhoout(k,2)
              !
              grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
              grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
              !
              grhoud = grho(1,k,1) * grho(1,k,2) + &
                       grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
              !
              CALL gcc_spin_more( rup, rdw, grhoup, grhodw, grhoud, &
                                  sc, v1cup, v1cdw, v2cup, v2cdw, v2cud )
              !
           ELSE
              !
              zeta = ( rhoout(k,1) - rhoout(k,2) ) / rh
              if (nspin.eq.4.and.domag) zeta=abs(zeta)*segni(k)
              !
              grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                     ( grho(2,k,1) + grho(2,k,2) )**2 + &
                     ( grho(3,k,1) + grho(3,k,2) )**2
              !
              CALL gcc_spin( rh, zeta, grh2, sc, v1cup, v1cdw, v2c )
              !
              v2cup = v2c
              v2cdw = v2c
              v2cud = v2c
              !
           END IF
           !
        ELSE
           !
           sc    = 0.D0
           v1cup = 0.D0
           v1cdw = 0.D0
           v2c   = 0.D0
           v2cup = 0.D0
           v2cdw = 0.D0
           v2cud = 0.D0
           !
        ENDIF
#if defined (EXX)
!        if (lexx) then
!           sc    = (1.d0-exxalfa)*sc
!           v1cup = (1.d0-exxalfa)*v1cup
!           v1cdw = (1.d0-exxalfa)*v1cdw
!           v2c   = (1.d0-exxalfa)*v2c
!           v2cup = (1.d0-exxalfa)*v2cup
!           v2cdw = (1.d0-exxalfa)*v2cdw
!           v2cud = (1.d0-exxalfa)*v2cud
!        end if
#endif
        !
        ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        v(k,1) = v(k,1) + e2 * ( v1xup + v1cup )
        v(k,2) = v(k,2) + e2 * ( v1xdw + v1cdw )
        !
        ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3
           !
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2 * ( ( v2xup + v2cup ) * grup + v2cud * grdw )
           h(ipol,k,2) = e2 * ( ( v2xdw + v2cdw ) * grdw + v2cud * grup )
           !
        END DO
        !
        vtxcgc = vtxcgc + &
                 e2 * ( v1xup + v1cup ) * ( rhoout(k,1) - rho_core(k) * fac )
        vtxcgc = vtxcgc + &
                 e2 * ( v1xdw + v1cdw ) * ( rhoout(k,2) - rho_core(k) * fac )
        etxcgc = etxcgc + e2 * ( sx + sc )
        !
     END DO
     !
  END IF
  !
  DO is = 1, nspin0
     !
     rhoout(:,is) = rhoout(:,is) - fac * rho_core(:)
     !
  END DO
  !
  DEALLOCATE( grho )
  !
  ALLOCATE( dh( nrxx ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
     CALL paw_grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                    nrxx, h(1,1,is), ngm, g, nl, alat, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     vtxcgc = vtxcgc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  !
  vtxc = vtxc + omega * vtxcgc / ( nr1 * nr2 * nr3 )
  etxc = etxc + omega * etxcgc / ( nr1 * nr2 * nr3 )

  IF (nspin==4.AND.domag) THEN
     DO is=1,nspin0
        vgg(:,is)=v(:,is)
     ENDDO
     v=vsave
     DO k=1,nrxx
        v(k,1)=v(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
        amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag.GT.1.d-12) THEN
           v(k,2)=v(k,2)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
           v(k,3)=v(k,3)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
           v(k,4)=v(k,4)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
        ENDIF
     ENDDO
  ENDIF
  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( rhoout )
  IF (nspin==4.and.domag) THEN
     DEALLOCATE( vgg )
     DEALLOCATE( vsave )
     DEALLOCATE( segni )
  ENDIF

  !
  RETURN
  !
END SUBROUTINE paw_gradcorr
!
!----------------------------------------------------------------------------
SUBROUTINE paw_gradient( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                     nrxx, a, ngm, g, nl, alat, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE wvfct,     ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN)     :: nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                                    nrxx, ngm, nl(ngm)
  REAL (DP), INTENT(IN)     :: a(nrxx), g(3,ngm), alat
  REAL (DP), INTENT(OUT)    :: ga(3,nrxx)
  !
  INTEGER                        :: n, ipol
  COMPLEX (DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 )
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.D0
  !
  DO ipol = 1, 3
     !
     gaux(:) = 0.D0
     !
     gaux(nl(:)) = g(ipol,:) * &
          CMPLX( -AIMAG( aux(nl(:)) ) , DBLE( aux(nl(:)) ) )
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( DBLE( gaux(nl(:)) ) , -AIMAG( gaux(nl(:)) ) )
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL cft3( gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = ga(ipol,:) + tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE paw_gradient
!
!----------------------------------------------------------------------------
SUBROUTINE paw_grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                     nrxx, a, ngm, g, nl, alat, da )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates da = \sum_i \grad_i a_i in R-space
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE wvfct,     ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN)     :: nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                                    nrxx, ngm, nl(ngm)
  REAL (DP), INTENT(IN)     :: a(3,nrxx), g(3,ngm), alat
  REAL (DP), INTENT(OUT)    :: da(nrxx)
  !
  INTEGER                        :: n, ipol
  COMPLEX (DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  gaux(:) = 0.D0
  !
  DO ipol = 1, 3
     !
     aux = CMPLX( a(ipol,:), 0.D0 )
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     gaux(nl(:)) = gaux(nl(:)) + g(ipol,:) * &
          CMPLX( -AIMAG( aux(nl(:)) ) , DBLE( aux(nl(:)) ) )
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     gaux(nlm(:)) = CMPLX( DBLE( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) )
     !
  END IF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  CALL cft3( gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  da(:) = tpiba * DBLE( gaux(:) )
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE paw_grad_dot

