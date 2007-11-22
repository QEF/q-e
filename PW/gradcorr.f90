!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define __OLD_NONCOLIN_GGA
!----------------------------------------------------------------------------
SUBROUTINE gradcorr( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   nl, ngm, g
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : omega, alat
  USE funct,                ONLY : gcxc, gcx_spin, gcc_spin, &
                                   gcc_spin_more, dft_is_gradient, get_igcc
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : ux
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rho(nrxx,nspin), rho_core(nrxx)
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm,nspin), rhog_core(ngm)
  REAL(DP),    INTENT(OUT)   :: v(nrxx,nspin)
  REAL(DP),    INTENT(INOUT) :: vtxc, etxc
  !
  INTEGER :: k, ipol, is, nspin0, ir, jpol
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: rhoout(:,:), segni(:), vgg(:,:), vsave(:,:)
  REAL(DP),    ALLOCATABLE :: gmag(:,:,:)

  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  !
  LOGICAL  :: igcc_is_lyp
  REAL(DP) :: grho2(2), sx, sc, v1x, v2x, v1c, v2c, &
              v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw , &
              etxcgc, vtxcgc, segno, arho, fac, zeta, rh, grh2, amag 
  REAL(DP) :: v2cup, v2cdw,  v2cud, rup, rdw, &
              grhoup, grhodw, grhoud, grup, grdw, seg
  !
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN

  igcc_is_lyp = (get_igcc() == 3)
  !
  etxcgc = 0.D0
  vtxcgc = 0.D0
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2
  fac = 1.D0 / DBLE( nspin0 )
  !
  ALLOCATE(    h( 3, nrxx, nspin0) )
  ALLOCATE( grho( 3, nrxx, nspin0) )
  ALLOCATE( rhoout( nrxx, nspin0) )
  IF (nspin==4.AND.domag) THEN
     ALLOCATE( vgg( nrxx, nspin0 ) )
     ALLOCATE( vsave( nrxx, nspin ) )
     ALLOCATE( segni( nrxx ) )
#ifdef __OLD_NONCOLIN_GGA
#else
     ALLOCATE( gmag( 3, nrxx, nspin) )
#endif
     vsave=v
     v=0.d0
  ENDIF
  !
  ALLOCATE( rhogsum( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
#ifdef __OLD_NONCOLIN_GGA
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho(rho,rhoout,segni,nrxx)
     !
     ! ... bring starting rhoout to G-space
     !
     DO is = 1, nspin0
        !
        psic(:) = rhoout(:,is)
        !
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
        !
        rhogsum(:,is) = psic(nl(:))
        !
     END DO
  ELSE
     !
     rhoout(:,1:nspin0)  = rho(:,1:nspin0)
     rhogsum(:,1:nspin0) = rhog(:,1:nspin0)
     !
  ENDIF
  DO is = 1, nspin0
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
     CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                   rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
     !
  END DO
#else
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho_new(rho,rhoout,segni,nrxx,ux)
     !
     rhogsum(:,1) =rhog_core(:) + rhog(:,1)
     !
     CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                   rhogsum, ngm, g, nl, gmag )
     DO is = 2, nspin
        rhogsum(:,1) = rhog(:,is)
        !
        CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                      rhogsum, ngm, g, nl, gmag(1,1,is) )
     END DO
     DO is=1,nspin0
        IF (is==1) seg=0.5_dp
        IF (is==2) seg=-0.5_dp
        DO ipol=1,3
           grho(ipol,:,is) = 0.5_dp*gmag(ipol,:,1)
        ENDDO
        DO ir=1,nrxx
           amag=sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
           rhoout(ir,is) = fac*rho_core(ir) + 0.5_dp*rho(ir,1)
           IF (amag>1.d-12) THEN
              rhoout(ir,is)= rhoout(ir,is) + segni(ir)*seg*amag
              DO ipol=1,3
                 DO jpol=2,4 
                    grho(ipol,ir,is)=grho(ipol,ir,is)+ segni(ir)*seg*rho(ir,jpol)* &
                                              gmag(ipol,ir,jpol)/amag
                 END DO
              END DO
           END IF
        END DO
     END DO
  ELSE
     DO is = 1, nspin0
        !
        rhoout(:,is)  = fac * rho_core(:)  + rho(:,is)
        rhogsum(:,is) = fac * rhog_core(:) + rhog(:,is)
        !
        CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                      rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
        !
     END DO
 END IF
#endif
  !
  DEALLOCATE( rhogsum )
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
              CALL gcxc( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
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
           ELSE
              h(:,k,1)=0.D0
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
        !
        IF ( rh > epsr ) THEN
           !
           IF ( igcc_is_lyp ) THEN
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
#ifdef __OLD_NONCOLIN_GGA
              !
              grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                     ( grho(2,k,1) + grho(2,k,2) )**2 + &
                     ( grho(3,k,1) + grho(3,k,2) )**2
#else
              if (nspin==4) then
                 grh2= gmag(1,k,1)**2+ gmag(2,k,1)**2+gmag(3,k,1)**2
              else
                 grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                        ( grho(2,k,1) + grho(2,k,2) )**2 + &
                        ( grho(3,k,1) + grho(3,k,2) )**2
              endif
#endif
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
     CALL grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
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
#ifdef __OLD_NONCOLIN_GGA
#else
     DEALLOCATE( gmag )
#endif
  ENDIF
  !
  RETURN
  !
END SUBROUTINE gradcorr
!
!----------------------------------------------------------------------------
SUBROUTINE gradrho( nrx1, nrx2, nrx3, &
                    nr1, nr2, nr3, nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is in G-space)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)  :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx
  INTEGER,     INTENT(IN)  :: ngm, nl(ngm)
  COMPLEX(DP), INTENT(IN)  :: a(ngm)
  REAL(DP),    INTENT(IN)  :: g(3,ngm)
  REAL(DP),    INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: gaux(:)
  !
  !
  ALLOCATE( gaux( nrxx ) )
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.D0
  !
  DO ipol = 1, 3
     !
     gaux(:) = 0.D0
     !
     gaux(nl(:)) = g(ipol,:) * CMPLX( -AIMAG( a(:) ), REAL( a(:) ) )
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) )
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL cft3( gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = ga(ipol,:) + tpiba * REAL( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  !
  RETURN
  !
END SUBROUTINE gradrho
!
!----------------------------------------------------------------------------
SUBROUTINE gradient( nrx1, nrx2, nrx3, &
                     nr1, nr2, nr3, nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
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
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) )
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) )
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
END SUBROUTINE gradient
!
!----------------------------------------------------------------------------
SUBROUTINE grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                     nrxx, a, ngm, g, nl, alat, da )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates da = \sum_i \grad_i a_i in R-space
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)     :: nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                              nrxx, ngm, nl(ngm)
  REAL(DP), INTENT(IN)     :: a(3,nrxx), g(3,ngm), alat
  REAL(DP), INTENT(OUT)    :: da(nrxx)
  !
  INTEGER                  :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE( aux( nrxx ), gaux( nrxx ) )
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
     DO n = 1, ngm
        !
        gaux(nl(n)) = gaux(nl(n)) + g(ipol,n) * &
                      CMPLX( -AIMAG( aux(nl(n)) ), REAL( aux(nl(n)) ) )
        !
     END DO
    !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     DO n = 1, ngm
        !
        gaux(nlm(n)) = CONJG( gaux(nl(n)) )
        !
     END DO
     !
  END IF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  CALL cft3( gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  da(:) = tpiba * REAL( gaux(:) )
  !
  DEALLOCATE( aux, gaux )
  !
  RETURN
  !
END SUBROUTINE grad_dot
