!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_vsgga( rhoout, grho, vsgga )
  !----------------------------------------------------------------------------
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl, ngm, g
  USE cell_base,            ONLY : alat
  USE noncollin_module,     ONLY : noncolin, nspin_gga
  USE funct,                ONLY : gcxc, gcx_spin, gcc_spin, &
                                   gcc_spin_more, dft_is_gradient, get_igcc
  USE spin_orb,             ONLY : domag
  USE fft_base,             ONLY : dfftp
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rhoout(dfftp%nnr,nspin_gga)
  REAL(DP),    INTENT(IN)    :: grho(3,dfftp%nnr,nspin_gga)
  REAL(DP),    INTENT(OUT)   :: vsgga(dfftp%nnr)
  !
  INTEGER :: k, ipol, is
  !
  REAL(DP),    ALLOCATABLE :: h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: vaux(:,:)
  !
  LOGICAL  :: igcc_is_lyp
  REAL(DP) :: grho2(2), sx, sc, v2c, &
              v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw , &
              arho, zeta, rh, grh2
  REAL(DP) :: v2cup, v2cdw,  v2cud, rup, rdw, &
              grhoup, grhodw, grhoud, grup, grdw
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-6, &
                         vanishing_mag    = 1.D-12
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  
  IF ( .NOT. (noncolin.and.domag) ) &
     call errore('compute_vsgga','routine called in the wrong case',1)

  igcc_is_lyp = (get_igcc() == 3)
  !
  ALLOCATE(    h( 3, dfftp%nnr, nspin_gga) )
  ALLOCATE( vaux( dfftp%nnr, nspin_gga ) )

  DO k = 1, dfftp%nnr
     !
     rh = rhoout(k,1) + rhoout(k,2)
     !
     arho=abs(rh)
     !
     IF ( arho > vanishing_charge ) THEN
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        IF ( grho2(1) > epsg .OR. grho2(2) > epsg ) THEN
           CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
                          grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
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
        ELSE
           !
           sc    = 0.D0
           sx    = 0.D0
           v1xup = 0.D0
           v1xdw = 0.D0
           v2xup = 0.D0
           v2xdw = 0.D0
           v1cup = 0.D0
           v1cdw = 0.D0
           v2c   = 0.D0
           v2cup = 0.D0
           v2cdw = 0.D0
           v2cud = 0.D0
        ENDIF
     ELSE
        !
        sc    = 0.D0
        sx    = 0.D0
        v1xup = 0.D0
        v1xdw = 0.D0
        v2xup = 0.D0
        v2xdw = 0.D0
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
     vaux(k,1) = e2 * ( v1xup + v1cup )
     vaux(k,2) = e2 * ( v1xdw + v1cdw )
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
  END DO
  !
  ALLOCATE( dh( dfftp%nnr ) )
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin_gga
     !
     CALL grad_dot( dfftp%nnr, h(1,1,is), ngm, g, nl, alat, dh )
     !
     vaux(:,is) = vaux(:,is) - dh(:)
     !
  END DO

  vsgga(:)=(vaux(:,1)-vaux(:,2))

  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE compute_vsgga
!
