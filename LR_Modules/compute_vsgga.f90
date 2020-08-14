!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
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
  USE gvect,                ONLY : ngm, g
  USE noncollin_module,     ONLY : noncolin, nspin_gga
  USE funct,                ONLY : dft_is_gradient, get_igcc
  USE xc_gga,               ONLY : xc_gcx
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
  REAL(DP), ALLOCATABLE :: h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: vaux(:,:)
  !
  LOGICAL  :: igcc_is_lyp
  REAL(DP), DIMENSION(dfftp%nnr) :: sx, sc, v2cud
  REAL(DP), DIMENSION(dfftp%nnr,2) :: v1x, v2x, v1c, v2c
  REAL(DP) :: gr(2)
  !
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
  ALLOCATE( h(3,dfftp%nnr,nspin_gga)  )
  ALLOCATE( vaux(dfftp%nnr,nspin_gga) )
  
  CALL xc_gcx( dfftp%nnr, 2, rhoout, grho, sx, sc, v1x, v2x, v1c, v2c, v2cud )
     
  DO k = 1, dfftp%nnr
     !
     ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
     !
     vaux(k,1) = e2 * ( v1x(k,1) + v1c(k,1) )
     vaux(k,2) = e2 * ( v1x(k,2) + v1c(k,2) )
     !
     ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
     !
     DO ipol = 1, 3
        !
        gr(:) = grho(ipol,k,:)
        h(ipol,k,1) = e2 * ( ( v2x(k,1) + v2c(k,1) ) * gr(1) + v2cud(k) * gr(2) )
        h(ipol,k,2) = e2 * ( ( v2x(k,2) + v2c(k,2) ) * gr(2) + v2cud(k) * gr(1) )
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
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
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
