!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod
  !
  ! ... *bec* contain <beta|psi> - used in h_psi, s_psi, many other places
  ! ... calbec( npw, beta, psi, betapsi [, nbnd ] ) is an interface calculating
  ! ...    betapsi(i,j)  = <beta(i)|psi(j)>   (the sum is over npw components)
  ! ... or betapsi(i,s,j)= <beta(i)|psi(s,j)> (s=polarization index)
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE gvect,            ONLY : gstart
  USE noncollin_module, ONLY : noncolin, npol
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: &
       rbecp(:,:)       !   <beta|psi> for real (at Gamma) wavefunctions 
  COMPLEX(DP), ALLOCATABLE ::  &
       becp (:,:), &    !  as above for complex wavefunctions
       becp_nc(:,:,:)   !  as above for spinors
  !
  INTERFACE calbec
     !
     MODULE PROCEDURE calbec_k, calbec_gamma, calbec_nc
     !
  END INTERFACE
  !
  PRIVATE
  PUBLIC :: rbecp, becp, becp_nc, allocate_bec, deallocate_bec, calbec
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on 
    ! ... half of the G-vectors or PWs - assuming k=0 is the G=0 component:
    ! ... betapsi(i,j) = 2Re(\sum_k beta^*(i,k)psi(k,j)) + beta^*(i,0)psi(0,j)
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (IN) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (OUT) :: betapsi(:,:)
    INTEGER, INTENT (IN) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
    nkb = SIZE (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    npwx= SIZE (beta, 1)
    IF ( npwx /= SIZE (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( PRESENT (nbnd) ) THEN
        m = nbnd
    ELSE
        m = SIZE ( psi, 2)
    END IF
    IF ( nkb /= SIZE (betapsi,1) .OR. m > SIZE (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
        !
        CALL DGEMV( 'C', 2*npw, nkb, 2.0_DP, beta, 2*npwx, psi, 1, 0.0_DP, &
                     betapsi, 1 )
        IF ( gstart == 2 ) betapsi(:,1) = betapsi(:,1) - beta(1,:)*psi(1,1)
        !
    ELSE
        !
        CALL DGEMM( 'C', 'N', nkb, m, 2*npw, 2.0_DP, beta, 2*npwx, psi, &
                    2*npwx, 0.0_DP, betapsi, nkb )
        IF ( gstart == 2 ) &
           CALL DGER( nkb, m, -1.0_DP, beta, 2*npwx, psi, 2*npwx, betapsi, nkb )
        !
     END IF
     !
     CALL reduce( nkb*m, betapsi )
     !
     CALL stop_clock( 'calbec' )
     !
    RETURN
    !
  END SUBROUTINE calbec_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on 
    ! ... G-vectors or PWs : betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (IN) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (OUT) :: betapsi(:,:)
    INTEGER, INTENT (IN) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
    nkb = SIZE (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    npwx= SIZE (beta, 1)
    IF ( npwx /= SIZE (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( PRESENT (nbnd) ) THEN
        m = nbnd
    ELSE
        m = SIZE ( psi, 2)
    END IF
    IF ( nkb /= SIZE (betapsi,1) .OR. m > SIZE (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
       !
       CALL ZGEMV( 'C', npw, nkb, (1.0_DP,0.0_DP), beta, npwx, psi, 1, &
                   (0.0_DP, 0.0_DP), betapsi, 1 )
       !
    ELSE
       !
       CALL ZGEMM( 'C', 'N', nkb, m, npw, (1.0_DP,0.0_DP), &
                 beta, npwx, psi, npwx, (0.0_DP,0.0_DP), betapsi, nkb )
       !
    END IF
    !
    CALL reduce( 2*nkb*m, betapsi )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k below) running on 
    ! ... G-vectors or PWs corresponding to two different polarizations:
    ! ... betapsi(i,1,j) = \sum_k=1,npw beta^*(i,k) psi(k,j)
    ! ... betapsi(i,2,j) = \sum_k=1,npw beta^*(i,k) psi(k+npwx,j)
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (IN) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (OUT) :: betapsi(:,:,:)
    INTEGER, INTENT (IN) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, npol, m
    !
    nkb = SIZE (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    call start_clock ('calbec')
    npwx= SIZE (beta, 1)
    IF ( 2*npwx /= SIZE (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( PRESENT (nbnd) ) THEN
        m = nbnd
    ELSE
        m = SIZE ( psi, 2)
    END IF
    npol= SIZE (betapsi, 2)
    IF ( nkb /= SIZE (betapsi,1) .OR. m > SIZE (betapsi, 3) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    call ZGEMM ('C', 'N', nkb, m*npol, npw, (1.0_DP, 0.0_DP), beta, &
              npwx, psi, npwx, (0.0_DP, 0.0_DP),  betapsi, nkb)
    !
    CALL reduce( 2*nkb*m*npol, betapsi )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_nc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec ( nkb, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nkb, nbnd
    !
    IF ( gamma_only ) THEN 
       !
       ALLOCATE( rbecp( nkb, nbnd ) )
       !
    ELSE IF ( noncolin) THEN
       !
       ALLOCATE( becp_nc( nkb, npol, nbnd ) )
       !
    ELSE
       !
       ALLOCATE( becp( nkb, nbnd ) )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE allocate_bec
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec ()
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF ( gamma_only ) THEN 
       !
       DEALLOCATE( rbecp )
       !
    ELSE IF ( noncolin) THEN
       !
       DEALLOCATE( becp_nc )
       !
    ELSE
       !
       DEALLOCATE( becp )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec
  !
END MODULE becmod
