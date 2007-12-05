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
  USE control_flags, ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, ONLY :  DP
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
    CALL pw_gemm ('Y', nkb, m, npw, beta, npwx, psi, npwx, betapsi, nkb)
    !
    RETURN
    !
  END SUBROUTINE calbec_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
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
    CALL ccalbec( nkb, npwx, npw, m, betapsi, beta, psi )
    !
    RETURN
    !
  END SUBROUTINE calbec_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
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
    CALL ccalbec_nc( nkb, npwx, npw, npol, m, betapsi, beta, psi )
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
