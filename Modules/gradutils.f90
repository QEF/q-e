!
! Copyright (C) 2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
! Routines computing gradient via FFT
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
SUBROUTINE external_gradient( a, grada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradients of a real function in real space,
  ! to be called by an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL ffT_gradient( dfftp, a, g, grada )

  RETURN

END SUBROUTINE external_gradient
!----------------------------------------------------------------------------
SUBROUTINE fft_gradient( dfft, a, g, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         g(3,:)   G-vectors, in 2pi/a units
  ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft, invfft 
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: a(dfft%nnr), g(3,dfft%ngm)
  REAL(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !
  INTEGER  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( gaux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.0_dp, 0.0_dp)
     !
     gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
                                             REAL( aux(dfft%nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
                                  -AIMAG( gaux(dfft%nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_gradient
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! Routines computing laplacian via FFT
!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE external_laplacian( a, lapla )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing laplacian in real space, to be called by 
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : gg
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: lapla( dfftp%nnr )

! A in real space, lapl(A) in real space
  CALL fft_laplacian( dfftp, a, gg, lapla )

  RETURN

END SUBROUTINE external_laplacian

!--------------------------------------------------------------------
SUBROUTINE fft_laplacian( dfft, a, gg, lapla )
!--------------------------------------------------------------------
  !
  ! ... Calculates lapla = laplacian(a)
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         gg(:)    square modules of G-vectors, in (2pi/a)^2 units
  ! ... output: lapla(:) \nabla^2 a, real, on the real-space FFT grid
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba2
  USE control_flags, ONLY : gamma_only
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: a(dfft%nnr), gg(dfft%ngm)
  REAL(DP), INTENT(OUT) :: lapla(dfft%nnr)
  !
  INTEGER                  :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
  !
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( laux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... Compute the laplacian
  !
  laux(:) = CMPLX(0.0_dp, 0.0_dp)
  !
  DO ig = 1, dfft%ngm
     !
     laux(dfft%nl(ig)) = -gg(ig)*aux(dfft%nl(ig))
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     laux(dfft%nlm(:)) = CMPLX( REAL(laux(dfft%nl(:)) ), &
                              -AIMAG(laux(dfft%nl(:)) ) ,kind=DP)
     !
  ENDIF
  !
  ! ... bring back to R-space, (\lapl a)(r) ...
  !
  CALL invfft ('Rho', laux, dfft)
  !
  ! ... add the missing factor (2\pi/a)^2 in G
  !
  lapla = tpiba2 * DBLE( laux )
  !
  DEALLOCATE( laux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_laplacian
