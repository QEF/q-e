!
! Copyright (C) 2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
! Various routines computing gradient and similar quantities via FFT
!--------------------------------------------------------------------
! FIXME: there is a dependency upon "cell_base" via variable tpiba
!        (2\pi/a) that maybe should be taken out from here?
!--------------------------------------------------------------------
SUBROUTINE external_gradient( a, grada )
  !--------------------------------------------------------------------
  !! Interface for computing gradients of a real function in real space,
  !! to be called by an external module.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  !! a real function on the real-space
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  !! grad(A) in real space
  !
  ! A in real space, grad(A) in real space
  CALL fft_gradient_r2r( dfftp, a, g, grada )
  !
  RETURN
  !
END SUBROUTINE external_gradient
!
!----------------------------------------------------------------------------
SUBROUTINE fft_gradient_r2r( dfft, a, g, ga )
  !----------------------------------------------------------------------------
  !! Calculates \({\bf ga}\), the gradient of \({\bf a}\).
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba
  USE fft_interfaces,  ONLY : fwfft, invfft 
  USE fft_types,       ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  REAL(DP), INTENT(IN)  :: a(dfft%nnr)
  !! a real function on the real-space FFT grid
  REAL(DP), INTENT(IN)  :: g(3,dfft%ngm)
  !! G-vectors, in 2\pi/a units
  REAL(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !! gradient of a, real, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( gaux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = (0.0_dp, 0.0_dp)
     !
     gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
                                             REAL( aux(dfft%nl(:)) ), kind=DP)
     !
     IF ( dfft%lgamma ) THEN
        !
        gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
                                  -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
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
END SUBROUTINE fft_gradient_r2r
!
!--------------------------------------------------------------------
SUBROUTINE fft_qgradient( dfft, a, xq, g, ga )
  !--------------------------------------------------------------------
  !! Like \texttt{fft\_gradient\_r2r}, for complex arrays having a 
  !! \(e^{iqr}\) behavior.
  !
  USE kinds,     ONLY: dp
  USE cell_base, ONLY: tpiba
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY: fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  COMPLEX(DP), INTENT(IN)  :: a(dfft%nnr)
  !! A complex function on the real-space FFT grid
  REAL(DP), INTENT(IN):: xq(3)
  !! q-vector, in 2\pi/a units
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  !! G-vectors, in 2\pi/a units
  COMPLEX(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !! The gradient of \({\bf a}\), complex, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER  :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)

  ALLOCATE (gaux(dfft%nnr))
  ALLOCATE (aux (dfft%nnr))

  ! bring a(r) to G-space, a(G) ...
  aux (:) = a(:)

  CALL fwfft ('Rho', aux, dfft)
  ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  DO ipol = 1, 3
     gaux (:) = (0.0_dp, 0.0_dp)
     DO n = 1, dfft%ngm
        gaux(dfft%nl(n)) = CMPLX( 0.0_dp, xq (ipol) + g(ipol,n), kind=DP ) * &
             aux (dfft%nl(n))
        IF ( dfft%lgamma ) gaux(dfft%nlm(n)) = CONJG( gaux (dfft%nl(n)) )
     END DO
     ! bring back to R-space, (\grad_ipol a)(r) ...

     CALL invfft ('Rho', gaux, dfft)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     DO n = 1, dfft%nnr
        ga (ipol,n) = gaux (n) * tpiba
     END DO
  END DO

  DEALLOCATE (aux)
  DEALLOCATE (gaux)

  RETURN

END SUBROUTINE fft_qgradient
!
!
!----------------------------------------------------------------------------
SUBROUTINE fft_gradient_g2r( dfft, a, g, ga )
  !----------------------------------------------------------------------------
  !! Calculates \({\bf ga}\), the gradient of \({\bf a}\) - like
  !! \(\textrm{fft_gradient}\), but with \({\bf a(G)}\) instead of \({\bf a(r)}\).
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  COMPLEX(DP), INTENT(IN) :: a(dfft%ngm)
  !! a(G), a complex function in G-space
  REAL(DP), INTENT(IN)  :: g(3,dfft%ngm)
  !! G-vectors, in \( 2\pi/a \) units
  REAL(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !! The gradient of \({\bf a}\), real, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER :: ipol, n, ip
  COMPLEX(DP), ALLOCATABLE :: gaux(:)
  !
#if defined(__CUDA) && defined(_OPENACC)
  INTEGER, POINTER, DEVICE :: nl_d(:), nlm_d(:)
  !
  nl_d  => dfft%nl_d
  nlm_d => dfft%nlm_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:), nlm_d(:)
  !
  ALLOCATE( nl_d(dfft%ngm) )
  nl_d = dfft%nl
  IF ( dfft%lgamma ) THEN
    ALLOCATE( nlm_d(dfft%ngm) )
    nlm_d = dfft%nlm
  ENDIF
  !$acc data copyin( nl_d, nlm_d )
#endif
  !
  !$acc data present_or_copyin( a, g ) present_or_copyout( ga )
  !
  ALLOCATE( gaux(dfft%nnr) )
  !$acc data create( gaux )
  !
  IF ( dfft%lgamma ) THEN
     !
     ! ... Gamma tricks: perform 2 FFT's in a single shot
     ! x and y
     ipol = 1
     !$acc parallel loop
     DO n = 1, dfft%nnr
       DO ip = 1, 3
         ga(ip,n) = 0.D0
       ENDDO
       gaux(n) = (0.0_dp,0.0_dp)
     ENDDO
     !
     ! ... multiply a(G) by iG to get the gradient in real space
     !
     !$acc parallel loop
     DO n = 1, dfft%ngm
        gaux(nl_d(n) ) = CMPLX( 0.0_dp, g(ipol,n),kind=DP)* a(n) - &
                                      CMPLX(g(ipol+1,n),kind=DP) * a(n)
        gaux(nlm_d(n)) = CMPLX( 0.0_dp,-g(ipol,n),kind=DP)*CONJG(a(n)) + &
                                      CMPLX(g(ipol+1,n),kind=DP) * CONJG(a(n))
     ENDDO
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     !$acc host_data use_device( gaux )
     CALL invfft( 'Rho', gaux, dfft )
     !$acc end host_data
     !
     ! ... bring back to R-space, (\grad_ipol a)(r)
     ! ... add the factor 2\pi/a  missing in the definition of q+G
     !
     !$acc parallel loop
     DO n = 1, dfft%nnr
        ga(ipol,n) = REAL( gaux(n), kind=DP ) * tpiba
        ga(ipol+1,n) = AIMAG( gaux(n) ) * tpiba
     ENDDO
     ! ... for z
     ipol = 3
     !$acc parallel loop
     DO n = 1, dfft%nnr
        gaux(n) = (0.0_dp,0.0_dp)
     ENDDO
     !
     ! ... multiply a(G) by iG to get the gradient in real space
     !
     !$acc parallel loop
     DO n = 1, dfft%ngm
        gaux(nl_d(n)) = CMPLX(g(ipol,n),kind=DP) * &
                        CMPLX( -AIMAG(a(n)), REAL(a(n)),kind=DP)
        gaux(nlm_d(n)) = CONJG( gaux(nl_d(n)) )
     ENDDO
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     !$acc host_data use_device( gaux )
     CALL invfft( 'Rho', gaux, dfft )
     !$acc end host_data
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     !$acc parallel loop
     DO n = 1, dfft%nnr
       ga(ipol,n) = tpiba * REAL( gaux(n), kind=DP )
     ENDDO
     !
  ELSE
     !
     DO ipol = 1, 3
        !
        !$acc parallel loop
        DO n = 1, dfft%nnr
           ga(ipol,n) = 0.D0
           gaux(n) = (0.0_dp,0.0_dp)
        ENDDO
        !
        ! ... multiply a(G) by iG to get the gradient in real space
        !
        !$acc parallel loop
        DO n = 1, dfft%ngm
          gaux(nl_d(n)) = CMPLX(g(ipol,n), kind=DP) * &
                          CMPLX( -AIMAG(a(n)), REAL(a(n)), kind=DP)
        ENDDO
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        !$acc host_data use_device( gaux )
        CALL invfft( 'Rho', gaux, dfft )
        !$acc end host_data
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        !$acc parallel loop
        DO n = 1, dfft%nnr
          ga(ipol,n) = tpiba * REAL( gaux(n), kind=DP )
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( gaux )
  !
  !$acc end data
  !
#if !defined(__CUDA) || !defined(_OPENACC)
  !$acc end data
  DEALLOCATE( nl_d )
  IF ( dfft%lgamma ) DEALLOCATE( nlm_d )
#endif
  !
  RETURN
  !
END SUBROUTINE fft_gradient_g2r
!
!
!----------------------------------------------------------------------------
SUBROUTINE fft_graddot( dfft, a, g, da )
  !---------------------------------------------------------------------------
  !! Calculates \( da = \sum_i \nabla_i a_i \) in R-space.
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  !! FFT descriptor
  REAL(DP), INTENT(IN) :: a(3,dfft%nnr)
  !! A real function on the real-space FFT grid
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  !! G-vectors, in \( 2\pi/a \) units
  REAL(DP), INTENT(OUT) :: da(dfft%nnr)
  !! \( \sum_i \nabla_i a_i \), real, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  COMPLEX(DP) :: fp, fm, aux1, aux2
  !
#if defined(__CUDA) && defined(_OPENACC)
  INTEGER, POINTER, DEVICE :: nl_d(:), nlm_d(:)
  !
  nl_d  => dfft%nl_d
  nlm_d => dfft%nlm_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:), nlm_d(:)
  !
  ALLOCATE( nl_d(dfft%ngm) )
  nl_d = dfft%nl
  IF ( dfft%lgamma ) THEN
    ALLOCATE( nlm_d(dfft%ngm) )
    nlm_d = dfft%nlm
  ENDIF
  !$acc data copyin( nl_d, nlm_d )
#endif
  !
  !$acc data present_or_copyin( a, g ) present_or_copyout( da )
  !
  ALLOCATE( aux(dfft%nnr), gaux(dfft%nnr) )
  !$acc data create( aux, gaux )
  !
  !$acc parallel loop
  DO n = 1, dfft%nnr
     gaux(n) = (0.0_dp,0.0_dp)
  ENDDO
  !
  IF ( dfft%lgamma ) THEN
     !
     ! Gamma tricks: perform 2 FFT's in a single shot
     ! x and y
     ipol = 1
     !
     !$acc parallel loop
     DO n = 1, dfft%nnr
       aux(n) = CMPLX( a(ipol,n), a(ipol+1,n), kind=DP)
     ENDDO
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     !$acc host_data use_device( aux )
     CALL fwfft( 'Rho', aux, dfft )
     !$acc end host_data
     !
     ! ... multiply by iG to get the gradient in G-space
     !
     !$acc parallel loop
     DO n = 1, dfft%ngm
        fp = (aux(nl_d(n)) + aux(nlm_d(n)))*0.5_dp
        fm = (aux(nl_d(n)) - aux(nlm_d(n)))*0.5_dp
        aux1 = CMPLX( REAL(fp), AIMAG(fm), kind=DP)
        aux2 = CMPLX(AIMAG(fp), -REAL(fm), kind=DP)
        gaux(nl_d(n)) = &
             CMPLX(0.0_dp, g(ipol  ,n),kind=DP) * aux1 + &
             CMPLX(0.0_dp, g(ipol+1,n),kind=DP) * aux2
     ENDDO
     ! ... for z
     ipol = 3
     !$acc parallel loop
     DO n = 1, dfft%nnr
       aux(n) = CMPLX( a(ipol,n), 0.0_dp, kind=DP)
     ENDDO
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     !$acc host_data use_device( aux )
     CALL fwfft( 'Rho', aux, dfft )
     !$acc end host_data
     !
     ! ... multiply by iG to get the gradient in G-space
     ! ... fill both gaux(G) and gaux(-G) = gaux*(G)
     !
     !$acc parallel loop
     DO n = 1, dfft%ngm
        gaux(nl_d(n)) = gaux(nl_d(n)) + CMPLX(g(ipol,n),kind=DP) * &
             CMPLX( -AIMAG( aux(nl_d(n)) ), &
                      REAL( aux(nl_d(n)) ), kind=DP)
        gaux(nlm_d(n)) = CONJG( gaux(nl_d(n)) )
     ENDDO
     !
  ELSE
     !
     DO ipol = 1, 3
        !
        !$acc parallel loop
        DO n = 1, dfft%nnr
          aux(n) = CMPLX( a(ipol,n), 0.0_dp, kind=DP)
        ENDDO
        !
        ! ... bring a(ipol,r) to G-space, a(G) ...
        !
        !$acc host_data use_device( aux )
        CALL fwfft( 'Rho', aux, dfft )
        !$acc end host_data
        !
        ! ... multiply by iG to get the gradient in G-space
        !
        !$acc parallel loop
        DO n = 1, dfft%ngm
           gaux(nl_d(n)) = gaux(nl_d(n)) + CMPLX(g(ipol,n),kind=DP) * &
                CMPLX( -AIMAG( aux(nl_d(n)) ), &
                         REAL( aux(nl_d(n)),kind=DP ), kind=DP)
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  !$acc host_data use_device( gaux )
  CALL invfft( 'Rho', gaux, dfft )
  !$acc end host_data
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  !$acc parallel loop
  DO n = 1, dfft%nnr
    da(n) = tpiba * REAL( gaux(n), kind=DP )
  ENDDO
  !
  !$acc end data
  DEALLOCATE( aux, gaux )
  !
  !$acc end data
  !
#if !defined(__CUDA) || !defined(_OPENACC)
  !$acc end data
  DEALLOCATE( nl_d )
  IF ( dfft%lgamma ) DEALLOCATE( nlm_d )
#endif  
  !
  RETURN
  !
END SUBROUTINE fft_graddot
!
!
!--------------------------------------------------------------------
SUBROUTINE fft_qgraddot ( dfft, a, xq, g, da)
  !--------------------------------------------------------------------
  !! Like \(\textrm{fft_graddot}\), for complex arrays having a \(e^{iqr}\)
  !! dependency.
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba
  USE fft_interfaces, ONLY : fwfft, invfft
  USE fft_types,      ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  !! FFT descriptor
  COMPLEX(DP), INTENT(IN)  :: a(3,dfft%nnr)
  !! a complex function on the real-space FFT grid
  REAL(DP), INTENT(IN) :: xq(3)
  !! q-vector, in \( 2\pi/a \) units
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  !! G-vectors, in \( 2\pi/a \) units
  COMPLEX(DP), INTENT(OUT) :: da(dfft%nnr)
  !! \( \sum_i \nabla_i a_i \), complex, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER :: n, ipol
  COMPLEX(DP), allocatable :: aux (:)

  ALLOCATE (aux (dfft%nnr))
  da(:) = (0.0_dp, 0.0_dp)
  DO ipol = 1, 3
     ! copy a(ipol,r) to a complex array...
     DO n = 1, dfft%nnr
        aux (n) = a (ipol, n)
     END DO
     ! bring a(ipol,r) to G-space, a(G) ...
     CALL fwfft ('Rho', aux, dfft)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     DO n = 1, dfft%ngm
        da (dfft%nl(n)) = da (dfft%nl(n)) + &
             CMPLX(0.0_dp, xq (ipol) + g (ipol, n),kind=DP) * aux(dfft%nl(n))
     END DO
  END DO
  IF ( dfft%lgamma ) THEN
     DO n = 1, dfft%ngm
        da (dfft%nlm(n)) = CONJG( da (dfft%nl(n)) )
     END DO
  END IF
  !  bring back to R-space, (\grad_ipol a)(r) ...
  CALL invfft ('Rho', da, dfft)
  ! ...add the factor 2\pi/a  missing in the definition of q+G
  da (:) = da (:) * tpiba

  DEALLOCATE(aux)

  RETURN
 
END SUBROUTINE fft_qgraddot

!--------------------------------------------------------------------
! Routines computing laplacian via FFT
!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE external_laplacian( a, lapla )
  !--------------------------------------------------------------------
  !! Interface for computing laplacian in real space, to be called by 
  !! an external module.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : gg
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: a( dfftp%nnr )
  !! A real function on the real-space FFT grid
  REAL(DP), INTENT(OUT) :: lapla( dfftp%nnr )
  !! Laplacian of \( {\bf a} \) in real space.
  !
  ! A in real space, lapl(A) in real space
  CALL fft_laplacian( dfftp, a, gg, lapla )

  RETURN

END SUBROUTINE external_laplacian

!--------------------------------------------------------------------
SUBROUTINE fft_laplacian( dfft, a, gg, lapla )
  !--------------------------------------------------------------------
  !! Calculates \(\text{lapla} = \nabla^2(a)\).
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba2
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  REAL(DP), INTENT(IN) :: a(dfft%nnr)
  !! A real function on the real-space FFT grid
  REAL(DP), INTENT(IN) :: gg(dfft%ngm)
  !! Square modules of G-vectors, in \( (2\pi/a)^2 \) units
  REAL(DP), INTENT(OUT) :: lapla(dfft%nnr)
  !! \(\text{lapla}(:)=\nabla^2(a)\), real, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
  !
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( laux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... Compute the laplacian
  !
  laux(:) = (0.0_dp, 0.0_dp)
  !
  DO ig = 1, dfft%ngm
     !
     laux(dfft%nl(ig)) = -gg(ig)*aux(dfft%nl(ig))
     !
  END DO
  !
  IF ( dfft%lgamma ) THEN
     !
     laux(dfft%nlm(:)) = CMPLX( REAL(laux(dfft%nl(:)) ), &
                              -AIMAG(laux(dfft%nl(:)) ), kind=DP)
     !
  ENDIF
  !
  ! ... bring back to R-space, (\lapl a)(r) ...
  !
  CALL invfft ('Rho', laux, dfft)
  !
  ! ... add the missing factor (2\pi/a)^2 in G
  !
  lapla = tpiba2 * REAL( laux )
  !
  DEALLOCATE( laux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_laplacian
!
!--------------------------------------------------------------------
! Routines computing hessian via FFT
!--------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE fft_hessian_g2r ( dfft, a, g, ha )
  !----------------------------------------------------------------------
  !! Calculates \( \text{ha} = \text{hessian}(a) \)
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : tpiba
  USE fft_types,              ONLY : fft_type_descriptor
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE fft_helper_subroutines, ONLY : fftx_oned2threed
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  !! G-vectors, in \( (2\pi/a)^2 \) units.
  COMPLEX(DP), INTENT(IN) :: a(dfft%ngm)
  !! A real function on the real-space FFT grid.
  REAL(DP), INTENT(OUT) :: ha( 6, dfft%nnr )
  !! \( \text{Hessian}(a) \), real, on the real-space FFT grid.  
  !! Lower-packed matrix indeces 1-6 correspond to:
  !! \(1=xx\), \(2=yx\), \(3=yy\), \(4=zx\), \(5=zy\), \(6=zz\).
  !
  ! ... local variables
  !
  INTEGER                  :: ig, ir
  COMPLEX(DP), ALLOCATABLE :: aux(:), haux(:,:)
  !
  IF ( .NOT. dfft%lgamma ) CALL errore ('fft_hessian_g2r',&
       'only gamma case is implemented',1)
  ALLOCATE ( aux(dfft%nnr))
  ALLOCATE (haux(dfft%ngm,2))
  ! xx, yx
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(1,ig)**2     *a(ig)
     haux(ig,2) = -tpiba**2*g(1,ig)*g(2,ig)*a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(1,ir) = DBLE(aux(ir))
     ha(2,ir) =AIMAG(aux(ir))
  END DO
  ! yy, zx
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(2,ig)**2     *a(ig)
     haux(ig,2) = -tpiba**2*g(1,ig)*g(3,ig)*a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(3,ir) = DBLE(aux(ir))
     ha(4,ir) =AIMAG(aux(ir))
  END DO
  ! zy, zz
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(2,ig)*g(3,ig)*a(ig)
     haux(ig,2) = -tpiba**2*g(3,ig)**2     *a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(5,ir) = DBLE(aux(ir))
     ha(6,ir) =AIMAG(aux(ir))
  END DO
  !
  DEALLOCATE(aux)
  DEALLOCATE(haux)
  
END SUBROUTINE fft_hessian_g2r
!--------------------------------------------------------------------
SUBROUTINE fft_hessian( dfft, a, g, ga, ha )
  !--------------------------------------------------------------------
  !! Calculates \( \text{ga} = \nabla a\) and \(\text{ha} = \text{hessian}(a)\)
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !! FFT descriptor
  REAL(DP), INTENT(IN) :: a(dfft%nnr)
  !! A real function on the real-space FFT grid
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  !! G-vectors, in \( (2\pi/a)^2 \) units
  REAL(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !! \(\nabla a\), real, on the real-space FFT grid
  REAL(DP), INTENT(OUT) :: ha(3,3,dfft%nnr)
  !! \(text{hessian}(a)\), real, on the real-space FFT grid
  !
  ! ... local variables
  !
  INTEGER :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
  !
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( gaux( dfft%nnr ) )
  ALLOCATE( haux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = (0.0_dp,0.0_dp)
     !
     gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
                                             REAL( aux(dfft%nl(:)) ), kind=DP )
     !
     IF ( dfft%lgamma ) THEN
        !
        gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
                                  -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * REAL( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        haux(:) = (0.0_dp,0.0_dp)
        !
        haux(dfft%nl(:)) = - g(ipol,:) * g(jpol,:) * &
             CMPLX( REAL( aux(dfft%nl(:)) ), &
                   AIMAG( aux(dfft%nl(:)) ), kind=DP)
        !
        IF ( dfft%lgamma ) THEN
           !
           haux(dfft%nlm(:)) = CMPLX(  REAL( haux(dfft%nl(:)) ), &
                                     -AIMAG( haux(dfft%nl(:)) ), kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Rho', haux, dfft)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        ha(ipol, jpol, :) = tpiba * tpiba * REAL( haux(:) )
        !
        ha(jpol, ipol, :) = ha(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( haux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_hessian
!--------------------------------------------------------------------
SUBROUTINE external_hessian( a, grada, hessa )
  !--------------------------------------------------------------------
  !! Interface for computing hessian in real space, to be called by
  !! an external module.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)   :: a( dfftp%nnr )
  !! A real function on the real-space FFT grid.
  REAL(DP), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  !! gradient in real space.
  REAL(DP), INTENT(OUT)  :: hessa( 3, 3, dfftp%nnr )
  !! Hessian in real space.
  !
  ! A in real space, grad(A) and hess(A) in real space
  CALL fft_hessian( dfftp, a, g, grada, hessa )

  RETURN

END SUBROUTINE external_hessian

!--------------------------------------------------------------------
