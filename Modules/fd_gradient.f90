!
! Copyright (C) 2006-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!
!=----------------------------------------------------------------------=!
   MODULE fd_gradient
!=----------------------------------------------------------------------=!
   !! Module to compute finite differences gradients on dense real space grid.  
   !! Written by Oliviero Andreussi.


      USE kinds, ONLY: DP

      IMPLICIT NONE

  CONTAINS
!=----------------------------------------------------------------------=!
SUBROUTINE calc_fd_gradient( nfdpoint, icfd, ncfd, nnr, f, grad )
!=----------------------------------------------------------------------=!
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg, alat
  USE fft_base,      ONLY : dfftp
  USE fft_types,     ONLY : fft_index_to_3d
  USE scatter_mod,   ONLY : scatter_grid
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nfdpoint
  INTEGER, INTENT(IN)  :: ncfd
  INTEGER, INTENT(IN)  :: icfd(-nfdpoint:nfdpoint)
  INTEGER, INTENT(IN)  :: nnr
  REAL( DP ), DIMENSION( nnr ), INTENT(IN) :: f
  REAL( DP ), DIMENSION( 3, nnr ), INTENT(OUT) :: grad

  INTEGER :: i, ir, ir_end, ipol, in
  INTEGER :: ix(-nfdpoint:nfdpoint),iy(-nfdpoint:nfdpoint),iz(-nfdpoint:nfdpoint)
  INTEGER :: ixc, iyc, izc, ixp, ixm, iyp, iym, izp, izm
  REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gradtmp, gradaux
  LOGICAL  :: offrange
  !
  grad = 0.D0
  !
  ALLOCATE( gradtmp( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3 ) )
  gradtmp = 0.D0
  !
#if defined (__MPI)
  ir_end = MIN(nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
  ir_end = nnr
#endif
  !
  DO ir = 1, ir_end
    !
    CALL fft_index_to_3d (ir, dfftp, ix(0), iy(0), iz(0), offrange)
    IF ( offrange ) CYCLE
    !
    DO in = 1, nfdpoint
      ix(in) = ix(in-1) + 1
      IF( ix(in) .GT. dfftp%nr1-1 ) ix(in) = 0
      ix(-in) = ix(-in+1) - 1
      IF( ix(-in) .LT. 0 ) ix(-in) = dfftp%nr1-1
      iy(in) = iy(in-1) + 1
      IF( iy(in) .GT. dfftp%nr2-1 ) iy(in) = 0
      iy(-in) = iy(-in+1) - 1
      IF( iy(-in) .LT. 0 ) iy(-in) = dfftp%nr2-1
      iz(in) = iz(in-1) + 1
      IF( iz(in) .GT. dfftp%nr3-1 ) iz(in) = 0
      iz(-in) = iz(-in+1) - 1
      IF( iz(-in) .LT. 0 ) iz(-in) = dfftp%nr3-1
    ENDDO
    !
    DO in = -nfdpoint, nfdpoint
      i = ix(in) + iy(0) * dfftp%nr1x + iz(0) * dfftp%nr1x * dfftp%nr2x + 1
      gradtmp(i,1) = gradtmp(i,1) - icfd(in)*f(ir)*dfftp%nr1
      i = ix(0) + iy(in) * dfftp%nr1x + iz(0) * dfftp%nr1x * dfftp%nr2x + 1
      gradtmp(i,2) = gradtmp(i,2) - icfd(in)*f(ir)*dfftp%nr2
      i = ix(0) + iy(0) * dfftp%nr1x + iz(in) * dfftp%nr1x * dfftp%nr2x + 1
      gradtmp(i,3) = gradtmp(i,3) - icfd(in)*f(ir)*dfftp%nr3
    ENDDO
    !
  ENDDO
  !
  ALLOCATE( gradaux(nnr,3) )
#if defined (__MPI)
  DO ipol = 1, 3
    CALL mp_sum( gradtmp(:,ipol), intra_bgrp_comm )
    CALL scatter_grid ( dfftp, gradtmp(:,ipol), gradaux(:,ipol) )
  ENDDO
#else
  gradaux(1:nnr,:) = gradtmp(1:nnr,:)
#endif
  !
  DEALLOCATE( gradtmp )
  !
  DO ir = 1,nnr
    grad(:,ir) = MATMUL( bg, gradaux(ir,:) )
  ENDDO
  DEALLOCATE( gradaux )
  !
  grad = grad / DBLE(ncfd) / alat
  !
  RETURN

END SUBROUTINE calc_fd_gradient

SUBROUTINE init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )

  USE kinds,         ONLY : DP

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: ifdtype, nfdpoint
  INTEGER, INTENT(OUT) :: ncfd
  INTEGER, INTENT(OUT) :: icfd(-nfdpoint:nfdpoint)
  !
  INTEGER :: in
  !
  ncfd = 0
  icfd = 0
  !
  SELECT CASE ( ifdtype )
    !
  CASE ( 1 )
    ! (2N+1)-point Central Differences
    IF ( nfdpoint .EQ. 1 ) THEN
      ncfd = 2
      icfd(  1 ) =   1
    ELSE IF ( nfdpoint .EQ. 2 ) THEN
      ncfd = 12
      icfd(  2 ) =  -1
      icfd(  1 ) =   8
    ELSE IF ( nfdpoint .EQ. 3 ) THEN
      ncfd = 60
      icfd(  3 ) =   1
      icfd(  2 ) =  -9
      icfd(  1 ) =  45
    ELSE IF ( nfdpoint .EQ. 4 ) THEN
      ncfd = 840
      icfd(  4 ) =  -3
      icfd(  3 ) =  32
      icfd(  2 ) =-168
      icfd(  1 ) = 672
    ELSE
      WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
      STOP
    ENDIF
    !
  CASE ( 2 )
    ! Low-Noise Lanczos Differentiators ( M = 2 )
    IF ( nfdpoint .GE. 2 ) THEN
      ncfd = (nfdpoint)*(nfdpoint+1)*(2*nfdpoint+1)/3
      DO in = 1,nfdpoint
        icfd( in ) = in
      ENDDO
    ELSE
      WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
      STOP
    END IF
    !
  CASE ( 3 )
    ! Super Lanczos Low-Noise Differentiators ( M = 4 )
    IF ( nfdpoint .EQ. 3 ) THEN
      ncfd = 252
      icfd(  3 ) = -22
      icfd(  2 ) =  67
      icfd(  1 ) =  58
    ELSE IF ( nfdpoint .EQ. 4 ) THEN
      ncfd = 1188
      icfd(  4 ) = -86
      icfd(  3 ) = 142
      icfd(  2 ) = 193
      icfd(  1 ) = 126
    ELSE IF ( nfdpoint .EQ. 5 ) THEN
      ncfd = 5148
      icfd(  5 ) =-300
      icfd(  4 ) = 294
      icfd(  3 ) = 532
      icfd(  2 ) = 503
      icfd(  1 ) = 296
    ELSE
      WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
      STOP
    ENDIF
    !
  CASE ( 4 )
    ! Smooth Noise-Robust Differentiators  ( n = 2 )
    IF ( nfdpoint .EQ. 2 ) THEN
      ncfd = 8
      icfd(  2 ) =   1
      icfd(  1 ) =   2
    ELSE IF ( nfdpoint .EQ. 3 ) THEN
      ncfd = 32
      icfd(  3 ) =   1
      icfd(  2 ) =   4
      icfd(  1 ) =   5
    ELSE IF ( nfdpoint .EQ. 4 ) THEN
      ncfd = 128
      icfd(  4 ) =   1
      icfd(  3 ) =   6
      icfd(  2 ) =  14
      icfd(  1 ) =  14
    ELSE IF ( nfdpoint .EQ. 5 ) THEN
      ncfd = 512
      icfd(  5 ) =   1
      icfd(  4 ) =   8
      icfd(  3 ) =  27
      icfd(  2 ) =  48
      icfd(  1 ) =  42
    ELSE
      WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
      STOP
    ENDIF
    !
  CASE ( 5 )
    ! Smooth Noise-Robust Differentiators  ( n = 4 )
    IF ( nfdpoint .EQ. 3 ) THEN
      ncfd = 96
      icfd(  3 ) =  -5
      icfd(  2 ) =  12
      icfd(  1 ) =  39
    ELSE IF ( nfdpoint .EQ. 4 ) THEN
      ncfd = 96
      icfd(  4 ) =  -2
      icfd(  3 ) =  -1
      icfd(  2 ) =  16
      icfd(  1 ) =  27
    ELSE IF ( nfdpoint .EQ. 5 ) THEN
      ncfd = 1536
      icfd(  5 ) = -11
      icfd(  4 ) = -32
      icfd(  3 ) =  39
      icfd(  2 ) = 256
      icfd(  1 ) = 322
    ELSE
      WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
      STOP
    ENDIF
    !
  CASE DEFAULT
    !
    WRITE(*,*)'ERROR: finite difference type unknown, ifdtype=',ifdtype
    STOP
    !
  END SELECT
  !
  DO in = 1,nfdpoint
    icfd( -in ) = - icfd( in )
  ENDDO
  !
  RETURN
  !
END SUBROUTINE init_fd_gradient

!=----------------------------------------------------------------------=!
   END MODULE fd_gradient
!=----------------------------------------------------------------------=!
