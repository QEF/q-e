!-----------------------------------------------------------------------
!OBM
! 150608 pfft replaced by fft_base :: dfftp
SUBROUTINE lr_set_boxes_density()
  !---------------------------------------------------------------------
  ! ... set boxes for the calculation of density response
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE lr_variables,         ONLY : cube_save
  !use pfft,                 only : npp
  USE fft_base,              ONLY : dfftp
  USE mp_global,            ONLY : me_pool
  USE lr_variables,   ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  INTEGER :: index0, index, ir
  INTEGER :: i, j, k, p, nr
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_set_boxes_density>")')
  ENDIF
  CALL start_clock( 'lr_set_boxes' )
  !
  ALLOCATE( cube_save( dfftp%nnr, 3 ) )
  cube_save = 0
  !
  index0 = 0
  !
#if defined (__MPI)
  !
  DO i = 1, me_pool
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  ENDDO
  !
#endif
  !
  DO ir = 1, dfftp%nnr
     !
     ! ... three dimensional indexes
     !
     index = index0 + ir - 1
     k     = index / (dfftp%nr1x*dfftp%nr2x)
     index = index - (dfftp%nr1x*dfftp%nr2x)*k
     j     = index / dfftp%nr1x
     index = index - dfftp%nr1x*j
     i     = index
     !
     IF ( i>=dfftp%nr1 .or. j>=dfftp%nr2 .or. k>=dfftp%nr3 ) CYCLE
     !
     cube_save(ir,1) = i
     cube_save(ir,2) = j
     cube_save(ir,3) = k
     !
  ENDDO
  !
  CALL stop_clock( 'lr_set_boxes' )
  !
  END SUBROUTINE lr_set_boxes_density
