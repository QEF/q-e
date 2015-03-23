!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_set_boxes_density()
  !---------------------------------------------------------------------
  !
  ! Set boxes for the calculation of the charge density response.
  ! Inspired by Modules/compute_dipole.f90
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2009)
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE lr_variables,         ONLY : cube_save
  USE fft_base,             ONLY : dfftp
  USE mp_global,            ONLY : me_bgrp
  USE lr_variables,         ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ir, nnr, ir_end, index0
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_set_boxes_density>")')
  ENDIF
  !
  CALL start_clock( 'lr_set_boxes' )
  !
  ALLOCATE( cube_save( dfftp%nnr, 3 ) )
  cube_save = 0
  !
  nnr = dfftp%nnr
  !
#if defined (__MPI)
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%npp(1:me_bgrp))
  ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  index0 = 0
  ir_end = nnr
#endif
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
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
  RETURN
  !
END SUBROUTINE lr_set_boxes_density
