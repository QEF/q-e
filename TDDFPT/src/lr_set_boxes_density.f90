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
  USE lr_variables,         ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ir, nnr, ir_end, idx, j0, k0
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_set_boxes_density>")')
  ENDIF
  !
  CALL start_clock( 'lr_set_boxes' )
  !
  nnr = dfftp%nnr
  !
  ALLOCATE( cube_save(nnr,3) )
  cube_save = 0
  !
#if defined (__MPI)
  j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
  j0 = 0 ; k0 = 0
  ir_end = dfftp%nnr
#endif  
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     idx = ir -1
     k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
     idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
     k   = k + k0
     IF ( k .GE. dfftp%nr3 ) CYCLE
     j   = idx / dfftp%nr1x
     idx = idx - dfftp%nr1x * j
     j   = j + j0
     IF ( j .GE. dfftp%nr2 ) CYCLE
     i   = idx
     IF ( i .GE. dfftp%nr1 ) CYCLE
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
