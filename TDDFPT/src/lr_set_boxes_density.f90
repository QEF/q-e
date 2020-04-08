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
  USE fft_types,            ONLY : fft_index_to_3d
  USE lr_variables,         ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ir, ir_end
  LOGICAL :: offrange
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_set_boxes_density>")')
  ENDIF
  !
  CALL start_clock( 'lr_set_boxes' )
  !
  ALLOCATE( cube_save(dfftp%nnr,3) )
  cube_save = 0
  !
#if defined (__MPI)
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
  ir_end = dfftp%nnr
#endif  
  !
  DO ir = 1, ir_end
     !
     ! ... three dimensional indexes
     !
     CALL fft_index_to_3d (ir, dfftp, i,j,k, offrange)
     IF ( offrange ) CYCLE
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
