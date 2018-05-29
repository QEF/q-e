!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! #define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE scf_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE :: vrs_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: vrs_d

     LOGICAL :: vrs_ood = .false.    ! used to flag out of date variables
     LOGICAL :: vrs_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_vrs(intento)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE scf, ONLY : vrs
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (vrs_ood) THEN
             IF (.not. allocated(vrs_d)) THEN
                CALL errore('using_vrs_d', 'PANIC: sync of vrs from vrs_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(vrs)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of vrs with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF

                ! IF (intento_ > 0)    vrs_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied vrs D->H"
                vrs = vrs_d
             END IF
             vrs_ood = .false.
         ENDIF
         IF (intento_ > 0)    vrs_d_ood = .true.
#endif
     END SUBROUTINE using_vrs
     !
     SUBROUTINE using_vrs_d(intento)
         !
         USE scf, ONLY : vrs
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         !
         IF (.not. allocated(vrs)) THEN
             IF (intento /= 2) print *, "WARNING: sync of vrs_d with unallocated array and intento /= 2?"
             IF (allocated(vrs_d)) DEALLOCATE(vrs_d)
             vrs_d_ood = .false.
             RETURN
         END IF
         ! here we know that vrs is allocated, check if size is 0 
         IF ( SIZE(vrs) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array vrs_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (vrs_d_ood) THEN
             IF ( allocated(vrs_d) .and. (SIZE(vrs_d)/=SIZE(vrs))) deallocate(vrs_d)
             IF (.not. allocated(vrs_d)) ALLOCATE(vrs_d, MOLD=vrs)  ! this copy may be avoided
             IF (intento < 2) THEN
                print *, "Really copied vrs H->D"
                vrs_d = vrs
             END IF
             vrs_d_ood = .false.
         ENDIF
         IF (intento > 0)    vrs_ood = .true.
#else
         CALL errore('using_vrs_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_vrs_d
     !     
     SUBROUTINE deallocate_scf_gpu
       IF( ALLOCATED( vrs_d ) ) DEALLOCATE( vrs_d )
     END SUBROUTINE deallocate_scf_gpu
!=----------------------------------------------------------------------------=!
   END MODULE scf_gpum
!=----------------------------------------------------------------------------=!