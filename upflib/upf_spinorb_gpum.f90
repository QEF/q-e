!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS5D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3),lbound(arr,4):ubound(arr,4),lbound(arr,5):ubound(arr,5)
!=----------------------------------------------------------------------------=!
   MODULE upf_spinorb_gpum
!=----------------------------------------------------------------------------=!
#if defined(__CUDA)
     USE cudafor
#endif
     USE upf_spinorb,  ONLY: fcoef_d
     IMPLICIT NONE
     SAVE
     !
     INTEGER :: iverbosity = 0
#if defined(__DEBUG)
     iverbosity = 1
#endif
     !
     LOGICAL :: fcoef_ood = .false.    ! used to flag out of date variables
     LOGICAL :: fcoef_d_ood = .false.    ! used to flag out of date variables
     !
     CONTAINS
     !
     SUBROUTINE using_fcoef(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE upf_spinorb, ONLY : fcoef, fcoef_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_fcoef ", debug_info, fcoef_ood
         !
         IF (fcoef_ood) THEN
             IF ((.not. allocated(fcoef_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_fcoef_d', 'PANIC: sync of fcoef from fcoef_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(fcoef)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of fcoef with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    fcoef_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied fcoef D->H"
                fcoef = fcoef_d
             END IF
             fcoef_ood = .false.
         ENDIF
         IF (intento_ > 0)    fcoef_d_ood = .true.
#endif
     END SUBROUTINE using_fcoef
     !
     SUBROUTINE using_fcoef_d(intento, debug_info)
         !
         USE upf_spinorb, ONLY : fcoef, fcoef_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_fcoef_d ", debug_info, fcoef_d_ood
         !
         IF (.not. allocated(fcoef)) THEN
             IF (intento /= 2) print *, "WARNING: sync of fcoef_d with unallocated array and intento /= 2?"
             IF (allocated(fcoef_d)) DEALLOCATE(fcoef_d)
             fcoef_d_ood = .false.
             RETURN
         END IF
         ! here we know that fcoef is allocated, check if size is 0
         IF ( SIZE(fcoef) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array fcoef_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (fcoef_d_ood) THEN
             IF ( allocated(fcoef_d) .and. (SIZE(fcoef_d)/=SIZE(fcoef))) deallocate(fcoef_d)
             IF (.not. allocated(fcoef_d)) ALLOCATE(fcoef_d(DIMS5D(fcoef)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied fcoef H->D"
                fcoef_d = fcoef
             END IF
             fcoef_d_ood = .false.
         ENDIF
         IF (intento > 0)    fcoef_ood = .true.
#else
         CALL errore('using_fcoef_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_fcoef_d
     !
     SUBROUTINE deallocate_spinorb_gpu
       fcoef_d_ood = .false.
     END SUBROUTINE deallocate_spinorb_gpu
!=----------------------------------------------------------------------------=!
   END MODULE upf_spinorb_gpum
!=----------------------------------------------------------------------------=!
