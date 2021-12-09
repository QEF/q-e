!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS1D(arr) lbound(arr,1):ubound(arr,1)
#define DIMS2D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2)
#define DIMS3D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3)
#define DIMS4D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3),lbound(arr,4):ubound(arr,4)
#define DIMS5D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3),lbound(arr,4):ubound(arr,4),lbound(arr,5):ubound(arr,5)
!=----------------------------------------------------------------------------=!
   MODULE scf_gpum
!=----------------------------------------------------------------------------=!
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
     INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
     INTEGER, PARAMETER :: i4b = selected_int_kind(9)
     INTEGER, PARAMETER :: i8b = selected_int_kind(18)
#if defined(__DEBUG)
     INTEGER :: iverbosity = 1
#else  
     INTEGER :: iverbosity = 0
#endif
     !
     REAL(DP), ALLOCATABLE :: vrs_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: vrs_d
#endif

     LOGICAL :: vrs_ood = .false.    ! used to flag out of date variables
     LOGICAL :: vrs_d_ood = .false.    ! used to flag out of date variables
     !
     CONTAINS
     !
     SUBROUTINE using_vrs(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE scf, ONLY : vrs
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_vrs ", debug_info, vrs_ood
         !
         IF (vrs_ood) THEN
             IF ((.not. allocated(vrs_d)) .and. (intento_ < 2)) THEN
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
                IF ( iverbosity > 0 ) print *, "Really copied vrs D->H"
                vrs = vrs_d
             END IF
             vrs_ood = .false.
         ENDIF
         IF (intento_ > 0)    vrs_d_ood = .true.
#endif
     END SUBROUTINE using_vrs
     !
     SUBROUTINE using_vrs_d(intento, debug_info)
         !
         USE scf, ONLY : vrs
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_vrs_d ", debug_info, vrs_d_ood
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
             IF (.not. allocated(vrs_d)) ALLOCATE(vrs_d(DIMS2D(vrs)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied vrs H->D"
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
       vrs_d_ood = .false.
     END SUBROUTINE deallocate_scf_gpu
!=----------------------------------------------------------------------------=!
   END MODULE scf_gpum
!=----------------------------------------------------------------------------=!
