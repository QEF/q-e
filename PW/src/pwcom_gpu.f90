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
   MODULE wvfct_gpum
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
     REAL(DP), ALLOCATABLE :: et_d(:, :)
     REAL(DP), ALLOCATABLE :: wg_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: et_d, wg_d
#endif

     LOGICAL :: et_ood = .false.    ! used to flag out of date variables
     LOGICAL :: et_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: wg_ood = .false.    ! used to flag out of date variables
     LOGICAL :: wg_d_ood = .false.    ! used to flag out of date variables
     !
     CONTAINS
     !
     SUBROUTINE using_et(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wvfct, ONLY : et
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_et ", debug_info, et_ood
         !
         IF (et_ood) THEN
             IF ((.not. allocated(et_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_et_d', 'PANIC: sync of et from et_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(et)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of et with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    et_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied et D->H"
                et = et_d
             END IF
             et_ood = .false.
         ENDIF
         IF (intento_ > 0)    et_d_ood = .true.
#endif
     END SUBROUTINE using_et
     !
     SUBROUTINE using_et_d(intento, debug_info)
         !
         USE wvfct, ONLY : et
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_et_d ", debug_info, et_d_ood
         !
         IF (.not. allocated(et)) THEN
             IF (intento /= 2) print *, "WARNING: sync of et_d with unallocated array and intento /= 2?"
             IF (allocated(et_d)) DEALLOCATE(et_d)
             et_d_ood = .false.
             RETURN
         END IF
         ! here we know that et is allocated, check if size is 0
         IF ( SIZE(et) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array et_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (et_d_ood) THEN
             IF ( allocated(et_d) .and. (SIZE(et_d)/=SIZE(et))) deallocate(et_d)
             IF (.not. allocated(et_d)) ALLOCATE(et_d(DIMS2D(et)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied et H->D"
                et_d = et
             END IF
             et_d_ood = .false.
         ENDIF
         IF (intento > 0)    et_ood = .true.
#else
         CALL errore('using_et_d', 'Trying to use device data without device compiled code!', 1)
#endif
     END SUBROUTINE using_et_d
     !
     SUBROUTINE using_wg(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wvfct, ONLY : wg
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_wg ", debug_info, wg_ood
         !
         IF (wg_ood) THEN
             IF ((.not. allocated(wg_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_wg_d', 'PANIC: sync of wg from wg_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(wg)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of wg with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    wg_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied wg D->H"
                wg = wg_d
             END IF
             wg_ood = .false.
         ENDIF
         IF (intento_ > 0)    wg_d_ood = .true.
#endif
     END SUBROUTINE using_wg
     !
     SUBROUTINE using_wg_d(intento, debug_info)
         !
         USE wvfct, ONLY : wg
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_wg_d ", debug_info, wg_d_ood
         !
         IF (.not. allocated(wg)) THEN
             IF (intento /= 2) print *, "WARNING: sync of wg_d with unallocated array and intento /= 2?"
             IF (allocated(wg_d)) DEALLOCATE(wg_d)
             wg_d_ood = .false.
             RETURN
         END IF
         ! here we know that wg is allocated, check if size is 0
         IF ( SIZE(wg) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array wg_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (wg_d_ood) THEN
             IF ( allocated(wg_d) .and. (SIZE(wg_d)/=SIZE(wg))) deallocate(wg_d)
             IF (.not. allocated(wg_d)) ALLOCATE(wg_d(DIMS2D(wg)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied wg H->D"
                wg_d = wg
             END IF
             wg_d_ood = .false.
         ENDIF
         IF (intento > 0)    wg_ood = .true.
#else
         CALL errore('using_wg_d', 'Trying to use device data without device compiled code!', 1)
#endif
     END SUBROUTINE using_wg_d
     !
     SUBROUTINE deallocate_wvfct_gpu
       IF( ALLOCATED( et_d ) ) DEALLOCATE( et_d )
       et_d_ood = .false.
       IF( ALLOCATED( wg_d ) ) DEALLOCATE( wg_d )
       wg_d_ood = .false.
     END SUBROUTINE deallocate_wvfct_gpu
!=----------------------------------------------------------------------------=!
   END MODULE wvfct_gpum
!=----------------------------------------------------------------------------=!!
