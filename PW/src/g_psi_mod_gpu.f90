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
   MODULE g_psi_mod_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE :: h_diag_d(:, :)
     REAL(DP), ALLOCATABLE :: s_diag_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: h_diag_d, s_diag_d

     LOGICAL :: h_diag_ood = .false.    ! used to flag out of date variables
     LOGICAL :: h_diag_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: s_diag_ood = .false.    ! used to flag out of date variables
     LOGICAL :: s_diag_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_h_diag(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE g_psi_mod, ONLY : h_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_h_diag ", debug_info, h_diag_ood
         !
         IF (h_diag_ood) THEN
             IF ((.not. allocated(h_diag_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_h_diag_d', 'PANIC: sync of h_diag from h_diag_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(h_diag)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of h_diag with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    h_diag_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied h_diag D->H"
                h_diag = h_diag_d
             END IF
             h_diag_ood = .false.
         ENDIF
         IF (intento_ > 0)    h_diag_d_ood = .true.
#endif
     END SUBROUTINE using_h_diag
     !
     SUBROUTINE using_h_diag_d(intento, debug_info)
         !
         USE g_psi_mod, ONLY : h_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_h_diag_d ", debug_info, h_diag_d_ood
         !
         IF (.not. allocated(h_diag)) THEN
             IF (intento /= 2) print *, "WARNING: sync of h_diag_d with unallocated array and intento /= 2?"
             IF (allocated(h_diag_d)) DEALLOCATE(h_diag_d)
             h_diag_d_ood = .false.
             RETURN
         END IF
         ! here we know that h_diag is allocated, check if size is 0 
         IF ( SIZE(h_diag) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array h_diag_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (h_diag_d_ood) THEN
             IF ( allocated(h_diag_d) .and. (SIZE(h_diag_d)/=SIZE(h_diag))) deallocate(h_diag_d)
             IF (.not. allocated(h_diag_d)) ALLOCATE(h_diag_d(DIMS2D(h_diag)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied h_diag H->D"
                h_diag_d = h_diag
             END IF
             h_diag_d_ood = .false.
         ENDIF
         IF (intento > 0)    h_diag_ood = .true.
#else
         CALL errore('using_h_diag_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_h_diag_d
     !
     SUBROUTINE using_s_diag(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE g_psi_mod, ONLY : s_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_s_diag ", debug_info, s_diag_ood
         !
         IF (s_diag_ood) THEN
             IF ((.not. allocated(s_diag_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_s_diag_d', 'PANIC: sync of s_diag from s_diag_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(s_diag)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of s_diag with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    s_diag_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied s_diag D->H"
                s_diag = s_diag_d
             END IF
             s_diag_ood = .false.
         ENDIF
         IF (intento_ > 0)    s_diag_d_ood = .true.
#endif
     END SUBROUTINE using_s_diag
     !
     SUBROUTINE using_s_diag_d(intento, debug_info)
         !
         USE g_psi_mod, ONLY : s_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_s_diag_d ", debug_info, s_diag_d_ood
         !
         IF (.not. allocated(s_diag)) THEN
             IF (intento /= 2) print *, "WARNING: sync of s_diag_d with unallocated array and intento /= 2?"
             IF (allocated(s_diag_d)) DEALLOCATE(s_diag_d)
             s_diag_d_ood = .false.
             RETURN
         END IF
         ! here we know that s_diag is allocated, check if size is 0 
         IF ( SIZE(s_diag) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array s_diag_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (s_diag_d_ood) THEN
             IF ( allocated(s_diag_d) .and. (SIZE(s_diag_d)/=SIZE(s_diag))) deallocate(s_diag_d)
             IF (.not. allocated(s_diag_d)) ALLOCATE(s_diag_d(DIMS2D(s_diag)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied s_diag H->D"
                s_diag_d = s_diag
             END IF
             s_diag_d_ood = .false.
         ENDIF
         IF (intento > 0)    s_diag_ood = .true.
#else
         CALL errore('using_s_diag_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_s_diag_d
     !
     SUBROUTINE deallocate_g_psi_mod_gpu
       IF( ALLOCATED( h_diag_d ) ) DEALLOCATE( h_diag_d )
       h_diag_d_ood = .false.
       IF( ALLOCATED( s_diag_d ) ) DEALLOCATE( s_diag_d )
       s_diag_d_ood = .false.
     END SUBROUTINE deallocate_g_psi_mod_gpu
!=----------------------------------------------------------------------------=!
   END MODULE g_psi_mod_gpum
!=----------------------------------------------------------------------------=!