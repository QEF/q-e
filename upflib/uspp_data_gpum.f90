!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS3D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3)
#define DIMS4D(arr) lbound(arr,1):ubound(arr,1),lbound(arr,2):ubound(arr,2),lbound(arr,3):ubound(arr,3),lbound(arr,4):ubound(arr,4)
!=----------------------------------------------------------------------------=!
   MODULE uspp_data_gpum
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
     INTEGER :: iverbosity = 0
#if defined(__DEBUG)
     iverbosity = 1
#endif
     !
     LOGICAL :: qrad_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qrad_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_at_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_at_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d2y_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d2y_d_ood = .false.    ! used to flag out of date variables
     !
     CONTAINS
     !
     SUBROUTINE using_qrad(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE uspp_data, ONLY : qrad, qrad_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_qrad ", debug_info, qrad_ood
         !
         IF (qrad_ood) THEN
             IF ((.not. allocated(qrad_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_qrad_d', 'PANIC: sync of qrad from qrad_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(qrad)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of qrad with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    qrad_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied qrad D->H"
                qrad = qrad_d
             END IF
             qrad_ood = .false.
         ENDIF
         IF (intento_ > 0)    qrad_d_ood = .true.
#endif
     END SUBROUTINE using_qrad
     !
     SUBROUTINE using_qrad_d(intento, debug_info)
         !
         USE uspp_data, ONLY : qrad, qrad_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_qrad_d ", debug_info, qrad_d_ood
         !
         IF (.not. allocated(qrad)) THEN
             IF (intento /= 2) print *, "WARNING: sync of qrad_d with unallocated array and intento /= 2?"
             IF (allocated(qrad_d)) DEALLOCATE(qrad_d)
             qrad_d_ood = .false.
             RETURN
         END IF
         ! here we know that qrad is allocated, check if size is 0
         IF ( SIZE(qrad) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array qrad_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (qrad_d_ood) THEN
             IF ( allocated(qrad_d) .and. (SIZE(qrad_d)/=SIZE(qrad))) deallocate(qrad_d)
             IF (.not. allocated(qrad_d)) ALLOCATE(qrad_d(DIMS4D(qrad)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied qrad H->D"
                qrad_d = qrad
             END IF
             qrad_d_ood = .false.
         ENDIF
         IF (intento > 0)    qrad_ood = .true.
#else
         CALL errore('using_qrad_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_qrad_d
     !
     SUBROUTINE using_tab(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE uspp_data, ONLY : tab, tab_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab ", debug_info, tab_ood
         !
         IF (tab_ood) THEN
             IF ((.not. allocated(tab_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_tab_d', 'PANIC: sync of tab from tab_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(tab)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of tab with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    tab_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab D->H"
                tab = tab_d
             END IF
             tab_ood = .false.
         ENDIF
         IF (intento_ > 0)    tab_d_ood = .true.
#endif
     END SUBROUTINE using_tab
     !
     SUBROUTINE using_tab_d(intento, debug_info)
         !
         USE uspp_data, ONLY : tab, tab_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_d ", debug_info, tab_d_ood
         !
         IF (.not. allocated(tab)) THEN
             IF (intento /= 2) print *, "WARNING: sync of tab_d with unallocated array and intento /= 2?"
             IF (allocated(tab_d)) DEALLOCATE(tab_d)
             tab_d_ood = .false.
             RETURN
         END IF
         ! here we know that tab is allocated, check if size is 0
         IF ( SIZE(tab) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array tab_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (tab_d_ood) THEN
             IF ( allocated(tab_d) .and. (SIZE(tab_d)/=SIZE(tab))) deallocate(tab_d)
             IF (.not. allocated(tab_d)) ALLOCATE(tab_d(DIMS3D(tab)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab H->D"
                tab_d = tab
             END IF
             tab_d_ood = .false.
         ENDIF
         IF (intento > 0)    tab_ood = .true.
#else
         CALL errore('using_tab_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_tab_d
     !
     SUBROUTINE using_tab_at(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE uspp_data, ONLY : tab_at, tab_at_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_at ", debug_info, tab_at_ood
         !
         IF (tab_at_ood) THEN
             IF ((.not. allocated(tab_at_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_tab_at_d', 'PANIC: sync of tab_at from tab_at_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(tab_at)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of tab_at with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    tab_at_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab_at D->H"
                tab_at = tab_at_d
             END IF
             tab_at_ood = .false.
         ENDIF
         IF (intento_ > 0)    tab_at_d_ood = .true.
#endif
     END SUBROUTINE using_tab_at
     !
     SUBROUTINE using_tab_at_d(intento, debug_info)
         !
         USE uspp_data, ONLY : tab_at, tab_at_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_at_d ", debug_info, tab_at_d_ood
         !
         IF (.not. allocated(tab_at)) THEN
             IF (intento /= 2) print *, "WARNING: sync of tab_at_d with unallocated array and intento /= 2?"
             IF (allocated(tab_at_d)) DEALLOCATE(tab_at_d)
             tab_at_d_ood = .false.
             RETURN
         END IF
         ! here we know that tab_at is allocated, check if size is 0
         IF ( SIZE(tab_at) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array tab_at_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (tab_at_d_ood) THEN
             IF ( allocated(tab_at_d) .and. (SIZE(tab_at_d)/=SIZE(tab_at))) deallocate(tab_at_d)
             IF (.not. allocated(tab_at_d)) ALLOCATE(tab_at_d(DIMS3D(tab_at)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab_at H->D"
                tab_at_d = tab_at
             END IF
             tab_at_d_ood = .false.
         ENDIF
         IF (intento > 0)    tab_at_ood = .true.
#else
         CALL errore('using_tab_at_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_tab_at_d
     !
     SUBROUTINE using_tab_d2y(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE uspp_data, ONLY : tab_d2y, tab_d2y_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)  || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_d2y ", debug_info, tab_d2y_ood
         !
         IF (tab_d2y_ood) THEN
             IF ((.not. allocated(tab_d2y_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_tab_d2y_d', 'PANIC: sync of tab_d2y from tab_d2y_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(tab_d2y)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of tab_d2y with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    tab_d2y_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab_d2y D->H"
                tab_d2y = tab_d2y_d
             END IF
             tab_d2y_ood = .false.
         ENDIF
         IF (intento_ > 0)    tab_d2y_d_ood = .true.
#endif
     END SUBROUTINE using_tab_d2y
     !
     SUBROUTINE using_tab_d2y_d(intento, debug_info)
         !
         USE uspp_data, ONLY : tab_d2y, tab_d2y_d
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_d2y_d ", debug_info, tab_d2y_d_ood
         !
         IF (.not. allocated(tab_d2y)) THEN
             IF (intento /= 2) print *, "WARNING: sync of tab_d2y_d with unallocated array and intento /= 2?"
             IF (allocated(tab_d2y_d)) DEALLOCATE(tab_d2y_d)
             tab_d2y_d_ood = .false.
             RETURN
         END IF
         ! here we know that tab_d2y is allocated, check if size is 0
         IF ( SIZE(tab_d2y) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array tab_d2y_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (tab_d2y_d_ood) THEN
             IF ( allocated(tab_d2y_d) .and. (SIZE(tab_d2y_d)/=SIZE(tab_d2y))) deallocate(tab_d2y_d)
             IF (.not. allocated(tab_d2y_d)) ALLOCATE(tab_d2y_d(DIMS3D(tab_d2y)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied tab_d2y H->D"
                tab_d2y_d = tab_d2y
             END IF
             tab_d2y_d_ood = .false.
         ENDIF
         IF (intento > 0)    tab_d2y_ood = .true.
#else
         CALL errore('using_tab_d2y_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_tab_d2y_d
     !
     SUBROUTINE deallocate_us_gpu
       qrad_d_ood = .false.
       tab_d_ood = .false.
       tab_at_d_ood = .false.
       tab_d2y_d_ood = .false.
     END SUBROUTINE deallocate_us_gpu
!=----------------------------------------------------------------------------=!
   END MODULE uspp_data_gpum
!=----------------------------------------------------------------------------=!!

