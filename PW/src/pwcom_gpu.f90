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
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE :: g2kin_d(:)
     REAL(DP), ALLOCATABLE :: et_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: g2kin_d, et_d

     LOGICAL :: g2kin_ood = .false.    ! used to flag out of date variables
     LOGICAL :: g2kin_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: et_ood = .false.    ! used to flag out of date variables
     LOGICAL :: et_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_g2kin(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wvfct, ONLY : g2kin
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_g2kin ", debug_info, g2kin_ood
         !
         IF (g2kin_ood) THEN
             IF (.not. allocated(g2kin_d)) THEN
                CALL errore('using_g2kin_d', 'PANIC: sync of g2kin from g2kin_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(g2kin)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of g2kin with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    g2kin_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied g2kin D->H"
                g2kin = g2kin_d
             END IF
             g2kin_ood = .false.
         ENDIF
         IF (intento_ > 0)    g2kin_d_ood = .true.
#endif
     END SUBROUTINE using_g2kin
     !
     SUBROUTINE using_g2kin_d(intento, debug_info)
         !
         USE wvfct, ONLY : g2kin
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_g2kin_d ", debug_info, g2kin_d_ood
         !
         IF (.not. allocated(g2kin)) THEN
             IF (intento /= 2) print *, "WARNING: sync of g2kin_d with unallocated array and intento /= 2?"
             IF (allocated(g2kin_d)) DEALLOCATE(g2kin_d)
             g2kin_d_ood = .false.
             RETURN
         END IF
         ! here we know that g2kin is allocated, check if size is 0 
         IF ( SIZE(g2kin) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array g2kin_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (g2kin_d_ood) THEN
             IF ( allocated(g2kin_d) .and. (SIZE(g2kin_d)/=SIZE(g2kin))) deallocate(g2kin_d)
             IF (.not. allocated(g2kin_d)) ALLOCATE(g2kin_d(DIMS1D(g2kin)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied g2kin H->D"
                g2kin_d = g2kin
             END IF
             g2kin_d_ood = .false.
         ENDIF
         IF (intento > 0)    g2kin_ood = .true.
#else
         CALL errore('using_g2kin_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_g2kin_d
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
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_et ", debug_info, et_ood
         !
         IF (et_ood) THEN
             IF (.not. allocated(et_d)) THEN
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
                print *, "Really copied et D->H"
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
#if defined(__CUDA)
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
                print *, "Really copied et H->D"
                et_d = et
             END IF
             et_d_ood = .false.
         ENDIF
         IF (intento > 0)    et_ood = .true.
#else
         CALL errore('using_et_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_et_d
     !
     SUBROUTINE deallocate_wvfct_gpu
       IF( ALLOCATED( g2kin_d ) ) DEALLOCATE( g2kin_d )
       IF( ALLOCATED( et_d ) ) DEALLOCATE( et_d )
     END SUBROUTINE deallocate_wvfct_gpu
!=----------------------------------------------------------------------------=!
   END MODULE wvfct_gpum
!=----------------------------------------------------------------------------=!!
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
   MODULE us_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE :: qrad_d(:, :, :, :)
     REAL(DP), ALLOCATABLE :: tab_d(:, :, :)
     REAL(DP), ALLOCATABLE :: tab_at_d(:, :, :)
     REAL(DP), ALLOCATABLE :: tab_d2y_d(:, :, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: qrad_d, tab_d, tab_at_d, tab_d2y_d

     LOGICAL :: qrad_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qrad_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_at_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_at_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d2y_ood = .false.    ! used to flag out of date variables
     LOGICAL :: tab_d2y_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_qrad(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE us, ONLY : qrad
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_qrad ", debug_info, qrad_ood
         !
         IF (qrad_ood) THEN
             IF (.not. allocated(qrad_d)) THEN
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
                print *, "Really copied qrad D->H"
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
         USE us, ONLY : qrad
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
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
                print *, "Really copied qrad H->D"
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
         USE us, ONLY : tab
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab ", debug_info, tab_ood
         !
         IF (tab_ood) THEN
             IF (.not. allocated(tab_d)) THEN
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
                print *, "Really copied tab D->H"
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
         USE us, ONLY : tab
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
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
                print *, "Really copied tab H->D"
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
         USE us, ONLY : tab_at
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_at ", debug_info, tab_at_ood
         !
         IF (tab_at_ood) THEN
             IF (.not. allocated(tab_at_d)) THEN
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
                print *, "Really copied tab_at D->H"
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
         USE us, ONLY : tab_at
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
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
                print *, "Really copied tab_at H->D"
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
         USE us, ONLY : tab_d2y
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_tab_d2y ", debug_info, tab_d2y_ood
         !
         IF (tab_d2y_ood) THEN
             IF (.not. allocated(tab_d2y_d)) THEN
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
                print *, "Really copied tab_d2y D->H"
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
         USE us, ONLY : tab_d2y
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
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
                print *, "Really copied tab_d2y H->D"
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
       IF( ALLOCATED( qrad_d ) ) DEALLOCATE( qrad_d )
       IF( ALLOCATED( tab_d ) ) DEALLOCATE( tab_d )
       IF( ALLOCATED( tab_at_d ) ) DEALLOCATE( tab_at_d )
       IF( ALLOCATED( tab_d2y_d ) ) DEALLOCATE( tab_d2y_d )
     END SUBROUTINE deallocate_us_gpu
!=----------------------------------------------------------------------------=!
   END MODULE us_gpum
!=----------------------------------------------------------------------------=!!
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
   MODULE spin_orb_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     COMPLEX(DP), ALLOCATABLE :: fcoef_d(:, :, :, :, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: fcoef_d

     LOGICAL :: fcoef_ood = .false.    ! used to flag out of date variables
     LOGICAL :: fcoef_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_fcoef(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE spin_orb, ONLY : fcoef
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_fcoef ", debug_info, fcoef_ood
         !
         IF (fcoef_ood) THEN
             IF (.not. allocated(fcoef_d)) THEN
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
                print *, "Really copied fcoef D->H"
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
         USE spin_orb, ONLY : fcoef
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
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
                print *, "Really copied fcoef H->D"
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
     SUBROUTINE deallocate_spin_orb_gpu
       IF( ALLOCATED( fcoef_d ) ) DEALLOCATE( fcoef_d )
     END SUBROUTINE deallocate_spin_orb_gpu
!=----------------------------------------------------------------------------=!
   END MODULE spin_orb_gpum
!=----------------------------------------------------------------------------=!