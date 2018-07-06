!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS1D(my_array) lbound(my_array,1):ubound(my_array,1)
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
#define DIMS3D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2),lbound(my_array,3):ubound(my_array,3)
#define DIMS4D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2),lbound(my_array,3):ubound(my_array,3),lbound(my_array,4):ubound(my_array,4)
!=----------------------------------------------------------------------------=!
   MODULE gvect_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE :: gg_d(:)
     REAL(DP), ALLOCATABLE :: g_d(:, :)
     INTEGER, ALLOCATABLE :: mill_d(:, :)
     COMPLEX(DP), ALLOCATABLE :: eigts1_d(:, :)
     COMPLEX(DP), ALLOCATABLE :: eigts2_d(:, :)
     COMPLEX(DP), ALLOCATABLE :: eigts3_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: gg_d, g_d, mill_d, eigts1_d, eigts2_d, eigts3_d

     LOGICAL :: gg_ood = .false.    ! used to flag out of date variables
     LOGICAL :: gg_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: g_ood = .false.    ! used to flag out of date variables
     LOGICAL :: g_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: mill_ood = .false.    ! used to flag out of date variables
     LOGICAL :: mill_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts1_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts1_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts2_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts2_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts3_ood = .false.    ! used to flag out of date variables
     LOGICAL :: eigts3_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_gg(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : gg
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_gg ", debug_info, gg_ood
         !
         IF (gg_ood) THEN
             IF (.not. allocated(gg_d)) THEN
                CALL errore('using_gg_d', 'PANIC: sync of gg from gg_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(gg)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of gg with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    gg_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied gg D->H"
                gg = gg_d
             END IF
             gg_ood = .false.
         ENDIF
         IF (intento_ > 0)    gg_d_ood = .true.
#endif
     END SUBROUTINE using_gg
     !
     SUBROUTINE using_gg_d(intento, debug_info)
         !
         USE gvect, ONLY : gg
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_gg_d ", debug_info, gg_d_ood
         !
         IF (.not. allocated(gg)) THEN
             IF (intento /= 2) print *, "WARNING: sync of gg_d with unallocated array and intento /= 2?"
             IF (allocated(gg_d)) DEALLOCATE(gg_d)
             gg_d_ood = .false.
             RETURN
         END IF
         ! here we know that gg is allocated, check if size is 0 
         IF ( SIZE(gg) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array gg_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (gg_d_ood) THEN
             IF ( allocated(gg_d) .and. (SIZE(gg_d)/=SIZE(gg))) deallocate(gg_d)
             IF (.not. allocated(gg_d)) ALLOCATE(gg_d(DIMS1D(gg)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied gg H->D"
                gg_d = gg
             END IF
             gg_d_ood = .false.
         ENDIF
         IF (intento > 0)    gg_ood = .true.
#else
         CALL errore('using_gg_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_gg_d
     !
     SUBROUTINE using_g(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : g
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_g ", debug_info, g_ood
         !
         IF (g_ood) THEN
             IF (.not. allocated(g_d)) THEN
                CALL errore('using_g_d', 'PANIC: sync of g from g_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(g)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of g with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    g_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied g D->H"
                g = g_d
             END IF
             g_ood = .false.
         ENDIF
         IF (intento_ > 0)    g_d_ood = .true.
#endif
     END SUBROUTINE using_g
     !
     SUBROUTINE using_g_d(intento, debug_info)
         !
         USE gvect, ONLY : g
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_g_d ", debug_info, g_d_ood
         !
         IF (.not. allocated(g)) THEN
             IF (intento /= 2) print *, "WARNING: sync of g_d with unallocated array and intento /= 2?"
             IF (allocated(g_d)) DEALLOCATE(g_d)
             g_d_ood = .false.
             RETURN
         END IF
         ! here we know that g is allocated, check if size is 0 
         IF ( SIZE(g) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array g_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (g_d_ood) THEN
             IF ( allocated(g_d) .and. (SIZE(g_d)/=SIZE(g))) deallocate(g_d)
             IF (.not. allocated(g_d)) ALLOCATE(g_d(DIMS2D(g)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied g H->D"
                g_d = g
             END IF
             g_d_ood = .false.
         ENDIF
         IF (intento > 0)    g_ood = .true.
#else
         CALL errore('using_g_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_g_d
     !
     SUBROUTINE using_mill(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : mill
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_mill ", debug_info, mill_ood
         !
         IF (mill_ood) THEN
             IF (.not. allocated(mill_d)) THEN
                CALL errore('using_mill_d', 'PANIC: sync of mill from mill_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(mill)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of mill with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    mill_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied mill D->H"
                mill = mill_d
             END IF
             mill_ood = .false.
         ENDIF
         IF (intento_ > 0)    mill_d_ood = .true.
#endif
     END SUBROUTINE using_mill
     !
     SUBROUTINE using_mill_d(intento, debug_info)
         !
         USE gvect, ONLY : mill
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_mill_d ", debug_info, mill_d_ood
         !
         IF (.not. allocated(mill)) THEN
             IF (intento /= 2) print *, "WARNING: sync of mill_d with unallocated array and intento /= 2?"
             IF (allocated(mill_d)) DEALLOCATE(mill_d)
             mill_d_ood = .false.
             RETURN
         END IF
         ! here we know that mill is allocated, check if size is 0 
         IF ( SIZE(mill) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array mill_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (mill_d_ood) THEN
             IF ( allocated(mill_d) .and. (SIZE(mill_d)/=SIZE(mill))) deallocate(mill_d)
             IF (.not. allocated(mill_d)) ALLOCATE(mill_d(DIMS2D(mill)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied mill H->D"
                mill_d = mill
             END IF
             mill_d_ood = .false.
         ENDIF
         IF (intento > 0)    mill_ood = .true.
#else
         CALL errore('using_mill_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_mill_d
     !
     SUBROUTINE using_eigts1(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : eigts1
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts1 ", debug_info, eigts1_ood
         !
         IF (eigts1_ood) THEN
             IF (.not. allocated(eigts1_d)) THEN
                CALL errore('using_eigts1_d', 'PANIC: sync of eigts1 from eigts1_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(eigts1)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of eigts1 with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    eigts1_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied eigts1 D->H"
                eigts1 = eigts1_d
             END IF
             eigts1_ood = .false.
         ENDIF
         IF (intento_ > 0)    eigts1_d_ood = .true.
#endif
     END SUBROUTINE using_eigts1
     !
     SUBROUTINE using_eigts1_d(intento, debug_info)
         !
         USE gvect, ONLY : eigts1
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts1_d ", debug_info, eigts1_d_ood
         !
         IF (.not. allocated(eigts1)) THEN
             IF (intento /= 2) print *, "WARNING: sync of eigts1_d with unallocated array and intento /= 2?"
             IF (allocated(eigts1_d)) DEALLOCATE(eigts1_d)
             eigts1_d_ood = .false.
             RETURN
         END IF
         ! here we know that eigts1 is allocated, check if size is 0 
         IF ( SIZE(eigts1) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array eigts1_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (eigts1_d_ood) THEN
             IF ( allocated(eigts1_d) .and. (SIZE(eigts1_d)/=SIZE(eigts1))) deallocate(eigts1_d)
             IF (.not. allocated(eigts1_d)) ALLOCATE(eigts1_d(DIMS2D(eigts1)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied eigts1 H->D"
                eigts1_d = eigts1
             END IF
             eigts1_d_ood = .false.
         ENDIF
         IF (intento > 0)    eigts1_ood = .true.
#else
         CALL errore('using_eigts1_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_eigts1_d
     !
     SUBROUTINE using_eigts2(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : eigts2
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts2 ", debug_info, eigts2_ood
         !
         IF (eigts2_ood) THEN
             IF (.not. allocated(eigts2_d)) THEN
                CALL errore('using_eigts2_d', 'PANIC: sync of eigts2 from eigts2_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(eigts2)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of eigts2 with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    eigts2_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied eigts2 D->H"
                eigts2 = eigts2_d
             END IF
             eigts2_ood = .false.
         ENDIF
         IF (intento_ > 0)    eigts2_d_ood = .true.
#endif
     END SUBROUTINE using_eigts2
     !
     SUBROUTINE using_eigts2_d(intento, debug_info)
         !
         USE gvect, ONLY : eigts2
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts2_d ", debug_info, eigts2_d_ood
         !
         IF (.not. allocated(eigts2)) THEN
             IF (intento /= 2) print *, "WARNING: sync of eigts2_d with unallocated array and intento /= 2?"
             IF (allocated(eigts2_d)) DEALLOCATE(eigts2_d)
             eigts2_d_ood = .false.
             RETURN
         END IF
         ! here we know that eigts2 is allocated, check if size is 0 
         IF ( SIZE(eigts2) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array eigts2_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (eigts2_d_ood) THEN
             IF ( allocated(eigts2_d) .and. (SIZE(eigts2_d)/=SIZE(eigts2))) deallocate(eigts2_d)
             IF (.not. allocated(eigts2_d)) ALLOCATE(eigts2_d(DIMS2D(eigts2)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied eigts2 H->D"
                eigts2_d = eigts2
             END IF
             eigts2_d_ood = .false.
         ENDIF
         IF (intento > 0)    eigts2_ood = .true.
#else
         CALL errore('using_eigts2_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_eigts2_d
     !
     SUBROUTINE using_eigts3(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE gvect, ONLY : eigts3
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts3 ", debug_info, eigts3_ood
         !
         IF (eigts3_ood) THEN
             IF (.not. allocated(eigts3_d)) THEN
                CALL errore('using_eigts3_d', 'PANIC: sync of eigts3 from eigts3_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(eigts3)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of eigts3 with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    eigts3_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied eigts3 D->H"
                eigts3 = eigts3_d
             END IF
             eigts3_ood = .false.
         ENDIF
         IF (intento_ > 0)    eigts3_d_ood = .true.
#endif
     END SUBROUTINE using_eigts3
     !
     SUBROUTINE using_eigts3_d(intento, debug_info)
         !
         USE gvect, ONLY : eigts3
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_eigts3_d ", debug_info, eigts3_d_ood
         !
         IF (.not. allocated(eigts3)) THEN
             IF (intento /= 2) print *, "WARNING: sync of eigts3_d with unallocated array and intento /= 2?"
             IF (allocated(eigts3_d)) DEALLOCATE(eigts3_d)
             eigts3_d_ood = .false.
             RETURN
         END IF
         ! here we know that eigts3 is allocated, check if size is 0 
         IF ( SIZE(eigts3) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array eigts3_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (eigts3_d_ood) THEN
             IF ( allocated(eigts3_d) .and. (SIZE(eigts3_d)/=SIZE(eigts3))) deallocate(eigts3_d)
             IF (.not. allocated(eigts3_d)) ALLOCATE(eigts3_d(DIMS2D(eigts3)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied eigts3 H->D"
                eigts3_d = eigts3
             END IF
             eigts3_d_ood = .false.
         ENDIF
         IF (intento > 0)    eigts3_ood = .true.
#else
         CALL errore('using_eigts3_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_eigts3_d
     !
     SUBROUTINE deallocate_gvect_gpu
       IF( ALLOCATED( gg_d ) ) DEALLOCATE( gg_d )
       IF( ALLOCATED( g_d ) ) DEALLOCATE( g_d )
       IF( ALLOCATED( mill_d ) ) DEALLOCATE( mill_d )
       IF( ALLOCATED( eigts1_d ) ) DEALLOCATE( eigts1_d )
       IF( ALLOCATED( eigts2_d ) ) DEALLOCATE( eigts2_d )
       IF( ALLOCATED( eigts3_d ) ) DEALLOCATE( eigts3_d )
     END SUBROUTINE deallocate_gvect_gpu
!=----------------------------------------------------------------------------=!
   END MODULE gvect_gpum
!=----------------------------------------------------------------------------=!