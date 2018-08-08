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
   MODULE wavefunctions_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP
#if defined(__CUDA)
     USE cudafor
#endif
     IMPLICIT NONE
     SAVE
     !
     COMPLEX(DP), ALLOCATABLE :: evc_d(:, :)
     COMPLEX(DP), ALLOCATABLE :: psic_d(:)
     COMPLEX(DP), ALLOCATABLE :: psic_nc_d(:, :)
     !
#if defined(__CUDA)
     attributes (DEVICE) :: evc_d, psic_d, psic_nc_d

     LOGICAL :: evc_ood = .false.    ! used to flag out of date variables
     LOGICAL :: evc_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: psic_ood = .false.    ! used to flag out of date variables
     LOGICAL :: psic_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: psic_nc_ood = .false.    ! used to flag out of date variables
     LOGICAL :: psic_nc_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_evc(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wavefunctions, ONLY : evc
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_evc ", debug_info, evc_ood
         !
         IF (evc_ood) THEN
             IF ((.not. allocated(evc_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_evc_d', 'PANIC: sync of evc from evc_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(evc)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of evc with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    evc_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied evc D->H"
                evc = evc_d
             END IF
             evc_ood = .false.
         ENDIF
         IF (intento_ > 0)    evc_d_ood = .true.
#endif
     END SUBROUTINE using_evc
     !
     SUBROUTINE using_evc_d(intento, debug_info)
         !
         USE wavefunctions, ONLY : evc
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_evc_d ", debug_info, evc_d_ood
         !
         IF (.not. allocated(evc)) THEN
             IF (intento /= 2) print *, "WARNING: sync of evc_d with unallocated array and intento /= 2?"
             IF (allocated(evc_d)) DEALLOCATE(evc_d)
             evc_d_ood = .false.
             RETURN
         END IF
         ! here we know that evc is allocated, check if size is 0 
         IF ( SIZE(evc) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array evc_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (evc_d_ood) THEN
             IF ( allocated(evc_d) .and. (SIZE(evc_d)/=SIZE(evc))) deallocate(evc_d)
             IF (.not. allocated(evc_d)) ALLOCATE(evc_d(DIMS2D(evc)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied evc H->D"
                evc_d = evc
             END IF
             evc_d_ood = .false.
         ENDIF
         IF (intento > 0)    evc_ood = .true.
#else
         CALL errore('using_evc_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_evc_d
     !
     SUBROUTINE using_psic(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wavefunctions, ONLY : psic
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_psic ", debug_info, psic_ood
         !
         IF (psic_ood) THEN
             IF ((.not. allocated(psic_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_psic_d', 'PANIC: sync of psic from psic_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(psic)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of psic with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    psic_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied psic D->H"
                psic = psic_d
             END IF
             psic_ood = .false.
         ENDIF
         IF (intento_ > 0)    psic_d_ood = .true.
#endif
     END SUBROUTINE using_psic
     !
     SUBROUTINE using_psic_d(intento, debug_info)
         !
         USE wavefunctions, ONLY : psic
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_psic_d ", debug_info, psic_d_ood
         !
         IF (.not. allocated(psic)) THEN
             IF (intento /= 2) print *, "WARNING: sync of psic_d with unallocated array and intento /= 2?"
             IF (allocated(psic_d)) DEALLOCATE(psic_d)
             psic_d_ood = .false.
             RETURN
         END IF
         ! here we know that psic is allocated, check if size is 0 
         IF ( SIZE(psic) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array psic_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (psic_d_ood) THEN
             IF ( allocated(psic_d) .and. (SIZE(psic_d)/=SIZE(psic))) deallocate(psic_d)
             IF (.not. allocated(psic_d)) ALLOCATE(psic_d(DIMS1D(psic)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied psic H->D"
                psic_d = psic
             END IF
             psic_d_ood = .false.
         ENDIF
         IF (intento > 0)    psic_ood = .true.
#else
         CALL errore('using_psic_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_psic_d
     !
     SUBROUTINE using_psic_nc(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE wavefunctions, ONLY : psic_nc
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, "using_psic_nc ", debug_info, psic_nc_ood
         !
         IF (psic_nc_ood) THEN
             IF ((.not. allocated(psic_nc_d)) .and. (intento_ < 2)) THEN
                CALL errore('using_psic_nc_d', 'PANIC: sync of psic_nc from psic_nc_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(psic_nc)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of psic_nc with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2
                END IF
                ! IF (intento_ > 0)    psic_nc_d_ood = .true.
             END IF
             IF (intento_ < 2) THEN
                print *, "Really copied psic_nc D->H"
                psic_nc = psic_nc_d
             END IF
             psic_nc_ood = .false.
         ENDIF
         IF (intento_ > 0)    psic_nc_d_ood = .true.
#endif
     END SUBROUTINE using_psic_nc
     !
     SUBROUTINE using_psic_nc_d(intento, debug_info)
         !
         USE wavefunctions, ONLY : psic_nc
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
#if defined(__CUDA)
         !
         IF (PRESENT(debug_info) ) print *, "using_psic_nc_d ", debug_info, psic_nc_d_ood
         !
         IF (.not. allocated(psic_nc)) THEN
             IF (intento /= 2) print *, "WARNING: sync of psic_nc_d with unallocated array and intento /= 2?"
             IF (allocated(psic_nc_d)) DEALLOCATE(psic_nc_d)
             psic_nc_d_ood = .false.
             RETURN
         END IF
         ! here we know that psic_nc is allocated, check if size is 0 
         IF ( SIZE(psic_nc) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array psic_nc_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (psic_nc_d_ood) THEN
             IF ( allocated(psic_nc_d) .and. (SIZE(psic_nc_d)/=SIZE(psic_nc))) deallocate(psic_nc_d)
             IF (.not. allocated(psic_nc_d)) ALLOCATE(psic_nc_d(DIMS2D(psic_nc)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                print *, "Really copied psic_nc H->D"
                psic_nc_d = psic_nc
             END IF
             psic_nc_d_ood = .false.
         ENDIF
         IF (intento > 0)    psic_nc_ood = .true.
#else
         CALL errore('using_psic_nc_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_psic_nc_d
     !
     SUBROUTINE deallocate_wavefunctions_gpu
       IF( ALLOCATED( evc_d ) ) DEALLOCATE( evc_d )
       evc_d_ood = .false.
       IF( ALLOCATED( psic_d ) ) DEALLOCATE( psic_d )
       psic_d_ood = .false.
       IF( ALLOCATED( psic_nc_d ) ) DEALLOCATE( psic_nc_d )
       psic_nc_d_ood = .false.
     END SUBROUTINE deallocate_wavefunctions_gpu
!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_gpum
!=----------------------------------------------------------------------------=!