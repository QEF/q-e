!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE g_psi_mod_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE, TARGET :: h_diag_d(:, :)

     REAL(DP), ALLOCATABLE, TARGET :: s_diag_d(:, :)

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
     SUBROUTINE using_h_diag(intento)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE g_psi_mod, ONLY : h_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         !
         IF (h_diag_ood) THEN
             IF (.not. allocated(h_diag_d)) THEN
                CALL errore('using_h_diag_d', 'PANIC: sync of h_diag from h_diag_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(h_diag)) THEN
                IF (intento /= 2) print *, "WARNING: sync of h_diag with unallocated array and intento /= 2?"
                IF (intento > 0)    h_diag_d_ood = .true.
                return
             END IF
             IF (intento < 2) THEN
                print *, "Really copied h_diag D->H"
                h_diag = h_diag_d
                h_diag_ood = .false.
             END IF
         ENDIF
         IF (intento > 0)    h_diag_d_ood = .true.
#endif
     END SUBROUTINE using_h_diag
     !
     SUBROUTINE using_h_diag_d(intento)
         !
         USE g_psi_mod, ONLY : h_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         !
         IF (.not. allocated(h_diag)) THEN
             IF (intento /= 2) print *, "WARNING: sync of h_diag_d with unallocated array and intento /= 2?"
             IF (allocated(h_diag_d)) DEALLOCATE(h_diag_d)
             h_diag_d_ood = .false.
             RETURN
         END IF
         ! here we know that h_diag is allocated, check if size if 0 
         IF ( SIZE(h_diag) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array h_diag_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (h_diag_d_ood) THEN
             IF ( allocated(h_diag_d) .and. (SIZE(h_diag_d)/=SIZE(h_diag))) deallocate(h_diag_d)
             IF (.not. allocated(h_diag_d)) ALLOCATE(h_diag_d, MOLD=h_diag)  ! this copy may be avoided
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
     SUBROUTINE using_s_diag(intento)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE g_psi_mod, ONLY : s_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         !
         IF (s_diag_ood) THEN
             IF (.not. allocated(s_diag_d)) THEN
                CALL errore('using_s_diag_d', 'PANIC: sync of s_diag from s_diag_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(s_diag)) THEN
                IF (intento /= 2) print *, "WARNING: sync of s_diag with unallocated array and intento /= 2?"
                IF (intento > 0)    s_diag_d_ood = .true.
                return
             END IF
             IF (intento < 2) THEN
                print *, "Really copied s_diag D->H"
                s_diag = s_diag_d
                s_diag_ood = .false.
             END IF
         ENDIF
         IF (intento > 0)    s_diag_d_ood = .true.
#endif
     END SUBROUTINE using_s_diag
     !
     SUBROUTINE using_s_diag_d(intento)
         !
         USE g_psi_mod, ONLY : s_diag
         implicit none
         INTEGER, INTENT(IN) :: intento
#if defined(__CUDA)
         !
         IF (.not. allocated(s_diag)) THEN
             IF (intento /= 2) print *, "WARNING: sync of s_diag_d with unallocated array and intento /= 2?"
             IF (allocated(s_diag_d)) DEALLOCATE(s_diag_d)
             s_diag_d_ood = .false.
             RETURN
         END IF
         ! here we know that s_diag is allocated, check if size if 0 
         IF ( SIZE(s_diag) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array s_diag_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (s_diag_d_ood) THEN
             IF ( allocated(s_diag_d) .and. (SIZE(s_diag_d)/=SIZE(s_diag))) deallocate(s_diag_d)
             IF (.not. allocated(s_diag_d)) ALLOCATE(s_diag_d, MOLD=s_diag)  ! this copy may be avoided
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
       IF( ALLOCATED( s_diag_d ) ) DEALLOCATE( s_diag_d )
     END SUBROUTINE deallocate_g_psi_mod_gpu
!=----------------------------------------------------------------------------=!
   END MODULE g_psi_mod_gpum
!=----------------------------------------------------------------------------=!