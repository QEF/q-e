!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE gvect_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE, TARGET :: gg_d(:)

     REAL(DP), ALLOCATABLE, TARGET :: g_d(:, :)

     !
#if defined(__CUDA)
     attributes (DEVICE) :: gg_d, g_d

     LOGICAL :: gg_ood = .false.    ! used to flag out of date variables
     LOGICAL :: gg_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: g_ood = .false.    ! used to flag out of date variables
     LOGICAL :: g_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_gg(changing)
         !
         USE gvect, ONLY : gg
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (gg_ood) THEN
             IF (.not. allocated(gg_d)) THEN
                CALL errore('using_gg_d', 'PANIC: sync of gg from gg_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(gg)) THEN
                print *, "WARNING: sync of gg with unallocated array. Bye!"
                IF (changing)    gg_d_ood = .true.
                return
             END IF
             print *, "Really copied gg D->H"
             gg = gg_d
             gg_ood = .false.
         ENDIF
         IF (changing)    gg_d_ood = .true.
#endif
     END SUBROUTINE using_gg
     !
     SUBROUTINE using_gg_d(changing)
         !
         USE gvect, ONLY : gg
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(gg)) THEN
             IF (allocated(gg_d)) DEALLOCATE(gg_d)
             gg_d_ood = .false.
             RETURN
         END IF
         IF (gg_d_ood) THEN
             IF ( allocated(gg_d) .and. (SIZE(gg_d)/=SIZE(gg))) deallocate(gg_d)
             IF (.not. allocated(gg_d)) THEN
                 ALLOCATE(gg_d, SOURCE=gg)
             ELSE
                 print *, "Really copied gg H->D"
                 gg_d = gg
             ENDIF
             gg_d_ood = .false.
         ENDIF
         IF (changing)    gg_ood = .true.
#else
         CALL errore('using_gg_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_gg_d
     !
     SUBROUTINE using_g(changing)
         !
         USE gvect, ONLY : g
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (g_ood) THEN
             IF (.not. allocated(g_d)) THEN
                CALL errore('using_g_d', 'PANIC: sync of g from g_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(g)) THEN
                !print *, "WARNING: sync of g with unallocated array. Bye!"
                IF (changing)    g_d_ood = .true.
                return
             END IF
             print *, "Really copied g D->H"
             g = g_d
             g_ood = .false.
         ENDIF
         IF (changing)    g_d_ood = .true.
#endif
     END SUBROUTINE using_g
     !
     SUBROUTINE using_g_d(changing)
         !
         USE gvect, ONLY : g
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(g)) THEN
             IF (allocated(g_d)) DEALLOCATE(g_d)
             g_d_ood = .false.
             RETURN
         END IF
         IF (g_d_ood) THEN
             IF ( allocated(g_d) .and. (SIZE(g_d)/=SIZE(g))) deallocate(g_d)
             IF (.not. allocated(g_d)) THEN
                 ALLOCATE(g_d, SOURCE=g)
             ELSE
                 print *, "Really copied g H->D"
                 g_d = g
             ENDIF
             g_d_ood = .false.
         ENDIF
         IF (changing)    g_ood = .true.
#else
         CALL errore('using_g_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_g_d
     !     
     SUBROUTINE deallocate_gvect_gpu
       IF( ALLOCATED( gg_d ) ) DEALLOCATE( gg_d )
       IF( ALLOCATED( g_d ) ) DEALLOCATE( g_d )
     END SUBROUTINE deallocate_gvect_gpu
!=----------------------------------------------------------------------------=!
   END MODULE gvect_gpum
!=----------------------------------------------------------------------------=!
