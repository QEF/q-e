!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE scf_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE
     !
     REAL(DP), ALLOCATABLE, TARGET :: vrs_d(:, :)

     !
#if defined(__CUDA)
     attributes (DEVICE) :: vrs_d

     LOGICAL :: vrs_ood = .false.    ! used to flag out of date variables
     LOGICAL :: vrs_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_vrs(changing)
         !
         USE scf, ONLY : vrs
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (vrs_ood) THEN
             IF (.not. allocated(vrs_d)) THEN
                CALL errore('using_vrs_d', 'PANIC: sync of vrs from vrs_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(vrs)) THEN
                !print *, "WARNING: sync of vrs with unallocated array. Bye!"
                IF (changing)    vrs_d_ood = .true.
                return
             END IF
             print *, "Really copied vrs D->H"
             vrs = vrs_d
             vrs_ood = .false.
         ENDIF
         IF (changing)    vrs_d_ood = .true.
#endif
     END SUBROUTINE using_vrs
     !
     SUBROUTINE using_vrs_d(changing)
         !
         USE scf, ONLY : vrs
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(vrs)) THEN
             IF (allocated(vrs_d)) DEALLOCATE(vrs_d)
             vrs_d_ood = .false.
             RETURN
         END IF
         IF (vrs_d_ood) THEN
             IF ( allocated(vrs_d) .and. (SIZE(vrs_d)/=SIZE(vrs))) deallocate(vrs_d)
             IF (.not. allocated(vrs_d)) THEN
                 ALLOCATE(vrs_d, SOURCE=vrs)
             ELSE
                 print *, "Really copied vrs H->D"
                 vrs_d = vrs
             ENDIF
             vrs_d_ood = .false.
         ENDIF
         IF (changing)    vrs_ood = .true.
#else
         CALL errore('using_vrs_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_vrs_d
     !     
     SUBROUTINE deallocate_scf_gpu
       IF( ALLOCATED( vrs_d ) ) DEALLOCATE( vrs_d )
     END SUBROUTINE deallocate_scf_gpu
!=----------------------------------------------------------------------------=!
   END MODULE scf_gpum
!=----------------------------------------------------------------------------=!
