!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE wavefunctions_module_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE
     !
     COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc_d(:,:)     ! wavefunctions in the PW basis set
                    ! noncolinear case: first index
                    ! is a combined PW + spin index
     !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: &
       psic_d(:), &      ! additional memory for FFT
       psic_nc_d(:,:)    ! as above for the noncolinear case
     !
#if defined(__CUDA)
     attributes (DEVICE) :: psic_d, psic_nc_d, evc_d
     
     LOGICAL :: evc_ood = .false.    ! used to flag out of date variables
     LOGICAL :: evc_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_evc(changing)
         !
         USE wavefunctions_module, ONLY : evc
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (evc_ood) THEN
             IF (.not. allocated(evc_d)) THEN
                CALL errore('using_evc_d', 'PANIC: sync of evc from evc_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(evc)) THEN
                print *, "WARNING: sync of evc with unallocated array. Bye!"
                IF (changing)    evc_d_ood = .true.
                return
             END IF
             print *, "Really copied D->H"
             evc = evc_d
             evc_ood = .false.
         ENDIF
         IF (changing)    evc_d_ood = .true.
#endif
     END SUBROUTINE using_evc
     !
     SUBROUTINE using_evc_d(changing)
         !
         USE wavefunctions_module, ONLY : evc
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(evc)) THEN
             IF (allocated(evc_d)) DEALLOCATE(evc_d)
             evc_d_ood = .false.
             RETURN
         END IF
         IF (evc_d_ood) THEN
             IF ( allocated(evc_d) .and. (SIZE(evc_d)/=SIZE(evc))) deallocate(evc_d)
             IF (.not. allocated(evc_d)) THEN
                 ALLOCATE(evc_d, SOURCE=evc)
             ELSE
                 print *, "Really copied H->D"
                 evc_d = evc
             ENDIF
             evc_d_ood = .false.
         ENDIF
         IF (changing)    evc_ood = .true.
#else
         CALL errore('using_evc_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_evc_d
     !
     SUBROUTINE deallocate_wavefunctions_gpu
       IF( ALLOCATED( psic_nc_d ) ) DEALLOCATE( psic_nc_D )
       IF( ALLOCATED( psic_d ) ) DEALLOCATE( psic_d )
       IF( ALLOCATED( evc_d ) ) DEALLOCATE( evc_d )
     END SUBROUTINE deallocate_wavefunctions_gpu
!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_module_gpum
!=----------------------------------------------------------------------------=!
