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
   MODULE becmod_gpum
   !! Routines for the device-host cuda memory management of \(\texttt{becmod}\)-
   !! related quantities.
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
     TYPE bec_type_d
#if defined(__CUDA)
        REAL(DP), ALLOCATABLE, DEVICE :: r_d(:, :)
#else
        REAL(DP), ALLOCATABLE :: r_d(:, :)
#endif
#if defined(__CUDA)
        COMPLEX(DP), ALLOCATABLE, DEVICE :: k_d(:, :)
#else
        COMPLEX(DP), ALLOCATABLE :: k_d(:, :)
#endif
#if defined(__CUDA)
        COMPLEX(DP), ALLOCATABLE, DEVICE :: nc_d(:, :, :)
#else
        COMPLEX(DP), ALLOCATABLE :: nc_d(:, :, :)
#endif
        INTEGER :: comm
        INTEGER :: nbnd
        INTEGER :: nproc
        INTEGER :: mype
        INTEGER :: nbnd_loc
        INTEGER :: ibnd_begin
     END TYPE bec_type_d
     !
     TYPE (bec_type_d), TARGET :: becp_d  ! <beta|psi>
     !

     LOGICAL :: becp_r_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becp_d_r_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becp_k_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becp_d_k_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becp_nc_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becp_d_nc_d_ood = .false.    ! used to flag out of date variables
     !
     CONTAINS
     !
     SUBROUTINE using_becp_r(intento, debug_info)
         !
         ! intento is used to specify what the variable will be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (becp_r_ood) THEN
             IF (.not. allocated(becp_d%r_d)) THEN
                CALL errore('using_r_d', 'PANIC: sync of becp%r from becp_d%r_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(becp%r)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of becp%r with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2 
                END IF
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%r D->H"
                becp%r = becp_d%r_d
             END IF
             ! ALWAYS update auxiliary variables
             !IF ( becp%comm /= becp_d%comm ) &
             !     print *, "WARNING: auxiliary variable becp%comm changed"
             becp%comm = becp_d%comm
             !IF ( becp%nbnd /= becp_d%nbnd ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd changed"
             becp%nbnd = becp_d%nbnd
             !IF ( becp%nproc /= becp_d%nproc ) &
             !     print *, "WARNING: auxiliary variable becp%nproc changed"
             becp%nproc = becp_d%nproc
             !IF ( becp%mype /= becp_d%mype ) &
             !     print *, "WARNING: auxiliary variable becp%mype changed"
             becp%mype = becp_d%mype
             !IF ( becp%nbnd_loc /= becp_d%nbnd_loc ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd_loc changed"
             becp%nbnd_loc = becp_d%nbnd_loc
             !IF ( becp%ibnd_begin /= becp_d%ibnd_begin ) &
             !     print *, "WARNING: auxiliary variable becp%ibnd_begin changed"
             becp%ibnd_begin = becp_d%ibnd_begin
             !
             becp_r_ood = .false.
         ENDIF
         IF (intento_ > 0)    becp_d_r_d_ood = .true.
#endif
     END SUBROUTINE using_becp_r
     !
     SUBROUTINE using_becp_r_d(intento, debug_info)
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (.not. allocated(becp%r)) THEN
             IF (intento /= 2) print *, "WARNING: sync of becp%r_d with unallocated array and intento /= 2?"
             IF (allocated(becp_d%r_d)) DEALLOCATE(becp_d%r_d)
             becp_d_r_d_ood = .false.
             RETURN
         END IF
         ! here we know that r is allocated, check if size is 0 
         IF ( SIZE(becp%r) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array becp_d%r_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (becp_d_r_d_ood) THEN
             IF ( allocated(becp_d%r_d) .and. (SIZE(becp_d%r_d)/=SIZE(becp%r))) deallocate(becp_d%r_d)
             IF (.not. allocated(becp_d%r_d)) ALLOCATE(becp_d%r_d(DIMS2D(becp%r)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%r H->D"
                becp_d%r_d = becp%r
             END IF
             ! ALWAYS update auxiliary variables
             becp_d%comm = becp%comm
             becp_d%nbnd = becp%nbnd
             becp_d%nproc = becp%nproc
             becp_d%mype = becp%mype
             becp_d%nbnd_loc = becp%nbnd_loc
             becp_d%ibnd_begin = becp%ibnd_begin
             !
             becp_d_r_d_ood = .false.
         ENDIF
         IF (intento > 0)    becp_r_ood = .true.
#else
         CALL errore('using_becp_d%r_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_becp_r_d
     !
     SUBROUTINE using_becp_k(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (becp_k_ood) THEN
             IF (.not. allocated(becp_d%k_d)) THEN
                CALL errore('using_k_d', 'PANIC: sync of becp%k from becp_d%k_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(becp%k)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of becp%k with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2 
                END IF
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%k D->H"
                becp%k = becp_d%k_d
             END IF
             ! ALWAYS update auxiliary variables
             !IF ( becp%comm /= becp_d%comm ) &
             !     print *, "WARNING: auxiliary variable becp%comm changed"
             becp%comm = becp_d%comm
             !IF ( becp%nbnd /= becp_d%nbnd ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd changed"
             becp%nbnd = becp_d%nbnd
             !IF ( becp%nproc /= becp_d%nproc ) &
             !     print *, "WARNING: auxiliary variable becp%nproc changed"
             becp%nproc = becp_d%nproc
             !IF ( becp%mype /= becp_d%mype ) &
             !     print *, "WARNING: auxiliary variable becp%mype changed"
             becp%mype = becp_d%mype
             !IF ( becp%nbnd_loc /= becp_d%nbnd_loc ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd_loc changed"
             becp%nbnd_loc = becp_d%nbnd_loc
             !IF ( becp%ibnd_begin /= becp_d%ibnd_begin ) &
             !     print *, "WARNING: auxiliary variable becp%ibnd_begin changed"
             becp%ibnd_begin = becp_d%ibnd_begin
             !
             becp_k_ood = .false.
         ENDIF
         IF (intento_ > 0)    becp_d_k_d_ood = .true.
#endif
     END SUBROUTINE using_becp_k
     !
     SUBROUTINE using_becp_k_d(intento, debug_info)
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (.not. allocated(becp%k)) THEN
             IF (intento /= 2) print *, "WARNING: sync of becp%k_d with unallocated array and intento /= 2?"
             IF (allocated(becp_d%k_d)) DEALLOCATE(becp_d%k_d)
             becp_d_k_d_ood = .false.
             RETURN
         END IF
         ! here we know that k is allocated, check if size is 0 
         IF ( SIZE(becp%k) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array becp_d%k_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (becp_d_k_d_ood) THEN
             IF ( allocated(becp_d%k_d) .and. (SIZE(becp_d%k_d)/=SIZE(becp%k))) deallocate(becp_d%k_d)
             IF (.not. allocated(becp_d%k_d)) ALLOCATE(becp_d%k_d(DIMS2D(becp%k)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%k H->D"
                becp_d%k_d = becp%k
             END IF
             ! ALWAYS update auxiliary variables
             becp_d%comm = becp%comm
             becp_d%nbnd = becp%nbnd
             becp_d%nproc = becp%nproc
             becp_d%mype = becp%mype
             becp_d%nbnd_loc = becp%nbnd_loc
             becp_d%ibnd_begin = becp%ibnd_begin
             !
             becp_d_k_d_ood = .false.
         ENDIF
         IF (intento > 0)    becp_k_ood = .true.
#else
         CALL errore('using_becp_d%k_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_becp_k_d
     !
     SUBROUTINE using_becp_nc(intento, debug_info)
         !
         ! intento is used to specify what the variable will  be used for :
         !  0 -> in , the variable needs to be synchronized but won't be changed
         !  1 -> inout , the variable needs to be synchronized AND will be changed
         !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         INTEGER :: intento_
         intento_ = intento
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (becp_nc_ood) THEN
             IF (.not. allocated(becp_d%nc_d)) THEN
                CALL errore('using_nc_d', 'PANIC: sync of becp%nc from becp_d%nc_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(becp%nc)) THEN
                IF (intento_ /= 2) THEN
                   print *, "WARNING: sync of becp%nc with unallocated array and intento /= 2? Changed to 2!"
                   intento_ = 2 
                END IF
             END IF
             IF (intento_ < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%nc D->H"
                becp%nc = becp_d%nc_d
             END IF
             ! ALWAYS update auxiliary variables
             !IF ( becp%comm /= becp_d%comm ) &
             !     print *, "WARNING: auxiliary variable becp%comm changed"
             becp%comm = becp_d%comm
             !IF ( becp%nbnd /= becp_d%nbnd ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd changed"
             becp%nbnd = becp_d%nbnd
             !IF ( becp%nproc /= becp_d%nproc ) &
             !     print *, "WARNING: auxiliary variable becp%nproc changed"
             becp%nproc = becp_d%nproc
             !IF ( becp%mype /= becp_d%mype ) &
             !     print *, "WARNING: auxiliary variable becp%mype changed"
             becp%mype = becp_d%mype
             !IF ( becp%nbnd_loc /= becp_d%nbnd_loc ) &
             !     print *, "WARNING: auxiliary variable becp%nbnd_loc changed"
             becp%nbnd_loc = becp_d%nbnd_loc
             !IF ( becp%ibnd_begin /= becp_d%ibnd_begin ) &
             !     print *, "WARNING: auxiliary variable becp%ibnd_begin changed"
             becp%ibnd_begin = becp_d%ibnd_begin
             !
             becp_nc_ood = .false.
         ENDIF
         IF (intento_ > 0)    becp_d_nc_d_ood = .true.
#endif
     END SUBROUTINE using_becp_nc
     !
     SUBROUTINE using_becp_nc_d(intento, debug_info)
         !
         USE becmod, ONLY : becp
         implicit none
         INTEGER, INTENT(IN) :: intento
         CHARACTER(len=*), INTENT(IN), OPTIONAL :: debug_info
         !
#if defined(__CUDA) || defined(__CUDA_GNU)
         !
         IF (PRESENT(debug_info) ) print *, debug_info
         !
         IF (.not. allocated(becp%nc)) THEN
             IF (intento /= 2) print *, "WARNING: sync of becp%nc_d with unallocated array and intento /= 2?"
             IF (allocated(becp_d%nc_d)) DEALLOCATE(becp_d%nc_d)
             becp_d_nc_d_ood = .false.
             RETURN
         END IF
         ! here we know that nc is allocated, check if size is 0 
         IF ( SIZE(becp%nc) == 0 ) THEN
             print *, "Refusing to allocate 0 dimensional array becp_d%nc_d. If used, code will crash."
             RETURN
         END IF
         !
         IF (becp_d_nc_d_ood) THEN
             IF ( allocated(becp_d%nc_d) .and. (SIZE(becp_d%nc_d)/=SIZE(becp%nc))) deallocate(becp_d%nc_d)
             IF (.not. allocated(becp_d%nc_d)) ALLOCATE(becp_d%nc_d(DIMS3D(becp%nc)))  ! MOLD does not work on all compilers
             IF (intento < 2) THEN
                IF ( iverbosity > 0 ) print *, "Really copied becp%nc H->D"
                becp_d%nc_d = becp%nc
             END IF
             ! ALWAYS update auxiliary variables
             becp_d%comm = becp%comm
             becp_d%nbnd = becp%nbnd
             becp_d%nproc = becp%nproc
             becp_d%mype = becp%mype
             becp_d%nbnd_loc = becp%nbnd_loc
             becp_d%ibnd_begin = becp%ibnd_begin
             !
             becp_d_nc_d_ood = .false.
         ENDIF
         IF (intento > 0)    becp_nc_ood = .true.
#else
         CALL errore('using_becp_d%nc_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_becp_nc_d
     !     
     SUBROUTINE deallocate_becmod_gpu
       IF( ALLOCATED( becp_d%r_d ) ) DEALLOCATE( becp_d%r_d )
       IF( ALLOCATED( becp_d%k_d ) ) DEALLOCATE( becp_d%k_d )
       IF( ALLOCATED( becp_d%nc_d ) ) DEALLOCATE( becp_d%nc_d )
     END SUBROUTINE deallocate_becmod_gpu
!=----------------------------------------------------------------------------=!
   END MODULE becmod_gpum
!=----------------------------------------------------------------------------=!
