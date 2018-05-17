!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)
!=----------------------------------------------------------------------------=!
   MODULE uspp_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE
     !
     INTEGER, ALLOCATABLE, TARGET :: indv_d(:, :)

     INTEGER, ALLOCATABLE, TARGET :: nhtol_d(:, :)

     INTEGER, ALLOCATABLE, TARGET :: nhtolm_d(:, :)

     INTEGER, ALLOCATABLE, TARGET :: ijtoh_d(:, :, :)

     INTEGER, ALLOCATABLE, TARGET :: indv_ijkb0_d(:)

     COMPLEX(DP), ALLOCATABLE, TARGET :: vkb_d(:, :)

     COMPLEX(DP), ALLOCATABLE, TARGET :: becsum_d(:, :, :)

     REAL(DP), ALLOCATABLE, TARGET :: dvan_d(:, :, :)

     REAL(DP), ALLOCATABLE, TARGET :: deeq_d(:, :, :, :)

     REAL(DP), ALLOCATABLE, TARGET :: qq_nt_d(:, :, :)

     REAL(DP), ALLOCATABLE, TARGET :: qq_at_d(:, :, :)

     REAL(DP), ALLOCATABLE, TARGET :: nhtoj_d(:, :)

     COMPLEX(DP), ALLOCATABLE, TARGET :: qq_so_d(:, :, :, :)

     COMPLEX(DP), ALLOCATABLE, TARGET :: dvan_so_d(:, :, :, :)

     COMPLEX(DP), ALLOCATABLE, TARGET :: deeq_nc_d(:, :, :, :)

     !
#if defined(__CUDA)
     attributes (DEVICE) :: indv_d, nhtol_d, nhtolm_d, ijtoh_d, indv_ijkb0_d, vkb_d, becsum_d, dvan_d, deeq_d, qq_nt_d, qq_at_d, nhtoj_d, qq_so_d, dvan_so_d, deeq_nc_d

     LOGICAL :: indv_ood = .false.    ! used to flag out of date variables
     LOGICAL :: indv_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtol_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtol_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtolm_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtolm_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: ijtoh_ood = .false.    ! used to flag out of date variables
     LOGICAL :: ijtoh_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: indv_ijkb0_ood = .false.    ! used to flag out of date variables
     LOGICAL :: indv_ijkb0_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: vkb_ood = .false.    ! used to flag out of date variables
     LOGICAL :: vkb_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becsum_ood = .false.    ! used to flag out of date variables
     LOGICAL :: becsum_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: dvan_ood = .false.    ! used to flag out of date variables
     LOGICAL :: dvan_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: deeq_ood = .false.    ! used to flag out of date variables
     LOGICAL :: deeq_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_nt_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_nt_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_at_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_at_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtoj_ood = .false.    ! used to flag out of date variables
     LOGICAL :: nhtoj_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_so_ood = .false.    ! used to flag out of date variables
     LOGICAL :: qq_so_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: dvan_so_ood = .false.    ! used to flag out of date variables
     LOGICAL :: dvan_so_d_ood = .false.    ! used to flag out of date variables
     LOGICAL :: deeq_nc_ood = .false.    ! used to flag out of date variables
     LOGICAL :: deeq_nc_d_ood = .false.    ! used to flag out of date variables
     !
#endif
     CONTAINS
     !
     SUBROUTINE using_indv(changing)
         !
         USE uspp, ONLY : indv
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (indv_ood) THEN
             IF (.not. allocated(indv_d)) THEN
                CALL errore('using_indv_d', 'PANIC: sync of indv from indv_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(indv)) THEN
                !print *, "WARNING: sync of indv with unallocated array. Bye!"
                IF (changing)    indv_d_ood = .true.
                return
             END IF
             print *, "Really copied indv D->H"
             indv = indv_d
             indv_ood = .false.
         ENDIF
         IF (changing)    indv_d_ood = .true.
#endif
     END SUBROUTINE using_indv
     !
     SUBROUTINE using_indv_d(changing)
         !
         USE uspp, ONLY : indv
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(indv)) THEN
             IF (allocated(indv_d)) DEALLOCATE(indv_d)
             indv_d_ood = .false.
             RETURN
         END IF
         IF (indv_d_ood) THEN
             IF ( allocated(indv_d) .and. (SIZE(indv_d)/=SIZE(indv))) deallocate(indv_d)
             IF (.not. allocated(indv_d)) THEN
                 ALLOCATE(indv_d, SOURCE=indv)
             ELSE
                 print *, "Really copied indv H->D"
                 indv_d = indv
             ENDIF
             indv_d_ood = .false.
         ENDIF
         IF (changing)    indv_ood = .true.
#else
         CALL errore('using_indv_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_indv_d
     !
     SUBROUTINE using_nhtol(changing)
         !
         USE uspp, ONLY : nhtol
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (nhtol_ood) THEN
             IF (.not. allocated(nhtol_d)) THEN
                CALL errore('using_nhtol_d', 'PANIC: sync of nhtol from nhtol_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(nhtol)) THEN
                !print *, "WARNING: sync of nhtol with unallocated array. Bye!"
                IF (changing)    nhtol_d_ood = .true.
                return
             END IF
             print *, "Really copied nhtol D->H"
             nhtol = nhtol_d
             nhtol_ood = .false.
         ENDIF
         IF (changing)    nhtol_d_ood = .true.
#endif
     END SUBROUTINE using_nhtol
     !
     SUBROUTINE using_nhtol_d(changing)
         !
         USE uspp, ONLY : nhtol
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(nhtol)) THEN
             IF (allocated(nhtol_d)) DEALLOCATE(nhtol_d)
             nhtol_d_ood = .false.
             RETURN
         END IF
         IF (nhtol_d_ood) THEN
             IF ( allocated(nhtol_d) .and. (SIZE(nhtol_d)/=SIZE(nhtol))) deallocate(nhtol_d)
             IF (.not. allocated(nhtol_d)) THEN
                 ALLOCATE(nhtol_d, SOURCE=nhtol)
             ELSE
                 print *, "Really copied nhtol H->D"
                 nhtol_d = nhtol
             ENDIF
             nhtol_d_ood = .false.
         ENDIF
         IF (changing)    nhtol_ood = .true.
#else
         CALL errore('using_nhtol_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_nhtol_d
     !
     SUBROUTINE using_nhtolm(changing)
         !
         USE uspp, ONLY : nhtolm
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (nhtolm_ood) THEN
             IF (.not. allocated(nhtolm_d)) THEN
                CALL errore('using_nhtolm_d', 'PANIC: sync of nhtolm from nhtolm_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(nhtolm)) THEN
                !print *, "WARNING: sync of nhtolm with unallocated array. Bye!"
                IF (changing)    nhtolm_d_ood = .true.
                return
             END IF
             print *, "Really copied nhtolm D->H"
             nhtolm = nhtolm_d
             nhtolm_ood = .false.
         ENDIF
         IF (changing)    nhtolm_d_ood = .true.
#endif
     END SUBROUTINE using_nhtolm
     !
     SUBROUTINE using_nhtolm_d(changing)
         !
         USE uspp, ONLY : nhtolm
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(nhtolm)) THEN
             IF (allocated(nhtolm_d)) DEALLOCATE(nhtolm_d)
             nhtolm_d_ood = .false.
             RETURN
         END IF
         IF (nhtolm_d_ood) THEN
             IF ( allocated(nhtolm_d) .and. (SIZE(nhtolm_d)/=SIZE(nhtolm))) deallocate(nhtolm_d)
             IF (.not. allocated(nhtolm_d)) THEN
                 ALLOCATE(nhtolm_d, SOURCE=nhtolm)
             ELSE
                 print *, "Really copied nhtolm H->D"
                 nhtolm_d = nhtolm
             ENDIF
             nhtolm_d_ood = .false.
         ENDIF
         IF (changing)    nhtolm_ood = .true.
#else
         CALL errore('using_nhtolm_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_nhtolm_d
     !
     SUBROUTINE using_ijtoh(changing)
         !
         USE uspp, ONLY : ijtoh
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (ijtoh_ood) THEN
             IF (.not. allocated(ijtoh_d)) THEN
                CALL errore('using_ijtoh_d', 'PANIC: sync of ijtoh from ijtoh_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(ijtoh)) THEN
                !print *, "WARNING: sync of ijtoh with unallocated array. Bye!"
                IF (changing)    ijtoh_d_ood = .true.
                return
             END IF
             print *, "Really copied ijtoh D->H"
             ijtoh = ijtoh_d
             ijtoh_ood = .false.
         ENDIF
         IF (changing)    ijtoh_d_ood = .true.
#endif
     END SUBROUTINE using_ijtoh
     !
     SUBROUTINE using_ijtoh_d(changing)
         !
         USE uspp, ONLY : ijtoh
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(ijtoh)) THEN
             IF (allocated(ijtoh_d)) DEALLOCATE(ijtoh_d)
             ijtoh_d_ood = .false.
             RETURN
         END IF
         IF (ijtoh_d_ood) THEN
             IF ( allocated(ijtoh_d) .and. (SIZE(ijtoh_d)/=SIZE(ijtoh))) deallocate(ijtoh_d)
             IF (.not. allocated(ijtoh_d)) THEN
                 ALLOCATE(ijtoh_d, SOURCE=ijtoh)
             ELSE
                 print *, "Really copied ijtoh H->D"
                 ijtoh_d = ijtoh
             ENDIF
             ijtoh_d_ood = .false.
         ENDIF
         IF (changing)    ijtoh_ood = .true.
#else
         CALL errore('using_ijtoh_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_ijtoh_d
     !
     SUBROUTINE using_indv_ijkb0(changing)
         !
         USE uspp, ONLY : indv_ijkb0
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (indv_ijkb0_ood) THEN
             IF (.not. allocated(indv_ijkb0_d)) THEN
                CALL errore('using_indv_ijkb0_d', 'PANIC: sync of indv_ijkb0 from indv_ijkb0_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(indv_ijkb0)) THEN
                !print *, "WARNING: sync of indv_ijkb0 with unallocated array. Bye!"
                IF (changing)    indv_ijkb0_d_ood = .true.
                return
             END IF
             print *, "Really copied indv_ijkb0 D->H"
             indv_ijkb0 = indv_ijkb0_d
             indv_ijkb0_ood = .false.
         ENDIF
         IF (changing)    indv_ijkb0_d_ood = .true.
#endif
     END SUBROUTINE using_indv_ijkb0
     !
     SUBROUTINE using_indv_ijkb0_d(changing)
         !
         USE uspp, ONLY : indv_ijkb0
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(indv_ijkb0)) THEN
             IF (allocated(indv_ijkb0_d)) DEALLOCATE(indv_ijkb0_d)
             indv_ijkb0_d_ood = .false.
             RETURN
         END IF
         IF (indv_ijkb0_d_ood) THEN
             IF ( allocated(indv_ijkb0_d) .and. (SIZE(indv_ijkb0_d)/=SIZE(indv_ijkb0))) deallocate(indv_ijkb0_d)
             IF (.not. allocated(indv_ijkb0_d)) THEN
                 ALLOCATE(indv_ijkb0_d, SOURCE=indv_ijkb0)
             ELSE
                 print *, "Really copied indv_ijkb0 H->D"
                 indv_ijkb0_d = indv_ijkb0
             ENDIF
             indv_ijkb0_d_ood = .false.
         ENDIF
         IF (changing)    indv_ijkb0_ood = .true.
#else
         CALL errore('using_indv_ijkb0_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_indv_ijkb0_d
     !
     SUBROUTINE using_vkb(changing)
         !
         USE uspp, ONLY : vkb
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (vkb_ood) THEN
             IF (.not. allocated(vkb_d)) THEN
                CALL errore('using_vkb_d', 'PANIC: sync of vkb from vkb_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(vkb)) THEN
                !print *, "WARNING: sync of vkb with unallocated array. Bye!"
                IF (changing)    vkb_d_ood = .true.
                return
             END IF
             print *, "Really copied vkb D->H"
             vkb = vkb_d
             vkb_ood = .false.
         ENDIF
         IF (changing)    vkb_d_ood = .true.
#endif
     END SUBROUTINE using_vkb
     !
     SUBROUTINE using_vkb_d(changing)
         !
         USE uspp, ONLY : vkb
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(vkb)) THEN
             IF (allocated(vkb_d)) DEALLOCATE(vkb_d)
             vkb_d_ood = .false.
             RETURN
         END IF
         IF (vkb_d_ood) THEN
             IF ( allocated(vkb_d) .and. (SIZE(vkb_d)/=SIZE(vkb))) deallocate(vkb_d)
             IF (.not. allocated(vkb_d)) THEN
                 ALLOCATE(vkb_d, SOURCE=vkb)
             ELSE
                 print *, "Really copied vkb H->D"
                 vkb_d = vkb
             ENDIF
             vkb_d_ood = .false.
         ENDIF
         IF (changing)    vkb_ood = .true.
#else
         CALL errore('using_vkb_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_vkb_d
     !
     SUBROUTINE using_becsum(changing)
         !
         USE uspp, ONLY : becsum
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (becsum_ood) THEN
             IF (.not. allocated(becsum_d)) THEN
                CALL errore('using_becsum_d', 'PANIC: sync of becsum from becsum_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(becsum)) THEN
                !print *, "WARNING: sync of becsum with unallocated array. Bye!"
                IF (changing)    becsum_d_ood = .true.
                return
             END IF
             print *, "Really copied becsum D->H"
             becsum = becsum_d
             becsum_ood = .false.
         ENDIF
         IF (changing)    becsum_d_ood = .true.
#endif
     END SUBROUTINE using_becsum
     !
     SUBROUTINE using_becsum_d(changing)
         !
         USE uspp, ONLY : becsum
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(becsum)) THEN
             IF (allocated(becsum_d)) DEALLOCATE(becsum_d)
             becsum_d_ood = .false.
             RETURN
         END IF
         IF (becsum_d_ood) THEN
             IF ( allocated(becsum_d) .and. (SIZE(becsum_d)/=SIZE(becsum))) deallocate(becsum_d)
             IF (.not. allocated(becsum_d)) THEN
                 ALLOCATE(becsum_d, SOURCE=becsum)
             ELSE
                 print *, "Really copied becsum H->D"
                 becsum_d = becsum
             ENDIF
             becsum_d_ood = .false.
         ENDIF
         IF (changing)    becsum_ood = .true.
#else
         CALL errore('using_becsum_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_becsum_d
     !
     SUBROUTINE using_dvan(changing)
         !
         USE uspp, ONLY : dvan
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (dvan_ood) THEN
             IF (.not. allocated(dvan_d)) THEN
                CALL errore('using_dvan_d', 'PANIC: sync of dvan from dvan_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(dvan)) THEN
                !print *, "WARNING: sync of dvan with unallocated array. Bye!"
                IF (changing)    dvan_d_ood = .true.
                return
             END IF
             print *, "Really copied dvan D->H"
             dvan = dvan_d
             dvan_ood = .false.
         ENDIF
         IF (changing)    dvan_d_ood = .true.
#endif
     END SUBROUTINE using_dvan
     !
     SUBROUTINE using_dvan_d(changing)
         !
         USE uspp, ONLY : dvan
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(dvan)) THEN
             IF (allocated(dvan_d)) DEALLOCATE(dvan_d)
             dvan_d_ood = .false.
             RETURN
         END IF
         IF (dvan_d_ood) THEN
             IF ( allocated(dvan_d) .and. (SIZE(dvan_d)/=SIZE(dvan))) deallocate(dvan_d)
             IF (.not. allocated(dvan_d)) THEN
                 ALLOCATE(dvan_d, SOURCE=dvan)
             ELSE
                 print *, "Really copied dvan H->D"
                 dvan_d = dvan
             ENDIF
             dvan_d_ood = .false.
         ENDIF
         IF (changing)    dvan_ood = .true.
#else
         CALL errore('using_dvan_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_dvan_d
     !
     SUBROUTINE using_deeq(changing)
         !
         USE uspp, ONLY : deeq
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (deeq_ood) THEN
             IF (.not. allocated(deeq_d)) THEN
                CALL errore('using_deeq_d', 'PANIC: sync of deeq from deeq_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(deeq)) THEN
                !print *, "WARNING: sync of deeq with unallocated array. Bye!"
                IF (changing)    deeq_d_ood = .true.
                return
             END IF
             print *, "Really copied deeq D->H"
             deeq = deeq_d
             deeq_ood = .false.
         ENDIF
         IF (changing)    deeq_d_ood = .true.
#endif
     END SUBROUTINE using_deeq
     !
     SUBROUTINE using_deeq_d(changing)
         !
         USE uspp, ONLY : deeq
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(deeq)) THEN
             IF (allocated(deeq_d)) DEALLOCATE(deeq_d)
             deeq_d_ood = .false.
             RETURN
         END IF
         IF (deeq_d_ood) THEN
             IF ( allocated(deeq_d) .and. (SIZE(deeq_d)/=SIZE(deeq))) deallocate(deeq_d)
             IF (.not. allocated(deeq_d)) THEN
                 ALLOCATE(deeq_d, SOURCE=deeq)
             ELSE
                 print *, "Really copied deeq H->D"
                 deeq_d = deeq
             ENDIF
             deeq_d_ood = .false.
         ENDIF
         IF (changing)    deeq_ood = .true.
#else
         CALL errore('using_deeq_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_deeq_d
     !
     SUBROUTINE using_qq_nt(changing)
         !
         USE uspp, ONLY : qq_nt
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (qq_nt_ood) THEN
             IF (.not. allocated(qq_nt_d)) THEN
                CALL errore('using_qq_nt_d', 'PANIC: sync of qq_nt from qq_nt_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(qq_nt)) THEN
                !print *, "WARNING: sync of qq_nt with unallocated array. Bye!"
                IF (changing)    qq_nt_d_ood = .true.
                return
             END IF
             print *, "Really copied qq_nt D->H"
             qq_nt = qq_nt_d
             qq_nt_ood = .false.
         ENDIF
         IF (changing)    qq_nt_d_ood = .true.
#endif
     END SUBROUTINE using_qq_nt
     !
     SUBROUTINE using_qq_nt_d(changing)
         !
         USE uspp, ONLY : qq_nt
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(qq_nt)) THEN
             IF (allocated(qq_nt_d)) DEALLOCATE(qq_nt_d)
             qq_nt_d_ood = .false.
             RETURN
         END IF
         IF (qq_nt_d_ood) THEN
             IF ( allocated(qq_nt_d) .and. (SIZE(qq_nt_d)/=SIZE(qq_nt))) deallocate(qq_nt_d)
             IF (.not. allocated(qq_nt_d)) THEN
                 ALLOCATE(qq_nt_d, SOURCE=qq_nt)
             ELSE
                 print *, "Really copied qq_nt H->D"
                 qq_nt_d = qq_nt
             ENDIF
             qq_nt_d_ood = .false.
         ENDIF
         IF (changing)    qq_nt_ood = .true.
#else
         CALL errore('using_qq_nt_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_qq_nt_d
     !
     SUBROUTINE using_qq_at(changing)
         !
         USE uspp, ONLY : qq_at
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (qq_at_ood) THEN
             IF (.not. allocated(qq_at_d)) THEN
                CALL errore('using_qq_at_d', 'PANIC: sync of qq_at from qq_at_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(qq_at)) THEN
                !print *, "WARNING: sync of qq_at with unallocated array. Bye!"
                IF (changing)    qq_at_d_ood = .true.
                return
             END IF
             print *, "Really copied qq_at D->H"
             qq_at = qq_at_d
             qq_at_ood = .false.
         ENDIF
         IF (changing)    qq_at_d_ood = .true.
#endif
     END SUBROUTINE using_qq_at
     !
     SUBROUTINE using_qq_at_d(changing)
         !
         USE uspp, ONLY : qq_at
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(qq_at)) THEN
             IF (allocated(qq_at_d)) DEALLOCATE(qq_at_d)
             qq_at_d_ood = .false.
             RETURN
         END IF
         IF (qq_at_d_ood) THEN
             IF ( allocated(qq_at_d) .and. (SIZE(qq_at_d)/=SIZE(qq_at))) deallocate(qq_at_d)
             IF (.not. allocated(qq_at_d)) THEN
                 ALLOCATE(qq_at_d, SOURCE=qq_at)
             ELSE
                 print *, "Really copied qq_at H->D"
                 qq_at_d = qq_at
             ENDIF
             qq_at_d_ood = .false.
         ENDIF
         IF (changing)    qq_at_ood = .true.
#else
         CALL errore('using_qq_at_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_qq_at_d
     !
     SUBROUTINE using_nhtoj(changing)
         !
         USE uspp, ONLY : nhtoj
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (nhtoj_ood) THEN
             IF (.not. allocated(nhtoj_d)) THEN
                CALL errore('using_nhtoj_d', 'PANIC: sync of nhtoj from nhtoj_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(nhtoj)) THEN
                !print *, "WARNING: sync of nhtoj with unallocated array. Bye!"
                IF (changing)    nhtoj_d_ood = .true.
                return
             END IF
             print *, "Really copied nhtoj D->H"
             nhtoj = nhtoj_d
             nhtoj_ood = .false.
         ENDIF
         IF (changing)    nhtoj_d_ood = .true.
#endif
     END SUBROUTINE using_nhtoj
     !
     SUBROUTINE using_nhtoj_d(changing)
         !
         USE uspp, ONLY : nhtoj
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(nhtoj)) THEN
             IF (allocated(nhtoj_d)) DEALLOCATE(nhtoj_d)
             nhtoj_d_ood = .false.
             RETURN
         END IF
         IF (nhtoj_d_ood) THEN
             IF ( allocated(nhtoj_d) .and. (SIZE(nhtoj_d)/=SIZE(nhtoj))) deallocate(nhtoj_d)
             IF (.not. allocated(nhtoj_d)) THEN
                 ALLOCATE(nhtoj_d, SOURCE=nhtoj)
             ELSE
                 print *, "Really copied nhtoj H->D"
                 nhtoj_d = nhtoj
             ENDIF
             nhtoj_d_ood = .false.
         ENDIF
         IF (changing)    nhtoj_ood = .true.
#else
         CALL errore('using_nhtoj_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_nhtoj_d
     !
     SUBROUTINE using_qq_so(changing)
         !
         USE uspp, ONLY : qq_so
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (qq_so_ood) THEN
             IF (.not. allocated(qq_so_d)) THEN
                CALL errore('using_qq_so_d', 'PANIC: sync of qq_so from qq_so_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(qq_so)) THEN
                !print *, "WARNING: sync of qq_so with unallocated array. Bye!"
                IF (changing)    qq_so_d_ood = .true.
                return
             END IF
             print *, "Really copied qq_so D->H"
             qq_so = qq_so_d
             qq_so_ood = .false.
         ENDIF
         IF (changing)    qq_so_d_ood = .true.
#endif
     END SUBROUTINE using_qq_so
     !
     SUBROUTINE using_qq_so_d(changing)
         !
         USE uspp, ONLY : qq_so
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(qq_so)) THEN
             IF (allocated(qq_so_d)) DEALLOCATE(qq_so_d)
             qq_so_d_ood = .false.
             RETURN
         END IF
         IF (qq_so_d_ood) THEN
             IF ( allocated(qq_so_d) .and. (SIZE(qq_so_d)/=SIZE(qq_so))) deallocate(qq_so_d)
             IF (.not. allocated(qq_so_d)) THEN
                 ALLOCATE(qq_so_d, SOURCE=qq_so)
             ELSE
                 print *, "Really copied qq_so H->D"
                 qq_so_d = qq_so
             ENDIF
             qq_so_d_ood = .false.
         ENDIF
         IF (changing)    qq_so_ood = .true.
#else
         CALL errore('using_qq_so_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_qq_so_d
     !
     SUBROUTINE using_dvan_so(changing)
         !
         USE uspp, ONLY : dvan_so
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (dvan_so_ood) THEN
             IF (.not. allocated(dvan_so_d)) THEN
                CALL errore('using_dvan_so_d', 'PANIC: sync of dvan_so from dvan_so_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(dvan_so)) THEN
                !print *, "WARNING: sync of dvan_so with unallocated array. Bye!"
                IF (changing)    dvan_so_d_ood = .true.
                return
             END IF
             print *, "Really copied dvan_so D->H"
             dvan_so = dvan_so_d
             dvan_so_ood = .false.
         ENDIF
         IF (changing)    dvan_so_d_ood = .true.
#endif
     END SUBROUTINE using_dvan_so
     !
     SUBROUTINE using_dvan_so_d(changing)
         !
         USE uspp, ONLY : dvan_so
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(dvan_so)) THEN
             IF (allocated(dvan_so_d)) DEALLOCATE(dvan_so_d)
             dvan_so_d_ood = .false.
             RETURN
         END IF
         IF (dvan_so_d_ood) THEN
             IF ( allocated(dvan_so_d) .and. (SIZE(dvan_so_d)/=SIZE(dvan_so))) deallocate(dvan_so_d)
             IF (.not. allocated(dvan_so_d)) THEN
                 ALLOCATE(dvan_so_d, SOURCE=dvan_so)
             ELSE
                 print *, "Really copied dvan_so H->D"
                 dvan_so_d = dvan_so
             ENDIF
             dvan_so_d_ood = .false.
         ENDIF
         IF (changing)    dvan_so_ood = .true.
#else
         CALL errore('using_dvan_so_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_dvan_so_d
     !
     SUBROUTINE using_deeq_nc(changing)
         !
         USE uspp, ONLY : deeq_nc
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (deeq_nc_ood) THEN
             IF (.not. allocated(deeq_nc_d)) THEN
                CALL errore('using_deeq_nc_d', 'PANIC: sync of deeq_nc from deeq_nc_d with unallocated array. Bye!!', 1)
                stop
             END IF
             IF (.not. allocated(deeq_nc)) THEN
                !print *, "WARNING: sync of deeq_nc with unallocated array. Bye!"
                IF (changing)    deeq_nc_d_ood = .true.
                return
             END IF
             print *, "Really copied deeq_nc D->H"
             deeq_nc = deeq_nc_d
             deeq_nc_ood = .false.
         ENDIF
         IF (changing)    deeq_nc_d_ood = .true.
#endif
     END SUBROUTINE using_deeq_nc
     !
     SUBROUTINE using_deeq_nc_d(changing)
         !
         USE uspp, ONLY : deeq_nc
         implicit none
         LOGICAL, INTENT(IN) :: changing
#if defined(__CUDA)
         !
         IF (.not. allocated(deeq_nc)) THEN
             IF (allocated(deeq_nc_d)) DEALLOCATE(deeq_nc_d)
             deeq_nc_d_ood = .false.
             RETURN
         END IF
         IF (deeq_nc_d_ood) THEN
             IF ( allocated(deeq_nc_d) .and. (SIZE(deeq_nc_d)/=SIZE(deeq_nc))) deallocate(deeq_nc_d)
             IF (.not. allocated(deeq_nc_d)) THEN
                 ALLOCATE(deeq_nc_d, SOURCE=deeq_nc)
             ELSE
                 print *, "Really copied deeq_nc H->D"
                 deeq_nc_d = deeq_nc
             ENDIF
             deeq_nc_d_ood = .false.
         ENDIF
         IF (changing)    deeq_nc_ood = .true.
#else
         CALL errore('using_deeq_nc_d', 'Trying to use device data without device compilated code!', 1)
#endif
     END SUBROUTINE using_deeq_nc_d
     !     
     SUBROUTINE deallocate_uspp_gpu
       IF( ALLOCATED( indv_d ) ) DEALLOCATE( indv_d )
       IF( ALLOCATED( nhtol_d ) ) DEALLOCATE( nhtol_d )
       IF( ALLOCATED( nhtolm_d ) ) DEALLOCATE( nhtolm_d )
       IF( ALLOCATED( ijtoh_d ) ) DEALLOCATE( ijtoh_d )
       IF( ALLOCATED( indv_ijkb0_d ) ) DEALLOCATE( indv_ijkb0_d )
       IF( ALLOCATED( vkb_d ) ) DEALLOCATE( vkb_d )
       IF( ALLOCATED( becsum_d ) ) DEALLOCATE( becsum_d )
       IF( ALLOCATED( dvan_d ) ) DEALLOCATE( dvan_d )
       IF( ALLOCATED( deeq_d ) ) DEALLOCATE( deeq_d )
       IF( ALLOCATED( qq_nt_d ) ) DEALLOCATE( qq_nt_d )
       IF( ALLOCATED( qq_at_d ) ) DEALLOCATE( qq_at_d )
       IF( ALLOCATED( nhtoj_d ) ) DEALLOCATE( nhtoj_d )
       IF( ALLOCATED( qq_so_d ) ) DEALLOCATE( qq_so_d )
       IF( ALLOCATED( dvan_so_d ) ) DEALLOCATE( dvan_so_d )
       IF( ALLOCATED( deeq_nc_d ) ) DEALLOCATE( deeq_nc_d )
     END SUBROUTINE deallocate_uspp_gpu
!=----------------------------------------------------------------------------=!
   END MODULE uspp_gpum
!=----------------------------------------------------------------------------=!