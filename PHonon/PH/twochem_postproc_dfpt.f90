!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE twochem_postproc_dfpt(npe, nsolv, imode0, lmetq0, convt, dos_ef, ldos, ldoss, &
                                 drhop, dbecsum, becsum1)
   USE kinds,                ONLY : DP
   USE mp,                   ONLY : mp_sum
   USE mp_pools,             ONLY : inter_pool_comm
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE fft_interfaces,       ONLY : fft_interpolate
   USE fft_base,             ONLY : dfftp, dffts
   USE noncollin_module,     ONLY : noncolin, nspin_mag, domag
   USE gvecs,                ONLY : doublegrid
   USE ions_base,            ONLY : nat
   USE uspp,                 ONLY : okvan
   USE uspp_param,           ONLY : nhm
   USE paw_variables,        ONLY : okpaw
   !
   USE two_chem,             ONLY : twochem
   USE lr_two_chem
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(in) :: npe
   INTEGER, INTENT(in) :: nsolv
   INTEGER, INTENT(in) :: imode0
   LOGICAL, INTENT(in) :: lmetq0
   LOGICAL, INTENT(in) :: convt
   REAL(DP), INTENT(in) :: dos_ef
   COMPLEX(DP), INTENT(in) :: ldos(dfftp%nnr, nspin_mag)
   COMPLEX(DP), INTENT(in) :: ldoss(dffts%nnr, nspin_mag)
   COMPLEX(DP), INTENT(inout) :: drhop(dfftp%nnr, nspin_mag, npe)
   COMPLEX(DP), INTENT(inout) :: dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag , npe)
   REAL(DP), INTENT(in), OPTIONAL :: becsum1((nhm * (nhm + 1))/2 , nat , nspin_mag)
   !
   INTEGER :: is, ipert
   !
   COMPLEX(DP), ALLOCATABLE :: dbecsum_aux(:,:,:,:)
   !
   IF (.NOT. twochem) CALL errore('twochem_postproc_dfpt', &
      'twochem_postproc_dfpt must be called only if twochem is true', 1)
   !
   IF (noncolin .AND. domag .AND. okvan) THEN
      ALLOCATE(dbecsum_aux((nhm * (nhm + 1))/2, nat, nspin_mag, npe))
   ENDIF
   !
   IF (nsolv==2) THEN
      drhos_cond = drhos_cond / 2.0_DP
      dbecsum_cond = dbecsum_cond / 2.0_DP
      dbecsum_cond_nc = dbecsum_cond_nc / 2.0_DP
   ENDIF
   !
   IF (noncolin) THEN
      CALL mp_sum (dbecsum_cond_nc, intra_bgrp_comm)
   ELSE
      CALL mp_sum (dbecsum_cond, intra_bgrp_comm)
   ENDIF
   !
   if (doublegrid) then
      do is = 1, nspin_mag
         do ipert = 1, npe
            call fft_interpolate (dffts, drhos_cond(:,is,ipert), dfftp, drhop_cond(:,is,ipert))
         enddo
      enddo
   else
      call zcopy (npe*nspin_mag*dfftp%nnr, drhos_cond, 1, drhop_cond, 1)
   endif
   !
   !rotate also dbecsum_cond in the twochem case
   IF (noncolin.and.okvan) THEN
      CALL set_dbecsum_nc(dbecsum_cond_nc, dbecsum_cond, npe)
      IF (nsolv==2) THEN
         dbecsum_aux=(0.0_DP,0.0_DP)
         CALL set_dbecsum_nc(dbecsum_cond_nc(1,1,1,1,1,2), dbecsum_aux, npe)
         dbecsum_cond(:,:,1,:)=dbecsum_cond(:,:,1,:)+dbecsum_aux(:,:,1,:)
         dbecsum_cond(:,:,2:4,:)=dbecsum_cond(:,:,2:4,:)-dbecsum_aux(:,:,2:4,:)
      ENDIF
   ENDIF
   !
   call addusddens_cond(drhop_cond, dbecsum_cond, imode0, npe, 0)
   !
   CALL mp_sum ( drhos_cond, inter_pool_comm )
   CALL mp_sum ( drhop_cond, inter_pool_comm )
   IF (okpaw) call mp_sum ( dbecsum_cond, inter_pool_comm )
   !
   ! q=0 in metallic case deserve special care (e_Fermi can shift)
   !
   IF (lmetq0) THEN
      IF (okpaw) THEN
         CALL ef_shift_twochem(npe, dos_ef, dos_ef_cond, ldos, ldos_cond, drhop,&
                           drhop_cond,dbecsum,dbecsum_cond, becsum1,becsum1_cond)
      ELSE
         CALL ef_shift_twochem(npe, dos_ef,dos_ef_cond, ldos,ldos_cond, drhop,&
                                         drhop_cond)
      ENDIF
   ENDIF
   !
   IF (lmetq0 .AND. convt) CALL ef_shift_wfc_twochem(npe, ldoss, ldoss_cond, drhop)
   !
END SUBROUTINE twochem_postproc_dfpt
