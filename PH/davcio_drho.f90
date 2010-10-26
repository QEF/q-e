!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE davcio_drho( drho, lrec, iunit, nrec, isw )
  !----------------------------------------------------------------------------
  !
  ! ... reads/writes variation of the charge with respect to a perturbation
  ! ... on a file.
  ! ... isw = +1 : gathers data from the processors, writes to a single file
  ! ... isw = -1 : reads data from a single file and distributes them
  !
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp, cgather_sym, cscatter_sym
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global, ONLY : inter_pool_comm, me_pool, intra_image_comm
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE noncollin_module,  ONLY : nspin_mag
  USE paw_variables, ONLY : okpaw
  USE phus,     ONLY : int3_paw
  USE units_ph, ONLY : iuint3paw, lint3paw
  !
  IMPLICIT NONE
  !
  INTEGER          :: iunit, lrec, nrec, isw
  COMPLEX(DP) :: drho(dfftp%nnr,nspin_mag)
  !
#ifdef __PARA
  !
  ! ... local variables
  !
  INTEGER                        :: is
  LOGICAL   :: exst
  COMPLEX(DP), ALLOCATABLE  :: ddrho(:,:)
  !
  !
  IF ( ionode ) INQUIRE (UNIT = iunit, OPENED = exst)
  CALL mp_bcast(exst,ionode_id, intra_image_comm)
  IF (.NOT.exst) RETURN

  ALLOCATE( ddrho( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x , nspin_mag) )
  !
  IF ( isw == 1 ) THEN
     !
     ! ... First task is the only task allowed to write the file
     !
     DO is = 1, nspin_mag
        !
        CALL cgather_sym( drho(:,is), ddrho(:,is) )
        !
     END DO
     !
     call mp_barrier(intra_image_comm)
     !
     IF ( ionode ) THEN
        CALL davcio( ddrho, lrec, iunit, nrec, + 1 )
        IF (okpaw) CALL davcio( int3_paw, lint3paw, iuint3paw, nrec, + 1 )
     END IF
     !
  ELSE IF ( isw < 0 ) THEN
     !
     ! ... First task reads and broadcasts ddrho to all pools
     !
     IF ( ionode ) THEN
        CALL davcio( ddrho, lrec, iunit, nrec, - 1 )
        IF (okpaw) CALL davcio( int3_paw, lint3paw, iuint3paw, nrec, - 1 )
     ENDIF
     !
     CALL mp_bcast( ddrho, ionode_id, inter_pool_comm )
     CALL mp_bcast( int3_paw, ionode_id, inter_pool_comm )
     !
     ! ... distributes ddrho between between the tasks of the pool
     !
     DO is = 1, nspin_mag
        !
        CALL cscatter_sym ( ddrho(:,is), drho(:,is) )
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( ddrho )
  !
#else
  !
  CALL davcio( drho, lrec, iunit, nrec, isw )
  IF (okpaw) CALL davcio( int3_paw, lint3paw, iuint3paw, nrec, isw )
! !
#endif
  !
  RETURN
  !
END SUBROUTINE davcio_drho
