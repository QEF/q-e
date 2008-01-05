!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE fft_base,  ONLY : dfftp
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global, ONLY : inter_pool_comm, me_pool
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE gvect,     ONLY : nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,  ONLY : nspin
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER          :: iunit, lrec, nrec, isw
  COMPLEX(DP) :: drho(nrxx,nspin)
  !
#ifdef __PARA
  !
  ! ... local variables
  !
  INTEGER                        :: itmp, proc, is, dim
  COMPLEX(DP), ALLOCATABLE  :: ddrho(:,:)
  !
  !
  ALLOCATE( ddrho( nrx1 * nrx2 * nrx3 , nspin) )
  !
  IF ( isw == 1 ) THEN
     !
     ! ... First task is the only task allowed to write the file
     !
     DO is = 1, nspin
        !   
        CALL cgather_sym( drho(1,is), ddrho(1,is) )
        !
     END DO
     !
     call mp_barrier()
     !
     IF ( ionode ) CALL davcio( ddrho, lrec, iunit, nrec, + 1 )
     !
  ELSE IF ( isw < 0 ) THEN
     !
     ! ... First task reads and broadcasts ddrho to all pools
     !
     IF ( ionode ) CALL davcio( ddrho, lrec, iunit, nrec, - 1 )
     !
     CALL mp_bcast( ddrho, ionode_id, inter_pool_comm )
     !
     ! ... distributes ddrho between between the tasks of the pool
     !
     DO is = 1, nspin
        !   
        CALL cscatter_sym ( ddrho(1,is), drho(1,is) )
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
! !
#endif
  !
  RETURN
  !
END SUBROUTINE davcio_drho
