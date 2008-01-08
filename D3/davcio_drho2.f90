!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE davcio_drho2 (drho, lrec, iunit, nrec, isw)
  !-----------------------------------------------------------------------
  !
  ! reads/writes variation of the charge with respect to a perturbation
  ! on a file.
  ! isw = +1 : gathers data from the nodes and writes on a single file
  ! isw = -1 : reads data from a single file and distributes them
  !
  USE pwcom
  USE kinds,     ONLY : DP
  USE phcom
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : intra_pool_comm, me_pool, root_pool
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE fft_base,  ONLY : dfftp, cgather_sym
  !
  IMPLICIT NONE
  !
  INTEGER :: iunit, lrec, nrec, isw
  COMPLEX(DP) :: drho (nrxx)
#ifdef __PARA
  !
  ! local variables
  !
  INTEGER :: root, errcode, itmp, proc
  COMPLEX(DP), ALLOCATABLE :: ddrho (:)

  ALLOCATE (ddrho( nrx1 * nrx2 * nrx3 ))    

  IF (isw == 1) THEN
     !
     ! First task of the pool gathers and writes in the file
     !
     CALL cgather_sym (drho, ddrho)
     root = 0
     CALL mp_barrier()
     IF ( me_pool == root_pool ) CALL davcio (ddrho, lrec, iunit, nrec, + 1)
  ELSEIF (isw < 0) THEN
     !
     ! First task of the pool reads ddrho, and broadcasts to all the
     ! processors of the pool
     !
     IF ( me_pool == root_pool ) CALL davcio (ddrho, lrec, iunit, nrec, - 1)
     CALL mp_bcast( ddrho, root_pool, intra_pool_comm )
     !
     ! Distributes ddrho between between the tasks of the pool
     !
     itmp = 1
     DO proc = 1, me_pool
        itmp = itmp + dfftp%nnp * dfftp%npp (proc)
     ENDDO
     drho (:) = (0.d0, 0.d0)
     CALL ZCOPY (dfftp%nnp * dfftp%npp (me_pool+1), ddrho (itmp), 1, drho, 1)
  ENDIF

  DEALLOCATE(ddrho)
#else
  CALL davcio (drho, lrec, iunit, nrec, isw)
#endif
  RETURN
END SUBROUTINE davcio_drho2
