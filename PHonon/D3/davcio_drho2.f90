!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  USE io_global, ONLY : ionode_id, ionode
  USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, me_pool, root_pool
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_world,  ONLY : world_comm
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : gather_grid
  !
  IMPLICIT NONE
  !
  INTEGER :: iunit, lrec, nrec, isw
  COMPLEX(DP) :: drho (dfftp%nnr)
#ifdef __MPI
  !
  ! local variables
  !
  INTEGER :: root, errcode, itmp, proc
  COMPLEX(DP), ALLOCATABLE :: ddrho (:)

  ALLOCATE (ddrho( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

  IF (isw == 1) THEN
     !
     ! First task of the pool gathers and writes in the file
     !
     CALL gather_grid (dfftp, drho, ddrho)
     root = 0
     CALL mp_barrier( world_comm )
     IF ( ionode ) CALL davcio (ddrho, lrec, iunit, nrec, + 1)
  ELSEIF (isw < 0) THEN
     !
     ! First task of the pool reads ddrho, and broadcasts to all the
     ! processors of the pool
     !
     IF ( ionode ) CALL davcio (ddrho, lrec, iunit, nrec, - 1)
     CALL mp_bcast( ddrho, ionode_id, inter_pool_comm )
     CALL mp_bcast( ddrho, root_pool, intra_pool_comm )
     !
     ! Distributes ddrho between between the tasks of the pool
     !
     itmp = 1
     DO proc = 1, me_pool
        itmp = itmp + dfftp%nnp * dfftp%npp (proc)
     ENDDO
     drho (:) = (0.d0, 0.d0)
     CALL zcopy (dfftp%nnp * dfftp%npp (me_pool+1), ddrho (itmp), 1, drho, 1)
  ENDIF

  DEALLOCATE(ddrho)
#else
  CALL davcio (drho, lrec, iunit, nrec, isw)
#endif
  RETURN
END SUBROUTINE davcio_drho2
