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
  ! ... isw = +1 : gathers data from the nodes and writes on a single file
  ! ... isw = -1 : reads data from a single file and distributes them
  !
  use pwcom
  USE kinds,     ONLY : DP
  use phcom
  USE pfft,      ONLY : npp, ncplane
  USE mp_global, ONLY : intra_pool_comm, me_pool, root_pool, my_image_id
  USE mp,        ONLY : mp_bcast, mp_barrier
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
     ! ... First task of each pool is the only task allowed to write
     ! ... the file
     !
     DO is = 1, nspin
        !   
        CALL cgather_sym( drho(1,is), ddrho(1,is) )
        !
     END DO
     !
     call mp_barrier()
     !
     IF ( me_pool == root_pool ) CALL davcio( ddrho, lrec, iunit, nrec, + 1 )
     !
  ELSE IF ( isw < 0 ) THEN
     !
     ! ... First task of the pool reads ddrho, and broadcasts to all the
     ! ... processors of the pool
     !
     IF ( me_pool == root_pool ) CALL davcio( ddrho, lrec, iunit, nrec, - 1 )
     !
     CALL mp_bcast( ddrho, root_pool, intra_pool_comm )
     !
     ! ... Distributes ddrho between between the tasks of the pool
     !
     itmp = 1
     !
     DO proc = 1, me_pool
        !
        itmp = itmp + ncplane * npp(proc)
        !
     END DO
     !
     dim = ncplane * npp(me_pool+1)
     !
     DO is = 1, nspin
        ! 
        drho(:,is) = ( 0.D0, 0.D0 )
        !
        drho(1:dim,is) = ddrho(itmp:itmp+dim,is)
        !CALL ZCOPY( ncplane*npp(me_pool+1), ddrho(itmp,is), 1, drho(1,is), 1 )
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
