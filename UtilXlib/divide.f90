!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE divide (comm, ntodiv, startn, lastn)
  !-----------------------------------------------------------------------
  !
  ! Given "ntodiv" objects, distribute index across a group of processors
  ! belonging to communicator "comm"
  ! Each processor gets index from "startn" to "lastn"
  ! If the number of processors nproc exceeds the number of objects,
  ! the last nproc-ntodiv processors return startn = ntodiv+1 > lastn = ntodiv
  !
  USE mp, ONLY : mp_size, mp_rank
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm
  ! communicator
  INTEGER, INTENT(in) :: ntodiv
  ! index to be distributed
  INTEGER, INTENT(out):: startn, lastn
  ! indices for this processor: from startn to lastn
  !
  INTEGER :: me_comm, nproc_comm
  ! identifier of current processor
  ! number of processors
  !
  INTEGER :: ndiv, rest
  ! number of points per processor
  ! number of processors having one more points
  !
  nproc_comm = mp_size(comm)
  me_comm = mp_rank(comm)
  !
  rest = mod ( ntodiv, nproc_comm )
  ndiv = int( ntodiv / nproc_comm ) 
  !
  IF (rest > me_comm) THEN 
     startn =  me_comm    * (ndiv+1) + 1
     lastn  = (me_comm+1) * (ndiv+1) 
  ELSE
     startn=  me_comm    * ndiv + rest + 1
     lastn = (me_comm+1) * ndiv + rest 
  ENDIF

  RETURN

END SUBROUTINE divide

SUBROUTINE divide_all (comm, ntodiv, startn, lastn, counts, displs)
  !-----------------------------------------------------------------------
  !
  ! Given "ntodiv" objects, distribute index across a group of processors
  ! belonging to communicator "comm"
  ! Each processor gets index from "startn" to "lastn"
  ! If the number of processors nproc exceeds the number of objects,
  ! the last nproc-ntodiv processors return startn = ntodiv+1 > lastn = ntodiv
  !
  USE mp, ONLY : mp_size, mp_rank
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm
  ! communicator
  INTEGER, INTENT(in) :: ntodiv
  ! index to be distributed
  INTEGER, INTENT(out):: startn, lastn
  ! indices for this processor: from startn to lastn
  INTEGER, INTENT(out):: counts(*), displs(*)
  ! indice counts and displacements of all ranks
  !
  INTEGER :: me_comm, nproc_comm
  ! identifier of current processor
  ! number of processors
  !
  INTEGER :: ndiv, rest
  ! number of points per processor
  ! number of processors having one more points
  INTEGER :: ip
  !
  nproc_comm = mp_size(comm)
  me_comm = mp_rank(comm)
  !
  rest = mod ( ntodiv, nproc_comm )
  ndiv = int( ntodiv / nproc_comm )
  !
  DO ip = 1, nproc_comm
     IF (rest >= ip) THEN
        counts(ip) = ndiv + 1
        displs(ip)  = (ip-1) * (ndiv+1)
     ELSE
        counts(ip) = ndiv
        displs(ip)  = (ip-1) * ndiv + rest
     ENDIF
  ENDDO
  ! seting startn and lastn
  startn =  displs(me_comm+1) + 1
  lastn = displs(me_comm+1) + counts(me_comm+1)

  RETURN

END SUBROUTINE divide_all
