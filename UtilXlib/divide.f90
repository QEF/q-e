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
  ! Divide ntodiv poins across processors belonging to communicator comm 
  ! Each processor gets points from startn to lastn
  !
#if defined(__MPI)
  !
  USE mp, ONLY : mp_size, mp_rank
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm
  INTEGER, INTENT(in) :: ntodiv
  INTEGER, INTENT(out):: startn, lastn
  !
  INTEGER :: me_comm, nproc_comm
  !
  INTEGER :: nb, resto, idx, ip
  ! number of bands per processor
  ! one additional band if me_pool+1 <= resto
  ! counter on bands
  ! counter on processors
  !
  nproc_comm = mp_size(comm)
  me_comm = mp_rank(comm)
  !
  nb = ntodiv / nproc_comm
  resto = ntodiv - nb * nproc_comm
  idx = 0
  DO ip = 1, nproc_comm
     IF (ip <= resto) THEN
        IF (me_comm+1 == ip) THEN
           startn = idx + 1
           lastn = startn + nb
        ENDIF
        idx = idx + nb + 1
     ELSE
        IF (me_comm+1 == ip) THEN
           startn = idx + 1
           lastn = startn + nb - 1
        ENDIF
        idx = idx + nb
     ENDIF
  ENDDO
#else

  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm
  INTEGER, INTENT(in) :: ntodiv
  INTEGER, INTENT(out):: startn, lastn
 
  startn = 1
  lastn = ntodiv

#endif
  RETURN

END SUBROUTINE divide

!
!-----------------------------------------------------------------------
SUBROUTINE divide2 (comm1, comm2, ntodiv, startn, lastn)
  !-----------------------------------------------------------------------
  ! Divide ntodiv points across processors belonging to two communicators 
  ! comm1 and comm2. The final quantity must be collected among the two
  ! Each processor gets points from startn to lastn
  !
#if defined(__MPI)
  !
  USE mp, ONLY : mp_size, mp_rank
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm1
  INTEGER, INTENT(in) :: comm2
  INTEGER, INTENT(in) :: ntodiv
  INTEGER, INTENT(out):: startn, lastn

  INTEGER :: ntodiv1, start_n1, end_n1, start_n2, end_n2
  !
  CALL divide( comm1, ntodiv, start_n1, end_n1 )
  ntodiv1=end_n1-start_n1+1
  CALL divide( comm2, ntodiv1, start_n2, end_n2 )
  startn=start_n1+start_n2-1
  lastn=start_n1+end_n2-1

#else

  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: comm1, comm2
  INTEGER, INTENT(in) :: ntodiv
  INTEGER, INTENT(out):: startn, lastn
 
  startn = 1
  lastn = ntodiv

#endif
  RETURN

END SUBROUTINE divide2

