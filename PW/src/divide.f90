!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine divide (ntodiv, startn, lastn, comm)
  !-----------------------------------------------------------------------
  ! Divide the bands among processors in the phonon code
  !
  USE mp, ONLY : mp_size, mp_rank
  !
#ifdef __MPI
  !
  implicit none
  !
  integer :: comm
  !
  integer :: ntodiv, startn, lastn
  ! input: the number to divide
  ! output: the first band of this processor
  ! output: the last band of this processor
  !
  integer :: me_comm, nproc_comm
  !
  integer :: nb, resto, idx, ip
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
  do ip = 1, nproc_comm
     if (ip.le.resto) then
        if (me_comm+1.eq.ip) then
           startn = idx + 1
           lastn = startn + nb
        endif
        idx = idx + nb + 1
     else
        if (me_comm+1.eq.ip) then
           startn = idx + 1
           lastn = startn + nb - 1
        endif
        idx = idx + nb
     endif
  enddo
#else

  integer :: ntodiv, startn, lastn
  ! input: the number to divide among proce
  ! output: the first band of this processo
  ! output: the last band of this processor
  startn = 1
  lastn = ntodiv
#endif
  return

end subroutine divide

