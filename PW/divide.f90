!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine divide (ntodiv, startn, lastn)
  !-----------------------------------------------------------------------
  ! Divide the bands among processors in the phonon code
  !
  USE mp_global, ONLY : me_pool, nproc_pool
  !
#ifdef __PARA
  !
  implicit none

  integer :: ntodiv, startn, lastn
  ! input: the number to divide
  ! output: the first band of this processor
  ! output: the last band of this processor

  integer :: nb, resto, index, ip
  ! number of bands per processor
  ! one additional band if me_pool+1 <= resto
  ! counter on bands
  ! counter on processors
  nb = ntodiv / nproc_pool
  resto = ntodiv - nb * nproc_pool
  index = 0
  do ip = 1, nproc_pool
     if (ip.le.resto) then
        if (me_pool+1.eq.ip) then
           startn = index + 1
           lastn = startn + nb
        endif
        index = index + nb + 1
     else
        if (me_pool+1.eq.ip) then
           startn = index + 1
           lastn = startn + nb - 1
        endif
        index = index + nb
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

