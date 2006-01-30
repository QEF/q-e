!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine n_plane_waves (ecutwfc, tpiba2, nks, nkstot, xk, g, ngm, &
     npwx, ngk)
  !-----------------------------------------------------------------------
  !
  ! Find number of plane waves for each k-point
  !
  USE kinds
  implicit none
  ! input
  integer :: nks, nkstot, ngm
  real(DP) :: ecutwfc, tpiba2, xk (3, nks), g (3, ngm)
  ! output
  integer :: npwx, ngk (nkstot)
  ! local
  integer :: nk, ng
  real(DP) :: q2
  !
  npwx = 0
  do nk = 1, nks
     ngk (nk) = 0
     do ng = 1, ngm
        q2 = (xk (1, nk) + g (1, ng) ) **2 + (xk (2, nk) + g (2, ng) ) ** &
             2 + (xk (3, nk) + g (3, ng) ) **2
        if (q2 <= ecutwfc / tpiba2) then
           !
           ! here if |k+G|^2 <= Ecut increase the number of G inside the sphere
           !
           ngk (nk) = ngk (nk) + 1
        else
           if (sqrt (g (1, ng) **2 + g (2, ng) **2 + g (3, ng) **2) &
                .gt.sqrt (xk (1, nk) **2 + xk (2, nk) **2 + xk (3, nk) **2) &
                + sqrt (ecutwfc / tpiba2) ) goto 100
           !
           ! if |G| > |k| + sqrt(Ecut)  stop search
           !
        endif
     enddo
100  npwx = max (npwx, ngk (nk) )
  enddo
  if (npwx <= 0) call errore ('n_plane_waves', 'No PWs !', 1)
  return
end subroutine n_plane_waves
