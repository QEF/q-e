!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_weight (wgg)
  !-----------------------------------------------------------------------
  !
  !     This routine implements Eq.15 (B) of the notes. It computes the
  !     weight to give to the v,v' terms in the orthogonality term
  !

  use pwcom
  USE kinds, only : DP
  use phcom
  implicit none

  real(DP) :: wgg (nbnd, nbnd, nksq)
  ! output: the weights

  integer :: ik, ikk, ikq, ibnd, jbnd
  ! counters
  real(DP) :: wg1, wg2, theta
  ! auxiliary variables
  real(DP), external :: wgauss
  real(DP), parameter :: eps = 1.0d-12
  !
  !     the weights are computed for each k point ...
  !
  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
        ikq = ik
     else
        ikk = 2 * ik - 1
        ikq = ikk + 1
     endif
     !
     !     each band v ...
     !
     do ibnd = 1, nbnd
        if (wk (ikk) .eq.0.d0) then
           wg1 = 0.d0
        else
           wg1 = wg (ibnd, ikk) / wk (ikk)
        endif
        !
        !     and each band v' ...
        !
        do jbnd = 1, nbnd
           if (lgauss) then
              theta = wgauss ( (et (jbnd,ikq) - et (ibnd,ikk) ) / degauss, 0)
              wg2 = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
           else
              IF (et (jbnd,ikq) > et (ibnd,ikk)) THEN
                 theta = 1.0d0
              ELSE
                 theta = 0.d0
              ENDIF
              IF (ABS(et (jbnd,ikq) - et (ibnd,ikk)) < 1.d-8) theta=0.5d0
              if (wk (ikk) .le.eps) then
                 wg2 = 0.d0
              else
                 wg2 = wg (jbnd, ikk) / wk (ikk)
              endif
           endif
           wgg (ibnd, jbnd, ik) = wg1 * (1.d0 - theta) + wg2 * theta
        enddo
     enddo
     !         do ibnd=1,nbnd
     !            do jbnd=1,nbnd
     !               WRITE( stdout,'(3i5,f20.10)') ibnd, jbnd, ik,wgg(ibnd,jbnd,ik)
     !            enddo
     !         enddo

  enddo
  !      call stop_ph(.true.)
  return

end subroutine compute_weight
