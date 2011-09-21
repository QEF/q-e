!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine symd0rho (npertx, nper, irr, d0rho, s, ftau, nsymq, &
     irgq, t, nat, nr1, nr2, nr3, nr1x, nr2x, nr3x)
  !---------------------------------------------------------------------
  !  symmetrizes q=0 drho
  !
  !
  USE kinds, only : DP
  implicit none
  integer :: nper, irr, s (3, 3, 48), ftau (3, 48), nsymq, irgq (48) &
       , nat, nr1, nr2, nr3, nr1x, nr2x, nr3x, npertx
  ! nper: the number of perturbations
  ! irr: the representation under consideration

  complex (DP) :: d0rho (nr1x, nr2x, nr3x, nper),        &
       t (npertx, npertx, 48, 3 * nat)
  ! charge variation to symmetrize

  integer :: ri, rj, rk, i, j, k, ipert, jpert, isym, irot
  !  ri, rj, rk: rotated points
  !  counters

  complex (DP), allocatable :: aux1 (:,:,:,:)
  ! the symmetrized charge


  call start_clock ('symd0rho')
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           do ipert = 1, nper
              d0rho (i, j, k, ipert) =  DBLE (d0rho (i, j, k, ipert) )
           enddo
        enddo
     enddo
  enddo

  if (nsymq == 1) return

  allocate  (aux1( nr1x, nr2x, nr3x, nper))
  !
  ! Here we symmetrize with respect to the group
  !
  aux1 (:,:,:,:) = (0.d0, 0.d0)
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           do isym = 1, nsymq
              irot = irgq (isym)
              ri = s (1, 1, irot) * (i - 1) + s (2, 1, irot) * (j - 1) + &
                   s (3, 1, irot) * (k - 1) - ftau (1, irot)
              ri = mod (ri, nr1) + 1
              if (ri < 1) ri = ri + nr1
              rj = s (1, 2, irot) * (i - 1) + s (2, 2, irot) * (j - 1) + &
                   s (3, 2, irot) * (k - 1) - ftau (2, irot)
              rj = mod (rj, nr2) + 1
              if (rj < 1) rj = rj + nr2
              rk = s (1, 3, irot) * (i - 1) + s (2, 3, irot) * (j - 1) + &
                   s (3, 3, irot) * (k - 1) - ftau (3, irot)
              rk = mod (rk, nr3) + 1
              if (rk < 1) rk = rk + nr3
              do ipert = 1, nper
                 do jpert = 1, nper
                    aux1 (i, j, k, ipert) = aux1 (i, j, k, ipert) + &
                         t(jpert, ipert, irot, irr) * d0rho (ri, rj, rk, jpert)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  d0rho (:,:,:,:) = aux1 (:,:,:,:) / DBLE (nsymq)

  deallocate (aux1)

  call stop_clock ('symd0rho')
  return
end subroutine symd0rho
