!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine setupkpoint (s, nrot, xk, wk, nks, npk, nk1, nk2, nk3, &
     k1, k2, k3, at, bg, tipo)
  !-----------------------------------------------------------------------
  !
  !     This routine generates the point of a monkhorst and pack mesh
  !     If tipo=1 it uses the standard generation,
  !     If tipo=2 it uses the mesh for hexagonal bravais lattice.
  !
  !
  use parameters, only : DP
  implicit none

  integer :: nrot, nks, npk, nk1, nk2, nk3, k1, k2, k3, tipo, s (3, &
       3, 48)
  ! input: the number of symmetries
  ! output: the number of k points
  ! input: the dimension of xk and wk
  ! input: the number of points in the three dire
  ! input: kpoint shift
  ! input: hexagonal or cubic type
  ! input: the symmetry matrices

  real(kind=DP) :: at (3, 3), bg (3, 3), xk (3, npk), wk (npk)
  ! input: direct and reciprocal lattice ve
  ! output: k points and weights

  integer :: npkk
  ! maximum number of unreduced points

  parameter (npkk = 100000)

  real(kind=DP) :: xknew (3, npkk), xkk (npkk), esort (npkk), d (3), &
       dnorm, dk (3), buffer (3), swap, ui, uj, uk, eps
  ! the unreduced points
  ! their modulus
  ! used for ordering
  ! used for ordering
  ! auxiliary
  ! save the rotated point
  ! used for exchange
  ! auxiliary
  ! a small parameter

  parameter (eps = 1e-5)

  integer :: start, igk (npkk), ikprova, irot, iswap, indsw, count, &
       contatore, l, i, j, k
  ! initial point in the list
  ! used to order the k points
  ! counters on k points and rotations
  ! used to exchange points
  ! used to count on k points
  ! counter on polarizations
  ! counter on k point

  integer :: kpoint

  logical :: fatto (npkk)
  ! if true this point has been found
  if (npk.gt.npkk) call error ('setupkpoint', 'npkk too small', &
       npkk)

  if ( (tipo.ne.1) .and. (tipo.ne.2) ) call error ('setupkpoint', &
       'wrong tipo', 1)
  call setv (npkk, 1.d20, esort, 1)
  !
  !     set d vector for unique ordering
  !
  d (1) = 0.25657642786d0
  d (2) = 0.35342818974d0
  d (3) = 0.56421652427d0
  dnorm = sqrt (d (1) * d (1) + d (2) * d (2) + d (3) * d (3) )
  d (1) = d (1) / dnorm
  d (2) = d (2) / dnorm

  d (3) = d (3) / dnorm
  !
  !       Here we generate the points with the MP method
  !
  count = 0
  do i = 1, nk1
     if (tipo.eq.1) then
        ui = (k1 + 2.d0 * i - nk1 - 1.d0) / (2.d0 * nk1)
     else
        ui = (k1 + i - 1.d0) / nk1
     endif
     do j = 1, nk2
        if (tipo.eq.1) then
           uj = (k2 + 2.d0 * j - nk2 - 1.d0) / (2.d0 * nk2)
        else
           uj = (k2 + j - 1.d0) / nk2
        endif
        do k = 1, nk3
           uk = (k3 + 2.d0 * k - nk3 - 1.d0) / (2.d0 * nk3)
           count = count + 1
           if (count.gt.npkk) call error ('setupkpt', 'mpmesh too large', &
                npkk)
           xknew (1, count) = ui
           xknew (2, count) = uj
           xknew (3, count) = uk
           call modulo2 (xknew (1, count), bg, xkk (count), d, esort (count) &
                )
           igk (count) = count
        enddo
     enddo
  enddo
  !
  !     Now order the point in increasing modulus order
  !
  call hpsort (count, esort, igk)
  !
  !     now order the k point
  !
  do i = 1, count - 1
62   indsw = igk (i)
     if (indsw.ne.i) then
        do l = 1, 3
           swap = xknew (l, indsw)
           xknew (l, indsw) = xknew (l, igk (indsw) )
           xknew (l, igk (indsw) ) = swap
        enddo
        iswap = igk (i)
        igk (i) = igk (iswap)
        igk (iswap) = iswap
        swap = xkk (i)
        xkk (i) = xkk (iswap)
        xkk (iswap) = swap
        goto 62
     endif
  enddo
  !
  !      The k points are divided in shells of equivalent points
  !
  do i = 1, npkk
     fatto (i) = .false.

  enddo
  contatore = 0
  do kpoint = 1, count
     if (.not. (fatto (kpoint) ) ) then
        !
        !     We found the first vector of a new shell. Now found all the equiva
        !     vectors
        !
        contatore = contatore+1
        !
        !    This k point has been found
        !
        start = kpoint
        fatto (kpoint) = .true.
        do l = 1, 3
           xk (l, contatore) = xknew (l, kpoint)
        enddo
        wk (contatore) = 1.d0
        !
        !       Ora applichiamo in successione tutte le simmetrie
        !
        do irot = 1, nrot
           call prodotto3dk (s (1, 1, irot), xknew (1, kpoint), buffer)
           !
           !             now we look in the xknew list if there is buffer
           !
           do ikprova = start, count
              if (.not. (fatto (ikprova) ) ) then
                 do l = 1, 3
                    dk (l) = abs (xknew (l, ikprova) - buffer (l) )
                 enddo
                 if ( (abs (dk (1) - int (dk (1) + eps) ) .lt.2.d0 * eps) &
                      .and. (abs (dk (2) - int (dk (2) + eps) ) .lt.2.d0 * eps) &
                      .and. (abs (dk (3) - int (dk (3) + eps) ) .lt.2.d0 * eps) ) &
                      then
                    !
                    !   we have found the equivalent vector in the list
                    !
                    fatto (ikprova) = .true.
                    wk (contatore) = wk (contatore) + 1.d0
                    goto 90
                 endif
              endif
           enddo
90         continue
        enddo
     endif
  enddo
  nks = contatore
  call cryst_to_cart (nks, xk, bg, 1)
  return

end subroutine setupkpoint
!
!-----------------------------------------------------------------------
subroutine modulo2 (vect, bg, modulo, d, esort)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the modulus of the vector vect
  !     written in the bg basis. bg could be non orthogonal
  !
  use parameters, only : DP
  implicit none

  real(kind=DP) :: vect (3), bg (3, 3), d (3), esort, modulo
  ! input: the vector
  ! input: the basis
  ! input: direction for ordering
  ! output: the projected modulus
  ! output: the modulus

  real(kind=DP) :: buffer (3)
  ! auxiliary


  integer :: l
  modulo = 0.d0
  do l = 1, 3
     buffer (l) = vect (1) * bg (l, 1) + vect (2) * bg (l, 2) + vect ( &
          3) * bg (l, 3)
     modulo = modulo + buffer (l) * buffer (l)
     if (modulo.gt.1.d-8) then
        esort = 1.d4 * modulo + (buffer (1) * d (1) + buffer (2) &
             * d (2) + buffer (3) * d (3) ) / sqrt (modulo)
     else
        esort = 0.d0
     endif

  enddo
  return

end subroutine modulo2
!
!-----------------------------------------------------------------------
subroutine prodotto3dk (a, v, w)
  !-----------------------------------------------------------------------
  !
  !       This subrutine computes w=A v where A is a 3*3 matrix
  !
  use parameters, only : DP
  implicit none
  integer :: a (3, 3)
  real(kind=DP) :: v (3)

  real(kind=DP) :: w (3)

  integer :: l
  do l = 1, 3
     w (l) = a (l, 1) * v (1) + a (l, 2) * v (2) + a (l, 3) * v (3)

  enddo
  return
end subroutine prodotto3dk
