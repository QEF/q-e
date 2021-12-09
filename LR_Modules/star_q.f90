!
! Copyright (C) 2001-2020 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine star_q1(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, verbosity,&
     t_rev )
  !-----------------------------------------------------------------------
  ! FIXME: MERGE WITH STAR_Q
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  real(DP), parameter :: accep=1.e-5_dp

  integer, intent(in) :: nsym, s (3, 3, 48), invs(48)
  integer, intent(in) :: t_rev(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  ! t_rev: if a given simmetry operation needs time reversal
  real(DP), intent(in) :: xq (3), at (3, 3), bg (3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  integer, intent(out) :: nq, isq (48), imq
  ! nq  : degeneracy of the star of q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  real(DP), intent(out) :: sxq (3, 48)
  ! list of vectors in the star of q
  logical, intent(in) :: verbosity
  ! if true prints several messages.
  !
  integer :: nsq (48), isym, ism1, iq, i
  ! number of symmetry ops. of bravais lattice
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counters
  real(DP) :: saq (3, 48), aq (3), raq (3), zero (3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect

  logical, external :: eqvect
  ! function used to compare two vectors
  !
  !
  zero(:) = 0.d0
  saq(:, :) = 0.0d0 
  !
  ! go to  crystal coordinates
  !
  do i = 1, 3
     aq(i) = xq(1) * at(1,i) + xq(2) * at(2,i) + xq(3) * at(3,i)
  enddo
  !
  ! create the list of rotated q
  !
  do i = 1, 48
     nsq (i) = 0
     isq (i) = 0
  enddo
  nq = 0
  do isym = 1, nsym
     ism1 = invs (isym)
     do i = 1, 3
        raq (i) = s (i, 1, ism1) * aq (1) &
                + s (i, 2, ism1) * aq (2) &
                + s (i, 3, ism1) * aq (3)
     enddo
     IF (t_rev(isym)==1) raq = -raq
     do i = 1, 3
        sxq (i, 48) = bg (i, 1) * raq (1) &
                    + bg (i, 2) * raq (2) &
                    + bg (i, 3) * raq (3)
     enddo
     do iq = 1, nq
        if (eqvect (raq, saq (1, iq), zero, accep) ) then
           isq (isym) = iq
           nsq (iq) = nsq (iq) + 1
        endif
     enddo
     if (isq (isym) == 0) then
        nq = nq + 1
        nsq (nq) = 1
        isq (isym) = nq
        saq(:,nq) = raq(:)
        do i = 1, 3
           sxq (i, nq) = bg (i, 1) * saq (1, nq) &
                       + bg (i, 2) * saq (2, nq) &
                       + bg (i, 3) * saq (3, nq)
        enddo
     endif
  enddo
  !
  ! set imq index if needed and check star degeneracy
  !
  raq (:) = - aq(:)
  imq = 0
  do iq = 1, nq
     if (eqvect (raq, saq (1, iq), zero, accep) ) imq = iq
     if (nsq(iq)*nq /= nsym) call errore ('star_q', 'wrong degeneracy', iq)
  enddo
  !
  ! writes star of q
  !
  IF (verbosity) THEN
  WRITE( stdout, * )
  WRITE( stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
  WRITE( stdout, '(5x,a)') 'List of q in the star:'
  WRITE( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
  if (imq == 0) then
     WRITE( stdout, '(5x,a)') 'In addition there is the -q list: '
     WRITE( stdout, '(7x,i4,3f14.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
  endif
  ENDIF
  return
end subroutine star_q1
!-----------------------------------------------------------------------
subroutine star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, verbosity)
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  real(DP), parameter :: accep=1.e-5_dp

  integer, intent(in) :: nsym, s (3, 3, 48), invs(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  real(DP), intent(in) :: xq (3), at (3, 3), bg (3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  integer, intent(out) :: nq, isq (48), imq
  ! nq  : degeneracy of the star of q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  real(DP), intent(out) :: sxq (3, 48)
  ! list of vectors in the star of q
  logical, intent(in) :: verbosity
  ! if true prints several messages.
  !
  integer :: nsq (48), isym, ism1, iq, i, t_rev(48)
  ! number of symmetry ops. of bravais lattice
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counters
  ! t_rev variable, says if a given simmetry operation
  ! needs the time reversal to be a symmetry of the crystal.
  real(DP) :: saq (3, 48), aq (3), raq (3), zero (3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect

  logical, external :: eqvect
  ! function used to compare two vectors
  !
  !
  zero(:) = 0.d0
  saq(:, :) = 0.0d0 
  !
  ! go to  crystal coordinates
  !
  do i = 1, 3
     aq(i) = xq(1) * at(1,i) + xq(2) * at(2,i) + xq(3) * at(3,i)
  enddo
  !
  ! create the list of rotated q
  !
  do i = 1, 48
     nsq (i) = 0
     isq (i) = 0
  enddo
  nq = 0
  do isym = 1, nsym
     ism1 = invs (isym)
     do i = 1, 3
        raq (i) = s (i, 1, ism1) * aq (1) &
                + s (i, 2, ism1) * aq (2) &
                + s (i, 3, ism1) * aq (3)
     enddo
     do i = 1, 3
        sxq (i, 48) = bg (i, 1) * raq (1) &
                    + bg (i, 2) * raq (2) &
                    + bg (i, 3) * raq (3)
     enddo
     do iq = 1, nq
        if (eqvect (raq, saq (1, iq), zero, accep) ) then
           isq (isym) = iq
           nsq (iq) = nsq (iq) + 1
        endif
     enddo
     if (isq (isym) == 0) then
        nq = nq + 1
        nsq (nq) = 1
        isq (isym) = nq
        saq(:,nq) = raq(:)
        do i = 1, 3
           sxq (i, nq) = bg (i, 1) * saq (1, nq) &
                       + bg (i, 2) * saq (2, nq) &
                       + bg (i, 3) * saq (3, nq)
        enddo
     endif
  enddo
  !
  ! set imq index if needed and check star degeneracy
  !
  raq (:) = - aq(:)
  imq = 0
  do iq = 1, nq
     if (eqvect (raq, saq (1, iq), zero, accep) ) imq = iq
     if (nsq(iq)*nq /= nsym) call errore ('star_q', 'wrong degeneracy', iq)
  enddo
  !
  ! writes star of q
  !
  IF (verbosity) THEN
  WRITE( stdout, * )
  WRITE( stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
  WRITE( stdout, '(5x,a)') 'List of q in the star:'
  WRITE( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
  if (imq == 0) then
     WRITE( stdout, '(5x,a)') 'In addition there is the -q list: '
     WRITE( stdout, '(7x,i4,3f14.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
  endif
  ENDIF
  return
end subroutine star_q
