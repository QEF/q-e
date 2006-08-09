!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine star_q (xq, at, bg, ibrav, symm_type, nat, tau, ityp, &
     nr1, nr2, nr3, nsym, s, invs, irt, rtau, nq, sxq, isq, imq, noinv, &
     modenum)
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! and return their list along with the symmetry ops. needed to obtain
  ! them.
  ! NB: output values of symmetry arrays (nsym, s, rtau, irt) are those
  ! appropriate to the crystal symmetry (not to the small-qroup of q).
  ! User is responsible for calling this routine with different array-name
  ! if information for both symmetry-groups needs to be kept
  !
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !-input variables
  integer :: ibrav, nat, ityp (nat), modenum, nr1, nr2, nr3
  ! input: bravais lattice index
  ! input: number of atom
  ! input: atomic type
  ! input: the mode to be done
  ! input: fft grid dimensions
  real(DP) :: xq (3), at (3, 3), bg (3, 3), tau (3, nat)
  ! input: q vector
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  ! input: coordinates of atomic positions

  character (len=9) :: symm_type
  ! input: 'cubic' or 'hexagonal' when ibrav=0

  logical :: noinv
  ! input: if true eliminates symmetries z <-> -z
  !-output variables
  integer :: nsym, s (3, 3, 48), invs (48), irt (48, nat), nq, isq (48), imq
  ! output: number of symmetry operations
  ! output: the first nq matrices are those that generate the star of q
  !         starting from it
  ! output: list of inverse operation indices
  ! output: for each atom gives the rotated atom
  ! output: degeneracy of the star of q
  ! output: index of q in the star for a given sym
  ! output: index of -q in the star (0 if not present)

  real(DP) :: rtau (3, 48, nat), sxq (3, 48)
  ! output: for each atom and rotation gives the R vector involved
  ! output: list of vectors in the star of q
  !
  ! Local variables
  !
  integer :: nsq (48), ftau(3,48), nrot, isym, jsym, ism1, table (48, 48), &
       iq, i, j, nks, npk
  ! number of symmetry ops. of bravais lattice.
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! group table
  ! counter on q-vectors
  ! generic counter
  ! number of dummy k-points
  ! maximum allowed number of dummy k-points
  integer :: t_rev(48) = 0
  ! for magnetic symetries - not actually used
  real(DP) :: saq (3, 48), aq (3), raq (3), xk0 (3), wk(1), zero (3), &
       mdum(3,nat)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! dummy k-points list
  ! a zero vector: used in eqvect and as dummy q-vector in sgama

  logical :: invsym, minus_q, nosym, sym (48)
  ! .t. if the crystal has inversion
  ! dummy output from sgama
  ! input for sgama

  character (len=45) :: sname (48)
  ! name of the rotat. part of each selected symmetry operation

  logical, external :: eqvect
  ! function used to compare two vectors
  !
  !  initialize dummy k-point list and zero vector
  !
  npk = 1
  nks = 1
  wk(:) = 1.d0
  xk0(:)= 0.d0
  zero(:) = 0.d0
  !
  !  generate transformation matrices for the bravais lattice
  !
  if (ibrav == 4 .or. ibrav == 5) then
     call hexsym (at, s, sname, nrot)
  elseif (ibrav >= 1 .and. ibrav <= 14) then
     call cubicsym (at, s, sname, nrot)
  elseif (ibrav == 0) then
     if (symm_type == 'cubic') call cubicsym (at, s, sname, nrot)
     if (symm_type == 'hexagonal') call hexsym (at, s, sname, nrot)
  else
     call errore ('star_q', 'wrong ibrav', 1)
  endif

  if (noinv) then
     jsym = 0
     do isym = 1, nrot
        if ( s (1, 3, isym) == 0 .and. s (3, 1, isym) == 0 .and. &
             s (2, 3, isym) == 0 .and. s (3, 2, isym) == 0 .and. &
             s (3, 3, isym) == 1) then
           jsym = jsym + 1
           do i = 1, 3
              do j = 1, 3
                 s (i, j, jsym) = s (i, j, isym)
              enddo
           enddo
           sname (jsym) = sname (isym)
        endif
     enddo
     nrot = jsym
  endif
  !
  ! extract from it the crystal symmetry group by calling sgama
  !
  nosym = .false.
  call sgama (nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym, nr1, &
       nr2, nr3, irt, ftau, npk, nks, xk0, wk, invsym, minus_q, zero, &
       modenum, .false., .false., mdum)
  do isym = 1, nsym
     sym (isym) = .true.
  enddo
  call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !
  ! computes the inverse of each matrix
  !
  call multable (nsym, s, table)
  do isym = 1, nsym
     do jsym = 1, nsym
        if (table (isym, jsym) .eq.1) invs (isym) = jsym
     enddo
  enddo
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
        if (eqvect (raq, saq (1, iq), zero) ) then
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
     if (eqvect (raq, saq (1, iq), zero) ) imq = iq
     if (nsq(iq)*nq /= nsym) call errore ('star_q', 'wrong degeneracy', iq)
  enddo
  !
  ! writes star of q
  !
  WRITE( stdout, * )
  WRITE( stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
  WRITE( stdout, '(5x,a)') 'List of q in the star:'
  WRITE( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
  if (imq == 0) then
     WRITE( stdout, '(5x,a)') 'In addition there is the -q list: '
     WRITE( stdout, '(7x,i4,3f12.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
  endif
  return
end subroutine star_q
