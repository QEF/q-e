!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, &
     irt, irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, &
     npert, nirr, gi, gimq, iverbosity, modenum)
  !---------------------------------------------------------------------
  !
  !    This routine computes the symmetry matrix of the mode defined
  !    by modenum. It sets also the modes u for all the other
  !    representation
  !
  !
  !
  USE kinds, only : DP
  USE constants, ONLY: tpi
  implicit none
  !
  !   first the dummy variables
  !

  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, modenum, npert (3 * nat), irgq (48), nsymq, irotmq, &
       nirr, npertx
  ! input: the number of atoms
  ! input: the number of symmetries
  ! input: the symmetry matrices
  ! input: the inverse of each matrix
  ! input: the rotated of each atom
  ! input: write control
  ! input: the mode to be done
  ! output: the dimension of each represe
  ! output: the small group of q
  ! output: the order of the small group
  ! output: the symmetry sending q -> -q+
  ! output: the number of irr. representa

  real(DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3), &
       gi (3, 48), gimq (3)
  ! input: the q point
  ! input: the R associated to each tau
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! output: [S(irotq)*q - q]
  ! output: [S(irotmq)*q + q]

  complex(DP) :: u(3*nat, 3*nat), t(npertx, npertx, 48, 3*nat),&
       tmq (npertx, npertx, 3 * nat)
  ! output: the pattern vectors
  ! output: the symmetry matrices
  ! output: the matrice sending q -> -q+G
  logical :: minus_q
  ! output: if true one symmetry send q -> -q+G
  !
  !   here the local variables
  !
  integer :: na, imode, jmode, ipert, jpert, nsymtot, imode0, irr, &
       ipol, jpol, isymq, irot, sna
  ! counters and auxilary variables

  real(DP) :: modul, arg
  ! the modulus of the mode
  ! the argument of the phase

  complex(DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase
  ! one pattern
  ! the rotated of one pattern
  ! the phase factor

  logical :: lgamma
  ! if true gamma point
  !
  !   Allocate the necessary quantities
  !
  lgamma = (xq (1) == 0.d0 .and. xq (2) == 0.d0 .and. xq (3) == 0.d0)
  !
  !   find the small group of q
  !
  call smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, minus_q, gi, gimq)
  !
  !    set the modes to be done
  !
  u (:, :) = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     u (imode, imode) = (1.d0, 0.d0)
  enddo
  !
  !  Here we count the irreducible representations and their dimensions
  !
  nirr = 3 * nat
  ! initialization
  npert (:) = 1
  !
  !   And we compute the matrices which represent the symmetry transformat
  !   in the basis of the displacements
  !

  t(:, :, :, :) = (0.d0, 0.d0)
  tmq (:, :, :) = (0.d0, 0.d0)
  if (minus_q) then
     nsymtot = nsymq + 1
  else
     nsymtot = nsymq
  endif

  do isymq = 1, nsymtot
     if (isymq.le.nsymq) then
        irot = irgq (isymq)
     else
        irot = irotmq
     endif
     imode0 = 0
     do irr = 1, nirr
        do ipert = 1, npert (irr)
           imode = imode0 + ipert
           do na = 1, nat
              do ipol = 1, 3
                 jmode = 3 * (na - 1) + ipol
                 wrk_u (ipol, na) = u (jmode, imode)
              enddo
           enddo
           !
           !     transform this pattern to crystal basis
           !
           do na = 1, nat
              call trnvecc (wrk_u (1, na), at, bg, - 1)
           enddo
           !
           !     the patterns are rotated with this symmetry
           !
           wrk_ru(:,:) = (0.d0, 0.d0)
           do na = 1, nat
              sna = irt (irot, na)
              arg = 0.d0
              do ipol = 1, 3
                 arg = arg + xq (ipol) * rtau (ipol, irot, na)
              enddo
              arg = arg * tpi
              if (isymq == nsymtot .and. minus_q) then
                 fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
              else
                 fase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
              endif
              do ipol = 1, 3
                 do jpol = 1, 3
                    wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + fase * &
                         s (jpol, ipol, irot) * wrk_u (jpol, na)
                 enddo
              enddo
           enddo
           !
           !    Transform back the rotated pattern
           !
           do na = 1, nat
              call trnvecc (wrk_ru (1, na), at, bg, 1)
           enddo
           !
           !     Computes the symmetry matrices on the basis of the pattern
           !
           do jpert = 1, npert (irr)
              imode = imode0 + jpert
              do na = 1, nat
                 do ipol = 1, 3
                    jmode = ipol + (na - 1) * 3
                    if (isymq == nsymtot .and. minus_q) then
                       tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + &
                            CONJG(u (jmode, imode) * wrk_ru (ipol, na) )
                    else
                       t (jpert, ipert, irot, irr) = t (jpert, ipert, irot, irr) &
                            + CONJG(u (jmode, imode) ) * wrk_ru (ipol, na)
                    endif
                 enddo
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo

  enddo
  !      WRITE( stdout,*) 'nsymq',nsymq
  !      do isymq=1,nsymq
  !        irot=irgq(isymq)
  !        WRITE( stdout,'("t(1,1,irot,modenum)",i5,2f10.5)')
  !     +                 irot,t(1,1,irot,modenum)
  !      enddo
  return
end subroutine set_irr_mode
