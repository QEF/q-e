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
     irt, irgq, nsymq, minus_q, irotmq, t, tmq, u, npert, nirr, gi, &
     gimq, iverbosity, modenum)
  !---------------------------------------------------------------------
  !
  !    This routine computes the symmetry matrix of the mode defined
  !    by modenum. It sets also the modes u for all the other
  !    representation
  !
  !
  !
#include "machine.h"
  use parameters, only : DP
#ifdef PARA
  use para  
#endif
  implicit none  
#ifdef PARA
  include 'mpif.h'  
#endif
  !
  !   first the dummy variables
  !

  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, modenum, npert (3 * nat), irgq (48), nsymq, irotmq, &
       nirr
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

  real(kind=DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3), &
       gi (3, 48), gimq (3)
  ! input: the q point
  ! input: the R associated to each tau
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! output: [S(irotq)*q - q]
  ! output: [S(irotmq)*q + q]

  complex(kind=DP) :: u (3 * nat, 3 * nat), t (3, 3, 48, 3 * nat), &
       tmq (3, 3, 3 * nat)
  ! output: the pattern vectors
  ! output: the symmetry matrices
  ! output: the matrice sending q -> -q+G
  logical :: minus_q  
  ! output: if true one symmetry send q -
  !
  !   here the local variables
  !
  real(kind=DP) :: tpi  

  parameter (tpi = 2.0d0 * 3.14159265358979d0)  


  integer :: na, imode, jmode, ipert, jpert, nsymtot, imode0, irr, &
       ipol, jpol, isymq, irot, sna
  ! counter on atoms
  ! counter on modes
  ! counter on modes
  ! counter on perturbations
  ! counter on perturbations
  ! total number of symmetries
  ! auxiliry variable for mode counting
  ! counter on irreducible representation
  ! counter on polarizations
  ! counter on polarizations
  ! counter on symmetries
  ! counter on rotations
  ! the rotated atom

  real(kind=DP) :: modul, arg  
  ! the modulus of the mode
  ! the argument of the phase

  complex(kind=DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase  
  ! one pattern
  ! the rotated of one pattern
  ! the phase factor

  logical :: lgamma  
  ! if true gamma point
  !
  !   Allocate the necessary quantities
  !
  lgamma = (xq (1) .eq.0.d0.and.xq (2) .eq.0.d0.and.xq (3) .eq.0.d0)  
  !
  !   find the small group of q
  !
  call smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, minus_q, &
       gi, gimq)
  !
  !    set the modes to be done
  !
  call setv (18 * nat * nat, 0.d0, u, 1)  
  do imode = 1, 3 * nat  
     u (imode, imode) = (1.d0, 0.d0)  
  enddo
  !
  !  Here we count the irreducible representations and their dimensions
  !
  nirr = 3 * nat  
  do imode = 1, 3 * nat  
     ! initialization
     npert (imode) = 1  
  enddo
  !
  !   And we compute the matrices which represent the symmetry transformat
  !   in the basis of the displacements
  !
  call setv (2 * 3 * 3 * 48 * 3 * nat, 0.d0, t, 1)  
  call setv (2 * 3 * 3 * 3 * nat, 0.d0, tmq, 1)  
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
           call setv (2 * 3 * nat, 0.d0, wrk_ru, 1)  
           do na = 1, nat  
              sna = irt (irot, na)  
              arg = 0.d0  
              do ipol = 1, 3  
                 arg = arg + xq (ipol) * rtau (ipol, irot, na)  
              enddo
              arg = arg * tpi  
              if (isymq.eq.nsymtot.and.minus_q) then  
                 fase = DCMPLX (cos (arg), sin (arg) )  
              else  
                 fase = DCMPLX (cos (arg), - sin (arg) )  
              endif
              do ipol = 1, 3  
                 do jpol = 1, 3  
                    wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + s (jpol, ipol, irot) &
                         * wrk_u (jpol, na) * fase
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
                    if (isymq.eq.nsymtot.and.minus_q) then  
                       tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + conjg (u ( &
                            jmode, imode) * wrk_ru (ipol, na) )
                    else  
                       t (jpert, ipert, irot, irr) = t (jpert, ipert, irot, irr) &
                            + conjg (u (jmode, imode) ) * wrk_ru (ipol, na)
                    endif
                 enddo
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)  
     enddo

  enddo
  !      write(6,*) 'nsymq',nsymq
  !      do isymq=1,nsymq
  !        irot=irgq(isymq)
  !        write(6,'("t(1,1,irot,modenum)",i5,2f10.5)')
  !     +                 irot,t(1,1,irot,modenum)
  !      enddo
  return  
end subroutine set_irr_mode
