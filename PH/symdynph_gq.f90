!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symdynph_gq (xq, phi, s, invs, rtau, irt, irgq, nsymq, &
     nat, irotmq, minus_q)
  !-----------------------------------------------------------------------
  !
  !     This routine receives as input an unsymmetrized dynamical
  !     matrix expressed on the crystal axes and imposes the symmetry
  !     of the small group of q. Furthermore it imposes also the symmetry
  !     q -> -q+G if present.
  !
  !
#include "machine.h"
  use parameters, only : DP
  implicit none  
  !
  !    The dummy variables
  !
  integer :: nat, s (3, 3, 48), irt (48, nat), irgq (48), invs (48), &
       nsymq, irotmq
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each vector
  ! input: the small group of q
  ! input: the inverse of each matrix
  ! input: the order of the small gro
  ! input: the rotation sending q ->
  real(kind=DP) :: xq (3), rtau (3, 48, nat)  
  ! input: the q point
  ! input: the R associated at each t

  logical :: minus_q  
  ! input: true if a symmetry q->-q+G
  complex(kind=DP) :: phi (3, 3, nat, nat)  
  ! inp/out: the matrix to symmetrize
  !
  !     One parameter
  !
  real(kind=DP) :: tpi  
  parameter (tpi = 2.0d0 * 3.14159265358979d0)  
  !
  !   and the local variables
  !

  integer :: isymq, sna, snb, irot, na, nb, ipol, jpol, lpol, kpol, &
       iflb (nat, nat)
  ! counter on symmetries
  ! the rotated of the a atom
  ! the rotated of the b atom
  ! counter on rotations
  ! counter on atoms
  ! counter on atoms
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! used to account for symmetrized e

  real(kind=DP) :: arg  
  ! the argument of the phase

  complex(kind=DP) :: phip (3, 3, nat, nat), work (3, 3), fase, faseq ( &
       48)
  ! working space
  ! working space
  ! the phase factor
  ! the phases for each symmetry
  !
  !    We start by imposing hermiticity
  !
  do na = 1, nat  
     do nb = 1, nat  
        do ipol = 1, 3  
           do jpol = 1, 3  
              phi (ipol, jpol, na, nb) = 0.5d0 * (phi (ipol, jpol, na, nb) &
                   + conjg (phi (jpol, ipol, nb, na) ) )
              phi (jpol, ipol, nb, na) = conjg (phi (ipol, jpol, na, nb) )  
           enddo
        enddo
     enddo
  enddo
  !
  !    If no other symmetry is present we quit here
  !

  if ( (nsymq.eq.1) .and. (.not.minus_q) ) return  
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
  if (minus_q) then  
     do na = 1, nat  
        do nb = 1, nat  
           do ipol = 1, 3  
              do jpol = 1, 3  
                 call setv (18, 0.d0, work, 1)  
                 sna = irt (irotmq, na)  
                 snb = irt (irotmq, nb)  
                 arg = 0.d0  
                 do kpol = 1, 3  
                    arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) - rtau (kpol, &
                         irotmq, nb) ) )
                 enddo
                 arg = arg * tpi  
                 fase = DCMPLX (cos (arg), sin (arg) )  
                 do kpol = 1, 3  
                    do lpol = 1, 3  
                       work (ipol, jpol) = work (ipol, jpol) + s (ipol, kpol, irotmq) &
                            * s (jpol, lpol, irotmq) * phi (kpol, lpol, sna, snb) * fase
                    enddo
                 enddo
                 phip (ipol, jpol, na, nb) = (phi (ipol, jpol, na, nb) + conjg ( &
                      work (ipol, jpol) ) ) * 0.5d0
              enddo
           enddo
        enddo
     enddo
     call ZCOPY (9 * nat * nat, phip, 1, phi, 1)  
  endif

  !
  !    Here we symmetrize with respect to the small group of q
  !
  if (nsymq.eq.1) return  

  do na = 1, nat  
     do nb = 1, nat  
        iflb (na, nb) = 0  
     enddo
  enddo
  do na = 1, nat  
     do nb = 1, nat  
        if (iflb (na, nb) .eq.0) then  
           call setv (18, 0.d0, work, 1)  
           do isymq = 1, nsymq  
              irot = irgq (isymq)  
              sna = irt (irot, na)  
              snb = irt (irot, nb)  
              arg = 0.d0  
              do ipol = 1, 3  
                 arg = arg + (xq (ipol) * (rtau (ipol, irot, na) - rtau (ipol, &
                      irot, nb) ) )
              enddo
              arg = arg * tpi  
              faseq (isymq) = DCMPLX (cos (arg), sin (arg) )  
              do ipol = 1, 3  
                 do jpol = 1, 3  
                    do kpol = 1, 3  
                       do lpol = 1, 3  
                          work (ipol, jpol) = work (ipol, jpol) + s (ipol, kpol, irot) &
                               * s (jpol, lpol, irot) * phi (kpol, lpol, sna, snb) * faseq ( &
                               isymq)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           do isymq = 1, nsymq  
              irot = irgq (isymq)  
              sna = irt (irot, na)  
              snb = irt (irot, nb)  
              do ipol = 1, 3  
                 do jpol = 1, 3  
                    phi (ipol, jpol, sna, snb) = (0.d0, 0.d0)  
                    do kpol = 1, 3  
                       do lpol = 1, 3  
                          phi (ipol, jpol, sna, snb) = phi (ipol, jpol, sna, snb) &
                               + s (ipol, kpol, invs (irot) ) * s (jpol, lpol, invs (irot) ) &
                               * work (kpol, lpol) * conjg (faseq (isymq) )
                       enddo
                    enddo
                 enddo
              enddo
              iflb (sna, snb) = 1  
           enddo
        endif
     enddo
  enddo
  call DSCAL (18 * nat * nat, 1.d0 / nsymq, phi, 1)  
  return  
end subroutine symdynph_gq
