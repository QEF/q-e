!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
     irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
     nirr, gi, gimq, iverbosity)
!---------------------------------------------------------------------
!
!     This subroutine computes a basis for all the irreducible
!     representations of the small group of q, which are contained
!     in the representation which has as basis the displacement vectors.
!     This is achieved by building a random hermitean matrix,
!     symmetrizing it and diagonalizing the result. The eigenvectors
!     give a basis for the irreducible representations of the
!     small group of q.
!
!     Furthermore it computes:
!     1) the small group of q
!     2) the possible G vectors associated to every symmetry operation
!     3) the matrices which represent the small group of q on the
!        pattern basis.
!
!     Original routine was from C. Bungaro.
!     Revised Oct. 1995 by Andrea Dal Corso.
!     April 1997: parallel stuff added (SdG)
!
#include "machine.h"
  use parameters, only : DP
#ifdef __PARA
  use mp, only: mp_bcast
#endif
  implicit none
!
!   first the dummy variables
!

  integer ::  nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, npert (3 * nat), irgq (48), nsymq, irotmq, nirr, max_irr_dim
! input: the number of atoms
! input: the number of symmetries
! input: the symmetry matrices
! input: the inverse of each matrix
! input: the rotated of each atom
! input: write control
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

  complex(kind=DP) :: u (3 * nat, 3 * nat), t (max_irr_dim, max_irr_dim, 48, 3 * nat), &
       tmq (max_irr_dim, max_irr_dim, 3 * nat)
! output: the pattern vectors
! output: the symmetry matrices
! output: the matrice sending q -> -q+G
  logical :: minus_q
! output: if true one symmetry send q -
!
!   here the local variables
!
  real(kind=DP), parameter :: tpi = 2.0d0 * 3.14159265358979d0

  integer :: na, nb, imode, jmode, ipert, jpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna
  ! counter on atoms
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

  integer :: info

  real(kind=DP) :: eigen (3 * nat), modul, arg
! the eigenvalues of dynamical ma
! the modulus of the mode
! the argument of the phase

  complex(kind=DP) :: wdyn (3, 3, nat, nat), phi (3 * nat, 3 * nat), &
       wrk_u (3, nat), wrk_ru (3, nat), fase
! the dynamical matrix
! the bi-dimensional dynamical ma
! one pattern
! the rotated of one pattern
! the phase factor

  logical :: lgamma
! if true gamma point
!
!   Allocate the necessary quantities
!
  lgamma = (xq(1).eq.0.d0 .and. xq(2).eq.0.d0 .and. xq(3).eq.0.d0)
!
!   find the small group of q
!
  call smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
!
!   then we generate a random hermitean matrix
!
  call set_rndm_seed(1)
  call random_matrix (irt,irgq,nsymq,minus_q,irotmq,nat,wdyn,lgamma)
!call write_matrix('random matrix',wdyn,nat)
!
! symmetrize the random matrix with the little group of q
!
  call symdynph_gq (xq,wdyn,s,invs,rtau,irt,irgq,nsymq,nat,irotmq,minus_q)
!call write_matrix('symmetrized matrix',wdyn,nat)
!
!  Diagonalize the symmetrized random matrix.
!  Transform the symmetryzed matrix, currently in crystal coordinates,
!  in cartesian coordinates.
!
  do na = 1, nat
     do nb = 1, nat
        call trntnsc( wdyn(1,1,na,nb), at, bg, 1 )
     enddo
  enddo
!
!     We copy the dynamical matrix in a bidimensional array
!
  do na = 1, nat
     do nb = 1, nat
        do ipol = 1, 3
           imode = ipol + 3 * (na - 1)
           do jpol = 1, 3
              jmode = jpol + 3 * (nb - 1)
              phi (imode, jmode) = wdyn (ipol, jpol, na, nb)

           enddo
        enddo
     enddo
  enddo
!
!   Diagonalize
!
  call cdiagh (3 * nat, phi, 3 * nat, eigen, u)
#ifdef __PARA
!
!  Make sure all nodes have the same patterns
!
  call check (3 * nat, eigen)
  call check (18 * nat * nat, u)
#endif
!
!   We adjust the phase of each mode in such a way that the first
!   non zero element is real
!
  do imode = 1, 3 * nat
     do na = 1, 3 * nat
        modul = abs (u(na, imode) )
        if (modul.gt.1d-9) then
           fase = u (na, imode) / modul
           goto 110
        endif
     enddo
     call errore ('set_irr', 'one mode is zero', imode)
110  do na = 1, 3 * nat
        u (na, imode) = - u (na, imode) * conjg (fase)
     enddo
  enddo
!
!  We have here a test which writes eigenvectors and eigenvalues
!
      if (iverbosity.eq.1) then
         do imode=1,3*nat
            write(6, '(2x,"autoval = ", e10.4)') eigen(imode)
            write(6, '(2x,"Real(aut_vet)= ( ",6f10.5,")")') &
                ( DREAL(u(na,imode)), na=1,3*nat )
            write(6, '(2x,"Imm(aut_vet)= ( ",6f10.5,")")') &
                ( DIMAG(u(na,imode)), na=1,3*nat )
         end do
      end if
!
!  Here we count the irreducible representations and their dimensions
!
  do imode = 1, 3 * nat
! initialization
     npert (imode) = 0
  enddo
  nirr = 1
  npert (1) = 1
  do imode = 2, 3 * nat
     if (abs (eigen (imode) - eigen (imode-1) ) / (abs (eigen (imode) ) &
          + abs (eigen (imode-1) ) ) .lt.1.d-4) then
        npert (nirr) = npert (nirr) + 1
        if (npert (nirr) .gt. max_irr_dim) call errore &
                         ('set_irr', 'npert > max_irr_dim ', nirr)
     else
        nirr = nirr + 1
        npert (nirr) = 1
     endif

  enddo
!
!   And we compute the matrices which represent the symmetry transformat
!   in the basis of the displacements
!
  t(:,:,:,:) = (0.d0, 0.d0)
  tmq(:,:,:) = (0.d0, 0.d0)
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
!
!    Note: the following lines are for testing purposes
!
!      nirr = 1
!      npert(1)=1
!      do na=1,3*nat/2
!        u(na,1)=(0.d0,0.d0)
!        u(na+3*nat/2,1)=(0.d0,0.d0)
!      enddo
!      u(1,1)=(-1.d0,0.d0)
!      write(6,'(" Setting mode for testing ")')
!      do na=1,3*nat
!         write(6,*) u(na,1)
!      enddo
!      nsymq=1
!      minus_q=.false.

#ifdef __PARA
!
! parallel stuff: first node broadcasts everything to all nodes
!
400 continue
  call mp_bcast (gi, 0)
  call mp_bcast (gimq, 0)
  call mp_bcast (t, 0)
  call mp_bcast (tmq, 0)
  call mp_bcast (u, 0)
  call mp_bcast (nsymq, 0)
  call mp_bcast (npert, 0)
  call mp_bcast (nirr, 0)
  call mp_bcast (irotmq, 0)
  call mp_bcast (irgq, 0)
  call mp_bcast (minus_q, 0)
#endif
  return
end subroutine set_irr
