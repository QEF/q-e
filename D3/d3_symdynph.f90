!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine d3_symdynph (xq, phi, s, invs, rtau, irt, irgq, nsymq, &
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
  real (8) :: xq (3), rtau (3, 48, nat)
  ! input: the q point
  ! input: the R associated at each t

  logical :: minus_q
  ! input: true if a symmetry q->-q+G
  complex (8) :: phi (3, 3, 3, nat, nat, nat)
  ! inp/out: the matrix to symmetr
  !
  !     One parameter
  !
  real (8) :: tpi
  parameter (tpi = 2.0d0 * 3.14159265358979d0)
  !
  !   and the local variables
  !

  integer :: isymq, sna, snb, snc, irot, na, nb, nc, ipol, jpol, &
       lpol, kpol, mpol, npol
  ! counter on symmetries
  ! the rotated of the a atom
  ! the rotated of the b atom
  ! the rotated of the b atom
  ! counter on rotations
  ! counter on atoms
  ! counter on atoms
  ! counter on atoms
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  integer, allocatable:: iflb (:,:,:)
  ! used to account for symmetrized elements

  real (8) :: arg
  ! the argument of the phase

  complex (8), allocatable :: phip (:,:,:,:,:,:)
  complex (8) :: work (3, 3, 3), fase, faseq (48)
  ! working space
  ! the phase factor
  ! the phases for each symmetry

  !
  !    We start by imposing hermiticity
  !
  do nc = 1, nat
     do na = 1, nat
        do nb = 1, nat
           do kpol = 1, 3
              do ipol = 1, 3
                 do jpol = 1, 3
                    phi (kpol, ipol, jpol, nc, na, nb) = 0.5d0 * (phi (kpol, ipol, &
                         jpol, nc, na, nb) + conjg (phi (kpol, jpol, ipol, nc, nb, na) ) )
                    phi (kpol, jpol, ipol, nc, nb, na) = conjg (phi (kpol, ipol, jpol, &
                         nc, na, nb) )
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  !    If no other symmetry is present we quit here
  !

  if ( (nsymq.eq.1) .and. (.not.minus_q) ) return
  allocate  (phip( 3, 3, 3, nat, nat, nat))    
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
  if (minus_q) then
     do nc = 1, nat
        do na = 1, nat
           do nb = 1, nat
              do mpol = 1, 3
                 do ipol = 1, 3
                    do jpol = 1, 3
                       call setv (54, 0.d0, work, 1)
                       snc = irt (irotmq, nc)
                       sna = irt (irotmq, na)
                       snb = irt (irotmq, nb)
                       arg = 0.d0
                       do kpol = 1, 3
                          arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) - rtau (kpol, &
                               irotmq, nb) ) )
                       enddo
                       arg = arg * tpi
                       fase = DCMPLX (cos (arg), sin (arg) )
                       do npol = 1, 3
                          do kpol = 1, 3
                             do lpol = 1, 3
                                work (mpol, ipol, jpol) = work (mpol, ipol, jpol) + fase * s ( &
                                     ipol, kpol, irotmq) * s (jpol, lpol, irotmq) * s (mpol, npol, &
                                     irotmq) * phi (npol, kpol, lpol, snc, sna, snb)
                             enddo
                          enddo
                       enddo
                       phip (mpol, ipol, jpol, nc, na, nb) = (phi (mpol, ipol, jpol, &
                            nc, na, nb) + conjg (work (mpol, ipol, jpol) ) ) * 0.5d0
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     call ZCOPY (27 * nat * nat * nat, phip, 1, phi, 1)
  endif

  deallocate (phip)
  !
  !    Here we symmetrize with respect to the small group of q
  !
  if (nsymq.eq.1) return

  allocate  (iflb( nat, nat, nat))    
  do na = 1, nat
     do nb = 1, nat
        do nc = 1, nat
           iflb (nc, na, nb) = 0
        enddo
     enddo

  enddo
  do nc = 1, nat
     do na = 1, nat
        do nb = 1, nat
           if (iflb (nc, na, nb) .eq.0) then
              call setv (54, 0.d0, work, 1)
              do isymq = 1, nsymq
                 irot = irgq (isymq)
                 snc = irt (irot, nc)
                 sna = irt (irot, na)
                 snb = irt (irot, nb)
                 arg = 0.d0
                 do ipol = 1, 3
                    arg = arg + (xq (ipol) * (rtau (ipol, irot, na) - rtau (ipol, &
                         irot, nb) ) )
                 enddo
                 arg = arg * tpi
                 faseq (isymq) = DCMPLX (cos (arg), sin (arg) )
                 do mpol = 1, 3
                    do ipol = 1, 3
                       do jpol = 1, 3
                          do npol = 1, 3
                             do kpol = 1, 3
                                do lpol = 1, 3
                                   work (mpol, ipol, jpol) = work (mpol, ipol, jpol) + s (ipol, &
                                        kpol, irot) * s (jpol, lpol, irot) * s (mpol, npol, irot) &
                                        * phi (npol, kpol, lpol, snc, sna, snb) * faseq (isymq)
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              do isymq = 1, nsymq
                 irot = irgq (isymq)
                 snc = irt (irot, nc)
                 sna = irt (irot, na)
                 snb = irt (irot, nb)
                 do mpol = 1, 3
                    do ipol = 1, 3
                       do jpol = 1, 3
                          phi (mpol, ipol, jpol, snc, sna, snb) = (0.d0, 0.d0)
                          do npol = 1, 3
                             do kpol = 1, 3
                                do lpol = 1, 3
                                   phi (mpol, ipol, jpol, snc, sna, snb) = phi (mpol, ipol, jpol, &
                                        snc, sna, snb) + s (mpol, npol, invs (irot) ) * s (ipol, kpol, &
                                        invs (irot) ) * s (jpol, lpol, invs (irot) ) * work (npol, &
                                        kpol, lpol) * conjg (faseq (isymq) )
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
                 iflb (snc, sna, snb) = 1
              enddo
           endif
        enddo
     enddo
  enddo
  call DSCAL (54 * nat * nat * nat, 1.d0 / nsymq, phi, 1)
  deallocate (iflb)
  return
end subroutine d3_symdynph
