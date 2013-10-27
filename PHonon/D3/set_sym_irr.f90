!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE set_sym_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
     irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u,         &
     npert, nirr, gi, gimq, iverbosity)
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

  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm, mpime, root
  !
  IMPLICIT NONE
!
!   first the dummy variables
!

  INTEGER ::  nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, npert (3 * nat), irgq (48), nsymq, irotmq, nirr, npertx
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

  REAL(DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3), &
       gi (3, 48), gimq (3)
! input: the q point
! input: the R associated to each tau
! input: the direct lattice vectors
! input: the reciprocal lattice vectors
! output: [S(irotq)*q - q]
! output: [S(irotmq)*q + q]

  COMPLEX(DP) :: u (3 * nat, 3 * nat),         &
       t (npertx, npertx, 48, 3 * nat),   &
       tmq (npertx, npertx, 3 * nat)
! output: the pattern vectors
! output: the symmetry matrices
! output: the matrice sending q -> -q+G
  LOGICAL :: minus_q
! output: if true one symmetry send q -
!
  INTEGER :: na, nb, imode, jmode, ipert, jpert, nsymtot, imode0, &
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

  REAL(DP) :: eigen (3 * nat), modul, arg
! the eigenvalues of dynamical ma
! the modulus of the mode
! the argument of the phase

  COMPLEX(DP) :: wdyn (3, 3, nat, nat), phi (3 * nat, 3 * nat), &
       wrk_u (3, nat), wrk_ru (3, nat), fase
! the dynamical matrix
! the bi-dimensional dynamical ma
! one pattern
! the rotated of one pattern
! the phase factor

  LOGICAL :: lgamma
! if true gamma point

  IF ( mpime == root ) THEN
     !
     !   Allocate the necessary quantities
     !
     lgamma = (xq(1).EQ.0.d0 .AND. xq(2).EQ.0.d0 .AND. xq(3).EQ.0.d0)
     !
     !   find the small group of q
     !
     CALL smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
     !
     !   And we compute the matrices which represent the symmetry transformat
     !   in the basis of the displacements
     !
     t(:,:,:,:) = (0.d0, 0.d0)
     tmq(:,:,:) = (0.d0, 0.d0)
     IF (minus_q) THEN
        nsymtot = nsymq + 1
     ELSE
        nsymtot = nsymq

     ENDIF
     DO isymq = 1, nsymtot
        IF (isymq.LE.nsymq) THEN
           irot = irgq (isymq)
        ELSE
           irot = irotmq
        ENDIF
        imode0 = 0
        DO irr = 1, nirr
           DO ipert = 1, npert (irr)
              imode = imode0 + ipert
              DO na = 1, nat
                 DO ipol = 1, 3
                    jmode = 3 * (na - 1) + ipol
                    wrk_u (ipol, na) = u (jmode, imode)
                 ENDDO
              ENDDO
              !
              !     transform this pattern to crystal basis
              !
              DO na = 1, nat
                 CALL trnvecc (wrk_u (1, na), at, bg, - 1)
              ENDDO
              !
              !     the patterns are rotated with this symmetry
              !
              wrk_ru(:,:) = (0.d0, 0.d0)
              DO na = 1, nat
                 sna = irt (irot, na)
                 arg = 0.d0
                 DO ipol = 1, 3
                    arg = arg + xq (ipol) * rtau (ipol, irot, na)
                 ENDDO
                 arg = arg * tpi
                 IF (isymq.EQ.nsymtot.AND.minus_q) THEN
                    fase = CMPLX(COS (arg), SIN (arg) ,kind=DP)
                 ELSE
                    fase = CMPLX(COS (arg), - SIN (arg) ,kind=DP)
                 ENDIF
                 DO ipol = 1, 3
                    DO jpol = 1, 3
                       wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + s (jpol, ipol, irot) &
                            * wrk_u (jpol, na) * fase
                    ENDDO
                 ENDDO
              ENDDO
              !
              !    Transform back the rotated pattern
              !
              DO na = 1, nat
                 CALL trnvecc (wrk_ru (1, na), at, bg, 1)
              ENDDO
              !
              !     Computes the symmetry matrices on the basis of the pattern
              !
              DO jpert = 1, npert (irr)
                 imode = imode0 + jpert
                 DO na = 1, nat
                    DO ipol = 1, 3
                       jmode = ipol + (na - 1) * 3
                       IF (isymq.EQ.nsymtot.AND.minus_q) THEN
                          tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + CONJG (u ( &
                               jmode, imode) * wrk_ru (ipol, na) )
                       ELSE
                          t (jpert, ipert, irot, irr) = t (jpert, ipert, irot, irr) &
                               + CONJG (u (jmode, imode) ) * wrk_ru (ipol, na)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           imode0 = imode0 + npert (irr)
        ENDDO

     ENDDO
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
     !      WRITE( stdout,'(" Setting mode for testing ")')
     !      do na=1,3*nat
     !         WRITE( stdout,*) u(na,1)
     !      enddo
     !      nsymq=1
     !      minus_q=.false.

     !
     ! parallel stuff: first node broadcasts everything to all nodes
     !
  END IF

  CALL mp_bcast (gi, root, world_comm)
  CALL mp_bcast (gimq, root, world_comm)
  CALL mp_bcast (t, root, world_comm)
  CALL mp_bcast (tmq, root, world_comm)
  CALL mp_bcast (u, root, world_comm)
  CALL mp_bcast (nsymq, root, world_comm)
  CALL mp_bcast (npert, root, world_comm)
  CALL mp_bcast (nirr, root, world_comm)
  CALL mp_bcast (irotmq, root, world_comm)
  CALL mp_bcast (irgq, root, world_comm)
  CALL mp_bcast (minus_q, root, world_comm)

  RETURN
END SUBROUTINE set_sym_irr
