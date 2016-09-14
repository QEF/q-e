!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  This file contains a series of obsolete routines that are here 
!  for compatibility. 
!
!
!-----------------------------------------------------------------------

subroutine smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, &
     minus_q, gi, gimq)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave q unchanged.
  ! Furthermore it checks if one of the matrices send q <-> -q+G. In
  ! this case minus_q is set true.
  !
  ! Revised   2 Sept. 1995 by Andrea Dal Corso
  ! Modified 22 April 1997 by SdG: minus_q is sought also among sym.op.
  !                such that Sq=q+G (i.e. the case q=-q+G is dealt with).
  !
  !
  !  The dummy variables
  !
  USE kinds, only : DP
  implicit none

  real(DP), parameter :: accep=1.e-5_dp
  real(DP) :: bg (3, 3), at (3, 3), xq (3), gi (3, 48), gimq (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal
  ! output: the G associated to a symmetry:[S(irotq)*q - q]
  ! output: the G associated to:  [S(irotmq)*q + q]

  integer :: s (3, 3, 48), irgq (48), irotmq, nsymq, nsym
  ! input: the symmetry matrices
  ! output: the symmetry of the small group
  ! output: op. symmetry: s_irotmq(q)=-q+G
  ! output: dimension of the small group of q
  ! input: dimension of the point group

  logical :: minus_q
  ! input: .t. if sym.ops. such that Sq=-q+G are searched for
  ! output: .t. if such a symmetry has been found

  real(DP) :: wrk (3), aq (3), raq (3), zero (3)
  ! additional space to compute gi and gimq
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: isym, ipol, jpol
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations

  logical :: look_for_minus_q, eqvect
  ! .t. if sym.ops. such that Sq=-q+G are searched for
  ! logical function, check if two vectors are equal
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  look_for_minus_q = minus_q
  !
  minus_q = .false.
  zero = 0.d0
  gi   = 0.d0
  gimq = 0.d0
  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  !
  !   test all symmetries to see if the operation S sends q in q+G ...
  !
  nsymq = 0
  do isym = 1, nsym
     raq = 0.d0
     do ipol = 1, 3
        do jpol = 1, 3
           raq (ipol) = raq (ipol) + DBLE (s (ipol, jpol, isym) ) * &
                aq (jpol)
        enddo
     enddo
     if (eqvect (raq, aq, zero, accep) ) then
        nsymq = nsymq + 1
        irgq (nsymq) = isym
        do ipol = 1, 3
           wrk (ipol) = raq (ipol) - aq (ipol)
        enddo
        call cryst_to_cart (1, wrk, bg, 1)
        gi (:, nsymq) = wrk (:)
        !
        !   ... and in -q+G
        !
        if (look_for_minus_q.and..not.minus_q) then
           raq (:) = - raq(:)
           if (eqvect (raq, aq, zero, accep) ) then
              minus_q = .true.
              irotmq = isym
              do ipol = 1, 3
                 wrk (ipol) = - raq (ipol) + aq (ipol)
              enddo
              call cryst_to_cart (1, wrk, bg, 1)
              gimq (:) = wrk (:)
           endif
        endif
     endif
  enddo
  !
  ! if xq=(0,0,0) minus_q always apply with the identity operation
  !
  if (xq (1) == 0.d0 .and. xq (2) == 0.d0 .and. xq (3) == 0.d0) then
     minus_q = .true.
     irotmq = 1
     gimq = 0.d0
  endif
  !
  return
end subroutine smallgq

!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine set_irr (nat, at, bg, xq, s, sr, tau, ntyp, ityp, ftau, invs, nsym, &
                    rtau, irt, irgq, nsymq, minus_q, irotmq, u, npert,   &
                    nirr, gi, gimq, iverbosity, u_from_file, eigen, search_sym,&
                    nspin_mag, t_rev, amass, num_rap_mode, name_rap_mode)
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
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  USE constants, ONLY: tpi
  USE random_numbers, ONLY : randy
  USE rap_point_group, ONLY : name_rap
#if defined(__MPI)
  use mp, only: mp_bcast
  use io_global, only : ionode_id
  use mp_global, only : intra_image_comm
#endif
  implicit none
!
!   first the dummy variables
!
  integer ::  nat, ntyp, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, npert (3 * nat), irgq (48), nsymq, irotmq, nirr, &
       ftau(3,48), nspin_mag, t_rev(48), ityp(nat), num_rap_mode(3*nat)
! input: the number of atoms
! input: the number of types of atoms
! input: the number of symmetries
! input: the symmetry matrices
! input: the inverse of each matrix
! input: the rotated of each atom
! input: write control
! output: the dimension of each representation
! output: the small group of q
! output: the order of the small group
! output: the symmetry sending q -> -q+
! output: the number of irr. representation
! input: the fractionary translations
! input: the number of spin components
! input: the time reversal symmetry
! input: the type of each atom
! output: the number of the representation of each mode

  real(DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3), &
       gi (3, 48), gimq (3), sr(3,3,48), tau(3,nat), amass(ntyp)
! input: the q point
! input: the R associated to each tau
! input: the direct lattice vectors
! input: the reciprocal lattice vectors
! output: [S(irotq)*q - q]
! output: [S(irotmq)*q + q]
! input: symmetry matrices in cartesian coordinates
! input: the atomic positions
! input: the mass of each atom (in amu)

  complex(DP) :: u(3*nat, 3*nat)
! output: the pattern vectors
  logical :: minus_q, u_from_file, search_sym
! output: if true one symmetry send q -
! input: if true the displacement patterns are not calculated here
! output: if true the symmetry of each mode has been calculated

  character(len=15) :: name_rap_mode( 3 * nat )
  ! output: the name of the representation for each group of modes
!
!   here the local variables
!
  integer :: na, nb, imode, jmode, ipert, jpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna, isym
  ! counters and auxiliary variables

  integer :: info, mode_per_rap(12), count_rap(12), rap, init, pos, irap, &
             num_rap_aux( 3 * nat )

  real(DP) :: eigen (3 * nat), modul, arg, eig(3*nat)
! the eigenvalues of dynamical matrix
! the modulus of the mode
! the argument of the phase

  complex(DP) :: wdyn (3, 3, nat, nat), phi (3 * nat, 3 * nat), &
       wrk_u (3, nat), wrk_ru (3, nat), fase
! the dynamical matrix
! the dynamical matrix with two indices
! pattern
! rotated pattern
! the phase factor

  logical :: lgamma, magnetic_sym
! if true gamma point
!
!   Allocate the necessary quantities
!
  lgamma = (xq(1) == 0.d0 .and. xq(2) == 0.d0 .and. xq(3) == 0.d0)
!
!   find the small group of q
!
  call smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
  ! are there non-symmorphic operations?
  ! note that in input search_sym should be initialized to=.true.
  IF ( ANY ( ftau(:,1:nsymq) /= 0 ) ) THEN
     DO isym=1,nsymq
        search_sym=( search_sym.and.(abs(gi(1,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(2,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(3,irgq(isym)))<1.d-8) )
     END DO
  END IF
  num_rap_mode=-1
  IF (search_sym) THEN
     magnetic_sym=(nspin_mag==4)
     CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
  ENDIF

  IF (.NOT. u_from_file) THEN
!
!   then we generate a random hermitean matrix
!
     arg = randy(0)
     call random_matrix (irt,irgq,nsymq,minus_q,irotmq,nat,wdyn,lgamma)
!call write_matrix('random matrix',wdyn,nat)
!
! symmetrize the random matrix with the little group of q
!
     call symdynph_gq (xq,wdyn,s,invs,rtau,irt,irgq,nsymq,nat,irotmq,minus_q)
!call write_matrix('symmetrized matrix',wdyn,nat)
!
!  Diagonalize the symmetrized random matrix.
!  Transform the symmetrized matrix, currently in crystal coordinates,
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
110     do na = 1, 3 * nat
           u (na, imode) = - u (na, imode) * CONJG(fase)
        enddo
     enddo
!
!  We have here a test which writes eigenvectors and eigenvalues
!
     if (iverbosity.eq.1) then
        do imode=1,3*nat
           WRITE( stdout, '(2x,"autoval = ", e10.4)') eigen(imode)
           WRITE( stdout, '(2x,"Real(aut_vet)= ( ",6f10.5,")")') &
               (  DBLE(u(na,imode)), na=1,3*nat )
           WRITE( stdout, '(2x,"Imm(aut_vet)= ( ",6f10.5,")")') &
               ( AIMAG(u(na,imode)), na=1,3*nat )
        end do
     end if

     IF (search_sym) THEN
        CALL find_mode_sym (u, eigen, at, bg, tau, nat, nsymq, &
                  sr, irt, xq, rtau, amass, ntyp, ityp, 0, lgamma, &
                  .FALSE., nspin_mag, name_rap_mode, num_rap_mode )
!
!   Order the modes so that we first make all those that belong to the first
!   representation, then the second ect.
!
!
!   First count, for each irreducible representation, how many modes
!   belong to that representation
!
        mode_per_rap=0
        DO imode=1,3*nat
           mode_per_rap(num_rap_mode(imode))=mode_per_rap(num_rap_mode(imode))+1
        ENDDO
!
!   The position of each mode on the list is the following:
!   The positions from 1 to nrap(1) contain the modes that transform according
!   to the first representation. From nrap(1)+1 to nrap(1)+nrap(2) the
!   mode that transform according to the second ecc.
!
        count_rap=1
        DO imode=1,3*nat
           rap=num_rap_mode(imode)
           IF (rap>12) call errore('set_irr',&
                                   'problem with the representation',1)
           init=0
           DO irap=1,rap-1
              init=init+mode_per_rap(irap)
           ENDDO
           pos=init+count_rap(rap)
           eig(pos)=eigen(imode)
           phi(:,pos)=u(:,imode)
           num_rap_aux(pos)=num_rap_mode(imode)
           count_rap(rap)=count_rap(rap)+1
        ENDDO
        eigen=eig
        u=phi
        num_rap_mode=num_rap_aux

! Modes with accidentally degenerate eigenvalues, or with eigenvalues 
! degenerate due to time reversal must be calculated together even if
! they belong to different irreducible representations. 
!
        DO imode=1,3*nat-1
           DO jmode = imode+1, 3*nat
              IF ((num_rap_mode(imode) /= num_rap_mode(jmode)).AND.  &
                  (ABS(eigen(imode) - eigen(jmode))/   &
                  (ABS(eigen(imode)) + ABS (eigen (jmode) )) < 1.d-4) ) THEN
                 eig(1)=eigen(jmode)
                 phi(:,1)=u(:,jmode)
                 num_rap_aux(1)=num_rap_mode(jmode)
                 eigen(jmode)=eigen(imode+1)
                 u(:,jmode)=u(:,imode+1)
                 num_rap_mode(jmode)=num_rap_mode(imode+1)
                 eigen(imode+1)=eig(1)
                 u(:,imode+1)=phi(:,1)
                 num_rap_mode(imode+1)=num_rap_aux(1)
              ENDIF
           ENDDO
        ENDDO
     ENDIF
!
!  Here we count the irreducible representations and their dimensions
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
        else
           nirr = nirr + 1
           npert (nirr) = 1
        endif
     enddo
     IF (search_sym) THEN
        imode=1
        DO irr=1,nirr
           name_rap_mode(irr)=name_rap(num_rap_mode(imode))
           imode=imode+npert(irr)
        ENDDO
     ENDIF
  endif
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

#if defined(__MPI)
!
! parallel stuff: first node broadcasts everything to all nodes
!
400 continue
  call mp_bcast (gi, ionode_id, intra_image_comm)
  call mp_bcast (gimq, ionode_id, intra_image_comm)
  call mp_bcast (u, ionode_id, intra_image_comm)
  call mp_bcast (nsymq, ionode_id, intra_image_comm)
  call mp_bcast (npert, ionode_id, intra_image_comm)
  call mp_bcast (nirr, ionode_id, intra_image_comm)
  call mp_bcast (irotmq, ionode_id, intra_image_comm)
  call mp_bcast (irgq, ionode_id, intra_image_comm)
  call mp_bcast (minus_q, ionode_id, intra_image_comm)
  call mp_bcast (num_rap_mode, ionode_id, intra_image_comm)
  call mp_bcast (name_rap_mode, ionode_id, intra_image_comm)
#endif
  return
end subroutine set_irr


!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, &
     irt, irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, &
     npert, nirr, gi, gimq, iverbosity)
  !---------------------------------------------------------------------
  !
  !     This routine substitute set_irr when there are no symmetries.
  !     The irreducible representations are all one dimensional and
  !     we set them to the displacement of a single atom in one direction
  !
  USE kinds, only : DP
  implicit none
  !
  !   first the dummy variables
  !

  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
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
  integer :: imode
  ! counter on modes
  !
  !    set the information on the symmetry group
  !
  call smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
  !
  !     set the modes
  !
  u (:,:) = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     u (imode, imode) = (1.d0, 0.d0)
  enddo
  nirr = 3 * nat
  do imode = 1, 3 * nat
     npert (imode) = 1
  enddo
  !
  !   And we compute the matrices which represent the symmetry transformat
  !   in the basis of the displacements
  !
  t(:, :, :, :) = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     t (1, 1, 1, imode) = (1.d0, 0.d0)
  enddo

  tmq (:, :, :) = (0.d0, 0.d0)
  if (minus_q) then
     tmq (1, 1, :) = (1.d0, 0.d0)
  end if

  return
end subroutine set_irr_nosym

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


!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine set_irr_sym (nat, at, bg, xq, s, rtau, irt, &
     irgq, nsymq, minus_q, irotmq, t, tmq, u, npert, nirr, npertx )
!---------------------------------------------------------------------
!
!     This subroutine computes: 
!     1) the matrices which represent the small group of q on the
!        pattern basis.
!
  USE kinds, ONLY : DP
  USE constants, ONLY: tpi

  USE mp, ONLY: mp_bcast
  USE mp_global, ONLY : intra_image_comm
  USE io_global, ONLY : ionode_id
  implicit none
!
!   first the dummy variables
!

  integer, intent(in) ::  nat, s (3, 3, 48), irt (48, nat), npert (3 * nat), &
                          irgq (48), nsymq, irotmq, nirr, npertx
! input: the number of atoms
! input: the symmetry matrices
! input: the rotated of each atom
! input: the dimension of each represe
! input: the small group of q
! input: the order of the small group
! input: the symmetry sending q -> -q+
! input: the number of irr. representa

  real(DP), intent(in) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3)
! input: the q point
! input: the R associated to each tau
! input: the direct lattice vectors
! input: the reciprocal lattice vectors

  complex(DP), intent(in) :: u(3*nat, 3*nat)
! input: the pattern vectors

  complex(DP), intent(out) :: t(npertx, npertx, 48, 3*nat), tmq (npertx, npertx, 3*nat)
! output: the symmetry matrices
! output: the matrice sending q -> -q+G
  logical :: minus_q
! output: if true one symmetry send q -
!
!   here the local variables
!
  integer :: na, imode, jmode, ipert, jpert, kpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna
  ! counters and auxiliary variables

  real(DP) :: arg
! the argument of the phase

  complex(DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase, wrk
! pattern
! rotated pattern
! the phase factor

!
!   We compute the matrices which represent the symmetry transformation
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
                 fase = CMPLX (cos (arg), sin (arg) )
              else
                 fase = CMPLX (cos (arg), - sin (arg) )
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
                       tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + CONJG(u ( &
                            jmode, imode) * wrk_ru (ipol, na) )
                    else
                       t (jpert, ipert, irot, irr) = t (jpert, ipert, irot, irr) &
                            + CONJG(u (jmode, imode) ) * wrk_ru (ipol, na)
                    endif
                 enddo
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)

!
! If the representations are irreducible, the rotations should be unitary matrices
! if this is not the case, the way the representations have been chosen has failed
! for some reasons (check set_irr.f90)
!
        do ipert = 1, npert (irr)
           do jpert = 1, npert (irr)
              wrk = cmplx(0.d0,0.d0)
              do kpert = 1, npert (irr)
                 wrk = wrk + t (ipert,kpert,irot,irr) * conjg( t(jpert,kpert,irot,irr))
              enddo
              if (jpert.ne.ipert .and. abs(wrk).gt. 1.d-6 ) &
                     call errore('set_irr_sym','wrong representation',100*irr+10*jpert+ipert)
              if (jpert.eq.ipert .and. abs(wrk-1.d0).gt. 1.d-6 ) &
                     call errore('set_irr_sym','wrong representation',100*irr+10*jpert+ipert)
           enddo
        enddo
     enddo

  enddo

#if defined(__MPI)
!
! parallel stuff: first node broadcasts everything to all nodes
!
  call mp_bcast (t, ionode_id, intra_image_comm)
  call mp_bcast (tmq, ionode_id, intra_image_comm)
#endif
  return
end subroutine set_irr_sym

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine dynmat0
  !-----------------------------------------------------------------------
  !
  !     This routine computes the part of the dynamical matrix which
  !     does not depend upon the change of the Bloch wavefunctions.
  !     It is a driver which calls the routines dynmat_## and d2ionq
  !     for computing respectively the electronic part and
  !     the ionic part
  !
  !
  !
  USE ions_base, ONLY : nat,ntyp => nsp, ityp, zv, tau
  USE cell_base, ONLY: alat, omega, at, bg
  USE gvect, ONLY: g, gg, ngm, gcutm
  USE symm_base, ONLY: irt, s, invs
  USE control_flags, ONLY : modenum
  USE kinds,         ONLY : DP
  USE ph_restart,    ONLY : ph_writefile
  USE control_ph,    ONLY : rec_code_read, current_iq
  USE modes,         ONLY : u, nmodes
  USE partial,       ONLY : done_irr, comp_irr
  USE dynmat,        ONLY : dyn, dyn00, dyn_rec

  USE lr_symm_base, ONLY : minus_q, irotmq, irgq, rtau, nsymq
  USE qpoint,       ONLY : xq

  implicit none

  integer :: nu_i, nu_j, na_icart, nb_jcart, ierr
  ! counters

  complex(DP) :: wrk, dynwrk (3 * nat, 3 * nat)
  ! auxiliary space

  IF ( .NOT.comp_irr(0) .or. done_irr(0) ) RETURN
  IF (rec_code_read > -30 ) RETURN

  call start_clock ('dynmat0')
  call zcopy (9 * nat * nat, dyn00, 1, dyn, 1)
  !
  ! first electronic contribution arising from the term  <psi|d2v|psi>
  !
  call dynmat_us()
  !
  !   Here the ionic contribution
  !
  call d2ionq (nat, ntyp, ityp, zv, tau, alat, omega, xq, at, bg, g, &
       gg, ngm, gcutm, nmodes, u, dyn)
  !
  !   Add non-linear core-correction (NLCC) contribution (if any)
  !
  call dynmatcc()
  !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q and of
  !   mode. This is done here, because this part of the dynmical matrix is
  !   saved with recover and in the other runs the symmetry group might change
  !
  if (modenum .ne. 0) then

     call symdyn_munu (dyn, u, xq, s, invs, rtau, irt, irgq, at, bg, &
          nsymq, nat, irotmq, minus_q)
     !
     ! rotate again in the pattern basis
     !
     call zcopy (9 * nat * nat, dyn, 1, dynwrk, 1)
     do nu_i = 1, 3 * nat
        do nu_j = 1, 3 * nat
           wrk = (0.d0, 0.d0)
           do nb_jcart = 1, 3 * nat
              do na_icart = 1, 3 * nat
                 wrk = wrk + CONJG(u (na_icart, nu_i) ) * &
                             dynwrk (na_icart, nb_jcart) * &
                             u (nb_jcart, nu_j)
              enddo
           enddo
           dyn (nu_i, nu_j) = wrk
        enddo
     enddo
  endif
  !      call tra_write_matrix('dynmat0 dyn',dyn,u,nat)
  dyn_rec(:,:)=dyn(:,:)
  done_irr(0) = .TRUE.
  CALL ph_writefile('data_dyn',current_iq,0,ierr)

  call stop_clock ('dynmat0')
  return
end subroutine dynmat0

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symdyn_munu (dyn, u, xq, s, invs, rtau, irt, irgq, at, &
     bg, nsymq, nat, irotmq, minus_q)
  !-----------------------------------------------------------------------
  !
  !    This routine symmetrize the dynamical matrix written in the basis
  !    of the modes
  !
  !
  USE kinds, only : DP
  implicit none
  integer :: nat, s (3, 3, 48), irt (48, nat), irgq (48), invs (48), &
       nsymq, irotmq
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each atom
  ! input: the small group of q
  ! input: the inverse of each matrix
  ! input: the order of the small gro
  ! input: the symmetry q -> -q+G

  real(DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3)
  ! input: the coordinates of q
  ! input: the R associated at each r
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors

  logical :: minus_q
  ! input: if true symmetry sends q->

  complex(DP) :: dyn (3 * nat, 3 * nat), u (3 * nat, 3 * nat)
  ! inp/out: matrix to symmetrize
  ! input: the patterns

  integer :: i, j, icart, jcart, na, nb, mu, nu
  ! counter on modes
  ! counter on modes
  ! counter on cartesian coordinates
  ! counter on cartesian coordinates
  ! counter on atoms
  ! counter on atoms
  ! counter on modes
  ! counter on modes

  complex(DP) :: work, phi (3, 3, nat, nat)
  ! auxiliary variable
  ! the dynamical matrix
  !
  ! First we transform in the cartesian coordinates
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        work = (0.d0, 0.d0)
        do mu = 1, 3 * nat
           do nu = 1, 3 * nat
              work = work + u (i, mu) * dyn (mu, nu) * CONJG(u (j, nu) )
           enddo
        enddo
        phi (icart, jcart, na, nb) = work
     enddo
  enddo
  !
  ! Then we transform to the crystal axis
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, - 1)
     enddo
  enddo
  !
  !   And we symmetrize in this basis
  !
  call symdynph_gq (xq, phi, s, invs, rtau, irt, irgq, nsymq, nat, &
       irotmq, minus_q)
  !
  !  Back to cartesian coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, + 1)
     enddo
  enddo
  !
  !  rewrite the dynamical matrix on the array dyn with dimension 3nat x 3
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        dyn (i, j) = phi (icart, jcart, na, nb)
     enddo

  enddo
  return
end subroutine symdyn_munu
!
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
  USE kinds, only : DP
  USE constants, ONLY: tpi
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
  real(DP) :: xq (3), rtau (3, 48, nat)
  ! input: the q point
  ! input: the R associated at each t

  logical :: minus_q
  ! input: true if a symmetry q->-q+G
  complex(DP) :: phi (3, 3, nat, nat)
  ! inp/out: the matrix to symmetrize
  !
  !   local variables
  !
  integer :: isymq, sna, snb, irot, na, nb, ipol, jpol, lpol, kpol, &
       iflb (nat, nat)
  ! counters, indices, work space

  real(DP) :: arg
  ! the argument of the phase

  complex(DP) :: phip (3, 3, nat, nat), work (3, 3), fase, faseq (48)
  ! work space, phase factors
  !
  !    We start by imposing hermiticity
  !
  do na = 1, nat
     do nb = 1, nat
        do ipol = 1, 3
           do jpol = 1, 3
              phi (ipol, jpol, na, nb) = 0.5d0 * (phi (ipol, jpol, na, nb) &
                   + CONJG(phi (jpol, ipol, nb, na) ) )
              phi (jpol, ipol, nb, na) = CONJG(phi (ipol, jpol, na, nb) )
           enddo
        enddo
     enddo
  enddo
  !
  !    If no other symmetry is present we quit here
  !
  if ( (nsymq == 1) .and. (.not.minus_q) ) return
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
  if (minus_q) then
     do na = 1, nat
        do nb = 1, nat
           do ipol = 1, 3
              do jpol = 1, 3
                 work(:,:) = (0.d0, 0.d0)
                 sna = irt (irotmq, na)
                 snb = irt (irotmq, nb)
                 arg = 0.d0
                 do kpol = 1, 3
                    arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) - &
                                              rtau (kpol, irotmq, nb) ) )
                 enddo
                 arg = arg * tpi
                 fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
                 do kpol = 1, 3
                    do lpol = 1, 3
                       work (ipol, jpol) = work (ipol, jpol) + &
                            s (ipol, kpol, irotmq) * s (jpol, lpol, irotmq) &
                            * phi (kpol, lpol, sna, snb) * fase
                    enddo
                 enddo
                 phip (ipol, jpol, na, nb) = (phi (ipol, jpol, na, nb) + &
                      CONJG( work (ipol, jpol) ) ) * 0.5d0
              enddo
           enddo
        enddo
     enddo
     phi = phip
  endif

  !
  !    Here we symmetrize with respect to the small group of q
  !
  if (nsymq == 1) return

  iflb (:, :) = 0
  do na = 1, nat
     do nb = 1, nat
        if (iflb (na, nb) == 0) then
           work(:,:) = (0.d0, 0.d0)
           do isymq = 1, nsymq
              irot = irgq (isymq)
              sna = irt (irot, na)
              snb = irt (irot, nb)
              arg = 0.d0
              do ipol = 1, 3
                 arg = arg + (xq (ipol) * (rtau (ipol, irot, na) - &
                                           rtau (ipol, irot, nb) ) )
              enddo
              arg = arg * tpi
              faseq (isymq) = CMPLX(cos (arg), sin (arg) ,kind=DP)
              do ipol = 1, 3
                 do jpol = 1, 3
                    do kpol = 1, 3
                       do lpol = 1, 3
                          work (ipol, jpol) = work (ipol, jpol) + &
                               s (ipol, kpol, irot) * s (jpol, lpol, irot) &
                               * phi (kpol, lpol, sna, snb) * faseq (isymq)
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
                               * work (kpol, lpol) * CONJG(faseq (isymq) )
                       enddo
                    enddo
                 enddo
              enddo
              iflb (sna, snb) = 1
           enddo
        endif
     enddo
  enddo
  phi (:, :, :, :) = phi (:, :, :, :) / DBLE(nsymq)
  return
end subroutine symdynph_gq

!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmatrix(iq_)
  !-----------------------------------------------------------------------
  !
  ! This routine is a driver which computes the symmetrized dynamical
  ! matrix at q (and in the star of q) and diagonalizes it.
  ! It writes the result on a iudyn file and writes the eigenvalues on
  ! output.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : FPI, BOHR_RADIUS_ANGS
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, amass, zv
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum
  USE cell_base,     ONLY : at, bg, celldm, ibrav, omega
  USE symm_base,     ONLY : s, sr, irt, nsym, time_reversal, invs
  USE run_info, ONLY : title
  USE dynmat,        ONLY : dyn, w2
  USE noncollin_module, ONLY : nspin_mag
  USE modes,         ONLY : u, nmodes, npert, nirr, name_rap_mode, num_rap_mode
  USE gamma_gamma,   ONLY : nasr, asr, equiv_atoms, has_equivalent, &
                            n_diff_sites
  USE efield_mod,    ONLY : epsilon, zstareu, zstarue0, zstarue
  USE control_ph,    ONLY : epsil, zue, lgamma_gamma, search_sym, ldisp, &
                            start_irr, last_irr, done_zue, where_rec, &
                            rec_code, ldiag, done_epsil, done_zeu, xmldyn, &
                            current_iq
  USE ph_restart,    ONLY : ph_writefile
  USE partial,       ONLY : all_comp, comp_irr, done_irr, nat_todo_input
  USE units_ph,      ONLY : iudyn
  USE noncollin_module, ONLY : m_loc, nspin_mag
  USE output,        ONLY : fildyn, fildrho, fildvscf
  USE io_dyn_mat,    ONLY : write_dyn_mat_header
  USE ramanm,        ONLY : lraman, ramtns
  USE dfile_star,    ONLY : write_dfile_star, drho_star, dvscf_star !write_dfile_mq
  USE units_ph,      ONLY : iudrho, iudvscf

  USE lr_symm_base, ONLY : minus_q, irotmq, nsymq, irgq, rtau
  USE qpoint,       ONLY : xq
  USE control_lr,   ONLY : lgamma

  implicit none
  INTEGER, INTENT(IN) :: iq_
  ! local variables
  !
  integer :: nq, isq (48), imq, na, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j, nqq, ierr
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real(DP) :: sxq (3, 48), work(3)
  ! list of vectors in the star of q
  real(DP), allocatable :: zstar(:,:,:)
  integer :: icart, jcart
  logical :: ldiag_loc, opnd
  !
  call start_clock('dynmatrix')
  ldiag_loc=ldiag.OR.(nat_todo_input > 0).OR.all_comp
  !
  !     set all noncomputed elements to zero
  !
  if (.not.lgamma_gamma) then
     imode0 = 0
     do irr = 1, nirr
        jmode0 = 0
        do jrr = 1, nirr
           if (.NOT.done_irr (irr) .and..NOT.done_irr (jrr) ) then
              do ipert = 1, npert (irr)
                 mu = imode0 + ipert
                 do jpert = 1, npert (jrr)
                    nu = jmode0 + jpert
                    dyn (mu, nu) = CMPLX(0.d0, 0.d0,kind=DP)
                 enddo
              enddo
           elseif (.NOT.done_irr (irr) .and.done_irr (jrr) ) then
              do ipert = 1, npert (irr)
                 mu = imode0 + ipert
                 do jpert = 1, npert (jrr)
                    nu = jmode0 + jpert
                    dyn (mu, nu) = CONJG(dyn (nu, mu) )
                 enddo
              enddo
           endif
           jmode0 = jmode0 + npert (jrr)
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  else
     do irr = 1, nirr
        if (.NOT.comp_irr(irr)) then
           do nu=1,3*nat
              dyn(irr,nu)=(0.d0,0.d0)
           enddo
        endif
     enddo
  endif

  !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q
  !
  IF (lgamma_gamma) THEN
     CALL generate_dynamical_matrix (nat, nsym, s, invs, irt, at, bg, &
                       n_diff_sites, equiv_atoms, has_equivalent, dyn)
     IF (asr) CALL set_asr_c(nat,nasr,dyn)
  ELSE
     CALL symdyn_munu (dyn, u, xq, s, invs, rtau, irt, irgq, at, bg, &
          nsymq, nat, irotmq, minus_q)
  ENDIF
  !
  !  if only one mode is computed write the dynamical matrix and stop
  !
  if (modenum .ne. 0) then
     WRITE( stdout, '(/,5x,"Dynamical matrix:")')
     do nu = 1, 3 * nat
        WRITE( stdout, '(5x,2i5,2f10.6)') modenum, nu, dyn (modenum, nu)
     enddo
     call stop_ph (.true.)
  endif

  IF ( .NOT. ldiag_loc ) THEN
     DO irr=0,nirr
        IF (.NOT.done_irr(irr)) THEN
           IF (.not.ldisp) THEN
              WRITE(stdout, '(/,5x,"Stopping because representation", &
                                 & i5, " is not done")') irr
              CALL close_phq(.TRUE.)
              CALL stop_ph(.TRUE.)
           ELSE
              WRITE(stdout, '(/5x,"Not diagonalizing because representation", &
                                 & i5, " is not done")') irr
           END IF
           RETURN
        ENDIF
     ENDDO
     ldiag_loc=.TRUE.
  ENDIF
  !
  !   Generates the star of q
  !
  call star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE. )
  !
  ! write on file information on the system
  !
  IF (xmldyn) THEN
     nqq=nq
     IF (imq==0) nqq=2*nq
     IF (lgamma.AND.done_epsil.AND.done_zeu) THEN
        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
             nqq, epsilon, zstareu, lraman, ramtns)
     ELSE
        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau,ityp,m_loc,nqq)
     ENDIF
  ELSE
     CALL write_old_dyn_mat_head(iudyn)
  ENDIF
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  call q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
       nq, sxq, isq, imq, iudyn)

  ! Rotates and write drho_q* (to be improved)
  IF(drho_star%open) THEN
    INQUIRE (UNIT = iudrho, OPENED = opnd)
    IF (opnd) CLOSE(UNIT = iudrho, STATUS='keep')
    CALL write_dfile_star(drho_star, fildrho, nsym, xq, u, nq, sxq, isq, &
                          s, sr, invs, irt, ntyp, ityp,(imq==0), -1 )
  ENDIF
  IF(dvscf_star%open) THEN
    INQUIRE (UNIT = iudvscf, OPENED = opnd)
    IF (opnd) CLOSE(UNIT = iudvscf, STATUS='keep')
    CALL write_dfile_star(dvscf_star, fildvscf, nsym, xq, u, nq, sxq, isq, &
                          s, sr, invs, irt, ntyp, ityp,(imq==0), iq_ )
  ENDIF
  !
  !   Writes (if the case) results for quantities involving electric field
  !
  if (epsil) call write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  IF (zue.AND..NOT.done_zue) THEN
     done_zue=.TRUE.
     IF (lgamma_gamma) THEN
        ALLOCATE(zstar(3,3,nat))
        zstar(:,:,:) = 0.d0
        DO jcart = 1, 3
           DO mu = 1, 3 * nat
              na = (mu - 1) / 3 + 1
              icart = mu - 3 * (na - 1)
              zstar(jcart, icart, na) = zstarue0 (mu, jcart)
           ENDDO
           DO na=1,nat
              work(:)=0.0_DP
              DO icart=1,3
                 work(icart)=zstar(jcart,1,na)*at(1,icart)+ &
                             zstar(jcart,2,na)*at(2,icart)+ &
                             zstar(jcart,3,na)*at(3,icart)
              ENDDO
              zstar(jcart,:,na)=work(:)
           ENDDO
        ENDDO
        CALL generate_effective_charges_c ( nat, nsym, s, invs, irt, at, bg, &
           n_diff_sites, equiv_atoms, has_equivalent, asr, nasr, zv, ityp, &
           ntyp, atm, zstar )
        DO na=1,nat
           do icart=1,3
              zstarue(:,na,icart)=zstar(:,icart,na)
           ENDDO
        ENDDO
        CALL summarize_zue()
        DEALLOCATE(zstar)
     ELSE
        CALL sym_and_write_zue
     ENDIF
  ELSEIF (lgamma) THEN
     IF (done_zue) CALL summarize_zue()
  ENDIF

  if (lraman) call write_ramtns (iudyn, ramtns)
  !
  !   Diagonalizes the dynamical matrix at q
  !
  IF (ldiag_loc) THEN
     call dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, w2)
     IF (search_sym) CALL find_mode_sym (dyn, w2, at, bg, tau, nat, nsymq, sr,&
              irt, xq, rtau, amass, ntyp, ityp, 1, lgamma, lgamma_gamma, &
                                        nspin_mag, name_rap_mode, num_rap_mode)
  END IF
!
! Here we save the dynamical matrix and the effective charges dP/du on
! the recover file. If a recover file with this very high recover code
! is found only the final result is rewritten on output.
!
  rec_code=30
  where_rec='dynmatrix.'
  CALL ph_writefile('status_ph',current_iq,0,ierr)

  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix
!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the vector rtau which contains for each
  !     atom and each rotation the vector S\tau_a - \tau_b, where
  !     b is the rotated a atom, given by the array irt. These rtau are
  !     non zero only if fractional translations are present.
  !
  USE kinds, ONLY : DP
  implicit none
  !
  !     first the dummy variables
  !
  integer, intent(in) :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! nsym: number of symmetries of the point group
  ! s:    matrices of symmetry operations
  ! nat : number of atoms in the unit cell
  ! irt(n,m) = transformed of atom m for symmetry n
  real(DP), intent(in) :: at (3, 3), bg (3, 3), tau (3, nat)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  ! tau: coordinates of the atoms
  logical, intent(in) :: sym (nsym)
  ! sym(n)=.true. if operation n is a symmetry
  real(DP), intent(out):: rtau (3, 48, nat)
  ! rtau: the direct translations
  !
  !    here the local variables
  !
  integer :: na, nb, isym, ipol
  ! counters on: atoms, symmetry operations, polarization
  real(DP) , allocatable :: xau (:,:)
  real(DP) :: ft (3)
  !
  allocate (xau(3,nat))
  !
  !   compute the atomic coordinates in crystal axis, xau
  !
  do na = 1, nat
     do ipol = 1, 3
        xau (ipol, na) = bg (1, ipol) * tau (1, na) + &
                         bg (2, ipol) * tau (2, na) + &
                         bg (3, ipol) * tau (3, na)
     enddo
  enddo
  !
  !    for each symmetry operation, compute the atomic coordinates
  !    of the rotated atom, ft, and calculate rtau = Stau'-tau
  !
  rtau(:,:,:) = 0.0_dp
  do isym = 1, nsym
     if (sym (isym) ) then
        do na = 1, nat
           nb = irt (isym, na)
           do ipol = 1, 3
              ft (ipol) = s (1, ipol, isym) * xau (1, na) + &
                          s (2, ipol, isym) * xau (2, na) + &
                          s (3, ipol, isym) * xau (3, na) - xau (ipol, nb)
           enddo
           do ipol = 1, 3
              rtau (ipol, isym, na) = at (ipol, 1) * ft (1) + &
                                      at (ipol, 2) * ft (2) + &
                                      at (ipol, 3) * ft (3)
           enddo
        enddo
     endif
  enddo
  !
  !    deallocate workspace
  !
  deallocate(xau)
  return
end subroutine sgam_ph
!
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine random_matrix (irt, irgq, nsymq, minus_q, irotmq, nat, &
     wdyn, lgamma)
  !----------------------------------------------------------------------
  !
  !   Create a random hermitian matrix with non zero elements similar to
  !   the dynamical matrix of the system
  !
  !
  USE kinds, only : DP
  USE random_numbers, ONLY : randy
  implicit none
  !
  !    The dummy variables
  !

  integer :: nat, irt (48, nat), irgq (48), nsymq, irotmq
  ! input: number of atoms
  ! input: index of the rotated atom
  ! input: the small group of q
  ! input: the order of the small group
  ! input: the rotation sending q -> -q

  complex(DP) :: wdyn (3, 3, nat, nat)
  ! output: random matrix
  logical :: lgamma, minus_q
  ! input: if true q=0
  ! input: if true there is a symmetry
  !
  !    The local variables
  !
  integer :: na, nb, ipol, jpol, isymq, irot, ira, iramq
  ! counters
  ! ira:   rotated atom
  ! iramq: rotated atom with the q->-q+G symmetry
  !
  !
  wdyn (:, :, :, :) = (0d0, 0d0)
  do na = 1, nat
     do ipol = 1, 3
        wdyn (ipol, ipol, na, na) = CMPLX(2 * randy () - 1, 0.d0,kind=DP)
        do jpol = ipol + 1, 3
           if (lgamma) then
              wdyn (ipol, jpol, na, na) = CMPLX(2 * randy () - 1, 0.d0,kind=DP)
           else
              wdyn (ipol, jpol, na, na) = &
                   CMPLX(2 * randy () - 1, 2 * randy () - 1,kind=DP)
           endif
           wdyn (jpol, ipol, na, na) = CONJG(wdyn (ipol, jpol, na, na) )
        enddo
        do nb = na + 1, nat
           do isymq = 1, nsymq
              irot = irgq (isymq)
              ira = irt (irot, na)
              if (minus_q) then
                 iramq = irt (irotmq, na)
              else
                 iramq = 0
              endif
              if ( (nb == ira) .or. (nb == iramq) ) then
                 do jpol = 1, 3
                    if (lgamma) then
                       wdyn (ipol, jpol, na, nb) = CMPLX(2*randy () - 1, 0.d0,kind=DP)
                    else
                       wdyn (ipol, jpol, na, nb) = &
                            CMPLX(2*randy () - 1, 2*randy () - 1,kind=DP)
                    endif
                    wdyn(jpol, ipol, nb, na) = CONJG(wdyn(ipol, jpol, na, nb))
                 enddo
                 goto 10
              endif
           enddo
10         continue
        enddo
     enddo
  enddo
  return
end subroutine random_matrix

!
! Copyright (C) 2006-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym (u, w2, at, bg, tau, nat, nsym, sr, irt, xq, &
     rtau, amass, ntyp, ityp, flag, lri, lmolecule, nspin_mag,        &
     name_rap_mode, num_rap_mode)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of eigenvectors of the dynamical
  !   matrix. It does NOT work at zone border in non symmorphic space groups.
  !   if flag=1 the true displacements are given in input, otherwise the
  !   eigenvalues of the dynamical matrix are given.
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, RY_TO_CMM1
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_is, ONLY : gname_is
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::             &
       nat,         &
       nsym,        &
       flag,        &
       ntyp,        &
       nspin_mag,   &
       ityp(nat),   &
       irt(48,nat)

  REAL(DP), INTENT(IN) ::   &
       at(3,3),        &
       bg(3,3),        &
       xq(3),          &
       tau(3,nat),     &
       rtau(3,48,nat), &
       amass(ntyp),    &
       w2(3*nat),      &
       sr(3,3,48)

  COMPLEX(DP), INTENT(IN) ::  &
       u(3*nat, 3*nat)       ! The eigenvectors or the displacement pattern
  LOGICAL, INTENT(IN) :: lri      ! if .true. print the Infrared/Raman flag
  LOGICAL, INTENT(IN) :: lmolecule ! if .true. these are eigenvalues of an
                                   ! isolated system

  CHARACTER(15), INTENT(OUT) :: name_rap_mode( 3 * nat )
  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ngroup, &   ! number of different frequencies groups
       nmodes, &   ! number of modes
       imode, imode1, igroup, dim_rap, nu_i,  &
       irot, irap, iclass, mu, na, i

  INTEGER, ALLOCATABLE :: istart(:)

  COMPLEX(DP) :: times              ! safe dimension
  ! in case of accidental degeneracy
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), ALLOCATABLE :: w1(:)
  COMPLEX(DP), ALLOCATABLE ::  rmode(:), trace(:,:), z(:,:)
  LOGICAL :: is_linear
  CHARACTER(3) :: cdum
  INTEGER :: counter, counter_s
  !
  !    Divide the modes on the basis of the mode degeneracy.
  !
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(w1(nmodes))
  ALLOCATE(rmode(nmodes))
  ALLOCATE(trace(48,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvalues of the dynmaical matrix
     !  Note that amass is in amu; amu_ry converts it to Ry au
     !
     DO nu_i = 1, nmodes
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           z (mu, nu_i) = u (mu, nu_i) * SQRT (amu_ry*amass (ityp (na) ) )
        END DO
     END DO
  ELSE
     z=u
  ENDIF

  DO imode=1,nmodes
     w1(imode)=SIGN(SQRT(ABS(w2(imode)))*RY_TO_CMM1,w2(imode))
  ENDDO

  ngroup=1
  istart(ngroup)=1
  imode1=1
  IF (lmolecule) THEN
     istart(ngroup)=7
     imode1=6
     IF(is_linear(nat,tau)) istart(ngroup)=6
  ENDIF
  DO imode=imode1+1,nmodes
     IF (ABS(w1(imode)-w1(imode-1)) > 5.0d-2) THEN
        ngroup=ngroup+1
        istart(ngroup)=imode
     END IF
  END DO
  istart(ngroup+1)=nmodes+1
  !
  !  Find the character of one symmetry operation per class
  !
  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     DO iclass=1,nclass
        irot=elem(1,iclass)
        trace(iclass,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap
           nu_i=istart(igroup)+i-1
           CALL rotate_mod(z(1,nu_i),rmode,sr(1,1,irot),irt,rtau,xq,nat,irot)
           trace(iclass,igroup)=trace(iclass,igroup) + &
                zdotc(3*nat,z(1,nu_i),1,rmode,1)
        END DO
!              write(6,*) igroup,iclass, trace(iclass,igroup)
     END DO
  END DO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  IF (nspin_mag==4) THEN
     IF (flag==1) WRITE(stdout,  &
          '(/,5x,"Mode symmetry, ",a11," [",a11,"] magnetic point group:",/)') &
          gname, gname_is
  ELSE
     IF (flag==1) WRITE(stdout,'(/,5x,"Mode symmetry, ",a11," point group:",/)') gname
  END IF
  num_rap_mode=-1
  counter=1
  DO igroup=1,ngroup
     IF (ABS(w1(istart(igroup)))<1.d-3) CYCLE
     DO irap=1,nclass
        times=(0.d0,0.d0)
        DO iclass=1,nclass
           times=times+CONJG(trace(iclass,igroup))*char_mat(irap, &
                which_irr(iclass))*nelem(iclass)
           !         write(6,*) igroup, irap, iclass, which_irr(iclass)
        ENDDO
        times=times/nsym
        cdum="   "
        IF (lri) cdum=ir_ram(irap)
        IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN
           IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') &
                istart(igroup), istart(igroup+1)-1, w1(istart(igroup))
        ENDIF

        IF (ABS(times) > eps) THEN
           IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
              IF (flag==1) WRITE(stdout,'(5x, "omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",a19)') &
                   istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3*nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           ELSE
              IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",i3,a19)') &
                   istart(igroup), istart(igroup+1)-1, &
                   w1(istart(igroup)), NINT(DBLE(times)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(times))*&
                                            NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3 * nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           END IF
        END IF
     END DO
  END DO
  IF (flag==1) WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(z)
  DEALLOCATE(w1)
  DEALLOCATE(rmode)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym
