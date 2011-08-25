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
                    nspin_mag, t_rev, pmass, num_rap_mode, name_rap_mode)
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
#ifdef __PARA
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
       gi (3, 48), gimq (3), sr(3,3,48), tau(3,nat), pmass(ntyp)
! input: the q point
! input: the R associated to each tau
! input: the direct lattice vectors
! input: the reciprocal lattice vectors
! output: [S(irotq)*q - q]
! output: [S(irotmq)*q + q]
! input: symmetry matrices in cartesian coordinates
! input: the atomic positions
! input: the mass of each atom

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
                  sr, irt, xq, rtau, pmass, ntyp, ityp, 0, lgamma, &
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

#ifdef __PARA
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
