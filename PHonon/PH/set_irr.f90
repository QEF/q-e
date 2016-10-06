!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_irr_new (xq, u, npert, nirr, eigen) 
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
!     Original routine was from C. Bungaro.
!     Revised Oct. 1995 by Andrea Dal Corso.
!     April 1997: parallel stuff added (SdG)
!
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  USE ions_base, ONLY : nat, tau, ntyp => nsp, ityp, amass
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, sr, ftau, invs, nsym, irt, t_rev
  USE modes,     ONLY : num_rap_mode, name_rap_mode
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE spin_orb,  ONLY : domag
  USE constants, ONLY: tpi
  USE control_ph, ONLY : search_sym
  USE control_flags, ONLY : iverbosity
  USE random_numbers, ONLY : randy
  USE rap_point_group, ONLY : name_rap

  use mp, only: mp_bcast
  use io_global, only : ionode_id
  use mp_images, only : intra_image_comm

  USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, gi, gimq, rtau
  USE control_lr,   ONLY : lgamma

  implicit none
!
!   first the dummy variables
!
  real(DP), INTENT(IN) :: xq (3)
! input: the q point

  complex(DP), INTENT(OUT) :: u(3*nat, 3*nat)

  INTEGER, INTENT(OUT) :: npert(3*nat), nirr

  REAL(DP), INTENT(OUT) :: eigen(3*nat)
  
!
!   here the local variables
!
  integer :: na, nb, imode, jmode, ipert, jpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna, isym
  ! counters and auxiliary variables

  integer :: info, mode_per_rap(0:12), count_rap(0:12), rap, init, pos, irap, &
             num_rap_aux( 3 * nat ), ierr

  real(DP) :: modul, arg, eig(3*nat)
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

  logical :: magnetic_sym

   magnetic_sym=noncolin.AND.domag
!
!   then we generate a random hermitean matrix
!
     arg = randy(0)
     call random_matrix_new (irt,nsymq,minus_q,irotmq,nat,wdyn,lgamma)
!call write_matrix('random matrix',wdyn,nat)
!
! symmetrize the random matrix with the little group of q
!
     call symdynph_gq_new (xq,wdyn,s,invs,rtau,irt,nsymq,nat,irotmq,minus_q)
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
     CALL compact_dyn(nat, phi, wdyn)
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
        npert=1
        do imode=1,3*nat
           WRITE( stdout, '(2x,"autoval = ", e10.4)') eigen(imode)
           CALL write_modes_out(imode,imode-1)
        end do
     end if

     IF (search_sym) THEN
        CALL find_mode_sym_new (u, eigen, tau, nat, nsymq, s, sr, irt, xq,    &
             rtau, amass, ntyp, ityp, 0, .FALSE., .TRUE., num_rap_mode, ierr)

!
!   Order the modes so that we first make all those that belong to the first
!   representation, then the second ect.
!
!
!   First count, for each irreducible representation, how many modes
!   belong to that representation. Modes that could not be classified
!   have num_rap_mode 0.
!
        mode_per_rap=0
        DO imode=1,3*nat
           mode_per_rap(num_rap_mode(imode))= &
                           mode_per_rap(num_rap_mode(imode))+1
        ENDDO
!
!   The position of each mode on the list is the following:
!   The positions from 1 to mode_per_rap(0) contain the modes that transform 
!   according to the first representation. From mode_per_rap(1)+1 to 
!   mode_per_rap(1) + mode_per_rap(2) the mode that transform according 
!   to the second ecc.
!
        count_rap=1
        DO imode=1,3*nat
           rap=num_rap_mode(imode)
           IF (rap>12) call errore('set_irr',&
                                   'problem with the representation',1)
!
!   Determine the first position for the representation rap
!
           init=0
           DO irap=0,rap-1
              init=init+mode_per_rap(irap)
           ENDDO
!
!   Determine in which position to put this mode. count_rap keep into
!   account how many modes of that representation we have already
!   assigned
!
           pos=init+count_rap(rap)
!
!   the eigenvalue, the mode and the number of its representation are
!   copied in the auxiliary list 
!
!
           eig(pos)=eigen(imode)
           phi(:,pos)=u(:,imode)
           num_rap_aux(pos)=num_rap_mode(imode)
!
!   Update the number of modes already found for a representation
!
           count_rap(rap)=count_rap(rap)+1
        ENDDO
!
!  Copy the new exchanged array in the old ones
!
        eigen=eig
        u=phi
        num_rap_mode=num_rap_aux
!
!  If two almost degenerate modes have been assigned to different
!  representations, we force them to be close in the list independently
!  from their representation in order not to change previous behaviour
!  of the code. These instructions should not be needed.
!
        DO imode=1,3*nat-1
           DO jmode = imode+1, 3*nat
              IF ((num_rap_mode(imode) /= num_rap_mode(jmode)).AND.  &
                  (ABS(eigen(imode) - eigen(jmode))/   &
                  (ABS(eigen(imode)) + ABS (eigen (jmode) )) < 1.d-4) ) THEN
                 WRITE(stdout,'("Eigenvectors exchange needed",2i5)') imode, &
                                                     jmode
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
!
!  Here we set the name of the representation for each mode
!
        name_rap_mode=' '
        DO imode = 1, 3*nat
           IF (num_rap_mode(imode) > 0 ) &
              name_rap_mode(imode)=name_rap(num_rap_mode(imode))
        ENDDO
     ENDIF
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
400 continue
  call mp_bcast (gi, ionode_id, intra_image_comm)
  call mp_bcast (gimq, ionode_id, intra_image_comm)
  call mp_bcast (u, ionode_id, intra_image_comm)
  call mp_bcast (nsymq, ionode_id, intra_image_comm)
  call mp_bcast (npert, ionode_id, intra_image_comm)
  call mp_bcast (nirr, ionode_id, intra_image_comm)
  call mp_bcast (irotmq, ionode_id, intra_image_comm)
  call mp_bcast (minus_q, ionode_id, intra_image_comm)
  call mp_bcast (num_rap_mode, ionode_id, intra_image_comm)
  call mp_bcast (name_rap_mode, ionode_id, intra_image_comm)

  return
end subroutine set_irr_new
