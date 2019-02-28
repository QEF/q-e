!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE lr_smallgq (xq)
  !----------------------------------------------------------------------
  !
  ! This subroutine finds the small group of q: S q = q + G,
  ! i.e. selects among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave q unchanged up to G.
  !
  ! Needed on input (read from data file):
  ! "nsym" crystal symmetries s, "nrot" lattice symetries s.
  !
  ! Inspired by PH/smallg_q, PH/smallgq, and PH/setup_nscf.
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,              ONLY : DP
  USE cell_base,          ONLY : at, bg
  USE noncollin_module,   ONLY : noncolin
  USE symm_base,          ONLY : s, nrot, nsym, sname, copy_sym, s_axis_to_cart
  USE control_flags,      ONLY : noinv

  USE lr_symm_base, ONLY : nsymq, invsymq, gi, minus_q
  !
#if defined(__MPI)
  USE mp,                 ONLY : mp_bcast
  USE io_global,          ONLY : ionode_id
  USE mp_global,          ONLY : intra_image_comm
#endif

  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: accep = 1.e-5_dp
  REAL(DP), INTENT(in) :: xq(3)
  ! input: the q point of the perturbation
  ! 
  LOGICAL :: sym(48)
  ! .true. if symm. op. S q = q + G
  !
  !  local variables
  !
  REAL(DP) :: aq(3), raq(3), wrk(3), zero(3)
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! additional space to compute gi
  ! the zero vector
  !
  INTEGER :: irot, isym, ipol, jpol
  ! counters
  !
  LOGICAL :: eqvect
  ! logical function, check if two vectors are equivalent
  ! 
  IF ( nsym == 1 ) THEN
     !
     nsymq = 1
     RETURN
     !
  ENDIF
  !
  CALL start_clock ('lr_smallgq')
  !
  zero(:) = 0.d0
  !
  ! Transform xq to the crystal basis
  !
  aq(:) = xq(:)
  CALL cryst_to_cart (1, aq, at, - 1)
  !
  ! Test all symmetries to see if this operation send Sq in q+G.
  !
  sym(1:nsym) = .true.
  sym(nsym+1:nrot) = .false.
  !
  DO irot = 1, nrot
     !
     IF (.not.sym(irot) ) goto 100
     !
     ! Rotate q : S*q
     !
     raq(:) = 0.0d0
     DO ipol = 1, 3
        DO jpol = 1, 3
           raq(ipol) = raq(ipol) + DBLE( s(ipol,jpol,irot) ) * aq(jpol)
        ENDDO
     ENDDO
     !
     sym(irot) = eqvect(raq, aq, zero, accep)
     !
     IF (sym(irot)) THEN
        !
        raq = - raq
        minus_q = eqvect (raq, aq, zero, accep)
        if (minus_q) CALL errore( 'lr_smalgq', 'minus_q=.true., &
                        & bug, do not use symmetry for this q!', 1 )
        !
     ENDIF
     !
100  CONTINUE
     !
  ENDDO
  !
  ! Re-order all rotations in such a way that true sym.ops.
  ! are the first nsymq; rotations that are not sym.ops. follow.
  !
  nsymq = copy_sym ( nsym, sym )
  !
  ! Determine the G-vectors : G = S*q - q 
  !
  gi(:,:) = 0.0d0
  !
  DO isym = 1, nsymq
     !
     ! Rotate q : S*q
     !
     raq(:) = 0.0d0
     DO ipol = 1, 3
        DO jpol = 1, 3
           raq(ipol) = raq(ipol) + DBLE( s(ipol,jpol,isym) ) * aq(jpol)
        ENDDO
     ENDDO
     !
     ! G = S*q - q
     !
     DO ipol = 1, 3
        wrk(ipol) = raq(ipol) - aq(ipol)
     ENDDO
     !
     ! Transform G to the cartesian basis
     !
     CALL cryst_to_cart (1, wrk, bg, 1)
     !
     gi(:,isym) = wrk(:)
     !
  ENDDO
  !
  ! Check if inversion (I) is a symmetry. If so, there should be nsymq/2
  ! symmetries without inversion, followed by nsymq/2 with inversion
  ! Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  ! IT: it seems that invsymq is useless (used nowhere)...
  !
  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
  !
  ! The order of the s matrices has changed.
  ! Transform symmetry matrices s from crystal to cartesian axes.
  !
  CALL s_axis_to_cart ()
  !
  minus_q = .false.
  !
  ! Parallel stuff: the first node broadcasts everything to all nodes.
  ! Copied from PH/set_irr.f90
  !
#if defined(__MPI)
  CALL mp_bcast (nsymq, ionode_id, intra_image_comm)
  CALL mp_bcast (gi, ionode_id, intra_image_comm)
  CALL mp_bcast (minus_q, ionode_id, intra_image_comm)
  CALL mp_bcast (invsymq, ionode_id, intra_image_comm)
#endif
  !
  CALL stop_clock ('lr_smallgq')
  !
  RETURN
  !
END SUBROUTINE lr_smallgq
