!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmat0_new
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
  USE kinds,         ONLY : DP
  USE ions_base, ONLY : nat,ntyp => nsp, ityp, zv, tau
  USE cell_base, ONLY: alat, omega, at, bg
  USE gvect, ONLY: g, gg, ngm, gcutm
  USE symm_base, ONLY: irt, s, invs
  USE control_flags, ONLY : modenum, llondon
  USE ph_restart,    ONLY : ph_writefile
  USE control_ph,    ONLY : rec_code_read, current_iq
  USE qpoint,        ONLY : xq
  USE modes,         ONLY : u, nmodes
  USE partial,       ONLY : done_irr, comp_irr
  USE dynmat,        ONLY : dyn, dyn00, dyn_rec
  USE london_module, ONLY : init_london, dealloca_london
  USE lr_symm_base,  ONLY : minus_q, irotmq, nsymq, rtau

  implicit none

  integer :: nu_i, nu_j, na_icart, nb_jcart, ierr
  ! counters

  complex(DP) :: wrk, dynwrk (3 * nat, 3 * nat)
  ! auxiliary space

  IF ( .NOT. comp_irr(0) .or. done_irr(0) ) RETURN
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
  !   Here the Grimme D2 dispersion contribution (if present)
  !
  IF ( llondon ) THEN
     CALL init_london ()
     CALL d2ionq_mm (alat, nat, ityp, at, bg, tau, xq, dynwrk )
     CALL rotate_pattern_add ( nat, u, dyn, dynwrk )
     CALL dealloca_london ()
  END IF
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

     call symdyn_munu_new (dyn, u, xq, s, invs, rtau, irt, at, bg, &
          nsymq, nat, irotmq, minus_q)
     !
     ! rotate again in the pattern basis
     !
     call zcopy (9 * nat * nat, dyn, 1, dynwrk, 1)

     dyn=(0.d0, 0.d0)

     CALL rotate_pattern_add(nat, u, dyn, dynwrk)

  endif
  !      call tra_write_matrix('dynmat0 dyn',dyn,u,nat)
  dyn_rec(:,:)=dyn(:,:)
  done_irr(0) = .TRUE.
  CALL ph_writefile('data_dyn',current_iq,0,ierr)

  call stop_clock ('dynmat0')
  return
end subroutine dynmat0_new
