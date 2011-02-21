!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmatrix
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
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, pmass, zv
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum
  USE cell_base,     ONLY : at, bg, celldm, ibrav, omega, symm_type
  USE symm_base,     ONLY : s, sr, irt, nsym, time_reversal, invs
  USE printout_base, ONLY : title
  USE dynmat,        ONLY : dyn, w2
  USE qpoint,        ONLY : xq
  USE noncollin_module, ONLY : nspin_mag
  USE modes,         ONLY : u, nmodes, minus_q, irotmq, nsymq, irgq, &
                            rtau, npert, nirr, name_rap_mode, num_rap_mode
  USE gamma_gamma,   ONLY : nasr, asr, equiv_atoms, has_equivalent, &
                            n_diff_sites
  USE efield_mod,    ONLY : epsilon, zstareu, zstarue0, zstarue
  USE control_ph,    ONLY : epsil, zue, lgamma, lgamma_gamma, search_sym, ldisp, &
                            start_irr, last_irr, done_zue, where_rec, &
                            rec_code, ldiag, done_epsil, done_zeu, xmldyn
  USE ph_restart,    ONLY : ph_writefile
  USE partial,       ONLY : all_comp, comp_irr, done_irr, nat_todo_input
  USE units_ph,      ONLY : iudyn
  USE noncollin_module, ONLY : m_loc, nspin_mag
  USE output,        ONLY : fildyn
  USE io_dyn_mat,    ONLY : write_dyn_mat_header
  USE ramanm,        ONLY: lraman, ramtns
  implicit none
  ! local variables
  !
  integer :: nq, isq (48), imq, na, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j, nqq
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real (DP), parameter ::   convfact = BOHR_RADIUS_ANGS**2
  real(DP) :: sxq (3, 48), work(3)
  ! list of vectors in the star of q
  real(DP), allocatable :: zstar(:,:,:)
  integer :: icart, jcart
  logical :: ldiag_loc
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
           if (done_irr (irr) .eq.0.and.done_irr (jrr) .eq.0) then
              do ipert = 1, npert (irr)
                 mu = imode0 + ipert
                 do jpert = 1, npert (jrr)
                    nu = jmode0 + jpert
                    dyn (mu, nu) = CMPLX(0.d0, 0.d0,kind=DP)
                 enddo
              enddo
           elseif (done_irr (irr) .eq.0.and.done_irr (jrr) .ne.0) then
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
        if (comp_irr(irr)==0) then
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
        IF (done_irr(irr)==0) THEN
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
             celldm, at, bg, omega, symm_type, atm, pmass, tau, ityp, m_loc, &
             nqq, epsilon, zstareu, lraman, ramtns*omega/fpi*convfact)
     ELSE
        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, symm_type, atm, pmass, tau,ityp,m_loc,nqq)
     ENDIF
  ELSE
     CALL write_old_dyn_mat(iudyn)
  ENDIF
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  call q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
       nq, sxq, isq, imq, iudyn)
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
     call dyndia (xq, nmodes, nat, ntyp, ityp, pmass, iudyn, dyn, w2)
     IF (search_sym) CALL find_mode_sym (dyn, w2, at, bg, tau, nat, nsymq, sr,&
              irt, xq, rtau, pmass, ntyp, ityp, 1, lgamma, lgamma_gamma, &
                                        nspin_mag, name_rap_mode, num_rap_mode)
  END IF
!
! Here we save the dynamical matrix and the effective charges dP/du on
! the recover file. If a recover file with this very high recover code
! is found only the final result is rewritten on output.
!
  rec_code=30
  where_rec='dynmatrix.'
  CALL ph_writefile('data',0)

  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix

  SUBROUTINE write_old_dyn_mat(iudyn)
!
!  This routine is here for compatibility with the old code.
!  It will be removed when the xml file format of the dynamical matrix
!  will be tested.
!
  USE ions_base, ONLY : ntyp => nsp, nat, ityp, tau, atm, pmass
  USE cell_base, ONLY : ibrav, celldm, at, symm_type
  USE printout_base, ONLY : title

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iudyn
  INTEGER :: nt, na, i, j

  WRITE (iudyn, '("Dynamical matrix file")')
  WRITE (iudyn, '(a)') title
  WRITE (iudyn, '(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  IF (ibrav==0) THEN
     WRITE (iudyn,'(a)') symm_type
     WRITE (iudyn,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  END IF
  DO nt = 1, ntyp
     WRITE (iudyn, * ) nt, ' ''', atm (nt) , ' '' ', pmass (nt)
  ENDDO
  DO na = 1, nat
     WRITE (iudyn, '(2i5,3f15.7)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  ENDDO

  RETURN
  END SUBROUTINE write_old_dyn_mat
!
