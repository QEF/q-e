!
! Copyright (C) 2001-2008 QUANTUM-ESPRESSO group
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
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, pmass, zv
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum
  USE cell_base,     ONLY : at, bg, celldm, ibrav, symm_type
  USE gvect,         ONLY : nr1, nr2, nr3
  USE symme,         ONLY : s, irt, nsym, time_reversal
  USE printout_base, ONLY : title
  USE dynmat,        ONLY : dyn, w2
  USE qpoint,        ONLY : xq
  USE modes,         ONLY : u, nmodes, minus_q, irotmq, nsymq, irgq, invs, &
                            rtau, npert, nirr
  USE gamma_gamma,   ONLY : nasr, asr, equiv_atoms, has_equivalent, &
                            n_diff_sites
  USE efield_mod,    ONLY : epsilon, zstareu, zstarue0
  USE control_ph,    ONLY : epsil, zue, lgamma_gamma, search_sym, ldisp
  USE partial,       ONLY : all_comp, comp_irr, done_irr
  USE units_ph,      ONLY : iudyn
  USE ramanm,        ONLY: lraman, ramtns
  implicit none
  ! local variables
  !
  integer :: nq, isq (48), imq, na, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real(DP) :: sxq (3, 48), work(3)
  ! list of vectors in the star of q
  real(DP), allocatable :: zstar(:,:,:)
  integer :: icart, jcart
  !
  call start_clock('dynmatrix')
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
                    dyn (mu, nu) = CMPLX (0.d0, 0.d0)
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
     CALL generate_dynamical_matrix &
          (nat,nsym,s,irt,at,bg, n_diff_sites,equiv_atoms,has_equivalent,dyn)
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

  DO irr=0,nirr
     IF (done_irr(irr)==0.and..not.ldisp) THEN
        WRITE(stdout,&
        '(/,5x,"Stopping because representation", i5, " is not done")') irr
        CALL close_phq(.TRUE.)
        CALL stop_ph(.TRUE.)
     ENDIF
     IF (done_irr(irr)==0.and.ldisp) THEN
        WRITE(stdout, '(/5x,"Not diagonalizing because representation", &
                         & i5, " is not done")') irr
        RETURN
     ENDIF
  ENDDO
  !
  !   Generates the star of q
  !
  call star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq )
  !
  ! write on file information on the system
  !
  write (iudyn, '("Dynamical matrix file")') 
  write (iudyn, '(a)') title
  write (iudyn, '(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  if (ibrav==0) then
     write (iudyn,'(a)') symm_type
     write (iudyn,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  end if
  do nt = 1, ntyp
     write (iudyn, * ) nt, ' ''', atm (nt) , ' '' ', pmass (nt)
  enddo
  do na = 1, nat
     write (iudyn, '(2i5,3f15.7)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  enddo
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  call q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
       nq, sxq, isq, imq, iudyn)
  !
  !   Writes (if the case) results for quantities involving electric field
  !
  if (epsil) call write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  IF (zue) THEN
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
        CALL generate_effective_charges_c (nat,nsym,s,irt,at,bg,n_diff_sites,&
                    equiv_atoms, has_equivalent,asr,nasr,zv,ityp,ntyp,atm,zstar)
        DEALLOCATE(zstar)
     ELSE
        CALL sym_and_write_zue
     ENDIF
  ENDIF

  if (lraman) call write_ramtns (iudyn, ramtns)
  !
  !   Diagonalizes the dynamical matrix at q
  !
  IF (all_comp) THEN
     call dyndia (xq, nmodes, nat, ntyp, ityp, pmass, iudyn, dyn, w2)
     IF (search_sym) CALL find_mode_sym (dyn, w2, at, bg, nat, nsymq, s, irt, &
                                     xq, rtau, pmass, ntyp, ityp)
  END IF
  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix
