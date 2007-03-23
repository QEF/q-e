!
! Copyright (C) 2001 PWSCF group
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
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, amass
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum, noinv
  USE cell_base,     ONLY : at, bg, celldm, ibrav, symm_type
  USE gvect,         ONLY : nr1, nr2, nr3
  USE symme,         ONLY : s, irt, nsym
  USE printout_base, ONLY : title
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : domag
  use phcom
  USE ramanm,        ONLY: lraman, ramtns
  implicit none
  ! local variables
  INTEGER :: s_ (3, 3, 48), invs_ (48), nsym_
  INTEGER, ALLOCATABLE :: irt_ (:,:)
  REAL (KIND=DP), ALLOCATABLE :: rtau_ (:,:,:)
  !
  ! s, irt, nsym, rtau, invs are recalculated by star_q because the full
  ! crystal symmetry is needed. Local variables are used to prevent trouble 
  ! with subsequent electron-phonon symmetrization. 
  ! FIXME: both the full symmetry group and the small group of q should be
  !        stored and clearly distinguishable - see also elphon
  integer :: nq, isq (48), imq, na, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real(DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  call start_clock('dynmatrix')
  ALLOCATE ( rtau_ (3, 48, nat), irt_ (48, nat) )
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
  !
  !   Generates the star of q
  !
  call star_q (xq, at, bg, ibrav, symm_type, nat, tau, ityp, nr1, &
       nr2, nr3, nsym_, s_, invs_, irt_, rtau_, nq, sxq, isq, imq,&
       noinv, modenum, noncolin, domag)
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
     write (iudyn, * ) nt, ' ''', atm (nt) , ' '' ', amass (nt)
  enddo
  do na = 1, nat
     write (iudyn, '(2i5,3f15.7)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  enddo
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  call q2qstar_ph (dyn, at, bg, nat, nsym_, s_, invs_, irt_, rtau_, &
       nq, sxq, isq, imq, iudyn)
  DEALLOCATE ( irt_, rtau_ )
  !
  !   Writes (if the case) results for quantities involving electric field
  !
  if (epsil) call write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  if (zue) call sym_and_write_zue
  if (lraman) call write_ramtns (iudyn, ramtns)
  !
  !   Diagonalizes the dynamical matrix at q
  !
  IF (all_comp) THEN
     call dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, w2)
     IF (search_sym) CALL find_mode_sym (dyn, w2, at, bg, nat, nsym, s, irt, &
                                     xq, rtau, amass, ntyp, ityp)
  END IF
  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix
