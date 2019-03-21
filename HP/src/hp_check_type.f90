!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine hp_check_type(na)
  !----------------------------------------------------------------------
  !
  ! If the symmetry is not used, then this routine does nothing.
  ! Instead if the symmetry is used, then this routine can 
  ! recompute the symmetries in the ground state, depending on
  ! the following condition:
  ! - if the perturbed Hubbard atom is unique by its type, then
  !   this routine will not recompute the symmetries and use
  !   the original nsym computed by PWscf;
  ! - if the perturbed Hubbard atom has the same type as some 
  !   other atoms in the system, then change the type of the 
  !   perturbed atom in order to make it inequivalent to all 
  !   other atoms (this is needed in order to mimic an isolated
  !   perturbation in a supercell approach).
  ! 
  USE ions_base,          ONLY : ityp, nat, ntyp => nsp, tau
  USE io_global,          ONLY : stdout
  USE symm_base,          ONLY : nsym, set_sym, ft, ftau
  USE noncollin_module,   ONLY : nspin_mag, m_loc
  USE fft_base,           ONLY : dfftp
  USE ldaU_hp,            ONLY : recalc_sym
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na ! the atom under consideration
  !
  INTEGER :: nt, ityp_save, nsym_old
  INTEGER :: na_, nt_, isym
  !
  IF (nsym==1) RETURN
  !
  nt = ityp(na)
  !
  recalc_sym = .false.
  !
  ! If there are at least two Hubbard atoms of the same type
  ! and we want to perturb one of them, then we have to change
  ! the type of the perturbed atom (in order to make them inequivalent)
  ! and to recompute the symmetries.
  !
  DO na_ = 1, nat 
     !
     nt_ = ityp(na_)
     !
     IF ( na_.NE.na .AND. nt_.EQ.nt ) THEN
        !
        WRITE( stdout, '(/5x,"The perturbed atom has a type which is not unique!")')
        WRITE( stdout, '(5x,"Changing the type of the perturbed atom and recomputing the symmetries...")')
        !
        ! Save the type of the perturbed atom na
        !
        ityp_save = ityp(na) 
        !
        ! Change temporarily the type of the perturbed atom na
        ! in order to recompute the symmetry
        !
        ityp(na) = ntyp + 1
        !
        recalc_sym = .true.
        !
        GO TO 100  
        !    
     ENDIF
     !
  ENDDO
  !
  WRITE( stdout, '(/5x,"The perturbed atom has a type which is unique!")')
  !
100 CONTINUE
  !
  ! Recalculate the symmetry of the unperturbed lattice (if needed)
  !
  IF (recalc_sym) THEN
     !
     nsym_old = nsym ! number of symmetries from PWscf
     !
     IF (.not.ALLOCATED(m_loc)) ALLOCATE(m_loc(3,nat))
     m_loc(:,:) = 0.0d0 
     !
     CALL set_sym (nat, tau, ityp, nspin_mag, m_loc)
     !
     DEALLOCATE(m_loc)
     !
     ! Since symmetries were recomputed, we need to reinitialize vectors
     ! of fractional translations
     !
     DO isym = 1, nsym
        ftau(1,isym) = NINT( ft(1,isym) * DBLE(dfftp%nr1) )
        ftau(2,isym) = NINT( ft(2,isym) * DBLE(dfftp%nr2) )
        ftau(3,isym) = NINT( ft(3,isym) * DBLE(dfftp%nr3) )
     ENDDO
     !
     IF ( nsym == nsym_old ) THEN
        WRITE( stdout, '(5x,"The number of symmetries is the same as in PWscf :")')
        recalc_sym = .false.
     ELSE
        WRITE( stdout, '(5x,"The number of symmetries is reduced :")')
     ENDIF
     WRITE( stdout, '(5x,"nsym =",1x,i2,2x,"nsym_PWscf =",1x,i2)') nsym, nsym_old
     WRITE( stdout, '(5x,"Changing the type of the perturbed atom back to its original type...")') 
     !
     ! Change back the type of the perturbed atom to its original type
     !
     ityp(na) = ityp_save
     !
  ENDIF
  !
  RETURN
  !
end subroutine hp_check_type
