!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE write_rism_type(rismt)
  !--------------------------------------------------------------------------
  !
  ! ... write rism_type for debug.
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_1DRISM, ITYPE_3DRISM, ITYPE_LAUERISM, &
                      & CLOSURE_HNC, CLOSURE_KH
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN) :: rismt
  !
  LOGICAL          :: laue
  INTEGER          :: isite
  CHARACTER(LEN=9) :: stype
  CHARACTER(LEN=3) :: sclosure
  !
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  IF (rismt%itype == ITYPE_1DRISM) THEN
    stype = '1D-RISM'
  ELSE IF (rismt%itype == ITYPE_3DRISM) THEN
    stype = '3D-RISM'
  ELSE IF (rismt%itype == ITYPE_LAUERISM) THEN
    stype = 'Laue-RISM'
  ELSE
    stype = '?????????'
  END IF
  !
  IF (rismt%closure == CLOSURE_HNC) THEN
    sclosure = 'HNC'
  ELSE IF (rismt%closure == CLOSURE_KH) THEN
    sclosure = 'KH '
  ELSE
    sclosure = '???'
  END IF
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"**** RISM data ****")')
  WRITE(stdout, '(5X,"Avairable ?      = ",L1)')    rismt%avail
  WRITE(stdout, '(5X,"Type of data     = ",A)')     TRIM(stype)
  WRITE(stdout, '(5X,"Closure eqn.     = ",A)')     TRIM(sclosure)
  WRITE(stdout, '(5X,"Temperature      = ",E16.8)') rismt%temp
  WRITE(stdout, '(5X,"Smearing radius  = ",E16.8)') rismt%tau
  WRITE(stdout, '(5X,"# solvent site   = ",I10)')   rismt%nsite
  WRITE(stdout, '(5X,"# R-space        = ",I10)')   rismt%nr
  WRITE(stdout, '(5X,"# Z-stick(short) = ",I10)')   rismt%nrzs
  WRITE(stdout, '(5X,"# Z-stick(long)  = ",I10)')   rismt%nrzl
  WRITE(stdout, '(5X,"# G-space        = ",I10)')   rismt%ng
  WRITE(stdout, '(5X,"# G-shell        = ",I10)')   rismt%ngs
  WRITE(stdout, '(5X,"# Gxy-plane      = ",I10)')   rismt%ngxy
  WRITE(stdout, '(5X,"Solvation erg.   = ",E16.8)') rismt%esol
  WRITE(stdout, '(5X,"Solvation pot.   = ",E16.8)') rismt%vsol
  !
  IF (ASSOCIATED(rismt%csr)) THEN
    WRITE(stdout, '(5X,"[Cs(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%csr)
  END IF
  !
  IF (ASSOCIATED(rismt%csg)) THEN
    WRITE(stdout, '(5X,"[Cs(g)]")')
    CALL write_rarray(rismt%ng, rismt%nsite, rismt%csg)
  END IF
  !
  IF (ASSOCIATED(rismt%csgz)) THEN
    IF (.NOT. laue) THEN
      WRITE(stdout, '(5X,"[Cs(g), complex]")')
      CALL write_carray(rismt%ng, rismt%nsite, rismt%csgz)
    ELSE
      WRITE(stdout, '(5X,"[Cs(z,gxy), complex]")')
      CALL write_carray(rismt%nrzs * rismt%ngxy, rismt%nsite, rismt%csgz)
    END IF
  END IF
  !
  IF (ASSOCIATED(rismt%uljr)) THEN
    WRITE(stdout, '(5X,"[Ulj(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%uljr)
  END IF
  !
  IF (ASSOCIATED(rismt%uwr)) THEN
    WRITE(stdout, '(5X,"[Uw(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%uwr)
  END IF
  !
  IF (ASSOCIATED(rismt%usr)) THEN
    WRITE(stdout, '(5X,"[Us(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%usr)
  END IF
  !
  IF (ASSOCIATED(rismt%ulr)) THEN
    WRITE(stdout, '(5X,"[Ul(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%ulr)
  END IF
  !
  IF (ASSOCIATED(rismt%ulg)) THEN
    WRITE(stdout, '(5X,"[Ul(g)]")')
    CALL write_rarray(rismt%ng, rismt%nsite, rismt%ulg)
  END IF
  !
  IF (ASSOCIATED(rismt%ulgz)) THEN
    IF (.NOT. laue) THEN
      WRITE(stdout, '(5X,"[Ul(g), complex]")')
      CALL write_carray(rismt%ng, rismt%nsite, rismt%ulgz)
    ELSE
      WRITE(stdout, '(5X,"[Ul(z,gxy), complex]")')
      CALL write_carray(rismt%nrzl * rismt%ngxy, rismt%nsite, rismt%ulgz)
    END IF
  END IF
  !
  IF (ASSOCIATED(rismt%vlgz)) THEN
    WRITE(stdout, '(5X,"[Vl(z,gxy)]")')
    CALL write_carray(rismt%nrzl * rismt%ngxy, 1, rismt%vlgz)
  END IF
  !
  IF (ASSOCIATED(rismt%vright)) THEN
    WRITE(stdout, '(5X,"[Vright(gxy), complex]")')
    CALL write_carray(rismt%ngxy, 1, rismt%vright)
  END IF
  !
  IF (ASSOCIATED(rismt%vleft)) THEN
    WRITE(stdout, '(5X,"[Vleft(gxy), complex]")')
    CALL write_carray(rismt%ngxy, 1, rismt%vleft)
  END IF
  !
  IF (ASSOCIATED(rismt%do_vright)) THEN
    WRITE(stdout, '(5X,"[DoVright(gxy)]")')
    CALL write_larray(rismt%ngxy, 1, rismt%do_vright)
  END IF
  !
  IF (ASSOCIATED(rismt%do_vleft)) THEN
    WRITE(stdout, '(5X,"[DoVleft(gxy)]")')
    CALL write_larray(rismt%ngxy, 1, rismt%do_vleft)
  END IF
  !
  IF (ASSOCIATED(rismt%hr)) THEN
    WRITE(stdout, '(5X,"[H(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%hr)
  END IF
  !
  IF (ASSOCIATED(rismt%hg)) THEN
    WRITE(stdout, '(5X,"[H(g)]")')
    CALL write_rarray(rismt%ng, rismt%nsite, rismt%hg)
  END IF
  !
  IF (ASSOCIATED(rismt%hgz)) THEN
    IF (.NOT. laue) THEN
      WRITE(stdout, '(5X,"[H(g), complex]")')
      CALL write_carray(rismt%ng, rismt%nsite, rismt%hgz)
    ELSE
      WRITE(stdout, '(5X,"[H(z,gxy), complex]")')
      CALL write_carray(rismt%nrzs * rismt%ngxy, rismt%nsite, rismt%hgz)
    END IF
  END IF
  !
  IF (ASSOCIATED(rismt%hsgz)) THEN
    WRITE(stdout, '(5X,"[Hs(z,gxy), complex]")')
    CALL write_carray(rismt%nrzl * rismt%ngxy, rismt%nsite, rismt%hsgz)
  END IF
  !
  IF (ASSOCIATED(rismt%hlgz)) THEN
    WRITE(stdout, '(5X,"[Hl(z,gxy), complex]")')
    CALL write_carray(rismt%nrzl * rismt%ngxy, rismt%nsite, rismt%hlgz)
  END IF
  !
  IF (ASSOCIATED(rismt%gr)) THEN
    WRITE(stdout, '(5X,"[G(r)]")')
    CALL write_rarray(rismt%nr, rismt%nsite, rismt%gr)
  END IF
  !
  IF (ASSOCIATED(rismt%wg)) THEN
    WRITE(stdout, '(5X,"[W(g)]")')
    CALL write_rarray(rismt%ng, rismt%nsite, rismt%wg)
  END IF
  !
  IF (ASSOCIATED(rismt%zg)) THEN
    WRITE(stdout, '(5X,"[Z(g)]")')
    CALL write_rarray(rismt%ng, rismt%nsite, rismt%zg)
  END IF
  !
  IF (ASSOCIATED(rismt%xgs)) THEN
    DO isite = 1, rismt%mp_site%nsite
      IF (.NOT. laue) THEN
        WRITE(stdout, '(5X,"[X(g)]",I10)') isite
        CALL write_rarray(rismt%ngs, rismt%nsite, rismt%xgs(:, :, isite))
      ELSE
        WRITE(stdout, '(5X,"[X(z,gxy)]",I10)') isite
        CALL write_rarray(rismt%nrzl * rismt%ngs, rismt%nsite, rismt%xgs(:, :, isite))
      END IF
    END DO
  END IF
  !
  IF (ASSOCIATED(rismt%ygs)) THEN
    DO isite = 1, rismt%mp_site%nsite
      IF (.NOT. laue) THEN
        WRITE(stdout, '(5X,"[Y(g)]",I10)') isite
        CALL write_rarray(rismt%ngs, rismt%nsite, rismt%ygs(:, :, isite))
      ELSE
        WRITE(stdout, '(5X,"[Y(z,gxy)]",I10)') isite
        CALL write_rarray(rismt%nrzl * rismt%ngs, rismt%nsite, rismt%ygs(:, :, isite))
      END IF
    END DO
  END IF
  !
  IF (ASSOCIATED(rismt%usol)) THEN
    WRITE(stdout, '(5X,"[Usol]")')
    WRITE(stdout, '(5X,5E16.8)') rismt%usol
  END IF
  !
  IF (ASSOCIATED(rismt%usol_GF)) THEN
    WRITE(stdout, '(5X,"[Usol(GF)]")')
    WRITE(stdout, '(5X,5E16.8)') rismt%usol_GF
  END IF
  !
  IF (ASSOCIATED(rismt%rhog)) THEN
    IF (.NOT. laue) THEN
      WRITE(stdout, '(5X,"[Rho(g)]")')
      CALL write_carray(rismt%ng, 1, rismt%rhog)
    ELSE
      WRITE(stdout, '(5X,"[Rho(z,gxy)]")')
      CALL write_carray(rismt%nrzl * rismt%ngxy, 1, rismt%rhog)
    END IF
  END IF
  !
  IF (ASSOCIATED(rismt%vpot)) THEN
    IF (.NOT. laue) THEN
      WRITE(stdout, '(5X,"[Vpot(g)]")')
      CALL write_carray(rismt%ng, 1, rismt%vpot)
    ELSE
      WRITE(stdout, '(5X,"[Vpot(z,gxy)]")')
      CALL write_carray(rismt%nrzl * rismt%ngxy, 1, rismt%vpot)
    END IF
  END IF
  !
  IF (ASSOCIATED(rismt%rhog_pbc)) THEN
    WRITE(stdout, '(5X,"[Rho_PBC(g)]")')
    CALL write_carray(rismt%ng, 1, rismt%rhog_pbc)
  END IF
  !
  IF (ASSOCIATED(rismt%vpot_pbc)) THEN
    WRITE(stdout, '(5X,"[Vpot_PBC(g)]")')
    CALL write_carray(rismt%ng, 1, rismt%vpot_pbc)
  END IF
  !
  WRITE(stdout, '()')
  !
CONTAINS
  !
  SUBROUTINE write_rarray(nt, nsite, rt)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: nt
    INTEGER,  INTENT(IN) :: nsite
    REAL(DP), INTENT(IN) :: rt(nt,1:*)
    !
    INTEGER :: isite
    INTEGER :: nt1
    INTEGER :: nt2
    !
    INTEGER, PARAMETER :: MT = 10
    !
    nt1 = MIN(MT, nt)
    nt2 = MAX(nt - MT + 1, nt1 + 1)
    !
    DO isite = 1, nsite
      !
      IF (nsite > 1) THEN
        WRITE(stdout, '(5X,"#site =",I10)') isite
      END IF
      !
      WRITE(stdout, '(5X,5E16.8)') rt(1:nt1, isite)
      IF (nt2 <= nt) THEN
        WRITE(stdout, '(5X,5("  .............."))')
        WRITE(stdout, '(5X,5E16.8)') rt(nt2:nt, isite)
      END IF
      !
    END DO
  END SUBROUTINE write_rarray
  !
  SUBROUTINE write_carray(nt, nsite, ct)
    IMPLICIT NONE
    INTEGER,     INTENT(IN) :: nt
    INTEGER,     INTENT(IN) :: nsite
    COMPLEX(DP), INTENT(IN) :: ct(nt,1:*)
    !
    INTEGER :: isite
    INTEGER :: nt1
    INTEGER :: nt2
    !
    INTEGER, PARAMETER :: MT = 10
    !
    nt1 = MIN(MT, nt)
    nt2 = MAX(nt - MT + 1, nt1 + 1)
    !
    DO isite = 1, nsite
      !
      IF (nsite > 1) THEN
        WRITE(stdout, '(5X,"#site =",I10)') isite
      END IF
      !
      WRITE(stdout, '(5X,"real:")')
      WRITE(stdout, '(5X,5E16.8)') DBLE(ct(1:nt1, isite))
      IF (nt2 <= nt) THEN
        WRITE(stdout, '(5X,5("  .............."))')
        WRITE(stdout, '(5X,5E16.8)') DBLE(ct(nt2:nt, isite))
      END IF
      !
      WRITE(stdout, '(5X,"imag:")')
      WRITE(stdout, '(5X,5E16.8)') AIMAG(ct(1:nt1, isite))
      IF (nt2 <= nt) THEN
        WRITE(stdout, '(5X,5("  .............."))')
        WRITE(stdout, '(5X,5E16.8)') AIMAG(ct(nt2:nt, isite))
      END IF
      !
    END DO
  END SUBROUTINE write_carray
  !
  SUBROUTINE write_larray(nt, nsite, lt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nt
    INTEGER, INTENT(IN) :: nsite
    LOGICAL, INTENT(IN) :: lt(nt,1:*)
    !
    INTEGER :: isite
    INTEGER :: nt1
    INTEGER :: nt2
    !
    INTEGER, PARAMETER :: MT = 10
    !
    nt1 = MIN(MT, nt)
    nt2 = MAX(nt - MT + 1, nt1 + 1)
    !
    DO isite = 1, nsite
      !
      IF (nsite > 1) THEN
        WRITE(stdout, '(5X,"#site =",I10)') isite
      END IF
      !
      WRITE(stdout, '(5X,5L4)') lt(1:nt1, isite)
      IF (nt2 <= nt) THEN
        WRITE(stdout, '(5X,5("  .."))')
        WRITE(stdout, '(5X,5L4)') lt(nt2:nt, isite)
      END IF
      !
    END DO
  END SUBROUTINE write_larray
  !
END SUBROUTINE write_rism_type
