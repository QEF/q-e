!
! Copyright (C) 2004 Tone Kokalj
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_whole_cell (alat, at, nat, tau, atm, ityp, &
     nr1, nr2, nr3, nr1x, nr2x, nr3x, rho, output_format, ounit)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER          :: nat, ityp (nat), output_format, ounit
  INTEGER          :: nr1x, nr2x, nr3x, nr1, nr2, nr3
  CHARACTER(len=3) :: atm(*)
  real(DP)    :: alat, tau (3, nat), at (3, 3), rho(2, nr1x,nr2x,nr3x)

  IF ( output_format == 3 ) THEN
     !
     ! XCRYSDEN FORMAT
     !
     CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
     CALL xsf_fast_datagrid_3d &
          (rho, nr1, nr2, nr3, nr1x, nr2x, nr3x, at, alat, ounit)

  ELSEIF ( output_format == 4 ) THEN
     !
     ! gOpenMol format
     !

     ! not yet implemented
     ! add code here ...
  ELSE
     CALL errore('plot_whole_cell', 'wrong output_format', 1)
  ENDIF
END SUBROUTINE plot_whole_cell
