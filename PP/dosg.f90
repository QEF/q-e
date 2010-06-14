!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE dos_g (et, nspin, nbnd, nks, wk, Degauss, ngauss, E, dosg)
  !--------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nspin, nks, nbnd, ngauss

  real(DP) :: wk (nks), et (nbnd, nks), Degauss, E, dosg (2)
  real(DP) :: w0gauss
  INTEGER :: n, ns, nk0, nk, ik
  INTEGER :: nspin0
  EXTERNAL w0gauss
  !
  IF (nspin == 1 .or. nspin == 4) THEN
     nk = nks
  ELSE
     nk = nks / 2
  ENDIF
  nspin0=nspin
  IF (nspin==4) nspin0=1
  !
  DO ns = 1, nspin0
     IF (ns==1) THEN
        nk0 = 1
     ELSE
        nk0 = nks / 2 + 1
     ENDIF
     dosg (ns) = 0.0d0
     DO ik = nk0, nk0 + nk-1
        DO n = 1, nbnd
           dosg (ns) = dosg (ns) + wk (ik) * w0gauss ( (E-et (n, ik) ) &
                / Degauss, ngauss)
        ENDDO
     ENDDO
     !
     dosg (ns) = dosg (ns) / Degauss
     !
  ENDDO
  !
  RETURN
END SUBROUTINE dos_g
