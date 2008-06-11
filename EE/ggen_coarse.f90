!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE ggen_coarse(g, gg, ngm)
  !----------------------------------------------------------------------
  !
  !     This routine generates all the reciprocal lattice vectors
  !     contained in the sphere of radius gcutm. Furthermore it
  !     computes the indices nl which give the correspondence
  !     between the fft mesh points and the array of g vectors.
  !
  USE kinds,              ONLY : DP
  USE cell_base,          ONLY : at
  USE gcoarse,            ONLY : ngmc, gcutmc, ngmc_g, nr1c, nr2c, nr3c, &
                                 nrx1c, nrx2c, nrx3c, nlc

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: ngm
  REAL(DP), INTENT(IN) :: g(3,ngm), gg(ngm)
  !
  !     here a few local variables
  !
  INTEGER :: n1c, n2c, n3c
  !
  INTEGER :: ng
  !
  ! counters
     ngmc = 0
     DO ng = 1, ngm
        !
        n1c = NINT ( g(1,ng)*at(1,1) + g(2,ng)*at(2,1) + g(3,ng)*at(3,1) ) + 1
        IF (n1c.LT.1) n1c = n1c + nr1c
        !
        n2c = NINT ( g(1,ng)*at(1,2) + g(2,ng)*at(2,2) + g(3,ng)*at(3,2) ) + 1
        IF (n2c.LT.1) n2c = n2c + nr2c
        !
        n3c = NINT ( g(1,ng)*at(1,3) + g(2,ng)*at(2,3) + g(3,ng)*at(3,3) ) + 1
        IF (n3c.LT.1) n3c = n3c + nr3c
        !
        IF (gg(ng).LE.gcutmc) THEN
           nlc (ng) = n1c + (n2c - 1) * nrx1c + (n3c - 1) * nrx1c * nrx2c
           ngmc = ngmc + 1
        END IF
        !
     ENDDO

     RETURN
   END SUBROUTINE ggen_coarse


