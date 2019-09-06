!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ewald_dipole( tens, dipole )
  !-----------------------------------------------------------------------
  !! Calculates the ewald field on each atom due to the presence of dipole, or
  !! the electic field on each atom due to the ionic charge of other atoms,
  !! with both G- and R-space terms. 
  !! Determines optimal alpha. Should hopefully work for any structure.
  !
  USE kinds,      ONLY: DP
  USE gvect,      ONLY: gcutm, gstart, ngm, g, gg
  USE constants,  ONLY: tpi, e2, fpi, pi
  USE cell_base,  ONLY: tpiba2, omega, alat, at, bg
  USE ions_base,  ONLY: nat, ntyp => nsp, ityp, tau
  USE vlocal,     ONLY: strf
  USE mp_bands,   ONLY: intra_bgrp_comm
  USE mp,         ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: tens(nat,3,3)
  !! the ewald field/electric field
  REAL(DP) :: dipole(ntyp)
  !! the dipole
  REAL(DP) :: charge, eta, arg, upperbound, temp
  COMPLEX(DP) :: rhon
  REAL(DP), EXTERNAL :: qe_erfc
  COMPLEX(DP), ALLOCATABLE:: ewaldg(:,:,:), ewaldr(:,:,:)
  INTEGER :: alpha, beta, na, ng, nt, ipol, nb, nrm, nr
  !
  INTEGER, PARAMETER :: mxr = 50
  REAL (DP)   :: r(3,mxr), r2(mxr), rmax, rr, dtau(3)
  REAL (DP)   :: expcoeff
  COMPLEX(DP) :: carg, recarg, recarg_dgg
  !
  ALLOCATE( ewaldg(nat,3,3) )
  ALLOCATE( ewaldr(nat,3,3) )
  !
  ewaldg = (0.d0,0.d0)  
  ewaldr = (0.d0,0.d0)
  !
  !  e2=1.d0 !hartree
  charge = 0.d0
  DO na = 1, nat
     charge = charge + dipole(ityp(na))
  ENDDO
  eta = 2.9d0
  DO
    eta = eta - 0.1d0
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    !
    IF ( eta <= 0.d0 ) CALL errore( 'ewald_dipole', 'optimal eta not found', 1 )
    upperbound = 2.d0 * charge**2 * SQRT(2.d0 * eta / tpi) &
                      * qe_erfc ( SQRT(tpiba2 * gcutm / 4.d0 / eta) )
    IF ( upperbound <= 1.0d-7 ) EXIT
  ENDDO
  !
  ! G-space sum here.
  !
  DO ng = gstart, ngm
     rhon = (0.d0, 0.d0)
     expcoeff = EXP( - gg(ng) * tpiba2 * 0.25d0 / eta )
     DO nt = 1, ntyp
        rhon = rhon + dipole(nt) * CONJG(strf(ng, nt))
     ENDDO
     DO na=1, nat
        arg = (g(1, ng) * tau(1, na) + g(2, ng) * tau(2, na) &
             + g(3, ng) * tau(3, na) ) * tpi
        carg = CMPLX(COS(arg), -SIN(arg),KIND=DP)
        recarg = rhon*expcoeff*carg
        recarg_dgg = recarg / gg(ng)
        DO alpha = 1,3
           DO beta=1,3
              ewaldg(na, alpha, beta) = ewaldg(na, alpha, beta) &
                - recarg_dgg * g(alpha,ng) * g(beta,ng)
           ENDDO
           !
           ewaldg(na, alpha, alpha) = ewaldg(na, alpha, alpha) &
                                      + 1.d0/3.d0 * recarg
        ENDDO
     ENDDO
  ENDDO
  ewaldg = e2 / 2.d0 * fpi / omega * ewaldg !Temp to compare with paratec
  !  ewaldg = e2 * fpi / omega * ewaldg
  !
  CALL mp_sum( ewaldg, intra_bgrp_comm )
  !
  ! R-space sum here (only for the processor that contains G=0)
  !
  ewaldr = 0.d0
  IF ( gstart==2 ) THEN
     rmax = 4.d0 / SQRT(eta) / alat
     !
     ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
     !
    DO na = 1, nat
       DO nb = 1, nat
          DO ipol = 1, 3
             dtau(ipol) = tau(ipol, na) - tau(ipol, nb)
          ENDDO
           !
           ! generates nearest-neighbors shells
           !
          CALL rgen(dtau, rmax, mxr, at, bg, r, r2, nrm)
           !
           ! and sum to the REAL space part
           !
          r = r * alat
          DO nr = 1, nrm
             rr = SQRT( r2(nr) ) * alat
             temp= dipole(ityp(na))  * ( 3.d0 / rr**3 * qe_erfc( SQRT(eta) * rr) &
                           + (6.d0 * SQRT(eta/pi) * 1.d0 / rr*2 + 4.d0 * SQRT(eta**3/pi)) &
                                * EXP(-eta* rr**2))
             DO alpha=1,3
                DO beta=1,3
                   ewaldr(na, alpha,beta) = ewaldr(na, alpha,beta) &
                                           + temp*r(alpha,nr)*r(beta,nr) / rr**2
                ENDDO
                ewaldr(na, alpha,alpha)= ewaldr(na, alpha,alpha) &
                                        - 1.d0/3.d0 * temp
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDIF
 ewaldr = e2 * ewaldr
 !
 CALL mp_sum( ewaldr, intra_bgrp_comm )
 ! 
 tens = ewaldg + ewaldr
 !
END SUBROUTINE ewald_dipole
