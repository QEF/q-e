!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                gg, ngm, gcutm, gstart, gamma_only, strf )
  !-----------------------------------------------------------------------
  !! Calculates Ewald energy with both G- and R-space terms. 
  !! Determines optimal alpha. Should hopefully work for any structure.
  !
  USE kinds
  USE constants,         ONLY : tpi, e2
  USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,                ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_ewald, do_comp_mt
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_ewald
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the unit cell
  INTEGER, INTENT(IN) :: ntyp
  !! number of different types of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! the type of each atom
  INTEGER, INTENT(IN) :: ngm
  !! number of plane waves for G sum
  INTEGER, INTENT(IN) :: gstart
  !! first non-zero G vector
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma only
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! the positions of the atoms in the cell
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! the coordinates of G vectors
  REAL(DP), INTENT(IN) :: gg(ngm)
  !! the square moduli of G vectors
  REAL(DP), INTENT(IN) :: zv(ntyp)
  !! the charge of each type of atoms
  REAL(DP), INTENT(IN) :: at(3,3)
  !! the direct lattice vectors
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! the reciprocal lattice vectors
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: alat
  !! lattice parameter
  REAL(DP), INTENT(IN) :: gcutm
  !! cut-off of g vectors
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  !! structure factor
  REAL(DP) :: ewald
  !! output: the ewald energy
  !
  !  ... local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r
  INTEGER :: ng, nr, na, nb, nt, nrm
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atomic types
  ! number of R vectors included in r sum
  INTEGER :: na_s, na_e, mykey
  ! for parallelization of real-space sums
  !
  REAL(DP) :: charge, tpiba2, ewaldg, ewaldr, dtau(3), alpha, &
              r(3,mxr), r2(mxr), rmax, rr, upperbound, fact
  ! total ionic charge in the cell
  ! length in reciprocal space
  ! ewald energy computed in reciprocal space
  ! ewald energy computed in real space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  COMPLEX(DP) :: rhon
  !
  tpiba2 = (tpi / alat) **2
  charge = 0.d0
  DO na = 1, nat
     charge = charge + zv( ityp(na) )
  ENDDO
  alpha = 2.9d0
100 alpha = alpha - 0.1d0
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  IF ( alpha <= 0.d0 ) CALL errore( 'ewald', 'optimal alpha not found', 1 )
  upperbound = 2.d0 * charge**2 * SQRT(2.d0 * alpha / tpi) * erfc( &
       SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
  IF (upperbound > 1.0d-7) GOTO 100
  !
  ! G-space sum here.
  ! Determine if this processor contains G=0 and set the constant term
  !
  IF ( do_cutoff_2D ) THEN ! cutoff ewald sums
     ewaldg = cutoff_ewald( gamma_only, alpha, omega )
  ELSE
     IF ( gstart==2 ) THEN
        ewaldg = - charge**2 / alpha / 4.0d0
     ELSE
        ewaldg = 0.0d0
     ENDIF
     IF ( gamma_only ) THEN
        fact = 2.d0
     ELSE
        fact = 1.d0
     ENDIF
     DO ng = gstart, ngm
        rhon = (0.d0, 0.d0)
        DO nt = 1, ntyp
           rhon = rhon + zv(nt) * CONJG(strf(ng, nt))
        ENDDO
        ewaldg = ewaldg + fact * ABS(rhon)**2 * EXP( - gg(ng) * tpiba2 / &
                 alpha / 4.d0) / gg(ng) / tpiba2
     ENDDO
     ewaldg = 2.d0 * tpi / omega * ewaldg
     !
     !  Here add the other constant term
     !
     IF (gstart==2) THEN
        DO na = 1, nat
           ewaldg = ewaldg - zv(ityp(na))**2 * SQRT(8.d0 / tpi * &
                alpha)
        ENDDO
     ENDIF
  ENDIF
  !
  ! R-space sum here
  !
  ! poor-man parallelization over atoms
  ! - if nproc_bgrp=1   : na_s=1, na_e=nat, mykey=0
  ! - if nproc_bgrp<=nat: each processor calculates atoms na_s to na_e; mykey=0
  ! - if nproc_bgrp>nat : each processor takes care of atom na_s=na_e;
  !   mykey labels how many times each atom appears (mykey=0 first time etc.)
  !
  CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
  ewaldr = 0.d0
  IF ( mykey == 0 ) THEN
     rmax = 4.d0 / SQRT(alpha) / alat
     !
     ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
     !
     DO na = na_s, na_e
        DO nb = 1, nat
           dtau(:) = tau(:,na) - tau(:,nb)
           !
           ! generates nearest-neighbors shells
           !
           CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
           !
           ! and sum to the real space part
           !
           DO nr = 1, nrm
              rr = SQRT(r2(nr)) * alat
              ewaldr = ewaldr + zv(ityp(na)) * zv(ityp(nb)) * &
                       erfc(SQRT(alpha) * rr) / rr
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  ewald = 0.5d0 * e2 * (ewaldg + ewaldr)
  IF ( do_comp_mt ) ewald =  ewald + wg_corr_ewald( omega, ntyp, ngm, zv, strf )
  !
  CALL mp_sum( ewald, intra_bgrp_comm )
  !      CALL mp_sum( ewaldr, intra_bgrp_comm )
  !      CALL mp_sum( ewaldg, intra_bgrp_comm )
  !      WRITE( stdout,'(/5x,"alpha used in ewald term: ",f4.2/
  !     + 5x,"R-space term: ",f12.7,5x,"G-space term: ",f12.7/)')
  !     + alpha, ewaldr, ewaldg
  RETURN
  !
END FUNCTION ewald

