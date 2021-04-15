!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_tab_atwfc( omega, intra_bgrp_comm)
  !-----------------------------------------------------------------------
  !! This routine computes a table with the radial Fourier transform 
  !! of the atomic wavefunctions.
  !
  USE upf_kinds,    ONLY : DP
  USE atom,         ONLY : rgrid, msh
  USE upf_const,    ONLY : fpi
  USE uspp_data,    ONLY : tab_at, tab_at_d, nqx, dq
  USE uspp_param,   ONLY : nsp, upf
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: omega
  INTEGER,  INTENT(IN) :: intra_bgrp_comm
  !
  INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm
  !
  REAL(DP), ALLOCATABLE :: aux(:), vchi(:)
  REAL(DP) :: vqint, pref, q
  !
  ndm = MAXVAL(msh(1:nsp))
  ALLOCATE( aux(ndm), vchi(ndm) )
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  pref = fpi / SQRT(omega)
  ! needed to normalize atomic wfcs (not a bad idea in general and 
  ! necessary to compute correctly lda+U projections)
  CALL divide( intra_bgrp_comm, nqx, startq, lastq )
  !
  tab_at(:,:,:) = 0.0_DP
  !
  DO nt = 1, nsp
     DO nb = 1, upf(nt)%nwfc
        !
        IF (upf(nt)%oc(nb) >= 0.0_DP) THEN
           l = upf(nt)%lchi (nb)
           !
           DO iq = startq, lastq
              q = dq * (iq - 1)
              CALL sph_bes( msh(nt), rgrid(nt)%r, q, l, aux )
              DO ir = 1, msh(nt)
                 vchi(ir) = upf(nt)%chi(ir,nb) * aux(ir) * rgrid(nt)%r(ir)
              ENDDO
              CALL simpson( msh(nt), vchi, rgrid(nt)%rab, vqint )
              tab_at( iq, nb, nt ) = vqint * pref
           ENDDO
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  CALL mp_sum( tab_at, intra_bgrp_comm )
  !
#if defined __CUDA
  ! update GPU memory (taking care of zero-dim allocations)
  if (SIZE(tab_at)>0) tab_at_d=tab_at
#endif
  !
  DEALLOCATE( aux, vchi )
  !
  RETURN
  !
END SUBROUTINE init_tab_atwfc

