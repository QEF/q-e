!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_tab_beta ( omega, intra_bgrp_comm ) 
  !----------------------------------------------------------------------
  !
  ! Compute interpolation table for beta(G) radial functions
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam, nsp
  USE uspp_data,    ONLY : nqx, dq, tab, tab_d
  USE mp,           ONLY : mp_sum
  USE m_gth,        ONLY : mk_ffnl_gth
  !
  IMPLICIT NONE
  !
  real(DP), intent(in) :: omega
  integer,  intent(in) :: intra_bgrp_comm
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, iq, ir
  REAL(dp) :: qi
  ! q-point grid for interpolation
  REAL(dp) :: pref
  ! the prefactor of the Q functions
  real(DP) ::  vqint, d1
  !
  REAL(dp), allocatable :: aux (:)
  ! work space
  REAL(dp), allocatable :: besr(:)
  ! work space
  !
  ndm = MAXVAL ( upf(:)%kkbeta )
  allocate( aux (ndm) )
  allocate (besr( ndm))
  pref = fpi / sqrt (omega)
  call divide (intra_bgrp_comm, nqx, startq, lastq)
  tab (:,:,:) = 0.d0
  do nt = 1, nsp
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        do iq = startq, lastq
           qi = (iq - 1) * dq
           if ( upf(nt)%is_gth ) then
              CALL mk_ffnl_gth( nt, nb, 1, omega, [ qi ] , tab(iq,nb,nt) )
           else
              call sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
              do ir = 1, upf(nt)%kkbeta
                 aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * rgrid(nt)%r(ir)
              enddo
              call simpson (upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
              tab (iq, nb, nt) = vqint * pref
           end if
        enddo
     enddo
  enddo
  deallocate (besr)
  deallocate (aux)
  !
  call mp_sum(  tab, intra_bgrp_comm )
  !
  ! update GPU memory (taking care of zero-dim allocations)
  !
#if defined __CUDA
  if ( nbetam > 0 ) tab_d=tab
#endif
  !
END SUBROUTINE init_tab_beta
