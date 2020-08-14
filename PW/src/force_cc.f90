!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE force_cc( forcecc )
  !----------------------------------------------------------------------
  !! Calculates the NLCC contribution to the force.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions,        ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcecc(3,nat)
  !! output: the NLCC forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms
  REAL(DP), ALLOCATABLE :: vxc(:,:), rhocg(:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  REAL(DP) ::  arg, fact
  !
  !
  forcecc(:,:) = 0.d0
  !
  IF ( ANY( upf(1:ntyp)%nlcc ) ) GOTO 15
  RETURN
  !
15 CONTINUE
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
  ! ... recalculate the exchange-correlation potential
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
  !
  psic = (0.0_DP,0.0_DP)
  IF (nspin == 1 .OR. nspin == 4) THEN
     DO ir = 1, dfftp%nnr
        psic(ir) = vxc (ir,1)
     ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     ENDDO
  ENDIF
  !
  DEALLOCATE( vxc )
  !
  CALL fwfft( 'Rho', psic, dfftp )
  !
  ! ... psic contains now Vxc(G)
  !
  ALLOCATE( rhocg(ngl) )
  !
  ! ... core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN
        !
        CALL drhoc( ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r, &
                    rgrid(nt)%rab, upf(nt)%rho_atc, rhocg )
!$omp parallel do private(arg)
        DO na = 1, nat
           IF (nt == ityp (na) ) THEN
              DO ig = gstart, ngm
                 arg = (g(1,ig) * tau(1,na) + g (2, ig) * tau (2, na) &
                      + g(3,ig) * tau(3,na) ) * tpi
                 forcecc (1:3, na) = forcecc(1:3, na) + tpiba * omega * &
                         rhocg(igtongl(ig)) * CONJG(psic(dfftp%nl(ig) ) ) * &
                         CMPLX( SIN(arg), COS(arg), KIND=DP) * g(1:3,ig) * fact
              ENDDO
           ENDIF
        ENDDO
!$omp end parallel do
     ENDIF
  ENDDO
  !
  CALL mp_sum( forcecc, intra_bgrp_comm )
  !
  DEALLOCATE( rhocg )
  !
  RETURN
  !
END SUBROUTINE force_cc
