!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE add_shift_cc (shift_cc)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY: tpi
  USE atom, ONLY: rgrid
  USE uspp_param, ONLY: upf
  USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau
  USE cell_base, ONLY: alat, omega, tpiba, tpiba2
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY: ngm, gstart, nl, g, gg, ngl, gl, igtongl
  USE ener, ONLY: etxc, vtxc
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho, rho_core, rhog_core
  USE control_flags, ONLY: gamma_only
  USE wavefunctions_module,    ONLY : psic
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE
  !
  !   first the dummy variable
  !

  real(DP) :: shift_cc (nat)
  ! output: the local forces on atoms

  INTEGER :: ig, ir, nt, na
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms


  real(DP), ALLOCATABLE :: vxc (:,:), rhocg (:), shift_(:)
  ! exchange-correlation potential
  ! radial fourier trasform of rho core
  real(DP)  ::  arg, fact
  !
  IF ( any (upf(1:ntyp)%nlcc) ) GOTO 15
  RETURN
  !
15 CONTINUE
   CALL infomsg ('add_shift_cc','BEWARE: shift with CC never tested !!!')
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
  ! recalculate the exchange-correlation potential
  !
  ALLOCATE ( vxc(dfftp%nnr,nspin), shift_(nat) )
  shift_(:) = 0.d0
  !
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  IF (nspin==1) THEN
     DO ir = 1, dfftp%nnr
        psic (ir) = vxc (ir, 1)
     ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     ENDDO
  ENDIF
  DEALLOCATE (vxc)
  CALL fwfft ('Dense', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  ALLOCATE ( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN

        CALL drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, rgrid(nt)%r, &
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        DO na = 1, nat
           IF (nt == ityp (na) ) THEN
              IF (gstart==2)  shift_(na) = omega * rhocg (igtongl (1) ) * &
                                                     conjg(psic (nl (1) ) )
              DO ig = gstart, ngm
                 arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                      + g (3, ig) * tau (3, na) ) * tpi
                 shift_ (na) = shift_( na) + omega * &
                         rhocg (igtongl (ig) ) * conjg(psic (nl (ig) ) ) * &
                         cmplx( cos(arg), -sin(arg),kind=DP) * fact
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
#if defined(__MPI)
  CALL mp_sum( shift_ , intra_pool_comm )
#endif
  shift_cc(:) = shift_cc(:) + shift_(:)
  DEALLOCATE (rhocg, shift_)
  !
  RETURN
END SUBROUTINE add_shift_cc

