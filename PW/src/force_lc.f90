!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, rho, nl, gstart, gamma_only, vloc, forcelc)
  !----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : tpi
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE esm,       ONLY : esm_force_lc, do_comp_esm, esm_bc
  USE Coul_cut_2D, ONLY : do_cutoff_2D, cutoff_force_lc
  implicit none
  !
  !   first the dummy variables
  !
  integer, intent(in) :: nat, ngm, ngl, gstart, &
                         igtongl (ngm), nl (ngm), ityp (nat)
  ! nat:    number of atoms in the cell
  ! ngm:    number of G vectors
  ! ngl:    number of shells
  ! igtongl correspondence G <-> shell of G
  ! nl:     correspondence fft mesh <-> G vec
  ! ityp:   types of atoms

  logical, intent(in) :: gamma_only

  real(DP), intent(in) :: tau (3, nat), g (3, ngm), vloc (ngl, * ), &
       rho (dfftp%nnr), alat, omega
  ! tau:  coordinates of the atoms
  ! g:    coordinates of G vectors
  ! vloc: local potential
  ! rho:  valence charge
  ! alat: lattice parameter
  ! omega: unit cell volume

  real(DP), intent(out) :: forcelc (3, nat)
  ! the local-potential contribution to forces on atoms

  integer :: ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  complex(DP), allocatable :: aux (:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  allocate (aux(dfftp%nnr))
  !
  aux(:) = CMPLX( rho(:), 0.0_dp, kind=dp )
  !
  CALL fwfft ('Rho', aux, dfftp)
  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
!$omp parallel do private(arg)
  do na = 1, nat
     !
     forcelc (1:3, na) = 0.d0
     ! contribution from G=0 is zero
     do ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        forcelc (1:3, na) = forcelc (1:3, na) + &
                g (1:3, ig) * vloc (igtongl (ig), ityp (na) ) * &
                (sin(arg)*DBLE(aux(nl(ig))) + cos(arg)*AIMAG(aux(nl(ig))) )
     enddo
     !
     forcelc (1:3, na) = fact * forcelc (1:3, na) * omega * tpi / alat
  enddo
!$omp end parallel do
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     !
     CALL esm_force_lc ( aux, forcelc )
  ENDIF
  !
  ! IN 2D calculations: re-add the erf/r contribution to the forces. It was substracted from
  ! vloc (in vloc_of_g) and readded to vltot only (in setlocal)
  IF ( do_cutoff_2D ) THEN
     !     !
      CALL cutoff_force_lc (aux, forcelc )
  ENDIF 
  !
  !
  call mp_sum(  forcelc, intra_bgrp_comm )
  !
  deallocate (aux)
  return
end subroutine force_lc
