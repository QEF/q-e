!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_lc_gpu (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl_d, g_d, rho, nl_d, gstart, gamma_only, vloc_d, forcelc)
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
  USE gbuffers,  ONLY : dev_buf
  implicit none
  !
  !   first the dummy variables
  !
  integer, intent(in) :: nat, ngm, ngl, gstart, &
                         igtongl_d (ngm), nl_d (ngm), ityp (nat)
#if defined(__CUDA)
  attributes(DEVICE) :: igtongl_d, nl_d
#endif
  ! nat:    number of atoms in the cell
  ! ngm:    number of G vectors
  ! ngl:    number of shells
  ! igtongl correspondence G <-> shell of G
  ! nl:     correspondence fft mesh <-> G vec
  ! ityp:   types of atoms

  logical, intent(in) :: gamma_only

  real(DP), intent(in) :: tau (3, nat), g_d (3, ngm), vloc_d (ngl, * ), &
       rho (dfftp%nnr), alat, omega
#if defined(__CUDA)
  attributes(DEVICE) :: g_d, vloc_d
#endif
  ! tau:  coordinates of the atoms
  ! g:    coordinates of G vectors
  ! vloc: local potential
  ! rho:  valence charge
  ! alat: lattice parameter
  ! omega: unit cell volume

  real(DP), intent(out) :: forcelc (3, nat)
  ! the local-potential contribution to forces on atoms

  integer :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  complex(DP), allocatable :: aux (:)
  complex(DP), pointer     :: aux_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d
#endif
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! auxiliary for cuf kernels and error checking
  integer :: ityp_na, ierr
  real(DP):: forcelc_x, forcelc_y, forcelc_z, tau1, tau2, tau3
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  allocate (aux(dfftp%nnr))
  CALL dev_buf%lock_buffer(aux_d, dfftp%nnr, ierr)
  !
  aux(:) = CMPLX( rho(:), 0.0_dp, kind=dp )
  aux_d = aux
  !
  CALL fwfft ('Rho', aux_d, dfftp)
  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nat
     ! This is no longher needed
     !do ipol = 1, 3
     !   forcelc (ipol, na) = 0.d0
     !enddo
     !
     ityp_na = ityp (na)
     tau1 = tau (1, na); tau2 = tau (2, na); tau3 = tau (3, na);
     forcelc_x = 0.d0 ; forcelc_y = 0.d0 ; forcelc_z = 0.d0
     ! contribution from G=0 is zero
     !$cuf kernel do(1)
     do ig = gstart, ngm
        arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 + &
               g_d (3, ig) * tau3 ) * tpi
        !
        ! arg is used as auxiliary variable here
        arg = vloc_d (igtongl_d (ig), ityp_na ) * &
                (sin(arg)*DBLE(aux_d(nl_d(ig))) + cos(arg)*AIMAG(aux_d(nl_d(ig))) )
        !
        forcelc_x = forcelc_x + g_d (1, ig) * arg
        forcelc_y = forcelc_y + g_d (2, ig) * arg
        forcelc_z = forcelc_z + g_d (3, ig) * arg
        !
     enddo
     forcelc (1, na) = forcelc_x; forcelc (2, na) = forcelc_y; forcelc (3, na) = forcelc_z; 
     do ipol = 1, 3
        forcelc (ipol, na) = fact * forcelc (ipol, na) * omega * tpi / alat
     enddo
  enddo
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
  CALL dev_buf%release_buffer(aux_d, ierr)
  return
end subroutine force_lc_gpu
