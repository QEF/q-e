!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE add_shift_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, nrxx, g, rho, nl, nspin, gstart, gamma_only, vloc, shift_lc)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY : fwfft
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE
  !
  !   first the dummy variables
  !
  INTEGER :: nat, ngm, nrxx, nspin, &
       ngl, gstart, igtongl (ngm), nl (ngm), ityp (nat)
  ! input: the number of atoms in the cell
  ! input: the number of G vectors
  ! input: number of spin polarizations
  ! input: the number of shells
  ! input: correspondence G <-> shell of G
  ! input: the correspondence fft mesh <-> G vec
  ! input: the types of atoms

  LOGICAL :: gamma_only

  real(DP) :: tau (3, nat), g (3, ngm), vloc (ngl, * ), &
       rho (nrxx, nspin), alat, omega
  ! input: the coordinates of the atoms
  ! input: the coordinates of G vectors
  ! input: the local potential
  ! input: the valence charge
  ! input: the length measure
  ! input: the volume of the cell

  real(DP) :: shift_lc ( nat)
  ! output: the local forces on atoms

  INTEGER :: ig, na
  ! counter on G vectors
  ! counter on atoms

  real(DP), ALLOCATABLE :: shift_(:)
  complex(DP), ALLOCATABLE :: aux (:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  ALLOCATE (aux(nrxx), shift_(nat) )
  shift_(:) = 0.d0

  IF (nspin==2) THEN
     aux(:) = CMPLX ( rho(:,1)+rho(:,2), 0.0_dp, KIND=dp )
  ELSE
     aux(:) = CMPLX ( rho(:,1), 0.0_dp, KIND=dp )
  END IF
  CALL fwfft ('Dense', aux, dfftp)
  !
  !    aux contains now  n(G)
  !
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  DO na = 1, nat
     ! contribution from G=0 is not zero but should be counted only once
     IF (gstart==2) shift_(na) = vloc(igtongl(1),ityp(na)) * DBLE (aux(nl(1))) / fact
     DO ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        shift_ ( na) = shift_ (na) + &
                vloc (igtongl (ig), ityp (na) ) * &
                (cos (arg) * DBLE (aux(nl(ig))) - sin (arg) * AIMAG (aux(nl(ig))) )
     ENDDO
     shift_ (na) = fact * shift_ (na) * omega
  ENDDO
#if defined(__MPI)
  CALL mp_sum(  shift_, intra_pool_comm )
#endif

  shift_lc(:) = shift_lc(:) + shift_(:)

  DEALLOCATE (aux,shift_)


  RETURN
END SUBROUTINE add_shift_lc
