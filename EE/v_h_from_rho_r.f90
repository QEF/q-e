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
!----------------------------------------------------------------------------
SUBROUTINE v_h_from_rho_r( rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm, &
                ngm, gg, gstart, alat, omega, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from n(r)
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : fpi, e2
  USE cell_base,      ONLY : tpiba2
  USE mp_global,      ONLY : me_pool, intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE control_flags,  ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                         nrxx, ngm, gstart, nl(ngm), nlm(ngm)
  !
  REAL (DP), INTENT(IN) :: rho(nrxx), gg(ngm), alat, omega
  !
  REAL (DP), INTENT(OUT) :: v(nrxx), ehart, charge
  !
  ! ... local variables
  !
  REAL (DP)              :: fac
  REAL (DP), ALLOCATABLE :: aux(:,:), aux1(:,:)
  INTEGER                     :: ir, is, ig
  !
  !
  ALLOCATE( aux( 2, nrxx ), aux1( 2, ngm ) )
  !
  ! ... copy total rho in aux
  !
  aux(2,:) = 0.D0
  aux(1,:) = rho(:)
  !
  ! ... bring rho (aux) to G space
  ! 
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  !
  charge = 0.D0
  !
  IF ( gstart == 2 ) charge = omega * aux(1,nl(1))
  !
  CALL mp_sum( charge, intra_pool_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  ehart     = 0.D0
  aux1(:,:) = 0.D0
  !
  DO ig = gstart, ngm
     !
     fac = 1.D0 / gg(ig)
     !
     ehart = ehart + ( aux(1,nl(ig))**2 + aux(2,nl(ig))**2 ) * fac
     !
     aux1(1,ig) = aux(1,nl(ig)) * fac
     aux1(2,ig) = aux(2,nl(ig)) * fac
     !
  ENDDO
  !
  fac = e2 * fpi / tpiba2
  !
  ehart = ehart * fac
  !
  aux1 = aux1 * fac
  !
  IF (gamma_only) THEN
     !
     ehart = ehart * omega
     !
  ELSE
     !
     ehart = ehart * 0.5D0 * omega
     !
  END IF
  !
  CALL mp_sum( ehart, intra_pool_comm )
  !
  aux(:,:) = 0.D0
  !
  DO ig = 1, ngm
     !
     aux(1,nl(ig)) = aux1(1,ig)
     aux(2,nl(ig)) = aux1(2,ig)
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     DO ig = 1, ngm
        !
        aux(1,nlm(ig)) =   aux1(1,ig)
        aux(2,nlm(ig)) = - aux1(2,ig)
        !
     END DO
     !
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  ! ... add hartree potential to the xc potential
  !
  v(:) = v(:) + aux(1,:)
  !
  DEALLOCATE( aux, aux1 )
  !
  RETURN
  !
END SUBROUTINE v_h_from_rho_r


