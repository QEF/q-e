!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
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
SUBROUTINE v_h_from_rho_r( rhotot, nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, nl, nlm, &
                ngm, gg, gstart, alat, omega, ehart, charge, vltot )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from n(r)
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : fpi, e2
  USE cell_base,      ONLY : tpiba2
  USE mp_global,      ONLY : me_pool, intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE fft_base,       ONLY : grid_gather, grid_scatter, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE control_flags,  ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                         nrxx, ngm, gstart, nl(ngm), nlm(ngm)
  !
  REAL (DP), INTENT(IN) :: rhotot(nr1x*nr2x*nr3x), gg(ngm), alat, omega
  !
  REAL (DP), INTENT(OUT) :: vltot(nr1x*nr2x*nr3x), ehart, charge


  REAL (DP),ALLOCATABLE  :: rho(:)
  !
  REAL (DP),ALLOCATABLE  :: v(:)
  !
  ! ... local variables
  !
  REAL (DP)              :: fac
  COMPLEX (DP), ALLOCATABLE :: aux(:), aux1(:)
  INTEGER                     :: ir, is, ig
  !

  ALLOCATE( rho( nrxx ) )
  ALLOCATE( v( nrxx ) )
  
#ifdef __PARA
  CALL grid_scatter(rhotot, rho)
  CALL grid_scatter(vltot,    v)
#else
  rho( : ) = rhotot( : )
  v( : ) = vltot( : )
#endif
  
  !
  ALLOCATE( aux( nrxx ), aux1( ngm ) )
  !
  ! ... copy total rho in aux
  !
  aux(:) = CMPLX ( rho(:), 0.0_dp, KIND=dp )
  !
  ! ... bring rho (aux) to G space
  ! 
  CALL fwfft ('Dense', aux, dfftp)

  !
  charge = 0.D0
  !

  IF ( gstart == 2 ) charge = omega * DBLE (aux(nl(1)))
  !
  CALL mp_sum( charge, intra_pool_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  ehart   = 0.D0
  aux1(:) = 0.D0
  !
  DO ig = gstart, ngm
     !
     fac = 1.D0 / gg(ig)
     !
     ehart = ehart + ( DBLE(aux(nl(ig)))**2 + AIMAG(aux(nl(ig)))**2 ) * fac
     !
     aux1(ig) = aux(nl(ig)) * fac
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
  aux(:) = 0.D0
  !
  DO ig = 1, ngm
     !
     aux(nl(ig)) = aux1(ig)
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     DO ig = 1, ngm
        !
        aux(nlm(ig)) =  CONJG( aux1(ig) )
        !
     END DO
     !
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL invfft ('Dense', aux, dfftp)
  !
  ! ... add hartree potential to the xc potential
  !
  v(:) = v(:) + DBLE ( aux(:) )
  !
  DEALLOCATE( aux, aux1 )
  !
#ifdef __PARA
  vltot(:)=0
  CALL grid_gather(v,    vltot)
  CALL mp_sum( vltot, intra_pool_comm )
#else
  vltot ( : ) = v( : )
#endif

  DEALLOCATE( rho )
  DEALLOCATE( v )


  RETURN
  !
END SUBROUTINE v_h_from_rho_r


