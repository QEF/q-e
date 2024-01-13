!
! Copyright (C) 2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dylmr2( nylm, ngy, g, gg, dylm, ipol )
  !-----------------------------------------------------------------------
  !! Compute \partial Y_lm(G) \over \partial (G)_ipol
  !! using simple numerical derivation (SdG).
  !! The spherical harmonics are calculated in ylmr2.
  !! ACC-enabled version
  !
  USE upf_kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nylm
  ! input: number of spherical harmonics
  INTEGER, INTENT(IN) :: ngy
  ! input: the number of g vectors to compute
  INTEGER, INTENT(IN) :: ipol
  ! input: desired polarization
  REAL(DP), INTENT(IN) :: g(3,ngy)
  !! input: the coordinates of g vectors
  REAL(DP), INTENT(IN) :: gg(ngy)
  !! input: the moduli of g vectors
  REAL(DP), INTENT(OUT) :: dylm(ngy,nylm)
  !! output: the spherical harmonics derivatives
  !
  ! ... local variables
  !
  INTEGER :: ig, lm, i
  ! counter on g vectors
  ! counter on l,m component
  !
  INTEGER :: apol, bpol
  !
  REAL(DP), PARAMETER :: delta = 1.e-6_dp, eps = 1.e-9_dp
  REAL(DP), ALLOCATABLE :: dg(:), gx(:,:), ggx(:), ylmaux(:,:)
  ! dg is the finite increment for numerical derivation:
  ! dg = delta |G| = delta * sqrt(gg) (later overwritten by 1/dg)
  ! gx = g +/- dg
  ! ggx = gx^2
  !
  IF ( ipol==1 ) THEN
    apol = 2 ; bpol = 3
  ELSEIF ( ipol==2 ) THEN
    apol = 1 ; bpol = 3
  ELSEIF ( ipol==3 ) THEN
    apol = 1 ; bpol = 2
  ENDIF
  !
  ALLOCATE( gx(3,ngy), ggx(ngy), dg(ngy), ylmaux(ngy,nylm) )
  !$acc data create(gx, ggx, dg, ylmaux) &
  !$acc      present_or_copyin(g, gg) present_or_copyout(dylm)
  !
  !$acc parallel loop
  DO ig = 1, ngy
     dg(ig) = delta * SQRT(gg(ig) )
  ENDDO
  !$acc parallel loop
  DO ig = 1, ngy
     gx(apol,ig) = g(apol,ig)
     gx(bpol,ig) = g(bpol,ig)
     gx(ipol,ig) = g(ipol,ig) + dg(ig)
     ggx(ig) = gx(1,ig) * gx(1,ig) + &
               gx(2,ig) * gx(2,ig) + &
               gx(3,ig) * gx(3,ig)
  ENDDO
  !
  CALL ylmr2( nylm, ngy, gx, ggx, dylm )
  !
  !$acc parallel loop
  DO ig = 1, ngy
     gx(ipol,ig) = g(ipol,ig) - dg(ig)
     ggx(ig) = gx(1,ig) * gx(1,ig) + &
               gx(2,ig) * gx(2,ig) + &
               gx(3,ig) * gx(3,ig)
  ENDDO
  !
  CALL ylmr2( nylm, ngy, gx, ggx, ylmaux )
  !
  !$acc parallel loop
  DO ig = 1, ngy
     IF (gg(ig) > eps) THEN
        dg(ig) = 1.0_dp / dg(ig)
     ELSE
        dg(ig) = 0.0_dp
     ENDIF
  ENDDO
  !$acc parallel loop collapse(2)
  DO lm = 1, nylm
     DO ig = 1, ngy
        dylm(ig,lm) = (dylm(ig,lm)-ylmaux(ig,lm)) * 0.5_dp * dg(ig)
     ENDDO
  ENDDO
  !
  !$acc end data
  DEALLOCATE( gx, ggx, dg, ylmaux )
  !
  RETURN
  !
END SUBROUTINE dylmr2
