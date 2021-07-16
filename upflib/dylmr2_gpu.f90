!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dylmr2_gpu( nylm, ngy, g_d, gg_d, dylm_d, ipol )
  !-----------------------------------------------------------------------
  !! Compute \partial Y_lm(G) \over \partial (G)_ipol
  !! using simple numerical derivation (SdG).
  !! The spherical harmonics are calculated in ylmr2.
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
  REAL(DP), INTENT(IN) :: g_d(3,ngy)
  !! input: the coordinates of g vectors
  REAL(DP), INTENT(IN) :: gg_d(ngy)
  !! input: the moduli of g vectors
  REAL(DP), INTENT(OUT) :: dylm_d(ngy,nylm)
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
  REAL(DP), PARAMETER :: delta = 1.E-6_DP
  REAL(DP), ALLOCATABLE :: dg_d(:), dgi_d(:), gx_d(:,:)
  REAL(DP), ALLOCATABLE :: ggx_d(:), ylmaux_d(:,:)
  ! dg is the finite increment for numerical derivation:
  ! dg = delta |G| = delta * sqrt(gg)
  ! dgi= 1 /(delta * sqrt(gg))
  ! gx = g +/- dg
  ! ggx = gx^2
  !
#if defined(__CUDA)
  attributes(DEVICE) :: g_d, gg_d, dylm_d, gx_d, ggx_d, dg_d, &
                        dgi_d, ylmaux_d
#endif  
  !
  
  !
  ALLOCATE( gx_d(3,ngy), ggx_d(ngy), dg_d(ngy) )
  ALLOCATE( dgi_d(ngy), ylmaux_d(ngy,nylm) )

  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, ngy
     dg_d(ig) = delta * SQRT(gg_d(ig) )
     IF (gg_d(ig) > 1.E-9_DP) THEN
        dgi_d(ig) = 1._DP / dg_d(ig)
     ELSE
        dgi_d(ig) = 0._DP
     ENDIF
  ENDDO
  !
  IF ( ipol==1 ) THEN
    apol = 2 ; bpol = 3
  ELSEIF ( ipol==2 ) THEN
    apol = 1 ; bpol = 3
  ELSEIF ( ipol==3 ) THEN
    apol = 1 ; bpol = 2
  ENDIF
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, ngy
     gx_d(apol,ig) = g_d(apol,ig)
     gx_d(bpol,ig) = g_d(bpol,ig)
     gx_d(ipol,ig) = g_d(ipol,ig) + dg_d(ig)
     ggx_d(ig) = gx_d(1,ig) * gx_d(1,ig) + &
                 gx_d(2,ig) * gx_d(2,ig) + &
                 gx_d(3,ig) * gx_d(3,ig)
  ENDDO
  !
  CALL ylmr2_gpu( nylm, ngy, gx_d, ggx_d, dylm_d )
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, ngy
     gx_d(ipol,ig) = g_d(ipol,ig) - dg_d(ig)
     ggx_d(ig) = gx_d(1,ig) * gx_d(1,ig) + &
                 gx_d(2,ig) * gx_d(2,ig) + &
                 gx_d(3,ig) * gx_d(3,ig)
  ENDDO
  !
  CALL ylmr2_gpu( nylm, ngy, gx_d, ggx_d, ylmaux_d )
  !
  !$cuf kernel do (2) <<<*,*>>>
  DO lm = 1, nylm
     DO ig = 1, ngy
        dylm_d(ig,lm) = (dylm_d(ig,lm)-ylmaux_d(ig,lm)) * 0.5_DP * dgi_d(ig)
     ENDDO
  ENDDO
  !
  DEALLOCATE( gx_d, ggx_d, dg_d, dgi_d, ylmaux_d )
  !
  RETURN
  !
END SUBROUTINE dylmr2_gpu

