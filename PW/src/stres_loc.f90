! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_loc( sigmaloc )
  !----------------------------------------------------------------------
  !! Calculate the local term of the stress.
  !
  USE kinds,                ONLY : DP
  USE vloc_mod,             ONLY : dvloc_of_g
  USE atom,                 ONLY : msh, rgrid
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : omega, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE scf,                  ONLY : rho
  USE vlocal,               ONLY : strf, vloc
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE esm,                  ONLY : do_comp_esm, esm_bc
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_stres_evloc, &
                                   cutoff_stres_sigmaloc 
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmaloc(3,3)
  REAL(DP), ALLOCATABLE :: dvloc(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP) :: evloc, fact
  INTEGER :: ng, nt, l, m
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components
  LOGICAL :: modified_coulomb
  REAL(DP) :: sigma11, sigma21, sigma22, spart, &
              sigma31, sigma32, sigma33
  !
  ALLOCATE( dvloc(ngl), rhog(dfftp%nnr,1) )
  sigmaloc(:,:) = 0.d0
  !
  !$acc data create(rhog)
  CALL rho_r2g( dfftp, rho%of_r(:,1), rhog )
  !
  !$acc data copyin(vloc,strf,gl) present(igtongl) create(dvloc)
  !
  modified_coulomb = do_cutoff_2D .OR. (do_comp_esm .and. ( esm_bc .ne. 'pbc' ))
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
  evloc = 0.0d0
  !
  IF (gstart==2) THEN
    !$acc parallel loop reduction(+:evloc)
    DO nt = 1, ntyp
       evloc = evloc + rhog(1,1)*strf(1,nt)*vloc(igtongl(1),nt)
    ENDDO
  ENDIF
  !
  !$acc parallel loop collapse(2) reduction(+:evloc)
  DO nt = 1, ntyp
    DO ng = gstart, ngm
      evloc = evloc + DBLE(CONJG(rhog(ng,1)) * strf(ng,nt)) &
                      * vloc(igtongl(ng),nt) * fact
    ENDDO
  ENDDO
  !
  ! ... 2D: add contribution from cutoff long-range part of Vloc
  IF (do_cutoff_2D)  CALL cutoff_stres_evloc( gamma_only, rhog(:,1), strf, evloc )
  !
  !      WRITE( 6,*) ' evloc ', evloc, evloc*omega   ! DEBUG
  !
  DO nt = 1, ntyp
     !
     CALL dvloc_of_g( nt, ngl, gl, tpiba2, modified_coulomb, omega, dvloc )
     !
     sigma11 = 0._DP ; sigma21 = 0._DP ; sigma22 = 0._DP
     sigma31 = 0._DP ; sigma32 = 0._DP ; sigma33 = 0._DP
     !
     !$acc parallel loop reduction(+:sigma11,sigma21,sigma22,sigma31,sigma32,&
     !$acc                           sigma33)
     DO ng = 1, ngm
       spart = DBLE(CONJG(rhog(ng,1)) * strf(ng,nt)) * 2.0_DP * &
               dvloc(igtongl(ng))
       sigma11 = sigma11 + spart * g(1,ng) * g(1,ng)
       sigma21 = sigma21 + spart * g(2,ng) * g(1,ng)
       sigma22 = sigma22 + spart * g(2,ng) * g(2,ng)
       sigma31 = sigma31 + spart * g(3,ng) * g(1,ng)
       sigma32 = sigma32 + spart * g(3,ng) * g(2,ng)
       sigma33 = sigma33 + spart * g(3,ng) * g(3,ng)
     ENDDO
     !
     sigmaloc(1,1) = sigmaloc(1,1) + sigma11 * fact * tpiba2
     sigmaloc(2,1) = sigmaloc(2,1) + sigma21 * fact * tpiba2
     sigmaloc(2,2) = sigmaloc(2,2) + sigma22 * fact * tpiba2
     sigmaloc(3,1) = sigmaloc(3,1) + sigma31 * fact * tpiba2
     sigmaloc(3,2) = sigmaloc(3,2) + sigma32 * fact * tpiba2
     sigmaloc(3,3) = sigmaloc(3,3) + sigma33 * fact * tpiba2  
     !
  ENDDO
  !
  ! ... 2D: re-add LR Vloc to sigma here
  IF (do_cutoff_2D)  CALL cutoff_stres_sigmaloc( gamma_only, rhog(:,1), strf, sigmaloc )
  !
  !$acc end data
  !$acc end data
  !
  DO l = 1, 3
     sigmaloc(l,l) = sigmaloc(l,l) + evloc
     DO m = 1, l-1
        sigmaloc(m,l) = sigmaloc(l,m)
     ENDDO
  ENDDO
  !
  CALL mp_sum( sigmaloc, intra_bgrp_comm )
  !
  DEALLOCATE( dvloc, rhog )
  !
  RETURN
  !
END SUBROUTINE stres_loc
