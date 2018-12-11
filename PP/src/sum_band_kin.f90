!-----------------------------------------------------------------------
! Program written by Yang Jiao,  Oct 2016, GPL, No warranties.
! 
! This subroutine is adapted from SUBROUTINE sum_band from PW/src/sum_band.f90
! ----------------------------------------------------------------------
SUBROUTINE sum_band_kin(kin_r)
  !----------------------------------------------------------------------------
  !
  ! ... Calculates the Kohn-Sham kinetic-energy density t_{KS}
  ! ... adapted from the original sum_band subroutine
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft, fft_interpolate
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE spin_orb,             ONLY : lspinorb, domag, fcoef
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : dft_is_meta
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: kin_r(dfftp%nnr,nspin)
  !
  ! ... local variables
  !
  INTEGER :: ir,   &! counter on 3D r points
             is,   &! counter on spin polarizations
             ig,   &! counter on g vectors
             ibnd, &! counter on bands
             ik,   &! counter on k points
             npol_,&! auxiliary dimension for noncolin case
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  COMPLEX(DP), ALLOCATABLE :: kin_g(:,:)          ! the kinetic energy density in G-space
  !
  !
  CALL start_clock( 'sum_band_kin' )
  !
  ALLOCATE(kin_g(ngm,nspin))
  kin_r=0.D0
  kin_g=(0.D0,0.D0)
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL weights ( )
  !
  IF (one_atom_occupations) CALL new_evc()
  !
  ! ... btype, used in diagonalization, is set here: a band is considered empty
  ! ... and computed with low accuracy only when its occupation is < 0.01, and
  ! ... only if option diago_full_acc is false; otherwise, use full accuracy
  !
!  ALLOCATE( btype( nbnd, nkstot ) )
!  btype(:,:) = 1
!  IF ( .NOT. diago_full_acc ) THEN
!     !
!     FORALL( ik = 1:nks, wk(ik) > 0.D0 )
!        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
!     END FORALL
!     !
!  END IF
  !
  ! ... Needed for LDA+U: compute occupations of Hubbard states
  !
!  IF (lda_plus_u) THEN
!     IF(noncolin) THEN
!        CALL new_ns_nc(rho%ns_nc)
!     ELSE
!        CALL new_ns(rho%ns)
!     ENDIF
!  ENDIF
  !
  ! ... for band parallelization: set band computed by this processor
  !
  call divide ( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type (nkb,nbnd, becp,intra_bgrp_comm)
  ALLOCATE (kplusg(npwx))
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  ! ... The band energy contribution eband is computed together with the charge
  !
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_kin_gamma(kin_r,kin_g)
     !
  ELSE
     !
     CALL sum_band_kin_k(kin_r,kin_g)
     !
  END IF
  !
  DEALLOCATE (kplusg)
  IF ( okvan ) CALL deallocate_bec_type ( becp )

  !
  ! ... kin_r: sum over bands, k-points, bring to G-space, symmetrize,
  ! ... synchronize with kin_g
  !
!  IF ( dft_is_meta() .OR. lxdm) THEN
     !
     CALL mp_sum( kin_r, inter_pool_comm )
     CALL mp_sum( kin_r, inter_bgrp_comm )
     DO is = 1, nspin
        psic(1:dffts%nnr) = kin_r(1:dffts%nnr,is)
        psic(dffts%nnr+1:) = 0.0_dp
        CALL fwfft ('Rho', psic, dffts)
        kin_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
        kin_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
     END DO
     !
     IF (.NOT. gamma_only) CALL sym_rho( nspin, kin_g )
     !
     DO is = 1, nspin
        psic(:) = ( 0.D0, 0.D0 )
        psic(dfftp%nl(:)) = kin_g(:,is)
        IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( kin_g(:,is) )
        CALL invfft ('Rho', psic, dfftp)
        kin_r(:,is) = psic(:)
     END DO
     !
!  END IF
  !
  CALL stop_clock( 'sum_band_kin' )
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_kin_gamma(kin_r,kin_g)
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       USE becmod,        ONLY : becp
       USE mp_bands,      ONLY : me_bgrp
       USE mp,            ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       !
       IMPLICIT NONE
       REAL(DP), INTENT(INOUT)    :: kin_r(dfftp%nnr,nspin)
       COMPLEX(DP), INTENT(INOUT) :: kin_g(ngm,nspin)
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, ioff_tg, nxyp, incr, v_siz, j, ir3
       LOGICAL :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = .False.
       !
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             !
          !   IF (dft_is_meta() .OR. lxdm) THEN
                DO j=1,3
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   kplusg (1:npw) = (xk(j,ik)+g(j,1:npw)) * tpiba

                   IF ( ibnd < ibnd_end ) THEN
                      ! ... two ffts at the same time
                      psic(dffts%nl (1:npw))=CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                            ( evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   ELSE
                      psic(dffts%nl(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) )
                   END IF
                   !
                   CALL invfft ('Wave', psic, dffts)
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, dffts%nnr
                      kin_r(ir,current_spin) = &
                                           kin_r(ir,current_spin) + &
                                           w1 *  DBLE( psic(ir) )**2 + &
                                           w2 * AIMAG( psic(ir) )**2
                   END DO
                   !
                END DO
          !   END IF
             !
             !
          END DO
          !
          !
       END DO k_loop
       !
       !
       RETURN
       !
     END SUBROUTINE sum_band_kin_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_kin_k(kin_r,kin_g)
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE mp_bands,     ONLY : me_bgrp
       USE mp,           ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       !
       IMPLICIT NONE
       REAL(DP), INTENT(INOUT)    :: kin_r(dfftp%nnr,nspin)
       COMPLEX(DP), INTENT(INOUT) :: kin_g(ngm,nspin)
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz, j, ir3
       LOGICAL  :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp
       !
       ! chunking parameters
       INTEGER, PARAMETER :: blocksize = 256
       INTEGER :: numblock
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = .False.
       !
       incr = 1
       !
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             ! ... the sum of eband and demet is the integral for e < ef of
             ! ... e n(e) which reduces for degauss=0 to the sum of the
             ! ... eigenvalues
             w1 = wg(ibnd,ik) / omega
             !
             IF (.NOT. noncolin) THEN
                !
!                IF (dft_is_meta() .OR. lxdm) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
                      psic(dffts%nl(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      !
                      CALL invfft ('Wave', psic, dffts)
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho(kin_r(:,current_spin), dffts%nnr, w1, psic)
                   END DO
!                END IF
                !
             END IF
             !
          END DO
          !
       END DO k_loop
       !
       RETURN
       !
     END SUBROUTINE sum_band_kin_k
     !
     !
     SUBROUTINE get_rho(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * ( DBLE( psic_loc(ir) )**2 + &
                                   AIMAG( psic_loc(ir) )**2 )
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho

     SUBROUTINE get_rho_gamma(rho_loc, nrxxs_loc, w1_loc, w2_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * DBLE( psic_loc(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc(ir) )**2
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho_gamma


     SUBROUTINE get_rho_domag(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(:, :)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir,2) = rho_loc(ir,2) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))* DBLE(psic_loc(ir,2)) + &
                          AIMAG(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)))
 
           rho_loc(ir,3) = rho_loc(ir,3) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)) - &
                           DBLE(psic_loc(ir,2))*AIMAG(psic_loc(ir,1)))

           rho_loc(ir,4) = rho_loc(ir,4) + w1_loc* &
                          (DBLE(psic_loc(ir,1))**2+AIMAG(psic_loc(ir,1))**2 &
                          -DBLE(psic_loc(ir,2))**2-AIMAG(psic_loc(ir,2))**2)
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho_domag

END SUBROUTINE sum_band_kin

