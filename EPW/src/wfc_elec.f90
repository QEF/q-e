!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!-----------------------------------------------------------------------
MODULE polaron
    USE kinds,  ONLY : dp
    IMPLICIT NONE
    ! Data block, try to keep it to minimal
    COMPLEX(KIND = DP), ALLOCATABLE :: epfall(:, :, :, :, :)
    !! el-ph element for all local k and all q
    !! epfall need to be filled in ephwann_shuffle
    COMPLEX(KIND = DP), ALLOCATABLE :: ufall(:, :, :)
    !! el-ph element for all local k and all q
    !! epfall need to be filled in ephwann_shuffle
    COMPLEX(KIND = DP), ALLOCATABLE :: Hamil(:, :)
    !! Hamil need to be passed to h_psi because the parameter space is fixed
    !! to meet the requirement of Davidson diagonalization. Ugly but workable.
    COMPLEX(KIND = DP), ALLOCATABLE :: eigVec(:, :)
    !! polaron eigenvector
    REAL(KIND = DP),    ALLOCATABLE :: etf_all(:, :)
    !! Gather all the eigenvalues
    INTEGER,            ALLOCATABLE :: ikq_all(:, :), kpg_map(:)
    !
    PUBLIC  :: wfc_elec, interp_plrn_wf, interp_plrn_bq, plot_plrn_wf
CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
        USE io_global,     ONLY : stdout, ionode

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: nrr_k, dims, ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
        REAL(DP), INTENT (IN) :: irvec_r(3, nrr_k)
        COMPLEX(DP), ALLOCATABLE :: eigvec_wan(:, :)
        INTEGER :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE interp_plrn_bq(nrr_q, ndegen_q, irvec_q)
        USE epwcom,        ONLY : nkf1, nkf2, nkf3, nbndsub
        USE elph2, only : xqf, wf, nqtotf
        USE modes,         ONLY : nmodes
        USE constants_epw, only : eps8, czero, one, two, twopi, ci
        USE ions_base,     ONLY : nat, amass, ityp, tau
        USE wan2bloch, only : dynwan2bloch

        IMPLICIT NONE
        INTEGER, INTENT (IN) :: nrr_q, ndegen_q(:,:,:) ! ! Added for polaron calculations by Chao Lian.
        INTEGER, INTENT (IN) :: irvec_q(3, nrr_q)

        INTEGER :: dtau_file
        INTEGER :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nat_p

        INTEGER :: iq, inu, ierr, imu, na, iatm, idir
        INTEGER :: icount, ix, iy, iz, bmat_file
        COMPLEX(DP) :: ctemp, shift(3)

        COMPLEX(DP), ALLOCATABLE :: uf(:, :), Bmat(:,:)
        COMPLEX(DP),  ALLOCATABLE :: dtau(:, :)
        REAL(DP),  ALLOCATABLE :: w2(:)
        REAL(KIND=dp) :: xxq(3)
        COMPLEX(KIND=dp) :: expTable(3)

    END SUBROUTINE

    SUBROUTINE wfc_elec (nrr_k, ndegen_k, irvec_r, dims)
        !
        ! Self consistency calculation of polaron wavefunction.
        ! Rewritten by Chao Lian based on the implementation by Danny Sio.
        !
        USE modes,         ONLY : nmodes
        USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero
        USE constants_epw, ONLY : czero, cone, pi, ci, twopi, eps6, eps8, eps5
        USE epwcom,        ONLY : num_cbands, polaron_type, sigma_plrn, full_diagon_plrn
        USE epwcom,        ONLY : r01, r02, r03, nPlrn, conv_thr_polaron, cb_shift
        USE epwcom,        ONLY : mixing_Plrn, init_plrn_wf, niterPlrn
        USE epwcom,        ONLY : nkf1, nkf2, nkf3, nbndsub
        USE io_global,     ONLY : stdout, ionode, meta_ionode_id
        USE elph2,         ONLY : etf, ibndmin, ibndmax, nbndfst
        USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
        USE elph2,         ONLY : xkf, xqf, wf, xkq, chw
        USE mp_global,     ONLY : inter_pool_comm
        USE mp_world,      ONLY : world_comm
        USE cell_base,     ONLY : bg
        USE mp,            ONLY : mp_sum, mp_bcast
        USE poolgathering, ONLY : poolgather2
        USE test_tools,    ONLY : para_write
        USE wan2bloch,     ONLY : hamwan2bloch
        USE ions_base,     ONLY : nat

        IMPLICIT NONE

        ! local variables
        LOGICAL :: debug
        INTEGER :: inu, iq, ik, ikk, jk, ibnd, jbnd, ikq, ik_global, iplrn, ierr
        INTEGER :: iter, icount, ix, iy, iz, start_mode, ik_bm, idos, iatm

        INTEGER, INTENT (IN) :: nrr_k, dims, ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
        REAL(DP), INTENT (IN) :: irvec_r(3, nrr_k)

        COMPLEX(DP),  ALLOCATABLE :: Bmat(:,:), Bmat_save(:,:)
        COMPLEX(DP),  ALLOCATABLE :: eigvec_wan(:, :), dtau(:, :)
        REAL(DP),     ALLOCATABLE :: rmat_tmp(:, :)

        COMPLEX(KIND=dp) :: cufkk ( nbndsub, nbndsub ), cfac(nrr_k, dims, dims)
        !! Rotation matrix, fine mesh, points k

        REAL(dp):: estmteRt(nPlrn),  eigVal(nPlrn), esterr

        REAL(KIND=dp) :: qcart(3), r0(3), xxk(3), xxq(3), prefac, norm
        REAL(KIND=dp) :: ef

        INTEGER :: band_pos, iqpg, ikpg, ikGamma, iqGamma
        INTEGER :: nkf1_p, nkf2_p, nkf3_p, nbndsub_p, nPlrn_p, nktotf_p

        REAL(DP) :: eb
        REAL(DP) :: xkf_all(3, nktotf), xkf_all_tmp(3, nktotf*2)
        REAL(DP) :: EPlrnTot, EPlrnElec, EPlrnPhon
        REAL(DP) :: disK, disK_t, shift(3)

        COMPLEX(DP) :: ctemp
        REAL(DP)    :: rtemp
        INTEGER     :: itemp, jtemp
        INTEGER     :: dos_file, wan_func_file, bloch_func_file, bmat_file, dtau_file
        !LOGICAL     :: SCF_run

    END SUBROUTINE
    SUBROUTINE plot_plrn_wf()
    END SUBROUTINE
END MODULE
