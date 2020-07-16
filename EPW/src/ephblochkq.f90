!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
! Feliciano Giustino
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano
! Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!-----------------------------------------------------------------------
MODULE ephblochkq
    PUBLIC :: gkg, phonon_eigvector, interpol_a_k, interpol_bq, get_cfac, compute_a_re
CONTAINS
    SUBROUTINE gkq ( iq, nrr_k, nrr_q, nrr_g, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g, &
            w2, uf, epmatwef, irvec_r, dims, dims2)
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
            et, xk, ef,  nelec
        USE cell_base,     ONLY : at, bg, omega, alat
        USE start_k,       ONLY : nk1, nk2, nk3
        USE ions_base,     ONLY : nat, amass, ityp, tau
        USE phcom,         ONLY : nq1, nq2, nq3
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange, &
            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
            etf_mem, scr_typ, &
            elecselfen, phonselfen, nest_fn, a2f, specfun_ph, &
            vme, eig_read, ephwrite, nkf1, nkf2, nkf3, &
            efermi_read, fermi_energy, specfun_el, band_plot, &
            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk, &
            plselfen, specfun_pl, wfcelec
        USE noncollin_module, ONLY : noncolin
        USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
        USE io_files,      ONLY : prefix, diropn
        USE io_global,     ONLY : stdout, ionode
        USE elph2,         ONLY : cu, cuq, lwin, lwinq,&
            chw, chw_ks, cvmew, cdmew, rdw, &
            epmatwp, epmatq, etf, etf_k, etf_ks, xqf, xkf, &
            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
            sigmai_all, sigmai_mode, gamma_all, epsi, zstar, &
            efnew, ifc, sigmar_all, zi_all, nkqtotf, eps_rpa, &
            g2_4, wf, nbndskip
#if defined(__NAG)
        USE f90_unix_io,   ONLY : FLUSH
#endif
        USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
        USE io_global,     ONLY : ionode_id
        USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool
        USE mp_world,      ONLY : mpime
        USE division,  ONLY : fkbounds
        USE wan2bloch, ONLY : dynwan2bloch, dynifc2blochf, hamwan2bloch, ephwan2blochp, ephwan2bloch
        USE rigid_epw,     ONLY : rpa_epsilon, tf_epsilon, compute_umn_f, rgd_blk_epw_fine
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT (IN) :: iq, nrr_k, nrr_q, nrr_g
        INTEGER, INTENT (IN) :: irvec_q(:,:), irvec_g(:,:)
        INTEGER, INTENT (IN) :: ndegen_k(:,:,:), ndegen_q(:,:,:), ndegen_g(:,:,:,:)
        REAL(KIND=dp), INTENT (INOUT) :: w2(3*nat)
        COMPLEX(KIND=dp), INTENT (INOUT) :: uf ( nmodes, nmodes), epmatwef( nbndsub, nbndsub, nrr_k, nmodes)
        REAL(KIND=dp), INTENT (IN) :: irvec_r(3,nrr_k)
        !
        ! Local  variables
        !! FIXME: dims should be nbnd_sub and intent(in)
        INTEGER, INTENT (IN) :: dims
        !! Dims is either nbndsub if use_ws or 1 if not
        INTEGER, INTENT (IN) :: dims2
        !! Dims is either nat if use_ws or 1 if not
        COMPLEX(KIND=dp) :: cfac(nrr_k, dims, dims)
        !! Used to store $e^{2\pi r \cdot k}$ exponential
        COMPLEX(KIND=dp) :: cfacq(nrr_k, dims, dims)
        !! Used to store $e^{2\pi r \cdot k+q}$ exponential
        COMPLEX(KIND=dp) :: cufkk ( nbndsub, nbndsub )
        !! Rotation matrix, fine mesh, points k
        COMPLEX(KIND=dp) :: cufkq (  nbndsub, nbndsub )
        !! the same, for points k+q
        COMPLEX(KIND=dp) :: epmatf( nbndsub, nbndsub, nmodes)
        !! e-p matrix  in smooth Bloch basis, fine mesh
        COMPLEX(KIND=dp) :: bmatf ( nbndsub, nbndsub)
        !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
        INTEGER :: nksqtotf
        !! half of the total number of k+q points (fine grid)
        INTEGER :: lower_bnd
        !! lower bound for the k-depend index among the mpi pools
        INTEGER :: upper_bnd
        !! lower bound for the k-depend index among the mpi pools
        INTEGER :: ik
        !! Counter on coarse k-point grid
        INTEGER :: ikk
        !! Counter on k-point when you have paired k and q
        INTEGER :: ikq
        !! Paired counter so that q is adjacent to its k
        INTEGER :: ibnd
        !! Counter on band
        INTEGER :: jbnd
        !! Counter on band
        INTEGER :: na
        !! Counter on atom
        INTEGER :: mu
        !! counter on mode
        INTEGER :: nu
        !! counter on mode
        INTEGER :: fermicount
        !! Number of states at the Fermi level
        INTEGER :: nrws
        !! Number of real-space Wigner-Seitz
        INTEGER, PARAMETER :: nrwsx=200
        !! Maximum number of real-space Wigner-Seitz
        REAL(KIND=dp) :: xxq(3)
        !! Current q-point
        REAL(KIND=dp) :: xxk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkq(3)
        !! Current k+q point on the fine grid
        REAL(KIND=dp) :: rws(0:3,nrwsx)
        !! Real-space wigner-Seitz vectors
        REAL(KIND=dp), PARAMETER :: eps = 0.01/ryd2mev
        !! Tolerence

    END SUBROUTINE gkq

    SUBROUTINE phonon_eigvector ( iq, nrr_k, nrr_q, irvec_q, ndegen_q, w2, uf, epmatwef)
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
            et, xk, ef,  nelec
        USE cell_base,     ONLY : at, bg, omega, alat
        USE start_k,       ONLY : nk1, nk2, nk3
        USE ions_base,     ONLY : nat, amass, ityp, tau
        USE phcom,         ONLY : nq1, nq2, nq3
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange, &
            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen,&
            etf_mem, scr_typ,&
            elecselfen, phonselfen, nest_fn, a2f, specfun_ph, &
            vme, eig_read, ephwrite, nkf1, nkf2, nkf3, &
            efermi_read, fermi_energy, specfun_el, band_plot, &
            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk, &
            plselfen, specfun_pl, wfcelec
        USE noncollin_module, ONLY : noncolin
        USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
        USE io_files,      ONLY : prefix, diropn
        USE io_global,     ONLY : stdout, ionode
        USE elph2,         ONLY : cu, cuq, lwin, lwinq,&
            chw, chw_ks, cvmew, cdmew, rdw, &
            epmatwp, epmatq, etf, etf_k, etf_ks, xqf, xkf, &
            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
            sigmai_all, sigmai_mode, gamma_all, epsi, zstar, &
            efnew, ifc, sigmar_all, zi_all, nkqtotf, eps_rpa, &
            g2_4, wf, nbndskip
#if defined(__NAG)
        USE f90_unix_io,   ONLY : FLUSH
#endif
        USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
        USE io_global,     ONLY : ionode_id
        USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool
        USE mp_world,      ONLY : mpime
        USE division,  ONLY : fkbounds
        USE wan2bloch, ONLY : dynwan2bloch, dynifc2blochf, hamwan2bloch, ephwan2blochp
        !
        IMPLICIT NONE
        !
        !INTEGER :: nrws
        !REAL(kind=DP) :: xxq(3), w2(3*nat), rws(0:3,200)
        !COMPLEX(kind=DP) :: uf ( nmodes, nmodes)
        !
        INTEGER, INTENT (IN) :: iq !, nrr_q!, ndegen_q(20*nq1*nq2*nq3)
        INTEGER, INTENT (IN) :: nrr_k, nrr_q, irvec_q(:,:), ndegen_q(:,:,:) ! ! Added for polaron calculations by Chao Lian.

        INTEGER :: nksqtotf, lower_bnd, upper_bnd
        REAL(KIND=dp), INTENT (INOUT) :: w2(3*nat)
        COMPLEX(KIND=dp), INTENT (INOUT) :: uf ( nmodes, nmodes), epmatwef( nbndsub,nbndsub, nrr_k, nmodes)
        COMPLEX(KIND=dp) :: cfac(nrr_k), cfacq(nrr_k), cufkk ( nbndsub, nbndsub),cufkq ( nbndsub, nbndsub), &
            epmatf( nbndsub, nbndsub, nmodes), bmatf( nbndsub, nbndsub)
        !
        ! Local  variables
        LOGICAL :: already_skipped
        !! Skipping band during the Wannierization
        LOGICAL :: exst
        !! If the file exist
        LOGICAL :: opnd
        !! Check whether the file is open.
        !
        CHARACTER (LEN=256) :: filint
        !! Name of the file to write/read
        CHARACTER (LEN=256) :: namef
        !! Name of the file
        CHARACTER (LEN=30)  :: myfmt
        !! Variable used for formatting output
        !
        INTEGER :: ios
        !! integer variable for I/O control
        INTEGER :: ik
        !! Counter on coarse k-point grid
        INTEGER :: ikk
        !! Counter on k-point when you have paired k and q
        INTEGER :: ikq
        !! Paired counter so that q is adjacent to its k
        INTEGER :: ibnd
        !! Counter on band
        INTEGER :: jbnd
        !! Counter on band
        INTEGER :: imode
        !! Counter on mode
        INTEGER :: na
        !! Counter on atom
        INTEGER :: mu
        !! counter on mode
        INTEGER :: nu
        !! counter on mode
        INTEGER :: fermicount
        !! Number of states at the Fermi level
        INTEGER :: nrec
        !! record index when reading file
        INTEGER :: lrepmatw
        !! record length while reading file
        INTEGER :: i,j
        !! Index when writing to file
        INTEGER :: ikx
        !! Counter on the coase k-grid
        INTEGER :: ikfx
        !! Counter on the fine k-grid.
        INTEGER :: xkk1, xkq1
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk2, xkq2
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk3, xkq3
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: ir
        !! Counter for WS loop
        INTEGER :: nrws
        !! Number of real-space Wigner-Seitz
        INTEGER :: valuerss(2)
        !! Return virtual and resisdent memory from system
        INTEGER, PARAMETER :: nrwsx=200
        !! Maximum number of real-space Wigner-Seitz
        !
        REAL(KIND=dp) :: rdotk_scal
        !! Real (instead of array) for $r\cdot k$
        REAL(KIND=dp) :: xxq(3)
        !! Current q-point
        REAL(KIND=dp) :: xxk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkq(3)
        !! Current k+q point on the fine grid
        REAL(KIND=dp) :: rws(0:3,nrwsx)
        !! Real-space wigner-Seitz vectors
        REAL(KIND=dp) :: atws(3,3)
        !! Maximum vector: at*nq
        REAL(KIND=dp), EXTERNAL :: efermig
        !! External function to calculate the fermi energy
        REAL(KIND=dp), EXTERNAL :: efermig_seq
        !! Same but in sequential
        REAL(KIND=dp), PARAMETER :: eps = 0.01/ryd2mev
        !! Tolerence
        !
        COMPLEX(KIND=dp) :: tableqx (4*nk1+1,2*nkf1+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqy (4*nk2+1,2*nkf2+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqz (4*nk3+1,2*nkf3+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe  (:,:,:,:,:)
        !! e-p matrix  in wannier basis - electrons
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
        !! e-p matrix  in wannier basis - electrons (written on disk)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
        !
        !
    END SUBROUTINE phonon_eigvector

    SUBROUTINE bubble_sort(array,sizes,output,repeat_list)
        USE kinds,         ONLY : dp
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: sizes
        REAL(KIND=dp), INTENT (IN) :: array(sizes)
        REAL(KIND=dp) :: temp, input_array(sizes)
        INTEGER :: bubble, j, lsup, degen_label
        LOGICAL, INTENT(OUT) :: output
        INTEGER, INTENT(OUT) :: repeat_list(sizes)

    END SUBROUTINE bubble_sort

    SUBROUTINE compute_a_re ( iq, nrr_k, ndegen_k, irvec_r, dims)
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
            et, xk, ef,  nelec
        USE cell_base,     ONLY : at, bg, omega, alat
        USE start_k,       ONLY : nk1, nk2, nk3
        USE ions_base,     ONLY : nat, amass, ityp, atm, ntyp => nsp, tau
        USE phcom,         ONLY : nq1, nq2, nq3
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange, &
            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
            etf_mem, scr_typ, &
            elecselfen, phonselfen, nest_fn, a2f, specfun_ph, &
            vme, eig_read, ephwrite, nkf1, nkf2, nkf3, &
            efermi_read, fermi_energy, specfun_el, band_plot, &
            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk, &
            plselfen, specfun_pl, wfcelec, num_cbands
        USE noncollin_module, ONLY : noncolin
        USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
        USE io_files,      ONLY : prefix, diropn
        USE io_global,     ONLY : stdout, ionode, meta_ionode_id
        USE elph2,         ONLY : cu, cuq, lwin, lwinq,&
            chw, chw_ks, cvmew, cdmew, rdw, &
            epmatwp, epmatq, etf, etf_k, etf_ks, xqf, xkf, &
            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
            sigmai_all, sigmai_mode, gamma_all, epsi, zstar, &
            efnew, ifc, sigmar_all, zi_all, nkqtotf, eps_rpa, &
            g2_4, wf, nbndskip
#if defined(__NAG)
        USE f90_unix_io,   ONLY : FLUSH
#endif
        USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
        USE io_global,     ONLY : ionode_id
        USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool, &
            world_comm
        USE mp_world,      ONLY : mpime
        USE division,  ONLY : fkbounds
        USE wan2bloch, ONLY : dynwan2bloch, dynifc2blochf, hamwan2bloch, ephwan2blochp
        !
        IMPLICIT NONE
        !
        !INTEGER :: nrws
        !REAL(kind=DP) :: xxq(3), w2(3*nat), rws(0:3,200)
        !COMPLEX(kind=DP) :: uf ( nmodes, nmodes)
        !
        INTEGER, INTENT (IN) :: iq, nrr_k
        INTEGER, INTENT (IN) :: ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
        INTEGER, INTENT (IN) :: dims
        INTEGER :: nksqtotf, lower_bnd, upper_bnd
        REAL(KIND=dp), INTENT (IN) :: irvec_r(3,nrr_k)
        !
        ! Local  variables
        COMPLEX(KIND=dp) :: cfac(nrr_k), cfacq(nrr_k), cufkk ( nbndsub, nbndsub),cufkq ( nbndsub, nbndsub), &
            epmatf( nbndsub, nbndsub, nmodes), bmatf( nbndsub, nbndsub)
        LOGICAL :: already_skipped
        !! Skipping band during the Wannierization
        LOGICAL :: exst
        !! If the file exist
        LOGICAL :: opnd
        !! Check whether the file is open.
        !
        CHARACTER (LEN=256) :: filint
        !! Name of the file to write/read
        CHARACTER (LEN=256) :: namef
        !! Name of the file
        CHARACTER (LEN=30)  :: myfmt
        !! Variable used for formatting output
        !
        INTEGER :: ios
        !! integer variable for I/O control
        INTEGER :: ik
        !! Counter on coarse k-point grid
        INTEGER :: ikk
        !! Counter on k-point when you have paired k and q
        INTEGER :: ikq
        !! Paired counter so that q is adjacent to its k
        INTEGER :: ibnd
        !! Counter on band
        INTEGER :: jbnd
        !! Counter on band
        INTEGER :: imode
        !! Counter on mode
        INTEGER :: na
        !! Counter on atom
        INTEGER :: mu
        !! counter on mode
        INTEGER :: nu
        !! counter on mode
        INTEGER :: fermicount
        !! Number of states at the Fermi level
        INTEGER :: nrec
        !! record index when reading file
        INTEGER :: lrepmatw
        !! record length while reading file
        INTEGER :: i,j
        !! Index when writing to file
        INTEGER :: ikx
        !! Counter on the coase k-grid
        INTEGER :: ikfx
        !! Counter on the fine k-grid.
        INTEGER :: xkk1, xkq1
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk2, xkq2
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk3, xkq3
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: ir
        !! Counter for WS loop
        INTEGER :: nrws
        !! Number of real-space Wigner-Seitz
        INTEGER :: valuerss(2)
        !! Return virtual and resisdent memory from system
        INTEGER, PARAMETER :: nrwsx=200
        !! Maximum number of real-space Wigner-Seitz
        !
        REAL(KIND=dp) :: rdotk_scal
        !! Real (instead of array) for $r\cdot k$
        REAL(KIND=dp) :: xxq(3)
        !! Current q-point
        REAL(KIND=dp) :: xxk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkq(3)
        !! Current k+q point on the fine grid
        REAL(KIND=dp) :: rws(0:3,nrwsx)
        !! Real-space wigner-Seitz vectors
        REAL(KIND=dp) :: atws(3,3)
        !! Maximum vector: at*nq
        REAL(KIND=dp), EXTERNAL :: efermig
        !! External function to calculate the fermi energy
        REAL(KIND=dp), EXTERNAL :: efermig_seq
        !! Same but in sequential
        REAL(KIND=dp), PARAMETER :: eps = 0.01/ryd2mev
        !! Tolerence
        !
        COMPLEX(KIND=dp) :: tableqx (4*nk1+1,2*nkf1+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqy (4*nk2+1,2*nkf2+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqz (4*nk3+1,2*nkf3+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe  (:,:,:,:,:)
        !! e-p matrix  in wannier basis - electrons
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
        !! e-p matrix  in wannier basis - electrons (written on disk)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
        INTEGER :: ounit, nx,ny,nz, np, hh, sort_indice_cutoff, n1_dim, k
        INTEGER :: ne, irr, ncb, ib, ii, jj, kk
        INTEGER :: n1_dim_x,n1_dim_y,n1_dim_z,hhx,hhy,hhz,npx,npy,npz
        INTEGER :: grid_nature_s(3), grid_infor(9)
        COMPLEX(KIND=dp) :: z1, z2, au
        COMPLEX(KIND=dp), ALLOCATABLE :: ac(:), carica(:,:,:), am(:,:)
        REAL(KIND=dp) :: wannier_func, iso_gaussian, xcart(3), rcart(3), rcoor(3)
        REAL(KIND=dp) :: e1(3), e2(3), e3(3), x0(3), m1, m2, m3, xkkf(3)
        REAL(KIND=dp) :: deltax, deltay, deltaz, rdk
        CHARACTER*12 :: aclist
        CHARACTER(*), PARAMETER :: fileplace = "./grid/"
        !
    END SUBROUTINE compute_a_re

    SUBROUTINE interpol_a_k ( iq, nrr_k, ndegen_k, irvec_r, dims)
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
            et, xk, ef,  nelec
        USE cell_base,     ONLY : at, bg, omega, alat
        USE start_k,       ONLY : nk1, nk2, nk3
        USE ions_base,     ONLY : nat, amass, ityp, atm, ntyp => nsp, tau
        USE phcom,         ONLY : nq1, nq2, nq3
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange, &
            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
            etf_mem, scr_typ, &
            elecselfen, phonselfen, nest_fn, a2f, specfun_ph, &
            vme, eig_read, ephwrite, nkf1, nkf2, nkf3, &
            efermi_read, fermi_energy, specfun_el, band_plot, &
            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk, &
            plselfen, specfun_pl, wfcelec, num_cbands
        USE noncollin_module, ONLY : noncolin
        USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
        USE io_files,      ONLY : prefix, diropn
        USE io_global,     ONLY : stdout, ionode, meta_ionode_id
        USE elph2,         ONLY : cu, cuq, lwin, lwinq,&
            chw, chw_ks, cvmew, cdmew, rdw, &
            epmatwp, epmatq, etf, etf_k, etf_ks, xqf, xkf, &
            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
            sigmai_all, sigmai_mode, gamma_all, epsi, zstar, &
            efnew, ifc, sigmar_all, zi_all, nkqtotf, eps_rpa, &
            g2_4, wf, nbndskip
#if defined(__NAG)
        USE f90_unix_io,   ONLY : FLUSH
#endif
        USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
        USE io_global,     ONLY : ionode_id
        USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool, &
            world_comm
        USE mp_world,      ONLY : mpime
        USE division,  ONLY : fkbounds
        USE wan2bloch, ONLY : dynwan2bloch, dynifc2blochf, hamwan2bloch, ephwan2blochp
        !
        IMPLICIT NONE
        !
        !INTEGER :: nrws
        !REAL(kind=DP) :: xxq(3), w2(3*nat), rws(0:3,200)
        !COMPLEX(kind=DP) :: uf ( nmodes, nmodes)
        !
        INTEGER, INTENT (IN) :: iq !,!, ndegen_q(20*nq1*nq2*nq3)
        INTEGER, INTENT (IN) :: dims
        INTEGER, INTENT (IN) :: nrr_k, ndegen_k(:,:,:)! ! Added for polaron calculations by Chao Lian.
        INTEGER :: nksqtotf, lower_bnd, upper_bnd
        REAL(KIND=dp), INTENT (IN) :: irvec_r(3,nrr_k)
        !
        ! Local  variables
        COMPLEX(KIND=dp) :: cfac(nrr_k), cfacq(nrr_k), cufkk ( nbndsub, nbndsub),cufkq( nbndsub, nbndsub), &
            epmatf( nbndsub, nbndsub, nmodes), bmatf( nbndsub,nbndsub)
        LOGICAL :: already_skipped
        !! Skipping band during the Wannierization
        LOGICAL :: exst
        !! If the file exist
        LOGICAL :: opnd
        !! Check whether the file is open.
        !
        CHARACTER (LEN=256) :: filint
        !! Name of the file to write/read
        CHARACTER (LEN=256) :: namef
        !! Name of the file
        CHARACTER (LEN=30)  :: myfmt
        !! Variable used for formatting output
        !
        INTEGER :: ios
        !! integer variable for I/O control
        INTEGER :: ik
        !! Counter on coarse k-point grid
        INTEGER :: ikk
        !! Counter on k-point when you have paired k and q
        INTEGER :: ikq
        !! Paired counter so that q is adjacent to its k
        INTEGER :: ibnd
        !! Counter on band
        INTEGER :: jbnd
        !! Counter on band
        INTEGER :: imode
        !! Counter on mode
        INTEGER :: na
        !! Counter on atom
        INTEGER :: mu
        !! counter on mode
        INTEGER :: nu
        !! counter on mode
        INTEGER :: fermicount
        !! Number of states at the Fermi level
        INTEGER :: nrec
        !! record index when reading file
        INTEGER :: lrepmatw
        !! record length while reading file
        INTEGER :: i,j
        !! Index when writing to file
        INTEGER :: ikx
        !! Counter on the coase k-grid
        INTEGER :: ikfx
        !! Counter on the fine k-grid.
        INTEGER :: xkk1, xkq1
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk2, xkq2
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk3, xkq3
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: ir
        !! Counter for WS loop
        INTEGER :: nrws
        !! Number of real-space Wigner-Seitz
        INTEGER :: valuerss(2)
        !! Return virtual and resisdent memory from system
        INTEGER, PARAMETER :: nrwsx=200
        !! Maximum number of real-space Wigner-Seitz
        !
        REAL(KIND=dp) :: rdotk_scal
        !! Real (instead of array) for $r\cdot k$
        REAL(KIND=dp) :: xxq(3)
        !! Current q-point
        REAL(KIND=dp) :: xxk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkq(3)
        !! Current k+q point on the fine grid
        REAL(KIND=dp) :: rws(0:3,nrwsx)
        !! Real-space wigner-Seitz vectors
        REAL(KIND=dp) :: atws(3,3)
        !! Maximum vector: at*nq
        REAL(KIND=dp), EXTERNAL :: efermig
        !! External function to calculate the fermi energy
        REAL(KIND=dp), EXTERNAL :: efermig_seq
        !! Same but in sequential
        REAL(KIND=dp), PARAMETER :: eps = 0.01/ryd2mev
        !! Tolerence
        !
        COMPLEX(KIND=dp) :: tableqx (4*nk1+1,2*nkf1+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqy (4*nk2+1,2*nkf2+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqz (4*nk3+1,2*nkf3+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe  (:,:,:,:,:)
        !! e-p matrix  in wannier basis - electrons
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
        !! e-p matrix  in wannier basis - electrons (written on disk)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
        INTEGER :: ounit, nx,ny,nz, hh, sort_indice_cutoff, n1_dim, k
        INTEGER :: ne, irr, ncb, ib, ii, jj, kk, nkf_global
        INTEGER :: n1_dim_x,n1_dim_y,n1_dim_z,hhx,hhy,hhz,npx,npy,npz
        INTEGER :: grid_nature_s(3), grid_infor(9)
        COMPLEX(KIND=dp) :: z1, z2, au
        COMPLEX(KIND=dp), ALLOCATABLE :: ac(:,:), carica(:,:,:), am(:), ac_full(:,:)
        REAL(KIND=dp) :: wannier_func, iso_gaussian, xcart(3), rcart(3), rcoor(3)
        REAL(KIND=dp) :: e1(3), e2(3), e3(3), x0(3), m1, m2, m3, xkkf(3)
        REAL(KIND=dp) :: deltax, deltay, deltaz, rdk
        CHARACTER*12 :: aclist
        CHARACTER(*), PARAMETER :: fileplace = "./grid/"
        !
    END SUBROUTINE interpol_a_k

    SUBROUTINE interpol_bq ( iq, nrr_k, nrr_q, nrr_g, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g, &
            w2, uf, epmatwef, irvec_r, dims, dims2)
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
            et, xk, ef,  nelec
        USE cell_base,     ONLY : at, bg, omega, alat
        USE start_k,       ONLY : nk1, nk2, nk3
        USE ions_base,     ONLY : nat, amass, ityp, atm, ntyp => nsp, tau
        USE phcom,         ONLY : nq1, nq2, nq3
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange, &
            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
            etf_mem, scr_typ, &
            elecselfen, phonselfen, nest_fn, a2f, specfun_ph, &
            vme, eig_read, ephwrite, nkf1, nkf2, nkf3, &
            efermi_read, fermi_energy, specfun_el, band_plot, &
            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk, &
            plselfen, specfun_pl, wfcelec, num_cbands
        USE noncollin_module, ONLY : noncolin
        USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
        USE io_files,      ONLY : prefix, diropn
        USE io_global,     ONLY : stdout, ionode, meta_ionode_id
        USE elph2,         ONLY : cu, cuq, lwin, lwinq,&
            chw, chw_ks, cvmew, cdmew, rdw, &
            epmatwp, epmatq, etf, etf_k, etf_ks, xqf, xkf, &
            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
            sigmai_all, sigmai_mode, gamma_all, epsi, zstar, &
            efnew, ifc, sigmar_all, zi_all, nkqtotf, eps_rpa, &
            g2_4, wf, nbndskip
#if defined(__NAG)
        USE f90_unix_io,   ONLY : FLUSH
#endif
        USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
        USE io_global,     ONLY : ionode_id
        USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool, &
            world_comm
        USE mp_world,      ONLY : mpime
        USE division,  ONLY : fkbounds
        USE wan2bloch, ONLY : dynwan2bloch, dynifc2blochf, hamwan2bloch, ephwan2blochp
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT (IN) :: iq, nrr_k, nrr_q, nrr_g
        INTEGER, INTENT (IN) :: irvec_q(:,:), irvec_g(:,:)
        INTEGER, INTENT (IN) :: ndegen_k(:,:,:), ndegen_q(:,:,:), ndegen_g(:,:,:,:)
        REAL(KIND=dp), INTENT (INOUT) :: w2(3*nat)
        COMPLEX(KIND=dp), INTENT (INOUT) :: uf ( nmodes, nmodes), epmatwef( nbndsub, nbndsub, nrr_k, nmodes)
        REAL(KIND=dp), INTENT (IN) :: irvec_r(3,nrr_k)
        !
        ! Local  variables
        !! FIXME: dims should be nbnd_sub and intent(in)
        INTEGER, INTENT (IN) :: dims
        !! Dims is either nbndsub if use_ws or 1 if not
        INTEGER, INTENT (IN) :: dims2
        !! Dims is either nat if use_ws or 1 if not

        INTEGER :: nksqtotf, lower_bnd, upper_bnd
        INTEGER :: n1_dim_x,n1_dim_y,n1_dim_z,hhx,hhy,hhz,npx,npy,npz
        INTEGER :: grid_nature_s(3), grid_infor(9)

        COMPLEX(KIND=dp) :: cfac(nrr_k), cfacq(nrr_k), cufkk ( nbndsub,nbndsub),cufkq( nbndsub, nbndsub), &
            epmatf( nbndsub, nbndsub, nmodes), bmatf( nbndsub,nbndsub)
        CHARACTER*12 :: aclist
        CHARACTER(*), PARAMETER :: fileplace = "./grid/"
        !
        ! Local  variables
        LOGICAL :: already_skipped
        !! Skipping band during the Wannierization
        LOGICAL :: exst
        !! If the file exist
        LOGICAL :: opnd
        !! Check whether the file is open.
        !
        CHARACTER (LEN=256) :: filint
        !! Name of the file to write/read
        CHARACTER (LEN=256) :: namef
        !! Name of the file
        CHARACTER (LEN=30)  :: myfmt
        !! Variable used for formatting output
        !
        INTEGER :: ios
        !! integer variable for I/O control
        INTEGER :: ik
        !! Counter on coarse k-point grid
        INTEGER :: ikk
        !! Counter on k-point when you have paired k and q
        INTEGER :: ikq
        !! Paired counter so that q is adjacent to its k
        INTEGER :: ibnd
        !! Counter on band
        INTEGER :: jbnd
        !! Counter on band
        INTEGER :: imode
        !! Counter on mode
        INTEGER :: na
        !! Counter on atom
        INTEGER :: mu
        !! counter on mode
        INTEGER :: nu
        !! counter on mode
        INTEGER :: fermicount
        !! Number of states at the Fermi level
        INTEGER :: nrec
        !! record index when reading file
        INTEGER :: lrepmatw
        !! record length while reading file
        INTEGER :: i,j
        !! Index when writing to file
        INTEGER :: ikx
        !! Counter on the coase k-grid
        INTEGER :: ikfx
        !! Counter on the fine k-grid.
        INTEGER :: xkk1, xkq1
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk2, xkq2
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: xkk3, xkq3
        !! Integer of xkk when multiplied by nkf/nk
        INTEGER :: ir
        !! Counter for WS loop
        INTEGER :: nrws
        !! Number of real-space Wigner-Seitz
        INTEGER :: valuerss(2)
        !! Return virtual and resisdent memory from system
        INTEGER, PARAMETER :: nrwsx=200
        !! Maximum number of real-space Wigner-Seitz
        !
        REAL(KIND=dp) :: rdotk_scal
        !! Real (instead of array) for $r\cdot k$
        REAL(KIND=dp) :: xxq(3)
        !! Current q-point
        REAL(KIND=dp) :: xxk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkk(3)
        !! Current k-point on the fine grid
        REAL(KIND=dp) :: xkq(3)
        !! Current k+q point on the fine grid
        REAL(KIND=dp) :: rws(0:3,nrwsx)
        !! Real-space wigner-Seitz vectors
        REAL(KIND=dp) :: atws(3,3)
        !! Maximum vector: at*nq
        REAL(KIND=dp), EXTERNAL :: efermig
        !! External function to calculate the fermi energy
        REAL(KIND=dp), EXTERNAL :: efermig_seq
        !! Same but in sequential
        REAL(KIND=dp), PARAMETER :: eps = 0.01/ryd2mev
        !! Tolerence
        !
        COMPLEX(KIND=dp) :: tableqx (4*nk1+1,2*nkf1+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqy (4*nk2+1,2*nkf2+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp) :: tableqz (4*nk3+1,2*nkf3+1)
        !! Look-up table for the exponential (speed optimization) in the case of
        !! homogeneous grids.
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe  (:,:,:,:,:)
        !! e-p matrix  in wannier basis - electrons
        COMPLEX(KIND=dp), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
        !! e-p matrix  in wannier basis - electrons (written on disk)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
        !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
        INTEGER :: ounit, nx,ny,nz,  hh, sort_indice_cutoff, n1_dim, k
        INTEGER :: ne, irr, ncb, ib, ii, jj, kk
        COMPLEX(KIND=dp) :: z1, z2, au
        COMPLEX(KIND=dp), ALLOCATABLE :: ac(:,:), carica(:,:,:), am(:), ac_full(:,:), &
            ac_read(:), bq(:)
        REAL(KIND=dp) :: wannier_func, iso_gaussian, xcart(3), rcart(3), rcoor(3)
        REAL(KIND=dp) :: e1(3), e2(3), e3(3), x0(3), m1, m2, m3, xkkf(3)
        REAL(KIND=dp) :: deltax, deltay, deltaz, rdk
        !
    END SUBROUTINE interpol_bq

    SUBROUTINE get_cfac(xk, nrr_k, ndegen_k, irvec_r, dims, cfac)
        USE epwcom, ONLY : use_ws
        USE constants_epw, ONLY : twopi, ci, czero
        USE kinds,         ONLY : dp, i4b

        IMPLICIT NONE


        INTEGER, INTENT(IN):: nrr_k, dims
        INTEGER, INTENT(IN):: ndegen_k(nrr_k, dims, dims)
        REAL(KIND=dp), INTENT (IN) :: xk(3), irvec_r(3, nrr_k)
        COMPLEX(KIND=dp), INTENT(OUT) :: cfac(nrr_k, dims, dims)
        ! Local Variables
        REAL(KIND=dp) :: rdotk(nrr_k)
        INTEGER:: ikk, ikq, iw, iw2, ir

        cfac = czero
        rdotk = czero

        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xk, 1, 0.0_dp, rdotk, 1 )
        !
        IF (use_ws) THEN
            DO iw=1, dims
                DO iw2=1, dims
                    DO ir = 1, nrr_k
                        IF (ndegen_k(ir,iw2,iw) > 0) THEN
                            cfac(ir,iw2,iw)  = EXP( ci*rdotk(ir) ) / ndegen_k(ir,iw2,iw)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            cfac(:,1,1)   = EXP( ci*rdotk(:) ) / ndegen_k(:,1,1)
        ENDIF
    END SUBROUTINE
SUBROUTINE ksstate_extract ( )
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : dp
    USE mp_global,     ONLY : my_pool_id, nproc_pool,    &
        intra_pool_comm, &
        inter_pool_comm, inter_image_comm, world_comm
    USE wavefunctions,  ONLY: evc
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE pwcom,         ONLY : nks, nkstot
    USE elph2,         ONLY : ngk_all, igk_k_all, xqf, nqf, xkf, nkf, nkqtotf
    USE cell_base,  ONLY : at, alat, celldm
    USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau
    USE gvect,         ONLY : g, ngm
    USE constants, ONLY :  pi
    USE constants_epw, ONLY : twopi
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE wvfct,     ONLY : npwx
    USE division,  ONLY : fkbounds
    USE io_epw,           ONLY : readwfc
    !
    INTEGER :: ipooltmp, ik, iq, npw, lower_bnd, upper_bnd, ig, igp, igpp, nksqtotf
    INTEGER :: ounit, nx,ny,nz, np, hh, sort_indice_cutoff, n1_dim, i, j, k
    REAL(KIND=dp) :: cg2, temp1(3), temp2(3)
    REAL(KIND=dp) :: e1(3), e2(3), e3(3), x0(3), m1, m2, m3
    REAL(KIND=dp) :: deltax, deltay, deltaz
    COMPLEX(KIND=dp), ALLOCATABLE :: carica (:,:,:)
    INTEGER, ALLOCATABLE :: igk(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: eigx (:), eigy (:), eigz (:)
    COMPLEX(KIND=dp), ALLOCATABLE :: ac(:), acp(:)

END SUBROUTINE ksstate_extract

!----------------------------------------------------------------------------
END MODULE

MODULE polaron_old
    PUBLIC  :: wfc_elec_old
CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE wfc_elec_old ( nrr_k, nrr_q, nrr_g, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g, &
            w2, uf, epmatwef, irvec_r, dims, dims2 )
        !-----------------------------------------------------------------------
        !
        !  Compute the polaron envelop function and formation energy (Pekar's continuum model,
        !  generalized Frohlich vertex )
        !  DS,CV
        !
        !  Use effective mass, static dielectric permittivity and Born effective
        !  charge as input variables
        !
        !  Use matrix elements, electronic eigenvalues and phonon frequencies
        !  from ep-wannier interpolation
        !
        !  This subroutine computes the contribution from phonon iq to all k-points
        !  The outer loop in ephwann_shuffle.f90 will loop over all iq points
        !  The contribution from each iq is summed at the end of this subroutine for iq=nqtotf
        !  to recover the per-ik electron wavefunction
        !
        !  Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
        !-----------------------------------------------------------------------
        USE kinds,         ONLY : dp
        use test_tools, only : para_write
        USE io_global,     ONLY : stdout,ionode_id, meta_ionode_id
        USE modes,         ONLY : nmodes
        USE epwcom,        ONLY : nbndsub, shortrange, restart_polaron,&
            fsthick, ngaussw, degaussw,spherical_cutoff,&
            eps_acustic, efermi_read, fermi_energy, lscreen, &
            model_vertex, nkf1, nkf2, nkf3, conv_thr_polaron, &
            r01, r02, r03, num_cbands, start_mode, cb_shift, &
            polaron_dos, polaron_type, &
            electron_dos, phonon_dos, diag_mode, restart_polaron_mode
        USE pwcom,         ONLY : ef !,nelec, isk
        USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
            nkf, nqf,epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf, &
            efnew, eps_rpa, g2_all, &
            ac, hkk, ec, ekf, gq, n1_dim,hh,np, &
            g2_4
        USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, twopi,eps6,&
            czero, cone
        USE mp,            ONLY : mp_barrier, mp_sum,mp_bcast
        USE mp_global,     ONLY : inter_pool_comm
        USE mp_world,      ONLY : mpime, world_comm
        USE ions_base,     ONLY : nat, tau
        USE start_k,       ONLY : nk1, nk2, nk3
        USE cell_base,     ONLY : at, bg, alat, omega
        USE parallel_include
        USE division,  ONLY : fkbounds
        USE ephblochkq, ONLY: gkq, phonon_eigvector
        USE poolgathering, ONLY : poolgather2
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT (IN) :: nrr_k, nrr_q, nrr_g ! Added for polaron calculations by Chao Lian.
        REAL(KIND=dp), INTENT (INOUT) :: w2(3*nat)
        INTEGER, INTENT (IN) :: irvec_q(:,:), irvec_g(:,:)
        REAL(KIND=dp), ALLOCATABLE :: irvec_r(:,:)
        INTEGER, INTENT (IN) :: ndegen_k(:,:,:), ndegen_q(:,:,:), ndegen_g(:,:,:,:)
        COMPLEX(KIND=dp), INTENT (INOUT) :: uf ( nmodes, nmodes), epmatwef( nbndsub, nbndsub, nrr_k, nmodes)
        INTEGER, INTENT (IN) :: dims
        !! Dims is either nbndsub if use_ws or 1 if not
        INTEGER, INTENT (IN) :: dims2
        !! Dims is either nat if use_ws or 1 if not
        !
        INTEGER :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ir
        COMPLEX(KIND=dp) :: cfac, weight, zdotu
        REAL(KIND=dp) :: g2, ekk, ekq, wq, ef0, wgkq, inv_eptemp0, w0g1, w0g2, &
            g2_tmp, inv_wq, inv_degaussw
        REAL(KIND=dp), EXTERNAL :: wgauss, w0gauss
        REAL(KIND=dp), PARAMETER :: eps2 = 2.d0/ryd2mev
        !
        ! variables for collecting data from all pools in parallel case
        !
        INTEGER :: nksqtotf, lower_bnd, upper_bnd
        REAL(KIND=dp), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:)
        !
        ! variables defined by DS
        INTEGER :: delta_function,counter,i,j,k,ii,n3,n4,sort_indice_even ! define the size of supercell
        INTEGER :: ncb, ikbnd, kbnd, delta, ikkk, sort_indice, band_pos
        INTEGER :: grid_infor(9)
        REAL(KIND=dp) :: diff_k(3) ,r0(3),a,c, omega_c, m_2, diff_kk(3),var, ecb0, tauu0(3)
        REAL(KIND=dp), ALLOCATABLE :: k_grid(:,:)  ! k grids used for matrix
        COMPLEX(KIND=dp), ALLOCATABLE :: z1(:), z2(:)
        COMPLEX(KIND=dp), ALLOCATABLE :: mt(:,:), vt(:), tm(:),outt(:),lac(:,:)

        !
        INTEGER :: sort_indice_cutoff, m1,m2,m3, na, ipol, grid_nature, igamma
        INTEGER :: lg(nqf)
        INTEGER, ALLOCATABLE :: dk_list(:), dkk_list(:), mygrow_list(:)
        INTEGER :: n1_dim_x,n1_dim_y,n1_dim_z,hhx,hhy,hhz,npx,npy,npz
        INTEGER :: grid_nature_s(3)
        !
        COMPLEX(KIND=dp) :: phase_f,zz, dtauu(3), apa
        COMPLEX(KIND=dp) :: tv(3),m(3,3), m_copy(3,3)!, g_frohlich0, g_frohlich, g_frohlich2 ! test variables, checking
        COMPLEX(KIND=dp) :: am(2,2), eig(2), am_copy(2), gv(nmodes,nqf)
        !
        INTEGER :: ierr !
        INTEGER :: n, nb  ! problem size and block size
        INTEGER :: myarows, myacols, mygrows, mygcols ! size of local subset of global matrix
        INTEGER :: myxrows, myxcols   ! size of local subset of global vector
        INTEGER :: myi, myj, rows, loc_i, loc_j, max_row, max_col
        INTEGER, ALLOCATABLE :: call_list(:), request_list(:,:)
        COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: mya,myx,myy, mymsg
        COMPLEX(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: myg


        INTEGER :: tar_p, ib, ibp, i_n, i_m, i_ki, i_kj
        INTEGER, EXTERNAL :: numroc   ! blacs routine
        INTEGER :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
        INTEGER :: p_col_i, p_row_i, pme, blacs_pnum, icaller, tar_prow, tar_pcol
        INTEGER :: p_row_i2, icaller2, tar_prow2, tar_pcol2, loc_i2, loc_j2
        INTEGER :: info    ! scalapack return value
        INTEGER, DIMENSION(9)   :: ides_a, ides_x, ides_y, ides_g ! matrix descriptors
        INTEGER, DIMENSION(2) :: dimsl
        CHARACTER*12 :: aclist, folder
        CHARACTER(*), PARAMETER :: fileplace = "./grid/"
        INTEGER :: ntot, mk
        REAL(KIND=dp) :: ds, he, sigma, dos
        REAL(KIND=dp) :: rp(3), qrp, qcart(3), e_formation2, e_formation1, wlo
        !!!
#if defined (__MPI)
        COMPLEX(KIND=dp), ALLOCATABLE :: array(:,:,:), VAL(:,:,:), nval(:,:,:)
        INTEGER (KIND = mpi_address_kind) :: SIZE,lowerbound, sizeofreal, disp_aint,&
            disp_int
        INTEGER :: win, tar
#endif
    END SUBROUTINE wfc_elec_old
END MODULE
