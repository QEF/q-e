  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !--------------------------------------------------------------------------
  MODULE supercond_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for real- and imag-axis Eliashberg equations
  !! SP: Sept. 2019 - Cleaning.
  !!
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: limag_fly
  !! .TRUE. if anisotropic Eliashberg eqns on imag-axis are solved on the fly
  LOGICAL :: lacon_fly
  !! .TRUE. if anisotropic Eliashberg eqns on real-axis are solved on the fly
  !
  INTEGER :: nsw
  !! Nr. of grid points between (0,wscut) for real-axis, analytical continuation and Pade approximants
  INTEGER :: ndos
  !! Nr. of energy bins in Fermi window for dos
  INTEGER :: nkfs
  !! nr. of irreducible k-points within the Fermi shell on the fine mesh
  INTEGER :: nkfs_all
  !! nr. of irreducible k-points within the Fermi shell on the fine mesh
  INTEGER :: nbndfs
  !! nr. of electronic bands within the Fermi shell
  INTEGER :: nbndfs_all
  !! nr. of bands of ekfs_all; it should be same with nbndfst
  INTEGER :: size_ir = 20
  !! size_ir is the size of arrays related to sparse-ir sampling. 
  !! Currently, the value of size_ir is assigned to siz_ir and siz_ir_cl in eliashberg_init.
  !! In our code, matrix-vector operations like y=Ax are frequently performed for each vector x, 
  !! where A remains constant. To optimize these linear transformations, the LAPACK subroutines 
  !! (ZGEMM or DGEMM) are utilized. 
  !! Therefore, it's more efficient to compute multiple vectors x together, as illustrated below:
  !!   [y1, y2, y3, ... , yn] = A [x1, x2, x3, ... , xn],
  !! where [x1, x2, x3, ... , xn] forms a matrix with x as its column vectors. 
  !! size_ir determines the number of vectors to be processed simultaneously. 
  !! Although vectors xn and yn are dependent on the final states (n', k+q), 
  !! memory limitations prevent allocating arrays as large as the number of final states.
  !! Test calculations have confirmed that the computation speed does not change that much even if the size is larger than 20.
  !! Depending on LAPACK and compiler, increasing size_ir beyond 20 might result in faster computation.
  INTEGER :: nlambda
  !! lambda = 10 ** (nlambda) determines wcut for each temperature: wcut = lambda * temperature
  INTEGER :: ndigit
  !! eps = 10 ** (-ndigit)
  INTEGER, ALLOCATABLE :: ixkff(:)
  !! index of k-point on the full k-grid ixkff(nkftot)
  INTEGER, ALLOCATABLE :: ixkf(:)
  !! index of k-point on the irreducible k-grid within the Fermi shell ixkf(nkf)
  INTEGER, ALLOCATABLE :: ixkf_inv(:)
  !! index of k-point on the irreducible k-grid ixkf_inv(nkfs) (for ekfs_all)
  INTEGER, ALLOCATABLE :: ixkqf(:, :)
  !! index k+q or k-q on the irreducilble k-grid within the Fermi shell ixkqf(nkfs,nqftot)
  INTEGER, ALLOCATABLE :: ixqfs(:, :)
  !! index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell ixqfs(nkfs,nqfs(ik))
  INTEGER, ALLOCATABLE :: nqfs(:)
  !! nr of q-points at each k-point for which k+sign*q is within the Fermi shell nqfs(nkfs)
  INTEGER, ALLOCATABLE :: ibnd_kfs_to_kfs_all(:, :)
  !! the function to convert an index of band of ekfs to the one of ekfs_all.
  INTEGER, ALLOCATABLE :: ibnd_kfs_all_to_kfs(:, :)
  !! the function to convert an index of band of ekf to the one of ekfs.
  !! ibnd_kfs_all_to_kfs(i, ik) == 0 means the corresponding the eigenvalue, ekfs_all(i, ik), is on the outside of fsthick window.
  INTEGER, ALLOCATABLE :: nsiw(:)
  !! Nr. of grid points at each temperature on imag-axis, nsiw(nstemp)
  INTEGER, ALLOCATABLE :: wsn(:)
  !! frequency "indices" on imag-axis at iw, wsn(nsiw(nstemp))
  INTEGER :: siz_ir
  !! siz_ir is the number of row of arrays related to sparse-ir (ir_giw, ...)
  INTEGER :: siz_ir_cl
  !! siz_ir_cl is the number of row of arrays related to sparse-ir (ir_giw_cl, ...)
  INTEGER, ALLOCATABLE :: num_js1(:, :)
  !! the number of states (n',k') inside the fsthick window for each state (n, k)
  INTEGER, ALLOCATABLE :: num_js2(:, :)
  !! the number of states (n',k') inside of [emin_coulomb: emax_coulomb] for each state (n, k)
  INTEGER, ALLOCATABLE :: num_js3(:)
  !! the number of states (n',k') inside of [emin_coulomb: emax_coulomb] for each state (n, k)
  INTEGER :: nk1_cl, nk2_cl, nk3_cl
  !! kx,ky,kz sizes of the uniform electron mesh of the data read from the nscf result
  INTEGER :: nbnd_cl
  !! nr. of bands of the data read from the nscf result
  INTEGER :: nkstot_cl
  !! total number of k points of the data read from the nscf result
  INTEGER :: nbnd_offset
  !! The band-index offset: we can obtain the index of ek_cl corresponding to ekf by adding the offset to the index of ekf.
  INTEGER, ALLOCATABLE :: ik_cl_to_fs(:)
  !! The function to convert an index of k points given by prefix.bands.out to an index of k points of the finer k grid.
  !! ik_cl_to_fs(ik_cl) == 0 means the corresponding k point does not have any bands within Fermi shell.
  INTEGER, ALLOCATABLE :: ik_fs_to_cl(:)
  !! the function to convert an index of k points of the finer k grid to an index of k points given by prefix.bands.out.
  !! ik_fs_to_cl(ik) == 0 means the corresponding k point is not on the coarser k grid given by the prefix.bands.out file.
  INTEGER, ALLOCATABLE :: ik_bz_to_ibz_cl(:)
  !! the function to convert an index of k points in the 1st BZ to an index of k points given by prefix.bands.out.
  INTEGER, ALLOCATABLE :: ik_ibz_to_bz_cl(:)
  !! the function to convert an index of k points given by prefix.bands.out to an index of k points in the 1st BZ.
  !
  REAL(KIND = DP) :: wsphmax
  !! maximum phonon frequency for evaluation of the integral over Omega (0, wsphmax)
  REAL(KIND = DP) :: dwsph
  !! frequency step for Eliashberg spectral function
  REAL(KIND = DP) :: gap0
  !! initial guess for delta
  REAL(KIND = DP) :: muintr
  !! superconducting (interacting) chemical potential
  REAL(KIND = DP) :: spin_fac
  !! If noncolin is true, spin_fac = 0.5; otherwise, spin_fac = 1.0.
  REAL(KIND = DP), ALLOCATABLE :: ws(:)
  !! frequency on real-axis, ws(nsw)
  REAL(KIND = DP), ALLOCATABLE :: wsph(:)
  !! frequency on real-axis, wsph(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: wsi(:)
  !! frequency on imag-axis at iw, wi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: en(:)
  !! Energy grid over Fermi window
  REAL(KIND = DP), ALLOCATABLE :: dosen(:)
  !! DOS (state/spin/eV/u.c.) over Fermi window
  REAL(KIND = DP), ALLOCATABLE :: a2f_tmp(:)
  !! Temporary isotropic Eliashberg spectral function a2f_tmp(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: fdwp(:)
  !! Fermi-Dirac distribution at frequency wp, fdwp(nsw)
  REAL(KIND = DP), ALLOCATABLE :: bewph(:)
  !! Bose-Einstein distribution at frequency wph, bewph(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: deltai(:)
  !! gap function on imag-axis at iw, deltai(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: deltaip(:)
  !! gap function on imag-axis at iwp, deltaip(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: znormi(:)
  !! renormalization function on the imag-axis at iw, znormi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: nznormi(:)
  !! normal state renormalization function on the imag-axis at iw, nznormi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: keri(:)
  !! phonon kernel on imag-axis, keri(2*nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: dsumi(:)
  !! contribution to delta eqn from the imaginary-axis in the analytic continuation dsumi(nsw)
  REAL(KIND = DP), ALLOCATABLE :: zsumi(:)
  !! contribution to znorm eqn from the imaginary-axis in the analytic continuation zsumi(nsw)
  REAL(KIND = DP), ALLOCATABLE :: gp(:, :)
  !! -bose(omegap)-fermi( omega+omegap) (eqn for delta and znorm analytic continuation)
  REAL(KIND = DP), ALLOCATABLE :: gm(:, :)
  !! bose(omegap)+fermi(-omega+omegap) (eqn for delta and znorm analytic continuation)
  REAL(KIND = DP), ALLOCATABLE :: znormip(:)
  !! renormalization function on imag-axis at iwp, znormip(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: shifti(:)
  !! energy shift on imag-axis at iw, shifti(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: shiftip(:)
  !! energy shift on imag-axis at iwp, shiftip(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: orderi(:)
  !! order paramter on the imag-axis at iw, orderi(nsiw(nstemp))
  REAL(KIND = DP) :: ef0
  !! Fermi energy
  REAL(KIND = DP) :: dosef
  !! density of states at the Fermi energy
  REAL(KIND = DP), ALLOCATABLE :: g2(:, :, :, :, :)
  !! e-ph matrix element squared |g_ji^nu(k,q)|^2, g2(nkfs_pool,nqftot,nbndfs,nbndfs,nmodes)
  REAL(KIND = DP), ALLOCATABLE :: ekfs(:, :)
  !! eigenvalues at E_i(k), etf(nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: xkff(:, :)
  !! coordintates of the k-points of the full k-grid xkff(3,nkftot)
  REAL(KIND = DP), ALLOCATABLE :: xkfs(:, :)
  !! coordintates of the k-points xkf(3,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: wkfs(:)
  !! weights of the irreducible k-points wkf(nkfs)
  REAL(KIND = DP), ALLOCATABLE :: ekfs_all(:, :)
  !! eigenvalues at E_i(k), etf(nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: xkfs_all(:, :)
  !! coordintates of the k-points xkf(3,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: wkfs_all(:)
  !! weights of the irreducible k-points wkf(nkfs)
  REAL(KIND = DP), ALLOCATABLE :: a2fij(:, :, :, :, :)
  !! spectral function a2fij(nqstep,nbndfs,nqftot,nbndfs,nkfs_pool)
  REAL(KIND = DP), ALLOCATABLE :: w0g(:, :)
  !! approximation for delta function w0g(nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: agap(:, :)
  !! superconducting gap edge agap(nkfs,nbndfs)
  REAL(KIND = DP), ALLOCATABLE :: adeltai(:, :, :)
  !! gap function on imag-axis at iw, adeltai(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: adeltaip(:, :, :)
  !! gap function on imag-axis at iwp, adeltaip(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: aznormi(:, :, :)
  !! renormalization function on imag-axis at iw, aznormi(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: naznormi(:, :, :)
  !! normal state renormalization function on imag-axis at iw, naznormi(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: akeri(:, :, :, :, :)
  !! phonon kernel on imag-axis, akeri(2*nsiw(nstemp),nbndfs,nqftot,nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: adsumi(:, :, :)
  !! contribution to delta eqn from the imaginary-axis in the analytic continuation adsumi(nsw,nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: azsumi(:, :, :)
  !! contribution to znorm eqn from the imaginary-axis in the analytic continuation azsumi(nsw,nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: memlt_pool(:)
  !! maximum allocatable memory per pool
  REAL(KIND = DP), ALLOCATABLE :: aznormip(:, :, :)
  !! renormalization function on imag-axis at iwp, aznormip(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: ashifti(:, :, :)
  !! energy shift on imag-axis at iw, ashifti(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: ashiftip(:, :, :)
  !! energy shift on imag-axis at iwp, ashiftip(nsiw(itemp),nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: gl_abs(:, :, :)
  !! absolute values of normal Green's functions within the fsthick window
  REAL(KIND = DP), ALLOCATABLE :: fl_abs(:, :, :)
  !! absolute values of anormalous Green's functions within the fsthick window
  REAL(KIND = DP), ALLOCATABLE :: knll_abs(:, :, :)
  !! absolute values of kernels within the fsthick window
  REAL(KIND = DP), ALLOCATABLE :: weight_q(:)
  !! To store wqf for each jx
  REAL(KIND = DP), ALLOCATABLE :: ir_gl_d(:,:)
  !! Green's function in the intermediate representations
  REAL(KIND = DP), ALLOCATABLE :: ir_gtau_d(:,:)
  !! Green's function of sampling imaginary time.
  REAL(KIND = DP), ALLOCATABLE :: ir_knll_d(:,:)
  !! kernel in the intermediate representations
  REAL(KIND = DP), ALLOCATABLE :: ir_knltau_d(:,:)
  !! kernel of sampling imaginary time.
  REAL(KIND = DP), ALLOCATABLE :: ir_cvll_d(:,:)
  !! convolution in the intermediate representations
  REAL(KIND = DP), ALLOCATABLE :: ir_cvltau_d(:,:)
  !! convolution of sampling imaginary time.
  REAL(KIND = DP), ALLOCATABLE :: ir_gl_cl_d(:,:)
  !! Green's function in the intermediate representations
  REAL(KIND = DP), ALLOCATABLE :: ir_gtau_cl_d(:,:)
  !! Green's function of sampling imaginary time.
  REAL(KIND = DP), ALLOCATABLE :: ek_cl(:, :)
  !! eigenvalues at E_i(k), ek_cl(nbnd_cl, nkstot_cl)
  REAL(KIND = DP), ALLOCATABLE :: xk_cl(:, :)
  !! coordintates of the k-points xk_cl(3, nkstot_cl)
  REAL(KIND = DP), ALLOCATABLE :: xk_bz_cl(:, :)
  !! coordintates of the k-points xk_bz_cl(3, nk1_cl * nk2_cl * nk3_cl)
  REAL(KIND = DP), ALLOCATABLE :: wk_cl(:)
  !! weights of the irreducible k-points wk_cl(nkstot_cl)
  REAL(KIND = DP), ALLOCATABLE :: adeltai_cl(:, :)
  !! gap function on imag-axis at iw, adeltai(nsiw(itemp),is_start:is_stop)
  REAL(KIND = DP), ALLOCATABLE :: adeltaip_cl(:, :, :)
  !! gap function on imag-axis at iwp, adeltaip(nsiw(itemp),nbnd_cl,nkstot_cl)
  REAL(KIND = DP), ALLOCATABLE :: w_stat(:, :)
  !! static part of coulomb interaction
  REAL(KIND = DP), ALLOCATABLE :: weight_cl(:)
  !! To store w_stat * (weight of q-point) for each jx
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: delta(:)
  !! gap function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: deltap(:)
  !! gap function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: znorm(:)
  !! renormalization function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: znormp(:)
  !! renormalization function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: kp(:, :)
  !! phonon kernel on real-axis (eqn for delta)
  COMPLEX(KIND = DP), ALLOCATABLE :: km(:, :)
  !! phonon kernel on real-axis (eqn for znorm)
  COMPLEX(KIND = DP), ALLOCATABLE :: shift(:)
  !! energy shift on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: aznorm(:, :, :)
  !! renormalization function on real-axis aznorm(nsw,nbndfs,nkfs)
  COMPLEX(KIND = DP), ALLOCATABLE :: aznormp(:, :, :)
  !! renormalization function on real-axis aznormkq(nsw,nbndfs,nkfs)
  COMPLEX(KIND = DP), ALLOCATABLE :: adelta(:, :, :)
  !! gap function on real-axis adelta(nsw,nbndfs,nkfs)
  COMPLEX(KIND = DP), ALLOCATABLE :: adeltap(:, :, :)
  !! gap function on real-axis adeltap(nsw,nbndfs,nkfs)
  COMPLEX(KIND = DP), ALLOCATABLE :: ashift(:, :, :)
  !! energy shift on real-axis ashift(nsw,nbndfs,nkfs)
  COMPLEX(KIND = DP), ALLOCATABLE :: fft_in1(:)
  !! inout array of FFT
  COMPLEX(KIND = DP), ALLOCATABLE :: fft_out1(:)
  !! outout array of FFT
  COMPLEX(KIND = DP), ALLOCATABLE :: fft_in2(:)
  !! inout array of FFT
  COMPLEX(KIND = DP), ALLOCATABLE :: fft_out2(:)
  !! outout array of FFT
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_giw(:,:)
  !! Green's function of sampling matsubara freq.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_gl(:,:)
  !! Green's function in the intermediate representations
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_gtau(:,:)
  !! Green's function of sampling imaginary time.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_knliw(:,:)
  !! kernel of sampling matsubara freq.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_knll(:,:)
  !! kernel in the intermediate representations
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_knltau(:,:)
  !! kernel of sampling imaginary time.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_cvliw(:,:)
  !! convolution of sampling matsubara freq.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_cvll(:,:)
  !! convolution in the intermediate representations
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_cvltau(:,:)
  !! convolution of sampling imaginary time.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_giw_cl(:,:)
  !! Green's function of sampling matsubara freq.
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_gl_cl(:,:)
  !! Green's function in the intermediate representations
  COMPLEX(KIND = DP), ALLOCATABLE :: ir_gtau_cl(:,:)
  !! Green's function of sampling imaginary time.
  !
  !--------------------------------------------------------------------------
  END MODULE supercond_common
  !--------------------------------------------------------------------------
