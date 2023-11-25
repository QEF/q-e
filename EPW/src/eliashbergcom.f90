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
  MODULE eliashberg_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for real- and imag-axis Eliashberg equations
  !! SP: Sept. 2019 - Cleaning.
  !!
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER :: nsw
  !! Nr. of grid points between (0,wscut) for real-axis, analytical continuation and Pade approximants
  INTEGER :: ndos
  !! Nr. of energy bins in Fermi window for dos
  INTEGER, ALLOCATABLE :: nsiw(:)
  !! Nr. of grid points at each temperature on imag-axis, nsiw(nstemp)
  INTEGER, ALLOCATABLE :: wsn(:)
  !! frequency "indices" on imag-axis at iw, wsn(nsiw(nstemp))
  !
  REAL(KIND = DP) :: wsphmax
  !! maximum phonon frequency for evaluation of the integral over Omega (0, wsphmax)
  REAL(KIND = DP) :: dwsph
  !! frequency step for Eliashberg spectral function
  REAL(KIND = DP) :: gap0
  !! initial guess for delta
  REAL(KIND = DP) :: muintr
  !! superconducting (interacting) chemical potential
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
  !
  !--------------------------------------------------------------------------
  END MODULE eliashberg_common
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  MODULE eliashberg_common_iso
  !--------------------------------------------------------------------------
  !!
  !! Global variables for real and imag-axis isotropic equations Eliashberg equations
  !!
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(KIND = DP), ALLOCATABLE :: a2f_iso(:)
  !! isotropic Eliashberg spectral function a2f_iso(nqstep)
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
  !
  !--------------------------------------------------------------------------
  END MODULE eliashberg_common_iso
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  MODULE eliashberg_common_aniso
  !--------------------------------------------------------------------------
  USE kinds, ONLY : DP
  !!
  !! Global variables for real and imag-axis anisotropic equations Eliashberg equations
  !! nkf = nr of irreducible k-points on the fine grid, if mp_mesh_k = .TRUE.
  !! nkf = total nr of k-points on the fine grid, otherwise
  !!
  SAVE
  !
  LOGICAL :: limag_fly
  !! .TRUE. if anisotropic Eliashberg eqns on imag-axis are solved on the fly
  LOGICAL :: lacon_fly
  !! .TRUE. if anisotropic Eliashberg eqns on real-axis are solved on the fly
  !
  INTEGER :: nkfs
  !! nr. of irreducible k-points within the Fermi shell on the fine mesh
  INTEGER :: nkfs_all
  !! nr. of irreducible k-points within the Fermi shell on the fine mesh
  INTEGER :: nbndfs
  !! nr. of electronic bands within the Fermi shell
  INTEGER :: nbndfs_all
  !! nr. of bands of ekfs_all; it should be same with nbndfst
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
  !
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
  !
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
  !
  !--------------------------------------------------------------------------
  END MODULE eliashberg_common_aniso
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  MODULE eliashbergcom
  !--------------------------------------------------------------------------
  !
  USE eliashberg_common
  USE eliashberg_common_iso
  USE eliashberg_common_aniso
  !
  !--------------------------------------------------------------------------
  END MODULE eliashbergcom
  !--------------------------------------------------------------------------
