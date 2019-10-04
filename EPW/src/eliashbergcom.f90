  !                                                                            
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
  INTEGER, ALLOCATABLE :: nsiw(:)
  !! Nr of grid points at each temperature on imag-axis, nsiw(nstemp)
  REAL(KIND = DP) :: wsphmax
  !! maximum phonon frequency for evaluation of the integral over Omega (0, wsphmax)
  REAL(KIND = DP) :: dwsph
  !! frequency step for Eliashberg spectral function
  REAL(KIND = DP) :: gap0
  !! initial guess for Delta
  REAL(KIND = DP), ALLOCATABLE :: dws(:)
  !! grid size at each bin dws(nsw)
  REAL(KIND = DP), ALLOCATABLE :: ws(:)
  !! frequency on real-axis, ws(nsw)
  REAL(KIND = DP), ALLOCATABLE :: wsph(:)
  !! frequency on real-axis, wsph(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: wsi(:)
  !! frequency on imag-axis at iw, wi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: estemp(:)
  !! temperature in eV entering in the Eliashberg equtions estemp(nstemp)
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
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  REAL(KIND = DP), ALLOCATABLE :: a2f_iso(:)
  !! isotropic Eliashberg spectral function a2f_iso(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: gap(:)
  !! superconducting gap edge gap(nstemp)
  REAL(KIND = DP), ALLOCATABLE :: fdwp(:)
  !! Fermi-Dirac distribution at frequency wp, fdwp(nsw)
  REAL(KIND = DP), ALLOCATABLE :: bewph(:)
  !! Bose-Einstein distribution at frequency wph, bewph(nqstep)
  REAL(KIND = DP), ALLOCATABLE :: Deltai(:)
  !! gap function on imag-axis at iw, Deltai(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: Deltaip(:)
  !! gap function on imag-axis at iwp, Deltaip(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: Znormi(:)
  !! renormalization function on the imag-axis at iw, Znormi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: NZnormi(:)
  !! normal state renormalization function on the imag-axis at iw, NZnormi(nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: Keri(:)
  !! phonon kernel on imag-axis, Keri(2*nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: Dsumi(:)
  !! contribution to Delta eqn from the imaginary-axis in the analytic continuation Dsumi(nsw)
  REAL(KIND = DP), ALLOCATABLE :: Zsumi(:)
  !! contribution to Znorm eqn from the imaginary-axis in the analytic continuation Zsumi(nsw)
  REAL(KIND = DP), ALLOCATABLE :: Gp(:, :)
  !! -bose(omegap)-fermi( omega+omegap) (eqn for Delta and Znorm analytic continuation) 
  REAL(KIND = DP), ALLOCATABLE :: Gm(:, :) 
  !! bose(omegap)+fermi(-omega+omegap) (eqn for Delta and Znorm analytic continuation) 
  COMPLEX(KIND = DP), ALLOCATABLE :: Delta(:)
  !! gap function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: Deltap(:) 
  !! gap function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: Znorm(:)
  !! renormalization function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: Znormp(:)
  !! renormalization function on real-axis at iw
  COMPLEX(KIND = DP), ALLOCATABLE :: Kp(:, :)
  !! phonon kernel on real-axis (eqn for Delta)
  COMPLEX(KIND = DP), ALLOCATABLE :: Km(:, :)
  !! phonon kernel on real-axis (eqn for Znorm)
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
  !! FIXME
  LOGICAL :: lacon_fly
  !! FIXME
  INTEGER :: nkfs
  !! nr. of irreducible k-points within the Fermi shell on the fine mesh
  INTEGER :: nbndfs
  !! nr. of electronic bands within the Fermi shell
  INTEGER, ALLOCATABLE :: ixkff(:)
  !! index of k-point on the full k-grid ixkff(nkftot)
  INTEGER, ALLOCATABLE :: ixkf(:)
  !! index of k-point on the irreducible k-grid within the Fermi shell ixkf(nkf)
  INTEGER, ALLOCATABLE :: ixkqf(:, :)
  !! index k+q or k-q on the irreducilble k-grid within the Fermi shell ixkqf(nkfs,nqftot)
  INTEGER, ALLOCATABLE :: ixqfs(:, :)
  !! index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell ixqfs(nkfs,nqfs(ik))
  INTEGER, ALLOCATABLE :: nqfs(:)
  !! nr of q-points at each k-point for which k+sign*q is within the Fermi shell nqfs(nkfs)
  REAL(KIND = DP) :: ef0
  !! Fermi energy
  REAL(KIND = DP) :: dosef
  !! density of states at the Fermi energy
  !
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
  REAL(KIND = DP), ALLOCATABLE :: a2fij(:, :, :, :, :)
  !! spectral function a2fij(nkfs_pool,nqftot,nbndfs,nbndfs,nqstep)
  REAL(KIND = DP), ALLOCATABLE :: w0g(:, :)
  !! approximation for delta function w0g(nbndfs,nkfs)
  REAL(KIND = DP), ALLOCATABLE :: Agap(:, :, :)
  !! superconducting gap edge Agap(nkfs,nbndfs,nstemp)
  REAL(KIND = DP), ALLOCATABLE :: ADeltai(:, :, :)
  !! gap function on imag-axis at iw, ADeltai(nbndfs,nkfs,nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: ADeltaip(:, :, :)
  !! gap function on imag-axis at iwp, ADeltaip(nbndfs,nkfs,nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: AZnormi(:, :, :)
  !! renormalization function on imag-axis at iw, AZnormi(nbndfs,nkfs,nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: NAZnormi(:, :, :)
  !! normal state renormalization function on imag-axis at iw, NAZnormi(nbndfs,nkfs,nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: AKeri(:, :, :, :, :)
  !! phonon kernel on imag-axis, AKeri(nkfs,nqftot,nbndfs,nbndfs,2*nsiw(nstemp))
  REAL(KIND = DP), ALLOCATABLE :: ADsumi(:, :, :)
  !! contribution to Delta eqn from the imaginary-axis in the analytic continuation ADsumi(nbndfs,nkfs,nsw)
  REAL(KIND = DP), ALLOCATABLE :: AZsumi(:, :, :)
  !! contribution to Znorm eqn from the imaginary-axis in the analytic continuation AZsumi(nbndfs,nkfs,nsw)
  REAL(KIND = DP), ALLOCATABLE :: memlt_pool(:)
  !! maximum allocatable memory per pool
  COMPLEX(KIND = DP), ALLOCATABLE :: AZnorm(:, :, :)
  !! renormalization function on real-axis AZnorm(nbndfs,nkfs,nsw)
  COMPLEX(KIND = DP), ALLOCATABLE :: AZnormp(:, :, :)
  !! renormalization function on real-axis AZnormkq(nbndfs,nkfs,nsw)
  COMPLEX(KIND = DP), ALLOCATABLE :: ADelta(:, :, :)
  !! gap function on real-axis ADelta(nbndfs,nkfs,nsw)
  COMPLEX(KIND = DP), ALLOCATABLE :: ADeltap(:, :, :)
  !! gap function on real-axis ADeltap(nbndfs,nkfs,nsw)
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
