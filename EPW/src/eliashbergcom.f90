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
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! Global variables for real- and imag-axis Eliashberg equations
  !
  INTEGER :: nsw
  !
  ! nsw : nr. of grid points between (0,wscut) for real-axis, analytical
  !       continuation and Pade approximants
  !
  INTEGER, ALLOCATABLE :: nsiw(:)
  !
  ! nsiw(:) : nr of grid points at each temperature on imag-axis, nsiw(nstemp)
  !
  REAL(DP) :: wsphmax, dwsph, gap0
  !
  ! wsphmax  : maximum phonon frequency for evaluation of the integral over Omega (0,wsphmax)
  ! dwsph : frequency step for Eliashberg spectral function 
  ! gap0  : initial guess for Delta
  !
  REAL(DP), ALLOCATABLE :: dws(:), ws(:), wsph(:), wsi(:), estemp(:)
  !
  ! dws    : grid size at each bin dws(nsw)
  ! ws     : frequency on real-axis, ws(nsw)
  ! wsph   : frequency on real-axis, wsph(nqstep)
  ! wsi(:) : frequency on imag-axis at iw, wi(nsiw(nstemp))
  ! estemp : temperature in eV entering in the Eliashberg equtions estemp(nstemp)
  !
  END MODULE eliashberg_common
  !
  !--------------------------------------------------------------------------
  MODULE eliashberg_common_iso
  !--------------------------------------------------------------------------
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! Global variables for real and imag-axis isotropic equations Eliashberg equations
  !
  REAL(DP), ALLOCATABLE :: a2f_iso(:), gap(:), fdwp(:), bewph(:)
  !
  ! a2f_iso: isotropic Eliashberg spectral function a2f_iso(nqstep)
  ! gap   : superconducting gap edge gap(nstemp)
  ! fdwp  : Fermi-Dirac distribution at frequency wp, fdwp(nsw)
  ! bewph : Boise-Einstein distribution at frequency wph, bewph(nqstep)
  !
  REAL(DP), ALLOCATABLE :: Deltai(:), Deltaip(:), Znormi(:), NZnormi(:), Keri(:), Dsumi(:), Zsumi(:)
  !
  ! Deltai  : gap function on imag-axis at iw, Deltai(nsiw(nstemp))
  ! Deltaip : gap function on imag-axis at iwp, Deltaip(nsiw(nstemp))
  ! Znormi  : renormalization function on the imag-axis at iw, Znormi(nsiw(nstemp))
  ! NZnormi : normal state renormalization function on the imag-axis at iw, NZnormi(nsiw(nstemp))
  ! Keri : phonon kernel on imag-axis, Keri(2*nsiw(nstemp))
  ! Dsumi : contribution to Delta eqn from the imaginary-axis in the analytic continuation Dsumi(nsw)
  ! Zsumi : contribution to Znorm eqn from the imaginary-axis in the analytic continuation Zsumi(nsw)
  !
  REAL(DP), ALLOCATABLE :: Gp(:,:), Gm(:,:)
  !
  ! Gp(nsw,nqstep)  : -bose(omegap)-fermi( omega+omegap) (eqn for Delta and Znorm analytic continuation)
  ! Gm(nsw,nqstep)  :  bose(omegap)+fermi(-omega+omegap) (eqn for Delta and Znorm analytic continuation)
  !
  COMPLEX(DP), ALLOCATABLE :: Delta(:), Deltap(:), Znorm(:), Znormp(:), Kp(:,:), Km(:,:)
  !
  ! Delta(nsw) : gap function on real-axis at iw
  ! Deltap(nsw): gap function on real-axis at iw
  ! Znorm(nsw) : renormalization function on real-axis at iw
  ! Znormp(nsw): renormalization function on real-axis at iw
  ! Kp(nsw,nsw) : phonon kernel on real-axis (eqn for Delta)
  ! Km(nsw,nsw) : phonon kernel on real-axis (eqn for Znorm)
  !
  END MODULE eliashberg_common_iso
  !
  !--------------------------------------------------------------------------
  MODULE eliashberg_common_aniso
  !--------------------------------------------------------------------------
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! Global variables for real and imag-axis anisotropic equations Eliashberg equations
  !
  LOGICAL :: limag_fly, lacon_fly
  !
  INTEGER :: nkfs, nbndfs
  !
  ! nkfs : nr. of irreducible k-points within the Fermi shell on the fine mesh
  ! nbndfs  : nr. of electronic bands within the Fermi shell
  !
  INTEGER, ALLOCATABLE :: equivk(:), ixkff(:), ixkf(:), ixkqf(:,:), ixqfs(:,:), nqfs(:)
  !
  ! nkf = nr of irreducible k-points on the fine grid, if mp_mesh_k = .true.
  ! nkf = total nr of k-points on the fine grid,       otherwise
  ! equivk : index of equivalent k-points on the fine k-grid equivk(nkf)
  ! ixkff : index of k-point on the full k-grid ixkff(nkftot)
  ! ixkf : index of k-point on the irreducible k-grid within the Fermi shell ixkf(nkf)
  ! ixkqf : index k+q or k-q on the irreducilble k-grid within the Fermi shell ixkqf(nkfs,nqftot)
  ! ixqfs : index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell 
  !         ixqfs(nkfs,nqfs(ik))
  ! nqfs : nr of q-points at each k-point for which k+sign*q is within the Fermi shell nqfs(nkfs)
  !
  REAL(DP) :: ef0, dosef
  !
  ! ef0     : Fermi energy 
  ! dosef   : density of states at the Fermi energy
  !
  REAL(DP), ALLOCATABLE :: g2(:,:,:,:,:), ekfs(:,:), xkff(:,:), xkfs(:,:), wkfs(:), & 
                           a2fij(:,:,:,:,:), w0g(:,:), Agap(:,:,:)
  !
  ! g2      : e-ph matrix element squared |g_ji^nu(k,q)|^2, g2(nkfs_pool,nqftot,nbndfs,nbndfs,nmodes)
  ! ekfs : eigenvalues at E_i(k), etf(nbndfs,nkfs)
  ! xkff : coordintates of the k-points of the full k-grid xkff(3,nkftot)
  ! xkfs : coordintates of the k-points xkf(3,nkfs)
  ! wkfs : weights of the irreducible k-points wkf(nkfs)
  ! a2fij  : spectral function a2fij(nkfs_pool,nqftot,nbndfs,nbndfs,nqstep) 
  ! w0g    : approximation for delta function w0g(nbndfs,nkfs)
  ! Agap   : superconducting gap edge Agap(nkfs,nbndfs,nstemp)
  !
  REAL(DP), ALLOCATABLE :: ADeltai(:,:,:), ADeltaip(:,:,:), AZnormi(:,:,:), NAZnormi(:,:,:), & 
                           AKeri(:,:,:,:,:), ADsumi(:,:,:), AZsumi(:,:,:)
  !
  ! ADeltai  : gap function on imag-axis at iw, ADeltai(nbndfs,nkfs,nsiw(nstemp))
  ! ADeltaip : gap function on imag-axis at iwp, ADeltaip(nbndfs,nkfs,nsiw(nstemp))
  ! AZnormi  : renormalization function on imag-axis at iw, AZnormi(nbndfs,nkfs,nsiw(nstemp))
  ! NAZnormi : normal state renormalization function on imag-axis at iw, NAZnormi(nbndfs,nkfs,nsiw(nstemp))
  ! AKeri : phonon kernel on imag-axis, AKeri(nkfs,nqftot,nbndfs,nbndfs,2*nsiw(nstemp))
  ! ADsumi : contribution to Delta eqn from the imaginary-axis in the analytic continuation ADsumi(nbndfs,nkfs,nsw)
  ! AZsumi : contribution to Znorm eqn from the imaginary-axis in the analytic continuation AZsumi(nbndfs,nkfs,nsw)
  !
  COMPLEX(DP), ALLOCATABLE :: AZnorm(:,:,:), AZnormp(:,:,:), ADelta(:,:,:), ADeltap(:,:,:)
  !
  ! AZnorm   : renormalization function on real-axis AZnorm(nbndfs,nkfs,nsw)
  ! AZnormp  : renormalization function on real-axis AZnormkq(nbndfs,nkfs,nsw)
  ! ADelta   : gap function on real-axis ADelta(nbndfs,nkfs,nsw)
  ! ADeltap  : gap function on real-axis ADeltap(nbndfs,nkfs,nsw)
  !
  REAL(DP), ALLOCATABLE :: memlt_pool(:)
  !
  ! memlt_pool : maximum allocatable memory per pool
  !
  END MODULE eliashberg_common_aniso
  !
  !--------------------------------------------------------------------------
  MODULE eliashbergcom
  !-------------------------------------------------------------------------- 
  !
  USE eliashberg_common
  USE eliashberg_common_iso
  USE eliashberg_common_aniso
  !
  END MODULE eliashbergcom

