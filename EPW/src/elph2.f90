  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  MODULE elph2
  !-------------------------------------------------------------------------- 
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       el_ph_mat  (:,:,:,:),   &!  e-p matrix  (nbnd, nbnd, nks, 3*nat)
       cu(:,:,:),              &!  rot matrix for wannier interpolation of k point, coarse mesh (nbnd*nbnd*nkstot)
       cuq(:,:,:),             &!  rot matrix for wannier interpolation of k+q point, coarse mesh (nbnd*nbnd*nkstot)
       chw (:,:,:),            &!  Hamiltonian in wannier basis 
       chw_ks (:,:,:),         &!  Hamiltonian in wannier basis  (Kohn-Sham eigenvalues if many-body eigenvalues are read in)
       cdmew (:,:,:,:),        &!  Dipole matrix in wannier basis 
       cvmew (:,:,:,:),        &!  Velocity matrix in wannier basis 
       rdw (:,:,:),            &!  dynamical matrix in wannier basis (real) 
       epmatwp (:,:,:,:,:),    &!  e-p matrix  in wannier basis - electrons and phonons
       umat(:,:,:),            &!  the rotation matrix for the unique setting of the wfs gauge -- on the local pool
       umatq(:,:,:),           &!  the rotation matrix for the unique setting of the wfs gauge for the k + q-- on the local pool
       umat_all(:,:,:),        &!  the rotation matrix for the unique setting of the wfs gauge -- for all k points
       umatq_all(:,:,:),       &!  the rotation matrix for the unique setting of the wfs gauge -- for all k+q points
       dynq  (:,:,:),          &!  dynamical matrix for every q (nmode, nmodes, nqtot)
       epmatq (:,:,:,:,:),     &!  e-p matrix for every q (nbnd, nbnd, nks, nmodes, nqtot)
       epf17 (:, :, :, :),  &!  full ep matrix in bloch rep stored in mem (nkqtotf, nbnd, nbnd, nmodes)-nbnd inside wndw 
       dmec(:,:,:,:),          &!  dipole matrix elements on the coarse mesh (ipol, nbnd, nbnd, nks)
       dmef(:,:,:,:),          &!  dipole matrix elements on the fine   mesh (ipol, nbnd, nbnd, nks)
       vmef(:,:,:,:),          &!  velocity matrix elements on the fine mesh (ipol, nbnd, nbnd, nks)
       bmat(:,:,:,:)            !   overlap U_k+q U_k^\dagger on the coarse mesh (nbnd, nbnd, nks, nqtot)
  REAL(KIND=DP), ALLOCATABLE ::&
       a_all(:,:),             &!  electronic spectral function du to electron-phonon interaction
       xk_all(:,:),            &!  full k point grid, coarse (3, nkstot)
       et_all(:,:),            &!  full eigenvalue list, coarse (nbnd, nkstot)
       et_ks(:,:),             &!  lda eigenvalues
       et_mb(:,:),             &!  gw eigenvalues
       xkq(:,:),               &!  local k+q grid, coarse (3, nks)
       etq(:,:),               &!  eigenvalues of k+q wavefunctions
       xkf(:,:),               &!  fine k point grid (3, nkqf)
       wkf(:),                 &!  weights on the fine grid (nkqf)
       xqf(:,:),               &!  fine q point grid 
       wqf(:),                 &!  weights on the fine q grid 
       etf(:,:),               &!  interpolated eigenvalues (nbnd, nkqf)
       etf_k(:,:),             &!  Saved interpolated KS eigenenergies for later used in q-parallelization (nbnd, nkqf)
       etf_ks(:,:),            &!  interpolated eigenvalues (nbnd, nkqf) KS eigenvalues in the case of eig_read
       wf(:,:),                &!  interpolated eigenfrequencies 
       wslen(:),               &!  length of the wigner seitz points in units of alat
       gamma_all(:,:,:),       &!
       gamma_nest(:,:),        &!  Nesting function in the case of q-parallelization
       gamma_v_all(:,:,:),     &!
       lambda_all(:,:,:),      &!
       lambda_v_all(:,:,:),    &!
       sigmar_all(:,:),        &!  To store sigmar, sigmai and zi globally
       sigmai_all(:,:),        &!
       sigmai_mode(:,:,:),     &! 
       zi_all(:,:),            &!
       esigmar_all(:,:,:),     &!
       esigmai_all(:,:,:),     &!   
       jdos(:),                &!
       spectra(:,:,:,:,:,:),   &!  dipole absorption spectra, polarizations, nomega, nsmear, dme/vme, absorption/emission
       sumr(:,:,:,:),          &!  to apply the ASR correction to every xq
       zstar(:,:,:),           &!  Born effective charges
       epsi(:,:),              &!  dielectric tensor
       inv_tau(:,:,:)           !  scattering rate
  REAL(KIND=DP) ::             &!
       efnew                    !  SP: Fermi level on the fine grid. Added globaly for efficiency reason 
  INTEGER ::                   &!
       nkqf,                   &!  number of k+q points per pool (fine grid)
       nkf,                    &!  number of k points per pool (fine grid)
       nqf,                    &!  number of q points per pool (fine grid)
       nkqtotf,                &!  total number of k+q points (fine grid)
       nqtotf,                 &!  total number of q points (fine grid)
       nrr,                    &!  number of wigner-seitz points (elec interp only)
       nrr_k,                  &!  number of wigner-seitz points for electrons
       nrr_q,                  &!  number of wigner-seitz points for phonons
       ibndmin,                &!  band bounds for slimming down electron-phonon matrix 
       ibndmax                  !
  INTEGER, ALLOCATABLE ::      & 
       irvec(:,:),             &!  crys coordinates of wigner-seitz vectors (both elec and phon)
       ndegen(:),              &!  corresponding degeneragy, electrons (old version)
       ndegen_k(:),            &!  corresponding degeneragy, electrons
       ndegen_q(:),            &!  corresponding degeneragy, phonons
       igk(:),                 &!  Index for k+G vector
       igkq(:),                &!  Index for k+q+G vector
       igk_k_all(:,:),         &!  Global index (in case of parallel)
       ngk_all(:)               !  Global number of plane wave for each global k-point
  INTEGER, allocatable ::      &
       shift (:),              &!  for every k+q, index of the G0 which folds k+q into k+q+G0 of the first BZ
       gmap(:)                  !  the map G -> G-G_0 in the large (density) G vectors set, for every G_0
  LOGICAL, allocatable ::      &
       lwin(:,:),              &!  identify bands within outer energy windows (when disentanglement is used)
       lwinq(:,:),             &!
       done_elph(:)
  LOGICAL ::                   &
       elph                   
  END MODULE elph2
