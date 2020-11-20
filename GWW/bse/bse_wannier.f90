!this module contains bse code input variables  
MODULE bse_wannier

  USE kinds, ONLY: DP

  INTEGER :: numw_prod!number of products w_i(r)*w_j(r) then of orthonormalized products
  INTEGER :: num_nbndv(2) !number of valence bands
  INTEGER :: nset!number of states to be read  written from/to file simultaneously
  INTEGER :: maxit!maximun number of iterarion for conj gradient minimiaztion
  INTEGER :: n_eig!number of excitoninc eigenvalues to be found
  INTEGER :: cg_nreset!conj gradient variable number of steps after which 
                      !the residual is effectively computed with the application
                      ! of the Hamiltonian  
  INTEGER :: plotn_min!start plotting from plot_nmin eigenstate 
  INTEGER :: plotn_max!end plotting at plot_nmax eigenstate
  INTEGER :: spectra_nstep! number of frequencies for which the spectrum will be
                          !computed 
  INTEGER :: l_restart! if 1 restart from the computation of the n_eig_start
                      ! vector; if 2 restart from the calculation of the
                      ! absorption spectrum
  INTEGER :: n_eig_start! eigenvector index where the calculation restarts from 
  INTEGER :: nlcz_restart! number of lanczos iteration where the calculation
                            !restarts from
  INTEGER :: nit_lcz! number of total lanczos iteration to be performed 
  INTEGER :: qpe_imin ! index of the lowest band for which GW energies are known
  INTEGER :: qpe_imax ! index of the highest band for which GW energies are known 
  LOGICAL :: l_truncated_coulomb!if true the Coulomb potential is truncated
  LOGICAL :: l_verbose!if true higher verbosity output
  LOGICAL :: l_cgrad! if true conjgradient diago
  LOGICAL :: l_plotexc !if true plot the excitonic wavefunctions
  LOGICAL :: l_plotaverage !if true plot the average excitonic wavefunction (for degenerate states)
  LOGICAL :: l_tspace !if true solve the BSE in transition space
  LOGICAL :: l_finite !if true computes the spectrum using the expectation
                       !value of the position operator, makes sense only for finite
                       !systems
  LOGICAL :: l_lanczos !if true computes absorption through a lanczos procedure
  LOGICAL :: l_restart_lcz! if true restart the lanczos calculation from
                          ! previous run
  LOGICAL :: l_fullbse ! if true perform full BSE calculation
  LOGICAL :: l_tdhf ! if true perform do not include Wc,i.e. perform
                    ! a td-hf like calculation
  LOGICAL :: l_lf   ! do not include the W=v+Wc term local fields only
  LOGICAL :: l_rpa  ! RPA non-local field calculation
  LOGICAL :: l_scissor! if true use scissor operator, else use computed QP energies
  REAL(kind=DP) :: truncation_radius!truncation radius for Coulomb potential
  REAL(kind=DP) :: dual_bse!dual factor for bse calculations
  REAL(kind=DP) :: lambda!small positive real for steepest descendent
  REAL(kind=DP) :: eps!small positive real for steepest descendent convergence
  REAL(kind=DP) :: eps_eig!small positive real for steepest descendent convergence check on eigenstate
  REAL(kind=DP) :: lm_delta! magnitude of the step along the search direction in
                           ! the line minimization subroutine  
  REAL(kind=DP) :: scissor!scissor operator opening the KS gap in eV
  REAL(kind=DP) :: r_hole(3)!coordinate of the hole (in alat units) for the excitonic wavefunction plot
  REAL(kind=DP) :: r_pola(3)!polarization direction
  REAL(kind=DP) :: spectra_e_min!minimum energy for the absorption spectrum
  REAL(kind=DP) :: spectra_e_max!max energy for the absorption spectrum
  REAL(kind=DP) :: spectra_broad!range of broadening for the spectra (in eV) 
  LOGICAL :: l_contraction!if true uses contraction for speeding-up
  LOGICAL :: l_gtrick!save arrays in G space for limiting memory usage
  LOGICAL :: l_dielectric!if true in case of l_lanczos=.false. calculates the dielectric function from excitonic amplitudes

END MODULE bse_wannier

