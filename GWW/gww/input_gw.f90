!P.Umari program GWW
! Modified by G. Stenuit
!
  MODULE input_gw
!thsi module contain the definition and routines to read
!input parameters from file
   USE kinds,  ONLY: DP

   TYPE input_options
!structure defining all the inputs required
      INTEGER :: n!number of time intervals in the positive or negative range, total of 2n+1 samples
      REAL(kind=DP) :: tau!imaginary positive time interval
      LOGICAL :: whole_s=.false.!if .true. also off-diagonal elements of self energy are calculated
      INTEGER :: max_i!maximum state to be calculated
      CHARACTER(len=256) :: prefix!prefix to designate the files same as in PW
      INTEGER :: num_rows=50!number of rows of polarization to be taken together for FFT
      LOGICAL :: xc_together=.false.!if true exchange and correlation treated together, usually NOT
      LOGICAL :: debug=.false.!if .true. debug output is considered
!the following options are for the fit
      INTEGER :: n_multipoles=2!number of poles used in the fit
      REAL(kind=DP) :: fit_dt!delta for verlet
      REAL(kind=DP) :: fit_thres=1.d-180!threshold for convergence
      REAL(kind=DP) :: fit_ma_0!mass paremeter for verlet
      REAL(kind=DP) :: fit_ma!mass paremeter for verlet
      REAL(kind=DP) :: fit_mb!mass paremeter for verlet
      REAL(kind=DP) :: fit_frice!frice for verlet
      INTEGER :: fit_maxiter=1000 !max number of iterations
      INTEGER :: starting_point=1!defines the starting point: 1-all 2-from polarization 3-from W 4-from Sigma
                               !5-from fit... HAS BEEN CHANGED A BIT
      INTEGER :: ending_point=7!defines the ending point for faster remainder calculation
      LOGICAL :: use_contractions=.true.!if true perform contraction for calculating the self energy
      LOGICAL :: lnonorthogonal=.false.!if true a non orthogonal basis set is considered
      LOGICAL :: l_hf_energies=.false.! if true uses perturbative HF energies for calculating G and W
      INTEGER :: n_fit!number of time intervals in the positive or negative range to be used for fit
      INTEGER :: remainder=0!if 0 not remainder if 1,2 our remainders , 3 Chelikowsky COH remainder, 4 Chelikowsky
                         !reaminder calculated in PW
      LOGICAL :: lconduction=.true.!if true accurate calculation of self_energy for conduction states
      LOGICAL :: l_remainder_cutoff=.false.!if true uses the reduced based for the remainder description
      LOGICAL :: l_contraction_single_state=.true.!if true write the contractions on disk state by state for saving memory
!the following options are for the treatment of time/frequency with user defined grids

      LOGICAL :: l_fft_timefreq=.false.!if true performs usual fft treatment of time/frequency
      INTEGER :: grid_time=2!0 = Gauss Legendre  1 = Gauss Laguerre 2 Gauss Legendre 0-T
      INTEGER :: grid_freq=2!0 = Gauss Legendre  1 = Gauss Laguerre 2 Gauss Legendre 0-T
      REAL(kind=DP) :: omega=20!max frequency if required
! the following for fitting tails
      LOGICAL :: l_fourier_fit_time=.false.!if true fits the tails in time
      LOGICAL :: l_fourier_fit_freq=.false.!if true fits the tails in frequency
      REAL(kind=DP) :: r_tau=1.3!ratio for finding outer time point
      REAL(kind=DP) :: r_omega=2.0!ratio for finding outer frequency point
      REAL(kind=DP) :: g_tau=10.0!ratio for treating bad cases in time
      REAL(kind=DP) :: g_omega=0.9!ration for treating bad cases in frequency
!the following for defining the grid for fitting the self energy

      INTEGER :: grid_fit=1!grid for self energy ON FREQUENCY: 0 uses the same as for P,W, 1 equally spaced, 2 GL
      REAL(kind=DP) :: omega_fit=20!max frequency to be considered
      INTEGER :: n_grid_fit=240!number of grid points on half-axes
!the following defines the range of states we want to studi
      INTEGER :: i_min=1!minimum state
      INTEGER :: i_max!maximum state

      INTEGER :: n_max_minpack=10000!maximum number of minpack iterations
      INTEGER :: cyc_minpack=1!number of minpack cycles

      LOGICAL :: l_lda_hartree=.true.!if true uses LDA previously  calculated hartree, otherwise recalculate
      LOGICAL :: l_lda_exchange=.true.!if true uses LDA previously  calculated exchange, otherwise recalculate

      LOGICAL :: l_read_exchange=.false.!force to read again exchange energies from disk

      LOGICAL :: l_symm_epsilon=.true.!if true calculates the symmetrized dielectric matrix for calculating W
      LOGICAL :: l_head_epsilon=.false.!if true the head of the symmetrized dielectric matrix is taken from
                               !a RPA external calculation, OTHERWISE it's set to zero
      INTEGER :: w_divergence=0!treatment of G=0,G=0 divergence of W: 0 no treatment potentially wrong
                             !1 W=v^1/2 Epsi^-1 v^1/2 with v truncated, 2 calculates terms from external file
      LOGICAL :: l_wing_epsilon=.false.!if true the wings of the symmetrized dielectric matrix are taken from
                               !a RPA external calculation, OTHERWISE it's set to zero
      INTEGER :: offset_fit=2!offset for fit from 0

      LOGICAL :: lcprim_file=.false. !if true read S_{ c' c i} terms from file
      LOGICAL :: lvcprim_file=.true. !if truem read S_{ i, v,c i} terms from file
      LOGICAL :: l_dft_xc_file=.false.!if true read dft exchange and correlation energies from file
      LOGICAL :: lpola_file=.true.!if true calculates the polarization directly from the psi_v psi_c \tilde{w}^P_i file

      LOGICAL :: l_self_from_pola=.false.!if true calculates the self energies terms directly from the dresses polarization
                                 !and not from dress interaction NOW ONLY FOR polarization analysis, it requires
                                 !the wp_v file
      LOGICAL :: l_pola_beta=.false.!if true calculate the polarization with beta strategy (optimized)
      LOGICAL :: l_self_beta=.false.!if true calculate the self-energy with beta strategy (optimized)

      LOGICAL :: l_pola_upper=.false.!if true uses reduced conduction states for calculating the polarization
      LOGICAL :: l_self_upper=.false.!if true uses reduced conduction states for calculating the self-energy
      INTEGER :: nc_minus=0 !number of conduction states not to be considered for alpha-strategy

   END TYPE  input_options

  CONTAINS

    SUBROUTINE  read_input_gww( ggwin )
!this subroutines reads the input structure  from file
!a namelist scheme is used
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE mp,                   ONLY : mp_bcast

    implicit none

    TYPE(input_options) :: ggwin!in output the input parameters
    INTEGER :: iun

    NAMELIST/inputgww/ ggwin

!read namelist

    if(ionode) then
      CALL input_from_file()
      WRITE(stdout, '(5x,a)') "Waiting for input..."
      READ(5, NML=inputgww)
      WRITE(stdout, '(5x,a)') "done"


!writes on screen

       write(stdout,*) 'Number of intervals: ', ggwin%n
       write(stdout,*) 'Number of intervals for fit:', ggwin%n_fit
       write(stdout,*) 'Maximum imaginary time: ',ggwin%tau
       write(stdout,*) 'Print whole Sigma: ', ggwin%whole_s
       write(stdout,*) 'Maximum state considered:', ggwin%max_i
       write(stdout,*) 'Prefix:', trim(ggwin%prefix)
       write(stdout,*) 'Rows together in FFT:',ggwin%num_rows
       if(ggwin%use_contractions) write(stdout,*) 'Uses contraction'
       if(ggwin%lnonorthogonal) write(stdout,*) 'Ultralocalized generalized Wanniers are assumed'
       if(ggwin%l_hf_energies) write(stdout, *) ' Uses perturbative HF energies for G and W'
       if(ggwin%remainder == 1) write(stdout, *) ' Uses remainder for self-energy'
       if(ggwin%lconduction) write(stdout, *) 'Accurate treatment of conduction states'
       if(ggwin%l_remainder_cutoff) write(stdout,*) 'Uses reduced set for remainder'
       if(ggwin%l_contraction_single_state) write(stdout,*) 'Uses contractions on single states'
       if(.not.ggwin%l_fft_timefreq) then
          write(stdout,*) 'Uses no fft grid for time/space integrations'
          write(stdout,*) 'Time grid :', ggwin%grid_time
          write(stdout,*) 'Frequency grid :', ggwin%grid_freq
          write(stdout,*) 'Max frequency :', ggwin%omega
          if(ggwin%l_fourier_fit_time) then
             write(stdout,*) 'Uses fit of long tails in time'
             write(stdout,*) 'Ratio: Time', ggwin%r_tau
             write(stdout,*) 'Ratio recover: Time', ggwin%g_tau
          endif
          if(ggwin%l_fourier_fit_freq) then
             write(stdout,*) 'Uses fit of long tails in frequency'
             write(stdout,*) 'Ratio: Freq', ggwin%r_omega
             write(stdout,*) 'Ratio recover: Freq', ggwin%g_omega
          endif
       endif
       write(stdout,*) 'Consider states from:', ggwin%i_min, ' to ', ggwin%i_max
       write(stdout,*) 'Maximum number of iterations in minpack: ', ggwin%n_max_minpack
       write(stdout,*) 'Number of cycles over minpack: ', ggwin%cyc_minpack
       write(stdout,*) 'Uses LDA Hartree :', ggwin%l_lda_hartree
       write(stdout,*) 'Uses LDA Exchange :', ggwin%l_lda_exchange
       if(ggwin%l_read_exchange) write(stdout,*) 'ReREAD EXCHANGE FROM DISK'
       if(ggwin%l_symm_epsilon) write(stdout,*) 'Uses SYMMETRIZED DIELECTRIC MATRIX'
       if(ggwin%l_head_epsilon) write(stdout,*) 'Head of dielectric matrix from file'
       write(stdout,*) 'Treatment of W divergence:', ggwin%w_divergence
       if(ggwin%l_wing_epsilon) write(stdout,*) 'Wings of dielectric matrix from file'
       write(stdout,*) 'Offset fit :', ggwin%offset_fit
       if(ggwin%lcprim_file) write(stdout,*) 'Read cprim terms from file'
       if(ggwin%lvcprim_file) write(stdout,*) 'Read vcprim terms from file'
       if(ggwin%l_dft_xc_file) write(stdout,*) 'Read DFT XC from file'
       if(ggwin%lpola_file)  write(stdout,*) 'Calculates the polarization directly from the overlaps'
       if(ggwin%l_self_from_pola) write(stdout,*) 'Calculates Self Energy from Dresses Polarization'
       if(ggwin%l_pola_beta) write(stdout,*) 'Calculate polarization with BETA strategy'
       if(ggwin%l_self_beta) write(stdout,*) 'Calculate self-energy with BETA strategy'
       if(ggwin%l_pola_upper) write(stdout,*) 'Reduced Upper states for Polarization'
       if(ggwin%l_self_upper) write(stdout,*) 'Reduced Upper states for Self-Energy'
       write(stdout,*) 'Number of excluded conduction states', ggwin%nc_minus
    endif

#ifdef __MPI

    call mp_bcast(ggwin%n,ionode_id)
    call mp_bcast(ggwin%tau,ionode_id)
    call mp_bcast(ggwin%whole_s,ionode_id)
    call mp_bcast(ggwin%max_i,ionode_id)
    call mp_bcast(ggwin%prefix,ionode_id)
    call mp_bcast(ggwin%num_rows,ionode_id)
    call mp_bcast(ggwin%xc_together,ionode_id)
    call mp_bcast(ggwin%debug, ionode_id)
    call mp_bcast(ggwin%n_multipoles,ionode_id)
    call mp_bcast(ggwin%fit_dt,ionode_id)
    call mp_bcast(ggwin%fit_thres,ionode_id)
    call mp_bcast(ggwin%fit_ma_0,ionode_id)
    call mp_bcast(ggwin%fit_ma,ionode_id)
    call mp_bcast(ggwin%fit_mb,ionode_id)
    call mp_bcast(ggwin%fit_frice,ionode_id)
    call mp_bcast(ggwin%fit_maxiter,ionode_id)
    call mp_bcast(ggwin%starting_point, ionode_id)
    call mp_bcast(ggwin%ending_point, ionode_id)
    call mp_bcast(ggwin%use_contractions,ionode_id)
    call mp_bcast(ggwin%lnonorthogonal, ionode_id)
    call mp_bcast(ggwin%n_fit, ionode_id)
    call mp_bcast(ggwin%l_hf_energies, ionode_id)
    call mp_bcast(ggwin%remainder, ionode_id)
    call mp_bcast(ggwin%lconduction, ionode_id)
    call mp_bcast(ggwin%l_remainder_cutoff, ionode_id)
    call mp_bcast(ggwin%l_contraction_single_state, ionode_id)
    call mp_bcast(ggwin%l_fft_timefreq, ionode_id)
    call mp_bcast(ggwin%grid_time, ionode_id)
    call mp_bcast(ggwin%grid_freq, ionode_id)
    call mp_bcast(ggwin%omega, ionode_id)
    call mp_bcast(ggwin%l_fourier_fit_time, ionode_id)
    call mp_bcast(ggwin%l_fourier_fit_freq, ionode_id)
    call mp_bcast(ggwin%r_tau, ionode_id)
    call mp_bcast(ggwin%r_omega, ionode_id)
    call mp_bcast(ggwin%g_tau, ionode_id)
    call mp_bcast(ggwin%g_omega, ionode_id)
    call mp_bcast(ggwin%grid_fit, ionode_id)
    call mp_bcast(ggwin%omega_fit, ionode_id)
    call mp_bcast(ggwin%n_grid_fit, ionode_id)
    call mp_bcast(ggwin%i_min, ionode_id)
    call mp_bcast(ggwin%i_max, ionode_id)
    call mp_bcast(ggwin%n_max_minpack, ionode_id)
    call mp_bcast(ggwin%cyc_minpack, ionode_id)
    call mp_bcast(ggwin%l_lda_hartree, ionode_id)
    call mp_bcast(ggwin%l_lda_exchange, ionode_id)
    call mp_bcast(ggwin%l_read_exchange, ionode_id)
    call mp_bcast(ggwin%l_symm_epsilon, ionode_id)
    call mp_bcast(ggwin%l_head_epsilon, ionode_id)
    call mp_bcast(ggwin%w_divergence, ionode_id)
    call mp_bcast(ggwin%l_wing_epsilon, ionode_id)
    call mp_bcast(ggwin%offset_fit, ionode_id)
    call mp_bcast(ggwin%lcprim_file, ionode_id)
    call mp_bcast(ggwin%lvcprim_file, ionode_id)
    call mp_bcast(ggwin%l_dft_xc_file, ionode_id)
    call mp_bcast(ggwin%lpola_file, ionode_id)
    call mp_bcast(ggwin%l_self_from_pola, ionode_id)
    call mp_bcast(ggwin%l_pola_beta, ionode_id)
    call mp_bcast(ggwin%l_self_beta, ionode_id)
    call mp_bcast(ggwin%l_pola_upper, ionode_id)
    call mp_bcast(ggwin%l_self_upper, ionode_id)
    call mp_bcast(ggwin%nc_minus, ionode_id)
#endif

    return
  END SUBROUTINE



  END MODULE

