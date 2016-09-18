!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


  MODULE input_gw
!thsi module contain the definition and routines to read
!input parameters from file
   USE kinds,  ONLY: DP

   TYPE input_options
!structure defining all the inputs required
      INTEGER :: n!number of time intervals in the positive or negative range, total of 2n+1 samples
      REAL(kind=DP) :: tau=0.d0!imaginary positive time interval
      LOGICAL :: whole_s=.false.!if .true. also off-diagonal elements of self energy are calculated
      INTEGER :: max_i!maximum state to be calculated
      CHARACTER(len=256) :: prefix = 'prefix'!prefix to designate the files same as in PW
      CHARACTER(len=256) :: outdir = './'!outdir to designate the files same as in PW
      INTEGER :: num_rows=50!number of rows of polarization to be taken together for FFT
      LOGICAL :: xc_together=.false.!if true exchange and correlation treated together, usually NOT
      LOGICAL :: debug=.false.!if .true. debug output is considered
!the following options are for the fit
      INTEGER :: n_multipoles=2!number of poles used in the fit
      REAL(kind=DP) :: fit_dt!delta for verlet
      REAL(kind=DP) :: fit_thres=1.d-5!threshold for convergence
      REAL(kind=DP) :: fit_ma_0!mass paremeter for verlet
      REAL(kind=DP) :: fit_ma!mass paremeter for verlet
      REAL(kind=DP) :: fit_mb!mass paremeter for verlet
      REAL(kind=DP) :: fit_frice!frice for verlet
      INTEGER :: fit_maxiter=5000 !max number of iterations
      INTEGER :: starting_point=1!defines the starting point: 1-all 2-from polarization 3-from W 4-from Sigma
                               !5-from fit... HAS BEEN CHANGED A BIT
      INTEGER :: ending_point=7!defines the ending point 
      LOGICAL :: use_contractions=.false.!if true perform contraction for calculating the self energy
      LOGICAL :: lnonorthogonal=.false.!if true a non orthogonal basis set is considered
      LOGICAL :: l_hf_energies=.false.! if true uses perturbative HF energies for calculating G and W
      INTEGER :: n_fit!number of time intervals in the positive or negative range to be used for fit
      LOGICAL :: lconduction=.true.!if true accurate calculation of self_energy for conduction states
      LOGICAL :: l_contraction_single_state=.true.!if true write the contractions on disk state by state for saving memory
!the following options are for the treatment of time/frequency with user defined grids

      LOGICAL :: l_fft_timefreq=.false.!if true performs usual fft treatment of time/frequency
      INTEGER :: grid_time=3!0 = Gauss Legendre  1 = Gauss Laguerre 2 Gauss Legendre 0-T 3 = equally spaced 0 centered 4 = equally spaced 
      INTEGER :: grid_freq=3!0 = Gauss Legendre  1 = Gauss Laguerre 2 Gauss Legendre 0-T 3 = equally spaced 0 centered 4 = equally spaced 5=double equally spaced grid
      REAL(kind=DP) :: omega!max frequency if required
! the following for fitting tails
      LOGICAL :: l_fourier_fit_time=.false.!if true fits the tails in time
      LOGICAL :: l_fourier_fit_freq=.false.!if true fits the tails in frequency
      REAL(kind=DP) :: r_tau=0.d0!ratio for finding outer time point
      REAL(kind=DP) :: r_omega=0.d0!ratio for finding outer frequency point
      REAL(kind=DP) :: g_tau=0.d0!ratio for treating bad cases in time
      REAL(kind=DP) :: g_omega=0.d0!ration for treating bad cases in frequency
!the following for defining the grid for fitting the self energy

      INTEGER :: grid_fit=1!grid for self energy ON FREQUENCY: 0 uses the same as for P,W, 1 equally spaced, 2 GL
      REAL(kind=DP) :: omega_fit!max frequency to be considered
      INTEGER :: n_grid_fit!number of grid points on half-axes
!the following defines the range of states we want to studi
      INTEGER :: i_min!minimum state
      INTEGER :: i_max!maximum state
      
      INTEGER :: n_max_minpack=20000!maximum number of minpack iterations
      INTEGER :: cyc_minpack=1!number of minpack cycles

      LOGICAL :: l_lda_hartree=.true.!if true uses LDA previously  calculated hartree, otherwise recalculate
      LOGICAL :: l_lda_exchange=.true.!if true uses LDA previously  calculated exchange, otherwise recalculate

      LOGICAL :: l_read_exchange=.false.!force to read again exchange energies from disk

      LOGICAL :: l_symm_epsilon=.true.!if true calculates the symmetrized dielectric matrix for calculating W
      LOGICAL :: l_head_epsilon=.true.!if true the head of the symmetrized dielectric matrix is taken from
                               !a RPA external calculation, OTHERWISE it's set to zero
      INTEGER :: w_divergence!treatment of G=0,G=0 divergence of W: 0 no treatment potentially wrong
                             !1 W=v^1/2 Epsi^-1 v^1/2 with v truncated, 2 calculates terms from external file
                             !3 Lanczos chain and Gygi Baldereschi algorithm
      LOGICAL :: l_wing_epsilon=.true.!if true the wings of the symmetrized dielectric matrix are taken from
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

      LOGICAL :: l_pola_lanczos=.true.!if true calculate the polarization through a lanczos scheme
      LOGICAL :: l_self_lanczos=.true.!if true calculate the self-energy through a lanczos scheme
      LOGICAL :: l_lanczos_conv=.true.!if the convolution is done analytically
      INTEGER :: n_set=100!block length over frequency in do_self_lanczos
      INTEGER :: n_set_ii=5!block length over KS states in do_self_lanczos
      INTEGER :: n_set_pola=4!block length over valence states in do_pola_lanczos
      INTEGER :: n_set_self=1!block length over frequency states in do_self_lanczos
      LOGICAL :: l_yet_pola=.false.!if true it assumes that (restarting at point 6) the 
                           !polaw's contain already the dress polarization
      LOGICAL :: l_reduce_io=.true. !if true reduces disk I/O
!options for grid_freq=5
      INTEGER :: second_grid_n=10!sub spacing for second grid
      INTEGER :: second_grid_i=1!max regular step using the second grid
      INTEGER :: grid_levels=1!for grid type = 6
      LOGICAL :: l_t_wannier=.true.!if true t vectors have been constructed from wannier valence functions
      LOGICAL :: l_truncated_coulomb=.true.!if true the system is finite otherwise is extended (with head and wings)
      LOGICAL :: l_self_time=.true.!if true calculates the self energy in imaginary time

      LOGICAL :: l_g_grid=.false.!if true use a dedicated grid on frequency for G 
      INTEGER :: grid_freq_g=3!for G grid:0 =  Gauss Legendre  1 = Gauss Laguerre 2 Gauss Legendre 0-T 3 = equally spaced 0 centered 4 = equally spaced 5=double spaced grid
      INTEGER :: n_g!for G grid: number of time intervals in the positive or negative range, total of 2n+1 samples
      INTEGER :: second_grid_n_g=10!for G grid: sub spacing for second grid
      INTEGER :: second_grid_i_g=1!for G grid: max regular step using the second grid
      REAL(kind=DP) :: omega_g! for G: max frequency

      INTEGER :: i_min_whole=0!set the minimum of range for off diagonal elements of self-energy
      INTEGER :: i_max_whole=0!set the maximum of range for off diagonal elements of self-energy
      INTEGER :: nspin=1!spin multiplicity

      LOGICAL :: l_frac_occ=.false.!if true uses routines for fractional occupancy
      LOGICAL :: l_semicore=.false.!if true add to the self energy semicore terms
      
      LOGICAL :: l_order=.false.!just for gww_fit if true from file order.dat takes the right order of the self_energies

      LOGICAL :: l_verbose=.false.
      REAL(kind=DP) :: real_energy_min=-1.0!minimum energy on real frequency axis
      REAL(kind=DP) :: real_energy_max=1.0!maximum energy on real frequency axis
      INTEGER :: n_real_axis=0!number of grid points on real frequency axis 
      LOGICAL :: l_big_system=.false.!if true calculate the self-energy state  by state only through local s vectors
      LOGICAL :: l_list=.false.!if true uses startegy for large systems from list of states included in i_min im_max
      LOGICAL :: l_full=.false.!if true in points 6 calls routine for full relativistic calculations EXPERTS only
      INTEGER :: n_full=0!number of KS states with explicit treatment in G

      LOGICAL :: l_scissor=.false.!if true shifts DFT occupied bands of scissor
      REAL(kind=DP) :: scissor(2)=0.d0!value for scissor of valence manifold (1) and conduction manifold (2) in eV
      

   END TYPE  input_options

  CONTAINS

    SUBROUTINE  read_input_gww( ggwin )
!this subroutines reads the input structure  from file
!a namelist scheme is used
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
    USE io_files,             ONLY :  tmp_dir,prefix
   
    implicit none

    CHARACTER(LEN=256), EXTERNAL :: trimcheck

    INTEGER, EXTERNAL :: find_free_unit

    TYPE(input_options) :: ggwin!in output the input parameters
    INTEGER :: iun

    NAMELIST/inputgww/ ggwin

!read namelist

    if(ionode) then
       iun=find_free_unit()
       !open(iun,file='inputgww')
       !read(iun, NML=inputgww)
       !close(iun)
       read(*, NML=inputgww)
!OBM: file handling in a more QE manner
     IF ( TRIM( ggwin%outdir ) == ' ' ) ggwin%outdir = './'
     tmp_dir = trimcheck(ggwin%outdir)
     prefix = trim(ggwin%prefix)
!set up parameter for calculation with Lanczos scheme
     
       if(ggwin%l_self_lanczos.and.ggwin%l_lanczos_conv.and. .not. ggwin%l_self_time) then
          ggwin%grid_fit=1
          ggwin%omega_fit=ggwin%omega
          ggwin%n_grid_fit=ggwin%n
       endif
!writes on screen
       
       write(stdout,*) 'Number of intervals: ', ggwin%n
       write(stdout,*) 'Number of intervals for fit:', ggwin%n_fit
       if(ggwin%tau==0.d0) ggwin%tau=2.d0/ggwin%omega*dble(ggwin%n)
       write(stdout,*) 'Maximum imaginary time: ',ggwin%tau
       write(stdout,*) 'Print whole Sigma: ', ggwin%whole_s
       write(stdout,*) 'Maximum state considered:', ggwin%max_i
       write(stdout,*) 'Prefix:', trim(ggwin%prefix)
       write(stdout,*) 'Rows together in FFT:',ggwin%num_rows
       if(ggwin%use_contractions) write(stdout,*) 'Uses contraction'
       if(ggwin%lnonorthogonal) write(stdout,*) 'Ultralocalized generalized Wanniers are assumed'
       if(ggwin%l_hf_energies) write(stdout, *) ' Uses perturbative HF energies for G and W'
       if(ggwin%lconduction) write(stdout, *) 'Accurate treatment of conduction states'
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
       if(ggwin%l_pola_lanczos) write(stdout,*) 'Polarization Calculated through Lanczos scheme'
       if(ggwin%l_self_lanczos) write(stdout,*) 'Self-energy Calculated through Lanczos scheme'
       if(ggwin%l_lanczos_conv) write(stdout,*) 'Convolution done analytically'
       if(ggwin%l_self_lanczos) write(stdout,*) 'Block length fequency', ggwin%n_set
       if(ggwin%l_self_lanczos) write(stdout,*) 'Block length states', ggwin%n_set_ii
       if(ggwin%l_pola_lanczos) write(stdout,*) 'Block length valence states', ggwin%n_set_pola
       if(ggwin%l_self_lanczos) write(stdout,*) 'Block length fequency lc', ggwin%n_set_self
       if(ggwin%starting_point < 6 ) then
          ggwin%l_yet_pola=.false.
       else
          if(ggwin%l_self_lanczos) write(stdout,*) 'Dressed polarization already calculated'
       endif
       if(ggwin%l_reduce_io) write(stdout,*) 'Reduced disk I/O'
       if(ggwin%grid_freq==5) then
          write(stdout,*) 'Uses double grid, subdivisions:', ggwin%second_grid_n,' till :', ggwin%second_grid_i
       endif
       if(ggwin%l_t_wannier) write(stdout,*) 't vectors from Wannier products'
       if(ggwin%l_truncated_coulomb) then
          write(stdout,*) 'Use truncated Coulomb interaction'
          ggwin%w_divergence=0
          ggwin%l_head_epsilon=.false.
          ggwin%l_wing_epsilon=.false.
       else
          write(stdout,*) 'Use truncated Extended interaction'
          ggwin%w_divergence=3
       endif
        if(ggwin%l_self_time) write(stdout,*) 'Calculate Self-energy through FT'
        if(ggwin%l_g_grid) write(stdout,*) 'Uses dedicated grid for G'
        if(ggwin%whole_s) then
           if(ggwin%i_min_whole==0) ggwin%i_min_whole=ggwin%i_min
           if(ggwin%i_max_whole==0) ggwin%i_max_whole=ggwin%i_max
           write(stdout,*) 'Calculate off-diagonal elements in range:', ggwin%i_min_whole,ggwin%i_max_whole
        endif
        write(stdout,*) 'Spin multiplicity:', ggwin%nspin
        write(stdout,*) 'Partiallly occpuied states:', ggwin%l_frac_occ
        if(ggwin%l_semicore) write(stdout,*) 'Add semicore terms to self-energy'
        if(ggwin%n_real_axis /=0) then
           write(stdout,*) 'N of grid points on REAL frequency:', ggwin%n_real_axis
           write(stdout,*) 'REAL frequency range:', ggwin%real_energy_min,ggwin%real_energy_max
        endif
        if(ggwin%l_big_system) write(stdout,*) 'USING ONLY LOCAL S VECTORS'
        if(ggwin%l_list) write(stdout,*) 'FROM LIST'
        if(ggwin%l_full) write(stdout,*) 'FULL RELATIVISTIC CALCULATION with:', ggwin%n_full
        if(ggwin%l_scissor) write(stdout,*) 'USE SCISSOR:', ggwin%scissor
     endif
#if defined(__MPI)
    CALL mp_bcast( tmp_dir,ionode_id, world_comm )
    CALL mp_bcast( prefix,ionode_id, world_comm )
    call mp_bcast(ggwin%n,ionode_id,world_comm)
    call mp_bcast(ggwin%tau,ionode_id,world_comm)
    call mp_bcast(ggwin%whole_s,ionode_id,world_comm)
    call mp_bcast(ggwin%max_i,ionode_id,world_comm)
    call mp_bcast(ggwin%prefix,ionode_id,world_comm)
    call mp_bcast(ggwin%num_rows,ionode_id,world_comm)
    call mp_bcast(ggwin%xc_together,ionode_id,world_comm)
    call mp_bcast(ggwin%debug, ionode_id,world_comm)
    call mp_bcast(ggwin%n_multipoles,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_dt,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_thres,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_ma_0,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_ma,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_mb,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_frice,ionode_id,world_comm)
    call mp_bcast(ggwin%fit_maxiter,ionode_id,world_comm)
    call mp_bcast(ggwin%starting_point, ionode_id,world_comm)
    call mp_bcast(ggwin%ending_point, ionode_id,world_comm)
    call mp_bcast(ggwin%use_contractions,ionode_id,world_comm)
    call mp_bcast(ggwin%lnonorthogonal, ionode_id,world_comm)
    call mp_bcast(ggwin%n_fit, ionode_id,world_comm)
    call mp_bcast(ggwin%l_hf_energies, ionode_id,world_comm)
    call mp_bcast(ggwin%lconduction, ionode_id,world_comm)
    call mp_bcast(ggwin%l_contraction_single_state, ionode_id,world_comm)
    call mp_bcast(ggwin%l_fft_timefreq, ionode_id,world_comm)
    call mp_bcast(ggwin%grid_time, ionode_id,world_comm)
    call mp_bcast(ggwin%grid_freq, ionode_id,world_comm)
    call mp_bcast(ggwin%omega, ionode_id,world_comm)
    call mp_bcast(ggwin%l_fourier_fit_time, ionode_id,world_comm)
    call mp_bcast(ggwin%l_fourier_fit_freq, ionode_id,world_comm)
    call mp_bcast(ggwin%r_tau, ionode_id,world_comm)
    call mp_bcast(ggwin%r_omega, ionode_id,world_comm)
    call mp_bcast(ggwin%g_tau, ionode_id,world_comm)
    call mp_bcast(ggwin%g_omega, ionode_id,world_comm)
    call mp_bcast(ggwin%grid_fit, ionode_id,world_comm)
    call mp_bcast(ggwin%omega_fit, ionode_id,world_comm)
    call mp_bcast(ggwin%n_grid_fit, ionode_id,world_comm)
    call mp_bcast(ggwin%i_min, ionode_id,world_comm)
    call mp_bcast(ggwin%i_max, ionode_id,world_comm)
    call mp_bcast(ggwin%n_max_minpack, ionode_id,world_comm)
    call mp_bcast(ggwin%cyc_minpack, ionode_id,world_comm)
    call mp_bcast(ggwin%l_lda_hartree, ionode_id,world_comm)
    call mp_bcast(ggwin%l_lda_exchange, ionode_id,world_comm)
    call mp_bcast(ggwin%l_read_exchange, ionode_id,world_comm)
    call mp_bcast(ggwin%l_symm_epsilon, ionode_id,world_comm)
    call mp_bcast(ggwin%l_head_epsilon, ionode_id,world_comm)
    call mp_bcast(ggwin%w_divergence, ionode_id,world_comm)
    call mp_bcast(ggwin%l_wing_epsilon, ionode_id,world_comm)
    call mp_bcast(ggwin%offset_fit, ionode_id,world_comm)
    call mp_bcast(ggwin%lcprim_file, ionode_id,world_comm)
    call mp_bcast(ggwin%lvcprim_file, ionode_id,world_comm)
    call mp_bcast(ggwin%l_dft_xc_file, ionode_id,world_comm)
    call mp_bcast(ggwin%lpola_file, ionode_id,world_comm)
    call mp_bcast(ggwin%l_self_from_pola, ionode_id,world_comm)
    call mp_bcast(ggwin%l_pola_beta, ionode_id,world_comm)
    call mp_bcast(ggwin%l_self_beta, ionode_id,world_comm)
    call mp_bcast(ggwin%l_pola_upper, ionode_id,world_comm)
    call mp_bcast(ggwin%l_self_upper, ionode_id,world_comm)
    call mp_bcast(ggwin%l_pola_lanczos, ionode_id,world_comm)
    call mp_bcast(ggwin%l_self_lanczos, ionode_id,world_comm)
    call mp_bcast(ggwin%l_lanczos_conv, ionode_id,world_comm)
    call mp_bcast(ggwin%n_set, ionode_id,world_comm)
    call mp_bcast(ggwin%n_set_ii, ionode_id,world_comm)
    call mp_bcast(ggwin%n_set_pola, ionode_id,world_comm)
    call mp_bcast(ggwin%n_set_self, ionode_id,world_comm)
    call mp_bcast(ggwin%l_yet_pola, ionode_id,world_comm)
    call mp_bcast(ggwin%l_reduce_io, ionode_id,world_comm)
    call mp_bcast(ggwin%second_grid_n, ionode_id,world_comm)
    call mp_bcast(ggwin%second_grid_i, ionode_id,world_comm)
    call mp_bcast(ggwin%l_t_wannier, ionode_id,world_comm)
    call mp_bcast(ggwin%l_truncated_coulomb, ionode_id,world_comm)
    call mp_bcast(ggwin%l_self_time, ionode_id,world_comm)
    call mp_bcast(ggwin%l_g_grid, ionode_id,world_comm)
    call mp_bcast(ggwin%grid_freq_g, ionode_id,world_comm)
    call mp_bcast(ggwin%n_g, ionode_id,world_comm)
    call mp_bcast(ggwin%second_grid_n_g, ionode_id,world_comm)
    call mp_bcast(ggwin%second_grid_i_g, ionode_id,world_comm)
    call mp_bcast(ggwin%omega_g, ionode_id,world_comm)
    call mp_bcast(ggwin%i_min_whole, ionode_id,world_comm)
    call mp_bcast(ggwin%i_max_whole, ionode_id,world_comm)
    call mp_bcast(ggwin%nspin, ionode_id,world_comm)
    call mp_bcast(ggwin%l_frac_occ, ionode_id,world_comm)
    call mp_bcast(ggwin%l_semicore, ionode_id,world_comm)
    call mp_bcast(ggwin%l_order, ionode_id,world_comm)
    call mp_bcast(ggwin%l_verbose, ionode_id,world_comm)
    call mp_bcast(ggwin%n_real_axis, ionode_id,world_comm)
    call mp_bcast(ggwin%real_energy_min, ionode_id,world_comm)
    call mp_bcast(ggwin%real_energy_max, ionode_id,world_comm)
    call mp_bcast(ggwin%grid_levels, ionode_id,world_comm)
    call mp_bcast(ggwin%l_big_system, ionode_id,world_comm)
    call mp_bcast(ggwin%l_list, ionode_id,world_comm)
    call mp_bcast(ggwin%l_full, ionode_id,world_comm)
    call mp_bcast(ggwin%n_full, ionode_id,world_comm)
    call mp_bcast(ggwin%l_scissor, ionode_id, world_comm)
    call mp_bcast(ggwin%scissor, ionode_id, world_comm)
#endif

    return
  END SUBROUTINE read_input_gww



  END MODULE

