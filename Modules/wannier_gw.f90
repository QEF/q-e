!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE wannier_gw
  !
  ! ... The variables needed for gww-gwl code (head.x)
  !
  USE kinds, ONLY: DP
  !
  SAVE
  !
  !From HEAD
  LOGICAL :: l_head=.false.!if true calculates the head of the symmetrized dielectric matrix -1
  INTEGER :: n_gauss!number of frequency steps for head calculation
  REAL(kind=DP) :: omega_gauss!period for frequency calculation
  INTEGER :: grid_type!0 GL -T,T 2 GL 0 T 3 Equally spaced 0 Omega
  INTEGER :: nsteps_lanczos!number of lanczos steps
    !options for grid_freq=5
  INTEGER :: second_grid_n!sub spacing for second grid
  INTEGER :: second_grid_i!max regular step using the second grid
  LOGICAL :: l_scissor!if true displaces occupied manifold of scissor(1) and unoccupied manifold of scissor(2)
  REAL(kind=DP) :: scissor(2)!see above
  !From pw4gww
  TYPE real_matrix_pointer
     REAL(kind=DP), DIMENSION(:,:), POINTER :: p
  END TYPE real_matrix_pointer

  TYPE complex_matrix_pointer
     COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: p
  END TYPE complex_matrix_pointer


  TYPE optimal_options!options for obtaining optimal basis sets
     LOGICAL :: l_complete!if true just do a diagonalization
     INTEGER :: idiago !kind of optimization: 0=Gram-Schmidt like
     INTEGER :: ithres!kind of threshold: 0=on modulus square
     REAL(kind=DP) :: thres!value of threshold
  END TYPE optimal_options


 
  REAL(kind=DP),  ALLOCATABLE :: wannier_centers(:,:,:)!wannier centers in a.u.
  REAL(kind=DP),  ALLOCATABLE :: wannier_radii(:,:)!wannier centers in a.u.
  INTEGER, ALLOCATABLE  :: w_centers(:,:,:)!wanier centers on the grid

  INTEGER, ALLOCATABLE  :: w_radii(:,:)!wannier lengths in grid units
  COMPLEX(kind=DP), ALLOCATABLE :: u_trans(:,:,:)!unitarian transformation from bloch wfcs to wannier'
  INTEGER :: numw_prod!number of products w_i(r)*w_j(r) then of orthonormalized products
  INTEGER :: num_nbndv(2) !number of valence bands
  INTEGER :: num_nbnds !number of studied bands valence plus  a part of conduction's
  REAL(kind=DP), ALLOCATABLE :: becp_gw(:,:,:)!to store projections of wfcs with us projectors
  REAL(kind=DP), ALLOCATABLE :: becp_gw_c(:,:,:)!to store projections of wfcs with us projectors for {c'} subspace
  COMPLEX(kind=DP), ALLOCATABLE :: expgsave(:,:,:,:) !to store exp_igx  on us augmentation functions
  INTEGER :: nset!number of states to be read  written from/to file simultaneously
  LOGICAL :: l_truncated_coulomb!if true the Coulomb potential is truncated
  REAL(kind=DP) :: truncation_radius!truncation radius for Coulomb potential
  INTEGER :: remainder!1-cutoff 2-distance 3-no remainder 4-postprocessing from W 5-postprocessing from dressed polarization P
  INTEGER :: restart_gww!for restarting the calculation of gww stuff, 0 begins from beginning


  LOGICAL :: l_gram!if true uses gram schmidt for orthonormalizing the products of wanniers
  !LOGICAL :: l_head!if true calculates the head of the symmetrized dielectric matrix -1
  !INTEGER :: n_gauss!number of frequency steps for head calculation
  !REAL(kind=DP) :: omega_gauss!period for frequency calculation
  LOGICAL :: l_exchange!if true calculate the exchange terms with k-points sampling


  LOGICAL :: l_zero!if .true. calculate also the v e v^1/2 operators with G=0,G'=0 put to 0

  LOGICAL :: l_wing!if .true. calculate also the wing terms, it requires the file .e_head


  !INTEGER :: grid_type!0 GL -T,T 2 GL 0 T 3 Equally spaced 0 Omega


  INTEGER :: nset_overlap!number of states to be read  written from/to file simultaneously, when
                         !calculating overlaps
  INTEGER :: nspace!space on grid for evalueation of exchange-type integrals


  REAL(kind=DP) :: ecutoff_global!cut off in Rydbergs for G sum on (dense charge grid)


  INTEGER :: maxiter2!max number of iteration for the genaralized maximally localized wannier
                      !of the second conduction manifold
  REAL(kind=DP) :: diago_thr2!thresold for electronic states used in c_bands for upper
                              !conduction manifold if any, if ==0 used same cutoff as for valence
  LOGICAL :: l_plot_mlwf!if true save the orthonormal wannier for plotting



  INTEGER :: max_ngm!max number of g vector for charge grid effctively stored

!variables for parallelization on matrices

  LOGICAL :: l_pmatrix !if true parallelize on  matrices
  INTEGER :: p_mpime!processor number
  INTEGER :: p_nproc!number of processors
  INTEGER :: npcol!number of processor columns
  INTEGER :: nprow!number of processor rows
  INTEGER :: icontxt!blacs descriptor
  INTEGER :: myrow!actual processor row
  INTEGER :: mycol!actual processor column



  LOGICAL :: l_coulomb_analysis!if true after polarization analysis consider eigenvalues of coulomb potential
  REAL(kind=DP) ::  cutoff_coulomb_analysis!cutoff for coulomb analysis


  INTEGER :: n_pola_lanczos!number of orthonormal states for polarization lanczos-style
  INTEGER :: n_self_lanczos!number of orthonormal states for self-energy lanczos-style
  INTEGER :: nsteps_lanczos_pola!number of lanczos steps for the polarizability
  INTEGER :: nsteps_lanczos_self!number of lanczos steps for the self_energy
  REAL(kind=DP) :: s_pola_lanczos!cutoff for lanczos basis for polarization
  REAL(kind=DP) :: s_self_lanczos!cutoff for lanczos basis for self-energy
  INTEGER :: nump_lanczos!dimension of basis for lanczos calculation of the polarization
  INTEGER :: nums_lanczos!dimension of basis for lanczos calculation of the self-energy
  REAL(kind=DP) :: s_g_lanczos!cutoff for absolute value of trial green function

  LOGICAL :: l_pmat_diago!if true find the basis for the polarization diagonalizing the O matrix
  REAL(kind=DP) :: pmat_ethr!threshold for diagonalizing the O matrix

  REAL(kind=DP) :: pmat_cutoff!cutoff (in Ryd) for polarization diagonalization
  INTEGER :: pmat_type!type of approximation 1 usual, 2 with wanniers, 3 with optimal representation,5 just plane waves
  INTEGER :: n_fast_pmat!number of states  for fast evaluation of conduction manifold if =0 disabled
  INTEGER :: n_pmat!number of orthonormal states for optimal representation O matrix
  REAL(kind=DP) :: s_pmat!cutoff for optimal basis for O matrix
  INTEGER :: lanczos_restart!restart point for lanczos
  INTEGER :: n_pola_lanczos_eff!effective number of pola states; if 0 equal to n_pola_lanczos
  INTEGER :: n_self_lanczos_eff!effective number of self states; if 0 equal to n_self_lanczos
  REAL(kind=DP) :: off_fast_pmat!offset in Ry for fast assessment of polarizability if =0 disabled
  LOGICAL :: l_fast_pola!if true fast assessment of polarizability for basis construction  
  LOGICAL :: l_v_basis!if true valuate the polarizability basis vectors as eigenstates of v operator
  REAL(kind=DP) :: v_cutoff!cutoff in Ryd for v operator
  LOGICAL :: l_iter_algorithm!if true uses iterative algorithms
  REAL(kind=DP) :: dual_pb!dual value till 4.d0 for defing the grid on which the polarizability basis is created
  
  REAL(kind=DP), ALLOCATABLE :: vg_q(:) ! contains the elements V(G) of the Coulomb potential obtained upon integration over q
  LOGICAL :: l_t_wannier!if true builds t verctors starting from KS valence wannier functions
  REAL(kind=DP) :: dual_vt!dual value till 4.d0 for defing the grid on which the t vectors created
  REAL(kind=DP) :: dual_vs!dual value till 4.d0 for defing the grid on which the s vectors created

  LOGICAL  :: lwannier!if true take advantage of localization of wannier functions
  REAL(kind=DP) :: wannier_thres!threshold for modulus of wannier function in a.u.
  

  INTEGER :: s_first_state!if different from 0, first KS state for calculatin s vectors (if last 1)
  INTEGER :: s_last_state!if different from 0, last KS state for calculatin s vectors (if last num_nbnds)

  LOGICAL :: l_selfconsistent!if true do selfconsistent GW calculation, requires file band.dat
  REAL(kind=DP), ALLOCATABLE :: ene_gw(:,:)!GW energies of previous iteration for selfconsistent calculation
  INTEGER :: n_gw_states!number of GW states for selfconsistent calculation
  REAL(kind=DP) :: delta_self!delta energy for selfconsistent calculation

  LOGICAL :: l_whole_s!if true calculates also the off-diagonal elemenets of V_xc for then 
                      !calculating the off-diagonal elements of sigma 
  LOGICAL :: l_ts_eigen!if true the t and global vectors are calculated considering also the eigenvalues of the partial basis (recommanded)

  LOGICAL :: l_frac_occ! if true consider fractional occupancies
  INTEGER :: num_nbndv_min(2)!limits for fully occupied states

  LOGICAL :: l_cond_pol_base!if true uses conduction states till num_nbnds for the construction of the polarizability bases


  LOGICAL :: l_semicore!if true evaluate semicore terms
  INTEGER :: n_semicore!number of semicore states staring from the bottom of valence states
  LOGICAL :: l_semicore_read!if true reads semicore file for calculating products for Green's function

  LOGICAL :: l_verbose!if true a lot of ouput for debug


  LOGICAL :: l_contour! if true calculates the terms for contour integration
  LOGICAL :: l_real!if true calculate the polarizability basis, s and t vectors avoiding ffts it requires more memory

  LOGICAL :: l_big_system!if true uses startegy for large systems: just local s vectors are used

  REAL(kind=DP) ::extra_pw_cutoff!cutoff to add to the optimal polarizability basis plane-waves (sin and cos functions) 
                                 !if 0 no plane waves is added
  
  !REAL(kind=DP) :: exchange_fast_dual!for defining the fast exchnage routines

  LOGICAL :: l_bse!if true computing quantities for bse calculation
  REAL(kind=DP) :: s_bse!threshold for wannier function overlap
  REAL(kind=DP) :: dual_bse!dual factor for bse calculations

  LOGICAL :: l_simple!if true writes on disk polarizability basis on real space for further post-processing

  LOGICAL :: l_list !if true uses startegy for large systems from list of states included in s_first_state, s_last_state
  INTEGER :: n_list(2)!number of states in list for the 2 spin channels
  INTEGER, ALLOCATABLE :: i_list(:,:) !list of KS states to be computed 

  LOGICAL :: l_full!if true prepare data for further post-processing for a full-relativistic calculation
  INTEGER :: n_full(2)!numeber of proper relativistic states in G of GW for collinear spin channel

!variables for splitting the head calculation in blocks in order to reduce memory usage
!no effect on results

  INTEGER :: len_head_block_freq!length of blocks on frequency
  INTEGER :: len_head_block_wfc!length of blocks on unperturbed occupied wfcs 

  

  INTERFACE free_memory

  MODULE PROCEDURE free_complex,free_real

  END INTERFACE

  CONTAINS
   
    subroutine free_complex( c)
      implicit none
      type(complex_matrix_pointer) :: c
      deallocate(c%p)
      return
    end subroutine


   subroutine free_real( r)
      implicit none
      type(real_matrix_pointer) :: r
      deallocate(r%p)
      return
    end subroutine


      
    subroutine  max_ngm_set
 !set the value of max_ngm
      use io_global, only : stdout
      use gvect, only :     ngm,gg
      use cell_base, only : tpiba2
      
      implicit none

      integer :: ig

      max_ngm=0
      do ig=1,ngm
         if(gg(ig)*tpiba2 >= ecutoff_global) exit
         max_ngm=max_ngm+1
      enddo

      write(stdout,*) 'MAX_NGM:', max_ngm, ngm

end subroutine max_ngm_set
          
 

END MODULE wannier_gw
