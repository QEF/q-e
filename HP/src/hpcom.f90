!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE ldaU_hp
  !
  ! Common variables for the HP program
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  SAVE
  !
  LOGICAL :: skip_type(ntypx),        &     ! If .true. skip the calculation for a specific type 
                                            ! (e.g. Ni_up can be used for Ni_down with some spin 
                                            ! considerations)
             perturb_only_atom(500),  &     ! If perturb_only_atom(i)=.true. perterb only i-th atom
             recalc_sym,              &     ! If .true. we recalculate the number of symmetries
                                            ! of the unperturbed lattice due to the change of 
                                            ! the atomic type of one of the atoms
             compute_hp,              &     ! If .true. collects all pieces of chi0 and chi
             sum_pertq,               &     ! If .true. collects dns0 and dnsscf for all q points
                                            ! (for the specific perturbed atom) and computes their 
                                            ! sum with the phase factor
             determine_num_pert_only, &     ! If .true. determine only which atoms must be perterbed
             skip_equivalence_q,      &     ! If .true. the full frid of q points will be used
             disable_type_analysis,   &     ! If .true. disable the algorithm which detects whether
                                            ! there are atoms of the same type but with different occupations
             skip_atom(500),          &     ! If .true. no LR calculation will be performed 
                                            ! for a selected atomic site.
                                            ! skip_atom(i), where i runs over atoms. If skip_atom(i)=.true.
                                            ! then no linear-response calculation will be performed for the
                                            ! i-th atom. This keyword cannot be used when find_atpert=1.
                                            ! Warning: Make sure you know what you are doing! This option might
                                            ! be useful in several cases:
                                            ! - Debugging purposes;
                                            ! - You know that the atom which you do not want to perturb
                                            !   is equivalent to some other atom (but the code does not recognizes
                                            !   this from the symmetry analysis). In this case check that there is
                                            !   at least one atom of the same type which was perturbed (this can
                                            !   happen only when find_atpert=3), otherwise the post-processing
                                            !   calculation of U will fail.
             search_sym                     ! If .TRUE. search for the symmetry of q
  !
  LOGICAL, ALLOCATABLE :: todo_atom(:),              &      ! Which atoms must be perturbed
                          perturbed_atom(:),         & ! Controls which atom is perturbed in the HP
                                                       ! calculation
                          this_pert_is_on_file(:),   & ! The perturbation is written on file or not
                          comp_iq(:)                   ! If .true. this q point has to be calculated
  !
  INTEGER :: nath,            &             ! Number of (real) atoms in the primitive cell
                                            ! with Hubbard_U \= 0
             nath_sc,         &             ! Total number of real+virtual Hubbard atoms
                                            ! in the virtual supercell
             nqsh,            &             ! Number of q points in the grid (without symmetry 
                                            ! reduction) = number of primitive cells in the
                                            ! virtual supercell
             nah_pert,        &             ! Site number of the perturbed Hubbard atom
             nath_pert,       &             ! Number of actual perturbed Hubbard atoms in the primitive cell
             find_atpert,     &             ! Method of searching for atoms to be perturbed
             iter_best,       &             ! Number of the iteration at which the best accuracy
                                            ! for chi was reached
             ntyp_new,        &             ! Maximum number of different types detected in the calculation
                                            ! (used only when find_atpert=3)
             num_neigh,       &             ! Used in the postprocessing: number of nearest neighbors of every atom 
                                            ! which will be written to the file parameters.out 
                                            ! (can be used only with lda_plus_u_kind = 2)      
             lmin,            &             ! Used in the postprocessing: minimum value of the orbital quantum number
                                            ! of the Hubbard atoms (Hubbard_l) starting from which (and up to the maximum
                                            ! Hubbard_l in the system) Hubbard V will be written to the file parameters.out  
             niter_max,       &             ! Maximum number of iterations to solve the linear system
             nmix,            &             ! Number of iterations used in potential mixing
             nq1, nq2, nq3,   &             ! Number of q points in each direction
             nqs,             &             ! Number of q points to be calculated
             start_q, last_q, &             ! Initial and final q in the list
             iudwfc,          &             ! Unit for response wavefunctions  
             lrdwfc,          &             ! Length of the record for response wavefunctions
             iudvwfc,         &             ! Unit for the perturbing potential * wavefunctions
             lrdvwfc                        ! Length of the record for the perturbing potential * wavefunctions
  !
  INTEGER :: equiv_type(ntypx)              ! equiv_type(i)=j, will merge type i to type j 
                                            ! (useful when nspin=2)
  !
  CHARACTER(LEN=16)  :: background          ! Background correction
  CHARACTER(LEN=256) :: tmp_dir_save,  &    ! Temprary directory  
                        tmp_dir_hp                  
  CHARACTER(LEN=5)   :: atm_new(500)        ! Label for atomic types
  CHARACTER(LEN=2)   :: code = 'HP'         ! Name of the code
  !
  INTEGER, PARAMETER :: maxter = 500        ! Maximum number of iterations
  !
  REAL(DP) :: conv_thr_chi,      &          ! Threshold for the calculation of chi 
              thresh_init,       &          ! Initial threshold for the solution of 
                                            ! the linear system
              conv_thr_chi_best, &          ! The best accuracy which was reached for chi 
                                            ! (useful for the analysis when the convergence 
                                            ! is not reached)
              docc_thr,          &          ! Threshold for the comparison of the unperturbed 
                                            ! occupations (used only with find_atpert=1 for
                                            ! determination of atoms which must be perturbed)
              rmax,              &          ! Maximum distance (in Bohr) between two atoms 
                                            ! to search for neighbors (used only at the 
                                            ! postprocessing step when lda_plus_u_kind = 2).
              alpha_mix(maxter)             ! The mixing parameter
  !
  REAL(DP), ALLOCATABLE :: ns(:),      &    ! Trace of unperturbed occupations (spin up + spin down)
                           magn(:),    &    ! Unperturbed magnetization
                           Rvect(:,:), &    ! Radius-vector of the primitive cell
                           chi0(:,:),  &    ! Bare response function (from 1st iteration)
                           chi(:,:),   &    ! SCF response function
                           x_q(:,:)         ! Coordinates of q points
  !
  COMPLEX(DP), ALLOCATABLE :: dns0(:,:,:,:,:),         & ! Bare response occupation matrix 
                                                         ! (from 1st iteration)
                              dnsscf(:,:,:,:,:),       & ! SCF response occupation matrix 
                              dns0_tot(:,:,:,:,:),     & ! Total bare response occupation matrix 
                                                         ! (summed over q)
                              dnsscf_tot(:,:,:,:,:),   & ! Total SCF  response occupation matrix 
                                                         ! (summed over q)
                              trace_dns_tot_old(:)       ! Trace of the response occupation matrix 
                                                         ! (for a convergence test) 
  !
  INTEGER, ALLOCATABLE :: ityp_new(:)       ! Types of atoms
  !
END MODULE ldaU_hp
