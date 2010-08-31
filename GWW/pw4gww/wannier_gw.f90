! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit
!
!#ifdef __GWW
!
!new module for wannier function support
MODULE wannier_gw

  USE kinds, ONLY: DP

  SAVE

  TYPE wannier_product!this structure described the product of 2 wannier functions w_i(r)*w_j(r)
     INTEGER :: i! i
     INTEGER :: j! j
     INTEGER :: center(3)!center of the product on the grid
     INTEGER :: radius(3)! radius length in grid units on three cartesian directions
  END TYPE wannier_product

  TYPE real_matrix_pointer
     REAL(kind=DP), DIMENSION(:,:), POINTER :: p
  END TYPE real_matrix_pointer

  TYPE complex_matrix_pointer
     COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: p
  END TYPE complex_matrix_pointer

  TYPE wannier_P!this structure described the localized and normalized products of wanniers  w_P
     INTEGER :: numij!number of (unordered)couples of wannier w_i w_j which overlap both with w_P
     INTEGER, DIMENSION(:,:), POINTER :: ij!array for i,j (ij,:)
     REAL(kind=DP),DIMENSION(:), POINTER :: o!overlap <w_P|w_i*w_j>
  END TYPE wannier_P

  LOGICAL  :: lwannier!if true calculates wannier functions
  REAL(kind=DP),  ALLOCATABLE :: wannier_centers(:,:)!wannier centers in a.u.
  REAL(kind=DP),  ALLOCATABLE :: wannier_radii(:)!wannier centers in a.u.
  INTEGER, ALLOCATABLE  :: w_centers(:,:)!wanier centers on the grid

  INTEGER, ALLOCATABLE  :: w_radii(:)!wannier lengths in grid units
  INTEGER, ALLOCATABLE  :: w_centers_c(:,:)!wanier centers on the grid for {c'} subspace
  INTEGER, ALLOCATABLE  :: w_radii_c(:)!wannier lengths in grid units for {c'} subspace
  COMPLEX(kind=DP), ALLOCATABLE :: u_trans(:,:)!unitarian transformation from bloch wfcs to wannier'
  REAL(kind=DP) :: cutoff_wsq!cutoff for |w|^2
  REAL(kind=DP) :: cutoff_wsq_c!cutoff for |w|^2 for conduction states
  REAL(kind=DP) :: cutoff_wpr!cutoff for w_i*w_j
  REAL(kind=DP) :: cutoff_overlap!for overlaps <w^{p}|w_iw_j>
  INTEGER :: numw_prod!number of products w_i(r)*w_j(r) then of orthonormalized products
  INTEGER :: numw_prod_c!number of products w_c'(r)*w_c(r)
  INTEGER :: numw_prod_vvc!!number of products w_i(r)*w_j(r) original
  INTEGER :: numw_prod_val_cond!number of products w_val(r)*w_cond(r) original
  INTEGER :: numw_prod_val_cond_sec!number of products w_val(r)*w_cond^2(r) original
  TYPE(wannier_product), POINTER :: w_prod(:)!array for describing products of wannier wfcs
  INTEGER :: num_nbndv !number of valence bands
  INTEGER :: num_nbnds !number of studied bands valence plus  a part of conduction's
  LOGICAL :: lsmallgrid!uses small wavefunction G grid, instead of denser R grid for wannier products
  LOGICAL :: lnonorthogonal!if true non orthogonal ultralocalized wanniers
  REAL(kind=DP) :: no_radius!radius to be used for nonorthogonal localization
  REAL(kind=DP), ALLOCATABLE :: becp_gw(:,:)!to store projections of wfcs with us projectors
  REAL(kind=DP), ALLOCATABLE :: becp_gw_c(:,:)!to store projections of wfcs with us projectors for {c'} subspace
  LOGICAL :: lggrid!if true calculates overlaps in g space
  COMPLEX(kind=DP), ALLOCATABLE :: expgsave(:,:,:,:) !to store exp_igx  on us augmentation functions
  INTEGER :: nset!number of states to be read  written from/to file simultaneously
  REAL(kind=DP) :: ultra_alpha_v!alpha factor for valence states
  REAL(kind=DP) :: ultra_alpha_c!alpha factor for conduction states
  REAL(kind=DP) :: ultra_alpha_c2!alpha factor for conduction states for second block
  REAL(kind=DP) :: ultra_alpha_c_prim!alpha factor for the subspace of observed conduction states
  INTEGER :: num_nbndc_set!number of conduction bands to be treated as a separate manifold, for sigma calculation
  LOGICAL :: l_truncated_coulomb!if true the Coulomb potential is truncated
  REAL(kind=DP) :: truncation_radius!truncation radius for Coulomb potential
  INTEGER :: numw_prodprod!total number of products of products of non-orthogonal wanniers
  REAL(kind=DP) :: cutoff_wpr_wpr!cutoff for  products of products of non-orthogonal wanniers
  REAL(kind=DP) :: r_cutoff_products!cutoff radius for distance of products of wanniers
  LOGICAL, ALLOCATABLE :: l_on_products(:,:)!if true the distance of two products of wanniers is less than r_cutoff_products
  INTEGER :: remainder!1-cutoff 2-distance 3-no remainder 4-postprocessing from W 5-postprocessing from dressed polarization P
  INTEGER :: restart_gww!for restarting the calculation of gww stuff, 0 begins from beginning

  REAL(kind=DP) :: cutoff_wpr_vc!cutoff for w_i*w_j for V-C first block
  REAL(kind=DP) :: cutoff_wpr_vc2!cutoff for w_i*w_j for V-C second block
  INTEGER :: num_nbnd_first!defines first block of conduction states

  REAL(kind=DP) :: cutoff_wpr_prim!cutoff for w_i*w_j for C-C' first block
  REAL(kind=DP) :: cutoff_wpr_prim2!cutoff for w_i*w_j for C-C' second block

  LOGICAL :: l_gram!if true uses gram schmidt for orthonormalizing the products of wanniers

  LOGICAL :: l_head!if true calculates the head of the symmetrized dielectric matrix -1
  INTEGER :: n_gauss!number of frequency steps for head calculation
  REAL(kind=DP) :: omega_gauss!period for frequency calculation
  LOGICAL :: l_exchange!if true calculate the exchange terms with k-points sampling
  REAL(kind=DP) :: tau_gauss!period for the calculation of the gauss legendre time grid

  LOGICAL :: l_zero!if .true. calculate also the v e v^1/2 operators with G=0,G'=0 put to 0

  LOGICAL :: l_wing!if .true. calculate also the wing terms, it requires the file .e_head


  INTEGER :: grid_type!0 GL -T,T 2 GL 0 T

  INTEGER :: cprim_type!if == 0 treats c' manifold through ultraorthogonalized wannier, if ==1 calculates all
                       !elemenst S_{c'c i}, if ==2 calculates S{c' v,c i}
  INTEGER :: cprim_first,cprim_last!define the range for c'

  LOGICAL :: l_vcw_overlap!if true calculates the  overlaps \int dr \Psi_v(r)\Psi_c(r) \tilde{w}^P_i

  INTEGER :: nset_overlap!number of states to be read  written from/to file simultaneously, when
                         !calculating overlaps
  INTEGER :: nspace!space on grid for evalueation of exchange-type integrals

  REAL(kind=DP) :: lambda_ene!if >0.  localizes wanniers also on energy
  REAL(kind=DP) :: e_min_cutoff!energy for Ev-Ec for which the cutoff is /= 0 (in Ry)
  REAL(kind=DP) :: e_max_cutoff!energy for Ev-Ec for which the cutoff is /= 0 and /= +inf (in Ry)
  REAL(kind=DP) :: v_min_cutoff!value of cutoff at e_min_cutoff (in a.u.)
  REAL(kind=DP) :: v_max_cutoff!value of cutoff at e_max_cutoff (in a.u.)

  LOGICAL :: l_orthonorm_products!if true orthonormalize directly the products of wannier
  LOGICAL :: l_wpwt_terms!if .true. (default) calculates  the overlaps <\tilde{w}_i| w_j>

  REAL(kind=DP) :: cutoff_products!cutoff for orthonormalization of wannier products

  REAL(kind=DP) :: ecutoff_global!cut off in Rydbergs for G sum on (dense charge grid)


  LOGICAL :: l_polarization_analysis!if true calculate the polarization and retain only significant eigenvectors
  REAL(kind=DP) :: cutoff_polarization!cutoff for polarization analysis
  INTEGER :: nc_polarization_analysis!number of conduction states for calculating the polarization
                                           !in polarization analysis
  LOGICAL :: l_only_val_cond!if true one calculating orthonormalization of wannier products
                            !considers only valence*conduction products
  LOGICAL :: l_no_val_cond_sec!if true one calculating orthonormalization of wannier products
                            !considers only valence*(second block of conduction) products

  INTEGER :: maxiter2!max number of iteration for the genaralized maximally localized wannier
                      !of the second conduction manifold
  REAL(kind=DP) :: diago_thr2!thresold for electronic states used in c_bands for upper
                              !conduction manifold if any, if ==0 used same cutoff as for valence
  LOGICAL :: l_plot_mlwf!if true save the orthonormal wannier for plotting
  LOGICAL :: l_plot_ulwf!if true save the ultralocalized wannier for plotting

  LOGICAL :: l_ultra_external!if true ultralocalize through the new procedure

  INTEGER :: nbnd_normal!max number of bands with normal treatment if == 0 sets nbnd_max = nbnd
  INTEGER :: num_nbnd_delta!number of upper bands for which the same energy is considered
  INTEGER :: num_nbnd_upper!number of upper bands AFTER THE REDUCTION for which the same energy is considered


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

  LOGICAL :: l_assume_ortho!during s_3 step assumes orthonormality (and uses s_2 cutoff)

  LOGICAL :: l_coulomb_analysis!if true after polarization analysis consider eigenvalues of coulomb potential
  REAL(kind=DP) ::  cutoff_coulomb_analysis!cutoff for coulomb analysis
  REAL(kind=DP) :: mem_per_core! memory (in Bytes) per core (not per CPU or per node !)

  INTERFACE free_memory

  MODULE PROCEDURE free_complex,free_real,free_wannier_P

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


   subroutine free_wannier_P( w_P)
      implicit none
      type(wannier_P) :: w_P
      deallocate(w_P%ij)
      deallocate(w_P%o)
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

!#endif

