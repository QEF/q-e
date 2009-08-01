! Copyright (C) 2006-2008 Dmitry Korotin - dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!

MODULE wannier_new
  !
  ! ... Variables to construct and store wannier functions
  !
  USE kinds,      ONLY : DP
  !
  SAVE
  ! 
  INTEGER, PARAMETER :: ningx = 10 ! max number of trial wavefunction ingredients

  LOGICAL :: &
     use_wannier,              &! if .TRUE. wannier functions are constructed
     rkmesh,                   &! if .TRUE. regular k-mesh without symmetry is used !now used in input_parameters_mod
     plot_wannier,             &! if .TRUE. wannier number plot_wan_num is plotted
     use_energy_int,           &! if .TRUE. uses energy interval for wannier generation, not band numbers   
        print_wannier_coeff           ! if .TRUE. computes and prints coefficients of wannier decomp. on atomic functions
  INTEGER :: &
     nwan,                     &! number of wannier functions
     plot_wan_num,             &! number of wannier for plotting  
     plot_wan_spin              ! spin of wannier for plotting
  REAL(kind=DP), allocatable ::  &
     wan_pot(:,:),             &! constrained potential
     wannier_energy(:,:),      &! energy of each wannier (of each spin)
     wannier_occ(:,:,:)         ! occupation matrix of wannier functions(of each spin)
  COMPLEX(kind=DP), allocatable :: &
     pp(:,:),                  &! <phi|S|psi> projections
     coef(:,:,:)                ! coefficients of wannier decomp. on atomic functions
        
  TYPE ingredient
         INTEGER :: l = 0, &           ! l value for atomic wfc
                    m = 0, &           ! m value for atomic wfc
                    iatomwfc = 0       ! number of corresponding atomic orbital
         REAL :: c = 0.d0              ! coefficient
  END TYPE ingredient
  
  TYPE wannier_data
          INTEGER :: iatom = 0,              &
                     ning = 0
          REAL ::        bands_from = 0.d0, &
                         bands_to = 0.d0
          TYPE (ingredient) :: ing(ningx)
  END TYPE wannier_data
  
  TYPE (wannier_data), allocatable :: wan_in(:,:)
END MODULE wannier_new
