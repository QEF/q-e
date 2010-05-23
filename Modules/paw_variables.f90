MODULE paw_variables
    !
    USE kinds,      ONLY : DP
    !
    IMPLICIT NONE
    PUBLIC
    SAVE

    !!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Control flags: !!!!

    ! Set to true after initialization, to prevent double allocs:
    LOGICAL              :: paw_is_init = .false.
    ! Analogous to okvan in  "uspp_param" (Modules/uspp.f90)
    LOGICAL :: &
         okpaw = .FALSE.          ! if .TRUE. at least one pseudo is PAW

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Pseudopotential data: !!!!

    ! There is (almost) no pseudopotential data here, it is all stored in the upf type.
    ! See files pseudo_types.f90 and read_uspp.f90

    ! Constant to be added to etot to get all-electron energy 
    REAL(DP) :: total_core_energy = 0._dp
    ! true if all the pseudopotentials are PAW
    LOGICAL  :: only_paw

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Initialization data: !!!!

    INTEGER,PARAMETER    :: lm_fact = 3   ! To converge E_xc integrate up to LM = lm_fact * lm_max
    INTEGER,PARAMETER    :: lm_fact_x = 3 ! As above, for gradient corrected functionals
    INTEGER,PARAMETER    :: xlm = 2       ! Additional factor to add to have a good grad.corr.
    INTEGER,PARAMETER    :: radial_grad_style = 0 ! = 0 or 1, algorithm to use for d/dr

    TYPE paw_radial_integrator
        ! the following variables are used to integrate radial sampling
        INTEGER          :: lmax        ! max l component that can be integrated correctly
        INTEGER          :: ladd        ! additional l max that have been added for grad.corr.
        INTEGER          :: lm_max      ! as above, but +1 and squared
        INTEGER          :: nx          ! number of integration directions
        REAL(DP),POINTER :: ww(:)       ! integration weights (one per direction)
        REAL(DP),POINTER :: ylm(:,:)    ! Y_lm(nx,lm_max)
        REAL(DP),POINTER :: wwylm(:,:)  ! ww(nx) * Y_lm(nx,lm_max)
        ! additional variables for gradient correction
        REAL(DP),POINTER :: dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:)  ! |d(ylm)/dphi|**2
        REAL(DP),POINTER :: cos_phi(:)  ! cos(phi)  
        REAL(DP),POINTER :: sin_phi(:)  ! sin(phi)  
        REAL(DP),POINTER :: cos_th(:)  ! cos(theta)  (for divergence)
        REAL(DP),POINTER :: sin_th(:)  ! sin(theta)  (for divergence)
        REAL(DP),POINTER :: cotg_th(:)  ! cos(theta)/sin(theta)  (for divergence)
    END TYPE
    TYPE(paw_radial_integrator), ALLOCATABLE :: &
        rad(:) ! information to integrate different atomic species

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! self-consistent variables: !!!!

    ! This type contains some useful data that has to be passed to all
    ! functions, but cannot stay in global variables for parallel:
    TYPE paw_info
        INTEGER :: a ! atom index
        INTEGER :: t ! atom type index = itype(a)
        INTEGER :: m ! atom mesh = g(t)%mesh
        INTEGER :: b ! number of beta functions = upf(t)%nbeta
        INTEGER :: l ! max angular index l+1 -> (l+1)**2 is max
                     ! lm index, it is used to allocate rho
        INTEGER :: ae ! tells if we are doing all-electron (1) or pseudo (2)
    END TYPE

    ! Analogous to deeq in "uspp_param" (Modules/uspp.f90)
    REAL(DP), ALLOCATABLE :: &
         ddd_paw(:,:,:)  ! D: D^1_{ij} - \tilde{D}^1_{ij} (only Hxc part)

    REAL(DP), ALLOCATABLE ::  vs_rad(:,:,:)

 END MODULE paw_variables
