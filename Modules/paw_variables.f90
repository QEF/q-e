MODULE paw_variables
    !
    USE kinds,      ONLY : DP
    USE parameters, ONLY : lqmax, npsx
    USE radial_grids, ONLY: ndmx
    !
    IMPLICIT NONE
    PUBLIC
    SAVE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Legacy (to be removed): !!!!
    INTEGER, PARAMETER :: nbrx = 14  ! max number of beta functions

    !!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Control flags: !!!!

    ! Set to true after initialization, to prevent double allocs:
    LOGICAL              :: is_init = .false.
    ! Analogous to okvan in  "uspp_param" (Modules/uspp.f90)
    LOGICAL :: &
         okpaw              ! if .TRUE. at least one pseudo is PAW

    ! Analogous to tvanp in "uspp_param" (Modules/uspp.f90)
!     LOGICAL :: &
!          tpawp(npsx) = .false.   ! if .TRUE. the atom is of PAW type

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Initialization data: !!!!

    ! the following variables are used to convert spherical harmonics expansion
    ! to radial sampling, they are initialized for an angular momentum up to
    ! l = l_max and (l+1)**2 = lm_max
    ! see function PAW_rad_init for details
    INTEGER              :: l_max  = 0
    INTEGER              :: lm_max = 0
    INTEGER              :: nx     = 0
    REAL(DP),ALLOCATABLE :: ww(:)
    REAL(DP),ALLOCATABLE :: ylm(:,:) ! Y_lm(nx,lm_max)

    ! additional variables for gradient correction
    INTEGER,PARAMETER    :: xlm = 2     ! Additional angular momentum to
                                        ! integrate to have a good GC
    REAL(DP),ALLOCATABLE :: dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:)  ! |d(ylm)/dphi|**2
    REAL(DP),ALLOCATABLE :: cos_th(:),& ! cos(theta) (for divergence)
                            sin_th(:)   ! sin(theta) (for divergence)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Pseudopotential data: !!!!

    ! Analogous to qfunc in "uspp_param" (Modules/uspp.f90)
!     REAL(DP), TARGET :: &
!          pfunc(ndmx,nbrx,nbrx,npsx), &! AE: \phi_{mu}(|r|)-\phi_{nu}(|r|)
!          ptfunc(ndmx,nbrx,nbrx,npsx)  ! PS: \tilde{\phi}_{mu}(|r|)-\tilde{\phi}_{nu}(|r|)

    ! Augmentation on radial grid:
!     TYPE augfun_t
!       REAL(DP), ALLOCATABLE :: fun(:,:,:,:)
!     END TYPE
!     TYPE(augfun_t) :: aug(npsx)
    ! Moments of the augmentation functions
    REAL (DP) :: &
         augmom(nbrx,nbrx,0:6,npsx)  ! moments of PAW augm. functions
    INTEGER :: &
         nraug(npsx)                 ! augm. functions cutoff parameter

    ! Analogous to rho_atc in "atom" (Modules/atom.f90)
!     REAL(DP), TARGET :: &
!          aerho_atc(ndmx,npsx),        &! radial AE core charge density
!          psrho_atc(ndmx,npsx)          ! radial PS core charge density          

    ! Analogous to vloc in "vlocal" (PW/pwcom.f90)
    REAL(DP), ALLOCATABLE, TARGET :: &
         psvloc(:,:)              ! PS local 1-c potential for each atom type

    ! Analogous to dion in "uspp_param" (Modules/uspp.f90)
    REAL(DP) :: &
         kdiff (nbrx,nbrx,npsx)                ! Kinetic energy differences

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! self-consistent variables: !!!!

    ! We need a place to store the radial AE and pseudo potential,
    ! as different atoms may have different max_lm, and different max(|r|)
    ! using a derived type is the way to go
    TYPE paw_saved_potential
        REAL(DP),ALLOCATABLE :: &
            v(:,:,:,:)  ! indexes: |r|, lm, spin, {AE|PS}
    END TYPE
    TYPE(paw_saved_potential),ALLOCATABLE :: &
         saved(:) ! allocated in PAW_rad_init

    ! This type contains some useful data that has to be passed to all
    ! functions, but cannot stay in global variables for parallel:
    TYPE paw_info
        INTEGER :: a ! atom index
        INTEGER :: t ! atom type index
        INTEGER :: m ! atom mesh = g(nt)%mesh
        INTEGER :: w ! number of atomic wavefunctions
        INTEGER :: l ! max angular index l
    END TYPE

    ! Analogous to deeq in "uspp_param" (Modules/uspp.f90)
    REAL(DP), ALLOCATABLE :: &
         dpaw_ae(:,:,:,:),         &! AE D: D^1_{ij}         (except for K.E.)
         dpaw_ps(:,:,:,:)           ! PS D: \tilde{D}^1_{ij} (except for K.E.)

    ! new vectors needed for mixing of augm. channel occupations
    REAL(DP), ALLOCATABLE :: &
         becnew(:,:,:)       ! new augmentation channel occupations


 CONTAINS
    ! From PW/init_paw_1.f90
    SUBROUTINE step_f(f2,f,r,nrs,nrc,pow,mesh)
      USE kinds , ONLY : dp
      !
      ! This routine apply a function which goes smoothly to zero from rs to rc
      ! 
      IMPLICIT NONE
      INTEGER :: mesh
      REAL(DP), INTENT(out):: f2(mesh)
      REAL(DP), INTENT(in) :: f(mesh), r(mesh)
      REAL(DP), INTENT(in) :: pow
      INTEGER :: nrs, nrc 

      INTEGER :: n,i
      REAL(DP) :: rcp, rsp

      rcp = r(nrc)
      rsp = r(nrs)

      DO i=1,mesh
         IF(r(i).LE.rsp) THEN
            f2(i) = f(i)
         ELSE
            IF(r(i).LE.rcp) THEN
               f2(i)=f(i)* (1.d0-3.d0*((r(i)-rsp)/(rcp-rsp))**2+ &
                    2.d0*((r(i)-rsp)/(rcp-rsp))**3)**pow
            ELSE
               f2(i)=0.d0
            ENDIF
         ENDIF

      END DO

    END SUBROUTINE step_f

END MODULE paw_variables
