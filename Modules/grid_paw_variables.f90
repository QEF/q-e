module grid_paw_variables
  !
  !   WARNINGS:
  !
  ! NO spin-orbit
  ! NO EXX
  ! NO rinner > 0
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lqmax, npsx
  USE radial_grids, ONLY: ndmx
  !
  implicit none
  public!              <===
  save

  LOGICAL, PARAMETER :: really_do_paw = .true.

  INTEGER, PARAMETER :: nbrx = 14  ! max number of beta functions

  ! Analogous to okvan in  "uspp_param" (Modules/uspp.f90)
  LOGICAL :: &
       okpaw              ! if .TRUE. at least one pseudo is PAW

  ! Analogous to tvanp in "uspp_param" (Modules/uspp.f90)
  LOGICAL :: &
       tpawp(npsx) = .false.   ! if .TRUE. the atom is of PAW type

  ! Analogous to qfunc in "uspp_param" (Modules/uspp.f90)
  REAL(DP), TARGET :: &
       pfunc(ndmx,nbrx,nbrx,npsx), &! AE: \phi_{mu}(|r|)-\phi_{nu}(|r|)
       ptfunc(ndmx,nbrx,nbrx,npsx)  ! PS: \tilde{\phi}_{mu}(|r|)-\tilde{\phi}_{nu}(|r|)

  ! Augmentation on radial grid:
  TYPE augfun_t
    REAL(DP), ALLOCATABLE :: fun(:,:,:,:)
  END TYPE
  TYPE(augfun_t) :: aug(npsx)
  ! Moments of the augmentation functions
  REAL (DP) :: &
       augmom(nbrx,nbrx,0:6,npsx)     ! moments of PAW augm. functions
  INTEGER :: &
       nraug(npsx)                 ! augm. functions cutoff parameter


#ifdef __GRID_PAW
  REAL(DP), TARGET :: &
       !augfun(ndmx,nbrx,nbrx,0:lqmax,npsx), & ! changed to aug of type augfun_t (see below)
       pmultipole(nbrx,nbrx,0:lqmax,npsx), &! AE multipoles
       ptmultipole(nbrx,nbrx,0:lqmax,npsx)  ! PS multipoles

  ! Analogous to qq in "uspp_param" (Modules/uspp.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       pp(:,:,:),             &! the integrals of p functions in the solid
       ppt(:,:,:)              ! the integrals of pt functions in the solid

  ! Analogous to qrad in "us" (PW/pwcom.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       prad(:,:,:,:),         &! radial FT of P functions
       ptrad(:,:,:,:)          ! radial FT of \tilde{P} functions

  ! Products \Sum_k (P_ij(k)*P_ij'(k))/k**2
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       prodp(:,:,:),              &! AE product in reciprocal space
       prodpt(:,:,:),             &! PS product in reciprocal space
       prod0p(:,:,:),             &! k=0 AE product in reciprocal space
       prod0pt(:,:,:)              ! k=0 PS product in reciprocal space

  ! Analogous to rho in "scf" (PW/pwcom.f90) + index scanning atoms
  REAL(DP), ALLOCATABLE, TARGET :: &
       rho1(:,:,:),             &! 1center AE charge density in real space
       rho1t(:,:,:)              ! 1center PS charge density in real space

  ! Analogous to vr in "scf" (PW/pwcom.f90) + index scanning atoms
  REAL(DP), ALLOCATABLE, TARGET :: &
       vr1(:,:,:),        &! the Hartree+XC potential in real space of rho1
       vr1t(:,:,:)         ! the Hartree+XC potential in real space of rho1t

  ! Analogous to qq in "uspp_param" (Modules/uspp.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       int_r2pfunc(:,:,:),   &! Integrals of r^2 * pfunc(r) (AE)
       int_r2ptfunc(:,:,:)    ! Integrals of r^2 * pfunc(r) (PS)
#endif

  ! Analogous to rho_atc in "atom" (Modules/atom.f90)
  REAL(DP), TARGET :: &
       aerho_atc(ndmx,npsx),        &! radial AE core charge density
       psrho_atc(ndmx,npsx)          ! radial PS core charge density          
  
  ! Analogous to vloc_at in "uspp_param" (Modules/uspp.f90)
  ! actually pseudopotential (AE and PS) on radial grid.
  REAL(DP), TARGET :: &
      aevloc_at(ndmx,npsx),               &! AE descreened potential
      psvloc_at(ndmx,npsx)                 ! PS descreened potential

  ! Analogous to vloc in "vlocal" (PW/pwcom.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
#ifdef __GRID_PAW
       aevloc(:,:),            &! AE local 1-c potential for each atom type
#endif
       psvloc(:,:)              ! PS local 1-c potential for each atom type

#ifdef __GRID_PAW
  ! Analogous to rho_core in "scf" (PW/pwcom.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       aerho_core(:,:),            &! AE core charge density in real space
       psrho_core(:,:)              ! PS core charge density in real space

  !
  REAL(DP), ALLOCATABLE :: &
       radial_distance(:), &   ! radial distance from origin (minimum image conv)
       radial_r(:,:)           ! radial r from origin (minimum image conv)

  ! Analogous to vltot in "scf" (PW/pwcom.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       aevloc_r(:,:),            &! AE local potential in real space
       psvloc_r(:,:)              ! PS local potential in real space

  ! One-center energies
  REAL(DP), ALLOCATABLE, TARGET :: &
       ehart1 (:),                & ! Hartree energy (AE)
       etxc1  (:),                & ! XC: energy (AE)
       vtxc1  (:),                & ! XC: Int V*rho (AE)
       ehart1t(:),                & ! Hartree energy (PS)
       etxc1t (:),                & ! XC: energy (PS)
       vtxc1t (:)                   ! XC: Int V*rho (PS)
#endif

  ! Analogous to dion in "uspp_param" (Modules/uspp.f90)
  REAL(DP) :: &
       kdiff (nbrx,nbrx,npsx)                ! Kinetic energy differences

  ! Analogous to deeq in "uspp_param" (Modules/uspp.f90)
  REAL(DP), ALLOCATABLE :: &
       dpaw_ae(:,:,:,:),         &! AE D: D^1_{ij}         (except for K.E.)
       dpaw_ps(:,:,:,:)           ! PS D: \tilde{D}^1_{ij} (except for K.E.)

#ifdef __GRID_PAW
  ! TMP analogous to rhonew in PW/electrons.f90
  REAL(DP), ALLOCATABLE, TARGET :: &
       rho1new(:,:,:),             &! new 1center AE charge density in real space
       rho1tnew(:,:,:)              ! new 1center PS charge density in real space

  ! analogous to deband and descf in PW/electrons.f90
  REAL(DP) ::  deband_1ae, deband_1ps, descf_1ae, descf_1ps
#endif

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

end module grid_paw_variables
