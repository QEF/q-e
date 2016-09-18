MODULE tsvdw_module
!
!----------------------------------------------------------------------------------------------------------------
! TS-vdW Code Version 14.0 (RAD/BS, Princeton University, February 2013)
!----------------------------------------------------------------------------------------------------------------
! All quantities necessary for the evaluation of the TS-vdW energy and forces are computed on the real-space
! mesh using linear interpolation of the atomic pseudo-densities and their first derivatives which have been
! mapped onto linear equispaced atomic grids from their original form computed on radial atomic grids via the
! ATOMIC code.
!----------------------------------------------------------------------------------------------------------------
! SYNOPSIS: radial form of rhoA & drhoA mapped onto linear grid;
!           atrho & rhosad on real-space mesh via linear interpolation;
!           integration on spherical atomic domains (subsets of real-space mesh);
!           quadratic veff derivatives computed linearly using sparse domain intersection algorithm.
!----------------------------------------------------------------------------------------------------------------
!
USE cell_base,          ONLY: h                  !h matrix for converting between r and s coordinates via r = h s
USE cell_base,          ONLY: ainv               !h^-1 matrix for converting between r and s coordinates via s = h^-1 r)
USE cell_base,          ONLY: omega              !cell volume (in au^3)
USE constants,          ONLY: pi                 !pi in double-precision
USE fft_base,           ONLY: dfftp              !FFT derived data type 
USE funct,              ONLY: get_iexch          !retrieves type of exchange utilized in functional
USE funct,              ONLY: get_icorr          !retrieves type of correlation utilized in functional
USE funct,              ONLY: get_igcx           !retrieves type of gradient correction to exchange utilized in functional
USE funct,              ONLY: get_igcc           !retrieves type of gradient correction to correlation utilized in functional
USE io_global,          ONLY: stdout             !print/write argument for standard output (to output file)
USE ions_base,          ONLY: nat                !number of total atoms (all atomic species)
USE ions_base,          ONLY: nsp                !number of unique atomic species
USE ions_base,          ONLY: na                 !number of atoms within each atomic species
USE ions_base,          ONLY: ityp               !ityp(i):=type/species of ith atom
USE ions_base,          ONLY: atm                !atm(j):=name of jth atomic species (3 characters)
USE kinds,              ONLY: DP                 !double-precision kind (selected_real_kind(14,200))
! the charge density is parallelized over the "band group" or processors
USE mp_bands,           ONLY: nproc_bgrp         !number of processors
USE mp_bands,           ONLY: me_bgrp            !processor number (0,1,...,nproc_bgrp-1)
USE mp_bands,           ONLY: intra_bgrp_comm    !standard MPI communicator
! atoms are parallelized over the "image group"
USE mp_images,          ONLY: nproc_image        !number of processors
USE mp_images,          ONLY: me_image           !processor number (0,1,...,nproc_image-1)
USE mp_images,          ONLY: intra_image_comm   !standard MPI communicator
USE mp,                 ONLY: mp_sum             !MPI collection with sum
USE parallel_include                             !MPI header
USE uspp_param,         ONLY: upf                !atomic pseudo-potential data
!
IMPLICIT NONE
!
SAVE
!
! PUBLIC variables 
!
LOGICAL, PUBLIC :: vdw_isolated    ! isolated system control
REAL(DP), PUBLIC:: vdw_econv_thr   ! energy convergence threshold for periodic systems
REAL(DP), PUBLIC :: EtsvdW                                   !the TS-vdW energy
REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: UtsvdW        !the TS-vdW wavefunction forces (dispersion potential)
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: FtsvdW      !the TS-vdW ionic forces (-dE/dR)
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: HtsvdW      !the TS-vdW cell forces (dE/dh)
REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: VefftsvdW     !the TS-vdW effective Hirshfeld volume
!
! PRIVATE variables 
!
INTEGER, PARAMETER, PRIVATE :: NgpA=1000                     !number of grid points for linear equispaced atomic grids (current value=1000pts)
INTEGER, PARAMETER, PRIVATE :: bsint=BIT_SIZE(NgpA)          !integer bit size (for use in bit array manipulation)
INTEGER, PRIVATE :: me                                       !processor number (1,2,...,nproc_image)
INTEGER, PRIVATE :: iproc                                    !processor dummy index
INTEGER, PRIVATE :: nr1,nr2,nr3                              !real space grid dimensions (global first, second, and third dimensions of the 3D grid)
INTEGER, PRIVATE :: nr1r,nr2r,nr3r                           !reduced real space grid dimensions (global first, second, and third dimensions of the 3D grid)
REAL(DP), PRIVATE :: ddamp                                   !damping function parameter #1
REAL(DP), PRIVATE :: sR                                      !damping function parameter #2
REAL(DP), PRIVATE :: spcutAmax                               !maximum radial cutoff for all atomic species
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: nstates       !number of atoms per processor
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: sdispls       !send displacement (offset) array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: rdispls       !receive displacement (offset) array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: sendcount     !send count array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: recvcount     !receive count array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: istatus       !MPI status array
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: NsomegaA      !number of points in the spherical atomic integration domain
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: NsomegaAr     !number of points in the reduced spherical atomic integration domain
INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: npair         !number of unique atom pairs
INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: pair        !unique atom pair overlap matrix
INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: gomegar     !precursor to spherical atomic integration domain (intersection bit array)
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: somegaA   !spherical atomic integration domain
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: somegaAr  !reduced spherical atomic integration domain
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: gomegaAr  !reduced spherical atomic integration domain (intersection bit array)
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE:: predveffAdn   !atomic dispersion potential prefactor
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: vfree        !free atomic volumes for each atomic species
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: dpfree       !free atomic static dipole polarizability for each atomic species
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: R0free       !free atomic vdW radius for each atomic species
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: C6AAfree     !free atomic homonuclear C6 coefficient for each atomic species
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: veff         !effective atomic volumes for each atom in the simulation cell
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: dpeff        !effective atomic static dipole polarizability for each atom in the simulation cell
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: R0eff        !effective atomic vdW radius for each atom in the simulation cell
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: C6AAeff      !effective atomic homonuclear C6 coefficient for each atom in the simulation cell
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: rhosad       !molecular pro-density (superposition of atomic densities) on real-space mesh
REAL(DP), DIMENSION(:), ALLOCATABLE, PRIVATE :: rhotot       !molecular charge density on real-space mesh
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE:: dveffAdn    !the local copy of the TS-vdW wavefunction forces (dispersion potential)
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: spgrd      !linear equispaced grid for each atomic species
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: sprho      !atomic pseudo-density for each atomic species
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: spdrho     !first derivative of atomic pseudo-density for each atomic species
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: spdata     !linear grid cutoff (is,1) and linear grid spacing (is,2) for each atomic species
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: LIA        !A coefficient for linear interpolation of rhoA
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: LIB        !B coefficient for linear interpolation of rhoA
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: dLIA       !A coefficient for linear interpolation of drhoA
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: dLIB       !B coefficient for linear interpolation of drhoA
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: atxyz      !Cartesian coordinates of ions adjusted according to PBC
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: C6ABfree   !free atomic heteronuclear C6 coefficient for each atom pair
REAL(DP), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: C6ABeff    !effective atomic heteronuclear C6 coefficient for each atom pair
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: dveffdR  !first derivative of effective volume wrt nuclear displacement
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: dveffdh  !first derivative of effective volume wrt cell displacement
!
! PUBLIC subroutines
!
PUBLIC :: tsvdw_initialize
PUBLIC :: tsvdw_calculate
PUBLIC :: tsvdw_finalize
!
! PRIVATE subroutines
!
PRIVATE :: tsvdw_para_init
PRIVATE :: tsvdw_pbc
PRIVATE :: tsvdw_unique_pair
PRIVATE :: tsvdw_rhotot
PRIVATE :: tsvdw_screen
PRIVATE :: tsvdw_veff
PRIVATE :: tsvdw_dveff
PRIVATE :: tsvdw_effqnts
PRIVATE :: tsvdw_energy
PRIVATE :: tsvdw_wfforce
PRIVATE :: tsvdw_cleanup
PRIVATE :: Num1stDer
PRIVATE :: CubSplCoeff
PRIVATE :: GetVdWParam
  !
  !
  CONTAINS
  !
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_initialize()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local Variables 
  !
  LOGICAL :: uniform_grid=.FALSE.
  INTEGER :: ip,iq,ir,is,it,NrgpA,NrgpintA,icutrA,Ndim
  REAL(DP) :: dxA,gfctrA,vref,eref,verr,d,dk1,dk2,dk3,num,den,drab,f1,f2,f3,L1,L2,L3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: atgrdr,atgrdrab,atrhor,datrhor,d2atrhor,CSA,CSB,CSC,CSD
  !
  ! Start of calculation banner...
  !
  WRITE(stdout,*)
  WRITE(stdout,'(3X,"TS-vdW initialization")')
  WRITE(stdout,'(3X,"---------------------")')
  WRITE(stdout,*)
  !
  ! Error messages for inconsistencies with current version of code...
  !
  !RAD: Have we missed any inconsistencies?
  !
  ! Setup variables for use in TS-vdW module...
  !
  nr1=dfftp%nr1; nr2=dfftp%nr2; nr3=dfftp%nr3 
  nr1r=nr1/2; nr2r=nr2/2; nr3r=nr3/2 
  IF(MOD(nr1,2).EQ.1) nr1r=(nr1+1)/2
  IF(MOD(nr2,2).EQ.1) nr2r=(nr2+1)/2
  IF(MOD(nr3,2).EQ.1) nr3r=(nr3+1)/2
  !
  ! Initialize the TS-vdW ionic forces, cell forces, and dispersion potential (wavefunction forces)...
  !
  ALLOCATE(FtsvdW(3,nat)); FtsvdW=0.0_DP
  ALLOCATE(HtsvdW(3,3)); HtsvdW=0.0_DP
  !
  ! Initialization of TS-vdW Hirshfeld effective volume public variable ... used in CP print_out.f90
  !
  ALLOCATE(VefftsvdW(nat)); VefftsvdW=0.0_DP
  !
  ALLOCATE(UtsvdW(dfftp%nnr)); UtsvdW=0.0_DP
  !
  ! Set ddamp damping function parameter (set to 20 and functional independent)...
  !
  WRITE(stdout,'(3X,"Determining TS-vdW damping function parameters...")')
  ddamp=20.0_DP
  WRITE(stdout,'(5X,"ddamp = ",F9.6)') ddamp
  !
  ! Set sR damping function parameter (functional dependent and currently only available for PBE & PBE0)...
  !
  IF (get_iexch().EQ.1.AND.get_icorr().EQ.4.AND.get_igcx().EQ.3.AND.get_igcc().EQ.4) THEN
    !
    sR=0.94_DP !PBE=sla+pw+pbx+pbc
    !
  ELSE IF (get_iexch().EQ.6.AND.get_icorr().EQ.4.AND.get_igcx().EQ.8.AND.get_igcc().EQ.4) THEN
    !
    sR=0.96_DP !PBE0=pb0x+pw+pb0x+pbc !RAD/BS: This line will not work in CP unless PBE0 code update funct.f90...
    !
  ELSE 
    !
    CALL errore('tsvdw','TS-vdW sR parameter only available for PBE and PBE0 functionals...',1)
    !
  END IF
  !
  WRITE(stdout,'(5X,"sR = ",F9.6)') sR
  !
  ! Allocate and initialize species-specific quantities...
  !
  ALLOCATE(vfree(nsp)); vfree=0.0_DP
  ALLOCATE(dpfree(nsp)); dpfree=0.0_DP
  ALLOCATE(R0free(nsp)); R0free=0.0_DP
  ALLOCATE(C6AAfree(nsp)); C6AAfree=0.0_DP
  ALLOCATE(C6ABfree(nsp,nsp)); C6ABfree=0.0_DP
  ALLOCATE(spdata(nsp,2)); spdata=0.0_DP
  ALLOCATE(spgrd(nsp,0:NgpA)); spgrd=0.0_DP
  ALLOCATE(sprho(nsp,0:NgpA)); sprho=0.0_DP
  ALLOCATE(spdrho(nsp,0:NgpA)); spdrho=0.0_DP
  ALLOCATE(LIA(nsp,0:NgpA)); LIA=0.0_DP
  ALLOCATE(LIB(nsp,0:NgpA)); LIB=0.0_DP
  ALLOCATE(dLIA(nsp,0:NgpA)); dLIA=0.0_DP
  ALLOCATE(dLIB(nsp,0:NgpA)); dLIB=0.0_DP
  !
  spcutAmax=0.0_DP
  !
  ! Loop over atomic species and extract species-dependent quantities to modular arrays...
  !
  DO is=1,nsp
    !
    ! Obtain the radial grid and radial atomic pseudo-density from pseudo-potential file (via upf module) for
    ! the given atomic species.  Convert the radial atomic pseudo-density to the real atomic pseudo-density using
    ! rho_real(r) = rho_radial(r) / (4*pi*r^2)...
    !
    WRITE(stdout,'(3X,"Initializing species # ",I3," with atomic symbol ",A3)') is,atm(is)
    !
    ! Read in the number of grid points in radial mesh from upf...
    !
    NrgpA=upf(is)%mesh
    !
    ! Transfer radial atomic grid (in upf) to local atgrdr array...
    !
    ALLOCATE(atgrdr(NrgpA)); atgrdr=0.0_DP
    !
    DO ir=1,NrgpA
      !
      atgrdr(ir)=upf(is)%r(ir)
      !
    END DO
    !
    ! Transfer radial atomic grid spacing (in upf) to local atgrdrab array...
    !
    ALLOCATE(atgrdrab(NrgpA)); atgrdrab=0.0_DP
    !
    DO ir=1,NrgpA
      !
      atgrdrab(ir)=upf(is)%rab(ir)
      !
    END DO
    !
    ! Determine whether radial grid is logarithmic/exponential or equispaced/uniform...
    !
    drab=atgrdrab(NrgpA)-atgrdrab(1)
    IF (DABS(drab).LT.(1.0E-6_DP)) uniform_grid=.TRUE.
    IF (uniform_grid) WRITE(stdout,'(5X,"Equispaced/Uniform radial atomic grid detected...")')
    !
    ! ----------------------------------------------------------------
    ! Logarithmic/Exponential grid (3 parameters: zmesh, xmin, dxA)
    ! ----------------------------------------------------------------
    !
    !   For i = 1,2,...,NrgpA: 
    !     r(i) = exp[xmin+(i-1)*dxA]/zmesh 
    !          = exp[xmin]/zmesh * exp[(i-1)*dxA] 
    !          = gfctrA * exp[(i-1)*dxA]
    !   rab(i) = r(i) * dxA 
    !
    !   Assumptions: grid does NOT start from zero (use simpson_cp90()).
    !
    ! ---------------------------------------------
    ! Equispaced/Uniform grid (1 parameter: dxA)
    ! ---------------------------------------------
    !
    !   For i = 1,2,...,NrgpA: 
    !     r(i) = (i-1) * dxA
    !   rab(i) = dxA
    !
    !   Assumptions: grid starts from zero (use simpson() for integration).
    !
    ! Determine atomic radial grid parameters...
    !
    IF (uniform_grid.EQV..TRUE.) THEN
      !
      gfctrA=1.0_DP
      dxA=atgrdrab(1)
      !
    ELSE
      !
      gfctrA=upf(is)%r(1)
      dxA=DLOG(upf(is)%r(2)/upf(is)%r(1))
      !
    END IF
    !
    WRITE(stdout,'(5X,"Radial grid parameter: NrgpA is ",I5,".")') NrgpA
    WRITE(stdout,'(5X,"Radial grid parameter: gfctrA is ",F9.6,".")') gfctrA
    WRITE(stdout,'(5X,"Radial grid parameter: dxA is ",F9.6,".")') dxA
    !
    ! Transfer radial atomic pseudo-density to atrhor array...
    ! Convert radial atomic pseudo-density to real atomic pseudo-density [n(r) = nrad(r)/(4*pi*r^2)]...
    !
    ALLOCATE(atrhor(NrgpA)); atrhor=0.0_DP
    !
    IF (uniform_grid.EQV..TRUE.) THEN
      !
      DO ir=2,NrgpA
        !
        atrhor(ir)=(upf(is)%rho_at(ir))/(4.0_DP*pi*atgrdr(ir)**(2.0_DP)) ! skip point at r=0...
        !
      END DO
      !
      ! Quadratic extrapolation of the atomic density to r=0...
      !
      L1=((0.0_DP-atgrdr(3))*(0.0_DP-atgrdr(4)))/((atgrdr(2)-atgrdr(3))*(atgrdr(2)-atgrdr(4)))
      L2=((0.0_DP-atgrdr(2))*(0.0_DP-atgrdr(4)))/((atgrdr(3)-atgrdr(2))*(atgrdr(3)-atgrdr(4)))
      L3=((0.0_DP-atgrdr(2))*(0.0_DP-atgrdr(3)))/((atgrdr(4)-atgrdr(2))*(atgrdr(4)-atgrdr(3)))
      atrhor(1)=L1*atrhor(2)+L2*atrhor(3)+L3*atrhor(4)
      !
    ELSE
      !
      DO ir=1,NrgpA
        !
        atrhor(ir)=(upf(is)%rho_at(ir))/(4.0_DP*pi*atgrdr(ir)**(2.0_DP))
        !
      END DO
      !
    END IF
    !
    ! Set NrgpintA as the number of grid points (which must be odd) used during numerical integration using Simpson's rule...
    !
    IF (IAND(NrgpA,1).EQ.1) THEN
      !
      NrgpintA=NrgpA
      !
    ELSE
      !
      NrgpintA=NrgpA-1
      !
    END IF
    !
    ! Compute the number of electrons (eref) for each atomic species via numerical integration
    ! of the atomic pseudo-density on the radial atomic grid using Simpson's rule...
    !
    eref=0.0_DP
    !
    DO ir=1,NrgpintA-2,2
      !  
      f1=atrhor(ir  )*atgrdrab(ir  )*atgrdr(ir  )**(2.0_DP)  ! integrated quantity is rho
      f2=atrhor(ir+1)*atgrdrab(ir+1)*atgrdr(ir+1)**(2.0_DP)
      f3=atrhor(ir+2)*atgrdrab(ir+2)*atgrdr(ir+2)**(2.0_DP)
      !
      eref=eref+(f1+4.0_DP*f2+f3) 
      !
    END DO
    !
    eref=(4.0_DP*pi/3.0_DP)*eref
    WRITE(stdout,'(5X,"The number of valence electrons, eref, is ",F25.15,".")') eref
    !
    ! Compute the reference free atom volume (vref) for each atomic species via numerical integration
    ! of the atomic pseudo-density on the radial atomic grid using Simpson's rule...
    !
    vref=0.0_DP
    !
    DO ir=1,NrgpintA-2,2
      !  
      f1=atrhor(ir  )*atgrdrab(ir  )*atgrdr(ir  )**(5.0_DP)  ! integrated quantity is rho * r^3
      f2=atrhor(ir+1)*atgrdrab(ir+1)*atgrdr(ir+1)**(5.0_DP)
      f3=atrhor(ir+2)*atgrdrab(ir+2)*atgrdr(ir+2)**(5.0_DP)
      !
      vref=vref+(f1+4.0_DP*f2+f3) 
      !
    END DO
    !
    vref=(4.0_DP*pi/3.0_DP)*vref
    WRITE(stdout,'(5X,"The reference free atom volume, vref, is ",F25.15," bohr^3.")') vref
    !
    ! Using the reference free atom volume, determine an acceptable radial grid cutoff value such that the 
    ! free atom volume obtained using this cutoff does not deviate from the reference value by more than 1.0%.
    !
    WRITE(stdout,'(5X,"Determining intial radial grid cutoff...")')
    !
    DO iq=5,NrgpintA,2
      !
      vfree(is)=0.0_DP
      verr=0.0_DP
      !
      DO ir=1,iq-2,2
        !
        f1=atrhor(ir  )*atgrdrab(ir  )*atgrdr(ir  )**(5.0_DP)  ! integrated quantity is rho * r^3
        f2=atrhor(ir+1)*atgrdrab(ir+1)*atgrdr(ir+1)**(5.0_DP)
        f3=atrhor(ir+2)*atgrdrab(ir+2)*atgrdr(ir+2)**(5.0_DP)
        !
        vfree(is)=vfree(is)+(f1+4.0_DP*f2+f3)
        !
      END DO
      !
      vfree(is)=(4.0_DP*pi/3.0_DP)*vfree(is)
      verr=(vref-vfree(is))/vref*100.0_DP
      !
      IF (verr.LE.1.0_DP) THEN
        !
        icutrA=iq
        !
        WRITE(stdout,'(5X,"An acceptable radial grid cutoff was determined by retaining ",I4," of ",I4," radial grid points.")') &
              icutrA,NrgpA
        !
        EXIT
        !
      END IF
      !
    END DO
    !
    WRITE(stdout,'(5X,"The magnitude of the atomic pseudo-density at the radial grid cutoff is ",ES13.6,".")') atrhor(icutrA)
    WRITE(stdout,'(5X,"Using this radial grid cutoff value of ",F25.15," au:")') atgrdr(icutrA)
    WRITE(stdout,'(5X,"The free atom volume computed with this cutoff is ",F25.15," bohr^3 with an error of ",F6.3,"%.")') &
         vfree(is),verr
    !
    ! Form 1st derivative of atrhor for input into cubic spline coefficient subroutine...
    !
    ALLOCATE(datrhor(NrgpA)); datrhor=0.0_DP
    CALL Num1stDer(atgrdr,atrhor,NrgpA,dxA,datrhor)
    !
    ! For logarithmic/exponential grid, transform linear derivative back to radial grid...
    !
    IF (.NOT.uniform_grid) THEN
      !
      DO ir=1,NrgpA
        !  
        datrhor(ir)=datrhor(ir)/atgrdr(ir)
        !
      END DO
      !
    END IF
    !
    ! Form the coefficients of the cubic spline interpolant (2nd derivatives) for the real atomic pseudo-density
    ! for use during cubic spline interpolation of the pseudo-density onto the linear equispaced atomic grid...
    !
    ALLOCATE(d2atrhor(NrgpA)); d2atrhor=0.0_DP
    CALL CubSplCoeff(atgrdr,atrhor,NrgpA,datrhor,d2atrhor)
    !
    ! Precompute cubic spline interpolation vectors (utilizing Taylor series form) via:
    !
    !            y(x) = CSA + CSB*(x-x(k)) + CSC*(x-x(k))**2 + CSD*(x-x(k))**3
    !
    ALLOCATE(CSA(NrgpA)); CSA=0.0_DP
    ALLOCATE(CSB(NrgpA)); CSB=0.0_DP
    ALLOCATE(CSC(NrgpA)); CSC=0.0_DP
    ALLOCATE(CSD(NrgpA)); CSD=0.0_DP
    !
    DO ir=1,NrgpA-1
      !
      ! CSA(k) := y(k)
      !
      CSA(ir)=atrhor(ir)
      !
      ! CSB(k) := delta(y)/delta(x) - 1/3*delta(x)*y''(k) - 1/6*delta(x)*y''(k+1)
      !
      CSB(ir)=(atrhor(ir+1)-atrhor(ir))/(atgrdr(ir+1)-atgrdr(ir))
      CSB(ir)=CSB(ir)-((1.0_DP/3.0_DP)*(atgrdr(ir+1)-atgrdr(ir))*d2atrhor(ir))
      CSB(ir)=CSB(ir)-((1.0_DP/6.0_DP)*(atgrdr(ir+1)-atgrdr(ir))*d2atrhor(ir+1))
      !
      ! CSC(k) := 1/2*y''(k)
      !
      CSC(ir)=(1.0_DP/2.0_DP)*d2atrhor(ir)
      !
      ! CSD(k) := 1/6*delta(y'')/delta(x)
      !
      CSD(ir)=((1.0_DP/6.0_DP)*(d2atrhor(ir+1)-d2atrhor(ir))/(atgrdr(ir+1)-atgrdr(ir)))
      !
    END DO
    !
    ! Pack species-specific radial cutoff into (is,1) of spdata array...
    !
    spdata(is,1)=atgrdr(icutrA)
    IF (spdata(is,1).GT.spcutAmax) spcutAmax=spdata(is,1) 
    !
    ! Compute and pack grid spacing of species-specific linear equispaced grid into (is,2) of spdata array...
    !
    spdata(is,2)=(atgrdr(icutrA)+1.0_DP)/DBLE(NgpA) !include additional buffer of 1 bohr...
    WRITE(stdout,'(5X,"Linear grid spacing was computed as: ",F25.15," bohr.")') spdata(is,2)
    !
    ! Form linear equispaced atomic grid (NOT including point at r=0) and pack into argument (is,:) of spgrd array...
    !
    DO ip=1,NgpA
      !
      spgrd(is,ip)=DBLE(ip)*spdata(is,2)   
      !
    END DO
    !
    ! Map atomic pseudo-density (currently on the radial atomic grid) onto linear equispaced atomic grid using
    ! cubic spline interpolation...Form first derivative of the atomic pseudo-density on the linear equispaced
    ! atomic grid via differentiation of the cubic spline interpolant...
    !
    DO ip=1,NgpA
      !  
      d=spgrd(is,ip)
      !
      IF (uniform_grid.EQV..TRUE.) THEN
        !
        ir=INT(d/dxA)+1 !since the equispaced/uniform grid first point is at r=0...
        !
      ELSE
        !
        ir=FLOOR(DLOG(d*EXP(dxA)/gfctrA)/dxA)
        !
      END IF
      !
      dk1=d-atgrdr(ir); dk2=dk1*dk1; dk3=dk2*dk1
      sprho(is,ip)=CSA(ir)+CSB(ir)*dk1+CSC(ir)*dk2+CSD(ir)*dk3      !Pack density into argument (is,:) of sprho array
      spdrho(is,ip)=CSB(ir)+2.0_DP*CSC(ir)*dk1+3.0_DP*CSD(ir)*dk2   !Pack density derivative into argument (is,:) of spdrho array
      !  
    END DO
    !
    ! For computational efficiency during the remainder of the calculation, extrapolate sprho and spdrho to
    ! include the point at r=0 (this eliminates an if statement in crucial inner loops)...
    ! Use quadratic extrapolation to obtain these points...Hence the 0:NgpA dimension above...
    !
    spgrd(is,0)=0.0_DP   !Extend linear grid to include point at r=0...
    !
    L1=((0.0_DP-spgrd(is,2))*(0.0_DP-spgrd(is,3)))/((spgrd(is,1)-spgrd(is,2))*(spgrd(is,1)-spgrd(is,3)))
    L2=((0.0_DP-spgrd(is,1))*(0.0_DP-spgrd(is,3)))/((spgrd(is,2)-spgrd(is,1))*(spgrd(is,2)-spgrd(is,3)))
    L3=((0.0_DP-spgrd(is,1))*(0.0_DP-spgrd(is,2)))/((spgrd(is,3)-spgrd(is,1))*(spgrd(is,3)-spgrd(is,2)))
    sprho(is,0)=L1*sprho(is,1)+L2*sprho(is,2)+L3*sprho(is,3) !Extend atomic pseudo-density to include point at r=0...
    spdrho(is,0)=L1*spdrho(is,1)+L2*spdrho(is,2)+L3*spdrho(is,3) !Extend atomic pseudo-density derivative to include point at r=0...
    !
    ! Throughout the remainder of the code, to map the atomic quantities onto the real-space mesh, we will be
    ! utilizing the Taylor series form of linear interpolation, given by:
    !
    !   y(x) = LIA + LIB*(x-x(k))    y'(x) = dLIA + dLIB*(x-x(k)) 
    !
    ! for x(k) <= x <= x(k+1)...
    !
    DO ip=0,NgpA-1
      !
      ! LIA(k) := y(k)
      !
      LIA(is,ip)=sprho(is,ip)
      dLIA(is,ip)=spdrho(is,ip)
      !
      ! LIB(k) := delta(y)/delta(x)
      !
      LIB(is,ip)=(sprho(is,ip+1)-sprho(is,ip))/(spgrd(is,ip+1)-spgrd(is,ip))
      dLIB(is,ip)=(spdrho(is,ip+1)-spdrho(is,ip))/(spgrd(is,ip+1)-spgrd(is,ip))
      !
    END DO
    !
    ! Populate reference free atom quantities...
    !
    CALL GetVdWParam(atm(is),C6AAfree(is),dpfree(is),R0free(is))
    !
    WRITE(stdout,'(5X,"The free atom static dipole polarizability is ",F13.6," bohr^3.")') dpfree(is)
    WRITE(stdout,'(5X,"The free atom homonuclear C6 coefficient is ",F13.6," Hartree bohr^6.")') C6AAfree(is)
    WRITE(stdout,'(5X,"The free atom vdW radius is ",F13.6," bohr.")') R0free(is)
    !
    ! Clean-up all species-specific temporary arrays
    !
    IF (ALLOCATED(atgrdr))   DEALLOCATE(atgrdr)
    IF (ALLOCATED(atgrdrab)) DEALLOCATE(atgrdrab)
    IF (ALLOCATED(atrhor))   DEALLOCATE(atrhor)
    IF (ALLOCATED(datrhor))  DEALLOCATE(datrhor)
    IF (ALLOCATED(d2atrhor)) DEALLOCATE(d2atrhor)
    IF (ALLOCATED(CSA))      DEALLOCATE(CSA)
    IF (ALLOCATED(CSB))      DEALLOCATE(CSB)
    IF (ALLOCATED(CSC))      DEALLOCATE(CSC)
    IF (ALLOCATED(CSD))      DEALLOCATE(CSD)
    !
  END DO !is
  !
  ! Compute free heteronuclear C6 coefficient matrix...
  !  C6ABfree(A,B)=[2*C6AAfree(A)*C6AAfree(B)]/[(dpfree(B)/dpfree(A))*C6AAfree(A)+(dpfree(A)/dpfree(B))*C6AAfree(B)]
  !
  DO is=1,nsp
    !
    DO it=1,nsp
      !
      num=2.0_DP*C6AAfree(is)*C6AAfree(it)
      den=(dpfree(it)/dpfree(is))*C6AAfree(is)+(dpfree(is)/dpfree(it))*C6AAfree(it)
      C6ABfree(is,it)=num/den
      !
    END DO
    !
  END DO
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_initialize
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_calculate(tauin, rhor)
  !--------------------------------------------------------------------------------------------------------------
  ! TS-vdW Management Code: Manages entire calculation of TS-vdW energy, wavefunction forces, and ion forces via
  ! calls to PRIVATE subroutines below (called in each MD step). The calls to tsvdw_initialize and tsvdw_finalize
  ! are done once at the beginning (init_run) and the end (terminate_run). 
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP), INTENT(IN) :: rhor(:,:)
  REAL(DP) :: tauin(3,nat)
  !
  ! Parallel initialization...
  !
  CALL tsvdw_para_init() 
  !
  ! Move all atoms into simulation cell by adjusting Cartesian coordinates according to PBCs...
  !
  CALL tsvdw_pbc(tauin)
  !
  ! Compute unique atom pair list... 
  !
  CALL tsvdw_unique_pair()
  !
  ! Obtain molecular charge density given on the real-space mesh...
  !
  CALL tsvdw_rhotot( rhor )
  !
  ! Determine spherical atomic integration domains and atom overlap (bit array)...
  ! Compute molecular pro-density (superposition of atomic densities) on the real-space mesh...
  ! Compute functional derivative of vdW energy wrt charge density (numerator only)...
  !
  CALL tsvdw_screen()
  !
  ! Compute effective volume for each atom in the simulation cell...
  ! Complete functional derivative of vdW energy wrt charge density... 
  !
  CALL tsvdw_veff()
  !
  ! Calculate first derivative of veff wrt nuclear and cell displacements...
  !
  CALL tsvdw_dveff()
  !
  ! Calculate effective quantities for each atom in the simulation cell...
  !
  CALL tsvdw_effqnts()
  !
  ! Calculate total TS-vdW energy, dispersion potential prefactor, ionic forces, and cell forces...
  !
  CALL tsvdw_energy()
  !
  ! Calculate total TS-vdW wavefunction forces (dispersion potential)...
  !
  CALL tsvdw_wfforce()
  !
  ! Deallocate all arrays specific to tsvdw_calculate...
  !
  CALL tsvdw_cleanup()
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_calculate
  !--------------------------------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------------------------------- 
  SUBROUTINE tsvdw_para_init()
  !-------------------------------------------------------------------------------------------------------------- 
  !
  IMPLICIT NONE
  !
  INTEGER :: i,j,k
  !
  me=me_image+1
  !
  ALLOCATE(nstates(nproc_image)); nstates=0
  ALLOCATE(sdispls(nproc_image)); sdispls=0
  ALLOCATE(sendcount(nproc_image)); sendcount=0
  ALLOCATE(rdispls(nproc_image)); rdispls=0
  ALLOCATE(recvcount(nproc_image)); recvcount=0
  ALLOCATE(istatus(nproc_image)); istatus=0
  !
  ! Assign workload of atoms over nproc_image processors
  !
  IF (nat.LE.nproc_image) THEN
    !
    DO i=1,nat
      !
      nstates(i)=1
      !
    END DO
    !
  ELSE
    !
    k=0
    !
10  DO j=1,nproc_image
      !
      nstates(j)=nstates(j)+1
      !
      k=k+1
      !
      IF (k.GE.nat) GO TO 20
      !
    END DO
    !
    IF (k.LT.nat) GO TO 10
    !
  END IF
  !
  20 CONTINUE
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_para_init
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_pbc(tauin)
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP) :: tauin(3,nat)
  !
  ! Local variables
  !
  INTEGER :: ia
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: atxyzs
  !
  ! Initialization of PBC-adjusted Cartesian coordinates...
  !
  ALLOCATE(atxyz(3,nat)); atxyz=0.0_DP
  ALLOCATE(atxyzs(3,nat)); atxyzs=0.0_DP
  !
  ! Adjust Cartesian coordinates of ions according to periodic boundary conditions...
  ! N.B.: PBC are imposed here in the range [0,1)... 
  !
  DO ia = 1, nat
    !
    atxyzs(1,ia)=ainv(1,1)*tauin(1,ia)+ainv(1,2)*tauin(2,ia)+ainv(1,3)*tauin(3,ia)   ! s = h^-1 r
    atxyzs(2,ia)=ainv(2,1)*tauin(1,ia)+ainv(2,2)*tauin(2,ia)+ainv(2,3)*tauin(3,ia)   ! s = h^-1 r
    atxyzs(3,ia)=ainv(3,1)*tauin(1,ia)+ainv(3,2)*tauin(2,ia)+ainv(3,3)*tauin(3,ia)   ! s = h^-1 r
    !
    atxyzs(1,ia)=atxyzs(1,ia)-FLOOR(atxyzs(1,ia))   ! impose PBC on s in range: [0,1)
    atxyzs(2,ia)=atxyzs(2,ia)-FLOOR(atxyzs(2,ia))   ! impose PBC on s in range: [0,1)
    atxyzs(3,ia)=atxyzs(3,ia)-FLOOR(atxyzs(3,ia))   ! impose PBC on s in range: [0,1)
    !
    atxyz(1,ia)=h(1,1)*atxyzs(1,ia)+h(1,2)*atxyzs(2,ia)+h(1,3)*atxyzs(3,ia)   ! r = h s
    atxyz(2,ia)=h(2,1)*atxyzs(1,ia)+h(2,2)*atxyzs(2,ia)+h(2,3)*atxyzs(3,ia)   ! r = h s
    atxyz(3,ia)=h(3,1)*atxyzs(1,ia)+h(3,2)*atxyzs(2,ia)+h(3,3)*atxyzs(3,ia)   ! r = h s
    !
  END DO
  !
  IF (ALLOCATED(atxyzs))   DEALLOCATE(atxyzs)
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_pbc
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_unique_pair()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,ib,ias,ibs,ip,ir,i,j,k,jj,nj_max,nbmax,num,num1,jj_neib_of_i
  REAL(DP) :: spcutA,spcutB,dAB(3),dAB2(3)
  INTEGER, DIMENSION(:), ALLOCATABLE :: nj,overlap2
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: overlap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dABmic
  !
  CALL start_clock('tsvdw_pair')
  !
  ! Allocate and initialize temporary arrays...
  !
  ALLOCATE(dABmic(nat)); dABmic=0.0_DP
  ALLOCATE(overlap(nat,nat)); overlap=0
  ALLOCATE(overlap2(nat)); overlap2=0
  ALLOCATE(nj(nat)); nj=0
  !
  ! Outer loop over atoms A to form non-unique atom pair overlap matrix...
  !
  DO ia=1,nat
    !
    nj(ia)=0; dABmic=0.0_DP
    !
    ! Connect atom type with species-dependent quantities...
    !
    ias=ityp(ia)
    !
    ! Transfer species-specific cutoff to spcutA...
    !
    spcutA=spdata(ias,1)
    !
    ! Inner loop over atoms B...
    !
    DO ib=1,nat 
      !
      IF(ib.NE.ia) THEN
        !
        ! Connect atom type with species-dependent quantities...
        !
        ibs=ityp(ib)
        !
        ! Transfer species-specific cutoff to spcutB...
        !
        spcutB=spdata(ibs,1)
        !
        ! Compute distance between atom A and atom B (according to the minimum image convention)...
        !
        dAB(1)=atxyz(1,ia)-atxyz(1,ib)   ! r_AB = r_A - r_B   
        dAB(2)=atxyz(2,ia)-atxyz(2,ib)   ! r_AB = r_A - r_B   
        dAB(3)=atxyz(3,ia)-atxyz(3,ib)   ! r_AB = r_A - r_B   
        !
        dAB2(1)=ainv(1,1)*dAB(1)+ainv(1,2)*dAB(2)+ainv(1,3)*dAB(3)   ! s_AB = h^-1 r_AB
        dAB2(2)=ainv(2,1)*dAB(1)+ainv(2,2)*dAB(2)+ainv(2,3)*dAB(3)   ! s_AB = h^-1 r_AB
        dAB2(3)=ainv(3,1)*dAB(1)+ainv(3,2)*dAB(2)+ainv(3,3)*dAB(3)   ! s_AB = h^-1 r_AB
        !
        dAB2(1)=dAB2(1)-IDNINT(dAB2(1))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        dAB2(2)=dAB2(2)-IDNINT(dAB2(2))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        dAB2(3)=dAB2(3)-IDNINT(dAB2(3))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        !
        dAB(1)=h(1,1)*dAB2(1)+h(1,2)*dAB2(2)+h(1,3)*dAB2(3)   ! r_AB = h s_AB (MIC)
        dAB(2)=h(2,1)*dAB2(1)+h(2,2)*dAB2(2)+h(2,3)*dAB2(3)   ! r_AB = h s_AB (MIC)
        dAB(3)=h(3,1)*dAB2(1)+h(3,2)*dAB2(2)+h(3,3)*dAB2(3)   ! r_AB = h s_AB (MIC)
        !
        dABmic(ib)=DSQRT(dAB(1)*dAB(1)+dAB(2)*dAB(2)+dAB(3)*dAB(3))   ! |r_A - r_B| (MIC)
        !
        IF(dABmic(ib).LT.(spcutA+spcutB)) THEN
          !
          nj(ia)=nj(ia)+1
          overlap(nj(ia),ia)=ib
          !
          IF(nj(ia).EQ.1) THEN
            !
            overlap(nj(ia),ia)=ib
            !
          ELSE IF(dABmic(overlap(nj(ia)-1,ia)).LE.dABmic(ib)) THEN
            !
            overlap(nj(ia),ia)=ib
            !
          ELSE
            !
            overlap2(:)=0
            !
            DO ir=1,nj(ia)-1
              !
              IF(dABmic(overlap(ir,ia)).LT.dABmic(ib)) THEN
                !
                overlap2(ir)=overlap(ir,ia)
                !
              ELSE
                !
                overlap2(ir)=ib
                !
                DO ip=ir+1,nj(ia)
                  !
                  overlap2(ip)=overlap(ip-1,ia)
                  !
                END DO
                !
                GO TO 30
                !
              END IF
              !
            END DO !ir
            !
30          CONTINUE
            !
            DO ir=1,nj(ia)
              !
              overlap(ir,ia)=overlap2(ir)
              !
            END DO
            !
          END IF !nj(ia)
          !
        END IF !dABmic(j)
        !
      END IF !ia/=ib
      !
    END DO !ib 
    !
  END DO !ia
  !
  IF (ALLOCATED(dABmic))   DEALLOCATE(dABmic)
  IF (ALLOCATED(overlap2)) DEALLOCATE(overlap2)
  !
  ! Now form unique atom pair overlap matrix...
  !
  nbmax=nat
  !
  ALLOCATE(pair(nbmax,nat)); pair=0
  ALLOCATE(npair(nat)); npair=0
  !
  num=0; num1=0
  !
  DO j=1,nbmax
    !
    DO ia=1,nat
      !
      DO jj=1,nj(ia)
        !
        jj_neib_of_i=overlap(jj,ia)
        !
        IF(jj_neib_of_i.GT.0) THEN
           !
           pair(j,ia)=jj_neib_of_i
           overlap(jj,ia)=0
           num=num+1
           !
           DO k=1,nj(jj_neib_of_i)
             !
             IF(overlap(k,jj_neib_of_i).EQ.ia) THEN
               !
               overlap(k,jj_neib_of_i)=0
               num1=num1+1
               !
               GO TO 40
               !
             END IF
             !
           END DO !k
           !
        END IF
        !
      END DO !jj
      !
40    CONTINUE
      !
    END DO !ia
    !
  END DO !j
  !
  IF(num.NE.num1) THEN
    !
    CALL errore('tsvdw','ERROR: num .NE. num1...',1)
    !
  END IF
  !
  ! Count number of unique atom pairs for each atom...
  !
  DO ia=1,nat
      !
      num=0
      !
      DO j=1,nbmax
        !
        IF(pair(j,ia).NE.0) num=num+1
        !
      END DO
      !
      npair(ia)=num
      !
  END DO
  !
  IF (ALLOCATED(overlap))   DEALLOCATE(overlap)
  IF (ALLOCATED(nj))        DEALLOCATE(nj)
  !
  CALL stop_clock('tsvdw_pair')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_unique_pair
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_rhotot( rhor )
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rhor(:,:)
  !
  ! Local variables
  !
  INTEGER :: ir,ierr,nspin
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rhor_tmp1,rhor_tmp2
  !  
  CALL start_clock('tsvdw_rhotot')
  !
  ! Initialization of rhotot array (local copy of the real-space charge density)...
  !
  ALLOCATE(rhotot(nr1*nr2*nr3)); rhotot=0.0_DP
  nspin = SIZE(rhor,2)
  IF ( nspin < 1 .OR.  nspin > 2 ) CALL errore ('tsvdw','invalid nspin',1)
#if defined(__MPI)
  !
  ! Initialization of rhor_tmp temporary buffers...
  !
  ALLOCATE(rhor_tmp1(nr1*nr2*nr3)); rhor_tmp1=0.0_DP
  !
  IF (nspin.EQ.2) THEN
    !
    ALLOCATE(rhor_tmp2(nr1*nr2*nr3)); rhor_tmp2=0.0_DP
    !
  END IF
  !
  ! Collect distributed rhor and broadcast to all processors...
  !
  DO iproc=1,nproc_bgrp
    !
    recvcount(iproc)=dfftp%npp(iproc)*nr1*nr2
    !
  END DO
  !
  rdispls(1) = 0
  !
  DO iproc=2,nproc_bgrp
    !
    rdispls(iproc)=rdispls(iproc-1)+recvcount(iproc-1)
    !
  END DO
  !
  CALL MPI_ALLGATHERV(rhor(1,1),dfftp%npp(me_bgrp+1)*nr1*nr2,&
      MPI_DOUBLE_PRECISION,rhor_tmp1(1),recvcount,rdispls,&
      MPI_DOUBLE_PRECISION,intra_bgrp_comm,ierr)
  !
  IF (nspin.EQ.2) THEN
    !
    CALL MPI_ALLGATHERV(rhor(1,2),dfftp%npp(me_bgrp+1)*nr1*nr2,&
        MPI_DOUBLE_PRECISION,rhor_tmp2(1),recvcount,rdispls,&
        MPI_DOUBLE_PRECISION,intra_bgrp_comm,ierr)
    !
  END IF
  ! 
  ! Transfer rhor temporary arrays to rhotot array...
  !
  rhotot=rhor_tmp1
  !
  IF (nspin.EQ.2) THEN
    !
    rhotot=rhotot+rhor_tmp2
    !
  END IF
  !
  ! Clean-up temporary arrays...
  !
  IF (ALLOCATED(rhor_tmp1))     DEALLOCATE(rhor_tmp1)
  IF (ALLOCATED(rhor_tmp2))     DEALLOCATE(rhor_tmp2)
  !
#else
  rhotot(:) = rhor(:,1)
  IF (nspin == 2) rhotot(:) = rhotot(:) + rhor(:,2)
#endif
  CALL stop_clock('tsvdw_rhotot')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_rhotot
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_screen()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,ias,Ntmp,Ntmpr,Npts,Nptsr,ir1,ir2,ir3,Ndim,Nints,off1,off1r,ioff,boff,iq,ir
  REAL(DP) :: spcutA,spdA,dq(3),dqA(3),dqAs(3),dk1,rhoA
  REAL(DP), ALLOCATABLE :: dqAmic(:,:,:),dveffAdntmp(:,:,:)
  !
  CALL start_clock('tsvdw_screen')
  !
  ! Allocate and initialize gomegar array which contains (in a bit array) which atoms contribute to a given point 
  ! on the reduced real-space grid for all atoms...
  !
  Nints=nat/bsint+1
  ALLOCATE(gomegar(nr1r*nr2r*nr3r,Nints)); gomegar=0
  !
  ! Allocate and initialize NsomegaA and NsomegaAr arrays which contains the number of points in the 
  ! full and reduced spherical atomic integration domains for all atoms...
  !
  ALLOCATE(NsomegaA(nat)); NsomegaA=0
  ALLOCATE(NsomegaAr(nat)); NsomegaAr=0
  !
  ! Allocate and initialize somegaA and somegaAr arrays which contains the grid indices (along nr1,nr2,nr3) for each point in the 
  ! full and reduced spherical atomic integration domains for each atom assigned to a given processor...
  !
  Ntmp=INT(1.10_DP*((4.0_DP/3.0_DP)*pi*spcutAmax**(3.0_DP)/omega)*(nr1*nr2*nr3))    ! Number of points in the full sphere + 10% buffer (to be safe)
  Ntmpr=INT(1.10_DP*((4.0_DP/3.0_DP)*pi*spcutAmax**(3.0_DP)/omega)*(nr1r*nr2r*nr3r)) ! Number of points in the reduced sphere + 10% buffer (to be safe)
  Ndim=MAX(1,nstates(me))
  ALLOCATE(somegaA(Ntmp,3,Ndim)); somegaA=0
  ALLOCATE(somegaAr(Ntmpr,3,Ndim)); somegaAr=0
  !
  ! Allocate and initialize gomegaAr array which contains in a bit array all of the atoms that intersect with each point in the 
  ! reduced spherical atomic integration domain for each atom assigned to a given processor...
  !
  ALLOCATE(gomegaAr(Ntmpr,Nints,Ndim)); gomegaAr=0
  !
  ! Initialization of rhosad(r)...
  !
  ALLOCATE(rhosad(nr1*nr2*nr3)); rhosad=0.0_DP
  !
  ! Initialization of dVA/dn(r)...
  !
  ALLOCATE(dveffAdn(Ntmp,Ndim)); dveffAdn=0.0_DP
  !
  DO iproc=1,nstates(me)
    !
    ! Connect processor number with atom...
    !
    ia=me+nproc_image*(iproc-1)
    !
    ! Connect atom type with species-dependent quantities...
    !
    ias=ityp(ia)
    !
    ! Transfer species-specific cutoff to spcutA...
    !
    spcutA=spdata(ias,1)
    !
    ! Precompute inverse of species-specific linear grid spacing (replaces / with * inside inner loop)...
    !
    spdA=1.0_DP/spdata(ias,2)
    !
    ! Loop over grid points and determine if they belong to spherical atomic integration domain (if r < RcutA)...
    !
    Npts=0; Nptsr=0
    !
    ALLOCATE(dqAmic(nr1,nr2,nr3)); dqAmic=0.0_DP
    ALLOCATE(dveffAdntmp(nr1,nr2,nr3)); dveffAdntmp=0.0_DP
    !
!$omp parallel do private(dq,dqA,dqAs,ir,dk1,rhoA,off1,ioff,boff,off1r)
    DO ir1=1,nr1
      !
      dq(1)=DBLE(ir1-1)/DBLE(nr1) ! s_i(1)
      !
      DO ir2=1,nr2
        !
        dq(2)=DBLE(ir2-1)/DBLE(nr2) ! s_i(2)
        !
        DO ir3=1,nr3
          !
          dq(3)=DBLE(ir3-1)/DBLE(nr3) ! s_i(3)
          !
          ! Compute distance between grid point and atom according to minimum image convention (MIC)...
          !
          dqA(1)=h(1,1)*dq(1)+h(1,2)*dq(2)+h(1,3)*dq(3)   ! r_i = h s_i
          dqA(2)=h(2,1)*dq(1)+h(2,2)*dq(2)+h(2,3)*dq(3)   ! r_i = h s_i
          dqA(3)=h(3,1)*dq(1)+h(3,2)*dq(2)+h(3,3)*dq(3)   ! r_i = h s_i
          !
          dqA(1)=dqA(1)-atxyz(1,ia)   ! r_iA = r_i - r_A
          dqA(2)=dqA(2)-atxyz(2,ia)   ! r_iA = r_i - r_A
          dqA(3)=dqA(3)-atxyz(3,ia)   ! r_iA = r_i - r_A
          !
          dqAs(1)=ainv(1,1)*dqA(1)+ainv(1,2)*dqA(2)+ainv(1,3)*dqA(3)   ! s_iA = h^-1 r_iA
          dqAs(2)=ainv(2,1)*dqA(1)+ainv(2,2)*dqA(2)+ainv(2,3)*dqA(3)   ! s_iA = h^-1 r_iA
          dqAs(3)=ainv(3,1)*dqA(1)+ainv(3,2)*dqA(2)+ainv(3,3)*dqA(3)   ! s_iA = h^-1 r_iA
          !
          dqAs(1)=dqAs(1)-IDNINT(dqAs(1))   ! impose MIC on s_iA in range: [-0.5,+0.5]
          dqAs(2)=dqAs(2)-IDNINT(dqAs(2))   ! impose MIC on s_iA in range: [-0.5,+0.5]
          dqAs(3)=dqAs(3)-IDNINT(dqAs(3))   ! impose MIC on s_iA in range: [-0.5,+0.5]
          !
          dqA(1)=h(1,1)*dqAs(1)+h(1,2)*dqAs(2)+h(1,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
          dqA(2)=h(2,1)*dqAs(1)+h(2,2)*dqAs(2)+h(2,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
          dqA(3)=h(3,1)*dqAs(1)+h(3,2)*dqAs(2)+h(3,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
          !
          dqAmic(ir1,ir2,ir3)=DSQRT(dqA(1)*dqA(1)+dqA(2)*dqA(2)+dqA(3)*dqA(3))   ! |r_i - r_A| (MIC)
          !
          ! Screen grid point according to atomic radial cutoff...
          !
          IF (dqAmic(ir1,ir2,ir3).LE.spcutA) THEN
            !
            ! Form rhosad(r) on the real-space mesh...
            ! N.B. This algorithm only works when the images of a given atom are greater than the radial grid cutoff values for ALL atomic species...
            !
            ! Determine the index in the atomic linear equispaced grid such that grd(ir) <= dqA <= grd(ir+1) and distance between dqA and grd(ir)...
            !
            ir=INT(dqAmic(ir1,ir2,ir3)*spdA)
            dk1=dqAmic(ir1,ir2,ir3)-spgrd(ias,ir)
            !
            ! Perform linear interpolation to obtain the value of the atomic pseudo-density at the given grid point...
            !
            rhoA=LIA(ias,ir)+LIB(ias,ir)*dk1
            !
            ! Increment contribution to rhosad(r)...
            !
            off1=ir1+(ir2-1)*nr1+(ir3-1)*nr1*nr2    !global offset [nr1,nr2,nr3]
            rhosad(off1)=rhosad(off1)+rhoA
            !
            ! Form numerator of dVA/dn(r) only...
            !
            dveffAdntmp(ir1,ir2,ir3)=dqAmic(ir1,ir2,ir3)**(3.0_DP)*rhoA
            !
            ! On reduced grid only, form screened somegaAr and gomegar...
            !
            IF ((MOD(ir1,2).EQ.1).AND.(MOD(ir2,2).EQ.1).AND.(MOD(ir3,2).EQ.1))  THEN
              !
              ioff=((ia-1)/bsint)+1                                    ! integer offset for gomegar bit array
              boff=(ia-((ioff-1)*bsint))-1                             ! bit offset for gomegar bit array
              off1r=(ir1+1)/2+((ir2-1)/2)*nr1r+((ir3-1)/2)*nr1r*nr2r   ! reduced global offset [nr1r,nr2r,nr3r]
              !
              gomegar(off1r,ioff)=IBSET(gomegar(off1r,ioff),boff)
              !
            END IF
            !
          END IF
          !
        END DO !ir3
        !
      END DO !ir2
      !
    END DO !ir1
!$omp end parallel do 
    !
    DO ir1=1,nr1
      !
      DO ir2=1,nr2
        !
        DO ir3=1,nr3
          !
          ! Screen grid point according to atomic radial cutoff...
          !
          IF (dqAmic(ir1,ir2,ir3).LE.spcutA) THEN
            !
            Npts=Npts+1
            !
            ! Form screened somegaA...
            !
            somegaA(Npts,1,iproc)=ir1
            somegaA(Npts,2,iproc)=ir2
            somegaA(Npts,3,iproc)=ir3
            !
            dveffAdn(Npts,iproc)=dveffAdntmp(ir1,ir2,ir3)
            !
            ! On reduced grid only, form screened somegaAr ... 
            !
            IF ((MOD(ir1,2).EQ.1).AND.(MOD(ir2,2).EQ.1).AND.(MOD(ir3,2).EQ.1))  THEN
              !
              Nptsr=Nptsr+1
              !
              ! Form reduced screened somegaAr...
              !
              somegaAr(Nptsr,1,iproc)=ir1
              somegaAr(Nptsr,2,iproc)=ir2
              somegaAr(Nptsr,3,iproc)=ir3
              !
            END IF
            !
          END IF

        END DO !ir3
        !
      END DO !ir2
      !
    END DO !ir1
    !
    NsomegaA(ia)=Npts
    NsomegaAr(ia)=Nptsr
    !
    IF (ALLOCATED(dqAmic))    DEALLOCATE(dqAmic)
    IF (ALLOCATED(dveffAdntmp))    DEALLOCATE(dveffAdntmp)
    ! 
  END DO ! iproc
  !
  ! Collect NsomegaA, NsomegaAr, gomegar, and rhosad over all processors and broadcast...
  !
  CALL mp_sum(NsomegaA,intra_image_comm)
  CALL mp_sum(NsomegaAr,intra_image_comm)
  CALL mp_sum(gomegar,intra_image_comm)
  CALL mp_sum(rhosad,intra_image_comm)
  !
  ! Decompose gomegar to gomegaAr to save on memory storage...
  !
  DO iproc=1,nstates(me)
    !
    ! Connect processor number with atom...
    !
    ia=me+nproc_image*(iproc-1)
    !
    ! Loop over points in the (pre-screened) reduced spherical atomic integration domain...
    !
!$omp parallel do private(off1r,ir)
    DO iq=1,NsomegaAr(ia)
      !
      DO ir=1,Nints
        !
        off1r=(somegaAr(iq,1,iproc)+1)/2+((somegaAr(iq,2,iproc)-1)/2)*nr1r+((somegaAr(iq,3,iproc)-1)/2)*nr1r*nr2r   ! reduced global offset [nr1r,nr2r,nr3r]
        gomegaAr(iq,ir,iproc)=gomegar(off1r,ir)
        !
      END DO
      !
    END DO !iq
!$omp end parallel do 
    !
  END DO ! iproc
  !
  ! Clean-up temporary arrays...
  !
  IF (ALLOCATED(gomegar))    DEALLOCATE(gomegar)
  !
  CALL stop_clock('tsvdw_screen')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_screen
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_veff()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,iq,off1
  REAL(DP) :: normr
  !
  CALL start_clock('tsvdw_veff')
  !
  ! Initialization of effective volume...
  !
  ALLOCATE(veff(nat)); veff=0.0_DP
  !
  ! Normalization factor for veff integral...
  !
  normr=omega/DBLE(nr1r*nr2r*nr3r)
  !
  ! Loop over atoms in the simulation cell...
  !
  DO iproc=1,nstates(me)
    !
    ! Connect processor number with atom...
    !
    ia=me+nproc_image*(iproc-1)
    !
    ! Loop over points in the (pre-screened) spherical atomic integration domain...
    !
!$omp parallel do private(off1),reduction(+:veff)
    DO iq=1,NsomegaA(ia)
      !
      ! Compute veff integrand and complete dispersion potential (functional derivative of veff(A) wrt charge density)...
      !
      !        veff(A) = INT [|r-rA|^3*rhoA(|r-rA|)*rhotot(r)/rhosad(r)]
      !
      !       dveff(A)/dn(r) = |r-rA|^3*rhoA(|r-rA|)/rhosad(r)
      !
      off1=somegaA(iq,1,iproc)+(somegaA(iq,2,iproc)-1)*nr1+(somegaA(iq,3,iproc)-1)*nr1*nr2    !global offset [nr1,nr2,nr3]
      dveffAdn(iq,iproc)=dveffAdn(iq,iproc)/rhosad(off1)
      !
      ! Increment veff...
      !
      IF ((MOD(somegaA(iq,1,iproc),2).EQ.1).AND.(MOD(somegaA(iq,2,iproc),2).EQ.1).AND.(MOD(somegaA(iq,3,iproc),2).EQ.1))  THEN
        !
        veff(ia)=veff(ia)+(dveffAdn(iq,iproc)*rhotot(off1))
        !
      END IF
      !
    END DO !iq
!$omp end parallel do 
    !
    ! Apply final normalization to veff integral...
    !
    veff(ia)=normr*veff(ia)
    !
  END DO !iproc
  !
  ! Collect veff over all processors and broadcast...
  !
  CALL mp_sum(veff,intra_image_comm)
  !
  VefftsvdW = veff
  !
  CALL stop_clock('tsvdw_veff')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_veff
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_dveff()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,ib,ias,ibs,iq,ir,i,j,ipair,off1,ioff,boff
  REAL(DP) :: spcutA,spcutB,spdA,spdB,dq(3),dqA(3),dqAs(3),dqB(3),dqBs(3),dqAmic,dqBmic,dABmic,normr
  REAL(DP) :: dk1,dptmp1,dptmp2,dptmp3,dptmp4,dptmp5,rhoA,rhoB,drhoA,drhoB,dVAdRA(3),dVAdRB(3),dVBdRA(3)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: predVAdRB
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dqxyzr,dqAxyzs,predVBdRA
  !
  CALL start_clock('tsvdw_dveff')
  !
  ! Initialization of the dveff/dR and dveff/dh arrays...
  !
  ALLOCATE(dveffdR(nat,nat,3)); dveffdR=0.0_DP
  ALLOCATE(dveffdh(nat,3,3)); dveffdh=0.0_DP
  !
  ! Normalization factor for dveff integrals...
  !
  normr=omega/DBLE(nr1r*nr2r*nr3r)
  !
  ! Loop over atoms A in the simulation cell and compute dveffdR and dveffdh...
  !
  DO iproc=1,nstates(me)
    !
    ! Connect processor number with atom...
    !
    ia=me+nproc_image*(iproc-1)
    !
    ! Connect atom type with species-dependent quantities...
    !
    ias=ityp(ia) 
    !
    ! Transfer species-specific cutoff to spcutA...
    !
    spcutA=spdata(ias,1)
    !
    ! Precompute inverse of species-specific linear grid spacing (replaces / with * during interpolation)...
    !
    spdA=1.0_DP/spdata(ias,2)
    !
    ! Allocate and initialize atom-specific arrays...
    !
    ALLOCATE(dqxyzr(NsomegaAr(ia),3)); dqxyzr=0.0_DP
    ALLOCATE(dqAxyzs(NsomegaAr(ia),3)); dqAxyzs=0.0_DP
    ALLOCATE(predVAdRB(NsomegaAr(ia))); predVAdRB=0.0_DP
    ALLOCATE(predVBdRA(NsomegaAr(ia),3)); predVBdRA=0.0_DP
    !
    ! Initial loop over points in the (pre-screened) reduced spherical atomic integration domain for atom A to compute 
    ! self-derivative (dV(A)/dR(A)), quantities necessary for dV(A)/dR(B) and dV(B)/dR(A), and self-contribution to dV(A)/dh...
    !
!$omp parallel do private(dq,dqA,dqAs,dqAmic,ir,dk1,rhoA,drhoA, &
!$omp off1,dptmp1,dptmp2,dptmp3,dptmp4,dVAdRA,dptmp5,i,j), &
!$omp reduction(-:dveffdh),reduction(+:dveffdR) 
    DO iq=1,NsomegaAr(ia)
      !
      ! Compute global/cell reference frame Cartesian coordinates of given real-space grid point...
      !
      dq(1)=DBLE(somegaAr(iq,1,iproc)-1)/DBLE(nr1)   ! s_i(1)
      dq(2)=DBLE(somegaAr(iq,2,iproc)-1)/DBLE(nr2)   ! s_i(2)
      dq(3)=DBLE(somegaAr(iq,3,iproc)-1)/DBLE(nr3)   ! s_i(3)
      !
      dqA(1)=h(1,1)*dq(1)+h(1,2)*dq(2)+h(1,3)*dq(3)   ! r_i = h s_i
      dqA(2)=h(2,1)*dq(1)+h(2,2)*dq(2)+h(2,3)*dq(3)   ! r_i = h s_i
      dqA(3)=h(3,1)*dq(1)+h(3,2)*dq(2)+h(3,3)*dq(3)   ! r_i = h s_i
      !
      ! Accumulate the Cartesian coordinates of the real-space grid point in the dqxyzr array:
      !
      !   dqxyzr(:,1) := x-coordinate of grid point (global/cell reference frame)
      !   dqxyzr(:,2) := y-coordinate of grid point (global/cell reference frame)
      !   dqxyzr(:,3) := z-coordinate of grid point (global/cell reference frame)
      !
      dqxyzr(iq,1)=dqA(1)
      dqxyzr(iq,2)=dqA(2)
      dqxyzr(iq,3)=dqA(3)
      !
      ! Compute distance between grid point and atom in scaled coordinates (s_iA) according to minimum image convention (MIC)...
      !
      dqA(1)=dqA(1)-atxyz(1,ia)   ! r_iA = r_i - r_A
      dqA(2)=dqA(2)-atxyz(2,ia)   ! r_iA = r_i - r_A
      dqA(3)=dqA(3)-atxyz(3,ia)   ! r_iA = r_i - r_A
      !
      dqAs(1)=ainv(1,1)*dqA(1)+ainv(1,2)*dqA(2)+ainv(1,3)*dqA(3)   ! s_iA = h^-1 r_iA
      dqAs(2)=ainv(2,1)*dqA(1)+ainv(2,2)*dqA(2)+ainv(2,3)*dqA(3)   ! s_iA = h^-1 r_iA
      dqAs(3)=ainv(3,1)*dqA(1)+ainv(3,2)*dqA(2)+ainv(3,3)*dqA(3)   ! s_iA = h^-1 r_iA
      !
      dqAs(1)=dqAs(1)-IDNINT(dqAs(1))   ! impose MIC on s_iA in range: [-0.5,+0.5]
      dqAs(2)=dqAs(2)-IDNINT(dqAs(2))   ! impose MIC on s_iA in range: [-0.5,+0.5]
      dqAs(3)=dqAs(3)-IDNINT(dqAs(3))   ! impose MIC on s_iA in range: [-0.5,+0.5]
      !
      ! Accumulate the components of the s_i - s_A vector in the dqAxyzs array:
      !
      !   dqAxyzs(:,1) := 1-coordinate of s_i - s_A vector (local/atom reference frame)
      !   dqAxyzs(:,2) := 2-coordinate of s_i - s_A vector (local/atom reference frame)
      !   dqAxyzs(:,3) := 3-coordinate of s_i - s_A vector (local/atom reference frame)
      !
      dqAxyzs(iq,1)=dqAs(1)
      dqAxyzs(iq,2)=dqAs(2)
      dqAxyzs(iq,3)=dqAs(3)
      !
      ! Convert MIC distance components from scaled coordinates to Cartesian coordinates (s_iA -> r_iA)...
      !
      dqA(1)=h(1,1)*dqAs(1)+h(1,2)*dqAs(2)+h(1,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
      dqA(2)=h(2,1)*dqAs(1)+h(2,2)*dqAs(2)+h(2,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
      dqA(3)=h(3,1)*dqAs(1)+h(3,2)*dqAs(2)+h(3,3)*dqAs(3)   ! r_iA = h s_iA (MIC)
      !
      dqAmic=DSQRT(dqA(1)*dqA(1)+dqA(2)*dqA(2)+dqA(3)*dqA(3))   ! |r_i - r_A| (MIC)
      !
      ! Determine the index in the atomic linear equispaced grid such that grd(ir) <= dqA <= grd(ir+1) and distance between dqA and grd(ir)...
      !
      ir=INT(dqAmic*spdA)
      dk1=dqAmic-spgrd(ias,ir)
      !
      ! Perform linear interpolation to obtain the value of the atomic pseudo-density and its derivative at the given grid point...
      !
      rhoA=LIA(ias,ir)+LIB(ias,ir)*dk1        !rhoA at grid point via linear interpolation
      drhoA=dLIA(ias,ir)+dLIB(ias,ir)*dk1     !drhoA at grid point via linear interpolation
      !
      ! Compute global offset for rhosad(r) and rhotot(r), both computed on the real-space mesh...
      !
      off1=somegaAr(iq,1,iproc)+(somegaAr(iq,2,iproc)-1)*nr1+(somegaAr(iq,3,iproc)-1)*nr1*nr2    !global offset [nr1,nr2,nr3]
      !
      ! Compute self-derivative dVA/dpA integrand for p={x,y,z}...
      !
      !   dVA/dpA = INT {(p-pA)*|r-rA|*rhotot(r)/rhosad(r)*[|r-rA|*rho(|r-rA|)*drho(|r-rA|)/rhosad(r)-|r-rA|*drho(|r-rA|)-3*rho(|r-rA|)]}
      !
      dptmp1=1.0_DP/rhosad(off1)
      dptmp2=dqAmic*drhoA
      dptmp3=dqAmic*rhotot(off1)*dptmp1
      dptmp4=((rhoA*dptmp1-1.0_DP)*dptmp2-3.0_DP*rhoA)*dptmp3
      !
      dVAdRA=dqA*dptmp4  !dVA/dpA integrand/contribution for the given grid point...
      !
      ! Increment self-derivative dVA/dpA for p={x,y,z}...
      !
      DO i=1,3
        !
        dveffdR(ia,ia,i)=dveffdR(ia,ia,i)+dVAdRA(i)
        !
      END DO !i
      !
      ! Increment self-contribution to dVA/dhpq for p,q={x,y,z}...
      !
      !   dVA/dhpq <-- INT {-(p-pA)*(qs-qsA)*|r-rA|*rhotot(r)/rhosad(r)*[|r-rA|*rho(|r-rA|)*drho(|r-rA|)/rhosad(r)-|r-rA|*drho(|r-rA|)-3*rho(|r-rA|)]}
      !
      DO i=1,3
        !
        DO j=1,3
          !
          dveffdh(ia,i,j)=dveffdh(ia,i,j)-dVAdRA(i)*dqAxyzs(iq,j)
          !
        END DO !j
        !
      END DO !i
      !
      ! Precompute quantities necessary for dV(A)/dR(B) and dV(B)/dR(A)...
      !
      predVAdRB(iq)=dptmp1*rhoA*dqAmic*dqAmic*dptmp3
      !
      dptmp5=dptmp1*dptmp1*drhoA*rhotot(off1)
      !
      IF (dqAmic.LT.(1.0E-12_DP)) THEN
        !
        predVBdRA(iq,:)=dptmp5
        !
      ELSE
        !
        predVBdRA(iq,:)=dptmp5*dqA(:)/dqAmic
        !
      END IF
      !
    END DO !iq
!$omp end parallel do 
    !
    ! Inner loop over unique atom pairs B in the simulation cell to compute pair contributions to dveffdR and dveffdh...
    !
!$omp parallel do private(dqB,dqBs,dqBmic,ir,dk1,rhoB,drhoB,dVAdRB,dVBdRA, &
!$omp i,j,ib,ibs,spcutB,spdB,ioff,boff), &
!$omp reduction(+:dveffdR),reduction(-:dveffdh) 
    DO ipair=1,npair(ia)
      !
      ! Connect pair number with atom...
      !
      ib=pair(ipair,ia)
      !
      ! Connect atom type with species-dependent quantities...
      !
      ibs=ityp(ib) 
      !
      ! Transfer species-specific cutoff to spcutB...
      !
      spcutB=spdata(ibs,1)
      !
      ! Precompute inverse of species-specific linear grid spacing (replaces / with * during interpolation)...
      !
      spdB=1.0_DP/spdata(ibs,2)
      !
      ! Determine atom B offsets for using gomegaAr bit array screening below...
      !
      ioff=((ib-1)/bsint)+1                     ! integer offset for gomegaAr bit array
      boff=(ib-((ioff-1)*bsint))-1              ! bit offset for gomegaAr bit array
      !
      ! Inner loop over points in the (pre-screened) reduced spherical atomic integration domain for atom A to compute 
      ! non-self-derivatives (dV(A)/dR(B) and dV(B)/dR(A)) and non-self-contributions to dV(A)/dh and dV(B)/dh in the overlapping integration domain...
      !
      DO iq=1,NsomegaAr(ia)
        !
        ! Determine if atom B contributes to the given point on the reduced spherical atomic integration domain on atom A (using gomegaAr bit array)...
        !
        IF (BTEST(gomegaAr(iq,ioff,iproc),boff)) THEN
          !
          ! Compute distance between grid point and atom B according to minimum image convention (MIC)...
          !
          dqB(1)=dqxyzr(iq,1)-atxyz(1,ib)   ! r_iB = r_i - r_B
          dqB(2)=dqxyzr(iq,2)-atxyz(2,ib)   ! r_iB = r_i - r_B
          dqB(3)=dqxyzr(iq,3)-atxyz(3,ib)   ! r_iB = r_i - r_B
          !
          dqBs(1)=ainv(1,1)*dqB(1)+ainv(1,2)*dqB(2)+ainv(1,3)*dqB(3)   ! s_iB = h^-1 r_iB
          dqBs(2)=ainv(2,1)*dqB(1)+ainv(2,2)*dqB(2)+ainv(2,3)*dqB(3)   ! s_iB = h^-1 r_iB
          dqBs(3)=ainv(3,1)*dqB(1)+ainv(3,2)*dqB(2)+ainv(3,3)*dqB(3)   ! s_iB = h^-1 r_iB
          !
          dqBs(1)=dqBs(1)-IDNINT(dqBs(1))   ! impose MIC on s_iB in range: [-0.5,+0.5]
          dqBs(2)=dqBs(2)-IDNINT(dqBs(2))   ! impose MIC on s_iB in range: [-0.5,+0.5]
          dqBs(3)=dqBs(3)-IDNINT(dqBs(3))   ! impose MIC on s_iB in range: [-0.5,+0.5]
          !
          dqB(1)=h(1,1)*dqBs(1)+h(1,2)*dqBs(2)+h(1,3)*dqBs(3)   ! r_iB = h s_iB (MIC)
          dqB(2)=h(2,1)*dqBs(1)+h(2,2)*dqBs(2)+h(2,3)*dqBs(3)   ! r_iB = h s_iB (MIC)
          dqB(3)=h(3,1)*dqBs(1)+h(3,2)*dqBs(2)+h(3,3)*dqBs(3)   ! r_iB = h s_iB (MIC)
          !
          dqBmic=DSQRT(dqB(1)*dqB(1)+dqB(2)*dqB(2)+dqB(3)*dqB(3))   ! |r_i - r_B| (MIC)
          !
          ! Final screening based on the (pre-screened) spherical atomic integration domain on atom B...
          !
          IF (dqBmic.LE.spcutB) THEN
            !
            ! Determine the index in the atomic linear equispaced grid such that grd(ir) <= dqB <= grd(ir+1) and distance between dqB and grd(ir)...
            !
            ir=INT(dqBmic*spdB)
            dk1=dqBmic-spgrd(ibs,ir)
            !
            ! Perform linear interpolation to obtain the value of the atomic pseudo-density and its derivative at the given grid point...
            !
            rhoB=LIA(ibs,ir)+LIB(ibs,ir)*dk1         !rhoB at grid point via linear interpolation
            drhoB=dLIA(ibs,ir)+dLIB(ibs,ir)*dk1      !drhoB at grid point via linear interpolation
            !
            ! Compute dVA/dpB integrand for p={x,y,z}...
            !
            !   dVA/dpB = INT {(p-pB)/|r-rB|*[drho(|r-rB|)*|r-rA|^3*rho(|r-rA|)*rhotot(r)/rhosad(r)^2]}
            !
            IF (dqBmic.LT.(1.0E-12_DP)) THEN
              !
              dVAdRB(:)=predVAdRB(iq)*drhoB
              !
            ELSE
              !
              dVAdRB(:)=predVAdRB(iq)*drhoB*dqB(:)/dqBmic
              !
            END IF
            !
            ! Increment non-self-derivative dVA/dpB for p={x,y,z}...
            !
            DO i=1,3
              !
              dveffdR(ia,ib,i)=dveffdR(ia,ib,i)+dVAdRB(i)
              !
            END DO !i
            !
            ! Increment non-self-contribution to dVA/dhpq for p,q={x,y,z} from atom B...
            !
            !   dVA/dhpq <-- INT {-(p-pB)*(qs-qsB)/|r-rB|*[drho(|r-rB|)*|r-rA|^3*rho(|r-rA|)*rhotot(r)/rhosad(r)^2]}
            !
            DO i=1,3
              !
              DO j=1,3
                !
                dveffdh(ia,i,j)=dveffdh(ia,i,j)-dVAdRB(i)*dqBs(j)
                !
              END DO !j
              !
            END DO !i
            !
            ! Compute dVB/dpA integrand for p={x,y,z}...
            !
            !   dVB/dpA = INT {(p-pA)/|r-rA|*[drho(|r-rA|)*|r-rB|^3*rho(|r-rB|)*rhotot(r)/rhosad(r)^2]}
            !
            dVBdRA(:)=predVBdRA(iq,:)*rhoB*dqBmic*dqBmic*dqBmic
            !
            ! Increment non-self-derivative dVB/dpA for p={x,y,z} from atom A...
            !
            DO i=1,3
              !
              dveffdR(ib,ia,i)=dveffdR(ib,ia,i)+dVBdRA(i)
              !
            END DO !i
            !
            ! Increment non-self-contribution to dVB/dhpq for p,q={x,y,z} from atom A...
            !
            !   dVB/dhpq <-- INT {-(p-pA)*(qs-qsA)/|r-rA|*[drho(|r-rA|)*|r-rB|^3*rho(|r-rB|)*rhotot(r)/rhosad(r)^2]}
            !
            DO i=1,3
              !
              DO j=1,3
                !
                dveffdh(ib,i,j)=dveffdh(ib,i,j)-dVBdRA(i)*dqAxyzs(iq,j)
                !
              END DO !j
              !
            END DO !i
            !
          END IF
          !
        END IF !BTEST
        !
      END DO !iq
      !
    END DO !ipair
!$omp end parallel do 
    !
    ! Deallocate temporary arrays...
    !
    IF (ALLOCATED(dqxyzr))     DEALLOCATE(dqxyzr)
    IF (ALLOCATED(dqAxyzs))    DEALLOCATE(dqAxyzs)
    IF (ALLOCATED(predVAdRB))  DEALLOCATE(predVAdRB)
    IF (ALLOCATED(predVBdRA))  DEALLOCATE(predVBdRA) 
    !
  END DO !iproc
  !
  ! Apply final normalization of dVA/dR integrals...
  !
  dveffdR=normr*dveffdR
  !
  ! Apply final normalization of dVA/dhab integrals...
  !
  dveffdh=normr*dveffdh
  !
  ! Collect dveffdR and dveffdh over all processors and broadcast...
  !
  CALL mp_sum(dveffdR,intra_image_comm)
  CALL mp_sum(dveffdh,intra_image_comm)
  !
  CALL stop_clock('tsvdw_dveff')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_dveff
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_effqnts()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,ib,ias,ibs
  REAL(DP) :: vA,vB,num,den
  !
  ! Initialization of base effective atomic quantities...
  !
  ALLOCATE(dpeff(nat)); dpeff=0.0_DP
  ALLOCATE(R0eff(nat)); R0eff=0.0_DP
  ALLOCATE(C6AAeff(nat)); C6AAeff=0.0_DP
  ALLOCATE(C6ABeff(nat,nat)); C6ABeff=0.0_DP
  !
  ! Population of base effective atomic quantities...
  !
  DO ia=1,nat
    !
    ! Connect atom type with species-dependent quantities...
    !
    ias=ityp(ia) 
    !
    ! Precompute veff(A)/vfree(A) ratio...
    !
    vA=(veff(ia)/vfree(ias))
    !
    ! Effective atomic static dipole polarizability array...
    ! dpeff(A)=[veff(A)/vfree(A)]*dpfree(A)
    !
    dpeff(ia)=vA*dpfree(ias)
    !
    ! Effective atomic vdW radius array...
    ! R0eff(A)=[veff(A)/vfree(A)]^1/3*R0free(A)
    !
    R0eff(ia)=(vA**(1.0_DP/3.0_DP))*R0free(ias)
    !
    ! Effective homonuclear C6 coefficient array...
    ! C6AAeff(A)=[veff(A)/vfree(A)]^2*C6AAfree(A)
    !
    C6AAeff(ia)=(vA**(2.0_DP))*C6AAfree(ias)
    !
    DO ib=1,nat
      !
      ! Connect atom type with species-dependent quantities...
      !
      ibs=ityp(ib) 
      !
      ! Precompute veff(B)/vfree(B) ratio...
      !
      vB=(veff(ib)/vfree(ibs))
      !
      ! Effective heteronuclear C6 coefficient matrix...
      ! C6ABeff(A,B)=(veff(A)/vfree(A))*(veff(B)/vfree(B))*C6ABfree(A,B)
      !
      C6ABeff(ia,ib)=(vA*vB)*C6ABfree(ias,ibs) 
      !
    END DO !ib
    !
  END DO !ia
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_effqnts
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_energy()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  LOGICAL :: periodic_converged
  INTEGER :: ia,ib,ic,ias,ibs,n_period,n1,n2,n3,i,j

  REAL(DP) :: dAB(3),dAB2(3),dsAB(3),dABimg,dABimg2,dABimgn1,dABimgn2,dABimgn5,dABimgn6
  REAL(DP) :: FDV0,FDR0,FCV0,FRR0,FDV1,FCVA1,FCVB1,FDVA2,FDVB2,FDR2,FRR2,FCVA2,FCVB2,FDVi,FDRi(3),FDRii(3,3),FCVi,FRRi(3),FRRii(3,3)
  REAL(DP) :: EtsvdW_period,RAB0,edamp,fdamp,fdamp2,D1A,D1B,D2A,D2B,D12A,D12B,dptmp1,dptmp2,vtmp1(3),vtmp2(3)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: predveffAdn_period
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: FtsvdW_period,HtsvdW_period
  !
  CALL start_clock('tsvdw_energy')
  !
  ! Initialize total TS-vdW energy, ion force, and cell force ...
  !
  EtsvdW=0.0_DP; FtsvdW=0.0_DP; HtsvdW=0.0_DP
  !
  ! Allocate and initialize TS-vdW dispersion potential prefactor...
  !
  ALLOCATE(predveffAdn(nat)); predveffAdn=0.0_DP
  !
  ! Allocate and initialize periodic contributions to TS-vdW ionic forces, cell forces, and dispersion potential prefactor...
  !
  ALLOCATE(FtsvdW_period(3,nat)); FtsvdW_period=0.0_DP
  ALLOCATE(HtsvdW_period(3,3)); HtsvdW_period=0.0_DP
  ALLOCATE(predveffAdn_period(nat)); predveffAdn_period=0.0_DP
  !
  ! Precompute quantities outside all loops...
  !
  FDR0=ddamp/(2.0_DP*sR)
  FDV0=-FDR0/3.0_DP
  FCV0=0.5_DP
  FRR0=-3.0_DP
  !
  ! For periodic systems, converge the energy with respect to neighboring images...
  !
  n_period=0
  periodic_converged=.FALSE.
  !
  DO WHILE (.NOT.periodic_converged)
    !
    EtsvdW_period=0.0_DP
    FtsvdW_period=0.0_DP
    HtsvdW_period=0.0_DP
    predveffAdn_period=0.0_DP
    !
    ! Outer loop over atoms A...
    !
    DO iproc=1,nstates(me)
      !
      ! Connect processor number with atom...
      !
      ia=me+nproc_image*(iproc-1)
      !
      ! Connect atom type with species-dependent quantities...
      !
      ias=ityp(ia)
      !
      ! Precompute quantities outside loop over B...
      !
      FDV1=R0free(ias)/(vfree(ias)**(1.0_DP/3.0_DP)*veff(ia)**(2.0_DP/3.0_DP))
      FCVA1=1.0_DP/vfree(ias)
      FCVB1=veff(ia)*FCVA1
      !
      ! Inner loop over atoms B...
      !
      
!$omp parallel private(ibs,RAB0,FRR2,FDR2,FDVA2,FDVB2,FCVB2,FCVA2, &
!$omp dAB,dAB2,FDVi,FDRi,FDRii,FCVi,FRRi,FRRii,n1,n2,n3,dsAB,dABimg2, &
!$omp dABimg,dABimgn1,dABimgn2,dABimgn5,dABimgn6,edamp,fdamp,fdamp2,dptmp1, &
!$omp dptmp2,i,j,vtmp1,vtmp2,D1A,D2A,D1B,D2B,D12A,D12B,ic), &
!$omp reduction(-:EtsvdW_period),reduction(+:FtsvdW_period), &
!$omp reduction(+:HtsvdW_period),reduction(-:predveffAdn_period)
!$omp do      
      DO ib=1,nat
        !
        ! Connect atom type with species-dependent quantities...
        !
        ibs=ityp(ib)
        !
        ! Compute RAB0 as the sum of the effective vdW radii of atoms A and B...
        !
        RAB0=R0eff(ia)+R0eff(ib)
        !
        ! Precompute quantities outside loop over image cells...
        !
        FRR2=C6ABeff(ia,ib)
        FDR2=FRR2/RAB0
        FDVA2=FDR2/RAB0
        FDVB2=FDVA2*R0free(ibs)/(vfree(ibs)**(1.0_DP/3.0_DP)*veff(ib)**(2.0_DP/3.0_DP))
        FCVB2=C6ABfree(ias,ibs)/vfree(ibs)
        FCVA2=FCVB2*veff(ib)
        !
        ! Compute distance between atom A and atom B (according to the minimum image convention)...
        !
        dAB(1)=atxyz(1,ia)-atxyz(1,ib)   ! r_AB = r_A - r_B   
        dAB(2)=atxyz(2,ia)-atxyz(2,ib)   ! r_AB = r_A - r_B   
        dAB(3)=atxyz(3,ia)-atxyz(3,ib)   ! r_AB = r_A - r_B   
        !
        dAB2(1)=ainv(1,1)*dAB(1)+ainv(1,2)*dAB(2)+ainv(1,3)*dAB(3)   ! s_AB = h^-1 r_AB
        dAB2(2)=ainv(2,1)*dAB(1)+ainv(2,2)*dAB(2)+ainv(2,3)*dAB(3)   ! s_AB = h^-1 r_AB
        dAB2(3)=ainv(3,1)*dAB(1)+ainv(3,2)*dAB(2)+ainv(3,3)*dAB(3)   ! s_AB = h^-1 r_AB
        !
        dAB2(1)=dAB2(1)-IDNINT(dAB2(1))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        dAB2(2)=dAB2(2)-IDNINT(dAB2(2))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        dAB2(3)=dAB2(3)-IDNINT(dAB2(3))   ! impose MIC on s_AB in range: [-0.5,+0.5]
        !
        ! Initialize image-summed matrix elements...
        !
        FDVi=0.0_DP; FDRi=0.0_DP; FDRii=0.0_DP; FCVi=0.0_DP; FRRi=0.0_DP; FRRii=0.0_DP
        !
        ! Loop over image cells...
        !
        DO n1=-n_period,n_period
          !
          DO n2=-n_period,n_period
            !
            DO n3=-n_period,n_period
              !
              IF ((ABS(n1).EQ.n_period).OR.(ABS(n2).EQ.n_period).OR.(ABS(n3).EQ.n_period)) THEN
                !
                ! Recover MIC distance between atom A and atom B in crystal coordinates...
                !
                dsAB(1)=dAB2(1)   ! s_AB (MIC)
                dsAB(2)=dAB2(2)   ! s_AB (MIC)
                dsAB(3)=dAB2(3)   ! s_AB (MIC)
                !
                ! Increment MIC distance in crystal coordinates...
                !
                dsAB(1)=dsAB(1)+DBLE(n1)   ! s_AB (incremented, MIC only if n_period == 0)
                dsAB(2)=dsAB(2)+DBLE(n2)   ! s_AB (incremented, MIC only if n_period == 0)
                dsAB(3)=dsAB(3)+DBLE(n3)   ! s_AB (incremented, MIC only if n_period == 0)
                !
                ! Convert incremented distance back into cartesian coordinates...
                !
                dAB(1)=h(1,1)*dsAB(1)+h(1,2)*dsAB(2)+h(1,3)*dsAB(3)   ! r_AB = h s_AB (MIC only if n_period == 0)
                dAB(2)=h(2,1)*dsAB(1)+h(2,2)*dsAB(2)+h(2,3)*dsAB(3)   ! r_AB = h s_AB (MIC only if n_period == 0)
                dAB(3)=h(3,1)*dsAB(1)+h(3,2)*dsAB(2)+h(3,3)*dsAB(3)   ! r_AB = h s_AB (MIC only if n_period == 0)
                !
                ! Compute incremented distance between atom A and atom B...
                !
                dABimg2=dAB(1)*dAB(1)+dAB(2)*dAB(2)+dAB(3)*dAB(3)
                dABimg=DSQRT(dABimg2)
                !
                ! Precompute inverse powers of incremented distance between atom A and atom B...
                !
                IF ( dABimg > 0.0_dp ) THEN
                   dABimgn1=1.0_DP/dABimg
                   dABimgn2=dABimgn1*dABimgn1
                   dABimgn5=dABimgn2*dABimgn2*dABimgn1
                   dABimgn6=dABimgn5*dABimgn1
                ELSE
                   dABimgn1=0.0_DP
                   dABimgn2=0.0_DP
                   dABimgn5=0.0_DP
                   dABimgn6=0.0_DP
                END IF
                !
                ! Precompute damping function (fdamp) and damping function exponential (edamp)...
                ! 
                edamp=EXP(-ddamp*(dABimg/(sR*RAB0)-1.0_DP))
                fdamp=1.0_DP/(1.0_DP+edamp)
                fdamp2=fdamp*fdamp
                !
                ! Apply delta[ia;ib] x delta[n1,n2,n3;0,0,0] conditional...
                !
                IF (n_period.EQ.0.AND.ia.EQ.ib) THEN
                  !
                  ! Do not include self-interaction in the simulation cell...
                  !
                  FDVi=FDVi+0.0_DP
                  FDRi=FDRi+0.0_DP
                  FDRii=FDRii+0.0_DP
                  FCVi=FCVi+0.0_DP
                  FRRi=FRRi+0.0_DP
                  FRRii=FRRii+0.0_DP
                  !
                ELSE
                  !
                  ! Increment image-summed matrix elements...
                  !
                  dptmp1=edamp*fdamp2*dABimgn5
                  FDVi=FDVi+dptmp1
                  !
                  dptmp2=fdamp*dABimgn6
                  FCVi=FCVi+dptmp2
                  !
                  dptmp1=dptmp1*dABimgn2
                  dptmp2=dptmp2*dABimgn2
                  !
                  DO i=1,3
                    !
                    vtmp1(i)=dptmp1*dAB(i)
                    FDRi(i)=FDRi(i)+vtmp1(i)
                    !
                    vtmp2(i)=dptmp2*dAB(i)
                    FRRi(i)=FRRi(i)+vtmp2(i)
                    !
                    DO j=1,3
                      !
                      FDRii(i,j)=FDRii(i,j)+vtmp1(i)*dsAB(j)
                      FRRii(i,j)=FRRii(i,j)+vtmp2(i)*dsAB(j)
                      !
                    END DO
                    !
                  END DO
                  !
                END IF
                !
              END IF !n_period conditional
              !
            END DO !n3
            !
          END DO !n2
          !
        END DO !n1
        !
        ! Increment period energy via EtsvdWAB = - 1/2 * C6ABeff * FAB...
        !
        EtsvdW_period=EtsvdW_period-(FCV0*FRR2*FCVi)
        !
        ! Increment dispersion potential (predveffAdn) prefactor...
        !
        ! predveffAdn(A) = (d * R0freeA * C6ABeff * edamp * fdamp^2) / (6 * sR * vfreeA^1/3 * veffA^2/3 * RAB0^2 * RAB^5)
        !                - (C6ABfree * veffB * fdamp) / (2 * vfreeA * vfreeB * RAB^6)
        !
        ! predveffAdn(B) = (d * R0freeB * C6ABeff * edamp * fdamp^2) / (6 * sR * vfreeB^1/3 * veffB^2/3 * RAB0^2 * RAB^5)
        !                - (C6ABfree * veffA * fdamp) / (2 * vfreeA * vfreeB * RAB^6)
        !
        predveffAdn_period(ia)=predveffAdn_period(ia)-(FDV0*FDV1*FDVA2*FDVi+FCV0*FCVA1*FCVA2*FCVi)
        predveffAdn_period(ib)=predveffAdn_period(ib)-(FDV0*FDVB2*FDVi+FCV0*FCVB1*FCVB2*FCVi)
        !
        ! Increment effective volume derivative contributions to ionic and cell forces...
        !
        ! (dfdamp/dVA) --> D1A = - (d * R0freeA * C6ABeff * edamp * fdamp^2) / (6 * sR * vfreeA^1/3 * veffA^2/3 * RAB0^2 * RAB^5)
        !
        ! (dfdamp/dVB) --> D1B = - (d * R0freeB * C6ABeff * edamp * fdamp^2) / (6 * sR * vfreeB^1/3 * veffB^2/3 * RAB0^2 * RAB^5)
        !
        ! (dC6AB/dVA)  --> D2A =   (C6ABfree * veffB * fdamp) / (2 * vfreeA * vfreeB * RAB^6)
        !
        ! (dC6AB/dVB)  --> D2B =   (C6ABfree * veffA * fdamp) / (2 * vfreeA * vfreeB * RAB^6) 
        !
        D1A=FDV0*FDV1*FDVA2*FDVi  ! (dfdamp/dVA)
        D2A=FCV0*FCVA1*FCVA2*FCVi ! (dC6AB/dVA)
        !
        D1B=FDV0*FDVB2*FDVi       ! (dfdamp/dVB)
        D2B=FCV0*FCVB1*FCVB2*FCVi ! (dC6AB/dVB)
        !
        D12A=D1A+D2A; D12B=D1B+D2B
        !
        DO i=1,3
          !
          DO ic=1,nat
            !
            FtsvdW_period(i,ic)=FtsvdW_period(i,ic)+(dveffdR(ia,ic,i)*D12A+dveffdR(ib,ic,i)*D12B)
            !
          END DO
          !                                            
          DO j=1,3
            !
            HtsvdW_period(i,j)=HtsvdW_period(i,j)+(dveffdh(ia,i,j)*D12A+dveffdh(ib,i,j)*D12B)
            !
          END DO
          !
        END DO
        !
        ! Increment RAB derivative contributions to ionic and cell forces...
        !
        ! (dfdamp/dRA)  --> D1A = (d * C6ABeff * edamp * fdamp^2) / (2 * sR * RAB0 * RAB^7)
        !
        ! (dfdamp/dRB)  --> D1B = - D1A
        !
        ! (dRAB^-6/dRA) --> D2A = - (3 * C6ABeff * fdamp) / (RAB^8)
        !
        ! (dRAB^-6/dRB) --> D2B = - D2A
        !
        D1A=FDR0*FDR2  ! (dfdamp/dRA)
        D2A=FRR0*FRR2  ! (dRAB^-6/dRA)
        !
        ! N.B.: Manually zero out the force contribution from an atom in the simulation cell and
        !       any of its images (this applies to distance derivatives NOT volume derivatives)...
        !
        IF (ia.NE.ib) THEN
          !
          DO i=1,3
            !
            FtsvdW_period(i,ia)=FtsvdW_period(i,ia)+(D1A*FDRi(i)+D2A*FRRi(i))
            FtsvdW_period(i,ib)=FtsvdW_period(i,ib)+(-D1A*FDRi(i)-D2A*FRRi(i))
            !
            DO j=1,3
              !
              HtsvdW_period(i,j)=HtsvdW_period(i,j)+(D1A*FDRii(i,j)+D2A*FRRii(i,j))
              !
            END DO
            !
          END DO
          !
        END IF
        !
      END DO !ib
!$omp end do 
!$omp end parallel
      !
    END DO !iproc
    !
    ! Synchronize n_period contribution from all processors...
    !
    CALL mp_sum(EtsvdW_period,intra_image_comm)
    CALL mp_sum(FtsvdW_period,intra_image_comm)
    CALL mp_sum(HtsvdW_period,intra_image_comm)
    CALL mp_sum(predveffAdn_period,intra_image_comm)
    !
    ! Increment total quantities...
    !
    EtsvdW=EtsvdW+EtsvdW_period                  !EvdW
    FtsvdW=FtsvdW+FtsvdW_period                  !(-dE/dR)
    HtsvdW=HtsvdW-HtsvdW_period                  !(dE/dh)
    predveffAdn=predveffAdn+predveffAdn_period   !(dE/dVA) & (dE/dVB)
    !
    ! DEBUGGING
    !WRITE(stdout,'(I10,2F25.12)') n_period,EtsvdW_period,EtsvdW
    ! DEBUGGING
    !
    ! Periodic convergence loop conditionals...
    !
    IF ( vdw_isolated .OR. (ABS(EtsvdW_period) <= vdw_econv_thr) ) periodic_converged=.TRUE.
    !
    n_period=n_period+1
    !
  END DO !convergence loop
  !
  ! Deallocate temporary arrays...
  !
  IF (ALLOCATED(FtsvdW_period))      DEALLOCATE(FtsvdW_period)
  IF (ALLOCATED(HtsvdW_period))      DEALLOCATE(HtsvdW_period)
  IF (ALLOCATED(predveffAdn_period)) DEALLOCATE(predveffAdn_period)
  !
  CALL stop_clock('tsvdw_energy')
  !
 !! DEBUGGING
 !WRITE(stdout,'(3X,"<START veff ARRAY>")') 
 !DO ia=1,nat
 !  WRITE(stdout,'(5X,I5,F25.12)') ia,veff(ia)
 !END DO
 !WRITE(stdout,'(3X,"<END veff ARRAY>")') 
 !!
 !WRITE(stdout,'(3X,"<START FtsvdW ARRAY>")') 
 !DO ia=1,nat
 !  WRITE(stdout,'(5X,I5,3F25.12)') ia,FtsvdW(:,ia)
 !END DO
 !WRITE(stdout,'(3X,"<END FtsvdW ARRAY>")') 
 !!
 !WRITE(stdout,'(3X,"<START HtsvdW ARRAY>")')
 !DO i=1,3
 !  WRITE(stdout,'(5X,3F25.12)') HtsvdW(i,:)
 !END DO
 !WRITE(stdout,'(3X,"<END HtsvdW ARRAY>")')
 !! DEBUGGING
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_energy
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_wfforce()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ia,ip,iq,off1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: UtsvdWA
  !
  CALL start_clock('tsvdw_wfforce')
  !
  ! Initialization of UtsvdwA array...
  !
  ALLOCATE(UtsvdWA(nr1*nr2*nr3)); UtsvdWA=0.0_DP
  !
  ! Loop over atoms and populate UtsvdWA from predveffAdn and dveffAdn...
  !
  DO iproc=1,nstates(me)
    !
    ! Connect processor number with atom...
    !
    ia=me+nproc_image*(iproc-1)
    !
    ! Loop over points in the (pre-screened) spherical atomic integration domain...
    !
!$omp parallel do private(off1)
    DO iq=1,NsomegaA(ia)
      !
      off1=somegaA(iq,1,iproc)+(somegaA(iq,2,iproc)-1)*nr1+(somegaA(iq,3,iproc)-1)*nr1*nr2    !global offset [nr1,nr2,nr3]
      UtsvdWA(off1)=UtsvdWA(off1)+predveffAdn(ia)*dveffAdn(iq,iproc)
      !
    END DO !iq
!$omp end parallel do
    !
  END DO !iproc
  !
  ! Collect UtsvdWA over all processors and broadcast...
  !
  CALL mp_sum(UtsvdWA,intra_image_comm)
  !
  ! Partition out dispersion potential consistent with slabs of the charge density...
  !
  IF (dfftp%npp(me_bgrp+1).NE.0) THEN
    !
!$omp parallel do 
    DO ip=1,dfftp%npp(me_bgrp+1)*nr1*nr2
      !
      UtsvdW(ip)=UtsvdWA(ip+rdispls(me_bgrp+1)) 
      !
    END DO
!$omp end parallel do 
    !
  END IF
  !
  ! Deallocate temporary arrays...
  !
  IF (ALLOCATED(UtsvdWA))      DEALLOCATE(UtsvdWA)
  !
  CALL stop_clock('tsvdw_wfforce')
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_wfforce
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_cleanup()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Deallocate tsvdw_calculate specific arrays...
  !
  IF (ALLOCATED(atxyz))       DEALLOCATE(atxyz)
  IF (ALLOCATED(rhosad))      DEALLOCATE(rhosad)
  IF (ALLOCATED(rhotot))      DEALLOCATE(rhotot)
  IF (ALLOCATED(veff))        DEALLOCATE(veff)
  IF (ALLOCATED(dpeff))       DEALLOCATE(dpeff)
  IF (ALLOCATED(R0eff))       DEALLOCATE(R0eff)
  IF (ALLOCATED(C6AAeff))     DEALLOCATE(C6AAeff)
  IF (ALLOCATED(C6ABeff))     DEALLOCATE(C6ABeff)
  IF (ALLOCATED(dveffdR))     DEALLOCATE(dveffdR)
  IF (ALLOCATED(dveffdh))     DEALLOCATE(dveffdh)
  IF (ALLOCATED(somegaA))     DEALLOCATE(somegaA)
  IF (ALLOCATED(somegaAr))    DEALLOCATE(somegaAr)
  IF (ALLOCATED(gomegaAr))    DEALLOCATE(gomegaAr)
  IF (ALLOCATED(NsomegaA))    DEALLOCATE(NsomegaA)
  IF (ALLOCATED(NsomegaAr))   DEALLOCATE(NsomegaAr)
  IF (ALLOCATED(nstates))     DEALLOCATE(nstates)
  IF (ALLOCATED(sdispls))     DEALLOCATE(sdispls)
  IF (ALLOCATED(sendcount))   DEALLOCATE(sendcount)
  IF (ALLOCATED(rdispls))     DEALLOCATE(rdispls)
  IF (ALLOCATED(recvcount))   DEALLOCATE(recvcount)
  IF (ALLOCATED(istatus))     DEALLOCATE(istatus)
  IF (ALLOCATED(npair))       DEALLOCATE(npair)
  IF (ALLOCATED(pair))        DEALLOCATE(pair)
  IF (ALLOCATED(dveffAdn))    DEALLOCATE(dveffAdn)
  IF (ALLOCATED(predveffAdn)) DEALLOCATE(predveffAdn)
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_cleanup
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE tsvdw_finalize()
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! Deallocate module-specific arrays...
  !
  IF (ALLOCATED(UtsvdW))   DEALLOCATE(UtsvdW)
  IF (ALLOCATED(FtsvdW))   DEALLOCATE(FtsvdW)
  IF (ALLOCATED(HtsvdW))   DEALLOCATE(HtsvdW)
  IF (ALLOCATED(VefftsvdW))DEALLOCATE(VefftsvdW)
  IF (ALLOCATED(vfree))    DEALLOCATE(vfree)
  IF (ALLOCATED(dpfree))   DEALLOCATE(dpfree)
  IF (ALLOCATED(R0free))   DEALLOCATE(R0free)
  IF (ALLOCATED(C6AAfree)) DEALLOCATE(C6AAfree)
  IF (ALLOCATED(C6ABfree)) DEALLOCATE(C6ABfree)
  IF (ALLOCATED(spgrd))    DEALLOCATE(spgrd)
  IF (ALLOCATED(sprho))    DEALLOCATE(sprho)
  IF (ALLOCATED(spdrho))   DEALLOCATE(spdrho)
  IF (ALLOCATED(spdata))   DEALLOCATE(spdata)
  IF (ALLOCATED(LIA))      DEALLOCATE(LIA)
  IF (ALLOCATED(LIB))      DEALLOCATE(LIB)
  IF (ALLOCATED(dLIA))     DEALLOCATE(dLIA)
  IF (ALLOCATED(dLIB))     DEALLOCATE(dLIB)
  !
  RETURN
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE tsvdw_finalize
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE Num1stDer(r,f,N,h,df)
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER :: N
  REAL(DP) :: h,r(N),f(N),df(N)
  !
  ! Local variables
  !
  INTEGER, PARAMETER :: Ndp=7 !using 7-point formulae...
  INTEGER :: ip,ir
  INTEGER :: A1(Ndp),A2(Ndp),A3(Ndp),A4(Ndp),A5(Ndp),A6(Ndp),A7(Ndp)
  REAL(DP) :: dsum
  !
  ! Populate Bickley coefficient vectors (Math. Gaz.,v25,p19-27,1941) according to cases...
  !
  DATA  A1/ -1764_DP, 4320_DP, -5400_DP,  4800_DP, -2700_DP,   864_DP, -120_DP /
  DATA  A2/  -120_DP, -924_DP,  1800_DP, -1200_DP,   600_DP,  -180_DP,   24_DP /
  DATA  A3/    24_DP, -288_DP,  -420_DP,   960_DP,  -360_DP,    96_DP,  -12_DP /
  DATA  A4/   -12_DP,  108_DP,  -540_DP,     0_DP,   540_DP,  -108_DP,   12_DP /
  DATA  A5/    12_DP,  -96_DP,   360_DP,  -960_DP,   420_DP,   288_DP,  -24_DP /
  DATA  A6/   -24_DP,  180_DP,  -600_DP,  1200_DP, -1800_DP,   924_DP,  120_DP /
  DATA  A7/   120_DP, -864_DP,  2700_DP, -4800_DP,  5400_DP, -4320_DP, 1764_DP /
  !
  ! Compute first derivative on linear mesh and then transform back to radial/exponential grid...
  !
  DO ir=1,N
    !
    dsum=0.0_DP
    !
    ! Deal with different cases one-by-one...
    !
    IF (ir.EQ.1) THEN
      DO ip=1,Ndp
        dsum=dsum+A1(ip)*f(ir-1+ip)     
      END DO
    ELSE IF (ir.EQ.2) THEN
      DO ip=1,Ndp
        dsum=dsum+A2(ip)*f(ir-2+ip)     
      END DO
    ELSE IF (ir.EQ.3) THEN
      DO ip=1,Ndp
        dsum=dsum+A3(ip)*f(ir-3+ip)     
      END DO
    ELSE IF (ir.GE.4.AND.ir.LE.N-3) THEN
      DO ip=1,Ndp
        dsum=dsum+A4(ip)*f(ir-4+ip)     
      END DO
    ELSE IF (ir.EQ.N-2) THEN
      DO ip=1,Ndp
        dsum=dsum+A5(ip)*f(ir-5+ip)     
      END DO
    ELSE IF (ir.EQ.N-1) THEN
      DO ip=1,Ndp
        dsum=dsum+A6(ip)*f(ir-6+ip)     
      END DO
    ELSE IF (ir.EQ.N) THEN
      DO ip=1,Ndp
        dsum=dsum+A7(ip)*f(ir-7+ip)     
      END DO
    ELSE
      WRITE(stdout,'("Error in Num1stDer subroutine...")')
    END IF
    !
    ! Final Normalization...
    !
    df(ir)=dsum/(720.0_DP*h)
    !
  END DO !ir
  ! 
  RETURN 
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE Num1stDer
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE CubSplCoeff(r,f,N,df,d2f)
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER :: N
  REAL(DP) :: r(N),f(N),df(N),d2f(N)
  !
  ! Local variables
  !
  INTEGER :: i,j
  REAL(DP) :: dy1,dyn,p,q,un,qn
  REAL(DP), DIMENSION(:), ALLOCATABLE :: work
  !
  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! SYNOPSIS: Compute second derivatives at each of the atomic radial grid points using the cubic spline methodology (i.e., smooth &
  ! continuous piecewise first and second derivatives). These second derivatives will be utilized during cubic spline interpolation 
  ! as a higher accuracy alternative to linear interpolation during the construction of the linear atomic grids. The two-parameter
  ! boundary conditions that will be utilized below are known as a clamped cubic spline in that the first derivative at both the 
  ! first and last grid point were computed numerically and provided as input...
  ! ----------------------------------------------------------------------------------------------------------------------------------
  !
  ALLOCATE(work(N)); work=0.0_DP
  !
  d2f=0.0_DP
  !
  ! Enforce 'clamped' boundary condition at the first radial grid point...
  !
  dy1=df(1)
  d2f(1)=-0.5_DP
  work(1)=(3.0_DP/(r(2)-r(1)))*((f(2)-f(1))/(r(2)-r(1))-dy1)
  !
  ! Decomposition loop of the tridiagonal algorithm for the second derivatives...
  !
  DO i=2,N-1
    p=(r(i)-r(i-1))/(r(i+1)-r(i-1))
    q=p*d2f(i-1)+2.0_DP
    d2f(i)=(p-1.0_DP)/q
    work(i)=(f(i+1)-f(i))/(r(i+1)-r(i))-(f(i)-f(i-1))/(r(i)-r(i-1))
    work(i)=(6.0_DP*work(i)/(r(i+1)-r(i-1))-p*work(i-1))/q
  END DO
  !
  ! Enforce 'clamped' boundary condition at the last radial grid point...
  !
  dyn=df(N)
  qn=0.5_DP
  un=(3.0_DP/(r(N)-r(N-1)))*(dyn-(f(N)-f(N-1))/(r(N)-r(N-1)))
  d2f(N)=(un-qn*work(N-1))/(qn*d2f(N-1)+1.0_DP)
  !
  ! Back substitution loop of the tridiagonal algorithm for the second derivatives...
  !
  DO j=N-1,1,-1
    d2f(j)=d2f(j)*d2f(j+1)+work(j)
  END DO
  !
  ! Clean-up and return home...
  ! 
  DEALLOCATE(work)
  ! 
  RETURN 
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE CubSplCoeff
  !--------------------------------------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE GetVdWParam(atom,C6,alpha,R0)
  !--------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  CHARACTER(LEN=3) :: atom
  REAL(DP) :: C6,alpha,R0
  !
  SELECT CASE (atom) 
  !
  CASE ('H')
  alpha=4.500000_DP
  C6=6.500000_DP
  R0=3.100000_DP
  !
  CASE ('He')
  alpha=1.380000_DP
  C6=1.460000_DP
  R0=2.650000_DP
  !
  CASE ('Li')
  alpha=164.200000_DP
  C6=1387.000000_DP
  R0=4.160000_DP
  !
  CASE ('Be')
  alpha=38.000000_DP
  C6=214.000000_DP
  R0=4.170000_DP
  !
  CASE ('B')
  alpha=21.000000_DP
  C6=99.500000_DP
  R0=3.890000_DP
  !
  CASE ('C')
  alpha=12.000000_DP
  C6=46.600000_DP
  R0=3.590000_DP
  !
  CASE ('N')
  alpha=7.400000_DP
  C6=24.200000_DP
  R0=3.340000_DP
  !
  CASE ('O')
  alpha=5.400000_DP
  C6=15.600000_DP
  R0=3.190000_DP
  !
  CASE ('F')
  alpha=3.800000_DP
  C6=9.520000_DP
  R0=3.040000_DP
  !
  CASE ('Ne')
  alpha=2.670000_DP
  C6=6.380000_DP
  R0=2.910000_DP
  !
  CASE ('Na')
  alpha=162.700000_DP
  C6=1556.000000_DP
  R0=3.730000_DP
  !
  CASE ('Mg')
  alpha=71.000000_DP
  C6=627.000000_DP
  R0=4.270000_DP
  !
  CASE ('Al')
  alpha=60.000000_DP
  C6=528.000000_DP
  R0=4.330000_DP
  !
  CASE ('Si')
  alpha=37.000000_DP
  C6=305.000000_DP
  R0=4.200000_DP
  !
  CASE ('P')
  alpha=25.000000_DP
  C6=185.000000_DP
  R0=4.010000_DP
  !
  CASE ('S')
  alpha=19.600000_DP
  C6=134.000000_DP
  R0=3.860000_DP
  !
  CASE ('Cl')
  alpha=15.000000_DP
  C6=94.600000_DP
  R0=3.710000_DP
  !
  CASE ('Ar')
  alpha=11.100000_DP
  C6=64.300000_DP
  R0=3.550000_DP
  !
  CASE ('K')
  alpha=292.900000_DP
  C6=3897.000000_DP
  R0=3.710000_DP
  !
  CASE ('Ca')
  alpha=160.000000_DP
  C6=2221.000000_DP
  R0=4.650000_DP
  !
  CASE ('Sc')
  alpha=120.000000_DP
  C6=1383.000000_DP
  R0=4.590000_DP
  !
  CASE ('Ti')
  alpha=98.000000_DP
  C6=1044.000000_DP
  R0=4.510000_DP
  !
  CASE ('V')
  alpha=84.000000_DP
  C6=832.000000_DP
  R0=4.440000_DP
  !
  CASE ('Cr')
  alpha=78.000000_DP
  C6=602.000000_DP
  R0=3.990000_DP
  !
  CASE ('Mn')
  alpha=63.000000_DP
  C6=552.000000_DP
  R0=3.970000_DP
  !
  CASE ('Fe')
  alpha=56.000000_DP
  C6=482.000000_DP
  R0=4.230000_DP
  !
  CASE ('Co')
  alpha=50.000000_DP
  C6=408.000000_DP
  R0=4.180000_DP
  !
  CASE ('Ni')
  alpha=48.000000_DP
  C6=373.000000_DP
  R0=3.820000_DP
  !
  CASE ('Cu')
  alpha=42.000000_DP
  C6=253.000000_DP
  R0=3.760000_DP
  !
  CASE ('Zn')
  alpha=40.000000_DP
  C6=284.000000_DP
  R0=4.020000_DP
  !
  CASE ('Ga')
  alpha=60.000000_DP
  C6=498.000000_DP
  R0=4.190000_DP
  !
  CASE ('Ge')
  alpha=41.000000_DP
  C6=354.000000_DP
  R0=4.200000_DP
  !
  CASE ('As')
  alpha=29.000000_DP
  C6=246.000000_DP
  R0=4.110000_DP
  !
  CASE ('Se')
  alpha=25.000000_DP
  C6=210.000000_DP
  R0=4.040000_DP
  !
  CASE ('Br')
  alpha=20.000000_DP
  C6=162.000000_DP
  R0=3.930000_DP
  !
  CASE ('Kr')
  alpha=16.800000_DP
  C6=129.600000_DP
  R0=3.820000_DP
  !
  CASE ('Rb')
  alpha=319.200000_DP
  C6=4691.000000_DP
  R0=3.720000_DP
  !
  CASE ('Sr')
  alpha=199.000000_DP
  C6=3170.000000_DP
  R0=4.540000_DP
  !
  CASE ('Y')
  alpha=126.7370_DP
  C6=1968.580_DP
  R0=4.81510_DP
  !
  CASE ('Zr')
  alpha=119.97_DP
  C6=1677.91_DP
  R0=4.53_DP
  !
  CASE ('Nb')
  alpha=101.603_DP
  C6=1263.61_DP
  R0=4.2365_DP
  !
  CASE ('Mo')
  alpha=88.4225785_DP
  C6=1028.73_DP
  R0=4.099_DP
  !
  CASE ('Tc')
  alpha=80.083_DP
  C6=1390.87_DP
  R0=4.076_DP
  !
  CASE ('Ru')
  alpha=65.8950_DP
  C6=609.754_DP
  R0=3.99530_DP
  !
  CASE ('Rh')
  alpha=56.1_DP
  C6=469.0_DP
  R0=3.95_DP
  !
  CASE ('Pd')
  alpha=23.680000_DP
  C6=157.500000_DP
  R0=3.66000_DP
  !
  CASE ('Ag')
  alpha=50.600000_DP
  C6=339.000000_DP
  R0=3.820000_DP
  !
  CASE ('Cd')
  alpha=39.7_DP
  C6=452.0_DP
  R0=3.99_DP
  !
  CASE ('In')
  alpha=70.22000_DP
  C6=707.046000_DP
  R0=4.23198000_DP
  !
  CASE ('Sn')
  alpha=55.9500_DP
  C6=587.41700_DP
  R0=4.303000_DP
  !
  CASE ('Sb')
  alpha=43.671970_DP
  C6=459.322_DP
  R0=4.2760_DP
  !
  CASE ('Te')
  alpha=37.65_DP
  C6=396.0_DP
  R0=4.22_DP
  !
  CASE ('I')
  alpha=35.000000_DP
  C6=385.000000_DP
  R0=4.170000_DP
  !
  CASE ('Xe')
  alpha=27.300000_DP
  C6=285.900000_DP
  R0=4.080000_DP
  !
  CASE ('Cs')
  alpha=427.12_DP
  C6=6582.08_DP
  R0=3.78_DP
  !
  CASE ('Ba')
  alpha=275.0_DP
  C6=5727.0_DP
  R0=4.77_DP
  !
  CASE ('Hf')
  alpha=99.52_DP
  C6=1274.8_DP
  R0=4.21_DP
  !
  CASE ('Ta')
  alpha=82.53_DP
  C6=1019.92_DP
  R0=4.15_DP
  !
  CASE ('W')
  alpha=71.041_DP
  C6=847.93_DP
  R0=4.08_DP
  !
  CASE ('Re')
  alpha=63.04_DP
  C6=710.2_DP
  R0=4.02_DP
  !
  CASE ('Os')
  alpha=55.055_DP
  C6=596.67_DP
  R0=3.84_DP
  !
  CASE ('Ir')
  alpha=42.51_DP
  C6=359.1_DP
  R0=4.00_DP
  !
  CASE ('Pt')
  alpha=39.68_DP
  C6=347.1_DP
  R0=3.92_DP
  !
  CASE ('Au')
  alpha=36.5_DP
  C6=298.0_DP
  R0=3.86_DP
  !
  CASE ('Hg')
  alpha=33.9_DP
  C6=392.0_DP
  R0=3.98_DP
  !
  CASE ('Tl')
  alpha=69.92_DP
  C6=717.44_DP
  R0=3.91_DP
  !
  CASE ('Pb')
  alpha=61.8_DP
  C6=697.0_DP
  R0=4.31_DP
  !
  CASE ('Bi')
  alpha=49.02_DP
  C6=571.0_DP
  R0=4.32_DP
  !
  CASE ('Po')
  alpha=45.013_DP
  C6=530.92_DP
  R0=4.097_DP
  !
  CASE ('At')
  alpha=38.93_DP
  C6=457.53_DP
  R0=4.07_DP
  !
  CASE ('Rn')
  alpha=33.54_DP
  C6=390.63_DP
  R0=4.23_DP
  !
  CASE DEFAULT
  !
  CALL errore('tsvdw','Reference free atom parameters not available for requested atom type...',1)
  !
  END SELECT
  !
  RETURN 
  !
  !--------------------------------------------------------------------------------------------------------------
  END SUBROUTINE GetVdWParam
  !--------------------------------------------------------------------------------------------------------------
  !
  !
END MODULE tsvdw_module

