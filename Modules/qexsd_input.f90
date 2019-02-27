!
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------
MODULE qexsd_input
!--------------------------------------------------------
  ! This module contains the routines needed to initialise the data-structures 
  ! contained in the PW XML input
  !----------------------------------------------------------------------------
  ! First version March 2016 
  !----------- ------------- --------------------------------------------------- 
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : input_xml_schema_file
  !
  USE constants,        ONLY : e2,bohr_radius_angs
  USE qes_types_module
  USE qes_libs_module
  !
  IMPLICIT NONE
  !
  
  PUBLIC                 
  SAVE
  !
  TYPE(input_type)               ::     input 
  CONTAINS
  !--------------------------------------------------------------------------------------------------------------------  
  SUBROUTINE  qexsd_init_control_variables(obj,title,calculation,restart_mode,&
                  prefix,pseudo_dir,outdir,stress,forces,wf_collect,disk_io,  &
                  max_seconds,etot_conv_thr,forc_conv_thr,press_conv_thr,verbosity, &
                  iprint, nstep) 
  !---------------------------------------------------------------------------------------------------------------------
  !
  TYPE(control_variables_type)         :: obj
  CHARACTER(LEN=*),INTENT(IN)          :: title,calculation,restart_mode,prefix,&
                                          pseudo_dir,outdir,disk_io,verbosity
  LOGICAL,INTENT(IN)                   :: stress,forces,wf_collect
  REAL(DP),INTENT(IN)                  :: max_seconds,etot_conv_thr,forc_conv_thr,&
                                          press_conv_thr   
  INTEGER,INTENT(IN)                   :: iprint, nstep
  OPTIONAL                             :: nstep
  !
  !
  CHARACTER(LEN=*),PARAMETER           :: TAGNAME='control_variables'
  CHARACTER(LEN=256)                   :: verbosity_value, disk_io_value
  INTEGER                              :: int_max_seconds
  LOGICAL                              :: nstep_ispresent

  int_max_seconds=nint(max_seconds)
  IF ( TRIM( verbosity ) .EQ. 'default' ) THEN 
        verbosity_value = "low"
  ELSE
        verbosity_value=TRIM(verbosity) 
  END IF
  IF ( TRIM(disk_io) .EQ. 'default' ) THEN 
     disk_io_value="low"
  ELSE
     disk_io_value=TRIM(disk_io)
  END IF
  !
  !
  CALL qes_init (obj,tagname,title=TRIM(title),calculation=TRIM(calculation),&
                                  restart_mode=TRIM(restart_mode),prefix=TRIM(prefix),        &
                                  pseudo_dir=TRIM(pseudo_dir),outdir=TRIM(outdir),disk_io=TRIM(disk_io_value),&
                                  verbosity=TRIM(verbosity_value),stress=stress,forces=forces,    &
                                  wf_collect=wf_collect,max_seconds=int_max_seconds,  &
                                  etot_conv_thr=etot_conv_thr,forc_conv_thr=forc_conv_thr, &
                                  press_conv_thr=press_conv_thr,print_every=iprint, NSTEP = nstep )

  END SUBROUTINE qexsd_init_control_variables
  !
  !
  !----------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_spin(obj,lsda,noncolin,spinorbit) 
  ! 
  IMPLICIT NONE
  ! 
  TYPE(spin_type)                 :: obj
  LOGICAL,INTENT(IN)              :: lsda,noncolin,spinorbit
  !
  CHARACTER(LEN=*),PARAMETER      :: TAGNAME="spin"
  
  CALL qes_init (obj,TAGNAME,lsda=lsda,noncolin=noncolin,spinorbit=spinorbit)
  
  END SUBROUTINE qexsd_init_spin  
  !
  !
  !-------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_bands(obj, nbnd, smearing, degauss, occupations, tot_charge, nspin, & 
                              input_occupations, input_occupations_minority, tot_mag)
  !
  IMPLICIT NONE
  ! 
  TYPE ( bands_type)                           :: obj
  INTEGER,OPTIONAL, INTENT(IN)                 :: nbnd 
  INTEGER,INTENT(IN)                           :: nspin
  CHARACTER(LEN=*),INTENT(IN)                  :: occupations,smearing
  REAL(DP),INTENT(IN)                          :: degauss 
  REAL(DP),DIMENSION(:),OPTIONAL,INTENT(IN)    :: input_occupations, input_occupations_minority
  REAL(DP),OPTIONAL,INTENT(IN)                 :: tot_mag, tot_charge 
  !
  CHARACTER(25)                                :: smearing_local
  INTEGER                                      :: spin_degeneracy, inpOcc_size = 0
  CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="bands"
  TYPE(smearing_type),POINTER                  :: smearing_obj => NULL()
  TYPE(occupations_type)                       :: occup_obj
  TYPE(inputoccupations_type),ALLOCATABLE      :: inpOcc_objs(:)
  LOGICAL                                      :: tot_mag_ispresent = .FALSE., &
                                                  inp_occ_arepresent = .FALSE.
  ! 
  IF (TRIM(occupations) .EQ. "smearing")  THEN
     ALLOCATE(smearing_obj) 
     SELECT CASE (TRIM  (smearing))
       CASE ("gaussian", "gauss")
           smearing_local="gaussian"
       CASE ('methfessel-paxton', 'm-p', 'mp')
           smearing_local="mp"
       CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv','Marzari-Vanderbilt') 
           smearing_local="mv"
       CASE ('fermi-dirac', 'f-d', 'fd') 
           smearing_local="fd"
     END SELECT 
     CALL qes_init (smearing_obj,"smearing",degauss=degauss,smearing=smearing_local)
  END IF
  CALL  qes_init (occup_obj, "occupations", occupations = TRIM(occupations))
  !
  IF (PRESENT(input_occupations) ) THEN 
     SELECT CASE ( nspin)
       CASE (2) 
          inpOcc_size=2
       CASE default
          inpOcc_size=1
     END SELECT
     ALLOCATE (inpOcc_objs(inpOcc_size))
     IF ( inpOcc_size .GT. 1) THEN 
        CALL qes_init ( inpOcc_objs(1),"input_occupations", ISPIN = 1, &
                  SPIN_FACTOR = 1._DP, INPUTOCCUPATIONS = input_occupations(2:nbnd) ) 
        CALL qes_init ( inpOcc_objs(2),"input_occupations", 2, & 
                  SPIN_FACTOR = 1._DP , INPUTOCCUPATIONS = input_occupations_minority(2:nbnd))
     ELSE 
        CALL qes_init ( inpOcc_objs(1),"input_occupations", ISPIN = 1, SPIN_FACTOR = 2._DP , &
                                                                 INPUTOCCUPATIONS = input_occupations(2:nbnd) )   
     END IF
  END IF
  ! 
  CALL qes_init (obj, TAGNAME, NBND = nbnd, SMEARING = smearing_obj, TOT_CHARGE = tot_charge, &
                      TOT_MAGNETIZATION = tot_mag, OCCUPATIONS=occup_obj, INPUTOCCUPATIONS = inpOcc_objs )
  IF (ASSOCIATED(smearing_obj)) THEN 
      CALL qes_reset (smearing_obj)
      DEALLOCATE ( smearing_obj) 
  END IF 
  CALL qes_reset (occup_obj)
  IF (ALLOCATED(inpOcc_objs)) THEN 
     CALL qes_reset (inpocc_objs(1))
     IF (inpOcc_size .GT. 1 ) CALL qes_reset (inpocc_objs(2))
     DEALLOCATE (inpocc_objs)
  END IF
  !
  END SUBROUTINE qexsd_init_bands
  !
  !
  !--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_basis(obj,k_points,ecutwfc,ecutrho,nr,nrs,nrb)
  !--------------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE (basis_type)                 :: obj
  CHARACTER(LEN=*),INTENT(IN)       :: k_points
  REAL(DP),INTENT(IN)               :: ecutwfc 
  REAL(DP),OPTIONAL,INTENT(IN)      :: ecutrho
  INTEGER,OPTIONAL,INTENT(IN)       :: nr(:), nrs(:), nrb(:) 
  ! 
  TYPE(basisSetItem_type),POINTER   :: grid_obj => NULL(), smooth_grid_obj => NULL(), box_obj => NULL()
  CHARACTER(LEN=*),PARAMETER        :: TAGNAME="basis",FFT_GRID="fft_grid",FFT_SMOOTH="fft_smooth", FFT_BOX="fft_box"
  LOGICAL                           :: gamma_only=.FALSE. 
  !
  IF ( PRESENT(nr)) THEN
    ALLOCATE(grid_obj) 
    CALL qes_init (grid_obj,FFT_GRID,nr(1),nr(2),nr(3),"grid set in input")
  END IF
  ! 
  IF( PRESENT(nrs)) THEN
    ALLOCATE(smooth_grid_obj)
    CALL qes_init (smooth_grid_obj,FFT_SMOOTH,nrs(1),nrs(2),nrs(3),"grid set in input")
  END IF
  ! 
  IF( PRESENT(nrb)) THEN
    ALLOCATE(box_obj) 
    CALL qes_init (box_obj,FFT_BOX,nrb(1),nrb(2),nrb(3),"grid set in input")
  END IF
  ! 
  IF (TRIM(k_points) .EQ. "gamma" ) gamma_only=.TRUE.

  CALL qes_init (obj,TAGNAME, GAMMA_ONLY=gamma_only,ECUTWFC=ecutwfc, ECUTRHO=ecutrho, FFT_GRID=grid_obj, &
                                                                    FFT_SMOOTH=smooth_grid_obj, FFT_BOX=box_obj)
  !
  IF (ASSOCIATED(grid_obj))        CALL  qes_reset( grid_obj )
  IF (ASSOCIATED(smooth_grid_obj)) CALL  qes_reset( smooth_grid_obj )
  IF (ASSOCIATED(box_obj))         CALL  qes_reset( box_obj )
  ! 
  !
  !
  END SUBROUTINE qexsd_init_basis
  !-------------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_electron_control( obj,diagonalization,mixing_mode,mixing_beta,&
                                          conv_thr, mixing_ndim, max_nstep, tqr, real_space, &
                                          tq_smoothing, tbeta_smoothing, & 
                                          diago_thr_init, diago_full_acc, &
                                          diago_cg_maxiter, diago_ppcg_maxiter, diago_david_ndim)
  !-------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE(electron_control_type)             ::  obj
  CHARACTER(LEN=*),INTENT(IN)             :: diagonalization,mixing_mode
  REAL(DP),INTENT(IN)                     :: mixing_beta, conv_thr, diago_thr_init
  INTEGER,INTENT(IN)                      :: mixing_ndim,max_nstep, diago_cg_maxiter, &
                                             diago_ppcg_maxiter, diago_david_ndim
  LOGICAL,OPTIONAL,INTENT(IN)             :: diago_full_acc,tqr, real_space, tq_smoothing, tbeta_smoothing
  !
  CHARACTER(LEN=*),PARAMETER              :: TAGNAME="electron_control"
  !
  CALL qes_init (obj,TAGNAME, DIAGONALIZATION=diagonalization,&
                                MIXING_MODE=mixing_mode,MIXING_BETA=mixing_beta,&
                                CONV_THR=conv_thr,MIXING_NDIM=mixing_ndim,MAX_NSTEP=max_nstep,&
                                TQ_SMOOTHING= tq_smoothing, TBETA_SMOOTHING = tbeta_smoothing,& 
                                REAL_SPACE_Q=tqr, REAL_SPACE_BETA = real_space, DIAGO_THR_INIT=diago_thr_init,& 
                                DIAGO_FULL_ACC=diago_full_acc,DIAGO_CG_MAXITER=diago_cg_maxiter, &
                                DIAGO_PPCG_MAXITER=diago_ppcg_maxiter)
   !
   END SUBROUTINE qexsd_init_electron_control
   !
   !
   !-------------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_k_points_ibz(obj,k_points,calculation,nk1,nk2,nk3,s1,s2,s3,nk,xk,wk,alat,a1, ibrav_lattice)
   ! 
   IMPLICIT NONE
   ! 
   TYPE (k_points_IBZ_type)             :: obj
   CHARACTER(LEN=*),INTENT(IN)          :: k_points,calculation
   INTEGER,INTENT(IN)                   :: nk1,nk2,nk3,s1,s2,s3,nk
   REAL(DP),INTENT(IN)                  :: xk(:,:),wk(:)
   REAL(DP),INTENT(IN)                  :: alat,a1(3)
   LOGICAL,INTENT(IN)                   :: ibrav_lattice
   !
   CHARACTER(LEN=*),PARAMETER           :: TAGNAME="k_points_IBZ"
   TYPE(monkhorst_pack_type),POINTER    :: mpack_obj_pt => NULL() 
   TYPE(monkhorst_pack_type),TARGET     :: mpack_obj_ 
   TYPE(k_point_type),ALLOCATABLE       :: kp_obj(:)
   LOGICAL                              :: mpack_ispresent,kp_ispresent
   CHARACTER(LEN=100)                   :: kind_of_grid
   INTEGER                              :: ik,jk,kcount
   REAL(DP),DIMENSION(3)                :: my_xk
   REAL(DP)                             :: scale_factor
   INTEGER, POINTER                     :: kdim_opt => NULL()
   INTEGER, TARGET                      :: kdim 
   !
  
   IF (TRIM(k_points).EQ."automatic") THEN 
      !
      IF ((s1+s2+s3).EQ.0) THEN
         kind_of_grid="Monkhorst-Pack"
      ELSE
         kind_of_grid="Uniform grid with offset"
      END IF
      CALL qes_init (mpack_obj_,"monkhorst_pack",nk1,nk2,nk3, s1,s2,s3,kind_of_grid)
      mpack_obj_pt => mpack_obj_
   ELSE
      kdim_opt => kdim 
      IF ( ibrav_lattice ) THEN 
         scale_factor = 1.d0
      ELSE 
         scale_factor=alat/sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
      END IF 
      !
      IF (TRIM(calculation).NE.'bands' .AND. (TRIM(k_points).EQ.'tpiba_b' .OR. &
                                              TRIM(k_points) .EQ. 'crystal_b')) THEN
          kdim=NINT(sum(wk(1:nk-1)))+1
          ALLOCATE (kp_obj(kdim))
          kcount=1
          CALL qes_init (kp_obj(kcount),"k_point", WEIGHT = 1.d0, K_POINT = xk(:,1))
          kcount=kcount+1
          DO ik=1,nk-1
             DO jk=1,NINT(wk(ik))
                my_xk=xk(:,ik)+(DBLE(jk)/wk(ik))*(xk(:,ik+1)-xk(:,ik))
                my_xk=my_xk*scale_factor
                CALL qes_init (kp_obj(kcount),"k_point",WEIGHT = 1.d0, K_POINT = my_xk)
                kcount=kcount+1
             END DO
          END DO
      ELSE
          kdim=nk
          ALLOCATE  (kp_obj(kdim))      
          DO ik=1,kdim
             my_xk=xk(:,ik)*scale_factor
             CALL qes_init (kp_obj(ik),"k_point", WEIGHT = wk(ik),K_POINT=my_xk)
          END DO
      END IF
   END IF    
   CALL qes_init (obj, TAGNAME, MONKHORST_PACK = mpack_obj_pt, NK = kdim_opt , K_POINT = kp_obj) 
   IF (ASSOCIATED (mpack_obj_pt)) THEN 
      CALL qes_reset (mpack_obj_)
      mpack_obj_pt => NULL() 
   ELSE  IF (ALLOCATED(kp_obj)) THEN 
      DO ik = 1, kdim 
         CALL qes_reset(kp_obj(ik))
      END DO 
      DEALLOCATE (kp_obj) ! this line is redundant because kp_obj is a local allocatable   
   END IF 
   
   END SUBROUTINE qexsd_init_k_points_ibz
   !
   ! 
   !--------------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_ion_control(obj,ion_dynamics,upscale,remove_rigid_rot,&
                                     refold_pos,pot_extrapolation,wfc_extrapolation,&
                                      ion_temperature,tempw,tolp,delta_t,nraise,dt,&
                                      bfgs_ndim,trust_radius_min,trust_radius_max,&
                                      trust_radius_init,w_1,w_2)
   !--------------------------------------------------------------------------------------------------
   !
   IMPLICIT NONE
   ! 
   TYPE (ion_control_type)                 :: obj
   CHARACTER(LEN=*),INTENT(IN)             :: ion_dynamics,pot_extrapolation,wfc_extrapolation,&
                                              ion_temperature
   REAL(DP),OPTIONAL,INTENT(IN)            :: upscale, tempw,tolp,delta_t,trust_radius_min,trust_radius_max,&
                                              trust_radius_init,w_1,w_2
   INTEGER,INTENT(IN)                      :: nraise,bfgs_ndim
   REAL(DP),INTENT(IN)                     :: dt
   LOGICAL,OPTIONAL,INTENT(IN)             :: remove_rigid_rot,refold_pos
   !
   !
   TYPE(md_type),POINTER                   :: md_obj =>NULL()
   TYPE(bfgs_type),POINTER                 :: bfgs_obj => NULL() 
   CHARACTER(LEN=*),PARAMETER              :: TAGNAME="ion_control"
   LOGICAL                                 :: bfgs_ispresent,md_ispresent
   ! 
   !
   IF (TRIM(ion_dynamics)=="bfgs") THEN
      ALLOCATE (bfgs_obj) 
      CALL qes_init (bfgs_obj,"bfgs",ndim=bfgs_ndim,trust_radius_min=trust_radius_min,&
                         trust_radius_max=trust_radius_max,trust_radius_init=trust_radius_init,&
                         w1=w_1,w2=w_2)
   ELSE IF(TRIM(ion_dynamics)=="verlet" .OR. TRIM(ion_dynamics)=="langevin" .OR. &
           TRIM(ion_dynamics) == "langevin-smc" ) THEN
      ALLOCATE(md_obj) 
      CALL qes_init (md_obj,"md",pot_extrapolation=pot_extrapolation,&
                      wfc_extrapolation=wfc_extrapolation,ion_temperature=ion_temperature,&
                      tolp=tolp,timestep=dt,deltaT=delta_t,nraise=nraise,tempw=tempw)
   END IF 
   CALL qes_init (obj,TAGNAME,ion_dynamics=TRIM(ion_dynamics), UPSCALE=upscale, REMOVE_RIGID_ROT=remove_rigid_rot,&
                  REFOLD_POS=refold_pos, BFGS=bfgs_obj, MD=md_obj)
   IF (ASSOCIATED(bfgs_obj)) THEN 
      CALL qes_reset (bfgs_obj)
      DEALLOCATE(bfgs_obj) 
   END IF 
   IF (ASSOCIATED(md_obj)) THEN   
      CALL qes_reset (md_obj)
      DEALLOCATE (md_obj) 
   END IF 
   !
   END SUBROUTINE qexsd_init_ion_control
   !
   !
   !------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_cell_control(obj,cell_dynamics, pressure, wmass,cell_factor,cell_dofree,iforceh)
   !------------------------------------------------------------------------------------------
   !
   IMPLICIT NONE
   ! 
   TYPE (cell_control_type)                     :: obj
   CHARACTER(LEN=*),INTENT(IN)                  :: cell_dynamics, cell_dofree
   REAL(DP),INTENT(IN)                          :: pressure, wmass, cell_factor
   INTEGER,DIMENSION(3,3),INTENT(IN)            :: iforceh
   ! 
   CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="cell_control"
   INTEGER,DIMENSION(3,3)                       :: my_forceh    
   !
   LOGICAL                                      :: fix_volume=.FALSE.,&
                                                   fix_area=.FALSE.,&
                                                   isotropic=.FALSE. 
   INTEGER                                      :: i,j
   TYPE(integerMatrix_type),TARGET              :: free_cell_obj
   TYPE(integerMatrix_type),POINTER             :: free_cell_ptr => NULL()  
   !
   IF (ANY(iforceh /= 1)) THEN 
      free_cell_ptr => free_cell_obj
      FORALL (i=1:3,j=1:3) my_forceh(i,j) = iforceh(i,j)
   END IF 
   SELECT CASE  (TRIM(cell_dofree))
      CASE ('all') 
         my_forceh = 1 
      CASE ('shape') 
         fix_volume = .TRUE.
      CASE ('2Dshape') 
         fix_area = .TRUE.
      CASE ('volume') 
         isotropic = .TRUE. 
      !CASE default 
         !NULLIFY ( free_cell_ptr) 
   END SELECT  
   IF (ASSOCIATED (free_cell_ptr)) CALL  qes_init (free_cell_obj,"free_cell",[3,3],my_forceh, ORDER = 'F' )
   !
   CALL qes_init (obj,TAGNAME, PRESSURE = pressure, CELL_DYNAMICS=cell_dynamics, WMASS=wmass, CELL_FACTOR=cell_factor,&
                  FIX_VOLUME=fix_volume, FIX_AREA=fix_area, ISOTROPIC=isotropic, FREE_CELL=free_cell_ptr)
   IF( ASSOCIATED(free_cell_ptr))   CALL qes_reset (free_cell_obj)
   END SUBROUTINE  qexsd_init_cell_control
   !
   !
   !-------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_symmetry_flags(obj,nosym,nosym_evc,noinv,no_t_rev,force_symmorphic,&
                                        use_all_frac)
   !-------------------------------------------------------------------------------------------
   !
   IMPLICIT NONE
   ! 
   TYPE ( symmetry_flags_type)                       :: obj
   LOGICAL,INTENT(IN)                                :: nosym,nosym_evc,noinv,no_t_rev,&
                                                        force_symmorphic,use_all_frac
   ! 
   CHARACTER(LEN=*),PARAMETER                        :: TAGNAME="symmetry_flags"
   CALL qes_init (obj,TAGNAME,nosym=nosym,nosym_evc=nosym_evc,noinv=noinv,&
                                no_t_rev=no_t_rev,force_symmorphic=force_symmorphic,&
                                use_all_frac=use_all_frac)
   ! 
   END SUBROUTINE qexsd_init_symmetry_flags
   !
   ! 
   !--------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_boundary_conditions(obj,assume_isolated,esm_bc, fcp_opt, fcp_mu, esm_nfit,esm_w, esm_efield)
   !--------------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE (boundary_conditions_type)              :: obj
   CHARACTER(LEN=*),INTENT(IN)                  :: assume_isolated
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)         :: esm_bc
   LOGICAL,OPTIONAL,INTENT(IN)                  :: fcp_opt
   REAL(DP),OPTIONAL,INTENT(IN)                 :: fcp_mu
   INTEGER,OPTIONAL,INTENT(IN)                  :: esm_nfit
   REAL(DP),OPTIONAL,INTENT(IN)                 :: esm_w,esm_efield
   ! 
   TYPE (esm_type),POINTER                      :: esm_obj => NULL() 
   LOGICAL                                      :: esm_ispresent = .FALSE.
   CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="boundary_conditions"
   !
   IF ( TRIM(assume_isolated) .EQ. "esm" ) THEN 
      esm_ispresent = .TRUE. 
      ALLOCATE(esm_obj) 
      CALL qes_init (esm_obj,"esm",bc=TRIM(esm_bc),nfit=esm_nfit,w=esm_w,efield=esm_efield)
   END IF 
   CALL qes_init (obj,TAGNAME,ASSUME_ISOLATED =assume_isolated, FCP_OPT= fcp_opt, FCP_MU = fcp_mu, ESM = esm_obj)
   IF ( esm_ispresent ) THEN
      CALL qes_reset (esm_obj)
      DEALLOCATE(esm_obj) 
   END IF 
   END SUBROUTINE qexsd_init_boundary_conditions
   ! 
   !
   !--------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_ekin_functional(obj,ecfixed,qcutz,q2sigma)
   !--------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE (ekin_functional_type)                   :: obj
   REAL(DP),INTENT(IN)                           :: ecfixed,qcutz,q2sigma
   ! 
   CHARACTER(LEN=*),PARAMETER                    :: TAGNAME="ekin_functional"
   CALL qes_init (obj,TAGNAME,ecfixed=ecfixed,qcutz=qcutz,q2sigma=q2sigma)
   END SUBROUTINE qexsd_init_ekin_functional
   !
   ! 
   !---------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_external_atomic_forces(obj,extfor,nat)
   !
   TYPE(matrix_type)                           :: obj
   REAL(DP),DIMENSION(:,:),INTENT(IN)          :: extfor
   INTEGER,INTENT(IN)                          :: nat
   ! 
   CHARACTER(LEN=*),PARAMETER                  :: TAGNAME="external_atomic_forces"
   !
   CALL qes_init (obj,TAGNAME,[3,nat],mat=extfor, order = 'F' )
   END SUBROUTINE qexsd_init_external_atomic_forces
   !
   !     
   !-------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_free_positions(obj,if_pos,nat)
   !
   IMPLICIT NONE
   !
   TYPE(integerMatrix_type)             :: obj
   INTEGER,DIMENSION(:,:),INTENT(IN)    :: if_pos
   INTEGER,INTENT(IN)                   :: nat
   !
   CHARACTER(LEN=*),PARAMETER           :: TAGNAME = "free_positions" 
   REAL(DP),DIMENSION(:,:),ALLOCATABLE  :: free_positions
   ! 
   CALL qes_init (obj,TAGNAME,DIMS = [3,nat], MAT = if_pos, ORDER = 'F' )
   END SUBROUTINE qexsd_init_free_positions
   ! 
   !----------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_starting_atomic_velocities(obj,tv0rd,rd_vel,nat)
   !----------------------------------------------------------------------------------
   !
   IMPLICIT NONE
   ! 
   TYPE (matrix_type)                   :: obj
   LOGICAL,INTENT(IN)                   :: tv0rd
   REAL(DP),DIMENSION(:,:),INTENT(IN)   :: rd_vel
   INTEGER,INTENT(IN)                   :: nat
   ! 
   CHARACTER(LEN=*),PARAMETER           :: TAGNAME="starting_atomic_velocities"
   INTEGER                              :: xdim=0,ydim=0
   IF (tv0rd) THEN
      xdim=3
      ydim=nat
   END IF
   CALL qes_init (obj,TAGNAME,[xdim,ydim],rd_vel )
   END SUBROUTINE qexsd_init_starting_atomic_velocities
   ! 
   !-------------------------------------------------------------------------------------
   SUBROUTINE  qexsd_init_spin_constraints(obj,constrained_magnetization,lambda,&
                                          fixed_magnetization)
   !-------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE(spin_constraints_type)                :: obj
   CHARACTER(LEN=*),INTENT(IN)                :: constrained_magnetization
   REAL(DP),INTENT(IN)                        :: lambda
   REAL(DP),DIMENSION(3),OPTIONAL,INTENT(IN)  :: fixed_magnetization
   ! 
   CHARACTER(LEN=*),PARAMETER                 :: TAGNAME="spin_constraints"
   REAL(DP),DIMENSION(3)                      :: target_magnetization=0.d0
   ! 
   IF (PRESENT(fixed_magnetization)) target_magnetization=fixed_magnetization
   CALL  qes_init (obj,TAGNAME,SPIN_CONSTRAINTS=TRIM(constrained_magnetization),&
                                   TARGET_MAGNETIZATION=fixed_magnetization ,LAGRANGE_MULTIPLIER=lambda)
   END SUBROUTINE qexsd_init_spin_constraints
   !
   ! 
   !-------------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_electric_field_input (obj,tefield,dipfield,lelfield,lberry,edir,gdir,emaxpos,eopreg,eamp,    &
         efield,efield_cart,nberrycyc,nppstr, gate, zgate, relaxz, block, block_1, block_2, block_height )
   !---------------------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE (electric_field_type)                   :: obj
   LOGICAL,INTENT(IN)                           :: tefield,lelfield, lberry
   LOGICAL,OPTIONAL,INTENT(IN)                  :: dipfield  
   INTEGER,INTENT(IN),OPTIONAL                  :: edir,gdir,nberrycyc,nppstr
   REAL(DP),INTENT(IN),OPTIONAL                 :: emaxpos,eopreg,eamp
   REAL(DP),INTENT(IN),OPTIONAL                 :: efield
   REAL(DP),INTENT(IN),OPTIONAL,DIMENSION(3)    :: efield_cart
   LOGICAL,INTENT(IN),OPTIONAL                  :: gate, block,relaxz 
   REAL(DP),INTENT(IN),OPTIONAL                 :: zgate,block_1, block_2, block_height
   ! 
   CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="electric_field",&
                                                   SAWTOOTH="sawtooth_potential",&
                                                   HOMOGENEOUS="homogenous_field",&
                                                   BERRYPHASE="Berry_Phase"
   REAL(DP),POINTER                             :: efield_cart_loc(:)=>NULL(), electric_field_amplitude=>NULL() 
   INTEGER,POINTER                              :: electric_field_direction => NULL() 
   CHARACTER(LEN=256)                           :: electric_potential
   LOGICAL                                      :: gate_, block_
   REAL(DP)                                     :: block_1_, block_2_, block_3_
   TYPE(gate_settings_type),TARGET              :: gata_settings_obj
   TYPE(gate_settings_type),POINTER             :: gata_settings_ptr => NULL() 
   TARGET                                       :: eamp, edir, efield, gdir 
   ! 
   electric_potential = "none"
   IF (tefield) THEN  
      electric_potential=SAWTOOTH
      electric_field_amplitude=>eamp
      electric_field_direction=>edir
   ELSE  IF (lelfield) THEN
      electric_potential=HOMOGENEOUS
      IF (PRESENT(efield)) electric_field_amplitude => efield
      IF ( gdir .GT. 0 )  electric_field_direction  => gdir
   ELSE IF (lberry) THEN
      electric_potential=BERRYPHASE
      IF ( gdir .GT. 0) electric_field_direction    =>  gdir
   END IF  
   IF (PRESENT (gate)) THEN 
      gata_settings_ptr => gata_settings_obj 
      CALL qes_init (gata_settings_obj, "gate_settings", gate, zgate, relaxz,&
         block, block_1, block_2, block_height ) 
   END IF 
   CALL  qes_init ( obj, TAGNAME, electric_potential=electric_potential, dipole_correction = dipfield,   &
                    electric_field_direction=electric_field_direction, potential_max_position = emaxpos, &
                    potential_decrease_width = eopreg, electric_field_amplitude=electric_field_amplitude,&                
                    electric_field_vector = efield_cart, n_berry_cycles=nberrycyc, nk_per_string=nppstr, &
                    gate_settings = gata_settings_obj)
   END SUBROUTINE qexsd_init_electric_field_input
   !
   !----------------------------------------------------------------------------------------------------------
   SUBROUTINE  qexsd_init_atomic_constraints(obj,ion_dynamics,lconstrain,nconstr,constr_type,constr_tol,     &
                                            constr_target,constr)
   !----------------------------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE (atomic_constraints_type)                           :: obj
   CHARACTER(LEN=*),INTENT(IN)                             :: ion_dynamics
   LOGICAL,INTENT(IN)                                      :: lconstrain
   INTEGER,OPTIONAL,INTENT(IN)                             :: nconstr
   REAL(DP),OPTIONAL,INTENT(IN)                            :: constr(:,:)
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)                    :: constr_type(:)
   REAL(DP),OPTIONAL,INTENT(IN)                            :: constr_target(:),constr_tol
   ! 
   CHARACTER(LEN=*),PARAMETER                             :: TAGNAME="atomic_constraints"
   TYPE(atomic_constraint_type),ALLOCATABLE                :: constr_objs(:)
   INTEGER                                                 :: iconstr
   !
   !
   ALLOCATE (constr_objs(nconstr))
   DO iconstr=1,nconstr
      CALL qes_init (constr_objs(iconstr),"atomic_constraint", constr_parms=constr(:,iconstr),&
                                     constr_type=TRIM(constr_type(iconstr)),constr_target=constr_target(iconstr))
   END DO
   CALL    qes_init (obj,TAGNAME, num_of_constraints=nconstr, atomic_constraint=constr_objs,tolerance=constr_tol)
   DO iconstr=1,nconstr
      CALL qes_reset (constr_objs(iconstr))
   END DO
   DEALLOCATE (constr_objs)
   END SUBROUTINE qexsd_init_atomic_constraints
   !
   !------------------------------------------------------------------------------------------------------------ 
   SUBROUTINE qexsd_init_occupations(obj, occupations, nspin) 
   !------------------------------------------------------------------------------------------------------------
     !
     IMPLICIT NONE
     TYPE(occupations_type),INTENT(OUT)        :: obj
     CHARACTER(LEN=*),INTENT(IN)               :: occupations
     INTEGER,INTENT(IN)                        :: nspin
     ! 
     INTEGER                                   :: spin_degeneracy
     !   
     IF (nspin .GT. 1) THEN 
        spin_degeneracy = 1
     ELSE 
        spin_degeneracy = 2
     END IF
     CALL  qes_init (obj, "occupations", occupations = TRIM(occupations))
   END SUBROUTINE qexsd_init_occupations
   !
   !---------------------------------------------------------
   SUBROUTINE qexsd_init_smearing(obj, smearing, degauss)
   !---------------------------------------------------------
      !
      IMPLICIT NONE
      TYPE(smearing_type),INTENT(OUT)      :: obj
      CHARACTER(LEN = * ), INTENT(IN)      :: smearing
      REAL(DP),INTENT(IN)                  :: degauss 
      ! 
      CHARACTER(LEN=256)                   :: smearing_local 

      SELECT CASE (TRIM  (smearing))
        CASE ("gaussian", "gauss")
            smearing_local="gaussian"
        CASE ('methfessel-paxton', 'm-p', 'mp')
            smearing_local="mp"
        CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt') 
            smearing_local="mv"
        CASE ('fermi-dirac', 'f-d', 'fd') 
            smearing_local="fd"
        CASE default
            smearing_local='not set'
      END SELECT 
      CALL qes_init (obj,"smearing",degauss=degauss,smearing=smearing_local)
      !
      END SUBROUTINE qexsd_init_smearing
      !--------------------------------------------------------------------------------------------
      !
END MODULE qexsd_input          
  
