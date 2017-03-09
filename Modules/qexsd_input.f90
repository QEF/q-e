!
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__XSD)
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
  USE iotk_base,        ONLY : iotk_indent, iotk_maxindent
  USE constants,        ONLY : e2,bohr_radius_angs
  USE iotk_module
  USE qes_module
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
  SELECT CASE ( TRIM (calculation)) 
    CASE ('scf', 'nscf', 'bands') 
       IF ( nstep == 1) THEN 
          nstep_ispresent = .FALSE. 
       ELSE 
          nstep_ispresent = .TRUE. 
       END IF 
   CASE DEFAULT 
       IF ( nstep == 50 ) THEN 
          nstep_ispresent = .FALSE. 
       ELSE 
          nstep_ispresent = .TRUE.
       END IF 
  END SELECT 
  !
  CALL qes_init_control_variables(obj,tagname,title=title,calculation=calculation,&
                                  restart_mode=restart_mode,prefix=prefix,        &
                                  pseudo_dir=pseudo_dir,outdir=outdir,disk_io=disk_io_value,&
                                  verbosity=TRIM(verbosity_value),stress=stress,forces=forces,    &
                                  wf_collect=wf_collect,max_seconds=int_max_seconds,  &
                                  etot_conv_thr=etot_conv_thr,forc_conv_thr=forc_conv_thr, &
                                  press_conv_thr=press_conv_thr,print_every=iprint, NSTEP = nstep, &
                                  NSTEP_ISPRESENT = nstep_ispresent )

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
  
  CALL qes_init_spin(obj,TAGNAME,lsda=lsda,noncolin=noncolin,spinorbit=spinorbit)
  
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
  INTEGER,INTENT(IN)                           :: nbnd,nspin
  CHARACTER(LEN=*),INTENT(IN)                  :: occupations,smearing
  REAL(DP),INTENT(IN)                          :: degauss,tot_charge
  REAL(DP),DIMENSION(:),OPTIONAL,INTENT(IN)    :: input_occupations, input_occupations_minority
  REAL(DP),OPTIONAL,INTENT(IN)                 :: tot_mag
  !
  CHARACTER(25)                                :: smearing_local
  INTEGER                                      :: spin_degeneracy, inpOcc_size = 0
  CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="bands"
  TYPE(smearing_type)                          :: smearing_obj
  TYPE(occupations_type)                       :: occup_obj
  TYPE(inputoccupations_type),ALLOCATABLE      :: inpOcc_objs(:)
  LOGICAL                                      :: tot_mag_ispresent = .FALSE., &
                                                  inp_occ_arepresent = .FALSE.
  ! 
  IF (TRIM(occupations) .NE. "smearing")  THEN
     CALL qes_init_smearing ( smearing_obj, "smearing", degauss=0.d0, smearing="")
     smearing_obj%lread  = .FALSE.
     smearing_obj%lwrite = .FALSE.
  ELSE
     SELECT CASE (TRIM  (smearing))
       CASE ("gaussian", "gauss")
           smearing_local="gaussian"
       CASE ('methfessel-paxton', 'm-p', 'mp')
           smearing_local="mp"
       CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv') 
           smearing_local="mv"
       CASE ('fermi-dirac', 'f-d', 'fd') 
           smearing_local="fd"
     END SELECT 
     CALL qes_init_smearing(smearing_obj,"smearing",degauss=degauss,smearing=smearing_local)
  END IF
  IF (nspin .GT. 1) THEN 
     spin_degeneracy = 1
  ELSE 
     spin_degeneracy = 2
  END IF
  CALL  qes_init_occupations(occup_obj, "occupations", spin= spin_degeneracy, & 
                              spin_ispresent =.FALSE., occupations = TRIM(occupations))
  !
  IF (PRESENT(input_occupations) ) THEN 
     inp_occ_arepresent = .TRUE.
     SELECT CASE ( nspin)
       CASE (2) 
          inpOcc_size=2
       CASE default
          inpOcc_size=1
     END SELECT
     ALLOCATE (inpOcc_objs(inpOcc_size))
     IF ( inpOcc_size .GT. 1) THEN 
        CALL qes_init_inputOccupations( inpOcc_objs(1),"input_occupations", 1, &
                  REAL(spin_degeneracy,KIND=DP),nbnd-1, input_occupations(2:nbnd) ) 
        CALL qes_init_inputOccupations( inpOcc_objs(2),"input_occupations", 2, & 
                  REAL(spin_degeneracy,KIND=DP) , nbnd-1,input_occupations_minority(2:nbnd))
     ELSE 
        CALL qes_init_inputOccupations( inpOcc_objs(1),"input_occupations", 1,            &
                                        REAL(spin_degeneracy,KIND=DP) , nbnd-1, input_occupations(2:nbnd) )   
     END IF
  END IF
  !
  IF (PRESENT ( tot_mag)) tot_mag_ispresent = .TRUE.
        
  CALL qes_init_bands(obj,TAGNAME,NBND_ISPRESENT=(nbnd .GT. 0), NBND = nbnd,&
                      SMEARING_ISPRESENT = smearing_obj%lread, SMEARING = smearing_obj,& 
                      TOT_CHARGE_ISPRESENT=.TRUE., TOT_CHARGE = tot_charge,  &
                      TOT_MAGNETIZATION_ISPRESENT = tot_mag_ispresent, TOT_MAGNETIZATION = tot_mag, &
                      OCCUPATIONS=occup_obj, INPUTOCCUPATIONS_ISPRESENT=inp_occ_arepresent, &
                      NDIM_INPUTOCCUPATIONS= inpOcc_size, INPUTOCCUPATIONS = inpOcc_objs)
  CALL qes_reset_smearing(smearing_obj)
  CALL qes_reset_occupations(occup_obj)
  IF (inp_occ_arepresent) THEN 
     CALL qes_reset_inputoccupations(inpocc_objs(1))
     IF (inpOcc_size .GT. 1 ) CALL qes_reset_inputoccupations(inpocc_objs(2))
     DEALLOCATE (inpocc_objs)
  END IF
  !
  END SUBROUTINE qexsd_init_bands
  !
  !
  !--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_basis(obj,k_points,ecutwfc,ecutrho,nr1,nr2,nr3,nr1s,nr2s,nr3s,nr1b,nr2b,nr3b)
  !--------------------------------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE (basis_type)                         :: obj
  CHARACTER(LEN=*),INTENT(IN)               :: k_points
  REAL(DP),INTENT(IN)                       :: ecutwfc,ecutrho
  INTEGER,INTENT(IN)                        :: nr1,nr2,nr3,nr1s,nr2s,nr3s,nr1b,nr2b,nr3b
  ! 
  TYPE(basisSetItem_type)                  :: grid_obj,smooth_grid_obj,box_obj
  CHARACTER(LEN=*),PARAMETER                :: TAGNAME="basis",FFT_GRID="fft_grid",FFT_SMOOTH="fft_smooth",&
                                               FFT_BOX="fft_box"
  LOGICAL                                   :: fft_grid_ispresent=.FALSE.,&
                                               fft_smooth_ispresent=.FALSE.,&
                                               fft_box_ispresent   = .FALSE., &
                                               gamma_only=.FALSE., ecutrho_ispresent=.FALSE.
  IF( ( nr1 .NE. 0 ) .AND. ( nr2 .NE. 0 ) .AND. ( nr3 .NE. 0 )) THEN
    fft_grid_ispresent=.TRUE.
    CALL qes_init_basisSetItem(grid_obj,FFT_GRID,nr1,nr2,nr3,"grid set in input")
  END IF
  ! 
  IF( ( nr1s .NE. 0 ) .AND. ( nr2s .NE. 0 ) .AND. ( nr3s .NE. 0 )) THEN
    fft_smooth_ispresent=.TRUE.
    CALL qes_init_basisSetItem(smooth_grid_obj,FFT_SMOOTH,nr1s,nr2s,nr3s,"grid set in input")
  END IF
  ! 
  IF( ( nr1b .NE. 0 ) .AND. ( nr2b .NE. 0 ) .AND. ( nr3b .NE. 0 )) THEN
    fft_box_ispresent=.TRUE.
    CALL qes_init_basisSetItem(box_obj,FFT_BOX,nr1b,nr2b,nr3b,"grid set in input")
  END IF
  ! 
  IF (TRIM(k_points) .EQ. "gamma" ) gamma_only=.TRUE.
  IF (ecutrho .GT. 4.d0*ecutwfc) ecutrho_ispresent=.TRUE. 

  CALL qes_init_basis(obj,TAGNAME,gamma_only_ispresent=gamma_only,gamma_only=gamma_only,ecutwfc=ecutwfc,               &
                      ecutrho=ecutrho,ecutrho_ispresent=ecutrho_ispresent,fft_grid_ispresent=fft_grid_ispresent,       &
                      fft_grid=grid_obj,fft_smooth_ispresent=fft_smooth_ispresent,fft_smooth=smooth_grid_obj,           &
                      fft_box=box_obj,fft_box_ispresent=fft_box_ispresent)
  !
  IF (fft_grid_ispresent)   CALL  qes_reset_basisSetItem( grid_obj )
  IF (fft_smooth_ispresent) CALL  qes_reset_basisSetItem( smooth_grid_obj )
  IF ( fft_box_ispresent )  CALL  qes_reset_basisSetItem( box_obj )
  ! 
  !
  !
  END SUBROUTINE qexsd_init_basis
  !-------------------------------------------------------------------------------------------
  SUBROUTINE qexsd_init_electron_control( obj,diagonalization,mixing_mode,mixing_beta,&
                                          conv_thr, mixing_ndim, max_nstep, tqr,tq_smoothing, &
                                          tbeta_smoothing, & 
                                          diago_thr_init, diago_full_acc, diago_cg_maxiter,&
                                          diago_david_ndim)
  !-------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE(electron_control_type)             ::  obj
  CHARACTER(LEN=*),INTENT(IN)             :: diagonalization,mixing_mode
  REAL(DP),INTENT(IN)                     :: mixing_beta, conv_thr, diago_thr_init
  INTEGER,INTENT(IN)                      :: mixing_ndim,max_nstep,diago_cg_maxiter,&
                                             diago_david_ndim
  LOGICAL,INTENT(IN)                      :: diago_full_acc,tqr, tq_smoothing, tbeta_smoothing
  !
  CHARACTER(LEN=*),PARAMETER              :: TAGNAME="electron_control"
  !
  CALL qes_init_electron_control(obj,TAGNAME,diagonalization=diagonalization,&
                                mixing_mode=mixing_mode,mixing_beta=mixing_beta,&
                                conv_thr=conv_thr,mixing_ndim=mixing_ndim,max_nstep=max_nstep,&
                                tq_smoothing= tq_smoothing, tbeta_smoothing = tbeta_smoothing,& 
                                real_space_q=tqr,diago_thr_init=diago_thr_init,& 
                                diago_full_acc=diago_full_acc,diago_cg_maxiter=diago_cg_maxiter)
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
   TYPE(monkhorst_pack_type)            :: mpack_obj
   TYPE(k_point_type),ALLOCATABLE       :: kp_obj(:)
   LOGICAL                              :: mpack_ispresent,kp_ispresent
   CHARACTER(LEN=100)                   :: kind_of_grid
   INTEGER                              :: kdim,ik,jk,kcount
   REAL(DP),DIMENSION(3)                :: my_xk
   REAL(DP)                             :: scale_factor
   !
  
   IF (TRIM(k_points).EQ."automatic") THEN 
      !
      IF ((s1+s2+s3).EQ.0) THEN
         kind_of_grid="Monkhorst-Pack"
      ELSE
         kind_of_grid="Uniform grid with offset"
      END IF
      CALL qes_init_monkhorst_pack(mpack_obj,"monkhorst_pack",nk1,nk2,nk3,&
                                   s1,s2,s3,kind_of_grid)
      CALL qes_init_k_points_IBZ(obj,TAGNAME,monkhorst_pack_ispresent=.TRUE.,&
                                 monkhorst_pack=mpack_obj,nk_ispresent=.FALSE.,&
                                 nk=0,k_point_ispresent=.FALSE.,ndim_k_point=0,k_point=kp_obj)
      CALL qes_reset_monkhorst_pack(mpack_obj)
   ELSE
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
          CALL qes_init_k_point(kp_obj(kcount),"k_point",1.d0,.TRUE.,LABEL= "", LABEL_ISPRESENT=.FALSE., &
                                K_POINT = xk(:,1))
          kcount=kcount+1
          DO ik=1,nk-1
             DO jk=1,NINT(wk(ik))
                my_xk=xk(:,ik)+(DBLE(jk)/wk(ik))*(xk(:,ik+1)-xk(:,ik))
                my_xk=my_xk*scale_factor
                CALL qes_init_k_point(kp_obj(kcount),"k_point",1.d0,.TRUE.,LABEL="", LABEL_ISPRESENT = .FALSE., &
                                      K_POINT = my_xk)
                kcount=kcount+1
             END DO
          END DO
      ELSE
          kdim=nk
          ALLOCATE  (kp_obj(kdim))      
          DO ik=1,kdim
             my_xk=xk(:,ik)*scale_factor
             CALL qes_init_k_point(kp_obj(ik),"k_point",wk(ik),.TRUE.,label="",label_ispresent=.FALSE.,K_POINT=my_xk)
          END DO
      END IF
      CALL qes_init_k_points_IBZ(obj,TAGNAME,monkhorst_pack_ispresent=.FALSE.,&
                                 monkhorst_pack=mpack_obj,nk_ispresent=.TRUE.,nk=kdim,&
                                 k_point_ispresent=.TRUE.,ndim_k_point=kdim,k_point=kp_obj)
      DO ik = 1,kdim
         CALL qes_reset_k_point(kp_obj(ik))
      END DO
      DEALLOCATE (kp_obj)
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
   REAL(DP),INTENT(IN)                     :: upscale,tempw,tolp,delta_t,trust_radius_min,trust_radius_max,&
                                              trust_radius_init,w_1,w_2
   INTEGER,INTENT(IN)                      :: nraise,bfgs_ndim
   REAL(DP),INTENT(IN)                     :: dt
   LOGICAL,INTENT(IN)                      :: remove_rigid_rot,refold_pos
   !
   !
   TYPE(md_type)                           :: md_obj
   TYPE(bfgs_type)                         :: bfgs_obj
   CHARACTER(LEN=*),PARAMETER              :: TAGNAME="ion_control"
   LOGICAL                                 :: bfgs_ispresent,md_ispresent
   ! 
   !
   IF (TRIM(ion_dynamics)=="bfgs") THEN
      bfgs_ispresent=.TRUE.
      md_ispresent=  .FALSE.
      CALL qes_init_bfgs(bfgs_obj,"bfgs",ndim=bfgs_ndim,trust_radius_min=trust_radius_min,&
                         trust_radius_max=trust_radius_max,trust_radius_init=trust_radius_init,&
                         w1=w_1,w2=w_2)
   ELSE IF(TRIM(ion_dynamics)=="verlet" .OR. TRIM(ion_dynamics)=="langevin" .OR. &
           TRIM(ion_dynamics) == "langevin-smc" ) THEN
      bfgs_ispresent=.FALSE.
      md_ispresent=.TRUE.
      CALL qes_init_md(md_obj,"md",pot_extrapolation=pot_extrapolation,&
                      wfc_extrapolation=wfc_extrapolation,ion_temperature=ion_temperature,&
                      tolp=tolp,timestep=dt,deltaT=delta_t,nraise=nraise,tempw=tempw)
   ELSE
      bfgs_ispresent=.FALSE.
      md_ispresent  =.FALSE.
   END IF 
   CALL qes_init_ion_control(obj,TAGNAME,ion_dynamics=TRIM(ion_dynamics),upscale_ispresent=bfgs_ispresent,&
                             upscale=upscale,remove_rigid_rot_ispresent=.true.,&
                             remove_rigid_rot=remove_rigid_rot,refold_pos_ispresent=.TRUE.,&
                             refold_pos=refold_pos,bfgs_ispresent=bfgs_ispresent,bfgs=bfgs_obj,&
                             md_ispresent=md_ispresent,md=md_obj)
   IF (bfgs_ispresent) CALL qes_reset_bfgs(bfgs_obj)
   IF (md_ispresent)   CALL qes_reset_md(md_obj)
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
                                                   fix_volume_ispresent=.FALSE.,&
                                                   fix_area=.FALSE.,&
                                                   fix_area_ispresent=.FALSE.,&
                                                   isotropic=.FALSE.,&
                                                   isotropic_ispresent=.FALSE.,&
                                                   free_cell_ispresent=.TRUE.
   INTEGER                                      :: i,j
   TYPE(integerMatrix_type)                     :: free_cell_obj
   !
   FORALL (i=1:3,j=1:3) my_forceh(i,j) = iforceh(i,j)
   IF (TRIM(cell_dofree)=='default') THEN
      free_cell_ispresent=.FALSE.
      my_forceh=1
   ELSE IF (TRIM(cell_dofree)=='all' ) THEN 
      my_forceh=1
   ELSE IF (TRIM(cell_dofree)=='shape') THEN 
      fix_volume=.TRUE.
      fix_volume_ispresent=.TRUE.
   ELSE IF ( TRIM(cell_dofree)=='2Dshape') THEN 
      fix_area = .TRUE.
      fix_area_ispresent=.TRUE.
   ELSE IF (TRIM(cell_dofree)=='volume') THEN
      isotropic=.TRUE.
      isotropic_ispresent=.TRUE.  
   END IF
   IF (free_cell_ispresent) CALL  qes_init_integerMatrix(free_cell_obj,"free_cell",3,3,my_forceh)
   !
   CALL qes_init_cell_control(obj,TAGNAME, PRESSURE = pressure, CELL_DYNAMICS=cell_dynamics, WMASS_ISPRESENT=.TRUE.,&
                              WMASS=wmass, CELL_FACTOR_ISPRESENT=.TRUE., CELL_FACTOR=cell_factor,&
                              FIX_VOLUME_ISPRESENT=fix_volume_ispresent,FIX_VOLUME=fix_volume,&
                              FIX_AREA_ISPRESENT=fix_area_ispresent, FIX_AREA=fix_area,& 
                              ISOTROPIC_ISPRESENT=isotropic_ispresent,ISOTROPIC=isotropic,&
                              FREE_CELL_ISPRESENT=free_cell_ispresent, FREE_CELL=free_cell_obj)
   IF( free_cell_ispresent ) CALL qes_reset_integerMatrix(free_cell_obj)
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
   CALL qes_init_symmetry_flags(obj,TAGNAME,nosym=nosym,nosym_evc=nosym_evc,noinv=noinv,&
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
   TYPE (esm_type)                              :: esm_obj
   LOGICAL                                      :: esm_ispresent = .FALSE., fcp_opt_ispresent = .TRUE., &
                                                   fcp_mu_ispresent = .FALSE. , fcp_opt_ = .FALSE.
   REAL(DP)                                     :: fcp_mu_ = 0.d0  
   CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="boundary_conditions"
   !
   IF ( TRIM(assume_isolated) .EQ. "esm" ) THEN 
      esm_ispresent = .TRUE. 
      CALL qes_init_esm(esm_obj,"esm",bc=TRIM(esm_bc),nfit=esm_nfit,w=esm_w,efield=esm_efield)
      IF ( PRESENT(fcp_opt) ) THEN 
          fcp_opt_ = fcp_opt
          fcp_mu_ispresent = .TRUE. 
          IF ( fcp_opt_ .AND. PRESENT ( fcp_mu)) fcp_mu_ = fcp_mu
      END IF 
   END IF 
   CALL qes_init_boundary_conditions(obj,TAGNAME,ASSUME_ISOLATED =assume_isolated, &
                                     FCP_OPT_ISPRESENT = fcp_opt_ispresent, FCP_OPT= fcp_opt_, &
                                     FCP_MU_ISPRESENT = fcp_mu_ispresent, FCP_MU = fcp_mu_, &  
                                     ESM_ISPRESENT = esm_ispresent, ESM = esm_obj)
   IF ( esm_ispresent ) CALL qes_reset_esm(esm_obj)
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
   CALL qes_init_ekin_functional(obj,TAGNAME,ecfixed=ecfixed,qcutz=qcutz,q2sigma=q2sigma)
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
   CALL qes_init_matrix(obj,TAGNAME,ndim1_mat=3,ndim2_mat=nat,mat=extfor)
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
   CALL qes_init_integerMatrix(obj,TAGNAME,ndim1_int_mat=3,ndim2_int_mat=nat,int_mat=if_pos(:,1:nat))
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
   CALL qes_init_matrix(obj,TAGNAME,xdim,ydim,rd_vel)
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
   CALL  qes_init_spin_constraints(obj,TAGNAME,spin_constraints=TRIM(constrained_magnetization),&
                                   target_magnetization_ispresent=PRESENT(fixed_magnetization), &
                                   target_magnetization=target_magnetization,lagrange_multiplier=lambda)
   END SUBROUTINE qexsd_init_spin_constraints
   !
   ! 
   !-------------------------------------------------------------------------------------------------
   SUBROUTINE qexsd_init_electric_field_input (obj,tefield,dipfield,lelfield,lberry,edir,gdir,emaxpos,eopreg,eamp,    &
                                               efield,efield_cart,nberrycyc,nppstr)
   !---------------------------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   TYPE (electric_field_type)                   :: obj
   LOGICAL,INTENT(IN)                           :: tefield,lelfield,dipfield,lberry
   INTEGER,INTENT(IN),OPTIONAL                  :: edir,gdir,nberrycyc,nppstr
   REAL(DP),INTENT(IN),OPTIONAL                 :: emaxpos,eopreg,eamp
   REAL(DP),INTENT(IN),OPTIONAL                 :: efield
   REAL(DP),INTENT(IN),OPTIONAL,DIMENSION(3)    :: efield_cart
   ! 
   CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="electric_field",&
                                                   SAWTOOTH="sawtooth_potential",&
                                                   HOMOGENEOUS="homogenous_field",&
                                                   BERRYPHASE="Berry_Phase"
   REAL(DP)                                     :: emaxpos_loc=0.d0,eopreg_loc=0.d0,electric_field_amplitude=0.d0
   REAL(DP),DIMENSION(3)                        :: efield_cart_loc=0.d0
   INTEGER                                      :: electric_field_direction,nberrycyc_loc=0,nppstr_loc=0
   CHARACTER(LEN=256)                           :: electric_potential
   LOGICAL                                      :: dir_ispresent=.FALSE., amp_ispresent= .FALSE.,&
                                                   nberrycyc_ispresent=.FALSE.,nppstr_ispresent=.FALSE., &
                                                   electric_field_ispresent = .FALSE.
   ! 
   IF (tefield) THEN  
      electric_potential=SAWTOOTH
      emaxpos_loc=emaxpos
      eopreg_loc=eopreg
      electric_field_amplitude=eamp
      electric_field_direction=edir
      dir_ispresent=.TRUE.
      amp_ispresent=.TRUE.
   ELSE  IF (lelfield) THEN
      electric_potential=HOMOGENEOUS
      nberrycyc_loc = nberrycyc
      nberrycyc_ispresent = .TRUE.
      nppstr_loc = nppstr
      nppstr_ispresent = .TRUE.
      IF (PRESENT(efield_cart)) THEN 
         efield_cart_loc=efield_cart
         electric_field_ispresent = .TRUE.
      END IF
      IF (PRESENT(efield)) THEN
         electric_field_amplitude = efield
         amp_ispresent = .TRUE.
      END IF
      IF ( gdir .GT. 0 ) THEN 
         dir_ispresent = .TRUE. 
         electric_field_direction = gdir
      END IF       
   ELSE IF (lberry) THEN
      electric_potential=BERRYPHASE
      nberrycyc_loc=nberrycyc
      nppstr_ispresent = .TRUE.
      nppstr_loc = nppstr
      IF ( gdir .GT. 0) THEN 
         dir_ispresent=.TRUE.
         electric_field_direction = gdir
      END IF
   END IF  
      
   CALL  qes_init_electric_field( obj, TAGNAME, electric_potential=electric_potential,        &
                                dipole_correction_ispresent=dipfield, dipole_correction = dipfield, &
                                electric_field_direction_ispresent= dir_ispresent, &
                                electric_field_direction=electric_field_direction,&
                                potential_max_position_ispresent=tefield, potential_max_position=emaxpos_loc,  &
                                potential_decrease_width = eopreg_loc, potential_decrease_width_ispresent=tefield,  &
                                electric_field_amplitude=electric_field_amplitude, &                
                                electric_field_amplitude_ispresent=amp_ispresent, &
                                electric_field_vector = efield_cart_loc,                                     &
                                electric_field_vector_ispresent= electric_field_ispresent, &
                                n_berry_cycles_ispresent=nberrycyc_ispresent,n_berry_cycles=nberrycyc_loc,&
                                nk_per_string_ispresent=nppstr_ispresent,nk_per_string=nppstr_loc  )
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
      CALL qes_init_atomic_constraint(constr_objs(iconstr),"atomic_constraint", constr_parms=constr(:,iconstr),&
                                     constr_type=TRIM(constr_type(iconstr)),constr_target=constr_target(iconstr))
   END DO
   CALL    qes_init_atomic_constraints(obj,TAGNAME,num_of_constraints=nconstr,ndim_atomic_constraint=nconstr,   &
                                          atomic_constraint=constr_objs,tolerance=constr_tol)
   DO iconstr=1,nconstr
      CALL qes_reset_atomic_constraint(constr_objs(iconstr))
   END DO
   DEALLOCATE (constr_objs)
   END SUBROUTINE qexsd_init_atomic_constraints
   !
   !------------------------------------------------------------------------------------------------------------ 
END MODULE qexsd_input          
         
#else
! 
MODULE qexsd_input
  IMPLICIT NONE
  INTEGER :: dummy__
END MODULE qexsd_input
! 
#endif
  
