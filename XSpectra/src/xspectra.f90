!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM X_Spectra
  USE kinds, ONLY : DP
  USE constants,       ONLY : rytoev,pi,fpi
  USE io_global,       ONLY : stdout,ionode,ionode_id   ! Modules/io_global.f90
  USE io_files,        ONLY : prefix, iunwfc, nwordwfc, tmp_dir, diropn
  USE cell_base,       ONLY : bg, at, celldm
  USE parameters,      ONLY : ntypx,lmaxx,lqmax
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
  USE ktetra,          ONLY : ltetra, ntetra, tetra
  USE start_k,         ONLY : nk1, nk2, nk3, k1, k2, k3
  USE wvfct,           ONLY : npwx,nbnd,npw,igk,et, wg! et(nbnd,nkstot)
  USE radial_grids,    ONLY : ndmx
  USE atom,            ONLY : rgrid
  USE becmod,          ONLY : becp
  USE uspp,            ONLY : vkb, nkb, okvan 
  USE uspp_param,      ONLY : upf
  USE xspectra
  USE ener,            ONLY : ef, ef_up, ef_dw !Fermi energy (ef in Ry)
  USE symm_base,       ONLY : nsym,s
  USE paw_gipaw,       ONLY : read_recon,  &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,          &
       paw_recon,           &
       set_paw_upf
  USE klist,           ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points for local pool
       nelec,nelup,neldw,             & !number of electrons
       xk,                & ! k-points coordinates
       wk ,               & ! k-points weight
       npk,               &
       degauss,lgauss,ngauss,    &
       two_fermi_energies
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE noncollin_module,ONLY : noncolin
  USE mp,              ONLY : mp_bcast, mp_sum             !parallelization
  USE mp_global,       ONLY : mp_startup, mp_global_end
  USE mp_pools,        ONLY : intra_pool_comm, npool
  USE mp_world,        ONLY : nproc, world_comm
  USE control_flags,   ONLY : gamma_only
  USE environment,     ONLY : environment_start

  USE cut_valence_green, ONLY :&
       cut_ierror, &    ! convergence tolerance for one step in the integral
       cut_stepu , &    ! integration initial step, upper side
       cut_stepl , &    ! integration initial step, lower side
       cut_startt, &    ! integration start value of the t variable
       cut_tinf  , &    ! maximum value of the lower integration boundary
       cut_tsup  , &    ! minimum value of the upper integration boudary
       cut_desmooth,&   ! size of the interval near the fermi energy 
                        ! in which cross section is smoothed
       cut_nmemu,&      ! size of the memory of the values of the green function, upper side
       cut_nmeml,&      ! size of the memory of the values of the green function, lower side
       cut_occ_states  ! true if you want tou remove occupied states from the spectrum

  USE control_flags,   ONLY : twfcollect
  !<CG>
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
                                 gamma_lines, gamma_tab, gamma_points, &
                                 gamma_mode, gamma_file
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm, init_xspectra_paw_nhm
  USE edge_energy, ONLY: getE
  !</CG>

  IMPLICIT NONE 
  !
  ! ... local variables
  !
  INTEGER :: nargs,iiarg,ierr,ios,il,ibnd,ibnd_up,ibnd_dw !,xm_r,nc_r,ncomp_max
  INTEGER :: iargc
  !INTEGER :: nt,nb,na,i,j,k,nrest,nline
  INTEGER :: nt,na,i,j,k !nrest,nline
  INTEGER, ALLOCATABLE :: ncalcv(:,:)
  INTEGER, ALLOCATABLE :: paw_iltonhb(:,:,:)
              ! corresp l, projector, type <--> cumulative over all the species
  REAL (DP) :: norm, core_energy
  REAL (DP) :: rc(ntypx,0:lmaxx),r_paw(0:lmaxx)
  REAL (DP) :: core_wfn(ndmx)
  REAL (DP) :: ehomo, elumo, middle_gap ! in eV 
  REAL (DP) :: e_core, e_core_ryd
  REAL (DP), ALLOCATABLE :: a(:,:,:),b(:,:,:),xnorm(:,:)      !lanczos vectors
  !REAL (DP), EXTERNAL   :: efermig,efermit
 
  LOGICAL :: exst !, loc_set

  CHARACTER (LEN=10)  :: dummy_char
  CHARACTER (LEN=256) :: filerecon(ntypx)
  !REAL (DP) :: gamma_energy(2), gamma_value(2)
  REAL (DP) :: xeps_dot_xk ! scalar product between xepsilon and xkvec
  !REAL(dp) :: auxrpol(3,2) non used variable


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Initialize MPI environment, clocks and a few other things
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'XSpectra' )

  CALL banner_xspectra()

  CALL read_input_and_bcast(filerecon, r_paw)

  call write_sym_param_to_stdout()

  call select_nl_init(edge, nl_init, two_edges, n_lanczos)     

  CALL start_clock( calculation  )

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Reads, initializes and writes several things in two cases :
  ! $   case 1: xonlyplot=.false.  (complete calc., i.e. Lanczos + plot)
  ! $   case 2: xonlyplot=.true.   (the plot only)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(.NOT.xonly_plot) THEN

     WRITE(stdout,1000) ! return+line
     WRITE(stdout,'(5x,2a)')  &
     '                 Reading SCF save directory: ',trim(prefix)//'.save'
     WRITE(stdout,1001) ! line+return

     CALL read_file()

     CALL calculate_and_write_homo_lumo_to_stdout(ehomo,elumo)

     call reset_k_points_and_reinit_nscf()

     IF(xread_wf) THEN
        IF (okvan) THEN
           WRITE(stdout,'(/,5x,a,f10.6,/)') &
                        'Approx. ram memory needed per proc in MB = ',&
                        (16.0*4.0*npwx*npool)/(nproc*1000000.0)
        ELSE
           WRITE(stdout,'(/,5x,a,f10.6,/)') &
                        'Approx. ram memory needed per proc in MB = ',&
                        (16.0*3.0*npwx*npool)/(nproc*1000000.0)
        ENDIF
     ENDIF

        ! ... normalize xepsilon

     IF ( xcoordcrys ) CALL cryst_to_cart(1,xepsilon,at,1)
     norm=DSQRT(xepsilon(1)**2+xepsilon(2)**2+xepsilon(3)**2)
     DO i=1,3
       xepsilon(i)=xepsilon(i)/norm
    ENDDO

    ! ... If needed normalize xkvec
    !     and check orthogonality with xepsilon
    
    IF(xang_mom.EQ.2) THEN
       
       IF ( xcoordcrys ) CALL cryst_to_cart(1,xkvec,at,1)
       norm=DSQRT(xkvec(1)**2+xkvec(2)**2+xkvec(3)**2)
       xkvec(:)=xkvec(:)/norm
       xeps_dot_xk=xkvec(1)*xepsilon(1)+&
                    xkvec(2)*xepsilon(2)+&
                    xkvec(3)*xepsilon(3)
       IF ((ABS(xeps_dot_xk)) >= 1.0d-6) THEN
          WRITE(stdout,'(5x,a)') &
               'ERROR: xkvec and xepsilon are not orthogonal'
          WRITE(stdout,'(12x,a,f10.6,/)') 'scalar product=', xeps_dot_xk
          WRITE(stdout,'(5x,a)') 'STOP'
          CALL stop_xspectra ()
       ENDIF
       
    ENDIF
    

 
     !... Is the type associated to xiabs existing ?
     
     i=0
     DO na=1,nat
        IF(ityp(na).EQ.xiabs) i=i+1
     ENDDO
     IF(i.NE.1) THEN
        CALL errore( 'main program', 'Wrong xiabs!!!',i)
     ENDIF

     !...  Reads reconstruction files
     WRITE(stdout,1000) ! return+line 
     WRITE(stdout,'(5x,a)')  &
     '          Reading core wavefunction file for the absorbing atom'
     WRITE(stdout,1001) ! line+return
     DO nt = 1, ntyp
        call set_paw_upf(nt, upf(nt))
     ENDDO

     CALL read_core_abs(filecore,core_wfn, nl_init)

     WRITE(stdout,'(5x,a," successfully read")') TRIM(ADJUSTL(filecore))
     
     IF ( .NOT. paw_recon(xiabs)%gipaw_data_in_upf_file ) &
          CALL read_recon ( filerecon(xiabs), xiabs, paw_recon(xiabs) ) !*apsi

     !... Assign paw radii to species (this will become soon obsolete)
     !
     !  read recon should be parallelized
     !
     WRITE(stdout,1000) ! return+line 
     WRITE(stdout,'(5x,a)')  &
     '                         Attributing the PAW radii '
     WRITE(stdout,'(5x,a)')  &
     '                for the absorbing atom [units: Bohr radius]'
     WRITE(stdout,1001) ! line+return

     DO nt=1,ntyp
        IF ((.NOT.paw_recon(nt)%gipaw_data_in_upf_file)) then
           paw_recon(nt)%paw_nbeta=0
        ENDIF
        DO  j=1,paw_recon(nt)%paw_nbeta
           il=paw_recon(nt)%psphi(j)%label%l
           IF(xiabs.EQ.nt.AND.DABS(r_paw(il)).lt.1.d-6) THEN
              !*apsi  if(r(kkbeta(nt),nt).GT.1.d-3) THEN
              !*apsi     rc(nt,il)=  r(kkbeta(nt),nt)
              !*apsi  ELSE
              !*apsi     WRITE(stdout,*) 'Warning-no-kkbeta r_paw(',il,')=1.0'
              !*apsi     rc(nt,il)=1.0
              !*apsi  ENDIF
              !<CG>  to be verified
              IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
                 rc(nt, il)= paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
                 WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                      'PAW proj', j,': r_paw(l=',il,')=',rc(nt,il),'  (1.5*r_cut)'
              ELSE
                 rc(nt, il)= 1.5d0
                 WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                      'PAW proj',j,': r_paw(l=',il,')=',rc(nt,il),'  (set to 1.5)'
              ENDIF
           ELSEIF(xiabs.EQ.nt.AND.DABS(r_paw(il)).GT.1.d-6) THEN
              rc(nt,il)=r_paw(il)
           ELSEIF(nt.NE.xiabs) THEN
              !*apsi  IF(r(kkbeta(nt),nt).GT.1.d-3) THEN
              !*apsi     rc(nt,il)=r(kkbeta(nt),nt)
              !*apsi  ELSE
              !*apsi     rc(nt,il)=1.0
              !*apsi  ENDIF
              !<CG> to be verified
              IF(paw_recon(nt)%psphi(j)%label%rc.GT.1.d-3) THEN
                 rc(nt,il)=paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
              ELSE
                 rc(nt,il)=1.5
              ENDIF
              !</CG> 
           ENDIF
           
        ENDDO
     ENDDO
     WRITE(stdout,*)
     WRITE(stdout,'(8x,a)')&
           'NB: The calculation will not necessary use all these r_paw values.'
     WRITE(stdout,'(8x,a)')&
           '    - For a edge in the electric-dipole approximation,'
     WRITE(stdout,'(8x,a)')&
           '      only the r_paw(l=1) values are used.'
     WRITE(stdout,'(8x,a)')&
           '    - For a K edge in the electric-quadrupole approximation,'
     WRITE(stdout,'(8x,a,/)')'      only the r_paw(l=2) values are used.'
     WRITE(stdout,'(8x,a,/)')     '    - For a L2 or L3 edge in the electric-quadrupole approximation,'
     WRITE(stdout,'(8x,a,/)')'      all projectors (s, p and d) are used.'

     !<CG>
     DO nt=1,ntyp
        DO j = 1,paw_recon(nt)%paw_nbeta
           paw_recon(nt)%psphi(j)%label%rc=rc(nt,paw_recon(nt)%psphi(j)%label%l)
           paw_recon(nt)%aephi(j)%label%rc=rc(nt,paw_recon(nt)%aephi(j)%label%l)
        ENDDO
     ENDDO
     !</CG>


     !...  write band energies if xread_wf=true

     IF(xread_wf.AND.TRIM(verbosity).EQ.'high') THEN

        WRITE(stdout,'(5x,a)') &
           'Band energies read from scf save file [units: eV]'
        WRITE(stdout,'(5x,a)') &
           '-------------------------------------------------'
        DO i=1,nkstot
           WRITE(stdout,'(5x,"k=[",3f14.8,"]   spin=",1i2)') &
                xk(1,i),xk(2,i),xk(3,i),isk(i)
           WRITE(stdout, '(8f9.4)') (et(j,i)*rytoev,j=1,nbnd)
        ENDDO

        WRITE(stdout,*)

     ENDIF

     CALL init_gipaw_1
  
          
 
     !        CALL mp_bcast( ef, ionode_id )  !Why should I need this ?


     !... Allocation of variables for paw projectors
     !... and initialization of Vanderbilt and Paw projectors
     !
     ! CG : becp allocated in the xanes_dipole and xanes_quadrupole subroutines
     !     call allocate_bec_type ( nkb, nbnd, becp )

     !CALL init_gipaw_1 ! Already called above

     !...  Definition of a specific indexation to avoid M. Profeta's crazy one 
     CALL init_xspectra_paw_nhm

     ALLOCATE (paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp))

     CALL define_index_arrays(paw_iltonhb)

     !... Allocates PAW projectors
     ALLOCATE (paw_vkb( npwx,  paw_nkb))
     ALLOCATE (paw_becp(paw_nkb, nbnd))

  ELSEIF(xonly_plot) THEN  

     CALL read_header_save_file(x_save_file)
     nks = nkstot
     IF(lsda) THEN
        isk(1:nkstot/2)=1
        isk(nkstot/2+1:nkstot)=2
     !  wk(1:nkstot)=2.d0/nkstot
     ELSEIF(.NOT.lsda) THEN
        isk(1:nkstot)=1
     !  wk(1:nkstot)=2.d0/nkstot
     ENDIF
     wk(1:nkstot)=2.d0/nkstot
     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )

  ENDIF   

  !<CG>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Computing gamma tabulated values
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !<MCB> THIS MUST BE CHANGED
  !<DC> why ? </DC>

  IF (TRIM(gamma_mode).EQ.'file') THEN
     CALL read_gamma_file
  ELSEIF (TRIM(gamma_mode).EQ.'variable') THEN
     gamma_lines=2
     ALLOCATE(gamma_points(2,2))
     gamma_points(1,1)=gamma_energy(1)
     gamma_points(2,1)=gamma_energy(2)
     gamma_points(1,2)=gamma_value(1)
     gamma_points(2,2)=gamma_value(2)
  ENDIF

  IF ((TRIM(gamma_mode).EQ.'file').OR.(TRIM(gamma_mode).EQ.'variable')) THEN
     ALLOCATE( gamma_tab(xnepoint))
     CALL initialize_gamma_tab
     DEALLOCATE(gamma_points)
  ENDIF

  !</CG>

  xnitermax=xniter

  !<DC> Block moved below into the previous "if xonlyplot" block
  !IF(xonly_plot) THEN
  !   CALL read_header_save_file(x_save_file)
  !   nks = nkstot
  !   WRITE(6,*) 'nks=',nks
  !   IF(lsda) THEN
  !      isk(1:nkstot/2)=1
  !      isk(nkstot/2+1:nkstot)=2
  !      wk(1:nkstot)=2.d0/nkstot
  !   ELSEIF(.NOT.lsda) THEN
  !      isk(1:nkstot)=1
  !      wk(1:nkstot)=2.d0/nkstot
  !   ENDIF
  !   CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $  Checks PAW relations between pseudo partial waves and projector 
  ! $  (radial parts)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(.NOT.xonly_plot.AND.TRIM(verbosity).EQ.'high') &
      CALL check_paw_projectors(xiabs)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Allocates and initializes Lanczos variables
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ALLOCATE(a(xnitermax,n_lanczos,nks))
  ALLOCATE(b(xnitermax,n_lanczos,nks))
  ALLOCATE(xnorm(n_lanczos,nks))
  ALLOCATE(ncalcv(n_lanczos,nks))

  a(:,:,:)=0.d0
  b(:,:,:)=0.d0
  xnorm(:,:)=0.d0
  ncalcv(:,:)=0

  ! for restart
  ALLOCATE(calculated(n_lanczos,nks))
  calculated(:,:)=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $  And now we go...  XANES CALCULATION
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(TRIM(calculation).EQ.'xanes') THEN
     IF(.NOT.xonly_plot) THEN ! complte calculation
        WRITE(stdout, 1000) ! line 
        WRITE(stdout,'(5x,a)')&
        '                     Starting XANES calculation'
        IF(nl_init(2).eq.0) then
           IF (xang_mom==1) WRITE(stdout,'(5x,a)')&
           '                in the electric dipole approximation'
           IF (xang_mom==2) WRITE(stdout,'(5x,a)')&
           '              in the electric quadrupole approximation'
           WRITE(stdout,1001)  ! line 
        ELSEIF(nl_init(2).eq.1) then
           WRITE(stdout,'(5x,a)')&
                '                in the electric dipole approximation'
        ENDIF
           
        IF(TRIM(ADJUSTL(restart_mode)).eq.'restart') THEN
          CALL read_header_save_file(x_save_file)
          CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
        ENDIF

        !...  Writes information about the method

        WRITE(stdout,'(7(5x,a,/))') &
        "Method of calculation based on the Lanczos recursion algorithm",&
        "--------------------------------------------------------------",&
        "   - STEP 1: Construction of a kpoint-dependent Lanczos basis,",&
        "     in which the Hamiltonian is tridiagonal (each 'iter' ",&
        "     corresponds to the calculation of one more Lanczos vector)",&
        "   - STEP 2: Calculation of the cross-section as a continued fraction",&
        "     averaged over the k-points."

        WRITE(stdout,'(5x,"... Begin STEP 1 ...",/)')

        save_file_version = 2 ! adds ef after core_energy

        IF(nl_init(2).EQ.0.AND.xang_mom .eq. 1) THEN

           save_file_kind='xanes_dipole'
           CALL xanes_dipole(a,b,ncalcv,xnorm,core_wfn,&
                             paw_iltonhb,terminator,verbosity)

        ELSEIF(nl_init(2).EQ.0.AND.xang_mom .eq. 2) THEN

           save_file_kind='xanes_quadrupole'
           CALL xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,&
                                 paw_iltonhb,terminator,verbosity)

        ELSEIF(nl_init(2).EQ.1.AND.xang_mom .eq. 1) then
           call xanes_dipole_general_edge(a,b,ncalcv,nl_init,xnorm,core_wfn,paw_iltonhb,terminator, verbosity)
        ENDIF

        CALL write_save_file(a,b,xnorm,ncalcv,x_save_file)
        WRITE(stdout,'(/,5x,a)') &
              'Results of STEP 1 successfully written in x_save_file'
        WRITE(stdout,'(5x,a18,/,5x,a2,2x,a65)') 'x_save_file name: ',&
                      '->', x_save_file 
        WRITE(stdout,'(5x,a21,i2)') 'x_save_file version: ', save_file_version

        WRITE(stdout,'(/,5x,"... End STEP 1 ...",/)')

     ELSE ! only the spectrum plot
        WRITE(stdout,1000) ! return+line
        WRITE(stdout,'(5x,a)')&
        '                          Reading x_save_file'
        WRITE(stdout,1001) ! line+return

        CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)

        IF (save_file_version .eq. 1 .AND. abs(xe0-xe0_default)<1.e-3) THEN
          WRITE(stdout,'(5x,3a)') &
            "STOP: the variable 'xe0' must be assigned in ", &
             'input file', ' (since save_file_version = 1)'
          CALL stop_xspectra()
        ENDIF
     ENDIF

     IF (TRIM(save_file_kind).eq.'unfinished') CALL stop_xspectra ()       
 
     IF (.NOT.xonly_plot) THEN
        !
        WRITE(stdout,'(5x,"... Begin STEP 2 ...",/)')
        WRITE(stdout,'(5x,a)') &
         'The spectrum is calculated using the following parameters:'
        !
        e_core = getE(upf(xiabs)%psd,edge)

     ELSE
        WRITE(stdout,1000) ! return+line
        WRITE(stdout,'(5x,a,a)')  &
        '               Starting the calculation of the spectrum'
        WRITE(stdout,1001) !line+return
        WRITE(stdout,'(5x,a)') 'Using the following parameters:'
        !
        e_core = core_energy
        !
     ENDIF

     e_core_ryd = e_core/rytoev

     IF (abs(xe0-xe0_default)<1.d-3) THEN ! xe0 not in input_file
        write(stdout,'(8x,a,f9.4,a)') 'energy-zero of the spectrum [eV]: ',&
          ef 
        ! ef is either read in x_save_file (version 2) if xonlyplot=T
        ! or determined from the scf calculation
        xe0_ry = ef/rytoev    ! ef est en eV ici
     ELSE ! xe0 is read in the input
        WRITE(stdout,'(8x,a,f9.4,a)') 'xe0 [eV]: ', xe0
        xe0_ry=xe0/rytoev
     ENDIF

     IF (cut_occ_states) THEN
        WRITE(stdout,'(8x,a)') 'the occupied states are cut'
     ELSE
        WRITE(stdout,'(8x,a)') 'the occupied states are NOT cut'
     ENDIF
     WRITE(stdout,'(8x,a,f6.2)') 'xemin [eV]: ', xemin
     WRITE(stdout,'(8x,a,f6.2)') 'xemax [eV]: ', xemax
     WRITE(stdout,'(8x,a,i4)')   'xnepoint: ', xnepoint
     IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
        WRITE(stdout,'(8x,a,f8.3)')'constant broadening parameter [eV]: ', xgamma
     ELSE
        WRITE(stdout,'(8x,a)') 'energy-dependent broadening parameter:'
        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
           WRITE(stdout,'(8x,a,a30)')' -> using gamma_file: ', gamma_file
        ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
           WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
                 ' -> first, constant up to point (', gamma_energy(1), &
                 ',', gamma_value(1), ') [eV]'
           WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
                 ' -> then, linear up to point (', gamma_energy(2),   &
                 ',', gamma_value(2), ') [eV]'
           WRITE(stdout,'(8x,a)') ' -> finally, constant up to xemax'
        ENDIF
     ENDIF
     WRITE(stdout,'(8x,"Core level energy [eV]:",1x,g11.4)') -e_core
     WRITE(stdout,'(8x,a,/)') &
     ' (from electron binding energy of neutral atoms in X-ray data booklet)'
     
     IF(xang_mom.EQ.1) THEN
        CALL plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e_core_ryd,1)
        if(two_edges) CALL plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e_core_ryd,2)
     ELSEIF(xang_mom.EQ.2) THEN
        CALL plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,e_core_ryd)
     ELSEIF(nl_init(2).eq.1) then
        
     ENDIF
     !
     WRITE(stdout,'(5x,"Cross-section successfully written in ",a,/)') &
      'xanes.dat'
     IF (.NOT. xonly_plot) WRITE(stdout,'(5x,"... End STEP 2 ...",/)')

  ELSEIF(TRIM(calculation).EQ.'rxes') THEN
     CALL errore( 'Main', 'rxes Not yet implemented',1)
  ELSEIF(TRIM(calculation).EQ.'bethe_salpeter') THEN
     CALL errore( 'Main', 'bethe_salpeter Not yet implemented',1)
  ENDIF
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  IF(.NOT.xonly_plot) THEN
  !     call deallocate_bec_type ( becp )
  !  ENDIF

  DEALLOCATE(a)
  DEALLOCATE(b)
  DEALLOCATE(xnorm)
  DEALLOCATE(ncalcv)


  !WRITE (stdout,*) 'End program ', TRIM(calculation)

  CALL stop_clock( calculation  )
  CALL print_clock( calculation )

  WRITE (stdout, 1000)
  WRITE (stdout,'(5x,a)') '                           END JOB XSpectra'
  WRITE (stdout, 1001)
  CALL stop_xspectra () 

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Formats 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 1000 FORMAT(/,5x,&
  '-------------------------------------------------------------------------')
 1001 FORMAT(5x,&
  '-------------------------------------------------------------------------',&
  /)

END program X_Spectra

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE stop_xspectra
  !----------------------------------------------------------------------------
  !
  ! Synchronize processes before stopping. This is a copy of stop_pp.
  !
  USE control_flags, ONLY: twfcollect
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  USE parallel_include
  !
#ifdef __MPI

  INTEGER :: info
  LOGICAL :: op

  INQUIRE ( iunwfc, opened = op )

  IF ( op ) THEN
     IF (twfcollect) THEN
        CLOSE (unit = iunwfc, status = 'delete')
     ELSE
        CLOSE (unit = iunwfc, status = 'keep')
     ENDIF
  ENDIF

  CALL mp_global_end()

#endif

  STOP
END SUBROUTINE stop_xspectra


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE read_core_abs(filename,core_wfn,nl_init)
  !--------------------------------------------------------------------------
  USE kinds,   ONLY: DP
  USE atom,    ONLY: rgrid
  !USE atom,        ONLY : mesh     !mesh(ntypx) number of mesh points
  USE xspectra,ONLY: xiabs
  USE io_global,       ONLY : ionode, stdout

  IMPLICIT NONE

  INTEGER :: i, ierr, nbp, iblind
  INTEGER, dimension(2) :: nl_init
  CHARACTER (LEN=80) :: filename
  REAL(KIND=dp):: x
  REAL(KIND=dp):: core_wfn(*)

  !WRITE(6,*) 'xmesh=',rgrid(xiabs)%mesh
  open(unit=33,file=filename,form='formatted',iostat=ierr,status='old',err=123)

123   if( ierr == 29 )  then
          if( ionode ) write(stdout, &
          '("ERROR: core wavefunction file ",A,&
         & " does not exist or is not located in the right folder !")') &
           trim(filename)
          call stop_xspectra
        else if( ierr /= 0 ) then
          if( ionode ) write(stdout, '("ERROR reading the core wavefunction file")')
          call stop_xspectra
        else
          continue
        end if

  rewind(33)

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !<OB>           Brute force identification of the core state:
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  nbp = rgrid(xiabs)%mesh
! Determine how many lines to skip before reading the good core function
  select case( nl_init(1) )
    case(1)
      iblind = 1                       !1s
    case(2)
      if( nl_init(2) == 0 ) then       !2s
        iblind = nbp + 2
      else                             !2p
        iblind = 2*nbp + 3
      end if
    case(3)
      if( nl_init(2) == 0 ) then       !3s
        iblind = 3*nbp + 4
      else if( nl_init(2) == 1 ) then  !3p
        iblind = 4*nbp + 5
      else
        iblind = 5*nbp + 6             !3d
      end if
  end select

  do i = 1, iblind
    READ(33,*) 
  end do

  DO i=1,nbp
   READ(33 ,*) x,core_wfn(i)
  ENDDO

  close(33)

END SUBROUTINE read_core_abs

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE define_index_arrays(paw_iltonhb)
  !----------------------------------------------------------------------------
  USE paw_gipaw,   ONLY: paw_lmaxkb, paw_recon
  USE ions_base,   ONLY: ntyp => nsp
  USE parameters,  ONLY: lmaxx
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm ! CG

  IMPLICIT NONE
  ! Arguments
  INTEGER :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm,ntyp) ! CG

  ! Local
  INTEGER :: nt,ih,nb,l
  INTEGER :: ip_per_l(0:lmaxx) 

  DO nt = 1, ntyp
     ih = 1
     ip_per_l(:) = 0
     DO nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        ip_per_l(l) = ip_per_l(l) + 1
        paw_iltonhb(l,ip_per_l(l),nt) = ih
        ih = ih + 2*l + 1
     ENDDO
  ENDDO

END SUBROUTINE define_index_arrays

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION lastterm(a,b,g)
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp)    :: a, b, g, y1, y2, z1, z2, r
  COMPLEX(dp) :: lastterm

  y1 = a*a - g*g - 4*b
  y2 = -2*a*g
  r  = 0.5*SQRT(y1*y1 + y2*y2)

  IF (g<0) THEN
     z1 =  a/2 + 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 + 0.5*SQRT(-y1/2 + r)
  ELSE
     z1 =  a/2 - 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 - 0.5*SQRT(-y1/2+r)
  ENDIF
 
 lastterm = CMPLX(z1,z2,kind=DP)

END FUNCTION lastterm

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION paste_fermi(e,e0,a,b,gamma,m,term,first)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY: dp
  USE constants, ONLY: tpi
  USE xspectra,  ONLY: xnitermax, xcheck_conv
  USE cut_valence_green, ONLY: cut_ierror, cut_stepu, cut_stepl, &
                               cut_startt, cut_tinf, cut_tsup,   &
                               cut_nmemu, cut_nmeml, memu, meml
  IMPLICIT NONE

  REAL(dp) :: paste_fermi
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e, e0
  LOGICAL  :: term, first
  ! Local
  COMPLEX(dp) :: green, y, dy, c1, c2, e1, e2
  REAL(dp)    :: t, dt, t1, ta, tb
  INTEGER, save :: n1
  INTEGER, save :: n2
  INTEGER :: nn1, nn2

  IF (first) THEN
     memu(:,:) = (0.d0,0.d0)
     meml(:,:) = (0.d0,0.d0) 
     n1 = 0
     n2 = 0
     first = .false.
  ENDIF

  dy = cut_ierror + 1.0
  y = 0.d0

  nn1 = 1
  nn2 = 1

  t1 = 0.5773502692

  t = cut_startt

  DO WHILE ((abs(dy)>cut_ierror).OR.(t<cut_tsup))
     dt = cut_stepu*t
     ta = t + dt*(1-t1)/2
     tb = t + dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn1>n1) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn1<cut_nmemu) THEN
           memu(nn1,1) = c1
           memu(nn1,2) = c2
           n1 = nn1
        ENDIF
     ELSE
        c1 = memu(nn1,1)
        c2 = memu(nn1,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)          &
           + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
           + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
           + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP) )
     y = y + dy
     t = t + dt
     nn1 = nn1 + 1
  ENDDO

  t = cut_startt
  dy = cut_ierror + 1

  DO WHILE((abs(dy)>cut_ierror).OR.(t>cut_tinf))
     dt = cut_stepl * t
     ta = t - dt*(1-t1)/2
     tb = t - dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn2>n2) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn2<cut_nmeml) THEN
           meml(nn2,1) = c1
           meml(nn2,2) = c2
           n2 = nn2
        ENDIF
     ELSE
        c1 = meml(nn2,1)
        c2 = meml(nn2,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)           &
            + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
            + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
            + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP))
     y = y + dy
     t = t - dt
     nn2 = nn2 + 1
  ENDDO

  paste_fermi = AIMAG(y)/tpi

END FUNCTION paste_fermi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION green(a,b,e,m,term)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv

  IMPLICIT NONE

  COMPLEX(dp) :: green
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  COMPLEX(dp) :: e
  LOGICAL :: term
  ! Local
  INTEGER :: i, p, q
  COMPLEX(dp) :: res, lastterm 
  REAL(dp) :: aa, bb

  q = xcheck_conv/2
  IF (term) THEN
     aa = 0.0
     bb = 0.0
     DO p = 1, q
        aa = aa + a(m-p)
        bb = bb + b(m-p)
     ENDDO
     aa = aa/q
     bb = bb/q

     res = lastterm(aa-REAL(e), bb*bb, AIMAG(e))
  ELSE
     res = CMPLX(a(m)-REAL(e),AIMAG(e),kind=DP)
  ENDIF
  DO i = 1, m-1
     res = a(m-i) - e - b(m-i)*b(m-i)/res
  ENDDO

  green = 1/res

END FUNCTION green

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE check_paw_projectors(xiabs)
  !----------------------------------------------------------------------------
  USE kinds,           ONLY: DP
  USE constants,       ONLY: pi
  USE paw_gipaw,       ONLY: paw_lmaxkb, paw_recon
  USE xspectra_paw_variables, ONLY: xspectra_paw_nhm
  USE atom,            ONLY: rgrid, msh
  !  USE atom,  ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points              
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration 
  !       r, rab
  USE ions_base,       ONLY: ntyp => nsp
  USE io_global,       ONLY: stdout
  USE radin_mod

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: xiabs
  ! Local  variables
  INTEGER  :: nr, nrc, ip, jp, lmax, l, ip_l, jtyp, n1, n2, nrs, ndm, ih, jh
  REAL(dp) :: overlap, rexx, overlap2
  REAL (dp), ALLOCATABLE :: aux(:), f(:,:)
  REAL(dp) , ALLOCATABLE :: s(:,:), e(:), v(:,:)

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE(f(rgrid(xiabs)%mesh,2)) !allocation too big, it needs only up to msh

  WRITE(stdout,'(/,5x,a)')&
  '-------------------------------------------------------------------------'
  WRITE(stdout,'(5x,a)')  &
  '                         Verification of the PAW relations '
  WRITE(stdout,'(5x,a,/)')&
  '-------------------------------------------------------------------------'

  WRITE(stdout,'(8x,a)') 'atom type   total number of projectors'
  DO jtyp=1,ntyp
     WRITE (stdout,'(13x,i4,3x,i4)') jtyp, paw_recon(jtyp)%paw_nbeta
  ENDDO
  WRITE(stdout,*)
  
  !... I calculate maximum l
  !
  lmax=0
  DO ip = 1, paw_recon(xiabs)%paw_nbeta
     IF(paw_recon(xiabs)%psphi(ip)%label%l.GT.lmax) &
          lmax = paw_recon(xiabs)%psphi(ip)%label%l
  ENDDO

  WRITE(stdout,'(8x,a)') 'atom type    l   number of projectors per ang. mom.'
  DO jtyp = 1, ntyp
     DO l = 0, lmax
        WRITE(stdout,'(13x,i4,3x,i2,3x,i4)') jtyp, l, paw_recon(jtyp)%paw_nl(l)
     ENDDO
  ENDDO
  WRITE(stdout,*)

  !... We calculate the overlaps between partial waves and projectors
  !    to see if they are equal to the Croneker delta.

  nr = msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.

  WRITE(stdout,'(8x,a)')'Overlaps between partial waves and projectors (radial)'
  WRITE(stdout,'(8x,a)')'------------------------------------------------------'
  WRITE(stdout,*)
  WRITE(stdout,'(8x,a)') &
                    '< \tilde{phi}_{l,n} | \tilde{p}_{l,nn} > = delta_{n,nn}  ?'
  WRITE(stdout,*) 

  DO ip = 1, paw_recon(xiabs)%paw_nbeta
     DO jp = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l .EQ. &
                        paw_recon(xiabs)%psphi(jp)%label%l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc > nr) THEN
              WRITE(stdout,'(8x,a,i8,a,i8)') 'STOP: nrc=', nrc,' > nr=', nr
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc) = paw_recon(xiabs)%psphi(ip)%psi(1:nrc) &
                        * paw_recon(xiabs)%paw_betar(1:nrc,jp)
           aux(nrc+1:nr) = 0.d0
           WRITE(stdout,'(8x,"<tilde{phi}_",2i2,10X,"|tilde{p}_",2i2,">=",1f14.8)')  &
                ip,paw_recon(xiabs)%psphi(ip)%label%l,jp, &
                paw_recon(xiabs)%psphi(jp)%label%l, &
                para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr)
        ENDIF
     ENDDO
  ENDDO

  WRITE(stdout,*)

  WRITE(stdout,'(8x,a)')'Checking normalization of pseudo,ae wfc and projectors'
  WRITE(stdout,'(8x,a)')'(radial part only, integral up to r_c)'
  WRITE(stdout,'(8x,a)')'------------------------------------------------------'
  WRITE(stdout,*)
  WRITE(stdout,'(8x,a)') 'l,   n, |proj|^2, |pswf|^2 , |aewf|^2'
  DO l = 0, lmax
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc =Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           aux(1:nrc) = paw_recon(xiabs)%paw_betar(1:nrc,ip) &
                      * paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = paw_recon(xiabs)%aephi(ip)%psi(1:nrc) &
                      * paw_recon(xiabs)%aephi(ip)%psi(1:nrc)
           overlap2 = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = paw_recon(xiabs)%psphi(ip)%psi(1:nrc) &
                      * paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           rexx = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'(8x,2i4,3f14.8)')l,ip,overlap,overlap2,rexx
        ENDIF
     ENDDO
  ENDDO
  WRITE(stdout,*)
  
  GOTO 323

  WRITE(stdout,'(8x,a)') '<phi|chi> = \sum_nl <phi|phi_l> <p_l|chi>_nrc ?'
  WRITE(stdout,'(8x,a)') '-----------------------------------------------'
  WRITE(stdout,'(8x,a)') &
        'WARNING: this test assumes a form of the phi/chi function'
  
  !  DO l=0,lmax
  DO l = 1, 1
     ip_l = 0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(ip_l.EQ.0.AND.paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) ip_l = ip
     ENDDO
     
     f(:,:) = 0.d0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           f(1:nr,1) = f(1:nr,1) + &
                       paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
           f(1:nr,2) = f(1:nr,2) + &
                       1.123*paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
        ENDIF
     ENDDO
     rexx = 0.d0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc > nr) THEN
              WRITE(stdout,'(8x,a,i8,a,i8)') 'STOP: nrc=', nrc,' > nr=', nr
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc) = f(1:nrc,1)*paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = f(1:nrc,2)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           overlap = overlap*para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'(8x,"overlap(l=",1i2,",n=",1i2,")= ",1f14.8)')&
                          l, ip, overlap
           rexx = rexx + overlap
        ENDIF
     ENDDO
     aux(1:nr) = f(1:nr,1)*f(1:nr,2)
     WRITE(stdout,'(8x,"sum/overlap=",1f14.8)') &
                    rexx/para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,'(8x,"sum projectors=",1f14.8," overlap=",1f14.8)') &
                    rexx,para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,*)
  ENDDO
  ! ENDDO
  !WRITE(stdout,*)
  !WRITE(stdout,*) '================================================================'
  
323 CONTINUE
  !
  !  Check linear dependence of projectors
  !

  WRITE(stdout,'(8x,a)') 'Checking linear dependence of projectors'
  WRITE(stdout,'(8x,a)') '----------------------------------------'
  WRITE(stdout,*)

  DEALLOCATE(aux)

  ndm = MAXVAL (msh(1:ntyp))

  ALLOCATE(aux(ndm))

  DO l = 0, paw_lmaxkb
     IF (paw_recon(xiabs)%paw_nl(l)>0) THEN
        ALLOCATE (s(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        ALLOCATE (e(paw_recon(xiabs)%paw_nl(l)),v(paw_recon(xiabs)%paw_nl(l),&
                    paw_recon(xiabs)%paw_nl(l)))
        DO ih = 1, paw_recon(xiabs)%paw_nl(l)
           n1 = paw_recon(xiabs)%paw_iltonh(l,ih)
           nrc = paw_recon(xiabs)%psphi(n1)%label%nrc
           nrs = paw_recon(xiabs)%psphi(n1)%label%nrs
           DO jh = 1, paw_recon(xiabs)%paw_nl(l)
              n2 = paw_recon(xiabs)%paw_iltonh(l,jh)
              CALL step_f(aux,paw_recon(xiabs)%psphi(n1)%psi(1:msh(xiabs)) * &
                   paw_recon(xiabs)%psphi(n2)%psi(1:msh(xiabs)), &
                   rgrid(xiabs)%r(1:msh(xiabs)),nrs,nrc, 1.d0, msh(xiabs) )
              CALL simpson( msh(xiabs), aux, rgrid(xiabs)%rab(1), s(ih,jh))
           ENDDO
        ENDDO
  
        WRITE(stdout,'(8x,"atom type:",1i4)') xiabs
        WRITE(stdout,'(8x,"number of projectors projector  =",1i3, " angular momentum=",1i4)') &
                    paw_recon(xiabs)%paw_nl(l),l

        DO ih = 1, paw_recon(xiabs)%paw_nl(l)
           WRITE(stdout,'(8x,10f14.8)') (s(ih,jh),jh=1,paw_recon(xiabs)%paw_nl(l)) 
        ENDDO
        WRITE(stdout,'(8x,a)') 'Eigenvalues S matrix:'

        IF(paw_recon(xiabs)%paw_nl(l).EQ.1) THEN
           WRITE(stdout,'(10x,1i4,1f14.8)') 1,s(1,1)
        ELSE 
           CALL rdiagh(paw_recon(xiabs)%paw_nl(l), s, &
                       paw_recon(xiabs)%paw_nl(l) , e, v )
           DO ih=1,paw_recon(xiabs)%paw_nl(l)
              WRITE(stdout,'(10x,1i4,1f14.8)') ih,e(ih)
           ENDDO
        ENDIF
        WRITE(stdout,*)
        call flush(stdout)
        DEALLOCATE(s,e,v)
     ENDIF
  ENDDO
  !WRITE(stdout,*) '================================================================'

  DEALLOCATE(aux,f)

END SUBROUTINE check_paw_projectors

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
  !----------------------------------------------------------------------------
  ! This routine reads the x_save_file (default name: xanes.sav)
  ! (automatically named prefix_K.save for a K edge, in the next version). 
  ! This routine is used only if xonly_plot=.true.
  !----------------------------------------------------------------------------
  USE kinds,       ONLY: DP
  USE constants,   ONLY: rytoev
  USE klist,       ONLY: nks, nkstot
  USE xspectra,    ONLY: xnitermax, xang_mom, xiabs,   &
                         n_lanczos, save_file_version, &
                         save_file_kind, calculated,   &
                         xe0, xe0_default, xe0_ry
  USE ener,        ONLY: ef
  USE io_global,   ONLY: stdout, ionode
  USE lsda_mod,    ONLY: nspin, lsda

  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (INOUT) ::  a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (INOUT) ::  b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (INOUT) ::  xnorm(n_lanczos,nks) 
  INTEGER, INTENT (INOUT)  ::  ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  REAL(dp), INTENT (OUT)   ::  core_energy
  ! Local variables
  ! NB: The '_r' suffix means 'read in x_save_file'
  INTEGER  :: ierr, nkstot_r
  INTEGER  :: xm_r, nc_r, ncomp_max
  INTEGER  :: i, j, k, ncalcv_max
  INTEGER  :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp) :: xepsilon_r(3), xkvec_r(3)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:), aux(:,:)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))
  ALLOCATE(aux(xnitermax,nks))

  a(:,:,:) = 0.d0
  b(:,:,:) = 0.d0
  xnorm(:,:) = 0.d0
  ncalcv(:,:) = 0
  calculated_all(:,:) = 0


  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  WRITE(stdout,'(5x,"x_save_file name: ",a)') TRIM( x_save_file )
  REWIND(10)

  IF (save_file_version == 0) then
     WRITE(stdout,'(5x,a)') 'x_save_file version: old'
  !ELSEIF(save_file_version == 1) then
  !   WRITE(stdout,'(5x,a)') 'x_save_file version: 1'
  ELSE
     WRITE(stdout,'(5x,a,i3)') 'x_save_file version: ', save_file_version
     DO i = 1, 6 
        READ(10,*)      
     ENDDO
  ENDIF
  WRITE(stdout,*)

  READ(10,*) lsda, nspin
  WRITE(stdout,'(5x,a,i2)') 'nspin:',nspin
  READ(10,*) xm_r, nkstot_r, xnitermax
  WRITE(stdout,'(5x,a,i4)') 'number of k-points:',nkstot
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a,i4)') 'final-state angular momentum (xm_r): ',xm_r
  IF (xm_r==1) THEN
     WRITE(stdout,'(5x,a)') ' => electric-dipole approximation'
  ELSEIF (xm_r==2) THEN
     WRITE(stdout,'(5x,a)') ' => electric-quadrupole approximation'
  ELSE
     WRITE(stdout,'(5x,a)') 'Wrong value of xm_r: STOP'
     CALL stop_xspectra
  ENDIF
  IF(xm_r.NE.xang_mom) & 
     CALL errore('read_save_file','xm_r is different from xang_mom=',xang_mom)

  READ(10,*) ncalcv_max
  IF(ncalcv_max.GT.xnitermax) THEN
     WRITE(stdout,'(5x,a,i5)') 'ncalcv_max=',ncalcv_max
     CALL errore('read_save_file','ncalcv_max is grater than xnitermax=', &
                 xnitermax)
  ENDIF

  WRITE(stdout,*)
  IF (save_file_version < 2) THEN
     READ(10,*) core_energy
  ELSE
     READ(10,*) core_energy, ef
     WRITE(stdout,'(5x,a,f9.4)') 'Fermi level [eV]:', ef 
  ENDIF
  WRITE(stdout,'(5x,a,f10.3,/)') 'core energy [eV]:', core_energy

  READ(10,*) (xkvec_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !WRITE(stdout,*) 'xkvec read from savefile'
  !WRITE(stdout,*) (xkvec_r(i),i=1,3)
  READ(10,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) 'xepsilon read from file'
  !WRITE(stdout,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !write(stdout,*) 'n_lanczos=', n_lanczos
  WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') 'xepsilon [Cartesian frame]:',&
       (xepsilon_r(i),i=1,3)
  IF (xm_r==2) WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') &
       'xkvec [Cartesian frame]:', (xkvec_r(i),i=1,3)
  WRITE(stdout,*)

  DO i = 1, n_lanczos
     IF (TRIM(save_file_kind).eq.'unfinished') THEN
       READ(10,*) (calculated_all(i,j),j=1,nkstot)
     ENDIF
     READ(10,*) (xnorm_all(i,j),j=1,nkstot)
     READ(10,*) (ncalcv_all(i,j),j=1,nkstot)
     READ(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     READ(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
#ifdef __MPI
     CALL poolscatter(xnitermax,nkstot,a_all,nks,aux)
     a(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
     CALL poolscatter(xnitermax,nkstot,b_all,nks,aux)
     b(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
#else
     a(1:xnitermax,i,1:nkstot)=a_all(1:ncalcv_max,1:nkstot)
     b(1:xnitermax,i,1:nkstot)=b_all(1:ncalcv_max,1:nkstot)
#endif

  ENDDO
  CLOSE(10)

#ifdef __MPI 
  CALL poolscatter(n_lanczos,nkstot,xnorm_all,nks,xnorm)
  CALL ipoolscatter(n_lanczos,nkstot,ncalcv_all,nks,ncalcv)
  CALL ipoolscatter(n_lanczos,nkstot,calculated_all,nks,calculated)
#else
  IF(nks.NE.nkstot) THEN
     CALL errore('read_save_file','nks\=nkstot',1)
  ENDIF
  xnorm(1:n_lanczos,1:nkstot) = xnorm_all(1:n_lanczos,1:nkstot)
  ncalcv(1:n_lanczos,1:nkstot) = ncalcv_all(1:n_lanczos,1:nkstot)
  calculated(1:n_lanczos,1:nkstot) = calculated_all(1:n_lanczos,1:nkstot)
#endif
  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)
  DEALLOCATE(aux)

END SUBROUTINE read_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_header_save_file(x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY: DP
  USE klist,    ONLY: nkstot
  USE lsda_mod, ONLY: nspin,lsda
  USE xspectra, ONLY: save_file_version, save_file_kind, n_lanczos
  USE io_global,ONLY: stdout

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  INTEGER    :: ierr, nkstot_r
  INTEGER    :: xm_r
  CHARACTER  :: c

  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  REWIND(10)

  READ(10, '(a1)') c
  REWIND(10)
  IF (c == '#') then
     READ(10, '(20x,i8)') save_file_version
     READ(10, '(20x,a32)') save_file_kind
     READ(10,*) 
     READ(10,'(27x,i4)') n_lanczos
     READ(10,*) 
     READ(10,*) 
  ELSE
     save_file_version=0
     save_file_kind='xanes_old'
     n_lanczos=1
  ENDIF

  READ(10,*) lsda,nspin
  READ(10,*) xm_r,nkstot_r
  nkstot=nkstot_r
  CLOSE(10)

  RETURN 
END SUBROUTINE read_header_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_save_file(a,b,xnorm,ncalcv,x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: rytoev
  USE klist,      ONLY: nks, nkstot
  USE xspectra,   ONLY: xnitermax, xang_mom, xkvec, xepsilon, xiabs, &
                        save_file_version, save_file_kind,           &
                        n_lanczos, calculated, edge
  USE ener,       ONLY: ef
  USE io_global,  ONLY: ionode
  !*apsi  USE uspp_param, ONLY : psd
  USE lsda_mod,   ONLY: nspin,lsda
  USE uspp_param, ONLY: upf
  USE edge_energy, ONLY: getE
  !
  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (IN) :: a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (IN) :: b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (IN) :: xnorm(n_lanczos,nks)
  INTEGER,  INTENT (IN) :: ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file ! could be removed, 
                                                 ! exists in module
  ! Local variables
  INTEGER :: ierr
  INTEGER :: i, j, k 
  INTEGER :: ncalcv_max
  INTEGER :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:)
  CHARACTER(LEN=8) :: dte

  IF (ionode)  CALL DATE_AND_TIME(date=dte)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))

  ncalcv_all(:,:) = 0
  xnorm_all(:,:)  = 0.d0

  ncalcv_all(1:n_lanczos,1:nks) = ncalcv(1:n_lanczos,1:nks)
  xnorm_all(1:n_lanczos,1:nks)  = xnorm(1:n_lanczos,1:nks)
  calculated_all(1:n_lanczos,1:nks) = calculated(1:n_lanczos,1:nks)

#ifdef __MPI
  CALL poolrecover(xnorm_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(ncalcv_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(calculated_all,n_lanczos,nkstot,nks)
#endif

  ncalcv_max = 0
  DO j = 1, n_lanczos
     DO i = 1, nkstot       
        IF (ncalcv_all(j,i).GT.ncalcv_max) ncalcv_max = ncalcv_all(j,i)
     ENDDO
  ENDDO

  IF ( ionode ) THEN
     OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
            STATUS = 'UNKNOWN', IOSTAT = ierr )
     REWIND(10)
     WRITE(10, '(a20,i8)') '# save_file_version=',save_file_version
     WRITE(10, '(a20,a32)') '# save_file_kind   =',save_file_kind
     WRITE(10,'(a7,a8)') '# date=', dte
     WRITE(10,'(a27,i4)') '# number of lanczos stored=',n_lanczos
     WRITE(10,'(a1)') '#'
     WRITE(10,'(a1)') '#'
     WRITE(10,*) lsda,nspin
     WRITE(10,*) xang_mom,nkstot,xnitermax
     WRITE(10,*) ncalcv_max
!     WRITE(10,*) mygetK(upf(xiabs)%psd), ef ! ef here is in eV
     WRITE(10,*) getE(upf(xiabs)%psd,edge), ef ! ef here is in eV
     WRITE(10,*) (xkvec(i),i=1,3)
     WRITE(10,*) (xepsilon(i),i=1,3)
  ENDIF
  DO i=1, n_lanczos
     a_all(:,:) = 0.d0
     b_all(:,:) = 0.d0
     a_all(1:xnitermax,1:nks) =  a(1:xnitermax,i,1:nks)
     b_all(1:xnitermax,1:nks) =  b(1:xnitermax,i,1:nks)
#ifdef __MPI
     CALL poolrecover(a_all,xnitermax,nkstot,nks)
     CALL poolrecover(b_all,xnitermax,nkstot,nks)
#endif
     IF ( ionode) THEN
        IF (TRIM(save_file_kind).EQ.'unfinished') THEN
          WRITE(10,*) (calculated_all(i,j),j=1,nkstot)
        ENDIF
        WRITE(10,*) (xnorm_all(i,j),j=1,nkstot)
        WRITE(10,*) (ncalcv_all(i,j),j=1,nkstot)
        WRITE(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
        WRITE(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     ENDIF
  ENDDO
  CLOSE(10)

  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)

END SUBROUTINE write_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_status_of_the_code
  !----------------------------------------------------------------------------
  USE io_global,       ONLY : stdout
  IMPLICIT NONE
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') &
   '                      STATUS OF THE CODE (22/04/2009) '
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') 'Working features (22/04/2009)'
  WRITE (stdout,'(5x,a)') '-----------------------------'
  WRITE (stdout,'(5x,a)') &
   '- XANES works both in the electric-dipole and -quadrupole approximation,'
  WRITE (stdout,'(5x,a)') '- Spin polarized works'
  WRITE (stdout,'(5x,a)') '- DFT+U implemented, validated'
  WRITE (stdout,'(5x,a)') '- Ultrasoft pseudo works'
  WRITE (stdout,'(5x,a)') '- Cut occupied states working, improved'
  WRITE (stdout,'(5x,a)') '- Terminator working'
  WRITE (stdout,'(5x,a)') '- Multiprojectors TM+USPP working (MCB,CG)'
  WRITE (stdout,'(5x,a)') '- New save file format, with version numbering'
  WRITE (stdout,'(5x,a)') &
   '- Time limit implemented, with restart, seems to work'
  WRITE (stdout,'(5x,a)') &
   '- DFT+U tested ONLY for non ortho wfc, but implemented'
  WRITE(stdout,*) 
  WRITE (stdout,'(5x,a)') 'TO DO'
  WRITE (stdout,'(5x,a)') '-----'
  WRITE (stdout,'(5x,a)') '- L2,3 edges [OB]'
  WRITE (stdout,'(5x,a)') '- Generalization to all edges [OB]'
  WRITE (stdout,'(5x,a)') '- XMCD [?]'
  WRITE (stdout,'(5x,a)') '- IXS [DC]'
  WRITE (stdout,'(5x,a)') '- EELS [DC]'
  WRITE (stdout,'(5x,a)') '- REXS [DC]'
  WRITE (stdout,'(5x,a)') '- Bethe-Salpeter [?] '
  WRITE (stdout,'(5x,a)') '- RXES [?]' 

END SUBROUTINE write_status_of_the_code

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_gamma_file
  !----------------------------------------------------------------------------
  USE gamma_variable_mod

  IMPLICIT NONE
  INTEGER :: nl, ierr, i

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  CALL errore('io ', 'gamma file '//TRIM(gamma_file)//' not found', abs (ierr))
  REWIND(21)

  nl = 0

  DO
     READ (21,'(a1)',iostat=ierr)
     IF (ierr.NE.0) EXIT
     nl = nl + 1
  ENDDO
  CLOSE(21)

  gamma_lines = nl
  ALLOCATE(gamma_points(nl,2))

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  REWIND(21)

  DO i=1,nl
     READ(21,*) gamma_points(i,1), gamma_points(i,2)
  ENDDO

  CLOSE(21)

END SUBROUTINE read_gamma_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE initialize_gamma_tab
  !----------------------------------------------------------------------------
  USE kinds,     ONLY: dp
  USE xspectra,  ONLY: xemin, xemax, xnepoint
  USE io_global, ONLY: stdout
  USE constants, ONLY: rytoev
  USE gamma_variable_mod

  IMPLICIT NONE
  REAL(dp) :: e, x, y, dx
  INTEGER  :: i, j, n

  dx = (xemax-xemin)/dfloat(xnepoint)

  DO n = 1, xnepoint
     x = xemin + (n-1)*dx
     i = 1
     DO j = 1, gamma_lines
        IF(x > gamma_points(j,1)) i = i + 1
     ENDDO

     IF (i == 1) THEN
        y = gamma_points(1,2)
     ELSEIF (i == (gamma_lines+1)) THEN
        y = gamma_points(gamma_lines,2)
     ELSE
        y = (  gamma_points(i-1,2) * (gamma_points(i,1)-x)     &
             + gamma_points(i,2)   * (x-gamma_points(i-1,1))  )&
            / ( gamma_points(i,1) - gamma_points(i-1,1) )
     ENDIF
     gamma_tab(n) = y/rytoev
  ENDDO

END SUBROUTINE initialize_gamma_tab

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  !----------------------------------------------------------------------------
  ! Calculates the interpolation polynome between 2 points - CG
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp

  IMPLICIT NONE
  REAL(dp) :: t1, t2, f1, f2, df1, df2, poly(4)

  poly(4) = ( (t2-t1) * (df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(3) = (df2-df1)/(2*(t2-t1)) - &
            1.5d0*(t2+t1)*( (t2-t1)*(df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(2) = df1 - 2*t1*poly(3) - 3*t1*t1*poly(4)
  poly(1) = f1 - poly(2)*t1 - poly(3)*t1**2 - poly(4)*t1**3

END SUBROUTINE determine_polycut
