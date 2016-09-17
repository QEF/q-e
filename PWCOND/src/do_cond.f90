!
! Copyright (C) 2003-2009 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE do_cond(done)
!
!   This is the main driver of the pwcond.x program.
!   It calculates the complex band structure, solves the
!   scattering problem and calculates the transmission coefficient.
!
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau, atm
  USE cell_base,  ONLY : at, bg, tpiba
  USE klist,      ONLY : npk, xk, two_fermi_energies
  USE ldaU,       ONLY : lda_plus_U, lda_plus_u_kind, U_projection, &
                         Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha, &
                         Hubbard_J0, Hubbard_beta
  USE spin_orb,   ONLY : lspinorb, domag
  USE uspp,       ONLY : okvan
  USE gvect,      ONLY : ecutrho
  USE gvecw,      ONLY : ecutwfc
  USE symm_base,  ONLY: nsym, s, t_rev, time_reversal
  USE cond
  USE io_files,   ONLY: tmp_dir, prefix
  !!! RECOVER
  USE cond_restart
  USE input_parameters, ONLY: max_seconds
  USE check_stop, ONLY: check_stop_init, check_stop_now
  !!!
  USE noncollin_module, ONLY : noncolin, i_cons
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp_global, ONLY : mp_startup, npool
  USE mp_world,  ONLY : world_comm, nproc
  USE paw_onecenter,      ONLY : PAW_potential
  USE paw_variables,      ONLY : okpaw, ddd_PAW
  USE mp
  USE environment,   ONLY : environment_start

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: outdir
  !
  LOGICAL, INTENT(OUT) :: done
  !
  REAL(DP) :: wtot, tk
  INTEGER :: ik, ien, ipol, ios, nt
  INTEGER :: loop1, loop2, loop1_in, loop1_fin, loop2_in, loop2_fin
  LOGICAL :: lso_l, lso_s, lso_r, skip_equivalence = .FALSE.
  REAL(DP) :: ecutwfc_l, ecutwfc_s, ecutwfc_r
  REAL(DP) :: ecutrho_l, ecutrho_s, ecutrho_r

  !!! RECOVER
  LOGICAL :: tran_save
  !!!

  NAMELIST /inputcond/ outdir, prefixt, prefixl, prefixs, prefixr,     &
                       band_file, tran_file, save_file, fil_loc,       &
                       lwrite_loc, lread_loc, lwrite_cond, lread_cond, &
                       tran_prefix, recover, max_seconds, loop_ek,     &
                       orbj_in,orbj_fin,ikind,iofspin,llocal,          &
                       bdl, bds, bdr, nz1, energy0, denergy, ecut2d,   &
                       start_e, last_e, start_k, last_k,               &
                       ewind, epsproj, delgep, cutplot,                &
                       tk_plot, lorb, lorb3d, lcharge
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PWCOND' )
  CALL start_clock('init')
!
!   set default values for variables in namelist
!
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefixt = ' '
  prefixl = ' '
  prefixs = ' '
  prefixr = ' '
  band_file = ' '
  tran_file = ' '
  save_file = ' '
  fil_loc = ' '
  loop_ek = .FALSE.
  lwrite_loc = .FALSE.
  lread_loc = .FALSE.
  lwrite_cond = .FALSE.
  lread_cond  = .FALSE.
  !!! RECOVER
  tran_prefix = ' '
  recover = .FALSE.
  !!!
  orbj_in = 0
  orbj_fin = 0
  iofspin = 1
  ikind = 0
  bdl = 0.d0
  bds = 0.d0
  bdr = 0.d0
  nz1 = 11
  energy0 = 0.d0
  denergy = 0.d0
  ecut2d = 0.d0
  start_e = 0
  last_e = 0
  start_k = 0
  last_k = 0
  ewind = 1.d0
  llocal = .FALSE.
  epsproj = 1.d-3
  delgep = 5.d-10
  cutplot = 2.d0
  tk_plot = 0
  lorb=.FALSE.
  lorb3d=.FALSE.
  lcharge=.FALSE.

  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     !     reading the namelist inputcond
     !
     READ (5, inputcond, err=200, iostat=ios )
200  CALL errore ('do_cond','reading inputcond namelist',ABS(ios))
     tmp_dir=trimcheck (outdir)
     !
     !     Reading 2D k-point
     READ(5, *, err=250, iostat=ios ) nkpts
250     CALL errore ('do_cond','reading number of kpoints',ABS(ios))
     IF (nkpts>0) THEN
        ALLOCATE( xyk(2,nkpts) )
        ALLOCATE( wkpt(nkpts) )
        wtot = 0.d0
        DO ik = 1, nkpts
           READ(5, *, err=300, iostat=ios) xyk(1,ik), xyk(2,ik), wkpt(ik)
           wtot = wtot + wkpt(ik)
        ENDDO
        DO ik = 1, nkpts
           wkpt(ik) = wkpt(ik)/wtot
        ENDDO
300     CALL errore ('do_cond','2D k-point',ABS(ios))
     ELSE
        ALLOCATE( xyk(2,npk) )
        ALLOCATE( wkpt(npk) )
        READ(5, *, err=350, iostat=ios) nk1ts, nk2ts, k1ts, k2ts
350     CALL errore ('do_cond','2D k-point size or shift',ABS(ios))
     ENDIF

     !
     !     To form the array of energies for calculation
     !
     READ(5, *, err=400, iostat=ios ) nenergy
     ALLOCATE( earr(nenergy) )
     ALLOCATE( tran_tot(nenergy) )
     IF(ABS(denergy).LE.1.d-8) THEN
        !   the list of energies is read
        DO ien = 1, nenergy
           READ(5, *, err=400, iostat=ios) earr(ien)
           tran_tot(ien) = 0.d0
        ENDDO
     ELSE
        !   the array of energies is automatically formed
        DO ien = 1, nenergy
           earr(ien) = energy0 + (ien-1)*denergy
           tran_tot(ien) = 0.d0
        ENDDO
     ENDIF

     IF (start_e .GT. 0) THEN
        IF (start_e .GT. last_e .OR. start_e .GT. nenergy) &
           CALL errore('do_cond','wrong value of start_e',1)
        IF (last_e .GT. nenergy) last_e = nenergy
     ELSE
        start_e = 1
        last_e = nenergy
     ENDIF

400  CALL errore ('do_cond','reading energy list',ABS(ios))
     !
  END IF

#if defined(__MPI)
   IF (npool > 1) CALL errore('pwcond','pools not implemented',npool)
   ik = IAND ( nproc, nproc-1 )
   IF ( nproc /= 1 .AND. ik /= 0 ) &
       CALL errore('pwcond','you should use 2^N number of CPUs',1)
#endif

!-- Some check and initialization for plotting the scattering states
  IF ( lorb .AND. ikind == 2 ) &
       CALL errore('do_cond','lorb not working with ikind = 2',1)
  IF (lorb3d) lorb = .TRUE.
  IF (lcharge) lorb = .TRUE.
!--

!
! ... Broadcast variables
!
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefixt, ionode_id, world_comm )
  CALL mp_bcast( prefixl, ionode_id, world_comm )
  CALL mp_bcast( prefixs, ionode_id, world_comm )
  CALL mp_bcast( prefixr, ionode_id, world_comm )
  CALL mp_bcast( band_file, ionode_id, world_comm )
  CALL mp_bcast( tran_file, ionode_id, world_comm )
  CALL mp_bcast( fil_loc, ionode_id, world_comm )
  CALL mp_bcast( save_file, ionode_id, world_comm )
  CALL mp_bcast( loop_ek, ionode_id, world_comm )
  CALL mp_bcast( lwrite_loc, ionode_id, world_comm )
  CALL mp_bcast( lread_loc, ionode_id, world_comm )
  CALL mp_bcast( lwrite_cond, ionode_id, world_comm )
  CALL mp_bcast( lread_cond, ionode_id, world_comm )
  !!! RECOVER
  CALL mp_bcast( tran_prefix, ionode_id, world_comm )
  CALL mp_bcast( max_seconds, ionode_id, world_comm )
  CALL mp_bcast( recover, ionode_id, world_comm )
  !!!
  CALL mp_bcast( ikind, ionode_id, world_comm )
  CALL mp_bcast( iofspin, ionode_id, world_comm )
  CALL mp_bcast( orbj_in, ionode_id, world_comm )
  CALL mp_bcast( orbj_fin, ionode_id, world_comm )
  CALL mp_bcast( llocal, ionode_id, world_comm )
  CALL mp_bcast( tk_plot, ionode_id, world_comm )
  CALL mp_bcast( lorb, ionode_id, world_comm )
  CALL mp_bcast( lorb3d, ionode_id, world_comm )
  CALL mp_bcast( lcharge, ionode_id, world_comm )
  CALL mp_bcast( bdl, ionode_id, world_comm )
  CALL mp_bcast( bds, ionode_id, world_comm )
  CALL mp_bcast( bdr, ionode_id, world_comm )
  CALL mp_bcast( nz1, ionode_id, world_comm )
  CALL mp_bcast( energy0, ionode_id, world_comm )
  CALL mp_bcast( denergy, ionode_id, world_comm )
  CALL mp_bcast( ecut2d, ionode_id, world_comm )
  CALL mp_bcast( start_e, ionode_id, world_comm )
  CALL mp_bcast( last_e, ionode_id, world_comm )
  CALL mp_bcast( ewind, ionode_id, world_comm )
  CALL mp_bcast( epsproj, ionode_id, world_comm )
  CALL mp_bcast( delgep, ionode_id, world_comm )
  CALL mp_bcast( cutplot, ionode_id, world_comm )
  CALL mp_bcast( nkpts, ionode_id, world_comm )
  CALL mp_bcast( nenergy, ionode_id, world_comm )
  CALL mp_bcast( nk1ts, ionode_id, world_comm )
  CALL mp_bcast( nk2ts, ionode_id, world_comm )
  CALL mp_bcast( k1ts, ionode_id, world_comm )
  CALL mp_bcast( k2ts, ionode_id, world_comm )

  IF ( .NOT. ionode ) THEN
     IF (nkpts>0) THEN
        ALLOCATE( xyk(2,nkpts) )
        ALLOCATE( wkpt(nkpts) )
     ELSE
        ALLOCATE( xyk(2,npk) )
        ALLOCATE( wkpt(npk) )
     ENDIF
     ALLOCATE( earr(nenergy) )
     ALLOCATE( tran_tot(nenergy) )
  ENDIF
  IF (nkpts>0) THEN
     CALL mp_bcast( xyk, ionode_id, world_comm )
     CALL mp_bcast( wkpt, ionode_id, world_comm )
  ENDIF
  CALL mp_bcast( earr, ionode_id, world_comm )
  CALL mp_bcast( tran_tot, ionode_id, world_comm )


!
! Now allocate space for pwscf variables, read and check them.
!

IF (lread_cond) THEN
  call save_cond (.false.,1,efl,nrzl,nocrosl,noinsl,   &
                  norbl,rl,rabl,betarl)
  if(ikind.eq.1) then
    call save_cond (.false.,2,efs,nrzs,ik,      &
                             noinss,norbs,rs,rabs,betars)
    norbr = norbl
    nocrosr = nocrosl
    noinsr = noinsl
  endif
  if(ikind.eq.2) then
    call save_cond (.false.,2,efs,nrzs,ik,      &
                             noinss,norbs,rs,rabs,betars)

    call save_cond (.false.,3,efr,nrzr,nocrosr,&
                             noinsr,norbr,rr,rabr,betarr)
  endif
ELSE
  lso_l=.false.
  lso_s=.false.
  lso_r=.false.
  ecutwfc_l=0.0_DP
  ecutwfc_s=0.0_DP
  ecutwfc_r=0.0_DP
  ecutrho_l=0.0_DP
  ecutrho_s=0.0_DP
  ecutrho_r=0.0_DP
  IF (prefixt.ne.' ') then
    prefix = prefixt

    call read_file
    IF (ikind.eq.0) then
      CALL init_cond(1,'t')
    ELSEIF (ikind.eq.1) then
      CALL init_cond(2,'t')
    ELSEIF (ikind.eq.2) then
      CALL init_cond(3,'t')
    ENDIF
    CALL clean_pw(.true.)
  ENDIF
  IF (prefixl.ne.' ') then
    prefix = prefixl
    call read_file
    lso_l=lspinorb
    ecutwfc_l=ecutwfc
    ecutrho_l=ecutrho
    CALL init_cond(1,'l')
  ENDIF

  IF (prefixr.ne.' ') then
    CALL clean_pw(.true.)
    prefix = prefixr
    call read_file
    lso_r=lspinorb
    ecutwfc_r=ecutwfc
    ecutrho_r=ecutrho
    CALL init_cond(1,'r')
  ENDIF
  IF (prefixs.ne.' ') then
    call clean_pw(.true.)
    prefix = prefixs
    call read_file
    lso_s=lspinorb
    ecutwfc_s=ecutwfc
    ecutrho_s=ecutrho
    CALL init_cond(1,'s')
  ENDIF

  IF (two_fermi_energies.or.i_cons /= 0) &
     CALL errore('pwcond',&
     'The pwcond code with constrained magnetization is not yet available',1)
  IF (ikind==1.and.(lso_l.neqv.lso_s)) &
     CALL errore('pwcond',&
     'Spin-orbit flag in left lead and scattering region do not match',1)
  IF (ikind==2.and.((lso_l.neqv.lso_s).or.(lso_r.neqv.lso_s))) &
     CALL errore('pwcond',&
     'Spin-orbit flag in left, right lead and scattering region do not match',1)
  IF (ikind>0.and.((ecutwfc_l.ne.ecutwfc_s).or.(ecutrho_l.ne.ecutrho_s)))  &
     CALL errore('do_cond',&
            'different cut-offs on left lead and scattering region',1)
  IF ((ecutwfc_r>0.0_DP)) THEN
     IF ((ecutwfc_r.ne.ecutwfc_s).or.(ecutrho_r.ne.ecutrho_s))  &
     CALL errore('do_cond',&
        'different cut-offs on right lead and scattering region',1)
  ENDIF
ENDIF

IF (lwrite_cond) then
  call save_cond (.true.,1,efl,nrzl,nocrosl,noinsl,         &
                  norbl,rl,rabl,betarl)
  if(ikind.gt.0) call save_cond (.true.,2,efs,nrzs,-1,      &
                             noinss,norbs,rs,rabs,betars)
  if(ikind.gt.1) call save_cond (.true.,3,efr,nrzr,nocrosr,&
                             noinsr,norbr,rr,rabr,betarr)
  write(stdout,*) 'information needed for PWCOND has been written in file'
  CALL stop_clock('init')
  CALL stop_clock('PWCOND')
  CALL print_clock_pwcond()
  return
endif

IF (nkpts==0) THEN
   time_reversal = .NOT. (noncolin .AND. domag)
   IF (ionode) THEN
      CALL kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, bg, &
                        npk, k1ts, k2ts, 0, nk1ts, nk2ts, 1, nkpts, xk, wkpt )
      call cryst_to_cart(nkpts,xk,at,-1)
      DO ik=1,nkpts
         xyk(1,ik)=xk(1,ik)
         xyk(2,ik)=xk(2,ik)
      ENDDO
   ENDIF
   CALL mp_bcast( nkpts, ionode_id, world_comm )
   CALL mp_bcast( xyk, ionode_id, world_comm )
   CALL mp_bcast( wkpt, ionode_id, world_comm )
ELSE
   tk_plot = 0
ENDIF

if (tk_plot.lt.0) CALL errore('do_cond','tk_plot should be > 0',1)
If (tk_plot.gt.0) loop_ek = .TRUE.
IF (ikind.ne.0.and.tk_plot.gt.0) ALLOCATE( tran_k(npk) )

IF (start_k .GT. 0) THEN
   IF (start_k .GT. last_k .OR. start_k .GT. nkpts) &
     CALL errore('do_cond','wrong value of start_k',1)
   IF (last_k .GT. nkpts) last_k = nkpts
ELSE
   start_k = 1
   last_k = nkpts
ENDIF
CALL mp_bcast( start_k, ionode_id, world_comm )
CALL mp_bcast( last_k, ionode_id, world_comm )

  !!! RECOVER
  ! Simple restart mechanism for transmission calculations
  ! (tran_prefix must be specified on input in order to enable restart)
  !!!
  ! Initialization of recover:
  IF (ikind.NE.0 .AND. tran_prefix.NE.' ') THEN
     !
     tran_save = .TRUE.
     CALL check_stop_init ()
     ! if recover flag is set to true, then check info file
     IF ( recover ) THEN
        ! read and check info file
        ! (lists of energies and k-points read from info file
        ! must coindice with those built from input parameters)
        CALL cond_readfile( 'init', ios )
     ELSE
        ! create restart directory and write info file
        CALL cond_writefile( 'init' )
     ENDIF

  ELSE
     !
     tran_save = .FALSE.
     IF (recover)  call errore('do_cond',&
        'you must specify tran_prefix in order to restart',1)
  ENDIF
  !!!

  CALL cond_out

  CALL stop_clock('init')

  IF (llocal) &
  CALL local_set(nocrosl,noinsl,norbl,noinss,norbs,nocrosr,noinsr,norbr)

!-- 
! Set up 2 loops over energies and over k-points 
  if (loop_ek) then
    loop1_in  = start_e
    loop1_fin = last_e
    loop2_in  = start_k
    loop2_fin = last_k
  else
    loop1_in  = start_k
    loop1_fin = last_k
    loop2_in  = start_e
    loop2_fin = last_e
  endif
!--

  DO loop1 = loop1_in, loop1_fin

    if (.not.loop_ek) then
      CALL init_gper(loop1)
      CALL local(1)
    endif

    DO loop2 = loop2_in, loop2_fin

      if (loop_ek) then
        ien = loop1
        ik  = loop2
      else
        ik  = loop1
        ien = loop2
      endif

!      write(6,*) loop1_in, loop1_fin, loop2_in, loop2_fin, loop1, loop2 
      WRITE(stdout,'("---  E-Ef = ",f12.7, "  k = ",2f12.7)') &
            earr(ien), (xyk (ipol, ik) , ipol = 1, 2)
      WRITE(stdout,'("---  ie = ",i10, "  ik = ",i10)') &
            ien, ik

      !!! RECOVER
      ! if recover mechanism is enabled
      IF (recover .AND. ikind.NE.0) THEN
         !
         WRITE(stdout,'(A)') 'Reading transmission from restart file:'
         ! check if the transmission has already been calculated for
         ! this specific k-point and energy value
         CALL cond_readfile( 'tran', ios, ik, ien, tk )
         ! if so, add it to the total transmission with the correct weight
         ! and skip to the next energy in the list
         IF ( ios .EQ. 0 ) THEN
            WRITE(stdout,'(a24, 2f12.7,/)') 'E-Ef(ev), T = ',earr(ien),tk
            tran_tot(ien) = tran_tot(ien) + wkpt(ik)*tk
            !CALL mp_bcast( tran_tot(ien), ionode_id, world_comm )
            CYCLE
         ! if not, do the actual calculation
         ELSE
            IF ( ios .EQ. 1 ) THEN
               write(stdout,'(" File not found...")')
            ELSE
               write(stdout,'(" FAILED: could not read from file...")')
            ENDIF
            write(stdout,'(" ... computing transmission",/)')
         ENDIF
      ENDIF
      !!!

      if (loop_ek) then
        CALL init_gper(ik)
        CALL local(ien)
      endif

      eryd = earr(ien)/rytoev + efl
      CALL form_zk(n2d, nrzpl, zkrl, zkl, eryd, tpiba)
      CALL scatter_forw(nrzl, nrzpl, zl, psiperl, zkl, norbl,     &
                        tblml, crosl, taunewl, rl, rabl, betarl)
      CALL compbs(1, nocrosl, norbl, nchanl, kvall, kfunl, kfundl, &
                  kintl, kcoefl, ik, ien)

      IF (ikind.EQ.2) THEN
        eryd = earr(ien)/rytoev + efr
        CALL form_zk(n2d, nrzpr, zkrr, zkr, eryd, tpiba)
        CALL scatter_forw(nrzr, nrzpr, zr, psiperr, zkr, norbr,    &
                          tblmr, crosr, taunewr, rr, rabr, betarr)
        CALL compbs(0, nocrosr, norbr, nchanr, kvalr, kfunr, kfundr,&
                     kintr, kcoefr, ik, ien)
      ENDIF

      CALL summary_band(ik,ien)

      IF (ikind.NE.0) THEN
         eryd = earr(ien)/rytoev + efs
         CALL form_zk(n2d, nrzps, zkrs, zks, eryd, tpiba)
         CALL scatter_forw(nrzs, nrzps, zs, psipers, zks, norbs,    &
                          tblms, cross, taunews, rs, rabs, betars)

         WRITE(stdout,*) 'to transmit'
         CALL transmit(ik, ien, tk, .true.)
         !
         ! To add T(k) to the total T
         !
         tran_tot(ien) = tran_tot(ien) + wkpt(ik)*tk
         if (tk_plot.gt.0) tran_k(ik) = tk
         if (lorb) CALL transmit(ik, ien, tk, .false.)

         !
         !!! RECOVER
         ! if recover is enabled, save the partial transmission on file,
         ! and then check for stopping condition
         IF ( tran_save ) THEN
            CALL cond_writefile( 'tran', ik, ien, tk )
            IF ( check_stop_now() ) THEN
               CALL free_mem
               CALL stop_clock('PWCOND')
               CALL print_clock_pwcond()
               done = .FALSE.
               RETURN
            ENDIF
         ENDIF
         !!!

      ENDIF

      if (loop_ek) CALL free_mem

    ENDDO

    if (ikind.ne.0.and.tk_plot.gt.0.and.ionode) &
             CALL summary_tran_k(ien,nk1ts,nk2ts,k1ts,k2ts)

    if (.not.loop_ek) CALL free_mem

  ENDDO

  IF(ikind.ne.0.and.ionode) CALL summary_tran_tot()

  CALL stop_clock('PWCOND')
  CALL print_clock_pwcond()

  done = .TRUE.

  DEALLOCATE( xyk )
  DEALLOCATE( wkpt )
  DEALLOCATE( earr )
  DEALLOCATE( tran_tot )
  IF (ikind.ne.0.and.tk_plot.gt.0) DEALLOCATE( tran_k )

  RETURN

END SUBROUTINE do_cond


