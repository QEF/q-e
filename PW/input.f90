!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine iosys
  !-----------------------------------------------------------------------
  !   this subroutine reads input data from standard input (unit 5)
  !     ------------------------------------------------------------------
  use pwcom, only: dt, iswitch, ltaucry, restart, lscf, iprint, &
       time_max, iverbosity, order, reduce_io, modenum, ibrav, &
       celldm, nat, input_drho, output_drho, &
       ntyp, nbnd, nelec, nr1,nr2,nr3, nr1s,nr2s,nr3s, lstres, &
       degauss, ngauss, nosym, ltetra, ecfixed, qcutz, q2sigma, &
       ecutwfc, lsda, nspin, dual, lxkcry, noinv, starting_magnetization, &
       lda_plus_U, Hubbard_U, Hubbard_alpha, niter_with_fixed_ns, rytoev, &
       niter, tr2, ethr, mixing_beta, nmix,&
       isolve, max_cg_iter, david, loverlap, diis_buff, diis_wfc_keep, &
       diis_start_dav, startingwfc, startingpot, startingconfig, &
       nstep, epse, epsf, amass, &
       temperature, lforce, ttol, delta_t, nraise, ntcheck, upscale, &
       press, cmass, calc, cell_factor, xqq, alat, ntypx, &
       lmovecell, imix, at, omega, ityp, tau, nks, xk, wk, uakbar, amconv, &
       force, at_old, omega_old, starting_scf_threshold, title, crystal,  &
       atm, nk1, nk2, nk3, k1, k2, k3
  use io, only : tmp_dir, prefix, pseudo_dir, pseudop
  use constants, only: pi
#ifdef PARA
  use para, only: me
  use mp
#endif
  !
  implicit none
  !
  !
  ! local variables
  !
  character (len=30) :: atomic_positions
  integer :: unit = 5, ionode_id = 0, i, ia, ios, is
  !
  ! CONTROL namelist

  character(len=256) :: outdir
  logical :: tstress, tprnfor
  real (kind=8) :: max_seconds
  real (kind=8) :: ekin_conv_thr, etot_conv_thr, forc_conv_thr
  character(len=80) :: restart_mode, disk_io, calculation, verbosity
  integer :: isave, ndr, ndw
  NAMELIST / control / title, calculation, verbosity, &
       restart_mode, nstep, iprint, isave, tstress, tprnfor, &
       dt, ndr, ndw, outdir, prefix, max_seconds, ekin_conv_thr,&
       etot_conv_thr, forc_conv_thr, pseudo_dir, disk_io

  ! SYSTEM namelist

  character(len=80) :: occupations, smearing, xc_type
  integer :: nelup, neldw, nr1b, nr2b, nr3b
  real (kind=8) :: ecutrho
  NAMELIST / system / ibrav, celldm, nat, ntyp, nbnd, nelec, &
       ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s, nr3s, &
       nr1b, nr2b, nr3b, nosym, starting_magnetization, &
       occupations, degauss, smearing, &
       nelup, neldw, nspin, ecfixed, qcutz, q2sigma, xc_type, &
       lda_plus_U, Hubbard_U, Hubbard_alpha

  ! ELECTRONS namelist

  character(len=80) :: orthogonalization, &
       electron_dynamics, electron_velocities, electron_temperature
  integer :: electron_maxstep
  real(kind=8) :: emass_cutoff, electron_damping, fnosee, &
       ortho_eps, ortho_max, emass, ekincw, ampre, grease
  integer :: empty_states_nbnd, empty_states_maxstep
  real(kind=8) :: empty_states_delt, empty_states_emass, empty_states_ethr
  integer :: diis_size, diis_nreset, diis_maxstep
  logical :: diis_rot, diis_chguess, twall
  real(kind=8) :: diis_hcut, diis_wthr, diis_delt
  real(kind=8) :: diis_fthr, diis_temp, diis_achmix, diis_g0chmix
  integer :: diis_nchmix, diis_nrot(3)
  real(kind=8) :: diis_g1chmix, diis_rothr(3), diis_ethr
  character(len=80) :: mixing_mode
  integer :: mixing_ndim, mixing_fixed_ns
  real (kind=8) :: conv_thr
  character(len=80) :: diagonalization
  integer :: diago_cg_maxiter, diago_david_ndim, diago_diis_buff, &
       diago_diis_start
  logical :: diago_diis_keep

  NAMELIST / electrons / emass, emass_cutoff, orthogonalization, &
       electron_maxstep, ortho_eps, ortho_max, electron_dynamics, &
       electron_damping, electron_velocities, electron_temperature,&
       ekincw, fnosee, ampre, grease, twall, &
       empty_states_nbnd, empty_states_maxstep, empty_states_delt, &
       empty_states_emass, empty_states_ethr, &
       diis_size, diis_nreset, diis_hcut, diis_wthr, diis_delt, &
       diis_maxstep, diis_rot, diis_fthr, diis_temp, diis_achmix, &
       diis_g0chmix, diis_g1chmix, diis_nchmix, diis_nrot, diis_rothr, &
       diis_ethr, diis_chguess, &
       mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns, &
       diago_cg_maxiter, diago_david_ndim, diago_diis_buff, &
       diago_diis_start, diago_diis_keep, diagonalization, &
       startingwfc, startingpot, conv_thr

  ! IONS namelist

  character(len=80) :: ion_dynamics, ion_positions, ion_velocities, &
       ion_temperature, potential_extrapolation
  integer :: ion_nstepe
  integer :: ion_maxstep
  real(kind=8) :: ion_radius(ntypx), ion_damping, fnosep, tempw, &
       amprp(ntypx), greasp, tolp
  logical :: tranp(ntypx)

  NAMELIST / ions / ion_dynamics, ion_radius, ion_damping, ion_positions, &
       ion_velocities, ion_temperature, &
       tempw, fnosep, tranp, amprp, greasp, tolp, &
       ion_nstepe, ion_maxstep, upscale, potential_extrapolation

  ! CELL namelist

  character(len=80) :: cell_parameters, cell_dynamics, cell_velocities, &
       cell_temperature, cell_dofree
  real(kind=8) :: cell_damping, fnoseh, wmass, temph, greash

  NAMELIST / cell / cell_parameters, cell_dynamics, cell_velocities, press, &
       wmass, cell_temperature, temph, fnoseh, cell_dofree, greash, &
       cell_factor

  ! PHONON namelist

  NAMELIST / phonon / modenum, xqq

  ! ...   Variables initialization for CONTROL
  !
  title = ' '
  calculation='scf'
  verbosity = 'default'
  max_seconds  = 1.d+6
  restart_mode = 'from_scratch'
  nstep  = 50
  iprint = 100000
  isave  = 100
  tstress = .FALSE.
  tprnfor = .FALSE.
  dt    = 1.0d0
  ndr = 50
  ndw = 50
  outdir = './'      ! use the path specified as Outdir and
  prefix = 'pwscf'   ! the filename prefix to store the output
  ekin_conv_thr = 1.d-6
  etot_conv_thr = 1.d-4
  forc_conv_thr = 1.d-3
  disk_io = 'default'
  noinv = .false.    ! not actually used
  !
#ifdef T3D
  call pxfgetenv('HOME',0,pseudo_dir,i,ios)
#else
  call getenv('HOME',pseudo_dir)
#endif
  pseudo_dir=trim(pseudo_dir)//'/pw/pseudo/'

  ! ...   Variables initialization for SYSTEM

  ibrav  =-1
  celldm = (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
  nat    = 0
  ntyp   = 0
  nbnd   = 0
  nelec  = 0
  ecutwfc= 0.d0
  ecutrho= 0.d0
  nr1  = 0
  nr2  = 0
  nr3  = 0
  nr1s = 0
  nr2s = 0
  nr3s = 0
  nr1b = 0
  nr2b = 0
  nr3b = 0
  occupations = 'fixed'
  smearing = 'gaussian'
  degauss = 0.d0
  nelup = 0
  neldw = 0
  nspin = 1
  nosym = .false.
  ecfixed = 0.d0
  qcutz   = 0.d0
  q2sigma = 0.01d0
  xc_type = 'PZ'
  lda_plus_U = .false.
  Hubbard_U(:) = 0.d0
  Hubbard_alpha(:) = 0.d0
  starting_magnetization(:) = 0.d0

  ! ...   Variables initialization for ELECTRONS
  !
  electron_maxstep = 100
  emass = 400.d0
  emass_cutoff = 2.5d0
  orthogonalization = 'ortho'
  ortho_eps = 1.d-8
  ortho_max = 20
  electron_dynamics = 'none'
  electron_damping = 0.1d0
  electron_velocities = 'default' ! ( 'zero' | 'default' )
  electron_temperature = 'not_controlled'
  ! ( 'nose' | 'not_controlled' | 'rescaling')
  ekincw = 0.001d0
  fnosee = 1.0d0
  ampre  = 0.0d0
  grease = 1.0d0
  twall  = .FALSE.
  empty_states_nbnd = 0
  empty_states_maxstep = 100
  empty_states_delt = 0.0d0
  empty_states_emass = 0.0d0
  empty_states_ethr = 0.0d0
  diis_size = 4
  diis_nreset = 3
  diis_hcut = 1.0d0
  diis_wthr = 0.0d0
  diis_delt = 0.0d0
  diis_maxstep = 0
  diis_rot = .FALSE.
  diis_fthr = 0.0d0
  diis_temp = 0.0d0
  diis_achmix = 0.0d0
  diis_g0chmix = 0.0d0
  diis_g1chmix = 0.0d0
  diis_nchmix = 3
  diis_nrot = 3
  diis_rothr  = 0.0d0
  diis_ethr   = 0.0d0
  diis_chguess = .FALSE.
  mixing_mode ='plain'
  mixing_fixed_ns = 0
  mixing_beta = 0.7
  mixing_ndim = 8

  diago_cg_maxiter = 20
  diago_david_ndim  =4
  diago_diis_buff = 200
  diago_diis_keep = .false.
  diago_diis_start= 0
  startingwfc = 'atomic'
  startingpot = 'atomic'
  conv_thr = 1.d-6

  ! ...   Variables initialization for IONS
  !
  ion_dynamics = 'none'  ! ( 'sd' | 'cg' | 'damp' | 'verlet' | 'beeman' |'bfgs' | 'none' )
  ion_radius = 0.5d0
  ion_damping = 0.1
  ion_positions = 'default' ! ( 'default' | 'from_input' )
  ion_velocities = 'default'
  ! ( 'zero' | 'default' | 'random' | 'from_input' )
  ion_temperature = 'not_controlled'
  ! ( 'nose' | 'not_controlled' | 'rescaling' )
  tempw = 300.0d0
  fnosep = 1.0d0
  tranp(:) = .FALSE.
  amprp(:) = 0.0d0
  greasp = 1.0d0
  tolp = 100.d0
  ion_nstepe = 1
  ion_maxstep = 100
  upscale = 10
  potential_extrapolation='default'
  !
  ! ...   Variables initialization for CELL
  !
  cell_parameters = 'default'
  cell_dynamics = 'none'      ! ( 'sd' | 'md' | 'damp' | 'md-w' | 'damp-w' | 'none' )
  cell_velocities = 'default' ! ( 'zero' | 'default' )
  press = 0.0d0
  wmass = 0.0d0
  cell_temperature = 'not_controlled' ! ( 'nose' | 'not_controlled' | 'rescaling' )
  temph = 0.0d0
  fnoseh = 0.0d0
  greash = 1.0d0
  cell_dofree = 'all'
  ! ('all'* | 'volume' | 'x' | 'y' | 'z' | 'xy' | 'xz' | 'yz' | 'xyz' )
  nraise = 100
  delta_t = 1.d0
  !
  ! ...   Variables initialization for PHONON
  !
  modenum = 0
  xqq = 0.d0
  !
#ifdef PARA
  if (me == 1)  then
#endif
     ios = 0
     READ (unit, control, iostat = ios )
     if (ios /= 0) call error ('reading','namelist &control',1)
     !
     ! reset default values for ion_dynamics according to definition
     ! of calculation in &control
     !
     if (TRIM(calculation) == 'relax' ) then
        ion_dynamics='bfgs'
     end if
     if (TRIM(calculation) == 'md' ) then
        ion_dynamics='verlet'
     end if
     if (TRIM(calculation) == 'vc-relax' ) then
        ion_dynamics ='damp'
     end if
     if (TRIM(calculation) == 'vc-md' ) then
        ion_dynamics ='beeman'
     end if
     ! reset defaulf value of startingpot according to calculation value

     if (TRIM(calculation) == 'nscf' .or. TRIM(calculation) == 'phonon' ) then
        startingpot = 'file'
     else
        startingpot = 'atomic'
     end if

     READ (unit, system, iostat = ios )
     if (ios /= 0) call error ('reading','namelist &system',2)
     READ (unit, electrons, iostat = ios )
     if (ios /= 0) call error ('reading','namelist &electrons',3)
     if ( TRIM(calculation) == 'relax'   .or. TRIM(calculation) == 'md'    .or.&
          TRIM(calculation) == 'vc-relax'.or. TRIM(calculation) == 'vc-md' .or.&
          TRIM(calculation) == 'cpmd' .or. TRIM(calculation) == 'vc-cpmd' ) then
        READ (unit, ions, iostat = ios )
        if (ios /= 0) call error ('reading','namelist &ions',4)
     end if
     if ( TRIM(calculation) == 'vc-relax'.or. TRIM(calculation) == 'vc-md' .or.&
          TRIM(calculation) == 'vc-cpmd' ) then
        READ (unit, cell, iostat = ios )
        if (ios /= 0) call error ('reading','namelist &cell',5)
     end if
     if (TRIM(calculation) == 'phonon' ) then
        READ (unit, phonon, iostat = ios )
        if (ios /= 0) call error ('reading','namelist &phonon',5)
     end if
#ifdef PARA
  end if

  !
  ! ...   CONTROL Variables Broadcast
  !
  CALL mp_bcast( title, ionode_id )
  CALL mp_bcast( calculation, ionode_id )
  CALL mp_bcast( verbosity, ionode_id )
  CALL mp_bcast( restart_mode, ionode_id )
  CALL mp_bcast( nstep, ionode_id )
  CALL mp_bcast( iprint, ionode_id )
  CALL mp_bcast( isave, ionode_id )
  CALL mp_bcast( tstress, ionode_id )
  CALL mp_bcast( tprnfor, ionode_id )
  CALL mp_bcast( dt, ionode_id )
  CALL mp_bcast( ndr, ionode_id )
  CALL mp_bcast( ndw, ionode_id )
  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( max_seconds, ionode_id )
  CALL mp_bcast( ekin_conv_thr, ionode_id )
  CALL mp_bcast( etot_conv_thr, ionode_id )
  CALL mp_bcast( forc_conv_thr, ionode_id )
  CALL mp_bcast( pseudo_dir, ionode_id )
  CALL mp_bcast( disk_io, ionode_id )
  !
  ! ...   SYSTEM Variables Broadcast
  !
  CALL mp_bcast( ibrav, ionode_id  )
  CALL mp_bcast( celldm, ionode_id  )
  CALL mp_bcast( nat, ionode_id  )
  CALL mp_bcast( ntyp, ionode_id  )
  CALL mp_bcast( nbnd, ionode_id  )
  CALL mp_bcast( nelec, ionode_id  )
  CALL mp_bcast( ecutwfc, ionode_id  )
  CALL mp_bcast( ecutrho, ionode_id  )
  CALL mp_bcast( nr1, ionode_id  )
  CALL mp_bcast( nr2, ionode_id  )
  CALL mp_bcast( nr3, ionode_id  )
  CALL mp_bcast( nr1s, ionode_id  )
  CALL mp_bcast( nr2s, ionode_id  )
  CALL mp_bcast( nr3s, ionode_id  )
  CALL mp_bcast( nr1b, ionode_id  )
  CALL mp_bcast( nr2b, ionode_id  )
  CALL mp_bcast( nr3b, ionode_id  )
  CALL mp_bcast( occupations, ionode_id  )
  CALL mp_bcast( smearing, ionode_id  )
  CALL mp_bcast( degauss, ionode_id  )
  CALL mp_bcast( nelup, ionode_id )
  CALL mp_bcast( neldw, ionode_id )
  CALL mp_bcast( nspin, ionode_id )
  CALL mp_bcast( nosym, ionode_id )
  CALL mp_bcast( ecfixed, ionode_id )
  CALL mp_bcast( qcutz, ionode_id )
  CALL mp_bcast( q2sigma, ionode_id )
  CALL mp_bcast( xc_type, ionode_id )
  CALL mp_bcast( lda_plus_U, ionode_id )
  CALL mp_bcast( Hubbard_U, ionode_id )
  CALL mp_bcast( Hubbard_alpha, ionode_id )
  CALL mp_bcast( starting_magnetization, ionode_id )

  ! ...   ELECTRONS Variables Broadcast
  !
  CALL mp_bcast( emass, ionode_id )
  CALL mp_bcast( emass_cutoff, ionode_id )
  CALL mp_bcast( orthogonalization, ionode_id )
  CALL mp_bcast( electron_maxstep, ionode_id )
  CALL mp_bcast( ortho_eps, ionode_id )
  CALL mp_bcast( ortho_max, ionode_id )
  CALL mp_bcast( electron_dynamics, ionode_id )
  CALL mp_bcast( electron_damping, ionode_id )
  CALL mp_bcast( electron_velocities, ionode_id )
  CALL mp_bcast( electron_temperature, ionode_id )
  CALL mp_bcast( ekincw, ionode_id )
  CALL mp_bcast( fnosee, ionode_id )
  CALL mp_bcast( ampre, ionode_id )
  CALL mp_bcast( grease, ionode_id )
  CALL mp_bcast( twall, ionode_id )
  CALL mp_bcast( empty_states_nbnd, ionode_id )
  CALL mp_bcast( empty_states_maxstep, ionode_id )
  CALL mp_bcast( empty_states_delt, ionode_id )
  CALL mp_bcast( empty_states_emass, ionode_id )
  CALL mp_bcast( empty_states_ethr, ionode_id )
  CALL mp_bcast( diis_size, ionode_id )
  CALL mp_bcast( diis_nreset, ionode_id )
  CALL mp_bcast( diis_hcut, ionode_id )
  CALL mp_bcast( diis_wthr, ionode_id )
  CALL mp_bcast( diis_delt, ionode_id )
  CALL mp_bcast( diis_maxstep, ionode_id )
  CALL mp_bcast( diis_rot, ionode_id )
  CALL mp_bcast( diis_fthr, ionode_id )
  CALL mp_bcast( diis_temp, ionode_id )
  CALL mp_bcast( diis_achmix, ionode_id )
  CALL mp_bcast( diis_g0chmix, ionode_id )
  CALL mp_bcast( diis_g1chmix, ionode_id )
  CALL mp_bcast( diis_nchmix, ionode_id )
  CALL mp_bcast( diis_nrot, ionode_id )
  CALL mp_bcast( diis_rothr, ionode_id )
  CALL mp_bcast( diis_ethr, ionode_id )
  CALL mp_bcast( diis_chguess, ionode_id )
  CALL mp_bcast( mixing_mode, ionode_id )
  CALL mp_bcast( mixing_beta, ionode_id )
  CALL mp_bcast( mixing_ndim, ionode_id )
  CALL mp_bcast( mixing_fixed_ns, ionode_id )
  CALL mp_bcast( diagonalization, ionode_id )
  CALL mp_bcast( diago_cg_maxiter, ionode_id )
  CALL mp_bcast( diago_david_ndim, ionode_id )
  CALL mp_bcast( diago_diis_buff, ionode_id )
  CALL mp_bcast( diago_diis_keep, ionode_id )
  CALL mp_bcast( diago_diis_start,ionode_id )
  CALL mp_bcast( startingwfc, ionode_id )
  CALL mp_bcast( startingpot, ionode_id )
  CALL mp_bcast( conv_thr, ionode_id )

  ! ...   IONS Variables Broadcast
  !
  CALL mp_bcast( ion_dynamics, ionode_id )
  CALL mp_bcast( ion_radius, ionode_id )
  CALL mp_bcast( ion_damping, ionode_id )
  CALL mp_bcast( ion_positions, ionode_id )
  CALL mp_bcast( ion_velocities, ionode_id )
  CALL mp_bcast( ion_temperature, ionode_id )
  CALL mp_bcast( tempw, ionode_id )
  CALL mp_bcast( fnosep, ionode_id )
  CALL mp_bcast( tranp, ionode_id )
  CALL mp_bcast( amprp, ionode_id )
  CALL mp_bcast( greasp, ionode_id )
  CALL mp_bcast( tolp, ionode_id )
  CALL mp_bcast( ion_nstepe, ionode_id )
  CALL mp_bcast( ion_maxstep, ionode_id )
  CALL mp_bcast( upscale, ionode_id )
  CALL mp_bcast( potential_extrapolation, ionode_id )

  ! ...   CELL Variables Broadcast
  !
  CALL mp_bcast( cell_parameters, ionode_id )
  CALL mp_bcast( cell_dynamics, ionode_id )
  CALL mp_bcast( cell_velocities, ionode_id )
  CALL mp_bcast( cell_dofree, ionode_id )
  CALL mp_bcast( press, ionode_id )
  CALL mp_bcast( wmass, ionode_id )
  CALL mp_bcast( cell_temperature, ionode_id )
  CALL mp_bcast( temph, ionode_id )
  CALL mp_bcast( fnoseh, ionode_id )
  CALL mp_bcast( cell_factor, ionode_id )

  ! ...   PHONON Variables Broadcast
  !
  CALL mp_bcast( modenum, ionode_id )
  CALL mp_bcast( xqq, ionode_id )
  !
#endif

  ! translate from input to internals of PWscf, various checks
  time_max = max_seconds

  ! ...   Set the number of species

  IF( ntyp < 1 .OR. ntyp > ntypx ) THEN
     CALL error(' iosys ',' ntyp out of range ', ntyp )
  END IF

  ! ...   IBRAV and CELLDM

  IF( ibrav /= 0 .and. celldm(1) == 0.d0 ) THEN
     CALL error(' iosys ',' invalid value in celldm ', 1 )
  END IF
  IF( ibrav < 0 .OR. ibrav > 14 ) THEN
     CALL error(' iosys ',' ibrav out of range ', 1 )
  END IF

  ! ...   Set Values for electron and bands

  SELECT CASE ( TRIM(occupations) )
     CASE ('fixed')
        ngauss = 0
        ltetra = .false.
        IF (degauss /= 0.d0) THEN
           CALL error(' iosys ',&
                ' fixed occupations, gauss. broadening ignored',-1 )
           degauss= 0.d0
        END IF
     CASE ('smearing')
        ltetra = .false.
        IF (degauss == 0.d0) THEN
           CALL error(' iosys ',&
                ' smearing requires gaussian broadening', 1 )
        END IF
        SELECT CASE ( TRIM(smearing) )
        CASE ('gaussian', 'gauss')
           ngauss = 0
        CASE ('methfessel-paxton', 'm-p', 'mp')
           ngauss = 1
        CASE ('marzari-vanderbilt', 'cold', 'm-v', 'mv')
           ngauss =-1
        CASE ('fermi-dirac', 'f-d', 'fd')
           ngauss = -99
        END SELECT
     CASE ('tetrahedra')
        ngauss = 0
        ltetra = .true.
     CASE DEFAULT
        CALL error(' iosys ',' occupations '//trim(occupations)// &
             & 'not implemented', 1 )
  END SELECT

  IF( nbnd < 1 ) THEN
     CALL error(' iosys ',' nbnd less than 1 ', nbnd )
  END IF
  IF( nelec < 0 ) THEN
     CALL error(' iosys ',' nelec less than 0 ', nelec )
  END IF
  IF( nspin < 1 .OR. nspin > 2 ) THEN
     CALL error(' iosys ',' nspin out of range ', nspin )
  END IF
  lsda = (nspin.eq.2)

  ! ...   Set Values for the cutoff

  IF( ecutwfc <= 0.d0 ) THEN
     CALL error(' iosys ',' invalid ecutwfc ', INT(ecutwfc) )
  END IF

  if (ecutrho <= 0.d0) then
     dual = 4.d0
  else
     dual = ecutrho/ecutwfc
     IF( dual <= 1.d0 ) THEN
        CALL error(' iosys ',' invalid dual? ', 1)
     end if
  end if

  SELECT CASE ( TRIM(restart_mode) )
  CASE ('from_scratch')
     restart = .false.
     startingconfig = 'input'
  CASE ('restart')
     restart = .true.
     startingpot = 'file'
     startingwfc = 'file'
     if (TRIM(ion_positions) .eq. 'from_input') then
        startingconfig = 'input'
     else
        startingconfig = 'file'
     end if
  CASE DEFAULT
     CALL error(' iosys ',' unknown restart_mode '//trim(restart_mode), 1 )
  END SELECT

  SELECT CASE ( TRIM(disk_io) )
  CASE ('high')
     reduce_io = .false.
  CASE DEFAULT
     reduce_io = .true.
     restart =.false.
  END SELECT

  Hubbard_U(:)    = Hubbard_U(:)/rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:)/rytoev

  SELECT CASE ( TRIM(calculation) )
  CASE ('scf' )
     lscf = .true.
     iswitch = 0
     lforce = tprnfor
     lmovecell=.false.
     nstep = 1
 CASE ('nscf')
     lscf = .false.
     iswitch = -1
     lforce = .false.
     lmovecell=.false.
     nstep = 1
  CASE ('relax')
     lscf = .true.
     iswitch = 1
     lforce = .true.
     lmovecell=.false.
  CASE ('md')
     lscf = .true.
     iswitch = 3
     lforce = .true.
     lmovecell=.false.
  CASE ('vc-relax','vc-md')
     lscf = .true.
     iswitch = 3
     lmovecell=.true.
     lforce = .true.
 CASE ('phonon')
     lscf = .false.
     iswitch = -2
     lforce = .false.
     lmovecell=.false.
     nstep = 1
  CASE DEFAULT
     CALL error(' iosys ',' calculation '//&
          trim(electron_dynamics)//' not implemented',1)
  END SELECT

  if (modenum /= 0) then
     iswitch = -4
  end if

  if ( startingpot.ne.'atomic' .and. startingpot.ne.'file' ) then
     call error(' iosys','wrong startingpot: use default',-1)
     if (      lscf ) startingpot = 'atomic'
     if ( .not.lscf ) startingpot = 'file'
  end if

  if ( .not.lscf .and. startingpot.ne.'file' ) then
     call error(' iosys','wrong startingpot: use default',-1)
     startingpot = 'file'
  end if

  if ( startingwfc.ne.'atomic' .and. startingwfc.ne.'random' .and. &
      startingwfc.ne.'file' ) then
     call error(' iosys','wrong startingwfc: use default',-1)
     startingwfc='atomic'
  end if

  SELECT CASE ( TRIM(electron_dynamics) )
  CASE ('none' )
     continue
  CASE DEFAULT
     CALL error(' iosys ',' unknown electron_dynamics '//&
          trim(electron_dynamics),1)
  END SELECT

  SELECT CASE ( TRIM(diagonalization) )
  CASE ('cg')
     isolve = 1
     max_cg_iter= diago_cg_maxiter
     loverlap =.false.
  CASE ('diis')
     isolve = 2
     diis_buff = diago_diis_buff
     diis_start_dav = diago_diis_start
     diis_wfc_keep  = diago_diis_keep
     loverlap =.true.
  CASE ('david', 'david_overlap')
     isolve = 0
     david = diago_david_ndim
     loverlap =.true.
  CASE ('david_nooverlap')
     isolve = 0
     david = diago_david_ndim
     loverlap =.false.
  CASE DEFAULT
     isolve = 0
     david = diago_david_ndim
     loverlap =.true.
  END SELECT

  ethr = 0.d0
  tr2  = conv_thr
  niter= electron_maxstep

  SELECT CASE ( TRIM(potential_extrapolation) )
  CASE ('none')
     order = 0
  CASE ('atomic')
     order = 1
  CASE ('wfc')
     order = 2
  CASE ('wfc2')
     order = 3
  CASE DEFAULT
     order = 1
  END SELECT

  IF ( occupations == 'fixed' .AND. nspin == 2  .AND. lscf ) THEN
     CALL error(' iosys ',&
                ' fixed occupations and lsda not implemented ', 1 )
  END IF

  calc = ' '
  if ( TRIM(calculation) == 'relax' ) then
     SELECT CASE ( TRIM(ion_dynamics) )
     CASE ('bfgs')
        iswitch = 1
        epse = etot_conv_thr
        epsf = forc_conv_thr
     CASE ('constrained-bfgs')
        iswitch = 2
        epse = etot_conv_thr
        epsf = forc_conv_thr
     CASE ('damp')
        iswitch = 3
        calc = 'mm'
        epse = etot_conv_thr
        epsf = forc_conv_thr
        ntcheck=nstep+1
     CASE DEFAULT
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dymanics='//trim(ion_dynamics)//' not supported', 1 )
     END SELECT
  endif
  if ( TRIM(calculation) == 'md' ) then
     SELECT CASE ( TRIM(ion_dynamics) )
     CASE ('verlet')
        iswitch = 3
     CASE ('constrained-verlet')
        iswitch = 4
     CASE ('beeman')
        iswitch = 3
        calc = 'md'
        ntcheck=nstep+1
     CASE DEFAULT
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dymanics='//trim(ion_dynamics)//' not supported', 1 )
     END SELECT
  endif
  if ( TRIM(calculation) == 'vc-relax' ) then
     SELECT CASE ( TRIM(cell_dynamics) )
     CASE ('none')
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'mm'
        ntcheck=nstep+1
     CASE ('damp-pr')
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'cm'
        ntcheck=nstep+1
     CASE ('damp-w')
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'nm'
        ntcheck=nstep+1
     CASE DEFAULT
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': cell_dymanics='//trim(cell_dynamics)//' not supported', 1 )
     END SELECT
     if ( TRIM(ion_dynamics) .ne. 'damp' ) then
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dymanics='//trim(ion_dynamics)//' not supported', 1 )
     end if
  endif
  if ( TRIM(calculation) == 'vc-md' ) then
     SELECT CASE ( TRIM(cell_dynamics) )
     CASE ('none')
        iswitch = 3
        calc = 'md'
        ntcheck=nstep+1
     CASE ('pr')
        iswitch = 3
        calc = 'cd'
        ntcheck=nstep+1
     CASE ('w')
        iswitch = 3
        calc = 'nd'
        ntcheck=nstep+1
     CASE DEFAULT
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dymanics='//trim(ion_dynamics)//' not supported', 1 )
     END SELECT
     if ( TRIM(ion_dynamics) .ne. 'beeman' ) then
        CALL error(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dymanics='//trim(ion_dynamics)//' not supported', 1 )
     end if
  endif

  !

  SELECT CASE ( TRIM(ion_temperature) )
  CASE ('not_controlled')
     continue
  CASE ('rescaling' )
     temperature = tempw
     ttol = tolp
  CASE DEFAULT
     CALL error(' iosys ',' unknown ion_temperature '//&
          trim(ion_temperature), 1 )
  END SELECT

  SELECT CASE ( TRIM(cell_temperature) )
  CASE ('not_controlled')
     continue
  CASE DEFAULT
     CALL error(' iosys ',' unknown cell_temperature '//&
          trim(cell_temperature), 1 )
  END SELECT

  SELECT CASE ( TRIM(cell_dofree) )
  CASE ('all')
     continue
  CASE DEFAULT
     CALL error(' iosys ',' unknown cell_dofree '//trim(cell_dofree), 1 )
  END SELECT

  SELECT CASE ( TRIM(mixing_mode) )
  CASE ('plain')
     imix = 0
     starting_scf_threshold = tr2
  CASE ('TF')
     imix = 1
     starting_scf_threshold = tr2
  CASE ('local-TF')
     imix = 2
     starting_scf_threshold = tr2
  CASE ('potential')
     imix = -1
     starting_scf_threshold = sqrt(tr2)
  CASE DEFAULT
     CALL error(' iosys ',' unknown mixing '//trim(mixing_mode), 1)
  END SELECT
  nmix = mixing_ndim
  niter_with_fixed_ns = mixing_fixed_ns

  SELECT CASE ( TRIM(verbosity) )
  CASE ('high')
     iverbosity = 1
  CASE DEFAULT
     iverbosity = 0
  END SELECT
  tmp_dir = trim(outdir)
  lstres = tstress .AND. lscf

  ! read following cards

  allocate (tau(3,nat))
  allocate (ityp(nat))
  allocate (force(3,nat)) ! compatibility with old readin
  !
  CALL read_cards (pseudop, atomic_positions)
  !
  ! set up atomic positions and crystal lattice
  !
  if (ibrav == 0) then
     if (celldm (1) == 0.d0) then
        celldm (1) = sqrt(at(1,1)**2+at(1,2)**2+at(1,3)**2)
        at(:,:) = at(:,:) / celldm(1)
     end if
  else
     CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3))
  end if
  alat = celldm(1)
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  SELECT CASE ( atomic_positions )
     !
     !  convert input atomic positions to internally used format:
     !  tau in a0 units
     !
  CASE ('alat')
     !
     !  input atomic positions are divided by a0: do nothing
     !
     ltaucry = .false.
  CASE ('bohr')
     !
     !  input atomic positions are in a.u.: divide by alat
     !
     tau = tau/alat
     ltaucry = .false.
  CASE ('crystal')
     !
     !  input atomic positions are in crystal axis
     !
     ltaucry = .true.
  CASE ('angstrom')
     !
     !  atomic positions in A: convert to a.u. and divide by alat
     !
     tau = tau/0.529177/alat
     ltaucry = .false.
  CASE DEFAULT
     CALL error(' iosys ',' atomic_positions='//trim(atomic_positions)// &
          ' not implemented ', 1 )
  END SELECT

  ! If  ltaucry = .true. , the input atomic positions in crystallographic
  ! units are transformed in cartesian coordinates
  !
  if (ltaucry) call cryst_to_cart (nat, tau, at, 1)
  !
  ! set default value of wmass
  !
  if (wmass.eq.0.d0) then
     if (calc.eq.'nd' .or. calc.eq.'nm') then
        do ia=1,nat
           wmass = wmass + amass(ityp(ia))
        end do
        wmass =  0.75d0 * wmass / pi / pi / omega**(2.d0/3.d0)
     end if
     if (calc.eq.'cd' .or. calc.eq.'cm') then
        do ia=1,nat
           wmass = wmass + amass(ityp(ia))
        end do
        wmass =  0.75d0 * wmass / pi / pi
     end if
  end if
  !
  !
  ! unit conversion for cell mass and pressure
  !
  cmass = wmass * amconv
  press = press / uakbar
  !
  !    read pseudopotentials
  !
  CALL readpp(pseudo_dir,pseudop)
  !
  !  Renata's dynamics uses masses in atomic units
  !
  if (calc.ne.' ') then
     amass = amass * amconv
  end if
  !
  !
  !    In the case of variable cell dynamics save old cell variables
  !    and initialize a few other variables
  !
  if (lmovecell) then
     at_old = at
     omega_old = omega
     lstres = .true.
     if (cell_factor.le.0.d0) cell_factor = 1.2d0
     if (cmass.le.0.d0) call error('readin',&
          &      'vcsmd: a positive value for cell mass is required',1)
  else
     cell_factor = 1.d0
  end if
  !
  call verify_tmpdir
  !
  write (6,'(/5x,"current restart_mode = ",a)') trim(restart_mode)
  write (6,'( 5x,"current disk_io mode = ",a)') trim(disk_io)
  call restart_from_file
  if (startingconfig.eq.'file') call read_config_from_file
  call write_config_to_file

  !
  ! Files
  !
  input_drho = ' '
  output_drho= ' '
  crystal = ' '
  !
  return
end subroutine iosys
!
!-----------------------------------------------------------------------
subroutine read_cards (pseudop, atomic_positions)
  !-----------------------------------------------------------------------
  !
  use parser
  use pwcom, only: at, ityp, tau, nat, ntyp, atm, amass, ibrav, &
       nks, nk1, nk2, nk3, k1, k2, k3, xk, wk, lxkcry, symm_type, &
       fixatom
  implicit none
  !
  character (len=30) :: pseudop (ntyp), atomic_positions
  !
  real(kind=8), allocatable :: tau_inp(:,:)
  integer, allocatable :: ityp_inp(:), iforce_inp(:,:)
  character(len=3) :: atom_label(ntyp), lb_pos
  character(len=256) :: line, input_line
  logical :: tcell=.false., tatms=.false., tatmp=.false., tend
  integer :: unit = 5, i, is, ns, ia, ios, ik, nf

  amass = 0

100 CALL read_line (line, end_of_file = tend)
  if (tend) go to 200
  if (line(1:1).eq.'#') go to 100
  do i=1,len_trim(line)
     line(i:i) = capital(line(i:i))
  end do
  if (matches('ATOMIC_SPECIES',line)) then

     do is = 1, ntyp

        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300

        read (input_line,*) atom_label(is), amass(is), pseudop(is)

        IF( amass(is) <= 0.d0 ) THEN
           CALL error(' iosys ',' invalid  mass ', is)
        END IF
        atm(is) = atom_label(is)

     end do
     tatms =.true.

  else if (matches('ATOMIC_POSITIONS',line)) then

     allocate ( ityp_inp(nat) )
     allocate ( tau_inp(3,nat) )
     allocate ( iforce_inp(3,nat) )
     tau_inp = 0.d0
     ityp_inp= 0
     do ia = 1, nat

        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        CALL field_count(nf,input_line)

        if (nf.eq.4) then
           READ(input_line,*) lb_pos, &
                tau_inp(1,ia), tau_inp(2,ia), tau_inp(3,ia)
           iforce_inp(:,ia) = 1
        else if (nf.eq.7) then
           READ (input_line, * ) lb_pos, &
                tau_inp(1,ia), tau_inp(2,ia), tau_inp(3,ia), &
                iforce_inp(1,ia), iforce_inp(2,ia), iforce_inp(3,ia)
        else
           call error (' cards','wrong number of tokens', ia)
        end if

        match_label: DO is = 1, ntyp
           IF( lb_pos == atom_label(is) ) THEN
              ityp_inp(ia) = is
              EXIT match_label
           END IF
        END DO match_label

        if (ityp_inp(ia) <= 0 .OR. ityp_inp(ia) > ntyp) &
             call error (' cards','wrong atomic positions', ia)

     end do
     !
     ! TEMP: calculate fixatom
     !
     fixatom = 0
     fix1: DO ia = nat, 1, -1
        if ( iforce_inp(1,ia) /= 0 .or. &
             iforce_inp(2,ia) /= 0 .or. &
             iforce_inp(3,ia) /= 0 ) EXIT fix1
        fixatom = fixatom + 1
     END DO fix1
     !
     !  read option to card ATOMIC_POSITIONS
     !
     if ( matches('ALAT', line) ) then
        atomic_positions = 'alat'
     else if ( matches('BOHR', line) ) then
        atomic_positions = 'bohr'
     else if ( matches('CRYSTAL', line) ) then
        atomic_positions = 'crystal'
     else if ( matches('ANGSTROM', line) ) then
        atomic_positions = 'angstrom'
     else
        atomic_positions = 'alat'
     end if
     !
     tau = tau_inp
     ityp= ityp_inp
     !
     deallocate (iforce_inp)
     deallocate (tau_inp)
     deallocate (ityp_inp)
     tatmp =.true.

  else if (matches('OCCUPATIONS',line)) then
!!!         allocate (f(n))
!!!         read (unit, *) (f(i), i=1,n)
     call error (' cards','card OCCUPATIONS not implemented', 1)

  else if (matches('K_POINTS',line)) then
     !
     !  read option to card K_POINTS
     !
     if ( matches('AUTOMATIC',line) ) then
        !
        !  automatic generation of k-points
        !
        lxkcry = .false.
        nks = 0
        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        read (input_line, *) nk1, nk2, nk3, k1, k2 ,k3
     else if ( matches('TPIBA',line) ) then
        !
        !  input k-points are in 2pi/a units
        !
        lxkcry = .false.
        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        read (input_line, *) nks
        do ik=1, nks
           CALL read_line (input_line, end_of_file = tend)
           if (tend) go to 300
           read (input_line, *) xk(1,ik), xk(2,ik), xk(3,ik), wk(ik)
        end do
     else if ( matches('CRYSTAL',line) ) then
        !
        !  input k-points are in crystal (reciprocal lattice) axis
        !
        lxkcry = .true.
        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        read (input_line, *) nks
        do ik=1, nks
           CALL read_line (input_line, end_of_file = tend)
           if (tend) go to 300
           read (input_line, *) xk(1,ik), xk(2,ik), xk(3,ik), wk(ik)
        end do
     else if ( matches('GAMMA',line) ) then
        !
        !  Only Gamma (k=0) is used
        !
        lxkcry = .false.
        nks = 1
        xk(:,1) = 0.0
        wk(1)   = 1.0
     else
        !
        !  default: input k-points are in 2pi/a units
        !
        lxkcry = .false.
        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        read (input_line, *) nks
        do ik=1, nks
           CALL read_line (input_line, end_of_file = tend)
           if (tend) go to 300
           read (input_line, *) xk(1,ik), xk(2,ik), xk(3,ik), wk(ik)
        end do
     end if

  else if (matches('CELL_PARAMETERS',line)) then

     if (matches('HEXAGONAL',line)) then
        symm_type = 'hexagonal'
     else
        symm_type = 'cubic'
     end if

     do i=1,3
        CALL read_line (input_line, end_of_file = tend)
        if (tend) go to 300
        read(input_line,*) at(1,i),at(2,i),at(3,i)
     enddo
     tcell = .true.

  else

     write (6,'(a)') 'Warning: card '//trim(line)//' ignored'

  end if

  go to 100

200 if (ibrav == 0 .and. .not.tcell ) &
       CALL error(' cards ',' ibrav=0: must read cell parameters', 1 )
  if (ibrav /= 0 .and. tcell ) &
       CALL error(' cards ',' redundant data for cell parameters', 2 )
  if (.not.tatms ) &
       CALL error(' cards ',' atomic species info missing', 1 )
  if (.not.tatmp ) &
       CALL error(' cards ',' atomic position info missing', 1 )

  return
300 CALL error(' cards ',' unexpected end of file', 1 )
end subroutine read_cards
!
!-----------------------------------------------------------------------
subroutine verify_tmpdir
  !-----------------------------------------------------------------------
  !
  use io, only: tmp_dir, nd_nmbr
  implicit none
  integer :: l, ios
  !
  !    verify if tmp_dir ends with /, add one if needed
  !
  l=len_trim(tmp_dir)
  if (tmp_dir(l:l).ne.'/') then
     if (l > 0 .and. l < len(tmp_dir)) then
        tmp_dir(l+1:l+1)='/'
     else
        call error('reading',tmp_dir//' truncated or empty',1)
     end if
  end if
  ios = 0
  open (unit=4,file=trim(tmp_dir)//'pwscf'// nd_nmbr,status='unknown',&
        form='unformatted', iostat=ios)
  if (ios /= 0 ) call error('reading', &
                    trim(tmp_dir)//' non existent or non writable',1)
  close (unit=4,status='delete')
  !
  return
  end
