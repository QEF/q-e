!
! Copyright (C) 2002-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "machine.h"
!-----------------------------------------------------------------------
subroutine iosys
  !-----------------------------------------------------------------------
  !   this subroutine reads input data from standard input (unit 5)
  !     ------------------------------------------------------------------

  !
  !    access the modules renaming the variables that have the same name 
  !    as the input parameters, this is required in order to use a conde
  !    independent input parser
  ! 

  use bp, only: &
     nppstr_ => nppstr, &
     gdir_   => gdir, &
     lberry_ => lberry
  use brilz, only: &
     at, alat, omega, &
     celldm_ => celldm, &
     ibrav_  => ibrav
  use basis, only: &
     nat_  => nat, &
     ntyp_ => ntyp, &
     ityp, tau, atomic_positions, atm, &
     startingwfc_ => startingwfc, &
     startingpot_ => startingpot, &
     startingconfig
  use char, only: &
     title_ => title, &
     crystal
  use cellmd, only: &
     cmass, ttol, omega_old, at_old, ntcheck, &
     cell_factor_ => cell_factor , &
     press_ => press, &
     calc, lmovecell
  use constants, only: &
     pi, rytoev, uakbar, amconv, bohr_radius_angs
  use dynam, only: &
     dt_ => dt, &
     temperature, amass, delta_t, nraise
  use extfield, only : &
     tefield_ => tefield, &
     dipfield_ => dipfield, &
     edir_ => edir, &
     emaxpos_ => emaxpos, &
     eopreg_ => eopreg, &
     eamp_ => eamp, &
     forcefield
  use filnam, only : &
     input_drho, output_drho
  use force_mod, only: &
     lforce, lstres, force
  use gvect, only: &
     dual, &
     nr1_ => nr1, &
     nr2_ => nr2, &
     nr3_ => nr3,  &
     ecutwfc_ => ecutwfc, &
     ecfixed_ => ecfixed, &
     qcutz_ => qcutz, &
     q2sigma_ => q2sigma
  use gsmooth, only : &
     nr1s_ => nr1s, &
     nr2s_ => nr2s, &
     nr3s_ => nr3s
  use klist, only: &
     xk, wk, nks, ngauss,&
     xqq_ => xqq, &
     degauss_ => degauss, &
     nelec_ => nelec
  use ktetra, only : &
     nk1, nk2, nk3, k1, k2, k3, ltetra
  use ldaU, only: &
     Hubbard_U_ => hubbard_u, &
     Hubbard_alpha_ => hubbard_alpha, &
     niter_with_fixed_ns, &
     lda_plus_u_ => lda_plus_u
  use lsda_mod, only: &
     nspin_ => nspin, &
     starting_magnetization_ => starting_magnetization, &
     lsda
  use io, only : &
     tmp_dir, &
     prefix_ => prefix, &
     pseudo_dir_ => pseudo_dir, &
     pseudop
  use relax, only: &
     epsf, starting_scf_threshold, restart_bfgs, epse
  use varie, only: &
     diis_start_cg, diis_ndim, diis_wfc_keep, isolve, max_cg_iter, &
     diis_buff, david, imix, nmix, iverbosity, tr2, niter, order, iswitch, &
     ntypx, &
     upscale_ => upscale, &
     mixing_beta_ => mixing_beta, &
     nstep_ => nstep, &
     iprint_ => iprint, &
     nosym_ => nosym, &
     modenum_ => modenum, &
     reduce_io, ethr, lscf, noinv, time_max, restart
  use wvfct, only: &
     nbnd_ => nbnd
  use fixed_occ, only : &
      tfixed_occ
  !
  ! CONTROL namelist

  use input_parameters, only: &
       title, calculation, verbosity, &
       restart_mode, nstep, iprint, tstress, tprnfor, &
       dt, outdir, prefix, max_seconds, &
       etot_conv_thr, forc_conv_thr, pseudo_dir, disk_io, tefield, &
       dipfield, lberry, gdir, nppstr

  ! SYSTEM namelist

  use input_parameters, only: &
       ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
       nat, ntyp, nbnd, nelec, ecutwfc, ecutrho, &
       nr1, nr2, nr3, nr1s, nr2s, nr3s, &
       nosym, starting_magnetization, &
       occupations, degauss, smearing, &
       nspin, ecfixed, qcutz, q2sigma, &
       lda_plus_U, Hubbard_U, Hubbard_alpha, &
       edir, emaxpos, eopreg, eamp

  ! ELECTRONS namelist

  use input_parameters, only: &
       electron_maxstep, electron_dynamics, &
       mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns, &
       diago_cg_maxiter, diago_david_ndim, diago_diis_buff, &
       diago_diis_start, diago_diis_ndim, diago_diis_keep, diagonalization, &
       startingwfc, startingpot, conv_thr

  ! IONS namelist

  use input_parameters, only: &
       ion_dynamics, ion_positions, ion_temperature, &
       tempw, tolp, upscale, potential_extrapolation

  ! CELL namelist

  use input_parameters, only: &
       cell_parameters, cell_dynamics, press, &
       wmass, cell_temperature, cell_dofree, cell_factor

  ! PHONON namelist

  use input_parameters, only: phonon, &
       modenum, xqq

  use read_namelists_module, only: &
       read_namelists

  !
  implicit none
  !

  !
  ! local variables
  !

  integer :: unit = 5, i, ia, ios, ierr, ilen, is

  !
  call getenv('HOME',pseudo_dir)
  pseudo_dir=trim(pseudo_dir)//'/pw/pseudo/'

  CALL read_namelists( 'PW' )

  nraise = 100
  delta_t = 1.d0
  !
  ! translate from input to internals of PWscf, various checks
  time_max = max_seconds

  if (tefield.and..not.nosym) then
     nosym=.true.
     write(6,'(5x,"Presently no symmetry can be used with electric field",/)')
  endif
  if (tefield.and.(tstress)) then
     tstress=.false.
     write(6,'(5x,"Presently stress not available with electric field",/)')
  endif
  if (tefield.and.(nspin.eq.2)) then
     call errore('input','LSDA not available with electric field',1)
  endif

  ! ...   Set Values for electron and bands

  tfixed_occ=.false.
  SELECT CASE ( TRIM(occupations) )
     CASE ('fixed')
        ngauss = 0
        ltetra = .false.
        IF (degauss /= 0.d0) THEN
           CALL errore(' iosys ',&
                ' fixed occupations, gauss. broadening ignored',-1 )
           degauss= 0.d0
        END IF
     CASE ('smearing')
        ltetra = .false.
        IF (degauss == 0.d0) THEN
           CALL errore(' iosys ',&
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
     CASE ('from_input')
        ngauss=0
        ltetra=.false.
        tfixed_occ=.true.
     CASE DEFAULT
        CALL errore(' iosys ',' occupations '//trim(occupations)// &
             & 'not implemented', 1 )
  END SELECT

  IF( nbnd < 1 ) THEN
     CALL errore(' iosys ',' nbnd less than 1 ', nbnd )
  END IF
  IF( nelec < 0 ) THEN
     CALL errore(' iosys ',' nelec less than 0 ', nelec )
  END IF
  lsda = (nspin.eq.2)

  if (ecutrho <= 0.d0) then
     dual = 4.d0
  else
     dual = ecutrho/ecutwfc
     IF( dual <= 1.d0 ) THEN
        CALL errore(' iosys ',' invalid dual? ', 1)
     end if
  end if

  SELECT CASE ( TRIM(restart_mode) )
  CASE ('from_scratch')
     restart = .false.
     restart_bfgs = .false.
     startingconfig = 'input'
  CASE ('restart')
     restart = .true.
     restart_bfgs = .true.
     startingpot = 'file'
     startingwfc = 'file'
     if (TRIM(ion_positions) .eq. 'from_input') then
        startingconfig = 'input'
     else
        startingconfig = 'file'
     end if
  CASE DEFAULT
     CALL errore(' iosys ',' unknown restart_mode '//trim(restart_mode), 1 )
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

  ethr = 0.d0
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
!     ethr = 1.d-6
! I think ethr should not be more strict than that in a simple band
! structure calculation but there is still something unsatisfactory 
! in the Davidson diagonalization convergence. SdG 20/03/2003
!
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
     CALL errore(' iosys ',' calculation '//&
          trim(calculation)//' not implemented',1)
  END SELECT

  if (modenum /= 0) then
     iswitch = -4
  end if

  if ( startingpot.ne.'atomic' .and. startingpot.ne.'file' ) then
     call errore(' iosys','wrong startingpot: use default',-1)
     if (      lscf ) startingpot = 'atomic'
     if ( .not.lscf ) startingpot = 'file'
  end if

  if ( .not.lscf .and. startingpot.ne.'file' ) then
     call errore(' iosys','wrong startingpot: use default',-1)
     startingpot = 'file'
  end if

  if ( startingwfc.ne.'atomic' .and. startingwfc.ne.'random' .and. &
      startingwfc.ne.'file' ) then
     call errore(' iosys','wrong startingwfc: use default',-1)
     startingwfc='atomic'
  end if

  SELECT CASE ( TRIM(electron_dynamics) )
  CASE ('none' )
     continue
  CASE DEFAULT
     CALL errore(' iosys ',' unknown electron_dynamics '//&
          trim(electron_dynamics),1)
  END SELECT

  SELECT CASE ( TRIM(diagonalization) )
  CASE ('cg')
     isolve = 1
     max_cg_iter= diago_cg_maxiter
  CASE ('diis')
     isolve = 2
     max_cg_iter= diago_cg_maxiter
     diis_buff = diago_diis_buff
     diis_start_cg = diago_diis_start   ! SF
     diis_ndim = diago_diis_ndim        ! SF
     diis_wfc_keep  = diago_diis_keep
  CASE ('david')
     isolve = 0
     david = diago_david_ndim
  CASE DEFAULT
     isolve = 0
     david = diago_david_ndim
  END SELECT

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
     CALL errore(' iosys ',&
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
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dynamics='//trim(ion_dynamics)//' not supported', 1 )
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
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dynamics='//trim(ion_dynamics)//' not supported', 1 )
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
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': cell_dynamics='//trim(cell_dynamics)//' not supported', 1 )
     END SELECT
     if ( TRIM(ion_dynamics) .ne. 'damp' ) then
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dynamics='//trim(ion_dynamics)//' not supported', 1 )
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
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dynamics='//trim(ion_dynamics)//' not supported', 1 )
     END SELECT
     if ( TRIM(ion_dynamics) .ne. 'beeman' ) then
        CALL errore(' iosys ','calculation='//trim(calculation)// &
&             ': ion_dynamics='//trim(ion_dynamics)//' not supported', 1 )
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
     CALL errore(' iosys ',' unknown ion_temperature '//&
          trim(ion_temperature), 1 )
  END SELECT

  SELECT CASE ( TRIM(cell_temperature) )
  CASE ('not_controlled')
     continue
  CASE DEFAULT
     CALL errore(' iosys ',' unknown cell_temperature '//&
          trim(cell_temperature), 1 )
  END SELECT

  SELECT CASE ( TRIM(cell_dofree) )
  CASE ('all')
     continue
  CASE DEFAULT
     CALL errore(' iosys ',' unknown cell_dofree '//trim(cell_dofree), 1 )
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
     CALL errore(' iosys ',' unknown mixing '//trim(mixing_mode), 1)
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

  !  Copy values from input module to PW internals

  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  title_      = title
  dt_         = dt
  tefield_    = tefield
  dipfield_   = dipfield
  prefix_     = TRIM(prefix)
  pseudo_dir_ = TRIM(pseudo_dir)
  nstep_      = nstep
  iprint_     = iprint

  celldm_ = celldm
  ibrav_ = ibrav
  nat_ = nat 
  ntyp_ = ntyp
  edir_ = edir
  emaxpos_ = emaxpos
  eopreg_ = eopreg
  eamp_ = eamp
  nr1_ = nr1
  nr2_ = nr2
  nr3_ = nr3
  ecutwfc_ = ecutwfc
  ecfixed_ = ecfixed
  qcutz_ = qcutz
  q2sigma_ = q2sigma
  nr1s_ = nr1s
  nr2s_ = nr2s
  nr3s_ = nr3s
  degauss_ = degauss
  nelec_ = nelec
  Hubbard_U_( 1 : ntyp ) = hubbard_u( 1 : ntyp )
  Hubbard_alpha_ ( 1 : ntyp )= hubbard_alpha( 1 : ntyp )
  lda_plus_u_ = lda_plus_u
  nspin_ = nspin
  starting_magnetization_ = starting_magnetization
  nosym_ = nosym
  nbnd_  = nbnd

  startingwfc_ = startingwfc
  startingpot_ = startingpot
  mixing_beta_ = mixing_beta

  upscale_ = upscale
  press_ = press
  cell_factor_ = cell_factor
  modenum_ = modenum
  xqq_ = xqq

  ! read following cards

  allocate ( tau( 3, nat_ ) )
  allocate ( ityp( nat_ ) )
  allocate ( force( 3, nat_ ) ) ! compatibility with old readin
  if ( tefield ) allocate( forcefield( 3, nat_ ) )
  !
  CALL read_cards (pseudop, atomic_positions)
  !
  ! set up atomic positions and crystal lattice
  !
  if (celldm_ (1) == 0.d0 .and. a /= 0.d0) then
     if (ibrav_ == 0) ibrav = 14 
     celldm_ (1) = a*bohr_radius_angs
     celldm_ (2) = b/a
     celldm_ (3) = c/a
     celldm_ (4) = cosab
     celldm_ (5) = cosac
     celldm_ (6) = cosbc 
  else if (celldm_ (1) /= 0.d0 .and. a /= 0.d0) then
     call errore('input', ' do not specify both celldm and a,b,c!',1)
  end if
  !
  if (ibrav_ == 0 .and. celldm_ (1) /= 0.d0) then
     ! input at are in units of alat
     alat = celldm_(1)
  else if (ibrav_ == 0 .and. celldm_ (1) == 0.d0) then
     ! input at are in atomic units: define alat
     celldm_ (1) = sqrt(at(1,1)**2+at(1,2)**2+at(1,3)**2)
     alat = celldm_(1)
     ! bring at to alat units
     at(:,:) = at(:,:) / alat
  else
     ! generate at (atomic units)
     CALL latgen(ibrav,celldm_,at(1,1),at(1,2),at(1,3),omega)
     alat = celldm_(1) 
     ! bring at to alat units
     at(:,:) = at(:,:) / alat
  end if
  !
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
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
  CASE ('bohr')
     !
     !  input atomic positions are in a.u.: divide by alat
     !
     tau = tau/alat
  CASE ('crystal')
     !
     !  input atomic positions are in crystal axis
     !
     call cryst_to_cart ( nat_ , tau, at, 1)
  CASE ('angstrom')
     !
     !  atomic positions in A: convert to a.u. and divide by alat
     !
     tau = tau/bohr_radius_angs/alat
  CASE DEFAULT
     CALL errore(' iosys ',' atomic_positions='//trim(atomic_positions)// &
          ' not implemented ', 1 )
  END SELECT

  !
  ! set default value of wmass
  !
  if ( wmass .eq. 0.d0 ) then
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
  cmass  = wmass * amconv
  press_ = press_ / uakbar
  !
  !    read pseudopotentials
  !
  CALL readpp
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
     if ( cell_factor_ .le. 0.d0 ) cell_factor_ = 1.2d0
     if (cmass.le.0.d0) call errore('readin',&
          &      'vcsmd: a positive value for cell mass is required',1)
  else
     cell_factor_ = 1.d0
  end if
  !
  call verify_tmpdir
  !
!  write (6,'(/5x,"current restart_mode = ",a)') trim(restart_mode)
!  write (6,'( 5x,"current disk_io mode = ",a)') trim(disk_io)
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

  !
  return
end subroutine iosys

!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
subroutine read_cards ( pseudop, atomic_positions_ )
  !-----------------------------------------------------------------------
  !
  use parser
  use pwcom, only: at, ityp, tau, nat, ntyp, atm, amass, ibrav, nks, &
       nk1_ => nk1, &
       nk2_ => nk2, &
       nk3_ => nk3, &
       k1_  => k1,  &
       k2_  => k2,  &
       k3_  => k3,  &
       xk_  => xk,  &
       wk_  => wk,  &
       f_inp_ => f_inp, &
       lxkcry, symm_type, fixatom, tfixed_occ
 
  use input_parameters, only: &
       atom_label, atom_pfile, atom_mass, atom_ptyp, taspc, &
       tapos, rd_pos, atomic_positions, if_pos, sp_pos, &
       k_points, xk, wk, nk1, nk2, nk3, k1, k2, k3, nkstot, &
       cell_symmetry, rd_ht, trd_ht, f_inp

  use read_cards_module, only: read_cards_base => read_cards

  implicit none
  !
  character (len=80) :: pseudop (ntyp)
  character (len=30) :: atomic_positions_
  !
  logical :: tcell=.false.
  integer :: i, is, ns, ia, ik

  amass = 0

  call read_cards_base( 'PW' )

  if (.not.taspc ) &
       CALL errore(' cards ',' atomic species info missing', 1 )
  if (.not.tapos ) &
       CALL errore(' cards ',' atomic position info missing', 1 )

  do is = 1, ntyp
    amass( is ) = atom_mass( is )
    pseudop( is ) = atom_pfile( is )
    atm(is) = atom_label(is)
    IF( amass(is) <= 0.d0 ) THEN
      CALL errore(' iosys ',' invalid  mass ', is)
    END IF
  end do

  do ia = 1, nat
    tau( : , ia ) = rd_pos( : , ia )
    ityp( ia ) = sp_pos( ia )
  end do

  !
  ! TEMP: calculate fixatom
  !
  fixatom = 0
  fix1: DO ia = nat, 1, -1
    if ( if_pos(1,ia) /= 0 .or. &
         if_pos(2,ia) /= 0 .or. &
         if_pos(3,ia) /= 0 ) EXIT fix1
    fixatom = fixatom + 1
  END DO fix1

  atomic_positions_ = TRIM(atomic_positions)

  if ( k_points == 'automatic' ) then
    !  automatic generation of k-points
    lxkcry = .false.
    nks = 0
    nk1_ = nk1
    nk2_ = nk2
    nk3_ = nk3
    k1_  = k1
    k2_  = k2
    k3_  = k3
  else if ( k_points == 'tpiba' ) then
    !  input k-points are in 2pi/a units
    lxkcry = .false.
    nks  = nkstot
    xk_( :, 1 : nks )  = xk( :, 1 : nks )
    wk_( 1 : nks )  = wk( 1 : nks )
  else if ( k_points == 'crystal' ) then
    !  input k-points are in crystal (reciprocal lattice) axis
    lxkcry = .true.
    nks  = nkstot
    xk_( :, 1 : nks )  = xk( :, 1 : nks )
    wk_( 1 : nks )  = wk( 1 : nks )
  else if ( k_points == 'gamma' ) then
    !  Only Gamma (k=0) is used
    lxkcry = .false.
    nks = 1
    xk_(:,1) = 0.0
    wk_(1)   = 1.0
  else
    !  default: input k-points are in 2pi/a units
    lxkcry = .false.
    nks  = nkstot
    xk_( :, 1 : nks )  = xk( :, 1 : nks )
    wk_( 1 : nks )  = wk( 1 : nks )
  end if

  if (tfixed_occ) then
     if (nks.gt.1.or.nk1*nk2*nk3.gt.1) call errore('read_cards', &
                          'only one k point with fixed occupations',1)
     f_inp_=f_inp
  endif

  
  if ( trd_ht ) then

    symm_type = cell_symmetry 
    at = TRANSPOSE( rd_ht )
    tcell = .true.

  end if


  if (ibrav == 0 .and. .not.tcell ) &
       CALL errore(' cards ',' ibrav=0: must read cell parameters', 1 )
  if (ibrav /= 0 .and. tcell ) &
       CALL errore(' cards ',' redundant data for cell parameters', 2 )

  return
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
        call errore('reading',tmp_dir//' truncated or empty',1)
     end if
  end if
  ios = 0
  open (unit=4,file=trim(tmp_dir)//'pwscf'// nd_nmbr,status='unknown',&
        form='unformatted', iostat=ios)
  if (ios /= 0 ) call errore('reading', &
                    trim(tmp_dir)//' non existent or non writable',1)
  close (unit=4,status='delete')
  !
  return
  end
