!
! Copyright (C) 2020-2020 Quantum ESPRESSO group
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!------------------------------------------------------------------------------
PROGRAM dvscf_q2r
!------------------------------------------------------------------------------
!!
!! dvscf_q2r.x
!!    Fourier transform the dvscf computed at coarse q points to real space.
!!    Originally proposed by [1]. For charge neutrality correction see [2].
!!    Dipole long-range part described in [3].
!!    Quadrupole long-range part described in [4] (not implemented here).
!!
!! Input data: Namelist "input"
!!    prefix  : Prepended to input/output filenames; must be the same used
!!              in the calculation of the phonon code.
!!              (character, Default: 'pwscf')
!!    outdir  : Directory containing input, output, and scratch files; must
!!              be the same as specified in the calculation of ph.x.
!!              (character, Default: value of the ESPRESSO_TMPDIR environment
!!               variable if set; current directory ('./') otherwise)
!!    fildyn  : File where the dynamical matrix is written. Normally, this
!!              should be the same as specified on the input to ph.x.
!!              Only "fildyn"0 is used here.
!!              (character, Must be specified)
!!    fildvscf : File where the potential variation is written. This should be
!!               the same as specified on the input to ph.x.
!!               (character, Default: 'dvscf')
!!    wpot_dir : Directory where the w_pot binary files are written. Real space
!!               potential files are stored in wpot_dir with names
!!               prefix.wpot.irc$irc//"1".
!!               (character, Default: outdir // 'w_pot/')
!!    do_long_range : If .true., subtract the long-range part of the potential
!!                    before interpolation. Requires epsilon and Born effective
!!                    charge data in _ph0/prefix.phsave/tensor.xml.
!!                    (logical, Default: .false.)
!!    do_charge_neutral : If .true., renormalize phonon potential to impose
!!                    neutrality of Born effective charges. See [2] for
!!                    details. Both the Hartree and exchange-correlation
!!                    parts are renormalized, while in [2] only the
!!                    Hartree part is renormalized.
!!                    (logical, Default: .false.)
!!    verbosity : If 'high', write more information to stdout. Used by the
!!                test-suite.
!!                (character, Default: 'default'. Only 'high' is allowed.)
!!
!! [1] A. Eiguren and C. Ambrosch-Draxl, PRB 78, 045124 (2008)
!! [2] S. Ponce et al, J. Chem. Phys. 143, 102813 (2015)
!! [3] Xavier Gonze et al, Comput. Phys. Commun. 248 107042 (2020)
!! [4] Guillaume Brunin et al, arXiv:2002.00628 (2020)
!!
!! Not implemented for the following cases:
!!    - PAW
!!    - DFPT+U
!!    - magnetism (both collinear and noncollinear magnetism)
!!    - 2d Coulomb cutoff
!!
!! dvscf = dvscf_ind + dvscf_bare (All are lattice periodic.)
!! dvscf_ind: computed in the ph.x run, read from files.
!! dvscf_bare: computed on the fly by subroutine calc_dvscf_bare.
!! All potentials are computed in the Cartesian basis of atomic displacements.
!!
!! * Charge neutrality correction
!! If the sum of Born effective charge zeu is not zero, the potential has
!! has non-physical divergence around q = Gamma. To correct this problem,
!! renormalize the Hartree term by a q-dependent constant factor.
!! Since the dvscf file contains the sum of Hartree and xc contribution,
!! we renormalize the xc term as well as the Hartree term.
!! To renormalize only the Hartree term, one need to read drho and
!! compute the corresponding dvscf.
!!
!! dvscf_ind_ren(q,iat,idir) = dvscf_ind(q,iat,idir) * coeff
!! coeff = ( Z*q_idir - sum_jdir (Z* - Z*avg)_{idir,jdir} * q_jdir / epsil_q )
!!       / ( Z*q_idir - sum_jdir (Z*)_{idir,jdir} * q_jdir / epsil_q )
!!
!! epsil_q = 1/q^2 * (q.T * epsil * q): dielectric constant
!! Z = upf(nt)%zp: bare valence charge of the atom iat
!! Z*(:,:) = zeu(:,:,iat): Born effective charge. Read from fildyn
!! Z*avg(:,:) = 1 / nat * sum_jatm zeu(:,:,jatm): average zeu
!!
!! * Long-range part correction
!! dvlong: long-range part. Subtracted for smooth Fourier interpolation.
!! Taken from Eq.(13) of Ref. [3]
!! dvlong(G,q)_{a,x} = 1j * 4pi / Omega
!!                   * [ (q+G)_y * Zstar_{a,yx} * exp(-i*(q+G)*tau_a)) ]
!!                   / [ (q+G)_y * epsilon_yz * (q+G)_z ]
!!  a: atom index, x, y: Cartesian direction index
!!
!! w_pot(r,R) = 1/N_q * sum_q exp(-iqR) exp(iqr) (dvscf(r,q) - dvlong(r,q))
!! w_pot are computed and written to file.
!!
!! Later, dvtot at fine q points can be computed as
!! dvscf(r,q) = exp(-iqr) (dvlong(r,q) + sum_R exp(iqR) w_pot(r,R))
!!
!! Only the dipole (Frohlich) potential is considered. The quadrupole
!! potential [4] is not implemented.
!!
!! * Parallelization
!! We use PW and pool parallelization.
!! Here, the pool parallelization is for the q points, not for the k points.
!!
!------------------------------------------------------------------------------
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi
  USE mp,          ONLY : mp_bcast, mp_sum
  USE mp_images,   ONLY : my_image_id
  USE mp_world,    ONLY : world_comm
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE mp_pools,    ONLY : root_pool, me_pool, my_pool_id, inter_pool_comm, &
                          npool, intra_pool_comm
  USE scatter_mod, ONLY : gather_grid
  USE io_global,   ONLY : ionode_id, ionode, stdout
  USE io_files,    ONLY : tmp_dir, prefix, diropn, create_directory
  USE environment, ONLY : environment_start, environment_end
  USE ions_base,   ONLY : ntyp => nsp, nat, ityp
  USE cell_base,   ONLY : at, bg
  USE fft_base,    ONLY : dfftp
  USE lsda_mod,    ONLY : nspin
  USE gvect,       ONLY : ngm
  USE noncollin_module, ONLY : nspin_mag
  USE uspp,        ONLY : nlcc_any
  USE paw_variables, ONLY : okpaw
  USE eqv,         ONLY : vlocq
  USE qpoint,      ONLY : eigqts
  USE output,      ONLY : fildyn
  USE uspp_param,  ONLY : upf
  USE control_ph,  ONLY : tmp_dir_ph
  USE ph_restart,  ONLY : ph_readfile
  USE efield_mod,  ONLY : zstareu, zstarue, zstarue0, zstareu0, epsilon
  USE dvscf_interpolate, ONLY : dvscf_shift_center, dvscf_bare_calc, &
                                dvscf_long_range, multiply_iqr
  !
  IMPLICIT NONE
  !
  ! Input variables
  CHARACTER(LEN=256) :: outdir
  !! Directory containing input, output, and scratch files
  CHARACTER(80) :: verbosity
  !! if 'high', verbose output
  LOGICAL :: do_charge_neutral
  !! Renormalize dvscf to impose neutrality of Born effective charges
  LOGICAL :: do_long_range
  !! Subtract the long-range part of the potential before interpolation
  CHARACTER(LEN=256) :: fildvscf
  !! file where the potential variation is written
  CHARACTER(LEN=256) :: wpot_dir
  !! folder where w_pot binary files are written
  !
  ! --------------------------------------------------------------------
  CHARACTER(LEN=256) :: wpot_file
  !! filename of the w_pot binary file
  LOGICAL :: verbose
  !! if verbosity == 'high', set to .true.
  LOGICAL :: exst
  !! on output of diropn, exst is .True. if opened file already exists
  LOGICAL :: needwf = .FALSE.
  !! do not need to read wavefunction data
  LOGICAL :: xmldyn
  !! is fildyn xml format
  INTEGER :: iq, iq1, iq2, iq3, irc, irc1, irc2, irc3, imode, lrwpot, rest, &
             ierr, iat, idir, nt, is
  !! indices
  INTEGER :: nq1, nq2, nq3
  !! Shape of the coarse q grid where ph.x calculation is done
  INTEGER :: nqirr
  !! number of irreducible q points computed in ph.x
  INTEGER :: nqtot
  !! nqtot = nq1 * nq2 * nq3
  INTEGER :: ios
  !! io status
  !
  ! Pool parallelization of q points
  INTEGER :: nqlocal
  !! Numebr of q points to compute in given pool
  INTEGER :: nqbase
  !! Given pool reads nqbase + 1 to nqbase + nqlocal q points from dfile_dir
  INTEGER, ALLOCATABLE :: iq_l2g(:)
  !! Index map from local q points (1 ~ nqlocal) to global q points
  !
  LOGICAL :: shift_half(3)
  !! true when the center of the supercell is at (0.5).
  INTEGER :: rlatt(3)
  !! Real space unit cell index for w_pot
  INTEGER :: iun
  !! Unit for reading files
  INTEGER :: iunrlatt
  !! Unit for writing rlatt.txt file
  INTEGER :: iunwpot
  !! Unit for writing w_pot binary file
  REAL(DP) :: epsil(3,3)
  !! dynamical matrix, read from fildyn
  REAL(DP) :: arg, xq_cart(3), w_pot_sum(3), coeff, epsil_q, zeu_avg(3, 3)
  !!
  COMPLEX(DP) :: phase
  !!
  REAL(DP), ALLOCATABLE :: xqirr(:, :)
  !! (3, nqirr) irreducible q points computed in ph.x, in Cartesian corrdinate
  REAL(DP), ALLOCATABLE :: zeu(:,:,:)
  !! Born effective charge tensor, read from fildyn
  REAL(DP), ALLOCATABLE :: xqs_cry_global(:, :)
  !! (3, nqtot) coarse grid for phonon calculations, in crystal coordinate
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  ! (dfftp%nnr) auxiliary variable for multiplying exp(iqr)
  COMPLEX(DP), ALLOCATABLE :: dvscf(:, :, :, :)
  !! (dfftp%nnr, nspin_mag, 3*nat, nqtot) Total (bare + induced) perturbed
  !! potential.
  COMPLEX(DP), ALLOCATABLE :: w_pot(:, :, :)
  !! (dfftp%nnr, nspin_mag, 3*nat) Fourier transformed potential in real space.
  !! Computed at one unit cell index R at a time to save memory.
  COMPLEX(DP), ALLOCATABLE :: w_pot_gathered(:, :)
  !! (dfftp%nnr, nspin_mag) Temporary storage for gathered w_pot
  COMPLEX(DP), ALLOCATABLE :: dvscf_bare(:, :, :)
  !! (dfftp%nnr, nspin_mag, 3*nat) Temporary storage foir the bare potential.
  COMPLEX(DP), ALLOCATABLE :: dvscf_long(:, :, :)
  !! (dfftp%nnr, nspin_mag, 3*nat) Temporary storage foir the long-range potential.
  !
  LOGICAL, EXTERNAL :: has_xml
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  NAMELIST / input / prefix, outdir, wpot_dir, do_long_range, fildvscf, &
                     do_charge_neutral, fildyn, verbosity
  !
  CALL mp_startup()
  CALL environment_start('DVSCF_Q2R')
  !
  ! ---------------------------------------------------------------------------
  !
  ! Reading input arguments
  !
  ! Default values of input arguments
  !
  IF (ionode) CALL input_from_file()
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  do_charge_neutral = .FALSE.
  do_long_range = .FALSE.
  fildyn = ''
  wpot_dir = ' '
  fildvscf = 'dvscf'
  verbosity = 'default'
  !
  ! Read input file
  !
  IF (ionode) READ (5, input, IOSTAT=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('dvscf_q2r','error reading input namelist', ABS(ios))
  !
  ! Broadcast input arguments
  !
  CALL mp_bcast(prefix, ionode_id, world_comm)
  CALL mp_bcast(outdir, ionode_id, world_comm)
  CALL mp_bcast(do_charge_neutral, ionode_id, world_comm)
  CALL mp_bcast(do_long_range, ionode_id, world_comm)
  CALL mp_bcast(fildyn, ionode_id, world_comm)
  CALL mp_bcast(wpot_dir, ionode_id, world_comm)
  CALL mp_bcast(fildvscf, ionode_id, world_comm)
  CALL mp_bcast(verbosity, ionode_id, world_comm)
  !
  ! Check input arguments validity
  !
  IF (fildyn == '') CALL errore('dvscf_q2r', 'fildyn must be specified', 1)
  !
  IF (wpot_dir == ' ') wpot_dir = TRIM(outdir) // "/w_pot/"
  !
  tmp_dir = trimcheck(outdir)
  wpot_dir = trimcheck(wpot_dir)
  tmp_dir_ph = trimcheck(TRIM(tmp_dir) // '_ph' // int_to_char(my_image_id))
  xmldyn = has_xml(fildyn)
  !
  IF (TRIM(verbosity) == 'high') THEN
    verbose = .TRUE.
  ELSE
    verbose = .FALSE.
  ENDIF
  !
  ! Create output directory
  !
  IF (ionode) INQUIRE(FILE=TRIM(wpot_dir), EXIST=exst)
  CALL mp_bcast(exst, ionode_id, world_comm)
  IF (.NOT. exst) CALL create_directory(wpot_dir)
  !
  ! ---------------------------------------------------------------------------
  !
  ! Read xml data. Do not need wavefunction information.
  !
  needwf = .FALSE.
  CALL read_file_new(needwf)
  !
  IF (nspin_mag /= 1) CALL errore('dvscf_q2r', 'magnetism not implemented', 1)
  IF (nspin == 2) CALL errore('dvscf_q2r', 'LSDA magnetism not implemented', 1)
  IF (okpaw) CALL errore('dvscf_q2r', 'PAW not implemented', 1)
  !
  ! ---------------------------------------------------------------------------
  !
  ! Read nq1, nq2, nq3, and irreducible q points from fildyn0
  !
  IF (ionode) THEN
    !
    WRITE(stdout,'(/,4x," Reading grid info from file ",a)') TRIM(fildyn)//'0'
    !
    iun = find_free_unit()
    !
    OPEN(UNIT=iun, FILE=TRIM(fildyn)//'0', STATUS='old', FORM='formatted', &
         IOSTAT=ios)
    IF (ios /= 0) CALL errore('dvscf_q2r', 'problem opening fildyn0', ios)
    !
    READ(iun, *) nq1, nq2, nq3
    READ(iun, *) nqirr
    !
    ALLOCATE(xqirr(3, nqirr))
    DO iq = 1, nqirr
      READ(iun, *) xqirr(:, iq)
    ENDDO
    !
    CLOSE(UNIT=iun, STATUS='keep')
  ENDIF
  !
  CALL mp_bcast(nq1, ionode_id, world_comm)
  CALL mp_bcast(nq2, ionode_id, world_comm)
  CALL mp_bcast(nq3, ionode_id, world_comm)
  CALL mp_bcast(nqirr, ionode_id, world_comm)
  IF (.NOT. ionode) ALLOCATE(xqirr(3, nqirr))
  CALL mp_bcast(xqirr, ionode_id, world_comm)
  !
  nqtot = nq1 * nq2 * nq3
  !
  ! ---------------------------------------------------------------------------
  ! Read tensors.xml file for epsilon and Born effective charge
  !
  IF (do_charge_neutral .OR. do_long_range) THEN
    !
    ALLOCATE(zeu(3, 3, nat))
    !
    WRITE(stdout, *)
    WRITE(stdout, '(5x,a)') 'Reading epsilon and Born effective charge from file'
    !
    ! These arrays must be set to call ph_readfile
    !
    ALLOCATE(zstareu (3, 3,  nat))
    ALLOCATE(zstareu0(3, 3 * nat))
    ALLOCATE(zstarue (3 , nat, 3))
    ALLOCATE(zstarue0(3 * nat, 3))
    !
    CALL ph_readfile('tensors', 0, 0, ierr)
    !
    IF (ierr == 0) THEN
      zeu = zstareu
      epsil = epsilon
    ELSE
      CALL errore('dvscf_q2r', &
          'problem reading epsilon and zeu from tensors.xml', ierr)
    ENDIF
    !
    DEALLOCATE(zstareu)
    DEALLOCATE(zstareu0)
    DEALLOCATE(zstarue)
    DEALLOCATE(zstarue0)
    !
    ! Write zeu and epsil to file. To be read in dvscf_r2q
    !
    IF (ionode) THEN
      iun = find_free_unit()
      OPEN(iun, FILE=TRIM(wpot_dir)//'tensors.dat', FORM='formatted', &
          ACTION='write', IOSTAT=ios)
      IF (ios /= 0) CALL errore('dvscf_q2r', &
          'problem opening tensors.dat file for writing zeu and epsil', ios)
      WRITE(iun, *) '# dielectric constant epsil and Born effective charge zeu'
      WRITE(iun, *) epsil
      WRITE(iun, *) zeu
      CLOSE(iun, STATUS='KEEP')
    ENDIF
    !
    WRITE(stdout, '(5x,a)') 'Done reading epsilon and Born effective charge'
    !
  ENDIF
  !
  ! ---------------------------------------------------------------------------
  !
  ! Determine q point parallelization parameters
  ! Adapted from PW/divide_et_impera.f90
  !
  IF (npool == 1) THEN
    nqbase = 0
    nqlocal = nqtot
  ELSE
    nqlocal = nqtot / npool
    rest = nqtot - nqlocal * npool
    IF ( my_pool_id < rest ) nqlocal = nqlocal + 1
    ! nqbase: the position in the list of the first point that belongs to this
    !         pool, minus one
    nqbase = nqlocal * my_pool_id
    IF ( my_pool_id >= rest ) nqbase = nqbase + rest
  ENDIF
  !
  ! Define map of local q point index to global q point index
  !
  ALLOCATE(iq_l2g(nqlocal))
  DO iq = 1, nqlocal
    iq_l2g(iq) = iq + nqbase
  ENDDO
  !
  ! ---------------------------------------------------------------------------
  !
  ! Set the coarse q grid points
  !
  ALLOCATE(xqs_cry_global(3, nqtot))
  iq = 0
  DO iq3 = 0, nq3 - 1
    DO iq2 = 0, nq2 - 1
      DO iq1 = 0, nq1 - 1
        iq = iq + 1
        xqs_cry_global(1, iq) = REAL(iq1, DP) / REAL(nq1, DP)
        xqs_cry_global(2, iq) = REAL(iq2, DP) / REAL(nq2, DP)
        xqs_cry_global(3, iq) = REAL(iq3, DP) / REAL(nq3, DP)
      ENDDO
    ENDDO
  ENDDO
  !
  ! ---------------------------------------------------------------------------
  !
  ! Allocate dvscf
  !
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a,3I8)') "Allocating dvscf, shape", dfftp%nnr, nat*3, nqlocal
  ALLOCATE(dvscf(dfftp%nnr, nspin_mag, nat*3, nqlocal))
  dvscf = (0.d0, 0.d0)
  WRITE(stdout, '(5x,a)') "Done allocating dvscf"
  !
  ! ---------------------------------------------------------------------------
  ! Read dvscf files (induced part)
  !
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') "Reading and rotating dvscf files (induced part of dvscf)"
  CALL read_rotate_dvscf()
  !
  ! ---------------------------------------------------------------------------
  ! Renormalize induced part of dvscf for charge neutrality
  IF (do_charge_neutral) THEN
    WRITE(stdout, *)
    WRITE(stdout, '(5x,a)') "Renormalize induced part of dvscf for charge neutrality"
    zeu_avg = 0.d0
    DO iat = 1, nat
      zeu_avg = zeu_avg + zeu(:, :, iat)
    ENDDO
    zeu_avg = zeu_avg / REAL(nat, DP)
    !
    DO iq = 1, nqlocal
      xq_cart = xqs_cry_global(:, iq_l2g(iq))
      CALL cryst_to_cart(1, xq_cart, bg, +1)
      !
      ! Skip q = Gamma
      IF (SUM(ABS(xq_cart)) < 1.d-5) CYCLE
      !
      ! epsil_q = 1/q^2 * (q.T * epsil * q)
      epsil_q = 0.d0
      DO idir = 1, 3
        epsil_q = epsil_q + xq_cart(idir) * SUM(epsil(idir,:) * xq_cart(:))
      ENDDO
      epsil_q = epsil_q / SUM(xq_cart * xq_cart)
      !
      DO imode = 1, 3 * nat
        iat = (imode - 1) / 3 + 1
        idir = imode - 3 * (iat - 1)
        nt = ityp(iat)
        !
        ! Skip if xq_cart is orthonormal to atomic displacement
        IF (ABS(xq_cart(idir)) < 1.d-5) CYCLE
        !
        !! coeff = ( Z*q_i - sum_j (Z* - Z*avg)_{i,j} * q(j) / epsil_q )
        !!       / ( Z*q_i - sum_j (Z*)_{i,j} * q_j / epsil_q )
        coeff = ( upf(nt)%zp * xq_cart(idir)  &
          - SUM((zeu(idir,:,iat)-zeu_avg(idir,:)) * xq_cart(:)) / epsil_q ) &
              / ( upf(nt)%zp * xq_cart(idir)  &
          - SUM(zeu(idir,:,iat) * xq_cart(:)) / epsil_q )
        !
        dvscf(:, :, imode, iq) = dvscf(:, :, imode, iq) * coeff
        !
      ENDDO ! imode
    ENDDO ! iq
    WRITE(stdout, '(5x,a)') "Charge neutrality done"
  ENDIF ! do_charge_neutral
  !
  ! ---------------------------------------------------------------------------
  ! Compute bare part
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') "Computing bare part of dvscf"
  CALL start_clock('dvscf_bare')
  !
  ALLOCATE(dvscf_bare(dfftp%nnr, nspin_mag, nat*3))
  DO iq = 1, nqlocal
    WRITE (stdout, '(i8)', ADVANCE='no') iq
    IF(MOD(iq, 10) == 0 ) WRITE(stdout,*)
    FLUSH(stdout)
    !
    xq_cart = xqs_cry_global(:, iq_l2g(iq))
    CALL cryst_to_cart(1, xq_cart, bg, +1)
    !
    ! Need to do some initializations before computing dvscf_bare
    CALL init_phq_dvscf_q2r(xq_cart)
    !
    CALL dvscf_bare_calc(xq_cart, dvscf_bare, .FALSE.)
    dvscf(:,:,:,iq) = dvscf(:,:,:,iq) + dvscf_bare(:,:,:)
    !
    DEALLOCATE(eigqts)
    DEALLOCATE(vlocq)
    !
  ENDDO
  DEALLOCATE(dvscf_bare)
  CALL stop_clock('dvscf_bare')
  !
  ! ---------------------------------------------------------------------------
  ! Subtract long-range part (dipole potential) from
  ! When charge neutrality renormalization is done, zeu must be modified
  IF (do_long_range) THEN
    !
    WRITE(stdout, *)
    WRITE(stdout, '(5x,a)') "Computing and subtracting long-range part of dvscf"
    !
    IF (do_charge_neutral) THEN
      DO iat = 1, nat
        zeu(:,:,iat) = zeu(:,:,iat) - zeu_avg
      ENDDO
    ENDIF
    !
    ALLOCATE(dvscf_long(dfftp%nnr, nspin_mag, nat*3))
    !
    DO iq = 1, nqlocal
      xq_cart = xqs_cry_global(:, iq_l2g(iq))
      CALL cryst_to_cart(1, xq_cart, bg, +1)
      CALL dvscf_long_range(xq_cart, zeu, epsil, dvscf_long)
      !
      dvscf(:,:,:,iq) = dvscf(:,:,:,iq) - dvscf_long(:,:,:)
    ENDDO
    !
    DEALLOCATE(dvscf_long)
    !
    WRITE(stdout, '(5x,a)') "Long-range part done"
    !
  ENDIF ! do_long_range
  ! ---------------------------------------------------------------------------
  !
  ! Shift center of the potential from tau to the center of the supercell
  !
  shift_half(1) = ( MOD(nq1, 2) == 1 )
  shift_half(2) = ( MOD(nq2, 2) == 1 )
  shift_half(3) = ( MOD(nq3, 2) == 1 )
  !
  DO iq = 1, nqlocal
    xq_cart = xqs_cry_global(:, iq_l2g(iq))
    CALL cryst_to_cart(1, xq_cart, bg, +1)
    !
    CALL dvscf_shift_center(dvscf(:, :, :, iq), xq_cart, shift_half, -1)
  ENDDO
  !
  ! ---------------------------------------------------------------------------
  ! Multiply exp(iqr) to dvscf
  CALL start_clock('mul_iqr')
  !
  ALLOCATE(aux(dfftp%nnr))
  !
  DO iq = 1, nqlocal
    xq_cart = xqs_cry_global(:, iq_l2g(iq))
    CALL cryst_to_cart(1, xq_cart, bg, +1)
    !
    aux = (1.d0, 0.d0)
    CALL multiply_iqr(dfftp, xq_cart, aux)
    !
    DO imode = 1, 3 * nat
      DO is = 1, nspin_mag
        dvscf(:, is, imode, iq) = dvscf(:, is, imode, iq) * aux
      ENDDO
    ENDDO ! imode
    !
  ENDDO ! iq
  !
  DEALLOCATE(aux)
  !
  CALL stop_clock('mul_iqr')
  !
  ! ---------------------------------------------------------------------------
  ! Fourier transfrom dvscf to w_pot, write to file
  !
  WRITE(stdout, *)
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') "Fourier transforming dvscf to w_pot, writing to file"
  !
  IF (ionode) THEN
    iunrlatt = find_free_unit()
    OPEN(iunrlatt, FILE=TRIM(wpot_dir)//'rlatt.txt', FORM='formatted', &
      ACTION='write')
    WRITE(iunrlatt, '(a)') "# Real space unit cell index for w_pot"
    WRITE(iunrlatt, '(a)') "#     ir     ir1     ir2     ir3"
    WRITE(iunrlatt, '(I8)') nqtot
  ENDIF
  !
  iunwpot = find_free_unit()
  lrwpot = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  CALL start_clock('w_pot')
  !
  ALLOCATE(w_pot(dfftp%nnr, nspin_mag, nat*3))
  ALLOCATE(w_pot_gathered(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag))
  !
  irc = 0
  DO irc1 = -(nq1 - 1) / 2, nq1 / 2
    DO irc2 = -(nq2 - 1) / 2, nq2 / 2
      DO irc3 = -(nq3 - 1) / 2, nq3 / 2
        !
        irc = irc + 1
        !
        w_pot = (0.d0, 0.d0)
        rlatt = (/ irc1, irc2, irc3 /)
        IF (ionode) WRITE(iunrlatt, '(4I8)') irc, rlatt
        !
        ! Calculate w_pot by slow Fourier transform
        CALL start_clock('calc_w_pot')
        DO iq = 1, nqlocal
          arg = tpi * SUM(xqs_cry_global(:, iq_l2g(iq)) * REAL(rlatt(:), DP))
          phase = CMPLX(COS(arg), -SIN(arg), KIND=DP)
          !
          ! w_pot(:,:) = w_pot(:,:) + dvscf(:, :, iq) * phase
          CALL ZAXPY(dfftp%nnr * 3 * nat * nspin_mag, phase, &
              dvscf(1, 1, 1, iq), 1, w_pot(1, 1, 1), 1)
          !
        ENDDO ! iq
        !
        ! sum w_pot over pools (q points)
        CALL mp_sum(w_pot, inter_pool_comm)
        !
        w_pot = w_pot / REAL(nqtot, DP)
        !
        CALL stop_clock('calc_w_pot')
        !
        ! Write SUM(ABS(w_pot)) to stdout for debugging and testing
        !
        IF (verbose) THEN
          WRITE(stdout, '(5x,a,I8)') 'unit_cell_index ', irc
          WRITE(stdout, '(5x,a,3I12)') 'rlatt_crys ', rlatt
          WRITE(stdout, '(5x,a,3F12.4)') 'rlatt_cart ', MATMUL(at, REAL(rlatt, DP))
          !
          DO iat = 1, nat
            DO idir = 1, 3
              imode = 3 * (iat - 1) + idir
              w_pot_sum(idir) = SUM(ABS(w_pot(:, :, imode)))
            ENDDO
            !
            CALL mp_sum(w_pot_sum, intra_pool_comm)
            WRITE(stdout, '(5x,a,3F16.6)') "sum_w_pot", w_pot_sum
          ENDDO ! imode
          WRITE(stdout, *)
        ENDIF
        !
        ! Write w_pot to file
        !
        CALL start_clock('write_w_pot')
        !
        WRITE(wpot_file, '(a,I0,a)') 'wpot.irc', irc, '.dat'
        IF (ionode) THEN
          CALL diropn(iunwpot, TRIM(wpot_file), lrwpot, exst, TRIM(wpot_dir))
        ENDIF
        !
        DO imode = 1, 3 * nat
#if defined(__MPI)
          ! gather w_pot
          DO is = 1, nspin_mag
            CALL gather_grid(dfftp, w_pot(:, is, imode), w_pot_gathered(:, is))
          ENDDO ! ispin
          !
          ! write gathered file
          IF (ionode) THEN
            CALL davcio(w_pot_gathered, lrwpot, iunwpot, imode, +1)
          ENDIF
#else
          CALL davcio(w_pot(:, :, imode), lrwpot, iunwpot, imode, +1)
#endif
        ENDDO ! imode
        !
        IF (ionode) THEN
          CLOSE(UNIT=iunwpot, STATUS='keep')
        ENDIF
        !
        CALL stop_clock('write_w_pot')
        !
      ENDDO ! irc3
    ENDDO ! irc2
  ENDDO ! irc1
  !
  IF (ionode) CLOSE(iunrlatt, STATUS='keep')
  !
  DEALLOCATE(dvscf)
  DEALLOCATE(xqs_cry_global)
  DEALLOCATE(w_pot)
  DEALLOCATE(w_pot_gathered)
  IF (do_long_range .OR. do_charge_neutral) DEALLOCATE(zeu)
  !
  WRITE(stdout, *)
  WRITE(stdout, *)
  CALL print_clock('read_dvscf')
  CALL print_clock('dvscf_bare')
  CALL print_clock('mul_iqr')
  CALL print_clock('dvscf_shift')
  CALL print_clock('calc_w_pot')
  CALL print_clock('write_w_pot')
  !
  CALL environment_end('DVSCF_Q2R')
  CALL mp_global_end()
  !
CONTAINS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE read_rotate_dvscf()
  !----------------------------------------------------------------------------
  !!
  !! Read _ph0/prefix.dvscf files for irreducible q points, which are written
  !! by ph.x. Rotate them to local coarse q points and add to dvscf.
  !!
  !! Rotation of dvscf is adapted from dfile_star.f90
  !!
  !----------------------------------------------------------------------------
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi
  USE io_global,   ONLY : ionode, ionode_id
  USE io_files,    ONLY : prefix, diropn
  USE mp,          ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_pools,    ONLY : me_pool, inter_pool_comm, root_pool
  USE mp_images,   ONLY : intra_image_comm
  USE scatter_mod, ONLY : scatter_grid, gather_grid
  USE ions_base,   ONLY : nat, tau
  USE cell_base,   ONLY : at, bg
  USE symm_base,   ONLY : time_reversal, nsym, s, invs, ft, irt
  USE fft_base,    ONLY : dfftp
  USE spin_orb,    ONLY : domag
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE control_flags, ONLY : noinv
  USE control_ph,  ONLY : tmp_dir_ph
  USE cryst_ph,    ONLY : magnetic_sym
  USE ph_restart,  ONLY : ph_readfile
  USE modes,       ONLY : u, npert
  USE modes, ONLY : nirr, name_rap_mode, num_rap_mode
  USE qpoint,      ONLY : xq
  USE disp,        ONLY : nqs
  USE lr_symm_base,ONLY : nsymq, invsymq, minus_q
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: tmp_dir_phq
  !! directory where dvscf file for iqirr is located
  LOGICAL :: exst
  !! on output of diropn, exst is .True. if opened file exists
  LOGICAL :: imq_iq
  !! apply time reversal
  INTEGER :: ierr
  !! error variable
  INTEGER :: iudvscf
  !! unit for reading dvscf
  INTEGER :: i, j, k, ri, rj, rk, n, nn
  !! index of real space lattice
  INTEGER :: is
  !! index of spin
  INTEGER :: iq
  !! index of q point
  INTEGER :: imode, jmode
  !! index of mode
  INTEGER :: iqirr
  !! index of irreducible q point
  INTEGER :: iat, iat_rot
  !! index of atoms
  INTEGER :: ipol, jpol
  !! index of directions
  INTEGER :: isym, isym_inv
  !! index of symmetry operation
  INTEGER :: imq
  !! index of -q in the star (0 if not present)
  INTEGER :: isq(48)
  !! index of q in the star for a given sym
  INTEGER :: lrdrho
  !! the length of the deltarho files ( = length of dvscf files)
  INTEGER :: ftau(3,nsym)
  !! fractional translation in fft grid
  INTEGER :: s_scaled(3,3,nsym)
  !! scaled rotations
  REAL(DP) :: xq_tau
  !! xq dot tau phase factor
  REAL(DP) :: xqtmp(3)
  !! q point vector
  REAL(DP) :: xq_diff(3)
  !! q point sxq - xqs_cry
  REAL(DP) :: sxq(3, 48)
  !! list of vectors in the star of q
  COMPLEX(DP) :: phase
  !! phase factor
  INTEGER, ALLOCATABLE :: nq_iqirr(:)
  !! (nqirr) number of local q points found for each iqirr
  INTEGER, ALLOCATABLE :: nqs_save(:)
  !! (nqirr) nqs for each iqirr
  INTEGER, ALLOCATABLE :: isym_iqloc(:)
  !! (nqlocal) isq index for local q point. -1 if not in the star
  LOGICAL, ALLOCATABLE :: imq_iqloc(:)
  !! (nqlocal) .TRUE. if symmetry index for local q point have time reversal
  REAL(DP), ALLOCATABLE :: g_residual(:, :)
  !! (3, nqlocal) sxq = xqs_cry_global + g_residual
  COMPLEX(DP), ALLOCATABLE :: rotmat(:, :)
  !! matrix for rotating pattern basis to crystal basis
  COMPLEX(DP), ALLOCATABLE :: dvscf_p_gathered(:)
  !! (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) dvscf in gathered form
  COMPLEX(DP), ALLOCATABLE :: dvscf_p_crys(:, :, :)
  !! (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag, 3*nat)
  !! dvscf for iqirr in crystal basis
  COMPLEX(DP), ALLOCATABLE :: dvscf_p_rotated(:, :, :)
  !! (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag, 3*nat)
  !! dvscf rotated from iqirr to iq
  COMPLEX(DP), ALLOCATABLE :: dvscf_p(:, :)
  ! (dfftp%nnr, nspin_mag) dvscf_p_gathered scattered to procs
  COMPLEX(DP), ALLOCATABLE :: dvscf_iq(:, :, :)
  ! (dfftp%nnr, nspin_mag, 3*nat) Temporary storage of read dvscf data
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  ! (dfftp%nnr) auxiliary variable for multiplying exp(iGr)
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  CALL start_clock('read_dvscf')
  !
  IF (me_pool == root_pool) THEN
    ALLOCATE(dvscf_p_crys(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag, 3*nat))
    ALLOCATE(dvscf_p_rotated(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag, 3*nat))
  ENDIF
  ALLOCATE(dvscf_p_gathered(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  ALLOCATE(aux(dfftp%nnr))
  ALLOCATE(dvscf_p(dfftp%nnr, nspin_mag))
  ALLOCATE(dvscf_iq(dfftp%nnr, nspin_mag, 3*nat))
  ALLOCATE(rotmat(3*nat, 3*nat))
  ALLOCATE(g_residual(3, nqlocal))
  ALLOCATE(isym_iqloc(nqlocal))
  ALLOCATE(imq_iqloc(nqlocal))
  ALLOCATE(nq_iqirr(nqirr))
  ALLOCATE(nqs_save(nqirr))
  ALLOCATE(npert(3 * nat))
  ALLOCATE(num_rap_mode(3 * nat))
  ALLOCATE(name_rap_mode(3 * nat))
  ALLOCATE(u(3 * nat, 3 * nat))
  !
  nq_iqirr = 0
  nqs_save = 0
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  DO iqirr = 1, nqirr
    !
    CALL mp_barrier(intra_image_comm)
    !
    CALL ph_readfile('data_u', iqirr, 0, ierr)
    !
    xq = xqirr(:, iqirr)
    !
    ! Rotation matrix to rotate dvscf from pattern (u) to crystal basis
    ! (Adapted from dfile_star.f90)
    !
    rotmat = (0.d0, 0.d0)
    DO iat = 1, nat
      DO ipol = 1, 3
        DO jpol = 1, 3
          imode = 3 * (iat - 1) + ipol
          jmode = 3 * (iat - 1) + jpol
          rotmat(jmode, :) = rotmat(jmode, :) + at(ipol, jpol) * CONJG(u(imode, :))
        ENDDO
      ENDDO
    ENDDO
    !
    ! Setup small group of q symmetry
    !
    CALL set_small_group_of_q(nsymq, invsymq, minus_q)
    !
    ! Rotate irreducible q points
    !
    CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, verbose)
    !
    nqs_save(iqirr) = nqs
    IF (time_reversal .AND. imq == 0) nqs_save(iqirr) = nqs * 2
    !
    isym_iqloc(:) = -1
    imq_iqloc(:) = .FALSE.
    !
    ! Find the coarse q points in xqs_cry_global from the star of xqirr
    ! sxq(:, isq(isym)) = xqs_cry + g_residual
    !
    ! If using time reversal symmetry, also find
    ! -sxq(:, isq(isym)) = xqs_cry + g_residual
    !
    DO iq = 1, nqlocal
      !
      ! check sxq = xqs_cry + g_residual
      !
      DO isym = 1, nsym
        xq_diff = sxq(:, isq(isym))
        CALL cryst_to_cart(1, xq_diff, at, -1)
        xq_diff = xq_diff - xqs_cry_global(:, iq + nqbase)
        !
        IF ( ALL( ABS(xq_diff(:) - NINT(xq_diff(:))) < 1.d-5 ) ) THEN
          nq_iqirr(iqirr) = nq_iqirr(iqirr) + 1
          g_residual(:, iq) = xq_diff(:)
          isym_iqloc(iq) = isym
          imq_iqloc(iq) = .FALSE.
          EXIT
        ENDIF
        !
      ENDDO
      !
      ! if not found, check time reversal -sxq = xqs_cry + g_residual
      !
      IF (isym_iqloc(iq) == -1 .AND. time_reversal .AND. imq == 0) THEN
        DO isym = 1, nsym
          xq_diff = -sxq(:, isq(isym))
          CALL cryst_to_cart(1, xq_diff, at, -1)
          xq_diff = xq_diff - xqs_cry_global(:, iq + nqbase)
          !
          IF ( ALL( ABS(xq_diff(:) - NINT(xq_diff(:))) < 1.d-5 ) ) THEN
            nq_iqirr(iqirr) = nq_iqirr(iqirr) + 1
            g_residual(:, iq) = xq_diff(:)
            isym_iqloc(iq) = isym
            imq_iqloc(iq) = .TRUE.
            EXIT
          ENDIF
          !
        ENDDO
      ENDIF
      !
    ENDDO ! iq
    !
    ! Open dvscf file
    !
    IF (iqirr == 1) THEN
      tmp_dir_phq = TRIM(tmp_dir_ph)
    ELSE
      tmp_dir_phq = TRIM(tmp_dir_ph) // TRIM(prefix) // '.q_' &
                  & // TRIM(int_to_char(iqirr)) // '/'
    ENDIF
    !
    iudvscf = find_free_unit()
    IF (ionode) THEN
      CALL diropn(iudvscf, fildvscf, lrdrho, exst, tmp_dir_phq)
    ENDIF
    CALL mp_bcast(exst, ionode_id, intra_image_comm)
    IF (.NOT. exst) CALL errore('read_rotate_dvscf', &
        'dvscf file does not exist', iqirr)
    !
    ! Read dvscf file ad rotate to crystal coordinate
    ! (Adapted from dfile_star.f90)
    !
    IF (me_pool == root_pool) dvscf_p_crys = (0.d0, 0.d0)
    !
    DO imode = 1, 3 * nat
      !
      ! read dvscf file to dvscf_p, and gather back to dvscf_p_gathered
      !
      CALL davcio_drho(dvscf_p, lrdrho, iudvscf, imode, -1)
      !
      DO is = 1, nspin_mag
#if defined(__MPI)
        CALL gather_grid(dfftp, dvscf_p(:, is), dvscf_p_gathered)
#else
        dvscf_p_gathered = dvscf_p(:, is)
#endif
        !
        ! Rotate from pattern to crystal coordiate
        !
        IF (me_pool == root_pool) THEN
          DO jmode = 1, 3 * nat
            dvscf_p_crys(:, is, jmode) = dvscf_p_crys(:, is, jmode) &
                + rotmat(jmode, imode) * dvscf_p_gathered(:)
          ENDDO
        ENDIF
      ENDDO
      !
    ENDDO ! imode
    !
    IF (ionode) CLOSE( UNIT = iudvscf, STATUS = 'KEEP' )
    !
    ! take away the phase due to the q-point
    !
    CALL scale_sym_ops( nsym, s, ft, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
         s_scaled, ftau )
    !
    IF (me_pool == root_pool) THEN
      DO iat = 1, nat
        !
        xq_tau = tpi * SUM(xq * tau(:, iat))
        phase = CMPLX(COS(xq_tau), SIN(xq_tau), KIND=DP)
        !
        DO ipol = 1, 3
           imode = (iat - 1) * 3 + ipol
           dvscf_p_crys(:, :, imode) = phase * dvscf_p_crys(:, :, imode)
        ENDDO
      ENDDO
    ENDIF ! root_pool
    !
    ! If no q point belongs to this star, skip reading dvscf
    !
    IF ( ALL(isym_iqloc == -1) ) CYCLE
    !
    ! Rotate dvscf from iqirr to iq (Adapted from dfile_star.f90)
    !
    DO iq = 1, nqlocal
      dvscf_iq = (0.d0, 0.d0)
      !
      IF (isym_iqloc(iq) == -1) CYCLE
      !
      isym = isym_iqloc(iq)
      isym_inv = invs(isym)
      imq_iq = imq_iqloc(iq)
      !
      IF (me_pool == root_pool) THEN
        dvscf_p_rotated = (0.d0, 0.d0)
        !
        DO is = 1, nspin_mag
          KLOOP : DO k = 1, dfftp%nr3
            JLOOP : DO j = 1, dfftp%nr2
              ILOOP : DO i = 1, dfftp%nr1
                !
                ! Here I rotate r
                !
                 CALL rotate_grid_point(s_scaled(1,1,isym_inv),ftau(1,isym_inv),&
                      i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                !
                n  = (i-1)  + (j-1)*dfftp%nr1  + (k-1)*dfftp%nr2*dfftp%nr1  + 1
                nn = (ri-1) + (rj-1)*dfftp%nr1 + (rk-1)*dfftp%nr2*dfftp%nr1 + 1
                !
                DO iat = 1, nat
                  iat_rot = irt(isym_inv, iat)
                  jmode = (iat_rot - 1) * 3
                  !
                  DO ipol = 1, 3
                      imode = (iat - 1) * 3 + ipol
                      !
                      dvscf_p_rotated(n, is, imode) = dvscf_p_rotated(n, is, imode) &
                        + ( s(ipol, 1, isym_inv) * dvscf_p_crys(nn, is, jmode + 1) + &
                            s(ipol, 2, isym_inv) * dvscf_p_crys(nn, is, jmode + 2) + &
                            s(ipol, 3, isym_inv) * dvscf_p_crys(nn, is, jmode + 3) )
                      !
                  ENDDO
                ENDDO
                !
              ENDDO ILOOP
            ENDDO JLOOP
          ENDDO KLOOP
          !
        ENDDO ! ispin
        !
        ! Take complex conjugate if using time reversal symmetry
        !
        IF (imq_iq) THEN
          dvscf_p_rotated = CONJG(dvscf_p_rotated)
          ! FIXME: magnetism
        ENDIF
        !
      ENDIF ! root_pool
      !
      DO is = 1, nspin_mag
        DO imode = 1, 3 * nat
          !
          ! Scatter dvscf_p_rotated to dvscf_p. dvscf_p_rotated is allocated
          ! only at root_pool, so we copy it to dvscf_p_gathered.
          !
#if defined(__MPI)
          dvscf_p_gathered = (0.d0, 0.d0)
          IF (me_pool == root_pool) THEN
            dvscf_p_gathered = dvscf_p_rotated(:, is, imode)
          ENDIF
          !
          CALL scatter_grid(dfftp, dvscf_p_gathered, dvscf_iq(:, is, imode))
#else
          dvscf_iq(1:dfftp%nnr, is, imode) = dvscf_p_rotated(1:dfftp%nnr, is, imode)
#endif
          !
        ENDDO ! imode
      ENDDO ! ispin
      !
      ! Add back the phase factor for the new q-point
      !
      xqtmp = sxq(:, isq(isym))
      if (imq_iq) xqtmp = -xqtmp
      !
      DO iat = 1, nat
        xq_tau = tpi * SUM(xqtmp * tau(:, iat))
        phase = CMPLX(COS(xq_tau), -SIN(xq_tau), KIND=DP)
        !
        DO ipol = 1, 3
          imode = (iat - 1) * 3 + ipol
          dvscf_iq(:, :, imode) = dvscf_iq(:, :, imode) * phase
        ENDDO
      ENDDO
      !
      ! Rotate from crystal to cartesian coordinates
      !
      DO is = 1, nspin_mag
        DO iat = 1, nat
          imode = (iat - 1) * 3
          DO ipol = 1, 3
            dvscf(:, is, imode + ipol, iq) = dvscf_iq(:, is, imode + 1) * bg(ipol, 1) &
                                           + dvscf_iq(:, is, imode + 2) * bg(ipol, 2) &
                                           + dvscf_iq(:, is, imode + 3) * bg(ipol, 3)
          ENDDO
        ENDDO
      ENDDO
      !
      ! Multiply phase factor exp(iGr) if necessary
      !
      IF ( ANY(ABS(g_residual(:, iq)) > 1.d-5) ) THEN
        xqtmp = g_residual(:, iq)
        CALL cryst_to_cart(1, xqtmp, bg, +1)
        !
        aux = (1.d0, 0.d0)
        CALL multiply_iqr(dfftp, xqtmp, aux)
        !
        DO is = 1, nspin_mag
          DO imode = 1, 3 * nat
            dvscf(:, is, imode, iq) = dvscf(:, is, imode, iq) * aux(:)
          ENDDO
        ENDDO
      ENDIF ! g_residual
      !
    ENDDO ! iq
    !
  ENDDO ! iqirr
  !
  ! Check whether all q points are found.
  !
  CALL mp_sum(nq_iqirr, inter_pool_comm)
  IF (ionode) THEN
    DO iqirr = 1, nqirr
      IF (nqs_save(iqirr) /= nq_iqirr(iqirr)) THEN
        CALL errore('read_rotate_dvscf', 'problem finding q points from star', 1)
      ENDIF
    ENDDO
  ENDIF
  !
  IF (SUM(nq_iqirr) < nqtot) CALL errore('read_rotate_dvscf', &
      'Not all coarse q points are found in the star', 1)
  IF (SUM(nq_iqirr) > nqtot) CALL errore('read_rotate_dvscf', &
      'More than nq1*nq2*nq3 q points are found in the star', 1)
  !
  IF (me_pool == root_pool) THEN
    DEALLOCATE(dvscf_p_crys)
    DEALLOCATE(dvscf_p_rotated)
  ENDIF
  DEALLOCATE(dvscf_p_gathered)
  DEALLOCATE(aux)
  DEALLOCATE(dvscf_p)
  DEALLOCATE(dvscf_iq)
  DEALLOCATE(rotmat)
  DEALLOCATE(g_residual)
  DEALLOCATE(isym_iqloc)
  DEALLOCATE(imq_iqloc)
  DEALLOCATE(nq_iqirr)
  DEALLOCATE(nqs_save)
  DEALLOCATE(npert)
  DEALLOCATE(num_rap_mode)
  DEALLOCATE(name_rap_mode)
  DEALLOCATE(u)
  !
  CALL stop_clock('read_dvscf')
  !
!------------------------------------------------------------------------------
END SUBROUTINE read_rotate_dvscf
!------------------------------------------------------------------------------
SUBROUTINE init_phq_dvscf_q2r(xq)
  !----------------------------------------------------------------------------
  !!
  !! Do the initializations needed in dvscf_bare_calc
  !!
  !----------------------------------------------------------------------------
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi
  USE cell_base,   ONLY : tpiba2, omega
  USE ions_base,   ONLY : ntyp => nsp, nat, tau
  USE eqv,         ONLY : vlocq
  USE qpoint,      ONLY : eigqts
  USE atom,        ONLY : rgrid, msh
  USE gvect,       ONLY : g, ngm
  USE uspp,        ONLY : nlcc_any
  USE uspp_param,  ONLY : upf
  USE m_gth,       ONLY : setlocq_gth
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xq(3)
  !! Input: q point (in cartesian coordinate)
  !
  INTEGER :: na, nt
  REAL(DP) :: arg
  !
  ! First, calculate vlocq
  ! Adapted from PH/phq_init.f90, step 0, a and b
  !
  ALLOCATE(eigqts(nat))
  !
  DO na = 1, nat
    !
    arg = ( xq(1) * tau(1,na) + &
            xq(2) * tau(2,na) + &
            xq(3) * tau(3,na) ) * tpi
    !
    eigqts(na) = CMPLX(COS(arg), -SIN(arg), KIND=DP)
    !
  END DO
  !
  ! ... b) the fourier components of the local potential at q+G
  !
  ALLOCATE(vlocq(ngm , ntyp))
  vlocq(:,:) = 0.D0
  !
  DO nt = 1, ntyp
     !
     IF (upf(nt)%tcoulombp) THEN
        CALL setlocq_coul( xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt))
     ELSE IF (upf(nt)%is_gth) THEN
        CALL setlocq_gth( nt, xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
     ELSE
        CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, &
                   vlocq(1,nt) )
     ENDIF
     !
  END DO
  !
!------------------------------------------------------------------------------
END SUBROUTINE init_phq_dvscf_q2r
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
END PROGRAM dvscf_q2r
!------------------------------------------------------------------------------
