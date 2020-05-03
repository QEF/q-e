!
! Copyright (C) 2020-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE dvscf_interpolate
  !----------------------------------------------------------------------------
  !!
  !! Module for Fourier interpolation of phonon potential dvscf.
  !! Uses the output of dvscf_q2r.x program to compute the induced part of
  !! dvscf at each q point.
  !!
  !! See header of dvscf_q2r.f90 for details.
  !!
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  !
  IMPLICIT NONE
  !
  ! Input parameters
  !
  LOGICAL, SAVE :: ldvscf_interpolate
  !! If .true., use Fourier interpolation of phonon potential to compute the
  !! induced part of phonon potential at each q point
  LOGICAL, SAVE :: do_long_range
  !! If .true., add the long-range part of the potential to dvscf
  LOGICAL, SAVE :: do_charge_neutral
  !! If .true., impose the neutrality condition on Born effective charges
  CHARACTER(LEN=256), SAVE :: wpot_dir
  !! folder where w_pot binary files are located
  !
  ! Local variables
  !
  LOGICAL, SAVE :: shift_half(3)
  !! true when the center of the supercell is at (0.5).
  INTEGER, SAVE :: nrtot
  !! Total number of unit cell grid points. Read from rlatt.txt
  INTEGER, SAVE :: nrlocal
  !! local number of R points to read in the given pool
  INTEGER, SAVE :: nrbase
  !! Each pool read R points ir = nrbase + 1, ..., nrbase + nrlocal
  INTEGER, SAVE :: lrwpot
  !! record length of wpot file
  INTEGER, ALLOCATABLE, SAVE :: iunwpot(:)
  !! units for reading w_pot file. The units are opened at start_q, and
  !! closed at last_q.
  INTEGER, ALLOCATABLE, SAVE :: rlatt(:, :)
  !! (3, nrtot) Unit cell grid indices. Read from rlatt.txt
  REAL(DP), SAVE :: epsilon_r2q(3,3)
  !! Dielectric matrix. Read from tensors.dat.
  REAL(DP), ALLOCATABLE, SAVE :: zeu_r2q(:, :, :)
  !! Born effective charge tensor. Read from tensors.dat.
  !
  CONTAINS
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_interpol_setup()
    !--------------------------------------------------------------------------
    !! Do setups for dvscf_r2q. Called in phq_setup at the beginning of each
    !! q point calculation.
    !--------------------------------------------------------------------------
    !
    USE kinds,       ONLY : DP
    USE mp,          ONLY : mp_bcast
    USE mp_pools,    ONLY : npool, my_pool_id
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : ionode_id, ionode, stdout
    USE io_files,    ONLY : prefix
    USE fft_base,    ONLY : dfftp
    USE ions_base,   ONLY : nat
    USE lsda_mod,    ONLY : nspin
    USE noncollin_module, ONLY : nspin_mag
    USE paw_variables,    ONLY : okpaw
    USE Coul_cut_2D,      ONLY : do_cutoff_2D
    USE ldaU,        ONLY : lda_plus_u
    !
    IMPLICIT NONE
    !
    INTEGER :: iun
    !! Index for reading tensors.dat and rlatt.txt files
    INTEGER :: i, j
    !! Index for directions
    INTEGER :: iatm
    !! Atom index
    INTEGER :: ir, ir_, ir1, ir2, ir3
    !! Real space unit cell index
    INTEGER :: nr
    !! Real space unit cell size
    INTEGER :: rest
    !! variable used for calculating nrlocal
    INTEGER :: ios
    !! iostat
    REAL(DP) :: zeu_r2q_avg(3,3)
    !! Average of Born effective charge. Used for charge neutrality correction.
    INTEGER, EXTERNAL :: find_free_unit
    !
    IF (nspin_mag /= 1) CALL errore(' dvscf_r2q', ' magnetism not implemented',1)
    IF (nspin == 2) CALL errore(' dvscf_r2q', ' LSDA magnetism not implemented',1)
    IF (okpaw) CALL errore(' dvscf_r2q', ' PAW not implemented',1)
    IF (do_cutoff_2D) CALL errore(' dvscf_r2q', ' 2D Coulomb cutoff not implemented',1)
    IF (lda_plus_u) CALL errore(' dvscf_r2q', ' lda_plus_u not implemented',1)
    !
    CALL start_clock('dvscf_setup')
    !
    IF (do_long_range) THEN
      !
      ! Read tensors.dat file for epsil and zeu
      !
      ALLOCATE(zeu_r2q(3, 3, nat))
      !
      IF (ionode) THEN
        iun = find_free_unit()
        OPEN(iun, FILE=TRIM(wpot_dir)//'tensors.dat', FORM='formatted', &
          STATUS='old', ACTION='read', IOSTAT=ios)
        IF (ios /= 0) CALL errore('dvscf_interpol_setup', &
            'problem opening tensors.dat', ios)
        !
        READ(iun, *, IOSTAT=ios)
        IF (ios /= 0) CALL errore('dvscf_interpol_setup', &
            'problem reading tensors.dat', ios)
        !
        READ(iun, *, IOSTAT=ios) epsilon_r2q
        IF (ios /= 0) CALL errore('dvscf_interpol_setup', &
            'problem reading epsil from tensors.dat', ios)
        !
        READ(iun, *, IOSTAT=ios) zeu_r2q
        IF (ios /= 0) CALL errore('dvscf_interpol_setup', &
            'problem reading zeu from tensors.dat', ios)
        !
        CLOSE(iun, STATUS='keep')
        !
        WRITE(stdout, *) ' Dielectric matrix epsilon'
        WRITE(stdout, '(3f10.5)') ( (epsilon_r2q(i,j), j=1,3), i=1,3 )
        WRITE(stdout, *) ' Born effective charge zeu'
        DO iatm = 1, nat
          WRITE(stdout, * ) ' na= ', iatm
          WRITE(stdout,'(3f10.5)') ( (zeu_r2q(i,j,iatm), j=1,3), i=1,3 )
        END DO
        !
      ENDIF
      !
      CALL mp_bcast(epsilon_r2q, ionode_id, intra_image_comm)
      CALL mp_bcast(zeu_r2q, ionode_id, intra_image_comm)
      !
      IF (do_charge_neutral) THEN
        ! Renormalize Born effective charge
        zeu_r2q_avg = 0.0
        DO iatm = 1, nat
          zeu_r2q_avg = zeu_r2q_avg + zeu_r2q(:,:, iatm)
        ENDDO
        zeu_r2q_avg = zeu_r2q_avg / nat
        !
        DO iatm = 1, nat
          zeu_r2q(:,:,iatm) = zeu_r2q(:,:,iatm) - zeu_r2q_avg
        ENDDO
        !
        WRITE(stdout, *) ' Born effective charge after charge neutrality correction'
        DO iatm = 1, nat
          WRITE(stdout, * ) ' na= ', iatm
          WRITE(stdout,'(3f10.5)') ( (zeu_r2q(i,j,iatm), j=1,3), i=1,3 )
        END DO
        !
      ENDIF ! do_charge_neutral
      !
    ENDIF ! do_long_range
    !
    ! read rlatt.txt file to parse real space unit cell grid
    !
    IF (ionode) THEN
      iun = find_free_unit()
      OPEN(iun, FILE=TRIM(wpot_dir)//'rlatt.txt', FORM='formatted', &
        STATUS='old', ACTION='read', IOSTAT=ios)
      IF (ios /= 0) CALL errore('dvscf_interpol_setup', &
          'problem reading rlatt.txt', ios)
      !
      READ(iun, *)
      READ(iun, *)
      READ(iun, *) nrtot
      !
      ALLOCATE(rlatt(3, nrtot))
      DO ir = 1, nrtot
        READ(iun, *) ir_, ir1, ir2, ir3
        rlatt(:, ir) = (/ ir1, ir2, ir3 /)
      ENDDO
      !
      CLOSE(iun, STATUS='keep')
    ENDIF
    !
    CALL mp_bcast(nrtot, ionode_id, intra_image_comm)
    IF (.NOT. ionode) ALLOCATE(rlatt(3, nrtot))
    CALL mp_bcast(rlatt, ionode_id, intra_image_comm)
    !
    ! We use pools to distribute R points. Determine which R points to read.
    ! Adapted from PW/divide_et_impera.f90
    !
    IF (npool == 1) THEN
      nrbase = 0
      nrlocal = nrtot
    ELSE
      nrlocal = nrtot / npool
      rest = nrtot - nrlocal * npool
      IF ( my_pool_id < rest ) nrlocal = nrlocal + 1
      !
      ! nrbase: the position in the list of the first point that belongs to
      !         this pool, minus one
      nrbase = nrlocal * my_pool_id
      IF ( my_pool_id >= rest ) nrbase = nrbase + rest
      !
    ENDIF ! npool
    !
    ! Set shift_half(ipol) == .true. if number of supercell along direction
    ! ipol is odd, because the center of the supercell is at (0.5) not (0.0).
    !
    DO i = 1, 3
      nr = MAXVAL(rlatt(i, :)) - MINVAL(rlatt(i, :)) + 1
      shift_half(i) = ( MOD(nr, 2) == 1 )
    ENDDO
    !
    ! iunwpot is allocated here, and its value is set in openfilq
    !
    ALLOCATE(iunwpot(nrlocal))
    !
    CALL stop_clock('dvscf_setup')
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_interpol_setup
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_r2q(xq, u_in, dvscf)
    !--------------------------------------------------------------------------
    !!
    !! Read inverse Fourier transformed potential w_pot (written by
    !! dvscf_q2r.x) and Fourier transform to compute dvscf at given q point.
    !!
    !! Originally proposed by [1], long-range part described in [2].
    !! [1] Eiguren and Ambrosch-Draxl, PRB 78, 045124 (2008)
    !! [2] Xavier Gonze et al, Comput. Phys. Commun., 107042 (2019)
    !!
    !! dvscf(r,q) = exp(-iqr) (dvlong(r,q) + sum_R exp(iqR) w_pot(r,R))
    !!
    !! In this subroutine, pool parallelization is used to distribute R points
    !! to nodes so that the root of each pool reads different w_pot file
    !! simultaneously, reducing the io time.
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,       ONLY : DP
    USE constants,   ONLY : tpi
    USE mp,          ONLY : mp_barrier, mp_sum
    USE mp_pools,    ONLY : root_pool, me_pool, inter_pool_comm
    USE mp_images,   ONLY : intra_image_comm
    USE io_global,   ONLY : ionode, stdout
    USE scatter_mod, ONLY : scatter_grid
    USE fft_base,    ONLY : dfftp
    USE cell_base,   ONLY : at
    USE modes,       ONLY : nmodes
    USE noncollin_module, ONLY : nspin_mag
    USE control_ph,  ONLY : current_iq, start_q, last_q
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xq(3)
    !! Input. q point to compute dvscf, in Cartesian coordinate
    COMPLEX(DP), INTENT(IN) :: u_in(nmodes, nmodes)
    !! Input. Basis of atomic displacement
    COMPLEX(DP) :: dvscf(dfftp%nnr, nspin_mag, nmodes)
    !! Output. Interpolated dvscf
    !
    INTEGER :: imode, jmode, ir, irlocal, is
    REAL(DP) :: xq_cry(3), arg
    COMPLEX(DP) :: phase
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !! (dfftp%nnr) temporary variable for multiplying exp(-iqr)
    COMPLEX(DP), ALLOCATABLE :: dvscf_cart(:, :, :)
    !! (dfftp%nnr, nmodes) dvscf in the Cartesian basis
    COMPLEX(DP), ALLOCATABLE :: dvscf_bare(:, :, :)
    !! (dfftp%nnr, nspin_mag, nmodes) bare part of dvscf. Computed
    !! and subtracted from dvscf
    COMPLEX(DP), ALLOCATABLE :: dvscf_long(:, :, :)
    !! (dfftp%nnr, nspin_mag, nmodes) long-range part of dvscf. Computed
    !! and added to dvscf
    COMPLEX(DP), ALLOCATABLE :: w_pot(:)
    !! (dfftp%nnr) w_pot in the Cartesian basis for each imode
    COMPLEX(DP), ALLOCATABLE :: w_pot_gathered(:, :)
    !! temporary storage of gathered w_pot
    !
    LOGICAL, EXTERNAL :: has_xml
    INTEGER, EXTERNAL :: find_free_unit
    !
    CALL mp_barrier(intra_image_comm)
    !
    CALL start_clock('dvscf_r2q')
    !
    xq_cry = xq
    CALL cryst_to_cart (1, xq_cry, at, -1)
    !
    ALLOCATE(w_pot_gathered(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin_mag))
    ALLOCATE(dvscf_cart(dfftp%nnr, nspin_mag, nmodes))
    ALLOCATE(w_pot(dfftp%nnr))
    !
    dvscf_cart = (0.d0, 0.d0)
    !
    !--------------------------------------------------------------------------
    !
    DO irlocal = 1, nrlocal
      ir = irlocal + nrbase
      !
      ! Phase factor for Fourier transformation exp(iqR)
      arg = tpi * SUM(xq_cry(1:3) * REAL(rlatt(1:3, ir), DP))
      phase = CMPLX(COS(arg), SIN(arg), KIND=DP)
      !
      !------------------------------------------------------------------------
      ! Read w_pot file
      !
      DO imode = 1, nmodes
        !
        CALL start_clock('dvscf_davcio')
        !
        w_pot = (0.d0, 0.d0)
        !
        IF (me_pool == root_pool) THEN
          ! read dvscf file to dvscf_p_gathered
          CALL davcio(w_pot_gathered, lrwpot, iunwpot(irlocal), imode, -1)
        ENDIF
        CALL stop_clock('dvscf_davcio')
        !
        DO is = 1, nspin_mag
          !
#if defined (__MPI)
          ! scatter to w_pot
          CALL start_clock('dvscf_scatgrid')
          CALL scatter_grid(dfftp, w_pot_gathered(:, is), w_pot)
          CALL stop_clock('dvscf_scatgrid')
#else
          w_pot(1:dfftp%nnr) = w_pot_gathered(1:dfftp%nnr, is)
#endif
          !
          ! Compute Fourier transform
          ! dvscf_cart(r,q) = exp(-iqr) (dvlong(r,q) + sum_R exp(iqR) w_pot(r,R))
          ! exp(-iqr) part is multiplied after the ir loop.
          !
          CALL start_clock('dvscf_fourier')
          CALL ZAXPY(dfftp%nnr, phase, w_pot(1), 1, dvscf_cart(1, is, imode), 1)
          CALL stop_clock('dvscf_fourier')
          !
        ENDDO ! ispin
        !
      ENDDO ! imode
      !
    ENDDO ! ir
    !
    ! sum dvscf over all pools as R points are distributed
    CALL start_clock('dvscf_mpsum')
    CALL mp_sum(dvscf_cart, inter_pool_comm)
    CALL stop_clock('dvscf_mpsum')
    !
    ! Multiply exp(-iqr) to dvscf
    !
    CALL start_clock('dvscf_iqr')
    !
    ALLOCATE(aux(dfftp%nnr))
    aux = (1.d0, 0.d0)
    CALL multiply_iqr(dfftp, -xq, aux)
    DO imode = 1, nmodes
      DO is = 1, nspin_mag
        dvscf_cart(:, is, imode) = dvscf_cart(:, is, imode) * aux(:)
      ENDDO
    ENDDO
    DEALLOCATE(aux)
    !
    CALL stop_clock('dvscf_iqr')
    !
    !--------------------------------------------------------------------------
    !
    ! In dvscf_q2r.x, center of dvscf is shifted from tau to (0,0,0). Here we
    ! shift it back to tau.
    ! Only the part interpolated from w_pot should be shifted. The center of
    ! dvscf_long_range and dvsc_bare is already at tau.
    !
    CALL dvscf_shift_center(dvscf_cart, xq, shift_half, +1)
    !
    !--------------------------------------------------------------------------
    ! For IR active materials with nonzero Born effective charge, add the
    ! long-range (dipole) potential after Fourier interpolation
    IF (do_long_range) THEN
      CALL start_clock('dvscf_long')
      !
      ALLOCATE(dvscf_long(dfftp%nnr, nspin_mag, nmodes))
      CALL dvscf_long_range(xq, zeu_r2q, epsilon_r2q, dvscf_long)
      dvscf_cart = dvscf_cart + dvscf_long
      DEALLOCATE(dvscf_long)
      !
      CALL stop_clock('dvscf_long')
    ENDIF ! do_long_range
    !
    !--------------------------------------------------------------------------
    ! In other parts of ph.x code, dvscf is assumed to contain only the
    ! induced part. induced part only.
    ! So, we subtract the bare part from dvscf.
    !
    CALL start_clock('dvscf_bare')
    !
    ALLOCATE(dvscf_bare(dfftp%nnr, nspin_mag, nmodes))
    CALL dvscf_bare_calc(xq, dvscf_bare, .FALSE.)
    dvscf_cart = dvscf_cart - dvscf_bare
    DEALLOCATE(dvscf_bare)
    !
    CALL stop_clock('dvscf_bare')
    !
    !--------------------------------------------------------------------------
    !
    ! Rotate dvscf_cart in Cartesian basis to dvscf in u_in basis
    ! dvscf(:, jmode) = u_in(imode, jmode) * dvscf_cart(:, imode)
    ! we assume nspin_mag == 1
    !
    CALL start_clock('dvscf_cart2u')
    dvscf = (0.d0, 0.d0)
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
        DO is = 1, nspin_mag
          CALL ZAXPY(dfftp%nnr, u_in(imode, jmode), &
            dvscf_cart(1, is, imode), 1, dvscf(1, is, jmode), 1)
        ENDDO ! is
      ENDDO ! jmode
    ENDDO ! imode
    CALL stop_clock('dvscf_cart2u')
    !
    DEALLOCATE(dvscf_cart)
    DEALLOCATE(w_pot)
    DEALLOCATE(w_pot_gathered)
    !
    CALL stop_clock('dvscf_r2q')
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_r2q
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_interpol_close()
    !--------------------------------------------------------------------------
    !! Close units and deallocate arrays. Called at the end of each q point
    !--------------------------------------------------------------------------
    !
    USE mp_pools,    ONLY : me_pool, root_pool
    !
    IMPLICIT NONE
    !
    INTEGER :: irlocal
    !! Real space unit cell index
    !
    DEALLOCATE(rlatt)
    !
    IF (me_pool == root_pool) THEN
      DO irlocal = 1, nrlocal
        CLOSE(UNIT=iunwpot(irlocal), STATUS='KEEP')
      ENDDO
    ENDIF
    !
    DEALLOCATE(iunwpot)
    !
    IF (do_long_range) DEALLOCATE(zeu_r2q)
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_interpol_close
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_shift_center(dvscf_q, xq, shift_half, sign)
  !----------------------------------------------------------------------------
  !! Shift center of phonon potential.
  !! If sign = +1, shift center from origin to tau (used in dvscf_r2q).
  !! If sign = -1, shift center from tau to origin (used in dvscf_q2r).
  !!
  !! For ipol = 1, 2, 3, if shift_half(ipol) is true, the origin is set at
  !! r = 0.5 for direction ipol. If false, the origin is set at r = 0.0.
  !! Setting shift_half = .true. is useful when the size of supercell is odd,
  !! becuase the center of the supercell is at (0.5, 0.5, 0.5).
  !!
  !----------------------------------------------------------------------------
    USE kinds, ONLY : DP
    USE constants, ONLY : tpi
    USE fft_interfaces, ONLY : invfft, fwfft
    USE fft_base, ONLY : dfftp
    USE cell_base, ONLY : bg, at
    USE ions_base, ONLY : nat, tau
    USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g, ngm
    USE noncollin_module, ONLY : nspin_mag
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: dvscf_q(dfftp%nnr, nspin_mag, 3*nat)
    REAL(DP), INTENT(IN) :: xq(3)
    LOGICAL, INTENT(IN) :: shift_half(3)
    INTEGER, INTENT(IN) :: sign
    !
    INTEGER :: imode, na, ig, ipol, is
    REAL(DP) :: arg
    COMPLEX(DP) :: gtau
    COMPLEX(DP), ALLOCATABLE :: aux1(:)
    COMPLEX(DP), ALLOCATABLE :: aux2(:)
    COMPLEX(DP), ALLOCATABLE :: eigqts(:)
    !
    CALL start_clock('dvscf_shift')
    !
    ALLOCATE(aux1(dfftp%nnr))
    ALLOCATE(aux2(dfftp%nnr))
    ALLOCATE(eigqts(nat))
    !
    DO na = 1, nat
      arg = SUM(xq(:) * tau(:, na)) * tpi
      eigqts(na) = CMPLX(COS(arg), -SIN(arg), KIND=DP)
    ENDDO
    !
    ! If shift_half(ipol) == .true., shift the center by at(:, ipol) * 0.5
    !
    DO ipol = 1, 3
      IF (shift_half(ipol)) THEN
        arg = SUM(xq(:) * at(:, ipol)) * tpi * 0.5d0
        eigqts(:) = eigqts(:) * CMPLX(COS(arg), SIN(arg), KIND=DP)
      ENDIF
    ENDDO
    !
    DO imode = 1, 3 * nat
      DO is = 1, nspin_mag
        !
        aux2 = (0.d0, 0.d0)
        aux1 = dvscf_q(:, is, imode)
        CALL fwfft('Rho', aux1, dfftp)
        !
        na = ((imode - 1) / 3) + 1
        !
        DO ig = 1, ngm
          gtau = eigts1(mill(1,ig), na) * eigts2(mill(2,ig), na) * &
                 eigts3(mill(3,ig), na) * eigqts(na)
          !
          ! shift by at(:, ipol) / 2 is multiplication of -1 when mill(ig) is odd
          !
          DO ipol = 1, 3
            IF (shift_half(ipol) .AND. MOD(mill(ipol, ig), 2) /= 0) THEN
              gtau = -gtau
            ENDIF
          ENDDO
          !
          IF (sign == -1) gtau = CONJG(gtau)
          !
          aux2(dfftp%nl(ig)) = aux1(dfftp%nl(ig)) * gtau
        ENDDO
        !
        CALL invfft('Rho', aux2, dfftp)
        dvscf_q(:, is, imode) = aux2
        !
      ENDDO ! nspin_mag
    ENDDO ! imode
    !
    DEALLOCATE(eigqts)
    DEALLOCATE(aux1)
    DEALLOCATE(aux2)
    !
    CALL stop_clock('dvscf_shift')
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_shift_center
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_long_range(xq, zeu, epsilon, dvscf_long)
  !----------------------------------------------------------------------------
  !! This subroutine calculates the long-range dipole potential for given
  !! xq and zeu.
  !! Input xq: q vector in Cartesian coordinate
  !!
  !! Currently, only the dipole part (Frohlich) is implemented. The quadrupole
  !! potential is not implemented.
  !!
  !! [1] Xavier Gonze et al, Comput. Phys. Commun., 107042 (2019)
  !!
  !! Taken from Eq.(13) of Ref. [1]
  !! dvlong(G,q)_{a,x} = 1j * 4pi / Omega * e^2
  !!                   * [ (q+G)_y * Zstar_{a,yx} * exp(-i*(q+G)*tau_a)) ]
  !!                   / [ (q+G)_y * epsilon_yz * (q+G)_z ]
  !!  a: atom index, x, y: Cartesian direction index
  !!
  !! Since QE uses Rydberg units, we multiply the e^2 = 2.0 factor.
  !! The units of q and G are 2pi/a, so we need to divide dvlong by tpiba.
  !!
  !----------------------------------------------------------------------------
    !
    USE kinds,          ONLY : DP
    USE constants,      ONLY : tpi, fpi, e2
    USE fft_base,       ONLY : dfftp
    USE fft_interfaces, ONLY : invfft
    USE gvect,          ONLY : g, ngm
    USE ions_base,      ONLY : nat, tau
    USE cell_base,      ONLY : omega, tpiba
    USE noncollin_module, ONLY : nspin_mag
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xq(3)
    !! Input: q point (in cartesian coordinate)
    REAL(DP), INTENT(IN) :: zeu(3, 3, nat)
    !! Input: Born effective charge
    REAL(DP), INTENT(IN) :: epsilon(3, 3)
    !! Input: Dielectric matrix
    COMPLEX(DP) :: dvscf_long(dfftp%nnr, nspin_mag, 3*nat)
    !! Output: long-range part of the potential, for all modes in Cartesian basis
    !
    INTEGER :: iatm, idir, imode, jdir, ig
    REAL(DP) :: xq_g(3), epsilon_denom, arg
    COMPLEX(DP) :: phase
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !! (dfftp%nnr) long-range part in G space
    !
    ALLOCATE(aux(dfftp%nnr))
    !
    dvscf_long = (0.d0, 0.d0)
    !
    DO imode = 1, 3 * nat
      iatm = (imode - 1) / 3 + 1
      idir = imode - 3 * (iatm - 1)
      !
      ! Compute aux(G) = 1j * 4pi / Omega
      !                * [ (q+G)_y * Zstar_{a,yx} * exp(-i*(q+G)*tau_a)) ]
      !                / [ (q+G)_y * epsilon_yz * (q+G)_z ]
      !
      aux(:) = (0.d0, 0.d0)
      !
      DO ig = 1, ngm
        xq_g(:) = xq(:) + g(:, ig)
        !
        ! skip if xq_g == 0
        IF (SUM(ABS(xq_g)) < 1.d-5) CYCLE
        !
        ! epsilon_denom = sum_yz (q+G)_y * epsilon_yz * (q+G)_z
        epsilon_denom = 0.d0
        DO jdir = 1, 3
          epsilon_denom = epsilon_denom + xq_g(jdir) * SUM(epsilon(jdir,:) * xq_g(:))
        ENDDO
        !
        ! phase = exp(-i*(q+G)*tau_a))
        arg = tpi * SUM(xq_g(:) * tau(:,iatm))
        phase = CMPLX(COS(arg), -SIN(arg), KIND=DP)
        !
        aux(dfftp%nl(ig)) = SUM(zeu(:, idir, iatm) * xq_g(:)) * phase / epsilon_denom
        !
      ENDDO
      !
      aux(:) = aux(:) * (0.d0, 1.d0) * fpi / omega * e2 / tpiba
      !
      ! Fourier transform aux to dvscf_long(:, imode)
      CALL invfft('Rho', aux, dfftp)
      dvscf_long(:, 1, imode) = aux(:)
      !
    ENDDO ! imode
    !
    DEALLOCATE(aux)
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_long_range
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dvscf_bare_calc(xq, dvscf_bare, addnlcc)
  !----------------------------------------------------------------------------
  !!
  !! This subroutine calculates the bare part of the perturbed potential
  !! in the cartesian basis.
  !! Input xq: q vector in Cartesian coordinate
  !!
  !! addnlcc is not implemented. nlcc need not be added here because it is
  !! included in the induced part of dvscf in solve_linter.
  !! (see solve_linter.f90, call to subroutine addcore).
  !!
  !! 2D Coulomb cutoff is not implemented.
  !!
  !! Adapted from PHonon/PH/dvqpsi_us.f90 by Jae-Mo Lihm
  !!
  !----------------------------------------------------------------------------
    !
    USE kinds,          ONLY : DP
    USE constants,      ONLY : tpi
    USE eqv,            ONLY : vlocq
    USE qpoint,         ONLY : eigqts
    USE fft_base,       ONLY : dfftp
    USE fft_interfaces, ONLY : invfft, fft_interpolate
    USE gvect,          ONLY : g, eigts1, eigts2, eigts3, mill, ngm
    USE atom,           ONLY : msh, rgrid
    USE m_gth,          ONLY : setlocq_gth
    USE cell_base,      ONLY : tpiba, tpiba2, omega
    USE ions_base,      ONLY : ntyp => nsp, tau, nat, ityp
    USE uspp_param,     ONLY : upf
    USE noncollin_module, ONLY : nspin_mag
    USE Coul_cut_2D,    ONLY : do_cutoff_2D
    ! USE Coul_cut_2D_ph, ONLY : cutoff_localq
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xq(3)
    ! Input: q point (in cartesian coordinate)
    COMPLEX(DP) :: dvscf_bare(dfftp%nnr, nspin_mag, 3*nat)
    ! Output: bare perturbation of the potential, for all modes
    LOGICAL, INTENT(IN) :: addnlcc
    ! Input: If true, add nlcc contribution
    !
    INTEGER :: nt, na, imode, ig
    COMPLEX(DP) :: gtau, gu, fact, u1, u2, u3, gu0
    COMPLEX(DP), ALLOCATABLE :: uact(:)
    COMPLEX(DP), ALLOCATABLE :: aux1(:)
    !
    IF (addnlcc) CALL errore('dvscf_bare_calc', 'addnlcc not implemented', 1)
    !
    IF (do_cutoff_2D) CALL errore('dvscf_bare_calc', &
        'do_cutoff_2D not implemented', 1)
    IF (nspin_mag /= 1) CALL errore('dvscf_bare_calc', &
        'magnetism not implemented', 1)
    !
    ! ! for 2d calculations, we need to initialize the fact for the q+G
    ! ! component of the cutoff of the COulomb interaction
    ! IF (do_cutoff_2D) call cutoff_fact_qg()
    ! !  in 2D calculations the long range part of vlocq(g) (erf/r part)
    ! ! was not re-added in g-space because everything is caclulated in
    ! ! radial coordinates, which is not compatible with 2D cutoff.
    ! ! It will be re-added each time vlocq(g) is used in the code.
    ! ! Here, this cutoff long-range part of vlocq(g) is computed only once
    ! ! by the routine below and stored
    ! IF (do_cutoff_2D) call cutoff_lr_Vlocq()
    !
    ! Below are adapted from PH/dvqpsi_us.f90
    ! Here, we consider all Cartesian perterbation
    ! uact(jmode) = delta_{imode, jmode}
    !
    ALLOCATE(uact(3*nat))
    ALLOCATE(aux1(dfftp%nnr))
    !
    DO imode = 1, 3 * nat
      uact = 0.d0
      uact(imode) = 1.d0
      !
      !    We start by computing the contribution of the local potential.
      !    The computation of the derivative of the local potential is done in
      !    reciprocal space, then transformed to real space.
      !
      aux1(:) = (0.d0, 0.d0)
      na = (imode - 1) / 3 + 1
      fact = tpiba * (0.d0, -1.d0) * eigqts (na)
      nt = ityp (na)
      u1 = uact (3*(na-1) + 1)
      u2 = uact (3*(na-1) + 2)
      u3 = uact (3*(na-1) + 3)
      gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
      DO ig = 1, ngm
        gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
               eigts3 (mill(3,ig), na)
        gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
        aux1 (dfftp%nl (ig) ) = aux1 (dfftp%nl (ig) ) + vlocq (ig, nt) * gu &
                                * fact * gtau
      ENDDO
      ! jml: do_cutoff_2D not implemented
      ! IF (do_cutoff_2D) then
      !    call cutoff_localq( aux1, fact, u1, u2, u3, gu0, nt, na)
      ! ENDIF
      !
      ! Now we transform dV_loc/dtau from G space to real space
      !
      CALL invfft ('Rho', aux1, dfftp)
      !
      dvscf_bare(:, 1, imode) = aux1(:)
      !
    ENDDO ! imode
    !
    DEALLOCATE(uact)
    DEALLOCATE(aux1)
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE dvscf_bare_calc
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE multiply_iqr(dfft, xq, func)
  !----------------------------------------------------------------------------
  !!
  !! Multiply exp(i*q*r) to func.
  !! Real space indexing is adapted from Modules/compute_dipole.f90
  !!
  !----------------------------------------------------------------------------
    USE kinds,               ONLY : DP
    USE constants,           ONLY : tpi
    USE cell_base,           ONLY : at
    USE fft_types,           ONLY : fft_type_descriptor, fft_index_to_3d
    !
    IMPLICIT NONE
    !
    TYPE (fft_type_descriptor), INTENT(IN) :: dfft
    ! fft_type_descriptor. dffts or dfftp
    REAL(DP), INTENT(IN) :: xq(3)
    ! q vector in cartesian coordinate
    COMPLEX(DP) :: func(dfft%nnr)
    ! input, output: real-space function
    !
    LOGICAL :: offrange
    INTEGER :: ir, i, j, k, ir_end
    REAL(DP) :: arg, xq_cry(3)
    COMPLEX(DP) :: phase
    !
    xq_cry = xq
    CALL cryst_to_cart(1, xq_cry, at, -1)
    !
#if defined (__MPI)
    ir_end = MIN(dfft%nnr, dfft%nr1x*dfft%my_nr2p*dfft%my_nr3p)
#else
    ir_end = dfft%nnr
#endif
    !
    DO ir = 1, ir_end
      !
      CALL fft_index_to_3d(ir, dfft, i, j, k, offrange)
      IF ( offrange ) CYCLE
      !
      ! (i,j,k) is the zero-based coordinate of the real-space grid
      arg = tpi * (  xq_cry(1) * REAL(i, DP) / REAL(dfft%nr1, DP) &
                   + xq_cry(2) * REAL(j, DP) / REAL(dfft%nr2, DP) &
                   + xq_cry(3) * REAL(k, DP) / REAL(dfft%nr3, DP)  )
      phase = CMPLX( COS(arg), SIN(arg), kind=DP )
      !
      func(ir) = func(ir) * phase
      !
    END DO ! ir
  !----------------------------------------------------------------------------
  END SUBROUTINE multiply_iqr
  !----------------------------------------------------------------------------
  !
!------------------------------------------------------------------------------
END MODULE dvscf_interpolate
!------------------------------------------------------------------------------
