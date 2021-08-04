  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE cum_mod
  !----------------------------------------------------------------------
  !!
  !! This module contains the various routines use for cumulant expansion
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_cumulant()
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the electron-
    !!  phonon interaction using the (retarded) cumulant expansion method.
    !!  Use the frequency-dependent self-energy from the one-shot Migdal approximation
    !!  read from specfun_sup.elself.
    !!
    !!  See e.g. PRL 77, 2268 (1996) and PRB 90, 085112 (2014):
    !!  the cumulant spectral function can be calculated using a series of convolutions
    !!  in frequency space, or using FFT.
    !!  To converge the FFT with a small dt, we use zero padding of ImSigma outside the
    !!  calculated frequency range. If the convergence is not satisfactory, one can
    !!  set 'fact' to a larger value (6 is used as default, see the SUBROUTINE cumulant_time).
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP, i4b
    USE constants_epw, ONLY : kelvin2eV, two, zero, ryd2ev, ryd2mev, ci
    USE constants,     ONLY : pi
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iospectral_sup, iospectral_cum
    USE epwcom,        ONLY : wmin_specfun, wmax_specfun, nw_specfun, &
                              bnd_cum, nstemp, eliashberg
    USE elph2,         ONLY : ibndmin, ibndmax, gtemp
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    CHARACTER(LEN = 256) :: filespecsup
    CHARACTER(LEN = 64) :: line
    !! Auxiliary string
    CHARACTER(LEN = 64) :: filespec
    !! Spectral function file
    INTEGER(KIND = i4b) :: maxrecs = 1.0E9
    !! Maximum number of lines in the specfun_sup.elself
    INTEGER :: ios
    !! Status of opening file
    INTEGER :: im
    !! Dummy counter when reading
    INTEGER :: i1, i2
    !! Auxiliary indices
    INTEGER :: ik
    !! K-point counter
    INTEGER :: iw
    !! Frequency counter
    INTEGER :: ibnd
    !! Band index
    INTEGER :: nk
    !! Total number of k-points
    INTEGER :: i0
    !! Energy index of Fermi level (w=0)
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: dw
    !! Freq. increment
    REAL(KIND = DP) :: e_thresh
    !! Do the cumulant only for states with energy below e_thresh
    REAL(KIND = DP) :: a1, a2, a3, a4
    !! Auxiliary variables
    REAL(KIND = DP) :: ekk
    !! K-S energy
    REAL(KIND = DP) :: zeta
    !! Z factor
    REAL(KIND = DP), ALLOCATABLE :: ww(:)
    !! Frequency variable
    REAL(KIND = DP), ALLOCATABLE :: ek(:)
    !! Eigenvalues for the band bnd_cum
    REAL(KIND = DP), ALLOCATABLE :: sigmar(:, :) ! kpt, omega
    !! Real self-energy for bnd_cum
    REAL(KIND = DP), ALLOCATABLE :: sigmai(:, :)
    !! Imaginary self-energy for bnd_cum
    REAL(KIND = DP), ALLOCATABLE :: a_mig(:, :)
    !! Migdal spectral function (same as in spectral_func.f90 )
    REAL(KIND = DP), ALLOCATABLE :: a_cw(:)
    !! Cumulant spectral function (convolutions)
    REAL(KIND = DP), ALLOCATABLE :: a_ct(:)
    !! Cumulant spectral function (FFT)
    REAL(KIND = DP), ALLOCATABLE :: a_tmp(:)
    !! Temporary spectral function
    !
    ! e_thresh can be changed if one needs e.g. only states below the Fermi level
    e_thresh = 10.d0 / ryd2ev ! referred to the Fermi level
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    !
    WRITE(stdout, '(/5x,a)') REPEAT('=',75)
    WRITE(stdout, '(5x,a)') 'Performing the CUMULANT EXPANSION of the retarded Greens function to obtain'
    WRITE(stdout, '(5x,a)') 'the electron spectral function.'
    WRITE(stdout, '(5x,a)') 'There is no spin degeneracy factor, as in spectral_func.f90'
    WRITE(stdout, '(5x,a)') 'Warning: the routine is sequential but very fast.'
    WRITE(stdout, '(5x,a/)') REPEAT('=',75)
    !
    DO itemp = 1, nstemp
      WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
      filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K'
      OPEN (UNIT = iospectral_sup, FILE = filespecsup, STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('spectral_cumulant', 'opening file specfun_sup.elself', ABS(ios))
      !
      ! determine number of k points, ibndmin, ibndmax
      DO im = 1, 6
        READ(iospectral_sup, '(a)') line
      ENDDO
      DO im = 1, maxrecs
        READ (iospectral_sup, *, IOSTAT = ios) i1, i2
        IF (im == 1) ibndmin = i2
        IF (ios /= 0) EXIT
        IF (im == maxrecs) CALL errore('spectral_cumulant', 'increase maxrecs', 1)
      ENDDO
      !
      REWIND(iospectral_sup)
      !
      nk = i1
      ibndmax = i2
      WRITE(stdout, '(5x,a/)') "Read self-energy from file specfun_sup.elself"
      WRITE(stdout, '(5x,a,i4,a,i4,a,i4,a,f12.6/)') "Check: nk = ", nk, &
             ", ibndmin = ", ibndmin, ", ibndmax = ", ibndmax, " kbT (eV) = ", gtemp(itemp) * ryd2ev
      !
      ALLOCATE(ww(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating ww', 1)
      ALLOCATE(ek(nk), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating ek', 1)
      ALLOCATE(sigmar(nk, nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating sigmar', 1)
      ALLOCATE(sigmai(nk, nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating sigmai', 1)
      ALLOCATE(a_mig(nw_specfun, nk), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_mig', 1)
      ALLOCATE(a_cw(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_cw', 1)
      ALLOCATE(a_ct(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_ct', 1)
      ALLOCATE(a_tmp(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_tmp', 1)
      !
      ! read and store Kohn-Sham energy, energy grid, real and im sigma for designated band
      DO im = 1,6
        READ(iospectral_sup, '(a)') line
      ENDDO
      DO ibnd = 1, ibndmax - ibndmin + 1
        DO ik = 1, nk
          DO iw = 1, nw_specfun
            READ(iospectral_sup,*) i1, i2, a1, a2, a3, a4
            IF (i2 == bnd_cum) THEN
              ! ek, w read in eV; Sigma read in meV
              ek(ik) = a1 / ryd2ev
              ww(iw) = a2 / ryd2ev
              sigmar(ik, iw) = a3 / ryd2mev ! / ( EXP(ww(iw)/eptemp )+1.d0 )
              sigmai(ik, iw) = a4 / ryd2mev ! / ( EXP(ww(iw)/eptemp )+1.d0 )
              ! spec func as in spectral_func.f90
              a_mig(iw, ik) = ABS(sigmai(ik, iw)) / pi / ((ww(iw) - ek(ik) - sigmar(ik, iw))**two + (sigmai(ik, iw) )**two)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE(iospectral_sup)
      !
      ! open file for cumulant spectral function
      WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
      IF (bnd_cum < 10) THEN
        WRITE(filespec, '(a,i1,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
      ELSEIF (bnd_cum > 9 .AND. bnd_cum < 100) THEN
        WRITE(filespec, '(a,i2,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
      ELSE
        WRITE(filespec, '(a,i3,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
      ENDIF
      OPEN(UNIT = iospectral_cum, FILE = filespec)
      !
      WRITE(iospectral_cum, '(a)') '#       k   Energy [eV]            A(k,w) [meV^-1]            Z-factor    '
      WRITE(iospectral_cum, '(a)') '#                       with convolutions  |    using FFT    '
      WRITE(stdout, '(8x,a)') 'k   Energy [eV]            A(k,w) [meV^-1]            Z-factor    '
      WRITE(stdout, '(8x,a)') '                with convolutions  |    using FFT '
      !
      ! define index corresponding to omega=0 (Fermi level)
      i0 = MINLOC(ABS(ww(:)), DIM = 1)
      IF (ABS(ww(i0)) > dw) CALL errore('spectral_cumulant', 'w=0 needs to be included in [wmin:wmax]', 1)
      a_cw = zero
      a_ct = zero
      a_tmp = zero
      !
      DO ik = 1, nk
        IF (ek(ik) < e_thresh) THEN
          !
          ekk = ek(ik)
          !
          ! Cumulant from convolutions in frequency space (ImSigma>0 in EPW)
          CALL cumulant_conv(ekk, ww, sigmar(ik, :), -sigmai(ik, :), a_mig(:, ik), a_cw, zeta)
          !
          ! Cumulant calculated in time domain + FFT (ImSigma>0 in EPW)
          CALL cumulant_time(ekk, ww, sigmar(ik, :), -sigmai(ik, :), a_mig(:, ik), a_tmp)
          !
          DO iw = 1, nw_specfun
            !
            ! map the indices of the FFT frequency grid onto the original one
            IF (iw >= i0) THEN
              a_ct(iw) = a_tmp(iw - i0 + 1)
            ELSE
              a_ct(iw) = a_tmp(iw + nw_specfun - i0 + 1)
            ENDIF
            !
            ! write cumulant spectral function on file (in meV^-1, as in spectral_func.f90)
            ! 3rd column: A_cum using convolutions; 4th column: A_cum using FFT
            IF (iw == 1) THEN
              WRITE(iospectral_cum, '(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,f8.4)') &
                    ik, ww(iw) * ryd2ev, a_cw(iw) / ryd2mev, a_ct(iw) / ryd2mev, zeta
              WRITE(stdout,'(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,f8.4)') &
                    ik, ww(iw) * ryd2ev, a_cw(iw) / ryd2mev, a_ct(iw) / ryd2mev, zeta
            ELSE
              WRITE(iospectral_cum, '(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7)') &
                    ik, ww(iw) * ryd2ev, a_cw(iw) / ryd2mev, a_ct(iw) / ryd2mev !/ ( EXP(ww(iw)/eptemp )+1.d0 )
              WRITE(stdout, '(2x,i7,2x,f10.5,3x,e16.7,3x,e16.7)') &
                    ik, ww(iw) * ryd2ev, a_cw(iw) / ryd2mev, a_ct(iw) / ryd2mev !/ ( EXP(ww(iw)/eptemp )+1.d0 )
              ! uncomment to multiply by Fermi occupation factor
            ENDIF
            !
          ENDDO
          !
          WRITE(iospectral_cum, '(a)') ' '
          WRITE(stdout, '(a)') ' '
          !
        ENDIF ! only states below energy e_thresh
        !
      ENDDO ! main loop k
      !
      WRITE(stdout, '(5x,a)') 'The file specfun_cum[BND].elself has been correctly written'
      !
      CLOSE(iospectral_cum)
      !
      DEALLOCATE(ww, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating ww', 1)
      DEALLOCATE(ek, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating ek', 1)
      DEALLOCATE(sigmar, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating sigmar', 1)
      DEALLOCATE(sigmai, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating sigmai', 1)
      DEALLOCATE(a_mig, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_mig', 1)
      DEALLOCATE(a_cw, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_cw', 1)
      DEALLOCATE(a_ct, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_ct', 1)
      DEALLOCATE(a_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_tmp', 1)
    ENDDO !itemp
    !
    ! Deallocate temperature when no  supercond
    IF (.NOT. eliashberg) THEN
      DEALLOCATE(gtemp, STAT = ierr)
      IF (ierr /= 0) CALL errore('cum_mod', 'Error deallocating gtemp', 1)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_cumulant
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE cumulant_conv(ek, ww, sigmar, sigmai, a_mig, a_cum, zeta)
    !-----------------------------------------------------------------------
    !!
    !! Convolution for the cumulant expansion.
    !!
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : two, zero, ci, ryd2ev
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : degaussw, wmin_specfun, wmax_specfun, nw_specfun
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: ek
    !! K-S energy
    REAL(KIND = DP), INTENT(in) :: ww(nw_specfun)
    !! Frequency variable
    REAL(KIND = DP), INTENT(in) :: sigmar(nw_specfun)
    !! Real self-energy
    REAL(KIND = DP), INTENT(in) :: sigmai(nw_specfun)
    !! Imaginary self-energy
    REAL(KIND = DP), INTENT(in) :: a_mig(nw_specfun)
    !! Migdal spectral function
    REAL(KIND = DP), INTENT(out) :: a_cum(nw_specfun)
    !! Cumulant spectral function
    REAL(KIND = DP), INTENT(out) :: zeta
    !! Z factor
    !
    ! local variables
    INTEGER :: iw
    !! Frequency counter
    INTEGER :: iw2
    !! Second freq. counter for the convolution
    INTEGER :: indw
    !! Auxillary indices for convolution
    INTEGER :: ind1
    !! Auxillary indices for convolution
    INTEGER :: isat
    !! Index to number the satellite
    INTEGER :: iqp
    !! Energy index of renormalized quasiparticle
    INTEGER :: iks
    !! Energy index for the KS quasiparticle
    INTEGER :: i0
    !! Energy index of Fermi level
    REAL(KIND = DP) :: dw
    !! Freq. increment
    REAL(KIND = DP) :: eqp
    !! Renpormalized quasiparticle energy
    REAL(KIND = DP) :: si_ks
    !! ImSigma at KS energy
    REAL(KIND = DP) :: si_qp
    !! ImSigma at renormalized energy
    REAL(KIND = DP) :: diS
    !! Derivative of ImSigma at KS energy
    REAL(KIND = DP) :: drS
    !! Derivative of ReSigma at qp energy
    REAL(KIND = DP) :: a_qp(nw_specfun)
    !! Quasiparticle contribution to the spectral function
    REAL(KIND = DP) :: a_s(nw_specfun)
    !! Temporary quantity needed to compute the satelite
    REAL(KIND = DP) :: conv(nw_specfun)
    !! Temporary quantity to compute the convolution
    REAL(KIND = DP) :: a_s1(nw_specfun), a_s2(nw_specfun), a_s3(nw_specfun)
    !! satellites
    ! 
    ! Initialization
    diS = zero
    drS = zero 
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    i0 = MINLOC(ABS(ww(:)), DIM = 1)
    ! define index and energy of renormalized qp (note: qp energy needs to be the absolute max of a_mig)
    iqp = MAXLOC(a_mig(:), DIM = 1)
    eqp = ww(iqp)
    !
    ! define index corresponding to unrenormalized (Kohn-Sham) energy
    iks = MINLOC(ABS(ww(:) - ek), DIM = 1)
    !
    si_qp = ABS(sigmai(iqp))
    si_ks = ABS(sigmai(iks))
    ! finite difference derivatives of Im and Re Sigma
    IF (iks - 1 > 0) diS = (ABS(sigmai(iks + 1)) - ABS(sigmai(iks - 1))) / (2.d0 * dw)
    IF (iqp - 1 > 0) drS = (sigmar(iqp + 1) - sigmar(iqp - 1)) / (2.d0 * dw)
    zeta = EXP(drS)
    !
    ! calculate Aqp and As1
    DO iw = 1, nw_specfun
      !
      a_qp(iw) = zeta * si_qp / pi / ((ww(iw) - eqp)**two + (si_qp)**two)
      conv(iw) = a_qp(iw)
      !
      ind1 = iks + iw - i0
      IF (ind1 > 0 .AND. ind1 < nw_specfun) THEN !RC
        a_s(iw) = (ABS(sigmai(ind1)) - si_ks - (ww(iw)) * diS) * REAL(1.d0 / (ww(iw) - ci * degaussw)**2.d0) / pi
      ELSE
        a_s(iw) = (-si_ks - (ww(iw)) * diS) * REAL(1.d0 / (ww(iw) - ci * degaussw)**2.d0) / pi
      ENDIF
      !
    ENDDO
    !
    DO isat = 1, 3
      !
      a_cum = zero
      !
      DO iw = 1, nw_specfun
        !
        DO iw2 = 1, nw_specfun
          !
          indw = i0 + iw - iw2
          IF (indw <= nw_specfun .AND. indw > 0) THEN
             a_cum(iw) = a_cum(iw) + ABS(a_s(iw2) * conv(indw)) * dw
          ENDIF
          !
        ENDDO
        !
      ENDDO
      !
      IF (isat == 1) a_s1 = a_cum
      IF (isat == 2) a_s2 = a_cum / 2.d0
      IF (isat == 3) a_s3 = a_cum / 6.d0
      conv = a_cum
      !
    ENDDO ! isat
    !
    DO iw = 1, nw_specfun
      !
      a_cum(iw) = a_qp(iw) + a_s1(iw) + a_s2(iw) + a_s3(iw)
      !
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cumulant_conv
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE cumulant_time(ek, ww, sigmar, sigmai, a_mig, a_cum)
    !-----------------------------------------------------------------------
    !!
    !! Cumulant calculated in the time domain
    !!
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : two, zero, czero, ci, ryd2ev
    USE constants,     ONLY : pi
    USE fft_scalar,    ONLY : cfft3d
    USE epwcom,        ONLY : degaussw, wmin_specfun, wmax_specfun, nw_specfun
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: ek
    !! K-S energy
    REAL(KIND = DP), INTENT(in) :: ww(nw_specfun)
    !! Frequency variable
    REAL(KIND = DP), INTENT(in) :: sigmar(nw_specfun)
    !! Real self-energy
    REAL(KIND = DP), INTENT(in) :: sigmai(nw_specfun)
    !! Imaginary self-energy
    REAL(KIND = DP), INTENT(in) :: a_mig(nw_specfun)
    !! Migdal spectral function
    REAL(KIND = DP), INTENT(out) :: a_cum(nw_specfun)
    !! Cumulant spectral function
    !
    ! local variables
    !
    INTEGER :: iw
    !! Frequency counters
    INTEGER :: it
    !! Time counters
    INTEGER :: iqp
    !! Energy index of renormalized quasiparticle
    INTEGER :: iks
    !! Energy index for the KS quasiparticle
    INTEGER :: i0
    !! Energy index of Fermi level
    INTEGER :: ind1, ind2
    !! aux indices
    INTEGER :: nw_new
    !! number of points used for FFT
    !! Zero padding is used for ImSigma outside the energy range computed.
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: fact
    !! factor used to increase the number of points for the FFT;
    !! fact=4 gives good convergence but it can be increased if needed.
    REAL(KIND = DP) :: dw
    !! Freq. increment
    REAL(KIND = DP) :: eqp
    !! Renormalized quasiparticle energy
    REAL(KIND = DP) :: drS, zeta
    !! Derivative of ReSigma at qp energy, Z factor
    REAL(KIND = DP) :: smeart
    !! small broadening
    REAL(KIND = DP) :: dt
    !! time increment
    REAL(KIND = DP) :: tmin, tmax
    !! min and max of time interval
    REAL(KIND = DP) :: tt
    !! time variable
    COMPLEX(KIND = DP), ALLOCATABLE :: cumS(:)
    !! satellite part of the cumulant function
    COMPLEX(KIND = DP), ALLOCATABLE :: cum(:)
    !! complete cumulant function
    COMPLEX(KIND = DP) :: qpfac
    !! quasiparticle factor
    !
    fact = 6.d0
    ! fact=1 corresponds to FFT with the original no. of points nw_specfun
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    smeart = -degaussw ! correct sign for RC
    !
    tmin = 0.d0
    tmax = two * pi / dw
    dt = two * pi / ((wmax_specfun - wmin_specfun) * fact)
    nw_new = INT(fact * (nw_specfun - 1) + 1) ! to be consistent with dt above
    !
    ALLOCATE(cumS(nw_new), STAT = ierr)
    IF (ierr /= 0) CALL errore('cumulant_time', 'Error allocating cumS', 1)
    ALLOCATE(cum(nw_new), STAT = ierr)
    IF (ierr /= 0) CALL errore('cumulant_time', 'Error allocating cum', 1)
    !
    i0 = MINLOC(ABS(ww(:)), DIM = 1)
    !
    ! define index and energy of renormalized qp (note: qp energy needs to be the absolute max of a_mig)
    iqp = MAXLOC(a_mig(:), DIM = 1)
    eqp = ww(iqp)
    !
    ! define index corresponding to unrenormalized (Kohn-Sham) energy
    iks = MINLOC(ABS(ww(:) - ek), DIM = 1)
    !
    drS = (sigmar(iqp + 1) - sigmar(iqp - 1)) / (two * dw)
    qpfac = EXP(-ci * (ek + sigmar(iqp)) + smeart * 0.5d0)
    zeta = EXP(drS)
    !
    cumS = czero
    DO iw = 1, nw_new
      !
      ! the w shift is needed because FFT uses positive w \in [0:Omega], Omega=wmax-wmin
      IF (iw <= (nw_specfun - i0 + 1)) THEN
        ind1 = iw + i0 - 1
        cumS(iw) = dw * ABS(sigmai(ind1)) / pi * REAL(1.d0 / (ek - ww(ind1) - ci * smeart)**two)
      ELSE IF (iw > (nw_new - i0 + 1)) THEN
        ind2 = iw - nw_new + i0 - 1
        cumS(iw) = dw * ABS(sigmai(ind2)) / pi * REAL(1.d0 / (ek - ww(ind2) - ci * smeart)**two)
      ENDIF
      !
    ENDDO
    !
    CALL cfft3d(cumS(:), nw_new, 1, 1, nw_new, 1, 1, 1, -1)
    !this is needed because cfft3d(...,-1) carries a renomalization factor 1/nw_new
    cumS = cumS * nw_new
    !
    DO it = 1, nw_new
      !
      tt = tmin + DBLE(it - 1) * dt
      cumS(it) = cumS(it) * EXP(ci * ek * tt)
      cum(it) = zeta * (qpfac**tt) * EXP(cumS(it))
      !
    ENDDO
    !
    CALL cfft3d(cum(:), nw_new, 1, 1, nw_new, 1, 1, 1, 1)
    cum = cum * dt / pi
    !
    ! extract the spectral function a_cum on the original w FFT grid (nw_specfun points)
    DO iw = 1, nw_specfun
      IF (iw <= (nw_specfun - i0 + 1)) THEN
        a_cum(iw) = REAL(REAL(cum(iw)))
      ELSE
        ind1 = iw + nw_new - nw_specfun
        a_cum(iw) = REAL(REAL(cum(ind1)))
      ENDIF
    ENDDO
    !
    DEALLOCATE(cumS, STAT = ierr)
    IF (ierr /= 0) CALL errore('cumulant_time', 'Error deallocating cumS', 1)
    DEALLOCATE(cum, STAT = ierr)
    IF (ierr /= 0) CALL errore('cumulant_time', 'Error deallocating cum', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cumulant_time
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE cum_mod
  !-----------------------------------------------------------------------
