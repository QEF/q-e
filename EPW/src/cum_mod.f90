  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
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
  !! In this module, the hole / electron cumulant is calculated using the
  !! convolution definition of the cumulant. This definition allows to 
  !! print out the quasi-particle and satellite contributions separately.
  !! 
  !! The implementation of the Fourier-transformed time-ordered cumulant
  !! has been retired.
  !!
  !! This module was initially implemented by C. Verdi and F. Caruso and
  !! later updated by S. Ponce'. The latest version from July 2021 is by
  !! Nikolaus Kandolf.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_cumulant()
    !-----------------------------------------------------------------------
    !!
    !!  Compute the second-order cumulant expansion spectral function using
    !!  the electron self-energy from the file specfun_sup.elself calculated
    !!  earlier. For occupied states e_{nk}<E_F, this routine calculates the
    !!  hole (or lesser) cumulant, for unoccupied states e_{nk}>E_F, the 
    !!  electron (or greater) cumulant. 
    !!
    !!  The implementation follows essentially the derivations given in 
    !!  Ref.[1]. For a simpler implementation, the quasi-particle function 
    !!  is written in terms of the full time-ordered self-energy as in 
    !!  Ref.[2], and not in terms of the lesser and greater self-energy as 
    !!  in Ref.[1]. Satellites are constructed from the imaginary part of 
    !!  the lesser or greater self-energy; see also documentation on 
    !!  https://epw-code.org
    !!
    !!  [1] Gumhalter et al., PRB 94, 035103 (2016)
    !!  [2] Aryasetiawan, in Strong coulomb correlations in electronic 
    !!      structure calculations edited by V.I. Anisimov (2000), 
    !!      chapter 1
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP, i4b
    USE constants_epw, ONLY : kelvin2eV, zero, two, ryd2ev, ryd2mev, ci
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
    REAL(KIND = DP) :: a1, a2, a3, a4
    !! Auxiliary variables
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
    REAL(KIND = DP), ALLOCATABLE :: a_ce(:)
    !! Full cumulant spectral function
    REAL(KIND = DP), ALLOCATABLE :: a_qp(:)
    !! Quasi-particle part of the cumulant spectrum
    REAL(KIND = DP), ALLOCATABLE :: a_s1(:)
    !! First satellite of the cumulant spectrum
    !
    WRITE(stdout, '(/5x,a)') REPEAT('=',75)
    WRITE(stdout, '(5x,a)') 'Calculating the cumulant expansion spectral function from the self-energy in specfun_sup.elself'
    WRITE(stdout, '(5x,a)') 'There is no spin degeneracy factor, as in spectral_func.f90'
    WRITE(stdout, '(5x,a)') 'Warning: the routine is sequential but very fast.'
    WRITE(stdout, '(5x,a/)') REPEAT('=',75)
    !
    !
    !  Read in self-energy from specfun_sup.elself
    !
    DO itemp = 1, nstemp
      !
      WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
      !
      filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K'
      OPEN (UNIT = iospectral_sup, FILE = filespecsup, STATUS = 'old', IOSTAT = ios)
      !
      IF (ios /= 0) CALL errore('spectral_cumulant', 'opening file specfun_sup.elself', ABS(ios))
      !
      !
      !  Determine number of k points, ibndmin, ibndmax
      !
      DO im = 1, 6
        !
        READ(iospectral_sup, '(a)') line
        !
      ENDDO
      !
      DO im = 1, maxrecs
        !
        READ (iospectral_sup, *, IOSTAT = ios) i1, i2
        !
        IF (im == 1) ibndmin = i2
        IF (ios /= 0) EXIT
        IF (im == maxrecs) CALL errore('spectral_cumulant', 'increase maxrecs', 1)
        !
      ENDDO
      !
      REWIND(iospectral_sup)
      !
      nk = i1
      ibndmax = i2
      !
      WRITE(stdout, '(5x,a/)') "Read self-energy from file specfun_sup.elself"
      WRITE(stdout, '(5x,a,i4,a,i4,a,i4,a,f12.6/)') "Check: nk = ", nk, &
             ", ibndmin = ", ibndmin, ", ibndmax = ", ibndmax, " kbT (eV) = ", gtemp(itemp) * ryd2ev
      !
      !
      !  Allocate arrays for frequency, KS energy, self-energy and spectral functions
      !
      ALLOCATE(ww(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating ww', 1)
      ALLOCATE(ek(nk), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating ek', 1)
      ALLOCATE(sigmar(nk, nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating sigmar', 1)
      ALLOCATE(sigmai(nk, nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating sigmai', 1)
      ALLOCATE(a_ce(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_ce', 1)
      ALLOCATE(a_qp(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_qp', 1)
      ALLOCATE(a_s1(nw_specfun), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error allocating a_s1', 1)
      !
      !
      !  Read in KS energy, frequency grid and self-energy for designated band
      !
      DO im = 1,6
        !
        READ(iospectral_sup, '(a)') line
        !
      ENDDO
      !
      !
      !
      !  Read in KS energy and frequency in eV, self-energy in meV
      !
      DO ibnd = 1, ibndmax - ibndmin + 1
        !
        DO ik = 1, nk
          !
          DO iw = 1, nw_specfun
            !
            READ(iospectral_sup,*) i1, i2, a1, a2, a3, a4
            !
            IF (i2 == bnd_cum) THEN
              !
              ek(ik) = DBLE(a1 / ryd2ev)
              ww(iw) = DBLE(a2 / ryd2ev)
              sigmar(ik, iw) = DBLE(a3 / ryd2mev) ! / ( EXP(ww(iw)/eptemp )+1.d0 ) 
              sigmai(ik, iw) = DBLE(a4 / ryd2mev) ! / ( EXP(ww(iw)/eptemp )+1.d0 )
              !
            ENDIF
            !
          ENDDO
          !
        ENDDO
        !
      ENDDO
      !
      CLOSE(iospectral_sup)
      !
      !
      !  Calculated self-energy is retarded: 
      !  Expecting Im Sigma_k(w) < 0 for all momenta and frequencies
      !
      IF( ANY(sigmai < zero ) ) THEN
        !
        WRITE(stdout,'(5x,a/)') 'WARNING: Some values of Im Sigma are larger than zero!'
        WRITE(stdout,'(5x,a/)') '         Check validity of spectral function result!'
        !
      ENDIF
      !
      !  Check if KS energy lies within the spectral function frequency window
      !
      DO ik = 1, nk
        !
        IF (ek(ik) < ww(1) .OR. ek(ik) > ww(nw_specfun-1)) THEN
          !
          WRITE(stdout, '(5x,a,i5/)') 'WARNING: Frequency window too small at ik = ', ik          
          WRITE(stdout, '(5x,a,i5/)') '         Frequency window for spectral function must contain KS energy.'
          WRITE(stdout, '(5x,a/)')    '         Skipping this k point'
          !
        ENDIF
        !
      ENDDO
      !
      ! Open file for cumulant spectral function output
      !
      WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
      !
      IF (bnd_cum < 10) THEN
        !
        WRITE(filespec, '(a,i1,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
        !
      ELSEIF (bnd_cum > 9 .AND. bnd_cum < 100) THEN
        !
        WRITE(filespec, '(a,i2,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
        !
      ELSE
        !
        WRITE(filespec, '(a,i3,a,a,a)') 'specfun_cum', bnd_cum, '.elself.', trim(adjustl(tp)), 'K'
        !
      ENDIF
      !
      OPEN(UNIT = iospectral_cum, FILE = filespec)
      !
      ! Prepare output to specfun_cum[bnd_cum].elself.[T]K and to standard output
      !
      WRITE(iospectral_cum, '(6x,a)') &
              '#k   ww-E_F[eV]     A(k,w)[meV^{-1}]   A^{QP}[meV^{-1}]   A^{S1}[meV^{-1}]  Z-factor'
      WRITE(stdout,'(6x,a)') &
              '#k   ww-E_F[eV]     A(k,w)[meV^{-1}]   A^{QP}[meV^{-1}]   A^{S1}[meV^{-1}]  Z-factor'
      !
      !
      DO ik = 1, nk ! main k loop
        !
        a_ce = zero
        a_qp = zero
        a_s1 = zero
        !
        ! Call cumulant convolution subroutine for k points whose 
        ! KS energy lies within the frequency range of the 
        ! spectral function.
        !
        IF (ek(ik) > ww(1) .AND. ek(ik) < ww(nw_specfun-1)) THEN
          CALL cumulant_conv(ek(ik), ww, sigmar(ik, :), sigmai(ik, :), a_ce, a_qp, a_s1, zeta)
        ENDIF
        !
        ! Write data to specfun_cum[bnd_cum].elself.[T]K and to standard output
        !
        DO iw = 1, nw_specfun
          !
          IF (iw == 1) THEN
            !
            WRITE(iospectral_cum, '(x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,e16.7,3x,f8.4)') &
                  ik, ww(iw) * ryd2ev, a_ce(iw) / ryd2mev, a_qp(iw) / ryd2mev, a_s1(iw) / ryd2mev, zeta
            WRITE(stdout,'(x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,e16.7,3x,f8.4)') &
                  ik, ww(iw) * ryd2ev, a_ce(iw) / ryd2mev, a_qp(iw) / ryd2mev, a_s1(iw) / ryd2mev, zeta
            !
          ELSE
            !
            WRITE(iospectral_cum, '(x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,e16.7)') &
                  ik, ww(iw) * ryd2ev, a_ce(iw) / ryd2mev, a_qp(iw) / ryd2mev, a_s1(iw) / ryd2mev
            WRITE(stdout, '(x,i7,2x,f10.5,3x,e16.7,3x,e16.7,3x,e16.7)') &
                  ik, ww(iw) * ryd2ev, a_ce(iw) / ryd2mev, a_qp(iw) / ryd2mev, a_s1(iw) / ryd2mev
            !
          ENDIF
          !
        ENDDO
        !
        WRITE(iospectral_cum, '(a)') ' '
        WRITE(stdout, '(a)') ' '
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
      DEALLOCATE(a_ce, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_cw', 1)
      DEALLOCATE(a_qp, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_ct', 1)
      DEALLOCATE(a_s1, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_cumulant', 'Error deallocating a_tmp', 1)
      !
    ENDDO ! Loop over temperatures
    !
    ! Deallocate temperature if supercond = .FALSE.
    !
    IF (.NOT. eliashberg) THEN
      !
      DEALLOCATE(gtemp, STAT = ierr)
      IF (ierr /= 0) CALL errore('cum_mod', 'Error deallocating gtemp', 1)
      !
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_cumulant
    !-----------------------------------------------------------------------
    !
    SUBROUTINE cumulant_conv(ek, ww, sigmar, sigmai, a_ce, a_qp, a_s1, zeta)
    !-----------------------------------------------------------------------
    !!
    !! Calculated cumulant spectral function from convolution
    !!
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero, one, two, ci, ryd2ev, ryd2mev
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : degaussw, wmin_specfun, wmax_specfun, nw_specfun
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: ek
    !! K-S energy
    REAL(KIND = DP), INTENT(in) :: ww(nw_specfun)
    !! Frequency variable
    REAL(KIND = DP), INTENT(in) :: sigmar(nw_specfun)
    !! Real self-energy
    REAL(KIND = DP), INTENT(INOUT) :: sigmai(nw_specfun)
    !! Imaginary self-energy
    REAL(KIND = DP), INTENT(out) :: a_ce(nw_specfun)
    !! Cumulant spectral function
    REAL(KIND = DP), INTENT(out) :: a_qp(nw_specfun)
    !! Quasi-particle part of spectral function
    REAL(KIND = DP), INTENT(out) :: a_s1(nw_specfun)
    !! First satellite of spectral function
    REAL(KIND = DP), INTENT(out) :: zeta
    !! Z factor
    !
    ! local variables
    INTEGER :: iw
    !! Frequency counter
    INTEGER :: iw2
    !! Second freq. counter for the convolution
    INTEGER :: icoul
    !! Auxillary indices for convolution
    INTEGER :: ind1
    !! Auxillary indices for convolution
    INTEGER :: iek
    !! Energy index for the KS quasiparticle
    INTEGER :: i0
    !! Energy index of Fermi level
    REAL(KIND = DP) :: dw
    !! Freq. increment
    REAL(KIND = DP) :: rek
    !! ReSigma at KS energy
    REAL(KIND = DP) :: imk
    !! ImSigma at KS energy
    REAL(KIND = DP) :: Gammak
    !! +/- ImSigma at KS energy
    REAL(KIND = DP) :: a_k
    !! alpha_k factor for A^{QP}
    REAL(KIND = DP) :: drek
    !! Derivative of ReSigma at ek
    REAL(KIND = DP) :: dimk
    !! Derivative of ImSigma at ek
    REAL(KIND = DP) :: d2imk
    !! Second-order derivative of ImSigma at ek
    REAL(KIND = DP) :: ww_as(nw_specfun)
    !! Auxiliary frequency grid to construct satellite function
    REAL(KIND = DP) :: a_s(nw_specfun)
    !! Temporary quantity needed to compute the satelite
    REAL(KIND = DP) :: conv(nw_specfun)
    !! Temporary quantity to compute the convolution
    REAL(KIND = DP) :: a_s2(nw_specfun), a_s3(nw_specfun)
    !! Second and third satellites
    !
    !  Extract indices for E_F (w=0) and the IP energy
    !
    i0 = MINLOC(ABS(ww(:)), DIM = 1)
    iek = MINLOC(ABS(ww(:) - ek), DIM = 1)
    !
    !  Construct auxiliary frequency grid for construction of satellite
    !  function, in which ww[iek] = ek is contained exactly.
    !
    IF ( ww(iek) > ek ) THEN
      !
      ww_as = ww - ww(iek) + ek
      !
    ELSE
      !
      ww_as = ww + ww(iek) - ek
      !
    ENDIF
    !
    !
    !  Initialise on-the-mass-shell self-energy: incoming self-energy
    !  is retarded ( Im(Sigma)>0 everywhere), need lesser and greater
    !  self-energy ( Im(Sigma)>0 for k>k_F, Im(Sigma)<0 for k<k_F)
    ! 
    rek = sigmar(iek)
    !
    IF (ek < 0) THEN            ! Lesser self-energy
      !
      sigmai(i0:) = zero
      imk = sigmai(iek)         ! This is Im Sigma^<(ek)
      Gammak = MAX(sigmai(iek),degaussw)      ! This is \Gamma in Eq.(51) in Ref.[1]
      !
    ELSE                        ! Greater self-energy
      !
      sigmai(i0:) = -sigmai(i0:)
      sigmai(1:i0) = zero      
      imk = sigmai(iek)         ! This is Im Sigma^>(ek)
      Gammak = MAX(-sigmai(iek),degaussw)     ! This is \Gamma in Eq.(51) in Ref.[1]
      !
    ENDIF
    !
    !
    !  Calculate finite difference derivatives of Re and Im Sigma
    !
    dw = ABS(ww_as(1)-ww_as(2))
    !
    IF (((iek-1) > 0) .AND. ((iek+1) <= nw_specfun)) THEN
      drek = (sigmar(iek+1) - sigmar(iek-1)) / (two * dw)
      dimk = (sigmai(iek+1) - sigmai(iek-1)) / (two * dw)
      d2imk = (sigmai(iek+1) - two * sigmai(iek) + sigmai(iek-1)) / &
              (dw**two)
    ELSE
      ! If iek+1 (or iek-1) is out of bounds, 
      ! skip calculating a_ce, a_qp, a_s1
      zeta = zero
      a_ce(:) = zero
      a_qp(:) = zero
      a_s1(:) = zero
      RETURN
    ENDIF
    !
    !
    !  Calculate Z_k and alpha_k factors
    !
    zeta = DEXP(drek)
    a_k = -dimk
    !
    !
    !  Calculate A^{QP}
    !
    a_qp = zero
    !
    DO iw = 1, nw_specfun
      !
      a_qp(iw) = zeta / pi * &
              (Gammak * COS(a_k) - (ww(iw) - ek - rek) * SIN(a_k)) / &
              ((ww(iw) - ek - rek)**two + Gammak**two)
      !
    ENDDO
    !
    !
    ! Calculate A^{S}
    !
    a_s = zero
    !
    DO iw = 1, nw_specfun
      !
      IF (iw == iek) THEN
        !
        a_s(iw) = one / two / pi * d2imk
        !
      ELSE
        !
        IF (ek < 0) THEN  ! Lesser self-energy
          !
          ! This is Eq.(67) in [1] for k<k_F
          !
          !S. Tiwari; Quick fix for out of bound error in GNU compiler
          ! 
          IF (((iw-iek+i0) > 0) .AND. ((iw-iek+i0) <= nw_specfun))THEN 
            !
            IF (sigmai(iw-iek+i0) > 0.d0) THEN
              a_s(iw) = (sigmai(iw) - (imk + (ww_as(iw) - ek) * dimk ) ) / &
                        (ww_as(iw) - ek)**two / pi
            ENDIF
            !
          ENDIF
          !
          ELSE ! Greater self-energy
          !
          ! This is Eq.(67) in [1] for k>k_F
          !
          IF (((iw-i0+iek) > 0) .AND. ((iw-i0+iek) <= nw_specfun))THEN           
            !
            IF (sigmai(iw-i0+iek) < 0.d0) THEN
              a_s(iw) = (-sigmai(iw) - (-imk - (ww_as(iw) - ek) * dimk ) ) / &
                      (ww_as(iw) - ek)**two / pi
            ENDIF
            ! 
          ENDIF
          !
        ENDIF
        !
      ENDIF
      !
    ENDDO
    !
    ! Calculate first satellite
    !
    a_s1 = zero
    !
    DO iw = 1, nw_specfun
      !
      DO iw2 = 1, nw_specfun
        !
        icoul = i0 + iw - iw2
        !
        IF ((icoul > 0) .AND. (icoul <= nw_specfun)) THEN
          a_s1(iw) = a_s1(iw) + a_qp(iw2) * a_s(icoul) * dw
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    ! Calculate second satellite
    !
    a_s2 = zero
    !
    DO iw = 1, nw_specfun
      !
      DO iw2 = 1, nw_specfun
        !
        icoul = i0 + iw - iw2
        !
        IF ((icoul > 0) .AND. (icoul <= nw_specfun)) THEN
          a_s2(iw) = a_s2(iw) + a_s1(iw2) * a_s(icoul) * dw
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    ! Calculate third satellite
    !
    a_s3 = zero
    !
    DO iw = 1, nw_specfun
      !
      DO iw2 = 1, nw_specfun
        !
        icoul = i0 + iw - iw2
        !
        IF ((icoul > 0) .AND. (icoul <= nw_specfun)) THEN
          a_s3(iw) = a_s3(iw) + a_s2(iw2) * a_s(icoul) * dw
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    a_ce = zero
    !
    a_ce = a_qp + a_s1 + a_s2 / two + a_s3 / 6.d0
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cumulant_conv
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE cum_mod
  !-----------------------------------------------------------------------
