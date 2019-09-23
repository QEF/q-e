  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  MODULE plot
  !----------------------------------------------------------------------
  !! 
  !! This module contains routine to plot data as well as DOS-like quantities. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE a2f_main()
    !-----------------------------------------------------------------------
    !!
    !! Compute the Eliasberg spectral function
    !! in the Migdal approximation. 
    !! 
    !! If the q-points are not on a uniform grid (i.e. a line)
    !! the function will not be correct
    !! 
    !! 02/2009 works in serial on ionode at the moment.  can be parallelized
    !! 03/2009 added transport spectral function -- this involves a v_k dot v_kq term 
    !!         in the quantities coming from selfen_phon.f90.  Not fully implemented  
    !! 10/2009 the code is transitioning to 'on-the-fly' phonon selfenergies
    !!         and this routine is not currently functional
    !! 10/2015 RM: added calcution of Tc based on Allen-Dynes formula 
    !! 09/2019 SP: Cleaning
    !!
    !
    USE kinds,     ONLY : DP
    USE phcom,     ONLY : nmodes
    USE cell_base, ONLY : omega
    USE epwcom,    ONLY : degaussq, delta_qsmear, nqsmear, nqstep, nsmear, eps_acustic, & 
                          delta_smear, degaussw, fsthick, nc
    USE elph2,     ONLY : nqtotf, wf, wqf, lambda_all, lambda_v_all
    USE constants_epw, ONLY : ryd2mev, ryd2ev, kelvin2eV, two, zero, kelvin2Ry, pi
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iua2ffil, iudosfil, iua2ftrfil, iures
    USE io_files,  ONLY : prefix
    ! 
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: fila2f_suffix
    !! FIXME
    CHARACTER(LEN = 256) :: fila2ftr
    !! FIXME
    CHARACTER(LEN = 256) :: fildos
    !! FIXME
    CHARACTER(LEN = 256) :: filres
    !! FIXME
    INTEGER :: imode
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: iw
    !! FIXME
    INTEGER :: ismear
    !! FIXME
    INTEGER :: isig
    !! FIXME
    INTEGER :: i
    !! FIXME
    INTEGER :: itemp
    !! FIXME
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: weight
    !! FIXME  
    REAL(KIND = DP) :: temp
    !! FIXME  
    REAL(KIND = DP) :: n
    !! FIXME  
    REAL(KIND = DP) :: be
    !! FIXME  
    REAL(KIND = DP) :: prefact
    !! FIXME  
    REAL(KIND = DP) :: lambda_tot
    !! FIXME  
    REAL(KIND = DP) :: lambda_tr_tot
    !! FIXME  
    REAL(KIND = DP) :: iomega
    !!
    REAL(KIND = DP) :: sigma
    !! 
    REAL(KIND = DP) :: a2f_tmp
    !! 
    REAL(KIND = DP) :: a2f_tr_tmp
    !! 
    REAL(KIND = DP) :: om_max
    !! 
    REAL(KIND = DP) :: dw
    !! 
    REAL(KIND = DP) :: w0
    !! 
    REAL(KIND = DP) :: l
    !! 
    REAL(KIND = DP) :: l_tr
    !! 
    REAL(KIND = DP) :: tc
    !! 
    REAL(KIND = DP) :: mu
    !! 
    REAL(KIND = DP), ALLOCATABLE :: a2fct(:, :), a2f_tr(:, :), l_a2f(:), l_a2f_tr(:), dosph(:, :), logavg(:), rho(:, :)
    !! 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! 
    !
    CALL start_clock('a2F')
    IF (mpime == ionode_id) THEN
      !
      ALLOCATE(a2fct(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating a2fct', 1)
      ALLOCATE(a2F_tr(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating a2f_tr', 1)
      ALLOCATE(dosph(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating dosph', 1)
      ALLOCATE(l_a2F(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating l_a2f', 1)
      ALLOCATE(l_a2F_tr(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating l_a2f_tr', 1)
      ALLOCATE(logavg(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating logavg', 1)
      ! The resitivity is computed for temperature between 0K-1000K by step of 10
      ! This is hardcoded and needs to be changed here if one wants to modify it
      ALLOCATE(rho(100, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f', 'Error allocating rho', 1)
      !  
      DO isig = 1, nsmear
        !
        IF (isig < 10) THEN
          WRITE(fila2f_suffix, '(a,a6,i1)') TRIM(prefix), '.a2f.0', isig
        ELSE 
          WRITE(fila2f_suffix, '(a,a5,i2)') TRIM(prefix), '.a2f.', isig
        ENDIF
        OPEN(UNIT = iua2ffil, FILE = fila2f_suffix, FORM = 'formatted')
        !
        IF (isig < 10) THEN
          WRITE(fila2ftr, '(a,a9,i1)') TRIM(prefix),'.a2f_tr.0', isig
        ELSE
          WRITE(fila2ftr, '(a,a8,i2)') TRIM(prefix),'.a2f_tr.', isig
        ENDIF
        OPEN(UNIT = iua2ftrfil, FILE = fila2ftr, FORM = 'formatted')
        !
        IF (isig < 10) THEN
          WRITE(filres, '(a,a6,i1)') TRIM(prefix), '.res.0', isig
        ELSE
          WRITE(filres, '(a,a5,i2)') TRIM(prefix), '.res.', isig
        ENDIF
        OPEN(UNIT = iures, FILE = filres, FORM = 'formatted')
        !
        IF (isig < 10) THEN
          WRITE(fildos, '(a,a8,i1)') TRIM(prefix), '.phdos.0', isig
        ELSE
          WRITE(fildos, '(a,a7,i2)') TRIM(prefix), '.phdos.', isig
        ENDIF
        OPEN(UNIT = iudosfil, FILE = fildos, FORM = 'formatted')
        !
        WRITE(stdout, '(/5x,a)') REPEAT('=',67)
        WRITE(stdout, '(5x,"Eliashberg Spectral Function in the Migdal Approximation")') 
        WRITE(stdout, '(5x,a/)') REPEAT('=',67)
        !
        om_max = 1.1d0 * MAXVAL(wf(:, :)) ! increase by 10%
        dw = om_max / DBLE(nqstep)
        !
        lambda_tot    = zero
        l_a2f(:)      = zero
        a2fct(:, :)     = zero
        lambda_tr_tot = zero
        l_a2f_tr(:)   = zero
        a2f_tr(:, :)  = zero
        dosph(:, :)   = zero
        logavg(:)     = zero
        !
        DO ismear = 1, nqsmear
          !
          sigma = degaussq + (ismear - 1) * delta_qsmear
          !
          DO iw = 1, nqstep  ! loop over points on the a2F(w) graph
            !
            iomega = DBLE(iw) * dw ! step through the frequncies we wish to plot
            !
            DO iq = 1, nqtotf ! loop over q-points 
              DO imode = 1, nmodes ! loop over modes
                w0 = wf(imode, iq)
                !
                IF (w0 > eps_acustic) THEN 
                  !
                  l  = lambda_all(imode, iq, isig)
                  IF (lambda_all(imode, iq, isig) < 0.d0)  l = 0.d0 ! sanity check
                  ! 
                  a2f_tmp    = wqf(iq) * w0 * l / two
                  !
                  weight = w0gauss((iomega - w0) / sigma, 0) / sigma
                  a2fct(iw, ismear) = a2fct(iw, ismear) + a2f_tmp * weight
                  dosph(iw, ismear)  = dosph(iw, ismear) + wqf(iq) * weight
                  !
                  l_tr = lambda_v_all(imode, iq, isig)
                  IF (lambda_v_all(imode, iq, isig) < 0.d0)  l_tr = 0.d0 !sanity check
                  ! 
                  a2f_tr_tmp = wqf(iq) * w0 * l_tr / two
                  !
                  a2f_tr(iw, ismear) = a2f_tr(iw, ismear) + a2f_tr_tmp * weight
                  !
                ENDIF
              ENDDO
            ENDDO
            !
            ! output a2F
            !
            IF (ismear == nqsmear) WRITE(iua2ffil,   '(f12.7, 15f12.7)') iomega * ryd2mev, a2fct(iw, :)
            IF (ismear == nqsmear) WRITE(iua2ftrfil, '(f12.7, 15f12.7)') iomega * ryd2mev, a2f_tr(iw, :)
            IF (ismear == nqsmear) WRITE(iudosfil,   '(f12.7, 15f12.7)') iomega * ryd2mev, dosph(iw, :) / ryd2mev
            !
            ! do the integral 2 int (a2F(w)/w dw)
            !
            l_a2f(ismear) = l_a2f(ismear) + two * a2fct(iw, ismear) / iomega * dw
            l_a2f_tr(ismear) = l_a2f_tr(ismear) + two * a2f_tr(iw, ismear) / iomega * dw
            logavg(ismear) = logavg(ismear) + two *  a2fct(iw, ismear) * LOG(iomega) / iomega * dw
            !
          ENDDO
          !
          logavg(ismear) = EXP(logavg(ismear) / l_a2f(ismear))
          !
        ENDDO
        !
        DO iq = 1, nqtotf ! loop over q-points 
          DO imode = 1, nmodes ! loop over modes
            IF (lambda_all(imode, iq, isig) > 0.d0 .AND. wf(imode, iq) > eps_acustic ) & 
               lambda_tot = lambda_tot + wqf(iq) * lambda_all(imode, iq, isig)
            IF (lambda_v_all(imode, iq, isig) > 0.d0 .AND. wf(imode, iq) > eps_acustic) &
               lambda_tr_tot = lambda_tr_tot + wqf(iq) * lambda_v_all(imode, iq, isig)
          ENDDO
        ENDDO
        WRITE(stdout, '(5x,a,f12.7)') "lambda : ", lambda_tot
        WRITE(stdout, '(5x,a,f12.7)') "lambda_tr : ", lambda_tr_tot
        WRITE(stdout, '(a)') " "
        !
        !
        ! Allen-Dynes estimate of Tc for ismear = 1
        !
        WRITE(stdout, '(5x,a,f12.7,a)') "Estimated Allen-Dynes Tc"
        WRITE(stdout, '(a)') " "
        WRITE(stdout, '(5x,a,f12.7,a,f12.7)') "logavg = ", logavg(1), " l_a2F = ", l_a2f(1)
        DO i = 1, 6
          !
          mu = 0.1d0 + 0.02d0 * DBLE(i - 1)
          tc = logavg(1) / 1.2d0 * EXP(-1.04d0 * (1.d0 + l_a2F(1)) / (l_a2f(1) - mu * ( 1.d0 + 0.62d0 * l_a2f(1))))
          ! tc in K
          !
          tc = tc * ryd2ev / kelvin2eV
          !SP: IF Tc is too big, it is not physical
          IF (tc < 1000.0) THEN
            WRITE(stdout, '(5x,a,f6.2,a,f22.12,a)') "mu = ", mu, " Tc = ", tc, " K"
          ENDIF 
          !
        ENDDO
        ! 
        rho(:, :) = zero
        ! Now compute the Resistivity of Metal using the Ziman formula
        ! rho(T,smearing) = 4 * pi * me/(n * e**2 * kb * T) int dw hbar w a2F_tr(w,smearing) n(w,T)(1+n(w,T))
        ! n is the number of electron per unit volume and n(w,T) is the Bose-Einstein distribution
        ! Usually this means "the number of electrons that contribute to the mobility" and so it is typically 8 (full shell)
        ! but not always. You might want to check this. 
        ! 
        n = nc / omega
        WRITE(iures, '(a)') '# Temperature [K]                &
                            &Resistivity [micro Ohm cm] for different Phonon smearing (meV)        '  
        WRITE(iures, '("#     ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1,nqsmear)
        DO ismear = 1, nqsmear
          DO itemp = 1, 100 ! Per step of 10K
            temp = itemp * 10 * kelvin2Ry
            ! omega is the volume of the primitive cell in a.u.  
            ! 
            prefact = 4.0 * pi / ( temp * n )
            DO iw = 1, nqstep  ! loop over points on the a2F(w)
              ! 
              iomega = DBLE(iw) * dw
              be = 1.0 / (EXP(iomega / temp) - 1); 
              ! Perform the integral with rectangle. 
              rho(itemp, ismear) = rho(itemp, ismear) + prefact * iomega * a2f_tr(iw, ismear) * be * (1.0 + be) * dw  
              ! 
            ENDDO
            ! From a.u. to micro Ohm cm
            ! Conductivity 1 a.u. = 2.2999241E6 S/m
            ! Now to go from Ohm*m to micro Ohm cm we need to multiply by 1E8 
            rho(itemp, ismear) = rho(itemp, ismear) * 1E8 / 2.2999241E6
            IF (ismear == nqsmear) WRITE (iures, '(i8, 15f12.7)') itemp * 10, rho(itemp, :)
          ENDDO
        ENDDO 
        CLOSE(iures)
        !
        WRITE(iua2ffil, *) "Integrated el-ph coupling"
        WRITE(iua2ffil, '("  #         ", 15f12.7)') l_a2f(:)
        WRITE(iua2ffil, *) "Phonon smearing (meV)"
        WRITE(iua2ffil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
        WRITE(iua2ffil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
        WRITE(iua2ffil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
        WRITE(iua2ffil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
        CLOSE(iua2ffil)
        !
        WRITE(iua2ftrfil, *) "Integrated el-ph coupling"
        WRITE(iua2ftrfil, '("  #         ", 15f12.7)') l_a2f_tr(:)
        WRITE(iua2ftrfil, *) "Phonon smearing (meV)"
        WRITE(iua2ftrfil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
        WRITE(iua2ftrfil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
        WRITE(iua2ftrfil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
        WRITE(iua2ftrfil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
        CLOSE(iua2ftrfil)
        !
        CLOSE(iudosfil)
        !
      ENDDO ! isig
      ! 
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2F', 1)
      DEALLOCATE(l_a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2F_tr', 1)
      DEALLOCATE(a2fct, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2F', 1)
      DEALLOCATE(a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2F_tr', 1)
      DEALLOCATE(rho, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating rho', 1)
      DEALLOCATE(dosph, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating dosph', 1)
      DEALLOCATE(logavg, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating logavg', 1)
      !
    ENDIF
    !
    CALL stop_clock('a2F')
    CALL print_clock('a2F')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE a2f_main
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE nesting_fn_q(iqq, iq)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation. This corresponds to 
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation. 
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : nbndsub, fsthick, &
                          eptemp, ngaussw, degaussw,     &
                          nsmear, delta_smear, efermi_read, fermi_energy
    USE pwcom,     ONLY : nelec, ef
    USE klist_epw, ONLY : isk_dummy
    USE elph2,     ONLY : ibndmax, ibndmin, etf, &
                          wkf, xqf, wqf, nkqf, nktotf, &
                          nkf, nkqtotf, xqf, nbndfst
    USE constants_epw, ONLY : ryd2ev, two
    USE mp,        ONLY : mp_barrier,mp_sum
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    ! 
    ! Local variables
    INTEGER :: ik
    !! Counter on the k-point index 
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index 
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: ismear
    !! Smearing for the Gaussian function 
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: weight
    !! Imaginary part of the phonhon self-energy factor 
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: w0g1
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0gauss
    !! 
    REAL(KIND = DP) :: dos_ef
    !! 
    REAL(KIND = DP) :: gamma
    !! 
    REAL(KIND = DP) :: degaussw0
    !! 
    REAL(KIND = DP), EXTERNAL :: efermig
    !! 
    !
    IF (iqq == 1) THEN
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Nesting Function in the double delta approx")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x,a,f10.6,a)' ) &
        'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
    ENDIF
    !
    ! SP: The Gamma function needs to be put to 0 for each q
    gamma = 0.0
    ! 
    ! Here we loop on smearing values
    DO ismear = 1, nsmear
      !
      degaussw0 = (ismear - 1) * delta_smear + degaussw
      !
      ! Fermi level and corresponding DOS
      !
      !   Note that the weights of k+q points must be set to zero here
      !   no spin-polarized calculation here
      IF (efermi_read) THEN
        ef0 = fermi_energy 
      ELSE
        ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk_dummy)
      ENDIF
      !
      dosef = dos_ef(ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
      !  N(Ef) in the equation for lambda is the DOS per spin
      dosef = dosef / two
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw0 * ryd2ev, ngaussw
        WRITE(stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
      ENDIF
      !
      !
      CALL start_clock('nesting')
      !
      fermicount = 0
      !
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) then
          !
          fermicount = fermicount + 1
          !
          DO ibnd = 1, nbndfst
            !
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            w0g1 = w0gauss(ekk / degaussw0, 0) / degaussw0
            !
            DO jbnd = 1, nbndfst
              !
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              w0g2 = w0gauss(ekq / degaussw0, 0) / degaussw0
              !
              ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
              ! This is the imaginary part of the phonon self-energy, sans the matrix elements
              !
              ! weight = wkf (ikk) * (wgkk - wgkq) * &
              !      aimag ( cone / ( ekq - ekk  - ci * degaussw ) ) 
              !
              ! the below expression is positive-definite, but also an approximation
              ! which neglects some fine features
              !
              weight = wkf(ikk) * w0g1 * w0g2
              !
              gamma  = gamma  + weight  
              !
            ENDDO ! jbnd
          ENDDO   ! ibnd
          !
        ENDIF ! endif fsthick
        !
      ENDDO ! loop on k
      !
      ! collect contributions from all pools (sum over k-points)
      ! this finishes the integral over the BZ  (k)
      !
      CALL mp_sum(gamma, inter_pool_comm) 
      CALL mp_sum(fermicount, inter_pool_comm)
      !
      WRITE(stdout, '(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:, iq) , wqf(iq)
      WRITE(stdout, '(5x,a)') REPEAT('-',67)
         ! 
      WRITE(stdout, 102) gamma
      WRITE(stdout, '(5x,a/)') REPEAT('-',67)
      !
      WRITE(stdout, '(/5x,a,i8,a,i8/)') &
        'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', nktotf
      !
      CALL stop_clock('nesting')
    ENDDO !smears
    !
100 FORMAT(5x, 'Gaussian Broadening: ',f7.3,' eV, ngauss=', i4)
101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
102 FORMAT(5x, 'Nesting function (q)=', e15.6, ' [Adimensional]')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE nesting_fn_q
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE plot_band()
    !-----------------------------------------------------------------------
    !!
    !! This routine writes output files for phonon dispersion and band structure 
    !! SP : Modified so that it works with the current plotband.x of QE 5
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, filqf, filkf
    USE elph2,         ONLY : etf, nkf, nqtotf, wf, xkf, xqf, nkqtotf, nktotf
    USE constants_epw, ONLY : ryd2mev, ryd2ev
    USE io_var,        ONLY : iufilfreq, iufileig
    USE elph2,         ONLY : nkqf
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! Global k-point index
    INTEGER :: ikk
    !! Index for the k-point
    INTEGER :: ikq
    !! Index for the q-point
    INTEGER :: ibnd
    !! Band index
    INTEGER :: imode
    !! Mode index
    INTEGER :: iq
    !! Global q-point index
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: dist
    !! Distance from Gamma
    REAL(KIND = DP) :: dprev
    !! Previous distance
    REAL(KIND = DP) :: dcurr
    !! Current distance
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! K-points on the full k grid (all pools)
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigenenergies on the full k grid (all pools)
    !
    IF (filqf /= ' ') THEN
      ! 
      IF (my_pool_id == ionode_id) THEN
        !
        OPEN(iufilfreq, FILE = "phband.freq", FORM = 'formatted')
        WRITE(iufilfreq, '(" &plot nbnd=",i4,", nks=",i6," /")') nmodes, nqtotf
        !
        ! crystal to cartesian coordinates
        CALL cryst_to_cart(nqtotf, xqf, bg, 1)
        !
        dist = 0.d0
        dprev = 0.d0
        dcurr = 0.d0
        DO iq = 1, nqtotf
          !
          IF (iq /= 1) THEN  
            dist = SQRT((xqf(1, iq) - xqf(1, iq - 1)) * (xqf(1, iq) - xqf(1, iq - 1)) & 
                      + (xqf(2, iq) - xqf(2, iq - 1)) * (xqf(2, iq) - xqf(2, iq - 1)) & 
                      + (xqf(3, iq) - xqf(3, iq - 1)) * (xqf(3, iq) - xqf(3, iq - 1)))
          ELSE 
            dist = 0.d0
          ENDIF
          dcurr = dprev + dist
          dprev = dcurr
          WRITE(iufilfreq, '(10x,3f10.6)') xqf(:, iq)
          WRITE(iufilfreq, '(1000f14.4)') (wf(imode, iq) * ryd2mev, imode = 1, nmodes)
          !
        ENDDO
        CLOSE(iufilfreq)
        !
        ! back from cartesian to crystal coordinates
        CALL cryst_to_cart(nqtotf, xqf, at, -1)
        !
      ENDIF
    ENDIF ! filqf
    ! 
    IF (filkf /= ' ') THEN
      !
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
      ENDDO
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error allocating etf_all', 1)
      !
#if defined(__MPI)
      CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_barrier(inter_pool_comm)
#else    
      !
      xkf_all = xkf
      etf_all = etf
#endif
      !
      IF (my_pool_id == ionode_id) THEN
        !
        OPEN(iufileig, FILE = "band.eig", FORM = 'formatted')
        WRITE(iufileig, '(" &plot nbnd=",i4,", nks=",i6," /")') nbndsub, nktotf
        !
        ! crystal to cartesian coordinates
        CALL cryst_to_cart(nkqtotf, xkf_all, bg, 1)
        !
        dist = 0.d0
        dprev = 0.d0
        dcurr = 0.d0
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF (ikk /= 1) THEN
            dist = SQRT((xkf_all(1, ikk) - xkf_all(1, ikk - 2)) * (xkf_all(1, ikk) - xkf_all(1, ikk - 2)) &
                      + (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) * (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) &
                      + (xkf_all(3, ikk) - xkf_all(3, ikk - 2)) * (xkf_all(3, ikk) - xkf_all(3, ikk - 2)))
          ELSE
            dist = 0.d0
          ENDIF
          dcurr = dprev + dist
          dprev = dcurr
          WRITE(iufileig, '(10x,3f10.6)') xkf_all(:, ikk)
          WRITE(iufileig, '(1000f20.12)') (etf_all(ibnd, ikk) * ryd2ev, ibnd = 1, nbndsub)
          !
        ENDDO
        CLOSE(iufileig)
        !
        ! back from cartesian to crystal coordinates
        CALL cryst_to_cart(nkqtotf, xkf_all, at, -1)
        !
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating etf_all', 1)
      !
    ENDIF ! filkf
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE plot_band
    !----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE plot
  !------------------------------------------------------------------------------
