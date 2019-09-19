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
  MODULE spectral_func
  !----------------------------------------------------------------------
  !! 
  !! This module contains the various spectral function routines
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_q(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation. 
    !!  
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically. 
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iospectral_sup ,iospectral
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, eps_acustic, &
                              fsthick, eptemp, ngaussw, degaussw, wmin_specfun, &
                              wmax_specfun, nw_specfun, shortrange, &
                              efermi_read, fermi_energy
    USE pwcom,         ONLY : nelec, ef
    USE klist_epw,     ONLY : isk_dummy
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, nktotf, &
                              epf17, wkf, nkf, wf, wqf, xkf, nkqtotf, adapt_smearing, &
                              esigmar_all, esigmai_all, a_all, nbndfst
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : me_pool, inter_pool_comm
    USE division,      ONLY : fkbounds
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index in selecq  
    INTEGER, INTENT(in) :: iq
    !! Current q-point index  
    INTEGER, INTENT(in) :: totq
    !! Total number of q-point in window
    ! 
    ! Local variables
    INTEGER :: iw
    !! Counter on the frequency
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
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgq
    !! Bose occupation factor $n_{q\nu}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Self-energy factor 
    !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$ 
    !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$
    REAL(KIND = DP) :: inv_wq
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: g2_tmp
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons  
    REAL(KIND = DP) :: ww
    !! Current frequency
    REAL(KIND = DP) :: dw 
    !! Frequency intervals
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP) :: esigmar0
    !! 
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! 
    REAL(KIND = DP), EXTERNAL :: efermig
    !! 
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! 
    REAL(KIND = DP), EXTERNAL :: wgauss
    !!
    ! 
    IF (adapt_smearing) CALL errore('spectral_func_q', 'adapt_smearing cannot be used with spectral functions ', 1)
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_eptemp0 = 1.0 / eptemp
    inv_degaussw = 1.0 / degaussw
    !   
    ! energy range and spacing for spectral function
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1.0d0)
    !
    IF (iqq == 1) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Electron Spectral Function in the Migdal Approximation")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
      !
    ENDIF
    !
    ! Fermi level and corresponding DOS
    !
    IF (efermi_read) THEN
      !
      ef0 = fermi_energy
      !
    ELSE
      !
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated in ephwann_shuffle
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
      !
    ENDIF
    !
    IF (iq == 1) THEN 
      WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    ! SP: Sum rule added to conserve the number of electron. 
    IF (iq == 1) THEN
      WRITE(stdout, '(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
      WRITE(stdout, '(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
      WRITE(stdout, '(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
    ! loop over all k points of the fine mesh
    !
    fermicount = 0 
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      ! (but in this case they are the same)
      !
      IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
          (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
          !
          ! the phonon frequency and Bose occupation
          wq = wf(imode, iq)
          ! SP: Define the inverse for efficiency
          inv_wq = 1.0 / (two * wq)
          wgq = wgauss(-wq * inv_eptemp0, -99)
          wgq = wgq / (one - two * wgq)
          !
          ! SP: Avoid if statement in inner loops
          IF (wq > eps_acustic) THEN
            g2_tmp = 1.0
          ELSE
            g2_tmp = 0.0
          ENDIF           
          !
          DO ibnd = 1, nbndfst
            !
            ! the energy of the electron at k (relative to Ef)
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !  
            DO jbnd = 1, nbndfst
              !
              ! the fermi occupation for k+q
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              wgkq = wgauss(-ekq / eptemp, -99)  
              !
              ! here we take into account the zero-point SQRT(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              IF (shortrange .AND. (ABS(xqf(1, iq))> eps8 .OR. ABS(xqf(2, iq))> eps8 &
                 .OR. ABS(xqf(3, iq))> eps8)) THEN
                ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                !     number, in which case its square will be a negative number. 
                g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq * g2_tmp)
              ELSE
                g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
              ENDIF
              !
              DO iw = 1, nw_specfun
                !
                ww = wmin_specfun + DBLE(iw - 1) * dw
                !
                weight = wqf(iq) * REAL(((wgkq + wgq) / (ww - (ekq - wq) - ci * degaussw)  +  &
                    (one - wgkq + wgq) / (ww - (ekq + wq) - ci * degaussw)))
                !
                esigmar_all(ibnd, ik + lower_bnd - 1, iw) = esigmar_all(ibnd, ik + lower_bnd - 1, iw) + g2 * weight 
                ! 
                ! SP : Application of the sum rule
                esigmar0 =  g2 *  wqf(iq) * REAL(((wgkq + wgq) / (-(ekq - wq) - ci * degaussw)  +  &
                    (one - wgkq + wgq) / (-(ekq + wq) - ci * degaussw)))
                esigmar_all(ibnd, ik + lower_bnd - 1, iw) = esigmar_all(ibnd, ik + lower_bnd - 1, iw) - esigmar0
                !
                weight = wqf(iq) * AIMAG(((wgkq + wgq) / (ww - (ekq - wq) - ci * degaussw)  +  &
                    (one - wgkq + wgq) / (ww - (ekq + wq) - ci * degaussw)))
                !
                esigmai_all(ibnd, ik + lower_bnd - 1, iw) = esigmai_all(ibnd, ik + lower_bnd - 1, iw) + g2 * weight
                !
              ENDDO
            ENDDO !jbnd
          ENDDO !ibnd
        ENDDO !imode
      ENDIF ! endif  fsthick
    ENDDO ! end loop on k
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      !
      ! Output electron spectral function here after looping over all q-points 
      ! (with their contributions summed in a etc.)
      !
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      ! construct the trace of the spectral function (assume diagonal selfenergy
      ! and constant matrix elements for dipole transitions)
      !
      IF (me_pool == 0) then
        OPEN(UNIT = iospectral, FILE = 'specfun.elself') 
        OPEN(UNIT = iospectral_sup, FILE = 'specfun_sup.elself') 
      ENDIF
      IF (me_pool == 0) then
        WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
        WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)' 
      ENDIF
      IF (me_pool == 0) then
        WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
        WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]      Real Sigma[meV]  Im Sigma[meV]'
      ENDIF
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout, '(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf_all (:, ikk)
        WRITE(stdout, '(5x,a)') REPEAT('-',67)
        !
        DO iw = 1, nw_specfun
          !
          ww = wmin_specfun + DBLE(iw - 1) * dw
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            a_all(iw, ik) = a_all(iw, ik) + ABS(esigmai_all(ibnd, ik, iw)) / pi / &
                 ((ww - ekk - esigmar_all(ibnd, ik, iw))**two + (esigmai_all(ibnd, ik, iw))**two)
            !
          ENDDO
          !
          WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw, ik) / ryd2mev
          !
        ENDDO
        !
        WRITE(stdout, '(5x,a/)') REPEAT('-',67)
        !
      ENDDO
      !
      DO ik = 1, nktotf
        !
        ! The spectral function should integrate to 1 for each k-point
        specfun_sum = 0.0
        ! 
        DO iw = 1, nw_specfun
          !
          ww = wmin_specfun + DBLE(iw - 1) * dw
          fermi(iw) = wgauss(-ww / eptemp, -99) 
          !
          specfun_sum = specfun_sum + a_all(iw, ik) * fermi(iw) * dw 
          !
        IF (me_pool == 0) WRITE(iospectral, '(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw, ik) / ryd2mev
           !
        ENDDO
        !
        IF (me_pool == 0) WRITE(iospectral, '(a)') ' '
        IF (me_pool == 0) WRITE(iospectral, '(2x,a,2x,e12.5)') '# Integrated spectral function ', specfun_sum
        !
      ENDDO
      !
      IF (me_pool == 0) CLOSE(iospectral)
      !
      DO ibnd = 1, nbndfst
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          !  the energy of the electron at k
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
          !
          DO iw = 1, nw_specfun
            !
            ww = wmin_specfun + DBLE(iw - 1) * dw
            WRITE(stdout, '(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik, &
              ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd, ik, iw), &
              ryd2mev * esigmai_all(ibnd, ik, iw)
            ! 
            IF (me_pool == 0) &
            WRITE(iospectral_sup, '(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
              ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd, ik, iw), &
              ryd2mev * esigmai_all(ibnd, ik, iw)
            !
          ENDDO
          !
        ENDDO
        !
        WRITE(stdout,*) ' '
        !
      ENDDO
      !
      IF (me_pool == 0) CLOSE(iospectral_sup)
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_q', 'Error deallocating etf_all', 1)
      !
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    103 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', e12.5, ' meV^-1')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_q
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_ph(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation. This corresponds to 
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation.  This routine is similar to the one above
    !! but it is ONLY called from within ephwann_shuffle and calculates 
    !! the selfenergy for one phonon at a time.  Much smaller footprint on the disk
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iospectral_sup, iospectral
    USE phcom,     ONLY : nmodes
    USE epwcom,    ONLY : nbndsub, fsthick, &
                          eptemp, ngaussw, degaussw, &
                          shortrange, nsmear, delta_smear, eps_acustic, &
                          efermi_read, fermi_energy, wmin_specfun,&
                          wmax_specfun, nw_specfun
    USE pwcom,     ONLY : nelec, ef
    USE klist_epw, ONLY : isk_dummy
    USE elph2,     ONLY : epf17, ibndmax, ibndmin, etf, nbndfst, &
                          wkf, xqf, nkqf, nkf, wf, a_all_ph, efnew
    USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi, cone, ci, eps8
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index 
    INTEGER, INTENT(in) :: totq
    !! Total q-points in selecq window
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
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: iw 
    !! Counter on frequency for the phonon spectra
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Self-energy factor 
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse Gaussian for efficiency reasons   
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse temperature
    REAL(KIND = DP) :: inv_wq
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons  
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: g2_tmp
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: gamma0(nmodes)
    !! Phonon self-energy
    REAL(KIND = DP) :: ww
    !! Current frequency
    REAL(KIND = DP) :: dw
    !! Frequency intervals 
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function to compute the Fermi energy 
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: gammai_all(nw_specfun, nmodes)
    !! Imaginary part of the frequency dependent spectral function
    REAL(KIND = DP) :: gammar_all(nw_specfun, nmodes)
    !!  Real part of the Phonon self-energy (freq. dependent for spectral function)
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    gammar_all(:, :) = zero
    gammai_all(:, :) = zero
    !
    ! Thomas-Fermi screening according to Resta PRB 1977
    ! Here specific case of Diamond
    !eps0   = 5.7
    !RTF    = 2.76 
    !qTF    = 1.36 
    !qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2) * tpiba2
    !epsTF =  (qTF**2 + qsquared) / (qTF**2/eps0 * sin (sqrt(qsquared)*RTF)/(sqrt(qsquared)*RTF)+qsquared)
    !
    IF (iqq == 1) THEN 
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Phonon Spectral Function Self-Energy in the Migdal Approximation (on the fly)")') 
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
    ENDIF
    !
    ! SP: Multiplication is faster than division ==> Important if called a lot in inner loops
    inv_degaussw = 1.0 / degaussw
    inv_eptemp   = 1.0 / eptemp
    !
    ! Fermi level and corresponding DOS
    !
    IF (efermi_read) THEN
      !
      ef0 = fermi_energy
      !
    ELSEIF (nsmear > 1) THEN
      !
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated in ephwann_shuffle
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
      !
    ELSE !SP: This is added for efficiency reason because the efermig routine is slow
      ef0 = efnew
    ENDIF
    !
    ! N(Ef) in the equation for lambda is the DOS per spin
    dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
    dosef = dosef / two
    !
    IF (iqq == 1) THEN 
      WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
      WRITE(stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
    ENDIF
    !
    CALL start_clock('PH SPECTRAL-FUNCTION')
    !
    fermicount = 0
    gamma0(:)  = zero
    !
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      ! Here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
          (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
        !
        fermicount = fermicount + 1
        !
        DO imode = 1, nmodes
          !
          ! the phonon frequency
          wq = wf(imode, iq)
          !
          ! SP : We should avoid branching statements (if statements) in innerloops. Therefore we do it here.
          inv_wq = 1.0 / (two * wq)
          ! the coupling from Gamma acoustic phonons is negligible
          IF (wq > eps_acustic) THEN
            g2_tmp = 1.0
          ELSE
            g2_tmp = 0.0
          ENDIF   
          !
          DO ibnd = 1, nbndfst
            !
            !  the fermi occupation for k
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            wgkk = wgauss(-ekk * inv_eptemp, -99)
            !
            DO jbnd = 1, nbndfst
              !
              !  the fermi occupation for k+q
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              wgkq = wgauss(-ekq * inv_eptemp, -99)  
              !
              ! here we take into account the zero-point SQRT(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              IF (shortrange .AND. (ABS(xqf(1, iq))> eps8 .OR. ABS(xqf(2, iq))> eps8 &
                 .OR. ABS(xqf(3, iq))> eps8 )) THEN
                ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                !     number, in which case its square will be a negative number. 
                g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq * g2_tmp) !* epsTF
              ELSE
                g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp !* epsTF
              ENDIF
              !
              ! SP - 03/2019 - Retarded phonon self-energy
              !                See Eq. 145 of RMP 89, 015003 (2017)
              ! \Pi^R = k-point weight * [ [f(E_k+q) - f(E_k)]/ [E_k+q - E_k -w_q - id] 
              !                           -[f(E_k+q) - f(E_k)]/ [E_k+q - E_k - id] ]
              ! The second term is gamma0 (static) 
              !
              weight = wkf(ikk) * (wgkq - wgkk) * REAL(cone / (ekq - ekk + ci * degaussw)) 
              !
              gamma0(imode) = gamma0(imode) + weight * g2
              ! 
              DO iw = 1, nw_specfun
                !
                ww = wmin_specfun + DBLE(iw - 1) * dw
                !
                weight = wkf(ikk) * (wgkq - wgkk) * REAL(cone / (ekq - ekk - ww + ci * degaussw)) 
                gammar_all(iw, imode) = gammar_all(iw, imode) + weight * g2
                !
                ! Normal implementation 
                !weight = wkf (ikk) * (wgkq - wgkk) * &
                !   aimag ( cone / ( ekq - ekk - ww + ci * degaussw0 ) ) 
                !  
                ! More stable:
                ! Analytical im. part 
                weight = pi * wkf(ikk) * (wgkq - wgkk) * w0gauss((ekq - ekk - ww) / degaussw, 0) / degaussw     
                !
                gammai_all(iw, imode) = gammai_all(iw, imode) + weight * g2 
                ! 
              ENDDO
            ENDDO ! jbnd
          ENDDO   ! ibnd
        ENDDO ! loop on q-modes
      ENDIF ! endif fsthick
    ENDDO ! loop on k
    !
    CALL stop_clock('PH SPECTRAL-FUNCTION')
    !
#if defined(__MPI)
    !
    ! collect contributions from all pools (sum over k-points) this finishes the integral over the BZ  (k)
    !
    CALL mp_sum(gammai_all,inter_pool_comm) 
    CALL mp_sum(gammar_all,inter_pool_comm) 
    CALL mp_sum(gamma0,inter_pool_comm) 
    CALL mp_sum(fermicount, inter_pool_comm)
    !
#endif
    !
    WRITE(stdout, '(5x,a)')
    !
    IF (iqq == 1) THEN
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iospectral, FILE = 'specfun.phon')
        OPEN(UNIT = iospectral_sup, FILE = 'specfun_sup.phon')
        WRITE(iospectral, '(/2x,a)') '#Phonon spectral function (meV)'
        WRITE(iospectral_sup, '(2x,a)') '#Phonon eigenenergies + real and im part of phonon self-energy (meV)'
        WRITE(iospectral, '(/2x,a)') '#Q-point    Energy[eV]     A(q,w)[meV^-1]'
        WRITE(iospectral_sup, '(2x,a)') '#Q-point    Mode       w_q[eV]        w[eV]    &
                                        &Real Sigma(w)[meV]   Real Sigma(w=0)[meV]     Im Sigma(w)[meV]'
      ENDIF
    ENDIF
    !
    ! Write to output file  
    WRITE(stdout, '(/5x,a)') 'Real and Imaginary part of the phonon self-energy (omega=0) without gamma0.'  
    DO imode = 1, nmodes
      wq = wf(imode, iq)
      ! Real and Im part of Phonon self-energy at 0 freq. 
      WRITE(stdout,105) imode, ryd2ev * wq, ryd2mev * gammar_all(1, imode), ryd2mev * gammai_all(1, imode)
    ENDDO 
    WRITE(stdout, '(5x,a,i8,a,i8)' ) 'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', totq
    !
    ! Write to support files
    DO iw = 1, nw_specfun
      !
      ww = wmin_specfun + DBLE(iw - 1) * dw
      !
      DO imode = 1, nmodes
        ! 
        wq = wf(imode, iq)
        !a_all(iw,iq) = a_all(iw,iq) + ABS(gammai_all(imode,iq,iw) ) / pi / &
        !      ( ( ww - wq - gammar_all (imode,iq,iw) + gamma0 (imode))**two + (gammai_all(imode,iq,iw) )**two )
        ! SP: From Eq. 16 of PRB 9, 4733 (1974)
        !    Also in Eq.2 of PRL 119, 017001 (2017). 
        a_all_ph(iw, iqq) = a_all_ph(iw, iqq) + (1.0d0 / pi) * ((2 * wq)**2) * ABS(gammai_all(iw, imode)) / &
              ((ww**2 - wq**2 - 2 * wq * (gammar_all(iw, imode) - gamma0(imode)))**two + (2.0d0 * wq * gammai_all(iw, imode))**two)
        !
        IF (mpime == ionode_id) THEN
          WRITE(iospectral_sup, '(2i9, 2x, f12.5, 2x, f12.5, 2x, E22.14, 2x, E22.14, 2x, E22.14)') iq, &
                imode, ryd2ev * wq, ryd2ev * ww, ryd2mev * gammar_all(iw, imode), ryd2mev * gamma0(imode), &
                ryd2mev * gammai_all(iw, imode)
        ENDIF
        !
      ENDDO 
      !
      IF (mpime == ionode_id) THEN 
        WRITE(iospectral, '(2x,i7,2x,f12.5,2x,E22.14)') iq, ryd2ev * ww, a_all_ph(iw, iqq) / ryd2mev ! print to file 
      ENDIF
      !
    ENDDO
    !
    IF (iqq == totq) THEN
      IF (mpime == ionode_id) THEN
        CLOSE(iospectral)
        CLOSE(iospectral_sup)
      ENDIF 
    ENDIF   
    WRITE(stdout, '(5x,a/)') REPEAT('-',67)
    ! 
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6,' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    105 FORMAT(5x, 'Omega( ', i3, ' )=', f9.4,' eV   Re[Pi]=', f15.6, ' meV Im[Pi]=', f15.6, ' meV')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_ph
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_pl_q(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation. 
    !!  
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically. 
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iospectral_sup, iospectral
    USE epwcom,        ONLY : nbndsub, &
                              fsthick, eptemp, ngaussw, degaussw, epsiHEG, &
                              wmax_specfun, nw_specfun, wmin_specfun,      &
                              efermi_read, fermi_energy, nel, meff
    USE pwcom,         ONLY : nelec, ef
    USE klist_epw,     ONLY : isk_dummy
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, nbndfst,   &
                              wkf, nkf, wqf, xkf, nkqtotf, xqf, dmef, &
                              esigmar_all, esigmai_all, a_all, nktotf
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : me_pool, inter_pool_comm 
    USE cell_base,     ONLY : omega, alat, bg
    USE division,      ONLY : fkbounds
    USE selfen,        ONLY : get_eps_mahan
    USE poolgathering, ONLY : poolgather2
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: iqq
    !! Q-point index in selecq
    INTEGER, INTENT(in) :: iq 
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in fsthick window
    !
    ! Local variables
    INTEGER :: iw
    !! 
    INTEGER :: ik 
    !! 
    INTEGER :: ikk 
    !! 
    INTEGER :: ikq
    !! 
    INTEGER :: ibnd
    !! 
    INTEGER :: jbnd
    !! 
    INTEGER :: fermicount
    !! 
    INTEGER :: lower_bnd
    !! 
    INTEGER :: upper_bnd
    !! 
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: g2
    !! 
    REAL(KIND = DP) :: ekk
    !! 
    REAL(KIND = DP) :: ekq
    !! 
    REAL(KIND = DP) :: wq
    !! 
    REAL(KIND = DP) :: ef0
    !! 
    REAL(KIND = DP) :: wgq
    !! 
    REAL(KIND = DP) :: wgkq
    !! 
    REAL(KIND = DP) :: ww
    !! Spectral frequency
    REAL(KIND = DP) :: dw
    !! Spectral frequency increment
    REAL(KIND = DP) :: weight
    !! 
    REAL(KIND = DP) :: specfun_sum
    !! 
    REAL(KIND = DP) :: esigmar0
    !! 
    REAL(KIND = DP) :: tpiba_new
    !! 
    REAL(KIND = DP) :: kF
    !! 
    REAL(KIND = DP) :: vF
    !! 
    REAL(KIND = DP) :: fermiHEG
    !! 
    REAL(KIND = DP) :: qin
    !! 
    REAL(KIND = DP) :: wpl0
    !! 
    REAL(KIND = DP) :: eps0
    !! 
    REAL(KIND = DP) :: deltaeps
    !! 
    REAL(KIND = DP) :: qcut
    !! 
    REAL(KIND = DP) :: qsquared
    !! 
    REAL(KIND = DP) :: qTF
    !! 
    REAL(KIND = DP) :: dipole
    !! 
    REAL(KIND = DP) :: rs
    !! 
    REAL(KIND = DP) :: ekk1
    !! 
    REAL(KIND = DP) :: degen
    !! 
    REAL(KIND = DP) :: fermi(nw_specfun)
    !!
    REAL(KIND = DP) :: q(3)
    !! The q-point in cartesian unit.
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !!  
    REAL(KIND = DP), EXTERNAL :: efermig
    !!
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! 
    REAL(KIND = DP), EXTERNAL :: wgauss
    !
    ! loop over temperatures can be introduced
    !
    ! energy range and spacing for spectral function
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    !
    IF (iqq == 1) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Electron Spectral Function in the Migdal Approximation")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
      !
    ENDIF
    !
    ! Fermi level and corresponding DOS
    !
    IF (efermi_read) THEN
      !
      ef0 = fermi_energy
      !
    ELSE
      !
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated in ephwann_shuffle
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
      !
    ENDIF
    !
    IF (iqq == 1) THEN 
      WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    ! SP: Sum rule added to conserve the number of electron. 
    IF (iqq == 1) THEN
      WRITE(stdout, '(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
      WRITE(stdout, '(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
      WRITE(stdout, '(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
!
    !nel      =  0.01    ! this should be read from input - # of doping electrons 
    !epsiHEG  =  12.d0   ! this should be read from input - # dielectric constant at zero doping  
    !meff     =  0.25    ! this should be read from input - effective mass 
    tpiba_new = 2.0d0 * pi / alat
    degen     = 1.0d0
    rs        = (3.d0 / (4.d0 * pi * nel / omega / degen))**(1.d0 / 3.d0) * meff * degen ! omega is the unit cell volume in Bohr^3
    kF        = (3.d0 * (pi**2) * nel / omega / degen)**(1.d0 / 3.d0)
    vF        = (1.d0 / meff) * (3.d0 * (pi**2) * nel / omega / degen)**(1.d0 / 3.d0)
    fermiHEG  = 1.d0 / (2.d0 * meff) * (3.d0 * (pi**2) * nel / omega / degen)**(2.d0 / 3.d0) * 2.d0 ! [Ryd] multiplication by 2 converts from Ha to Ry
    qTF       = (6.d0 * pi * nel / omega / degen / (fermiHEG / 2.d0))**(1.d0 / 2.d0)    ! [a.u.]
    wpl0      = SQRT(4.d0 * pi * nel / omega / meff / epsiHEG) * 2.d0         ! [Ryd] multiplication by 2 converts from Ha to Ryd
    wq        = wpl0 ! [Ryd] 
    q(:)      = xqf(:, iq) 
    CALL cryst_to_cart(1, q, bg, 1)
    qsquared  =  (q(1)**2 + q(2)**2 + q(3)**2)
    qin       =  SQRT(qsquared) * tpiba_new
    qcut      =  wpl0 / vF  / tpiba_new / 2.d0 ! 1/2 converts from Ryd to Ha
    ! 
    CALL get_eps_mahan(qin, rs, kF, eps0) ! qin should be in atomic units for Mahan formula
    deltaeps  = -(1.d0 / (epsiHEG + eps0 - 1.d0) - 1.d0 / epsiHEG)
    !
    IF (iq == 1) THEN
      WRITE(stdout, '(12x," nel       = ", E15.10)') nel
      WRITE(stdout, '(12x," meff      = ", E15.10)') meff
      WRITE(stdout, '(12x," rs        = ", E15.10)') rs
      WRITE(stdout, '(12x," kF        = ", E15.10)') kF
      WRITE(stdout, '(12x," vF        = ", E15.10)') vF
      WRITE(stdout, '(12x," fermi_en  = ", E15.10)') fermiHEG
      WRITE(stdout, '(12x," qTF       = ", E15.10)') qTF
      WRITE(stdout, '(12x," wpl       = ", E15.10)') wpl0
      WRITE(stdout, '(12x," qcut      = ", E15.10)') qcut
      WRITE(stdout, '(12x," eps0      = ", E15.10)') eps0
      WRITE(stdout, '(12x," epsiHEG   = ", E15.10)') epsiHEG
      WRITE(stdout, '(12x," deltaeps  = ", E15.10)') deltaeps
    ENDIF
    !
    IF (SQRT(qsquared) < qcut) THEN
      !
      ! loop over all k points of the fine mesh
      !
      fermicount = 0 
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle (but in this case they are the same)
        !
        IF ((MINVAL(ABS(etf (:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL( ABS(etf (:, ikq) - ef)) < fsthick)) THEN
          !
          fermicount = fermicount + 1
          !
          ! Bose occupation
          wgq = wgauss(-wq / eptemp, -99)
          wgq = wgq / (one - two * wgq)
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k (relative to Ef)
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !  
            DO jbnd = 1, nbndfst
              !
              !  the fermi occupation for k+q
              ekk1 = etf(ibndmin - 1 + jbnd, ikk) - ef0
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              wgkq = wgauss(-ekq / eptemp, -99)  
              !
              ! Computation of the dipole
              IF (ibnd == jbnd) THEN
                IF (SQRT(qsquared) > eps6) THEN
                  dipole = 1.0d0 / (qsquared * tpiba_new * tpiba_new)
                ELSE
                  dipole = 0.d0 
                ENDIF
              ELSE
                IF (ABS(ekq - ekk1) > eps6) THEN
                  dipole = REAL(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) * &
                               CONJG(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk)) / ((ekk1 - ekk)**2 + degaussw**2))
                ELSE 
                  dipole = 0.d0
                ENDIF
              ENDIF
              !
              g2 = dipole * 4.d0 * pi * (wq * deltaeps / 2.d0) / omega * 2.d0 ! The q^-2 is cancelled by the q->0 limit of the dipole. See e.g., pg. 258 of Grosso Parravicini. 
              !
              DO iw = 1, nw_specfun
                !
                ww = wmin_specfun + DBLE(iw - 1) * dw
                !
                weight = wqf(iq) * REAL(((wgkq + wgq) / (ww - (ekq - wq) - ci * degaussw)  +  &
                         (one - wgkq + wgq) / (ww - (ekq + wq) - ci * degaussw)))
                !
                esigmar_all(ibnd, ik + lower_bnd - 1, iw) = esigmar_all(ibnd, ik + lower_bnd - 1, iw) + g2 * weight 
                ! 
                ! SP : Application of the sum rule
                esigmar0 =  g2 *  wqf(iq) * REAL(((wgkq + wgq) / (-(ekq - wq) - ci * degaussw)  +  &
                    (one - wgkq + wgq) / (-(ekq + wq) - ci * degaussw)))
                esigmar_all(ibnd, ik + lower_bnd - 1, iw) = esigmar_all(ibnd, ik + lower_bnd - 1, iw) - esigmar0
                !
                weight = wqf(iq) * AIMAG((wgkq + wgq  / (ww - (ekq - wq) - ci * degaussw)  +  &
                         (one - wgkq + wgq) / (ww - (ekq + wq) - ci * degaussw)))
                !
                esigmai_all(ibnd, ik + lower_bnd - 1, iw) = esigmai_all(ibnd, ik + lower_bnd - 1, iw) + g2 * weight
                !
              ENDDO
            ENDDO !jbnd
          ENDDO !ibnd
        ENDIF ! endif  fsthick
      ENDDO ! end loop on k
    ENDIF
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! Note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      !
      ! Output electron spectral function here after looping over all q-points (with their contributions summed in a etc.)
      !
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      ! construct the trace of the spectral function (assume diagonal selfenergy
      ! and constant matrix elements for dipole transitions)
      !
      IF (me_pool == 0) then
        OPEN(UNIT = iospectral, FILE = 'specfun.plself') 
        OPEN(UNIT = iospectral_sup, FILE = 'specfun_sup.plself') 
      ENDIF
      IF (me_pool == 0) then
        WRITE(iospectral, '(/2x,a/)') '#Electron-plasmon spectral function (meV)'
        WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electron-plasmon self-energy (meV)' 
      ENDIF
      IF (me_pool == 0) then
        WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
        WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]       Real Sigma[meV]  Im Sigma[meV]'
      ENDIF
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout, '(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
        WRITE(stdout, '(5x,a)') REPEAT('-',67)
        !
        DO iw = 1, nw_specfun
          !
          ww = wmin_specfun + DBLE(iw - 1) * dw
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            a_all(iw, ik) = a_all(iw, ik) + ABS(esigmai_all(ibnd, ik, iw) ) / pi / &
                 ((ww - ekk - esigmar_all(ibnd, ik, iw))**two + (esigmai_all(ibnd, ik, iw))**two)
            !
          ENDDO
          !
          WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw, ik) / ryd2mev
          !
        ENDDO
        !
        WRITE(stdout, '(5x,a/)') REPEAT('-',67)
        !
      ENDDO
      !
      DO ik = 1, nktotf
        !
        ! The spectral function should integrate to 1 for each k-point
        specfun_sum = 0.0
        ! 
        DO iw = 1, nw_specfun
          !
          ww = wmin_specfun + DBLE(iw-1) * dw
          fermi(iw) = wgauss(-ww / eptemp, -99) 
          specfun_sum = specfun_sum + a_all(iw, ik) * fermi(iw) * dw !/ ryd2mev
          !
         IF (me_pool == 0) WRITE(iospectral, '(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw, ik) / ryd2mev
        ENDDO
        !
        IF (me_pool == 0) WRITE(iospectral, '(a)') ' '
        IF (me_pool == 0) WRITE(iospectral, '(2x,a,2x,e12.5)') '# Integrated spectral function ', specfun_sum
      ENDDO
      !
      IF (me_pool == 0) CLOSE(iospectral)
      !
      DO ibnd = 1, nbndfst
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          !  the energy of the electron at k
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
          !
          DO iw = 1, nw_specfun
            !
            ww = wmin_specfun + DBLE(iw - 1) * dw
            WRITE(stdout, '(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik, &
              ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd, ik, iw), &
              ryd2mev * esigmai_all(ibnd, ik, iw)
            ! 
            IF (me_pool == 0) &
            WRITE(iospectral_sup, '(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik, &
              ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd, ik, iw), &
              ryd2mev * esigmai_all(ibnd, ik, iw)
            !
          ENDDO
          !
        ENDDO
        !
        WRITE(stdout,*) ' '
        !
      ENDDO
      !
      IF (me_pool == 0) CLOSE(iospectral_sup)
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating etf_all', 1)
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    103 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', e12.5, ' meV^-1')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_pl_q
    !-----------------------------------------------------------------------
    ! 
  !-----------------------------------------------------------------------
  END MODULE spectral_func
  !-----------------------------------------------------------------------
