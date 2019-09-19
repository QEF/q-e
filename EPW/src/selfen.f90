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
  MODULE selfen
  !----------------------------------------------------------------------
  !! 
  !! This module contains the various self-energy routines
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_elec_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !! 
    !!  Compute the imaginary part of the electron self energy due to electron-
    !!  phonon interaction in the Migdal approximation. This corresponds to 
    !!  the electron linewidth (half width). The phonon frequency is taken into
    !!  account in the energy selection rule.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !!  This routines computes the contribution from phonon iq to all k-points
    !!  The outer loop in ephwann_shuffle.f90 will loop over all iq points
    !!  The contribution from each iq is summed at the end of this SUBROUTINE for iqq=totq
    !!  to recover the per-ik electron self energy
    !!
    !!  RM 24/02/2014
    !!  Redefined the size of sigmar_all, sigmai_all, and zi_all within the fermi windwow
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : linewidth_elself
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, shortrange, &
                              fsthick, eptemp, ngaussw, degaussw, &
                              eps_acustic, efermi_read, fermi_energy,&
                              restart, restart_freq
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, eta, nbndfst, &
                              nkf, epf17, wf, wqf, xkf, nkqtotf, adapt_smearing, &
                              sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew, &
                              nktotf, lower_bnd
    USE control_flags, ONLY : iverbosity
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_epw, ONLY : electron_write
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart 
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points from the selecq.fmt grid. 
    !
    ! Local variables 
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
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
    INTEGER :: ierr 
    !! Error status
    REAL(KIND = DP) :: tmp
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp2
    !! Temporary variable to store imag part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp3
    !! Temporary variable to store Z for the degenerate average
    REAL(KIND = DP) :: ekk2
    !! Temporary variable to the eigenenergies for the degenerate average
    REAL(KIND = DP) :: sigmar_tmp(nbndfst)
    !! Temporary array to store the real-part of Sigma 
    REAL(KIND = DP) :: sigmai_tmp(nbndfst)
    !! Temporary array to store the imag-part of Sigma 
    REAL(KIND = DP) :: zi_tmp(nbndfst)
    !! Temporary array to store the Z
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose occupation factor $n_{q\nu}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Self-energy factor 
    !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$ 
    !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$ 
    REAL(KIND = DP) :: w0g1
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: inv_eta(nbndfst, nmodes, nktotf)
    !! Inverse of the eta for speed purposes
    REAL(KIND = DP) :: eta2(nbndfst, nmodes, nktotf)
    !! Current eta 
    REAL(KIND = DP) :: inv_eta_tmp
    !! Temporary inv_eta
    REAL(KIND = DP) :: eta_tmp
    !! Temporrary eta2
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !  
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_eptemp = 1.0 / eptemp
    ! To avoid if branching in the loop
    inv_eta(:, :, :) = zero
    IF (adapt_smearing) THEN
      DO ik = 1, nkf
        DO ibnd = 1, nbndfst
          DO imode = 1, nmodes
            inv_eta(ibnd, imode, ik) = 1.0d0 / (SQRT(2d0) * eta(imode, ibnd, ik))
            eta2(ibnd, imode, ik) = SQRT(2d0) * eta(imode, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nkf
        DO ibnd = 1, nbndfst
          DO imode = 1, nmodes
            inv_eta(ibnd, imode, ik) = 1.0d0 / degaussw
            eta2(ibnd, imode, ik) = degaussw
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Now pre-treat phonon modes for efficiency
    ! Treat phonon frequency and Bose occupation
    wq(:) = zero
    DO imode = 1, nmodes
      wq(imode) = wf(imode, iq)
      IF (wq(imode) > eps_acustic) THEN
        g2_tmp(imode) = 1.0d0
        wgq(imode)    = wgauss(-wq(imode) * inv_eptemp, -99)
        wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
        inv_wq(imode) =  1.0d0 / (two * wq(imode))
      ELSE
        g2_tmp(imode) = 0.0
        wgq(imode)    = 0.0
        inv_wq(imode) = 0.0
      ENDIF
    ENDDO
    !
    IF (iqq == 1) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)') 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
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
      ef0 = efnew
      !
    ENDIF
    !
    IF ((iqq == 1) .AND. .NOT. adapt_smearing) THEN 
      WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
      WRITE (stdout, '(a)') ' '
    ENDIF
    !
    IF (restart) THEN
      ! Make everythin 0 except the range of k-points we are working on
      sigmar_all(:, 1:lower_bnd-1) = zero
      sigmar_all(:, lower_bnd + nkf:nktotf) = zero
      sigmai_all(:, 1:lower_bnd-1) = zero
      sigmai_all(:, lower_bnd + nkf:nktotf) = zero
      zi_all(:, 1:lower_bnd-1) = zero
      zi_all(:, lower_bnd + nkf:nktotf) = zero
      ! 
    ENDIF
    !
    ! In the case of a restart do not add the first step
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ELSE
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
            DO ibnd = 1, nbndfst
              !
              !  the energy of the electron at k (relative to Ef)
              ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
              ! 
              eta_tmp     = eta2(ibnd, imode, ik) 
              inv_eta_tmp = inv_eta(ibnd, imode, ik)  
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
                ! SP: Shortrange is disabled for efficiency reasons 
                !IF (shortrange .AND. ( ABS(xqf (1, iq))> eps8 .OR. ABS(xqf (2, iq))> eps8 &
                !   .OR. ABS(xqf (3, iq))> eps8 )) THEN                         
                !  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                !  !     number, in which case its square will be a negative number. 
                !  g2 = REAL((epf17 (jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode))
                !ELSE
                !  g2 = (ABS(epf17 (jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                !ENDIF        
                ! 
                g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                !
                ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                ! (Many-Particle Physics, 3rd edition)
                ! 
                weight = wqf(iq) * REAL(((wgkq + wgq(imode)) / (ekk - (ekq - wq(imode)) - ci * eta_tmp)  +  &
                                        (one - wgkq + wgq(imode)) / (ekk - (ekq + wq(imode)) - ci * eta_tmp)))
                !
                sigmar_all(ibnd, ik + lower_bnd - 1) = sigmar_all(ibnd, ik + lower_bnd - 1) + g2 * weight
                !
                ! Logical implementation
                ! weight = wqf(iq) * aimag (                                                  &
                !         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                !           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
                !
                ! Delta implementation 
                w0g1   = w0gauss((ekk - ekq + wq(imode)) * inv_eta_tmp, 0) * inv_eta_tmp
                w0g2   = w0gauss((ekk - ekq - wq(imode)) * inv_eta_tmp, 0) * inv_eta_tmp
                weight = pi * wqf(iq) * ((wgkq + wgq(imode)) * w0g1 + (one - wgkq + wgq(imode)) * w0g2)
                !
                sigmai_all(ibnd, ik + lower_bnd - 1) = sigmai_all(ibnd, ik + lower_bnd - 1) + g2 * weight
                !
                ! Mode-resolved
                IF (iverbosity == 3) THEN
                  sigmai_mode(ibnd, imode, ik + lower_bnd - 1) = sigmai_mode(ibnd, imode, ik + lower_bnd - 1) + g2 * weight
                ENDIF
                !
                ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                !
                weight = wqf(iq) * &
                         ((      wgkq + wgq(imode)) * ((ekk - (ekq - wq(imode)))**two - eta_tmp**two) /       &
                                                      ((ekk - (ekq - wq(imode)))**two + eta_tmp**two)**two +  &
                          (one - wgkq + wgq(imode)) * ((ekk - (ekq + wq(imode)))**two - eta_tmp**two) /       &
                                                      ((ekk - (ekq + wq(imode)))**two + eta_tmp**two)**two )  
                !
                zi_all(ibnd, ik + lower_bnd - 1) = zi_all(ibnd, ik + lower_bnd - 1) + g2 * weight
                ! 
              ENDDO !jbnd
            ENDDO !ibnd
          ENDDO !imode
        ENDIF ! endif  fsthick
      ENDDO ! end loop on k
      !
      ! Creation of a restart point
      IF (restart) THEN
        IF (MOD(iqq, restart_freq) == 0) THEN
          WRITE(stdout, '(a,i10)' ) '     Creation of a restart point at ', iqq
          CALL mp_sum(sigmar_all, inter_pool_comm)
          CALL mp_sum(sigmai_all, inter_pool_comm)
          CALL mp_sum(zi_all, inter_pool_comm)
          CALL mp_sum(fermicount, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          CALL electron_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
        ENDIF
      ENDIF 
    ENDIF ! in case of restart, do not do the first one
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(sigmar_all, inter_pool_comm)
      CALL mp_sum(sigmai_all, inter_pool_comm)
      IF (iverbosity == 3) CALL mp_sum(sigmai_mode, inter_pool_comm)
      CALL mp_sum(zi_all, inter_pool_comm)
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
      ! Average over degenerate eigenstates:
      WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
      ! 
      DO ik = 1, nktotf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk)
          n = 0
          tmp = 0.0_DP
          tmp2 = 0.0_DP
          tmp3 = 0.0_DP
          DO jbnd = 1, nbndfst
            ekk2 = etf_all(ibndmin - 1 + jbnd, ikk) 
            IF (ABS(ekk2-ekk) < eps6) THEN
              n    = n + 1
              tmp  =  tmp + sigmar_all(jbnd, ik)
              tmp2 = tmp2 + sigmai_all(jbnd, ik)
              tmp3 = tmp3 + zi_all(jbnd, ik)
            ENDIF
            !
          ENDDO ! jbnd
          sigmar_tmp(ibnd) = tmp / FLOAT(n)
          sigmai_tmp(ibnd) = tmp2 / FLOAT(n)
          zi_tmp(ibnd) = tmp3 / FLOAT(n)
          !
        ENDDO ! ibnd
        sigmar_all(:, ik) = sigmar_tmp(:) 
        sigmai_all(:, ik) = sigmai_tmp(:)
        zi_all(:, ik)     = zi_tmp(:)
        ! 
      ENDDO ! nktotf
      !  
      ! Output electron SE here after looping over all q-points (with their contributions 
      ! summed in sigmar_all, etc.)
      !
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      IF (mpime == ionode_id) THEN
        ! Write to file
        OPEN(UNIT = linewidth_elself, FILE = 'linewidth.elself')
        WRITE(linewidth_elself, '(a)') '# Electron linewidth = 2*Im(Sigma) (meV)'
        IF (iverbosity == 3) THEN
          WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      imode          Im(Sigma)(meV)'
        ELSE
          WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sigma)(meV)'
        ENDIF
        ! 
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout, '(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
          WRITE(stdout, '(5x,a)') REPEAT('-',67)
          !
          DO ibnd = 1, nbndfst
            !
            ! note that ekk does not depend on q 
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
            zi_all(ibnd, ik) = one / (one + zi_all(ibnd, ik))
            !
            WRITE(stdout, 102) ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik), &
                               ryd2mev * sigmai_all(ibnd,ik), zi_all(ibnd, ik), one / zi_all(ibnd, ik) - one
            IF (iverbosity == 3) THEN
              DO imode = 1, nmodes
                WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ik
                WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
                WRITE(linewidth_elself, '(E22.14,2x)', ADVANCE = 'no') ryd2ev * ekk
                WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') imode
                WRITE(linewidth_elself, '(E22.14,2x)') ryd2mev * sigmai_mode(ibnd, imode, ik)
              ENDDO
            ELSE
              WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ik
              WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
              WRITE(linewidth_elself, '(E22.14,2x)', ADVANCE = 'no') ryd2ev * ekk
              WRITE(linewidth_elself, '(E22.14,2x)') ryd2mev * sigmai_all(ibnd, ik)
            ENDIF
            !
          ENDDO
          WRITE(stdout, '(5x,a/)') REPEAT('-',67)
          !
        ENDDO
      ENDIF
      !
      DO ibnd = 1, nbndfst
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          ! note that ekk does not depend on q 
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
          !
          ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
          !zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
          !
          WRITE(stdout, '(2i9,5f12.4)') ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik), &
                         ryd2mev * sigmai_all(ibnd, ik), zi_all(ibnd, ik), one / zi_all(ibnd, ik) - one
          ! 
        ENDDO
        !
        WRITE(stdout, '(a)') '  '
        !
      ENDDO
      !
      CLOSE(linewidth_elself)
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error deallocating etf_all', 1)
      !
    ENDIF 
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    102 FORMAT(5x, 'E( ', i3, ' )=', f9.4, ' eV   Re[Sigma]=', f15.6, ' meV Im[Sigma]=', &
               f15.6, ' meV     Z=', f15.6, ' lam=', f15.6)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_elec_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_phon_q(iqq, iq, totq)
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
    !! the selfenergy for one phonon at a time.  Much smaller footprint on
    !! the disk
    !!
    !! RM 24/02/2014
    !! redefined the size of coskkq, vkk, vkq within the fermi windwow
    !! cleaned up the subroutine
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    use phcom,      ONLY : nmodes
    use epwcom,     ONLY : nbndsub, fsthick, efermi_read, fermi_energy,  &
                           eptemp, ngaussw, degaussw, shortrange,        &
                           nsmear, delta_smear, eps_acustic, specfun_ph, &
                           delta_approx, vme
    use pwcom,      ONLY : nelec, ef
    USE klist_epw,  ONLY : isk_dummy
    use elph2,      ONLY : epf17, ibndmax, ibndmin, etf, wkf, xqf, wqf, nkqf,   &
                           nkf, wf, nkqtotf, xqf, lambda_all, lambda_v_all,     &
                           dmef, vmef, gamma_all, gamma_v_all, efnew, nbndfst, &
                           nktotf, adapt_smearing
    USE mp,         ONLY : mp_barrier, mp_sum
    USE mp_global,  ONLY : inter_pool_comm
    USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi, eps4, eps6, eps8
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from the selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in selecq.fmt
    ! 
    ! Local variables 
    !
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
    INTEGER :: jmode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: ismear
    !! Upper bounds index after k or q paral
    !! Smearing for the Gaussian function 
    INTEGER :: n
    !! Counter on number of mode degeneracies
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: wq_tmp
    !! Temporary Phonon frequency on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Imaginary part of the phonhon self-energy factor 
    !!$$ \pi N_q \Im(\frac{f_{nk}(T) - f_{mk+q(T)}}{\varepsilon_{nk}-\varepsilon_{mk+q}-\omega_{q\nu}+i\delta}) $$
    !! In practice the imaginary is performed with a delta Dirac
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: w0g1
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: inv_wq
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: g2_tmp
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: gamma(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy 
    REAL(KIND = DP) :: gamma_v(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy multiplied by (1-coskkq)
    REAL(KIND = DP) :: coskkq(nbndfst, nbndfst)
    !! $$(v_k \cdot v_{k+q}) / |v_k|^2$$
    REAL(KIND = DP) :: DDOT
    !! Dot product function
    REAL(KIND = DP) :: degaussw0
    !! degaussw0 = (ismear-1) * delta_smear + degaussw
    REAL(KIND = DP) :: inv_degaussw0
    !! Inverse degaussw0 for efficiency reasons
    REAL(KIND = DP) :: lambda_tot
    !! Integrated lambda function
    REAL(KIND = DP) :: lambda_tr_tot
    !! Integrated transport lambda function
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: eptemp0
    !!eptemp0   = (ismear-1) * delta_smear + eptem
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electronic velocity $v_{nk}$
    REAL(KIND = DP) :: vkq(3, nbndfst)
    !! Electronic velocity $v_{nk+q}$
    REAL(KIND = DP) :: tmp
    !! Temporary value of lambda for av.
    REAL(KIND = DP) :: tmp2
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: tmp3
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: tmp4
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: lambda_tmp(nmodes)
    !! Temporary value of lambda for av.  
    REAL(KIND = DP) :: lambda_v_tmp(nmodes)
    !! Temporary value of lambda v for av.  
    REAL(KIND = DP) :: gamma_tmp(nmodes)
    !! Temporary value of gamma for av.  
    REAL(KIND = DP) :: gamma_v_tmp(nmodes)
    !! Temporary value of gamma v for av.  
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    !  
    IF (adapt_smearing) CALL errore('selfen_phon_q', 'adapt_smearing cannot be used with phonon self-energy ', 1) 
    ! 
    IF (iq == 1) THEN 
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x,a,f10.6,a)' ) &
           'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Golden Rule strictly enforced with T = ', eptemp * ryd2ev, ' eV'
      !
    ENDIF
    !
    DO ismear = 1, nsmear
      !
      degaussw0 = (ismear - 1) * delta_smear + degaussw
      eptemp0   = (ismear - 1) * delta_smear + eptemp
      ! 
      ! SP: Multiplication is faster than division ==> Important if called a lot
      !     in inner loops
      inv_degaussw0 = 1.0 / degaussw0
      inv_eptemp0   = 1.0 / eptemp0
      !
      ! Fermi level and corresponding DOS
      !
      IF (efermi_read) THEN
        !
        ef0 = fermi_energy
        !
      ELSE IF (nsmear > 1) THEN
        !
        ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk_dummy)
        ! if some bands are skipped (nbndskip /= 0), nelec has already been
        ! recalculated 
        ! in ephwann_shuffle
        !
      ELSE !SP: This is added for efficiency reason because the efermig routine is slow
        ef0 = efnew
      ENDIF
      !
      dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
      !  N(Ef) in the equation for lambda is the DOS per spin
      dosef = dosef / two
      !
      IF (iq == 1) THEN 
        WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
        WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
      ENDIF
      !
      CALL start_clock('PH SELF-ENERGY')
      !
      fermicount = 0
      wgkk = 0.0_DP
      w0g1 = 0.0_DP
      gamma(:)   = zero
      gamma_v(:) = zero
      !
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        coskkq = zero
        ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
        ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
        ! In principle the only coskkq contributing to lambda_tr are both near the
        ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
        ! we may implement the approximation to the angle between k and k+q
        ! vectors also listed in Grimvall
        !
        IF (vme) THEN 
          DO ibnd = 1, nbndfst
            DO jbnd = 1, nbndfst
              !
              ! vmef is in units of Ryd * bohr
              !
              vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
              vkq(:, jbnd) = REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
              IF (ABS(vkk(1, ibnd)**2 + vkk(2, ibnd)**2 + vkk(3, ibnd)**2) > eps4) &
                  coskkq(ibnd, jbnd ) = DDOT(3, vkk(:, ibnd ), 1, vkq(:, jbnd), 1) / &
                  DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd), 1)
            ENDDO
          ENDDO
        ELSE
          DO ibnd = 1, nbndfst
            DO jbnd = 1, nbndfst
              !
              ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
              ! 1/m  = 2 in Rydberg atomic units
              !
              vkk(:, ibnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
              vkq(:, jbnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
              IF (ABS(vkk(1, ibnd)**2 + vkk(2, ibnd)**2 + vkk(3, ibnd)**2) > eps4) &
                  coskkq(ibnd, jbnd) = DDOT(3, vkk(:, ibnd), 1, vkq(:, jbnd),1)  / &
                  DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd),1)
            ENDDO
          ENDDO
        ENDIF
        !
        ! Here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
          !
          fermicount = fermicount + 1
          DO imode = 1, nmodes
            !
            ! the phonon frequency
            wq = wf(imode, iq)
            !
            ! SP : We should avoid branching statements (if statements) in
            !      innerloops. Therefore we do it here.
            inv_wq = 1.0d0 / (two * wq)
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
              IF (delta_approx) THEN
                w0g1 = w0gauss(ekk / degaussw0, 0) / degaussw0
              ELSE
                wgkk = wgauss(-ekk * inv_eptemp0, -99)
              ENDIF
              !
              DO jbnd = 1, nbndfst
                !
                !  the fermi occupation for k+q
                ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                !
                ! here we take into account the zero-point SQRT(hbar/2M\omega)
                ! with hbar = 1 and M already contained in the eigenmodes
                ! g2 is Ry^2, wkf must already account for the spin factor
                !
                IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                   .OR. ABS(xqf(3, iq)) > eps8)) THEN              
                  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                  !     number, in which case its square will be a negative number. 
                  g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq * g2_tmp) 
                ELSE
                  g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                ENDIF
                !
                IF (delta_approx) THEN 
                  !
                  w0g2 = w0gauss(ekq / degaussw0, 0) / degaussw0
                  ! the expression below is positive-definite, but also an
                  ! approximation which neglects some fine features
                  weight = pi * wq * wkf(ikk) * w0g1 * w0g2
                  !
                ELSE
                  !
                  wgkq = wgauss(-ekq * inv_eptemp0, -99)
                  !
                  ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q + id]
                  ! This is the imaginary part of minus the phonon self-energy, sans
                  ! the matrix elements
                  !
                  !weight = wkf (ikk) * (wgkk - wgkq) * &
                  !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) )
                  !
                  ! SP: The expression below (minus phonon self-energy) corresponds to
                  !  = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q]
                  weight = pi * wkf(ikk) * (wgkk - wgkq) * w0gauss((ekq - ekk - wq) / degaussw0, 0) / degaussw0
                  !
                ENDIF  
                !
                gamma(imode) = gamma(imode) + weight * g2 
                gamma_v(imode) = gamma_v(imode) + weight * g2 * (1.0d0 - coskkq(ibnd, jbnd)) 
                !
              ENDDO ! jbnd
            ENDDO   ! ibnd
          ENDDO ! loop on q-modes
        ENDIF ! endif fsthick
      ENDDO ! loop on k
      !
      CALL stop_clock('PH SELF-ENERGY')
      !
      ! collect contributions from all pools (sum over k-points)
      ! this finishes the integral over the BZ  (k)
      !
      CALL mp_sum(gamma, inter_pool_comm) 
      CALL mp_sum(gamma_v, inter_pool_comm) 
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      ! 
      ! An average over degenerate phonon-mode is performed. 
      DO imode = 1, nmodes
        n = 0
        tmp = 0.0_DP
        tmp2 = 0.0_DP
        tmp3 = 0.0_DP
        tmp4 = 0.0_DP
        wq = wf(imode, iq)
        DO jmode = 1, nmodes
          wq_tmp = wf(jmode, iq)
          IF (ABS(wq - wq_tmp) < eps6) THEN
            n = n + 1
            IF (wq_tmp > eps_acustic) THEN 
              tmp  =  tmp  + gamma(jmode) / pi / wq**two / dosef
              tmp2 =  tmp2 + gamma_v(jmode) / pi / wq**two / dosef
            ENDIF
            tmp3 =  tmp3 + gamma(jmode)
            tmp4 =  tmp4 + gamma_v(jmode)
          ENDIF
        ENDDO ! jbnd
        lambda_tmp(imode)   = tmp / FLOAT(n)
        lambda_v_tmp(imode) = tmp2 / FLOAT(n)
        gamma_tmp(imode)    = tmp3 / FLOAT(n)
        gamma_v_tmp(imode)  = tmp4 / FLOAT(n)
      ENDDO
      lambda_all(:, iq, ismear)   = lambda_tmp(:)
      lambda_v_all(:, iq, ismear) = lambda_v_tmp(:)
      gamma_all(:, iq, ismear)    = gamma_tmp(:)
      gamma_v_all(:, iq, ismear)  = gamma_v_tmp(:)
      lambda_tot = SUM(lambda_all(:, iq, ismear))
      lambda_tr_tot = SUM(lambda_v_all(:, iq, ismear))
      !
      WRITE(stdout, '(/5x,"ismear = ",i5," iq = ",i7," coord.: ", 3f9.5, " wt: ", f9.5)') ismear, iq, xqf(:, iq), wqf(iq)
      WRITE(stdout, '(5x,a)') REPEAT('-',67)
      !
      DO imode = 1, nmodes
        ! 
        wq = wf(imode, iq)
        WRITE(stdout, 102) imode, lambda_all(imode, iq, ismear), ryd2mev * gamma_all(imode, iq, ismear), ryd2mev * wq
        WRITE(stdout, 104) imode, lambda_v_all(imode, iq, ismear), ryd2mev * gamma_v_all(imode, iq, ismear), ryd2mev * wq
        !
      ENDDO
      !
      WRITE(stdout, 103) lambda_tot
      WRITE(stdout, 105) lambda_tr_tot
      ! 
      IF (.NOT. specfun_ph) THEN
        WRITE(stdout, '(5x,a/)') REPEAT('-',67)
        WRITE(stdout, '(/5x,a,i8,a,i8/)' ) 'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', nktotf
      ENDIF
      !
    ENDDO !smears
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    102 FORMAT(5x, 'lambda___( ', i3, ' )=', f15.6, '   gamma___=', f15.6, ' meV', '   omega=', f12.4, ' meV')
    103 FORMAT(5x, 'lambda___( tot )=', f15.6)
    104 FORMAT(5x, 'lambda_tr( ',i3,' )=', f15.6, '   gamma_tr=', f15.6, ' meV', '   omega=', f12.4, ' meV')
    105 FORMAT(5x, 'lambda_tr( tot )=', f15.6)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_phon_q
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_pl_q(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the imaginary part of the electron self energy due to electron-
    !!  plasmon interaction. 
    !! 
    !!  The coupling coefficients have been evaluated analytically employing a 
    !!  Lindhard function model for the dielectric function contribution due to 
    !!  the extrinsic carriers.
    !! 
    !!  There are 3 parameters that the users should provide in the input: 
    !!    - DOS effective mass; 
    !!    - carrier concentration (Only for doped semiconductors, it shouldnt be used for insulators);
    !!    - epsilon_infinity (e.g, from exp. or from RPA).
    !!
    !!  F. Caruso and S. Ponce - 2017 
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : linewidth_elself
    USE epwcom,        ONLY : nbndsub, fsthick, eptemp, ngaussw,  &
                              efermi_read, fermi_energy, degaussw,& 
                              nel, meff, epsiHEG 
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, dmef, adapt_smearing, &
                              nkf, wqf, xkf, nkqtotf, efnew, nbndfst, nktotf,  &
                              sigmar_all, sigmai_all, sigmai_mode, zi_all
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm 
    USE cell_base,     ONLY : omega, alat, bg
    USE division,      ONLY : fkbounds
    USE poolgathering, ONLY : poolgather2
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Q-index from the selected q
    INTEGER, INTENT(in) :: iq 
    !! Q-index from the global q
    INTEGER, INTENT(in) :: totq
    !! Number of q-points in selecq window
    ! 
    ! Local varialbes
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
    INTEGER :: n
    !! 
    INTEGER :: ierr
    !! Error status
    !! Integer for the degenerate average over eigenstates
    REAL(KIND = DP) :: tmp
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp2
    !! Temporary variable to store imag part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp3
    !! Temporary variable to store Z for the degenerate average
    REAL(KIND = DP) :: ekk2
    !! Temporary variable to the eigenenergies for the degenerate average
    REAL(KIND = DP) :: sigmar_tmp(nbndfst)
    !! Temporary array to store the real-part of Sigma 
    REAL(KIND = DP) :: sigmai_tmp(nbndfst)
    !! Temporary array to store the imag-part of Sigma 
    REAL(KIND = DP) :: zi_tmp(nbndfst)
    !! Temporary array to store the Z
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
    REAL(KIND = DP) :: weight
    !! 
    REAL(KIND = DP) :: w0g1
    !! 
    REAL(KIND = DP) :: w0g2
    !! 
    REAL(KIND = DP) :: inv_eptemp0
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
    !!
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !!
    !
    IF (adapt_smearing) CALL errore('selfen_pl_q', 'adapt_smearing cannot be used with plasmon self-energy ', 1) 
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_eptemp0  = 1.0 / eptemp
    !
    IF (iqq == 1) THEN
      !
      WRITE(stdout,'(/5x,a)') REPEAT('=',67)
      WRITE(stdout,'(5x,"Electron-plasmon Self-Energy in the Migdal Approximation")')
      WRITE(stdout,'(5x,a/)') REPEAT('=',67)
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
      ef0 = efnew ! Fermi energy is recalculated on the fine mesh!! 
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated  in ephwann_shuffle
      !
    ENDIF
    !
    IF (iqq == 1) THEN 
      WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
      WRITE (stdout,'(a)') ' '
    ENDIF
    !
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    !nel      =  0.01    ! this should be read from input - # of doping electrons 
    !epsiHEG  =  12.d0   ! this should be read from input - # dielectric constant at zero doping  
    !meff     =  0.25    ! this should be read from input - effective mass 
    !
    tpiba_new = 2.0d0 * pi / alat
    degen     = 1.0d0
    rs        = (3.d0 / (4.d0 * pi * nel / omega / degen))**(1.d0 / 3.d0) * meff * degen ! omega is the unit cell volume in Bohr^3
    kF        = (3.d0 * (pi**2) * nel / omega / degen)**(1.d0 / 3.d0) 
    vF        = (1.d0 / meff) * (3.d0 * (pi**2) * nel / omega / degen)**(1.d0 / 3.d0) 
    fermiHEG  = (1.d0 / (2.d0 * meff)) * (3.d0 * (pi**2) * nel / omega / degen)**(2.d0 / 3.d0) * 2.d0 ! [Ryd] multiplication by 2 converts from Ha to Ry
    qTF       = (6.d0 * pi * nel / omega / degen / (fermiHEG / 2.d0))**(1.d0 / 2.d0)    ! [a.u.]
    wpl0      = SQRT(4.d0 * pi * nel / omega / meff / epsiHEG) * 2.d0 ! [Ryd] multiplication by 2 converts from Ha to Ryd
    wq        = wpl0 ! [Ryd] 
    q(:)      = xqf(:, iq)
    CALL cryst_to_cart(1, q, bg, 1)
    qsquared  = (q(1)**2 + q(2)**2 + q(3)**2)
    qin       = SQRT(qsquared) * tpiba_new
    qcut      = wpl0 / vF / tpiba_new / 2.d0    ! 1/2 converts from Ryd to Ha
    !
    !if (.TRUE.) qcut = qcut / 2.d0 ! renorm to account for Landau damping 
    !
    CALL get_eps_mahan(qin, rs, kF, eps0) ! qin should be in atomic units for Mahan formula
    deltaeps = -(1.d0 / (epsiHEG + eps0 - 1.d0) - 1.d0 / epsiHEG)
    !
    IF (iqq == 1) THEN 
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
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        ! (but in this case they are the same)
        !
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
          !
          fermicount = fermicount + 1
          wgq = wgauss(-wq * inv_eptemp0, -99)
          wgq = wgq / (one - two * wgq)
          !
          DO ibnd = 1, nbndfst
            !
            !  the energy of the electron at k (relative to Ef)
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !
            DO jbnd = 1, nbndfst
              !
              ekk1 = etf(ibndmin - 1 + jbnd, ikk) - ef0
              ekq  = etf(ibndmin - 1 + jbnd, ikq) - ef0
              wgkq = wgauss(-ekq * inv_eptemp0, -99)  
              !
              ! Computation of the dipole
              IF (ibnd == jbnd) THEN
                IF (SQRT(qsquared) > eps8) THEN
                  dipole = 1.0d0 / (qsquared * tpiba_new * tpiba_new)
                ELSE
                  dipole = 0.d0 
                ENDIF
              ELSE
                IF (ABS(ekk - ekk1) > eps8) THEN
                  dipole = REAL(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) * &
                               CONJG(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk)) / ((ekk1 - ekk)**2 + degaussw**2)) 
                ELSE 
                  dipole = 0.d0
                ENDIF
              ENDIF
              !
              IF (ABS(dipole * (qsquared * tpiba_new * tpiba_new)) > 1) THEN
                dipole = 1.0d0 / (qsquared * tpiba_new * tpiba_new)
              ENDIF
              !
              g2 = dipole * 4.d0 * pi * (wq * deltaeps / 2.d0) / omega * 2.d0 ! The q^-2 is cancelled by the q->0 limit of the dipole. See e.g., pg. 258 of Grosso Parravicini. 
              !
              ! The fermi occupation for k+q
              !
              ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
              ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
              ! (Many-Particle Physics, 3rd edition)
              ! 
              weight = wqf(iq) * REAL(((wgkq + wgq) / (ekk - (ekq - wq) - ci * degaussw)  +  &
                        (one - wgkq + wgq) / (ekk - (ekq + wq) - ci * degaussw)))
              !
              sigmar_all(ibnd, ik + lower_bnd - 1) = sigmar_all(ibnd, ik + lower_bnd - 1) + g2 * weight
              !
              ! Delta implementation 
              w0g1 = w0gauss((ekk - ekq + wq) / degaussw, 0) / degaussw
              w0g2 = w0gauss((ekk - ekq - wq) / degaussw, 0) / degaussw
              weight = pi * wqf(iq) * ((wgkq + wgq) * w0g1 + (one - wgkq + wgq) * w0g2)
              !
              sigmai_all(ibnd, ik + lower_bnd - 1) = sigmai_all(ibnd, ik + lower_bnd - 1) + g2 * weight
              !
              ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
              !
              weight = wqf(iq) * &
                      ((       wgkq + wgq) * ((ekk - (ekq - wq))**two - degaussw**two) /       &
                                             ((ekk - (ekq - wq))**two + degaussw**two)**two +  &
                        (one - wgkq + wgq) * ((ekk - (ekq + wq))**two - degaussw**two) /       &
                                             ((ekk - (ekq + wq))**two + degaussw**two)**two )  
              !
              zi_all(ibnd, ik + lower_bnd - 1) = zi_all(ibnd, ik + lower_bnd - 1) + g2 * weight
              !
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
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(sigmar_all, inter_pool_comm)
      CALL mp_sum(sigmai_all, inter_pool_comm)
      CALL mp_sum(zi_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      !
      ! Average over degenerate eigenstates:
      WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
      ! 
      DO ik = 1, nktotf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk)
          n = 0
          tmp = 0.0_DP
          tmp2 = 0.0_DP
          tmp3 = 0.0_DP
          DO jbnd = 1, nbndfst
            ekk2 = etf_all(ibndmin - 1 + jbnd, ikk)
            IF (ABS(ekk2 - ekk) < eps6) THEN
              n = n + 1
              tmp =  tmp + sigmar_all(jbnd, ik)
              tmp2 =  tmp2 + sigmai_all(jbnd, ik)
              tmp3 =  tmp3 + zi_all(jbnd, ik)
            ENDIF
            ! 
          ENDDO ! jbnd
          sigmar_tmp(ibnd) = tmp / FLOAT(n)
          sigmai_tmp(ibnd) = tmp2 / FLOAT(n)
          zi_tmp(ibnd) = tmp3 / FLOAT(n)
          !
        ENDDO ! ibnd
        sigmar_all(:, ik) = sigmar_tmp(:)
        sigmai_all(:, ik) = sigmai_tmp(:)
        zi_all(:, ik) = zi_tmp(:)
        ! 
      ENDDO ! nktotf
      !
      ! Output electron SE here after looping over all q-points (with their contributions summed in sigmar_all, etc.)
      !
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      ! Write to file
      OPEN(UNIT = linewidth_elself, FILE = 'linewidth.plself')
      WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
      WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sgima)(meV)'
      ! 
      DO ik = 1, nktotf
         !
         ikk = 2 * ik - 1
         ikq = ikk + 1
         !
         WRITE(stdout, '(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
         WRITE(stdout, '(5x,a)') REPEAT('-',67)
         !
         DO ibnd = 1, nbndfst
           !
           ! note that ekk does not depend on q 
           ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
           !
           ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
           zi_all(ibnd, ik) = one / (one + zi_all(ibnd, ik))
           !
           WRITE(stdout, 102) ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik), &
                              ryd2mev * sigmai_all(ibnd, ik), zi_all(ibnd, ik), one / zi_all(ibnd, ik) - one
           WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ik
           WRITE(linewidth_elself, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
           WRITE(linewidth_elself, '(E22.14,2x)', ADVANCE = 'no') ryd2ev * ekk
           WRITE(linewidth_elself, '(E22.14,2x)') ryd2mev * sigmai_all(ibnd, ik)
           !
         ENDDO
         WRITE(stdout, '(5x,a/)') REPEAT('-',67)
      ENDDO
      !
      CLOSE(linewidth_elself)
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error deallocating etf_all', 1)
      !
    ENDIF 
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    102 FORMAT(5x, 'E( ', i3, ' )=', f9.4, ' eV   Re[Sigma]=', f15.6, ' meV Im[Sigma]=', &
               f15.6, ' meV     Z=', f15.6, ' lam=', f15.6)
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_pl_q
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE get_eps_mahan(q, rs, kF, eps0)
    !--------------------------------------------------------------------------
    !!
    !! Based on Eq. 5.166 of Mahan 2000. 
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi, eps6, eps10
    ! 
    IMPLICIT NONE
    ! 
    REAL(KIND = DP), INTENT(out) :: eps0
    !! Output dielectric function
    REAL(KIND = DP), INTENT(in) ::  q
    !! Q-point
    REAL(KIND = DP), INTENT(in) :: rs
    !! 
    REAL(KIND = DP), INTENT(in) :: kF
    !! Fermi wave-vector
    ! 
    !Local variable
    REAL(KIND = DP) :: x 
    !! 
    REAL(KIND = DP) :: alpha
    !!  
    ! 
    alpha = (4.d0 / (9.d0 * pi))**(1.d0/3.d0)
    ! 
    IF (ABS(q) > eps10) THEN
      x    = q / (2.d0 * kF) 
      eps0 = 1.d0 + (1.d0 - x**2) / (2.d0 * x) * LOG(ABS((1.d0 + x)/(1.d0 - x))) 
      eps0 = 1.d0 + alpha * rs * eps0 / (2.d0 * pi * (x**2))
    ELSE
      x    = (q + eps6) / 2.d0 / kF 
      eps0 = 1.d0 + (1.d0 - x**2) / (2.d0 * x) * LOG(ABS((1.d0 + x) / (1.d0 - x))) 
      eps0 = 1.d0 + alpha * rs / 2.d0 / pi / x**2 * eps0
    ENDIF
    ! 
    !--------------------------------------------------------------------------
    END SUBROUTINE get_eps_mahan
    !--------------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    FUNCTION dos_ef_seq(ngauss, degauss, ef, et, wk, nks, nbnd)
    !-----------------------------------------------------------------------
    !
    USE kinds, ONLY : DP
    USE mp,    ONLY : mp_sum
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: ngauss
    !! Number of smearing
    INTEGER, INTENT(in) :: nbnd
    !! Total number of bands considered
    INTEGER, INTENT(in) :: nks
    !!  Number of kpoints
    REAL(KIND = DP), INTENT(in) :: et(nbnd, nks)
    !! Eigenenergies
    REAL(KIND = DP), INTENT(in) :: wk(nks)
    !! K-point weights
    REAL(KIND = DP), INTENT(in) :: ef
    !! Fermi level  
    REAL(KIND = DP), INTENT(in) :: degauss
    !! Smearing value
    ! 
    REAL(KIND = DP) :: dos_ef_seq
    !! Output of the function
    !
    ! Local variables
    INTEGER :: ik
    !! K-point value
    INTEGER :: ibnd
    !! Band number
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! Fermi-Dirac function
    !
    ! Compute DOS at E_F (states per Ry per unit cell)
    !
    dos_ef_seq = 0.0d0
    DO ik = 1, nks
      DO ibnd = 1, nbnd
        dos_ef_seq = dos_ef_seq + wk(ik) * w0gauss((et(ibnd, ik) - ef) / degauss, ngauss) / degauss
      ENDDO
    ENDDO
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dos_ef_seq
    !-----------------------------------------------------------------------
    ! 
  !-----------------------------------------------------------------------
  END MODULE selfen
  !-----------------------------------------------------------------------
