  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE indabs
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines related to phonon assisted absorption
  !! 12/03/2018 Kyle and E. Kioupakis: First implementation
  !! 08/04/2018 S. Ponce: Cleaning
  !! 17/09/2019 S. Ponce: Modularization and cleaning
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE indabs_main(iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !! Main routine for phonon assisted absorption
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nstemp, fsthick, degaussw, &
                              eps_acustic, efermi_read, fermi_energy,&
                              vme, omegamin, omegamax, omegastep, indabs_fca, &
                              nomega, neta, restart, restart_step
    USE elph2,         ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf, &
                              sigmar_all, efnew, gtemp, &
                              dmef, omegap, epsilon2_abs, epsilon2_abs_lorenz, vmef, &
                              nbndfst, nktotf, ef0_fca
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, czero
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega
    USE io_indabs,     ONLY : indabs_write
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points from the selecq.fmt grid.
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp
    !! Temperature, in string format
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
    INTEGER :: nksqtotf
    !! Total number of k+q points
    INTEGER :: iw
    !! Index for frequency
!    INTEGER :: nomega
    !! Number of points on the photon energy axis
    INTEGER :: mbnd
    !! Index for summation over intermediate bands
    INTEGER :: ipol
    !! Polarization direction
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
!    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: wq(nmodes), nqv(nmodes)
    !! Phonon frequencies and phonon occupations on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgkk, wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
    REAL(KIND = DP) :: weighta, weighte
    !!- delta function for absorption, emission
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: pfac
    !! Occupation prefactors
    REAL(KIND = DP) :: pface
    !! Occupation prefactors
    REAL(KIND = DP) :: cfac
    !! Absorption prefactor
    REAL(KIND = DP) :: eta(9) = (/ 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 /) / ryd2eV
    !! Imaginary broadening of matrix element denominators
    REAL(KIND = DP) :: etemp_fca
    !! Temperature for fermi level calculation
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    COMPLEX(KIND = DP) :: vkk(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: vkq(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: s1a(3), s1e(3), s2a(3), s2e(3)
    !! Transition probability function
    COMPLEX(KIND = DP) :: epf(nbndfst, nbndfst, nmodes)
    !! Generalized matrix elements for phonon-assisted absorption
    !
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    !
    inv_degaussw = 1.0 / degaussw
    !
    ! Done outside
!    nomega = INT((omegamax - omegamin) / omegastep) + 1
    !
    ! 300 K
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 for spin (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    !
    IF (iq == 1) THEN
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Phonon-assisted absorption")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a)') 'The following temperatures are calculated:'
      DO itemp = 1, nstemp
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Temperature T = ', gtemp(itemp) * ryd2ev, ' eV'
      ENDDO
      !
      !IF (.NOT. ALLOCATED (omegap) )    ALLOCATE(omegap(nomega))
      !IF (.NOT. ALLOCATED (epsilon2_abs) ) ALLOCATE(epsilon2_abs(3, nomega, neta))
      !IF (.NOT. ALLOCATED (epsilon2_abs_lorenz) ) ALLOCATE(epsilon2_abs_lorenz(3, nomega, neta))
      ALLOCATE(omegap(nomega), STAT = ierr)
      IF (ierr /= 0) CALL errore('indabs', 'Error allocating omegap', 1)
      ! Now move the allocation into ephwann instead of indabs
!      ALLOCATE(epsilon2_abs(3, nomega, neta, nstemp), STAT = ierr)
!      IF (ierr /= 0) CALL errore('indabs', 'Error allocating epsilon2_abs', 1)
!      ALLOCATE(epsilon2_abs_lorenz(3, nomega, neta, nstemp), STAT = ierr)
!      IF (ierr /= 0) CALL errore('indabs', 'Error allocating epsilon2_abs_lorenz', 1)
      !
      epsilon2_abs = 0.d0
      epsilon2_abs_lorenz = 0.d0
      DO iw = 1, nomega
        omegap(iw) = omegamin + (iw - 1) * omegastep
      ENDDO
!      IF (indabs_fca) THEN
        ! Calculates free carrier fermi level
!        DO itemp = 1, nstemp
!          etemp_fca = gtemp(itemp)
!          CALL fermi_carrier_indabs(itemp, etemp_fca, ef0_fca)
!        ENDDO
!      ENDIF
    ENDIF
    !
    ! The total number of k points
    !
    nksqtotf = nktotf ! odd-even for k,k+q
    !
    DO itemp = 1, nstemp
      IF (first_cycle .and. itemp == nstemp) THEN
        first_cycle = .false.
      ELSE
        !        
        IF (indabs_fca) THEN
          !
          ef0 = ef0_fca(itemp)
        ELSEIF (efermi_read) THEN
          !
          ef0 = fermi_energy
        ELSE
          !
          ef0 = efnew
        ENDIF
        !
        inv_eptemp0 = 1.0 / gtemp(itemp)
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          DO imode = 1, nmodes
            !
            ! the phonon frequency at this q and nu
            wq(imode) = wf(imode, iq)
            !
            epf(:, :, imode) = epf17(:, :, imode,ik)
            IF (wq(imode) > eps_acustic) THEN
              nqv(imode) = wgauss(-wq(imode) / gtemp(itemp), -99)
              nqv(imode) = nqv(imode) / (one - two * nqv(imode))
            ENDIF
          ENDDO
          !
          ! RM - vme version should be checked
          IF (vme) THEN
            DO ibnd = 1, nbndfst
              DO jbnd = 1, nbndfst
                ! vmef is in units of Ryd * bohr
                vkk(:, ibnd, jbnd) = vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikk)
                vkq(:, ibnd, jbnd) = vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikq)
              ENDDO
            ENDDO
          ELSE
            DO ibnd = 1, nbndfst
              DO jbnd = 1, nbndfst
                ! Dme's already corrected for GW corrections in wan2bloch.f90
                vkk(:, ibnd, jbnd) = 2.0 * dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikk)
                vkq(:, ibnd, jbnd) = 2.0 * dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikq)
              ENDDO
            ENDDO
          ENDIF
          !
          DO ibnd = 1, nbndfst
            !  the energy of the electron at k (relative to Ef)
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !
            IF (ABS(ekk) < fsthick) THEN
              !
              wgkk = wgauss(-ekk * inv_eptemp0, -99)
              !
              DO jbnd = 1, nbndfst
                !
                ! The fermi occupation for k+q
                ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                !
                IF (ABS(ekq) < fsthick .AND. ekq < ekk + wq(nmodes) + omegamax + 6.0 * degaussw) THEN
                  !
                  wgkq = wgauss(-ekq * inv_eptemp0, -99)
                  !
                  IF (ekq - ekk - wq(nmodes) - omegamax > 6.0 * degaussw) CYCLE
                  IF (ekq - ekk + wq(nmodes) - omegamin < - 6.0 * degaussw) CYCLE
                  !
                  DO imode = 1, nmodes
                    !
                    IF (wq(imode) > eps_acustic) THEN
                      !
                      DO m = 1, neta
                        s1a = czero
                        s1e = czero
                        s2a = czero
                        s2e = czero
                        !
                        DO mbnd = 1, nbndfst
                          !
                          ! The energy of the electron at k (relative to Ef)
                          ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                          ! The energy of the electron at k+q (relative to Ef)
                          ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                          !
                          s1a(:) = s1a(:) + epf(mbnd, jbnd,imode) * vkk(:, ibnd, mbnd) / &
                                   (ekmk  - ekq + wq(imode) + ci * eta(m))
                          s1e(:) = s1e(:) + epf(mbnd, jbnd,imode) * vkk(:, ibnd, mbnd) / &
                                   (ekmk  - ekq - wq(imode) + ci * eta(m))
                          s2a(:) =  s2a(:) + epf(ibnd, mbnd,imode) * vkq(:, mbnd, jbnd) / &
                                   (ekmq  - ekk - wq(imode)+ ci * eta(m))
                          s2e(:) =  s2e(:) + epf(ibnd, mbnd,imode) * vkq(:, mbnd, jbnd) / &
                                   (ekmq  - ekk + wq(imode)+ ci * eta(m))
                        ENDDO
                        !
                        pfac  =  nqv(imode)      * wgkk * (one - wgkq) - (nqv(imode) + one) * (one - wgkk) * wgkq
                        pface = (nqv(imode) + one) * wgkk * (one - wgkq) -  nqv(imode) * (one - wgkk) * wgkq
                        !
                        DO iw = 1, nomega
                          !
                          IF (ABS(ekq - ekk - wq(imode) - omegap(iw)) > 6.0 * degaussw .AND. &
                              ABS(ekq - ekk + wq(imode) - omegap(iw)) > 6.0 * degaussw) CYCLE
                          !
                          weighte = w0gauss((ekq - ekk - omegap(iw) + wq(imode)) / degaussw, 0) / degaussw
                          weighta = w0gauss((ekq - ekk - omegap(iw) - wq(imode)) / degaussw, 0) / degaussw
                          !
                          DO ipol = 1, 3
                            epsilon2_abs(ipol, iw, m, itemp) = epsilon2_abs(ipol, iw, m, itemp) + (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pfac  * weighta * ABS(s1a(ipol) + s2a(ipol))**2 / (2 * wq(imode) * omega)
                            epsilon2_abs(ipol, iw, m, itemp) = epsilon2_abs(ipol, iw, m, itemp) + (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pface * weighte * ABS(s1e(ipol) + s2e(ipol))**2 / (2 * wq(imode) * omega)
                            epsilon2_abs_lorenz(ipol, iw, m, itemp) = epsilon2_abs_lorenz(ipol, iw, m, itemp) + &
                                  (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pfac  * ABS(s1a(ipol) + s2a(ipol))**2 / (2 * wq(imode) * omega) * &
                                 (degaussw / (degaussw**2 + (ekq - ekk - omegap(iw) - wq(imode))**2)) / pi
                            epsilon2_abs_lorenz(ipol, iw, m, itemp) = epsilon2_abs_lorenz(ipol, iw, m, itemp)  + &
                                 (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pface * ABS(s1e(ipol) + s2e(ipol))**2 / (2 * wq(imode) * omega) * &
                                 (degaussw / (degaussw**2 + (ekq - ekk - omegap(iw) + wq(imode))**2 )) / pi
                          ENDDO ! ipol
                        ENDDO ! iw
                      ENDDO ! neta
                    ENDIF ! if wq > acoustic
                  ENDDO ! imode
                ENDIF ! endif  ekq in fsthick
              ENDDO ! jbnd
            ENDIF  ! endif  ekk in fsthick
          ENDDO ! ibnd
        ENDDO ! ik
        IF (restart) THEN
          IF (MOD(iq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iq
            CALL mp_sum(epsilon2_abs, inter_pool_comm)
            CALL mp_sum(epsilon2_abs_lorenz, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL indabs_write(iq, totq, epsilon2_abs, epsilon2_abs_lorenz)
          ENDIF
        ENDIF
      ENDIF ! Skip first step in restart
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iq == nqtotf) THEN
      !
#if defined(__MPI)
      !
      ! Note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(epsilon2_abs, inter_pool_comm)
      CALL mp_sum(epsilon2_abs_lorenz, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      ! Output to stdout
      c = 'X'
      WRITE(c,"(i0)") neta
      format_string = "(5x,f15.6," // TRIM(c) // "E22.14)"
      !
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'Phonon-assisted absorption versus energy'
      WRITE(stdout, '(5x,a,4f15.6)') 'Broadenings: ', eta(1:4)
      WRITE(stdout, '(5x,a,5f15.6)') '   ', eta(5:9)
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'For the first Broadening and Temperature we have:'
      ! For test-farm checking purposes, only show m=1
      WRITE(stdout, '(5x,a)') 'Photon energy (eV), Imaginary dielectric function along x,y,z'
      DO iw = 1, nomega
        WRITE(stdout, '(5x,f15.6,3E22.14)') omegap(iw) * ryd2ev, (epsilon2_abs(ipol, iw, 1, 1), ipol = 1, 3)
      ENDDO
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'Values with other broadenings for temperature X are reported in the files epsilon2_indabs_X.dat'
      WRITE(stdout, '(5x,a)')
      !
      ! Output to file
      DO itemp = 1,nstemp
        WRITE(c,"(i0)") neta + 1
        WRITE(tp,"(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
        format_string = "("//TRIM(c) // "E22.14)"
        nameF = 'epsilon2_indabs_' // trim(adjustl(tp)) // 'K.dat'
        OPEN(UNIT = iuindabs, FILE = nameF)
        WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
        WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
        DO iw = 1, nomega
          WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
        ENDDO
        CLOSE(iuindabs)
        ! 
        nameF = 'epsilon2_indabs_lorenz' // trim(adjustl(tp)) // 'K.dat'
        OPEN(UNIT = iuindabs, FILE = nameF)
        WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
        WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
        DO iw = 1, nomega
          WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs_lorenz(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
        ENDDO
        CLOSE(iuindabs)
      ENDDO
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE indabs_main
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE renorm_eig(ikk, ikq, nrr_k, dims, ndegen_k, irvec_k, irvec_r, cufkk, cufkq, cfac, cfacq)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the renormalization of the eigenenergies to be applied
    !! in case one read external eigenvalues.
    !! The implementation follows Eq. 30 of  Phys. Rev. B 62, 4927 (2000)
    !! Samuel Ponce, Kyle and Emmanouil Kioupakis
    !!
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : xkfd, chw, chw_ks, etf_ks, etf, vmef, nkqf
    USE epwcom,        ONLY : use_ws, nbndsub
    USE constants_epw, ONLY : eps40, ryd2mev, twopi, zero, eps6, ci, czero
    USE wan2bloch,     ONLY : hamwan2bloch, vmewan2bloch
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ikk
    !! Current k point on that core (ikk = 2 * ik + 1)
    INTEGER, INTENT(in) :: ikq
    !! k+q point on that core
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors for electrons
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !!  Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(inout) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), INTENT(inout) :: cufkq(nbndsub, nbndsub)
    !! the same, for points k+q
    COMPLEX(KIND = DP)  :: cufkkd(nbndsub,nbndsub)
    !! Rotation matrix, shifted mesh, points k
    COMPLEX(KIND = DP)  :: cufkqd(nbndsub,nbndsub)
    !! Rotation matrix, shifted mesh, k+q
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr_k, dims, dims)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: cfacq(nrr_k, dims, dims)
    !! Exponential factor
    !
    ! Local variables
    INTEGER :: icounter
    !! Integer counter for displaced points
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: ir
    !! Counter for WS loop
    INTEGER :: iw
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: iw2
    !! Counter on bands when use_ws == .TRUE.
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: rdotk2(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: etfd(nbndsub, nkqf, 6)
    !! interpolated eigenvalues (nbnd, nkqf) eigenvalues for shifted grid in the case of eig_read
    REAL(KIND = DP) :: etfd_ks(nbndsub, nkqf, 6)
    !! interpolated eigenvalues (nbnd, nkqf) KS eigenvalues for shifted grid in the case of eig_read
    COMPLEX(KIND = DP) :: cfacd(nrr_k, dims, dims, 6)
    !! Used to store $e^{2\pi r \cdot k}$ exponential of displaced vector
    COMPLEX(KIND = DP) :: cfacqd(nrr_k, dims, dims, 6)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential of dispaced vector
    !
    cfacd(:, :, :, :) = czero
    cfacqd(:, :, :, :)= czero
    DO icounter = 1, 6
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkfd(:, ikk, icounter), 1, 0.0_DP, rdotk, 1)
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkfd(:, ikq, icounter), 1, 0.0_DP, rdotk2, 1)
      IF (use_ws) THEN
        DO iw = 1, dims
          DO iw2 = 1, dims
            DO ir = 1, nrr_k
              IF (ndegen_k(ir, iw2, iw) > 0) THEN
                cfacd(ir, iw2, iw, icounter)  = EXP(ci * rdotk(ir))  / ndegen_k(ir, iw2, iw)
                cfacqd(ir, iw2, iw, icounter) = EXP(ci * rdotk2(ir)) / ndegen_k(ir, iw2, iw)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        cfacd(:, 1, 1, icounter)  = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
        cfacqd(:, 1, 1, icounter) = EXP(ci * rdotk2(:)) / ndegen_k(:, 1, 1)
      ENDIF
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cufkkd, etfd(:, ikk, icounter), chw, cfacd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkqd, etfd(:, ikq, icounter), chw, cfacqd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkkd, etfd_ks(:, ikk, icounter), chw_ks, cfacd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkqd, etfd_ks(:, ikq, icounter), chw_ks, cfacqd, dims)
    ENDDO ! icounter
    ! -----------------------------------------------------------------------------------------
    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), &
                      etf(:, ikk), etf_ks(:, ikk), chw_ks, cfac, dims)
    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkq, vmef(:, :, :, ikq), &
                      etf(:, ikq), etf_ks(:, ikq), chw_ks, cfacq, dims)
    ! To Satisfy Phys. Rev. B 62, 4927-4944 (2000) , Eq. (30)
    DO ibnd = 1, nbndsub
      DO jbnd = 1, nbndsub
        IF (ABS(etfd_ks(ibnd, ikk, 1) - etfd_ks(jbnd, ikk, 2)) > eps6) THEN
          vmef(1, ibnd, jbnd, ikk) = vmef(1, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 1) - etfd(jbnd, ikk, 2)) / &
              (etfd_ks(ibnd, ikk, 1) - etfd_ks(jbnd, ikk, 2))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikk, 3) - etfd_ks(jbnd, ikk, 4)) > eps6) THEN
          vmef(2, ibnd, jbnd, ikk) = vmef(2, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 3) - etfd(jbnd, ikk, 4)) / &
              (etfd_ks(ibnd, ikk, 3) - etfd_ks(jbnd, ikk, 4))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikk, 5) - etfd_ks(jbnd, ikk, 6)) > eps6) THEN
          vmef(3, ibnd, jbnd, ikk) = vmef(3, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 5) - etfd(jbnd, ikk, 6)) / &
              (etfd_ks(ibnd, ikk, 5) - etfd_ks(jbnd, ikk, 6))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 1) - etfd_ks(jbnd, ikq, 2)) > eps6) THEN
          vmef(1, ibnd, jbnd, ikq) = vmef(1, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 1) - etfd(jbnd, ikq, 2)) / &
              (etfd_ks(ibnd, ikq, 1) - etfd_ks(jbnd, ikq, 2))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 3) - etfd_ks(jbnd, ikq, 4)) > eps6) THEN
          vmef(2, ibnd, jbnd, ikq) = vmef(2, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 3) - etfd(jbnd, ikq, 4)) / &
              (etfd_ks(ibnd, ikq, 3) - etfd_ks(jbnd, ikq, 4))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 5) - etfd_ks(jbnd, ikq, 6)) > eps6) THEN
          vmef(3, ibnd, jbnd, ikq) = vmef(3, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 5) - etfd(jbnd, ikq, 6) ) / &
              (etfd_ks(ibnd, ikq, 5) - etfd_ks(jbnd, ikq, 6))
        ENDIF
      ENDDO
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE renorm_eig
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fermi_carrier_indabs(itemp, etemp_fca, ef0_fca)
    !-----------------------------------------------------------------------
    !!
    !! This is a slightly modified version of subroutine fermicarrier
    !! in Utility.f90. This is used to calculate the fermi energy associated
    !! with a certain carrier density in when the user want to calculate
    !! indirect optical absorption for doped semiconductor
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : omega, alat, at
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : etf, nkf, wkf, efnew, nkqf
    USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm, eps5, kelvin2eV, &
                              zero, eps80, eps6
    USE noncollin_module, ONLY : noncolin
    USE pwcom,     ONLY : nelec
    USE epwcom,    ONLY : int_mob, nbndsub, nc_indabs, nstemp, fermi_energy, &
                          system_2d, carrier, efermi_read, assume_metal, ngaussw
    USE klist_epw, ONLY : isk_dummy
    USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
!    INTEGER, INTENT(out) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: etemp_fca
    !! Temperature in kBT [Ry] unit.
    REAL(KIND = DP), INTENT(inout) :: ef0_fca(nstemp)
    !! Fermi level for the temperature itemp
!    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Second fermi level for the temperature itemp
    REAL(KIND = DP), EXTERNAL :: efermig
    !! External function to calculate the fermi energy
    !
    ! Local variables
    INTEGER :: i
    !! Index for the bisection iteration
    INTEGER :: ik
    !! k-point index per pool
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ivbm
    !! Index of the VBM
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER, PARAMETER :: maxiter = 500 ! 300
    !! Maximum interation
    REAL(KIND = DP) :: fermi
    !! Fermi level returned
    REAL(KIND = DP) :: fermicb
    !! Fermi level returned for second Fermi level
    REAL(KIND = DP) :: fnk
    !! Fermi-Diract occupation
    REAL(KIND = DP) :: ks_exp(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT
    REAL(KIND = DP) :: ks_expcb(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT for CB
    REAL(KIND = DP) :: fermi_exp
    !! Fermi level in exponential format
    REAL(KIND = DP) :: rel_err
    !! Relative error
    REAL(KIND = DP) :: factor
    !! Factor that goes from number of carrier per unit cell to number of
    !! carrier per cm^-3
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: evbm
    !! Energy of the VBM
    REAL(KIND = DP) :: ecbm
    !! Energy of the CBM
    REAL(KIND = DP) :: ef_tmp
    !! Energy of the current Fermi level for the bisection method
    REAL(KIND = DP) :: elw
    !! Energy lower bound for the bisection method
    REAL(KIND = DP) :: eup
    !! Energy upper bound for the bisection method
    REAL(KIND = DP) :: hole_density
    !! Hole carrier density
    REAL(KIND = DP) :: electron_density
    !! Electron carrier density
    REAL(KIND = DP) :: intrinsic_density
    !! Intrinsic carrier density
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    !
    IF (assume_metal) THEN
      !! set conduction band chemical potential to 0 since it is irrelevent
 !     efcb(itemp) = zero
      ef0_fca(itemp) = efermig(etf, nbndsub, nkqf, nelec, wkf, etemp_fca, ngaussw, 0, isk_dummy)
      RETURN
    ENDIF
    ef_tmp  = zero
    fermi   = zero
    fermicb = zero
    inv_cell = 1.0d0 / omega
    !
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF (system_2d) inv_cell = inv_cell * at(3, 3) * alat
    ! vbm index
    IF (noncolin) THEN
      ivbm = FLOOR(nelec / 1.0d0)
    ELSE
      ivbm = FLOOR(nelec / 2.0d0)
    ENDIF
    icbm = ivbm + 1 ! Nb of bands
    !
    ! Initialization value. Should be large enough ...
    evbm = -10000d0
    ecbm = 10000d0 ! In Ry
    !
    WRITE(stdout, '(5x, "calculating band extrema...")')
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      DO ibnd = 1, nbndsub
        IF (ibnd < ivbm + 1) THEN
          IF (etf(ibnd, ikk) > evbm) THEN
            evbm = etf(ibnd, ikk)
          ENDIF
        ENDIF
       ! Find cbm index
        IF (ibnd > ivbm) THEN
          IF (etf(ibnd, ikk) < ecbm) THEN
            ecbm = etf(ibnd, ikk)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    ! Find max and min across pools
    !
    CALL mp_max(evbm, inter_pool_comm)
    CALL mp_min(ecbm, inter_pool_comm)
    !
    IF (itemp == 1) THEN
      WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
      WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
    ENDIF
    !
    ! We first calculate intrinsic carrier density
    ! Using this we can determine if the user input make sense
    WRITE(stdout, '(5x, "Calculating intrinsic carrier density")')
    ef_tmp = (ecbm + evbm) / 2.d0
    eup = ecbm
    elw = evbm
    factor = inv_cell * (bohr2ang * ang2cm)**(-3.d0)
    Do i = 1, maxiter
      !
      ef_tmp = (eup + elw) / 2.d0
      hole_density = zero
      electron_density = zero
      !
      DO ik = 1, nkf
        ikk = 2 * ik -1
!        WRITE(stdout, '(5x, i, i)') icbm, nbndsub
        DO ibnd = icbm, nbndsub
          ekk = etf(ibnd, ikk) - ef_tmp
          fnk = wgauss(-ekk / etemp_fca, -99)
          ! The wkf(ikk) already include a factor 2
          electron_density = electron_density + wkf(ikk) * fnk * factor
!          WRITE(stdout, '(5x, f, f)') ekk, electron_density
        ENDDO ! ibnd
      ENDDO ! ik
      DO ik = 1, nkf
        ikk = 2 * ik -1
        DO ibnd = 1, ivbm
          ekk = etf(ibnd, ikk) - ef_tmp
          fnk = wgauss(-ekk / etemp_fca, -99)
          ! The wkf(ikk) already include a factor 2
          hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
        ENDDO ! ibnd
      ENDDO ! ik
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(hole_density, inter_pool_comm)
      CALL mp_sum(electron_density, inter_pool_comm)
      IF (ABS(hole_density) < eps80) THEN
        rel_err = -1.d0
      ELSE
        rel_err = (hole_density - electron_density) / hole_density
      ENDIF
      IF (ABS(rel_err) < eps5) THEN
        intrinsic_density = hole_density
        EXIT
      ELSEIF ((rel_err) > eps5) THEN
        elw = ef_tmp
      ELSE
        eup = ef_tmp
      ENDIF
      intrinsic_density = hole_density
    ENDDO ! maxiter
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x, a)') "Too many iterations when calculating intrinsic density"
      WRITE(stdout, '(5x, a, f8.5)') "Relative error of hole density versus electron density: ", rel_err
      WRITE(stdout, '(5x, a)') "Something likely wrong unless we are dealing with a metallic system"
    ENDIF
    !
    WRITE(stdout, '(5x, "Intrinsic density = ", E18.6, "Cm^-3")' ) intrinsic_density
    WRITE(stdout, '(/5x, "calculating fermi level...")')
    IF (ABS(nc_indabs) < intrinsic_density) THEN
      WRITE(stdout, '(5x, a)') 'nc_indabs not given, or smaller than intrinsic density'
      WRITE(stdout, '(5x, a)') 'Setting the fermi level to mig-gap'
      fermi = (evbm + ecbm) / 2.d0
    ELSEIF ( nc_indabs > intrinsic_density ) THEN
      ! assuming free electron density
      eup = 10000d0 + ecbm
      elw = evbm - 10000d0
      ef_tmp = (ecbm + ecbm) / 2.d0
      DO i = 1, maxiter
!        WRITE(stdout, '(5x, i)') nkf
        electron_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik -1
!          WRITE(stdout, '(5x, i, i)') icbm, nbndsub
          DO ibnd = icbm, nbndsub
            ekk = etf(ibnd, ikk) - ef_tmp
            fnk = wgauss(-ekk / etemp_fca, -99)
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk * factor
!            WRITE(stdout, '(5x, f, f)') ekk, electron_density
          ENDDO ! ibnd
        ENDDO ! ik
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(electron_density, inter_pool_comm)
!        WRITE(stdout, '(5x, f)') electron_density
!        IF (ABS(electron_density) < eps80) THEN
!          rel_err = -1.0d0
!        ELSE
          ! In this case ncarrier is a negative number
         rel_err = (electron_density - nc_indabs) / electron_density
!        ENDIF
        !
        IF (ABS(rel_err) < eps6) THEN
          fermi = ef_tmp
          EXIT
        ELSEIF ((rel_err) > eps6) THEN
          eup = ef_tmp
        ELSE
          elw = ef_tmp
        ENDIF
        IF (MOD(i, 10) == 0) THEN
          WRITE(stdout, '(5x, a, i, a, f8.3)') 'Iteration #', i, " Fermi level = ", ef_tmp
        ENDIF
        ef_tmp = (eup + elw) / 2.0d0
      ENDDO ! maxiter
      fermi = ef_tmp
    ELSEIF ( nc_indabs < - intrinsic_density) THEN
      ! assuming free hole density
      eup = 10000d0 + ecbm
      elw = evbm - 10000d0
      ef_tmp = (ecbm + evbm) / 2.d0
      DO i = 1, maxiter
        hole_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik -1
          DO ibnd = 1, ivbm
            ekk = etf(ibnd, ikk) - ef_tmp
            fnk = wgauss(-ekk / etemp_fca, -99)
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
          ENDDO ! ibnd
        ENDDO ! ik
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(hole_density, inter_pool_comm)
!        IF (ABS(hole_density) < eps80) THEN
!          rel_err = -1000.0d0
!        ELSE
          ! In this case ncarrier is a negative number
        rel_err = (hole_density - ABS(nc_indabs)) / hole_density
!        ENDIF
        !
        IF (ABS(rel_err) < eps6) THEN
          fermi = ef_tmp
          EXIT
        ELSEIF ((rel_err) > eps6) THEN
          elw = ef_tmp
        ELSE
          eup = ef_tmp
        ENDIF
        IF (MOD(i, 10) == 0) THEN
          WRITE(stdout, '(5x, a, i, a, f8.3)') 'Iteration #', i, " Fermi level = ", ef_tmp
        ENDIF
        ef_tmp = (eup + elw) / 2.0d0
      ENDDO ! maxiter
      fermi = ef_tmp
    ENDIF !ncarrier
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x, "Warning: too many iterations in bisection"/ &
                    5x, "ef_tmp = ", f10.6)' ) fermi * ryd2ev
    ENDIF
    WRITE(stdout, '(/5x, "Temperature ", f8.3, " K")' ) etemp_fca * ryd2ev / kelvin2eV
    IF (nc_indabs > intrinsic_density) THEN
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Electron density = ", E18.6, "Cm^-3")' ) electron_density
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ELSEIF (nc_indabs < - intrinsic_density) THEN
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Hole density = ", E18.6, "Cm^-3")' ) hole_density
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ELSE
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ENDIF
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE fermi_carrier_indabs
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE indabs
  !-----------------------------------------------------------------------


