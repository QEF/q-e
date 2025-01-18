  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
  ! Feliciano Giustino
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE screening
  !----------------------------------------------------------------------
  !! 
  !! This module contains routines for screening effect
  !!
  USE kinds, ONLY: DP
  !!
  IMPLICIT NONE
  CONTAINS
    !
    !--------------------------------------------------------------------------
    COMPLEX(KIND = DP) FUNCTION H_eps(z)
    !--------------------------------------------------------------------------
    !!
    !! Function used in the Lindhard function. See Eq.(56) of Hedin (1965)
    !!
    USE ep_constants,  ONLY : eps10
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: z
    !! Argument of the Lindhard function
    !
    IF (ABS(z - 1.d0) > eps10) THEN
      IF (ABS((z + 1.d0) / (z - 1.d0)) > eps10) THEN
        H_eps = 2.d0 * z + (1.d0 - z**2) * LOG((z + 1.d0) / (z - 1.d0))
      ENDIF
    ENDIF
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END FUNCTION H_eps
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------------
    SUBROUTINE rpa_epsilon(q, w, nmodes, epsil, eps_rpa)
    !-----------------------------------------------------------------------------
    !!
    !!  Compute the Lindhard dielectric function for the homogeneous electron gas
    !!
    USE cell_base,     ONLY : at, bg, omega, alat
    USE ep_constants,  ONLY : twopi, ha2ev, cone, ci, eps5, eps10
    USE ep_constants,  ONLY : pi
    USE input,         ONLY : meff, fermi_diff, nel, smear_rpa, system_2d
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Number of phonon modes
    REAL (KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates
    REAL (KIND = DP), INTENT(inout) :: w(nmodes)
    !! phonon frequencies associated with q
    REAL (KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    COMPLEX(KIND = DP), INTENT(out) :: eps_rpa(nmodes)
    !! electronic screening
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    INTEGER :: im
    !! Mode counter
    REAL(KIND = DP) :: n
    !! Electron density in atomic units
    REAL(KIND = DP) :: rs
    !! Prefactor for the dielectric screening
    REAL(KIND = DP) :: EF
    !! Fermi-level in eV
    REAL(KIND = DP) :: kF
    !! Fermi wavevector
    REAL(KIND = DP) :: pref
    !! Prefactor for the dielectric function
    REAL(KIND = DP) :: eta
    !! Broadening for the dielectric function
    REAL(KIND = DP) :: q2
    !! q-point square
    REAL(KIND = DP) :: qm
    !! Internal units for Hedin's formula
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    COMPLEX(KIND = DP) :: u
    !! Complex frequency argument
    !
    IF (system_2d == 'no') THEN
      n = nel / omega
    ELSE
      n = (at(3, 3) * alat) * nel / omega
    ENDIF
    EF = fermi_diff / ha2ev
    kF = (3.d0 * pi**2 * n)**(1.d0 / 3.d0)
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    rs  = (3.d0 / ( 4.d0 * pi * n ) )**(1.d0 / 3.d0) * meff / eps_ave
    w = w * 0.5d0 / EF / 4.d0 !Ha&internal units for Hedin's formula
    pref = (4.d0 / 9.d0 / pi )**(1.d0 / 3.0) * (rs / 8.d0 / pi)
    eta = smear_rpa / ha2ev / EF / 4.d0
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of Lindhard screening: use with care")')
      WRITE(stdout, '(5x,"Warning: current implementation for doubly degenerate band, one valley")')
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'Nel = ', nel, ', n = ', n, ' au^-3, meff = ', meff
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'EF = ', EF * ha2ev, ' eV, kF = ', kF, ' au^-1, rs = ', rs
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    IF (eps_ave < eps5) eps_ave = 1.d0
    !
    CALL cryst_to_cart(1, q, bg, 1)
    q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) * (twopi / alat) / kF / 2.d0 ! internal units for Hedin's formula
    !
    IF (ABS(qm) > eps10) THEN
      DO im = 1, nmodes
        u = w(im) + SIGN(eta, w(im)) * ci
        eps_rpa(im) = 1.d0 + pref * (H_eps(qm + u / qm) + H_eps(qm - u / qm)) / qm**3
      ENDDO
    ELSE
      eps_rpa = cone
    ENDIF
    !
    w = w / (0.5d0 / EF / 4.d0)
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE rpa_epsilon
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE tf_epsilon(q, nmodes, epsil, eps_tf)
    !--------------------------------------------------------------------------
    !!
    !!  Compute the Thomas-Fermi dielectric screening
    !!
    USE cell_base,     ONLY : at, bg, omega, alat
    USE ep_constants,  ONLY : twopi, ha2ev, cone, eps5, eps10
    USE ep_constants,  ONLY : pi
    USE input,         ONLY : fermi_diff, nel
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Number of phonon modes
    REAL(KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates)
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    COMPLEX(KIND = DP), INTENT(out) :: eps_tf(nmodes)
    !! electronic screening
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    REAL(KIND = DP) :: n
    !! Electron density in atomic units
    REAL(KIND = DP) :: EF
    !! Fermi-level in eV
    REAL(KIND = DP) :: q2
    !! q-point square
    REAL(KIND = DP) :: qtf
    !! Thomas-Fermi wavector
    REAL(KIND = DP) :: qtfc
    !! Thomas-Fermi wavector in unit of 2pi/a
    REAL(KIND = DP) :: qm
    !! Modulus of q
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    !
    n = nel / omega
    EF = fermi_diff / ha2ev
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    qtf = (6.d0 * pi * n / EF / eps_ave )**(1.d0 / 2.d0)
    qtfc = qtf / (twopi / alat)
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of Thomas-Fermi screening: use with care")')
      WRITE(stdout, '(5x,"Warning: current implementation for doubly degenerate band, one valley")')
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'Nel = ', nel, ', n = ', n, ' au^-3, EF (eV) = ', EF*ha2ev
      WRITE(stdout, '(5x,a,f12.8,a,f12.8)') 'q_tf (au-1) = ', qtf, ', q_tf (tpiba) = ', qtfc
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    IF (eps_ave < eps5) eps_ave = 1.d0
    !
    CALL cryst_to_cart(1, q, bg, 1)
    q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) ! in tpiba
    IF (ABS(qm) > eps10) THEN
      eps_tf = 1.d0 + qtfc**2 / q2
    ELSE
      eps_tf = cone
    ENDIF
    !
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE tf_epsilon
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE calc_qtf2_therm(itemp, etemp, ef0, efcb, ctype)
    !--------------------------------------------------------------------------
    !!
    !! JL 01/2022
    !! Compute thermal Thomas fermi wave vector for doped Semiconductors
    !! ctype = +/1 1 only!! assume_metal = .FALSE. only (for now...)
    !! q_{tf}^{2}(T) = \frac{4\pi}{\omega*\varepsilon_\infty}\sum_{n,\mathbf{k}}
    !! w_{n\mathbf k}\frac{\partial f_{n,\mathbf k}}{\partial \epsilon_{n\mathbf k}}
    !!
    !
    USE cell_base,     ONLY : omega
    USE input,         ONLY : nstemp, nbndsub
    USE global_var,    ONLY : etf, nkf, wkf, qtf2_therm, epsi_s
    USE ep_constants,  ONLY : pi
    USE ep_constants,  ONLY : eps12, one, thre
    USE io_global,     ONLY : stdout
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! temp index
    INTEGER, INTENT(in) :: ctype
    !! -1 for holes, 1 for electrons, 0 for both (forbidden)
    REAL(KIND = DP), INTENT(in)  :: etemp
    !! energy of temperatre at index itemp
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! fermi energy of p-type or metal material
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! fermi_energy of n-type materials
    !
    ! local variables
    INTEGER :: ik
    !! counter on k ind
    INTEGER :: ikk
    !! counter over k on fine grid
    INTEGER :: ibnd
    !! counter over bands
    REAL(KIND=DP) :: ekk
    !! Energy on fine grid at ekk in fermi window
    REAL(KIND=DP) :: dfnk
    !! derivative of the fermi-dirac distribution
    REAL(KIND=DP) :: inv_etemp
    !! beta = 1/(kB*T), inverse thermal temperature in Ryd
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    !
    eps_ave = (epsi_s(1,1) + epsi_s(2,2) + epsi_s(3,3)) / thre
    !
    IF (ctype==0) THEN
      WRITE(stdout, '(5x,"ctype=0 not supported, keeping epf_tf_therm(:) == 1.0d0")')
      RETURN
    ENDIF
    !
    !IF (my_pool_id == 0) THEN
    !  IF (ctype == -1) THEN
    !    WRITE(stdout, '(5x,"Using Fermi energy of: ")')
    !    WRITE(stdout,*) ef0(itemp)
    !  ELSEIF (ctype == 1) THEN
    !    WRITE(stdout, '(5x,"Using Fermi energy of: ")')
    !    WRITE(stdout,*) efcb(itemp)
    !  ENDIF
    !ENDIF
    !
    inv_etemp = one / etemp
    !
    ! calc qtf^2(T) for holes in case of ctype=-1
    IF (ctype==-1) THEN
      DO ik = 1, nkf
        ikk = 2*ik-1
        DO ibnd = 1, nbndsub
          ekk  = etf(ibnd, ikk)
          dfnk = w0gauss((ekk-ef0(itemp))*inv_etemp,-99)*inv_etemp
          qtf2_therm(itemp) = qtf2_therm(itemp) + (1.0d0/(eps_ave*omega))*8.0d0*pi*( wkf(ikk)*dfnk)
        ENDDO
      ENDDO
    ENDIF
    !
    ! calc qtf^2(T) for holes in case of ctype=-1
    IF (ctype==1) THEN
      DO ik = 1, nkf
        ikk = 2*ik-1
        DO ibnd = 1, nbndsub
          ekk  = etf(ibnd, ikk)
          dfnk = w0gauss((ekk-efcb(itemp))*inv_etemp,-99)*inv_etemp
          qtf2_therm(itemp) = qtf2_therm(itemp) + (1.0d0/(omega*eps_ave))*8.0d0*pi*( wkf(ikk)*dfnk  )
        ENDDO
      ENDDO
    ENDIF
    !
    CALL mp_sum(qtf2_therm(itemp), inter_pool_comm)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_qtf2_therm
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE calc_epstf_therm(q, nstemp, epsil)
    !--------------------------------------------------------------------------
    !!
    !! JL: 01/2022
    !! Calculate the thermal thomas fermi dielectric function (1 + qtf^2(T)/q^2)
    !!
    USE cell_base,     ONLY : at, bg, alat
    USE ep_constants,  ONLY : twopi, ha2ev, cone, eps5, eps10, one
    USE global_var,    ONLY : qtf2_therm, epstf_therm
    !
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nstemp
    !! number of temperature points
    REAL(KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates)
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    INTEGER :: itemp
    !! counter of temperature index
    INTEGER :: m1, m2, m3
    !! counter for G vecs
    INTEGER :: m1f, m2f, m3f
    !! index of G that minimizes q+G
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    REAL(KIND = DP) :: q2
    !! squared phonon wavevector
    REAL(KIND = DP) :: q2inv
    !! inverse of q2, which is summed over G
    REAL(KIND = DP) :: qm
    !! sqrt of q2
    REAL(KIND = DP) :: qtfc2
    !! square Thomas fermi wv in units of (two*pi/alat)
    !! qtfc = qtf / (twopi / alat)
    REAL(KIND = DP) :: val
    !! test value for finding G that minimizes q+G
    REAL(KIND = DP) :: g1, g2, g3, qeq, qmG
    !!
    REAL(KIND = DP) :: eps_loc(3,3)
    !
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    eps_loc(1,1) = 1.0d0
    eps_loc(2,2) = 1.0d0
    eps_loc(3,3) = 1.0d0
    eps_loc(1,2) = 0.0d0
    eps_loc(2,1) = 0.0d0
    eps_loc(1,3) = 0.0d0
    eps_loc(3,1) = 0.0d0
    eps_loc(2,3) = 0.0d0
    eps_loc(3,2) = 0.0d0
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of thermal Thomas-Fermi screening: use with care")')
      DO itemp = 1, nstemp
        WRITE(stdout, '(5x,a,i5)') 'itemp=', itemp
        WRITE(stdout, '(5x,a,f22.16)') 'q_tf (au^-1) = ', SQRT(qtf2_therm(itemp))
      ENDDO
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    !
    CALL cryst_to_cart(1, q, bg, 1)
    !
    ! JL: Look for G that minimizes q+G --> enforce periodicity
    val = 1.0d24
    m1f = 0
    m2f = 0
    m3f = 0
    !
    DO m1 = -5, 5
      DO m2 = -5, 5
        DO m3 = -5, 5
          g1 = q(1) + (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3))
          g2 = q(2) + (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3))
          g3 = q(3) + (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3))
          qmG = SQRT(g1**2+g2**2+g3**2)
          IF (qmG < val) THEN
            val = qmG
            m1f = m1
            m2f = m2
            m3f = m3
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    q2inv = 0.0d0
    !
    DO m1 = -1, 1 !-nqc1, nqc1
      DO m2 = -1, 1 !-nqc2, nqc2
        DO m3 = -1, 1 !-nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1) + m1f * bg(1, 1) + &
                  m2f * bg(1, 2) + m3f * bg(1, 3)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2) + m1f * bg(2, 1) + &
                  m2f * bg(2, 2) + m3f * bg(2, 3)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3) + m1f * bg(3, 1) + &
                  m2f * bg(3, 2) + m3f * bg(3, 3)
          !
          !qeq = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3)
          !+ &
          !       g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3)
          !       + &
          !       g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))
          !       !*twopi/alat
          qeq = (g1 * (eps_loc(1, 1) * g1 + eps_loc(1, 2) * g2 + eps_loc(1, 3) * g3) + &
                 g2 * (eps_loc(2, 1) * g1 + eps_loc(2, 2) * g2 + eps_loc(2, 3) * g3) + &
                 g3 * (eps_loc(3, 1) * g1 + eps_loc(3, 2) * g2 + eps_loc(3, 3) * g3)) !*twopi/alat
          IF (qeq > 0.0d0) THEN
            q2inv = q2inv + (1.0d0 / qeq)
          ELSE
            q2inv = q2inv + 0.0d0
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! JL now calculate contributions to epstf_therm(itmep)
    !
    IF (q2inv > 0.0d0) THEN
      q2 = 1.0d0 / q2inv
    ELSE
      q2 = 0.0
    ENDIF
    !q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) ! in tpiba
    DO itemp = 1, nstemp
      IF (ABS(qm) > eps10) THEN
        qtfc2 = qtf2_therm(itemp) / ((twopi / alat)**2.0d0)
        epstf_therm(itemp) = 1.d0 + (qtfc2 / q2)
      ELSE
        epstf_therm(itemp) = one
      ENDIF
    ENDDO
    !
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_epstf_therm
    !--------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE calc_eps_static(nrr_q, irvec_q, ndegen_q, rws, nrws)
    !-------------------------------------------------------------------------------
    !!
    !! Compute low-frequency dielectric tensor which includes the contributions
    !! from both electronic and ionic polarizations.
    !! The calculations are done based on equations (1)-(3) in the article
    !! [Rignanese, Gonze and Pasquarello, PRB 63, 104305 (2001)] (Ryd unit).
    !! In SI unit, this quantity is calculated as
    !! $$\epsilon_{\alpha\beta}^{0} = \epsilon_{\alpha\beta}^{\infty} +
    !!                           \sum_{\nu} \epsilon_{\alpha\beta\nu}^{\rm ion}$$
    !! where,
    !! $$\epsilon_{\alpha\beta\nu}^{\rm ion} = \frac{e^2}{4\pi^2\epsilon_0\Omega}
    !!              \frac{1}{\omega_{\nu}^2}
    !!              \sum\limits_{\gamma\kappa} Z^{\ast}_{\alpha\gamma\kappa}
    !!              \frac{{\bf e}_{\nu\gamma\kappa}}{\sqrt{M_{\kappa}}}
    !!              \sum\limits_{\gamma'\kappa'} Z^{\ast}_{\beta\gamma'\kappa'}
    !!              \frac{{\bf e}_{\nu\gamma'\kappa'}}{\sqrt{M_{\kappa'}}}$$
    !! Added by V.-A. Ha (2024)
    !!
    USE global_var,     ONLY : epsi, epsi_s, zstar
    USE ep_constants,   ONLY : zero, czero, fpi, eps8, echg, eps_vac, bohr, amu, &
                               AMU_RY, sm1toryd, bohr_radius_si, pi
    USE cell_base,      ONLY : omega
    USE ions_base,      ONLY : amass, nat, ityp
    USE modes,          ONLY : nmodes
    USE input,          ONLY : lpolar, lifc
    USE io_global,      ONLY : stdout
    USE wannier2bloch,  ONLY : dynwan2bloch, dynifc2blochf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)           :: nrr_q
    !! Number of WS points for phonons
    INTEGER, INTENT(in)           :: irvec_q(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, INTENT(in)           :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that
    !! depend on atomic positions $R + \tau(nb) - \tau(na)$
    REAL(KIND = DP), INTENT(in)   :: rws(0:3, 200)
    !! Number of real-space Wigner-Seitz vectors
    INTEGER, INTENT(in)           :: nrws
    !! Number of real-space Wigner-Seitz
    INTEGER :: ierr
    !! Error index when reading/writing a file or allocating array
    INTEGER :: i, j, k
    !! Cartesian direction indices
    INTEGER :: mu, nu
    !! Phonon mode indices
    INTEGER :: na
    !! Atom index
    REAL(KIND = DP)  :: mass
    !! mass of atom
    REAL(KIND = DP)  :: xxq_r(3)
    !! q-point coordinate
    REAL(KIND = DP)  :: fac
    !! factor e^2/(4pi*epsilon_0*omega*M_0) with M_0 being amu
    COMPLEX(KIND = DP)  :: zstar_ph_a
    !! Average Born effective charge for each phonon mode
    COMPLEX(KIND = DP)  :: zstar_ph_b
    !! Average Born effective charge for each phonon mode
    REAL(KIND = DP)  :: epsi_ion(3, 3)
    !! ionic part of dielectric tensor
    REAL(KIND = DP)  :: w2(nmodes)
    !! Square of phonon frequencies at Gamma point
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvector
    !
    ALLOCATE(epsi_s(3,3), STAT=ierr)
    IF (ierr /= 0) CALL errore('calc_eps_static', 'Error allocating epss', 1)
    epsi_s(:,:) = epsi(:,:)
    IF (.NOT. lpolar) GOTO 425
    ! --------------------------------------------------------------
    ! V.-A. Ha, we need phonon energy and eigenvector at Gamma point
    ! to compute ionic part of low-frequency dielectric tensor
    ! dynamical matrix : Wannier -> Bloch
    ! phonon energies w2 in Ry2, eigenvector uf is dimensionless
    ! --------------------------------------------------------------
    !
    xxq_r = (/zero, zero, zero/)
    uf = czero
    w2 = zero
    IF (.NOT. lifc) THEN
      CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, .false.)
    ELSE
      CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, .false.)
    ENDIF
    !
    epsi_ion = zero
    ! omega in bohr^3, M_{\kappa} in atomic mass unit (Atomic Rydberg unit)
    fac = echg**2 / (pi * fpi * eps_vac * omega * bohr_radius_si**3 * amu) ! in s-2
    ! convert frequencies square in s-2 to Ry2 for consistency with w2
    fac = fac * sm1toryd**2
    DO nu = 1, nmodes
      ! calculate average Born effective charge for each phonon mode
      ! Note: The implementation is more general than in "subroutine
      ! polar_mode_permittivity"
      ! of "LR_Modules/dynmat_sub.f90" when runing dynmat.x in which
      ! eigenvectors "uf" are real 
      DO i = 1, 3
        ! Born effective charge for each phonon mode \nu along direction i
        zstar_ph_a = czero
        DO na = 1, nat
          mass = amass(ityp(na))/AMU_RY  ! Atomic Rydberg unit to atomic mass unit
          DO k = 1, 3
            mu = (na - 1) * 3 + k
            zstar_ph_a = zstar_ph_a + zstar(i, k, na) * CONJG(uf(mu, nu)) / DSQRT(mass)
          ENDDO
        ENDDO
        !
        DO j = 1, 3
          ! Born effective charge for each phonon mode \nu along direction j
          zstar_ph_b = czero
          DO na = 1, nat
            mass = amass(ityp(na))/AMU_RY  ! Atomic Rydberg unit to atomic mass unit
            DO k = 1, 3
              mu = (na - 1) * 3 + k
              zstar_ph_b = zstar_ph_b + zstar(j, k, na) * uf(mu, nu) / DSQRT(mass)
            ENDDO
          ENDDO
          !
          IF (w2(nu) > eps8) epsi_ion(i, j) = epsi_ion(i, j) + fac / w2(nu) * zstar_ph_a * zstar_ph_b
        ENDDO ! j
      ENDDO ! i
      !
      ENDDO ! nmodes
    !
    epsi_s(:,:) = epsi_s(:,:) + epsi_ion(:,:)
    !
425 WRITE(stdout, '(5x,a)') ' '
    WRITE(stdout, '(5x,a)') 'Low-frequency (static) dielectric tensor for ionized-impurity scattering'
    WRITE(stdout, '(10x,3f16.10)') epsi_s(1,:)
    WRITE(stdout, '(10x,3f16.10)') epsi_s(2,:)
    WRITE(stdout, '(10x,3f16.10)') epsi_s(3,:)
    WRITE(stdout, '(5x,a)') ' '
    !
    RETURN
    !-------------------------------------------------------------------------------
    END SUBROUTINE calc_eps_static
    !-------------------------------------------------------------------------------
    !
  !---------------------------------------------------------------------------------
  END MODULE screening
  !---------------------------------------------------------------------------------
