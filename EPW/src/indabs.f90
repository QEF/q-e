  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE indabs(iq)
  !-----------------------------------------------------------------------
  !! 
  !!  Phonon assisted absorption
  !!  12/03/2018 E. Kioupakis: First implementation
  !!  08/04/2018 S. Ponce: Cleaning 
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode_id
  USE io_epw,        ONLY : iuindabs
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy,&
                            vme, omegamin, omegamax, omegastep, n_r, scissor, eig_read
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            nkf, epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, efnew, &
                            dmef, omegap, epsilon2_abs, epsilon2_abs_lorenz, vmef, etf_ks
  USE transportcom,  ONLY : lower_bnd, upper_bnd
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, czero
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE cell_base,     ONLY : omega
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: iq
  !! Q-point index    
  !
  ! Local variables 
  CHARACTER (len=256) :: nameF='epsilon2_indabs.dat'
  !! Name of the file
  CHARACTER (len=10) :: c
  !! Number of eta values, in string format
  CHARACTER (len=256) :: format_string
  !! Format string
  !
  LOGICAL :: opnd
  !! Check whether the file is open. 
  ! 
  INTEGER :: ios
  !! integer variable for I/O control
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
  INTEGER :: nrec
  !! Record index for reading the e-f matrix
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: nksqtotf
  !! Total number of k+q points 
  INTEGER :: i
  !! Index for reading files
  INTEGER :: iw
  !! Index for frequency
  INTEGER :: nomega
  !! Number of points on the photon energy axis
  INTEGER :: mbnd
  !! Index for summation over intermediate bands
  INTEGER :: ipol
  !! Polarization direction
  INTEGER :: m
  !! Counter on denominator imaginary broadening values
  INTEGER, PARAMETER :: neta = 9
  !! Broadening parameter
  REAL(kind=DP) :: tmp
  !! Temporary variable to store real part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp2
  !! Temporary variable to store imag part of Sigma for the degenerate average
  REAL(kind=DP) :: tmp3
  !! Temporary variable to store Z for the degenerate average
  REAL(kind=DP) :: ekk2
  !! Temporary variable to the eigenenergies for the degenerate average
  REAL(kind=DP) :: sigmar_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the real-part of Sigma 
  REAL(kind=DP) :: sigmai_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the imag-part of Sigma 
  REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
  !! Temporary array to store the Z
  REAL(kind=DP) :: g2
  !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekmk
  !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
  REAL(kind=DP) :: ekmq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
  REAL(kind=DP) :: wq(nmodes), nqv(nmodes)
  !! Phonon frequencies and phonon occupations on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgq
  !! Bose occupation factor $n_{q\nu}(T)$
  REAL(kind=DP) :: wgkk, wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
  REAL(kind=DP) :: weighta, weighte
  !!- delta function for absorption, emission
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: inv_degaussw
  !! Inverse of the smearing for efficiency reasons
  REAL(KIND=DP) :: pfac
  !! Occupation prefactors
  REAL(KIND=DP) :: pface
  !! Occupation prefactors
  REAL(KIND=DP) :: cfac
  !! Absorption prefactor
  REAL(kind=DP) :: eta(neta) = (/ 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 /)/ryd2eV
  !! Imaginary broadening of matrix element denominators
  REAL(kind=DP), EXTERNAL :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), EXTERNAL :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), EXTERNAL :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:)
  !! Collect k-point coordinate from all pools in parallel case
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Collect eigenenergies from all pools in parallel case
  COMPLEX(KIND=DP) :: vkk(3,ibndmax-ibndmin+1,ibndmax-ibndmin+1)
  !!- Velocity matrix elements at k, k+q
  COMPLEX(KIND=DP) :: vkq(3,ibndmax-ibndmin+1,ibndmax-ibndmin+1)
  !!- Velocity matrix elements at k, k+q
  COMPLEX (KIND=DP) :: s1a(3), s1e(3), s2a(3), s2e(3)
  !! Transition probability function    
  COMPLEX (KIND=DP) :: epf(ibndmax-ibndmin+1, ibndmax-ibndmin+1,nmodes)
  !! Generalized matrix elements for phonon-assisted absorption
  ! 
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  nomega = INT((omegamax - omegamin)/omegastep) + 1
  ! 
  ! 300 K  
  ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 for spin (weights for k-points are divided by 2 to be normalized to 1)
  ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
  ! 
  cfac = 16.d0*pi**2
  ! 
  IF (iq == 1) THEN
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Phonon-assisted absorption")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick .lt. 1.d3 ) &
         WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Temperature T = ',eptemp * ryd2ev, ' eV'
    ! 
    IF ( .not. ALLOCATED (omegap) )    ALLOCATE(omegap(nomega))
    IF ( .not. ALLOCATED (epsilon2_abs) ) ALLOCATE(epsilon2_abs(3,nomega,neta))
    IF ( .not. ALLOCATED (epsilon2_abs_lorenz) ) ALLOCATE(epsilon2_abs_lorenz(3,nomega,neta))
    ! 
    epsilon2_abs = 0.d0
    epsilon2_abs_lorenz = 0.d0
    DO iw = 1, nomega
      omegap(iw) = omegamin + (iw-1) * omegastep
    ENDDO
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  IF ( efermi_read ) THEN
     !
     ef0 = fermi_energy
  ELSE
     !
     ef0 = efnew
  ENDIF
  ! 
  DO ik = 1, nkf
    !
    ikk = 2 * ik - 1
    ikq = ikk + 1
    !
    DO imode = 1, nmodes
      !
      ! the phonon frequency at this q and nu 
      wq(imode) = wf (imode, iq)
      ! 
      epf(:,:,imode) = epf17(:, :, imode,ik)
      IF ( wq(imode) .gt. eps_acustic ) THEN
        nqv(imode) = wgauss( -wq(imode)/(eptemp), -99)
        nqv(imode) = nqv(imode) / ( one - two * nqv(imode) )
      ENDIF
    ENDDO
    !
    ! RM - vme version should be checked
    IF ( vme ) THEN 
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          ! vmef is in units of Ryd * bohr
          vkk(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk)
          vkq(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikq)
        ENDDO
      ENDDO
    ELSE
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
           ! Dme's already corrected for GW corrections in wan2bloch.f90
          vkk(:,ibnd,jbnd) = 2.0 * dmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk) 
          vkq(:,ibnd,jbnd) = 2.0 * dmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikq) 
        ENDDO
      ENDDO
    ENDIF
    ! 
    DO ibnd = 1, ibndmax-ibndmin+1
      !  the energy of the electron at k (relative to Ef)
      ekk = etf (ibndmin-1+ibnd, ikk) - ef0
      !
      IF ( abs(ekk) .lt. fsthick ) THEN
        !
        wgkk = wgauss( -ekk*inv_eptemp0, -99)  
        ! 
        DO jbnd = 1, ibndmax-ibndmin+1
          ! 
          ! The fermi occupation for k+q
          ekq = etf (ibndmin-1+jbnd, ikq) - ef0
          !  
          IF ( abs(ekq) < fsthick .AND. ekq < ekk+wq(nmodes)+omegamax + 6.0*degaussw ) THEN
            !
            wgkq = wgauss ( -ekq*inv_eptemp0, -99)  
            !
            IF ( ekq-ekk-wq(nmodes)-omegamax > 6.0*degaussw ) CYCLE 
            IF ( ekq-ekk+wq(nmodes)-omegamin < 6.0*degaussw ) CYCLE 
            ! 
            DO imode = 1, nmodes
              !
              IF ( wq(imode) > eps_acustic ) THEN
                !
                DO m = 1, neta
                  s1a = czero
                  s1e = czero
                  s2a = czero
                  s2e = czero
                  !
                  DO mbnd = 1, ibndmax-ibndmin+1
                    !
                    ! The energy of the electron at k (relative to Ef)
                    ekmk = etf (ibndmin-1+mbnd, ikk) - ef0
                    ! The energy of the electron at k+q (relative to Ef)
                    ekmq = etf (ibndmin-1+mbnd, ikq) - ef0
                    !
                    s1a(:)  = s1a(:) + epf(mbnd, jbnd,imode) * 0.5 * vkk(:,ibnd, mbnd)  / &
                         (  ekmk  - ekq + wq(imode) + ci * eta(m) )
                    s1e(:)  = s1e(:) + epf(mbnd, jbnd,imode) * 0.5 * vkk(:,ibnd, mbnd)  / &
                         (  ekmk  - ekq - wq(imode) + ci * eta(m) )
                    s2a(:) =  s2a(:) + epf(ibnd, mbnd,imode) * 0.5 * vkq(:,mbnd, jbnd)   / &
                         (  ekmq  - ekk - wq(imode)+ ci * eta(m))
                    s2e(:) =  s2e(:) + epf(ibnd, mbnd,imode) * 0.5 * vkq(:,mbnd, jbnd)   / &
                         (  ekmq  - ekk + wq(imode)+ ci * eta(m))
                  ENDDO
                  ! 
                  pfac  =  nqv(imode)      * wgkk *(one- wgkq ) - (nqv(imode)+one)*(one-wgkk) * wgkq 
                  pface = (nqv(imode)+one) * wgkk *(one- wgkq ) -  nqv(imode)     *(one-wgkk) * wgkq
                  ! 
                  DO iw = 1, nomega
                    !
                    IF ( ABS(ekq-ekk-wq(imode)-omegap(iw)) > 6.0*degaussw .AND. &
                         ABS(ekq-ekk+wq(imode)-omegap(iw)) > 6.0*degaussw) CYCLE
                    ! 
                    weighte = w0gauss( ( ekq - ekk - omegap(iw) + wq(imode))  / degaussw, 0) / degaussw
                    weighta = w0gauss( ( ekq - ekk - omegap(iw) - wq(imode))  / degaussw, 0) / degaussw
                    ! 
                    DO ipol = 1, 3
                      epsilon2_abs(ipol,iw,m) = epsilon2_abs(ipol,iw,m)  + &
                           (wkf(ikk)/2.0) * wqf(iq) * &
                           cfac / omegap(iw)**2 * pfac  * weighta * abs( s1a(ipol) + s2a(ipol) )**2 / (2 * wq(imode) * omega )
                      epsilon2_abs(ipol,iw,m) = epsilon2_abs(ipol,iw,m)  + &
                           (wkf(ikk)/2.0) * wqf(iq) * &
                           cfac / omegap(iw)**2 * pface * weighte * abs( s1e(ipol) + s2e(ipol) )**2 / (2 * wq(imode) * omega )
                      epsilon2_abs_lorenz(ipol,iw,m) = epsilon2_abs_lorenz(ipol,iw,m)  + &
                           (wkf(ikk)/2.0) * wqf(iq) * &
                           cfac / omegap(iw)**2 * pfac  * abs( s1a(ipol) + s2a(ipol) )**2 / (2 * wq(imode) * omega ) * &
                           ( degaussw / ( degaussw**2 + (ekq - ekk - omegap(iw) - wq(imode))**2 ))/pi
                      epsilon2_abs_lorenz(ipol,iw,m) = epsilon2_abs_lorenz(ipol,iw,m)  + &
                           (wkf(ikk)/2.0) * wqf(iq) * &
                           cfac / omegap(iw)**2 * pface * abs( s1e(ipol) + s2e(ipol) )**2 / (2 * wq(imode) * omega ) * &
                           ( degaussw / ( degaussw**2 + (ekq - ekk - omegap(iw) + wq(imode))**2 ))/pi
                    ENDDO ! ipol
                    !
                  ENDDO ! iw
                  !
                ENDDO ! neta
                ! 
              ENDIF ! if wq > acoustic
              ! 
            ENDDO ! imode
            !
          ENDIF ! endif  ekq in fsthick
          !
        ENDDO ! jbnd
        !
      ENDIF  ! endif  ekk in fsthick
      !
    ENDDO ! ibnd
    ! 
  ENDDO ! ik
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF ( iq .eq. nqtotf ) THEN
    !
#if defined(__MPI)
    !
    ! Note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_sum( epsilon2_abs, inter_pool_comm )
    CALL mp_sum( epsilon2_abs_lorenz, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
#endif
    ! 
    ! Output to stdout
    c = 'X'
    WRITE(c,"(i0)") neta
    format_string = "(5x,f15.6," // TRIM(c) // "E22.14)"
    ! 
    WRITE(stdout,'(5x,a)')
    WRITE(stdout,'(5x,a)') 'Phonon-assisted absorption versus energy'
    WRITE(stdout,'(5x,a,4f15.6)') 'Broadenings: ', eta(1:4)
    WRITE(stdout,'(5x,a,5f15.6)') '   ', eta(5:9)
    WRITE(stdout,'(5x,a)')
    WRITE(stdout,'(5x,a)') 'For the first Broadening we have:'
    ! For test-farm checking purposes, only show m=1
    WRITE(stdout,'(5x,a)') 'Photon energy (eV), Imaginary dielectric function along x,y,z'
    DO iw = 1, nomega
      WRITE(stdout,'(5x,f15.6,3E22.14)') omegap(iw)*ryd2ev, (epsilon2_abs(ipol,iw,1), ipol=1,3)
    ENDDO
    WRITE(stdout,'(5x,a)')
    WRITE(stdout,'(5x,a)') 'Values with other broadening are reported in the files epsilon2_indabs.dat'
    WRITE(stdout,'(5x,a)')
    ! 
    ! Output to file
    WRITE(c,"(i0)") neta+1
    format_string = "("//TRIM(c) // "E22.14)"

    OPEN(unit=iuindabs,file=nameF)
    WRITE(iuindabs,'(a)') '# Phonon-assisted absorption versus energy'
    WRITE(iuindabs,'(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
    DO iw = 1, nomega
      WRITE(iuindabs, format_string) omegap(iw)*ryd2ev, ( SUM(epsilon2_abs(:,iw,m))/3.0d0, m=1,neta)
    ENDDO
    CLOSE(iuindabs)
    ! 
    OPEN(unit=iuindabs,file='epsilon2_indabs_lorenz.dat')
    WRITE(iuindabs,'(a)') '# Phonon-assisted absorption versus energy'
    WRITE(iuindabs,'(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
    DO iw = 1, nomega
      WRITE(iuindabs, format_string) omegap(iw)*ryd2ev, ( SUM(epsilon2_abs_lorenz(:,iw,m))/3.0d0, m=1,neta )
    ENDDO
    CLOSE(iuindabs)
  ENDIF
  !
  END SUBROUTINE indabs
