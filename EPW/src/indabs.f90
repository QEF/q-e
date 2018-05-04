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
  USE epwcom,        ONLY : nbndsub, lrepmatf, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy,&
                            vme, omegamin, omegamax, omegastep, n_r, scissor
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            nkf, epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, efnew, &
                            dmef, omegap, alpha_abs, vmef, etf_ks
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
  CHARACTER (len=256) :: nameF='indabs.dat'
  !! Name of the file
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
  !!-- polarization direction
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
  !!-- Occupation prefactors
  REAL(KIND=DP) :: pface
  !!-- Occupation prefactors
  REAL(KIND=DP) :: cfac
  !!- Absorption prefactor
  REAL(kind=DP), EXTERNAL :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), EXTERNAL :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), EXTERNAL :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerence  
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
  ! C = 4*pi^2*e^2 / (n_r c m_e^2) * 2 = 4 * pi^2 * 2 *2^2 *2 / (nr * 137*2) = 32/ (nr*137) = 2*1.15235180919/n_r
  ! 
  cfac = two*1.15235180919/n_r 
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
    IF ( .not. ALLOCATED (alpha_abs) ) ALLOCATE(alpha_abs(3,nomega))
    ! 
    alpha_abs = 0.d0
    DO iw = 1, nomega
      omegap(iw) = omegamin + (iw-1) * omegastep
    END DO
  END IF
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
      END IF
    END DO
    !
    ! RM - vme version should be checked
    IF ( vme ) THEN 
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          ! vmef is in units of Ryd * bohr
          IF (ABS(scissor) > eps6 .AND. &
              ABS( etf_ks(ibndmin-1+ibnd,ikk)-etf_ks(ibndmin-1+jbnd,ikk)) > eps6 .AND. &
              ABS( etf_ks(ibndmin-1+ibnd,ikq)-etf_ks(ibndmin-1+jbnd,ikq)) > eps6 ) THEN
            vkk(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk) &
                             * ( etf(ibndmin-1+ibnd,ikk)    - etf(ibndmin-1+jbnd,ikk) ) & 
                             / ( etf_ks(ibndmin-1+ibnd,ikk) - etf_ks(ibndmin-1+jbnd,ikk) )
            vkq(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd,ikq) &
                             * ( etf(ibndmin-1+ibnd,ikq)    - etf(ibndmin-1+jbnd,ikq) ) & 
                             / ( etf_ks(ibndmin-1+ibnd,ikq) - etf_ks(ibndmin-1+jbnd,ikq) )
          ELSE 
            vkk(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk)
            vkq(:,ibnd,jbnd) = vmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikq)
          END IF
        END DO
      END DO
    ELSE
      DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
          IF (ABS(scissor) > eps6 .AND. &
              ABS( etf_ks(ibndmin-1+ibnd,ikk)-etf_ks(ibndmin-1+jbnd,ikk)) > eps6 .AND. &
              ABS( etf_ks(ibndmin-1+ibnd,ikq)-etf_ks(ibndmin-1+jbnd,ikq)) > eps6 ) THEN
            vkk(:,ibnd,jbnd) = 2.0 * dmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk) &
                             * ( etf(ibndmin-1+ibnd,ikk)    - etf(ibndmin-1+jbnd,ikk) ) & 
                             / ( etf_ks(ibndmin-1+ibnd,ikk) - etf_ks(ibndmin-1+jbnd,ikk) )
            vkq(:,ibnd,jbnd) = 2.0 * dmef (:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikq) &
                             * ( etf(ibndmin-1+ibnd,ikq)    - etf(ibndmin-1+jbnd,ikq) ) & 
                             / ( etf_ks(ibndmin-1+ibnd,ikq) - etf_ks(ibndmin-1+jbnd,ikq) )
          ELSE
            vkk(:,ibnd,jbnd) = 2.0 * dmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikk) 
            vkq(:,ibnd,jbnd) = 2.0 * dmef(:, ibndmin-1+ibnd, ibndmin-1+jbnd, ikq) 
          END IF
        END DO
      END DO
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
                s1a = czero
                s1e = czero
                s2a = czero
                s2e = czero
                !
                DO mbnd = 1, ibndmax-ibndmin+1
                  !
                  ! The energy of the electron at k (relative to Ef)
                  ekmk = etf (ibndmin-1+mbnd, ikk) - ef0
                  ! The fermi occupation for k+q
                  ekmq = etf (ibndmin-1+mbnd, ikq) - ef0
                  !
                  s1a(:)  = s1a(:) + epf(mbnd, jbnd,imode) * 0.5 * vkk(:,ibnd, mbnd)  / &
                       (  ekmk  - ekq + wq(imode) + ci * degaussw )
                  s1e(:)  = s1e(:) + epf(mbnd, jbnd,imode) * 0.5 * vkk(:,ibnd, mbnd)  / &
                       (  ekmk  - ekq - wq(imode) + ci * degaussw )
                  s2a(:) =  s2a(:) + epf(ibnd, mbnd,imode) * 0.5 * vkq(:,mbnd, jbnd)   / &
                       (  ekmq  - ekk - wq(imode)+ ci * degaussw)
                  s2e(:) =  s2e(:) + epf(ibnd, mbnd,imode) * 0.5 * vkq(:,mbnd, jbnd)   / &
                       (  ekmq  - ekk + wq(imode)+ ci * degaussw)
                END DO
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
                    alpha_abs(ipol,iw) = alpha_abs(ipol,iw)  + &
                         (wkf(ikk)/2.0) * wqf(iq) * &
                         cfac / omegap(iw) * pfac  * weighta * abs( s1a(ipol) + s2a(ipol) )**2 / (2 * wq(imode) * omega )
                    alpha_abs(ipol,iw) = alpha_abs(ipol,iw)  + &
                         (wkf(ikk)/2.0) * wqf(iq) * &
                         cfac / omegap(iw) * pface * weighte * abs( s1e(ipol) + s2e(ipol) )**2 / (2 * wq(imode) * omega )
                  END DO ! ipol
                  !
                END DO ! iw
                !
              END IF ! if wq > acoustic
              ! 
            END DO ! imode
            !
          END IF ! endif  ekq in fsthick
          !
        END DO ! jbnd
        !
      END IF  ! endif  ekk in fsthick
      !
    END DO ! ibnd
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
    CALL mp_sum( alpha_abs, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
#endif
    ! 
    ! Output to stdout
    WRITE(stdout,'(5x,a)')
    WRITE(stdout,'(5x,a)') 'Phonon-assisted absorption coefficient versus energy'
    WRITE(stdout,'(5x,a)') 'Photon energy (eV), absorption coefficient (cm-1) along x,y,z'
    DO iw = 1, nomega
      WRITE(stdout, '(5x,f15.6,3E22.14,a)') omegap(iw)*ryd2ev, (alpha_abs(ipol,iw)/0.529177E-8,ipol=1,3), '  (cm-1)'
    ENDDO 
    ! 
    ! Output to file
    OPEN(unit=iuindabs,file=nameF)
    WRITE(iuindabs,'(a)') '# Phonon-assisted absorption coefficient versus energy'
    WRITE(iuindabs,'(a)') '# Photon energy (eV), absorption coefficient (cm-1) along x,y,z'
    DO iw = 1, nomega
      !WRITE(iuindabs, '(4f12.7)') omegap(iw)*ryd2ev, (alpha_abs(ipol,iw)/0.529177E-8,ipol=1,3)
      WRITE(iuindabs, '(4E22.14)') omegap(iw)*ryd2ev, (alpha_abs(ipol,iw)/0.529177E-8,ipol=1,3)
    END DO
    CLOSE(iuindabs)
  END IF
  ! 
  !RETURN
  !
  END SUBROUTINE indabs
