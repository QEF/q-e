  !                                                                            
  ! Copyright (C) 2010-2017 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_ph ( iqq, iq, totq )
  !-----------------------------------------------------------------------
  !
  !  Compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation.  This routine is similar to the one above
  !  but it is ONLY called from within ephwann_shuffle and calculates 
  !  the selfenergy for one phonon at a time.  Much smaller footprint on
  !  the disk
  !
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iospectral_sup, iospectral
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw, &
                        shortrange, nsmear, delta_smear, eps_acustic, &
                        efermi_read, fermi_energy, wmin_specfun,&
                        wmax_specfun, nw_specfun
  USE pwcom,     ONLY : nelec, ef, isk
  USE elph2,     ONLY : epf17, ibndmax, ibndmin, etf, &
                        wkf, xqf, nkqf, nkf, wf, a_all, efnew
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi, cone, ci, eps8
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm, ionode_id
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iqq
  !! Current q-point index from selecq
  INTEGER, INTENT (in) :: iq
  !! Current q-point index 
  INTEGER, INTENT (in) :: totq
  !! Total q-points in selecq window
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
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: iw 
  !! Counter on frequency for the phonon spectra
  !
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: wq
  !! Phonon frequency on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgkk
  !! Fermi-Dirac occupation factor $f_{nk}(T)$
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: inv_degaussw
  !! Inverse Gaussian for efficiency reasons   
  REAL(kind=DP) :: inv_eptemp
  !! Inverse temperature
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons  
  REAL(kind=DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: g2_tmp
    !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: gamma0(nmodes)
  !! Phonon self-energy
  REAL(kind=DP) :: ww
  !! Current frequency
  REAL(kind=DP) :: dw
  !! Frequency intervals 
  REAL(kind=DP), EXTERNAL :: efermig
  !! Function to compute the Fermi energy 
  REAL(kind=DP), external :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), external :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), external :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  REAL(kind=DP) :: gammai_all(nw_specfun, nmodes)
  !! Imaginary part of the frequency dependent spectral function
  REAL(kind=DP) :: gammar_all(nw_specfun, nmodes)
  !!  Real part of the Phonon self-energy (freq. dependent for spectral function)
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  gammar_all(:,:) = zero
  gammai_all(:,:) = zero
  !
  ! Thomas-Fermi screening according to Resta PRB 1977
  ! Here specific case of Diamond
  !eps0   = 5.7
  !RTF    = 2.76 
  !qTF    = 1.36 
  !qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2) * tpiba2
  !epsTF =  (qTF**2 + qsquared) / (qTF**2/eps0 * sin (sqrt(qsquared)*RTF)/(sqrt(qsquared)*RTF)+qsquared)
  !
  IF ( iqq == 1 ) THEN 
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Phonon Spectral Function Self-Energy in the Migdal Approximation (on the fly)")') 
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick.lt.1.d3 ) &
         WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
  ENDIF
  !
  !
  ! SP: Multiplication is faster than division ==> Important if called a lot
  !     in inner loops
  inv_degaussw = 1.0/degaussw
  inv_eptemp   = 1.0/eptemp
  !
  ! Fermi level and corresponding DOS
  !
  IF ( efermi_read ) THEN
    !
    ef0 = fermi_energy
    !
  ELSE IF (nsmear > 1) THEN
    !
    ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
    ! if some bands are skipped (nbndskip.neq.0), nelec has already been
    ! recalculated in ephwann_shuffle
    !
  ELSE !SP: This is added for efficiency reason because the efermig routine is slow
    ef0 = efnew
  ENDIF
  !
  dosef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  dosef = dosef / two
  !
  IF ( iqq .eq. 1 ) THEN 
    WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
    WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
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
    IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
         ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
      !
      fermicount = fermicount + 1
      !
      DO imode = 1, nmodes
        !
        ! the phonon frequency
        wq = wf (imode, iq)
        !
        ! SP : We should avoid branching statements (if statements) in
        !      innerloops. Therefore we do it here.
        inv_wq =  1.0/(two * wq)
        ! the coupling from Gamma acoustic phonons is negligible
        IF ( wq .gt. eps_acustic ) THEN
          g2_tmp = 1.0
        ELSE
          g2_tmp = 0.0
        ENDIF   
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          !  the fermi occupation for k
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          wgkk = wgauss( -ekk*inv_eptemp, -99)
          !w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
          !
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            !  the fermi occupation for k+q
            ekq = etf (ibndmin-1+jbnd, ikq) - ef0
            wgkq = wgauss( -ekq*inv_eptemp, -99)  
            !w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
            !
            ! here we take into account the zero-point sqrt(hbar/2M\omega)
            ! with hbar = 1 and M already contained in the eigenmodes
            ! g2 is Ry^2, wkf must already account for the spin factor
            !
            IF ( shortrange .AND. ( abs(xqf (1, iq))> eps8 .OR. abs(xqf (2, iq))> eps8 &
               .OR. abs(xqf (3, iq))> eps8 )) THEN
              ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
              !     number, in which case its square will be a negative number. 
              g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp ) !* epsTF
            ELSE
              g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp !* epsTF
            ENDIF
            !
            ! SP - 03/2019 - Retarded phonon self-energy
            !                See Eq. 145 of RMP 89, 015003 (2017)
            ! \Pi^R = k-point weight * [ [f(E_k+q) - f(E_k)]/ [E_k+q - E_k -w_q - id] 
            !                           -[f(E_k+q) - f(E_k)]/ [E_k+q - E_k - id] ]
            ! The second term is gamma0 (static) 
            !
            weight = wkf (ikk) * (wgkq - wgkk) * &
               REAL ( cone / ( ekq - ekk + ci * degaussw )) 
            !
            gamma0 ( imode ) = gamma0 ( imode ) + weight * g2
            ! 
            DO iw = 1, nw_specfun
              !
              ww = wmin_specfun + dble (iw-1) * dw
              !
              weight = wkf (ikk) * (wgkq - wgkk) * &
                 REAL ( cone / ( ekq - ekk - ww + ci * degaussw )) 
              gammar_all (iw, imode) = gammar_all (iw, imode) + weight * g2
              !
              ! Normal implementation 
              !weight = wkf (ikk) * (wgkq - wgkk) * &
              !   aimag ( cone / ( ekq - ekk - ww + ci * degaussw0 ) ) 
              !  
              ! More stable:
              ! Analytical im. part 
              weight = pi * wkf (ikk) * (wgkq - wgkk) * &
                 w0gauss ( (ekq - ekk - ww) / degaussw, 0) / degaussw     
              !
              gammai_all (iw, imode) = gammai_all(iw, imode) + weight * g2 
              ! 
            ENDDO
            !
          ENDDO ! jbnd
          !
        ENDDO   ! ibnd
        !
      ENDDO ! loop on q-modes
      !
    ENDIF ! endif fsthick
    !
  ENDDO ! loop on k
  !
  CALL stop_clock('PH SPECTRAL-FUNCTION')
  !
#if defined(__MPI)
  !
  ! collect contributions from all pools (sum over k-points)
  ! this finishes the integral over the BZ  (k)
  !
  CALL mp_sum(gammai_all,inter_pool_comm) 
  CALL mp_sum(gammar_all,inter_pool_comm) 
  CALL mp_sum(gamma0,inter_pool_comm) 
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
#endif
  !
  WRITE(stdout,'(5x,a)')
  IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun, totq) )
  a_all(:,iqq) = zero
  !
  IF (iqq == 1 ) THEN
    IF (mpime.eq.ionode_id) THEN
      OPEN(unit=iospectral,file='specfun.phon')
      OPEN(unit=iospectral_sup,file='specfun_sup.phon')
      WRITE(iospectral, '(/2x,a)') '#Phonon spectral function (meV)'
      WRITE(iospectral_sup, '(2x,a)') '#Phonon eigenenergies + real and im part of phonon self-energy (meV)'
      WRITE(iospectral, '(/2x,a)') '#Q-point    Energy[eV]     A(q,w)[meV^-1]'
      WRITE(iospectral_sup, '(2x,a)') '#Q-point    Mode       w_q[eV]        w[eV]   &
& Real Sigma(w)[meV]   Real Sigma(w=0)[meV]     Im Sigma(w)[meV]'
    ENDIF
  ENDIF
  !
  ! Write to output file  
  !WRITE(stdout,'(/5x,"iq = ",i7," coord.: ", 3f12.7)') iq, xqf(:,iq)
  WRITE(stdout,'(/5x,a)') 'Real and Imaginary part of the phonon self-energy (omega=0) without gamma0.'  
  DO imode = 1, nmodes
    wq = wf (imode, iq)
    ! Real and Im part of Phonon self-energy at 0 freq. 
    WRITE(stdout,105) imode, ryd2ev * wq, ryd2mev * gammar_all(1,imode), ryd2mev * gammai_all(1,imode)
  ENDDO 
  WRITE( stdout, '(5x,a,i8,a,i8)' ) &
   'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', totq
  !
  ! Write to support files
  DO iw = 1, nw_specfun
    !
    ww = wmin_specfun + dble (iw-1) * dw
    !
    DO imode = 1, nmodes
      ! 
      wq = wf (imode, iq)
      !a_all(iw,iq) = a_all(iw,iq) + abs( gammai_all(imode,iq,iw) ) / pi / &
      !      ( ( ww - wq - gammar_all (imode,iq,iw) + gamma0 (imode))**two + (gammai_all(imode,iq,iw) )**two )
      ! SP: From Eq. 16 of PRB 9, 4733 (1974)
      !    Also in Eq.2 of PRL 119, 017001 (2017). 
      a_all(iw,iqq) = a_all(iw,iqq) + (1.0d0/pi) * ((2*wq)**2) * ABS( gammai_all(iw, imode) ) / &
            ( ( ww**2 - wq**2 - 2 * wq * ( gammar_all (iw, imode) - gamma0 (imode) ) )**two +&
              (2 * wq * gammai_all(iw, imode) )**two )
      !
      IF (mpime.eq.ionode_id) THEN
        WRITE(iospectral_sup,'(2i9,2x,f12.5,2x,f12.5,2x,E22.14,2x,E22.14,2x,E22.14)') iq,&
             imode, ryd2ev * wq, ryd2ev * ww, ryd2mev * gammar_all(iw, imode), ryd2mev * gamma0(imode),&
             ryd2mev * gammai_all(iw, imode)
      ENDIF
      !
    ENDDO 
    !
    IF (mpime.eq.ionode_id) THEN 
      WRITE(iospectral,'(2x,i7,2x,f12.5,2x,E22.14)') iq, ryd2ev * ww, a_all(iw,iqq) / ryd2mev ! print to file 
    ENDIF
    !
  ENDDO
  !
  IF (iqq == totq ) THEN
    IF (mpime == ionode_id) THEN
      CLOSE(iospectral)
      CLOSE(iospectral_sup)
    ENDIF 
  ENDIF   
  WRITE(stdout,'(5x,a/)') repeat('-',67)
  ! 
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
105 FORMAT(5x,'Omega( ',i3,' )=',f9.4,' eV   Re[Pi]=',f15.6,' meV Im[Pi]=',f15.6,' meV')
  !
  RETURN
  !
END SUBROUTINE spectral_func_ph 
