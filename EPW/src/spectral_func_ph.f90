  !                                                                            
  ! Copyright (C) 2010-2017 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_ph (iq )
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
  USE epwcom,    ONLY : nbndsub, lrepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, &
                        shortrange, nsmear, delta_smear, eps_acustic, &
                        efermi_read, fermi_energy, wmin_specfun,&
                        wmax_specfun, nw_specfun
  USE pwcom,     ONLY : nelec, ef, isk
  USE cell_base, ONLY : tpiba2
  USE elph2,     ONLY : epf17, ibndmax, ibndmin, etf, &
                        wkf, xqf, wqf, nkqf, nqtotf,   &
                        nkf, wf, a_all, &
                        dmef, gammai_all,gammar_all, efnew
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi, cone, ci
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm, ionode_id
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iq
  !! Current q-point index 
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
  INTEGER :: nrec
  !! Record index for reading the e-f matrix
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: ismear
  !! Smearing index
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
  REAL(kind=DP) :: wgq
  !! Bose occupation factor $n_{q\nu}(T)$
  REAL(kind=DP) :: wgkk
  !! Fermi-Dirac occupation factor $f_{nk}(T)$
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: degaussw0
  !! Gaussian smearing parameter
  REAL(kind=DP) :: inv_degaussw0
  !! Inverse Gaussian for efficiency reasons   
  REAL(kind=DP) :: eptemp0
  !! Temperature
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse temperature
  REAL(kind=DP) :: lambda_tot
  !! El-ph coupling strength
  REAL(kind=DP) :: lambda_tr_tot
  !! Transport el-ph coupling strength
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
  REAL(kind=DP) :: qsquared
  !! q-point squared
  REAL(kind=DP) :: eps0
  !! Dielectric function \varepsilon_\inf
  REAL(kind=DP) :: RTF
  !! Resta Thomas-Fermi
  REAL(kind=DP) :: qTF
  !! q Thomas-Fermi
  REAL(kind=DP) :: epsTF
  !! Thomas-Fermi dielectric function
  REAL(kind=DP), EXTERNAL :: efermig
  !! Function to compute the Fermi energy 
  REAL(kind=DP), external :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), external :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), external :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerence  
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  ! Thomas-Fermi screening according to Resta PRB 1977
  ! Here specific case of Diamond
  !eps0   = 5.7
  !RTF    = 2.76 
  !qTF    = 1.36 
  !qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2) * tpiba2
  !epsTF =  (qTF**2 + qsquared) / (qTF**2/eps0 * sin (sqrt(qsquared)*RTF)/(sqrt(qsquared)*RTF)+qsquared)
  !
  IF ( iq .eq. 1 ) THEN 
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Phonon Spectral Function Self-Energy in the Migdal Approximation (on the fly)")') 
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick.lt.1.d3 ) &
         WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
    !
    IF ( .not. ALLOCATED (gammai_all)  )  ALLOCATE( gammai_all (nmodes,nqtotf,nw_specfun) )
    IF ( .not. ALLOCATED (gammar_all)  )  ALLOCATE( gammar_all (nmodes,nqtotf,nw_specfun) )
    gammar_all(:,:,:)  = zero
    gammai_all(:,:,:)  = zero
    !
  ENDIF
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0   = (ismear-1) * delta_smear + eptemp
     ! 
     ! SP: Multiplication is faster than division ==> Important if called a lot
     !     in inner loops
     inv_degaussw0 = 1.0/degaussw0
     inv_eptemp0   = 1.0/eptemp0
     !
     ! Fermi level and corresponding DOS
     !
     IF ( efermi_read ) THEN
       !
       ef0 = fermi_energy
       !
     ELSE IF (nsmear > 1) THEN
       !
       ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw0,ngaussw,0,isk)
       ! if some bands are skipped (nbndskip.neq.0), nelec has already been
       ! recalculated 
       ! in ephwann_shuffle
       !
     ELSE !SP: This is added for efficiency reason because the efermig routine is slow
       ef0 = efnew
     ENDIF
     !
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     IF ( iq .eq. 1 ) THEN 
       WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
       WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
     ENDIF
     !
     CALL start_clock('PH SELF-ENERGY')
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
             wgkk = wgauss( -ekk*inv_eptemp0, -99)
             !w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
             !
             DO jbnd = 1, ibndmax-ibndmin+1
               !
               !  the fermi occupation for k+q
               ekq = etf (ibndmin-1+jbnd, ikq) - ef0
               wgkq = wgauss( -ekq*inv_eptemp0, -99)  
               !w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
               !
               ! here we take into account the zero-point sqrt(hbar/2M\omega)
               ! with hbar = 1 and M already contained in the eigenmodes
               ! g2 is Ry^2, wkf must already account for the spin factor
               !
               IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                  .OR. abs(xqf (3, iq))> eps2 )) THEN
                 ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                 !     number, in which case its square will be a negative number. 
                 g2 = (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp !* epsTF
               ELSE
                 g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp !* epsTF
               ENDIF
               !
               ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
               ! This is the imaginary part of the phonon self-energy, sans the matrix elements
               !
               !weight = wkf (ikk) * (wgkk - wgkq) * &
               !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
               !
               ! SP: The expression below (phonon self-energy) corresponds to 
               !  = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k -w_q] 
               !
               weight = wkf (ikk) * (wgkk - wgkq) * &
                  real ( cone / ( ekq - ekk + ci * degaussw0 )) 
               !
               gamma0  ( imode )  = gamma0   ( imode ) + weight * g2
               ! 
               DO iw = 1, nw_specfun
                 !
                 ww = wmin_specfun + dble (iw-1) * dw
                 !
                 weight = wkf (ikk) * (wgkk - wgkq) * &
                    real ( cone / ( ekq - ekk - ww + ci * degaussw0 )) 
                 gammar_all  (imode,iq,iw)  = gammar_all  (imode,iq,iw) + weight * g2
                 !
                 ! Normal implementation 
                 !weight = wkf (ikk) * (wgkk - wgkq) * &
                 !   aimag ( cone / ( ekq - ekk - ww + ci * degaussw0 ) ) 
                 !  
                 ! More stable:
                 ! Analytical im. part 
                 weight = pi * wkf (ikk) * (wgkk - wgkq) * &
                    w0gauss ( (ekq - ekk - ww) / degaussw0, 0) / degaussw0     
                 !
                 gammai_all (imode,iq,iw)  = gammai_all(imode,iq,iw)+ weight * g2 
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
     CALL stop_clock('PH SELF-ENERGY')
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
     IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun,nqtotf) )
     a_all(:,:) = zero
     !
     IF (iq == 1 ) THEN
       IF (mpime.eq.ionode_id) THEN
         OPEN(unit=iospectral,file='specfun.phon')
         OPEN(unit=iospectral_sup,file='specfun_sup.phon')
         WRITE(iospectral, '(/2x,a)') '#Phonon spectral function (meV)'
         WRITE(iospectral_sup, '(2x,a)') '#Phonon eigenenergies + real and im part of phonon self-energy (meV)'
         WRITE(iospectral, '(/2x,a)') '#K-point    Energy[meV]     A(q,w)[meV^-1]'
         WRITE(iospectral_sup, '(2x,a)') '#Q-point    Mode       w_q[eV]        w[eV]   &
&    Real Sigma(w)[meV]    Im Sigma(w=0)[meV]     Im Sigma(w)[meV]'
         WRITE(stdout,'(/5x,a)') 'Real and Imaginary part of the phonon self-energy (omega=0).'  
       ENDIF
     ENDIF
     !
     ! Write to output file  
     !WRITE(stdout,'(/5x,"iq = ",i7," coord.: ", 3f12.7)') iq, xqf(:,iq)
     DO imode = 1, nmodes
       wq = wf (imode, iq)
       ! Real and Im part of Phonon self-energy at 0 freq. 
       WRITE(stdout,105) imode, ryd2ev * wq, ryd2mev * gammar_all(imode,iq,1), ryd2mev * gammai_all(imode,iq,1)
     ENDDO 
     WRITE( stdout, '(5x,a,i8,a,i8)' ) &
      'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nqtotf

     !
     ! Write to support files
     DO iw = 1, nw_specfun
       !
       ww = wmin_specfun + dble (iw-1) * dw
       !
       DO imode = 1, nmodes
         ! 
         wq = wf (imode, iq)
         a_all(iw,iq) = a_all(iw,iq) + abs( gammai_all(imode,iq,iw) ) / pi / &
               ( ( ww - wq - gammar_all (imode,iq,iw) + gamma0 (imode))**two + (gammai_all(imode,iq,iw) )**two )
         !
         IF (mpime.eq.ionode_id) THEN
           WRITE(iospectral_sup,'(2i9,2x,f12.5,2x,f12.5,2x,E22.14,2x,E22.14,2x,E22.14)') iq,&
                imode, ryd2ev * wq, ryd2ev * ww, ryd2mev * gammar_all(imode,iq,iw), ryd2mev * gamma0(imode),&
                ryd2mev * gammai_all(imode,iq,iw)
         ENDIF
         !
       ENDDO 
       !
       IF (mpime.eq.ionode_id) THEN 
         WRITE(iospectral,'(2x,i7,2x,f12.5,2x,E22.14)') iq, ryd2ev * ww, a_all(iw,iq) / ryd2mev ! print to file 
       ENDIF
       !
     ENDDO
     !
     IF (iq == nqtotf ) THEN
       IF (mpime.eq.ionode_id) THEN
         CLOSE(iospectral)
         CLOSE(iospectral_sup)
       ENDIF 
     ENDIF   
     WRITE(stdout,'(5x,a/)') repeat('-',67)
     ! 
  ENDDO !smears
  !
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
103 FORMAT(5x,'iq = ',i7,'  w = ',f9.4,' eV   A(k,w) = ',e12.5,' meV^-1')
104 FORMAT(5x,'E( ',i3,' )=',f12.4,' meV   Re[Pi]=',f12.4,' - ', f12.4, ' meV Im[Pi]=',f12.4,' meV Im[Pi-o]=',f12.4,' meV ')
105 FORMAT(5x,'Omega( ',i3,' )=',f9.4,' eV   Re[Pi]=',f15.6,' meV Im[Pi]=',f15.6,' meV')
  !
  RETURN
  !
END SUBROUTINE spectral_func_ph 
