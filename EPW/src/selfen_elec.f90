  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_q ( iq )
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
  !!  02/06/2013 Modified by Roxana Margine 
  !!
  !!  This subroutine computes the contribution from phonon iq to all k-points
  !!  The outer loop in ephwann_shuffle.f90 will loop over all iq points
  !!  The contribution from each iq is summed at the end of this subroutine for iq=nqtotf 
  !!  to recover the per-ik electron self energy
  !!
  !!  RM 24/02/2014
  !!  redefined the size of sigmar_all, sigmai_all, and zi_all within the fermi windwow
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : linewidth_elself
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy,&
                            restart, restart_freq
  USE pwcom,         ONLY : ef !, nelec, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            nkf, epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew
  USE control_flags, ONLY : iverbosity
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iq
  !! Current q-point index 
  !
  ! Local variables 
  CHARACTER (len=256) :: nameF
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
  INTEGER :: lower_bnd
  !! Lower bounds index after k or q paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k or q paral
  INTEGER :: i
  !! Index for reading files
  ! 
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
  REAL(kind=DP) :: wq
  !! Phonon frequency on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgq
  !! Bose occupation factor $n_{q\nu}(T)$
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Self-energy factor 
  !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$ 
  !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$ 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
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
  !REAL(kind=DP), external :: efermig
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
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:)
  !! Collect k-point coordinate from all pools in parallel case
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Collect eigenenergies from all pools in parallel case
  !  
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  IF ( iq .eq. 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")')
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick .lt. 1.d3 ) &
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
           'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
     !
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  IF ( efermi_read ) THEN
    !
    ef0 = fermi_energy
    !
  ELSE
    !
    ef0 = efnew
    !ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
    ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
    ! in ephwann_shuffle
    !
  ENDIF
  !
  dosef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  dosef = dosef / two
  !
  IF ( iq .eq. 1 ) THEN 
     WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
     WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
     !WRITE (stdout, 101) dosef / ryd2ev, ef  * ryd2ev
     WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds( nksqtotf, lower_bnd, upper_bnd )
  ! 
  IF (restart) THEN
    ! Make everythin 0 except the range of k-points we are working on
    sigmar_all(:,1:lower_bnd-1) = zero
    sigmar_all(:,lower_bnd+nkf:nkqtotf/2) = zero
    sigmai_all(:,1:lower_bnd-1) = zero
    sigmai_all(:,lower_bnd+nkf:nkqtotf/2) = zero
    zi_all(:,1:lower_bnd-1) = zero
    zi_all(:,lower_bnd+nkf:nkqtotf/2) = zero
  ENDIF
  !
  IF ( iq .eq. 1 ) THEN 
     IF ( .not. ALLOCATED (sigmar_all) ) ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (sigmai_all) ) ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (zi_all) )     ALLOCATE( zi_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( iverbosity == 3 ) THEN
       IF ( .not. ALLOCATED (sigmai_mode) ) ALLOCATE( sigmai_mode(ibndmax-ibndmin+1, nmodes, nksqtotf) )
       sigmai_mode(:,:,:) = zero
     ENDIF
     sigmar_all(:,:) = zero
     sigmai_all(:,:) = zero
     zi_all(:,:) = zero
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
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
           !
           ! the phonon frequency and Bose occupation
           wq = wf (imode, iq)
           ! SP: Define the inverse for efficiency
           inv_wq = 1.0/( two * wq )
           wgq = wgauss( -wq*inv_eptemp0, -99)
           wgq = wgq / ( one - two * wgq )
           !
           ! SP: Avoid if statement in inner loops
           IF (wq .gt. eps_acustic) THEN
             g2_tmp = 1.0
           ELSE
             g2_tmp = 0.0
           ENDIF
           !
           DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  the energy of the electron at k (relative to Ef)
              ekk = etf (ibndmin-1+ibnd, ikk) - ef0
              !
              DO jbnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k+q
                 ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                 wgkq = wgauss( -ekq*inv_eptemp0, -99)  
                 !
                 ! here we take into account the zero-point sqrt(hbar/2M\omega)
                 ! with hbar = 1 and M already contained in the eigenmodes
                 ! g2 is Ry^2, wkf must already account for the spin factor
                 !
                 IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                    .OR. abs(xqf (3, iq))> eps2 )) THEN                         
                   ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                   !     number, in which case its square will be a negative number. 
                   g2 = (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp
                 ELSE
                   g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                 ENDIF        
                 !
                 ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                 ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                 ! (Many-Particle Physics, 3rd edition)
                 ! 
                 weight = wqf(iq) * real (                                                   &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
!                   ecutse needs to be defined if it's used 
!@                  if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 sigmar_all(ibnd,ik+lower_bnd-1) = sigmar_all(ibnd,ik+lower_bnd-1) + g2 * weight
                 !
                 ! Logical implementation
!                 weight = wqf(iq) * aimag (                                                  &
!                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
!                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 ! Delta implementation 
                 w0g1=w0gauss( (ekk-ekq+wq)/degaussw, 0) /degaussw
                 w0g2=w0gauss( (ekk-ekq-wq)/degaussw, 0) /degaussw
                 weight = pi * wqf(iq) * ( (wgkq+wgq)*w0g1 + (one-wgkq+wgq)*w0g2 )
                 !
                 sigmai_all(ibnd,ik+lower_bnd-1) = sigmai_all(ibnd,ik+lower_bnd-1) + g2 * weight
                 !
                 ! Mode-resolved
                 IF (iverbosity == 3) THEN
                   sigmai_mode(ibnd,imode,ik+lower_bnd-1) = sigmai_mode(ibnd,imode,ik+lower_bnd-1) + g2 * weight
                 ENDIF
                 !
                 ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                 !
                 weight = wqf(iq) * &
                         ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                           ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )  
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 zi_all(ibnd,ik+lower_bnd-1) = zi_all(ibnd,ik+lower_bnd-1) + g2 * weight
                 ! 
              ENDDO !jbnd
              !
           ENDDO !ibnd
           !
        ENDDO !imode
        !
     ENDIF ! endif  fsthick
     !
  ENDDO ! end loop on k
  !
  ! Creation of a restart point
  IF (restart) THEN
    IF (MOD(iq,restart_freq) == 0) THEN
      WRITE(stdout, '(a)' ) '     Creation of a restart point'
      ! 
      CALL mp_sum( sigmar_all, inter_pool_comm )
      CALL mp_sum( sigmai_all, inter_pool_comm )
      CALL mp_sum( zi_all, inter_pool_comm )
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      CALL electron_write(iq,nqtotf,nksqtotf,sigmar_all,sigmai_all,zi_all)
      ! 
    ENDIF
  ENDIF 
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF ( iq .eq. nqtotf ) THEN
     !
     ALLOCATE ( xkf_all      ( 3,       nkqtotf ), &
                etf_all      ( nbndsub, nkqtotf ) )
     xkf_all(:,:) = zero
     etf_all(:,:) = zero
     !
#if defined(__MPI)
     !
     ! note that poolgather2 works with the doubled grid (k and k+q)
     !
     CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
     CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
     CALL mp_sum( sigmar_all, inter_pool_comm )
     CALL mp_sum( sigmai_all, inter_pool_comm )
     IF (iverbosity == 3) CALL mp_sum( sigmai_mode, inter_pool_comm )
     CALL mp_sum( zi_all, inter_pool_comm )
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
     WRITE(stdout,'(5x,"Average over degenerate eigenstates is performed")')
     ! 
     DO ik = 1, nksqtotf
       ikk = 2 * ik - 1
       ikq = ikk + 1
       ! 
       DO ibnd = 1, ibndmax-ibndmin+1
         ekk = etf_all (ibndmin-1+ibnd, ikk)
         n = 0
         tmp = 0.0_DP
         tmp2 = 0.0_DP
         tmp3 = 0.0_DP
         !sigmar_tmp(:) = zero
         DO jbnd = 1, ibndmax-ibndmin+1
           ekk2 = etf_all (ibndmin-1+jbnd, ikk) 
           IF ( ABS(ekk2-ekk) < eps6 ) THEN
             n = n + 1
             tmp =  tmp + sigmar_all (jbnd,ik)
             tmp2 =  tmp2 + sigmai_all (jbnd,ik)
             tmp3 =  tmp3 + zi_all (jbnd,ik)
           ENDIF
           ! 
         ENDDO ! jbnd
         sigmar_tmp(ibnd) = tmp / float(n)
         sigmai_tmp(ibnd) = tmp2 / float(n)
         zi_tmp(ibnd) = tmp3 / float(n)
         !
       ENDDO ! ibnd
       sigmar_all (:,ik) = sigmar_tmp(:) 
       sigmai_all (:,ik) = sigmai_tmp(:)
       zi_all (:,ik)  = zi_tmp(:)
       ! 
     ENDDO ! nksqtotf
     !  
     ! Output electron SE here after looping over all q-points (with their contributions 
     ! summed in sigmar_all, etc.)
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     !
     IF (mpime.eq.ionode_id) THEN
       ! Write to file
       OPEN(unit=linewidth_elself,file='linewidth.elself')
       WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
       IF ( iverbosity == 3 ) THEN
         WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      imode          Im(Sgima)(meV)'
       ELSE
         WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sgima)(meV)'
       ENDIF
       ! 
       DO ik = 1, nksqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout,'(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:,ikk)
          WRITE(stdout,'(5x,a)') repeat('-',67)
          !
          DO ibnd = 1, ibndmax-ibndmin+1
            !
            ! note that ekk does not depend on q 
            ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
            !
            ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
            zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
            !
            WRITE(stdout, 102) ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
                               ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik), one/zi_all(ibnd,ik)-one
!            WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
!                               ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik)
            IF ( iverbosity == 3 ) THEN
              DO imode=1, nmodes
                WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
                WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
                WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
                WRITE(linewidth_elself,'(i9,2x)',advance='no') imode
                WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_mode(ibnd,imode,ik)
              ENDDO
            ELSE
              WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
              WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
              WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
              WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_all(ibnd,ik)
            ENDIF
            !
          ENDDO
          WRITE(stdout,'(5x,a/)') repeat('-',67)
          !
       ENDDO
     ENDIF
     !
     DO ibnd = 1, ibndmax-ibndmin+1
        !
        DO ik = 1, nksqtotf
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           ! note that ekk does not depend on q 
           ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
           !
           ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
           !zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
           !
           WRITE(stdout,'(2i9,5f12.4)') ik, ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd,ik), &
                                        ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik),  one/zi_all(ibnd,ik)-one
           ! 
        ENDDO
        !
        WRITE(stdout,'(a)') '  '
        !
     ENDDO
     !
     CLOSE(linewidth_elself)
     !
     IF ( ALLOCATED(xkf_all) )      DEALLOCATE( xkf_all )
     IF ( ALLOCATED(etf_all) )      DEALLOCATE( etf_all )
     IF ( ALLOCATED(sigmar_all) )   DEALLOCATE( sigmar_all )
     IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
     IF ( ALLOCATED(zi_all) )       DEALLOCATE( zi_all )
     IF ( ALLOCATED(sigmai_mode) )   DEALLOCATE( sigmai_mode )
     !
  ENDIF 
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_q
  !
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_k ( ik )
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
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : linewidth_elself
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy
  USE pwcom,         ONLY : ef !, nelec, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, etf_k, xqf, &
                            epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew, nqf
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6
  USE control_flags, ONLY : iverbosity
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
  !
  implicit none
  !
  CHARACTER (len=256) :: nameF
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
  INTEGER :: iq
  !! Counter on the q-point index 
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
  !! variables for collecting data from all pools in parallel case
  ! 
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
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, weight, dosef, &
                   w0g1, w0g2, inv_wq, inv_eptemp0, g2_tmp,&
                   inv_degaussw
  REAL(kind=DP), external :: efermig, dos_ef, wgauss, w0gauss, dos_ef_seq
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerence 
  !
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  IF ( ik .eq. 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")')
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick .lt. 1.d3 ) &
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
           'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     IF (mpime.eq.ionode_id) THEN
       OPEN(unit=linewidth_elself,file='linewidth.elself')
       WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
       IF ( iverbosity == 3 ) THEN
         WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      imode          Im(Sgima)(meV)'
       ELSE
         WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sgima)(meV)'
       ENDIF
     endif
     !
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  IF ( efermi_read ) THEN
     !
     ef0 = fermi_energy
     !
  ELSE
     !
     ef0 = efnew
     !ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
     ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
     ! in ephwann_shuffle
     !
  ENDIF
  !  
  IF (mpime .eq. ionode_id) THEN
    !   N(Ef) in the equation for lambda is the DOS per spin
    !dosef = dosef / two
    dosef = dos_ef_seq (ngaussw, degaussw, ef0, etf_k, wkf, nkqf, nbndsub)/2
    !
  ENDIF
  CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
  !
  !dosef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  !dosef = dosef / two
  !
  IF ( ik .eq. 1 ) THEN 
     WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
     WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
     !WRITE (stdout, 101) dosef / ryd2ev, ef  * ryd2ev
     WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  IF ( ik .eq. 1 ) THEN 
     IF ( .not. ALLOCATED (sigmar_all) ) ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (sigmai_all) ) ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (zi_all) )     ALLOCATE( zi_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( iverbosity == 3 ) THEN
       IF ( .not. ALLOCATED (sigmai_mode) ) ALLOCATE(sigmai_mode(ibndmax-ibndmin+1, nmodes, nksqtotf) )
       sigmai_mode(:,:,:) = zero
     ENDIF
     sigmar_all(:,:) = zero
     sigmai_all(:,:) = zero
     zi_all(:,:) = zero
  ENDIF
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
  DO iq = 1, nqf
     !
     ikq = 2 * iq
     ikk = ikq - 1
     !
     ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
     ! (but in this case they are the same)
     !
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
           !
           ! the phonon frequency and Bose occupation
           wq = wf (imode, iq)
           ! SP: Define the inverse for efficiency
           inv_wq = 1.0/( two * wq )
           wgq = wgauss( -wq*inv_eptemp0, -99)
           wgq = wgq / ( one - two * wgq )
           !
           ! SP: Avoid if statement in inner loops
           IF (wq .gt. eps_acustic) THEN
             g2_tmp = 1.0
           ELSE
             g2_tmp = 0.0
           ENDIF
           !
           DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  the energy of the electron at k (relative to Ef)
              ekk = etf (ibndmin-1+ibnd, ikk) - ef0
              !
              DO jbnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k+q
                 ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                 wgkq = wgauss( -ekq*inv_eptemp0, -99)  
                 !
                 ! here we take into account the zero-point sqrt(hbar/2M\omega)
                 ! with hbar = 1 and M already contained in the eigenmodes
                 ! g2 is Ry^2, wkf must already account for the spin factor
                 !
                 IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                    .OR. abs(xqf (3, iq))> eps2 )) THEN
                   ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                   !     number, in which case its square will be a negative number. 
                   g2 = (epf17 (jbnd, ibnd, imode, iq)**two)*inv_wq*g2_tmp
                 ELSE
                   g2 = (abs(epf17 (jbnd, ibnd, imode, iq))**two)*inv_wq*g2_tmp
                 ENDIF
                 !
                 ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                 ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                 ! (Many-Particle Physics, 3rd edition)
                 ! 
                 weight = wqf(iq) * real (                                                   &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
!                   ecutse needs to be defined if it's used 
!@                  if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 sigmar_all(ibnd,ik) = sigmar_all(ibnd,ik) + g2 * weight
                 !
                 ! Logical implementation
!                 weight = wqf(iq) * aimag (                                                  &
!                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
!                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 ! Delta implementation 
                 w0g1=w0gauss( (ekk-ekq+wq)/degaussw, 0) /degaussw
                 w0g2=w0gauss( (ekk-ekq-wq)/degaussw, 0) /degaussw
                 weight = pi * wqf(iq) * ( (wgkq+wgq)*w0g1 + (one-wgkq+wgq)*w0g2 )
                 !
                 sigmai_all(ibnd,ik) = sigmai_all(ibnd,ik) + g2 * weight
                 !
                 ! Mode-resolved
                 IF (iverbosity == 3) THEN
                   sigmai_mode(ibnd,imode,ik) = sigmai_mode(ibnd,imode,ik) + g2 * weight
                 ENDIF
                 !
                 ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                 !
                 weight = wqf(iq) * &
                         ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                           ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )  
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 zi_all(ibnd,ik) = zi_all(ibnd,ik) + g2 * weight
                 !
              ENDDO !jbnd
              !
           ENDDO !ibnd
           !
        ENDDO !imode
        !
     ENDIF ! endif  fsthick
     !
  ENDDO ! end loop on q
  !
  ! collect contributions from all pools (sum over k-points)
  ! this finishes the integral over the BZ  (k)
  !
  CALL mp_sum(sigmar_all,inter_pool_comm)
  CALL mp_sum(sigmai_all,inter_pool_comm)
  IF (iverbosity == 3) CALL mp_sum(sigmai_mode,inter_pool_comm)
  CALL mp_sum(zi_all,inter_pool_comm)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  ! Average over degenerate eigenstates:
  WRITE(stdout,'(5x,"Average over degenerate eigenstates is performed")')
  ! 
  ikk = 2 * ik - 1
  ikq = ikk + 1
  ! 
  DO ibnd = 1, ibndmax-ibndmin+1
    ekk = etf_k (ibndmin-1+ibnd, ikk)
    n = 0
    tmp = 0.0_DP
    tmp2 = 0.0_DP
    tmp3 = 0.0_DP
    DO jbnd = 1, ibndmax-ibndmin+1
      ekk2 = etf_k (ibndmin-1+jbnd, ikk)
      IF ( ABS(ekk2-ekk) < eps6 ) THEN
        n = n + 1
        tmp =  tmp + sigmar_all (jbnd,ik)
        tmp2 =  tmp2 + sigmai_all (jbnd,ik)
        tmp3 =  tmp3 + zi_all (jbnd,ik)
      ENDIF
      ! 
    ENDDO ! jbnd
    sigmar_tmp(ibnd) = tmp / float(n)
    sigmai_tmp(ibnd) = tmp2 / float(n)
    zi_tmp(ibnd) = tmp3 / float(n)
    !
  ENDDO ! ibnd
  sigmar_all (:,ik) = sigmar_tmp(:)
  sigmai_all (:,ik) = sigmai_tmp(:)
  zi_all (:,ik)  = zi_tmp(:)
  ! 
  ! ---
  !
  WRITE(stdout,'(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf(:,ikk)
  WRITE(stdout,'(5x,a)') repeat('-',67)
  !
  DO ibnd = 1, ibndmax-ibndmin+1
    !
    ! note that ekk does not depend on q 
    ekk = etf_k (ibndmin-1+ibnd, ikk) - ef0
    !
    ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
    zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
    !
    WRITE(stdout, 102) ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
                       ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik), one/zi_all(ibnd,ik)-one
!    WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
!                       ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik)

    IF (mpime.eq.ionode_id) THEN
      IF ( iverbosity == 3 ) THEN
        DO imode=1, nmodes
          WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
          WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
          WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
          WRITE(linewidth_elself,'(i9,2x)',advance='no') imode
          WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_mode(ibnd,imode,ik)
        ENDDO
      ELSE
        WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
        WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
        WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
        WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_all(ibnd,ik)
      ENDIF
    ENDIF
    !
  ENDDO
  WRITE(stdout,'(5x,a/)') repeat('-',67)
  !
  IF ( ik .eq. (nkqtotf - nqtotf)) THEN
    IF (mpime.eq.ionode_id) CLOSE(linewidth_elself)
    IF ( ALLOCATED(sigmar_all) )   DEALLOCATE( sigmar_all )
    IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
    IF ( ALLOCATED(zi_all) )       DEALLOCATE( zi_all )
    IF ( ALLOCATED(sigmai_mode) )   DEALLOCATE( sigmai_mode )
  ENDIF
  !
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  !
  RETURN
  !
END SUBROUTINE selfen_elec_k
