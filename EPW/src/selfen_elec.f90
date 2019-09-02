  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
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
  USE io_epw,        ONLY : linewidth_elself
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, shortrange, &
                            fsthick, eptemp, ngaussw, degaussw, &
                            eps_acustic, efermi_read, fermi_energy,&
                            restart, restart_freq
  USE pwcom,         ONLY : ef
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, eta, nbndfst, &
                            nkf, epf17, wf, wqf, xkf, nkqtotf, adapt_smearing, &
                            sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew
  USE transportcom,  ONLY : lower_bnd
  USE control_flags, ONLY : iverbosity
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, eps8
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
  USE io_scattering, ONLY : electron_write
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
  !
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
  INTEGER :: nksqtotf
  !! Total number of k+q points 
  ! 
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
    WRITE(stdout,'(/5x,a)') REPEAT('=',67)
    WRITE(stdout,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")')
    WRITE(stdout,'(5x,a/)') REPEAT('=',67)
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
     WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  nksqtotf = nktotf ! odd-even for k,k+q
  !
  IF (restart) THEN
    ! Make everythin 0 except the range of k-points we are working on
    sigmar_all(:,1:lower_bnd-1) = zero
    sigmar_all(:,lower_bnd+nkf:nktotf) = zero
    sigmai_all(:,1:lower_bnd-1) = zero
    sigmai_all(:,lower_bnd+nkf:nktotf) = zero
    zi_all(:,1:lower_bnd-1) = zero
    zi_all(:,lower_bnd+nkf:nktotf) = zero
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
             !
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
                   weight = wqf(iq) * REAL (                                                   &
                           ((      wgkq + wgq(imode)) / (ekk - (ekq - wq(imode)) - ci * eta_tmp)  +  &
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
      IF (MOD(iqq,restart_freq) == 0) THEN
        WRITE(stdout, '(a,i10)' ) '     Creation of a restart point at ',iqq
        ! 
        CALL mp_sum(sigmar_all, inter_pool_comm)
        CALL mp_sum(sigmai_all, inter_pool_comm)
        CALL mp_sum(zi_all, inter_pool_comm)
        CALL mp_sum(fermicount, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
        CALL electron_write(iqq, totq, nksqtotf, sigmar_all, sigmai_all, zi_all)
        ! 
      ENDIF
    ENDIF 
  ENDIF ! in case of restart, do not do the first one
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF (iqq == totq) THEN
    !
    ALLOCATE(xkf_all(3,       nkqtotf))
    ALLOCATE(etf_all(nbndsub, nkqtotf))
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
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
      DO ibnd = 1, nbndfst
        ekk = etf_all(ibndmin - 1 + ibnd, ikk)
        n = 0
        tmp = 0.0_DP
        tmp2 = 0.0_DP
        tmp3 = 0.0_DP
        DO jbnd = 1, nbndfst
          ekk2 = etf_all (ibndmin-1+jbnd, ikk) 
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
    ENDDO ! nksqtotf
    !  
    ! Output electron SE here after looping over all q-points (with their contributions 
    ! summed in sigmar_all, etc.)
    !
    WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
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
      DO ik = 1, nksqtotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout,'(/5x,"ik = ",i7," coord.: ", 3f12.7)') ik, xkf_all(:,ikk)
        WRITE(stdout,'(5x,a)') REPEAT('-',67)
        !
        DO ibnd = 1, nbndfst
          !
          ! note that ekk does not depend on q 
          ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
          !
          ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
          zi_all (ibnd,ik) = one / ( one + zi_all (ibnd,ik) )
          !
          WRITE(stdout, 102) ibndmin-1+ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
                             ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik), one/zi_all(ibnd,ik)-one
!          WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
!                             ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik)
          IF (iverbosity == 3) THEN
            DO imode = 1, nmodes
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
        WRITE(stdout,'(5x,a/)') REPEAT('-',67)
        !
      ENDDO
    ENDIF
    !
    DO ibnd = 1, nbndfst
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
    DEALLOCATE(xkf_all)
    DEALLOCATE(etf_all)
    !
  ENDIF 
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE selfen_elec_q
  !-----------------------------------------------------------------------
