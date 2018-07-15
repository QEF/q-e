  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2018 Samuel Ponce'
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE transport
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutine linked with electronic transport  
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE scattering_rate_q( iq, ef0, efcb, first_cycle ) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine computes the scattering rate (inv_tau)
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, fsthick, eps_acustic, degaussw, & 
                              nstemp, scattering_serta, scattering_0rta, shortrange,&
                              restart, restart_freq, restart_filq, vme
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, dmef, vmef, wf, wqf, & 
                              epf17, nqtotf, nkqtotf, inv_tau_all, inv_tau_allcb, &
                              xqf, zi_allvb, zi_allcb
    USE transportcom,  ONLY : transp_temp, lower_bnd
    USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, & 
                              eps6
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : world_comm
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT (INOUT) :: first_cycle
    !! Use to determine weather this is the first cycle after restart 
    INTEGER, INTENT(IN) :: iq
    !! Q-point inde
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp. Could be unused (0).
    !
    ! Local variables
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates  
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: imode
    !! Local mode index
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: nqtotf_new
    !! Number of q-point in the new dataset
    !
    REAL(kind=DP) :: tmp
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(kind=DP) :: tmp2
    !! Temporary variable for zi_all
    REAL(kind=DP) :: ekk2
    !! Temporary variable to the eigenenergies for the degenerate average  
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND=DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND=DP) :: g2
    !! Electron-phonon matrix elements squared (g2 is Ry^2) 
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: w0g1
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$ 
    REAL(KIND=DP) :: w0g2 
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
    REAL(KIND=DP) :: inv_wq 
    !! Inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND=DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
    REAL(KIND=DP) :: g2_tmp 
    !! Used to set component to 0 if the phonon freq. is too low. This is defined
    !! for efficiency reasons as if statement should be avoided in inner-most loops.
    REAL(KIND=DP) :: inv_degaussw
    !! 1.0/degaussw. Defined for efficiency reasons. 
    REAL(KIND=DP) :: wq
    !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.  
    REAL(KIND=DP) :: wgq
    !! Bose-Einstein occupation function $$n_{q\nu}$$
    REAL(kind=DP) :: weight
    !! Self-energy factor 
    REAL(KIND=DP) :: fmkq
    !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
    REAL(KIND=DP) :: trans_prob
    !! Transition probability function
    REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(KIND=DP) :: vkq(3,ibndmax-ibndmin+1)
    !! Electronic velocity $$v_{m\mathbf{k+q}}$$
    REAL(KIND=DP) :: vel_factor(ibndmax-ibndmin+1,ibndmax-ibndmin+1)
    !! Velocity factor  $$ 1 - \frac{(v_{nk} \cdot v_{mk+q})}{ |v_{nk}|^2} $$
    REAL(kind=DP) :: inv_tau_tmp(ibndmax-ibndmin+1)
    !! Temporary array to store the scattering rates
    REAL(kind=DP) :: zi_tmp(ibndmax-ibndmin+1)
    !! Temporary array to store the zi
    REAL(KIND=DP), ALLOCATABLE :: inv_tau_all_new (:,:,:)
    !! New scattering rates to be merged
    !
    REAL(KIND=DP), ALLOCATABLE :: etf_all(:,:)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND=DP), EXTERNAL :: DDOT
    !! Dot product function
    REAL(KIND=DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function  
    REAL(kind=DP), PARAMETER :: eps = 1.d-4
    !! Tolerence parameter for the velocity
    REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
    !! Tolerence
    ! 
    CALL start_clock ( 'SCAT' )
    ! 
    IF ( iq .eq. 1 ) THEN
      !
      WRITE(stdout,'(/5x,a)') repeat('=',67)
      WRITE(stdout,'(5x,"Scattering rate")')
      WRITE(stdout,'(5x,a/)') repeat('=',67)
      !
      IF ( fsthick .lt. 1.d3 ) &
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ',(ef-fsthick) * ryd2ev, ' eV and ',&
                (ef+fsthick) * ryd2ev, ' eV will be included'
        WRITE(stdout,'(5x,a/)')
      !
      !IF ( .not. ALLOCATED (inv_tau_all) ) ALLOCATE( inv_tau_all(nstemp,ibndmax-ibndmin+1,nkqtotf/2) )
      !inv_tau_all(:,:,:) = zero
      !
    ENDIF
    ! 
    ! In the case of a restart do not add the first step
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ELSE
      ! loop over temperatures
      DO itemp = 1, nstemp
        !
        etemp = transp_temp(itemp)
        !
        ! SP: Define the inverse so that we can efficiently multiply instead of
        ! dividing
        !
        inv_etemp = 1.0/etemp
        inv_degaussw = 1.0/degaussw
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF ( scattering_0rta ) THEN 
            !vel_factor = 1 - (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
            vel_factor(:,:) = zero
            IF ( vme ) THEN 
              DO ibnd = 1, ibndmax-ibndmin+1
                !
                ! vkk(3,nbnd) - velocity for k
                vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                !
                DO jbnd = 1, ibndmax-ibndmin+1
                  !
                  ! vkq(3,nbnd) - velocity for k + q
                  vkq(:,jbnd) = REAL (vmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
                  !
                  IF ( abs( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2 ) > eps) &
                    vel_factor(ibnd,jbnd) = DDOT(3, vkk(:,ibnd), 1, vkq(:,jbnd), 1) / &
                                            DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd), 1)
                ENDDO
              ENDDO
            ELSE
              DO ibnd = 1, ibndmax-ibndmin+1
                !
                ! vkk(3,nbnd) - velocity for k
                vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                !
                DO jbnd = 1, ibndmax-ibndmin+1
                  ! 
                  ! vkq(3,nbnd) - velocity for k + q
                  vkq(:,jbnd) = 2.0 * REAL (dmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
                  !
                  IF ( abs( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2 ) > eps) &
                    vel_factor(ibnd,jbnd) = DDOT(3, vkk(:,ibnd), 1, vkq(:,jbnd), 1) / &
                                            DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd), 1)
                ENDDO  
              ENDDO
            ENDIF
            vel_factor(:,:) = one - vel_factor(:,:)
          ENDIF
          !
          ! We are not consistent with ef from ephwann_shuffle but it should not 
          ! matter if fstick is large enough.
          !IF ( ( minval ( abs(etf (:, ikk) - ef0(itemp)) ) .lt. fsthick ) .AND. &
          !     ( minval ( abs(etf (:, ikq) - ef0(itemp)) ) .lt. fsthick ) ) THEN
          ! If scissor = 0 then 
          IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
               ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
            !
            DO imode = 1, nmodes
              !
              ! the phonon frequency and bose occupation
              wq = wf (imode, iq)
              !
              ! SP : Avoid if statement in inner loops
              ! the coupling from Gamma acoustic phonons is negligible
              IF ( wq .gt. eps_acustic ) THEN
                g2_tmp = 1.0
                wgq = wgauss( -wq*inv_etemp, -99)
                wgq = wgq / ( one - two * wgq )
                ! SP : Define the inverse for efficiency
                inv_wq =  1.0/(two * wq) 
              ELSE
                g2_tmp = 0.0
                wgq = 0.0
                inv_wq = 0.0
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                !
                !  energy at k (relative to Ef)
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                !
                DO jbnd = 1, ibndmax-ibndmin+1
                  !
                  !  energy and fermi occupation at k+q
                  ekq = etf (ibndmin-1+jbnd, ikq) - ef0(itemp)
                  fmkq = wgauss( -ekq*inv_etemp, -99)
                  !
                  ! here we take into account the zero-point sqrt(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                  ! need to treat it like the normal g with abs(g).
                  IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                     .OR. abs(xqf (3, iq))> eps2 )) THEN
                    ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                    !     number, in which case its square will be a negative number. 
                    g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP )
                  ELSE
                    g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                  ENDIF
                  !
                  ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                  w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                  w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                  !
                  ! transition probability 
                  ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                  ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                  !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                  !
                  ! DBSP Just to try
                  !trans_prob = pi *  g2 * & 
                  !             ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                  trans_prob = pi * wqf(iq) * g2 * & 
                               ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                  !
                  !if ((ik == 6) .and. (ibnd == 4) .and. jbnd==1 .and. imode==6) then
                  !  print*,'wqf(iq) ',wqf(iq)
                  !  print*,'wq ',wq
                  !  print*,'inv_etemp ',inv_etemp
                  !  print*,'g2 ',g2
                  !  print*,'fmkq ',fmkq
                  !  print*,'wgq ',wgq
                  !  print*,'w0g1 ',w0g1
                  !  print*,'trans_prob ik, ibnd, jbnd, imode ', trans_prob, ik, ibnd, jbnd, imode
                  !  print*,'inv_tau',SUM(inv_tau_all(1,5:8,27))
                  !end if
                  !  
                  IF ( scattering_serta ) THEN 
                    ! energy relaxation time approximation 
                    inv_tau_all(itemp,ibnd,ik+lower_bnd-1) = inv_tau_all(itemp,ibnd,ik+lower_bnd-1) + two * trans_prob
                    
                  ELSEIF ( scattering_0rta ) THEN 
                    ! momentum relaxation time approximation
                    inv_tau_all(itemp,ibnd,ik+lower_bnd-1) = inv_tau_all(itemp,ibnd,ik+lower_bnd-1) &
                                           + two * trans_prob * vel_factor(ibnd,jbnd)
                  ENDIF
                  !
                  ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                  !
                  weight = wqf(iq) * &
                          ( (       fmkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                   ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                            ( one - fmkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                   ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )
                  !
                  zi_allvb(itemp,ibnd,ik+lower_bnd-1) = zi_allvb(itemp,ibnd,ik+lower_bnd-1) + g2 * weight
                  ! 
                ENDDO !jbnd
                !
              ENDDO !ibnd
              !
              ! In this case we are also computing the scattering rate for another Fermi level position
              ! This is used to compute both the electron and hole mobility at the same time.  
              IF ( ABS(efcb(itemp)) > eps ) THEN
                ! 
                DO ibnd = 1, ibndmax-ibndmin+1
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                  !
                  DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  energy and fermi occupation at k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - efcb(itemp)
                    fmkq = wgauss( -ekq*inv_etemp, -99)
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                    ! need to treat it like the normal g with abs(g).
                    IF ( shortrange .AND. ( abs(xqf (1, iq))> eps2 .OR. abs(xqf (2, iq))> eps2 &
                       .OR. abs(xqf (3, iq))> eps2 )) THEN
                      ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                      !     number, in which case its square will be a negative number. 
                      g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp, KIND=DP)
                    ELSE
                      g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                    ENDIF
                    !
                    ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                    w0g1 = w0gauss( (ekk-ekq+wq) * inv_degaussw, 0) * inv_degaussw
                    w0g2 = w0gauss( (ekk-ekq-wq) * inv_degaussw, 0) * inv_degaussw
                    !
                    ! transition probability 
                    ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                    ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                    !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                    !
                    trans_prob = pi * wqf(iq) * g2 * &
                                 ( (fmkq+wgq)*w0g1 + (one-fmkq+wgq)*w0g2 )
                    !
                    IF ( scattering_serta ) THEN
                      ! energy relaxation time approximation 
                      inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) = inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) + two * trans_prob
                      !
                    ELSEIF ( scattering_0rta ) THEN
                      ! momentum relaxation time approximation
                      inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) = inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) &
                                             + two * trans_prob * vel_factor(ibnd,jbnd)
                    ENDIF
                    !
                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                    !
                    weight = wqf(iq) * &
                            ( (       fmkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                     ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                              ( one - fmkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                     ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )
                    !
                    zi_allcb(itemp,ibnd,ik+lower_bnd-1) = zi_allcb(itemp,ibnd,ik+lower_bnd-1) + g2 * weight
                    ! 
                  ENDDO !jbnd
                  !
                ENDDO !ibnd
                ! 
              ENDIF ! ABS(efcb) < eps
              !
            ENDDO !imode
            !
          ENDIF ! endif  fsthick
          !
        ENDDO ! end loop on k
      ENDDO ! itemp
      !
      ! Creation of a restart point
      IF (restart) THEN
        IF (MOD(iq,restart_freq) == 0 ) THEN
          WRITE(stdout, '(a)' ) '     Creation of a restart point'
          ! 
          ! The mp_sum will aggreage the results on each k-points. 
          CALL mp_sum( inv_tau_all, world_comm )
          CALL mp_sum( zi_allvb,    world_comm )
          !
          IF ( ABS(efcb(1)) > eps ) THEN
            ! 
            CALL mp_sum( inv_tau_allcb, world_comm ) 
            CALL mp_sum( zi_allcb,      world_comm ) 
            ! 
          ENDIF
          ! 
          IF ( ABS(efcb(1)) > eps ) THEN
            CALL tau_write(iq,nqtotf,nkqtotf/2,.TRUE.)
          ELSE
            CALL tau_write(iq,nqtotf,nkqtotf/2,.FALSE.)
          ENDIF
          ! 
          ! Now show intermediate mobility with that amount of q-points
          CALL transport_coeffs(ef0,efcb)
          ! 
        ENDIF
      ENDIF
      ! 
    ENDIF ! first_cycle
    ! 
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF ( iq .eq. nqtotf ) THEN
      !
      ! The total number of k points
      !
      ALLOCATE ( etf_all ( nbndsub, nkqtotf ))
      !
      CALL mp_sum( inv_tau_all, world_comm )
      IF (ABS(efcb(1)) > eps) CALL mp_sum( inv_tau_allcb, world_comm )
      CALL mp_sum( zi_allvb, world_comm )
      IF (ABS(efcb(1)) > eps) CALL mp_sum( zi_allcb, world_comm )
      !print*,'zi_allvb SUM ',SUM(zi_allvb)
      !print*,'inv_tau_all SUM ',SUM(inv_tau_all)
      !
#ifdef __MPI
      !
      ! collect contributions from all pools (sum over k-points)
      ! this finishes the integral over the BZ (k)
      !
      CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf, etf_all )
#else
      !
      etf_all = etf
#endif
      !
      DO itemp = 1, nstemp  
        ! 
        etemp = transp_temp(itemp)
        WRITE(stdout, '(a,f8.3,a)' ) '     Temperature ',etemp * ryd2ev / kelvin2eV,' K'
        !
        ! In case we read another q-file, merge the scattering here
        IF (restart_filq .ne. '') THEN
          ! 
          ALLOCATE( inv_tau_all_new(nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
          inv_tau_all_new(:,:,:) = zero
          ! 
          CALL merge_read( nkqtotf/2, nqtotf_new, inv_tau_all_new ) 
          ! 
          inv_tau_all(:,:,:) = ( inv_tau_all(:,:,:) * nqtotf &
                              + inv_tau_all_new(:,:,:) * nqtotf_new ) / (nqtotf+nqtotf_new)
          !
          WRITE(stdout, '(a)' ) '     '
          WRITE(stdout, '(a,i10,a)' ) '     Merge scattering for a total of ',nqtotf+nqtotf_new,' q-points'
          ! 
          CALL tau_write(iq+nqtotf_new,nqtotf+nqtotf_new,nkqtotf/2)
          WRITE(stdout, '(a)' ) '     Write to restart file the sum'
          WRITE(stdout, '(a)' ) '     '
          !
          ! 
        ENDIF
        ! Average over degenerate eigenstates:
        WRITE(stdout,'(5x,"Average over degenerate eigenstates is performed")')
        ! 
        DO ik = 1, nkqtotf/2
          ikk = 2 * ik - 1
          ikq = ikk + 1
          ! 
          DO ibnd = 1, ibndmax-ibndmin+1
            ekk = etf_all (ibndmin-1+ibnd, ikk)
            n = 0
            tmp = 0.0_DP
            tmp2 = 0.0_DP
            DO jbnd = 1, ibndmax-ibndmin+1
              ekk2 = etf_all (ibndmin-1+jbnd, ikk)
              IF ( ABS(ekk2-ekk) < eps6 ) THEN
                n = n + 1
                tmp =  tmp + inv_tau_all (itemp,jbnd,ik)
                tmp2 =  tmp2 + zi_allvb(itemp,jbnd,ik)
              ENDIF
              ! 
            ENDDO ! jbnd
            inv_tau_tmp(ibnd) = tmp / float(n)
            zi_tmp(ibnd) = tmp2 / float(n)
            !
          ENDDO ! ibnd
          inv_tau_all (itemp,:,ik) = inv_tau_tmp(:)
          zi_allvb (itemp,:,ik) = zi_tmp(:)
          ! 
        ENDDO ! nkqtotf
        !
        IF (ABS(efcb(itemp)) > eps) THEN 
          ! Average over degenerate eigenstates:
          WRITE(stdout,'(5x,"Average over degenerate eigenstates in CB is performed")')
          ! 
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1 
            ikq = ikk + 1 
            ! 
            DO ibnd = 1, ibndmax-ibndmin+1
              ekk = etf_all (ibndmin-1+ibnd, ikk)
              n = 0 
              tmp = 0.0_DP
              tmp2 = 0.0_DP
              DO jbnd = 1, ibndmax-ibndmin+1
                ekk2 = etf_all (ibndmin-1+jbnd, ikk)
                IF ( ABS(ekk2-ekk) < eps6 ) THEN
                  n = n + 1 
                  tmp =  tmp + inv_tau_allcb (itemp,jbnd,ik)
                  tmp2 =  tmp2 + zi_allcb (itemp,jbnd,ik)
                ENDIF
                ! 
              ENDDO ! jbnd
              inv_tau_tmp(ibnd) = tmp / float(n)
              zi_tmp(ibnd) = tmp2 / float(n)
              !
            ENDDO ! ibnd
            inv_tau_allcb (itemp,:,ik) = inv_tau_tmp(:)
            zi_allcb (itemp,:,ik) = zi_tmp(:)
            ! 
          ENDDO ! nkqtotf
        ENDIF
        ! 
        ! Output scattering rates here after looping over all q-points
        ! (with their contributions summed in inv_tau_all, etc.)
        !print*,'inv_tau_all(1,1,1) ',inv_tau_all(1,1,1)
        CALL scattering_write(itemp, etemp, ef0, etf_all)
        !
      ENDDO !nstemp 
      !
      IF ( ALLOCATED(etf_all) )     DEALLOCATE( etf_all )
      ! 
      ! Creation of a restart point at the end
      IF (restart) THEN
        WRITE(stdout, '(a)' ) '     Creation of the final restart point'
        ! 
        IF ( ABS(efcb(1)) > eps ) THEN
          CALL tau_write(iq,nqtotf,nkqtotf/2,.TRUE.)
        ELSE
          CALL tau_write(iq,nqtotf,nkqtotf/2,.FALSE.)
        ENDIF
        ! 
      ENDIF ! restart
      ! 
    ENDIF ! iq 
    !
    CALL stop_clock ( 'SCAT' )
    ! DBSP
    !write(stdout,*),'iq ',iq
    !print*,shape(inv_tau_all)
    !write(stdout,*),'inv_tau_all(1,5:8,21) ',SUM(inv_tau_all(3,5:8,1))
    !write(stdout,*),'inv_tau_all(1,5:8,:) ',SUM(inv_tau_all(3,5:8,:))
    !write(stdout,*),'SUM(inv_tau_all) ',SUM(inv_tau_all(3,:,:))
    !write(stdout,*),'first_cycle ',first_cycle
    !
    RETURN
    !
    END SUBROUTINE scattering_rate_q
    !-----------------------------------------------------------------------
    !       
    !-----------------------------------------------------------------------
    SUBROUTINE transport_coeffs (ef0,efcb)
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine computes the transport coefficients
    !!  SP - June 2018 - Update for symmetries in velocities when using homogeneous grids. 
    !!       This is currently commented out since it is ONLY needed if we are
    !!       interested in the off-diagonal mobility_\alpha\beta terms. 
    !!       At the moment we just want mobility_\alpha\alpha so it makes no difference.
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, meta_ionode_id
    USE cell_base, ONLY : alat, at, omega
    USE io_files,  ONLY : prefix 
    USE io_epw,    ONLY : iufilsigma 
    USE epwcom,    ONLY : nbndsub, fsthick, & 
                          system_2d, nstemp, &
                          int_mob, ncarrier, scatread, &
                          iterative_bte, vme
    USE pwcom,     ONLY : ef 
    USE elph2,     ONLY : ibndmax, ibndmin, etf, nkf, wkf, dmef, vmef, & 
                          inv_tau_all, nkqtotf, Fi_all, inv_tau_allcb, &
                          zi_allvb, zi_allcb, Fi_allcb, BZtoIBZ,       &
                          s_BZtoIBZ, map_rebal
    USE transportcom,  ONLY : transp_temp
    USE constants_epw, ONLY : zero, one, bohr2ang, ryd2ev, electron_SI, &
                              kelvin2eV, hbar, Ang2m, hbarJ, ang2cm, czero
    USE mp,        ONLY : mp_sum
    USE mp_global, ONLY : world_comm
    USE mp_world,  ONLY : mpime
    ! SP - Uncomment to use symmetries on velocities
    USE symm_base,     ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE io_global,     ONLY : ionode_id
    USE cell_base, ONLY : bg
    USE mp,        ONLY : mp_bcast
    USE mp_global, ONLY : inter_pool_comm
    USE epwcom,    ONLY : mp_mesh_k, nkf1, nkf2, nkf3
    USE constants_epw, ONLY : eps6
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    ! 
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND=DP), INTENT(IN) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp (could be 0)
    !
    ! Local variables
    INTEGER :: i
    !! Cartesian direction index 
    INTEGER :: j
    !! Cartesian direction index 
    INTEGER :: ij
    !! Cartesian coupled index for matrix. 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    !  SP - Uncomment to use symmetries on velocities
    INTEGER :: nkqtotf_tmp
    !! Temporary k-q points.
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: nb
    !! Number of points in the BZ corresponding to a point in IBZ    
    INTEGER :: BZtoIBZ_tmp(nkf1*nkf2*nkf3)
    !! Temporary mapping
    ! 
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
   ! REAL(KIND=DP) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND=DP) :: dfnk
    !! Derivative Fermi distribution $$-df_{nk}/dE_{nk}$$
    REAL(KIND=DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND=DP) :: tau 
    !! Relaxation time
    REAL(KIND=DP) :: conv_factor1
    !! Conversion factor for the conductivity 
    REAL(KIND=DP) :: inv_cell 
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND=DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND=DP) :: carrier_density_prt
    !! Carrier density [nb of carrier per unit cell] in cm^-3 unit
    REAL(KIND=DP) :: fnk
    !! Fermi-Dirac occupation function 
    REAL(KIND=DP) :: mobility
    !! Sum of the diagonalized mobilities [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND=DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs] 
    REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
    !! Electron velocity vector for a band. 
    REAL(KIND=DP) :: Sigma(9,nstemp)
    !! Conductivity matrix in vector form
    REAL(KIND=DP) :: SigmaZ(9,nstemp)
    !! Conductivity matrix in vector form with Znk
    REAL(KIND=DP) :: Sigma_m(3,3,nstemp)
    !! Conductivity matrix
    REAL(KIND=DP) :: sigma_up(3,3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND=DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND=DP) :: sigma_vect(3,3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND=DP) :: Znk
    !! Real Znk from \lambda_nk (called zi_allvb or zi_allcb)
    REAL(KIND=DP) :: tdf_sigma(9)
    !! Temporary file
    REAL(kind=DP), PARAMETER :: eps = 1.d-4
    !! Tolerence
    REAL(KIND=DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND=DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND=DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    CHARACTER (len=256) :: filsigma
    !! File for the conductivity  
    REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    COMPLEX(kind=DP), ALLOCATABLE :: dmef_all(:,:,:,:)
    !! dipole matrix elements on the fine mesh among all pools
    COMPLEX(kind=DP), ALLOCATABLE :: vmef_all(:,:,:,:)
    !! velocity matrix elements on the fine mesh among all pools
    REAL(DP), ALLOCATABLE :: tdf_sigma_m(:,:,:,:)
    !! transport distribution function
    REAL(DP), ALLOCATABLE :: wkf_all(:)
    !! k-point weight on the full grid across all pools
    !  SP - Uncomment to use symmetries on velocities
    REAL(kind=DP) :: xkf_tmp (3, nkqtotf)
    !! Temporary k-point coordinate (dummy variable)
    REAL(kind=DP) :: wkf_tmp(nkqtotf)
    !! Temporary k-weights (dummy variable)
    REAL(kind=DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(kind=DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(kind=DP) :: sa(3,3), sb(3,3), sr(3,3)
    CALL start_clock('MOB')
    !
    inv_cell = 1.0d0/omega
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF ( system_2d ) &
       inv_cell = inv_cell * at(3,3) * alat
  
    ! 
    ! We can read the scattering rate from files. 
    IF ( scatread ) THEN
      conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
      !
      ! Compute the Fermi level 
      DO itemp = 1, nstemp
        ! 
        etemp = transp_temp(itemp)
        ! 
        ! Lets gather the velocities from all pools
#ifdef __MPI
        IF ( vme ) THEN 
          IF ( .not. ALLOCATED(vmef_all) )  ALLOCATE( vmef_all(3,nbndsub,nbndsub,nkqtotf) )
          vmef_all(:,:,:,:) = czero
          CALL poolgatherc4 ( 3, nbndsub, nbndsub, nkqtotf, 2*nkf, vmef, vmef_all )
        ELSE
          IF ( .not. ALLOCATED(dmef_all) )  ALLOCATE( dmef_all(3,nbndsub,nbndsub,nkqtotf) )
          dmef_all(:,:,:,:) = czero
          CALL poolgatherc4 ( 3, nbndsub, nbndsub, nkqtotf, 2*nkf, dmef, dmef_all )
        ENDIF
        IF ( .not. ALLOCATED(wkf_all) )  ALLOCATE( wkf_all(nkqtotf) )
        wkf_all(:) = zero
        CALL poolgather2 ( 1, nkqtotf, 2*nkf, wkf, wkf_all  )
#else
        IF ( vme ) THEN
          vmef_all = vmef
        ELSE
          dmef_all = dmef
        ENDIF
#endif     
        ! 
        ! In this case, the sum over q has already been done. It should therefore be ok 
        ! to do the mobility in sequential. Each cpu does the same thing below
        ALLOCATE ( etf_all ( nbndsub, nkqtotf/2 ) )
        !
        CALL scattering_read(etemp, ef0(itemp), etf_all, inv_tau_all)
        ! 
        ! This is hole mobility. ----------------------------------------------------
        IF (int_mob .OR. (ncarrier < -1E5)) THEN
          IF (itemp == 1) THEN
            WRITE(stdout,'(/5x,a)') repeat('=',67)
            WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
            WRITE(stdout,'(5x,a/)') repeat('=',67)
          ENDIF
          !      
          IF ( itemp .eq. 1 ) THEN        
            IF ( .not. ALLOCATED(tdf_sigma_m) )  ALLOCATE( tdf_sigma_m(3,3,ibndmax-ibndmin+1,nkqtotf) )
            tdf_sigma_m(:,:,:,:) = zero
            Sigma_m(:,:,:)   = zero
          ENDIF
          !
          DO ik = 1, nkqtotf/2 
            ikk = 2 * ik - 1
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF ( minval ( abs(etf_all (:, ik) - ef ) ) < fsthick ) THEN
              DO ibnd = 1, ibndmax-ibndmin+1
                ! This selects only valence bands for hole conduction
                IF (etf_all (ibndmin-1+ibnd, ik) < ef0(itemp)  ) THEN
                  IF ( vme ) THEN
                    vkk(:,ibnd) = REAL (vmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk)) 
                  ELSE
                    vkk(:,ibnd) = 2.0 * REAL (dmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ENDIF
                  tau = one / inv_tau_all(itemp,ibnd,ik)
                  ekk = etf_all (ibndmin-1+ibnd, ik) -  ef0(itemp)
                  ! 
                  DO j = 1, 3
                    DO i = 1, 3
                      tdf_sigma_m(i,j,ibnd,ik) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  !
                  ! derivative Fermi distribution
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                  !
                  ! electrical conductivity matrix
                  Sigma_m(:,:,itemp) = Sigma_m(:,:,itemp) +  wkf_all(ikk) * dfnk * tdf_sigma_m(:,:,ibnd,ik)
                  !
                ENDIF ! valence bands
              ENDDO ! ibnd
            ENDIF ! fstick
          ENDDO ! ik
          ! 
          carrier_density = 0.0
          ! 
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only valence bands for hole conduction
              IF (etf_all (ibndmin-1+ibnd, ik) < ef0(itemp)  ) THEN
                !  energy at k (relative to Ef)
                ekk = etf_all (ibndmin-1+ibnd, ik) - ef0(itemp)
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ikk) * (1.0d0 - fnk )
              ENDIF
            ENDDO
          ENDDO
          ! 
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3,Sigma_m(:,:,itemp),3,sigma_eig,sigma_vect)
          ! 
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                  ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
          !
        ENDIF ! int_mob .OR. (ncarrier < -1E5)
        ! 
        ! This is electron mobility. ----------------------------------------------------
        IF (int_mob .OR. (ncarrier > 1E5)) THEN
          IF (itemp == 1) THEN
            WRITE(stdout,'(/5x,a)') repeat('=',67)
            WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Electron density [cm^-3]  Electron mobility [cm^2/Vs]")')
            WRITE(stdout,'(5x,a/)') repeat('=',67)
          ENDIF
          !      
          IF ( itemp .eq. 1 ) THEN
            IF ( .not. ALLOCATED(tdf_sigma_m) )  ALLOCATE( tdf_sigma_m(3,3,ibndmax-ibndmin+1,nkqtotf) )
            tdf_sigma_m(:,:,:,:) = zero
            Sigma_m(:,:,:)   = zero
          ENDIF
          !
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF ( minval ( abs(etf_all (:, ik) - ef ) ) < fsthick ) THEN
              DO ibnd = 1, ibndmax-ibndmin+1
                ! This selects only conduction bands for electron conduction
                IF (etf_all (ibndmin-1+ibnd, ik) > ef0(itemp)  ) THEN
                  IF ( vme ) THEN 
                    vkk(:,ibnd) = REAL (vmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ELSE 
                    vkk(:,ibnd) = 2.0 * REAL (dmef_all (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ENDIF
                  tau = one / inv_tau_all(itemp,ibnd,ik)
                  ekk = etf_all (ibndmin-1+ibnd, ik) -  ef0(itemp)
                  ! 
                  DO j = 1, 3
                    DO i = 1, 3
                      tdf_sigma_m(i,j,ibnd,ik) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    ENDDO
                  ENDDO
                  !
                  ! derivative Fermi distribution
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                  !
                  ! electrical conductivity matrix
                  Sigma_m(:,:,itemp) = Sigma_m(:,:,itemp) +  wkf_all(ikk) * dfnk * tdf_sigma_m(:,:,ibnd,ik)
                  !
                ENDIF ! valence bands
              ENDDO ! ibnd
            ENDIF ! fstick
          ENDDO ! ik
          ! 
          carrier_density = 0.0
          ! 
          DO ik = 1, nkqtotf/2
            ikk = 2 * ik - 1
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only conduction bands for electron conduction
              IF (etf_all (ibndmin-1+ibnd, ik) > ef0(itemp)  ) THEN
                !  energy at k (relative to Ef)
                ekk = etf_all (ibndmin-1+ibnd, ik) - ef0(itemp)
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ikk) * fnk 
              ENDIF
            ENDDO
          ENDDO
          ! 
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3,Sigma_m(:,:,itemp),3,sigma_eig,sigma_vect)
          ! 
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                  ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
          !
        ENDIF ! int_mob .OR. (ncarrier > 1E5)
        ! 
      ENDDO ! itemp
      !
    ELSE ! Case without reading the scattering rates from files.
      !
      !  SP - Uncomment to use symmetries on velocities
      IF (mp_mesh_k) THEN
        IF ( mpime .eq. ionode_id ) THEN
          ! 
          CALL set_sym_bl( )
          BZtoIBZ(:) = 0
          s_BZtoIBZ(:,:,:) = 0
          ! What we get from this call is BZtoIBZ
          CALL kpoint_grid_epw ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
                     nkf1,nkf2,nkf3, nkqtotf_tmp, xkf_tmp, wkf_tmp,BZtoIBZ,s_BZtoIBZ)
          ! 
          IF (iterative_bte) THEN
            ! Now we have to remap the points because the IBZ k-points have been
            ! changed to to load balancing. 
            BZtoIBZ_tmp(:) = 0
            DO ikbz=1, nkf1*nkf2*nkf3
              BZtoIBZ_tmp(ikbz) = map_rebal( BZtoIBZ( ikbz ) )
            ENDDO
            BZtoIBZ(:) = BZtoIBZ_tmp(:) 
          ENDIF
          !
        ENDIF
        CALL mp_bcast( s_BZtoIBZ, ionode_id, inter_pool_comm )
        CALL mp_bcast( BZtoIBZ, ionode_id, inter_pool_comm )
        ! 
        !WRITE(stdout,*)'s ',s
        !WRITE(stdout,*)'BZtoIBZ ',BZtoIBZ
        !WRITE(stdout,*)'s_BZtoIBZ ',s_BZtoIBZ
        ! 
      ENDIF
      !
      ! This is hole mobility. In the case of intrinsic mobilities we can do both
      ! electron and hole mobility because the Fermi level is the same. This is not
      ! the case for doped mobilities.
      ! 
      ! find the bounds of k-dependent arrays in the parallel case in each pool
      CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
      ! 
      IF (int_mob .OR. (ncarrier < -1E5)) THEN
        ! 
        DO itemp = 1, nstemp
          !
          etemp = transp_temp(itemp)
          !DBSP
          !write(stdout,*)'etemp ',etemp 
          !write(stdout,*)'inv_tau_all ', SUM(inv_tau_all(itemp,:,:))
          !write(stdout,*)'inv_tau_all ', SUM(inv_tau_all(:,:,:))
          !
          IF ( itemp .eq. 1 ) THEN 
            !
            ! tdf_sigma_ij(ibnd,ik) = v_i(ik,ibnd) * v_j(ik,ibnd) * tau(ik,ibnd)
            ! i,j - cartesian components and ij combined (i,j) index
            ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
            ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
            ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
            ! this can be reduced to 6 if we take into account symmetry xy=yx, ...
            tdf_sigma(:)  = zero
            Sigma(:,:)    = zero
            SigmaZ(:,:)   = zero
            !
          ENDIF
          !
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            !
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                !
                ! This selects only valence bands for hole conduction
                IF (etf (ibndmin-1+ibnd, ikk) < ef0(itemp) ) THEN 
                  !
                  ! vkk(3,nbnd) - velocity for k
                  tdf_sigma(:) = zero
                  IF ( vme ) THEN
                    ! vmef is in units of Ryd * bohr
                    vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ELSE 
                    ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                    ! 1/m  = 2 in Rydberg atomic units
                    ! dmef is in units of 1/a.u. (where a.u. is bohr)
                    ! v_(k,i) is in units of Ryd * a.u.
                    vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                  ENDIF
                  ! Use symmetries on k-point (from Homogeneous grid only)
                  IF (mp_mesh_k) THEN
                    !
                    vk_cart(:) = vkk(:,ibnd)
                    ! 
                    ! Loop on full BZ 
                    nb = 0
                    DO ikbz=1, nkf1*nkf2*nkf3
                      ! If the k-point from the full BZ is related by a symmetry operation 
                      ! to the current k-point, then take it.  
                      IF (BZtoIBZ(ikbz) == ik+lower_bnd-1) THEN
                        nb = nb + 1
                        ! Transform the symmetry matrix from Crystal to cartesian
                        sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                        sb = matmul ( bg, sa )
                        sr (:,:) = matmul ( at, transpose (sb) )
                        CALL dgemv( 'n', 3, 3, 1.d0,&
                          sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                        ij = 0
                        DO j = 1, 3
                          DO i = 1, 3
                            ij = ij + 1
                            ! The factor two in the weight at the end is to account for spin
                            IF (noncolin) THEN
                              tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                            ELSE
                              tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDIF
                    ENDDO ! ikbz 
                    IF (noncolin) THEN
                      IF (ABS(nb*1.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                        CALL errore ('transport', &
                               &' The number of kpoint in the IBZ is not equal to the weight', 1)
                      ENDIF
                    ELSE
                      IF (ABS(nb*2.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                        CALL errore ('transport', &
                               &' The number of kpoint in the IBZ is not equal to the weight', 1)
                      ENDIF
                    ENDIF 
                  ! withtout symmetries
                  ELSE
                    ! 
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * wkf(ikk)
                      ENDDO
                    ENDDO
                    !
                  ENDIF ! mp_mesh_k
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
                  ! derivative Fermi distribution
                  ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                  dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                  ! electrical conductivity
                  Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:) * tau
                  !
                  ! Now do the same but with Znk multiplied
                  ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                  Znk = one / ( one + zi_allvb (itemp,ibnd,ik+lower_bnd-1) )
                  tau = one / ( Znk * inv_tau_all(itemp,ibnd,ik+lower_bnd-1) )
                  SigmaZ(:,itemp) = SigmaZ(:,itemp) + dfnk * tdf_sigma(:) * tau
   
                  !print*,'itemp ik ibnd ',itemp, ik, ibnd
                  !print*,'Sigma ',Sigma(:,itemp)
                  !print*,'SigmaZ ',SigmaZ(:,itemp)
                  !print*,'Znk ',Znk
                  !
                ENDIF
                !
              ENDDO ! ibnd
              !
            ENDIF ! endif  fsthick
            !
          ENDDO ! end loop on k
          !
          ! The k points are distributed among pools: here we collect them
          !
          CALL mp_sum( Sigma(:,itemp),  world_comm )
          CALL mp_sum( SigmaZ(:,itemp), world_comm )
          !DBSP
          !write(stdout,*) 'ef0(itemp) ',ef0(itemp)    
          !write(stdout,*) 'Sigma ',SUM(Sigma(:,itemp))    
          !
        ENDDO ! nstemp
        !
        IF (mpime .eq. meta_ionode_id) THEN
          filsigma = TRIM(prefix) // '_elcond_h'
          OPEN(iufilsigma, file = filsigma, form = 'formatted')
          WRITE(iufilsigma,'(a)') "# Electrical conductivity in 1/(Ohm * m)"
          WRITE(iufilsigma,'(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // & 
                                                   "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                   "        Sigma_xz        Sigma_yz        Sigma_zz"
        ENDIF
        !
        conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
        !
        WRITE(stdout,'(/5x,a)') repeat('=',67)
        WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
        WRITE(stdout,'(5x,a/)') repeat('=',67)
        ! 
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          ! Sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
          IF (mpime.eq. meta_ionode_id) THEN
            WRITE(iufilsigma,'(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                         conv_factor1 * Sigma(:,itemp) * inv_cell
          ENDIF
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            DO ibnd = 1, ibndmax-ibndmin+1
              ! This selects only valence bands for hole conduction
              IF (etf (ibndmin-1+ibnd, ikk) < ef0(itemp) ) THEN
                !  energy at k (relative to Ef)
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)      
                fnk = wgauss( -ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf(ikk) * (1.0d0 - fnk ) 
              ENDIF
            ENDDO
          ENDDO 
          ! 
          CALL mp_sum( carrier_density, world_comm )
          !
          ! Diagonalize the conductivity matrix
          ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
          ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
          ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
          sigma_up(:,:) = zero
          sigma_up(1,1) = Sigma(1,itemp)
          sigma_up(1,2) = Sigma(2,itemp)
          sigma_up(1,3) = Sigma(3,itemp)
          sigma_up(2,1) = Sigma(4,itemp)
          sigma_up(2,2) = Sigma(5,itemp)
          sigma_up(2,3) = Sigma(6,itemp)
          sigma_up(3,1) = Sigma(7,itemp)
          sigma_up(3,2) = Sigma(8,itemp)
          sigma_up(3,3) = Sigma(9,itemp)
          ! 
          CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
          ! 
          !Sigma_diag = (Sigma(1,itemp)+Sigma(5,itemp)+Sigma(9,itemp))/3
          !Sigma_offdiag = (Sigma(2,itemp)+Sigma(3,itemp)+Sigma(4,itemp)+&
          !                 Sigma(6,itemp)+Sigma(7,itemp)+Sigma(8,itemp))/6
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
                  ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg' 
          ! 
  !        ! Now do Znk ----------------------------------------------------------
  !        sigma_up(:,:) = zero
  !        sigma_up(1,1) = SigmaZ(1,itemp)
  !        sigma_up(1,2) = SigmaZ(2,itemp)
  !        sigma_up(1,3) = SigmaZ(3,itemp)
  !        sigma_up(2,1) = SigmaZ(4,itemp)
  !        sigma_up(2,2) = SigmaZ(5,itemp)
  !        sigma_up(2,3) = SigmaZ(6,itemp)
  !        sigma_up(3,1) = SigmaZ(7,itemp)
  !        sigma_up(3,2) = SigmaZ(8,itemp)
  !        sigma_up(3,3) = SigmaZ(9,itemp)
  !        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
  !        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
  !        ! carrier_density in cm^-1
  ! DBSP - Z-factor
  !        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
  !                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'
  
          ! 
        ENDDO ! nstemp
        !
        IF (mpime .eq. meta_ionode_id) CLOSE(iufilsigma)
        !
      ENDIF ! Hole mob
      ! 
      ! Now the electron conduction and mobilities
      ! 
      IF (int_mob .OR. (ncarrier > 1E5)) THEN
        DO itemp = 1, nstemp
          !
          etemp = transp_temp(itemp)
          IF ( itemp .eq. 1 ) THEN
            tdf_sigma(:)  = zero
            Sigma(:,:)    = zero
            SigmaZ(:,:)   = zero
          ENDIF
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
              IF ( ABS(efcb(itemp)) < eps ) THEN  
                DO ibnd = 1, ibndmax-ibndmin+1
                  ! This selects only cond bands for electron conduction
                  IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                    ! vkk(3,nbnd) - velocity for k
                    tdf_sigma(:) = zero
                    IF ( vme ) THEN
                      ! vmef is in units of Ryd * bohr
                      vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                    ELSE
                      ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                      ! 1/m  = 2 in Rydberg atomic units
                      ! dmef is in units of 1/a.u. (where a.u. is bohr)
                      ! v_(k,i) is in units of Ryd * a.u.
                      vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                    ENDIF
                    IF (mp_mesh_k) THEN
                      !
                      vk_cart(:) = vkk(:,ibnd)
                      ! 
                      ! Loop on full BZ 
                      nb = 0
                      DO ikbz=1, nkf1*nkf2*nkf3
                        ! If the k-point from the full BZ is related by a symmetry operation 
                        ! to the current k-point, then take it.  
                        IF (BZtoIBZ(ikbz) == ik+lower_bnd-1) THEN
                          nb = nb + 1
                          ! Transform the symmetry matrix from Crystal to cartesian
                          sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                          sb = matmul ( bg, sa )
                          sr (:,:) = matmul ( at, transpose (sb) )
                          CALL dgemv( 'n', 3, 3, 1.d0,&
                            sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                          ij = 0
                          DO j = 1, 3
                            DO i = 1, 3
                              ij = ij + 1
                              ! The factor two in the weight at the end is to account for spin
                              IF (noncolin) THEN
                                tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                              ELSE
                                tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                      ENDDO ! ikbz 
                      IF (noncolin) THEN
                        IF (ABS(nb*1.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore ('transport', &
                                 &' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ELSE
                        IF (ABS(nb*2.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore ('transport', &
                                 &' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ENDIF
                    ! withtout symmetries
                    ELSE
                      ! 
                      ij = 0
                      DO j = 1, 3
                        DO i = 1, 3
                          ij = ij + 1
                          tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * wkf(ikk)
                        ENDDO
                      ENDDO
                      !
                    ENDIF ! mp_mesh_k
                    !
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                    tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
                    dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                    Sigma(:,itemp) = Sigma(:,itemp) + dfnk * tdf_sigma(:) * tau
                    !
                    ! Now do the same but with Znk multiplied
                    ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                    Znk = one / ( one + zi_allvb (itemp,ibnd,ik+lower_bnd-1) )
                    tau = one / ( Znk * inv_tau_all(itemp,ibnd,ik+lower_bnd-1) )
                    SigmaZ(:,itemp) = SigmaZ(:,itemp) + dfnk * tdf_sigma(:) * tau
                  ENDIF
                ENDDO 
              ELSE ! In this case we have 2 Fermi levels
                DO ibnd = 1, ibndmax-ibndmin+1
                  ! This selects only cond bands for hole conduction
                  IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                    ! 
                    !  SP - Uncomment to use symmetries on velocities
                    tdf_sigma(:) = zero
                    IF ( vme ) THEN
                      vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd,ibndmin-1+ibnd, ikk))
                    ELSE
                      vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd,ibndmin-1+ibnd, ikk))
                    ENDIF
                    ! 
                    IF (mp_mesh_k) THEN
                      !
                      vk_cart(:) = vkk(:,ibnd)
                      ! 
                      ! Loop on full BZ 
                      nb = 0
                      DO ikbz=1, nkf1*nkf2*nkf3
                        ! If the k-point from the full BZ is related by a symmetry operation 
                        ! to the current k-point, then take it.  
                        IF (BZtoIBZ(ikbz) == ik+lower_bnd-1) THEN
                          nb = nb + 1
                          ! Transform the symmetry matrix from Crystal to cartesian
                          sa (:,:) = dble ( s_BZtoIBZ(:,:,ikbz) )
                          sb = matmul ( bg, sa )
                          sr (:,:) = matmul ( at, transpose (sb) )
                          CALL dgemv( 'n', 3, 3, 1.d0,&
                            sr, 3, vk_cart(:),1 ,0.d0 , v_rot(:), 1 )
                          ij = 0
                          DO j = 1, 3
                            DO i = 1, 3
                              ij = ij + 1
                              ! The factor two in the weight at the end is to account for spin
                              IF (noncolin) THEN
                                tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 1.0 / (nkf1*nkf2*nkf3)
                              ELSE
                                tdf_sigma(ij) = tdf_sigma(ij) + ( v_rot(i) * v_rot(j) ) * 2.0 / (nkf1*nkf2*nkf3)
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                      ENDDO ! ikbz 
                      IF (noncolin) THEN
                        IF (ABS(nb*1.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore ('transport', &
                                 &' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ELSE
                        IF (ABS(nb*2.0/(nkf1*nkf2*nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore ('transport', &
                                 &' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ENDIF
                    ! withtout symmetries
                    ELSE
                      ! 
                      ij = 0
                      DO j = 1, 3
                        DO i = 1, 3
                          ij = ij + 1
                          tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * wkf(ikk)
                        ENDDO
                      ENDDO
                      !
                    ENDIF ! mp_mesh_k
                    ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                    tau = one / inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1)
                    dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                    Sigma(:,itemp) = Sigma(:,itemp) +  dfnk * tdf_sigma(:) * tau                      


                    !IF ( vme ) THEN
                    !  vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                    !ELSE
                    !  vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
                    !ENDIF
                    !ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                    !tau = one / inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1)
                    !ij = 0
                    !DO j = 1, 3
                    !  DO i = 1, 3
                    !    ij = ij + 1
                    !    tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                    !  ENDDO
                    !ENDDO
                    !dfnk = w0gauss( ekk / etemp, -99 ) / etemp
                    !Sigma(:,itemp) = Sigma(:,itemp) +  wkf(ikk) * dfnk * tdf_sigma(:)
                    !
                    ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                    Znk = one / ( one + zi_allcb (itemp,ibnd,ik+lower_bnd-1) )
                    tau = one / ( Znk * inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1) )
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i,ibnd) * vkk(j,ibnd) * tau
                      ENDDO
                    ENDDO
                    SigmaZ(:,itemp) = SigmaZ(:,itemp) + wkf(ikk) * dfnk * tdf_sigma(:)                  
                  ENDIF 
                ENDDO ! ibnd
              ENDIF ! etcb
            ENDIF ! endif  fsthick
          ENDDO ! end loop on k
          CALL mp_sum( Sigma(:,itemp),  world_comm )
          CALL mp_sum( SigmaZ(:,itemp), world_comm )
          ! 
        ENDDO ! nstemp
        IF (mpime .eq. meta_ionode_id) THEN
          filsigma = TRIM(prefix) // '_elcond_e'
          OPEN(iufilsigma, file = filsigma, form = 'formatted')
          WRITE(iufilsigma,'(a)') "# Electrical conductivity in 1/(Ohm * m)"
          WRITE(iufilsigma,'(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // &
                                                   "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                  "        Sigma_xz        Sigma_yz        Sigma_zz"
        ENDIF
        !
        conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
        WRITE(stdout,'(/5x,a)') repeat('=',67)
        WRITE(stdout,'(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
        WRITE(stdout,'(5x,a/)') repeat('=',67)
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          IF (mpime .eq. meta_ionode_id) THEN
            ! Sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
            IF ( ABS(efcb(itemp)) < eps ) THEN 
              WRITE(iufilsigma,'(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                           conv_factor1 * Sigma(:,itemp) * inv_cell
            ELSE
              WRITE(iufilsigma,'(11E16.8)') efcb(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                           conv_factor1 * Sigma(:,itemp) * inv_cell
            ENDIF
          ENDIF
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            DO ibnd = 1, ibndmax-ibndmin+1
              ikk = 2 * ik - 1
              ! This selects only conduction bands for electron conduction
              IF ( ABS(efcb(itemp)) < eps ) THEN 
                IF (etf (ibndmin-1+ibnd, ikk) > ef0(itemp) ) THEN
                  ekk = etf (ibndmin-1+ibnd, ikk) - ef0(itemp)
                  fnk = wgauss( -ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ELSE
                IF (etf (ibndmin-1+ibnd, ikk) > efcb(itemp) ) THEN
                  ekk = etf (ibndmin-1+ibnd, ikk) - efcb(itemp)
                  fnk = wgauss( -ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          CALL mp_sum( carrier_density, world_comm )
          ! Diagonalize the conductivity matrix
          ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
          ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
          ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
          sigma_up(:,:) = zero
          sigma_up(1,1) = Sigma(1,itemp)
          sigma_up(1,2) = Sigma(2,itemp)
          sigma_up(1,3) = Sigma(3,itemp)
          sigma_up(2,1) = Sigma(4,itemp)
          sigma_up(2,2) = Sigma(5,itemp)
          sigma_up(2,3) = Sigma(6,itemp)
          sigma_up(3,1) = Sigma(7,itemp)
          sigma_up(3,2) = Sigma(8,itemp)
          sigma_up(3,3) = Sigma(9,itemp)
          ! 
          CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
          ! 
          !Sigma_diag = (Sigma(1,itemp)+Sigma(5,itemp)+Sigma(9,itemp))/3
          !Sigma_offdiag = (Sigma(2,itemp)+Sigma(3,itemp)+Sigma(4,itemp)+&
          !                 Sigma(6,itemp)+Sigma(7,itemp)+Sigma(8,itemp))/6
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          !
       
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * ( bohr2ang * ang2cm  )**(-3)
          IF ( ABS(efcb(itemp)) < eps ) THEN
            WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                     ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          ELSE
            WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                     efcb(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          ENDIF
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg'
          ! Issue warning if the material is anisotropic
         ! IF (Sigma_offdiag > 0.1*Sigma_diag) THEN
         !   WRITE(stdout,'(5x,a,1f10.5,a)') 'Warning: Sigma_offdiag = ',(Sigma_offdiag*100)/Sigma_diag, '% of Sigma_diag'
         ! ENDIF
          ! Now do the mobility with Znk factor ----------------------------------------------------------
          sigma_up(:,:) = zero
          sigma_up(1,1) = SigmaZ(1,itemp)
          sigma_up(1,2) = SigmaZ(2,itemp)
          sigma_up(1,3) = SigmaZ(3,itemp)
          sigma_up(2,1) = SigmaZ(4,itemp)
          sigma_up(2,2) = SigmaZ(5,itemp)
          sigma_up(2,3) = SigmaZ(6,itemp)
          sigma_up(3,1) = SigmaZ(7,itemp)
          sigma_up(3,2) = SigmaZ(8,itemp)
          sigma_up(3,3) = SigmaZ(9,itemp)
          CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
          mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  / ( carrier_density * hbarJ)
          mobility = (mobility_xx+mobility_yy+mobility_zz)/3
          !
  ! DBSP - Z-factor
  !        IF ( ABS(efcb(itemp)) < eps ) THEN
  !          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
  !                                                   ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        ELSE
  !          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
  !                                                   efcb(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        ENDIF
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'
  
          ! 
        ENDDO ! nstemp
        WRITE(stdout,'(5x)')
        WRITE(stdout,'(5x,"Note: Mobility are sorted by ascending values and might not correspond to the expected (x,y,z) axis.")')
        !
        IF (mpime .eq. meta_ionode_id) CLOSE(iufilsigma)
        ! 
      ENDIF ! Electron mobilities
    ENDIF ! scatread
    !
    CALL stop_clock('MOB')
    ! 
    WRITE( stdout,  * ) '    Total time so far'
    CALL print_clock ('SCAT')
    CALL print_clock ('MOB')
    WRITE(stdout,'(5x)')
    ! 
    !! IF IBTE we want the SRTA solution to be the first iteration of IBTE
    !IF (iterative_bte) THEN
    !  Fi_all(:,:,:,:) = zero
    !  IF ( ABS(efcb(1)) > eps ) Fi_allcb(:,:,:,:) = zero
    !  ! 
    !  DO itemp = 1, nstemp
    !    DO ik = 1, nkf
    !      ikk = 2 * ik - 1
    !      IF ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) THEN
    !        DO ibnd = 1, ibndmax-ibndmin+1
    !          IF ( vme ) THEN
    !            vkk(:,ibnd) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
    !          ELSE
    !            vkk(:,ibnd) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk))
    !          ENDIF
    !          tau = one / inv_tau_all(itemp,ibnd,ik+lower_bnd-1)
    !          Fi_all(:,ibnd,ik+lower_bnd-1,itemp) = vkk(:,ibnd) * tau
    !          IF ( ABS(efcb(itemp)) > eps ) THEN
    !            tau = one / inv_tau_allcb(itemp,ibnd,ik+lower_bnd-1)
    !            Fi_allcb(:,ibnd,ik+lower_bnd-1,itemp) = vkk(:,ibnd) * tau 
    !          ENDIF
    !          !
    !        ENDDO
    !      ENDIF
    !    ENDDO ! kpoints
    !    !DBSP
    !    print*,'itemp ',itemp,' ',sum(Fi_all(:,:,:,itemp)), sum(inv_tau_all)
    !    print*,'itemp cb',itemp,' ',sum(Fi_allcb(:,:,:,itemp)), sum(inv_tau_allcb)
    !  ENDDO ! itemp
    !  CALL mp_sum( Fi_all, world_comm )
    !  IF ( ABS(efcb(1)) > eps ) CALL mp_sum( Fi_allcb, world_comm )
    !ENDIF
    !
    RETURN
    !
    END SUBROUTINE transport_coeffs
    !--------------------------------------------------------------------------
    ! 
  END MODULE transport
