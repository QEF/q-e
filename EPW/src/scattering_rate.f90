  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !       
  !-----------------------------------------------------------------------
  SUBROUTINE scattering_rate_q( iq, ef0, efcb, first_cycle ) 
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the scattering rate (inv_tau)
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, meta_ionode_id
  USE io_epw,        ONLY : iufilscatt_rate
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, etf_mem, efermi_read, lrepmatf, & 
                            eps_acustic, fermi_energy, ngaussw, degaussw, & 
                            nstemp, scattering_serta, scattering_0rta, shortrange,&
                            restart, restart_freq, restart_filq
  USE pwcom,         ONLY : ef, nelec, isk
  USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, wkf, dmef, wf, wqf, xkf, & 
                            epf17, efnew, nqtotf, nkqtotf, inv_tau_all, inv_tau_allcb, &
                            xqf, zi_allvb, zi_allcb
  USE transportcom,  ONLY : transp_temp, lower_bnd, upper_bnd
  USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, & 
                            meV2invps, eps6
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : world_comm
  USE mp_world,      ONLY : mpime
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
  INTEGER :: ikq
  !! Even k+q index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: jbnd
  !! Local band index
  INTEGER :: imode
  !! Local mode index
  INTEGER :: icbm
  !! Index of the CBM
  INTEGER :: nrec
  !! Record index
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
  REAL(KIND=DP) :: temp
  !! Temporary file name used to write scattering rate to file. 
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
  CHARACTER (len=256) :: name1
  !! Name used to write scattering rates to file. 
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
            wgq = wgauss( -wq*inv_etemp, -99)
            wgq = wgq / ( one - two * wgq )
            !
            ! SP : Define the inverse for efficiency
            inv_wq =  1.0/(two * wq)
            ! SP : Avoid if statement in inner loops
            ! the coupling from Gamma acoustic phonons is negligible
            IF ( wq .gt. eps_acustic ) THEN
              g2_tmp = 1.0
            ELSE
              g2_tmp = 0.0
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
                  g2 = (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp
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
               ! if ((ik .eq. 27) .and. (ibnd > 4) .and. trans_prob > 1E-7) then
               !    print*,'trans_prob ', trans_prob
               !    print*,'inv_tau',SUM(inv_tau(1,5:8,27))
               ! end if
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
                    g2 = (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp
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
      IF (MOD(iq,restart_freq) == 0) THEN
        WRITE(stdout, '(a)' ) '     Creation of a restart point'
        ! 
        ! The mp_sum will aggreage the results on each k-points. 
        CALL mp_sum( inv_tau_all, world_comm )
        CALL mp_sum( zi_allvb,    world_comm )
        !
        CALL tau_write(iq,nqtotf,nkqtotf/2,.FALSE.)
        ! 
        IF ( ABS(efcb(1)) > eps ) THEN
          ! 
          CALL mp_sum( inv_tau_allcb, world_comm ) 
          CALL mp_sum( zi_allcb,      world_comm ) 
          ! 
          CALL tau_write(iq,nqtotf,nkqtotf/2,.TRUE.) 
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
&                            + inv_tau_all_new(:,:,:) * nqtotf_new ) / (nqtotf+nqtotf_new)
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
  ENDIF
  ! DBSP
  !write(stdout,*),'iq ',iq
  !write(stdout,*),'inv_tau_all(1,5:8,21) ',SUM(inv_tau_all(3,5:8,21))
  !write(stdout,*),'inv_tau_all(1,5:8,:) ',SUM(inv_tau_all(3,5:8,:))
  !write(stdout,*),'SUM(inv_tau_all) ',SUM(inv_tau_all(3,:,:))
  !write(stdout,*),'first_cycle ',first_cycle
  !
  RETURN
  !
  END SUBROUTINE scattering_rate_q
  !-----------------------------------------------------------------------
