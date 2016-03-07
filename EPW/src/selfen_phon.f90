  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_phon_q (iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
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
  !  RM 24/02/2014
  !  redefined the size of coskkq, vkk, vkq within the fermi windwow
  !  cleaned up the subroutine
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_epw,     ONLY : iunepmatf, iuetf
  use phcom,      ONLY : nmodes
  USE control_lr, ONLY : lgamma
  use epwcom,     ONLY : nbndsub, lrepmatf, fsthick, &
                         eptemp, ngaussw, degaussw, &
                         etf_mem, nsmear, delta_smear, eps_acustic, &
                         efermi_read, fermi_energy, delta_approx
  use pwcom,      ONLY : nelec, ef, isk
  use elph2,      ONLY : epf17, ibndmax, ibndmin, etf, &
                         wkf, xqf, wqf, nkqf, nqtotf,   &
                         nkf, wf, nkqtotf, xqf, nqf, &
                         lambda_all, lambda_v_all, &
                         dmef, gamma_all,gamma_v_all, efnew
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi
#ifdef __PARA
  use mp,         ONLY : mp_barrier,mp_sum
  use mp_global,  ONLY : me_pool,inter_pool_comm,my_pool_id,npool
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, weight, dosef, &
                   degaussw0, eptemp0, lambda_tot, lambda_tr_tot,&
                   inv_wq, inv_degaussw0, g2_tmp, inv_eptemp0, w0g1, w0g2
  !
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  real(kind=DP) :: gamma(nmodes),gamma_v(nmodes)
  real(kind=DP) :: coskkq(ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: DDOT, vkk(3,ibndmax-ibndmin+1), vkq(3,ibndmax-ibndmin+1)
  !
  !
  IF ( iq .eq. 1 ) THEN 
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
     !
     IF ( .not. ALLOCATED (lambda_all) )    ALLOCATE( lambda_all  (nmodes, nqtotf, nsmear) )
     IF ( .not. ALLOCATED (lambda_v_all) )  ALLOCATE( lambda_v_all(nmodes, nqtotf, nsmear) )
     lambda_all(:,:,:)   = zero
     lambda_v_all(:,:,:) = zero
     IF ( .not. ALLOCATED (gamma_all) )    ALLOCATE( gamma_all  (nmodes,nqtotf,nsmear) )
     IF ( .not. ALLOCATED (gamma_v_all) )  ALLOCATE( gamma_v_all(nmodes,nqtotf,nsmear) )
     gamma_all(:,:,:)   = zero
     gamma_v_all(:,:,:) = zero
     !
  ENDIF
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0   = (ismear-1) * delta_smear + eptemp(1)
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
        !WRITE (stdout, 101) dosef / ryd2ev, ef  * ryd2ev
     ENDIF
     !
     CALL start_clock('PH SELF-ENERGY')
     !
     fermicount = 0
     gamma(:)   = zero
     gamma_v(:) = zero
     !
     DO ik = 1, nkf
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        ! 
        coskkq = 0.d0
        DO ibnd = 1, ibndmax-ibndmin+1
           DO jbnd = 1, ibndmax-ibndmin+1
              ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
              ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
              ! In principle the only coskkq contributing to lambda_tr are both near the
              ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
              ! we may implement the approximation to the angle between k and k+q vectors also 
              ! listed in Grimvall
              !
              ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
              ! 1/m  = 2 in Rydberg atomic units
              !
              vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk ) )
              vkq(:, jbnd ) = 2.0 * REAL (dmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
              IF ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) .gt. 1.d-4) &
                   coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                   DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
           ENDDO
        ENDDO
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        ! when we see references to iq for file readinq, it is always = 1 
        IF (.not. etf_mem) then
           nrec = ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
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
              !  we read the e-p matrix from disk / memory
              !
              IF (etf_mem) then
                 epf(:,:) = epf17 ( ik, :, :, imode)
              ELSE
                 nrec = (imode-1) * nkf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k
                 ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                 IF (delta_approx) THEN
                    w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                 ELSE
                    wgkk = wgauss( -ekk*inv_eptemp0, -99)
                 ENDIF
                 !
                 DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    g2 = (abs(epf (jbnd, ibnd))**two)*inv_wq*g2_tmp
                    !
                    IF (delta_approx) THEN 
                       !
                       w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                       ! the expression below is positive-definite, but also an
                       ! approximation which neglects some fine features
                       weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                       !
                    ELSE
                       !
                       wgkq = wgauss( -ekq*inv_eptemp0, -99)
                       !
                       ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q + id]
                       ! This is the imaginary part of the phonon self-energy, sans
                       ! the matrix elements
                       !
                       !weight = wkf (ikk) * (wgkk - wgkq) * &
                       !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) )
                       !
                       ! SP: The expression below (phonon self-energy) corresponds to
                       !  = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q]
                       weight = pi * wkf (ikk) * (wgkk - wgkq)* &
                            w0gauss ( (ekq - ekk - wq) / degaussw0, 0) / degaussw0
                       !
                    ENDIF  
                    !
                    gamma   (imode) = gamma   (imode) + weight * g2 
                    gamma_v (imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
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
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     ! this finishes the integral over the BZ  (k)
     !
     CALL mp_sum(gamma,inter_pool_comm) 
     CALL mp_sum(gamma_v,inter_pool_comm) 
     CALL mp_sum(fermicount, inter_pool_comm)
     CALL mp_barrier(inter_pool_comm)
     !
#endif
     !
     WRITE(stdout,'(/5x,"ismear = ",i5," iq = ",i7," coord.: ", 3f9.5, " wt: ", f9.5)') ismear, iq, xqf(:,iq), wqf(iq)
     WRITE(stdout,'(5x,a)') repeat('-',67)
     !
     lambda_tot = 0.d0
     lambda_tr_tot = 0.d0
     !
     DO imode = 1, nmodes
        ! 
        wq = wf (imode, iq)
        IF ( wq .gt. eps_acustic ) THEN 
           lambda_all  ( imode, iq, ismear ) = gamma  ( imode ) / pi / wq**two / dosef
           lambda_v_all( imode, iq, ismear ) = gamma_v( imode ) / pi / wq**two / dosef
        ENDIF
        gamma_all  ( imode, iq, ismear ) = gamma  ( imode )
        gamma_v_all( imode, iq, ismear ) = gamma_v( imode )
        lambda_tot    = lambda_tot    + lambda_all  ( imode, iq, ismear )
        lambda_tr_tot = lambda_tr_tot + lambda_v_all( imode, iq, ismear )
        !
        WRITE(stdout, 102) imode, lambda_all(imode,iq,ismear),ryd2mev*gamma_all(imode,iq,ismear), ryd2mev*wq
!        WRITE(stdout, 102) imode, lambda_v_all(imode,iq,ismear),ryd2mev*gamma_v(imode), ryd2mev*wq
        !
     ENDDO
     !
     WRITE(stdout, 103) lambda_tot
 !    WRITE(stdout, 105) lambda_tr_tot
     WRITE(stdout,'(5x,a/)') repeat('-',67)
     ! 
     ! test ONLY
#ifdef __PARA
     IF (me_pool == 0) &
#endif
     !     
     WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
           'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
     !
  ENDDO !smears
  !
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'lambda( ',i3,' )=',f15.6,'   gamma=',f15.6,' meV','   omega=',f12.4,' meV')
103 FORMAT(5x,'lambda( tot )=',f15.6)
104 FORMAT(5x,'lambda_tr( ',i3,' )=',f15.6,'   gamma_tr=',f15.6,' meV','   omega=',f12.4,' meV')
105 FORMAT(5x,'lambda_tr( tot )=',f15.6)
  !
  RETURN
  !
END SUBROUTINE selfen_phon_q
  !                                                                            
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_phon_k (ik )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
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
  !  k-point paralellization
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_epw,     ONLY : iunepmatf, iuetf, lambda_phself, linewidth_phself
  use phcom,      ONLY : nmodes
  USE control_lr, ONLY : lgamma
  use epwcom,     ONLY : nbndsub, lrepmatf, fsthick, &
                         eptemp, ngaussw, degaussw, &
                         etf_mem, nsmear, delta_smear, eps_acustic, &
                         efermi_read, fermi_energy, delta_approx
  use pwcom,      ONLY : nelec, ef, isk
  use elph2,      ONLY : epf17, ibndmax, ibndmin, etf, etf_k, &
                         wkf, xqf, wqf, nkqf, nqtotf,   &
                         nkf, wf, nkqtotf, xqf, nqf, &
                         lambda_all, lambda_v_all, &
                         dmef, gamma_all,gamma_v_all, efnew
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi
#ifdef __PARA
  use mp,         ONLY : mp_barrier, mp_sum, mp_bcast
  use mp_global,  ONLY : me_pool,inter_pool_comm,my_pool_id,npool
  USE mp_world,   ONLY : mpime
  USE io_global,  ONLY : ionode_id
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, weight, dosef, &
                   degaussw0, eptemp0, lambda_tot, lambda_tr_tot,&
                   inv_wq, inv_degaussw0, g2_tmp, inv_eptemp0, w0g1, w0g2
  !
  real(kind=DP), external :: efermig_seq, dos_ef_seq, w0gauss, wgauss
  real(kind=DP) :: coskkq(ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: DDOT, vkk(3,ibndmax-ibndmin+1), vkq(3,ibndmax-ibndmin+1)
  character (len=30)  :: myfmt
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf, lower_bnd, upper_bnd
  REAL(kind=DP), ALLOCATABLE :: xqf_all(:,:), etf_all(:,:), wqf_all(:,:), wf_all(:,:)
  !
  IF ( ik .eq. 1 ) THEN 
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
     !
     IF ( .not. ALLOCATED (lambda_all) )    ALLOCATE( lambda_all  (nmodes, nqtotf, nsmear) )
     IF ( .not. ALLOCATED (lambda_v_all) )  ALLOCATE( lambda_v_all(nmodes, nqtotf, nsmear) )
     lambda_all(:,:,:)   = zero
     lambda_v_all(:,:,:) = zero
     IF ( .not. ALLOCATED (gamma_all) )    ALLOCATE( gamma_all  (nmodes,nqtotf,nsmear) )
     IF ( .not. ALLOCATED (gamma_v_all) )  ALLOCATE( gamma_v_all(nmodes,nqtotf,nsmear) )
     gamma_all(:,:,:)   = zero
     gamma_v_all(:,:,:) = zero
     !
  ENDIF
  !
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1) * delta_smear + degaussw
    eptemp0   = (ismear-1) * delta_smear + eptemp(1)
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
#ifdef __PARA
       IF (mpime .eq. ionode_id) THEN
#endif
         ef0 = efermig_seq(etf_k,nbndsub,nkqf,nelec,wkf,degaussw0,ngaussw,0,isk)
#ifdef __PARA
       ENDIF
       CALL mp_bcast (ef0, ionode_id, inter_pool_comm)
#endif 
       ! if some bands are skipped (nbndskip.neq.0), nelec has already been
       ! recalculated in ephwann_shuffle
       !
    ELSE !SP: This is added for efficiency reason because the efermig routine is slow
       ef0 = efnew
    ENDIF
    !
#ifdef __PARA
    IF (mpime .eq. ionode_id) THEN
#endif
      !   N(Ef) in the equation for lambda is the DOS per spin
      !dosef = dosef / two
      dosef = dos_ef_seq (ngaussw, degaussw0, ef0, etf_k, wkf, nkqf, nbndsub) / two
      !
#ifdef __PARA
    ENDIF
    CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
#endif
    !   N(Ef) in the equation for lambda is the DOS per spin
    !dosef = dosef / two
    !
    IF ( ik .eq. 1 ) THEN 
       WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
       WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
       !WRITE (stdout, 101) dosef / ryd2ev, ef  * ryd2ev
    ENDIF
    !
    ! find the bounds of q-dependent arrays in the parallel case in each pool
    CALL fkbounds( nqtotf, lower_bnd, upper_bnd )
    !
    CALL start_clock('PH SELF-ENERGY')
    !
    fermicount = 0
    !
    DO iq = 1, nqf
       !
       IF (lgamma) THEN
         ikq = iq
         ikk = iq
       ELSE
         ikq = 2 * iq
         ikk = ikq - 1
       ENDIF
       ! 
       coskkq = 0.d0
       DO ibnd = 1, ibndmax-ibndmin+1
          DO jbnd = 1, ibndmax-ibndmin+1
             ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
             ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
             ! In principle the only coskkq contributing to lambda_tr are both near the
             ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
             ! we may implement the approximation to the angle between k and k+q vectors also 
             ! listed in Grimvall
             !
             ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
             ! 1/m  = 2 in Rydberg atomic units
             !
             vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk ) )
             vkq(:, jbnd ) = 2.0 * REAL (dmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
             IF ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) .gt. 1.d-4) &
                  coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                  DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
          ENDDO
       ENDDO
       !
       ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
       !
       ! when we see references to iq for file readinq, it is always = 1 
       IF (.not. etf_mem) then
          nrec = ikk
          CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
          nrec = ikq
          CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
       ENDIF
       !
       ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
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
             !  we read the e-p matrix from disk / memory
             !
             IF (etf_mem) then
                epf(:,:) = epf17 ( iq, :, :, imode)
             ELSE
                nrec = (imode-1) * nqf + iq
                CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
             ENDIF
             !
             DO ibnd = 1, ibndmax-ibndmin+1
                !
                !  the fermi occupation for k
                ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                IF (delta_approx) THEN
                   w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                ELSE
                   wgkk = wgauss( -ekk*inv_eptemp0, -99)
                ENDIF
                !
                DO jbnd = 1, ibndmax-ibndmin+1
                   !
                   !  the fermi occupation for k+q
                   ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                   !
                   ! here we take into account the zero-point sqrt(hbar/2M\omega)
                   ! with hbar = 1 and M already contained in the eigenmodes
                   ! g2 is Ry^2, wkf must already account for the spin factor
                   !
                   g2 = (abs(epf (jbnd, ibnd))**two)*inv_wq*g2_tmp
                   !
                   IF (delta_approx) THEN 
                      !
                      w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                      ! the expression below is positive-definite, but also an
                      ! approximation which neglects some fine features
                      weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                      !
                   ELSE
                      !
                      wgkq = wgauss( -ekq*inv_eptemp0, -99)
                      !
                      ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q + id]
                      ! This is the imaginary part of the phonon self-energy,
                      ! sans the matrix elements
                      !
                      !weight = wkf (ikk) * (wgkk - wgkq) * &
                      !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) )
                      !
                      ! SP: The expression below (phonon self-energy) corresponds to
                      !  = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q]
                      weight = pi * wkf (ikk) * (wgkk - wgkq)* &
                           w0gauss ( (ekq - ekk - wq) / degaussw0, 0) / degaussw0
                      !
                   ENDIF  

                   gamma_all(imode,iq+lower_bnd-1,ismear) = gamma_all(imode,iq+lower_bnd-1,ismear) + weight * g2 
                   gamma_v_all(imode,iq+lower_bnd-1,ismear) = gamma_v_all(imode,iq+lower_bnd-1,ismear) &
                                                              + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                   ! 
                   IF ( wq .gt. eps_acustic ) THEN
                     lambda_all  ( imode, iq+lower_bnd-1, ismear ) = gamma_all(imode,iq+lower_bnd-1,ismear)&
                                                                     / pi / wq**two / dosef
                     lambda_v_all( imode, iq+lower_bnd-1, ismear ) = gamma_v_all(imode,iq+lower_bnd-1,ismear)&
                                                                     / pi / wq**two / dosef
                   ENDIF
                   !
                ENDDO ! jbnd
                !
             ENDDO   ! ibnd
             !
          ENDDO ! loop on q-modes
          !
       ENDIF ! endif fsthick
       !
    ENDDO ! loop on q
    !
  ENDDO !smears
  !
  CALL stop_clock('PH SELF-ENERGY')
  ! 
  ! The q points are distributed among pools: here we collect them 
  ! when we reach the last k
  !
  IF ( ik .eq. (nkqtotf - nqtotf)) THEN
    !
    ALLOCATE ( xqf_all ( 3,       nqtotf ), &
               wqf_all (1,nqtotf), wf_all(nmodes,nqtotf) )
    xqf_all(:,:) = zero
    wqf_all(:,:) = zero
    !
    wf_all(:,:) = zero
    DO iq = 1, nqf
      DO imode = 1, nmodes
        wf_all(imode,iq+lower_bnd-1) = wf (imode, iq)
      ENDDO
    ENDDO  
    !
#ifdef __PARA
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2 ( 3,       nqtotf, nqf, xqf,    xqf_all  )
    CALL poolgather2 ( 1,       nqtotf, nqf, wqf,    wqf_all  )
    CALL mp_sum( gamma_all, inter_pool_comm )
    CALL mp_sum( gamma_v_all, inter_pool_comm )
    CALL mp_sum( lambda_all, inter_pool_comm )
    CALL mp_sum( lambda_v_all, inter_pool_comm )
    CALL mp_sum(fermicount, inter_pool_comm)
    CALL mp_sum(wf_all, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    xqf_all = xqf
    wqf_all = wqf
    !
#endif
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      OPEN(unit=lambda_phself,file='lambda.phself')
      WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
      WRITE(lambda_phself, *) '#Modes     ',(imode, imode=1,nmodes)
      DO iq = 1, nqtotf
          !
        myfmt = "(*(3x,E15.5))"
        WRITE(lambda_phself,'(i9,4x)',advance='no') iq
        WRITE(lambda_phself, fmt=myfmt) (REAL(lambda_all(imode,iq,1)),imode=1,nmodes)
          !
      ENDDO
      CLOSE(lambda_phself)
      OPEN(unit=linewidth_phself,file='linewidth.phself')
      WRITE(linewidth_phself, '(/2x,a/)') '#Phonon lifetime (meV) '
      WRITE(linewidth_phself,'(2x,a)',advance='no') '#Q-point     '
      Do imode=1, nmodes
        WRITE(linewidth_phself, '(a)',advance='no') '      Mode'
        WRITE(linewidth_phself, '(i3)',advance='no') imode
      enddo
      WRITE(linewidth_phself, '(/2x,a/)') '  '
      DO iq = 1, nqtotf
        !
        myfmt = "(*(3x,E15.5))"
        WRITE(linewidth_phself,'(i9,4x)',advance='no') iq
        WRITE(linewidth_phself, fmt=myfmt) (ryd2mev*REAL(gamma_all(imode,iq,1)), imode=1,nmodes)
        !
      ENDDO
      CLOSE(linewidth_phself)
#ifdef __PARA
    ENDIF
#endif
    !
    DO ismear = 1, nsmear
      lambda_tot = 0.d0
      lambda_tr_tot = 0.d0
      !
      DO iq = 1, nqtotf
        !
        WRITE(stdout,'(/5x,"ismear = ",i5," iq = ",i7," coord.: ", 3f9.5, " wt: ", f9.5)') ismear, iq, xqf_all(:,iq), wqf_all(1,iq)
        WRITE(stdout,'(5x,a)') repeat('-',67)
        !
        lambda_tot = 0.d0
        lambda_tr_tot = 0.d0
        !
        DO imode = 1, nmodes
          ! 
          wq = wf_all (imode, iq)
          lambda_tot    = lambda_tot    + lambda_all  ( imode, iq, ismear )
          lambda_tr_tot = lambda_tr_tot + lambda_v_all( imode, iq, ismear )
          !
          WRITE(stdout, 102) imode, lambda_all(imode,iq,ismear),ryd2mev*gamma_all(imode,iq,ismear), ryd2mev*wq
!          WRITE(stdout, 102) imode, lambda_v_all(imode,iq,ismear),ryd2mev*gamma_v(imode), ryd2mev*wq
          !
        ENDDO
        !
        WRITE(stdout, 103) lambda_tot
 !       WRITE(stdout, 105) lambda_tr_tot
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        ! 
        ! test ONLY
#ifdef __PARA
        IF (me_pool == 0) &
#endif
        !     
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
              'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
        !
      ENDDO 
      !
    ENDDO
  ENDIF  
  ! 
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'lambda( ',i3,' )=',f15.6,'   gamma=',f15.6,' meV','   omega=',f12.4,' meV')
103 FORMAT(5x,'lambda( tot )=',f15.6)
104 FORMAT(5x,'lambda_tr( ',i3,' )=',f15.6,'   gamma_tr=',f15.6,' meV','   omega=',f12.4,' meV')
105 FORMAT(5x,'lambda_tr( tot )=',f15.6)
  !
END SUBROUTINE selfen_phon_k
!
!-----------------------------------------------------------------------
FUNCTION dos_ef_seq (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE mp_pools, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_sum
  IMPLICIT NONE
  REAL(DP) :: dos_ef_seq
  INTEGER :: ngauss, nbnd, nks
  REAL(DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  INTEGER :: ik, ibnd
  REAL(DP), EXTERNAL :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef_seq = 0.0d0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef_seq = dos_ef_seq + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  RETURN
END FUNCTION dos_ef_seq


