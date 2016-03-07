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
  !
  !  Compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  02/06/2013 Modified by Roxana Margine 
  !
  !  This subroutine computes the contribution from phonon iq to all k-points
  !  The outer loop in ephwann_shuffle.f90 will loop over all iq points
  !  The contribution from each iq is summed at the end of this subroutine for iq=nqtotf 
  !  to recover the per-ik electron self energy
  !
  !  RM 24/02/2014
  !  redefined the size of sigmar_all, sigmai_all, and zi_all within the fermi windwow
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iunepmatf, iuetf, linewidth_elself
  USE phcom,         ONLY : nmodes
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, &
                           fsthick, eptemp, ngaussw, degaussw, &
                           etf_mem, eps_acustic, efermi_read, fermi_energy
  USE pwcom,         ONLY : ef !, nelec, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, zi_all, efnew, nqf
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, weight, dosef, eptemp0,&
                   w0g1, w0g2, inv_wq, inv_eptemp0, g2_tmp,&
                   inv_degaussw
  REAL(kind=DP), external :: efermig, dos_ef, wgauss, w0gauss
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf, lower_bnd, upper_bnd
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:)
  !
  ! loop over temperatures can be introduced
  !
  eptemp0 = eptemp(1)
  !
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp0
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
           'Golden Rule strictly enforced with T = ',eptemp0 * ryd2ev, ' eV'
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
  IF ( iq .eq. 1 ) THEN 
     IF ( .not. ALLOCATED (sigmar_all) ) ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (sigmai_all) ) ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nksqtotf) )
     IF ( .not. ALLOCATED (zi_all) )     ALLOCATE( zi_all(ibndmax-ibndmin+1, nksqtotf) )
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
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
     ! when we see references to iq, it is always = 1 
     !
     IF (.not. etf_mem) THEN
        nrec = ikk
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = ikq
        CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
     ENDIF
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
           !  we read the e-p matrix
           !
           IF (etf_mem) THEN
              epf(:,:) = epf17 ( ik, :, :, imode)
           ELSE
              nrec = (imode-1) * nkf + ik
              CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
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
                 g2 = (abs(epf (jbnd, ibnd))**two)*inv_wq*g2_tmp
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
                 ! Delta implementation 
                 w0g1=w0gauss( (ekk-ekq+wq)/degaussw, 0) /degaussw
                 w0g2=w0gauss( (ekk-ekq-wq)/degaussw, 0) /degaussw
                 weight = pi * wqf(iq) * ( (wgkq+wgq)*w0g1 + (one-wgkq+wgq)*w0g2 )
                 !
                 sigmai_all(ibnd,ik+lower_bnd-1) = sigmai_all(ibnd,ik+lower_bnd-1) + g2 * weight
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
  ! The k points are distributed among pools: here we collect them
  !
  IF ( iq .eq. nqtotf ) THEN
     !
     ALLOCATE ( xkf_all      ( 3,       nkqtotf ), &
                etf_all      ( nbndsub, nkqtotf ) )
     xkf_all(:,:) = zero
     etf_all(:,:) = zero
     !
#ifdef __PARA
     !
     ! note that poolgather2 works with the doubled grid (k and k+q)
     !
     CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
     CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
     CALL mp_sum( sigmar_all, inter_pool_comm )
     CALL mp_sum( sigmai_all, inter_pool_comm )
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
     ! Output electron SE here after looping over all q-points (with their contributions 
     ! summed in sigmar_all, etc.)
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     !
     ! Write to file
     OPEN(unit=linewidth_elself,file='linewidth.elself')
     WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
     WRITE(linewidth_elself, '(a)') '#    ik     E(ibnd)    Im(Sgima)(meV)'
     ! 
     DO ik = 1, nksqtotf
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
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
!          WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ik), &
!                             ryd2mev * sigmai_all (ibnd,ik), zi_all (ibnd,ik)

          WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
          WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
          WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
          WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_all(ibnd,ik)
          !
        ENDDO
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        !
     ENDDO
     !
     DO ibnd = 1, ibndmax-ibndmin+1
        !
        DO ik = 1, nksqtotf
           !
           IF (lgamma) THEN
              ikk = ik
              ikq = ik
           ELSE
              ikk = 2 * ik - 1
              ikq = ikk + 1
           ENDIF
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
     !
  ENDIF 
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  103 FORMAT(5x,'k( ',i7,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_q
  !
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_k ( ik )
  !-----------------------------------------------------------------------
  !
  !  Compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iunepmatf, iuetf, linewidth_elself
  USE phcom,         ONLY : nmodes
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, &
                           fsthick, eptemp, ngaussw, degaussw, &
                           etf_mem, eps_acustic, efermi_read, fermi_energy
  USE pwcom,         ONLY : ef !, nelec, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, etf_k, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf, &
                            sigmar_all, sigmai_all, zi_all, efnew, nqf
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id
#endif
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, weight, dosef, eptemp0,&
                   w0g1, w0g2, inv_wq, inv_eptemp0, g2_tmp,&
                   inv_degaussw
  REAL(kind=DP), external :: efermig, dos_ef, wgauss, w0gauss, dos_ef_seq
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf, lower_bnd, upper_bnd
  !
  ! loop over temperatures can be introduced
  !
  eptemp0 = eptemp(1)
  !
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp0
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
           'Golden Rule strictly enforced with T = ',eptemp0 * ryd2ev, ' eV'
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     OPEN(unit=linewidth_elself,file='linewidth.elself')
     WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
     WRITE(linewidth_elself, '(a)') '#    ik     E(ibnd)    Im(Sgima)(meV)'
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
#ifdef __PARA
  IF (mpime .eq. ionode_id) THEN
#endif
    !   N(Ef) in the equation for lambda is the DOS per spin
    !dosef = dosef / two
    dosef = dos_ef_seq (ngaussw, degaussw, ef0, etf_k, wkf, nkqf, nbndsub)/2
    !
#ifdef __PARA
  ENDIF
  CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
#endif  
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
     IF (lgamma) THEN
        ikk = iq
        ikq = iq
     ELSE
        ikq = 2 * iq
        ikk = ikq - 1
     ENDIF
     !
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
     ! when we see references to iq, it is always = 1 
     !
     IF (.not. etf_mem) THEN
        nrec = ikk
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = ikq
        CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
     ENDIF
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
           !  we read the e-p matrix
           !
           IF (etf_mem) THEN
              epf(:,:) = epf17 ( iq, :, :, imode)
           ELSE
              nrec = (imode-1) * nqf + iq
              CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
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
                 g2 = (abs(epf (jbnd, ibnd))**two)*inv_wq*g2_tmp
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
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  ! this finishes the integral over the BZ  (k)
  !
  CALL mp_sum(sigmar_all,inter_pool_comm)
  CALL mp_sum(sigmai_all,inter_pool_comm)
  CALL mp_sum(zi_all,inter_pool_comm)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
#endif  
  !
  IF (lgamma) THEN
     ikk = ik
     ikq = ik
  ELSE
     ikk = 2 * ik - 1
     ikq = ikk + 1
  ENDIF  
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

    WRITE(linewidth_elself,'(i9,2x)',advance='no') ik
    WRITE(linewidth_elself,'(i9,2x)',advance='no') ibndmin-1+ibnd
    WRITE(linewidth_elself,'(E22.14,2x)',advance='no') ryd2ev * ekk
    WRITE(linewidth_elself,'(E22.14,2x)') ryd2mev*sigmai_all(ibnd,ik)
    !
  ENDDO
  WRITE(stdout,'(5x,a/)') repeat('-',67)
  !
  IF ( ik .eq. (nkqtotf - nqtotf)) THEN
    CLOSE(linewidth_elself)
    IF ( ALLOCATED(sigmar_all) )   DEALLOCATE( sigmar_all )
    IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
    IF ( ALLOCATED(zi_all) )       DEALLOCATE( zi_all )
  ENDIF
  !
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  102 FORMAT(5x,'E( ',i3,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6,' lam=',f15.6)
  103 FORMAT(5x,'k( ',i7,' )=',f9.4,' eV   Re[Sigma]=',f15.6,' meV Im[Sigma]=',f15.6,' meV     Z=',f15.6)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_k






