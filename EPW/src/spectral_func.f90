  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_q ( iq )
  !-----------------------------------------------------------------------
  !
  !  Compute the electron spectral function including the  electron-
  !  phonon interaction in the Migdal approximation. 
  !  
  !  We take the trace of the spectral function to simulate the photoemission
  !  intensity. I do not consider the c-axis average for the time being.
  !  The main approximation is constant dipole matrix element and diagonal
  !  selfenergy. The diagonality can be checked numerically. 
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  01/2014 Modified by Roxana Margine 
  !
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iunepmatf, iospectral_sup ,iospectral
  USE phcom,         ONLY : nmodes
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, etf_mem, eps_acustic, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, &
                            efermi_read, fermi_energy
  USE pwcom,         ONLY : nelec, ef, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  real(kind=DP), external :: efermig, dos_ef, wgauss
  integer :: iw, ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, ww, dw, weight
  real(kind=DP) :: dosef, specfun_sum, esigmar0
  real(kind=DP) :: fermi(nw_specfun)
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf, lower_bnd, upper_bnd
  real(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:)
  ! 
  ! energy range and spacing for spectral function
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  IF ( iq .eq. 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron Spectral Function in the Migdal Approximation")')
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
     ef0 = efermig(etf,nbndsub,nkqf,nelec,wkf,degaussw,ngaussw,0,isk)
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
     IF ( .not. ALLOCATED(esigmar_all) ) ALLOCATE( esigmar_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
     IF ( .not. ALLOCATED(esigmai_all) ) ALLOCATE( esigmai_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
     esigmar_all(:,:,:) = zero
     esigmai_all(:,:,:) = zero
  ENDIF 
  !
  ! SP: Sum rule added to conserve the number of electron. 
  IF ( iq .eq. 1 ) THEN
    WRITE (stdout,'(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
    WRITE (stdout,'(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
    WRITE (stdout,'(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
    WRITE (stdout,'(a)') ' '
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
     ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
     ! (but in this case they are the same)
     !
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
           !
           ! the phonon frequency and Bose occupation
           wq = wf (imode, iq)
           wgq = wgauss( -wq/eptemp, -99)
           wgq = wgq / ( one - two * wgq )
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
                 wgkq = wgauss( -ekq/eptemp, -99)  
                 !
                 ! here we take into account the zero-point sqrt(hbar/2M\omega)
                 ! with hbar = 1 and M already contained in the eigenmodes
                 ! g2 is Ry^2, wkf must already account for the spin factor
                 !
                 IF (wq .gt. eps_acustic) THEN
                    g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                 ELSE
                    g2 = 0.d0
                 ENDIF
                 !
                 DO iw = 1, nw_specfun
                    !
                    ww = wmin_specfun + dble (iw-1) * dw
                    !
                    weight = wqf(iq) * real (                                            &
                      ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
                    !
                    esigmar_all(ibnd,ik+lower_bnd-1,iw) = esigmar_all(ibnd,ik+lower_bnd-1,iw) + g2 * weight 
                    ! 
                    ! SP : Application of the sum rule
                    esigmar0 =  g2 *  wqf(iq) * real (                                   &
                      ( (       wgkq + wgq ) / ( -( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( -( ekq + wq ) - ci * degaussw ) ) )
                    esigmar_all(ibnd,ik+lower_bnd-1,iw)=esigmar_all(ibnd,ik+lower_bnd-1,iw)-esigmar0
                    !
                    weight = wqf(iq) * aimag (                                           &
                      ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
                    !
                    esigmai_all(ibnd,ik+lower_bnd-1,iw) = esigmai_all(ibnd,ik+lower_bnd-1,iw) + g2 * weight
                    !
                 ENDDO
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
     CALL mp_sum( esigmar_all, inter_pool_comm )
     CALL mp_sum( esigmai_all, inter_pool_comm )
     CALL mp_sum( fermicount, inter_pool_comm )
     CALL mp_barrier(inter_pool_comm)
     !
#else
     !
     xkf_all = xkf
     etf_all = etf
     !
#endif
     !
     ! Output electron spectral function here after looping over all q-points 
     ! (with their contributions summed in a etc.)
     !
     WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     !
     ! construct the trace of the spectral function (assume diagonal selfenergy
     ! and constant matrix elements for dipole transitions)
     !
     IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun, nksqtotf) )
     a_all(:,:) = zero
     !
#ifdef __PARA
     IF (me_pool == 0) then
#endif      
     OPEN(unit=iospectral,file='specfun.elself') 
     OPEN(unit=iospectral_sup,file='specfun_sup.elself') 
#ifdef __PARA
     ENDIF
     IF (me_pool == 0) then
#endif
     WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
     WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)' 
#ifdef __PARA
     ENDIF
     IF (me_pool == 0) then
#endif
     WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
     WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
&          Real Sigma[meV]  Im Sigma[meV]'
#ifdef __PARA
     ENDIF
#endif

     !
     DO ik = 1, nksqtotf
        !
        IF (lgamma) then
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        !
        WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf_all (:,ikk)
        WRITE(stdout,'(5x,a)') repeat('-',67)
        !
        DO iw = 1, nw_specfun
           !
           ww = wmin_specfun + dble (iw-1) * dw
           !
           DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  the energy of the electron at k
              ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
              !
              a_all(iw,ik) = a_all(iw,ik) + abs( esigmai_all(ibnd,ik,iw) ) / pi / &
                 ( ( ww - ekk - esigmar_all(ibnd,ik,iw) )**two + (esigmai_all(ibnd,ik,iw) )**two )
              !
           ENDDO
           !
           WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
        !
        ENDDO
        !
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        !
     ENDDO
     !
     DO ik = 1, nksqtotf
        !
        ! The spectral function should integrate to 1 for each k-point
        specfun_sum = 0.0
        ! 
        DO iw = 1, nw_specfun
           !
           ww = wmin_specfun + dble (iw-1) * dw
           fermi(iw) = wgauss(-ww/eptemp, -99) 
           WRITE(stdout,'(2x,i7,2x,f12.4,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
           !
           specfun_sum = specfun_sum + a_all(iw,ik)* fermi(iw) * dw !/ ryd2mev
           !
#ifdef __PARA
        IF (me_pool == 0) &
#endif
        WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
           !
        ENDDO
        !
#ifdef __PARA
        IF (me_pool == 0) &
#endif
        WRITE(iospectral,'(a)') ' '
#ifdef __PARA
        IF (me_pool == 0) &
#endif
        WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
        !
     ENDDO
     !
#ifdef __PARA
     IF (me_pool == 0) &
#endif      
     CLOSE(iospectral)
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
           !  the energy of the electron at k
           ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
           !
           DO iw = 1, nw_specfun
              !
              ww = wmin_specfun + dble (iw-1) * dw
              WRITE(stdout,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
                ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
                ryd2mev * esigmai_all(ibnd,ik,iw)
              ! 
#ifdef __PARA
              IF (me_pool == 0) &
#endif
              WRITE(iospectral_sup,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
                ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
                ryd2mev * esigmai_all(ibnd,ik,iw)
              !
           ENDDO
           !
        ENDDO
        !
        WRITE(stdout,*) ' '
        !
     ENDDO
     !
#ifdef __PARA
     IF (me_pool == 0) &
#endif      
     CLOSE(iospectral_sup)
     !
     IF ( ALLOCATED(xkf_all) )      DEALLOCATE( xkf_all )
     IF ( ALLOCATED(etf_all) )      DEALLOCATE( etf_all )
     IF ( ALLOCATED(esigmar_all) )  DEALLOCATE( esigmar_all )
     IF ( ALLOCATED(esigmai_all) )  DEALLOCATE( esigmai_all )
     IF ( ALLOCATED(a_all) )        DEALLOCATE( a_all )
     !
  ENDIF
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  103 FORMAT(5x,'ik = ',i7,'  w = ',f9.4,' eV   A(k,w) = ',e12.5,' meV^-1')
  !
  RETURN
  !
  END SUBROUTINE spectral_func_q
  !
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_k ( ik )
  !-----------------------------------------------------------------------
  !
  !  Compute the electron spectral function including the  electron-
  !  phonon interaction in the Migdal approximation. 
  !  
  !  We take the trace of the spectral function to simulate the photoemission
  !  intensity. I do not consider the c-axis average for the time being.
  !  The main approximation is constant dipole matrix element and diagonal
  !  selfenergy. The diagonality can be checked numerically. 
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  01/2014 Modified by Roxana Margine 
  !
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iunepmatf, iospectral_sup, iospectral
  USE phcom,         ONLY : nmodes
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, etf_mem, eps_acustic, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, &
                            efermi_read, fermi_energy
  USE pwcom,         ONLY : nelec, ef, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, nqf, etf_k, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all, efnew
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id  
#endif
  implicit none
  !
  real(kind=DP), external :: efermig, dos_ef, wgauss
  integer :: iw, ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, ww, dw, weight
  real(kind=DP) :: dosef, eptemp0, specfun_sum, esigmar0
  real(kind=DP) :: fermi(nw_specfun)
   REAL(kind=DP), external ::  dos_ef_seq
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf, lower_bnd, upper_bnd
  ! 
  ! energy range and spacing for spectral function
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  IF ( ik .eq. 1 ) THEN
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Electron Spectral Function in the Migdal Approximation")')
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
     !
  ENDIF
  !
#ifdef __PARA
  IF (mpime .eq. ionode_id) THEN
#endif
    !
    dosef = dos_ef_seq (ngaussw, degaussw, ef0, etf_k, wkf, nkqf, nbndsub)/2
#ifdef __PARA
  ENDIF
  CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
#endif
  !
  IF ( ik .eq. 1 ) THEN 
     WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
     WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
     WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  IF ( ik .eq. 1 ) THEN 
     IF ( .not. ALLOCATED(esigmar_all) ) ALLOCATE( esigmar_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
     IF ( .not. ALLOCATED(esigmai_all) ) ALLOCATE( esigmai_all(ibndmax-ibndmin+1, nksqtotf, nw_specfun) )
     esigmar_all(:,:,:) = zero
     esigmai_all(:,:,:) = zero
    !
    ! SP: Sum rule added to conserve the number of electron. 
    !
    WRITE (stdout,'(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
    WRITE (stdout,'(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
    WRITE (stdout,'(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
    WRITE (stdout,'(a)') ' '
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
     ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
     ! (but in this case they are the same)
     !
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
           !
           ! the phonon frequency and Bose occupation
           wq = wf (imode, iq)
           wgq = wgauss( -wq/eptemp, -99)
           wgq = wgq / ( one - two * wgq )
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
                 wgkq = wgauss( -ekq/eptemp, -99)  
                 !
                 ! here we take into account the zero-point sqrt(hbar/2M\omega)
                 ! with hbar = 1 and M already contained in the eigenmodes
                 ! g2 is Ry^2, wkf must already account for the spin factor
                 !
                 IF (wq .gt. eps_acustic) THEN
                    g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                 ELSE
                    g2 = 0.d0
                 ENDIF
                 !
                 DO iw = 1, nw_specfun
                    !
                    ww = wmin_specfun + dble (iw-1) * dw
                    !
                    weight = wqf(iq) * real (                                            &
                      ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
                    !
                    esigmar_all(ibnd,ik,iw) = esigmar_all(ibnd,ik,iw) + g2 * weight 
                    ! 
                    ! SP : Application of the sum rule
                    esigmar0 =  g2 *  wqf(iq) * real (                                   &
                      ( (       wgkq + wgq ) / ( -( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( -( ekq + wq ) - ci * degaussw ) ) )
                    esigmar_all(ibnd,ik,iw)=esigmar_all(ibnd,ik,iw)-esigmar0
                    !
                    weight = wqf(iq) * aimag (                                           &
                      ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                        ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
                    !
                    esigmai_all(ibnd,ik,iw) = esigmai_all(ibnd,ik,iw) + g2 * weight
                    !
                 ENDDO
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
  ! collect contributions from all pools (sum over q-points)
  ! this finishes the integral over the BZ  (q)
  !
  CALL mp_sum(esigmar_all,inter_pool_comm)
  CALL mp_sum(esigmai_all,inter_pool_comm)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
#endif 
  !
  IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun, nksqtotf) )
  a_all(:,:) = zero  
  !
  ! Output electron spectral function here after looping over all q-points 
  ! (with their contributions summed in a etc.)
  IF ( ik .eq. 1 ) THEN
#ifdef __PARA
    IF (me_pool == 0) then
#endif      
      !
      OPEN(unit=iospectral,file='specfun.elself')
      OPEN(unit=iospectral_sup,file='specfun_sup.elself')
      WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
      WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
      WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
      WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
&         Real Sigma[meV]  Im Sigma[meV]'
      !
#ifdef __PARA
    ENDIF
#endif
  ENDIF
  !
  WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
  !
  IF (lgamma) then
     ikk = ik
     ikq = ik
  ELSE
     ikk = 2 * ik - 1
     ikq = ikk + 1
  ENDIF
  !
  WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf(:,ikk)
  WRITE(stdout,'(5x,a)') repeat('-',67)
  !
  DO iw = 1, nw_specfun
     !
     ww = wmin_specfun + dble (iw-1) * dw
     !
     DO ibnd = 1, ibndmax-ibndmin+1
        !
        !  the energy of the electron at k
        ekk = etf_k (ibndmin-1+ibnd, ikk) - ef0
        !
        a_all(iw,ik) = a_all(iw,ik) + abs( esigmai_all(ibnd,ik,iw) ) / pi / &
           ( ( ww - ekk - esigmar_all(ibnd,ik,iw) )**two + (esigmai_all(ibnd,ik,iw) )**two )
        !
     ENDDO
     !
     WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
     !
  ENDDO
  !
  WRITE(stdout,'(5x,a/)') repeat('-',67)
  !
  !
  ! The spectral function should integrate to 1 for each k-point
  specfun_sum = 0.0
  ! 
  DO iw = 1, nw_specfun
    !
    ww = wmin_specfun + dble (iw-1) * dw
    fermi(iw) = wgauss(-ww/eptemp, -99) 
    WRITE(stdout,'(2x,i7,2x,f12.4,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
    !
    specfun_sum = specfun_sum + a_all(iw,ik)* fermi(iw) * dw !/ ryd2mev
    !
#ifdef __PARA
    IF (me_pool == 0) &
#endif
      WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
      !
  ENDDO
  !
#ifdef __PARA
  IF (me_pool == 0) &
#endif
    WRITE(iospectral,'(a)') ' '
#ifdef __PARA
  IF (me_pool == 0) &
#endif
    WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
    !
  DO ibnd = 1, ibndmax-ibndmin+1
    !
    !
    IF (lgamma) THEN
      ikk = ik
      ikq = ik
    ELSE
      ikk = 2 * ik - 1
      ikq = ikk + 1
    ENDIF
    !
    !  the energy of the electron at k
    ekk = etf_k (ibndmin-1+ibnd, ikk) - ef0
    !
    DO iw = 1, nw_specfun
      !
      ww = wmin_specfun + dble (iw-1) * dw
      WRITE(stdout,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
        ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
        ryd2mev * esigmai_all(ibnd,ik,iw)
      ! 
#ifdef __PARA
      IF (me_pool == 0) &
#endif
        WRITE(iospectral_sup,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
          ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all(ibnd,ik,iw),&
          ryd2mev * esigmai_all(ibnd,ik,iw)
        !
    ENDDO
    !
    !
    WRITE(stdout,*) ' '
    !
  ENDDO ! ibnd
  !
  IF ( ik .eq. (nkqtotf - nqtotf)) THEN
    IF ( ALLOCATED(esigmar_all) )  DEALLOCATE( esigmar_all )
    IF ( ALLOCATED(esigmai_all) )  DEALLOCATE( esigmai_all )
    IF ( ALLOCATED(a_all) )        DEALLOCATE( a_all )
    !
#ifdef __PARA
    IF (me_pool == 0) THEN
#endif      
      CLOSE(iospectral_sup)
      CLOSE(iospectral)
#ifdef __PARA
    ENDIF
#endif    
    !
  ENDIF
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
  103 FORMAT(5x,'ik = ',i7,'  w = ',f9.4,' eV   A(k,w) = ',e12.5,' meV^-1')
  !
  RETURN
  !
  END SUBROUTINE spectral_func_k


