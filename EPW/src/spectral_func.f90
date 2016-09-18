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
  !!
  !!  Compute the electron spectral function including the  electron-
  !!  phonon interaction in the Migdal approximation. 
  !!  
  !!  We take the trace of the spectral function to simulate the photoemission
  !!  intensity. I do not consider the c-axis average for the time being.
  !!  The main approximation is constant dipole matrix element and diagonal
  !!  selfenergy. The diagonality can be checked numerically. 
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup ,iospectral
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, eps_acustic, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, shortrange, &
                            efermi_read, fermi_energy
  USE pwcom,         ONLY : nelec, ef, isk
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            epf17, wkf, nkf, nqtotf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm
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
  INTEGER :: iw
  !! Counter on the frequency
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
  ! 
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
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: inv_degaussw
  !! Inverse of the smearing for efficiency reasons  
  REAL(kind=DP) :: ww
  !! Current frequency
  REAL(kind=DP) :: dw 
  !! Frequency intervals
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerenc  
  real(kind=DP) :: specfun_sum, esigmar0
  real(kind=DP) :: fermi(nw_specfun)
  real(kind=DP), external :: efermig, dos_ef, wgauss
  !
  ! variables for collecting data from all pools in parallel case 
  !
  real(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:)
  ! 
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
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
    ikk = 2 * ik - 1
    ikq = ikk + 1
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
            wgkq = wgauss( -ekq/eptemp, -99)  
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
#if defined(__MPI)
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
    IF (me_pool == 0) then
      OPEN(unit=iospectral,file='specfun.elself') 
      OPEN(unit=iospectral_sup,file='specfun_sup.elself') 
    ENDIF
    IF (me_pool == 0) then
      WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
      WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)' 
    ENDIF
    IF (me_pool == 0) then
      WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
      WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
&         Real Sigma[meV]  Im Sigma[meV]'
    ENDIF
    !
    DO ik = 1, nksqtotf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
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
      IF (me_pool == 0) &
        WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
         !
      ENDDO
      !
      IF (me_pool == 0) &
        WRITE(iospectral,'(a)') ' '
      IF (me_pool == 0) &
        WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
      !
    ENDDO
    !
    IF (me_pool == 0)  CLOSE(iospectral)
    !
    DO ibnd = 1, ibndmax-ibndmin+1
      !
      DO ik = 1, nksqtotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
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
          IF (me_pool == 0) &
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
    IF (me_pool == 0)  CLOSE(iospectral_sup)
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
  !!
  !!  Compute the electron spectral function including the  electron-
  !!  phonon interaction in the Migdal approximation. 
  !!  
  !!  We take the trace of the spectral function to simulate the photoemission
  !!  intensity. I do not consider the c-axis average for the time being.
  !!  The main approximation is constant dipole matrix element and diagonal
  !!  selfenergy. The diagonality can be checked numerically. 
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation
  !!
  !!-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup, iospectral
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, eps_acustic, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, shortrange, &
                            efermi_read, fermi_energy
  USE pwcom,         ONLY : ef
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, nqf, etf_k, xqf, &
                            epf17, wkf, nqtotf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all, esigmai_all, a_all, efnew
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci
  USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_global,     ONLY : me_pool, inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE io_global,     ONLY : ionode_id  
  ! 
  implicit none
  !
  INTEGER, INTENT (in) :: ik
  !! Current k-point
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
  INTEGER :: iw
  !! Counter on the frequency
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
  ! 
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: inv_degaussw
  !! Inverse of the smearing for efficiency reasons   
  REAL(kind=DP), PARAMETER :: eps2 = 0.01/ryd2mev
  !! Tolerenc  
  real(kind=DP), external :: efermig, dos_ef, wgauss
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, ww, dw, weight
  real(kind=DP) :: dosef, specfun_sum, esigmar0
  real(kind=DP) :: fermi(nw_specfun)
  REAL(kind=DP), external ::  dos_ef_seq
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  ! 
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  ! 
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
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
  IF (mpime .eq. ionode_id) THEN
    !
    dosef = dos_ef_seq (ngaussw, degaussw, ef0, etf_k, wkf, nkqf, nbndsub)/2
  ENDIF
  CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
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
    ikq = 2 * iq
    ikk = ikq - 1
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
            wgkq = wgauss( -ekq/eptemp, -99)  
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
  ! Collect contributions from all pools (sum over q-points)
  ! this finishes the integral over the BZ  (q)
  !
  CALL mp_sum(esigmar_all,inter_pool_comm)
  CALL mp_sum(esigmai_all,inter_pool_comm)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  IF (.not. ALLOCATED (a_all)) ALLOCATE ( a_all(nw_specfun, nksqtotf) )
  a_all(:,:) = zero  
  !
  ! Output electron spectral function here after looping over all q-points 
  ! (with their contributions summed in a etc.)
  IF ( ik .eq. 1 ) THEN
    IF (me_pool == 0) then
      !
      OPEN(unit=iospectral,file='specfun.elself')
      OPEN(unit=iospectral_sup,file='specfun_sup.elself')
      WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
      WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
      WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
      WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
&         Real Sigma[meV]  Im Sigma[meV]'
      !
    ENDIF
  ENDIF
  !
  WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
  !
  ikk = 2 * ik - 1
  ikq = ikk + 1
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
    IF (me_pool == 0) &
      WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ik) / ryd2mev
      !
  ENDDO
  !
  IF (me_pool == 0) &
    WRITE(iospectral,'(a)') ' '
  IF (me_pool == 0) &
    WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
    !
  DO ibnd = 1, ibndmax-ibndmin+1
    !
    ikk = 2 * ik - 1
    ikq = ikk + 1
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
      IF (me_pool == 0) &
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
    IF (me_pool == 0) THEN
      CLOSE(iospectral_sup)
      CLOSE(iospectral)
    ENDIF
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
