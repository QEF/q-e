  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine nesting_fn_q (iq )
  !-----------------------------------------------------------------------
  !!
  !!  compute the imaginary part of the phonon self energy due to electron-
  !!  phonon interaction in the Migdal approximation. This corresponds to 
  !!  the phonon linewidth (half width). The phonon frequency is taken into
  !!  account in the energy selection rule.
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation. 
  !!
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw,     &
                        nsmear, delta_smear, efermi_read, fermi_energy
  USE pwcom,     ONLY : nelec, ef, isk
  USE elph2,     ONLY : ibndmax, ibndmin, etf, &
                        wkf, xqf, wqf, nkqf, &
                        nkf, nkqtotf, xqf
  USE constants_epw, ONLY : ryd2ev, two, pi
#if defined(__NAG)
  USE f90_unix_io,  ONLY : flush
#endif
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_global, ONLY : inter_pool_comm
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
  !! Upper bounds index after k or q paral
  !! Smearing for the Gaussian function 
  ! 
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: weight
  !! Imaginary part of the phonhon self-energy factor 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  real(kind=DP) :: w0gauss, dos_ef, gamma, degaussw0
  real(kind=DP), external :: efermig
  !
  !
  IF (iq.eq.1) then 
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout,'(5x,"Nesting Function in the double delta approx")')
     WRITE(stdout,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
  ENDIF
  !
  ! SP: The Gamma function needs to be put to 0 for each q
  gamma = 0.0
  ! 
  ! Here we loop on smearing values
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1)*delta_smear+degaussw
    !
    ! Fermi level and corresponding DOS
    !
    !   Note that the weights of k+q points must be set to zero here
    !   no spin-polarized calculation here
    IF ( efermi_read ) THEN
      ef0 = fermi_energy 
    ELSE
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk)
    ENDIF
    !
    dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
    !  N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    IF (iq.eq.1) then
      WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
      WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
    ENDIF
    !
    !
    CALL start_clock('nesting')
    !
    fermicount = 0
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
        !
        fermicount = fermicount + 1
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
          !
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ekq = etf (ibndmin-1+jbnd, ikq) - ef0
            w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
            !
            ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
            ! This is the imaginary part of the phonon self-energy, sans the matrix elements
            !
            ! weight = wkf (ikk) * (wgkk - wgkq) * &
            !      aimag ( cone / ( ekq - ekk  - ci * degaussw ) ) 
            !
            ! the below expression is positive-definite, but also an approximation
            ! which neglects some fine features
            !
            weight = wkf (ikk) * w0g1 * w0g2
            !
            gamma  =   gamma  + weight  
            !
          ENDDO ! jbnd
        ENDDO   ! ibnd
        !
      ENDIF ! endif fsthick
      !
    ENDDO ! loop on k
    !
    ! collect contributions from all pools (sum over k-points)
    ! this finishes the integral over the BZ  (k)
    !
    CALL mp_sum(gamma,inter_pool_comm) 
    CALL mp_sum(fermicount, inter_pool_comm)
    !
    WRITE(stdout,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
    WRITE(stdout,'(5x,a)') repeat('-',67)
       ! 
    WRITE(stdout, 102)  gamma
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    CALL flush(6)
    !
    WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
      'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
    !
    !
    CALL stop_clock('nesting')
  ENDDO !smears
  !
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' eV, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,' Nesting function (q)=',e15.6,' [Adimensional]')
  !
  end subroutine nesting_fn_q
  !
  !-----------------------------------------------------------------------
  subroutine nesting_fn_k ( ik )
  !-----------------------------------------------------------------------
  !!
  !!  Compute the imaginary part of the phonon self energy due to electron-
  !!  phonon interaction in the Migdal approximation. This corresponds to 
  !!  the phonon linewidth (half width). The phonon frequency is taken into
  !!  account in the energy selection rule.
  !!
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw,     &
                        nsmear, delta_smear, efermi_read, fermi_energy
  USE pwcom,     ONLY : nelec, ef, isk
  USE elph2,     ONLY : ibndmax, ibndmin, etf, etf_k, &
                        wkf, xqf, wqf, nkqf, nqf, nqtotf, &
                        nkqtotf, xqf, gamma_nest
  USE constants_epw, ONLY : ryd2ev, two, pi, zero
#if defined(__NAG)
  USE f90_unix_io,  ONLY : flush
#endif
  USE mp,        ONLY : mp_barrier,mp_sum, mp_bcast
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  implicit none
  !
  INTEGER, INTENT (in) :: ik
  !! Current k-point
  !
  ! Work variables
  INTEGER :: iq
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
  INTEGER :: ismear
  !! Smearing for the Gaussian function 
  INTEGER :: lower_bnd
  !! Upper bounds index after k or q paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k or q paral
  ! 
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: weight
  !! Imaginary part of the phonhon self-energy factor 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  ! 
  REAL(kind=DP) :: w0gauss, degaussw0
  REAL(kind=DP), external :: efermig_seq, dos_ef_seq
  REAL(kind=DP), ALLOCATABLE :: xqf_all(:,:), wqf_all(:,:)
  REAL(kind=DP), external :: efermig
  !
  !
  IF (ik.eq.1) then 
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Nesting Function in the double delta approx")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick.lt.1.d3 ) &
         WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
    IF ( .not. ALLOCATED (gamma_nest) )    ALLOCATE( gamma_nest (nqtotf,nsmear) )
    gamma_nest(:,:)   = zero
  ENDIF

! here we loop on smearing values
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1)*delta_smear+degaussw
    !
    ! Fermi level and corresponding DOS
    !
    !   Note that the weights of k+q points must be set to zero here
    !   no spin-polarized calculation here
    IF ( efermi_read ) THEN
       ef0 = fermi_energy 
    ELSE
       IF (mpime .eq. ionode_id) THEN
         ef0 = efermig_seq(etf_k, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk)
       ENDIF
       CALL mp_bcast (ef0, ionode_id, inter_pool_comm)
    ENDIF
    !
    IF (mpime .eq. ionode_id) THEN
      dosef = dos_ef_seq (ngaussw, degaussw0, ef0, etf_k, wkf, nkqf, nbndsub)
      !   N(Ef) in the equation for lambda is the DOS per spin
      dosef = dosef / two
    ENDIF
    CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
    !
    IF (ik.eq.1) then
      WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
      WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
    ENDIF
    !
    ! Find the bounds of q-dependent arrays in the parallel case in each pool
    CALL fkbounds( nqtotf, lower_bnd, upper_bnd )
    !
    CALL start_clock('nesting')
    !
    fermicount = 0
    !
    DO iq = 1, nqf
      !
      ikq = 2 * iq
      ikk = ikq - 1
      ! 
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
        !
        fermicount = fermicount + 1
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
          !
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ekq = etf (ibndmin-1+jbnd, ikq) - ef0
            w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
            !
            ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
            ! This is the imaginary part of the phonon self-energy, sans the matrix elements
            !
            ! weight = wkf (ikk) * (wgkk - wgkq) * &
            !      aimag ( cone / ( ekq - ekk  - ci * degaussw ) ) 
            !
            ! the below expression is positive-definite, but also an approximation
            ! which neglects some fine features
            !
            weight = wkf (ikk) * w0g1 * w0g2
            !
            gamma_nest(iq+lower_bnd-1,ismear) = gamma_nest(iq+lower_bnd-1,ismear) + weight 
            !
          ENDDO ! jbnd
        ENDDO   ! ibnd
        !
      ENDIF ! endif fsthick
      !
    ENDDO ! loop on q
    !
    CALL stop_clock('nesting')
    !
  ENDDO !smears
  !
  IF ( ik .eq. (nkqtotf - nqtotf)) THEN
    !
    ALLOCATE ( xqf_all (3,nqtotf ),wqf_all (1,nqtotf) )
    xqf_all(:,:) = zero
    wqf_all(:,:) = zero
    ! 
#if defined(__MPI)
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2 ( 3,       nqtotf, nqf, xqf,    xqf_all  )
    CALL poolgather2 ( 1,       nqtotf, nqf, wqf,    wqf_all  )
    CALL mp_sum(gamma_nest, inter_pool_comm )
    CALL mp_sum(fermicount, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    xqf_all = xqf
    DO iq = 1, nqtotf
      wqf_all(1,iq) = wqf(iq)
    ENDDO
    !
#endif
    DO ismear = 1, nsmear
      degaussw0 = (ismear-1)*delta_smear+degaussw
      WRITE(stdout,'(/5x,"smearing = ",f9.5)') degaussw0
      DO iq = 1, nqtotf
        WRITE(stdout,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf_all(:,iq), wqf_all(1,iq)
        WRITE(stdout,'(5x,a)') repeat('-',67)
           ! 
        WRITE(stdout, 102)  gamma_nest(iq,ismear)
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        CALL flush(6)
        !
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
           'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
        !  
      ENDDO
      ! 
    ENDDO
    !
  endif
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' eV, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,'Nesting function (q)=',e15.6,' [Adimensional]')
  !
  end subroutine nesting_fn_k
