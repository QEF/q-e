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
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation. 
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     only : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iuetf
  use epwcom,    only : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw,     &
                        etf_mem, &
                        nsmear, delta_smear, efermi_read, fermi_energy
  use pwcom,     only : nelec, ef, isk
  use elph2,     only : ibndmax, ibndmin, etf, &
                        wkf, xqf, wqf, nkqf, &
                        nkf, nkqtotf, xqf
  USE constants_epw, ONLY : ryd2ev, two, pi
#ifdef __NAG
  USE f90_unix_io,  ONLY : flush
#endif
#ifdef __PARA
  use mp,        only : mp_barrier,mp_sum
  use mp_global, only : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, nrec, iq, fermicount, ismear
  real(kind=DP) :: ekk, ekq, ef0, &
     weight, w0g1, w0g2, w0gauss, wgauss, dosef, dos_ef, gamma, &
     degaussw0, eptemp0
  !
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
          'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
  ENDIF
  !
  ! SP: The Gamma function needs to be put to 0 for each q
  gamma = 0.0
  ! 
  ! Here we loop on smearing values
  DO ismear = 1, nsmear
  !
  degaussw0 = (ismear-1)*delta_smear+degaussw
  eptemp0 = (ismear-1)*delta_smear+eptemp(1)
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
  !   N(Ef) in the equation for lambda is the DOS per spin
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
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
     !
     ! when we see references to iq for file readins, it is always = 1 
     IF (.not. etf_mem) then
        nrec = ikk
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = ikq
        CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
     ENDIF
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
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  ! this finishes the integral over the BZ  (k)
  !
  CALL mp_sum(gamma,inter_pool_comm) 
  CALL mp_sum(fermicount, inter_pool_comm)
  !
#endif
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
  !
  !  Compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     only : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iuetf
  use epwcom,    only : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw,     &
                        etf_mem, &
                        nsmear, delta_smear, efermi_read, fermi_energy
  use pwcom,     only : nelec, ef, isk
  use elph2,     only : ibndmax, ibndmin, etf, etf_k, &
                        wkf, xqf, wqf, nkqf, nqf, nqtotf, &
                        nkf, nkqtotf, xqf, gamma_nest
  USE constants_epw, ONLY : ryd2ev, two, pi, zero
#ifdef __NAG
  USE f90_unix_io,  ONLY : flush
#endif
#ifdef __PARA
  use mp,        only : mp_barrier,mp_sum, mp_bcast
  use mp_global, only : me_pool,inter_pool_comm,my_pool_id
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, nrec, iq, fermicount, ismear, &
             lower_bnd, upper_bnd
  real(kind=DP) :: ekk, ekq, ef0, &
     weight, w0g1, w0g2, w0gauss, wgauss, dosef, dos_ef,  &
     degaussw0, eptemp0
  real(kind=DP), external :: efermig_seq, dos_ef_seq
  REAL(kind=DP), ALLOCATABLE :: xqf_all(:,:), wqf_all(:,:)
  !
  real(kind=DP), external :: efermig
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
          'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
     IF ( .not. ALLOCATED (gamma_nest) )    ALLOCATE( gamma_nest (nqtotf,nsmear) )
     gamma_nest(:,:)   = zero
  ENDIF

! here we loop on smearing values
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1)*delta_smear+degaussw
    eptemp0 = (ismear-1)*delta_smear+eptemp(1)
    !
    ! Fermi level and corresponding DOS
    !
    !   Note that the weights of k+q points must be set to zero here
    !   no spin-polarized calculation here
    IF ( efermi_read ) THEN
       ef0 = fermi_energy 
    ELSE
#ifdef __PARA
       IF (mpime .eq. ionode_id) THEN
#endif
       ef0 = efermig_seq(etf_k, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk)
#ifdef __PARA
       ENDIF
       CALL mp_bcast (ef0, ionode_id, inter_pool_comm)
#endif 
    ENDIF
    !
#ifdef __PARA
    IF (mpime .eq. ionode_id) THEN
#endif
      dosef = dos_ef_seq (ngaussw, degaussw0, ef0, etf_k, wkf, nkqf, nbndsub)
      !   N(Ef) in the equation for lambda is the DOS per spin
      dosef = dosef / two
#ifdef __PARA
    ENDIF
    CALL mp_bcast (dosef, ionode_id, inter_pool_comm)
#endif
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
       ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
       !
       ! when we see references to iq for file readins, it is always = 1 
       IF (.not. etf_mem) then
          nrec = ikk
          CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
          nrec = ikq
          CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
       ENDIF
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
#ifdef __PARA
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
    wqf_all = wqf
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
