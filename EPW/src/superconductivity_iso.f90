  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE superconductivity_iso
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutines linked with superconductivity 
  !! using the isotropic Eliashberg formalism. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_iaxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic 
    !! Eliashberg equations on the imaginary-axis.  
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                              limag, lpade, lacon
    USE eliashbergcom, ONLY : nsw, nsiw, Deltai, Deltaip, Delta, Deltap, estemp
    USE constants_epw, ONLY : kelvin2eV, ci
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE superconductivity, ONLY : free_energy, dos_quasiparticle, gen_freqgrid_iaxis, & 
                                  deallocate_eliashberg_iso_iaxis, deallocate_eliashberg
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: N
    !! Maximum nr. frequency points in Pade approx
    REAL(DP) :: rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(DP) :: tcpu
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv
    !! True if calculation is converged
    !
    CALL start_clock( 'iso_iaxis' )
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a,i3,a,f12.5,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a)') 'Solve isotropic Eliashberg equations on imaginary-axis' 
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a,i6,a,i6)') 'Total number of frequency points nsiw ( ', itemp, ' ) = ', nsiw(itemp)
      WRITE(stdout,'(a)') '    '
      CALL start_clock( 'iaxis_imag' )
      CALL gen_freqgrid_iaxis( itemp )
      !
      IF ( limag ) THEN
        iter = 1
        conv = .false.
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL sum_eliashberg_iso_iaxis( itemp, iter, conv )
          CALL mix_broyden( nsiw(itemp), Deltai, Deltaip, broyden_beta, iter, broyden_ndim, conv )
          iter = iter + 1
        ENDDO ! iter
        !
        IF ( conv ) THEN
          !
          ! SP : Only print the Free energy if the user want it
          !
          IF ( iverbosity .eq. 2 ) THEN
            CALL free_energy( itemp )
          ENDIF
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'iaxis_imag' )
          CALL print_clock( 'iaxis_imag' )
        ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') 'not converged  '
          CALL stop_clock( 'iaxis_imag' )
          CALL print_clock( 'iaxis_imag' )
          CALL errore('sum_eliashberg_iso_iaxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF
      !
      IF ( lpade ) THEN 
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a)') 'Pade approximant of isotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '    '
        CALL start_clock( 'raxis_pade' )
        !
        iter = 1
        conv = .false.
        N = 80 * nsiw(itemp) / 100
        IF ( mod(N,2) .ne. 0 ) N = N + 1
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL pade_cont_iso_iaxis_to_raxis( itemp, N, conv )
          N = N - 2
          iter = iter + 1
        ENDDO
        !
        IF ( conv ) THEN
          CALL dos_quasiparticle( itemp )
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_pade' )
          CALL print_clock( 'raxis_pade' )
          WRITE(stdout,'(a)') '  '
        ELSEIF ( .not. conv  .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_pade' )
          CALL print_clock( 'raxis_pade' )
          CALL errore('pade_cont_iso_iaxis_to_raxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF 
      !
      IF ( lacon ) THEN 
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a)') 'Analytic continuation of isotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a,i6)') 'Total number of frequency points nsw = ', nsw
        WRITE(stdout,'(a)') '    '
        CALL start_clock( 'raxis_acon' )
        !
        iter = 1
        conv = .false.
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL analytic_cont_iso_iaxis_to_raxis( itemp, iter, conv )
          rdeltain(:)  = real(Deltap(:))
          cdeltain(:)  = aimag(Deltap(:))
          rdeltaout(:) = real(Delta(:))
          cdeltaout(:) = aimag(Delta(:))
          CALL mix_broyden ( nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
          CALL mix_broyden2( nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
          Deltap(:) = rdeltain(:) + ci * cdeltain(:)
          iter = iter + 1
        ENDDO ! iter
        !
        IF ( conv ) THEN
          CALL dos_quasiparticle( itemp )
          WRITE(stdout,'(a)') ' '
          CALL stop_clock( 'raxis_acon' )
          CALL print_clock( 'raxis_acon' )
          WRITE(stdout,'(a)') ' '
        ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_acon' )
          CALL print_clock( 'raxis_acon' )
          CALL errore('analytic_cont_iso_iaxis_to_raxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF
      !
      CALL deallocate_eliashberg_iso_iaxis
      !
      tcpu = get_clock('iso_iaxis')
      WRITE(stdout,'(5x,a,i3,a,f8.1,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      !
    ENDDO ! itemp
    !
    CALL stop_clock( 'iso_iaxis' )
    !
    RETURN
    !
    END SUBROUTINE eliashberg_iso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_iso_iaxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations on the imaginary-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, wsi, NZnormi, Znormi, Deltai, Deltaip, Keri
    USE constants_epw, ONLY : pi, zero
    USE io_eliashberg, ONLY : eliashberg_write_iaxis
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    ! 
    ! Local variables
    INTEGER :: iw, iwp 
    !! Counter on frequency imag-axis
    REAL(DP) :: kernelp, kernelm, lambdap, lambdam
    !! lambdam - K_{-}(n,n'T))
    !! lambdap - K_{+}(n,n',T) 
    !! kernelm = lambdam - lambdap
    !! kernelp = lambdam + lambdap
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    REAL(DP) :: esqrt
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    !! Temporary variables
    REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !! gap
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsiw(itemp)) )
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(gap) )      ALLOCATE( gap(nstemp) )
      IF ( .not. ALLOCATED(Deltai) )   ALLOCATE( Deltai(nsiw(itemp)) )
      IF ( .not. ALLOCATED(Deltaip) )  ALLOCATE( Deltaip(nsiw(itemp)) )
      IF ( .not. ALLOCATED(Znormi) )   ALLOCATE( Znormi(nsiw(itemp)) )
      IF ( .not. ALLOCATED(NZnormi) )  ALLOCATE( NZnormi(nsiw(itemp)) )
      gap(itemp) = zero
      Deltaip(:) = zero
      Deltaip(:) = gap0
      !
      CALL kernel_iso_iaxis( itemp )
    ENDIF
    Znormi(:) = zero
    NZnormi(:) = zero
    Deltai(:) = zero
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsiw(itemp)) )
      Deltaold(:) = gap0
    ENDIF
    absdelta = zero 
    reldelta = zero 
    DO iw = 1, nsiw(itemp) ! loop over omega
      DO iwp = 1, nsiw(itemp) ! loop over omega_prime
        ! this step is performed at each iter step only for iw=1 since it is independ of wsi(iw)
        IF ( iw .eq. 1 ) THEN
          esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + Deltaip(iwp)**2.d0 )
          wesqrt(iwp) =  wsi(iwp) * esqrt 
          desqrt(iwp) =  Deltaip(iwp) * esqrt 
        ENDIF
        lambdam = Keri( abs(iw-iwp)+1 )
        lambdap = Keri( abs(iw+iwp) )
        kernelm = lambdam - lambdap
        kernelp = lambdam + lambdap
        NZnormi(iw) = NZnormi(iw) + kernelm
        Znormi(iw) = Znormi(iw) + wesqrt(iwp) * kernelm 
        Deltai(iw) = Deltai(iw) + desqrt(iwp) * ( kernelp - 2.d0 * muc ) 
      ENDDO ! iwp
      Znormi(iw) = 1.d0 + pi * estemp(itemp) * Znormi(iw) / wsi(iw)
      NZnormi(iw) = 1.d0 + pi * estemp(itemp) * NZnormi(iw) / wsi(iw)
      Deltai(iw) = pi * estemp(itemp) * Deltai(iw) / Znormi(iw)
      reldelta = reldelta + abs( Deltai(iw) - Deltaold(iw) )
      absdelta = absdelta + abs( Deltai(iw) ) 
    ENDDO ! iw 
    errdelta = reldelta / absdelta
    Deltaold(:) = Deltai(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, '   error = ', errdelta, &
                                           '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
    !
    IF ( errdelta .lt. conv_thr_iaxis ) conv = .true.
    IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
      gap(itemp) = Deltai(1)
      gap0 = gap(itemp)
      CALL eliashberg_write_iaxis( itemp )
    ENDIF
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
      IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
      WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
      WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
      CALL errore('sum_eliashberg_iso_iaxis','increase nsiter or reduce conv_thr_iaxis',1)
      !CALL deallocate_eliashberg
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE sum_eliashberg_iso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_iso_iaxis_to_raxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine does the analyic continuation of the isotropic Eliashberg equations 
    !! from the imaginary-axis to the real axis
    !! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, nsiter, conv_thr_racon, lpade
    USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, gap, a2f_iso, Dsumi, Zsumi, & 
                              Delta, Deltap, Znorm, Znormp, Gp, Gm
    USE constants_epw, ONLY : pi, ci, zero, czero, cone
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT (in) :: iter
    !! Counter on the iteration number
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged
    ! 
    ! Local variables
    INTEGER :: i, iw, iwp 
    !! Counter on frequency real-axis
    REAL(kind=DP) :: rgammap
    !! -bose_einstein( w' ) - fermi_dirac(  w + w' )
    REAL(kind=DP) :: rgammam
    !! bose_einstein( w' ) + fermi_dirac( -w + w' )
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    COMPLEX(kind=DP) :: esqrt, root
    !! Temporary variables
    COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !! gap
    CHARACTER (len=256) :: cname
    !! character in file name
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
      IF ( .not. ALLOCATED(Deltap) )   ALLOCATE( Deltap(nsw) )
      IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
      IF ( .not. ALLOCATED(Znormp) )   ALLOCATE( Znormp(nsw) )
      IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
      Deltap(:) = czero
      Deltaold(:) = czero
      IF ( lpade ) THEN
         Deltap(:) = Delta(:)
         Deltaold(:) = Delta(:)
      ELSE 
         Deltap(:) = gap(itemp)
         Deltaold(:) = gap(itemp)
      ENDIF
      Znormp(:) = cone
      IF ( .not. ALLOCATED(Gp) ) ALLOCATE( Gp(nsw,nqstep) )
      IF ( .not. ALLOCATED(Gm) ) ALLOCATE( Gm(nsw,nqstep) )
      IF ( .not. ALLOCATED(Dsumi) ) ALLOCATE( Dsumi(nsw) )
      IF ( .not. ALLOCATED(Zsumi) ) ALLOCATE( Zsumi(nsw) )
      CALL kernel_iso_iaxis_analytic_cont( itemp )
    ENDIF
    Znorm(:) = czero
    Delta(:) = czero
    !
    absdelta = zero 
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nqstep ! loop over omega_prime
        IF ( iter .eq. 1 ) THEN 
          CALL gamma_acont( ws(iw), ws(iwp), estemp(itemp), rgammap, rgammam )
          Gp(iw,iwp) = rgammap
          Gm(iw,iwp) = rgammam
        ENDIF
        !
        i = iw + iwp - 1
        IF ( i .le. nsw ) THEN
          root = sqrt( Znormp(i)**2.d0 * ( ws(i)**2.d0 - Deltap(i)**2.d0 ) )
          IF ( aimag(root) .lt. zero ) THEN 
             esqrt = Znormp(i) / conjg(root)
          ELSE  
             esqrt = Znormp(i) / root
          ENDIF
          esqrt = esqrt * Gp(iw,iwp) * a2f_iso(iwp) 
          Znorm(iw) = Znorm(iw) - ws(i) * esqrt 
          Delta(iw) = Delta(iw) - Deltap(i) * esqrt 
        ENDIF
        ! 
        i = abs(iw - iwp) + 1
        root = sqrt( Znormp(i)**2.d0 * ( ws(i)**2.d0 - Deltap(i)**2.d0 ) )
        IF ( aimag(root) .lt. zero ) THEN 
           esqrt = Znormp(i) / conjg(root)
        ELSE  
           esqrt = Znormp(i) / root
        ENDIF
        esqrt = esqrt * Gm(iw,iwp) * a2f_iso(iwp) 
        IF ( iw .lt. iwp ) THEN 
           Znorm(iw) = Znorm(iw) - ws(i) * esqrt 
        ELSE
           Znorm(iw) = Znorm(iw) + ws(i) * esqrt 
        ENDIF
        Delta(iw) = Delta(iw) + Deltap(i) * esqrt
      ENDDO ! iwp
      Znorm(iw) = 1.d0 + pi * ( - estemp(itemp) * Zsumi(iw) + ci * Znorm(iw) * dwsph ) / ws(iw)
      Delta(iw) = pi * ( estemp(itemp) * Dsumi(iw) + ci * Delta(iw) * dwsph ) / Znorm(iw)
      reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
      absdelta = absdelta + abs( Delta(iw) ) 
    ENDDO ! iw 
    errdelta = reldelta / absdelta
    Deltaold(:) = Delta(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
                 '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
                 '   Re[Delta(1)] = ', real(Delta(1))
    !
    IF ( errdelta .lt. conv_thr_racon ) conv = .true.
    IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
       cname = 'acon'
       CALL eliashberg_write_raxis( itemp, cname )
    ENDIF
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
       IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
       WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
       CALL errore('analytic_cont_iso_iaxis_to_raxis','increase nsiter or reduce conv_thr_racon',-1)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE analytic_cont_iso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_iso_iaxis_to_raxis( itemp, N, conv )
    !-----------------------------------------------------------------------
    !
    ! This routine uses pade approximants to continue the isotropic Eliashberg equations 
    ! from the imaginary-axis to the real-axis
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Delta, Znorm, Deltai, Znormi
    USE constants_epw, ONLY : cone, ci, zero, czero
    USE superconductivity, ONLY : pade_coeff, pade_eval
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT (in) :: N
    !! Maximum number of frequency 
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged
    ! 
    ! Local variable
    INTEGER :: iw
    !! Counter on frequency imag- and real-axis
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta 
    ! arrays used in pade_coeff and pade_eval
    COMPLEX(DP) :: a(N)
    !! a - pade coeff for Deltai
    COMPLEX(DP) :: b(N)
    !! b - pade coeff for Znormi
    COMPLEX(DP) :: z(N)
    !! z - frequency imag-axis
    COMPLEX(DP) :: u(N)
    !! u - Deltai
    COMPLEX(DP) :: v(N)
    !! v - Znormi 
    COMPLEX(DP) :: omega
    !! frequency real-axis
    COMPLEX(DP) :: padapp
    !! Znorm or Delta on real-axis after pade_eval
    COMPLEX(DP) :: Deltaold(nsw)
    !! gap
    CHARACTER (len=256) :: cname
    !! character in file name 
    !
    Deltaold(:) = gap(itemp)
    absdelta = zero
    reldelta = zero
    !
    IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
    IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
    Znorm(:) = czero
    Delta(:) = czero
    !
    DO iw = 1, N
       z(iw) = ci * wsi(iw)
       u(iw) = cone * Deltai(iw)
       v(iw) = cone * Znormi(iw)
    ENDDO
    !
    CALL pade_coeff( N, z, u, a )
    CALL pade_coeff( N, z, v, b )
    !
    DO iw = 1, nsw
       omega = cone * ws(iw)
       CALL pade_eval( N, z, a, omega, padapp )
       Delta(iw) = padapp
       CALL pade_eval( N, z, b, omega, padapp )
       Znorm(iw) = padapp
       reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
       absdelta = absdelta + abs( Delta(iw) )
    ENDDO
    errdelta = reldelta / absdelta
    !
    IF ( errdelta .gt. zero ) THEN 
       conv = .true.
       WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'pade = ', N, & 
              '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
              '   Re[Delta(1)] = ', real(Delta(1))
       cname = 'pade'
       CALL eliashberg_write_raxis( itemp, cname )
    ENDIF
    !
  !  IF ( .not. conv ) THEN
  !     WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached pade = ', N
  !     CALL errore('pade_cont_iso_iaxis_to_raxis','decrease number of Pade approximants',-1)
  !  ENDIF
    !
    RETURN
    !
    END SUBROUTINE pade_cont_iso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !  
    ! computes kernels K_{+}(n,n',T) and K_{-}(n,n'T)
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    !
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : pi, zero
    USE eliashbergcom, ONLY : nsiw, estemp, Keri
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw, n
    !! Counter on frequency imag-axis
    REAL(DP) :: omega
    !! frequency imag-axis
    REAL(DP) :: lambda_eph
    !! electron-phonon coupling
    !
    IF ( .not. ALLOCATED(Keri) ) ALLOCATE( Keri(2*nsiw(itemp)) )
    Keri(:) = zero
    !
    DO iw = 1, 2*nsiw(itemp)
       n = iw - 1
       omega = dble(2*n) * pi * estemp(itemp)
       CALL lambdar_iso( omega, lambda_eph )
       Keri(iw) = lambda_eph
    ENDDO 
    !
    RETURN
    !
    END SUBROUTINE kernel_iso_iaxis                                                       
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_iso( omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(n-n')   
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nqstep
    USE constants_epw, ONLY : zero
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    ! 
    IMPLICIT NONE
    !
    REAL(DP), INTENT (in) :: omega
    !! frequency on imag-axis
    REAL(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda(n-n')
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency 
    !
    lambda_eph = zero
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       lambda_eph = lambda_eph + wsph(iwph) * a2f_iso(iwph) & 
                  / ( wsph(iwph)**2.d0 + omega**2.d0 )
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph 
    !
    RETURN
    !
    END SUBROUTINE lambdar_iso
  
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_iaxis_analytic_cont( itemp )
    !-----------------------------------------------------------------------
    !  
    ! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : muc
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, Deltai, Dsumi, Zsumi
    USE constants_epw, ONLY : zero
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: iwp
    !! Counter on frequency imag-axis
    !
    REAL(DP) :: esqrt, kernelp, kernelm
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    !! Temporary working variables
    !
    COMPLEX(DP) :: lambda_eph
    !! electron-phonon coupling lambda(w-iw_n)
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(Dsumi) )  ALLOCATE( Dsumi(nsw) )
    IF ( .not. ALLOCATED(Zsumi) )  ALLOCATE( Zsumi(nsw) )
    Dsumi(:) = zero
    Zsumi(:) = zero
    !
    DO iw = 1, nsw ! loop over omega
       DO iwp = 1, nsiw(itemp) ! loop over iw_n
          CALL lambdai_iso( ws(iw), wsi(iwp), lambda_eph )
          kernelp = 2.d0 * real(lambda_eph)
          kernelm = 2.d0 * aimag(lambda_eph) 
          IF ( iw .eq. 1 ) THEN
             esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + Deltai(iwp)**2.d0 )
             wesqrt(iwp) =  wsi(iwp) * esqrt
             desqrt(iwp) =  Deltai(iwp) * esqrt
          ENDIF
          Zsumi(iw) = Zsumi(iw) + kernelm * wesqrt(iwp)
          Dsumi(iw) = Dsumi(iw) + ( kernelp - 2.d0 * muc ) * desqrt(iwp)
       ENDDO
    ENDDO
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE (wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE (desqrt)
    !   
    RETURN
    !
    END SUBROUTINE kernel_iso_iaxis_analytic_cont      
    !                                                
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_iso( omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(w-iw_n)   
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    USE kinds, ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    USE constants_epw, ONLY : ci, czero
    ! 
    IMPLICIT NONE
    !
    REAL(DP), INTENT (in) :: omega
    !! frequency w at point iw on real-axis
    REAL(DP), INTENT (in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda(w-iw_n)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !    
    lambda_eph = czero
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       lambda_eph = lambda_eph & 
                  + wsph(iwph) * a2f_iso(iwph) / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
    ENDDO ! iwph
    lambda_eph = lambda_eph * 2.d0 * dwsph 
    !
    RETURN
    !
    END SUBROUTINE lambdai_iso
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gamma_acont( omega, omegap, temp, rgammap, rgammam )
    !-----------------------------------------------------------------------
    !!
    !! computes gammam(w,wp)  (notes RM)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : eps6, zero, one
    !
    IMPLICIT NONE
    !
    REAL(kind=DP), INTENT (in) :: omega
    !! frequency w at point iw on the real-axis
    REAL(kind=DP), INTENT (in) :: omegap
    !! frequency w' at point iwp on the real-axis
    REAL(kind=DP), INTENT (in) :: temp
    !! temperature in eV
    REAL(kind=DP), INTENT (out) :: rgammap
    !! -bose_einstein( w' ) - fermi_dirac(  w + w' )
    REAL(kind=DP), INTENT (out) :: rgammam
    !! bose_einstein( w' ) + fermi_dirac( -w + w' )
    !
    rgammap = zero
    rgammam = zero
    IF ( ABS(temp) < eps6 ) THEN
       rgammap = zero
       rgammam = one
    ELSEIF ( omegap .gt. zero ) THEN
       rgammap = 0.5d0 * (   tanh( 0.5d0 * ( omega + omegap ) / temp ) &
                           - 1.d0 / tanh( 0.5d0 * omegap / temp ) )
       rgammam = 0.5d0 * (   tanh( 0.5d0 * ( omega - omegap ) / temp ) &
                           + 1.d0 / tanh( 0.5d0 * omegap / temp ) )
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE gamma_acont
    !
    !-----------------------------------------------------------------------
    !          
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                              
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic 
    !! Eliashberg equations directly on the real-axis.  
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim
    USE eliashbergcom, ONLY : nsw, Delta, Deltap, gap, estemp
    USE constants_epw, ONLY : kelvin2eV, ci
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE superconductivity, ONLY : deallocate_eliashberg, gen_freqgrid_raxis
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    REAL(DP) :: rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(DP) :: tcpu
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv
    !
    CALL start_clock( 'iso_raxis' ) 
    !
    WRITE(stdout,'(5x,a)') 'Solve isotropic Eliashberg equations on real-axis'
    !
    CALL gen_freqgrid_raxis
    !
    DO itemp = 1, nstemp ! loop over temperature
       !
       WRITE(stdout,'(a)') '    '
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
       WRITE(stdout,'(a)') '    '
       iter = 1
       conv = .false.
       DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL integrate_eliashberg_iso_raxis( itemp, iter, conv )
          rdeltain(:) = real(Deltap(:))
          cdeltain(:) = aimag(Deltap(:))
          rdeltaout(:) = real(Delta(:))
          cdeltaout(:) = aimag(Delta(:))
          CALL mix_broyden ( nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
          CALL mix_broyden2( nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
          Deltap(:) = rdeltain(:) + ci * cdeltain(:)
          iter = iter + 1
       ENDDO ! iter
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a,f10.6,a,a,f10.6,a)') &
                    'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K ', &
                    '  gap_edge(', itemp, ') = ', gap(itemp), ' eV ', &
                    '  Re[Delta(1)] = ', real(Delta(1)), ' eV '
       WRITE(stdout,'(a)') '    '
       tcpu = get_clock( 'iso_raxis' )
       WRITE( stdout,'(5x,a,i3,a,f8.1,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
       !
       IF ( conv ) THEN
          WRITE(stdout,'(a)') '    '
          CALL print_clock( 'iso_raxis' )
          WRITE(stdout,'(a)') '    '
       ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'iso_raxis' )
          CALL print_clock( 'iso_raxis' )
          CALL errore('integrate_eliashberg_iso_raxis','converged was not reached',1)
          RETURN
       ENDIF
       !
    ENDDO ! itemp
    !    
    CALL stop_clock( 'iso_raxis' )
    ! 
    RETURN
    !
    END SUBROUTINE eliashberg_iso_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE integrate_eliashberg_iso_raxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations directly on the real-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilker
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : nswfc, nqstep, nsiter, muc, conv_thr_raxis, &
                              kerwrite, kerread, nstemp
    USE eliashbergcom, ONLY : nsw, estemp, ws, dws, gap0, gap, fdwp, Kp, Km, & 
                              Delta, Deltap, Znorm
    USE constants_epw, ONLY : kelvin2eV, ci, eps6, zero, czero
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT (in) :: iter
    !! Counter on iteration steps
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged  
    ! 
    ! Local variables
    INTEGER :: iw, iwp
    !! Counter on frequency real-axis
    REAL(DP) :: dstep
    !!
    REAL(DP) :: temp
    !! Temperature in K
    REAL(DP) :: a, b, c, d
    !! Temporary variables for reading kernelp and kernelm from file
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    COMPLEX(DP) :: esqrt
    !! Temporary working variables
    COMPLEX(DP) :: kernelp, kernelm
    !! Temporary arrays for kernels K_{+}(w,w',T) and K_{-}(w,w'T)
    COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !! gap
    REAL(DP), EXTERNAL :: wgauss
    CHARACTER(len=256) :: name1, cname
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsw) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsw) )
    !
    IF ( iter .eq. 1 ) THEN 
       IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
       IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
       IF ( .not. ALLOCATED(Deltap) ) ALLOCATE( Deltap(nsw) )
       IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
       gap(itemp) = zero
       Deltap(:)  = czero
       Deltap(:)  = gap0 
       IF ( .not. ALLOCATED(fdwp) ) ALLOCATE( fdwp(nsw) )
       IF ( .not. ALLOCATED(Kp) )   ALLOCATE( Kp(nsw,nsw) )
       IF ( .not. ALLOCATED(Km) )   ALLOCATE( Km(nsw,nsw) )
    ENDIF
    Delta(:) = czero
    Znorm(:) = czero
    !
    temp = estemp(itemp) / kelvin2eV
    IF ( temp .lt. 10.d0 ) THEN  
       WRITE(name1,'(a,a7,f4.2)') TRIM(prefix),'.ker_00', temp
    ELSEIF ( temp .ge. 10.d0 ) THEN 
       WRITE(name1,'(a,a6,f5.2)') TRIM(prefix),'.ker_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN 
       WRITE(name1,'(a,a5,f6.2)') TRIM(prefix),'.ker_', temp
    ENDIF
    OPEN(iufilker, file=name1, form='unformatted')
    !
    IF ( iter .eq. 1 ) THEN
       IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
       Deltaold(:) = gap0
    ENDIF          
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nsw ! loop over omega_prime
        IF ( iter .eq. 1 ) THEN
          IF ( iw .eq. 1 ) THEN
            IF ( ABS(estemp(itemp)) <  eps6 ) THEN
               fdwp(iwp) = zero
            ELSE
               fdwp(iwp) = wgauss( -ws(iwp) / estemp(itemp), -99 )
            ENDIF
          ENDIF
          !
          ! read the kernels from file if they were calculated before otherwise calculate them
          IF ( kerread ) THEN 
            READ(iufilker) a, b, c, d
            Kp(iw,iwp) = a + ci*b
            Km(iw,iwp) = c + ci*d
          ENDIF
          IF ( kerwrite ) THEN 
            CALL kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
            Kp(iw,iwp) = kernelp
            Km(iw,iwp) = kernelm
            WRITE(iufilker) real(Kp(iw,iwp)), aimag(Kp(iw,iwp)), &
                            real(Km(iw,iwp)), aimag(Km(iw,iwp))
          ENDIF
        ENDIF
        !
        ! this step is performed at each iter step only for iw=1 since it is independent of w(iw)
        IF ( iw .eq. 1 ) THEN
           esqrt = 1.d0 / sqrt( ws(iwp)**2.d0 - Deltap(iwp)**2.d0 )
           wesqrt(iwp) =  real( ws(iwp) * esqrt )
           desqrt(iwp) =  real( Deltap(iwp) * esqrt )
        ENDIF
        !
        ! end points contribute only half ( trapezoidal integration rule )
        IF ( (iwp .eq. 1) .OR. (iwp .eq. nsw) ) THEN
           dstep = 0.5d0 * dws(iwp) 
        ! boundary points contribute half from left and half from right side
        ELSEIF ( iwp .eq. nswfc ) THEN
           dstep = 0.5d0 * ( dws(iwp) + dws(iwp+1) )
        ELSE
           dstep = dws(iwp)
        ENDIF 
        Znorm(iw) = Znorm(iw) + dstep * wesqrt(iwp) * Km(iw,iwp)
        Delta(iw) = Delta(iw) + dstep * desqrt(iwp) &
                  * ( Kp(iw,iwp) - muc*( 1.d0 - 2.d0*fdwp(iwp) ) )
      ENDDO ! iwp
      Znorm(iw) = 1.d0 - Znorm(iw) / ws(iw)
      Delta(iw) = Delta(iw) / Znorm(iw)
      reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) * dws(iw)
      absdelta = absdelta + abs( Delta(iw) ) * dws(iw)
    ENDDO ! iw 
    CLOSE(iufilker)
    errdelta = reldelta / absdelta
    Deltaold(:) = Delta(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, '   error = ', errdelta, & 
                                  '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1)) 
    !
    IF ( errdelta .lt. conv_thr_raxis) conv = .true.
    IF ( errdelta .lt. conv_thr_raxis .OR. iter .eq. nsiter ) THEN
      cname = 'real'
      CALL eliashberg_write_raxis( itemp, cname )
      gap0 = gap(itemp)
    ENDIF
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
       IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
       CALL errore('integrate_eliashberg_iso_raxis','increase nsiter or reduce conv_thr_raxis',-1)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE integrate_eliashberg_iso_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
    !-----------------------------------------------------------------------
    !
    ! computes kernels K_{+}(w,w',T) and K_{-}(w,w'T)  
    ! reference M. J. Holcomb, PRB 54, 6648 (1996)   
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi, ci, eps6, zero, czero
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, bewph, wsph, dwsph, ws, fdwp, estemp
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: iw
    !! index frequency w : ws(iw)
    INTEGER, INTENT (in) :: iwp
    !! index frequency w' : ws(iwp)
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    COMPLEX(DP), INTENT (out) :: kernelp
    !! phonon kernel K_{+}(w,w',T) 
    COMPLEX(DP), INTENT (out) :: kernelm
    !! phonon kernel K_{-}(w,w',T)
    !
    ! Local variables 
    INTEGER :: iwph
    !! Counter on frequency
    INTEGER :: ngaussw0
    REAL(DP) :: degaussw0, f1, f2, f3, f4
    COMPLEX(DP) :: e1, e2, e3, e4
    REAL(DP), EXTERNAL :: wgauss, w0gauss
    !
    degaussw0 = 1.d0 * dwsph
    ngaussw0 = 0
    !
    f1 = zero
    f2 = zero
    f3 = zero
    f4 = zero
    kernelp = czero
    kernelm = czero
    e1 = czero
    e2 = czero
    e3 = czero
    e4 = czero
    !
    IF ( .not. ALLOCATED(bewph) ) ALLOCATE( bewph(nqstep) )
    ! Bose-Einstein distribution
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       IF ( iw .eq. 1 .AND. iwp .eq. 1 ) THEN
          IF ( ABS(estemp(itemp)) <  eps6 ) THEN
             bewph(iwph)  = zero
          ELSE
             bewph(iwph) = wgauss( -wsph(iwph) / estemp(itemp), -99 )
             bewph(iwph) = bewph(iwph) / (1.d0 - 2.d0 * bewph(iwph))
          ENDIF
       ENDIF
       !
       ! a small complex number is added to denominator to move the pole away from the real-axis
       !
       ! in order to reduce the numerical noise at very small frequencies coming from the complex number 
       ! added in the denominator, the contribution of the imaginary part is reestimated using 
       ! delta function (RM notes) 
       !
       ! subtract the imaginary part coming from e1 to e4 and add instead the imaginary part 
       ! coming from f1 to f4
       !
       e1 = 1.d0 / ( wsph(iwph) + ws(iwp) + ws(iw) + ci*degaussw0 ) 
       e2 = 1.d0 / ( wsph(iwph) + ws(iwp) - ws(iw) - ci*degaussw0 ) 
       e3 = 1.d0 / ( wsph(iwph) - ws(iwp) + ws(iw) + ci*degaussw0 ) 
       e4 = 1.d0 / ( wsph(iwph) - ws(iwp) - ws(iw) - ci*degaussw0 ) 
       !
       ! estimate of the imaginary part using delta function
       f1 = w0gauss( ( wsph(iwph) + ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f2 = w0gauss( ( wsph(iwph) + ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f3 = w0gauss( ( wsph(iwph) - ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f4 = w0gauss( ( wsph(iwph) - ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       !
       kernelp = kernelp + a2f_iso(iwph) &
               * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 + e2 - ci*aimag(e2) + ci*pi*f2 ) & 
                  - (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 + e4 - ci*aimag(e4) + ci*pi*f4 ) )
       kernelm = kernelm + a2f_iso(iwph) &
               * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 - e2 + ci*aimag(e2) - ci*pi*f2 ) &
                  + (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 - e4 + ci*aimag(e4) - ci*pi*f4 ) )
    ENDDO ! iwph
    kernelp = kernelp * dwsph 
    kernelm = kernelm * dwsph 
    !
    RETURN
    !
    END SUBROUTINE kernel_raxis
    !
    !-----------------------------------------------------------------------
    ! 
  END MODULE superconductivity_iso
