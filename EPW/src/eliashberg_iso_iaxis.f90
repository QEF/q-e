  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
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
  ! 
  IMPLICIT NONE
  !
  ! Local variables
  INTEGER :: itemp
  !! Counter on temperature index
  INTEGER :: iter
  !! Counter on iteration steps
  INTEGER :: N
  !! Maximum frequency 
  REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
  REAL(DP), EXTERNAL :: get_clock
  LOGICAL :: conv
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
  USE constants_epw, ONLY : pi, kelvin2eV
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: itemp
  !! Counter on temperature
  INTEGER, INTENT(in) :: iter
  !! Counter on iteration steps
  LOGICAL, INTENT(inout) :: conv
  !! True if the calculation is converged
  ! 
  ! Local variables
  INTEGER :: iw, iwp 
  REAL(DP) :: esqrt, kernelp, kernelm, lambdap, lambdam, absdelta, reldelta, errdelta
  REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
  REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
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
     gap(itemp) = 0.d0
     Deltaip(:) = 0.d0
     Deltaip(:) = gap0
     !
     CALL kernel_iso_iaxis( itemp )
  ENDIF
  Znormi(:) = 0.d0
  NZnormi(:) = 0.d0
  Deltai(:) = 0.d0
  !
  IF ( iter .eq. 1 ) THEN
     IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsiw(itemp)) )
     Deltaold(:) = gap0
  ENDIF
  absdelta = 0.d0   
  reldelta = 0.d0 
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
     CALL errore('sum_eliashberg_iso_iaxis','increase nsiter or reduce conv_thr_iaxis',-1)
     CALL deallocate_eliashberg
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE sum_eliashberg_iso_iaxis
  !
  !-----------------------------------------------------------------------
