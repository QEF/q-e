  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino   
  ! Copyright (C) 2007-2009 Roxana Margine, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_iso_raxis
  !-----------------------------------------------------------------------
  !!
  !! This routine is the driver of the self-consistent cycle for the isotropic 
  !! Eliashberg equations on the real-axis.  
  !!
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim
  USE eliashbergcom, ONLY : nsw, Delta, Deltap, gap, estemp
  USE constants_epw, ONLY : kelvin2eV, ci
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  ! Local variables
  INTEGER :: itemp, iter
  REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
  REAL(DP), EXTERNAL :: get_clock
  LOGICAL :: conv 
  CHARACTER (len=256) :: filgap
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
  filgap = TRIM(prefix) // '.gap'
  OPEN(iufilgap, file=filgap, status='unknown')
  DO itemp = 1, nstemp ! loop over temperature
     WRITE(iufilgap,'(2f12.6)') estemp(itemp)/kelvin2eV, gap(itemp)
  ENDDO
  CLOSE(iufilgap)
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
  !! This routine solves the isotropic Eliashberg equations on the real-axis
  !!
  !
  USE kinds,         ONLY : DP
  USE io_epw,        ONLY : iufilker, iufilgap
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE epwcom,        ONLY : nswfc, nqstep, nsiter, muc, conv_thr_raxis, &
                            kerwrite, kerread, nstemp
  USE eliashbergcom, ONLY : nsw, estemp, ws, dws, gap0, gap, fdwp, Kp, Km, & 
                            Delta, Deltap, Znorm
  USE constants_epw, ONLY : kelvin2eV, ci
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
  REAL(DP) :: dstep, a, b, c, d, absdelta, reldelta, errdelta, temp
  REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
  REAL(DP) :: eps=1.0d-6
  COMPLEX(DP) :: kernelp, kernelm, esqrt
  COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  REAL(DP), EXTERNAL :: wgauss
  LOGICAL :: lgap
  CHARACTER(len=256) :: name1, name2
  !
  IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsw) )
  IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsw) )
  !
  IF ( iter .eq. 1 ) THEN 
     IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
     IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
     IF ( .not. ALLOCATED(Deltap) ) ALLOCATE( Deltap(nsw) )
     IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
     gap(itemp) = 0.d0
     Deltap(:)  = (0.d0, 0.d0)
     Deltap(:)  = gap0 
     IF ( .not. ALLOCATED(fdwp) ) ALLOCATE( fdwp(nsw) )
     IF ( .not. ALLOCATED(Kp) )   ALLOCATE( Kp(nsw,nsw) )
     IF ( .not. ALLOCATED(Km) )   ALLOCATE( Km(nsw,nsw) )
  ENDIF
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  !
  temp = estemp(itemp) / kelvin2eV
  IF ( temp .lt. 10.d0 ) THEN  
     WRITE(name2,'(a,a6,f4.2)') TRIM(prefix),'.ker_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN 
     WRITE(name2,'(a,a5,f5.2)') TRIM(prefix),'.ker_', temp
  ENDIF
  !OPEN(iufilker, file=name2, form='formatted')
  OPEN(iufilker, file=name2, form='unformatted')
  !
  IF ( iter .eq. 1 ) THEN
     IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
     Deltaold(:) = gap0
  ENDIF          
  absdelta = 0.d0
  reldelta = 0.d0
  DO iw = 1, nsw ! loop over omega
     DO iwp = 1, nsw ! loop over omega_prime
        IF ( iter .eq. 1 ) THEN
           IF ( iw .eq. 1 ) THEN
              IF ( ABS(estemp(itemp)) <  eps ) THEN
                 fdwp(iwp) = 0.d0
              ELSE
                 fdwp(iwp) = wgauss( -ws(iwp) / estemp(itemp), -99 )
              ENDIF
           ENDIF
           !
           ! read the kernels from file if they were calculated before otherwise calculate them
           IF ( kerread ) THEN 
              !READ(iufilker,'(4ES20.10)') a, b, c, d
              READ(iufilker) a, b, c, d
              Kp(iw,iwp) = a + ci*b
              Km(iw,iwp) = c + ci*d
           ENDIF
           IF ( kerwrite ) THEN 
              CALL kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
              Kp(iw,iwp) = kernelp
              Km(iw,iwp) = kernelm
              !WRITE(iufilker,'(4ES20.10)') real(Kp(iw,iwp)), aimag(Kp(iw,iwp)), &
              !                             real(Km(iw,iwp)), aimag(Km(iw,iwp))
              WRITE(iufilker) real(Kp(iw,iwp)), aimag(Kp(iw,iwp)), &
                              real(Km(iw,iwp)), aimag(Km(iw,iwp))
           ENDIF
        ENDIF
        !
        ! this step is performed at each iter step only for iw=1 since it is independ of w(iw)
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
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a7,f4.2)') TRIM(prefix),'.gapr_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a6,f5.2)') TRIM(prefix),'.gapr_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     !
     WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     lgap = .true.
     ! DO iw = 1, nsw
     DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
        IF ( lgap .AND. iw .lt. (nqstep) .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
             ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN
           gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                      / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )
           lgap = .false.
        ENDIF
        WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                    real(Delta(iw)), aimag(Delta(iw))
     ENDDO
     CLOSE(iufilgap)
     IF ( lgap ) & 
        gap(itemp) =  real(Delta(1))
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
