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
  !
  ! This routine is the driver of the self-consistent cycle for the isotropic 
  ! Eliashberg equations on the imaginary-axis.  
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE control_flags, ONLY : iverbosity
  USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                            limag, lpade, lacon, imag_read, pade_read 
  USE eliashbergcom, ONLY : nsw, nsiw, Deltai, Deltaip, Delta, Deltap, estemp, gap
  USE constants_epw, ONLY : kelvin2eV, ci
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: itemp, iter, N
  REAL(DP) :: dFE, tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
  REAL(DP), EXTERNAL :: get_clock
  LOGICAL :: conv
  !
  CALL start_clock( 'iso_iaxis' )
!#ifdef __PARA 
!  IF (mpime .eq. ionode_id) THEN
!#endif
  DO itemp = 1, nstemp ! loop over temperature
     !
     WRITE(stdout,'(a)') '    '
     WRITE(stdout,'(5x,a,i3,a,f12.5,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
     WRITE(stdout,'(a)') '    '
     IF ( limag .AND. .not. imag_read ) THEN
        WRITE(stdout,'(5x,a)') 'Solve isotropic Eliashberg equations on imaginary-axis' 
     ELSEIF ( limag .AND. imag_read ) THEN
        WRITE(stdout,'(5x,a)') 'Read from file Delta and Znorm on imaginary-axis '
     ENDIF
     WRITE(stdout,'(a)') '    '
     WRITE(stdout,'(5x,a,i6,a,i6)') 'Total number of frequency points nsiw ( ', itemp, ' ) = ', nsiw(itemp)
     WRITE(stdout,'(a)') '    '
     CALL start_clock( 'iaxis_imag' )
     CALL gen_freqgrid_iaxis( itemp )
     !
     IF ( ( limag .AND. .not. imag_read ) .OR. ( limag .AND. imag_read .AND. itemp .ne. 1 ) ) THEN
!     IF ( limag .AND. .not. imag_read ) THEN
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
            CALL free_energy( itemp, dFE )
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
     ELSEIF ( limag .AND. imag_read .AND. itemp .eq. 1 ) THEN
!     ELSEIF ( limag .AND. imag_read ) THEN
        CALL read_eliashberg_iso_iaxis( itemp )
     ENDIF
     !
     If ( ( lpade .AND. .not. pade_read ) .OR. ( lpade .AND. pade_read .AND. itemp .ne. 1 ) ) THEN
!     IF ( lpade .AND. .not. pade_read ) THEN 
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
     ELSEIF ( lpade .AND. pade_read .AND. itemp .eq. 1 ) THEN
!     ELSEIF ( lpade .AND. pade_read ) THEN
        WRITE(stdout,'(5x,a)') 'Read from file Delta and Znorm on real-axis '
        CALL read_eliashberg_iso_raxis_pade( itemp )
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
!#ifdef __PARA 
!  ENDIF
!  CALL mp_barrier()
!#endif 
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
  !
  ! This routine solves the isotropic Eliashberg equations on the imaginary-axis
  !
  ! input
  !
  ! itemp  - temperature point
  ! iter   - iteration number
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iufilgap
  USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis
  USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, wsi, NZnormi, Znormi, Deltai, Deltaip, Keri
  USE constants_epw, ONLY : pi, kelvin2eV
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, iwp, itemp, iter
  REAL(DP) :: esqrt, kernelp, kernelm, lambdap, lambdam, absdelta, reldelta, errdelta, temp
  REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
  REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  LOGICAL :: conv
  CHARACTER(len=256) :: name1
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
     IF ( itemp .eq. 1 ) THEN 
        Deltaip(:) = gap0
     ELSE
        Deltaip(:) = gap(itemp-1)
        !Deltaip(:) = 2.d0*gap(1)*sqrt(1.d0 - 3.52d0*estemp(itemp)/2.d0/gap(1) )
     ENDIF
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
  WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9)') 'iter = ', iter, '   error = ', errdelta, &
                                              '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
  !
  IF ( errdelta .lt. conv_thr_iaxis ) conv = .true.
  IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
     temp = estemp(itemp) / kelvin2eV 
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a11,f10.5)') TRIM(prefix),'.imag_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a10,f10.5)') TRIM(prefix),'.imag_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(4a24)') 'w', 'Znorm(w)', 'Delta(w)', 'NZnorm(w)'
     DO iw = 1, nsiw(itemp) ! loop over omega
        WRITE(iufilgap,'(4ES20.10)') wsi(iw), Znormi(iw), Deltai(iw), NZnormi(iw)
     ENDDO
     gap(itemp) = Deltai(1)
     CLOSE(iufilgap)
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
  SUBROUTINE analytic_cont_iso_iaxis_to_raxis( itemp, iter, conv ) 
  !-----------------------------------------------------------------------
  !
  ! This routine does the analyic continuation of the isotropic Eliashberg equations 
  ! from the imaginary-axis to the real axis
  ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
  !
  ! input
  !
  ! itemp  - temperature point
  ! iter   - iteration number
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iufilgap
  USE epwcom,        ONLY : nqstep, nsiter, conv_thr_racon, lpade
  USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, gap, a2f_iso, Dsumi, Zsumi, & 
                            Delta, Deltap, Znorm, Znormp, Gp, Gm
  USE constants_epw, ONLY : pi, kelvin2eV, ci
  ! 
  IMPLICIT NONE
  !
  INTEGER :: i, iw, iwp, itemp, iter 
  REAL(DP) :: rgammap, rgammam, absdelta, reldelta, errdelta, temp
  COMPLEX(DP) :: esqrt, root
  COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  LOGICAL :: conv, lgap
  CHARACTER(len=256) :: name1
  !
  IF ( iter .eq. 1 ) THEN
     IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
     IF ( .not. ALLOCATED(Deltap) )   ALLOCATE( Deltap(nsw) )
     IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
     IF ( .not. ALLOCATED(Znormp) )   ALLOCATE( Znormp(nsw) )
     IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
     Deltap(:) = (0.d0, 0.d0)
     Deltaold(:) = (0.d0, 0.d0)
     IF ( lpade ) THEN
        Deltap(:) = Delta(:)
        Deltaold(:) = Delta(:)
     ELSE 
        Deltap(:) = gap(itemp)
        Deltaold(:) = gap(itemp)
     ENDIF
     Znormp(:) = (1.d0, 0.d0)
     IF ( .not. ALLOCATED(Gp) ) ALLOCATE( Gp(nsw,nqstep) )
     IF ( .not. ALLOCATED(Gm) ) ALLOCATE( Gm(nsw,nqstep) )
     IF ( .not. ALLOCATED(Dsumi) ) ALLOCATE( Dsumi(nsw) )
     IF ( .not. ALLOCATED(Zsumi) ) ALLOCATE( Zsumi(nsw) )
     CALL kernel_iso_iaxis_analytic_cont( itemp, Zsumi, Dsumi )
  ENDIF
  Znorm(:) = (0.d0, 0.d0)
  Delta(:) = (0.d0, 0.d0)
  !
  absdelta = 0.d0   
  reldelta = 0.d0
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
           IF ( aimag(root) .lt. 0.d0 ) THEN 
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
        IF ( aimag(root) .lt. 0.d0 ) THEN 
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
  WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9)') 'iter = ', iter, '   error = ', errdelta, &
                              '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1))
  !
  IF ( errdelta .lt. conv_thr_racon ) conv = .true.
  IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
     temp = estemp(itemp) / kelvin2eV
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a11,f10.5)') TRIM(prefix),'.acon_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a10,f10.5)') TRIM(prefix),'.acon_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     lgap = .true.
     DO iw = 1, nsw
        IF ( lgap .AND. iw .lt. nqstep .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. & 
             ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN 
           gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                      / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )  
           !WRITE(stdout,'(a)') '   '
           !WRITE(stdout,'(5x,a)') 'gap_edge'
           !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw, ws(iw), real(Delta(iw)), aimag(Delta(iw)), & 
           !                                 gap(itemp)
           !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw+1, ws(iw+1), real(Delta(iw+1)), aimag(Delta(iw+1)), &
           !                                 gap(itemp)
           !WRITE(stdout,'(a)') '   '
           !
           lgap = .false. 
           !
        ENDIF
        WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                   real(Delta(iw)), aimag(Delta(iw))
     ENDDO
     CLOSE(iufilgap)
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
  ! input
  !
  ! itemp  - temperature point
  ! iter   - iteration number
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iufilgap
  USE epwcom,        ONLY : nqstep
  USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Delta, Znorm, Deltai, Znormi, estemp
  USE constants_epw, ONLY : kelvin2eV, cone, ci
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, N
  REAL(DP) :: absdelta, reldelta, errdelta, temp
  COMPLEX(DP) :: a(N), b(N), z(N), u(N), v(N)
  COMPLEX(DP) :: omega, padapp, Deltaold(nsw)
  LOGICAL :: lgap, conv
  CHARACTER(len=256) :: name1
  !
  Deltaold(:) = gap(itemp)
  absdelta = 0.d0
  reldelta = 0.d0
  !
  IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
  Znorm(:) = (0.d0, 0.d0)
  Delta(:) = (0.d0, 0.d0)
  !
  DO iw = 1, N
     z(iw) = ci * wsi(iw)
     u(iw) = cone * Deltai(iw)
     v(iw) = cone * Znormi(iw)
!     u(iw) = cone * Deltai(iw) / wsi(iw)
!     v(iw) = cone * Znormi(iw) * Deltai(iw)
  ENDDO
  !
  CALL pade_coeff( N, z, u, a )
  CALL pade_coeff( N, z, v, b )
  !
  DO iw = 1, nsw
     omega = cone * ws(iw)
     CALL pade_eval( N, z, a, omega, padapp )
     Delta(iw) = padapp
!     Delta(iw) = padapp * omega
     CALL pade_eval( N, z, b, omega, padapp )
     Znorm(iw) = padapp
!     Znorm(iw) = padapp / Delta(iw)
     reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
     absdelta = absdelta + abs( Delta(iw) )
  ENDDO
  errdelta = reldelta / absdelta
  !
  IF ( errdelta .gt. 0.d0 ) THEN 
     conv = .true.
     WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9)') 'pade = ', N, '   error = ', errdelta, &
                  '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1))
     temp = estemp(itemp) / kelvin2eV
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a11,f10.5)') TRIM(prefix),'.pade_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a10,f10.5)') TRIM(prefix),'.pade_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     lgap = .true.
     DO iw = 1, nsw
        WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                  real(Delta(iw)), aimag(Delta(iw))
        IF ( lgap .AND. iw .lt. nqstep .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
             ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN 
           gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                   / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )  
           !WRITE(stdout,'(a)') '   '
           !WRITE(stdout,'(5x,a)') 'gap_edge'
           !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw, ws(iw), real(Delta(iw)), aimag(Delta(iw)), &
           !                                 gap(itemp)
           !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw+1, ws(iw+1), real(Delta(iw+1)), aimag(Delta(iw+1)), &
           !                                 gap(itemp)
           !WRITE(stdout,'(a)') '   '        
           !
           lgap = .false. 
           ! 
        ENDIF
     ENDDO
     CLOSE(iufilgap)
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
  SUBROUTINE read_eliashberg_iso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the anisotropic Delta and Znorm on the imaginary-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iufilgap
  USE epwcom,        ONLY : nstemp
  USE eliashbergcom, ONLY : nsiw, estemp, gap, wsi, Znormi, NZnormi, Deltai
  USE constants_epw, ONLY : kelvin2eV
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ios
  REAL(DP) :: temp, tmp
  CHARACTER (len=256) :: name1, word
  !
  IF ( .not. ALLOCATED(gap) )     ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Deltai) )  ALLOCATE( Deltai(nsiw(itemp)) )
  IF ( .not. ALLOCATED(Znormi) )  ALLOCATE( Znormi(nsiw(itemp)) )
  IF ( .not. ALLOCATED(NZnormi) ) ALLOCATE( NZnormi(nsiw(itemp)) )
  Deltai(:) = 0.d0
  Znormi(:) = 0.d0
  NZnormi(:) = 0.d0
  !
  temp = estemp(itemp) / kelvin2eV
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f10.5)') TRIM(prefix),'.imag_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f10.5)') TRIM(prefix),'.imag_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_iso_iaxis','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO iw = 1, nsiw(itemp) ! loop over omega
     READ(iufilgap,'(4ES20.10)') tmp, Znormi(iw), Deltai(iw), NZnormi(iw)
     IF ( abs(wsi(iw)-tmp) .gt. 10.d-8 ) & 
        CALL errore('read_eliashberg_iso_iaxis','temperature not the same with the input',1)        
  ENDDO
  gap(itemp) = Deltai(1)
  CLOSE(iufilgap)
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_iso_iaxis
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_eliashberg_iso_raxis_analytic_cont( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the isotropic Delta and Znorm on the real-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : prefix 
  USE epwcom,        ONLY : nstemp
  USE io_epw,        ONLY : iufilgap
  USE eliashbergcom, ONLY : nsw, estemp, gap, ws, Znorm, Delta
  USE constants_epw, ONLY : kelvin2eV, ci
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ios
  REAL(DP) :: a, b, c, d, temp, omega
  CHARACTER (len=256) :: name1, word
  !
  IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  !
  temp = estemp(itemp) / kelvin2eV
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.acon_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.acon_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_iso_raxis_analytic_cont','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO iw = 1, nsw ! loop over omega
     READ(iufilgap,'(5ES20.10)') omega, a, b, c, d
           Znorm(iw) = a + ci*b
           Delta(iw) = c + ci*d
     IF ( abs(ws(iw)-omega) .gt. 10.d-8 ) &
        CALL errore('read_eliashberg_iso_raxis_analytic_cont','temperature not the same with the input',1)
  ENDDO
  gap(itemp) = REAL(Delta(1))
  CLOSE(iufilgap)
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_iso_raxis_analytic_cont
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_eliashberg_iso_raxis_pade( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the isotropic Delta and Znorm on the real-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : prefix
  USE io_epw,        ONLY : iufilgap
  USE epwcom,        ONLY : nstemp
  USE eliashbergcom, ONLY : nsw, estemp, gap, ws, Znorm, Delta
  USE constants_epw, ONLY : kelvin2eV, ci
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ios
  REAL(DP) :: a, b, c, d, temp, omega
  CHARACTER (len=256) :: name1, word
  !
  IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  !   
  temp = estemp(itemp) / kelvin2eV
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f10.5)') TRIM(prefix),'.pade_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f10.5)') TRIM(prefix),'.pade_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_iso_raxis_pade','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO iw = 1, nsw ! loop over omega
     READ(iufilgap,'(5ES20.10)') omega, a, b, c, d
           Znorm(iw) = a + ci*b
           Delta(iw) = c + ci*d
     IF ( abs(ws(iw)-omega) .gt. 10.d-8 ) &
        CALL errore('read_eliashberg_iso_raxis_pade','temperature not the same with the input',1)
  ENDDO
  gap(itemp) = REAL(Delta(1))
  CLOSE(iufilgap)
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_iso_raxis_pade
  !
  !-----------------------------------------------------------------------


