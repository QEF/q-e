!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE stop_lr( full_run  )
  !----------------------------------------------------------------------------
  !
  ! This subroutine synchronizes processes before stopping.
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : mp_global_end
  USE parallel_include
  USE lr_variables,         ONLY : n_ipol, LR_polarization, beta_store,          &
                                 & gamma_store, zeta_store, norm0, code1,code2,  &
                                 & lr_verbosity, itermax, bgz_suffix,            &
                                   eels, q1, q2, q3              
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : tmp_dir, prefix, iunwfc
  USE environment,          ONLY : environment_end
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin
  USE ions_base,            ONLY : nat, ityp, atm, ntyp => nsp, tau
  USE cell_base,            ONLY : celldm, at, bg, alat, omega
  USE klist,                ONLY : nelec
  USE buffers,              ONLY : close_buffer
#if defined(__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE solvent_tddfpt,       ONLY : solvent_clean_tddfpt
#endif  
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: full_run
  CHARACTER(len=6), EXTERNAL :: int_to_char
  CHARACTER(len=256) :: filename
  INTEGER :: ip,i,j
  REAL(kind=dp) :: degspin
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<stop_lr>")')
  ENDIF
  !
  ! Write beta, gamma, and z coefficents to output directory for
  ! easier post processing. These can also be read from the output log file.
  !
  IF (full_run) THEN
#if defined(__MPI)
  IF (ionode) THEN
#endif
  !
  DO ip = 1,n_ipol
     !
     IF (eels) THEN
       filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
     ELSE
       IF (n_ipol==3) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(ip))
       IF (n_ipol==1) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
     ENDIF
     !
     filename = trim(tmp_dir) // trim(filename)
     !
     OPEN (158, file = filename, form = 'formatted', status = 'replace')
     !
     ! Write the number of iterations
     !
     WRITE(158,*) itermax
     !
     ! Write the norm of the starting Lanczos vectors
     !
     norm0(ip) = beta_store(ip,1)
     WRITE(158,*) norm0(ip)
     !
     IF (nspin==2) THEN
        degspin = 1.0d0
     ELSE
        degspin = 2.0d0
     ENDIF
     IF (noncolin) degspin = 1.0d0
     !
     ! Write the degenaracy wrt spin
     !
     WRITE(158,*) degspin
     !
     ! ------ Needed for EELS ----------
     !
     ! Write the lattice parameter
     !
     WRITE(158,*) alat
     !
     ! Write the unit-cell volume
     !
     WRITE(158,*) omega
     !
     ! Write the number of valence (and semicore electrons) in the unit cell
     !
     WRITE(158,*) nelec
     !
     ! Write the components of the transferred momentum
     !
     WRITE(158,*) q1
     WRITE(158,*) q2
     WRITE(158,*) q3
     !
     !-----------------------------------
     !
     DO i = 1, itermax-1
        !
        WRITE(158,*) beta_store(ip,i+1)
        WRITE(158,*) gamma_store(ip,i+1)
        !
        ! This is absolutely necessary for cross platform compatibilty
        !
        DO j = 1, n_ipol
           WRITE(158,*) zeta_store (ip,j,i)
        ENDDO
        !
     ENDDO
     !
     ! X. Ge: These two faked values will not be
     ! really used in the spectrum calculation.
     !
     WRITE(158,*) beta_store(ip,itermax) 
     WRITE(158,*) gamma_store(ip,itermax) 
     DO j=1,n_ipol                                               
        WRITE(158,*) zeta_store (ip,j,itermax)                    
     ENDDO
     !
     CLOSE(158)
     !
  ENDDO
  !
#if defined(__MPI)
  ENDIF
#endif
  !
  ENDIF
  !
  ! Deallocate lr variables
  !
  CALL lr_dealloc()
  !
#if defined(__ENVIRON)
  !
  ! Deallocate Environ related arrays
  !
  IF (use_environ) CALL solvent_clean_tddfpt()
  !
#endif
  !
  IF (eels) THEN
     CALL environment_end(code2)
  ELSE
     CALL environment_end(code1) 
  ENDIF
  !
  CALL mp_global_end( )
  !
#if defined (__T3E)
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
  !
#endif
  !
  ! EELS: Close the file where it read the wavefunctions at k and k+q.
  !
  IF (eels) CALL close_buffer(iunwfc, 'keep')
  !
  STOP
  !
END SUBROUTINE stop_lr
