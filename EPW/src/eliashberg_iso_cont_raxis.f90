  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
  ! Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
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
  USE constants_epw, ONLY : pi, ci
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: itemp
  !! Counter on iteration
  INTEGER, INTENT(in) :: iter
  !! Counter on the iteration number
  LOGICAL, INTENT(inout) :: conv
  !! True if the calculation is converged
  ! 
  ! Local variables
  INTEGER :: i, iw, iwp 
  REAL(kind=DP) :: rgammap, rgammam, absdelta, reldelta, errdelta
  COMPLEX(kind=DP) :: esqrt, root
  COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  CHARACTER (len=256) :: cname
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
  WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
               '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
               '   Re[Delta(1)] = ', real(Delta(1))
  !
  IF ( errdelta .lt. conv_thr_racon ) conv = .true.
  IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
     cname = 'acon'
     CALL eliashberg_write_cont_raxis( itemp, cname )
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
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Delta, Znorm, Deltai, Znormi
  USE constants_epw, ONLY : cone, ci
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: itemp
  !! Counter on temperature
  INTEGER, INTENT (in) :: N
  !! Maximum number of frequency 
  LOGICAL, INTENT (inout) :: conv
  !! True if the calculation is converged
  ! 
  ! Local variable
  INTEGER :: iw
  REAL(DP) :: absdelta, reldelta, errdelta
  COMPLEX(DP) :: a(N), b(N), z(N), u(N), v(N)
  COMPLEX(DP) :: omega, padapp, Deltaold(nsw)
  CHARACTER (len=256) :: cname
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
     WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'pade = ', N, & 
            '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
            '   Re[Delta(1)] = ', real(Delta(1))
     cname = 'pade'
     CALL eliashberg_write_cont_raxis( itemp, cname )
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
