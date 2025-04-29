!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE ns_adj
!-----------------------------------------------------------------------
   !
   !! This routine adjusts (modifies) the eigenvalues of the atomic
   !! occupation matrix using starting_ns eigenvalues as suggested
   !! by the user from the input. 
   !
   USE kinds,     ONLY : DP
   USE ions_base, ONLY : nat, ntyp => nsp, ityp
   USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, starting_ns,&
                         is_hubbard
   USE scf,       ONLY : rho
   USE lsda_mod,  ONLY : nspin
   USE noncollin_module, ONLY : noncolin, npol
   USE io_global, ONLY : stdout
#if defined (__OSCDFT)
   USE plugin_flags,     ONLY : use_oscdft
   USE oscdft_base,      ONLY : oscdft_ctx
   USE oscdft_functions, ONLY : oscdft_ns_set
#endif
   ! 
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: ldmx=7
   INTEGER :: na, nt, is, m1, m2, ldim, i, j, l 
   REAL(DP) :: lambda(npol*ldmx) 
   COMPLEX(DP) :: vet(npol*ldmx,npol*ldmx), f(npol*ldmx,npol*ldmx), temp
   LOGICAL :: print_new_ns
   !
   print_new_ns = .FALSE.
   !
   IF (ANY(starting_ns /= -1.0_dp)) THEN
     !
     WRITE( stdout, '(/5X,"WARNING!!! Modifying starting ns matrices according to input values")')
     !
     IF (2*Hubbard_lmax+1>ldmx) CALL errore('ns_adj',' ldmx too small',ldmx) 
     !
     DO na = 1, nat
        !
        nt = ityp(na)
        !
        IF (is_hubbard(nt)) THEN
           !
           ldim = 2 * Hubbard_l(nt) + 1 
           !
           IF (noncolin) THEN
              !
              ! Noncollinear case
              !
              DO m1 = 1, ldim
                 DO m2 = 1, ldim
                    f(m1, m2)           = rho%ns_nc(m1, m2, 1, na)
                    f(m1, ldim+m2)      = rho%ns_nc(m1, m2, 2, na)
                    f(ldim+m1, m2)      = rho%ns_nc(m1, m2, 3, na)
                    f(ldim+m1, ldim+m2) = rho%ns_nc(m1, m2, 4, na)
                 ENDDO
              ENDDO
              !
              ! Diagonalize the ocupation matrix
              CALL cdiagh( npol*ldim, f, npol*ldmx, lambda, vet)
              !
              ! Change the eigenvalues as requested from the input
              j = 0 
              DO is = 1, npol
                 DO i = 1, ldim
                    j = j + 1
                    IF (starting_ns(i,is,nt) >= 0.d0) &
                            lambda(j) = starting_ns(i,is,nt)
                 ENDDO
              ENDDO
              !
              ! Reconstruct back the occupation matrix from the
              ! modified eignevalues and the original eigenvectors
              DO m1 = 1, npol*ldim
                 DO m2 = m1, npol*ldim
                    temp = 0.d0
                    DO i = 1, npol*ldim
                       temp = temp + vet(m1,i)*lambda(i)*CONJG(vet(m2,i))     
                    ENDDO
                    f(m1,m2) = temp
                    f(m2,m1) = CONJG(temp) 
                 ENDDO
              ENDDO
              DO m1 = 1, ldim
                 DO m2 = 1, ldim
                    rho%ns_nc(m1, m2, 1, na) = f(m1, m2) 
                    rho%ns_nc(m1, m2, 2, na) = f(m1, ldim+m2) 
                    rho%ns_nc(m1, m2, 3, na) = f(ldim+m1, m2) 
                    rho%ns_nc(m1, m2, 4, na) = f(ldim+m1, ldim+m2) 
                 ENDDO
              ENDDO
              !
           ELSE
              !
              ! Collinear case
              !
              DO is = 1, nspin
                 ! 
                 DO m1 = 1, ldim
                    DO m2 = 1, ldim 
                       f(m1,m2) = rho%ns(m1,m2,is,na)
                    ENDDO
                 ENDDO
                 !
                 ! Diagonalize the ocupation matrix
                 CALL cdiagh(ldim, f, ldmx, lambda, vet)
                 !
                 ! Change the eigenvalues as requested from the input
                 DO i = 1, ldim
                    IF (starting_ns(i,is,nt) >= 0.d0) &
                            lambda(i) = starting_ns(i,is,nt)
                 ENDDO
                 !
                 ! Reconstruct back the occupation matrix from the
                 ! modified eignevalues and the original eigenvectors
                 DO m1 = 1, ldim
                    DO m2 = m1, ldim
                       temp = 0.d0
                       DO i = 1, ldim
                          temp = temp + CONJG(vet(m1,i))*lambda(i)*vet(m2,i)
                       ENDDO
                       rho%ns(m1,m2,is,na) =  DBLE(temp)
                       rho%ns(m2,m1,is,na) = rho%ns(m1,m2,is,na)
                    ENDDO
                 ENDDO
                 !
              ENDDO ! is
              !
           ENDIF ! noncolin
           !
        ENDIF ! is_hubbard
        !
     ENDDO ! na
     !
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
        ! Here we are copying the occupation matrix rho%ns to the 
        ! array oscdft_ctx%inp%occupation because it will be used in
        ! the constrained calculation
        CALL oscdft_ns_set (oscdft_ctx, Hubbard_lmax, Hubbard_l, rho%ns, 2)
     ENDIF
#endif
     !
     print_new_ns = .TRUE.
     !
   ELSE
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
        ! Here we are copying the target occupation matrix to the
        ! current/working occupation matrix rho%ns
        CALL oscdft_ns_set (oscdft_ctx, Hubbard_lmax, Hubbard_l, rho%ns, 1)
     ENDIF
     print_new_ns = .TRUE.
#endif
   ENDIF
   !
   ! Write the updated occupation matrices
   IF (print_new_ns) THEN
      IF (noncolin) THEN
         CALL write_ns_nc
      ELSE
         CALL write_ns
      ENDIF
   ENDIF
   !
   ! Reset starting_ns so that this step is not repeated
   starting_ns = -1.0_dp
   !
   RETURN
   !
END SUBROUTINE ns_adj
