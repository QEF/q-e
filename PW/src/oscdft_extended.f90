! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE oscdft_nsg (lflag)
   !
   !! This routine adjusts (modifies) the nsg based on constraints
   !! If lflag=1, then copy ctx%inp%occupation to nsgnew
   !! If lflag=2, then copy nsgnew to ctx%inp%occupation
   !! If lflag=3, then copy the diagonal components of the complex 
   !! If lflag=4, like lflag=1 but it does not nullify the occupations for 
   !!             Hubbard atoms atoms to which we do not apply the constraints
   !! generalized occupation matrix nsgnew to a real array nsnew that is used 
   !! to build the contraint 
   !
   USE kinds,           ONLY : DP
   USE parameters,      ONLY : ntypx
   USE io_global,       ONLY : stdout
   USE ions_base,       ONLY : nat, ityp
   USE lsda_mod,        ONLY : nspin
   USE upf_params,      ONLY : lqmax
   USE ldaU,            ONLY : max_num_neighbors, ldmx_tot, nsgnew, neighood, &
                               Hubbard_l, Hubbard_lmax, nsgnew, nsnew
#if defined (__OSCDFT)
   USE oscdft_base,     ONLY : oscdft_ctx
#endif
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: lflag
   !
   INTEGER :: na, na1, nt, viz, ldim, is, m1, m2
   LOGICAL :: found
   !
#if defined (__OSCDFT)
   IF (.NOT.(oscdft_ctx%inp%oscdft_type==2) .OR. .NOT.oscdft_ctx%is_constraint) RETURN
   !
   IF (lflag==1 .OR. lflag==2 .OR. lflag==4) &
   WRITE(stdout, '(/5x,"Modifying starting occupation matrices according to input constrained values")')
   !
   found = .true.
   !
   IF (lflag==1) nsgnew = (0.0d0, 0.0d0)
   !
   DO na = 1, nat
      IF (oscdft_ctx%constraining(na)) THEN
         nt = ityp(na)
         ldim = 2*Hubbard_l(nt) + 1
         DO is = 1, nspin
            DO viz = 1, neighood(na)%num_neigh
               na1 = neighood(na)%neigh(viz)
               IF (na1.EQ.na) THEN
                  DO m1 = 1, ldim
                     DO m2 = 1, ldim
                        IF (lflag==1 .OR. lflag==4) THEN
                           IF (oscdft_ctx%inp%occupation(m1,m2,is,na)==-2.d0) THEN
                              WRITE(stdout, '(/5x,"Warning!!! Missing element: ",4(1x,i4))') na, is, m1, m2
                              found = .false.
                           ELSE
                              nsgnew(m1,m2,viz,na,is) = oscdft_ctx%inp%occupation(m1,m2,is,na)
                           ENDIF
                        ELSEIF (lflag==2) THEN
                           oscdft_ctx%inp%occupation(m1,m2,is,na) = DBLE(nsgnew(m1,m2,viz,na,is))
                        ELSEIF (lflag==3) THEN
                           nsnew(m1,m2,is,na) = DBLE(nsgnew(m1,m2,viz,na,is))
                        ENDIF
                     ENDDO
                  ENDDO
                  GO TO 7
               ENDIF
            ENDDO
7           CONTINUE
         ENDDO
      ENDIF
   ENDDO
   !
   IF (.NOT.found) CALL errore( 'oscdft_nsg', 'Missing target occupation matrix element', 1 )
   !
#endif
   RETURN
   !
END SUBROUTINE oscdft_nsg

SUBROUTINE oscdft_v_constraint_extended (nsg, vtot, etot)
   !
   !! Computes the contribution to the potential from the Lagrange multipliers used
   !! to constrain the occupation matrix to the target (DFT+U+V case).
   !! Here we are applying the constraint only using the diagonal (onsite) component 
   !! of the generalized occupaion matrix nsg.
   !
   USE kinds,           ONLY : DP
   USE parameters,      ONLY : ntypx
   USE ions_base,       ONLY : nat, ityp
   USE lsda_mod,        ONLY : nspin
   USE io_global,       ONLY : stdout
   USE control_flags,   ONLY : iverbosity
   USE ldaU,            ONLY : ldim_u, ldmx_tot, max_num_neighbors, neighood 
#if defined (__OSCDFT)
   USE oscdft_base,     ONLY : oscdft_ctx
#endif

   IMPLICIT NONE

   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(INOUT) :: vtot(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   REAL(DP), INTENT(INOUT) :: etot
   REAL(DP) :: ec
   INTEGER :: is, na, na1, na2, viz, nt1, m1, m2
   INTEGER, EXTERNAL :: find_viz
   !
#if defined (__OSCDFT)
   IF (.NOT.(oscdft_ctx%inp%oscdft_type==2) .OR. .NOT.oscdft_ctx%is_constraint .OR. oscdft_ctx%conv) RETURN
   !
   ec = 0
   !
   DO na1 = 1, nat
      IF (oscdft_ctx%constraining(na1)) THEN
         nt1 = ityp(na1)
         DO is = 1, nspin
            DO viz = 1, neighood(na1)%num_neigh
               na2 = neighood(na1)%neigh(viz)
               IF ( na1.EQ.na2 ) THEN
                  DO m1 = 1, ldim_u(nt1)
                     DO m2 = 1, ldim_u(nt1)
                        vtot(m1,m2,viz,na1,is) = vtot(m1,m2,viz,na1,is) + &
                                                 oscdft_ctx%inp%constraint_strength * &
                                                 oscdft_ctx%constraint(m1,m2,is,na1)
                        ec = ec + oscdft_ctx%inp%constraint_strength * &
                                  oscdft_ctx%constraint(m1,m2,is,na1) * &
                                  (DBLE(nsg(m1,m2,viz,na1,is)) - oscdft_ctx%inp%occupation(m1,m2,is,na1))
                     ENDDO 
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
   ENDDO 
   !
   IF (nspin == 1) THEN
       ec = 2.d0 * ec
   ENDIF
   etot = etot + ec
   !
   IF (iverbosity > 0) THEN
      WRITE(stdout, '(/5x,"CONSTRAINT ENERGY = ", f9.4, 1x," (Ry)")') ec    
   ENDIF
   !
#endif
   RETURN
   !
END SUBROUTINE oscdft_v_constraint_extended
