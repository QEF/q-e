!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE nsg_adj
!-----------------------------------------------------------------------
   !
   ! This routine tries to suggest to the code the right atomic orbital to 
   ! localize the charge on.
   !
   USE kinds,            ONLY : DP
   USE ions_base,        ONLY : nat, ntyp => nsp, ityp
   USE ldaU,             ONLY : Hubbard_lmax, Hubbard_l, starting_ns, &
                                nsgnew, neighood, is_hubbard
   USE lsda_mod,         ONLY : nspin
   USE io_global,        ONLY : stdout
 
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: ldmx = 7
   INTEGER :: na, na1, nt, is, m1, m2, i, viz, ldim
   REAL(DP) :: lambda(ldmx)
   COMPLEX(DP) :: vet(ldmx,ldmx), f(ldmx,ldmx), temp
   !
   IF (ALL(starting_ns == -1.d0)) RETURN
   !
   WRITE(stdout,*) "Modify starting ns matrices according to input values"
   !
   IF (2*Hubbard_lmax+1 > ldmx) CALL errore('ns_adj',' ldmx is too small',ldmx)
   ! 
   DO na = 1, nat
      !
      nt = ityp(na)
      !
      IF (is_hubbard(nt)) THEN
         !
         ldim = 2*Hubbard_l(nt) + 1
         !
         DO is = 1, nspin
            !
            DO viz = 1, neighood(na)%num_neigh
               na1 = neighood(na)%neigh(viz)
               IF (na1.EQ.na) THEN
                  f(:,:) = (0.d0, 0.d0)
                  DO m1 = 1, ldim
                     DO m2 = 1, ldim
                        f(m1,m2) = nsgnew(m2,m1,viz,na,is)
                     ENDDO
                  ENDDO
                  GO TO 7
               ENDIF
            ENDDO
            !
7           CONTINUE
            !
            CALL cdiagh(ldim, f, ldmx, lambda, vet)
            !
            DO i = 1, ldim
              IF (starting_ns(i,is,nt) >= 0.d0) lambda(i) = starting_ns(i,is,nt)
            ENDDO
            !
            DO m1 = 1,ldim
               DO m2 = m1, ldim
                  temp = 0.d0
                  DO i = 1,ldim
                     temp = temp + CONJG(vet(m1,i))*lambda(i)*vet(m2,i)
                  ENDDO
                  nsgnew(m2,m1,viz,na,is) = DBLE(temp)
                  nsgnew(m1,m2,viz,na,is) = nsgnew(m2,m1,viz,na,is)
               ENDDO
            ENDDO
            !
         ENDDO
         !
      ENDIF
      !
   ENDDO
   !
   ! Uncomment the line below if needed (useful for debugging purposes)
   !CALL write_nsg
   !
   RETURN
   !
END SUBROUTINE nsg_adj
!-----------------------------------------------------------------------
