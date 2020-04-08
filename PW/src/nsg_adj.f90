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
   USE ldaU,             ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, starting_ns, &
                                Hubbard_l_back, Hubbard_U_back, starting_ns_back, ldim_u, &
                                nsg, nsgnew
   USE scf,              ONLY : rho
   USE lsda_mod,         ONLY : nspin
   USE noncollin_module, ONLY : noncolin, npol
   USE io_global,        ONLY : stdout
 
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: ldmx = 7
   INTEGER :: na, nt, is, m1, m2, majs, mins, adjs, mol(ldmx), &
              nel, i, j, l, index(ldmx), viz, ldim
   REAL(DP) :: totoc, delta,lambda(ldmx)
   COMPLEX(DP) :: vet(ldmx,ldmx), f(ldmx,ldmx), temp
   LOGICAL :: adjust
   INTEGER, EXTERNAL :: find_viz
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
      ldim = 2*Hubbard_l(nt) + 1
      !
      viz = find_viz(na,na)
      !
      IF (ldim_u(nt).GT.0) THEN
         !
         DO is = 1, nspin
            !
            DO m1 = 1, ldim
               DO m2 = 1, ldim
                 f(m1,m2) = nsgnew(na,m1,viz,m2,is)
               ENDDO
            ENDDO
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
                  nsgnew(na,m1,viz,m2,is) = DBLE(temp)
                  nsgnew(na,m2,viz,m1,is) = nsgnew(na,m1,viz,m2,is)
               ENDDO
            ENDDO
            !
         ENDDO
         !
      ENDIF
      !
   ENDDO
   !
   RETURN
   !
END SUBROUTINE nsg_adj
!-----------------------------------------------------------------------
