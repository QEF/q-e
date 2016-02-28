!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine qdipol_cryst()
!
! This subroutine puts the dipole of Q on the crystal basis
!

USE kinds,     ONLY : DP
USE lsda_mod,  ONLY : nspin
USE uspp_param, ONLY : nh
USE spin_orb,  ONLY : lspinorb
USE cell_base, ONLY : at
USE ions_base, ONLY : nat, ityp, ntyp => nsp
USE lrus,      ONLY : dpqq, dpqq_so

IMPLICIT NONE

REAL(DP) :: fact(3)
COMPLEX(DP) :: fact_so(3)
INTEGER :: nt, na, ih, jh, ipol, is

DO nt = 1, ntyp
   DO ih = 1, nh (nt)
      DO jh = 1, nh (nt)
         IF (lspinorb) THEN
            DO is=1,nspin
               DO ipol=1,3
                  fact_so(ipol)=at(1,ipol)*dpqq_so(ih,jh,is,1,nt)+  &
                                at(2,ipol)*dpqq_so(ih,jh,is,2,nt)+  &
                                at(3,ipol)*dpqq_so(ih,jh,is,3,nt)
               ENDDO
               dpqq_so(ih,jh,is,:,nt)=fact_so(:)
            ENDDO
         END IF
         DO ipol=1,3
            fact(ipol)=at(1,ipol)*dpqq(ih,jh,1,nt)+  &
                       at(2,ipol)*dpqq(ih,jh,2,nt)+  &
                       at(3,ipol)*dpqq(ih,jh,3,nt)
         ENDDO
         dpqq(ih,jh,:,nt)=fact(:)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE qdipol_cryst
