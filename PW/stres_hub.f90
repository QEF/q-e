!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE stres_hub ( sigmah )
   !----------------------------------------------------------------------
   !
   ! This routines computes the Hubbard contribution to the internal stress
   ! tensor. It gives in output the array sigmah(i,j) which corresponds to
   ! the quantity -(1/\Omega)dE_{h}/d\epsilon_{i,j}
   !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : omega, at, bg
  USE ldaU,      ONLY : hubbard_lmax, hubbard_l, hubbard_u, &
                        hubbard_alpha, ns, U_projection
  USE lsda_mod,  ONLY : nspin
  USE symme,     ONLY : s, nsym
  USE io_files,  ONLY : prefix, iunocc
  USE wvfct,     ONLY : gamma_only   
  USE io_global, ONLY : stdout, ionode
   !
   IMPLICIT NONE
   !
   REAL (DP) :: sigmah(3,3)        ! output: the Hubbard stresses

   INTEGER :: ipol, jpol, na, nt, is,isi, m1,m2,m3,m4
   INTEGER :: ldim
   REAL (DP) :: omin1, current_sum, inverse_sum, sum, temp, flag
   LOGICAL :: exst
   REAL (DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat), ! the derivative of the atomic occupations
 
   IF (U_projection .NE. "atomic") CALL errore("stres_hub", &
                   " stress for this U_projection_type not implemented",1)

   IF (gamma_only) CALL errore('stres_hub',&
                   ' LDA+U, stress AND gamma-only not implemented yet',1)

   sigmah(:,:) = 0.d0

   ldim = 2 * Hubbard_lmax + 1
   ALLOCATE (dns(ldim,ldim,nspin,nat))
   dns(:,:,:,:) = 0.d0

   IF ( ionode ) THEN

      CALL seqopn(iunocc,'occup','formatted',exst)
      READ(iunocc,*) ns
      CLOSE(unit=iunocc,status='keep')

   END IF


#ifdef DEBUG
   DO na=1,nat
      DO is=1,nspin
         nt = ityp(na)
         IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN
            WRITE( stdout,'(a,2i3)') 'NS(NA,IS) ', na,is
            DO m1=1,ldim
               WRITE( stdout,'(7f10.4)') (ns(m1,m2,is,na),m2=1,ldim)
            END DO
         END IF
      END DO
   END DO
#endif
   omin1 = 1.d0/omega
   DO ipol = 1,3
      DO jpol = 1,ipol
         CALL dndepsilon(dns,ldim,ipol,jpol)
         DO na = 1,nat                 
            nt = ityp(na)
            IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN
               DO is = 1,nspin
#ifdef DEBUG
                  WRITE( stdout,'(a,4i3)') 'DNS(IPOL,JPOL,NA,IS) ', ipol,jpol,na,is
                  WRITE( stdout,'(5f10.4)') ((dns(m1,m2,is,na),m2=1,5),m1=1,5)
#endif
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1
                     sigmah(ipol,jpol) = sigmah(ipol,jpol) - omin1 * &
                           Hubbard_U(nt) * 0.5d0 * dns(m2,m2,is,na) 
                     DO m1 = 1, 2 * Hubbard_l(nt) + 1
                        sigmah(ipol,jpol) = sigmah(ipol,jpol) + omin1 * &
                           Hubbard_U(nt) * ns(m2,m1,is,na) * dns(m1,m2,is,na)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
   IF (nspin.EQ.1) sigmah(:,:) = 2.d0 * sigmah(:,:)

   !
   ! Symmetryze the stress tensor
   !
   DO ipol = 1,3
      DO jpol = ipol,3
         sigmah(ipol,jpol) = sigmah(jpol,ipol)
      END DO
   END DO
   
   CALL trntns(sigmah,at,bg,-1)
   CALL symtns(sigmah,nsym,s)
   CALL trntns(sigmah,at,bg,1)

   DEALLOCATE (dns)

   RETURN
END  SUBROUTINE stres_hub
