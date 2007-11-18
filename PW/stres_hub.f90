!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
#undef TIMING
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
                        hubbard_alpha, U_projection
  USE scf,       ONLY : v
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
#ifdef TIMING
   CALL start_clock( 'stres_hub' )
#endif

 
   IF (U_projection .NE. "atomic") CALL errore("stres_hub", &
                   " stress for this U_projection_type not implemented",1)

   IF (gamma_only) CALL errore('stres_hub',&
                   ' LDA+U, stress AND gamma-only not implemented yet',1)

   sigmah(:,:) = 0.d0

   ldim = 2 * Hubbard_lmax + 1
   ALLOCATE (dns(ldim,ldim,nspin,nat))
   dns(:,:,:,:) = 0.d0

#ifdef DEBUG
   DO na=1,nat
      DO is=1,nspin
         nt = ityp(na)
         IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN
            WRITE( stdout,'(a,2i3)') 'NS(NA,IS) ', na,is
            DO m1=1,ldim
               WRITE( stdout,'(7f10.4)') (rho%ns(m1,m2,is,na),m2=1,ldim)
            END DO
         END IF
      END DO
   END DO
#endif
   omin1 = 1.d0/omega
!
!  NB: both ipol and jpol must run from 1 to 3 because this stress 
!      contribution is not in general symmetric when computed only 
!      from k-points in the irreducible wedge of the BZ. 
!      It is (must be) symmetric after symmetrization but this requires 
!      the full stress tensor not only its upper triangular part.
!
   DO ipol = 1,3
      DO jpol = 1,3
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
                     DO m1 = 1, 2 * Hubbard_l(nt) + 1
                        sigmah(ipol,jpol) = sigmah(ipol,jpol) - omin1 * &
                           v%ns(m2,m1,is,na) * dns(m1,m2,is,na)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
   IF (nspin.EQ.1) sigmah(:,:) = 2.d0 * sigmah(:,:)

   CALL trntns(sigmah,at,bg,-1)
   CALL symtns(sigmah,nsym,s)
   CALL trntns(sigmah,at,bg,1)
!
! Symmetryze the stress tensor with respect to cartesian coordinates
! it should NOT be needed, let's do it for safety.
!
   DO ipol = 1,3
      DO jpol = ipol,3
         if ( abs( sigmah(ipol,jpol)-sigmah(jpol,ipol) )  > 1.d-6 ) then
             write (stdout,'(2i3,2f12.7)') ipol,jpol,sigmah(ipol,jpol), &
                                                     sigmah(jpol,ipol)
            call errore('stres_hub',' non-symmetric stress contribution',1)
         end if
         sigmah(ipol,jpol) = 0.5d0* ( sigmah(ipol,jpol) + sigmah(jpol,ipol) )
         sigmah(jpol,ipol) = sigmah(ipol,jpol)
      END DO
   END DO
   
   DEALLOCATE (dns)
#ifdef TIMING
   CALL stop_clock( 'stres_hub' )
   CALL print_clock( 'stres_hub' )
#endif

   RETURN
END  SUBROUTINE stres_hub
