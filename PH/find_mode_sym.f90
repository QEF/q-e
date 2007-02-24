!
! Copyright (C) 2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym (dyn, w2, at, bg, nat, nsym, s, irt, xq, rtau, &
           amass, ntyp, ityp)
!
!   This subroutine finds the irreducible representations which give
!   the transformation properties of eigenvectors of the dynamical 
!   matrix. It does NOT work at zone border in non symmorphic space groups.
!  
!
#include "f_defs.h"
USE io_global,  ONLY : stdout
USE kinds, ONLY : DP
USE noncollin_module, ONLY : noncolin
USE spin_orb, ONLY : domag
USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
                            char_mat, name_rap, name_class, gname, ir_ram
USE rap_point_group_is, ONLY : gname_is
USE control_ph, ONLY : lgamma
IMPLICIT NONE
INTEGER ::                  &
          nat, nsym,        & 
          ntyp, ityp(nat),  &
          irt(48,nat),      &
          s(3,3,48)  

REAL(DP) ::                 &
            at(3,3),        &
            bg(3,3),        &
            xq(3),          &
            rtau(3,48,nat), &
            amass(ntyp),    &
            w2(3*nat)       

COMPLEX(DP) ::  &
            dyn(3*nat, 3*nat)       

REAL(DP), PARAMETER :: eps=1.d-5,  &
                       rydcm1 = 13.6058d0 * 8065.5d0

INTEGER ::      &
        ngroup, &   ! number of different frequencies groups
        nmodes, &   ! number of modes
        imode, igroup, dim_rap, nu_i, nu_j, irot, irap, iclass, mu, na, i, j

INTEGER, ALLOCATABLE :: istart(:)

REAL(DP) :: sr(3,3,48)
COMPLEX(DP) :: ZDOTC, times              ! safe dimension 
                                         ! in case of accidental degeneracy 
REAL(DP), ALLOCATABLE :: w1(:)
COMPLEX(DP), ALLOCATABLE ::  rmode(:), trace(:,:), z(:,:)
CHARACTER(3) :: cdum
!
!    Divide the modes on the basis of the mode degeneracy.
!
nmodes=3*nat

ALLOCATE(istart(nmodes+1))
ALLOCATE(z(nmodes,nmodes))
ALLOCATE(w1(nmodes))
ALLOCATE(rmode(nmodes))
ALLOCATE(trace(48,nmodes))

DO nu_i = 1, nmodes
   DO mu = 1, nmodes
      na = (mu - 1) / 3 + 1
      z (mu, nu_i) = dyn (mu, nu_i) * SQRT (amass (ityp (na) ) )
   END DO
END DO

DO imode=1,nmodes
   w1(imode)=SIGN(SQRT(ABS(w2(imode)))*rydcm1,w2(imode))
ENDDO

DO irot=1,nsym
   CALL s_axis_to_cart (s(1,1,irot), sr(1,1,irot), at, bg)
END DO

ngroup=1
istart(ngroup)=1
DO imode=2,nmodes
   IF (ABS(w1(imode)-w1(imode-1)) > 5.0d-2) THEN
      ngroup=ngroup+1
      istart(ngroup)=imode
   END IF
END DO
istart(ngroup+1)=nmodes+1
!
!  Find the character of one symmetry operation per class
!
DO igroup=1,ngroup
   dim_rap=istart(igroup+1)-istart(igroup)
   DO iclass=1,nclass
      irot=elem(1,iclass)
      trace(iclass,igroup)=(0.d0,0.d0)
      DO i=1,dim_rap
         nu_i=istart(igroup)+i-1
         CALL rotate_mod(z(1,nu_i),rmode,sr(1,1,irot),irt,rtau,xq,nat,irot)
         trace(iclass,igroup)=trace(iclass,igroup) + &
                            ZDOTC(3*nat,z(1,nu_i),1,rmode,1)
      END DO
!      write(6,*) igroup,iclass, trace(iclass,igroup)
   END DO
END DO
!
!  And now use the character table to identify the symmetry representation
!  of each group of modes
!
IF (noncolin.and.domag) THEN
   WRITE(stdout,  &
    '(/,5x,"Mode symmetry, ",a11," [",a11,"] magnetic point group:",/)') &
                   gname, gname_is
ELSE
   WRITE(stdout,'(/,5x,"Mode symmetry, ",a11," point group:",/)') gname
END IF

DO igroup=1,ngroup
   DO irap=1,nclass
      times=(0.d0,0.d0)
      DO iclass=1,nclass
         times=times+CONJG(trace(iclass,igroup))*char_mat(irap, &
                     which_irr(iclass))*nelem(iclass)
!         write(6,*) igroup, irap, iclass, which_irr(iclass)
      ENDDO
      times=times/nsym
      cdum="   "
      IF (lgamma) cdum=ir_ram(irap)
      IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
          (ABS(AIMAG(times)) > eps) ) THEN
            WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup))
      ENDIF

      IF (ABS(times) > eps) THEN
         IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
            WRITE(stdout,'(5x, "omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",a19)') &
              istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                                name_rap(irap)//" "//cdum
         ELSE
            WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",i3,a19)') &
              istart(igroup), istart(igroup+1)-1, &
              w1(istart(igroup)), NINT(DBLE(times)), &
                                  name_rap(irap)//" "//cdum
         END IF
      END IF
   END DO
END DO
WRITE( stdout, '(/,1x,74("*"))')

DEALLOCATE(trace)
DEALLOCATE(z)
DEALLOCATE(w1)
DEALLOCATE(rmode)
DEALLOCATE(istart)

RETURN
END SUBROUTINE find_mode_sym

SUBROUTINE rotate_mod(mode,rmode,sr,irt,rtau,xq,nat,irot)
USE kinds, ONLY : DP
USE constants, ONLY: tpi
IMPLICIT NONE

INTEGER :: nat, irot, irt(48,nat)
COMPLEX(DP) :: mode(3*nat), rmode(3*nat), phase
REAL(DP)  :: sr(3,3), rtau(3,48,nat), xq(3), arg
INTEGER :: na, nb, ipol, kpol, mu_i, mu_k

rmode=(0.d0,0.d0)
DO na=1,nat
   nb=irt(irot,na)
   arg = ( xq(1)*rtau(1,irot,na) + xq(2)*rtau(2,irot,na)+  &
           xq(3)*rtau(3,irot,na) ) * tpi
   phase = cmplx(cos(arg), sin(arg))
   DO ipol=1,3
      mu_i=3*(na-1)+ipol
      DO kpol=1,3
         mu_k=3*(nb-1)+kpol
         rmode(mu_i)=rmode(mu_i) + sr(kpol,ipol)*mode(mu_k)*phase
      END DO
   END DO
END DO

RETURN
END SUBROUTINE rotate_mod
