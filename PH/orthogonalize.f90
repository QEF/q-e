!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE orthogonalize(dvpsi, evq, ikk, ikq, dpsi)
!------------------------------------------------------------------------
  !
  ! This routine ortogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
  ! It should be quite general. It works for metals and insulators, with 
  ! NC as well as with US PP, both SR or FR.
  ! Note that on output it changes sign. So it applies -P^+_c.
  !
USE kinds, ONLY : DP
USE klist, ONLY : lgauss, degauss, ngauss
USE noncollin_module, ONLY : noncolin, npol
USE wvfct, ONLY : npwx, nbnd, et
USE ener, ONLY : ef
USE qpoint, ONLY : npwq
USE control_ph,  ONLY : alpha_pv, nbnd_occ
USE becmod,      ONLY : becp, becp_nc, calbec
USE uspp,        ONLY : vkb, okvan
USE mp_global,   ONLY : intra_pool_comm
USE mp,          ONLY : mp_sum

!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
COMPLEX(DP), INTENT(IN) :: evq(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbnd) ! work space allocated by
                                                   ! the calling routine

COMPLEX(DP), ALLOCATABLE :: ps(:,:)
INTEGER :: ibnd, jbnd, nbnd_eff
REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
REAL(DP), EXTERNAL :: w0gauss, wgauss
! functions computing the delta and theta function

CALL start_clock ('ortho')
ALLOCATE(ps(nbnd,nbnd))
!
if (lgauss) then
   !
   !  metallic case
   !
   ps = (0.d0, 0.d0)
   IF (noncolin) THEN
      CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
                 evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
   ELSE
      CALL ZGEMM( 'C', 'N', nbnd, nbnd_occ (ikk), npwq, (1.d0,0.d0), &
                 evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
   END IF
   !
   DO ibnd = 1, nbnd_occ (ikk)
      wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
      w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
      DO jbnd = 1, nbnd
         wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
         deltae = et (jbnd, ikq) - et (ibnd, ikk)
         theta = wgauss (deltae / degauss, 0)
         wwg = wg1 * (1.d0 - theta) + wgp * theta
         IF (jbnd <= nbnd_occ (ikq) ) THEN
            IF (abs (deltae) > 1.0d-5) THEN
                wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
            ELSE
               !
               !  if the two energies are too close takes the limit
               !  of the 0/0 ratio
               !
               wwg = wwg - alpha_pv * theta * w0g
            ENDIF
         ENDIF
         !
         ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
         !
      ENDDO
      IF (noncolin) THEN
         CALL DSCAL (2*npwx*npol, wg1, dvpsi(1,ibnd), 1)
      ELSE
         call DSCAL (2*npwq, wg1, dvpsi(1,ibnd), 1)
      END IF
   END DO
   nbnd_eff=nbnd
ELSE
   !
   !  insulators
   !
   ps = (0.d0, 0.d0)
   IF (noncolin) THEN
      CALL ZGEMM( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
             (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, &
             (0.d0,0.d0), ps, nbnd )
   ELSE
      CALL ZGEMM( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npwq, &
             (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
             (0.d0,0.d0), ps, nbnd )
   END IF
   nbnd_eff=nbnd_occ(ikk)
END IF
#ifdef __PARA
   call mp_sum(ps(:,1:nbnd_eff),intra_pool_comm)
#endif
!
! dpsi is used as work space to store S|evc>
!
IF (noncolin) THEN
   IF (okvan) CALL calbec ( npwq, vkb, evq, becp_nc, nbnd_eff )
ELSE
   IF (okvan) CALL calbec ( npwq, vkb, evq, becp, nbnd_eff)
ENDIF
CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)
!
! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
!
if (lgauss) then
   !
   !  metallic case
   !
   IF (noncolin) THEN
      CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
                (1.d0,0.d0), dpsi, npwx*npol, ps, nbnd, (-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSE
      CALL ZGEMM( 'N', 'N', npwq, nbnd_occ(ikk), nbnd, &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ELSE
   !
   !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
   !
   IF (noncolin) THEN
      CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
                (1.d0,0.d0),dpsi,npwx*npol,ps,nbnd,(-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSE
      CALL ZGEMM( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ENDIF

DEALLOCATE(ps)
CALL stop_clock ('ortho')
RETURN
END SUBROUTINE orthogonalize
