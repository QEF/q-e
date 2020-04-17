!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthogonalize2(dvpsi, evq, dpsi)
!------------------------------------------------------------------------
   !
   ! This routine ortogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
   ! It should be quite general. It works for metals and insulators, with
   ! NC as well as with US PP, both SR or FR.
   ! Note that on output it changes sign. So it applies -P^+_c.
   !
   ! NB: IN/OUT is dvpsi ; dpsi is used as work_space
   !
   USE kinds, ONLY: DP
   USE klist, ONLY: lgauss, degauss, ngauss
   USE noncollin_module, ONLY: noncolin, npol
   USE wvfct, ONLY: npwx, nbnd, et, npw
   USE ener, ONLY: ef
   USE control_lr, ONLY: alpha_pv, nbnd_occ
   USE becmod, ONLY: bec_type, becp, calbec
   USE uspp, ONLY: vkb, okvan
   USE mp_global, ONLY: intra_pool_comm
   USE mp, ONLY: mp_sum
   USE control_flags, ONLY: gamma_only
!USE realus,      ONLY : npw_k
   USE gvect, ONLY: gstart
!
   IMPLICIT NONE
!INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
!INTEGER, INTENT(IN) :: npwq       ! the number of plane waves for q
   COMPLEX(DP), INTENT(IN) :: evq(npwx*npol, nbnd)
   COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol, nbnd)
   COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol, nbnd) ! work space allocated by
   ! the calling routine

   COMPLEX(DP), ALLOCATABLE :: ps(:, :)
   REAL(DP), ALLOCATABLE :: ps_r(:, :)

   INTEGER :: ibnd, jbnd, nbnd_eff
   REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
   REAL(DP), EXTERNAL :: w0gauss, wgauss
   logical ::l_test
!
   integer ::npwq
! functions computing the delta and theta function

   CALL start_clock('ortho')

   npwq = npw
   IF (gamma_only) THEN
      ALLOCATE (ps_r(nbnd, nbnd))
      ps_r = 0.0_DP
   ENDIF

   ALLOCATE (ps(nbnd, nbnd))
   ps = (0.0_DP, 0.0_DP)

   IF (gamma_only) THEN
      CALL dgemm('C', 'N', nbnd, nbnd, 2*npwq, &
                 2.0_DP, evq, 2*npwx, dvpsi, 2*npwx, &
                 0.0_DP, ps_r, nbnd)
      IF (gstart == 2) THEN
         CALL DGER(nbnd, nbnd, -1.0_DP, evq, &
         & 2*npwq, dvpsi, 2*npwx, ps_r, nbnd)
      ENDIF
   END IF
#ifdef __MPI
   IF (gamma_only) THEN
      call mp_sum(ps_r(:, :), intra_pool_comm)
   END IF
#endif
!
   l_test = .false.
   if (l_test) then
      IF (gamma_only) THEN
         ps = CMPLX(ps_r, 0.0_DP, KIND=DP)
         CALL ZGEMM('N', 'N', npwq, nbnd, nbnd, &
                    (1.d0, 0.d0), dpsi, npwx, ps, nbnd, (-1.0d0, 0.d0), &
                    dvpsi, npwx)
      END IF
   else
      IF (gamma_only) THEN
         do ibnd = 1, nbnd
            dvpsi(1:npw, ibnd) = dvpsi(1:npw, ibnd) - evq(1:npw, ibnd)*ps_r(ibnd, ibnd)
         end do
      END IF
   end if

   DEALLOCATE (ps_r)
   CALL stop_clock('ortho')
   RETURN
END SUBROUTINE orthogonalize2
