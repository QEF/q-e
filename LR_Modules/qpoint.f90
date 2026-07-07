!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
! TODO: Remove qpoint_aux module.
!
! TODO: Use qpoint_setup_k_plus_q_indices in HP/src/hp_load_q.f90, KCW/src/kcw_initialize_ph.f90,
!       and PHonon/PH/initialize_ph.f90
!
!       USE qpoint,           ONLY : qpoint_setup_k_plus_q_indices
!       !
!       ! Setup k and k+q point indices (also for -k and -k-q for magnetic calculations)
!       !
!       CALL qpoint_setup_k_plus_q_indices()
!
! TODO: Remove xq = xk(ikq) - xk(ikk) check in HP/src/hp_init_q.f90, KCW/src/kcw_init_q.f90,
!       and PHonon/PH/phq_init.f90
!
! TODO: Use qpoint_setup_eigqts in HP/src/hp_init_q.f90, KCW/src/kcw_init_q.f90, PHonon/PH/phq_init.f90
!
!       USE qpoint,               ONLY : qpoint_setup_eigqts
!       ! 1) USPP: Compute the phase factor exp(-i q*\tau)
!       !    The result is stored in variable eigqts in module qpoints.
!       !
!       IF (okvan) CALL qpoint_setup_eigqts()
!
! TODO: Remove PHonon/PH/phcom.f90 MODULE disp, use qpoint instead.
!
MODULE qpoint
   !---------------------------------------------------------------------
   !! The variables needed to specify various indices,
   !! number of plane waves and k and k+q points and their coordiantes.
   !!
   !! In a LR calculation with finite q, we deal with k and k+q points.
   !! For magnetic calculations, we also need the -k and -k-q points.
   !! The number of the original, "real" k points in this core is nksq, whereas the
   !! number of full k points, including k+q, -k, and -k-q, is nks (in module klist).
   !! nksqtot and nkstot are the sum of nksq and nks over pools (k parallelization).
   !!
   !! For ik = 1, ..., nksq, the index of k, k+q, -k, and -k-q in the full list of k points
   !! is given by ikks(ik), ikqs(ik), ikmks(ik), and ikmkmqs(ik) respectively.
   !! For example, a typical loop over k points in a LR calculation with finite q looks like:
   !!     DO ik = 1, nksq
   !!        ikk = ikks(ik)
   !!        ikq = ikqs(ik)
   !!        npw = ngk(ikk)
   !!        npwq= ngk(ikq)
   !!        ...
   !!     ENDDO
   !!
   !! Use subroutine qpoint_setup_k_plus_q_points to initialize these parameters.
   !
   USE kinds,      ONLY : DP
   !
   SAVE
   !
   ! Information about the k points (k+q, -k, and -k-q points as well)
   !
   INTEGER, POINTER :: igkq(:)     ! npwx)
   ! correspondence k+q+G <-> G
   INTEGER :: nksq
   !! the real number of k points
   INTEGER :: nksqtot
   !! the total number of q points
   INTEGER, ALLOCATABLE :: ikks(:)
   !! the index of k in the full list of k points
   INTEGER, ALLOCATABLE :: ikqs(:)
   !! the index of k+q in the full list of k points
   INTEGER, ALLOCATABLE :: ikmks(:)
   !! the index of -k in the full list of k points (for magnetic calculations)
   INTEGER, ALLOCATABLE :: ikmkmqs(:)
   !! the index of -k-q in the full list of k points (for magnetic calculations)
   !
   INTEGER :: npwq
   !! the number of plane waves for k+q
   !! FIXME: In most cases, npwq is defined locally. In a few cases, it is used from this module.
   !! FIXME: Make this consistent (by removing npwq from this module).
   !! FIXME: The same problem is with npw in module klist.
   !
   REAL (DP) :: xq(3)
   ! the coordinates of the q point
   COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
   ! the phases associated to the q
   REAL (DP), ALLOCATABLE :: xk_col(:,:)
   !
   ! Information about the list of q points
   !
   LOGICAL, ALLOCATABLE :: lgamma_iq(:)
   !! if TRUE this q is gamma.
   LOGICAL, ALLOCATABLE :: done_iq(:)
   !! if TRUE this q point has been already calculated
   LOGICAL, ALLOCATABLE :: comp_iq(:)
   !! if TRUE this q point has to be calculated
   INTEGER :: nq1, nq2, nq3
   !! Number of q points in each direction
   INTEGER :: nqs
   !! Number of q points to be calculated
   INTEGER :: start_q
   !! Initial q point to be computed
   INTEGER :: last_q
   !! Final q point to be computed
   REAL(DP), ALLOCATABLE :: x_q(:, :)
   !! Coordinates of q points (in Cartesian, tpiba units)
   REAL(DP), ALLOCATABLE :: wq(:)
   !! for plotting
   !
CONTAINS
   !
   !------------------------------------------------------------------------------
   SUBROUTINE qpoint_setup_k_plus_q_indices()
      !----------------------------------------------------------------------------
      !! This subroutine sets up the k+q points and their indices.
      !! It is called at the beginning of the LR calculation.
      !! It allocates the arrays ikks, ikqs, ikmks, and ikmkmqs.
      !----------------------------------------------------------------------------
      !
      USE constants,        ONLY : eps8
      USE io_global,        ONLY : stdout
      USE klist,            ONLY : nks, nkstot, xk
      USE control_lr,       ONLY : lgamma
      USE noncollin_module, ONLY : noncolin, domag
      !
      IMPLICIT NONE
      !
      INTEGER :: ik, ikk, ikq, ikmk, ikmkmq, ipol
      !
      IF ( lgamma ) THEN
         !
         IF (noncolin .AND. domag) THEN
            ! q = 0, magnetic. Need k and -k points. (k+q is the same as k, -k-q is the same as -k)
            nksq = nks/2
            nksqtot = nkstot/2
            ALLOCATE(ikks(nksq))
            ALLOCATE(ikqs(nksq))
            ALLOCATE(ikmks(nksq))
            ALLOCATE(ikmkmqs(nksq))
            DO ik = 1, nksq
               ikks(ik)    = 2 * ik - 1
               ikqs(ik)    = 2 * ik - 1
               ikmks(ik)   = 2 * ik
               ikmkmqs(ik) = 2 * ik
            ENDDO
            !
         ELSE
            ! q = 0, nonmagnetic. Need only k points. (k+q is the same as k)
            nksq = nks
            nksqtot = nkstot
            ALLOCATE(ikks(nksq))
            ALLOCATE(ikqs(nksq))
            DO ik = 1, nksq
               ikks(ik) = ik
               ikqs(ik) = ik
            ENDDO
         ENDIF
         !
      ELSE
         !
         IF (noncolin .AND. domag) THEN
            ! q /= 0, magnetic. Need k, k+q, -k, and -k-q points.
            nksq = nks / 4
            nksqtot = nkstot / 4
            ALLOCATE(ikks(nksq))
            ALLOCATE(ikqs(nksq))
            ALLOCATE(ikmks(nksq))
            ALLOCATE(ikmkmqs(nksq))
            DO ik = 1, nksq
               ikks(ik)    = 4 * ik - 3
               ikqs(ik)    = 4 * ik - 2
               ikmks(ik)   = 4 * ik - 1
               ikmkmqs(ik) = 4 * ik
            ENDDO
            !
         ELSE
            ! q /= 0, nonmagnetic. Need k and k+q points.
            nksq = nks / 2
            nksqtot = nkstot / 2
            ALLOCATE(ikks(nksq))
            ALLOCATE(ikqs(nksq))
            DO ik = 1, nksq
               ikks(ik) = 2 * ik - 1
               ikqs(ik) = 2 * ik
            ENDDO
         ENDIF
         !
      ENDIF
      !
      ! Check order of k points
      !
      DO ik = 1, nksq
         !
         ! Check if xk(ikq) = xk(ikk) + xq
         !
         ikk = ikks(ik)
         ikq = ikqs(ik)
         !
         IF ( ANY( ABS( xk(:, ikq) - xk(:, ikk) - xq ) > eps8 )) THEN
            WRITE(stdout, '(/,5x,"k points #", i6, " and ", i6, 5x," total number ", i6)') ikk, ikq, nksq
            WRITE(stdout, '(  5x,"Current q    point  ", 3f10.7)') (xq(ipol), ipol = 1, 3)
            WRITE(stdout, '(  5x,"Found k      point  ", 3f10.7)') (xk(ipol,ikk), ipol = 1, 3)
            WRITE(stdout, '(  5x,"Found k+q    point  ", 3f10.7)') (xk(ipol,ikq), ipol = 1, 3)
            WRITE(stdout, '(  5x,"Expected k+q point  ", 3f10.7)') (xk(ipol,ikk) + xq(ipol), ipol = 1, 3)
            CALL errore('qpoint_setup_k_plus_q_indices', 'wrong order of k+q point', 1)
         ENDIF
         !
         IF (noncolin .AND. domag) THEN
            !
            ! Check if xk(ikmq) = -xk(ikk) and xk(ikmkmq) = -xk(ikk) - xq
            !
            ikmk = ikmks(ik)
            ikmkmq = ikmkmqs(ik)
            !
            IF ( ANY( ABS( xk(:, ikmk) + xk(:, ikk) ) > eps8 )) THEN
               WRITE(stdout, '(/,5x,"k points #", i6, " and ", i6, 5x," total number ", i6)') ikk, ikmk, nksq
               WRITE(stdout, '(  5x,"Current q   point  ", 3f10.7)') (xq(ipol), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Found k     point  ", 3f10.7)') (xk(ipol,ikk), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Found -k    point  ", 3f10.7)') (xk(ipol,ikmk), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Expected -k point  ", 3f10.7)') (-xk(ipol,ikk), ipol = 1, 3)
               CALL errore('qpoint_setup_k_plus_q_indices', 'wrong order of -k point', 1)
            ENDIF
            !
            IF ( ANY( ABS( xk(:, ikmkmq) + xk(:, ikk) + xq ) > eps8 )) THEN
               WRITE(stdout, '(/,5x,"k points #", i6, " and ", i6, 5x," total number ", i6)') ikk, ikmkmq, nksq
               WRITE(stdout, '(  5x,"Current q      point  ", 3f10.7)') (xq(ipol), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Found k        point  ", 3f10.7)') (xk(ipol,ikk), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Found -k-q     point  ", 3f10.7)') (xk(ipol,ikmkmq), ipol = 1, 3)
               WRITE(stdout, '(  5x,"Expected -k-q  point  ", 3f10.7)') (-xk(ipol,ikk) - xq(ipol), ipol = 1, 3)
               CALL errore('qpoint_setup_k_plus_q_indices', 'wrong order of -k-q point', 1)
            ENDIF
            !
         ENDIF ! noncolin .AND. domag
         !
      ENDDO ! ik
      !
   END SUBROUTINE qpoint_setup_k_plus_q_indices
   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------
   SUBROUTINE qpoint_setup_eigqts()
      !----------------------------------------------------------------------------
      !! Compute the phase factor eigqts = exp(-i q*\tau)
      !----------------------------------------------------------------------------
      !
      USE kinds,            ONLY : DP
      USE constants,        ONLY : tpi
      USE ions_base,        ONLY : nat, tau
      !
      IMPLICIT NONE
      !
      INTEGER :: na
      !! Atom index
      REAL(DP) :: arg
      !! argument of the phase factor exponent
      !
      DO na = 1, nat
         !
         arg = (  xq(1) * tau(1,na) + &
                  xq(2) * tau(2,na) + &
                  xq(3) * tau(3,na) ) * tpi
         !
         eigqts(na) = CMPLX(COS(arg), - SIN(arg), KIND = DP)
         !
      ENDDO ! na
      !
   END SUBROUTINE qpoint_setup_eigqts
   !--------------------------------------------------------------------------------
   !
END MODULE qpoint


MODULE qpoint_aux
   ! TODO: This module is deprecated and should be removed. It is here only to avoid
   ! TODO: breaking PH, HP, and KCW codes that use it.
   USE qpoint, ONLY : ikmks, ikmkmqs
   USE lrus, ONLY : becpt, alphapt
END MODULE qpoint_aux
