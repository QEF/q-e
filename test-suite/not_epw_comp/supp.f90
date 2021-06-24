!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE el_phon2
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  !
  LOGICAL :: elph_epw
  REAL(DP) :: kx, ky, kz
  !
END MODULE el_phon2

!------------------------------------------------------------------------
SUBROUTINE elphsum2( )
  !-----------------------------------------------------------------------
  !
  ! Print the |g| vertex for all n,n' and modes in meV and do average
  ! on degenerate states.
  !
  ! Rewritten by H. Lee based on the previous version
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : rytoev
  USE ions_base,   ONLY : nat
  USE klist,       ONLY : xk
  USE wvfct,       ONLY : nbnd, et
  USE el_phon,     ONLY : el_ph_mat
  USE modes,       ONLY : u, nmodes
  USE dynmat,      ONLY : dyn,w2
  USE io_global,   ONLY : stdout, ionode
  USE mp,          ONLY : mp_sum
  USE mp_pools,    ONLY : npool
  USE mp_images,   ONLY : intra_image_comm, me_image, nproc_image, root_image
  USE mp_bands,    ONLY : root_bgrp, me_bgrp
  USE qpoint,      ONLY : nksq, ikks, ikqs, xq
  USE el_phon2,    ONLY : kx, ky, kz
  USE parallel_include
  !
  IMPLICIT NONE
  !
  LOGICAL :: found
  !
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, pbnd, nu, mu, vu, ierr, istatus
  INTEGER :: nksq2, ikk2, ikq2, nkq2, ik1, ik2, ipert, jpert, n
  !
  REAL(DP), PARAMETER :: ryd2mev  = rytoev * 1.0E3_DP
  REAL(DP), PARAMETER :: eps = 0.01/ryd2mev
  REAL(DP) :: gamma, g2, w, w_1, w_2
  REAL(DP) :: kpoint(3)
  REAL(DP) :: epc(nbnd, nbnd, 3 * nat)
  REAL(DP) :: epc_sym(nbnd, nbnd, 3 * nat)
  REAL(DP), ALLOCATABLE :: et2(:, :)
  !
  COMPLEX(DP) :: el_ph_sum(3 * nat, 3 * nat)
  COMPLEX(DP) :: el_ph_sum_aux(3 * nat, 3 * nat)
  COMPLEX(DP), ALLOCATABLE :: el_ph_mat2(:, :, :, :)
  !
  CALL start_clock('elphsum2')
  !
  kpoint = (/ kx, ky, kz /)
  !
  WRITE(stdout, '(5x,/"electron-phonon interaction  ..."/)')
  !
  found = .FALSE.
  DO ik = 1, nksq
    IF (ANY(ABS(xk(:, ikks(ik))-kpoint(:)) > 1.0E-6_DP)) CYCLE
    found = .TRUE.
    ik1 = ik
    ikk = ikks(ik)
    ikq = ikqs(ik)
  ENDDO
  !
  ierr = 0
  IF (found) THEN
    ierr = 1
  ENDIF
  CALL mp_sum(ierr, intra_image_comm)
  IF (ierr == 0) CALL errore('elphsum2', 'kpoint not found', 1)  
  !
#if defined(__MPI)
  IF ((npool == 1) .AND. ionode) THEN
    ik2 = ik1
    ikk2 = ikk
    ikq2 = ikq
    nksq2 = nksq
    ALLOCATE(el_ph_mat2(nbnd, nbnd, nksq2, 3 * nat))
    el_ph_mat2 = el_ph_mat
    nkq2 = SIZE(et, 2)
    ALLOCATE(et2(nbnd, nkq2))
    et2 = et
  ELSE
    IF (found .AND. (me_bgrp == root_bgrp) .AND. (.NOT. ionode)) THEN
      CALL MPI_SEND(ik1, 1, MPI_INTEGER, root_image, 0, intra_image_comm, ierr)
      CALL MPI_SEND(ikk, 1, MPI_INTEGER, root_image, nproc_image, &
                    intra_image_comm, ierr)
      CALL MPI_SEND(ikq, 1, MPI_INTEGER, root_image, 2 * nproc_image, &
                    intra_image_comm, ierr)
      CALL MPI_SEND(nksq, 1, MPI_INTEGER, root_image, 3 * nproc_image, &
                    intra_image_comm, ierr)
      CALL MPI_SEND(el_ph_mat, nbnd * nbnd * nksq * 3 * nat, MPI_DOUBLE_COMPLEX, &
                    root_image, 4 * nproc_image, intra_image_comm, ierr)
      CALL MPI_SEND(SIZE(et, 2), 1, MPI_INTEGER, root_image, 5 * nproc_image, &
                    intra_image_comm, ierr)
      CALL MPI_SEND(et, nbnd * SIZE(et, 2), MPI_DOUBLE_PRECISION, root_image, &
                    6 * nproc_image, intra_image_comm, ierr)
    ELSEIF (ionode) THEN
      IF (.NOT. found) THEN
        CALL MPI_RECV(ik2, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, intra_image_comm, &
                      istatus, ierr )
        CALL MPI_RECV(ikk2, 1, MPI_INTEGER, MPI_ANY_SOURCE, nproc_image, &
                      intra_image_comm, istatus, ierr )
        CALL MPI_RECV(ikq2, 1, MPI_INTEGER, MPI_ANY_SOURCE, 2 * nproc_image, &
                      intra_image_comm, istatus, ierr )
        CALL MPI_RECV(nksq2, 1, MPI_INTEGER, MPI_ANY_SOURCE, 3 * nproc_image, &
                      intra_image_comm, istatus, ierr )
        ALLOCATE(el_ph_mat2(nbnd, nbnd, nksq2, 3 * nat))
        CALL MPI_RECV(el_ph_mat2, nbnd * nbnd * nksq2 * 3 * nat, MPI_DOUBLE_COMPLEX, &
                      MPI_ANY_SOURCE, 4 * nproc_image, intra_image_comm, istatus, ierr ) 
        CALL MPI_RECV(nkq2, 1, MPI_INTEGER, MPI_ANY_SOURCE, 5 * nproc_image, &
                      intra_image_comm, istatus, ierr )
        ALLOCATE(et2(nbnd, nkq2))
        CALL MPI_RECV(et2, nbnd * nkq2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      6 * nproc_image, intra_image_comm, istatus, ierr ) 
      ELSE
        ik2 = ik1
        ikk2 = ikk
        ikq2 = ikq
        nksq2 = nksq
        ALLOCATE(el_ph_mat2(nbnd, nbnd, nksq2, 3 * nat))
        el_ph_mat2 = el_ph_mat
        nkq2 = SIZE(et, 2)
        ALLOCATE(et2(nbnd, nkq2))
        et2 = et
      ENDIF
    ENDIF
  ENDIF
#else
  ik2 = ik1
  ikk2 = ikk
  ikq2 = ikq
  nksq2 = nksq
  ALLOCATE(el_ph_mat2(nbnd, nbnd, nksq2, 3 * nat))
  el_ph_mat2 = el_ph_mat
  nkq2 = SIZE(et, 2)
  ALLOCATE(et2(nbnd, nkq2))
  et2 = et
#endif
  !
  IF (ionode) THEN
    !
    DO ibnd = 1, nbnd
      DO jbnd = 1, nbnd
        !
        DO jpert = 1, 3 * nat
          DO ipert = 1, 3 * nat
            el_ph_sum(ipert, jpert) = CONJG(el_ph_mat2(jbnd, ibnd, ik2, ipert)) * &
                                            el_ph_mat2(jbnd, ibnd, ik2, jpert)
          ENDDO
        ENDDO
        !
        ! from pert to cart
        !
        CALL dyn_pattern_to_cart(nat, u, el_ph_sum, el_ph_sum_aux)
        CALL compact_dyn(nat, el_ph_sum, el_ph_sum_aux)
        !
        DO nu = 1, nmodes
          gamma = 0.d0
          DO vu = 1, 3 * nat
            DO mu = 1, 3 * nat
              gamma = gamma + REAL(CONJG(dyn(mu, nu)) * el_ph_sum(mu, vu) &
                      * dyn(vu, nu))
            ENDDO
          ENDDO
          gamma = gamma / 2.d0
          !
          ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that
          ! appears
          ! in the definition of the electron-phonon matrix element g
          ! The sqrt(1/M) factor is actually hidden into the normal modes
          ! we still need to divide by the phonon frequency in Ry
          !
          IF (w2(nu) .GT. 0.d0) THEN
            w = DSQRT(w2(nu))
            gamma = gamma / w
          ELSE
            w = DSQRT(-w2(nu))
            gamma = 0.d0
          ENDIF
          !
          IF (gamma .LT. 0.d0) gamma = 0.d0
          !
          gamma = DSQRT(gamma)
          !
          ! gamma = |g| [Ry]
          !
          epc(ibnd, jbnd, nu) = gamma
          !
        ENDDO
        !
      ENDDO
    ENDDO
    !
    ! HERE WE "SYMMETRIZE": actually we simply take the averages over
    ! degenerate states, it is only a convention because g is gauge-dependent!
    !
    ! first the phonons 
    !
    DO jbnd = 1, nbnd
      DO ibnd = 1, nbnd
        !
        DO nu = 1, nmodes
          !
          w_1 = DSQRT(ABS(w2(nu)))
          g2 = 0.d0
          n = 0
          !
          DO mu = 1, nmodes
            !
            w_2 = DSQRT(ABS(w2(mu)))
            !
            IF (ABS(w_2 - w_1) .LT. eps) THEN
              n = n + 1
              g2 = g2 + epc(ibnd, jbnd, mu) * epc(ibnd, jbnd, mu)
            ENDIF
            !
          ENDDO 
          !
          g2 = g2 / FLOAT(n)
          epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          !
        ENDDO
        !
      ENDDO
    ENDDO
    epc = epc_sym
    !
    ! then the k electrons
    !
    DO nu  = 1, nmodes
      DO jbnd = 1, nbnd
        !
        DO ibnd = 1, nbnd
          !
          w_1 = et2(ibnd, ikk2)
          g2 = 0.d0
          n  = 0
          !
          DO pbnd = 1, nbnd
            !
            w_2 = et2(pbnd, ikk2)
            !
            IF (ABS(w_2 - w_1) .LT. eps) THEN
              n = n + 1
              g2 = g2 + epc(pbnd, jbnd, nu) * epc(pbnd, jbnd, nu)
            ENDIF
            !
          ENDDO 
          !
          g2 = g2 / FLOAT(n)
          epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          !
        ENDDO
        !
      ENDDO
    ENDDO
    epc = epc_sym
    !
    ! and finally the k+q electrons
    !
    DO nu = 1, nmodes
      DO ibnd = 1, nbnd
        !
        DO jbnd = 1, nbnd
          !
          w_1 = et2(jbnd, ikq2)
          g2 = 0.d0
          n  = 0
          !
          DO pbnd = 1, nbnd
            !
            w_2 = et2(pbnd, ikq2)
            !
            IF (ABS(w_2 - w_1) .LT. eps) THEN
              n = n + 1
              g2 = g2 + epc(ibnd, pbnd, nu) * epc(ibnd, pbnd, nu)
            ENDIF
            !
          ENDDO 
          !
          g2 = g2 / FLOAT(n)
          epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          !
        ENDDO
        !
      ENDDO
    ENDDO
    epc = epc_sym
    !
    WRITE(stdout, '(5x, a)') ' Electron-phonon vertex |g| (meV)'
    WRITE(stdout, '(/5x, "q coord.: ", 3f12.7)') xq
    WRITE(stdout, '(5x, "k coord.: ", 3f12.7)') kpoint
    WRITE(stdout, '(5x, a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
    WRITE(stdout, '(5x, a)') REPEAT('-', 78)
    !
    DO ibnd = 1, nbnd
      DO jbnd = 1, nbnd
        DO nu = 1, nmodes
          !
          IF (w2(nu) .GT. 0.d0) THEN
            w = DSQRT(w2(nu))
          ELSE
            w = DSQRT(-w2(nu))
          ENDIF
          !
          WRITE(stdout, '(3i9, 2f12.4, 1f20.10, 1e20.10)') ibnd, jbnd, nu, &
                rytoev * et2(ibnd, ikk2), rytoev * et2(jbnd, ikq2), &
                ryd2mev * w, ryd2mev * epc(ibnd, jbnd, nu)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    WRITE(stdout, '(5x, a/)') REPEAT('-', 78)
    !
  ENDIF
  !
  CALL stop_clock('elphsum2')
  !
  RETURN
  !
END SUBROUTINE elphsum2
