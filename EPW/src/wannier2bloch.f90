  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE wannier2bloch
  !----------------------------------------------------------------------
  !!
  !! Modules that contains all the routines that transforms quantities from Wannier
  !! space to Bloch space.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------
    SUBROUTINE hamwan2bloch(nbnd, nrr, cuf, eig, chw, cfac, is_mirror)
    !--------------------------------------------------------------------------
    !!
    !!  From the Hamiltonian in Wannier representation, find the corresponding
    !!  Hamiltonian in Bloch representation for a given k point
    !!
    !!  input  : number of bands nbnd
    !!           number of WS vectors, coordinates and degeneracy
    !!           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
    !!           kpoint coordinate xk(3)
    !!
    !!  output : rotation matrix cuf(nbnd, nbnd)
    !!           interpolated hamiltonian eigenvalues eig(nbnd)
    !!
    !!  2021: CL : Replace the random perturbation matrix with prime number matrix
    !!        in Lifting of degeneracies; control tag : lphase
    !!  2021: CL : Rotate the the largest element in eigenvector to real axis. (lrot)
    !!  2019: Weng Hong Sio and SP: Lifting of degeneracies. control tag: lrot
    !!        P_prime = U^dag P U where P is a random perturbation matrix
    !!        cuf = (eigvector of P_prime) * U
    !!        P_prime spans the degenenrate subspace.
    !!  2016: SP: optimization
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero, cone, zero, one, eps12, eps16
    USE input,         ONLY : lphase, lrot
    USE low_lvl,       ONLY : utility_zdotu, degen_sort, prime_number_matrix
    USE rigid,         ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    REAL(KIND = DP), INTENT(out) :: eig (nbnd)
    !! interpolated hamiltonian eigenvalues for this kpoint
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: chw(nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(KIND = DP), INTENT(out) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh
    LOGICAL, INTENT(in), OPTIONAL :: is_mirror
    !! .true. if k-point is a time-reversal invariant point
    !
    ! Local variables
    LOGICAL :: duplicates
    !! Returns if the bands contains degeneracices for that k-point.
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: info
    !! Infor for lapack DSYEV
    INTEGER :: list_dup(nbnd)
    !! List of degenerate eigenvalues
    INTEGER :: ndeg
    !! Number of degeneracies
    INTEGER :: starting
    !! Starting position
    INTEGER :: ending
    !! Ending position
    INTEGER :: length
    !! Size of the degenerate subspace
    INTEGER :: ig
    !! Degenerate group index
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd_max(1)
    !! Index of the maximum element
    INTEGER, ALLOCATABLE :: iwork(:)
    !! IWORK(1) returns the optimal LIWORK.
    INTEGER, ALLOCATABLE :: degen_group(:, :)
    !! Index of degenerate subspace
    REAL(KIND = DP) :: w(nbnd)
    !! Eigenvalues
    REAL(KIND = DP) :: norm_vec(nbnd)
    !! Real Hamiltonian matrix in Bloch basis for TRI k points
    REAL(KIND = DP), ALLOCATABLE :: P_prime(:, :)
    !! Perturbation matrix on the subspace
    REAL(KIND = DP), ALLOCATABLE :: wp(:)
    !! Perturbed eigenvalues on the degenerate subspace
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! RWORK(1) returns the optimal LRWORK.
    COMPLEX(KIND = DP) :: chf(nbnd, nbnd)
    !! Hamiltonian in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: cz(nbnd, nbnd)
    !! Eigenvectors from diag of Hamiltonian
    COMPLEX(KIND = DP), ALLOCATABLE :: cwork(:)
    !! Complex work variable
    COMPLEX(KIND = DP), ALLOCATABLE :: Uk(:, :)
    !! Rotation matrix on the full space
    COMPLEX(KIND = DP) :: cfac2(nrr)
    !! cfac for the time-reversial point if pointed
    COMPLEX(KIND = DP) :: fac_max
    !! Maximal normalized eigenvector.
    !
    !$omp master
    CALL start_clock('HamW2B')
    !$omp end master
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 31 of PRB 76, 165108 (2007)]
    !  H~(k')    = 1/ndegen(R) sum_R e^{ik'R     } H(R)
    !  H~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} H(R)
    !  Note H~(k') is H^(W)(k') in PRB 74, 195118 (2006) notations
    !
    !  H~(k')    is chf( nbnd, nbnd, 2*ik-1 )
    !  H~(k'+q') is chf( nbnd, nbnd, 2*ik   )
    !
    chf(:, :) = czero
    !
    ! Calculate the -k if this is a mirror point
    ! Do nothing if not specified
    cfac2 = cfac
    IF (PRESENT(is_mirror)) THEN
      IF (is_mirror) cfac2 = CONJG(cfac)
    ENDIF
    !
    CALL ZGEMV('n', nbnd**2, nrr, cone, chw, nbnd**2, cfac2, 1, cone, chf, 1)
    !
    !---------------------------------------------------------------------
    !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
    !---------------------------------------------------------------------
    !
    ! Hermitization
    DO jbnd = 1, nbnd
      DO ibnd = 1, jbnd
        chf(ibnd, jbnd) = (chf(ibnd, jbnd) + &
                             CONJG(chf(jbnd, ibnd))) * 0.5d0
      ENDDO
    ENDDO
    !
    ! Diagonalization routine
    CALL cdiagh2(nbnd, chf, nbnd, w, cz)
    !
    ! Find the degenerate eigenvalues w
    CALL degen_sort(w, SIZE(w), duplicates, list_dup)
    !
    ndeg = MAXVAL(list_dup)
    ALLOCATE(degen_group(2, ndeg), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating degen_group', 1)
    degen_group(:, :) = 0
    !
    ! degen_group contains the starting and ending position of each group
    ! degen_group(1,1) = starting position of group 1
    ! degen_group(2,1) = ending position of group 1
    ! degen_group(1,2) = starting position of group 2 ...
    DO ig = 1, ndeg
      degen_group(2, ig) = 0
      DO jbnd = 1, nbnd
        IF (list_dup(jbnd) == ig) THEN
          IF (jbnd == 1) THEN
            degen_group(1, ig) = jbnd
          ELSE
            IF (list_dup(jbnd) - list_dup(jbnd - 1) /= 0) degen_group(1, ig) = jbnd
          ENDIF
          degen_group(2, ig) = degen_group(2, ig) + 1
        ENDIF
      ENDDO
      degen_group(2, ig) = degen_group(1, ig) + degen_group(2, ig) -1
    ENDDO
    !
    DO ig = 1, ndeg
      starting = degen_group(1, ig)
      ending   = degen_group(2, ig)
      ! Size of the degenerate subspace
      length   = ending - starting + 1
      !
      ALLOCATE(rwork(length**2 + 2 * length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating rwork', 1)
      ALLOCATE(iwork(3 + 5 * length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating iwork', 1)
      ALLOCATE(cwork(length**2 + 2 * length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating cwork', 1)
      ALLOCATE(Uk(nbnd, length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating Uk', 1)
      ALLOCATE(P_prime(length, length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating P_prime', 1)
      ALLOCATE(wp(length), STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error allocating wp', 1)
      !
      Uk(:, :) = cz(:, starting:ending)
      ! Create a matrix filled with prime numbers
      !
      ! P_prime = MATMUL(TRANSPOSE(CONJG(Uk)), MATMUL(P, Uk))
      CALL prime_number_matrix(P_prime, length)
      !
      ! Diagonalization of P_prime
      CALL DSYEV('V', 'L', length, P_prime, length, wp, rwork, &
                length**2 + 2 * length, info)
      !
      ! On exiting P_prime is the eigenvector of the P_prime matrix and wp the eigenvector.
      !
      IF (lphase) cz(:, starting:ending) = MATMUL(Uk, P_prime)
      !
      DEALLOCATE(rwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating rwork', 1)
      DEALLOCATE(iwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating iwork', 1)
      DEALLOCATE(cwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating cwork', 1)
      DEALLOCATE(Uk, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating Uk', 1)
      DEALLOCATE(P_prime, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating P_prime', 1)
      DEALLOCATE(wp, STAT = ierr)
      IF (ierr /= 0) CALL errore('hamwan2bloch', 'Error deallocating wp', 1)
    ENDDO ! ig
    !
    ! Find the largest element and set it to pure real
    IF (lrot) THEN
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          norm_vec(jbnd) = ABS(cz(jbnd, ibnd))
        ENDDO
        ibnd_max(:) = MAXLOC(norm_vec(1:nbnd))
        fac_max = cz(ibnd_max(1), ibnd) / norm_vec(ibnd_max(1))
        cz(1:nbnd, ibnd) = cz(1:nbnd, ibnd) * CONJG(fac_max)
      ENDDO
    ENDIF
    !
    DO jbnd = 1, nbnd
      INNER : DO ibnd = 1, nbnd
        IF (ABS(cz(ibnd, jbnd)) > eps12) THEN
          cz(:, jbnd) = cz(:, jbnd) * CONJG(cz(ibnd, jbnd))
          cz(:, jbnd) = cz(:, jbnd) / SQRT(utility_zdotu(CONJG(cz(:, jbnd)), cz(:, jbnd)))
          EXIT INNER
        ENDIF
      ENDDO INNER
    ENDDO
    !
    ! Rotation matrix and Ham eigenvalues in Ryd [mind when comparing with wannier code (eV units)]
    ! U^\dagger is cuf(nbnd, nbnd)
    !
    cuf = CONJG(TRANSPOSE(cz))
    eig = w
    !
    ! Do the conjugate on eigenvector
    IF (PRESENT(is_mirror)) THEN
      IF(is_mirror) cuf = TRANSPOSE(cz)
    ENDIF
    !
    !$omp master
    CALL stop_clock('HamW2B')
    !$omp end master
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE hamwan2bloch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE hamwan2bloch_batch(nbnd, nrr, cuf, eig, chw, cfac, nkf, iq, chf, wb)
    !--------------------------------------------------------------------------
    !!
    !!  From the Hamiltonian in Wannier representation, find the corresponding
    !!  Hamiltonian in Bloch representation for a given k point, batched
    !!
    !!  input  : number of bands nbnd
    !!           number of WS vectors, coordinates and degeneracy
    !!           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
    !!           kpoint coordinate xk(3)
    !!
    !!  output : rotation matrix cuf(nbnd, nbnd)
    !!           interpolated hamiltonian eigenvalues eig(nbnd)
    !!
    !!  2021: CL : Replace the random perturbation matrix with prime number matrix
    !!        in Lifting of degeneracies; control tag : lphase
    !!  2021: CL : Rotate the the largest element in eigenvector to real axis. (lrot)
    !!  2019: Weng Hong Sio and SP: Lifting of degeneracies. control tag: lrot
    !!        P_prime = U^dag P U where P is a random perturbation matrix
    !!        cuf = (eigvector of P_prime) * U
    !!        P_prime spans the degenenrate subspace.
    !!  2016: SP: optimization
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero, cone, zero, one, eps12, eps16
    USE input,         ONLY : lphase, lrot, fsthick
    USE low_lvl,       ONLY : utility_zdotu, degen_sort, prime_number_matrix
    USE ep_constants,  ONLY : cone, czero, ci, twopi
    USE pwcom,         ONLY : ef
#if defined (__CUDA)
    USE global_var,    ONLY : cublas_h, cusolverdn_h
    USE cublas
    USE cusolverdn
#endif
    USE iso_c_binding
#if defined (__ONEMKL)
    USE epw_onemkl
#endif
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    REAL(KIND = DP), INTENT(out) :: eig (nbnd, 2*nkf)
    !! interpolated hamiltonian eigenvalues for this kpoint
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, 2*nkf)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: chw(nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(KIND = DP), INTENT(out) :: cuf(nbnd, nbnd, 2*nkf)
    !! Rotation matrix U^\dagger, fine mesh
    INTEGER, INTENT(in) :: nkf
    !! Number of k points in the fine mesh (distributed over pools)
    INTEGER, INTENT(in) :: iq
    !! Index of the q point
    COMPLEX(KIND = DP), INTENT(inout) :: chf(nbnd, nbnd, 2*nkf)
    !! Electron Hamiltonian in fine mesh
    REAL(KIND = DP), INTENT(inout) :: wb(nbnd, 2*nkf)
    !! Electron energies in fine mesh
    !
    ! Local variables
    LOGICAL :: duplicates
    !! Returns if the bands contains degeneracices for that k-point.
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: info
    !! Infor for lapack ZHPEVX
    INTEGER :: list_dup(nbnd)
    !! List of degenerate eigenvalues
    INTEGER :: ndeg
    !! Number of degeneracies
    INTEGER :: starting
    !! Starting position
    INTEGER :: ending
    !! Ending position
    INTEGER :: length
    !! Size of the degenerate subspace
    INTEGER :: ig
    !! Degenerate group index
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd_max(1)
    !! Index of the maximum element
    INTEGER, ALLOCATABLE :: iwork(:)
    !! IWORK(1) returns the optimal LIWORK.
    INTEGER, ALLOCATABLE :: degen_group(:, :)
    !! Index of degenerate subspace
    REAL(KIND = DP) :: w(nbnd)
    !! Eigenvalues
    REAL(KIND = DP) :: norm_vec(nbnd)
    !! Real Hamiltonian matrix in Bloch basis for TRI k points
    REAL(KIND = DP), ALLOCATABLE :: P_prime(:, :)
    !! Perturbation matrix on the subspace
    REAL(KIND = DP), ALLOCATABLE :: wp(:)
    !! Perturbed eigenvalues on the degenerate subspace
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! RWORK(1) returns the optimal LRWORK.
    COMPLEX(KIND = DP) :: cz(nbnd, nbnd)
    !! Eigenvectors from diag of Hamiltonian
    COMPLEX(KIND = DP), ALLOCATABLE :: cwork(:)
    !! Complex work variable
    COMPLEX(KIND = DP), ALLOCATABLE :: Uk(:, :)
    !! Rotation matrix on the full space
    COMPLEX(KIND = DP) :: cfac2(nrr)
    !! cfac for the time-reversial point if pointed
    COMPLEX(KIND = DP) :: fac_max
    !! Maximal normalized eigenvector.
    INTEGER :: ik
    !! Counter on k point index
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ikq
    !! Paired counter so that q is adjacent to its k
    LOGICAL :: mirror_k
    !! Mirror symmetry for k-points
    LOGICAL :: mirror_kpq
    !! Mirror symmetry for (k+q)-points
    INTEGER :: lwork
    !! Size of the work array for diagonalization
    COMPLEX(KIND = DP), ALLOCATABLE :: work(:)
    !! Work array for diagonalization
    INTEGER, ALLOCATABLE :: infob(:)
    !! Info array for diagonalization
#if defined (__CUDA)
    INTERFACE
      INTEGER(C_INT) FUNCTION CUSOLVERDN_CREATESYEVJINFO(info) &
        BIND(C, NAME='cusolverDnCreateSyevjInfo')
        IMPORT
        TYPE(C_PTR) :: info
      END FUNCTION
      !! create and initialize the syevj info object
      INTEGER(C_INT) FUNCTION CUSOLVERDN_DESTROYSYEVJINFO(info) &
        BIND(C, NAME='cusolverDnDestroySyevjInfo')
        IMPORT
        TYPE(C_PTR), VALUE :: info
      END FUNCTION
      !! destroy the syevj info object
      INTEGER(C_INT) FUNCTION CUSOLVERDN_ZHEEVJBATCHED_BUFFERSIZE( &
        handle, jobz, uplo, n, A, lda, W, lwork, params, batchSize) &
        BIND(C, NAME='cusolverDnZheevjBatched_bufferSize')
        IMPORT
        TYPE(CUSOLVERDNHANDLE), VALUE :: handle
        INTEGER(C_INT), VALUE :: jobz, uplo, n
        COMPLEX(C_DOUBLE), DEVICE :: A(*)
        REAL(C_DOUBLE), DEVICE :: W(*)
        INTEGER(C_INT), VALUE :: lda
        INTEGER(C_INT) :: lwork
        TYPE(C_PTR), VALUE :: params
        INTEGER(C_INT), VALUE :: batchSize
      END FUNCTION
      !! calculate the size of workspace for the syevj
      INTEGER(C_INT) FUNCTION CUSOLVERDN_ZHEEVJBATCHED( &
        handle, jobz, uplo, n, A, lda, W, work, lwork, info, params, batchSize) &
        BIND(C, NAME='cusolverDnZheevjBatched')
        IMPORT
        TYPE(CUSOLVERDNHANDLE), VALUE :: handle
        INTEGER(C_INT), VALUE :: jobz, uplo, n
        COMPLEX(C_DOUBLE), DEVICE :: A(*)
        REAL(C_DOUBLE), DEVICE :: W(*)
        COMPLEX(C_DOUBLE), DEVICE :: work(*)
        INTEGER(C_INT), DEVICE :: info(*)
        INTEGER(C_INT), VALUE :: lda
        INTEGER(C_INT), VALUE :: lwork
        TYPE(C_PTR), VALUE :: params
        INTEGER(C_INT), VALUE :: batchSize
      END FUNCTION
      !! compute the eigenvalues and eigenvectors of a batch of hermitian matrices
    END INTERFACE
    !! custom Fortran interface for CUSOLVERDN subroutines
    TYPE(C_PTR) :: syevj_info
    !! syevj info object
#endif
    !
    !$omp master
    CALL start_clock('HamW2B')
    !$omp end master
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 31 of PRB 76, 165108 (2007)]
    !  H~(k')    = 1/ndegen(R) sum_R e^{ik'R     } H(R)
    !  H~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} H(R)
    !  Note H~(k') is H^(W)(k') in PRB 74, 195118 (2006) notations
    !
    !  H~(k')    is chf( nbnd, nbnd, 2*ik-1 )
    !  H~(k'+q') is chf( nbnd, nbnd, 2*ik   )
    !
    ! Calculate the -k if this is a mirror point
    ! Do nothing if not specified
    !
#if defined (__CUDA)
    !$acc enter data copyin(chw)
    !$acc host_data use_device(chw,cfac,chf)
    ! ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 'n', 'n', nbnd**2, 1, nrr, cone, &
    !   chw, nbnd**2, 0, cfac, nrr, nrr, czero, chf, nbnd**2, nbnd**2, 2*nkf)
    CALL CUBLASZGEMM('n', 'n', nbnd**2, 2*nkf, nrr, cone, chw, nbnd**2, cfac, nrr, czero, chf, nbnd**2)
    !$acc end host_data
    !$acc exit data delete(chw)
#elif defined(__ONEMKL)
    !$omp target enter data map(to:chw)
    !$omp dispatch
    CALL ZGEMM('n', 'n', nbnd**2, 2*nkf, nrr, cone, chw, nbnd**2, cfac, nrr, czero, chf, nbnd**2)
    !$omp target exit data map(delete:chw)
    !$omp target update from(chf)
#else
    CALL ZGEMM('n', 'n', nbnd**2, 2*nkf, nrr, cone, chw, nbnd**2, cfac, nrr, czero, chf, nbnd**2)
#endif
    !
#if defined (__CUDA)
    !$acc host_data use_device(chf,wb)
    ierr = CUSOLVERDN_CREATESYEVJINFO(syevj_info)
    ierr = CUSOLVERDN_ZHEEVJBATCHED_BUFFERSIZE(cusolverdn_h, 1, 0, nbnd, &
      chf, nbnd, wb, lwork, syevj_info, 2*nkf)
    !$acc end host_data
    !
    ALLOCATE(work(lwork), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating work', 1)
    ALLOCATE(infob(2*nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating infob', 1)
    !$acc enter data create(work,infob)
    !
    !$acc host_data use_device(chf,wb,work,infob)
    ierr = CUSOLVERDN_ZHEEVJBATCHED(cusolverdn_h, 1, 0, nbnd, &
      chf, nbnd, wb, work, lwork, infob, syevj_info, 2*nkf)
    ierr = CUSOLVERDN_DESTROYSYEVJINFO(syevj_info)
    !$acc end host_data

    !$acc update self(chf,wb)
    !$acc exit data delete(work,infob)
    DEALLOCATE(work, infob)
#else
    !
    !$omp parallel private(cz,w,rwork,iwork,cwork,info,ierr,ikk,ikq)
    ALLOCATE(rwork(1 + 5 * nbnd + 2 * (nbnd**2)), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating rwork', 1)
    ALLOCATE(iwork(3 + 5 * nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating iwork', 1)
    ALLOCATE(cwork(nbnd**2 + 2 * nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating cwork', 1)
    !
    !$omp do
    DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! NOTE: we don't need to slice chf here
        cz(:, :) = chf(:, :, ikk)
        CALL ZHEEVD('V', 'L', nbnd, cz, nbnd, w, cwork, 2 * nbnd + nbnd**2, &
                rwork, 1 + 5 * nbnd + 2 * (nbnd**2), iwork, 3 + 5 * nbnd, info)
        !
        chf(:, :, ikk) = cz
        wb(:, ikk) = w
        !
        cz(:, :) = chf(:, :, ikq)
        CALL ZHEEVD('V', 'L', nbnd, cz, nbnd, w, cwork, 2 * nbnd + nbnd**2, &
                rwork, 1 + 5 * nbnd + 2 * (nbnd**2), iwork, 3 + 5 * nbnd, info)
        !
        chf(:, :, ikq) = cz
        wb(:, ikq) = w
    ENDDO
    !$omp end do
    !
    DEALLOCATE(rwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating rwork', 1)
    DEALLOCATE(iwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating iwork', 1)
    DEALLOCATE(cwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating cwork', 1)
    !$omp end parallel
#endif
    !
    !$omp parallel do schedule(dynamic) &
    !$omp private(ikk,ikq,cz,w,duplicates,list_dup,ndeg,degen_group,ierr) &
    !$omp private(ig,jbnd,starting,ending,length,rwork,iwork,cwork,Uk,P_prime) &
    !$omp private(wp,info,ibnd,norm_vec,ibnd_max,fac_max,mirror_k,mirror_kpq)
    DO ik = 1, nkf
      !
      ! xkf is assumed to be in crys coord
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      IF ((MINVAL(ABS(wb(:, ikk) - ef)) < fsthick) .AND. (MINVAL(ABS(wb(:, ikq) - ef)) < fsthick)) THEN
        !
        !---------------------------------------------------------------------
        !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
        !---------------------------------------------------------------------
        !
        ! Hermitization
        ! chf(:, :, ikk) = 0.5d0 * (chf(:, :, ikk) + TRANSPOSE(CONJG(chf(:, :, ikk))))
        !
        ! Diagonalized Hamiltonian
        cz(:, :) = chf(:, :, ikk)
        w(:) = wb(:, ikk)
        !
        ! Find the degenerate eigenvalues w
        CALL degen_sort(w, SIZE(w), duplicates, list_dup)
        !
        ndeg = MAXVAL(list_dup)
        ALLOCATE(degen_group(2, ndeg), STAT = ierr)
        IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating degen_group for k', 1)
        degen_group(:, :) = 0
        !
        ! degen_group contains the starting and ending position of each group
        ! degen_group(1,1) = starting position of group 1
        ! degen_group(2,1) = ending position of group 1
        ! degen_group(1,2) = starting position of group 2 ...
        DO ig = 1, ndeg
          degen_group(2, ig) = 0
          DO jbnd = 1, nbnd
            IF (list_dup(jbnd) == ig) THEN
              IF (jbnd == 1) THEN
                degen_group(1, ig) = jbnd
              ELSE
                IF (list_dup(jbnd) - list_dup(jbnd - 1) /= 0) degen_group(1, ig) = jbnd
              ENDIF
              degen_group(2, ig) = degen_group(2, ig) + 1
            ENDIF
          ENDDO
          degen_group(2, ig) = degen_group(1, ig) + degen_group(2, ig) -1
        ENDDO
        !
        DO ig = 1, ndeg
          starting = degen_group(1, ig)
          ending   = degen_group(2, ig)
          ! Size of the degenerate subspace
          length   = ending - starting + 1
          !
          ALLOCATE(rwork(length**2 + 2 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating rwork', 1)
          ALLOCATE(iwork(3 + 5 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating iwork', 1)
          ALLOCATE(cwork(length**2 + 2 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating cwork', 1)
          ALLOCATE(Uk(nbnd, length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating Uk', 1)
          ALLOCATE(P_prime(length, length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating P_prime', 1)
          ALLOCATE(wp(length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating wp', 1)
          !
          Uk(:, :) = cz(:, starting:ending)
          ! Create a matrix filled with prime numbers
          !
          ! P_prime = MATMUL(TRANSPOSE(CONJG(Uk)), MATMUL(P, Uk))
          CALL prime_number_matrix(P_prime, length)
          !
          ! Diagonalization of P_prime
          CALL DSYEV('V', 'L', length, P_prime, length, wp, rwork, &
                    length**2 + 2 * length, info)
          !
          ! On exiting P_prime is the eigenvector of the P_prime matrix and wp the eigenvector.
          !
          IF (lphase) cz(:, starting:ending) = MATMUL(Uk, P_prime)
          !
          DEALLOCATE(rwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating rwork', 1)
          DEALLOCATE(iwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating iwork', 1)
          DEALLOCATE(cwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating cwork', 1)
          DEALLOCATE(Uk, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating Uk', 1)
          DEALLOCATE(P_prime, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating P_prime', 1)
          DEALLOCATE(wp, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating wp', 1)
        ENDDO ! ig
        DEALLOCATE(degen_group)
        !
        ! Find the largest element and set it to pure real
        IF (lrot) THEN
          DO ibnd = 1, nbnd
            DO jbnd = 1, nbnd
              norm_vec(jbnd) = ABS(cz(jbnd, ibnd))
            ENDDO
            ibnd_max(:) = MAXLOC(norm_vec)
            fac_max = cz(ibnd_max(1), ibnd) / norm_vec(ibnd_max(1))
            cz(1:nbnd, ibnd) = cz(1:nbnd, ibnd) * CONJG(fac_max)
          ENDDO
        ENDIF
        !
        DO jbnd = 1, nbnd
          INNER : DO ibnd = 1, nbnd
            IF (ABS(cz(ibnd, jbnd)) > eps12) THEN
              cz(:, jbnd) = cz(:, jbnd) * CONJG(cz(ibnd, jbnd))
              cz(:, jbnd) = cz(:, jbnd) / SQRT(utility_zdotu(CONJG(cz(:, jbnd)), cz(:, jbnd)))
              EXIT INNER
            ENDIF
          ENDDO INNER
        ENDDO
        !
        ! Rotation matrix and Ham eigenvalues in Ryd [mind when comparing with wannier code (eV units)]
        ! U^\dagger is cuf(nbnd, nbnd)
        !
        cuf(:, :, ikk) = CONJG(TRANSPOSE(cz))
        eig(:, ikk) = w
        !
        ! ! Do the conjugate on eigenvector
        ! IF (plrn .AND. time_rev_U_plrn) THEN
        !   mirror_k = is_mirror_k(ik)
        ! ELSE
        !   mirror_k = .FALSE.
        ! ENDIF
        ! IF(mirror_k) cuf(:, :, ikk) = TRANSPOSE(cz)
        !
        !---------------------------------------------------------------------
        !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
        !---------------------------------------------------------------------
        !
        ! Hermitization
        ! chf(:, :, ikq) = 0.5d0 * (chf(:, :, ikq) + TRANSPOSE(CONJG(chf(:, :, ikq))))
        !
        ! Diagonalized Hamiltonian
        cz(:, :) = chf(:, :, ikq)
        w(:) = wb(:, ikq)
        !
        ! Find the degenerate eigenvalues w
        CALL degen_sort(w, SIZE(w), duplicates, list_dup)
        !
        ndeg = MAXVAL(list_dup)
        ALLOCATE(degen_group(2, ndeg), STAT = ierr)
        IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating degen_group for k+q', 1)
        degen_group(:, :) = 0
        !
        ! degen_group contains the starting and ending position of each group
        ! degen_group(1,1) = starting position of group 1
        ! degen_group(2,1) = ending position of group 1
        ! degen_group(1,2) = starting position of group 2 ...
        DO ig = 1, ndeg
          degen_group(2, ig) = 0
          DO jbnd = 1, nbnd
            IF (list_dup(jbnd) == ig) THEN
              IF (jbnd == 1) THEN
                degen_group(1, ig) = jbnd
              ELSE
                IF (list_dup(jbnd) - list_dup(jbnd - 1) /= 0) degen_group(1, ig) = jbnd
              ENDIF
              degen_group(2, ig) = degen_group(2, ig) + 1
            ENDIF
          ENDDO
          degen_group(2, ig) = degen_group(1, ig) + degen_group(2, ig) -1
        ENDDO
        !
        DO ig = 1, ndeg
          starting = degen_group(1, ig)
          ending   = degen_group(2, ig)
          ! Size of the degenerate subspace
          length   = ending - starting + 1
          !
          ALLOCATE(rwork(length**2 + 2 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating rwork', 1)
          ALLOCATE(iwork(3 + 5 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating iwork', 1)
          ALLOCATE(cwork(length**2 + 2 * length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating cwork', 1)
          ALLOCATE(Uk(nbnd, length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating Uk', 1)
          ALLOCATE(P_prime(length, length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating P_prime', 1)
          ALLOCATE(wp(length), STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error allocating wp', 1)
          !
          Uk(:, :) = cz(:, starting:ending)
          ! Create a matrix filled with prime numbers
          !
          ! P_prime = MATMUL(TRANSPOSE(CONJG(Uk)), MATMUL(P, Uk))
          CALL prime_number_matrix(P_prime, length)
          !
          ! Diagonalization of P_prime
          CALL DSYEV('V', 'L', length, P_prime, length, wp, rwork, &
                    length**2 + 2 * length, info)
          !
          ! On exiting P_prime is the eigenvector of the P_prime matrix and wp the eigenvector.
          !
          IF (lphase) cz(:, starting:ending) = MATMUL(Uk, P_prime)
          !
          DEALLOCATE(rwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating rwork', 1)
          DEALLOCATE(iwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating iwork', 1)
          DEALLOCATE(cwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating cwork', 1)
          DEALLOCATE(Uk, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating Uk', 1)
          DEALLOCATE(P_prime, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating P_prime', 1)
          DEALLOCATE(wp, STAT = ierr)
          IF (ierr /= 0) CALL errore('hamwan2bloch_batch', 'Error deallocating wp', 1)
        ENDDO ! ig
        !
        DEALLOCATE(degen_group)
        !
        ! Find the largest element and set it to pure real
        IF (lrot) THEN
          DO ibnd = 1, nbnd
            DO jbnd = 1, nbnd
              norm_vec(jbnd) = ABS(cz(jbnd, ibnd))
            ENDDO
            ibnd_max(:) = MAXLOC(norm_vec)
            fac_max = cz(ibnd_max(1), ibnd) / norm_vec(ibnd_max(1))
            cz(1:nbnd, ibnd) = cz(1:nbnd, ibnd) * CONJG(fac_max)
          ENDDO
        ENDIF
        !
        DO jbnd = 1, nbnd
          INNER2 : DO ibnd = 1, nbnd
            IF (ABS(cz(ibnd, jbnd)) > eps12) THEN
              cz(:, jbnd) = cz(:, jbnd) * CONJG(cz(ibnd, jbnd))
              cz(:, jbnd) = cz(:, jbnd) / SQRT(utility_zdotu(CONJG(cz(:, jbnd)), cz(:, jbnd)))
              EXIT INNER2
            ENDIF
          ENDDO INNER2
        ENDDO
        !
        ! Rotation matrix and Ham eigenvalues in Ryd [mind when comparing with wannier code (eV units)]
        ! U^\dagger is cuf(nbnd, nbnd)
        !
        cuf(:, :, ikq) = CONJG(TRANSPOSE(cz))
        eig(:, ikq) = w
        !
        ! ! Do the conjugate on eigenvector
        ! IF (plrn .AND. time_rev_U_plrn) THEN
        !   ! kpg_map return global index, thus xkf_all is needed
        !   mirror_kpq = is_mirror_q(ikq_all(ik, iq))
        ! ELSE
        !   mirror_kpq = .FALSE.
        ! ENDIF
        ! IF(mirror_kpq) cuf(:, :, ikq) = TRANSPOSE(cz)
      ELSE
        eig(:, ikk) = wb(:, ikk)
        eig(:, ikq) = wb(:, ikq)
      ENDIF
    ENDDO
    !$omp end parallel do
    !
    !$acc update device(cuf)
    !$omp target update to(cuf)
    !
    !$omp master
    CALL stop_clock('HamW2B')
    !$omp end master
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE hamwan2bloch_batch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, cuf, eig, is_mirror)
    !--------------------------------------------------------------------------
    !!
    !! From the Hamiltonian in Wannier representation, find the corresponding
    !! Hamiltonian in Bloch representation for a given k point
    !!
    !! This SUBROUTINE is identical to hamwan2bloch.f90, except
    !! that here rdw is a real array, not a complex one. This is
    !! required to obtain proper phonon dispersion interpolation
    !! and corresponds to the reality of the interatomic force constants
    !!
    !!  2021: CL : Lifting of degeneracies using random perturbation matrix
    !!             with prime number matrix control tag : lphase
    !!             Rotate the the largest element in eigenvector to real axis. (lrot)
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE ions_base,     ONLY : amass, tau, nat, ityp
    USE global_var,    ONLY : rdw, epsi, zstar, qrpl, L
    USE input,         ONLY : lpolar, lphase, lrot, use_ws, nqc1, nqc2, nqc3
    USE ep_constants,  ONLY : twopi, ci, czero, zero, one, eps12
    USE rigid,         ONLY : cdiagh2
    USE low_lvl,       ONLY : utility_zdotu, prime_number_matrix, degen_sort
    USE longrange,     ONLY : rgd_blk
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of modes (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr_q
    !! number of WS points
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! coordinates of phononic WS points
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, nat, nat)
    !! degeneracy of WS points
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! kpoint coordinates for the interpolation
    REAL(KIND = DP), INTENT(out) :: eig(nmodes)
    !! interpolated dynamical matrix eigenvalues for this kpoint
    COMPLEX(KIND = DP), INTENT(out) :: cuf(nmodes, nmodes)
    !! Rotation matrix, fine mesh
    LOGICAL, INTENT(in), OPTIONAL :: is_mirror
    !! .true. if q-point is a the mirror point of some original point
    !
    ! Local variables
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: jmode
    !! Counter on modes
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: ierr
    !! Error status
    INTEGER :: list_dup(nmodes)
    !! List of degenerate eigenvalues
    INTEGER :: ndeg
    !! Number of degeneracies
    LOGICAL :: duplicates
    !! Returns if the bands contains degeneracices for that k-point.
    INTEGER :: ig
    !! Counter on real-space index
    INTEGER :: info
    !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
    INTEGER :: starting
    !! Starting position
    INTEGER :: ending
    !! Ending position
    INTEGER :: length
    !! Size of the degenerate subspace
    INTEGER :: imode_max(1)
    !! Size of the degenerate subspace
    INTEGER, ALLOCATABLE :: degen_group(:, :)
    !! Index of degenerate subspace
    INTEGER, ALLOCATABLE :: iwork(:)
    !! IWORK(1) returns the optimal LIWORK.
    REAL(KIND = DP) :: w(nmodes)
    !! Eigenvalues
    REAL(KIND = DP) :: xq(3)
    !! Coordinates q-point
    REAL(KIND = DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    REAL(KIND = DP) :: massfac
    !! inverse square root of masses
    REAL(KIND = DP) :: norm_vec(nmodes)
    !! Real Dynamical matrix in Bloch basis for TRI q points
    REAL(KIND = DP), ALLOCATABLE :: wp(:)
    !! Perturbed eigenvalues on the degenerate subspace
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: P_prime(:, :)
    !! Perturbation matrix on the subspace
    COMPLEX(KIND = DP) :: cz(nmodes, nmodes)
    !! Eigenvectors
    COMPLEX(KIND = DP) :: chf(nmodes, nmodes)
    ! Dynamical matrix in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: cfac
    !! Complex prefactor for Fourier transform.
    COMPLEX(KIND = DP) :: fac_max(1)
    !! max factor
    COMPLEX(KIND = DP), ALLOCATABLE :: Uk(:, :)
    !! Rotation matrix on the full space
    !
    CALL start_clock ('DynW2B')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 32 of PRB 76, 165108 (2007)]
    !  D~(k')    = 1/ndegen(R) sum_R e^{ik'R     } D(R)
    !  D~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} D(R)
    !
    !  D~(k')    is chf ( nmodes, nmodes, 2*ik-1 )
    !  D~(k'+q') is chf ( nmodes, nmodes, 2*ik   )
    !
    xq = xxq
    chf(:, :) = czero
    !
    IF (use_ws) THEN
      DO ir = 1, nrr_q
        rdotk = twopi * DOT_PRODUCT(xq, DBLE(irvec_q(:, ir)))
        DO na = 1, nat
          DO nb = 1, nat
            IF (ndegen_q(ir, na, nb) > 0) THEN
              cfac = EXP(ci * rdotk) / DBLE(ndegen_q(ir, na, nb))
              ! To map atom coordinate to mode basis.
              chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
              chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) &
              + cfac * rdw(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, ir)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE ! use_ws
      DO ir = 1, nrr_q
        rdotk = twopi * DOT_PRODUCT(xq, DBLE(irvec_q(:, ir)))
        cfac = EXP(ci * rdotk) / DBLE(ndegen_q(ir, 1, 1))
        chf = chf + cfac * rdw(:, :, ir)
      ENDDO
    ENDIF
    !
    ! bring xq in cart. coordinates (needed for rgd_blk call)
    CALL cryst_to_cart(1, xq, bg, 1)
    !
    !  add the long-range term to D(q)
    IF (lpolar .OR. qrpl) THEN
      ! xq has to be in 2pi/a
      CALL rgd_blk(L, nqc1, nqc2, nqc3, nat, chf, xq, tau, epsi, zstar, +1.d0)
      !
    ENDIF
    !
    !  divide by the square root of masses
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)) )
        !
        chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
        chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) * massfac
        !
      ENDDO
    ENDDO
    !
    ! bring xq back to crystal coordinates
    CALL cryst_to_cart(1, xq, at, -1)
    !
    !---------------------------------------------------------------------
    !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
    !---------------------------------------------------------------------
    !
    ! Here, we hermitianize only the upper triangular part for cdiagh2
    !
    DO jmode = 1, nmodes
      DO imode = 1, jmode
        chf(imode, jmode) = (chf(imode, jmode) + &
                             CONJG(chf(jmode, imode))) * 0.5d0
      ENDDO
    ENDDO
    !
    CALL cdiagh2(nmodes, chf, nmodes, w, cz)
    !
    ! Find the degenerate eigenvalues w
    CALL degen_sort(w, SIZE(w), duplicates, list_dup)
    !
    ndeg = MAXVAL(list_dup)
    ALLOCATE(degen_group(2, ndeg), STAT = ierr)
    IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating degen_group', 1)
    degen_group(:, :) = 0
    !
    ! degen_group contains the starting and ending position of each group
    ! degen_group(1,1) = starting position of group 1
    ! degen_group(2,1) = ending position of group 1
    ! degen_group(1,2) = starting position of group 2 ...
    DO ig = 1, ndeg
      degen_group(2, ig) = 0
      DO jmode = 1, nmodes
        IF (list_dup(jmode) == ig) THEN
          IF (jmode == 1) THEN
            degen_group(1, ig) = jmode
          ELSE
            IF (list_dup(jmode) - list_dup(jmode - 1) /= 0) degen_group(1, ig) = jmode
          ENDIF
          degen_group(2, ig) = degen_group(2, ig) + 1
        ENDIF
      ENDDO
      degen_group(2, ig) = degen_group(1, ig) + degen_group(2, ig) -1
    ENDDO
    !
    ! Generate a pertubation matrix of size (nbnd x nbnd) made of random number
    !
    DO ig = 1, ndeg
      starting = degen_group(1, ig)
      ending   = degen_group(2, ig)
      ! Size of the degenerate subspace
      length   = ending - starting + 1
      !
      ALLOCATE(rwork(length**2 + 2 * length), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating rwork', 1)
      ALLOCATE(iwork(3 + 5 * length), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating iwork', 1)
      ALLOCATE(Uk(nmodes, length), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating Uk', 1)
      ALLOCATE(P_prime(length, length), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating P_prime', 1)
      ALLOCATE(wp(length), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error allocating wp', 1)
      !
      Uk(:, 1:length) = cz(:, starting:ending)
      CALL prime_number_matrix(P_prime, length)
      ! Diagonalization of P_prime
      CALL DSYEV('V', 'L', length, P_prime, length, wp, rwork, &
                length**2 + 2 * length, info)
      ! On exiting P_prime is the eigenvector of the P_prime matrix and wp the eigenvector.
      !
      IF(lphase) cz(:, starting:ending) = MATMUL(Uk, P_prime)
      !
      DEALLOCATE(rwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error deallocating rwork', 1)
      DEALLOCATE(iwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error deallocating iwork', 1)
      DEALLOCATE(Uk, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error deallocating Uk', 1)
      DEALLOCATE(P_prime, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error deallocating P_prime', 1)
      DEALLOCATE(wp, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynwan2bloch', 'Error deallocating wp', 1)
    ENDDO ! ig
    !
    ! clean noise
    DO jmode = 1, nmodes
      DO imode = 1, nmodes
        IF (ABS(cz(imode, jmode)) < eps12) cz(imode, jmode) = czero
      ENDDO
    ENDDO
    !
    ! Find the largest element and set it to pure real
    IF(lrot) THEN
       DO imode = 1, nmodes
         DO jmode = 1, nmodes
            norm_vec(jmode) = ABS(cz(jmode, imode))
         END DO
         imode_max(:) = MAXLOC(norm_vec(1:nmodes))
         fac_max(1) = cz(imode_max(1), imode)/norm_vec(imode_max(1))
         cz(1:nmodes, imode) = cz(1:nmodes, imode) * CONJG(fac_max(1))
       ENDDO
     ENDIF
    !
    ! DS - Impose phase
    IF (lphase) THEN
      DO jmode = 1,nmodes
        INNER : DO imode = 1, nmodes
          IF (ABS(cz(imode, jmode)) > eps12) THEN
            cz(:, jmode) = cz(:, jmode) * CONJG(cz(imode,jmode))
            cz(:, jmode) = cz(:, jmode) / SQRT(utility_zdotu(CONJG(cz(:, jmode)),cz(:, jmode)))
            EXIT INNER
          ENDIF
        ENDDO INNER
      ENDDO
    ENDIF
    !
    ! cuf(nmodes,nmodes) is rotation matrix (eigenmodes e_k)
    !
    cuf = cz
    eig = w
    !
    IF(PRESENT(is_mirror)) THEN
      IF(is_mirror) cuf = CONJG(cz)
    ENDIF
    !
    CALL stop_clock('DynW2B')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE dynwan2bloch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dynifc2blochf(nmodes, rws, nrws, xxq, cuf, eig, is_mirror)
    !--------------------------------------------------------------------------
    !!
    !!  From the IFCs in the format of q2r, find the corresponding
    !!  dynamical matrix for a given q point (as in matdyn.x) on the fine grid
    !!  2021: CL : Lifting of degeneracies using random perturbation matrix
    !!             with prime number matrix control tag : lphase
    !!             Rotate the the largest element in eigenvector to real axis. (lrot)
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE ions_base,     ONLY : amass, tau, nat, ityp
    USE global_var,    ONLY : ifc, epsi, zstar, wscache, qrpl, L
    USE input,         ONLY : lpolar, nqc1, nqc2, nqc3, lphase, lrot, fd
    USE io_global,     ONLY : stdout
    USE longrange,     ONLY : rgd_blk
    USE low_lvl,       ONLY : utility_zdotu, degen_sort
    USE ep_constants,  ONLY : twopi, czero, zero, one, eps8, eps12
    USE rigid,         ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of modes
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Real space Wigner-Seitz vector
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! qpoint coordinates for the interpolation
    REAL(KIND = DP), INTENT(out) :: eig(nmodes)
    !! interpolated phonon eigenvalues for this qpoint
    COMPLEX(KIND = DP), INTENT(out) :: cuf(nmodes, nmodes)
    !! Rotation matrix, fine mesh
    LOGICAL, INTENT(in), OPTIONAL :: is_mirror
    !! .true. if q-point is a time-reversal point
    !
    ! Local variables
    LOGICAL, SAVE :: first = .TRUE.
    !! First time
    INTEGER :: n1, n2, n3
    !! Q-point grid dimensions
    INTEGER :: m1, m2, m3
    !! Mod of the grid dim
    INTEGER :: i
    !! Cartesian direction
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: jmode
    !! Counter on modes
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: jpol
    !! Counter on polarizations
    INTEGER :: imode_max(1)
    !! index of the max element
    REAL(KIND = DP) :: xq(3)
    !! Coordinates q-point
    REAL(KIND = DP) :: massfac
    !! inverse square root of masses
    REAL(KIND = DP) :: w(nmodes)
    !! Eigenvalues
    REAL(KIND = DP) :: total_weight
    !! Sum of the weigths
    REAL(KIND = DP) :: weight
    !! WS weights
    REAL(KIND = DP) :: arg
    !! 2 * pi * r
    REAL(KIND = DP) :: r(3)
    !! Real-space vector
    REAL(KIND = DP) :: r_ws(3)
    !! Real space vector including fractional translation
    REAL(KIND = DP) :: norm_vec(nmodes)
    !! Vector for one eigenmode
    REAL(KIND = DP), EXTERNAL :: wsweight
    !! Wigner-Seitz weights
    COMPLEX(KIND = DP) :: cz(nmodes, nmodes)
    !! Eigenvectors
    COMPLEX(KIND = DP) :: chf(nmodes, nmodes)
    !! Dynamical matrix in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: dyn(3, 3, nat, nat)
    !! Dynamical matrix
    COMPLEX(KIND = DP) :: fac_max(1)
    !! value of the max element
    !
    CALL start_clock('DynW2B')
    !
    xq = xxq
    ! bring xq in cart. coordinates
    CALL cryst_to_cart(1, xq, bg, 1)
    !
    IF (first) THEN
      first = .FALSE.
      DO na = 1, nat
        DO nb = 1, nat
          DO n1 = -2 * nqc1, 2 * nqc1
            DO n2 = -2 * nqc2, 2 * nqc2
              DO n3 = -2 * nqc3, 2 * nqc3
                DO i = 1, 3
                  r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
                  ! SM: IFC's  are from finite displacement calculation (see PH/matdyn.f90)
                  IF (fd) THEN
                    r_ws(i) = r(i)
                  ELSE
                    r_ws(i) = r(i) + tau(i, na) - tau(i, nb)
                  ENDIF
                END DO
                wscache(n3, n2, n1, nb, na) = wsweight(r_ws, rws, nrws)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 32 of PRB 76, 165108 (2007)]
    !  D~(k')    = 1/ndegen(R) sum_R e^{ik'R     } D(R)
    !  D~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} D(R)
    !
    !  D~(k')    is chf ( nmodes, nmodes, 2*ik-1 )
    !  D~(k'+q') is chf ( nmodes, nmodes, 2*ik   )
    !
    chf = czero
    dyn = czero
    !
    DO na = 1, nat
      DO nb = 1, nat
        total_weight = zero
        DO n1 = -2 * nqc1, 2 * nqc1
          DO n2 = -2 * nqc2, 2 * nqc2
            DO n3 = -2 * nqc3, 2 * nqc3
              !
              ! Sum over r-vectors in the supercell - safe range
              !
              DO i = 1, 3
                r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
              ENDDO
              !
              weight = wscache(n3, n2, n1, nb, na)
              IF (weight > zero) THEN
                !
                ! Find the vector corresponding to r in the original cell
                !
                m1 = MOD(n1 + 1, nqc1)
                IF (m1 <= 0) m1 = m1 + nqc1
                m2 = MOD(n2 + 1, nqc2)
                IF (m2 <= 0) m2 = m2 + nqc2
                m3 = MOD(n3 + 1, nqc3)
                IF (m3 <= 0) m3 = m3 + nqc3
                !
                arg = twopi * (xq(1) * r(1) + xq(2) * r(2) + xq(3) * r(3))
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    dyn(ipol, jpol, na, nb) = dyn(ipol, jpol, na, nb) +   &
                      ifc(m1, m2, m3, ipol, jpol, na, nb) * CMPLX(COS(arg), -SIN(arg), KIND = DP) * weight
                  ENDDO
                ENDDO
              ENDIF
              total_weight = total_weight + weight
            ENDDO
          ENDDO
        ENDDO
        IF (ABS(total_weight - nqc1 * nqc2 * nqc3) > eps8) THEN
          WRITE(stdout, *) total_weight
          CALL errore('dynifc2bloch', 'wrong total_weight', 1)
        ENDIF
      ENDDO
    ENDDO
    !
    DO na = 1, nat
      DO nb = 1, nat
        DO ipol = 1, 3
          DO jpol = 1, 3
            chf((na - 1) * 3 + ipol, (nb - 1) * 3 + jpol) = dyn(ipol, jpol, na, nb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    IF (lpolar .OR. qrpl) THEN
      ! xq has to be in 2pi/a
      CALL rgd_blk(L, nqc1, nqc2, nqc3, nat, chf, xq, tau, epsi, zstar, +1.d0)
      !
    ENDIF
    !
    !  divide by the square root of masses
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)))
        !
        chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
           chf(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) * massfac
        !
      ENDDO
    ENDDO
    !
    !
    !---------------------------------------------------------------------
    !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
    !---------------------------------------------------------------------
    !
    ! Here, we hermitianize only the upper triangular part for cdiagh2
    !
    DO jmode = 1, nmodes
      DO imode = 1, jmode
        chf(imode, jmode) = &
             (chf(imode, jmode) + CONJG(chf(jmode, imode))) * 0.5d0
      ENDDO
    ENDDO
    !
    CALL cdiagh2(nmodes, chf, nmodes, w, cz)
    !
    ! Find the largest element and set it to pure real
    IF (lrot) THEN
      DO imode = 1, nmodes
        DO jmode = 1, nmodes
          norm_vec(jmode) = ABS(cz(jmode, imode))
        END DO
        imode_max(:) = MAXLOC(norm_vec(1:nmodes))
        fac_max(1) = cz(imode_max(1), imode) / norm_vec(imode_max(1))
        cz(1:nmodes, imode) = cz(1:nmodes, imode) * CONJG(fac_max(1))
      ENDDO
    ENDIF
    ! clean noise
    DO jmode = 1,nmodes
      DO imode = 1,nmodes
        IF (ABS(cz(imode, jmode)) < eps12) cz(imode, jmode) = czero
      ENDDO
    ENDDO
    !
    ! DS - Impose phase
    IF (lphase) THEN
      DO jmode = 1,nmodes
        INNER : DO imode = 1,nmodes
          IF (ABS(cz(imode, jmode)) > eps12) THEN
            cz(:, jmode) = cz(:, jmode) * CONJG(cz(imode,jmode))
            cz(:, jmode) = cz(:, jmode) / SQRT(utility_zdotu(CONJG(cz(:, jmode)),cz(:, jmode)))
            EXIT INNER
          ENDIF
        ENDDO INNER
      ENDDO
    ENDIF
    !
    ! cuf(nmodes, nmodes) is rotation matrix (eigenmodes e_k)
    !
    cuf = cz
    eig = w
    !
    CALL stop_clock('DynW2B')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE dynifc2blochf
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dynifc2blochc(nmodes, rws, nrws, xxq, chf)
    !--------------------------------------------------------------------------
    !!
    !! From the IFCs in the format of q2r, find the corresponding
    !! dynamical matrix for a given q point (as in matdyn.x) on the coarse grid
    !!
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg
    USE ions_base, ONLY : tau, nat
    USE global_var,ONLY : ifc, epsi, zstar, wscache, qrpl, L
    USE input,     ONLY : lpolar, nqc1, nqc2, nqc3, fd
    USE ep_constants,  ONLY : twopi, czero, zero, eps8
    USE io_global, ONLY : stdout
    USE longrange, ONLY : rgd_blk
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of modes
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Wigner-Seitz radius
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! Crystal q-point coordinates for the interpolation
    COMPLEX(KIND = DP), INTENT(out) :: chf(nmodes, nmodes)
    !! dyn mat (not divided by the masses)
    !
    ! Local variables
    LOGICAL, SAVE :: first = .TRUE.
    !! First time in routine
    INTEGER :: n1, n2, n3
    !! Q-point grid dimensions
    INTEGER :: m1, m2, m3
    !! Mod of the grid dim
    INTEGER :: i
    !! Cartesian direction
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: jpol
    !! Counter on polarizations
    REAL(KIND = DP) :: total_weight
    !! Sum of the weigths
    REAL(KIND = DP) :: weight
    !! WS weights
    REAL(KIND = DP) :: arg
    !! 2 * pi * r
    REAL(KIND = DP) :: r(3)
    !! Real-space vector
    REAL(KIND = DP) :: xq(3)
    !! Q-point coordinate in Cartesian coordinate.
    REAL(KIND = DP) :: r_ws(3)
    !! Real space vector including fractional translation
    REAL(KIND = DP), EXTERNAL :: wsweight
    !! Wigner-Seitz weights
    COMPLEX(KIND = DP) :: dyn(3, 3, nat, nat)
    !! Dynamical matrix
    !
    xq = xxq
    ! bring xq in cart. coordinates
    CALL cryst_to_cart(1, xq, bg, 1)
    !
    IF (first) THEN
      first = .FALSE.
      DO na = 1, nat
        DO nb = 1, nat
          DO n1 = -2 * nqc1, 2 * nqc1
            DO n2 = -2 * nqc2, 2 * nqc2
              DO n3 = -2 * nqc3, 2 * nqc3
                DO i = 1, 3
                  r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
                  IF (fd) THEN
                    r_ws(i) = r(i)
                  ELSE
                    r_ws(i) = r(i) + tau(i, na) - tau(i, nb)
                  ENDIF
                END DO
                wscache(n3, n2, n1, nb, na) = wsweight(r_ws, rws, nrws)
              ENDDO
            ENDDO
          ENDDO
       ENDDO
      ENDDO
    ENDIF
    !
    chf = czero
    dyn = czero
    !
    DO na = 1, nat
      DO nb = 1, nat
        total_weight = zero
        DO n1 = -2 * nqc1, 2 * nqc1
          DO n2= -2 * nqc2, 2 * nqc2
            DO n3 = -2 * nqc3, 2 * nqc3
              !
              ! Sum over R vectors in the supercell - safe range
              !
              DO i = 1, 3
                r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
              ENDDO
              !
              weight = wscache(n3, n2, n1, nb, na)
              IF (weight > 0.0d0) THEN
                !
                ! Find the vector corresponding to R in the original cell
                !
                m1 = MOD(n1 + 1, nqc1)
                IF (m1 <= 0) m1 = m1 + nqc1
                m2 = MOD(n2 + 1, nqc2)
                IF (m2 <= 0) m2 = m2 + nqc2
                m3 = MOD(n3 + 1, nqc3)
                IF (m3 <= 0) m3 = m3 + nqc3
                !
                arg = twopi * (xq(1) * r(1) + xq(2) * r(2) + xq(3) * r(3))
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    dyn(ipol, jpol, na, nb) = dyn(ipol, jpol, na, nb) +  &
                        ifc(m1, m2, m3, ipol, jpol, na, nb) * CMPLX(COS(arg), -SIN(arg), KIND = DP) * weight
                  ENDDO
                ENDDO
              ENDIF
              total_weight = total_weight + weight
            ENDDO
          ENDDO
        ENDDO
        IF (ABS(total_weight - nqc1 * nqc2 * nqc3) > eps8) THEN
          WRITE(stdout,*) total_weight
          CALL errore ('dynifc2bloch', 'wrong total_weight', 1)
        ENDIF
      ENDDO
    ENDDO
    !
    DO na = 1, nat
      DO nb = 1, nat
        DO ipol = 1, 3
          DO jpol = 1, 3
            chf((na - 1) * 3 + ipol, (nb - 1) * 3 + jpol) = dyn(ipol, jpol, na, nb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    IF (lpolar .OR. qrpl) THEN
      ! xq has to be in 2pi/a
      CALL rgd_blk(L, nqc1, nqc2, nqc3, nat, chf, xq, tau, epsi, zstar, +1.d0)
      !
    ENDIF
    !
    !-----------------------------------------------------------------------------------------
    END SUBROUTINE dynifc2blochc
    !-----------------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dmewan2bloch(nbnd, nrr, cuf, dmef, etf, etf_ks, cfac, use_momentum)
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Wannier representation, find the corresponding
    !!  Dipole in Bloch representation for a given k point
    !!
    !!  input  : number of bands nbnd
    !!           number of WS vectors, coordinates and degeneracy
    !!           Dipole in Wannier representation  cdmew(3,nbnd,nbnd,nrr)
    !!           use_momentum: if .TRUE., use cpmew (momentum) instead of (default: .FALSE.)
    !!
    !!  output : interpolated dipole (v = dH/dk) matrix elements (dmef)
    !!
    !!  SP 09/2016: Optimization
    !!  JN, EK 09/2010
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : cdmew, cpmew
    USE input,         ONLY : eig_read
    USE ep_constants,  ONLY : cone, czero, eps4
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in), OPTIONAL :: use_momentum
    !! If .TRUE., use cpmew instead of cdmew (default: .FALSE.)
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! kpoint number for the interpolation
    REAL(KIND = DP), INTENT(in) :: etf(nbnd)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP), INTENT(in) :: etf_ks(nbnd)
    !! Kohn-Sham eigenvalues
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(out) :: dmef(3, nbnd, nbnd)
    !! interpolated dipole matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr)
    !! Exponential factor
    !
    ! Local variables
    LOGICAL :: use_momentum_
    !! Value of use_momentum
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    !
    COMPLEX( kind=DP ) :: cdmef(3, nbnd, nbnd)
    !! dipole matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: cdmef_tmp(nbnd, nbnd)
    !! dipole matrix elements in Bloch basis, fine mesh
    !
    ! Initialization
    use_momentum_ = .FALSE.
    IF (PRESENT(use_momentum)) use_momentum_ = use_momentum
    cdmef_tmp(:, :) = czero
    dmef(:, :, :)   = czero
    cdmef(:, :, :)  = czero
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  p~(k')    = 1/ndegen(R) sum_R e^{ik'R     } p(R)
    !  p~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} p(R)
    !  Note p~(k') is p^(W)(k') in PRB 74, 195118 (2006) notations
    !
    !  p~(k' )   is chf( nbnd, nbnd, 2*ik-1 )
    !  p~(k'+q') is chf( nbnd, nbnd, 2*ik   )
    !
    ! SUM on ir of: cdmef(1,ibnd,jbnd) = cdmef(1,ibnd,jbnd) + cfac(ir) * cdmew(1,ibnd,jbnd,ir)
    IF (use_momentum_) THEN
      CALL ZGEMV('n', 3 * (nbnd**2), nrr, cone, cpmew, 3 * (nbnd**2), cfac, 1, cone, cdmef, 1)
    ELSE
      CALL ZGEMV('n', 3 * (nbnd**2), nrr, cone, cdmew, 3 * (nbnd**2), cfac, 1, cone, cdmef, 1)
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! p(k') = U(k')^\dagger * p~(k') * U(k')
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    ! Note p(k') is p^(H)(k') in PRB 74, 195118 (2006) notations
    !
    DO ipol = 1, 3
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, cdmef(ipol, :, :), &
                 nbnd, cuf(:, :), nbnd, czero, cdmef_tmp(:, :), nbnd)
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:, :), &
                 nbnd, cdmef_tmp(:, :), nbnd, czero, dmef(ipol, :, :), nbnd)
    ENDDO
    !
    ! Satisfy Phys. Rev. B 62, 4927-4944 (2000) , Eq. (30)
    !
    IF (eig_read) THEN
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          IF (ABS(etf_ks(ibnd) - etf_ks(jbnd)) > eps4) THEN
            dmef(:, ibnd, jbnd) = dmef(:, ibnd, jbnd) * (etf(ibnd) - etf(jbnd)) / &
                                                        (etf_ks(ibnd) - etf_ks(jbnd))
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !--------------------------------------------------------------------------
    END SUBROUTINE dmewan2bloch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE vmewan2bloch(nbnd, nrr, irvec, cuf, vmef, etf, etf_ks, chw, cfac)
    !--------------------------------------------------------------------------
    !!
    !!  From the Velocity matrix elements in Wannier representation, find the corresponding
    !!  MEs in Bloch representation for a given k point
    !!
    !!  input  : nbnd, nrr, irvec, ndegen, xk, cuf, et
    !!
    !!  output : vmef; velocity matrix elements on the fine mesh
    !!
    !!  Adapted from hamwan2bloch by Jesse Noffsinger and Emmanouil Kioupakis
    !!  RM 04/2018: optimized
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : cvmew
    USE cell_base,     ONLY : at, alat
    USE input,         ONLY : eig_read
    USE ep_constants,  ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang, one
    USE low_lvl,       ONLY : degen_sort
    USE rigid,         ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! coordinates of WS points
    REAL(KIND = DP), INTENT(in) :: etf(nbnd)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP), INTENT(in) :: etf_ks(nbnd)
    !! Kohn-Sham eigenvalues
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: chw(nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(out) :: vmef(3, nbnd, nbnd)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: vmef_deg(:, :, :)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh, in the degenerate subspaces
    !
    ! local variables
    LOGICAL :: duplicates
    !! Returns if the bands contains degeneracices for that k-point.
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ideg
    !! Index on degenerations
    INTEGER :: ibndc
    !! Index to count degenerate iband
    INTEGER :: jbndc
    !! Index to count degenerate jband
    INTEGER :: ijbndc
    !! Index to deduce ibndc and jbndc
    INTEGER :: list_dup(nbnd)
    !! List of degenerate eigenvalues
    INTEGER :: info
    !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: deg_dim(:)
    !! Index that keeps track of degeneracies and their dimensionality
    REAL(KIND = DP) :: irvec_tmp(3)
    !! coordinates of WS points for the interpolation, cartesian coordinates
    REAL(KIND = DP), ALLOCATABLE :: w(:)
    !! Eigenvalues
    COMPLEX(KIND = DP) :: chf_a(3, nbnd, nbnd)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP) :: chf_a_tmp(nbnd, nbnd)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP) :: cvmef(3, nbnd, nbnd)
    !! velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: cvmef_tmp(nbnd, nbnd)
    !! velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: cz(:, :)
    !! Eigenvectors
    COMPLEX(KIND = DP) :: cfac_a(nrr, 3)
    !
    ! Initialization
    !
    !$omp master
    CALL start_clock('vmeW2B')
    !$omp end master
    !
    cvmef_tmp(:, :) = czero
    cvmef(:, :, :)  = czero
    vmef(:, :, :)   = czero
    chf_a_tmp(:, :) = czero
    chf_a(:, :, :)  = czero
    irvec_tmp(:)    = zero
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    ! [Eqn. 39 of PRB 74, 195118 (2006)]
    ! A^(W)_{mn,\alpha}(k') = 1/ndegen(R) sum_R e^{ik'R} r_{\alpha}(R)
    !
    ! SUM on ir of: cvmef(1,ibnd,jbnd) = cvmef(1,ibnd,jbnd) + cfac(ir) * cvmew(1,ibnd,jbnd,ir)
    CALL ZGEMV('n', 3 * (nbnd**2), nrr, cone, cvmew, 3 * (nbnd**2), cfac, 1, cone, cvmef, 1)
    !
    ! k-derivative of the Hamiltonian in the Wannier gauge
    ! [Eqn. 38 of PRB 74, 195118 (2006)] or [Eq. 29 of PRB 75, 195121 (2007)]
    !
    ! dH~_{\alpha}(k) = 1/ndegen(R) sum_R i * R_{\alpha} e^{ikR} H(R)
    !
    ! Note that chw stores H(R) / ndegen(R).
    ! Note dH~_{\alpha}(k') is H^(W)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations.
    !
    DO ipol = 1, 3
      DO ir = 1, nrr
        cfac_a(ir, ipol) = alat * DOT_PRODUCT(at(ipol, :), irvec(:, ir)) * cfac(ir)
      ENDDO
    ENDDO
    !
    DO ipol = 1, 3
      CALL ZGEMV('n', nbnd**2, nrr, ci, chw, nbnd**2, cfac_a(:, ipol), 1, czero, chf_a_tmp, 1)
      chf_a(ipol, :, :) = chf_a_tmp(:, :)
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! A^(H)_{mn,\alpha}(k') = U(k')^\dagger A^(W)_{mn,\alpha}(k') U(k')
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    DO ipol = 1, 3
      !
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, cvmef(ipol, :, :), &
                 nbnd, cuf(:, :), nbnd, czero, cvmef_tmp(:, :), nbnd)
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:, :), &
                 nbnd, cvmef_tmp(:, :), nbnd, czero, vmef(ipol, :, :), nbnd)
    ENDDO
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! dH_{\alpha}(k') = U(k')^\dagger dH~_{\alpha}(k') U(k')
    !
    ! Note dH_{\alpha}(k') is H^(H)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    DO ipol = 1, 3
      !
      ! chf_a_tmp(:, :) = matmul( chf_a(ipol,:,:), CONJG(transpose(cuf(:, :))) )
      ! chf_a(ipol,:,:) = matmul(cuf(:, :), chf_a_tmp(:, :) )
      !
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, chf_a(ipol, :, :), &
                 nbnd, cuf(:, :), nbnd, czero, chf_a_tmp(:, :), nbnd)
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:, :), &
                 nbnd, chf_a_tmp(:, :), nbnd, czero, chf_a(ipol, :, :), nbnd)
    ENDDO
    !
    ! velocity matrix elements
    ! [Eqn. 31 of PRB 74, 195118 (2006)]
    ! \hbar v_{mn,\alpha}(k') = H^(H)_{mn,\alpha}(k') &
    !                         - (E^(H)_nk'-E^(H)_mk') * A^(H)_{mn,\alpha}(k')
    !
    ! RM - use etf instead of etf_ks when eig_read=.FALSE.
    IF (eig_read) THEN
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          vmef(:, ibnd, jbnd) = chf_a(:, ibnd, jbnd) - ci * (etf_ks(jbnd) - etf_ks(ibnd)) * vmef(:, ibnd, jbnd)
        ENDDO
      ENDDO
    ELSE
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          vmef(:, ibnd, jbnd) = chf_a(:, ibnd, jbnd) - ci * (etf(jbnd) - etf(ibnd)) * vmef(:, ibnd, jbnd)
        ENDDO
      ENDDO
    ENDIF
    !
    ! Now find and sort degeneracies
    !
    CALL degen_sort(etf, SIZE(etf), duplicates, list_dup)
    !
    ! Count degeneracies and their dimensionality
    IF (duplicates .eqv. .TRUE.) THEN
      ALLOCATE(deg_dim(MAXVAL(list_dup)), STAT = ierr)
      IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error allocating deg_dim(MAXVAL', 1)
      deg_dim = 0
      DO ideg = 1, SIZE(deg_dim)
        DO ibnd = 1, nbnd
          IF (list_dup(ibnd) == ideg) THEN
            deg_dim(ideg) = deg_dim(ideg) + 1
          ENDIF
        ENDDO
      ENDDO
      ! Now allocate matrixes for each degenerate subspace
      DO ideg = 1, SIZE(deg_dim)
        ALLOCATE(vmef_deg(3, deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error allocating vmef_deg(3, deg_dim(ideg), deg_dim', 1)
        ALLOCATE(w(deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error allocating w(deg_dim', 1)
        ALLOCATE(cz(deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error allocating cz(deg_dim(ideg), deg_dim', 1)
        ijbndc = 0
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
            IF ((list_dup(ibnd) == ideg) .AND. (list_dup(jbnd) == ideg)) THEN
              ijbndc = ijbndc + 1
              jbndc = deg_dim(ideg) - MOD(ijbndc, deg_dim(ideg))
              ibndc = INT((ijbndc - 1) / deg_dim(ideg)) + 1
              vmef_deg(:, ibndc, jbndc) = vmef(:, ibnd, jbnd)
            ENDIF
          ENDDO
        ENDDO
        !
        DO ipol = 1, 3
          DO jbndc = 1, deg_dim(ideg)
            DO ibndc = 1, jbndc
              vmef_deg(ipol, ibndc, jbndc) = &
                   (vmef_deg(ipol, ibndc, jbndc) + CONJG(vmef_deg(ipol, jbndc, ibndc))) * 0.5d0
            ENDDO
          ENDDO
          !
          CALL cdiagh2(deg_dim(ideg), vmef_deg(ipol, :, :), deg_dim(ideg), w, cz)
          !
          vmef_deg(ipol, :, :) = zero
          DO ibndc = 1, deg_dim(ideg)
            vmef_deg(ipol, ibndc, ibndc) = w(ibndc)
          ENDDO
        ENDDO !ipol
        !
        ! Now store the values back in vmef
        !
        ijbndc = 0
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
            IF ((list_dup(ibnd) == ideg) .AND. (list_dup(jbnd) == ideg)) THEN
              ijbndc = ijbndc + 1
              jbndc = deg_dim(ideg) - MOD(ijbndc, deg_dim(ideg))
              ibndc = INT((ijbndc - 1) / deg_dim(ideg)) + 1
              vmef(:, ibnd, jbnd) = vmef_deg(:, ibndc, jbndc)
            ENDIF
          ENDDO
        ENDDO
        !
        DEALLOCATE(vmef_deg, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error deallocating vmef_deg', 1)
        DEALLOCATE(w, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error deallocating w', 1)
        DEALLOCATE(cz, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch', 'Error deallocating cz', 1)
        !
      ENDDO !ideg
      !
    ENDIF
    !
    !$omp master
    CALL stop_clock('vmeW2B')
    !$omp end master
    !--------------------------------------------------------------------------
    END SUBROUTINE vmewan2bloch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE vmewan2bloch_batch(nbnd, nrr, irvec, nkf, cuf, vmef, etf, etf_ks, &
      chw, cfac, vmef_tmp, chf_a)
    !--------------------------------------------------------------------------
    !!
    !!  From the Velocity matrix elements in Wannier representation, find the corresponding
    !!  MEs in Bloch representation for a given k point, batched
    !!
    !!  input  : nbnd, nrr, irvec, ndegen, xk, cuf, et
    !!
    !!  output : vmef; velocity matrix elements on the fine mesh
    !!
    !!  Adapted from hamwan2bloch by Jesse Noffsinger and Emmanouil Kioupakis
    !!  RM 04/2018: optimized
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : cvmew
    USE cell_base,     ONLY : at, alat
    USE input,         ONLY : eig_read
    USE ep_constants,  ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang, one
    USE low_lvl,       ONLY : degen_sort
#if defined(__CUDA)
    USE global_var,    ONLY : cublas_h
    USE cublas
#endif
#if defined (__ONEMKL)
    USE epw_onemkl
#endif
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! coordinates of WS points
    INTEGER, INTENT(in) :: nkf
    !! number of k points
    REAL(KIND = DP), INTENT(in) :: etf(nbnd, 2*nkf)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP), INTENT(in) :: etf_ks(nbnd, 2*nkf)
    !! Kohn-Sham eigenvalues
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, 2*nkf)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: chw(nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nbnd, nbnd, 2*nkf)
    !! Rotation matrix U^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(out) :: vmef(3, nbnd, nbnd, 2*nkf)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), INTENT(inout) :: vmef_tmp(nbnd, nbnd, 3, 2*nkf)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), INTENT(inout) :: chf_a(nbnd, nbnd, 2*nkf, 3)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: vmef_deg(:, :, :)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh, in the degenerate subspaces
    !
    ! local variables
    LOGICAL :: duplicates
    !! Returns if the bands contains degeneracices for that k-point.
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ideg
    !! Index on degenerations
    INTEGER :: ibndc
    !! Index to count degenerate iband
    INTEGER :: jbndc
    !! Index to count degenerate jband
    INTEGER :: ijbndc
    !! Index to deduce ibndc and jbndc
    INTEGER :: list_dup(nbnd)
    !! List of degenerate eigenvalues
    INTEGER :: neig
    !! The total number of eigenvalues found
    INTEGER :: info
    !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: ifail(:)
    !! Contains the indices of the eigenvectors that failed to converge
    INTEGER, ALLOCATABLE :: iwork(:)
    !! Integer work array
    INTEGER, ALLOCATABLE :: deg_dim(:)
    !! Index that keeps track of degeneracies and their dimensionality
    REAL(KIND = DP) :: irvec_tmp(3)
    !! coordinates of WS points for the interpolation, cartesian coordinates
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! Real work array
    REAL(KIND = DP), ALLOCATABLE :: w(:)
    !! Eigenvalues
    COMPLEX(KIND = DP) :: chf_a_tmp(nbnd, nbnd, 2*nkf, 3)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP) :: cvmef(3, nbnd, nbnd, 2*nkf)
    !! velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP) :: cvmef_tmp(nbnd, nbnd, 3, 2*nkf)
    !! velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: champ(:)
    !! Complex Hamiltonian packed in upper triangle
    COMPLEX(KIND = DP), ALLOCATABLE :: cwork(:)
    !! Complex work array
    COMPLEX(KIND = DP), ALLOCATABLE :: cz(:, :)
    !! Eigenvectors
    INTEGER :: ik
    !! Counter on coarse k-point grid
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ikq
    !! Paired counter so that q is adjacent to its k
    COMPLEX(KIND = DP), ALLOCATABLE :: ir_cfac(:, :, :)
    !! k derivative of Fourier factor
    REAL(KIND = DP) :: irvec_tmp2
    !! coordinates of WS points for the interpolation, cartesian coordinates
    !
    ! Initialization
    !
    !$omp master
    CALL start_clock('vmeW2B')
    !$omp end master
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    ! [Eqn. 39 of PRB 74, 195118 (2006)]
    ! A^(W)_{mn,\alpha}(k') = 1/ndegen(R) sum_R e^{ik'R} r_{\alpha}(R)
    !
    ! SUM on ir of: cvmef(1,ibnd,jbnd) = cvmef(1,ibnd,jbnd) + cfac(ir) * cvmew(1,ibnd,jbnd,ir)    
    !
    ! k-derivative of the Hamiltonian in the Wannier gauge
    ! [Eqn. 38 of PRB 74, 195118 (2006)] or [Eq. 29 of PRB 75, 195121 (2007)]
    !
    ! dH~_{\alpha}(k) = 1/ndegen(R) sum_R i * R_{\alpha} e^{ikR} H(R)
    !
    ! Note that chw stores H(R) / ndegen(R).
    ! Note dH~_{\alpha}(k') is H^(W)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations.
    !
#if defined(__CUDA)
    !$acc enter data copyin(cvmew) create(cvmef)
    !$acc host_data use_device(cvmew,cfac,cvmef)
    CALL CUBLASZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, cvmew, 3*(nbnd**2), &
      cfac, nrr, czero, cvmef, 3*(nbnd**2))
    !$acc end host_data
    !$acc exit data delete(cvmew)
#elif defined(__ONEMKL)
    !$omp target enter data map(alloc:cvmef) map(to:cvmew)
    !$omp dispatch
    CALL ZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, cvmew, 3*(nbnd**2), &
      cfac, nrr, czero, cvmef, 3*(nbnd**2))
    !$omp target exit data map(delete:cvmew)
#else
    CALL ZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, cvmew, 3*(nbnd**2), &
      cfac, nrr, czero, cvmef, 3*(nbnd**2))
#endif
    !
    ALLOCATE(ir_cfac(nrr, 2*nkf, 3), STAT = ierr)
#if defined(__CUDA)
    !$acc enter data create(ir_cfac)
    !$acc enter data copyin(at,irvec)
    !$acc parallel loop present(ir_cfac,at,irvec) private(irvec_tmp2) collapse(3)
#elif defined(__ONEMKL)
    ! TODO: temporary workaround using OpenMP threads
    !$omp parallel do private(irvec_tmp2) collapse(3)
#endif
    DO ipol = 1, 3
      DO ik = 1, 2*nkf
        DO ir = 1, nrr
          irvec_tmp2 = alat * (at(ipol, 1)*irvec(1, ir) + at(ipol, 2)*irvec(2, ir) &
            + at(ipol, 3)*irvec(3, ir))
          ir_cfac(ir, ik, ipol) = irvec_tmp2 * cfac(ir, ik)
        ENDDO
      ENDDO
    ENDDO
#if defined(__CUDA)
    !$acc end parallel loop
    !$acc exit data delete(at,irvec)
#elif defined(__ONEMKL)
    !$omp end parallel do
    !$omp target update to(ir_cfac)
#endif
    !
#if defined(__CUDA)
    !$acc enter data create(chf_a_tmp) copyin(chw)
    !$acc host_data use_device(chw,ir_cfac,chf_a)
    ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 0, nbnd**2, 1, nrr, ci, &
      chw, nbnd**2, 0, &
      ir_cfac, nrr, nrr, czero, &
      chf_a, nbnd**2, nbnd**2, 6*nkf)
    !$acc end host_data
    !$acc exit data delete(chw)
    !$acc exit data delete(ir_cfac)
#elif defined(__ONEMKL)
    !$omp target enter data map(alloc:chf_a_tmp) map(to:chw)
    !$omp dispatch
    CALL ZGEMV_BATCH_STRIDED('n', nbnd**2, nrr, ci, chw, nbnd**2, 0, ir_cfac, &
      1, nrr, czero, chf_a, 1, nbnd**2, 6*nkf)
    !$omp target exit data map(delete:chw)
    !$omp target exit data map(delete:ir_cfac)
#else
    DO ipol = 1, 3
      DO ik = 1, 2*nkf
        CALL ZGEMV('n', nbnd**2, nrr, ci, chw, nbnd**2, ir_cfac(:, ik, ipol), 1, &
          czero, chf_a(:, :, ik, ipol), 1)
      ENDDO
    ENDDO
#endif
    !
    DEALLOCATE(ir_cfac, STAT = ierr)
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! A^(H)_{mn,\alpha}(k') = U(k')^\dagger A^(W)_{mn,\alpha}(k') U(k')
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    !$acc enter data create(cvmef_tmp)
    !$omp target enter data map(alloc:cvmef_tmp)
    !
#if defined(__CUDA)
    !$acc host_data use_device(cuf,cvmef,cvmef_tmp)
    ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 2, nbnd, 3*nbnd, nbnd, cone, &
      cuf, nbnd, nbnd**2, &
      cvmef, 3*nbnd, 3*(nbnd**2), czero, &
      cvmef_tmp, nbnd, 3*(nbnd**2), 2*nkf)
    !$acc end host_data
#elif defined(__ONEMKL)
    ! Syntax:
    ! call zgemm_batch_strided(transa, transb, m, n, k, alpha, a, lda, stridea, &
    !                          b, ldb, strideb, beta, c, ldc, stridec, batch_size)
    !$omp dispatch
    CALL zgemm_batch_strided('n', 'c', nbnd, 3*nbnd, nbnd, cone, cuf, nbnd, &
      nbnd**2, cvmef, 3*nbnd, 3*nbnd*nbnd, czero, cvmef_tmp, nbnd, 3*nbnd*nbnd, 2*nkf)
#else
    DO ik = 1, 2*nkf
      CALL ZGEMM('n', 'c', nbnd, 3*nbnd, nbnd, cone, cuf(:, :, ik), nbnd, &
        cvmef(:, :, :, ik), 3*nbnd, &
        czero, cvmef_tmp(:, :, :, ik), nbnd)
    ENDDO
#endif
    !
    DO ipol = 1, 3
#if defined(__CUDA)
      !$acc host_data use_device(cuf,cvmef_tmp,vmef_tmp)
      ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 2, nbnd, nbnd, nbnd, cone, &
        cuf, nbnd, nbnd**2, &
        cvmef_tmp(1, 1, ipol, 1), nbnd, 3*(nbnd**2), czero, &
        vmef_tmp(1, 1, ipol, 1), nbnd, 3*(nbnd**2), 2*nkf)
      !$acc end host_data
#elif defined(__ONEMKL)
      !$omp dispatch
      CALL zgemm_batch_strided('n', 'c', nbnd, nbnd, nbnd, cone, cuf, nbnd, &
        nbnd*nbnd, cvmef_tmp(1, 1, ipol, 1), nbnd, 3*(nbnd**2), czero, &
        vmef_tmp(1, 1, ipol, 1), nbnd, 3*(nbnd**2), 2*nkf)
#else
      DO ik = 1, 2*nkf
        CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, cuf(:, :, ik), nbnd, &
          cvmef_tmp(:, :, ipol, ik), nbnd, czero, &
          vmef_tmp(:, :, ipol, ik), nbnd)
      ENDDO
#endif
    ENDDO
    !
    !$acc exit data delete(cvmef,cvmef_tmp)
    !$omp target exit data map(delete:cvmef,cvmef_tmp)
    !
    !$acc update self(vmef_tmp)
    !$omp target update from(vmef_tmp)
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! dH_{\alpha}(k') = U(k')^\dagger dH~_{\alpha}(k') U(k')
    !
    ! Note dH_{\alpha}(k') is H^(H)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    DO ipol = 1, 3
#if defined(__CUDA)
      !$acc host_data use_device(cuf,chf_a,chf_a_tmp)
      ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 2, nbnd, nbnd, nbnd, cone, &
        chf_a(1, 1, 1, ipol), nbnd, nbnd**2, &
        cuf, nbnd, nbnd**2, czero, &
        chf_a_tmp(1, 1, 1, ipol), nbnd, nbnd**2, 2*nkf)
      !$acc end host_data
#elif defined(__ONEMKL)
      !$omp dispatch
      CALL zgemm_batch_strided('n', 'c', nbnd, nbnd, nbnd, cone, &
        chf_a(1, 1, 1, ipol), nbnd, nbnd**2, cuf, nbnd, nbnd**2, czero, &
        chf_a_tmp(1, 1, 1, ipol), nbnd, nbnd**2, 2*nkf)
#else
      DO ik = 1, 2*nkf
        CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, chf_a(:, :, ik, ipol), &
          nbnd, cuf(:, :, ik), nbnd, czero, chf_a_tmp(:, :, ik, ipol), nbnd)
      ENDDO
#endif
      !
#if defined(__CUDA)
      !$acc host_data use_device(cuf,chf_a,chf_a_tmp)
      ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 0, nbnd, nbnd, nbnd, cone, &
        cuf, nbnd, nbnd**2, &
        chf_a_tmp(1, 1, 1, ipol), nbnd, nbnd**2, czero, &
        chf_a(1, 1, 1, ipol), nbnd, nbnd**2, 2*nkf)
      !$acc end host_data
#elif defined(__ONEMKL)
      !$omp dispatch
      CALL zgemm_batch_strided('n', 'n', nbnd, nbnd, nbnd, cone, &
        cuf, nbnd, nbnd**2, chf_a_tmp(1, 1, 1, ipol), nbnd, nbnd**2, czero, &
        chf_a(1, 1, 1, ipol), nbnd, nbnd**2, 2*nkf)
#else
      DO ik = 1, 2*nkf
        CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:, :, ik), &
          nbnd, chf_a_tmp(:, :, ik, ipol), nbnd, czero, chf_a(:, :, ik, ipol), nbnd)
      ENDDO
#endif
    ENDDO
    !$acc exit data delete(chf_a_tmp)
    !$omp target exit data map(delete:chf_a_tmp)
    !
    !$acc update self(chf_a)
    !$omp target update from(chf_a)
    !
    !$omp parallel do schedule(dynamic) &
    !$omp private(ipol,ibnd,jbnd,duplicates,list_dup) &
    !$omp private(deg_dim,ideg,vmef_deg,ifail,iwork,w,rwork,champ,cwork,cz,ierr) &
    !$omp private(ijbndc,ibndc,jbndc,info,neig)
    DO ik = 1, 2*nkf
      !
      ! velocity matrix elements
      ! [Eqn. 31 of PRB 74, 195118 (2006)]
      ! \hbar v_{mn,\alpha}(k') = H^(H)_{mn,\alpha}(k') &
      !                         - (E^(H)_nk'-E^(H)_mk') * A^(H)_{mn,\alpha}(k')
      !
      ! RM - use etf instead of etf_ks when eig_read=.FALSE.
      IF (eig_read) THEN
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
            vmef(:, ibnd, jbnd, ik) = vmef_tmp(ibnd, jbnd, :, ik)
            vmef(:, ibnd, jbnd, ik) = chf_a(ibnd, jbnd, ik, :) &
              - ci * (etf_ks(jbnd, ik) - etf_ks(ibnd, ik)) * vmef(:, ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ELSE
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
            vmef(:, ibnd, jbnd, ik) = vmef_tmp(ibnd, jbnd, :, ik)
            vmef(:, ibnd, jbnd, ik) = chf_a(ibnd, jbnd, ik, :) &
              - ci * (etf(jbnd, ik) - etf(ibnd, ik)) * vmef(:, ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ENDIF
      !
      ! Now find and sort degeneracies
      !
      CALL degen_sort(etf(:, ik), SIZE(etf(:, ik)), duplicates, list_dup)
      !
      ! Count degeneracies and their dimensionality
      IF (duplicates .eqv. .TRUE.) THEN
        ALLOCATE(deg_dim(MAXVAL(list_dup)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating deg_dim(MAXVAL)', 1)
        deg_dim = 0
        DO ideg = 1, SIZE(deg_dim)
          DO ibnd = 1, nbnd
            IF (list_dup(ibnd) == ideg) THEN
              deg_dim(ideg) = deg_dim(ideg) + 1
            ENDIF
          ENDDO
        ENDDO
        ! Now allocate matrixes for each degenerate subspace
        DO ideg = 1, SIZE(deg_dim)
          ALLOCATE(vmef_deg(3, deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating vmef_deg(3, deg_dim(ideg), deg_dim', 1)
          ALLOCATE(ifail(deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating ifail(deg_dim', 1)
          ALLOCATE(iwork(5 * deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating iwork(5 * deg_dim', 1)
          ALLOCATE(w(deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating w(deg_dim', 1)
          ALLOCATE(rwork(7 * deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating rwork(7 * deg_dim', 1)
          ALLOCATE(champ(deg_dim(ideg) * (deg_dim(ideg) + 1) / 2), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating champ(deg_dim(ideg) * (deg_dim', 1)
          ALLOCATE(cwork(2 * deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating cwork(2 * deg_dim', 1)
          ALLOCATE(cz(deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error allocating cz(deg_dim(ideg), deg_dim', 1)
          ijbndc = 0
          DO ibnd = 1, nbnd
            DO jbnd = 1, nbnd
              IF ((list_dup(ibnd) == ideg) .AND. (list_dup(jbnd) == ideg)) THEN
                ijbndc = ijbndc + 1
                jbndc = deg_dim(ideg) - MOD(ijbndc, deg_dim(ideg))
                ibndc = INT((ijbndc - 1) / deg_dim(ideg)) + 1
                vmef_deg(:, ibndc, jbndc) = vmef(:, ibnd, jbnd, ik)
              ENDIF
            ENDDO
          ENDDO
          !
          DO ipol = 1, 3
            DO jbndc = 1, deg_dim(ideg)
              DO ibndc = 1, jbndc
                champ(ibndc + (jbndc - 1) * jbndc / 2) = &
                      (vmef_deg(ipol, ibndc, jbndc) + CONJG(vmef_deg(ipol, jbndc, ibndc))) * 0.5d0
              ENDDO
            ENDDO
            !
            CALL ZHPEVX('V', 'A', 'U', deg_dim(ideg), champ , zero, zero, &
                      0, 0, -one, neig, w, cz, deg_dim(ideg), cwork, rwork, iwork, ifail, info)
            !
            vmef_deg(ipol, :, :) = zero
            DO ibndc = 1, deg_dim(ideg)
              vmef_deg(ipol, ibndc, ibndc) = w(ibndc)
            ENDDO
          ENDDO !ipol
          !
          ! Now store the values back in vmef
          !
          ijbndc = 0
          DO ibnd = 1, nbnd
            DO jbnd = 1, nbnd
              IF ((list_dup(ibnd) == ideg) .AND. (list_dup(jbnd) == ideg)) THEN
                ijbndc = ijbndc + 1
                jbndc = deg_dim(ideg) - MOD(ijbndc, deg_dim(ideg))
                ibndc = INT((ijbndc - 1) / deg_dim(ideg)) + 1
                vmef(:, ibnd, jbnd, ik) = vmef_deg(:, ibndc, jbndc)
              ENDIF
            ENDDO
          ENDDO
          !
          DEALLOCATE(vmef_deg, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating vmef_deg', 1)
          DEALLOCATE(ifail, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating ifail', 1)
          DEALLOCATE(iwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating iwork', 1)
          DEALLOCATE(w, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating w', 1)
          DEALLOCATE(rwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating rwork', 1)
          DEALLOCATE(champ, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating champ', 1)
          DEALLOCATE(cwork, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating cwork', 1)
          DEALLOCATE(cz, STAT = ierr)
          IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating cz', 1)
          !
        ENDDO !ideg
        !
        DEALLOCATE(deg_dim, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2bloch_batch', 'Error deallocating deg_dim', 1)
        !
      ENDIF
    ENDDO
    !$omp end parallel do
    !
    !$omp master
    CALL stop_clock('vmeW2B')
    !$omp end master
    !--------------------------------------------------------------------------
    END SUBROUTINE vmewan2bloch_batch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE vmewan2blochp(xxq, nmodes, nrr_q, irvec_q, ndegen_q, cuf, vmefp, wf, rws, nrws)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the phonon velocity by computing the q derivative of
    !! the dynamical matrix.
    !! This routines is required for adaptative broadening.
    !! Samuel Ponce & Francesco Macheda
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : rdw, epsi, zstar, wscache, ifc
    USE cell_base,     ONLY : at, alat, bg
    USE input,         ONLY : use_ws, lpolar, lifc, nqc1, nqc2, nqc3, fd
    USE ep_constants,  ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang, one, eps8
    USE ions_base,     ONLY : tau, nat
    USE io_global,     ONLY : stdout
    USE low_lvl,       ONLY : degen_sort
    USE longrange,     ONLY : rgd_blk_der
    USE rigid,         ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of modes (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr_q
    !! number of WS points
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! coordinates of WS points
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, nat, nat)
    !! degeneracy of WS points
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors when lifc == .TRUE.
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! q-point coordinates for the interpolation
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nmodes, nmodes)
    !! Rotation matrix e^\dagger, fine mesh
    REAL(KIND = DP), INTENT(in) :: wf(nmodes)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Real space Wigner-Seitz vector when lifc == .TRUE.
    COMPLEX(KIND = DP), INTENT(out) :: vmefp(3, nmodes, nmodes)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: vmef_deg(:, :, :)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh, in the degenerate subspaces
    !
    ! Local variables
    LOGICAL, SAVE :: first = .TRUE.
    !! First entrance [used when lifc == .TRUE.]
    LOGICAL :: duplicates
    !! Returns if the modes contains degeneracices for that q-point.
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: imode
    !! Counter on band index
    INTEGER :: jmode
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: jpol
    !! Counter on polarization
    INTEGER :: na
    !! Local index
    INTEGER :: nb
    !! Local index
    INTEGER :: list_dup(nmodes)
    !! List of degenerate modes
    INTEGER :: ideg
    !! Index on degenerations
    INTEGER :: imodec
    !! Index to count degenerate imode
    INTEGER :: jmodec
    !! Index to count degenerate jmode
    INTEGER :: ijmodec
    !! Index to deduce imodec and jmodec
    INTEGER :: n1, n2, n3
    !! WS dimensions [used when lifc == .TRUE.]
    INTEGER :: m1, m2, m3
    !! Find corresponding vector [used when lifc == .TRUE.]
    INTEGER :: i
    !! Cartesian direction
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: deg_dim(:)
    !! Index that keeps track of degeneracies and their dimensionality
    REAL(KIND = DP) :: irvec_tmp(3)
    !! coordinates of WS points for the interpolation, cartesian coordinates
    REAL(KIND = DP) :: xq(3)
    !! Coordinates q-point
    REAL(KIND = DP) :: rdotq
    !! $$\mathbf{r}\cdot\mathbf{q}
    REAL(KIND = DP) :: total_weight
    !! Total WS weight [used when lifc == .TRUE.]
    REAL(KIND = DP) :: weight
    !! WS weight [used when lifc == .TRUE.]
    REAL(KIND = DP) :: arg
    !! Argument of the exp. [used when lifc == .TRUE.]
    REAL(KIND = DP) :: r(3)
    !! Real space vector [used when lifc == .TRUE.]
    REAL(KIND = DP) :: r_ws(3)
    !! Real WS point [used when lifc == .TRUE.]
    REAL(KIND = DP), EXTERNAL :: wsweight
    !! WS weight [used when lifc == .TRUE.]
    COMPLEX(KIND = DP) :: cfac
    !! Complex prefactor for Fourier transform.
    COMPLEX(KIND = DP) :: chf_a(3, nmodes, nmodes)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP) :: chf_a_tmp(nmodes, nmodes)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(KIND = DP) :: dyn_a(3, 3, 3, nat, nat)
    !! Temp dyn mat. [used when lifc == .TRUE.]
    REAL(KIND = DP), ALLOCATABLE :: w(:)
    !! Workspace for storing the eigenvalue of the velocity matrix
    COMPLEX(KIND = DP), ALLOCATABLE :: cz(:, :)
    !! Workspace for storing the eigenvector of the velocity matrix
    !
    ! Initialization
    !
    CALL start_clock('vmeW2Bp')
    !
    vmefp(:, :, :)  = czero
    chf_a_tmp(:, :) = czero
    chf_a(:, :, :)  = czero
    irvec_tmp(:)    = zero
    xq = xxq
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine q meshes
    !----------------------------------------------------------
    ! Adapted for D Eqn. 38 of PRB 74, 195118 (2006)
    ! d_{\alpha} D_{\mu\nu}(q) = 1/ndegen(R) sum_R i*R_{\alpha} e^{iqR} D(R)
    !
    IF (lifc) THEN
      ! bring xq in cart. coordinates
      CALL cryst_to_cart(1, xq, bg, 1)
      !
      IF (first) THEN
        first = .FALSE.
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = -2 * nqc1, 2 * nqc1
              DO n2 = -2 * nqc2, 2 * nqc2
                DO n3 = -2 * nqc3, 2 * nqc3
                  DO i = 1, 3
                    r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
                    ! SM: Added for reading IFC's are from finite displacement
                    IF (fd) THEN
                      r_ws(i) = r(i)
                    ELSE
                      r_ws(i) = r(i) + tau(i, na) - tau(i, nb)
                    ENDIF
                  ENDDO
                  wscache(n3, n2, n1, nb, na) = wsweight(r_ws, rws, nrws)
                ENDDO
              ENDDO
            ENDDO
          ENDDO ! nb
        ENDDO ! na
      ENDIF
      !
      dyn_a(:, :, :, :, :) = czero
      !
      DO na = 1, nat
        DO nb = 1, nat
          total_weight = zero
          DO n1 = -2 * nqc1, 2 * nqc1
            DO n2 = -2 * nqc2, 2 * nqc2
              DO n3 = -2 * nqc3, 2 * nqc3
                !
                ! Sum over R vectors in the supercell - safe range
                DO i = 1, 3
                  r(i) = n1 * at(i, 1) + n2 * at(i, 2) + n3 * at(i, 3)
                ENDDO
                !
                weight = wscache(n3, n2, n1, nb, na)
                IF (weight > zero) THEN
                  !
                  ! Find the vector corresponding to R in the orginial cell
                  m1 = MOD(n1 + 1, nqc1)
                  IF (m1 <= 0) m1 = m1 + nqc1
                  m2 = MOD(n2 + 1, nqc2)
                  IF (m2 <= 0) m2 = m2 + nqc2
                  m3 = MOD(n3 + 1, nqc3)
                  IF (m3 <= 0) m3 = m3 + nqc3
                  !
                  arg = twopi * (xq(1) * r(1) + xq(2) * r(2) + xq(3) * r(3))
                  DO ipol = 1, 3
                    DO jpol = 1, 3
                      dyn_a(:, ipol, jpol, na, nb) = dyn_a(:, ipol, jpol, na, nb) -  &
                        ifc(m1, m2, m3, ipol, jpol, na, nb) * ci * alat * r(:) * CMPLX(COS(arg), -SIN(arg), KIND = DP) * weight
                    ENDDO
                  ENDDO
                ENDIF
                total_weight = total_weight + weight
              ENDDO
            ENDDO
          ENDDO
          IF (ABS(total_weight - nqc1 * nqc2 * nqc3) > eps8) THEN
            WRITE(stdout,*) total_weight
            CALL errore ('vmewan2blochp', 'wrong total_weight', 1)
          END IF
        ENDDO ! nb
      ENDDO ! na
      !
      DO na = 1, nat
        DO nb = 1, nat
          DO ipol = 1, 3
            DO jpol = 1, 3
              chf_a(:, (na - 1) * 3 + ipol, (nb - 1) * 3 + jpol) = dyn_a(:, ipol, jpol, na, nb)
            ENDDO
          ENDDO
        ENDDO ! nb
      ENDDO ! na
      !
    ELSE ! lifc
      IF (use_ws) THEN
        DO ir = 1, nrr_q
          rdotq = twopi * DOT_PRODUCT(xq, DBLE(irvec_q(:, ir)))
          irvec_tmp(:) = alat * MATMUL(at, DBLE(irvec_q(:, ir)))
          DO na = 1, nat
            DO nb = 1, nat
              IF (ndegen_q(ir, na, nb) > 0) THEN
                cfac = EXP(ci * rdotq) / DBLE(ndegen_q(ir, na, nb))
                ! To map atom coordinate to mode basis.
                DO ipol = 1, 3
                  chf_a(ipol, 3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
                  chf_a(ipol, 3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) &
                  + ci * irvec_tmp(ipol) * cfac * rdw(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, ir)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE ! use_ws
        DO ir = 1, nrr_q
          rdotq = twopi * DOT_PRODUCT(xq, DBLE(irvec_q(:, ir)))
          cfac = EXP(ci * rdotq) / DBLE(ndegen_q(ir, 1, 1))
          irvec_tmp(:) = alat * MATMUL(at, DBLE(irvec_q(:, ir)))
          DO ipol = 1, 3
            chf_a(ipol, :, :) = chf_a(ipol, :, :) + &
                ci * irvec_tmp(ipol) * cfac * rdw(:, :, ir)
          ENDDO
        ENDDO
      ENDIF
      ! bring xq in cart. coordinates (needed for rgd_blk call)
      CALL cryst_to_cart(1, xq, bg, 1)
    ENDIF ! lifc
    !
    ! add the long-range term to D(q)
    ! FIXME: add support for quadrupoles
    IF (lpolar) THEN
      ! xq has to be in 2pi/a
      CALL rgd_blk_der(nqc1, nqc2, nqc3, nat, chf_a, xq, tau, epsi, zstar, +1.d0)
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! d_{\alpha}D_{\mu\nu}^{(H)}(q) = e_{\mu\nu}(q)^\dagger d_{\alpha}D_{\mu\nu}(q) e_{\mu\nu}(q)
    ! Note that the e_{\mu\nu}(q) = cuf are already mass scaled with 1.0/sqrt(amass(ityp(na)))
    !
    DO ipol = 1, 3
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, chf_a(ipol, :, :), &
                 nmodes, cuf(:, :), nmodes, czero, chf_a_tmp(:, :), nmodes)
      CALL ZGEMM('c', 'n', nmodes, nmodes, nmodes, cone, cuf(:, :), &
                 nmodes, chf_a_tmp(:, :), nmodes, czero, chf_a(ipol, :, :), nmodes)
    ENDDO
    !
    ! velocity matrix elements
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
        vmefp(:, imode, jmode) = chf_a(:, imode, jmode)
      ENDDO
    ENDDO
    !
    ! Now find and sort degeneracies
    !
    CALL degen_sort(wf, SIZE(wf), duplicates, list_dup)
    !
    ! Count degeneracies and their dimensionality
    !
    IF (duplicates .eqv. .TRUE.) THEN
      ALLOCATE(deg_dim(MAXVAL(list_dup)), STAT = ierr)
      IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error allocating deg_dim(MAXVAL', 1)
      deg_dim = 0
      DO ideg = 1, size(deg_dim)
        DO imode = 1, nmodes
          IF (list_dup(imode) == ideg) deg_dim(ideg) = deg_dim(ideg) + 1
        ENDDO
      ENDDO
      ! Now allocate matrixes for each degenerate subspace
      DO ideg = 1, SIZE(deg_dim)
        ALLOCATE(vmef_deg(3, deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error allocating vmef_deg(3, deg_dim(ideg), deg_dim', 1)
        ALLOCATE(w(deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error allocating w(deg_dim', 1)
        ALLOCATE(cz(deg_dim(ideg), deg_dim(ideg)), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error allocating cz(deg_dim(ideg), deg_dim', 1)
        ijmodec = 0
        DO imode = 1, nmodes
          DO jmode = 1, nmodes
            IF ((list_dup(imode) == ideg) .AND. (list_dup(jmode) == ideg)) THEN
              ijmodec = ijmodec+1
              jmodec = deg_dim(ideg) - MOD(ijmodec, deg_dim(ideg))
              imodec = INT((ijmodec - 1) / deg_dim(ideg)) + 1
              vmef_deg(:, imodec, jmodec) = vmefp(:, imode, jmode)
            ENDIF
          ENDDO
        ENDDO
        !
        DO ipol = 1, 3
          DO jmodec = 1, deg_dim(ideg)
            DO imodec = 1, jmodec
              vmef_deg(ipol, imodec, jmodec) = &
               (vmef_deg(ipol, imodec, jmodec) + CONJG(vmef_deg(ipol, jmodec, imodec))) * 0.5d0
            ENDDO
          ENDDO
          !
          CALL cdiagh2(deg_dim(ideg), vmef_deg(ipol, :, :), deg_dim(ideg), w, cz)
          !
          vmef_deg(ipol, :, :) = zero
          DO imodec = 1, deg_dim(ideg)
            vmef_deg(ipol, imodec, imodec) = w(imodec)
          ENDDO
        ENDDO !ipol
        !
        ! Now store the values back in vmefp
        !
        ijmodec = 0
        DO imode = 1, nmodes
          DO jmode = 1, nmodes
            IF ((list_dup(imode) == ideg) .AND. (list_dup(jmode) == ideg)) THEN
              ijmodec = ijmodec + 1
              jmodec = deg_dim(ideg) - MOD(ijmodec, deg_dim(ideg))
              imodec = INT((ijmodec - 1) / deg_dim(ideg)) + 1
              vmefp(:, imode, jmode) = vmef_deg(:, imodec, jmodec)
            ENDIF
          ENDDO
        ENDDO
        !
        DEALLOCATE(vmef_deg, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error deallocating vmef_deg', 1)
        DEALLOCATE(w, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error deallocating w', 1)
        DEALLOCATE(cz, STAT = ierr)
        IF (ierr /= 0) CALL errore('vmewan2blochp', 'Error deallocating cz', 1)
        !
      ENDDO !ideg
    ENDIF
    !
    CALL stop_clock('vmeW2Bp')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE vmewan2blochp
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE rrwan2bloch(nbnd, nrr, cfac, rrf)
    !--------------------------------------------------------------------------
    !!
    !! Fourier interpolation of the position matrix elements
    !! SP - 06/22
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : crrw
    USE ep_constants,  ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang, one
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(out) :: rrf(3, nbnd, nbnd)
    !! interpolated position matrix element in Bloch basis, fine mesh
    !
    ! Initialization
    !
    !----------------------------------------------------------
    !  STEP 1: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !CALL cryst_to_cart(1, vector, at, 1)
    !
    !$omp master
    CALL start_clock('rrW2B')
    !$omp end master
    !
    CALL ZGEMV('n', 3 * (nbnd**2), nrr, cone, crrw, 3 * (nbnd**2), cfac, 1, czero, rrf, 1)
    !
    !DO ibnd = 1, nbnd
    !  rrf(:, ibnd, ibnd) = rrf(:, ibnd, ibnd) - vector * alat
    !  print*,'rrf ',rrf(:, ibnd, ibnd)
    !  print*,'vector alat ',vector * alat
    !ENDDO
    !CALL cryst_to_cart(1, vector, bg, -1)
    !
    !$omp master
    CALL stop_clock('rrW2B')
    !$omp end master
    !--------------------------------------------------------------------------
    END SUBROUTINE rrwan2bloch
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE rrwan2bloch_batch(nbnd, nrr, nkf, cfac, rrf)
    !--------------------------------------------------------------------------
    !!
    !! Fourier interpolation of the position matrix elements, batched
    !! SP - 06/22
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : crrw
    USE ep_constants,  ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang, one
#if defined(__CUDA)
    USE global_var,    ONLY : cublas_h
    USE cublas
#endif 
#if defined (__ONEMKL)
    USE epw_onemkl
#endif
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: nkf
    !! number of k-points, fine mesh
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, 2*nkf)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(out) :: rrf(3, nbnd, nbnd, 2*nkf)
    !! interpolated position matrix element in Bloch basis, fine mesh
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ierr
    !! Error status
    !
    ! Initialization
    !
    !----------------------------------------------------------
    !  STEP 1: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !$omp master
    CALL start_clock('rrW2B')
    !$omp end master
    !
#if defined(__CUDA)
    !$acc data copyin(crrw) present(cfac,rrf)
    !$acc host_data use_device(crrw,cfac,rrf)
    CALL CUBLASZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, crrw, 3*(nbnd**2), &
      cfac, nrr, czero, rrf, 3*(nbnd**2))
    !$acc end host_data
    !$acc end data
    !$acc update self(rrf)
#elif defined(__ONEMKL)
    !$omp target enter data map(to:crrw)
    !$omp dispatch
    CALL ZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, crrw, 3*(nbnd**2), &
      cfac, nrr, czero, rrf, 3*(nbnd**2))
    !$omp target exit data map(delete:crrw)
    !$omp target update from(rrf)
#else
    CALL ZGEMM('n', 'n', 3*(nbnd**2), 2*nkf, nrr, cone, crrw, 3*(nbnd**2), &
      cfac, nrr, czero, rrf, 3*(nbnd**2))
#endif
    !
    !$omp master
    CALL stop_clock('rrW2B')
    !$omp end master
    !--------------------------------------------------------------------------
    END SUBROUTINE rrwan2bloch_batch
    !--------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatf, nbnd, nrr_k, dims, nat)
    !---------------------------------------------------------------------------
    !!
    !! Even though this is for phonons, we use the same notations
    !! adopted for the electronic case (nmodes -> nmodes etc)
    !!
    USE kinds,            ONLY : DP
    USE input,            ONLY : etf_mem, use_ws, lopt_w2b
    USE global_var,       ONLY : epmatwp, irn_start, irn_stop, irg_start, &
                                 irg_stop, nirg_loc
    USE ep_constants,     ONLY : twopi, ci, czero, cone
    USE io_var,           ONLY : iunepmatwp, iunepmatwp2
    USE mp,               ONLY : mp_sum, mp_bcast
    USE mp_world,         ONLY : world_comm, mpime
    USE mp_pools,         ONLY : inter_pool_comm
    USE io_global,        ONLY : ionode_id, ionode
    USE io,               ONLY : rwepmatw
    USE wannier2bloch_opt,ONLY : ephwan2blochp_opt
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, &
                                 MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, &
                                 MPI_IN_PLACE, MPI_SUM
#endif
#if defined(__CUDA)
    USE cublas
#endif
#if defined (__ONEMKL)
    USE epw_onemkl
#endif    
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT(in) :: nrr_g
    !! Number of phononic WS points
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of WS points
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: nat
    !! Is equal to the number of atoms if use_ws == .TRUE. or 1 otherwise
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Number of degeneracy of WS points
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands
    INTEGER, INTENT(in) ::  nrr_k
    !! Number of electronic WS points
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(KIND = DP), INTENT(inout) :: epmatf(nbnd, nbnd, nrr_k, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! Local variables
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: iw
    !! Wannier function index
    INTEGER :: iw2
    !! Wannier function index
    INTEGER :: irn
    !! Combined WS and atom index
    INTEGER :: ir_start
    !! Starting ir for this pool
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: na
    !! Atom index
    INTEGER :: imode
    !! Number of modes
    INTEGER :: ipol
    !! Polarization index
    INTEGER :: ir_loc
    !! Local counter for ir in this pool
    INTEGER :: irr_k
    !! Counter on k-point
    INTEGER :: diff
    !! Difference between starting and ending on master core
    INTEGER :: add
    !! Additional element
    INTEGER :: ierr
    !! Error status
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER :: epmat_block_dtype
    !! MPI block datatype for parallel reading epmatwp
#endif
    INTEGER :: lsize_epmatw
    !! Number of e-ph matrix elements to be read from the file at a time
    REAL(KIND = DP) :: rdotk
    !! Exponential for the FT
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Factor for the FT
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatw(:, :, :, :)
    !! Buffer for el-ph matrix elements
    INTEGER :: ndim
    !! Last dimension of epmatw
    !
    IF (lopt_w2b) THEN
      ! Use optimized implementation in wannier2bloch_opt.f90
      CALL ephwan2blochp_opt(nmodes, xxq, irvec_g, nrr_g, epmatf, nbnd, nrr_k)
      RETURN
    ENDIF
    !
    CALL start_clock('ephW2Bp')
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(R_e,q') = 1/ndegen(R_p) sum_R_p e^{iq'R_p} g(R_e,R_p)
    !
    !  g~(R_e,q') is epmatf(nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    ! Assuming that ws_indexes_distribution was called before
    ! so that the irg_start and irg_stop are set
    ALLOCATE(cfac(irg_stop - irg_start + 1), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp', 'Error allocating cfac', 1)
    !
    cfac(:) = czero
    DO ir = irg_start, irg_stop
      ! note xxq is assumed to be already in cryst coord
      rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
      !
      ir_loc = ir - irg_start + 1
      cfac(ir_loc) = EXP(ci * rdotk)
    ENDDO
    !
    IF (etf_mem == 0) then
      ! Fourier transform of g to fine k mesh
#if defined(__CUDA)
      !$acc enter data copyin(cfac)
      !$acc host_data use_device(epmatwp,epmatf,cfac)
      CALL CUBLASZGEMM('n', 'n', nbnd * nbnd * nrr_k * nmodes, 1, nirg_loc, cone, &
                      epmatwp, nbnd * nbnd * nrr_k * nmodes, &
                      cfac, nirg_loc, czero, epmatf, nbnd * nbnd * nrr_k * nmodes)
      !$acc end host_data
      !$acc exit data delete(cfac)
      !$acc update self(epmatf)
#elif defined(__ONEMKL)
      !$omp target data map(to:cfac)
      !$omp dispatch
      CALL ZGEMV('n', nbnd * nbnd * nrr_k * nmodes, nirg_loc, cone, &
                epmatwp, nbnd * nbnd * nrr_k * nmodes, cfac, 1, czero, epmatf,1)
      !$omp end target data
      !$omp target update from(epmatf)
#else
      CALL ZGEMV('n', nbnd * nbnd * nrr_k * nmodes, nirg_loc, cone, &
                epmatwp, nbnd * nbnd * nrr_k * nmodes, cfac, 1, czero, epmatf, 1)
#endif
    ELSE ! etf_mem == 1
      !
      ! 1-1. Allocate epmatw as a buffer for reading the epmatwp file
      !
      ! Although this should almost never be problematic (see explaination below)
      IF (use_ws) THEN
        ndim = 3
      ELSE
        ndim = 1
      ENDIF
      !
      lsize_epmatw = nbnd * nbnd * nrr_k * ndim
      ALLOCATE(epmatw(nbnd, nbnd, nrr_k, ndim), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwan2blochp', 'Error allocating epmatw', 1)
      epmatw(:, :, :, :) = czero
      !
      ! 1-2. Since MPI_FILE_READ_AT_ALL is a collective operation, 
      !      we set the add variable to make the ranks that have no data to read
      !      do nothing but still be part of the MPI_FILE_READ_AT_ALL call.
      IF (ionode) THEN
        diff = irn_stop - irn_start
      ENDIF
      CALL mp_bcast(diff, ionode_id, inter_pool_comm)
      !
      ! If you are the last cpu with less element
      IF (irn_stop - irn_start /= diff) THEN
        add = 1
      ELSE
        add = 0
      ENDIF
      !
      ! 1-3. Zeros epmatf. This is important because ZAXPY will be used later. 
      epmatf(:, :, :, :) = czero
      !
      DO irn = irn_start, irn_stop + add
        ! 2-1. Read the irn-th slice of the epmatwp file
#if defined(__MPI)
        !
        ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or
        !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
        !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
        lrepmatw = 2_MPI_OFFSET_KIND * 8_MPI_OFFSET_KIND * &
                                       INT(lsize_epmatw, KIND = MPI_OFFSET_KIND) * &
                                       INT(irn - 1, KIND = MPI_OFFSET_KIND)
        !
        ! SP: mpi seek is used to set the position at which we should start
        ! reading the file. It is given in bits.
        ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
        !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ.
        !        Here we want non blocking because not all the process have the same nb of irn.
        !
        !CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
        !IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
        !CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
        !IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
        !
        IF (add == 1 .AND. irn == irn_stop + add) THEN
          CALL MPI_FILE_READ_AT_ALL(iunepmatwp2, 0_MPI_OFFSET_KIND, epmatw, 0, MPI_DOUBLE_COMPLEX, &
                                MPI_STATUS_IGNORE, ierr)
        ELSE
          CALL MPI_FILE_READ_AT_ALL(iunepmatwp2, lrepmatw, epmatw, lsize_epmatw, MPI_DOUBLE_COMPLEX, &
                                MPI_STATUS_IGNORE, ierr)
        ENDIF
        IF (ierr /= 0) CALL errore('ephwan2blochp', 'error in MPI_FILE_READ_AT_ALL', 1)
        !
        IF (add == 1 .AND. irn == irn_stop + add) CYCLE
#else
        CALL rwepmatw(epmatw, nbnd, nrr_k, ndim, irn, iunepmatwp, -1)
#endif
        !
        ! 2-2. Dividing epmatw by the degeneracy factor
        IF (use_ws) THEN
          ir = (irn - 1) / nat + 1
          na = MOD(irn - 1, nat) + 1
          !
          DO ipol = 1, 3
            DO irr_k = 1, nrr_k
              DO iw2 = 1, nbnd
                DO iw = 1, nbnd
                  IF (ndegen_g(iw, ir, na) > 0) THEN
                    epmatw(iw, iw2, irr_k, ipol) = epmatw(iw, iw2, irr_k, ipol) / &
                      REAL(ndegen_g(iw, ir, na), KIND = DP)
                  ELSE
                    epmatw(iw, iw2, irr_k, ipol) = czero
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ir = (irn - 1) / nmodes + 1
          !
          DO irr_k = 1, nrr_k
            DO iw2 = 1, nbnd
              DO iw = 1, nbnd
                epmatw(iw, iw2, irr_k, 1) = epmatw(iw, iw2, irr_k, 1) / &
                  REAL(ndegen_g(1, ir, 1), KIND = DP)
              ENDDO
            ENDDO
          ENDDO
        ENDIF ! use_ws
        !
        ! 2-3. Fourier transform of g to fine k mesh
        ir_loc = ir - irg_start + 1
        !
        IF (use_ws) THEN
          imode = 3 * (na - 1) + 1
        ELSE
          imode = MOD(irn - 1, nmodes) + 1
        ENDIF
        !
        CALL ZAXPY(lsize_epmatw, cfac(ir_loc), epmatw, 1, epmatf(1, 1, 1, imode), 1)
      ENDDO ! irn
      !
      DEALLOCATE(epmatw, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwan2blochp', 'Error deallocating epmatw', 1)
      !
    ENDIF ! etf_mem
    !
    !
    ! NOTE: MPI_ALLREDUCE with MPI_IN_PLACE and custum datatype seem to have an issue with MPI from NVHPC
    ! CALL mp_sum(epmatf, inter_pool_comm)
    DO imode = 1, nmodes
      CALL mp_sum(epmatf(:, :, :, imode), inter_pool_comm)
    ENDDO
    !
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp', 'Error deallocating cfac', 1)
    !
    CALL stop_clock('ephW2Bp')
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE ephwan2blochp
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp_dist(nmodes, xxq, irvec_g, nrr_g, epmatf, nbnd, nrr_k)
    !---------------------------------------------------------------------------
    !!
    !! This subroutine optimizes memory usage by distributing the memory 
    !! into the processes. It only works with etf_mem == 0.
    !!
    !! Zhe Liu 08/2024
    !!
    USE kinds,            ONLY : DP
    USE input,            ONLY : etf_mem
    USE global_var,       ONLY : epmatwp_dist, epmatwp
    USE ep_constants,     ONLY : twopi, ci, czero, cone
    USE mp,               ONLY : mp_sum
    USE mp_pools,         ONLY : inter_pool_comm
    USE parallelism,      ONLY : para_bounds
#if defined(__CUDA)
    USE cublas
#endif
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, MPI_SUM
#endif
#if defined (__ONEMKL)
    USE epw_onemkl
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT(in) :: nrr_g
    !! Number of phononic WS points
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of WS points
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd, nrr_k, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! Local variables
    INTEGER :: irn
    !! Combined WS and atom index
    INTEGER :: irn_start
    !! Starting irn for this pool
    INTEGER :: irn_stop
    !! Ending irn for this pool
    INTEGER :: imode
    !! Number of modes
    REAL(KIND = DP) :: rdotk
    !! Exponential for the FT
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: ir_start
    !! Starting ir in this pool
    INTEGER :: ir_stop
    !! Ending ir in this pool
    INTEGER :: nir_loc
    !! Number of ir in this pool
    INTEGER :: ir_loc
    !! Local ir in this pool
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Factor for the FT
    INTEGER :: ierr
    !! Error status
#if defined(__MPI)
    INTEGER :: epmat_block_dtype
    !! MPI block datatype for parallel reading epmatwp
#endif
    !
    IF (etf_mem /= 0) THEN
      CALL errore("ephwan2blochp_dist", "ephwan2blochp_dist must be called only if etf_mem == 0", 1)
    ENDIF
    !
    CALL start_clock('ephW2Bp')
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(R_e,q') = 1/ndegen(R_p) sum_R_p e^{iq'R_p} g(R_e,R_p)
    !
    !  g~(R_e,q') is epmatf(nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    ! SP: Because nrr_g can be quite small, we do a combined parallelization on WS vector and atoms
    !
    CALL para_bounds(irn_start, irn_stop, nrr_g * nmodes)
    !
    ! Courtesy of Zhe Liu, we can still use ZGEMV with dividing nrr_g * nmodes.
    !
    ir_start = (irn_start - 1) / nmodes + 1
    ir_stop = (irn_stop - 1) / nmodes + 1
    !
    nir_loc = ir_stop - ir_start + 1
    ALLOCATE(cfac(nir_loc), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_dist', 'Error allocating cfac', 1)
    cfac = czero
    !
    DO ir = ir_start, ir_stop
      ! note xxq is assumed to be already in cryst coord
      rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
      !
      ir_loc = ir - ir_start + 1
      cfac(ir_loc) = EXP(ci * rdotk)
    ENDDO
    !
#if defined(__CUDA)
    !$acc enter data copyin(cfac)
    !$acc host_data use_device(epmatwp,epmatf,cfac)
    CALL CUBLASZGEMM('n', 'n', nbnd * nbnd * nrr_k * nmodes, 1, nir_loc, cone, &
                     epmatwp, nbnd * nbnd * nrr_k * nmodes, &
                     cfac, nir_loc, czero, epmatf, nbnd * nbnd * nrr_k * nmodes)
    !$acc end host_data
    !$acc exit data delete(cfac)
    !$acc update self(epmatf)
#elif defined(__ONEMKL)
    !$omp target data map(to:cfac)
    !$omp dispatch
    CALL ZGEMV('n', nbnd * nbnd * nrr_k * nmodes, nir_loc, cone, &
               epmatwp, nbnd * nbnd * nrr_k * nmodes, cfac, 1, czero, epmatf, 1)
    !$omp end target data
    !$omp target update from(epmatf)
    !CALL ZGEMM('n', 'n', nbnd * nbnd * nrr_k * nmodes, 1, nir_loc, cone, &
    !  epmatwp, nbnd * nbnd * nrr_k * nmodes, &
    !  cfac, nir_loc, czero, epmatf, nbnd * nbnd * nrr_k * nmodes)
#else
    CALL ZGEMV('n', nbnd * nbnd * nrr_k * nmodes, nir_loc, cone, &
               epmatwp, nbnd * nbnd * nrr_k * nmodes, cfac, 1, czero, epmatf, 1)
#endif
    !
    ! NOTE: MPI_ALLREDUCE with MPI_IN_PLACE and custum datatype seem to have an issue with MPI from NVHPC
    ! CALL mp_sum(epmatf, inter_pool_comm)
    DO imode = 1, nmodes
      CALL mp_sum(epmatf(:, :, :, imode), inter_pool_comm)
    ENDDO
    !
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_dist', 'Error deallocating cfac', 1)
    !
    CALL stop_clock('ephW2Bp')
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE ephwan2blochp_dist
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch(nbnd, nrr, epmatw, cufkk, cufkq, epmatf, nmodes, cfac, q, rrf, long)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the electron-phonon
    !! matrix elements
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : twopi, ci, czero, cone, one
    USE input,         ONLY : lpolar, nqc1, nqc2, nqc3
    USE global_var,    ONLY : epsi, zstar, qrpl
    USE longrange,     ONLY : rgd_blk_epw
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: epmatw(nbnd, nbnd, nrr, nmodes)
    !! e-p matrix in Wannier representation
    COMPLEX(KIND = DP), INTENT(in) :: cufkk(nbnd, nbnd)
    !! rotation matrix U(k)^\dagger, fine k mesh
    COMPLEX(KIND = DP), INTENT(in) :: cufkq(nbnd, nbnd)
    !! rotation matrix U(k+q)^\dagger, fine q mesh
    COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid in Cartesian coordinates
    COMPLEX(KIND = DP), INTENT(inout) :: rrf(3, nbnd, nbnd)
    !! Position operator in Wannier basis.
    LOGICAL, INTENT(IN), OPTIONAL :: long
    !! If true, add long-range contribution. (Default: true)
    !
    ! Local variables
    LOGICAL :: long_
    !! Value of long
    INTEGER :: imode
    !! Counter on  phonon modes
    COMPLEX(KIND = DP) :: eptmp(nbnd, nbnd)
    !! Temporary variable
    !
    !$omp master
    CALL start_clock('ephW2B')
    !$omp end master
    !
    long_ = .TRUE.
    IF (PRESENT(long)) THEN
      long_ = long
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(k',q') = 1/ndegen(R_e) sum_R_e e^{ik'R_e} g(R_e,q')
    !
    !  g~(k',q') is epmatf(nmodes, nmodes, ik)
    !  every pool works with its own subset of k points on the fine grid
    !
    epmatf(:, :, :) = czero
    !
    DO imode = 1, nmodes
      CALL zgemv('n', nbnd**2, nrr, cone, epmatw(:, :, :, imode), nbnd**2, cfac, 1, cone, epmatf(:, :, imode), 1)
    ENDDO
    !
    IF (long_) THEN
      IF (lpolar .OR. qrpl) THEN
        CALL rgd_blk_epw(nqc1, nqc2, nqc3, q, epmatf, nbnd, nmodes, epsi, zstar, rrf, one)
      ENDIF
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g(k',q') = U(k'+q') * g~(k',q') * U(k')^\dagger
    !
    !  the two zgemm calls perform the following ops:
    !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
    !
    DO imode = 1, nmodes
      !
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
           nbnd, epmatf(:, :, imode), nbnd, czero, eptmp, nbnd)
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, eptmp, &
           nbnd, cufkk, nbnd, czero, epmatf(:, :, imode), nbnd)
      !
    ENDDO
    !
    !$omp master
    CALL stop_clock('ephW2B')
    !$omp end master
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE ephwan2bloch
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch_batch(nbnd, nrr, epmatw, nkf, cuf, epmatf, &
      nmodes, cfac, q, rrf, long)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the electron-phonon
    !! matrix elements, batched
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : twopi, ci, czero, cone, one
    USE input,         ONLY : lpolar, nqc1, nqc2, nqc3
    USE global_var,    ONLY : epsi, zstar, qrpl
    USE longrange,     ONLY : rgd_blk_epw
#if defined(__CUDA)
    USE global_var,    ONLY : cublas_h
    USE cublas
#endif
#if defined (__ONEMKL)
    USE epw_onemkl
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, 2*nkf)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: epmatw(nbnd, nbnd, nrr, nmodes)
    !! e-p matrix in Wannier representation
    INTEGER, INTENT(in) :: nkf
    !! number of k points in fine grid
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nbnd, nbnd, 2*nkf)
    !! rotation matrix U(k)^\dagger and U(k+q)^\dagger fine k and k+q mesh
    COMPLEX(KIND = DP), INTENT(inout) :: epmatf(nbnd, nbnd, 2*nkf, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid
    COMPLEX(KIND = DP), INTENT(inout) :: rrf(3, nbnd, nbnd, 2*nkf)
    !! Position operator in Wannier basis.
    LOGICAL, INTENT(IN), OPTIONAL :: long
    !! If true, add long-range contribution. (Default: true)
    !
    ! Local variables
    LOGICAL :: long_
    !! Value of long
    INTEGER :: imode
    !! Counter on phonon modes
    COMPLEX(KIND = DP) :: epmatf_k(nbnd, nbnd, nmodes)
    !! Temporary array for the k slice of epmatf
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatf_mode(:, :, :)
    !! Temporary array for the mode slice of epmatf
    INTEGER :: ik
    !! Counter on coarse k-point grid
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ikq
    !! Paired counter so that q is adjacent to its k
    INTEGER :: ierr
    !! Error status
    !
    !$omp master
    CALL start_clock('ephW2B')
    !$omp end master
    !
    long_ = .TRUE.
    IF (PRESENT(long)) THEN
      long_ = long
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(k',q') = 1/ndegen(R_e) sum_R_e e^{ik'R_e} g(R_e,q')
    !
    !  g~(k',q') is epmatf(nmodes, nmodes, ik)
    !  every pool works with its own subset of k points on the fine grid
    !
#if defined(__CUDA)
    !$acc update device(epmatw)
    !$acc host_data use_device(epmatw,cfac,epmatf)
    ierr = CUBLASZGEMMSTRIDEDBATCHED(cublas_h, 0, 0, nbnd**2, 2*nkf, nrr, &
      cone, epmatw, nbnd**2, (nbnd**2)*nrr, cfac, nrr, 0, czero, &
      epmatf, nbnd**2, (nbnd**2)*2*nkf, nmodes)
    !$acc end host_data
    !$acc update self(epmatf)
    ! TODO: add error checks for CUDA routines
#elif defined(__ONEMKL)
    !$omp target update to(epmatw)
    !$omp dispatch
    CALL zgemm_batch_strided('n', 'n', nbnd**2, 2*nkf, nrr, cone, epmatw, nbnd**2, &
      (nbnd**2)*nrr, cfac, nrr, 0, czero, epmatf, nbnd**2, (nbnd**2)*2*nkf, nmodes)
    !$omp target update from(epmatf)
#else
    epmatf = czero
    !$omp parallel do private(imode,ik,ikk,ikq) collapse(2)
    DO imode = 1, nmodes
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        CALL ZGEMV('n', nbnd**2, nrr, cone, epmatw(:, :, :, imode), nbnd**2, &
          cfac(:, ikk), 1, cone, epmatf(:, :, ikk, imode), 1)
        !
        CALL ZGEMV('n', nbnd**2, nrr, cone, epmatw(:, :, :, imode), nbnd**2, &
          cfac(:, ikq), 1, cone, epmatf(:, :, ikq, imode), 1)
      ENDDO
    ENDDO
    !$omp end parallel do
#endif
    !
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      IF (long_) THEN
        IF (lpolar .OR. qrpl) THEN
          epmatf_k = epmatf(:, :, ikk, :)
          CALL rgd_blk_epw(nqc1, nqc2, nqc3, q, epmatf_k, nbnd, nmodes, epsi, zstar, rrf(:, :, :, ikk), one)
          epmatf(:, :, ikk, :) = epmatf_k
        ENDIF
      ENDIF
    ENDDO
    !
    !$acc update device(epmatf)
    !$omp target update to(epmatf)
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g(k',q') = U(k'+q') * g~(k',q') * U(k')^\dagger
    !
    !  the two zgemm calls perform the following ops:
    !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
    !
    ALLOCATE(epmatf_mode(nbnd, nbnd, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_batch', 'Error allocating epmatf_mode', 1)
    !$acc enter data create(epmatf_mode)
    !$omp target enter data map(alloc:epmatf_mode)
    !
    DO imode = 1, nmodes
#if defined(__CUDA)
      !$acc host_data use_device(cuf,epmatf,epmatf_mode)
      ierr = CUBLASZGEMMSTRIDEDBATCHED( &
        cublas_h, 0, 0, nbnd, nbnd, nbnd, cone, &
        cuf(1, 1, 2), nbnd, 2*(nbnd**2), &
        epmatf(1, 1, 1, imode), nbnd, 2*(nbnd**2), czero, &
        epmatf_mode, nbnd, nbnd**2, nkf)
      !$acc end host_data
      !
      !$acc host_data use_device(epmatf_mode,cuf,epmatf)
      ierr = CUBLASZGEMMSTRIDEDBATCHED( &
        cublas_h, 0, 2, nbnd, nbnd, nbnd, cone, &
        epmatf_mode, nbnd, nbnd**2, &
        cuf(1, 1, 1), nbnd, 2*(nbnd**2), czero, &
        epmatf(1, 1, 1, imode), nbnd, nbnd**2, nkf)
      !$acc end host_data
      ! !$acc update self(epmatf)
#elif defined(__ONEMKL)
      !$omp dispatch
      CALL zgemm_batch_strided( &
        'n', 'n', nbnd, nbnd, nbnd, cone, &
        cuf(1, 1, 2), nbnd, 2*(nbnd**2), &
        epmatf(1, 1, 1, imode), nbnd, 2*(nbnd**2), czero, &
        epmatf_mode, nbnd, nbnd**2, nkf)
      !
      !$omp dispatch
      CALL zgemm_batch_strided( &
        'n', 'c', nbnd, nbnd, nbnd, cone, &
        epmatf_mode, nbnd, nbnd**2, &
        cuf(1, 1, 1), nbnd, 2*(nbnd**2), czero, &
        epmatf(1, 1, 1, imode), nbnd, nbnd**2, nkf)
#else
      !$omp parallel do private(ikk,ikq)
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:, :, ikq), nbnd, &
          epmatf(:, :, ikk, imode), nbnd, czero, epmatf_mode(:, :, ik), nbnd)
        !
        CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, epmatf_mode(:, :, ik), nbnd, &
          cuf(:, :, ikk), nbnd, czero, epmatf(:, :, ik, imode), nbnd)
      ENDDO
      !$omp end parallel do
#endif
    ENDDO
    !
    !$acc update self(epmatf)
    !$omp target update from(epmatf)
    !
    !$acc exit data delete(epmatf_mode)
    !$omp target exit data map(delete:epmatf_mode)
    DEALLOCATE(epmatf_mode)
    !
    !$omp master
    CALL stop_clock('ephW2B')
    !$omp end master
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE ephwan2bloch_batch
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch_mem(nbnd, nrr, epmatw, cufkk, cufkq, epmatf, cfac, dims)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the electron-phonon
    !! matrix elements
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : twopi, ci, czero, cone
    USE input,         ONLY : use_ws
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE. Is equal to 1 otherwise.
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, dims, dims)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: epmatw(nbnd, nbnd, nrr)
    !! e-p matrix in Wannier representation
    COMPLEX(KIND = DP), INTENT(in) :: cufkk(nbnd, nbnd)
    !! rotation matrix U(k)^\dagger, fine k mesh
    COMPLEX(KIND = DP), INTENT(in) :: cufkq(nbnd, nbnd)
    !! rotation matrix U(k+q)^\dagger, fine k mesh
    COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! Local variables
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: iw
    !! Counter on Wannier functions
    INTEGER :: iw2
    !! Counter on Wannier functions
    COMPLEX(KIND = DP) :: eptmp(nbnd, nbnd)
    !! Temporary variable
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(k',q') = 1/ndegen(R_e) sum_R_e e^{ik'R_e} g(R_e,q')
    !
    !  g~(k',q') is epmatf(nmodes, nmodes, ik)
    !  every pool works with its own subset of k points on the fine grid
    !
    epmatf(:, :) = czero
    !
    IF (use_ws) THEN
      DO iw2 = 1, dims
        DO iw = 1, dims
          DO ir = 1, nrr
           epmatf(iw, iw2) = epmatf(iw, iw2) + epmatw(iw, iw2, ir) * cfac(ir, iw, iw2)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      CALL ZGEMV('n', nbnd**2, nrr, cone, epmatw(:, :, :), nbnd**2, cfac(:, 1, 1), 1, cone, epmatf(:, :), 1)
    ENDIF
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g(k',q') = U(k'+q') * g~(k',q') * U(k')^\dagger
    !
    !  the two zgemm calls perform the following ops:
    !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
    !
    !
    CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
               nbnd, epmatf(:, :), nbnd, czero, eptmp, nbnd)
    CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, eptmp, &
               nbnd, cufkk, nbnd, czero, epmatf(:, :), nbnd)
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE ephwan2bloch_mem
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp_mem(imode, nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatf, nbnd, nrr_k, dims, nat)
    !---------------------------------------------------------------------------
    !!
    !! Even though this is for phonons, I use the same notations
    !! adopted for the electronic case (nmodes->nmodes etc)
    !!
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : twopi, ci, czero
    USE io_files,         ONLY : prefix, tmp_dir
    USE mp,               ONLY : mp_sum
    USE mp_world,         ONLY : world_comm
    USE mp_pools,         ONLY : inter_pool_comm
    USE input,            ONLY : use_ws
    USE parallelism,      ONLY : para_bounds
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, &
                                 MPI_MODE_RDONLY,MPI_INFO_NULL
#endif
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: imode
    !! Current mode
    INTEGER, INTENT(in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT(in) :: nrr_g
    !! Number of phononic WS points
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE.
    !! Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: nat
    !! Is equal to the number of atoms if use_ws == .TRUE. or 1 otherwise.
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of WS points
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Number of degeneracy of WS points
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands
    INTEGER, INTENT(in) ::  nrr_k
    !! Number of electronic WS points
    REAL(KIND = DP) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd, nrr_k)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! File name
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: ir_start
    !! Starting ir for this pool
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: iw
    !! Counter on Wannier functions
    INTEGER :: iunepmatwp2
    !! Return the file unit
    INTEGER :: na
    !! Index on atom
    INTEGER :: ierr
    !! Error status
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw2
    !! Offset to tell where to start reading the file
#else
    INTEGER(KIND = 8) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER(KIND = 8) :: lrepmatw2
    !! Offset to tell where to start reading the file
#endif
    !
    REAL(KIND = DP) :: rdotk
    !! Exponential for the FT
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:, :)
    !! Factor for the FT
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatw(:, :, :)
    !! El-ph matrix elements
    !
    CALL start_clock('ephW2Bp')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(R_e,q') = 1/ndegen(R_p) sum_R_p e^{iq'R_p} g(R_e,R_p)
    !
    !  g~(R_e,q') is epmatf(nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    CALL para_bounds(ir_start, ir_stop, nrr_g)
    !
#if defined(__MPI)
    filint = TRIM(tmp_dir) // TRIM(prefix) // '.epmatwp'
    CALL MPI_FILE_OPEN(inter_pool_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatwp2, ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'error in MPI_FILE_OPEN', 1)
#endif
    !
    ALLOCATE(cfac(nrr_g, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'Error allocating cfac', 1)
    cfac(:, :) = czero
    !
    IF (use_ws) THEN
      na = (imode - 1) / 3 + 1
      !
      DO ir = ir_start, ir_stop
        !
        ! note xxq is assumed to be already in cryst coord
        rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
        DO iw = 1, dims
          IF (ndegen_g(iw, ir, na) > 0) &
            cfac(ir, iw) = EXP(ci * rdotk) / DBLE(ndegen_g(iw, ir, na))
        ENDDO
      ENDDO
    ELSE
      DO ir = ir_start, ir_stop
        !
        ! note xxq is assumed to be already in cryst coord
        rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
        cfac(ir, 1) = EXP(ci * rdotk) / DBLE(ndegen_g(1, ir, 1))
      ENDDO
    ENDIF
    !
    ALLOCATE(epmatw(nbnd, nbnd, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'Error allocating epmatw', 1)
    epmatw(:, :, :) = czero
    !
#if defined(__MPI)
    lrepmatw2 = 2_MPI_OFFSET_KIND * INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                    INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                    INT(nrr_k, KIND = MPI_OFFSET_KIND)
#endif
    !
    DO ir = ir_start, ir_stop
      !
      ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or
      !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
      !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
#if defined(__MPI)
      lrepmatw = 2_MPI_OFFSET_KIND * INT(nbnd  ,KIND = MPI_OFFSET_KIND) * &
                                     INT(nbnd  ,KIND = MPI_OFFSET_KIND) * &
                                     INT(nrr_k ,KIND = MPI_OFFSET_KIND) * &
                                     INT(nmodes,KIND = MPI_OFFSET_KIND) * &
               8_MPI_OFFSET_KIND *  (INT(ir    ,KIND = MPI_OFFSET_KIND) - 1_MPI_OFFSET_KIND) + &
               2_MPI_OFFSET_KIND *   INT(nbnd  ,KIND = MPI_OFFSET_KIND) * &
                                     INT(nbnd  ,KIND = MPI_OFFSET_KIND) * &
                                     INT(nrr_k ,KIND = MPI_OFFSET_KIND) * &
               8_MPI_OFFSET_KIND *  (INT(imode ,KIND = MPI_OFFSET_KIND) - 1_MPI_OFFSET_KIND)
      !
      ! SP: mpi seek is used to set the position at which we should start
      ! reading the file. It is given in bits.
      ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
      !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ.
      !        Here we want non blocking because not all the process have the same nb of ir.
      !
      !CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
      !IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
      !CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
      !IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
      CALL MPI_FILE_READ_AT_ALL(iunepmatwp2, lrepmatw, epmatw, lrepmatw2, &
        MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'error in MPI_FILE_READ_AT_ALL', 1)
#endif
      !
      IF (use_ws) THEN
        DO iw = 1, dims
          CALL ZAXPY(nbnd * nrr_k, cfac(ir, iw), epmatw(iw, :, :), 1, epmatf(iw, :, :), 1)
        ENDDO
      ELSE
        CALL ZAXPY(nbnd * nbnd * nrr_k, cfac(ir, 1), epmatw, 1, epmatf, 1)
      ENDIF
      !
    ENDDO
    DEALLOCATE(epmatw, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'Error deallocating epmatw', 1)
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'Error deallocating cfac', 1)
    !
    CALL mp_sum(epmatf, inter_pool_comm)
    !
#if defined(__MPI)
    CALL MPI_FILE_CLOSE(iunepmatwp2, ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_mem', 'error in MPI_FILE_CLOSE', 1)
#endif
    !
    CALL stop_clock('ephW2Bp')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ephwan2blochp_mem
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dwwan2blochp(nbnd, dwmatef, cuf, dwmatpf, nmodes)
    !--------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the Debye-Waller
    !! matrix elements
    !!
    !! dwmatef(i, j, jdir, imode)
    !! = i * <\psi_i(k)| [dV_SCF^{imode}(q=0), p_jdir] |\psi_j(k)>
    !! where imode = 3 * (iatm - 1) + idir is in Cartesian basis.
    !!
    !! For q == Gamma, the output dwmatpf should be
    !! dwmatpf(i,j,kmode)
    !! = sum_{iatm,idir,jdir} i * <\psi_i(k)| [dV_SCF^{iatm, idir}, p_jdir] |\psi_j(k)>
    !!                        * CONJG(cuf(imode, kmode)) * cuf(jmode, kmode)
    !! = sum_{iatm,idir,jdir} dwmatef(i, j, jdir, imode)
    !!                        * CONJG(cuf(imode, kmode)) * cuf(jmode, kmode)
    !! where jmode = 3 * (iatm - 1) + jdir
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ci, czero, cone
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    COMPLEX(KIND = DP), INTENT(in) :: dwmatef(nbnd, nbnd, 3, nmodes)
    !! e-p matrix in electron Bloch, phonon Cartesian representation
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nmodes, nmodes)
    !! phonon eigenmode
    COMPLEX(KIND = DP), INTENT(out) :: dwmatpf(nbnd, nbnd, nmodes)
    !! e-p matrix in electron Bloch, phonon eigenmode representation
    !
    ! Local variables
    INTEGER :: iatm
    !! Counter on atoms
    INTEGER :: idir
    !! Counter on Cartesian direction
    INTEGER :: jdir
    !! Counter on Cartesian direction
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: jmode
    !! Counter on phonon modes
    INTEGER :: kmode
    !! Counter on phonon modes
    !
    CALL start_clock('dwW2Bp')
    !
    !----------------------------------------------------------
    !  STEP 1: rotate phonon basis from Cartesian to eigenmode
    !----------------------------------------------------------
    !
    ! dwmatpf(:,:,kmode) = dwmatef(:,:,jdir,imode)
    !                    * CONJG(cuf(imode, kmode)) * cuf(jmode, kmode)
    !
    ! imode, jmode : Cartesian (Wannier) representation
    ! kmode        : eigenmode (Bloch) representation
    !
    dwmatpf = czero
    !
    DO imode = 1, nmodes
      idir = MOD(imode - 1, 3) + 1
      iatm = (imode - 1) / 3 + 1
      DO jdir = 1, 3
        jmode = 3 * (iatm - 1) + jdir
        !
        DO kmode = 1, nmodes
          dwmatpf(:, :, kmode) = dwmatpf(:, :, kmode) &
            + dwmatef(:, :, jdir, imode) * CONJG(cuf(imode, kmode)) * cuf(jmode, kmode)
        ENDDO
        !
      ENDDO
    ENDDO
    !
    CALL stop_clock('dwW2Bp')
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE dwwan2blochp
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE sthwan2blochp(nbnd, sthmatef, cuf, sthmatpf, nmodes)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the Sternheimer
    !! matrix elements
    !!
    !! sthmatpf(:, :, imode) = sum_{jmode, kmode} sthmatef(:, :, kmode, jmode)
    !!                       * CONJG(cuf(kmode, imode)) * cuf(jmode, imode)
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ci, czero, cone
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    COMPLEX(KIND = DP), INTENT(in) :: sthmatef(nbnd, nbnd, nmodes, nmodes)
    !! Sterheimer matrix in electron Bloch, phonon Cartesian representation
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nmodes, nmodes)
    !! phonon eigenmode
    COMPLEX(KIND = DP), INTENT(out) :: sthmatpf(nbnd, nbnd, nmodes)
    !! Sterheimer matrix in electron Bloch, phonon eigenmode representation
    !
    ! Local variables
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: jmode
    !! Counter on phonon modes
    INTEGER :: kmode
    !! Counter on phonon modes
    CALL start_clock('sthW2Bp')
    !
    sthmatpf = czero
    !
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
        DO kmode = 1, nmodes
          sthmatpf(:, :, imode) = sthmatpf(:, :, imode) &
            + sthmatef(:, :, kmode, jmode) &
            * CONJG(cuf(kmode, imode)) * cuf(jmode, imode)
        ENDDO
      ENDDO
    ENDDO
    !
    CALL stop_clock('sthW2Bp')
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE sthwan2blochp
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE dgwan2blochp(nbnd, dgmatef, cuf, dgmatpf, nmodes)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the delta g
    !! matrix elements
    !!
    !! dgmatpf(:, :, imode) = sum_{jmode} dgmatef(:, :, jmode)
    !!                      * cuf(jmode, imode)
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ci, czero, cone
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    COMPLEX(KIND = DP), INTENT(in) :: dgmatef(nbnd, nbnd, nmodes)
    !! delta g matrix in electron Bloch, phonon Cartesian representation
    COMPLEX(KIND = DP), INTENT(in) :: cuf(nmodes, nmodes)
    !! phonon eigenmode
    COMPLEX(KIND = DP), INTENT(out) :: dgmatpf(nbnd, nbnd, nmodes)
    !! delta g matrix in electron Bloch, phonon eigenmode representation
    !
    ! Local variables
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: jmode
    !! Counter on phonon modes
    CALL start_clock('dgW2Bp')
    !
    dgmatpf = czero
    !
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
          dgmatpf(:, :, imode) = dgmatpf(:, :, imode) &
            + dgmatef(:, :, jmode) * cuf(jmode, imode)
      ENDDO
    ENDDO
    !
    CALL stop_clock('dgW2Bp')
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE dgwan2blochp
    !---------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  END MODULE wannier2bloch
  !--------------------------------------------------------------------------
