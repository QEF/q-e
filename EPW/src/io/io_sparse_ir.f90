  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE io_sparse_ir
  !----------------------------------------------------------------------
  !!
  !! HM: This is adapted from sparse_ir_fortran/sparse_ir_io.f90,
  !!     which is in sparse_ir_fortran repository released under the CC0-1.0 license.
  !!     See also the GitHub repository:
  !!     https://github.com/SpM-lab/sparse-ir-fortran
  !!
  !!     PLEASE do not carelessly rename arrays and do not change the orders of 
  !!     declarations. We should keep the original names and orders of arrays 
  !!     to keep it easy to compare subroutines with the ones in 
  !!     sparse_ir_fortran repository. 
  !!
  !! 07/2022 adapted from   (SHA) 6ecef444884c48e4a155d6d66966ad37322be249
  !! 09/2022 updated according to b257776a0441b6755add1610c5826dbef675e386
  !! 05/2024 updated according to babfec6783c614253bcaa80347e35071b21835a9
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : ionode_id
  USE mp_world,        ONLY : mpime
  USE mp,              ONLY : mp_barrier, mp_bcast
  USE mp_global,       ONLY : inter_pool_comm
  USE supercond_common,ONLY : nlambda, ndigit
  USE sparse_ir,       ONLY : IR, init_ir
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: read_ir_epw
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  FUNCTION read_ir_epw(unit, beta, positive_only) RESULT(obj)
  !-----------------------------------------------------------------------
  !!
  !! This function calls read_v1 
  !! to read the file including the ir-basis objects.
  !!
  !
  INTEGER, INTENT(IN) :: unit
  !! Unit number
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  LOGICAL, INTENT(IN), OPTIONAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  !
  TYPE(IR) :: obj
  !! to contain all the ir-basis objects from the file
  CHARACTER(LEN = 100) :: tmp_str
  !! dummy for characters
  INTEGER :: version
  !! version number
  !
  IF (mpime == ionode_id) THEN
    READ(unit,*) tmp_str, version
    !write(*, *) "Invalid version number", version
    !stop "Stopping..."
    IF (version .NE. 1) &
    CALL errore('read_ir_epw', 'Error while reading data of ir objects: Invalid version number',version)
  ENDIF
  CALL mp_bcast(version, ionode_id, inter_pool_comm)
  IF (version == 1) then
    IF ((.NOT. PRESENT(positive_only))) then
      obj = read_v1(unit, beta)
    ELSE
      obj = read_v1(unit, beta, positive_only)
    ENDIF
  ENDIF
  !-----------------------------------------------------------------------
  END FUNCTION read_ir_epw
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION read_v1(unit, beta, positive_only) RESULT(obj)
  !-----------------------------------------------------------------------
  !!
  !! This function reads the file to get the ir-basis objects.
  !! (version 1)
  !!
  !
  INTEGER, INTENT(IN) :: unit
  !! Unit number
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  LOGICAL, INTENT(IN), OPTIONAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  !
  REAL(KIND = DP), PARAMETER :: rtol = 1e-20
  !! 
  TYPE(IR) :: obj
  !! to contain all the ir-basis objects from the file
  CHARACTER(LEN = 100) :: tmp_str
  !! dummy for characters
  INTEGER :: i
  !! counter
  INTEGER :: l
  !! counter
  INTEGER :: t
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: size
  !! total number of IR basis functions (size of s)
  INTEGER :: ntau
  !! total number of sampling points of imaginary time
  INTEGER :: nfreq_f
  !! total number of sampling Matsubara freqs (Fermionic)
  INTEGER :: nfreq_b
  !! total number of sampling Matsubara freqs (Bosonic)
  INTEGER :: nomega
  !! total number of sampling points of real frequency
  INTEGER :: ierr
  !! Error status
  INTEGER, ALLOCATABLE :: freq_f(:)
  !! integer part of sampling Matsubara freqs (Fermion)
  INTEGER, ALLOCATABLE :: freq_b(:)
  !! integer part of sampling Matsubara freqs (Boson)
  REAL(KIND = DP) :: rtmp
  !! dummy for real variables
  REAL(KIND = DP) :: rtmp2
  !! dummy for real variables
  REAL(KIND = DP) :: lambda
  !! lambda = 10^{nlambda}, 
  !! which determines maximum sampling point of real frequency
  REAL(KIND = DP) :: eps
  !! eps = 10^{-ndigit}
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  REAL(KIND = DP), ALLOCATABLE :: tau(:)
  !! sampling points of imaginary time (dimenisonless)
  REAL(KIND = DP), ALLOCATABLE :: omega(:)
  !! sampling points of real frequency (dimenisonless)
  REAL(KIND = DP), ALLOCATABLE :: u(:, :)
  !! dimensionless ir-basis functions of tau
  REAL(KIND = DP), ALLOCATABLE :: v(:, :)
  !! this may be not used after getting dlr
  REAL(KIND = DP), ALLOCATABLE :: dlr(:, :)
  !! change-of-basis matrix from IR basis to DLR basis
  !! dlr(i, l) = - s(l) * v(i, l)
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_f(:, :)
  !! dimensionless ir-basis functions of Matsubara freqs
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_b(:, :)
  !! dimensionless ir-basis functions of Matsubara freqs
  !
  IF (mpime == ionode_id) THEN
    READ(unit, *) tmp_str, lambda
    READ(unit, *) tmp_str, eps
    !
    rtmp = LOG(lambda) / LOG(1.0E1_DP)
    nlambda = NINT(rtmp)
    !
    rtmp = LOG(eps) / LOG(1.0E-1_DP)
    ndigit = NINT(rtmp)
    !
    ! Singular values
    READ(unit, *)
    READ(unit, *) size
    ALLOCATE(s(size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating s', 1)
    DO i = 1, size
      READ(unit, *) s(i)
    ENDDO
    !
    ! Sampling times
    READ(unit, *)
    READ(unit, *) ntau
    ALLOCATE(tau(ntau), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating tau', 1)
    DO i = 1, ntau
      READ(unit, *) tau(i)
    ENDDO
    !
    ! Basis functions on sampling times
    READ(unit, *)
    ALLOCATE(u(ntau, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating u', 1)
    DO l = 1, size
      DO t = 1, ntau
        READ(unit, *) rtmp
        u(t, l) = rtmp
      ENDDO
    ENDDO
    !
    ! Sampling frequencies (F)
    READ(unit, *)
    READ(unit, *) nfreq_f
    ALLOCATE(freq_f(nfreq_f), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating freq_f', 1)
    DO i = 1, nfreq_f
      READ(unit, *) freq_f(i)
    ENDDO
    !
    READ(unit, *)
    ALLOCATE(uhat_f(nfreq_f, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating uhat_f', 1)
    DO l = 1, size
      DO n = 1, nfreq_f
        READ(unit, *) rtmp, rtmp2
        uhat_f(n, l) = CMPLX(rtmp, rtmp2, KIND = DP)
      ENDDO
    ENDDO
    !
    ! Sampling frequencies (B)
    READ(unit, *)
    READ(unit, *) nfreq_b
    ALLOCATE(freq_b(nfreq_b), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating freq_b', 1)
    DO i = 1, nfreq_b
      READ(unit, *) freq_b(i)
    ENDDO
    !
    READ(unit, *)
    ALLOCATE(uhat_b(nfreq_b, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating uhat_b', 1)
    DO l = 1, size
      DO n = 1, nfreq_b
        READ(unit, *) rtmp, rtmp2
        uhat_b(n, l) = CMPLX(rtmp, rtmp2, KIND = DP)
      ENDDO
    ENDDO
    !
    ! Sampling poles on real frequencies
    READ(unit, *)
    READ(unit, *) nomega
    ALLOCATE(omega(nomega), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating omega', 1)
    DO i = 1, nomega
      READ(unit, *) omega(i)
    ENDDO
    !
    ! Right singular functions on sampling poles
    READ(unit, *)
    ALLOCATE(v(nomega, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating v', 1)
    ALLOCATE(dlr(nomega, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating dlr', 1)
    DO l = 1, size
      DO i = 1, nomega
        READ(unit, *) rtmp
        v(i, l) = rtmp
        dlr(i, l) = - s(l) * v(i, l)
      ENDDO
    ENDDO
  ENDIF ! mpime
  !
  CALL mp_bcast(size   , ionode_id, inter_pool_comm)
  CALL mp_bcast(ntau   , ionode_id, inter_pool_comm)
  CALL mp_bcast(nfreq_f, ionode_id, inter_pool_comm)
  CALL mp_bcast(nfreq_b, ionode_id, inter_pool_comm)
  CALL mp_bcast(nomega , ionode_id, inter_pool_comm)
  !
  IF (mpime .NE. ionode_id) THEN
    ALLOCATE(s(size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating s', 1)
    ALLOCATE(tau(ntau), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating tau', 1)
    ALLOCATE(u(ntau, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating u', 1)
    ALLOCATE(freq_f(nfreq_f), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating freq_f', 1)
    ALLOCATE(uhat_f(nfreq_f, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating uhat_f', 1)
    ALLOCATE(freq_b(nfreq_b), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating freq_b', 1)
    ALLOCATE(uhat_b(nfreq_b, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating uhat_b', 1)
    ALLOCATE(omega(nomega), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating omega', 1)
    ALLOCATE(v(nomega, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating v', 1)
    ALLOCATE(dlr(nomega, size), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_v1', 'Error allocating dlr', 1)
  ENDIF
  !
  CALL mp_bcast(s     , ionode_id, inter_pool_comm)
  CALL mp_bcast(tau   , ionode_id, inter_pool_comm)
  CALL mp_bcast(u     , ionode_id, inter_pool_comm)
  CALL mp_bcast(freq_f, ionode_id, inter_pool_comm)
  CALL mp_bcast(uhat_f, ionode_id, inter_pool_comm)
  CALL mp_bcast(freq_b, ionode_id, inter_pool_comm)
  CALL mp_bcast(uhat_b, ionode_id, inter_pool_comm)
  CALL mp_bcast(omega , ionode_id, inter_pool_comm)
  CALL mp_bcast(v     , ionode_id, inter_pool_comm)
  CALL mp_bcast(dlr   , ionode_id, inter_pool_comm)
  !
  IF ((.NOT. PRESENT(positive_only))) then
    CALL init_ir(obj, beta, lambda, eps, s, tau, freq_f, freq_b, u, uhat_f, uhat_b, omega, v, dlr, 1d-16)
  ELSE
    CALL init_ir(obj, beta, lambda, eps, s, tau, freq_f, freq_b, u, uhat_f, uhat_b, omega, v, dlr, 1d-16, positive_only)
  ENDIF
  !
  DEALLOCATE(s, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating s', 1)
  DEALLOCATE(tau, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating tau', 1)
  DEALLOCATE(u, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating u', 1)
  DEALLOCATE(freq_f, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating freq_f', 1)
  DEALLOCATE(uhat_f, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating uhat_f', 1)
  DEALLOCATE(freq_b, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating freq_b', 1)
  DEALLOCATE(uhat_b, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating uhat_b', 1)
  DEALLOCATE(omega, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating omega', 1)
  DEALLOCATE(v, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating v', 1)
  DEALLOCATE(dlr, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_v1', 'Error deallocating dlr', 1)
  !
  !-----------------------------------------------------------------------
  END FUNCTION read_v1
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  END MODULE io_sparse_ir
  !-----------------------------------------------------------------------
