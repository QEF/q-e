  !
  ! Copyright (C) 2024 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE sparse_ir
  !----------------------------------------------------------------------
  !!
  !! HM: This is adapted from sparse_ir_fortran/sparse_ir.f90,
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
  USE ep_constants,    ONLY : one, zero, cone, czero, ci, pi
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  INTERFACE decompose
    MODULE PROCEDURE decompose_z, decompose_d
  END INTERFACE decompose
  !
  INTERFACE evaluate_tau
    MODULE PROCEDURE evaluate_tau_zz, evaluate_tau_dd, evaluate_tau_dz, evaluate_tau_zd
  END INTERFACE evaluate_tau
  !
  INTERFACE evaluate_matsubara_f
    MODULE PROCEDURE evaluate_matsubara_f_zz, evaluate_matsubara_f_dz
  END INTERFACE evaluate_matsubara_f
  !
  INTERFACE evaluate_matsubara_b
    MODULE PROCEDURE evaluate_matsubara_b_zz, evaluate_matsubara_b_dz
  END INTERFACE evaluate_matsubara_b
  !
  INTERFACE fit_tau
    MODULE PROCEDURE fit_tau_zz, fit_tau_dd, fit_tau_dz, fit_tau_zd
  END INTERFACE fit_tau
  !
  INTERFACE fit_matsubara_f
    MODULE PROCEDURE fit_matsubara_f_zz, fit_matsubara_f_zd
  END INTERFACE fit_matsubara_f
  !
  INTERFACE fit_matsubara_b
    MODULE PROCEDURE fit_matsubara_b_zz, fit_matsubara_b_zd
  END INTERFACE fit_matsubara_b
  !
  INTERFACE to_dlr
    MODULE PROCEDURE to_dlr_zz, to_dlr_dd, to_dlr_dz, to_dlr_zd
  END INTERFACE to_dlr
  !
  INTERFACE evaluate_tau_from_dlr
    MODULE PROCEDURE evaluate_tau_from_dlr_zz, evaluate_tau_from_dlr_dz, &
             evaluate_tau_from_dlr_zd, evaluate_tau_from_dlr_dd
  END INTERFACE evaluate_tau_from_dlr
  !
  INTERFACE evaluate_matsubara_f_from_dlr
    MODULE PROCEDURE evaluate_matsubara_f_from_dlr_zz, evaluate_matsubara_f_from_dlr_dz
  END INTERFACE evaluate_matsubara_f_from_dlr
  !
  INTERFACE evaluate_matsubara_b_from_dlr
    MODULE PROCEDURE evaluate_matsubara_b_from_dlr_zz, evaluate_matsubara_b_from_dlr_dz
  END INTERFACE evaluate_matsubara_b_from_dlr
  !
  PUBLIC :: DecomposedMatrix_z, DecomposedMatrix_d, IR
  PUBLIC :: init_ir, set_beta, finalize_ir
  PUBLIC :: evaluate_tau, evaluate_matsubara_f, evaluate_matsubara_b
  PUBLIC :: fit_tau, fit_matsubara_f, fit_matsubara_b
  PUBLIC :: to_dlr
  PUBLIC :: evaluate_tau_from_dlr, evaluate_matsubara_f_from_dlr, evaluate_matsubara_b_from_dlr
  !
  !-----------------------------------------------------------------------
  TYPE DecomposedMatrix_z
  !-----------------------------------------------------------------------
  !!
  !! This contains a given transformation matrix and the resultants from SVD of it
  !!
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: a(:, :)
  !! original matrix A
  REAL(KIND = DP), ALLOCATABLE :: a_real(:, :)
  !! Real part of original matrix A
  REAL(KIND = DP), ALLOCATABLE :: a_imag(:, :)
  !! Imaginary part of original matrix A
  REAL(KIND = DP), ALLOCATABLE :: a_odd(:, :)
  !! having the elements from odd rows of original matrix A
  REAL(KIND = DP), ALLOCATABLE :: a_even(:, :)
  !! having the elements from even rows of original matrix
  REAL(KIND = DP), ALLOCATABLE :: inv_s_dl(:)
  !! Inverse of dimensionless singular values
  REAL(KIND = DP), ALLOCATABLE :: inv_s(:)
  !! Inverse of singular values
  COMPLEX(KIND = DP), ALLOCATABLE :: ut(:, :)
  !! transposed matrix of U which appears in A = U \Sigma V^T
  COMPLEX(KIND = DP), ALLOCATABLE :: v(:, :)
  !! transposed matrix of V^T which appears in A = U \Sigma V^T
  REAL(KIND = DP), ALLOCATABLE :: ut_real(:, :)
  !! Real part of ut
  REAL(KIND = DP), ALLOCATABLE :: ut_imag(:, :)
  !! Imaginary part of ut
  REAL(KIND = DP), ALLOCATABLE :: v_real(:, :)
  !! Real part of v
  REAL(KIND = DP), ALLOCATABLE :: v_imag(:, :)
  !! Imaginary part of v
  INTEGER :: m
  !! number of rows of A
  INTEGER :: n
  !! number of columns of A
  INTEGER :: ns
  !! number of non-zero singular values
  !
  !-----------------------------------------------------------------------
  END TYPE DecomposedMatrix_z
  !-----------------------------------------------------------------------
  ! 
  !-----------------------------------------------------------------------
  TYPE DecomposedMatrix_d
  !-----------------------------------------------------------------------
  !!
  !! This contains a given transformation matrix and the resultants from SVD of it
  !!
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: a(:, :)
  !! original matrix A
  REAL(KIND = DP), ALLOCATABLE :: a_real(:, :)
  !! Real part of original matrix A
  REAL(KIND = DP), ALLOCATABLE :: inv_s_dl(:)
  !! Inverse of dimensionless singular values
  REAL(KIND = DP), ALLOCATABLE :: inv_s(:)
  !! Inverse of singular values
  COMPLEX(KIND = DP), ALLOCATABLE :: ut(:, :)
  !! transposed matrix of U which appears in A = U \Sigma V^T
  COMPLEX(KIND = DP), ALLOCATABLE :: v(:, :)
  !! transposed matrix of V^T which appears in A = U \Sigma V^T
  REAL(KIND = DP), ALLOCATABLE :: ut_real(:, :)
  !! Real part of ut
  REAL(KIND = DP), ALLOCATABLE :: v_real(:, :)
  !! Real part of v
  INTEGER :: m
  !! number of rows of A
  INTEGER :: n
  !! number of columns of A
  INTEGER :: ns
  !! number of non-zero singular values
  !
  !-----------------------------------------------------------------------
  END TYPE DecomposedMatrix_d
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  TYPE IR
  !-----------------------------------------------------------------------
  !!
  !! This contains all the IR-basis objects, 
  !! such as sampling points, and basis functions
  !!
  !
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
  REAL(KIND = DP) :: beta
  !! inverse temperature
  REAL(KIND = DP) :: lambda
  !! lambda = 10^{nlambda}, 
  !! which determines maximum sampling point of real frequency
  REAL(KIND = DP) :: wmax
  !! maximum real frequency: wmax = lambda / beta
  REAL(KIND = DP) :: eps
  !! eps = 10^{-ndigit}
  REAL(KIND = DP) :: eps_svd
  !! This is used in the SVD fitting.
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  REAL(KIND = DP), ALLOCATABLE :: tau(:)
  !! sampling points of imaginary time
  REAL(KIND = DP), ALLOCATABLE :: x(:)
  !! This is used to get tau: tau = 5.0d-1 * beta * (x + 1.d0)
  REAL(KIND = DP), ALLOCATABLE :: omega(:)
  !! sampling points of real frequency
  REAL(KIND = DP), ALLOCATABLE :: y(:)
  !! This is used to get omega: omega = y * wmax
  INTEGER, ALLOCATABLE :: freq_f(:)
  !! integer part of sampling Matsubara freqs (Fermion)
  INTEGER, ALLOCATABLE :: freq_b(:)
  !! integer part of sampling Matsubara freqs (Boson)
  REAL(KIND = DP), ALLOCATABLE :: u_data(:,:)
  !! dimensionless IR-basis functions of tau
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_f_data(:,:)
  !! dimensionless IR-basis functions of Matsubara freqs
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_b_data(:,:)
  !! dimensionless IR-basis functions of Matsubara freqs
  REAL(KIND = DP), ALLOCATABLE :: v_data(:,:)
  !! this may be not used after getting dlr_data
  REAL(KIND = DP), ALLOCATABLE :: dlr_data(:,:)
  !! change-of-basis matrix from IR basis to DLR basis
  !! dlr_data(i, l) = - s(l) * v_data(i, l)
  TYPE(DecomposedMatrix_d) :: u
  !! stores IR-basis functions of tau 
  !! and the resultants from SVD of it
  TYPE(DecomposedMatrix_z) :: uhat_f
  !! stores IR-basis functions of Matsubara freqs (F)
  !! and the resultants from SVD of it
  TYPE(DecomposedMatrix_z) :: uhat_b
  !! stores IR-basis functions of Matsubara freqs (B)
  !! and the resultants from SVD of it
  TYPE(DecomposedMatrix_d) :: dlr
  !! stores change-of-basis matrix 
  !! and the resultants from SVD of it
  LOGICAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  !-----------------------------------------------------------------------
  END TYPE IR
  !-----------------------------------------------------------------------
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE init_ir(obj, beta, lambda, eps, s, x, freq_f, freq_b, u, uhat_f, uhat_b, y, v, dlr, eps_svd, positive_only)
  !-----------------------------------------------------------------------
  !!
  !! This routine initializes arrays related to the IR-basis objects.
  !! This routine should be called by read_ir or mk_ir_preset. 
  !! Do not call in other routines directly.
  !!
  !
  TYPE(IR), INTENT(INOUT) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  REAL(KIND = DP), INTENT(IN) :: lambda
  !! lambda = 10^{nlambda}
  REAL(KIND = DP), INTENT(IN) :: eps
  !! eps = 10^{-ndigit}
  REAL(KIND = DP), INTENT(IN) :: s(:)
  !! singular values
  REAL(KIND = DP), INTENT(IN) :: x(:)
  !! This is used to get tau: tau = 5.0d-1 * beta * (x + 1.d0)
  REAL(KIND = DP), INTENT(IN) :: y(:)
  !! This is used to get omega: omega = y * wmax
  REAL(KIND = DP), INTENT(IN) :: eps_svd
  !! This is used in the SVD fitting.
  REAL(KIND = DP), INTENT(IN) :: u(:,:)
  !! IR-basis functions of tau for beta
  COMPLEX(KIND = DP), INTENT(IN) :: uhat_f(:, :)
  !! stores IR-basis functions of Matsubara freqs (F) for beta
  COMPLEX(KIND = DP), INTENT(IN) :: uhat_b(:, :)
  !! stores IR-basis functions of Matsubara freqs (B) for beta
  REAL(KIND = DP), INTENT(IN) :: v(:, :)
  !! this may be not used
  REAL(KIND = DP), INTENT(IN) :: dlr(:, :)
  !! change-of-basis matrix from IR basis to DLR basis for beta
  INTEGER, INTENT(IN) :: freq_f(:)
  !! integer part of sampling Matsubara freqs (Fermion)
  INTEGER, INTENT(IN) :: freq_b(:)
  !! integer part of sampling Matsubara freqs (Boson)
  LOGICAL, INTENT(IN), OPTIONAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  INTEGER :: ierr
  !! Error status
  !
  IF ((.NOT. PRESENT(positive_only)) .OR. (.NOT. positive_only)) THEN
    obj%positive_only = .false.
  ELSE
    obj%positive_only = .true.
  ENDIF
  !
  IF (ALLOCATED(obj%x)) THEN
    CALL errore('init_ir', 'IR%x is already allocated. &
               &You should call finalize_ir before recalling init_ir.', 1)
  ENDIF
  !
  obj%size = SIZE(s)
  obj%ntau = SIZE(x)
  IF (.NOT. obj%positive_only) THEN
    obj%nfreq_f = SIZE(freq_f)
    obj%nfreq_b = SIZE(freq_b)
  ELSE
    obj%nfreq_f = SIZE(freq_f) / 2
    obj%nfreq_b = (SIZE(freq_b)-1) / 2 + 1
  ENDIF
  obj%lambda = lambda
  obj%eps = eps
  obj%eps_svd = eps_svd
  obj%nomega = SIZE(y)
  !
  ALLOCATE(obj%x(obj%ntau), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%x', 1)
  obj%x = x
  !
  ALLOCATE(obj%tau(obj%ntau), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%tau', 1)
  !
  ALLOCATE(obj%s(obj%size), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%s', 1)
  obj%s = SQRT(5.0d-1*obj%lambda) * s
  !
  ALLOCATE(obj%freq_f(obj%nfreq_f), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%freq_f', 1)
  !
  ALLOCATE(obj%freq_b(obj%nfreq_b), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%freq_b', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    obj%freq_f = freq_f
    obj%freq_b = freq_b
  ELSE
    obj%freq_f(1:obj%nfreq_f) = freq_f((SIZE(freq_f) / 2 + 1):SIZE(freq_f))
    obj%freq_b(1:obj%nfreq_b) = freq_b(((SIZE(freq_b) + 1) / 2):SIZE(freq_b))
  ENDIF
  !
  ALLOCATE(obj%u_data(obj%ntau, obj%size), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%u_data', 1)
  obj%u_data = u
  !
  ALLOCATE(obj%uhat_f_data(obj%nfreq_f, obj%size), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%uhat_f_data', 1)
  !
  ALLOCATE(obj%uhat_b_data(obj%nfreq_b, obj%size), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%uhat_b_data', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    obj%uhat_f_data = uhat_f
    obj%uhat_b_data = uhat_b
  ELSE
    obj%uhat_f_data(1:obj%nfreq_f, :) = uhat_f((SIZE(freq_f) / 2 + 1):SIZE(freq_f), :)
    obj%uhat_b_data(1:obj%nfreq_b, :) = uhat_b(((SIZE(freq_b) + 1) / 2):SIZE(freq_b), :)
  ENDIF
  !
  ALLOCATE(obj%y(obj%nomega), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%y', 1)
  obj%y = y
  !
  ALLOCATE(obj%omega(obj%nomega), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%omega', 1)
  !
  ALLOCATE(obj%v_data(obj%nomega, obj%size), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%v_data', 1)
  obj%v_data = v
  !
  ALLOCATE(obj%dlr_data(obj%size, obj%nomega), STAT = ierr)
  IF (ierr /= 0) CALL errore('init_ir', 'Error allocating IR%dlr_data', 1)
  obj%dlr_data = TRANSPOSE(dlr)
  !
  obj%u = decompose(obj%u_data, obj%eps_svd, .false.)
  IF (.NOT. obj%positive_only) THEN
    obj%uhat_f = decompose(obj%uhat_f_data, obj%eps_svd, .false.)
    obj%uhat_b = decompose(obj%uhat_b_data, obj%eps_svd, .false.)
  ELSE
    obj%uhat_f = split_decompose(obj%uhat_f_data, .false., obj%eps_svd, .false.)
    obj%uhat_b = split_decompose(obj%uhat_b_data, .true., obj%eps_svd, .false.)
  ENDIF
  !
  obj%dlr = decompose(obj%dlr_data, obj%eps_svd, .true.)
  !
  ! Here we define basis sets for the input value of beta.
  CALL set_beta(obj, beta)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE init_ir
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_beta(obj, beta)
  !-----------------------------------------------------------------------
  !!
  !! This routine update IR-basis objects for a given beta.
  !! This is called in init_ir at once.
  !!
  !
  TYPE(IR), INTENT(INOUT) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  !
  obj%beta = beta
  obj%wmax = obj%lambda / beta
  !
  obj%tau = 5.0d-1 * beta * (obj%x + 1.d0)
  obj%omega = obj%y * obj%wmax
  !
  obj%u%a_real(:, :) = SQRT(2.0D0/beta)*obj%u_data(:, :)
  obj%uhat_f%a(:, :) = SQRT(beta) * obj%uhat_f_data(:, :)
  obj%uhat_b%a(:, :) = SQRT(beta) * obj%uhat_b_data(:, :)
  obj%dlr%a(:, :) = SQRT(5.0d-1*beta)*obj%dlr_data(:, :)
  !
  obj%u%a = CMPLX(obj%u%a_real, zero, KIND = DP)
  obj%uhat_f%a_real = REAL(obj%uhat_f%a, KIND = DP)
  obj%uhat_f%a_imag = AIMAG(obj%uhat_f%a)
  obj%uhat_b%a_real = REAL(obj%uhat_b%a, KIND = DP)
  obj%uhat_b%a_imag = AIMAG(obj%uhat_b%a)
  obj%uhat_f%a_odd(:, :) = obj%uhat_f%a_imag(:, 1:(obj%uhat_f%n - 1):2)
  obj%uhat_f%a_even(:, :) = obj%uhat_f%a_real(:, 2:obj%uhat_f%n:2)
  obj%uhat_b%a_odd(:, :) = obj%uhat_b%a_real(:, 1:(obj%uhat_b%n - 1):2)
  obj%uhat_b%a_even(:, :) = obj%uhat_b%a_imag(:, 2:obj%uhat_b%n:2)
  !
  obj%u%inv_s(:) = SQRT(5.0d-1*beta) * obj%u%inv_s_dl(:)
  obj%uhat_f%inv_s(:) = (1.0D0 / SQRT(beta)) * obj%uhat_f%inv_s_dl(:)
  obj%uhat_b%inv_s(:) = (1.0D0 / SQRT(beta)) * obj%uhat_b%inv_s_dl(:)
  obj%dlr%inv_s(:) = SQRT(2.0D0 / beta) * obj%dlr%inv_s_dl(:)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE set_beta
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE finalize_ir(obj)
  !-----------------------------------------------------------------------
  !!
  !! This routine deallocates IR-basis objects contained in obj
  !!
  !
  TYPE(IR) :: obj
  !! contains all the IR-basis objects
  INTEGER :: ierr
  !! Error status
  !
  DEALLOCATE(obj%x, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%x', 1)
  DEALLOCATE(obj%tau, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%tau', 1)
  DEALLOCATE(obj%s, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%s', 1)
  DEALLOCATE(obj%freq_f, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%freq_f', 1)
  DEALLOCATE(obj%freq_b, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%freq_b', 1)
  DEALLOCATE(obj%u_data, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%u_data', 1)
  DEALLOCATE(obj%uhat_f_data, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%uhat_f_data', 1)
  DEALLOCATE(obj%uhat_b_data, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%uhat_b_data', 1)
  DEALLOCATE(obj%y, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%y', 1)
  DEALLOCATE(obj%omega, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%omega', 1)
  DEALLOCATE(obj%v_data, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%v_data', 1)
  DEALLOCATE(obj%dlr_data, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_ir', 'Error deallocating IR%dlr_data', 1)
  !
  CALL finalize_dmat_d(obj%u)
  CALL finalize_dmat_z(obj%uhat_f)
  CALL finalize_dmat_z(obj%uhat_b)
  CALL finalize_dmat_d(obj%dlr)
  !-----------------------------------------------------------------------
  END SUBROUTINE finalize_ir
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE finalize_dmat_z(dmat)
  !-----------------------------------------------------------------------
  !!
  !! This routine deallocates arrays contained in dmat
  !!
  !
  TYPE(DecomposedMatrix_z) :: dmat
  !! contains a given transformation matrix and the resultants from SVD of it
  INTEGER :: ierr
  !! Error status
  !
  DEALLOCATE(dmat%a, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a', 1)
  DEALLOCATE(dmat%a_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a_real', 1)
  DEALLOCATE(dmat%a_imag, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a_imag', 1)
  DEALLOCATE(dmat%a_odd, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a_odd', 1)
  DEALLOCATE(dmat%a_even, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a_even', 1)
  DEALLOCATE(dmat%inv_s, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%inv_s', 1)
  DEALLOCATE(dmat%inv_s_dl, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%inv_s_dl', 1)
  DEALLOCATE(dmat%ut, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%ut', 1)
  DEALLOCATE(dmat%v, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%v', 1)
  DEALLOCATE(dmat%ut_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%ut_real', 1)
  DEALLOCATE(dmat%v_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%v_real', 1)
  DEALLOCATE(dmat%ut_imag, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%ut_imag', 1)
  DEALLOCATE(dmat%v_imag, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%v_imag', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE finalize_dmat_z
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE finalize_dmat_d(dmat)
  !-----------------------------------------------------------------------
  !!
  !! This routine deallocates arrays contained in dmat
  !!
  !
  TYPE(DecomposedMatrix_d) :: dmat
  !! contains a given transformation matrix and the resultants from SVD of it
  INTEGER :: ierr
  !! Error status
  !
  DEALLOCATE(dmat%a, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a', 1)
  DEALLOCATE(dmat%a_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%a_real', 1)
  DEALLOCATE(dmat%inv_s, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%inv_s', 1)
  DEALLOCATE(dmat%inv_s_dl, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%inv_s_dl', 1)
  DEALLOCATE(dmat%ut, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%ut', 1)
  DEALLOCATE(dmat%v, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%v', 1)
  DEALLOCATE(dmat%ut_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%ut_real', 1)
  DEALLOCATE(dmat%v_real, STAT = ierr)
  IF (ierr /= 0) CALL errore('finalize_dmat_z', 'Error deallocating dmat%v_real', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE finalize_dmat_d
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION decompose_z(a, eps, ill_conditioned) result(dmat)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs SVD of complex-valued matrix a. 
  !! Singular values smaller than eps * the largest one are dropped.
  !!
  !
  COMPLEX(KIND = DP), INTENT(IN) :: a(:, :)
  !! original matrix a
  REAL(KIND = DP), INTENT(IN) :: eps
  !! threshold for singular values
  LOGICAL, INTENT(IN) :: ill_conditioned
  !! If true, use ZGESVD instead of ZGESDD
  !
  INTEGER :: i
  !! Counter
  INTEGER :: info
  !! Error status
  INTEGER :: lda
  !! lda = m
  INTEGER :: ldu
  !! ldu = m
  INTEGER :: ldvt
  !! ldvt = n
  INTEGER :: lwork
  !! lwork = mn*mn + 3*mn (ill_conditioned = false) or 2*mn + m + n (true)
  INTEGER :: m
  !! number of rows of a
  INTEGER :: n
  !! number of columns of a
  INTEGER :: mn
  !! m times n
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: ierr
  !! Error status
  COMPLEX(KIND = DP), ALLOCATABLE :: a_copy(:, :)
  !! copy of a
  COMPLEX(KIND = DP), ALLOCATABLE :: u(:, :)
  !! U in A=U \Sigma V^T 
  COMPLEX(KIND = DP), ALLOCATABLE :: vt(:, :)
  !! V^T in A=U \Sigma V^T
  COMPLEX(KIND = DP), ALLOCATABLE :: work(:)
  !! arrays for ZGESDD
  REAL(KIND = DP), ALLOCATABLE :: rwork(:)
  !! arrays for ZGESDD
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  INTEGER, ALLOCATABLE :: iwork(:)
  !! arrays for ZGESDD
  TYPE(DecomposedMatrix_z)::dmat
  !! to contain matrix a and the resultants from SVD of it
  !
  IF (ALLOCATED(dmat%a)) THEN
    CALL errore('decompose_z', 'dmat%a is already allocated.  &
               &You should call finalize_dmat before recalling decompose.', 1)
  ENDIF
  !
  m = SIZE(a, 1)
  n = SIZE(a, 2)
  mn = min(m, n)
  lda = m
  ldu = m
  ldvt = n
  IF (.NOT. ill_conditioned) THEN
    lwork = mn*mn + 3*mn
    ALLOCATE(rwork((5*mn+7)*mn), STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating rwork', 1)
    ALLOCATE(iwork(8*mn), STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating iwork', 1)
  ELSE
    lwork = 2*mn + m + n
    ALLOCATE(rwork(5*n), STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating rwork', 1)
  ENDIF
  ALLOCATE(work(lwork), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating work', 1)
  ALLOCATE(a_copy(m,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating a_copy', 1)
  ALLOCATE(s(m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating s', 1)
  ALLOCATE(u(ldu,m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating u', 1)
  ALLOCATE(vt(ldvt,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating vt', 1)
  !
  a_copy(1:m, 1:n) = a(1:m, 1:n)
  IF (.NOT. ill_conditioned) THEN
    lwork = mn*mn + 3*mn
    call ZGESDD('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose_z', 'Failure in ZGESDD', info)
    ENDIF
  ELSE
    lwork = 2*mn + m + n
    call ZGESVD('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose_z', 'Failure in ZGESVD', info)
    ENDIF
  ENDIF
  !
  ! Number of relevant singular values s(i)/s(1) >= eps
  ns = 0
  DO i = 1, mn
    IF (s(i)/s(1) < eps) THEN
      exit
    ENDIF
    ns = ns + 1
  ENDDO
  !
  ALLOCATE(dmat%a(m, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%a', 1)
  ALLOCATE(dmat%a_real(m, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%a_real', 1)
  ALLOCATE(dmat%a_imag(m, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%a_imag', 1)
  ALLOCATE(dmat%inv_s_dl(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%inv_s_dl', 1)
  ALLOCATE(dmat%inv_s(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%inv_s', 1)
  ALLOCATE(dmat%ut(ns, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%ut', 1)
  ALLOCATE(dmat%ut_real(ns, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%ut_real', 1)
  ALLOCATE(dmat%ut_imag(ns, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%ut_imag', 1)
  ALLOCATE(dmat%v(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%v', 1)
  ALLOCATE(dmat%v_real(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%v_real', 1)
  ALLOCATE(dmat%v_imag(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%v_imag', 1)
  ALLOCATE(dmat%a_odd(m, n/2), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%a_odd', 1)
  ALLOCATE(dmat%a_even(m, n/2), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error allocating dmat%a_even', 1)
  !
  ! dmat%a temporarily stores the same data of input a
  dmat%a = a
  dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
  ! inv_s temporarily stores the same data of inv_s_dl
  dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
  dmat%ut(1:ns, 1:m) = CONJG(TRANSPOSE(u(1:m, 1:ns)))
  dmat%v(1:n, 1:ns) = CONJG(TRANSPOSE(vt(1:ns, 1:n)))
  dmat%m = SIZE(a, 1)
  dmat%n = SIZE(a, 2)
  dmat%ns = ns
  dmat%a_odd = zero
  dmat%a_even = zero
  !
  dmat%a_real = REAL(dmat%a, KIND = DP)
  dmat%a_imag = AIMAG(dmat%a)
  dmat%ut_real = REAL(dmat%ut, KIND = DP)
  dmat%ut_imag = AIMAG(dmat%ut)
  dmat%v_real = REAL(dmat%v, KIND = DP)
  dmat%v_imag = AIMAG(dmat%v)
  !
  DEALLOCATE(work, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating work', 1)
  DEALLOCATE(a_copy, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating a_copy', 1)
  DEALLOCATE(s, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating s', 1)
  DEALLOCATE(u, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating u', 1)
  DEALLOCATE(vt, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating vt', 1)
  DEALLOCATE(rwork, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_z', 'Error deallocating rwork', 1)
  !
  IF (.NOT. ill_conditioned) THEN
    DEALLOCATE(iwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose', 'Error deallocating iwork', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END FUNCTION decompose_z
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION decompose_d(a, eps, ill_conditioned) result(dmat)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs SVD of real-valued matrix a. 
  !! Singular values smaller than eps * the largest one are dropped.
  !!
  !
  REAL(KIND = DP), INTENT(IN) :: a(:, :)
  !! original matrix a
  REAL(KIND = DP), INTENT(IN) :: eps
  !! threshold for singular values
  LOGICAL, INTENT(IN) :: ill_conditioned
  !! If true, use DGESVD instead of DGESDD
  !
  INTEGER :: i
  !! Counter
  INTEGER :: info
  !! Error status
  INTEGER :: lda
  !! lda = m
  INTEGER :: ldu
  !! ldu = m
  INTEGER :: ldvt
  !! ldvt = n
  INTEGER :: lwork
  !! lwork = 4*mn*mn + 7*mn (ill_conditioned = false) or 5*mn + m + n (true)
  INTEGER :: m
  !! number of rows of a
  INTEGER :: n
  !! number of columns of a
  INTEGER :: mn
  !! m times n
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP), ALLOCATABLE :: a_copy(:, :)
  !! copy of a
  REAL(KIND = DP), ALLOCATABLE :: u(:, :)
  !! U in A=U \Sigma V^T 
  REAL(KIND = DP), ALLOCATABLE :: vt(:, :)
  !! V^T in A=U \Sigma V^T
  REAL(KIND = DP), ALLOCATABLE :: work(:)
  !! arrays for ZGESDD
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  INTEGER, ALLOCATABLE :: iwork(:)
  !! arrays for DGESDD
  TYPE(DecomposedMatrix_d)::dmat
  !! to contain matrix a and the resultants from SVD of it
  !
  IF (ALLOCATED(dmat%a)) THEN
    CALL errore('decompose_d', 'dmat%a is already allocated.  &
               &You should call finalize_dmat before recalling decompose.', 1)
  ENDIF
  !
  m = SIZE(a, 1)
  n = SIZE(a, 2)
  mn = min(m, n)
  lda = m
  ldu = m
  ldvt = n
  IF (.NOT. ill_conditioned) THEN
    lwork = 4*mn*mn + 7*mn
    ALLOCATE(iwork(8*mn), STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating iwork', 1)
  ELSE
    lwork = 5*mn + m + n
  ENDIF
  ALLOCATE(work(lwork), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating work', 1)
  ALLOCATE(a_copy(m,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating a_copy', 1)
  ALLOCATE(s(m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating s', 1)
  ALLOCATE(u(ldu,m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating u', 1)
  ALLOCATE(vt(ldvt,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating vt', 1)
  !
  a_copy(1:m, 1:n) = a(1:m, 1:n)
  IF (.NOT. ill_conditioned) THEN
    lwork = 4*mn*mn + 7*mn
    call DGESDD('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose_d', 'Failure in DGESDD', info)
    ENDIF
  ELSE
    lwork = 5*mn + m + n
    call DGESVD('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose_d', 'Failure in DGESVD', info)
    ENDIF
  ENDIF
  !
  ! Number of relevant singular values s(i)/s(1) >= eps
  ns = 0
  DO i = 1, mn
    IF (s(i)/s(1) < eps) THEN
      exit
    ENDIF
    ns = ns + 1
  ENDDO
  !
  ALLOCATE(dmat%a(m, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%a', 1)
  ALLOCATE(dmat%a_real(m, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%a_real', 1)
  ALLOCATE(dmat%inv_s_dl(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%inv_s_dl', 1)
  ALLOCATE(dmat%inv_s(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%inv_s', 1)
  ALLOCATE(dmat%ut(ns, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%ut', 1)
  ALLOCATE(dmat%ut_real(ns, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%ut_real', 1)
  ALLOCATE(dmat%v(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%v', 1)
  ALLOCATE(dmat%v_real(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating dmat%v_real', 1)
  !
  ! dmat%a temporarily stores the same data of input a
  dmat%a_real = a
  dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
  ! inv_s temporarily stores the same data of inv_s_dl
  dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
  dmat%ut_real(1:ns, 1:m) = TRANSPOSE(u(1:m, 1:ns))
  dmat%v_real(1:n, 1:ns) = TRANSPOSE(vt(1:ns, 1:n))
  dmat%m = SIZE(a, 1)
  dmat%n = SIZE(a, 2)
  dmat%ns = ns 
  !
  dmat%a = dmat%a_real
  dmat%ut = dmat%ut_real
  dmat%v = dmat%v_real
  !
  DEALLOCATE(work, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating work', 1)
  DEALLOCATE(a_copy, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating a_copy', 1)
  DEALLOCATE(s, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating s', 1)
  DEALLOCATE(u, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating u', 1)
  DEALLOCATE(vt, STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating vt', 1)
  !
  IF (.NOT. ill_conditioned) THEN
    DEALLOCATE(iwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_d', 'Error deallocating iwork', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END FUNCTION decompose_d
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION split_decompose(a, has_zero, eps, ill_conditioned) RESULT(dmat)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs "split SVD" of uhat_f or uhat_b.
  !! Singular values smaller than eps * the largest one are dropped.
  !!
  !! Stores a matrix "A" together with its "split SVD" form::
  !!
  !!     A == u * s @ vT
  !!
  !! where "vT" is a real matrix and "u" is a complex matrix.  The "split" SVD
  !! form differs from the SVD in that the least squares fit has to be
  !! constructed as follows:
  !!
  !!     fit(A, x) == v / s @ REAL(CONJG(TRANSPOSE(u)) * x)
  !!
  !! This again allows for fast and accurate least squares fits.
  !!
  !! This subroutine is called if positive_only = true.
  !! If s, u, and vT are obtained in this subroutine, the imaginary
  !! part of (CONJG(TRANSPOSE(u)) * x) returns a non-zero
  !! value indeed, but that value should be ignored. 
  !! Only the real part will be used in fit_matsubara_{f,b}
  !! if positive_only = true.
  !!
  !
  COMPLEX(KIND = DP), INTENT(IN) :: a(:, :)
  !! original matrix a
  LOGICAL, INTENT(IN) :: has_zero
  !! if a is uhat_b, has_zero should be true
  REAL(KIND = DP), INTENT(IN) :: eps
  !! threshold for singular values
  LOGICAL, INTENT(IN) :: ill_conditioned
  !! If true, use DGESVD instead of DGESDD
  !
  INTEGER :: i
  !! Counter
  INTEGER :: info
  !! Error status
  INTEGER :: lda
  !! lda = m
  INTEGER :: ldu
  !! ldu = m
  INTEGER :: ldvt
  !! ldvt = n
  INTEGER :: lwork
  !! lwork = 4*mn*mn + 7*mn (ill_conditioned = false) or 5*mn + m + n (true)
  INTEGER :: m_half
  !! number of rows of a
  INTEGER :: m
  !! m = 2 * m_half - 1 (has_zero = true) or 2 * m_half (false)
  INTEGER :: n
  !! number of columns of a
  INTEGER :: mn
  !! m times n
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP), ALLOCATABLE :: a_copy(:, :)
  !! copy of a
  REAL(KIND = DP), ALLOCATABLE :: u(:, :)
  !! U in A=U \Sigma V^T 
  REAL(KIND = DP), ALLOCATABLE :: vt(:, :)
  !! V^T in A=U \Sigma V^T
  REAL(KIND = DP), ALLOCATABLE :: work(:)
  !! arrays for ZGESDD
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  COMPLEX(KIND = DP), ALLOCATABLE :: u_copy(:, :)
  !! copy of u
  INTEGER, ALLOCATABLE :: iwork(:)
  !! arrays for DGESDD
  TYPE(DecomposedMatrix_z)::dmat
  !! to contain matrix a and the resultants from SVD of it
  !
  IF (ALLOCATED(dmat%a)) THEN
    CALL errore('split_decompose', 'dmat%a is already allocated.  &
               &You should call finalize_dmat before recalling split_decompose.', 1)
  ENDIF
  !
  m_half = SIZE(a, 1)
  IF (has_zero) THEN
    m = 2 * m_half - 1
  ELSE
    m = 2 * m_half
  ENDIF
  !
  n = SIZE(a, 2)
  mn = min(m, n)
  lda = m
  ldu = m
  ldvt = n
  IF (.NOT. ill_conditioned) THEN
    lwork = 4*mn*mn + 7*mn
    ALLOCATE(iwork(8*mn), STAT = ierr)
    IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating iwork', 1)
  ELSE
    lwork = 5*mn + m + n
  ENDIF
  ALLOCATE(work(lwork), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating work', 1)
  ALLOCATE(a_copy(m,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating a_copy', 1)
  ALLOCATE(s(m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating s', 1)
  ALLOCATE(u(ldu,m), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating u', 1)
  ALLOCATE(vt(ldvt,n), STAT = ierr)
  IF (ierr /= 0) CALL errore('decompose_d', 'Error allocating vt', 1)
  !
  a_copy(1:m_half, 1:n) = REAL(a(1:m_half, 1:n))
  !
  IF (has_zero) THEN
    a_copy(m_half+1:m, 1:n) = AIMAG(a(2:m_half, 1:n))
  ELSE
    a_copy(m_half+1:m, 1:n) = AIMAG(a(1:m_half, 1:n))
  ENDIF
  !
  IF (.NOT. ill_conditioned) THEN
    lwork = 4*mn*mn + 7*mn
    call DGESDD('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose', 'Failure in DGESDD', info)
    ENDIF
  ELSE
    lwork = 5*mn + m + n
    call DGESVD('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    IF (info /= 0) THEN
      CALL errore('decompose', 'Failure in DGESVD', info)
    ENDIF
  ENDIF
  !
  ! Number of relevant singular values s(i)/s(1) >= eps
  ns = 0
  DO i = 1, mn
    IF (s(i)/s(1) < eps) THEN
      exit
    ENDIF
    ns = ns + 1
  ENDDO
  !
  ALLOCATE(u_copy(m_half, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating u_copy', 1)
  !
  IF (has_zero) THEN
    u_copy(1, 1:ns) = CMPLX(u(1, 1:ns), 0.0D0, KIND = DP)
    u_copy(2:m_half, 1:ns) = CMPLX(u(2:m_half, 1:ns), u(m_half+1:m, 1:ns), KIND = DP)
  ELSE
    u_copy(1:m_half, 1:ns) = CMPLX(u(1:m_half, 1:ns), u(m_half+1:m, 1:ns), KIND = DP)
  ENDIF
  !
  ALLOCATE(dmat%a(m_half, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating a', 1)
  ALLOCATE(dmat%a_real(m_half, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating a_real', 1)
  ALLOCATE(dmat%a_imag(m_half, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating a_imag', 1)
  ALLOCATE(dmat%inv_s_dl(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating inv_s_dl', 1)
  ALLOCATE(dmat%inv_s(ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating inv_s', 1)
  ALLOCATE(dmat%ut(ns, m_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating ut', 1)
  ALLOCATE(dmat%ut_real(ns, m_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating ut_real', 1)
  ALLOCATE(dmat%ut_imag(ns, m_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating ut_imag', 1)
  ALLOCATE(dmat%v(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating v', 1)
  ALLOCATE(dmat%v_real(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating v_real', 1)
  ALLOCATE(dmat%v_imag(n, ns), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating v_imag', 1)
  ALLOCATE(dmat%a_odd(m_half, n/2), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating a_odd', 1)
  ALLOCATE(dmat%a_even(m_half, n/2), STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error allocating a_even', 1)
  !
  ! dmat%a temporarily stores the same data of input a
  dmat%a = a
  dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
  ! inv_s temporarily stores the same data of inv_s_dl
  dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
  dmat%ut(1:ns, 1:m_half) = CONJG(TRANSPOSE(u_copy(1:m_half, 1:ns)))
  dmat%v_real(1:n, 1:ns) = TRANSPOSE(vt(1:ns, 1:n))
  dmat%m = SIZE(a, 1)
  dmat%n = SIZE(a, 2)
  dmat%ns = ns
  dmat%a_odd = zero
  dmat%a_even = zero
  !
  dmat%a_real = REAL(dmat%a, KIND = DP)
  dmat%a_imag = AIMAG(dmat%a)
  dmat%ut_real = REAL(dmat%ut, KIND = DP)
  dmat%ut_imag = AIMAG(dmat%ut)
  dmat%v_imag = zero
  dmat%v = CMPLX(dmat%v_real, zero, KIND = DP)
  !
  DEALLOCATE(work, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating work', 1)
  DEALLOCATE(a_copy, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating a_copy', 1)
  DEALLOCATE(s, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating s', 1)
  DEALLOCATE(u, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating u', 1)
  DEALLOCATE(vt, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating vt', 1)
  DEALLOCATE(u_copy, STAT = ierr)
  IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating u_copy', 1)
  !
  IF (.NOT. ill_conditioned) THEN
    DEALLOCATE(iwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('split_decompose', 'Error deallocating iwork', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END FUNCTION split_decompose
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto sampling points of imaginary time
  !! from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_zz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('evaluate_tau_zz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('evaluate_tau_zz', 'wrong number of columns of output array.', 1)
  res(:, :) = czero
  !
  IF (.NOT. obj%positive_only) THEN
    !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
    !      l2, obj%u%a, m, czero, res(:, :), l2)
    !
    ALLOCATE(arr_tmp(l1, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error allocating res_r', 1)
    ALLOCATE(res_i(l2, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error allocating res_i', 1)
    ! calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
          l2, obj%u%a_real, m, zero, res_r(:, :), l2)
    !
    ! calculate the imaginary part
    arr_tmp = AIMAG(arr)
    CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
          l2, obj%u%a_real, m, zero, res_i(:, :), l2)
    res = CMPLX(res_r, res_i, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error deallocating res_r', 1)
    DEALLOCATE(res_i, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error deallocating res_i', 1)
  ELSE
    ALLOCATE(arr_tmp(l1, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error allocating res_r', 1)
    ! only calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
          l2, obj%u%a_real, m, zero, res_r(:, :), l2)
    res = CMPLX(res_r, zero, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_tau_zz', 'Error deallocating res_r', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_dd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto sampling points of imaginary time
  !! from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! If positive_only is false, evaluate_tau_zz 
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! both the input and output arrays are real arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from IR expansion coefficients
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_dd', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('evaluate_tau_dd', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('evaluate_tau_dd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_dd', 'input and output arrays should be complex arrays.', 1)
  res(:, :) = zero
  !
  CALL DGEMM('n', 't', l1, m, n, one, arr(:,:), &
        l2, obj%u%a_real, m, zero, res(:, :), l2)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_dd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_dz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto sampling points of imaginary time
  !! from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! If positive_only is false, evaluate_tau_zz 
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: res_tmp(:, :)
  !! dummy for res
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_dz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('evaluate_tau_dz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('evaluate_tau_dz', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_dz', 'input array should be a complex array.', 1)
  res(:, :) = czero
  !
  ALLOCATE(res_tmp(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_tau_dz', 'Error allocating res_tmp', 1)
  CALL DGEMM('n', 't', l1, m, n, one, arr(:,:), &
        l2, obj%u%a_real, m, zero, res_tmp(:, :), l2)
  res = CMPLX(res_tmp, zero, KIND = DP)
  DEALLOCATE(res_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_tau_dz', 'Error deallocating res_tmp', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_zd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto sampling points of imaginary time
  !! from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! If positive_only is false, evaluate_tau_zz 
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_zd', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('evaluate_tau_zd', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('evaluate_tau_zd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_zd', 'output array should be a complex array.', 1)
  res(:, :) = zero
  !
  ALLOCATE(arr_tmp(l1, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_tau_zd', 'Error allocating arr_tmp', 1)
  arr_tmp = REAL(arr, KIND = DP)
  CALL DGEMM('n', 't', l1, m, n, one, arr_tmp(:,:), &
        l2, obj%u%a_real, m, zero, res(:, :), l2)
  DEALLOCATE(arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_tau_zd', 'Error deallocating arr_tmp', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_f_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a fermionic function onto sampling 
  !! Matsubara frequencies from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! A_{mn} is the n-th basis function of m,
  !! where m is a sampling Matsubara frequency.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a fermionic function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: arr_half(:, :)
  !! dummy for arr
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%uhat_f%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%uhat_f%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: n_half
  !! n_half = n / 2
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_f_zz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%uhat_f%n) &
    CALL errore('evaluate_matsubara_f_zz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%uhat_f%m) &
    CALL errore('evaluate_matsubara_f_zz', 'wrong number of columns of output array.', 1)
  res(:, :) = czero
  !
  !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
  !      l2, obj%uhat_f%a, m, czero, res(:, :), l2)
  !
  n_half = n / 2
  !
  ALLOCATE(res_r(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error allocating res_r', 1)
  ALLOCATE(res_i(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error allocating res_i', 1)
  ALLOCATE(arr_half(l1, n_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error allocating arr_half', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ! take partial summations over even l
    arr_half(:,:) = REAL(arr(:,2:n:2))
    ! calculate the real part
    ! Re(arr) * Re(uhat_f)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
    !
    arr_half(:,:) = AIMAG(arr(:,2:n:2))
    ! calculate the imaginary part
    ! Im(arr) * Re(uhat_f)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_even, m, zero, res_i(:, :), l2)
    !
    res = CMPLX(res_r, res_i, KIND = DP)
    !
    ! take partial summations over odd l
    arr_half(:,:) = AIMAG(arr(:,1:(n-1):2))
    ! calculate the real part
    ! Im(arr) * Im(uhat_f)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_odd, m, zero, res_r(:, :), l2)
    !
    arr_half(:,:) = REAL(arr(:,1:(n-1):2))
    ! calculate the imaginary part
    ! Re(arr) * Im(uhat_f)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
    !
    res = res + CMPLX(-res_r, res_i, KIND = DP)
  ELSE
    ! take partial summations over even l
    arr_half(:,:) = REAL(arr(:,2:n:2))
    ! calculate the real part
    ! Re(arr) * Re(uhat_f)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
    !
    ! take partial summations over odd l
    arr_half(:,:) = REAL(arr(:,1:(n-1):2))
    ! calculate the imaginary part
    ! Re(arr) * Im(uhat_f)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
    !
    res = CMPLX(res_r, res_i, KIND = DP)
  ENDIF
  DEALLOCATE(res_r, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error deallocating res_r', 1)
  DEALLOCATE(res_i, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error deallocating res_i', 1)
  DEALLOCATE(arr_half, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_zz', 'Error deallocating arr_half', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_f_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_f_dz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a fermionic function onto sampling 
  !! Matsubara frequencies from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! A_{mn} is the n-th basis function of m,
  !! where m is a sampling Matsubara frequency.
  !!
  !! If positive_only is false, evaluate_matsubara_f_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a fermionic function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: arr_half(:, :)
  !! dummy for arr
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%uhat_f%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%uhat_f%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: n_half
  !! n_half = n / 2
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_f_dz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%uhat_f%n) &
    CALL errore('evaluate_matsubara_f_dz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%uhat_f%m) &
    CALL errore('evaluate_matsubara_f_dz', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_matsubara_f_dz', 'input array should be a complex array.', 1)
  res(:, :) = czero
  !
  !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
  !      l2, obj%uhat_f%a, m, czero, res(:, :), l2)
  !
  n_half = n / 2
  !
  ALLOCATE(res_r(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error allocating res_r', 1)
  ALLOCATE(res_i(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error allocating res_i', 1)
  ALLOCATE(arr_half(l1, n_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error allocating arr_half', 1)
  !
  ! take partial summations over even l
  arr_half(:,:) = arr(:,2:n:2)
  ! calculate the real part
  ! Re(arr) * Re(uhat_f)**T
  res_r(:, :) = zero
  CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
        l2, obj%uhat_f%a_even, m, zero, res_r(:, :), l2)
  !
  ! take partial summations over odd l
  arr_half(:,:) = arr(:,1:(n-1):2)
  ! calculate the imaginary part
  ! Re(arr) * Im(uhat_f)**T
  res_i(:, :) = zero
  CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
        l2, obj%uhat_f%a_odd, m, zero, res_i(:, :), l2)
  !
  res = CMPLX(res_r, res_i, KIND = DP)
  DEALLOCATE(res_r, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error deallocating res_r', 1)
  DEALLOCATE(res_i, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error deallocating res_i', 1)
  DEALLOCATE(arr_half, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_f_dz', 'Error deallocating arr_half', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_f_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_b_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a bosonic function onto sampling 
  !! Matsubara frequencies from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! A_{mn} is the n-th basis function of m,
  !! where m is a sampling Matsubara frequency.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a bosonic function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: arr_half(:, :)
  !! dummy for arr
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%uhat_b%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%uhat_b%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: n_half
  !! n_half = n / 2
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_b_zz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%uhat_b%n) &
    CALL errore('evaluate_matsubara_b_zz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%uhat_b%m) &
    CALL errore('evaluate_matsubara_b_zz', 'wrong number of columns of output array.', 1)
  res(:, :) = czero
  !
  !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
  !      l2, obj%uhat_b%a, m, czero, res(:, :), l2)
  !
  n_half = n / 2
  !
  ALLOCATE(res_r(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error allocating res_r', 1)
  ALLOCATE(res_i(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error allocating res_i', 1)
  ALLOCATE(arr_half(l1, n_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error allocating arr_half', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ! take partial summations over odd l
    arr_half(:,:) = REAL(arr(:,1:(n-1):2))
    ! calculate the real part
    ! Re(arr) * Re(uhat_b)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
    !
    arr_half(:,:) = AIMAG(arr(:,1:(n-1):2))
    ! calculate the imaginary part
    ! Im(arr) * Re(uhat_b)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_odd, m, zero, res_i(:, :), l2)
    !
    res = CMPLX(res_r, res_i, KIND = DP)
    !
    ! take partial summations over even l
    arr_half(:,:) = AIMAG(arr(:,2:n:2))
    ! calculate the real part
    ! Im(arr) * Im(uhat_b)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_even, m, zero, res_r(:, :), l2)
    !
    arr_half(:,:) = REAL(arr(:,2:n:2))
    ! calculate the imaginary part
    ! Re(arr) * Im(uhat_b)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
    !
    res = res + CMPLX(-res_r, res_i, KIND = DP)
  ELSE
    ! take partial summations over odd l
    arr_half(:,:) = REAL(arr(:,1:(n-1):2))
    ! calculate the real part
    ! Re(arr) * Re(uhat_b)**T
    res_r(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
    !
    ! take partial summations over even l
    arr_half(:,:) = REAL(arr(:,2:n:2))
    ! calculate the imaginary part
    ! Re(arr) * Im(uhat_b)**T
    res_i(:, :) = zero
    CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
          l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
    !
    res = CMPLX(res_r, res_i, KIND = DP)
  ENDIF
  DEALLOCATE(res_r, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error deallocating res_r', 1)
  DEALLOCATE(res_i, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error deallocating res_i', 1)
  DEALLOCATE(arr_half, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_zz', 'Error deallocating arr_half', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_b_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_b_dz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a bosonic function onto sampling 
  !! Matsubara frequencies from the IR expansion coefficients, G_n.
  !!     G_m = \sum_n A_{mn} * G_n.
  !! A_{mn} is the n-th basis function of m,
  !! where m is a sampling Matsubara frequency.
  !!
  !! If positive_only is false, evaluate_matsubara_b_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a fermionic function constructed from IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: arr_half(:, :)
  !! dummy for arr
  INTEGER :: m
  !! number of columns of res, 
  !! which should be number of rows of A = obj%uhat_b%a
  INTEGER :: n
  !! number of columns of arr,
  !! which should be number of columns of A = obj%uhat_b%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: n_half
  !! n_half = n / 2
  INTEGER :: ierr
  !! Error status
  !
  l1 = SIZE(arr, 1)
  n = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  m = SIZE(res, 2)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_b_dz', 'wrong number of rows of input array.', 1)
  IF (n .NE. obj%uhat_b%n) &
    CALL errore('evaluate_matsubara_b_dz', 'wrong number of columns of input array.', 1)
  IF (m .NE. obj%uhat_b%m) &
    CALL errore('evaluate_matsubara_b_dz', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_matsubara_b_dz', 'input array should be a complex array.', 1)
  res(:, :) = czero
  !
  !CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
  !      l2, obj%uhat_b%a, m, czero, res(:, :), l2)
  !
  n_half = n / 2
  !
  ALLOCATE(res_r(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error allocating res_r', 1)
  ALLOCATE(res_i(l2, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error allocating res_i', 1)
  ALLOCATE(arr_half(l1, n_half), STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error allocating arr_half', 1)
  !
  ! take partial summations over odd l
  arr_half(:,:) = arr(:,1:(n-1):2)
  ! calculate the real part
  ! Re(arr) * Re(uhat_b)**T
  res_r(:, :) = zero
  CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
        l2, obj%uhat_b%a_odd, m, zero, res_r(:, :), l2)
  !
  ! take partial summations over even l
  arr_half(:,:) = arr(:,2:n:2)
  ! calculate the imaginary part
  ! Re(arr) * Im(uhat_b)**T
  res_i(:, :) = zero
  CALL DGEMM('n', 't', l1, m, n_half, one, arr_half(:,:), &
        l2, obj%uhat_b%a_even, m, zero, res_i(:, :), l2)
  !
  res = CMPLX(res_r, res_i, KIND = DP)
  DEALLOCATE(res_r, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error deallocating res_r', 1)
  DEALLOCATE(res_i, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error deallocating res_i', 1)
  DEALLOCATE(arr_half, STAT = ierr)
  IF (ierr /= 0) CALL errore('evaluate_matsubara_b_dz', 'Error deallocating arr_half', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_b_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_tau_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! where, A, U^T, V, s^{-1}, and G_m are 
  !! obj%u%a, obj%u%ut, obj%u%v, obj%u%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input function of imaginary time
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%u%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_tau_zz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('fit_tau_zz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('fit_tau_zz', 'wrong number of columns of output array.', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating ut_arr', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating res_r', 1)
    ALLOCATE(res_i(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating res_i', 1)
    ! calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_r(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_r, l2)
    !
    ! calculate the imaginary part
    arr_tmp = AIMAG(arr)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_i(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_i, l2)
    res = CMPLX(res_r, res_i, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating res_r', 1)
    DEALLOCATE(res_i, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating res_i', 1)
  ELSE
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error allocating res_r', 1)
    ! calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_r(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_r, l2)
    res = CMPLX(res_r, zero, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating res_r', 1)
  ENDIF
  !
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zz', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_tau_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_tau_dz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! where, A, U^T, V, s^{-1}, and G_m are 
  !! obj%u%a, obj%u%ut, obj%u%v, obj%u%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, fit_tau_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input function of imaginary time
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: res_tmp(:, :)
  !! dummy for res
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%u%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_tau_dz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('fit_tau_dz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('fit_tau_dz', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('fit_tau_dz', 'input array should be a complex array.', 1)
  ALLOCATE(res_tmp(l2, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dz', 'Error allocating res_tmp', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dz', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  res_tmp(:, :) = zero
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res_tmp, l2)
  !
  res(:, :) = CMPLX(res_tmp(:, :), zero, KIND = DP)
  DEALLOCATE(res_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dz', 'Error deallocating res_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dz', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_tau_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_tau_zd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! where, A, U^T, V, s^{-1}, and G_m are 
  !! obj%u%a, obj%u%ut, obj%u%v, obj%u%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, fit_tau_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input function of imaginary time
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%u%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_tau_zd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('fit_tau_zd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('fit_tau_zd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('fit_tau_zd', 'output array should be a complex array.', 1)
  ALLOCATE(arr_tmp(l1, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zd', 'Error allocating arr_tmp', 1)
  arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zd', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res, l2)
  !
  DEALLOCATE(arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zd', 'Error deallocating arr_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_zd', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_tau_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_tau_dd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a function of imaginary time.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling point of imaginary time.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! where, A, U^T, V, s^{-1}, and G_m are 
  !! obj%u%a, obj%u%ut, obj%u%v, obj%u%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, fit_tau_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! both the input and output arrays are real arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input function of imaginary time
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%u%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%u%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%u%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_tau_dd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%u%m) &
    CALL errore('fit_tau_dd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%u%n) &
    CALL errore('fit_tau_dd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('fit_tau_dd', 'input and output arrays should be complex arrays.', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dd', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%u%ut_real, ns, arr, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%u%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%u%v_real, n, zero, res, l2)
  !
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_tau_dd', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_tau_dd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_matsubara_f_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a fermionic function of Matsubara frequency.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling Matsubara frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! If positive_only = true, we can find them using
  !! U^T, \Sigma^-1, and V obtained from split_decompose:
  !! 
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * REAL[(U^T)_{im} * G_m].
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%uhat_f%a, obj%uhat_f%ut, 
  !! obj%uhat_f%v, obj%uhat_f%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input fermionic function of Matsubara frequency
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  REAL(KIND = DP), ALLOCATABLE :: res_tmp(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_r(:, :)
  !! dummy for ut_arr
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_tmp(:, :)
  !! dummy for ut_arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%uhat_f%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%uhat_f%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%uhat_f%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_matsubara_f_zz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%uhat_f%m) &
    CALL errore('fit_matsubara_f_zz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%uhat_f%n) &
    CALL errore('fit_matsubara_f_zz', 'wrong number of columns of output array.', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ALLOCATE(ut_arr(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error allocating ut_arr', 1)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = czero
    CALL ZGEMM("n", "t", ns, l1, m, cone, obj%uhat_f%ut, ns, arr, l1, czero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%uhat_f%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res(:, :) = czero
    CALL ZGEMM("t", "t", l1, n, ns, cone, ut_arr, ns, obj%uhat_f%v, n, czero, res, l2)
    DEALLOCATE(ut_arr, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error deallocating ut_arr', 1)
  ELSE
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_tmp(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error allocating res_tmp', 1)
    ALLOCATE(ut_arr_r(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error allocating ut_arr_r', 1)
    ALLOCATE(ut_arr_tmp(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error allocating ut_arr_tmp', 1)
    arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_r(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_f%ut_real, ns, arr_tmp, l1, zero, ut_arr_r, ns)
    !
    arr_tmp(:, :) = AIMAG(arr(:, :))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_tmp(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_f%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    !
    ut_arr_r(:, :) = ut_arr_r(:, :) - ut_arr_tmp(:, :)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr_r(j, i) = ut_arr_r(j, i) * obj%uhat_f%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_tmp(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr_r, ns, obj%uhat_f%v_real, n, zero, res_tmp, l2)
    !
    res = CMPLX(res_tmp, zero, KIND = DP)
    !
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error deallocating res_tmp', 1)
    DEALLOCATE(ut_arr_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error deallocating ut_arr_r', 1)
    DEALLOCATE(ut_arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_f_zz', 'Error deallocating ut_arr_tmp', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_matsubara_f_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_matsubara_f_zd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a fermionic function of Matsubara frequency.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling Matsubara frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! If positive_only = true, we can find them using
  !! U^T, \Sigma^-1, and V obtained from split_decompose:
  !! 
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * REAL[(U^T)_{im} * G_m].
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%uhat_f%a, obj%uhat_f%ut, 
  !! obj%uhat_f%v, obj%uhat_f%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, fit_matsubara_f_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input fermionic function of Matsubara frequency
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_tmp(:, :)
  !! dummy for ut_arr
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%uhat_f%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%uhat_f%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%uhat_f%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_matsubara_f_zd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%uhat_f%m) &
    CALL errore('fit_matsubara_f_zd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%uhat_f%n) &
    CALL errore('fit_matsubara_f_zd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('fit_matsubara_f_zd', 'output array should be a complex array.', 1)
  ALLOCATE(arr_tmp(l1, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error allocating arr_tmp', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error allocating ut_arr', 1)
  ALLOCATE(ut_arr_tmp(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error allocating ut_arr_tmp', 1)
  !
  arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_f%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
  !
  arr_tmp(:, :) = AIMAG(arr(:, :))
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr_tmp(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_f%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
  !
  ut_arr(:, :) = ut_arr(:, :) - ut_arr_tmp(:, :)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%uhat_f%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  res(:, :) = zero
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_f%v_real, n, zero, res, l2)
  !
  DEALLOCATE(arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error deallocating arr_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error deallocating ut_arr', 1)
  DEALLOCATE(ut_arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_f_zd', 'Error deallocating ut_arr_tmp', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_matsubara_f_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_matsubara_b_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a bosonic function of Matsubara frequency.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling Matsubara frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! If positive_only = true, we can find them using
  !! U^T, \Sigma^-1, and V obtained from split_decompose:
  !! 
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * REAL[(U^T)_{im} * G_m].
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%uhat_b%a, obj%uhat_b%ut, 
  !! obj%uhat_b%v, obj%uhat_b%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input bosonic function of Matsubara frequency
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  REAL(KIND = DP), ALLOCATABLE :: res_tmp(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_r(:, :)
  !! dummy for ut_arr
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_tmp(:, :)
  !! dummy for ut_arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%uhat_b%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%uhat_b%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%uhat_b%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_matsubara_b_zz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%uhat_b%m) &
    CALL errore('fit_matsubara_b_zz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%uhat_b%n) &
    CALL errore('fit_matsubara_b_zz', 'wrong number of columns of output array.', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ALLOCATE(ut_arr(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error allocating ut_arr', 1)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = czero
    CALL ZGEMM("n", "t", ns, l1, m, cone, obj%uhat_b%ut, ns, arr, l1, czero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%uhat_b%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res(:, :) = czero
    CALL ZGEMM("t", "t", l1, n, ns, cone, ut_arr, ns, obj%uhat_b%v, n, czero, res, l2)
    DEALLOCATE(ut_arr, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error deallocating ut_arr', 1)
  ELSE
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_tmp(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error allocating res_tmp', 1)
    ALLOCATE(ut_arr_r(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error allocating ut_arr_r', 1)
    ALLOCATE(ut_arr_tmp(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error allocating ut_arr_tmp', 1)
    arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_r(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_b%ut_real, ns, arr_tmp, l1, zero, ut_arr_r, ns)
    !
    arr_tmp(:, :) = AIMAG(arr(:, :))
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr_tmp(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_b%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
    !
    ut_arr_r(:, :) = ut_arr_r(:, :) - ut_arr_tmp(:, :)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr_r(j, i) = ut_arr_r(j, i) * obj%uhat_b%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_tmp(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr_r, ns, obj%uhat_b%v_real, n, zero, res_tmp, l2)
    !
    res = CMPLX(res_tmp, zero, KIND = DP)
    !
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error deallocating res_tmp', 1)
    DEALLOCATE(ut_arr_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error deallocating ut_arr_r', 1)
    DEALLOCATE(ut_arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zz', 'Error deallocating ut_arr_tmp', 1)
  ENDIF
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_matsubara_b_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fit_matsubara_b_zd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine finds IR expansion coefficients G_n satisfying 
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is a bosonic function of Matsubara frequency.
  !! A_{mn} is the n-th basis function of m, 
  !! where m is a sampling Matsubara frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! If positive_only = true, we can find them using
  !! U^T, \Sigma^-1, and V obtained from split_decompose:
  !! 
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * REAL[(U^T)_{im} * G_m].
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%uhat_b%a, obj%uhat_b%ut, 
  !! obj%uhat_b%v, obj%uhat_b%inv_s, and arr, respectively.
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, fit_matsubara_b_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! input fermionic function of Matsubara frequency
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: ut_arr_tmp(:, :)
  !! dummy for ut_arr
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%uhat_f%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%uhat_f%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%uhat_b%ns
  IF (l1 .NE. l2) &
    CALL errore('fit_matsubara_b_zd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%uhat_b%m) &
    CALL errore('fit_matsubara_b_zd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%uhat_b%n) &
    CALL errore('fit_matsubara_b_zd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('fit_matsubara_b_zd', 'output array should be a complex array.', 1)
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error allocating arr_tmp', 1)
    ALLOCATE(ut_arr(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error allocating ut_arr', 1)
    ALLOCATE(ut_arr_tmp(ns, l1), STAT = ierr)
    IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error allocating ut_arr_tmp', 1)
  !
  arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_b%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
  !
  arr_tmp(:, :) = AIMAG(arr(:, :))
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr_tmp(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%uhat_b%ut_imag, ns, arr_tmp, l1, zero, ut_arr_tmp, ns)
  !
  ut_arr(:, :) = ut_arr(:, :) - ut_arr_tmp(:, :)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%uhat_b%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  res(:, :) = zero
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%uhat_b%v_real, n, zero, res, l2)
  !
  DEALLOCATE(arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error deallocating arr_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error deallocating ut_arr', 1)
  DEALLOCATE(ut_arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('fit_matsubara_b_zd', 'Error deallocating ut_arr_tmp', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE fit_matsubara_b_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE to_dlr_zz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine transfers expansion coefficients in the IR basis
  !! to that in the DLR basis.
  !! DLR expansion coefficients G_n satisfy
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is an IR expansion coefficient.
  !! A_{mn} is defined as:
  !!     A_{mn} = - s_{m} * V_{mn}
  !! where V is right-singular functions, which is obtained by
  !! the following singular value expansion:
  !!     K(tau_s, omega_m) : 
  !!     K_{sn} = \sum_m U_{sm} * s_{m} * V_{mn},
  !! where K is an integral kernel and omega_m is a real frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%dlr%a, obj%dlr%ut, 
  !! obj%dlr%v, obj%dlr%inv_s, and arr, respectively.
  !! (Do not confuse these V and s with the above V and s.)
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  REAL(KIND = DP), ALLOCATABLE :: res_r(:, :)
  !! dummy for res
  REAL(KIND = DP), ALLOCATABLE :: res_i(:, :)
  !! dummy for res
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%dlr%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%dlr%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%dlr%ns
  IF (l1 .NE. l2) &
    CALL errore('to_dlr_zz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%dlr%m) &
    CALL errore('to_dlr_zz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%dlr%n) &
    CALL errore('to_dlr_zz', 'wrong number of columns of output array.', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating ut_arr', 1)
  !
  IF (.NOT. obj%positive_only) THEN
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating res_r', 1)
    ALLOCATE(res_i(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating res_i', 1)
    ! calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_r(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res_r, l2)
    !
    ! calculate the imaginary part
    arr_tmp = AIMAG(arr)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_i(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res_i, l2)
    res = CMPLX(res_r, res_i, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating res_r', 1)
    DEALLOCATE(res_i, STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating res_i', 1)
  ELSE
    ALLOCATE(arr_tmp(l1, m), STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating arr_tmp', 1)
    ALLOCATE(res_r(l2, n), STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error allocating res_r', 1)
    ! calculate the real part
    arr_tmp = REAL(arr, KIND = DP)
    !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
    ut_arr(:, :) = zero
    CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
    DO j = 1, ns
      DO i = 1, l1
        ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
      ENDDO
    ENDDO
    !
    ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
    res_r(:, :) = zero
    CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res_r, l2)
    res = CMPLX(res_r, zero, KIND = DP)
    DEALLOCATE(arr_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating arr_tmp', 1)
    DEALLOCATE(res_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating res_r', 1)
  ENDIF
  !
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zz', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE to_dlr_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE to_dlr_dz(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine transfers expansion coefficients in the IR basis
  !! to that in the DLR basis.
  !! DLR expansion coefficients G_n satisfy
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is an IR expansion coefficient.
  !! A_{mn} is defined as:
  !!     A_{mn} = V_{mn} / s_{m}
  !! where V is right-singular functions, which is obtained by
  !! the following singular value expansion:
  !!     K(tau_s, omega_m) : 
  !!     K_{sn} = \sum_m U_{sm} * s_{m} * V_{mn},
  !! where K is an integral kernel and omega_m is a real frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%dlr%a, obj%dlr%ut, 
  !! obj%dlr%v, obj%dlr%inv_s, and arr, respectively.
  !! (Do not confuse these V and s with the above V and s.)
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, to_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: res_tmp(:, :)
  !! dummy for res
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%dlr%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%dlr%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%dlr%ns
  IF (l1 .NE. l2) &
    CALL errore('to_dlr_dz', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%dlr%m) &
    CALL errore('to_dlr_dz', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%dlr%n) &
    CALL errore('to_dlr_dz', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('to_dlr_dz', 'input array should be a complex array.', 1)
  ALLOCATE(res_tmp(l2, n), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dz', 'Error allocating res_tmp', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dz', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  res_tmp(:, :) = zero
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res_tmp, l2)
  !
  res(:, :) = CMPLX(res_tmp(:, :), zero, KIND = DP)
  DEALLOCATE(res_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dz', 'Error deallocating res_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dz', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE to_dlr_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE to_dlr_zd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine transfers expansion coefficients in the IR basis
  !! to that in the DLR basis.
  !! DLR expansion coefficients G_n satisfy
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is an IR expansion coefficient.
  !! A_{mn} is defined as:
  !!     A_{mn} = V_{mn} / s_{m}
  !! where V is right-singular functions, which is obtained by
  !! the following singular value expansion:
  !!     K(tau_s, omega_m) : 
  !!     K_{sn} = \sum_m U_{sm} * s_{m} * V_{mn},
  !! where K is an integral kernel and omega_m is a real frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%dlr%a, obj%dlr%ut, 
  !! obj%dlr%v, obj%dlr%inv_s, and arr, respectively.
  !! (Do not confuse these V and s with the above V and s.)
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, to_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  REAL(KIND = DP), ALLOCATABLE :: arr_tmp(:, :)
  !! dummy for arr
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%dlr%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%dlr%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%dlr%ns
  IF (l1 .NE. l2) &
    CALL errore('to_dlr_zd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%dlr%m) &
    CALL errore('to_dlr_zd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%dlr%n) &
    CALL errore('to_dlr_zd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('to_dlr_zd', 'output array should be a complex array.', 1)
  ALLOCATE(arr_tmp(l1, m), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zd', 'Error allocating arr_tmp', 1)
  arr_tmp(:, :) = REAL(arr(:, :), KIND = DP)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zd', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr_tmp, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res, l2)
  !
  DEALLOCATE(arr_tmp, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zd', 'Error deallocating arr_tmp', 1)
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_zd', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE to_dlr_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE to_dlr_dd(obj, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine transfers expansion coefficients in the IR basis
  !! to that in the DLR basis.
  !! DLR expansion coefficients G_n satisfy
  !! the following relation:
  !!     G_m = \sum_n A_{mn} * G_n.
  !! G_m is an IR expansion coefficient.
  !! A_{mn} is defined as:
  !!     A_{mn} = V_{mn} / s_{m}
  !! where V is right-singular functions, which is obtained by
  !! the following singular value expansion:
  !!     K(tau_s, omega_m) : 
  !!     K_{sn} = \sum_m U_{sm} * s_{m} * V_{mn},
  !! where K is an integral kernel and omega_m is a real frequency.
  !!
  !! Assuming the matrix A is already SVDecomposed,
  !!
  !!     A_{mn} = \sum_i U_{mi} * s_{i} * (V^T)_{in},
  !!
  !! we can easily find the expansion coefficients G_n as follows:
  !!
  !!     G_n = \sum_i \sum_m V_{ni} * s^{-1}_{i} * (U^T)_{im} * G_m.
  !!
  !! A, U^T, V, s^{-1}, and G_m are obj%dlr%a, obj%dlr%ut, 
  !! obj%dlr%v, obj%dlr%inv_s, and arr, respectively.
  !! (Do not confuse these V and s with the above V and s.)
  !!
  !! Please note: Numerical stability is achieved by multiplying G_m by
  !! U^T, \Sigma^-1, and V in sequence, rather than multiplying by
  !! the pseudo-inverse matrix A^-1 = U^T \Sigma^-1 V.
  !!
  !! If positive_only is false, to_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! both the input and output arrays are real arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! IR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), ALLOCATABLE :: ut_arr(:, :)
  !! temporarily stores U^T * G
  !
  INTEGER :: m
  !! number of columns of arr,
  !! which should be number of rows of A = obj%dlr%a
  INTEGER :: n
  !! number of columns of res,
  !! which should be number of columns of A = obj%dlr%a
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  INTEGER :: ns
  !! number of non-zero singular values
  INTEGER :: i
  !! counter
  INTEGER :: j
  !! counter
  INTEGER :: ierr
  !! Error status
  !
  ! ut(ns, m)
  ! v(n, ns)
  ! arr(l1, m)
  ! mat(m, n)
  ! ut_arr(ns, l1)
  ! res(l2, n)
  l1 = SIZE(arr, 1)
  m = SIZE(arr, 2)
  l2 = SIZE(res, 1)
  n = SIZE(res, 2)
  ns = obj%dlr%ns
  IF (l1 .NE. l2) &
    CALL errore('to_dlr_dd', 'wrong number of rows of input array.', 1)
  IF (m .NE. obj%dlr%m) &
    CALL errore('to_dlr_dd', 'wrong number of columns of input array.', 1)
  IF (n .NE. obj%dlr%n) &
    CALL errore('to_dlr_dd', 'wrong number of columns of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('to_dlr_dd', 'input and output arrays should be complex arrays.', 1)
  ALLOCATE(ut_arr(ns, l1), STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dd', 'Error allocating ut_arr', 1)
  !
  !ut(ns, m) * arr(l1, m) -> ut_arr(ns, l1)
  ut_arr(:, :) = zero
  CALL DGEMM("n", "t", ns, l1, m, one, obj%dlr%ut_real, ns, arr, l1, zero, ut_arr, ns)
  DO j = 1, ns
    DO i = 1, l1
      ut_arr(j, i) = ut_arr(j, i) * obj%dlr%inv_s(j)
    ENDDO
  ENDDO
  !
  ! ut_arr(ns, l1) * v(n, ns) -> res(l2, n)
  CALL DGEMM("t", "t", l1, n, ns, one, ut_arr, ns, obj%dlr%v_real, n, zero, res, l2)
  !
  DEALLOCATE(ut_arr, STAT = ierr)
  IF (ierr /= 0) CALL errore('to_dlr_dd', 'Error deallocating ut_arr', 1)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE to_dlr_dd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_from_dlr_zz(obj, tau, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto given points of imaginary time
  !! from the DLR expansion coefficients.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: tau(:)
  !! points of imaginary time, which must be in [0, beta]
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from DLR expansion coefficients
  REAL(KIND = DP) :: kernel
  !! logistic kernel, K(tau, omega)
  INTEGER :: ntau
  !! number of points of imaginary time, tau
  INTEGER :: nt
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: t
  !! counter on points of imaginary time
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  ntau = SIZE(tau)
  nt = SIZE(res, 2)
  IF (ntau .NE. nt) &
    CALL errore('evaluate_tau_from_dlr_zz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_from_dlr_zz', 'wrong number of rows of output array.', 1)
  !
  DO t = 1, ntau
   IF ((tau(t) < 0.0D0) .OR. (obj%beta < tau(t))) THEN
    CALL errore('evaluate_tau_from_dlr_zz', 'tau must be in [0, beta].', 1)
   ENDIF
  ENDDO
  !
  res(:, :) = czero
  !
  DO t = 1, ntau
    DO p = 1, obj%nomega
      IF (obj%omega(p) < 0) THEN
        IF ( (obj%beta - tau(t)) * obj%omega(p) < -100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) < -30.d0) THEN
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p))
        ELSE
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p)) / (exp(obj%beta * obj%omega(p)) + 1.0D0)
        ENDIF
      ELSE
        IF (tau(t) * obj%omega(p) > 100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) > 30.d0) THEN
          kernel = exp(- tau(t) * obj%omega(p))
        ELSE
          kernel = exp(- tau(t) * obj%omega(p)) / (1.0D0 + exp(- obj%beta * obj%omega(p)))
        ENDIF
      ENDIF
      res(:, t) = res(:, t) - kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_from_dlr_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_from_dlr_dz(obj, tau, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto given points of imaginary time
  !! from the DLR expansion coefficients.
  !!
  !! If positive_only is false, evaluate_tau_from_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: tau(:)
  !! points of imaginary time, which must be in [0, beta]
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from DLR expansion coefficients
  REAL(KIND = DP) :: kernel
  !! logistic kernel, K(tau, omega)
  INTEGER :: ntau
  !! number of points of imaginary time, tau
  INTEGER :: nt
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: t
  !! counter on points of imaginary time
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  ntau = SIZE(tau)
  nt = SIZE(res, 2)
  IF (ntau .NE. nt) &
    CALL errore('evaluate_tau_from_dlr_dz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_from_dlr_dz', 'wrong number of rows of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_from_dlr_dz', 'input array should be a complex array.', 1)
  !
  DO t = 1, ntau
   IF ((tau(t) < 0.0D0) .OR. (obj%beta < tau(t))) THEN
    CALL errore('evaluate_tau_from_dlr_dz', 'tau must be in [0, beta].', 1)
   ENDIF
  ENDDO
  !
  res(:, :) = czero
  !
  DO t = 1, ntau
    DO p = 1, obj%nomega
      IF (obj%omega(p) < 0) THEN
        IF ( (obj%beta - tau(t)) * obj%omega(p) < -100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) < -30.d0) THEN
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p))
        ELSE
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p)) / (exp(obj%beta * obj%omega(p)) + 1.0D0)
        ENDIF
      ELSE
        IF (tau(t) * obj%omega(p) > 100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) > 30.d0) THEN
          kernel = exp(- tau(t) * obj%omega(p))
        ELSE
          kernel = exp(- tau(t) * obj%omega(p)) / (1.0D0 + exp(- obj%beta * obj%omega(p)))
        ENDIF
      ENDIF
      res(:, t) = res(:, t) - kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_from_dlr_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_from_dlr_zd(obj, tau, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto given points of imaginary time
  !! from the DLR expansion coefficients.
  !!
  !! If positive_only is false, evaluate_tau_from_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is complex and the output array is real.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: tau(:)
  !! points of imaginary time, which must be in [0, beta]
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from DLR expansion coefficients
  REAL(KIND = DP) :: kernel
  !! logistic kernel, K(tau, omega)
  INTEGER :: ntau
  !! number of points of imaginary time, tau
  INTEGER :: nt
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: t
  !! counter on points of imaginary time
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  ntau = SIZE(tau)
  nt = SIZE(res, 2)
  IF (ntau .NE. nt) &
    CALL errore('evaluate_tau_from_dlr_zd', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_from_dlr_zd', 'wrong number of rows of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_from_dlr_zd', 'output array should be a complex array.', 1)
  !
  DO t = 1, ntau
   IF ((tau(t) < 0.0D0) .OR. (obj%beta < tau(t))) THEN
    CALL errore('evaluate_tau_from_dlr_zd', 'tau must be in [0, beta].', 1)
   ENDIF
  ENDDO
  !
  res(:, :) = czero
  !
  DO t = 1, ntau
    DO p = 1, obj%nomega
      IF (obj%omega(p) < 0) THEN
        IF ( (obj%beta - tau(t)) * obj%omega(p) < -100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) < -30.d0) THEN
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p))
        ELSE
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p)) / (exp(obj%beta * obj%omega(p)) + 1.0D0)
        ENDIF
      ELSE
        IF (tau(t) * obj%omega(p) > 100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) > 30.d0) THEN
          kernel = exp(- tau(t) * obj%omega(p))
        ELSE
          kernel = exp(- tau(t) * obj%omega(p)) / (1.0D0 + exp(- obj%beta * obj%omega(p)))
        ENDIF
      ENDIF
      res(:, t) = res(:, t) - kernel * REAL(arr(:, p), KIND = DP)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_from_dlr_zd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_tau_from_dlr_dd(obj, tau, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a function onto given points of imaginary time
  !! from the DLR expansion coefficients.
  !!
  !! If positive_only is false, evaluate_tau_from_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! both the input and output arrays are real arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  REAL(KIND = DP), INTENT(IN) :: tau(:)
  !! points of imaginary time, which must be in [0, beta]
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  REAL(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a function constructed from DLR expansion coefficients
  REAL(KIND = DP) :: kernel
  !! logistic kernel, K(tau, omega)
  INTEGER :: ntau
  !! number of points of imaginary time, tau
  INTEGER :: nt
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: t
  !! counter on points of imaginary time
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  ntau = SIZE(tau)
  nt = SIZE(res, 2)
  IF (ntau .NE. nt) &
    CALL errore('evaluate_tau_from_dlr_dd', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_tau_from_dlr_dd', 'wrong number of rows of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_tau_from_dlr_dd', 'input and output arrays should be complex arrays.', 1)
  !
  DO t = 1, ntau
   IF ((tau(t) < 0.0D0) .OR. (obj%beta < tau(t))) THEN
    CALL errore('evaluate_tau_from_dlr_dd', 'tau must be in [0, beta].', 1)
   ENDIF
  ENDDO
  !
  res(:, :) = czero
  !
  DO t = 1, ntau
    DO p = 1, obj%nomega
      IF (obj%omega(p) < 0) THEN
        IF ( (obj%beta - tau(t)) * obj%omega(p) < -100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) < -30.d0) THEN
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p))
        ELSE
          kernel = exp( (obj%beta - tau(t)) * obj%omega(p)) / (exp(obj%beta * obj%omega(p)) + 1.0D0)
        ENDIF
      ELSE
        IF (tau(t) * obj%omega(p) > 100.d0) THEN
          kernel = 0.0D0
        ELSEIF (obj%beta * obj%omega(p) > 30.d0) THEN
          kernel = exp(- tau(t) * obj%omega(p))
        ELSE
          kernel = exp(- tau(t) * obj%omega(p)) / (1.0D0 + exp(- obj%beta * obj%omega(p)))
        ENDIF
      ENDIF
      res(:, t) = res(:, t) - kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_tau_from_dlr_dd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_f_from_dlr_zz(obj, freq, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a fermionic function onto given Matsubara 
  !! frequencies from the DLR expansion coefficients.
  !! Matsubara frequencies must be given as odd integers.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  INTEGER, INTENT(IN) :: freq(:)
  !! integer parts of Matsubara freqs, which must be odd integers
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a fermionic function constructed from DLR expansion coefficients
  COMPLEX(KIND = DP) :: cfreq
  !! Matsubara frequencies
  COMPLEX(KIND = DP) :: kernel
  !! logistic kernel, K(iw_n, omega)
  INTEGER :: nfreq
  !! number of points of Matsubara freqs
  INTEGER :: nf
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  nfreq = SIZE(freq)
  nf = SIZE(res, 2)
  IF (nfreq .NE. nf) &
    CALL errore('evaluate_matsubara_f_from_dlr_zz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_f_from_dlr_zz', 'wrong number of rows of output array.', 1)
  !
  DO n = 1, nfreq
   IF (MOD(freq(n), 2) == 0) &
     CALL errore('evaluate_matsubara_f_from_dlr_zz', 'one of input integers is not odd.', 1)
  ENDDO
  !
  res(:, :) = czero
  !
  DO n = 1, nfreq
    cfreq = ci * pi * REAL(freq(n), KIND = DP) / obj%beta
    DO p = 1, obj%nomega
      kernel = 1.0D0 / (cfreq - obj%omega(p))
      res(:, n) = res(:, n) + kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_f_from_dlr_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_f_from_dlr_dz(obj, freq, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a fermionic function onto given Matsubara 
  !! frequencies from the DLR expansion coefficients.
  !! Matsubara frequencies must be given as odd integers.
  !!
  !! If positive_only is false, evaluate_matsubara_f_from_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  INTEGER, INTENT(IN) :: freq(:)
  !! integer parts of Matsubara freqs, which must be odd integers
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a fermionic function constructed from DLR expansion coefficients
  COMPLEX(KIND = DP) :: cfreq
  !! Matsubara frequencies
  COMPLEX(KIND = DP) :: kernel
  !! logistic kernel, K(iw_n, omega)
  INTEGER :: nfreq
  !! number of points of Matsubara freqs
  INTEGER :: nf
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  nfreq = SIZE(freq)
  nf = SIZE(res, 2)
  IF (nfreq .NE. nf) &
    CALL errore('evaluate_matsubara_f_from_dlr_dz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_f_from_dlr_dz', 'wrong number of rows of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_matsubara_f_from_dlr_dz', 'input array should be a complex array.', 1)
  !
  DO n = 1, nfreq
   IF (MOD(freq(n), 2) == 0) &
    CALL errore('evaluate_matsubara_f_from_dlr_dz', 'one of input integers is not odd.', 1)
  ENDDO
  !
  res(:, :) = czero
  !
  DO n = 1, nfreq
    cfreq = ci * pi * REAL(freq(n), KIND = DP) / obj%beta
    DO p = 1, obj%nomega
      kernel = 1.0D0 / (cfreq - obj%omega(p))
      res(:, n) = res(:, n) + kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_f_from_dlr_dz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_b_from_dlr_zz(obj, freq, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a bosonic function onto given Matsubara 
  !! frequencies from the DLR expansion coefficients.
  !! Matsubara frequencies must be given as even integers.
  !!
  !! This subroutine is called when
  !! both the input and output arrays are complex arrays.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  INTEGER, INTENT(IN) :: freq(:)
  !! integer parts of Matsubara freqs, which must be even integers
  COMPLEX(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a bosonic function constructed from DLR expansion coefficients
  COMPLEX(KIND = DP) :: cfreq
  !! Matsubara frequencies
  COMPLEX(KIND = DP) :: kernel
  !! logistic kernel, K(iw_n, omega)
  INTEGER :: nfreq
  !! number of points of Matsubara freqs
  INTEGER :: nf
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  nfreq = SIZE(freq)
  nf = SIZE(res, 2)
  IF (nfreq .NE. nf) &
    CALL errore('evaluate_matsubara_b_from_dlr_zz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_b_from_dlr_zz', 'wrong number of rows of output array.', 1)
  !
  DO n = 1, nfreq
   IF (MOD(freq(n), 2) .ne. 0) &
    CALL errore('evaluate_matsubara_b_from_dlr_zz', 'one of input integers is not even.', 1)
  ENDDO
  !
  res(:, :) = czero
  !
  DO n = 1, nfreq
    cfreq = ci * pi * REAL(freq(n), KIND = DP) / obj%beta
    DO p = 1, obj%nomega
      kernel = TANH(5.0d-1 * obj%beta * obj%omega(p)) / (cfreq - obj%omega(p))
      res(:, n) = res(:, n) + kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_b_from_dlr_zz
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_matsubara_b_from_dlr_dz(obj, freq, arr, res)
  !-----------------------------------------------------------------------
  !!
  !! This routine constructs a bosonic function onto given Matsubara 
  !! frequencies from the DLR expansion coefficients.
  !! Matsubara frequencies must be given as even integers.
  !!
  !! If positive_only is false, evaluate_matsubara_b_from_dlr_zz
  !! should be called instead of this subroutine.
  !! This subroutine is called when
  !! the input array is real and the output array is complex.
  !!
  !
  TYPE(IR), INTENT(IN) :: obj
  !! contains all the IR-basis objects
  INTEGER, INTENT(IN) :: freq(:)
  !! integer parts of Matsubara freqs, which must be even integers
  REAL(KIND = DP), INTENT(IN) :: arr(:, :)
  !! DLR expansion coefficients
  COMPLEX(KIND = DP), INTENT(OUT) :: res(:, :)
  !! a bosonic function constructed from DLR expansion coefficients
  COMPLEX(KIND = DP) :: cfreq
  !! Matsubara frequencies
  COMPLEX(KIND = DP) :: kernel
  !! logistic kernel, K(iw_n, omega)
  INTEGER :: nfreq
  !! number of points of Matsubara freqs
  INTEGER :: nf
  !! number of columns of res
  INTEGER :: p
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: l1
  !! number of rows of arr
  INTEGER :: l2
  !! number of rows of res
  !
  nfreq = SIZE(freq)
  nf = SIZE(res, 2)
  IF (nfreq .NE. nf) &
    CALL errore('evaluate_matsubara_b_from_dlr_dz', 'wrong number of columns of output array.', 1)
  l1 = SIZE(arr, 1)
  l2 = SIZE(res, 1)
  IF (l1 .NE. l2) &
    CALL errore('evaluate_matsubara_b_from_dlr_dz', 'wrong number of rows of output array.', 1)
  IF (.NOT. obj%positive_only) &
    CALL errore('evaluate_matsubara_b_from_dlr_dz', 'input array should be a complex array.', 1)
  !
  DO n = 1, nfreq
   IF (MOD(freq(n), 2) .ne. 0) &
    CALL errore('evaluate_matsubara_b_from_dlr_dz', 'one of input integers is not even.', 1)
  ENDDO
  !
  res(:, :) = czero
  !
  DO n = 1, nfreq
    cfreq = ci * pi * REAL(freq(n), KIND = DP) / obj%beta
    DO p = 1, obj%nomega
      kernel = TANH(5.0d-1 * obj%beta * obj%omega(p)) / (cfreq - obj%omega(p))
      res(:, n) = res(:, n) + kernel * arr(:, p)
    ENDDO
  ENDDO
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE evaluate_matsubara_b_from_dlr_dz
  !-----------------------------------------------------------------------
  !
  !----------------------------------------------------------------------
  END MODULE sparse_ir
  !----------------------------------------------------------------------