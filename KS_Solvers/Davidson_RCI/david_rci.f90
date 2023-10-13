!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
!               2017 M. Oliveira
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )

module david_rci_m
  use util_param
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id, &
                            me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_bcast
  implicit none
  private

  public :: david_rci_run,        &
            david_rci_work_alloc, &
            david_rci_work_free,  &
            david_rci_work_t

  integer, parameter, public :: &
    ESL_TASK_NONE = 0,          &
    ESL_TASK_EXIT = 1,          &
    ESL_TASK_HPSI = 2,          &
    ESL_TASK_SPSI = 4,          &
    ESL_TASK_GPSI = 8

  type david_rci_work_t
    integer :: ivec_start, ivec_size
    complex(DP), allocatable :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
    real(DP), allocatable :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  end type david_rci_work_t
  
contains

  subroutine david_rci_work_alloc(npw, npwx, nvec, nvecx, npol, uspp, work)
    integer, intent(in) :: npw, npwx, nvec, nvecx, npol
    logical, intent(in) :: uspp
    type(david_rci_work_t), intent(inout) :: work
    
    integer :: ierr

    ALLOCATE(  work%psi( npwx, npol, nvecx ), STAT=ierr )
    IF( ierr /= 0 ) &
      CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )

    ALLOCATE( work%hpsi( npwx, npol, nvecx ), STAT=ierr )
    IF( ierr /= 0 ) &
      CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
    !
    IF ( uspp ) THEN
      ALLOCATE( work%spsi( npwx, npol, nvecx ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
    END IF
    !
    ALLOCATE( work%ew( nvecx ), STAT=ierr )
    IF( ierr /= 0 ) &
      CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )

  end subroutine david_rci_work_alloc

  subroutine david_rci_work_free(work)
    type(david_rci_work_t), intent(inout) :: work

    if (allocated(work%psi)) then
      deallocate(work%psi)
    end if

    if (allocated(work%hpsi)) then
      deallocate(work%hpsi)
    end if

    if (allocated(work%spsi)) then
      deallocate(work%spsi)
    end if

    if (allocated(work%ew)) then
      deallocate(work%ew)
    end if

  end subroutine david_rci_work_free
  
  subroutine david_rci_run(npw, npwx, nvec, nvecx, npol, evc, ethr, uspp, e, btype, notcnv, lrot, dav_iter, work, task)
    !----------------------------------------------------------------------------
    !
    ! ... iterative solution of the eigenvalue problem:
    !
    ! ... ( H - e S ) * evc = 0
    !
    ! ... where H is an hermitean operator, e is a real scalar,
    ! ... S is an overlap matrix, evc is a complex vector
    !
    INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
    COMPLEX(DP), INTENT(INOUT) :: evc(npwx,npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
    REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
    LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
    INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
    LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
    REAL(DP), INTENT(INOUT) :: e(nvec)
    ! contains the estimated roots.
    INTEGER, INTENT(INOUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
    type(david_rci_work_t), intent(inout) :: work
    ! work space
    integer, intent(out) :: task
    ! Next task to be performed by the calling program
    !
    !
    include 'laxlib.fh'
    !
    ! ... LOCAL variables
    !
    INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
    !
    INTEGER, save :: nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
    INTEGER :: ierr
    COMPLEX(DP), ALLOCATABLE, save :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
    LOGICAL, ALLOCATABLE, save  :: conv(:)
    ! true if the root is converged
    REAL(DP), save :: empty_ethr 
    ! threshold for empty bands
    integer, save :: step = 1
    ! current calculation step
    REAL(DP), EXTERNAL :: ddot
    
    select case (step)
    case(1)
      CALL start_clock( 'cegterg' )
      !
      IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
      !
      ! ... threshold for empty bands
      !
      empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
      !
      IF ( npol == 1 ) THEN
        !
        kdim = npw
        kdmx = npwx
        !
      ELSE
        !
        kdim = npwx*npol
        kdmx = npwx*npol
        !
      END IF
      !
      ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate sc ', ABS(ierr) )
      ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate hc ', ABS(ierr) )
      ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
      ALLOCATE( conv( nvec ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
      !
      notcnv = nvec
      nbase  = nvec
      conv   = .FALSE.
      !
      IF ( uspp ) work%spsi = ZERO
      !
      work%hpsi = ZERO
      work%psi  = ZERO
      work%psi(:,:,1:nvec) = evc(:,:,1:nvec)

      dav_iter = 0

      work%ivec_start = 1
      work%ivec_size  = nvec
      task = ESL_TASK_HPSI
      if (uspp) task = task + ESL_TASK_SPSI
      
      step = step + 1

    case(2)
      ! ... hc contains the projection of the hamiltonian onto the reduced 
      ! ... space vc contains the eigenvectors of hc
      !
      hc(:,:) = ZERO
      sc(:,:) = ZERO
      vc(:,:) = ZERO
      !
      CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
        work%psi, kdmx, work%hpsi, kdmx, ZERO, hc, nvecx )
      !
      CALL mp_sum( hc( :, 1:nbase ), intra_bgrp_comm )
      !
      IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
          work%psi, kdmx, work%spsi, kdmx, ZERO, sc, nvecx )
        !     
      ELSE
        !
        CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
          work%psi, kdmx, work%psi, kdmx, ZERO, sc, nvecx )
        !
      END IF
      !
      CALL mp_sum( sc( :, 1:nbase ), intra_bgrp_comm )
      !
      IF ( lrot ) THEN
        !
        DO n = 1, nbase
          !
          e(n) = REAL( hc(n,n) )
          !
          vc(n,n) = ONE
          !
        END DO
        !
      ELSE
        !
        ! ... diagonalize the reduced hamiltonian
        !
        IF( my_bgrp_id == root_bgrp_id ) THEN
          CALL diaghg( nbase, nvec, hc, sc, nvecx, work%ew, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
        END IF
        IF( nbgrp > 1 ) THEN
          CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
          CALL mp_bcast( work%ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
        !
        e(1:nvec) = work%ew(1:nvec)
        !
      END IF

      dav_iter = 0

      step = step + 1
      task = ESL_TASK_NONE
      
    case (3)
      !
      ! ... iterate
      !
      dav_iter = dav_iter + 1
      !
      CALL start_clock( 'cegterg:update' )
      !
      np = 0
      !
      DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
          !
          ! ... this root not yet converged ... 
          !
          np = np + 1
          !
          ! ... reorder eigenvectors so that coefficients for unconverged
          ! ... roots come first. This allows to use quick matrix-matrix 
          ! ... multiplications to set a new basis vector (see below)
          !
          IF ( np /= n ) vc(:,np) = vc(:,n)
          !
          ! ... for use in g_psi_ptr
          !
          work%ew(nbase+np) = e(n)
          !
        END IF
        !
      END DO
      !
      nb1 = nbase + 1
      !
      ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
      !
      IF ( uspp ) THEN
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, work%spsi, &
          kdmx, vc, nvecx, ZERO, work%psi(1,1,nb1), kdmx )
        !     
      ELSE
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, work%psi, &
          kdmx, vc, nvecx, ZERO, work%psi(1,1,nb1), kdmx )
        !
      END IF
      !
      DO np = 1, notcnv
        !
        work%psi(:,:,nbase+np) = - work%ew(nbase+np)*work%psi(:,:,nbase+np)
        !
      END DO
      !
      CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, work%hpsi, &
        kdmx, vc, nvecx, ONE, work%psi(1,1,nb1), kdmx )
      !
      CALL stop_clock( 'cegterg:update' )
      !
      ! ... approximate inverse iteration
      !
      
      work%ivec_start = nb1
      work%ivec_size  = notcnv
      step = step + 1
      task = ESL_TASK_GPSI

    case (4)
      !
      ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
      ! ... order to improve numerical stability of subspace diagonalization
      ! ... (diaghg) ew is used as work array :
      !
      ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
      !
      DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
          !
          work%ew(n) = ddot( 2*npw, work%psi(1,1,nbn), 1, work%psi(1,1,nbn), 1 )
          !
        ELSE
          !
          work%ew(n) = ddot( 2*npw, work%psi(1,1,nbn), 1, work%psi(1,1,nbn), 1 ) + &
            ddot( 2*npw, work%psi(1,2,nbn), 1, work%psi(1,2,nbn), 1 )
          !
        END IF
        !
      END DO
      !
      CALL mp_sum( work%ew( 1:notcnv ), intra_bgrp_comm )
      !
      DO n = 1, notcnv
        !
        work%psi(:,:,nbase+n) = work%psi(:,:,nbase+n) / SQRT( work%ew(n) )
        !
      END DO
      !
      ! ... here compute the hpsi and spsi of the new functions
      !
      !
      work%ivec_start = nb1
      work%ivec_size  = notcnv
      step = step + 1
      task = ESL_TASK_HPSI
      if (uspp) task = task + ESL_TASK_SPSI

    case (5)
      !
      ! ... update the reduced hamiltonian
      !
      CALL start_clock( 'cegterg:overlap' )
      !
      CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, work%psi, &
        kdmx, work%hpsi(1,1,nb1), kdmx, ZERO, hc(1,nb1), nvecx )
      !
      CALL mp_sum( hc( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
      !
      IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, work%psi, &
          kdmx, work%spsi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !     
      ELSE
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, work%psi, &
          kdmx, work%psi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !
      END IF
      !
      CALL mp_sum( sc( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
      !
      CALL stop_clock( 'cegterg:overlap' )
      !
      nbase = nbase + notcnv
      !
      DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real 
        !
        hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
        sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        !
        DO m = n + 1, nbase
          !
          hc(m,n) = CONJG( hc(n,m) )
          sc(m,n) = CONJG( sc(n,m) )
          !
        END DO
        !
      END DO
      !
      ! ... diagonalize the reduced hamiltonian
      !
      IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc, sc, nvecx, work%ew, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
      END IF
      IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( work%ew, root_bgrp_id, inter_bgrp_comm )
      ENDIF
      !
      ! ... test for convergence
      !
      WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( work%ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
      ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( work%ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
      END WHERE
      ! ... next line useful for band parallelization of exact exchange
      IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
      !
      notcnv = COUNT( .NOT. conv(:) )
      !
      e(1:nvec) = work%ew(1:nvec)
      !
      ! ... if overall convergence has been achieved, or the dimension of
      ! ... the reduced basis set is becoming too large, or in any case if
      ! ... we are at the last iteration refresh the basis set. i.e. replace
      ! ... the first nvec elements with the current estimate of the
      ! ... eigenvectors;  set the basis dimension to nvec.
      !
      IF ( notcnv == 0 .or. dav_iter == maxter ) THEN
        ! ... all roots converged: return
        ! or
        ! ... last iteration, some roots not converged: return
        CALL start_clock( 'cegterg:last' )
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
          work%psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
        !
        CALL stop_clock( 'cegterg:last' )
      END IF
      
      IF ( nbase+notcnv > nvecx) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
          work%psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
        !
        ! ... refresh psi, H*psi and S*psi
        !
        work%psi(:,:,1:nvec) = evc(:,:,1:nvec)
        !
        IF ( uspp ) THEN
          !
          CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, work%spsi, &
            kdmx, vc, nvecx, ZERO, work%psi(1,1,nvec+1), kdmx )
          !
          work%spsi(:,:,1:nvec) = work%psi(:,:,nvec+1:nvec+nvec)
          !
        END IF
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, work%hpsi, &
          kdmx, vc, nvecx, ZERO, work%psi(1,1,nvec+1), kdmx )
        !
        work%hpsi(:,:,1:nvec) = work%psi(:,:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
        !
        DO n = 1, nbase
          !
          !           hc(n,n) = REAL( e(n) )
          hc(n,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
          !
          sc(n,n) = ONE
          vc(n,n) = ONE
          !
        END DO
        !
        CALL stop_clock( 'cegterg:last' )
        !
      END IF

      IF ( notcnv == 0 .or. nbase+notcnv > nvecx .or. dav_iter == maxter ) THEN              
        !
        DEALLOCATE( conv )
        DEALLOCATE( vc )
        DEALLOCATE( hc )
        DEALLOCATE( sc )
        !
        CALL stop_clock( 'cegterg' )
        task = ESL_TASK_EXIT
        step = 1
      ELSE
        ! go back to step 3 without doing anything
        task = ESL_TASK_NONE
        step = 3
      END IF
        
    end select

  end subroutine david_rci_run
  
end module david_rci_m
