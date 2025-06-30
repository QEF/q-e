!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE laxlib_end()
  use laxlib_processors_grid
  CALL laxlib_end_drv ( )
END SUBROUTINE laxlib_end


SUBROUTINE laxlib_getval_ ( nproc_ortho, leg_ortho, np_ortho, me_ortho, ortho_comm, ortho_row_comm, ortho_col_comm, &
  ortho_comm_id, ortho_parent_comm, ortho_cntx, do_distr_diag_inside_bgrp  )
  use laxlib_processors_grid, ONLY : &
    nproc_ortho_ => nproc_ortho, &
    leg_ortho_   => leg_ortho, &
    np_ortho_    => np_ortho, &
    me_ortho_    => me_ortho, &
    ortho_comm_  => ortho_comm, & 
    ortho_row_comm_ => ortho_row_comm, &
    ortho_col_comm_ => ortho_col_comm, & 
    ortho_comm_id_  => ortho_comm_id, &
    ortho_parent_comm_ => ortho_parent_comm, &
    ortho_cntx_  => ortho_cntx, &
    do_distr_diag_inside_bgrp_ => do_distr_diag_inside_bgrp
  IMPLICIT NONE
  INTEGER, OPTIONAL, INTENT(OUT) :: nproc_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: leg_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: np_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: me_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_row_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_col_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm_id
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_parent_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_cntx
  LOGICAL, OPTIONAL, INTENT(OUT) :: do_distr_diag_inside_bgrp
  IF( PRESENT(nproc_ortho) ) nproc_ortho = nproc_ortho_
  IF( PRESENT(leg_ortho) ) leg_ortho = leg_ortho_
  IF( PRESENT(np_ortho) ) np_ortho = np_ortho_
  IF( PRESENT(me_ortho) ) me_ortho = me_ortho_
  IF( PRESENT(ortho_comm) ) ortho_comm = ortho_comm_
  IF( PRESENT(ortho_row_comm) ) ortho_row_comm = ortho_row_comm_
  IF( PRESENT(ortho_col_comm) ) ortho_col_comm = ortho_col_comm_
  IF( PRESENT(ortho_comm_id) ) ortho_comm_id = ortho_comm_id_
  IF( PRESENT(ortho_parent_comm) ) ortho_parent_comm = ortho_parent_comm_
  IF( PRESENT(ortho_cntx) ) ortho_cntx = ortho_cntx_
  IF( PRESENT(do_distr_diag_inside_bgrp) ) do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp_
END SUBROUTINE
!
SUBROUTINE laxlib_get_status_x ( lax_status )
  use laxlib_processors_grid, ONLY : &
    nproc_ortho_ => nproc_ortho, &
    leg_ortho_   => leg_ortho, &
    np_ortho_    => np_ortho, &
    me_ortho_    => me_ortho, &
    ortho_comm_  => ortho_comm, & 
    ortho_row_comm_ => ortho_row_comm, &
    ortho_col_comm_ => ortho_col_comm, & 
    ortho_comm_id_  => ortho_comm_id, &
    ortho_parent_comm_ => ortho_parent_comm, &
    ortho_cntx_  => ortho_cntx, &
    do_distr_diag_inside_bgrp_ => do_distr_diag_inside_bgrp
  IMPLICIT NONE
  include 'laxlib_param.fh'
  INTEGER, INTENT(OUT) :: LAX_STATUS(:)
  lax_status(LAX_STATUS_NPROC)= nproc_ortho_
  lax_status(LAX_STATUS_LEG)= leg_ortho_
  lax_status(LAX_STATUS_NP1)= np_ortho_( 1 )
  lax_status(LAX_STATUS_NP2)= np_ortho_( 2 )
  lax_status(LAX_STATUS_ME1)= me_ortho_( 1 )
  lax_status(LAX_STATUS_ME2)= me_ortho_( 2 )
  lax_status(LAX_STATUS_COMM)= ortho_comm_
  lax_status(LAX_STATUS_ROWCOMM)= ortho_row_comm_
  lax_status(LAX_STATUS_COLCOMM)= ortho_col_comm_
  lax_status(LAX_STATUS_COMMID)= ortho_comm_id_
  lax_status(LAX_STATUS_PARENTCOMM)= ortho_parent_comm_
  lax_status(LAX_STATUS_ORTHOCNTX)= ortho_cntx_
  IF( do_distr_diag_inside_bgrp_ ) THEN
     lax_status(LAX_STATUS_DISTDIAG)= 1
  ELSE
     lax_status(LAX_STATUS_DISTDIAG)= 2
  END IF
END SUBROUTINE

!----------------------------------------------------------------------------

SUBROUTINE laxlib_start_drv( ndiag_, parent_comm, do_distr_diag_inside_bgrp_  )
    !
    use laxlib_processors_grid
    USE laxlib_parallel_include
    !
    !
    ! ... Ortho/diag/linear algebra group initialization
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT) :: ndiag_  ! (IN) input number of procs in the diag group, (OUT) actual number
    INTEGER, INTENT(IN) :: parent_comm ! parallel communicator inside which the distributed linear algebra group
                                       ! communicators are created
    LOGICAL, INTENT(IN) :: do_distr_diag_inside_bgrp_  ! comme son nom l'indique
    !
    INTEGER :: nproc_ortho_try
    INTEGER :: parent_nproc ! nproc of the parent group
    INTEGER :: my_parent_id ! id of the parent communicator 
    INTEGER :: ierr = 0
    !
    IF( lax_is_initialized ) &
       CALL laxlib_end_drv ( )

    parent_nproc = laxlib_size( parent_comm ) ! the number of processors in the current parent communicator
    my_parent_id = laxlib_rank( parent_comm ) ! set the index of the current parent communicator

    ! save input value inside the module
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp_ 

    !
    IF( ndiag_ > 0 ) THEN
       ! The caller suggested a diag group size
       ! Ensuring that it falls in the proper range
       nproc_ortho_try = MIN( ndiag_ , parent_nproc )
    ELSE 
       ! The caller didn't suggest a diag group size
       ! insert here custom architecture specific default definitions
#if defined(__SCALAPACK) && !defined(__CUDA)
       nproc_ortho_try = MAX( parent_nproc, 1 )
#else
       nproc_ortho_try = 1
#endif
    END IF
    !
    ! the ortho group for parallel linear algebra is a sub-group of the pool,
    ! then there are as many ortho groups as pools.
    !
    CALL init_ortho_group ( nproc_ortho_try, parent_comm )
    !
    ! set the number of processors in the diag group to the actual number used
    !
    ndiag_ = nproc_ortho
    !
    lax_is_initialized = .true.
    !  
    RETURN
    !
CONTAINS

  SUBROUTINE init_ortho_group ( nproc_try_in, comm_all )
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nproc_try_in, comm_all

    INTEGER :: ierr, color, key, me_all, nproc_all, nproc_try

#if defined __SCALAPACK
    INTEGER, ALLOCATABLE :: blacsmap(:,:)
    INTEGER, ALLOCATABLE :: ortho_cntx_pe(:)
    INTEGER :: nprow, npcol, myrow, mycol, i, j, k
    INTEGER, EXTERNAL :: BLACS_PNUM
#endif

#if defined __MPI

    me_all    = laxlib_rank( comm_all )
    nproc_all = laxlib_size( comm_all )
    !
    nproc_try = MIN( nproc_try_in, nproc_all )
    nproc_try = MAX( nproc_try, 1 )

    !  find the square closer (but lower) to nproc_try
    !
    CALL grid2d_dims( 'S', nproc_try, np_ortho(1), np_ortho(2) )
    !
    !  now, and only now, it is possible to define the number of tasks
    !  in the ortho group for parallel linear algebra
    !
    nproc_ortho = np_ortho(1) * np_ortho(2)
    !
    IF( nproc_all >= 4*nproc_ortho ) THEN
       !
       !  here we choose a processor every 4, in order not to stress memory BW
       !  on multi core procs, for which further performance enhancements are
       !  possible using OpenMP BLAS inside regter/cegter/rdiaghg/cdiaghg
       !  (to be implemented)
       !
       color = 0
       IF( me_all < 4*nproc_ortho .AND. MOD( me_all, 4 ) == 0 ) color = 1
       !
       leg_ortho = 4
       !
    ELSE IF( nproc_all >= 2*nproc_ortho ) THEN
       !
       !  here we choose a processor every 2, in order not to stress memory BW
       !
       color = 0
       IF( me_all < 2*nproc_ortho .AND. MOD( me_all, 2 ) == 0 ) color = 1
       !
       leg_ortho = 2
       !
    ELSE
       !
       !  here we choose the first processors
       !
       color = 0
       IF( me_all < nproc_ortho ) color = 1
       !
       leg_ortho = 1
       !
    END IF
    !
    key   = me_all
    !
    !  initialize the communicator for the new group by splitting the input communicator
    !
    CALL laxlib_comm_split ( comm_all, color, key, ortho_comm )
    !
    ! and remember where it comes from
    !
    ortho_parent_comm = comm_all
    !
    !  Computes coordinates of the processors, in row maior order
    !
    me_ortho1   = laxlib_rank( ortho_comm )
    !
    IF( me_all == 0 .AND. me_ortho1 /= 0 ) &
         CALL lax_error__( " init_ortho_group ", " wrong root task in ortho group ", ierr )
    !
    if( color == 1 ) then
       ortho_comm_id = 1
       CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
       CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
       IF( ierr /= me_ortho1 ) &
            CALL lax_error__( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
       IF( me_ortho1*leg_ortho /= me_all ) &
            CALL lax_error__( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )

       CALL laxlib_comm_split( ortho_comm, me_ortho(2), me_ortho(1), ortho_col_comm)
       CALL laxlib_comm_split( ortho_comm, me_ortho(1), me_ortho(2), ortho_row_comm)

#if defined __SCALAPACK
       !
       ortho_cntx = ortho_comm
       ! ortho_cntx is both an input and output. input is a system context. output is a BLACS context.
       ! In Fortran, a system context is just a MPI communicator. To avoid any unintended behavior,
       ! make sure ortho_comm matches exactly the BLACS processor grid.
       call BLACS_GRIDINIT(ortho_cntx, 'R', np_ortho(1), np_ortho(2))
#endif

    else
       ortho_comm_id = 0
       me_ortho(1) = me_ortho1
       me_ortho(2) = me_ortho1
    endif

#else

    ortho_comm_id = 1

#endif

    RETURN
  END SUBROUTINE init_ortho_group

END SUBROUTINE laxlib_start_drv

!------------------------------------------------------------------------------!

SUBROUTINE print_lambda_x( lambda, idesc, n, nshow, nudx, ccc, ionode, iunit )
    IMPLICIT NONE
    include 'laxlib_low.fh'
    include 'laxlib_kinds.fh'
    real(DP), intent(in) :: lambda(:,:,:), ccc
    INTEGER, INTENT(IN) :: idesc(:,:)
    integer, intent(in) :: n, nshow, nudx
    logical, intent(in) :: ionode
    integer, intent(in) :: iunit
    !
    integer :: nnn, j, i, is
    real(DP), allocatable :: lambda_repl(:,:)
    nnn = min( nudx, nshow )
    ALLOCATE( lambda_repl( nudx, nudx ) )
    IF( ionode ) WRITE( iunit,*)
    DO is = 1, SIZE( lambda, 3 )
       CALL collect_lambda( lambda_repl, lambda(:,:,is), idesc(:,is) )
       IF( ionode ) THEN
          WRITE( iunit,3370) '    lambda   nudx, spin = ', nudx, is
          IF( nnn < n ) WRITE( iunit,3370) '    print only first ', nnn
          DO i=1,nnn
             WRITE( iunit,3380) (lambda_repl(i,j)*ccc,j=1,nnn)
          END DO
       END IF
    END DO
    DEALLOCATE( lambda_repl )
3370   FORMAT(26x,a,2i4)
3380   FORMAT(9f8.4)
    RETURN
END SUBROUTINE print_lambda_x

  SUBROUTINE laxlib_desc_init1( nsiz, nx, la_proc, idesc, rank_ip, idesc_ip )
     !
     IMPLICIT NONE
     include 'laxlib_low.fh'
     include 'laxlib_param.fh'
     include 'laxlib_kinds.fh'
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: nx
     LOGICAL, INTENT(OUT) :: la_proc
     INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
     INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip( :, : )
     INTEGER, INTENT(OUT), ALLOCATABLE :: idesc_ip(:,:,:)
     !
     INTEGER :: ortho_comm, np_ortho(2), me_ortho(2), ortho_comm_id, &
          leg_ortho
     !
     CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, &
          ortho_comm = ortho_comm, leg_ortho = leg_ortho, &
          ortho_comm_id = ortho_comm_id )
     !
     IF ( .NOT. ALLOCATED (idesc_ip) ) THEN
        ALLOCATE( idesc_ip( LAX_DESC_SIZE, np_ortho(1), np_ortho(2) ) )
     ELSE
        IF ( SIZE (idesc_ip,2) /= np_ortho(1) .OR. &
             SIZE (idesc_ip,3) /= np_ortho(2) ) &
             CALL lax_error__( " desc_init ", " inconsistent dimension ", 2 )
     END IF
     IF ( .NOT. ALLOCATED (rank_ip) ) &
          ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
     !
     CALL laxlib_init_desc( idesc, idesc_ip, rank_ip, nsiz, nsiz )
     ! 
     nx = idesc(LAX_DESC_NRCX)
     !
     la_proc = .FALSE.
     IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) la_proc = .TRUE.
     !
     RETURN
   END SUBROUTINE laxlib_desc_init1
   !
   SUBROUTINE laxlib_desc_init2( nsiz, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
     !
     IMPLICIT NONE
     include 'laxlib_low.fh'
     include 'laxlib_param.fh'
     include 'laxlib_kinds.fh'
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: nx
     LOGICAL, INTENT(OUT) :: la_proc
     INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
     INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip(:,:)
     INTEGER, INTENT(OUT), ALLOCATABLE :: irc_ip(:)
     INTEGER, INTENT(OUT), ALLOCATABLE :: nrc_ip(:)

     INTEGER :: i, j, rank
     INTEGER :: ortho_comm, np_ortho(2), me_ortho(2), ortho_comm_id, &
          leg_ortho, ortho_cntx
     !
     CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, &
          ortho_comm = ortho_comm, leg_ortho = leg_ortho, &
          ortho_comm_id = ortho_comm_id, ortho_cntx = ortho_cntx )
     !
     CALL laxlib_init_desc( idesc, nsiz, nsiz, np_ortho, me_ortho, &
          ortho_comm, ortho_cntx, ortho_comm_id )
     !
     nx = idesc(LAX_DESC_NRCX)
     !
     IF ( .NOT. ALLOCATED (rank_ip) ) THEN
        ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
        ALLOCATE( irc_ip( np_ortho(1) ), nrc_ip (np_ortho(1) ) )
     ELSE
        IF ( SIZE (rank_ip,1) /= np_ortho(1) .OR. &
             SIZE (rank_ip,2) /= np_ortho(2) ) &
             CALL lax_error__( " desc_init ", " inconsistent dimension ", 1 )
     END IF
     DO j = 0, idesc(LAX_DESC_NPC) - 1
        CALL laxlib_local_dims( irc_ip( j + 1 ), nrc_ip( j + 1 ), &
             idesc(LAX_DESC_N), idesc(LAX_DESC_NX), np_ortho(1), j )
        DO i = 0, idesc(LAX_DESC_NPR) - 1
           CALL GRID2D_RANK( 'R', idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) la_proc = .TRUE.
     !
     RETURN
   END SUBROUTINE laxlib_desc_init2
  !
  !
SUBROUTINE laxlib_init_desc_x( idesc, n, nx, np, me, comm, cntx, comm_id )
    USE laxlib_descriptor,       ONLY: la_descriptor, descla_init, laxlib_desc_to_intarray
    IMPLICIT NONE
    include 'laxlib_param.fh'
    INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
    INTEGER, INTENT(IN)  :: n   !  the size of this matrix
    INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution
    INTEGER, INTENT(IN)  :: np(2), me(2), comm, cntx
    INTEGER, INTENT(IN)  :: comm_id
    !
    TYPE(la_descriptor) :: descla
    !
    CALL descla_init( descla, n, nx, np, me, comm, cntx, comm_id )
    CALL laxlib_desc_to_intarray( idesc, descla )
    RETURN
END SUBROUTINE laxlib_init_desc_x


SUBROUTINE laxlib_multi_init_desc_x( idesc, idesc_ip, rank_ip, n, nx  )
    USE laxlib_descriptor,       ONLY: la_descriptor, descla_init, laxlib_desc_to_intarray
    use laxlib_processors_grid,  ONLY: leg_ortho, np_ortho, me_ortho, ortho_comm, ortho_comm_id, ortho_cntx
    IMPLICIT NONE
    include 'laxlib_param.fh'
    INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
    INTEGER, INTENT(OUT) :: idesc_ip(:,:,:)
    INTEGER, INTENT(OUT) :: rank_ip(:,:)
    INTEGER, INTENT(IN)  :: n   !  the size of this matrix
    INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution

    INTEGER :: i, j, rank, includeme
    INTEGER :: coor_ip( 2 )
    !
    TYPE(la_descriptor) :: descla
    !
    CALL descla_init( descla, n, nx, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
    !
    CALL laxlib_desc_to_intarray( idesc, descla )
    !
    includeme = 1
    !
    DO j = 0, idesc(LAX_DESC_NPC) - 1
       DO i = 0, idesc(LAX_DESC_NPR) - 1
          coor_ip( 1 ) = i
          coor_ip( 2 ) = j
          CALL descla_init( descla, idesc(LAX_DESC_N), idesc(LAX_DESC_NX), &
               np_ortho, coor_ip, ortho_comm, ortho_cntx, includeme )
          CALL laxlib_desc_to_intarray( idesc_ip(:,i+1,j+1), descla )
          CALL GRID2D_RANK( 'R', idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC), i, j, rank )
          rank_ip( i+1, j+1 ) = rank * leg_ortho
       END DO
    END DO
    !
    RETURN
END SUBROUTINE laxlib_multi_init_desc_x

   SUBROUTINE descla_local_dims( i2g, nl, n, nx, np, me )
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: i2g  !  global index of the first local element
      INTEGER, INTENT(OUT) :: nl   !  local number of elements
      INTEGER, INTENT(IN)  :: n    !  number of actual element in the global array
      INTEGER, INTENT(IN)  :: nx   !  dimension of the global array (nx>=n) to be distributed
      INTEGER, INTENT(IN)  :: np   !  number of processors
      INTEGER, INTENT(IN)  :: me   !  taskid for which i2g and nl are computed
      !
      !  note that we can distribute a global array larger than the
      !  number of actual elements. This could be required for performance
      !  reasons, and to have an equal partition of matrix having different size
      !  like matrixes of spin-up and spin-down
      !
      INTEGER, EXTERNAL ::  ldim_block, ldim_cyclic, ldim_block_sca
      INTEGER, EXTERNAL ::  gind_block, gind_cyclic, gind_block_sca
      !
#if __SCALAPACK
      nl  = ldim_block_sca( nx, np, me )
      i2g = gind_block_sca( 1, nx, np, me )
#else
      nl  = ldim_block( nx, np, me )
      i2g = gind_block( 1, nx, np, me )
#endif
      ! This is to try to keep a matrix N * N into the same
      ! distribution of a matrix NX * NX, useful to have
      ! the matrix of spin-up distributed in the same way
      ! of the matrix of spin-down
      !
      IF( i2g + nl - 1 > n ) nl = n - i2g + 1
      IF( nl < 0 ) nl = 0
      RETURN
      !
   END SUBROUTINE descla_local_dims


!   ----------------------------------------------
!   Simplified driver 

   SUBROUTINE diagonalize_parallel_x( n, rhos, rhod, s, idesc )

      USE dspev_module

      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      include 'laxlib_param.fh'
      include 'laxlib_mid.fh'
      include 'laxlib_low.fh'
      REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
      REAL(DP)              :: rhod(:)   !  output eigenvalues
      REAL(DP)              :: s(:,:)    !  output eigenvectors
      INTEGER,  INTENT(IN) :: n         !  size of the global matrix
      INTEGER,  INTENT(IN) :: idesc(LAX_DESC_SIZE)

      IF( n < 1 ) RETURN

      !  Matrix is distributed on the same processors group
      !  used for parallel matrix multiplication
      !
      IF( SIZE(s,1) /= SIZE(rhos,1) .OR. SIZE(s,2) /= SIZE(rhos,2) ) &
         CALL lax_error__( " diagonalize_parallel ", " inconsistent dimension for s and rhos ", 1 )

      IF ( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         IF( SIZE(s,1) /= idesc(LAX_DESC_NRCX) ) &
            CALL lax_error__( " diagonalize_parallel ", " inconsistent dimension ", 1)
         !
         !  Compute local dimension of the cyclically distributed matrix
         !
         s = rhos
         !
#if defined(__SCALAPACK)
         CALL pdsyevd_drv( .true. , n, idesc(LAX_DESC_NRCX), s, SIZE(s,1), rhod, idesc(LAX_DESC_CNTX), idesc(LAX_DESC_COMM) )
#else
         CALL laxlib_pdsyevd( .true., n, idesc, s, SIZE(s,1), rhod )
#endif
         !
      END IF

      RETURN

   END SUBROUTINE diagonalize_parallel_x


   SUBROUTINE diagonalize_serial_x( n, rhos, rhod )
      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      include 'laxlib_low.fh'
      INTEGER,  INTENT(IN)  :: n
      REAL(DP)              :: rhos(:,:)
      REAL(DP)              :: rhod(:)
      !
      ! inputs:
      ! n     size of the eigenproblem
      ! rhos  the symmetric matrix
      ! outputs:
      ! rhos  eigenvectors
      ! rhod  eigenvalues
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER :: i, j, k

      IF( n < 1 ) RETURN

      ALLOCATE( aux( n * ( n + 1 ) / 2 ) )

      !  pack lower triangle of rho into aux
      !
      k = 0
      DO j = 1, n
         DO i = j, n
            k = k + 1
            aux( k ) = rhos( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, aux, rhod, rhos, SIZE(rhos,1) )

      DEALLOCATE( aux )

      RETURN
   END SUBROUTINE diagonalize_serial_x

   SUBROUTINE diagonalize_serial_gpu( m, rhos, rhod, s, info )
#if defined(__CUDA)
      use cudafor
      USE cusolverDn
      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      INTEGER, INTENT(IN) :: m
      REAL(DP), DEVICE, INTENT(IN) :: rhos(:,:)
      REAL(DP), DEVICE, INTENT(OUT) :: rhod(:)
      REAL(DP), DEVICE, INTENT(OUT) :: s(:,:)
      INTEGER, INTENT(OUT) :: info
      !
      INTEGER :: lwork_d
      INTEGER :: i, j, lda
      !
      !
      INTEGER, DEVICE        :: devInfo
      TYPE(cusolverDnHandle) :: cuSolverHandle
      REAL(DP), ALLOCATABLE, DEVICE :: work_d(:)
      !
      ! .... Subroutine Body
      !
      !
      s = rhos
      lda = SIZE( rhos, 1 )
      !
      info = cusolverDnCreate(cuSolverHandle)
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) &
         CALL lax_error__( ' diagonalize_serial_gpu ', 'cusolverDnCreate',  ABS( info ) )

      info = cusolverDnDsyevd_bufferSize( &
             cuSolverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, m, s, lda, rhod, lwork_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' error in solver 1 ', ABS( info ) )

      ALLOCATE( work_d ( lwork_d ), STAT=info )
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate work_d ', ABS( info ) )

      info = cusolverDnDsyevd( &
             cuSolverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, m, s, lda, rhod, work_d, lwork_d, devInfo)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' error in solver 2 ', ABS( info ) )

      info = cudaDeviceSynchronize()

      info = cusolverDnDestroy(cuSolverHandle)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' diagonalize_serial_gpu ', ' cusolverDnDestroy failed ', ABS( info ) )

      DEALLOCATE( work_d )

      !
#else
      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      INTEGER  :: m
      REAL(DP) :: rhos(:,:)
      REAL(DP) :: rhod(:)
      REAL(DP) :: s(:,:)
      INTEGER  :: info
      !
      CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' not compiled in this version ', 0 )
#endif
   END SUBROUTINE

