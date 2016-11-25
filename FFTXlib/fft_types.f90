!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
MODULE fft_types
!=----------------------------------------------------------------------------=!

  USE fft_support, ONLY : good_fft_order, good_fft_dimension
  USE fft_param

  IMPLICIT NONE
  PRIVATE
  SAVE

  TYPE fft_type_descriptor

    INTEGER :: nst      ! total number of sticks
    INTEGER, ALLOCATABLE :: nsp(:)   ! number of sticks per processor ( potential )
                                 ! using proc index starting from 1 !!
                                 ! on proc mpime -> nsp( mpime + 1 )
    INTEGER, ALLOCATABLE :: nsw(:)   ! number of sticks per processor ( wave func )
                                 ! using proc index as above
    INTEGER :: nr1    = 0  !
    INTEGER :: nr2    = 0  ! effective FFT dimensions of the 3D grid (global)
    INTEGER :: nr3    = 0  ! 
    INTEGER :: nr1x   = 0  ! FFT grids leading dimensions
    INTEGER :: nr2x   = 0  ! dimensions of the arrays for the 3D grid (global)
    INTEGER :: nr3x   = 0  ! may differ from nr1 ,nr2 ,nr3 in order to boost performances

    INTEGER :: npl    = 0  ! number of "Z" planes for this processor = npp( mpime + 1 )
    INTEGER :: nnp    = 0  ! number of 0 and non 0 sticks in a plane ( ~nr1*nr2/nproc )
    INTEGER :: nnr    = 0  ! local number of FFT grid elements  ( ~nr1*nr2*nr3/proc )
                           ! size of the arrays allocated for the FFT, local to each processor:
                           ! in parallel execution may differ from nr1x*nr2x*nr3x
                           ! Not to be confused either with nr1*nr2*nr3 
    INTEGER, ALLOCATABLE :: ngl(:)   ! per proc. no. of non zero charge density/potential components
    INTEGER, ALLOCATABLE :: nwl(:)   ! per proc. no. of non zero wave function plane components
    INTEGER, ALLOCATABLE :: npp(:)   ! number of "Z" planes per processor
    INTEGER, ALLOCATABLE :: ipp(:)   ! offset of the first "Z" plane on each proc ( 0 on the first proc!!!)
    INTEGER, ALLOCATABLE :: iss(:)   ! index of the first rho stick on each proc
    INTEGER, ALLOCATABLE :: isind(:) ! for each position in the plane indicate the stick index
    INTEGER, ALLOCATABLE :: ismap(:) ! for each stick in the plane indicate the position
    INTEGER, ALLOCATABLE :: iplp(:)  ! indicate which "Y" plane should be FFTed ( potential )
    INTEGER, ALLOCATABLE :: iplw(:)  ! indicate which "Y" plane should be FFTed ( wave func )
    !
    !  fft parallelization
    !
    INTEGER :: mype     = 0          ! my processor id (starting from 0) in the fft group
#if defined(__MPI)
    INTEGER :: comm     = MPI_COMM_NULL
#else
    INTEGER :: comm     = 0          ! communicator of the fft gruop 
#endif
    INTEGER :: nproc    = 1          ! number of processor in the fft group
    INTEGER :: root     = 0          ! root processor
    LOGICAL :: lpara    = .FALSE.
    !
  END TYPE


  REAL(DP) :: fft_dual = 4.0d0

  PUBLIC :: fft_type_descriptor, fft_type_init
  PUBLIC :: fft_type_allocate, fft_type_deallocate

CONTAINS

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_setdim( desc, nr1, nr2, nr3 )
     TYPE (fft_type_descriptor) :: desc
     INTEGER, INTENT(IN) :: nr1, nr2, nr3
     IF (desc%nr1 /= 0 .OR. desc%nr1 /= 0 .OR. desc%nr1 /= 0 ) &
        CALL fftx_error__(' fft_type_setdim ', ' fft dimensions already set ', 1 )
     desc%nr1 = nr1
     desc%nr2 = nr2
     desc%nr3 = nr3
     desc%nr1 = good_fft_order( desc%nr1 )
     desc%nr2 = good_fft_order( desc%nr2 )
     desc%nr3 = good_fft_order( desc%nr3 )
     desc%nr1x  = good_fft_dimension( desc%nr1 )
     desc%nr2x  = desc%nr2
     desc%nr3x  = good_fft_dimension( desc%nr3 )
  END SUBROUTINE

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_allocate( desc, at, bg, gcutm, comm, fft_fact  )
  !
  ! routine that allocate arrays of fft_type_descriptor
  ! must be called before fft_type_set
  !
    TYPE (fft_type_descriptor) :: desc
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
    REAL(DP), INTENT(IN) :: gcutm
    INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
    INTEGER, INTENT(in) :: comm ! mype starting from 0
    INTEGER :: nx, ny, ierr
    INTEGER :: mype, root, nproc ! mype starting from 0

    IF ( ALLOCATED( desc%nsp ) ) &
        CALL fftx_error__(' fft_type_allocate ', ' fft arrays already allocated ', 1 )

    desc%comm = comm 

#if defined(__MPI)
    IF( desc%comm == MPI_COMM_NULL ) THEN
       CALL fftx_error__( ' realspace_grid_init ', ' fft communicator is null ', 1 )
    END IF
#endif
       !
    mype = 0
    nproc = 1
    root = 0
#if defined(__MPI)
    CALL MPI_COMM_RANK( comm, mype, ierr )
    CALL MPI_COMM_SIZE( comm, nproc, ierr )
#endif

    CALL realspace_grid_init( desc, at, bg, gcutm, fft_fact )

    nx = desc%nr1x 
    ny = desc%nr2x

    ALLOCATE( desc%nsp( nproc ) )
    ALLOCATE( desc%nsw( nproc ) )
    ALLOCATE( desc%ngl( nproc ) )
    ALLOCATE( desc%nwl( nproc ) )
    ALLOCATE( desc%npp( nproc ) )
    ALLOCATE( desc%ipp( nproc ) )
    ALLOCATE( desc%iss( nproc ) )
    ALLOCATE( desc%isind( nx * ny ) )
    ALLOCATE( desc%ismap( nx * ny ) )
    ALLOCATE( desc%iplp( nx ) )
    ALLOCATE( desc%iplw( nx ) )

    desc%nsp   = 0
    desc%nsw   = 0
    desc%ngl   = 0
    desc%nwl   = 0
    desc%npp   = 0
    desc%ipp   = 0
    desc%iss   = 0
    desc%isind = 0
    desc%ismap = 0
    desc%iplp  = 0
    desc%iplw  = 0

    desc%mype  = mype
    desc%comm  = comm
    desc%nproc = nproc
    desc%root  = root

  END SUBROUTINE fft_type_allocate

  SUBROUTINE fft_type_deallocate( desc )
    TYPE (fft_type_descriptor) :: desc
    IF ( ALLOCATED( desc%nsp ) )    DEALLOCATE( desc%nsp )
    IF ( ALLOCATED( desc%nsw ) )    DEALLOCATE( desc%nsw )
    IF ( ALLOCATED( desc%ngl ) )    DEALLOCATE( desc%ngl )
    IF ( ALLOCATED( desc%nwl ) )    DEALLOCATE( desc%nwl )
    IF ( ALLOCATED( desc%npp ) )    DEALLOCATE( desc%npp )
    IF ( ALLOCATED( desc%ipp ) )    DEALLOCATE( desc%ipp )
    IF ( ALLOCATED( desc%iss ) )    DEALLOCATE( desc%iss )
    IF ( ALLOCATED( desc%isind ) )  DEALLOCATE( desc%isind )
    IF ( ALLOCATED( desc%ismap ) )  DEALLOCATE( desc%ismap )
    IF ( ALLOCATED( desc%iplp ) )   DEALLOCATE( desc%iplp )
    IF ( ALLOCATED( desc%iplw ) )   DEALLOCATE( desc%iplw )
#if defined(__MPI)
    desc%comm = MPI_COMM_NULL 
#endif
    desc%nr1    = 0 
    desc%nr2    = 0  
    desc%nr3    = 0  
    desc%nr1x   = 0  
    desc%nr2x   = 0  
    desc%nr3x   = 0  
  END SUBROUTINE fft_type_deallocate

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_set( desc, tk, lpara, ntg, nst, ub, lb, idx, in1, in2, ncp, ncpw, ngp, ngpw, st, stw )

    TYPE (fft_type_descriptor) :: desc

    LOGICAL, INTENT(in) :: tk               ! gamma/not-gamma logical
    LOGICAL, INTENT(in) :: lpara            ! set array for parallel or serial FFT drivers
    INTEGER, INTENT(in) :: ntg              ! number of task groups (optimal spacing for residual-plane distribution)
    INTEGER, INTENT(in) :: nst              ! total number of stiks 
    INTEGER, INTENT(in) :: ub(3), lb(3)     ! upper and lower bound of real space indices
    INTEGER, INTENT(in) :: idx(:)           ! sorting index of the sticks
    INTEGER, INTENT(in) :: in1(:)           ! x-index of a stick
    INTEGER, INTENT(in) :: in2(:)           ! y-index of a stick
    INTEGER, INTENT(in) :: ncp(:)           ! number of rho  columns per processor
    INTEGER, INTENT(in) :: ncpw(:)          ! number of wave columns per processor
    INTEGER, INTENT(in) :: ngp(:)           ! number of rho  G-vectors per processor
    INTEGER, INTENT(in) :: ngpw(:)          ! number of wave G-vectors per processor
    INTEGER, INTENT(in) :: st( lb(1) : ub(1), lb(2) : ub(2) )   ! stick owner of a given rho stick
    INTEGER, INTENT(in) :: stw( lb(1) : ub(1), lb(2) : ub(2) )  ! stick owner of a given wave stick

    INTEGER :: npp( desc%nproc ), n3( desc%nproc ), nsp( desc%nproc )
    INTEGER :: np, nq, i, is, iss, itg, i1, i2, m1, m2, n1, n2, ip
    INTEGER :: ncpx, nppx
    INTEGER :: nr1, nr2, nr3    ! size of real space grid 
    INTEGER :: nr1x, nr2x, nr3x ! padded size of real space grid
    !
    IF (.NOT. ALLOCATED( desc%nsp ) ) &
        CALL fftx_error__(' fft_type_allocate ', ' fft arrays not yet allocated ', 1 )

    IF ( desc%nr1 == 0 .OR. desc%nr2 == 0 .OR. desc%nr3 == 0 ) &
        CALL fftx_error__(' fft_type_set ', ' fft dimensions not yet set ', 1 )

    !  Set fft actual and leading dimensions to be used internally

    nr1  = desc%nr1
    nr2  = desc%nr2
    nr3  = desc%nr3
    nr1x = desc%nr1x
    nr2x = desc%nr2x
    nr3x = desc%nr3x

    IF( ( nr1 > nr1x ) .or. ( nr2 > nr2x ) .or. ( nr3 > nr3x ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong fft dimensions ', 1 )

    IF( ( size( desc%ngl ) < desc%nproc ) .or. ( size( desc%npp ) < desc%nproc ) .or.  &
        ( size( desc%ipp ) < desc%nproc ) .or. ( size( desc%iss ) < desc%nproc ) )     &
      CALL fftx_error__( ' fft_type_set ', ' wrong descriptor dimensions ', 2 )

    IF( ( size( idx ) < nst ) .or. ( size( in1 ) < nst ) .or. ( size( in2 ) < nst ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong number of stick dimensions ', 3 )

    IF( ( size( ncp ) < desc%nproc ) .or. ( size( ngp ) < desc%nproc ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong stick dimensions ', 4 )

    !  Set the number of "xy" planes for each processor
    !  in other word do a slab partition along the z axis

    npp = 0
    IF ( desc%nproc == 1 ) THEN   ! sigle processor: npp(1)=nr3
      npp(1) = nr3
    ELSE
      np = nr3 / desc%nproc
      nq = nr3 - np * desc%nproc
      npp(1:desc%nproc) = np      ! assign a base value to all processors
      IF( np == 0 ) THEN
        DO i = 1, desc%nproc
         IF ( i <= nq ) npp(i) = 1
        ENDDO
      ELSE
        reminder_loop : &           ! assign an extra plane to processors so that they are spaced by ntg
        DO itg = 1, ntg  
          DO i = itg, desc%nproc, ntg
             IF (nq==0) EXIT reminder_loop
             nq = nq - 1
             npp(i) = np + 1
          ENDDO
        ENDDO reminder_loop
      END IF
    ENDIF

    !-- npp(1:nproc) is the number of planes per processor
    !-- npl is the number of planes per processor of this processor

    desc%npp( 1:desc%nproc )  = npp    
    desc%npl = npp( desc%mype + 1 )    

    !  Find out the index of the starting plane on each proc

    n3 = 0
    DO i = 2, desc%nproc
      n3(i) = n3(i-1) + npp(i-1)
    ENDDO

    !-- ipp(1:proc) is the index-offset of the starting plane of each processor 

    desc%ipp( 1:desc%nproc )  = n3

    ! dimension of the xy plane. see ncplane

    desc%nnp  = nr1x * nr2x  

    !  Set fft local workspace dimension

    nppx = 0
    ncpx = 0
    DO i = 1, desc%nproc
       nppx = MAX( nppx, npp( i ) )  ! maximum number of planes per processor
       ncpx = MAX( ncpx, ncp( i ) )  ! maximum number of columns per processor 
    END DO

    IF ( desc%nproc == 1 ) THEN
      desc%nnr  = nr1x * nr2x * nr3x
    ELSE
      desc%nnr  = max( nr3x * ncpx, nr1x * nr2x * nppx )  ! this is required to contain the local data in R and G space
      desc%nnr  = max( desc%nnr, ncpx * nppx * desc%nproc )  ! this is required to use ALLTOALL instead of ALLTOALLV
      desc%nnr  = max( 1, desc%nnr ) ! ensure that desc%nrr > 0 ( for extreme parallelism )
    ENDIF

    desc%ngl( 1:desc%nproc )  = ngp( 1:desc%nproc )  ! local number of g vectors (rho) per processor
    desc%nwl( 1:desc%nproc )  = ngpw( 1:desc%nproc ) ! local number of g vectors (wave) per processor

    IF( size( desc%isind ) < ( nr1x * nr2x ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong descriptor dimensions, isind ', 5 )

    IF( size( desc%iplp ) < ( nr1x ) .or. size( desc%iplw ) < ( nr1x ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong descriptor dimensions, ipl ', 5 )

    !
    !  1. Temporarily store in the array "desc%isind" the index of the processor
    !     that own the corresponding stick (index of proc starting from 1)
    !  2. Set the array elements of  "desc%iplw" and "desc%iplp" to one
    !     for that index corresponding to YZ planes containing at least one stick
    !     this are used in the FFT transform along Y
    !

    desc%isind = 0    ! will contain the +ve or -ve of the processor number, if any, that owns the stick
    desc%iplp  = 0    ! 1 if the given x-plane location is active in rho  y-fft
    desc%iplw  = 0    ! 1 if the given x-plane location is active in wave y-fft

    !  Set nst the proper number of sticks (the total number of 1d fft along z to be done)
    !
    desc%nst = 0
    DO iss = 1, SIZE( idx )
      is = idx( iss )
      IF( is < 1 ) CYCLE
      i1 = in1( is )
      i2 = in2( is )
      IF( st( i1, i2 ) > 0 ) THEN
        desc%nst = desc%nst + 1
        m1 = i1 + 1; IF ( m1 < 1 ) m1 = m1 + nr1
        m2 = i2 + 1; IF ( m2 < 1 ) m2 = m2 + nr2
        IF( stw( i1, i2 ) > 0 ) THEN
          desc%isind( m1 + ( m2 - 1 ) * nr1x ) =  st( i1, i2 )
          desc%iplw( m1 ) = 1
        ELSE
          desc%isind( m1 + ( m2 - 1 ) * nr1x ) = -st( i1, i2 )
        ENDIF
        desc%iplp( m1 ) = 1
        IF( .not. tk ) THEN
          IF( i1 /= 0 .OR. i2 /= 0 ) desc%nst = desc%nst + 1
          n1 = -i1 + 1; IF ( n1 < 1 ) n1 = n1 + nr1
          n2 = -i2 + 1; IF ( n2 < 1 ) n2 = n2 + nr2
          IF( stw( -i1, -i2 ) > 0 ) THEN
            desc%isind( n1 + ( n2 - 1 ) * nr1x ) =  st( -i1, -i2 )
            desc%iplw( n1 ) = 1
          ELSE
            desc%isind( n1 + ( n2 - 1 ) * nr1x ) = -st( -i1, -i2 )
          ENDIF
          desc%iplp( n1 ) = 1
        ENDIF
      ENDIF
    ENDDO

    !
    !  Compute for each proc the global index ( starting from 0 ) of the first
    !  local stick ( desc%iss )
    !

    DO i = 1, desc%nproc
      IF( i == 1 ) THEN
        desc%iss( i ) = 0
      ELSE
        desc%iss( i ) = desc%iss( i - 1 ) + ncp( i - 1 )
      ENDIF
    ENDDO

    ! iss(1:nproc) is the index offset of the first column of a given processor

    IF( size( desc%ismap ) < ( nst ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong descriptor dimensions ', 6 )

    !
    !  1. Set the array desc%ismap which maps stick indexes to
    !     position in the plane  ( iss )
    !  2. Re-set the array "desc%isind",  that maps position
    !     in the plane with stick indexes (it is the inverse of desc%ismap )
    !

    !  wave function sticks first

    desc%ismap = 0     ! will be the global xy stick index in the global list of processor-ordered sticks
    nsp        = 0     ! will be the number of sticks of a given processor
    DO iss = 1, size( desc%isind )
      ip = desc%isind( iss ) ! processor that owns iss wave stick. if it's a rho stick it's negative !
      IF( ip > 0 ) THEN ! only operates on wave sticks
        nsp( ip ) = nsp( ip ) + 1
        desc%ismap( nsp( ip ) + desc%iss( ip ) ) = iss
        IF( ip == ( desc%mype + 1 ) ) THEN
          desc%isind( iss ) = nsp( ip ) ! isind is OVERWRITTEN as the ordered index in this processor stick list
        ELSE
          desc%isind( iss ) = 0         ! zero otherwise...
        ENDIF
      ENDIF
    ENDDO

    !  check number of stick against the input value

    IF( any( nsp( 1:desc%nproc ) /= ncpw( 1:desc%nproc ) ) ) THEN
      DO ip = 1, desc%nproc
        WRITE( stdout,*)  ' * ', ip, ' * ', nsp( ip ), ' /= ', ncpw( ip )
      ENDDO
      CALL fftx_error__( ' fft_type_set ', ' inconsistent number of sticks ', 7 )
    ENDIF

    desc%nsw( 1:desc%nproc ) = nsp( 1:desc%nproc )  ! -- number of wave sticks per porcessor

    !  then add pseudopotential stick

    DO iss = 1, size( desc%isind )
      ip = desc%isind( iss ) ! -ve of processor that owns iss rho stick. if it was a wave stick it's something non negative !
      IF( ip < 0 ) THEN
        nsp( -ip ) = nsp( -ip ) + 1
        desc%ismap( nsp( -ip ) + desc%iss( -ip ) ) = iss
        IF( -ip == ( desc%mype + 1 ) ) THEN
          desc%isind( iss ) = nsp( -ip ) ! isind is OVERWRITTEN as the ordered index in this processor stick list
        ELSE
          desc%isind( iss ) = 0         ! zero otherwise...
        ENDIF
      ENDIF
    ENDDO

    !  check number of stick against the input value

    IF( any( nsp( 1:desc%nproc ) /= ncp( 1:desc%nproc ) ) ) THEN
      DO ip = 1, desc%nproc
        WRITE( stdout,*)  ' * ', ip, ' * ', nsp( ip ), ' /= ', ncp( ip )
      ENDDO
      CALL fftx_error__( ' fft_type_set ', ' inconsistent number of sticks ', 8 )
    ENDIF

    desc%nsp( 1:desc%nproc ) = nsp( 1:desc%nproc ) ! -- number of rho sticks per processor

    IF( .NOT. lpara ) THEN

       desc%isind = 0
       desc%iplw  = 0
       desc%iplp  = 1

       ! here we are setting parameter as if we were
       ! in a serial code, sticks are along X dimension
       ! and not along Z
       desc%nsp(1) = 0
       desc%nsw(1) = 0
       DO i1 = lb( 1 ), ub( 1 )
         DO i2 = lb( 2 ), ub( 2 )
           m1 = i1 + 1; IF ( m1 < 1 ) m1 = m1 + nr1
           m2 = i2 + 1; IF ( m2 < 1 ) m2 = m2 + nr2
           IF( st( i1, i2 ) > 0 ) THEN
             desc%nsp(1) = desc%nsp(1) + 1
           END IF
           IF( stw( i1, i2 ) > 0 ) THEN
             desc%nsw(1) = desc%nsw(1) + 1
             desc%isind( m1 + ( m2 - 1 ) * nr1x ) =  1  ! st( i1, i2 )
             desc%iplw( m1 ) = 1
           ENDIF
         ENDDO
       ENDDO
       !
       ! if we are in a parallel run, but would like to use serial FFT, all
       ! tasks must have the same parameters as if serial run.
       !
       desc%nsw = desc%nsw(1)
       desc%nsp = desc%nsp(1)
       desc%nnr  = nr1x * nr2x * nr3x
       desc%npl  = nr3
       desc%nnp  = nr1x * nr2x
       desc%npp  = nr3
       desc%ipp  = 0
       desc%ngl  = SUM(ngp)  
       desc%nwl  = SUM(ngpw) 
       !
    END IF

    RETURN
  END SUBROUTINE fft_type_set
!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_init( dfft, smap, pers, lgamma, lpara, comm, at, bg, gcut_in, dual_in, ntask_groups )

     USE stick_base

     TYPE (fft_type_descriptor), INTENT(INOUT) :: dfft 
     TYPE (sticks_map), INTENT(INOUT) :: smap
     CHARACTER(LEN=*), INTENT(IN) :: pers ! fft personality
     LOGICAL, INTENT(IN) :: lpara
     LOGICAL, INTENT(IN) :: lgamma
     INTEGER, INTENT(IN) :: comm
     REAL(DP), INTENT(IN) :: gcut_in
     REAL(DP), INTENT(IN) :: bg(3,3)
     REAL(DP), INTENT(IN) :: at(3,3)
     REAL(DP), OPTIONAL, INTENT(IN) :: dual_in
     INTEGER,  OPTIONAL, INTENT(IN) :: ntask_groups
!
!    Potential or dual
!
     INTEGER, ALLOCATABLE :: st(:,:)
! ...   stick map, st(i,j) = number of G-vector in the
! ...   stick whose x and y miller index are i and j
     INTEGER, ALLOCATABLE :: nstp(:)
! ...   number of sticks, nstp(ip) = number of stick for processor ip
     INTEGER, ALLOCATABLE :: sstp(:)
! ...   number of G-vectors, sstp(ip) = sum of the
! ...   sticks length for processor ip = number of G-vectors owned by the processor ip
     INTEGER :: nst
! ...   nst      local number of sticks
!
! ...     Plane wave
!
     INTEGER, ALLOCATABLE :: stw(:,:)
! ...   stick map (wave functions), stw(i,j) = number of G-vector in the
! ...   stick whose x and y miller index are i and j
     INTEGER, ALLOCATABLE :: nstpw(:)
! ...   number of sticks (wave functions), nstpw(ip) = number of stick for processor ip
     INTEGER, ALLOCATABLE :: sstpw(:)
! ...   number of G-vectors (wave functions), sstpw(ip) = sum of the
! ...   sticks length for processor ip = number of G-vectors owned by the processor ip
     INTEGER :: nstw
! ...   nstw     local number of sticks (wave functions)
     INTEGER :: ntg
! ...   ntg       number of task groups (assigned from input if any)


     REAL(DP) :: gcut, gkcut, dual
     INTEGER  :: ngm, ngw

     ntg = 1
     IF( PRESENT( ntask_groups ) ) ntg = ntask_groups
     dual = fft_dual
     IF( PRESENT( dual_in ) ) dual = dual_in

     IF( pers == 'rho' ) THEN
        gcut = gcut_in
        gkcut = gcut / dual
     ELSE IF ( pers == 'wave' ) THEN
        gkcut = gcut_in
        gcut = gkcut * dual
     ELSE
        CALL fftx_error__(' fft_type_init ', ' unknown FFT personality ', 1 )
     END IF

     IF( .NOT. ALLOCATED( dfft%nsp ) ) THEN
        CALL fft_type_allocate( dfft, at, bg, gcut, comm  )
     ELSE
        IF( dfft%comm /= comm ) THEN
           CALL fftx_error__(' fft_type_init ', ' FFT already allocated with a different communicator ', 1 )
        END IF
     END IF

     dfft%lpara = lpara  !  this descriptor can be either a descriptor for a
                         !  parallel FFT or a serial FFT even in parallel build

     CALL sticks_map_allocate( smap, lgamma, lpara, dfft%nr1, dfft%nr2, dfft%nr3, bg, dfft%comm )

     ALLOCATE( stw ( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) )
     ALLOCATE( st  ( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) )
     ALLOCATE( nstp(smap%nproc) )
     ALLOCATE( sstp(smap%nproc) )
     ALLOCATE( nstpw(smap%nproc) )
     ALLOCATE( sstpw(smap%nproc) )

     CALL get_sticks(  smap, gkcut, nstpw, sstpw, stw, nstw, ngw )
     CALL get_sticks(  smap, gcut,  nstp, sstp, st, nst, ngm )

     CALL fft_type_set( dfft, .not.smap%lgamma, lpara, ntg, nst, smap%ub, smap%lb, smap%idx, &
                             smap%ist(:,1), smap%ist(:,2), nstp, nstpw, sstp, sstpw, st, stw )

     DEALLOCATE( st )
     DEALLOCATE( stw )
     DEALLOCATE( nstp )
     DEALLOCATE( sstp )
     DEALLOCATE( nstpw )
     DEALLOCATE( sstpw )


  END SUBROUTINE fft_type_init

!=----------------------------------------------------------------------------=!
!=----------------------------------------------------------------------------=!



     SUBROUTINE realspace_grid_init( dfft, at, bg, gcutm, fft_fact )
       !
       ! ... Sets optimal values for dfft%nr[123] and dfft%nr[123]x
       ! ... If fft_fact is present, force nr[123] to be multiple of fft_fac([123])
       !
       USE fft_support, only: good_fft_dimension, good_fft_order
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
       REAL(DP), INTENT(IN) :: gcutm
       INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
       TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
       !
       IF( dfft%nr1 == 0 .OR. dfft%nr2 == 0 .OR. dfft%nr3 == 0 ) THEN
         !
         ! ... calculate the size of the real-space dense grid for FFT
         ! ... first, an estimate of nr1,nr2,nr3, based on the max values
         ! ... of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3
         ! ... We use G*a_i = n_i => n_i .le. |Gmax||a_i|
         !
         dfft%nr1 = int ( sqrt (gcutm) * &
               sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
         dfft%nr2 = int ( sqrt (gcutm) * &
               sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
         dfft%nr3 = int ( sqrt (gcutm) * &
               sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
         !
         CALL grid_set( dfft, bg, gcutm, dfft%nr1, dfft%nr2, dfft%nr3 )
         !
       ELSE
          WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
       END IF

       IF (PRESENT(fft_fact)) THEN
          dfft%nr1 = good_fft_order( dfft%nr1, fft_fact(1) )
          dfft%nr2 = good_fft_order( dfft%nr2, fft_fact(2) )
          dfft%nr3 = good_fft_order( dfft%nr3, fft_fact(3) )
       ELSE
          dfft%nr1 = good_fft_order( dfft%nr1 )
          dfft%nr2 = good_fft_order( dfft%nr2 )
          dfft%nr3 = good_fft_order( dfft%nr3 )
       END IF

       dfft%nr1x  = good_fft_dimension( dfft%nr1 )
       dfft%nr2x  = dfft%nr2
       dfft%nr3x  = good_fft_dimension( dfft%nr3 )

     END SUBROUTINE realspace_grid_init

!=----------------------------------------------------------------------------=!

   SUBROUTINE grid_set( dfft, bg, gcut, nr1, nr2, nr3 )

!  this routine returns in nr1, nr2, nr3 the minimal 3D real-space FFT 
!  grid required to fit the G-vector sphere with G^2 <= gcut
!  On input, nr1,nr2,nr3 must be set to values that match or exceed
!  the largest i,j,k (Miller) indices in G(i,j,k) = i*b1 + j*b2 + k*b3
!  ----------------------------------------------

      IMPLICIT NONE

! ... declare arguments
      TYPE(fft_type_descriptor), INTENT(IN) :: dfft
      INTEGER, INTENT(INOUT) :: nr1, nr2, nr3
      REAL(DP), INTENT(IN) :: bg(3,3), gcut

! ... declare other variables
      INTEGER :: i, j, k, nr, nb(3)
      REAL(DP) :: gsq, g(3)

!  ----------------------------------------------

      nb     = 0

! ... calculate moduli of G vectors and the range of indices where
! ... |G|^2 < gcut (in parallel whenever possible)

      DO k = -nr3, nr3
        !
        ! ... me_image = processor number, starting from 0
        !
        IF( MOD( k + nr3, dfft%nproc ) == dfft%mype ) THEN
          DO j = -nr2, nr2
            DO i = -nr1, nr1

              g( 1 ) = DBLE(i)*bg(1,1) + DBLE(j)*bg(1,2) + DBLE(k)*bg(1,3)
              g( 2 ) = DBLE(i)*bg(2,1) + DBLE(j)*bg(2,2) + DBLE(k)*bg(2,3)
              g( 3 ) = DBLE(i)*bg(3,1) + DBLE(j)*bg(3,2) + DBLE(k)*bg(3,3)

! ...         calculate modulus

              gsq =  g( 1 )**2 + g( 2 )**2 + g( 3 )**2 

              IF( gsq < gcut ) THEN

! ...           calculate maximum index
                nb(1) = MAX( nb(1), ABS( i ) )
                nb(2) = MAX( nb(2), ABS( j ) )
                nb(3) = MAX( nb(3), ABS( k ) )
              END IF

            END DO
          END DO
        END IF
      END DO

#if defined(__MPI)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, nb, 3, MPI_INTEGER, MPI_MAX, dfft%comm, i )
#endif

! ... the size of the 3d FFT matrix depends upon the maximum indices. With
! ... the following choice, the sphere in G-space "touches" its periodic image

      nr1 = 2 * nb(1) + 1
      nr2 = 2 * nb(2) + 1
      nr3 = 2 * nb(3) + 1

      RETURN
   
   END SUBROUTINE grid_set

!=----------------------------------------------------------------------------=!
END MODULE fft_types
!=----------------------------------------------------------------------------=!
