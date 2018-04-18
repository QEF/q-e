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

  ! 
  !  Data type for FFT descriptor. 
  !
  TYPE fft_type_descriptor
    !
    ! FFT dimensions
    !
    INTEGER :: nr1    = 0  !
    INTEGER :: nr2    = 0  ! effective FFT dimensions of the 3D grid (global)
    INTEGER :: nr3    = 0  ! 
    INTEGER :: nr1x   = 0  ! FFT grids leading dimensions
    INTEGER :: nr2x   = 0  ! dimensions of the arrays for the 3D grid (global)
    INTEGER :: nr3x   = 0  ! may differ from nr1 ,nr2 ,nr3 in order to boost performances
    !
    !  Parallel layout: in reciprocal space data are organized in columns (sticks) along
    !                   the third direction and distributed across nproc processors.
    !                   In real space data are distributed in blocks comprising sections
    !                   of the Y and Z directions and complete rows in the X direction in 
    !                   a matrix of  nproc2 x nproc3  processors.
    !                   nproc = nproc2 x nproc3 and additional communicators are introduced
    !                   for data redistribution across matrix columns and rows.
    !
    ! communicators and processor coordinates
    !
    LOGICAL :: lpara  = .FALSE. ! .TRUE. if parallel FFT is active
    LOGICAL :: lgamma = .FALSE. ! .TRUE. if the grid has Gamma symmetry
    INTEGER :: root   = 0 ! root processor
    INTEGER :: comm   = MPI_COMM_NULL ! communicator for the main fft group 
    INTEGER :: comm2  = MPI_COMM_NULL ! communicator for the fft group along the second direction
    INTEGER :: comm3  = MPI_COMM_NULL ! communicator for the fft group along the third direction 
    INTEGER :: nproc  = 1 ! number of processor in the main fft group
    INTEGER :: nproc2 = 1 ! number of processor in the fft group along the second direction
    INTEGER :: nproc3 = 1 ! number of processor in the fft group along the third direction
    INTEGER :: mype   = 0 ! my processor id (starting from 0) in the fft main communicator
    INTEGER :: mype2  = 0 ! my processor id (starting from 0) in the fft communicator along the second direction (nproc2)
    INTEGER :: mype3  = 0 ! my processor id (starting from 0) in the fft communicator along the third direction (nproc3)

    INTEGER, ALLOCATABLE :: iproc(:,:) , iproc2(:), iproc3(:) ! subcommunicators proc mapping (starting from 1)
    !
    ! FFT distributed data dimensions and indices
    !
    INTEGER :: my_nr3p = 0 ! size of the "Z" section for this processor = nr3p( mype3 + 1 )    ~ nr3/nproc3
    INTEGER :: my_nr2p = 0 ! size of the "Y" section for this processor = nr2p( mype2 + 1 )    ~ nr2/nproc2

    INTEGER :: my_i0r3p = 0 ! offset of the first "Z" element of this proc in the nproc3 group = i0r3p( mype3 + 1 )
    INTEGER :: my_i0r2p = 0 ! offset of the first "Y" element of this proc in the nproc2 group = i0r2p( mype2 + 1 )

    INTEGER, ALLOCATABLE :: nr3p(:)  ! size of the "Z" section of each processor in the nproc3 group along Z
    INTEGER, ALLOCATABLE :: nr2p(:)  ! size of the "Y" section of each processor in the nproc2 group along Y
    INTEGER, ALLOCATABLE :: nr1p(:)  ! number of active "X" values ( potential ) for a given proc in the nproc2 group 
    INTEGER, ALLOCATABLE :: nr1w(:)  ! number of active "X" values ( wave func ) for a given proc in the nproc2 group
    INTEGER              :: nr1w_tg  ! total number of active "X" values ( wave func ). used in task group ffts

    INTEGER, ALLOCATABLE :: i0r3p(:) ! offset of the first "Z" element of each proc in the nproc3 group (starting from 0)
    INTEGER, ALLOCATABLE :: i0r2p(:) ! offset of the first "Y" element of each proc in the nproc2 group (starting from 0)

    INTEGER, ALLOCATABLE :: ir1p(:)  ! if >0 ir1p(m1) is the incremental index of the active ( potential ) X value of this proc
    INTEGER, ALLOCATABLE :: indp(:,:)! is the inverse of ir1p
    INTEGER, ALLOCATABLE :: ir1w(:)  ! if >0 ir1w(m1) is the incremental index of the active ( wave func ) X value of this proc
    INTEGER, ALLOCATABLE :: indw(:,:)! is the inverse of ir1w
    INTEGER, ALLOCATABLE :: ir1w_tg(:)! if >0 ir1w_tg(m1) is the incremental index of the active ( wfc ) X value in task group
    INTEGER, ALLOCATABLE :: indw_tg(:)! is the inverse of ir1w_tg

#if defined(__CUDA)
    INTEGER, POINTER, CONTIGUOUS, DEVICE :: ir1p_d(:), ir1w_d(:), ir1w_tg_d(:)
#endif

    INTEGER :: nst      ! total number of sticks ( potential )

    INTEGER, ALLOCATABLE :: nsp(:)   ! number of sticks per processor ( potential ) using proc index starting from 1 
                                     !                                              ... that is on proc mype -> nsp( mype + 1 )
    INTEGER, ALLOCATABLE :: nsw(:)   ! number of sticks per processor ( wave func ) using proc index as above
    INTEGER, ALLOCATABLE :: nsw_tg(:)! number of sticks per processor ( wave func ) using proc index as above. task group version

    INTEGER, ALLOCATABLE :: ngl(:) ! per proc. no. of non zero charge density/potential components
    INTEGER, ALLOCATABLE :: nwl(:) ! per proc. no. of non zero wave function plane components

    INTEGER :: ngm  ! my no. of non zero charge density/potential components
                    !    ngm = dfftp%ngl( dfftp%mype + 1 )
                    ! with gamma sym.   
                    !    ngm = ( dfftp%ngl( dfftp%mype + 1 ) + 1 ) / 2

    INTEGER :: ngw  ! my no. of non zero wave function plane components
                    !    ngw = dffts%nwl( dffts%mype + 1 )
                    ! with gamma sym.   
                    !    ngw = ( dffts%nwl( dffts%mype + 1 ) + 1 ) / 2

    INTEGER, ALLOCATABLE :: iplp(:) ! if > 0 is the iproc2 processor owning the active "X" value ( potential )
    INTEGER, ALLOCATABLE :: iplw(:) ! if > 0 is the iproc2 processor owning the active "X" value ( wave func )

    INTEGER :: nnp    = 0  ! number of 0 and non 0 sticks in a plane ( ~nr1*nr2/nproc )
    INTEGER :: nnr    = 0  ! local number of FFT grid elements  ( ~nr1*nr2*nr3/nproc )
                           ! size of the arrays allocated for the FFT, local to each processor:
                           ! in parallel execution may differ from nr1x*nr2x*nr3x
                           ! Not to be confused either with nr1*nr2*nr3 
    INTEGER :: nnr_tg = 0  ! local number of grid elements for task group FFT ( ~nr1*nr2*nr3/proc3 )
    INTEGER, ALLOCATABLE :: iss(:)   ! index of the first rho stick on each proc
    INTEGER, ALLOCATABLE :: isind(:) ! for each position in the plane indicate the stick index
    INTEGER, ALLOCATABLE :: ismap(:) ! for each stick in the plane indicate the position
#if defined(__CUDA)
    INTEGER, POINTER, CONTIGUOUS, DEVICE :: ismap_d(:)
#endif
    INTEGER, ALLOCATABLE :: nl(:)    ! position of the G vec in the FFT grid
    INTEGER, ALLOCATABLE :: nlm(:)   ! with gamma sym. position of -G vec in the FFT grid
    !
    ! task group ALLTOALL communication layout
    INTEGER, ALLOCATABLE :: tg_snd(:) ! number of elements to be sent in task group redistribution
    INTEGER, ALLOCATABLE :: tg_rcv(:) ! number of elements to be received in task group redistribution
    INTEGER, ALLOCATABLE :: tg_sdsp(:)! send displacement for task group A2A communication
    INTEGER, ALLOCATABLE :: tg_rdsp(:)! receive displacement for task group A2A communicattion
    !
    LOGICAL :: has_task_groups = .FALSE.
    !
    CHARACTER(len=12):: rho_clock_label  = ' '
    CHARACTER(len=12):: wave_clock_label = ' '

    INTEGER :: grid_id

  END TYPE

  REAL(DP) :: fft_dual = 4.0d0
  INTEGER  :: incremental_grid_identifier = 0

  PUBLIC :: fft_type_descriptor, fft_type_init
  PUBLIC :: fft_type_allocate, fft_type_deallocate
  PUBLIC :: fft_stick_index

CONTAINS

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_allocate( desc, at, bg, gcutm, comm, fft_fact, nyfft  )
  !
  ! routine that allocate arrays of fft_type_descriptor
  ! must be called before fft_type_init
  !
    TYPE (fft_type_descriptor) :: desc
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
    REAL(DP), INTENT(IN) :: gcutm
    INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
    INTEGER, INTENT(IN), OPTIONAL :: nyfft
    INTEGER, INTENT(in) :: comm ! mype starting from 0
    INTEGER :: nx, ny, ierr, nzfft
    INTEGER :: mype, root, nproc, nproc2, nproc3, iproc, iproc2, iproc3 ! mype starting from 0
    INTEGER :: color, key, comm2, comm3
     !write (6,*) ' inside fft_type_allocate' ; FLUSH(6)

    IF ( ALLOCATED( desc%nsp ) ) &
        CALL fftx_error__(' fft_type_allocate ', ' fft arrays already allocated ', 1 )

    desc%comm = comm 
#if defined(__MPI)
    IF( desc%comm == MPI_COMM_NULL ) THEN
       CALL fftx_error__( ' fft_type_allocate ', ' fft communicator is null ', 1 )
    END IF
#endif
    !
    root = 0 ; mype = 0 ; nproc = 1
#if defined(__MPI)
    CALL MPI_COMM_RANK( comm, mype, ierr )
    CALL MPI_COMM_SIZE( comm, nproc, ierr )
#endif
    desc%root = root ; desc%mype = mype ; desc%nproc   = nproc

    IF ( present(nyfft) ) THEN
      ! check on yfft group dimension
      CALL fftx_error__( ' fft_type_allocate ', ' MOD(nproc,nyfft) .ne. 0 ', MOD(nproc,nyfft) )

!#define ZCOMPACT
#if defined(__MPI)
#if defined(ZCOMPACT)
      !write (6,*) ' FFT IS ZCOMPACT '
      nzfft = nproc / nyfft
      color = mype / nzfft       ;      key   = MOD( mype, nzfft )
#else
      !write (6,*) ' FFT IS YCOMPACT '
      color = MOD( mype, nyfft ) ;      key   = mype / nyfft
#endif

      ! processes with the same key are in the same group along Y
      CALL MPI_COMM_SPLIT( comm, key, color, desc%comm2, ierr )
      CALL MPI_COMM_RANK( desc%comm2, desc%mype2, ierr )
      CALL MPI_COMM_SIZE( desc%comm2, desc%nproc2, ierr )

      ! processes with the same color are in the same group along Z
      CALL MPI_COMM_SPLIT( comm, color, key, desc%comm3, ierr )
      CALL MPI_COMM_RANK( desc%comm3, desc%mype3, ierr )
      CALL MPI_COMM_SIZE( desc%comm3, desc%nproc3, ierr )
#else
      desc%comm2 = desc%comm ; desc%mype2 = desc%mype ; desc%nproc2 = desc%nproc
      desc%comm3 = desc%comm ; desc%mype3 = desc%mype ; desc%nproc3 = desc%nproc
#endif

    ENDIF
    !write (6,*) '  nproc and  mype  '
    !write (6,*) desc%nproc, desc%nproc2, desc%nproc3
    !write (6,*) desc%mype, desc%mype2, desc%mype3

    ALLOCATE ( desc%iproc(desc%nproc2,desc%nproc3), desc%iproc2(desc%nproc), desc%iproc3(desc%nproc) )
    do iproc = 1, desc%nproc
#if defined(ZCOMPACT)
       iproc3 = MOD(iproc-1, desc%nproc3) + 1 ; iproc2 = (iproc-1)/desc%nproc3 + 1 
#else
       iproc2 = MOD(iproc-1, desc%nproc2) + 1 ; iproc3 = (iproc-1)/desc%nproc2 + 1
#endif
       desc%iproc2(iproc) = iproc2 ; desc%iproc3(iproc) = iproc3
       desc%iproc(iproc2,iproc3) = iproc
    end do

    CALL realspace_grid_init( desc, at, bg, gcutm, fft_fact )

    ALLOCATE( desc%nr2p ( desc%nproc2 ), desc%i0r2p( desc%nproc2 ) ) ; desc%nr2p = 0 ; desc%i0r2p = 0
    ALLOCATE( desc%nr3p ( desc%nproc3 ), desc%i0r3p( desc%nproc3 ) ) ; desc%nr3p = 0 ; desc%i0r3p = 0

    nx = desc%nr1x 
    ny = desc%nr2x

    ALLOCATE( desc%nsp( desc%nproc ) ) ; desc%nsp   = 0
    ALLOCATE( desc%nsw( desc%nproc ) ) ; desc%nsw   = 0
    ALLOCATE( desc%nsw_tg( desc%nproc ) ) ; desc%nsw_tg   = 0
    ALLOCATE( desc%ngl( desc%nproc ) ) ; desc%ngl   = 0
    ALLOCATE( desc%nwl( desc%nproc ) ) ; desc%nwl   = 0
    ALLOCATE( desc%iss( desc%nproc ) ) ; desc%iss   = 0
    ALLOCATE( desc%isind( nx * ny ) ) ; desc%isind = 0
    ALLOCATE( desc%ismap( nx * ny ) ) ; desc%ismap = 0
    ALLOCATE( desc%nr1p( desc%nproc2 ) ) ; desc%nr1p  = 0
    ALLOCATE( desc%nr1w( desc%nproc2 ) ) ; desc%nr1w  = 0
    ALLOCATE( desc%ir1p( desc%nr1x ) ) ; desc%ir1p  = 0
    ALLOCATE( desc%indp( desc%nr1x,desc%nproc2 ) ) ; desc%indp  = 0
    ALLOCATE( desc%ir1w( desc%nr1x ) ) ; desc%ir1w  = 0
    ALLOCATE( desc%ir1w_tg( desc%nr1x ) ) ; desc%ir1w_tg  = 0
    ALLOCATE( desc%indw( desc%nr1x, desc%nproc2 ) ) ; desc%indw  = 0
    ALLOCATE( desc%indw_tg( desc%nr1x ) ) ; desc%indw_tg  = 0
    ALLOCATE( desc%iplp( nx ) ) ; desc%iplp  = 0
    ALLOCATE( desc%iplw( nx ) ) ; desc%iplw  = 0

    ALLOCATE( desc%tg_snd( desc%nproc2) ) ; desc%tg_snd = 0
    ALLOCATE( desc%tg_rcv( desc%nproc2) ) ; desc%tg_rcv = 0
    ALLOCATE( desc%tg_sdsp( desc%nproc2) ) ; desc%tg_sdsp = 0
    ALLOCATE( desc%tg_rdsp( desc%nproc2) ) ; desc%tg_rdsp = 0

#if defined(__CUDA)
    ALLOCATE( desc%ir1p_d( desc%nr1x ) ) ; desc%ir1p_d  = 0
    ALLOCATE( desc%ir1w_d( desc%nr1x ) ) ; desc%ir1w_d  = 0
    ALLOCATE( desc%ir1w_tg_d( desc%nr1x ) ) ; desc%ir1w_tg_d  = 0
    ALLOCATE( desc%ismap_d( nx * ny ) ) ; desc%ismap_d = 0
#endif

    incremental_grid_identifier = incremental_grid_identifier + 1
    desc%grid_id = incremental_grid_identifier

  END SUBROUTINE fft_type_allocate

  SUBROUTINE fft_type_deallocate( desc )
    TYPE (fft_type_descriptor) :: desc
    INTEGER :: ierr
     !write (6,*) ' inside fft_type_deallocate' ; FLUSH(6)
    IF ( ALLOCATED( desc%nr2p ) )   DEALLOCATE( desc%nr2p )
    IF ( ALLOCATED( desc%i0r2p ) )  DEALLOCATE( desc%i0r2p )
    IF ( ALLOCATED( desc%nr3p ) )   DEALLOCATE( desc%nr3p )
    IF ( ALLOCATED( desc%i0r3p ) )  DEALLOCATE( desc%i0r3p )
    IF ( ALLOCATED( desc%nsp ) )    DEALLOCATE( desc%nsp )
    IF ( ALLOCATED( desc%nsw ) )    DEALLOCATE( desc%nsw )
    IF ( ALLOCATED( desc%nsw_tg ) ) DEALLOCATE( desc%nsw_tg )
    IF ( ALLOCATED( desc%ngl ) )    DEALLOCATE( desc%ngl )
    IF ( ALLOCATED( desc%nwl ) )    DEALLOCATE( desc%nwl )
    IF ( ALLOCATED( desc%iss ) )    DEALLOCATE( desc%iss )
    IF ( ALLOCATED( desc%isind ) )  DEALLOCATE( desc%isind )
    IF ( ALLOCATED( desc%ismap ) )  DEALLOCATE( desc%ismap )
    IF ( ALLOCATED( desc%nr1p ) )   DEALLOCATE( desc%nr1p )
    IF ( ALLOCATED( desc%nr1w ) )   DEALLOCATE( desc%nr1w )
    IF ( ALLOCATED( desc%ir1p ) )   DEALLOCATE( desc%ir1p )
    IF ( ALLOCATED( desc%indp ) )   DEALLOCATE( desc%indp )
    IF ( ALLOCATED( desc%ir1w ) )   DEALLOCATE( desc%ir1w )
    IF ( ALLOCATED( desc%ir1w_tg ) )DEALLOCATE( desc%ir1w_tg )
    IF ( ALLOCATED( desc%indw ) )   DEALLOCATE( desc%indw )
    IF ( ALLOCATED( desc%indw_tg ) )DEALLOCATE( desc%indw_tg )
    IF ( ALLOCATED( desc%iplp ) )   DEALLOCATE( desc%iplp )
    IF ( ALLOCATED( desc%iplw ) )   DEALLOCATE( desc%iplw )
    IF ( ALLOCATED( desc%iproc ) )  DEALLOCATE( desc%iproc )
    IF ( ALLOCATED( desc%iproc2 ) ) DEALLOCATE( desc%iproc2 )
    IF ( ALLOCATED( desc%iproc3 ) ) DEALLOCATE( desc%iproc3 )

    IF ( ALLOCATED( desc%tg_snd ) ) DEALLOCATE( desc%tg_snd )
    IF ( ALLOCATED( desc%tg_rcv ) ) DEALLOCATE( desc%tg_rcv )
    IF ( ALLOCATED( desc%tg_sdsp ) )DEALLOCATE( desc%tg_sdsp )
    IF ( ALLOCATED( desc%tg_rdsp ) )DEALLOCATE( desc%tg_rdsp )

    IF ( ALLOCATED( desc%nl ) )  DEALLOCATE( desc%nl )
    IF ( ALLOCATED( desc%nlm ) ) DEALLOCATE( desc%nlm )

#if defined(__CUDA)
    IF ( ALLOCATED( desc%ismap_d ) )   DEALLOCATE( desc%ismap_d )
    IF ( ALLOCATED( desc%ir1p_d ) )    DEALLOCATE( desc%ir1p_d )
    IF ( ALLOCATED( desc%ir1w_d ) )    DEALLOCATE( desc%ir1w_d )
    IF ( ALLOCATED( desc%ir1w_tg_d ) ) DEALLOCATE( desc%ir1w_tg_d )
#endif

    desc%comm  = MPI_COMM_NULL 
#if defined(__MPI)
    IF (desc%comm2 /= MPI_COMM_NULL) CALL MPI_COMM_FREE( desc%comm2, ierr )
    IF (desc%comm3 /= MPI_COMM_NULL) CALL MPI_COMM_FREE( desc%comm3, ierr )
#else
    desc%comm2 = MPI_COMM_NULL 
    desc%comm3 = MPI_COMM_NULL 
#endif

    desc%nr1    = 0 ; desc%nr2    = 0 ; desc%nr3    = 0  
    desc%nr1x   = 0 ; desc%nr2x   = 0 ; desc%nr3x   = 0  

    desc%grid_id = 0

  END SUBROUTINE fft_type_deallocate

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_set( desc, nst, ub, lb, idx, in1, in2, ncp, ncpw, ngp, ngpw, st, stw )

    TYPE (fft_type_descriptor) :: desc

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

    INTEGER :: nsp( desc%nproc ), nsw_tg, nr1w_tg
    INTEGER :: np, nq, i, is, iss, itg, i1, i2, m1, m2, ip
    INTEGER :: ncpx, nr1px, nr2px, nr3px
    INTEGER :: nr1, nr2, nr3    ! size of real space grid 
    INTEGER :: nr1x, nr2x, nr3x ! padded size of real space grid
     !write (6,*) ' inside fft_type_set' ; FLUSH(6)
    !
    IF (.NOT. ALLOCATED( desc%nsp ) ) &
        CALL fftx_error__(' fft_type_set ', ' fft arrays not yet allocated ', 1 )

    IF ( desc%nr1 == 0 .OR. desc%nr2 == 0 .OR. desc%nr3 == 0 ) &
        CALL fftx_error__(' fft_type_set ', ' fft dimensions not yet set ', 1 )

    !  Set fft actual and leading dimensions to be used internally

    nr1  = desc%nr1  ; nr2  = desc%nr2  ; nr3  = desc%nr3
    nr1x = desc%nr1x ; nr2x = desc%nr2x ; nr3x = desc%nr3x

    IF( ( nr1 > nr1x ) .or. ( nr2 > nr2x ) .or. ( nr3 > nr3x ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong fft dimensions ', 1 )

    IF( ( size( desc%ngl ) < desc%nproc ) .or.  ( size( desc%iss ) < desc%nproc ) .or. &
        ( size( desc%nr2p ) < desc%nproc2 ) .or. ( size( desc%i0r2p ) < desc%nproc2 ) .or. &
        ( size( desc%nr3p ) < desc%nproc3 ) .or. ( size( desc%i0r3p ) < desc%nproc3 ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong descriptor dimensions ', 2 )

    IF( ( size( idx ) < nst ) .or. ( size( in1 ) < nst ) .or. ( size( in2 ) < nst ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong number of stick dimensions ', 3 )

    IF( ( size( ncp ) < desc%nproc ) .or. ( size( ngp ) < desc%nproc ) ) &
      CALL fftx_error__( ' fft_type_set ', ' wrong stick dimensions ', 4 )

    !  Set the number of "Y" values for each processor in the nproc2 group
    np = nr2 / desc%nproc2
    nq = nr2 - np * desc%nproc2
    desc%nr2p(1:desc%nproc2) = np    ! assign a base value to all processors of the nproc2 group
    DO i =1, nq ! assign an extra unit to the first nq processors of the nproc2 group
       desc%nr2p(i) = np + 1
    ENDDO
    !-- my_nr2p is the number of planes per processor of this processor   in the Y group
    desc%my_nr2p = desc%nr2p( desc%mype2 + 1 )    

    !  Find out the index of the starting plane on each proc
    desc%i0r2p = 0
    DO i = 2, desc%nproc2
       desc%i0r2p(i) = desc%i0r2p(i-1) + desc%nr2p(i-1)
    ENDDO
    !-- my_i0r2p is the index-offset of the starting plane of this processor  in the Y group
    desc%my_i0r2p = desc%i0r2p( desc%mype2 + 1 ) 

    !  Set the number of "Z" values for each processor in the nproc3 group
    np = nr3 / desc%nproc3
    nq = nr3 - np * desc%nproc3
    desc%nr3p(1:desc%nproc3) = np    ! assign a base value to all processors
    DO i =1, nq ! assign an extra unit to the first nq processors of the nproc3 group
       desc%nr3p(i) = np + 1
    END DO
    !-- my_nr3p is the number of planes per processor of this processor   in the Z group
    desc%my_nr3p = desc%nr3p( desc%mype3 + 1 ) 

    !  Find out the index of the starting plane on each proc
    desc%i0r3p  = 0
    DO i = 2, desc%nproc3
       desc%i0r3p( i )  = desc%i0r3p( i-1 ) + desc%nr3p ( i-1 )  
    ENDDO
    !-- my_i0r3p is the index-offset of the starting plane of this processor  in the Z group
    desc%my_i0r3p = desc%i0r3p( desc%mype3 + 1 ) 

    ! dimension of the xy plane. see ncplane

    desc%nnp  = nr1x * nr2x  

!!!!!!!

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

    desc%isind = 0  ! will contain the +ve or -ve of the processor number, if any, that owns the stick
    desc%iplp  = 0  ! if > 0 is the nporc2 processor owning this ( potential ) X active plane
    desc%iplw  = 0  ! if > 0 is the nproc2 processor owning this ( wave func ) X active value 

    !  Set nst to the proper number of sticks (the total number of 1d fft along z to be done)
    
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
          desc%iplw( m1 ) = desc%iproc2(st(i1,i2))
        ELSE
          desc%isind( m1 + ( m2 - 1 ) * nr1x ) = -st( i1, i2 )
        ENDIF
        desc%iplp( m1 ) = desc%iproc2(st(i1,i2))
        IF( desc%lgamma ) THEN
          IF( i1 /= 0 .OR. i2 /= 0 ) desc%nst = desc%nst + 1
          m1 = -i1 + 1; IF ( m1 < 1 ) m1 = m1 + nr1
          m2 = -i2 + 1; IF ( m2 < 1 ) m2 = m2 + nr2
          IF( stw( -i1, -i2 ) > 0 ) THEN
            desc%isind( m1 + ( m2 - 1 ) * nr1x ) =  st( -i1, -i2 )
            desc%iplw( m1 ) = desc%iproc2(st(-i1,-i2))
          ELSE
            desc%isind( m1 + ( m2 - 1 ) * nr1x ) = -st( -i1, -i2 )
          ENDIF
          desc%iplp( m1 ) = desc%iproc2(st(-i1,-i2))
        ENDIF
      ENDIF
    ENDDO
    do m1=1,desc%nr1x
       if (desc%iplw(m1)>0) then
          if (desc%iplp(m1) /= desc%iplw(m1) )  then
             write (6,*) 'WRONG iplp/iplw arrays'
             write (6,*) desc%iplp
             write (6,*) desc%iplw
             CALL fftx_error__( ' fft_type_set ', ' iplp is wrong ', m1 )
          end if
       end if
    end do
    ! count how many active X values per each nproc2 processor and set the incremental index of this one

    ! wave func X values first
    desc%nr1w = 0 ; desc%ir1w = 0 ; desc%indw = 0 
    nr1w_tg = 0 ; desc%ir1w_tg = 0 ; desc%indw_tg = 0 
    do i1 = 1, nr1
       if (desc%iplw(i1) > 0 ) then
          desc%nr1w(desc%iplw(i1)) =  desc%nr1w(desc%iplw(i1)) + 1 
          desc%indw(desc%nr1w(desc%iplw(i1)),desc%iplw(i1)) = i1
          nr1w_tg = nr1w_tg + 1 ; desc%ir1w_tg(i1) = nr1w_tg ; desc%indw_tg(nr1w_tg) = i1
       end if
       if (desc%iplw(i1) == desc%mype2 +1) desc%ir1w(i1) = desc%nr1w(desc%iplw(i1)) 
    end do
    desc%nr1w_tg = nr1w_tg ! this is useful in task group ffts

    ! then potential X values 
    desc%nr1p = desc%nr1w ; desc%ir1p=desc%ir1w ; desc%indp = desc%indw
    do i1 = 1, nr1
       if ( (desc%iplw(i1) > 0) .and. (desc%iplp(i1) == 0) ) &
             CALL fftx_error__( ' fft_type_set ', ' bad distribution of X values ', i1 )
       if ( (desc%iplw(i1) > 0) ) cycle ! this X value has already been taken care of

       if (desc%iplp(i1) > 0 ) then
          desc%nr1p(desc%iplp(i1)) =  desc%nr1p(desc%iplp(i1)) + 1 
          desc%indp(desc%nr1p(desc%iplp(i1)),desc%iplp(i1)) = i1
       end if
       if (desc%iplp(i1) == desc%mype2+1) desc%ir1p(i1) = desc%nr1p(desc%iplp(i1)) 
    end do

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

    !  check number of sticks against the input value

    IF( any( nsp( 1:desc%nproc ) /= ncpw( 1:desc%nproc ) ) ) THEN
      DO ip = 1, desc%nproc
        WRITE( stdout,*)  ' * ', ip, ' * ', nsp( ip ), ' /= ', ncpw( ip )
      ENDDO
      CALL fftx_error__( ' fft_type_set ', ' inconsistent number of sticks ', 7 )
    ENDIF

    desc%nsw( 1:desc%nproc ) = nsp( 1:desc%nproc )  ! -- number of wave sticks per processor

    ! -- number of wave sticks per processor for task group ffts
    desc%nsw_tg( 1:desc%nproc ) = 0
    do ip =1, desc%nproc3
       nsw_tg = sum(desc%nsw(desc%iproc(1:desc%nproc2,ip)))
       desc%nsw_tg(desc%iproc(1:desc%nproc2,ip)) = nsw_tg 
    end do

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

    !  check number of sticks against the input value

    IF( any( nsp( 1:desc%nproc ) /= ncp( 1:desc%nproc ) ) ) THEN
      DO ip = 1, desc%nproc
        WRITE( stdout,*)  ' * ', ip, ' * ', nsp( ip ), ' /= ', ncp( ip )
      ENDDO
      CALL fftx_error__( ' fft_type_set ', ' inconsistent number of sticks ', 8 )
    ENDIF

    desc%nsp( 1:desc%nproc ) = nsp( 1:desc%nproc ) ! -- number of rho sticks per processor

    IF( .NOT. desc%lpara ) THEN

       desc%isind = 0
       desc%iplw  = 0
       desc%iplp  = 1

       ! here we are setting parameter as if we were in a serial code, 
       ! sticks are along X dimension and not along Z
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
       desc%nnr  = nr1x * nr2x * nr3x
       desc%nnp  = nr1x * nr2x
       desc%my_nr2p = nr2 ;  desc%nr2p = nr2 ;  desc%i0r2p = 0
       desc%my_nr3p = nr3 ;  desc%nr3p = nr3 ;  desc%i0r3p = 0
       desc%nsw = desc%nsw(1)
       desc%nsp = desc%nsp(1)
       desc%ngl  = SUM(ngp)  
       desc%nwl  = SUM(ngpw) 
       !
    END IF

    !write (6,*) 'fft_type_set SUMMARY'
    !write (6,*) 'desc%mype ', desc%mype
    !write (6,*) 'desc%mype2', desc%mype2
    !write (6,*) 'desc%mype3', desc%mype3
    !write (6,*) 'nr1  nr2  nr3  dimensions : ', desc%nr1, desc%nr2, desc%nr3
    !write (6,*) 'nr1x nr2x nr3x dimensions : ', desc%nr1x, desc%nr2x, desc%nr3x
    !write (6,*) 'nr3p arrays'
    !write (6,*) desc%nr3p
    !write (6,*) 'i0r3p arrays'
    !write (6,*) desc%i0r3p
    !write (6,*) 'nr2p arrays'
    !write (6,*) desc%nr2p
    !write (6,*) 'i0r2p arrays'
    !write (6,*) desc%i0r2p
    !write (6,*) 'nsp/nsw arrays'
    !write (6,*) desc%nsp
    !write (6,*) desc%nsw
    !write (6,*) 'nr1p/nr1w arrays'
    !write (6,*) desc%nr1p
    !write (6,*) desc%nr1w
    !write (6,*) 'ir1p/ir1w arrays'
    !write (6,*) desc%ir1p
    !write (6,*) desc%ir1w
    !write (6,*) 'indp/indw arrays'
    !write (6,*) desc%indp
    !write (6,*) desc%indw
    !write (6,*) 'iplp/iplw arrays'
    !write (6,*) desc%iplp
    !write (6,*) desc%iplw

    !  Finally set fft local workspace dimension

    nr1px = MAXVAL( desc%nr1p( 1:desc%nproc2 ) )  ! maximum number of X values per processor in the nproc2 group
    nr2px = MAXVAL( desc%nr2p( 1:desc%nproc2 ) )  ! maximum number of planes per processor in the nproc2 group
    nr3px = MAXVAL( desc%nr3p( 1:desc%nproc3 ) )  ! maximum number of planes per processor in the nproc3 group
    ncpx  = MAXVAL( ncp( 1:desc%nproc ) ) ! maximum number of columns per processor (use potential sticks to be safe)

    IF ( desc%nproc == 1 ) THEN
      desc%nnr  = nr1x * nr2x * nr3x
      desc%nnr_tg = desc%nnr * desc%nproc2 
    ELSE
      desc%nnr  = max( ncpx * nr3x, nr1x * nr2px * nr3px )  ! this is required to contain the local data in R and G space
      desc%nnr  = max( desc%nnr, ncpx*nr3px*desc%nproc3, nr1px*nr2px*nr3px*desc%nproc2)  ! this is required to use ALLTOALL instead of ALLTOALLV
      desc%nnr  = max( 1, desc%nnr ) ! ensure that desc%nrr > 0 ( for extreme parallelism )
      desc%nnr_tg = desc%nnr * desc%nproc2 
    ENDIF

    !write (6,*) ' nnr bounds'
    !write (6,*) ' nr1x ',nr1x,' nr2x ', nr2x, ' nr3x ', nr3x
    !write (6,*) ' nr1x * nr2x * nr3x',nr1x * nr2x * nr3x
    !write (6,*) ' ncpx ',ncpx,' nr3px ', nr3px, ' desc%nproc3 ', desc%nproc3
    !write (6,*) ' ncpx * nr3x ',ncpx * nr3x 
    !write (6,*) ' ncpx * nr3px * desc%nproc3 ',ncpx*nr3px*desc%nproc3
    !write (6,*) ' nr1px ', nr1px,' nr2px ',nr2px,' desc%nproc2 ', desc%nproc2
    !write (6,*) ' nr1x * nr2px * nr3px ',nr1x * nr2px * nr3px 
    !write (6,*) ' nr1px * nr2px *nr3px * desc%nproc2 ',nr1px*nr2px*nr3px*desc%nproc2
    !write (6,*) ' desc%nnr ', desc%nnr

    IF( desc%nr3x * desc%nsw( desc%mype + 1 ) > desc%nnr ) &
        CALL fftx_error__( ' task_groups_init ', ' inconsistent desc%nnr ', 1 )
    desc%tg_snd(1)  = desc%nr3x * desc%nsw( desc%mype + 1 )
    desc%tg_rcv(1)  = desc%nr3x * desc%nsw( desc%iproc(1,desc%mype3+1) )
    desc%tg_sdsp(1) = 0
    desc%tg_rdsp(1) = 0
    DO i = 2, desc%nproc2
       desc%tg_snd(i)  = desc%nr3x * desc%nsw( desc%mype + 1 )
       desc%tg_rcv(i)  = desc%nr3x * desc%nsw( desc%iproc(i,desc%mype3+1) )
       desc%tg_sdsp(i) = desc%tg_sdsp(i-1) + desc%nnr
       desc%tg_rdsp(i) = desc%tg_rdsp(i-1) + desc%tg_rcv(i-1)
    ENDDO

#if defined(__CUDA)
    desc%ismap_d = desc%ismap
    desc%ir1p_d = desc%ir1p
    desc%ir1w_d  = desc%ir1w
    desc%ir1w_tg_d = desc%ir1w_tg
#endif

    RETURN

  END SUBROUTINE fft_type_set

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_type_init( dfft, smap, pers, lgamma, lpara, comm, at, bg, gcut_in, dual_in, nyfft )

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
     INTEGER, INTENT(IN) :: nyfft
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

     REAL(DP) :: gcut, gkcut, dual
     INTEGER  :: ngm, ngw
     !write (6,*) ' inside fft_type_init' ; FLUSH(6)

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
     !write (*,*) 'FFT_TYPE_INIT pers, gkcut,gcut', pers, gkcut, gcut

     IF( .NOT. ALLOCATED( dfft%nsp ) ) THEN
        CALL fft_type_allocate( dfft, at, bg, gcut, comm, nyfft=nyfft )
     ELSE
        IF( dfft%comm /= comm ) THEN
           CALL fftx_error__(' fft_type_init ', ' FFT already allocated with a different communicator ', 1 )
        END IF
     END IF


     dfft%lpara = lpara  !  this descriptor can be either a descriptor for a
                         !  parallel FFT or a serial FFT even in parallel build

     CALL sticks_map_allocate( smap, lgamma, dfft%lpara, dfft%nproc2, dfft%iproc, dfft%iproc2, &
                                             dfft%nr1, dfft%nr2, dfft%nr3, bg, dfft%comm )

     dfft%lgamma = smap%lgamma ! .TRUE. if the grid has Gamma symmetry

     ALLOCATE( stw ( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) )
     ALLOCATE( st  ( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) )
     ALLOCATE( nstp(smap%nproc) )
     ALLOCATE( sstp(smap%nproc) )
     ALLOCATE( nstpw(smap%nproc) )
     ALLOCATE( sstpw(smap%nproc) )

     !write(*,*) 'calling get_sticks with gkcut =',gkcut
     CALL get_sticks(  smap, gkcut, nstpw, sstpw, stw, nstw, ngw )
     !write(*,*) 'calling get_sticks with gcut =',gcut
     CALL get_sticks(  smap, gcut,  nstp, sstp, st, nst, ngm )

     CALL fft_type_set( dfft, nst, smap%ub, smap%lb, smap%idx, &
                             smap%ist(:,1), smap%ist(:,2), nstp, nstpw, sstp, sstpw, st, stw )

     dfft%ngw = dfft%nwl( dfft%mype + 1 )
     dfft%ngm = dfft%ngl( dfft%mype + 1 )
     IF( dfft%lgamma ) THEN
        dfft%ngw = (dfft%ngw + 1)/2
        dfft%ngm = (dfft%ngm + 1)/2
     END IF

     IF( dfft%ngw /= ngw ) THEN
        CALL fftx_error__(' fft_type_init ', ' wrong ngw ', 1 )
     END IF
     IF( dfft%ngm /= ngm ) THEN
        CALL fftx_error__(' fft_type_init ', ' wrong ngm ', 1 )
     END IF


     DEALLOCATE( st )
     DEALLOCATE( stw )
     DEALLOCATE( nstp )
     DEALLOCATE( sstp )
     DEALLOCATE( nstpw )
     DEALLOCATE( sstpw )

  END SUBROUTINE fft_type_init

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
     !write (6,*) ' inside realspace_grid_init' ; FLUSH(6)
       !
       IF( dfft%nr1 == 0 .OR. dfft%nr2 == 0 .OR. dfft%nr3 == 0 ) THEN
         !
         ! ... calculate the size of the real-space dense grid for FFT
         ! ... first, an estimate of nr1,nr2,nr3, based on the max values
         ! ... of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3
         ! ... We use G*a_i = n_i => n_i .le. |Gmax||a_i|
         !
         dfft%nr1 = int ( sqrt (gcutm) * sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
         dfft%nr2 = int ( sqrt (gcutm) * sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
         dfft%nr3 = int ( sqrt (gcutm) * sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1

         !write (6,*) sqrt(gcutm)*sqrt(at(1,1)**2 + at(2,1)**2 + at(3,1)**2) , dfft%nr1
         !write (6,*) sqrt(gcutm)*sqrt(at(1,2)**2 + at(2,2)**2 + at(3,2)**2) , dfft%nr2
         !write (6,*) sqrt(gcutm)*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2) , dfft%nr3
         !
         CALL grid_set( dfft, bg, gcutm, dfft%nr1, dfft%nr2, dfft%nr3 )
         !
#if defined (__DEBUG)
       ELSE
          WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
#endif
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
     !write (6,*) ' inside grid_set' ; FLUSH(6)

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


   PURE FUNCTION fft_stick_index( desc, i, j )
      IMPLICIT NONE
      TYPE(fft_type_descriptor), INTENT(IN) :: desc
      INTEGER :: fft_stick_index
      INTEGER, INTENT(IN) :: i, j
      INTEGER :: mc, m1, m2
      m1 = mod (i, desc%nr1) + 1
      IF (m1 < 1) m1 = m1 + desc%nr1
      m2 = mod (j, desc%nr2) + 1
      IF (m2 < 1) m2 = m2 + desc%nr2
      mc = m1 + (m2 - 1) * desc%nr1x
      fft_stick_index = desc%isind ( mc ) 
   END FUNCTION

!=----------------------------------------------------------------------------=!
END MODULE fft_types
!=----------------------------------------------------------------------------=!
