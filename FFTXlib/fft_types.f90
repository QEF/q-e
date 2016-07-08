!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE fft_types

  IMPLICIT NONE

  SAVE

  INTEGER :: stdout = 6

  TYPE fft_dlay_descriptor

    INTEGER :: nst      ! total number of sticks
    INTEGER, POINTER :: nsp(:)   ! number of sticks per processor ( potential )
                                 ! using proc index starting from 1 !!
                                 ! on proc mpime -> nsp( mpime + 1 )
    INTEGER, POINTER :: nsw(:)   ! number of sticks per processor ( wave func )
                                 ! using proc index as above
    INTEGER :: nr1    = 0  !
    INTEGER :: nr2    = 0  ! effective FFT dimensions of the 3D grid (global)
    INTEGER :: nr3    = 0  ! 
    INTEGER :: nr1x   = 0  ! FFT grids leading dimensions
    INTEGER :: nr2x   = 0  ! dimensions of the arrays for the 3D grid (global)
    INTEGER :: nr3x   = 0  ! may differ from nr1 ,nr2 ,nr3 in order to boost performances
    LOGICAL :: dimensions_have_been_set = .FALSE.
    LOGICAL :: arrays_have_been_allocated = .FALSE.
    LOGICAL :: arrays_have_been_initialized = .FALSE.

    INTEGER :: npl    = 0  ! number of "Z" planes for this processor = npp( mpime + 1 )
    INTEGER :: nnp    = 0  ! number of 0 and non 0 sticks in a plane ( ~nr1*nr2/nproc )
    INTEGER :: nnr    = 0  ! local number of FFT grid elements  ( ~nr1*nr2*nr3/proc )
                           ! size of the arrays allocated for the FFT, local to each processor:
                           ! in parallel execution may differ from nr1x*nr2x*nr3x
                           ! Not to be confused either with nr1*nr2*nr3 
    INTEGER, POINTER :: ngl(:)   ! per proc. no. of non zero charge density/potential components
    INTEGER, POINTER :: nwl(:)   ! per proc. no. of non zero wave function plane components
    INTEGER, POINTER :: npp(:)   ! number of "Z" planes per processor
    INTEGER, POINTER :: ipp(:)   ! offset of the first "Z" plane on each proc ( 0 on the first proc!!!)
    INTEGER, POINTER :: iss(:)   ! index of the first rho stick on each proc
    INTEGER, POINTER :: isind(:) ! for each position in the plane indicate the stick index
    INTEGER, POINTER :: ismap(:) ! for each stick in the plane indicate the position
    INTEGER, POINTER :: iplp(:)  ! indicate which "Y" plane should be FFTed ( potential )
    INTEGER, POINTER :: iplw(:)  ! indicate which "Y" plane should be FFTed ( wave func )

    !
    !  descriptor id and pointer, for future use
    !
    INTEGER :: id
    INTEGER :: tptr
    !
    !  Sub (box) grid descriptor
    !
    INTEGER, POINTER :: irb(:,:)  ! the offset of the box corner
    INTEGER, POINTER :: imin3(:)  ! the starting local plane
    INTEGER, POINTER :: imax3(:)  ! the last local plane
    INTEGER, POINTER :: np3(:)    ! number of local plane for the box fft
    !
    !  fft parallelization
    !
    INTEGER :: mype     = 0          ! my processor id (starting from 0) in the fft group
    INTEGER :: comm     = 0          ! communicator of the fft gruop 
    INTEGER :: nproc    = 1          ! number of processor in the fft group
    INTEGER :: root     = 0          ! root processor
    !
    !  task groups
    !
    LOGICAL :: have_task_groups
    !
    INTEGER :: me_pgrp   = 0          ! task id for plane wave task group
    INTEGER :: nogrp     = 1          ! number of proc. in an orbital "task group"
    INTEGER :: npgrp     = 1          ! number of proc. in a plane-wave "task group"
    INTEGER :: ogrp_comm = 0          ! orbital group communicator
    INTEGER :: pgrp_comm = 0          ! plane-wave group communicator
    INTEGER, POINTER :: nolist(:) ! list of pes in orbital group
    INTEGER, POINTER :: nplist(:) ! list of pes in pw group
    !
    INTEGER :: tg_nnr = 0            ! maximum among nnr
    INTEGER, POINTER :: tg_nsw(:) ! number of sticks per task group ( wave func )
    INTEGER, POINTER :: tg_npp(:) ! number of "Z" planes per task group
    INTEGER, POINTER :: tg_snd(:) ! number of element to be sent in group redist
    INTEGER, POINTER :: tg_rcv(:) ! number of element to be received in group redist
    INTEGER, POINTER :: tg_psdsp(:)! send displacement for all to all (pack)
    INTEGER, POINTER :: tg_usdsp(:)! send displacement for all to all (unpack)
    INTEGER, POINTER :: tg_rdsp(:)! receive displacement for all to all
    INTEGER :: tg_nppx = 0  ! max of tg_npp
    INTEGER :: tg_ncpx = 0  ! max of tg_ncpx
    !
  END TYPE


  INTEGER, PRIVATE :: icount = 0


CONTAINS

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_dlay_set_dims( desc, nr1, nr2, nr3, nr1x, nr2x, nr3x)
  !
  ! routine that defines the dimensions of fft_dlay_descriptor
  ! must be called before fft_dlay_allocate and fft_dlay_set
  !
    TYPE (fft_dlay_descriptor) :: desc
    INTEGER, INTENT(in) :: nr1, nr2, nr3    ! size of real space grid
    INTEGER, INTENT(in) :: nr1x, nr2x, nr3x ! padded size of real space grid

    IF (desc%dimensions_have_been_set ) &
        CALL fftx_error__(' fft_dlay_set_dims ', ' fft dimensions already set ', 1 )

    !  Set fft actual and leading dimensions of fft_dlay_descriptor from input

    IF( nr1 > nr1x ) CALL fftx_error__( ' fft_dlay_set_dims ', ' wrong fft dimensions ', 1 )
    IF( nr2 > nr2x ) CALL fftx_error__( ' fft_dlay_set_dims ', ' wrong fft dimensions ', 2 )
    IF( nr3 > nr3x ) CALL fftx_error__( ' fft_dlay_set_dims ', ' wrong fft dimensions ', 3 )

    desc%nr1  = nr1
    desc%nr2  = nr2
    desc%nr3  = nr3
    desc%nr1x = nr1x
    desc%nr2x = nr2x
    desc%nr3x = nr3x

    desc%dimensions_have_been_set = .true.

  END SUBROUTINE fft_dlay_set_dims

!-------------------------------------------------------
  SUBROUTINE fft_dlay_allocate( desc, mype, root, nproc, comm, nogrp )
  !
  ! routine that allocate arrays of fft_dlay_descriptor
  ! must be called before fft_dlay_set
  !
    TYPE (fft_dlay_descriptor) :: desc
    INTEGER, INTENT(in) :: mype, root, nproc, comm ! mype starting from 0
    INTEGER, INTENT(in) :: nogrp   ! number of task groups
    INTEGER :: nx, ny

    IF (desc%arrays_have_been_allocated ) &
        CALL fftx_error__(' fft_dlay_allocate ', ' fft arrays already allocated ', 1 )

    IF (.NOT. desc%dimensions_have_been_set ) &
        CALL fftx_error__(' fft_dlay_allocate ', ' fft dimensions not yet set ', 1 )

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

    desc%id    = 0

    desc%mype  = mype
    desc%comm  = comm
    desc%nproc = nproc
    desc%root  = root
    desc%have_task_groups = ( nogrp > 1 )
    desc%me_pgrp = 0
    !
    IF( MOD( nproc, MAX( 1, nogrp ) ) /= 0 ) &
       CALL fftx_error__( " fft_dlay_allocate ", "the number of task groups should be a divisor of the number of MPI task ", 1 )
    IF( nogrp > nproc ) &
       CALL fftx_error__( " fft_dlay_allocate ", "the number of task groups should be less than the number of MPI task ", 1 )

    desc%nogrp = MAX( 1, nogrp )
    desc%npgrp = nproc / MAX( 1, nogrp )
    desc%ogrp_comm = 0
    desc%pgrp_comm = 0
    ALLOCATE( desc%nolist( desc%nogrp ) )
    ALLOCATE( desc%nplist( desc%npgrp ) )
    desc%nolist = 0
    desc%nplist = 0
    NULLIFY( desc%tg_nsw )
    NULLIFY( desc%tg_npp )
    NULLIFY( desc%tg_snd )
    NULLIFY( desc%tg_rcv )
    NULLIFY( desc%tg_psdsp )
    NULLIFY( desc%tg_usdsp )
    NULLIFY( desc%tg_rdsp )

    desc%arrays_have_been_allocated = .TRUE.

  END SUBROUTINE fft_dlay_allocate

  SUBROUTINE fft_dlay_deallocate( desc )
    TYPE (fft_dlay_descriptor) :: desc
    IF ( associated( desc%nsp ) )    DEALLOCATE( desc%nsp )
    IF ( associated( desc%nsw ) )    DEALLOCATE( desc%nsw )
    IF ( associated( desc%ngl ) )    DEALLOCATE( desc%ngl )
    IF ( associated( desc%nwl ) )    DEALLOCATE( desc%nwl )
    IF ( associated( desc%npp ) )    DEALLOCATE( desc%npp )
    IF ( associated( desc%ipp ) )    DEALLOCATE( desc%ipp )
    IF ( associated( desc%iss ) )    DEALLOCATE( desc%iss )
    IF ( associated( desc%isind ) )  DEALLOCATE( desc%isind )
    IF ( associated( desc%ismap ) )  DEALLOCATE( desc%ismap )
    IF ( associated( desc%iplp ) )   DEALLOCATE( desc%iplp )
    IF ( associated( desc%iplw ) )   DEALLOCATE( desc%iplw )
    IF ( associated( desc%nolist ) ) DEALLOCATE( desc%nolist )
    IF ( associated( desc%nplist ) ) DEALLOCATE( desc%nplist )
    desc%id = 0
    IF( desc%have_task_groups ) THEN
       IF ( associated( desc%tg_nsw ) )   DEALLOCATE( desc%tg_nsw )
       IF ( associated( desc%tg_npp ) )   DEALLOCATE( desc%tg_npp )
       IF ( associated( desc%tg_snd ) )   DEALLOCATE( desc%tg_snd )
       IF ( associated( desc%tg_rcv ) )   DEALLOCATE( desc%tg_rcv )
       IF ( associated( desc%tg_psdsp ) )   DEALLOCATE( desc%tg_psdsp )
       IF ( associated( desc%tg_usdsp ) )   DEALLOCATE( desc%tg_usdsp )
       IF ( associated( desc%tg_rdsp ) )   DEALLOCATE( desc%tg_rdsp )
    ENDIF
    desc%have_task_groups = .FALSE.

    desc%arrays_have_been_allocated = .FALSE.

    desc%dimensions_have_been_set = .FALSE.

  END SUBROUTINE fft_dlay_deallocate

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_allocate( desc, mype, root, nproc, comm, nat )
    TYPE (fft_dlay_descriptor) :: desc
    INTEGER, INTENT(in) :: nat, nproc, mype, root, comm  ! mype starting from 0
    ALLOCATE( desc%irb( 3, nat ) )
    ALLOCATE( desc%imin3( nat ) )
    ALLOCATE( desc%imax3( nat ) )
    ALLOCATE( desc%npp( nproc ) )
    ALLOCATE( desc%ipp( nproc ) )
    ALLOCATE( desc%np3( nat ) )
    desc%irb = 0
    desc%imin3 = 0
    desc%imax3 = 0
    desc%npp = 0
    desc%ipp = 0
    desc%np3 = 0
    desc%mype = mype
    desc%nproc = nproc
    desc%comm = comm
    desc%root = root
    desc%have_task_groups = .false.
  END SUBROUTINE fft_box_allocate

  SUBROUTINE fft_box_deallocate( desc )
    TYPE (fft_dlay_descriptor) :: desc
    IF( associated( desc%irb ) ) DEALLOCATE( desc%irb )
    IF( associated( desc%imin3 ) ) DEALLOCATE( desc%imin3 )
    IF( associated( desc%imax3 ) ) DEALLOCATE( desc%imax3 )
    IF( associated( desc%npp ) ) DEALLOCATE( desc%npp )
    IF( associated( desc%ipp ) ) DEALLOCATE( desc%ipp )
    IF( associated( desc%np3 ) ) DEALLOCATE( desc%np3 )
    desc%have_task_groups = .false.
  END SUBROUTINE fft_box_deallocate

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_dlay_set( desc, tk, nst, ub, lb, idx, in1, in2, ncp, ncpw, ngp, ngpw, st, stw )

    TYPE (fft_dlay_descriptor) :: desc

    LOGICAL, INTENT(in) :: tk               ! gamma/not-gamma logical
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
    INTEGER :: np, nq, i, is, iss, i1, i2, m1, m2, n1, n2, ip
    INTEGER :: ncpx, nppx
    INTEGER :: nr1, nr2, nr3    ! size of real space grid 
    INTEGER :: nr1x, nr2x, nr3x ! padded size of real space grid

    !  Task-grouping C. Bekas
    !
    INTEGER :: sm

    IF (.NOT. desc%arrays_have_been_allocated ) &
        CALL fftx_error__(' fft_dlay_allocate ', ' fft arrays not yet allocated ', 1 )

    IF (.NOT. desc%dimensions_have_been_set ) &
        CALL fftx_error__(' fft_dlay_set ', ' fft dimensions not yet set ', 1 )

    !  Set fft actual and leading dimensions to be used internally

    nr1  = desc%nr1
    nr2  = desc%nr2
    nr3  = desc%nr3
    nr1x = desc%nr1x
    nr2x = desc%nr2x
    nr3x = desc%nr3x

    IF( ( nr1 > nr1x ) .or. ( nr2 > nr2x ) .or. ( nr3 > nr3x ) ) &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong fft dimensions ', 1 )

    IF( ( size( desc%ngl ) < desc%nproc ) .or. ( size( desc%npp ) < desc%nproc ) .or.  &
        ( size( desc%ipp ) < desc%nproc ) .or. ( size( desc%iss ) < desc%nproc ) )     &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong descriptor dimensions ', 2 )

    IF( ( size( idx ) < nst ) .or. ( size( in1 ) < nst ) .or. ( size( in2 ) < nst ) ) &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong number of stick dimensions ', 3 )

    IF( ( size( ncp ) < desc%nproc ) .or. ( size( ngp ) < desc%nproc ) ) &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong stick dimensions ', 4 )

    !  Set the number of "xy" planes for each processor
    !  in other word do a slab partition along the z axis

    sm  = 0
    npp = 0
    IF ( desc%nproc == 1 ) THEN      ! sigle processor: npp(1)=nr3
      npp(1) = nr3
    ELSEIF( desc%nproc <= nr3 ) THEN ! normal case, more planes than processors.
      np = nr3 / desc%nproc
      nq = nr3 - np * desc%nproc
      DO i = 1, desc%nproc
        npp(i) = np
        IF ( i <= nq ) npp(i) = np + 1
      ENDDO
    ELSE                             ! less planes than processors. distribute them round robin with step nogrp...
      DO ip = 1, nr3  !  some compiler complains for empty DO loops
        DO i = 1, desc%nproc, desc%nogrp
             npp(i) = npp(i) + 1
             sm = sm + 1
             IF ( sm == nr3 ) exit
        ENDDO
        IF ( sm == nr3 ) exit
      ENDDO
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

    !  Set the proper number of sticks (the total number of 1d fft along z to be done)

    IF( .not. tk ) THEN
      desc%nst  = 2*nst - 1
    ELSE
      desc%nst  = nst
    ENDIF

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
      desc%tg_nnr = desc%nnr
    ELSE
      desc%nnr  = max( nr3x * ncpx, nr1x * nr2x * nppx )  ! this is required to contain the local data in R and G space
      desc%nnr  = max( desc%nnr, ncpx * nppx * desc%nproc )  ! this is required to use ALLTOALL instead of ALLTOALLV
      desc%nnr  = max( 1, desc%nnr ) ! ensure that desc%nrr > 0 ( for extreme parallelism )
      desc%tg_nnr = desc%nnr
      desc%tg_nnr = max( desc%tg_nnr, nr3x * ncpx ) ! this is required to contain the local data in G space (should be already granted!)
      desc%tg_nnr = max( desc%tg_nnr, nr1x * nr2x * nppx ) ! this is required to contain the local data in R space (should be already granted!)
      desc%tg_nnr = max( 1, desc%tg_nnr ) ! ensure that desc%nrr > 0 ( for extreme parallelism )
    ENDIF

    desc%ngl( 1:desc%nproc )  = ngp( 1:desc%nproc )  ! local number of g vectors (rho) per processor
    desc%nwl( 1:desc%nproc )  = ngpw( 1:desc%nproc ) ! local number of g vectors (wave) per processor

    IF( size( desc%isind ) < ( nr1x * nr2x ) ) &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong descriptor dimensions, isind ', 5 )

    IF( size( desc%iplp ) < ( nr1x ) .or. size( desc%iplw ) < ( nr1x ) ) &
      CALL fftx_error__( ' fft_dlay_set ', ' wrong descriptor dimensions, ipl ', 5 )

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

    DO iss = 1, nst
      is = idx( iss )
      i1 = in1( is )
      i2 = in2( is )
      IF( st( i1, i2 ) > 0 ) THEN
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
      CALL fftx_error__( ' fft_dlay_set ', ' wrong descriptor dimensions ', 6 )

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
      CALL fftx_error__( ' fft_dlay_set ', ' inconsistent number of sticks ', 7 )
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
      CALL fftx_error__( ' fft_dlay_set ', ' inconsistent number of sticks ', 8 )
    ENDIF

    desc%nsp( 1:desc%nproc ) = nsp( 1:desc%nproc ) ! -- number of rho sticks per processor

    icount    = icount + 1
    desc%id   = icount

    !  Initialize the pointer to the fft tables

    desc%tptr = icount

    RETURN
  END SUBROUTINE fft_dlay_set

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_set( desc, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, nat, &
                          irb, npp, ipp )

    IMPLICIT NONE

    TYPE (fft_dlay_descriptor) :: desc

    INTEGER, INTENT(in) :: nat
    INTEGER, INTENT(in) :: irb( :, : )
    INTEGER, INTENT(in) :: npp( : )
    INTEGER, INTENT(in) :: ipp( : )
    INTEGER, INTENT(in) :: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx

    INTEGER :: ir3, ibig3, irb3, imin3, imax3, nr3, isa

    IF( nat > size( desc%irb, 2 ) ) THEN
       WRITE( stdout, fmt="( ///,'NAT, SIZE = ',2I10)" ) nat, size( desc%irb, 2 )
       CALL fftx_error__(" fft_box_set ", " inconsistent dimensions ", 1 )
    ENDIF

    IF( desc%nproc > size( desc%npp ) ) &
       CALL fftx_error__(" fft_box_set ", " inconsistent dimensions ", 2 )

    desc%nr1 = nr1b
    desc%nr2 = nr2b
    desc%nr3 = nr3b
    desc%nr1x = nr1bx
    desc%nr2x = nr2bx
    desc%nr3x = nr3bx

    desc%irb( 1:3, 1:nat ) = irb( 1:3, 1:nat )
    desc%npp( 1:desc%nproc )    = npp( 1:desc%nproc )
    desc%ipp( 1:desc%nproc )    = ipp( 1:desc%nproc )

    nr3   = sum( npp( 1:desc%nproc ) )

    DO isa = 1, nat

       imin3 = nr3b
       imax3 = 1
       irb3  = irb( 3, isa )

       DO ir3 = 1, nr3b
          ibig3 = 1 + mod( irb3 + ir3 - 2, nr3 )
          IF( ibig3 < 1 .or. ibig3 > nr3 )   &
        &        CALL fftx_error__(' fft_box_set ',' ibig3 wrong ', ibig3 )
          ibig3 = ibig3 - ipp( desc%mype + 1 )
          IF ( ibig3 > 0 .and. ibig3 <= npp(desc%mype + 1) ) THEN
               imin3 = min( imin3, ir3 )
               imax3 = max( imax3, ir3 )
          ENDIF
       ENDDO

       desc%imin3( isa ) = imin3
       desc%imax3( isa ) = imax3
       desc%np3( isa )   = imax3 - imin3 + 1

    ENDDO

    desc%have_task_groups = .false.

  END SUBROUTINE fft_box_set


!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_dlay_scalar( desc, ub, lb, stw )

    IMPLICIT NONE

    TYPE (fft_dlay_descriptor) :: desc
    INTEGER, INTENT(in) :: lb(:), ub(:)
    INTEGER, INTENT(in) :: stw( lb(1) : ub(1), lb(2) : ub(2) )

    INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x
    INTEGER :: m1, m2, i1, i2

    IF (.NOT. desc%dimensions_have_been_set ) &
        CALL fftx_error__(' fft_dlay_scalar ', ' fft dimensions not yet set ', 1 )

    nr1  = desc%nr1
    nr2  = desc%nr2
    nr3  = desc%nr3
    nr1x = desc%nr1x
    nr2x = desc%nr2x
    nr3x = desc%nr3x

    IF( size( desc%iplw ) < nr1x .or. size( desc%isind ) < nr1x * nr2x ) &
      CALL fftx_error__(' fft_dlay_scalar ', ' wrong dimensions ', 1 )

    desc%isind = 0
    desc%iplw  = 0
    desc%iplp  = 1

    ! here we are setting parameter as if we were
    ! in a serial code, sticks are along X dimension
    ! and not along Z

    DO i1 = lb( 1 ), ub( 1 )
      DO i2 = lb( 2 ), ub( 2 )
        m1 = i1 + 1; IF ( m1 < 1 ) m1 = m1 + nr1
        m2 = i2 + 1; IF ( m2 < 1 ) m2 = m2 + nr2
        IF( stw( i1, i2 ) > 0 ) THEN
          desc%isind( m1 + ( m2 - 1 ) * nr1x ) =  1  ! st( i1, i2 )
          desc%iplw( m1 ) = 1
        ENDIF
      ENDDO
    ENDDO

    desc%nnr  = nr1x * nr2x * nr3x
    desc%npl  = nr3
    desc%nnp  = nr1x * nr2x
    desc%npp  = nr3
    desc%ipp  = 0
    desc%tg_nnr = desc%nnr
    !
    desc%have_task_groups = .false.

    RETURN
  END SUBROUTINE fft_dlay_scalar

END MODULE fft_types
