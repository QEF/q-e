!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE fft_smallbox_type

  IMPLICIT NONE

  SAVE

  INTEGER :: stdout = 6

  TYPE fft_box_descriptor
    !
    !  Sub (box) grid descriptor
    !
    INTEGER, ALLOCATABLE :: irb(:,:)  ! the offset of the box corner
    INTEGER, ALLOCATABLE :: imin3(:)  ! the starting local plane
    INTEGER, ALLOCATABLE :: imax3(:)  ! the last local plane
    INTEGER, ALLOCATABLE :: np3(:)    ! number of local plane for the box fft
    !
    INTEGER :: nr1    = 0  !
    INTEGER :: nr2    = 0  ! effective FFT dimensions of the 3D grid (global)
    INTEGER :: nr3    = 0  ! 
    INTEGER :: nr1x   = 0  ! FFT grids leading dimensions
    INTEGER :: nr2x   = 0  ! dimensions of the arrays for the 3D grid (global)
    INTEGER :: nr3x   = 0  ! may differ from nr1 ,nr2 ,nr3 in order to boost performances
    INTEGER :: nnr    = 0 
    !
    INTEGER, ALLOCATABLE :: npp(:)   ! number of "Z" planes per processor
    INTEGER, ALLOCATABLE :: ipp(:)   ! offset of the first "Z" plane on each proc ( 0 on the first proc!!!)
    !
    !
    !  fft parallelization
    !
    INTEGER :: mype     = 0          ! my processor id (starting from 0) in the fft group
    INTEGER :: comm     = 0          ! communicator of the fft gruop 
    INTEGER :: nproc    = 1          ! number of processor in the fft group
    INTEGER :: root     = 0          ! root processor

  END TYPE



CONTAINS

!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_allocate( desc, mype, root, nproc, comm, nat )
    TYPE (fft_box_descriptor) :: desc
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
  END SUBROUTINE fft_box_allocate

  SUBROUTINE fft_box_deallocate( desc )
    TYPE (fft_box_descriptor) :: desc
    IF( ALLOCATED( desc%irb ) ) DEALLOCATE( desc%irb )
    IF( ALLOCATED( desc%imin3 ) ) DEALLOCATE( desc%imin3 )
    IF( ALLOCATED( desc%imax3 ) ) DEALLOCATE( desc%imax3 )
    IF( ALLOCATED( desc%npp ) ) DEALLOCATE( desc%npp )
    IF( ALLOCATED( desc%ipp ) ) DEALLOCATE( desc%ipp )
    IF( ALLOCATED( desc%np3 ) ) DEALLOCATE( desc%np3 )
  END SUBROUTINE fft_box_deallocate

!=----------------------------------------------------------------------------=!
!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_set( desc, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, nat, &
                          irb, npp, ipp )

    IMPLICIT NONE

    TYPE (fft_box_descriptor) :: desc

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

    desc%nnr  = desc%nr1x * desc%nr2x * desc%nr3x

  END SUBROUTINE fft_box_set

END MODULE fft_smallbox_type
