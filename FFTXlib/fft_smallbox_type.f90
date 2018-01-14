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
    INTEGER, ALLOCATABLE :: imin2(:),imin3(:)  ! the starting local plane
    INTEGER, ALLOCATABLE :: imax2(:),imax3(:)  ! the last local plane
    INTEGER, ALLOCATABLE :: np2(:),np3(:)    ! number of local plane for the box fft
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
    ALLOCATE( desc%imin2( nat ) )
    ALLOCATE( desc%imax2( nat ) )
    ALLOCATE( desc%imin3( nat ) )
    ALLOCATE( desc%imax3( nat ) )
    ALLOCATE( desc%npp( nproc ) )
    ALLOCATE( desc%ipp( nproc ) )
    ALLOCATE( desc%np2( nat ) )
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
    IF( ALLOCATED( desc%imin2 ) ) DEALLOCATE( desc%imin2 )
    IF( ALLOCATED( desc%imax2 ) ) DEALLOCATE( desc%imax2 )
    IF( ALLOCATED( desc%imin3 ) ) DEALLOCATE( desc%imin3 )
    IF( ALLOCATED( desc%imax3 ) ) DEALLOCATE( desc%imax3 )
    IF( ALLOCATED( desc%npp ) ) DEALLOCATE( desc%npp )
    IF( ALLOCATED( desc%ipp ) ) DEALLOCATE( desc%ipp )
    IF( ALLOCATED( desc%np3 ) ) DEALLOCATE( desc%np3 )
    IF( ALLOCATED( desc%np2 ) ) DEALLOCATE( desc%np2 )
  END SUBROUTINE fft_box_deallocate

!=----------------------------------------------------------------------------=!
!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_set( desc, nat, irb, dfftp )

    USE fft_types

    IMPLICIT NONE

    TYPE (fft_box_descriptor), INTENT(INOUT) :: desc
    TYPE (fft_type_descriptor), INTENT(IN) :: dfftp

    INTEGER, INTENT(in) :: nat
    INTEGER, INTENT(in) :: irb( :, : )

    INTEGER :: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
    INTEGER :: isa
    INTEGER :: ir2, ibig2, irb2, imin2, imax2, nr2
    INTEGER :: ir3, ibig3, irb3, imin3, imax3, nr3

    IF( nat > size( desc%irb, 2 ) ) THEN
       WRITE( stdout, fmt="( ///,'NAT, SIZE = ',2I10)" ) nat, size( desc%irb, 2)
       CALL fftx_error__(" fft_box_set ", " inconsistent dimensions ", 1 )
    ENDIF

    if ( (desc%nr1.eq.0)  .OR. (desc%nr2.eq.0)  .OR. (desc%nr3.eq.0) .OR. &
         (desc%nr1x.eq.0) .OR. (desc%nr2x.eq.0) .OR. (desc%nr3x.eq.0) ) &
         call fftx_error__(" fft_box_set ", "descriptor dimensions must be already initialized", 1)
    nr1b  = desc%nr1  ; nr2b  = desc%nr2  ; nr3b  = desc%nr3
    nr1bx = desc%nr1x ; nr2bx = desc%nr2x ; nr3bx = desc%nr3x

    desc%irb( 1:3, 1:nat ) = irb( 1:3, 1:nat )


    nr3   = sum( dfftp%npp( 1:desc%nproc ) )

    DO isa = 1, nat

       imin3 = nr3b
       imax3 = 1
       irb3  = irb( 3, isa )

       DO ir3 = 1, nr3b
          ibig3 = 1 + mod( irb3 + ir3 - 2, nr3 )
          IF( ibig3 < 1 .or. ibig3 > nr3 )   &
        &        CALL fftx_error__(' fft_box_set ',' ibig3 wrong ', ibig3 )
          ibig3 = ibig3 - dfftp%ipp( desc%mype + 1 )
          IF ( ibig3 > 0 .and. ibig3 <= dfftp%npp(desc%mype + 1) ) THEN
               imin3 = min( imin3, ir3 )
               imax3 = max( imax3, ir3 )
          ENDIF
       ENDDO

       desc%imin3( isa ) = imin3
       desc%imax3( isa ) = imax3
       desc%np3( isa )   = imax3 - imin3 + 1

    ENDDO

    desc%imin2 = 1
    desc%imax2 = desc%nr2
    desc%np2(:) = desc%nr2

    desc%nnr  = desc%nr1x * desc%nr2x * desc%nr3x

  END SUBROUTINE fft_box_set

END MODULE fft_smallbox_type
