!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE fft_smallbox_type

  USE fft_types, ONLY: fft_type_descriptor

  IMPLICIT NONE

  SAVE

  INTEGER :: stdout = 6

  TYPE fft_box_descriptor
    !
    !  Sub (box) grid descriptor
    !
    INTEGER, ALLOCATABLE :: irb(:,:)  ! the offset of the box corner
    INTEGER, ALLOCATABLE :: imin2(:),imin3(:)  ! the starting index of local yz-plane section
    INTEGER, ALLOCATABLE :: imax2(:),imax3(:)  ! the last index of local yz-plane section
    INTEGER, ALLOCATABLE :: np2(:),np3(:)    ! number of local yz-plane section for the box fft
    !
    INTEGER :: nr1    = 0  !
    INTEGER :: nr2    = 0  ! effective FFT dimensions of the 3D grid (global)
    INTEGER :: nr3    = 0  ! 
    INTEGER :: nr1x   = 0  ! FFT grids leading dimensions
    INTEGER :: nr2x   = 0  ! dimensions of the arrays for the 3D grid (global)
    INTEGER :: nr3x   = 0  ! may differ from nr1 ,nr2 ,nr3 in order to boost performances
    INTEGER :: nnr    = 0 
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
    ALLOCATE( desc%imin2( nat ), desc%imin3( nat ) )
    ALLOCATE( desc%imax2( nat ), desc%imax3( nat ) )
    ALLOCATE( desc%np2( nat ), desc%np3( nat ) )
    desc%irb = 0
    desc%imin2 = 0; desc%imin3 = 0
    desc%imax2 = 0; desc%imax3 = 0
    desc%np2 = 0; desc%np3 = 0
    desc%mype = mype
    desc%nproc = nproc
    desc%comm = comm
    desc%root = root
  END SUBROUTINE fft_box_allocate

  SUBROUTINE fft_box_deallocate( desc )
    TYPE (fft_box_descriptor) :: desc
    IF( ALLOCATED( desc%irb ) ) DEALLOCATE( desc%irb )
    IF( ALLOCATED( desc%imin2 ) ) DEALLOCATE( desc%imin2 )
    IF( ALLOCATED( desc%imin3 ) ) DEALLOCATE( desc%imin3 )
    IF( ALLOCATED( desc%imax2 ) ) DEALLOCATE( desc%imax2 )
    IF( ALLOCATED( desc%imax3 ) ) DEALLOCATE( desc%imax3 )
    IF( ALLOCATED( desc%np2 ) ) DEALLOCATE( desc%np2 )
    IF( ALLOCATED( desc%np3 ) ) DEALLOCATE( desc%np3 )
  END SUBROUTINE fft_box_deallocate

!=----------------------------------------------------------------------------=!
!=----------------------------------------------------------------------------=!

  SUBROUTINE fft_box_set( desc, nat, irb, dfftp )

    IMPLICIT NONE

    TYPE (fft_box_descriptor), INTENT(INOUT) :: desc
    INTEGER, INTENT(in) :: nat
    INTEGER, INTENT(in) :: irb( :, : )
    TYPE (fft_type_descriptor), INTENT(IN) :: dfftp

    INTEGER :: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
    INTEGER :: isa
    INTEGER :: ir2, ibig2, irb2, imin2, imax2, nr2
    INTEGER :: ir3, ibig3, irb3, imin3, imax3, nr3

    IF( nat > size( desc%irb, 2 ) ) THEN
       WRITE( stdout, fmt="( ///,'NAT, SIZE = ',2I10)" ) nat, size( desc%irb, 2 )
       CALL fftx_error__(" fft_box_set ", " inconsistent dimensions ", 1 )
    ENDIF

    if ( (desc%nr1.eq.0)  .OR. (desc%nr2.eq.0)  .OR. (desc%nr3.eq.0) .OR. &
         (desc%nr1x.eq.0) .OR. (desc%nr2x.eq.0) .OR. (desc%nr3x.eq.0) ) &
         call fftx_error__(" fft_box_set ", "descriptor dimensions must be already initialized", 1)
    nr1b  = desc%nr1  ; nr2b  = desc%nr2  ; nr3b  = desc%nr3
    nr1bx = desc%nr1x ; nr2bx = desc%nr2x ; nr3bx = desc%nr3x

    desc%irb( 1:3, 1:nat ) = irb( 1:3, 1:nat )

    DO isa = 1, nat

       imin3 = nr3b
       imax3 = 1
       irb3  = irb( 3, isa )

       DO ir3 = 1, nr3b
          ibig3 = 1 + mod( irb3 + ir3 - 2, dfftp%nr3 )
          IF( ibig3 < 1 .or. ibig3 > dfftp%nr3 )   &
        &        CALL fftx_error__(' fft_box_set ',' ibig3 wrong ', ibig3 )
          ibig3 = ibig3 - dfftp%my_i0r3p
          IF ( ibig3 > 0 .and. ibig3 <= dfftp%my_nr3p ) THEN
               imin3 = min( imin3, ir3 )
               imax3 = max( imax3, ir3 )
          ENDIF
       ENDDO

       desc%imin3( isa ) = imin3
       desc%imax3( isa ) = imax3
       desc%np3( isa )   = imax3 - imin3 + 1

       imin2 = nr2b
       imax2 = 1
       irb2  = irb( 2, isa )

       DO ir2 = 1, nr2b
          ibig2 = 1 + mod( irb2 + ir2 - 2, dfftp%nr2 )
          IF( ibig2 < 1 .or. ibig2 > dfftp%nr2 )   &
        &        CALL fftx_error__(' fft_box_set ',' ibig2 wrong ', ibig2 )
          ibig2 = ibig2 - dfftp%my_i0r2p
          IF ( ibig2 > 0 .and. ibig2 <= dfftp%my_nr2p ) THEN
               imin2 = min( imin2, ir2 )
               imax2 = max( imax2, ir2 )
          ENDIF
       ENDDO

       desc%imin2( isa ) = imin2
       desc%imax2( isa ) = imax2
       desc%np2( isa )   = imax2 - imin2 + 1

    ENDDO

    desc%nnr  = desc%nr1x * desc%nr2x * desc%nr3x

  END SUBROUTINE fft_box_set

END MODULE fft_smallbox_type
