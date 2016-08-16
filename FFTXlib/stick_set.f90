!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=
   MODULE stick_set
!=----------------------------------------------------------------------=

!  ... Distribute G-vectors across processors as sticks and planes,
!  ... initialize FFT descriptors for both dense and smooth grids

      USE stick_base
      USE fft_types

      IMPLICIT NONE

      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

      PRIVATE
      SAVE

      TYPE(sticks_map) :: smap

      PUBLIC :: pstickset,  pstickdealloc
      PUBLIC :: smap

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE pstickdealloc()
         CALL sticks_map_deallocate( smap )
      END SUBROUTINE pstickdealloc


!=----------------------------------------------------------------------=

      SUBROUTINE pstickset( lgamma, dfftp, dffts, ngw, ngm, ngs, ionode, stdout )

          TYPE(fft_type_descriptor), INTENT(inout) :: dfftp, dffts
          INTEGER, INTENT(out) :: ngw, ngm, ngs
          LOGICAL, INTENT(IN) :: ionode, lgamma
          INTEGER, INTENT(IN) :: stdout
          !
          !  Initialize task groups descriptor
          !
          ngw = dffts%nwl( dffts%mype + 1 )
          ngs = dffts%ngl( dffts%mype + 1 )
          ngm = dfftp%ngl( dfftp%mype + 1 )
          !
          IF( lgamma ) THEN
             ngw = (ngw + 1)/2
             ngs = (ngs + 1)/2
             ngm = (ngm + 1)/2
          END IF
          !
          IF (ionode) THEN
             WRITE( stdout,*)
             IF ( dfftp%nproc > 1 ) THEN
                WRITE( stdout, '(5X,"Parallelization info")')
             ELSE
                WRITE( stdout, '(5X,"G-vector sticks info")')
             ENDIF
             WRITE( stdout, '(5X,"--------------------")')
             WRITE( stdout, '(5X,"sticks:   dense  smooth     PW", &
                            & 5X,"G-vecs:    dense   smooth      PW")') 
             IF ( dfftp%nproc > 1 ) THEN
                WRITE( stdout,'(5X,"Min",4X,2I8,I7,12X,2I9,I8)') &
                   minval(dfftp%nsp), minval(dffts%nsp), minval(dffts%nsw), &
                   minval(dfftp%ngl), minval(dffts%ngl), minval(dffts%nwl)
                WRITE( stdout,'(5X,"Max",4X,2I8,I7,12X,2I9,I8)') &
                   maxval(dfftp%nsp), maxval(dffts%nsp), maxval(dffts%nsw), &
                   maxval(dfftp%ngl), maxval(dffts%ngl), maxval(dffts%nwl)
             END IF
             WRITE( stdout,'(5X,"Sum",4X,2I8,I7,12X,2I9,I8)') &
                sum(dfftp%nsp), sum(dffts%nsp), sum(dffts%nsw), &
                sum(dfftp%ngl), sum(dffts%ngl), sum(dffts%nwl)
          ENDIF

          IF(ionode) WRITE( stdout,*)

          RETURN
        END SUBROUTINE pstickset



!=----------------------------------------------------------------------=
   END MODULE stick_set
!=----------------------------------------------------------------------=
