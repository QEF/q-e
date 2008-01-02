!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!=---------------------------------------------------------------------==!
!
!
!     Parallel 3D FFT high level Driver 
!     ( Charge density and Wave Functions )
!
!     Carlo Cavazzoni, Dec. 2007
!
!=---------------------------------------------------------------------==!
!
MODULE fft_parallel
!
   IMPLICIT NONE
   SAVE
!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn, use_task_groups )
  !----------------------------------------------------------------------------
  !
  ! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  ! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  ! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  ! ...              fft along z using pencils        (cft_1z)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code 
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE fft_base,   ONLY : fft_scatter
  USE kinds,      ONLY : DP
  USE mp_global,  only : me_pool, nproc_pool, ogrp_comm, npgrp, nogrp, &
                         intra_pool_comm, nolist, nplist
  USE fft_types,  ONLY : fft_dlay_descriptor
  USE parallel_include

  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT)    :: f( : )  ! array containing data to be transformed
  type (fft_dlay_descriptor), intent(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(IN)           :: isgn    ! fft direction
  LOGICAL, OPTIONAL, INTENT(IN) :: use_task_groups 
                                           ! specify if you want to use task groups parallelization
  !
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: yf(:), aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  LOGICAL                    :: use_tg
  !
  CALL start_clock( 'cft3s' )
  !
  IF( PRESENT( use_task_groups ) ) THEN
     use_tg = use_task_groups
  ELSE
     use_tg = .FALSE.
  END IF
  !
  IF( use_tg .AND. .NOT. dfft%have_task_groups ) &
     CALL errore( ' tg_cft3s ', ' call requiring task groups for a descriptor without task groups ', 1 )
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( use_tg ) THEN
     ALLOCATE( aux( nogrp * dfft%nnrx ) )
     ALLOCATE( YF ( nogrp * dfft%nnrx ) )
  ELSE
     ALLOCATE( aux( dfft%nnrx ) )
  END IF
  !
  me_p = me_pool + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        IF( use_tg ) &
           CALL errore( ' tg_cfft ', ' task groups on large mesh not implemented ', 1 )
        !
        call cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL pack_group_sticks()
        !
        IF( use_tg ) THEN
           CALL cft_1z( yf, dfft%tg_nsw( me_p ), n3, nx3, isgn, aux )
        ELSE
           call cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        END IF
        !
        planes = dfft%iplw
        !
     END IF
     !
     CALL fw_scatter( isgn ) ! forwart scatter from stick to planes
     !
     IF( use_tg ) THEN
        CALL cft_2xy( f, dfft%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        CALL cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     END IF
     !
  ELSE
     !
     if ( isgn /= -2 ) then
        !
        IF( use_tg ) &
           CALL errore( ' tg_cfft ', ' task groups on large mesh not implemented ', 1 )
        !
        planes = dfft%iplp
        !
     else
        !
        planes = dfft%iplw
        !
     endif

     IF( use_tg ) THEN
        CALL cft_2xy( f, dfft%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        call cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes)
     END IF
     !
     CALL bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        call cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
        !
     ELSE
        !
        IF( use_tg ) THEN
           CALL cft_1z( aux, dfft%tg_nsw( me_p ), n3, nx3, isgn, yf )
        ELSE
           call cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        END IF
        !
        CALL unpack_group_sticks()
        !
     END IF
     !
  END IF
  !
  DEALLOCATE( aux )
  !
  IF( use_tg ) THEN
     DEALLOCATE( yf )
  END IF
  !
  CALL stop_clock( 'cft3s' )
  !
  RETURN
  !
CONTAINS
  !

  SUBROUTINE pack_group_sticks()

     INTEGER                     :: idx, ierr
     INTEGER, DIMENSION(nogrp+1) :: send_cnt, send_displ, recv_cnt, recv_displ
     !
     if( .NOT. use_tg ) return
     !
     send_cnt(1)   = nx3 * dfft%nsw( me_p )
     IF( nx3 * dfft%nsw( me_p ) > dfft%nnrx ) THEN
        CALL errore( ' tg_cfft ', ' inconsistent dfft%nnrx ', 1 )
     END IF
     send_displ(1) = 0
     recv_cnt(1)   = nx3 * dfft%nsw( nolist(1) + 1 )
     recv_displ(1) = 0
     DO idx = 2, nogrp
        send_cnt(idx)   = nx3 * dfft%nsw( me_p )
        send_displ(idx) = send_displ(idx-1) + dfft%nnrx
        recv_cnt(idx)   = nx3 * dfft%nsw( nolist(idx) + 1 )
        recv_displ(idx) = recv_displ(idx-1) + recv_cnt(idx-1)
     ENDDO

     IF( recv_displ(nogrp) + recv_cnt(nogrp) > SIZE( yf ) ) THEN
        CALL errore( ' tg_cfft ', ' inconsistent size ', 1 )
     END IF
     IF( send_displ(nogrp) + send_cnt(nogrp) > SIZE( f ) ) THEN
        CALL errore( ' tg_cfft ', ' inconsistent size ', 2 )
     END IF

     CALL start_clock( 'ALLTOALL' )
     !
     !  Collect all the sticks of the different states,
     !  in "yf" processors will have all the sticks of the OGRP

#if defined __MPI

     CALL MPI_ALLTOALLV( f(1), send_cnt, send_displ, MPI_DOUBLE_COMPLEX, yf(1), recv_cnt, &
      &                     recv_displ, MPI_DOUBLE_COMPLEX, ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL errore( ' tg_cfft ', ' alltoall error 1 ', ABS(ierr) )
     END IF

#endif

     CALL stop_clock( 'ALLTOALL' )
     !
     !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
     !
     RETURN
  END SUBROUTINE pack_group_sticks

  !

  SUBROUTINE unpack_group_sticks()
     !
     !  Bring pencils back to their original distribution
     !
     INTEGER                     :: idx, ierr
     INTEGER, DIMENSION(nogrp+1) :: send_cnt, send_displ, recv_cnt, recv_displ
     !
     if( .NOT. use_tg ) return
     !
     send_cnt  (1) = nx3 * dfft%nsw( nolist(1) + 1 )
     send_displ(1) = 0
     recv_cnt  (1) = nx3 * dfft%nsw( me_p )
     recv_displ(1) = 0
     DO idx = 2, NOGRP
        send_cnt  (idx) = nx3 * dfft%nsw( nolist(idx) + 1 )
        send_displ(idx) = send_displ(idx-1) + send_cnt(idx-1)
        recv_cnt  (idx) = nx3 * dfft%nsw( me_p )
        recv_displ(idx) = recv_displ(idx-1) + recv_cnt(idx-1)
     ENDDO

     IF( recv_displ(nogrp) + recv_cnt(nogrp) > SIZE( f ) ) THEN
        CALL errore( ' tg_cfft ', ' inconsistent size ', 3 )
     END IF
     IF( send_displ(nogrp) + send_cnt(nogrp) > SIZE( yf ) ) THEN
        CALL errore( ' tg_cfft ', ' inconsistent size ', 4 )
     END IF

     CALL start_clock( 'ALLTOALL' )

#if defined __MPI
     CALL MPI_Alltoallv( yf(1), &
          send_cnt, send_displ, MPI_DOUBLE_COMPLEX, f(1), &
          recv_cnt, recv_displ, MPI_DOUBLE_COMPLEX, ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL errore( ' tg_cfft ', ' alltoall error 2 ', ABS(ierr) )
     END IF
#endif   

     CALL stop_clock( 'ALLTOALL' )

     RETURN
  END SUBROUTINE unpack_group_sticks

  !

  SUBROUTINE fw_scatter( iopt )

        !Transpose data for the 2-D FFT on the x-y plane
        !
        !NOGRP*dfft%nnr: The length of aux and f
        !nr3x: The length of each Z-stick 
        !aux: input - output
        !f: working space
        !isgn: type of scatter
        !dfft%nsw(me) holds the number of Z-sticks proc. me has.
        !dfft%npp: number of planes per processor
        !
     !
     use fft_base, only: fft_scatter
     !
     INTEGER, INTENT(IN) :: iopt
     INTEGER :: nppx, ip, nnp, npp, ii, i, mc, j
     !
     IF( iopt == 2 ) THEN
        !
        IF( use_tg ) THEN
           !
           nppx = dfft%tg_npp( me_p )
           npp  = dfft%tg_npp( me_p )
           nnp  = nx1*nx2
           !
           CALL fft_scatter( aux, nx3, nogrp*dfft%nnrx, f, dfft%tg_nsw, dfft%tg_npp, iopt, use_tg )
           !
        ELSE
           !
           nppx = dfft%npp( me_p )
           IF( nproc_pool == 1 ) nppx = dfft%nr3x
           npp  = dfft%npp( me_p )
           nnp  = dfft%nnp
           !
           call fft_scatter( aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        END IF
        !
        f(:) = (0.d0, 0.d0)
        ii = 0
        !
        do ip = 1, nproc_pool
           !
           do i = 1, dfft%nsw( ip )
              !
              mc = dfft%ismap( i + dfft%iss( ip ) )
              !
              ii = ii + 1
              !
              do j = 1, npp
                 f( mc + ( j - 1 ) * nnp ) = aux( j + ( ii - 1 ) * nppx )
              end do
              !
           end do
           !
        end do
        !
     ELSE IF( iopt == 1 ) THEN
        !
        if ( nproc_pool == 1 ) then
           nppx = dfft%nr3x
        else
           nppx = dfft%npp( me_p )
        end if
        !
        call fft_scatter( aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
        f(:) = (0.d0, 0.d0)
        !
        do i = 1, dfft%nst
           mc = dfft%ismap( i )
           do j = 1, dfft%npp( me_p )
              f( mc + ( j - 1 ) * dfft%nnp ) = aux( j + ( i - 1 ) * nppx )
           end do
        end do
        !
     END IF
     !
     RETURN
  END SUBROUTINE fw_scatter

  !

  SUBROUTINE bw_scatter( iopt )
     !
     use fft_base, only: fft_scatter
     !
     INTEGER, INTENT(IN) :: iopt
     INTEGER :: nppx, ip, nnp, npp, ii, i, mc, j
     !
     !  
     IF( iopt == -2 ) THEN
        !  
        IF( use_tg ) THEN
           !
           nppx = dfft%tg_npp( me_p )
           npp  = dfft%tg_npp( me_p )
           nnp  = nx1*nx2
           !
        ELSE
           !
           nppx = dfft%npp( me_p )
           IF( nproc_pool == 1 ) nppx = dfft%nr3x
           npp  = dfft%npp( me_p )
           nnp  = dfft%nnp
           !
        END IF

        ii = 0

        do ip = 1, nproc_pool
           do i = 1, dfft%nsw( ip )
              mc = dfft%ismap( i + dfft%iss( ip ) )
              ii = ii + 1
              do j = 1, npp
                 aux( j + ( ii - 1 ) * nppx ) = f( mc + ( j - 1 ) * nnp )
              end do
           end do
        end do
        !
        IF( use_tg ) THEN
           !
           CALL fft_scatter( aux, nx3, nogrp*dfft%nnrx, f, dfft%tg_nsw, dfft%tg_npp, iopt, use_tg )
           !
        ELSE
           !
           call fft_scatter( aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        END IF
        ! 
     ELSE IF( iopt == -1 ) THEN
        !
        if ( nproc_pool == 1 ) then
           nppx = dfft%nr3x
        else
           nppx = dfft%npp( me_p )
        end if
        do i = 1, dfft%nst
           mc = dfft%ismap( i )
           do j = 1, dfft%npp( me_p )
              aux( j + ( i - 1 ) * nppx ) = f( mc + ( j - 1 ) * dfft%nnp )
           end do
        end do
        call fft_scatter( aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
     END IF
     !
     RETURN
  END SUBROUTINE bw_scatter
  !
END SUBROUTINE tg_cft3s
!
END MODULE fft_parallel
