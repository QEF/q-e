!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni
!     Last update Apr. 2009
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel
!
   IMPLICIT NONE
   SAVE

   INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

   PRIVATE :: DP
!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn, dtgs )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
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
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,    ONLY : task_groups_descriptor

  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f( : )  ! array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  TYPE (task_groups_descriptor), OPTIONAL, INTENT(in) :: dtgs
                                           ! specify if you want to use task groups parallelization
  !
  ! the following ifdef prevents usage of directive in older ifort versions
  !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
  ! the following is a workaround for Intel 12.1 bug
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: yf, aux
#endif
#endif
#endif
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: yf(:), aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  LOGICAL                    :: use_tg
  !
  !
  IF( present( dtgs ) ) THEN
     use_tg = dtgs%have_task_groups
  ELSE
     use_tg = .false.
  ENDIF
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( use_tg ) THEN
     ALLOCATE( aux( dtgs%nogrp * dtgs%tg_nnr ) )
     ALLOCATE( YF ( dtgs%nogrp * dtgs%tg_nnr ) )
  ELSE
     ALLOCATE( aux( dfft%nnr ) )
  ENDIF
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        IF( use_tg ) &
           CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
        !
        CALL cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        IF( use_tg ) THEN
           CALL pack_group_sticks( f, yf, dtgs )
           CALL cft_1z( yf, dtgs%tg_nsw( me_p ), n3, nx3, isgn, aux )
        ELSE
           CALL cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        ENDIF
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL fw_scatter( isgn ) ! forward scatter from stick to planes
     !
     IF( use_tg ) THEN
        CALL cft_2xy( f, dtgs%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        CALL cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ENDIF
     !
  ELSE
     !
     IF ( isgn /= -2 ) THEN
        !
        IF( use_tg ) &
           CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        planes = dfft%iplw
        !
     ENDIF

     IF( use_tg ) THEN
        CALL cft_2xy( f, dtgs%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        CALL cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes)
     ENDIF
     !
     CALL bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
         !
     ELSE
        !
        IF( use_tg ) THEN
           CALL cft_1z( aux, dtgs%tg_nsw( me_p ), n3, nx3, isgn, yf )
           CALL unpack_group_sticks( yf, f, dtgs )
        ELSE
           CALL cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  IF( use_tg ) THEN
     DEALLOCATE( yf )
  ENDIF

  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
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
     USE scatter_mod, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        IF( use_tg ) THEN
           !
           CALL fft_scatter( dfft, aux, nx3, dtgs%nogrp*dtgs%tg_nnr, f, dtgs%tg_nsw, dtgs%tg_npp, iopt, dtgs )
           !
        ELSE
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        ENDIF
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE fw_scatter

  !

  SUBROUTINE bw_scatter( iopt )
     !
     USE scatter_mod, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        IF( use_tg ) THEN
           !
           CALL fft_scatter( dfft, aux, nx3, dtgs%nogrp*dtgs%tg_nnr, f, dtgs%tg_nsw, dtgs%tg_npp, iopt, dtgs )
           !
        ELSE
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        ENDIF
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE bw_scatter
  !
END SUBROUTINE tg_cft3s
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_z( f_in, dfft, f_out, dtgs )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_1z
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f_in( : )  ! INPUT array containing data to be transformed
  COMPLEX(DP), INTENT(inout)   :: f_out (:)  ! OUTPUT
  TYPE (fft_type_descriptor), INTENT(in) :: dfft ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  !
  CALL cft_1z( f_in, dtgs%tg_nsw( dtgs%mype + 1 ), dfft%nr3, dfft%nr3x, 2, f_out )
  !
END SUBROUTINE fw_tg_cft3_z
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_z( f_out, dfft, f_in, dtgs )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_1z
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f_out( : ) ! OUTPUT
  COMPLEX(DP), INTENT(inout)   :: f_in (:) ! INPUT array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  !
  CALL cft_1z( f_in, dtgs%tg_nsw( dtgs%mype + 1 ), dfft%nr3, dfft%nr3x, -2, f_out )
  !
END SUBROUTINE bw_tg_cft3_z
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_scatter( f, dfft, aux, dtgs )
  !----------------------------------------------------------------------------
  !
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ), aux( : )  ! array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft     ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  !
  CALL fft_scatter( dfft, aux, dfft%nr3x, dtgs%nogrp*dtgs%tg_nnr, f, dtgs%tg_nsw, dtgs%tg_npp, 2, dtgs )
  !
END SUBROUTINE fw_tg_cft3_scatter
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_scatter( f, dfft, aux, dtgs )
  !----------------------------------------------------------------------------
  !
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ), aux( : )  ! array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft     ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  !
  CALL fft_scatter( dfft, aux, dfft%nr3x, dtgs%nogrp*dtgs%tg_nnr, f, dtgs%tg_nsw, dtgs%tg_npp, -2, dtgs )
  !
END SUBROUTINE bw_tg_cft3_scatter
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_xy( f, dfft, dtgs )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_2xy
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ) ! INPUT/OUTPUT array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  INTEGER                    :: planes( dfft%nr1x )
  !
  planes = dfft%iplw
  CALL cft_2xy( f, dtgs%tg_npp( dtgs%mype + 1 ), dfft%nr1, dfft%nr2, dfft%nr1x, dfft%nr2x, 2, planes )
  !
END SUBROUTINE fw_tg_cft3_xy
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_xy( f, dfft, dtgs )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_2xy
  USE fft_types,  ONLY : fft_type_descriptor
  USE task_groups,  ONLY : task_groups_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ) ! INPUT/OUTPUT  array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft ! descriptor of fft data layout
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs ! descriptor of fft data layout
  INTEGER                    :: planes( dfft%nr1x )
  !
  planes = dfft%iplw
  CALL cft_2xy( f, dtgs%tg_npp( dtgs%mype + 1 ), dfft%nr1, dfft%nr2, dfft%nr1x, dfft%nr2x, -2, planes )
  !
END SUBROUTINE bw_tg_cft3_xy

#if defined(__DOUBLE_BUFFER)
  SUBROUTINE pack_group_sticks_i( f, yf, dtgs, req)

     USE task_groups,  ONLY : task_groups_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(in)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(out)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (task_groups_descriptor), INTENT(in) :: dtgs
     INTEGER                     :: ierr,req
     !
     CALL start_clock( 'IALLTOALL' )
     !
     !  Collect all the sticks of the different states,
     !  in "yf" processors will have all the sticks of the OGRP

#if defined(__MPI)

     CALL MPI_IALLTOALLV( f(1), dtgs%tg_snd, dtgs%tg_psdsp, MPI_DOUBLE_COMPLEX, yf(1), dtgs%tg_rcv, &
      &                     dtgs%tg_rdsp, MPI_DOUBLE_COMPLEX, dtgs%ogrp_comm, req, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'pack_group_sticks_i', ' alltoall error 1 ', abs(ierr) )
     ENDIF

#else

     IF( dfft%tg_rcv(dfft%nogrp) /= dfft%tg_snd(dfft%nogrp) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 3 )
     ENDIF

     yf( 1 : dfft%tg_rcv(dfft%nogrp) ) =  f( 1 : dfft%tg_snd(dfft%nogrp) )

#endif

     CALL stop_clock( 'IALLTOALL' )
     !
     !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
     !
     RETURN
  END SUBROUTINE pack_group_sticks_i
#endif

!----------------------------------------------------------------------------
  SUBROUTINE pack_group_sticks( f, yf, dtgs )

     USE task_groups,  ONLY : task_groups_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(in)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(out)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (task_groups_descriptor), INTENT(in) :: dtgs
     INTEGER                     :: ierr
     !
     IF( dtgs%tg_rdsp(dtgs%nogrp) + dtgs%tg_rcv(dtgs%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'pack_group_sticks' , ' inconsistent size ', 1 )
     ENDIF
     IF( dtgs%tg_psdsp(dtgs%nogrp) + dtgs%tg_snd(dtgs%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 2 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )
     !
     !  Collect all the sticks of the different states,
     !  in "yf" processors will have all the sticks of the OGRP

#if defined(__MPI)

     CALL MPI_ALLTOALLV( f(1), dtgs%tg_snd, dtgs%tg_psdsp, MPI_DOUBLE_COMPLEX, yf(1), dtgs%tg_rcv, &
      &                     dtgs%tg_rdsp, MPI_DOUBLE_COMPLEX, dtgs%ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' alltoall error 1 ', abs(ierr) )
     ENDIF

#else

     IF( dtgs%tg_rcv(dtgs%nogrp) /= dtgs%tg_snd(dtgs%nogrp) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 3 )
     ENDIF

     yf( 1 : dtgs%tg_rcv(dtgs%nogrp) ) =  f( 1 : dtgs%tg_snd(dtgs%nogrp) )

#endif

     CALL stop_clock( 'ALLTOALL' )
     !
     !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
     !
     RETURN
  END SUBROUTINE pack_group_sticks

  !
  SUBROUTINE unpack_group_sticks( yf, f, dtgs )

     USE task_groups,  ONLY : task_groups_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(out)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(in)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (task_groups_descriptor), INTENT(in) :: dtgs
     !
     !  Bring pencils back to their original distribution
     !
     INTEGER                     :: ierr
     !
     IF( dtgs%tg_usdsp(dtgs%nogrp) + dtgs%tg_snd(dtgs%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 3 )
     ENDIF
     IF( dtgs%tg_rdsp(dtgs%nogrp) + dtgs%tg_rcv(dtgs%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 4 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )

#if defined(__MPI)
     CALL MPI_Alltoallv( yf(1), &
          dtgs%tg_rcv, dtgs%tg_rdsp, MPI_DOUBLE_COMPLEX, f(1), &
          dtgs%tg_snd, dtgs%tg_usdsp, MPI_DOUBLE_COMPLEX, dtgs%ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' alltoall error 2 ', abs(ierr) )
     ENDIF
#endif

     CALL stop_clock( 'ALLTOALL' )

     RETURN
  END SUBROUTINE unpack_group_sticks


SUBROUTINE tg_gather( dffts, dtgs, v, tg_v )
   !
   USE fft_types,      ONLY : fft_type_descriptor
   USE task_groups,    ONLY : task_groups_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

   TYPE(fft_type_descriptor), INTENT(in) :: dffts
   TYPE(task_groups_descriptor), INTENT(in) :: dtgs

   REAL(DP) :: v(:)
   REAL(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dtgs%nogrp ), recv_displ( dtgs%nogrp )

   nsiz_tg = dtgs%tg_nnr * dtgs%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dtgs%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dtgs%nogrp
      recv_cnt(i) = dffts%npp( dtgs%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dtgs%nogrp) + recv_cnt( dtgs%nogrp ) + 1, size( tg_v )
      tg_v( i ) = 0.0d0
   ENDDO

#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dtgs%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_gather
!
!  Complex version of previous routine
!
SUBROUTINE tg_cgather( dffts, dtgs, v, tg_v )
   !
   USE fft_types,      ONLY : fft_type_descriptor
   USE task_groups,    ONLY : task_groups_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif

   TYPE(fft_type_descriptor), INTENT(in) :: dffts
   TYPE(task_groups_descriptor), INTENT(in) :: dtgs

   COMPLEX(DP) :: v(:)
   COMPLEX(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dtgs%nogrp ), recv_displ( dtgs%nogrp )

   nsiz_tg = dtgs%tg_nnr * dtgs%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dtgs%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dtgs%nogrp
      recv_cnt(i) = dffts%npp( dtgs%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dtgs%nogrp) + recv_cnt( dtgs%nogrp ) + 1, size( tg_v )
      tg_v( i ) = (0.0d0,0.0d0)
   ENDDO
   !
   ! The quantities are complex, multiply the cunters by 2 and gather
   ! real numbers
   !
   nsiz = 2 * nsiz
   recv_cnt = 2 * recv_cnt
   recv_displ = 2 * recv_displ
#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dtgs%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_cgather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_cgather

!--------------------------------------------------------------------------------
!   Auxiliary routines to read/write from/to a distributed array
!   NOT optimized for efficiency .... just to show how one can access the data
!--------------------------------------------------------------------------------
!
COMPLEX (DP) FUNCTION get_f_of_R (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  INTEGER :: kk, ii, jj, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do jj = 1, dfft%nproc
     if ( dfft%ipp(jj) < k ) kk = jj
  end do
  ii  = i + dfft%nr1x * ( j - 1 ) + dfft%nr1x * dfft%nr2x * ( k - dfft%ipp(kk) - 1 )
  f_aux = (0.d0,0.d0)
  if (kk == dfft%mype +1) f_aux = f(ii)
  CALL MPI_ALLREDUCE( f_aux, get_f_of_R,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  get_f_of_R = f(ii)
#endif
END FUNCTION get_f_of_R

SUBROUTINE put_f_of_R (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f_in
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  INTEGER :: kk, ii, jj, ierr

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do jj = 1, dfft%nproc
     if ( dfft%ipp(jj) < k ) kk = jj
  end do
  ii  = i + dfft%nr1x * ( j - 1 ) + dfft%nr1x * dfft%nr2x * ( k - dfft%ipp(kk) - 1 )
  if (kk == dfft%mype +1) f(ii) = f_in
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  f(ii) = f_in
#endif

END SUBROUTINE put_f_of_R

COMPLEX (DP) FUNCTION get_f_of_G (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  f_aux = (0.d0,0.d0)
  if ( jj > 0 ) f_aux = f( k + dfft%nr3x * (jj -1))
  CALL MPI_ALLREDUCE( f_aux, get_f_of_G,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  get_f_of_G = f(ii)
#endif
END FUNCTION get_f_of_G

SUBROUTINE put_f_of_G (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  COMPLEX(DP), INTENT (IN) :: f_in
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  if ( jj > 0 )   f( k + dfft%nr3x * (jj -1)) = f_in
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  f(ii) = f_in
#endif
END SUBROUTINE put_f_of_G

END MODULE fft_parallel
