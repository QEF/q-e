!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  ----------------------------------------------
!  This Module written by Carlo Cavazzoni 
!  Last modified April 2003
!  ----------------------------------------------

#include "f_defs.h"

!=---------------------------------------------------------------------==!
!
!
!     FFT high level Driver 
!     ( Charge density and Wave Functions )
!
!
!=---------------------------------------------------------------------==!
!

MODULE fft_cp

  USE fft_types, ONLY: fft_dlay_descriptor

  IMPLICIT NONE
  SAVE

  INTEGER, PRIVATE :: what_scatter = 1

CONTAINS

!----------------------------------------------------------------------
   SUBROUTINE cfft_cp ( f, nr1, nr2, nr3, nr1x, nr2x, nr3x, sign, dfft )
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!   sign = +-2 : parallel 3d fft for wavefunctions
!
!   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!              fft along z using pencils        (cft_1z)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along y (using planes) and x (cft_2xy)
!   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega 
!              fft along x and y(using planes)  (cft_2xy)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along z using pencils        (cft_1z)
!
!   The array "planes" signals whether a fft is needed along y :
!     planes(i)=0 : column f(i,*,*) empty , don't do fft along y
!     planes(i)=1 : column f(i,*,*) filled, fft along y needed
!   "empty" = no active components are present in f(i,*,*) 
!             after (sign>0) or before (sign<0) the fft on z direction
!
!   Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
!   and all planes(i) are set to 1
!
!   based on code written by Stefano de Gironcoli for PWSCF
!
      use mp_global, only: mpime, nproc, group
      use fft_scalar, only: cft_1z, cft_2xy
!
      implicit none
      !
      integer, intent(in) :: nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      type (fft_dlay_descriptor), intent(in) :: dfft
      complex(8) :: f( dfft%nnr )
      complex(8), allocatable :: aux( : )

      integer  mc, i, j, ii, proc, k, me
      integer planes(nr1x)
!
! see comments in cfftp for the logic (or lack of it) of the following
!
      if ( nr1  /= dfft%nr1  ) call errore(' cfft ',' wrong dims ', 1)
      if ( nr2  /= dfft%nr2  ) call errore(' cfft ',' wrong dims ', 2)
      if ( nr3  /= dfft%nr3  ) call errore(' cfft ',' wrong dims ', 3)
      if ( nr1x /= dfft%nr1x ) call errore(' cfft ',' wrong dims ', 4)
      if ( nr2x /= dfft%nr2x ) call errore(' cfft ',' wrong dims ', 5)
      if ( nr3x /= dfft%nr3x ) call errore(' cfft ',' wrong dims ', 6)

      me = mpime + 1

      allocate( aux( dfft%nnr ) )

      if ( sign > 0 ) then
         !
         if ( sign /= 2 ) then
            call cft_1z( f, dfft%nsp(me), nr3, nr3x, sign, aux )
            CALL fw_scatter( sign ) ! forwart scatter from stick to planes
            planes = dfft%iplp
         else
            call cft_1z( f, dfft%nsw(me), nr3, nr3x, sign, aux )
            CALL fw_scatter( sign ) ! forwart scatter from stick to planes
            planes = dfft%iplw
         end if
         !
         call cft_2xy( f, dfft%npp( me ), nr1, nr2, nr1x, nr2x, sign, planes )
         !
      else
         !
         if ( sign .ne. -2 ) then
            planes = dfft%iplp
         else
            planes = dfft%iplw
         endif
         !
         call cft_2xy( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
         !
         if ( sign /= -2 ) then
            call bw_scatter( sign )
            call cft_1z( aux, dfft%nsp( me ), nr3, nr3x, sign, f )
         else
            call bw_scatter( sign )
            call cft_1z( aux, dfft%nsw( me ), nr3, nr3x, sign, f )
         end if
         !
      end if
!
      deallocate( aux )

      RETURN
      !
      !
   CONTAINS
      !
      !
      SUBROUTINE fw_scatter( iopt )
         !
         use fft_base, only: fft_scatter, fft_transpose, fft_itranspose
         !
         INTEGER, INTENT(IN) :: iopt
         INTEGER :: nppx
         !
         !
         IF( iopt == 2 ) THEN
            !
            IF( what_scatter == 1 ) THEN
               call fft_transpose ( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, -2)
            ELSE IF( what_scatter == 2 ) THEN
               call fft_itranspose( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, -2)
            ELSE 
               if ( nproc == 1 ) then
                  nppx = dfft%nr3x
               else
                  nppx = dfft%npp( me )
               end if
               call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
               f(:) = (0.d0, 0.d0)
               ii = 0
               do proc = 1, nproc
                  do i = 1, dfft%nsw( proc )
                     mc = dfft%ismap( i + dfft%iss( proc ) )
                     ii = ii + 1 
                     do j = 1, dfft%npp( me )
                        f( mc + ( j - 1 ) * dfft%nnp ) = aux( j + ( ii - 1 ) * nppx )
                     end do
                  end do
               end do
            END IF
            !
         ELSE IF( iopt == 1 ) THEN
            !
            if ( nproc == 1 ) then
               nppx = dfft%nr3x
            else
               nppx = dfft%npp( me )
            end if
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
            f(:) = (0.d0, 0.d0)
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp( me )
                  f( mc + ( j - 1 ) * dfft%nnp ) = aux( j + ( i - 1 ) * nppx )
               end do
            end do
            !
         END IF
         !
         RETURN 
      END SUBROUTINE fw_scatter
      !
      !  
      !
      SUBROUTINE bw_scatter( iopt )
         !
         use fft_base, only: fft_scatter, fft_transpose, fft_itranspose
         !
         INTEGER, INTENT(IN) :: iopt
         INTEGER :: nppx
         !
         !
         IF( iopt == -2 ) THEN
            !
            IF( what_scatter == 1 ) THEN
               call fft_transpose ( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, 2)
            ELSE IF( what_scatter == 2 ) THEN
               call fft_itranspose( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, 2)
            ELSE 
               if ( nproc == 1 ) then
                  nppx = dfft%nr3x
               else
                  nppx = dfft%npp( me )
               end if
               ii = 0
               do proc = 1, nproc
                  do i = 1, dfft%nsw( proc )
                     mc = dfft%ismap( i + dfft%iss( proc ) )
                     ii = ii + 1 
                     do j = 1, dfft%npp( me )
                        aux( j + ( ii - 1 ) * nppx ) = f( mc + ( j - 1 ) * dfft%nnp )
                     end do
                  end do
               end do
               call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
            END IF
            !
         ELSE IF( iopt == -1 ) THEN
            !
            if ( nproc == 1 ) then
               nppx = dfft%nr3x
            else
               nppx = dfft%npp( me )
            end if
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp( me )
                  aux( j + ( i - 1 ) * nppx ) = f( mc + ( j - 1 ) * dfft%nnp )
               end do
            end do
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
            !
         END IF
         !
         RETURN 
      END SUBROUTINE bw_scatter
      !
      !  
      !
   END SUBROUTINE cfft_cp
   !
   !
   !
   !=======================================================================
   ! ADDED TASK GROUP FFT DRIVERS
   !=======================================================================

   !----------------------------------------------------------------------
   !TASK GROUPS FFT ROUTINE.
   !Added: C. Bekas, Oct. 2005. Adopted from the CPMD code (A. Curioni)
   !----------------------------------------------------------------------

!----------------------------------------------------------------------
    SUBROUTINE tg_cfft_cp ( f, nr1, nr2, nr3, nr1x, nr2x, nr3x, sign, dfft )
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!   sign = +-2 : parallel 3d fft for wavefunctions
!
!   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!              fft along z using pencils        (cft_1)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along y (using planes) and x (cft_2)
!   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
!              fft along x and y(using planes)  (cft_2)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along z using pencils        (cft_1)
!
!   The array "planes" signals whether a fft is needed along y :
!     planes(i)=0 : column f(i,*,*) empty , don't do fft along y
!     planes(i)=1 : column f(i,*,*) filled, fft along y needed
!   "empty" = no active components are present in f(i,*,*)
!             after (sign>0) or before (sign<0) the fft on z direction
!
!   Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
!   and all planes(i) are set to 1
!
!   based on code written by Stefano de Gironcoli for PWSCF
!
      USE mp_global,  only: mpime, nproc, me_ogrp, me_pgrp, npgrp, nogrp
      USE fft_base,   only: fft_scatter, group_fft_scatter
      USE fft_scalar, only: cft_1z, cft_2xy
      USE groups_module
      USE parallel_include

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      TYPE (fft_dlay_descriptor), INTENT(IN) :: dfft

      COMPLEX*16, DIMENSION((NOGRP+1)*strd), INTENT(INOUT) :: f

      integer  mc, i, j, ii, proc, k, nppx, me
      integer planes(nr1x)
      integer group_index, nnz_index, offset, ierr

      !--------------
      !C. Bekas
      !--------------
      INTEGER :: HWMN, IT1, FLAG, num_planes, num_sticks
      COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: XF, YF, aux 
      INTEGER, DIMENSION(NOGRP) :: local_send_cnt, local_send_displ, local_recv_cnt, local_recv_displ
      INTEGER  :: index

      ALLOCATE(XF((NOGRP+1)*strd))
      ALLOCATE(YF((NOGRP+1)*strd))
      ALLOCATE(aux((NOGRP+1)*strd))
 
      !
      ! see comments in cfftp for the logic (or lack of it) of the following
      !
      if ( nr1  /= dfft%nr1  ) call errore(' cfft ',' wrong dims ', 1)
      if ( nr2  /= dfft%nr2  ) call errore(' cfft ',' wrong dims ', 2)
      if ( nr3  /= dfft%nr3  ) call errore(' cfft ',' wrong dims ', 3)
      if ( nr1x /= dfft%nr1x ) call errore(' cfft ',' wrong dims ', 4)
      if ( nr2x /= dfft%nr2x ) call errore(' cfft ',' wrong dims ', 5)
      if ( nr3x /= dfft%nr3x ) call errore(' cfft ',' wrong dims ', 6)

      me = mpime + 1

      if ( nproc == 1 ) then
         nppx = dfft%nr3x
      else
         nppx = dfft%npp(me)
      end if

      !-----------------------
      !Inverse FFT
      !-----------------------

      if ( sign > 0 ) then

         if ( sign /= 2 ) then

            !-----------------------------------
            !Density - Potential FFT calculation
            !-----------------------------------

            !-----------------------------------------
            !ALL TO ALL IN THE ORBITAL GROUP (ME_OGRP)
            !-----------------------------------------
            !
            !Find out how many elements to exchange: Each processor holds dfft%nnr complex fourier
            !coefficients for each eigenvalue. The exchange will move NOGRP*dfft%nnr coefficients
            !between tasks so that each on of the NOGRP groups will hold all necessary coeficients
            !for 1 (one) eigenvalue

#if defined __MPI
            call MPI_ALLTOALL(f, dfft%nnr*16, MPI_BYTE, XF, dfft%nnr*16, MPI_BYTE, ME_OGRP, IERR)
#endif

            !-----------------------------------------------------------------------------------------
            !ADDED COMMENTS: C. Bekas, Oct. 2005
            !-----------------------------------SUBROUTINE: CFT_1Z: Sequence of FFTs
            !XF           : holds the data to be transformed: Pencils in the z-direction of the G-mesh
            !dfft%nsp(me) : Number of different z-pencils for current processor. Multiply this by NOGRP
            !nr3          : Length of each pencil
            !aux          : Output of results
            !sign         : Type of transform...+(forward)...-(inverse)
            !nr3x         : The length of the out vecs (stride between starts of vectors in the output)
            !------------------------------------------------------------------------------------------

            call cft_1z( XF, NOGRP*dfft%nsp(me), nr3, nr3x, sign, aux)

            !-----------------------------------------------------------------------------------------
            !ADDED COMMENTS: C. Bekas, Oct. 2005
            !-----------------------------------SUBROUTINE: FFT_SCATTER: Scatter data accros x-y planes
            !aux          : Input data (sequence of vectors) and output (overwritten)
            !nr3x         : Length of data across transformed z-direction (can be up to nr3+1)
            !dfft%nnr     : FFT data size?
            !f            : WORK space
            !dfft%nsp     : Control sizes of contigious slices along Z direction...
            !dfft%npp     : ...and mapping for the communication between processors
            !sign         : Type of transform...+(pencils to planes)...-(planes to pencils)
            !------------------------------------------------------------------------------------------

            call fft_scatter( aux, nr3x, NOGRP*dfft%nnr, f, NOGRP*dfft%nsp, NOGRP*dfft%npp, sign)

            !------------------------------------------------------------------------------------------
            !ADDED COMMENTS: C. Bekas, Oct. 2005
            !-----------------------------------Rearrange data in x-y planes:
            f(:) = (0.d0, 0.d0)
            do group_index = 1, NOGRP
               do i = 1, dfft%nst! FOR EACH ONE OF THE PENCILS
                  mc = dfft%ismap( i ) + (group_index-1)*dfft%nnr! THE POSITION OF THE PENCIL IN THE PLANE
                  do j = 1, dfft%npp(me)! FOR EACH ONE OF MY PLANES
                     f( mc + (j-1) * dfft%nnp ) = aux( j + (i-1) * nppx + (group_index-1)*dfft%nnr)
                  end do
               end do
            end do
            planes = dfft%iplp*NOGRP !THESE ARE THE Y PLANES TO BE FFTed
            !
            !
         else
            !
            !-----------------------------------
            !WAVE - FUNCTION FFT calculation
            !-----------------------------------

            local_send_cnt(1) = nr3x*dfft%nsw(mpime+1)
            local_send_displ(1) = 0
            local_recv_cnt(1) = nr3x*dfft%nsw(NOLIST(1)+1)
            local_recv_displ(1) = 0
            DO index=2, NOGRP
               local_send_cnt(index) = nr3x*dfft%nsw(mpime+1)
               local_send_displ(index) = local_send_displ(index-1) + strd ! local_send_cnt(index-1)

               local_recv_cnt(index) = nr3x*dfft%nsw(NOLIST(index)+1)
               local_recv_displ(index)  = local_recv_displ(index-1) + local_recv_cnt(index-1)
            ENDDO

            !------------------------------------------------------------------------------------------
            !ADDED COMMENTS: C. Bekas, Oct. 2005
            !-----------------------------------ANALOGOUS TO THE COMMENTS ABOVE

            CALL start_clock( 'ALLTOALL' )

#if defined __MPI
            CALL MPI_Alltoallv(f, local_send_cnt, local_send_displ, MPI_DOUBLE_COMPLEX, YF, local_recv_cnt, & 
         &                     local_recv_displ, MPI_DOUBLE_COMPLEX, ME_OGRP, IERR)  
#endif

            !-----------------------------------------
            !We need to get rid of all the zeros in XF
            !-----------------------------------------
            num_sticks = 0
            num_planes = 0
            DO ii=1, NOGRP
               num_sticks = num_sticks + dfft%nsw(NOLIST(ii)+1)
               num_planes = num_planes + dfft%npp(NOLIST(ii)+1)
            ENDDO

        !    IT1 = 1
        !    DO group_index = 1, NOGRP
               !----------------------------------------------------------
               !First find the global index of group_index
               !Then look for how many sticks there are
               !Then multiply by nr3: the number of coefficients per stick
               !----------------------------------------------------------
        !       HWMN = nr3 * ALL_Z_STICKS(NOLIST(group_index)+1)
        !       CALL DCOPY(2*HWMN, XF((group_index-1)*strd+1), 1, YF(IT1), 1)
        !       IT1 = IT1 + HWMN
        !    ENDDO


            CALL stop_clock( 'ALLTOALL' )
            !-------------------------------------------------------------
            !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
            !-------------------------------------------------------------
            !Do all decoupled FFTs across the Z-sticks
            !-------------------------------------------------------------
            CALL start_clock( '1D' )
            call cft_1z(YF, num_sticks, nr3, nr3x, sign, aux)
            CALL stop_clock( '1D' )

            !-------------------------------------------------------------------------------------
            !Transpose data for the 2-D FFT on the x-y plane
            !-----------------------------------------------
            !NOGRP*dfft%nnr: The length of aux and f
            !nr3x: The length of each Z-stick
            !aux: input - output
            !f: working space
            !sign: type of scatter
            !dfft%nsw(me) holds the number of Z-sticks proc. me has.
            !dfft%npp: number of planes per processor
            !-------------------------------------------------------------------------------------

            IF (tmp_npp(1).EQ.-1) THEN
#if defined __MPI
               CALL MPI_ALLGATHER(num_sticks, 1, MPI_INTEGER, tmp_nsw, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
               CALL MPI_ALLGATHER(num_planes, 1, MPI_INTEGER, tmp_npp, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
#endif
            ENDIF

            CALL start_clock( 'SCATTER' ) 
            call group_fft_scatter( aux, nr3x, (NOGRP+1)*strd, f, tmp_nsw, tmp_npp, sign)
            CALL stop_clock( 'SCATTER' )

            f(:) = (0.d0, 0.d0)

            ii=0
            do proc=1,nproc
                  do i=1,dfft%nsw(proc)
                     mc = dfft%ismap( i + dfft%iss(proc))
                     ii = ii + 1
                     do j=1,tmp_npp(me)
                        f(mc+(j-1)*nr1x*nr2x) = aux(j + (ii-1)*tmp_npp(me))
                     end do
                  end do
            end do

            planes = dfft%iplw

         end if

         !------------------------------------------------------------------------------------------
         !ADDED COMMENTS: C. Bekas, Oct. 2005
         !-----------------------------------DO THE 2-D FFT ON THE PLANES
         CALL start_clock( '2D' )
         call cft_2xy( f, tmp_npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
         CALL stop_clock( '2D' )

     !-----------------------
     !FORWARD FFT
     !-----------------------
      else
         if (sign.ne.-2) then
            planes = dfft%iplp
         else
            planes = dfft%iplw
         endif

!---------------------------------------------------------------------------
!        call cft_2xy( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
!---------------------------------------------------------------------------

         call cft_2xy( f, tmp_npp(me), nr1, nr2, nr1x, nr2x, sign, planes)


         !-----------------------------------
         !Density - Potential FFT calculation
         !-----------------------------------
         !NOT IMPLEMENTED YET
         !-----------------------------------
         if (sign.ne.-2) then
            do i=1,dfft%nst
               mc = dfft%ismap( i )
               do j=1,dfft%npp(me)
                  aux(j + (i-1)*nppx) = f(mc+(j-1)*dfft%nnp)
               end do
            end do
            call fft_scatter(aux,nr3x,dfft%nnr,f,dfft%nsp,dfft%npp,sign)
            call cft_1z( aux, dfft%nsp(me), nr3, nr3x, sign, f)

         !-----------------------------------
         !WAVE - FUNCTION FFT calculation
         !-----------------------------------
         else
            ii = 0
            do proc=1,nproc
               do i=1,dfft%nsw(proc)
                  mc = dfft%ismap( i + dfft%iss(proc) )
                  ii = ii + 1
                  do j=1,tmp_npp(me)
                     aux(j + (ii-1)*tmp_npp(me)) = f(mc+(j-1)*nr1x*nr2x)
                  end do
               end do
            end do

            call group_fft_scatter(aux,nr3x,(NOGRP+1)*strd,f,tmp_nsw,tmp_npp,sign)
            call cft_1z(aux,tmp_nsw(me),nr3,nr3x,sign,f)

         end if
      end if
!

      DEALLOCATE(XF)
      DEALLOCATE(YF)
      DEALLOCATE(aux)

      return
      !--------------
      !END
      !--------------
      end subroutine tg_cfft_cp
   !
   !
   !
END MODULE fft_cp
