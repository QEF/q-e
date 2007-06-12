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


!-----------------------------------------------------------------------
  subroutine invfft_x( grid_type, f, nr1, nr2, nr3, nr1x, nr2x, nr3x, ia )
!-----------------------------------------------------------------------
! grid_type = 'Dense'
!   inverse fourier transform of potentials and charge density
!   on the dense grid . On output, f is overwritten
! grid_type = 'Smooth'
!   inverse fourier transform of  potentials and charge density
!   on the smooth grid . On output, f is overwritten
! grid_type = 'Wave'
!   inverse fourier transform of  wave functions
!   on the smooth grid . On output, f is overwritten
! grid_type = 'Box'
!   not-so-parallel 3d fft for box grid, implemented only for sign=1
!   G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!   The array f (overwritten on output) is NOT distributed:
!   a copy is present on each processor.
!   The fft along z  is done on the entire grid.
!   The fft along xy is done only on planes that have components on the
!   dense grid for each processor. Note that the final array will no
!   longer be the same on all processors.
!

!
      USE kinds,         ONLY: DP
      use fft_cp,        only: cfft_cp, tg_cfft_cp
      use fft_base,      only: dfftp, dffts, dfftb
      use fft_scalar,    only: cfft3d, cfft3ds, cft_b
      USE control_flags, ONLY: use_task_groups

      IMPLICIT none

      INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
      INTEGER, OPTIONAL, INTENT(IN) :: ia
      CHARACTER(LEN=*), INTENT(IN) :: grid_type
      COMPLEX(DP) :: f(:)
      !
      INTEGER :: imin3, imax3, np3

      IF( grid_type == 'Dense' ) THEN
         call start_clock( 'fft' )
      ELSE IF( grid_type == 'Smooth' ) THEN
         call start_clock( 'ffts' )
      ELSE IF( grid_type == 'Wave' ) THEN
         call start_clock('fftw')
      ELSE IF( grid_type == 'Box' ) THEN
         call start_clock( 'fftb' )
      ELSE 
         call errore( ' invfft ', ' unknown grid: '//grid_type , 1 )
      END IF

#if defined __PARA && !defined __USE_3D_FFT

      IF( grid_type == 'Box' ) THEN
         imin3 = dfftb%imin3( ia )
         imax3 = dfftb%imax3( ia )
         np3   = dfftb%np3( ia )   ! imax3 - imin3 + 1
      END IF
      
      IF( grid_type == 'Dense' ) THEN
         call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1,dfftp)
      ELSE IF( grid_type == 'Smooth' ) THEN
         call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1,dffts)
      ELSE IF( grid_type == 'Wave' ) THEN
         IF( use_task_groups ) THEN
            call tg_cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,2,dffts)
         ELSE
            call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,2,dffts)
         END IF
      ELSE IF( grid_type == 'Box' .AND. np3 > 0 ) THEN
         call cft_b(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,imin3,imax3,1)
      END IF

#else


# if defined __COMPLIB || __SCSL || __SX6 || __USE_3D_FFT

      call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1)

# elif defined __ESSL || __LINUX_ESSL || __FFTW  || __FFTW3

      IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
          grid_type == 'Box' ) THEN
         call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1)
      ELSE IF( grid_type == 'Wave' ) THEN
         call cfft3ds(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1,dffts%isind, dffts%iplw)
      END IF
# endif

#endif

      IF( grid_type == 'Dense' ) THEN
         call stop_clock( 'fft' )
      ELSE IF( grid_type == 'Smooth' ) THEN
         call stop_clock( 'ffts' )
      ELSE IF( grid_type == 'Wave' ) THEN
         call stop_clock('fftw')
      ELSE IF( grid_type == 'Box' ) THEN
         call stop_clock( 'fftb' )
      END IF
!
      return
      end subroutine invfft_x



!-----------------------------------------------------------------------
      subroutine fwfft_x(grid_type,f,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!-----------------------------------------------------------------------
! grid_type = 'Dense'
!   forward fourier transform of potentials and charge density 
!   on the dense grid . On output, f is overwritten
! grid_type = 'Smooth'
!   forward fourier transform of potentials and charge density
!   on the smooth grid . On output, f is overwritten
! grid_type = 'Wave'
!   forward fourier transform of  wave functions
!   on the smooth grid . On output, f is overwritten
! 
      USE kinds,         ONLY: DP
      use fft_cp,        only: cfft_cp, tg_cfft_cp
      use fft_base,      only: dfftp, dffts
      use fft_scalar,    only: cfft3d, cfft3ds
      USE control_flags, ONLY: use_task_groups

      implicit none

      INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
      CHARACTER(LEN=*), INTENT(IN) :: grid_type
      COMPLEX(DP) :: f(:)

      IF( grid_type == 'Dense' ) THEN
         call start_clock( 'fft' )
      ELSE IF( grid_type == 'Smooth' ) THEN
         call start_clock( 'ffts' )
      ELSE IF( grid_type == 'Wave' ) THEN
         call start_clock( 'fftw' )
      ELSE
         call errore( ' fwfft ', ' unknown grid: '//grid_type , 1 )
      END IF

#if defined __PARA && !defined __USE_3D_FFT

      IF( grid_type == 'Dense' ) THEN
         call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1,dfftp)
      ELSE IF( grid_type == 'Smooth' ) THEN
         call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1,dffts)
      ELSE IF( grid_type == 'Wave' ) THEN
         IF( use_task_groups ) THEN
            call tg_cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-2,dffts)
         ELSE
            call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-2,dffts)
         END IF
      END IF

#else 

# if defined __COMPLIB || __SCSL || __SX6 || __USE_3D_FFT

      call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1)

# elif defined __ESSL || __LINUX_ESSL || __FFTW  || __FFTW3

      IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' ) THEN
         call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1)
      ELSE IF( grid_type == 'Wave' ) THEN
         call cfft3ds(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1,dffts%isind, dffts%iplw)
      END IF

# endif

#endif

      IF( grid_type == 'Dense' ) THEN
         call stop_clock( 'fft' )
      ELSE IF( grid_type == 'Smooth' ) THEN
         call stop_clock( 'ffts' )
      ELSE IF( grid_type == 'Wave' ) THEN
         call stop_clock( 'fftw' )
      END IF

      return
      end subroutine fwfft_x



!-----------------------------------------------------------------------


    SUBROUTINE c2psi( psi, nnr, c, ca, ng, iflg )
       !
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), c(*), ca(*)
       integer, intent(in) :: nnr, ng, iflg

       complex(DP), parameter :: ci=(0.0d0,1.0d0)
       integer :: ig
       
         psi( 1 : nnr ) = 0.0d0

         !
         !  iflg "cases"
         !
         !  0   Do not use gamma symmetry
         !
         !  1   set psi using a wf with Gamma symmetry
         
         !  2   set psi combining two wf with Gamma symmetry
         !

         SELECT CASE ( iflg )
           !
           !  Case 0, 1 and 2  SMOOTH MESH
           !
           CASE ( 0 )
             !
             do ig = 1, ng
               psi( nps( ig ) ) = c( ig )
             end do
             !
           CASE ( 1 )
             !
             do ig = 1, ng
               psi( nms( ig ) ) = CONJG( c( ig ) )
               psi( nps( ig ) ) = c( ig )
             end do
             !
           CASE ( 2 )
             !
             do ig = 1, ng
               psi( nms( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ) )
               psi( nps( ig ) ) = c( ig ) + ci * ca( ig )
             end do

           CASE DEFAULT
             !
             CALL errore(" c2psi "," wrong value for iflg ", ABS( iflg ) )

         END SELECT

        return
     END SUBROUTINE c2psi

!
!
!

     SUBROUTINE rho2psi( grid_type, psi, nnr, rho, ng )
       !
       use recvecs_indexes, only: nm, np
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), rho(*)
       integer, intent(in) :: nnr, ng
       character(len=*), intent(in) :: grid_type

       integer :: ig
       
         psi( 1 : nnr ) = 0.0d0

         SELECT CASE ( grid_type )
           !
           !  Case 0, 1 and 2  SMOOTH MESH
           !
           CASE ( 'Smooth' )
             !
             ! without gamma sym
             ! do ig = 1, ng
             !   psi( nps( ig ) ) = rho( ig )
             ! end do
             !
             do ig = 1, ng
               psi( nms( ig ) ) = CONJG( rho( ig ) )
               psi( nps( ig ) ) = rho( ig )
             end do
             !
           CASE ( 'Dense' )
             !
             ! do ig = 1, ng
             !   psi( np( ig ) ) = rho( ig )
             ! end do
             !
             do ig = 1, ng
               psi( nm( ig ) ) = CONJG( rho( ig ) )
               psi( np( ig ) ) = rho( ig )
             end do
             !
           CASE DEFAULT
             !
             CALL errore(" rho2psi "," wrong grid "//grid_type , 1 )

         END SELECT

       return
     END SUBROUTINE rho2psi

!-----------------------------------------------------------------------

     SUBROUTINE psi2c( psi, nnr, c, ca, ng, iflg )

       use recvecs_indexes, only: nm, np
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), c(*), ca(*)
       integer, intent(in) :: nnr, ng, iflg

       complex(DP), parameter :: ci=(0.0d0,1.0d0)
       integer :: ig

         !
         !  iflg "cases"
         !
         !  0, 10   Do not use gamma symmetry
         !
         !  1, 11   set psi using a wf with Gamma symmetry
         !
         !  2, 12   set psi combining two wf with Gamma symmetry
         !
       
         SELECT CASE ( iflg )

           !
           !  Case 0, 1 and 2  SMOOTH MESH
           !
           CASE ( 0 )
             !
             do ig = 1, ng
               c( ig ) = psi( nps( ig ) )
             end do
             !
           CASE ( 1 )
             !
             CALL errore(" psi2c "," wrong value for iflg ", 11 )
             !
           CASE ( 2 )
             !
             DO ig = 1, ng
               ca(ig) = psi( nms( ig ) )
               c (ig) = psi( nps( ig ) )
             END DO

           !
           !  Case 10, 11 and 12  DENSE MESH
           !
           CASE ( 10 )
             !
             do ig = 1, ng
               c( ig ) = psi( np( ig ) )
             end do
             !
           CASE ( 11 )
             !
             CALL errore(" psi2c "," wrong value for iflg ", 1 )
             !
           CASE ( 12 )
             !
             DO ig = 1, ng
               ca(ig) = psi( nm( ig ) )
               c (ig) = psi( np( ig ) )
             END DO

           CASE DEFAULT
             !
             CALL errore(" psi2c "," wrong value for iflg ", ABS( iflg ) )

         END SELECT
       
       return
     END SUBROUTINE psi2c

!-----------------------------------------------------------------------

     SUBROUTINE psi2rho( grid_type, psi, nnr, rho, ng )

       use recvecs_indexes, only: nm, np
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), rho(*)
       integer, intent(in) :: nnr, ng
       character(len=*), intent(in) :: grid_type

       integer :: ig

       SELECT CASE ( grid_type )
          !
          CASE ( 'Smooth' )
             !
             do ig = 1, ng
                rho( ig ) = psi( nps( ig ) )
             end do
             !
          CASE ( 'Dense' )
             !
             do ig = 1, ng
               rho( ig ) = psi( np( ig ) )
             end do
             !
          CASE DEFAULT
             !
             CALL errore(" psi2rho "," wrong grid "//grid_type , 1 )

         END SELECT
       
       return
     END SUBROUTINE psi2rho



!-----------------------------------------------------------------------
      SUBROUTINE box2grid(irb,nfft,qv,vr)
!-----------------------------------------------------------------------
!
! add array qv(r) on box grid to array vr(r) on dense grid
! irb   : position of the box in the dense grid
! nfft=1  add      real part of qv(r) to real part of array vr(r)
! nfft=2  add imaginary part of qv(r) to real part of array vr(r)
!
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image

      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(8), INTENT(in):: qv(2,nnrb)
      COMPLEX(8), INTENT(inout):: vr(nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me

      IF(nfft.LE.0.OR.nfft.GT.2) CALL errore('box2grid','wrong data',nfft)

      me = me_image + 1

      DO ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,nr3)
         IF(ibig3.LT.1.OR.ibig3.GT.nr3)                                 &
     &        CALL errore('box2grid','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         IF ( ibig3 .GT. 0 .AND. ibig3 .LE. ( dfftp%npp(me) ) ) THEN
            DO ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,nr2)
               IF(ibig2.LT.1.OR.ibig2.GT.nr2)                           &
     &              CALL errore('box2grid','ibig2 wrong',ibig2)
               DO ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+MOD(ibig1-1,nr1)
                  IF(ibig1.LT.1.OR.ibig1.GT.nr1)                        &
     &                 CALL errore('box2grid','ibig1 wrong',ibig1)
                  ibig=ibig1+(ibig2-1)*nr1x+(ibig3-1)*nr1x*nr2x
                  ir=ir1+(ir2-1)*nr1bx+(ir3-1)*nr1bx*nr2bx
                  vr(ibig) = vr(ibig)+qv(nfft,ir)
               END DO
            END DO
         END IF
      END DO
!
      RETURN
      END SUBROUTINE box2grid


!-----------------------------------------------------------------------
      SUBROUTINE box2grid2(irb,qv,v)
!-----------------------------------------------------------------------
!
! add array qv(r) on box grid to array v(r) on dense grid
! irb   : position of the box in the dense grid
!
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in):: irb(3)
      COMPLEX(8), INTENT(in):: qv(nnrb)
      COMPLEX(8), INTENT(inout):: v(nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me

      me = me_image + 1

      DO ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,nr3)
         IF(ibig3.LT.1.OR.ibig3.GT.nr3)                                 &
     &        CALL errore('box2grid2','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         IF (ibig3.GT.0.AND.ibig3.LE. dfftp%npp(me) ) THEN
            DO ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,nr2)
               IF(ibig2.LT.1.OR.ibig2.GT.nr2)                           &
     &              CALL errore('box2grid2','ibig2 wrong',ibig2)
               DO ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+MOD(ibig1-1,nr1)
                  IF(ibig1.LT.1.OR.ibig1.GT.nr1)                        &
     &                 CALL errore('box2grid2','ibig1 wrong',ibig1)
                  ibig=ibig1+(ibig2-1)*nr1x+(ibig3-1)*nr1x*nr2x
                  ir=ir1+(ir2-1)*nr1bx+(ir3-1)*nr1bx*nr2bx
                  v(ibig) = v(ibig)+qv(ir)
               END DO
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE box2grid2


!-----------------------------------------------------------------------
      REAL(8) FUNCTION boxdotgrid(irb,nfft,qv,vr)
!-----------------------------------------------------------------------
!
! Calculate \sum_i qv(r_i)*vr(r_i)  with r_i on box grid
! array qv(r) is defined on box grid, array vr(r)on dense grid
! irb   : position of the box in the dense grid
! nfft=1 (2): use real (imaginary) part of qv(r)
! Parallel execution: remember to sum the contributions from other nodes
!
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image
      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(8), INTENT(in):: qv(2,nnrb), vr(nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me
!
!
      IF(nfft.LE.0.OR.nfft.GT.2) CALL errore('boxdotgrid','wrong data',nfft)

      me = me_image + 1

      boxdotgrid=0.d0

      DO ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,nr3)
         ibig3=ibig3-dfftp%ipp(me)
         IF (ibig3.GT.0.AND.ibig3.LE. dfftp%npp(me) ) THEN
            DO ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,nr2)
               DO ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+MOD(ibig1-1,nr1)
                  ibig=ibig1 + (ibig2-1)*nr1x + (ibig3-1)*nr1x*nr2x
                  ir  =ir1 + (ir2-1)*nr1bx + (ir3-1)*nr1bx*nr2bx
                  boxdotgrid = boxdotgrid + qv(nfft,ir)*vr(ibig)
               END DO
            END DO
         ENDIF
      END DO

      RETURN
      END FUNCTION boxdotgrid


!
!----------------------------------------------------------------------
      subroutine parabox(nr3b,irb3,nr3,imin3,imax3)
!----------------------------------------------------------------------
!
! find if box grid planes in the z direction have component on the dense
! grid on this processor, and if, which range imin3-imax3
!
      use mp_global, only: me_image
      use fft_base, only: dfftp
! input
      integer nr3b,irb3,nr3
! output
      integer imin3,imax3
! local
      integer ir3, ibig3, me
!
      me = me_image + 1
      imin3=nr3b
      imax3=1
      do ir3=1,nr3b
         ibig3=1+mod(irb3+ir3-2,nr3)
         if(ibig3.lt.1.or.ibig3.gt.nr3)                                 &
     &        call errore('cfftpb','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         if (ibig3.gt.0.and.ibig3.le.dfftp%npp(me)) then
            imin3=min(imin3,ir3)
            imax3=max(imax3,ir3)
         end if
      end do
!
      return
      end subroutine parabox

