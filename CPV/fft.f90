!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  ----------------------------------------------
!  These subroutines written by Carlo Cavazzoni 
!  ----------------------------------------------

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
      USE kinds, ONLY: dp
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nrxx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image

      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(dp), INTENT(in):: qv(2,nnrb)
      COMPLEX(dp), INTENT(inout):: vr(nnr)
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
!$omp critical
                  vr(ibig) = vr(ibig)+qv(nfft,ir)
!$omp end critical
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
      USE kinds, ONLY: dp
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nrxx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in):: irb(3)
      COMPLEX(dp), INTENT(in):: qv(nnrb)
      COMPLEX(dp), INTENT(inout):: v(nnr)
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
      USE kinds, ONLY: dp
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nrxx
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      USE fft_base, ONLY: dfftp
      USE mp_global, ONLY: me_image
      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(dp), INTENT(in):: qv(2,nnrb), vr(nnr)
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

