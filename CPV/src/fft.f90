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
      SUBROUTINE box2grid(irb,nfft,qv,vr)
!-----------------------------------------------------------------------
!
! add array qv(r) on box grid to array vr(r) on dense grid
! irb   : position of the box in the dense grid
! nfft=1  add      real part of qv(r) to real part of array vr(r)
! nfft=2  add imaginary part of qv(r) to real part of array vr(r)
!
      USE kinds, ONLY: dp
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp

      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(dp), INTENT(in):: qv(2,dfftb%nnr)
      COMPLEX(dp), INTENT(inout):: vr(dfftp%nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me

      IF(nfft.LE.0.OR.nfft.GT.2) CALL errore('box2grid','wrong data',nfft)

      me = me_bgrp + 1

      DO ir3=1,dfftb%nr3
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,dfftp%nr3)
         IF(ibig3.LT.1.OR.ibig3.GT.dfftp%nr3)                                 &
     &        CALL errore('box2grid','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%my_i0r3p
         IF ( ibig3 .GT. 0 .AND. ibig3 .LE. ( dfftp%my_nr3p ) ) THEN

            DO ir2=1,dfftb%nr2
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,dfftp%nr2)
               IF(ibig2.LT.1.OR.ibig2.GT.dfftp%nr2)                           &
     &              CALL errore('box2grid','ibig2 wrong',ibig2)
               ibig2=ibig2-dfftp%my_i0r2p
               IF ( ibig2 .GT. 0 .AND. ibig2 .LE. ( dfftp%my_nr2p ) ) THEN

                  DO ir1=1,dfftb%nr1
                     ibig1=irb(1)+ir1-1
                     ibig1=1+MOD(ibig1-1,dfftp%nr1)
                     IF(ibig1.LT.1.OR.ibig1.GT.dfftp%nr1)                        &
     &                    CALL errore('box2grid','ibig1 wrong',ibig1)
                     ibig=ibig1+(ibig2-1)*dfftp%nr1x+(ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                     ir=ir1+(ir2-1)*dfftb%nr1x+(ir3-1)*dfftb%nr1x*dfftb%nr2x
!$omp critical
                     vr(ibig) = vr(ibig)+qv(nfft,ir)
!$omp end critical
                  END DO
               END IF
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
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in):: irb(3)
      COMPLEX(dp), INTENT(in):: qv(dfftb%nnr)
      COMPLEX(dp), INTENT(inout):: v(dfftp%nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me

      me = me_bgrp + 1

      DO ir3=1,dfftb%nr3
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,dfftp%nr3)
         IF(ibig3.LT.1.OR.ibig3.GT.dfftp%nr3)                                 &
     &        CALL errore('box2grid2','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%my_i0r3p
         IF (ibig3.GT.0.AND.ibig3.LE. dfftp%my_nr3p ) THEN

            DO ir2=1,dfftb%nr2
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,dfftp%nr2)
               IF(ibig2.LT.1.OR.ibig2.GT.dfftp%nr2)                           &
     &              CALL errore('box2grid2','ibig2 wrong',ibig2)
               ibig2=ibig2-dfftp%my_i0r2p
               IF (ibig2.GT.0.AND.ibig2.LE. dfftp%my_nr2p ) THEN

                  DO ir1=1,dfftb%nr1
                     ibig1=irb(1)+ir1-1
                     ibig1=1+MOD(ibig1-1,dfftp%nr1)
                     IF(ibig1.LT.1.OR.ibig1.GT.dfftp%nr1)                        &
     &                    CALL errore('box2grid2','ibig1 wrong',ibig1)
                     ibig=ibig1+(ibig2-1)*dfftp%nr1x+(ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                     ir=ir1+(ir2-1)*dfftb%nr1x+(ir3-1)*dfftb%nr1x*dfftb%nr2x
                     v(ibig) = v(ibig)+qv(ir)
                  END DO
               END IF
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
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp
      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft, irb(3)
      REAL(dp), INTENT(in):: qv(2,dfftb%nnr), vr(dfftp%nnr)
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me
!
!
      IF(nfft.LE.0.OR.nfft.GT.2) CALL errore('boxdotgrid','wrong data',nfft)

      me = me_bgrp + 1

      boxdotgrid=0.d0

      DO ir3=1,dfftb%nr3
         ibig3=irb(3)+ir3-1
         ibig3=1+MOD(ibig3-1,dfftp%nr3)
         ibig3=ibig3-dfftp%my_i0r3p
         IF (ibig3.GT.0.AND.ibig3.LE. dfftp%my_nr3p ) THEN
            DO ir2=1,dfftb%nr2
               ibig2=irb(2)+ir2-1
               ibig2=1+MOD(ibig2-1,dfftp%nr2)
               ibig2=ibig2-dfftp%my_i0r2p
               IF (ibig2.GT.0.AND.ibig2.LE. dfftp%my_nr2p ) THEN
                  DO ir1=1,dfftb%nr1
                     ibig1=irb(1)+ir1-1
                     ibig1=1+MOD(ibig1-1,dfftp%nr1)
                     ibig=ibig1 + (ibig2-1)*dfftp%nr1x + (ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                     ir  =ir1 + (ir2-1)*dfftb%nr1x + (ir3-1)*dfftb%nr1x*dfftb%nr2x
                     boxdotgrid = boxdotgrid + qv(nfft,ir)*vr(ibig)
                  END DO
               ENDIF
            END DO
         ENDIF
      END DO

      RETURN
      END FUNCTION boxdotgrid

