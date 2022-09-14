!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  -------------------------------------------------
!  These subroutines were written by Carlo Cavazzoni 
!  -------------------------------------------------
!
!=----------------------------------------------------------------------=
MODULE smallbox_subs
!=----------------------------------------------------------------------=
   !! Subroutines generating G-vectors and variables needed to map
   !! G-vector components onto the FFT grid(s) in reciprocal space
   !! Small-Box grid.

   USE small_box,     ONLY :  bgb, tpibab
   USE smallbox_gvec, ONLY :  ngb, ngbl, gb, gxb, glb, npb, nmb, mill_b, gcutb
   USE fft_base,      ONLY : dfftb

   PRIVATE
   SAVE

   INTERFACE fft_oned2box
      MODULE PROCEDURE fft_oned2box_dp
   END INTERFACE
   INTERFACE fft_add_oned2box
      MODULE PROCEDURE fft_add_oned2box_dp
   END INTERFACE
   INTERFACE boxdotgrid
      MODULE PROCEDURE boxdotgrid_dp, boxdotgridcplx_dp
   END INTERFACE
   INTERFACE box2grid
      MODULE PROCEDURE box2grid_dp, box2grid2_dp
   END INTERFACE

   PUBLIC :: ggenb, gcalb, fft_oned2box, fft_add_oned2box
   PUBLIC :: boxdotgrid, box2grid

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
  !
  SUBROUTINE ggenb( ecutrho, iprsta )
    !-----------------------------------------------------------------------
    !! As \(\texttt{ggen}\), for the box grid. A "b" is appended to box variables.
    !! The documentation for \(\texttt{ggen}\) applies.
    !
    USE kinds, ONLY: DP
    USE io_global, ONLY: stdout, ionode
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(in) :: ecutrho
    INTEGER, INTENT (in) :: iprsta
    !
    INTEGER, ALLOCATABLE:: idx(:), iglb(:)
    INTEGER n1pb, n2pb, n3pb, n1mb, n2mb, n3mb
    INTEGER it, icurr, nr1m1, nr2m1, nr3m1, ir, ig, i,j,k, itv(3), ip
    REAL(DP) t(3), g2
    !
    !   gcutb is the effective cut-off for G-vectors of the small box grid
    !
    gcutb = ecutrho / tpibab**2
    !    
    nr1m1=dfftb%nr1-1
    nr2m1=dfftb%nr2-1
    nr3m1=dfftb%nr3-1
    ngb=0
    !
    !     first step : count the number of vectors with g2 < gcutb
    !
    !     exclude space with x<0
    !
    DO i= 0,nr1m1
       DO j=-nr2m1,nr2m1
          !
          !     exclude plane with x=0, y<0
          !
          IF(i==0.and.j<0) GOTO 10
          !
          DO k=-nr3m1,nr3m1
             !
             !     exclude line with x=0, y=0, z<0
             !
             IF(i==0.and.j==0.and.k<0) GOTO 20
             g2=0.d0
             DO ir=1,3
                t(ir) = dble(i)*bgb(ir,1)+dble(j)*bgb(ir,2)+dble(k)*bgb(ir,3)
                g2=g2+t(ir)*t(ir)
             ENDDO
             IF(g2>gcutb) GOTO 20
             ngb=ngb+1
20           CONTINUE
          ENDDO
10        CONTINUE
       ENDDO
    ENDDO
    !
    !     second step: allocate space
    !
    ALLOCATE(gxb(3,ngb))
    ALLOCATE(gb(ngb))
    ALLOCATE(npb(ngb))
    ALLOCATE(nmb(ngb))
    ALLOCATE(iglb(ngb))
    ALLOCATE(mill_b(3,ngb))
    ALLOCATE(idx(ngb))
    !
    !     third step : find the vectors with g2 < gcutb
    !
    ngb=0
    !
    !     exclude space with x<0
    !
    DO i= 0,nr1m1
       DO j=-nr2m1,nr2m1
          !
          !     exclude plane with x=0, y<0
          !
          IF(i==0.and.j<0) GOTO 15
          !
          DO k=-nr3m1,nr3m1
             !
             !     exclude line with x=0, y=0, z<0
             !
             IF(i==0.and.j==0.and.k<0) GOTO 25
             g2=0.d0
             DO ir=1,3
                t(ir) = dble(i)*bgb(ir,1)+dble(j)*bgb(ir,2)+dble(k)*bgb(ir,3)
                g2=g2+t(ir)*t(ir)
             ENDDO
             IF(g2>gcutb) GOTO 25
             ngb=ngb+1
             gb(ngb)=g2
             mill_b(1,ngb)=i
             mill_b(2,ngb)=j
             mill_b(3,ngb)=k
25           CONTINUE
          ENDDO
15        CONTINUE
       ENDDO
    ENDDO

    IF( iprsta > 3 ) THEN
       WRITE( stdout,*)
       WRITE( stdout,170) ngb
170    FORMAT(' ggenb: # of gb vectors < gcutb ngb = ',i6)
    ENDIF

    idx(1)=0
    CALL hpsort (ngb,gb,idx)

    DO ig=1,ngb-1
       icurr=ig
30     IF(idx(icurr)/=ig) THEN
          itv=mill_b(:,icurr)
          mill_b(:,icurr)=mill_b(:,idx(icurr))
          mill_b(:,idx(icurr))=itv

          it=icurr
          icurr=idx(icurr)
          idx(it)=it
          IF(idx(icurr)==ig) THEN
             idx(icurr)=icurr
             GOTO 35
          ENDIF
          GOTO 30
       ENDIF
35     CONTINUE
    ENDDO
    !
    DEALLOCATE(idx)
    !
    ! costruct fft indexes (n1b,n2b,n3b) for the box grid
    !
    DO ig=1,ngb
       i=mill_b(1,ig)
       j=mill_b(2,ig)
       k=mill_b(3,ig)
       n1pb=i+1
       n2pb=j+1
       n3pb=k+1
       !
       ! n1pb,n2pb,n3pb: indexes of G
       ! negative indexes are refolded (note that by construction i.ge.0)
       !
       IF(i<0) n1pb=n1pb+dfftb%nr1
       IF(j<0) n2pb=n2pb+dfftb%nr2
       IF(k<0) n3pb=n3pb+dfftb%nr3
       !
       ! n1mb,n2mb,n3mb: indexes of -G
       !
       IF(i==0) THEN
          n1mb=1
       ELSE
          n1mb=dfftb%nr1-n1pb+2
       ENDIF
       IF(j==0) THEN
          n2mb=1
       ELSE
          n2mb=dfftb%nr2-n2pb+2
       ENDIF
       IF(k==0) THEN
          n3mb=1
       ELSE
          n3mb=dfftb%nr3-n3pb+2
       ENDIF
       !
       ! conversion from (i,j,k) index to combined 1-d ijk index:
       ! ijk = 1 + (i-1)+(j-1)*ix+(k-1)*ix*jx
       ! where the (i,j,k) array is assumed to be dimensioned (ix,jx,kx)
       !
       npb(ig) = n1pb+(n2pb-1)*dfftb%nr1x+(n3pb-1)*dfftb%nr1x*dfftb%nr2x
       nmb(ig) = n1mb+(n2mb-1)*dfftb%nr1x+(n3mb-1)*dfftb%nr1x*dfftb%nr2x
    ENDDO
    !
    ! shells of G - first calculate their number and position
    !
    CALL gshcount( ngb, gb, ngbl, iglb ) 
    !
    IF( iprsta > 3 ) THEN
       WRITE( stdout,180) ngbl
180    FORMAT(' ggenb: # of gb shells  < gcutb ngbl= ',i6)
    ENDIF
    !
    ! then allocate the array glb
    !
    ALLOCATE(glb(ngbl))
    !
    ! and finally fill glb with the values of the shells
    !
    glb(iglb(1))=gb(1)
    DO ig=2,ngb
       IF(iglb(ig)/=iglb(ig-1)) glb(iglb(ig))=gb(ig)
    ENDDO
    !
    ! calculation of G-vectors
    !
    DO ig=1,ngb
       i=mill_b(1,ig)
       j=mill_b(2,ig)
       k=mill_b(3,ig)
       gxb(:,ig)=i*bgb(:,1)+j*bgb(:,2)+k*bgb(:,3)
    ENDDO
    !
    DEALLOCATE (iglb)
    !
    RETURN
    !
  END SUBROUTINE ggenb
  !
  !-------------------------------------------------------------------------
  SUBROUTINE gshcount( ng, gg, ngl, igl )
    !-------------------------------------------------------------------------
    !
    USE kinds, ONLY: DP
    !
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: ng
    REAL(DP),INTENT (IN) :: gg(ng)
    INTEGER, INTENT (OUT) :: ngl, igl(ng)

    INTEGER :: ig

    ngl=1
    igl(1)=ngl
    DO ig=2,ng
       IF(abs(gg(ig)-gg(ig-1))>1.e-6) THEN
          ngl=ngl+1
       ENDIF
       igl(ig)=ngl
    ENDDO

    RETURN
    !
  END SUBROUTINE gshcount
  !
  !
  SUBROUTINE gcalb( )
    !! Re-generation of little box g-vectors
    !
    USE kinds, ONLY: DP
    !
    IMPLICIT NONE
    !
    INTEGER :: ig, i1,i2,i3

    IF ( dfftb%nr1 == 0 .OR. dfftb%nr2 == 0 .OR. dfftb%nr3 == 0 ) return
    !
    do ig=1,ngb
       i1=mill_b(1,ig)
       i2=mill_b(2,ig)
       i3=mill_b(3,ig)
       gxb(:,ig)=i1*bgb(:,1)+i2*bgb(:,2)+i3*bgb(:,3)
       gb(ig)=gxb(1,ig)**2 + gxb(2,ig)**2 + gxb(3,ig)**2
    enddo
    !
    RETURN
  END SUBROUTINE gcalb
  !
  SUBROUTINE fft_oned2box_dp( qv, fg1, fg2 )
    USE kinds, ONLY: DP
    IMPLICIT NONE
    COMPLEX(DP), INTENT(OUT) :: qv(:)
    COMPLEX(DP), INTENT(IN) :: fg1(:)
    COMPLEX(DP), OPTIONAL, INTENT(IN) :: fg2(:)
    INTEGER :: ig
    COMPLEX(DP) :: ci
    ci = ( 0.0d0, 1.0d0 )
    !
    qv = ( 0.0d0, 0.0d0 )
    IF( PRESENT( fg2 ) ) THEN
       DO ig=1,ngb
          qv(npb(ig)) = fg1(ig) + ci * fg2(ig)
          qv(nmb(ig))=  CONJG(fg1(ig)) + ci * CONJG(fg2(ig))
       END DO
    ELSE
       DO ig=1,ngb
          qv(npb(ig)) =       fg1(ig)
          qv(nmb(ig)) = CONJG(fg1(ig))
       END DO
    END IF
    RETURN
  END SUBROUTINE fft_oned2box_dp

  SUBROUTINE fft_add_oned2box_dp( qv, fg1, fg2 )
    USE kinds, ONLY: DP
    IMPLICIT NONE
    COMPLEX(DP), INTENT(INOUT) :: qv(:)
    COMPLEX(DP), INTENT(IN) :: fg1(:)
    COMPLEX(DP), OPTIONAL, INTENT(IN) :: fg2(:)
    INTEGER :: ig
    COMPLEX(DP) :: ci
    ci = ( 0.0d0, 1.0d0 )
    !
    IF( PRESENT( fg2 ) ) THEN
       DO ig=1,ngb
          qv(npb(ig)) = qv(npb(ig)) + fg1(ig) + ci * fg2(ig)
          qv(nmb(ig)) = qv(nmb(ig)) + CONJG(fg1(ig)) + ci * CONJG(fg2(ig))
       END DO
    ELSE
       DO ig=1,ngb
          qv(npb(ig)) = qv(npb(ig)) +       fg1(ig)
          qv(nmb(ig)) = qv(nmb(ig)) + CONJG(fg1(ig))
       END DO
    END IF
    RETURN
  END SUBROUTINE fft_add_oned2box_dp


!-----------------------------------------------------------------------
      SUBROUTINE box2grid_dp( irb, nfft, qv, vr )
!-----------------------------------------------------------------------
      !! Add array \(\text{qv}(r)\) on box grid to array \(\text{vr}(r)\)
      !! on dense grid.
!
      USE kinds, ONLY: dp
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp
!
      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft
      !! nfft=1 : add real part of qv(r) to real part of array vr(r);  
      !! nfft=2 : add imaginary part of qv(r) to real part of array vr(r).
      INTEGER, INTENT(in):: irb(3)
      !! position of the box in the dense grid
      COMPLEX(dp), INTENT(in):: qv(dfftb%nnr)
      !! input array on box grid
      COMPLEX(dp), INTENT(inout):: vr(dfftp%nnr)
      !! array on dense grid
!
      ! ... local variables
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
!$omp critical
                  DO ir1=1,dfftb%nr1
                     ibig1=irb(1)+ir1-1
                     ibig1=1+MOD(ibig1-1,dfftp%nr1)
                     IF(ibig1.LT.1.OR.ibig1.GT.dfftp%nr1)                        &
     &                    CALL errore('box2grid','ibig1 wrong',ibig1)
                     ibig=ibig1+(ibig2-1)*dfftp%nr1x+(ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                     ir=ir1+(ir2-1)*dfftb%nr1x+(ir3-1)*dfftb%nr1x*dfftb%nr2x
                     IF( nfft == 1 ) THEN
                        vr(ibig) = vr(ibig)+REAL(qv(ir))
                     ELSE
                        vr(ibig) = vr(ibig)+AIMAG(qv(ir))
                     END IF
                  END DO
!$omp end critical
               END IF
            END DO
         END IF
      END DO
!
      RETURN
      END SUBROUTINE box2grid_dp


!-----------------------------------------------------------------------
      SUBROUTINE box2grid2_dp(irb,qv,v)
!-----------------------------------------------------------------------
      !! Add array \(\text{qv}(r)\) on box grid to array \(\text{v}(r)\)
      !! on dense grid.
!
      USE kinds, ONLY: dp
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in):: irb(3)
      !! position of the box in the dense grid
      COMPLEX(dp), INTENT(in):: qv(dfftb%nnr)
      !! input array on box grid
      COMPLEX(dp), INTENT(inout):: v(dfftp%nnr)
      !! array on dense grid
!
      ! ... local variables
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
!$omp critical
                  DO ir1=1,dfftb%nr1
                     ibig1=irb(1)+ir1-1
                     ibig1=1+MOD(ibig1-1,dfftp%nr1)
                     IF(ibig1.LT.1.OR.ibig1.GT.dfftp%nr1)                        &
     &                    CALL errore('box2grid2','ibig1 wrong',ibig1)
                     ibig=ibig1+(ibig2-1)*dfftp%nr1x+(ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                     ir=ir1+(ir2-1)*dfftb%nr1x+(ir3-1)*dfftb%nr1x*dfftb%nr2x
                     v(ibig) = v(ibig)+qv(ir)
                  END DO
!$omp end critical
               END IF
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE box2grid2_dp


!-----------------------------------------------------------------------
      REAL(8) FUNCTION boxdotgrid_dp( irb, nfft, qv, vr )
!-----------------------------------------------------------------------
      !! Calculates \( \sum_i \text{qv}(r_i)\cdot\text{vr}(r_i) \) with \(r_i\)
      !! on box grid, array \(\text{qv}(r)\) is defined on box grid, array
      !! \(\text{vr}(r)\) on dense grid.  
      !! Parallel execution: remember to sum the contributions from other
      !! nodes.
!
      USE kinds, ONLY: dp
      USE fft_base, ONLY: dfftp, dfftb
      USE mp_global, ONLY: me_bgrp
      IMPLICIT NONE
      INTEGER, INTENT(in):: nfft
      !! nfft=1 (2): use real (imaginary) part of qv(r)
      INTEGER, INTENT(in):: irb(3)
      !! position of the box in the dense grid
      COMPLEX(dp), INTENT(in):: qv(dfftb%nnr)
      !! box grid array
      REAL(dp), INTENT(in):: vr(dfftp%nnr)
      !! dense grid array
!
      ! ... local variables
!
      INTEGER ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
      INTEGER me
!
!
      IF(nfft.LE.0.OR.nfft.GT.2) CALL errore('boxdotgrid','wrong data',nfft)

      me = me_bgrp + 1

      boxdotgrid_dp=0.d0

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
                     IF( nfft == 1 ) THEN
                        boxdotgrid_dp = boxdotgrid_dp + REAL(qv(ir))*vr(ibig)
                     ELSE
                        boxdotgrid_dp = boxdotgrid_dp + AIMAG(qv(ir))*vr(ibig)
                     END IF
                  END DO
               ENDIF
            END DO
         ENDIF
      END DO

      RETURN
      END FUNCTION boxdotgrid_dp


!-----------------------------------------------------------------------
FUNCTION boxdotgridcplx_dp(irb,qv,vr)
  !-----------------------------------------------------------------------
  !! Calculate \(\sum_i \text{qv}(r_i)\cdot\text{vr}(r_i)\)  with \(r_i\)
  !! on box grid. Array \(\text{qv}(r)\) is defined on box grid, array 
  !! \(\text{vr}(r)\) on dense grid.
  !
  ! irb   : position of the box in the dense grid
  ! Parallel execution: remember to sum the contributions from other nodes
  !
  !      use ion_parameters
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp, dfftb
  USE mp_global,       ONLY : me_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER,           INTENT(IN):: irb(3)
  COMPLEX(DP), INTENT(IN):: qv(dfftb%nnr), vr(dfftp%nnr)
  COMPLEX(DP)            :: boxdotgridcplx_dp
  !
  INTEGER :: ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig, me
  !
  me = me_bgrp + 1
  !
  boxdotgridcplx_dp = 0.0_DP

  DO ir3=1,dfftb%nr3
     ibig3=irb(3)+ir3-1
     ibig3=1+MOD(ibig3-1,dfftp%nr3)
#if defined(__MPI)
     ibig3 = ibig3 - dfftp%my_i0r3p
     IF (ibig3.GT.0.AND.ibig3.LE.dfftp%my_nr3p) THEN
#endif
        DO ir2=1,dfftb%nr2
           ibig2=irb(2)+ir2-1
           ibig2=1+MOD(ibig2-1,dfftp%nr2)
#if defined(__MPI)
           ibig2 = ibig2 - dfftp%my_i0r2p
           IF (ibig2.GT.0.AND.ibig2.LE.dfftp%my_nr2p) THEN
#endif
              DO ir1=1,dfftb%nr1
                 ibig1=irb(1)+ir1-1
                 ibig1=1+MOD(ibig1-1,dfftp%nr1)
                 ibig=ibig1 + (ibig2-1)*dfftp%nr1x + (ibig3-1)*dfftp%nr1x*dfftp%my_nr2p
                 ir  =ir1 + (ir2-1)*dfftb%nr1x + (ir3-1)*dfftb%nr1x*dfftb%nr2x
                 boxdotgridcplx_dp = boxdotgridcplx_dp + qv(ir)*vr(ibig)
              END DO
#if defined(__MPI)
           ENDIF
#endif
        END DO
#if defined(__MPI)
     ENDIF
#endif
  END DO
  !
  RETURN
  !
END FUNCTION boxdotgridcplx_dp
!
!=----------------------------------------------------------------------=
END MODULE smallbox_subs
!=----------------------------------------------------------------------=
