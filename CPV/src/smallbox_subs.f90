!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------=
MODULE smallbox_subs
!=----------------------------------------------------------------------=

!  ... subroutines generating G-vectors and variables needed to map
!  ... G-vector components onto the FFT grid(s) in reciprocal space
!  ... Small-Box grid

   USE small_box,     ONLY :  bgb, tpibab
   USE smallbox_gvec, ONLY :  ngb, ngbl, gb, gxb, glb, npb, nmb, mill_b, gcutb
   USE fft_base,      ONLY : dfftb

   PRIVATE
   SAVE

   PUBLIC :: ggenb, gcalb

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
  !
  SUBROUTINE ggenb ( ecutrho, iprsta )
    !-----------------------------------------------------------------------
    !
    ! As ggen, for the box grid. A "b" is appended to box variables.
    ! The documentation for ggen applies
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
  SUBROUTINE gcalb ( )
    !
    !     re-generation of little box g-vectors
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
!=----------------------------------------------------------------------=
END MODULE smallbox_subs
!=----------------------------------------------------------------------=
