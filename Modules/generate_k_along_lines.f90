SUBROUTINE generate_k_along_lines(nkaux, xkaux, wkaux, xk, wk, nkstot)
!
!  This routine recieves as input a set of k point (xkaux) and integer weights
!  (wkaux) and generates a set of k points along the lines 
!  xkaux(:,i+1)-xkaux(:,i). Each line contains wkaux(i) points.
!  The weights of each k point wk(i) is the length of the path from xk(:,1)
!  to xk(i). Points with wkaux=0 do not increase the path length.
!  The total number of output points must be nkstot, and xk and wk must
!  be array of length nkstot.
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN) :: nkaux, nkstot, wkaux(nkaux)
REAL(DP), INTENT(IN) :: xkaux(3,nkaux)
REAL(DP), INTENT(OUT) :: xk(3,nkstot), wk(nkstot)

INTEGER :: nkstot_, i, j
REAL(DP) :: delta, xkmod

nkstot_=1
wk(1)=0.0_DP
xk(:,1)=xkaux(:,1)
DO i=2,nkaux
   IF (wkaux(i-1)>0) THEN
      delta=1.0_DP/wkaux(i-1)
      DO j=1,wkaux(i-1)
         nkstot_=nkstot_+1
         IF (nkstot_ > nkstot) CALL errore ('generate_k_along_lines', &
                                            'internal error 1: wrong nkstot',i)
         xk(:,nkstot_)=xkaux(:,i-1)+delta*j*(xkaux(:,i)-xkaux(:,i-1))
         xkmod=SQRT( (xk(1,nkstot_)-xk(1,nkstot_-1))**2 +   &
                     (xk(2,nkstot_)-xk(2,nkstot_-1))**2 +   &
                     (xk(3,nkstot_)-xk(3,nkstot_-1))**2 )
         wk(nkstot_)=wk(nkstot_-1) + xkmod
      ENDDO
   ELSEIF (wkaux(i-1)==0) THEN
      nkstot_=nkstot_+1
      IF (nkstot_ > nkstot) CALL errore ('generate_k_along_lines', &
                                           'internal error 2: wrong nkstot',i)
      IF (nkstot_ ==1 ) CALL errore ('generate_k_along_lines', &
                                            'problems with weights',i)
      xk(:,nkstot_)=xkaux(:,i)
      wk(nkstot_)=wk(nkstot_-1) 
   ELSE
      CALL errore ('generate_k_along_lines', 'wrong number of points',i)
   ENDIF
ENDDO
IF (nkstot_ /= nkstot) CALL errore ('generate_k_along_lines', &
                                    'internal error 3: wrong nkstot',nkstot_)

RETURN
END SUBROUTINE generate_k_along_lines

SUBROUTINE generate_k_in_plane(nkaux, xkaux, wkaux, xk, wk, nkstot)
!
!   Generate a uniform mesh of k points on the plane defined by
!   the origin xkaux(:,1), and two vectors xkaux(:,2) and xkaux(:,3).
!   The size of the mesh is wkaux(2)*wkaux(3).
!
 
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN) :: nkaux, nkstot, wkaux(nkaux)
REAL(DP), INTENT(IN) :: xkaux(3,nkaux)
REAL(DP), INTENT(OUT) :: xk(3,nkstot), wk(nkstot)

REAL(DP) :: dkx(3), dky(3), wk0
INTEGER :: ijk, i, j

dkx(:)=(xkaux(:,2)-xkaux(:,1))/(wkaux(2)-1.0_DP)
dky(:)=(xkaux(:,3)-xkaux(:,1))/(wkaux(3)-1.0_DP)
wk0=1.0_DP/nkstot
ijk=0
DO i=1, wkaux(2)
   DO j = 1, wkaux(3)
      ijk=ijk+1
      IF (ijk > nkstot) CALL errore ('generate_k_in_plane', &
                                            'internal error : wrong nstot',i)
      
      xk(:,ijk) = xkaux(:,1) + dkx(:)*(i-1) + dky(:) * (j-1)
      wk(ijk) = wk0
   ENDDO
ENDDO

RETURN

END SUBROUTINE generate_k_in_plane
