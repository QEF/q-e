!                                                                            
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid_epw ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                         bg, npk, nk1, nk2, nk3, nks, xk, wk, BZtoIBZ, s_BZtoIBZ)
!-----------------------------------------------------------------------
!!
!!  Automatic generation of a uniform grid of k-points with symmetry. 
!!  Routine copied from PW/src/kpoint_grid.f90.
!!  We had to duplicate because the BZtoIBZ array was deallocated and is needed in
!!  EPW 
!!
USE kinds, ONLY: DP
IMPLICIT NONE
!
INTEGER, INTENT(in) :: nrot
!! Number of Bravais symmetry
INTEGER, INTENT(in) :: npk, nk1, nk2, nk3, t_rev(48), s(3,3,48)
INTEGER, INTENT(inout) :: s_BZtoIBZ(3,3,nk1*nk2*nk3)
!! Symeetry matrix that links an point to its IBZ friend.
INTEGER, INTENT(inout) :: BZtoIBZ (nk1*nk2*nk3)
!! Number of rotation
LOGICAL, INTENT(in) :: time_reversal
!! True if time reversal
LOGICAL, INTENT(in) :: skip_equivalence
!! True if equivalent point

REAL(DP), INTENT(in) :: bg(3,3)
!! Reciprocal space vectors
!
INTEGER, INTENT(out) :: nks
real(DP), INTENT(out):: xk(3,npk)
real(DP), INTENT(out):: wk(npk)
! LOCAL:

INTEGER :: s_save(3,3,nk1*nk2*nk3)

real(DP), PARAMETER :: eps=1.0d-5
real(DP) :: xkr(3), fact, xx, yy, zz
real(DP), ALLOCATABLE:: xkg(:,:), wkk(:)
INTEGER :: nkr, i,j,k, ns, n, nk, equiv(nk1*nk2*nk3), ik
LOGICAL :: in_the_list
!
nkr=nk1*nk2*nk3
ALLOCATE (xkg( 3,nkr),wkk(nkr))
equiv(:) = 0
s_save(:,:,:) = 0
!
DO i=1,nk1
   DO j=1,nk2
      DO k=1,nk3
         !  this is nothing but consecutive ordering
         n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
         !  xkg are the components of the complete grid in crystal axis
         xkg(1,n) = dble(i-1)/nk1 
         xkg(2,n) = dble(j-1)/nk2 
         xkg(3,n) = dble(k-1)/nk3 
      ENDDO
   ENDDO
ENDDO
!  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
!  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
DO nk=1,nkr
   equiv(nk)=nk
ENDDO

IF ( skip_equivalence ) THEN
  CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
  wkk = 1.d0
ELSE
  DO nk=1,nkr
  !  check if this k-point has already been found equivalent to another
    IF (equiv(nk) == nk) THEN
      wkk(nk)   = 1.0d0
      !  check if there are equivalent k-point to this in the list
      !  (excepted those previously found to be equivalent to another)
      !  check both k and -k
      DO ns=1,nrot
         DO i=1,3
            xkr(i) = s(i,1,ns) * xkg(1,nk) &
                   + s(i,2,ns) * xkg(2,nk) &
                   + s(i,3,ns) * xkg(3,nk)
            xkr(i) = xkr(i) - nint( xkr(i) )
         ENDDO
         IF(t_rev(ns)==1) xkr = -xkr
         xx = xkr(1)*nk1 
         yy = xkr(2)*nk2 
         zz = xkr(3)*nk3 
         in_the_list = abs(xx-nint(xx))<=eps .and. &
                       abs(yy-nint(yy))<=eps .and. &
                       abs(zz-nint(zz))<=eps
         IF (in_the_list) THEN
            i = mod ( nint ( xkr(1)*nk1 + 2*nk1), nk1 ) + 1
            j = mod ( nint ( xkr(2)*nk2 + 2*nk2), nk2 ) + 1
            k = mod ( nint ( xkr(3)*nk3 + 2*nk3), nk3 ) + 1
            n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
            IF (n>nk .and. equiv(n)==n) THEN
               equiv(n) = nk
               wkk(nk)=wkk(nk)+1.0d0
               s_save(:,:,n) = s(:,:,ns)
            ELSE
               IF (equiv(n)/=nk .or. n<nk ) CALL errore('kpoint_grid', &
                  'something wrong in the checking algorithm',1)
            ENDIF
         ENDIF
         IF ( time_reversal ) THEN
            xx =-xkr(1)*nk1 
            yy =-xkr(2)*nk2 
            zz =-xkr(3)*nk3 
            in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                               .and. abs(zz-nint(zz))<=eps
            IF (in_the_list) THEN
               i = mod ( nint (-xkr(1)*nk1  + 2*nk1), nk1 ) + 1
               j = mod ( nint (-xkr(2)*nk2  + 2*nk2), nk2 ) + 1
               k = mod ( nint (-xkr(3)*nk3  + 2*nk3), nk3 ) + 1
               n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
               IF (n>nk .and. equiv(n)==n) THEN
                  equiv(n) = nk
                  wkk(nk)=wkk(nk)+1.0d0
                  s_save(:,:,n) = -s(:,:,ns)
               ELSE
                  IF (equiv(n)/=nk.or.n<nk) CALL errore('kpoint_grid', &
                  'something wrong in the checking algorithm',2)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF

!  count irreducible points and order them

nks=0
fact=0.0d0
! 
DO nk=1,nkr
  BZtoIBZ(nk) = equiv(nk)
ENDDO
!
DO nk=1,nkr
  IF (equiv(nk) == nk) THEN
    nks=nks+1
    IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
    wk(nks) = wkk(nk)
    fact    = fact+wk(nks)
    !  bring back into to the first BZ
    DO i=1,3
       xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
    ENDDO
    ! DBSP
    BZtoIBZ(nk) = nks  
    ! Change all the one above
    DO ik=nk,nkr
      IF (equiv(ik) == nk) THEN
        BZtoIBZ(ik) = nks
      ENDIF 
    ENDDO
  ENDIF
ENDDO

! Now do the symmetry mapping. 
DO nk=1,nkr
  ! If its an irreducible point 
  IF ( equiv(nk) == nk  ) THEN
    ! Then you have the identity matrix
    s_BZtoIBZ(:,:,nk) = s(:,:,1)
  ELSE
    s_BZtoIBZ(:,:,nk) = s_save(:,:,nk)  
  ENDIF
ENDDO
! 
! Store mapping from all the point to IBZ
 
!  go to cartesian axis (in units 2pi/a0)
CALL cryst_to_cart(nks,xk,bg,1)
!  normalize weights to one
DO nk=1,nks
   wk(nk) = wk(nk)/fact
ENDDO

DEALLOCATE(xkg,wkk)

RETURN
END SUBROUTINE kpoint_grid_epw

