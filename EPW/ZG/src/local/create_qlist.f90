program qpoint_gen
!
IMPLICIT NONE
 
 integer :: i,j,k,n,ctr, ctrAB
 integer :: nk1, nk2, nk3, nk_tot !, k1, k2, k3
 DOUBLE PRECISION, ALLOCATABLE   :: xkg(:,:), xkg_AB(:,:)
 DOUBLE PRECISION                :: q_B(3), q_A(3), eps
 CHARACTER(len=256)              :: filename
 !
 eps = 1.0E-06
 !
 read(*,*) nk1, nk2, nk3
 !
 nk_tot=nk1*nk2*nk3
 !
 ALLOCATE(xkg(3,nk_tot))
 DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 ! + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 ! + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 ! + dble(k3)/2/nk3 ! k1 , k2 , k3 is for the shift
        ENDDO
     ENDDO
 ENDDO
 !
 ctr= 0
 ctrAB = 0
 DO i=1,nk_tot
   q_A = xkg(:,i) +  xkg(:,i) ! q_A to find if q belongs in A 
   IF (((abs(q_A(1)) .LT. eps) .OR. (abs(abs(q_A(1))-1) .LT. eps)) .AND. &
       ((abs(q_A(2)) .LT. eps) .OR. (abs(abs(q_A(2))-1) .LT. eps)) .AND. &
       ((abs(q_A(3)) .LT. eps) .OR. (abs(abs(q_A(3))-1) .LT. eps))) THEN
       ctrAB=ctrAB+1
   ELSE
    DO j=i+1,nk_tot
       q_B = xkg(:,i) +  xkg(:,j)       
      IF (((abs(q_B(1)) .LT. eps) .OR. (abs(abs(q_B(1))-1) .LT. eps)) .AND. &
          ((abs(q_B(2)) .LT. eps) .OR. (abs(abs(q_B(2))-1) .LT. eps)) .AND. &
          ((abs(q_B(3)) .LT. eps) .OR. (abs(abs(q_B(3))-1) .LT. eps))) THEN
          ctr = ctr+1
          ctrAB=ctrAB+1
      END IF 
     END DO
   END IF
 END DO
 ! 
 ALLOCATE(xkg_AB(3,ctrAB))
 ctr=0
 ctrAB = 0
 DO i=1,nk_tot
   q_A = xkg(:,i) +  xkg(:,i) ! q_A to find if q belongs in A 
   IF (((abs(q_A(1)) .LT. eps) .OR. (abs(abs(q_A(1))-1) .LT. eps)) .AND. &
       ((abs(q_A(2)) .LT. eps) .OR. (abs(abs(q_A(2))-1) .LT. eps)) .AND. &
       ((abs(q_A(3)) .LT. eps) .OR. (abs(abs(q_A(3))-1) .LT. eps))) THEN
       ctrAB=ctrAB+1
       xkg_AB(:,ctrAB) = xkg(:,i)
 !      write(*,*) "A", xkg(:,i)
   ELSE
    DO j=i+1,nk_tot
       q_B = xkg(:,i) +  xkg(:,j)       
      IF (((abs(q_B(1)) .LT. eps) .OR. (abs(abs(q_B(1))-1) .LT. eps)) .AND. &
          ((abs(q_B(2)) .LT. eps) .OR. (abs(abs(q_B(2))-1) .LT. eps)) .AND. &
          ((abs(q_B(3)) .LT. eps) .OR. (abs(abs(q_B(3))-1) .LT. eps))) THEN
          ctr = ctr+1
          ctrAB=ctrAB+1
          xkg_AB(:,ctrAB) = xkg(:,i)
 !         write(*,*) "B", xkg(:,i)
      END IF 
     END DO
   END IF
 END DO
 !
 !
 filename = 'qlist_AB.txt'
 OPEN (unit = 80, file = filename, status = 'unknown', form = 'formatted')
 WRITE(*,*) "points in set B and set AB:", ctr, ctrAB
 WRITE(80,'(10i10)') ctrAB
 DO i = 1, ctrAB 
    WRITE(80,'(3F10.6,4i2)') xkg_AB(:,i), 1
 END DO
  close(80)
 !
 DEALLOCATE(xkg,xkg_AB)
 !
end program qpoint_gen
