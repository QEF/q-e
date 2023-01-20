program kpoints_unfold

 integer :: nk,sdim1,sdim2,sdim3
 ! counter
 integer :: i, y_n
 DOUBLE PRECISION, ALLOCATABLE   :: kpts_mat(:,:),v2(:), v(:)
 real :: step
 !
 !   
 WRITE(*,*) "Write number of high-symmetry kpts" 
 read(*,*) nk
 WRITE(*,*) "Write high-symmetry kpts (3 columns per row)"
  ALLOCATE(kpts_mat(nk,3))
 DO i=1,nk
      read(*,*) kpts_mat(i,:)
 ENDDO
  !
  !
   WRITE(*,*) "Write x-positions of high-sym kpts"
  ALLOCATE(v2(nk))
  DO i=1,nk
       read(*,*) v2(i)
  ENDDO
  !
  ! The default for the intermediate points  between  high-sym kpts is 1 kpt every 2/300= 0.006667
  ! 
  step = 2/dble(70)
  WRITE(*,*) "Step is (default):", step
  ! 
  ! Find intermediate points
  ALLOCATE(v(nk))
  DO i=1,nk-1
         v(i)=(v2(i+1)-v2(i))/step
  ENDDO
  v(nk)=1 ! this should be always one because is the last one
  WRITE(*,*) "Supercell size (n * m * p) ?"
  read(*,*) sdim1, sdim2, sdim3
  !
  !
   WRITE(*,*) "kpts to use for Supercell calculation:"
    kpts_mat(:,1)=sdim1*kpts_mat(:,1)
    kpts_mat(:,2)=sdim2*kpts_mat(:,2)
    kpts_mat(:,3)=sdim3*kpts_mat(:,3)
    DO i=1,nk
        WRITE(*,'(3F11.6,I4)')  kpts_mat(i,:), abs(nint(v(i)))
    ENDDO

   WRITE(*,*) "Write every single kpt? (0=no, 1=yes)"
   read(*,*) y_n   
   IF (y_n == 1) THEN
     DO i=1,nk-1
       DO j=0,abs(nint(v(i)))-1
        WRITE(*,'(3F11.6,I4)')  kpts_mat(i,1)+ j*(kpts_mat(i+1,1) - kpts_mat(i,1))/abs(nint(v(i))), &
                                kpts_mat(i,2)+ j*(kpts_mat(i+1,2) - kpts_mat(i,2))/abs(nint(v(i))), &
                                kpts_mat(i,3)+ j*(kpts_mat(i+1,3) - kpts_mat(i,3))/abs(nint(v(i))), 1
       ENDDO
     ENDDO
   ENDIF

END program kpoints_unfold
