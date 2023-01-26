MODULE test_tools
   USE kinds,     ONLY : DP
   INTERFACE para_write
      MODULE PROCEDURE &
         para_write_i, para_write_i1, para_write_i2, para_write_i3, para_write_i4, para_write_i5, para_write_i6, &
         para_write_r, para_write_r1, para_write_r2, para_write_r3, para_write_r4, para_write_r5, para_write_r6,&
         para_write_c, para_write_c1, para_write_c2, para_write_c3, para_write_c4, para_write_c5, para_write_c6
   END INTERFACE
CONTAINS
   SUBROUTINE para_write_i(A, filename)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN) :: filename
   END SUBROUTINE
   SUBROUTINE para_write_i1(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A(:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(1), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      WRITE(23747806, '(I5)') (A(j), j=1, bounds(1))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_i2(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      integer, INTENT(IN) :: A(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(2), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      !print *, bounds
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      WRITE(23747806, '(I10)') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_i3(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A(:,:,:)
   END SUBROUTINE
   SUBROUTINE para_write_i4(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A(:,:,:,:)

   END SUBROUTINE
   SUBROUTINE para_write_i5(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A(:,:,:,:,:)

   END SUBROUTINE
   SUBROUTINE para_write_i6(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A(:,:,:,:,:,:)
   END SUBROUTINE

   SUBROUTINE para_write_r(A)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A
   END SUBROUTINE
   SUBROUTINE para_write_r1(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(1), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds(1), 1
      WRITE(23747806, '(G20.6E4)') (A(j), j=1, bounds(1))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_r2(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(2), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      WRITE(23747806, '(G20.6E4)') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_r3(A)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:,:,:)
   END SUBROUTINE
   SUBROUTINE para_write_r4(A)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:,:,:,:)

   END SUBROUTINE
   SUBROUTINE para_write_r5(A)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:,:,:,:,:)

   END SUBROUTINE
   SUBROUTINE para_write_r6(A)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: A(:,:,:,:,:,:)
   END SUBROUTINE

   SUBROUTINE para_write_c(A)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A
   END SUBROUTINE
   SUBROUTINE para_write_c1(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(1), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds(1), 1
      WRITE(23747806, '(2G20.6E4)') (A(j), j=1, bounds(1))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_c2(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(2), i, j

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      !WRITE(23747806, '("(",G20.6E4,"+",G20.6E4,"j",")")') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      WRITE(23747806, '(2G20.6E4)') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_c3(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:,:,:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(3), i, j, k, l

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      !WRITE(23747806, '("(",G20.6E4,"+",G20.6E4,"j",")")') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      WRITE(23747806, '(2G20.6E4)') (((A(j,i, k), j=1, bounds(1)), i=1, bounds(2)), k=1, bounds(3))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_c4(A, varName)
      USE mp_world,      ONLY : mpime
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:,:,:,:)
      CHARACTER(LEN=*), INTENT(IN) :: varName
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: lineLen

      INTEGER :: bounds(4), i, j, k, l

      bounds = SHAPE(A)
      WRITE(filename,'(a, i0)') varName, mpime
      OPEN(UNIT=23747806, FILE='./test_out/'//filename, POSITION="append", FORM='formatted')
      WRITE(23747806, *) REPEAT('-',10), bounds
      !WRITE(23747806, '("(",G20.6E4,"+",G20.6E4,"j",")")') ((A(j,i), j=1, bounds(1)), i=1, bounds(2))
      WRITE(23747806, '(2G20.6E4)') ((((A(j,i, k, l), j=1, bounds(1)), i=1, bounds(2)), k=1, bounds(3)),l=1, bounds(4))
      !WRITE(lineLen,'("(",i0,"G20.6E4)")') bounds(1) ! x2 due to complex
      !do i = 1, bounds(2)
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")') (A(j,i), j=1, bounds(1))
      !WRITE(23747806, '(G20.6E4,"+",G20.6E4,"j")', advance="no") A(:,i)
      !WRITE(23747806, lineLen) (A(j,i), j=1, bounds(1))
      !end do
      CLOSE(23747806)
   END SUBROUTINE
   SUBROUTINE para_write_c5(A)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:,:,:,:,:)

   END SUBROUTINE
   SUBROUTINE para_write_c6(A)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: A(:,:,:,:,:,:)
   END SUBROUTINE
END MODULE



