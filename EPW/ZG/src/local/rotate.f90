PROGRAM rotate
  !
IMPLICIT NONE

 INTEGER                         :: i, p, sym, nk_tot, ctr
 DOUBLE PRECISION, ALLOCATABLE   :: str_f(:,:), str_f_new(:,:), str_fp(:, :)
 CHARACTER(len=256)              :: filename
 !
 DOUBLE PRECISION                :: Rmat(2,2), pi, theta, eps
 ! 
 pi         = 3.1415926536
 eps        = 1d-5
 !
 sym        = 6
 WRITE(*,*) "Write number of rotations based on symmetry (integer)"
 READ(*,*)  sym
 theta      = 2.0 * pi / FLOAT(sym)
 WRITE(*,*) "Write number of entries (integer)"
 READ(*,*)  nk_tot
 !
 ALLOCATE(str_f(nk_tot, 4), str_fp(nk_tot, 4), str_f_new(nk_tot * sym, 4))
 str_f = 0.d0
 !
 filename = 'structure_factor_all-phonon.dat'
 OPEN (unit = 80, file = filename, status = 'unknown', form = 'formatted')
 DO p = 1, nk_tot
    READ(80,*) str_f(p, :)
 END DO
 CLOSE(80)
 !
 ! To remove double contribution upon rotation
 str_fp = 0.d0
 str_fp(1, :) = str_f(1, :)
 DO p = 2, nk_tot
   IF (ATAN(str_f(p, 1) / str_f(p, 2)) .LT. (2.d0 * pi / float(sym) - eps) ) THEN
   str_fp(p, :) = str_f(p, :)
   ENDIF
 ENDDO
 !
 OPEN (unit = 80, file = filename, status = 'unknown', form = 'formatted')
 !
 str_f_new = 0.d0
 ctr = 1
 DO i = 0, sym - 1
   Rmat(1, :) = (/ COS(i * theta), -SIN(i * theta) /)
   Rmat(2, :) = (/ SIN(i * theta),  COS(i * theta) /)
   DO p = 1, nk_tot
     str_f_new(ctr, 1 : 2) = MATMUL(Rmat, str_fp(p, 1 : 2))
     str_f_new(ctr, 3) = str_fp(p, 3)
     str_f_new(ctr, 4) = str_fp(p, 4)
     ctr = ctr + 1
   END DO
 ENDDO
 CLOSE(80)
 ! 
 filename = 'structure_factor_all-phonon_rot.dat'
 OPEN (unit = 80, file = filename, status = 'unknown', form = 'formatted')
 !
 ctr = 1
 DO i = 1, sym
   DO p = 1, nk_tot
     WRITE(80,'(4f26.6)') str_f_new(ctr, :)
     ctr = ctr + 1
   END DO
 ENDDO
CLOSE(80)
! 
DEALLOCATE(str_f, str_f_new)
END PROGRAM
