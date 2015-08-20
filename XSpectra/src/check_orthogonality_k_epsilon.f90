subroutine check_orthogonality_k_epsilon( xcoordcrys, xang_mom )
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : bg, at
  USE xspectra,        ONLY : xepsilon, xkvec
  USE io_global,       ONLY : stdout
  implicit none
  LOGICAL, INTENT(IN) :: xcoordcrys
  INTEGER, INTENT(IN) :: xang_mom
! internal
  INTEGER :: i
  REAL(kind=dp)       :: norm,xeps_dot_xk
  

  IF ( xcoordcrys ) CALL cryst_to_cart(1,xepsilon,at,1)
  norm=DSQRT(xepsilon(1)**2+xepsilon(2)**2+xepsilon(3)**2)
  DO i=1,3
     xepsilon(i)=xepsilon(i)/norm
  ENDDO
  
  ! ... If needed normalize xkvec
  !     and check orthogonality with xepsilon
  
  IF(xang_mom.EQ.2) THEN
     
     IF ( xcoordcrys ) CALL cryst_to_cart(1,xkvec,at,1)
     norm=DSQRT(xkvec(1)**2+xkvec(2)**2+xkvec(3)**2)
     xkvec(:)=xkvec(:)/norm
     xeps_dot_xk=xkvec(1)*xepsilon(1)+&
          xkvec(2)*xepsilon(2)+&
          xkvec(3)*xepsilon(3)
     IF ((ABS(xeps_dot_xk)) .gt. 1.0d-6) THEN
        WRITE(stdout,'(5x,a)') &
             'ERROR: xkvec and xepsilon are not orthogonal'
        WRITE(stdout,'(12x,a,f10.6,/)') 'scalar product=', xeps_dot_xk
        WRITE(stdout,'(5x,a)') 'STOP'
        CALL stop_xspectra ()
     ENDIF
     
  ENDIF
  
  
  
end subroutine check_orthogonality_k_epsilon
