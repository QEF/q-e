MODULE roctx
  USE ISO_C_BINDING,   ONLY : c_char, c_null_char
  IMPLICIT NONE
  character(len=256),private :: tempName
  INTEGER :: ret

#if defined(__PROFILE_ROCTX)
 INTERFACE 
     SUBROUTINE roctxMarkA(message) BIND(c, name="roctxMarkA")
       USE ISO_C_BINDING,   ONLY : C_CHAR
       IMPLICIT NONE
       CHARACTER(C_CHAR) :: message(*)
     END SUBROUTINE roctxMarkA

     FUNCTION roctxRangePushA(message) BIND(c, name="roctxRangePushA")
       USE ISO_C_BINDING,   ONLY: C_INT,C_CHAR
       IMPLICIT NONE
       INTEGER(C_INT) :: roctxRangePushA
       CHARACTER(C_CHAR) :: message(*)
     END FUNCTION roctxRangePushA

     SUBROUTINE roctxRangePop() BIND(c, name="roctxRangePop")
       IMPLICIT NONE
     END SUBROUTINE roctxRangePop
 END INTERFACE
#endif

 CONTAINS
 
      SUBROUTINE roctxStartRange(name)
         character(kind=c_char,len=*) :: name

#if defined(__PROFILE_ROCTX)
         tempName=trim(name)//c_null_char
         ret = roctxRangePushA(tempname)
#endif
      END SUBROUTINE 

      SUBROUTINE roctxEndRange()
#if defined(__PROFILE_ROCTX)
         CALL roctxRangePop()
#endif
      END SUBROUTINE

END MODULE roctx
