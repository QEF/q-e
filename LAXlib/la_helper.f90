!----------------------------------------------------------------------------
SUBROUTINE laxlib_free_ortho_group()
  !----------------------------------------------------------------------------
  !
  use mp_diag
  IMPLICIT NONE
#if defined (__MPI)
  CALL clean_ortho_group ( )
#endif
  RETURN
  !
END SUBROUTINE laxlib_free_ortho_group

