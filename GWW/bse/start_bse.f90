subroutine start_bse
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N]
  !
  !  Wrapper routine for postprocessing initialization
  !
  USE mp_global,     ONLY: mp_startup
  USE environment,   ONLY: environment_start
  implicit none
  character(len=9) :: code = 'BSE'
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
  ! 
  return
end subroutine start_bse
