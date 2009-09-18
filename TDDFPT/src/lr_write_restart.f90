!-----------------------------------------------------------------------
subroutine lr_write_restart()
  !---------------------------------------------------------------------
  ! ... reads in and stores the vectors necessary to 
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  use io_files,             only : tmp_dir, prefix
  use lr_variables,         only : beta_store, gamma_store, zeta_store, norm0, &
                                   LR_polarization, LR_iteration, n_ipol
  USE io_global,            ONLY : ionode
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout
  !
  implicit none
  !
  character(len=6), external :: int_to_char
  !
  !integer, intent(in) :: pol, iter
  !
  ! local variables
  !
  integer :: i, pol_index
  character(len=256) :: tempfile, filename
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_write_restart>")')
  endif
  pol_index=1
  if ( n_ipol /= 1 ) pol_index=LR_polarization

#ifdef __PARA
  if (ionode) then
#endif
  !
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
  !
  !
  open (158, file = tempfile, form = 'formatted', status = 'unknown')
  !
  write(158,*) LR_iteration
  !
  write(158,*) norm0(pol_index)
  !
  do i=1,LR_iteration
     !
     write(158,*) beta_store(pol_index,i)
     write(158,*) gamma_store(pol_index,i)
     write(158,*) zeta_store (pol_index,:,i)
     !
  end do
  !
  close(158)
  !
#ifdef __PARA
  end if
#endif
  !
  return
end subroutine lr_write_restart
!-----------------------------------------------------------------------
