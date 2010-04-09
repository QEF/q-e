!-----------------------------------------------------------------------
subroutine lr_read_d0psi()
  !---------------------------------------------------------------------
  ! ... reads in and stores the vectors necessary to 
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  ! Modified by Osman Baris Malcioglu (2009)
  !
#include "f_defs.h"
  !
  use klist,                only : nks,degauss
  use io_files,             only : prefix
  use lr_variables,         only : d0psi, n_ipol,LR_polarization
  use lr_variables,         only : nwordd0psi, iund0psi
  use wvfct,                only : nbnd, npwx,et
  USE lr_variables,   ONLY : lr_verbosity,restart
   USE io_global,      ONLY : stdout
 !
  implicit none
  !
  ! local variables
  integer :: ip
  character(len=6), external :: int_to_char
  logical :: exst
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_read_d0psi>")')
  endif
  nwordd0psi = 2 * nbnd * npwx * nks
  !
  do ip=1,n_ipol
     !
     if (n_ipol==1) call diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     if (n_ipol==3) call diropn ( iund0psi, 'd0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
     !
     if (.not.exst) then 
      call errore('lr_read_d0psi', TRIM( prefix )//'.d0psi.'//trim(int_to_char(ip))//' not found',1)
     endif
     !
     call davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,-1)
     !
     CLOSE( UNIT = iund0psi)
     !
  end do
  !
end subroutine lr_read_d0psi
!-----------------------------------------------------------------------
