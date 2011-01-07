!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_lr( )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  ! Modified by O. Baris Malcioglu (2009)
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_end, mp_barrier
  !
  USE parallel_include
  use lr_variables,         only : n_ipol, LR_polarization, beta_store
  use lr_variables,         only :  gamma_store, zeta_store, norm0, rho_1_tot
  use lr_variables,         only : lr_verbosity, itermax, bgz_suffix
  USE io_global,            ONLY : ionode
  use io_files,             only : tmp_dir, prefix
  USE io_global,      ONLY : stdout
  ! For gaussian cube file
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau
  USE cell_base,  ONLY : celldm, at, bg
  !
  IMPLICIT NONE
  !
  character(len=6), external :: int_to_char
  !
  character(len=256) :: filename
  !
  integer :: ip,i
  !
  !
  
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<stop_lr>")')
  endif
  ! I write the beta gamma and z coefficents to output directory for
  ! easier post processing. These can also be read from the output log file
#ifdef __PARA
  if (ionode) then
#endif
  !
  do ip=1,n_ipol
   if (n_ipol==3) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(ip))
   if (n_ipol==1) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
   filename = trim(tmp_dir) // trim(filename)
  !
  !
  open (158, file = filename, form = 'formatted', status = 'replace')
  !
  write(158,*) itermax
  !
  write(158,*) norm0(ip)
  !
  do i=1,itermax
     !
     write(158,*) beta_store(ip,i)
     write(158,*) gamma_store(ip,i)
     !This is absolutely necessary for cross platform compatibilty
     do j=1,n_ipol                                    
      write(158,*) zeta_store (ip,j,i)
     end do
     !
  end do
  !
  close(158)
  !
  enddo
#ifdef __PARA
  end if
#endif
  !
  !   Deallocate lr variables
  !
  CALL lr_dealloc()


  CALL mp_barrier()
  !
  CALL mp_end()
  !
#if defined (__T3E)
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
  !
#endif
  !
  STOP
  !
  !
END SUBROUTINE stop_lr
