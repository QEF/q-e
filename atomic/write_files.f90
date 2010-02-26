!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
subroutine write_wfcfile(filename,wfc,elaux,num)
!-------------------------------------------------------------------------
!
!  This subroutine writes a formatted wavefunction file
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id
USE ld1inc, ONLY  : grid
USE radial_grids, ONLY  : ndmx
USE mp,      ONLY : mp_bcast

implicit none
integer, intent(in) :: num   ! the number of wavefunctions to write  
real(DP), intent(in) :: wfc(ndmx,num)
character(len=2), intent(in) :: elaux(num)
character(len=256), intent(in) :: filename
integer :: ios, n, ns

if (filename .eq. ' ') return

if (ionode) &
   open(unit=19,file=filename, status='unknown', iostat=ios, err=800)
800  call mp_bcast(ios, ionode_id)
call errore('write_wfcfile','opening file '//TRIM(filename),abs(ios))
if (ionode) then
   write(19,'("#     r  ",38(14x,a2))') (elaux(n),n=1,num)
   do n=1,grid%mesh
      write(19,'(38f20.12)') grid%r(n), (wfc(n,ns), ns=1,num)
   enddo
   close(19)
endif

return
end subroutine write_wfcfile

!-------------------------------------------------------------------------
subroutine write_wfcfile_ft(filename_in,wfc,num)
!-------------------------------------------------------------------------
!
! calculate and print the fourier transform of wfc and the  
! convergence of their norm with cutoff
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id
USE ld1inc, ONLY  : grid, lls
USE radial_grids, ONLY  : ndmx
USE mp,      ONLY : mp_bcast

implicit none
integer, intent(in) :: num   ! the number of wavefunctions to write  
real(DP), intent(in) :: wfc(ndmx,num)
character(len=256), intent(in) :: filename_in

integer :: ios, n, ns, nmax
character(len=256) :: filename
real(DP) :: q, fac, pi, wrk(ndmx), jlq(ndmx), norm(num), normr(num), work(num)
real(DP), external :: int_0_inf_dr   ! the function calculating the integral 

IF (filename_in .eq. ' ') RETURN

IF (ionode) THEN
   filename = trim(filename_in)//'.q'
   open(unit=19,file=filename, status='unknown', iostat=ios)
   filename = trim(filename_in)//'.norm_q'
   open(unit=29,file=filename, status='unknown', iostat=ios)
   pi = 4._dp*atan(1._dp)
   do ns=1,num
      wrk(1:grid%mesh)=(wfc(1:grid%mesh,ns)*exp(-0.04*grid%r2(1:grid%mesh)))**2
      normr(ns) = int_0_inf_dr ( wrk, grid, grid%mesh, 2*lls(ns)+2 )
   end do
   !write (*,*) normr(1:num)
   fac = pi /grid%r(grid%mesh) * 2._dp/pi
   norm(:) = 0._dp
   nmax = int (10 * grid%r(grid%mesh) /pi)
   do n=1,nmax 
      q=n * pi /grid%r(grid%mesh)
      do ns=1,num
         call sph_bes ( grid%mesh, grid%r, q, lls(ns), jlq )
         wrk(1:grid%mesh)=jlq(1:grid%mesh)*wfc(1:grid%mesh,ns)*  &
               exp(-0.04*grid%r2(1:grid%mesh))*grid%r(1:grid%mesh)
         work(ns) = int_0_inf_dr ( wrk, grid, grid%mesh, 2*lls(ns)+2 )
      end do
      norm(1:num) = norm(1:num) + work(1:num)*work(1:num)*q*q*fac
      write (19,'(15f12.6)') q, (work(ns), ns=1,num)
      write (29,'(15f12.6)') q, (norm(ns)/normr(ns), ns=1,num)
   end do
   close (29)
   close (19)
   !
   ! end of Fourier analysis
   !
endif
return
end subroutine write_wfcfile_ft
!
!-------------------------------------------------------------------------
subroutine write_efun(filename,fun,x,npte,num)
!-------------------------------------------------------------------------
!
!  This subroutine writes a functions defined on the energy grid
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id
USE mp,      ONLY : mp_bcast

implicit none
integer, intent(in) :: num   ! the number of wavefunctions to write  
integer, intent(in) :: npte   ! the number of energies  
real(DP), intent(in) :: fun(npte,num), x(npte)
character(len=256), intent(in) :: filename
integer :: ios, n, ns

if (filename .eq. ' ') return

if (ionode) &
   open(unit=19,file=filename, status='unknown', iostat=ios, err=800)
800  call mp_bcast(ios, ionode_id)
call errore('write_wfcfile','opening file '//TRIM(filename),abs(ios))
if (ionode) then
   do n=1,npte
      write(19,'(38f20.12)') x(n), (max(min(fun(n,ns),9.d4),-9d4), ns=1,num)
   enddo
   close(19)
endif

return
end subroutine write_efun

