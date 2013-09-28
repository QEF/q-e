!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE write_wfcfile(filename,wfc,elaux,num)
!-------------------------------------------------------------------------
!
!  This subroutine writes a formatted wavefunction file
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, ionode_id
  USE ld1inc, ONLY  : grid
  USE radial_grids, ONLY  : ndmx
  USE mp,      ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm

  IMPLICIT NONE
  INTEGER, INTENT(in) :: num   ! the number of wavefunctions to write
  REAL(DP), INTENT(in) :: wfc(ndmx,num)
  CHARACTER(len=2), INTENT(in) :: elaux(num)
  CHARACTER(len=256), INTENT(in) :: filename
  INTEGER :: ios, n, ns

  IF (filename == ' ') RETURN

  IF (ionode) &
     OPEN(unit=19,file=filename, status='unknown', iostat=ios, err=80)
80 CALL mp_bcast(ios, ionode_id,world_comm)
  CALL errore('write_wfcfile','opening file '//trim(filename),abs(ios))
  IF (ionode) THEN
     WRITE(19,'("#",12x,"r",38(18x,a2))') (elaux(n),n=1,num)
     DO n=1,grid%mesh
        WRITE(19,'(38f20.12)') grid%r(n), (wfc(n,ns), ns=1,num)
     ENDDO
     CLOSE(19)
  ENDIF
  RETURN
END SUBROUTINE write_wfcfile

!-------------------------------------------------------------------------
SUBROUTINE write_wfcfile_ft(filename_in,wfc,num)
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
USE mp_world,      ONLY : world_comm

IMPLICIT NONE
INTEGER, INTENT(in) :: num   ! the number of wavefunctions to write
real(DP), INTENT(in) :: wfc(ndmx,num)
CHARACTER(len=256), INTENT(in) :: filename_in

INTEGER :: ios, n, ns, nmax
CHARACTER(len=256) :: filename
real(DP) :: q, fac, pi, wrk(ndmx), jlq(ndmx), norm(num), normr(num), work(num)
real(DP), EXTERNAL :: int_0_inf_dr   ! the function calculating the integral

IF (filename_in == ' ') RETURN

IF (ionode) THEN
   filename = trim(filename_in)//'.q'
   OPEN(unit=19,file=filename, status='unknown', iostat=ios)
   filename = trim(filename_in)//'.norm_q'
   OPEN(unit=29,file=filename, status='unknown', iostat=ios)
   pi = 4._dp*atan(1._dp)
   DO ns=1,num
      wrk(1:grid%mesh)=(wfc(1:grid%mesh,ns)*exp(-0.04*grid%r2(1:grid%mesh)))**2
      normr(ns) = int_0_inf_dr ( wrk, grid, grid%mesh, 2*lls(ns)+2 )
   ENDDO
   !write (*,*) normr(1:num)
   fac = pi /grid%r(grid%mesh) * 2._dp/pi
   norm(:) = 0._dp
   nmax = int (10 * grid%r(grid%mesh) /pi)
   DO n=1,nmax
      q=n * pi /grid%r(grid%mesh)
      DO ns=1,num
         CALL sph_bes ( grid%mesh, grid%r, q, lls(ns), jlq )
         wrk(1:grid%mesh)=jlq(1:grid%mesh)*wfc(1:grid%mesh,ns)*  &
               exp(-0.04*grid%r2(1:grid%mesh))*grid%r(1:grid%mesh)
         work(ns) = int_0_inf_dr ( wrk, grid, grid%mesh, 2*lls(ns)+2 )
      ENDDO
      norm(1:num) = norm(1:num) + work(1:num)*work(1:num)*q*q*fac
      WRITE (19,'(15f12.6)') q, (work(ns), ns=1,num)
      WRITE (29,'(15f12.6)') q, (norm(ns)/normr(ns), ns=1,num)
   ENDDO
   CLOSE (29)
   CLOSE (19)
   !
   ! end of Fourier analysis
   !
ENDIF
RETURN
END SUBROUTINE write_wfcfile_ft
!
!-------------------------------------------------------------------------
SUBROUTINE write_efun(filename,fun,x,npte,num)
!-------------------------------------------------------------------------
!
!  This subroutine writes a functions defined on the energy grid
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id
USE mp,      ONLY : mp_bcast
USE mp_world,      ONLY : world_comm

IMPLICIT NONE
INTEGER, INTENT(in) :: num   ! the number of wavefunctions to write
INTEGER, INTENT(in) :: npte   ! the number of energies
real(DP), INTENT(in) :: fun(npte,num), x(npte)
CHARACTER(len=256), INTENT(in) :: filename
INTEGER :: ios, n, ns

IF (filename == ' ') RETURN

IF (ionode) &
   OPEN(unit=19,file=filename, status='unknown', iostat=ios, err=800)
800  CALL mp_bcast(ios, ionode_id, world_comm)
CALL errore('write_wfcfile','opening file '//trim(filename),abs(ios))
IF (ionode) THEN
   DO n=1,npte
      WRITE(19,'(38f20.12)') x(n), (max(min(fun(n,ns),9.d4),-9d4), ns=1,num)
   ENDDO
   CLOSE(19)
ENDIF

RETURN
END SUBROUTINE write_efun

