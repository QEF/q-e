! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

#include "f_defs.h"
!----------------------------------------------------------------------- 
SUBROUTINE wannier_hamiltonian_JK(nwan,hamk,outfile)
!----------------------------------------------------------------------- 
! for Jan Kunis code
  use io_global, only: stdout
  use kinds, only: DP 
  use constants,  ONLY : rytoev
  use klist, only: nks, wk

  implicit none
  integer, intent(in) :: nwan, outfile
  complex(DP) :: hamk(nwan,nwan,nks)

  integer :: i,j, ik
  complex(DP), allocatable :: hamk2(:,:)
  real(DP) :: eps = 1.d-8, hr,hi

  write(stdout,'(/5x,a32,i5,a9)') 'Hamiltonian is in the JK format,', nks, 'k-points'
  write(stdout,'(5x,a48/)') 'ATTENTION: All k-points weights are real weights'

  allocate(hamk2(nwan,nwan))

  write(outfile,*) nks,nwan
  do ik = 1, nks

!		if(ik.eq.43) then
!			write(stdout,*) 'Omitting point', ik
!			CYCLE
!		end if

  	write(outfile,'(f15.12)') wk(ik)

		! eg-orbitals should be the first
	hamk2 = ZERO
!		hamk2(1,:) = hamk(1,:,ik)
!		hamk2(2,:) = hamk(4,:,ik)
!		hamk2(3,:) = hamk(2,:,ik)
!		hamk2(4,:) = hamk(3,:,ik)
!		hamk2(5:nwan,:) = hamk(5:nwan,:,ik)
!		hamk(:,:,ik) = hamk2
!		hamk2(:,1) = hamk(:,1,ik)
!		hamk2(:,2) = hamk(:,4,ik)
!		hamk2(:,3) = hamk(:,2,ik)
!		hamk2(:,4) = hamk(:,3,ik)
!		hamk2(:,5:nwan) = hamk(:,5:nwan,ik)

!rearrange
!		hamk2(1,:) = hamk(5,:,ik)
!		hamk2(2,:) = hamk(3,:,ik)
!		hamk2(3,:) = hamk(1,:,ik)
!		hamk2(4,:) = hamk(2,:,ik)
!		hamk2(5,:) = hamk(4,:,ik)
!		hamk2(6,:) = hamk(10,:,ik)
!		hamk2(7,:) = hamk(8,:,ik)
!		hamk2(8,:) = hamk(6,:,ik)
!		hamk2(9,:) = hamk(7,:,ik)
!		hamk2(10,:) = hamk(9,:,ik)
!		hamk2(11:nwan,:) = hamk(11:nwan,:,ik)
!		hamk(:,:,ik) = hamk2
!		hamk2(:,1) = hamk(:,5,ik)
!		hamk2(:,2) = hamk(:,3,ik)
!		hamk2(:,3) = hamk(:,1,ik)
!		hamk2(:,4) = hamk(:,2,ik)
!		hamk2(:,5) = hamk(:,4,ik)
!		hamk2(:,6) = hamk(:,10,ik)
!		hamk2(:,7) = hamk(:,8,ik)
!		hamk2(:,8) = hamk(:,6,ik)
!		hamk2(:,9) = hamk(:,7,ik)
!		hamk2(:,10) = hamk(:,9,ik)
!		hamk2(:,11:nwan) = hamk(:,11:nwan,ik)

		hamk2 = hamk2 * rytoev
		
		hamk2 = hamk(:,:,ik) * rytoev
		do i=1, nwan
			do j=1, nwan
				hr = ABS(dreal(hamk2(i,j)))
				hi = ABS(aimag(hamk2(i,j)))
				if((hr.ge.eps).AND.(hi.ge.eps)) write(outfile,'(2f12.8)') dreal(hamk2(i,j)), aimag(hamk2(i,j))
				if ((hr.lt.eps).AND.(hi.ge.eps)) write(outfile,'(f3.0,f12.8)') 0., aimag(hamk2(i,j))
				if ((hr.ge.eps).AND.(hi.lt.eps)) write(outfile,'(f12.8,f3.0)') dreal(hamk2(i,j)), 0.
				if ((hr.lt.eps).AND.(hi.lt.eps)) write(outfile,'(2f3.0)') 0., 0.
			end do
		end do

	end do 

!for debug
!	write(stdout,*) 'Real part of first 5x5 block in Gamma'
!	do i=1,5
!			write(stdout,'(5f17.12)') (dreal(hamk2(i,j,1)), j=1,5)
!	end do
!	write(stdout,*) 'Imag part of first 5x5 block in Gamma'
!	do i=1,5
!			write(stdout,'(5f17.12)') (aimag(hamk2(i,j,1)), j=1,5)
!	end do
!end for debug

	deallocate(hamk2)

END SUBROUTINE wannier_hamiltonian_JK

!----------------------------------------------------------------------- 
SUBROUTINE wannier_hamiltonian_IL(nwan,hamk,outfile)
!----------------------------------------------------------------------- 
! Ivan Leonov's code
  use io_global, only: stdout
  use kinds, only: DP 
  use constants,  ONLY : rytoev
	use klist, only: nks
	use ktetra
  use klist, only: xk, wk
	use lsda_mod, only: nspin

	implicit none
	integer, intent(in) :: nwan, outfile
	complex(DP) :: hamk(nwan,nwan,nks)

	integer :: i,j, ik

	write(stdout,*) 'Hamiltonian is in the IL format,', nks, 'k-points'

  write(outfile,*) nks, ntetra
  write(outfile,*) nspin, nwan
  write(outfile,*) (wk(ik), ik=1,nks)
  write(outfile,*) ((xk(i,ik), i=1,3), ik=1,nks)
  write(outfile,*) (1, (tetra(i,j), i=1,4), j=1,ntetra)

	do ik = 1, nks
    write(outfile,*) ((dreal(hamk(i,j,ik)),j=i,nwan),i=1,nwan)
    write(outfile,*) ((dimag(hamk(i,j,ik)),j=i,nwan),i=1,nwan)
	end do 

END SUBROUTINE wannier_hamiltonian_IL
