! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!-----------------------------------------------------------------------
SUBROUTINE wannier_hamiltonian_JK(nwan,hamk,outfile)
!-----------------------------------------------------------------------
! for Jan Kunis code
  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE constants,  ONLY : rytoev
  USE klist, ONLY: nks, wk

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nwan, outfile
  COMPLEX(DP) :: hamk(nwan,nwan,nks)

  INTEGER :: i,j, ik
  COMPLEX(DP), ALLOCATABLE :: hamk2(:,:)
  real(DP) :: eps = 1.d-8, hr,hi

  WRITE(stdout,'(/5x,a32,i5,a9)') 'Hamiltonian is in the JK format,', nks, 'k-points'
  WRITE(stdout,'(5x,a48/)') 'ATTENTION: All k-points weights are real weights'

  ALLOCATE(hamk2(nwan,nwan))

  WRITE(outfile,*) nks,nwan
  DO ik = 1, nks

     ! if(ik.eq.43) then
     !    write(stdout,*) 'Omitting point', ik
     !    CYCLE
     ! end if

     WRITE(outfile,'(f15.12)') wk(ik)

     ! eg-orbitals should be the first
     hamk2 = ZERO
     !hamk2(1,:) = hamk(1,:,ik)
     !hamk2(2,:) = hamk(4,:,ik)
     !hamk2(3,:) = hamk(2,:,ik)
     !hamk2(4,:) = hamk(3,:,ik)
     !hamk2(5:nwan,:) = hamk(5:nwan,:,ik)
     !hamk(:,:,ik) = hamk2
     !hamk2(:,1) = hamk(:,1,ik)
     !hamk2(:,2) = hamk(:,4,ik)
     !hamk2(:,3) = hamk(:,2,ik)
     !hamk2(:,4) = hamk(:,3,ik)
     !hamk2(:,5:nwan) = hamk(:,5:nwan,ik)

!rearrange
!hamk2(1,:) = hamk(5,:,ik)
!hamk2(2,:) = hamk(3,:,ik)
!hamk2(3,:) = hamk(1,:,ik)
!hamk2(4,:) = hamk(2,:,ik)
!hamk2(5,:) = hamk(4,:,ik)
!hamk2(6,:) = hamk(10,:,ik)
!hamk2(7,:) = hamk(8,:,ik)
!hamk2(8,:) = hamk(6,:,ik)
!hamk2(9,:) = hamk(7,:,ik)
!hamk2(10,:) = hamk(9,:,ik)
!hamk2(11:nwan,:) = hamk(11:nwan,:,ik)
!hamk(:,:,ik) = hamk2
!hamk2(:,1) = hamk(:,5,ik)
!hamk2(:,2) = hamk(:,3,ik)
!hamk2(:,3) = hamk(:,1,ik)
!hamk2(:,4) = hamk(:,2,ik)
!hamk2(:,5) = hamk(:,4,ik)
!hamk2(:,6) = hamk(:,10,ik)
!hamk2(:,7) = hamk(:,8,ik)
!hamk2(:,8) = hamk(:,6,ik)
!hamk2(:,9) = hamk(:,7,ik)
!hamk2(:,10) = hamk(:,9,ik)
!hamk2(:,11:nwan) = hamk(:,11:nwan,ik)

     hamk2 = hamk2 * rytoev

     hamk2 = hamk(:,:,ik) * rytoev
     DO i=1, nwan
        DO j=1, nwan
           hr = abs(dreal(hamk2(i,j)))
           hi = abs(aimag(hamk2(i,j)))
           IF((hr>=eps).and.(hi>=eps)) WRITE(outfile,'(2f12.8)') dreal(hamk2(i,j)), aimag(hamk2(i,j))
           IF ((hr<eps).and.(hi>=eps)) WRITE(outfile,'(f3.0,f12.8)') 0., aimag(hamk2(i,j))
           IF ((hr>=eps).and.(hi<eps)) WRITE(outfile,'(f12.8,f3.0)') dreal(hamk2(i,j)), 0.
           IF ((hr<eps).and.(hi<eps)) WRITE(outfile,'(2f3.0)') 0., 0.
        ENDDO
     ENDDO

  ENDDO

  !for debug
  ! write(stdout,*) 'Real part of first 5x5 block in Gamma'
  ! do i=1,5
  !    write(stdout,'(5f17.12)') (dreal(hamk2(i,j,1)), j=1,5)
  ! end do
  ! write(stdout,*) 'Imag part of first 5x5 block in Gamma'
  ! do i=1,5
  !    write(stdout,'(5f17.12)') (aimag(hamk2(i,j,1)), j=1,5)
  ! end do
  !end for debug

  DEALLOCATE(hamk2)

END SUBROUTINE wannier_hamiltonian_JK

!-----------------------------------------------------------------------
SUBROUTINE wannier_hamiltonian_IL(nwan,hamk,outfile)
!-----------------------------------------------------------------------
! Ivan Leonov's code
  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE constants,  ONLY : rytoev
  USE klist, ONLY: nks
  USE ktetra
  USE klist, ONLY: xk, wk
  USE lsda_mod, ONLY: nspin

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nwan, outfile
  COMPLEX(DP) :: hamk(nwan,nwan,nks)

  INTEGER :: i,j, ik

  WRITE(stdout,*) 'Hamiltonian is in the IL format,', nks, 'k-points'

  WRITE(outfile,*) nks, ntetra
  WRITE(outfile,*) nspin, nwan
  WRITE(outfile,*) (wk(ik), ik=1,nks)
  WRITE(outfile,*) ((xk(i,ik), i=1,3), ik=1,nks)
  WRITE(outfile,*) (1, (tetra(i,j), i=1,4), j=1,ntetra)

  DO ik = 1, nks
     WRITE(outfile,*) ((dreal(hamk(i,j,ik)),j=i,nwan),i=1,nwan)
     WRITE(outfile,*) ((dimag(hamk(i,j,ik)),j=i,nwan),i=1,nwan)
  ENDDO

END SUBROUTINE wannier_hamiltonian_IL
