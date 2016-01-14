!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutine contains the routines dedicated to obtaining optimal basis sets


SUBROUTINE optimal_driver(num_in,wfcs,lda,options,num_out, info)
!this routine is a driver for performing the calculation of the optimal basis set
!using the appropriate method

  USE kinds,                ONLY : DP
  USE wannier_gw,           ONLY : optimal_options
  USE io_global,            ONLY : stdout

  implicit none

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(lda,num_in)!in input non-orthonormal in output optimal basis
  INTEGER, INTENT(in) :: lda!leading dimension of wfcs, essentially npw or npwx
  TYPE(optimal_options), INTENT(in) :: options!options to be used
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions
  INTEGER, INTENT(out) :: info!final outcome status 0== OK
  
  REAL(kind=DP) :: tr

  !select routine
  select case (options%idiago)

  case(0)
!Gram_Schmidt like
     if(options%l_complete) then
        tr=0.d0
     else
        tr=options%thres
     endif
     call optimal_gram_schmidt(num_in,wfcs,lda,options%ithres,tr,num_out)
     
  case default
     write(stdout,*) 'optimal driver: NOT IMPLEMENTED YET'
     FLUSH(stdout)
     stop
  end select

  info=0
  return
END SUBROUTINE optimal_driver






SUBROUTINE optimal_gram_schmidt(num_in,wfcs,lda,ithres,thres,num_out)
!this subroutine performs a gram_schmidt orthonormalization and retains
!vectors which are above the give threshold

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE wvfct,                ONLY : npwx, npw
  USE gvect,                ONLY : gstart
 
 implicit none

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(lda,num_in)!in input non-orthonormal in output optimal basis
  INTEGER, INTENT(in) :: lda!leading dimension of wfcs, essetally npw or npwx
  INTEGER, INTENT(in) :: ithres!kind of threshold
  REAL(kind=DP), INTENT(in) :: thres!thrshold for the optimal basis
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions



  INTEGER :: i,j
  REAL(kind=DP), ALLOCATABLE :: prod(:)
  REAL(kind=DP) :: sca
  REAL(kind=DP), EXTERNAL :: ddot


  allocate(prod(num_in))
  num_out=0

  do i=1,num_in
     if(num_out >0) then
        call dgemv('T',2*npw,num_out,2.d0, wfcs,2*lda,wfcs(:,i),1,0.d0,prod,1)
        if(gstart==2) then
           prod(1:num_out)=prod(1:num_out) - dble(wfcs(1,1:num_out)*conjg(wfcs(1,i)))
        endif
        call mp_sum(prod(1:num_out),world_comm)
        call dgemm('N','N',2*npw,1,num_out,-1.d0,wfcs,2*lda,prod,num_in,1.d0,wfcs(:,i),2*lda)
     endif
     sca = 2.d0*ddot(2*npw,wfcs(:,i),1,wfcs(:,i),1)
     if(gstart==2) then
        sca=sca-dble((wfcs(1,i)*conjg(wfcs(1,i))))
     endif
     call mp_sum(sca,world_comm)
     if(sca >= thres) then
        num_out=num_out+1
        sca=dsqrt(sca)
        call dcopy(2*npw,wfcs(:,i),1,wfcs(:,num_out),1)
        wfcs(1:npw,num_out)=wfcs(1:npw,num_out)/sca
     endif
  enddo

  deallocate(prod)
  return
END SUBROUTINE optimal_gram_schmidt
