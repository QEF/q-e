!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


SUBROUTINE v_basis(numpw,o_basis,cutoff)

!this subroutine calculate the coulomb v operator on the space of polarizability basis
!functions and retains only eigenstate larger than cutoff


  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE cell_base,            ONLY : tpiba2,tpiba
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec
  USE gvect,                ONLY : g, gstart
  USE wvfct,                ONLY : npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE mp_world,             ONLY : world_comm
  USE klist,                ONLY : xk
  USE wannier_gw,           ONLY : l_truncated_coulomb,truncation_radius
  USE exx,                  ONLY : exx_divergence, exx_grid_init, yukawa
  USE constants,            ONLY : e2,fpi

  implicit none

  
  INTEGER, INTENT(inout) :: numpw!dimension of polarization basis      
  COMPLEX(kind=DP), INTENT(inout) :: o_basis(npw,numpw)
  REAL(kind=DP), INTENT(in) :: cutoff!cutoff for v eigen values


  REAL(kind=DP), ALLOCATABLE :: fac(:)
  INTEGER :: ig,ii,jj
  REAL(kind=DP) :: qq
  REAL(kind=DP), ALLOCATABLE :: vmat(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: vo_basis(:,:)
  REAL(kind=DP), ALLOCATABLE :: eigen(:),vectors(:,:)
  INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
  INTEGER, ALLOCATABLE :: isuppz(:)
  INTEGER :: n_found
  REAL(kind=DP), ALLOCATABLE :: work(:)
  INTEGER :: lwork,info,liwork


  allocate(fac(npw))

  do ig=1,npw
      qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0

      if(.not.l_truncated_coulomb) then

         if (qq > 1.d-8) then
            fac(ig)=e2*fpi/(tpiba2*qq + yukawa )
         else
            fac(ig)= 0.d0
            if (yukawa .gt. 1.d-8 ) then
               fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
            end if
         end if
      else
         if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
         endif
         !if(abs(fac(ig)) <=1.d-8) fac(ig)=1.d-8!ATTENZIONE
      endif

   end do
   allocate(vmat(numpw,numpw))
   allocate(vo_basis(npw,numpw))
   
   do ii=1,numpw
      vo_basis(:,ii)=fac(:)*o_basis(:,ii)
   enddo
   deallocate(fac)

   call dgemm('T','N',numpw,numpw,2*npw,2.d0,o_basis,2*npw,vo_basis,2*npw,0.d0,vmat,numpw)
   if(gstart==2) then
      do ii=1,numpw
         do jj=1,numpw
            vmat(jj,ii)=vmat(jj,ii)-dble(conjg(o_basis(1,jj))*vo_basis(1,ii))
         enddo
      enddo
   endif
   do ii=1,numpw
      call mp_sum(vmat(:,ii),world_comm)
   enddo

  
   allocate(eigen(numpw))
   allocate(vectors(numpw,numpw))
   if(ionode) then
      allocate(isuppz(2*numpw))
      allocate(work(1),iwork(1))
      call DSYEVR('V','V','U',numpw,vmat,numpw,cutoff,1.d5,1,1,0.d0,n_found,eigen,&
           & vectors,numpw,isuppz,work, -1,iwork,-1, info)
      lwork=work(1)
      liwork=iwork(1)
      deallocate(work,iwork)
      allocate(work(lwork))
      allocate(iwork(liwork))
      call DSYEVR('V','V','U',numpw,vmat,numpw, cutoff,1.d5,1,1,0.d0,n_found,eigen,&
           & vectors,numpw,isuppz,work,lwork,iwork,liwork, info)
      if(info/=0) then
         write(stdout,*) 'ROUTINE v_basis DSYEVR, INFO:', info
         stop
      endif
      deallocate(isuppz)
      deallocate(work,iwork)
   else
      eigen(:)=0.d0
      vectors(:,:)=0.d0
      n_found=0
   endif
   call mp_sum(n_found,world_comm)
   call mp_sum(eigen(1:n_found),world_comm)
   do ii=1,n_found
      write(stdout,*) 'v_basis:',ii,eigen(ii)
      call mp_sum(vectors(:,ii),world_comm)
   enddo
   deallocate(vmat)

   vo_basis(:,:)=o_basis(:,:)
   call dgemm('N','N',2*npw,n_found,numpw,1.d0,vo_basis,2*npw,vectors,numpw,0.d0,o_basis,2*npw)
   numpw=n_found

   deallocate(eigen,vectors,vo_basis)
 
  return
END SUBROUTINE v_basis
