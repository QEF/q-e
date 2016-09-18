!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE create_vcprim(cpp, istate,uu, qm)
!this subroutine creates the (v)cprim structure
!starting from the strategy beta
! Z_{\alpha,iv}=v_{\alpha,i'j'}U_{i,i'}U_{j, j'}
!it is working with parallelized qm

  USE kinds, ONLY : DP
  USE basic_structures,     ONLY : cprim_prod,q_mat,wannier_u,wannier_P, free_memory

  implicit none
  TYPE(cprim_prod), INTENT(out) :: cpp!the structure to be created
  INTEGER, INTENT(in)  :: istate !state i to be calculated
  TYPE(wannier_u), INTENT(in) :: uu!for the energies and trasformation matrix   
  TYPE(q_mat), INTENT(in) :: qm ! for S matrices 

  INTEGER :: i,j,k,ip, iw
  REAL(kind=DP), ALLOCATABLE :: v_mat(:,:),u_mat(:,:)


!initialization

  call free_memory(cpp)

  cpp%cprim=istate
  cpp%nums=uu%nums
  cpp%nums_occ=uu%nums_occ(1)
  cpp%nums_cond=cpp%nums-cpp%nums_occ
  cpp%numpw=qm%numpw
  cpp%is_parallel=qm%is_parallel
  cpp%numpw_para=qm%numpw_para
  cpp%first_para=qm%first_para

  allocate(cpp%cpmat(cpp%numpw_para,cpp%nums))
  cpp%lda=cpp%numpw_para

!calculate V_{\alpha,v'}=\sum_i' v_{\alpha,i'v'}U_{i,i'}
  allocate(v_mat(cpp%numpw_para,cpp%nums))
  v_mat(:,:)=0.d0
  if(cpp%cprim <= cpp%nums_occ)then
!i==valence state
!WE ASSUME SAME ORDER
      do ip=1,qm%wp(1)%numij
         i=qm%wp(1)%ij(1,ip)!valence only                                             
         j=qm%wp(1)%ij(2,ip)!valence and conduction         
         do iw=1,cpp%numpw_para
            v_mat(iw,j)=v_mat(iw,j)+qm%wp(iw)%o(ip)*dble(uu%umat(cpp%cprim,i,1))
         enddo
         if( j<=cpp%nums_occ .and. i/=j ) then
            do iw=1,cpp%numpw_para
               v_mat(iw,i)=v_mat(iw,i)+qm%wp(iw)%o(ip)*dble(uu%umat(cpp%cprim,j,1))
            enddo
         endif
      enddo
  else
!i==conduction state
     do ip=1,qm%wp(1)%numij
!WE ASSUME SAME ORDER
        i=qm%wp(1)%ij(1,ip)!valence only                                                                                    
        j=qm%wp(1)%ij(2,ip)!valence and conduction                                                                          
        if(j>uu%nums_occ(1)) then
           do iw=1,cpp%numpw_para
              v_mat(iw,i)=v_mat(iw,i)+qm%wp(iw)%o(ip)*dble(uu%umat(cpp%cprim,j,1))
           enddo
        endif
     enddo
  endif

!calculate Z_{\alpha j}=v_{\alpha, j'}U_{j,j'}

  allocate(u_mat(uu%nums,uu%nums))
  u_mat(:,:)=dble(uu%umat(:,:,1))
  call dgemm('N','T',cpp%numpw_para,uu%nums,uu%nums,1.d0,v_mat,cpp%numpw_para,u_mat,uu%nums,0.d0,cpp%cpmat,cpp%numpw_para)


  deallocate(u_mat,v_mat)


  return

END SUBROUTINE create_vcprim


SUBROUTINE add_vcprim_conduction(cpp, uu, up, vp)
!this subroutine adds to the (v)cprim structure 
!the contribution from conduction states                                                                              
!starting from the strategy beta        
                                                                                      
! Z_{\alpha,iv}+=v_{\alpha,i'j'}U'_{i,i'}U_{j, j'}  
                                                                            

  USE kinds, ONLY : DP
  USE basic_structures,     ONLY : cprim_prod,wannier_u,wannier_u_prim,v_pot_prim
  USE io_global,  ONLY : stdout

  implicit none
  TYPE(cprim_prod), INTENT(inout) :: cpp!the structure to be calaculated 
  TYPE(wannier_u), INTENT(in) :: uu!for the energies and trasformation matrix
  TYPE(wannier_u_prim), INTENT(in) :: up!for the U'_{cc'} trasform
  TYPE(v_pot_prim), INTENT(in) :: vp!


  INTEGER :: i,j,k,ip, iw
  REAL(kind=DP), ALLOCATABLE :: v_mat(:,:),u_mat(:,:)

  if(cpp%cprim <= cpp%nums_occ) return
  if(cpp%numpw_para /= vp%numpw_para) then
     write(stdout,*) 'add_vcprim_conduction NOT CORRESPONDING'
     FLUSH(stdout)
     stop
  endif

!calculate V_{\alpha,v'}=\sum_i' v_{\alpha,i'v'}U_{i,i'}                                                                    
  allocate(v_mat(cpp%numpw_para,cpp%nums))
  v_mat(:,:)=0.d0
  
!WE ASSUME SAME ORDER                                                     
                                                   
  do ip=1,vp%numpw_prim
     i=vp%ij(1,ip)! on c'   
     j=vp%ij(2,ip)! on c         
     do iw=1,cpp%numpw_para
        v_mat(iw,j)=v_mat(iw,j)+vp%vmat(ip,iw)*dble(up%umat(cpp%cprim-cpp%nums_occ,i))
     enddo
     
  enddo
  

!calculate Z_{\alpha j}+=v_{\alpha, j'}U_{j,j'}  
                                                                              

  allocate(u_mat(uu%nums,uu%nums))
  u_mat(:,:)=dble(uu%umat(:,:,1))
  call dgemm('N','T',cpp%numpw_para,uu%nums,uu%nums,1.d0,v_mat,cpp%numpw_para,u_mat,uu%nums,1.d0,cpp%cpmat,cpp%numpw_para)


  deallocate(u_mat,v_mat)

  return

END SUBROUTINE add_vcprim_conduction


SUBROUTINE distribute_qmat(qm,qmd)
!this subroutine distributes q_mat on parallel processors
  USE kinds, ONLY : DP
  USE basic_structures,    ONLY : q_mat,wannier_P
  USE mp_world,            ONLY : nproc,mpime

  implicit none

  TYPE(q_mat), INTENT(in) :: qm ! input q_mat
  TYPE(q_mat), INTENT(out) :: qmd ! output distributed q_mat

  INTEGER :: l_blk,nbegin,nend,ii

!set up qmd
  qmd%numpw=qm%numpw
  qmd%is_parallel=.true.


  l_blk= qm%numpw/nproc
  if(l_blk*nproc < qm%numpw) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > qm%numpw) nend = qm%numpw

  qmd%numpw_para=nend-nbegin+1
  qmd%first_para=nbegin

  if(qmd%numpw_para>1) then
     allocate(qmd%wp(qmd%numpw_para))
     do ii=1,qmd%numpw_para
        qmd%wp(ii)%numij=qm%wp(ii+qmd%first_para-1)%numij
        allocate(qmd%wp(ii)%ij(2,qmd%wp(ii)%numij))
        qmd%wp(ii)%ij(:,:)=qm%wp(ii+qmd%first_para-1)%ij(:,:)
        allocate(qmd%wp(ii)%o(qmd%wp(ii)%numij))
        qmd%wp(ii)%o(:)=qm%wp(ii+qmd%first_para-1)%o(:)
     enddo

  endif
  return

END SUBROUTINE distribute_qmat


SUBROUTINE collect_cprim_prod(cpp,cppd)
!this subroutine collects the cprim structures from parallel processor

  USE kinds, ONLY : DP
  USE basic_structures,     ONLY : cprim_prod, free_memory
  USE mp_world,             ONLY : nproc,mpime, world_comm!group
  USE io_global,            ONLY : stdout
  USE parallel_include

  implicit none

  TYPE(cprim_prod), INTENT(out) :: cpp!structure to be collected
  TYPE(cprim_prod), INTENT(in) :: cppd!distributed structures

  REAL(kind=DP), ALLOCATABLE :: sndbuf(:)
  INTEGER :: l_blk,ii, ierr

!initializations

  call free_memory(cpp)

  cpp%cprim=cppd%cprim
  cpp%nums=cppd%nums
  cpp%nums_occ=cppd%nums_occ
  cpp%nums_cond=cpp%nums-cpp%nums_occ
  cpp%numpw=cppd%numpw
  cpp%is_parallel=.false.
  cpp%numpw_para=cpp%numpw
  cpp%first_para=1
  
  if(.not.cppd%is_parallel) then
     write(stdout,*) 'collect_cprim_prod: NOT CORRESPONDING'
     FLUSH(stdout)
     stop
  endif

  l_blk= cpp%numpw/nproc
  if(l_blk*nproc < cpp%numpw) l_blk = l_blk+1


  allocate(cpp%cpmat(nproc*l_blk,cpp%nums))
  cpp%lda=nproc*l_blk

  allocate(sndbuf(l_blk))


  do ii=1,cpp%nums
      
     sndbuf(:)=0.d0
     sndbuf(1:cppd%numpw_para)=cppd%cpmat(1:cppd%numpw_para,ii)

#if defined(__MPI)
     call MPI_ALLGATHER(sndbuf,l_blk,MPI_DOUBLE_PRECISION,cpp%cpmat(:,ii),l_blk,MPI_DOUBLE_PRECISION,&
          world_comm, ierr)

#else
     cpp%cpmat(:,ii)=cppd%cpmat(:,ii)

#endif  
  enddo
  deallocate(sndbuf)

return

END SUBROUTINE collect_cprim_prod


SUBROUTINE distribute_v_pot_prim(vp,vpd)
!this subroutine distributes the structure v_pot_prim among processors
!on the basis of the polarization

  USE kinds, ONLY : DP
  USE basic_structures,     ONLY : v_pot_prim, free_memory
  USE mp_world,             ONLY : nproc,mpime, world_comm!group
  USE io_global,            ONLY : stdout
  USE parallel_include

  implicit none

  TYPE(v_pot_prim), INTENT(in) :: vp!structure to be distribute
  TYPE(v_pot_prim), INTENT(out) :: vpd!distribute structure

  INTEGER :: l_blk,nbegin,nend,ii

!initializations

  vpd%numpw=vp%numpw
  vpd%numpw_prim=vp%numpw_prim
  vpd%is_parallel=.true.
 

  l_blk= vp%numpw/nproc
  if(l_blk*nproc < vp%numpw) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > vp%numpw) nend = vp%numpw
 
  vpd%numpw_para=nend-nbegin+1
  vpd%first_para=nbegin

  allocate(vpd%ij(2,vp%numpw_prim))
  vpd%ij(:,:)=vp%ij(:,:)
   
  allocate(vpd%vmat(vpd%numpw_prim,vpd%numpw_para))
  do ii=1,vpd%numpw_para
     vpd%vmat(:,ii)=vp%vmat(:,ii+vpd%first_para-1)
  enddo

  return

END SUBROUTINE distribute_v_pot_prim
