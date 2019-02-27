!this module contains objected derived and contracted 
MODULE derived_objects

  USE kinds, ONLY : DP

  TYPE prod_proj
!terms <\mathcal{E}_\alpha|u_{k,c}u_{k,v}^*>
!k ponits distributed on MPI tasks as bands object
     INTEGER :: numv !number of valence states (those considered for excitons only) 
     INTEGER :: numc !number of conduction states 
     INTEGER :: nk!total number of k, points 
     INTEGER :: nk_loc!local number of k points 
     INTEGER :: ik_first!first local k point 
     INTEGER :: ik_last!last local k point 
     INTEGER :: ntot_e!dimension of global to all k, basis for KS states
     INTEGER :: nprod_e!number of product terms 
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: javc! (nprod_e,numv,numc,nk_loc)
  END TYPE prod_proj

  TYPE prod_mix
!terms <\mathcal{E}_\alpha|(u_{k,v}u_{k',v'}^*)>
!and terms <\mathcal{E}_\alpha|(u_{k,c}u_{k',c'}^*)>
!k' distributed over MPI tasks
!k NOT distributed
     INTEGER :: numv !number of valence states (those considered for excitons only) 
     INTEGER :: numc !number of conduction states
     INTEGER :: nk!total number of k, points
     INTEGER :: nk_loc!local number of k points
     INTEGER :: ik_first!first local k point
     INTEGER :: ik_last!last local k point
     INTEGER :: ntot_e!dimension of global to all k, basis for KS states
     INTEGER :: nprod_e!number of product terms
     COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: gvv! (nprod_e,numv,nk,numv',nk_loc)! ' means relative to nk_loc
     COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: gcc! (nprod_e,numc,nk,numc',nk_loc)


  END TYPE prod_mix


  CONTAINS

    SUBROUTINE initialize_prod_proj(pp)
      implicit none
      TYPE(prod_proj) :: pp
      
      nullify(pp%javc)

      return
    END SUBROUTINE initialize_prod_proj



    SUBROUTINE deallocate_prod_proj(pp)
      implicit none
      TYPE(prod_proj) :: pp
      if(associated(pp%javc)) deallocate(pp%javc)
      nullify(pp%javc)
      return
    END SUBROUTINE deallocate_prod_proj






    SUBROUTINE initialize_prod_mix(pm)
      implicit none
      TYPE(prod_mix) :: pm

      nullify(pm%gvv)
      nullify(pm%gcc)

      return
    END SUBROUTINE initialize_prod_mix



    SUBROUTINE deallocate_prod_mix(pm)
      implicit none
      TYPE(prod_mix) :: pm
      if(associated(pm%gvv)) deallocate(pm%gvv)
      nullify(pm%gvv)
      if(associated(pm%gcc)) deallocate(pm%gcc)
      nullify(pm%gcc)
      return
    END SUBROUTINE deallocate_prod_mix
    

    SUBROUTINE build_prod_proj(bd,pd,pp)
!this subroutine constructs the prod_proj object
      USE simple_objects, ONLY : bands,product

      

      implicit none
      TYPE(bands), INTENT(in) :: bd
      TYPE(product), INTENT(in) :: pd
      TYPE(prod_proj), INTENT(out) :: pp

      
      INTEGER :: ik,iv,ic
      COMPLEX(kind=DP), ALLOCATABLE :: tmp_mat(:,:,:),tmp_mat2(:,:),tmp_mat3(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: zmat(:,:,:,:),zmat0(:,:),zmat1(:,:,:,:)
      LOGICAL, parameter :: debug = .false.


      pp%numv=bd%numv
      pp%numc=bd%numc
      pp%ntot_e=bd%ntot_e
      pp%nk=bd%nk
      pp%nk_loc=bd%nk_loc
      pp%ik_first=bd%ik_first
      pp%ik_last=bd%ik_last
      pp%nprod_e=pd%nprod_e
      
      if(pp%nk_loc>0) then
         allocate(pp%javc(pp%nprod_e,pp%numv, pp%numc,pp%nk_loc))
         allocate(tmp_mat(pp%nprod_e,pp%ntot_e,pp%numc))
         allocate(tmp_mat2(pp%ntot_e,pp%numv))
         allocate(tmp_mat3(pp%ntot_e,pp%numc))
         do ik=1,pp%nk_loc
            tmp_mat2(1:pp%ntot_e,1:pp%numv)=conjg(bd%omat(1:pp%ntot_e,1:pp%numv,ik))
            tmp_mat3(1:pp%ntot_e,1:pp%numc)=bd%omat(1:pp%ntot_e,pp%numv+1:pp%numv+pp%numc,ik)
            call ZGEMM('N','N',pp%nprod_e*pp%ntot_e,pp%numc,pp%ntot_e,(1.d0,0.d0),pd%fij,&
                 &pp%nprod_e*pp%ntot_e,tmp_mat3,pp%ntot_e,(0.d0,0.d0),tmp_mat,pp%nprod_e*pp%ntot_e)
            do ic=1,pp%numc
               call ZGEMM('N','N',pp%nprod_e,pp%numv,pp%ntot_e,(1.d0,0.d0),tmp_mat(1,1,ic),pp%nprod_e,tmp_mat2,pp%ntot_e,&
                   &(0.d0,0.d0),pp%javc(1,1,ic,ik),pp%nprod_e)
            enddo
         enddo
         deallocate(tmp_mat,tmp_mat2,tmp_mat3)
      else
         nullify(pp%javc)
      endif
      if(debug) then
         !test for consistency
         allocate(zmat(pp%numc,pp%numv,pp%numv,pp%numc))
         allocate(zmat0(bd%num,bd%num))
         allocate(zmat1(pd%ntot_e,pd%ntot_e,pd%ntot_e,pd%ntot_e))
         do ik=1,bd%nk_loc
            call ZGEMM('C','N',bd%num,bd%num,bd%ntot_e,(1.d0,0.d0),bd%omat(1,1,ik),bd%ntot_e,bd%omat(1,1,ik),bd%ntot_e,&
                 &(0.d0,0.d0),zmat0,bd%num)
            
            do iv=1,bd%num
               do ic=1,bd%num
                  write(*,*) 'CHECK OMAT ik, ic,iv', ik, ic, iv, zmat0(ic,iv)
               enddo
            enddo
         enddo
         call ZGEMM('C','N',pd%ntot_e*pd%ntot_e,pd%ntot_e*pd%ntot_e,pd%nprod_e,(1.d0,0.d0),pd%fij,&
              pd%nprod_e,pd%fij,pd%nprod_e,(0.d0,0.d0),zmat1,pd%ntot_e*pd%ntot_e)
         do iv=1,pd%ntot_e
            do ic=1,pd%ntot_e
               write(*,*) 'CHECK  FIJ ic,iv', ic, iv, zmat1(iv,ic,iv,ic)
            enddo
         enddo
        
         do ik=1,pp%nk_loc
            call ZGEMM('C','N',pp%numc*pp%numv,pp%numc*pp%numv,pp%nprod_e,(1.d0,0.d0),pp%javc(1,1,1,ik),pp%nprod_e,&
                 & pp%javc(1,1,1,ik),pp%nprod_e,(0.d0,0.0),zmat, pp%numc*pp%numv)
            do iv=1,pp%numv
               do ic=1,pp%numc
                  write(*,*) 'CHECK JAVC ik, ic,iv', ik, ic, iv, zmat(iv,ic,iv,ic)
               enddo
            enddo
         enddo
          deallocate(zmat,zmat0,zmat1)
       endif
      return
    END SUBROUTINE build_prod_proj

  SUBROUTINE build_prod_mix(sin,bd,pd,pm,pt)
!this subroutine constructs the prod_mix object
!COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: gvv! (nprod_e,numv,nk,numv,nk_loc)                                     
!COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: gcc! (nprod_e,numc,nk,numc,nk_loc)
 
     USE input_simple_exc
     USE simple_objects, ONLY : bands,product,potential
     USE mp_world, ONLY : mpime, world_comm
     USE mp, ONLY : mp_sum, mp_bcast
     USE io_global, ONLY : stdout

      implicit none

      TYPE(input_options) :: sin
      TYPE(bands), INTENT(in) :: bd
      TYPE(product), INTENT(in) :: pd
      TYPE(prod_mix), INTENT(out) :: pm
      TYPE(potential) :: pt

      INTEGER :: ik,iv,ic
      COMPLEX(kind=DP), ALLOCATABLE :: tmp_mat(:,:,:),tmp_mat2(:,:),tmp_mat3(:,:,:)
      LOGICAL, parameter :: debug = .true.
      COMPLEX(kind=DP), ALLOCATABLE :: emat(:,:)
      INTEGER :: is_mine
      INTEGER :: jk
      COMPLEX(kind=DP), ALLOCATABLE :: tmp_pot(:,:),tmp_fij(:,:,:)
      INTEGER :: ii,jj,kk

      pm%numv=bd%numv
      pm%numc=bd%numc
      pm%ntot_e=bd%ntot_e
      pm%nk=bd%nk
      pm%nk_loc=bd%nk_loc
      pm%ik_first=bd%ik_first
      pm%ik_last=bd%ik_last
      pm%nprod_e=pd%nprod_e
      if(pm%nk_loc>0) then
         allocate( pm%gvv(pm%nprod_e,pm%numv,pm%nk,pm%numv,pm%nk_loc))
         allocate( pm%gcc(pm%nprod_e,pm%numc,pm%nk,pm%numc,pm%nk_loc))
      else
         nullify(pm%gcc)
         nullify(pm%gvv)
      endif
!now do gvv
!loop on nk
!if ik is my distribute to others
!loop on nk_loc
!calculate terms
      allocate(emat(pm%ntot_e,pm%numv))
      allocate(tmp_mat2(pm%ntot_e,pm%numv))
      allocate(tmp_mat(pm%nprod_e,pm%ntot_e,pm%numv))
      allocate(tmp_mat3(pm%nprod_e,pm%numv,pm%numv))
      allocate(tmp_pot(pm%nprod_e,pm%nprod_e))
      allocate(tmp_fij(pm%nprod_e,pm%ntot_e,pm%ntot_e))

      



      do ik=1,pm%nk
         if(ik>=pm%ik_first.and.ik<=pm%ik_last) then
            emat(1:pm%ntot_e,1:pm%numv)=bd%omat(1:pm%ntot_e,1:pm%numv,ik-pm%ik_first+1)
            is_mine=mpime+1
         else
            is_mine=0
         endif
         call mp_sum(is_mine,world_comm)
         is_mine=is_mine-1
         call mp_bcast( emat,is_mine, world_comm )

         do jk=1,pm%nk_loc!on k' local

!find out q=k_i-k_j
            ii=pt%ijk(1,ik,jk+pm%ik_first-1)+1
            jj=pt%ijk(2,ik,jk+pm%ik_first-1)+1
            kk=pt%ijk(3,ik,jk+pm%ik_first-1)+1
            tmp_pot(1:pm%nprod_e,1:pm%nprod_e)= pt%vpotq(1:pm%nprod_e,1:pm%nprod_e,ii,jj,kk)
         
           
            
            !tmp_pot(1:pm%nprod_e,1:pm%nprod_e)=1.d0!DEBUG
            if(sin%h_level >= 3) then
               tmp_pot(1:pm%nprod_e,1:pm%nprod_e)= tmp_pot(1:pm%nprod_e,1:pm%nprod_e) +&
              &pt%wpotq(1:pm%nprod_e,1:pm%nprod_e,ii,jj,kk)
            endif
            call ZGEMM('N','N',pm%nprod_e,pm%ntot_e*pm%ntot_e,pm%nprod_e,(-1.d0,0.d0),&
                 &tmp_pot,pm%nprod_e,pd%fij,pm%nprod_e,(0.d0,0.d0),tmp_fij,pm%nprod_e)

    
            tmp_mat2(1:pm%ntot_e,1:pm%numv)=conjg(bd%omat(1:pm%ntot_e,1:pm%numv,jk))
            !call ZGEMM('N','N',pm%nprod_e*pm%ntot_e,pm%numv,pm%ntot_e,(1.d0,0.d0),pd%fij,&
            !      &pm%nprod_e*pm%ntot_e,emat,pm%ntot_e,(0.d0,0.d0),tmp_mat,pm%nprod_e*pm%ntot_e)

            call ZGEMM('N','N',pm%nprod_e*pm%ntot_e,pm%numv,pm%ntot_e,(1.d0,0.d0),tmp_fij,&
                  &pm%nprod_e*pm%ntot_e,emat,pm%ntot_e,(0.d0,0.d0),tmp_mat,pm%nprod_e*pm%ntot_e)


            do iv=1,pm%numv
                call ZGEMM('N','N',pm%nprod_e,pm%numv,pm%ntot_e,(1.d0,0.d0),tmp_mat(1,1,iv),pm%nprod_e,tmp_mat2,pm%ntot_e,&
                   &(0.d0,0.d0),tmp_mat3(1,1,iv),pm%nprod_e)
               
            enddo
            do iv=1,pm%numv
               pm%gvv(1:pm%nprod_e,1:pm%numv,ik,iv,jk)=tmp_mat3(1:pm%nprod_e,iv,1:pm%numv)
            enddo
         enddo
      enddo

      deallocate(emat)
      deallocate(tmp_mat2)
      deallocate(tmp_mat)
      deallocate(tmp_mat3)
      deallocate(tmp_pot)
      deallocate(tmp_fij)

      allocate(emat(pm%ntot_e,pm%numc))
      allocate(tmp_mat2(pm%ntot_e,pm%numc))
      allocate(tmp_mat(pm%nprod_e,pm%ntot_e,pm%numc))
      allocate(tmp_mat3(pm%nprod_e,pm%numc,pm%numc))

      do ik=1,pm%nk
         if(ik>=pm%ik_first.and.ik<=pm%ik_last) then
            emat(1:pm%ntot_e,1:pm%numc)=bd%omat(1:pm%ntot_e,pm%numv+1:pm%numv+pm%numc,ik-pm%ik_first+1)
            is_mine=mpime+1
         else
            is_mine=0
         endif
         call mp_sum(is_mine,world_comm)
         is_mine=is_mine-1
         call mp_bcast( emat,is_mine, world_comm )

         do jk=1,pm%nk_loc!on k' local                                                                                         
            tmp_mat2(1:pm%ntot_e,1:pm%numc)=conjg(bd%omat(1:pm%ntot_e,pm%numv+1:pm%numv+pm%numc,jk))
            call ZGEMM('N','N',pm%nprod_e*pm%ntot_e,pm%numc,pm%ntot_e,(1.d0,0.d0),pd%fij,&
                  &pm%nprod_e*pm%ntot_e,emat,pm%ntot_e,(0.d0,0.d0),tmp_mat,pm%nprod_e*pm%ntot_e)
            do ic=1,pm%numc
                call ZGEMM('N','N',pm%nprod_e,pm%numc,pm%ntot_e,(1.d0,0.d0),tmp_mat(1,1,ic),pm%nprod_e,tmp_mat2,pm%ntot_e,&
                   &(0.d0,0.d0),tmp_mat3(1,1,ic),pm%nprod_e)

            enddo
            do ic=1,pm%numc
               pm%gcc(1:pm%nprod_e,1:pm%numc,ik,ic,jk)=tmp_mat3(1:pm%nprod_e,ic,1:pm%numc)
            enddo
         enddo
      enddo

      deallocate(emat)
      deallocate(tmp_mat2)
      deallocate(tmp_mat)
      deallocate(tmp_mat3)




    END SUBROUTINE build_prod_mix



END MODULE derived_objects
