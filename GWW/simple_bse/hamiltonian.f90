SUBROUTINE hamiltonian(sin,n,bd,pp,pt,pm,a,ha,ilevel)
!this subroutine applies the excitonic hamiltonian

  USE simple_objects
  USE derived_objects 
  USE input_simple_exc
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global, ONLY : stdout

  implicit none

  TYPE(input_options), INTENT(in) :: sin
  INTEGER :: n!number of vectors
  TYPE(bands), INTENT(in) :: bd
  TYPE(prod_proj), INTENT(in) :: pp
  TYPE(potential), INTENT(in) :: pt
  TYPE(prod_mix), INTENT(in) :: pm
  TYPE(exc), INTENT(in) :: a(n)
  TYPE(exc), INTENT(inout) :: ha(n)
  INTEGER :: ilevel!for analysis:0, All, 1 Diagonal, 2 Exchange ,3 Direct

  INTEGER :: i,ik,jk,ii,jj,kk
  COMPLEX(kind=DP), ALLOCATABLE :: emat(:,:),vemat(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: hmat(:,:,:), imat(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: bvc(:,:),bvc_t(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: gcc_t(:,:,:),imat_t(:,:,:), gvv_t(:,:,:)
  COMPLEX(kind=DP) :: fact

  INTEGER :: ia,ikk,iv,ic,icp,ivp

  
!Diagonal part  
  do i=1,n
     ha(i)= bd .hd. a(i)
  end do
  if(ilevel>1) ha(1)%avc(1:ha(1)%numv,1:ha(1)%numc,1:ha(1)%nk_loc)=(0.d0,0.d0)
  if(sin%h_level >= 1 .and. sin%spin_state>=1 .and. (ilevel/=1)) then
!Exchange term
     allocate(emat(pp%nprod_e,n),vemat(pp%nprod_e,n))
     emat=(0.d0,0.d0)
     if(sin%spin_state==1) then
        fact=(2.d0,0.d0)
     else
        fact=(1.d0,0.d0)
     endif

     do i=1,n
        if(pp%nk_loc>0) then
           call ZGEMM('N','N',pp%nprod_e,1,pp%numc*pp%numv*pp%nk_loc,(1.d0,0.d0),pp%javc,pp%nprod_e,&
                &a(i)%avc,pp%numc*pp%numv*pp%nk_loc,(0.d0,0.d0),emat(1,i),pp%nprod_e)
        else
           emat(1:pp%nprod_e,i)=(0.d0,0.d0)
        endif
     end do
     call mp_sum(emat,world_comm)
     call ZGEMM('N','N',pp%nprod_e,n,pp%nprod_e,fact,pt%vpot,pt%nprod_e,emat,pp%nprod_e,(0.d0,0.d0),&
          &vemat, pp%nprod_e)
     do i=1,n
        if(pp%nk_loc>0) then
           call ZGEMM('C','N',pp%numc*pp%numv*pp%nk_loc,1,pp%nprod_e,(1.d0,0.d0),pp%javc,pp%nprod_e,vemat(1,i),pp%nprod_e,&
                (1.d0,0.d0),ha(i)%avc,pp%numc*pp%numv*pp%nk_loc)
        endif
     enddo

     deallocate(emat,vemat)
  endif
  if(ilevel==3) ha(1)%avc=0.d0
  if(sin%h_level >= 2 .and. ( ilevel==0 .or. ilevel==3)) then
!TD-HF term
     allocate(hmat(pm%nprod_e,pm%numv,pm%numv))  
     allocate(imat(pm%nprod_e,pm%numv,pm%numc))
     allocate(bvc(pm%numv,pm%numc))
     allocate(imat_t(pm%nprod_e,pm%numc, pm%numv))
     allocate(gcc_t(pm%nprod_e,pm%numc, pm%numc))
     allocate(bvc_t(pm%numc,pm%numv))
     allocate(gvv_t(pm%nprod_e,pm%numv,pm%numv))
     !write(stdout,*) 'DEBUG',pm%nk,pm%nk_loc

     do i=1,n

!loop on k
        do ik=1,pm%nk
           !loop on k'
           bvc=(0.d0,0.d0)
           do jk=1,pm%nk_loc
!multiply V*Gvv'
              ii=pt%ijk(1,ik,jk+pm%ik_first-1)+1
              jj=pt%ijk(2,ik,jk+pm%ik_first-1)+1
              kk=pt%ijk(3,ik,jk+pm%ik_first-1)+1

              

!change indices
              do ia=1,pm%nprod_e
                 do iv=1,pm%numv
                    do ivp=1,pm%numv
                      ! gvv_t(ia,iv,ivp)=pm%gvv(ia,iv,ik,ivp,jk)
                       hmat(ia,iv,ivp)=pm%gvv(ia,iv,ik,ivp,jk)
                    enddo
                 enddo
              enddo
             ! hmat=1.d0!DEBUG

  !            call ZGEMM('N','N',pm%nprod_e,pm%numv*pm%numv, pm%nprod_e,&
  !   &(-1.d0,0.d0),pt%vpotq(1,1,ii,jj,kk),pm%nprod_e,gvv_t, pm%nprod_e,(0.d0,0.d0),hmat,pm%nprod_e)
!in case add W_c term
  !            if(sin%h_level >= 3) then
  !               call ZGEMM('N','N',pm%nprod_e,pm%numv*pm%numv, pm%nprod_e,&
  !   &(-1.d0,0.d0),pt%wpotq(1,1,ii,jj,kk),pm%nprod_e,gvv_t, pm%nprod_e,(1.d0,0.d0),hmat,pm%nprod_e)
  !            endif
!multiply *Akv
!multiply Gcc'*
              call ZGEMM('N','N',pm%nprod_e*pm%numv,pm%numc,pm%numv,(1.d0,0.d0)&
     &,hmat,pm%nprod_e*pm%numv,a(i)%avc(1,1,jk),pm%numv,(0.d0,0.d0),imat,pm%nprod_e*pm%numv)
!now invert indices
              do ia=1,pm%nprod_e
                 do icp=1,pm%numc
                    do iv=1,pm%numv
                       imat_t(ia,icp,iv)=imat(ia,iv,icp)
                    enddo
                 enddo
              enddo

              do ia=1,pm%nprod_e
                 do icp=1,pm%numc
                    do ic=1,pm%numc
                       gcc_t(ia,icp,ic)=pm%gcc(ia,ic,ik,icp,jk)
                    enddo
                 enddo
              enddo

              call ZGEMM('C','N',pm%numc,pm%numv,pm%nprod_e*pm%numc,(1.d0,0.d0),&
        &gcc_t,pm%nprod_e*pm%numc,imat_t,pm%nprod_e*pm%numc,(0.d0,0.d0),&
  &bvc_t,pm%numc)


!put indices in correct order
              do iv=1,pm%numv
                 do ic=1,pm%numc
                    bvc(iv,ic)=bvc(iv,ic)+bvc_t(ic,iv)
                 enddo
              enddo
           enddo!on jk
           call mp_sum(bvc,world_comm)
           !DEBUG START
           !bvc=bvc/dble(pm%nk)
           !DEBUG END
           if(ik>=pm%ik_first .and. ik<=pm%ik_last) then
              ha(i)%avc(1:pm%numv,1:pm%numc,ik-pm%ik_first+1)=ha(i)%avc(1:pm%numv,1:pm%numc,ik-pm%ik_first+1)+&
     &bvc(1:pm%numv,1:pm%numc)
           endif
           
        enddo!ok ik
     
     enddo!on i

     deallocate(hmat,imat,bvc)
     deallocate(gcc_t,imat_t,bvc_t,gvv_t)
  endif

  return

END SUBROUTINE hamiltonian
