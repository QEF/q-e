!this module contains basic types and routine for dealimg with lanczos chains
 
MODULE lanczos
 USE kinds, ONLY : DP

  SAVE

  TYPE lanczos_chain
     INTEGER :: nc!number of different chains
     INTEGER :: ns!order of tri-diagonal matrix
     REAL(kind=DP), POINTER :: a(:,:)!diagonal part a(ns,nc)
     REAL(kind=DP), POINTER :: b(:,:)!lower and upper parts b(ns,nc)
     LOGICAL, PRIVATE :: l_wfc!if true contains  also wavefunction basis
     INTEGER ::npw_lanczos
     COMPLEX(kind=DP), POINTER :: wfc(:,:,:) !wfc basis (npw,ns,nc) , where (1:npw,nc) contains seeds
     COMPLEX(kind=DP), POINTER :: wfc0(:,:,:) !wfc basis (npw,ns,nc) , where (1:npw,nc) contains seeds the same as above untouch for multiplications with valence sates
     LOGICAL :: l_first!if true, the chain has not yet been created
     LOGICAL :: l_krylov!if true a more generic krylov basis is stored
     COMPLEX(kind=DP), POINTER  :: h(:,:,:) !h matrix on krylov basis (ns,ns,nc)
     INTEGER, POINTER :: ns_chain(:)!number of steps for each chain 
  END  TYPE lanczos_chain

  INTERFACE free_memory
     MODULE PROCEDURE free_lanczos_chain
  END INTERFACE

  INTERFACE initialize_memory
     MODULE PROCEDURE initialize_lanczos_chain
  END  INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_lanczos_chain
  END INTERFACE copy



  CONTAINS

    SUBROUTINE  free_lanczos_chain(lc)
      IMPLICIT NONE
      TYPE(lanczos_chain) :: lc
      if(associated(lc%a)) deallocate(lc%a)
      if(associated(lc%b)) deallocate(lc%b)
      if(associated(lc%wfc)) deallocate(lc%wfc)
      if(associated(lc%wfc0)) deallocate(lc%wfc0)
      if(associated(lc%h)) deallocate(lc%h)
      if(associated(lc%ns_chain)) deallocate(lc%ns_chain)
      nullify(lc%a,lc%b,lc%wfc,lc%wfc0,lc%h,lc%ns_chain)
      lc%l_first=.true.
      RETURN
    END SUBROUTINE free_lanczos_chain

    SUBROUTINE  initialize_lanczos_chain(lc)
      IMPLICIT NONE
      TYPE(lanczos_chain) :: lc
      nullify(lc%a,lc%b,lc%wfc,lc%wfc0,lc%h,lc%ns_chain)
      lc%l_first=.true.
      RETURN
    END SUBROUTINE initialize_lanczos_chain



    SUBROUTINE copy_lanczos_chain(lc_a,lc_b)
!copies A into B
      IMPLICIT NONE
      TYPE(lanczos_chain), INTENT(in) :: lc_a
      TYPE(lanczos_chain), INTENT(out) :: lc_b
      call free_memory(lc_b)
      lc_b%nc=lc_a%nc
      lc_b%ns=lc_a%ns
      lc_b%npw_lanczos=lc_a%npw_lanczos
      allocate(lc_b%a(lc_b%ns,lc_b%nc),lc_b%b(lc_b%ns,lc_b%nc))
      allocate(lc_b%wfc(lc_b%npw_lanczos,lc_b%ns,lc_b%nc))
      allocate(lc_b%wfc0(lc_b%npw_lanczos,lc_b%ns,lc_b%nc))
      lc_b%a(1:lc_b%ns,1:lc_b%nc)=lc_a%a(1:lc_b%ns,1:lc_b%nc)
      lc_b%b(1:lc_b%ns,1:lc_b%nc)=lc_a%b(1:lc_b%ns,1:lc_b%nc)
      lc_b%wfc(1:lc_b%npw_lanczos,1:lc_b%ns,1:lc_b%nc)=lc_a%wfc(1:lc_b%npw_lanczos,1:lc_b%ns,1:lc_b%nc)
      lc_b%wfc0(1:lc_b%npw_lanczos,1:lc_b%ns,1:lc_b%nc)=lc_a%wfc0(1:lc_b%npw_lanczos,1:lc_b%ns,1:lc_b%nc)
      lc_b%l_first=lc_a%l_first
      lc_b%l_krylov=lc_a%l_krylov
      if(lc_b%l_krylov) then
         allocate(lc_b%h(lc_b%ns,lc_b%ns,lc_b%nc))
         lc_b%h(1:lc_b%ns,1:lc_b%ns,1:lc_b%nc)=lc_a%h(1:lc_b%ns,1:lc_b%ns,1:lc_b%nc)
      else
         nullify(lc_b%h)
      endif
      allocate(lc_b%ns_chain(lc_b%nc))
      lc_b%ns_chain(1:lc_b%nc)=lc_a%ns_chain(1:lc_b%nc)

      RETURN
    END SUBROUTINE copy_lanczos_chain

    



    SUBROUTINE create_lanczos_chain(lc,op,nc,ns,seed,l_wfc,param, ic_update)
!this subroutine calculates all the chains AT ONCE for operator ip
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum, mp_barrier
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout 
      USE wannier_gw, ONLY : l_verbose

      IMPLICIT NONE
      TYPE(lanczos_chain) :: lc!lanczos chain to be created
      EXTERNAL :: op!operator for creating the lanczos chains
      INTEGER  :: nc !number of chains
      INTEGER  :: ns!number of steps in each chain
      COMPLEX(kind=DP) :: seed(npw,nc)!seeds for starting the chains
      LOGICAL :: l_wfc!if true allocate and stores the chain vectors
      REAL(kind=DP) :: param(nc)!energies or parameters for op routine
      INTEGER :: ic_update!if 0 create new chains if /= 0 updated ic_update chain

      COMPLEX(kind=DP), ALLOCATABLE :: psi_1(:,:),psi_2(:,:),psi_3(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: u_0(:,:),u_1(:,:)
      REAL(kind=DP), ALLOCATABLE :: alpha(:),beta(:),gamma(:), n_1(:)
      REAL(kind=DP), ALLOCATABLE :: c(:)
      INTEGER :: it,ic,ig
      REAL(kind=DP)::sca
      INTEGER :: ii,jj
      INTEGER :: ic_first,ic_last


      if(ic_update==0) then
         allocate(lc%a(ns,nc),lc%b(ns,nc))
         lc%a=0.d0
         lc%b=0.d0
         if(l_wfc) then
            lc%l_wfc=.true.
            allocate(lc%wfc(npw,ns,nc))
            allocate(lc%wfc0(npw,ns,nc))
         endif

         lc%nc=nc
         lc%ns=ns
         lc%npw_lanczos=npw

         ic_first=1
         ic_last=lc%nc
      else
         lc%a(:,ic_update)=0.d0
         lc%b(:,ic_update)=0.d0
         ic_first=ic_update
         ic_last=ic_update
         lc%npw_lanczos=npw
      endif
      lc%l_krylov=.false.

      allocate(psi_1(npw,nc),u_0(npw,nc),psi_2(npw,nc),psi_3(npw,nc),u_1(npw,nc))
      allocate(n_1(nc),alpha(nc),beta(nc),gamma(nc))
!set first state and NORMALISE IT

      psi_1(1:npw,ic_first:ic_last)=seed(1:npw,ic_first:ic_last)

      do ic=ic_first,ic_last
         sca=0.d0
         do ig=1,npw
            sca=sca+2.d0*dble(conjg(psi_1(ig,ic))*psi_1(ig,ic))
         enddo
         if(gstart==2) sca=sca-dble(conjg(psi_1(1,ic))*psi_1(1,ic))
         call mp_sum(sca,world_comm)
         psi_1(1:npw,ic)=(1.d0/sqrt(sca))*psi_1(1:npw,ic)
      enddo

      call op( npw, psi_1(1,ic_first), u_0(1,ic_first),param,1,ic_last-ic_first+1 )

      !calculate n_1
      n_1(ic_first:ic_last)=0.d0
      do ic=ic_first,ic_last
         do ig=1,npw
            n_1(ic)=n_1(ic)+2.d0*dble(conjg(u_0(ig,ic))*u_0(ig,ic))
         enddo
         if(gstart==2) n_1(ic)=n_1(ic)-dble(conjg(u_0(1,ic))*u_0(1,ic))
      enddo
      call mp_sum(n_1, world_comm)
      n_1(ic_first:ic_last)=dsqrt(n_1(ic_first:ic_last))
!calculate alpha     
      alpha(ic_first:ic_last)=0.d0
     do ic=ic_first,ic_last
        do ig=1,npw
           alpha(ic)=alpha(ic)+2.d0*dble(conjg(psi_1(ig,ic))*u_0(ig,ic))
        enddo
        if(gstart==2) alpha(ic)=alpha(ic)-dble(conjg(psi_1(1,ic))*u_0(1,ic))
     enddo
     call mp_sum(alpha, world_comm)
     alpha(ic_first:ic_last)=alpha(ic_first:ic_last)/n_1(ic_first:ic_last)

!calculate psi_2 and beta                                                                                                      
     do ic=ic_first,ic_last
        psi_2(1:npw,ic)=u_0(1:npw,ic)/n_1(ic)-alpha(ic)*psi_1(1:npw,ic)
     enddo
     beta=0.d0
     do ic=ic_first,ic_last
        do ig=1,npw
           beta(ic)=beta(ic)+2.d0*dble(conjg(psi_2(ig,ic))*psi_2(ig,ic))
        enddo
        if(gstart==2) beta(ic)=beta(ic)-dble(conjg(psi_2(1,ic))*psi_2(1,ic))
     enddo
     call mp_sum(beta(ic_first:ic_last), world_comm)
     beta(ic_first:ic_last)=dsqrt(beta(ic_first:ic_last))
     do ic=ic_first,ic_last
        psi_2(1:npw,ic)=psi_2(1:npw,ic)/beta(ic)
     enddo

!calculate d
     do ic=ic_first,ic_last
        do ig=1,npw
           lc%a(1,ic)=lc%a(1,ic)+2.d0*dble(conjg(psi_1(ig,ic))*u_0(ig,ic))
        enddo
        if(gstart==2) lc%a(1,ic)=lc%a(1,ic)-dble(conjg(psi_1(1,ic))*u_0(1,ic))
     enddo
     call mp_sum(lc%a(1,ic_first:ic_last), world_comm)

!calculate f

     do ic=ic_first,ic_last
        do ig=1,npw
           lc%b(1,ic)=lc%b(1,ic)+2.d0*dble(conjg(psi_2(ig,ic))*u_0(ig,ic))
        enddo
        if(gstart==2) lc%b(1,ic)=lc%b(1,ic)-dble(conjg(psi_2(1,ic))*u_0(1,ic))
     enddo
     call mp_sum(lc%b(1,ic_first:ic_last), world_comm)

     do ic=ic_first,ic_last
        lc%wfc(1:npw,1,ic)=psi_1(1:npw,ic)
        lc%wfc(1:npw,2,ic)=psi_2(1:npw,ic)
     enddo

     do it=2,ns
        if(l_verbose) write(stdout,*) 'LANCZOS STEP: ', it
        call op( npw, psi_2(1,ic_first), u_1(1,ic_first),param,1,ic_last-ic_first+1 )
!calculate n_1

        n_1(ic_first:ic_last)=0.d0
        do ic=ic_first,ic_last
           do ig=1,npw
              n_1(ic)=n_1(ic)+2.d0*dble(conjg(u_1(ig,ic))*u_1(ig,ic))
           enddo
           if(gstart==2) n_1(ic)=n_1(ic)-dble(conjg(u_1(1,ic))*u_1(1,ic))
        enddo
        call mp_sum(n_1(ic_first:ic_last), world_comm)
        n_1(ic_first:ic_last)=dsqrt(n_1(ic_first:ic_last))

!calculate alpha

        alpha(ic_first:ic_last)=0.d0
        do ic=ic_first,ic_last
           do ig=1,npw
              alpha(ic)=alpha(ic)+2.d0*dble(conjg(psi_1(ig,ic))*u_1(ig,ic))
           enddo
           if(gstart==2) alpha(ic)=alpha(ic)-dble(conjg(psi_1(1,ic))*u_1(1,ic))
        enddo
        call mp_sum(alpha(ic_first:ic_last), world_comm)
        alpha(ic_first:ic_last)=alpha(ic_first:ic_last)/n_1(ic_first:ic_last)

!calculate beta                 
      
        beta(ic_first:ic_last)=0.d0
        do ic=ic_first,ic_last
           do ig=1,npw
              beta(ic)=beta(ic)+2.d0*dble(conjg(psi_2(ig,ic))*u_1(ig,ic))
           enddo
           if(gstart==2) beta(ic)=beta(ic)-dble(conjg(psi_2(1,ic))*u_1(1,ic))
        enddo
        call mp_sum(beta, world_comm)
        beta(ic_first:ic_last)=beta(ic_first:ic_last)/n_1(ic_first:ic_last)

!calculate psi_3 and gamma                                                                                                     
        do ic=ic_first,ic_last
           psi_3(1:npw,ic)=u_1(1:npw,ic)/n_1(ic)-alpha(ic)*psi_1(1:npw,ic)-beta(ic)*psi_2(1:npw,ic)
           if(l_verbose) write(stdout,*) 'LANCZOS: ',n_1(ic),alpha(ic),beta(ic)
        enddo
        gamma(ic_first:ic_last)=0.d0
        do ic=ic_first,ic_last
           do ig=1,npw
              gamma(ic)=gamma(ic)+2.d0*dble(conjg(psi_3(ig,ic))*psi_3(ig,ic))
           enddo
           if(gstart==2) gamma(ic)=gamma(ic)-dble(conjg(psi_3(1,ic))*psi_3(1,ic))
        enddo
        call mp_sum(gamma(ic_first:ic_last), world_comm)
        gamma(ic_first:ic_last)=dsqrt(gamma(ic_first:ic_last))
        do ic=ic_first,ic_last
           psi_3(1:npw,ic)=psi_3(1:npw,ic)/gamma(ic)
        enddo

!calculate diagonal

        do ic=ic_first,ic_last
           do ig=1,npw
              lc%a(it,ic)=lc%a(it,ic)+2.d0*dble(conjg(psi_2(ig,ic))*u_1(ig,ic))
           enddo
           if(gstart==2) lc%a(it,ic)=lc%a(it,ic)-dble(conjg(psi_2(1,ic))*u_1(1,ic))
        enddo
        call mp_sum(lc%a(it,ic_first:ic_last), world_comm)

!calculate upper/lower diagonal 
        do ic=ic_first,ic_last
           do ig=1,npw
              lc%b(it,ic)=lc%b(it,ic)+2.d0*dble(conjg(psi_3(ig,ic))*u_1(ig,ic))
           enddo
           if(gstart==2) lc%b(it,ic)=lc%b(it,ic)-dble(conjg(psi_3(1,ic))*u_1(1,ic))
        enddo
        call mp_sum(lc%b(it,ic_first:ic_last), world_comm)

        do ic=ic_first,ic_last
           lc%wfc(1:npw,it,ic)=psi_2(1:npw,ic)
        enddo

!update arrays                                                                                                                 
        psi_1(1:npw,ic_first:ic_last)=psi_2(1:npw,ic_first:ic_last)
        psi_2(1:npw,ic_first:ic_last)=psi_3(1:npw,ic_first:ic_last)
        u_0(1:npw,ic_first:ic_last)=u_1(1:npw,ic_first:ic_last)







     enddo
     do ii=1,ns
        do jj=1,ns
           sca=0.d0
           do ig=1,npw
              sca=sca+2.d0*dble(conjg(lc%wfc(ig,ii,1))*lc%wfc(ig,jj,1))
           enddo
           if(gstart==2) sca=sca-dble(conjg(lc%wfc(1,ii,1))*lc%wfc(1,jj,1))
           call mp_sum(sca,world_comm)
           if(l_verbose) write(stdout,*) 'ORTHONORMALITY LANCZOS', ii, jj, sca
        enddo
     enddo
     if(lc%l_wfc) then
        lc%wfc0(1:npw,1:lc%ns,1:lc%nc)=lc%wfc(1:npw,1:lc%ns,1:lc%nc)
        
     endif
     allocate(lc%ns_chain(lc%nc))
     lc%ns_chain=lc%ns

     deallocate(psi_1,u_0,psi_2,psi_3,u_1)
     deallocate(n_1,alpha,beta,gamma)
      RETURN
    END SUBROUTINE create_lanczos_chain



    SUBROUTINE product_chain(lc,v_states,ic_update)
!multiply on real space grid the lanczos basis with the vstate vectors
!this is used for calculating W
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout
      USE fft_base,             ONLY : dfftp, dffts
      USE fft_interfaces,       ONLY : fwfft, invfft
      USE wavefunctions, ONLY : psic



      IMPLICIT NONE

      TYPE(lanczos_chain), INTENT(inout) :: lc !data to be modified
      REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,lc%nc)!the states to be multiplied with, namely valence states
      INTEGER :: ic_update!if 0 create new chains if /= 0 updated ic_update chain 
      
      REAL(kind=DP), ALLOCATABLE :: psi_tmp(:,:)
      INTEGER :: it,ic
      INTEGER :: ic_first,ic_last

      allocate(psi_tmp(dffts%nnr,2))

      if(ic_update==0) then
         ic_first=1
         ic_last=lc%nc
      else
          ic_first=ic_update
          ic_last=ic_update
      endif


      do ic=ic_first,ic_last
         do it=1,lc%ns,2
            call pc_operator(lc%wfc(1,it,ic),1,.false.)
            if(it/=lc%ns)  call pc_operator(lc%wfc(1,it+1,ic),1,.false.)
!put lanczos basis on real space
            psic(:)=(0.d0,0.d0)
            if(it/=lc%ns) then
               psic(dffts%nl(1:npw))  = lc%wfc(1:npw,it,ic)+(0.d0,1.d0)*lc%wfc(1:npw,it+1,ic)
               psic(dffts%nlm(1:npw)) = CONJG( lc%wfc(1:npw,it,ic) )+(0.d0,1.d0)*conjg(lc%wfc(1:npw,it+1,ic))
            else
               psic(dffts%nl(1:npw))  = lc%wfc(1:npw,it,ic)
               psic(dffts%nlm(1:npw)) = CONJG( lc%wfc(1:npw,it,ic) )
            endif
            CALL invfft ('Wave', psic, dffts)
            if(it/=lc%ns) then
               psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
               psi_tmp(1:dffts%nnr,2)= aimag(psic(1:dffts%nnr))
            else
               psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
            endif
!!product with valence states                     
            if(it/=lc%ns) then
               psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,ic)
               psi_tmp(1:dffts%nnr,2)=psi_tmp(1:dffts%nnr,2)*v_states(1:dffts%nnr,ic)
            else
               psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,ic)
            endif
!bring back to G space

            if(it/=lc%ns) then
               psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),psi_tmp(1:dffts%nnr,2))
            else
               psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),0.d0)
            endif
            CALL fwfft ('Wave', psic, dffts)
            if(it/=lc%ns) then
                lc%wfc(1:npw,it,ic)=0.5d0*(psic(dffts%nl(1:npw))+conjg(psic(dffts%nlm(1:npw))))
                lc%wfc(1:npw,it+1,ic)=(0.d0,-0.5d0)*(psic(dffts%nl(1:npw))-conjg(psic(dffts%nlm(1:npw))))
                if(gstart==2) lc%wfc(1,it,ic)=dble(lc%wfc(1,it,ic))
                if(gstart==2) lc%wfc(1,it+1,ic)=dble(lc%wfc(1,it+1,ic))
             else
                lc%wfc(1:npw,it,ic)=psic(dffts%nl(1:npw))
                if(gstart==2) lc%wfc(1,it,ic)=dble(lc%wfc(1,it,ic))
             endif

          enddo
       enddo

      deallocate(psi_tmp)
      RETURN

    END SUBROUTINE product_chain



    SUBROUTINE solve_lanczos(lc,wfc_in,freq,wfc_out_re,wfc_out_im, l_extrapolate,wfc_in_mult0,wfc_in_mult,op, param)
      
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout



      IMPLICIT NONE
      TYPE(lanczos_chain), INTENT(in) :: lc!date for the lanczos chain previously calcaulted
      COMPLEX(kind=DP), INTENT(in) :: wfc_in(npw,lc%nc)!wavefunction (real) on which the inverted operator should be applied
      COMPLEX(kind=DP),INTENT(in) :: freq!complex number to be added to the diagonal
      COMPLEX(kind=DP), INTENT(out) :: wfc_out_re(npw,lc%nc)!solution (real part)
      COMPLEX(kind=DP), INTENT(out) :: wfc_out_im(npw,lc%nc)!solution (imaginary part)  
      LOGICAL :: l_extrapolate!if true estrapolate solution
      COMPLEX(kind=DP), INTENT(in) :: wfc_in_mult0(npw,lc%nc)!wavefunction on input evetually multiplied with \psi_v 
      COMPLEX(kind=DP), INTENT(in) :: wfc_in_mult(npw,lc%nc)!wavefunction on input evetually multiplied with \psi_v \_psi_v
      EXTERNAL :: op!operator for creating the lanczos chains
      REAL(kind=DP) :: param(lc%nc)!energies or parameters for op routine


      COMPLEX(kind=DP), ALLOCATABLE :: dl(:),du(:),d(:),t(:)
      REAL(kind=DP), ALLOCATABLE  :: t_re(:),t_im(:)
      INTEGER :: ic,ig,it
      INTEGER :: info
      REAL(kind=DP) :: sca,sca1,ene
      COMPLEX(kind=DP), ALLOCATABLE :: wfc_remainder(:),h_wfc_remainder(:)

      allocate(dl(lc%ns-1),du(lc%ns-1),d(lc%ns),t(lc%ns),t_re(lc%ns),t_im(lc%ns))
      allocate(wfc_remainder(npw),h_wfc_remainder(npw))


      do ic=1,lc%nc
         

         dl(1:lc%ns-1)=cmplx(lc%b(1:lc%ns-1,ic),0.d0)
         du(1:lc%ns-1)=cmplx(lc%b(1:lc%ns-1,ic),0.d0)
         d(1:lc%ns)=cmplx(lc%a(1:lc%ns,ic),0.d0)-freq
         sca1=0.d0
         do ig=1,npw
            sca1=sca1+2.d0*dble(conjg(wfc_in(ig,ic)*wfc_in(ig,ic)))
         enddo
         if(gstart==2) sca1=sca1-dble(conjg(wfc_in(1,ic)*wfc_in(1,ic)))
         call mp_sum(sca1,world_comm)

!calculate t vectors projecting wfc_in
         CALL ZGEMM('C','N',lc%ns,1,npw,(1.d0,0.d0),lc%wfc(1,1,ic),npw,wfc_in(1,ic),npw,(0.d0,0.d0),t,lc%ns)
         t(1:lc%ns)=t(1:lc%ns)+conjg(t(1:lc%ns))
         if(gstart==2) then
            t(1:lc%ns)=t(1:lc%ns)-conjg(lc%wfc(1,1:lc%ns,ic))*wfc_in(1,ic)
         endif
         call mp_sum(t,world_comm)
         wfc_remainder(1:npw)=wfc_in_mult0(1:npw,ic)
         CALL ZGEMM('N','N',npw,1,lc%ns,(-1.d0,0.d0),lc%wfc0(1,1,ic),npw,t,lc%ns,(1.d0,0.d0),wfc_remainder,npw)
!calculate remainder energy
         call op( npw, wfc_remainder, h_wfc_remainder,param,1,1 )
         ene=0.d0
         sca=0.d0
         do ig=1,npw
            ene=ene+2.d0*dble(conjg(wfc_remainder(ig))*h_wfc_remainder(ig))
            sca=sca+2.d0*dble(conjg(wfc_remainder(ig))*wfc_remainder(ig))
         enddo
         if(gstart==2) then
            ene=ene-dble(conjg(wfc_remainder(1))*h_wfc_remainder(1))
            sca=sca-dble(conjg(wfc_remainder(1))*wfc_remainder(1))
         endif
         call mp_sum(ene,world_comm)
         call mp_sum(sca,world_comm)
         ene=ene/sca

         wfc_remainder(1:npw)=wfc_in_mult(1:npw,ic)
         CALL ZGEMM('N','N',npw,1,lc%ns,(-1.d0,0.d0),lc%wfc(1,1,ic),npw,t,lc%ns,(1.d0,0.d0),wfc_remainder,npw)
         

         sca=0.d0
         do it=1,lc%ns
            sca=sca+dble(conjg(t(it))*t(it))
         enddo
       
        ! write(stdout,*) 'PROIEZIONE', ic, sca
       
!!call lapack 

        call ZGTSV(lc%ns,1,dl,d,du,t,lc%ns,info)
        if(info /= 0) then
           write(stdout,*) 'DGTSV info:', info
           FLUSH(stdout)
           stop
        endif
     
        t_re(1:lc%ns)=dble(t(1:lc%ns))
        t_im(1:lc%ns)=aimag(t(1:lc%ns))

        CALL DGEMM('N','N',2*npw,1,lc%ns,1.d0,lc%wfc(1,1,ic),2*npw,t_re,lc%ns,0.d0,wfc_out_re(1,ic),2*npw)
        CALL DGEMM('N','N',2*npw,1,lc%ns,1.d0,lc%wfc(1,1,ic),2*npw,t_im,lc%ns,0.d0,wfc_out_im(1,ic),2*npw)


!in case add extrapolation
        if(l_extrapolate) then
           !wfc_out_re(1:npw,ic)=wfc_out_re(1:npw,ic)+dble(1.d0/(lc%a(lc%ns,ic)-freq))*wfc_remainder(1:npw)
           !wfc_out_im(1:npw,ic)=wfc_out_im(1:npw,ic)+aimag(1.d0/(lc%a(lc%ns,ic)-freq))*wfc_remainder(1:npw)
           wfc_out_re(1:npw,ic)=wfc_out_re(1:npw,ic)+dble(1.d0/(ene-freq))*wfc_remainder(1:npw)
           wfc_out_im(1:npw,ic)=wfc_out_im(1:npw,ic)+aimag(1.d0/(ene-freq))*wfc_remainder(1:npw)
        endif
     enddo


     deallocate(dl,du,d,t,t_re,t_im)
     deallocate(wfc_remainder,h_wfc_remainder)

     RETURN

    END SUBROUTINE solve_lanczos

    SUBROUTINE norms_lanczos(lc,wfc_in, norms)
!calculate the norms of wfc_in projected on the lanczos bais

      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout



      IMPLICIT NONE


      TYPE(lanczos_chain), INTENT(in) :: lc!date for the lanczos chain previously calcaulted
      COMPLEX(kind=DP), INTENT(in) :: wfc_in(npw,lc%nc)!wavefunction (real) on which the inverted operator should be applied
      REAL(kind=DP), INTENT(out) :: norms(lc%nc)
      
      
      INTEGER :: ic,it 
      COMPLEX(kind=DP), ALLOCATABLE :: t(:)
      allocate(t(lc%ns))
      do ic=1,lc%nc
         CALL ZGEMM('C','N',lc%ns,1,npw,(1.d0,0.d0),lc%wfc(1,1,ic),npw,wfc_in(1,ic),npw,(0.d0,0.d0),t,lc%ns)
         t(1:lc%ns)=t(1:lc%ns)+conjg(t(1:lc%ns))
         if(gstart==2) then
            t(1:lc%ns)=t(1:lc%ns)-conjg(lc%wfc(1,1:lc%ns,ic))*wfc_in(1,ic)
         endif
         call mp_sum(t,world_comm)
         norms(ic)=0.d0
         do it=1,lc%ns
            norms(ic)=norms(ic)+dble(conjg(t(it))*t(it))
         enddo
      enddo

      deallocate(t)
    END SUBROUTINE norms_lanczos

    SUBROUTINE product_wfc(v_states,n_states, wfc_in,wfc_out)
!multiply on real space grid the wave-functrionswith the vstate vectors 
!this is used for calculating W 

      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout
      USE fft_base,             ONLY : dfftp, dffts
      USE fft_interfaces,       ONLY : fwfft, invfft
      USE wavefunctions, ONLY : psic



      IMPLICIT NONE


      REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,n_states)!the states to be multiplied with, namely valence states                  
      INTEGER :: n_states
      COMPLEX(kind=DP), INTENT(in) :: wfc_in(npw,n_states)!input states
      COMPLEX(kind=DP), INTENT(out) :: wfc_out(npw,n_states)!ouput states

      REAL(kind=DP), ALLOCATABLE :: psi_tmp(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: psi_g_tmp(:,:)
      INTEGER :: ic
    

      allocate(psi_tmp(dffts%nnr,2))
      allocate(psi_g_tmp(npw,2))

    

      do ic=1,n_states,2
         psi_g_tmp(1:npw,1)=wfc_in(1:npw,ic)
         call pc_operator(psi_g_tmp(1,1),1,.false.)
         if(ic/=n_states)  then
            psi_g_tmp(1:npw,2)=wfc_in(1:npw,ic+1)
            call pc_operator(psi_g_tmp(1,2),1,.false.)
         endif
!put lanczos basis on real space                                                                                                                                                                                                                                                                                                                                                    
            psic(:)=(0.d0,0.d0)
            if(ic/=n_states) then
               psic(dffts%nl(1:npw))  = psi_g_tmp(1:npw,1)+(0.d0,1.d0)*psi_g_tmp(1:npw,2)
               psic(dffts%nlm(1:npw)) = CONJG( psi_g_tmp(1:npw,1) )+(0.d0,1.d0)*conjg(psi_g_tmp(1:npw,2))
            else
               psic(dffts%nl(1:npw))  = psi_g_tmp(1:npw,1)
               psic(dffts%nlm(1:npw)) = CONJG( psi_g_tmp(1:npw,2 ))
            endif
            CALL invfft ('Wave', psic, dffts)
            if(ic/=n_states) then
               psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
               psi_tmp(1:dffts%nnr,2)= aimag(psic(1:dffts%nnr))
            else
               psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
            endif
!!product with valence states                                                                                                                                                                                                                                                                                                                                                       
            if(ic/=n_states) then
               psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,ic)
               psi_tmp(1:dffts%nnr,2)=psi_tmp(1:dffts%nnr,2)*v_states(1:dffts%nnr,ic+1)
            else
               psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,ic)
            endif
!bring back to G space                                                                                                                                                                                                                                                                                                                                                              

            if(ic/=n_states) then
               psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),psi_tmp(1:dffts%nnr,2))
            else
               psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),0.d0)
            endif
            CALL fwfft ('Wave', psic, dffts)
            if(ic/=n_states) then
                wfc_out(1:npw,ic)=0.5d0*(psic(dffts%nl(1:npw))+conjg(psic(dffts%nlm(1:npw))))
                wfc_out(1:npw,ic+1)=(0.d0,-0.5d0)*(psic(dffts%nl(1:npw))-conjg(psic(dffts%nlm(1:npw))))
                if(gstart==2) wfc_out(1,ic)=dble(wfc_out(1,ic))
                if(gstart==2) wfc_out(1,ic+1)=dble(wfc_out(1,ic+1))
             else
                wfc_out(1:npw,ic)=psic(dffts%nl(1:npw))
                if(gstart==2) wfc_out(1,ic)=dble(wfc_out(1,ic))
             endif

         
       enddo

      deallocate(psi_tmp)
      deallocate(psi_g_tmp)
      RETURN

    END SUBROUTINE product_wfc
    


    SUBROUTINE create_krylov(lc,op,nc,ns,seed,param,ispin)!,ks_wfcs)
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum, mp_barrier
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout
      USE wannier_gw, ONLY : l_verbose
      USE lsda_mod, ONLY : nspin

      IMPLICIT NONE
      TYPE(lanczos_chain) :: lc!lanczos chain to be created
      EXTERNAL :: op!operator for creating the lanczos chains
      COMPLEX(kind=DP) :: seed(npw,nc)!seeds for updating the chains
      REAL(kind=DP) :: param(nc)!energies or parameters for op routine
      INTEGER, INTENT(in) :: ns!number of basis terms 
      INTEGER, INTENT(in) :: nc!number of chains
      INTEGER, INTENT(in) :: ispin!spin index
     ! COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd, nspin)!KS wavefunctions

      INTEGER :: ns_new
      INTEGER :: ii,jj,is,ig,iv
      COMPLEX(kind=DP), ALLOCATABLE :: h_tmp(:,:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: wfc_tmp(:),wfc_tmp_all(:,:,:),h_wfc_tmp(:)
      REAL(kind=DP), ALLOCATABLE :: omat(:,:), param_tmp(:)
      REAL(kind=DP) ::sca

      lc%l_krylov=.true.
      lc%ns=ns
      lc%nc=nc
      lc%npw_lanczos=npw
      allocate(lc%ns_chain(lc%nc))
      lc%ns_chain=lc%ns
      allocate(lc%h(ns,ns,lc%nc))
      allocate(lc%wfc0(npw,ns,lc%nc))
      allocate(lc%wfc(npw,ns,lc%nc))
      allocate(lc%a(ns,lc%nc),lc%b(ns,lc%nc))
        
      allocate(wfc_tmp(npw),h_wfc_tmp(npw))

      do iv=1,lc%nc
         wfc_tmp(1:npw)=seed(1:npw,iv)
         ns_new=0
      !loop on steps new 
         do is=1,ns

            if(ns_new>0) then
               allocate(omat(ns_new,1))
         !project out old basis vectors 
               CALL DGEMM( 'T','N',ns_new,1,2*npw,2.d0,lc%wfc0(1,1,iv),2*npw,&
                    &wfc_tmp,2*npw,0.d0,omat,ns_new)

               if(gstart==2) then
                  do ii=1,ns_new
                     omat(ii,1)=omat(ii,1)-dble(conjg(lc%wfc0(1,ii,iv))*wfc_tmp(1))
                  enddo
               endif
               call mp_sum(omat,world_comm)

               CALL DGEMM('N','N',2*npw,1,ns_new,-1.d0,lc%wfc0(1,1,iv),2*npw,omat,ns_new,1.d0,wfc_tmp,2*npw)
               deallocate(omat)
            endif

         !normalize
            sca=0.d0
            do ig=1,npw
               sca=sca+2.d0*dble(conjg(wfc_tmp(ig))*wfc_tmp(ig))
            enddo
            if(gstart==2) sca=sca-dble(conjg(wfc_tmp(1))*wfc_tmp(1))
            call mp_sum(sca,world_comm)
            wfc_tmp(1:npw)=(1.d0/sqrt(sca))*wfc_tmp(1:npw)
         
            allocate(omat(ns_new+1,1))
          !add to basis
            lc%wfc0(1:npw,ns_new+1,iv)=wfc_tmp(1:npw)
         !apply Op
            allocate(param_tmp(lc%nc))
            param_tmp(1)=param(iv)
            call op( npw,wfc_tmp,h_wfc_tmp,param_tmp,ispin,1)!,ks_wfcs )
            deallocate(param_tmp)
                    
            CALL DGEMM('T','N',ns_new+1,1,2*npw,2.d0,lc%wfc0(1,1,iv),2*npw,h_wfc_tmp,2*npw,0.d0,omat,ns_new+1)


            if(gstart==2) then
               do ii=1,ns_new+1
                  omat(ii,1)=omat(ii,1)-dble(conjg(lc%wfc0(1,ii,iv))*h_wfc_tmp(1))
               enddo
            endif
            call mp_sum(omat,world_comm)
         !update h  
            lc%h(1:ns_new+1,ns_new+1,iv)=omat(1:ns_new+1,1)
            lc%h(ns_new+1,1:ns_new+1,iv)=omat(1:ns_new+1,1)



            deallocate(omat)
        !update gradient state  
            wfc_tmp(1:npw)=h_wfc_tmp(1:npw)
            ns_new=ns_new+1
         enddo
      enddo
      lc%wfc(1:npw,1:lc%ns,1:lc%nc)=lc%wfc0(1:npw,1:lc%ns,1:lc%nc)
     
      deallocate(wfc_tmp,h_wfc_tmp)

    END SUBROUTINE create_krylov



    SUBROUTINE update_krylov(lc,op,seed,param, ns_update,iv,ispin,ks_wfcs)
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum, mp_barrier
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout 
      USE wannier_gw, ONLY : l_verbose
      USE lsda_mod, ONLY : nspin

      IMPLICIT NONE
      TYPE(lanczos_chain) :: lc!lanczos chain to be increased
      EXTERNAL :: op!operator for creating the lanczos chains
      COMPLEX(kind=DP) :: seed(npw,lc%nc)!seeds for updating the chains
      REAL(kind=DP) :: param(lc%nc)!energies or parameters for op routine
      INTEGER, INTENT(in) :: ns_update!number of basis terms to be added
      INTEGER, INTENT(in) :: iv!chain to be updated
      INTEGER, INTENT(in) :: ispin!spin index  
      COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npw,nbnd, nspin)!KS wavefunctions      

      INTEGER :: ns_old,ns,ns_new,ns_recent
      INTEGER :: ii,jj,is,ig
      COMPLEX(kind=DP), ALLOCATABLE :: h_tmp(:,:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: wfc_tmp(:),wfc_tmp_all(:,:,:),h_wfc_tmp(:)
      REAL(kind=DP), ALLOCATABLE :: omat(:,:), param_tmp(:)
      REAL(kind=DP) ::sca


      if(.not.lc%l_krylov) then
         ns_old=lc%ns
         ns=lc%ns+ns_update
         lc%ns_chain(iv)=ns
         lc%l_krylov=.true.
         allocate(lc%h(ns,ns,lc%nc))
         lc%h(1:ns_old,1:ns_old,1:lc%nc)=0.d0
         do ii=1,ns_old
            lc%h(ii,ii,1:lc%nc)=lc%a(ii,1:lc%nc)
         enddo
         do ii=1,ns_old-1
            lc%h(ii,ii+1,1:lc%nc)=lc%b(ii,1:lc%nc)
         enddo
         do ii=2,ns_old
            lc%h(ii,ii-1,1:lc%nc)=lc%b(ii,1:lc%nc)
         enddo
      else
         ns_old=lc%ns_chain(iv)
         ns=lc%ns_chain(iv)+ns_update
         lc%ns_chain(iv)=ns
      endif
      if(ns>lc%ns) then
         allocate(h_tmp(lc%ns,lc%ns,lc%nc))
         h_tmp(1:lc%ns,1:lc%ns,1:lc%nc)=lc%h(1:lc%ns,1:lc%ns,1:lc%nc)
         deallocate(lc%h)
         allocate(lc%h(ns,ns,lc%nc))
         lc%h(1:lc%ns,1:lc%ns,1:lc%nc)=h_tmp(1:lc%ns,1:lc%ns,1:lc%nc)
         deallocate(h_tmp)
         
         allocate(wfc_tmp_all(npw,lc%ns,lc%nc))
         wfc_tmp_all(1:npw,1:lc%ns,1:lc%nc)=lc%wfc0(1:npw,1:lc%ns,1:lc%nc)
         deallocate(lc%wfc0)
         allocate(lc%wfc0(npw,ns,lc%nc))
         lc%wfc0(1:npw,1:lc%ns,1:lc%nc)=wfc_tmp_all(1:npw,1:lc%ns,1:lc%nc)
         wfc_tmp_all(1:npw,1:lc%ns,1:lc%nc)=lc%wfc(1:npw,1:lc%ns,1:lc%nc)
         deallocate(lc%wfc)
         allocate(lc%wfc(npw,ns,lc%nc))
         lc%wfc(1:npw,1:lc%ns,1:lc%nc)=wfc_tmp_all(1:npw,1:lc%ns,1:lc%nc)
         deallocate(wfc_tmp_all)
         lc%ns=ns
         deallocate(lc%a)
         deallocate(lc%b)
         allocate(lc%a(lc%ns,lc%nc),lc%b(lc%ns,lc%nc))
         lc%a=0.d0
         lc%b=0.d0
      endif

      allocate(wfc_tmp(npw),h_wfc_tmp(npw))
      wfc_tmp(1:npw)=seed(1:npw,iv)
      ns_new=ns_old
      !loop on steps new
      do is=1,ns_update

         allocate(omat(ns_new,1))
         !project out old basis vectors
         CALL DGEMM( 'T','N',ns_new,1,2*npw,2.d0,lc%wfc0(1,1,iv),2*npw,&
              &wfc_tmp,2*npw,0.d0,omat,ns_new)
        
         if(gstart==2) then
            do ii=1,ns_new
               omat(ii,1)=omat(ii,1)-dble(conjg(lc%wfc0(1,ii,iv))*wfc_tmp(1))
            enddo
         endif
         call mp_sum(omat,world_comm)
         
         CALL DGEMM('N','N',2*npw,1,ns_new,-1.d0,lc%wfc0(1,1,iv),2*npw,omat,ns_new,1.d0,wfc_tmp,2*npw)
         !normalize
         sca=0.d0
         do ig=1,npw
            sca=sca+2.d0*dble(conjg(wfc_tmp(ig))*wfc_tmp(ig))
         enddo
         if(gstart==2) sca=sca-dble(conjg(wfc_tmp(1))*wfc_tmp(1))
         call mp_sum(sca,world_comm)
         wfc_tmp(1:npw)=(1.d0/sqrt(sca))*wfc_tmp(1:npw)
         
         deallocate(omat)
         allocate(omat(ns_new+1,1))
          !add to basis
         lc%wfc0(1:npw,ns_new+1,iv)=wfc_tmp(1:npw)
         !apply Op
         allocate(param_tmp(lc%nc))
         param_tmp(1)=param(iv)
         call op( npw,wfc_tmp,h_wfc_tmp,param_tmp,ispin,1 ,ks_wfcs)
         deallocate(param_tmp)

         ns_recent=ns_new+1
         CALL DGEMM('T','N',ns_recent,1,2*npw,2.d0,lc%wfc0(1,1,iv),2*npw,h_wfc_tmp,2*npw,0.d0,omat,ns_recent)
  
        
         if(gstart==2) then
            do ii=1,ns_new+1
               omat(ii,1)=omat(ii,1)-dble(conjg(lc%wfc0(1,ii,iv))*h_wfc_tmp(1))
            enddo
         endif
         call mp_sum(omat,world_comm)
         !update h
         lc%h(1:ns_new+1,ns_new+1,iv)=omat(1:ns_new+1,1)
         lc%h(ns_new+1,1:ns_new+1,iv)=omat(1:ns_new+1,1)
            
               
            
         deallocate(omat)
        !update gradient state
         wfc_tmp(1:npw)=h_wfc_tmp(1:npw)
         ns_new=ns_new+1
      enddo

      
      deallocate(wfc_tmp,h_wfc_tmp)

    END SUBROUTINE update_krylov


    SUBROUTINE solve_krylov(lc,wfc_in,freq,wfc_out_re,wfc_out_im)
!solve generic hamiltonian 
     
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum,mp_bcast
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout, ionode,ionode_id



      IMPLICIT NONE
      TYPE(lanczos_chain), INTENT(in) :: lc!date for the lanczos chain previously calcaulted
      COMPLEX(kind=DP), INTENT(in) :: wfc_in(npw,lc%nc)!wavefunction (real) on which the inverted operator should be applied
      COMPLEX(kind=DP),INTENT(in) :: freq!complex number to be added to the diagonal
      COMPLEX(kind=DP), INTENT(out) :: wfc_out_re(npw,lc%nc)!solution (real part)
      COMPLEX(kind=DP), INTENT(out) :: wfc_out_im(npw,lc%nc)!solution (imaginary part)  

      COMPLEX(kind=DP), ALLOCATABLE :: hmat(:,:),t(:)
      INTEGER, ALLOCATABLE :: ipiv(:)
      INTEGER :: ic,is
      REAL(kind=DP), ALLOCATABLE  :: t_re(:),t_im(:)
      INTEGER :: info


      allocate(hmat(lc%ns,lc%ns))
      allocate(ipiv(lc%ns))
      allocate(t(lc%ns),t_re(lc%ns),t_im(lc%ns))
      do ic=1,lc%nc

         hmat(1:lc%ns_chain(ic),1:lc%ns_chain(ic))=lc%h(1:lc%ns_chain(ic),1:lc%ns_chain(ic),ic)
         do is=1,lc%ns_chain(ic)
            hmat(is,is)=hmat(is,is)-freq
         enddo

!calculate t vectors projecting wfc_in
         CALL ZGEMM('C','N',lc%ns_chain(ic),1,npw,(1.d0,0.d0),lc%wfc(1,1,ic),npw,wfc_in(1,ic),npw,(0.d0,0.d0),t,lc%ns)
         t(1:lc%ns_chain(ic))=t(1:lc%ns_chain(ic))+conjg(t(1:lc%ns_chain(ic)))
         if(gstart==2) then
            t(1:lc%ns_chain(ic))=t(1:lc%ns_chain(ic))-conjg(lc%wfc(1,1:lc%ns_chain(ic),ic))*wfc_in(1,ic)
         endif
         call mp_sum(t,world_comm)


         
         CALL ZGESV(lc%ns_chain(ic),1,hmat,lc%ns,ipiv,t,lc%ns,info)
         if(info /= 0) then
           write(stdout,*) 'ZGESV info:', info
           FLUSH(stdout)
           stop
        endif
        if(.not.ionode) t=0.d0
        call mp_bcast(t,ionode_id, world_comm)
           
        t_re(1:lc%ns_chain(ic))=dble(t(1:lc%ns_chain(ic)))
        t_im(1:lc%ns_chain(ic))=aimag(t(1:lc%ns_chain(ic)))

        CALL DGEMM('N','N',2*npw,1,lc%ns_chain(ic),1.d0,lc%wfc(1,1,ic),2*npw,t_re,lc%ns,0.d0,wfc_out_re(1,ic),2*npw)
        CALL DGEMM('N','N',2*npw,1,lc%ns_chain(ic),1.d0,lc%wfc(1,1,ic),2*npw,t_im,lc%ns,0.d0,wfc_out_im(1,ic),2*npw)

      enddo


      deallocate(hmat,ipiv)
      deallocate(t,t_re,t_im)
    END SUBROUTINE solve_krylov

    SUBROUTINE product_chain_krylov(lc,v_states,ns_new,iv)
!multiply on real space grid the lanczos basis with the vstate vectors
!this is used for calculating W
      USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE gvect, ONLY : gstart
      USE io_global, ONLY : stdout
      USE fft_base,             ONLY : dfftp, dffts
      USE fft_interfaces,       ONLY : fwfft, invfft
      USE wavefunctions, ONLY : psic



      IMPLICIT NONE

      TYPE(lanczos_chain), INTENT(inout) :: lc !data to be modified
      REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,lc%nc)!the states to be multiplied with, namely valence states
      INTEGER :: ns_new!number of states to be update
      INTEGER, INTENT(in) :: iv!chain to be updated
      
      REAL(kind=DP), ALLOCATABLE :: psi_tmp(:,:)
      INTEGER :: it


      allocate(psi_tmp(dffts%nnr,2))

      lc%wfc(1:npw,lc%ns_chain(iv)-ns_new:lc%ns_chain(iv),iv)=lc%wfc0(1:npw,lc%ns_chain(iv)-ns_new:lc%ns_chain(iv),iv)

      do it=lc%ns_chain(iv)-ns_new,lc%ns_chain(iv),2
         call pc_operator(lc%wfc(1,it,iv),1,.false.)
         if(it/=lc%ns)  call pc_operator(lc%wfc(1,it+1,iv),1,.false.)
!put lanczos basis on real space
         psic(:)=(0.d0,0.d0)
         if(it/=lc%ns_chain(iv)) then
            psic(dffts%nl(1:npw))  = lc%wfc(1:npw,it,iv)+(0.d0,1.d0)*lc%wfc(1:npw,it+1,iv)
            psic(dffts%nlm(1:npw)) = CONJG( lc%wfc(1:npw,it,iv) )+(0.d0,1.d0)*conjg(lc%wfc(1:npw,it+1,iv))
         else
            psic(dffts%nl(1:npw))  = lc%wfc(1:npw,it,iv)
            psic(dffts%nlm(1:npw)) = CONJG( lc%wfc(1:npw,it,iv) )
         endif
         CALL invfft ('Wave', psic, dffts)
         if(it/=lc%ns_chain(iv)) then
            psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
            psi_tmp(1:dffts%nnr,2)= aimag(psic(1:dffts%nnr))
         else
            psi_tmp(1:dffts%nnr,1)= DBLE(psic(1:dffts%nnr))
         endif
!!product with valence states                     
         if(it/=lc%ns_chain(iv)) then
            psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
            psi_tmp(1:dffts%nnr,2)=psi_tmp(1:dffts%nnr,2)*v_states(1:dffts%nnr,iv)
         else
            psi_tmp(1:dffts%nnr,1)=psi_tmp(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
         endif
!bring back to G space

         if(it/=lc%ns_chain(iv)) then
            psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),psi_tmp(1:dffts%nnr,2))
         else
            psic(1:dffts%nnr)=cmplx(psi_tmp(1:dffts%nnr,1),0.d0)
         endif
         CALL fwfft ('Wave', psic, dffts)
         if(it/=lc%ns_chain(iv)) then
            lc%wfc(1:npw,it,iv)=0.5d0*(psic(dffts%nl(1:npw))+conjg(psic(dffts%nlm(1:npw))))
            lc%wfc(1:npw,it+1,iv)=(0.d0,-0.5d0)*(psic(dffts%nl(1:npw))-conjg(psic(dffts%nlm(1:npw))))
            if(gstart==2) lc%wfc(1,it,iv)=dble(lc%wfc(1,it,iv))
            if(gstart==2) lc%wfc(1,it+1,iv)=dble(lc%wfc(1,it+1,iv))
         else
            lc%wfc(1:npw,it,iv)=psic(dffts%nl(1:npw))
            if(gstart==2) lc%wfc(1,it,iv)=dble(lc%wfc(1,it,iv))
         endif

         
      enddo

      deallocate(psi_tmp)
      RETURN
      
    END SUBROUTINE product_chain_krylov


  END MODULE lanczos

