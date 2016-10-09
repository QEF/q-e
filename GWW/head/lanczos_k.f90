!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

subroutine lanczos_state_k(ik,nstates, nsteps,in_states,d,f,omat,dpsi_ipol, t_out)
!this subroutine perform nsteps collective lanczos iterations
!on orthonormal zstates state
! k points version

  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix
  USE kinds,    ONLY : DP
  USE wannier_gw
  USE gvect
  USE constants, ONLY : e2, pi, tpi, fpi
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wvfct,    ONLY : g2kin, npwx, nbnd
  USE wavefunctions_module, ONLY : evc, psic
  USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world, ONLY : mpime, world_comm
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE klist,                ONLY : xk,igk_k, ngk
  USE noncollin_module,     ONLY : noncolin, npol

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  INTEGER, INTENT(in) :: ik!k point
  INTEGER, INTENT(in) :: nstates!number of states
  INTEGER, INTENT(in) :: nsteps!number of Lanczos iteration to be performed
  COMPLEX(kind=DP), INTENT(in) :: in_states(npwx,nstates)!states for starting lanczos chains
  COMPLEX(kind=DP), INTENT(out) :: d(nsteps,nstates)!diagonal part
  COMPLEX(kind=DP), INTENT(out) :: f(nsteps,nstates)!off-diagonal part
  COMPLEX(kind=DP), INTENT(out) :: omat(nsteps,3,nstates)!overlaps
  COMPLEX(kind=DP), INTENT(in) :: dpsi_ipol(npwx,nstates,3)!other r|\Psi_v> states
  COMPLEX(kind=DP), INTENT(out) :: t_out(npwx,nsteps,nstates)!complete orthonormal basis

  COMPLEX(kind=DP), ALLOCATABLE :: psi_1(:,:),psi_2(:,:),psi_3(:,:),spsi(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: u_0(:,:),u_1(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: alpha(:),  delta(:)
  REAL(kind=DP), ALLOCATABLE :: gamma(:), n_1(:), beta(:)
  REAL(kind=DP), ALLOCATABLE :: c(:)
  

  INTEGER :: npw, is,ig,ii,jj,it,ipol
  INTEGER :: iunlan
  COMPLEX(kind=DP) :: csca

  allocate(psi_1(npwx,nstates),psi_2(npwx,nstates),psi_3(npwx,nstates))
  allocate(u_0(npwx,nstates),u_1(npwx,nstates))
  allocate(alpha(nstates),beta(nstates),gamma(nstates),n_1(nstates),delta(nstates))
  allocate(c(nstates))
  allocate(spsi(npwx,nstates))
 
  npw = ngk(ik)
  t_out(:,:,:)=(0.d0,0.d0)

  !first step
  psi_1(1:npw,1:nstates)=in_states(1:npw,1:nstates)
!calculate n_1
  n_1(:)=0.d0
  do is=1,nstates
     do ig=1,npw
        n_1(is)=n_1(is)+dble(conjg(psi_1(ig,is))*psi_1(ig,is))
     enddo
  enddo
  call mp_sum(n_1(:),world_comm)
  n_1(:)=dsqrt(n_1(:))
  do is=1,nstates
     psi_1(1:npw,is)=psi_1(1:npw,is)/n_1(is)
  enddo
!for h_psi allocations are required

  ALLOCATE( h_diag( npwx, npol ) )
  ALLOCATE( s_diag( npwx, npol ) )
  !
    !

!npw and igk should already been read!!

 
  IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
  do ig = 1, npw
     g2kin (ig) = ( (xk (1,ik ) + g (1,igk_k (ig,ik)) ) **2 + &
          (xk (2,ik ) + g (2,igk_k (ig,ik)) ) **2 + &
          (xk (3,ik ) + g (3,igk_k (ig,ik)) ) **2 ) * tpiba2
  enddo


!calculate H|\phi_i>
   call h_psi( npwx, npw, nstates, psi_1, u_0 )
   call s_psi( npwx, npw, nstates, psi_1, spsi ) ! is this needed?
   if(l_scissor) then
       call h_psi_scissor( ik,npwx, npw, nstates, psi_1, u_0 )
   endif
!calculate n_1
  n_1(:)=0.d0
  do is=1,nstates
     do ig=1,npw
        n_1(is)=n_1(is)+dble(conjg(u_0(ig,is))*u_0(ig,is))
     enddo
  enddo
  call mp_sum(n_1(:),world_comm)
  n_1(:)=dsqrt(n_1(:))
   write(stdout,*) 'Lanczos N1', n_1(:)
  FLUSH(stdout)


!calculate alpha
  alpha(:)=(0.d0,0.d0)
  do is=1,nstates
     do ig=1,npw
        alpha(is)=alpha(is)+conjg(psi_1(ig,is))*u_0(ig,is)
     enddo
  enddo
  call mp_sum(alpha(:),world_comm)
  alpha(:)=alpha(:)/n_1(:)
  write(stdout,*) 'Lanczos alpha', alpha(:)
  FLUSH(stdout)

!calculate psi_2 and beta
  do is=1,nstates
     psi_2(1:npw,is)=u_0(1:npw,is)/n_1(is)-alpha(is)*psi_1(1:npw,is)
  enddo
  
  beta(:)=0.d0
  do is=1,nstates
     do ig=1,npw
        beta(is)=beta(is)+dble(conjg(psi_2(ig,is))*psi_2(ig,is))
     enddo
  enddo
  call mp_sum(beta(:),world_comm)
  beta(:)=dsqrt(beta(:))
  write(stdout,*) 'Lanczos beta', beta(:)
  FLUSH(stdout)

  do is=1,nstates
     psi_2(:,is)=psi_2(:,is)/beta(is)
  enddo



!calculate d
  d(:,:)=0.d0
  do is=1,nstates
     do ig=1,npw
        d(1,is)=d(1,is)+conjg(psi_1(ig,is))*u_0(ig,is)
     enddo
  enddo
  do is=1,nstates
     call mp_sum(d(1,is),world_comm)
  enddo
  write(stdout,*) 'Lanczos Diagonal 1', d(1,:)
    FLUSH(stdout)

!calculate f

  f(:,:)=(0.d0,0.d0)
  do is=1,nstates
     do ig=1,npw
        f(1,is)=f(1,is)+conjg(psi_2(ig,is))*u_0(ig,is)
     enddo
     call mp_sum(f(1,is),world_comm)
  enddo

  write(stdout,*) 'ATTENZIONE1'
  FLUSH(stdout)
  omat(:,:,:)=(0.d0,0.d0)
  
  do is=1,nstates
     do ipol=1,3
        call zgemm('C','N',1,1,npw,(1.d0,0.d0),dpsi_ipol(:,is,ipol),npwx,psi_1(:,is),npwx,(0.d0,0.d0),omat(1,ipol,is),1)
        call zgemm('C','N',1,1,npw,(1.d0,0.d0),dpsi_ipol(:,is,ipol),npwx,psi_2(:,is),npwx,(0.d0,0.d0),omat(2,ipol,is),1)
     enddo
     t_out(1:npw,1,is)=psi_1(1:npw,is)
     t_out(1:npw,2,is)=psi_2(1:npw,is)
  end do
  call mp_sum(omat(1:2,1:3,1:nstates),world_comm)

   
  !do iterate
  do it=2,nsteps
     write(stdout,*) 'lanczos h_psi'
     FLUSH(stdout)

!calculate H|\phi_i+1>
     call h_psi( npwx, npw, nstates, psi_2, u_1 )
     call s_psi (npwx, npw, nstates, psi_2, spsi) ! is this needed?
     if(l_scissor) then
        call h_psi_scissor( ik,npwx, npw, nstates, psi_2, u_1 )
     endif

     write(stdout,*) 'lanczos alfa beta gamma'
     FLUSH(stdout)
!calculate n_1
     n_1(:)=0.d0
     do is=1,nstates
        do ig=1,npw
           n_1(is)=n_1(is)+dble(conjg(u_1(ig,is))*u_1(ig,is))
        enddo
     enddo
     call mp_sum(n_1(:),world_comm)
     n_1(:)=dsqrt(n_1(:))

!calculate alpha
     alpha(:)=(0.d0,0.d0)
     do is=1,nstates
        do ig=1,npw
           alpha(is)=alpha(is)+conjg(psi_1(ig,is))*u_1(ig,is)
        enddo
     enddo
     call mp_sum(alpha(:),world_comm)
     alpha(:)=alpha(:)/n_1(:)

!calculate beta
     delta(:)=(0.d0,0.d0)
     do is=1,nstates
        do ig=1,npw
           delta(is)=delta(is)+conjg(psi_2(ig,is))*u_1(ig,is)
        enddo
     enddo
     call mp_sum(delta(:),world_comm)
     delta(:)=delta(:)/n_1(:)

!calculate psi_3 and gamma
     do is=1,nstates
        psi_3(1:npw,is)=u_1(1:npw,is)/n_1(is)-alpha(is)*psi_1(1:npw,is)-delta(is)*psi_2(1:npw,is)
     enddo





     gamma(:)=0.d0
     do is=1,nstates
        do ig=1,npw
           gamma(is)=gamma(is)+dble(conjg(psi_3(ig,is))*psi_3(ig,is))
        enddo
     enddo
     call mp_sum(gamma(:),world_comm)
     gamma(:)=dsqrt(gamma(:))
     do is=1,nstates
        psi_3(:,is)=psi_3(:,is)/gamma(is)
     enddo
     write(stdout,*) 'lanczos d f omat'
     FLUSH(stdout)


!calculate d
     do is=1,nstates
        do ig=1,npw
           d(it,is)=d(it,is)+dble(conjg(psi_2(ig,is))*u_1(ig,is))
        enddo
         call mp_sum(d(it,is),world_comm)
     enddo
    

!calculate f
     do is=1,nstates
        do ig=1,npw
           f(it,is)=f(it,is)+conjg(psi_3(ig,is))*u_1(ig,is)
        enddo
        call mp_sum(f(it,is),world_comm)
     enddo
  

     if(it/=nsteps) then
        do is=1,nstates
           do ipol=1,3
              call zgemm('C','N',1,1,npw,(1.d0,0.d0),dpsi_ipol(:,is,ipol),npwx,psi_3(:,is),npwx,(0.d0,0.d0),&
                   &omat(it+1,ipol,is),1)
           enddo
           t_out(1:npw,it+1,is)=psi_3(1:npw,is)
        end do
        call mp_sum(omat(it+1,1:3,1:nstates),world_comm)
     endif
        


!update arrays
     psi_1(:,:)=psi_2(:,:)
     psi_2(:,:)=psi_3(:,:)
     u_0(:,:)=u_1(:,:)


  enddo
  
  deallocate(psi_1,psi_2,psi_3)
  deallocate(u_0,u_1)
  deallocate(alpha,beta,gamma,n_1)
  deallocate(c)
  deallocate(delta)
  deallocate(h_diag,s_diag)
  deallocate(spsi)

  return
end subroutine lanczos_state_k


subroutine  h_psi_scissor( ik,lda, n, m, psi, hpsi )
!NOT_TO_BE_INCLUDED_START
!add to hpsi part dur to self-consistent GW calculation
! ... input
! ...    lda   leading dimension of arrays psi, spsi, hpsi
! ...    n     true dimension of psi, spsi, hpsi
! ...    m     number of states psi
! ...    psi
! ... output:
! ...    hpsi  H*psi   

  USE kinds,    ONLY : DP
  USE gvect,    ONLY : gstart
  USE wvfct,    ONLY : npwx, npw, nbnd,et
  USE wavefunctions_module, ONLY : evc
  USE wannier_gw, ONLY : scissor
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE control_lr,           ONLY : nbnd_occ
  USE constants, ONLY : rytoev

  implicit none

  INTEGER, INTENT(in)  :: ik!k point
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(kind=DP), INTENT(IN)  :: psi(lda,m) 
  COMPLEX(kind=DP), INTENT(OUT) :: hpsi(lda,m)   

  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: prod(:,:)



  allocate(prod(nbnd_occ(ik),m))
  prod=0.d0
  call dgemm('T','N', nbnd_occ(ik),m,2*npw,2.d0,evc,2*npwx,psi,2*lda,0.d0,prod,nbnd_occ(ik))
  do ii=1,nbnd_occ(ik)
     do jj=1,m
        if(gstart==2) prod(ii,jj)=prod(ii,jj)-dble(conjg(evc(1,ii))*psi(1,jj))
     enddo
  enddo
  call mp_sum(prod,world_comm)

  do jj=1,m
     do ii=1,nbnd_occ(ik)
        prod(ii,jj)=prod(ii,jj)*(scissor(1)-scissor(2))/rytoev
     enddo
  enddo
  call dgemm('N','N',2*npw,m,nbnd_occ(ik),1.d0,evc,2*npwx,prod,nbnd_occ(ik),1.d0+scissor(2)/rytoev,hpsi,2*lda)


  deallocate(prod)
  return
!NOT_TO_BE_INCLUDED_END
end subroutine h_psi_scissor
