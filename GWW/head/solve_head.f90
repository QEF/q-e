!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!-----------------------------------------------------------------------
subroutine solve_head
  !-----------------------------------------------------------------------
  !
  !calculates the head and wings of the dielectric matrix
  !
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode,ionode_id
  USE io_files,              ONLY : diropn,prefix, tmp_dir
  use pwcom
  USE check_stop,            ONLY : max_seconds
  USE wavefunctions_module,  ONLY : evc
  USE kinds,                 ONLY : DP
  USE becmod,                ONLY : becp,calbec
  USE uspp_param,            ONLY : nhm
  use phcom
  USE wannier_gw,           ONLY : n_gauss, omega_gauss, grid_type,&
                                   nsteps_lanczos,second_grid_n,second_grid_i,&
                                   l_scissor,scissor,len_head_block_freq, &
                                   len_head_block_wfc
  USE control_ph,           ONLY : tr2_ph
  USE gvect,                ONLY : ig_l2g
  USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world,             ONLY : world_comm, mpime, nproc
  USE uspp,                 ONLY : nkb, vkb
!  USE symme, ONLY: s
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm
  USE symme, only : crys_to_cart, symmatrix
  USE mp_wave, ONLY : mergewf,splitwf
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE buffers,              ONLY : get_buffer
  USE constants,            ONLY : rytoev

  use qpoint,                ONLY : npwq, nksq
  use control_lr,            ONLY : nbnd_occ, lgamma

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
 
 

 
  complex(DP) , allocatable ::    ps (:,:)

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol,jpol, ibnd, jbnd, iter, lter, &
       ik, ig, irr, ir, is, nrec, ios
  ! counters
  integer :: ltaver, lintercall

  real(DP) :: tcpu, get_clock
  ! timing variables
  
  REAL(kind=DP), ALLOCATABLE :: head(:,:),head_tmp(:)
  COMPLEX(kind=DP) :: sca, sca2
  REAL(kind=DP), ALLOCATABLE :: x(:),w(:), freqs(:)
 ! COMPLEX(kind=DP), ALLOCATABLE :: e_head(:,:)!wing of symmetric dielectric matrix (for G of local processor)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head_g(:),e_head_g_tmp(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head_pol(:,:,:)
  INTEGER :: i, j,k,iun
  REAL(kind=DP) :: ww, weight
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:)
  COMPLEX(kind=DP), ALLOCATABLE :: psi_v(:,:), prod(:)
  COMPLEX(kind=DP), ALLOCATABLE :: pola_charge(:,:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: dpsi_ipol(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: epsilon_g(:,:,:)
  INTEGER :: i_start,idumm,idumm1,idumm2,idumm3,ii
  REAL(kind=DP) :: rdumm
  COMPLEX(kind=DP), ALLOCATABLE :: d(:,:),f(:,:),omat(:,:,:)
  INTEGER :: iv, info
  COMPLEX(kind=DP), ALLOCATABLE :: z_dl(:),z_d(:),z_du(:),z_b(:)
  COMPLEX(kind=DP) :: csca, csca1
  COMPLEX(kind=DP), ALLOCATABLE :: t_out(:,:,:), psi_tmp(:)
  INTEGER :: n
  INTEGER :: npwx_g

  INTEGER :: ib,lenb,first_b,last_b,n_block,nbnd_block

  INTEGER :: freq_block, m_block,first_f,last_f,lenf,im


  write(stdout,*) 'Routine solve_head'
  FLUSH(stdout)

  if(grid_type==5) then
     n=n_gauss
     n_gauss=n+second_grid_n*(1+second_grid_i*2)
  endif

  !allocate(e_head(npw,n_gauss+1))
  !allocate(e_head_pol(ngm,n_gauss+1,3))
  !e_head(:,:) =(0.d0,0.d0)  
  allocate(x(2*n_gauss+1),w(2*n_gauss+1), freqs(n_gauss+1))
  allocate(head(n_gauss+1,3),head_tmp(n_gauss+1))
  head(:,:)=0.d0
  allocate(prod(dfftp%nnr))
  allocate (tmp_g(ngm))
 ! allocate( pola_charge(dfftp%nnr,nspin,3,n_gauss+1))
  allocate(epsilon_g(3,3,n_gauss+1))
  allocate(psi_tmp(npwx))




  epsilon_g(:,:,:)=0.d0
!  e_head_pol(:,:,:)=0.d0
!  pola_charge(:,:,:,:)=0.d0


!setup Gauss Legendre frequency grid
!IT'S OF CAPITAL IMPORTANCE TO NULLIFY THE FOLLOWING ARRAYS
  x(:)=0.d0
  w(:)=0.d0
  if(grid_type==0) then
     call legzo(n_gauss*2+1,x,w)
     freqs(1:n_gauss+1)=-x(n_gauss+1:2*n_gauss+1)*omega_gauss
  else if(grid_type==2) then
     call legzo(n_gauss,x,w)
     freqs(1) = 0.d0
     freqs(2:n_gauss+1)=(1.d0-x(1:n_gauss))*omega_gauss/2.d0
  else if(grid_type==3) then!equally spaced grid
     freqs(1) = 0.d0
     do i=1,n_gauss
        freqs(1+i)=omega_gauss*dble(i)/dble(n_gauss)
     enddo
  else  if(grid_type==4) then!equally spaced grid shifted of 1/2
     freqs(1) = 0.d0
     do i=1,n_gauss
        freqs(i+1)=(omega_gauss/dble(n_gauss))*dble(i)-(0.5d0*omega_gauss/dble(n_gauss))
     enddo
  else!equally spaced grid more dense at -1 , 0 and 1
     freqs(1)=0.d0
          
     ii=2
     do i=1,second_grid_n
        freqs(ii)=(omega_gauss/dble(2*second_grid_n*n))*dble(i)-0.5d0*omega_gauss/dble(2*second_grid_n*n)
        ii=ii+1
     enddo
     do j=1,second_grid_i
        do i=1,second_grid_n
           freqs(ii)=(omega_gauss/dble(2*second_grid_n*n))*dble(i+second_grid_n+2*second_grid_n*(j-1))&
      &-0.5d0*omega_gauss/dble(2*second_grid_n*n)
           ii=ii+1
        enddo
        freqs(ii)=omega_gauss/dble(n)*dble(j)
        ii=ii+1
        do i=1,second_grid_n
           freqs(ii)=(omega_gauss/dble(2*second_grid_n*n))*dble(i+2*second_grid_n*j)&
    &-0.5d0*omega_gauss/dble(2*second_grid_n*n)
           ii=ii+1
        enddo
     enddo
     do i=second_grid_i+1,n
        freqs(ii)=omega_gauss/dble(n)*dble(i)
        ii=ii+1
     enddo
     


  endif
  do i=1,n_gauss+1
     write(stdout,*) 'Freq',i,freqs(i)
  enddo
  FLUSH( stdout )
  
  deallocate(x,w)
  head(:,:)=0.d0





  !if (lsda) call errore ('solve_head', ' LSDA not implemented', 1)

  call start_clock ('solve_head')
 
 
 
  allocate (ps  (nbnd,nbnd))    
  ps (:,:) = (0.d0, 0.d0)
 
  
 
  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL DIROPN (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  !

  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if (degauss.ne.0.d0.or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)
  !

  !
  !   only one iteration is required
  !

  if(.not.l_scissor) scissor=0.d0

!loop on frequency blocks
  if(len_head_block_freq==0) then
     freq_block=n_gauss+1
     m_block=1
  else
     m_block=(n_gauss+1)/len_head_block_freq
     if(m_block*len_head_block_freq < (n_gauss+1)) m_block=m_block+1
     if(len_head_block_freq>(n_gauss+1)) then
        freq_block=n_gauss+1
     else
        freq_block=len_head_block_freq
     endif
  endif

!loop on frequency blocks
  do im=1,m_block

     epsilon_g=0.d0

   write(stdout,*) 'FREQUENCY BLOCK : ', im
   
     first_f=(im-1)*freq_block+1
     last_f=im*freq_block
     if(last_f>(n_gauss+1)) last_f=n_gauss+1
     lenf=last_f-first_f+1
     
     allocate( pola_charge(dfftp%nnr,nspin,3,lenf))  
     allocate(e_head_pol(ngm,lenf,3))

     e_head_pol(:,:,:)=0.d0
     pola_charge(:,:,:,:)=0.d0



!loop on k points
     do ik=1, nksq
        if(len_head_block_wfc==0) then
           nbnd_block=nbnd_occ(ik)
           n_block=1
        else
           n_block=nbnd_occ(ik)/len_head_block_wfc
           if(n_block*len_head_block_wfc < nbnd_occ(ik)) n_block=n_block+1
           if(len_head_block_wfc>nbnd_occ(ik)) then
              nbnd_block=nbnd_occ(ik)
           else
              nbnd_block=len_head_block_wfc
           endif
        endif



        write(stdout,*) 'ik:', ik
        FLUSH(stdout)
        weight = wk (ik)
        ww = fpi * weight / omega
      
        if (lsda) current_spin = isk (ik)
        npw = ngk(ik)
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
        if (nksq.gt.1)  call get_buffer(evc, lrwfc, iuwfc, ik)
        npwq = npw
        call init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)



     !
     ! compute the kinetic energy
     !
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ik ) + g (1,igk_k (ig,ik)) ) **2 + &
                (xk (2,ik ) + g (2,igk_k (ig,ik)) ) **2 + &
                (xk (3,ik ) + g (3,igk_k (ig,ik)) ) **2 ) * tpiba2
        enddo
      !
        do ib=1,n_block

           write(stdout,*) 'BLOCK : ', ib
        
           first_b=(ib-1)*nbnd_block+1
           last_b=ib*nbnd_block
           if(last_b>nbnd_occ(ik)) last_b=nbnd_occ(ik)
           lenb=last_b-first_b+1
           
           allocate (dpsi_ipol(npwx,lenb,3))
           allocate(t_out(npwx,nsteps_lanczos,lenb))
           allocate(psi_v(dffts%nnr, lenb))

           dpsi_ipol(:,:,:)=(0.d0,0.d0)
        

           !trasform valence wavefunctions to real space                                                                                        
           do ibnd=1,lenb
              psi_v(:,ibnd) = ( 0.D0, 0.D0 )
              psi_v(nls(igk_k(1:npw,ik)),ibnd) = evc(1:npw,first_b+ibnd-1)
              CALL invfft ('Wave',  psi_v(:,ibnd), dffts)
           enddo



!loop on carthesian directions
           do ipol = 1,3
              write(stdout,*) 'ipol:', ipol
              FLUSH(stdout)
        !
        ! computes/reads P_c^+ x psi_kpoint into dvpsi array
        !

              do jpol=1,3
         
                 call dvpsi_e (ik, jpol)
              
          !
        ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
        !
                 CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
                      (1.d0,0.d0), evc(1,1), npwx, dvpsi(1,1), npwx, (0.d0,0.d0), &
                      ps(1,1), nbnd )
#if defined(__MPI)
           !call reduce (2 * nbnd * nbnd_occ (ik), ps)
                 call mp_sum(ps(1:nbnd_occ (ik),1:nbnd_occ (ik)),world_comm)
#endif
        ! dpsi is used as work space to store S|evc>
        !

                 CALL calbec(npw,vkb,evc,becp,nbnd_occ(ik))
                 CALL s_psi (npwx, npw, nbnd_occ(ik), evc, dpsi)
           !
        
        ! |dvpsi> = - (|dvpsi> - S|evc><evc|dvpsi>)
        ! note the change of sign!
           !
                 CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
                    (1.d0,0.d0), dpsi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
                    dvpsi(1,1), npwx )
!create lanczos chain for dvpsi
                 dpsi_ipol(1:npw,1:lenb,jpol)=dvpsi(1:npw,first_b:last_b)
              enddo
              dvpsi(1:npw,1:lenb)=dpsi_ipol(1:npw,1:lenb,ipol)
            
            
              allocate(d(nsteps_lanczos,lenb),f(nsteps_lanczos,lenb))
              allocate(omat(nsteps_lanczos,3,lenb))
              write(stdout,*) 'before lanczos_state_k'
     



              call lanczos_state_k(ik,lenb, nsteps_lanczos ,dvpsi,d,f,omat,dpsi_ipol,t_out)
              write(stdout,*) 'after lanczos_state_k'
!loop on frequency
              allocate(z_dl(nsteps_lanczos-1),z_d(nsteps_lanczos),z_du(nsteps_lanczos-1),z_b(nsteps_lanczos))
              do i=1,lenf
!loop on valence states
                 do iv=1,lenb
!invert Hamiltonian
                    z_dl(1:nsteps_lanczos-1)=conjg(f(1:nsteps_lanczos-1,iv))
                    z_du(1:nsteps_lanczos-1)=f(1:nsteps_lanczos-1,iv)
                    z_d(1:nsteps_lanczos)=d(1:nsteps_lanczos,iv)+dcmplx(-et(first_b+iv-1,ik)-scissor(1)/rytoev,freqs(first_f+i-1))
                    z_b(:)=(0.d0,0.d0)
                    z_b(1)=dble(omat(1,ipol,iv))
                    call zgtsv(nsteps_lanczos,1,z_dl,z_d,z_du,z_b,nsteps_lanczos,info)
                    if(info/=0) then
                       write(stdout,*) 'problems with ZGTSV'
                       FLUSH(stdout)
                       stop
                    endif
                    do jpol=1,3
!multiply with overlap factors
                       call zgemm('T','N',1,1,nsteps_lanczos,(1.d0,0.d0),omat(:,jpol,iv),nsteps_lanczos&
     &,z_b,nsteps_lanczos,(0.d0,0.d0),csca,1)
!update epsilon array NO SYMMETRIES for the moment
                       epsilon_g(jpol,ipol,i)=epsilon_g(jpol,ipol,i)+4.d0*ww*dble(csca)
                    enddo
!update part for wing calculation 
                    call zgemm('N','N',npw,1,nsteps_lanczos,(1.d0,0.d0),t_out(:,:,iv),npwx,z_b,nsteps_lanczos,&
                         &(0.d0,0.d0),psi_tmp,npwx) 
!fourier trasform
                    prod(:) = ( 0.D0, 0.D0 )
                    prod(nls(igk_k(1:npw,ik))) = psi_tmp(1:npw)
                    CALL invfft ('Wave', prod, dffts)
           

!      product dpsi * psi_v
                    prod(1:dffts%nnr)=conjg(prod(1:dffts%nnr))*psi_v(1:dffts%nnr,iv)
                    if(doublegrid) then
                       call cinterpolate(prod,prod,1)
                    endif

!US part STLL TO BE ADDED!!
                    pola_charge(1:dffts%nnr,1,ipol,i)=pola_charge(1:dffts%nnr,1,ipol,i)-prod(1:dffts%nnr)*ww


                 enddo
              enddo
              deallocate(z_dl,z_d,z_du,z_b)
              deallocate(d,f,omat)
           enddo
           deallocate(dpsi_ipol)
           deallocate(t_out)
           deallocate(psi_v)
        end do!on ib
     enddo! on ik

!print out results



!
!      symmetrize
!

     do i=1,lenf
        WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
        WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon_g(ipol,jpol,i),&
          &                                ipol=1,3),jpol=1,3)

  !   call symtns (epsilon_g(:,:,i), nsym, s)
  !
  !    pass to cartesian axis
  !
        WRITE( stdout,'(/,10x,"Symmetrized in crystal axis ",/)')
        WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon_g(ipol,jpol,i),&
       &                                ipol=1,3),jpol=1,3)
  !   call trntns (epsilon_g(:,:,i), at, bg, 1)

        call crys_to_cart ( epsilon_g(:,:,i) )
        call symmatrix ( epsilon_g(:,:,i))
  !
  ! add the diagonal part
  !

  !
  !  and print the result
  !
        WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')
        
        WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon_g(ipol,jpol,i), ipol=1,3), jpol=1,3)

        head(first_f+i-1,1)=epsilon_g(1,1,i)
        head(first_f+i-1,2)=epsilon_g(2,2,i)
        head(first_f+i-1,3)=epsilon_g(3,3,i)


#if defined(__MPI)
        call mp_sum ( pola_charge(:,:,:,i) , inter_pool_comm )
        call psyme (pola_charge(:,:,:,i))
#else
        call syme (pola_charge(:,:,:,i))
#endif
     
        do ipol=1,3
           CALL fwfft ('Dense',  pola_charge(1:dfftp%nnr,1,ipol,i), dfftp)
           tmp_g(:)=(0.d0,0.d0)
           tmp_g(gstart:ngm)=pola_charge(nl(gstart:ngm),1,ipol,i) 
          
!loop on frequency
           do ig=gstart,ngm
              e_head_pol(ig,i,ipol)=-4.d0*tmp_g(ig)
           enddo
        enddo
       

     enddo


!writes on file wings

!collect data
 
  


!calculate total number of G for wave function
     npwx_g=ngm
     call mp_sum(npwx_g,world_comm)
     allocate(e_head_g(ngm_g))

     if(ionode) then
        iun =  find_free_unit()
        if(im==1) then
           open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.e_head', status='unknown',form='unformatted')
           write(iun) n_gauss
           write(iun) omega_gauss
           write(iun) freqs(1:n_gauss+1)
           write(iun) npwx_g
        else
           open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.e_head', status='old',form='unformatted',&
                &position='append')
        endif
     endif

     call mp_barrier( world_comm )
     
     
     do i=1,lenf
        do ipol=1,3
           e_head_g(:)=(0.d0,0.d0)
           call mergewf(e_head_pol(:,i,ipol),e_head_g ,ngm,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
           if(ionode) then
              write(iun) e_head_g(1:npwx_g)
           endif
        enddo
     enddo
     call mp_barrier( world_comm )
     write(stdout,*) 'ATT02'
     if(ionode) close(iun)

     call mp_barrier( world_comm )
     write(stdout,*) 'ATT1'
     deallocate(pola_charge)
     deallocate(e_head_pol)
     deallocate(e_head_g)
  end do!im

!writes on file head

     if(ionode) then
        iun =  find_free_unit()
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.head', status='unknown',form='unformatted')
        write(iun) n_gauss
        write(iun) omega_gauss
        write(iun) freqs(1:n_gauss+1)
        write(iun) head(1:n_gauss+1,1)
        write(iun) head(1:n_gauss+1,2)
        write(iun) head(1:n_gauss+1,3)
        close(iun)
     endif



  deallocate(psi_tmp)
  deallocate(prod)
  
  deallocate (ps)

 
  
    

  deallocate(head,head_tmp,freqs)
  deallocate( tmp_g)
  deallocate(epsilon_g)

  

   call mp_barrier( world_comm )
   write(stdout,*) 'ATT2'

  call stop_clock ('solve_head')
  return
end subroutine solve_head

