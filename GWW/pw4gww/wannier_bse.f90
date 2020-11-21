!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Modified by Joshua Elliott November 2020 as JDE
!
! This subroutine computes the overlap between Wannier orbitals Ovv' and
! computes the (v*w_v*w_v')(r) term for each vv' such that Ovv'>s_bse, and
! writes to disk Ovv' and (v*w_v*w_v')(r) ,and z_beta_v_v'=v*phi_beta*wv*wv'

subroutine wannier_bse(ispin,w_wfcs,o_mat)

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY :  prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw, ONLY : num_nbndv,dual_bse,s_bse,l_truncated_coulomb,truncation_radius,vg_q,&
                          max_ngm,numw_prod, pmat_type, l_no_GW_just_screening, l_no_GW_bare_coulomb,& ! JDE
                          no_GW_cg_maxit, no_GW_cg_threshold, ewvc                                     ! JDE
   USE fft_custom_gwl
   USE wvfct,    ONLY : npwx, npw, nbnd, nbndx  ! JDE
   USE gvecw,    ONLY : ecutwfc
   USE mp_pools, ONLY : intra_pool_comm
   USE mp_world, ONLY : mpime, nproc, world_comm
   USE mp,             ONLY : mp_sum, mp_bcast ! JDE
   USE gvect
   USE wavefunctions, ONLY :  psic
   USE constants, ONLY : e2, fpi
   USE cell_base, ONLY: tpiba,tpiba2,omega
   USE mp_wave, ONLY : mergewf,splitwf






   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   INTEGER, intent(in) :: ispin
   COMPLEX(kind=DP), intent(in) :: w_wfcs(npw,num_nbndv(ispin))
   REAL(kind=DP), intent(out) :: o_mat(num_nbndv(ispin),num_nbndv(ispin))
   TYPE(fft_cus) :: fc

   COMPLEX(kind=DP), allocatable :: w_wfcs_t(:,:)
   REAL(kind=DP), ALLOCATABLE :: w_wfcs_r(:,:)
   REAL(kind=DP), ALLOCATABLE :: w_wfcs_2(:,:)
   REAL(kind=DP), ALLOCATABLE :: ww_prod(:)
   COMPLEX(kind=DP), allocatable :: ww_prodg(:),ww_prodg2(:)
   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
   COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
   REAL(kind=DP), ALLOCATABLE :: fac(:)
   REAL(kind=DP), ALLOCATABLE :: z(:)
   INTEGER, ALLOCATABLE ::iww(:)


   INTEGER :: ii,ig, iv,jv,np,np_max,i
   INTEGER :: iunu, iungprod,iunz,iuni
   REAL(kind=DP) :: qq
   LOGICAL :: exst
   logical :: debug, debug_operator  ! JDE 
   integer :: iundebug

! JDE
   INTEGER :: wannier_product_unit, loop_wannier_prod, &
              wannier_write_product_unit, wwwdebug, iunwww, &
              start_loop_wannier
   EXTERNAL :: operator_1_vp, operator_debug
   COMPLEX(kind=DP), allocatable :: ww_read_prodg2(:), ww_solved_prodg2(:)
   REAL(KIND=dp), ALLOCATABLE :: overlap(:,:), diagonal(:), e(:)
   COMPLEX(KIND=dp), ALLOCATABLE :: evc_jsh(:,:)
   COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: x, b, u, h, Ah, evc, pu
   INTEGER :: cg_iter
   LOGICAL :: cgsolve_conv
! JDE


   fc%ecutt=ecutwfc
   fc%dual_t=dual_bse
    
   debug=.true.
   debug_operator=.false.

! FFT the wannier function to r-space (dual grid)

   write(stdout,*) 'Call initialize_fft_custom'
   call initialize_fft_custom(fc)
   allocate(w_wfcs_t(fc%npwt,num_nbndv(ispin)))
   allocate( evc_g(fc%ngmt_g ) )
   allocate(w_wfcs_r(fc%nrxxt,num_nbndv(ispin)))
   allocate(w_wfcs_2(fc%nrxxt,num_nbndv(ispin)))


   if(fc%dual_t==4.d0) then
           w_wfcs_t(1:fc%npwt,1:num_nbndv(ispin))= w_wfcs(1:fc%npwt,1:num_nbndv(ispin))
   else
           do ii=1, num_nbndv(ispin)
              call mergewf(w_wfcs(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
              call splitwf(w_wfcs_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           enddo
   endif

   do ii=1,num_nbndv(ispin),2
      psic(1:fc%nrxxt)=(0.d0,0.d0)
      if (ii==num_nbndv(ispin)) then
          psic(fc%nlt(1:fc%npwt))  = w_wfcs_t(1:fc%npwt,ii)
          psic(fc%nltm(1:fc%npwt)) = CONJG( w_wfcs_t(1:fc%npwt,ii) )
      else
          psic(fc%nlt(1:fc%npwt))=w_wfcs_t(1:fc%npwt,ii)+(0.d0,1.d0)*w_wfcs_t(1:fc%npwt,ii+1)
          psic(fc%nltm(1:fc%npwt)) =CONJG(w_wfcs_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(w_wfcs_t(1:fc%npwt,ii+1))
      endif
      
      CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
      
      w_wfcs_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
      if(ii/=num_nbndv(ispin)) w_wfcs_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))

      w_wfcs_2(1:fc%nrxxt,ii)=w_wfcs_r(1:fc%nrxxt,ii)**2
      if(ii/=num_nbndv(ispin)) w_wfcs_2(1:fc%nrxxt,ii+1)=w_wfcs_r(1:fc%nrxxt,ii+1)**2
   enddo

! compute the overlap matrix o_vv'
   call dgemm('T','N',num_nbndv(ispin),num_nbndv(ispin),fc%nrxxt,1.d0, w_wfcs_2, &
        & fc%nrxxt,w_wfcs_2, fc%nrxxt, 0.d0, o_mat,num_nbndv(ispin))        
   call mp_sum(o_mat,world_comm)
   
   o_mat(1:num_nbndv(ispin),1:num_nbndv(ispin))= &
        & o_mat(1:num_nbndv(ispin),1:num_nbndv(ispin))/(fc%nr1t*fc%nr2t*fc%nr3t)

   do iv=1,num_nbndv(ispin) 
     do jv=1,num_nbndv(ispin)
        write(stdout,*) 'iv,jv,o_mat',iv,jv,o_mat(iv,jv)
     enddo
   enddo

   FLUSH(stdout)
! write it on disk
   if(ionode) then 
     iunu = find_free_unit()
  
     if (ispin==1) open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wbse1',status='unknown',form='unformatted')
     if (ispin==2) open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wbse2',status='unknown',form='unformatted')

     write(iunu) num_nbndv(ispin)
     write(iunu) s_bse
     do ii=1,num_nbndv(ispin)
       write(iunu) o_mat(1:num_nbndv(ispin),ii)
     enddo
     close(iunu)
   
   endif
! if Ovv'> s_bse compute in G-space the v*w_v*w_v' product

   iungprod = find_free_unit()
   if (ispin==1)  CALL diropn( iungprod, 'vww_bse1.',npw*2, exst)
   if (ispin==2)  CALL diropn( iungprod, 'vww_bse2.',npw*2, exst)

   if(debug) then
      iundebug = find_free_unit()
      open(iundebug,file='vww_pw4gww.dat')
   endif
! compute V(G)

   allocate(fac(npw))
   if(l_truncated_coulomb) then
     do ig=1,npw
        qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
        if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
        else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
        endif
     enddo
     fac(:)=fac(:)/omega
   else
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif

   allocate(ww_prod(fc%nrxxt))
   allocate(ww_prodg(fc%npwt))
   allocate(ww_prodg2(npw))

   ii=0
   do iv=1, num_nbndv(ispin)
     do jv=1, num_nbndv(ispin)
       if (o_mat(iv,jv)>=s_bse) then
        ii=ii+1 
        ww_prod(1:fc%nrxxt)= w_wfcs_r(1:fc%nrxxt,iv)* w_wfcs_r(1:fc%nrxxt,jv)
        psic(1:fc%nrxxt)=ww_prod(1:fc%nrxxt)
        CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
        ww_prodg(1:fc%npwt) = psic(fc%nlt(1:fc%npwt))
      

        if(fc%dual_t==4.d0) then
           ww_prodg2(:)=ww_prodg(:)
        else
           call mergewf(ww_prodg,evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           call splitwf(ww_prodg2,evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
        endif 

        ww_prodg2(1:npw)=ww_prodg2(1:npw)*fac(1:npw)
        call davcio(ww_prodg2,npw*2,iungprod,ii,1)

        if(debug) then
          if(ionode) then
             write(iundebug,*) npw
             do i=1,npw
                write(iundebug,*) ww_prodg2(i)
             enddo
          endif
        endif

       endif
    enddo
 enddo
   write(stdout,*) 'bse ii found=',ii 

  write(stdout,*) 'max_ngm=',max_ngm  
  write(stdout,*) 'npw=',npw  
  FLUSH(stdout)

  close(iungprod)
! JDE start
! Obtain polarization basis
! Iterative method
  IF (l_no_GW_just_screening) THEN
     allocate(ww_read_prodg2(npw))
     allocate(ww_solved_prodg2(npw))

     if (debug) then
        wwwdebug = find_free_unit()
        open(wwwdebug,file='Www_pw4gww.dat')
     endif

     ! Open files for wannier products
     wannier_product_unit = find_free_unit()
     IF (ispin .EQ. 1) CALL diropn( wannier_product_unit, 'vww_bse1.',npw*2, exst)
     IF (ispin .EQ. 2) CALL diropn( wannier_product_unit, 'vww_bse2.',npw*2, exst)
     wannier_write_product_unit = find_free_unit()
     IF (ispin .EQ. 1) CALL diropn( wannier_write_product_unit, 'Www_bse1.',npw*2, exst)
     IF (ispin .EQ. 2) CALL diropn( wannier_write_product_unit, 'Www_bse2.',npw*2, exst)

     ! Allocate variables needed for cgsolve routine 
     IF (.NOT. l_no_GW_bare_coulomb) THEN
        ALLOCATE(overlap(nbndx,1))
        ALLOCATE(diagonal(npw))
        ALLOCATE(evc_jsh(npw,1))
        ALLOCATE(e(1))
        ALLOCATE(b(npw,1))
        ALLOCATE(u(npw,1))
        ALLOCATE(h(npw,1))
        ALLOCATE(Ah(npw,1))
        ALLOCATE(pu(npw,1))

        evc_jsh= 0.d0
        diagonal= 0.d0
        overlap= 0.d0
! BE VERY CAREFUL PMAT_TYPE IS A GLOBAL VARIABLE
WRITE(stdout,*) 'Note: we have changed pmat_type on the fly!'
pmat_type = 0
! BE VERY CAREFUL PMAT_TYPE IS A GLOBAL VARIABLE
     END IF

     INQUIRE(FILE=TRIM(tmp_dir)//TRIM(prefix)//'.restart_Www_stat', EXIST=exst)
     IF (exst) THEN
        IF (ionode) THEN
           iunwww = find_free_unit()
           OPEN(UNIT=iunwww,FILE=TRIM(tmp_dir)//TRIM(prefix)//'.restart_Www_stat')
           READ(iunwww,*) start_loop_wannier
           CLOSE(iunwww)
        END IF
        start_loop_wannier = start_loop_wannier + 1 
        CALL mp_bcast(start_loop_wannier, ionode_id, world_comm)
     ELSE
        start_loop_wannier = 1 
     END IF

     DO loop_wannier_prod = start_loop_wannier, ii
        CALL davcio(ww_read_prodg2,npw*2,wannier_product_unit,loop_wannier_prod,-1)
        ww_solved_prodg2(:) = 0.d0
        ! NOTE: read Vww :: do not need to apply Coulomb operator
        !ww_read_prodg2(1:npw) =  ww_read_prodg2(1:npw) * fac(1:npw) 

        IF (.NOT. l_no_GW_bare_coulomb) THEN
           IF (debug_operator) THEN
              WRITE(*,*) 'before cg'
              CALL cgsolve(operator_debug, npw, evc_jsh, npw, 1, overlap, &
                           1, .TRUE., .FALSE., diagonal, &
                           .FALSE.,e,ww_read_prodg2,u,h,Ah,pu,no_GW_cg_maxit,no_GW_cg_threshold,cg_iter,ww_solved_prodg2 )
           ELSE 

              WRITE(stdout,*) 'CGSOLVE ww PAIR:', loop_wannier_prod
              cgsolve_conv = .false.
              DO WHILE (.NOT. cgsolve_conv)
                 overlap=0.d0
                 diagonal=0.d0
                 evc_jsh=0.d0
                 CALL cgsolve(operator_1_vp, npw, evc_jsh, npw, 1, overlap, 1, .TRUE., .FALSE., diagonal, .FALSE.,& 
                              e, ww_read_prodg2, u, h, Ah, pu, no_GW_cg_maxit, no_GW_cg_threshold, cg_iter, ww_solved_prodg2 )
                 IF (cg_iter .LT. no_GW_cg_maxit) THEN
                    cgsolve_conv = .true.
                 END IF
              END DO
           END IF
        ELSE
           ww_solved_prodg2(1:npw) = ww_read_prodg2(1:npw)
        END IF

        CALL davcio(ww_solved_prodg2,npw*2,wannier_write_product_unit,loop_wannier_prod,1)

        if(debug) then
          if(ionode) then
             write(wwwdebug,*) npw
             do i=1,npw
                write(wwwdebug,*) ww_solved_prodg2(i)
             enddo
          endif
        endif

        IF (ionode) THEN ! Store the last wannier product index
           iunwww = find_free_unit()
           OPEN(UNIT=iunwww,FILE=TRIM(tmp_dir)//TRIM(prefix)//'.restart_Www_stat')
           WRITE(iunwww,*) loop_wannier_prod
           CLOSE(iunwww)
        END IF
     END DO

     IF (.NOT. l_no_GW_bare_coulomb) THEN
        DEALLOCATE(overlap, diagonal, e, b, u, h, Ah, pu)
     END IF
     
     close(wannier_product_unit)
     close(wannier_write_product_unit)
     if (debug) close(wwwdebug)
     write(stdout,*) 'files closed' 

  ! Or read from file
  ELSE
     allocate(p_basis(npw,numw_prod))
     CALL diropn( iungprod, 'wiwjwfc_red', npw*2, exst )

     do ii=1,numw_prod
        call davcio(p_basis(:,ii),npw*2,iungprod,ii,-1)
        p_basis(1:npw,ii)=p_basis(1:npw,ii)*fac(1:npw)
     enddo

     close(iungprod)
  END IF
! JDE end

 ! maximum number of non-zero overlap 
   np_max=0
   do iv=1, num_nbndv(ispin)
     np=0
     do jv=1, num_nbndv(ispin)
        if (o_mat(iv,jv)>=s_bse) np=np+1 
     enddo
     if (np>np_max) np_max=np
   enddo

   IF (.NOT. l_no_GW_just_screening) THEN ! JDE
      if(ionode) then 
        iunz = find_free_unit()
        if (ispin==1) open(unit=iunz,file=trim(tmp_dir)//trim(prefix)//'.zbse1',status='unknown',form='unformatted')
        if (ispin==2) open(unit=iunz,file=trim(tmp_dir)//trim(prefix)//'.zbse2',status='unknown',form='unformatted')
        write(iunz) num_nbndv(ispin)
        write(iunz) s_bse 
        write (iunz) np_max
        write (iunz) numw_prod
      endif  
   
      allocate(z(numw_prod))
      z(1:numw_prod)=0.d0
   
   
      do iv=1, num_nbndv(ispin)
        do jv=1, num_nbndv(ispin)
          if (o_mat(jv,iv)>=s_bse) then
           ww_prod(1:fc%nrxxt)= w_wfcs_r(1:fc%nrxxt,iv)* w_wfcs_r(1:fc%nrxxt,jv)
           psic(1:fc%nrxxt)=dcmplx(ww_prod(1:fc%nrxxt),0.d0)
           CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
           ww_prodg(1:fc%npwt) = psic(fc%nlt(1:fc%npwt))
   
           call mergewf(ww_prodg,evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           call splitwf(ww_prodg2,evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
   
           call dgemm('T','N',numw_prod,1,2*npw,2.d0,p_basis,&
                       &2*npw,ww_prodg2,2*npw,0.d0,z,numw_prod)
   
           if(gstart==2) then
             do ii=1,numw_prod
                z(ii)=z(ii)-dble(p_basis(1,ii)*conjg(ww_prodg2(1)))
             enddo
           endif 
           call mp_sum(z,world_comm)
           if(ionode) then
              write(iunz) z
           endif
          endif
         enddo
      enddo
      close(iunz)
   END IF ! JDE

   allocate(iww(np_max)) 

!  in file iwwbse1 we write, for each iv, the set of jv for which
!  o_mat(iv,jv)>=s_bse

   if(ionode) then 
     iuni = find_free_unit()
     if (ispin==1) open(unit=iuni,file=trim(tmp_dir)//trim(prefix)//'.iwwbse1',status='unknown',form='unformatted')
     if (ispin==2) open(unit=iuni,file=trim(tmp_dir)//trim(prefix)//'.iwwbse2',status='unknown',form='unformatted')
     write(iuni) num_nbndv(ispin)
     write(iuni) s_bse 
     write (iuni) np_max
   endif 

  
   do iv=1, num_nbndv(ispin)
     ii=0
     iww(1:np_max)=0
     do jv=1, num_nbndv(ispin)
       if (o_mat(jv,iv)>=s_bse) then
         ii=ii+1
         iww(ii)=jv
       endif 
     enddo
     if(ionode) write(iuni) iww 
   enddo 
   if (ionode) close(iuni) ! JDE fixes a seg fault when fcheck=bounds
   
   if(debug) close(iundebug) 
  
 
  deallocate(iww) 
  IF (.NOT. l_no_GW_just_screening) THEN ! JDE
     deallocate(z)
     deallocate(p_basis)
  END IF ! JDE
  deallocate(ww_prod)
  deallocate(ww_prodg)
  deallocate(ww_prodg2)
   
  deallocate(fac)


   call deallocate_fft_custom(fc)

   deallocate(w_wfcs_t)
   deallocate(w_wfcs_r)
   deallocate(w_wfcs_2)
   deallocate(evc_g)


end subroutine
