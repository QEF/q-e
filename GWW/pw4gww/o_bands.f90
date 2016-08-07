!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE o_bands(numv, v_states,numpw,o_basis,ethr,cutoff,ptype)
!this subroutines find the lowest eigenstates of the O matrix
!GAMMA ONLY VERSION
!ONLY FOR NORM_CONSERVING PSEUDOPOTENTIALS

  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps4
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec, igk_k
  USE gvect,                ONLY : g, gstart
  USE wvfct,                ONLY : g2kin, wg, et, nbnd, npwx, npw, current_k
  USE control_flags,        ONLY : max_cg_iter, david
  USE g_psi_mod,            ONLY : h_diag
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE klist,                ONLY : xk
  USE control_flags,        ONLY : isolve
  USE io_files, ONLY : prefix, tmp_dir, diropn
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE fft_base,             ONLY : dfftp, dffts

  implicit none


  INTEGER, EXTERNAL :: find_free_unit
  INTEGER, INTENT(in) :: numv!number of valence states
  REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space
  INTEGER, INTENT(inout) :: numpw!dimension of polarization basis
  COMPLEX(kind=DP), INTENT(inout) :: o_basis(npw,numpw)
  REAL(kind=DP), INTENT(in) :: ethr!threshold for diagonalization
  REAL(kind=DP), INTENT(in) :: cutoff!cutoff for plane waves
  INTEGER, INTENT(in) :: ptype!type of approximation for O operator


  REAL(kind=DP), ALLOCATABLE :: precondition(:)
  REAL(kind=DP), ALLOCATABLE :: o_values(:)!eigenvalues of O operator
  COMPLEX(kind=DP), ALLOCATABLE :: psi_test(:)
  INTEGER :: notconv
  REAL(kind=DP) :: avg_iter
  INTEGER :: ig, iw, il,ii,ip
  REAL(kind=DP)::sca
  REAL(kind=DP), ALLOCATABLE :: hdiag(:)
  INTEGER :: dav_iter

  REAL(kind=DP), ALLOCATABLE :: ovec(:,:)
  INTEGER :: nfound
  INTEGER :: iunwfc
  INTEGER :: l_blk,nbegin,nend,nsize,iunfcw
  LOGICAL :: exst

  INTEGER :: fcw_number!number of "fake conduction" states for O matrix method
  COMPLEX(kind=DP), ALLOCATABLE  :: fcw_state(:,:)! "fake conduction" states for O matrix method
  REAL(kind=DP), ALLOCATABLE :: fcw_mat(:,:)! "fake conduction" matrix

  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  INTEGER :: lwork,info,liwork
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL(kind=DP), ALLOCATABLE :: omat(:,:)
  INTEGER, ALLOCATABLE :: isuppz(:)
  INTEGER :: n_out,num_fc
  INTEGER :: istep
  REAL(kind=DP), PARAMETER :: lambda=0.001d0


  isolve=1
  allocate(o_values(numpw))

  if(ptype==0.or.ptype==1.or.ptype==2) then
     write(stdout,*) 'PTYPE = 1 or 2 NOT IMPLEMENTED YET'
     FLUSH(stdout)
    ! stop
!#ifdef __NOTIMPLEMENTED

     allocate(psi_test(npw))
     allocate(hdiag(npw))
     allocate(precondition(npw))


     write(stdout,*) 'setting preconditioning'
     FLUSH(stdout)

!the following is for calling h_1psi routine
     g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
          ( g(2,igk_k(1:npw,1)) )**2 + &
          ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2
     
     do ig=1,npw
        if(g2kin(ig) <= cutoff) then
           hdiag(ig)=1.d0!/(1.d0+g2kin(ig))
        else
           hdiag(ig)=0.d0
        endif
     enddo
     hdiag=1.d0
  !precondition(1:npw) = 1.D0 + g2kin(1:npw) + &
  !     SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
                  !
     precondition(:)=1.d0
   
!simple steepest descent only for debug purposes
       do istep=1,0!100
           call o_1psi_gamma( numv, v_states, o_basis(:,1), psi_test,.false.,hdiag,ptype,fcw_number,fcw_state,fcw_mat,ethr)
           sca=0.d0
           do ig=1,npw
              sca=sca+2.d0*dble(conjg(o_basis(ig,1))*psi_test(ig))
           enddo
           if(gstart==2) sca=sca-dble(conjg(o_basis(1,1))*psi_test(1))
           call mp_sum(sca,world_comm)
           write(stdout,*) 'Steepest:',sca
           
           o_basis(1:npw,1)= o_basis(1:npw,1)-lambda*psi_test(1:npw)
           
           sca=0.d0
           do ig=1,npw
              sca=sca+2.d0*dble(conjg(o_basis(ig,1))*o_basis(ig,1))
           enddo
           if(gstart==2) sca=sca-dble(conjg(o_basis(1,1))*o_basis(1,1))
           call mp_sum(sca,world_comm)
           sca=dsqrt(sca)
           o_basis(:,1)=o_basis(:,1)/sca


        enddo

        if(isolve==1) then
  
        write(stdout,*) 'call o_rcgdiagg',max_cg_iter
        FLUSH(stdout)
     !precondition(:)=hdiag(:)
        do il=1,1!ATTENZIONE DEBUG
           
           call  o_rcgdiagg( npw, npw, numpw, o_basis, o_values, precondition, &
                ethr, 100, .true., notconv, avg_iter, numv, v_states,hdiag,ptype,fcw_number,fcw_state,fcw_mat)
           
           if(notconv <= 0) exit
        enddo

        do iw=1,numpw
           write(stdout,*) 'Eigen:', iw, o_values(iw)
        enddo
     else
!davidson strategy
     !precondition(:)=hdiag(:)
        do il=1,50
          ! call  o_regterg( npw, npw, numpw, david*numpw, o_basis, ethr, &
          !      .true., gstart, o_values, notconv, .true., dav_iter,&
          !      numv, v_states,hdiag, precondition)
           if(notconv <= 0) exit
        enddo

        do iw=1,numpw
           write(stdout,*) 'Eigen:', iw, o_values(iw)
        enddo
     

     endif
!controlla se sono autovalori
     do iw=1,numpw
        call o_1psi_gamma( numv, v_states, o_basis(:,iw), psi_test,.false.,hdiag,ptype,fcw_number,fcw_state,fcw_mat,ethr)
        sca=0.d0
        do ig=1,npw
           sca=sca+2.d0*dble(conjg(psi_test(ig))*psi_test(ig))
        enddo
        if(gstart==2) sca=sca-dble(conjg(psi_test(1))*psi_test(1))
        call mp_sum(sca,world_comm)
        sca=dsqrt(sca)
        psi_test(:)=psi_test(:)/sca
        sca=0.d0
        do ig=1,npw
           sca=sca+2.d0*dble(conjg(o_basis(ig,iw))*psi_test(ig))
        enddo
        if(gstart==2) sca=sca-dble(conjg(o_basis(1,iw))*psi_test(1))
        call mp_sum(sca,world_comm)
        write(stdout,*) 'Eig prod:',iw,sca

     enddo

     FLUSH(stdout)

     deallocate(precondition)
     deallocate(psi_test)
   !  deallocate (becp%r)
    ! call deallocate_bec_type ( becp)
     deallocate(hdiag)
!#endif __NOTIMPLEMENTED
  else if(ptype==3 .or. ptype==4) then
     !read from file
     if(ionode) then
        iunfcw = find_free_unit()
        open(unit=iunfcw,file=trim(tmp_dir)//trim(prefix)//'.nfcws',status='old')
        read(iunfcw,*) fcw_number
        close(iunfcw)
     endif
 
     call mp_bcast(fcw_number, ionode_id,world_comm)

     if(numpw>fcw_number) then
        numpw=fcw_number
        write(stdout,*) 'Set polarizability basis dimension:', numpw
        FLUSH(stdout)
     endif

     allocate(ovec(fcw_number,numpw))

     write(stdout,*) 'ATT1', fcw_number
     FLUSH(stdout)

     l_blk= (fcw_number)/nproc
     if(l_blk*nproc < (fcw_number)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > fcw_number) nend=fcw_number
     nsize=nend-nbegin+1

     write(stdout,*) 'ATT2', fcw_number
     FLUSH(stdout)

     if(nsize>0) then
        allocate(fcw_mat(fcw_number,nsize))
     else
        allocate(fcw_mat(fcw_number,1))
     endif
     allocate(fcw_state(npw,fcw_number))

     write(stdout,*) 'ATT3', fcw_number
     FLUSH(stdout)

     iunfcw = find_free_unit()
     CALL diropn( iunfcw, 'fcw', npw*2, exst )
     do ii=1,fcw_number
        CALL davcio(fcw_state(:,ii), 2*npw,iunfcw,ii,-1)
     enddo
     close(iunfcw)

     write(stdout,*) 'ATT4', fcw_number
     FLUSH(stdout)


     CALL diropn( iunfcw, 'fmat',fcw_number, exst )
     do ii=1,nsize
        CALL davcio(fcw_mat(:,ii), fcw_number,iunfcw,ii,-1)
     enddo
     close(iunfcw)

     write(stdout,*) 'ATT5', fcw_number
     FLUSH(stdout)


     if(ptype==3) then
        call  diago_cg(fcw_number,fcw_mat,1000,numpw,o_values,ovec,0.d0,1d-8,nfound,.true.)
     else
        allocate(omat(fcw_number,fcw_number))
        omat(:,:)=0.d0
        if(nsize>0) omat(1:fcw_number,nbegin:nend)=fcw_mat(1:fcw_number,1:nsize)
        do iw=1,fcw_number
           call mp_sum(omat(:,iw),world_comm)
        enddo
             
        if(ionode) then
           allocate(isuppz(fcw_number))
           allocate(work(1),iwork(1))
           call DSYEVR('V','I','U',fcw_number,omat,fcw_number,0.d0,0.d0,&
        &fcw_number-numpw+1,fcw_number,1d-8,n_out,o_values,ovec,fcw_number,isuppz,work,-1,iwork,-1,info)
           lwork=work(1)
           liwork=iwork(1)
           deallocate(work,iwork)
           allocate(work(lwork))
           allocate(iwork(liwork))
         
           call DSYEVR('V','I','U',fcw_number,omat,fcw_number,0.d0,0.d0,&
     &fcw_number-numpw+1,fcw_number,1d-8,n_out,o_values,ovec,fcw_number,isuppz,work,lwork,iwork,liwork,info)
           if(info/=0) then
              write(stdout,*) 'ROUTINE fake_conduction_wannier, INFO:', info
              stop
           endif
         
           deallocate(isuppz)
           deallocate(work,iwork)
        else
           o_values(:)=0.d0
           ovec(:,:)=0.d0
        endif
        do iw=1,numpw
           call mp_sum(ovec(:,iw),world_comm)
        enddo
        call mp_sum(o_values(:),world_comm)
        deallocate(omat)
        do iw=1,numpw
           write(stdout,*) 'POLARIZABILITY eigen:', iw, o_values(iw)
        enddo
        FLUSH(stdout)
     endif
     call dgemm('N','N',2*npw,numpw,fcw_number,1.d0,fcw_state,2*npw,ovec,fcw_number,0.d0,o_basis,2*npw)

     deallocate(ovec)
     deallocate(fcw_mat,fcw_state)
  else if(ptype==5) then
!just real plane waves
     g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
          ( g(2,igk_k(1:npw,1)) )**2 + &
          ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2
     
     num_fc=0
     do ig=1,npw
        if(g2kin(ig) <= cutoff) num_fc=num_fc+1
     enddo
     call mp_sum(num_fc,world_comm)
     num_fc=(num_fc-1)*2!G=0 excluded
  
     o_basis(:,:)=(0.d0,0.d0)
     write(stdout,*) 'Number of G states', num_fc
     if(num_fc>numpw) then 
        write(stdout,*) 'numw_prod too small:', num_fc
        FLUSH(stdout)
        stop
     endif
     numpw=num_fc
     ii=0
     do ip=0,nproc-1
        if(mpime==ip) then
           do ig=gstart,npw
              if(g2kin(ig) <= cutoff) then
                 ii=ii+1
                 o_basis(ig,ii)=cmplx(dsqrt(0.5d0),0.d0)
                 ii=ii+1
                 o_basis(ig,ii)=cmplx(0.d0,dsqrt(0.5d0))
              endif
           enddo
        else
           ii=0
        endif
        call mp_sum(ii,world_comm)
     enddo
     if(ii/=num_fc) then
        write(stdout,*) 'ERRORE G STATES',ii
        FLUSH(stdout)
        stop
        return
     endif
  

  endif
  if(allocated(o_values) )deallocate(o_values)
  return

end SUBROUTINE o_bands


subroutine o_extra_pw( p_basis, numwp, numwp_max,cutoff)

  !this subroutines add to the polarizability basis at the end, plane waves (sin and cos) up to the specified cutoff


  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec, igk_k
  USE gvect,                ONLY : g, gstart
  USE wvfct,                ONLY : g2kin, wg, et, nbnd, npwx, npw, current_k
  USE control_flags,        ONLY : max_cg_iter, david
  USE g_psi_mod,            ONLY : h_diag
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE klist,                ONLY : xk
  USE mp_world,             ONLY : mpime, nproc, world_comm
  USE wannier_gw,           ONLY : optimal_options

  implicit none

  INTEGER, INTENT(inout) :: numwp!dimension of polarization basis
  COMPLEX(kind=DP), INTENT(inout) :: p_basis(npw,numwp_max)
  REAL(kind=DP), INTENT(in) :: cutoff
  INTEGER, INTENT(in) :: numwp_max!max dimension of polarizability basis
  INTEGER :: num_fc, ii, ip, ig, info,numwp2
  TYPE(optimal_options) :: options


  !just real plane waves                        
  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2

  num_fc=0
  do ig=1,npw
     if(g2kin(ig) <= cutoff) num_fc=num_fc+1
  enddo
  call mp_sum(num_fc,world_comm)
  num_fc=(num_fc-1)*2!G=0 excluded           
  p_basis(:,numwp+1:numwp+num_fc)=(0.d0,0.d0)
  write(stdout,*) 'Number of G states added to the polarizability basis', num_fc
 
  ii=numwp
  do ip=0,nproc-1
     if(mpime==ip) then
        do ig=gstart,npw
           if(g2kin(ig) <= cutoff) then
              ii=ii+1
              p_basis(ig,ii)=cmplx(dsqrt(0.5d0),0.d0)
              ii=ii+1
              p_basis(ig,ii)=cmplx(0.d0,dsqrt(0.5d0))
           endif
        enddo
     else
        ii=0
     endif
     call mp_sum(ii,world_comm)
  enddo
  if(ii/=num_fc+numwp) then
     write(stdout,*) 'ERRORE G STATES',ii
     FLUSH(stdout)
     stop
     return
  endif
  numwp=numwp+num_fc
  write(stdout,*) 'UPDATED DIMESION OF POLARIZABILITY BASIS: ', numwp
!now re-orthonormalize
  options%l_complete=.true.
  options%idiago=0
  options%ithres=0
  options%thres=0.d0
  call optimal_driver(numwp,p_basis,npw,options,numwp2, info)
  write(stdout,*) 'UPDATED DIMESION OF POLARIZABILITY BASIS: ', numwp
  if(info/=0) then
     write(stdout,*) 'PROBLEM WITH OPTIMAL_DRIVER'
     FLUSH(stdout)
     stop
     return
  endif


  return

end subroutine o_extra_pw


subroutine  update_numwp(numwp, cutoff)
!this subroutine adds to numwp the number of plane-waves differnt from 0 till the specified cutoff

  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec, igk_k
  USE gvect,                ONLY : g, gstart
  USE wvfct,                ONLY : g2kin, wg, et, nbnd, npwx, npw, current_k
  USE control_flags,        ONLY : max_cg_iter, david
  USE g_psi_mod,            ONLY : h_diag
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE klist,                ONLY : xk
  USE mp_world,             ONLY : mpime, nproc, world_comm




  implicit none

  INTEGER, INTENT(inout) :: numwp!dimension of polarization basis
  REAL(kind=DP), INTENT(in) :: cutoff

  INTEGER :: num_fc, ig

  !just real plane waves 

  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2

  num_fc=0
  do ig=1,npw
     if(g2kin(ig) <= cutoff) num_fc=num_fc+1
  enddo
  call mp_sum(num_fc,world_comm)
  num_fc=(num_fc-1)*2!G=0 excluded              

  numwp=numwp+num_fc

  return

end subroutine update_numwp
