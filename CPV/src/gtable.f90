!
! Copyright (C) 2002-2010 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine gtable( ipol, ctable)

  ! this subroutine prepares the correspondence array to
  ! compute the operator exp(iG_ipol.r)

  !   ctable : output correspondence table
  !    in (ig,1) correspondence for g+1
  !    in (ig,2) correspondence for (-g)+1
  !    we use the rule: if non point ngw+1
  !    if found positive = normal
  !    negative = conjugate
  !   ipol   : input polarization direction
  !            a orthorombic primitive cell is supposed
  use kinds, only: dp
  use gvecw, only: ngw  
  use gvect, only: mill
  use mp, only: mp_sum
  use io_global, only: ionode, stdout
  use mp_global, only: intra_bgrp_comm 

  implicit none
  integer :: ipol, ctable(ngw,2)
  !local variables
  integer :: i,j,k, ig, jg
  logical :: found
  real(dp) :: test

  test=0.d0
  do ig=1,ngw!loop on g vectors
     ! first +g
     i = mill(1,ig)
     j = mill(2,ig)
     k = mill(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1
     
     found = .false.
     
     do jg=1,ngw
        if(mill(1,jg).eq.i .and. mill(2,jg).eq.j .and. mill(3,jg).eq.k) then
           found=.true.
           ctable(ig,1)=jg
        endif
     enddo

     if(.not. found) then        
        do jg=1,ngw
           if(-mill(1,jg).eq.i .and. -mill(2,jg).eq.j .and. -mill(3,jg).eq.k) then
              found=.true.
              ctable(ig,1)=-jg
           endif
        enddo
        if(.not. found)  then
           ctable(ig,1)= ngw+1
           test=test+1.d0
        endif
     endif

     ! now -g
     i = -mill(1,ig)
     j = -mill(2,ig)
     k = -mill(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1
     
     found = .false.

     do jg=1,ngw
        if (-mill(1,jg).eq.i .and. -mill(2,jg).eq.j .and. -mill(3,jg).eq.k)then
           found=.true.
           ctable(ig,2)=-jg
        endif
     enddo
     
     if(.not.found) then
        do jg=1,ngw
           if(mill(1,jg).eq.i .and. mill(2,jg).eq.j .and. mill(3,jg).eq.k)then
              found=.true.
              ctable(ig,2)=jg
           endif
        enddo
        if(.not.found) then
           ctable(ig,2)=ngw+1
           test=test+1.d0
        endif
     endif
  enddo

  call mp_sum(test, intra_bgrp_comm)
  if(ionode) write(stdout,*) '#not found, gtable: ', test

  return
end subroutine gtable

subroutine gtablein( ipol, ctabin)
  
  ! this subroutine prepare the inverse correspondence array to
  ! compute the operator exp(iG_ipol.r)
  
  !   ctabin(ngw,2) : output correspondence table
  !   if negative to take complex conjugate, 1 g'+1, 2 g' -1
  !   if not found = ngw+1
  !   ipol   : input polarization direction
  !            a orthorombic primitive cell is supposed

  use kinds, only: dp
  use gvecw, only: ngw  
  use gvect, only: mill
  use mp, only: mp_sum
  use io_global, only: ionode, stdout
  use mp_global, only: intra_bgrp_comm

  implicit none

  integer :: ipol, ctabin(ngw,2)

  !local variables
  integer :: i,j,k, ig, jg
  logical :: found
  real(dp) :: test

  test=0.d0
  
  do ig=1,ngw!loop on g vectors
     i = mill(1,ig)
     j = mill(2,ig)
     k = mill(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1
     found = .false.

     do jg=1,ngw
        if(i.eq.mill(1,jg).and. j.eq.mill(2,jg) .and. k.eq.mill(3,jg))then
           found = .true.
           ctabin(ig,1)=jg
        else if(i.eq.-mill(1,jg).and. j.eq.-mill(2,jg) .and. k.eq.-mill(3,jg))then
           found=.true.
           ctabin(ig,1)=-jg
        endif
     enddo
     if(.not.found) then
        ctabin(ig,1)=ngw+1
        test=test+1
     endif
  enddo
  
  do ig=1,ngw!loop on g vectors
     i = mill(1,ig)
     j = mill(2,ig)
     k = mill(3,ig)
     if(ipol.eq.1) i=i-1
     if(ipol.eq.2) j=j-1
     if(ipol.eq.3) k=k-1
     found = .false.

     do jg=1,ngw
        if(i.eq.mill(1,jg).and. j.eq.mill(2,jg) .and. k.eq.mill(3,jg))then
           found = .true.
           ctabin(ig,2)=jg
        else if(i.eq.-mill(1,jg).and. j.eq.-mill(2,jg) .and. k.eq.-mill(3,jg))then
           found=.true.
           ctabin(ig,2)=-jg
        endif
     enddo
     if(.not.found) then
        ctabin(ig,2)=ngw+1
        test=test+1
     endif
  enddo

  call mp_sum(test, intra_bgrp_comm)
  if(ionode) write(stdout,*) '#not found, gtabin: ', test

  return
  
end subroutine gtablein



subroutine find_whose_is_g
!this subroutine set the correspondence G-->Proc

  USE gvecw,              ONLY : ngw, ngw_g
  USE gvect,              ONLY : ig_l2g, mill_g, mill
  USE mp,                 ONLY : mp_sum
  USE io_global,          ONLY : stdout
  USE mp_global,          ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
  USE efield_module,      ONLY : whose_is_g

  implicit none

  INTEGER :: ig

  whose_is_g(:)=0


  do ig=1,ngw
     if(ig_l2g(ig) > ngw_g) then
        write(stdout,*) 'find_whose_is_g: too large'
        stop
     endif
     whose_is_g(ig_l2g(ig))=me_bgrp+1
  enddo
  call mp_sum(whose_is_g,intra_bgrp_comm)
  whose_is_g(:)=whose_is_g(:)-1
  ! mill_g is used in gtable_missing and re-initialized here
  ! workaround by PG to avoid a large array like mill_g allocated all the time

  allocate ( mill_g(3,ngw_g) )
  do ig=1,ngw
     mill_g(:,ig_l2g(ig)) = mill(:,ig)
  end do
  call mp_sum(mill_g,intra_bgrp_comm)

  
return
end subroutine find_whose_is_g


subroutine gtable_missing

  USE efield_module,      ONLY : ctable_missing_1, ctable_missing_2,  &
                                 whose_is_g,n_g_missing_p, &
                                 ctable_missing_rev_1,ctable_missing_rev_2 
  USE gvecw,              ONLY : ngw, ngw_g
  USE gvect, ONLY : ig_l2g, mill_g, mill, gstart
  USE mp,                 ONLY : mp_sum, mp_max, mp_alltoall
  USE io_global,          ONLY : stdout
  USE mp_global,          ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
  USE parallel_include

  implicit none

  INTEGER :: ipol, i,j,k,ig,igg, nfound_max, ip
  LOGICAL :: found
  INTEGER :: nfound_proc(nproc_bgrp,2)
  INTEGER, ALLOCATABLE :: igg_found(:,:,:), ig_send(:,:,:), igg_found_snd(:,:,:)
  INTEGER, ALLOCATABLE :: igg_found_rcv(:,:,:)
  INTEGER :: ierr,sndint,rcvint

  allocate( igg_found(ngw_g,2,nproc_bgrp), ig_send(ngw_g,2,nproc_bgrp) )
  do ipol=1,2
    
     nfound_max=0
     nfound_proc(:,:)=0
     ig_send(:,:,:)=0

     do ig=1,ngw!loop on g vectors
        ! first +g
        i = mill(1,ig)
        j = mill(2,ig)
        k = mill(3,ig)
        if(ipol.eq.1) i=i+1
        if(ipol.eq.2) j=j+1
        if(ipol.eq.3) k=k+1
        do igg=1,ngw_g
           if( i==mill_g(1,igg) .and. j==mill_g(2,igg) .and. k==mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,1)=nfound_proc(whose_is_g(igg)+1,1)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=igg
              endif
              
           else if( i==-mill_g(1,igg) .and. j==-mill_g(2,igg) .and. k==-mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,1)=nfound_proc(whose_is_g(igg)+1,1)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=-igg
              endif
           endif
              
        enddo
     enddo
  
     do ig=gstart,ngw!loop on g vectors
        ! first +g
        i = -mill(1,ig)
        j = -mill(2,ig)
        k = -mill(3,ig)
        if(ipol.eq.1) i=i+1
        if(ipol.eq.2) j=j+1
        if(ipol.eq.3) k=k+1
        do igg=1,ngw_g
           if( i==mill_g(1,igg) .and. j==mill_g(2,igg) .and. k==mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,2)=nfound_proc(whose_is_g(igg)+1,2)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=igg
              endif

           else if( i==-mill_g(1,igg) .and. j==-mill_g(2,igg) .and. k==-mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                  nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,2)=nfound_proc(whose_is_g(igg)+1,2)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=-igg
              endif
           endif

        enddo
     enddo
      
! determine the largest nfound for processor and set it as dimension
! for ctable_missing and ctable_missing_rev
! copy ig_send to ctable_missing

     call mp_sum(nfound_max, intra_bgrp_comm)
     write(stdout,*) 'Additional found:', nfound_max

     n_g_missing_p(ipol)=maxval(nfound_proc(:,:))
     
     call mp_max(n_g_missing_p(ipol), intra_bgrp_comm)
         

    if(ipol==1) then
        allocate(ctable_missing_1(n_g_missing_p(ipol),2,nproc_bgrp))
        ctable_missing_1(:,:,:)=0
        do ip=1,nproc_bgrp
           ctable_missing_1(1:nfound_proc(ip,1),1,ip)=ig_send(1:nfound_proc(ip,1),1,ip)
           ctable_missing_1(1:nfound_proc(ip,2),2,ip)=ig_send(1:nfound_proc(ip,2),2,ip)
        enddo
     else
        allocate(ctable_missing_2(n_g_missing_p(ipol),2,nproc_bgrp))
        ctable_missing_2(:,:,:)=0
        do ip=1,nproc_bgrp
           ctable_missing_2(1:nfound_proc(ip,1),1,ip)=ig_send(1:nfound_proc(ip,1),1,ip)
           ctable_missing_2(1:nfound_proc(ip,2),2,ip)=ig_send(1:nfound_proc(ip,2),2,ip)
        enddo
     endif
    

!mpi all to all for igg_found

     allocate(igg_found_snd(n_g_missing_p(ipol),2,nproc_bgrp))
     allocate(igg_found_rcv(n_g_missing_p(ipol),2,nproc_bgrp))
     igg_found_snd(:,:,:)=0
     do ip=1,nproc_bgrp
        igg_found_snd(1:nfound_proc(ip,1),1,ip)=igg_found(1:nfound_proc(ip,1),1,ip)
        igg_found_snd(1:nfound_proc(ip,2),2,ip)=igg_found(1:nfound_proc(ip,2),2,ip)
     enddo

  
     call mp_alltoall( igg_found_snd, igg_found_rcv, intra_bgrp_comm )
   
     if(ipol==1) then
        allocate(ctable_missing_rev_1(n_g_missing_p(ipol),2,nproc_bgrp))
        ctable_missing_rev_1(:,:,:)=0
     else
        allocate(ctable_missing_rev_2(n_g_missing_p(ipol),2,nproc_bgrp))
        ctable_missing_rev_2(:,:,:)=0
     endif


  
     nfound_max=0

     do ip=1,nproc_bgrp
        do igg=1, n_g_missing_p(ipol)
           if(igg_found_rcv(igg,1,ip) /= 0 ) then
              found=.false.
              do ig=1,ngw
                 if(igg_found_rcv(igg,1,ip)>0) then
                    if(ig_l2g(ig)==igg_found_rcv(igg,1,ip)) then
                       nfound_max=nfound_max+1
                     if(ipol==1) then
                          ctable_missing_rev_1(igg,1,ip)=ig
                       else
                          ctable_missing_rev_2(igg,1,ip)=ig
                       endif
                       found=.true.
                    endif
                 else
                    if(ig_l2g(ig)==-igg_found_rcv(igg,1,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctable_missing_rev_1(igg,1,ip)=-ig
                       else
                          ctable_missing_rev_2(igg,1,ip)=-ig
                       endif
                       found=.true.
                    endif
                 endif
              enddo
              if(.not.found) write(stdout,*) 'NOT FOUND:', igg_found_rcv(igg,1,ip)
           endif
        enddo
         do igg=1, n_g_missing_p(ipol)
           if(igg_found_rcv(igg,2,ip) /= 0 ) then
              found=.false.
              do ig=1,ngw
                 if(igg_found_rcv(igg,2,ip)>0) then
                    if(ig_l2g(ig)==igg_found_rcv(igg,2,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctable_missing_rev_1(igg,2,ip)=ig
                       else
                          ctable_missing_rev_2(igg,2,ip)=ig
                       endif
                       found=.true.
                    endif
                 else
                    if(ig_l2g(ig)==-igg_found_rcv(igg,2,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctable_missing_rev_1(igg,2,ip)=-ig
                       else
                          ctable_missing_rev_2(igg,2,ip)=-ig
                       endif
                       found=.true.
                    endif
                 endif
              enddo
              if(.not.found) write(stdout,*) 'NOT FOUND:', igg_found_rcv(igg,2,ip)
           endif
        enddo

     enddo
     call mp_sum(nfound_max, intra_bgrp_comm)
     !write(stdout,*) 'Found check', nfound_max
     deallocate(igg_found_snd,igg_found_rcv)
  enddo



  deallocate(igg_found, ig_send)
return

end subroutine gtable_missing




subroutine gtable_missing_inv

  USE efield_module, ONLY : ctabin_missing_1,ctabin_missing_2, whose_is_g,n_g_missing_m,&
                          &      ctabin_missing_rev_1,ctabin_missing_rev_2 
  USE gvecw,              ONLY : ngw, ngw_g
  USE gvect, ONLY : ig_l2g, mill_g, mill, gstart
  USE mp,                 ONLY : mp_sum, mp_max, mp_alltoall
  USE io_global,          ONLY : stdout
  USE mp_global,          ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
  USE parallel_include


  implicit none

  INTEGER :: ipol, i,j,k,ig,igg, nfound_max, ip
  LOGICAL :: found
  INTEGER :: nfound_proc(nproc_bgrp,2)
  INTEGER, ALLOCATABLE :: igg_found(:,:,:), ig_send(:,:,:), igg_found_snd(:,:,:)
  INTEGER, ALLOCATABLE :: igg_found_rcv(:,:,:)
  INTEGER :: ierr,sndint,rcvint



  allocate( igg_found(ngw_g,2,nproc_bgrp), ig_send(ngw_g,2,nproc_bgrp))
  do ipol=1,2



     nfound_max=0
     nfound_proc(:,:)=0
     ig_send(:,:,:)=0

     do ig=1,ngw!loop on g vectors
        ! first +g
        i = mill(1,ig)
        j = mill(2,ig)
        k = mill(3,ig)
        if(ipol.eq.1) i=i+1
        if(ipol.eq.2) j=j+1
        if(ipol.eq.3) k=k+1
        do igg=1,ngw_g
           if( i==mill_g(1,igg) .and. j==mill_g(2,igg) .and. k==mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,1)=nfound_proc(whose_is_g(igg)+1,1)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=igg
              endif
              
           else if( i==-mill_g(1,igg) .and. j==-mill_g(2,igg) .and. k==-mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,1)=nfound_proc(whose_is_g(igg)+1,1)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,1),1,whose_is_g(igg)+1)=-igg
              endif
           endif
              
        enddo
     enddo
  
     do ig=1,ngw!loop on g vectors
        ! first +g
        i = mill(1,ig)
        j = mill(2,ig)
        k = mill(3,ig)
        if(ipol.eq.1) i=i-1
        if(ipol.eq.2) j=j-1
        if(ipol.eq.3) k=k-1
        do igg=1,ngw_g
           if( i==mill_g(1,igg) .and. j==mill_g(2,igg) .and. k==mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                 nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,2)=nfound_proc(whose_is_g(igg)+1,2)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=igg
              endif

           else if( i==-mill_g(1,igg) .and. j==-mill_g(2,igg) .and. k==-mill_g(3,igg)) then
              if(whose_is_g(igg) /= -1 .and. whose_is_g(igg) /= me_bgrp) then
                  nfound_max=nfound_max+1
                 nfound_proc(whose_is_g(igg)+1,2)=nfound_proc(whose_is_g(igg)+1,2)+1
                 ig_send(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=ig
                 igg_found(nfound_proc(whose_is_g(igg)+1,2),2,whose_is_g(igg)+1)=-igg
              endif
           endif

        enddo
     enddo
        

!determine the largest nfound for processor and set it as dimensione for ctabin_missing and ctabin_missing_rev
!copy ig_send to ctabin_missing

     call mp_sum(nfound_max, intra_bgrp_comm)
     write(stdout,*) 'Additional found:', nfound_max

     
     n_g_missing_m(ipol)=maxval(nfound_proc(:,:))
     call mp_max(n_g_missing_m(ipol), intra_bgrp_comm)
          

    if(ipol==1) then
        allocate(ctabin_missing_1(n_g_missing_m(ipol),2,nproc_bgrp))
        ctabin_missing_1(:,:,:)=0
        do ip=1,nproc_bgrp
           ctabin_missing_1(1:nfound_proc(ip,1),1,ip)=ig_send(1:nfound_proc(ip,1),1,ip)
           ctabin_missing_1(1:nfound_proc(ip,2),2,ip)=ig_send(1:nfound_proc(ip,2),2,ip)
        enddo
     else
        allocate(ctabin_missing_2(n_g_missing_m(ipol),2,nproc_bgrp))
        ctabin_missing_2(:,:,:)=0
        do ip=1,nproc_bgrp
           ctabin_missing_2(1:nfound_proc(ip,1),1,ip)=ig_send(1:nfound_proc(ip,1),1,ip)
           ctabin_missing_2(1:nfound_proc(ip,2),2,ip)=ig_send(1:nfound_proc(ip,2),2,ip)
        enddo
     endif
    

!mpi all to all for igg_found

     allocate(igg_found_snd(n_g_missing_m(ipol),2,nproc_bgrp))
     allocate(igg_found_rcv(n_g_missing_m(ipol),2,nproc_bgrp))
     igg_found_snd(:,:,:)=0
     do ip=1,nproc_bgrp
        igg_found_snd(1:nfound_proc(ip,1),1,ip)=igg_found(1:nfound_proc(ip,1),1,ip)
        igg_found_snd(1:nfound_proc(ip,2),2,ip)=igg_found(1:nfound_proc(ip,2),2,ip)
     enddo

   
     CALL mp_alltoall( igg_found_snd, igg_found_rcv, intra_bgrp_comm )

     if(ipol==1) then
        allocate(ctabin_missing_rev_1(n_g_missing_m(ipol),2,nproc_bgrp))
        ctabin_missing_rev_1(:,:,:)=0
     else
        allocate(ctabin_missing_rev_2(n_g_missing_m(ipol),2,nproc_bgrp))
        ctabin_missing_rev_2(:,:,:)=0
     endif


   
     nfound_max=0

     do ip=1,nproc_bgrp
        do igg=1, n_g_missing_m(ipol)
           if(igg_found_rcv(igg,1,ip) /= 0 ) then
              found=.false.
              do ig=1,ngw
                 if(igg_found_rcv(igg,1,ip)>0) then
                    if(ig_l2g(ig)==igg_found_rcv(igg,1,ip)) then
                       nfound_max=nfound_max+1
                     if(ipol==1) then
                          ctabin_missing_rev_1(igg,1,ip)=ig
                       else
                          ctabin_missing_rev_2(igg,1,ip)=ig
                       endif
                       found=.true.
                    endif
                 else
                    if(ig_l2g(ig)==-igg_found_rcv(igg,1,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctabin_missing_rev_1(igg,1,ip)=-ig
                       else
                          ctabin_missing_rev_2(igg,1,ip)=-ig
                       endif
                       found=.true.
                    endif
                 endif
              enddo
              if(.not.found) write(stdout,*) 'NOT FOUND:', igg_found_rcv(igg,1,ip)
           endif
        enddo
         do igg=1, n_g_missing_m(ipol)
           if(igg_found_rcv(igg,2,ip) /= 0 ) then
              found=.false.
              do ig=1,ngw
                 if(igg_found_rcv(igg,2,ip)>0) then
                    if(ig_l2g(ig)==igg_found_rcv(igg,2,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctabin_missing_rev_1(igg,2,ip)=ig
                       else
                          ctabin_missing_rev_2(igg,2,ip)=ig
                       endif
                       found=.true.
                    endif
                 else
                    if(ig_l2g(ig)==-igg_found_rcv(igg,2,ip)) then
                       nfound_max=nfound_max+1
                       if(ipol==1) then
                          ctabin_missing_rev_1(igg,2,ip)=-ig
                       else
                          ctabin_missing_rev_2(igg,2,ip)=-ig
                       endif
                       found=.true.
                    endif
                 endif
              enddo
              if(.not.found) write(stdout,*) 'NOT FOUND:', igg_found_rcv(igg,2,ip)
           endif
        enddo

     enddo
     call mp_sum(nfound_max, intra_bgrp_comm)
     !write(stdout,*) 'Found check', nfound_max
     deallocate(igg_found_snd,igg_found_rcv)
  enddo

  deallocate(igg_found, ig_send)
  ! workaround by PG to avoid a large array like mill_g allocated all the time
  deallocate ( mill_g )

return

end subroutine gtable_missing_inv
