!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine gtable( ipol, ctable)

!respect to vectorial (serial) program I changed ngwx to ngw :-)

! this subroutine prepare the correspondence array to
! compute the operator exp(iG_ipol.r)

!   ctable : output coorespondence table
!    in (ig,1) correspondence for g+1
!    in (ig,2) correspondence for (-g)+1
!    we use the rule: if non point ngw+1
!    if found positive = normal
!    negative = conjugate
!   ipol   : input polarization direction
!            a orthorombic primitive cell is supposed


  use gvecw, only: ngw  
  use reciprocal_vectors, only: mill_l


  implicit none
  integer ctable(ngw,2), ipol

!local variables
  integer i,j,k, ig, jg
  logical found
  real(8) test

  test=0.

  do ig=1,ngw!loop on g vectors
! first +g
     i = mill_l(1,ig)
     j = mill_l(2,ig)
     k = mill_l(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1

     found = .false.
!#ifdef NEC
!     *vdir nodep(ctable)
!#endif 
     do jg=1,ngw
        if(mill_l(1,jg).eq.i .and. mill_l(2,jg).eq.j .and. mill_l(3,jg).eq.k) then
           found=.true.
           ctable(ig,1)=jg
!               jg=ngw
        endif
     enddo
     if(.not. found) then
!#ifdef NEC
!        *vdir nodep
!#endif 
        do jg=1,ngw
           if(-mill_l(1,jg).eq.i .and. -mill_l(2,jg).eq.j .and. -mill_l(3,jg).eq.k) then
              found=.true.
              ctable(ig,1)=-jg
!     jg=ngw 
           endif
        enddo
        if(.not. found)  then
           ctable(ig,1)= ngw+1
           test=test+1.
        endif
     endif


    
! now -g
     i = -mill_l(1,ig)
     j = -mill_l(2,ig)
     k = -mill_l(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1

     found = .false.

!#ifdef NEC
!     *vdir nodep
!#endif 
     do jg=1,ngw
        if (-mill_l(1,jg).eq.i .and. -mill_l(2,jg).eq.j .and. -mill_l(3,jg).eq.k)then
           found=.true.
           ctable(ig,2)=-jg
!               jg=ngw
        endif
     enddo
     if(.not.found) then

!#ifdef NEC
!        *vdir nodep
!#endif 
        do jg=1,ngw
           if(mill_l(1,jg).eq.i .and. mill_l(2,jg).eq.j .and. mill_l(3,jg).eq.k)then
              found=.true.
              ctable(ig,2)=jg
!                  jg=ngw
           endif
        enddo
        if(.not.found) then
           ctable(ig,2)=ngw+1
           test=test+1.
        endif
     endif
  enddo


#ifdef __PARA
  call reduce(1,test)
#endif
  write(6,*) '#not found, gtable: ', test

  return
end subroutine gtable




! this subroutine prepare the inverse correspondence array to
! compute the operator exp(iG_ipol.r)

subroutine gtablein( ipol, ctabin)

!   ctabin(ngw,2) : output coorespondence table
!   if negative to take complex conjugate, 1 g'+1, 2 g' -1
!   if not found = ngw+1
!   ipol   : input polarization direction
!            a orthorombic primitive cell is supposed

  
  use gvecw, only: ngw  
  use reciprocal_vectors, only: mill_l

  implicit none

      
  integer ctabin(ngw,2), ipol

!local variables
  integer i,j,k, ig, jg
  logical found
  real(8) test

  test=0.

  do ig=1,ngw!loop on g vectors
     i = mill_l(1,ig)
     j = mill_l(2,ig)
     k = mill_l(3,ig)
     if(ipol.eq.1) i=i+1
     if(ipol.eq.2) j=j+1
     if(ipol.eq.3) k=k+1
     found = .false.

!#ifdef NEC
!     *vdir nodep
!#endif 
     do jg=1,ngw
        if(i.eq.mill_l(1,jg).and. j.eq.mill_l(2,jg) .and. k.eq.mill_l(3,jg))then
           found = .true.
           ctabin(ig,1)=jg
!               jg=ngw
        else if(i.eq.-mill_l(1,jg).and. j.eq.-mill_l(2,jg) .and. k.eq.-mill_l(3,jg))then
           found=.true.
           ctabin(ig,1)=-jg
!               jg=ngw
        endif
     enddo
     if(.not.found) then
        ctabin(ig,1)=ngw+1
        test=test+1
     endif
  enddo

  do ig=1,ngw!loop on g vectors
     i = mill_l(1,ig)
     j = mill_l(2,ig)
     k = mill_l(3,ig)
     if(ipol.eq.1) i=i-1
     if(ipol.eq.2) j=j-1
     if(ipol.eq.3) k=k-1
     found = .false.

!#ifdef NEC
!     *vdir nodep
!#endif 
     do jg=1,ngw
        if(i.eq.mill_l(1,jg).and. j.eq.mill_l(2,jg) .and. k.eq.mill_l(3,jg))then
           found = .true.
           ctabin(ig,2)=jg
!               jg=ngw
        else if(i.eq.-mill_l(1,jg).and. j.eq.-mill_l(2,jg) .and. k.eq.-mill_l(3,jg))then
           found=.true.
           ctabin(ig,2)=-jg
!               jg=ngw
        endif
     enddo
     if(.not.found) then
        ctabin(ig,2)=ngw+1
        test=test+1
     endif



  enddo
#ifdef __PARA
  call reduce(1,test)
#endif
  write(6,*) '#not found, gtabin: ', test

  return

end subroutine gtablein







