!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine compbs(lright, zin, zfin, nocros, norbnow, orbin,     &
                  nchan, kval, kfun, kfund, kint, kcoef)
!
! Using the basis functions obtained by scatter_forw it computes
! the complex band structure (CBS) of the tip ( zin<z<zfin ). 
! Some variables needed for wave-function matching in transmission
! calculation are constructed and saved.
!
#include "machine.h"
  use pwcom
  use cond
  implicit none
  integer ::   &
     nocros,   & ! number of orbitals crossing the boundary
     noins,    & ! number of interior orbitals 
     norbnow,  & ! total number of orbitals
     orbin,    & ! number of 1-st orbital in full list 
     lright      ! 1/0 if it is right/left tip    
  integer :: ik, i, j, kin, kfin, ig, info, lam, n, iorb, iorb1, &
             iorb2, aorb, borb, nt, startk, lastk, nchan, nb, ih, &
             ih1 
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: raux, zin, zfin, dz, DDOT
  real(kind=DP), allocatable :: zps(:,:)
  complex(kind=DP), parameter :: cim=(0.d0,1.d0) 
  complex(kind=DP) :: x1, x2, x3, x4, y1, y2, y3, y4,                  &
            kval(2*(n2d+nocros)), kfun(n2d,2*(n2d+nocros)),            &
            kint(nocros,2*(n2d+nocros)), kcoef(nocros,2*(n2d+nocros)), &
            kfund(n2d,2*(n2d+nocros))   
  complex(kind=DP), allocatable :: amat(:,:), bmat(:,:), vec(:,:) 

  noins = norbnow-2*nocros

  allocate( amat( (2*n2d+norbnow), (2*n2d+norbnow) ) )
  allocate( bmat( (2*n2d+norbnow), (2*n2d+norbnow) ) )
  allocate( vec( (2*n2d+norbnow), 2*(n2d+nocros) ) )
  allocate( zps( norbnow, norbnow ) )

!
! To find indeces of initial and final slab
!
  do ik=1, nrz
    if (z(ik).le.zin+eps)  kin=ik
    if (z(ik).le.zfin-eps) kfin=ik
  enddo
  dz=z(kin+1)-z(kin)

  amat=(0.d0,0.d0)
  bmat=(0.d0,0.d0)
!
! zps=zpseu-e*qq for US-PP and zps=zpseu for norm-conserv. PP
!
  orbin=orbin-1
  do iorb=1, norbnow
    do iorb1=1, norbnow 
      zps(iorb, iorb1)=zpseu(orbin+iorb,orbin+iorb1,iofspin)
    enddo
  enddo 
  do iorb=1, nocros+noins
    nt=itnew(orbin+iorb) 
    if(tvanp(nt)) then
      do iorb1=1, nocros+noins
        ih=natih(orbin+iorb,2)
        if (natih(orbin+iorb,1).eq.natih(orbin+iorb1,1)) then
          ih1=natih(orbin+iorb1,2)               
          zps(iorb,iorb1)=zps(iorb,iorb1)-eryd*qq(ih,ih1,nt)
        endif 
      enddo
    endif
  enddo 
  do iorb=1, nocros
    do iorb1=1, nocros
      zps(iorb+noins+nocros,iorb1+noins+nocros)=zps(iorb,iorb1)       
    enddo
  enddo 

!
! Forming the matrices A and B for generalized eigenvalue problem

!
!   1
!

  do n=1, 2*n2d
    do ig=1, n2d
      amat(ig, n)=fun1(ig, n)
      amat(ig+n2d,n)=fund1(ig, n)
      bmat(ig,n)=fun0(ig, n)
      bmat(ig+n2d,n)=fund0(ig, n)
     enddo
  enddo
!
!   2 
!
  do iorb=1, norbnow
    do ig=1, n2d
      amat(ig, 2*n2d+iorb)=funl1(ig, iorb)
      amat(n2d+ig, 2*n2d+iorb)=fundl1(ig, iorb)      
      bmat(ig, 2*n2d+iorb)=funl0(ig, iorb)
      bmat(n2d+ig, 2*n2d+iorb)=fundl0(ig, iorb) 
    enddo
  enddo 
!
!  3
!
  do iorb=1, norbnow
    aorb=iorb
    borb=iorb
    if (iorb.le.nocros) aorb=iorb+noins+nocros 
    if (iorb.gt.nocros) borb=iorb-noins-nocros 
    do n=1, 2*n2d
      do iorb1=1, norbnow 
        amat(2*n2d+iorb,n)=amat(2*n2d+iorb,n)+                  &
                           zps(iorb1,aorb)*intw1(iorb1,n)         
        if (borb.gt.0) bmat(2*n2d+iorb,n)=                      &
         bmat(2*n2d+iorb,n)-zps(iorb1,borb)*intw1(iorb1,n)
      enddo
    enddo
  enddo
!
!  4
!
  do iorb=1, nocros
    do iorb1=1, norbnow
      do iorb2=1, norbnow
        bmat(2*n2d+iorb,2*n2d+iorb1)=bmat(2*n2d+iorb,2*n2d+iorb1) &
               -zps(iorb2,iorb)*intw2(iorb2, iorb1)
      enddo
      bmat(2*n2d+iorb+noins+nocros,2*n2d+iorb1)=                  &
                              bmat(2*n2d+iorb,2*n2d+iorb1)
    enddo
    bmat(2*n2d+iorb,2*n2d+iorb)=                                  &
                         bmat(2*n2d+iorb,2*n2d+iorb)+(1.d0,0.d0)
  enddo
!
!   5 
!
  do iorb=1, norbnow
    aorb=iorb
    if (iorb.le.nocros) aorb=iorb+noins+nocros
    do iorb1=1, norbnow
      do iorb2=1, norbnow
        amat(2*n2d+iorb,2*n2d+iorb1)=amat(2*n2d+iorb,2*n2d+iorb1)+  &
                 zps(iorb2,aorb)*intw2(iorb2, iorb1)
      enddo  
    enddo
    if (aorb.eq.iorb) amat(2*n2d+iorb,2*n2d+iorb)=                  &
                        amat(2*n2d+iorb,2*n2d+iorb)-(1.d0,0.d0)             
  enddo  

!
! To reduce matrices and solve GEP A X = c B X; X = {a_n, a_\alpha}
! 
  call compbs_2(nocros, norbnow, n2d, 2*(n2d+nocros),               &
                amat, bmat, vec, kval, llapack) 

!
! To normalize (over XY plane) all the states
!

  kfun=(0.d0,0.d0)
  kfund=(0.d0,0.d0)
  do ik=1, 2*(n2d+nocros)
    do ig=1, n2d
      do j=1,  2*n2d+norbnow
        kfun(ig,ik)= kfun(ig,ik)+amat(ig,j)*vec(j,ik)
        kfund(ig,ik)= kfund(ig,ik)+amat(n2d+ig,j)*vec(j,ik) 
      enddo
    enddo
  enddo
  do ik=1, 2*(n2d+nocros)
    raux=DDOT(2*n2d,kfun(1,ik),1,kfun(1,ik),1)*sarea
    raux=1.d0/sqrt(raux)
    call DSCAL(2*(2*n2d+norbnow),raux,vec(1,ik),1)
    call DSCAL(2*n2d,raux,kfun(1,ik),1)
    call DSCAL(2*n2d,raux,kfund(1,ik),1)
  enddo

!
! To find k-vector and the current of Bloch states
!

  call kbloch (2*(n2d+nocros), kval)

  call jbloch(2*(n2d+nocros), n2d, norbf, norbnow, nocros,  &
              kfun, kfund, vec, kval, intw1, intw2, sarea, nchan)

!
! To save band structure result
!
  kfun=(0.d0,0.d0)
  kfund=(0.d0,0.d0)
  kint=(0.d0,0.d0)
!
! To account for the case of the right lead  
!
  if (lright.eq.1) then
    do i=1, 2*n2d
      do j=1, 2*n2d+norbnow
        amat(i,j)=bmat(i,j) 
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+nocros
      do j=1, 2*n2d+norbnow 
        amat(i,j)=-bmat(i+nocros+noins,j) 
      enddo
    enddo  
  endif
!
! psi_k and psi'_k on the scattering region boundary
!
  do ik=1, 2*(n2d+nocros)
    do ig=1, n2d
      do j=1,  2*n2d+norbnow   
        kfun(ig,ik)= kfun(ig,ik)+amat(ig,j)*vec(j,ik)
        kfund(ig,ik)= kfund(ig,ik)+amat(n2d+ig,j)*vec(j,ik)
      enddo
    enddo
  enddo

!
! kint(iorb, ik)=\sum_{iorb1} D_{iorb,iorb1} 
!            \int_{cell} W_iorb1^* psi_ik )
! for the orbitals crossing the boundary
!                              
  do ik=1,  2*(n2d+nocros)
    do iorb=1, nocros
      do j=1, 2*n2d+norbnow 
        kint(iorb,ik)=kint(iorb,ik)+amat(2*n2d+iorb,j)*vec(j,ik)
      enddo
    enddo
  enddo 
!
! a_iorb = kcoef(iorb,ik) = \sum_{iorb1} D_{iorb,iorb1}
!           \int_{all space} W_iorb1^* psi_ik )
! for the orbitals crossing the boundary
!
  do ik=1, 2*(n2d+nocros)
    do iorb=1, nocros
      if (lright.eq.1) then 
        kcoef(iorb,ik)=vec(2*n2d+iorb,ik)
      else
        kcoef(iorb,ik)=vec(2*n2d+nocros+noins+iorb,ik)
      endif   
    enddo
  enddo 

!
! to set up B.S. for the right lead in the case of identical tips
!
  if(lright.eq.0.and.ikind.eq.1) then
    nchanr=nchan 
    call DCOPY(2*(n2d+nocros), kval, 1, kvalr, 1)
    kfunr=(0.d0,0.d0)
    kfundr=(0.d0,0.d0)
    kintr=(0.d0,0.d0)

    do i=1, 2*n2d
      do j=1, 2*n2d+norbnow
        amat(i,j)=bmat(i,j)
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+nocros
      do j=1, 2*n2d+norbnow
        amat(i,j)=-bmat(i+nocros+noins,j)
      enddo
    enddo               

    do ik=1, 2*(n2d+nocros)
      do ig=1, n2d
        do j=1,  2*n2d+norbnow
          kfunr(ig,ik)= kfunr(ig,ik)+amat(ig,j)*vec(j,ik)
          kfundr(ig,ik)= kfundr(ig,ik)+amat(n2d+ig,j)*vec(j,ik)
        enddo
      enddo
    enddo        
    do ik=1,  2*(n2d+nocros)
      do iorb=1, nocros
        do j=1, 2*n2d+norbnow
          kintr(iorb,ik)=kintr(iorb,ik)+amat(2*n2d+iorb,j)*vec(j,ik)
        enddo
      enddo
    enddo
    do ik=1, 2*(n2d+nocros)
      do iorb=1, nocros
        kcoefr(iorb,ik)=vec(2*n2d+iorb,ik)
      enddo
    enddo                     

  endif

  deallocate(amat)
  deallocate(bmat)
  deallocate(vec)  
  deallocate(zps)

  return
end subroutine compbs  
