

!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
!
subroutine compbs(lright, zin, zfin, nocros, norbnow, orbin,     &
                  nchan, kval, kfun, kfund, kint, kcoef)
!
! Using the basis functions obtained by scatter_forw it computes
! the complex band structure (CBS) of the tip ( zin<z<zfin ). 
! Some variables needed for wave-function matching in transmission
! calculation are constructed and saved.
!
#include "f_defs.h"
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param, ONLY: tvanp
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
             ih1, ij, is, js 
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: raux, zin, zfin, dz, DDOT
  real(kind=DP), allocatable :: zps(:,:)
  complex(kind=DP), parameter :: cim=(0.d0,1.d0) 
  complex(kind=DP) :: x1, x2, x3, x4, y1, y2, y3, y4,                  &
            kval(2*(n2d+npol*nocros)), kfun(n2d,2*(n2d+npol*nocros)),  &
            kint(nocros*npol,2*(n2d+npol*nocros)),  &
            kcoef(nocros*npol,2*(n2d+npol*nocros)), &
            kfund(n2d,2*(n2d+npol*nocros))   
  complex(kind=DP), allocatable :: amat(:,:), bmat(:,:), vec(:,:), &
            zps_nc(:,:), aux(:,:) 
  complex(kind=DP), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)

  call start_clock('compbs')
  noins = norbnow-2*nocros

  allocate( amat( (2*n2d+npol*norbnow), (2*n2d+npol*norbnow) ) )
  allocate( bmat( (2*n2d+npol*norbnow), (2*n2d+npol*norbnow) ) )
  allocate( vec( (2*n2d+npol*norbnow), 2*(n2d+npol*nocros) ) )
  allocate( zps( norbnow, norbnow ) )
  allocate( aux( n2d, 2*n2d+npol*norbnow))
  if (noncolin) allocate( zps_nc( norbnow*npol, norbnow*npol ) )

!
! To find indices of initial and final slab
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
      if (noncolin) then
        ij=0
        do is=1,npol
          do js=1,npol
            ij=ij+1
            zps_nc(npol*(iorb-1)+is, npol*(iorb1-1)+js)=  &
                 zpseu_nc(orbin+iorb, orbin+iorb1,ij)
          enddo
        enddo
      else
        zps(iorb, iorb1)=zpseu(orbin+iorb,orbin+iorb1,iofspin)
      endif
    enddo
  enddo 
  do iorb=1, nocros+noins
    nt=itnew(orbin+iorb) 
    if(tvanp(nt).or.lspinorb) then
      ih=natih(orbin+iorb,2)
      do iorb1=1, nocros+noins
        if (natih(orbin+iorb,1).eq.natih(orbin+iorb1,1)) then
          ih1=natih(orbin+iorb1,2)               
          if (noncolin) then
            ij=0
            do is=1,npol
              do js=1,npol
                ij=ij+1
                if (lspinorb) then
                  zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)= &
                        zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                                -eryd*qq_so(ih,ih1,ij,nt)
                else
                  if (is.eq.js) &
                      zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)= &
                          zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                                 -eryd*qq(ih,ih1,nt)
                endif
              enddo
            enddo
          else
            zps(iorb,iorb1)=zps(iorb,iorb1)-eryd*qq(ih,ih1,nt)
          endif
        endif 
      enddo
    endif
  enddo 

  do iorb=1, nocros*npol
    do iorb1=1, nocros*npol
      if (noncolin) then
         zps_nc(iorb+npol*(noins+nocros), iorb1+npol*(noins+nocros))= &
                   zps_nc(iorb,iorb1)       
      else
        zps(iorb+noins+nocros,iorb1+noins+nocros)=zps(iorb,iorb1)       
      endif
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
  do iorb=1, norbnow*npol
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
  do iorb=1, norbnow*npol
    aorb=iorb
    borb=iorb
    if (iorb.le.npol*nocros) aorb=iorb+npol*(noins+nocros) 
    if (iorb.gt.npol*nocros) borb=iorb-npol*(noins+nocros) 
    do n=1, 2*n2d
      do iorb1=1, norbnow*npol 
        if (noncolin) then
           amat(2*n2d+iorb,n)=amat(2*n2d+iorb,n)+                  &
                             zps_nc(aorb, iorb1)*intw1(iorb1,n)         
           if (borb.gt.0) bmat(2*n2d+iorb,n)=                      &
              bmat(2*n2d+iorb,n)-zps_nc(borb,iorb1)*intw1(iorb1,n)
        else
           amat(2*n2d+iorb,n)=amat(2*n2d+iorb,n)+                  &
                             zps(aorb,iorb1)*intw1(iorb1,n)         
           if (borb.gt.0) bmat(2*n2d+iorb,n)=                      &
              bmat(2*n2d+iorb,n)-zps(borb,iorb1)*intw1(iorb1,n)
        endif
      enddo
    enddo
  enddo
!
!  4
!
  do iorb=1, nocros*npol
    do iorb1=1, norbnow*npol
      do iorb2=1, norbnow*npol
        if (noncolin) then
           bmat(2*n2d+iorb,2*n2d+iorb1)=bmat(2*n2d+iorb,2*n2d+iorb1) &
                 -zps_nc(iorb,iorb2)*intw2(iorb2, iorb1)
        else
           bmat(2*n2d+iorb,2*n2d+iorb1)=bmat(2*n2d+iorb,2*n2d+iorb1) &
                 -zps(iorb,iorb2)*intw2(iorb2, iorb1)
        endif
      enddo
      bmat(2*n2d+iorb+npol*(noins+nocros),2*n2d+iorb1)=                  &
                              bmat(2*n2d+iorb,2*n2d+iorb1)
    enddo
    bmat(2*n2d+iorb,2*n2d+iorb)=                                  &
                         bmat(2*n2d+iorb,2*n2d+iorb)+(1.d0,0.d0)
  enddo
!
!   5 
!
  do iorb=1, norbnow*npol
    aorb=iorb
    if (iorb.le.npol*nocros) aorb=iorb+npol*(noins+nocros)
    do iorb1=1, norbnow*npol
      do iorb2=1, norbnow*npol
        if (noncolin) then
           amat(2*n2d+iorb,2*n2d+iorb1)=amat(2*n2d+iorb,2*n2d+iorb1)+  &
                   zps_nc(aorb,iorb2)*intw2(iorb2, iorb1)
        else
           amat(2*n2d+iorb,2*n2d+iorb1)=amat(2*n2d+iorb,2*n2d+iorb1)+  &
                   zps(aorb,iorb2)*intw2(iorb2, iorb1)
        endif
      enddo  
    enddo
    if (aorb.eq.iorb) amat(2*n2d+iorb,2*n2d+iorb)=                  &
                        amat(2*n2d+iorb,2*n2d+iorb)-(1.d0,0.d0)             
  enddo  

!
! To reduce matrices and solve GEP A X = c B X; X = {a_n, a_\alpha}
! 
  call compbs_2(npol*nocros, npol*norbnow, n2d, 2*(n2d+npol*nocros),   &
                amat, bmat, vec, kval, llapack) 

!
! To normalize (over XY plane) all the states
!

!  kfun=(0.d0,0.d0)
!  kfund=(0.d0,0.d0)
!  do ik=1, 2*(n2d+npol*nocros)
!    do ig=1, n2d
!      do j=1,  2*n2d+npol*norbnow
!        kfun(ig,ik)= kfun(ig,ik)+amat(ig,j)*vec(j,ik)
!        kfund(ig,ik)= kfund(ig,ik)+amat(n2d+ig,j)*vec(j,ik) 
!      enddo
!    enddo
!  enddo
!
   call ZGEMM('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norbnow,   &
               one, amat, 2*n2d+npol*norbnow, vec, 2*n2d+npol*norbnow,  &
               zero, kfun, n2d)
   do ig=1,n2d
      do ik=1, 2*n2d+npol*norbnow
         aux(ig,ik)=amat(n2d+ig,ik)
      enddo
   enddo
   call ZGEMM('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norbnow,   &
               one, aux, n2d, vec, 2*n2d+npol*norbnow, zero, kfund, n2d)


  do ik=1, 2*(n2d+npol*nocros)
    raux=DDOT(2*n2d,kfun(1,ik),1,kfun(1,ik),1)*sarea
    raux=1.d0/sqrt(raux)
    call DSCAL(2*(2*n2d+npol*norbnow),raux,vec(1,ik),1)
    call DSCAL(2*n2d,raux,kfun(1,ik),1)
    call DSCAL(2*n2d,raux,kfund(1,ik),1)
  enddo

!
! To find k-vector and the current of Bloch states
!

  call kbloch (2*(n2d+npol*nocros), kval)

  call jbloch(2*(n2d+npol*nocros), n2d, norbf, norbnow, nocros,  &
              kfun, kfund, vec, kval, intw1, intw2, sarea, nchan, npol)
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
      do j=1, 2*n2d+npol*norbnow
        amat(i,j)=bmat(i,j) 
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+npol*nocros
      do j=1, 2*n2d+npol*norbnow 
        amat(i,j)=-bmat(i+npol*(nocros+noins),j) 
      enddo
    enddo  
  endif
!
! psi_k and psi'_k on the scattering region boundary
!
!  do ik=1, 2*(n2d+npol*nocros)
!    do ig=1, n2d
!      do j=1,  2*n2d+npol*norbnow   
!        kfun(ig,ik)= kfun(ig,ik)+amat(ig,j)*vec(j,ik)
!        kfund(ig,ik)= kfund(ig,ik)+amat(n2d+ig,j)*vec(j,ik)
!      enddo
!    enddo
!  enddo

   call ZGEMM('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norbnow,   &
               one, amat, 2*n2d+npol*norbnow, vec, 2*n2d+npol*norbnow,  &
               zero, kfun, n2d)
   do ig=1,n2d
      do ik=1, 2*n2d+npol*norbnow
         aux(ig,ik)=amat(n2d+ig,ik)
      enddo
   enddo
   call ZGEMM('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norbnow,   &
               one, aux, n2d, vec, 2*n2d+npol*norbnow, zero, kfund, n2d)

!  do j=1,2*nchan
!     write(6,'("------------------------------",i5,2f15.6)') j, kval(j)
!     do ik=1,n2d
!        if (ik.eq.n2d/2+1.and.noncolin) write(6,'("-------")')
!        write(6,'(i5,f15.7)') ik, abs(kfun(ik,j))
!     enddo
!  enddo
!  stop

!
! kint(iorb, ik)=\sum_{iorb1} D_{iorb,iorb1} 
!            \int_{cell} W_iorb1^* psi_ik )
! for the orbitals crossing the boundary
!                              
  do ik=1,  2*(n2d+npol*nocros)
    do iorb=1, nocros*npol
      do j=1, 2*n2d+npol*norbnow 
        kint(iorb,ik)=kint(iorb,ik)+amat(2*n2d+iorb,j)*vec(j,ik)
      enddo
    enddo
  enddo 
!
! a_iorb = kcoef(iorb,ik) = \sum_{iorb1} D_{iorb,iorb1}
!           \int_{all space} W_iorb1^* psi_ik )
! for the orbitals crossing the boundary
!
  do ik=1, 2*(n2d+npol*nocros)
    do iorb=1, nocros*npol
      if (lright.eq.1) then 
        kcoef(iorb,ik)=vec(2*n2d+iorb,ik)
      else
        kcoef(iorb,ik)=vec(2*n2d+npol*(nocros+noins)+iorb,ik)
      endif   
    enddo
  enddo 

!
! to set up B.S. for the right lead in the case of identical tips
!
  if(lright.eq.0.and.ikind.eq.1) then
    nchanr=nchan 
    call DCOPY(2*(n2d+npol*nocros), kval, 1, kvalr, 1)
    kfunr=(0.d0,0.d0)
    kfundr=(0.d0,0.d0)
    kintr=(0.d0,0.d0)

    do i=1, 2*n2d
      do j=1, 2*n2d+npol*norbnow
        amat(i,j)=bmat(i,j)
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+npol*nocros
      do j=1, 2*n2d+npol*norbnow
        amat(i,j)=-bmat(i+npol*(nocros+noins),j)
      enddo
    enddo               

    do ik=1, 2*(n2d+npol*nocros)
      do ig=1, n2d
        do j=1,  2*n2d+npol*norbnow
          kfunr(ig,ik)= kfunr(ig,ik)+amat(ig,j)*vec(j,ik)
          kfundr(ig,ik)= kfundr(ig,ik)+amat(n2d+ig,j)*vec(j,ik)
        enddo
      enddo
    enddo        
    do ik=1,  2*(n2d+npol*nocros)
      do iorb=1, nocros*npol
        do j=1, 2*n2d+npol*norbnow
          kintr(iorb,ik)=kintr(iorb,ik)+amat(2*n2d+iorb,j)*vec(j,ik)
        enddo
      enddo
    enddo
    do ik=1, 2*(n2d+npol*nocros)
      do iorb=1, nocros*npol
        kcoefr(iorb,ik)=vec(2*n2d+iorb,ik)
      enddo
    enddo                     

  endif

  deallocate(amat)
  deallocate(bmat)
  deallocate(vec)
  deallocate(zps)
  deallocate(aux)
  if (noncolin) deallocate(zps_nc)
  call stop_clock('compbs')

  return
end subroutine compbs  
