!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine transmit(ik, ien)
!
! This subroutine constructs the scattering states
! using the CBS of the left and right tips (saved by compbs)
! and the functions and integrals computed by scatter_forw in
! the scattering region.
!
#include "machine.h"
  use pwcom
  use cond
implicit none

  integer :: ik, ien, n, iorb, iorb1, iorb2, nt, &
             ih, ih1, orbin, ig, ntran, info
  integer, allocatable :: ipiv(:)
  real(kind=DP) :: tk, tj, tij, eev
  real(kind=DP), allocatable :: zps(:,:), eigen(:)
  complex(kind=DP) :: x1, x2
  complex(kind=DP), allocatable :: amat(:,:), vec1(:,:), &
                     tchan(:,:), veceig(:,:)

  eev = earr(ien)
  ntran=4*n2d+norbs+nocrosl+nocrosr

  if(nchanl*nchanr.eq.0) then
    tk = 0.d0
    write(6,'(a20, 2f12.7)') 'E-Ef(ev), TOTAL T = ',   &
                           eev, tk
    return
  endif

  orbin=nocrosl+noinsl

  allocate( ipiv( ntran ) )
  allocate( zps( norbs, norbs ) )
  allocate( amat( ntran, ntran ) )
  allocate( vec1( ntran, nchanl ) )

  allocate( tchan( nchanr, nchanl ) )
  allocate( veceig( nchanl, nchanl ) )
  allocate( eigen( nchanl ) ) 

  amat=(0.d0,0.d0)
!
! To form  zps=zpseu-e*qq
!                                                                            
  do iorb=1, norbs
    do iorb1=1, norbs
      zps(iorb, iorb1)=zpseu(orbin+iorb,orbin+iorb1,iofspin)
    enddo
  enddo
  do iorb=1, norbs-nocrosr
    nt=itnew(orbin+iorb)
    if(tvanp(nt)) then
      do iorb1=1, norbs-nocrosr
        ih=natih(orbin+iorb,2)
        if (natih(orbin+iorb,1).eq.natih(orbin+iorb1,1)) then
          ih1=natih(orbin+iorb1,2)
          zps(iorb,iorb1)=zps(iorb,iorb1)-                &
                                    eryd*qq(ih,ih1,nt)
        endif
      enddo
    endif
  enddo
  do iorb=norbs-nocrosr+1, norbs
    nt=itnew(orbin+iorb)
    if(tvanp(nt)) then
      do iorb1=norbs-nocrosr+1, norbs
        ih=natih(orbin+iorb,2)
        if (natih(orbin+iorb,1).eq.natih(orbin+iorb1,1)) then
          ih1=natih(orbin+iorb1,2)
          zps(iorb,iorb1)=zps(iorb,iorb1)-              &
                                   eryd*qq(ih,ih1,nt)
        endif
      enddo
    endif
  enddo
!
! Compute the part of amat which comes from the matching of
! the wave function on the boundaries
!
!     1) local functions
  do n=1, n2d
    do ig=1, n2d
      amat(ig,n)=fun0(ig,n)
      amat(ig+n2d,n)=fund0(ig,n)
      amat(ig+2*n2d,n)=fun1(ig,n)
      amat(ig+3*n2d,n)=fund1(ig,n)
      amat(ig,n+n2d)=fun0(ig,n2d+n)
      amat(ig+n2d,n+n2d)=fund0(ig,n2d+n)
      amat(ig+2*n2d,n+n2d)=fun1(ig,n2d+n)
      amat(ig+3*n2d,n+n2d)=fund1(ig,n2d+n)
    enddo
  enddo
!     2) nonlocal functions
  do iorb=1, norbs
    do ig=1, n2d
      amat(ig,2*n2d+iorb)=funl0(ig, iorb)
      amat(ig+n2d,2*n2d+iorb)=fundl0(ig, iorb)
      amat(ig+2*n2d,2*n2d+iorb)=funl1(ig, iorb)
      amat(ig+3*n2d,2*n2d+iorb)=fundl1(ig, iorb)
    enddo
  enddo                                                                 
!     4) to add reflection and transmission parts
  do ig=1, n2d
    do n=1, n2d+nocrosl
      amat(ig,2*n2d+norbs+n)=-kfunl(ig,n2d+nocrosl+n)
      amat(ig+n2d,2*n2d+norbs+n)=-kfundl(ig,n2d+nocrosl+n)
    enddo
    do n=1, n2d+nocrosr
      amat(ig+2*n2d,3*n2d+norbs+nocrosl+n)=-kfunr(ig,n)
      amat(ig+3*n2d,3*n2d+norbs+nocrosl+n)=-kfundr(ig,n)
    enddo
  enddo
!
!     3) Part coming from the definion of C_alpha 
!
  do iorb=1, norbs
    do n=1, 2*n2d
      do iorb1=1, norbs
        amat(4*n2d+iorb,n)=amat(4*n2d+iorb,n)-             &
                    zps(iorb,iorb1)*intw1(iorb1,n)
      enddo
    enddo
    do iorb1=1, norbs
      do iorb2=1, norbs
        amat(4*n2d+iorb,2*n2d+iorb1)=amat(4*n2d+iorb,2*n2d+iorb1)-     &
                     zps(iorb,iorb2)*intw2(iorb2, iorb1)
      enddo
    enddo
    amat(4*n2d+iorb,2*n2d+iorb)=amat(4*n2d+iorb,2*n2d+iorb)+(1.d0,0.d0)
  enddo
!     To set up terms depending on kint and kcoef of
!     5) left tip
  do ig=1, nocrosl
    do n=1, n2d+nocrosl
      amat(4*n2d+ig,2*n2d+norbs+n)=-kintl(ig,n2d+nocrosl+n)
      amat(4*n2d+norbs+ig,2*n2d+norbs+n)=-kcoefl(ig,n2d+nocrosl+n)
    enddo
  enddo
!     6) right tip
  do ig=1, nocrosr
    do n=1, n2d+nocrosr
      amat(4*n2d+nocrosl+noinss+ig,3*n2d+norbs+nocrosl+n)=-kintr(ig,n)
      amat(4*n2d+norbs+nocrosl+ig,3*n2d+norbs+nocrosl+n)=-kcoefr(ig,n)
    enddo
  enddo                                                                 
!     7) to match C_alpha for crossing orbitals with the tips ones
  do ig=1, nocrosl
    amat(4*n2d+norbs+ig,2*n2d+ig)=(1.d0,0.d0)
  enddo
  do ig=1, nocrosr
    amat(4*n2d+norbs+nocrosl+ig,2*n2d+nocrosl+noinss+ig)=(1.d0,0.d0)
  enddo
!
! Form the vector of free coefficients
!
  vec1=(0.d0,0.d0)
  do n=1, nchanl
    do ig=1, n2d
      vec1(ig,n)=kfunl(ig,n)
      vec1(n2d+ig,n)=kfundl(ig,n)
    enddo
    do ig=1, nocrosl
      vec1(4*n2d+ig,n)=kintl(ig,n)
      vec1(4*n2d+norbs+ig,n)=kcoefl(ig,n)
    enddo
  enddo
!
! Solve the system on the coefficiens vec1 of scattering states
!
  call ZGESV(ntran, nchanl, amat, ntran, ipiv, vec1,  &
             ntran, info)
!
! transmission coeff. k --> i
!
  do n=1, nchanl
    do ig=1, nchanr
      tchan(ig,n)=vec1(ntran-n2d-nocrosr+ig,n)
    enddo
  enddo
!
! transmission of each band
!
  write(6,*) 'T_ij for propagating states:'
  do n=1, nchanl
    tj=0.d0
    do ig=1, nchanr
      tij=DREAL(tchan(ig,n))**2+DIMAG(tchan(ig,n))**2
      tj=tj+tij
      write(6,'(i5,'' --> '',i5,f12.7)') n, ig, tij
    enddo
    write(6,'(15x,f9.5)') tj
  enddo
!
! eigenchannel decomposition
!
  call eigenchnl(nchanl, nchanr, tchan, veceig, eigen)
  write(6,*) 'Eigenchannel decomposition:'
  tk=0
  do n=1, nchanl
    write(6,'(''@'',i5, 2f9.5)') n, eev, eigen(n)
    do ig=1, nchanl
      tj=DREAL(veceig(ig,n))**2+DIMAG(veceig(ig,n))**2
      write(6,'(20x, f9.5)') tj
    enddo
    tk=tk+eigen(n)
  enddo
  write(6,'(a20, 2f12.7)') 'E-Ef(ev), TOTAL T = ', eev, tk 
!
! To add to the total transmission 
!
  tran_tot(ien) = tran_tot(ien) + wkpt(ik)*tk

  deallocate(ipiv)
  deallocate(zps)
  deallocate(amat)
  deallocate(vec1)

  deallocate(tchan)
  deallocate(veceig)
  deallocate(eigen)

  return
end subroutine transmit                                     

