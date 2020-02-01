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
subroutine jbloch (nst, n2d, norbf, norb, nocros, kfun, kfund, &
                     vec, kval, c1, c2, nchan, npol)
!
! This routine computes the current I carrying by the Bloch state
! so that <psi_k|I|psi_k'> = \delta_{kk'} I_k and does some
! rearrangements.
!
  USE kinds, only : DP
  USE noncollin_module, only : noncolin
  USE cond, only : sarea
  implicit none

  real(DP), parameter :: eps=1.d-4
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  integer ::                                     &
      n2d,     &  ! 2D-dimension
      noins,   &  ! number of interior orbitals
      nocros,  &  ! number of the orbitals crossing the boundary
      norb, &  ! total number of orbitals =noins+2*nocros
      norbf,   &  ! max number of orbitals
      npol,    &  ! number of wave-functions components
      nst,     &  ! number of Bloch states =2*(n2d+nocros)
      nchan,   &  ! number of propagating channels
      info, j, k, n, iorb, il, ir, in
  real(DP) ::    &
      k1
  complex(DP) ::                           &
      kfun(n2d, nst),                           & ! phi_k(d)
      kfund(n2d, nst),                          & ! phi'_k(d)
      vec(2*n2d+npol*norb, nst),             & ! exp. coeff. for phi_k
      kval(nst),                                & ! k of phi_k
      c1(norbf*npol,2*n2d), c2(norbf*npol,norbf*npol),   & ! nonlocal integrals
      z, zdotc
  integer, allocatable ::   &
      ncond(:)    ! channel --> Bloch state correspondence
  real(DP), allocatable :: ej(:), kcur(:)
  complex(DP), allocatable :: kval1(:), vec1(:,:), kcoef(:,:), &
                               kcuroff(:,:), valj(:,:)

  noins = norb-2*nocros
!
! To compute the number of channels and ncond(+-,>,<) --> (...)
!
  allocate( ncond( nst ) )
  nchan=0
  do k=1, nst
    if (abs(AIMAG(kval(k))).le.eps) then
      nchan=nchan+1
      ncond(nchan)=k
    endif
  enddo
  nchan=nchan/2

  allocate( kcoef( npol*nocros, 2*nchan ) )
  allocate( kcuroff( 2*nchan, 2*nchan ) )
  allocate( ej( 2*nchan ) )
  allocate( valj( 2*nchan, 2*nchan ) )
  allocate( kval1( nst ) )
  allocate( kcur( nst ) )
  allocate( vec1( (2*n2d+npol*norb), nst ) )

  kcur=0.d0
  kcoef=(0.d0,0.d0)

  do k=1, 2*nchan
    ir=ncond(k)
    do iorb=1, nocros*npol
      do j=1, 2*n2d
        kcoef(iorb,k)=kcoef(iorb,k)+       &
                      vec(j,ir)*c1(npol*(nocros+noins)+iorb,j)
      enddo
      do j=1, norb*npol
        kcoef(iorb,k)=kcoef(iorb,k)+       &
           vec(2*n2d+j,ir)*c2(npol*(nocros+noins)+iorb,j)
      enddo
    enddo
  enddo

!
! Calculation of the current matrix <phi_k|I|phi_n>
!

  do k=1, 2*nchan
    do n=1, 2*nchan
      ir=ncond(k)
      il=ncond(n)
      z=zdotc(n2d,kfun(1,ir),1,kfund(1,il),1)
      kcuroff(k,n)=-cim*          &
              (z-zdotc(n2d,kfund(1,ir),1,kfun(1,il),1))*sarea
!     ---------------------------------------------
      do iorb=1, nocros*npol
        kcuroff(k,n)=kcuroff(k,n)-cim*(                               &
            CONJG(vec(2*n2d+npol*(nocros+noins)+iorb,ir))*kcoef(iorb,n)-   &
            vec(2*n2d+npol*(nocros+noins)+iorb,il)*CONJG(kcoef(iorb,k)))
      enddo
!   WRITE( 6,'(2i5, 2f12.6)') k,n, DBLE(kcuroff(k,n)),AIMAG(kcuroff(k,n))
    enddo
  enddo

!
! Diagonalizing of the current matrix
!
  if(nchan.gt.0) then
   info=-1
   call hev_ab(2*nchan, kcuroff, 2*nchan, ej, valj, 0.d0, 0.d0, info)
  endif
!
! Right ordering (+, >, -, <)
!
!-- propagating modes
    ir=0
    il=nst/2
    do in=1, 2*nchan
      if (ej(in).gt.0.d0) then
        ir=ir+1
        do n=1, 2*n2d+norb*npol
          vec1(n,ir)=0.d0
          do j=1, 2*nchan
            vec1(n,ir)=vec1(n,ir)+vec(n,ncond(j))*valj(j,in)
          enddo
        enddo
        kcur(ir)=ej(in)
        do n=1, 2*nchan
          k1= DBLE(valj(n,in))**2+AIMAG(valj(n,in))**2
          if(abs(k1).gt.eps) kval1(ir)=kval(ncond(n))
        enddo
      else
        il=il+1
        do n=1, 2*n2d+npol*norb
          vec1(n,il)=0.d0
          do j=1, 2*nchan
            vec1(n,il)=vec1(n,il)+vec(n,ncond(j))*valj(j,in)
          enddo
        enddo
        kcur(il)=ej(in)
        do n=1, 2*nchan
          k1= DBLE(valj(n,in))**2+AIMAG(valj(n,in))**2
          if(abs(k1).gt.eps) kval1(il)=kval(ncond(n))
        enddo
      endif
    enddo

!--  decaying states
    do in=1, nst
      if (AIMAG(kval(in)).gt.eps) then
        ir=ir+1
        kval1(ir)=kval(in)
        call dcopy(2*(2*n2d+npol*norb),vec(1,in),1,vec1(1,ir),1)
      endif
      if (-AIMAG(kval(in)).gt.eps) then
        il=il+1
        kval1(il)=kval(in)
        call dcopy(2*(2*n2d+npol*norb),vec(1,in),1,vec1(1,il),1)
      endif
    enddo

!
! Normalization to the unit current
!
  do k=1, nchan
    k1=1.d0/abs(kcur(k))
    call dscal(2*(2*n2d+npol*norb),sqrt(k1),vec1(1,k),1)
  enddo
  do k=nst/2+1, nst/2+nchan
    k1=1.d0/abs(kcur(k))
    call dscal(2*(2*n2d+npol*norb),sqrt(k1),vec1(1,k),1)
  enddo

!
! Final rearrangement to ascending |Re(k)|
!
    do k=1, nst
      ncond(k)=k
    enddo
    do k=1, nchan
      kval(k)=1.d2
      do n=1, nchan
        if(abs( DBLE(kval1(n))).le.abs( DBLE(kval(k)))) then
          kval(k)=kval1(n)
          ncond(k)=n
        endif
      enddo
      kval1(ncond(k))=1.d3
    enddo
    do k=nchan+1, nst/2
      kval(k)=kval1(k)
    enddo
    do k=nst/2+1, nst/2+nchan
      kval(k)=1.d2
      do n=nst/2+1, nst/2+nchan
        if(abs( DBLE(kval1(n))).le.abs( DBLE(kval(k)))) then
          kval(k)=kval1(n)
          ncond(k)=n
        endif
      enddo
      kval1(ncond(k))=1.d3
    enddo
    do k=nst/2+nchan+1, nst
      kval(k)=kval1(k)
    enddo
    do k=1, nst
      call dcopy(2*(2*n2d+npol*norb),vec1(1,ncond(k)),1,vec(1,k),1)
    enddo

  deallocate(kcoef)
  deallocate(kval1)
  deallocate(kcur)
  deallocate(vec1)
  deallocate(ncond)
  deallocate(kcuroff)
  deallocate(ej)
  deallocate(valj)

  return
end subroutine jbloch
