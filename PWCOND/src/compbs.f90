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
subroutine compbs(lleft, nocros, norb, nchan, kval, kfun,  &
                  kfund, kint, kcoef, ikk, ien)
!
! Using the basis functions obtained by scatter_forw it computes
! the complex band structure (CBS) of the lead.
! Some variables needed for wave-function matching in transmission
! calculation are constructed and saved.
!
  USE constants, ONLY : tpi
  USE noncollin_module, ONLY : noncolin, npol
  USE spin_orb,  ONLY : lspinorb
  USE lsda_mod,  ONLY : nspin
  USE cond
  USE cell_base, ONLY : alat, at, omega
  USE ions_base, ONLY : nat, ityp

  implicit none
  integer ::   &
     nocros,   & ! number of orbitals crossing the boundary
     noins,    & ! number of interior orbitals
     norb,     & ! total number of orbitals
     lleft       ! 1/0 if it is left/right tip
  integer :: ik, ikk, i, j, ig, n, iorb, iorb1, &
             iorb2, aorb, borb, nchan,     &
             ij, is, js, ichan, ien
  REAL(DP), PARAMETER :: eps=1.d-8
  REAL(DP) :: raux, ddot
  REAL(DP), ALLOCATABLE :: zpseu(:,:,:), zps(:,:)
  COMPLEX(DP), PARAMETER :: cim=(0.d0,1.d0)
  COMPLEX(DP) :: x1,          &
            kval(2*(n2d+npol*nocros)), kfun(n2d,2*(n2d+npol*nocros)),  &
            kint(nocros*npol,2*(n2d+npol*nocros)),                     &
            kcoef(nocros*npol,2*(n2d+npol*nocros)),                    &
            kfund(n2d,2*(n2d+npol*nocros))
  COMPLEX(DP), ALLOCATABLE :: amat(:,:), bmat(:,:), vec(:,:), &
                              zpseu_nc(:,:,:,:), zps_nc(:,:), &
                              aux(:,:), veceig(:,:), korb(:,:)
  COMPLEX(DP), PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)


  call start_clock('compbs')

  noins = norb-2*nocros
  IF (norb>0) THEN
     if(lleft.eq.1) then
       if (noncolin) then
         allocate( zpseu_nc(2,norb,norb,nspin) )
         zpseu_nc = zpseul_nc
       else
         allocate( zpseu(2,norb,norb) )
         zpseu = zpseul
       endif
     else
       if (noncolin) then
          allocate( zpseu_nc(2,norb,norb,nspin) )
          zpseu_nc = zpseur_nc
       else
          allocate( zpseu(2,norb,norb) )
          zpseu = zpseur
       endif
     endif
     if (noncolin) then
        allocate( zps_nc( norb*npol, norb*npol ) )
     else
        allocate( zps( norb, norb ) )
     endif
  END IF

  allocate( amat( (2*n2d+npol*norb), (2*n2d+npol*norb) ) )
  allocate( bmat( (2*n2d+npol*norb), (2*n2d+npol*norb) ) )
  allocate( vec( (2*n2d+npol*norb), 2*(n2d+npol*nocros) ) )
  allocate( aux( n2d, 2*n2d+npol*norb))
  IF (lorb) allocate( korb(npol*(nocros+noins),2*(n2d+npol*nocros)) )

  amat=(0.d0,0.d0)
  bmat=(0.d0,0.d0)

!
! zps=zpseu-e*qq for US-PP and zps=zpseu for norm-conserv. PP
!
  do iorb=1, norb
    do iorb1=1, norb
      if (noncolin) then
        ij=0
        do is=1,npol
          do js=1,npol
            ij=ij+1
            zps_nc(npol*(iorb-1)+is, npol*(iorb1-1)+js)=         &
                zpseu_nc(1,iorb,iorb1,ij)
            if (lspinorb) then
              zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)=        &
                    zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                            -eryd*zpseu_nc(2,iorb,iorb1,ij)
            else
              if (is.eq.js)                                      &
                  zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)=    &
                      zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                             -eryd*zpseu_nc(2,iorb,iorb1,ij)
            endif
          enddo
        enddo
      else
        zps(iorb,iorb1)=zpseu(1,iorb,iorb1)-eryd*zpseu(2,iorb,iorb1)
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
  do iorb=1, norb*npol
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
  do iorb=1, norb*npol
    aorb=iorb
    borb=iorb
    if (iorb.le.npol*nocros) aorb=iorb+npol*(noins+nocros)
    if (iorb.gt.npol*nocros) borb=iorb-npol*(noins+nocros)
    do n=1, 2*n2d
      do iorb1=1, norb*npol
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
    do iorb1=1, norb*npol
      do iorb2=1, norb*npol
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
  do iorb=1, norb*npol
    aorb=iorb
    if (iorb.le.npol*nocros) aorb=iorb+npol*(noins+nocros)
    do iorb1=1, norb*npol
      do iorb2=1, norb*npol
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

  call compbs_2(npol*nocros, npol*norb, n2d, 2*(n2d+npol*nocros),   &
                amat, bmat, vec, kval)

!
! To normalize (over XY plane) all the states
!
   call zgemm('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norb,   &
               one, amat, 2*n2d+npol*norb, vec, 2*n2d+npol*norb,     &
               zero, kfun, n2d)
   do ig=1,n2d
      do ik=1, 2*n2d+npol*norb
         aux(ig,ik)=amat(n2d+ig,ik)
      enddo
   enddo
   call zgemm('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norb,   &
               one, aux, n2d, vec, 2*n2d+npol*norb, zero, kfund, n2d)

  do ik=1, 2*(n2d+npol*nocros)
    raux=ddot(2*n2d,kfun(1,ik),1,kfun(1,ik),1)*sarea
    raux=1.d0/sqrt(raux)
    call dscal(2*(2*n2d+npol*norb),raux,vec(1,ik),1)
    call dscal(2*n2d,raux,kfun(1,ik),1)
    call dscal(2*n2d,raux,kfund(1,ik),1)
  enddo

!
! To find k-vector and the current of Bloch states
!

  call kbloch (2*(n2d+npol*nocros), kval)

  call jbloch(2*(n2d+npol*nocros), n2d, norbf, norb, nocros,  &
              kfun, kfund, vec, kval, intw1, intw2, nchan, npol)
!
! To save band structure result
!
  kfun=(0.d0,0.d0)
  kfund=(0.d0,0.d0)
  kint=(0.d0,0.d0)
!
! To account for the case of the right lead
!
  if (lleft.eq.0) then
    do i=1, 2*n2d
      do j=1, 2*n2d+npol*norb
        amat(i,j)=bmat(i,j)
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+npol*nocros
      do j=1, 2*n2d+npol*norb
        amat(i,j)=-bmat(i+npol*(nocros+noins),j)
      enddo
    enddo
  endif
!
! psi_k and psi'_k on the scattering region boundary
!
   call zgemm('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norb,&
               one, amat, 2*n2d+npol*norb, vec, 2*n2d+npol*norb,  &
               zero, kfun, n2d)
   do ig=1,n2d
      do ik=1, 2*n2d+npol*norb
         aux(ig,ik)=amat(n2d+ig,ik)
      enddo
   enddo
   call zgemm('n', 'n', n2d, 2*(n2d+npol*nocros), 2*n2d+npol*norb,&
               one, aux, n2d, vec, 2*n2d+npol*norb, zero, kfund, n2d)

!
! kint(iorb, ik)=\sum_{iorb1} D_{iorb,iorb1}
!            \int_{cell} W_iorb1^* psi_ik )
! for the orbitals crossing the boundary
!
  do ik=1,  2*(n2d+npol*nocros)
    do iorb=1, nocros*npol
      do j=1, 2*n2d+npol*norb
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
      if (lleft.eq.0) then
        kcoef(iorb,ik)=vec(2*n2d+iorb,ik)
      else
        kcoef(iorb,ik)=vec(2*n2d+npol*(nocros+noins)+iorb,ik)
      endif
    enddo
  enddo

!
! to set up B.S. for the right lead in the case of identical tips
!
  if(lleft.eq.1.and.ikind.eq.1) then
    nchanr=nchan
    call dcopy(2*(n2d+npol*nocros), kval, 1, kvalr, 1)
    kfunr=(0.d0,0.d0)
    kfundr=(0.d0,0.d0)
    kintr=(0.d0,0.d0)

    do i=1, 2*n2d
      do j=1, 2*n2d+npol*norb
        amat(i,j)=bmat(i,j)
      enddo
    enddo
    do i=2*n2d+1, 2*n2d+npol*nocros
      do j=1, 2*n2d+npol*norb
        amat(i,j)=-bmat(i+npol*(nocros+noins),j)
      enddo
    enddo

    do ik=1, 2*(n2d+npol*nocros)
      do ig=1, n2d
        do j=1,  2*n2d+npol*norb
          kfunr(ig,ik)= kfunr(ig,ik)+amat(ig,j)*vec(j,ik)
          kfundr(ig,ik)= kfundr(ig,ik)+amat(n2d+ig,j)*vec(j,ik)
        enddo
      enddo
    enddo
    do ik=1,  2*(n2d+npol*nocros)
      do iorb=1, nocros*npol
        do j=1, 2*n2d+npol*norb
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

!--
!  integrals of Bloch states with boundary orbitals for left/right leads
!
 if (lorb) then
   korb = 0.d0
   do ik = 1, 2*(n2d+npol*nocros)

     do iorb = 1, npol*nocros
      iorb1 = iorb + npol*(nocros+noins)
      do ig = 1, 2*n2d
        korb(iorb,ik) = korb(iorb,ik)+       &
                      intw1(iorb1,ig)*vec(ig,ik)
      enddo
      do ig = 1, norb*npol
        korb(iorb,ik) = korb(iorb,ik)+       &
               intw2(iorb1,ig)*vec(2*n2d+ig,ik)
      enddo
     enddo

     do iorb = 1, npol*nocros
      x1 = 0.d0
      do ig = 1, 2*n2d
        x1 = x1 + intw1(iorb,ig)*vec(ig,ik)
      enddo
      do ig = 1, norb*npol
        x1 = x1 + intw2(iorb,ig)*vec(2*n2d+ig,ik)
      enddo
      korb(iorb,ik) = korb(iorb,ik) + x1* &
                   exp(kval(ik)*(0.d0,1.d0)*tpi)
     enddo
   enddo

   if (ikind.ne.2.or.lleft.ne.0) korbl(:,:) = korb(:,:)
   if (ikind.ne.2.or.lleft.ne.1) then
     do ik = 1, 2*(n2d+npol*nocros)
       x1 = exp(-kval(ik)*(0.d0,1.d0)*tpi)
       do iorb = 1, npol*nocros
        korbr(iorb,ik) = korb(iorb,ik) * x1
       enddo
     enddo
   endif

 endif
!--

!--
! Computes and writes the propagating Bloch states
!
  if (lorb.and.ikind.eq.0.and.nchan /= 0) then
    allocate( veceig(nchan, nchan) )
    deallocate( aux )
    allocate( aux(4*n2d+npol*(norb+2*nocros),nchan) )

!-- right moving states
    veceig = 0.d0
    aux = 0.d0
    do ichan = 1, nchan
      do ig = 1, 2*n2d + npol*norb
        aux(ig, ichan) = vec(ig,ichan)
      enddo
    enddo
    CALL scat_states_plot(ikk,ien,norb,nocros,nchan,aux,veceig,.true.)

!-- left moving states
    veceig = 0.d0
    aux = 0.d0
    do ichan = 1, nchan
      do ig = 1, 2*n2d + npol*norb
        aux(ig, ichan) = vec(ig,n2d+npol*nocros+ichan)
      enddo
    enddo
    CALL scat_states_plot(ikk,ien,norb,nocros,nchan,aux,veceig,.false.)

    deallocate( veceig )
  endif
!--

  deallocate(amat)
  deallocate(bmat)
  deallocate(vec)
  deallocate(aux)
  IF (norb>0) THEN
     if (noncolin) then
        deallocate(zpseu_nc)
        deallocate(zps_nc)
     else
        deallocate(zpseu)
        deallocate(zps)
     endif
  ENDIF
  if (lorb) deallocate(korb)
  call stop_clock('compbs')

  return
end subroutine compbs
