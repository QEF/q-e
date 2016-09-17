
!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Optimized Aug. 2004 (ADC)
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
! Optimized Oct. 2006 (A. Smogunov)
!
!
subroutine scatter_forw(nrz, nrzp, z, psiper, zk, norb, tblm, cros, &
                        taunew, r, rab, betar)
!
! This subroutine computes local Phi_n and partial nonlocal Phi_alpha
! solutions of the Schrodinger equation in the region zin<z<zfin
! on which the general solution may be constructed:
!        \sum a_n Phi_n + \sum a_alpha Phi_alpha
! It computes also the integrals (intw1, intw2) of  Phi_n and
! Phi_alpha over beta-functions inside the unit cell.
!
  USE constants, ONLY : tpi
  USE parameters, only : npsx
  use radial_grids, only: ndmx
  USE cell_base, ONLY : tpiba
  USE noncollin_module, ONLY : npol
  USE mp_global, ONLY : intra_pool_comm
  USE cond
  !
  IMPLICIT NONE

  INTEGER ::      &
        nrz, nrzp, norb,&
        cros(norb,nrz), &
        tblm(4,norb),   &
        k, kz, n, lam, ig, lam1, mdim, itt, nbb, iorb, iorb1,   &
        iorba, iorb1a, is, kp, k1, nt, nb, kin, kfin
  INTEGER :: i, j, kp1, info
  INTEGER, ALLOCATABLE :: ipiv(:), inslab(:)
  real(DP) :: z(nrz+1), r(1:ndmx,npsx), rab(1:ndmx,npsx),        &
                   betar(1:ndmx,nbrx,npsx), taunew(4,norb)
  REAL(DP), PARAMETER :: eps=1.d-8
  REAL(DP) :: ddot, dz, tr, dz1
  COMPLEX(DP), PARAMETER :: cim=(0.d0,1.d0), one=(1.d0, 0.d0), &
                                 zero=(0.d0,0.d0)
  COMPLEX(DP) :: int1d, int2d, c, d, e, f, arg,&
                      zdotc, fact, factm, psiper(n2d,n2d,nrzp), &
                      zk(n2d,nrzp)
  COMPLEX(DP), ALLOCATABLE ::   &
     psigper(:,:), & ! psigper(g,lam)=newbg(g,lam1) psiper(lam1,lam)
     w0(:,:,:),    & ! w0(z,g,m) are 2D Fourier components (see four.f)
     w(:,:,:),     & ! w(z,lam,m)=psigper(g,lam)^* \exp{-igr^m_perp}
                     !                            w0(z,g,m)
     hmat(:,:),    & ! ci(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z-z(k))}
     amat(:,:),    & ! di(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z(k+1)-z)}
     xmat(:,:),    & !
     ci(:,:),      & ! ci(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z-z(k))}
     di(:,:),      & ! di(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z(k+1)-z)}
     cix(:,:,:),   & !
     dix(:,:,:),   & !
     f0(:,:), f1(:,:), f2(:,:), f_aux(:,:), &
     zkk(:), ezk(:), emzk(:), zk2(:), ezk1(:,:), emzk1(:,:)

  CALL start_clock('scatter_forw')

!
! Find first and last slab (kin, kfin)
!
  dz = z(2) - z(1)
  dz1 = dz/nz1
!
! Divide the slabs among CPU
!
  call divide(intra_pool_comm,nrz,kin,kfin)

  ALLOCATE( psigper( ngper*npol, n2d ) )
  ALLOCATE( w0( nz1, ngper, 5 ) )
  ALLOCATE( zkk( n2d ) )
  ALLOCATE( ezk( n2d ) )
  ALLOCATE( emzk( n2d ) )
  ALLOCATE( ezk1( nz1, n2d ) )
  ALLOCATE( emzk1( nz1, n2d ) )
  ALLOCATE( zk2( n2d ) )
  ALLOCATE( amat( 2*n2d, 2*n2d ) )
  ALLOCATE( xmat( 2*n2d, 2*n2d+norb*npol ) )
  ALLOCATE( ipiv(2*n2d) )

  IF (norb>0) THEN
     ALLOCATE( w( nz1, n2d, norb*npol ) )
     ALLOCATE( cix( nz1, n2d, norb*npol ) )
     ALLOCATE( dix( nz1, n2d, norb*npol ) )
     ALLOCATE( ci( norb*npol, n2d ) )
     ALLOCATE( di( norb*npol, n2d ) )
     ALLOCATE( inslab( norb ) )
     ALLOCATE( f0(n2d,norb*npol) )
     ALLOCATE( f1(n2d,norb*npol) )
     ALLOCATE( f2(n2d,norb*npol) )
     ALLOCATE( f_aux(norb*npol,n2d) )
     intw1=(0.d0,0.d0)
     intw2=(0.d0,0.d0)
  END IF

!
! some orbitals relations
!
  do iorb=1, norb
    inslab(iorb) = 0
  enddo
  do iorb = 1, norb
    if(inslab(iorb).eq.0.and.tblm(4,iorb).eq.1) then
      itt = tblm(1,iorb)
      nbb = tblm(2,iorb)
      do iorb1 = iorb, norb
       if(tblm(4,iorb1).eq.1.and.tblm(2,iorb1).eq.nbb) then
         tr = abs(taunew(3,iorb)-taunew(3,iorb1))
         if(tblm(1,iorb1).eq.itt.and.tr.le.eps) then
           inslab(iorb1) = iorb
         endif
       endif
      enddo
    endif
  enddo

!--- initial conditions for a_n coefficients
  xmat = 0.d0
  do lam = n2d+1, 2*n2d
    xmat(lam,lam) = 1.d0
  enddo
!---

  do k = kin, kfin
    kp = k-kin+1

!------
!  Start of 2D Fourier components calculations and depending
!  variables
!
    do lam=1,n2d
      arg=cim*tpi*zk(lam, kp)*dz
      zkk(lam)=cim*zk(lam,kp)*tpiba
      ezk(lam)=exp(arg)
      emzk(lam)=exp(-arg)
      zk2(lam)=cim/(2.d0*zk(lam,kp)*tpiba)

      arg=cim*tpi*zk(lam,kp)*dz1
      fact=exp(arg)
      factm=exp(-arg)
      ezk1(1,lam)=fact
      emzk1(1,lam)=factm
      do k1=2,nz1
        ezk1(k1,lam)=ezk1(k1-1,lam)*fact
        emzk1(k1,lam)=emzk1(k1-1,lam)*factm
      enddo
    enddo

    if(ewind.le.100.d0) then
        CALL zgemm('n', 'n', ngper*npol, n2d, n2d, one, newbg,   &
                   ngper*npol, psiper(1,1,kp), n2d, zero,        &
                   psigper, ngper*npol)
    else
        psigper(:,:) = psiper(:,:,kp)
    endif

    IF (norb>0) THEN
       w = 0.d0
       ci = 0.d0
       di = 0.d0
    ENDIF
    DO iorb=1, norb
      IF(cros(iorb,k).EQ.1.AND.inslab(iorb).EQ.iorb) THEN
       mdim = 2*tblm(3,iorb)+1
       nt = tblm(1,iorb)
       nb = tblm(2,iorb)
       call four(w0, z(k), dz, tblm(1,iorb), taunew(1,iorb),  &
                 r(1,nt), rab(1,nt), betar(1,nb,nt))
       DO iorb1=1, norb
        IF(inslab(iorb1).EQ.iorb) THEN
         DO ig=1, ngper
          tr=-tpi*ddot(2,gper(1,ig),1,taunew(1,iorb1),1)
          c=CMPLX(COS(tr),SIN(tr),kind=DP)
          DO lam=1, n2d
           DO is=1,npol
            d=CONJG(psigper(ngper*(is-1)+ig,lam))*c
            DO n=1, mdim
             DO kz=1, nz1
              w(kz,lam,npol*(iorb1-2+n)+is)= &
                 w(kz,lam,npol*(iorb1-2+n)+is)+d*w0(kz,ig,n)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDIF
    ENDDO
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF (cros(iorba,k).EQ.1) THEN
          DO lam=1, n2d
             CALL setint(w(1,lam,iorb),cix(1,lam,iorb),dix(1,lam,iorb),  &
                                ezk1(1,lam), emzk1(1,lam), nz1)
          ENDDO
       ENDIF
    ENDDO
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF (cros(iorba,k).EQ.1) THEN
          DO lam=1, n2d
             ci(iorb,lam)=int1d(w(1,lam,iorb),   &
                          zk(lam,kp),dz,dz1,nz1,tpiba,1)
             di(iorb,lam)=int1d(w(1,lam,iorb),   &
                          zk(lam,kp),dz,dz1,nz1,tpiba,-1)
          ENDDO
          DO iorb1=1, norb*npol
             iorb1a=iorb1
             IF (npol.EQ.2) iorb1a=(iorb1+1)/2
             IF (cros(iorb1a,k).EQ.1) THEN
                DO lam=1, n2d
                   intw2(iorb,iorb1)=intw2(iorb,iorb1)-                     &
                      int2d(w(1,lam,iorb),w(1,lam,iorb1),cix(1,lam,iorb1),  &
                            dix(1,lam,iorb1),ezk1(1,lam),emzk1(1,lam),      &
                            zk(lam,kp),dz1,tpiba,nz1)*zk2(lam)
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
!----

!-----------
! computes f1 and f2
!
    IF (norb>0) THEN
       f1 = 0.d0
       f2 = 0.d0
    ENDIF
    DO iorb=1, norb*npol
      iorba=iorb
      IF (npol.EQ.2) iorba=(iorb+1)/2
      IF (cros(iorba,k).EQ.1) THEN
        DO lam=1, n2d
          IF (ABS(AIMAG(zk(lam, kp))).LT.eps) THEN
           f1(lam,iorb)=-ezk(lam)*CONJG(di(iorb,lam))*zk2(lam)
           f2(lam,iorb)=-ezk(lam)*CONJG(ci(iorb,lam))*zk2(lam)
          ELSE
           f1(lam,iorb)=-CONJG(ci(iorb,lam))*zk2(lam)
           f2(lam,iorb)=-CONJG(di(iorb,lam))*zk2(lam)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
!------------

    if(kp.eq.1) then
!-------
!     b coeff. and fp on the left boundary
      fun0 = 0.d0
      funl0 = 0.d0
      IF (norb>0) f0 = f1
      do n = 1, n2d
        fun0(n,n) = 1.d0
      enddo
!-----------
      goto 11
    endif

!-------
!    adding nonlocal part
    IF (norb>0) THEN
       CALL zgemm('n','n',n2d,norb*npol,n2d,-one,psiper(1,1,kp), &
                  n2d,f1,n2d,one,funl1,n2d)
       do i = 1, norb*npol
          do lam = 1, n2d
             f1(lam,i) = -zkk(lam)*f1(lam,i)
          enddo
       enddo
       CALL zgemm('n','n',n2d,norb*npol,n2d,-one,psiper(1,1,kp), &
                  n2d,f1,n2d,one,fundl1,n2d)
    END IF
!-------

!------
!    constructs matrices
    do i = 1, n2d
      do j = 1, n2d
        amat(i,j) = fun1(i,j)
        amat(n2d+i,j) = fund1(i,j)
        amat(i,n2d+j) = -psiper(i,j,kp)
        amat(n2d+i,n2d+j) = amat(i,n2d+j)*zkk(j)
        xmat(i,j) = psiper(i,j,kp)*ezk(j)
        xmat(n2d+i,j) = amat(n2d+i,n2d+j)*ezk(j)
      enddo
      do j = n2d+1, 2*n2d
        xmat(i,j) = -fun1(i,j)
        xmat(n2d+i,j) = -fund1(i,j)
      enddo
      do j = 1, norb*npol
        xmat(i,2*n2d+j) = -funl1(i,j)
        xmat(n2d+i,2*n2d+j) = -fundl1(i,j)
      enddo
    enddo
!------

!    Solve the system of linear equations
    call ZGESV(2*n2d,2*n2d+norb*npol,amat,2*n2d,ipiv,xmat,2*n2d,info)

!-------
!  rotates integrals
!
    IF (norb>0) THEN
       do i = 1, norb*npol
          do j = 1, n2d
             f_aux(i,j) = intw1(i,j)
          enddo
       enddo
       call zgemm('n','n',norb*npol,n2d,n2d,one,f_aux,norb*npol,      &
                  xmat(1,n2d+1),2*n2d,one,intw1(1,n2d+1),norbf*npol)
       call zgemm('n','n',norb*npol,norb*npol,n2d,one,f_aux,norb*npol,&
                  xmat(1,2*n2d+1),2*n2d,one,intw2,norbf*npol)
       call zgemm('n','n',norb*npol,n2d,n2d,one,f_aux,norb*npol,xmat, &
                  2*n2d,zero,intw1,norbf*npol)
    ENDIF
!--------

!-------
! rotates b coeff. on the left boundary
!
    do i = 1, n2d
      do j = 1, n2d
       amat(i,j) = fun0(i,j)
      enddo
    enddo
    call zgemm('n','n',n2d,n2d,n2d,one,amat,2*n2d,xmat,          &
               2*n2d,zero,fun0,n2d)
    call zgemm('n','n',n2d,n2d,n2d,one,amat,2*n2d,xmat(1,n2d+1), &
               2*n2d,one,fun0(1,n2d+1),n2d)
    IF (norb>0) &
       call zgemm('n','n',n2d,norb*npol,n2d,one,amat,2*n2d,         &
                  xmat(1,2*n2d+1),2*n2d,one,funl0,n2d)
!---------------

!-------
! rotates all previous functions if lorb is .true.
!
    IF (lorb) THEN
       DO kp1=kp,2,-1
          DO i = 1, n2d
             DO j = 1, n2d
                amat(i,j) = funz0(i,j,kp1)
             END DO
          END DO
          CALL zgemm('n','n',n2d,n2d,n2d,one,amat,2*n2d,xmat,          &
                      2*n2d,zero,funz0(1,1,kp1),n2d)
          CALL zgemm('n','n',n2d,n2d+norb*npol,n2d,one,amat,2*n2d,     &
                      xmat(1,n2d+1), 2*n2d,one,funz0(1,n2d+1,kp1),n2d)
       END DO
    END IF

11  continue

!------
!   Add to the integrals
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF(cros(iorba,k).EQ.1) THEN
          DO n=1, n2d
             DO lam=1, n2d
              intw1(iorb,n) = intw1(iorb,n) + xmat(n2d+lam,n)*ci(iorb,lam)
              intw1(iorb,n2d+n) = intw1(iorb,n2d+n) + &
                                  xmat(n2d+lam,n2d+n)*ci(iorb,lam)
             ENDDO
             intw1(iorb,n)=intw1(iorb,n)+di(iorb,n)
          ENDDO
          DO iorb1=1, norb*npol
             iorb1a=iorb1
             IF (npol.EQ.2) iorb1a=(iorb1+1)/2
             tr=taunew(3,iorb1a)-taunew(4,iorb1a)
             IF (z(k)+dz.GT.tr) THEN
               c=(0.d0, 0.d0)
               DO lam=1, n2d
                  c=c+xmat(n2d+lam,2*n2d+iorb1)*ci(iorb,lam)
               ENDDO
               intw2(iorb,iorb1)=intw2(iorb,iorb1)+c
             ENDIF
          ENDDO
       ENDIF
    ENDDO
!---------------


!-------
!    wave functions on the right boundary
    do i = 1, n2d
     do j = 1, 2*n2d
       amat(i,j) = xmat(n2d+i,j)*ezk(i)
       amat(n2d+i,j) = amat(i,j)*zkk(i)
     enddo
     amat(i,i) = amat(i,i) + 1.d0
     amat(n2d+i,i) = amat(n2d+i,i) - zkk(i)
     do j = 1, norb*npol
       f1(i,j)=xmat(n2d+i,2*n2d+j)*ezk(i)+f2(i,j)
       f2(i,j)=f1(i,j)*zkk(i)
     enddo
    enddo
    CALL zgemm('n','n',n2d,2*n2d,n2d,one,psiper(1,1,kp),     &
               n2d,amat,2*n2d,zero,fun1,n2d)
    CALL zgemm('n','n',n2d,2*n2d,n2d,one,psiper(1,1,kp),     &
               n2d,amat(n2d+1,1),2*n2d,zero,fund1,n2d)
    IF (norb>0) THEN
       CALL zgemm('n','n',n2d,norb*npol,n2d,one,psiper(1,1,kp), &
                  n2d,f1,n2d,zero,funl1,n2d)
       CALL zgemm('n','n',n2d,norb*npol,n2d,one,psiper(1,1,kp), &
                  n2d,f2,n2d,zero,fundl1,n2d)
    END IF

    IF (lorb.and.kp<nrzp) THEN
       DO i=1,n2d
          DO j=1,2*n2d
             funz0(i,j,kp+1)=fun1(i,j)
          END DO
          DO j=1,norb*npol
             funz0(i,2*n2d+j,kp+1)=funl1(i,j)
          END DO
       END DO
    END IF
!---------

  enddo

!-------
!    wave functions on the left boundary
  do lam = 1, n2d
    arg = cim*tpi*zk(lam,1)*dz
    zkk(lam) = cim*zk(lam,1)*tpiba
    ezk(lam) = exp(arg)
  enddo

  do i = 1, n2d
    do j = 1, 2*n2d
      amat(i,j) = fun0(i,j)*ezk(i)
      amat(n2d+i,j) = -amat(i,j)*zkk(i)
    enddo
    amat(i,n2d+i) = amat(i,n2d+i) + 1.d0
    amat(n2d+i,n2d+i) = amat(n2d+i,n2d+i) + zkk(i)
    do j = 1, norb*npol
      f1(i,j) = funl0(i,j)*ezk(i)+f0(i,j)
      f2(i,j) = -f1(i,j)*zkk(i)
    enddo
  enddo
  CALL zgemm('n','n',n2d,2*n2d,n2d,one,psiper(1,1,1),     &
             n2d,amat,2*n2d,zero,fun0,n2d)
  CALL zgemm('n','n',n2d,2*n2d,n2d,one,psiper(1,1,1),     &
             n2d,amat(n2d+1,1),2*n2d,zero,fund0,n2d)
  IF (norb > 0) THEN
     CALL zgemm('n','n',n2d,norb*npol,n2d,one,psiper(1,1,1), &
                n2d,f1,n2d,zero,funl0,n2d)
     CALL zgemm('n','n',n2d,norb*npol,n2d,one,psiper(1,1,1), &
                n2d,f2,n2d,zero,fundl0,n2d)
  ENDIF
!---------
  IF (lorb) THEN
     DO i=1,n2d
        DO j=1,2*n2d
           funz0(i,j,1)=fun0(i,j)
        END DO
        DO j=1,norb*npol
           funz0(i,2*n2d+j,1)=funl0(i,j)
        END DO
     END DO
  END IF
!---------

! scaling the integrals
  IF (norbf > 0) THEN
     CALL dscal(2*norbf*npol*2*n2d, sarea, intw1, 1)
     CALL dscal(2*norbf*npol*norbf*npol, sarea, intw2, 1)
  ENDIF

!
! To construct the functions in the whole rigion zin<z<zfin in the
! case of multiparallel running
!
#if defined(__MPI)
  CALL rotproc(fun0, fund0, fun1, fund1, funl0, fundl0, funl1, &
               fundl1, intw1, intw2, n2d, norbf, norb, nrzp)
#endif

  deallocate( psigper )
  deallocate( w0 )
  deallocate( zkk )
  deallocate( ezk )
  deallocate( emzk )
  deallocate( ezk1 )
  deallocate( emzk1 )
  deallocate( zk2 )
  deallocate( amat )
  deallocate( xmat )
  deallocate( ipiv )

  IF (norb>0) THEN
     deallocate( w )
     deallocate( cix )
     deallocate( dix )
     deallocate( ci )
     deallocate( di )
     deallocate( inslab )
     deallocate( f0 )
     deallocate( f1 )
     deallocate( f2 )
     deallocate( f_aux )
  END IF
  CALL stop_clock('scatter_forw')

  return
END SUBROUTINE scatter_forw

