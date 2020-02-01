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
SUBROUTINE local (ien)
!
! This subroutine computes 2D eigenfunctions and eigenvalues for
! the local potential in each slab and performs 2D reduction of
! the plane wave basis set (local_1). Using this reduced basis
! set it solves again 2D EV problem (local_2).
!
  USE constants, ONLY : rytoev
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : npol
  USE io_files
  USE cond
  !
  USE mp_pools, ONLY : intra_pool_comm
  !

  IMPLICIT NONE

  INTEGER :: ien, ig, il, k, kin, kfin
  REAL(DP) :: edummy
  COMPLEX(DP), ALLOCATABLE :: psibase(:,:)
  LOGICAL :: exst

!
! To divide the slabs between CPU
!
  CALL start_clock('local')

!
! If all the information is already contained in the file it reads it.
!
  IF (lread_loc) THEN
    CALL seqopn(4,fil_loc,'unformatted',exst)
    IF(.NOT.exst) CALL errore ('local','fil_loc not found',1)
    READ(4) n2d
    READ(4) nrzpl, nrzps, nrzpr
!   Allocate variables depending on n2d
    CALL allocate_cond

    READ(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d)
    READ(4) (((psiperl(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzpl)
    READ(4) ((zkrl(il,k),il=1,n2d),k=1,nrzpl)
    if(ikind.gt.0) then
     READ(4) (((psipers(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzps)
     READ(4) ((zkrs(il,k),il=1,n2d),k=1,nrzps)
    endif
    if(ikind.gt.1) then
     READ(4) (((psiperr(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzpr)
     READ(4) ((zkrr(il,k),il=1,n2d),k=1,nrzpr)
    endif
    CLOSE(unit=4)
    RETURN
  ENDIF

  allocate( psibase( ngper*npol, ngper*npol ) )
  psibase = 0.d0

  if(ewind.le.100.d0) then
    n2d = 0
    edummy = earr(ien)/rytoev + efl
    call local_1(edummy,nrzl,vppotl,n2d,psibase)
    if(ikind.gt.0) then
     edummy = earr(ien)/rytoev + efs
     call local_1(edummy,nrzs,vppots,n2d,psibase)
    endif
    if(ikind.eq.2) then
     edummy = earr(ien)/rytoev + efr
     call local_1(edummy,nrzr,vppotr,n2d,psibase)
    endif
  else
    n2d = ngper*npol
  endif

!
! Allocate variables depending on n2d
!
  nrzps = 0
  nrzpr = 0
  call divide(intra_pool_comm, nrzl, kin, kfin)
  nrzpl = kfin - kin + 1
  if(ikind.gt.0) then
    call divide(intra_pool_comm, nrzs, kin, kfin)
    nrzps = kfin - kin + 1
  endif
  if(ikind.gt.1) then
    call divide(intra_pool_comm, nrzr, kin, kfin)
    nrzpr = kfin - kin + 1
  endif

  CALL allocate_cond

  IF (npol.EQ.2) THEN
     WRITE( stdout,*) 'ngper, ngper*npol, n2d = ', ngper, ngper*npol, n2d
  ELSE
     WRITE( stdout,*) 'ngper, n2d = ', ngper, n2d
  ENDIF

!
! Construct components of basis vector set on G_per
!
  if(ewind.le.100.d0) then
    CALL dcopy(2*ngper*npol*n2d,psibase,1,newbg,1)
  else
    newbg = 0.d0
    do ig=1, n2d
      newbg(ig,ig) = 1.d0
    enddo
  endif

  deallocate( psibase )

  call local_2(nrzl,nrzpl,vppotl,psiperl,zkrl)
  if(ikind.gt.0) call local_2(nrzs,nrzps,vppots,psipers,zkrs)
  if(ikind.gt.1) call local_2(nrzr,nrzpr,vppotr,psiperr,zkrr)

!
! saving the 2D data on the file if lwrite_loc=.t.
!
  IF (lwrite_loc) THEN
    IF(fil_loc.EQ.' ') CALL errore ('local','fil_loc no name',1)
    CALL seqopn(4,fil_loc,'unformatted',exst)
    WRITE(4) n2d
    WRITE(4) nrzpl, nrzps, nrzpr
    WRITE(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d)
    WRITE(4) (((psiperl(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzpl)
    WRITE(4) ((zkrl(il,k),il=1,n2d),k=1,nrzpl)
    if(ikind.gt.0) then
     WRITE(4) (((psipers(ig,il,k),ig=1,n2d),il=1,n2d), &
                                             k=1,nrzps)
     WRITE(4) ((zkrs(il,k),il=1,n2d),k=1,nrzps)
    endif
    if(ikind.gt.1) then
     WRITE(4) (((psiperr(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzpr)
     WRITE(4) ((zkrr(il,k),il=1,n2d),k=1,nrzpr)
    endif
    CLOSE(unit=4)
  ENDIF

  CALL stop_clock('local')
  RETURN
END SUBROUTINE local
!-----------------------------------

subroutine local_1 (edummy, nrz, vppot, n2d, psibase)

  USE kinds, only : DP
  USE cell_base, ONLY : at, tpiba2
  USE noncollin_module, ONLY : npol
  USE mp_world,        ONLY : world_comm, nproc
  USE mp_pools,        ONLY : me_pool, root_pool, intra_pool_comm
  USE mp,         ONLY : mp_barrier, mp_bcast
  USE io_global, ONLY : ionode, ionode_id
  USE parallel_include
  use cond, only : nrx, nry, ngper, gper, ewind, epsproj
  !
  !
  IMPLICIT NONE

  INTEGER :: nrz, n2d
  INTEGER :: i, il, j, jl, ig, jg, ipol,      &
             idx, number, nprob, nteam, nteamnow,     &
             info, kin, kfin, is, js
#if defined(__MPI)
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER, ALLOCATABLE :: fftxy(:,:)
  REAL(DP) :: edummy
  REAL(DP), ALLOCATABLE :: el(:), gp(:)
  complex(DP) :: psibase(ngper*npol,ngper*npol),   &
                      vppot(nrz,nrx*nry,npol,npol)
  COMPLEX(DP), ALLOCATABLE :: amat(:,:), psiprob(:,:)
  COMPLEX(DP),PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)


  ALLOCATE( gp( 2 ) )
  ALLOCATE( el( ngper * npol ) )
  ALLOCATE( amat( ngper * npol, ngper * npol ) )
  ALLOCATE( psiprob( ngper * npol, ngper * npol ) )
  ALLOCATE( fftxy(-nrx:nrx,-nry:nry) )

!
! To form fftxy correspondence
!
  fftxy = 0
  DO i=1, nrx
    il=i-1
    IF (il.GT.nrx/2) il=il-nrx
    DO j=1, nry
       jl=j-1
       IF (jl.GT.nry/2) jl=jl-nry
       fftxy(il,jl)=i+(j-1)*nrx
    ENDDO
  ENDDO

!
! Starting k and number of CPU
!
  kin = 1
  kfin = nrz
  kin = kin + me_pool
  nteam = nproc
  nprob=0
!
! set and solve the eigenvalue equation for each slab
!
  DO WHILE(kin.LE.kfin)

     amat=(0.d0,0.d0)
     DO ig=1, ngper
        DO jg=ig, ngper
           DO ipol=1, 2
              gp(ipol)=gper(ipol,ig)-gper(ipol,jg)
           ENDDO
           idx=number(gp, at, fftxy, nrx, nry)
           IF (idx.GT.0) THEN
              DO is=1,npol
                 DO js=1,npol
                    amat(ig+(is-1)*ngper,jg+(js-1)*ngper)=vppot(kin,idx,is,js)
                    amat(jg+(js-1)*ngper,ig+(is-1)*ngper)= &
                        CONJG(amat(ig+(is-1)*ngper,jg+(js-1)*ngper))
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        DO is=1,npol
           amat(ig+(is-1)*ngper,ig+(is-1)*ngper)=    &
             amat(ig+(is-1)*ngper,ig+(is-1)*ngper)+ &
                   (gper(1,ig)**2 + gper(2,ig)**2)*tpiba2
        ENDDO
     ENDDO
     CALL hev_ab(ngper*npol, amat, ngper*npol, el, psiprob,      &
                     -1.d1, edummy+ewind, nprob)

#if defined(__MPI)
    IF ( me_pool.ne.root_pool ) THEN
        CALL mpi_send(nprob,1,MPI_INTEGER,0,17,     &
                                    MPI_COMM_WORLD,info )
      CALL errore ('n2d reduction','info<>0 in send',info)
      CALL mpi_send(psiprob,2*ngper*npol*ngper*npol,MPI_DOUBLE_PRECISION,0,18,&
                                    MPI_COMM_WORLD,info )
      CALL errore ('n2d reduction','info<>0 in send',info)
    ELSE
      CALL gramsh(ngper*npol,nprob,1,nprob,           &
                 psibase,psiprob,n2d,epsproj)

      nteamnow=kfin-kin+1
      IF(nteamnow.GT.nteam) nteamnow=nteam

      DO ig=1, nteamnow-1
        CALL mpi_recv(nprob,1,MPI_INTEGER,       &
                       ig,17,MPI_COMM_WORLD,status,info )
        CALL errore ('n2d reduction','info<>0 in recv',info)
        CALL mpi_recv(psiprob,2*ngper*npol*ngper*npol,MPI_DOUBLE_PRECISION,  &
                       ig,18,MPI_COMM_WORLD,status,info )
        CALL errore ('n2d reduction','info<>0 in recv',info)
        CALL gramsh(ngper*npol,nprob,1,nprob,         &
                 psibase,psiprob,n2d,epsproj)
      ENDDO
    ENDIF
#else
    CALL gramsh(ngper*npol,nprob,1,nprob,psibase,psiprob,n2d,epsproj)
#endif
    kin=kin+nteam

  ENDDO

#if defined(__MPI)
  CALL mp_barrier(world_comm)
  CALL mp_bcast(n2d,ionode_id, world_comm)
  CALL mp_bcast(psibase,ionode_id, world_comm)
#endif

  deallocate( gp )
  deallocate( el )
  deallocate( amat )
  deallocate( psiprob )
  deallocate( fftxy )

  return
end subroutine local_1

subroutine local_2(nrz, nrzp, vppot, psiper, zkr)

  USE kinds, only : DP
  USE cell_base, ONLY : at, tpiba2
  USE mp,         ONLY : mp_barrier
  USE noncollin_module, ONLY : npol
  use cond, only : nrx, nry, ngper, n2d, gper, newbg
  !
  USE mp_pools, ONLY : intra_pool_comm
  USE mp_world, ONLY : world_comm
  !

  IMPLICIT NONE

  INTEGER :: nrz, nrzp
  INTEGER :: i, il, j, jl, ig, jg, ipol, k, kp, &
             info, idx, number, kin, kfin, is, js
  INTEGER, ALLOCATABLE :: fftxy(:,:)
  REAL(DP) :: zkr(n2d,nrzp)
  REAL(DP), ALLOCATABLE :: gp(:)
  complex(DP) :: psiper(n2d,n2d,nrzp),   &
                      vppot(nrz,nrx*nry,npol,npol), aij, zdotc
  COMPLEX(DP), ALLOCATABLE :: amat(:,:), amat1(:,:), ymat(:,:)
  COMPLEX(DP),PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)

  allocate( gp( 2 ) )
  allocate( fftxy(-nrx:nrx, -nry:nry) )
  ALLOCATE( amat( ngper * npol, ngper * npol ) )
  ALLOCATE( amat1( n2d, n2d ) )
  ALLOCATE( ymat( ngper*npol, n2d ) )

!
! To form fftxy correspondence
!
  fftxy(:,:) = 0
  do i = 1, nrx
    il = i-1
    if (il.gt.nrx/2) il=il-nrx
    do j = 1, nry
       jl = j-1
       if (jl.gt.nry/2) jl = jl-nry
       fftxy(il,jl) = i+(j-1)*nrx
    enddo
  enddo

  call divide(intra_pool_comm, nrz, kin, kfin)

! for reduced basis set H'_{ab}=e*^i_aH_{ij}e^j_b
  do k = kin, kfin
      kp = k - kin + 1
      ymat=(0.d0,0.d0)
!
!     First compute y_{ib}=H_{ij}e_{jb}
!
      DO ig=1, ngper
         DO jg=1, ngper
            DO ipol=1, 2
               gp(ipol) = gper(ipol,ig) - gper(ipol,jg)
            ENDDO
            idx=number(gp, at, fftxy, nrx, nry)
            DO is=1,npol
               DO js=1,npol
                  IF (idx.GT.0) THEN
                     aij=vppot(k,idx,is,js)
                  ELSE
                     aij=(0.d0,0.d0)
                  ENDIF
                  IF ((ig.EQ.jg).AND.(is.EQ.js))          &
                     aij=aij+(gper(1,ig)**2+              &
                              gper(2,ig)**2)*tpiba2
                     amat(ig+(is-1)*ngper,jg+(js-1)*ngper)= aij
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL zgemm('n','n',ngper*npol,n2d,ngper*npol,one,amat,ngper*npol, &
                         newbg,ngper*npol,zero,ymat,ngper*npol)
!
!     and construct H'_{ab}=<e_a|y_b>
!
      DO il=1, n2d
        DO jl=il, n2d
          amat1(il,jl)=zdotc(ngper*npol,newbg(1,il),1,ymat(1,jl),1)
          amat1(jl,il)=CONJG(amat1(il,jl))
        ENDDO
      ENDDO
!
!     Solving the eigenvalue problem and construction zk
!
      info=-1
      CALL hev_ab(n2d, amat1, n2d, zkr(1,kp),       &
                  psiper(1,1,kp), 0.d0, 0.d0, info)


  ENDDO

#if defined(__MPI)
  CALL mp_barrier(world_comm)
#endif

  deallocate(amat)
  deallocate(amat1)
  deallocate(ymat)
  deallocate(gp)
  deallocate(fftxy)

  return
end subroutine local_2

FUNCTION number(gp, at, fftxy, nrx, nry)
!
! This function receives as input the coordinates of 2D g vector
! and write on output its fft position.
!
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER :: nrx, nry, fftxy(-nrx:nrx, -nry:nry), &
             number, n1, n2
  REAL(DP) :: gp(2), at(3,3), x1, x2

  x1=gp(1)*at(1,1)+gp(2)*at(2,1)
  x2=gp(1)*at(1,2)+gp(2)*at(2,2)
  n1=NINT(x1)
  n2=NINT(x2)
  IF (n1.LE.nrx/2.AND.n1.GE.-(nrx-1)/2.AND.    &
      n2.LE.nry/2.AND.n2.GE.-(nry-1)/2) THEN
    number=fftxy(n1,n2)
  ELSE
!
! The g vector is not inside the 2D mesh
!
    number=-1
  ENDIF

  RETURN
END FUNCTION number

