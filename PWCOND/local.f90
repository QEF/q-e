!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
#include "f_defs.h"
!
SUBROUTINE local
!
! This subroutine computes 2D eigenfunctions and eigenvalues for
! the local potential in each slab and performs 2D reduction of 
! the plane wave basis set.
!

  USE io_global,        ONLY : stdout, ionode
  USE pwcom
  USE noncollin_module, ONLY : npol
  USE io_files
  USE cond
  USE mp_global,        ONLY : nproc, me_pool, root_pool
  USE parallel_include

  IMPLICIT NONE 

  INTEGER :: i, il, j, jl, ixy, ig, jg, ipol, igper, k,      &
             ios, index, number, nprob, nteam, nteamnow,     &
             status, info, kin, kfin, is, js
  INTEGER, ALLOCATABLE :: fftxy(:,:) 
  REAL(kind=DP), PARAMETER :: eps=1.d-6
  REAL(kind=DP), ALLOCATABLE :: el(:), gp(:)
  COMPLEX(kind=DP), ALLOCATABLE :: amat(:,:), amat1(:,:), ymat(:,:),     &
                                   psibase(:,:), psiprob(:,:)
  COMPLEX(kind=DP),PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  COMPLEX(kind=DP) :: aij, xfact, ZDOTC
  LOGICAL :: exst

!
! To divide the slabs between CPU
!
  CALL start_clock('local')
  CALL slabcpu(nrz, nrzp, nkofz, bdl1, bdl2, bdr1, bdr2, z)

!
! If all the information is already contained in the file it reads it.
!
  IF (lread_loc) THEN 
    CALL seqopn(4,fil_loc,'unformatted',exst)
    IF(.NOT.exst) CALL errore ('local','fil_loc not found',1)
    READ(4) n2d
!   Allocate variables depending on n2d
    CALL allocate_cond_2
    READ(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d) 
!    WRITE( stdout,*) 'ngper, n2d = ', ngper, n2d
    READ(4) (((psiper(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzp)
    READ(4) ((zkr(il,k),il=1,n2d),k=1,nrzp)

    CLOSE(unit=4)
    RETURN
  ENDIF

  ALLOCATE( gp( 2 ) )
  ALLOCATE( el( ngper * npol ) )  
  ALLOCATE( amat( ngper * npol, ngper * npol ) )
  ALLOCATE( psibase( ngper * npol, ngper * npol ) )
  ALLOCATE( psiprob( ngper * npol, ngper * npol ) )
  ALLOCATE( fftxy(-(nrx-1)/2:nrx/2,-(nry-1)/2:nry/2) )  

!
! To form fftxy correspondence
!      
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
! to find kin and kfin
!
  DO k=1, nrz
    IF (z(k).LE.bdl1+eps)  kin=k
    IF (z(k).LE.bdr2-eps)  kfin=k
  ENDDO             

!
! Starting k and number of CPU
!
  nteam=1

  kin = kin + me_pool
  nteam = nproc

!
! set and solve the eigenvalue equation for each slab
!                                                       
  n2d=0
  nprob=0
  psibase=(0.d0,0.d0)

  DO WHILE(kin.LE.kfin) 
     amat=(0.d0,0.d0)
     DO ig=1, ngper
        DO jg=ig, ngper
           DO ipol=1, 2
              gp(ipol)=gper(ipol,ig)-gper(ipol,jg)
           ENDDO
           index=number(gp, at, fftxy, nrx, nry)
           IF (index.GT.0) THEN
              DO is=1,npol
                 DO js=1,npol
                    amat(ig+(is-1)*ngper,jg+(js-1)*ngper)=vppot(kin,index,is,js)
                    amat(jg+(js-1)*ngper,ig+(is-1)*ngper)= &
                        DCONJG(amat(ig+(is-1)*ngper,jg+(js-1)*ngper))
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
                     -1.d1, eryd+ewind, nprob)

!     do is=1,ngper*npol
!        write(6,'("------------------------------",i5,f15.7)') is, el(is)
!        do ig=1,ngper*npol
!           if (ig.eq.ngper+1) write(6,'("----")')
!           write(6,'(i5,2f15.7)') ig, psiprob(ig,is)
!        enddo
!     enddo
!     stop

#ifdef __PARA
     IF ( me_pool == root_pool ) THEN
        CALL mpi_send(nprob,1,MPI_INTEGER,0,17,     &
                                    MPI_COMM_WORLD,info )
      CALL errore ('n2d reduction','info<>0 in send',info)       
      CALL mpi_send(psiprob,2*ngper*npol*ngper*npol,MPI_REAL8,0,18, &
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
        CALL mpi_recv(psiprob,2*ngper*npol*ngper*npol,MPI_REAL8,  &
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

#ifdef __PARA
  CALL mpi_barrier( )
  CALL mpi_bcast(n2d,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  CALL errore ('reduction','mpi_bcast 1',info) 
  CALL mpi_bcast(psibase,2*ngper*npol*ngper*npol,MPI_REAL8,0,  &
                  MPI_COMM_WORLD,info)
  CALL errore ('reduction','mpi_bcast 1',info)       
#endif

!
! Allocate variables depending on n2d
!
  CALL allocate_cond_2
  IF (npol.EQ.2) THEN
     WRITE( stdout,*) 'ngper, ngper*npol, n2d = ', ngper, ngper*npol, n2d
  ELSE
     WRITE( stdout,*) 'ngper, n2d = ', ngper, n2d
  ENDIF
!
! Construct components of basis vector set on G_per
!
  CALL DCOPY(2*ngper*npol*n2d,psibase,1,newbg,1)

!
! set and solve the eigenvalue equation for each slab
!
  ALLOCATE( amat1( n2d, n2d ) )
  ALLOCATE( ymat( ngper*npol, n2d ) )

! for reduced basis set H'_{ab}=e*^i_aH_{ij}e^j_b
  DO k=1, nrz
    IF(nkofz(k).NE.0) THEN
      ymat=(0.d0,0.d0)
!
!     First compute y_{ib}=H_{ij}e_{jb}
!
      DO ig=1, ngper
         DO jg=1, ngper
            DO ipol=1, 2
               gp(ipol) = gper(ipol,ig) - gper(ipol,jg)
            ENDDO
            index=number(gp, at, fftxy, nrx, nry)
            DO is=1,npol
               DO js=1,npol
                  IF (index.GT.0) THEN
                     aij=vppot(k,index,is,js)
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
      CALL ZGEMM('n','n',ngper*npol,n2d,ngper*npol,one,amat,ngper*npol, &
                         newbg,ngper*npol,zero,ymat,ngper*npol)
!
!     and construct H'_{ab}=<e_a|y_b>
!
      DO il=1, n2d
        DO jl=il, n2d
          amat1(il,jl)=ZDOTC(ngper*npol,newbg(1,il),1,ymat(1,jl),1)
          amat1(jl,il)=CONJG(amat1(il,jl))
        ENDDO
      ENDDO
!
!     Solving the eigenvalue problem and construction zk
!
      info=-1
      CALL hev_ab(n2d, amat1, n2d, zkr(1,nkofz(k)),       &
                  psiper(1,1,nkofz(k)), 0.d0, 0.d0, info)

    ENDIF
  ENDDO

#ifdef __PARA
  CALL mpi_barrier()
#endif

!
! saving the 2D data on the file if lwrite_loc=.t. 
!
  IF (lwrite_loc) THEN
    IF(fil_loc.EQ.' ') CALL errore ('local','fil_loc no name',1)
    CALL seqopn(4,fil_loc,'unformatted',exst)
    WRITE(4) n2d
    WRITE(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d)

    WRITE(4) (((psiper(ig,il,k),ig=1,n2d),il=1,n2d), &
                                                k=1,nrzp)
    WRITE(4) ((zkr(il,k),il=1,n2d),k=1,nrzp)
    CLOSE(unit=4)
  ENDIF

  DEALLOCATE(amat)
  DEALLOCATE(amat1)
  DEALLOCATE(ymat)
  DEALLOCATE(gp)
  DEALLOCATE(psibase)
  DEALLOCATE(psiprob)
  DEALLOCATE(el)
  DEALLOCATE(fftxy)

  CALL stop_clock('local')
  RETURN 
END SUBROUTINE local
!-----------------------------------

FUNCTION number(gp, at, fftxy, nrx, nry)
!
! This function receives as input the coordinates of 2D g vector 
! and write on output its fft position. 
!
  IMPLICIT NONE
  INTEGER :: nrx, nry, fftxy(-(nrx-1)/2:nrx/2, -(nry-1)/2:nry/2), &
             number, n1, n2
  REAL(kind=KIND(0.d0)) :: gp(2), at(3,3), x1, x2 
  REAL(kind=KIND(0.d0)), PARAMETER :: eps=1.d-4

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
      
