
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
subroutine local
!
! This subroutine computes 2D eigenfunctions and eigenvalues for
! the local potential in each slab and performs 2D reduction of 
! the plane wave basis set.
!
#include "f_defs.h"
  USE io_global,  ONLY :  stdout
  USE pwcom
  USE noncollin_module, ONLY : npol
  USE io_files
  USE cond
#ifdef __PARA
  USE mp_global, ONLY: nproc
  use para
  USE parallel_include
#endif      

  implicit none 

  integer :: i, il, j, jl, ixy, ig, jg, ipol, igper, k,      &
             ios, index, number, nprob, nteam, nteamnow,     &
             status, info, kin, kfin, is, js
  integer, allocatable :: fftxy(:,:) 
  real(kind=DP), parameter :: eps=1.d-6
  real(kind=DP), allocatable :: el(:), gp(:)
  complex(kind=DP), allocatable :: amat(:,:), amat1(:,:), ymat(:,:),     &
                                   psibase(:,:), psiprob(:,:)
  complex(kind=DP),parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  complex(kind=DP) :: aij, xfact, ZDOTC
  logical :: exst

!
! To divide the slabs between CPU
!
  call start_clock('local')
  call slabcpu(nrz, nrzp, nkofz, bdl1, bdl2, bdr1, bdr2, z)

!
! If all the information is already contained in the file it reads it.
!
  if (lread_loc) then 
    call seqopn(4,fil_loc,'unformatted',exst)
    if(.not.exst) call errore ('local','fil_loc not found',1)
    read(4) n2d
!   Allocate variables depending on n2d
    call allocate_cond_2
    read(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d) 
!    WRITE( stdout,*) 'ngper, n2d = ', ngper, n2d
    read(4) (((psiper(ig,il,k),ig=1,n2d),il=1,n2d), &
                                               k=1,nrzp)
    read(4) ((zkr(il,k),il=1,n2d),k=1,nrzp)

    close(unit=4)
    return
  endif

  allocate( gp( 2 ) )
  allocate( el( ngper * npol ) )  
  allocate( amat( ngper * npol, ngper * npol ) )
  allocate( psibase( ngper * npol, ngper * npol ) )
  allocate( psiprob( ngper * npol, ngper * npol ) )
  allocate( fftxy(-(nrx-1)/2:nrx/2,-(nry-1)/2:nry/2) )  

!
! To form fftxy correspondence
!      
  do i=1, nrx
    il=i-1
    if (il.gt.nrx/2) il=il-nrx 
    do j=1, nry
       jl=j-1
       if (jl.gt.nry/2) jl=jl-nry  
       fftxy(il,jl)=i+(j-1)*nrx   
    enddo
  enddo     

!
! to find kin and kfin
!
  do k=1, nrz
    if (z(k).le.bdl1+eps)  kin=k
    if (z(k).le.bdr2-eps)  kfin=k
  enddo             

!
! Starting k and number of CPU
!
  nteam=1
#ifdef __PARA
  kin=kin+me-1
  nteam=nproc
#endif

!
! set and solve the eigenvalue equation for each slab
!                                                       
  n2d=0
  nprob=0
  psibase=(0.d0,0.d0)

  do while(kin.le.kfin) 
     amat=(0.d0,0.d0)
     do ig=1, ngper
        do jg=ig, ngper
           do ipol=1, 2
              gp(ipol)=gper(ipol,ig)-gper(ipol,jg)
           enddo
           index=number(gp, at, fftxy, nrx, nry)
           if (index.gt.0) then
              do is=1,npol
                 do js=1,npol
                    amat(ig+(is-1)*ngper,jg+(js-1)*ngper)=vppot(kin,index,is,js)
                    amat(jg+(js-1)*ngper,ig+(is-1)*ngper)= &
                        DCONJG(amat(ig+(is-1)*ngper,jg+(js-1)*ngper))
                 enddo
              enddo
           endif
        enddo
        do is=1,npol
           amat(ig+(is-1)*ngper,ig+(is-1)*ngper)=    &
             amat(ig+(is-1)*ngper,ig+(is-1)*ngper)+ &
                   (gper(1,ig)**2 + gper(2,ig)**2)*tpiba2
        enddo
     enddo       
     call hev_ab(ngper*npol, amat, ngper*npol, el, psiprob,      &
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
     if ( me.ne.1 ) then
        call mpi_send(nprob,1,MPI_INTEGER,0,17,     &
                                    MPI_COMM_WORLD,info )
      call errore ('n2d reduction','info<>0 in send',info)       
      call mpi_send(psiprob,2*ngper*npol*ngper*npol,MPI_REAL8,0,18, &
                                    MPI_COMM_WORLD,info )           
      call errore ('n2d reduction','info<>0 in send',info)
    else
      call gramsh(ngper*npol,nprob,1,nprob,           &
                 psibase,psiprob,n2d,epsproj) 

      nteamnow=kfin-kin+1
      if(nteamnow.gt.nteam) nteamnow=nteam

      do ig=1, nteamnow-1 
        call mpi_recv(nprob,1,MPI_INTEGER,       &
                       ig,17,MPI_COMM_WORLD,status,info )
        call errore ('n2d reduction','info<>0 in recv',info)
        call mpi_recv(psiprob,2*ngper*npol*ngper*npol,MPI_REAL8,  &
                       ig,18,MPI_COMM_WORLD,status,info )         
        call errore ('n2d reduction','info<>0 in recv',info)
        call gramsh(ngper*npol,nprob,1,nprob,         &
                 psibase,psiprob,n2d,epsproj)  
      enddo
    endif
#else
    call gramsh(ngper*npol,nprob,1,nprob,psibase,psiprob,n2d,epsproj) 
#endif
    kin=kin+nteam
  enddo

#ifdef __PARA
  call mpi_barrier( MPI_COMM_WORLD, info )
  call mpi_bcast(n2d,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call errore ('reduction','mpi_bcast 1',info) 
  call mpi_bcast(psibase,2*ngper*npol*ngper*npol,MPI_REAL8,0,  &
                  MPI_COMM_WORLD,info)
  call errore ('reduction','mpi_bcast 1',info)       
#endif

!
! Allocate variables depending on n2d
!
  call allocate_cond_2
  if (npol.eq.2) then
     WRITE( stdout,*) 'ngper, ngper*npol, n2d = ', ngper, ngper*npol, n2d
  else
     WRITE( stdout,*) 'ngper, n2d = ', ngper, n2d
  endif
!
! Construct components of basis vector set on G_per
!
  call DCOPY(2*ngper*npol*n2d,psibase,1,newbg,1)

!
! set and solve the eigenvalue equation for each slab
!
  allocate( amat1( n2d, n2d ) )
  allocate( ymat( ngper*npol, n2d ) )

! for reduced basis set H'_{ab}=e*^i_aH_{ij}e^j_b
  do k=1, nrz
    if(nkofz(k).ne.0) then
      ymat=(0.d0,0.d0)
!
!     First compute y_{ib}=H_{ij}e_{jb}
!
      do ig=1, ngper
         do jg=1, ngper
            do ipol=1, 2
               gp(ipol) = gper(ipol,ig) - gper(ipol,jg)
            enddo
            index=number(gp, at, fftxy, nrx, nry)
            do is=1,npol
               do js=1,npol
                  if (index.gt.0) then
                     aij=vppot(k,index,is,js)
                  else
                     aij=(0.d0,0.d0)
                  endif
                  if ((ig.eq.jg).and.(is.eq.js))          &
                     aij=aij+(gper(1,ig)**2+              &
                              gper(2,ig)**2)*tpiba2
                     amat(ig+(is-1)*ngper,jg+(js-1)*ngper)= aij
               enddo
            enddo
         enddo
      enddo
      call ZGEMM('n','n',ngper*npol,n2d,ngper*npol,one,amat,ngper*npol, &
                         newbg,ngper*npol,zero,ymat,ngper*npol)
!
!     and construct H'_{ab}=<e_a|y_b>
!
      do il=1, n2d
        do jl=il, n2d
          amat1(il,jl)=ZDOTC(ngper*npol,newbg(1,il),1,ymat(1,jl),1)
          amat1(jl,il)=conjg(amat1(il,jl))
        enddo
      enddo
!
!     Solving the eigenvalue problem and construction zk
!
      info=-1
      call hev_ab(n2d, amat1, n2d, zkr(1,nkofz(k)),       &
                  psiper(1,1,nkofz(k)), 0.d0, 0.d0, info)

    endif
  enddo

#ifdef __PARA
  call mpi_barrier( MPI_COMM_WORLD, info )
#endif

!
! saving the 2D data on the file if lwrite_loc=.t. 
!
  if (lwrite_loc) then
    if(fil_loc.eq.' ') call errore ('local','fil_loc no name',1)
    call seqopn(4,fil_loc,'unformatted',exst)
    write(4) n2d
    write(4) ((newbg(ig,il), ig=1, ngper*npol), il=1, n2d)

    write(4) (((psiper(ig,il,k),ig=1,n2d),il=1,n2d), &
                                                k=1,nrzp)
    write(4) ((zkr(il,k),il=1,n2d),k=1,nrzp)
    close(unit=4)
  endif

  deallocate(amat)
  deallocate(amat1)
  deallocate(ymat)
  deallocate(gp)
  deallocate(psibase)
  deallocate(psiprob)
  deallocate(el)
  deallocate(fftxy)

  call stop_clock('local')
  return 
end subroutine local
!-----------------------------------

function number(gp, at, fftxy, nrx, nry)
!
! This function receives as input the coordinates of 2D g vector 
! and write on output its fft position. 
!
  implicit none
  integer :: nrx, nry, fftxy(-(nrx-1)/2:nrx/2, -(nry-1)/2:nry/2), &
             number, n1, n2
  real(kind=kind(0.d0)) :: gp(2), at(3,3), x1, x2 
  real(kind=kind(0.d0)), parameter :: eps=1.d-4

  x1=gp(1)*at(1,1)+gp(2)*at(2,1)
  x2=gp(1)*at(1,2)+gp(2)*at(2,2) 
  n1=nint(x1)
  n2=nint(x2)
  if (n1.le.nrx/2.and.n1.ge.-(nrx-1)/2.and.    &
      n2.le.nry/2.and.n2.ge.-(nry-1)/2) then
    number=fftxy(n1,n2)
  else
!
! The g vector is not inside the 2D mesh
!
    number=-1
  endif

  return
end function number 
      
