!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

module para_mod

  USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, &
     fft_dlay_deallocate, fft_dlay_set
  USE fft_base, ONLY: dfftp, dffts

  integer maxproc, ncplanex
  parameter (maxproc=64, ncplanex=37000)
  
  character(len=3) :: node

! node:    node number, useful for opening files

  integer nproc, me, mygroup

! nproc:   number of processors
! me:      number of this processor
!

! parallel fft information for the dense grid
!
! npp:     number of plane per processor                      
! n3:      n3(me)+1 = first  plane on proc. me            
! ncp:     number of (density) columns per proc  
! ncp0:    starting column for each processor
! ncplane: number of columns in a plane                   
! nct:     total number of non-zero columns               
! nnr_:    local fft data size                            
! ipc:     index saying which proc owns columns in a plane
! icpl:    index relating columns and pos. in the plane   
!
! n3 -> dfftp%ipp
! ncplane -> dfftp%nnp
! ncp  -> dfftp%nsp
! ncp0 -> dfftp%iss
! npp  -> dfftp%npp
! ipc  -> dfftp%isind
! icpl -> dfftp%ismap
! nnr_ -> dfftp%nnr
!
!  integer  npp(maxproc), n3(maxproc), ncp(maxproc), ncp0(maxproc), &
!           ncplane, nct, nnr_, ipc(ncplanex), icpl(ncplanex)
!
! parallel fft information for the smooth mesh
!
! npps:    number of plane per processor
! ncps:    number of (density) columns per proc 
! ncpw:    number of (wfs) columns per processor
! ncps0:   starting column for each processor
! ncplanes:number of columns in a plane (smooth)
! ncts:    total number of non-zero columns
! nnrs_:   local fft data size
! ipcs:    saying which proc owns columns in a plane
! icpls:   index relating columns and pos. in the plane 
!
! ncpw -> dffts%ncpw
! n3s -> dffts%ipp
! ncplanes -> dffts%nnp
! ncps  -> dffts%nsp
! ncps0 -> dffts%iss
! npps  -> dffts%npp
! ipcs  -> dffts%isind
! icpls -> dffts%ismap
! nnrs_ -> dffts%nnr
!

contains

  subroutine deallocate_para_mod
    use stick_base, only: sticks_deallocate
    call fft_dlay_deallocate( dfftp )
    call fft_dlay_deallocate( dffts )
    call sticks_deallocate()
  end subroutine

end module para_mod
!
!-----------------------------------------------------------------------
      subroutine startup
!-----------------------------------------------------------------------
!
!  This subroutine initializes the Message Passing environment. 
!  The number of processors "nproc" returned by this routine is 
!  determined by the way the process was started. 
!  This is configuration- and machine-dependent. For inter-
!  active execution it is determined by environment variable MP_PROCS
!
      use para_mod
      use mp, only: mp_start, mp_env
      use io_global, only: stdout
      use global_version
!
      implicit none

      integer ierr
!
      call mp_start()
      call mp_env( nproc, me, mygroup )
      me = me + 1
!
!
! parent process (source) will have me=1 - child process me=2,...,NPROC
! (for historical reasons: MPI uses 0,...,NPROC-1 instead )
!
      if ( nproc > maxproc) &
     &   call errore('startup',' too many processors ',nproc)
!
      if (me < 10) then
         write(node,'(i1,2x)') me
      else if (me < 100) then
         write(node,'(i2,1x)') me
      else if (me < 1000) then
         write(node,'(i3)') me
      else
         call errore('startup','wow, >1000 nodes !!',nproc)
      end if
!
! only the first processor writes
!
      if ( me == 1 ) then
         WRITE( stdout,'(72("*"))')
         WRITE( stdout,'(4("*"),64x,4("*"))')
         WRITE( stdout,'(4("*"),"  CPV: variable-cell Car-Parrinello ", &
              &  "molecular dynamics          ",4("*"))') 
         WRITE( stdout,'(4("*"),"  using ultrasoft Vanderbilt ", &
              &  "pseudopotentials - v.",a6,8x,4("*"))') version_number
         WRITE( stdout,'(4("*"),64x,4("*"))')
         WRITE( stdout,'(72("*"))')
         WRITE( stdout,'(/5x,''Parallel version (MPI)'')')
         WRITE( stdout,'(5x,''Number of processors in use:   '',i4)') nproc
      else
         open(6,file='/dev/null',status='unknown')
!
! useful for debugging purposes 
!         open(6,file='out.'//node,status='unknown')
      end if
!
      return
      end


!
!----------------------------------------------------------------------
      subroutine read_rho(unit,nspin,rhor)
!----------------------------------------------------------------------
!
! read from file rhor(nnr,nspin) on first node and distribute to other nodes
!
      use para_mod
      use parallel_include
      use grid_dimensions, only: nr1x, nr2x, nr3x, nnr => nnrx
      implicit none
      integer unit, nspin
      real(kind=8) rhor(nnr,nspin)
!
      integer ir, is
      integer root, proc, ierr, n, displs(nproc), sendcount(nproc)
      real(kind=8), allocatable:: rhodist(:)
!
!
      if (me.eq.1) allocate(rhodist(nr1x*nr2x*nr3x))
      root = 0
      do proc=1,nproc
         sendcount(proc) =  dfftp%nnp * ( dfftp%npp(proc) )
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + sendcount(proc-1)
         end if
      end do
      do is=1,nspin
!
! read the charge density from unit "unit" on first node only
!
         if (me.eq.1) read(unit) (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
!
! distribute the charge density to the other nodes
!
#if defined __PARA
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         call mpi_scatterv(rhodist, sendcount, displs, MPI_DOUBLE_PRECISION,       &
     &                     rhor(1,is),sendcount(me),   MPI_DOUBLE_PRECISION,       &
     &                     root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('mpi_scatterv','ierr<>0',ierr)
#endif
!
! just in case: set to zero unread elements (if any)
!
         do ir=sendcount(me)+1,nnr
            rhor(ir,is)=0.d0
         end do
      end do
      if (me.eq.1) deallocate(rhodist)
!
      return
      end subroutine read_rho
!
!----------------------------------------------------------------------
      subroutine write_rho(unit,nspin,rhor)
!----------------------------------------------------------------------
!
! collect rhor(nnr,nspin) on first node and write to file
!
      use para_mod
      use parallel_include
      use grid_dimensions, only: nr1x, nr2x, nr3x, nnr => nnrx
      use gvecw , only : ngw
      implicit none
      integer unit, nspin
      real(kind=8) rhor(nnr,nspin)
!
      integer ir, is
      integer root, proc, ierr, displs(nproc), recvcount(nproc)
      real(kind=8), allocatable:: rhodist(:)
!
!
      if (me.eq.1) allocate(rhodist(nr1x*nr2x*nr3x))
!
      root = 0
      do proc=1,nproc
         recvcount(proc) =  dfftp%nnp  * ( dfftp%npp(proc) )
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
      end do
!
      do is=1,nspin
!
! gather the charge density on the first node
!

#if defined __PARA
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         call mpi_gatherv (rhor(1,is), recvcount(me), MPI_DOUBLE_PRECISION,        &
     &                     rhodist,recvcount, displs, MPI_DOUBLE_PRECISION,        &
     &                     root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('mpi_gatherv','ierr<>0',ierr)
#endif
!
! write the charge density to unit "unit" from first node only
!
         if (me.eq.1) write(unit) (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
         ! if (me.eq.1) write(unit,'(f12.7)') (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
      end do
      if (me.eq.1) deallocate(rhodist)
!
      return
      end subroutine write_rho
!
!
!----------------------------------------------------------------------
      subroutine cfftpb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign)
!----------------------------------------------------------------------
!
!   not-so-parallel 3d fft for box grid, implemented only for sign=1
!   G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!   The array f (overwritten on output) is NOT distributed:
!   a copy is present on each processor.
!   The fft along z  is done on the entire grid.
!   The fft along xy is done only on planes that have components on the
!   dense grid for each processor. Note that the final array will no
!   longer be the same on all processors.
!
      use para_mod
      use grid_dimensions, only: nr3
      use fft_scalar, only: cft_b
!
      implicit none
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign
      complex(kind=8) f(nr1bx*nr2bx*nr3bx)
!
      integer ir3, ibig3, imin3, imax3, np3
!
      call parabox(nr3b,irb3,nr3,imin3,imax3)
      np3=imax3-imin3+1
! np3 is the number of planes to be transformed
      if (np3.le.0) return
      call cft_b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,imin3,imax3,sign)
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine parabox(nr3b,irb3,nr3,imin3,imax3)
!----------------------------------------------------------------------
!
! find if box grid planes in the z direction have component on the dense
! grid on this processor, and if, which range imin3-imax3
!
      use para_mod
! input
      integer nr3b,irb3,nr3
! output
      integer imin3,imax3
! local
      integer ir3, ibig3
!
      imin3=nr3b
      imax3=1
      do ir3=1,nr3b
         ibig3=1+mod(irb3+ir3-2,nr3)
         if(ibig3.lt.1.or.ibig3.gt.nr3)                                 &
     &        call errore('cfftpb','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         if (ibig3.gt.0.and.ibig3.le.dfftp%npp(me)) then
            imin3=min(imin3,ir3)
            imax3=max(imax3,ir3)
         end if
      end do
!
      return
      end
!
!-----------------------------------------------------------------------
      subroutine reduce(size,ps)
!-----------------------------------------------------------------------
!
!     sums a distributed variable s(size) over the processors.
!     This version uses a fixed-length buffer of appropriate (?) size
!
      use para_mod
      use parallel_include
!
      implicit none
      integer size
      real(kind=8)  ps(size)
!
      integer ierr, n, nbuf
      integer, parameter:: MAXB=10000
      real(kind=8) buff(MAXB)
!
      if (nproc.le.1) return
      if (size.le.0) return
      call start_clock( 'reduce' )
!
!  syncronize processes
!
#if defined __PARA
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) call errore('reduce','error in barrier',ierr)
!
      nbuf=size/MAXB
!
      do n=1,nbuf
         call mpi_allreduce (ps(1+(n-1)*MAXB), buff, MAXB,              &
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0)                                                 &
     &        call errore('reduce','error in allreduce1',ierr)
         call DCOPY(MAXB,buff,1,ps(1+(n-1)*MAXB),1)
      end do
!
!    possible remaining elements < maxb
!
      if (size-nbuf*MAXB.gt.0) then
          call mpi_allreduce (ps(1+nbuf*MAXB), buff, size-nbuf*MAXB,    &
     &          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
          if (ierr.ne.0)                                                &
     &         call errore('reduce','error in allreduce2',ierr)
          call DCOPY(size-nbuf*MAXB,buff,1,ps(1+nbuf*MAXB),1)
      endif
#endif

      call stop_clock( 'reduce' )
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine nrbounds(ngw,nr1s,nr2s,nr3s,mill,nmin,nmax)
!----------------------------------------------------------------------
!
! find the bounds for (i,j,k) indexes of all wavefunction G-vectors
! The (i,j,k) indexes are defined as: G=i*g(1)+j*g(2)+k*g(3)
! where g(1), g(2), g(3) are basis vectors of the reciprocal lattice
!
      use parallel_include
      use mp, only: mp_min, mp_max
      implicit none
! input
      integer ngw,nr1s,nr2s,nr3s,mill(3,*)
! output
      integer nmin(3), nmax(3)
! local
      integer nmin0(3), nmax0(3), ig, ierr
!
!
      nmin0(1)=  nr1s
      nmax0(1)= -nr1s
      nmin0(2)=  nr2s
      nmax0(2)= -nr2s
      nmin0(3)=  nr3s
      nmax0(3)= -nr3s
!
      do ig=1,ngw
         nmin0(1) = min(nmin0(1),mill(1,ig))
         nmin0(2) = min(nmin0(2),mill(2,ig))
         nmin0(3) = min(nmin0(3),mill(3,ig))
         nmax0(1) = max(nmax0(1),mill(1,ig))
         nmax0(2) = max(nmax0(2),mill(2,ig))
         nmax0(3) = max(nmax0(3),mill(3,ig))
      end do
!
! find minima and maxima for the FFT box across all nodes
!
      CALL mp_min( nmin0 )
      CALL mp_max( nmax0 )
      nmin = nmin0
      nmax = nmax0

      return
      end subroutine nrbounds

!----------------------------------------------------------------------
      subroutine write_pot(unit,rhos2)
!     - To write the hartree potential
!        M.S
!----------------------------------------------------------------------
!
! collect rhos2(nnrs) on first node and write to file
!
      use para_mod
      use smooth_grid_dimensions , nnrs => nnrsx
      use parallel_include

      implicit none

      integer unit, nspin
      real(kind=8) rhos2(nnrs)
!
      integer ir, is
      integer root, proc, ierr, displs(nproc), recvcount(nproc)
      real(kind=8), allocatable:: rhodist(:)
!
!
      if (me.eq.1) allocate(rhodist(nr1sx*nr2sx*nr3sx))
!
      root = 0
      do proc=1,nproc
         recvcount(proc) =   dffts%nnp * dffts%npp(proc)
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
      end do
!
!      do is=1,nspin
!
! gather the charge density on the first node
#if defined __PARA
         call mpi_barrier ( MPI_COMM_WORLD, ierr)
         call mpi_gatherv (rhos2, recvcount(me), MPI_DOUBLE_PRECISION,        &
     &                     rhodist,recvcount, displs, MPI_DOUBLE_PRECISION,        &
     &                     root, MPI_COMM_WORLD, ierr)
         if (ierr.ne.0) call errore('mpi_gatherv','ierr<>0',ierr)
#endif
!
! write the charge density to unit "unit" from first node only
!
         if (me.eq.1) write(unit,'(f12.6)') (rhodist(ir),ir=1,nr1sx*nr2sx*nr3sx)
!      end do
      if (me.eq.1) deallocate(rhodist)
!
      return
      end subroutine write_pot
