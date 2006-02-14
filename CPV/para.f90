!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

! nproc:   number of processors
! me:      number of this processor ( starting from one )
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
!  integer maxproc, ncplanex
!  parameter (maxproc=64, ncplanex=37000)
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
!
!
!----------------------------------------------------------------------
      subroutine read_rho( unit, nspin, rhor )
!----------------------------------------------------------------------
      !
      ! read rhor(nnr,nspin) from file
      !
      use kinds,           only: DP
      USE fft_base,        ONLY: dfftp
      use grid_dimensions, only: nr1, nr2, nr3, nr1x, nr2x, nnr => nnrx
      use xml_io_base,     only: read_rho_xml
      !
      implicit none
      !
      integer  :: unit, nspin
      real(DP) :: rhor( nnr, nspin )
      !
      integer            :: is
      CHARACTER(LEN=256) :: filename
      !
      do is=1,nspin
        IF( nspin == 2 .AND. is == 1 ) THEN
           filename = 'rho.1.xml'
        ELSE IF( nspin == 2 .AND. is == 2 ) THEN
           filename = 'rho.2.xml'
        ELSE
           filename = 'rho.xml'
        END IF
        CALL read_rho_xml( filename, rhor(:,is), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
      END DO
      return
      end subroutine read_rho
!
!----------------------------------------------------------------------
      subroutine write_rho( unit, nspin, rhor )
!----------------------------------------------------------------------
      !
      ! write rho to file
      !
      use kinds,           only: DP
      USE fft_base,        ONLY: dfftp
      use grid_dimensions, only: nr1, nr2, nr3, nr1x, nr2x, nnr => nnrx
      use xml_io_base,     only: write_rho_xml
      !
      implicit none
      !
      integer  :: unit, nspin
      real(DP) :: rhor( nnr, nspin )
      !
      integer            :: is
      CHARACTER(LEN=256) :: filename
      !
      do is=1,nspin
        IF( nspin == 2 .AND. is == 1 ) THEN
           filename = 'rho.1.xml'
        ELSE IF( nspin == 2 .AND. is == 2 ) THEN
           filename = 'rho.2.xml'
        ELSE
           filename = 'rho.xml'
        END IF
        CALL write_rho_xml( filename, rhor(:,is), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
      END DO
      return
      end subroutine write_rho
!
!
!-----------------------------------------------------------------------
      subroutine reduce(size,ps)
!-----------------------------------------------------------------------
!
!     sums a distributed variable s(size) over the processors.
!     This version uses a fixed-length buffer of appropriate (?) size
!
      use parallel_include
      use mp_global, only: nproc
!
      implicit none
      integer size
      real(8)  ps(size)
!
      integer ierr, n, nbuf
      integer, parameter:: MAXB=10000
      real(8) buff(MAXB)
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
      end subroutine reduce
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
      subroutine write_rho_xsf(tau0,h,rho)
!----------------------------------------------------------------------
      
      use ions_base, only: nsp, na, pmass
      use parameters
      use grid_dimensions, only: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr => nnrx
      use io_global, only: ionode
      use mp_global, only: nproc, mpime
      use fft_base, only: dfftp
      use mp, only: mp_bcast

      implicit none

#ifdef __PARA
      include 'mpif.h'
#endif

      integer specie(100)!specie atomica da modificare 
      integer i, j, k, ir, ia, is, ityp(natx), nat00, tp(3)
      integer ir1, ir2, ir3, ip1, ip2, ip3, ipn
      real(8) tau0(3,natx), tau00(3,natx), rho(nnr)
      real(8) h(3,3), l1, l2, l3, shift(3), maxr, minr, cm(3)
      real(8) ll1, ll2, ll3, tot_m
      real(8), allocatable:: rho_aux(:)

      integer :: isa, me

#ifdef __PARA
      integer ip, ierr, incr(nproc), displs(nproc)
      real(8), allocatable:: rhow(:)
#endif

      me = mpime + 1

#ifdef __PARA
! in parallel execution, only the first nodes writes
      if (me.eq.1) then
#endif
      open(40,file='plot.xsf',status='unknown')

      write(40,*) ' DIM-GROUP'
      write(40,*) ' 3 1'
      write(40,*) ' PRIMVEC'
      do i = 1,3
         write(40,'(3(1x,f12.6))') (h(j,i)*0.529177d0,j=1,3)
      end do
      write(40,*) ' PRIMCOORD'

      do ir = 1,nr1*10
         write(26,*) ir, rho(ir)
      end do

#ifdef __PARA
      end if
#endif




      do i=1,100
         specie(i)=mod(10*ityp(i),70)+ityp(i)
      enddo

      specie(1)=1
      specie(2)=78
      specie(3)=8





      l1 = h(1,1) + h(2,1) + h(3,1)
      l2 = h(1,2) + h(2,2) + h(3,2)
      l3 = h(1,3) + h(2,3) + h(3,3)
      ll1 = dsqrt(h(1,1)**2+h(1,2)**2+h(1,3)**2)
      ll2 = dsqrt(h(2,1)**2+h(2,2)**2+h(2,3)**2)
      ll3 = dsqrt(h(3,1)**2+h(3,2)**2+h(3,3)**2)
      nat00 = 0
      tot_m = 0.d0
      do i = 1,3
         cm(i) = 0.d0
      end do
      isa = 0
      do is = 1,nsp
         do ia = 1,na(is)
            tot_m = tot_m + pmass(is)
            nat00 = nat00 + 1
            isa = isa + 1
            do i = 1,3
               tau00(i,nat00) = tau0(i,isa)
               cm(i) = cm(i) + tau00(i,nat00)*pmass(is)
            end do
            ityp(nat00) = is
         end do
      end do
      do i = 1,3
         cm(i) = cm(i)/tot_m
      end do

! to center the plot of the charge density at the center of the unit cell where also
! the center of mass of the system is moved

      shift(1) = 0.5d0*l1 - cm(1)
      tp(1) = nint(shift(1)*DBLE(nr1)/ll1)
      shift(1) = 0.d0 !DBLE(tp(1))*ll1/DBLE(nr1)
      shift(2) = 0.5d0*l2 - cm(2)
      tp(2) = nint(shift(2)*DBLE(nr2)/ll2)
      shift(2) = 0.d0 !DBLE(tp(2))*ll2/DBLE(nr2)
      shift(3) = 0.5d0*l3 - cm(3)
      tp(3) = nint(shift(3)*DBLE(nr3)/ll3)
      shift(3) = 0.d0 !DBLE(tp(3))*ll3/DBLE(nr3)

#ifdef __PARA
! in parallel execution, only the first nodes writes
      if (me.eq.1) then
#endif
      write(40,*) ' ',nat00,' 1'
      do i = 1,nat00
         write(40,'(2x,i2,2x,3(f12.6,1x))')  mod(10*ityp(i),70)+ityp(i), &
     &        (((tau00(j,i)+shift(j))*0.529177d0),j=1,3)
      end do
      write(40,*) ' ATOMS'
      do i = 1,nat00
         write(40,'(2x,i2,2x,3(f12.6,1x))')  mod(10*ityp(i),70)+ityp(i), &
     &        (((tau00(j,i)+shift(j))*0.529177d0),j=1,3)
      end do
      write(40,*) ' BEGIN_BLOCK_DATAGRID3D'
      write(40,*) ' 3D_PWSCF'
      write(40,*) ' DATAGRID_3D_UNKNOWN'
      write(40,*) nr1, nr2, nr3
      write(40,'(3(1x,f10.6))') -0.0d0*l1*0.529177d0,-0.0d0*l2*0.529177d0,  &
     &                         -0.0d0*l3*0.529177d0
      do i = 1,3
         write(40,'(3(1x,f10.6))')  (h(j,i)*0.529177d0,j=1,3)
      end do
#ifdef __PARA
      end if

      if (me.eq.1) allocate(rhow(nr1x*nr2x*nr3x))
      do ip=1,nproc
         incr(ip) = dfftp%nnp * ( dfftp%npp(ip) )
         if (ip.eq.1) then
            displs(ip)=0
         else
            displs(ip)=displs(ip-1) + incr(ip)
         end if
      end do
      call mpi_barrier ( MPI_COMM_WORLD, ierr)
      call mpi_gatherv (rho, incr(me), MPI_REAL8,                       &
     &                  rhow,incr, displs, MPI_REAL8,                   &
     &                     0, MPI_COMM_WORLD, ierr)
      if (ierr.ne.0) call errore('mpi_gatherv','ierr<>0',ierr)

! in parallel execution, only the first nodes writes

      if (me.eq.1) then
#endif

      allocate(rho_aux(nr1x*nr2x*nr3x))

      maxr = 0.d0
      minr = 10.d0
      do ir3 = 1,nr3
         if ((ir3-tp(3)).le.0) then
            ip3 = (ir3-tp(3))+nr3
         else
            ip3 = (ir3-tp(3))
         end if
         do ir2 = 1,nr2
            if ((ir2-tp(2)).le.0) then
               ip2 = (ir2-tp(2))+nr2
            else
               ip2 = (ir2-tp(2))
            end if
            do ir1 = 1,nr1
               if ((ir1-tp(1)).le.0) then
                  ip1 = (ir1-tp(1))+nr1
               else
                  ip1 = (ir1-tp(1))
               end if
               ir = ir1 + (ir2-1)*nr1 + (ir3-1)*nr2*nr1
               ipn = ip1 + (ip2-1)*nr1 + (ip3-1)*nr2*nr1
#ifdef __PARA
               rho_aux(ir) = rhow(ir) !rhow(ipn)
#else
               rho_aux(ir) = rho(ir) !rho(ipn)
#endif
               maxr = max(rho_aux(ir),maxr)
               minr = min(rho_aux(ir),minr)
            end do
         end do
      end do
#ifdef __PARA
      deallocate (rhow)
#endif
      write(6,*) 'minr1 = ', minr
      write(6,*) 'maxr1 = ', maxr

      write(40,'(6e13.5)')                                              &
     &        (((rho_aux((k-1)*nr1*nr2+(j-1)*nr1+i),                    &
     &                                 i=1,nr1),j=1,nr2),k=1,nr3)
      deallocate (rho_aux)
      write(40,*) ' END_DATAGRID_3D'
      write(40,*) ' END_BLOCK_DATAGRID3D'

      close(40)
#ifdef __PARA
      end if
#endif

      return
    end subroutine write_rho_xsf


