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
    SUBROUTINE read_rho( nspin, rhor )
!----------------------------------------------------------------------
      !
      ! read rhor(nnr,nspin) from file
      !
      use kinds,           ONLY: DP
      USE fft_base,        ONLY: dfftp
      use grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nnrx
      use xml_io_base,     ONLY: read_rho_xml, restart_dir
      use control_flags,   ONLY: ndr
      USE io_files,        ONLY: outdir
      !
      implicit none
      !
      integer  :: nspin
      real(DP) :: rhor( nnrx, nspin )
      !
      integer            :: is
      CHARACTER(LEN=256) :: filename
      !
      filename = restart_dir( outdir, ndr )
      !
      filename = TRIM(filename) // '/' // 'charge-density'
      !
      CALL read_rho_xml( filename, rhor(:,1), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
      !
      IF( nspin == 2 ) THEN
         !
         filename = TRIM(filename) // '/' // 'spin-polarization'
         !
         CALL read_rho_xml( filename, rhor(:,2), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
         !
         !  Convert rho_tot, spin_pol back to rho_up, rho_down
         !
         rhor(:,2) = 0.5d0 * ( rhor(:,1) - rhor(:,2) )
         rhor(:,1) = rhor(:,1) - rhor(:,2)
         !
      END IF

      RETURN
    END SUBROUTINE read_rho
!
!----------------------------------------------------------------------
      subroutine old_write_rho( unit, nspin, rhor )
!----------------------------------------------------------------------
!
! collect rhor(nnrx,nspin) on first node and write to file
!
      use parallel_include
      use grid_dimensions, only : nr1x, nr2x, nr3x, nnrx
      use gvecw ,          only : ngw
      USE mp_global,       ONLY : me_image, nproc_image, intra_image_comm
      USE io_global,       ONLY : ionode, ionode_id
      USE fft_base,        ONLY : dfftp
      USE mp,              ONLY : mp_barrier
      !
      implicit none
      !
      integer,       INTENT(IN) :: unit, nspin
      real(kind=DP), INTENT(IN) :: rhor(nnrx,nspin)
      !
      integer :: ir, is
      integer :: proc, ierr
      integer, allocatable:: displs(:), recvcount(:)
      real(kind=DP), allocatable:: rhodist(:)
      !
      ALLOCATE( displs( nproc_image ), recvcount( nproc_image ) )
      !
      if (ionode) allocate(rhodist(nr1x*nr2x*nr3x))
      !
      do proc=1,nproc_image
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
         call mp_barrier()
         call mpi_gatherv( rhor(1,is), recvcount(me_image), MPI_DOUBLE_PRECISION,        &
     &                     rhodist,recvcount, displs, MPI_DOUBLE_PRECISION,        &
     &                     ionode_id, intra_image_comm, ierr)
         call errore('mpi_gatherv','ierr<>0',ierr)
#endif
!
! write the charge density to unit "unit" from first node only
!
         if ( ionode ) &
            write( unit, '(F12.7)' ) (rhodist(ir),ir=1,nr1x*nr2x*nr3x)
         !
      end do
      
      DEALLOCATE( displs, recvcount )
      if (ionode) deallocate(rhodist)
!
      return
      end subroutine old_write_rho
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
      use mp_global, only: intra_image_comm
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
      CALL mp_min( nmin0, intra_image_comm )
      CALL mp_max( nmax0, intra_image_comm )
      nmin = nmin0
      nmax = nmax0

      return
      end subroutine nrbounds
