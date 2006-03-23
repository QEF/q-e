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
      subroutine read_rho( nspin, rhor )
!----------------------------------------------------------------------
      !
      ! read rhor(nnr,nspin) from file
      !
      use kinds,           ONLY: DP
      USE fft_base,        ONLY: dfftp
      use grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nnr => nnrx
      use xml_io_base,     ONLY: read_rho_xml, restart_dir
      use control_flags,   ONLY: ndr
      USE io_files,        ONLY: scradir
      !
      implicit none
      !
      integer  :: nspin
      real(DP) :: rhor( nnr, nspin )
      !
      integer            :: is
      CHARACTER(LEN=256) :: filename
      !
      filename = restart_dir( scradir, ndr )
      !
      do is=1,nspin
        IF( nspin == 2 .AND. is == 1 ) THEN
           filename = TRIM(filename) // '/' // 'charge-density-up'
        ELSE IF( nspin == 2 .AND. is == 2 ) THEN
           filename = TRIM(filename) // '/' // 'charge-density-dw'
        ELSE
           filename = TRIM(filename) // '/' // 'charge-density'
        END IF
        CALL read_rho_xml( filename, rhor(:,is), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
      END DO
      return
      end subroutine read_rho
!
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
