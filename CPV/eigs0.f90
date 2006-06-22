!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine eigs0( ei, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nudx )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds, only            : DP
      use io_global, only        : stdout
      use constants, only        : autoev
      use parallel_toolkit, only : dspev_drv
      USE sic_module, only       : self_interaction

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin)
      real(DP), intent(in) :: lambda( nudx, nudx, nspin ), f( nx )
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: lambdar(:), wr(:)
      real(DP) zr(1)
      integer :: iss, j, i, ierr, k, n, nspin_eig, npaired
      logical :: tsic
!
      tsic = ( ABS( self_interaction) /= 0 )

      IF( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      END IF

      do iss = 1, nspin_eig

         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF

         allocate( lambdar( n * ( n + 1 ) / 2 ), wr(n) )

         k = 0

         do i = 1, n
            do j = i, n
               k = k + 1
               lambdar( k ) = lambda( j, i, iss )
            end do
         end do

         CALL dspev_drv( 'N', 'L', n, lambdar, wr, zr, 1 )

         if( lf ) then
            do i = 1, n
               if ( f(iupdwn(iss)-1+i).gt.1.e-6) then
                  wr(i)=wr(i)/f(iupdwn(iss)-1+i)
               else
                  wr(i)=wr(i)/2.0d0 * nspin  ! fake occupation factor to print empty states
               end if
            end do
         end if
         !
         !     store eigenvalues
         !
         IF( SIZE( ei, 1 ) < nupdwn(iss) ) &
            CALL errore( ' eigs0 ', ' wrong dimension array ei ', 1 )

         ei( 1:n, iss ) = wr( 1:n )

         IF( tsic ) THEN
            !
            !  store unpaired state
            !
            ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
            ei( nupdwn(1), 1 ) = lambda( nupdwn(1), nupdwn(1), 1 )
            !
         END IF

         ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
         ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug

         deallocate( lambdar, wr )

      end do
      !
      !
      do iss = 1, nspin

         IF( tsic .AND. iss == 2 ) THEN
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
         END IF

         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0., 0., 0.
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )

         ENDIF

      end do

      IF( tprint ) WRITE( stdout,*)

   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
!
      return
      end subroutine eigs0
