!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine eigs0( tprint, nspin, nx, nupdwn, iupdwn, f, lambda )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use io_global, only        : stdout
      use constants, only        : au
      use electrons_module, only : ei
      implicit none
! input
      logical, intent(in) :: tprint
      integer, intent(in) :: nspin, nx, nupdwn(nspin), iupdwn(nspin)
      real(8), intent(in) :: lambda(nx,nx), f(nx)
! local variables
      real(8), allocatable :: lambdar(:,:)
      real(8) wr(nx), fv1(nx),fm1(2,nx), zr
      integer iss, j, i, ierr
!
      do iss=1,nspin
         allocate( lambdar( nupdwn(iss), nupdwn(iss) ) )
         do i=1,nupdwn(iss)
            do j=1,nupdwn(iss)
               lambdar(j,i)=lambda(iupdwn(iss)-1+j,iupdwn(iss)-1+i)
            end do
         end do
         call rs( nupdwn(iss), nupdwn(iss), lambdar, wr, 0, zr, fv1, fm1, ierr)
         do i=1,nupdwn(iss)
            if (f(iupdwn(iss)-1+i).gt.1.e-6) then
               wr(i)=wr(i)/f(iupdwn(iss)-1+i)
            else
               wr(i)=0.0
            end if
         end do
         !
         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0., 0., 0.
            WRITE( stdout,14) (wr(i)*au,i=1,nupdwn(iss))

         ELSE
            !
            !     store eigenvalues
            !
            IF( SIZE( ei, 1 ) < nupdwn(iss) ) &
               CALL errore( ' eigs0 ', ' wrong dimension array ei ', 1 )

            ei( 1:nupdwn(iss), 1, iss ) = wr( 1:nupdwn(iss) )

         END IF


         deallocate( lambdar )

      end do

      IF( tprint ) WRITE( stdout,*)

   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
!
      return
      end subroutine eigs0
