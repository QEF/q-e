!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
         subroutine intref(lam,e,mesh,grid,vpot,ze2,chi)
!---------------------------------------------------------------
!
!  numerical integration of the radial schroedinger equation 
!  computing logarithmic derivatives
!  thresh dermines the absolute accuracy for the eigenvalue
!
!
!
      use kinds, only : DP
      use radial_grids, only: radial_grid_type, series
      implicit none
      type(radial_grid_type), intent(in):: grid
      integer :: &
              mesh,    &     ! the mesh size
              lam           ! the angular momentum
      real(DP) :: &
              vpot(mesh), &   ! the local potential
              chi(mesh),  &   ! the solution
              ze2,       &   ! the nuclear charge in Ry units
              e             ! the eigenvalue

      integer :: &
              ierr,  &      ! used to control allocation
              n             ! generic counter 

      real(DP) :: &
              lamsq,   &     ! combined angular momentum
              b(0:3),c(4)      ! used for starting guess of the solution 

      real(DP),allocatable :: &
              al(:)      ! the known part of the differential equation

      if (mesh.gt.grid%mesh) &
          call errore('intref','mesh dimension is too large',1)

      allocate(al(mesh),stat=ierr)

      do n=1,4
         al(n)=vpot(n)-ze2/grid%r(n)
      enddo
      call series(al,grid%r,grid%r2,b)


      lamsq=(lam+0.5_DP)**2
!
!     b) find the value of solution s in the first two points
!
      call start_scheq( lam, e, b, grid, ze2, chi )

      do n=1,mesh
         al(n)=( (vpot(n)-e)*grid%r2(n) + lamsq )*grid%dx**2/12.0_DP
         al(n)=1.0_DP-al(n)
      enddo
!
!     Integrate forward the equation:
!     c) integrate the equation from 0 to matching radius
!
      do n=2,mesh-1
         chi(n+1)=((12.0_DP-10.0_DP*al(n))*chi(n) -al(n-1)*chi(n-1))/al(n+1)
      enddo
!
!     set the correct r dependence
!            
      do n=1,mesh
         chi(n)= chi(n)*grid%sqr(n)
      enddo


      deallocate(al)
      return
      end subroutine intref
