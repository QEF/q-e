      subroutine cryst_to_cart (nvec, vec, trmat, iflag)
      !-----------------------------------------------------------------------
      !
      !     This routine transforms the atomic positions or the k-point
      !     components from crystallographic to cartesian coordinates 
      !     ( iflag=1 ) and viceversa ( iflag=-1 ).
      !     Output cartesian coordinates are stored in the input ('vec') array
      !
      !
      implicit none
      !
      integer, intent(in) :: nvec, iflag
      ! nvec:  number of vectors (atomic positions or k-points)
      !        to be transformed from crystal to cartesian and vice versa
      ! iflag: gives the direction of the transformation
      double precision, intent(in) :: trmat (3, 3)
      ! trmat: transformation matrix
      ! if iflag=1:
      !    trmat = at ,  basis of the real-space lattice,       for atoms   or
      !          = bg ,  basis of the reciprocal-space lattice, for k-points
      ! if iflag=-1: the opposite
      double precision, intent(inout) :: vec (3, nvec)
      ! coordinates of the vector (atomic positions or k-points) to be
      ! transformed - overwritten on output
      !
      !    local variables
      !
      integer :: nv, kpol
      ! counter on vectors
      ! counter on polarizations
      double precision :: vau (3)
      ! workspace
      !
      !     Compute the cartesian coordinates of each vectors
      !     (atomic positions or k-points components)
      !
      do nv = 1, nvec
         if (iflag.eq.1) then
            do kpol = 1, 3
               vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
                    * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
            end do
         else
            do kpol = 1, 3
               vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
                    * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
            end do
         endif
         do kpol = 1, 3
            vec (kpol, nv) = vau (kpol)
         end do
      end do
      !
      return

      end subroutine cryst_to_cart

