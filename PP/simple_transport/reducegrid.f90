!     This subroutine generates a reduced grid around the Fermi level
!     for efficient calculation of transport integrals 
!
      subroutine reducegrid (nktot,nks,nbnd,etk,ekq,ef,cut,deg,nkeff,iflag)

      implicit none

      integer, intent(in) :: nktot, nks, nbnd
                             ! nktot : total number of k-points (full grid)
                             ! nks : total number of k-points (IBZ)
                             ! nbnd : number of bands of interest

      double precision, intent(in) :: etk(nbnd,nks),ef,deg,cut
                                      ! etk :  energy eigenvalues
                                      ! ef, deg : Fermi level and smearing (or temperature)
                                      ! cut : scanning parameter (recommended value is ~ 10)

      integer, intent(in) :: ekq(nktot)
                             ! Map of IBZ to full-grid      

      integer, intent(out) :: nkeff(nbnd), iflag(nbnd,nktot)
                              ! nkeff : effective reduced grid size for each band of interest
                              ! iflag : index of the reduced grid point

      ! LOCAL parameters
      integer :: ik,ikk,ibnd, ind_k


      ! Initiate nkeff
      nkeff = 0
      iflag = 0

      do ibnd = 1, nbnd
         !
         ikk = 1
         do ik = 1, nktot
            !
            ind_k = ekq(ik) ! map between IBZ and full grid
            !
            if ( abs(etk(ibnd,ind_k)-ef) < cut*deg) then
               !
               iflag(ibnd,ikk) = ik
               nkeff(ibnd) = nkeff(ibnd) + 1
               ikk = ikk + 1
               !
            end if
         end do ! ik
      end do ! ibnd


      end subroutine reducegrid
