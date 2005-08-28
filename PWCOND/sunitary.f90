!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine sunitary(nchanl, nchanr, smat, info)
!
! It performs a check of unitarity of the scattering matrix
! S = {R_ij, T_ij}
! \sum R_ki^* R_kj + \sum T_ki^* R_kj = \delta_ij
!
  USE kinds, ONLY : DP
  USE io_global,  ONLY :  stdout
implicit none

  integer ::  info, i, j, k
  integer ::  nchanl, &  ! number of channels in the left tip
              nchanr     !     ------------         right tip
  real(DP), parameter :: eps=1.d-4
  real(DP) :: raux, raux1
  complex(DP) ::   &
              s, smat(nchanl+nchanr, nchanl)    ! S matrix

  info = 0
  do i = 1, nchanl
    do j = 1, nchanl
      s = 0.d0
      do k = 1, nchanl
       s = s + CONJG(smat(k,i))*smat(k,j)
      enddo
      do k = 1, nchanr
       s = s + CONJG(smat(nchanl+k,i))*smat(nchanl+k,j)
      enddo
      raux = sqrt(DBLE(s)**2 + AIMAG(s)**2)
      raux1 = raux
      if(i.eq.j) raux1 = abs(raux1-1.d0)
      if(raux1.gt.eps) then
        write(stdout, '(2i5,f12.7)') i, j, raux
        info = 1
      endif
    enddo
  enddo

  return
end subroutine sunitary

