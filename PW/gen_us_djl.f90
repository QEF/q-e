!
! Copyright (C) 2010 Layla Martin-Samos
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine gen_us_djl (ik,npw,djl,size_tab,vec_tab, spline_ps, vec_tab_d2y)
  !----------------------------------------------------------------------
  !
  !  Calculates the kleinman-bylander pseudopotentials with the
  !  derivative of the spherical harmonics projected on vector u
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : tpiba
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : igk
  USE gvect,      ONLY : g
  USE us,         ONLY : nqx, dq
  USE splinelib,  ONLY : splint_deriv
  USE uspp_param, ONLY : upf
  !
  implicit none
  !
  real(DP), intent(inout) ::djl(1:npw)
  integer, intent(in) :: ik, npw
  integer, intent(in) :: size_tab
  real(DP), intent(in) :: vec_tab(1:size_tab)
  real(DP), optional, intent(in) :: vec_tab_d2y(1:size_tab)
  logical :: spline_ps
  !
  integer :: i0, i1, i2, &
       i3, ig
  real(DP), allocatable :: gk(:,:), q (:)
  real(DP) :: px, ux, vx, wx
  complex(DP), allocatable :: sk (:)

  integer :: iq
  real(DP), allocatable :: xdata(:)
  real(DP) :: qt


  allocate ( gk(3,npw) )
  allocate ( q(npw) )

  do ig = 1, npw
     gk (1, ig) = xk (1, ik) + g (1, igk (ig) )
     gk (2, ig) = xk (2, ik) + g (2, igk (ig) )
     gk (3, ig) = xk (3, ik) + g (3, igk (ig) )
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  do ig = 1, npw
     q (ig) = sqrt ( q(ig) ) * tpiba
  end do

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  ! calculate beta in G-space using an interpolation table
  do ig = 1, npw
    qt = sqrt(q(ig)) * tpiba
    if (spline_ps) then
        djl(ig) = splint_deriv(xdata, vec_tab(:), &
                                vec_tab_d2y(:), qt)
    else
        px = qt / dq - int (qt / dq)
        ux = 1.d0 - px
        vx = 2.d0 - px
        wx = 3.d0 - px
        i0 = qt / dq + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        djl (ig) = vec_tab (i0) * (-vx*wx-ux*wx-ux*vx) / 6.d0 + &
                          vec_tab (i1) * (+vx*wx-px*wx-px*vx) / 2.d0 - &
                          vec_tab (i2) * (+ux*wx-px*wx-px*ux) / 2.d0 + &
                          vec_tab (i3) * (+ux*vx-px*vx-px*ux) / 6.d0
    endif
  enddo

  deallocate (q)
  deallocate ( gk )
  if (spline_ps) deallocate(xdata)
  return
end subroutine gen_us_djl
