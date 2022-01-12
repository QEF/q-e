!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE shepardmod
implicit none
save
  !
  integer, parameter :: dp = selected_real_kind(14,200)  
  !
  integer :: PMetric   ! metric for the (inverse) distance in the Shepard formula
  !
  real(dp) :: ScaleSphere  ! scaling factor for the radius of the sphere for the modified method
  !
CONTAINS
!----------------------------------------------------------------------------
subroutine shepard(iwhat, Nb, Nq, q, eq, Nk, k, ek)
!
! compute the band structure with Shepard's interpolation
! iwhat = 1 ... basic Shepard method
!         2 ... modified method with the sphere radius
!
implicit none
  integer, intent(in) :: iwhat 
  integer, intent(in) :: Nq, Nk, Nb
  real(dp), intent(in) :: q(3,Nq), k(3,Nk), eq(Nq,Nb)
  real(dp), intent(out) :: ek(Nk,Nb)
  ! local variables
  real(dp) :: w, d, dsum, esum, dthr, R
  integer :: ib, iq, ik
  !
  if (iwhat.ne.1.and.iwhat.ne.2) then
    write(*,*) 'wrong iwhat in shepard'
    stop
  else
    write(*,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,'(A)') 'Shepard interpolation method'
    write(*,'(4(A,I5))') 'iwhat: ',iwhat, ' Nb: ',Nb, ' Nq: ',Nq, ' Nk: ',Nk 
  end if 
  !
  PMetric =  6  ! set a metric for distance
  R = ScaleSphere * sqrt( (q(1,1)-q(1,2))**2 + (q(2,1)-q(2,2))**2 + (q(3,1)-q(3,2))**2 ) ! the radius of the search sphere is proportional to the uniform grid spacing
  dthr = 0.00010d0**PMetric
  !
  ek = 0.0d0
  !
  do ib = 1, Nb
    do ik = 1, Nk
      !
      dsum = 0.0d0
      esum = 0.0d0
      do iq = 1, Nq
        d = abs(k(1,ik)-q(1,iq))**PMetric + abs(k(2,ik)-q(2,iq))**PMetric + abs(k(3,ik)-q(3,iq))**PMetric
        if(d.gt.dthr) then
          if(iwhat.eq.1) then 
            ! basic Shepard method
            w = 1.0d0/d
          elseif(iwhat.eq.2) then
            ! search only inside the sphere R
            w = (max(0.0d0, (R-d))/(R*d))**2
          end if
          dsum = dsum + w
          esum = esum + w * eq(iq,ib)
        else
          ek(ik,ib) = eq(iq,ib)
          write(*,*) ib, ik, iq, ' found', d
          go to 10  
        end if
      end do 
      ek(ik,ib) = esum / dsum
      !
10    continue
      !
    end do      
  end do 
  !
  return
  !
end subroutine shepard
!----------------------------------------------------------------------------
END MODULE
!----------------------------------------------------------------------------
