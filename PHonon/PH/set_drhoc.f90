!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_drhoc (q, drc)
  !---------------------------------------------------------------------
  !! Calculate the Fourier transform of the core charge for all pseudo
  !! without structure factor and put it in \(\text{drc}\), at \(\text{q}\)
  !! point used to calculate derivatives of the core charge.
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba2
  USE gvect,     ONLY : g, ngm
  USE ions_base, ONLY : ntyp => nsp
  USE uspp_param,ONLY : upf
  USE uspp,      ONLY : nlcc_any
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(in) :: q(3)
  !! the q-point used for structure factor
  REAL(DP),INTENT(inout) :: drc(ngm,ntyp)
  !! Fourier-transform of core charge at q+G
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: qg2(:)
  !! (q+G)^2 in (2\pi/a)^2 units
  INTEGER :: ng, nt
  ! counters

  IF ( .not. nlcc_any ) RETURN

  CALL start_clock('set_drhoc')
  !
  ALLOCATE ( qg2(ngm) )
  do ng = 1, ngm
     qg2(ng) = ( g(1,ng) + q(1) )**2 + ( g(2,ng) + q (2) )**2 + &
               ( g(3,ng) + q(3) )**2
  end do
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        call interp_rhc( nt, ngm, qg2, tpiba2, drc(1,nt) )
     else
        drc (:,nt) = 0.0_dp
     end if
  end do
  !
  DEALLOCATE( qg2 )
  CALL stop_clock('set_drhoc')

  RETURN
END SUBROUTINE set_drhoc
