!
! Copyright (C) 2001-2016 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE dnonloccorr( rho, drho, xq, dvaux )
  !------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE fft_base,         ONLY : dfftp
  USE funct,            ONLY : dft_is_gradient, dft_is_nonlocc, get_inlc
  USE ph_vdW_DF,        ONLY : dv_drho_vdwdf
  USE ph_rVV10,         ONLY : dv_drho_rvv10
  !
  implicit none
  !
  real(DP),    INTENT(IN)    :: rho (dfftp%nnr, nspin_mag), xq(3)
  complex(DP), INTENT(IN)    :: drho (dfftp%nnr, nspin_mag)
  complex(DP), INTENT(INOUT) :: dvaux (dfftp%nnr, nspin_mag)
  complex(DP), ALLOCATABLE   :: dvnonloc(:,:)
  integer                    :: inlc


  allocate( dvnonloc(dfftp%nnr, nspin_mag) )
  dvnonloc(:,:) = (0.D0,0.D0)
  inlc = get_inlc()

  IF ( inlc > 0 .AND. inlc < 26 ) THEN
     CALL dv_drho_vdwdf(rho, drho, nspin_mag, xq, dvnonloc)
  ELSEIF ( inlc == 26) THEN
     CALL dv_drho_rvv10(rho, drho, nspin_mag, xq, dvnonloc)
  ENDIF

  dvaux = dvaux + dvnonloc

  DEALLOCATE( dvnonloc )
 
END SUBROUTINE dnonloccorr
