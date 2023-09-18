!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------
SUBROUTINE swfc (npw_, nbnd_, vkb_, wfc_, swfc_) 
  !---------------------------------------------------------------
  !! This routine applies the S operator to the function \(\text{wfc_}\)
  !! and puts the result in \(\text{swfc_}\), i.e.
  !! \(\text{swfc_} = S * \text{wfc_}\).
  !
  !! Important notice: here, the global array \(\text{vkb}\) is used as 
  !! a workspace, because the routine \(\texttt{s_psi}\) uses 
  !! \(\text{vkb}\) internally.  
  !! \(\text{vkb_}\) can be a beta function at k or k+q, therefore, 
  !! before changing the global array \(\text{vkb}\), we need to be very 
  !! careful: save \(\text{vkb}\) to a temporary array \(\text{vkb_save}\).
  !! Then we copy \(\text{vkb_save}\) back to \(\text{vkb}\), such that
  !! the meaning of \(\text{vkb}\) outside of this routine is restored
  !! (whatever it is).
  !
  !! Written by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,  ONLY : DP
  USE becmod, ONLY : becp, calbec
  USE uspp,   ONLY : vkb, nkb
  USE wvfct,  ONLY : npwx
  
  IMPLICIT NONE
  
  INTEGER,     INTENT(IN)  :: npw_
  INTEGER,     INTENT(IN)  :: nbnd_
  COMPLEX(DP), INTENT(IN)  :: vkb_ (npwx, nkb)
  COMPLEX(DP), INTENT(IN)  :: wfc_ (npwx, nbnd_) 
  COMPLEX(DP), INTENT(OUT) :: swfc_(npwx, nbnd_) 
  !
  COMPLEX(DP), ALLOCATABLE :: vkb_save(:,:)
  !
  CALL start_clock( 'swfc' )
  !
  swfc_ = (0.d0, 0.d0)
  !
  ALLOCATE(vkb_save(npwx,nkb))
  !
  vkb_save = vkb
  vkb = vkb_
  ! 
  CALL calbec (npw_, vkb, wfc_, becp)
  !
  ! s_psi uses vkb
  !
  CALL s_psi (npwx, npw_, nbnd_, wfc_, swfc_)
  !
  vkb = vkb_save
  ! 
  CALL stop_clock( 'swfc' )
  !
  RETURN
  !
END SUBROUTINE swfc
            
