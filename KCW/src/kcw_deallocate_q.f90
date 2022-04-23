!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE kcw_deallocate_q()
!------------------------------------------------------------
  !
  !  Deallocates variables allocated in hp_allocate_q 
  !  (and in some other routines)
  !
  USE noncollin_module,    ONLY : m_loc
  USE becmod,              ONLY : deallocate_bec_type, becp
  USE uspp,                ONLY : okvan
  USE qpoint,              ONLY : eigqts, ikks, ikqs
  USE lrus,                ONLY : becp1
  USE gc_lr,               ONLY : grho, gmag, dvxc_rr, dvxc_sr, dvxc_ss, &
                                  & dvxc_s, vsgga, segni
  USE eqv,                 ONLY : dmuxc, dpsi, dvpsi, evq
  USE control_lr,          ONLY : lgamma, nbnd_occ
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  IF (lgamma) THEN
     if (associated(evq))  nullify(evq)
  ELSE
     if (associated(evq))  deallocate(evq)
  ENDIF
  !
  if (allocated(dvpsi))     deallocate (dvpsi)
  if (allocated(dpsi))      deallocate ( dpsi)
  if (allocated(dmuxc))     deallocate (dmuxc)
  if (allocated(nbnd_occ))  deallocate (nbnd_occ)
  if (allocated(ikks))      deallocate (ikks)
  if (allocated(ikqs))      deallocate (ikqs)
  if (allocated(m_loc))     deallocate (m_loc)
  !
  IF (okvan) THEN 
     if (allocated(eigqts)) deallocate (eigqts)
     if (allocated(becp1))  then
        do ik=1,size(becp1)
           call deallocate_bec_type ( becp1(ik) )
        enddo
        deallocate(becp1)
     endif
  ENDIF
  !
  CALL deallocate_bec_type ( becp )
  !
  ! GGA-specific arrays
  !
  if (allocated(dvxc_rr))         deallocate (dvxc_rr)
  if (allocated(dvxc_sr))         deallocate (dvxc_sr)
  if (allocated(dvxc_ss))         deallocate (dvxc_ss)
  if (allocated(dvxc_s))          deallocate (dvxc_s)
  if (allocated(grho))            deallocate (grho)
  if (allocated(segni))           deallocate (segni)
  if (allocated(vsgga))           deallocate (vsgga)
  if (allocated(gmag))            deallocate (gmag)
  !
  RETURN
  !
END SUBROUTINE kcw_deallocate_q
