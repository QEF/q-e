!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_dealloc_q()
!------------------------------------------------------------
  !
  !  Deallocates variables allocated in hp_allocate_q 
  !  (and in some other routines)
  !
  USE noncollin_module,    ONLY : m_loc
  USE becmod,              ONLY : deallocate_bec_type
  USE uspp,                ONLY : okvan
  USE qpoint,              ONLY : eigqts, ikks, ikqs
  USE lrus,                ONLY : becp1
  USE gc_lr,               ONLY : grho, gmag, dvxc_rr, dvxc_sr, dvxc_ss, &
                                  & dvxc_s, vsgga, segni
  USE eqv,                 ONLY : dmuxc, dpsi, dvpsi, evq
  USE control_lr,          ONLY : lgamma, nbnd_occ
  USE ldaU_lr,             ONLY : swfcatomk, swfcatomkpq, vh_u_save, vh_uv_save
  USE qpoint_aux,          ONLY : ikmks, ikmkmqs, becpt
  USE lr_nc_mag,           ONLY : deeq_nc_save
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  IF (lgamma) THEN
     if (associated(evq))  nullify(evq)
  ELSE
     !$acc exit data delete(evq)
     if (associated(evq))  deallocate(evq)
  ENDIF
  !
  if (allocated(dvpsi))     deallocate (dvpsi)
  if (allocated(dpsi))      deallocate ( dpsi)
  if (allocated(dmuxc))     deallocate (dmuxc)
  if (allocated(nbnd_occ))  deallocate (nbnd_occ)
  if (allocated(ikks))      deallocate (ikks)
  if (allocated(ikqs))      deallocate (ikqs)
  if (allocated(ikmks))     deallocate (ikmks)
  if (allocated(ikmkmqs))   deallocate (ikmkmqs)
  if (allocated(m_loc))     deallocate (m_loc)
  !
  if (allocated(eigqts)) deallocate (eigqts)
  !
  if (allocated(becp1))  then
     do ik=1,size(becp1)
        call deallocate_bec_type ( becp1(ik) )
     enddo
     deallocate(becp1)
  endif
  if (allocated(becpt)) then  
      do ik=1,size(becpt)
         call deallocate_bec_type ( becpt(ik) )
      enddo
      deallocate(becpt)
  endif
  if (allocated(deeq_nc_save)) deallocate(deeq_nc_save)
  !
  if (allocated(vh_u_save))     deallocate (vh_u_save)
  if (allocated(vh_uv_save))    deallocate (vh_uv_save)
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
  if (allocated(swfcatomk))       deallocate (swfcatomk)
  !
  if (lgamma) then
     if (associated(swfcatomkpq)) nullify (swfcatomkpq)
  else
     if (associated(swfcatomkpq)) deallocate (swfcatomkpq)
  endif
  !
  RETURN
  !
END SUBROUTINE hp_dealloc_q
