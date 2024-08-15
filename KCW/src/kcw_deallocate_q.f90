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
  USE qpoint,              ONLY : eigqts, ikks, ikqs, nksq
  USE lrus,                ONLY : becp1
  USE gc_lr,               ONLY : grho, gmag, dvxc_rr, dvxc_sr, dvxc_ss, &
                                  & dvxc_s, vsgga, segni
  USE eqv,                 ONLY : dmuxc, dpsi, dvpsi, evq
  USE control_lr,          ONLY : lgamma, nbnd_occ
  USE qpoint_aux,          ONLY : ikmks, ikmkmqs, becpt, alphapt
  !
  IMPLICIT NONE
  INTEGER :: ik, ipol
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
  IF(ALLOCATED(ikmks))      DEALLOCATE(ikmks)
  IF(ALLOCATED(ikmkmqs))    DEALLOCATE(ikmkmqs)
  if (allocated(m_loc))     deallocate (m_loc)
  !
  if (allocated(eigqts)) deallocate (eigqts)
  if (allocated(becp1))  then
     do ik=1,size(becp1)
        call deallocate_bec_type ( becp1(ik) )
     enddo
     deallocate(becp1)
  endif
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
  IF (ALLOCATED(alphapt)) THEN
     DO ik=1,nksq
        DO ipol=1,3
           CALL deallocate_bec_type ( alphapt(ipol,ik) )
        ENDDO
     ENDDO
     DEALLOCATE (alphapt)
  ENDIF
  IF (ALLOCATED(becpt))  THEN
     DO ik=1, nksq
        CALL deallocate_bec_type ( becpt(ik) )
     ENDDO
     DEALLOCATE(becpt)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE kcw_deallocate_q
