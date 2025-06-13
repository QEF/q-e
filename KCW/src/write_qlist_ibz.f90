SUBROUTINE write_qlist_ibz()
  USE control_kcw,         ONLY : tmp_dir_kcw
  USE control_kcw,         ONLY : num_wann, spin_component
  USE io_global,           ONLY : ionode
  USE klist,               ONLY : nkstot, xk
  USE lsda_mod,            ONLY : lsda, isk
  USE control_kcw,         ONLY : nqstot_ibz, fbz2ibz, wq_ibz, nqstot
  !
  implicit none
  !
  INTEGER             :: ik, iwann, iq_ibz
  INTEGER             :: iun_qlist_ibz
  character(len=1024) :: filename  
  !
  DO iwann = 1, num_wann
    iun_qlist_ibz = 155 + iwann
    WRITE (filename, "(A,I0.3,A)") TRIM(tmp_dir_kcw)//'qlist_ibz_iwann_', iwann, '.txt'
    OPEN (iun_qlist_ibz, file = filename)
    IF (ionode) THEN 
      WRITE(iun_qlist_ibz,*)  nqstot_ibz( iwann )
      DO ik = 1, nqstot
        !IF (lsda .AND. isk(ik) /= spin_component) CYCLE
        !
        IF( fbz2ibz(ik, iwann) .eq. -1 ) THEN
          WRITE(iun_qlist_ibz, *) -1, -1
        ELSE 
          iq_ibz = fbz2ibz(ik, iwann)
          WRITE(iun_qlist_ibz, *) fbz2ibz(ik, iwann), &
                                  wq_ibz(iq_ibz, iwann)
        END IF
      ENDDO!ik 
    ENDIF
  END DO!iwann    
  !
END SUBROUTINE
