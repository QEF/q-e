SUBROUTINE read_qlist_ibz()
  USE constants,           ONLY : DP
  USE control_kcw,         ONLY : tmp_dir_kcw
  USE control_kcw,         ONLY : num_wann
  USE klist,               ONLY : xk
  USE control_kcw,         ONLY : nqstot_ibz, fbz2ibz, ibz2fbz, wq_ibz, nqstot
  USE control_kcw,         ONLY : xq_ibz
  USE control_kcw,         ONLY : mp1, mp2, mp3
  USE mp_world, only: mpime
  !
  implicit none
  !
  INTEGER             :: iq, iwann, iq_ibz
  INTEGER             :: iun_qlist_ibz
  character(len=1024) :: filename
  REAL(DP) :: wq
  !
  ALLOCATE( nqstot_ibz( num_wann ) )
  ALLOCATE( xq_ibz( 3, mp1*mp2*mp3, num_wann ) )
  ALLOCATE( wq_ibz( mp1*mp2*mp3, num_wann ) )
  ALLOCATE( ibz2fbz( mp1*mp2*mp3, num_wann))
  ALLOCATE( fbz2ibz( mp1*mp2*mp3, num_wann))  
  !
  
  DO iwann = 1, num_wann
    iun_qlist_ibz = 155 + iwann
    WRITE (filename, "(A,I0.3,A)") TRIM(tmp_dir_kcw)//'qlist_ibz_iwann_', iwann, '.txt'
    OPEN (iun_qlist_ibz, file = filename)
    READ(iun_qlist_ibz,*)  nqstot_ibz( iwann )
    DO iq = 1, nqstot
      !
      !READ(iun_qlist_ibz,*) iq_ibz, &
      !                    wq_ibz(iq_ibz, iwann)
      READ(iun_qlist_ibz,*) iq_ibz, wq
      !
      fbz2ibz(iq, iwann) = iq_ibz
      IF( fbz2ibz(iq, iwann) .ne. -1 ) THEN 
        wq_ibz(iq_ibz, iwann) = wq
        ibz2fbz(iq_ibz, iwann) = iq
        xq_ibz(:, iq_ibz, iwann) = xk(:, iq)
      END IF
    ENDDO!ik 
  END DO!iwann    
  !
  END SUBROUTINE
