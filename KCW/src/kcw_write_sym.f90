SUBROUTINE write_symmetry_op(iwann)
!this function reset s, ft, nsym with the symmetries respect by wf iwann, 
!computed at step wann2kcw
  USE control_kcw,         ONLY : tmp_dir_kcw
  USE symm_base,           ONLY : s, ft, nsym
  USE control_kcw,         ONLY : nsym_w, s_w, ft_w
  !
  implicit none
  !
  INTEGER, intent(in) :: iwann
  INTEGER             :: isym
  INTEGER             :: iun_sym
  character(len=1024) :: filename  
  !
  !
  iun_sym = 555 + iwann
  WRITE (filename, "(A,I0.3,A)") TRIM(tmp_dir_kcw)//'sym_iwann_', iwann, '.txt'
  OPEN (iun_sym, file = filename)
  WRITE(iun_sym,*)  nsym_w(iwann)
  DO isym = 1, nsym_w(iwann)
    WRITE(iun_sym,*) s_w(:,:,isym, iwann)
    WRITE(iun_sym,*) ft_w(:, isym, iwann)
  END DO
  !
END SUBROUTINE

