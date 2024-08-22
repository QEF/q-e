SUBROUTINE read_symmetry_op(iwann, nsym_aux, s_aux, ft_aux)
!this function reset s, ft, nsym with the symmetries respect by wf iwann, 
!computed at step wann2kcw
  USE kinds,               ONLY : DP
  USE control_kcw,         ONLY : tmp_dir_kcw
  !
  implicit none
  !
  INTEGER, intent(in) :: iwann
  INTEGER, INTENT(OUT) :: nsym_aux
  INTEGER, INTENT(OUT) :: s_aux(3,3,48)
  REAL(DP), INTENT(OUT):: ft_aux(3,48)
  INTEGER             :: isym
  INTEGER             :: iun_sym
  character(len=1024) :: filename  
  !
  !
  iun_sym = 555 + iwann
  WRITE (filename, "(A,I0.3,A)") TRIM(tmp_dir_kcw)//'sym_iwann_', iwann, '.txt'
  OPEN (iun_sym, file = filename)
  READ(iun_sym,*)  nsym_aux
  DO isym = 1, nsym_aux
    READ(iun_sym,*) s_aux(:,:,isym)
    READ(iun_sym,*) ft_aux(:, isym)
  END DO
  CLOSE(iun_sym)
  !
END SUBROUTINE


SUBROUTINE reset_symmetry_op(iwann)
!this function reset all the quantities needed from symm_base
  USE kinds,               ONLY : DP
  USE symm_base,           ONLY : invsym, s, ft, nsym, inverse_s, nrot
  USE symm_base,           ONLY : s_axis_to_cart
  USE symm_base,           ONLY : copy_sym, time_reversal
  !USE symm_base,           ONLY : is_group
  USE io_global,            ONLY : stdout, ionode
  USE control_kcw,         ONLY : nsym_old, mp1, mp2, mp3
  USE cell_base,          ONLY: bg
  USE klist,                ONLY : xk, wk, nkstot
  USE start_k,              ONLY : nks_start, xk_start, wk_start
  !
  implicit none
  !
  INTEGER, intent(in) :: iwann
  INTEGER             :: isym, jsym, ik
  INTEGER             :: nsym_aux
  INTEGER             :: s_aux(3,3,48)
  REAL(DP)            :: ft_aux(3,48)
  LOGICAL             :: sym(48)
  INTEGER             :: t_rev_eff(48)
  !
  CALL read_symmetry_op(iwann, nsym_aux, s_aux, ft_aux)
  !
  ! find index of the symmetry on the total symmetries
  !
  sym = .false.
  nsym = nsym_old
  DO isym = 1, nsym
    DO jsym = 1, nsym_aux
      IF ( ALL ( s(:,:,isym) - s_aux(:,:, jsym) .eq. 0 ) &
    .and. SUM( ABS( ft(:,isym) - ft_aux(:, jsym ) ) ) .lt. 1.D-03 ) THEN
        sym(isym) = .true.
        EXIT
      END IF 
    END DO
  END DO
  !
  ! from now on, 
  ! adapted from symm_base/find_sym
  !
  !
  ! ... Here we re-order all rotations in such a way that true sym.ops
  ! are the first nsym; rotations that are not sym.ops. follow
  nsym = copy_sym( nrot, sym )
  WRITE(stdout,'(8X, "SYM : number of symmetry for iwann =", I5, " :", I5, 3x "(out of ", I5, " )")') iwann, nsym, nrot
  !
  !IF ( .NOT. is_group( nsym ) ) THEN
  !  CALL infomsg( 'find_sym', 'not a group! symmetry disabled' )
  !  nsym = 1
  !END IF
  !
  ! ... check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = ALL( s(:,:,nsym/2+1) == -s(:,:,1) )
  !
  CALL inverse_s()
  !
  CALL s_axis_to_cart()
  !
  ! calculate IBZ of k points
  !
  !t_rev_eff=0 
  !CALL kpoint_grid ( nsym, time_reversal, .FALSE., s, t_rev_eff, &
  !bg, mp1*mp2*mp3, 0.D0, 0.D0, 0.D0, mp1,mp2,mp3, nkstot, xk, wk)
  !nks_start = nkstot
  !xk_start(:,1:nkstot) = xk(:,1:nkstot) 
  !wk_start(1:nkstot) = wk(1:nkstot)  
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*) 
  !WRITE(*,*) ( xk(:,ik), ik=1, nkstot)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !WRITE(*,*)
  !
END SUBROUTINE
