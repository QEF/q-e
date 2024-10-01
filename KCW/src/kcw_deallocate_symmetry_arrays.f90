SUBROUTINE kcw_deallocate_symmetry_arrays()
  !
  ! deallocate arrays related to calculation of symmetries
  !
  USE control_kcw,     ONLY: rir, r, s_w, ft_w, nsym_w, nqstot_ibz
  USE control_kcw,     ONLY: xq_ibz, wq_ibz, ibz2fbz, fbz2ibz
  !
  IMPLICIT NONE
  !
  if( allocated(rir)) deallocate(rir)
  if( allocated(r)) deallocate(r)
  if( allocated(s_w)) deallocate(s_w)
  if( allocated(ft_w)) deallocate(ft_w)
  if( allocated(nsym_w)) deallocate(nsym_w)
  if( allocated(nqstot_ibz)) deallocate(nqstot_ibz)
  if( allocated(xq_ibz)) deallocate(xq_ibz)
  if( allocated(wq_ibz)) deallocate(wq_ibz)
  if( allocated(ibz2fbz)) deallocate(ibz2fbz)
  if( allocated(fbz2ibz)) deallocate(fbz2ibz)
  !
END SUBROUTINE kcw_deallocate_symmetry_arrays
  