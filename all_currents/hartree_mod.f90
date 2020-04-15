MODULE hartree_mod

  USE kinds, ONLY: DP

  SAVE

  CHARACTER(len=256) :: init_linear, file_output, trajdir
   real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)

  real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
  complex(kind=DP), allocatable :: evc_uno(:,:), evc_due(:, :) ! TODO: maybe save one allocation of npwx*nb
   !evc_uno will be a copy of last evc
   !evc_due will be a copy of evc of the first electrons() call

  integer :: first_step, last_step

END MODULE hartree_mod


