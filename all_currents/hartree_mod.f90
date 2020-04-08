MODULE hartree_mod

  USE kinds, ONLY: DP

  SAVE

  CHARACTER(len=256) :: init_linear, file_output, trajdir

  real(kind=DP) ::delta_t 
  complex(kind=DP), allocatable :: evc_uno(:,:), evc_due(:, :) ! TODO: maybe save one allocation of npwx*nb
   !evc_uno will be a copy of last evc
   !evc_due will be a copy of evc of the first electrons() call

END MODULE hartree_mod


