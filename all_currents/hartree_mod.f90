MODULE hartree_mod

  USE kinds, ONLY: DP

  SAVE

  CHARACTER(len=256) :: prefix_uno, prefix_due, &
       init_linear, file_output, file_dativel, thermodir

  real(kind=DP) ::delta_t 
  complex(kind=DP), allocatable :: evc_uno(:,:), evc_due(:, :)
   !evc_uno will be last evc, that contains the FIRST position in the input
   !evc_due will be a copy of evc of the first electrons() call

END MODULE hartree_mod


