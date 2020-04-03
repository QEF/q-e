MODULE hartree_mod

  USE kinds, ONLY: DP

  SAVE

  CHARACTER(len=256) :: prefix_uno, prefix_due, prefix_thermal, &
       init_linear, file_output, file_dativel

  integer :: blocco,passo,calcolato

  real(kind=DP) ::delta_t 
  

END MODULE hartree_mod


