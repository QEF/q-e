MODULE hartree_mod

   USE kinds, ONLY: DP

   SAVE

   CHARACTER(len=256) :: file_output, trajdir = '' !, init_linear
   real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)
   real(kind=DP), allocatable :: v_cm(:, :)

   real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
   complex(kind=DP), allocatable :: evc_uno(:, :), evc_due(:, :) ! TODO: maybe save one allocation of npwx*nb
   !evc_uno will be a copy of last evc
   !evc_due will be a copy of evc of the first electrons() call
   complex(kind=DP), allocatable :: dvpsi_save(:,:,:) ! to save the solution of the system between iterations
   logical :: save_dvpsi = .true. ! if true dvpsi_save is allocated and used

   integer :: first_step, last_step, step_mul, step_rem
   logical :: restart ! if true try to read last calculated step from output and set first_step
   logical :: subtract_cm_vel ! if true do velocity renormalization
   logical :: ec_test ! activates tests for debugging purposes

END MODULE hartree_mod

