!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine print_clock_ph
!-----------------------------------------------------------------------


use pwcom
use parameters, only : DP
use phcom
implicit none
write (6, * )
call print_clock ('PHONON')
write (6,  * ) '    INITIALIZATION: '
call print_clock ('phq_setup')
call print_clock ('phq_init')
write (6, * )
call print_clock ('phq_init')
if (nlcc_any) call print_clock ('set_drhoc')
call print_clock ('init_vloc')
call print_clock ('init_us_1')
!call print_clock ('init_us_2')
call print_clock ('newd')
call print_clock ('dvanqq')
call print_clock ('drho')
if ((epsil.or.zue).and.okvan) call print_clock ('cmpt_qdipol')

if(epsil) then 
  write (6, * )
  write (6,  * ) '    DIELECTRIC CONSTANT AND EFFECTIVE CHARGES:'
  call print_clock ('solve_e')
  call print_clock ('dielec')
  call print_clock ('zstar_eu')
#ifdef TIMING_ZSTAR_US
  write (6, * )
  call print_clock ('zstar_eu_us')
  call print_clock ('zstar_us_1')
  call print_clock ('zstar_us_2')
  call print_clock ('zstar_us_3')
  call print_clock ('zstar_us_4')
  call print_clock ('zstar_us_5')
#endif
#ifdef TIMING_ADD_DKMDS
  write (6, * )
  call print_clock ('add_dkmds')
  call print_clock ('add_dkmds1')
  call print_clock ('add_dkmds2')
  call print_clock ('add_dkmds3')
  call print_clock ('add_dkmds4')
  call print_clock ('add_dkmds5')
  call print_clock ('add_dkmds6')
#endif
endif
if(trans) then 
  write (6, * )
  write (6,  * ) '    DYNAMICAL MATRIX:'
  call print_clock ('dynmat0')
  call print_clock ('phqscf')
  call print_clock ('dynmatrix')
  write (6, * )
  call print_clock ('phqscf')
  call print_clock ('solve_linter')
  call print_clock ('drhodv')
  if (zue) call print_clock('add_zstar_ue')
  if (zue) call print_clock('add_zstar_1')
  if (zue.and.okvan) call print_clock('add_zstar_us')
endif
write (6, * )
call print_clock ('dynmat0')
call print_clock ('dynmat_us')
call print_clock ('addusdynmat1')
call print_clock ('d2ionq')
if (nlcc_any) call print_clock ('dynmatcc')
write (6, * )
call print_clock ('dynmat_us')
call print_clock ('addusdynmat')
write (6, * )
call print_clock ('phqscf')
call print_clock ('solve_linter')
write (6, * )
call print_clock ('solve_linter')
call print_clock ('dvqpsi_us')
call print_clock ('ortho')
call print_clock ('cgsolve')
call print_clock ('incdrhoscf')
call print_clock ('addusddens')
call print_clock ('vpsifft')
call print_clock ('dv_of_drho')
call print_clock ('mix_pot')
call print_clock ('ef_shift')
call print_clock ('localdos')
#ifdef __PARA
call print_clock ('psymdvscf')
#else
call print_clock ('symdvscf')
#endif
call print_clock ('newdq')
call print_clock ('adddvscf')


call print_clock ('drhodvus')
write (6, * )
call print_clock ('dvqpsi_us')

call print_clock ('dvqpsi_us_on')
write (6, * )
call print_clock ('cgsolve')

call print_clock ('ch_psi')
write (6, * )
call print_clock ('ch_psi')
call print_clock ('first')
call print_clock ('h_psiq')

call print_clock ('last')
write (6, * )
call print_clock ('h_psiq')
call print_clock ('firstfft')
call print_clock ('product')
call print_clock ('secondfft')

call print_clock ('add_vuspsi')
write (6, * )
call print_clock ('incdrhoscf')

call print_clock ('addusdbec')
write (6, * )
call print_clock ('drhodvus')

call print_clock ('addusddort')
write (6, * )
write (6,  * ) '     General routines'
call print_clock ('ccalbec')
call print_clock ('cft3')
call print_clock ('cft3s')
call print_clock ('cinterpolate')
call print_clock ('davcio')
call print_clock ('write_rec')
write (6, * )
#ifdef __PARA
write (6,  * ) '     Parallel routines'
call print_clock ('reduce')
call print_clock ('poolreduce')
#endif
return
end subroutine print_clock_ph
