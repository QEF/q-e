!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine print_clock_d3  

  use d3com
  implicit none
  write (6, * )  
  call print_clock ('D3TOTEN')  
  call print_clock ('d3_setup')  

  call print_clock ('phq_init')  
  write (6, * )  
  call print_clock ('solve_linter')  
  call print_clock ('ortho')  
  call print_clock ('cgsolve')  
  call print_clock ('incdrhoscf')  
  call print_clock ('dv_of_drho')  
#ifdef PARA
  call print_clock ('psymdvscf')  
  call print_clock ('psymd0rho')  
#else

  call print_clock ('symdvscf')  
#endif
  write (6, * )  
  call print_clock ('cgsolve')  

  call print_clock ('ch_psi')  
  write (6, * )  
  call print_clock ('ch_psi')  
  call print_clock ('h_psiq')  

  call print_clock ('last')  
  write (6, * )  
  call print_clock ('h_psiq')  
  call print_clock ('firstfft')  
  call print_clock ('product')  

  call print_clock ('secondfft')  
  write (6, * )  
  write (6,  * ) '     General routines'  
  call print_clock ('ccalbec')  
  call print_clock ('cft3')  
  call print_clock ('cft3s')  
  call print_clock ('cinterpolate')  
  call print_clock ('davcio')  
  write (6, * )  
#ifdef PARA
  write (6,  * ) '     Parallel routines'  
  call print_clock ('reduce')  
  call print_clock ('poolreduce')  
#endif
  return  
end subroutine print_clock_d3
