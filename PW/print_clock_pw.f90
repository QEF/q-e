!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

subroutine print_clock_pw
   !
   ! this routine prints out the clocks at the end of the run
   ! it tries to construct the calling tree of the program..

   use pwcom
   implicit none

   write (6, * )
   call print_clock ('PWSCF')
   call print_clock ('init_run')
   call print_clock ('electrons')
   if (lforce) call print_clock ('forces')
   if (lstres) call print_clock ('stress')

   write (6, * )
   call print_clock ('electrons')
   call print_clock ('c_bands')
   call print_clock ('sum_band')
   call print_clock ('v_of_rho')
   call print_clock ('newd')
   if (imix.ge.0) then
      call print_clock ('mix_rho')
   else
      call print_clock ('mix_pot')
   endif
   write (6, * )
   call print_clock ('c_bands')
   call print_clock ('init_us_2')
   if (isolve.eq.0) then
      call print_clock ('cegterg')
   elseif (isolve.eq.1) then
      call print_clock ('ccgdiagg')
   else
      call print_clock ('cdiisg')
   endif
   write (6, * )
   call print_clock ('sum_band')
   call print_clock ('sumbec')

   call print_clock ('addusdens')
   write (6, * )
   if (isolve.eq.0) then
      call print_clock ('cegterg')
   elseif (isolve.eq.1) then
      call print_clock ('ccdiagg')
   else
      call print_clock ('cdiisg')
   endif
   if (isolve.eq.0.or.isolve.eq.2) then
      call print_clock ('h_psi')
      call print_clock ('g_psi')
      if (loverlap) then
         call print_clock ('overlap')
         call print_clock ('cdiaghg')
      else
         call print_clock ('cdiagh')
         call print_clock ('cgramg1')
      endif
      call print_clock ('update')
      call print_clock ('last')
      write (6, * )
      call print_clock ('h_psi')
      call print_clock ('init')
      call print_clock ('firstfft')
      call print_clock ('secondfft')
      call print_clock ('add_vuspsi')
      call print_clock ('s_psi')
   else
      call print_clock ('h_1psi')
      call print_clock ('s_1psi')
      write (6, * )
      call print_clock ('h_1psi')
      call print_clock ('init')
      call print_clock ('firstfft')
      call print_clock ('secondfft')
      call print_clock ('add_vuspsi')
   endif
   write (6, * )
   write (6, * ) '     General routines'
   call print_clock ('ccalbec')
   call print_clock ('cft3')
   call print_clock ('cft3s')
   call print_clock ('interpolate')
   call print_clock ('davcio')
   write (6, * )
#ifdef PARA
   write (6,  * ) '     Parallel routines'
   call print_clock ('reduce')
   call print_clock ('fft_scatter')
!   call print_clock('poolreduce')
#endif
   return

end subroutine print_clock_pw

