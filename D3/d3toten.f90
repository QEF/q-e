!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program d3toten
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use io
  implicit none
#ifdef __PARA
  include 'mpif.h'
#endif
  character :: cdate * 9, ctime * 9, version * 12
  integer :: nu_i, nu_i0, irecv
  real (8) :: t0, t1, get_clock

  external date_and_tim
  !      call sigcatch( )
  ! use ".false." to disable all clocks except the total cpu time clock
  ! use ".true."  to enable clocks
  !      call init_clocks (.false.)

  call init_clocks (.true.)
  call start_clock ('D3TOTEN')
  version = 'D3TOTEN1.2.0'
#ifdef __PARA
  call startup (nd_nmbr, version)
#else
  nd_nmbr = '   '
  call date_and_tim (cdate, ctime)
  write (6, 9000) version, cdate, ctime

9000 format (/5x,'Program ',a12,' starts ...',/5x, &
       &            'Today is ',a9,' at ',a9)
#endif
  write (6, '(/5x,"UltraSoft (Vanderbilt) ", &
       &                "Pseudopotentials")')
  !
  ! Initialization routines
  !
  call d3_readin
  call allocate_d3
  call d3_setup
  call d3_summary
  call openfild3
  call d3_init
  call show_memory ()
  call print_clock ('D3TOTEN')
  !
  ! Used for testing purposes: if wraux=.true. it writes
  ! different terms of the third derivative matrix in different files.
  !
  if (wraux) call write_aux (1)

  call setv (54 * nat * nat * nat, 0.d0, d3dyn, 1)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  nu_i0 = 1
  if (recv) then
     !
     ! If recv.eq.true. this is a recover run
     !
     call d3_recover (irecv, - 1)

     write (6,  * ) ' Recover Run             index:', irecv
     if (irecv.ge.401.and.irecv.lt.499) then
        nu_i0 = irecv - 400
        goto 304
     else
        goto (301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, &
             312, 313) irecv
     endif
  endif
  !
  ! It calculates the variation of wavefunctions | d/du(q) psi(k) >
  !
  t0 = get_clock ('D3TOTEN')
  if (.not.lgamma) then
     write (6, '(/,5x,"calling gen_dwfc(1)")')
     call gen_dwfc (1)
     call d3_recover (1, + 1)
     t1 = get_clock ('D3TOTEN') - t0
     t0 = get_clock ('D3TOTEN')
     write (6, '(5x,"gen_dwfc(1)   time: ",f12.2, &
          &         " sec    Total time:",f12.2," sec")') t1, t0
  endif
  if (istop.eq.1) stop
  !
  ! It calculates the variation of wavefunctions | d/du(q=0) psi(k) >
  !
301 continue
  write (6, '(/,5x,"calling gen_dwfc(3)")')
  call gen_dwfc (3)
  call d3_recover (2, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"gen_dwfc(3)   time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.2) stop
  !
  ! It calculates the variation of wavefunctions | d/du(q=0) psi(k+q) >
  !  to be used for the terms < dpsi | dpsi ><psi| dH |psi>
  !
302 continue
  if (.not.lgamma) then
     write (6, '(/,5x,"calling gen_dwfc(2)")')
     call gen_dwfc (2)
     call d3_recover (3, + 1)
     t1 = get_clock ('D3TOTEN') - t0
     t0 = get_clock ('D3TOTEN')
     write (6, '(5x,"gen_dwfc(2)   time: ",f12.2, &
          &          " sec    Total time:",f12.2," sec")') t1, t0
  endif
  if (istop.eq.3) stop
  !
  ! It writes on files terms of the type: <dpsi| dH | psi>, that
  ! will be used for the metallic case
  !
303 continue
  write (6, '(/,5x,"calling gen_dpdvp")')
  call gen_dpdvp
  call d3_recover (4, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"gen_dpdvp     time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.4) stop
  !
  ! It calculates the term < dpsi| dH | dpsi >
  !
304 continue
  do nu_i = nu_i0, 3 * nat
     if (q0mode (nu_i) ) then
        write (6, '(/,5x,"calling dpsidvdpsi:",i3)') nu_i
        call dpsidvdpsi (nu_i)
        call d3_recover (401 + nu_i, + 1)
        t1 = get_clock ('D3TOTEN') - t0
        t0 = get_clock ('D3TOTEN')

        write (6, '(5x,"dpsidvdpsi",i3," time: ",f12.2, &
             &   " sec    Total time:",f12.2," sec")') nu_i, t1, t0

        if (istop.gt.400.and.nu_i.eq. (istop - 400) ) stop
     endif
  enddo
  call d3_recover (5, + 1)
  if (istop.eq.5) stop
  !
  ! It calculates the term < dpsi| dpsi > < psi | dH | psi>
  !
305 continue
  write (6, '(/,5x,"calling dpsidpsidv")')
  call dpsidpsidv
  call d3_recover (6, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"dpsidpsidv    time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.6) stop
  !
  ! It calculates the term   drho * d2V
  !
306 continue
  write (6, '(/,5x,"calling drhod2v")')
  call drhod2v
  call d3_recover (7, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"drhod2v       time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.7) stop
  !
  ! It calculates the term   rho * d3V
  !
307 continue
  write (6, '(/,5x,"calling d3vrho")')
  call d3vrho
  call d3_recover (8, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"d3vrho        time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.8) stop
  !
  ! It calculates the contribution due to ionic term
  !
308 continue
  write (6, '(/,5x,"calling d3ionq")')
  call d3ionq (nat, ntyp, ityp, zv, tau, alat, omega, xq, at, bg, g, &
       gg, ngm, gcutm, nmodes, u, ug0, npert_i, npert_f, q0mode, d3dyn)
  call d3_recover (9, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"d3ionq        time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.9) stop
  !
  ! In the metallic case some additional terms are needed
  !
309 continue
  write (6, '(/,5x,"calling d3_valence")')
  call d3_valence
  call d3_recover (10, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"d3_valence    time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (istop.eq.10) stop
  !
  ! drho_cc(+1) adds to the variation or the charge -written on a file-
  ! the variation of the core charge. The variation of the charge,
  ! modified this way is used by the routines d3_exc and d3dyn_cc.
  ! drho_cc(-1) restores drho as it was before (useless)
  !
310 continue
  write (6, '(/,5x,"calling drho_cc(+1)")')
  call drho_cc ( + 1)
  call d3_recover (11, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"drho_cc(+1)   time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! It calculates d3Ei * drho * drho * drho, where drho is the variation
  ! of the charge and d3Ei is the third derivative of the
  ! Kohn-Sham-Energy term depending on the charge density.
  !
311 continue
  write (6, '(/,5x,"calling d3_exc")')
  call d3_exc
  call d3_recover (12, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')
  write (6, '(5x,"d3_exc        time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! It calculates additional terms due to non_linear-core-corrections
  !
312 continue
  write (6, '(/,5x,"calling d3dyn_cc")')
  call d3dyn_cc
  call d3_recover (13, + 1)
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')

  write (6, '(5x,"d3dyn_cc      time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  !
  ! drho is restored as it was before
  !
  !      write(6,'(/,5x,"calling drho_cc(-1)")')
  !      call drho_cc(-1)
  !      t1 = get_clock('D3TOTEN') - t0
  !      t0 = get_clock('D3TOTEN')
  !      write(6,'(5x,"drho_cc(-1)   time: ",f12.2,
  !     +       " sec    Total time:",f12.2," sec")') t1,t0
  if (wraux) call write_aux (2)
  !
  ! Symmetrizes d3dyn, calculates the q in the star and writes the result
  ! for every q on a file.
  !
313 continue
  write (6, '(/,5x,"calling d3matrix")')
  call d3matrix
  t1 = get_clock ('D3TOTEN') - t0
  t0 = get_clock ('D3TOTEN')

  write (6, '(5x,"d3matrix      time: ",f12.2, &
       &         " sec    Total time:",f12.2," sec")') t1, t0
  if (wraux) call write_aux (3)

  call stop_d3 (.true.)
end program d3toten
