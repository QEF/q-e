!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine d3_readin
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq. T
  !    input is read from unit 5. A namelist is used on the machine which
  !    allows it. A second routine readfile reads the variables saved
  !    on the filpun file by the self-consistent program.
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use io
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: ios, ipol, iter, na, it, ii
  ! integer variable for I/O control
  ! counter on polarizations
  ! counter on iterations
  ! counter on atoms
  ! counter on types
  ! counter
  ! eigenvalues convergence thresho

  namelist / inputph / ethr_ph, amass, iverbosity, tmp_dir, filpun, &
       fildyn, fildrho, fild0rho, q0mode_todo, wraux, recv, istop, &
       testflag, testint, testreal
  ! atomic masses
  ! write control
  ! directory for temporary files
  ! computed
  ! the punch file produced by pwsc
  ! the file with the dynamical mat
  ! the file with the deltarho
  ! the file with q=0 deltarho
  ! list of the q=0 modes to be com
  ! .true.==> writes some auxiliary
  ! .true.==> this is a recover run
  ! to stop the program at a given
  ! variables used for testing purp
  ! variables used for testing purp
  ! variables used for testing purp
#ifdef __PARA

  if (me.ne.1) goto 400
#endif
  !
  !    Read the first line of the input file
  !
  read (5, '(a)', err = 100, iostat = ios) title_ph
100 call errore ('d3_readin', 'reading title ', abs (ios) )
  !
  !   set default values for variables in namelist
  !
  ethr_ph = 1.d-5
  iverbosity = 0
  tmp_dir = './'
  filpun = ' '
  fildyn = 'd3dyn'
  fildrho = ' '
  fild0rho = ' '
  do ii = 1, 300
     q0mode_todo (ii) = 0
  enddo
  wraux = .false.
  recv = .false.
  istop = 0
  do ii = 1, 50
     testflag (ii) = .false.


  enddo
  !
  !     reading the namelist inputph
  !
#ifdef CRAYY
  !   The Cray does not accept "err" and "iostat" together with a namelist
  read (5, inputph)
  ios = 0
#else
  !
  !   Note: for AIX machine (xlf compiler version 3.0 or higher):
  !   The variable XLFRTEOPTS must be set to "namelist=old"
  !   in order to have "&end" to end the namelist
  !
  read (5, inputph, err = 200, iostat = ios)
#endif

200 call errore ('d3_readin', 'reading inputph namelist', abs (ios) )
  !
  !     Check all namelist variables
  !
  if (ethr_ph.le.0.d0) call errore (' d3_readin', ' Wrong ethr_ph ', &
       1)
  if (iverbosity.ne.0.and.iverbosity.ne.1) call errore ('d3_readin', ' Wrong &
       &iverbosity ', 1)
  if (fildyn.eq.' ') call errore ('d3_readin', ' Wrong fildyn ', 1)
  if (filpun.eq.' ') call errore ('d3_readin', ' Wrong filpun ', 1)
  if (fildrho.eq.' ') call errore ('d3_readin', ' Wrong fildrho ', 1)
  if (fild0rho.eq.' ') call errore ('d3_readin', ' Wrong fild0rho ', &
       1)
  !
  !    reads the q point
  !
  read (5, *, err = 300, iostat = ios) (xq (ipol), ipol = 1, 3)
300 call errore ('d3_readin', 'reading xq', abs (ios) )

  lgamma = xq (1) .eq.0.d0.and.xq (2) .eq.0.d0.and.xq (3) .eq.0.d0
#ifdef __PARA
400 continue

  call bcast_d3_input
  call init_pool
#endif
  call DCOPY (3, xq, 1, xqq, 1)
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file
  !
  if (lgamma) then
     nksq = nks
  else
     nksq = nks / 2
  endif
  !
  if (lsda) call errore ('d3_readin', 'lsda not implemented', 1)
  if (okvan) call errore ('d3_readin', 'US not implemented', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  call allocate_part
#ifdef __PARA

  if (me.ne.1.or.mypool.ne.1) goto 800
#endif
#ifdef __PARA

800 continue
#endif

  if (iswitch.ne. - 2.and.iswitch.ne. - 3.and.iswitch.ne. - &
       4.and..not.lgamma) call errore ('d3_readin', ' Wrong iswitch ', 1 + &
       abs (iswitch) )
  do it = 1, ntyp
     if (amass (it) .le.0.d0) call errore ('d3_readin', 'Wrong masses', &
          it)

  enddo
  if (mod (nks, 2) .ne.0.and..not.lgamma) call errore ('d3_readin', &
       'k-points are odd', nks)
  !
  ! q0mode, and q0mode_todo are not allocated dynamically. Their
  ! dimension is fixed to 300
  !

  if (3 * nat.gt.300) call errore ('d3_readin', 'wrong dimension of q &
       &0mode variable', 1)
  do ii = 1, 3 * nat
     if (q0mode_todo (ii) .gt.3 * nat) call errore ('d3_readin', ' wrong &
          & q0mode_todo ', 1)

  enddo
  return
end subroutine d3_readin
