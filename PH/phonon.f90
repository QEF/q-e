!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program phonon
  !-----------------------------------------------------------------------
  !
  !   This is the main driver of the phonon program. It controls
  !   the initialization routines and the self-consistent cycle.
  !   At the end of the self-consistent run the dynamical matrix is
  !   computed. In the case q=0 the dielectric constant and the effective
  !   charges are computed.
  !

  use pwcom
  use parameters, only : DP
  use phcom
  use io
  implicit none

  character :: cdate * 9, ctime * 9, version * 12

  external date_and_tim
  !      call sigcatch( )
  ! use ".false." to disable all clocks except the total cpu time clock
  ! use ".true."  to enable clocks
  !      call init_clocks (.false.)

  call init_clocks (.true.)
  call start_clock ('PHONON')
  version = 'PHONON 1.2.0'
#ifdef PARA
  call startup (nd_nmbr, version)
#else
  nd_nmbr = '   '
  call date_and_tim (cdate, ctime)
  write (6, 9000) version, cdate, ctime
9000 format (/5x,'Program ',a12,' starts ...',/5x, &
       &            'Today is ',a9,' at ',a9)
#endif
  write (6, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')
  !
  !   and begin with the initialization part
  !
  call phq_readin

  call allocate_phq
  call phq_setup
  call phq_recover
  call phq_summary

  call openfilq

  call phq_init
  call show_memory ()

  call print_clock ('PHONON')
  if (epsil.and.irr0.le.0) then

     write (6, '(/,5x," Computing electric fields")')

     call solve_e
     if (convt) then
        !
        ! calculate the dielectric tensor epsilon
        !
        call dielec
        !
        ! calculate the effective charges Z(E,Us) (E=scf,Us=bare)
        !
        call zstar_eu
        if (fildrho.ne.' ') call punch_plot_e
     else
        call stop_ph (.false.)
     endif

  endif

  if (trans.and.irr0.le.0) call dynmat0
  if (trans) then
     call phqscf
     call dynmatrix
     if (fildrho.ne.' ') call punch_plot_ph
  endif
  if (elph) then
     if (.not.trans) call elphon
     call elphsum
  endif
  call stop_ph (.true.)
end program phonon
