!
! Copyright (C) 2001-2003 PWSCF group
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
  USE io_global,  ONLY : stdout
  use pwcom
  USE kinds, only : DP
  use phcom
  use io_files
  use global_version
  implicit none

  character (len=9) :: cdate, ctime, code = 'PHONON'

  external date_and_tim
  !      call sigcatch( )
  ! use ".false." to disable all clocks except the total cpu time clock
  ! use ".true."  to enable clocks
  !      call init_clocks (.false.)

  call init_clocks (.true.)
  call start_clock ('PHONON')
  gamma_only = .false.
  call startup (nd_nmbr, code, version_number)
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')
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
  if (trans.and..not.recover) call dynmat0

  if (epsil.and.irr0.le.0) then

     WRITE( stdout, '(/,5x," Computing electric fields")')

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

  if (trans) then
     call phqscf
     call dynmatrix
     if (fildrho.ne.' ') call punch_plot_ph
  endif
  if (elph) then
     if (.not.trans) then 
         call dvanqq
         call elphon
     end if
     call elphsum
  endif
  call stop_ph (.true.)
end program phonon
