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
  use io_files
  USE io_global, ONLY: ionode_id
  USE mp, ONLY: mp_bcast
  USE ions_base,  ONLY : nat
  USE kinds, only : DP
  USE parser,           ONLY : int_to_char
  USE control_flags, ONLY : iswitch,  restart, lphonon, tr2, &
       mixing_beta, lscf, ldisp, david, isolve
  use phcom
  use global_version
  implicit none

  integer :: iq, iq_start, iustat

  integer :: nks_start
  ! number of initial k points

  real(kind = dp), dimension(:), allocatable :: wk_start
  ! initial weight of k points
  real(kind = dp), dimension(:,:), allocatable :: xk_start
  ! initial coordinates of k points

  LOGICAL :: exst

  character (len=9) :: cdate, ctime, code = 'PHONON'
  CHARACTER (LEN=80) :: auxdyn
  
  character (len=256) :: filname

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
  !
  !  Checking the status of the calculation
  !
  iustat = 98
  filname = trim(tmp_dir) // trim(prefix) //'.stat'
  inquire(file=filname, exist= exst)
  if(exst) then
     open (unit = iustat, file = filname)
     read(iustat, *) iq_start
     write(stdout, "(//5x' STARTING FROM AND OLD RUN')")
     write(stdout, "(5x' Doing now the calculation for q point nr',i3)")  &
          iq_start
  else
     open (unit = iustat, file = filname)
     iq_start = 1
  end if
  !
  if (ldisp) then
     !
     ! Calculate the q-points for the dispersion
     !
     call q_points
     !
     ! Store the name of the matdyn file in auxdyn
     !
     auxdyn = fildyn
     !
     ! Save the starting k points 
     !
     nks_start = nks
     if(.not. allocated(xk_start)) allocate(xk_start(3,nks_start))
     if(.not. allocated(wk_start)) allocate(wk_start(nks_start))
     xk_start(:,1:nks_start) = xk(:,1:nks_start)
     wk_start(1:nks_start) = wk(1:nks_start)
     !
     ! do always a non scf calculation
     !
     lnscf = .true.
  else
     nqs = 1
  end if

  if (lnscf) CALL start_clock('PWSCF')

  do iq = iq_start, nqs

     REWIND (IUSTAT)
     write(iustat,*) iq

     if (ldisp) then
        !
        ! set the name for the output file
        !
        fildyn = trim(auxdyn) // TRIM( int_to_char( iq ))
        !
        ! set the q point
        !
        xqq(1:3) = x_q(1:3,iq)
        xq(1:3)  = x_q(1:3,iq)
        lgamma = xqq (1) .eq.0.d0.and.xqq (2) .eq.0.d0.and.xqq (3) .eq.0.d0
        !
        ! in the case of an insulator one has to calculate 
        ! the dielectric constant and the Born eff. charges
        !
        if (lgamma .and. degauss.eq.0.d0) then
           epsil = .true.
           zue = .true.
        end if
        !
        ! for q != 0 no calculation of the dielectric tensor 
        !            and Born eff. charges
        !
        if(.not. lgamma) then
           epsil = .false.
           zue = .false.
        end if

        call mp_bcast (epsil, ionode_id)
        call mp_bcast (zue, ionode_id)
        call mp_bcast (lgamma, ionode_id)

        nks = nks_start
        xk(:,1:nks_start) = xk_start(:,1:nks_start)
        wk(1:nks_start) = wk_start(1:nks_start)

     end if
     !
     !  In the case of q != 0, we make first an non selfconsistent run
     !
     if (.not. lgamma .and. lnscf) then

        WRITE(stdout,'("Calculation of q = "3f8.4)') XQQ

        call clean_pw(.false.)
        CALL close_files()
     
        !
        !  Setting the values for the nscf run
        !
        lphonon   = .TRUE.
        lscf = .false.
     
        restart        = .FALSE.
        restart_bfgs   = .FALSE.
        startingconfig = 'input'

        startingpot = 'file'
        startingwfc = 'atomic'
        tr2 = 1.0e-8_dp
     

        if(.not. allocated(force)) allocate(force(3,nat))
        !
        ! Set the value for the david
        !
        if (isolve == 0)  david = 4

        CALL init_run()
        
        call electrons()
     
        CALL hinit1()

        call sum_band
    
        call close_files()

     end if
     !
     ! Setting nksq
     !
     if (lgamma) then
        nksq = nks
     else
        nksq = nks / 2
     endif
     !
     ! Calculation of the dispersion: do all modes 
     !
     maxirr = 0
     !
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
     !
     ! cleanup of the variables
     !
     call clean_pw(.false.)
     call deallocate_phq()
     !
     !  Close the files
     !
     call close_phq(.true.)

  end do

  close (unit = iustat, status='delete')

  if (allocated(xk_start)) deallocate (xk_start)
  if (allocated(wk_start)) deallocate (wk_start)
  
  if (lnscf) call print_clock_pw()

  call stop_ph (.true.)
end program phonon
