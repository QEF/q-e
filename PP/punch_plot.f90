!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine punch_plot (filplot, plot_num, sample_bias, z, dz, &
     stm_wfc_matching, emin, emax, kpoint, kband, spin_component, lsign)
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output several quantities
  !     in a real space 3D mesh which can be plotted using ch.x
  !     The integer variable plot_num is used to chose the output quantity
  !
  !           plot_num                  quantity
  !              0                 self consistent charge density
  !              1                 the total potential V_bare+V_H + V_xc
  !              2                 the local ionic pseudopotential
  !              3                 the local density of states at e_fermi
  !              4                 the local density of electronic entropy
  !              5                 stm images
  !              6                 spin polarisation (rho(up)-rho(down)
  !              7                 square of a wavefunction
  !              8                 electron localization function (ELF)
  !              9                 planar averages of each wavefunction
  !             10                 integrated local dos from emin to emax
  !             11                 the V_bare + V_H potential
  !             12                 The electric field potential
  !
  !     The output quantity is written (formatted) on file filplot.
  !
#include "machine.h"

  use pwcom
  USE io_global,  ONLY : stdout
  
#ifdef __PARA
  use para
#endif
  implicit none
  character(len=*) :: filplot
  integer :: kpoint, kband, spin_component, plot_num
  logical :: stm_wfc_matching, lsign
  real(kind=DP) :: sample_bias, z, dz
  real(kind=DP) :: emin, emax, wf, charge

  integer :: is, ik, ibnd, ir, ninter
#ifdef __PARA
  ! auxiliary vector (parallel case)
  real(kind=DP), allocatable :: raux1 (:)

#endif
  ! auxiliary vector
  real(kind=DP), allocatable :: raux (:)
  real(kind=DP), allocatable :: averag (:,:,:), plan (:,:,:)


  if (filplot == ' ') return
#ifdef __PARA
  allocate (raux1( nrx1 * nrx2 * nrx3))    
#endif

  WRITE( stdout, '(/5x,"Calling punch_plot, plot_num = ",i3)') plot_num
  if ( ( plot_num == 8 .OR. plot_num == 9 ) .AND. gamma_only) &
       call errore('punch_plot', &
      ' gamma_only not implemented for this plot ',1)
  !
  allocate (raux( nrxx))    
  !
  !     Here we decide which quantity to plot
  !
  if (plot_num == 0) then
     !
     !      plot of the charge density
     !
     if (spin_component == 0) then
        call DCOPY (nrxx, rho (1, 1), 1, raux, 1)
        do is = 2, nspin
           call DAXPY (nrxx, 1.d0, rho (1, is), 1, raux, 1)
        enddo
     else
        if (nspin == 2) current_spin = spin_component
        call DCOPY (nrxx, rho (1, current_spin), 1, raux, 1)
        call DSCAL (nrxx, 0.5d0 * nspin, raux, 1)
     endif

  elseif (plot_num == 1) then
     !
     !       The total self-consistent potential V_H+V_xc on output
     !
     if (spin_component == 0) then
        call DCOPY (nrxx, vr (1, 1), 1, raux, 1)
        do is = 2, nspin
           call DAXPY (nrxx, 1.0d0, vr (1, is), 1, raux, 1)
        enddo
        call DSCAL (nrxx, 1.d0 / nspin, raux, 1)
     else
        if (nspin == 2) current_spin = spin_component
        call DCOPY (nrxx, vr (1, current_spin), 1, raux, 1)
     endif
     call DAXPY (nrxx, 1.0d0, vltot, 1, raux, 1)

  elseif (plot_num == 2) then
     !
     !       The local pseudopotential on output
     !
     call DCOPY (nrxx, vltot, 1, raux, 1)

  elseif (plot_num == 3) then
     !
     !       The local density of states at e_fermi on output
     !
     call local_dos (1, lsign, kpoint, kband, emin, emax, raux)

  elseif (plot_num == 4) then
     !
     !       The local density of electronic entropy on output
     !
     call local_dos (2, lsign, kpoint, kband, emin, emax, raux)

  elseif (plot_num == 5) then

     call work_function (wf)
#ifdef __PARA
     call stm (wf, sample_bias, z, dz, stm_wfc_matching, raux1)
#else
     call stm (wf, sample_bias, z, dz, stm_wfc_matching, raux)
#endif
     if (stm_wfc_matching) then
        write (title, '("Matching z = ",f4.2," dz in a.u. = ", &
             &       f4.2," Bias in eV = ",f10.4," # states",i4)') &
             z, dz * alat, sample_bias * rytoev, nint (wf)
     else
        write (title, '("No matching, scf wave-functions. ", &
             &       " Bias in eV = ",f10.4," # states",i4)') &
             sample_bias * rytoev, nint (wf)
     endif

  elseif (plot_num == 6) then
     !
     !      plot of the spin polarisation
     !
     if (nspin == 2) then
        call DCOPY (nrxx, rho (1, 1), 1, raux, 1)
        call DAXPY (nrxx, - 1.d0, rho (1, 2), 1, raux, 1)
     else
        raux(:) = 0.d0
     endif

  elseif (plot_num == 7) then

     call local_dos (0, lsign, kpoint, kband, emin, emax, raux)

  elseif (plot_num == 8) then

     call do_elf (raux)

  elseif (plot_num == 9) then

     allocate (averag( nat, nbnd, nkstot))    
     allocate (plan(nr3, nbnd, nkstot))    
     call plan_avg (averag, plan, ninter)

  elseif (plot_num == 10) then

     call local_dos (3, lsign, kpoint, kband, emin, emax, raux)

  elseif (plot_num == 11) then

     raux(:) = vltot(:) 
     if (nspin == 2) then
        rho(:,1) =  rho(:,1) +  rho(:,2)
        nspin = 1
     end if
     call v_h (rho(1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
       nl, ngm, gg, gstart, nspin, alat, omega, ehart, charge, raux)
     if (tefield.and.dipfield) call add_efield(raux)
  elseif (plot_num == 12) then
     raux=0.d0
     if (tefield) then
         call add_efield(raux)
     else
         call errore('punch_plot','e_field is not calculated',-1)
     endif
  else

     call errore ('punch_plot', 'plot_num not implemented', - 1)

  endif
#ifdef __PARA
  if (.not. (plot_num == 5 .or. plot_num == 9) ) &
       call gather (raux, raux1)
  if (me == 1.and.mypool == 1) call plot_io (filplot, title, nrx1, &
       nrx2, nrx3, nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, &
       dual, ecutwfc, plot_num, atm, ityp, zv, tau, raux1, + 1)
  deallocate (raux1)
#else

  call plot_io (filplot, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num, &
       atm, ityp, zv, tau, raux, + 1)
#endif
  if (plot_num == 9) then
#ifdef __PARA
     if (me == 1.and.mypool == 1) then
#endif
        write (4, '(3i8)') ninter, nkstot, nbnd
        do ik = 1, nkstot
           do ibnd = 1, nbnd
              write (4, '(3f15.9,i5)') xk (1, ik) , xk (2, ik) , xk (3, &
                   ik) , ibnd
              write (4, '(4(1pe17.9))') (averag (ir, ibnd, ik) , ir = 1, &
                   ninter)
              do ir = 1, nr3
                 write (4, * ) ir, plan (ir, ibnd, ik)
              enddo
           enddo
        enddo
#ifdef __PARA
     endif
#endif
     deallocate (plan)
     deallocate (averag)
  endif

  deallocate (raux)
  return
end subroutine punch_plot
