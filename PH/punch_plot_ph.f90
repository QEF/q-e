!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine punch_plot_ph
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output the change of the charge density,
  !     due to a perturbation ubare defined in the routine phq_setup. It
  !     can be read by chdens.f which cuts a bidimensional plane
  !     to plot contour levels,
  !     or selects a line for a usual line plot. The routine produces
  !     a file with the name in the variable fildrho# given in input.
  !
#include "machine.h"


  use pwcom
  use parameters, only : DP
  use phcom
#ifdef __PARA
  use para
#endif
  implicit none

  integer :: iunplot, ios, ipert, irr, na, ir, imode0, plot_num, jpol
  ! unit of the plot file
  ! integer variable for I/O contr
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on mesh points
  ! the starting mode
  ! compatibility variable
  ! counter on polarization

  character :: filin * 42
  ! complete name of the file

  real(kind=DP), allocatable :: raux (:)
  ! auxiliary vector

  complex(kind=DP) :: ps, ZDOTC
  complex(kind=DP), allocatable :: aux (:,:,:), aux1 (:,:)
  ! the scalar product
  ! scalar product function
  ! auxiliary space to rotate the
  ! induced charge
#ifdef __PARA
  ! auxiliary vector
  real(kind=DP), allocatable :: raux1 (:)
#endif

  if (fildrho.eq.' ') return
  write (6, '(/5x,"Calling punch_plot_ph" )')
  write (6, '(5x,"Writing on file  ",a)') fildrho
  !
  !    reads drho from the file
  !
  allocate (aux  (  nrxx,nspin,3))    
  allocate (aux1 (  nrxx,nspin))    
  allocate (raux (  nrxx))    
  !
  !
  !     reads the delta_rho on the aux variable
  !
  call setv (2 * nrxx * nspin, 0.d0, aux1, 1)
  imode0 = 0
  do irr = 1, nirr
     if (comp_irr (irr) .eq.1) then
        do ipert = 1, npert (irr)
           call davcio_drho (aux (1, 1, ipert), lrdrho, iudrho, imode0 + &
                ipert, - 1)
        enddo
#ifdef __PARA
        call psymdvscf (npert (irr), irr, aux)
#else
        call symdvscf (npert (irr), irr, aux)
#endif
        do ipert = 1, npert (irr)
           ps = ZDOTC (3 * nat, ubar, 1, u (1, imode0 + ipert), 1)
           call ZAXPY (nrxx * nspin, ps, aux (1, 1, ipert), 1, aux1, 1)
        enddo
     endif
     imode0 = imode0 + npert (irr)
  enddo
  !
  !     write on output the change of the charge
  !
  iunplot = 4
  filin = trim(fildrho)
#ifdef __PARA
  if (me.eq.1.and.mypool.eq.1) then
#endif
     open (unit = iunplot, file = filin, status = 'unknown', err = &
          100, iostat = ios)

100  call errore ('plotout', 'opening file'//filin, abs (ios) )
     rewind (iunplot)
     !
     !       Here we write some information quantity which are always necessa
     !
     plot_num = 0
     write (iunplot, '(a)') title
     write (iunplot, '(8i8)') nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, &
          ntyp
     write (iunplot, '(i6,6f12.8)') ibrav, celldm
     write (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, plot_num
     write (iunplot, 200) (na, atm (ityp (na) ), zv (ityp (na) ), &
          (tau (jpol, na), jpol = 1, 3), na = 1, nat)
200  format   (3x,i2,3x,a6,3x,f5.2,3x,3f14.10)
#ifdef __PARA
  endif
#endif
  !
  !      plot of the charge density
  !
  call DCOPY (nrxx, aux1 (1, 1), 2, raux, 1)

  if (lsda) call DAXPY (nrxx, 1.d0, aux1 (1, 2), 2, raux, 1)
#ifdef __PARA
  allocate (raux1( nrx1 * nrx2 * nrx3))    
  call gather (raux, raux1)
  if (me.eq.1.and.mypool.eq.1) write (iunplot, * ) (raux1 (ir), &
       ir = 1, nrx1 * nrx2 * nrx3)
  deallocate (raux1)
#else
  write (iunplot, * ) (raux (ir), ir = 1, nrxx)
#endif

  close (unit = iunplot)
  deallocate (raux)
  deallocate (aux1)
  deallocate (aux)
  return
end subroutine punch_plot_ph
