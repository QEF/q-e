!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine punch_plot_e
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output the change of the charge density,
  !     due to an electric field in a real space mesh which can be read
  !     by chdens.f which cuts a bidimensional plane to plot contour level
  !     or selects a line for a usual line plot. The routine produces
  !     3 files with the change of charge density due to perturbations in
  !     three cartesian directions. The names of the files are
  !     in the variable fildrho given in input.
  !
#include "machine.h"

  USE io_global,  ONLY : stdout
  use pwcom
  use parameters, only : DP
  use phcom

#ifdef __PARA
  use para
#endif
  implicit none

  integer :: plot_num, iunplot, ios, ipol, jpol, na, ir, nt
  ! type of plot (not used)
  ! unit of the plot file
  ! integer variable for I/O contr
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on mesh points

  character :: caux * 1, filin * 80, which*2
  ! used to compose the name
  ! complete name of the file

  real(kind=DP), allocatable  :: raux (:)
  ! auxiliary vector

  complex(kind=DP), allocatable :: aux (:,:), aux1 (:,:)
  ! auxiliary space to rotate the
  ! induced charge

#ifdef __PARA
  ! auxiliary vector
  real(kind=DP), allocatable :: raux1 (:)
#endif

  if (fildrho.eq.' ') return
  WRITE( stdout, '(/5x,"Calling punch_plot_e" )')
  WRITE( stdout, '(5x,"Writing on file  ",a)') fildrho
  !
  !    reads drho from the file
  !
  allocate (aux  (  nrxx,3))    
  allocate (aux1 (  nrxx,3))    
  allocate (raux (  nrxx))    
  !
  !     reads the delta_rho on the aux variable
  !
  do ipol = 1, 3
     call davcio_drho (aux (1, ipol), lrdrho, iudrho, ipol, - 1)
  enddo
  !
  !     symmetrize
  !
#ifdef __PARA
     call psyme (aux)
#else
     call syme (aux)
#endif
  !
  !     rotate the charge and transform to cartesian coordinates
  !
  call setv (6 * nrxx, 0.0d0, aux1, 1)
  do ipol = 1, 3
     do jpol = 1, 3
        call DAXPY (2 * nrxx, bg (ipol, jpol), aux (1, jpol), 1, aux1 (1, &
             ipol), 1)
     enddo
  enddo
  !
  !     write on output the change of the charge
  !
  iunplot = 4
  which='_e'
  do ipol = 1, 3
     write (caux, '(i1)') ipol
     filin = trim(fildrho) //which//caux
#ifdef __PARA
     if (me.eq.1.and.mypool.eq.1) then
#endif
        open (unit = iunplot, file = filin, status = 'unknown', err = &
             100, iostat = ios)

100     call errore ('plotout', 'opening file'//filin, abs (ios) )
        rewind (iunplot)
        !
        !       Here we write some information quantity which are always necessa
        !
        ! not used
        plot_num = - 1
        write (iunplot, '(a)') title
        write (iunplot, '(8i8)') nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, &
             ntyp
        write (iunplot, '(i6,6f12.8)') ibrav, celldm
        write (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, plot_num
        write (iunplot, '(i4,3x,a2,3x,f5.2)') &
                                (nt, atm (nt), zv (nt), nt=1, ntyp)
        write (iunplot, '(i4,3x,3f14.10,3x,i2)') (na, &
          (tau (jpol, na), jpol = 1, 3), ityp (na), na = 1, nat)

#ifdef __PARA
     endif
#endif
     !
     !      plot of the charge density
     !

     call DCOPY (nrxx, aux1 (1, ipol), 2, raux, 1)
#ifdef __PARA
     allocate (raux1( nrx1 * nrx2 * nrx3))    
     call gather (raux, raux1)
     if (me.eq.1.and.mypool.eq.1) write (iunplot, '(5(1pe17.9))') &
          (raux1 (ir) , ir = 1, nrx1 * nrx2 * nrx3)
     deallocate (raux1)
#else
     write (iunplot, '( 5( 1pe17.9 ) )') (raux (ir) , ir = 1, nrxx)
#endif
     close (unit = iunplot)
  enddo
  deallocate (raux)
  deallocate (aux1)
  deallocate (aux)
  return

end subroutine punch_plot_e
