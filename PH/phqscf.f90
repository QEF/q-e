!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine phqscf
  !-----------------------------------------------------------------------
  !
  !     This subroutine is the main driver of the self consistent cycle
  !     which gives as output the change of the wavefunctions and the
  !     change of the self-consistent potential due to a phonon of
  !     a fixed q or to an electric field.
  !
#include "machine.h"

  USE io_global,  ONLY : stdout
  use pwcom
  use parameters, only : DP
  use phcom
  implicit none

  integer :: irr, irr1, irrc, imode0
  ! counter on the representations
  ! counter on the representations
  ! number of representation computed
  ! counter on the modes

  real(kind=DP) :: tcpu, get_clock
  ! timing variables

  logical :: exst
  ! used to test the recover file

  external get_clock
  ! the change of density due to perturbations

  complex(kind=DP), allocatable :: drhoscf (:,:,:)

  call start_clock ('phqscf')
  !
  !    For each irreducible representation we compute the change
  !    of the wavefunctions
  !
  irrc = 0
  if (irr0 < 0) then
     irr0 = 0
     iter0= 0
  endif

 allocate (drhoscf( nrxx , nspin, npertx))    
  do irr = 1, nirr
     if ( (comp_irr (irr) == 1) .and. (done_irr (irr) == 0) ) then
        irrc = irrc + 1
        imode0 = 0
        do irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        enddo
        if (npert (irr) == 1) then
           WRITE( stdout, '(//,5x,"Representation #", i3," mode # ",i3)') &
                              irr, imode0 + 1
        else
           WRITE( stdout, '(//,5x,"Representation #", i3," modes # ",4i3)') &
                              irr, (imode0+irr1, irr1=1,npert(irr))
        endif
        !
        !    then for this irreducible representation we solve the linear system
        !
        call solve_linter (irr, imode0, npert (irr), drhoscf)
        !
        !   Add the contribution of this mode to the dynamical matrix
        !
        if (convt) then
           call drhodv (imode0, npert (irr), drhoscf)
           !
           !   add the contribution of the modes imode0+1 -> imode+npe
           !   to the effective charges Z(Us,E) (Us=scf,E=bare)
           !
           if (zue) call add_zstar_ue (imode0, npert (irr) )
           if (zue.and. okvan) call add_zstar_ue_us(imode0, npert (irr) )

        endif
        if (convt) then
           WRITE( stdout, '(/,5x,"Convergence has been achieved ")')
           done_irr (irr) = 1
           iter0 = 0
        else
           WRITE( stdout, '(/,5x,"No convergence has been achieved ")')
           call stop_ph (.false.)
        endif
        !
        !   We test here if we have done the appropriate number of
        !   representation
        !
        call seqopn (iunrec, 'recover', 'unformatted', exst)
        if (okvan) write (iunrec) int1, int2

        write(iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, zstarue0

        write(iunrec) irr, 0, convt, done_irr, comp_irr, ifat

        close (unit = iunrec, status = 'keep')
        tcpu = get_clock ('PHONON')
        if (tcpu > 1000000000) then
           !
           !            if (tcpu.gt.time_max) then
           ! temporary fix: recover does not work if program stop here
           !
           WRITE( stdout, '(/,5x,"Stopping for time limit ",2f10.0)') &
                tcpu, time_max
           call stop_ph (.false.)

        endif

        if (irrc >= maxirr) then
           WRITE( stdout, '(/,5x,"Stopping at Representation #",i6)') irr
#ifdef DEBUG
           if (me /= 1 .or. mypool /= 1) close (6)
#endif
           call stop_ph (.false.)
        endif
     endif

  enddo
  deallocate (drhoscf)

  call stop_clock ('phqscf')
  return
end subroutine phqscf
