!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine potinit
  !-----------------------------------------------------------------------
  !
  !     This routine initializes the self consistent potential in the array
  !     vr. There are three possible cases:
  !
  !     a) In this run the code is restarting from a broken run
  !     b) The potential (or rho) is read from file
  !     c) if a and b are both false, the total charge is computed
  !        as a sum of atomic charges, and the corresponding potential
  !        is saved in vr
  !
#include "machine.h"
  use pwcom
  use io, only: prefix
#ifdef __PARA
  use para
  use mp
#endif
  implicit none

  real(kind=DP) :: charge
  ! the starting charge
  integer :: ios, ionode_id=0 
  integer :: ldim
  ! integer variable for I/O control
  logical :: exst
  !
#ifdef __PARA
  if (me.eq.1.and.mypool.eq.1) then
#endif
     if (imix.ge.0) then
        call seqopn (4, trim(prefix)//'.rho', 'unformatted', exst)
     else
        call seqopn (4, trim(prefix)//'.pot', 'unformatted', exst)
     end if
     if (exst) then
        close (unit =4, status = 'keep')
     else
        close (unit =4, status = 'delete')
     end if
#ifdef __PARA
  endif
  call mp_bcast( exst, ionode_id )
#endif
  if (startingpot=='file' .and. .not.exst) then
     write (6, '(5x,"Cannot read pot/rho file: not found")')
     startingpot='atomic'
  end if

  ! First case, the potential is read from file
  ! NB: this case applies also for a restarting run, in which case
  !     potential and rho files have been read from the restart file
  !
  if (startingpot=='file') then
     if (imix.ge.0) then
        call io_pot ( -1, trim(prefix)//'.rho', rho, nspin)
        write (6, '(/5x,"The initial density is read from file ", &
                   &    a14)') trim(prefix)//'.rho'
        !
        ! here we compute the potential which correspond to the initial charge
        !
        call v_of_rho (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
             nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
             ehart, etxc, vtxc, charge, vr)
        !
        if (abs (charge-nelec)  / charge.gt.1.0d-4) &
             write (6, '(/5x,"starting charge =",f10.5)') charge
     else
        call io_pot ( - 1,  trim(prefix)//'.pot', vr, nspin)
        write (6, '(/5x,"The initial potential is read from file ", &
                   &     a14)') trim(prefix)//'.pot'
     end if
     !
     ! The occupations ns also need to be read in order to build up the poten
     !
     if (lda_plus_u) then  
        ldim = 2 * Hubbard_lmax + 1
#ifdef __PARA
        if (me.eq.1.and.mypool.eq.1) then
#endif
           call seqopn (iunocc, trim(prefix)//'.occup', 'formatted', exst)
           read (iunocc, * ) ns
           close (unit = iunocc, status = 'keep')
#ifdef __PARA
        else  
           call setv (nat * nspin * ldim * ldim, 0.d0, ns, 1)  
        endif
        call reduce (nat * nspin * ldim * ldim, ns)  
        call poolreduce (nat * nspin * ldim * ldim, ns)  
#endif
        call DCOPY(nat*nspin*ldim*ldim,ns,1,nsnew,1)
     endif
  else
     !
     ! Second case, the potential is built from a superposition of atomic
     ! charges contained in the array rho_at and already set in readin-readva
     !
     write (6, '(/5x,"Initial potential from superposition", &
          &                " of free atoms")')
     !
     ! in the lda+U case set the initial value of ns
     !

     if (lda_plus_u) then
        ldim = 2 * Hubbard_lmax + 1
        call init_ns  
        call DCOPY(nat*nspin*ldim*ldim,ns,1,nsnew,1)
     end if

     call atomic_rho (rho, nspin)
     if (input_drho.ne.' ') then
        if (lsda) call errore ('potinit', ' lsda not allowed in drho', 1)
        call io_pot ( - 1, input_drho, vr, nspin)
        write (6, '(/5x,"a scf correction to at. rho is read from", &
             &          a14)') input_drho
        call DAXPY (nrxx, 1.d0, vr, 1, rho, 1)
     endif
     !
     ! here we compute the potential which correspond to the initial charge
     !
     call v_of_rho (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
          ehart, etxc, vtxc, charge, vr)
     !
     if (abs (charge-nelec)  / charge.gt.1.0d-4) &
          write (6, '(/5x,"starting charge =",f10.5)') charge

  endif
  !
  ! define the total local potential (external+scf)
  !

  call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)
  !
  ! write on output the parameters used in the lda+U calculation
  !
  if (lda_plus_u) then
     write (6, '(/5x,"Parameters of the lda+U calculation:")')
     write (6, '(5x,"Number of iteration with fixed ns =",i3)') &
          niter_with_fixed_ns
     write (6, '(5x,"Starting ns and Hubbard U :")')
     call write_ns

  endif

  if (imix.ge.0) call io_pot ( +1,  trim(prefix)//'.rho', rho, nspin)
  call io_pot ( +1,  trim(prefix)//'.pot', vr, nspin)

  return
20 call errore ('potinit', 'error reading '//trim(prefix)//'.pot', abs(ios) )
end subroutine potinit

