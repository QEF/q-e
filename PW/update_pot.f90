!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine update_pot
  !-----------------------------------------------------------------------
  !
  !     update potential, use the integer variable order to decide the way
  !
  !     order = 0       copy the old potential (nothing is done)
  !
  !     order = 1       subtract old atomic charge density and sum the new
  !                     if dynamics is done the routine extrapolates also
  !                     the difference between the the scf charge and the
  !                     atomic one,
  !
  !     order = 2       extrapolate the wavefunctions:
  !                       |psi(t+dt)> = 2*|psi(t)> - |psi(t-dt)>
  !
  !     order = 3       extrapolate the wavefunctions with the second-order
  !                     formula:
  !                       |psi(t+dt)> = |psi(t) +
  !                                   + alpha0*(|psi(t)> -    |psi(t-dt)>
  !                                   + beta0* (|psi(t-dt)> - |psi(t-2*dt)>
  !
  !                     where alpha0 and beta0 are calculated in "dynamics" so
  !                     that |tau'-tau(t+dt)| is minimum; tau' and tau(t+dt)
  !                     are respectively the atomic positions at time t+dt
  !                     and  the extrapolated one:
  !                       tau(t+dt) = tau(t) +
  !                                    + alpha0*( tau(t) - tau(t-dt) )
  !                                    + beta0*( tau(t-dt) -tau(t-2*dt) )
  !
  !

  use pwcom
  implicit none

  call start_clock ('update_pot')
  if (order.eq.0) return
  if (order.gt.2.and.iswitch.le.2) then
     order = 2
     call errore ('update_pot', 'order > 2 not allowed in bfgs', - 1)

  endif
  call extrapolate_charge

  if (order.ge.2) call extrapolate_wfcs

  call stop_clock ('update_pot')
  return

end subroutine update_pot
!-----------------------------------------------------------------------
subroutine extrapolate_charge
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  !

  use pwcom
  use io_files, only: prefix
  implicit none
  integer :: ir
  ! do-loop variable on FFT grid

  real(kind=DP), allocatable :: work (:), work1 (:)
  ! work is the difference between charge density and atomic charge at time t
  ! work1 is the same thing at time t-dt
  real(kind=DP) :: charge

  allocate(work(nrxx))
  work(:) = 0.d0
  !
  !     if order = 1 update the potential subtracting to the charge density
  !     the "old" atomic charge and summing the new one
  !
  write (6,'(/5x,"NEW-OLD atomic charge density approx. for the potential")')
  !
  ! in the lsda case the magnetization will follow rigidly the density kee
  ! fixed the value of zeta=mag/rho_tot. zeta is set here and put in rho(*
  ! while rho(*,1) will contain the total valence charge
  !
  if (lsda) call rho2zeta (rho, rho_core, nrxx, nspin, + 1)
  !
  !     subtract the old atomic charge density
  !
  call atomic_rho (work, 1)

  call DAXPY (nrxx, - 1.0d0, work, 1, rho, 1)
  if (lmovecell) call DSCAL (nrxx, omega_old, rho, 1)
  !
  !     if dynamics extrapolate the difference between the atomic charge a
  !     the self-consistent one
  !
  if (iswitch.gt.2) then
     if (istep.eq.1) then
        call io_pot ( + 1,trim(prefix)//'.oldrho', rho, 1)
     else
        allocate(work1(nrxx))
        work1(:) = 0.d0
        call io_pot ( - 1,trim(prefix)//'.oldrho', work, 1)
        call io_pot ( + 1,trim(prefix)//'.oldrho', rho, 1)
        if (istep.eq.2) then
           call io_pot ( + 1,trim(prefix)//'.oldrho2', work, 1)
        endif
        call io_pot ( - 1,trim(prefix)//'.oldrho2', work1, 1)
        call io_pot ( + 1,trim(prefix)//'.oldrho2', work, 1)
        !
        ! alpha0 and beta0 have been calculated in dynamics or in vcsmd subs.
        !
        do ir = 1, nrxx
           rho(ir,1) = rho(ir,1) + alpha0 * ( rho(ir,1) - work(ir) ) + &
                                    beta0 * ( work(ir) - work1(ir) )
        enddo
        deallocate(work1)
     endif

  endif
  if (lmovecell) call DSCAL (nrxx, 1.0d0 / omega, rho, 1)
  !
  !     calculate structure factors for the new positions
  !
  if (lmovecell) call scale_h
  call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)
  !
  !     add atomic charges in the new positions
  !
  call atomic_rho (work, 1)
  call DAXPY (nrxx, 1.0d0, work, 1, rho, 1)
  call set_rhoc
  !
  ! reset up and down charge densities in the LSDA case
  !

  if (lsda) call rho2zeta (rho, rho_core, nrxx, nspin, -1)

  call v_of_rho (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
       ehart, etxc, vtxc, charge, vr)
  !
  !     write potential (and rho) on file
  !

  if (imix.ge.0) call io_pot(+1,trim(prefix)//'.rho',rho,nspin)
  call io_pot(+1,trim(prefix)//'.pot',vr,nspin)

  deallocate(work)

  return


end subroutine extrapolate_charge
!-----------------------------------------------------------------------
subroutine extrapolate_wfcs
  !-----------------------------------------------------------------------
  !
  !     This routine extrapolate the wfc's after a "parallel alignment"
  !     of the basis of the t-dt and t time steps, according to the Mead
  !     recipe, see Rev. Mod. Phys., vol 64, pag. 51 (1992), eqs. 3.20-3.2
  !
  !
#include "machine.h"
  use pwcom
  USE wavefunctions,    ONLY : evc
  implicit none
#define ONE (1.d0,0.d0)
#define ZERO (0.d0,0.d0)

  integer :: j, i, ik
  ! do-loop variables
  ! counter on k-points

  complex(kind=DP), allocatable:: u_m (:,:), s_m (:,:), sp_m (:,:), temp (:,:)
  ! the unitary matrix (eq. 3.21)
  ! the overlap matrix s (eq. 3.24)
  ! its dagger
  ! workspace
  complex(kind=DP), allocatable:: evcold(:,:)
  ! wavefunctions at previous iteration

  real(kind=DP), allocatable :: ew (:)
  ! the eigenvalues of sp_m*s_m

  logical :: first
  ! Used for initialization
  data first / .true. /


  save first
  if (first) then
     first = .false.
     if (isolve.eq.1.and.startingwfc.eq.'atomic') then
        deallocate(evc)
        allocate(evc(npwx,nbnd))
     endif
  endif
  allocate(evcold(npwx,nbnd))
  if (istep.eq.1) then
     if (nks.gt.1) rewind (iunigk)
     do ik = 1, nks
        if (nks.gt.1) read (iunigk) npw, igk
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call ZCOPY (npwx * nbnd, evc, 1, evcold, 1)
        call davcio (evcold, nwordwfc, iunoldwfc, ik, 1)
     enddo
  else
     if (order.eq.2) then
        write (6, '(5x,"Extrapolating wave-functions (first order) ...")')
     else
        write (6, '(5x,"Extrapolating wave-functions (second order) ...")')
     endif

     allocate ( u_m(nbnd,nbnd), s_m(nbnd,nbnd), sp_m(nbnd,nbnd), &
                temp(nbnd,nbnd), ew(nbnd) )

     if (nks.gt.1) rewind (iunigk)
     do ik = 1, nks
        if (nks.gt.1) read (iunigk) npw, igk
        call davcio (evcold, nwordwfc, iunoldwfc, ik, - 1)
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        if (istep.eq.2.and.order.gt.2) then
           call davcio (evcold, nwordwfc, iunoldwfc2, ik, 1)
        endif
        !
        !     construct s_m = <evcold|evc>
        !
        call ZGEMM ('c', 'n', nbnd, nbnd, npw, ONE, evcold, npwx, evc, &
             npwx, ZERO, s_m, nbnd)
#ifdef __PARA
        call reduce (2 * nbnd * nbnd, s_m)
#endif
        !
        !     temp = sp_m*s_m
        !
        call ZGEMM ('c', 'n', nbnd, nbnd, nbnd, ONE, s_m, nbnd, s_m, &
             nbnd, ZERO, temp, nbnd)
        !
        !     diagonalize temp, use u_m as workspace to accomodate the eigenvect
        !     matrix which diagonalizes temp, sp_m is its hermitean conjugate
        !
        call cdiagh (nbnd, temp, nbnd, ew, u_m)
        do i = 1, nbnd
           do j = 1, nbnd
              sp_m (j, i) = conjg (u_m (i, j) ) / sqrt (ew (j) )
           enddo
        enddo
        call ZGEMM ('n', 'n', nbnd, nbnd, nbnd, ONE, u_m, nbnd, sp_m, &
             nbnd, ZERO, temp, nbnd)
        !
        !     temp = [ sp_m * s_m ]^(-1/2)
        !
        call ZGEMM ('n', 'c', nbnd, nbnd, nbnd, ONE, temp, nbnd, s_m, &
             nbnd, ZERO, u_m, nbnd)
        !
        !     and u_m is the unitary matrix [ sp_m * s_m ]^(-1/2)*sp_m (eq.3.29)
        !     now use evcold as workspace to calculate
        !
        !                        evcold_i = sum_j evc_j*u_m_ji
        !
        call ZGEMM ('n', 'n', npw, nbnd, nbnd, ONE, evc, npwx, u_m, &
             nbnd, ZERO, evcold, npwx)
        !
        !     and copy evcold in evc
        !
        call ZCOPY (npwx * nbnd, evcold, 1, evc, 1)
        !
        !     save on file evc
        !
        call davcio (evc, nwordwfc, iunwfc, ik, 1)
        !
        !     re-read from file the right evcold
        !
        call davcio (evcold, nwordwfc, iunoldwfc, ik, - 1)
        !
        !     extrapolate the wfc's, if order=3 use the second order extrapolati
        !     formula, alpha0 and beta0 are calculated in "dynamics"
        !
        if (order.gt.2) then
           do j = 1, nbnd
              do i = 1, npw
                 evc (i, j) = (1 + alpha0) * evc (i, j) + (beta0 - alpha0) &
                      * evcold (i, j)
              enddo
           enddo
           call davcio (evcold, nwordwfc, iunoldwfc2, ik, - 1)
           do j = 1, nbnd
              do i = 1, npw
                 evc (i, j) = evc (i, j) - beta0 * evcold (i, j)
              enddo
           enddo
        else
           do j = 1, nbnd
              do i = 1, npw
                 evc (i, j) = 2 * evc (i, j) - evcold (i, j)
              enddo
           enddo
        endif
        !
        !     move the files: "old" -> "old1" and "now" -> "old"
        !
        if (order.gt.2) then
           call davcio (evcold, nwordwfc, iunoldwfc, ik, - 1)
           call davcio (evcold, nwordwfc, iunoldwfc2, ik, 1)
        endif
        call davcio (evcold, nwordwfc, iunwfc, ik, - 1)
        call davcio (evcold, nwordwfc, iunoldwfc, ik, 1)
        !
        !     save evc on file iunwfc
        !
        call davcio (evc, nwordwfc, iunwfc, ik, 1)
     enddo
     deallocate(u_m, s_m, sp_m, temp, ew)
  endif
  deallocate (evcold)
  return
end subroutine extrapolate_wfcs

