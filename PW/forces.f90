!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine forces
  !----------------------------------------------------------------------
  !
  !    This routine is a driver routine which computes the forces
  !    acting on the atoms. The complete expression of the forces
  !    contains four parts which are computed by different routines:
  !    a)   force_lc,   local contribution to the forces
  !    b)   force_cc,   contribution due to NLCC
  !    c)   force_ew,   contribution due to the electrostatic ewald term
  !    d)   force_us ,  contribution due to the non-local potential
  !    e)   force_corr  correction term for incomplete self-consistency
  !
#include "machine.h"

  use pwcom
  implicit none

  real(kind=DP), allocatable :: forcenl (:,:), forcelc (:,:), forcecc (:,:), &
       forceion (:,:), forcescc (:,:), forceh(:,:)
  ! nonlocal, local, core-correction, ewald, and scf correction terms
  real(kind=DP) :: sum, sumscf

  integer :: ipol, na
  ! counter on polarization
  ! counter on atoms

  call start_clock ('forces')

  allocate ( forcenl(3,nat), forcelc(3,nat), forcecc(3,nat), forceh(3,nat), &
             forceion(3,nat), forcescc(3,nat) )
  forcescc(:,:) = 0.d0
  forceh(:,:) = 0.d0

  write (6, '(/,5x,"Forces acting on atoms (Ry/au):", / )')
  !
  !    The  nonlocal contribution is computed here
  !
  call force_us (forcenl)
  !
  !    The local contribution
  !
  call force_lc (nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
       nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, g, rho, nl, nspin, &
       gstart, gamma_only, vloc, forcelc)
  !
  !    The NLCC contribution
  !
  call force_cc (forcecc)
  !
  !    The Hubbard contribution
  !
  if (lda_plus_u) call force_hub ( forceh )
  !
  !    The ionic contribution is computed here
  !
  call force_ew (alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
       gg, ngm, gstart, gamma_only, gcutm, strf, forceion)
  !
  !    The SCF contribution
  !
  call force_corr (forcescc)
  !
  !  here we sum all the contributions and compute the total force acting
  !  on the crstal
  !
  do ipol = 1, 3
     sum = 0.d0
     do na = 1, nat
        force(ipol,na) = forcenl (ipol, na) + &
                         forceion(ipol, na) + &
                         forcelc (ipol, na) + &
                         forcecc (ipol, na) + &
                         forceh (ipol, na) + &
                         forcescc(ipol, na)
        if (tefield) force(ipol,na)=force(ipol,na)+forcefield(ipol,na)
        sum = sum + force (ipol, na)
     enddo
     !         write(6,*) 'sum = ', sum
     !
     ! impose total force = 0
     !
     do na = 1, nat
        force (ipol, na) = force (ipol, na) - sum / nat
     enddo
  enddo
  !
  ! resymmetrize (should not be needed, but...)
  !
  if (nsym.gt.1) then
     do na = 1, nat
        call trnvect (force (1, na), at, bg, - 1)
     enddo
     call symvect (nat, force, nsym, s, irt)
     do na = 1, nat
        call trnvect (force (1, na), at, bg, 1)
     enddo
  endif
  !
  !   write on output the forces
  !
  do na = 1, nat
     write (6, 9035) na, ityp (na), (force (ipol, na), ipol = 1, 3)
  enddo
#ifdef DEBUG
  write (6, '(5x,"The SCF correction term to forces")')
  do na = 1, nat
     write (6, 9035) na, ityp (na), (forcescc (ipol, na), ipol = 1, 3)
  enddo
  write (6, '(5x,"The Hubbard contribution to forces")')
  do na = 1, nat
     write (6, 9035) na, ityp (na), (forceh(ipol, na), ipol = 1, 3)
  enddo
#endif
  sum = 0.d0
  sumscf = 0.d0
  do ipol = 1, 3
     do na = 1, nat
        sum = sum + abs (force (ipol, na) )
        sumscf = sumscf + abs (forcescc (ipol, na) )
     enddo
  enddo

  write (6, '(/5x,"Total force = ",f12.6,5x, &
       &                "Total SCF correction = ",f12.6)') sum, sumscf
#ifdef __PARA

  call check (3 * nat, force)
#endif

  deallocate (forcenl, forcelc, forcecc, forceh, forceion, forcescc)

  lforce = .true.
  call stop_clock ('forces')

  return
9035 format (5x,'atom ',i3,' type ',i2,'   force = ',3f14.8)
end subroutine forces

