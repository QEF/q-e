!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, &
     w2)
  !-----------------------------------------------------------------------
  !
  !   This routine diagonalizes the dynamical matrix and returns
  !   displacement patterns in "dyn". The frequencies are written
  !   on output from this routine.
  !
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  !   first the dummy variables
  !

  integer :: nmodes, nat, ntyp, ityp (nat), iudyn
  ! input: the total number of modes
  ! input: the number of atoms
  ! input: the number of types
  ! input: the types of atoms
  ! input: the unit with the dynamical matrix

  real(kind=DP) :: xq (3), amass (ntyp), w2 (3 * nat)
  ! input: q vector
  ! input: the masses
  ! output: the frequencies square
  complex(kind=DP) :: dyn (3 * nat, nmodes)
  ! input: the dynamical matrix
  !
  !   here the local variables
  !


  integer :: nta, ntb, nu_i, nu_j, mu, na, nb, i
  ! counter on types
  ! counter on types
  ! counter on modes
  ! counter on modes
  ! counter on 3*nat
  ! counter on atoms
  ! counter on atoms
  ! counter on modes

  real(kind=DP) :: rydthz, rydcm1, w1, unorm
  ! conversion from a.u. to terahertz
  ! conversion from a.u. to cm^-1
  ! the frequency
  ! norm of u

  complex(kind=DP) :: z (3 * nat, 3 * nat)
  ! the eigenvectors

  !
  ! fill the second half of the matrix (imposing hermiticity !)
  !
  do nu_i = 1, nmodes
     do nu_j = 1, nu_i
        dyn (nu_i, nu_j) = 0.5d0 * (dyn (nu_i, nu_j) + conjg (dyn (nu_j, &
             nu_i) ) )
        dyn (nu_j, nu_i) = conjg (dyn (nu_i, nu_j) )
     enddo
  enddo
  !
  ! divide the dynamical matrix by the masses
  !
  do nu_i = 1, nmodes
     na = (nu_i - 1) / 3 + 1
     nta = ityp (na)
     do nu_j = 1, nmodes
        nb = (nu_j - 1) / 3 + 1
        ntb = ityp (nb)
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) / sqrt (amass (nta) * amass ( &
             ntb) )
     enddo
  enddo
  !
  ! solve the eigenvalue problem
  !
  call cdiagh (nmodes, dyn, 3 * nat, w2, z)
  !
  !     conversion factors  ryd=>thz e ryd=>1/cm
  !
  rydthz = 13.6058d0 * 241.796d0
  rydcm1 = 13.6058d0 * 8065.5d0
  !
  !    Writes on output the displacements and the normalized frequencies.
  !
  WRITE( stdout, 9000) (xq (i), i = 1, 3)
  if (iudyn.ne.0) write (iudyn, 9000) (xq (i), i = 1, 3)

9000 format(/,5x,'Diagonalizing the dynamical matrix', &
       &       //,5x,'q = ( ',3f14.9,' ) ',//,1x,74('*'))
  call setv (2 * 3 * nat * nmodes, 0.d0, dyn, 1)
  do nu_i = 1, nmodes
     w1 = sqrt (abs (w2 (nu_i) ) )
     if (w2 (nu_i) .lt.0.d0) w1 = - w1
     WRITE( stdout, 9010) nu_i, w1 * rydthz, w1 * rydcm1
     if (iudyn.ne.0) write (iudyn, 9010) nu_i, w1 * rydthz, w1 * &
          rydcm1
9010 format   (5x,'omega(',i2,') =',f15.6,' [THz] =',f15.6,' [cm-1]')
     !
     ! write displacements onto matrix dyn
     !
     unorm = 0.d0
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu_i) = z (mu, nu_i) / sqrt (amass (ityp (na) ) )
        unorm = unorm + dyn (mu, nu_i) * conjg (dyn (mu, nu_i) )
     enddo
     if (iudyn.ne.0) write (iudyn, '(" (",6f10.6," ) ")') (dyn (mu, &
          nu_i) / sqrt (unorm) , mu = 1, 3 * nat)
  enddo
  WRITE( stdout, '(1x,74("*"))')
  if (iudyn.ne.0) write (iudyn, '(1x,74("*"))')
  return

end subroutine dyndia
