!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, w2)
  !-----------------------------------------------------------------------
  !
  !   This routine diagonalizes the dynamical matrix and returns
  !   displacement patterns in "dyn". The frequencies are written
  !   on output from this routine.
  !
  !
  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  USE constants, ONLY : amu_ry, RY_TO_THZ, RY_TO_CMM1
  USE io_dyn_mat, ONLY : write_dyn_mat_tail
  USE control_ph, ONLY : xmldyn
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

  real(DP) :: xq (3), amass (ntyp), w2 (3 * nat)
  ! input: q vector
  ! input: the masses
  ! output: the frequencies squared
  complex(DP) :: dyn (3 * nat, nmodes)
  ! input: the dynamical matrix
  !
  !   here the local variables
  !
  integer :: nta, ntb, nu_i, nu_j, mu, na, nb, i
  ! counters

  real(DP) :: w1, unorm
  ! the frequency
  ! norm of u

  complex(DP) :: z (3 * nat, 3 * nat)
  ! the eigenvectors

  !
  ! fill the second half of the matrix (imposing hermiticity !)
  !
  do nu_i = 1, nmodes
     do nu_j = 1, nu_i
        dyn (nu_i, nu_j) = 0.5d0 * (dyn (nu_i, nu_j) + &
                             CONJG(dyn (nu_j, nu_i) ) )
        dyn (nu_j, nu_i) = CONJG(dyn (nu_i, nu_j) )
     enddo
  enddo
  !
  ! divide the dynamical matrix by the masses (beware: amass is in amu)
  !
  do nu_i = 1, nmodes
     na = (nu_i - 1) / 3 + 1
     nta = ityp (na)
     do nu_j = 1, nmodes
        nb = (nu_j - 1) / 3 + 1
        ntb = ityp (nb)
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) / sqrt (amass (nta)*amass (ntb)) &
                                            / amu_ry
     enddo
  enddo
  !
  ! solve the eigenvalue problem
  !
  call cdiagh (nmodes, dyn, 3 * nat, w2, z)
  !
  !    Writes on output the displacements and the normalized frequencies.
  !
  WRITE( stdout, 9000) (xq (i), i = 1, 3)
  if (iudyn /= 0) write (iudyn, 9000) (xq (i), i = 1, 3)

9000 format(/,5x,'Diagonalizing the dynamical matrix', &
       &       //,5x,'q = ( ',3f14.9,' ) ',//,1x,74('*'))
  dyn (:,:) = (0.d0, 0.d0)
  do nu_i = 1, nmodes
     w1 = sqrt (abs (w2 (nu_i) ) )
     if (w2 (nu_i) < 0.d0) w1 = - w1
     WRITE( stdout, 9010) nu_i, w1 * RY_TO_THZ, w1 * RY_TO_CMM1
     if (iudyn /= 0) write (iudyn, 9010) nu_i, w1 * RY_TO_THZ, w1 * RY_TO_CMM1
9010 format   (5x,'freq (',i5,') =',f15.6,' [THz] =',f15.6,' [cm-1]')
     !
     ! write displacements onto matrix dyn
     !
     unorm = 0.d0
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu_i) = z (mu, nu_i) / sqrt (amu_ry * amass (ityp (na) ) )
        unorm = unorm + dyn (mu, nu_i) * CONJG(dyn (mu, nu_i) )
     enddo
     if (iudyn /= 0) then
         write (iudyn, '(" (",6f10.6," ) ")') &
          (dyn (mu, nu_i) / sqrt (unorm) , mu = 1, 3 * nat)
     else
         z(:,nu_i)=dyn (:, nu_i) / sqrt (unorm)
     endif
  enddo
  WRITE( stdout, '(1x,74("*"))')
  if (iudyn /= 0) write (iudyn, '(1x,74("*"))')
  IF (xmldyn) CALL write_dyn_mat_tail(nat, w2, z)
  return

end subroutine dyndia
