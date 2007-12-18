!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine pcgreen (avg_iter, thresh, ik, et_ )
  !-----------------------------------------------------------------------
  !
  ! Solve the linear system which defines the change of the wavefunctions
  ! due to the electric field for a given k_point in a non self-consistent
  ! way. The self-consistent variation of the potential has been computed
  ! previously and is in the common variable dvscfs
  !
#include "f_defs.h"
  use kinds, only : DP
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  use phcom
  implicit none

  !
  ! Input variables
  !
  integer :: ik
  ! input: k-point under consideration

  real(DP) :: avg_iter, thresh, et_ (nbnd)
  ! in/out: # of diagonalization iterations
  ! input: convergence threshold
  ! input: eigenvalues of the hamiltonian

  !
  ! Local variables
  !
  logical :: conv_root
  ! .true. if linter is converged

  integer :: ibnd, jbnd, ig, lter
  ! counters on bands
  ! counter on G-points
  ! # of diagonalization iterations

  real(DP) :: anorm
  ! the norm of the error

  real(DP) , allocatable :: h_diag(:,:), eprec(:)
  ! the diagonal part of the Hamiltonian
  ! cut-off for preconditioning

  complex(DP) :: ZDOTC
  ! the scalar product function

  complex(DP) , allocatable :: ps(:,:), auxg (:)
  ! auxiliary work space

  external   ch_psi_all, cg_psi

  allocate (h_diag ( npwx, nbnd ))
  allocate (auxg   ( npwx ))
  allocate (eprec  ( nbnd ))
  allocate (ps     ( nbnd, nbnd ))

  !
  ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
  !
  CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
       (1.d0,0.d0), evc(1,1), npwx, dvpsi(1,1), npwx, (0.d0,0.d0), &
       ps(1,1), nbnd )
#ifdef __PARA
  call reduce (2 * nbnd * nbnd_occ (ik), ps)
#endif
  !
  ! |dvspi> = - (|dvpsi> - S|evc><evc|dvpsi>)
  ! note the change of sign!
  !
  ! uncomment for ultrasoft PPs
  ! note that spsi is used as work space to store S|evc>
  ! CALL calbec ( npw, vkb, evc, becp, nbnd_occ(ik) )
  ! CALL s_psi (npwx, npw, nbnd_occ(ik), evc, spsi)
  ! CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
  !     (1.d0,0.d0), spsi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
  !     dvpsi(1,1), npwx )
  !
  ! comment  for ultrasoft PPs
  CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
       (1.d0,0.d0), evc(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
       dvpsi(1,1), npwx )
  !
  ! iterative solution of the linear system (H-e)*dpsi=dvpsi
  ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
  !
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        auxg (ig) = g2kin (ig) * evc(ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0 * ZDOTC (npw, evc(1, ibnd), 1, auxg, 1)
  enddo
#ifdef __PARA
  call reduce (nbnd_occ (ik), eprec)
#endif
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     enddo
  enddo
  conv_root = .true.

  call cgsolve_all( ch_psi_all, cg_psi, et_, dvpsi, dpsi, h_diag, &
                    npwx, npw, thresh, ik, lter, conv_root, anorm, &
                    nbnd_occ(ik), 1 )

  avg_iter = avg_iter + DBLE (lter)
  if (.not.conv_root) write(6, &
      "(5x,'kpoint',i4,' ibnd',i4, ' pcgreen: root not converged',e10.3)") &
      ik,ibnd,anorm

  deallocate (ps)
  deallocate (eprec)
  deallocate (auxg)
  deallocate (h_diag)

  return
end subroutine pcgreen
