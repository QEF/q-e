!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine pcgreen (avg_iter, thresh, ik, et_, auxg, spsi )
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

  real(kind=DP) :: avg_iter, thresh, et_ (nbnd)
  ! in/out: # of diagonalization iterations
  ! input: convergence threshold
  ! input: eigenvalues of the hamiltonian

  complex(kind=DP) :: auxg (npwx), spsi (npwx)
  ! auxiliary space

  !
  ! Local variables
  !
  logical :: conv_root
  ! .true. if linter is converged

  integer :: ibnd, jbnd, ig, lter
  ! counters on bands
  ! counter on G-points
  ! # of diagonalization iterations

  real(kind=DP) :: anorm
  ! the norm of the error

  real(kind=DP) , allocatable :: h_diag(:,:), eprec(:)
  ! the diagonal part of the Hamiltonian
  ! cut-off for preconditioning

  complex(kind=DP) :: ZDOTC
  ! the scalar product function

  complex(kind=DP) , allocatable :: ps(:)
  ! the scalar product

  external   ch_psi_all, cg_psi

  allocate (h_diag ( npwx, nbnd ))
  allocate (eprec  ( nbnd ))
  allocate (ps     ( nbnd ))

  !
  ! Ortogonalize dvpsi 
  !
  do ibnd = 1, nbnd_occ (ik)
     auxg (:) = (0.d0, 0.d0)
     do jbnd = 1, nbnd_occ (ik)
        ps (jbnd) = -ZDOTC (npw, evc (1, jbnd), 1, dvpsi (1, ibnd), 1)
     enddo
#ifdef __PARA
     call reduce (2 * nbnd, ps) 
#endif
     do jbnd = 1, nbnd_occ (ik)
        call ZAXPY (npw, ps (jbnd), evc (1, jbnd), 1, auxg, 1)
     enddo
! If you uncomment the following three lines you should comment the fourth
!     call ccalbec (nkb, npwx, npw, 1, becp, vkb, auxg)
!     call s_psi (npwx, npw, 1, auxg, spsi)
!     call DAXPY ( 2 * npw, 1.0d0, spsi, 1, dvpsi (1, ibnd), 1)
     call DAXPY (2 * npw, 1.0d0, auxg, 1, dvpsi (1, ibnd), 1)

  enddo
  !
  !    Here we change the sign of the known term
  !
  call DSCAL ( 2 * npwx * nbnd, -1.d0, dvpsi, 1)
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
                    nbnd_occ(ik) )

  avg_iter = avg_iter + dfloat (lter)
  if (.not.conv_root) write(6, &
      "(5x,'kpoint',i4,' ibnd',i4, ' pcgreen: root not converged',e10.3)") &
      ik,ibnd,anorm

  deallocate (h_diag)
  deallocate (eprec)
  deallocate (ps)

  return
end subroutine pcgreen
