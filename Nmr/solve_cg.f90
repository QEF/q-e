!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*******************************************************************

subroutine solve_cg(dev,ik, dpsi)

!*******************************************************************

#include "f_defs.h"
  use kinds,            only : dp
  use wvfct,            only :npwx, nbnd, npw, igk, g2kin, et
  use uspp,             only : nkb, vkb
  use wavefunctions_module, only: evc
  use phcom,            only: evq
  use control_ph,       ONLY : nbnd_occ, lgamma

  implicit none
  integer :: ik
  integer :: lter, ibnd, ig, jbnd
  real(kind=DP) :: thresh, anorm, emin, emax, alpha_pv
  logical :: conv_root
  complex(kind=dp) :: dpsi(npwx,nbnd), dev(npwx,nbnd)
  real(kind=DP), allocatable :: h_diag (:,:), eprec(:)
  complex(kind=dp), allocatable ::  auxg(:), ps(:)
  complex(kind=DP) :: ZDOTC
  external ch_psi_all, cg_psi


  allocate (auxg(npwx))
  allocate (h_diag(npwx, nbnd))    
  allocate (eprec(nbnd))
  allocate (ps(nbnd))
  lgamma=.true. !trick to use ch_psi_all as is, check side effect...

 ! starting guess 

  dpsi(:,:)=(0.d0,0.d0)



 !Set the threshold of the diagonalisation

  do ibnd = 1, nbnd_occ (ik)
     auxg(:) = (0.d0, 0.d0)
     do jbnd = 1, nbnd_occ (ik)
        ps(jbnd)=-ZDOTC(npw,evc(1,jbnd),1,dev(1,ibnd),1)
     enddo
#ifdef __PARA
     call reduce (2 * nbnd, ps)
#endif
     do jbnd = 1, nbnd_occ (ik)
        call ZAXPY (npw, ps (jbnd), evc (1, jbnd), 1, auxg, 1)
     enddo
     call DAXPY (2*npw, 1.0d0, auxg, 1, dev (1, ibnd), 1)
  enddo
 !
 !    Here we change the sign of the known term
 !
  call DSCAL (2*npwx*nbnd, -1.d0, dev, 1)

 ! do preconditionning
  thresh=1.d-12
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        auxg (ig) = g2kin (ig) * evc (ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0*ZDOTC(npw,evc(1,ibnd),1,auxg,1)
  enddo
#ifdef __PARA
  call reduce (nbnd_occ (ik), eprec)
#endif
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
     enddo
  enddo

  print *, nbnd_occ(ik)
  
  call cgsolve_all(ch_psi_all, cg_psi, et, dev, dpsi, h_diag, &
      npwx, npw, thresh, ik, lter, conv_root, anorm, &
      nbnd_occ(ik)  )

! print *,'et',et
! print *, 'dpsi',dpsi(1,1)
  print *,'anorm',anorm

 return
  
end subroutine solve_cg

