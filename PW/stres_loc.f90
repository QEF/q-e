!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_loc (sigmaloc)  
  !----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom  
  use allocate
  implicit none  !
  real(kind=DP) :: sigmaloc (3, 3)
  real(kind=DP) , allocatable :: dvloc(:)
  real(kind=DP) :: evloc, fact
  integer :: ng, nt, l, m, is  
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components

  allocate(dvloc(ngl))
  sigmaloc(:,:) = 0.d0
  call setv (2 * nrxx, 0.d0, psic, 1)  
  do is = 1, nspin  
     call DAXPY (nrxx, 1.d0, rho (1, is), 1, psic, 2)  
  enddo

  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)  
  ! psic contains now the charge density in G space
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  evloc = 0.0  
  do nt = 1, ntyp  
     if (gstart==2) evloc = evloc + &
          psic (nl (1) ) * strf (1, nt) * vloc (igtongl (1), nt)
     do ng = gstart, ngm  
        evloc = evloc + DREAL (conjg (psic (nl (ng) ) ) * strf (ng, nt) ) &
             * vloc (igtongl (ng), nt) * fact
     enddo
  enddo
  !      write (6,*) ' evloc ', evloc, evloc*omega
  !
  do nt = 1, ntyp  
     ! dvloc contains dV_loc(G)/dG
     call dvloc_of_g (lloc (nt), lmax (nt), numeric (nt), mesh (nt), &
          msh (nt), rab (1, nt), r (1, nt), vnl (1, lloc (nt), nt), &
          cc (1, nt), alpc (1, nt), nlc (nt), nnl (nt), zp (nt), &
          aps (1, 0, nt), alps (1, 0, nt), tpiba2, ngl, gl, omega, dvloc)
     ! no G=0 contribution
     do ng = 1, ngm  
        do l = 1, 3  
           do m = 1, l  
              sigmaloc(l, m) = sigmaloc(l, m) + DREAL( conjg( psic(nl(ng) ) ) &
                    * strf (ng, nt) ) * 2.0 * dvloc (igtongl (ng) ) &
                    * tpiba2 * g (l, ng) * g (m, ng) * fact
           enddo
        enddo
     enddo
  enddo
  !
  do l = 1, 3  
     sigmaloc (l, l) = sigmaloc (l, l) + evloc  
     do m = 1, l - 1  
        sigmaloc (m, l) = sigmaloc (l, m)  
     enddo
  enddo
#ifdef PARA
  call reduce (9, sigmaloc)  
#endif
  deallocate(dvloc)  
  return  
end subroutine stres_loc

