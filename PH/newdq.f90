!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine newdq (dvscf, irr, npe)  
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
#include "machine.h"

  use pwcom 
  use allocate 
  use parameters, only : DP 
  use phcom  
  implicit none 
  !
  !   The dummy variables
  !

  integer :: npe  
  ! input: the number of perturbations

  complex(kind=DP) :: dvscf (nrxx, nspin, npe)  
  ! input: the change of the self
  !consistent pot.
  integer :: irr  
  ! input: the irreducible representation
  !
  !   And the local variables
  !

  integer :: na, ig, nt, ir, ipert, is, ih, jh  
  ! counter on atoms
  ! counter on G vectors
  ! counter on atomic types
  ! counter on real mesh
  ! counter on perturbations
  ! counter on spin
  ! counter on beta functions
  ! counter on beta functions

  real(kind=DP), pointer :: qmod (:), qg (:,:), ylmk0 (:,:)  
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(kind=DP) :: ZDOTC  
  ! the scalar product function

  complex(kind=DP), pointer :: aux1 (:), aux2 (:,:), veff (:)  
  ! space for several quantities
  ! space for veff
  ! a mesh space for the FFT of the V_eff

  if (.not.okvan) return  
  call setv (2 * nhm * nhm * 3 * nat * nspin, 0.0d0, int3, 1)  

  call start_clock ('newdq')  
  call mallocate(aux1 ,  ngm)  
  call mallocate(aux2 ,  ngm , nspin)  
  call mallocate(veff ,  nrxx)  
  call mallocate(ylmk0 ,  ngm , lqx * lqx)  
  call mallocate(qmod ,  ngm)  

  if (.not.lgamma) call mallocate(qg ,3,  ngm)  
  !
  !    first compute the spherical harmonics
  !
  if (.not.lgamma) then  
     call setqmod (ngm, xq, g, qmod, qg)  
     call ylmr2 (lqx * lqx, ngm, qg, qmod, ylmk0)  
     do ig = 1, ngm  
        qmod (ig) = sqrt (qmod (ig) )  
     enddo
  else  
     call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)  
     do ig = 1, ngm  
        qmod (ig) = sqrt (gg (ig) )  
     enddo

  endif
  !
  !     and for each perturbation of this irreducible representation
  !     integrate the change of the self consistent potential and
  !     the Q functions
  !
  do ipert = 1, npert (irr)  
     do is = 1, nspin  
        do ir = 1, nrxx  
           veff (ir) = dvscf (ir, is, ipert)  
        enddo
        call cft3 (veff, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)  
        do ig = 1, ngm  
           aux2 (ig, is) = veff (nl (ig) )  
        enddo

     enddo
     do nt = 1, ntyp  
        if (tvanp (nt) ) then  
           do ih = 1, nh (nt)  
              do jh = ih, nh (nt)  
                 call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)  
                 do na = 1, nat  
                    if (ityp (na) .eq.nt) then  
                       do ig = 1, ngm  
                          aux1 (ig) = qgm (ig) * eigts1 (ig1 (ig), na) * eigts2 (ig2 ( &
                               ig), na) * eigts3 (ig3 (ig), na) * eigqts (na)

                       enddo
                       do is = 1, nspin  
                          int3 (ih, jh, ipert, na, is) = omega * ZDOTC (ngm, aux1, 1, &
                               aux2 (1, is), 1)
                       enddo
                       !
                       !  ps contain the integral of V_loc and Q_nm
                       !
                    endif
                 enddo
              enddo

           enddo
           do na = 1, nat  
              if (ityp (na) .eq.nt) then  
                 !
                 !    We use the symmetry properties of the ps factor
                 !
                 do ih = 1, nh (nt)  
                    do jh = ih, nh (nt)  
                       do is = 1, nspin  
                          int3 (jh, ih, ipert, na, is) = int3 (ih, jh, ipert, na, is)  
                       enddo
                    enddo
                 enddo

              endif
           enddo
        endif
     enddo
  enddo
#ifdef PARA
  call reduce (2 * nhm * nhm * 3 * nat * nspin, int3)  
#endif
  !      do ih = 1,nh(1)
  !         do jh=1,nh(1)
  !            write(6,*) int3(jh,ih,1,1,1)
  !         enddo
  !      enddo
  !      call stop_ph(.true.)
  if (.not.lgamma) call mfree (qg)  
  call mfree (qmod)  
  call mfree (ylmk0)  
  call mfree (veff)  
  call mfree (aux2)  
  call mfree (aux1)  

  call stop_clock ('newdq')  
  return  
end subroutine newdq
