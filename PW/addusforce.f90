!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusforce (forcenl)  
  !----------------------------------------------------------------------
  !
  !   This routine computes the contribution to atomic forces due
  !   to the dependence of the Q function on the atomic position.
  !   On output: the contribution is added to forcenl
  !
#include "machine.h"

  use pwcom
  use allocate
  implicit none
  real(kind=DP) :: forcenl (3, nat)  
  ! local variables
  integer :: ig, ir, dim, nt, ih, jh, ijh, ipol, is, na  
  complex(kind=DP):: cfac 
  real(kind=DP) :: fact, DDOT
  ! work space
  complex(kind=DP), pointer :: aux(:,:), aux1(:,:), vg(:)
  real(kind=DP) , pointer :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)
  !
  if (.not.okvan) return
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  call mallocate(aux,ngm,nspin)
  call mallocate(aux1,ngm,3)
  call mallocate(ddeeq, 3, (nhm*(nhm+1))/2,nat,nspin)
  call mallocate(qmod, ngm)  
  call mallocate(ylmk0,ngm,lqx*lqx)
  !
  ddeeq(:,:,:,:) = 0.d0
  !
  call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)  
  !
  qmod (:) = sqrt (gg (:) )  
  !
  ! fourier transform of the total effective potential
  !
  call mallocate(vg,nrxx)
  do is = 1, nspin  
     vg (:) = vltot (:) + vr (:, is)  
     call cft3 (vg, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)  
     aux (:, is) = vg (nl (:) ) * tpiba * (0.d0, - 1.d0)  
  enddo
  call mfree (vg)  
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  !
  do nt = 1, ntyp  
     if (tvanp (nt) ) then  
        ijh = 1  
        do ih = 1, nh (nt)  
           do jh = ih, nh (nt)  
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)  
              do na = 1, nat  
                 if (ityp (na) .eq.nt) then  
                    !
                    ! The product of potential, structure factor and iG
                    !
                    do is = 1, nspin  
                       do ig = 1, ngm  
                          cfac = aux (ig, is) * conjg (eigts1 (ig1 (ig), na) *&
                                                       eigts2 (ig2 (ig), na) *&
                                                       eigts3 (ig3 (ig), na) )
                          aux1 (ig, 1) = g (1, ig) * cfac  
                          aux1 (ig, 2) = g (2, ig) * cfac  
                          aux1 (ig, 3) = g (3, ig) * cfac  
                       enddo
                       !
                       !    and the product with the Q functions
                       !    G=0 term gives no contribution
                       !
                       do ipol = 1, 3  
                          ddeeq (ipol, ijh, na, is) = omega * fact * &
                               DDOT (2 * ngm, aux1 (1, ipol), 1, qgm, 1)
                       enddo
                    enddo
                 endif
              enddo
              ijh = ijh + 1  
           enddo
        enddo
     endif

  enddo
#ifdef PARA
  call reduce (3 * nhm * (nhm + 1) * nat * nspin / 2, ddeeq)  
#endif
  !            write(6,'( "dmatrix atom ",i4)') na
  !            do ih = 1, nh(nt)
  !               write(6,'(8f9.4)') (ddeeq(ipol,ih,jh,na),jh=1,nh(nt))
  !            end do
  !            write(6,'( "dion pseudo ",i4)') nt
  !            do ih = 1, nh(nt)
  !               write(6,'(8f9.4)') (dvan(ih,jh,nt),jh=1,nh(nt))
  !            end do
  do is = 1, nspin  
     do na = 1, nat  
        nt = ityp (na)  
        dim = (nh (nt) * (nh (nt) + 1) ) / 2  
        do ipol = 1, 3  
           do ir = 1, dim  
              forcenl (ipol, na) = forcenl (ipol, na) + &
                   ddeeq (ipol, ir, na, is) * becsum (ir, na, is)
           enddo
        enddo
     enddo
  enddo
  call mfree (ylmk0)  
  call mfree (qmod)  
  call mfree (ddeeq)  
  call mfree (aux1)  
  call mfree (aux)  

  return  
end subroutine addusforce

