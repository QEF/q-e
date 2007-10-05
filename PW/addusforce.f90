!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine addusforce (forcenl)
  !----------------------------------------------------------------------
  !
  !   This routine computes the contribution to atomic forces due
  !   to the dependence of the Q function on the atomic position.
  !   On output: the contribution is added to forcenl
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, &
                         nl, nlm, gg, g, eigts1, eigts2, eigts3, ig1, ig2, ig3
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : vr, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE wvfct,      ONLY : gamma_only
  !
  implicit none
  !
  real(DP) :: forcenl (3, nat)
  ! local variables
  integer :: ig, ir, dim, nt, ih, jh, ijh, ipol, is, na
  complex(DP):: cfac
  real(DP) :: fact, DDOT
  ! work space
  complex(DP), allocatable :: aux(:,:), aux1(:,:), vg(:), qgm(:)
  real(DP) , allocatable :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)

  !
  if (.not.okvan) return
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  allocate (aux(ngm,nspin))    
  !
  ! fourier transform of the total effective potential
  !
  allocate (vg(nrxx))    
  do is = 1, nspin
     if (nspin.eq.4.and.is.ne.1) then
        vg (:) = vr(:,is)
     else
        vg (:) = vltot (:) + vr (:, is)
     endif
     call cft3 (vg, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     aux (:, is) = vg (nl (:) ) * tpiba * (0.d0, -1.d0)
  enddo
  deallocate (vg)
  !
  allocate (aux1(ngm,3))    
  allocate (ddeeq( 3, (nhm*(nhm+1))/2,nat,nspin))    
  allocate (qgm( ngm))
  allocate (qmod( ngm))    
  allocate (ylmk0(ngm,lmaxq*lmaxq))    
  !
  ddeeq(:,:,:,:) = 0.d0
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  !
  qmod (:) = sqrt (gg (:) )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  !
  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        ijh = 1
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              do na = 1, nat
                 if (ityp (na) == nt) then
                    !
                    ! The product of potential, structure factor and iG
                    !
                    do is = 1, nspin
                       do ig = 1, ngm
                          cfac = aux (ig, is) * CONJG(eigts1 (ig1 (ig), na) *&
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
#ifdef __PARA
  call reduce (3 * nhm * (nhm + 1) * nat * nspin / 2, ddeeq)
#endif
  !            WRITE( stdout,'( "dmatrix atom ",i4)') na
  !            do ih = 1, nh(nt)
  !               WRITE( stdout,'(8f9.4)') (ddeeq(ipol,ih,jh,na),jh=1,nh(nt))
  !            end do
  !            WRITE( stdout,'( "dion pseudo ",i4)') nt
  !            do ih = 1, nh(nt)
  !               WRITE( stdout,'(8f9.4)') (dvan(ih,jh,nt),jh=1,nh(nt))
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
  deallocate (ylmk0)
  deallocate (qgm)
  deallocate (qmod)
  deallocate (ddeeq)
  deallocate (aux1)
  deallocate (aux)

  return
end subroutine addusforce

