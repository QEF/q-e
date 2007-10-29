!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us (ik, mode, uact, addnlcc)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp
  use pwcom
  USE noncollin_module, ONLY : npol
  use uspp_param, only: upf
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  complex(DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  logical :: addnlcc
  !
  !   And the local variables
  !

  integer :: na, mu, ikk, ig, nt, ibnd, ir, is, ip
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh

  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable, target :: aux (:)
  complex(DP) , allocatable :: aux1 (:), aux2 (:)
  complex(DP) , pointer :: auxs (:)
  ! work space

  call start_clock ('dvqpsi_us')
  if (nlcc_any) then
     allocate (aux( nrxx))    
     if (doublegrid) then
        allocate (auxs( nrxxs))    
     else
        auxs => aux
     endif
  endif
  allocate (aux1( nrxxs))    
  allocate (aux2( nrxxs))    
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  if (lgamma) then
     ikk = ik
  else
     ikk = 2 * ik - 1
  endif
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) then
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) * eigts3 ( &
                ig3 (ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (nls (ig) ) = aux1 (nls (ig) ) + vlocq (ig, nt) * gu * &
                fact * gtau
        enddo
     endif
  enddo
  !
  ! add NLCC when present
  !
   if (nlcc_any.and.addnlcc) then
      aux(:) = (0.d0, 0.d0)
      do na = 1,nat
         fact = tpiba*(0.d0,-1.d0)*eigqts(na)
         mu = 3*(na-1)
         if (abs(uact(mu+1))+abs(uact(mu+2))  &
                         +abs(uact(mu+3)).gt.1.0d-12) then
            nt=ityp(na)
            u1 = uact(mu+1)
            u2 = uact(mu+2)
            u3 = uact(mu+3)
            gu0 = xq(1)*u1 +xq(2)*u2+xq(3)*u3
            if (upf(nt)%nlcc) then
               do ig = 1,ngm
                  gtau = eigts1(ig1(ig),na)*   &
                         eigts2(ig2(ig),na)*   &
                         eigts3(ig3(ig),na)
                  gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
                  aux(nl(ig))=aux(nl(ig))+drc(ig,nt)*gu*fact*gtau
               enddo
            endif
         endif
      enddo
      call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
      if (.not.lsda) then
         do ir=1,nrxx
            aux(ir) = aux(ir) * dmuxc(ir,1,1)
         end do
      else
         is=isk(ikk)
         do ir=1,nrxx
            aux(ir) = aux(ir) * 0.5d0 *  &
                 (dmuxc(ir,is,1)+dmuxc(ir,is,2))
         enddo
      endif
      call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
      if (doublegrid) then
         auxs(:) = (0.d0, 0.d0)
         do ig=1,ngms
            auxs(nls(ig)) = aux(nl(ig))
         enddo
      endif
      aux1(:) = aux1(:) + auxs(:)
   endif
  !
  ! Now we compute dV_loc/dtau in real space
  !
  call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
  do ibnd = 1, nbnd
     do ip=1,npol
        aux2(:) = (0.d0, 0.d0)
        if (ip==1) then
           do ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig, ibnd)
           enddo
        else
           do ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig+npwx, ibnd)
           enddo
        end if
        !
        !  This wavefunction is computed in real space
        !
        call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
        do ir = 1, nrxxs
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        enddo
        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        if (ip==1) then
           do ig = 1, npwq
              dvpsi (ig, ibnd) = aux2 (nls (igkq (ig) ) )
           enddo
        else
           do ig = 1, npwq
              dvpsi (ig+npwx, ibnd) = aux2 (nls (igkq (ig) ) )
           enddo
        end if
     enddo
  enddo
  !
  deallocate (aux2)
  deallocate (aux1)
  if (nlcc_any) then
     deallocate (aux)
     if (doublegrid) deallocate (auxs)
  endif
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients in the perturbat
  !
  call dvqpsi_us_only (ik, mode, uact)

  call stop_clock ('dvqpsi_us')
  return
end subroutine dvqpsi_us
