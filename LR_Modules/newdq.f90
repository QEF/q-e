!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine newdq (dvscf, npe)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, mill, eigts1, eigts2, eigts3
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE paw_variables,        ONLY : okpaw
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  USE lrus,                 ONLY : int3, int3_paw
  USE qpoint,               ONLY : xq, eigqts
  USE control_lr,           ONLY : lgamma

  implicit none
  !
  !   The dummy variables
  !
  integer, intent(in) :: npe
  ! input: the number of perturbations

  complex(DP), intent(in) :: dvscf (dfftp%nnr, nspin_mag, npe)
  ! input: the change of the selfconsistent pot.
  !
  !   And the local variables
  !
  integer :: na, nb, nab, ig, nt, ir, ipert, is, ih, jh
  ! counters

  real(DP), allocatable :: qmod (:), qg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP), allocatable :: aux1 (:,:), veff (:), qgm(:), sk(:,:)
  ! work space

  COMPLEX(DP) :: tmp
  ! temporary scalar to hold accumulated value

  if (.not.okvan) return
  !
  call start_clock ('newdq')
  !
  allocate (aux1 (ngm , nspin_mag))
  allocate (veff (dfftp%nnr))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm  (ngm))
  allocate (qmod (ngm))
  !
  if (.not.lgamma) allocate (qg (3,  ngm))
  !
  !    first compute the spherical harmonics
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qg)
     call ylmr2 (lmaxq * lmaxq, ngm, qg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) ) * tpiba
     enddo
  else
     call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) ) * tpiba
     enddo
  endif
  !
  !     and for each perturbation of this irreducible representation
  !     integrate the change of the self consistent potential and
  !     the Q functions
  !
  !$acc data create(qgm) copyin(ylmk0, qmod, eigqts) copyout(int3)
  !
  !$acc kernels
  int3 (:,:,:,:,:) = (0.d0, 0.0d0)
  !$acc end kernels
  do ipert = 1, npe

     do is = 1, nspin_mag
        do ir = 1, dfftp%nnr
           veff (ir) = dvscf (ir, is, ipert)
        enddo
        CALL fwfft ('Rho', veff, dfftp)
        do ig = 1, ngm
           aux1 (ig, is) = veff (dfftp%nl (ig) )
        enddo
     enddo
     !$acc data copyin(aux1)

     do nt = 1, ntyp
        if (upf(nt)%tvanp ) then
           !
           ! Count the number of atoms of type nt and allocate sk accordingly
           !
           nab = 0
           do na = 1, nat
              if (ityp(na) == nt) nab = nab + 1
           enddo
           !
           allocate(sk(ngm, nab))
           !$acc data create(sk)
           !
           nb = 0
           do na = 1, nat
              if (ityp(na) == nt) then
                 nb = nb + 1
                 !
                 ! Compute the structure factor for all atoms of type nt
                 !
                 !$acc parallel loop present(eigts1,eigts2,eigts3,mill,eigqts)
                 do ig = 1, ngm
                    sk(ig,nb) = eigts1(mill(1,ig),na) * &
                                eigts2(mill(2,ig),na) * &
                                eigts3(mill(3,ig),na) * &
                                eigqts(na)
                 enddo
              endif
           enddo
           !
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                 nb = 0
                 do na = 1, nat
                    if (ityp (na) == nt) then
                       nb = nb + 1
                       do is = 1, nspin_mag
                          tmp = (0.d0, 0.d0)
                          !$acc parallel loop reduction(+:tmp) present(qgm, sk, aux1)
                          do ig = 1, ngm
                             tmp = tmp + conjg(qgm(ig) * sk(ig,nb)) * aux1(ig,is)
                          enddo
                          !
                          ! Compute and assign int3 on the device. The assignment is a single
                          ! operation hence the serial region below (executed by 1 thread only).
                          !
                          ! The alternative is to keep int3 on the host, and send tmp (scalar)
                          ! from device to host at every (ih,jh,na,is,ipert) iteration.
                          !
                          ! The choice is between many small (scalar) transfers vs a serial region
                          ! and a single copyout of int3 at the end (this is what is done below).
                          !
                          !$acc serial present(int3)
                          int3(ih,jh,na,is,ipert) = omega * tmp
                          !
                          !    We use the symmetry properties of the ps factor:
                          !    int3(jh,ih,...) = int3(ih,jh,...)
                          !
                          int3(jh,ih,na,is,ipert) = omega * tmp
                          !$acc end serial
                       enddo
                    endif
                 enddo
              enddo
           enddo
           !
           !$acc end data ! sk
           deallocate(sk)
           !
        endif
     enddo
     !$acc end data ! aux1
  enddo
  !$acc end data ! int3
#if defined(__MPI)
  call mp_sum ( int3, intra_bgrp_comm )
#endif
  !
  ! Sum of the USPP and PAW terms 
  ! (see last two terms in Eq.(12) in PRB 81, 075123 (2010))
  !
  IF (okpaw) int3 = int3 + int3_paw
  !
  IF (noncolin) CALL set_int3_nc(npe)
  !
  if (.not.lgamma) deallocate (qg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (veff)
  deallocate (aux1)
  !
  call stop_clock ('newdq')
  !
  return
  !
end subroutine newdq
