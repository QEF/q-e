!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, mill, eigts1, eigts2, eigts3, nl
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE paw_variables,        ONLY : okpaw

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  USE lrus,       ONLY : int3, int3_paw
  USE qpoint,     ONLY : xq, eigqts
  USE control_lr, ONLY : lgamma

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
  integer :: na, ig, nt, ir, ipert, is, ih, jh
  ! countera

  real(DP), allocatable :: qmod (:), qg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP), external :: zdotc
  ! the scalar product function

  complex(DP), allocatable :: aux1 (:), aux2 (:,:), veff (:), qgm(:)
  ! work space

  if (.not.okvan) return
  call start_clock ('newdq')

  int3 (:,:,:,:,:) = (0.d0, 0.0d0)
  allocate (aux1 (ngm))
  allocate (aux2 (ngm , nspin_mag))
  allocate (veff (dfftp%nnr))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm  (ngm))
  allocate (qmod (ngm))

  if (.not.lgamma) allocate (qg (3,  ngm))
  !
  !    first compute the spherical harmonics
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qg)
     call ylmr2 (lmaxq * lmaxq, ngm, qg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  else
     call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) )
     enddo
  endif
  !
  !     and for each perturbation of this irreducible representation
  !     integrate the change of the self consistent potential and
  !     the Q functions
  !
  do ipert = 1, npe

     do is = 1, nspin_mag
        do ir = 1, dfftp%nnr
           veff (ir) = dvscf (ir, is, ipert)
        enddo
        CALL fwfft ('Dense', veff, dfftp)
        do ig = 1, ngm
           aux2 (ig, is) = veff (nl (ig) )
        enddo
     enddo

     do nt = 1, ntyp
        if (upf(nt)%tvanp ) then
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                 do na = 1, nat
                    if (ityp (na) == nt) then
                       do ig = 1, ngm
                          aux1(ig) = qgm(ig) * eigts1(mill(1,ig),na) * &
                                               eigts2(mill(2,ig),na) * &
                                               eigts3(mill(3,ig),na) * &
                                               eigqts(na)
                       enddo
                       do is = 1, nspin_mag
                          int3(ih,jh,na,is,ipert) = omega * &
                                             zdotc(ngm,aux1,1,aux2(1,is),1)
                       enddo
                    endif
                 enddo
              enddo
           enddo
           do na = 1, nat
              if (ityp(na) == nt) then
                 !
                 !    We use the symmetry properties of the ps factor
                 !
                 do ih = 1, nh (nt)
                    do jh = ih, nh (nt)
                       do is = 1, nspin_mag
                          int3(jh,ih,na,is,ipert) = int3(ih,jh,na,is,ipert)
                       enddo
                    enddo
                 enddo
              endif
           enddo
        endif
     enddo

  enddo
#if defined(__MPI)
  call mp_sum ( int3, intra_bgrp_comm )
#endif
  IF (noncolin) CALL set_int3_nc(npe)
  IF (okpaw) int3=int3+int3_paw

  if (.not.lgamma) deallocate (qg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (veff)
  deallocate (aux2)
  deallocate (aux1)

  call stop_clock ('newdq')
  return
end subroutine newdq
