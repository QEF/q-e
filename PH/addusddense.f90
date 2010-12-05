!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddense (drhoscf, dbecsum)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the change of the charge and magnetization
  !  densities due to an electric field perturbation
  !  the part due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term.
  !  The expression implemented is given in Eq. B32 of PRB 64, 235118
  !  (2001) with b=c=0.
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  USE gvect, ONLY : nl, g, gg, ngm, eigts1, eigts2, eigts3, mill
  USE uspp, ONLY: okvan
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE noncollin_module, ONLY : nspin_mag

  USE qpoint, ONLY : eigqts
  implicit none
  !
  !   the dummy variables
  !

  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(DP) :: drhoscf(dfftp%nnr,nspin_mag,3), &
                 dbecsum(nhm*(nhm+1)/2,nat,nspin_mag,3)

  ! inp/out: change of the charge density
  ! input: sum over kv of bec
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ipert, ijh, is

  ! counters

  real(DP), allocatable  :: qmod(:), ylmk0(:,:)
  ! the modulus of q+G
  ! the spherical harmonics

  complex(DP) :: zsum
  complex(DP), allocatable ::  sk (:), qg (:), qgm (:), aux (:,:,:)
  ! the structure factor
  ! work space

  if (.not.okvan) return
  call start_clock ('addusddense')
  allocate (aux(  ngm, nspin_mag, 3))
  allocate (sk (  ngm))
  allocate (qg (  dfftp%nnr))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm  (ngm))
  allocate (qmod (ngm))

  !
  !  And then we compute the additional charge in reciprocal space
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo

  aux (:,:,:) = (0.d0, 0.d0)
  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) == nt) then
                    !
                    ! calculate the structure factor
                    !
                    do ig = 1, ngm
                       sk(ig)=eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na) &
                             *eigts3(mill(3,ig),na)*eigqts(na)*qgm(ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do is=1,nspin_mag
                       do ipert = 1, 3
                          zsum = dbecsum (ijh, na, is,ipert)
                          call zaxpy(ngm,zsum,sk,1,aux(1,is,ipert),1)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space
  !
  do is=1,nspin_mag
     do ipert = 1, 3
        qg (:) = (0.d0, 0.d0)
        qg (nl (:) ) = aux (:, is, ipert)
        CALL invfft ('Dense', qg, dfftp)
        drhoscf(:,is,ipert) = drhoscf(:,is,ipert) + 2.d0*qg(:)
     enddo
  enddo
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (qg)
  deallocate (sk)
  deallocate (aux)

  call stop_clock ('addusddense')
  return
end subroutine addusddense
