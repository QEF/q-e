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
  !  This routine adds to the change of the charge density the part
  !  which is due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term.
  !
#include "machine.h"

  use pwcom
  use phcom
  USE kinds, only : DP
  implicit none
  !
  !   the dummy variables
  !

  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(kind=DP) :: drhoscf(nrxx,nspin,3), dbecsum(nhm*(nhm+1)/2,nat,nspin,3)

  ! inp/out: change of the charge density
  ! input: sum over kv of bec
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ir, mode, ipert, ijh, is

  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! pointer on modes
  ! pointer on the mode
  ! counter on perturbations
  ! counter on spin
  ! counter on combined beta functions

  real(kind=DP), allocatable  :: qmod(:), ylmk0(:,:)
  ! the modulus of q+G
  ! the spherical harmonics

  complex(kind=DP) :: zsum
  complex(kind=DP), allocatable ::  sk (:), qg (:), aux (:,:,:)
  ! the structure factor
  ! work space

  if (.not.okvan) return
  call start_clock ('addusddense')
  allocate (aux(  ngm, nspin, 3))    
  allocate (sk (  ngm))    
  allocate (qg (  nrxx))    
  allocate (ylmk0(  ngm , lqx * lqx))    
  allocate (qmod (  ngm))    
  !
  !  And then we compute the additional charge in reciprocal space
  !
  call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo

  aux (:,:,:) = (0.d0, 0.d0)
  do nt = 1, ntyp
     if (tvanp (nt) ) then
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
                       sk(ig)=eigts1(ig1(ig),na)*eigts2(ig2(ig),na) &
                             *eigts3(ig3(ig),na)*eigqts(na)*qgm(ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do is=1,nspin
                       do ipert = 1, 3
                          zsum = dbecsum (ijh, na, is,ipert)
                          call ZAXPY(ngm,zsum,sk,1,aux(1,is,ipert),1)
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
  do is=1,nspin
     do ipert = 1, 3
        qg (:) = (0.d0, 0.d0)
        qg (nl (:) ) = aux (:, is, ipert)
        call cft3 (qg, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        drhoscf(:,is,ipert) = drhoscf(:,is,ipert) + 2.d0*qg(:)
     enddo
  enddo
  deallocate (ylmk0)
  deallocate (qmod)
  deallocate (qg)
  deallocate (sk)
  deallocate (aux)

  call stop_clock ('addusddense')
  return
end subroutine addusddense
