!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddens (drhoscf, dbecsum, irr, mode0, npe, iflag)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the change of the charge and of the
  !  magnetization densities the part due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE wavefunctions_module,  ONLY: psic
  use phcom
  USE kinds, only : DP
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm

  implicit none
  !
  !   the dummy variables
  !

  integer :: iflag, npe
  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(DP) :: drhoscf (nrxx, nspin, npe), &
                      dbecsum (nhm*(nhm+1)/2, nat, nspin, npe)
  ! inp/out: change of the charge density
  !input: sum over kv of bec
  integer :: irr, mode0
  ! input:the index of the irreducible repr.
  ! input:the mode of the representation
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ir, mu, mode, ipert, is, ijh, nspin0
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

  real(DP), allocatable  :: qmod (:), qpg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP) :: fact, zsum, bb, alpha, alpha_0, u1, u2, u3
  ! auxiliary variables
  complex(DP), allocatable ::  sk (:), qgm(:), drhous (:,:), aux (:,:,:)
  ! the structure factor
  ! q_lm(G)
  ! contain the charge of drho
  ! auxiliary variable for drho(G)

  if (.not.okvan) return
  nspin0=nspin
  if (nspin==4.and..not.domag) nspin0=1
  call start_clock ('addusddens')
  allocate (aux(  ngm , nspin , npertx))    
  allocate (sk (  ngm))    
  allocate (ylmk0(ngm , lmaxq * lmaxq))    
  allocate (qgm(  ngm))    
  allocate (qmod( ngm))    
  if (.not.lgamma) allocate (qpg( 3  , ngm))    
  !      WRITE( stdout,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  else
     call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) )
     enddo
  endif
  fact = 0.5d0 * CMPLX (0.d0, - tpiba)
  aux(:,:,:) = (0.d0, 0.d0)
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    mu = 3 * (na - 1)
                    !
                    ! calculate the structure factor
                    !
                    do ig = 1, ngm
                       sk (ig) = eigts1 (ig1 (ig), na) * &
                                 eigts2 (ig2 (ig), na) * &
                                 eigts3 (ig3 (ig), na) * &
                                 eigqts (na) * qgm (ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do ipert = 1, npert (irr)
                       do is = 1, nspin0
                          mode = mode0 + ipert
                          zsum = dbecsum (ijh, na, is, ipert)
                          u1 = u (mu + 1, mode)
                          u2 = u (mu + 2, mode)
                          u3 = u (mu + 3, mode)
                          if (abs(u1) + abs(u2) + abs(u3) .gt.1d-12 .and. &
                              iflag.eq.1) then
                             bb = becsum (ijh, na, is)
                             zsum = zsum + 0.5d0 * &
                                  ( alphasum (ijh, 1, na, is) * u1 &
                                  + alphasum (ijh, 2, na, is) * u2 &
                                  + alphasum (ijh, 3, na, is) * u3)
                             u1 = u1 * fact
                             u2 = u2 * fact
                             u3 = u3 * fact
                             alpha_0 = xq(1)*u1 + xq(2)*u2 + xq(3)*u3
                             do ig = 1, ngm
                                alpha = alpha_0 + &
                                        g(1,ig)*u1 + g(2,ig)*u2 + g(3,ig)*u3
                                aux(ig,is,ipert) = aux(ig,is,ipert) + &
                                                   (zsum + alpha*bb) * sk(ig)
                             enddo
                          else
                             call ZAXPY (ngm, zsum, sk, 1, aux(1,is,ipert), 1)
                          endif
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
  do ipert = 1, npert (irr)
     mu = mode0 + ipert
     do is = 1, nspin0
        psic(:) = (0.d0, 0.d0)
        do ig = 1, ngm
           psic (nl (ig) ) = aux (ig, is, ipert)
        enddo
        call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        call DAXPY (2*nrxx, 2.d0, psic, 1, drhoscf(1,is,ipert), 1)
     enddo
  enddo
  if (.not.lgamma) deallocate (qpg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (sk)
  deallocate (aux)

  if (iflag == 0) then
     allocate (drhous( nrxx, nspin))    
     do ipert = 1, npert (irr)
        mu = mode0 + ipert
        call davcio (drhous, lrdrhous, iudrhous, mu, -1)
        call DAXPY (2*nrxx*nspin, 1.d0, drhous, 1, drhoscf(1,1,ipert), 1)
     end do
     deallocate (drhous)
  end if

  call stop_clock ('addusddens')
  return
end subroutine addusddens
