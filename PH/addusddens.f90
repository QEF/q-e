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
  !  This routine adds to the change of the charge density the part
  !  which is due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term.
  !
#include "machine.h"

  use pwcom
  use phcom
  use parameters, only : DP
  implicit none
  !
  !   the dummy variables
  !

  integer :: iflag, npe
  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(kind=DP) :: drhoscf (nrxx, nspin, npe), dbecsum (nhm * &
       (nhm + 1) / 2, nat, nspin, npe)
  ! inp/out: change of the charge density
  !input: sum over kv of bec
  integer :: irr, mode0
  ! input:the index of the irreducible repr.
  ! input:the mode of the representation
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ir, mu, mode, ipert, is, ijh
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

  real(kind=DP), allocatable  :: qmod (:),&
       qpg (:,:),&
       ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(kind=DP) :: fact, zsum, bb, alpha, alpha_0, u1, u2, u3
  complex(kind=DP), allocatable ::  sk (:),&
       qg (:),&
       drhous (:,:),&
       aux (:,:,:)
  ! auxiliary variables
  ! auxiliary variable
  ! auxiliary variable
  ! auxiliary variables
  ! the structure factor
  ! auxiliary variable for FFT
  ! contain the charge of drho
  ! auxiliary variable for drho(G)

  if (.not.okvan) return
  call start_clock ('addusddens')
  allocate (aux(   ngm , nspin , 3))    
  allocate (sk (  ngm))    
  allocate (qg (  nrxx))    
  allocate (ylmk0(  ngm , lqx * lqx))    
  allocate (qmod (  ngm))    
  if (.not.lgamma) allocate (qpg( 3  , ngm))    
  if (iflag.eq.0) allocate (drhous( nrxx, nspin))    
  !      write(6,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lqx * lqx, ngm, qpg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  else
     call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) )
     enddo

  endif
  fact = 0.5d0 * DCMPLX (0.d0, - tpiba)
  call setv (6 * ngm * nspin, 0.d0, aux, 1)
  do nt = 1, ntyp
     if (tvanp (nt) ) then
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
                       sk (ig) = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) &
                            * eigts3 (ig3 (ig), na) * eigqts (na) * qgm (ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do ipert = 1, npert (irr)
                       do is = 1, nspin
                          mode = mode0 + ipert
                          zsum = dbecsum (ijh, na, is, ipert)
                          u1 = u (mu + 1, mode)
                          u2 = u (mu + 2, mode)
                          u3 = u (mu + 3, mode)
                          if (abs (u1) + abs (u2) + abs (u3) .gt.1d-12.and.iflag.eq.1) &
                               then
                             bb = becsum (ijh, na, is)
                             zsum = zsum + 0.5d0 * (alphasum (ijh, 1, na, is) * u1 + &
                                  alphasum (ijh, 2, na, is) * u2 + alphasum (ijh, 3, na, &
                                  is) * u3)
                             u1 = u1 * fact
                             u2 = u2 * fact
                             u3 = u3 * fact
                             alpha_0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
                             do ig = 1, ngm

                                alpha = alpha_0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, &
                                     ig) * u3

                                aux (ig, is, ipert) = aux (ig, is, ipert) + (zsum + &
                                     alpha * bb) * sk (ig)
                             enddo
                          else
                             call ZAXPY (ngm, zsum, sk, 1, aux (1, is, ipert), &
                                  1)
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
     if (iflag.eq.0) call davcio (drhous, lrdrhous, iudrhous, mu, &
          - 1)
     do is = 1, nspin
        call setv (2 * nrxx, 0.d0, qg, 1)
        do ig = 1, ngm
           qg (nl (ig) ) = aux (ig, is, ipert)
        enddo
        call cft3 (qg, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        if (iflag.eq.0) then
           do ir = 1, nrxx
              drhoscf (ir, is, ipert) = drhoscf (ir, is, ipert) + 2.d0 * qg ( &
                   ir) + drhous (ir, is)
           enddo
        else
           do ir = 1, nrxx
              drhoscf (ir, is, ipert) = drhoscf (ir, is, ipert) + 2.d0 * qg ( &
                   ir)
           enddo
        endif
     enddo

  enddo
  if (iflag.eq.0) deallocate (drhous)
  if (.not.lgamma) deallocate (qpg)
  deallocate (ylmk0)
  deallocate (qmod)
  deallocate (qg)
  deallocate (sk)
  deallocate (aux)

  call stop_clock ('addusddens')
  return
end subroutine addusddens
