!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine newd
  !----------------------------------------------------------------------
  !
  !   This routine computes the integral of the effective potential with
  !   the Q function and adds it to the bare ionic D term which is used
  !   to compute the non-local term in the US scheme.
  !
#include "machine.h"

  use pwcom
  USE wavefunctions,    ONLY : psic
  implicit none
  integer :: ig, nt, ih, jh, na, is
  ! counters on g vectors, atom type, beta functions x 2, atoms, spin
  complex(kind=DP), allocatable :: aux (:,:), qgm_na (:)
  ! work space
  real(kind=DP), allocatable  :: ylmk0 (:,:), qmod (:)
  ! spherical harmonics, modulus of G
  real(kind=DP) :: fact, DDOT
  !
  !
  if (.not.okvan) then
     ! no ultrasoft potentials: use bare coefficients for projectors
     do na = 1, nat
        nt = ityp (na)
        do is = 1, nspin
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 deeq (ih, jh, na, is) = dvan (ih, jh, nt)
                 deeq (jh, ih, na, is) = deeq (ih, jh, na, is)
              enddo
           enddo
        enddo
     end do
     return
  end if

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  call start_clock ('newd')
  allocate ( aux(ngm,nspin), qgm_na(ngm), qmod(ngm), ylmk0(ngm, lqx*lqx) )
  !
  deeq(:,:,:,:) = 0.d0
  !
  call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  !
  ! fourier transform of the total effective potential
  !
#ifdef DEBUG_NEWD
   call start_clock ('newd:fftvg')
#endif
  do is = 1, nspin
     psic (:) = vltot (:) + vr (:, is)
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     do ig = 1, ngm
        aux (ig, is) = psic (nl (ig) )
     enddo
  enddo
#ifdef DEBUG_NEWD
   call stop_clock ('newd:fftvg')
#endif
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G exp(-iR.G) Q_nm v^*
  !
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
#ifdef DEBUG_NEWD
   call start_clock ('newd:qvan2')
#endif
              !
              ! The Q(r) for this atomic species without structure factor
              !
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
#ifdef DEBUG_NEWD
   call stop_clock ('newd:qvan2')
#endif
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
#ifdef DEBUG_NEWD
   call start_clock ('newd:int1')
#endif
                    !
                    ! The Q(r) for this specific atom
                    !
                    do ig = 1, ngm
                       qgm_na (ig) = qgm(ig) * eigts1 (ig1(ig), na) &
                                             * eigts2 (ig2(ig), na) &
                                             * eigts3 (ig3(ig), na) 
                     enddo
#ifdef DEBUG_NEWD
   call stop_clock ('newd:int1')
#endif
#ifdef DEBUG_NEWD
   call start_clock ('newd:int2')
#endif
                    !
                    !  and the product with the Q functions
                    !
                    do is = 1, nspin
                       deeq (ih, jh, na, is) = fact * omega * &
                            DDOT (2 * ngm, aux(1,is), 1, qgm_na, 1)
                       if (gamma_only .and. gstart==2) &
                           deeq (ih, jh, na, is) = deeq (ih, jh, na, is) - &
                                             omega*real(aux(1,is)*qgm_na(1))
                    enddo
#ifdef DEBUG_NEWD
   call stop_clock ('newd:int2')
#endif
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
#ifdef __PARA
  call reduce (nhm * nhm * nat * nspin, deeq)
#endif
  do na = 1, nat
     nt = ityp (na)
     do is = 1, nspin
        !           WRITE( stdout,'( "dmatrix atom ",i4, " spin",i4)') na,is
        !           do ih = 1, nh(nt)
        !              WRITE( stdout,'(8f9.4)') (deeq(ih,jh,na,is),jh=1,nh(nt))
        !           end do
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              deeq (ih, jh, na, is) = deeq (ih, jh, na, is) + dvan (ih, jh, nt)
              deeq (jh, ih, na, is) = deeq (ih, jh, na, is)
           enddo
        enddo
     enddo
     !        WRITE( stdout,'( "dion pseudo ",i4)') nt
     !        do ih = 1, nh(nt)
     !           WRITE( stdout,'(8f9.4)') (dvan(ih,jh,nt),jh=1,nh(nt))
     !        end do

  enddo
  deallocate ( aux, qgm_na, qmod, ylmk0 )
  call stop_clock ('newd')

  return
end subroutine newd

