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
subroutine newd
  !----------------------------------------------------------------------
  !
  !   This routine computes the integral of the effective potential with
  !   the Q function and adds it to the bare ionic D term which is used
  !   to compute the non-local term in the US scheme.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, g, gg, &
                                   ngm, gstart, ig1, ig2, ig3, &
                                   eigts1, eigts2, eigts3, nl
  USE lsda_mod,             ONLY : nspin
  USE spin_orb,             ONLY : lspinorb
  USE noncollin_module,     ONLY : noncolin
  USE scf,                  ONLY : vr, vltot
  USE us,                   ONLY : okvan
  USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so
  USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none
  !
  integer :: ig, nt, ih, jh, na, is, nht
  ! counters on g vectors, atom type, beta functions x 2, atoms, spin
  complex(kind=DP), allocatable :: aux (:,:), qgm(:), qgm_na (:)
  ! work space
  real(kind=DP), allocatable  :: ylmk0 (:,:), qmod (:)
  ! spherical harmonics, modulus of G
  real(kind=DP) :: fact, DDOT
  !
  !
  if (.not.okvan) then
     ! no ultrasoft potentials: use bare coefficients for projectors
     do na = 1, nat
        nt=ityp(na)
        nht = nh (nt)
        if (lspinorb) then
           deeq_nc (1:nht,1:nht, na, 1:nspin) = &
                dvan_so (1:nht, 1:nht, 1:nspin, nt)
        else if (noncolin) then
           deeq_nc (1:nht,1:nht, na, 1) = dvan (1:nht, 1:nht, nt)
           deeq_nc (1:nht,1:nht, na, 2) = (0.d0, 0.d0)
           deeq_nc (1:nht,1:nht, na, 3) = (0.d0, 0.d0)
           deeq_nc (1:nht,1:nht, na, 4) = dvan (1:nht, 1:nht, nt)
        else
           do is=1, nspin
              deeq (1:nht,1:nht, na, is) = dvan (1:nht, 1:nht, nt)
           end do
        end if
     end do
     return
  end if

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  call start_clock ('newd')
  allocate ( aux(ngm,nspin), qgm_na(ngm), qgm(ngm), qmod(ngm), &
       ylmk0(ngm, lmaxq*lmaxq) )
  !
  deeq(:,:,:,:) = 0.d0
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  !
  ! fourier transform of the total effective potential
  !
  do is = 1, nspin
     if (nspin == 4 .and. is /= 1) then 
        psic (:) = vr (:, is)
     else
        psic (:) = vltot (:) + vr (:, is)
     endif
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     do ig = 1, ngm
        aux (ig, is) = psic (nl (ig) )
     enddo
  enddo
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G exp(-iR.G) Q_nm v^*
  !
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              !
              ! The Q(r) for this atomic species without structure factor
              !
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              !
              do na = 1, nat
                 if (ityp (na) == nt) then
                    !
                    ! The Q(r) for this specific atom
                    !
                    do ig = 1, ngm
                       qgm_na (ig) = qgm(ig) * eigts1 (ig1(ig), na) &
                                             * eigts2 (ig2(ig), na) &
                                             * eigts3 (ig3(ig), na) 
                    enddo
                    !
                    !  and the product with the Q functions
                    !
                    do is = 1, nspin
                       deeq (ih, jh, na, is) = fact * omega * &
                            DDOT (2 * ngm, aux(1,is), 1, qgm_na, 1)
                       if (gamma_only .and. gstart==2) &
                           deeq (ih, jh, na, is) = deeq (ih, jh, na, is) - &
                                             omega*real(aux(1,is)*qgm_na(1))
                       deeq (jh,ih,na,is) = deeq (ih,jh,na,is)
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo

  deallocate (aux, qgm_na, qgm, qmod, ylmk0)

#ifdef __PARA
  call reduce (nhm * nhm * nat * nspin, deeq)
#endif

  if (lspinorb) then
     call newd_so ( )
  else if (noncolin) then 
     call newd_nc ( )
  else
      do na = 1, nat
        nt = ityp (na)
        do is = 1, nspin
        !           WRITE( stdout,'( "dmatrix atom ",i4, " spin",i4)') na,is
        !           do ih = 1, nh(nt)
        !              WRITE( stdout,'(8f9.4)') (deeq(ih,jh,na,is),jh=1,nh(nt))
        !           end do
           do ih = 1, nh (nt)
              do jh = 1, nh (nt)
                   deeq (ih, jh, na, is) = deeq (ih, jh, na, is) &
                &                      + dvan (ih, jh, nt)
              enddo
           enddo
        enddo
     end do
  end if

  call stop_clock ('newd')

  return

CONTAINS
  !
  SUBROUTINE newd_so ( )
    !
    USE spin_orb, ONLY : fcoef
    implicit none
    integer :: ijs, is1, is2, kh, lh
    
    do na = 1, nat
       nt = ityp (na)
       ijs=0
       do is1=1,2
          do is2=1,2
             ijs=ijs+1
             do ih = 1, nh(nt)
                do jh = 1, nh(nt)
                   deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                   do kh= 1, nh(nt)
                      do lh= 1, nh(nt)
                         deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                             deeq (kh,lh,na,1)*            &
                           (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                            fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                                             deeq (kh,lh,na,2)*            &
                           (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  + &
                            fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                (0.d0,-1.d0)*deeq (kh,lh,na,3)*            &
                           (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  - &
                            fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                             deeq (kh,lh,na,4)*            &
                           (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  - &
                            fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))   
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    return
  end subroutine newd_so
  !
  SUBROUTINE newd_nc ( )
    !
    implicit none

    do na = 1, nat
       nt = ityp (na)
       do is = 1, nspin
          do ih = 1, nh (nt)
             do jh = 1, nh (nt)
                deeq_nc(ih, jh, na, 1) = dvan(ih,jh,nt) + &
                     deeq (ih,jh,na,1) + deeq (ih,jh,na,4)
                deeq_nc(ih, jh, na, 2) = &
                     deeq (ih,jh,na,2) - (0.d0,1.d0)*deeq (ih,jh,na,3)
                deeq_nc(ih, jh, na, 3) = &
                     deeq (ih,jh,na,2) + (0.d0,1.d0)*deeq (ih,jh,na,3)
                deeq_nc(ih, jh, na, 4) = dvan(ih,jh,nt) + &
                     deeq (ih,jh,na,1) - deeq (ih,jh,na,4)
             enddo
          enddo
       enddo
    enddo
  end subroutine newd_nc

end subroutine newd
