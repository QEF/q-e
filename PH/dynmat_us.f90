!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dynmat_us
  !-----------------------------------------------------------------------
  !
  !  This routine calculates the electronic term: <psi|V"|psi>
  !  of the dynamical matrix.
  !
#include "machine.h"

  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE io_files, ONLY: iunigk
  USE kinds, only : DP
  use phcom
#ifdef __PARA
  use para
#endif

  implicit none
  integer :: icart, jcart, na_icart, na_jcart, na, nb, ng, nt, ik, &
       ig, ir, is, ibnd, nu_i, nu_j, ijkb0, ikb, jkb, ih, jh, ikk
  ! counter on polarizations
  ! counter on polarizations
  ! counter on modes
  ! counter on modes
  ! counter on atoms
  ! counter on atoms
  ! counter on G vectors
  ! counter on atomic types
  ! counter on k points
  ! counter on G vectors
  ! counter on real space mesh
  ! counter on spin
  ! counter on bands
  ! counter on modes
  ! counter on modes
  ! initial value of ikb
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! record position of wfc at k

  real(kind=DP) :: gtau, fac, wgg
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  complex(kind=DP) :: work, dynwrk (3 * nat, 3 * nat)
  ! auxiliary space
  ! work space
  complex(kind=DP), allocatable :: rhog (:), &
       gammap (:,:,:,:), aux1 (:,:), work1 (:), work2 (:)
  ! fourier transform of rho
  ! the second derivative of the beta
  ! auxiliary space
  ! auxiliary space
  ! auxiliary space

  call start_clock ('dynmat_us')
  allocate (rhog  ( nrxx))    
  allocate (work1 ( npwx))    
  allocate (work2 ( npwx))    
  allocate (aux1  ( npwx , nbnd))    
  allocate (gammap(  nkb, nbnd , 3 , 3))    

  call setv (2 * 3 * 3 * nat * nat, 0.0d0, dynwrk, 1)
  !
  !   We first compute the part of the dynamical matrix due to the local
  !   potential
#ifdef __PARA
  !   ... only the first pool does the calculation (no sum over k needed)
  if (mypool.ne.1) goto 100
#endif
  !
  call setv (2 * nrxx, 0.d0, rhog, 1)
  do is = 1, nspin
     do ir = 1, nrxx
        rhog (ir) = rhog (ir) + DCMPLX (rho (ir, is), 0.d0)
     enddo
  enddo

  call cft3 (rhog, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  ! there is a delta ss'
  !
  do na = 1, nat
     do icart = 1, 3
        na_icart = 3 * (na - 1) + icart
        do jcart = 1, 3
           na_jcart = 3 * (na - 1) + jcart
           do ng = 1, ngm
              gtau = tpi * (g (1, ng) * tau (1, na) + &
                            g (2, ng) * tau (2, na) + &
                            g (3, ng) * tau (3, na) )
              fac = omega * vloc (igtongl (ng), ityp (na) ) * tpiba2 * &
                   (DREAL (rhog (nl (ng) ) ) * cos (gtau) - &
                    DIMAG (rhog (nl (ng) ) ) * sin (gtau) )
              dynwrk (na_icart, na_jcart) = dynwrk (na_icart, na_jcart) - &
                   fac * g (icart, ng) * g (jcart, ng)
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  call reduce (18 * nat * nat, dynwrk)
  !
  ! each pool contributes to next term
  !
100 continue
#endif
  !      goto 500
  !
  ! Here we compute  the nonlocal Ultra-soft contribution
  !
  if (nksq.gt.1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
     else
        ikk = 2 * ik - 1
     endif
     if (lsda) current_spin = isk (ikk)
     if (nksq.gt.1) read (iunigk) npw, igk
     ! npwq and igkq are not actually used

     if (nksq.gt.1.and..not.lgamma) read (iunigk) npwq, igkq

     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ikk, - 1)
     call init_us_2 (npw, igk, xk (1, ikk), vkb)
     !
     !    We first prepare the gamma terms, which are the second derivatives
     !    becp terms.
     !
     do icart = 1, 3
        do jcart = 1, icart
           do ibnd = 1, nbnd
              do ig = 1, npw
                 aux1 (ig, ibnd) = - evc (ig, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk (ig) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk (ig) ) )
              enddo
           enddo

           call ccalbec(nkb,npwx,npw,nbnd,gammap(1,1,icart,jcart),vkb,aux1)
           if (jcart.lt.icart) &
               call ZCOPY (nkb * nbnd, gammap (1, 1, icart, jcart), 1, &
                                       gammap (1, 1, jcart, icart), 1)
        enddo
     enddo
     !
     !   And then compute the contribution from the US pseudopotential
     !   which is  similar to the KB one
     !
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do icart = 1, 3
                 na_icart = 3 * (na - 1) + icart
                 do jcart = 1, 3
                    na_jcart = 3 * (na - 1) + jcart
                    do ibnd = 1, nbnd_occ (ikk)
                       wgg = wg (ibnd, ikk)
                       do ih = 1, nh (nt)
                          ikb = ijkb0 + ih
                          do jh = 1, nh (nt)
                             jkb = ijkb0 + jh
                             dynwrk(na_icart,na_jcart) = &
                                  dynwrk(na_icart,na_jcart) + &
                                  (deeq (ih, jh, na, current_spin) - &
                                   et (ibnd, ikk) * qq (ih,jh,nt) ) * wgg * &
                                  (conjg(gammap(ikb, ibnd, icart, jcart)) *&
                                   becp1 (jkb, ibnd, ik) + &
                                   conjg (becp1 (ikb, ibnd, ik) ) * &
                                   gammap (jkb, ibnd, icart, jcart) + &
                                   conjg (alphap (ikb, ibnd, icart, ik) ) * &
                                   alphap (jkb, ibnd, jcart, ik) + &
                                   conjg (alphap (ikb, ibnd, jcart, ik) ) * &
                                   alphap (jkb, ibnd, icart, ik) )
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     enddo
  enddo
  !
  !   For true US pseudopotentials there is an additional term in the second
  !   derivative which is due to the change of the self consistent D part
  !   when the atom moves. We compute these terms in an additional routine
  !
  call addusdynmat (dynwrk)
#ifdef __PARA
  call poolreduce (18 * nat * nat, dynwrk)
#endif
  !      do na = 1,nat
  !         do nb = 1,nat
  !           WRITE( stdout, '(2i3)') na,nb
  !            do icart = 1,3
  !              na_icart = 3*(na-1)+icart
  !               WRITE( stdout,'(6f13.8)')
  !     +               (dynwrk(na_icart,3*(nb-1)+jcart), jcart=1,3)
  !            end do
  !         end do
  !      end do
  !      call stop_ph(.false.)
  !
  !  We rotate the dynamical matrix on the basis of patterns
  !
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        work = (0.0d0, 0.0d0)
        do na_jcart = 1, 3 * nat
           do na_icart = 1, 3 * nat
              work = work + conjg (u (na_icart, nu_i) ) * &
                            dynwrk (na_icart, na_jcart) * &
                            u (na_jcart, nu_j)
           enddo
        enddo
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) + work
     enddo

  enddo
  deallocate (gammap)
  deallocate (aux1)
  deallocate (work2)
  deallocate (work1)
  deallocate (rhog)

  call stop_clock ('dynmat_us')
  return
end subroutine dynmat_us
