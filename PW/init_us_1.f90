!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_1
  !----------------------------------------------------------------------
  !
  !   This routine performs the following tasks:
  !   a) For each non vanderbilt pseudopotential it computes the D and
  !      the betar in the same form of the Vanderbilt pseudopotential.
  !   b) It computes the indices indv which establish the correspondence
  !      nh <-> beta in the atom
  !   c) It computes the indices nhtol which establish the correspondence
  !      nh <-> angular momentum of the beta function
  !   d) It computes the indices nhtom which establish the correspondence
  !      nh <-> magnetic angular momentum of the beta function.
  !   e) It computes the coefficients c_{LM}^{nm} which relates the
  !      spherical harmonics in the Q expansion
  !   f) It computes the radial fourier transform of the Q function on
  !      all the g vectors
  !   g) It computes the q terms which define the S matrix.
  !   h) It fills the interpolation table for the beta functions
  !
#include "machine.h"
  USE kinds, ONLY: DP
  USE parameters, ONLY: lmaxx, ndm, nbrx, lqmax
  USE constants, ONLY: fpi
  USE atom, ONLY: r, rab, vnl, chi, msh
  USE basis, ONLY: ntyp
  USE cell_base, ONLY: omega, tpiba
  USE gvect, ONLY: g, gg
  USE pseud, ONLY: lloc, lmax
  USE us, ONLY: lqx, dion, betar, qfunc, qfcoef, rinner, nh, nbeta, &
       kkbeta, nqf, nqlc, lll, tvanp, okvan, newpseudo, lmaxkb, mx, &
       nlx, ap, lpx, lpl, nqxq, dq, qgm, nqx, tab, dvan, qq, nhtol, &
       nhtom, qrad, indv, nhm
  implicit none
  !
  !     here a few local variables
  !

  integer :: nt, ih, jh, nb, mb, nmb, l, m, ir, iq, is, startq, &
       lastq, ilast
  ! various counters
  real(kind=DP), allocatable :: aux (:), aux1 (:), besr (:), qtot (:,:,:)
  ! various work space
  real(kind=DP) :: prefr, pref, q, qi
  ! the prefactor of the q functions
  ! the prefactor of the beta functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(kind=DP), allocatable :: ylmk0 (:)
  ! the spherical harmonics
  real(kind=DP) ::  vll (0:lmaxx),vqint
  ! the denominator in KB case
  ! interpolated value

  call start_clock ('init_us_1')
  !
  !    Initialization of the variables
  !
  allocate (aux ( ndm))    
  allocate (aux1( ndm))    
  allocate (besr( ndm))    
  allocate (qtot( ndm , nbrx , nbrx))    
  allocate (ylmk0( lqx * lqx))    
  dvan (:,:,:) = 0.d0
  qq (:,:,:)   = 0.d0
  ap (:,:,:)   = 0.d0
  if (lqx > 0) qrad(:,:,:,:)= 0.d0

  prefr = fpi / omega
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtom,
  !   indv, and if the pseudopotential is of KB type we initialize the
  !   atomic D terms
  !
  do nt = 1, ntyp
     if (tvanp (nt) .or.newpseudo (nt) ) then
        ih = 1
        do nb = 1, nbeta (nt)
           l = lll (nb, nt)
           do m = 1, 2 * l + 1
              nhtol (ih, nt) = l
              nhtom (ih, nt) = m
              indv (ih, nt) = nb
              ih = ih + 1
           enddo
        enddo
     else
        !
        !     In this case we simulate a KB pseudopotential as a US one
        !     We set the important dimensions,
        !
        nbeta (nt) = lmax (nt)
        kkbeta (nt) = msh (nt)
        !
        !     And we compute the equivalent of D which is the energy parameter
        !
        do l = 0, lmax (nt)
           if (l.ne.lloc (nt) ) then
              !
              !      simpson integration of the chi functions
              !
              do ir = 1, msh (nt)
                 aux (ir) = chi (ir, l + 1, nt) **2 * vnl (ir, l, nt)
              enddo
              call simpson (msh (nt), aux, rab (1, nt), vll (l) )
           endif

        enddo
        !
        !    Here we set the D functions and the beta functions
        !
        dion (:,:,nt) = 0.d0
        nb = 1
        ih = 1
        do l = 0, lmax (nt)
           if (l.ne.lloc (nt) ) then
              dion (nb, nb, nt) = 1.d0 / vll (l)
              betar (0, nb, nt) = 0.d0
              do ir = 1, kkbeta (nt)
                 betar (ir, nb, nt) = vnl (ir, l, nt) * chi (ir, l + 1, nt)
              enddo
              lll (nb, nt) = l
              do m = 1, 2 * l + 1
                 nhtol (ih, nt) = l
                 nhtom (ih, nt) = m
                 indv (ih, nt) = nb
                 ih = ih + 1
              enddo
              nb = nb + 1
           endif
        enddo
     endif
     !
     !    From now on the only difference between KB and US pseudopotentials
     !    is in the presence of the q and Q functions.
     !
     !    Here we initialize the D of the solid
     !
     do ih = 1, nh (nt)
        do jh = 1, nh (nt)
           if (nhtol (ih, nt) .eq.nhtol (jh, nt) .and.nhtom (ih, nt).eq.nhtom (jh, nt) ) then
              ir = indv (ih, nt)
              is = indv (jh, nt)
              dvan (ih, jh, nt) = dion (ir, is, nt)
           endif
        enddo
     enddo

  enddo
  !
  !  compute Clebsch-Gordan coefficients
  !
  if (okvan) call aainit (lmaxkb + 1, lqmax, mx, nlx, ap, lpx,lpl)
  !
  !   here for the US types we compute the Fourier transform of the
  !   Q functions.
  !
  call divide (nqxq, startq, lastq)
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do l = 0, nqlc (nt) - 1
           !
           !     first we build for each nb,mb,l the total Q(|r|) function
           !     note that l is the true angular momentum, and the arrays
           !     have dimensions 1..l+1
           !
           do nb = 1, nbeta (nt)
              do mb = nb, nbeta (nt)
                 if ( (l.ge.abs (lll (nb, nt) - lll (mb, nt) ) ) .and. &
                      (l.le.lll (nb, nt) + lll (mb, nt) )        .and. &
                      (mod (l + lll (nb, nt) + lll (mb, nt), 2) .eq.0) ) then
                    do ir = 1, kkbeta (nt)
                       if (r (ir, nt) .ge.rinner (l + 1, nt) ) then
                          qtot (ir, nb, mb) = qfunc (ir, nb, mb, nt)
                       else
                          ilast = ir
                       endif
                    enddo
                    if (rinner (l + 1, nt) .gt.0.d0) &
                         call setqf(qfcoef (1, l+1, nb, mb, nt), &
                                    qtot(1,nb,mb), r(1,nt), nqf(nt),l,ilast)
                 endif
              enddo
           enddo
           !
           !     here we compute the spherical bessel function for each |g|
           !
           do iq = startq, lastq
              q = (iq - 1) * dq * tpiba
              call sph_bes (kkbeta (nt), r (1, nt), q, l, aux)
              !
              !   and then we integrate with all the Q functions
              !
              do nb = 1, nbeta (nt)
                 !
                 !    the Q are symmetric with respect to indices
                 !
                 do mb = nb, nbeta (nt)
                    nmb = mb * (mb - 1) / 2 + nb
                    if ( (l.ge.abs (lll (nb, nt) - lll (mb, nt) ) ) .and. &
                         (l.le.lll (nb, nt) + lll (mb, nt) )        .and. &
                         (mod (l + lll(nb, nt) + lll(mb, nt), 2) .eq.0) ) then
                       do ir = 1, kkbeta (nt)
                          aux1 (ir) = aux (ir) * qtot (ir, nb, mb)
                       enddo
                       call simpson (kkbeta(nt), aux1, rab(1, nt), &
                                     qrad(iq,nmb,l + 1, nt) )
                    endif
                 enddo
              enddo
              ! igl
           enddo
           ! l
        enddo
        call DSCAL (nqxq * nbrx * (nbrx + 1) / 2 * lqx, prefr, &
                    qrad (1, 1, 1, nt), 1)
#ifdef __PARA
        call reduce (nqxq * nbrx * (nbrx + 1) / 2 * lqx, qrad (1, 1, 1, nt) )
#endif
     endif
     ! ntyp

  enddo
  !
  !   and finally we compute the qq coefficients by integrating the Q.
  !   q are the g=0 components of Q.
  !
#ifdef __PARA
  if (gg (1) .gt.1.0d-8) goto 100
#endif
  call ylmr2 (lqx * lqx, 1, g, gg, ylmk0)
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (1, ih, jh, nt, gg, qgm, ylmk0)
              qq (ih, jh, nt) = omega * DREAL (qgm (1) )
              qq (jh, ih, nt) = qq (ih, jh, nt)
              !                  WRITE( stdout,*) ih,jh,nt,qq(ih,jh,nt)
           enddo
        enddo
     endif
  enddo
#ifdef __PARA
100 continue
  call reduce (nhm * nhm * ntyp, qq)
#endif
  !
  !     fill the interpolation table tab
  !
  pref = fpi / sqrt (omega)
  call divide (nqx, startq, lastq)
  tab (:,:,:) = 0.d0
  do nt = 1, ntyp
     do nb = 1, nbeta (nt)
        l = lll (nb, nt)
        do iq = startq, lastq
           qi = (iq - 1) * dq
           call sph_bes (kkbeta (nt), r (1, nt), qi, l, besr)
           do ir = 1, kkbeta (nt)
              aux (ir) = betar (ir, nb, nt) * besr (ir) * r (ir, nt)
           enddo
           call simpson (kkbeta (nt), aux, rab (1, nt), vqint)
           tab (iq, nb, nt) = vqint * pref
        enddo
     enddo
  enddo
#ifdef __PARA
  call reduce (nqx * nbrx * ntyp, tab)
#endif
  deallocate (ylmk0)
  deallocate (qtot)
  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)

  call stop_clock ('init_us_1')
  return
end subroutine init_us_1

