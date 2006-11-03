!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  !   d) It computes the indices nhtolm which establish the correspondence
  !      nh <-> combined (l,m) index for the beta function.
  !   e) It computes the coefficients c_{LM}^{nm} which relates the
  !      spherical harmonics in the Q expansion
  !   f) It computes the radial fourier transform of the Q function on
  !      all the g vectors
  !   g) It computes the q terms which define the S matrix.
  !   h) It fills the interpolation table for the beta functions
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lmaxx
  USE constants,  ONLY : fpi
  USE atom,       ONLY : r, rab
  USE ions_base,  ONLY : ntyp => nsp
  USE cell_base,  ONLY : omega, tpiba
  USE gvect,      ONLY : g, gg
  USE lsda_mod,   ONLY : nspin
#ifdef USE_SPLINES
  USE us,         ONLY : nqxq, dq, nqx, tab, tab_d2y, qrad
  USE splinelib
#else
  USE us,         ONLY : nqxq, dq, nqx, tab, qrad
#endif
  USE uspp,       ONLY : nhtol, nhtoj, nhtolm, dvan, qq, indv, ap, aainit, &
                         qq_so, dvan_so, okvan
  USE uspp_param, ONLY : lmaxq, dion, betar, qfunc, qfcoef, rinner, nbeta, &
                         kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, &
                         nbetam, nhm
  USE spin_orb,   ONLY : lspinorb, so, rot_ylm, fcoef
  !
  implicit none
  !
  !     here a few local variables
  !

  integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, iq, is, startq, &
       lastq, ilast, ndm
  ! various counters
  real(DP), allocatable :: aux (:), aux1 (:), besr (:), qtot (:,:)
  ! various work space
  real(DP) :: prefr, pref, q, qi
  ! the prefactor of the q functions
  ! the prefactor of the beta functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP), allocatable :: ylmk0 (:)
  ! the spherical harmonics
  real(DP) ::  vll (0:lmaxx), vqint, sqrt2, j
  ! the denominator in KB case
  ! interpolated value
  integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
             lk, mk, vk, kh, lh
  integer, external :: sph_ind
  complex(DP) :: coeff, qgm(1)
  real(DP) :: spinor, ji, jk
#ifdef USE_SPLINES
  real(DP), allocatable :: xdata(:)
  real(DP) :: d1
#endif

  call start_clock ('init_us_1')
  !
  !    Initialization of the variables
  !
  ndm = MAXVAL (kkbeta(1:ntyp))
  allocate (aux ( ndm))    
  allocate (aux1( ndm))    
  allocate (besr( ndm))    
  allocate (qtot( ndm , nbetam*(nbetam+1)/2 ))    
  allocate (ylmk0( lmaxq * lmaxq))    
  ap (:,:,:)   = 0.d0
  if (lmaxq > 0) qrad(:,:,:,:)= 0.d0
  !
  ! the following prevents an out-of-bound error: nqlc(is)=2*lmax+1
  ! but in some versions of the PP files lmax is not set to the maximum
  ! l of the beta functions but includes the l of the local potential
  !
  do nt=1,ntyp
     nqlc(nt) = MIN ( nqlc(nt), lmaxq )
     IF ( nqlc(nt) < 0 )  nqlc(nt) = 0
  end do

  prefr = fpi / omega
  if (lspinorb) then
!
!  In the spin-orbit case we need the unitary matrix u which rotates the
!  real spherical harmonics and yields the complex ones.
!
     sqrt2=1.d0/dsqrt(2.d0)
     rot_ylm=(0.d0,0.d0)
     l=lmaxx
     rot_ylm(l+1,1)=(1.d0,0.d0)
     do n1=2,2*l+1,2
       m=n1/2
       n=l+1-m
       rot_ylm(n,n1)=CMPLX((-1.d0)**m*sqrt2,0.d0)
       rot_ylm(n,n1+1)=CMPLX(0.d0,-(-1.d0)**m*sqrt2)
       n=l+1+m
       rot_ylm(n,n1)=CMPLX(sqrt2,0.d0)
       rot_ylm(n,n1+1)=CMPLX(0.d0, sqrt2)
     enddo
     fcoef=(0.d0,0.d0)
     dvan_so = (0.d0,0.d0)
     qq_so=(0.d0,0.d0)
     qq  = 0.d0
  else
     qq  = 0.d0
     dvan = 0.d0
  endif
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtolm,
  !   nhtoj, indv, and if the pseudopotential is of KB type we initialize the
  !   atomic D terms
  !
  do nt = 1, ntyp
     ih = 1
     do nb = 1, nbeta (nt)
        l = lll (nb, nt)
        j = jjj (nb, nt)
        do m = 1, 2 * l + 1
           nhtol (ih, nt) = l
           nhtolm(ih, nt) = l*l+m
           nhtoj (ih, nt) = j
           indv  (ih, nt) = nb
           ih = ih + 1
        enddo
     enddo
     !
     !    From now on the only difference between KB and US pseudopotentials
     !    is in the presence of the q and Q functions.
     !
     !    Here we initialize the D of the solid
     !
     if (so(nt)) then
     !
     !  first calculate the fcoef coefficients
     !
       do ih = 1, nh (nt)
          li = nhtol(ih, nt)
          ji = nhtoj(ih, nt)
          mi = nhtolm(ih, nt)-li*li
          vi = indv (ih, nt)
          do kh=1,nh(nt)
            lk = nhtol(kh, nt)
            jk = nhtoj(kh, nt)
            mk = nhtolm(kh, nt)-lk*lk
            vk = indv (kh, nt)
            if (li.eq.lk.and.abs(ji-jk).lt.1.d-7) then
              do is1=1,2
                do is2=1,2
                  coeff = (0.d0, 0.d0)
                  do m=-li-1, li
                    m0= sph_ind(li,ji,m,is1) + lmaxx + 1
                    m1= sph_ind(lk,jk,m,is2) + lmaxx + 1
                    coeff=coeff + rot_ylm(m0,mi)*spinor(li,ji,m,is1)* &
                            CONJG(rot_ylm(m1,mk))*spinor(lk,jk,m,is2)
                  enddo
                  fcoef(ih,kh,is1,is2,nt)=coeff
                enddo
              enddo
            endif
          enddo
        enddo
!
!   and calculate the bare coefficients
!
        do ih = 1, nh (nt)
           vi = indv (ih, nt)
           do jh = 1, nh (nt)
              vj = indv (jh, nt)
              ijs=0
              do is1=1,2
                 do is2=1,2
                    ijs=ijs+1
                    dvan_so(ih,jh,ijs,nt) = dion(vi,vj,nt) * &
                                            fcoef(ih,jh,is1,is2,nt)
                    if (vi.ne.vj) fcoef(ih,jh,is1,is2,nt)=(0.d0,0.d0)
                 enddo
              enddo
           enddo
        enddo
     else
        do ih = 1, nh (nt)
          do jh = 1, nh (nt)
            if (nhtol (ih, nt) == nhtol (jh, nt) .and. &
              nhtolm(ih, nt) == nhtolm(jh, nt) ) then
              ir = indv (ih, nt)
              is = indv (jh, nt)
              if (lspinorb) then
                 dvan_so (ih, jh, 1, nt) = dion (ir, is, nt)
                 dvan_so (ih, jh, 4, nt) = dion (ir, is, nt)
              else
                 dvan (ih, jh, nt) = dion (ir, is, nt)
              endif
            endif
          enddo
        enddo
     endif
  enddo
  !
  !  compute Clebsch-Gordan coefficients
  !
  if (okvan) call aainit (lmaxkb + 1)
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
           !     note that l is the true (combined) angular momentum
           !     and that the arrays have dimensions 1..l+1
           !
           do nb = 1, nbeta (nt)
              do mb = nb, nbeta (nt)
                 ijv = mb * (mb-1) / 2 + nb
                 if ( (l >= abs (lll (nb, nt) - lll (mb, nt) ) ) .and. &
                      (l <= lll (nb, nt) + lll (mb, nt) )        .and. &
                      (mod (l + lll (nb, nt) + lll (mb, nt), 2) == 0) ) then
                    do ir = 1, kkbeta (nt)
                       if (r (ir, nt) >= rinner (l + 1, nt) ) then
                          qtot (ir, ijv) = qfunc (ir, ijv, nt)
                       else
                          ilast = ir
                       endif
                    enddo
                   if (rinner (l + 1, nt) > 0.d0) &
                         call setqf(qfcoef (1, l+1, nb, mb, nt), &
                                    qtot(1,ijv), r(1,nt), nqf(nt),l,ilast)
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
                    ijv = mb * (mb - 1) / 2 + nb
                    if ( (l >= abs (lll (nb, nt) - lll (mb, nt) ) ) .and. &
                         (l <= lll (nb, nt) + lll (mb, nt) )        .and. &
                         (mod (l + lll(nb, nt) + lll(mb, nt), 2) == 0) ) then
                       do ir = 1, kkbeta (nt)
                          aux1 (ir) = aux (ir) * qtot (ir, ijv)
                       enddo
                       call simpson (kkbeta(nt), aux1, rab(1, nt), &
                                     qrad(iq,ijv,l + 1, nt) )
                    endif
                 enddo
              enddo
              ! igl
           enddo
           ! l
        enddo
        qrad (:, :, :, nt) = qrad (:, :, :, nt)*prefr
#ifdef __PARA
        call reduce (nqxq * nbetam*(nbetam+1) / 2 * lmaxq, qrad (1, 1, 1, nt) )
#endif
     endif
     ! ntyp

  enddo
  !
  !   and finally we compute the qq coefficients by integrating the Q.
  !   q are the g=0 components of Q.
  !
#ifdef __PARA
  if (gg (1) > 1.0d-8) goto 100
#endif
  call ylmr2 (lmaxq * lmaxq, 1, g, gg, ylmk0)
  do nt = 1, ntyp
    if (tvanp (nt) ) then
      if (so(nt)) then
        do ih=1,nh(nt)
          do jh=1,nh(nt)
            call qvan2 (1, ih, jh, nt, gg, qgm, ylmk0)
            qq (ih, jh, nt) = omega *  DBLE (qgm (1) )
            do kh=1,nh(nt)
              do lh=1,nh(nt)
                ijs=0
                do is1=1,2
                  do is2=1,2
                    ijs=ijs+1
                    do is=1,2
                      qq_so(kh,lh,ijs,nt) = qq_so(kh,lh,ijs,nt)       &
                          + omega* DBLE(qgm(1))*fcoef(kh,ih,is1,is,nt)&
                                               *fcoef(jh,lh,is,is2,nt)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      else
        do ih = 1, nh (nt)
          do jh = ih, nh (nt)
             call qvan2 (1, ih, jh, nt, gg, qgm, ylmk0)
             if (lspinorb) then
                 qq_so (ih, jh, 1, nt) = omega *  DBLE (qgm (1) )
                 qq_so (jh, ih, 1, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (ih, jh, 4, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (jh, ih, 4, nt) = qq_so (ih, jh, 4, nt)
             endif
             qq (ih, jh, nt) = omega *  DBLE (qgm (1) )
             qq (jh, ih, nt) = qq (ih, jh, nt)
          enddo
        enddo
      endif
    endif
  enddo
#ifdef __PARA
100 continue
  if (lspinorb) then
    call reduce ( nhm * nhm * ntyp * 8, qq_so )
    call reduce ( nhm * nhm * ntyp, qq )
  else
    call reduce ( nhm * nhm * ntyp, qq )
  endif
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
  call reduce (nqx * nbetam * ntyp, tab)
#endif

#ifdef USE_SPLINES
  ! initialize spline interpolation
  allocate(xdata(lastq-startq+1))
  do iq = startq, lastq
    xdata(iq) = (iq - 1) * dq
  enddo
  do nt = 1, ntyp
     do nb = 1, nbeta (nt)
        l = lll (nb, nt)
        d1 = (tab(2,nb,nt) - tab(1,nb,nt)) / dq
        call spline(xdata, tab(:,nb,nt), 0.d0, d1, tab_d2y(:,nb,nt))
     enddo
  enddo
  deallocate(xdata)
#endif

  deallocate (ylmk0)
  deallocate (qtot)
  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)

  call stop_clock ('init_us_1')
  return
end subroutine init_us_1

