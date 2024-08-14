!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_1( nat, ityp, omega, qmax, intra_bgrp_comm )
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
  !   f) It computes the interpolation table "qrad" for Q(G)
  !   g) It computes the qq terms which define the S matrix.
  !
  USE upf_kinds,    ONLY : DP
  USE upf_const,    ONLY : fpi, sqrt2
  USE uspp,         ONLY : nhtol, nhtoj, nhtolm, ijtoh, dvan, qq_at, qq_nt, indv, &
                           ap, aainit, qq_so, dvan_so, okvan, ofsbeta
  USE uspp_param,   ONLY : upf, lmaxq, nh, nhm, lmaxkb, nsp
  USE upf_spinorb,  ONLY : is_spinorbit, rot_ylm, fcoef, lmaxx, &
                           transform_qq_so
  USE qrad_mod,     ONLY : init_tab_qrad
  USE paw_variables,ONLY : okpaw
  USE mp,           ONLY : mp_sum
  implicit none
  !
  integer,  intent(in) :: nat
  integer,  intent(in) :: ityp(nat)
  real(DP), intent(in) :: qmax
  real(DP), intent(in) :: omega
  integer,  intent(in) :: intra_bgrp_comm
  !
  !     here a few local variables
  !
  integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, iq, is, ia
  ! various counters
  real(DP) ::  j
  ! J=L+S (noninteger!)
  integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
             lk, mk, kh, lh, ijkb0, na
  integer, external :: sph_ind
  complex(DP) :: coeff
  real(DP) :: ji, jk
  real(DP), EXTERNAL :: spinor
  !
  call start_clock ('init_us_1')
  !
  !    NB: duplicated modules' variables are syncronized at the end. This
  !        may lead to problems if these variables are using during function
  !        calls in this subroutines. However this should never happen.
  !
  !    Initialization of the variables
  !
  ap (:,:,:)   = 0.d0
  !
  ! the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
  ! but in some versions of the PP files lmax is not set to the maximum
  ! l of the beta functions but includes the l of the local potential
  !
  do nt=1,nsp
     upf(nt)%nqlc = MIN ( upf(nt)%nqlc, lmaxq )
     IF ( upf(nt)%nqlc < 0 )  upf(nt)%nqlc = 0
  end do

  if (is_spinorbit) then
     !
     !  In the spin-orbit case we need the unitary matrix u which rotates the
     !  real spherical harmonics and yields the complex ones.
     !
     rot_ylm=(0.d0,0.d0)
     l=lmaxx
     rot_ylm(l+1,1)=(1.d0,0.d0)
     do n1=2,2*l+1,2
       m=n1/2
       n=l+1-m
       rot_ylm(n,n1)=CMPLX((-1.d0)**m/sqrt2,0.0_dp,kind=DP)
       rot_ylm(n,n1+1)=CMPLX(0.d0,-(-1.d0)**m/sqrt2,kind=DP)
       n=l+1+m
       rot_ylm(n,n1)=CMPLX(1.0_dp/sqrt2,0.d0,kind=DP)
       rot_ylm(n,n1+1)=CMPLX(0.d0, 1.0_dp/sqrt2,kind=DP)
     enddo
  endif
  if ( nhm > 0 ) then
     if (is_spinorbit) then
        fcoef=(0.d0,0.d0)
        dvan_so = (0.d0,0.d0)
        qq_so=(0.d0,0.d0)
     else
        dvan = 0.d0
     end if
     qq_nt = 0.d0
     qq_at = 0.d0
  endif
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtolm,
  !   nhtoj, indv, and if the pseudopotential is of KB type we initialize the
  !   atomic D terms
  !
  ijkb0 = 0
  do nt = 1, nsp
     ih = 1
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        do m = 1, 2 * l + 1
           nhtol (ih, nt) = l
           nhtolm(ih, nt) = l*l+m
           indv  (ih, nt) = nb
           ih = ih + 1
        enddo
     enddo
     if ( upf(nt)%has_so ) then
        ih = 1
        do nb = 1, upf(nt)%nbeta
           l = upf(nt)%lll (nb)
           j = upf(nt)%jjj (nb)
           do m = 1, 2 * l + 1
              nhtoj (ih, nt) = j
              ih = ih + 1
           enddo
        enddo
     endif
     !
     ! ijtoh map augmentation channel indexes ih and jh to composite
     ! "triangular" index ijh
     if ( nhm > 0 ) ijtoh(:,:,nt) = -1
     ijv = 0
     do ih = 1,nh(nt)
         do jh = ih,nh(nt)
             ijv = ijv+1
             ijtoh(ih,jh,nt) = ijv
             ijtoh(jh,ih,nt) = ijv
         enddo
     enddo
     !
     ! ijkb0 points to the last beta "in the solid" for atom ia-1
     ! i.e. ijkb0+1,.. ijkb0+nh(ityp(ia)) are the nh beta functions of
     !      atom ia in the global list of beta functions (ijkb0=0 for ia=1)
     do ia = 1,nat
       IF ( ityp(ia) == nt ) THEN
          ofsbeta(ia) = ijkb0
          ijkb0 = ijkb0 + nh(nt)
        END IF
     enddo
     !
     !    From now on the only difference between KB and US pseudopotentials
     !    is in the presence of the qq and Q functions.
     !
     !    Here we initialize the D of the solid
     !
     if (upf(nt)%has_so) then
       !
       !  first calculate the fcoef coefficients
       !
       do ih = 1, nh (nt)
          li = nhtol(ih, nt)
          ji = nhtoj(ih, nt)
          mi = nhtolm(ih, nt)-li*li
          do kh=1,nh(nt)
            lk = nhtol(kh, nt)
            jk = nhtoj(kh, nt)
            mk = nhtolm(kh, nt)-lk*lk
            if (li == lk .and. abs(ji-jk) < 1.d-7) then
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
                    dvan_so(ih,jh,ijs,nt) = upf(nt)%dion(vi,vj) * &
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
              if (is_spinorbit) then
                 dvan_so (ih, jh, 1, nt) = upf(nt)%dion (ir, is)
                 dvan_so (ih, jh, 4, nt) = upf(nt)%dion (ir, is)
              else
                 dvan (ih, jh, nt) = upf(nt)%dion (ir, is)
              endif
            endif
          enddo
        enddo
     endif
  enddo
  !
  !  compute Clebsch-Gordan coefficients
  !
  if (okvan .or. okpaw) call aainit (lmaxkb + 1)
  !
  !   here for the US types we compute the Fourier transform of the
  !   Q functions.
  !
  IF ( lmaxq > 0 ) CALL init_tab_qrad(qmax, omega, intra_bgrp_comm, ir)
  !
  !   and finally we compute the qq coefficients by integrating the Q.
  !   The qq are the g=0 components of Q
  !
  call compute_qqr ( 1.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], omega, qq_nt )
  if ( is_spinorbit ) call transform_qq_so( qq_nt, qq_so )
  !
  ! finally we set the atomic specific qq_at matrices
  if ( nhm > 0 ) then
     do na=1, nat
        qq_at(:,:, na) = qq_nt(:,:,ityp(na))
     end do
  end if
  !
  ! update GPU memory (taking care of zero-dim allocations)
  !
  if (nhm>0) then
     !$acc update device(indv)
     !$acc update device(nhtol)
     !$acc update device(nhtoj)
     !$acc update device(ijtoh)
     !$acc update device(qq_at)
     if (is_spinorbit) then
        !$acc update device(dvan_so)
        !$acc update device(fcoef)
        !$acc update device(qq_so)
     else 
        !$acc update device(dvan)
     endif
  endif
  !
  call stop_clock ('init_us_1')
  return
  !
end subroutine init_us_1
