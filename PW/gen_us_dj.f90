!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine gen_us_dj (ik, dvkb)  
  !----------------------------------------------------------------------
  !
  !   Calculates the beta function pseudopotentials with
  !   the derivative of the Bessel functions
  !
#include "machine.h"
  use pwcom  
  use allocate
  implicit none
  !
  integer :: ik
  complex(kind=DP) :: dvkb (npwx, nkb)  
  !
  ! local variables
  !
  integer :: ikb, nb, ih, ig, i0, nt  
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on G vectors
  ! index of the first nonzero point in the r
  ! counter on atomic type

  real(kind=DP) :: arg  
  ! argument of the atomic phase factor

  complex(kind=DP) :: phase, pref  
  ! atomic phase factor
  ! prefactor

  integer :: na, i, m, l, iig, lm  
  real(kind=DP), pointer :: djl (:,:,:), ylm (:,:), q (:), gk (:,:)
  real(kind=DP) ::  jl (ndm), jlm1 (ndm), qt, dv, eps
  parameter (eps = 1.0e-8)  

  complex(kind=DP), pointer :: sk (:)  

  if (nkb.eq.0) return  
  call mallocate(djl, npw , nbrx , ntyp)  
  call mallocate(ylm, npw ,(lmaxkb + 1) **2)  
  call mallocate(gk, 3, npw)  
  call mallocate(q, npw)  
  do ig = 1, npw  
     gk (1,ig) = xk (1, ik) + g(1, igk(ig) )
     gk (2,ig) = xk (2, ik) + g(2, igk(ig) )
     gk (3,ig) = xk (3, ik) + g(3, igk(ig) )
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  call ylmr2 ((lmaxkb+1)**2, npw, gk, q, ylm)  

  do nt = 1, ntyp  
     do nb = 1, nbeta (nt)  
        l = lll (nb, nt)  
        do ig = 1, npw  
           qt = sqrt(q (ig)) * tpiba  
           if (qt.lt.eps) then  
              if (l.ne.1) then  
                 do i = 1, kkbeta (nt)  
                    jl (i) = 0.0d0  
                 enddo
              else  
                 ! Note that dj_1/dx (x=0) = 1/3
                 do i = 1, kkbeta (nt)  
                    jl (i) = 1.0d0 / 3.d0  
                 enddo
              endif
           else  
              !
              !    in order to avoid a division by zero i0 is defined as
              !    the first point in the radial mesh such that q*r>eps
              !    NB: eps value must be consistent with its value in sph_bes
              !
              i0 = 1  
              do while ( qt * r(i0,nt) .lt. eps )
                 i0 = i0 + 1
              end do
              call sph_bes (kkbeta(nt)+1-i0, r(i0, nt), qt, l, jl (i0))
              call sph_bes (kkbeta(nt)+1-i0, r(i0, nt), qt, l - 1, jlm1(i0) )
              ! recurrence relation for jl
              do i = i0, kkbeta (nt)  
                 jl (i) = jlm1 (i) - (l + 1) / (qt * r (i, nt) ) * jl (i)  
              enddo
              if (i0.eq.2) jl (1) = jl (2)  
           endif
           ! jl is now the derivative of the Bessel functions
           do i = 1, kkbeta (nt)  
              jlm1 (i) = jl (i) * betar (i, nb, nt) * r (i, nt) **2  
           enddo
           call simpson (kkbeta (nt), jlm1, rab (1, nt), dv)  
           djl (ig, nb, nt) = dv * fpi / sqrt (omega)  
        enddo
     enddo
  enddo

  call mfree (q)  
  call mfree (gk)  

  call mallocate(sk, npw)  
  ikb = 0  
  do nt = 1, ntyp  
     do na = 1, nat  
        if (ityp (na) .eq.nt) then  
           arg = (xk (1, ik) * tau(1,na) + &
                  xk (2, ik) * tau(2,na) + &
                  xk (3, ik) * tau(3,na) ) * tpi
           phase = DCMPLX (cos (arg), - sin (arg) )  
           do ig = 1, npw  
              iig = igk (ig)  
              sk (ig) = eigts1 (ig1 (iig), na) * &
                        eigts2 (ig2 (iig), na) * &
                        eigts3 (ig3 (iig), na) * phase
           enddo
           do ih = 1, nh (nt)  
              nb = indv (ih, nt)  
              l = nhtol (ih, nt)  
              m = nhtom (ih, nt)  
              lm = l * l + m  
              ikb = ikb + 1  
              pref = (0.d0, - 1.d0) **l  
              !
              do ig = 1, npw  
                 dvkb (ig, ikb) = djl (ig, nb, nt) * sk (ig) * ylm (ig, lm) &
                      * pref
              enddo
           enddo
        endif
     enddo

  enddo

  if (ikb.ne.nkb) call error ('gen_us_dj', 'unexpected error', 1)  
  call mfree (sk)  
  call mfree (ylm)  
  call mfree (djl)  
  return  
end subroutine gen_us_dj

