!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*******************************************************************

subroutine take_nloc_k_kq(kpoint, kpoint_q, phi, jpol, duk)

!*******************************************************************

  use ions_base, only: nat, ntyp => nsp, ityp
  use cell_base, only: tpiba, at, alat
  use lsda_mod, only: isk, current_spin, lsda
  use kinds, only: DP
  use uspp_param, only: nh
  use uspp, only: nkb, vkb, deeq
  use wvfct, only:npwx , npw, igk, nbnd
  USE wavefunctions_module,  ONLY: evc
  use klist, only: xk
  use gvect, only: g
  use nmr_mod, only: igk_q, vkb_q
  use phcom , only : igkq
  use control_ph,       ONLY : nbnd_occ
  use becmod , only : becp
  implicit none

  integer, intent(in) :: kpoint,kpoint_q
  real(kind=DP) :: ej(3)
  complex(kind=DP), intent(in) :: phi (npwx,nbnd)
  complex(kind=DP), intent(out) :: duk (npwx,nbnd)

  integer :: i, jpol,  nt, na, ikb, jkb, ibnd, jbnd, ig, ijkb0, ih, jh

  real(kind=DP), allocatable  :: gk (:,:), gk_q (:,:), g2kin(:), g2kin_q(:)

  complex(kind=DP), allocatable :: ps(:,:), work (:,:), work_q (:,:)
  complex(kind=DP), allocatable :: becp2(:,:), dvkb(:,:), dvkb1(:,:)

  complex(kind=DP), external :: ZDOTC


  allocate (work(npwx,nkb))    
  allocate (work_q(npwx,nkb))
  allocate (gk(3, npwx))
  allocate (gk_q(3, npwx))
  allocate (g2kin(npwx), g2kin_q(npwx))
  allocate (ps(2,nbnd))
  !     allocate (becp(nkb,nbnd))
  allocate (becp2(nkb,nbnd))
  allocate (dvkb(npwx,nkb))
  allocate (dvkb1(npwx,nkb))


  duk=(0.d0,0.d0)
  work=(0.d0,0.d0)
  work_q=(0.d0,0.d0)
  ps = (0.d0,0.d0)
  
!
! direction of the derivative

  ej=0.d0
  ej(jpol)=1.d0

  call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, phi)

! compute derivative for kpoint

  call gen_us_dj (kpoint, dvkb)
  call gen_us_dy (kpoint, ej, dvkb1)

! compute (k+G)/|k+G|

  do ig = 1, npw
     gk (1, ig) = (xk (1, kpoint) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, kpoint) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, kpoint) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
     if (g2kin (ig) < 1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     else
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     endif
  enddo


  jkb = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
              do ig = 1, npw
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * gk(jpol,ig)
              enddo
           enddo
        endif
     enddo
  enddo

  call ccalbec (nkb, npwx, npw, nbnd, becp2, work, phi)
  
  if (kpoint .ne. kpoint_q) then


! compute derivative for kpoint_q

     call gen_us_dj (kpoint_q, dvkb)
     call gen_us_dy (kpoint_q, ej, dvkb1)
 

! compute (k+q+G)/|k+G|

     do ig = 1, npw
        gk_q (1, ig) = (xk (1, kpoint_q) + g (1, igkq (ig) ) ) * tpiba
        gk_q (2, ig) = (xk (2, kpoint_q) + g (2, igkq (ig) ) ) * tpiba
        gk_q (3, ig) = (xk (3, kpoint_q) + g (3, igkq (ig) ) ) * tpiba
        g2kin_q (ig) = gk_q (1, ig) **2 + gk_q (2, ig) **2 + gk_q (3, ig) **2
        if (g2kin_q (ig) < 1.0d-10) then
           gk_q (1, ig) = 0.d0
           gk_q (2, ig) = 0.d0
           gk_q (3, ig) = 0.d0
        else
           gk_q (1, ig) = gk_q (1, ig) / sqrt (g2kin_q (ig) )
           gk_q (2, ig) = gk_q (2, ig) / sqrt (g2kin_q (ig) )
           gk_q (3, ig) = gk_q (3, ig) / sqrt (g2kin_q (ig) )
        endif
     enddo
     

     jkb = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (nt == ityp (na)) then
              do ikb = 1, nh (nt)
                 jkb = jkb + 1
                 do ig = 1, npw
                    work_q (ig,jkb) = dvkb1(ig, jkb) + dvkb(ig, jkb) &
                         * gk_q(jpol,ig) 
                 enddo
              enddo
           endif
        enddo
     enddo

  else
     work_q=work
  endif
 

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              ps(:,:)=(0.d0,0.d0)
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ibnd = 1, nbnd_occ (kpoint)
                    ps (1, ibnd) = ps(1,ibnd)+ becp2(jkb,ibnd)* &
                         (0.d0,-1.d0)*deeq(ih,jh,na,current_spin)
                    ps (2, ibnd) = ps(2,ibnd) +becp(jkb,ibnd) * &
                         (0.d0,-1.d0)*deeq(ih,jh,na,current_spin)
                 enddo
              enddo
              do ibnd = 1, nbnd_occ (kpoint)
                 call ZAXPY(npw,ps(1,ibnd),vkb_q(1,ikb),1,duk(1,ibnd),1)
                 call ZAXPY(npw,ps(2,ibnd),work_q(1,ikb),1,duk(1,ibnd),1)
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        end if
     end do
  end do
  if (jkb /= nkb) call errore ('take_nloc_k_kq', 'unexpected error', 1)

  call grad(kpoint, phi, jpol, duk)

  return

end subroutine take_nloc_k_kq
