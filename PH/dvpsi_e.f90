!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvpsi_e (kpoint, ipol)  
  !----------------------------------------------------------------------
  !
  ! Calculates x * psi_k  for each k-point and polarization
  ! Requires on input: vkb, evc, igk
  ! On output: dvpsi is in crystal axis (projected on at(*,ipol) )
  !
#include "machine.h"
  !

  use pwcom 
  use allocate 
  use parameters, only : DP 
  use becmod
  use phcom  
  implicit none 
  !
  integer :: ipol, ig, na, ibnd, jbnd, ikb, jkb, nt, lter, kpoint, &
                   ih, jh, ijkb0 
  ! input: the polarization
  ! counter on G vectors
  ! counter on atoms
  ! counters on bands
  ! counter on beta functions
  ! counter on beta functions
  ! counter on pseudopotentials
  ! iterations of linter
  ! input: the kpoint

  real(kind=DP), pointer  :: gk (:,:), h_diag (:,:),  eprec (:) 
  ! the derivative of |k+G|
  real(kind=DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter

  logical :: conv_root  
  ! true if convergence has been achieved

  complex(kind=DP), pointer :: ps (:,:), dvkb (:,:), dvkb1 (:,:), work (:,:), &
                               becp2(:,:), spsi(:,:)
  complex(kind=DP) :: ZDOTC
  ! the scalar products
  external ch_psi_all, cg_psi  
  !
  !
  call start_clock ('dvpsi_e')  
  call mallocate(work , npwx, nkb)  
  call mallocate(gk , 3, npwx)  
  call mallocate(h_diag, npwx , nbnd)  
  call mallocate(ps , 2 , nbnd)  
  call mallocate(spsi , npwx, nbnd)  
  call mallocate(eprec , nbnd)  
  if (nkb.gt.0) call mallocate(becp2 , nkb, nbnd)  
  if (nkb.gt.0) call mallocate(dvkb , npwx , nkb)  
  if (nkb.gt.0) call mallocate(dvkb1, npwx , nkb)  
  call setv (2 * npwx * nkb, 0.d0, dvkb, 1)  
  call setv (2 * npwx * nkb, 0.d0, dvkb1, 1)  
  call setv (2 * npwx * nbnd, 0.d0, dvpsi, 1)  
  do ig = 1, npw  
     gk (1, ig) = (xk (1, kpoint) + g (1, igk (ig) ) ) * tpiba  
     gk (2, ig) = (xk (2, kpoint) + g (2, igk (ig) ) ) * tpiba  
     gk (3, ig) = (xk (3, kpoint) + g (3, igk (ig) ) ) * tpiba  
!     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2  
  enddo
  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  do ibnd = 1, nbnd_occ (kpoint)  
     do ig = 1, npw  
        dpsi (ig, ibnd) = (at(1, ipol) * gk(1, ig) + &
                           at(2, ipol) * gk(2, ig) + &
                           at(3, ipol) * gk(3, ig) ) &
                           *(0.d0,-2.d0)*evc (ig, ibnd)
     enddo
  enddo
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  call gen_us_dj (kpoint, dvkb)  
  call gen_us_dy (kpoint, at (1, ipol), dvkb1)  
  do ig = 1, npw  
     if (g2kin (ig) .lt.1.0d-10) then  
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
        if (nt.eq.ityp (na)) then
           do ikb = 1, nh (nt)  
              jkb = jkb + 1  
              do ig = 1, npw  
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                            (at (1, ipol) * gk (1, ig) + &
                             at (2, ipol) * gk (2, ig) + &
                             at (3, ipol) * gk (3, ig) )
              enddo
           enddo
        endif
     enddo
  enddo
  call ccalbec (nkb, npwx, npw, nbnd, becp2, work, evc)

  ijkb0 = 0
  do nt = 1, ntyp  
     do na = 1, nat  
        if (nt.eq.ityp (na)) then
           do ih = 1, nh (nt)  
              ikb = ijkb0 + ih  
              ps(:,:)=(0.d0,0.d0)
              do jh = 1, nh (nt)  
                 jkb = ijkb0 + jh  
                 do ibnd = 1, nbnd_occ (kpoint)  
                    ps (1, ibnd) = ps(1,ibnd)+ becp2(jkb,ibnd)* &
                           (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                           -et(ibnd,kpoint)*qq(ih,jh,nt))
                    ps (2, ibnd) = ps(2,ibnd) +becp1(jkb,ibnd,kpoint) * &
                           (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                           -et(ibnd,kpoint)*qq(ih,jh,nt))
                 enddo
              enddo
              do ibnd = 1, nbnd_occ (kpoint)  
                 call ZAXPY(npw,ps(1,ibnd),vkb(1,ikb),1,dpsi(1,ibnd),1)
                 call ZAXPY(npw,ps(2,ibnd),work(1,ikb),1,dpsi(1,ibnd),1)  
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        end if
     end do
  end do
  if (jkb.ne.nkb) call error ('dvpsi_e', 'unexpected error', 1)  
  !
  !    orthogonalize dpsi to the valence subspace
  !
  do ibnd = 1, nbnd_occ (kpoint)  
     call setv (2 * npwx, 0.d0, work, 1)  
     do jbnd = 1, nbnd_occ (kpoint)  
        ps (1, jbnd) = - ZDOTC(npw,evc(1,jbnd),1,dpsi(1, ibnd),1)
     enddo
#ifdef PARA
     call reduce (4 * nbnd, ps)  
#endif
     do jbnd = 1, nbnd_occ (kpoint)  
        call ZAXPY (npw, ps (1, jbnd), evc (1, jbnd), 1, work, 1)  
     enddo
     call ccalbec (nkb, npwx, npw, 1, becp, vkb, work)
     call s_psi (npwx, npw, 1, work, spsi)
     call DAXPY (2 * npw, 1.0d0, spsi, 1, dpsi (1, ibnd), 1)  
  enddo
  !
  !   dpsi contains now P^+_c [H-eS,x] psi_v  for the three crystal 
  !   polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  thresh = 1.d-5  
  do ibnd = 1, nbnd_occ (kpoint)  
     conv_root = .true.  
     do ig = 1, npwq  
        work (ig,1) = g2kin (ig) * evc (ig, ibnd)  
     enddo
     eprec (ibnd) = 1.35d0 * ZDOTC (npwq, evc (1, ibnd), 1, work, 1)  
  enddo
#ifdef PARA
  call reduce (nbnd_occ (kpoint), eprec)  
#endif
  do ibnd = 1, nbnd_occ (kpoint)  
     do ig = 1, npwq  
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )  
     enddo
  enddo
  !
  call cgsolve_all (ch_psi_all, cg_psi, et (1, kpoint), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, kpoint, lter, conv_root, anorm, &
       nbnd_occ (kpoint) )

  if (.not.conv_root) write (6, '(5x,"kpoint",i4," ibnd",i4, &
       &                         " linter: root not converged ",e10.3)') &
       & kpoint, ibnd, anorm
#ifdef FLUSH
  call flush (6)  
#endif
!
! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>, 
! therefore we apply S again, and then subtract the additional term
! furthermore we add the term due to dipole of the augmentation charges.
!  
  if (okvan) then
     call ccalbec(nkb,npwx,npw,nbnd,becp,vkb,dvpsi)
     call s_psi(npwx,npw,nbnd,dvpsi,spsi)
     call DCOPY(2*npwx*nbnd,spsi,1,dvpsi,1)
     call adddvepsi_us(becp2,ipol,kpoint)
  endif

  if (nkb.gt.0) call mfree (dvkb1)  
  if (nkb.gt.0) call mfree (dvkb)  
  if (nkb.gt.0) call mfree (becp2)  
  call mfree (eprec)  
  call mfree (spsi)  
  call mfree (ps)  
  call mfree (h_diag)  
  call mfree (gk)  
  call mfree (work)  

  call stop_clock ('dvpsi_e')  
  return  
end subroutine dvpsi_e
