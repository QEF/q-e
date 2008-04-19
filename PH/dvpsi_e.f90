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
subroutine dvpsi_e (ik, ipol)
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis 
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is READ from file if this_pcxpsi_is_on_file(ik,ipol)=.true. 
  ! otherwise dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set)
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE io_global,      ONLY : stdout
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE noncollin_module,  ONLY : noncolin, npol
  USE kinds, only : DP
  USE becmod, ONLY: becp, becp_nc, calbec
  USE uspp_param, ONLY: nh
  USE ramanm, ONLY: eth_rps
  use phcom
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  implicit none
  !
  integer, intent(IN) :: ipol, ik
  !
  ! Local variables
  !
  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, nrec
  ! counters

  real(DP), allocatable  :: gk (:,:), h_diag (:,:),  eprec (:)
  ! the derivative of |k+G|
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter

  logical :: conv_root
  ! true if convergence has been achieved

  complex(DP), allocatable :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),  &
       work (:,:), becp2(:,:), becp2_nc(:,:,:), spsi(:,:), ps(:,:), &
       psc(:,:,:,:), aux(:)
  complex(DP), external :: ZDOTC
  ! the scalar products
  external ch_psi_all, cg_psi
  !
  !
  call start_clock ('dvpsi_e')
  dpsi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  if (this_pcxpsi_is_on_file(ik,ipol)) then
     nrec = (ipol - 1)*nksq + ik
     call davcio(dvpsi, lrebar, iuebar, nrec, -1)
     call stop_clock ('dvpsi_e')
     return
  end if
  !
  allocate (work ( npwx, MAX(nkb,1)))
  allocate (aux ( npwx*npol ))
  allocate (gk ( 3, npwx))    
  allocate (h_diag( npwx*npol, nbnd))    
  allocate (eprec( nbnd))    
  if (nkb > 0) then
     IF (noncolin) THEN
        allocate (becp2_nc (nkb, npol, nbnd))
     ELSE
        allocate (becp2 (nkb, nbnd))
     END IF
     allocate (dvkb (npwx, nkb), dvkb1(npwx, nkb))
     dvkb (:,:) = (0.d0, 0.d0)
     dvkb1(:,:) = (0.d0, 0.d0)
  end if
  do ig = 1, npw
     gk (1, ig) = (xk (1, ik) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, ik) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, ik) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  enddo
  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        dpsi (ig, ibnd) = (at(1, ipol) * gk(1, ig) + &
             at(2, ipol) * gk(2, ig) + &
             at(3, ipol) * gk(3, ig) ) &
             *(0.d0,-2.d0)*evc (ig, ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, npw
           dpsi (ig+npwx, ibnd) = (at(1, ipol) * gk(1, ig) + &
                at(2, ipol) * gk(2, ig) + &
                at(3, ipol) * gk(3, ig) ) &
                 *(0.d0,-2.d0)*evc (ig+npwx, ibnd)
        end do
     END IF
  enddo
!
! Uncomment this goto and the continue below to calculate 
! the matrix elements of p without the commutator with the
! nonlocal potential.
!
!  goto 111
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  call gen_us_dj (ik, dvkb)
  call gen_us_dy (ik, at (1, ipol), dvkb1)
  do ig = 1, npw
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
  work=(0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
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
  deallocate (gk)

  IF ( noncolin ) THEN 
     call calbec (npw, work, evc, becp2_nc)
  ELSE
     call calbec (npw, work, evc, becp2)
  END IF

  ijkb0 = 0
  IF (noncolin) THEN
     allocate (psc ( nkb, npol, nbnd, 2))
     psc=(0.d0,0.d0)
  ELSE
     allocate (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
  END IF
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ibnd = 1, nbnd_occ (ik)
                    IF (noncolin) THEN
                       IF (lspinorb) THEN 
                          psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+(0.d0,-1.d0)* &
                             (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+       &
                              becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)- &
                                       et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                          psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+(0.d0,-1.d0)*  &
                             (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+       &
                              becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)- &
                                       et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                          psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+(0.d0,-1.d0)* &
                             (becp1_nc(jkb,1,ibnd,ik)*(deeq_nc(ih,jh,na,1)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+      &
                             becp1_nc(jkb,2,ibnd,ik)*(deeq_nc(ih,jh,na,2)-  &
                                       et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                          psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+(0.d0,-1.d0)*  &
                             (becp1_nc(jkb,1,ibnd,ik)*(deeq_nc(ih,jh,na,3)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+      &
                             becp1_nc(jkb,2,ibnd,ik)*(deeq_nc(ih,jh,na,4)-  &
                                       et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                       ELSE
                          psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+ (0.d0,-1.d0)* &
                              ( becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                             -et(ibnd,ik)*qq(ih,jh,nt)) + &
                                becp2_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                          psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+ (0.d0,-1.d0)* &
                              ( becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+  &
                                becp2_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                          psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+ (0.d0,-1.d0)* &
                              ( becp1_nc(jkb,1,ibnd,ik)*(deeq_nc(ih,jh,na,1) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                becp1_nc(jkb,2,ibnd,ik)*deeq_nc(ih,jh,na,2) )
                          psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+ (0.d0,-1.d0)* &
                              ( becp1_nc(jkb,2,ibnd,ik)*(deeq_nc(ih,jh,na,4) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                becp1_nc(jkb,1,ibnd,ik)*deeq_nc(ih,jh,na,3) )
                       END IF
                    ELSE
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                           (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                           -et(ibnd,ik)*qq(ih,jh,nt))
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                           (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                           -et(ibnd,ik)*qq(ih,jh,nt))
                    END IF
                 enddo
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        end if
     end do
  end do
  if (ikb /= nkb .OR. jkb /= nkb) call errore ('dvpsi_e', 'unexpected error',1)
  IF (noncolin) THEN
     CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik)*npol, nkb, &
          (1.d0,0.d0), vkb(1,1), npwx, psc(1,1,1,1), nkb, (1.d0,0.d0), &
          dpsi, npwx )
     CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik)*npol, nkb, &
          (1.d0,0.d0),work(1,1), npwx, psc(1,1,1,2), nkb, (1.d0,0.d0), &
          dpsi, npwx )
  ELSE
     CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
          (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
          dpsi(1,1), npwx )
     CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
          (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
          dpsi(1,1), npwx )
  ENDIF

  IF (noncolin) THEN
     deallocate (psc)
  ELSE
     deallocate (ps2)
  END IF

!  111 continue
  !
  !    orthogonalize dpsi to the valence subspace: ps = <evc|dpsi>
  !
  allocate (ps ( nbnd, nbnd ))
  IF (noncolin) THEN
     CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwx*npol, &
          (1.d0,0.d0), evc, npwx*npol, dpsi, npwx*npol, (0.d0,0.d0), &
          ps, nbnd )
  ELSE
     CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
          (1.d0,0.d0), evc(1,1), npwx, dpsi(1,1), npwx, (0.d0,0.d0), &
          ps(1,1), nbnd )
  END IF
#ifdef __PARA
  call mp_sum ( ps(:, 1:nbnd_occ(ik) ), intra_pool_comm )
#endif
  ! dvpsi is used as work space to store S|evc>
  !
  IF (noncolin) THEN
     CALL calbec ( npw, vkb, evc, becp_nc, nbnd_occ(ik) )
     CALL s_psi_nc (npwx, npw, nbnd_occ(ik), evc, dvpsi)
  ELSE
     CALL calbec ( npw, vkb, evc, becp, nbnd_occ(ik) )
     CALL s_psi (npwx, npw, nbnd_occ(ik), evc, dvpsi)
  END IF
  !
  ! |dpsi> = |dpsi> - S|evc><evc|dpsi>)
  !
  IF (noncolin) THEN
     CALL ZGEMM( 'N', 'N', npwx*npol, nbnd_occ(ik), nbnd_occ(ik), &
          (-1.d0,0.d0), dvpsi(1,1), npwx*npol, ps(1,1), nbnd, (1.d0,0.d0), &
          dpsi(1,1), npwx*npol )
  ELSE
     CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
          (-1.d0,0.d0), dvpsi(1,1), npwx, ps(1,1), nbnd, (1.d0,0.d0), &
          dpsi(1,1), npwx )
  END IF
  deallocate (ps)
  !
  !   dpsi contains P^+_c [H-eS,x] psi_v for the three crystal polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  thresh = eth_rps
  do ibnd = 1, nbnd_occ (ik)
     conv_root = .true.
     aux=(0.d0,0.d0)
     do ig = 1, npwq
        aux(ig) = g2kin (ig) * evc (ig, ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, npwq
           aux(ig+npwx) = g2kin (ig) * evc (ig+npwx, ibnd)
        enddo
     END IF
     eprec (ibnd) = 1.35d0 * ZDOTC (npwx*npol, evc (1, ibnd), 1, aux, 1)
  enddo
#ifdef __PARA
  call mp_sum ( eprec( 1:nbnd_occ(ik) ), intra_pool_comm )
#endif
  h_diag=0.d0
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npwq
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     enddo
     IF (noncolin) THEN
        do ig = 1, npwq
           h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
        enddo
     END IF
  enddo
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  call cgsolve_all (ch_psi_all, cg_psi, et (1, ik), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ (ik), npol)

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm
  !
  CALL flush_unit( stdout )
  !
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born 
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  if (okvan) then
     !
     ! for effective charges
     !
     nrec = (ipol - 1) * nksq + ik
     call davcio (dvpsi, lrcom, iucom, nrec, 1)
     !
     allocate (spsi ( npwx*npol, nbnd))    
     IF (noncolin) THEN
        CALL calbec (npw, vkb, dvpsi, becp_nc )
        CALL s_psi_nc(npwx,npw,nbnd,dvpsi,spsi)
     ELSE
        CALL calbec (npw, vkb, dvpsi, becp )
        CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
     END IF
     call DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
     deallocate (spsi)
     IF (noncolin) THEN
        call adddvepsi_us(becp2_nc,ipol,ik)
     ELSE
        call adddvepsi_us(becp2,ipol,ik)
     END IF
  endif


  IF (nkb > 0) THEN
     deallocate (dvkb1, dvkb)
     IF (noncolin) THEN
        deallocate(becp2_nc)
     ELSE
        deallocate(becp2)
     ENDIF
  END IF

  deallocate (eprec)
  deallocate (h_diag)
  deallocate (work)
  deallocate (aux)

  nrec = (ipol - 1)*nksq + ik
  call davcio(dvpsi, lrebar, iuebar, nrec, 1)
  this_pcxpsi_is_on_file(ik,ipol) = .true.
  call stop_clock ('dvpsi_e')
  return
end subroutine dvpsi_e
