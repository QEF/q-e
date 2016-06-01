!
! Copyright (C) 2003-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine psidspsi (ik, uact, pdsp)
!----------========----------------------------------------------------
  !
  ! This routine calculates <psi_v'|ds/du|psi_v>
  ! at q=0. The displacements are described by a vector uact.
  ! The result is stored in pdsp. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : g
  USE klist,     ONLY : xk, ngk, igk_k
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,  ONLY : lsda, current_spin, isk
  USE spin_orb,  ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions_module,    ONLY : evc
  USE wvfct,     ONLY : nbnd, npwx
  USE uspp,      ONLY: nkb, vkb, qq, qq_so
  USE uspp_param,ONLY : nh
  USE phus,      ONLY : alphap

  USE lrus,       ONLY : becp1
  USE control_lr, ONLY : lgamma

  implicit none
  !
  !   The dummy variables
  !

  integer, intent(in) :: ik
  ! input: the k point
  complex(DP) :: uact (3 * nat), pdsp(nbnd,nbnd)
  ! input: the pattern of displacements
  ! output: <psi|ds/du|psi>
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, jbnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, npw
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! the point k+q
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(DP), parameter :: eps = 1.d-12

  complex(DP), ALLOCATABLE :: ps1 (:,:), ps2 (:,:,:), aux (:), aux1(:,:), &
                              dspsi(:,:)
  complex(DP), ALLOCATABLE :: ps1_nc(:,:,:), ps2_nc(:,:,:,:)
  ! the scalar product
  ! the scalar product
  ! a mesh space for psi
  ! the matrix dspsi

  logical :: ok
  ! used to save time

  if (noncolin) then
     allocate (ps1_nc ( nkb, npol, nbnd ))
     allocate (ps2_nc ( nkb, npol, 3, nbnd))
  else
     allocate (ps1 ( nkb, nbnd ))
     allocate (ps2 ( nkb, 3, nbnd))
  endif
  allocate (dspsi (npwx*npol, nbnd))
  allocate (aux ( npwx*npol ))

  if (lgamma) then
     ikk = ik
     ikq = ik
     npw = ngk(ik)
  else
     call infomsg ('psidspsi', 'called for lgamma .eq. false')
  endif
  if (lsda) current_spin = isk (ikk)

  if (noncolin) then
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
  else
     ps1(:,:)   = (0.d0, 0.d0)
     ps2(:,:,:) = (0.d0, 0.d0)
  endif
  pdsp(:,:)   = (0.d0, 0.d0)
  dspsi = (0.d0,0.d0)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           mu = 3 * (na - 1)
           if ( abs (uact (mu + 1) ) + &
                abs (uact (mu + 2) ) + &
                abs (uact (mu + 3) ) > eps) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ipol = 1, 3
                       do ibnd = 1, nbnd
                          if (noncolin) then
                             if (lspinorb) then
                                ps1_nc(ikb,1,ibnd)=ps1_nc(ikb,1,ibnd) +     &
                                  (qq_so(ih,jh,1,nt)*                    &
                                  alphap(ipol,ik)%nc(jkb,1,ibnd)+         &
                                  qq_so(ih,jh,2,nt)*                     &
                                  alphap(ipol,ik)%nc(jkb,2,ibnd) )*       &
                                  uact (mu + ipol)
                                ps1_nc(ikb,2,ibnd)=ps1_nc(ikb,2,ibnd) +     &
                                  (qq_so(ih,jh,3,nt)*                    &
                                  alphap(ipol,ik)%nc(jkb,1,ibnd)+         &
                                  qq_so(ih,jh,4,nt)*                     &
                                  alphap(ipol,ik)%nc(jkb,2,ibnd) )*       &
                                  uact (mu + ipol)
                                ps2_nc(ikb,1,ipol,ibnd)=                 &
                                      ps2_nc(ikb,1,ipol,ibnd) +          &
                                  (qq_so (ih, jh, 1, nt) *               &
                                  becp1(ik)%nc (jkb, 1, ibnd) +          &
                                   qq_so (ih, jh, 2, nt) *               &
                                  becp1(ik)%nc (jkb, 2, ibnd) )*         &
                                  (0.d0, -1.d0)* uact (mu + ipol) * tpiba
                                ps2_nc(ikb,2,ipol,ibnd)=                 &
                                      ps2_nc(ikb,2,ipol,ibnd) +          &
                                  (qq_so (ih, jh, 3, nt) *               &
                                  becp1(ik)%nc (jkb, 1, ibnd) +          &
                                   qq_so (ih, jh, 4, nt) *               &
                                  becp1(ik)%nc (jkb, 2, ibnd) )*         &
                                  (0.d0, -1.d0)* uact (mu + ipol) * tpiba
                             else
                                do is=1,npol
                                   ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd) +   &
                                     qq(ih,jh,nt)*                          &
                                     alphap(ipol,ik)%nc(jkb,is,ibnd)*        &
                                     uact (mu + ipol)
                                   ps2_nc(ikb,is,ipol,ibnd)=                 &
                                         ps2_nc(ikb,is,ipol,ibnd) +          &
                                     qq (ih, jh, nt) *(0.d0, -1.d0)*         &
                                     becp1(ik)%nc (jkb,is,ibnd) *            &
                                     uact (mu + ipol) * tpiba
                                enddo
                             endif
                          else
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +    &
                               qq (ih, jh, nt) *                 &
                               alphap(ipol,ik)%k(jkb,ibnd) *     &
                               uact (mu + ipol)
                             ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) + &
                               qq (ih, jh, nt) *                          &
                               (0.d0, -1.d0) *                            &
                               becp1(ik)%k (jkb, ibnd) *                  &
                               uact (mu + ipol) * tpiba
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0= ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dspsi, npwx)
     else
        call zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dspsi, npwx)
!        dspsi = matmul(vkb,ps1)+ dspsi
     endif
  endif
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           if (noncolin) then
              ok = ok.or. (ABS (ps2_nc (ikb, 1, ipol, ibnd) ) .gt.eps) &
                     .or. (ABS (ps2_nc (ikb, 2, ipol, ibnd) ) .gt.eps)
           else
              ok = ok.or. (ABS (ps2 (ikb, ipol, ibnd) ) .gt.eps)
           endif
        enddo
        if (ok) then
           do ig = 1, npw
              igg = igk_k (ig,ikk)
              aux (ig) =  vkb(ig, ikb) *    &
                   (xk(ipol, ik) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              if (noncolin) then
                 dspsi(1:npw,ibnd) = ps2_nc(ikb,1,ipol,ibnd) * aux(1:npw) &
                      + dspsi(1:npw,ibnd)
                 dspsi(1+npwx:npw+npwx,ibnd) = ps2_nc(ikb,2,ipol,ibnd)* &
                                aux(1:npw) + dspsi(1+npwx:npw+npwx,ibnd)
              else
                 dspsi(1:npw,ibnd) = ps2(ikb,ipol,ibnd) * aux(1:npw) &
                      + dspsi(1:npw,ibnd)
              endif
           enddo
        endif
     enddo

  enddo
  do ibnd = 1, nbnd
     do jbnd=1, nbnd
        pdsp(ibnd,jbnd) =  &
                dot_product(evc(1:npwx*npol,ibnd),dspsi(1:npwx*npol,jbnd))
     enddo
  enddo

  if (allocated(aux)) deallocate (aux)
  if (noncolin) then
     if (allocated(ps2_nc)) deallocate (ps2_nc)
     if (allocated(ps1_nc)) deallocate (ps1_nc)
  else
     if (allocated(ps2)) deallocate (ps2)
     if (allocated(ps1)) deallocate (ps1)
  endif
  if (allocated(dspsi)) deallocate (dspsi)

  return
end subroutine psidspsi
