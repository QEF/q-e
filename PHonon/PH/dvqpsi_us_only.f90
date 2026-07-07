!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvqpsi_us_only (ik, uact, becp1, alphap)
   !----------------------------------------------------------------------
   !! This routine calculates \(\text{dV_bare}/\text{dtau}\cdot\text{psi}\)
   !! for one perturbation with a given q.  
   !! The displacements are described by a vector uact.  
   !! The result is stored in \(\text{dvpsi}\). The routine is called for
   !! each k-point and for each pattern u. It computes simultaneously all
   !! the bands.  
   !! This routine implements Eq. (B29) of PRB 64, 235118 (2001).
   !! Only the contribution of the nonlocal potential is calculated here;
   !! both norm-conserving term and ultrasoft correction are calculated here.
   !
   !
   USE kinds, only : DP
   USE cell_base, ONLY : tpiba
   USE gvect,     ONLY : g
   USE klist,     ONLY : xk, ngk, igk_k
   USE ions_base, ONLY : nat, ityp, ntyp => nsp
   USE lsda_mod,  ONLY : lsda, current_spin, isk, nspin
   USE wvfct,     ONLY : nbnd, npwx, et
   USE noncollin_module, ONLY : noncolin, npol, lspinorb
   USE uspp, ONLY: okvan, nkb, vkb, ofsbeta
   USE uspp_param, ONLY: nh, nhm
   USE phus,      ONLY : int1, int1_nc, int2, int2_so

   USE qpoint,     ONLY : nksq, ikks, ikqs
   USE becmod,     ONLY : bec_type
   USE eqv,        ONLY : dvpsi
   USE control_lr, ONLY : lgamma

   IMPLICIT NONE
   !
   integer :: ik
   !! input: the k point
   complex(DP) :: uact(3*nat)
   !! input: the pattern of displacements
   TYPE(bec_type) :: becp1(nksq), alphap(3,nksq)
   !
   !   ... local variables
   !

   integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
        ikb, jkb, ih, jh, ipol, is, js, ijs, npwq
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

   complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:), deff_nc(:,:,:,:)
   real(DP), allocatable :: deff(:,:,:)
   complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
   ! work space
   complex(DP), allocatable :: temp_ps1 (:), temp_ps2 (:,:)
   complex(DP), allocatable :: temp_ps1_nc (:,:), temp_ps2_nc (:,:,:)
   ! temporary buffers
   complex(DP), allocatable :: temp_alphap_nc (:,:), temp_becp1_nc (:,:)
   complex(DP), allocatable :: temp_alphap_k (:), temp_becp1_k (:)
   ! temporary buffers for alphap and becp1
   real(DP), allocatable :: temp_reduction (:)
   ! temporary buffer for reduction operations
   complex(DP), allocatable :: temp_ps2_values (:), temp_ps2_nc_values (:,:)
   complex(DP), allocatable :: temp_dvpsi (:)
   ! temporary buffers for ps2 values and dvpsi

   logical :: ok

   call start_clock ('dvqpsi_us_on')
   IF (noncolin) THEN
      allocate (ps1_nc(nkb , npol, nbnd))
      allocate (ps2_nc(nkb , npol, nbnd , 3))
      allocate (deff_nc(nhm, nhm, nat, nspin))
      allocate (temp_ps1_nc(nkb, npol))
      allocate (temp_ps2_nc(nkb, npol, 3))
      allocate (temp_alphap_nc(nkb, npol))
      allocate (temp_becp1_nc(nkb, npol))
   ELSE
      allocate (ps1 ( nkb , nbnd))
      allocate (ps2 ( nkb , nbnd , 3))
      allocate (deff(nhm, nhm, nat))
      allocate (temp_ps1(nkb))
      allocate (temp_ps2(nkb, 3))
      allocate (temp_alphap_k(nkb))
      allocate (temp_becp1_k(nkb))
   END IF
   allocate (aux ( npwx))
   allocate (temp_reduction(nbnd))
   IF (noncolin) THEN
      allocate (temp_ps2_nc_values(nbnd, npol))
   ELSE
      allocate (temp_ps2_values(nbnd))
   END IF
   allocate (temp_dvpsi(npwx))
   ikk = ikks(ik)
   ikq = ikqs(ik)
   IF (lsda) current_spin = isk (ikk)
   !
   !   we first compute the coefficients of the vectors
   !
   IF (noncolin) THEN
      ps1_nc(:,:,:)   = (0.d0, 0.d0)
      ps2_nc(:,:,:,:) = (0.d0, 0.d0)
   ELSE
      ps1(:,:)   = (0.d0, 0.d0)
      ps2(:,:,:) = (0.d0, 0.d0)
   END IF
   DO ibnd = 1, nbnd
      ! Initialize temporary buffers
      IF (noncolin) THEN
         temp_ps1_nc(:,:) = (0.d0, 0.d0)
         temp_ps2_nc(:,:,:) = (0.d0, 0.d0)
         temp_becp1_nc(:,:) = becp1(ik)%nc(:,:,ibnd)
      ELSE
         temp_ps1(:) = (0.d0, 0.d0)
         temp_ps2(:,:) = (0.d0, 0.d0)
         temp_becp1_k(:) = becp1(ik)%k(:,ibnd)
      END IF
      
      IF (noncolin) THEN
         CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
      ELSE
         CALL compute_deff(deff,et(ibnd,ikk))
      ENDIF
      DO na = 1, nat
         ijkb0 = ofsbeta(na) 
         nt = ityp(na) 
         mu = 3 * (na - 1)
         ! First loop: deff calculation
         IF ( abs (uact (mu + 1) ) + &
              abs (uact (mu + 2) ) + &
              abs (uact (mu + 3) ) > eps) THEN
            DO ih = 1, nh (nt)
               ikb = ijkb0 + ih
               DO jh = 1, nh (nt)
                  jkb = ijkb0 + jh
                  DO ipol = 1, 3
                     ! Copy alphap for current ipol and ibnd
                     IF (noncolin) THEN
                        temp_alphap_nc(:,:) = alphap(ipol, ik)%nc(:,:,ibnd)
                     ELSE
                        temp_alphap_k(:) = alphap(ipol, ik)%k(:,ibnd)
                     END IF
                     
                     IF (noncolin) THEN
                        ijs=0
                        DO is=1,npol
                           DO js=1,npol
                              ijs=ijs+1
                              temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is) +  &
                                 deff_nc(ih,jh,na,ijs) * &
                                 temp_alphap_nc(jkb,js)* &
                                  uact(mu + ipol)
                              temp_ps2_nc(ikb,is,ipol)=               &
                                     temp_ps2_nc(ikb,is,ipol)+        &
                                     deff_nc(ih,jh,na,ijs) *          &
                                     temp_becp1_nc(jkb,js) *      &
                                     (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                           END DO
                        END DO
                     ELSE
                        temp_ps1 (ikb) = temp_ps1 (ikb) +      &
                                   deff(ih, jh, na) *            &
                           temp_alphap_k(jkb) * uact (mu + ipol)
                        temp_ps2 (ikb, ipol) = temp_ps2 (ikb, ipol) +&
                             deff(ih,jh,na)*temp_becp1_k(jkb) *  &
                             (0.0_DP,-1.0_DP) * uact (mu + ipol) * tpiba
                     ENDIF
                     IF (okvan) THEN
                        IF (noncolin) THEN
                           ijs=0
                           DO is=1,npol
                              DO js=1,npol
                                 ijs=ijs+1
                                 temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is)+ &
                                    int1_nc(ih,jh,ipol,na,ijs) *     &
                                    temp_becp1_nc(jkb,js)*uact(mu+ipol)
                              END DO
                           END DO
                        ELSE
                           temp_ps1 (ikb) = temp_ps1 (ikb) + &
                             (int1 (ih, jh, ipol,na, current_spin) * &
                             temp_becp1_k(jkb) ) * uact (mu +ipol)
                        END IF
                     END IF
                  ENDDO ! ipol
               ENDDO ! jh
            ENDDO ! ih
         END IF  ! uact>0
      enddo  ! na
      ! Second loop: okvan nb loop
      IF (okvan) THEN
         DO nb = 1, nat
            nu = 3 * (nb - 1)
            IF ( abs (uact (nu + 1) ) + &
                 abs (uact (nu + 2) ) + &
                 abs (uact (nu + 3) ) > eps) THEN
               DO na = 1, nat
                  ijkb0 = ofsbeta(na) 
                  nt = ityp(na) 
                  DO ih = 1, nh (nt)
                     ikb = ijkb0 + ih
                     DO jh = 1, nh (nt)
                        jkb = ijkb0 + jh
                        DO ipol = 1, 3
                           IF (noncolin) THEN
                              IF (lspinorb) THEN
                                 ijs=0
                                 DO is=1,npol
                                    DO js=1,npol
                                       ijs=ijs+1
                                       temp_ps1_nc(ikb,is)= &
                                                 temp_ps1_nc(ikb,is)+ &
                                       int2_so(ih,jh,ipol,nb,na,ijs)* &
                                        temp_becp1_nc(jkb,js)*uact(nu+ipol)
                                    END DO
                                 END DO
                              ELSE
                                 DO is=1,npol
                                    temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is)+ &
                                       int2(ih,jh,ipol,nb,na) * &
                                       temp_becp1_nc(jkb,is)*uact(nu+ipol)
                                 END DO
                              END IF
                           ELSE
                              temp_ps1 (ikb) = temp_ps1 (ikb) + &
                                  (int2 (ih, jh, ipol, nb, na) * &
                                   temp_becp1_k(jkb) ) * uact (nu + ipol)
                           END IF
                        ENDDO ! ipol
                     ENDDO ! jh
                  ENDDO ! ih
               ENDDO  ! na
            END IF  ! uact>0
         ENDDO  ! nb
      ENDIF  ! okvan
      
      ! Assign temporary buffers to original arrays
      IF (noncolin) THEN
         ps1_nc(:,:,ibnd) = temp_ps1_nc(:,:)
         ps2_nc(:,:,ibnd,:) = temp_ps2_nc(:,:,:)
      ELSE
         ps1(:,ibnd) = temp_ps1(:)
         ps2(:,ibnd,:) = temp_ps2(:,:)
      END IF
   ENDDO ! nbnd
   !
   !      This term is proportional to beta(k+q+G)
   ! 
   npwq = ngk(ikq)
   IF (nkb.gt.0) THEN
      IF (noncolin) THEN
         call zgemm ('N', 'N', npwq, nbnd*npol, nkb, &
          (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
      ELSE
         call zgemm ('N', 'N', npwq, nbnd, nkb, &
          (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
      END IF
   END IF
   !
   !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
   !
   DO ipol = 1, 3
      DO ikb = 1, nkb
         ok = .false.
         IF (noncolin) THEN
            ! Copy data to contiguous buffer for reduction
            temp_reduction(1:nbnd) = abs(ps2_nc(ikb, 1, 1:nbnd, ipol))
            ok = any(temp_reduction(1:nbnd) .gt. eps)
            if (.not. ok) then
               temp_reduction(1:nbnd) = abs(ps2_nc(ikb, 2, 1:nbnd, ipol))
               ok = any(temp_reduction(1:nbnd) .gt. eps)
            ENDIF
         ELSE
            ! Copy data to contiguous buffer for reduction
            temp_reduction(1:nbnd) = abs(ps2(ikb, 1:nbnd, ipol))
            ok = any(temp_reduction(1:nbnd) .gt. eps)
         ENDIF
         IF (ok) THEN
            DO ig = 1, npwq
               igg = igk_k (ig,ikq)
               aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
            ENDDO
            
            ! Copy ps2 values to contiguous buffer
            IF (noncolin) THEN
               temp_ps2_nc_values(1:nbnd, 1) = ps2_nc(ikb, 1, 1:nbnd, ipol)
               temp_ps2_nc_values(1:nbnd, 2) = ps2_nc(ikb, 2, 1:nbnd, ipol)
            ELSE
               temp_ps2_values(1:nbnd) = ps2(ikb, 1:nbnd, ipol)
            END IF
            
            DO ibnd = 1, nbnd
               ! Initialize temp_dvpsi to zero for each band
               temp_dvpsi(1:npwq) = (0.0_DP, 0.0_DP)
               
               IF (noncolin) THEN
                  call zaxpy(npwq, temp_ps2_nc_values(ibnd,1), aux, 1, temp_dvpsi, 1)
                  ! Copy temp_dvpsi to dvpsi for first spin component
                  dvpsi(1:npwq, ibnd) = dvpsi(1:npwq, ibnd) + temp_dvpsi(1:npwq)
                  
                  ! Reset temp_dvpsi for second spin component
                  temp_dvpsi(1:npwq) = (0.0_DP, 0.0_DP)
                  call zaxpy(npwq, temp_ps2_nc_values(ibnd,2), aux, 1, temp_dvpsi, 1)
                  ! Copy temp_dvpsi to dvpsi for second spin component
                  dvpsi(1+npwx:npwx+npwq, ibnd) = dvpsi(1+npwx:npwx+npwq, ibnd) + temp_dvpsi(1:npwq)
               ELSE
                  call zaxpy(npwq, temp_ps2_values(ibnd), aux, 1, temp_dvpsi, 1)
                  ! Copy temp_dvpsi to dvpsi
                  dvpsi(1:npwq, ibnd) = dvpsi(1:npwq, ibnd) + temp_dvpsi(1:npwq)
               END IF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   deallocate (aux)
   deallocate (temp_reduction)
   if (noncolin) then
      deallocate (temp_ps2_nc_values)
   else
      deallocate (temp_ps2_values)
   end if
   deallocate (temp_dvpsi)
   IF (noncolin) THEN
      deallocate (ps2_nc)
      deallocate (ps1_nc)
      deallocate (deff_nc)
      deallocate (temp_ps1_nc)
      deallocate (temp_ps2_nc)
      deallocate (temp_alphap_nc)
      deallocate (temp_becp1_nc)
   ELSE
      deallocate (ps2)
      deallocate (ps1)
      deallocate (deff)
      deallocate (temp_ps1)
      deallocate (temp_ps2)
      deallocate (temp_alphap_k)
      deallocate (temp_becp1_k)
   END IF

   call stop_clock ('dvqpsi_us_on')
   RETURN
END SUBROUTINE dvqpsi_us_only