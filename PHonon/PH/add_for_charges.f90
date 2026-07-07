!
! Copyright (C) 2001-2007 Quantum ESPRESSO PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE add_for_charges (ik, uact)
  !----------===============-----------------------------------------------
  !! Applies \(\frac{dS}{du}\) to \(\text{dpsi}\) and accumulates the result into \(\text{dvpsi}\).
  !! On input \(\text{dpsi}\) is expected to contain \(P_c x |\psi_{ik}\rangle\), as computed by
  !! \(\texttt{dvpsi\_e}\) and stored to disk via the \(\texttt{iucom}\) buffer.
  !

  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE cell_base, ONLY : tpiba
  USE gvect, ONLY : g
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE klist, ONLY : xk, ngk, igk_k
  USE uspp, ONLY : nkb, qq_nt, qq_so, vkb, ofsbeta
  USE wvfct, ONLY : npwx, nbnd
  USE becmod, ONLY: calbec, bec_type, allocate_bec_type, deallocate_bec_type
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE uspp_param, ONLY: nh, nhm
  USE eqv, ONLY : dvpsi, dpsi
  USE control_lr, ONLY : lgamma
  USE qpoint,  ONLY : ikks
  IMPLICIT NONE
  !
  !   The dummy variables
  !

  INTEGER :: ik
  !! input: the k point
  INTEGER :: mode
  ! input: the actual perturbation
  COMPLEX(DP) :: uact (3 * nat)
  !! input: the pattern of displacements
  !
  !   And the local variables
  !

  INTEGER :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, js, ijs, npw
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

  REAL(DP), PARAMETER :: eps = 1.d-12

  COMPLEX(DP), ALLOCATABLE :: ps1 (:,:), ps2 (:,:,:), aux (:)
  COMPLEX(DP), ALLOCATABLE :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! temporary arrays for optimization
  COMPLEX(DP) :: temp_ps1, temp_ps2
  COMPLEX(DP) :: temp_ps1_nc(npol), temp_ps2_nc(npol)
  COMPLEX(DP) :: ps2_col(nbnd)
  ! small buffers for loop optimization (allocated once with max size)
  COMPLEX(DP) :: alphapp_buf_nc(nhm, npol), bedp_buf_nc(nhm, npol)
  COMPLEX(DP) :: alphapp_buf_k(nhm), bedp_buf_k(nhm)
  COMPLEX(DP) :: qq_so_buf(nhm, npol*npol), qq_nt_buf(nhm)
  ! the scalar product
  ! the scalar product
  ! a mesh space for psi
  TYPE(bec_type) :: bedp, alphapp(3)
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)

  LOGICAL :: ok
  ! used to save time

  ALLOCATE (aux ( npwx))
  ALLOCATE (aux1( npwx*npol, nbnd))
  CALL allocate_bec_type(nkb,nbnd,bedp)
  DO ipol=1,3
     CALL allocate_bec_type(nkb,nbnd,alphapp(ipol))
  ENDDO
  IF (noncolin) THEN
     ALLOCATE (ps1_nc ( nkb, npol, nbnd))
     ALLOCATE (ps2_nc ( nkb, npol, nbnd , 3))
  ELSE
     ALLOCATE (ps1 ( nkb , nbnd))
     ALLOCATE (ps2 ( nkb , nbnd , 3))
  ENDIF
  IF (lgamma) THEN
     ikk = ikks(ik)
     ikq = ikk
     npw =ngk(ikk)
  ELSE
     CALL infomsg ('add_for_charges', 'called for lgamma .eq. false')
  ENDIF
  IF (lsda) current_spin = isk (ikk)
  !
  !   we first compute the coefficients of the vectors
  !
  IF (noncolin) THEN
     ps1_nc   = (0.d0, 0.d0)
     ps2_nc   = (0.d0, 0.d0)
     bedp%nc = (0.d0,0.d0)
     DO ipol=1,3
        alphapp(ipol)%nc = (0.d0,0.d0)
     END DO
  ELSE
     ps1   = (0.d0, 0.d0)
     ps2   = (0.d0, 0.d0)
     bedp%k = (0.d0,0.d0)
     DO ipol=1,3
        alphapp(ipol)%k = (0.d0,0.d0)
     END DO
  ENDIF
  aux1  = (0.d0, 0.d0)

  !
  ! first we calculate the products of the beta functions with dpsi
  !
  CALL calbec (npw, vkb, dpsi, bedp)
  DO ipol = 1, 3
     aux1=(0.d0,0.d0)
     DO ibnd = 1, nbnd
        DO ig = 1, npw
           aux1 (ig, ibnd) = dpsi(ig,ibnd) *           &
                tpiba * (0.d0,1.d0) *                  &
                ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
        ENDDO
        IF (noncolin) THEN
           DO ig = 1, npw
              aux1 (ig+npwx, ibnd) = dpsi(ig+npwx,ibnd) *           &
                   tpiba * (0.d0,1.d0) *                  &
                  ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
           ENDDO
        ENDIF
     ENDDO
     CALL calbec ( npw, vkb, aux1, alphapp(ipol) )
  ENDDO
  DEALLOCATE (aux1)

  DO na = 1, nat
     nt = ityp (na)
     ijkb0=ofsbeta(na)
     mu = 3 * (na - 1)
     IF ( ABS (uact (mu + 1) ) + &
          ABS (uact (mu + 2) ) + &
          ABS (uact (mu + 3) ) > eps) THEN
        DO ih = 1, nh (nt)
           ikb = ijkb0 + ih
           DO ipol = 1, 3
              DO ibnd = 1, nbnd
                 ! Populate small buffers for this iteration
                 qq_nt_buf(1:nh(nt)) = qq_nt(ih, 1:nh(nt), nt)

                 IF (noncolin) THEN
                    DO jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphapp_buf_nc(jh, 1:npol) = alphapp(ipol)%nc(jkb, 1:npol, ibnd)
                       bedp_buf_nc(jh, 1:npol) = bedp%nc(jkb, 1:npol, ibnd)
                    ENDDO
                    IF (lspinorb) THEN
                       DO jh = 1, nh(nt)
                          qq_so_buf(jh, 1:npol*npol) = qq_so(ih, jh, 1:npol*npol, nt)
                       ENDDO
                    ENDIF
                 ELSE
                    DO jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphapp_buf_k(jh) = alphapp(ipol)%k(jkb, ibnd)
                       bedp_buf_k(jh) = bedp%k(jkb, ibnd)
                    ENDDO
                 ENDIF

                 ! Initialize temp arrays
                 IF (noncolin) THEN
                    temp_ps1_nc = (0.d0, 0.d0)
                    temp_ps2_nc = (0.d0, 0.d0)
                 ELSE
                    temp_ps1 = (0.d0, 0.d0)
                    temp_ps2 = (0.d0, 0.d0)
                 ENDIF

                 DO jh = 1, nh (nt)
                    IF (noncolin) THEN
                       IF (lspinorb) THEN
                          ijs=0
                          DO is=1,npol
                             DO js=1,npol
                                ijs=ijs+1
                                temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                (qq_so_buf(jh,ijs) *              &
                                alphapp_buf_nc(jh,js))*         &
                                uact (mu + ipol)
                                temp_ps2_nc(is) = temp_ps2_nc(is) + &
                                (qq_so_buf(jh,ijs) *              &
                                 bedp_buf_nc(jh,js)) *            &
                                 uact (mu + ipol)
                             ENDDO
                          ENDDO
                       ELSE
                          DO is=1,npol
                             temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                 qq_nt_buf(jh) *                     &
                                 alphapp_buf_nc(jh,is) *     &
                                 uact (mu + ipol)
                             temp_ps2_nc(is) = temp_ps2_nc(is) + &
                                 qq_nt_buf(jh) *                  &
                                 bedp_buf_nc(jh,is) *             &
                                 uact (mu + ipol)
                          END DO
                       ENDIF
                    ELSE
                       temp_ps1 = temp_ps1 + qq_nt_buf(jh)*alphapp_buf_k(jh)* &
                            uact (mu + ipol)
                       temp_ps2 = temp_ps2 + qq_nt_buf(jh) * &
                             bedp_buf_k(jh) * uact (mu + ipol)
                    ENDIF
                 ENDDO

                 ! ps1(ikb,ibnd)      = sum_{jh,ipol} q_{ih,jh} <d(beta_{jh})/d(tau)*u | dpsi_ibnd>
                 ! ps2(ikb,ibnd,ipol) = sum_{jh}      q_{ih,jh} <beta_{jh} | dpsi_ibnd> * uact(ipol)
                 IF (noncolin) THEN
                    DO is=1,npol
                       ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + temp_ps1_nc(is)
                       ps2_nc(ikb,is,ibnd,ipol) = temp_ps2_nc(is)
                    ENDDO
                 ELSE
                    ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + temp_ps1
                    ps2 (ikb, ibnd, ipol) =  temp_ps2
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  CALL deallocate_bec_type(bedp)
  DO ipol=1,3
     CALL deallocate_bec_type(alphapp(ipol))
  END DO
  !
  !  ps1(ikb,ibnd)   = sum_{jh,ipol} q_{ih,jh} <d(beta_{jh})/d(tau)*u | dpsi_ibnd>
  !  adds sum_{ih,jh} |beta_{ih}> q_{ih,jh} <d(beta_{jh})/d(tau)*u | dpsi> to dvpsi
  !
  IF (nkb.GT.0) THEN
     IF (noncolin) THEN
        CALL zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     ELSE
        CALL zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     ENDIF
  ENDIF
  IF (noncolin) THEN
     DEALLOCATE (ps1_nc)
  ELSE
     DEALLOCATE (ps1)
  END IF
  !
  !  ps2(ikb,ibnd,ipol) = sum_{jh} q_{ih,jh} <beta_{jh} | dpsi_ibnd> * uact(ipol)
  !  adds sum_{ih,jh} |d(beta_{ih})/d(tau)*u> q_{ih,jh} <beta_{jh} | dpsi> to dvpsi
  !
  DO ikb = 1, nkb
     DO ipol = 1, 3
        IF (noncolin) THEN
           ok = ANY(ABS(ps2_nc(ikb, 1, 1:nbnd, ipol)) > eps) .OR. &
                ANY(ABS(ps2_nc(ikb, 2, 1:nbnd, ipol)) > eps)
        ELSE
           ok = ANY(ABS(ps2(ikb, 1:nbnd, ipol)) > eps)
        ENDIF
        IF (ok) THEN
           DO ig = 1, npw
              igg = igk_k (ig,ikq)
              aux (ig) = vkb(ig, ikb) * (0.d0,-1.d0) * tpiba * (xk(ipol, ikq) + g(ipol, igg) )
           ENDDO
           IF (noncolin) THEN
              ps2_col(1:nbnd) = ps2_nc(ikb, 1, 1:nbnd, ipol)
              CALL zgemm('N','T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, ps2_col, nbnd, (1.d0,0.d0), dvpsi(1,1), npwx*npol)
              ps2_col(1:nbnd) = ps2_nc(ikb, 2, 1:nbnd, ipol)
              CALL zgemm('N','T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, ps2_col, nbnd, (1.d0,0.d0), dvpsi(npwx+1,1), npwx*npol)
           ELSE
              ps2_col(1:nbnd) = ps2(ikb, 1:nbnd, ipol)
              CALL zgemm('N','T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, ps2_col, nbnd, (1.d0,0.d0), dvpsi, npwx*npol)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE (aux)
  IF (noncolin) THEN
     DEALLOCATE (ps2_nc)
  ELSE
     DEALLOCATE (ps2)
  END IF
!
!    Now dvpsi contains dS/du applied to dpsi, i.e. (dS/du) P_c x |psi>
!
  RETURN
END SUBROUTINE add_for_charges
