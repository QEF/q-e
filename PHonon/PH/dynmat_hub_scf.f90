!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_scf (irr, nu_i0, nper)
  !---------------------------------------------------------------------
  !
  !  DFPT+U: This routine adds to the dynamical matrix the scf+orth Hubbard terms.   
  !  It adds the scf derivative of the Hubbard energy (terms 1 and 2). Moreover,
  !  it adds 3 terms due to the orthogonality constraints in the USPP formalism (terms 3,4,5). 
  !  Note, the orthogonality terms 4 and 5 DO NOT follow simply from the 2nd derivative of 
  !  the Hubbard energy E_U; They stem from the coupling of the -\epsilon\dS\dlambda in 
  !  the USPP forces with the part of d\psi\d\mu projected onto the occupied manifold 
  !  (the variable dpsi stems from the scf linear system and lives in the unoccupied manifold).   
  !
  !  Terms implemented:  
  !   dyn_hub_scf (ipert, imode)  
  !  1) =  +\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi(ipert,ibnd,k+q, is)|dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)| psi(ibnd,k,is> + 
  !                                   <dpsi(ipert,ibnd,k+q, is)|S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)| psi(ibnd,k,is> +  m1<=>m2]   
  !  2) = -\sum_{I,m1,m2,is} Hubbard_U(I) * CONJG(dnsscf(m1,m2,is,I,ipert)) * dnsbare(m2,m1,is,I,imode)  
  ! 
  ! Orthogonality terms present only in USPP:
  !  3) = -\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi_orth(ipert,ibnd,k+q, is)|dqsphi(imode,I,k+q,m2)><S\phi(I,k,m1)| psi(ibnd,k,is> + 
  !                                   <dpsi_orth(ipert,ibnd,k+q, is)|S\phi(I,k+q,m2)><dmqsphi(imode,I,k,m1)| psi(ibnd,k,is>] 
  !  4) = -\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi_orth(imode,ibnd,k+q, is)|dqsphi(ipert,I,k+q,m1)><S\phi(I,k,m2)| psi(ibnd,k,is> + 
  !                                   <dpsi_orth(imode,ibnd,k+q, is)|S\phi(I,k+q,m1)><dmqsphi(ipert,I,k,m2)| psi(ibnd,k,is>]  
  !  5) = -\sum_{I,m1,m2,is} Hubbard_U(I) * [CONJG(dnsscf(m1,m2,is,I,ipert) + CONJG(dnsbare(m2,m1,is,I,ipert))] 
  !                                       * dnsorth(m1,m2,is,I,imode)  
  !
  ! Note: dnsscf includes the orthogonality part (dnsorth)
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018) 
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, Hubbard_J0
  USE ldaU_ph,       ONLY : dnsbare, dnsbare_all_modes, dnsscf, &
                            dnsorth_cart, effU
  USE lsda_mod,      ONLY : lsda, current_spin, isk, nspin
  USE modes,         ONLY : u, nmodes
  USE dynmat,        ONLY : dyn, dyn_rec, dyn_hub_scf
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE eqv,           ONLY : dpsi
  USE wvfct,         ONLY : npwx, nbnd
  USE control_lr,    ONLY : lgamma
  USE control_ph,    ONLY : rec_code_read
  USE units_ph,      ONLY : iudwf, lrdwf
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE wavefunctions, ONLY : evc
  USE klist,         ONLY : wk, lgauss, ltetra, ngk, igk_k
  USE uspp,          ONLY : okvan
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp_pools,      ONLY : inter_pool_comm       
  USE mp,            ONLY : mp_sum
  USE buffers,       ONLY : get_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: irr, nu_i0, nper
  ! the irreducible representation
  ! the initial position of the mode
  ! number of perturbations in the current irreducible representation
  !
  ! Local variables
  !
  INTEGER :: icar, jcar, na, nap, nah, ihubst1, ihubst2, nt, counter, n, l, is,    & 
             na_icar, nap_jcar, m1, m2, imode, jmode, ik, ikk, ikq, ipert, nrec1,  &
             ibnd, ig, isi, op_spin, npw, npwq
  COMPLEX(DP), ALLOCATABLE :: dyn1 (:,:), dyn_orth(:,:), term(:,:), term_aux(:,:), &
                              dpsi_orth_cart(:,:,:,:), dyn_orth_cart(:,:),         &
                              dvqhbar(:,:,:,:), dvqhbar_orth(:,:,:,:),             &
                              dvqhbar_orth_lm(:,:,:,:), aux1(:,:,:), dyn1_test(:,:) 
  REAL(DP), ALLOCATABLE :: wgg(:,:,:)
  COMPLEX(DP) :: dvi, prj, prj_orth
  COMPLEX(DP), EXTERNAL :: ZDOTC   
  LOGICAL :: lmetq0   ! .true. if q=0 for a metal
  ! 
  CALL start_clock( 'dynmat_hub_scf' )
  !
  ALLOCATE (dyn1(nper,nmodes))
  ALLOCATE (dyn_orth(nper,nmodes))
  ALLOCATE (term(nat,3))
  ALLOCATE (term_aux(nper,nmodes))
  ALLOCATE (dpsi_orth_cart(npwx,nbnd,3,nat))
  ALLOCATE (dyn_orth_cart(3*nat,3*nat))
  ALLOCATE (dvqhbar(npwx,nbnd,3,nat))
  ALLOCATE (aux1(npwx,nbnd,nmodes))
  ALLOCATE (dyn1_test(nmodes,nmodes))
  ALLOCATE (wgg(nbnd,nbnd,nksq)) 
  ALLOCATE (dvqhbar_orth(npwx,nbnd,3,nat))
  ALLOCATE (dvqhbar_orth_lm(npwx,nbnd,3,nat))
  ! 
  dyn1          = (0.d0, 0.d0)  
  dyn_orth_cart = (0.d0, 0.d0)  
  dyn1_test     = (0.d0, 0.d0)  
  ! 
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  !
  ! USPP: compute the weights as in square bracket 
  ! of Eq. (27) of A. Dal Corso PRB 64, 235118 (2001).
  ! Or see Eq. (B19) in the same reference.
  !
  IF (okvan) CALL compute_weight (wgg)
  ! 
  ! If recover=.true. recompute dnsscf for the current irr 
  ! (e.g. when convt=.true. in solve_linter but the code is 
  ! interrupted before the call of this routine)
  !
  IF (rec_code_read==10) THEN
     !WRITE(stdout,*) 'rec_code_read', rec_code_read
     CALL dnsq_scf (nper, lmetq0, nu_i0, irr, .true.)  
  ENDIF
  !
  ! We need a sum over all k points 
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     IF (lsda) THEN
        current_spin = isk(ikk)
        IF (current_spin .eq. 1) THEN
           op_spin = 2
        ELSE
           op_spin = 1
        END IF
     ELSE        
        op_spin = 1
     ENDIF
     !
     ! Read unperturbed KS wavefuctions psi(k)     
     !
     IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     !
     ! Calculate d^{na,icar}V_tilde(q) |psi(ibnd,k)> =  
     ! = |ket> Bloch vector at k+q for all cartesian coordinates (na,icar)
     !
     CALL dvqhub_barepsi_us2 (ik, dvqhbar, dvqhbar_orth, dvqhbar_orth_lm)
     !
     ! Compute terms 3 and 4.
     ! Compute the orthogonality term in cartesian coordinates.
     ! Hubbard-potential analogue to the two terms in 
     ! Eq. (42) of A. Dal Corso PRB 64, 235118 (2001).
     !
     IF (okvan)  THEN
        !
        ! Calculate dpsi_orth_cart which is the valence component 
        ! of the response KS wavefunction (non-zero only in the USPP case).
        !
        CALL dpsi_orth (ik, wgg, dpsi_orth_cart)
        !
        DO na = 1, nat
           DO icar = 1, 3
              na_icar = 3*(na-1)+icar
              DO nap = 1, nat
                 DO jcar = 1, 3
                    nap_jcar = 3*(nap-1)+jcar
                    DO ibnd = 1, nbnd
                       !
                       !          from E_Hub alone
                       prj_orth = ZDOTC (npwq, dpsi_orth_cart(:,ibnd,jcar,nap), 1,  &
                                               dvqhbar_orth(:,ibnd,icar,na), 1)   + & 
                       !          extra term     
                                  ZDOTC (npwq, dvqhbar_orth_lm(:,ibnd,jcar,nap), 1, &
                                               dpsi_orth_cart(:,ibnd,icar,na), 1)
                       ! 
                       CALL mp_sum(prj_orth, intra_bgrp_comm) 
                       ! 
                       dyn_orth_cart (nap_jcar,na_icar) = dyn_orth_cart (nap_jcar,na_icar)  &
                                                          - wk(ikk) * prj_orth 
                       !
                    ENDDO 
                 ENDDO 
              ENDDO 
           ENDDO
        ENDDO
        ! 
     ENDIF
     !
     ! Compute term 1.
     ! Transform dvqhbar from the Cartesian to the pattern basis u. 
     ! aux1 is dvqhbar in the pattern basis.
     ! 
     aux1 = (0.d0, 0.d0)
     !
     DO imode = 1, nmodes 
        DO na = 1, nat
           DO icar = 1, 3 
              na_icar = 3*(na-1)+icar
              DO ibnd = 1, nbnd           
                 DO ig = 1, npwq
                    aux1(ig,ibnd,imode) = aux1(ig,ibnd,imode) + & 
                         dvqhbar(ig,ibnd,icar,na) * u(na_icar,imode)
                 ENDDO 
              ENDDO 
           ENDDO 
        ENDDO 
     ENDDO 
     !
     DO ipert = 1, nper
        !                
        nrec1 = (ipert-1)*nksq+ik
        !
        ! Read the response KS wave functions dpsi from file 
        !
        CALL get_buffer (dpsi, lrdwf, iudwf, nrec1)
        ! 
        DO imode = 1, nmodes 
           ! 
           ! Calculate prj = <d^(ipert)dpsi_k+q | d^{imode}V_tilde(q) | psi(ibnd,k)>
           !
           DO ibnd = 1, nbnd
              !
              prj = ZDOTC (npwq, dpsi(:,ibnd), 1, aux1(:,ibnd,imode), 1)
                          
              CALL mp_sum (prj, intra_bgrp_comm) 
              !
              ! dyn1 is a sum over all the k and ibnd
              !
              dyn1(ipert,imode) =  dyn1(ipert,imode) + wk(ikk)*prj 
              !
           ENDDO 
           !
        ENDDO 
        ! 
     ENDDO 
     !
  ENDDO ! ik
  !
  CALL mp_sum (dyn1, inter_pool_comm) 
  CALL mp_sum (dyn_orth_cart, inter_pool_comm)
  !
  ! Compute terms 2 and 5 (and the equivalent for Hubbard_J0)
  !
  DO ipert = 1, nper
     !
     term = (0.d0, 0.d0)   
     !    
     DO na = 1, nat
        !
        DO icar = 1, 3
           !
           DO nah = 1, nat
              !
              nt = ityp(nah)
              !
              ! For effU = Hubbard_U - Hubbard_J0
              !
              IF (is_hubbard(nt)) THEN
                 !
                 dvi = (0.d0, 0.d0)
                 !           
                 DO is = 1, nspin
                    !
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       !
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          ! 
                          ! Term 2
                          !
                          dvi = dvi + CONJG(dnsscf(m1,m2,is,nah,ipert)) * &
                                      dnsbare(m2,m1,is,nah,icar,na)   
                          ! 
                          ! Term 5
                          !
                          IF (okvan) dvi = dvi + CONJG( dnsscf(m1,m2,is,nah,ipert) + &
                                                        dnsbare_all_modes(m1,m2,is,nah,nu_i0+ipert) ) &
                                                      * dnsorth_cart(m1,m2,is,nah,icar,na)   
                          !
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDDO
                 !                    
                 ! term is implicitly defined for each ipert
                 !
                 term(na,icar) = term(na,icar) - dvi*effU(nt)  
                 ! 
              ENDIF
              !
              ! For Hubbard_J0
              !
              IF (Hubbard_J0(nt).NE.0.d0) THEN
                 !
                 dvi = (0.d0, 0.d0)   
                 !        
                 DO is = 1, nspin
                    !                    
                    IF ((is==1).AND.(nspin==2)) THEN
                       isi = 2
                    ELSE 
                       isi = 1
                    ENDIF
                    !
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       !
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          !
                          ! Term 2
                          !
                          dvi = dvi - CONJG(dnsscf(m1,m2,isi,nah,ipert)) * & ! sign change
                                      dnsbare(m2,m1,is,nah,icar,na)   
                          !
                          ! Term 5
                          !
                          IF (okvan) dvi = dvi - CONJG( dnsscf(m1,m2,isi,nah,ipert) + & ! sign change
                                                        dnsbare_all_modes(m1,m2,isi,nah,nu_i0+ipert) ) &
                                                      * dnsorth_cart(m1,m2,is,nah,icar,na)  
                          ! 
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDDO
                 !                    
                 ! term is implicitely defined for each ipert
                 !
                 term(na,icar) = term(na,icar) - dvi*Hubbard_J0(nt)  
                 !
              ENDIF
              !
           ENDDO ! nah
           !
        ENDDO ! icar
        !
     ENDDO ! na
     !
     ! transform term from the Cartesian to the pattern basis u. 
     ! The result is stored in term_aux.
     ! 
     DO imode = 1, nmodes
        term_aux(ipert,imode) = (0.d0, 0.d0)
        DO na = 1, nat
           DO icar = 1, 3
              na_icar = 3*(na-1)+icar
              term_aux(ipert,imode) = term_aux(ipert,imode) + &
                                      term(na,icar) * u(na_icar,imode)         
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO ! ipert
  !
  ! Factor of 2 taking into account that
  ! if nspin=1 the spin sum above does not occur   
  !
  IF (nspin==1) term_aux = 2.d0 * term_aux
  !
  ! Add the contribution to the dyn_hub_scf matrix
  !
  dyn1 = dyn1 + term_aux    
  !  
  dyn_orth = (0.d0, 0.d0)
  !
  ! Adding the orthogonality terms 3 and 4, 
  ! which were computed above (and summed up over k), to dyn1.  
  !
  IF (okvan) THEN
     DO ipert = 1, nper
        DO imode = 1, nmodes
           DO na_icar = 1, 3*nat
              DO nap_jcar = 1, 3*nat 
                 dyn_orth(ipert,imode) = dyn_orth(ipert,imode) &
                                         + u(na_icar,imode)    &  
                                         * dyn_orth_cart(nap_jcar,na_icar) &
                                         * CONJG(u(nap_jcar,nu_i0+ipert)) 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     dyn1 = dyn1 + dyn_orth
  ENDIF
  ! 
  DO imode = 1, nmodes
     DO ipert = 1, nper
        !
        jmode = nu_i0 + ipert
        !
        dyn_hub_scf(jmode,imode) = dyn1(ipert,imode)
        !
        dyn_rec(jmode,imode) = dyn_rec(jmode,imode) + dyn1(ipert,imode)
        ! 
        dyn1_test(jmode,imode) = dyn1 (ipert, imode)
        !
     ENDDO
  ENDDO
  !
  ! Write the symmetrized dyn_hub_scf (upgraded for the new ipert) and total 
  ! dynamical matrix dyn in  cartesian coordinates.
  ! dyn1 is instead the contribution of the particular irrep only.
  !
  IF (iverbosity==1) THEN
     CALL tra_write_matrix_no_sym('dyn1 NOT SYMMETRIZED',dyn1_test,nat)
     CALL tra_write_matrix('dyn1 SYMMETRIZED',dyn1_test,u,nat)
     CALL tra_write_matrix('dyn',dyn,u,nat)
  ENDIF
  !
  ! Add Hubbard terms in the dynamical matrix to the total
  ! dynamical matrix.
  !
  DO imode = 1, nmodes
     DO ipert = 1, nper
        jmode = nu_i0 + ipert
        dyn(jmode,imode) = dyn(jmode,imode) + dyn1(ipert,imode)
     ENDDO
  ENDDO
  !
  DEALLOCATE (dyn1)
  DEALLOCATE (dyn_orth)
  DEALLOCATE (term)
  DEALLOCATE (term_aux)
  DEALLOCATE (dpsi_orth_cart)
  DEALLOCATE (dyn_orth_cart)
  DEALLOCATE (dvqhbar)
  DEALLOCATE (aux1)
  DEALLOCATE (dyn1_test)
  DEALLOCATE (wgg)
  DEALLOCATE (dvqhbar_orth)
  DEALLOCATE (dvqhbar_orth_lm)
  !   
  CALL stop_clock ( 'dynmat_hub_scf' )
  !
  RETURN
  !
END SUBROUTINE dynmat_hub_scf
