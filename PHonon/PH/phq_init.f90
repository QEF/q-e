!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE phq_init()
  !----------------------------------------------------------------------------
  !! This subroutine computes the quantities necessary to describe the
  !! local and nonlocal pseudopotential in the phononq program.
  !! In detail it computes:  
  !! 0)  initialize the structure factors;  
  !! a)  compute rhocore for each atomic-type if needed for nlcc;  
  !! b)  the local potential at q+G-G'. Needed for the second
  !!     second part of the dynamical matrix;  
  !! c)  the D coefficients for the US pseudopotential or the E_l parame
  !!     of the KB pseudo. In the US case it prepares also the integrals
  !!     qrad and qradq which are needed for computing \(Q_{nm}(G)\) and
  !!     \(Q_{nm}(q+G)\);  
  !! d)  the functions \(\text{vkb}(k+G)\) needed for the part of the 
  !!     dynamical matrix independent of deltapsi;  
  !! e)  The becp functions for the k points;  
  !! e') The derivative of the becp term with respect to a displacement;  
  !! f)  The functions \(\text{vkb}(k+q+G)\), needed for the linear system
  !!     and the second part of the dynamical matrix.
  !
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : bg, tpiba
  USE ions_base,            ONLY : nat, ityp, tau
  USE becmod,               ONLY : calbec, becp, becupdate, bec_type, &
                                   allocate_bec_type_acc, deallocate_bec_type_acc
  USE control_flags,        ONLY : offload_type
  USE constants,            ONLY : eps8, tpi
  USE gvect,                ONLY : g
  USE klist,                ONLY : xk, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : npwx, nbnd
  USE gvecw,                ONLY : gcutw
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc_d
#endif
  USE noncollin_module,     ONLY : noncolin, domag, npol, lspinorb
  USE uspp,                 ONLY : okvan, vkb, nlcc_any, nkb
  USE phus,                 ONLY : alphap
  USE nlcc_ph,              ONLY : drc
  USE control_ph,           ONLY : trans, zue, epsil, all_done
  USE units_lr,             ONLY : lrwfc, iuwfc
  USE mp,                   ONLY : mp_sum
  USE acfdtest,             ONLY : acfdt_is_active, acfdt_num_der
  USE el_phon,              ONLY : elph_mat, iunwfcwann, npwq_refolded, &
                                   kpq,g_kpq,igqg,xk_gamma, lrwfcr
  USE wannier_gw,           ONLY : l_head
  USE lrus,                 ONLY : becp1, dpqq, dpqq_so
  USE qpoint,               ONLY : xq, nksq, eigqts, ikks, ikqs
  USE qpoint_aux,           ONLY : becpt, alphapt, ikmks
  USE eqv,                  ONLY : evq
  USE control_lr,           ONLY : nbnd_occ, lgamma
  USE ldaU,                 ONLY : lda_plus_u
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: nt, ik, ikq, ipol, ibnd, ikk, na, ig, irr, imode0, itmp
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  INTEGER :: ikqg         !for the case elph_mat=.true.
  INTEGER :: npw, npwq
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:), tevc(:,:)
    ! used to compute alphap
#if defined(__CUDA)
  TYPE(bec_type) :: bectmp
    ! temporary buffer to work with offload of arrays of derived types
#endif
  !
  !
  IF (all_done) RETURN
  !
  CALL start_clock( 'phq_init' )
  !
#if defined(__CUDA)
  Call allocate_bec_type_acc ( nkb, nbnd, bectmp )
#endif
  !
  DO na = 1, nat
     !
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     !
     eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
     !
  END DO
  !
  ! ... a) compute rhocore for each atomic-type if needed for nlcc
  !
  IF ( nlcc_any ) CALL set_drhoc( xq, drc )
  !
  ! ... b) the fourier components of the local potential at q+G
  !
  CALL init_vlocq ( xq ) 
  !
  ! only for electron-phonon coupling with wannier functions
  ! 
  if(elph_mat) then
    ALLOCATE(kpq(nksq),g_kpq(3,nksq),igqg(nksq))
    ALLOCATE (xk_gamma(3,nksq))
 
    do ik=1,nksq
      xk_gamma(1:3,ik)=xk(1:3,ikks(ik))
    enddo
      !
      !first of all I identify q' in the list of xk such that
      !   (i) q' is in the set of xk
      !   (ii) k+q'+G=k+q
      !  and G is a G vector depending on k and q.
      !
    call get_equivalent_kpq(xk_gamma,xq,kpq,g_kpq,igqg)

  endif
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )
  !
  IF (noncolin.AND.domag) THEN
          ALLOCATE(tevc(npwx*npol,nbnd))
          !$acc enter data create(tevc(1:npwx*npol,1:nbnd))
  ENDIF
  !
  !$acc data copyin(xk) create(aux1( 1:npwx*npol, 1:nbnd )) 
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     IF ( lsda ) current_spin = isk( ikk )
     !
     IF ( .NOT. lgamma ) THEN
        !
        IF ( ABS( xq(1) - ( xk(1,ikq) - xk(1,ikk) ) ) > eps8 .OR. &
             ABS( xq(2) - ( xk(2,ikq) - xk(2,ikk) ) ) > eps8 .OR. &
             ABS( xq(3) - ( xk(3,ikq) - xk(3,ikk) ) ) > eps8 ) THEN
           WRITE( stdout,'(/,5x,"k points #",i6," and ", &
                  & i6,5x," total number ",i6)') ikk, ikq, nksq
           WRITE( stdout, '(  5x,"Expected q ",3f10.7)')(xq(ipol), ipol=1,3)
           WRITE( stdout, '(  5x,"Found      ",3f10.7)')((xk(ipol,ikq) &
                                                -xk(ipol,ikk)), ipol = 1, 3)
           CALL errore( 'phq_init', 'wrong order of k points', 1 )
        END IF
        !
     END IF
     !
     ! ... d) The functions vkb(k+G)
     !
     CALL init_us_2( npw, igk_k(1,ikk), xk(1,ikk), vkb, .true. )
     !
     ! ... read the wavefunctions at k
     !
     if(elph_mat) then
        call read_wfc_rspace_and_fwfft( evc, ik, lrwfcr, iunwfcwann, npw, igk_k(1,ikk) )
        !       CALL davcio (evc, lrwfc, iunwfcwann, ik, - 1)
     else
        CALL get_buffer( evc, lrwfc, iuwfc, ikk )
        IF (noncolin.AND.domag) THEN
           CALL get_buffer( tevc, lrwfc, iuwfc, ikmks(ik) )
#if defined(__CUDA)
           !$acc update device(tevc)
           Call calbec ( offload_type, npw, vkb, tevc, bectmp )
           Call becupdate( offload_type, becpt, ik, nksq, bectmp )
#else
           Call calbec ( offload_type, npw, vkb, tevc, becpt(ik) )
#endif
        ENDIF
     endif
     !
     ! ... e) we compute the becp terms which are used in the rest of
     ! ...    the code
     !
#if defined(__CUDA)
     evc_d = evc
     !$acc data present_or_copyin(evc)
     Call calbec( offload_type, npw, vkb, evc, bectmp )
     !$acc end data
     Call becupdate( offload_type, becp1, ik, nksq, bectmp ) 
#else
     Call calbec( offload_type, npw, vkb, evc, becp1(ik) )
#endif
     !
     ! ... e') we compute the derivative of the becp term with respect to an
     !         atomic displacement
     !
     DO ipol = 1, 3
#if defined(__CUDA)
        !$acc parallel loop collapse(2)
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              aux1(ig,ibnd)=(0.d0,0.d0)
           END DO
        END DO
#else
        aux1=(0.d0,0.d0)
#endif
        !$acc parallel loop collapse(2) 
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              itmp = igk_k(ig,ikk)
#if defined(__CUDA)
              aux1(ig,ibnd) = evc_d(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * &
                   ( xk(ipol,ikk) + g(ipol,itmp) )
#else
              aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * &
                   ( xk(ipol,ikk) + g(ipol,itmp) )
#endif
           END DO
        END DO
        IF (noncolin) THEN
           !$acc parallel loop collapse(2) 
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 itmp = igk_k(ig,ikk)
#if defined(__CUDA)
                 aux1(ig+npwx,ibnd)=evc_d(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*&
                      ( xk(ipol,ikk) + g(ipol,itmp) )
#else
                 aux1(ig+npwx,ibnd)=evc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*&
                      ( xk(ipol,ikk) + g(ipol,itmp) )
#endif
              END DO
           END DO
        END IF
#if defined(__CUDA)
        Call calbec ( offload_type, npw, vkb, aux1, bectmp )
        Call becupdate( offload_type, alphap, ipol, 3, ik, nksq, bectmp )
#else
        Call calbec ( offload_type, npw, vkb, aux1, alphap(ipol,ik) )
#endif
     END DO
     !
     IF (noncolin.AND.domag) THEN
        DO ipol = 1, 3
#if defined(__CUDA)
           !$acc parallel loop collapse(2)
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 aux1(ig,ibnd)=(0.d0,0.d0)
              END DO
           END DO
#else
           aux1=(0.d0,0.d0)
#endif
           !$acc parallel loop collapse(2) 
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 itmp = igk_k(ig,ikk)
                 aux1(ig,ibnd) = tevc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * &
                      ( xk(ipol,ikk) + g(ipol,itmp) )
              END DO
           END DO
           IF (noncolin) THEN
              !$acc parallel loop collapse(2) 
              DO ibnd = 1, nbnd
                 DO ig = 1, npw
                    itmp = igk_k(ig,ikk)
                    aux1(ig+npwx,ibnd)=tevc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*&
                         ( xk(ipol,ikk) + g(ipol,itmp) )
                 END DO
              END DO
           END IF
#if defined(__CUDA)
           Call calbec( offload_type, npw, vkb, aux1, bectmp )
           Call becupdate( offload_type, alphapt, ipol, 3, ik, nksq, bectmp )
#else
           Call calbec( offload_type, npw, vkb, aux1, alphapt(ipol,ik) )
#endif
         END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
     IF (acfdt_is_active) THEN
        ! ACFDT -test always read calculated wcf from non_scf calculation
        IF(acfdt_num_der) then 
           CALL get_buffer( evq, lrwfc, iuwfc, ikq )
        ELSE
           IF ( .NOT. lgamma ) &
                CALL get_buffer( evq, lrwfc, iuwfc, ikq )
        ENDIF
     ELSE
        ! this is the standard treatment
        IF ( .NOT. lgamma .and..not. elph_mat )then 
           CALL get_buffer( evq, lrwfc, iuwfc, ikq )
        ELSEIF(.NOT. lgamma .and. elph_mat) then
           !
           ! I read the wavefunction in real space and fwfft it
           !
           ikqg = kpq(ik)
           call read_wfc_rspace_and_fwfft( evq, ikqg, lrwfcr, iunwfcwann, npwq, &
                igk_k(1,ikq) )
           !        CALL davcio (evq, lrwfc, iunwfcwann, ikqg, - 1)
           call calculate_and_apply_phase(ik, ikqg, igqg, &
                npwq_refolded, g_kpq, xk_gamma, evq, .false.)
        ENDIF
     ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! END OF ACFDT TEST !!!!!!!!!!!!!!!!
     !

  END DO
  !$acc end data
  !
  DEALLOCATE( aux1 )
  IF (noncolin.AND.domag) THEN
          !$acc exit data delete(tevc)
          DEALLOCATE( tevc )
  ENDIF
  !
  CALL dvanqq()
  CALL drho()
  !
  IF ( ( epsil .OR. zue .OR. l_head) .AND. okvan ) THEN
     CALL compute_qdipol(dpqq)
     IF (lspinorb) CALL compute_qdipol_so(dpqq, dpqq_so)
     CALL qdipol_cryst()
  END IF
  !
  ! DFPT+U
  ! 
  IF (lda_plus_u)  THEN 
     !
     ! Calculate and write to file the atomic orbitals 
     ! \phi and S\phi at k and k+q
     ! Note: the array becp will be temporarily used
     ! in the routine lr_orthoUwfc.
     !
     CALL deallocate_bec_type_acc(becp)
     CALL lr_orthoUwfc (.TRUE.)
     CALL allocate_bec_type_acc(nkb,nbnd,becp)   
     !
     ! Calculate dnsbare, i.e. the bare variation of ns, 
     ! for all cartesian coordinates
     !
     CALL dnsq_bare()
     !
     ! Calculate the orthogonality term in the USPP case
     !
     IF (okvan) CALL dnsq_orth()
     !
  ENDIF
  !
  IF ( trans ) CALL dynmat0_new()
  !
#if defined(__CUDA)
  Call deallocate_bec_type_acc ( bectmp )
#endif
  !
  CALL stop_clock( 'phq_init' )
  !
  RETURN
  !
END SUBROUTINE phq_init
