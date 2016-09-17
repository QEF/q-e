!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_sm1_psi (recalculate, ik, lda, n, m, psi, spsi)
  !----------------------------------------------------------------------------
  !
  ! This subroutine applies the S^{-1} matrix to m wavefunctions psi
  ! and puts the results in spsi.
  ! See Eq.(13) in B. Walker and R. Gebauer, JCP 127, 164106 (2007).
  ! Requires the products of psi with all beta functions
  ! in array becp(nkb,m) (calculated in h_psi or by ccalbec).
  !
  ! INPUT: recalculate   Decides if the overlap of beta 
  !                      functions is recalculated or not.
  !                      This is needed e.g. if ions are moved 
  !                      and the overlap changes accordingly.
  !        ik            k point under consideration
  !        lda           leading dimension of arrays psi, spsi
  !        n             true dimension of psi, spsi
  !        m             number of states psi
  !        psi           the wavefunction to which the S^{-1} 
  !                      matrix is applied
  ! OUTPUT: spsi = S^{-1}*psi
  !
  ! Original routine written by R. Gebauer
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE uspp,             ONLY : okvan, vkb, nkb, qq
  USE uspp_param,       ONLY : nh, upf
  USE ions_base,        ONLY : ityp,nat,ntyp=>nsp
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : noncolin, npol
  USE matrix_inversion
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in)      :: recalculate
  INTEGER, INTENT(in)      :: lda, n, m, ik
  COMPLEX(DP), INTENT(in)  :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(out) :: spsi(lda*npol,m)
  !
  LOGICAL :: recalc
  !
  CALL start_clock( 'lr_sm1_psi' )
  !
  recalc = recalculate
  !
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSEIF (noncolin) THEN
     CALL errore( 'lr_sm1_psi', 'Noncollinear case is not implemented', 1 )
  ELSE
     CALL sm1_psi_k()
  ENDIF
  !
  CALL stop_clock( 'lr_sm1_psi' )
  !
  RETURN
  !
CONTAINS
  ! 
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    ! gamma_only version
    !
    ! Note: 1) vkb must be computed before calling this routine
    !       2) the array bbg must be deallocated somewhere
    !          outside of this routine.
    !
    USE becmod,   ONLY : bec_type,becp,calbec
    USE realus,   ONLY : real_space, invfft_orbital_gamma, initialisation_level, &
                         fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                         v_loc_psir, s_psir_gamma, real_space_debug
    USE lrus,     ONLY : bbg
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    ! counters
    REAL(DP), ALLOCATABLE :: ps(:,:)
    LOGICAL, SAVE :: first_entry = .true.
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda * npol * m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B but
    ! use the ones which were already calculated and saved (if recalc=.false.).
    !
    IF (first_entry) THEN
       !
       IF (allocated(bbg)) DEALLOCATE(bbg)
       first_entry = .false.
       recalc = .true.
       !
    ENDIF
    !
    IF (.not.allocated(bbg)) recalc = .true.
    IF (recalc .and. allocated(bbg)) DEALLOCATE(bbg)
    !
    IF (recalc) THEN
       !
       ALLOCATE(bbg(nkb,nkb))
       bbg = 0.d0
       !
       ! The beta-functions vkb must have beed calculated elsewhere. 
       ! So there is no need to call the routine init_us_2.
       !
       ! Calculate the coefficients B_ij defined by Eq.(15).
       ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
       !
       CALL calbec (n, vkb, vkb, bbg(:,:), nkb)
       !
       ALLOCATE(ps(nkb,nkb))
       !
       ps(:,:) = 0.d0
       !
       ! Calculate the product of q_nm and B_ij of Eq.(16).
       !
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + &
                               & qq(ih,jh,nt) * bbg(jkb,ii)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO
       !
       ! Add identity to q_nm * B_ij [see Eq.(16)].
       ! ps = (1 + q*B)
       !
       DO ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + 1.d0
       ENDDO
       !
       ! Invert matrix: (1 + q*B)^{-1}
       !
       CALL invmat( nkb, ps )
       !
       ! Use the array bbg as a workspace in order to save memory
       !
       bbg(:,:) = 0.d0
       !
       ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
       !
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            bbg(ii,jkb) = bbg(ii,jkb) &
                                    & - ps(ii,ikb) * qq(ih,jh,nt)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO
       DEALLOCATE(ps)
    ENDIF
    !
    ! Compute the product of the beta-functions vkb with the functions psi,
    ! and put the result in becp.
    ! becp(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    IF (real_space_debug>3) THEN 
       !
       DO ibnd=1,m,2
          CALL invfft_orbital_gamma(psi,ibnd,m)
          CALL calbec_rs_gamma(ibnd,m,becp%r)
       ENDDO
       !
    ELSE
       CALL calbec(n,vkb,psi,becp,m)
    ENDIF
    !
    ! Use the array ps as a workspace
    ALLOCATE(ps(nkb,m))
    ps(:,:) = 0.D0
    !
    ! Now let us apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : ps = lambda * <beta|psi>
    ! Here, lambda = bbg, and <beta|psi> = becp%r
    !
    call DGEMM( 'N','N',nkb,m,nkb,1.d0,bbg(:,:),nkb,becp%r,nkb,0.d0,ps,nkb)
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta>
    !
    call DGEMM('N','N',2*n,m,nkb,1.d0,vkb,2*lda,ps,nkb,1.d0,spsi,2*lda)
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_gamma
  
  SUBROUTINE sm1_psi_k()
    !-----------------------------------------------------------------------
    !
    ! k-points version
    ! Note: the array bbk must be deallocated somewhere
    ! outside of this routine.
    !
    USE becmod,   ONLY : bec_type,becp,calbec
    USE klist,    ONLY : nks, xk, ngk, igk_k
    USE lrus,     ONLY : bbk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ik1
    ! counters
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*npol*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B 
    ! (they are already allocated) but use the ones which were already calculated 
    ! and saved (if recalc=.false.).
    !
    IF (.NOT.ALLOCATED(bbk)) recalc = .true.
    IF (recalc .AND. ALLOCATED(bbk)) DEALLOCATE(bbk)
    !
    ! If recalc=.true. we (re)calculate the coefficients lambda_nm defined by Eq.(16).
    !
    IF (recalc) THEN
       !
       ALLOCATE(bbk(nkb,nkb,nks))
       bbk = (0.d0,0.d0)
       !
       ALLOCATE(ps(nkb,nkb))
       !
       DO ik1 = 1, nks
          !
          ! Calculate beta-functions vkb for a given k point.
          !
          CALL init_us_2(ngk(ik1),igk_k(:,ik1),xk(1,ik1),vkb)
          !
          ! Calculate the coefficients B_ij defined by Eq.(15).
          ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
          !
          CALL zgemm('C','N',nkb,nkb,ngk(ik1),(1.d0,0.d0),vkb, &
                 & lda,vkb,lda,(0.d0,0.d0),bbk(1,1,ik1),nkb)
          !
#if defined(__MPI)
          CALL mp_sum(bbk(:,:,ik1), intra_bgrp_comm)
#endif
          !
          ps(:,:) = (0.d0,0.d0)
          !
          ! Calculate the product of q_nm and B_ij of Eq.(16).
          !
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               ps(ikb,ii) = ps(ikb,ii) + &
                                & bbk(jkb,ii, ik1) * qq(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO
          !
          ! Add an identity to q_nm * B_ij [see Eq.(16)].
          ! ps = (1 + q*B)
          !
          DO ii = 1, nkb
             ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
          ENDDO
          !
          ! Invert matrix: (1 + q*B)^{-1}
          !
          CALL invmat( nkb, ps )
          ! 
          ! Use the array bbk as a work space in order to save memory.
          !
          bbk(:,:,ik1) = (0.d0,0.d0)
          !
          ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
          !
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               bbk(ii,jkb,ik1) = bbk(ii,jkb,ik1) &
                                           & - ps(ii,ikb) * qq(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO
          !
       ENDDO ! loop on k points
       DEALLOCATE(ps)
       !
    ENDIF
    !
    IF (n.NE.ngk(ik)) CALL errore( 'sm1_psi_k', &
                    & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k point 'ik'.
    !
    CALL init_us_2(n,igk_k(:,ik),xk(1,ik),vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n,vkb,psi,becp,m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb,m))
    ps(:,:) = (0.d0,0.d0)
    !
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi> 
    !
    DO ibnd=1,m
       DO jkb=1,nkb
          DO ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd) + bbk(jkb,ii,ik) * becp%k(ii,ibnd)
          ENDDO
       ENDDO
    ENDDO
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta>
    !
    CALL ZGEMM( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
               & lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_k

END SUBROUTINE lr_sm1_psi


!----------------------------------------------------------------------------
SUBROUTINE lr_sm1_psiq (recalculate, ik, lda, n, m, psi, spsi)
  !----------------------------------------------------------------------------
  !
  ! This subroutine applies the S^{-1} matrix to m wavefunctions psi
  ! and puts the results in spsi.
  ! See Eq.(13) in B. Walker and R. Gebauer, JCP 127, 164106 (2007).
  ! Requires the products of psi with all beta functions
  ! in array becp(nkb,m) (calculated in h_psi or by ccalbec)
  !
  ! INPUT: recalculate   Decides if the overlap of beta 
  !                      functions is recalculated or not.
  !                      This is needed e.g. if ions are moved 
  !                      and the overlap changes accordingly.
  !        ik            k point under consideration
  !        lda           leading dimension of arrays psi, spsi
  !        n             true dimension of psi, spsi
  !        m             number of states psi
  !        psi           the wavefunction to which the S^{-1} 
  !                      matrix is applied
  ! OUTPUT: spsi = S^{-1}*psi  
  !
  ! Written by Iurii Timrov (2016) on the basis of lr_sm1_psi.
  ! Note: The difference of this routine from lr_sm1_psi is that 
  ! it can be used when there is a perturbation with the finite
  ! (transferred) momentum q (e.g. by PHonon and turboEELS codes). 
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE klist,            ONLY : xk, igk_k, ngk
  USE qpoint,           ONLY : ikks, ikqs, nksq
  USE uspp,             ONLY : okvan, vkb, nkb, qq
  USE uspp_param,       ONLY : nh, upf
  USE ions_base,        ONLY : ityp,nat,ntyp=>nsp
  USE becmod,           ONLY : bec_type, becp, calbec
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE matrix_inversion
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in)      :: recalculate
  INTEGER, INTENT(in)      :: lda, n, m, ik
  COMPLEX(DP), INTENT(in)  :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(out) :: spsi(lda*npol,m)
  !
  LOGICAL :: recalc
  !
  CALL start_clock( 'lr_sm1_psiq' )
  !
  recalc = recalculate
  !
  IF ( gamma_only ) THEN
     CALL errore( 'lr_sm1_psiq', 'gamma_only is not supported', 1 )
  ELSEIF (noncolin) THEN
     CALL sm1_psiq_nc()
  ELSE
     CALL sm1_psiq_k()
  ENDIF
  !
  CALL stop_clock( 'lr_sm1_psiq' )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE sm1_psiq_k()
    !-----------------------------------------------------------------------
    !
    ! k-points version
    ! Note: the array bbk must be deallocated somewhere
    ! outside of this routine.
    !
    USE lrus,      ONLY : bbk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    INTEGER :: ik1, & ! dummy index for k points
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B 
    ! (they are already allocated) but use the ones which were already
    ! calculated and saved (if recalc=.false.).
    !
    IF (.NOT.ALLOCATED(bbk)) recalc = .true.
    IF (recalc .AND. ALLOCATED(bbk)) DEALLOCATE(bbk)
    !
    ! If recalc=.true. we (re)calculate the coefficients lambda_nm defined by Eq.(16).
    !
    IF ( recalc ) THEN
       !
       ALLOCATE(bbk(nkb,nkb,nksq))
       bbk = (0.0d0,0.0d0)
       ! 
       ALLOCATE(ps(nkb,nkb))
       !
       DO ik1 = 1, nksq
          !
          ikk  = ikks(ik1)
          ikq  = ikqs(ik1)
          npwq = ngk(ikq)
          !
          ! Calculate beta-functions vkb for a given k+q point.
          !
          CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
          !
          ! Calculate the coefficients B_ij defined by Eq.(15).
          ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
          !
          CALL zgemm('C','N',nkb,nkb,npwq,(1.d0,0.d0),vkb, &
                & lda,vkb,lda,(0.d0,0.d0),bbk(1,1,ik1),nkb)
          !
#if defined(__MPI)
          CALL mp_sum(bbk(:,:,ik1), intra_bgrp_comm)
#endif
          ! 
          ps(:,:) = (0.d0,0.d0)
          !
          ! Calculate the product of q_nm and B_ij of Eq.(16).
          !
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            !
                            jkb=ijkb0 + jh
                            !
                            do ih=1,nh(nt)
                               !
                               ikb = ijkb0 + ih
                               !
                               ps(ikb,ii) = ps(ikb,ii) + bbk(jkb,ii,ik1)*qq(ih,jh,nt)
                               !
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
          !
          ! Add an identity to q_nm * B_ij [see Eq.(16)].
          ! ps = (1 + q*B)
          ! 
          DO ii=1,nkb
             ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
          ENDDO
          !
          ! Invert matrix: (1 + q*B)^{-1}
          ! 
          CALL invmat( nkb, ps )
          !
          ! Use the array bbk as a work space in order to save the memory.
          !
          bbk(:,:,ik1) = (0.d0,0.d0)
          !
          ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
          ! Let us use the array bbk and put there the values of the lambda
          ! coefficients.
          !
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            !
                            jkb=ijkb0 + jh
                            !
                            do ih=1,nh(nt)
                               !
                               ikb = ijkb0 + ih
                               !
                               bbk(ii,jkb,ik1) = bbk(ii,jkb,ik1) - ps(ii,ikb) * qq(ih,jh,nt)
                               !
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
          !
       ENDDO ! loop on k points
       !
       DEALLOCATE(ps)
       !
    ENDIF
    !
    ! Now set up the indices ikk and ikq such that they
    ! correspond to the points k and k+q using the index ik,
    ! which was passed to this routine as an input.
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    !
    IF (n.NE.ngk(ikq)) CALL errore( 'sm1_psiq_k', &
                     & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k+q point.
    !
    CALL init_us_2 (n, igk_k(1,ikq), xk(1,ikq), vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k+q, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n, vkb, psi, becp, m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb,m))
    ps(:,:) = (0.d0,0.d0)
    !
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi>  
    !
    DO ibnd=1,m
       DO jkb=1,nkb
          DO ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd) + bbk(jkb,ii,ik) * becp%k(ii,ibnd)
          ENDDO
       ENDDO
    ENDDO
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta> 
    !
    CALL zgemm( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
                & lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
END SUBROUTINE sm1_psiq_k

SUBROUTINE sm1_psiq_nc()
    !-----------------------------------------------------------------------
    !
    ! Noncollinear case
    ! Note: 1) the implementation of this routine is not finished...
    !       2) the array bbnc must be deallocated somewhere
    !          outside of this routine.
    !
    USE uspp,       ONLY : qq_so
    USE spin_orb,   ONLY : lspinorb
    USE lrus,       ONLY : bbnc
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ipol
    INTEGER :: ik1, & ! dummy index for k points
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    COMPLEX(DP), ALLOCATABLE :: ps(:,:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*npol*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    CALL errore( 'sm1_psiq_nc', 'USPP + noncolin is not implemented', 1 )
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B 
    ! (they are already allocated) but use the ones which were already calculated 
    ! and saved (if recalc=.false.).
    !
    IF (.NOT.ALLOCATED(bbnc)) recalc = .true.
    IF (recalc .AND. ALLOCATED(bbnc)) DEALLOCATE(bbnc)
    !
    ! If recalc=.true. we (re)calculate the coefficients lambda_nm defined by Eq.(16).
    !
    IF ( recalc ) THEN
       !
       ALLOCATE(bbnc(nkb,nkb,nspin_mag,nksq))
       bbnc = (0.0d0,0.0d0)
       ! 
       ALLOCATE(ps(nkb,nkb,nspin_mag))
       !
       DO ik1 = 1, nksq
          !
          ikk  = ikks(ik1)
          ikq  = ikqs(ik1)
          npwq = ngk(ikq)
          !
          ! Calculate beta-functions vkb for a given k+q point.
          !
          CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
          !
          ! Calculate the coefficients B_ij defined by Eq.(15).
          ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
          !
          CALL ZGEMM('C','N',nkb,nkb,npwq,(1.d0,0.d0),vkb, &
               & lda,vkb,lda,(0.d0,0.d0),bbnc(1,1,1,ik1),nkb)
          !
          IF (lspinorb) THEN
             bbnc(:,:,2,ik1) = bbnc(:,:,1,ik1)
             bbnc(:,:,3,ik1) = bbnc(:,:,1,ik1)
             bbnc(:,:,4,ik1) = bbnc(:,:,1,ik1)
          ENDIF
          !
#if defined(__MPI)
          CALL mp_sum(bbnc(:,:,:,ik1), intra_bgrp_comm)
#endif
          !
          ps(:,:,:) = (0.d0,0.d0)
          !
          ! Calculate the product of q_nm and B_ij of Eq.(16).
          !
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            !
                            jkb=ijkb0 + jh
                            !
                            do ih=1,nh(nt)
                               !
                               ikb = ijkb0 + ih
                               !
                               if (lspinorb) then
                                  do ipol=1,4
                                     ps(ikb,ii,ipol) = ps(ikb,ii,ipol) + &
                                        & bbnc(jkb,ii,ipol,ik1) * qq_so(ih,jh,ipol,nt)
                                  enddo
                               else
                                  ps(ikb,ii,1) = ps(ikb,ii,1) + &
                                        & bbnc(jkb,ii,1,ik1) * qq(ih,jh,nt) 
                               endif
                               !
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
          !
          ! Add an identity to q_nm * B_ij [see Eq.(16)].
          ! ps = (1 + q*B)
          ! 
          DO ii=1,nkb
             IF (lspinorb) THEN
                DO ipol=1,4
                   ps(ii,ii,ipol) = ps(ii,ii,ipol) + (1.d0,0.d0)
                ENDDO
             ELSE
                ps(ii,ii,1) = ps(ii,ii,1) + (1.d0,0.d0)
             ENDIF
          ENDDO
          !
          ! Invert matrix: (1 + q*B)^{-1}
          ! WARNING: How to do the invertion of the matrix in the spin-orbit case,
          ! when there are 4 components?
          !
          IF (lspinorb) THEN
             DO ipol=1,4
                CALL invmat( nkb, ps(:,:,ipol) )  
             ENDDO
          ELSE
             CALL invmat ( nkb, ps(:,:,1) )
          ENDIF
          !
          ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
          !
          ! Use the array bbnc as a workspace and put there 
          ! the values of the lambda coefficients.
          !
          bbnc(:,:,:,ik1) = (0.d0,0.d0)
          !
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            do ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               if (lspinorb) then
                                  bbnc(ii,jkb,1,ik1) = bbnc(ii,jkb,1,ik1) &
                                     & - ps(ii,ikb,1) * qq_so(ih,jh,1,nt) &
                                     & - ps(ii,ikb,2) * qq_so(ih,jh,3,nt)
                                  !
                                  bbnc(ii,jkb,2,ik1) = bbnc(ii,jkb,2,ik1) &
                                     & - ps(ii,ikb,1) * qq_so(ih,jh,2,nt) &
                                     & - ps(ii,ikb,2) * qq_so(ih,jh,4,nt)
                                  !
                                  bbnc(ii,jkb,3,ik1) = bbnc(ii,jkb,3,ik1) &
                                     & - ps(ii,ikb,3) * qq_so(ih,jh,1,nt) &
                                     & - ps(ii,ikb,4) * qq_so(ih,jh,3,nt)
                                  !
                                  bbnc(ii,jkb,4,ik1) = bbnc(ii,jkb,4,ik1) &
                                     & - ps(ii,ikb,3) * qq_so(ih,jh,2,nt) &
                                     & - ps(ii,ikb,4) * qq_so(ih,jh,4,nt)
                                else
                                  bbnc(ii,jkb,1,ik1) = bbnc(ii,jkb,1,ik1) &
                                            & - ps(ii,ikb,1)*qq(ih,jh,nt)
                                endif
                                !
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
          !
       enddo ! loop on k points
       !
       deallocate(ps)
       !
    endif
    !
    ! Now set up the indices ikk and ikq such that they
    ! correspond to the points k and k+q using the index ik,
    ! which was passed to this routine as an input.
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    !
    IF (n.NE.ngk(ikq)) CALL errore( 'sm1_psiq_nc', &
                     & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k+q point.
    !
    CALL init_us_2 (n, igk_k(1,ikq), xk(1,ikq), vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k+q, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n, vkb, psi, becp, m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb,npol,m))
    ps(:,:,:) = (0.d0,0.d0)
    !    
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi>  
    !
    DO ibnd=1,m
       DO jkb=1,nkb
          DO ii=1,nkb
             IF (lspinorb) THEN
                !
                ps(jkb,1,ibnd) = ps(jkb,1,ibnd) &
                    & + bbnc(jkb,ii,1,ik) * becp%nc(ii,1,ibnd) &
                    & + bbnc(jkb,ii,2,ik) * becp%nc(ii,2,ibnd)
                !
                ps(jkb,2,ibnd) = ps(jkb,2,ibnd) &
                    & + bbnc(jkb,ii,3,ik) * becp%nc(ii,1,ibnd) &
                    & + bbnc(jkb,ii,4,ik) * becp%nc(ii,2,ibnd)
                !
             ELSE
                DO ipol=1,npol
                   ps(jkb,ipol,ibnd) = ps(jkb,ipol,ibnd) &
                   & + bbnc(jkb,ii,1,ik) * becp%nc(ii,ipol,ibnd) 
                ENDDO
             ENDIF
           ENDDO
       ENDDO
    ENDDO
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta> 
    !
    CALL ZGEMM( 'N', 'N', n, m*npol, nkb, (1.D0, 0.D0), vkb, &
                     & lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
END SUBROUTINE sm1_psiq_nc

END SUBROUTINE lr_sm1_psiq
