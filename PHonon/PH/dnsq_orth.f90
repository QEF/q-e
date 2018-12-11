!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE dnsq_orth() 
  !-----------------------------------------------------------------------
  !
  ! DFPT+U : This routine calculates, in case of USPP, the bare variation 
  ! of the occupation matrix due to orthogonality contraints.
  !
  ! dnsorth_cart(m1,m2,ispin,I,icart,na) = 
  !   - \sum_{k,n,n'} wgg(n,n',k) * <psi(n,k,ispin)| S_{k}\phi_(k,I,m2)>  
  !                             * < S_{k+q}\phi_(k+q,I,m1)| psi(n',k+q,ispin)> 
  !    * \sum_{l1,l2} [ <psi(n',k+q,ispin)| d_{L,icart}beta(k+q,L,l1)> * q(L,l1,l2) *
  !                     <beta(k,L,l2)| psi(n,k,ispin)> 
  !                   + <psi(n',k+q,ispin)| beta(k+q,L,l1)> * q(L,l1,l2) *
  !                     <d_{L,icart}beta(k,L,l2)| psi(n,k,ispin)> ]                                !
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, offsetU, nwfcU
  USE ldaU_ph,       ONLY : swfcatomk, swfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            proj1, proj2, dnsorth_cart, &
                            read_dns_bare, dnsorth
  USE klist,         ONLY : xk, wk,  ngk, igk_k
  USE wvfct,         ONLY : npwx, wg, nbnd 
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma, ofsbeta
  USE units_lr,      ONLY : iuatswfc
  USE uspp_param,    ONLY : nh
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : evq
  USE uspp,          ONLY : okvan, nkb, vkb
  USE control_flags, ONLY : iverbosity
  USE mp_global,     ONLY : intra_pool_comm, inter_pool_comm
  USE mp,            ONLY : mp_sum, mp_bcast 
  USE io_files,      ONLY : seqopn 
  USE buffers,       ONLY : get_buffer
  USE mp_world,      ONLY : world_comm
  USE mp_images,     ONLY : intra_image_comm
  USE doubleprojqq_module
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: i, j, k, ios, icart, nt, na, l, ina, ih, n,               &
             ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, jbnd, is, &
             iat, ic, nti, ibeta, imode, na_icart, ldim, iundnsorth,   &
             npw, npwq, ik, ikk, ikq
  COMPLEX(DP), ALLOCATABLE :: dpqq(:), dpqq1(:), sum_dpqq(:,:)
  REAL(DP), ALLOCATABLE :: wgg(:,:,:)
  LOGICAL :: exst 
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock( 'dnsq_orth' )
  !
  ios = 0 
  !
  ALLOCATE (dpqq(nbnd))
  ALLOCATE (dpqq1(nbnd))
  ALLOCATE (sum_dpqq(nbnd,nbnd))
  ALLOCATE (wgg(nbnd,nbnd,nksq))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  !  
  ldim = 2 * Hubbard_lmax + 1
  !
  ! USPP: compute the weights as in square bracket 
  ! of Eq. (27) of A. Dal Corso PRB 64, 235118 (2001).
  ! Or see Eq. (B19) in the same reference.
  !
  CALL compute_weight (wgg)
  ! 
  ALLOCATE (dnsorth_cart(ldim,ldim,nspin,nat,3,nat))
  dnsorth_cart = (0.d0, 0.d0)  
  !
  exst = .FALSE.
  !
  !The unit number 
  iundnsorth = 38
  !
  IF (ionode) CALL seqopn (iundnsorth, 'dnsorth_cart', 'formatted', exst)
  !
  ! Read (unsymmetrized) dnsorth_cart from file (if it was already computed)
  ! 
  IF (read_dns_bare) THEN 
     ! 
     WRITE(stdout,*) 'Does the file dnsorth_cart exist in tmp_dir and is it read correctly?'
     ! 
     IF (ionode) THEN 
        IF (exst) THEN
           READ(iundnsorth,*,iostat=ios) dnsorth_cart
           IF (ios==0) THEN 
              WRITE(stdout,*) '...yes. the dnsorth_cart matrix was read correctly from file'
           ELSE
              WRITE(stdout,*) '...no. the file exists but it seems to be corrupted'
           ENDIF
        ELSE
           WRITE(stdout,*) '...no. the file does not exist. '
        ENDIF
     ENDIF
     ! 
     CALL mp_bcast(exst, ionode_id, world_comm)           
     CALL mp_bcast(ios,  ionode_id, world_comm)    
     IF (exst.and.ios==0) CALL mp_bcast(dnsorth_cart, ionode_id, world_comm)    
     !
     ! IT: Is it needed to broadcast for intra_image_comm?
     CALL mp_bcast(ios,  ionode_id, intra_image_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsorth_cart, ionode_id, intra_image_comm)
     ! 
  ENDIF
  ! 
  ! Compute dnsorth_cart (if it was not already done)
  ! 
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !
     WRITE(stdout,'(/5x,"Calculating the dnsorth_cart matrix...")')
     dnsorth_cart = (0.d0, 0.d0)  
     ! 
     DO ik = 1, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        IF (lsda) current_spin = isk (ikk)
        !
        ! Reads unperturbed KS wavefuctions psi(k) and psi(k+q)
        !
        IF (nksq.GT.1) THEN
           CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           IF (.NOT.lgamma) CALL get_buffer (evq, lrwfc, iuwfc, ikq)
        ENDIF
        !
        ! Read the atomic orbitals S*\phi at k and k+q from file (unit iuatswfc)
        !  
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)        
        IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
        ! 
        ! Compute beta functions at k (vkb) and at k+q (vkbkpq)
        ! 
        CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(:,ikq), vkbkpq)
        !
        ! Calculate:
        ! proj1 (ibnd, ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
        ! proj2 (ibnd, ihubst) = < S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k+q) >
        !
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO m = 1, 2*Hubbard_l(nt)+1
                 ihubst = offsetU(nah) + m   ! I m index
                 DO ibnd = 1, nbnd
                    proj1(ibnd,ihubst) = ZDOTC (npw,  swfcatomk(:,ihubst),   1, evc(:,ibnd), 1)
                    proj2(ibnd,ihubst) = ZDOTC (npwq, swfcatomkpq(:,ihubst), 1, evq(:,ibnd), 1)
                 ENDDO
              ENDDO
           ENDIF 
        ENDDO  
        !
        CALL mp_sum (proj1, intra_pool_comm)  
        CALL mp_sum (proj2, intra_pool_comm)
        ! 
        DO na = 1, nat
           !
           nt = ityp(na)
           !
           DO icart = 1, 3 
              !   
              ! Calculates the derivatives \delta(na,icart) of beta functions
              ! only for j=na, at k and k+q (for all the states l,l')
              !
              DO ih = 1, nh(nt)
                 !
                 ibeta = ofsbeta(na) + ih
                 !
                 CALL dwfc (npw, igk_k(:,ikk), ikk, icart, &
                            vkb(:,ibeta), dvkb(:,ibeta,icart))
                 IF (.NOT.lgamma) &
                 CALL dwfc (npwq, igk_k(:,ikq), ikq, icart, &
                            vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
                 !     
              ENDDO
              !
              sum_dpqq = (0.d0, 0.d0) 
              !
              ! Calculate for all bands jbnd (in evq) and for the band ibnd:
              ! dpqq  = \sum_{l1,l2} <psi(n',k+q,ispin) | d_{L,icart}beta(k+q,L,l1)> 
              !                      * qq_nt(L,l1,l2) * <beta(k,L,l2)| psi(n,k)>
              ! dpqq1 = \sum_{l1,l2} <psi(n',k+q,ispin) | beta(k+q,L,l1)>
              !                      * qq_nt(L,l1,l2) * <d_{L,icart}beta(k,L,l2)| psi(n,k)>
              !
              DO ibnd = 1, nbnd
                 !      
                 CALL doubleprojqq (na, evq, dvkbkpq(:,:,icart), vkb, evc(:,ibnd),  &
                                    npwq, npw, dpqq)
                 CALL doubleprojqq (na, evq, vkbkpq, dvkb(:,:,icart), evc(:,ibnd), &
                                    npwq, npw, dpqq1)
                 sum_dpqq(ibnd,:) = dpqq + dpqq1
                 !
              ENDDO
              ! 
              ! Finally calculate dnsorth_cart
              ! 
              DO nah = 1, nat
                 !
                 nt = ityp(nah)
                 !
                 IF (is_hubbard(nt)) THEN
                    !
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       ! 
                       ihubst1 = offsetU(nah) + m1
                       !
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          !
                          ihubst2 = offsetU(nah) + m2
                          ! 
                          DO ibnd = 1, nbnd
                             DO jbnd = 1, nbnd
                                dnsorth_cart (m1,m2,current_spin,nah,icart,na) = &
                                dnsorth_cart (m1,m2,current_spin,nah,icart,na) - &
                                    wgg(ibnd,jbnd,ik) * wk(ikk) * &              
                                    CONJG(proj1(ibnd,ihubst2)) * proj2(jbnd,ihubst1) * &
                                    sum_dpqq(ibnd,jbnd)              
                             ENDDO
                          ENDDO
                          !
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDIF
                 !
              ENDDO
              !
           ENDDO ! icart 
           !
        ENDDO ! na
        !
     ENDDO ! ik
     !
     CALL mp_sum (dnsorth_cart, inter_pool_comm)
     !
     ! If nspin=1 k-point weight is normalized to 2 el/band 
     ! in the whole BZ but we are interested in dns of one spin component
     !
     IF (nspin.EQ.1) dnsorth_cart = 0.5d0 * dnsorth_cart
     ! 
     ! Write (unsymmetrized) dnsorth_cart to file
     !
     IF (ionode) WRITE(iundnsorth,*) dnsorth_cart
     ! 
  ENDIF
  !
  ! Close the unit iundnsorth
  ! 
  IF (ionode) CLOSE (unit=iundnsorth,status='keep')  
  ! 
  ! Symmetrize dnsorth_cart, and also save the result
  ! in the pattern basis dnsorth
  !
  CALL sym_dns_wrapper (ldim, dnsorth_cart, dnsorth)
  !
  ! Write symmetryzed dnsorth_cart in cartesian coordinates
  ! to the standard output
  !
  IF (iverbosity==1) THEN 
     WRITE(stdout,*) 'DNS_ORTH SYMMETRIZED IN CARTESIAN COORDINATES'
     DO na = 1, nat
        DO icart = 1, 3 
           WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
           DO nah = 1, nat
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 DO is = 1, nspin
                    WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       WRITE( stdout,'(14(f15.10,1x))') dnsorth_cart (m1,:,is,nah,icart,na) 
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     WRITE(stdout,*)
  ENDIF

  DEALLOCATE (dpqq)
  DEALLOCATE (dpqq1)
  DEALLOCATE (sum_dpqq)
  DEALLOCATE (wgg)  
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  !
  CALL stop_clock( 'dnsq_orth' )
  ! 
  RETURN
  ! 
END SUBROUTINE dnsq_orth
!---------------------------------------------------------------------
