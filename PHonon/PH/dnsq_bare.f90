!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dnsq_bare 
  !-----------------------------------------------------------------------
  !
  ! DFPT+U : This routine calculates dnsbare, i.e. the bare variation 
  ! of the occupation matrix ns under the perturbation q, due to a 
  ! displacement of the atom L in the direction icart. 
  !
  ! dnsbare(m1,m2,ispin,I,icart,L)= 
  !  = \sum_{k,n} [ <psi(n,k,ispin)| S_{k}\phi_(k,I,m1)> * 
  !                 <\Delta^{L icart}_{-q}(S_{k+q}\phi_(k+q,I,m2))|psi(n,k,ispin)>  
  !               + <psi(n,k,ispin)| S_{k}\phi_(k,I,m2)> * 
  !                 <\Delta^{L icart}_{-q}(S_{k+q}\phi_(k+q,I,m1))|psi(n,k,ispin)> ]
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  ! 
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, swfcatomk, swfcatomkpq, dwfcatomkpq, &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            dnsbare, dnsbare_all_modes, proj1, proj2, read_dns_bare
  USE wvfct,         ONLY : npwx, wg, nbnd 
  USE uspp,          ONLY : vkb, nkb 
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma, ofsbeta
  USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh, nhm 
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_bands,      ONLY : intra_bgrp_comm 
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : seqopn  
  USE buffers,       ONLY : get_buffer
  USE control_flags, ONLY : iverbosity
  USE wavefunctions, ONLY : evc  
  USE mp_images,     ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: i, j, k, icart, ios, is, na, nt, l, ih,        &
             ibeta, n, ihubst, ihubst1, ihubst2, nah, m, m1, m2,  &
             ibnd, ldim, npw, npwq, ik, ikk, ikq, iundnsbare
  COMPLEX(DP), ALLOCATABLE :: dqsphi(:,:), dmqsphi(:,:), dwfcatom_(:) 
  COMPLEX(DP), EXTERNAL :: ZDOTC
  LOGICAL :: exst
  ! 
  CALL start_clock( 'dnsq_bare' )
  ! 
  ios = 0
  !
  ALLOCATE (dqsphi(npwx,nwfcU)) 
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dwfcatom_(npwx))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ! 
  proj1   = (0.d0, 0.d0) 
  proj2   = (0.d0, 0.d0) 
  dnsbare = (0.d0, 0.d0) 
  !
  ! Note: sdwfcatomkpq is not needed in this routine when we call delta_sphi,
  ! so we set it to zero and do not calculate.
  !
  sdwfcatomk   = (0.d0, 0.d0)
  sdwfcatomkpq = (0.d0, 0.d0)
  !
  ldim = 2 * Hubbard_lmax + 1
  !
  exst = .FALSE.
  !
  ! The unit number
  iundnsbare = 37
  !
  IF (ionode) CALL seqopn (iundnsbare, 'dnsbare', 'formatted', exst)
  !
  ! Read (unsymmetrized) dnsbare from file (if it was already computed)
  !
  IF (read_dns_bare) THEN 
     !   
     WRITE(stdout,*) 'Does the file dnsbare exist in tmp_dir and is it read correctly?'
     ! 
     IF (ionode) THEN 
        IF (exst) THEN
           READ(iundnsbare,*,iostat=ios) dnsbare
           IF (ios==0) THEN 
              WRITE(stdout,*) '...yes. the dnsbare matrix was read correctly from file'
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
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, world_comm) 
     !
     ! IT: Is it needed to broadcast for intra_image_comm?
     CALL mp_bcast(ios,  ionode_id, intra_image_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, intra_image_comm)
     !
  ENDIF
  !
  ! Compute dnsbare (if it was not already done)
  !
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !
     WRITE(stdout,'(/5x,"Calculating the dnsbare matrix...")')
     dnsbare = (0.d0, 0.d0)
     ! 
     DO ik = 1, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        IF (lsda) current_spin = isk(ikk)
        !
        ! Read unperturbed KS wavefuctions psi(k)
        !
        IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! Compute beta functions at k (vkb) and at k+q (vkbkpq)
        ! 
        CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(:,ikq), vkbkpq)
        ! 
        ! Compute the derivatives of beta functions d^{icart}beta at k and k+q
        ! for all bands and for the 3 cartesian directions
        !
        DO icart = 1, 3
           DO na = 1, nat 
              nt = ityp(na)
              DO ih = 1, nh(nt)
                 !
                 ibeta = ofsbeta(na) + ih
                 !
                 CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                            vkb(:,ibeta), dvkb(:,ibeta,icart))
                 !
                 IF (.NOT.lgamma) &
                 CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, & 
                            vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
                 !
              ENDDO
           ENDDO
        ENDDO
        ! 
        ! Read the atomic orbitals \phi at k and k+q from file (unit iuatwfc)
        ! 
        CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
        IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
        ! 
        ! Read S*\phi at k from file (unit iuatswfc)
        !
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)      
        !
        DO na = 1, nat
           ! 
           DO icart = 1, 3 
              ! 
              dqsphi  = (0.d0, 0.d0)  
              dmqsphi = (0.d0, 0.d0)
              !
              DO nah = 1, nat
                 !
                 nt = ityp(nah)
                 !
                 IF (is_hubbard(nt)) THEN
                    !    
                    DO m = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst = offsetU(nah) + m ! I m index
                       ! 
                       IF (nah==na) THEN 
                          !  
                          ! Calculate |d_icart\phi_(k,I,m)) >
                          !
                          CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                                     wfcatomk(:,ihubst), dwfcatom_) 
                          !
                          ! Apply the S matrix: | S_{k} d_^(I,icart)\phi_(k,I,m) > 
                          !
                          CALL swfc (npw, 1, vkb, dwfcatom_, sdwfcatomk(:,ihubst))
                          !
                       ENDIF
                       !
                       ! Calculate dmqsphi = |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                       !
                       CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq, &
                                        sdwfcatomk, sdwfcatomkpq, vkb, vkbkpq, dvkb(:,:,icart), &
                                        dvkbkpq(:,:,icart), dqsphi, dmqsphi, 0)  
                       !
                       ! Calculate: 
                       ! proj1 (ibnd,ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
                       ! proj2 (ibnd,ihubst) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k)>
                       !
                       DO ibnd = 1, nbnd
                          proj1 (ibnd,ihubst) = ZDOTC (npw, swfcatomk(1,ihubst), 1, evc(1,ibnd), 1)
                          proj2 (ibnd,ihubst) = ZDOTC (npw, dmqsphi(1,ihubst ), 1, evc(1,ibnd), 1)
                       ENDDO
                       !  
                    ENDDO 
                    !
                 ENDIF
                 !
              ENDDO
              !
              CALL mp_sum (proj1, intra_bgrp_comm)
              CALL mp_sum (proj2, intra_bgrp_comm)
              !
              ! Finally calculate dnsbare.
              ! dnsbare is symmetric: dnsbare(m1,m2) = dnsbare(m2,m1)
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
                       DO m2 = m1, 2*Hubbard_l(nt)+1
                          !
                          ihubst2 = offsetU(nah) + m2
                          !
                          DO ibnd = 1, nbnd
                             dnsbare (m1,m2,current_spin,nah,icart,na) = &
                             dnsbare (m1,m2,current_spin,nah,icart,na) + wg(ibnd,ikk) * &
                                  ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2)  + &
                                    CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
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
     CALL mp_sum (dnsbare, inter_pool_comm)
     !
     ! Filling the m1 m2 lower triangular part 
     !
     DO m1 = 2, ldim
        DO m2 = 1, m1-1
           dnsbare(m1,m2,:,:,:,:) = dnsbare(m2,m1,:,:,:,:)
        ENDDO
     ENDDO
     !
     ! In nspin.eq.1 k-point weight is normalized to 2 el/band 
     ! in the whole BZ but we are interested in dns of one spin component
     !
     IF (nspin.EQ.1) dnsbare = 0.5d0 * dnsbare
     ! 
     ! Write (unsymmetrized) dnsbare to file
     !  
     IF (ionode) WRITE(iundnsbare,*) dnsbare
     !
  ENDIF
  !
  ! Close the unit iundnsbare
  !
  IF (ionode) CLOSE (unit=iundnsbare,status='keep')   
  !
  ! Symmetrize dnsbare, and also save the result
  ! in the pattern basis dnsbare_all_modes
  !
  CALL sym_dns_wrapper (ldim, dnsbare, dnsbare_all_modes) 
  !
  ! Write symmetryzed dnsbare in cartesian coordinates 
  ! to the standard output
  ! 
  IF (iverbosity==1) THEN 
     WRITE(stdout,*) 'DNS_BARE SYMMETRIZED IN CARTESIAN COORDINATES'
     DO na = 1, nat        
        DO icart = 1, 3
           WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
           DO nah = 1, nat  
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 DO is = 1, nspin
                    WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                    DO m1 = 1, ldim
                       WRITE(stdout,'(10(f15.10,1x))') dnsbare (m1,:,is,nah,icart,na)
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     WRITE(stdout,*)
  ENDIF
  ! 
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dwfcatom_)
  DEALLOCATE (proj1 )
  DEALLOCATE (proj2 )
  ! 
  CALL stop_clock ('dnsq_bare')
  ! 
  RETURN
  !
END SUBROUTINE dnsq_bare
!----------------------------------------------------------------------------------------
