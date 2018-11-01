!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_bare
  !---------------------------------------------------------------------
  !
  ! DFPT+U: This routine does two tasks:
  ! 1) Computes d2ns_bare (very slow) or reads it
  ! 2) Adds to the dynamical matrix the bare part due the Hubbard term:  
  !
  !  \sum_{nah} Hubbard_U(nah) \sum_{is, m1, m2} [ 
  !      (0.5*delta_m1m2 - ns(m1,m2,is,nah)) * d2ns_bare 
  !      - conjg(dnsbare ( m1, m2, is, nah, icart,na)) 
  !            * dnsbare ( m1, m2, is, nah, jcart,nap) ]
  !
  !  Addition of the J0 term:
  !  
  !  \sum_{nah} Hubbard_J0(nah) \sum_{is, m1, m2} [ 
  !      ns(m1,m2,-is,nah) * d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart)  
  !      + conjg(dnsbare ( m1, m2, is, nah, icart,na))*
  !              dnsbare ( m1, m2, -is, nah, jcart,nap) ]
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, Hubbard_J0, &
                            Hubbard_lmax, offsetU, nwfcU
  USE ldaU_ph,       ONLY : dnsbare, wfcatomk, wfcatomkpq, dvkb, vkbkpq, &
                            dvkbkpq, proj1, projpb, projpdb, swfcatomk, &
                            effU, read_dns_bare, d2ns_type
  USE wavefunctions, ONLY : evc
  USE units_lr,      ONLY : iuwfc, lrwfc, iuatwfc, iuatswfc
  USE uspp,          ONLY : vkb, nkb
  USE uspp_param,    ONLY : nh
  USE klist,         ONLY : xk, ngk, igk_k
  USE control_lr,    ONLY : lgamma, ofsbeta
  USE io_files,      ONLY : nwordwfcU, tmp_dir 
  USE lsda_mod,      ONLY : lsda, current_spin, isk, nspin
  USE modes,         ONLY : u, nmodes
  USE dynmat,        ONLY : dyn
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE wvfct,         ONLY : npwx, nbnd
  USE control_flags, ONLY : iverbosity
  USE d2nsq_bare_module
  USE scf,           ONLY : rho
  USE mp_global,     ONLY : intra_pool_comm, inter_pool_comm       
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : seqopn
  USE buffers,       ONLY : get_buffer
  ! test: variables useful to symmetrize dynmat_hub_bare, 
  ! before (and only for) printing it
  USE lr_symm_base,  ONLY : rtau, nsymq, irotmq, minus_q
  USE modes,         ONLY : u
  USE qpoint,        ONLY : xq
  USE symm_base,     ONLY : irt, s, invs
  USE cell_base,     ONLY : at, bg
  USE dynmat,        ONLY : dyn_hub_bare 
  USE io_global,     ONLY : ionode, ionode_id, stdout
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: icart, jcart, na, nap, nah, &
             ihubst1, ihubst2, ibeta, nt, n, l, is, ina, ih, & 
             na_icart, nap_jcart, na_icar, m1, m2, op_spin, isi, &
             imode, jmode, ldim, ik, ikk, ikq, icar, ibnd, &
             ios, iund2nsbare, npw, npwq
  REAL(DP)    :: nsaux
  COMPLEX(DP) :: dnsaux1, dnsaux2, d2ns_bare_aux, d2ns_bare_k, work
  LOGICAL :: exst
  COMPLEX(DP), ALLOCATABLE :: d2ns_bare(:,:,:,:,:,:), dynwrk(:,:)  
  COMPLEX(DP), EXTERNAL :: ZDOTC  
  !
  CALL start_clock ( 'dynmat_hub_bare' )
  !
  ldim = 2*Hubbard_lmax + 1
  !
  ALLOCATE (dyn_hub_bare(3*nat,3*nat))  
  ALLOCATE (d2ns_bare(ldim,ldim,nspin,nat,3*nat,3*nat) )
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (projpb(nbnd,nkb))
  ALLOCATE (projpdb(nbnd,nkb,3))
  ALLOCATE (dynwrk(3*nat,3*nat))
  !
  dyn_hub_bare = (0.d0, 0.d0)
  dynwrk       = (0.d0, 0.d0) 
  d2ns_bare    = (0.d0, 0.d0)
  !
  exst = .false.
  iund2nsbare = 39
  !
  IF (ionode) CALL seqopn(iund2nsbare, 'd2nsbare', 'formatted', exst)
  !
  ! Read dnsbare
  !
  IF (read_dns_bare) THEN 
     !
     WRITE( stdout,*) 'D2NS_BARE FILE EXISTS IN TMP_DIR AND IS READ CORRECTLY?'
     !   
     IF (ionode) THEN
        !
        IF (exst) THEN
           !
           ios=0
           READ(iund2nsbare,*,iostat=ios) d2ns_bare
           !
           IF (ios==0) THEN 
              WRITE( stdout,*) '...YES. THE D2NS_BARE MATRIX WAS CORRECTLY READ FROM FILE'
           ELSE
              WRITE( stdout,*) '...NO. THE FILE EXISTS BUT IT SEEMS TO BE CORRUPTED'
           ENDIF
           !
        ELSE
           WRITE( stdout,*) 'THE FILE DOES NOT EXIST'
        ENDIF
        !
     ENDIF
     !
     CALL mp_bcast(exst, ionode_id, world_comm)       
     CALL mp_bcast(ios, ionode_id, world_comm)   
     !
     IF (exst .AND. ios==0) CALL mp_bcast(d2ns_bare, ionode_id, world_comm)     
     !           
  ENDIF
  !
  ! If dnsbare was not read, then compute it
  !
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !     
     d2ns_bare = (0.d0, 0.d0)
     ! 
     WRITE(stdout,'(/5x,"Calculating the d2ns_bare matrix. It might take a while!")')
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
           IF ((current_spin==1).and.(nspin==2)) THEN
              op_spin = 2
           ELSE
              op_spin = 1
           ENDIF
        ELSE        
           op_spin = 1     
        ENDIF
        !
        ! Read unperturbed wavefuctions psi(k)
        !
        IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! Compute the beta function vkb at k, and vkbkpq at k+q
        !
        CALL init_us_2 (npw, igk_k(1,ikk), xk(1,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkbkpq)
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
        ! Calculate proj1(ibnd) = <S(k)*phi(k,I,m')|psi(ibnd,k)> 
        !
        proj1 = (0.d0, 0.d0)
        !
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO m1 = 1, 2*Hubbard_l(nt)+1
                 ihubst1 = offsetU(nah) + m1   ! I m index
                 DO ibnd = 1, nbnd
                    proj1(ibnd, ihubst1) = ZDOTC (npw, swfcatomk(:,ihubst1), 1, evc(:,ibnd), 1)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        !
        ! Calculate: projpb(nbnd,nkb) = <psi|vkb> for all nbnd and vkb states, 
        ! useful in d2ns_bare_k 
        !
        projpb = (0.d0, 0.d0)
        !
        DO na = 1, nat    
           nt = ityp(na)  
           DO ih = 1, nh(nt)
              ibeta = ofsbeta(na) + ih
              DO ibnd = 1, nbnd
                 projpb(ibnd, ibeta) = ZDOTC (npw, evc(:,ibnd), 1, vkb(:,ibeta), 1)
              ENDDO
           ENDDO
        ENDDO
        !
#if defined(__MPI)
        CALL mp_sum(proj1, intra_pool_comm) 
        CALL mp_sum(projpb, intra_pool_comm)
#endif
        !
        ! Calculate the m1 m2 upper triangular part (UPT) of the UPT matrices
        !
        projpdb = (0.d0, 0.d0)
        !
        DO na = 1, nat  ! the displaced atom
           !
           DO icart = 1, 3
              !
              na_icart = 3*(na-1) + icart
              !
              DO nap = 1, nat ! the displaced atom
                 !
                 DO jcart = 1, 3
                    !
                    nap_jcart = 3*(nap-1) + jcart
                    !
                    IF (nap_jcart.GE.na_icart) THEN
                       !
                       ! Calculate the derivatives of the beta functions  
                       ! at k and k+q only at icart, jcart
                       ! and only afor atoms na, nap (for all the states)
                       !                  
                       DO icar = 1, 3
                          !
                          IF ((icar==icart) .OR. (icar==jcart)) THEN 
                             ! we want only icart jcart
                             !
                             DO nt = 1, ntyp
                                !
                                DO ina = 1, nat 
                                   !
                                   IF (ityp(ina).EQ.nt) THEN
                                      !
                                      IF ((ina==na) .OR. (ina==nap)) THEN 
                                         ! we want only dbeta at na or nap
                                         !
                                         DO ih = 1, nh(nt)
                                            !
                                            ibeta = ofsbeta(ina) + ih
                                            !
                                            CALL dwfc(npw, igk_k(:,ikk), ikk, icar, &
                                                      vkb(:,ibeta), dvkb(:,ibeta,icar))
                                            IF (.NOT.lgamma) &
                                            CALL dwfc(npwq, igk_k(:,ikq), ikq, icar, &
                                                      vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icar))
                                            ! 
                                            DO ibnd = 1, nbnd
                                               projpdb(ibnd, ibeta, icar) = &
                                                 ZDOTC (npw, evc(:,ibnd), 1, dvkb(:,ibeta,icar), 1)
                                            ENDDO
                                            !
                                         ENDDO ! ih
                                         !
                                      ENDIF
                                      !
                                   ENDIF
                                   !
                                ENDDO ! ina
                                !
                             ENDDO ! nt
                             !
                          ENDIF
                          !
                       ENDDO ! icar
                       !
#if defined(__MPI)
     CALL mp_sum(projpdb, intra_pool_comm) 
#endif
                       !
                       DO nah = 1, nat ! the Hubbard atom
                          !
                          nt = ityp(nah)
                          !
                          IF (is_hubbard(nt)) THEN
                             !                          
                             DO m1 = 1, 2*Hubbard_l(nt)+1
                                !
                                ihubst1 = offsetU(nah) + m1  ! I m index
                                !
                                DO m2 = m1, 2*Hubbard_l(nt)+1
                                   !
                                   ihubst2 = offsetU(nah) + m2 ! I m index
                                   !
                                   ! Calculate d2ns_bare_k  = the bare part of 
                                   ! the 2nd derivative of the occupation matrix ns 
                                   ! it is an object calculated at 
                                   ! (ik,na,icart,nap,jcart,ihubst1,ihubst2)
                                   !
                                   IF ((d2ns_type=='full') .OR. (d2ns_type=='fmmp')) THEN
                                      !           
                                      CALL d2nsq_bare_k (ik, icart, jcart, na, nap, nah, &
                                                         ihubst1, ihubst2, d2ns_bare_k)
                                      ! 
                                   ELSEIF ((d2ns_type=='diag') .OR. (d2ns_type=='dmmp')) THEN
                                      !
                                      CALL d2nsq_bare_k_diag (ik, icart, jcart, na, nap, nah, &
                                                              ihubst1, ihubst2, d2ns_bare_k)
                                      !
                                   ELSE
                                      CALL errore( 'dynmat_hub_bare', &
                                                    & 'wrong d2ns_type value', 1 )
                                   ENDIF
                                   !
                                   d2ns_bare (m1, m2, current_spin, nah, nap_jcart, na_icart) = &
                                   d2ns_bare (m1, m2, current_spin, nah, nap_jcart, na_icart) + d2ns_bare_k    
                                   ! 
                                ENDDO ! m2
                                !
                             ENDDO ! m1
                             !
                          ENDIF
                          !
                       ENDDO ! nah
                       !
                    ENDIF  ! IF  (nap_jcart >= na_icart) 
                    ! 
                 ENDDO ! jcart
                 !
              ENDDO ! nap
              !
           ENDDO ! icart
           !
        ENDDO ! na
        !
        ! Print information about how much the calculation has advanced
        ! because it may take really long...
        !
        !WRITE(stdout,'(4x,i4,1x,"%")') INT(ik*100/nksq)
        WRITE(stdout,'(5x,"k point #",1x,i5,3x,"out of",1x,i5)') ik, nksq
        !
     ENDDO ! ik
     !
#ifdef __MPI
     CALL mp_sum(d2ns_bare, inter_pool_comm) 
#endif
     !
     ! If nspin=1, k point weight is normalized to 2 el/band 
     !
     IF (nspin.EQ.1) d2ns_bare = 0.5d0 * d2ns_bare  
     !   
     ! Filling the m1 m2 lower triangular part (LTP) of the calculated UPT matrices
     !   
     DO na_icart = 1, 3*nat
        DO nap_jcart = na_icart, 3*nat
           DO m1 = 2, ldim
              DO m2 = 1, m1-1
                 DO nah = 1, nat ! the Hubbard atom
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       DO is = 1, nspin
                          d2ns_bare (m1, m2, is, nah, nap_jcart, na_icart) = &
                          d2ns_bare (m2, m1, is, nah, nap_jcart, na_icart)   
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !   
     ! Filling the LTP matrices
     !   
     DO na_icart = 2, 3*nat
        DO nap_jcart = 1, na_icart-1
           d2ns_bare (:,:,:,:, nap_jcart, na_icart) = &
           d2ns_bare (:,:,:,:, na_icart, nap_jcart)   
        ENDDO
     ENDDO
     !
     ! Write d2ns_bare to file
     !
     IF (ionode) WRITE(iund2nsbare,*) d2ns_bare
     !
  ENDIF
  !
  IF (ionode) CLOSE(unit=iund2nsbare,status='KEEP')
  !   
  ! Writing the full d2ns_bare matrix to the standard output file
  !
  IF (iverbosity==1) THEN 
     WRITE(stdout,*)
     WRITE(stdout,*) 'd2ns_bare NOT SYMMETRIZED IN CARTESIAN COORDINATES' 
     DO na = 1, nat  ! the displaced atom
        DO icart = 1, 3
           na_icart = 3*(na-1) + icart
           DO nap = 1, nat ! the displaced atom
              DO jcart = 1, 3
                 nap_jcart = 3*(nap-1) + jcart
                 DO nah = 1, nat ! the Hubbard atom
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       DO is = 1, nspin
                          WRITE(stdout,'(2(a,1x,i2,2x,a,1x,i2,2x))') &
                                  'displaced atom L =', na, 'ipol=', &
                                   icart,'displaced atom Lp =', nap, 'ipol=', jcart
                          WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                          DO m1 = 1, ldim
                             WRITE(stdout,'(14(f15.10,1x))') &
                                  (d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart), m2 = 1, ldim)
                          ENDDO
                       ENDDO
                    ENDIF 
                 ENDDO 
              ENDDO 
           ENDDO  
        ENDDO 
     ENDDO
     WRITE( stdout,*)
  ENDIF
  !
  ! And now we compute the bare part of the dynamical matrix and add it to dyn
  !      
  DO na = 1, nat  ! the displaced atom
     !
     DO icart = 1, 3
        !
        na_icart = 3*(na-1) + icart
        !      
        DO nap = 1, nat ! the displaced atom
           !
           DO jcart = 1, 3
              !
              nap_jcart = 3*(nap-1) + jcart
              !      
              DO nah = 1, nat ! the Hubbard atom
                 !
                 nt = ityp(nah)
                 !
                 ! The effective Hubbard_U part
                 ! 
                 IF (is_hubbard(nt)) THEN
                    !
                    ! work = \sum_{is, m1, m2} [ 
                    ! (0.5*delta_m1m2- ns(m1,m2,is,nah))*d2ns_bare - 
                    ! - conjg(dnsbare ( m1, m2, is, nah, icart,na))*
                    !         dnsbare ( m1, m2, is, nah, jcart,nap) ]
                    !      
                    work = (0.d0, 0.d0)
                    !       
                    DO is = 1, nspin
                       !      
                       DO m1 = 1, 2*Hubbard_l(nt) + 1
                          !
                          DO m2 = 1, 2*Hubbard_l(nt) + 1
                             !      
                             nsaux = rho%ns(m1, m2, is, nah)
                             !      
                             ! dnsbare matrix is symmetric i.e. dnsbare(m1,m2) = dnsbare(m2,m1) 
                             ! (when keeping only the terms in u and neglecting the ones in u*)
                             !      
                             dnsaux1 = dnsbare(m1, m2, is, nah, icart, na)
                             dnsaux2 = dnsbare(m1, m2, is, nah, jcart, nap)
                             d2ns_bare_aux = d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart) 
                             !                              
                             ! Include the delta_m1m2 contribution              
                             !
                             IF (m1==m2) work = work + 0.5d0*d2ns_bare_aux
                             !      
                             work = work - nsaux   * d2ns_bare_aux &
                                         - dnsaux1 * CONJG(dnsaux2)
                             ! 
                          ENDDO ! m2
                          !
                       ENDDO ! m1
                       !      
                    ENDDO ! is
                    ! 
                    ! Summing the contributions for each Hubbard atom                    
                    !
                    dynwrk(na_icart, nap_jcart) = dynwrk(na_icart, nap_jcart) & 
                                                  + work * effU(nt)
                    !      
                 ENDIF
                 !
                 ! The Hubbard_J0 part (if present)
                 !
                 IF (Hubbard_J0(nt).NE.0.d0) THEN
                    !                     
                    ! work = \sum_{is, m1, m2} [ 
                    !    ns(m1,m2,-is,nah) * d2ns_bare - 
                    !    dnsbare( m1, m2, is, nah, icart,na) *
                    !    conjg(dnsbare( m1, m2, -is, nah, jcart,nap)) ]
                    !      
                    work = (0.d0, 0.d0)
                    !      
                    DO is = 1, nspin
                       !      
                       IF ((is==1) .AND. (nspin==2)) THEN
                          isi = 2
                       ELSE 
                          isi = 1
                       ENDIF
                       !
                       DO m1 = 1, 2*Hubbard_l(nt) + 1
                          !
                          DO m2 = 1, 2*Hubbard_l(nt) + 1
                             !      
                             nsaux = rho%ns(m1, m2, isi, nah)
                             !      
                             ! dnsbare matrix is symmetric i.e. dnsbare(m1,m2) = dnsbare(m2,m1) 
                             ! (when keeping only the terms in u and neglecting the ones in u*)
                             !      
                             dnsaux1 = dnsbare(m1, m2, is,  nah, icart, na)
                             dnsaux2 = dnsbare(m1, m2, isi, nah, jcart, nap) 
                             d2ns_bare_aux = d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart) 
                             !                              
                             ! DO NOT include the delta_m1m2 contribution  
                             ! Note the sign change            
                             !
                             work = work + nsaux *  d2ns_bare_aux + &
                                         + dnsaux1 * CONJG(dnsaux2) 
                             ! 
                          ENDDO ! m2
                          !
                       ENDDO ! m1
                       !      
                    ENDDO ! is
                    ! 
                    ! Summing the contributions for each Hubbard atom                    
                    !
                    dynwrk(na_icart, nap_jcart) = dynwrk(na_icart, nap_jcart) & 
                                                  + work * Hubbard_J0(nt)
                    !
                 ENDIF
                 !
              ENDDO ! nah
              !
           ENDDO ! jcart
           !
        ENDDO  ! nap
        !
     ENDDO ! icart
     !
  ENDDO ! na
  !      
  IF (nspin.EQ.1) dynwrk = 2.d0 * dynwrk  
  !      
  ! Transform from the Cartesian to the pattern basis u. 
  ! We rotate the dynamical matrix on the basis of patterns u.
  !
  DO imode = 1, nmodes
     DO jmode = 1, nmodes
        !
        work = (0.d0, 0.d0)
        !
        DO nap_jcart = 1, 3*nat
           DO na_icart = 1, 3*nat
              work = work + CONJG (u(na_icart, imode) ) * &
                            dynwrk(na_icart, nap_jcart) * &
                            u(nap_jcart, jmode)
           ENDDO
        ENDDO
        !
        ! Adding the contribution to the total dynamical matrix dyn
        !        
        dyn(imode, jmode) = dyn(imode, jmode) + work
        !     
        ! Keep separately the bare part of the dynamical matrix
        ! 
        dyn_hub_bare(imode, jmode) = dyn_hub_bare(imode, jmode) + work
        !
     ENDDO
  ENDDO
  !
  IF (iverbosity==1) THEN
     !
     ! Write the UNSYMMETRIZED bare Hubbard dynamical matrix 
     ! in the pattern basis
     !
     CALL tra_write_matrix_no_sym ('dyn_hub_bare NOT SYMMETRIZED',dyn_hub_bare,nat) 
     !
     ! Write the SYMMETRIZED bare Hubbard dynamical matrix 
     ! in carthesian coordinates
     !
     CALL tra_write_matrix ('dyn_hub_bare SYMMETRIZED',dyn_hub_bare,u,nat)                   
     !
  ENDIF
  !
  DEALLOCATE (proj1)
  DEALLOCATE (projpb)
  DEALLOCATE (projpdb)
  DEALLOCATE (d2ns_bare)
  DEALLOCATE (dynwrk)     
  !
  CALL stop_clock ( 'dynmat_hub_bare' )
  ! 
  RETURN
  !      
END SUBROUTINE dynmat_hub_bare
      
