!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_nsb( ns )
  !-----------------------------------------------------------------------
  !! This routine computes the new value for ns (the occupation numbers of
  !! ortogonalized atomic wfcs) for background states.
  !! These quantities are defined as follows: 
  !! 
  !! ns_{I,s,m1,m2} = \sum_{k,v}
  !! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !!
  !! ns is symmetric (and real for collinear case, due to time reversal),
  !! hence the order of m1 and m2 does not matter.
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : Hubbard_l, q_ae, wfcU, U_projection, backall, &
                                   is_hubbard_back, nwfcU, offsetU_back, ldmx_b, &
                                   ldim_back, Hubbard_l_back, Hubbard_l1_back,   &
                                   offsetU_back1
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type

  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: ns(ldmx_b,ldmx_b,nspin,nat)
  !
  TYPE (bec_type) :: proj     
  ! proj(nwfcU,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m0, m00, ldim
  INTEGER :: off, off1, off2, off3, m11, m22, ldim_std
  INTEGER, ALLOCATABLE :: lhub(:)
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  REAL(DP) , ALLOCATABLE :: nr (:,:,:,:)
  REAL(DP) :: psum
  !
  CALL start_clock('new_nsb')
  !
  ALLOCATE( nr(ldmx_b, ldmx_b, nspin, nat) )  
  !
  CALL allocate_bec_type ( nwfcU, nbnd, proj ) 
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  nr (:,:,:,:) = 0.d0
  !
  ! we start a loop on k points
  !
  DO ik = 1, nks
     ! 
     IF (lsda) current_spin = isk(ik)
     !
     npw = ngk(ik)
     !
     IF (nks > 1) CALL get_buffer (evc, nwordwfc, iunwfc, ik)
     !
     ! make the projection
     !
     IF ( U_projection == 'pseudo' ) THEN
        CALL errore('new_nsb', 'U_projection = pseudo is not supported',1)
        !CALL compute_pproj( ik, q_ae, proj )
     ELSE
        IF (nks > 1) CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
        CALL calbec ( npw, wfcU, evc, proj )
     END IF
     !
     ! compute the occupation matrix (ns_{I,s,m1,m2}) of the
     ! atomic orbitals for the background states
     !
     DO na = 1, nat  
        nt = ityp (na)  
        IF ( is_hubbard_back(nt) ) THEN  
           DO m1 = 1, ldim_back(nt)  
              off1 = offsetU_back(na)+m1
              IF (backall(nt) .AND. m1.GT.2*Hubbard_l_back(nt)+1) THEN
                 off1 = offsetU_back1(na)+m1-2*Hubbard_l_back(nt)-1
              ENDIF
              DO m2 = m1, ldim_back(nt) 
                 off2 = offsetU_back(na)+m2
                 IF (backall(nt) .AND. m2.GT.2*Hubbard_l_back(nt)+1) THEN
                    off2 = offsetU_back1(na)+m2-2*Hubbard_l_back(nt)-1
                 ENDIF
                 IF ( gamma_only ) THEN
                    DO ibnd = 1, nbnd  
                       nr(m1,m2,current_spin,na) = nr(m1,m2,current_spin,na) + &
                            proj%r(off2,ibnd) * &
                            proj%r(off1,ibnd) * wg(ibnd,ik) 
                    ENDDO
                 ELSE
                    DO ibnd = 1, nbnd  
                       nr(m1,m2,current_spin,na) = nr(m1,m2,current_spin,na) + &
                            DBLE( proj%k(off2,ibnd) * &
                            CONJG(proj%k(off1,ibnd)) ) * wg(ibnd,ik) 
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     ! 
  ENDDO
  !
  CALL deallocate_bec_type (proj) 
  !
  CALL mp_sum( nr, inter_pool_comm )
  !
  IF (nspin.EQ.1) nr = 0.5d0 * nr
  !
  ! impose hermiticity of n_{m1,m2}
  !
  DO na = 1, nat  
     nt = ityp(na)
     DO is = 1, nspin  
        DO m1 = 1, ldim_back(nt) 
           DO m2 = m1 + 1, ldim_back(nt)
              nr(m2,m1,is,na) = nr(m1,m2,is,na)  
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! symmetrize the quantities nr -> ns
  !
  ns (:,:,:,:) = 0.d0
  !
  ALLOCATE (lhub(ldmx_b))
  !
  DO na = 1, nat  
     nt = ityp (na)  
     IF ( is_hubbard_back(nt) ) THEN 
        ! 
        DO m1 = 1, ldim_back(nt)
           lhub(m1) = Hubbard_l_back(nt)
           IF (backall(nt) .AND. m1.GT.(2*Hubbard_l_back(nt)+1)) THEN
              lhub(m1) = Hubbard_l1_back(nt)
           ENDIF
        ENDDO
        !
        ldim_std = 2*Hubbard_l(nt)+1
        !
        DO is = 1, nspin  
           off = 0 !2*Hubbard_l(nt)+1
           off2 = 2*Hubbard_l_back(nt)+1 !2*(Hubbard_l(nt)+Hubbard_l_back(nt)+1)
           DO m1 = 1, ldim_back(nt)
              m11 = m1
              !
              IF (backall(nt) .AND. m1.GT.(2*Hubbard_l_back(nt)+1)) THEN
                 off = 2*Hubbard_l_back(nt)+1 !2*(Hubbard_l(nt)+Hubbard_l_back(nt)+1)
                 off2 = 2*(Hubbard_l_back(nt)+Hubbard_l1_back(nt)+1) !ldim_u(nt)
                 m11 = m1-2*Hubbard_l_back(nt)-1
              ENDIF
              !
              off1 = 0 !2*Hubbard_l(nt)+1
              off3 = 2*Hubbard_l_back(nt)+1 !2*(Hubbard_l(nt)+Hubbard_l_back(nt)+1)
              !
              DO m2 = 1, ldim_back(nt) 
                 m22 = m2
                 !
                 IF (backall(nt) .AND. m2.GT.(2*Hubbard_l_back(nt)+1)) THEN
                    off1 = 2*Hubbard_l_back(nt)+1 !2*(Hubbard_l(nt)+Hubbard_l_back(nt)+1)
                    off3 = 2*(Hubbard_l_back(nt)+Hubbard_l1_back(nt)+1) !ldim_u(nt)
                    m22 = m2 - 2*Hubbard_l_back(nt)-1
                 ENDIF
                 !
                 DO isym = 1, nsym  
                    nb = irt (isym, na)  
                    DO m0 = 1, 2*lhub(m1)+1     ! 2 * Hubbard_l(nt) + 1  
                       DO m00 = 1, 2*lhub(m2)+1 ! 2 * Hubbard_l(nt) + 1  
                          !
                          if (lhub(m1).eq.0.and.lhub(m2).eq.0) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                  nr(m0+off,m00+off1,is,nb) / nsym
                          else if (lhub(m1).eq.0.and.lhub(m2).eq.1) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                  nr(m0+off,m00+off1,is,nb) * &
                                  d1(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.0.and.lhub(m2).eq.2) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                  nr(m0+off,m00+off1,is,nb) * &
                                  d2(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.0.and.lhub(m2).eq.3) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                  nr(m0+off,m00+off1,is,nb) * &
                                  d3(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.1.and.lhub(m2).eq.0) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                             d1(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) / nsym
                          else if (lhub(m1).eq.1.and.lhub(m2).eq.1) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d1(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d1(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.1.and.lhub(m2).eq.2) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d1(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d2(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.1.and.lhub(m2).eq.3) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d1(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d3(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.2.and.lhub(m2).eq.0) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                             d2(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) / nsym
                          else if (lhub(m1).eq.2.and.lhub(m2).eq.1) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d2(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d1(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.2.and.lhub(m2).eq.2) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d2(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d2(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.2.and.lhub(m2).eq.3) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d2(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d3(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.3.and.lhub(m2).eq.0) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                             d3(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) / nsym
                          else if (lhub(m1).eq.3.and.lhub(m2).eq.1) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d3(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d1(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.3.and.lhub(m2).eq.2) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d3(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d2(m00,m22,isym) / nsym
                          else if (lhub(m1).eq.3.and.lhub(m2).eq.3) then
                             ns(m1,m2,is,na) = &
                                  ns(m1,m2,is,na) +  &
                                d3(m0,m11,isym)*nr(m0+off,m00+off1,is,nb) * &
                                  d3(m00,m22,isym) / nsym
                          else
                             CALL errore ('new_ns', &
                                         'angular momentum not implemented', &
                                          ABS(Hubbard_l(nt)) )
                          endif
                          !
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  DEALLOCATE ( lhub )
  DEALLOCATE ( nr )
  !
  ! Now we make the matrix ns(m1,m2) strictly hermitean
  !
  DO na = 1, nat  
     nt = ityp (na)  
     IF ( is_hubbard_back(nt) ) THEN  
        DO is = 1, nspin  
           DO m1 = 1, ldim_back(nt)  
              DO m2 = 1, ldim_back(nt) 
                 psum = ABS ( ns(m1,m2,is,na) - ns(m2,m1,is,na) )  
                 IF (psum.GT.1.d-10) THEN  
                    WRITE( stdout, * ) na, is, m1, m2  
                    WRITE( stdout, * ) ns (m1, m2, is, na)  
                    WRITE( stdout, * ) ns (m2, m1, is, na)  
                    CALL errore ('new_nsb', 'non hermitean matrix', 1)  
                 ELSE  
                    ns(m1,m2,is,na) = 0.5d0 * (ns(m1,m2,is,na) + &
                                               ns(m2,m1,is,na) )
                    ns(m2,m1,is,na) = ns(m1,m2,is,na)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF 
  ENDDO
  !
  CALL stop_clock('new_nsb')
  !
  RETURN
  !
END SUBROUTINE new_nsb
