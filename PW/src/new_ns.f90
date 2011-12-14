!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_ns(ns)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the new value for ns (the occupation numbers of
  ! ortogonalized atomic wfcs).
  ! These quantities are defined as follows: ns_{I,s,m1,m2} = \sum_{k,v}
  ! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE basis,                ONLY : natomwfc
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, oatwfc, &
                                   Hubbard_U, Hubbard_alpha, swfcatom
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum


  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)

  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, &
       m0, m00, ldim
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  ! in the natomwfc ordering
  REAL(DP) , ALLOCATABLE :: nr (:,:,:,:)
  REAL(DP) ::  t0, scnds
  ! cpu time spent

  REAL(DP), EXTERNAL :: ddot
  COMPLEX(DP) :: zdotc
  COMPLEX(DP) , ALLOCATABLE :: proj(:,:)

  REAL(DP) :: psum

  t0 = scnds ()  
  ldim = 2 * Hubbard_lmax + 1
  ALLOCATE( proj(natomwfc,nbnd), nr(ldim,ldim,nspin,nat) )  
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
  !
  nr (:,:,:,:) = 0.d0
  ns (:,:,:,:) = 0.d0
  !
  !    we start a loop on k points
  !

  IF (nks > 1) REWIND (iunigk)

  DO ik = 1, nks
     IF (lsda) current_spin = isk(ik)
     npw = ngk (ik)
     IF (nks > 1) THEN
        READ (iunigk) igk
        CALL get_buffer  (evc, nwordwfc, iunwfc, ik)
     END IF
     CALL davcio (swfcatom, nwordatwfc, iunsat, ik, - 1)
     !
     ! make the projection
     !
        DO ibnd = 1, nbnd
           DO i = 1, natomwfc
              IF ( gamma_only ) THEN 
                 proj (i, ibnd) = 2.d0 * &
                      ddot(2*npw, swfcatom (1, i), 1, evc (1, ibnd), 1) 
                 IF (gstart.EQ.2) proj (i, ibnd) = proj (i, ibnd) - &
                      swfcatom (1, i) * evc (1, ibnd)
              ELSE 
                 proj (i, ibnd) = zdotc (npw, swfcatom (1, i), 1, evc (1, ibnd), 1)
              ENDIF
           ENDDO
        ENDDO

#ifdef __PARA
        CALL mp_sum ( proj, intra_pool_comm )
#endif

     !
     ! compute the occupation numbers (the quantities n(m1,m2)) of the
     ! atomic orbitals
     !
     DO na = 1, nat  
        nt = ityp (na)  
        IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = m1, 2 * Hubbard_l(nt) + 1
                 DO ibnd = 1, nbnd  
                    nr(m1,m2,current_spin,na) = nr(m1,m2,current_spin,na) + &
                         wg(ibnd,ik) *  DBLE( proj(oatwfc(na)+m2,ibnd) * &
                         CONJG(proj(oatwfc(na)+m1,ibnd)) )
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
    ENDDO
! on k-points

  ENDDO
#ifdef __PARA
  CALL mp_sum( nr, inter_pool_comm )
#endif
  IF (nspin.EQ.1) nr = 0.5d0 * nr
  !
  ! impose hermiticity of n_{m1,m2}
  !
  DO na = 1, nat  
     nt = ityp(na)
     DO is = 1, nspin  
        DO m1 = 1, 2 * Hubbard_l(nt) + 1
           DO m2 = m1 + 1, 2 * Hubbard_l(nt) + 1  
              nr (m2, m1, is, na) = nr (m1, m2, is, na)  
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! symmetrize the quantities nr -> ns
  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 DO isym = 1, nsym  
                    nb = irt (isym, na)  
                    DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                       DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                          IF (Hubbard_l(nt).EQ.0) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                   nr(m0,m00,is,nb) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.1) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                   d1(m0 ,m1,isym) * nr(m0,m00,is,nb) * &
                                   d1(m00,m2,isym) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.2) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                   d2(m0 ,m1,isym) * nr(m0,m00,is,nb) * &
                                   d2(m00,m2,isym) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.3) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                   d3(m0 ,m1,isym) * nr(m0,m00,is,nb) * &
                                   d3(m00,m2,isym) / nsym
                          ELSE
                             CALL errore ('new_ns', &
                                         'angular momentum not implemented', &
                                          ABS(Hubbard_l(nt)) )
                          END IF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  ! Now we make the matrix ns(m1,m2) strictly hermitean
  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = m1, 2 * Hubbard_l(nt) + 1  
                 psum = ABS ( ns(m1,m2,is,na) - ns(m2,m1,is,na) )  
                 IF (psum.GT.1.d-10) THEN  
                    WRITE( stdout, * ) na, is, m1, m2  
                    WRITE( stdout, * ) ns (m1, m2, is, na)  
                    WRITE( stdout, * ) ns (m2, m1, is, na)  
                    CALL errore ('new_ns', 'non hermitean matrix', 1)  
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

  DEALLOCATE ( proj, nr )

  RETURN

END SUBROUTINE new_ns
