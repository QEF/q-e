!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE new_ns()
  !-----------------------------------------------------------------------
  !
  ! This routine computes the new value for ns (the occupation numbers of
  ! ortogonalized atomic wfcs).
  ! These quantities are defined as follows: ns_{I,s,m1,m2} = \sum_{k,v}
  ! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : lchi, nchi, oc
  USE ions_base,            ONLY : nat, ityp
  USE basis,                ONLY : natomwfc
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : ns, nsnew, Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_alpha, swfcatom, &
                                   eth, d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symme,                ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg, gamma_only
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat

  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, n, counter, m1, m2, &
       m0, m00, l, ldim
  INTEGER, ALLOCATABLE ::  offset (:)
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  ! offset of d electrons of atom d
  ! in the natomwfc ordering
  REAL(DP) , ALLOCATABLE :: nr (:,:,:,:)
  REAL(DP) ::  t0, scnds
  ! cpu time spent

  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP) :: ZDOTC
  COMPLEX(DP) , ALLOCATABLE :: proj(:,:)

  REAL(DP) :: psum

  t0 = scnds ()  
  ldim = 2 * Hubbard_lmax + 1
  ALLOCATE( offset(nat), proj(natomwfc,nbnd), nr(ldim,ldim,nspin,nat) )  
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  counter = 0  
  DO na = 1, nat  
     nt = ityp (na)  
     DO n = 1, nchi (nt)  
        IF (oc (n, nt) >= 0.d0) THEN  
           l = lchi (n, nt)  
           IF (l == Hubbard_l(nt)) offset (na) = counter  
           counter = counter + 2 * l + 1  
        ENDIF
     ENDDO

  ENDDO

  IF (counter.NE.natomwfc) CALL errore ('new_ns', 'nstart<>counter', 1)
  nr    (:,:,:,:) = 0.d0
  nsnew (:,:,:,:) = 0.d0
  !
  !    we start a loop on k points
  !

  IF (nks.GT.1) REWIND (iunigk)

  DO ik = 1, nks
     IF (lsda) current_spin = isk(ik)
     npw = ngk (ik)
     IF (nks > 1) THEN
        READ (iunigk) igk
        CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     END IF
     CALL davcio (swfcatom, nwordatwfc, iunsat, ik, - 1)
     !
     ! make the projection
     !
        DO ibnd = 1, nbnd
           DO i = 1, natomwfc
              IF ( gamma_only ) THEN 
                 proj (i, ibnd) = 2.d0 * &
                      DDOT(2*npw, swfcatom (1, i), 1, evc (1, ibnd), 1) 
                 IF (gstart.EQ.2) proj (i, ibnd) = proj (i, ibnd) - &
                      swfcatom (1, i) * evc (1, ibnd)
              ELSE 
                 proj (i, ibnd) = ZDOTC (npw, swfcatom (1, i), 1, evc (1, ibnd), 1)
              ENDIF
           ENDDO
        ENDDO

#ifdef __PARA
        CALL reduce (2 * natomwfc * nbnd, proj)
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
                         wg(ibnd,ik) *  DBLE( proj(offset(na)+m2,ibnd) * &
                         CONJG(proj(offset(na)+m1,ibnd)) )
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
    ENDDO
! on k-points

  ENDDO
#ifdef __PARA
  CALL poolreduce (ldim * ldim * nspin * nat , nr)  
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

  ! symmetryze the quantities nr -> nsnew
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
                             nsnew(m1,m2,is,na) = nsnew(m1,m2,is,na) +  &
                                   nr(m0,m00,is,nb) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.1) THEN
                             nsnew(m1,m2,is,na) = nsnew(m1,m2,is,na) +  &
                                   d1(m0 ,m1,isym) * nr(m0,m00,is,nb) * &
                                   d1(m00,m2,isym) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.2) THEN
                             nsnew(m1,m2,is,na) = nsnew(m1,m2,is,na) +  &
                                   d2(m0 ,m1,isym) * nr(m0,m00,is,nb) * &
                                   d2(m00,m2,isym) / nsym
                          ELSE IF (Hubbard_l(nt).EQ.3) THEN
                             nsnew(m1,m2,is,na) = nsnew(m1,m2,is,na) +  &
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
                 psum = ABS ( nsnew(m1,m2,is,na) - nsnew(m1,m2,is,na) )  
                 IF (psum.GT.1.d-10) THEN  
                    WRITE( stdout, * ) na, is, m1, m2  
                    WRITE( stdout, * ) nsnew (m1, m2, is, na)  
                    WRITE( stdout, * ) nsnew (m2, m1, is, na)  
                    CALL errore ('new_ns', 'non hermitean matrix', 1)  
                 ELSE  
                    nsnew(m1,m2,is,na) = 0.5d0 * (nsnew(m1,m2,is,na) + &
                                                  nsnew(m2,m1,is,na) )
                    nsnew(m2,m1,is,na) = nsnew(m1,m2,is,na)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF 
  ENDDO

  !
  ! Now the contribution to the total energy is computed. The corrections
  ! needed to obtain a variational expression are already included
  !
  eth = 0.d0  
  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 eth = eth + Hubbard_U(nt) * nsnew(m1,m2,is,na) * &
                          (ns(m2,m1,is,na) - nsnew(m2,m1,is,na) * 0.5d0)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  DEALLOCATE ( offset, proj, nr )
  IF (nspin.EQ.1) eth = 2.d0 * eth

  RETURN

END SUBROUTINE new_ns
