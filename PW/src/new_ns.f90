!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
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

  !  It seems that the order of {m1, m2} in the definition should be opposite. 
  !  Hovewer, since ns is symmetric (and real for collinear case, due to time reversal) 
  !  it does not matter.  
  !  (A.Smogunov) 


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
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type

  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !
  TYPE (bec_type) :: proj     ! proj(natomwfc,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m0, m00, ldim
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  ! in the natomwfc ordering
  REAL(DP) , ALLOCATABLE :: nr (:,:,:,:)
  REAL(DP) :: psum

  CALL start_clock('new_ns')
  ldim = 2 * Hubbard_lmax + 1
  ALLOCATE( nr(ldim,ldim,nspin,nat) )  
  CALL allocate_bec_type ( natomwfc, nbnd, proj ) 
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
     CALL calbec ( npw, swfcatom, evc, proj )
     !
     ! compute the occupation numbers (the quantities n(m1,m2)) of the
     ! atomic orbitals
     !
     DO na = 1, nat  
        nt = ityp (na)  
        IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = m1, 2 * Hubbard_l(nt) + 1
                 IF ( gamma_only ) THEN
                    DO ibnd = 1, nbnd  
                       nr(m1,m2,current_spin,na) = nr(m1,m2,current_spin,na) + &
                            proj%r(oatwfc(na)+m2,ibnd) * &
                            proj%r(oatwfc(na)+m1,ibnd) * wg(ibnd,ik) 
                    ENDDO
                 ELSE
                    DO ibnd = 1, nbnd  
                       nr(m1,m2,current_spin,na) = nr(m1,m2,current_spin,na) + &
                            DBLE( proj%k(oatwfc(na)+m2,ibnd) * &
                            CONJG(proj%k(oatwfc(na)+m1,ibnd)) ) * wg(ibnd,ik) 
                    ENDDO
                 END IF
              ENDDO
           ENDDO
        ENDIF
    ENDDO
! on k-points

  ENDDO
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

  DEALLOCATE ( nr )
  CALL stop_clock('new_ns')

  RETURN

END SUBROUTINE new_ns 


SUBROUTINE new_ns_nc(ns)
  !-----------------------------------------------------------------------
  !
  ! Noncollinear version (A. Smogunov). 
  !

  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE basis,                ONLY : natomwfc
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, oatwfc, &
                                   Hubbard_U, Hubbard_alpha, swfcatom, d_spin_ldau
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE noncollin_module, ONLY : noncolin, npol
  USE symm_base,            ONLY : nsym, irt, time_reversal, t_rev
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : intra_bgrp_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum

  IMPLICIT NONE
  !
  ! I/O variables
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  INTEGER :: ik, ibnd, is, js, i, j, sigmay2, na, nb, nt, isym,  &
             m1, m2, m3, m4, is1, is2, is3, is4, m0, m00, ldim

  COMPLEX(DP) , ALLOCATABLE :: nr (:,:,:,:,:), nr1 (:,:,:,:,:), proj(:,:) 

  COMPLEX(DP) :: z, zdotc

  REAL(DP) :: psum


  CALL start_clock('new_ns')
  ldim = 2 * Hubbard_lmax + 1

  ALLOCATE( nr(ldim,ldim,npol,npol,nat), nr1(ldim,ldim,npol,npol,nat) )  
  ALLOCATE( proj(natomwfc,nbnd) )

  nr  (:,:,:,:,:) = 0.d0
  nr1 (:,:,:,:,:) = 0.d0
  ns  (:,:,:,:)   = 0.d0

!--
!    loop on k points
!
  IF (nks > 1) REWIND (iunigk)
  DO ik = 1, nks

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
         proj(i, ibnd) = zdotc (npwx*npol, swfcatom (1, i), 1, evc (1, ibnd), 1)
       ENDDO
     ENDDO
#ifdef __MPI
     CALL mp_sum ( proj, intra_bgrp_comm )
#endif
     !
     ! compute the occupation matrix
     !
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0) then  
          ldim = 2 * Hubbard_l(nt) + 1

          do m1 = 1, 2 * Hubbard_l(nt) + 1
            do m2 = 1, 2 * Hubbard_l(nt) + 1
              do is1 = 1, npol
                do is2 = 1, npol

                  do ibnd = 1, nbnd
                    nr(m1,m2,is1,is2,na) = nr(m1,m2,is1,is2,na) + &
                            wg(ibnd,ik) * CONJG( proj(oatwfc(na)+m1+ldim*(is1-1),ibnd) ) * &
                                           proj(oatwfc(na)+m2+ldim*(is2-1),ibnd)
                  enddo

                enddo
              enddo
            enddo
          enddo

        endif
     enddo

  enddo
!---

  CALL mp_sum( nr, inter_pool_comm )

!--  symmetrize: nr  -->  nr1
!
  do na = 1, nat  
    nt = ityp (na)  
    if (Hubbard_U(nt).ne.0.d0) then  

      do m1 = 1, 2 * Hubbard_l(nt) + 1  
        do m2 = 1, 2 * Hubbard_l(nt) + 1  
          do is1 = 1, npol
            do is2 = 1, npol

loopisym:     do isym = 1, nsym  
                nb = irt (isym, na)  

                do m3 = 1, 2 * Hubbard_l(nt) + 1
                  do m4 = 1, 2 * Hubbard_l(nt) + 1
                    do is3 = 1, npol
                      do is4 = 1, npol
  
                        if (Hubbard_l(nt).eq.0) then
                          if (t_rev(isym).eq.1) then
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*                &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  
                          else
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*                &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  
                          endif
                        elseif (Hubbard_l(nt).eq.1) then
                          if (t_rev(isym).eq.1) then
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d1(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d1(m2,m4,isym)
                          else
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d1(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d1(m2,m4,isym)
                          endif
                        elseif (Hubbard_l(nt).eq.2) then
                          if (t_rev(isym).eq.1) then
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d2(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d2(m2,m4,isym)
                          else
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d2(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d2(m2,m4,isym)
                          endif
                        elseif (Hubbard_l(nt).eq.3) then
                          if (t_rev(isym).eq.1) then
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d3(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d3(m2,m4,isym)
                          else
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d3(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d3(m2,m4,isym)
                          endif
                        else
                          CALL errore ('new_ns', &
                                         'angular momentum not implemented', &
                                          ABS(Hubbard_l(nt)) )
                        endif

                      enddo
                    enddo
                  enddo
                enddo

              enddo  loopisym 

            enddo
          enddo
        enddo
      enddo 

    endif
  enddo
!--

!-- Setup the output matrix ns with combined spin index 
!
  DO na = 1, nat
     nt = ityp (na)
     IF (Hubbard_U(nt).NE.0.d0) THEN
        DO is1 = 1, npol
         do is2 = 1, npol
           i = npol*(is1-1) + is2
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                ns(m1,m2,i,na) = nr1(m1,m2,is1,is2,na)
              ENDDO
           ENDDO
         enddo
        ENDDO
     ENDIF
  ENDDO
!--

!-- make the matrix ns strictly hermitean
!
  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
        DO is1 = 1, npol  
         do is2 = 1, npol
           i = npol*(is1-1) + is2
           j = is1 + npol*(is2-1)
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 psum = ABS ( ns(m1,m2,i,na) - CONJG(ns(m2,m1,j,na)) )  
                 IF (psum.GT.1.d-10) THEN  
                    WRITE( stdout, * ) na, m1, m2, is1, is2  
                    WRITE( stdout, * ) ns (m1, m2, i, na)  
                    WRITE( stdout, * ) ns (m2, m1, j, na)  
                    CALL errore ('new_ns', 'non hermitean matrix', 1)  
                 ELSE  
                    ns (m2, m1, j, na) = CONJG( ns(m1, m2, i, na))
                 ENDIF
              ENDDO
           ENDDO
         enddo
        ENDDO
     ENDIF 
  ENDDO
!--


  DEALLOCATE ( nr, nr1 )
  CALL stop_clock('new_ns')

  RETURN

END SUBROUTINE new_ns_nc



