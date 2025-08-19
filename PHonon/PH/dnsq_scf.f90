!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dnsq_scf (npe, lmetq0)
  !-----------------------------------------------------------------------
  !! DFPT+U: This routine calculates, for each SCF iteration, 
  !! the SCF variation of the occupation matrix ns, for npe perturbations.
  !! The result is stored in the variable dnsscf:
  !
  !! $$ \text{dnsscf}(m1,m2,\text{ispin},I,\text{ipert}) = 
  !!   \sum_{k,n} [ \langle \psi(n,k,\text{ispin})| S_{k}\phi(k,I,m1)\rangle  
  !!                \cdot \langle S_{k+q}\phi(k+q,I,m2)|
  !!                       d\psi(\text{ipert},n,k+q,\text{ispin})\rangle
  !!                + \langle\psi(n,k,\text{ispin})| S_{k}\phi(k,I,m2)\rangle  
  !!                \cdot \langle S_{k+q}\phi(k+q,I,m1)|
  !!                 d\psi(\text{ipert},n,k+q,\text{ispin})\rangle ] $$
  !
  !! Written  by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc, iudwf, lrdwf
  USE ions_base,     ONLY : nat, ityp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, offsetU, nwfcU
  USE ldaU_ph,       ONLY : proj1, proj2
  USE ldaU_lr,       ONLY : swfcatomk, swfcatomkpq, dnsscf, lr_has_dnsorth, lr_dnsorth
  USE klist,         ONLY : wk, degauss, ngauss, ngk, degauss_cond
  USE wvfct,         ONLY : npwx, nbnd, et, nbnd_cond
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma, nbnd_occ
  USE units_lr,      ONLY : iuatswfc
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE ener,          ONLY : ef, ef_cond
  USE uspp,          ONLY : okvan 
  USE mp_pools,      ONLY : inter_pool_comm 
  USE mp_bands,      ONLY : intra_bgrp_comm   
  USE mp,            ONLY : mp_sum
  USE io_global,     ONLY : stdout
  USE buffers,       ONLY : get_buffer
  USE efermi_shift,  ONLY : def
  USE lr_two_chem,   ONLY : def_val,def_cond
  USE two_chem,      ONLY : twochem
  USE control_flags, ONLY : iverbosity
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: npe 
  !! the number of perturbations
  LOGICAL,  INTENT(IN) :: lmetq0
  !! TRUE if \(xq=(0,0,0)\) in a metal
  !
  ! ... local variables
  !
  INTEGER  :: nt, ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, is, &
              ipert, nrec, ldim, npw, npwq, ik , ikk, ikq
  COMPLEX(DP), ALLOCATABLE :: dpsi(:,:)
  REAL(DP) :: weight, wdelta, w1
  REAL(DP), EXTERNAL :: w0gauss 
  !
  CALL start_clock( 'dnsq_scf' )
  ! 
  ALLOCATE (dpsi(npwx,nbnd))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  !
  ldim = 2 * Hubbard_lmax + 1
  !   
  ! At each iteration dnsscf is set to zero, as it is recomputed
  ! 
  dnsscf = (0.d0, 0.d0) 
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
     ! Read the atomic orbitals S*\phi at k and k+q from file (unit iuatswfc)
     !  
     CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
     IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
     !
     DO ipert = 1, npe
        !   
        nrec = (ipert - 1) * nksq + ik
        !
        ! At each SCF iteration for each ik and ipert read dpsi on iunit iudwf
        ! iudwf contains data for the actual perturbation ipert.
        !
        CALL get_buffer (dpsi, lrdwf, iudwf, nrec)
        !
        ! Calculate:
        ! proj1 (ibnd, ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
        ! proj2 (ibnd, ihubst) = < S_{k+q}\phi_(k+q,I,m)) | dpsi(ipert,ibnd,k+q) >
        ! 
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           IF (is_hubbard(nt)) THEN
              !   
              DO m = 1, 2*Hubbard_l(nt)+1
                 !
                 ihubst = offsetU(nah) + m   ! I m index
                 !
                 DO ibnd = 1, nbnd_occ(ikk)
                    !
                    proj1(ibnd,ihubst) = dot_product (swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
                    proj2(ibnd,ihubst) = dot_product (swfcatomkpq(1:npwq,ihubst), dpsi(1:npwq,ibnd))
                    !
                 ENDDO
                 ! 
              ENDDO ! m
              !
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum (proj1, intra_bgrp_comm)  
        CALL mp_sum (proj2, intra_bgrp_comm)
        !
        ! Calculate dnsscf. It is in the pattern basis, because dpsi is.                       
        ! The weights (theta functions) are already contained in dpsi, 
        ! according to Eq. (75) in Rev. Mod. Phys. 73, 515, (2001).
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
                    DO ibnd = 1, nbnd_occ(ikk)
                       !    
                       dnsscf(m1,m2,current_spin,nah,ipert) = dnsscf(m1,m2,current_spin,nah,ipert) + &
                                        wk(ikk) * ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2) &
                                                  + CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
                       !
                       ! Correction for metals at q=0 to be added when dpsi is NOT converged.
                       ! After the convergence the correction is built into 
                       ! dpsi by ef_fermi.f90 with flag==.true. (called in solve_linter)
                       ! wdelta is a smeared delta function at the Fermi level.  
                       !                       
                       IF (lmetq0) THEN
                          !
                          wdelta = w0gauss ( (ef-et(ibnd,ikk))/degauss, ngauss) / degauss
                          weight = wk(ikk)  
                          w1 = weight * wdelta   
                          !
                          IF (twochem) THEN
                           !two chem case, valence and condcution must be treated separately for the two fermi shifts
                             IF (ibnd.le.nbnd-nbnd_cond) then
                                  dnsscf(m1,m2,current_spin,nah,ipert) = &
                                  dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def_val(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                             ELSE
                                        wdelta = w0gauss ( (ef_cond-et(ibnd,ikk))/degauss_cond, ngauss) / degauss_cond
                                        dnsscf(m1,m2,current_spin,nah,ipert) = &
                                        dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def_cond(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                             !  
                             END IF
                          !
                          ELSE
                                          dnsscf(m1,m2,current_spin,nah,ipert) = &
                               dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                          END IF
                          !
                       ENDIF
                       !
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
     ENDDO ! ipert
     !
  ENDDO ! ik 
  !
  CALL mp_sum (dnsscf, inter_pool_comm) 
  !
  ! Filling the m1 m2 lower triangular part 
  !
  DO m1 = 2, ldim
     DO m2 = 1, m1-1
        dnsscf(m1,m2,:,:,:) = dnsscf(m2,m1,:,:,:)
     ENDDO
  ENDDO
  !
  ! In nspin.eq.1 k-point weight wk is normalized to 2 el/band 
  ! in the whole BZ but we are interested in dns of one spin component
  !
  IF (nspin.EQ.1) dnsscf = 0.5d0 * dnsscf
  !
  ! USPP case: add to the dnsscf calculated with the P_c^+dpsi 
  ! the non-scf part coming from the orthogonality contraints. 
  ! The orthogonality correction comes only for an atomic displacement. 
  ! In the case of the electric field calculation, lr_has_dnsorth is false.
  !
  IF (okvan .AND. lr_has_dnsorth) THEN
     dnsscf = dnsscf + lr_dnsorth
  ENDIF
  !
  ! Symmetrize dnsscf
  !
  CALL sym_dns(ldim, npe, dnsscf)
  !
  ! Write symmetrized dnsscf in the pattern basis 
  ! to the standard output
  !
  IF (iverbosity==1) THEN
     WRITE(stdout,*) 'DNSSCF SYMMETRIZED IN THE PATTERN BASIS'
     DO ipert = 1, npe
        WRITE(stdout,'(a,1x,i2)') ' partner # ', ipert
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                 DO m1 = 1, ldim
                    WRITE(stdout,'(10(f15.10,1x))') dnsscf(m1,:,is,nah,ipert)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO     
  ENDIF
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  DEALLOCATE (dpsi)
  !
  CALL stop_clock( 'dnsq_scf' )
  ! 
  RETURN
  ! 
END SUBROUTINE dnsq_scf
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE dnsq_store(npe, imode0)
!----------------------------------------------------------------------------
  !! Store the computed dnsscf in the full matrix dnsscf_all_modes
  !! (i.e for all modes and not only for the npe irreducible representations)
  !
  USE ions_base,     ONLY : nat, ityp
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : is_hubbard, Hubbard_l
  USE ldaU_ph,       ONLY : dnsscf_all_modes
  USE ldaU_lr,       ONLY : dnsscf
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: npe
  !! the number of perturbations
  INTEGER , INTENT(IN) :: imode0
  !! the position of the modes
  !
  INTEGER :: ipert, nah, nt, is, m1, m2
  !
  DO ipert = 1, npe
   DO nah = 1, nat
      nt = ityp(nah)
      IF (is_hubbard(nt)) THEN
         DO is = 1, nspin
            DO m1 = 1, 2*Hubbard_l(nt)+1
               DO m2 = 1, 2*Hubbard_l(nt)+1
                  dnsscf_all_modes(m1,m2,is,nah,imode0+ipert) = &
                                   dnsscf(m1,m2,is,nah,ipert)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE dnsq_store
!----------------------------------------------------------------------------
