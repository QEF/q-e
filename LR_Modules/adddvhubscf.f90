!                                         
! Copyright (C) 2001-2015 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE adddvhubscf (ipert, ik)
  !--------------------------------------------------------------------------------------
  ! 
  !! DFPT+U  
  !! This routine calculates the SCF derivative of the Hubbard potential times \(\psi\):
  !! \begin{equation}\notag
  !! \begin{split}
  !!   |\Delta V_{SCF}(k+q,is) \psi(\text{ibnd},k,is)\rangle = 
  !!   - &\sum_{I,m1,m2} \text{Hubbard_U}(I)\cdot\text{dnsscf}(m1,m2,is,I,\text{imode})\cdot \\
  !!     &|S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2|\psi(\text{ibnd},k,is)\rangle
  !! \end{split}
  !! \end{equation}
  !
  !! Addition of the \(\text{J0}\) terms:
  !
  !! $$ + \sum_{I,m1,m2} \text{Hubbard_J0}(I)\cdot \text{dnsscf}(m1,m2,\text{isi},I,\text{imode})\cdot 
  !! |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
  !
  !! Where:  
  !! \(\text{is}\)  = current_spin;  
  !! \(\text{isi}\) = opposite of the current_spin.
  !
  !! Written by A. Floris. Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                            Hubbard_J0, nwfcU
  USE ldaU_lr,       ONLY : effU, swfcatomk, swfcatomkpq, dnsscf
  USE wvfct,         ONLY : npwx, nbnd 
  USE control_lr,    ONLY : lgamma, nbnd_occ
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : dvpsi  
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  USE klist,         ONLY : ngk
  USE buffers,       ONLY : get_buffer
  USE qpoint,        ONLY : xq, ikks, ikqs
  USE units_lr,      ONLY : iuatswfc
  USE noncollin_module,     ONLY : npol
  !  
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k point under consideration
  INTEGER, INTENT(IN) :: ipert  
  !! the index of perturbation
  !  
  ! Local variables
  !
  INTEGER :: npw, npwq, ikk, ikq, op_spin
  INTEGER ::  i, j, k, nt, l, ih, n, ig, ihubst, ihubst1, ihubst2, &
              nah, m, m1, m2, ibnd, ldim
  COMPLEX(DP), ALLOCATABLE :: dvhubscf(:,:), dvqi(:,:), proj1(:, :)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock ( 'adddvhubscf' )
  !
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (dvhubscf(npwx,nbnd))
  ALLOCATE (dvqi(npwx,nbnd))
  !  
  ldim = 2 * Hubbard_lmax + 1
  !
  ! At each iteration dvhubscf is set to zero, 
  ! because it is recomputed through dnsscf 
  !
  dvhubscf = (0.d0, 0.d0) 
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  IF (lsda) THEN
     current_spin = isk (ikk)
     IF (current_spin==1) THEN   
        op_spin = 2
     ELSE
        op_spin = 1
     END IF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Read S * atomic wfc's at k and k+q on unit iuatswfc
  !  
  CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     IF (is_hubbard(nt)) THEN        
        !   
        ! Calculate proj1(ibnd,ihubst) = < psi(inbd,k) | S_{k}\phi_(k,I,m) >
        ldim = 2*Hubbard_l(nt) + 1
        ihubst = offsetU(nah) + 1   ! I m index
        CALL ZGEMM('C', 'N', nbnd, ldim, npw, (1.0d0, 0.0d0), &
                evc, npwx*npol, swfcatomk(1, ihubst), npwx*npol, &
                (0.0d0, 0.0d0), proj1(1, ihubst), nbnd)
        !
     ENDIF
     !
  ENDDO
  ! 
  CALL mp_sum (proj1, intra_bgrp_comm)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     dvqi = (0.d0, 0.d0)
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
                 DO ig = 1, npwq
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) - effU(nt) * &
                                    dnsscf(m1,m2,current_spin,nah,ipert) * &
                                    swfcatomkpq(ig,ihubst1) * CONJG(proj1(ibnd,ihubst2))
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nbnd
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
     !  Hubbard_J0
     !
     dvqi = (0.d0, 0.d0)
     !
     IF (Hubbard_J0(nt).NE.0.d0) THEN
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
                 DO ig = 1, npwq
                    ! Notice that sign change here
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) + Hubbard_J0(nt) * &
                                    dnsscf(m1,m2,op_spin,nah,ipert) * &
                                    swfcatomkpq(ig,ihubst1) * CONJG(proj1(ibnd,ihubst2))
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nbnd
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !  
  ! Add the result to dvpsi
  !
  DO ibnd = 1, nbnd
     DO ig = 1, npwq
        dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + dvhubscf(ig,ibnd)  
     ENDDO
  ENDDO
  !
  DEALLOCATE (proj1)
  DEALLOCATE (dvhubscf)
  DEALLOCATE (dvqi)
  !
  CALL stop_clock ( 'adddvhubscf' )
  !  
  RETURN
  ! 
END SUBROUTINE adddvhubscf

SUBROUTINE revert_mag_u (ns_nc)
   !
   ! This routine reverts the sign of the Hubbard magnetization,
   ! to be used in the noncollinear magnetic case.
   ! Written by L. Binci (2023)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin
   USE ldaU,                 ONLY : lda_plus_u_kind, Hubbard_l, is_hubbard, &
                                    Hubbard_lmax
   USE noncollin_module,     ONLY : npol
   !
   IMPLICIT NONE
   INTEGER :: na, nt, m1, m2, is, is1, is2, ldim
   COMPLEX(DP) :: ns, ms(3)
   COMPLEX(DP), INTENT(INOUT) :: ns_nc(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
   !
   DO na = 1, nat
      nt = ityp(na)
      IF (is_hubbard(nt)) THEN
          ldim = 2 * Hubbard_l(nt) + 1
          DO m1 = 1, ldim
              DO m2 = 1, ldim
                   !
                   ! charge
                   ns = (ns_nc(m1,m2,1,na) + ns_nc(m1,m2,4,na))
                   ! magnetization
                   ms(1) =  (ns_nc(m1,m2,2,na)  + ns_nc(m1,m2,3,na))
                   ms(2) =  -(0.d0,1.0d0)*(ns_nc(m1,m2,2,na)  - ns_nc(m1,m2,3,na))
                   ms(3) =  (ns_nc(m1,m2,1,na)  - ns_nc(m1,m2,4,na))
                   !
                   ms(:) = -ms(:)
                   !
                   ns_nc(m1,m2,1,na) = 0.5d0*( ns + ms(3))
                   ns_nc(m1,m2,2,na) = 0.5d0*( ms(1) + (0.d0,1.0d0)*ms(2) )
                   ns_nc(m1,m2,3,na) = 0.5d0*( ms(1) - (0.d0,1.0d0)*ms(2) )
                   ns_nc(m1,m2,4,na) = 0.5d0*( ns - ms(3) )
              ENDDO
          ENDDO
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE revert_mag_u

 SUBROUTINE calc_vh_u( ns, v_hub)
   !
   !! Recomputes the noncollinear Hubbard potential.
   !! Used for a (spin-)Hubbard matrix with reversed magnetization
   !! Similar to v_hubbard_nc within PW/src/v_of_rho.f90
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U
   USE lsda_mod,             ONLY : nspin
   USE control_flags,        ONLY : iverbosity, dfpt_hub
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
   !! Occupation matrix
   COMPLEX(DP), INTENT(INOUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
   !! Hubbard potential
   INTEGER  :: is, is1, na, nt, m1, m2
   ! ... local variables
   !
   v_hub(:,:,:,:) = 0.d0
   !
   DO na = 1, nat
      nt = ityp (na)
      IF (Hubbard_U(nt) /= 0.d0) THEN
         !
         ! is=1 and is=4 are diagonal components
         ! is=2 and is=3 are off-diagonal components
         DO is = 1, nspin
            IF (is == 2) THEN
             is1 = 3
            ELSEIF (is == 3) THEN
             is1 = 2
            ELSE
             is1 = is ! is=1 or is=4
            ENDIF
            !
            IF (is1 == is) THEN
               !
               ! Non spin-flip contribution (is=1 and is=4)
               ! (diagonal [spin indexes] occupancy matrices)
               DO m1 = 1, 2*Hubbard_l(nt) + 1
                  v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + 0.5D0*Hubbard_U(nt)
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1
                     v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - Hubbard_U(nt) * ns(m2,m1,is,na)
                  ENDDO
               ENDDO
            ELSE
               !
               ! Spin-flip contribution (is=2 and is=3)
               ! (NON-diagonal [spin indexes] occupancy matrices)
               DO m1 = 1, 2*Hubbard_l(nt) + 1
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1
                     v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - Hubbard_U(nt)*ns(m2,m1,is1,na)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO !is
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE calc_vh_u
 !--------------------------------------------------------------------------

 SUBROUTINE revert_mag_uv (nsg_nc)
   !
   ! This routine reverts the sign of the Hubbard magnetization,
   ! to be used in the noncollinear magnetic case.
   ! Written by L. Binci (2023)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin
   USE ldaU,                 ONLY : lda_plus_u_kind, Hubbard_l, is_hubbard, &
                                    Hubbard_lmax, ldmx_tot, max_num_neighbors,ldim_u,&
                                    neighood, at_sc
   USE noncollin_module,     ONLY : npol
   USE control_lr,           ONLY : nbnd_occ
   !
   IMPLICIT NONE
   INTEGER :: na, nt, m1, m2, is, is1, is2, ldim, ibnd, nt1, nt2, na1, na2, viz,&
              ldim1, ldim2, eq_na2, isp1, isp2
   INTEGER, EXTERNAL :: find_viz
   COMPLEX(DP) :: ns, ms(3)
   COMPLEX(DP), INTENT(INOUT) :: nsg_nc(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   !
   DO na1 = 1, nat
      nt1 = ityp (na1)
      IF ( ldim_u(nt1).GT.0 ) THEN
         ldim1 = ldim_u(nt1)
         DO viz = 1, neighood(na1)%num_neigh
            na2 = neighood(na1)%neigh(viz)
            eq_na2 = at_sc(na2)%at
            nt2 = ityp (eq_na2)
            ldim2 = ldim_u(nt2)
            !IF (na1.GT.na2) THEN
            !   DO m1 = 1, ldim1
            !      DO m2 = 1, ldim2
            !         DO is1 = 1, npol
            !            DO is2 = 1, npol
            !               isp1 = is2 + npol*(is1-1)
            !               isp2 = is1 + npol*(is2-1)
            !               nsg_nc(m2,m1,viz,na1,isp1) = &
            !                  CONJG(nsg_nc(m1,m2,find_viz(na2,na1),na2,isp2))
            !            ENDDO
            !         ENDDO
            !      ENDDO
            !   ENDDO
            !ELSE
               DO m1 = 1, ldim1
                  DO m2 = 1, ldim2
                     ! charge
                     ns = (nsg_nc(m2,m1,viz,na1,1) + nsg_nc(m2,m1,viz,na1,4))
                     ! magnetization
                     ms(1) = (nsg_nc(m2,m1,viz,na1,2) + nsg_nc(m2,m1,viz,na1,3))
                     ms(2) = -(0.d0,1.0d0)*(nsg_nc(m2,m1,viz,na1,2) - nsg_nc(m2,m1,viz,na1,3))
                     ms(3) = (nsg_nc(m2,m1,viz,na1,1) - nsg_nc(m2,m1,viz,na1,4))
                     !
                     ms(:) = -ms(:)
                     !
                     nsg_nc(m2,m1,viz,na1,1) = 0.5d0*( ns + ms(3))
                     nsg_nc(m2,m1,viz,na1,2) = 0.5d0*( ms(1) + (0.d0,1.0d0)*ms(2) )
                     nsg_nc(m2,m1,viz,na1,3) = 0.5d0*( ms(1) - (0.d0,1.0d0)*ms(2) )
                     nsg_nc(m2,m1,viz,na1,4) = 0.5d0*( ns - ms(3) )
                  ENDDO
               ENDDO
            !ENDIF
         ENDDO
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE revert_mag_uv

 SUBROUTINE calc_vh_uv( nsg, v_hub)
   !
   !! Recomputes the noncollinear Hubbard potential.
   !! Used for a (spin-)Hubbard matrix with reversed magnetization
   !! Similar to v_hubbard_nc within PW/src/v_of_rho.f90
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                 ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                 Hubbard_V, Hubbard_alpha_back, is_hubbard, is_hubbard_back
   USE lsda_mod,          ONLY : nspin
   USE control_flags,     ONLY : iverbosity, dfpt_hub
   USE io_global,         ONLY : stdout
   USE noncollin_module,  ONLY : npol
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   !
   ! Local variables
   !
   INTEGER :: is, is1, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2
   INTEGER, EXTERNAL :: find_viz
   !
   v_hub(:,:,:,:,:) = (0.d0, 0.d0)
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      IF ( is_hubbard(nt1) ) THEN
         DO is = 1, nspin
            IF (is == 2) THEN
               is1 = 3
            ELSEIF (is == 3) THEN
               is1 = 2
            ELSE
               is1 = is
            ENDIF
            DO viz = 1, neighood(na1)%num_neigh
               na2 = neighood(na1)%neigh(viz)
               equiv_na2 = at_sc(na2)%at
               nt2 = ityp(equiv_na2)
               !
               IF (is_hubbard(nt2) .AND. &
                   (Hubbard_V(na1,na2,1).NE.0.d0) ) THEN
                   !
                   ! Here no need to use is1: complex conjugation is enough
                   ! For both standard and background states of a center atom
                   DO m1 = 1, ldim_u(nt1)
                      ! For both standard and background states of the neighbor atom
                      DO m2 = 1, ldim_u(nt2)
                         v_hub(m2,m1,viz,na1,is) = - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                      ENDDO
                   ENDDO
                   !
                   IF ( na1.EQ.na2 .AND. is1.EQ.is) THEN
                      !
                      na = find_viz(na1,na1)
                      ! This is the diagonal term (like in the DFT+U only case)
                      DO m1 = 1, ldim_u(nt1)
                         v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                        + Hubbard_V(na1,na1,1) * 0.5d0
                      ENDDO
                      !
                   ENDIF
                   !
               ENDIF
               !
            ENDDO ! viz
            !
         ENDDO ! is
         !
      ENDIF
      !
   ENDDO ! na1
   !
   !
   RETURN
   !
 END SUBROUTINE calc_vh_uv
 !--------------------------------------------------------------------------
