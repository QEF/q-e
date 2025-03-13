!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_becsum( iflag )
  !----------------------------------------------------------------------------
  !! Compute "becsum" = \sum_i w_i <psi_i|beta_l><beta_m|\psi_i> term.
  !! Output in module uspp and (PAW only) in rho%bec (symmetrized)
  !! if iflag = 1, weights w_k are re-computed.
  !! This is basically a stripped-down version of sum_band,
  !! without the calculation of the charge density
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : nkb, vkb, becsum, okvan
  USE uspp_param,           ONLY : nhm
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm 
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type_acc, &
                                   deallocate_bec_type_acc, becp
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iflag
  !
  INTEGER :: ik, & ! counter on k points
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  !
  !
  IF ( .NOT. okvan ) RETURN
  !
  CALL start_clock( 'compute_becsum' )
  !
  ! ... calculates weights of Kohn-Sham orbitals
  !
  IF ( iflag == 1) CALL weights( )
  !
  !$acc kernels
  becsum(:,:,:) = 0.D0
  !$acc end kernels
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  CALL allocate_bec_type_acc( nkb, this_bgrp_nbnd, becp,intra_bgrp_comm )
  !
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
     !$acc update device(evc)
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
     !
     ! ... actual calculation is performed (on GPU) inside routine "sum_bec"
     !
     CALL sum_bec( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd )
     !
  ENDDO k_loop
  ! ... Use host copy to do the communications
  !$acc update host(becsum)
  !
  CALL deallocate_bec_type_acc( becp )
  !
  ! ... becsums must be summed over bands (with bgrp parallelization)
  ! ... and over k-points (unsymmetrized).
  !
  CALL mp_sum(becsum, inter_bgrp_comm )
  CALL mp_sum(becsum, inter_pool_comm )
  !
  ! ... Needed for PAW: becsums are stored into rho%bec and symmetrized so that they reflect
  ! ... a real integral in k-space, not only on the irreducible zone. 
  ! ... No need to symmetrize becsums or to align GPU and CPU copies: they are used only here.
  !
  IF ( okpaw )  THEN
     rho%bec(:,:,:) = becsum(:,:,:)
     CALL PAW_symmetrize( rho%bec )
  ENDIF
  !
  CALL stop_clock( 'compute_becsum' )
  !
END SUBROUTINE compute_becsum

!----------------------------------------------------------------------------
SUBROUTINE sum_bec ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
  !----------------------------------------------------------------------------
  !
  !! This routine computes the sum over bands:
  !
  !! \[ \sum_i \langle\psi_i|\beta_l\rangle w_i \langle\beta_m|\psi_i\rangle \]
  !
  !! for point "ik" and, for LSDA, spin "current_spin".  
  !! Calls calbec to compute \(\text{"becp"}=\langle \beta_m|\psi_i \rangle\).  
  !! Output is accumulated (unsymmetrized) into "becsum", module "uspp";
  !! for real-space algorithm, the needed correction is computed and stored into "ebecsum"
  !! Routine used in sum_band (if okvan) and in compute_becsum, called by hinit1 (if okpaw).
  !
  USE kinds,              ONLY : DP
  USE becmod,             ONLY : becp, calbec
  USE control_flags,      ONLY : gamma_only, tqr, offload_type 
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp,               ONLY : nkb, becsum, ebecsum, ofsbeta, vkb
  USE uspp_param,         ONLY : upf, nh, nhm
  USE wvfct,              ONLY : nbnd, wg, et, current_k
  USE klist,              ONLY : ngk, nkstot
  USE noncollin_module,   ONLY : noncolin, npol
  USE wavefunctions,      ONLY : evc
  USE realus,             ONLY : real_space, &
                                 invfft_orbital_gamma, calbec_rs_gamma, &
                                 invfft_orbital_k, calbec_rs_k
  USE us_exx,             ONLY : store_becxx0
  USE mp_bands,           ONLY : nbgrp,inter_bgrp_comm
  USE mp,                 ONLY : mp_sum
  !
  ! Used to avoid unnecessary memcopy
  USE xc_lib,             ONLY : xclib_dft_is
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1(:,:), auxk2(:,:), aux_nc(:,:)
  !$acc declare device_resident (auxk1, auxk2)
  REAL(DP), ALLOCATABLE    :: auxg1(:,:), auxg2(:,:), aux_gk(:,:), aux_egk(:,:)
  !$acc declare device_resident (auxg1, auxg2, aux_gk, aux_egk) 
  INTEGER :: ibnd, kbnd, ibnd_loc, nbnd_loc  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js, nhnt, offset
  ! counters on beta functions, atoms, atom types, spin, and auxiliary vars
  !
  CALL start_clock( 'sum_band:calbec' )
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     !$acc data present(evc) 
     CAll calbec(offload_type, npw, vkb, evc(:,ibnd_start:ibnd_end), becp )
     !$acc end data
  ELSE
     if (gamma_only) then
        do kbnd = 1, this_bgrp_nbnd, 2
           ibnd = ibnd_start + kbnd -1 
           call invfft_orbital_gamma(evc,ibnd,ibnd_end) 
           call calbec_rs_gamma(kbnd,this_bgrp_nbnd,becp%r)
        enddo
        !call mp_sum(becp%r,inter_bgrp_comm)
        !$acc update device(becp%r)
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do kbnd = 1, this_bgrp_nbnd
           ibnd = ibnd_start + kbnd -1
           call invfft_orbital_k(evc,ibnd,ibnd_end)
           call calbec_rs_k(kbnd,this_bgrp_nbnd)
        enddo
        !call mp_sum(becp%k,inter_bgrp_comm)
        !$acc update device(becp%k)
     endif
  ENDIF
  CALL stop_clock( 'sum_band:calbec' )
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  ! FIXME: is this needed? is this a wise thing to do?
  IF(xclib_dft_is('hybrid')) THEN
     if(allocated(becp%r)) then
       !$acc update self(becp%r)
     elseif(allocated(becp%k)) then
       !$acc update self(becp%k)
     elseif(allocated(becp%nc)) then
       !$acc update self(becp%nc)
     endif
     CALL store_becxx0(ik, becp)
  ENDIF
  !
  CALL start_clock( 'sum_band:becsum' )
  !
  !$acc data copyin(wg)
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
           nbnd_loc = becp%nbnd
           ALLOCATE( auxg1( nbnd_loc, nh(np) ) )
           ALLOCATE( auxg2( nbnd_loc, nh(np) ) )
        ELSE
           ALLOCATE( auxk1( ibnd_start:ibnd_end, nh(np)*npol ), &
                     auxk2( ibnd_start:ibnd_end, nh(np)*npol ) )
        END IF
        IF ( noncolin ) THEN
           ALLOCATE ( aux_nc( nh(np)*npol,nh(np)*npol ) ) 
           !$acc enter data create(aux_nc)
        ELSE
           ALLOCATE ( aux_gk( nh(np),nh(np) ) ) 
           if (tqr) ALLOCATE ( aux_egk( nh(np),nh(np) ) ) 
        END IF
        !
        !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
        !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
        !
        nhnt = nh(np)
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              offset = ofsbeta(na)
              IF ( noncolin ) THEN
                 !
                 !$acc parallel loop collapse(2)
                 DO is = 1, npol
                    DO ih = 1, nhnt
                       ikb = offset + ih
                       DO kbnd = 1, this_bgrp_nbnd 
                          ibnd = ibnd_start + kbnd -1 
                          auxk1(ibnd,ih+(is-1)*nhnt)= becp%nc(ikb,is,kbnd)
                          auxk2(ibnd,ih+(is-1)*nhnt)= wg(ibnd,ik) * &
                                                        becp%nc(ikb,is,kbnd)
                       END DO
                    END DO
                 END DO
                 !
                 !$acc host_data use_device(auxk1, auxk2, aux_nc)
                 CALL MYZGEMM ( 'C', 'N', npol*nhnt, npol*nhnt, this_bgrp_nbnd, &
                      (1.0_dp,0.0_dp), auxk1, this_bgrp_nbnd, auxk2, this_bgrp_nbnd, &
                      (0.0_dp,0.0_dp), aux_nc, npol*nhnt )
                 !$acc end host_data
                 !
              ELSE IF ( gamma_only ) THEN
                 !
                 !$acc parallel loop collapse(2)
                 DO ih = 1, nhnt
                    DO ibnd_loc = 1, nbnd_loc
                       ikb = offset + ih
                       ibnd = (ibnd_start -1) + ibnd_loc
                       auxg1(ibnd_loc,ih) = becp%r(ikb,ibnd_loc)
                       auxg2(ibnd_loc,ih) = becp%r(ikb,ibnd_loc) * wg(ibnd,ik)
                    END DO
                 END DO
                 !$acc host_data use_device(auxg1, auxg2, aux_gk)
                 CALL MYDGEMM ( 'T', 'N', nhnt, nhnt, nbnd_loc, &
                      1.0_dp, auxg1, nbnd_loc,    &
                      auxg2, nbnd_loc, 0.0_dp, aux_gk, nhnt )
                 !$acc end host_data
                 !
                 if (tqr) then
                   !$acc parallel loop collapse(2)
                   DO ih = 1, nhnt
                      DO ibnd_loc = 1, nbnd_loc
                      ibnd = (ibnd_start -1) + ibnd_loc
                      auxg2(ibnd_loc,ih) = et(ibnd,ik) * auxg2(ibnd_loc,ih)
                      END DO
                   END DO
                   !$acc host_data use_device(auxg1, auxg2, aux_egk)
                   CALL MYDGEMM ( 'T', 'N', nhnt, nhnt, nbnd_loc, &
                        1.0_dp, auxg1, nbnd_loc,    &
                        auxg2, nbnd_loc, 0.0_dp, aux_egk, nhnt )
                   !$acc end host_data
                 end if
                 !
              ELSE
                 !
                 !$acc parallel loop collapse(2)
                 DO ih = 1, nhnt
                    DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                       ibnd = ibnd_start + kbnd -1 
                       ikb = offset + ih
                       auxk1(ibnd,ih) = becp%k(ikb,kbnd) 
                       auxk2(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
                    END DO
                 END DO
                 !
                 ! only the real part is computed
                 !
                 !$acc host_data use_device(auxk1, auxk2, aux_gk)
                 CALL MYDGEMM ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk, nhnt )
                 !$acc end host_data
                 !
                 if (tqr) then
                   !$acc parallel loop collapse(2)
                   DO ih = 1, nhnt
                      DO ibnd = ibnd_start, ibnd_end
                         auxk2(ibnd,ih) = et(ibnd,ik)*auxk2(ibnd,ih)
                      END DO
                   END DO

                   !$acc host_data use_device(auxk1, auxk2, aux_egk)
                   CALL MYDGEMM ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                        1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                        0.0_dp, aux_egk, nhnt )
                   !$acc end host_data
                 end if

              END IF
              !
              ! copy output from GEMM into desired format
              !
              IF (noncolin .AND. .NOT. upf(np)%has_so) THEN
                 CALL add_becsum_nc(na, np, aux_nc, becsum )
              ELSE IF (noncolin .AND. upf(np)%has_so) THEN
                 CALL add_becsum_so (na, np, aux_nc, becsum )
              ELSE
                 !
                 !$acc parallel loop collapse(2) present(becsum)
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2  ! or use  ijtoh(ih,jh,np) ?  
                       !
                       ! nondiagonal terms summed and collapsed into a
                       ! single index (matrix is symmetric wrt (ih,jh))
                       !
                       IF ( jh == ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk (ih,jh)
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                             ebecsum(ijh,na,current_spin) + aux_egk (ih,jh)
                       ELSE IF ( jh > ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk(ih,jh)*2.0_dp
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                             ebecsum(ijh,na,current_spin) + aux_egk(ih,jh)*2.0_dp
                       END IF
                    END DO
                 END DO
                 !
              END IF
           END IF
           !
        END DO
        !
        IF ( noncolin ) THEN
           !$acc exit data delete(aux_nc)
           DEALLOCATE ( aux_nc )
        ELSE
           DEALLOCATE ( aux_gk  ) 
           if (tqr) DEALLOCATE ( aux_egk  ) 
        END IF
        IF ( gamma_only ) THEN
           DEALLOCATE( auxg2, auxg1 )
        ELSE
           DEALLOCATE( auxk2, auxk1 )
        END IF
        !
     END IF
     !
  END DO
  !$acc end data
  !
  CALL stop_clock( 'sum_band:becsum' )
  !
END SUBROUTINE sum_bec
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_nc ( na, np, becsum_nc, becsum )
!----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the
  !! Pauli matrices, saves it in \(\text{becsum}\) for the calculation of 
  !! augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE uspp,                 ONLY : ijtoh
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, ijh, ipol, jpol, nhnp
  REAL(DP) :: fac
  !
  nhnp = nh(np)

  !$acc parallel loop collapse(2) present(ijtoh, becsum, becsum_nc)
  DO ih = 1, nhnp
     DO jh = 1, nhnp
        IF ( jh >= ih ) THEN
           !ijh = jh + ((ih-1)*(2*nhnp-ih))/2  is this faster? Does it matter?
           ijh=ijtoh(ih,jh,np)
           IF ( ih == jh ) THEN
              fac = 1.0_dp
           ELSE
              fac = 2.0_dp
           END IF
           becsum(ijh,na,1)= becsum(ijh,na,1) + fac * &
                   DBLE( becsum_nc(ih,1,jh,1) + becsum_nc(ih,2,jh,2) )
           IF (domag) THEN
              becsum(ijh,na,2)= becsum(ijh,na,2) + fac *  &
                   DBLE( becsum_nc(ih,1,jh,2) + becsum_nc(ih,2,jh,1) )
              becsum(ijh,na,3)= becsum(ijh,na,3) + fac * DBLE( (0.d0,-1.d0)* &
                  (becsum_nc(ih,1,jh,2) - becsum_nc(ih,2,jh,1)) )
              becsum(ijh,na,4)= becsum(ijh,na,4) + fac * &
                   DBLE( becsum_nc(ih,1,jh,1) - becsum_nc(ih,2,jh,2) )
           END IF
        END IF
     END DO
  END DO
  
END SUBROUTINE add_becsum_nc
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_so( na, np, becsum_nc, becsum )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the Pauli
  !! matrices, rotates it as appropriate for the spin-orbit case, saves it in 
  !! \(\text{becsum}\) for the calculation of augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  USE uspp,                 ONLY : ijtoh, nhtol, nhtoj, indv
  USE upf_spinorb,          ONLY : fcoef
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2, nhnt
  COMPLEX(DP) :: fac
  !
  nhnt = nh(np)
  !
  ! For an obscure reason, if you collapse the first two loops into one
  ! with collapse(2) in the line below, the calculation also collapses!
  !$acc parallel loop present(fcoef, ijtoh, nhtoj, nhtol, indv, becsum, becsum_nc)
  DO ih = 1, nhnt
     DO jh = 1, nhnt
        ijh=ijtoh(ih,jh,np)
        DO kh = 1, nhnt
           IF ( (nhtol(kh,np)==nhtol(ih,np)).AND. &
                (ABS(nhtoj(kh,np)-nhtoj(ih,np))<1.d8).AND. &
                (indv(kh,np)==indv(ih,np)) ) THEN ! same_lj(kh,ih,np)
              DO lh=1,nhnt
                 IF ( (nhtol(lh,np)==nhtol(jh,np)).AND. &
                      (ABS(nhtoj(lh,np)-nhtoj(jh,np))<1.d8).AND. &
                      (indv(lh,np)==indv(jh,np)) ) THEN   !same_lj(lh,jh,np)) THEN
                    DO is1=1,npol
                       DO is2=1,npol
                          fac=becsum_nc(kh,is1,lh,is2)
                          becsum(ijh,na,1)=becsum(ijh,na,1) + DBLE( fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  ) )
                          IF (domag) THEN
                            becsum(ijh,na,2)=becsum(ijh,na,2) + DBLE( fac * &
                                (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) +&
                                 fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  ) )
                            becsum(ijh,na,3)=becsum(ijh,na,3) + DBLE( fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  ))
                            becsum(ijh,na,4)=becsum(ijh,na,4) + DBLE(fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  ) )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO

END SUBROUTINE add_becsum_so
