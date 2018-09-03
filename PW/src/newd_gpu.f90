!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dfunct_gpum

CONTAINS
!---------------------------------------

SUBROUTINE newq_gpu(vr,deeq_d,skip_vltot)
  !
  !   This routine computes the integral of the perturbed potential with
  !   the Q function
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasZgemm Zgemm
#define cublasDGEMM Dgemm
#define cudaDGER    Dger
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, gstart, mill, &
                                   eigts1, eigts2, eigts3
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vltot
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE wavefunctions_gpum,   ONLY : psic_d
  USE spin_orb,             ONLY : lspinorb, domag
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  USE gvect_gpum,           ONLY : g_d, gg_d, mill_d, &
                                   eigts1_d, eigts2_d, eigts3_d
  USE gbuffers,             ONLY : buffer=>dev_buf
  !
  IMPLICIT NONE
  !
  !
  ! Input: potential , output: contribution to integral
  REAL(kind=dp), intent(in)  :: vr(dfftp%nnr,nspin)
  REAL(kind=dp), intent(out) :: deeq_d( nhm, nhm, nat, nspin )
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_d
#endif
  LOGICAL, intent(in) :: skip_vltot !If .false. vltot is added to vr when necessary
  ! INTERNAL
  INTEGER :: ngm_s, ngm_e, ngm_l
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, nt, ih, jh, na, is, ijh, nij, nb, nab, nhnt, ierr
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  COMPLEX(DP), ALLOCATABLE :: vaux_d(:,:), aux_d(:,:), qgm_d(:,:)
    ! work space
  REAL(DP), ALLOCATABLE :: ylmk0_d(:,:), qmod_d(:), deeaux_d(:,:)
    ! spherical harmonics, modulus of G
  REAL(DP) :: fact
  !
  INTEGER, POINTER :: dfftp_nl_d(:)
  ! workaround for cuf kernel limitations
  !
#if defined(__CUDA)
  attributes(DEVICE) :: vaux_d, aux_d, qgm_d, ylmk0_d, qmod_d, deeaux_d, dfftp_nl_d
#endif
  ! variable to map index of atoms of the same type
  INTEGER, ALLOCATABLE :: na_to_nab_h(:)
  INTEGER, POINTER     :: na_to_nab_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: na_to_nab_d
#endif
  !
  ALLOCATE(na_to_nab_h(nat))
  CALL buffer%lock_buffer(na_to_nab_d, nat, ierr)
  !
  IF ( gamma_only ) THEN
     fact = 2.0_dp
  ELSE
     fact = 1.0_dp
  END IF
  !
  deeq_d(:,:,:,:) = 0.D0
  !
  ! With k-point parallelization, distribute G-vectors across processors
  ! ngm_s = index of first G-vector for this processor
  ! ngm_e = index of last  G-vector for this processor
  ! ngm_l = local number of G-vectors 
  !
  CALL divide (inter_pool_comm, ngm, ngm_s, ngm_e)
  ngm_l = ngm_e-ngm_s+1
  ! for the extraordinary unlikely case of more processors than G-vectors
  !
  IF ( ngm_l > 0 ) THEN
     ALLOCATE( vaux_d(ngm_l,nspin_mag), qmod_d(ngm_l), ylmk0_d( ngm_l, lmaxq*lmaxq ) )
     !
     CALL ylmr2_gpu (lmaxq * lmaxq, ngm_l, g_d(1,ngm_s), gg_d(ngm_s), ylmk0_d)
!$cuf kernel do
     DO ig = 1, ngm_l
        qmod_d (ig) = sqrt (gg_d (ngm_s+ig-1) )
     ENDDO
  END IF
  !
  ! ... fourier transform of the total effective potential
  !
  dfftp_nl_d => dfftp%nl_d
  DO is = 1, nspin_mag
     !
     IF ( (nspin_mag == 4 .AND. is /= 1) .or. skip_vltot ) THEN 

        psic_d(1:dfftp%nnr) = vr(1:dfftp%nnr,is)

     ELSE
!$omp parallel do default(shared) private(ig)
        do ig=1,dfftp%nnr
           psic(ig) = vltot(ig) + vr(ig,is)
        end do
!$omp end parallel do
        psic_d(1:dfftp%nnr) = psic(1:dfftp%nnr)
     END IF
     CALL fwfft ('Rho', psic_d, dfftp)
!$cuf kernel do
     do ig=1,ngm_l
        vaux_d(ig, is) = psic_d(dfftp_nl_d(ngm_s+ig-1))
     end do
     !
  END DO

  DO nt = 1, ntyp
     !
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nhnt = nh(nt)
        nij = nh(nt)*(nh(nt)+1)/2
        ALLOCATE ( qgm_d(ngm_l,nij) )
        !
        ! ... Compute and store Q(G) for this atomic species 
        ! ... (without structure factor)
        !
        ijh = 0
        DO ih = 1, nhnt
           DO jh = ih, nhnt
              ijh = ijh + 1
              CALL qvan2_gpu ( ngm_l, ih, jh, nt, qmod_d, qgm_d(1,ijh), ylmk0_d )
           END DO
        END DO
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
           IF ( ityp(na) == nt ) THEN
              na_to_nab_h(na)           = nab
           ELSE
              na_to_nab_h(na)           = -1
           END IF
        END DO
        na_to_nab_d(1:nat) = na_to_nab_h(1:nat)
        
        ALLOCATE ( aux_d (ngm_l, nab ), deeaux_d(nij, nab) )
        !
        ! ... Compute and store V(G) times the structure factor e^(-iG*tau)
        !
        DO is = 1, nspin_mag
!$cuf kernel do(2)
           DO na = 1, nat
              DO ig=1,ngm_l
                 nb = na_to_nab_d(na)
                 IF (nb > 0) &
                    aux_d(ig, nb) = vaux_d(ig,is) * CONJG ( &
                      eigts1_d(mill_d(1,ngm_s+ig-1),na) * &
                      eigts2_d(mill_d(2,ngm_s+ig-1),na) * &
                      eigts3_d(mill_d(3,ngm_s+ig-1),na) )
              END DO
           END DO
           !
           ! ... here we compute the integral Q*V for all atoms of this kind
           !
           CALL cublasDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm_d, 2*ngm_l, aux_d, &
                    2*ngm_l, 0.0_dp, deeaux_d, nij )
           IF ( gamma_only .AND. gstart == 2 ) &
                CALL cudaDGER(nij, nab,-1.0_dp, qgm_d, 2*ngm_l,aux_d,2*ngm_l,deeaux_d,nij)
           !
           nhnt = nh(nt)
!$cuf kernel do(3)
           DO na = 1, nat
              DO ih = 1, nhnt
                 DO jh = 1, nhnt
                    nb = na_to_nab_d(na)
                    IF (nb > 0) THEN
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                       IF (jh >= ih) deeq_d(ih,jh,na,is) = omega * deeaux_d(ijh,nb)
                       IF (jh > ih) deeq_d(jh,ih,na,is) = deeq_d(ih,jh,na,is)
                    END IF
                 END DO
              END DO
           END DO
           !
        END DO
        !
        DEALLOCATE ( deeaux_d, aux_d, qgm_d )
        !
     END IF
     !
  END DO
  !
  DEALLOCATE( qmod_d, ylmk0_d, vaux_d )
  DEALLOCATE(na_to_nab_h); CALL buffer%release_buffer(na_to_nab_d, ierr)
  !
  ! REPLACE THIS WITH THE NEW allgather with type or use CPU variable! OPTIMIZE HERE
  CALL mp_sum( deeq_d( :, :, :, 1:nspin_mag ), inter_pool_comm )
  CALL mp_sum( deeq_d( :, :, :, 1:nspin_mag ), intra_bgrp_comm )
  !
END SUBROUTINE newq_gpu
  !
!----------------------------------------------------------------------------
SUBROUTINE newd_gpu( ) 
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the integral of the effective potential with
  ! ... the Q function and adds it to the bare ionic D term which is used
  ! ... to compute the non-local term in the US scheme.
  !
#if defined(__CUDA)
  use cudafor
  use cublas
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE spin_orb,             ONLY : lspinorb, domag
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : nhtol, nhtolm
  USE scf,                  ONLY : v
  USE realus,        ONLY : newq_r
  USE control_flags, ONLY : tqr
  USE ldaU,          ONLY : lda_plus_U, U_projection
  !
  USE uspp_gpum,     ONLY : deeq_d, deeq_nc_d, &
                            using_deeq, using_deeq_nc, &
                            using_deeq_d, using_deeq_nc_d, &
                            dvan_d, dvan_so_d
                            
  USE gbuffers,      ONLY : buffer=>dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht, nb, mb, ierr
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  REAL(kind=dp), allocatable :: deeq_h( :,:,:,: )
  INTEGER, POINTER :: ityp_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: ityp_d
#endif
  !
  IF ( .NOT. okvan ) THEN
     ! Sync
     IF ( lspinorb .or. noncolin ) THEN
       CALL using_deeq_nc_d(1)
     ELSE
       CALL using_deeq_d(1)    ! Can these be changed from INOUT to OUT?
     END IF
     !
     ! ... no ultrasoft potentials: use bare coefficients for projectors
     !
     CALL buffer%lock_buffer(ityp_d, nat, ierr)
     ityp_d(1:nat)=ityp(1:nat)
     DO nt = 1, ntyp
        !
        nht = nh(nt)
        !
        IF ( lspinorb ) THEN
           !
           !$cuf kernel do(4)
           DO is =  1, nspin
              DO na = 1, nat
                 DO jh = 1, nht
                    DO ih = 1, nht
                       IF ( ityp_d(na) == nt ) deeq_nc_d(ih,jh,na,is) = dvan_so_d(ih,jh,is,nt)
                    END DO
                 END DO
              END DO
           END DO
           !
        ELSE IF ( noncolin ) THEN
           !
           !$cuf kernel do(3)
           DO na = 1, nat
              DO jh = 1, nht
                 DO ih = 1, nht
                    IF ( ityp_d(na) == nt ) THEN
                       deeq_nc_d(ih,jh,na,1) = dvan_d(ih,jh,nt)
                       deeq_nc_d(ih,jh,na,2) = ( 0.D0, 0.D0 )
                       deeq_nc_d(ih,jh,na,3) = ( 0.D0, 0.D0 )
                       deeq_nc_d(ih,jh,na,4) = dvan_d(ih,jh,nt)
                    END IF
                 END DO
              END DO
           END DO
           !
        ELSE
           !
           !$cuf kernel do(4)
           DO is = 1, nspin
              DO na = 1, nat
                 DO jh = 1, nht
                    DO ih = 1, nht
                       !
                       IF ( ityp_d(na) == nt ) deeq_d(ih,jh,na,is) = dvan_d(ih,jh,nt)
                       !
                    END DO
                 END DO
              END DO
           END DO
           !
        END IF
        !
     END DO
     !
     ! ... early return
     !
     CALL buffer%release_buffer(ityp_d, ierr)
     RETURN
     !
  END IF
  !
  CALL start_clock( 'newd' )
  allocate(deeq_h( nhm, nhm, nat, nspin ))
  !
  ! move atom type info to GPU
  CALL buffer%lock_buffer(ityp_d, nat, ierr)
  ityp_d(1:nat)=ityp(1:nat)
  !
  ! Sync
  CALL using_deeq_d(2)
  IF ( lspinorb .or. noncolin ) CALL using_deeq_nc_d(1) ! lspinorb implies noncolin 
  !
  IF (tqr) THEN
     CALL using_deeq(2)
     CALL newq_r(v%of_r,deeq,.false.)
     CALL using_deeq_d(1)
  ELSE
     CALL newq_gpu(v%of_r,deeq_d,.false.)
  END IF
  !
  IF (noncolin) call add_paw_to_deeq_gpu(deeq_d)
  !
  types : &
  DO nt = 1, ntyp
     !
     if_noncolin:&
     IF ( noncolin ) THEN
        !
        IF (upf(nt)%has_so) THEN
           !
           CALL newd_so_gpu(nt)
           !
        ELSE
           !
           CALL newd_nc_gpu(nt)
           !
        END IF
        !
     ELSE if_noncolin
        !
        nht = nh(nt)
        !$cuf kernel do(4)
        DO is = 1, nspin
           DO na = 1, nat
              DO ih = 1, nht
                 DO jh = 1, nht
                    IF ( ityp_d(na) == nt ) THEN
                       deeq_d(ih,jh,na,is) = deeq_d(ih,jh,na,is) + dvan_d(ih,jh,nt)
                    END IF
                 END DO
              END DO
           END DO
        END DO
        !
     END IF if_noncolin
     !
  END DO types
  !
  IF (.NOT.noncolin) CALL add_paw_to_deeq_gpu(deeq_d)
  !
  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL add_vhub_to_deeq_gpu(deeq_d)
  !
  CALL buffer%release_buffer(ityp_d, ierr)
  CALL stop_clock( 'newd' )
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_so_gpu(nt)
      !------------------------------------------------------------------------
      !
      USE spin_orb_gpum, ONLY : fcoef_d, using_fcoef_d
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER :: nt

      INTEGER :: ijs, is1, is2, kh, lh, nhnt, ih, jh, na
      !
      CALL using_fcoef_d(0)
      !
      nhnt = nh(nt)
      ijs = 0
      !
      DO is1 = 1, 2
         !
         DO is2 =1, 2
            !
            ijs = ijs + 1
            !
            IF (domag) THEN
!$cuf kernel do(3)
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp_d(na) == nt ) THEN
                           !
                           deeq_nc_d(ih,jh,na,ijs) = dvan_so_d(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc_d(ih,jh,na,ijs) = deeq_nc_d(ih,jh,na,ijs) +   &
                                      deeq_d (kh,lh,na,1)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt)) + &
                                   deeq_d (kh,lh,na,2)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,2,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,1,is2,nt)) + &
                                   (0.D0,-1.D0)*deeq_d (kh,lh,na,3)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,2,is2,nt)  - &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,1,is2,nt)) + &
                                   deeq_d (kh,lh,na,4)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  - &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt))   
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                  END DO
                  !
               END DO
               !
            ELSE
               !
!$cuf kernel do(3) <<<*,*>>>
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp_d(na) == nt ) THEN
                           !
                           deeq_nc_d(ih,jh,na,ijs) = dvan_so_d(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc_d(ih,jh,na,ijs) = deeq_nc_d(ih,jh,na,ijs) +   &
                                      deeq_d (kh,lh,na,1)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt) ) 
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
               !
            END IF
            !
         END DO
         !
      END DO
      !
    RETURN
      !
    END SUBROUTINE newd_so_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_nc_gpu(nt)
      !------------------------------------------------------------------------
      !
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nt
      INTEGER :: nhnt, na, ih, jh
      !
      nhnt = nh(nt)
      !
!$cuf kernel do(3)
      DO na = 1, nat
         !
         DO ih = 1, nhnt
            !
            DO jh = 1, nhnt
               !
               IF ( ityp_d(na) == nt ) THEN
                  !
                  IF (lspinorb) THEN
                     deeq_nc_d(ih,jh,na,1) = dvan_so_d(ih,jh,1,nt) + &
                                           deeq_d(ih,jh,na,1) + deeq_d(ih,jh,na,4)
                     !                      
                     deeq_nc_d(ih,jh,na,4) = dvan_so_d(ih,jh,4,nt) + &
                                           deeq_d(ih,jh,na,1) - deeq_d(ih,jh,na,4)
                     !
                  ELSE
                     deeq_nc_d(ih,jh,na,1) = dvan_d(ih,jh,nt) + &
                                           deeq_d(ih,jh,na,1) + deeq_d(ih,jh,na,4)
                     !                      
                     deeq_nc_d(ih,jh,na,4) = dvan_d(ih,jh,nt) + &
                                           deeq_d(ih,jh,na,1) - deeq_d(ih,jh,na,4)
                     !
                  END IF
                  deeq_nc_d(ih,jh,na,2) = deeq_d(ih,jh,na,2) - &
                                        ( 0.D0, 1.D0 ) * deeq_d(ih,jh,na,3)
                  !                      
                  deeq_nc_d(ih,jh,na,3) = deeq_d(ih,jh,na,2) + &
                                        ( 0.D0, 1.D0 ) * deeq_d(ih,jh,na,3)
                  !
               END IF
               !
            END DO
            !
         END DO
         !
      END DO
      !
    RETURN
    END SUBROUTINE newd_nc_gpu
    !
END SUBROUTINE newd_gpu

END MODULE dfunct_gpum
