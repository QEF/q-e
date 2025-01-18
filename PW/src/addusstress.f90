!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusstress( sigmanlc )
  !----------------------------------------------------------------------
  !! Driver routine to compute the part of the crystal stress which is due
  !! to the dependence of the Q function on the atomic position.
  !
  USE kinds,          ONLY : dp
  USE control_flags,  ONLY : tqr
  USE realus,         ONLY : addusstress_r
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: sigmanlc(3, 3)
  !! the nonlocal stress
  !
  ! ... local variables
  !
  REAL(DP) :: sigma_r(3,3), sigma_g(3,3)
  INTEGER  :: na, ijh, ipol, jpol
  !
  IF ( tqr ) THEN
     sigma_r(:,:) = 0.d0
     CALL addusstress_r( sigma_r )
     !WRITE (6,'(A)') 'addusstress_r'
     !WRITE (6,'(3f13.8)') sigma_r
     sigmanlc = sigmanlc + sigma_r
     sigma_g(:,:) = 0.d0
     CALL addusstress_g( sigma_g )
!     sigmanlc = sigmanlc + sigma_g
     !WRITE (6,'(A)') 'addusstress_g'
     !WRITE (6,'(3f13.8)') sigma_g
  ELSE
     sigma_g(:,:) = 0.d0
     CALL addusstress_g( sigma_g )
     sigmanlc = sigmanlc + sigma_g
     !WRITE (6,'(A)') 'addusstress_g'
     !WRITE (6,'(3f13.8)') sigma_g
  END IF
  !
END SUBROUTINE addusstress
!
!----------------------------------------------------------------------
SUBROUTINE addusstress_g( sigmanlc )
  !----------------------------------------------------------------------
  !! This routine computes the part of the crystal stress which is due
  !! to the dependence of the Q function on the atomic position.  
  !! It adds contribution to input \(\text{sigmanlc}\), it does not sum 
  !! contributions from various processors (sum is performed by calling
  !! routine).
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp=>nsp, ityp
  USE cell_base,      ONLY : omega, tpiba
  USE fft_base,       ONLY : dfftp
  USE fft_rho,        ONLY : rho_r2g
  USE gvect,          ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE lsda_mod,       ONLY : nspin
  USE scf,            ONLY : v, vltot
  USE uspp,           ONLY : becsum, okvan
  USE uspp_param,     ONLY : upf, lmaxq, nh, nhm
  USE control_flags,  ONLY : gamma_only
  USE mp_pools,       ONLY : inter_pool_comm
  USE mp,             ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! the nonlocal stress
  !
  ! ... local variables
  !
  INTEGER :: ngm_s, ngm_e, ngm_l
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, igm, nt, ih, jh, ijh, ipol, jpol, is, na, nij
  ! counters
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:), aux2(:,:), vg(:,:), qgm(:,:)
  ! work space (complex)
  COMPLEX(DP) :: cfac 
  REAL(DP) :: fac(3,nspin), sus(3,3)
  ! auxiliary variables
  REAL(DP) , ALLOCATABLE :: qmod(:), ylmk0(:,:), dylmk0(:,:), tbecsum(:,:)
  ! work space (real)
  !
  sus(:,:) = 0.d0
  !
  ! ... Fourier transform of the total effective potential
  !
  ALLOCATE( vg(ngm,nspin) )
  !$acc data create( vg )
  !
  DO is = 1, nspin
    IF ( nspin == 4 .and. is /= 1 ) THEN
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is) )
    ELSE
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is), v=vltot(:) )
    ENDIF
  ENDDO
  !
  ! ... With k-point parallelization, distribute G-vectors across processors:
  !     ngm_s = index of first G-vector for this processor
  !     ngm_e = index of last  G-vector for this processor
  !     ngm_l = local number of G-vectors 
  !
  CALL divide( inter_pool_comm, ngm, ngm_s, ngm_e )
  ngm_l = ngm_e-ngm_s+1
  ! for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
  ALLOCATE( aux1(ngm_l,3), aux2(ngm_l,nspin), qmod(ngm_l) )
  ALLOCATE( ylmk0(ngm_l,lmaxq*lmaxq), dylmk0(ngm_l,lmaxq*lmaxq) )
  !$acc data create(aux1,aux2,qmod)
  !$acc data create(ylmk0,dylmk0)
  !
  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  !
  !$acc parallel loop
  DO ig = 1, ngm_l
     qmod(ig) = SQRT( gg(ngm_s+ig-1) ) * tpiba
  ENDDO
  !
  ! ... here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  !     (no contribution from G=0)
  !
  DO ipol = 1, 3
     !
     CALL dylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), dylmk0, ipol )
     !
     DO nt = 1, ntyp
        !
        IF ( upf(nt)%tvanp ) THEN
           nij = nh(nt)*(nh(nt)+1)/2
           ALLOCATE( qgm(ngm_l,nij), tbecsum(nij,nspin) )
           !$acc data create(qgm,tbecsum)
           ijh = 0
           DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                 ijh = ijh + 1
                 CALL dqvan2( ih, jh, nt, ipol, ngm_l, g(1,ngm_s), tpiba, &
                              qmod, ylmk0, dylmk0, qgm(1,ijh) )
              ENDDO
           ENDDO
           !
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 !
                 !$acc kernels
                 tbecsum(:,:) = becsum(1:nij,na,1:nspin)
                 !$acc end kernels
                 !
                 !$acc host_data use_device(qgm,tbecsum,aux2)
                 CALL MYDGEMM( 'N', 'N', 2*ngm_l, nspin, nij, 1.0_dp, &
                               qgm, 2*ngm_l, tbecsum, nij, 0.0_dp, aux2, 2*ngm_l )
                 !$acc end host_data
                 !
#if defined(_OPENACC)
!$acc parallel loop collapse(2)
#else
!$omp parallel do default(shared) private(is, ig)
#endif
                 DO is = 1, nspin
                    DO ig = 1, ngm_l
                       aux2(ig,is) = aux2(ig,is) * CONJG(vg (ngm_s+ig-1, is))
                    ENDDO
                 ENDDO
#if defined(_OPENACC)
!$acc parallel loop present(eigts1,eigts2,eigts3,mill,g)
#else
!$omp end parallel do
!$omp parallel do default(shared) private(ig,igm,cfac)
#endif
                 DO ig = 1, ngm_l
                    igm = ngm_s+ig-1
                    cfac = CONJG( eigts1(mill(1,igm),na) * &
                                  eigts2(mill(2,igm),na) * &
                                  eigts3(mill(3,igm),na) ) * tpiba
                    aux1(ig,1) = cfac * g(1,igm)
                    aux1(ig,2) = cfac * g(2,igm)
                    aux1(ig,3) = cfac * g(3,igm)
                 ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
                 !
                 !$acc data copyout( fac )
                 !$acc host_data use_device(aux1,aux2,fac)
                 CALL MYDGEMM( 'T','N', 3, nspin, 2*ngm_l, 1.0_dp, aux1, 2*ngm_l, &
                               aux2, 2*ngm_l, 0.0_dp, fac, 3 )
                 !$acc end host_data
                 !$acc end data
                 !
                 DO is = 1, nspin
                    DO jpol = 1, 3
                       sus(ipol,jpol) = sus(ipol,jpol) - omega * fac(jpol,is)
                    ENDDO
                 ENDDO
                 !
              ENDIF
           ENDDO
           !$acc end data
           DEALLOCATE( tbecsum, qgm )
        ENDIF
     ENDDO
     !
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( ylmk0, dylmk0 )
  DEALLOCATE( aux1, aux2, qmod )
  !
10 CONTINUE
  !
  CALL mp_sum( sus, inter_pool_comm )
  !
  IF (gamma_only) THEN
     sigmanlc(:,:) = sigmanlc(:,:) + 2.0_dp*sus(:,:)
  ELSE
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  ENDIF
  !
  !$acc end data
  DEALLOCATE( vg )
  !
  RETURN
  !
END SUBROUTINE addusstress_g
