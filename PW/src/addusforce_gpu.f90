!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusforce_gpu( forcenl )
  !----------------------------------------------------------------------
  !! Wrapper to \(\texttt{addusforce_g}\) or \(\texttt{addusforce_r}\).
  !
  USE kinds,         ONLY : dp
  USE ions_base,     ONLY : nat
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : addusforce_r
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  IF ( tqr ) THEN
     CALL addusforce_r( forcenl )
  ELSE
     CALL addusforce_g_gpu( forcenl )
  ENDIF
  !
END SUBROUTINE addusforce_gpu 
!
!----------------------------------------------------------------------
SUBROUTINE addusforce_g_gpu( forcenl )
  !----------------------------------------------------------------------
  !! This routine computes the contribution to atomic forces due
  !! to the dependence of the Q function on the atomic position.
  !! \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
  !!    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
  !! where:
  !! \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
  !!    w_i\langle \beta_m|\psi_i\rangle \]
  !! On output: the contribution is added to \(\text{forcenl}\).
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE cell_base,          ONLY : omega, tpiba
  USE fft_base,           ONLY : dfftp
  USE gvect,              ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,                ONLY : v, vltot
  USE uspp,               ONLY : becsum, okvan
  USE uspp_param,         ONLY : upf, lmaxq, nh, nhm
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp_pools,           ONLY : inter_pool_comm
  USE mp,                 ONLY : mp_sum
  USE control_flags,      ONLY : gamma_only
  USE fft_interfaces,     ONLY : fwfft
  !
  USE gvect_gpum,         ONLY : gg_d, g_d, eigts1_d, eigts2_d, eigts3_d, mill_d
  !
  USE uspp_gpum,     ONLY : using_becsum_d, becsum_d
  USE device_fbuff_m,      ONLY : dev_buf
#if defined(__CUDA) 
  USE cudafor 
  USE cublas
#else
#define cublasDGEMM Dgemm
#endif
 
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl (3, nat)
  !! the non-local contribution to the force
  !
  ! ... local variables
  !
  INTEGER  :: ngm_s, ngm_e, ngm_l
  ! starting/ending indices, local number of G-vectors
  INTEGER  :: ig, nt, ih, jh, ijh, nij, ipol, is, na, nb, nab, ir
  REAL(DP) :: fact
  COMPLEX(DP) :: cfac
  ! work space
  COMPLEX(DP), POINTER :: aux_d(:), aux1_d(:,:,:), vg_d(:,:), qgm_d(:,:)
  REAL(DP),    POINTER :: ddeeq_d(:,:,:,:), qmod_d(:), ylmk0_d(:,:) 
  REAL(DP),    ALLOCATABLE ::  forceq(:,:)
  INTEGER,POINTER          :: nl_d(:)  
  INTEGER                  :: ierr 
  !
  REAL(DP)                 :: forceqx, forceqy, forceqz
#if defined(__CUDA) 
ATTRIBUTES (DEVICE) aux_d, aux1_d, vg_d, qgm_d, ddeeq_d, qmod_d, ylmk0_d,nl_d 
#endif 
  nl_d => dfftp%nl_d
  IF (.NOT.okvan) RETURN
  !
  ALLOCATE( forceq(3,nat) )
  forceq(:,:) = 0.0_dp
  IF ( gamma_only ) THEN
     fact = 2.d0*omega
  ELSE
     fact = 1.d0*omega
  ENDIF
  !
  ! fourier transform of the total effective potential
  !
  CALL dev_buf%lock_buffer( vg_d, [ngm, nspin_mag] , ierr )
  IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
  !
  CALL dev_buf%lock_buffer( aux_d, dfftp%nnr, ierr )
  IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
  !
  DO is = 1, nspin_mag
     IF (nspin_mag==4.AND.is/=1) THEN
        aux_d(:) = v%of_r(:,is)
     ELSE
        aux_d(:) = vltot (:) + v%of_r (:, is)
     ENDIF
     CALL fwfft( 'Rho', aux_d, dfftp )
     ! Note the factors -i and 2pi/a *units of G) here in V(G) !
     !
     !$cuf kernel do
     do ir=1, ngm  
        vg_d(ir, is) = aux_d(nl_d(ir)) * tpiba * (0.d0, -1.d0)
     end do
  ENDDO
  CALL dev_buf%release_buffer( aux_d, ierr )
  !
  ! With k-point parallelization, distribute G-vectors across processors
  ! ngm_s = index of first G-vector for this processor
  ! ngm_e = index of last  G-vector for this processor
  ! ngm_l = local number of G-vectors 
  !
  CALL divide( inter_pool_comm, ngm, ngm_s, ngm_e )
  ngm_l = ngm_e-ngm_s+1
  ! for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
  CALL dev_buf%lock_buffer( ylmk0_d, [ngm_l,lmaxq*lmaxq], ierr )
  IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
  !
  CALL ylmr2_gpu( lmaxq * lmaxq, ngm_l, g_d(1,ngm_s), gg_d(ngm_s), ylmk0_d )
  !
  CALL dev_buf%lock_buffer( qmod_d, ngm_l, ierr  )
  IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
  !
  !$cuf kernel do 
  DO ig = 1, ngm_l
     qmod_d(ig) = SQRT( gg_d(ngm_s+ig-1) )*tpiba
  ENDDO
  !
  ! Sync if needed
  CALL using_becsum_d(0)
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        ! qgm contains the Q functions in G space
        !
        nij = nh(nt)*(nh(nt)+1)/2
        CALL dev_buf%lock_buffer(qgm_d, [ngm_l, nij],ierr) 
        IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
        !
        ijh = 0
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              ijh = ijh + 1
              CALL qvan2_gpu( ngm_l, ih, jh, nt, qmod_d, qgm_d(1,ijh), ylmk0_d )
           ENDDO
        ENDDO
        !
        ! nab = number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        CALL dev_buf%lock_buffer(aux1_d, [ngm_l, nab, 3], ierr )
        IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
        !
        CALL dev_buf%lock_buffer( ddeeq_d, [nij, nab, 3, nspin_mag],ierr )
        IF (ierr /= 0) CALL errore( 'addusforce_gpu', 'cannot allocate buffers', -1 )
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 !
                 ! aux1 = product of potential, structure factor and iG
                 !
                 !$cuf kernel do 
                 do ig = 1, ngm_l
                    cfac = vg_d(ngm_s+ig-1, is) * &
                         CONJG(eigts1_d(mill_d(1,ngm_s+ig-1),na) * &
                               eigts2_d(mill_d(2,ngm_s+ig-1),na) * &
                               eigts3_d(mill_d(3,ngm_s+ig-1),na) )
                    aux1_d(ig, nb, 1) = g_d(1,ngm_s+ig-1) * cfac
                    aux1_d(ig, nb, 2) = g_d(2,ngm_s+ig-1) * cfac
                    aux1_d(ig, nb, 3) = g_d(3,ngm_s+ig-1) * cfac
                 enddo
                 !
              ENDIF
           ENDDO
           !
           !    ddeeq = dot product of aux1 with the Q functions
           !    No need for special treatment of the G=0 term (is zero)
           !
           DO ipol = 1, 3
              CALL cublasDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm_d, 2*ngm_l, &
                   aux1_d(1,1,ipol), 2*ngm_l, 0.0_dp, ddeeq_d(1,1,ipol,is), nij )
           ENDDO
           !
        ENDDO
        !
        CALL dev_buf%release_buffer(aux1_d, ierr)
        CALL dev_buf%release_buffer(qgm_d,  ierr)
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 !DO ipol = 1, 3
                 forceqx = 0
                 forceqy = 0
                 forceqz = 0
                 !$cuf kernel do
                    DO ijh = 1, nij 
                       forceqx = forceqx + ddeeq_d(ijh, nb, 1, is) * becsum_d(ijh, na, is)
                       forceqy = forceqy + ddeeq_d(ijh, nb, 2, is) * becsum_d(ijh, na, is)                    
                       forceqz = forceqz + ddeeq_d(ijh, nb, 3, is) * becsum_d(ijh, na, is)   
                    ENDDO
                forceq(1,na) = forceq(1,na)+forceqx
                forceq(2,na) = forceq(2,na)+forceqy
                forceq(3,na) = forceq(3,na)+forceqz
              ENDIF
           ENDDO
        ENDDO
        CALL dev_buf%release_buffer( ddeeq_d , ierr)
     ENDIF
  ENDDO
  !
  10 CONTINUE
  CALL mp_sum( forceq, inter_pool_comm )
  CALL mp_sum( forceq, intra_bgrp_comm )
  !
  forcenl(:,:) = forcenl(:,:) + forceq(:,:)
  !
  CALL dev_buf%release_buffer(qmod_d, ierr)
  CALL dev_buf%release_buffer(ylmk0_d, ierr)
  CALL dev_buf%release_buffer (vg_d, ierr )
  DEALLOCATE(forceq)
  !
  RETURN
  !
END SUBROUTINE addusforce_g_gpu
