!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
#define DEVICEATTR ,DEVICE
#else
#define DEVICEATTR
#endif
!
#if defined(__CUDA)
#define PINMEM 
#else
#define PINMEM
#endif
!
!-----------------------------------------------------------------------
   subroutine nlsm1us_x ( n, betae, c, becp )
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_bgrp, intra_bgrp_comm
      USE gvecw,      only : ngw
      USE uspp,       only : nkb
      USE gvect,      ONLY : gstart
#if defined(__CUDA)
      USE cudafor
      USE cublas
#endif
!
      implicit none

      integer,     intent(in)  :: n
      complex(DP) DEVICEATTR , intent(in)  :: c( :, : )
      complex(DP) DEVICEATTR , intent(inout)  :: betae( :, : )
      real(DP)    DEVICEATTR , intent(out) :: becp( :, : )
      INTEGER :: i
      !
      call start_clock( 'nlsm1us' )

      IF( ngw > 0 .AND. nkb > 0 ) THEN
         IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, nkb
               betae( 1, i ) = 0.5d0 * betae( 1, i )
            END DO
         END IF
         CALL MYDGEMM( 'T', 'N', nkb, n, 2*ngw, 2.0d0, betae, 2*ngw, c, 2*ngw, 0.0d0, becp, nkb )
         IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, nkb
               betae( 1, i ) = 2.0d0 * betae( 1, i )
            END DO
         END IF
      END IF

      IF( nproc_bgrp > 1 ) THEN
        CALL mp_sum( becp, intra_bgrp_comm )
      END IF

      call stop_clock( 'nlsm1us' )

      return
   end subroutine nlsm1us_x
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   subroutine nlsm1_x ( n, betae, c, becp, pptype_ )
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_bgrp, intra_bgrp_comm
      USE ions_base,  only : nat, nsp, ityp
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta, ofsbeta
      USE uspp_param, only : nh, upf, nhm
      USE gvect,      ONLY : gstart
!
      implicit none

      integer,     intent(in)  :: n
      complex(DP), intent(in)  :: c( :, : )
      complex(DP), intent(inout)  :: betae( :, : )
      real(DP),    intent(out) :: becp( :, : )
      INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      ! pptype_: pseudo type to process: 0 = all, 1 = norm-cons, 2 = ultra-soft
      !
      integer :: ig, is, iv, ia, l, inl
      integer :: pptype
      real(DP), allocatable :: becps( :, : )
      LOGICAL :: nothing_to_do
      !
      call start_clock( 'nlsm1' )

      IF( PRESENT( pptype_ ) ) THEN
         pptype = pptype_
      ELSE
         pptype = 0
      END IF

      IF( pptype == 1 ) THEN
         nothing_to_do = .TRUE.
         do is = 1, nsp
            IF( .NOT. upf(is)%tvanp ) THEN
               nothing_to_do = .FALSE.
            END IF
         END DO
         IF( nothing_to_do ) GO TO 100
      END IF

      IF( pptype == 0 ) THEN
         IF( ngw > 0 .AND. nkb > 0 ) THEN
            IF( gstart > 1 ) betae( 1, : ) = 0.5d0 * betae( 1, : )
            CALL dgemm( 'T', 'N', nkb, n, 2*ngw, 2.0d0, betae, 2*ngw, c, 2*ngw, 0.0d0, becp, SIZE(becp,1) )
            IF( gstart > 1 ) betae( 1, : ) = 2.0d0 * betae( 1, : )
         END IF
         IF( nproc_bgrp > 1 ) THEN
            CALL mp_sum( becp, intra_bgrp_comm )
         END IF
         GO TO 100
      END IF
      
      allocate( becps( SIZE(becp,1), SIZE(becp,2) ) ) 
 
      IF( ngw > 0 .AND. nkb > 0 ) THEN
         IF( gstart > 1 ) betae( 1, : ) = 0.5d0 * betae( 1, : )
         CALL dgemm( 'T', 'N', nkb, n, 2*ngw, 2.0d0, betae, 2*ngw, c, 2*ngw, 0.0d0, becps, nkb )
         IF( gstart > 1 ) betae( 1, : ) = 2.0d0 * betae( 1, : )
      END IF

      IF( nproc_bgrp > 1 ) THEN
        CALL mp_sum( becps, intra_bgrp_comm )
      END IF
      do is = 1, nsp
        IF( pptype == 2 .AND. .NOT. upf(is)%tvanp ) CYCLE
        IF( pptype == 1 .AND. upf(is)%tvanp ) CYCLE
          DO ia = 1, nat
            IF( ityp(ia) == is ) THEN
              inl = ofsbeta(ia)
              do iv = 1, nh( is )
                becp(inl+iv,:) = becps( inl+iv, : )
              end do
            END IF
         end do
      end do
              !
      DEALLOCATE( becps )

100   CONTINUE

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_x
!-----------------------------------------------------------------------

#if defined (__CUDA)
!-----------------------------------------------------------------------
   subroutine nlsm1_gpu_x ( n, betae, c, becp, pptype_ )
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_bgrp, intra_bgrp_comm
      USE ions_base,  only : nat, nsp, ityp
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta, ofsbeta
      USE uspp_param, only : nh, upf, nhm
      USE gvect,      ONLY : gstart
      USE cp_main_variables, ONLY : nlsm1_wrk_d
!
      implicit none

      integer,     intent(in)  :: n
      complex(DP), intent(in)  :: c( :, : )
      complex(DP), intent(inout)  :: betae( :, : )
      real(DP),    intent(out) :: becp( :, : )
#if defined(__CUDA)
      attributes(DEVICE) :: c, betae, becp
#endif
      INTEGER,     INTENT(IN), OPTIONAL  :: pptype_
      ! pptype_: pseudo type to process: 0 = all, 1 = norm-cons, 2 = ultra-soft
      !
      integer :: ig, is, iv, ia, l, inl
      integer :: pptype
      integer :: i
      LOGICAL :: nothing_to_do
      !
      call start_clock( 'nlsm1' )

      IF( PRESENT( pptype_ ) ) THEN
         pptype = pptype_
      ELSE
         pptype = 0
      END IF

      IF( pptype == 1 ) THEN
         nothing_to_do = .TRUE.
         do is = 1, nsp
            IF( .NOT. upf(is)%tvanp ) THEN
               nothing_to_do = .FALSE.
            END IF
         END DO
         IF( nothing_to_do ) GO TO 100
      END IF

      IF( pptype == 0 ) THEN

         IF( ngw > 0 .AND. nkb > 0 ) THEN
            IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, SIZE(betae,2)
                  betae( 1, i ) = 0.5d0 * betae( 1, i )
               END DO
            END IF
            CALL MYDGEMM( 'T', 'N', nkb, n, 2*ngw, 2.0d0, betae, 2*ngw, c, 2*ngw, 0.0d0, becp, SIZE(becp,1) )
            IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, SIZE(betae,2)
                  betae( 1, i ) = 2.0d0 * betae( 1, i )
               END DO
            END IF
         END IF

         IF( nproc_bgrp > 1 ) THEN
            CALL mp_sum( becp, intra_bgrp_comm )
         END IF

         GO TO 100

      END IF

      IF( ngw > 0 .AND. nkb > 0 ) THEN
         IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, SIZE(betae,2)
               betae( 1, i ) = 0.5d0 * betae( 1, i )
            END DO
         END IF
         CALL MYDGEMM( 'T', 'N', nkb, n, 2*ngw, 2.0d0, betae, 2*ngw, c, 2*ngw, 0.0d0, nlsm1_wrk_d, nkb )
         IF( gstart > 1 ) THEN
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, SIZE(betae,2)
               betae( 1, i ) = 2.0d0 * betae( 1, i )
            END DO
         END IF
      END IF

      IF( nproc_bgrp > 1 ) THEN
        CALL mp_sum( nlsm1_wrk_d, intra_bgrp_comm )
      END IF

      do is = 1, nsp
        IF( pptype == 2 .AND. .NOT. upf(is)%tvanp ) CYCLE
        IF( pptype == 1 .AND. upf(is)%tvanp ) CYCLE
          DO ia = 1, nat
            IF( ityp(ia) == is ) THEN
              inl = ofsbeta(ia)
!$cuf kernel do(2) <<<*,*>>>
              do i = 1, SIZE(becp,2)
                do iv = 1, nh( is )
                  becp(inl+iv,i) = nlsm1_wrk_d(inl+iv,i)
                end do
              end do
            END IF
         end do
      end do
              !
100   CONTINUE

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_gpu_x
!-----------------------------------------------------------------------
#endif

!-------------------------------------------------------------------------
   subroutine nlsm2_bgrp_x( ngw, nkb, betae, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
!-----------------------------------------------------------------------

      !     computes: the array becdr
      !     becdr(ia,n,iv,is,k)
      !      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
      !     input : eigr, c
      !     output: becdr
      !
 
      USE kinds,      ONLY : DP
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_bgrp, intra_bgrp_comm
      use gvect,      only : g, gstart
!
      implicit none
    
      integer,     intent(in)  :: ngw, nkb, nbspx_bgrp, nbsp_bgrp
      complex(DP), intent(in)  :: betae(:,:), c_bgrp(:,:)
      real(DP),    intent(out) :: becdr_bgrp(:,:,:)
      !
      complex(DP), allocatable :: wrk2(:,:)
      !
      integer  :: ig, iv, k, info
      complex(DP) :: cfact1, cfact2
!
      call start_clock( 'nlsm2' )

      cfact2 = - cmplx( 0.0_dp , 1.0_dp ) * tpiba * 2.0d0
      cfact1 = - cmplx( 0.0_dp , 1.0_dp ) * tpiba 

      allocate( wrk2, MOLD = betae, STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' nlsm2 ', ' allocating wrk2', ABS( info ) )
!
      DO k = 1, 3
         do iv=1,nkb
            wrk2(1,iv) = cfact1 * g(k,1) * betae(1,iv)
         end do
!$omp parallel do default(shared) private(iv,ig) collapse(2)
         do iv=1,nkb
            do ig=gstart,ngw
               wrk2(ig,iv) = cfact2 * g(k,ig) * betae(ig,iv)
            end do
         end do
!$omp end parallel do
         IF( ngw > 0 .AND. nkb > 0 ) THEN
            CALL dgemm( 'T', 'N', nkb, nbsp_bgrp, 2*ngw, 1.0d0, wrk2(1,1), 2*ngw, &
                 c_bgrp, 2*ngw, 0.0d0, becdr_bgrp( 1, 1, k ), nkb )
         END IF
      end do

      deallocate( wrk2 )

      IF( nproc_bgrp > 1 ) THEN
         CALL mp_sum( becdr_bgrp, intra_bgrp_comm )
      END IF

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2_bgrp_x
!-----------------------------------------------------------------------

#if defined (__CUDA)
!-------------------------------------------------------------------------
   subroutine nlsm2_bgrp_gpu_x( ngw, nkb, betae, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
!-----------------------------------------------------------------------

      !     computes: the array becdr
      !     becdr(ia,n,iv,is,k)
      !      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
      !     input : eigr, c
      !     output: becdr
      !
 
      USE kinds,      ONLY : DP
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_bgrp, intra_bgrp_comm
      use gvect,      only : gstart, g
      USE device_memcpy_m, ONLY : dev_memcpy
      USE cudafor
      USE cublas
!
      implicit none
    
      integer,     intent(in)  :: ngw, nkb, nbspx_bgrp, nbsp_bgrp
      complex(DP), intent(in), DEVICE :: c_bgrp(:,:)
      complex(DP), intent(in), DEVICE :: betae(:,:)
      real(DP),    intent(out) PINMEM :: becdr_bgrp(:,:,:)
      !
      complex(DP), allocatable, DEVICE :: wrk2(:,:)
      real(DP), allocatable, DEVICE :: becdr_d(:,:)
      !
      integer  :: ig, iv, k, info
      complex(DP) :: cfact1, cfact2
!
      call start_clock( 'nlsm2' )

      ALLOCATE( wrk2, MOLD=betae, STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' nlsm2 ', ' allocating wrk2', ABS( info ) )
      ALLOCATE( becdr_d( SIZE( becdr_bgrp, 1 ), SIZE( becdr_bgrp, 2 ) ), STAT=info ) 
      IF( info /= 0 ) &
         CALL errore( ' nlsm2 ', ' allocating becdr_d ', ABS( info ) )

      cfact2 = - cmplx( 0.0_dp , 1.0_dp ) * tpiba * 2.0d0
      cfact1 = - cmplx( 0.0_dp , 1.0_dp ) * tpiba 

      DO k = 1, 3
!$acc parallel present(g(:,:)) deviceptr(wrk2, betae)
!$acc loop gang 
         do iv=1,nkb
            wrk2(1,iv) = cfact1 * g(k,1) * betae(1,iv)
!$acc loop vector
            do ig=gstart,ngw
               wrk2(ig,iv) = cfact2 * g(k,ig) * betae(ig,iv)
            end do
         end do
!$acc end parallel
         IF( ngw > 0 .AND. nkb > 0 ) THEN
            CALL MYDGEMM( 'T', 'N', nkb, nbsp_bgrp, 2*ngw, 1.0d0, wrk2(1,1), 2*ngw, &
                 c_bgrp, 2*ngw, 0.0d0, becdr_d, nkb )
            CALL dev_memcpy( becdr_bgrp(:,:,k), becdr_d )
         END IF
      end do

      DEALLOCATE( becdr_d )
      deallocate( wrk2 )

      IF( nproc_bgrp > 1 ) THEN
         CALL mp_sum( becdr_bgrp, intra_bgrp_comm )
      END IF

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2_bgrp_gpu_x
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
   SUBROUTINE ennl_x( ennl_val, rhovan, bec_bgrp )
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use uspp_param,     only : nh, upf
      use uspp,           only : dvan, ofsbeta
      use electrons_base, only : nbsp_bgrp, nspin, ispin_bgrp, f_bgrp, nbspx_bgrp
      use ions_base,      only : nsp, nat, ityp
      !
      implicit none
      !
      ! input
      !
      real(DP), intent(out) :: ennl_val
      real(DP), intent(out) :: rhovan( :, :, : )
      real(DP), intent(in)  :: bec_bgrp( :, : )
      !
      ! local
      !
      real(DP) :: sumt, sums(2), ennl_t
      integer  :: is, iv, jv, ijv, inl, jnl, ia, iss, i, indv
      INTEGER  :: omp_get_num_threads
      !
      ennl_t = 0.d0  
      !
!$omp parallel num_threads(min(4,omp_get_num_threads())) default(none) &
!$omp shared(nat,ityp,ofsbeta,nh,nbsp_bgrp,ispin_bgrp,f_bgrp,bec_bgrp,rhovan,dvan,nspin,ennl_t) &
!$omp private(ia,is,indv,iv,inl,jv,ijv,jnl,sums,iss,i,sumt)
!$omp do reduction(+:ennl_t)
      do ia = 1, nat
         is   = ityp(ia)
         indv = ofsbeta(ia)
         do iv = 1, nh(is)
            inl = indv + iv
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               jnl = indv + jv
               sums = 0.d0
               do i = 1, nbsp_bgrp
                  iss = ispin_bgrp(i)
                  sums(iss) = sums(iss) + f_bgrp(i) * bec_bgrp(inl,i) * bec_bgrp(jnl,i)
               end do
               sumt = 0.d0
               do iss = 1, nspin
                  rhovan( ijv, ia, iss ) = sums( iss )
                  sumt = sumt + sums( iss )
               end do
               if( iv .ne. jv ) sumt = 2.d0 * sumt
               ennl_t = ennl_t + sumt * dvan( jv, iv, is)
            end do
         end do
      end do
!$omp end do
!$omp end parallel
      !
      ennl_val = ennl_t
      !
      return
   end subroutine ennl_x
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   subroutine calrhovan_x( rhovan, bec, iwf )
!-----------------------------------------------------------------------
      !
      ! calculation of rhovan relative to state iwf
      !
      use kinds,          only : DP
      use uspp_param,     only : nh
      use uspp,           only : ofsbeta
      use electrons_base, only : ispin, f
      use ions_base,      only : nat, ityp
      !
      implicit none
      !
      ! input
      !
      real(DP), intent(out) :: rhovan( :, :, : )
      real(DP), intent(in) :: bec( :, : )
      integer, intent(in) :: iwf
      !
      ! local
      !
      integer   :: is, iv, jv, ijv, inl, jnl, ia, iss
      !
      iss = ispin(iwf)
      !
      do ia = 1, nat
         is = ityp(ia)
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               inl = ofsbeta(ia) + iv
               jnl = ofsbeta(ia) + jv
               rhovan( ijv, ia, iss ) = f(iwf) * bec(inl,iwf) * bec(jnl,iwf)
            end do
         end do
      end do
      !
      return
   end subroutine calrhovan_x
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine calbec_x ( n, betae, c, bec, pptype_ )
!-----------------------------------------------------------------------

      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !
      
      USE kinds,          ONLY : DP
      use cp_interfaces,  only : nlsm1
!
      implicit none
      !
      INTEGER,     INTENT(IN)    :: n
      real(DP),    intent(out)   :: bec( :, : )
      complex(DP), intent(in)    :: c( :, : )
      complex(DP), intent(inout) :: betae( :, : )
      INTEGER,     INTENT(IN), OPTIONAL  :: pptype_

      call start_clock( 'calbec' )
      call nlsm1( n, betae, c, bec, pptype_ )
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec_x
!-----------------------------------------------------------------------

#if defined (__CUDA)
!-----------------------------------------------------------------------
   subroutine calbec_gpu_x ( n, betae, c, bec, pptype_ )
!-----------------------------------------------------------------------

      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !

      USE kinds,          ONLY : DP
      use cp_interfaces,  only : nlsm1
!
      implicit none
      !
      INTEGER,     INTENT(IN)    :: n
      real(DP),    intent(out), device   :: bec( :, : )
      complex(DP), intent(in), device    :: c( :, : )
      complex(DP), intent(inout), device :: betae( :, : )
      INTEGER,     INTENT(IN), OPTIONAL  :: pptype_

      call start_clock( 'calbec' )
      call nlsm1( n, betae, c, bec, pptype_ )
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec_gpu_x
!-----------------------------------------------------------------------
#endif


!-----------------------------------------------------------------------
SUBROUTINE dbeta_eigr_x( dbeigr, eigr )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  use ions_base,  only : nat, ityp
  use uspp,       only : nhtol, nkb, dbeta, ofsbeta
  use uspp_param, only : nh, nhm
  use gvect,      only : gstart
  use gvecw,      only : ngw
  !
  implicit none
  !
  include 'laxlib.fh'
  !
  complex(DP), intent(out) :: dbeigr( :, :, :, : )
  complex(DP), intent(in)  :: eigr(:,:)
  !
  integer   :: ig, is, iv, ia, l, inl, i, j
  complex(DP) :: cfact(4)
  !
  call start_clock( 'dbeta_eigr' )
  !
  !if (l == 0) then
  cfact(1) =   cmplx( 1.0_dp , 0.0_dp )
  !else if (l == 1) then
  cfact(2) = - cmplx( 0.0_dp , 1.0_dp )
  !else if (l == 2) then
  cfact(3) = - cmplx( 0.0_dp , 1.0_dp )
  cfact(3) = cfact(3) * cfact(3)
  !else if (l == 3) then
  cfact(4) = - cmplx( 0.0_dp , 1.0_dp )
  cfact(4) = cfact(4) * cfact(4) * cfact(4)
  !endif

  do j=1,3
     do i=1,3
        do ia = 1, nat
           is = ityp(ia) 
           inl = ofsbeta(ia)
           do iv=1,nh(is)
              l=nhtol(iv,is)
              !     q = 0   component (with weight 1.0)
              dbeigr(1,iv+inl,i,j)= cfact(l+1)*dbeta(1,iv,is,i,j)*eigr(1,ia)
              !     q > 0   components (with weight 2.0)
              do ig = gstart, ngw
                 dbeigr(ig,iv+inl,i,j) = 2.0d0*cfact(l+1)*dbeta(ig,iv,is,i,j)*eigr(ig,ia)
              end do
           end do
        end do
     end do
  end do
  !
  call stop_clock( 'dbeta_eigr' )
  !
  return
end subroutine dbeta_eigr_x
!-----------------------------------------------------------------------

#if defined (__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE dbeta_eigr_gpu_x( dbeigr, eigr )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  use ions_base,  only : nat, ityp
  use uspp,       only : nhtol, nkb, dbeta, ofsbeta
  use uspp_param, only : nh, nhm
  use gvect,      only : gstart
  use gvecw,      only : ngw
  use device_memcpy_m, only : dev_memcpy
  !
  implicit none
  !
  include 'laxlib.fh'
  !
  complex(DP), device, intent(out) :: dbeigr( :, :, :, : )
  complex(DP), device, intent(in)  :: eigr(:,:)
  !
  integer   :: ig, is, iv, ia, l, inl, i, j
  complex(DP) :: cfact(4)
  !
  complex(DP), device :: cfact_d(4)
  integer, allocatable, device :: nhtol_d(:,:)
  real(DP), allocatable, device :: dbeta_d(:,:,:,:,:)
  !
  call start_clock( 'dbeta_eigr' )
  !
  !if (l == 0) then
  cfact(1) =   cmplx( 1.0_dp , 0.0_dp )
  !else if (l == 1) then
  cfact(2) = - cmplx( 0.0_dp , 1.0_dp )
  !else if (l == 2) then
  cfact(3) = - cmplx( 0.0_dp , 1.0_dp )
  cfact(3) = cfact(3) * cfact(3)
  !else if (l == 3) then
  cfact(4) = - cmplx( 0.0_dp , 1.0_dp )
  cfact(4) = cfact(4) * cfact(4) * cfact(4)
  !endif

  call dev_memcpy(cfact_d, cfact)

  allocate(nhtol_d, MOLD=nhtol)
  call dev_memcpy(nhtol_d, nhtol)

  allocate(dbeta_d, SOURCE=dbeta)

  do ia = 1, nat
     is = ityp(ia)
     inl = ofsbeta(ia)

!$cuf kernel do(3) <<<*,*>>>
     do j=1,3
        do i=1,3
           do iv=1,nh(is)
              l=nhtol_d(iv,is)
              ! q = 0   component (with weight 1.0)
              dbeigr(1,iv+inl,i,j)= cfact_d(l+1)*dbeta_d(1,iv,is,i,j)*eigr(1,ia)
           end do
        end do
     end do

!$cuf kernel do(4) <<<*,*>>>
     do j=1,3
        do i=1,3
           do iv=1,nh(is)
              ! q > 0   components (with weight 2.0)
              do ig = gstart, ngw
                 l=nhtol_d(iv,is)
                 dbeigr(ig,iv+inl,i,j) = 2.0d0*cfact_d(l+1)*dbeta_d(ig,iv,is,i,j)*eigr(ig,ia)
              end do
           end do
        end do
     end do

  end do
  !
  deallocate(nhtol_d)
  deallocate(dbeta_d)
  !
  call stop_clock( 'dbeta_eigr' )
  !
  return
end subroutine dbeta_eigr_gpu_x
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
SUBROUTINE caldbec_bgrp_x( eigr, c_bgrp, dbec, idesc )
  !-----------------------------------------------------------------------
  !
  !     this routine calculates array dbec, derivative of bec:
  !
  !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
  !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
  !
  !     with respect to cell parameters h
  !
  !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,      ONLY : DP
  use mp,         only : mp_sum
  use mp_global,  only : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm, nbgrp
  use ions_base,  only : nat, ityp
  use uspp,       only : nhtol, nkb, dbeta, ofsbeta
  use uspp_param, only : nh, nhm
  use gvect,      only : gstart
  use gvecw,      only : ngw
  use electrons_base, only : nspin, iupdwn, nupdwn, nbspx_bgrp, iupdwn_bgrp, nupdwn_bgrp, &
                             ibgrp_g2l, i2gupdwn_bgrp, nbspx, nbsp_bgrp
  use cp_interfaces,  only : dbeta_eigr
  !
  implicit none
  !
  include 'laxlib.fh'
  !
  complex(DP), intent(in)  :: c_bgrp( :, : )
  complex(DP), intent(in)  :: eigr(:,:)
  real(DP),    intent(out) :: dbec( :, :, :, : )
  integer, intent(in) :: idesc( :, : )
  !
  complex(DP), allocatable :: wrk2(:,:,:,:)
  real(DP),    allocatable :: dwrk_bgrp(:,:)
  !
  integer   :: ig, is, iv, ia, l, inl, i, j, ii, iw, iss, nr, ir, istart, nss
  integer   :: n1, n2, m1, m2, ibgrp_i, nrcx
  complex(DP) :: cfact
  !
  call start_clock( 'caldbec_bgrp' )
  !
  nrcx = MAXVAL(idesc(LAX_DESC_NRCX,:))
  !
  dbec = 0.0d0
  !
  allocate( wrk2( ngw, nkb, 3, 3 ) )
  allocate( dwrk_bgrp( nkb, nbspx_bgrp ) )
  !
  CALL dbeta_eigr( wrk2, eigr )
  !
  do j=1,3
     do i=1,3
        IF( ngw > 0 .AND. nkb > 0 ) THEN
           CALL dgemm( 'T', 'N', nkb, nbsp_bgrp, 2*ngw, 1.0d0, wrk2(1,1,i,j), 2*ngw, &
                             c_bgrp, 2*ngw, 0.0d0, dwrk_bgrp(1,1), nkb )
        END IF
        if( nproc_bgrp > 1 ) then
           call mp_sum( dwrk_bgrp, intra_bgrp_comm )
        end if
        do ia = 1, nat
           is = ityp(ia) 
           inl = ofsbeta(ia)
           do iss=1,nspin
              IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
                 nr = idesc( LAX_DESC_NR, iss )
                 ir = idesc( LAX_DESC_IR, iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do ii = 1, nr
                    ibgrp_i = ibgrp_g2l( ii + ir - 1 + istart - 1 )
                    IF( ibgrp_i > 0 ) THEN
                       do iw = 1, nh(is)
                          dbec( inl + iw, ii + (iss-1)*nrcx, i, j ) = dwrk_bgrp( inl + iw, ibgrp_i )
                       end do
                    END IF
                 end do
              END IF
           end do
        end do
     end do
  end do

  deallocate( wrk2 )
  deallocate( dwrk_bgrp )
  if( nbgrp > 1 ) then
     CALL mp_sum( dbec, inter_bgrp_comm )
  end if
  !
  call stop_clock( 'caldbec_bgrp' )
  !
  return
end subroutine caldbec_bgrp_x
!-----------------------------------------------------------------------

#if defined (__CUDA)
!-----------------------------------------------------------------------
SUBROUTINE caldbec_bgrp_gpu_x( eigr, c_bgrp, dbec, idesc )
  !-----------------------------------------------------------------------
  !
  !     this routine calculates array dbec, derivative of bec:
  !
  !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
  !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
  !
  !     with respect to cell parameters h
  !
  !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,      ONLY : DP
  use mp,         only : mp_sum
  use mp_global,  only : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm, nbgrp
  use ions_base,  only : nat, ityp
  use uspp,       only : nhtol, nkb, dbeta, ofsbeta
  use uspp_param, only : nh, nhm
  use gvect,      only : gstart
  use gvecw,      only : ngw
  use electrons_base, only : nspin, iupdwn, nupdwn, nbspx_bgrp, iupdwn_bgrp, nupdwn_bgrp, &
                             ibgrp_g2l, i2gupdwn_bgrp, nbspx, nbsp_bgrp
  use cp_interfaces,  only : dbeta_eigr
  use device_memcpy_m, only : dev_memcpy
  use cp_main_variables, only : dbec_d, caldbec_wrk_d, caldbec_dwrk_d
  !
  implicit none
  !
  include 'laxlib.fh'
  !
  complex(DP), intent(in), device :: c_bgrp( :, : )
  complex(DP), intent(in), device :: eigr(:,:)
  real(DP),    intent(out) PINMEM :: dbec( :, :, :, : )
  integer, intent(in) :: idesc( :, : )
  !
  integer   :: ig, is, iv, ia, l, inl, i, j, ii, iw, iss, nr, ir, istart, nss
  integer   :: n1, n2, m1, m2, ibgrp_i, nrcx
  complex(DP) :: cfact
  !

  integer :: ibgrp
  integer, allocatable, device :: ibgrp_g2l_d(:)
  !
  call start_clock( 'caldbec_bgrp' )
  !
  dbec_d = 0.0d0
  !
  allocate(ibgrp_g2l_d, MOLD=ibgrp_g2l)
  call dev_memcpy(ibgrp_g2l_d, ibgrp_g2l)
  !
  caldbec_wrk_d = (0.D0, 0.D0)
  !
  nrcx = MAXVAL(idesc(LAX_DESC_NRCX,:))
  !
  CALL dbeta_eigr( caldbec_wrk_d, eigr )
  !
  caldbec_dwrk_d = 0.D0
  !
  do j=1,3
     do i=1,3
        IF( ngw > 0 .AND. nkb > 0 ) THEN
           CALL MYDGEMM( 'T', 'N', nkb, nbsp_bgrp, 2*ngw, 1.0d0, caldbec_wrk_d(1,1,i,j), 2*ngw, &
                             c_bgrp, 2*ngw, 0.0d0, caldbec_dwrk_d(1,1), nkb )
        END IF
        if( nproc_bgrp > 1 ) then
           call mp_sum( caldbec_dwrk_d, intra_bgrp_comm )
        end if
        do ia = 1, nat
           is = ityp(ia)
           inl = ofsbeta(ia)
           do iss=1,nspin
              IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
                 nr = idesc( LAX_DESC_NR, iss )
                 ir = idesc( LAX_DESC_IR, iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
!$cuf kernel do(2) <<<*,*>>>
                 do ii = 1, nr
                    do iw = 1, nh(is)
                       ibgrp = ibgrp_g2l_d( ii + ir - 1 + istart - 1 )
                       IF( ibgrp > 0 ) THEN
                          dbec_d( inl + iw, ii + (iss-1)*nrcx, i, j ) = caldbec_dwrk_d( inl + iw, ibgrp )
                       END IF
                    end do
                 end do
              END IF
           end do
        end do
     end do
  end do
  !
  CALL dev_memcpy( dbec, dbec_d )
  !
  if( nbgrp > 1 ) then
     CALL mp_sum( dbec, inter_bgrp_comm )
  end if
  !
  call stop_clock( 'caldbec_bgrp' )
  !
  return
end subroutine caldbec_bgrp_gpu_x
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
subroutine dennl_x( bec_bgrp, dbec, drhovan, denl, idesc )
  !-----------------------------------------------------------------------
  !
  !  compute the contribution of the non local part of the
  !  pseudopotentials to the derivative of E with respect to h
  !
  USE kinds,      ONLY : DP
  use uspp_param, only : nh
  use uspp,       only : nkb, dvan, deeq, ofsbeta
  use ions_base,  only : nat, ityp
  use cell_base,  only : h
  use io_global,  only : stdout
  use mp,         only : mp_sum
  use mp_global,  only : intra_bgrp_comm
  use electrons_base,     only : nbspx_bgrp, nbsp_bgrp, ispin_bgrp, f_bgrp, nspin, iupdwn, nupdwn, ibgrp_g2l
  use gvect, only : gstart

  implicit none

  include 'laxlib.fh'

  real(DP), intent(in)  :: dbec( :, :, :, : )
  real(DP), intent(in)  :: bec_bgrp( :, : )
  real(DP), intent(out) :: drhovan( :, :, :, :, : )
  real(DP), intent(out) :: denl( 3, 3 )
  INTEGER, intent(in) :: idesc( :, : )

  real(DP) :: dsum(3,3),dsums(2,3,3), detmp(3,3)
  integer   :: is, iv, jv, ijv, inl, jnl, ia, iss, i,j,k
  integer   :: istart, nss, ii, ir, nr, ibgrp, nrcx
  !
  nrcx = MAXVAL(idesc(LAX_DESC_NRCX,:))
  !
  denl=0.d0
  drhovan=0.0d0

!$omp parallel default(none) &
!$omp shared(nat,ityp,ofsbeta,nh,nbsp_bgrp,ispin_bgrp,f_bgrp,bec_bgrp,drhovan,dvan,nspin,denl) &
!$omp shared(idesc,iupdwn,nupdwn,ibgrp_g2l,nrcx,dbec) &
!$omp private(ia,is,iv,inl,jv,ijv,jnl,dsums,iss,i,dsum,ii,ir,k,j,nr,istart,nss,ibgrp)
!$omp do reduction(+:denl)
  do ia=1,nat
     is = ityp(ia) 
     do iv=1,nh(is)
        do jv=iv,nh(is)
           ijv = (jv-1)*jv/2 + iv
           inl = ofsbeta(ia) + iv
           jnl = ofsbeta(ia) + jv
           dsums=0.d0
           do iss=1,nspin
              IF( ( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) .AND. &
                  ( idesc( LAX_DESC_MYR, iss ) == idesc( LAX_DESC_MYC, iss ) ) ) THEN
                 nr = idesc( LAX_DESC_NR, iss )
                 ir = idesc( LAX_DESC_IR, iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do i=1,nr
                    ii = i+istart-1+ir-1
                    ibgrp = ibgrp_g2l( ii )
                    IF( ibgrp > 0 ) THEN
                       do k=1,3
                          do j=1,3
                             dsums(iss,k,j) = dsums(iss,k,j) + f_bgrp(ibgrp) *       &
 &                          ( dbec(inl,i+(iss-1)*nrcx,k,j)*bec_bgrp(jnl,ibgrp)          &
 &                          + bec_bgrp(inl,ibgrp)*dbec(jnl,i+(iss-1)*nrcx,k,j) )
                          enddo
                       enddo
                    END IF
                 end do
                 dsum=0.d0
                 do k=1,3
                    do j=1,3
                       drhovan(ijv,ia,iss,j,k)=dsums(iss,j,k)
                       dsum(j,k)=dsum(j,k)+dsums(iss,j,k)
                    enddo
                 enddo
                 if(iv.ne.jv) dsum=2.d0*dsum
                 denl = denl + dsum * dvan(jv,iv,is)
              END IF
           end do
        end do
     end do
  end do
!$omp end do
!$omp end parallel

  CALL mp_sum( denl,    intra_bgrp_comm )
  CALL mp_sum( drhovan, intra_bgrp_comm )

!  WRITE(6,*) 'DEBUG enl (CP) = '
!  detmp = denl
!  detmp = MATMUL( detmp(:,:), TRANSPOSE( h ) )
!  WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
!5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
!     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
!     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!
  !
  return
end subroutine dennl_x
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine nlfq_bgrp_x( c_bgrp, betae, bec_bgrp, becdr_bgrp, fion )
  !-----------------------------------------------------------------------
  !
  !     contribution to fion due to nonlocal part
  !
  USE kinds,          ONLY : DP
  use uspp,           only : nkb, dvan, deeq, ofsbeta
  use uspp_param,     only : nhm, nh
  use ions_base,      only : nax, nat, ityp
  use electrons_base, only : nbsp_bgrp, f_bgrp, nbspx_bgrp, ispin_bgrp
  use gvecw,          only : ngw
  use constants,      only : pi, fpi
  use mp_global,      only : intra_bgrp_comm, nbgrp, inter_bgrp_comm, world_comm
  use mp_global,      only : me_bgrp, nproc_bgrp
  use mp,             only : mp_sum
  use cp_interfaces,  only : nlsm2_bgrp
  !
  implicit none
  !
  COMPLEX(DP), INTENT(IN) :: c_bgrp( :, : )
  COMPLEX(DP), INTENT(IN) ::  betae( :, : )
#if defined (__CUDA)
  ATTRIBUTES( DEVICE ) :: c_bgrp, betae
#endif
  REAL(DP),    INTENT(IN)  ::  bec_bgrp( :, : )
  REAL(DP),    INTENT(OUT)  ::  becdr_bgrp( :, :, : )
  REAL(DP),    INTENT(OUT) ::  fion( :, : )
  !
  integer  :: k, is, ia, inl, jnl, iv, jv, i
  real(DP) :: temp
  real(DP) :: sum_tmpdr
  !
  real(DP), allocatable :: tmpbec(:,:), tmpdr(:,:) 
  real(DP), allocatable :: fion_loc(:,:)
#if defined(_OPENMP) 
  INTEGER :: mytid, ntids, omp_get_thread_num, omp_get_num_threads
#endif  
  !
  call start_clock( 'nlfq' )
  !
  !     nlsm2 fills becdr
  !
  call nlsm2_bgrp( ngw, nkb, betae, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
  !
  allocate ( fion_loc( 3, nat ) )
  !
  fion_loc = 0.0d0
  !
!$omp parallel default(none), &
!$omp shared(becdr_bgrp,bec_bgrp,fion_loc,f_bgrp,deeq,dvan,nbsp_bgrp,ofsbeta,nh, &
!$omp        nat,nhm,nbspx_bgrp,ispin_bgrp,nproc_bgrp,me_bgrp,ityp), &
!$omp private(tmpbec,tmpdr,is,ia,iv,jv,k,inl,jnl,temp,i,mytid,ntids,sum_tmpdr)

#if defined(_OPENMP)
  mytid = omp_get_thread_num()  ! take the thread ID
  ntids = omp_get_num_threads() ! take the number of threads
#endif

  allocate ( tmpbec( nbspx_bgrp, nhm ), tmpdr( nbspx_bgrp, nhm ) )

  DO k = 1, 3
     DO ia = 1, nat
        is = ityp(ia)

        ! better if we distribute to MPI tasks too!
        !
        IF( MOD( ia + (k-1)*nat, nproc_bgrp ) /= me_bgrp ) CYCLE

#if defined(_OPENMP)
        ! distribute atoms round robin to threads
        !
        IF( MOD( ( ia + (k-1)*nat ) / nproc_bgrp, ntids ) /= mytid ) CYCLE
#endif  
        tmpbec = 0.d0
        do jv=1,nh(is)
           jnl = ofsbeta(ia) + jv
           do iv=1,nh(is)
              do i = 1, nbsp_bgrp
                 temp = dvan(iv,jv,is) + deeq(jv,iv,ia,ispin_bgrp( i ) )
                 tmpbec(i,iv) = tmpbec(i,iv) + temp * bec_bgrp(jnl,i)
              end do
           end do
        end do

        do iv = 1, nh(is)
           inl = ofsbeta(ia) + iv
           do i = 1, nbsp_bgrp
              tmpdr(i,iv) = f_bgrp( i ) * becdr_bgrp( inl, i, k )
           end do
        end do

        sum_tmpdr = 0.0d0
        do iv = 1, nh(is)
           do i = 1, nbsp_bgrp
              sum_tmpdr = sum_tmpdr + tmpdr(i,iv)*tmpbec(i,iv)
           end do
        end do

        fion_loc(k,ia) = fion_loc(k,ia)-2.d0*sum_tmpdr

     END DO
  END DO
  deallocate ( tmpbec, tmpdr )

!$omp end parallel
  !
  CALL mp_sum( fion_loc, intra_bgrp_comm )
  IF( nbgrp > 1 ) THEN
     CALL mp_sum( fion_loc, inter_bgrp_comm )
  END IF
  !
  fion = fion + fion_loc
  !
  !     end of x/y/z loop
  !
  deallocate ( fion_loc )
  !
  call stop_clock( 'nlfq' )
  !
  return
end subroutine nlfq_bgrp_x
