!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
   subroutine nlsm1_x ( n, nspmn, nspmx, eigr, c, becp )
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
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE uspp_param, only : nh, ish
      !
      USE gvect, ONLY : gstart
!
      implicit none

      integer,     intent(in)  :: n, nspmn, nspmx
      complex(DP), intent(in)  :: eigr( :, : ), c( :, : )
      real(DP), intent(out) :: becp( :, : )
      !
      integer   :: isa, ig, is, iv, ia, l, inl, i, nhx
      real(DP), allocatable :: becps( :, : )
      complex(DP), allocatable :: wrk2( :, : )
      complex(DP) :: cfact
      !
      call start_clock( 'nlsm1' )

      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         allocate( wrk2( ngw, na( is ) ) )
         !
         IF( nproc_bgrp > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) == 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
!$omp parallel default(shared), private(l,ig,ia,cfact)
            l = nhtol( iv, is )
            !
            if (l == 0) then
               cfact =   cmplx( 1.0_dp , 0.0_dp )
            else if (l == 1) then
               cfact = - cmplx( 0.0_dp , 1.0_dp )
            else if (l == 2) then
               cfact = - cmplx( 0.0_dp , 1.0_dp )
               cfact = cfact * cfact
            else if (l == 3) then
               cfact = - cmplx( 0.0_dp , 1.0_dp )
               cfact = cfact * cfact * cfact
            endif
!
!$omp do
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2( 1, ia ) = cfact * beta(1,iv,is) * eigr(1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  wrk2( ig, ia ) = 2.0d0 * cfact * beta(ig,iv,is) * eigr(ig,ia+isa)
               end do
               !
            end do
!$omp end do
            
!$omp end parallel
            !
            IF( nproc_bgrp > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do

         deallocate( wrk2 )


         IF( nproc_bgrp > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_bgrp_comm )

            do i = 1, n
               do iv = inl , ( inl + na(is) * nh(is) - 1 )
                  becp( iv, i ) = becps( iv - inl + 1, i )
               end do
            end do

            DEALLOCATE( becps )

         END IF

         isa = isa + na(is)

      end do

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_x
!-----------------------------------------------------------------------


!-------------------------------------------------------------------------
   subroutine nlsm2_bgrp_x( ngw, nkb, eigr, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
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
      use ions_base,  only : nsp, na, nat
      use uspp,       only : nhtol, beta
      use uspp_param, only : nh, ish
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_bgrp, intra_bgrp_comm
      use gvect,      only : g, gstart
!
      implicit none
    
      integer,     intent(in)  :: ngw, nkb, nbspx_bgrp, nbsp_bgrp
      complex(DP), intent(in)  :: eigr(:,:), c_bgrp(:,:)
      real(DP),    intent(out) :: becdr_bgrp(:,:,:)
      !
      real(DP),    allocatable :: gk(:)
      complex(DP), allocatable :: wrk2(:,:)
      !
      integer  :: ig, is, iv, ia, k, l, inl, isa, i
      real(DP) :: arg
      complex(DP) :: cfact
!
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )

      becdr_bgrp = 0.d0
!
      do k = 1, 3

         do ig=1,ngw
            gk(ig)=g(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            allocate( wrk2( ngw, na( is ) ) )

            do iv=1,nh(is)
               !
               !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
               !
!$omp parallel default(none), shared(na,nhtol,gstart,wrk2,gk,beta,eigr,ngw,iv,is,isa), private(l,cfact,ig,arg,ia)
               l=nhtol(iv,is)

               ! compute (-i)^(l+1)
               !
               if (l == 0) then
                  cfact = - cmplx( 0.0_dp , 1.0_dp )
               else if (l == 1) then
                  cfact = - cmplx( 0.0_dp , 1.0_dp )
                  cfact = cfact * cfact
               else if (l == 2) then
                  cfact = - cmplx( 0.0_dp , 1.0_dp )
                  cfact = cfact * cfact * cfact
               else if (l == 3) then
                  cfact =   cmplx( 1.0_dp , 0.0_dp )
               endif

!$omp do
               do ia=1,na(is)
                  !    q = 0   component (with weight 1.0)
                  if (gstart == 2) then
                     wrk2(1,ia) = cfact*gk(1)*beta(1,iv,is)*eigr(1,ia+isa)
                  end if
                  !    q > 0   components (with weight 2.0)
                  do ig=gstart,ngw
                     arg = 2.0d0*gk(ig)*beta(ig,iv,is)
                     wrk2(ig,ia) = cfact * arg * eigr(ig,ia+isa)
                  end do
               end do
!$omp end do
!$omp end parallel 
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), nbsp_bgrp, 2*ngw, 1.0d0, wrk2, 2*ngw, &
                           c_bgrp, 2*ngw, 0.0d0, becdr_bgrp( inl, 1, k ), nkb )
            end do

            deallocate( wrk2 )

            isa = isa + na(is)

         end do

      end do

      deallocate( gk )

      IF( nproc_bgrp > 1 ) THEN
         CALL mp_sum( becdr_bgrp, intra_bgrp_comm )
      END IF

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2_bgrp_x
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
   SUBROUTINE ennl_x( ennl_val, rhovan, bec_bgrp )
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use uspp_param,     only : nh, ish
      use uspp,           only : dvan
      use electrons_base, only : nbsp_bgrp, nspin, ispin_bgrp, f_bgrp, nbspx_bgrp
      use ions_base,      only : nsp, na
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
      integer  :: is, iv, jv, ijv, inl, jnl, isa, isat, ism, ia, iss, i
      !
      ennl_t = 0.d0  
      !
      !  xlf does not like name of function used for OpenMP reduction
      !
!$omp parallel default(shared), &
!$omp private(is,iv,jv,ijv,isa,isat,ism,ia,inl,jnl,sums,i,iss,sumt), reduction(+:ennl_t)
      do is = 1, nsp
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               isa = 0
               do ism = 1, is - 1
                  isa = isa + na(ism)
               end do
!$omp do
               do ia = 1, na(is)
                  inl = ish(is)+(iv-1)*na(is)+ia
                  jnl = ish(is)+(jv-1)*na(is)+ia
                  isat = isa+ia
                  sums = 0.d0
                  do i = 1, nbsp_bgrp
                     iss = ispin_bgrp(i)
                     sums(iss) = sums(iss) + f_bgrp(i) * bec_bgrp(inl,i) * bec_bgrp(jnl,i)
                  end do
                  sumt = 0.d0
                  do iss = 1, nspin
                     rhovan( ijv, isat, iss ) = sums( iss )
                     sumt = sumt + sums( iss )
                  end do
                  if( iv .ne. jv ) sumt = 2.d0 * sumt
                  ennl_t = ennl_t + sumt * dvan( jv, iv, is)
               end do
!$omp end do
            end do
         end do
      end do
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
      use uspp_param,     only : nhm, nh, ish
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
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
      integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss
      !
      do is = 1, nsp
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               isa = 0
               do ism = 1, is - 1
                  isa = isa + na(ism)
               end do
               do ia = 1, na(is)
                  inl = ish(is)+(iv-1)*na(is)+ia
                  jnl = ish(is)+(jv-1)*na(is)+ia
                  isa = isa+1
                  iss = ispin(iwf)
                  rhovan( ijv, isa, iss ) = f(iwf) * bec(inl,iwf) * bec(jnl,iwf)
               end do
            end do
         end do
      end do
      !
      return
   end subroutine calrhovan_x
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
   subroutine calbec_x ( nspmn, nspmx, eigr, c, bec )
!-----------------------------------------------------------------------

      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !
      
      USE kinds,          ONLY : DP
      use ions_base,      only : na
      use io_global,      only : stdout
      use electrons_base, only : nbsp
      use control_flags,  only : iverbosity
      use uspp_param,     only : nh, ish
      use cp_interfaces,  only : nlsm1
!
      implicit none
      !
      integer,     intent(in)  :: nspmn, nspmx
      real(DP),    intent(out) :: bec( :, : )
      complex(DP), intent(in)  :: c( :, : ), eigr( :, : )

      ! local variables

      integer :: is, ia, i , iv
!
      call start_clock( 'calbec' )
      !
      call nlsm1( nbsp, nspmn, nspmx, eigr, c, bec )
!
      if ( iverbosity > 1 ) then
         WRITE( stdout,*)
         do is=1,nspmx
            WRITE( stdout,'(33x,a,i4)') ' calbec: bec (is)',is
            do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calbec: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &             ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,nbsp)
            end do
         end do
      endif
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec_x
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine calbec_bgrp_x ( nspmn, nspmx, eigr, c_bgrp, bec_bgrp )
!-----------------------------------------------------------------------

      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !

      USE kinds,          ONLY : DP
      use ions_base,      only : na, nat
      use electrons_base, only : nbsp_bgrp, nbspx_bgrp
      use gvecw,          only : ngw
      use uspp_param,     only : nh, ish
      use uspp,           only : nkb
      use cp_interfaces,  only : nlsm1
!
      implicit none
      !
      integer,     intent(in)  :: nspmn, nspmx
      real(DP),    intent(out) :: bec_bgrp( :, : )
      complex(DP), intent(in)  :: c_bgrp( :, : ), eigr( :, : )
!
      call start_clock( 'calbec' )
      !
      call nlsm1( nbsp_bgrp, nspmn, nspmx, eigr, c_bgrp, bec_bgrp )
      !
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec_bgrp_x


!-----------------------------------------------------------------------
SUBROUTINE caldbec_bgrp_x( eigr, c_bgrp, dbec, descla )
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
  use ions_base,  only : na, nat, nsp
  use uspp,       only : nhtol, nkb, dbeta
  use uspp_param, only : nh, nhm, ish
  use gvect,      only : gstart
  use gvecw,      only : ngw
  USE descriptors,        ONLY : la_descriptor
  use electrons_base,     only : nspin, iupdwn, nupdwn, nbspx_bgrp, iupdwn_bgrp, nupdwn_bgrp, &
                                 ibgrp_g2l, i2gupdwn_bgrp, nbspx, nbsp_bgrp
  !
  implicit none
  !
  complex(DP), intent(in)  :: c_bgrp( :, : )
  complex(DP), intent(in)  :: eigr(:,:)
  real(DP),    intent(out) :: dbec( :, :, :, : )
  TYPE(la_descriptor), intent(in) :: descla( : )
  !
  complex(DP), allocatable :: wrk2(:,:)
  real(DP),    allocatable :: dwrk_bgrp(:,:)
  !
  integer   :: ig, is, iv, ia, l, inl, i, j, ii, isa, nanh, iw, iss, nr, ir, istart, nss
  integer   :: n1, n2, m1, m2, ibgrp_i, nrcx
  complex(DP) :: cfact
  !
  nrcx = MAXVAL(descla(:)%nrcx)
  !
  dbec = 0.0d0
  !
  do j=1,3
     do i=1,3

        isa = 0

        do is=1,nsp
           allocate( wrk2( ngw, na(is) ) )
           nanh = na(is)*nh(is)
           allocate( dwrk_bgrp( nanh, nbspx_bgrp ) )
           do iv=1,nh(is)
              l=nhtol(iv,is)
              if (l == 0) then
                 cfact =   cmplx( 1.0_dp , 0.0_dp )
              else if (l == 1) then
                 cfact = - cmplx( 0.0_dp , 1.0_dp )
              else if (l == 2) then
                 cfact = - cmplx( 0.0_dp , 1.0_dp )
                 cfact = cfact * cfact
              else if (l == 3) then
                 cfact = - cmplx( 0.0_dp , 1.0_dp )
                 cfact = cfact * cfact * cfact
              else
                 CALL errore(' caldbec  ', ' l not implemented ', ABS( l ) )
              endif
              !
              do ia=1,na(is)
                 if (gstart == 2) then
                    !     q = 0   component (with weight 1.0)
                    wrk2(1,ia)= cfact*dbeta(1,iv,is,i,j)*eigr(1,ia+isa)
                 end if
                 !     q > 0   components (with weight 2.0)
                 do ig = gstart, ngw
                    wrk2(ig,ia) = 2.0d0*cfact*dbeta(ig,iv,is,i,j)*eigr(ig,ia+isa)
                 end do
              end do
              inl=(iv-1)*na(is)+1
              CALL dgemm( 'T', 'N', na(is), nbsp_bgrp, 2*ngw, 1.0d0, wrk2, 2*ngw, c_bgrp, 2*ngw, 0.0d0, dwrk_bgrp(inl,1), nanh )
           end do
           deallocate( wrk2 )

           if( nproc_bgrp > 1 ) then
              call mp_sum( dwrk_bgrp, intra_bgrp_comm )
           end if

           inl=ish(is)+1
           do iss=1,nspin
              IF( descla( iss )%active_node > 0 ) THEN
                 nr = descla( iss )%nr
                 ir = descla( iss )%ir
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do ii = 1, nr
                    ibgrp_i = ibgrp_g2l( ii + ir - 1 + istart - 1 )
                    IF( ibgrp_i > 0 ) THEN
                       do iw = 1, nanh
                          dbec( iw + inl - 1, ii + (iss-1)*nrcx, i, j ) = dwrk_bgrp( iw, ibgrp_i )
                       end do
                    END IF
                 end do
              END IF
           end do
           deallocate( dwrk_bgrp )
           isa = isa + na(is)
        end do
     end do
  end do

  if( nbgrp > 1 ) then
     CALL mp_sum( dbec, inter_bgrp_comm )
  end if
  !
  return
end subroutine caldbec_bgrp_x
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dennl_x( bec_bgrp, dbec, drhovan, denl, descla )
  !-----------------------------------------------------------------------
  !
  !  compute the contribution of the non local part of the
  !  pseudopotentials to the derivative of E with respect to h
  !
  USE kinds,      ONLY : DP
  use uspp_param, only : nh, nhm, ish
  use uspp,       only : nkb, dvan, deeq
  use ions_base,  only : nsp, na, nat
  use cell_base,  only : h
  use io_global,  only : stdout
  use mp,         only : mp_sum
  use mp_global,  only : intra_bgrp_comm
  USE descriptors,        ONLY : la_descriptor
  use electrons_base,     only : nbspx_bgrp, nbsp_bgrp, ispin_bgrp, f_bgrp, nspin, iupdwn, nupdwn, ibgrp_g2l
  use gvect, only : gstart

  implicit none

  real(DP), intent(in)  :: dbec( :, :, :, : )
  real(DP), intent(in)  :: bec_bgrp( :, : )
  real(DP), intent(out) :: drhovan( :, :, :, :, : )
  real(DP), intent(out) :: denl( 3, 3 )
  TYPE(la_descriptor), intent(in) :: descla( : )

  real(DP) :: dsum(3,3),dsums(2,3,3), detmp(3,3)
  integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i,j,k
  integer   :: istart, nss, ii, ir, nr, ibgrp, nrcx
  !
  nrcx = MAXVAL(descla(:)%nrcx)
  !
  denl=0.d0
  drhovan=0.0d0

  do is=1,nsp
     do iv=1,nh(is)
        do jv=iv,nh(is)
           ijv = (jv-1)*jv/2 + iv
           isa=0
           do ism=1,is-1
              isa=isa+na(ism)
           end do
           do ia=1,na(is)
              inl=ish(is)+(iv-1)*na(is)+ia
              jnl=ish(is)+(jv-1)*na(is)+ia
              isa=isa+1
              dsums=0.d0
              do iss=1,nspin
                 IF( ( descla( iss )%active_node > 0 ) .AND. ( descla( iss )%myr == descla( iss )%myc ) ) THEN
                    nr = descla( iss )%nr
                    ir = descla( iss )%ir
                    istart = iupdwn( iss )
                    nss    = nupdwn( iss )
                    do i=1,nr
                       ii = i+istart-1+ir-1
                       ibgrp = ibgrp_g2l( ii )
                       IF( ibgrp > 0 ) THEN
                          do k=1,3
                             do j=1,3
                                dsums(iss,k,j)=dsums(iss,k,j)+f_bgrp(ibgrp)*       &
 &                          (dbec(inl,i+(iss-1)*nrcx,k,j)*bec_bgrp(jnl,ibgrp)          &
 &                          + bec_bgrp(inl,ibgrp)*dbec(jnl,i+(iss-1)*nrcx,k,j))
                             enddo
                          enddo
                       END IF
                    end do
                    dsum=0.d0
                    do k=1,3
                       do j=1,3
                          drhovan(ijv,isa,iss,j,k)=dsums(iss,j,k)
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
  end do

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
subroutine nlfq_bgrp_x( c_bgrp, eigr, bec_bgrp, becdr_bgrp, fion )
  !-----------------------------------------------------------------------
  !
  !     contribution to fion due to nonlocal part
  !
  USE kinds,          ONLY : DP
  use uspp,           only : nkb, dvan, deeq
  use uspp_param,     only : nhm, nh, ish, nvb
  use ions_base,      only : nax, nat, nsp, na
  use electrons_base, only : nbsp_bgrp, f_bgrp, nbspx_bgrp, ispin_bgrp
  use gvecw,          only : ngw
  use constants,      only : pi, fpi
  use mp_global,      only : intra_bgrp_comm, nbgrp, inter_bgrp_comm
  use mp,             only : mp_sum
  use cp_interfaces,  only : nlsm2_bgrp
  !
  implicit none
  !
  COMPLEX(DP), INTENT(IN)  ::  c_bgrp( :, : ), eigr( :, : )
  REAL(DP),    INTENT(IN)  ::  bec_bgrp( :, : )
  REAL(DP),    INTENT(OUT)  ::  becdr_bgrp( :, :, : )
  REAL(DP),    INTENT(OUT) ::  fion( :, : )
  !
  integer  :: k, is, ia, isa, inl, iv, jv, i
  real(DP) :: temp
  real(DP) :: sum_tmpdr
  !
  real(DP), allocatable :: tmpbec(:,:), tmpdr(:,:) 
  real(DP), allocatable :: fion_loc(:,:)
#if defined(__OPENMP) 
  INTEGER :: mytid, ntids, omp_get_thread_num, omp_get_num_threads
#endif  
  !
  call start_clock( 'nlfq' )
  !
  !     nlsm2 fills becdr
  !
  call nlsm2_bgrp( ngw, nkb, eigr, c_bgrp, becdr_bgrp, nbspx_bgrp, nbsp_bgrp )
  !
  allocate ( fion_loc( 3, nat ) )
  !
  fion_loc = 0.0d0
  !
  DO k = 1, 3

!$omp parallel default(none), &
!$omp shared(becdr_bgrp,bec_bgrp,fion_loc,k,f_bgrp,deeq,dvan,nbsp_bgrp,ish,nh,na,nsp,nhm,nbspx_bgrp,ispin_bgrp), &
!$omp private(tmpbec,tmpdr,isa,is,ia,iv,jv,inl,temp,i,mytid,ntids,sum_tmpdr)

#if defined(__OPENMP)
     mytid = omp_get_thread_num()  ! take the thread ID
     ntids = omp_get_num_threads() ! take the number of threads
#endif

     allocate ( tmpbec( nhm, nbspx_bgrp ), tmpdr( nhm, nbspx_bgrp ) )

     isa = 0
     !
     DO is=1,nsp
        DO ia=1,na(is)

           isa=isa+1

#if defined(__OPENMP)
           ! distribute atoms round robin to threads
           !
           IF( MOD( isa, ntids ) /= mytid ) CYCLE
#endif  

                 tmpbec = 0.d0
                 tmpdr  = 0.d0

                 do iv=1,nh(is)
                    do jv=1,nh(is)
                       inl=ish(is)+(jv-1)*na(is)+ia
                       do i = 1, nbsp_bgrp
                          temp = dvan(iv,jv,is) + deeq(jv,iv,isa,ispin_bgrp( i ) )
                          tmpbec(iv,i) = tmpbec(iv,i) + temp * bec_bgrp(inl,i)
                       end do
                    end do
                 end do

                 do iv=1,nh(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    do i = 1, nbsp_bgrp
                       tmpdr(iv,i) = f_bgrp( i ) * becdr_bgrp( inl, i, k )
                    end do
                 end do

                 sum_tmpdr = 0.0d0
                 do i = 1, nbsp_bgrp
                    do iv = 1, nh(is)
                       sum_tmpdr = sum_tmpdr + tmpdr(iv,i)*tmpbec(iv,i)
                    end do
                 end do

                 fion_loc(k,isa) = fion_loc(k,isa)-2.d0*sum_tmpdr

        END DO
     END DO
     deallocate ( tmpbec, tmpdr )
!$omp end parallel

  END DO
  !
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
