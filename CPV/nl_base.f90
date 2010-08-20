!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
   subroutine nlsm1 ( n, nspmn, nspmx, eigr, c, becp )
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
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
!
      implicit none

      integer,   intent(in)  :: n, nspmn, nspmx
      real(DP), intent(in)  :: eigr( 2, ngw, nat ), c( 2, ngw, n )
      real(DP), intent(out) :: becp( nkb, n )
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      real(DP), allocatable :: wrk2( :, :, : )
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
         allocate( wrk2( 2, ngw, na( is ) ) )
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) == 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
            l = nhtol( iv, is )
            !
            if (l == 0) then
               ixr = 1
               ixi = 2
               signre =  1.0d0
               signim =  1.0d0
            else if (l == 1) then
               ixr = 2
               ixi = 1
               signre =  1.0d0
               signim = -1.0d0
            else if (l == 2) then
               ixr = 1
               ixi = 2
               signre = -1.0d0
               signim = -1.0d0
            else if (l == 3) then
               ixr = 2
               ixi = 1
               signre = -1.0d0
               signim =  1.0d0
            endif
!
!$omp do
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2( 1, 1, ia ) = signre*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                  wrk2( 2, 1, ia ) = signim*beta(1,iv,is)*eigr(ixi,1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  arg = 2.0d0 * beta(ig,iv,is)
                  wrk2( 1, ig, ia ) = signre*arg*eigr(ixr,ig,ia+isa)
                  wrk2( 2, ig, ia ) = signim*arg*eigr(ixi,ig,ia+isa)
               end do
               !
            end do
!$omp end do
            
!$omp end parallel
            !
            IF( nproc_image > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do

         deallocate( wrk2 )


         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_image_comm )

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
   end subroutine nlsm1
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
   subroutine nlsm1_dist ( n, nspmn, nspmx, eigr, c, becp, nlax, nspin, desc )
!-----------------------------------------------------------------------
      !  
      ! This version is for becp distributed over procs  
      !  
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
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
      USE descriptors,        ONLY : descla_siz_ , lambda_node_ , nlar_ , ilar_ , la_n_
!
      implicit none

      integer,  intent(in)  :: n, nspmn, nspmx, nlax, nspin
      integer,  intent(in)  :: desc( descla_siz_ , nspin )
      real(DP), intent(in)  :: eigr( 2, ngw, nat ), c( 2, ngw, n )
      real(DP), intent(out) :: becp( nkb, nlax*nspin )
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      integer   :: nr, ir, nup
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      real(DP), allocatable :: wrk2( :, :, : )
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
         allocate( wrk2( 2, ngw, na( is ) ) )
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
            l = nhtol( iv, is )
            !
            if (l == 0) then
               ixr = 1
               ixi = 2
               signre =  1.0d0
               signim =  1.0d0
            else if (l == 1) then
               ixr = 2
               ixi = 1
               signre =  1.0d0
               signim = -1.0d0
            else if (l == 2) then
               ixr = 1
               ixi = 2
               signre = -1.0d0
               signim = -1.0d0
            else if (l == 3) then
               ixr = 2
               ixi = 1
               signre = -1.0d0
               signim =  1.0d0
            endif
!
!$omp do
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2( 1, 1, ia ) = signre*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                  wrk2( 2, 1, ia ) = signim*beta(1,iv,is)*eigr(ixi,1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  arg = 2.0d0 * beta(ig,iv,is)
                  wrk2( 1, ig, ia ) = signre*arg*eigr(ixr,ig,ia+isa)
                  wrk2( 2, ig, ia ) = signim*arg*eigr(ixi,ig,ia+isa)
               end do
               !
            end do
!$omp end do
            
!$omp end parallel
            
            !
            IF( nproc_image > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do

         deallocate( wrk2 )


         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_image_comm )

            IF( desc( lambda_node_ , 1 ) > 0 ) THEN
               ir = desc( ilar_ , 1 )
               nr = desc( nlar_ , 1 )
               do i = 1, nr
                  do iv = inl , ( inl + na(is) * nh(is) - 1 )
                     becp( iv, i ) = becps( iv - inl + 1, i + ir - 1 )
                  end do
               end do
            END IF
            !
            IF( nspin == 2 ) THEN
               IF( desc( lambda_node_ , 2 ) > 0 ) THEN
                  nup = desc( la_n_ , 1 )
                  ir = desc( ilar_ , 2 )
                  nr = desc( nlar_ , 2 )
                  do i = 1, nr
                     do iv = inl , ( inl + na(is) * nh(is) - 1 )
                        becp( iv, i + nlax ) = becps( iv - inl + 1, i + ir - 1 + nup )
                     end do
                  end do
               END IF
            END IF

            DEALLOCATE( becps )

         END IF

         isa = isa + na(is)

      end do

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_dist
!-----------------------------------------------------------------------


!-------------------------------------------------------------------------
   subroutine nlsm2( ngw, nkb, n, nspin, eigr, c, becdr )
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
      use uspp,       only : nhtol, beta  !, nkb
      use cvan,       only : ish
      use uspp_param, only : nh
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_image, intra_image_comm
      use cp_main_variables,  only : nlax, descla, distribute_bec
      use reciprocal_vectors, only : gx, gstart
!
      implicit none
    
      integer,  intent(in)  :: ngw, nkb, n, nspin
      real(DP), intent(in)  :: eigr(2,ngw,nat), c(2,ngw,n)
      real(DP), intent(out) :: becdr(nkb,nspin*nlax,3)
      !
      real(DP), allocatable :: gk(:)
      real(DP), allocatable :: wrk2(:,:,:)
      real(DP), allocatable :: becdr_repl(:,:)
      !
      integer   :: ig, is, iv, ia, k, l, ixr, ixi, inl, isa, i
      real(DP) :: signre, signim, arg
!
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )
      allocate( becdr_repl( nkb, n ) )

      becdr = 0.d0
!
      do k = 1, 3

         becdr_repl = 0.d0

         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            allocate( wrk2( 2, ngw, na( is ) ) )

            do iv=1,nh(is)
               !
               !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
               !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
               l=nhtol(iv,is)
               if (l.eq.0) then
                  ixr = 2
                  ixi = 1
                  signre =  1.0d0
                  signim = -1.0d0
               else if (l.eq.1) then
                  ixr = 1
                  ixi = 2
                  signre = -1.0d0
                  signim = -1.0d0
               else if (l.eq.2) then
                  ixr = 2
                  ixi = 1
                  signre = -1.0d0
                  signim =  1.0d0
               else if (l == 3) then
                  ixr = 1
                  ixi = 2
                  signre =  1.0d0
                  signim =  1.0d0
               endif
!    
!$omp do
               do ia=1,na(is)
                  !    q = 0   component (with weight 1.0)
                  if (gstart == 2) then
                     wrk2(1,1,ia) = signre*gk(1)*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                     wrk2(2,1,ia) = signim*gk(1)*beta(1,iv,is)*eigr(ixi,1,ia+isa)
                  end if
                  !    q > 0   components (with weight 2.0)
                  do ig=gstart,ngw
                     arg = 2.0d0*gk(ig)*beta(ig,iv,is)
                     wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                     wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                  end do
               end do
!$omp end do
!$omp end parallel 
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becdr_repl( inl, 1 ), nkb )
            end do

            deallocate( wrk2 )

            isa = isa + na(is)

         end do

         IF( nproc_image > 1 ) THEN
            CALL mp_sum( becdr_repl(:,:), intra_image_comm )
         END IF
         CALL distribute_bec( becdr_repl, becdr(:,:,k), descla, nspin )
      end do

      deallocate( gk )
      deallocate( becdr_repl )

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
   subroutine nlsm2_repl( ngw, nkb, n, eigr, c, becdr )
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
      use uspp,       only : nhtol, beta  !, nkb
      use cvan,       only : ish
      use uspp_param, only : nh
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_image, intra_image_comm
      use reciprocal_vectors, only : gx, gstart
!
      implicit none
    
      integer,  intent(in)  :: ngw, nkb, n
      real(DP), intent(in)  :: eigr(2,ngw,nat), c(2,ngw,n)
      real(DP), intent(out) :: becdr(nkb,n,3)
      !
      real(DP), allocatable :: gk(:)
      real(DP), allocatable :: wrk2(:,:,:)
      !
      integer   :: ig, is, iv, ia, k, l, ixr, ixi, inl, isa, i
      real(DP) :: signre, signim, arg
!
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )

      becdr = 0.d0
!
      do k = 1, 3

         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            allocate( wrk2( 2, ngw, na( is ) ) )

            do iv=1,nh(is)
               !
               !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
               !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
               l=nhtol(iv,is)
               if (l.eq.0) then
                  ixr = 2
                  ixi = 1
                  signre =  1.0d0
                  signim = -1.0d0
               else if (l.eq.1) then
                  ixr = 1
                  ixi = 2
                  signre = -1.0d0
                  signim = -1.0d0
               else if (l.eq.2) then
                  ixr = 2
                  ixi = 1
                  signre = -1.0d0
                  signim =  1.0d0
               else if (l == 3) then
                  ixr = 1
                  ixi = 2
                  signre =  1.0d0
                  signim =  1.0d0
               endif
!    
!$omp do
               do ia=1,na(is)
                  !    q = 0   component (with weight 1.0)
                  if (gstart == 2) then
                     wrk2(1,1,ia) = signre*gk(1)*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                     wrk2(2,1,ia) = signim*gk(1)*beta(1,iv,is)*eigr(ixi,1,ia+isa)
                  end if
                  !    q > 0   components (with weight 2.0)
                  do ig=gstart,ngw
                     arg = 2.0d0*gk(ig)*beta(ig,iv,is)
                     wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                     wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                  end do
               end do
!$omp end do
!$omp end parallel 
               inl=ish(is)+(iv-1)*na(is)+1
               CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becdr( inl, 1, k ), nkb )
            end do

            deallocate( wrk2 )

            isa = isa + na(is)

         end do

         IF( nproc_image > 1 ) THEN
            CALL mp_sum( becdr(:,:,k), intra_image_comm )
         END IF
      end do

      deallocate( gk )

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2_repl
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
   real(8) function ennl( rhovan, bec )
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      !
      implicit none
      !
      ! input
      !
      real(DP) :: bec( nkb, n )
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
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
                  do i = 1, n
                     iss = ispin(i)
                     sums(iss) = sums(iss) + f(i) * bec(inl,i) * bec(jnl,i)
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
      ennl = ennl_t
      !
      return
   end function ennl
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   subroutine calrhovan( rhovan, bec, iwf )
!-----------------------------------------------------------------------
      !
      ! calculation of rhovan relative to state iwf
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      !
      implicit none
      !
      ! input
      !
      real(DP) :: bec( nkb, n )
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
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
   end subroutine calrhovan
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
   subroutine calbec ( nspmn, nspmx, eigr, c, bec )
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
      use io_global,      only : stdout
      use cvan,           only : ish
      use electrons_base, only : n => nbsp
      use gvecw,          only : ngw
      use control_flags,  only : iprint, iprsta
      use uspp_param,     only : nh
      use uspp,           only : nkb
!
      implicit none
      !
      integer,     intent(in)  :: nspmn, nspmx
      real(DP),    intent(out) :: bec( nkb, n )
      complex(DP), intent(in)  :: c( ngw, n ), eigr( ngw,nat )

      ! local variables

      integer :: is, ia, i , iv
!
      call start_clock( 'calbec' )
      !
      call nlsm1( n, nspmn, nspmx, eigr, c, bec )
!
      if ( iprsta > 2 ) then
         WRITE( stdout,*)
         do is=1,nspmx
            WRITE( stdout,'(33x,a,i4)') ' calbec: bec (is)',is
            do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calbec: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &             ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
            end do
         end do
      endif
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE caldbec( ngw, nkb, n, nspmn, nspmx, eigr, c, dbec )
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
  use mp_global,  only : nproc_image, intra_image_comm
  use ions_base,  only : na, nat
  use cvan,       only : ish
  use cdvan,      only : dbeta
  use uspp,       only : nhtol
  use uspp_param, only : nh, nhm
  use reciprocal_vectors, only : gstart
  USE cp_main_variables,  ONLY : descla, la_proc, nlax, nlam
  USE descriptors,        ONLY : nlar_ , nlac_ , ilar_ , ilac_ , nlax_ , la_myr_ , la_myc_
  use electrons_base,     only : nspin, iupdwn, nupdwn
  !
  implicit none
  !
  integer,      intent(in)  :: ngw, nkb, n
  integer,      intent(in)  :: nspmn, nspmx
  complex(DP), intent(in)  :: c(ngw,n)
  real(DP),    intent(in)  :: eigr(2,ngw,nat)
  real(DP),    intent(out) :: dbec( nkb, 2*nlam, 3, 3 )
  !
  real(DP), allocatable :: wrk2(:,:,:), dwrk(:,:)
  !
  integer   :: ig, is, iv, ia, l, ixr, ixi, inl, i, j, ii, isa, nanh, iw, iss, nr, ir, istart, nss
  real(DP) :: signre, signim, arg
  !
  !
  !
  do j=1,3
     do i=1,3

        isa = 0
        do is = 1, nspmn - 1
          isa = isa + na(is)
        end do

        do is=nspmn,nspmx
           allocate( wrk2( 2, ngw, na(is) ) )
           nanh = na(is)*nh(is)
           allocate( dwrk( nanh, n ) )
           do iv=1,nh(is)
              l=nhtol(iv,is)
              if (l == 0) then
                 ixr = 1
                 ixi = 2
                 signre =  1.0d0
                 signim =  1.0d0
              else if (l == 1) then
                 ixr = 2
                 ixi = 1
                 signre =  1.0d0
                 signim = -1.0d0
              else if (l == 2) then
                 ixr = 1
                 ixi = 2
                 signre = -1.0d0
                 signim = -1.0d0
              else if (l == 3) then
                 ixr = 2
                 ixi = 1
                 signre = -1.0d0
                 signim =  1.0d0
              else
                 CALL errore(' caldbec  ', ' l not implemented ', ABS( l ) )
              endif
              !
              do ia=1,na(is)
                 if (gstart == 2) then
                    !     q = 0   component (with weight 1.0)
                    wrk2(1,1,ia)= signre*dbeta(1,iv,is,i,j)*eigr(ixr,1,ia+isa)
                    wrk2(2,1,ia)= signim*dbeta(1,iv,is,i,j)*eigr(ixi,1,ia+isa)
                 end if
                 !     q > 0   components (with weight 2.0)
                 do ig = gstart, ngw
                    arg = 2.0d0*dbeta(ig,iv,is,i,j)
                    wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                    wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                 end do
              end do
              inl=(iv-1)*na(is)+1
              CALL dgemm( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, dwrk(inl,1), nanh )
           end do
           deallocate( wrk2 )
           if( nproc_image > 1 ) then
              call mp_sum( dwrk, intra_image_comm )
           end if
           inl=ish(is)+1
           do iss=1,nspin
              IF( la_proc ) THEN
                 nr = descla( nlar_ , iss )
                 ir = descla( ilar_ , iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do ii = 1, nr
                    do iw = 1, nanh
                       dbec( iw + inl - 1, ii + (iss-1)*nlam, i, j ) = dwrk( iw, ii + ir - 1 + istart - 1 )
                       !dbec( iw + inl - 1, ii + (iss-1)*nspin, i, j ) = dwrk( iw, ii + ir - 1 + istart - 1 )
                    end do
                 end do
              END IF
           end do
           deallocate( dwrk )
           isa = isa + na(is)
        end do
     end do
  end do

  !
  return
end subroutine caldbec
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine dennl( bec, dbec, drhovan, denl )
  !-----------------------------------------------------------------------
  !
  !  compute the contribution of the non local part of the
  !  pseudopotentials to the derivative of E with respect to h
  !
  USE kinds,      ONLY : DP
  use cvan,       only : ish
  use uspp_param, only : nh, nhm
  use uspp,       only : nkb, dvan, deeq
  use ions_base,  only : nsp, na, nat
  use cell_base,  only : h
  use io_global,  only : stdout
  use mp,         only : mp_sum
  use mp_global,  only : intra_image_comm
  USE cp_main_variables,  ONLY : descla, la_proc, nlax, nlam
  USE descriptors,        ONLY : nlar_ , nlac_ , ilar_ , ilac_ , nlax_ , la_myr_ , la_myc_
  use electrons_base,     only : n => nbsp, ispin, f, nspin, iupdwn, nupdwn
  use reciprocal_vectors, only : gstart

  implicit none

  real(DP), intent(in)  :: dbec( nkb, 2*nlam, 3, 3 )
  real(DP), intent(in)  :: bec( nkb, n )
  real(DP), intent(out) :: drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 )
  real(DP), intent(out) :: denl( 3, 3 )

  real(DP) :: dsum(3,3),dsums(2,3,3), detmp(3,3)
  integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i,j,k
  integer   :: istart, nss, ii, ir, nr
  !
  denl=0.d0
  drhovan=0.0d0

  IF( la_proc ) THEN


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
                 IF( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) THEN
                 nr = descla( nlar_ , iss )
                 ir = descla( ilar_ , iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do i=1,nr
                    ii = i+istart-1+ir-1
                    do k=1,3
                       do j=1,3
                          dsums(iss,k,j)=dsums(iss,k,j)+f(ii)*       &
 &                          (dbec(inl,i+(iss-1)*nlam,k,j)*bec(jnl,ii)          &
 &                          + bec(inl,ii)*dbec(jnl,i+(iss-1)*nlam,k,j))
                       enddo
                    enddo
                 end do
                 END IF
              end do
              !
              do iss=1,nspin
                 IF( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) THEN
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

  END IF

  CALL mp_sum( denl,    intra_image_comm )
  do k=1,3
     do j=1,3
        CALL mp_sum( drhovan(:,:,:,j,k), intra_image_comm )
     end do
  end do

!  WRITE(6,*) 'DEBUG enl (CP) = '
!  detmp = denl
!  detmp = MATMUL( detmp(:,:), TRANSPOSE( h ) )
!  WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

  !
  return
end subroutine dennl
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine nlfq( c, eigr, bec, becdr, fion )
  !-----------------------------------------------------------------------
  !
  !     contribution to fion due to nonlocal part
  !
  USE kinds,          ONLY : DP
  use uspp,           only : nkb, dvan, deeq
  use uspp_param,     only : nhm, nh
  use cvan,           only : ish, nvb
  use ions_base,      only : nax, nat, nsp, na
  use electrons_base, only : n => nbsp, ispin, f, nspin, iupdwn, nupdwn
  use gvecw,          only : ngw
  use constants,      only : pi, fpi
  use mp_global,      only : me_image, intra_image_comm, nproc_image
  use mp,             only : mp_sum
  USE cp_main_variables, ONLY: nlax, descla, la_proc
  USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , &
                               la_myr_ , la_myc_
  !
  implicit none
  !
  real(DP),    intent(in)  :: bec( nkb, n ), c( 2, ngw, n )
  real(DP),    intent(out) :: becdr( nkb, nspin*nlax, 3 )
  complex(DP), intent(in)  :: eigr( ngw, nat )
  real(DP),    intent(out) :: fion( 3, nat )
  !
  integer   :: k, is, ia, isa, iss, inl, iv, jv, i, ir, nr, nss, istart, ioff
  real(DP) :: temp
  !
  real(DP), allocatable :: tmpbec(:,:), tmpdr(:,:) 
  real(DP), allocatable :: fion_loc(:,:)
#ifdef __OPENMP 
  INTEGER :: mytid, ntids, omp_get_thread_num, omp_get_num_threads
#endif  
  !
  call start_clock( 'nlfq' )
  !
  !
  !     nlsm2 fills becdr
  !
  call nlsm2( ngw, nkb, n, nspin, eigr, c, becdr )
  !
  allocate ( fion_loc( 3, nat ) )
  !
  fion_loc = 0.0d0

  !
  DO k = 1, 3

!$omp parallel default(shared), &
!$omp private(tmpbec,tmpdr,isa,is,ia,iss,nss,istart,ir,nr,ioff,iv,jv,inl,temp,i,mytid,ntids)

#ifdef __OPENMP
     mytid = omp_get_thread_num()  ! take the thread ID
     ntids = omp_get_num_threads() ! take the number of threads
#endif

     allocate ( tmpbec( nhm, nlax ), tmpdr( nhm, nlax ) )

     isa = 0
     !
     DO is=1,nsp
        DO ia=1,na(is)

           isa=isa+1

#ifdef __OPENMP
           ! distribute atoms round robin to threads
           !
           IF( MOD( isa, ntids ) /= mytid ) CYCLE
#endif  
           DO iss = 1, nspin

              nss = nupdwn( iss )
              istart = iupdwn( iss )

              IF( la_proc .AND. &
                  ( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) ) THEN

                 ! only processors on the diagonal of the square proc grid enter here.
                 ! This is to distribute the load among different multi-core nodes,
                 ! and maximize the memory bandwith per core.

                 tmpbec = 0.d0
                 tmpdr  = 0.d0

                 ir = descla( ilar_ , iss )
                 nr = descla( nlar_ , iss )

                 ioff = istart-1+ir-1

                 do iv=1,nh(is)
                    do jv=1,nh(is)
                       inl=ish(is)+(jv-1)*na(is)+ia
                       temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
                       do i=1,nr
                          tmpbec(iv,i)=tmpbec(iv,i)+temp*bec(inl,i+ioff)
                       end do
                    end do
                 end do

                 do iv=1,nh(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    do i=1,nr
                       tmpdr(iv,i)=f(i+ioff)*becdr( inl, i+(iss-1)*nlax, k )
                    end do
                 end do

                 do i=1,nr
                    do iv=1,nh(is)
                       tmpdr(iv,i)=tmpdr(iv,i)*tmpbec(iv,i)
                    end do
                 end do

                 fion_loc(k,isa) = fion_loc(k,isa)-2.d0*SUM(tmpdr)

              END IF
           END DO
        END DO
     END DO
     deallocate ( tmpbec, tmpdr )
!$omp end parallel
  END DO
  !
  CALL mp_sum( fion_loc, intra_image_comm )
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
end subroutine nlfq

