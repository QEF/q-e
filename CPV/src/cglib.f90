!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
   subroutine calcmt( nrlx, fdiag, zmat, fmat )
!-----------------------------------------------------------------------
!
!  constructs fmat=z0^t.fdiag.z0    zmat = z0^t
!
      USE kinds,             ONLY: DP
      use electrons_base,    ONLY: nudx, nspin, nupdwn, iupdwn, nx => nbspx
      USE cp_main_variables, ONLY: descla
      USE descriptors,       ONLY: la_descriptor , ldim_cyclic
      USE mp,                ONLY: mp_sum, mp_bcast

      implicit none

      integer, intent(in)  :: nrlx
      real(DP) :: zmat( nrlx, nudx, nspin ), fmat( nrlx, nudx, nspin )
      !  NOTE: zmat and fmat are distributed by row across processors
      real(dp) :: fdiag( nx ) !  fdiag is replicated

      integer  :: iss, nss, istart, i, j, k, ii, jj, kk
      integer  :: np_rot, me_rot, nrl, comm_rot, ip, nrl_ip

      real(DP), ALLOCATABLE :: mtmp(:,:)
      real(DP) :: f_z0t


      call start_clock('calcmt')

      fmat = 0.0d0

      DO iss = 1, nspin

         nss      = nupdwn( iss )
         istart   = iupdwn( iss )
         np_rot   = descla( iss )%npr * descla( iss )%npc
         me_rot   = descla( iss )%mype
         nrl      = descla( iss )%nrl
         comm_rot = descla( iss )%comm

         IF( descla( iss )%active_node > 0 ) THEN

            ALLOCATE( mtmp( MAXVAL(descla(:)%nrlx), nudx ) )

            DO ip = 1, np_rot

               IF( me_rot == ( ip - 1 ) ) THEN
                  mtmp = zmat(:,:,iss)
               END IF
               nrl_ip = ldim_cyclic( nss, np_rot, ip - 1 )
               CALL mp_bcast( mtmp , ip - 1 , comm_rot )

               DO j = 1, nss
                  ii = ip
                  DO i = 1, nrl_ip
                     f_z0t = fdiag( j + istart - 1 ) * mtmp( i, j )
                     DO k = 1, nrl
                        fmat( k, ii, iss ) = fmat( k, ii, iss )+ zmat( k, j, iss ) * f_z0t 
                     END DO
                     ii = ii + np_rot
                  END DO
               END DO

            END DO

            DEALLOCATE( mtmp )

         END IF

      END DO

      call stop_clock('calcmt')

      RETURN
      END SUBROUTINE calcmt


!-----------------------------------------------------------------------
      subroutine rotate( nrlx, z0, c0, bec, c0diag, becdiag )
!-----------------------------------------------------------------------
      use kinds, only: dp
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx, n => nbsp
      use uspp_param, only: nh, ish, nvb
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use gvecw, only: ngw
      use ions_base, only: nsp, na
      USE cp_main_variables, ONLY: descla
      USE descriptors,       ONLY: la_descriptor
      USE cp_interfaces,     ONLY: protate

      implicit none
      integer, intent(in) :: nrlx
      real(kind=DP)    z0( nrlx, nudx, nspin )
      real(kind=DP)    bec( nhsa, n ), becdiag( nhsa, n )
      complex(kind=DP) c0( ngw, nx ), c0diag( ngw, nx )
      integer :: np_rot, me_rot, nrl, comm_rot
      integer iss, nss, istart
  
      CALL start_clock( 'rotate' )

      DO iss = 1, nspin
         istart   = iupdwn( iss )
         nss      = nupdwn( iss )
         np_rot   = descla( iss )%npr * descla( iss )%npc
         me_rot   = descla( iss )%mype
         nrl      = descla( iss )%nrl
         comm_rot = descla( iss )%comm
         CALL protate ( c0, bec, c0diag, becdiag, ngw, nss, istart, z0(:,:,iss), nrl, &
                        na, nsp, ish, nh, np_rot, me_rot, comm_rot )
      END DO

      CALL stop_clock( 'rotate' )
      return
      end subroutine rotate


!-----------------------------------------------------------------------
      subroutine ddiag(nx,n,amat,dval,dvec,iflag)
!-----------------------------------------------------------------------
!
      use dspev_module, only: dspev_drv
      use kinds , only : dp

      implicit none

      integer nx,n,ndim,iflag,k,i,j
      real(dp)   dval(n)
      real(dp) amat(nx,n), dvec(nx,n)
      real(dp), allocatable::  ap(:)

      ndim=(n*(n+1))/2
      allocate(ap(ndim))
      ap(:)=0.d0

      k=0
      do j=1,n
       do i=1,j
        k=k+1
        ap(k)=amat(i,j)
       end do
      end do

      CALL dspev_drv( 'V', 'U', n, ap, dval, dvec, nx )

      deallocate(ap)

      return
    end subroutine ddiag

    subroutine minparabola(ene0,dene0,ene1,passop,passo,stima)
!this subroutines finds the minimum of a quadratic real function
      
      use kinds, only : dp

      implicit none
      real(dp) ene0,dene0,ene1,passop,passo,stima
      real(dp) a,b,c!a*x^2+b*x+c
      
      c=ene0
      b=dene0
      a=(ene1-b*passop-c)/(passop**2.d0)
      
      passo = -b/(2.d0*a)
      if( a.lt.0.d0) then
         if(ene1.lt.ene0) then
            passo=passop
         else 
            passo=0.5d0*passop
         endif
      endif


      stima=a*passo**2.d0+b*passo+c


      return
    end subroutine minparabola


subroutine pc2(a,beca,b,becb)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
    
      use kinds, only: dp 
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin,  nupdwn, iupdwn, nspin
      use uspp_param, only: nh, nvb, ish
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      
                           
      implicit none        
                           
      complex(kind=DP) a(ngw,n), b(ngw,n)
                     
      real(kind=DP)    beca(nhsa,n),becb(nhsa,n)
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(kind=DP) sca
      real(DP), allocatable :: bectmp(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc2' )

      do iss= 1, nspin
         nss= nupdwn( iss )
         istart= iupdwn( iss )

         allocate(bectmp(nss,nss))
         allocate(zbectmp(nss,nss))
         bectmp(:,:)=0.d0

         call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart),ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

         do j=1,nss
            do i=1,nss
               bectmp(i,j)=2.d0*dble(zbectmp(i,j))
               if(gstart==2) bectmp(i,j)=bectmp(i,j)-DBLE(a(1,j))*DBLE(b(1,i))
               
            enddo
         enddo
         deallocate(zbectmp)
         call mp_sum( bectmp(:,:), intra_bgrp_comm)
         if(nvb >= 0) then

            nl_max=0
            do is=1,nvb
               nl_max=nl_max+nh(is)*na(is)
            enddo
            allocate (qq_tmp(nl_max,nl_max))
            allocate (qqb_tmp(nl_max,nss))
            qq_tmp(:,:)=0.d0
            do is=1,nvb
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        qq_tmp(inl,jnl)=qq(iv,jv,is)
                     enddo
                  enddo
               enddo
            enddo
            if(.not. mat_par)  then
               call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
               call dgemm('T','N',nss,nss,nl_max,1.d0,beca(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
            else
               call para_dgemm &
& ('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_bgrp_comm)
               call para_dgemm &
&('T','N',nss,nss,nl_max,1.d0,beca(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss, intra_bgrp_comm)
            endif
            deallocate(qq_tmp,qqb_tmp)
         endif
         allocate(zbectmp(nss,nss))
         do i=1,nss
            do j=1,nss
               zbectmp(i,j)=CMPLX(bectmp(i,j),0.d0,kind=dp)
            enddo
         enddo
         call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),a(:,istart),ngw,zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
         deallocate(zbectmp)
         call dgemm('N','N',nhsa,nss,nss,1.0d0,beca(:,istart),nhsa,bectmp,nss,1.0d0,becb(:,istart),nhsa)
         deallocate(bectmp)
      enddo!on spin
      CALL stop_clock( 'pc2' )
      return
    end subroutine pc2


    subroutine pcdaga2(a,as ,b )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use kinds
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin
      use uspp_param, only: nh, ish, nvb
      use uspp, only :nhsa=>nkb


      implicit none

      complex(dp) a(ngw,n), b(ngw,n), as(ngw,n)
      ! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(dp) sca
      real(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate(scar(n))
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               if (gstart==2) b(1,i) = CMPLX(dble(b(1,i)),0.0d0,kind=dp)
               do  ig=1,ngw           !loop on g vectors
                  sca=sca+DBLE(CONJG(a(ig,j))*b(ig,i))
               enddo
               sca = sca*2.0d0  !2. for real weavefunctions
               if (gstart==2) sca = sca - dble(a(1,j))*dble(b(1,i))
            endif
            scar(i) = sca
         enddo
         
          
         call mp_sum( scar, intra_bgrp_comm )


         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca*as(ig,j)
               enddo
               ! this to prevent numerical errors
               if (gstart==2) b(1,i) = CMPLX(dble(b(1,i)),0.0d0,kind=dp)
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2')
      return
      end subroutine pcdaga2

     subroutine set_x_minus1(betae,m_minus1,ema0bg,use_ema)

! this function calculates the factors for the inverse of the US K  matrix
! it takes care of the preconditioning

      use kinds, only: dp
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum, mp_bcast
      use electrons_base, only: n => nbsp, ispin
      use uspp_param, only: nh, ish, nvb
      use uspp, only :nhsa=>nkb,qq,nhsavb=>nkbus
      use io_global, ONLY: ionode, ionode_id

      implicit none

      complex(DP) :: betae(ngw,nhsa)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      real(DP)    :: ema0bg(ngw)
      logical     :: use_ema


! local variables
      real(DP),allocatable :: q_matrix(:,:), b_matrix(:,:),c_matrix(:,:)
      integer is, iv, jv, ia, inl, jnl, i, j, k,ig, js, ja
      real(DP) sca
      integer info, lwork
      integer, allocatable :: ipiv(:)
      real(dp),allocatable :: work(:)

      call start_clock('set_x_minus1')
      allocate(ipiv(nhsavb))
      allocate(work(nhsavb))

      lwork=nhsavb

      allocate(q_matrix(nhsavb,nhsavb),c_matrix(nhsavb,nhsavb))
!construct q matrix
      q_matrix(:,:) = 0.d0

      do is=1,nvb
         do iv=1,nh(is)
            do jv=1,nh(is)
               do ia=1,na(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    jnl=ish(is)+(jv-1)*na(is)+ia
                    q_matrix(inl,jnl)= qq(iv,jv,is)
               enddo
            enddo
         enddo
      enddo

!construct b matrix
! m_minus1 used to be b matrix
      m_minus1(:,:) = 0.d0
      do is=1,nvb
         do ia=1,na(is)
            do iv=1,nh(is)
               do js=1,nvb
                  do ja=1,na(js)
                     do jv=1,nh(js)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(js)+(jv-1)*na(js)+ja
                        sca=0.d0
                        if (use_ema) then
                           ! k_minus case
                        do  ig=1,ngw           !loop on g vectors
                           sca=sca+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
                        enddo
                        sca = sca*2.0d0  !2. for real weavefunctions
                        if (gstart==2) sca = sca - ema0bg(1)*DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                        else
                           ! s_minus case
                        do  ig=1,ngw           !loop on g vectors
                           sca=sca+DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
                        enddo
                        sca = sca*2.0d0  !2. for real weavefunctions
                        if (gstart==2) sca = sca - DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                        endif
                        m_minus1(inl,jnl)=sca
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      call mp_sum( m_minus1, intra_bgrp_comm )

!calculate -(1+QB)**(-1) * Q
      CALL dgemm('N','N',nhsavb,nhsavb,nhsavb,1.0d0,q_matrix,nhsavb,m_minus1,nhsavb,0.0d0,c_matrix,nhsavb)

      do i=1,nhsavb
         c_matrix(i,i)=c_matrix(i,i)+1.d0
      enddo

      if(ionode) then
        call dgetrf(nhsavb,nhsavb,c_matrix,nhsavb,ipiv,info)
        if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetrf :', info
        call dgetri(nhsavb,c_matrix,nhsavb,ipiv,work,lwork,info)
        if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetri :', info
      endif
      call mp_bcast( c_matrix, ionode_id, intra_bgrp_comm )


      CALL dgemm('N','N',nhsavb,nhsavb,nhsavb,-1.0d0,c_matrix,nhsavb,q_matrix,nhsavb,0.0d0,m_minus1,nhsavb)

      deallocate(q_matrix,c_matrix)
      deallocate(ipiv,work)
      call stop_clock('set_x_minus1')
      return
    end subroutine set_x_minus1
!
      subroutine xminus1(c0,betae,ema0bg,beck,m_minus1,do_k)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0> 
! endif
      use kinds, only: dp
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use uspp_param, only: nh, nvb, ish
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use gvect, only: gstart
!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp)    beck(nhsa,n), ema0bg(ngw)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      logical :: do_k
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, jv, ia, inl, jnl, i, j, js, ja,ig
      real(dp) becktmp

      
      logical :: mat_par=.true.!if true uses parallel routines      

      call start_clock('xminus1')
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0

            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (gstart==2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i)) 
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               enddo
            enddo
            call mp_sum( beck, intra_bgrp_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nhsavb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      if(.not.mat_par) then
         call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb )
      else
         call para_dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb,intra_bgrp_comm )
      endif

!NB  nhsavb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones

         CALL dgemm( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
                    qtemp, nhsavb, 0.0d0, phi, 2*ngw )
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
         else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
         endif
      deallocate(qtemp,phi)

      else
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=c0(ig,j)*ema0bg(ig)
               end do
            end do
         endif
      endif
      call stop_clock('xminus1')
      return
     end subroutine xminus1

      SUBROUTINE emass_precond_tpa( ema0bg, tpiba2, emaec )
       use kinds, ONLY : dp
       use gvecw, ONLY : g2kin,ngw
       IMPLICIT NONE
       REAL(DP), INTENT(OUT) :: ema0bg(ngw)
       REAL(DP), INTENT(IN) ::  tpiba2, emaec
       INTEGER :: i

       real(DP) :: x

       call start_clock('emass_p_tpa')
       do i = 1, ngw

          x=0.5d0*tpiba2*g2kin(i)/emaec
          ema0bg(i) = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
       end do
       call stop_clock('emass_p_tpa')
      RETURN
      END SUBROUTINE emass_precond_tpa

      subroutine ave_kin( c, ngwx, n, ene_ave )
!this subroutine calculates the average kinetic energy of
!each state , to be used for preconditioning


      USE kinds,              ONLY: DP
      USE constants,          ONLY: pi, fpi
      USE gvecw,              ONLY: ngw
      USE gvect, ONLY: gstart
      USE gvecw,              ONLY: g2kin
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE cell_base,          ONLY: tpiba2
                                                                                                                             
      IMPLICIT NONE
                                                                                                                            
                                                                                                                             
      ! input
                                                                                                                             
      INTEGER,     INTENT(IN) :: ngwx, n
      COMPLEX(kind=DP), INTENT(IN) :: c( ngwx, n )
      REAL(kind=DP), INTENT(out) :: ene_ave(n)!average kinetic energy to be calculated
      !
      ! local
                                                                                                                             
      INTEGER  :: ig, i

      !
      DO i=1,n
         ene_ave(i)=0.d0
         DO ig=gstart,ngw
            ene_ave(i)=ene_ave(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*g2kin(ig)
         END DO
      END DO


      CALL mp_sum( ene_ave(1:n), intra_bgrp_comm )
      ene_ave(:)=ene_ave(:)*tpiba2
                                                                                                                             
      RETURN
    END subroutine ave_kin



      subroutine xminus1_state(c0,betae,ema0bg,beck,m_minus1,do_k,ave_kin)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0>
! endif
!adapted for state by state
      use kinds, only: dp
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use uspp_param, only: nh, nvb, ish
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use gvect, only: gstart
      USE gvecw,              ONLY: g2kin
      USE cell_base,          ONLY: tpiba2


!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp)    beck(nhsa,n), ema0bg(ngw)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      logical :: do_k
      real(kind=DP) :: ave_kin(n)!average kinetic energy per state
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, jv, ia, inl, jnl, i, j, js, ja,ig
      real(dp) becktmp
      real(kind=DP) :: prec_fact, x

 
      call start_clock('xminus1')
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0
 
            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (gstart==2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i))
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               enddo
            enddo
            call mp_sum( beck, intra_bgrp_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nhsavb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb )
 
 
 
!NB nhsavb is the total number of US projectors, it works because the first pseudos are the vanderbilt's ones

      CALL dgemm( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
           qtemp, nhsavb, 0.0d0, phi, 2*ngw )
      if (do_k) then
         do j=1,n
            do ig=1,ngw
               x=tpiba2*g2kin(i)/ave_kin(j)
               prec_fact = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
                c0(ig,j)=c0(ig,j)*prec_fact
               !c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
            end do
         end do
      else
         do j=1,n
            do i=1,ngw
               c0(i,j)=(phi(i,j)+c0(i,j))
            end do
         end do
      endif
      deallocate(qtemp,phi)
      
   else
      if (do_k) then
         do j=1,n
            do ig=1,ngw
               x=tpiba2*g2kin(ig)/ave_kin(j)
               prec_fact = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
               c0(ig,j)=c0(ig,j)*prec_fact
            end do
         end do
      endif
   endif
   call stop_clock('xminus1')
   return
 end subroutine xminus1_state
!
! ... some simple routines for parallel linear algebra (the matrices are
! ... always replicated on all the cpus)
!
! ... written by carlo sbraccia ( 2006 )
!
!----------------------------------------------------------------------------
SUBROUTINE para_dgemm( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix B by columns) of dgemm 
  !
  USE kinds, ONLY : DP
  USE parallel_toolkit
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
  INTEGER,          INTENT(IN)    :: m, n, k
  REAL(DP),         INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, ldb, ldc
  REAL(DP),         INTENT(INOUT) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER,          INTENT(IN)    :: comm
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( ( alpha == 0.0_DP .OR. k == 0 ) .AND. beta == 1.0_DP ) ) RETURN
  !
!write(*,*) 'DEBUG: para_dgemm'
  !
  CALL rep_matmul_drv( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  RETURN
  !
END SUBROUTINE para_dgemm
