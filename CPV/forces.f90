!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


    SUBROUTINE dforce1_x( co, ce, dco, dce, fio, fie, hg, v, psi_stored )

      USE kinds,      ONLY: DP
      USE fft_base,   ONLY: dffts
      USE gvecw,      ONLY: ngw
      USE cp_interfaces, ONLY: fwfft, invfft
      USE cell_base,  ONLY: tpiba2

      IMPLICIT NONE

      ! ... declare subroutine arguments
      COMPLEX(DP), INTENT(OUT) :: dco(:), dce(:)
      COMPLEX(DP), INTENT(IN)  :: co(:), ce(:)
      REAL(DP),    INTENT(IN)  :: fio, fie
      REAL(DP),    INTENT(IN)  :: v(:)
      REAL(DP),    INTENT(IN)  :: hg(:)
      COMPLEX(DP), OPTIONAL    :: psi_stored(:)

      ! ... declare other variables
      !
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP) :: fp, fm, aro, are
      REAL(DP)    :: fioby2, fieby2, arg
      INTEGER      :: ig

      !  end of declarations

      ALLOCATE( psi( SIZE(v) ) )

      IF( PRESENT( psi_stored ) ) THEN
        psi = psi_stored * CMPLX(v, 0.0d0)
      ELSE
        CALL c2psi( psi, dffts%nnr, co, ce, ngw, 2 )
        CALL invfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )
        psi = psi * CMPLX(v, 0.0d0)
      END IF

      CALL fwfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )
      CALL psi2c( psi, dffts%nnr, dco, dce, ngw, 2 )

      DEALLOCATE(psi)

      fioby2   = fio * 0.5
      fieby2   = fie * 0.5

      DO ig = 1, SIZE(co)
        fp = dco(ig) + dce(ig)
        fm = dco(ig) - dce(ig)
        aro = CMPLX(  DBLE(fp), AIMAG(fm) )
        are = CMPLX( AIMAG(fp), -DBLE(fm))
        arg = tpiba2 * hg(ig)
        dco(ig) = -fioby2 * (arg * co(ig) + aro)
        dce(ig) = -fieby2 * (arg * ce(ig) + are)
      END DO

    RETURN
    END SUBROUTINE dforce1_x


!=----------------------------------------------------------------------------=!


    SUBROUTINE dforce2_x( fio, fie, df, da, vkb, beco, bece )

        !  this routine computes:
        !  the generalized force df=CMPLX(dfr,dfi) acting on the i-th
        !  electron state at the ik-th point of the Brillouin zone
        !  represented by the vector c=CMPLX(cr,ci)
        !  ----------------------------------------------

      USE kinds,                   ONLY: DP
      USE ions_base,               ONLY: na
      USE read_pseudo_module_fpmd, ONLY: nspnl
      USE uspp_param,              ONLY: nh
      USE uspp,                    ONLY: nhtol, nhtolm, indv, beta, dvan, nkb
      use cvan,                    ONLY: ish
      !

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(DP), INTENT(IN) :: vkb(:,:)
      REAL(DP), INTENT(IN) :: fio, fie
      COMPLEX(DP)  :: df(:), da(:)
      REAL(DP), INTENT(IN) :: beco(:)
      REAL(DP), INTENT(IN) :: bece(:)

! ... declare other variables
      REAL(DP) :: to, te
      INTEGER  :: l, is, ig, ngw, iv, inl, isa

      ! ----------------------------------------------

      ngw  = SIZE(df)

      isa = 1
      
      DO is = 1, nspnl
        !
        DO iv = 1, nh( is )
          !
          inl = ish(is) + (iv-1) * na(is) + 1

          to = - fio * dvan( iv, iv, is )
          !
          te = - fie * dvan( iv, iv, is )

          CALL DGEMV('N', 2*ngw, na(is), to, vkb( 1, inl ), &
               2*SIZE(vkb,1), beco( inl ), 1, 1.0d0, df, 1)
          !
          CALL DGEMV('N', 2*ngw, na(is), te, vkb( 1, inl ), &
               2*SIZE(vkb,1), bece( inl ), 1, 1.0d0, da, 1)
          !
        END DO
        !
        isa = isa + na( is )
        !
      END DO
      !

      RETURN
    END SUBROUTINE dforce2_x



!=----------------------------------------------------------------------------=!

     

    SUBROUTINE dforce_fpmd( ib, c, f, df, da, v, vkb, bec, n, noffset )
       !
       USE kinds,              ONLY: DP
       USE reciprocal_vectors, ONLY: ggp, g, gx
       USE cp_interfaces
       !
       IMPLICIT NONE
       !
       INTEGER,     INTENT(IN)  :: ib     ! band index
       COMPLEX(DP), INTENT(IN)  :: c(:,:)
       COMPLEX(DP), INTENT(OUT) :: df(:), da(:)
       REAL (DP),   INTENT(IN)  :: v(:), bec(:,:), f(:)
       COMPLEX(DP), INTENT(IN)  :: vkb(:,:)
       INTEGER,     INTENT(IN)  :: n, noffset  ! number of bands, and band index offset
       !
       COMPLEX(DP), ALLOCATABLE :: dum( : )   
       !
       INTEGER :: in
       !
       IF( ib > n ) CALL errore( ' dforce ', ' ib out of range ', 1 )
       !
       in = noffset + ib - 1 
       !
       IF( ib == n ) THEN
          !
          ALLOCATE( dum( SIZE( df ) ) )
          !
          CALL dforce1( c( :, in ), c( :, in ), df, dum, f( in ), f( in ), ggp, v )
          !
          CALL dforce2( f( in ), f( in ), df , dum , vkb, bec( :, in ), bec( :, in ) )
          !
          DEALLOCATE( dum )
          !
       ELSE
          !
          CALL dforce1( c( :, in ), c( :, in+1 ), df, da, f( in ), f(in+1), ggp, v )
          !
          CALL dforce2( f(in), f(in+1), df, da, vkb, bec( :, in ), bec( :, in+1 ) )
          !
       END IF
       !
       RETURN
    END SUBROUTINE dforce_fpmd


!  ----------------------------------------------
  

    SUBROUTINE dforce_all( c, f, cgrad, vpot, vkb, bec, n, noffset )
       !
       USE kinds,              ONLY: DP
       USE cp_interfaces

       IMPLICIT NONE

       COMPLEX(DP),           INTENT(INOUT) :: c(:,:)
       REAL(DP),              INTENT(IN)    :: vpot(:), f(:)
       COMPLEX(DP),           INTENT(OUT)   :: cgrad(:,:)
       COMPLEX(DP),           INTENT(IN)    :: vkb(:,:)
       REAL(DP),              INTENT(IN)    :: bec(:,:)
       INTEGER,               INTENT(IN)    :: n, noffset
       
       INTEGER :: ib, in
       !
       IF( n > 0 ) THEN
          !
          !   Process two states at the same time
          !
          DO ib = 1, n-1, 2
             !
             in = ib + noffset - 1
             !
             CALL dforce( ib, c, f, cgrad(:,in), cgrad(:,in+1), vpot, vkb, bec, n, noffset )
             !
          END DO
          !
          !   and now process the last state in case that n is odd
          !
          IF( MOD( n, 2 ) /= 0 ) THEN
             !
             in = n + noffset - 1
             !
             CALL dforce( n, c, f, cgrad(:,in), cgrad(:,in), vpot, vkb, bec, n, noffset )
             !
          END IF
          !
       END IF
       !
       RETURN
    END SUBROUTINE dforce_all



!
!-------------------------------------------------------------------------
      SUBROUTINE dforce ( bec, betae, i, c, ca, df, da, v, ispin, f, n, nspin )
!-----------------------------------------------------------------------
!computes: the generalized force df=CMPLX(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=CMPLX(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
      USE kinds, ONLY: dp
      USE control_flags, ONLY: iprint
      USE gvecs
      USE gvecw, ONLY: ngw
      USE cvan, ONLY: ish
      USE uspp, ONLY: nhsa=>nkb, dvan, deeq
      USE uspp_param, ONLY: nhm, nh
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      USE constants, ONLY: pi, fpi
      USE ions_base, ONLY: nsp, na, nat
      USE gvecw, ONLY: ggp
      USE cell_base, ONLY: tpiba2
      USE ensemble_dft, ONLY: tens
      USE funct, ONLY: dft_is_meta
      USE cp_interfaces, ONLY: fwfft, invfft
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: n, nspin
      INTEGER, INTENT(IN) :: ispin( n )
      REAL(DP), INTENT(IN) :: f( n )
      !
      COMPLEX(DP) betae(ngw,nhsa), c(ngw), ca(ngw), df(ngw), da(ngw)
      REAL(DP) bec(nhsa,n), v(nnrsx,nspin)
      !
      ! local variables
      !
      INTEGER iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      REAL(DP) fi, fip, dd
      COMPLEX(DP) fp,fm,ci
      REAL(DP) af(nhsa), aa(nhsa) ! automatic arrays
      COMPLEX(DP), ALLOCATABLE :: dtemp(:)
      COMPLEX(DP), ALLOCATABLE :: psi(:)
!
!
      CALL start_clock( 'dforce' ) 
      !
      ALLOCATE( psi( nnrsx ) )
      ALLOCATE( dtemp( ngw ) )
!
!     important: if n is odd => c(*,n+1)=0.
! 
      IF (MOD(n,2).NE.0.AND.i.EQ.n) THEN
         DO ig=1,ngw
            ca(ig)=(0.,0.)
         END DO
      ENDIF
!
      ci=(0.0,1.0)
!
      psi (:) = (0.d0, 0.d0)
      DO ig=1,ngw
            psi(nms(ig))=CONJG(c(ig)-ci*ca(ig))
            psi(nps(ig))=c(ig)+ci*ca(ig)
      END DO

      CALL invfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!     
      iss1=ispin(i)
!
! the following avoids a potential out-of-bounds error
!
      IF (i.NE.n) THEN
         iss2=ispin(i+1)
      ELSE
         iss2=iss1
      END IF
!
      DO ir=1,nnrsx
         psi(ir)=CMPLX(v(ir,iss1)* DBLE(psi(ir)), v(ir,iss2)*AIMAG(psi(ir)) )
      END DO
!
      CALL fwfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
!     note : the factor 0.5 appears 
!       in the kinetic energy because it is defined as 0.5*g**2
!       in the potential part because of the logics
!
   
     IF (tens) THEN
        fi =-0.5
        fip=-0.5
      ELSE
        fi =-  f(i)*0.5
        fip=-f(i+1)*0.5
      END IF

      DO ig=1,ngw
         fp= psi(nps(ig)) + psi(nms(ig))
         fm= psi(nps(ig)) - psi(nms(ig))
         df(ig)= fi*(tpiba2*ggp(ig)* c(ig)+CMPLX(DBLE(fp), AIMAG(fm)))
         da(ig)=fip*(tpiba2*ggp(ig)*ca(ig)+CMPLX(AIMAG(fp),-DBLE(fm)))
      END DO

      IF(dft_is_meta()) CALL dforce_meta(c,ca,df,da,psi,iss1,iss2,fi,fip) !METAGGA
!
!     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
! 
      IF(nhsa.GT.0)THEN
         DO inl=1,nhsa
            af(inl)=0.
            aa(inl)=0.
         END DO
!
         DO is=1,nsp
            DO iv=1,nh(is)
               DO jv=1,nh(is)
                  isa=0
                  DO ism=1,is-1
                     isa=isa+na(ism)
                  END DO
                  DO ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     isa=isa+1
                     dd = deeq(iv,jv,isa,iss1)+dvan(iv,jv,is)
                     IF(tens) THEN
                      af(inl)=af(inl)-dd*bec(jnl,  i)
                     ELSE
                      af(inl)=af(inl)- f(i)*dd*bec(jnl,  i)
                     END IF
                     dd = deeq(iv,jv,isa,iss2)+dvan(iv,jv,is)
                     IF(tens) THEN
                      IF (i.NE.n) aa(inl)=aa(inl)-dd*bec(jnl,i+1)
                     ELSE
                      IF (i.NE.n) aa(inl)=aa(inl)-f(i+1)*dd*bec(jnl,i+1)
                     END IF
                  END DO
               END DO
            END DO
         END DO
!
         dtemp = 0.0d0
         CALL MXMA                                                      &
     &        (betae,1,2*ngw,af,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         DO ig=1,ngw
            df(ig)=df(ig)+dtemp(ig)
         END DO
!
         dtemp = 0.0d0
         CALL MXMA                                                      &
     &        (betae,1,2*ngw,aa,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         DO ig=1,ngw
            da(ig)=da(ig)+dtemp(ig)
         END DO
      ENDIF

      DEALLOCATE( dtemp )
      DEALLOCATE( psi )
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
      END SUBROUTINE dforce
!
!=========================================================================
!C. Bekas, IBM Research Zurich.
! dforce with Task Groups parallelization
!=========================================================================
!-------------------------------------------------------------------------
      subroutine dforce_bgl (bec,betae,i,c,ca,df,da,v)
!-----------------------------------------------------------------------
!computes: the generalized force df=CMPLX(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=CMPLX(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
      use kinds, only: dp
      use control_flags, only: iprint
      use gvecs
      use gvecw, only: ngw
      use cvan, only: ish
      use uspp, only: nhsa=>nkb, dvan, deeq
      use uspp_param, only: nhm, nh
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: n => nbsp, ispin, f, nspin
      use constants, only: pi, fpi
      use ions_base, only: nsp, na, nat
      use gvecw, only: ggp
      use cell_base, only: tpiba2
      use ensemble_dft, only: tens
      use funct, only: dft_is_meta
      USE task_groups
      use fft_base,  only : dffts
      use mp_global, only : nogrp, me_image, me_ogrp
      USE cp_interfaces, ONLY: fwfft, invfft
      use parallel_include
!
      implicit none
!
      !--------
      !C. Bekas
      !   c and ca hold the coefficients for the input eigenvalues
      !   originaly they are vectors of length ngw
      !   In the task-groups version they are matrices with
      !   ngw rows and NOGRP columns
      !-----------------------------------------------------------

      !--------
      !C. Bekas
      !--------
      !Observe the increased sizes for Task Groups
      !C. Bekas: Increased size for matrix v
      complex(DP) :: betae(ngw,nhsa), c(ngw,2*NOGRP), ca(ngw,2*NOGRP), df(ngw*(NOGRP+1)), da(ngw*(NOGRP+1)) 
      real(DP) :: bec(nhsa,n), v( ( NOGRP + 1 ) * nr1sx * nr2sx * nr3sx, nspin ) 
      integer i
! local variables
      integer iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      real(DP) fi, fip, dd
      complex(DP) fp,fm,ci
      real(DP),    ALLOCATABLE :: af(:,:), aa(:,:) 
                               ! C. Bekas: increased size for automatic arrays
      complex(DP), ALLOCATABLE :: dtemp(:,:)    
      complex(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: temp_psi( : )

      !--------
      !C. Bekas
      !--------
      INTEGER  ::  eig_index, index, index_df_da, ierr
      INTEGER, DIMENSION(NOGRP) :: local_send_cnt, local_send_displ, local_recv_cnt, local_recv_displ
!
      call start_clock( 'dforce' ) 
      !
#ifdef __BGL

      ALLOCATE( psi( strd * ( NOGRP+1 ) ))
      ALLOCATE( temp_psi( 2 * (NOGRP+1) * dffts%nsw(1) * nr3sx ) )
      ALLOCATE( af( nhsa, NOGRP ), aa( nhsa, NOGRP ), dtemp( ngw, NOGRP ) )
!
!     important: if n is odd => c(*,n+1)=0.
! 
      if ( MOD(n,2) .ne. 0 .and. i .eq. n ) then
         do ig = 1, ngw
            ca(ig,:) = 0.0d0
         end do
      end if
!
      ci = ( 0.0d0, 1.0d0 )
!

         psi(:) = (0D0,0D0)
         index = 1
         eig_offset = 0
         do eig_index = 1, 2*NOGRP, 2! Outer loop for eigenvalues
            !The  eig_index loop is executed only ONCE when NOGRP=1.
            !Equivalent to the case with no task-groups
            !dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
            !We can either send these in the group with an mpi_allgather...or put the
            !in the PSIS vector (in special positions) and send them with them.
            !Otherwise we can do this once at the beginning, before the loop.
            !we choose to do the latter one.

            !---------------------------------------------
            !strd is defined earlier in the rhoofr routine
            !---------------------------------------------

            do ig=1,ngw
               psi(nms(ig)+eig_offset*strd)=conjg( c(ig,index) - ci*ca(ig,index) )
               psi(nps(ig)+eig_offset*strd)=c(ig,index)+ci*ca(ig,index)
            end do
            eig_offset = eig_offset + 1
            index = index + 2
         end do

         CALL invfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)  

         ! task group is managed inside the fft driver

! 
      !==================================================================
      !C. Bekas
      !This logic is altered in the TG case, see below
      !------------------------------------------------------------------
      iss1=ispin(i)

!
      if (i.ne.n) then
         iss2=ispin(i+1)
      else
         iss2=iss1
      end if
      !==================================================================


      !------------------------------------------------------------------
      !Each wave function is multiplied term - to - term by the local
      !potential, which is always the same for all eigenvalues
      !The length of psi is so that it holds all parts of it in the
      !plane-wave group
      !------------------------------------------------------------------
      do ir=1, nr1sx*nr2sx*tmp_npp(me_image+1)
         psi(ir)=cmplx(v(ir,iss1)* DBLE(psi(ir)),                       &
     &                 v(ir,iss2)*AIMAG(psi(ir)) )
      end do
!

      !-----------------------------------------------
      !CALL TASK GROUP PARALLEL FORWARD FFT
      !Note that the wavefunctions are already
      !distributed according to the TASK-GROUPS
      !scheme
      !-----------------------------------------------

      CALL fwfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)


      !-------------------------------------------------
      !Bring pencils back to their original distribution
      !-------------------------------------------------
      local_send_cnt(1) = nr3sx*dffts%nsw(NOLIST(1)+1)
      local_send_displ(1) = 0
      local_recv_cnt(1) = nr3sx*dffts%nsw(me_image+1)
      local_recv_displ(1) = 0
      DO index=2, NOGRP
         local_send_cnt(index) = nr3sx*dffts%nsw(NOLIST(index)+1)
         local_send_displ(index) = local_send_displ(index-1) + local_send_cnt(index-1)

         local_recv_cnt(index) = nr3sx*dffts%nsw(me_image+1)
         local_recv_displ(index)  = local_recv_displ(index-1) + local_recv_cnt(index-1)
      ENDDO

      CALL start_clock('DFORCE_ALL')
#if defined __MPI
      CALL MPI_Alltoallv(psi, &
           local_send_cnt, local_send_displ, MPI_DOUBLE_COMPLEX, temp_psi, &
           local_recv_cnt, local_recv_displ, MPI_DOUBLE_COMPLEX, ME_OGRP, IERR)
#endif
      CALL stop_clock('DFORCE_ALL')


!
!     note : the factor 0.5 appears 
!       in the kinetic energy because it is defined as 0.5*g**2
!       in the potential part because of the logics
!
   

      !--------------------------------------------------------------
      !Each processor will treat its own part of the eigenstate
      !assigned to its ORBITAL group
      !--------------------------------------------------------------
      eig_offset = 0
      index_df_da = 1
      DO index = 1, 2*NOGRP, 2
         do ig=1,ngw
            if (tens) then
               fi = -0.5
               fip = -0.5
            else
               fi = -0.5*f(i+index-1)
               fip = -0.5*f(i+index)
            endif
            fp= temp_psi(nps(ig)+eig_offset) +  temp_psi(nms(ig)+eig_offset)
            fm= temp_psi(nps(ig)+eig_offset) -  temp_psi(nms(ig)+eig_offset)
            df(index_df_da)= fi*(tpiba2 * ggp(ig) * c(ig,index)+cmplx(real(fp), aimag(fm)))
            da(index_df_da)= fip*(tpiba2 * ggp(ig) * ca(ig,index)+cmplx(aimag(fp),-real(fm)))
            index_df_da = index_df_da + 1
         enddo
         eig_offset = eig_offset + nr3sx * dffts%nsw(me_image+1)
         !We take into account the number of elements received from other members of the orbital group
      ENDDO

      !--------------------------------------------------------------------------------------------
      !C. Bekas: I am not sure whether this is implemented correctly...need to check this carefully 
      if(dft_is_meta()) call dforce_meta(c,ca,df,da,psi,iss1,iss2,fi,fip) !METAGGA
      !--------------------------------------------------------------------------------------------
!
!     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
! 
      IF( nhsa > 0 ) THEN

         do inl=1,nhsa
            af(inl,:)=0.0d0
            aa(inl,:)=0.0d0
         end do
!
         do is=1,nsp
            do iv=1,nh(is)
               do jv=1,nh(is)
                  isa=0
                  do ism=1,is-1
                     isa=isa+na(ism)
                  end do
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     isa=isa+1
                     dd = deeq(iv,jv,isa,iss1)+dvan(iv,jv,is)
                    
                     !-------------------------------------------------
                     !C. Bekas
                     !Work on all currently treated (NOGRP) eigenvalues
                     !-------------------------------------------------
                     ig = 1
                     DO index = 1, 2*NOGRP, 2
                        if (tens) then 
                           af(inl,ig) = af(inl,ig) -  dd*bec(jnl,  i+index-1 )
                        else
                           af(inl,ig) = af(inl,ig) -  f(i+index-1)*dd*bec(jnl,  i+index-1 )
                        endif

                        dd = deeq(iv,jv,isa,iss2)+dvan(iv,jv,is)
                      
                        if (tens) then
                           if ((i+index-1).ne.n) aa(inl,ig) = aa(inl,ig) - dd*bec(jnl,i+index)
                        else
                           if ((i+index-1).ne.n) aa(inl,ig) = aa(inl,ig) - f(i+index)*dd*bec(jnl,i+index)
                        endif    
                        ig = ig + 1
                     ENDDO
                  end do
               end do
            end do
         end do
!
         dtemp(:,:) = 0.0d0

         call MXMA (betae, 1, 2*ngw, af, 1, nhsa, dtemp, 1, 2*ngw, 2*ngw, nhsa, NOGRP)

!         DO index = 1, NOGRP
!            DO ig = 1+(index-1)*ngw, index*ngw
!               df(ig) = df(ig) + dtemp(ig,index)
!            END DO
!         ENDDO

         DO ig = 1, NOGRP
            df(1+(ig-1)*ngw:ig*ngw) = df(1+(ig-1)*ngw:ig*ngw) + dtemp(:,ig)
         ENDDO


         dtemp(:,:) = 0.0d0

         call MXMA (betae, 1, 2*ngw, aa, 1, nhsa, dtemp, 1, 2*ngw, 2*ngw, nhsa, NOGRP)

!         DO index = 1, NOGRP
!            DO ig = 1+(index-1)*ngw, index*ngw
!               da(ig) = da(ig) + dtemp(ig,index)
!            ENDDO
!         ENDDO

         DO ig = 1, NOGRP
            da(1+(ig-1)*ngw:ig*ngw) = da(1+(ig-1)*ngw:ig*ngw) + dtemp(:,ig)
         ENDDO

      END IF


      DEALLOCATE( psi )
      DEALLOCATE( temp_psi ) 
      DEALLOCATE( af, aa, dtemp )

#endif
!
      call stop_clock( 'dforce' ) 
!
      return
      end subroutine dforce_bgl
!
