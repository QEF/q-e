!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written and revised by Carlo Cavazzoni
! Task Groups parallelization by C. Bekas (IBM Research Zurich).
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

      fioby2   = fio * 0.5d0
      fieby2   = fie * 0.5d0

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

     

    SUBROUTINE dforce_fpmd_x( ib, c, f, df, da, v, vkb, bec, n, noffset )
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
    END SUBROUTINE dforce_fpmd_x


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
             CALL dforce_fpmd( ib, c, f, cgrad(:,in), cgrad(:,in+1), vpot, vkb, bec, n, noffset )
             !
          END DO
          !
          !   and now process the last state in case that n is odd
          !
          IF( MOD( n, 2 ) /= 0 ) THEN
             !
             in = n + noffset - 1
             !
             CALL dforce_fpmd( n, c, f, cgrad(:,in), cgrad(:,in), vpot, vkb, bec, n, noffset )
             !
          END IF
          !
       END IF
       !
       RETURN
    END SUBROUTINE dforce_all



!
!-------------------------------------------------------------------------
      SUBROUTINE dforce_x ( i, bec, vkb, c, df, da, v, ldv, ispin, f, n, nspin, v1 )
!-----------------------------------------------------------------------
!computes: the generalized force df=CMPLX(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=CMPLX(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
!
      USE parallel_include
      USE kinds,                  ONLY: dp
      USE control_flags,          ONLY: iprint, use_task_groups
      USE gvecs,                  ONLY: nms, nps
      USE cvan,                   ONLY: ish
      USE uspp,                   ONLY: nhsa=>nkb, dvan, deeq
      USE uspp_param,             ONLY: nhm, nh
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
                                        nr1sx, nr2sx, nr3sx, nnrsx
      USE constants,              ONLY: pi, fpi
      USE ions_base,              ONLY: nsp, na, nat
      USE gvecw,                  ONLY: ngw, ggp
      USE cell_base,              ONLY: tpiba2
      USE ensemble_dft,           ONLY: tens
      USE fft_base,               ONLY: dffts
      USE funct,                  ONLY: dft_is_meta
      USE cp_interfaces,          ONLY: fwfft, invfft
      USE mp_global,              ONLY: nogrp, me_image, ogrp_comm
      USE task_groups,            ONLY: tmp_npp, nolist, strd
!
      IMPLICIT NONE
!
      INTEGER,     INTENT(IN)    :: i
      REAL(DP)                   :: bec(:,:)
      COMPLEX(DP)                :: vkb(:,:)
      COMPLEX(DP)                :: c(:,:)
      COMPLEX(DP)                :: df(:), da(:)
      INTEGER,     INTENT(IN)    :: ldv
      REAL(DP)                   :: v( ldv, * )
      INTEGER                    :: ispin( : )
      REAL(DP)                   :: f( : )
      INTEGER,     INTENT(IN)    :: n, nspin
      REAL(DP),    OPTIONAL      :: v1( ldv, * )
      !
      !
      ! local variables
      !
      INTEGER     :: iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      INTEGER     :: ivoff, jvoff, igoff, igno, igrp, ierr
      INTEGER     :: idx, eig_offset, eig_index, nogrp_
      REAL(DP)    :: fi, fip, dd
      COMPLEX(DP) :: fp, fm, ci
      REAL(DP),    ALLOCATABLE :: af( :, : ), aa( :, : ) ! automatic arrays
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: temp_psi(:)
      INTEGER,     ALLOCATABLE :: send_cnt(:), recv_cnt(:), send_displ(:), recv_displ(:)
!
!
      CALL start_clock( 'dforce' ) 
      !
      IF( use_task_groups ) THEN
         nogrp_ = nogrp
         ALLOCATE( psi( strd * ( nogrp + 1 ) ) )
         ALLOCATE( temp_psi( 2 * ( nogrp + 1 ) * dffts%nsw(1) * nr3sx ) )
      ELSE
         nogrp_ = 1
         ALLOCATE( psi( nnrsx ) )
         ALLOCATE( temp_psi( 1 ) )
      END IF
      !
      ci=(0.0d0,1.0d0)
      !
      psi (:) = (0.d0, 0.d0)

      igoff = 0

      DO idx = 1, 2*nogrp_ , 2

         !This loop is executed only ONCE when NOGRP=1.
         !Equivalent to the case with no task-groups
         !dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
         !We can either send these in the group with an mpi_allgather...or put the
         !in the PSIS vector (in special positions) and send them with them.
         !Otherwise we can do this once at the beginning, before the loop.
         !we choose to do the latter one.

         !     important: if n is odd => c(*,n+1)=0.
         ! 
         IF ( ( idx + i - 1 ) == n ) c( : , idx + i ) = 0.0d0

         IF( idx + i - 1 <= n ) THEN
            DO ig=1,ngw
               psi(nms(ig)+igoff) = conjg( c(ig,idx+i-1) - ci * c(ig,idx+i) )
               psi(nps(ig)+igoff) =        c(ig,idx+i-1) + ci * c(ig,idx+i)
            END DO
         END IF

         igoff = igoff + strd

      END DO

      CALL invfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!     
      !
      ! the following avoids a potential out-of-bounds error
      !
      IF ( i < n ) THEN
         iss1 = ispin(i)
         iss2 = ispin(i+1)
      ELSE
         iss1 = ispin(i)
         iss2 = iss1
      END IF
!
      IF( use_task_groups ) THEN
         !
         DO ir = 1, nr1sx * nr2sx * tmp_npp( me_image + 1 )
            psi(ir) = CMPLX( v(ir,iss1) * DBLE( psi(ir) ), v(ir,iss2) * AIMAG( psi(ir) ) )
         END DO
         !
      ELSE
         !
         IF( PRESENT( v1 ) ) THEN
            DO ir=1,nnrsx
               psi(ir)=CMPLX(v(ir,iss1)* DBLE(psi(ir)), v1(ir,iss2)*AIMAG(psi(ir)) )
            END DO
         ELSE
            DO ir=1,nnrsx
               psi(ir)=CMPLX(v(ir,iss1)* DBLE(psi(ir)), v(ir,iss2)*AIMAG(psi(ir)) )
            END DO
         END IF
         !
      END IF
      !
      CALL fwfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
      !
      IF( use_task_groups ) THEN
         !
         !  Bring pencils back to their original distribution
         !
         ALLOCATE( send_cnt( nogrp ), recv_cnt( nogrp ) )
         ALLOCATE( send_displ( nogrp ), recv_displ( nogrp ) )
         send_cnt  (1) = nr3sx*dffts%nsw(nolist(1)+1)
         send_displ(1) = 0
         recv_cnt  (1) = nr3sx*dffts%nsw(me_image+1)
         recv_displ(1) = 0
         DO idx = 2, NOGRP
            send_cnt  (idx) = nr3sx*dffts%nsw(nolist(idx)+1)
            send_displ(idx) = send_displ(idx-1) + send_cnt(idx-1)

            recv_cnt  (idx) = nr3sx*dffts%nsw(me_image+1)
            recv_displ(idx) = recv_displ(idx-1) + recv_cnt(idx-1)
         ENDDO

         CALL start_clock('DFORCE_ALL')
#if defined __MPI
         CALL MPI_Alltoallv( psi(1), &
              send_cnt, send_displ, MPI_DOUBLE_COMPLEX, temp_psi(1), &
              recv_cnt, recv_displ, MPI_DOUBLE_COMPLEX, ogrp_comm, IERR)
#endif
         CALL stop_clock('DFORCE_ALL')
         !
         DEALLOCATE( send_cnt, recv_cnt )
         DEALLOCATE( send_displ, recv_displ )
         !
      END IF
!
      !     note : the factor 0.5 appears 
      !       in the kinetic energy because it is defined as 0.5*g**2
      !       in the potential part because of the logics
!
   
      !
      !--------------------------------------------------------------
      !Each processor will treat its own part of the eigenstate
      !assigned to its ORBITAL group
      !--------------------------------------------------------------
      !
      eig_offset = 0
      igno = 1

      DO idx = 1, 2*nogrp_ , 2

         IF( idx + i - 1 <= n ) THEN
            if (tens) then
               fi = -0.5d0
               fip = -0.5d0
            else
               fi = -0.5d0*f(i+idx-1)
               fip = -0.5d0*f(i+idx)
            endif
            IF( use_task_groups ) THEN
               DO ig=1,ngw
                  fp= temp_psi(nps(ig)+eig_offset) +  temp_psi(nms(ig)+eig_offset)
                  fm= temp_psi(nps(ig)+eig_offset) -  temp_psi(nms(ig)+eig_offset)
                  df(igno)= fi *(tpiba2 * ggp(ig) * c(ig,idx+i-1)+cmplx(real (fp), aimag(fm)))
                  da(igno)= fip*(tpiba2 * ggp(ig) * c(ig,idx+i  )+cmplx(aimag(fp),-real (fm)))
                  igno = igno + 1
               END DO
            ELSE
               DO ig=1,ngw
                  fp= psi(nps(ig)) + psi(nms(ig))
                  fm= psi(nps(ig)) - psi(nms(ig))
                  df(ig)= fi*(tpiba2*ggp(ig)* c(ig,idx+i-1)+CMPLX(DBLE(fp), AIMAG(fm)))
                  da(ig)=fip*(tpiba2*ggp(ig)* c(ig,idx+i  )+CMPLX(AIMAG(fp),-DBLE(fm)))
               END DO
            END IF
         END IF

         eig_offset = eig_offset + nr3sx * dffts%nsw(me_image+1)
         !We take into account the number of elements received from other members of the orbital group

      ENDDO

      !

      IF(dft_is_meta()) THEN
         CALL dforce_meta(c(1,i),c(1,i+1),df,da,psi,iss1,iss2,fi,fip) !METAGGA
      END IF


      IF( nhsa > 0 ) THEN
         !
         !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
         ! 
         ALLOCATE( af( nhsa, nogrp_ ), aa( nhsa, nogrp_ ) )

         af = 0.0d0
         aa = 0.0d0
         !
         igrp = 1

         DO idx = 1, 2*nogrp_ , 2

            IF( idx + i - 1 <= n ) THEN

               IF (tens) THEN
                  fi = 1.0d0
                  fip= 1.0d0
               ELSE
                  fi = f(i+idx-1)
                  fip= f(i+idx)
               END IF
               !
               DO is = 1, nsp
                  DO iv = 1, nh(is)
                     DO jv = 1, nh(is)
                        isa = 0
                        DO ism = 1, is-1
                           isa = isa + na( ism )
                        END DO
                        ivoff = ish(is)+(iv-1)*na(is)
                        jvoff = ish(is)+(jv-1)*na(is)
                        IF( i + idx - 1 /= n ) THEN
                           DO ia=1,na(is)
                              inl = ivoff + ia
                              jnl = jvoff + ia
                              isa = isa + 1
                              dd = deeq(iv,jv,isa,iss1) + dvan(iv,jv,is)
                              af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                              dd = deeq(iv,jv,isa,iss2) + dvan(iv,jv,is)
                              aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                           END DO
                        ELSE
                           DO ia=1,na(is)
                              inl = ivoff + ia
                              jnl = jvoff + ia
                              isa = isa + 1
                              dd = deeq(iv,jv,isa,iss1) + dvan(iv,jv,is)
                              af(inl,igrp) = af(inl,igrp) - fi * dd * bec(jnl,i+idx-1)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
      
            END IF

            igrp = igrp + 1

         END DO
!
         CALL DGEMM ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, af, nhsa, 1.0d0, df, 2*ngw)

         CALL DGEMM ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, aa, nhsa, 1.0d0, da, 2*ngw)
         !
         DEALLOCATE( aa, af )
         !
      ENDIF

      DEALLOCATE( psi )
      DEALLOCATE( temp_psi )
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
   END SUBROUTINE dforce_x

