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
      USE control_flags,          ONLY: iprint, use_task_groups, program_name
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
      USE task_groups,            ONLY: tmp_npp, nolist, strd, nswx
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
      REAL(DP)    :: fi, fip, dd, dv
      COMPLEX(DP) :: fp, fm, ci
      REAL(DP),    ALLOCATABLE :: af( :, : ), aa( :, : )
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      !
      CALL start_clock( 'dforce' ) 
      !
      IF( use_task_groups ) THEN
         nogrp_ = nogrp
         ALLOCATE( psi( strd * ( nogrp + 1 ) ) )
      ELSE
         nogrp_ = 1
         ALLOCATE( psi( nnrsx ) )
      END IF
      !
      ci = ( 0.0d0, 1.0d0 )
      !
      psi( : ) = (0.d0, 0.d0)

      igoff = 0

      DO idx = 1, 2*nogrp_ , 2
         !
         !  This loop is executed only ONCE when NOGRP=1.
         !  Equivalent to the case with no task-groups
         !  dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
         !  We can either send these in the group with an mpi_allgather...or put the
         !  in the PSIS vector (in special positions) and send them with them.
         !  Otherwise we can do this once at the beginning, before the loop.
         !  we choose to do the latter one.
         !
         !  important: if n is odd => c(*,n+1)=0.
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

      CALL invfft( 'Wave', psi, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx )
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
      CALL fwfft( 'Wave', psi, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx )
      !
      !   note : the factor 0.5 appears 
      !       in the kinetic energy because it is defined as 0.5*g**2
      !       in the potential part because of the logics
      !
      !   Each processor will treat its own part of the eigenstate
      !   assigned to its ORBITAL group
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
                  fp= psi(nps(ig)+eig_offset) +  psi(nms(ig)+eig_offset)
                  fm= psi(nps(ig)+eig_offset) -  psi(nms(ig)+eig_offset)
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

         ! We take into account the number of elements received from other members of the orbital group

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
                     IF( program_name == 'FPMD' ) THEN
                        ivoff = ish(is) + (iv-1) * na(is)
                        dd = dvan( iv, iv, is )
                        DO inl = ivoff + 1, ivoff + na(is)
                           af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(inl,i+idx-1)
                        END DO
                        IF( i + idx - 1 /= n ) THEN
                           DO inl = ivoff + 1, ivoff + na(is)
                              aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(inl,i+idx)
                           END DO
                        END IF
                     ELSE
                     DO jv = 1, nh(is)
                        isa = 0
                        DO ism = 1, is-1
                           isa = isa + na( ism )
                        END DO
                        dv = dvan(iv,jv,is)
                        ivoff = ish(is)+(iv-1)*na(is)
                        jvoff = ish(is)+(jv-1)*na(is)
                        IF( i + idx - 1 /= n ) THEN
                           DO ia=1,na(is)
                              inl = ivoff + ia
                              jnl = jvoff + ia
                              isa = isa + 1
                              dd = deeq(iv,jv,isa,iss1) + dv
                              af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                              dd = deeq(iv,jv,isa,iss2) + dv
                              aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                           END DO
                        ELSE
                           DO ia=1,na(is)
                              inl = ivoff + ia
                              jnl = jvoff + ia
                              isa = isa + 1
                              dd = deeq(iv,jv,isa,iss1) + dv
                              af(inl,igrp) = af(inl,igrp) - fi * dd * bec(jnl,i+idx-1)
                           END DO
                        END IF
                     END DO
                     END IF
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
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
   END SUBROUTINE dforce_x

