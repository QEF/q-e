!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!
!-----------------------------------------------------------------------
      SUBROUTINE atomic_wfc(eigr,n_atomic_wfc,wfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
!
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart, g, gx
      USE ions_base, ONLY: nsp, na, nat
      USE cell_base, ONLY: tpiba
      USE atom, ONLY: nchi, lchi, mesh, r, chi, rab
!
      IMPLICIT NONE
      INTEGER, INTENT(in) :: n_atomic_wfc
      COMPLEX(8), INTENT(in) ::  eigr(ngw,nat)
      COMPLEX(8), INTENT(out):: wfc(ngw,n_atomic_wfc)
!
      INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
      REAL(8), ALLOCATABLE::  ylm(:,:), q(:), jl(:), vchi(:),      &
     &     chiq(:)
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         DO nb = 1, nchi(is)
            lmax_wfc = MAX (lmax_wfc, lchi (nb, is) )
         ENDDO
      ENDDO
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, gx, g, ylm)
      ndm = MAXVAL(mesh(1:nsp))
      ALLOCATE(jl(ndm), vchi(ndm))
      ALLOCATE(q(ngw), chiq(ngw))
!
      DO i=1,ngw
         q(i) = SQRT(g(i))*tpiba
      END DO
!
      natwfc=0
      isa   = 0
      DO is=1,nsp
!
!   radial fourier transform of the chi functions
!   NOTA BENE: chi is r times the radial part of the atomic wavefunction
!
         DO nb = 1,nchi(is)
            l = lchi(nb,is)
            DO i=1,ngw
               CALL sph_bes (mesh(is), r(1,is), q(i), l, jl)
               DO ir=1,mesh(is)
                  vchi(ir) = chi(ir,nb,is)*r(ir,is)*jl(ir)
               ENDDO
               CALL simpson_cp90(mesh(is),vchi,rab(1,is),chiq(i))
            ENDDO
!
!   multiply by angular part and structure factor
!   NOTA BENE: the factor i^l MUST be present!!!
!
            DO m = 1,2*l+1
               lm = l**2 + m
               DO ia = 1 + isa, na(is) + isa
                  natwfc = natwfc + 1
                  wfc(:,natwfc) = (0.d0,1.d0)**l * eigr(:,ia)* ylm(:,lm)*chiq(:)
               ENDDO
            ENDDO
         ENDDO
         isa = isa + na(is)
      ENDDO
!
      IF (natwfc.NE.n_atomic_wfc)                                       &
     &     CALL errore('atomic_wfc','unexpected error',natwfc)
!
      DEALLOCATE(q, chiq, vchi, jl, ylm)
!
      RETURN
      END SUBROUTINE atomic_wfc
!
!

!-----------------------------------------------------------------------
      REAL(8) FUNCTION cscnorm( bec, nkbx, cp, ngwx, i, n )
!-----------------------------------------------------------------------
!     requires in input the updated bec(i)
!
      USE ions_base,  ONLY: na
      USE gvecw,      ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE cvan,       ONLY: ish, nvb
      USE uspp_param, ONLY: nh
      USE uspp,       ONLY: qq
      USE mp,         ONLY: mp_sum
      USE kinds,      ONLY: DP
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, n
      INTEGER, INTENT(IN) :: ngwx, nkbx
      REAL(DP)    :: bec( nkbx, n )
      COMPLEX(DP) :: cp( ngwx, n )
!
      INTEGER ig, is, iv, jv, ia, inl, jnl
      REAL(8) rsum
      REAL(8), ALLOCATABLE:: temp(:)
!
!
      ALLOCATE(temp(ngw))
      DO ig=1,ngw
         temp(ig)=DBLE(CONJG(cp(ig,i))*cp(ig,i))
      END DO
      rsum=2.*SUM(temp)
      IF (gstart == 2) rsum=rsum-temp(1)

      CALL mp_sum( rsum )

      DEALLOCATE(temp)
!
      DO is=1,nvb
         DO iv=1,nh(is)
            DO jv=1,nh(is)
               IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN 
                  DO ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     rsum = rsum +                                        &
     &                    qq(iv,jv,is)*bec(inl,i)*bec(jnl,i)
                  END DO
               ENDIF
            END DO
         END DO
      END DO
!
      cscnorm=SQRT(rsum)
!
      RETURN
      END FUNCTION cscnorm
!
!-----------------------------------------------------------------------
      SUBROUTINE denkin(c,dekin)
!-----------------------------------------------------------------------
!
      USE constants, ONLY: pi, fpi
      USE electrons_base, ONLY: n => nbsp, nx => nbspx, f
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart, g, gx
      USE cell_base, ONLY: ainv, tpiba2
      USE gvecw, ONLY: ggp, ecutz, ecsig, ecfix
      USE mp, ONLY: mp_sum
!
      IMPLICIT NONE
! input
      COMPLEX(8) c(ngw,nx)
! output
      REAL(8) dekin(3,3)
! local
      INTEGER j, k, ig, i
      REAL(8), ALLOCATABLE:: gtmp(:)
      REAL(8) sk(n)  ! automatic array
      REAL(8) :: ga, dggp, efac
!
      ALLOCATE (gtmp(ngw))
      dekin=0.d0
      DO j=1,3
         DO k=1,3
            DO ig=1,ngw
               efac     = 2.d0 * ecutz / ecsig / SQRT(pi)
               dggp     = 1.d0 + efac * EXP( - ( tpiba2 * g(ig) - ecfix ) * ( tpiba2 * g(ig) - ecfix ) / ecsig / ecsig )
               ga       = gx(1,ig) * ainv(k,1) + gx(2,ig) * ainv(k,2) + gx(3,ig) * ainv(k,3)
               gtmp(ig) = gx(j,ig) * ga * dggp
            END DO
            DO i=1,n
               sk(i)=0.d0
               DO ig=gstart,ngw
                  sk(i)=sk(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*gtmp(ig)
               END DO
            END DO
            DO i=1,n
               dekin(j,k)=dekin(j,k)-2.d0*tpiba2*(f(i)*sk(i))
            END DO
         END DO
      END DO
      DEALLOCATE (gtmp)

      CALL mp_sum( dekin( 1:3, 1:3 ) )
!
      RETURN
      END SUBROUTINE denkin
!
!-----------------------------------------------------------------------
      SUBROUTINE denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
!-----------------------------------------------------------------------
!
! derivative of hartree energy wrt cell parameters h
! Output in dh
!
! rhotmp input : total electronic + ionic broadened charge (G)
! drhotmp input and work space
! sfac   input : structure factors
! wtemp work space
! eh input: hartree energy
!
      USE constants, ONLY: pi, fpi
      USE ions_base, ONLY: nsp
      USE gvecs
      USE gvecp, ONLY: ng => ngm
      USE reciprocal_vectors, ONLY: gstart, gx, g
      USE cell_base, ONLY: omega
      USE cell_base, ONLY: ainv, tpiba2
      USE local_pseudo, ONLY: rhops, drhops
      USE mp, ONLY: mp_sum

      IMPLICIT NONE
! input
      COMPLEX(8) rhotmp(ng), drhotmp(ng,3,3), vtemp(ng), sfac(ngs,nsp)
      REAL(8) eh
! output
      REAL(8) dh(3,3)
! local
      INTEGER i, j, ig, is
      REAL(8) wz
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.d0
      DO j=1,3
         DO i=1,3
            DO is=1,nsp
               DO ig=1,ngs
                  drhotmp(ig,i,j) = drhotmp(ig,i,j) -                   &
     &                    sfac(ig,is)*drhops(ig,is)*                    &
     &                    2.d0*tpiba2*gx(i,ig)*(gx(1,ig)*ainv(j,1)+     &
     &                     gx(2,ig)*ainv(j,2)+gx(3,ig)*ainv(j,3))-      &
     &                    sfac(ig,is)*rhops(ig,is)*ainv(j,i)
               ENDDO
            ENDDO
            IF (gstart == 2) vtemp(1)=(0.d0,0.d0)
            DO ig=gstart,ng
               vtemp(ig)=CONJG(rhotmp(ig))*rhotmp(ig)/(tpiba2*g(ig))**2 &
     &                 * tpiba2*gx(i,ig)*(gx(1,ig)*ainv(j,1)+           &
     &                   gx(2,ig)*ainv(j,2)+gx(3,ig)*ainv(j,3)) +       &
     &                 CONJG(rhotmp(ig))/(tpiba2*g(ig))*drhotmp(ig,i,j)
            ENDDO
            dh(i,j)=fpi*omega*DBLE(SUM(vtemp))*wz
         ENDDO
      ENDDO

      CALL mp_sum( dh( 1:3, 1:3 ) )

      DO i=1,3
         DO j=1,3
            dh(i,j)=dh(i,j)+omega*eh*ainv(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE denh
!
!-----------------------------------------------------------------------
      SUBROUTINE denps(rhotmp,drhotmp,sfac,vtemp,dps)
!-----------------------------------------------------------------------
!
! derivative of local potential energy wrt cell parameters h
! Output in dps
!
! rhotmp input : rho(G) (up and down spin components summed)
! drhotmp input
! sfac   input : structure factors
! wtemp work space
!
      USE ions_base, ONLY: nsp
      USE gvecs, ONLY: ngs
      USE gvecp, ONLY: ng => ngm
      USE reciprocal_vectors, ONLY: gstart, gx
      USE cell_base, ONLY: omega
      USE cell_base, ONLY: ainv, tpiba2
      USE local_pseudo, ONLY: vps, dvps
      USE mp, ONLY: mp_sum

      IMPLICIT NONE
! input
      COMPLEX(8) rhotmp(ng), drhotmp(ng,3,3), vtemp(ng), sfac(ngs,nsp)
! output
      REAL(8) dps(3,3)
! local
      INTEGER i, j, ig, is
      REAL(8) wz
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.d0
      DO i=1,3
         DO j=1,3
            DO ig=1,ngs
               vtemp(ig)=(0.,0.)
            ENDDO
            DO is=1,nsp
               DO ig=1,ngs
                  vtemp(ig)=vtemp(ig)-CONJG(rhotmp(ig))*sfac(ig,is)*    &
     &                    dvps(ig,is)*2.d0*tpiba2*gx(i,ig)*             &
     &                    (gx(1,ig)*ainv(j,1) +                         &
     &                     gx(2,ig)*ainv(j,2) +                         &
     &                     gx(3,ig)*ainv(j,3) ) +                       &
     &                    CONJG(drhotmp(ig,i,j))*sfac(ig,is)*vps(ig,is)
               ENDDO
            ENDDO
            dps(i,j)=omega*DBLE(wz*SUM(vtemp))
            IF (gstart == 2) dps(i,j)=dps(i,j)-omega*DBLE(vtemp(1))
         ENDDO
      ENDDO

      CALL mp_sum( dps( 1:3, 1:3 ) )

      RETURN
      END SUBROUTINE denps


!-----------------------------------------------------------------------
      SUBROUTINE denlcc( nnr, nspin, vxcr, sfac, drhocg, dcc )
!-----------------------------------------------------------------------
!
! derivative of non linear core correction exchange energy wrt cell 
! parameters h 
! Output in dcc
!
      USE kinds, ONLY: DP
      USE ions_base, ONLY: nsp
      USE reciprocal_vectors, ONLY: gstart, gx, ngs, g, ngm
      USE recvecs_indexes, ONLY: np
      USE cell_base, ONLY: omega, ainv, tpiba2
      USE mp, ONLY: mp_sum
      USE atom, ONLY: nlcc
      USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x
      USE fft_module, ONLY: fwfft

      IMPLICIT NONE

      ! input

      INTEGER, INTENT(IN)   :: nnr, nspin
      REAL(DP)              :: vxcr( nnr, nspin )
      COMPLEX(DP)           :: sfac( ngs, nsp )
      REAL(DP)              :: drhocg( ngm, nsp )

      ! output

      REAL(DP), INTENT(OUT) ::  dcc(3,3)

      ! local

      INTEGER     :: i, j, ig, is
      COMPLEX(DP) :: srhoc
      REAL(DP)    :: vxcc
      !
      COMPLEX(DP), ALLOCATABLE :: vxc( : )
!
      dcc = 0.0d0
      !
      ALLOCATE( vxc( nnr ) )
      !
      vxc(:) = vxcr(:,1)
      !
      IF( nspin > 1 ) vxc(:) = vxc(:) + vxcr(:,2)
      !
      CALL fwfft( 'Dense', vxc, nr1, nr2, nr3, nr1x, nr2x, nr3x )
      !
      DO i=1,3
         DO j=1,3
            DO ig = gstart, ngs
               srhoc = 0.0d0
               DO is = 1, nsp
                 IF( nlcc( is ) ) srhoc = srhoc + sfac( ig, is ) * drhocg( ig, is )
               ENDDO
               vxcc = DBLE( CONJG( vxc( np( ig ) ) ) * srhoc ) / SQRT( g( ig ) * tpiba2 )
               dcc(i,j) = dcc(i,j) + vxcc * &
     &                      2.d0 * tpiba2 * gx(i,ig) *                  &
     &                    (gx(1,ig)*ainv(j,1) +                         &
     &                     gx(2,ig)*ainv(j,2) +                         &
     &                     gx(3,ig)*ainv(j,3) )
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( vxc )

      dcc = dcc * omega

      CALL mp_sum( dcc( 1:3, 1:3 ) )

      RETURN
      END SUBROUTINE denlcc



!
!-------------------------------------------------------------------------
      SUBROUTINE dforce (bec,betae,i,c,ca,df,da,v)
!-----------------------------------------------------------------------
!computes: the generalized force df=CMPLX(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=CMPLX(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
      USE kinds, ONLY: dp
      USE control_flags, ONLY: iprint, tbuff
      USE gvecs
      USE gvecw, ONLY: ngw
      USE cvan, ONLY: ish
      USE uspp, ONLY: nhsa=>nkb, dvan, deeq
      USE uspp_param, ONLY: nhm, nh
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      USE electrons_base, ONLY: n => nbsp, ispin, f, nspin
      USE constants, ONLY: pi, fpi
      USE ions_base, ONLY: nsp, na, nat
      USE gvecw, ONLY: ggp
      USE cell_base, ONLY: tpiba2
      USE ensemble_dft, ONLY: tens
      USE funct, ONLY: dft_is_meta
      USE fft_module, ONLY: fwfft, invfft
!
      IMPLICIT NONE
!
      COMPLEX(8) betae(ngw,nhsa), c(ngw), ca(ngw), df(ngw), da(ngw)
      REAL(8) bec(nhsa,n), v(nnrsx,nspin)
      INTEGER i
! local variables
      INTEGER iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      REAL(8) fi, fip, dd
      COMPLEX(8) fp,fm,ci
      REAL(8) af(nhsa), aa(nhsa) ! automatic arrays
      COMPLEX(8)  dtemp(ngw)    !
      COMPLEX(8), ALLOCATABLE :: psi(:)
!
!
      CALL start_clock( 'dforce' ) 
      !
      ALLOCATE( psi( nnrsx ) )
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
      IF (.NOT.tbuff) THEN
!
         psi (:) = (0.d0, 0.d0)
         DO ig=1,ngw
            psi(nms(ig))=CONJG(c(ig)-ci*ca(ig))
            psi(nps(ig))=c(ig)+ci*ca(ig)
         END DO
!
         CALL invfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!     
      ELSE
!
!     read psi from buffer 21
!
         READ(21,iostat=ios) psi
         IF(ios.NE.0) CALL errore                                        &
     &       (' dforce',' error in reading unit 21',ios)
!
      ENDIF
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
         DO ig=1,ngw
            dtemp(ig)=(0.,0.)
         END DO
         CALL MXMA                                                      &
     &        (betae,1,2*ngw,af,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         DO ig=1,ngw
            df(ig)=df(ig)+dtemp(ig)
         END DO
!
         DO ig=1,ngw
            dtemp(ig)=(0.,0.)
         END DO
         CALL MXMA                                                      &
     &        (betae,1,2*ngw,aa,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         DO ig=1,ngw
            da(ig)=da(ig)+dtemp(ig)
         END DO
      ENDIF

      DEALLOCATE( psi )
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
      END SUBROUTINE dforce
!
!-----------------------------------------------------------------------
      SUBROUTINE dotcsc(eigr,cp)
!-----------------------------------------------------------------------
!
      USE ions_base, ONLY: nas => nax, na, nsp, nat
      USE io_global, ONLY: stdout
      USE gvecw, ONLY: ngw
      USE electrons_base, ONLY: n => nbsp
      USE reciprocal_vectors, ONLY: gstart
      USE cvan, ONLY: ish, nvb
      USE uspp, ONLY: nhsa=>nkb, qq
      USE uspp_param, ONLY: nh
      USE mp, ONLY: mp_sum
!
      IMPLICIT NONE
!
      COMPLEX(8)  eigr(ngw,nat), cp(ngw,n)
! local variables
      REAL(8) rsum, csc(n) ! automatic array
      COMPLEX(8) temp(ngw) ! automatic array
 
      REAL(8), ALLOCATABLE::  becp(:,:)
      INTEGER i,kmax,nnn,k,ig,is,ia,iv,jv,inl,jnl
!
      ALLOCATE(becp(nhsa,n))
!
!     < beta | phi > is real. only the i lowest:
!
      nnn=MIN(12,n)
      DO i=nnn,1,-1
         kmax=i
         CALL nlsm1(i,1,nvb,eigr,cp,becp)
!
         DO k=1,kmax
            DO ig=1,ngw
               temp(ig)=CONJG(cp(ig,k))*cp(ig,i)
            END DO
            csc(k)=2.*DBLE(SUM(temp))
            IF (gstart == 2) csc(k)=csc(k)-DBLE(temp(1))
         END DO

         CALL mp_sum( csc( 1:kmax ) )

         DO k=1,kmax
            rsum=0.
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     DO ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        rsum = rsum +                                    &
     &                   qq(iv,jv,is)*becp(inl,i)*becp(jnl,k)
                     END DO
                  END DO
               END DO
            END DO
            csc(k)=csc(k)+rsum
         END DO
!
         WRITE( stdout,'(a,12f18.15)')' dotcsc = ',(csc(k),k=1,i)
!
      END DO
      WRITE( stdout,*)
!
      DEALLOCATE(becp)
!
      RETURN
      END SUBROUTINE dotcsc
!-----------------------------------------------------------------------
      SUBROUTINE drhov(irb,eigrb,rhovan,rhog,rhor)
!-----------------------------------------------------------------------
!     this routine calculates arrays drhog drhor, derivatives wrt h of:
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     Same logic as in routine rhov.
!     On input rhor and rhog must contain the smooth part only !!!
!     Output in module derho (drhor, drhog)
!
      USE kinds, ONLY: dp
      USE control_flags, ONLY: iprint
      USE ions_base, ONLY: na, nsp, nat, nas => nax
      USE cvan
      USE uspp_param, ONLY: nhm, nh
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      USE electrons_base, ONLY: nspin
      USE gvecb
      USE gvecp, ONLY: ng => ngm
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE cell_base, ONLY: ainv
      USE qgb_mod
      USE cdvan
      USE derho
      USE dqgb_mod
      USE recvecs_indexes, ONLY: nm, np
      USE fft_module, ONLY: fwfft, invfft
      USE fft_base, ONLY: dfftb

      IMPLICIT NONE
! input
      INTEGER, INTENT(in) ::  irb(3,nat)
      REAL(8), INTENT(in)::  rhor(nnr,nspin)
      REAL(8) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      COMPLEX(8), INTENT(in)::  eigrb(ngb,nat), rhog(ng,nspin)
! local
      INTEGER i, j, isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss,   &
     &     isa, ia, ir
      REAL(8) sum, dsum
      COMPLEX(8) fp, fm, ci
      COMPLEX(8), ALLOCATABLE :: v(:)
      COMPLEX(8), ALLOCATABLE:: dqgbt(:,:)
      COMPLEX(8), ALLOCATABLE :: qv(:)
!
!
      DO j=1,3
         DO i=1,3
            DO iss=1,nspin
               DO ir=1,nnr
                  drhor(ir,iss,i,j)=-rhor(ir,iss)*ainv(j,i)
               END DO
               DO ig=1,ng
                  drhog(ig,iss,i,j)=-rhog(ig,iss)*ainv(j,i)
               END DO
            END DO
         END DO
      END DO
!
      IF (nvb.EQ.0) RETURN
!
      ALLOCATE( v( nnr ) )
      ALLOCATE( qv( nnrb ) )
      ALLOCATE( dqgbt( ngb, 2 ) )

      ci=(0.,1.)
!
      IF(nspin.EQ.1) THEN
!     ------------------------------------------------------------------
!     nspin=1 : two fft at a time, one per atom, if possible
!     ------------------------------------------------------------------
         DO i=1,3
            DO j=1,3
!
               v(:) = (0.d0, 0.d0)
!
               iss=1
               isa=1
               DO is=1,nvb
#ifdef __PARA
                  DO ia=1,na(is)
                     nfft=1
                     IF ( dfftb%np3( isa ) <= 0 ) go to 15
#else
                  DO ia=1,na(is),2
                     nfft=2
#endif
                     dqgbt(:,:) = (0.d0, 0.d0) 
                     IF (ia.EQ.na(is)) nfft=1
!
!  nfft=2 if two ffts at the same time are performed
!
                     DO ifft=1,nfft
                        DO iv=1,nh(is)
                           DO jv=iv,nh(is)
                              ijv = (jv-1)*jv/2 + iv
                              sum = rhovan(ijv,isa+ifft-1,iss)
                              dsum=drhovan(ijv,isa+ifft-1,iss,i,j)
                              IF(iv.NE.jv) THEN
                                 sum =2.*sum
                                 dsum=2.*dsum
                              ENDIF
                              DO ig=1,ngb
                                 dqgbt(ig,ifft)=dqgbt(ig,ifft) +        &
     &                                (sum*dqgb(ig,ijv,is,i,j) +        &
     &                                dsum*qgb(ig,ijv,is) )
                              END DO
                           END DO
                        END DO
                     END DO
!     
! add structure factor
!
                     qv(:) = (0.d0, 0.d0)
                     IF(nfft.EQ.2) THEN
                        DO ig=1,ngb
                           qv(npb(ig)) = eigrb(ig,isa   )*dqgbt(ig,1)  &
     &                        + ci*      eigrb(ig,isa+1 )*dqgbt(ig,2)
                           qv(nmb(ig))=                                 &
     &                             CONJG(eigrb(ig,isa  )*dqgbt(ig,1)) &
     &                        + ci*CONJG(eigrb(ig,isa+1)*dqgbt(ig,2))
                        END DO
                     ELSE
                        DO ig=1,ngb
                           qv(npb(ig)) = eigrb(ig,isa)*dqgbt(ig,1)
                           qv(nmb(ig)) =                                &
     &                             CONJG(eigrb(ig,isa)*dqgbt(ig,1))
                        END DO
                     ENDIF
!
                     CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
!  qv = US contribution in real space on box grid
!       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
!
!  add qv(r) to v(r), in real space on the dense grid
!
                     CALL box2grid(irb(1,isa),1,qv,v)
                     IF (nfft.EQ.2) CALL box2grid(irb(1,isa+1),2,qv,v)
  15                 isa=isa+nfft
!
                  END DO
               END DO
!
               DO ir=1,nnr
                  drhor(ir,iss,i,j)=drhor(ir,iss,i,j)+DBLE(v(ir))
               END DO
!
               CALL fwfft('Dense', v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
               DO ig=1,ng
                  drhog(ig,iss,i,j)=drhog(ig,iss,i,j)+v(np(ig))
               END DO
!
            ENDDO
         ENDDO
!
      ELSE
!     ------------------------------------------------------------------
!     nspin=2: two fft at a time, one for spin up and one for spin down
!     ------------------------------------------------------------------
         isup=1
         isdw=2
         DO i=1,3
            DO j=1,3
               v(:) = (0.d0, 0.d0)
               isa=1
               DO is=1,nvb
                  DO ia=1,na(is)
#ifdef __PARA
                     IF ( dfftb%np3( isa ) <= 0 ) go to 25
#endif
                     DO iss=1,2
                        dqgbt(:,iss) = (0.d0, 0.d0)
                        DO iv= 1,nh(is)
                           DO jv=iv,nh(is)
                              ijv = (jv-1)*jv/2 + iv
                              sum=rhovan(ijv,isa,iss)
                              dsum =drhovan(ijv,isa,iss,i,j)
                              IF(iv.NE.jv) THEN
                                 sum =2.*sum
                                 dsum=2.*dsum
                              ENDIF
                              DO ig=1,ngb
                                 dqgbt(ig,iss)=dqgbt(ig,iss)  +         &
     &                               (sum*dqgb(ig,ijv,is,i,j) +         &
     &                               dsum*qgb(ig,ijv,is))
                              END DO
                           END DO
                        END DO
                     END DO
!     
! add structure factor
!
                     qv(:) = (0.d0, 0.d0)
                     DO ig=1,ngb
                        qv(npb(ig))= eigrb(ig,isa)*dqgbt(ig,1)        &
     &                    + ci*      eigrb(ig,isa)*dqgbt(ig,2)
                        qv(nmb(ig))= CONJG(eigrb(ig,isa)*dqgbt(ig,1)) &
     &                    +       ci*CONJG(eigrb(ig,isa)*dqgbt(ig,2))
                     END DO
!
                     CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
!  qv is the now the US augmentation charge for atomic species is
!  and atom ia: real(qv)=spin up, imag(qv)=spin down
!
!  add qv(r) to v(r), in real space on the dense grid
!
                     CALL box2grid2(irb(1,isa),qv,v)
  25                 isa=isa+1
                  END DO
               END DO
!
               DO ir=1,nnr
                  drhor(ir,isup,i,j) = drhor(ir,isup,i,j) + DBLE(v(ir))
                  drhor(ir,isdw,i,j) = drhor(ir,isdw,i,j) +AIMAG(v(ir))
               ENDDO
!
               CALL fwfft('Dense', v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
               DO ig=1,ng
                  fp=v(np(ig))+v(nm(ig))
                  fm=v(np(ig))-v(nm(ig))
                  drhog(ig,isup,i,j) = drhog(ig,isup,i,j) +             &
     &                 0.5*CMPLX( DBLE(fp),AIMAG(fm))
                  drhog(ig,isdw,i,j) = drhog(ig,isdw,i,j) +             &
     &                 0.5*CMPLX(AIMAG(fp),-DBLE(fm))
               END DO
!
            END DO
         END DO
      ENDIF
      DEALLOCATE(dqgbt)
      DEALLOCATE( v )
      DEALLOCATE( qv )
!
      RETURN
      END SUBROUTINE drhov

!
!-----------------------------------------------------------------------
   FUNCTION enkin( c, ngwx, f, n )
!-----------------------------------------------------------------------
      !
      ! calculation of kinetic energy term
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: pi, fpi
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE gvecw,              ONLY: ggp
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: group
      USE cell_base,          ONLY: tpiba2

      IMPLICIT NONE

      REAL(DP)                :: enkin

      ! input

      INTEGER,     INTENT(IN) :: ngwx, n
      COMPLEX(DP), INTENT(IN) :: c( ngwx, n )
      REAL(DP),    INTENT(IN) :: f( n )
      !
      ! local

      INTEGER  :: ig, i
      REAL(DP) :: sk(n)  ! automatic array
      !
      DO i=1,n
         sk(i)=0.0
         DO ig=gstart,ngw
            sk(i)=sk(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*ggp(ig)
         END DO
      END DO

      CALL mp_sum( sk(1:n), group )

      enkin=0.0
      DO i=1,n
         enkin=enkin+f(i)*sk(i)
      END DO

      ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
      ! ... multiplicative factor (2 pi/alat)**2 is required

      enkin = enkin * tpiba2
!
      RETURN
   END FUNCTION enkin
!
!
!-----------------------------------------------------------------------
      SUBROUTINE force_ion(tau0,esr,fion,dsr)
!-----------------------------------------------------------------------
!
!     forces on ions, ionic term in real space (also stress if requested)
!
      USE parameters, ONLY: natx
      USE control_flags, ONLY: iprint, tpre
      USE constants, ONLY: pi, fpi
      USE cell_base, ONLY: ainv, a1, a2, a3
      USE ions_base, ONLY: nsp, na, rcmax, zv
      IMPLICIT NONE
! input
      REAL(8) tau0(3,natx)
! output
      REAL(8) fion(3,natx), dsr(3,3), esr
! local variables
      INTEGER i,j,k,l,m, ii, lax, inf, isak, isaj
      REAL(8) rlm(3), rckj, rlmn, arg, addesr, addpre, repand, fxx
      REAL(8), EXTERNAL :: erfc
!
!
      esr=0.d0
      IF(tpre) dsr=0.d0
!
      isak = 0
      DO k=1,nsp
         isaj = 0
         DO j = 1, k-1
           isaj = isaj + na(j)
         END DO
         DO j=k,nsp
            rckj=SQRT(rcmax(k)**2+rcmax(j)**2)
            lax=na(k)
            IF(k.EQ.j) lax=lax-1
!
            DO l=1,lax
               inf=1
               IF(k.EQ.j) inf=l+1
!
               DO m=inf,na(j)
                  rlm(1) = tau0(1,l + isak) - tau0(1,m + isaj)
                  rlm(2) = tau0(2,l + isak) - tau0(2,m + isaj)
                  rlm(3) = tau0(3,l + isak) - tau0(3,m + isaj)
                  CALL pbc(rlm,a1,a2,a3,ainv,rlm)
!
                  rlmn=SQRT(rlm(1)**2+rlm(2)**2+rlm(3)**2)
!
                  arg=rlmn/rckj
                  addesr=zv(k)*zv(j)*erfc(arg)/rlmn
                  esr=esr+addesr
                  addpre=2.d0*zv(k)*zv(j)*EXP(-arg*arg)/rckj/SQRT(pi)
                  repand=(addesr+addpre)/rlmn/rlmn
!
                  DO i=1,3
                     fxx=repand*rlm(i)
                     fion(i,l+isak)=fion(i,l+isak)+fxx
                     fion(i,m+isaj)=fion(i,m+isaj)-fxx
                     IF(tpre)THEN
                        DO ii=1,3
                           dsr(i,ii)=dsr(i,ii)-                         &
     &                             repand*rlm(i)*rlm(1)*ainv(ii,1)-     &
     &                             repand*rlm(i)*rlm(2)*ainv(ii,2)-     &
     &                             repand*rlm(i)*rlm(3)*ainv(ii,3)
                        END DO
                     ENDIF
                  END DO
               END DO
            END DO
            isaj = isaj + na(j)
         END DO
         isak = isak + na(k)
      END DO

      RETURN
      END SUBROUTINE force_ion
!
!-----------------------------------------------------------------------
      SUBROUTINE force_ps(rhotemp,rhog,vtemp,ei1,ei2,ei3,fion1)
!-----------------------------------------------------------------------
!
! Contribution to ionic forces from local pseudopotential
!
      USE kinds, ONLY: dp
      USE constants, ONLY: pi, fpi
      USE electrons_base, ONLY: nspin
      USE gvecs
      USE gvecp, ONLY: ng => ngm
      USE reciprocal_vectors, ONLY: gstart, gx, mill_l, g
      USE cell_base, ONLY: omega, tpiba, tpiba2
      USE ions_base, ONLY: nsp, na, nas => nax, nat
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE parameters, ONLY: natx
      USE local_pseudo, ONLY: vps, rhops
!
      IMPLICIT NONE
! input
      COMPLEX(8) rhotemp(ng), rhog(ng,nspin), vtemp(ng),           &
     &           ei1(-nr1:nr1,nat),                                 &
     &           ei2(-nr2:nr2,nat),                                 &
     &           ei3(-nr3:nr3,nat)
! output
      REAL(8) fion1(3,natx)
! local
      INTEGER ig, is, isa, ism, ia, ix, iss, isup, isdw
      INTEGER i, j, k
      REAL(8)  wz
      COMPLEX(8) eigrx, vcgs, cnvg, cvn
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.0
      DO is=1,nsp
         isa=0
         DO ism=1,is-1
            isa=isa+na(ism)
         END DO
         DO ia=1,na(is)
            isa=isa+1
            DO ix=1,3
               IF(nspin.EQ.1)THEN
                  iss=1
                  IF (gstart == 2) vtemp(1)=0.0
                  DO ig=gstart,ngs
                     vcgs=CONJG(rhotemp(ig))*fpi/(tpiba2*g(ig))
                     cnvg=rhops(ig,is)*vcgs
                     cvn=vps(ig,is)*CONJG(rhog(ig,iss))
                     i = mill_l(1,ig)
                     j = mill_l(2,ig)
                     k = mill_l(3,ig)
                     eigrx=ei1(i,isa)*ei2(j,isa)*ei3(k,isa)
                     vtemp(ig)=eigrx*(cnvg+cvn)*CMPLX(0.d0,gx(ix,ig)) 
                  END DO
               ELSE
                  isup=1
                  isdw=2
                  IF (gstart == 2) vtemp(1)=0.0
                  DO ig=gstart,ngs
                     vcgs=CONJG(rhotemp(ig))*fpi/(tpiba2*g(ig))
                     cnvg=rhops(ig,is)*vcgs
                     cvn=vps(ig,is)*CONJG(rhog(ig,isup)                 &
     &                                   +rhog(ig,isdw))
                     i = mill_l(1,ig)
                     j = mill_l(2,ig)
                     k = mill_l(3,ig)
                     eigrx=ei1(i,isa)*ei2(j,isa)*ei3(k,isa)
                     vtemp(ig)=eigrx*(cnvg+cvn)*CMPLX(0.d0,gx(ix,ig)) 
                  END DO
               ENDIF
               fion1(ix,isa) = fion1(ix,isa) + tpiba*omega* wz*DBLE(SUM(vtemp))
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE force_ps
!
!-----------------------------------------------------------------------
      SUBROUTINE gausin(eigr,cm)
!-----------------------------------------------------------------------
!
! initialize wavefunctions with gaussians - edit to fit your system
!
      USE ions_base, ONLY: nas => nax, na, nsp, nat
      USE electrons_base, ONLY: n => nbsp
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gx, g
!
      IMPLICIT NONE
!
      COMPLEX(8) eigr(ngw,nat), cm(ngw,n)
      REAL(8)    sigma, auxf
      INTEGER nband, is, ia, ig, isa
!
      sigma=12.0
      nband=0
!!!      do is=1,nsp
      isa = 0
      is=1
         DO ia=1,na(is)
! s-like gaussians
            nband=nband+1
            DO ig=1,ngw
               auxf=EXP(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)
            END DO
! px-like gaussians
            nband=nband+1
            DO ig=1,ngw
               auxf=EXP(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)
            END DO
! py-like gaussians
            nband=nband+1
            DO ig=1,ngw
               auxf=EXP(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)
            END DO
! pz-like gaussians
            nband=nband+1
            DO ig=1,ngw
               auxf=EXP(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(3,ig)
            END DO
         END DO
      isa = isa + na(is)
      is=2
         DO ia=1,na(is)
! s-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)
!            end do
! px-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)
!            end do
! py-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)
!            end do
! pz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(2,ig)
!            end do
! dxz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)*gx(3,ig)
!            end do
! dx2-y2-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*                        &
!     &              (gx(1,ig)**2-gx(2,ig)**2)
!            end do
         END DO
!!!      end do
      RETURN
      END SUBROUTINE gausin
!            

!-------------------------------------------------------------------------
      SUBROUTINE gracsc( bec, nkbx, betae, cp, ngwx, i, csc, n )
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
      USE ions_base,      ONLY: na
      USE cvan,           ONLY :nvb, ish
      USE uspp,           ONLY : nkb, nhsavb=>nkbus, qq
      USE uspp_param,     ONLY:  nh
      USE electrons_base, ONLY: ispin
      USE gvecw,          ONLY: ngw
      USE mp,             ONLY: mp_sum
      USE kinds,          ONLY: DP
      USE reciprocal_vectors, ONLY: gstart
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: i, nkbx, ngwx, n
      COMPLEX(DP) :: betae( ngwx, nkb )
      REAL(DP)    :: bec( nkbx, n ), cp( 2, ngwx, n )
      REAL(DP)    :: csc( n )
      INTEGER     :: k, kmax,ig, is, iv, jv, ia, inl, jnl
      REAL(DP)    :: rsum
      REAL(DP), ALLOCATABLE :: temp(:) 

      !
      !     calculate csc(k)=<cp(i)|cp(k)>,  k<i
      !
      ALLOCATE( temp( ngw ) )

      kmax = i - 1

      DO k = 1, kmax
         csc(k) = 0.0d0
         IF ( ispin(i) .EQ. ispin(k) ) THEN
            DO ig = 1, ngw
               temp(ig) = cp(1,ig,k) * cp(1,ig,i) + cp(2,ig,k) * cp(2,ig,i)
            END DO
            csc(k) = 2.0d0 * SUM(temp)
            IF (gstart == 2) csc(k) = csc(k) - temp(1)
         ENDIF
      END DO

      CALL mp_sum( csc( 1:kmax ) )

      !
      !     calculate bec(i)=<cp(i)|beta>
      !
      DO inl=1,nhsavb
         DO ig=1,ngw
            temp(ig)=cp(1,ig,i)* DBLE(betae(ig,inl))+             &
     &               cp(2,ig,i)*AIMAG(betae(ig,inl))
         END DO
         bec(inl,i)=2.*SUM(temp)
         IF (gstart == 2) bec(inl,i)= bec(inl,i)-temp(1)
      END DO

      CALL mp_sum( bec( 1:nhsavb, i ) )
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
      DO k=1,kmax
         IF (ispin(i).EQ.ispin(k)) THEN
            rsum=0.
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN 
                        DO ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           rsum = rsum + qq(iv,jv,is)*bec(inl,i)*bec(jnl,k)
                        END DO
                     ENDIF
                  END DO
               END DO
            END DO
            csc(k)=csc(k)+rsum
         ENDIF
      END DO
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
!
      DO k=1,kmax
         DO inl=1,nhsavb
            bec(inl,i)=bec(inl,i)-csc(k)*bec(inl,k)
         END DO
      END DO

      DEALLOCATE( temp )
!
      RETURN
      END SUBROUTINE gracsc

!-------------------------------------------------------------------------
      SUBROUTINE gram( betae, bec, nkbx, cp, ngwx, n )
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the set of wavefunctions cp
!
      USE uspp,           ONLY : nkb, nhsavb=> nkbus
      USE gvecw,          ONLY : ngw
      USE kinds,          ONLY : DP
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx, n
      REAL(DP)      :: bec( nkbx, n )
      COMPLEX(DP)   :: cp( ngwx, n ), betae( ngwx, nkb )
!
      REAL(DP) :: anorm, cscnorm
      REAL(DP), ALLOCATABLE :: csc( : )
      INTEGER :: i,k
      EXTERNAL cscnorm
!
      CALL start_clock( 'gram' )

      ALLOCATE( csc( n ) )
!
      DO i = 1, n
         !
         CALL gracsc( bec, nkbx, betae, cp, ngwx, i, csc, n )
         !
         ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
         !
         DO k = 1, i - 1
            CALL DAXPY( 2*ngw, -csc(k), cp(1,k), 1, cp(1,i), 1 )
         END DO
         anorm = cscnorm( bec, nkbx, cp, ngwx, i, n )
         CALL DSCAL( 2*ngw, 1.0/anorm, cp(1,i), 1 )
         !
         !         these are the final bec's
         !
         CALL DSCAL( nhsavb, 1.0/anorm, bec(1,i), 1 )
      END DO
!
      DEALLOCATE( csc )

      CALL stop_clock( 'gram' )
!
      RETURN
      END SUBROUTINE gram
!
!-----------------------------------------------------------------------
      SUBROUTINE initbox ( tau0, taub, irb )
!-----------------------------------------------------------------------
!
!     sets the indexes irb and positions taub for the small boxes 
!     around atoms
!
      USE ions_base, ONLY: nsp, na, nat
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE cell_base, ONLY: ainv, a1, a2, a3
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
      USE control_flags, ONLY: iprsta
      USE io_global, ONLY: stdout
      USE mp_global, ONLY: nproc, mpime
      USE fft_base,  ONLY: dfftb, dfftp, fft_dlay_descriptor
      USE fft_types, ONLY: fft_box_set

      IMPLICIT NONE
! input
      REAL(8), INTENT(in):: tau0(3,nat)
! output
      INTEGER, INTENT(out):: irb(3,nat)
      REAL(8), INTENT(out):: taub(3,nat)
! local
      REAL(8) x(3), xmod
      INTEGER nr(3), nrb(3), xint, is, ia, i, isa
!
      nr (1)=nr1
      nr (2)=nr2
      nr (3)=nr3
      nrb(1)=nr1b
      nrb(2)=nr2b
      nrb(3)=nr3b
!
      isa = 0
      DO is=1,nsp
         DO ia=1,na(is)
           isa = isa + 1
!
            DO i=1,3
!
! bring atomic positions to crystal axis
!
               x(i) = ainv(i,1)*tau0(1,isa) +                         &
     &                ainv(i,2)*tau0(2,isa) +                         &
     &                ainv(i,3)*tau0(3,isa)
!
! bring x in the range between 0 and 1
!
               x(i) = MOD(x(i),1.d0)
               IF (x(i).LT.0.d0) x(i)=x(i)+1.d0
!
! case of nrb(i) even
!
               IF (MOD(nrb(i),2).EQ.0) THEN
!
! find irb = index of the grid point at the corner of the small box
!           (the indices of the small box run from irb to irb+nrb-1)
!
                  xint=INT(x(i)*nr(i))
                  irb (i,isa)=xint+1-nrb(i)/2+1
                  IF(irb(i,isa).LT.1) irb(i,isa)=irb(i,isa)+nr(i)
!
! x(i) are the atomic positions in crystal coordinates, where the
! "crystal lattice" is the small box lattice and the origin is at
! the corner of the small box. Used to calculate phases exp(iG*taub)
!
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+nrb(i)/2-1)/nr(i)
               ELSE
!
! case of nrb(i) odd - see above for comments
!
                  xint=NINT(x(i)*nr(i))
                  irb (i,isa)=xint+1-(nrb(i)-1)/2
                  IF(irb(i,isa).LT.1) irb(i,isa)=irb(i,isa)+nr(i)
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+(nrb(i)-1)/2)/nr(i)
               END IF
            END DO
!
! bring back taub in cartesian coordinates
!
            DO i=1,3
               taub(i,isa)= x(1)*a1(i) + x(2)*a2(i) + x(3)*a3(i)
            END DO
         END DO
      END DO

      CALL fft_box_set( dfftb, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, &
                        nat, irb, mpime+1, nproc, dfftp%npp, dfftp%ipp )

      IF( iprsta > 2 ) THEN
           isa = 1
           DO is=1,nsp
              WRITE( stdout,'(/,2x,''species= '',i2)') is
              DO ia=1,na(is)
                 WRITE( stdout,2000) ia, (irb(i,isa),i=1,3)
2000             FORMAT(2x,'atom= ',i3,' irb1= ',i3,' irb2= ',i3,' irb3= ',i3)
                 isa = isa + 1
               END DO
            END DO
      ENDIF

!
      RETURN
      END SUBROUTINE initbox
!
!-------------------------------------------------------------------------
      SUBROUTINE newd(vr,irb,eigrb,rhovan,fion)
!-----------------------------------------------------------------------
!     this routine calculates array deeq:
!         deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!     and the corresponding term in forces
!         fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!     where
!         rho_lm = \sum_j f_j <psi_j|beta_l><beta_m|psi_j>
!
      USE kinds, ONLY: dp
      USE uspp_param, ONLY: nh, nhm
      USE uspp, ONLY: deeq
      USE cvan, ONLY: nvb
      USE ions_base, ONLY: nas => nax, nat, nsp, na
      USE parameters, ONLY: natx, nsx
      USE constants, ONLY: pi, fpi
      USE grid_dimensions, ONLY: nr3, nnr => nnrx
      USE gvecb
      USE small_box, ONLY: omegab, tpibab
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE qgb_mod
      USE electrons_base, ONLY: nspin
      USE control_flags, ONLY: iprint, thdyn, tfor, tprnfor
      USE mp, ONLY: mp_sum
      USE fft_module, ONLY: invfft
      USE fft_base, ONLY: dfftb
!
      IMPLICIT NONE
! input
      INTEGER irb(3,nat)
      REAL(8) rhovan(nhm*(nhm+1)/2,nat,nspin)
      COMPLEX(8) eigrb(ngb,nat)
      REAL(8)  vr(nnr,nspin)
! output
      REAL(8)  fion(3,natx)
! local
      INTEGER isup,isdw,iss, iv,ijv,jv, ik, nfft, isa, ia, is, ig
      REAL(8)  fvan(3,natx,nsx), fac, fac1, fac2, boxdotgrid
      COMPLEX(8) ci, facg1, facg2
      COMPLEX(8), ALLOCATABLE :: qv(:)
      EXTERNAL boxdotgrid
!
      CALL start_clock( 'newd' )
      ci=(0.d0,1.d0)
      fac=omegab/DBLE(nr1b*nr2b*nr3b)
      deeq (:,:,:,:) = 0.d0
      fvan (:,:,:) = 0.d0

      ALLOCATE( qv( nnrb ) )
!
! calculation of deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!
      isa=1
      DO is=1,nvb
#ifdef __PARA
         DO ia=1,na(is)
            nfft=1
            IF ( dfftb%np3( isa ) <= 0 ) go to 15
#else
         DO ia=1,na(is),2
            nfft=2
#endif
            IF( ia .EQ. na(is) ) nfft=1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
            DO iv=1,nh(is)
               DO jv=iv,nh(is)
                  ijv = (jv-1)*jv/2 + iv
                  qv(:) = (0.d0, 0.d0)
                  IF (nfft.EQ.2) THEN
                     DO ig=1,ngb
                        qv(npb(ig))= eigrb(ig,isa  )*qgb(ig,ijv,is)   &
     &                          + ci*eigrb(ig,isa+1)*qgb(ig,ijv,is)
                        qv(nmb(ig))= CONJG(                             &
     &                               eigrb(ig,isa  )*qgb(ig,ijv,is))  &
     &                          + ci*CONJG(                             &
     &                               eigrb(ig,isa+1)*qgb(ig,ijv,is))
                     END DO
                  ELSE
                     DO ig=1,ngb
                        qv(npb(ig)) = eigrb(ig,isa)*qgb(ig,ijv,is)
                        qv(nmb(ig)) = CONJG(                            &
     &                                eigrb(ig,isa)*qgb(ig,ijv,is))
                     END DO
                  END IF
!
                  CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
                  DO iss=1,nspin
                     deeq(iv,jv,isa,iss) = fac *                        &
     &                    boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
                     IF (iv.NE.jv)                                      &
     &                    deeq(jv,iv,isa,iss)=deeq(iv,jv,isa,iss)
!
                     IF (nfft.EQ.2) THEN
                        deeq(iv,jv,isa+1,iss) = fac*                    &
     &                       boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
                        IF (iv.NE.jv)                                   &
     &                       deeq(jv,iv,isa+1,iss)=deeq(iv,jv,isa+1,iss)
                     END IF
                  END DO
               END DO
            END DO
  15        isa=isa+nfft
         END DO
      END DO

      CALL reduce(nat*nhm*nhm*nspin,deeq)

      IF (.NOT.( tfor .OR. thdyn .OR. tprnfor ) ) go to 10
!
! calculation of fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!
      isa=1
      IF(nspin.EQ.1) THEN
!     =================================================================
!     case nspin=1: two ffts at the same time, on two atoms (if possible)
!     -----------------------------------------------------------------
         iss=1
         isa=1
         DO is=1,nvb
#ifdef __PARA
            DO ia=1,na(is)
               nfft=1
               IF ( dfftb%np3( isa ) <= 0 ) go to 20
#else
            DO ia=1,na(is),2
               nfft=2
#endif
               IF( ia.EQ.na(is)) nfft=1
               DO ik=1,3
                  qv(:) = (0.d0, 0.d0)
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        IF(iv.NE.jv) THEN
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,iss)
                           IF (nfft.EQ.2) fac2=2.d0*fac*tpibab*         &
     &                                           rhovan(ijv,isa+1,iss)
                        ELSE
                           fac1=     fac*tpibab*rhovan(ijv,isa,iss)
                           IF (nfft.EQ.2) fac2=     fac*tpibab*        &
     &                                           rhovan(ijv,isa+1,iss)
                        ENDIF
                        IF (nfft.EQ.2) THEN
                           DO ig=1,ngb
                              facg1 = CMPLX(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is) * fac1
                              facg2 = CMPLX(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is) * fac2
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa  )*facg1  &
     &                                    + ci*eigrb(ig,isa+1)*facg2
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                                +   CONJG(eigrb(ig,isa  )*facg1)&
     &                                +ci*CONJG(eigrb(ig,isa+1)*facg2)
                           END DO
                        ELSE
                           DO ig=1,ngb
                              facg1 = CMPLX(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is)*fac1
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa)*facg1
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                               +  CONJG( eigrb(ig,isa)*facg1)
                           END DO
                        END IF
                     END DO
                  END DO
!
                  CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
                  fvan(ik,ia,is) =                                      &
     &                    boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
!
                  IF (nfft.EQ.2) fvan(ik,ia+1,is) =                     &
     &                    boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
               END DO
 20            isa = isa+nfft
            END DO
         END DO
      ELSE
!     =================================================================
!     case nspin=2: up and down spin fft's combined into a single fft
!     -----------------------------------------------------------------
         isup=1
         isdw=2
         isa=1
         DO is=1,nvb
            DO ia=1,na(is)
#ifdef __PARA
               IF ( dfftb%np3( isa ) <= 0 ) go to 25
#endif
               DO ik=1,3
                  qv(:) = (0.d0, 0.d0)
!
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        IF(iv.NE.jv) THEN
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=2.d0*fac*tpibab*rhovan(ijv,isa,isdw)
                        ELSE
                           fac1=     fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=     fac*tpibab*rhovan(ijv,isa,isdw)
                        END IF
                        DO ig=1,ngb
                           facg1 = fac1 * CMPLX(0.d0,-gxb(ik,ig)) *     &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           facg2 = fac2 * CMPLX(0.d0,-gxb(ik,ig)) *     &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           qv(npb(ig)) = qv(npb(ig))                    &
     &                                    + facg1 + ci*facg2
                           qv(nmb(ig)) = qv(nmb(ig))                    &
     &                                    +CONJG(facg1)+ci*CONJG(facg2)
                        END DO
                     END DO
                  END DO
!
                  CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
                  fvan(ik,ia,is) =                                      &
     &                    boxdotgrid(irb(1,isa),isup,qv,vr(1,isup)) + &
     &                    boxdotgrid(irb(1,isa),isdw,qv,vr(1,isdw))
               END DO
25             isa = isa+1
            END DO
         END DO
      END IF

      CALL reduce(3*natx*nvb,fvan)

      isa = 0
      DO is = 1, nvb
        DO ia = 1, na(is)
          isa = isa + 1
          fion(:,isa) = fion(:,isa) - fvan(:,ia,is)
        END DO
      END DO

      DEALLOCATE( qv )
!
  10  CALL stop_clock( 'newd' )
!
      RETURN
      END SUBROUTINE newd
!-------------------------------------------------------------------------
      SUBROUTINE nlfl(bec,becdr,lambda,fion)
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
! 
!
      USE io_global, ONLY: stdout
      USE ions_base, ONLY: na, nsp
      USE parameters, ONLY: natx
      USE uspp, ONLY :nhsa=>nkb, qq
      USE uspp_param, ONLY: nhm, nh
      USE cvan, ONLY: ish, nvb
      USE electrons_base, ONLY: nbspx, nbsp, nudx, nspin, iupdwn, nupdwn
      USE constants, ONLY: pi, fpi
!
      IMPLICIT NONE
      REAL(8) bec(nhsa,nbsp), becdr(nhsa,nbsp,3), lambda(nudx,nudx,nspin)
      REAL(8) fion(3,natx)
!
      INTEGER k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart
      REAL(8), ALLOCATABLE :: temp(:,:), tmpbec(:,:),tmpdr(:,:) 
      !
      CALL start_clock( 'nlfl' )
      !
      ALLOCATE( temp( nudx, nudx ), tmpbec( nhm, nudx ), tmpdr( nudx, nhm ) )

      DO k=1,3
         isa = 0
         DO is=1,nvb
            DO ia=1,na(is)
               isa = isa + 1
               !
               DO iss = 1, nspin
                  !
                  nss = nupdwn( iss )
                  istart = iupdwn( iss )
                  !
                  tmpbec = 0.d0
                  tmpdr  = 0.d0
!
                  DO iv=1,nh(is)
                     DO jv=1,nh(is)
                        inl=ish(is)+(jv-1)*na(is)+ia
                        IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                           DO i=1,nss
                              tmpbec(iv,i)=tmpbec(iv,i)                    &
     &                             + qq(iv,jv,is)*bec(inl,i+istart-1)
                           END DO
                        ENDIF
                     END DO
                  END DO
!
                  DO iv=1,nh(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     DO i=1,nss
                        tmpdr(i,iv)=becdr(inl,i+istart-1,k)
                     END DO
                  END DO
!
                  IF(nh(is).GT.0)THEN
                     !
                     temp = 0.d0
!
                     CALL MXMA                                             &
     &                 (tmpdr,1,nudx,tmpbec,1,nhm,temp,1,nudx,nss,nh(is),nss)
!
                     DO j=1,nss
                        DO i=1,nss
                           temp(i,j)=temp(i,j)*lambda(i,j,iss)
                        END DO
                     END DO
!
                     fion(k,isa)=fion(k,isa)+2.*SUM(temp)
                  ENDIF

               END DO
!
            END DO
         END DO
      END DO
      !
      !     end of x/y/z loop

      DEALLOCATE( temp, tmpbec, tmpdr )
      !
      CALL stop_clock( 'nlfl' )
      !
      RETURN
      END SUBROUTINE nlfl

!
!-----------------------------------------------------------------------
      SUBROUTINE pbc(rin,a1,a2,a3,ainv,rout)
!-----------------------------------------------------------------------
!
!     brings atoms inside the unit cell
!
      IMPLICIT NONE
! input
      REAL(8) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
! output
      REAL(8) rout(3)
! local
      REAL(8) x,y,z
!
! bring atomic positions to crystal axis
!
      x = ainv(1,1)*rin(1)+ainv(1,2)*rin(2)+ainv(1,3)*rin(3)
      y = ainv(2,1)*rin(1)+ainv(2,2)*rin(2)+ainv(2,3)*rin(3)
      z = ainv(3,1)*rin(1)+ainv(3,2)*rin(2)+ainv(3,3)*rin(3)
!
! bring x,y,z in the range between -0.5 and 0.5
!
      x = x - NINT(x)
      y = y - NINT(y)
      z = z - NINT(z)
!
! bring atomic positions back in cartesian axis
!
      rout(1) = x*a1(1)+y*a2(1)+z*a3(1)
      rout(2) = x*a1(2)+y*a2(2)+z*a3(2)
      rout(3) = x*a1(3)+y*a2(3)+z*a3(3)
!
      RETURN
      END SUBROUTINE pbc

!
!-------------------------------------------------------------------------
      SUBROUTINE prefor(eigr,betae)
!-----------------------------------------------------------------------
!
!     input :        eigr =  e^-ig.r_i
!     output:        betae_i,i(g) = (-i)**l beta_i,i(g) e^-ig.r_i 
!
      USE ions_base, ONLY: nas => nax, nsp, na, nat
      USE gvecw, ONLY: ngw
      USE cvan, ONLY: ish
      USE uspp, ONLY :nhsa=>nkb, beta, nhtol
      USE uspp_param, ONLY: nh
!
      IMPLICIT NONE
      COMPLEX(8) eigr(ngw,nat)
      COMPLEX(8) betae(ngw,nhsa)
!
      INTEGER is, iv, ia, inl, ig, isa
      COMPLEX(8) ci
!
      CALL start_clock( 'prefor' )
      isa = 0
      DO is=1,nsp
         DO iv=1,nh(is)
            ci=(0.,-1.)**nhtol(iv,is)
            DO ia=1,na(is)
               inl=ish(is)+(iv-1)*na(is)+ia
               DO ig=1,ngw
                  betae(ig,inl)=ci*beta(ig,iv,is)*eigr(ig,ia+isa)
               END DO
            END DO
         END DO
         isa = isa + na(is)
      END DO
      CALL stop_clock( 'prefor' )
!
      RETURN
      END SUBROUTINE prefor
!
!-----------------------------------------------------------------------
      SUBROUTINE projwfc(c,eigr,betae)
!-----------------------------------------------------------------------
!
! Projection on atomic wavefunctions
!
      USE io_global, ONLY: stdout
      USE electrons_base, ONLY: nx => nbspx, n => nbsp
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE ions_base, ONLY: nsp, na, nas => nax, nat
      USE uspp, ONLY: nhsa => nkb
      USE atom, ONLY: nchi, lchi
!
      IMPLICIT NONE
      COMPLEX(8), INTENT(in) :: c(ngw,nx), eigr(ngw,nat),      &
     &                               betae(ngw,nhsa)
!
      COMPLEX(8), ALLOCATABLE:: wfc(:,:), swfc(:,:), becwfc(:,:)
      REAL(8), ALLOCATABLE   :: overlap(:,:), e(:), z(:,:),        &
     &                               proj(:,:), temp(:)
      REAL(8)                :: somma
      INTEGER n_atomic_wfc
      INTEGER is, ia, nb, l, m, k, i
!
! calculate number of atomic states
!
      n_atomic_wfc=0
      DO is=1,nsp
         DO nb = 1,nchi(is)
            l = lchi(nb,is)
            n_atomic_wfc = n_atomic_wfc + (2*l+1)*na(is)
         END DO
      END DO
      IF (n_atomic_wfc.EQ.0) RETURN
!
      ALLOCATE(wfc(ngw,n_atomic_wfc))
!
! calculate wfc = atomic states
!
      CALL atomic_wfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
      ALLOCATE(becwfc(nhsa,n_atomic_wfc))
      CALL nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)

! calculate swfc = S|wfc>
!
      ALLOCATE(swfc(ngw,n_atomic_wfc))
      CALL s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate overlap(i,j) = <wfc_i|S|wfc_j> 
!
      ALLOCATE(overlap(n_atomic_wfc,n_atomic_wfc))
!
      CALL MXMA(wfc,2*ngw,1,swfc,1,2*ngw,overlap,1,                     &
     &          n_atomic_wfc,n_atomic_wfc,2*ngw,n_atomic_wfc)

      CALL reduce(n_atomic_wfc**2,overlap)

      overlap=overlap*2.d0
      IF (gstart == 2) THEN
         DO l=1,n_atomic_wfc
            DO m=1,n_atomic_wfc
               overlap(m,l)=overlap(m,l)-DBLE(wfc(1,m))*DBLE(swfc(1,l))
            END DO
         END DO
      END IF
!
! calculate (overlap)^(-1/2)(i,j). An orthonormal set of vectors |wfc_i>
! is obtained by introducing |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>
!
      ALLOCATE(z(n_atomic_wfc,n_atomic_wfc))
      ALLOCATE(e(n_atomic_wfc))
      CALL rdiag(n_atomic_wfc,overlap,n_atomic_wfc,e,z)
      overlap=0.d0
      DO l=1,n_atomic_wfc
         DO m=1,n_atomic_wfc
            DO k=1,n_atomic_wfc
               overlap(l,m)=overlap(l,m)+z(m,k)*z(l,k)/SQRT(e(k))
            END DO
         END DO
      END DO
      DEALLOCATE(e)
      DEALLOCATE(z)
!
! calculate |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>   (note the S matrix!)
!
      wfc=0.d0
      DO m=1,n_atomic_wfc
         DO l=1,n_atomic_wfc
            wfc(:,m)=wfc(:,m)+overlap(l,m)*swfc(:,l)
         END DO
      END DO
      DEALLOCATE(overlap)
      DEALLOCATE(swfc)
      DEALLOCATE(becwfc)
!
! calculate proj = <c|S|wfc> 
!
      ALLOCATE(proj(n,n_atomic_wfc))
      ALLOCATE(temp(ngw))
      DO m=1,n
         DO l=1,n_atomic_wfc
            temp(:)=DBLE(CONJG(c(:,m))*wfc(:,l))
            proj(m,l)=2.d0*SUM(temp)
            IF (gstart == 2) proj(m,l)=proj(m,l)-temp(1)
         END DO
      END DO
      DEALLOCATE(temp)

      CALL reduce(n*n_atomic_wfc,proj)

      i=0
      WRITE( stdout,'(/''Projection on atomic states:'')')
      DO is=1,nsp
         DO nb = 1,nchi(is)
            l=lchi(nb,is)
            DO m = -l,l
               DO ia=1,na(is)
                  i=i+1
                  WRITE( stdout,'(''atomic state # '',i3,'': atom # '',i3,    &
     &                      ''  species # '',i2,''  wfc # '',i2,        &
     &                      '' (l='',i1,'' m='',i2,'')'')')             &
     &                 i, ia, is, nb, l, m
               END DO
            END DO
         END DO
      END DO

      WRITE( stdout,*)
      DO m=1,n
         somma=0.d0
         DO l=1,n_atomic_wfc
            somma=somma+proj(m,l)**2
         END DO
         WRITE( stdout,'(''state # '',i4,''    sum c^2 ='',f7.4)') m,somma
         WRITE( stdout,'(10f7.4)') (ABS(proj(m,l)),l=1,n_atomic_wfc)
      END DO
!
      DEALLOCATE(proj)
      DEALLOCATE(wfc)
      RETURN
      END SUBROUTINE projwfc
!---------------------------------------------------------------------
      SUBROUTINE randin(nmin,nmax,gstart,ngw,ampre,c)
!---------------------------------------------------------------------
!
      USE wave_functions, ONLY: wave_rand_init
      IMPLICIT NONE

! input
      INTEGER nmin, nmax, gstart, ngw
      REAL(8) ampre
! output
      COMPLEX(8) c(ngw,nmax)
! local
      INTEGER i,j
      REAL(8) ranf1, randy, ranf2, ampexp
!
      CALL wave_rand_init( c )
!      do i=nmin,nmax
!         do j=gstart,ngw
!            ranf1=.5-randy()
!            ranf2=.5-randy()
!            ampexp=ampre*exp(-(4.*j)/ngw)
!            c(j,i)=c(j,i)+ampexp*CMPLX(ranf1,ranf2)
!         end do
!      end do
!
      RETURN
      END SUBROUTINE randin
!
!-----------------------------------------------------------------------
      SUBROUTINE rdiag (n,h,ldh,e,v)
!-----------------------------------------------------------------------
!
!   calculates all the eigenvalues and eigenvectors of a complex
!   hermitean matrix H . On output, the matrix H is destroyed
!
      USE kinds,            ONLY: DP
      USE parallel_toolkit, ONLY: zhpev_drv
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in)           :: n, ldh
      COMPLEX(DP), INTENT(inout):: h(ldh,n)
      REAL   (DP), INTENT(out)  :: e(n)
      COMPLEX(DP), INTENT(out)  :: v(ldh,n)
!
      INTEGER :: i, j, k
      COMPLEX(DP), ALLOCATABLE :: ap( : )
!
      ALLOCATE( ap( n * ( n + 1 ) / 2 ) )

      K = 0
      DO J = 1, n
         DO I = J, n
            K = K + 1
            ap( k ) = h( i, j )
         END DO
      END DO

      CALL zhpev_drv( 'V', 'L', n, ap, e, v, ldh )

      DEALLOCATE( ap )
!
      RETURN
      END SUBROUTINE rdiag

!
!-----------------------------------------------------------------------
      SUBROUTINE rhov(irb,eigrb,rhovan,rhog,rhor)
!-----------------------------------------------------------------------
!     Add Vanderbilt contribution to rho(r) and rho(g)
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      USE kinds, ONLY: dp
      USE ions_base, ONLY: nas => nax, nat, na, nsp
      USE io_global, ONLY: stdout
      USE cvan, ONLY: nvb
      USE uspp_param, ONLY: nh, nhm
      USE uspp, ONLY: deeq
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      USE electrons_base, ONLY: nspin
      USE gvecb
      USE gvecp, ONLY: ng => ngm
      USE cell_base, ONLY: omega
      USE small_box, ONLY: omegab
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE control_flags, ONLY: iprint, iprsta
      USE qgb_mod
      USE recvecs_indexes, ONLY: np, nm
      USE fft_module, ONLY: fwfft, invfft
      USE fft_base, ONLY: dfftb
!
      IMPLICIT NONE
!
      REAL(8) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      INTEGER, INTENT(in) :: irb(3,nat)
      COMPLEX(8), INTENT(in):: eigrb(ngb,nat)
      REAL(8), INTENT(inout):: rhor(nnr,nspin)
      COMPLEX(8),  INTENT(inout):: rhog(ng,nspin)
!
      INTEGER isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss,           &
     &     isa, ia, ir
      REAL(8) sumrho
      COMPLEX(8) ci, fp, fm, ca
      COMPLEX(8), ALLOCATABLE::  qgbt(:,:)
      COMPLEX(8), ALLOCATABLE:: v(:)
      COMPLEX(8), ALLOCATABLE:: qv(:)
!
      IF (nvb.EQ.0) RETURN
      CALL start_clock( 'rhov' )
      ci=(0.,1.)
!
!
      ALLOCATE( v( nnr ) )
      ALLOCATE( qv( nnrb ) )
      v (:) = (0.d0, 0.d0)
      ALLOCATE( qgbt( ngb, 2 ) )

!
      IF(nspin.EQ.1) THEN
         ! 
         !     nspin=1 : two fft at a time, one per atom, if possible
         !
         iss=1
         isa=1

         DO is = 1, nvb

#ifdef __PARA

            DO ia = 1, na(is)
               nfft = 1
               IF ( dfftb%np3( isa ) <= 0 ) go to 15
#else

            DO ia = 1, na(is), 2
               nfft = 2
#endif

               IF( ia .EQ. na(is) ) nfft = 1

               !
               !  nfft=2 if two ffts at the same time are performed
               !
               DO ifft=1,nfft
                  qgbt(:,ifft) = (0.d0, 0.d0)
                  DO iv= 1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        sumrho=rhovan(ijv,isa+ifft-1,iss)
                        IF(iv.NE.jv) sumrho=2.*sumrho
                        DO ig=1,ngb
                           qgbt(ig,ifft)=qgbt(ig,ifft) + sumrho*qgb(ig,ijv,is)
                        END DO
                     END DO
                  END DO
               END DO
               !
               ! add structure factor
               !
               qv(:) = (0.d0, 0.d0)
               IF(nfft.EQ.2)THEN
                  DO ig=1,ngb
                     qv(npb(ig))=  &
                                   eigrb(ig,isa  )*qgbt(ig,1)  &
                        + ci*      eigrb(ig,isa+1)*qgbt(ig,2)
                     qv(nmb(ig))=                                       &
                             CONJG(eigrb(ig,isa  )*qgbt(ig,1))        &
                        + ci*CONJG(eigrb(ig,isa+1)*qgbt(ig,2))
                  END DO
               ELSE
                  DO ig=1,ngb
                     qv(npb(ig)) = eigrb(ig,isa)*qgbt(ig,1)
                     qv(nmb(ig)) = CONJG(eigrb(ig,isa)*qgbt(ig,1))
                  END DO
               ENDIF

               CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)

               !
               !  qv = US augmentation charge in real space on box grid
               !       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
 
               IF(iprsta.GT.2) THEN
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom g-sp = ',         &
     &                 omegab*DBLE(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom r-sp = ',         &
     &                 omegab*DBLE(ca)/(nr1b*nr2b*nr3b)
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom g-sp = ',         &
     &                 omegab*DBLE(qgbt(1,2))
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom r-sp = ',         &
     &                 omegab*AIMAG(ca)/(nr1b*nr2b*nr3b)
               ENDIF
               !
               !  add qv(r) to v(r), in real space on the dense grid
               !
               CALL  box2grid(irb(1,isa),1,qv,v)
               IF (nfft.EQ.2) CALL  box2grid(irb(1,isa+1),2,qv,v)
  15           isa=isa+nfft
!
            END DO
         END DO
         !
         !  rhor(r) = total (smooth + US) charge density in real space
         !
         DO ir=1,nnr
            rhor(ir,iss)=rhor(ir,iss)+DBLE(v(ir))        
         END DO
!
         IF(iprsta.GT.2) THEN
            ca = SUM(v)

            CALL reduce(2,ca)

            WRITE( stdout,'(a,2f12.8)')                                  &
     &           ' rhov: int  n_v(r)  dr = ',omega*ca/(nr1*nr2*nr3)
         ENDIF
!
         CALL fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         IF(iprsta.GT.2) THEN
            WRITE( stdout,*) ' rhov: smooth ',omega*rhog(1,iss)
            WRITE( stdout,*) ' rhov: vander ',omega*v(1)
            WRITE( stdout,*) ' rhov: all    ',omega*(rhog(1,iss)+v(1))
         ENDIF
         !
         !  rhog(g) = total (smooth + US) charge density in G-space
         !
         DO ig=1,ng
            rhog(ig,iss)=rhog(ig,iss)+v(np(ig))
         END DO
!
         IF(iprsta.GT.1) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) = ',omega*DBLE(rhog(1,iss))
!
      ELSE
         !
         !     nspin=2: two fft at a time, one for spin up and one for spin down
         !
         isup=1
         isdw=2
         isa=1
         DO is=1,nvb
            DO ia=1,na(is)
#ifdef __PARA
               IF ( dfftb%np3( isa ) <= 0 ) go to 25
#endif
               DO iss=1,2
                  qgbt(:,iss) = (0.d0, 0.d0)
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        sumrho=rhovan(ijv,isa,iss)
                        IF(iv.NE.jv) sumrho=2.*sumrho
                        DO ig=1,ngb
                           qgbt(ig,iss)=qgbt(ig,iss)+sumrho*qgb(ig,ijv,is)
                        END DO
                     END DO
                  END DO
               END DO
!     
! add structure factor
!
               qv(:) = (0.d0, 0.d0)
               DO ig=1,ngb
                  qv(npb(ig)) =    eigrb(ig,isa)*qgbt(ig,1)           &
     &                  + ci*      eigrb(ig,isa)*qgbt(ig,2)
                  qv(nmb(ig)) = CONJG(eigrb(ig,isa)*qgbt(ig,1))       &
     &                  + ci*   CONJG(eigrb(ig,isa)*qgbt(ig,2))
               END DO
!
               CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
!
!  qv is the now the US augmentation charge for atomic species is
!  and atom ia: real(qv)=spin up, imag(qv)=spin down
!
               IF(iprsta.GT.2) THEN
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: up   g-space = ',        &
     &                 omegab*DBLE(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: up r-sp = ',             &
     &                 omegab*DBLE(ca)/(nr1b*nr2b*nr3b)
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw g-space = ',          &
     &                 omegab*DBLE(qgbt(1,2))
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw r-sp = ',             &
     &                 omegab*AIMAG(ca)/(nr1b*nr2b*nr3b)
               ENDIF
!
!  add qv(r) to v(r), in real space on the dense grid
!
               CALL box2grid2(irb(1,isa),qv,v)
  25           isa=isa+1
!
            END DO
         END DO
!
         DO ir=1,nnr
            rhor(ir,isup)=rhor(ir,isup)+DBLE(v(ir)) 
            rhor(ir,isdw)=rhor(ir,isdw)+AIMAG(v(ir)) 
         END DO
!
         IF(iprsta.GT.2) THEN
            ca = SUM(v)
            CALL reduce(2,ca)
            WRITE( stdout,'(a,2f12.8)') 'rhov:in n_v  ',omega*ca/(nr1*nr2*nr3)
         ENDIF
!
         CALL fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         IF(iprsta.GT.2) THEN
            WRITE( stdout,*) 'rhov: smooth up',omega*rhog(1,isup)
            WRITE( stdout,*) 'rhov: smooth dw',omega*rhog(1,isdw)
            WRITE( stdout,*) 'rhov: vander up',omega*DBLE(v(1))
            WRITE( stdout,*) 'rhov: vander dw',omega*AIMAG(v(1))
            WRITE( stdout,*) 'rhov: all up',                                  &
     &           omega*(rhog(1,isup)+DBLE(v(1)))
            WRITE( stdout,*) 'rhov: all dw',                                  &
     &           omega*(rhog(1,isdw)+AIMAG(v(1)))
         ENDIF
!
         DO ig=1,ng
            fp=  v(np(ig)) + v(nm(ig))
            fm=  v(np(ig)) - v(nm(ig))
            rhog(ig,isup)=rhog(ig,isup) + 0.5*CMPLX(DBLE(fp),AIMAG(fm))
            rhog(ig,isdw)=rhog(ig,isdw) + 0.5*CMPLX(AIMAG(fp),-DBLE(fm))
         END DO
!
         IF(iprsta.GT.2) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) up   = ',omega*DBLE (rhog(1,isup))
         IF(iprsta.GT.2) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) down = ',omega*DBLE(rhog(1,isdw))
!
      ENDIF

      DEALLOCATE(qgbt)
      DEALLOCATE( v )
      DEALLOCATE( qv )

      CALL stop_clock( 'rhov' )
!
      RETURN
      END SUBROUTINE rhov
!
!
!-------------------------------------------------------------------------
      SUBROUTINE s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
      USE ions_base, ONLY: na
      USE cvan, ONLY: nvb, ish
      USE uspp, ONLY: nhsa => nkb, nhsavb=>nkbus, qq
      USE uspp_param, ONLY: nh
      USE gvecw, ONLY: ngw
      !use parm
      USE constants, ONLY: pi, fpi
      IMPLICIT NONE
! input
      INTEGER, INTENT(in)         :: n_atomic_wfc
      COMPLEX(8), INTENT(in) :: betae(ngw,nhsa),                   &
     &                               wfc(ngw,n_atomic_wfc)
      REAL(8), INTENT(in)    :: becwfc(nhsa,n_atomic_wfc)
! output
      COMPLEX(8), INTENT(out):: swfc(ngw,n_atomic_wfc)
! local
      INTEGER is, iv, jv, ia, inl, jnl, i
      REAL(8) qtemp(nhsavb,n_atomic_wfc)
!
      swfc=0.d0
!
      IF (nvb.GT.0) THEN
         qtemp=0.d0
         DO is=1,nvb
            DO iv=1,nh(is)
               DO jv=1,nh(is)
                  IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                     DO ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        DO i=1,n_atomic_wfc
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*becwfc(jnl,i)
                        END DO
                     END DO
                  ENDIF
               END DO
            END DO
         END DO
!
         CALL MXMA (betae,1,2*ngw,qtemp,1,nhsavb,swfc,1,                &
     &              2*ngw,2*ngw,nhsavb,n_atomic_wfc)
      END IF
!
      swfc=swfc+wfc
!
      RETURN
      END SUBROUTINE s_wfc

!-----------------------------------------------------------------------
      SUBROUTINE spinsq (c,bec,rhor)
!-----------------------------------------------------------------------
!
!     estimate of <S^2>=s(s+1) in two different ways.
!     1) using as many-body wavefunction a single Slater determinant
!        constructed with Kohn-Sham orbitals:
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw - 
!                \sum_up\sum_dw < psi_up | psi_dw >
!
!        where Nup, Ndw = number of up and down states, the sum is over 
!        occupied states. Not suitable for fractionary occupancy.
!        In the ultrasoft scheme (c is the smooth part of \psi): 
!
!        < psi_up | psi_dw > = \sum_G c*_up(G) c_dw(G) +
!                              \int Q_ij <c_up|beta_i><beta_j|c_dw>
!
!        This is the usual formula, unsuitable for fractionary occupancy.
!     2) using the "LSD model" of Wang, Becke, Smith, JCP 102, 3477 (1995):
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw -
!                \int max(rhoup(r),rhodw(r)) dr
!
!     Requires on input: c=psi, bec=<c|beta>, rhoup(r), rhodw(r)
!     Assumes real psi, with only half G vectors.
!
      USE electrons_base, ONLY: nx => nbspx, n => nbsp, iupdwn, nupdwn, f, nel, nspin
      USE io_global, ONLY: stdout
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nnr => nnrx
      USE cell_base, ONLY: omega
      USE cvan, ONLY: nvb, ish
      USE uspp, ONLY: nhsa => nkb, nhsavb=>nkbus, qq
      USE uspp_param, ONLY: nh
      USE ions_base, ONLY: na
!
      IMPLICIT NONE
! input
      REAL(8) bec(nhsa,n), rhor(nnr,nspin)
      COMPLEX(8) c(ngw,nx)
! local variables
      INTEGER nup, ndw, ir, i, j, jj, ig, ia, is, iv, jv, inl, jnl
      REAL(8) spin0, spin1, spin2, fup, fdw
      REAL(8), ALLOCATABLE:: overlap(:,:), temp(:)
      LOGICAL frac
!
!
      IF (nspin.EQ.1) RETURN
!
! find spin-up and spin-down states
!
      fup = 0.0
      DO i=iupdwn(1),nupdwn(1)
         fup = fup + f(i)
      END DO
      nup = NINT(fup)
      ndw = nel(1)+nel(2) - nup
!
! paranoid checks
!
      frac= ABS(fup-nup).GT.1.0e-6
      fup = 0.0
      DO i=1,nup
         fup = fup + f(i)
      END DO
      frac=frac.OR.ABS(fup-nup).GT.1.0e-6
      fdw = 0.0
      DO j=iupdwn(2),iupdwn(2)-1+ndw
         fdw = fdw + f(j)
      END DO
      frac=frac.OR.ABS(fdw-ndw).GT.1.0e-6
!
      spin0 = ABS(fup-fdw)/2.d0 * ( ABS(fup-fdw)/2.d0 + 1.d0 ) + fdw
!
!     Becke's formula for spin polarization
!
      spin1 = 0.0
      DO ir=1,nnr
         spin1 = spin1 - MIN(rhor(ir,1),rhor(ir,2))
      END DO
      CALL reduce(1,spin1)
      spin1 = spin0 + omega/(nr1*nr2*nr3)*spin1
      IF (frac) THEN
         WRITE( stdout,'(/'' Spin contamination: s(s+1)='',f5.2,'' (Becke) '',&
     &                             f5.2,'' (expected)'')')              &
     &          spin1, ABS(fup-fdw)/2.d0*(ABS(fup-fdw)/2.d0+1.d0)
         RETURN
      END IF
!
!     Slater formula, smooth contribution to  < psi_up | psi_dw >
!
      ALLOCATE (overlap(nup,ndw))
      ALLOCATE (temp(ngw))
      DO j=1,ndw
         jj=j+iupdwn(2)-1
         DO i=1,nup
            overlap(i,j)=0.d0
            DO ig=1,ngw
               temp(ig)=DBLE(CONJG(c(ig,i))*c(ig,jj))
            END DO
            overlap(i,j) = 2.d0*SUM(temp)
            IF (gstart == 2) overlap(i,j) = overlap(i,j) - temp(1)
         END DO
      END DO
      DEALLOCATE (temp)
      CALL reduce(nup*ndw,overlap)
      DO j=1,ndw
         jj=j+iupdwn(2)-1
         DO i=1,nup
!
!     vanderbilt contribution to  < psi_up | psi_dw >
!
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN 
                        DO ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           overlap(i,j) = overlap(i,j) +                &
     &                          qq(iv,jv,is)*bec(inl,i)*bec(jnl,jj)
                        END DO
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
!
      spin2 = spin0
      DO j=1,ndw
         DO i=1,nup
            spin2 = spin2 - overlap(i,j)**2
         END DO
      END DO
!
      DEALLOCATE (overlap)
!
      WRITE( stdout,'(/" Spin contamination: s(s+1)=",f5.2," (Slater) ",  &
     &          f5.2," (Becke) ",f5.2," (expected)")')              &
     &     spin2,spin1, ABS(fup-fdw)/2.d0*(ABS(fup-fdw)/2.d0+1.d0)
!
      RETURN
      END SUBROUTINE spinsq

!
!-----------------------------------------------------------------------
      SUBROUTINE vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast,           &
     &     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot,
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!
      USE kinds, ONLY: dp
      USE control_flags, ONLY: iprint, iprsta, thdyn, tpre, tfor, tprnfor
      USE io_global, ONLY: stdout
      USE parameters, ONLY: natx
      USE ions_base, ONLY: nas => nax, nsp, na, nat
      USE gvecs
      USE gvecp, ONLY: ng => ngm
      USE cell_base, ONLY: omega
      USE cell_base, ONLY: a1, a2, a3, tpiba2
      USE reciprocal_vectors, ONLY: gstart, g
      USE recvecs_indexes, ONLY: np, nm
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      USE electrons_base, ONLY: nspin
      USE constants, ONLY: pi, fpi
      USE energies, ONLY: etot, eself, enl, ekin, epseu, esr, eht, exc 
      USE local_pseudo, ONLY: vps, rhops
      USE core, ONLY: nlcc_any
      USE gvecb
      USE dener
      USE derho
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: group
      USE funct, ONLY: dft_is_meta
      USE fft_module, ONLY: fwfft, invfft
      USE sic_module, ONLY: self_interaction, sic_epsilon, sic_alpha
      USE energies,   ONLY: self_sxc, self_ehte
!
      IMPLICIT NONE
!
      LOGICAL :: tlast, tfirst
      INTEGER :: nfi
      REAL(DP)  rhor(nnr,nspin), rhos(nnrsx,nspin), fion(3,natx)
      REAL(DP)  rhoc(nnr), tau0(3,natx)
      COMPLEX(DP) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),     &
     &                ei3(-nr3:nr3,nat), eigrb(ngb,nat),        &
     &                rhog(ng,nspin), sfac(ngs,nsp)
!
      INTEGER irb(3,nat), iss, isup, isdw, ig, ir,i,j,k,is, ia
      REAL(DP) fion1(3,natx), vave, ebac, wz, eh
      COMPLEX(DP)  fp, fm, ci
      COMPLEX(DP), ALLOCATABLE :: v(:), vs(:)
      COMPLEX(DP), ALLOCATABLE :: rhotmp(:), vtemp(:), drhotmp(:,:,:)
!
      COMPLEX(DP), ALLOCATABLE :: self_vloc(:)
      COMPLEX(DP)              :: self_rhoeg
      REAL(DP)                 :: self_ehtet , fpibg
      LOGICAL                  :: ttsic
!

      CALL start_clock( 'vofrho' )
      ci=(0.,1.)
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz = 2.0
      ALLOCATE( v( nnr ) )
      ALLOCATE( vs( nnrsx ) )
      ALLOCATE(vtemp(ng))
      ALLOCATE(rhotmp(ng))
      IF (tpre) ALLOCATE(drhotmp(ng,3,3))
!
      ttsic = ( ABS(self_interaction) /= 0 )
      IF( tpre .AND. ttsic ) &
      &   CALL errore('in cplib tpre and ttsic are not possible', 1)

      IF( ttsic ) ALLOCATE(self_vloc(ng))
!
!     first routine in which fion is calculated: annihilation
!
      fion =0.d0
      fion1=0.d0
!
!     ===================================================================
!     forces on ions, ionic term in real space
!     -------------------------------------------------------------------
      IF( tprnfor .OR. tfor .OR. tfirst .OR. tpre ) THEN
        CALL force_ion( tau0, esr, fion, dsr )
      END IF
!

      IF( nspin == 1 ) THEN
         iss=1
         DO ig=1,ng
            rhotmp(ig)=rhog(ig,iss)
         END DO
         IF(tpre)THEN
            DO j=1,3
               DO i=1,3
                  DO ig=1,ng
                     drhotmp(ig,i,j)=drhog(ig,iss,i,j)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ELSE
         isup=1
         isdw=2
         DO ig=1,ng
            rhotmp(ig)=rhog(ig,isup)+rhog(ig,isdw)
         END DO
         IF(tpre)THEN
            DO i=1,3
               DO j=1,3
                  DO ig=1,ng
                     drhotmp(ig,i,j) = drhog(ig,isup,i,j) +           &
     &                                 drhog(ig,isdw,i,j)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      END IF
!     ===================================================================
!     calculation local potential energy
!     -------------------------------------------------------------------
      vtemp=(0.,0.)
      DO is=1,nsp
         DO ig=1,ngs
            vtemp(ig)=vtemp(ig)+CONJG(rhotmp(ig))*sfac(ig,is)*vps(ig,is)
         END DO
      END DO

      epseu=wz*DBLE(SUM(vtemp))
      IF (gstart == 2) epseu=epseu-vtemp(1)
      CALL mp_sum( epseu, group )
      epseu=epseu*omega
!
      IF(tpre) CALL denps(rhotmp,drhotmp,sfac,vtemp,dps)
!
!     ===================================================================
!     calculation hartree energy
!     -------------------------------------------------------------------
     !
      self_ehtet = 0.d0
      IF( ALLOCATED( self_vloc ) ) self_vloc = 0.d0

      DO is=1,nsp
         DO ig=1,ngs
            rhotmp(ig)=rhotmp(ig)+sfac(ig,is)*rhops(ig,is)
         END DO
      END DO
      IF (gstart == 2) vtemp(1)=0.0
      DO ig=gstart,ng
         vtemp(ig)=CONJG(rhotmp(ig))*rhotmp(ig)/g(ig)
      END DO
!
      eh=DBLE(SUM(vtemp))*wz*0.5*fpi/tpiba2
!
      IF ( ttsic ) THEN
       DO ig = gstart,ng
!
         fpibg = fpi /(tpiba2 *g(ig))
!
         self_rhoeg    = rhog(ig,1) - rhog(ig,2)
         self_ehtet    = self_ehtet +  fpibg * DBLE(self_rhoeg * CONJG(self_rhoeg))
         self_vloc(ig) = sic_epsilon * fpibg * self_rhoeg
!
       ENDDO

        IF(gstart == 2) self_vloc(1) = 0.D0
         self_ehte = sic_epsilon * self_ehtet * wz * 0.5d0
         eh = eh - self_ehte

         CALL mp_sum( self_ehte, group )

      ENDIF
!
      CALL mp_sum( eh, group )
      !
      IF(tpre) CALL denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
      IF(tpre) DEALLOCATE(drhotmp)

!     ===================================================================
!     forces on ions, ionic term in reciprocal space
!     -------------------------------------------------------------------
      IF( tprnfor .OR. tfor .OR. tpre)                                                  &
     &    CALL force_ps(rhotmp,rhog,vtemp,ei1,ei2,ei3,fion1)
!     ===================================================================
!     calculation hartree + local pseudo potential
!     -------------------------------------------------------------------
!
      IF (gstart == 2) vtemp(1)=(0.,0.)
      DO ig=gstart,ng
         vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      END DO
!
      DO is=1,nsp
         DO ig=1,ngs
            vtemp(ig)=vtemp(ig)+sfac(ig,is)*vps(ig,is)
         END DO
      END DO
!
!     vtemp = v_loc(g) + v_h(g)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
      IF (nlcc_any) CALL add_cc(rhoc,rhog,rhor)
!
      CALL exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_sxc )

!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
      IF( nspin == 1 ) THEN
         iss = 1
         DO ir = 1, nnr
            v(ir) = CMPLX( rhor( ir, iss ), 0.d0 )
         END DO
         !
         !     v_xc(r) --> v_xc(g)
         !
         CALL fwfft( 'Dense', v, nr1, nr2, nr3, nr1x, nr2x, nr3x )
!
         DO ig = 1, ng
            rhog( ig, iss ) = vtemp(ig) + v( np( ig ) )
         END DO
         !
         !     v_tot(g) = (v_tot(g) - v_xc(g)) +v_xc(g)
         !     rhog contains the total potential in g-space
         !
      ELSE
         isup=1
         isdw=2
         DO ir=1,nnr
            v(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw))
         END DO
         CALL fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         DO ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            IF( ttsic ) THEN
             rhog(ig,isup)=vtemp(ig)-self_vloc(ig) +0.5*CMPLX( DBLE(fp),AIMAG(fm))
             rhog(ig,isdw)=vtemp(ig)+self_vloc(ig) +0.5*CMPLX(AIMAG(fp),-DBLE(fm))
            ELSE
             rhog(ig,isup)=vtemp(ig)+0.5*CMPLX( DBLE(fp),AIMAG(fm))
             rhog(ig,isdw)=vtemp(ig)+0.5*CMPLX(AIMAG(fp),-DBLE(fm))
            ENDIF
         END DO
      ENDIF

!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
      IF( tprnfor .OR. tfor ) THEN

         IF ( nlcc_any ) CALL force_cc( irb, eigrb, rhor, fion1 )

         CALL mp_sum( fion1 )
         !
         !    add g-space ionic and core correction contributions to fion
         !
         fion = fion + fion1

      END IF
!
      IF( ALLOCATED( self_vloc ) ) DEALLOCATE( self_vloc )
!
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
      v(:) = (0.d0, 0.d0)
      IF(nspin.EQ.1) THEN
         iss=1
         DO ig=1,ng
            v(np(ig))=rhog(ig,iss)
            v(nm(ig))=CONJG(rhog(ig,iss))
         END DO
!
!     v(g) --> v(r)
!
         CALL invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         DO ir=1,nnr
            rhor(ir,iss)=DBLE(v(ir))
         END DO
!
!     calculation of average potential
!
         vave=SUM(rhor(:,iss))/DBLE(nr1*nr2*nr3)
      ELSE
         isup=1
         isdw=2
         DO ig=1,ng
            v(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            v(nm(ig))=CONJG(rhog(ig,isup)) +ci*CONJG(rhog(ig,isdw))
         END DO
!
         CALL invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         DO ir=1,nnr
            rhor(ir,isup)= DBLE(v(ir))
            rhor(ir,isdw)=AIMAG(v(ir))
         END DO
!
!     calculation of average potential
!
         vave=(SUM(rhor(:,isup))+SUM(rhor(:,isdw)))       &
     &        /2.0/DBLE(nr1*nr2*nr3)
      ENDIF

      CALL mp_sum( vave, group )

      !
      !     fourier transform of total potential to r-space (smooth grid)
      !
      vs (:) = (0.d0, 0.d0)
      !
      IF(nspin.EQ.1)THEN
         iss=1
         DO ig=1,ngs
            vs(nms(ig))=CONJG(rhog(ig,iss))
            vs(nps(ig))=rhog(ig,iss)
         END DO
!
         CALL invfft('Smooth',vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
         DO ir=1,nnrsx
            rhos(ir,iss)=DBLE(vs(ir))
         END DO
      ELSE
         isup=1
         isdw=2
         DO ig=1,ngs
            vs(nps(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            vs(nms(ig))=CONJG(rhog(ig,isup)) +ci*CONJG(rhog(ig,isdw))
         END DO 
         CALL invfft('Smooth',vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
         DO ir=1,nnrsx
            rhos(ir,isup)= DBLE(vs(ir))
            rhos(ir,isdw)=AIMAG(vs(ir))
         END DO
      ENDIF

      IF( dft_is_meta() ) CALL vofrho_meta( v, vs )  !METAGGA

      ebac=0.0
!
      eht=eh*omega+esr-eself
!
!     etot is the total energy ; ekin, enl were calculated in rhoofr
!
      etot=ekin+eht+epseu+enl+exc+ebac
      IF(tpre) detot=dekin+dh+dps+denl+dxc+dsr
!
      DEALLOCATE(rhotmp)
      DEALLOCATE(vtemp)
      DEALLOCATE( v )
      DEALLOCATE( vs )
!
      CALL stop_clock( 'vofrho' )

      IF( ( nfi == 0 ) .OR. tfirst .OR. tlast ) GOTO 999
      IF( MOD( nfi - 1, iprint) /= 0 ) RETURN
!
 999  IF ( tpre ) THEN
         IF( iprsta >= 2 ) THEN  
            WRITE( stdout,*)
            WRITE( stdout,*) "From vofrho:"
            WRITE( stdout,*) "cell parameters h"
            WRITE( stdout,5555) (a1(i),a2(i),a3(i),i=1,3)
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(tot)"
            WRITE( stdout,5555) ((detot(i,j),j=1,3),i=1,3)
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(kin)"
            WRITE( stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(electrostatic)"
            WRITE( stdout,5555) (((dh(i,j)+dsr(i,j)),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(h)"
            WRITE( stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(sr)"
            WRITE( stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(ps)"
            WRITE( stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(nl)"
            WRITE( stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(xc)"
            WRITE( stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         ENDIF
      ENDIF
5555  FORMAT(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!
      RETURN
      END SUBROUTINE vofrho

!
!----------------------------------------------------------------------
      SUBROUTINE checkrho(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
!----------------------------------------------------------------------
!
!     check \int rho(r)dr and the negative part of rho
!
      IMPLICIT NONE
      INTEGER nnr, nspin
      REAL(8) rhor(nnr,nspin), rmin, rmax, rsum, rnegsum
!
      REAL(8) roe
      INTEGER ir, iss
!
      rsum   =0.0
      rnegsum=0.0
      rmin   =100.
      rmax   =0.0
      DO iss=1,nspin
         DO ir=1,nnr
            roe=rhor(ir,iss)
            rsum=rsum+roe
            IF (roe.LT.0.0) rnegsum=rnegsum+roe
            rmax=MAX(rmax,roe)
            rmin=MIN(rmin,roe)
         END DO
      END DO
      CALL reduce(1,rsum)
      CALL reduce(1,rnegsum)
      RETURN
    END SUBROUTINE checkrho
!______________________________________________________________________

