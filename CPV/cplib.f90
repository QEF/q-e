!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!
!-----------------------------------------------------------------------
      SUBROUTINE atomic_wfc( eigr, n_atomic_wfc, wfc )
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
!
      USE kinds,              ONLY: DP
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart, g, gx
      USE ions_base,          ONLY: nsp, na, nat
      USE cell_base,          ONLY: tpiba
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
!
      IMPLICIT NONE
      INTEGER,     INTENT(in) :: n_atomic_wfc
      COMPLEX(DP), INTENT(in) :: eigr( ngw, nat )
      COMPLEX(DP), INTENT(out):: wfc( ngw, n_atomic_wfc )
!
      INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
      REAL(DP), ALLOCATABLE ::  ylm(:,:), q(:), jl(:), vchi(:), chiq(:)
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, gx, g, ylm)
      ndm = MAXVAL(rgrid(1:nsp)%mesh)
      !
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
         DO nb = 1,upf(is)%nwfc
            l = upf(is)%lchi(nb)
            DO i=1,ngw
               CALL sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
               DO ir=1,rgrid(is)%mesh
                  vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
               ENDDO
               CALL simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
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

!-----------------------------------------------------------------------
FUNCTION n_atom_wfc_x( )
!----------------------------------------------------------------------------
  !
  ! ... Find max number of bands needed
  !
  USE ions_base,        ONLY : na, nsp
  USE kinds,            ONLY : DP
  USE uspp_param,       ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER  :: n_atom_wfc_x
  INTEGER  :: is, n
  !
  n_atom_wfc_x = 0
  !
  DO is = 1, nsp
     !
     DO n = 1, upf(is)%nwfc
        !
        IF ( upf(is)%oc(n) >= 0.D0 ) THEN
           !
           n_atom_wfc_x = n_atom_wfc_x + na(is) * ( 2*upf(is)%lchi(n) + 1 )
           !
        END IF
        !
     END DO
     !
  END DO
  !
  RETURN
END FUNCTION

!

!-----------------------------------------------------------------------
   FUNCTION cscnorm( bec, nkbx, cp, ngwx, i, n )
!-----------------------------------------------------------------------
!     requires in input the updated bec(i)
!
      USE ions_base,          ONLY: na
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE cvan,               ONLY: ish, nvb
      USE uspp_param,         ONLY: nh
      USE uspp,               ONLY: qq
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
      USE kinds,              ONLY: DP
!
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: i, n
      INTEGER, INTENT(IN) :: ngwx, nkbx
      REAL(DP)    :: bec( nkbx, n )
      COMPLEX(DP) :: cp( ngwx, n )
      !
      REAL(DP) :: cscnorm
      !
      INTEGER ig, is, iv, jv, ia, inl, jnl
      REAL(DP) rsum
      REAL(DP), ALLOCATABLE:: temp(:)
!
!
      ALLOCATE(temp(ngw))
      DO ig=1,ngw
         temp(ig)=DBLE(CONJG(cp(ig,i))*cp(ig,i))
      END DO
      rsum=2.d0*SUM(temp)
      IF (gstart == 2) rsum=rsum-temp(1)

      CALL mp_sum( rsum, intra_image_comm )

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
!
!-----------------------------------------------------------------------
      SUBROUTINE denlcc( nnr, nspin, vxcr, sfac, drhocg, dcc )
!-----------------------------------------------------------------------
!
! derivative of non linear core correction exchange energy wrt cell 
! parameters h 
! Output in dcc
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE reciprocal_vectors, ONLY: gstart, gx, ngs, g, ngm
      USE recvecs_indexes,    ONLY: np
      USE cell_base,          ONLY: omega, ainv, tpiba2
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
      USE uspp_param,         ONLY: upf
      USE grid_dimensions,    ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x
      USE cp_interfaces,      ONLY: fwfft

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
                 IF( upf(is)%nlcc ) srhoc = srhoc + sfac( ig, is ) * drhocg( ig, is )
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

      CALL mp_sum( dcc( 1:3, 1:3 ), intra_image_comm )

      RETURN
      END SUBROUTINE denlcc



!-----------------------------------------------------------------------
      SUBROUTINE dotcsc( eigr, cp, ngw, n )
!-----------------------------------------------------------------------
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: na, nsp, nat
      USE io_global,          ONLY: stdout
      USE reciprocal_vectors, ONLY: gstart
      USE cvan,               ONLY: ish, nvb
      USE uspp,               ONLY: nkb, qq
      USE uspp_param,         ONLY: nh
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ngw, n
      COMPLEX(DP) ::  eigr(ngw,nat), cp(ngw,n)
! local variables
      REAL(DP) rsum, csc(n) ! automatic array
      COMPLEX(DP) temp(ngw) ! automatic array
 
      REAL(DP), ALLOCATABLE::  becp(:,:)
      INTEGER i,kmax,nnn,k,ig,is,ia,iv,jv,inl,jnl
!
      ALLOCATE(becp(nkb,n))
!
!     < beta | phi > is real. only the i lowest:
!
      nnn = MIN( 12, n )

      DO i = nnn, 1, -1
         kmax = i
         CALL nlsm1(i,1,nvb,eigr,cp,becp)
!
         DO k=1,kmax
            DO ig=1,ngw
               temp(ig)=CONJG(cp(ig,k))*cp(ig,i)
            END DO
            csc(k)=2.d0*DBLE(SUM(temp))
            IF (gstart == 2) csc(k)=csc(k)-DBLE(temp(1))
         END DO

         CALL mp_sum( csc( 1:kmax ), intra_image_comm )

         DO k=1,kmax
            rsum=0.d0
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
         WRITE( stdout,'("dotcsc =",12f18.15)') (csc(k),k=1,i)
!
      END DO
      WRITE( stdout,*)
!
      DEALLOCATE(becp)
!
      RETURN
      END SUBROUTINE dotcsc


!-----------------------------------------------------------------------
   SUBROUTINE dotcsv( csv, eigr, c, v, ngw )
!-----------------------------------------------------------------------
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: na, nsp, nat
      USE io_global,          ONLY: stdout
      USE reciprocal_vectors, ONLY: gstart
      USE cvan,               ONLY: ish, nvb
      USE uspp,               ONLY: nkb, qq
      USE uspp_param,         ONLY: nh
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ngw
      COMPLEX(DP) ::  eigr(ngw,nat), c(ngw), v(ngw)
      REAL(DP), INTENT(OUT) ::  csv

      ! local variables
      COMPLEX(DP) temp(ngw) ! automatic array
 
      REAL(DP), ALLOCATABLE ::  bec(:), bev(:)
      INTEGER ig,is,ia,iv,jv,inl,jnl
!
      ALLOCATE(bec(nkb))
      ALLOCATE(bev(nkb))
!
!     < beta | c > is real. only the i lowest:
!
      CALL nlsm1(1,1,nvb,eigr,c,bec)
      CALL nlsm1(1,1,nvb,eigr,v,bev)
!
      DO ig=1,ngw
         temp(ig)=CONJG(c(ig))*v(ig)
      END DO
      csv = 2.0d0 * DBLE(SUM(temp))
      IF (gstart == 2) csv = csv - DBLE(temp(1))

      CALL mp_sum( csv, intra_image_comm )

      DO is=1,nvb
         DO iv=1,nh(is)
            DO jv=1,nh(is)
               DO ia=1,na(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  jnl=ish(is)+(jv-1)*na(is)+ia
                  csv = csv + qq(iv,jv,is)*bec(inl)*bev(jnl)
               END DO
            END DO
         END DO
      END DO
      !
      DEALLOCATE(bec)
      DEALLOCATE(bev)
      !
      RETURN
   END SUBROUTINE dotcsv



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
      USE kinds,                    ONLY: DP
      USE control_flags,            ONLY: iprint
      USE ions_base,                ONLY: na, nsp, nat
      USE cvan,                     ONLY: nvb
      USE uspp_param,               ONLY: nhm, nh
      USE grid_dimensions,          ONLY: nr1, nr2, nr3, &
                                          nr1x, nr2x, nr3x, nnr => nnrx
      USE electrons_base,           ONLY: nspin
      USE gvecb,                    ONLY: ngb, npb, nmb
      USE gvecp,                    ONLY: ng => ngm
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
                                          nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE cell_base,                ONLY: ainv
      USE qgb_mod,                  ONLY: qgb
      USE cdvan,                    ONLY: drhovan
      USE derho,                    ONLY: drhor, drhog
      USE dqgb_mod,                 ONLY: dqgb
      USE recvecs_indexes,          ONLY: nm, np
      USE cp_interfaces,            ONLY: fwfft, invfft
      USE fft_base,                 ONLY: dfftb

      IMPLICIT NONE
! input
      INTEGER, INTENT(in) ::  irb(3,nat)
      REAL(DP), INTENT(in)::  rhor(nnr,nspin)
      REAL(DP) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      COMPLEX(DP), INTENT(in)::  eigrb(ngb,nat), rhog(ng,nspin)
! local
      INTEGER i, j, isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss,   &
     &     isa, ia, ir
      REAL(DP) sum, dsum
      COMPLEX(DP) fp, fm, ci
      COMPLEX(DP), ALLOCATABLE :: v(:)
      COMPLEX(DP), ALLOCATABLE:: dqgbt(:,:)
      COMPLEX(DP), ALLOCATABLE :: qv(:)
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

      IF ( nvb == 0 ) RETURN

      ALLOCATE( v( nnr ) )
      ALLOCATE( qv( nnrb ) )
      ALLOCATE( dqgbt( ngb, 2 ) )

      ci =( 0.0d0, 1.0d0 )

      IF( nspin == 1 ) THEN
         !  
         !  nspin=1 : two fft at a time, one per atom, if possible
         ! 
         DO i=1,3
            DO j=1,3

               v(:) = (0.d0, 0.d0)

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
                                 sum =2.d0*sum
                                 dsum=2.d0*dsum
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
                     CALL box2grid( irb(1,isa), 1, qv, v )
                     IF (nfft.EQ.2) CALL box2grid(irb(1,isa+1),2,qv,v)

  15                 isa = isa + nfft
!
                  END DO
               END DO
!
               DO ir=1,nnr
                  drhor(ir,iss,i,j) = drhor(ir,iss,i,j) + DBLE(v(ir))
               END DO
!
               CALL fwfft( 'Dense', v, nr1, nr2, nr3, nr1x, nr2x, nr3x )
!
               DO ig=1,ng
                  drhog(ig,iss,i,j) = drhog(ig,iss,i,j) + v(np(ig))
               END DO
!
            ENDDO
         ENDDO
!
      ELSE
         !
         !     nspin=2: two fft at a time, one for spin up and one for spin down
         ! 
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
                                 sum =2.d0*sum
                                 dsum=2.d0*dsum
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
                     !
  25                 isa = isa + 1
                     !
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
     &                 0.5d0*CMPLX( DBLE(fp),AIMAG(fm))
                  drhog(ig,isdw,i,j) = drhog(ig,isdw,i,j) +             &
     &                 0.5d0*CMPLX(AIMAG(fp),-DBLE(fm))
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
      USE mp_global,          ONLY: intra_image_comm
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
         sk(i)=0.0d0
         DO ig=gstart,ngw
            sk(i)=sk(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*ggp(ig)
         END DO
      END DO

      CALL mp_sum( sk(1:n), intra_image_comm )

      enkin=0.0d0
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
      SUBROUTINE gausin(eigr,cm)
!-----------------------------------------------------------------------
!
! initialize wavefunctions with gaussians - edit to fit your system
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: na, nsp, nat
      USE electrons_base,     ONLY: n => nbsp
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gx, g
!
      IMPLICIT NONE
!
      COMPLEX(DP) eigr(ngw,nat), cm(ngw,n)
      REAL(DP)    sigma, auxf
      INTEGER nband, is, ia, ig, isa
!
      sigma=12.0d0
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
      USE mp_global,      ONLY: intra_image_comm
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

      CALL mp_sum( csc( 1:kmax ), intra_image_comm )

      !
      !     calculate bec(i)=<cp(i)|beta>
      !
      DO inl=1,nhsavb
         DO ig=1,ngw
            temp(ig)=cp(1,ig,i)* DBLE(betae(ig,inl))+             &
     &               cp(2,ig,i)*AIMAG(betae(ig,inl))
         END DO
         bec(inl,i)=2.d0*SUM(temp)
         IF (gstart == 2) bec(inl,i)= bec(inl,i)-temp(1)
      END DO

      CALL mp_sum( bec( 1:nhsavb, i ), intra_image_comm )
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
      DO k=1,kmax
         IF (ispin(i).EQ.ispin(k)) THEN
            rsum=0.d0
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
          DO inl=1,nkbx
            bec(inl,i)=bec(inl,i)-csc(k)*bec(inl,k)
         END DO
      END DO

      DEALLOCATE( temp )
!
      RETURN
      END SUBROUTINE gracsc


!-------------------------------------------------------------------------
      SUBROUTINE smooth_csv( c, v, ngwx, csv, n )
!-----------------------------------------------------------------------

      USE gvecw,              ONLY: ngw
      USE kinds,              ONLY: DP
      USE reciprocal_vectors, ONLY: gstart
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ngwx, n
      REAL(DP)    :: c( 2, ngwx )
      REAL(DP)    :: v( 2, ngwx, n )
      REAL(DP)    :: csv( n )
      INTEGER     :: k, ig
      REAL(DP), ALLOCATABLE :: temp(:) 

      !
      !     calculate csv(k)=<c|v(k)>
      !
      ALLOCATE( temp( ngw ) )

      DO k = 1, n
         DO ig = 1, ngw
            temp(ig) = v(1,ig,k) * c(1,ig) + v(2,ig,k) * c(2,ig)
         END DO
         csv(k) = 2.0d0 * SUM(temp)
         IF (gstart == 2) csv(k) = csv(k) - temp(1)
      END DO

      DEALLOCATE( temp )
!
      RETURN
      END SUBROUTINE smooth_csv


!-------------------------------------------------------------------------
      SUBROUTINE grabec( becc, nkbx, betae, c, ngwx )
!-----------------------------------------------------------------------
      !
      !     on output: bec(i) is recalculated
      !
      USE uspp,           ONLY : nkb, nhsavb=>nkbus
      USE gvecw,          ONLY: ngw
      USE kinds,          ONLY: DP
      USE reciprocal_vectors, ONLY: gstart
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx
      COMPLEX(DP) :: betae( ngwx, nkb )
      REAL(DP)    :: becc( nkbx ), c( 2, ngwx )
      INTEGER     :: ig, inl
      REAL(DP), ALLOCATABLE :: temp(:) 
      !
      ALLOCATE( temp( ngw ) )
      !
      !     calculate becc=<c|beta>
      !
      DO inl=1,nhsavb
         DO ig=1,ngw
            temp(ig)=c(1,ig)* DBLE(betae(ig,inl))+             &
     &               c(2,ig)*AIMAG(betae(ig,inl))
         END DO
         becc(inl)=2.d0*SUM(temp)
         IF (gstart == 2) becc(inl)= becc(inl)-temp(1)
      END DO

      DEALLOCATE( temp )

      RETURN
      END SUBROUTINE grabec


!-------------------------------------------------------------------------
      SUBROUTINE bec_csv( becc, becv, nkbx, csv, n )
!-----------------------------------------------------------------------
!     requires in input the updated becc and becv(k)
!     on output: csv is updated
!
      USE ions_base,      ONLY: na
      USE cvan,           ONLY :nvb, ish
      USE uspp,           ONLY : nkb, nhsavb=>nkbus, qq
      USE uspp_param,     ONLY:  nh
      USE kinds,          ONLY: DP
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, n
      REAL(DP)    :: becc( nkbx )
      REAL(DP)    :: becv( nkbx, n )
      REAL(DP)    :: csv( n )
      INTEGER     :: k, is, iv, jv, ia, inl, jnl
      REAL(DP)    :: rsum

!     calculate csv(k) = csv(k) + <c| SUM_nm |beta(n)><beta(m)|v(k)>,  k<i
!
      DO k=1,n
            rsum=0.d0
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN 
                        DO ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           rsum = rsum + qq(iv,jv,is)*becc(inl)*becv(jnl,k)
                        END DO
                     ENDIF
                  END DO
               END DO
            END DO
            csv(k)=csv(k)+rsum
      END DO
!
      RETURN
      END SUBROUTINE bec_csv



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
         CALL DSCAL( 2*ngw, 1.0d0/anorm, cp(1,i), 1 )
         !
         !         these are the final bec's
         !
         CALL DSCAL( nkbx, 1.0d0/anorm, bec(1,i), 1 )
      END DO
!
      DEALLOCATE( csc )

      CALL stop_clock( 'gram' )
!
      RETURN
      END SUBROUTINE gram
!
!-----------------------------------------------------------------------
      SUBROUTINE initbox ( tau0, taub, irb, ainv, a1, a2, a3 )
!-----------------------------------------------------------------------
!
!     sets the indexes irb and positions taub for the small boxes 
!     around atoms
!
      USE kinds,                    ONLY: DP
      USE ions_base,                ONLY: nsp, na, nat
      USE grid_dimensions,          ONLY: nr1, nr2, nr3
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
      USE control_flags,            ONLY: iprsta
      USE io_global,                ONLY: stdout
      USE mp_global,                ONLY: nproc_image, me_image
      USE fft_base,                 ONLY: dfftb, dfftp, fft_dlay_descriptor
      USE fft_types,                ONLY: fft_box_set
      USE cvan,                     ONLY: nvb

      IMPLICIT NONE
! input
      REAL(DP), INTENT(in)  :: tau0(3,nat)
! output
      INTEGER,  INTENT(out) :: irb(3,nat)
      REAL(DP), INTENT(out) :: taub(3,nat)
! input
      REAL(DP), INTENT(in)  :: ainv(3,3)
      REAL(DP), INTENT(in)  :: a1(3)
      REAL(DP), INTENT(in)  :: a2(3)
      REAL(DP), INTENT(in)  :: a3(3)
! local
      REAL(DP) :: x(3), xmod
      INTEGER  :: nr(3), nrb(3), xint, is, ia, i, isa
!
      IF ( nr1b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 1)
      IF ( nr2b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 2)
      IF ( nr3b < 1) CALL errore &
         ('initbox', 'incorrect value for box grid dimensions', 3)

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
                        nat, irb, me_image+1, nproc_image, dfftp%npp, dfftp%ipp )

      IF( iprsta > 2 ) THEN
           isa = 1
           DO is=1,nsp
              WRITE( stdout, '( /, 2x, "species= ", i2 )' ) is
              DO ia=1,na(is)
                 WRITE( stdout,2000) ia, (irb(i,isa),i=1,3)
2000             FORMAT(2x, 'atom= ', i3, ' irb1= ', i3, ' irb2= ', i3, ' irb3= ', i3)
                 isa = isa + 1
               END DO
            END DO
      ENDIF

#ifdef __PARA
      ! 
      ! for processor that do not call fft on the box
      ! artificially start the clock
      ! 
      CALL start_clock( 'fftb' )
      CALL stop_clock( 'fftb' )
      !
#endif
!
      RETURN
   END SUBROUTINE initbox
!
!-------------------------------------------------------------------------
   SUBROUTINE newd(vr,irb,eigrb,rhovan,fion)
!-----------------------------------------------------------------------
!
!     this routine calculates array deeq:
!         deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!     and the corresponding term in forces
!         fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!     where
!         rho_lm = \sum_j f_j <psi_j|beta_l><beta_m|psi_j>
!
      USE kinds,                    ONLY: dp
      USE uspp_param,               ONLY: nh, nhm
      USE uspp,                     ONLY: deeq
      USE cvan,                     ONLY: nvb
      USE ions_base,                ONLY: nat, nsp, na
      USE parameters,               ONLY: nsx
      USE constants,                ONLY: pi, fpi
      USE grid_dimensions,          ONLY: nr3, nnr => nnrx
      USE gvecb,                    ONLY: ngb, npb, nmb, gxb
      USE small_box,                ONLY: omegab, tpibab
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
                                          nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE qgb_mod,                  ONLY: qgb
      USE electrons_base,           ONLY: nspin
      USE control_flags,            ONLY: iprint, thdyn, tfor, tprnfor
      USE mp,                       ONLY: mp_sum
      USE mp_global,                ONLY: intra_image_comm
      USE cp_interfaces,            ONLY: invfft
      USE fft_base,                 ONLY: dfftb
!
      IMPLICIT NONE
! input
      INTEGER irb(3,nat)
      REAL(DP) rhovan(nhm*(nhm+1)/2,nat,nspin)
      COMPLEX(DP) eigrb(ngb,nat)
      REAL(DP)  vr(nnr,nspin)
! output
      REAL(DP)  fion(3,nat)
! local
      INTEGER isup,isdw,iss, iv,ijv,jv, ik, nfft, isa, ia, is, ig
      REAL(DP)  fvan(3,nat,nsx), fac, fac1, fac2, boxdotgrid
      COMPLEX(DP) ci, facg1, facg2
      COMPLEX(DP), ALLOCATABLE :: qv(:)
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

      CALL mp_sum( deeq, intra_image_comm )

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

      CALL mp_sum( fvan, intra_image_comm )

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
      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE ions_base,         ONLY: na, nsp, nat
      USE uspp,              ONLY: nhsa=>nkb, qq
      USE uspp_param,        ONLY: nhm, nh
      USE cvan,              ONLY: ish, nvb
      USE electrons_base,    ONLY: nbspx, nbsp, nudx, nspin, iupdwn, nupdwn
      USE constants,         ONLY: pi, fpi
      USE cp_main_variables, ONLY: nlam, nlax, descla, la_proc
      USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ 
      USE mp,                ONLY: mp_sum
      USE mp_global,         ONLY: intra_image_comm
!
      IMPLICIT NONE
      REAL(DP) bec(nhsa,nbsp), becdr(nhsa,nbsp,3), lambda(nlam,nlam,nspin)
      REAL(DP) fion(3,nat)
!
      INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc
      REAL(DP), ALLOCATABLE :: temp(:,:), tmpbec(:,:),tmpdr(:,:) 
      REAL(DP), ALLOCATABLE :: fion_tmp(:,:)
      !
      CALL start_clock( 'nlfl' )
      !
      ALLOCATE( fion_tmp( 3, nat ) )
      !
      fion_tmp = 0.0d0
      !

      ALLOCATE( temp( nlax, nlax ), tmpbec( nhm, nlax ), tmpdr( nlax, nhm ) )

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
                  IF( la_proc ) THEN
                     ! tmpbec distributed by columns
                     ic = descla( ilac_ , iss )
                     nc = descla( nlac_ , iss )
                     DO iv=1,nh(is)
                        DO jv=1,nh(is)
                           inl=ish(is)+(jv-1)*na(is)+ia
                           IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                              DO i=1,nc
                                 tmpbec(iv,i)=tmpbec(iv,i) + qq(iv,jv,is)*bec(inl,i+istart-1+ic-1)
                              END DO
                           ENDIF
                        END DO
                     END DO
                     ! tmpdr distributed by rows
                     ir = descla( ilar_ , iss )
                     nr = descla( nlar_ , iss )
                     DO iv=1,nh(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        DO i=1,nr
                           tmpdr(i,iv)=becdr(inl,i+istart-1+ir-1,k)
                        END DO
                     END DO
                  END IF
!
                  IF(nh(is).GT.0)THEN
                     !
                     !  CALL DGEMM &
                     !  ( 'N', 'N', nss, nss, nh(is), 1.0d0, tmpdr, nudx, tmpbec, nhm, 0.0d0, temp, nudx )
                     !
                     IF( la_proc ) THEN
                        ir = descla( ilar_ , iss )
                        ic = descla( ilac_ , iss )
                        nr = descla( nlar_ , iss )
                        nc = descla( nlac_ , iss )
                        CALL DGEMM( 'N', 'N', nr, nc, nh(is), 1.0d0, tmpdr, nlax, tmpbec, nhm, 0.0d0, temp, nlax )
                        DO j = 1, nc
                           DO i = 1, nr
                              fion_tmp(k,isa) = fion_tmp(k,isa) + 2D0 * temp( i, j ) * lambda( i, j, iss )
                           END DO
                        END DO
                     END IF
!
                  ENDIF

               END DO
!
            END DO
         END DO
      END DO
      !
      DEALLOCATE( temp, tmpbec, tmpdr )
      !
      CALL mp_sum( fion_tmp, intra_image_comm )
      !
      fion = fion + fion_tmp
      !
      DEALLOCATE( fion_tmp )
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
      USE kinds,  ONLY: DP

      IMPLICIT NONE
! input
      REAL(DP) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
! output
      REAL(DP) rout(3)
! local
      REAL(DP) x,y,z
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
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nsp, na, nat
      USE gvecw,      ONLY : ngw
      USE cvan,       ONLY : ish
      USE uspp,       ONLY : nkb, beta, nhtol
      USE uspp_param, ONLY : nh
!
      IMPLICIT NONE
      COMPLEX(DP) :: eigr( ngw, nat )
      COMPLEX(DP) :: betae( ngw, nkb )
!
      INTEGER     :: is, iv, ia, inl, ig, isa
      COMPLEX(DP) :: ci
!
      CALL start_clock( 'prefor' )
      isa = 0
      DO is=1,nsp
         DO iv=1,nh(is)
            ci=(0.0d0,-1.0d0)**nhtol(iv,is)
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
   SUBROUTINE projwfc( c, nx, eigr, betae, n, ei  )
!-----------------------------------------------------------------------
      !
      ! Projection on atomic wavefunctions
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: autoev
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE ions_base,          ONLY: nsp, na, nat
      USE uspp,               ONLY: nhsa => nkb
      USE uspp_param,         ONLY: upf
!
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: nx, n
      COMPLEX(DP), INTENT(IN) :: c( ngw, nx ), eigr(ngw,nat), betae(ngw,nhsa)
      REAL(DP),    INTENT(IN) :: ei( nx )
!
      COMPLEX(DP), ALLOCATABLE :: wfc(:,:), swfc(:,:), becwfc(:,:)
      REAL(DP),    ALLOCATABLE :: overlap(:,:), e(:), z(:,:)
      REAL(DP),    ALLOCATABLE :: proj(:,:), temp(:)
      REAL(DP)                 :: somma

      INTEGER :: n_atomic_wfc
      INTEGER :: is, ia, nb, l, m, k, i
      !
      ! calculate number of atomic states
      !
      n_atomic_wfc = 0
      !
      DO is=1,nsp
         DO nb = 1,upf(is)%nwfc
            l = upf(is)%lchi(nb)
            n_atomic_wfc = n_atomic_wfc + (2*l+1)*na(is)
         END DO
      END DO
      IF ( n_atomic_wfc .EQ. 0 ) RETURN
      !
      ALLOCATE( wfc( ngw, n_atomic_wfc ) )
      !
      ! calculate wfc = atomic states
      !
      CALL atomic_wfc( eigr, n_atomic_wfc, wfc )
      !
      ! calculate bec = <beta|wfc>
      !
      ALLOCATE( becwfc( nhsa, n_atomic_wfc ) )
      !
      CALL nlsm1( n_atomic_wfc, 1, nsp, eigr, wfc, becwfc )
      !
      ! calculate swfc = S|wfc>
      !
      ALLOCATE( swfc( ngw, n_atomic_wfc ) )
      !
      CALL s_wfc( n_atomic_wfc, becwfc, betae, wfc, swfc )
      !
      ! calculate overlap(i,j) = <wfc_i|S|wfc_j> 
      !
      ALLOCATE( overlap( n_atomic_wfc, n_atomic_wfc ) )

      CALL DGEMM &
           ( 'T', 'N', n_atomic_wfc, n_atomic_wfc, 2*ngw, 1.0d0, wfc, 2*ngw, &
             swfc, 2*ngw, 0.0d0, overlap, n_atomic_wfc )

!      CALL MXMA( wfc, 2*ngw, 1, swfc, 1, 2*ngw, overlap, 1,                     &
!     &          n_atomic_wfc, n_atomic_wfc, 2*ngw, n_atomic_wfc )

      CALL mp_sum( overlap, intra_image_comm )

      overlap = overlap * 2.d0
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
      !
      CALL rdiag( n_atomic_wfc, overlap, n_atomic_wfc, e, z )
      !
      overlap=0.d0
      !
      DO l=1,n_atomic_wfc
         DO m=1,n_atomic_wfc
            DO k=1,n_atomic_wfc
               overlap(l,m)=overlap(l,m)+z(m,k)*z(l,k)/SQRT(e(k))
            END DO
         END DO
      END DO
      !
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

      CALL mp_sum( proj, intra_image_comm )

      i=0
      WRITE( stdout,  90 ) 
      WRITE( stdout, 100 )
      DO is=1,nsp
         DO nb = 1,upf(is)%nwfc
            l=upf(is)%lchi(nb)
            DO m = -l,l
               DO ia=1,na(is)
                  i=i+1
               END DO
               WRITE( stdout, 110 ) i-na(is)+1, i, na(is), is, nb, l, m
            END DO
         END DO
      END DO

      WRITE( stdout,*)
      DO m=1,n
         somma=0.d0
         DO l=1,n_atomic_wfc
            somma=somma+proj(m,l)**2
         END DO
         WRITE( stdout, 120 ) m, somma, ei(m)*autoev
         WRITE( stdout, 130 ) (ABS(proj(m,l)),l=1,n_atomic_wfc)
      END DO

  90  FORMAT( 3X,'Projection on atomic states')
 100  FORMAT( 3X,'atomic state    atom   specie  wfc  l  m')
 110  FORMAT( 3X, I4, ' - ', I4, 4X, I4, 6X, I3,   I5, I4,I3)
 120  FORMAT( 3X,'state # ',i4,'    sum c^2 = ',f7.4, ' eV = ', F7.2 )
 130  FORMAT( 3X, 10f7.4)

!
      DEALLOCATE(proj)
      DEALLOCATE(wfc)
      RETURN
      END SUBROUTINE projwfc

!
!-----------------------------------------------------------------------
      SUBROUTINE rdiag ( n, h, ldh, e, v )
!-----------------------------------------------------------------------
!
!   calculates all the eigenvalues and eigenvectors of a complex
!   hermitean matrix H . On output, the matrix H is destroyed
!
      USE kinds,            ONLY: DP
      USE dspev_module,     ONLY: dspev_drv
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(in)   :: n, ldh
      REAL(DP), INTENT(inout):: h(ldh,n)
      REAL(DP), INTENT(out)  :: e(n)
      REAL(DP), INTENT(out)  :: v(ldh,n)
!
      INTEGER :: i, j, k
      REAL(DP), ALLOCATABLE :: ap( : )
!
      ALLOCATE( ap( n * ( n + 1 ) / 2 ) )

      K = 0
      DO J = 1, n
         DO I = J, n
            K = K + 1
            ap( k ) = h( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, ap, e, v, ldh )

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
      USE ions_base, ONLY: nat, na, nsp
      USE io_global, ONLY: stdout
      USE mp_global, ONLY: intra_image_comm
      USE mp, ONLY: mp_sum
      USE cvan, ONLY: nvb
      USE uspp_param, ONLY: nh, nhm
      USE uspp, ONLY: deeq
      USE grid_dimensions, ONLY: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      USE electrons_base, ONLY: nspin
      USE gvecb
      USE gvecp, ONLY: ng => ngm
      USE cell_base, ONLY: omega, ainv
      USE small_box, ONLY: omegab
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      USE control_flags, ONLY: iprint, iprsta, tpre
      USE qgb_mod
      USE recvecs_indexes, ONLY: np, nm
      USE cp_interfaces, ONLY: fwfft, invfft
      USE fft_base, ONLY: dfftb
!
      IMPLICIT NONE
!
      REAL(DP) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      INTEGER, INTENT(in) :: irb(3,nat)
      COMPLEX(DP), INTENT(in):: eigrb(ngb,nat)
      REAL(DP), INTENT(inout):: rhor(nnr,nspin)
      COMPLEX(DP),  INTENT(inout):: rhog(ng,nspin)
!
      INTEGER     :: isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss, isa, ia, ir, i, j
      REAL(DP)    :: sumrho
      COMPLEX(DP) :: ci, fp, fm, ca
      COMPLEX(DP), ALLOCATABLE :: qgbt(:,:)
      COMPLEX(DP), ALLOCATABLE :: v(:)
      COMPLEX(DP), ALLOCATABLE :: qv(:)

      !  Quick return if this sub is not needed
      !
      IF ( nvb == 0 ) RETURN

      CALL start_clock( 'rhov' )
      ci=(0.d0,1.d0)
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
                        IF(iv.NE.jv) sumrho=2.d0*sumrho
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

            CALL mp_sum( ca, intra_image_comm )

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
         DO ig = 1, ng
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
                        IF(iv.NE.jv) sumrho=2.d0*sumrho
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
            CALL mp_sum( ca, intra_image_comm )
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
            rhog(ig,isup)=rhog(ig,isup) + 0.5d0*CMPLX(DBLE(fp),AIMAG(fm))
            rhog(ig,isdw)=rhog(ig,isdw) + 0.5d0*CMPLX(AIMAG(fp),-DBLE(fm))
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
      SUBROUTINE s_wfc(n_atomic_wfc1,becwfc,betae,wfc,swfc) !@@@@ Changed n_atomic_wfc to n_atomic_wfc1
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
      USE kinds, ONLY: DP
      USE ions_base, ONLY: na
      USE cvan, ONLY: nvb, ish
      USE uspp, ONLY: nhsa => nkb, nhsavb=>nkbus, qq
      USE uspp_param, ONLY: nh
      USE gvecw, ONLY: ngw
      USE constants, ONLY: pi, fpi
      IMPLICIT NONE
! input
      INTEGER, INTENT(in)         :: n_atomic_wfc1
      COMPLEX(DP), INTENT(in) :: betae(ngw,nhsa),                   &
     &                               wfc(ngw,n_atomic_wfc1)
      REAL(DP), INTENT(in)    :: becwfc(nhsa,n_atomic_wfc1)
! output
      COMPLEX(DP), INTENT(out):: swfc(ngw,n_atomic_wfc1)
! local
      INTEGER is, iv, jv, ia, inl, jnl, i
      REAL(DP) qtemp(nhsavb,n_atomic_wfc1)
!
      swfc = wfc
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
                        DO i=1,n_atomic_wfc1
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*becwfc(jnl,i)
                        END DO
                     END DO
                  ENDIF
               END DO
            END DO
         END DO
!
!         CALL MXMA (betae,1,2*ngw,qtemp,1,nhsavb,swfc,1,                &
!     &              2*ngw,2*ngw,nhsavb,n_atomic_wfc1)
!
         CALL DGEMM &
              ('N','N',2*ngw,n_atomic_wfc1,nhsavb,1.0d0,betae,2*ngw,&
               qtemp,nhsavb,1.0d0,swfc,2*ngw)
!
      END IF
!
!      swfc=swfc+wfc
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
      USE mp_global, ONLY: intra_image_comm
      USE mp, ONLY: mp_sum
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
      fup = 0.0d0
      DO i=iupdwn(1),nupdwn(1)
         fup = fup + f(i)
      END DO
      nup = NINT(fup)
      ndw = nel(1)+nel(2) - nup
!
! paranoid checks
!
      frac= ABS(fup-nup).GT.1.0d-6
      fup = 0.0d0
      DO i=1,nup
         fup = fup + f(i)
      END DO
      frac=frac.OR.ABS(fup-nup).GT.1.0d-6
      fdw = 0.0d0
      DO j=iupdwn(2),iupdwn(2)-1+ndw
         fdw = fdw + f(j)
      END DO
      frac=frac.OR.ABS(fdw-ndw).GT.1.0d-6
!
      spin0 = ABS(fup-fdw)/2.d0 * ( ABS(fup-fdw)/2.d0 + 1.d0 ) + fdw
!
!     Becke's formula for spin polarization
!
      spin1 = 0.0d0
      DO ir=1,nnr
         spin1 = spin1 - MIN(rhor(ir,1),rhor(ir,2))
      END DO
      CALL mp_sum( spin1, intra_image_comm )
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
      CALL mp_sum( overlap, intra_image_comm )
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
      USE kinds,              ONLY: dp
      USE control_flags,      ONLY: iprint, iprsta, thdyn, tpre, tfor, &
                                    tprnfor, iesr
      USE io_global,          ONLY: stdout
      USE ions_base,          ONLY: nsp, na, nat, rcmax
      USE gvecs
      USE gvecp,              ONLY: ng => ngm
      USE cell_base,          ONLY: omega, r_to_s
      USE cell_base,          ONLY: a1, a2, a3, tpiba2, h, ainv
      USE reciprocal_vectors, ONLY: gstart, g, gx
      USE recvecs_indexes,    ONLY: np, nm
      USE grid_dimensions,    ONLY: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnr => nnrx
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
                                        nr1sx, nr2sx, nr3sx, nnrsx
      USE electrons_base,   ONLY: nspin
      USE constants,        ONLY: pi, fpi, au_gpa
      USE energies,         ONLY: etot, eself, enl, ekin, epseu, esr, eht, exc 
      USE local_pseudo,     ONLY: vps, dvps, rhops
      USE core,             ONLY: nlcc_any
      USE gvecb
      USE dener,            ONLY: detot, dekin, dps, dh, dsr, dxc, denl, &
                                  detot6, dekin6, dps6, dh6, dsr6, dxc6, denl6
      USE derho
      USE mp,               ONLY: mp_sum
      USE mp_global,        ONLY: intra_image_comm
      USE funct,            ONLY: dft_is_meta
      USE pres_ai_mod,      ONLY: abivol, abisur, v_vol, P_ext, volclu,  &
                                  Surf_t, surfclu
      USE cp_interfaces,    ONLY: fwfft, invfft, self_vofhar
      USE sic_module,       ONLY: self_interaction, sic_epsilon, sic_alpha
      USE energies,         ONLY: self_exc, self_ehte
      USE cp_interfaces,    ONLY: pseudo_stress, compute_gagb, stress_hartree, &
                                  add_drhoph, stress_local, force_loc
!@@@@@
      USE ldaU,             ONLY: e_hubbard
!@@@@@
      IMPLICIT NONE
!
      LOGICAL :: tlast, tfirst
      INTEGER :: nfi
      REAL(DP)  rhor(nnr,nspin), rhos(nnrsx,nspin), fion(3,nat)
      REAL(DP)  rhoc(nnr), tau0(3,nat)
      COMPLEX(DP) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),     &
     &                ei3(-nr3:nr3,nat), eigrb(ngb,nat),        &
     &                rhog(ng,nspin), sfac(ngs,nsp)
      !
      INTEGER irb(3,nat)
      !
      INTEGER iss, isup, isdw, ig, ir, i, j, k, ij, is, ia
      REAL(DP) vave, ebac, wz, eh, ehpre
      COMPLEX(DP)  fp, fm, ci, drhop
      COMPLEX(DP), ALLOCATABLE :: rhotmp(:), vtemp(:)
      ! COMPLEX(DP), ALLOCATABLE :: drhotmp(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: drhot(:,:)
      COMPLEX(DP), ALLOCATABLE :: v(:), vs(:)
      REAL(DP), ALLOCATABLE    :: gagb(:,:)
      !
      REAL(DP), ALLOCATABLE :: fion1( :, : )
      REAL(DP), ALLOCATABLE :: stmp( :, : )
      !
      COMPLEX(DP), ALLOCATABLE :: self_vloc(:)
      COMPLEX(DP)              :: self_rhoeg
      REAL(DP)                 :: self_ehtet, fpibg
      LOGICAL                  :: ttsic
      REAL(DP)                 :: detmp( 3, 3 ), desr( 6 ), deps( 6 )
      REAL(DP)                 :: detmp2( 3, 3 )
      REAL(DP)                 :: ht( 3, 3 )
      REAL(DP)                 :: deht( 6 )
      COMPLEX(DP)              :: screen_coul( 1 )
!
      INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)

      ! ...  dalbe(:) = delta( alpha(:), beta(:) )
      REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)



      CALL start_clock( 'vofrho' )

      ci = ( 0.0d0, 1.0d0 )
      !
      !     wz = factor for g.neq.0 because of c*(g)=c(-g)
      !
      wz = 2.0d0
      !
      ht = TRANSPOSE( h )
      !
      ALLOCATE( v( nnr ) )
      ALLOCATE( vs( nnrsx ) )
      ALLOCATE( vtemp( ng ) )
      ALLOCATE( rhotmp( ng ) )
      !
      IF ( tpre ) THEN
         ALLOCATE( drhot( ng, 6 ) )
         ALLOCATE( gagb( 6, ng ) )
         CALL compute_gagb( gagb, gx, ng, tpiba2 )
      END IF
!
!     ab-initio pressure and surface tension contributions to the potential
!
      if (abivol.or.abisur) call vol_clu(rhor,rhog,sfac,nfi)
      !
      ttsic = ( ABS( self_interaction ) /= 0 )
      !
      IF( ttsic ) ALLOCATE( self_vloc( ng ) )
      !
      !     first routine in which fion is calculated: annihilation
      !
      fion  = 0.d0
      !
      !     forces on ions, ionic term in real space
      !
      IF( tprnfor .OR. tfor .OR. tfirst .OR. tpre ) THEN
         !
         ALLOCATE( stmp( 3, nat ) )
         !
         CALL r_to_s( tau0, stmp, na, nsp, ainv )
         !
         CALL vofesr( iesr, esr, dsr6, fion, stmp, tpre, h )
         !
         call mp_sum( fion, intra_image_comm )
         !
         DEALLOCATE( stmp )
         !
      END IF
!
      rhotmp( 1:ng ) = rhog( 1:ng, 1 )
      !

      IF( tpre ) THEN
         DO ij = 1, 6
            i = alpha( ij )
            j = beta( ij )
            drhot( :, ij ) = 0.0d0
            DO k = 1, 3
               drhot( :, ij ) = drhot( :, ij ) +  drhog( :, 1, i, k ) * ht( k, j )
            END DO
         END DO
      END IF
      !
      IF( nspin == 2 ) THEN
         rhotmp( 1:ng ) = rhotmp( 1:ng ) + rhog( 1:ng, 2 )
         IF(tpre)THEN
            DO ij = 1, 6
               i = alpha( ij )
               j = beta( ij )
               DO k = 1, 3
                  drhot( :, ij ) = drhot( :, ij ) +  drhog( :, 2, i, k ) * ht( k, j )
               END DO
            END DO
         ENDIF
      END IF
      !
      !     calculation local potential energy
      !
      vtemp=(0.d0,0.d0)
      DO is=1,nsp
         DO ig=1,ngs
            vtemp(ig)=vtemp(ig)+CONJG(rhotmp(ig))*sfac(ig,is)*vps(ig,is)
         END DO
      END DO

      epseu = wz * DBLE(SUM(vtemp))
      IF (gstart == 2) epseu=epseu-vtemp(1)
      CALL mp_sum( epseu, intra_image_comm )

      epseu = epseu * omega

!
      IF( tpre ) THEN
         !
         CALL stress_local( dps6, epseu, gagb, sfac, rhotmp, drhot, omega )
         !
      END IF
      !
      !     
      !     calculation hartree energy
      !    
      !
      self_ehtet = 0.d0  
      !
      IF( ttsic ) self_vloc = 0.d0 

      DO is=1,nsp
         DO ig=1,ngs
            rhotmp(ig)=rhotmp(ig)+sfac(ig,is)*rhops(ig,is)
         END DO
      END DO
      !
      IF (gstart == 2) vtemp(1)=0.0d0
      DO ig = gstart, ng
         vtemp(ig) = CONJG( rhotmp( ig ) ) * rhotmp( ig ) / g( ig )
      END DO
!
      eh = DBLE( SUM( vtemp ) ) * wz * 0.5d0 * fpi / tpiba2
!
      IF ( ttsic ) THEN
         !
         CALL self_vofhar( .false., self_ehte, self_vloc, rhog, omega, h )
         !
         eh = eh - self_ehte / omega

         CALL mp_sum( self_ehte, intra_image_comm )
         !
      END IF
      !
      CALL mp_sum( eh, intra_image_comm )
      !
      IF(tpre) THEN
         !
         CALL add_drhoph( drhot, sfac, gagb )
         !
         CALL stress_hartree(dh6, eh*omega, sfac, rhotmp, drhot, gagb, omega )
         !
      END IF
      !
      IF(tpre) THEN
         DEALLOCATE( drhot )
      END IF
      !    
      !     forces on ions, ionic term in reciprocal space
      !     
      ALLOCATE( fion1( 3, nat ) )
      !
      fion1 = 0.d0
      !
      IF( tprnfor .OR. tfor .OR. tpre) THEN
          vtemp( 1:ng ) = rhog( 1:ng, 1 )
          IF( nspin == 2 ) THEN
             vtemp( 1:ng ) = vtemp(1:ng) + rhog( 1:ng, 2 )
          END IF
          CALL force_loc( .false., vtemp, fion1, rhops, vps, ei1, ei2, ei3, sfac, omega, screen_coul )
      END IF
      !
      !     calculation hartree + local pseudo potential
      !
      !
      IF (gstart == 2) vtemp(1)=(0.d0,0.d0)
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
      IF ( nlcc_any ) CALL add_cc( rhoc, rhog, rhor )
!
      CALL exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_exc )


!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
      IF( nspin == 1 ) THEN
         iss = 1
         if (abivol.or.abisur) then
            do ir=1,nnr
               v(ir)=CMPLX( rhor( ir, iss ) + v_vol( ir ), 0.d0 )
            end do           
         else
            do ir=1,nnr
               v(ir)=CMPLX( rhor( ir, iss ), 0.d0 )
            end do
         end if
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
         if (abivol.or.abisur) then
            do ir=1,nnr
               v(ir)=CMPLX(rhor(ir,isup)+v_vol(ir),rhor(ir,isdw)+v_vol(ir))
            end do
         else
            do ir=1,nnr
               v(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw))
            end do
         end if
         CALL fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         DO ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            IF( ttsic ) THEN
             rhog(ig,isup)=vtemp(ig)-self_vloc(ig) +0.5d0*CMPLX( DBLE(fp),AIMAG(fm))
             rhog(ig,isdw)=vtemp(ig)+self_vloc(ig) +0.5d0*CMPLX(AIMAG(fp),-DBLE(fm))
            ELSE
             rhog(ig,isup)=vtemp(ig)+0.5d0*CMPLX( DBLE(fp),AIMAG(fm))
             rhog(ig,isdw)=vtemp(ig)+0.5d0*CMPLX(AIMAG(fp),-DBLE(fm))
            ENDIF
         END DO
      ENDIF

!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
      IF( tprnfor .OR. tfor ) THEN

         IF ( nlcc_any ) CALL force_cc( irb, eigrb, rhor, fion1 )

         CALL mp_sum( fion1, intra_image_comm )
         !
         !    add g-space ionic and core correction contributions to fion
         !
         fion = fion + fion1

      END IF

      DEALLOCATE( fion1 )
!
      IF( ttsic ) DEALLOCATE( self_vloc )
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
         vave=(SUM(rhor(:,isup))+SUM(rhor(:,isdw))) / 2.0d0 / DBLE( nr1 * nr2 * nr3 )
      ENDIF

      CALL mp_sum( vave, intra_image_comm )

      !
      !     fourier transform of total potential to r-space (smooth grid)
      !
      vs (:) = (0.d0, 0.d0)
      !
      IF(nspin.EQ.1)THEN
         !
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
         !
      ELSE
         !
         isup=1
         isdw=2
         DO ig=1,ngs
            vs(nps(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            vs(nms(ig))=CONJG(rhog(ig,isup)) +ci*CONJG(rhog(ig,isdw))
         END DO 
         !
         CALL invfft('Smooth',vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
         !
         DO ir=1,nnrsx
            rhos(ir,isup)= DBLE(vs(ir))
            rhos(ir,isdw)=AIMAG(vs(ir))
         END DO
         !
      ENDIF

      IF( dft_is_meta() ) CALL vofrho_meta( v, vs )  !METAGGA

      ebac = 0.0d0
      !
      eht = eh * omega + esr - eself
      !
      !     etot is the total energy ; ekin, enl were calculated in rhoofr
      !
      etot = ekin + eht + epseu + enl + exc + ebac +e_hubbard
      !
      if (abivol) etot = etot + P_ext*volclu
      if (abisur) etot = etot + Surf_t*surfclu
      !
      IF( tpre ) THEN
         !
         detot6 = dekin6 + dh6 + dps6 + dsr6
         !
         call mp_sum( detot6, intra_image_comm )
         !
         DO k = 1, 6
            detmp( alpha(k), beta(k) ) = detot6(k)
            detmp( beta(k), alpha(k) ) = detmp( alpha(k), beta(k) )
         END DO
         !
         detot = MATMUL( detmp(:,:), TRANSPOSE( ainv(:,:) ) )
         !
         detot = detot + denl + dxc
         !
      END IF
      !
      !
      CALL stop_clock( 'vofrho' )
      !
      !
      IF ( tpre ) THEN
         !
         DEALLOCATE( gagb )
         !
         IF( ( iprsta >= 2 ) .AND. ( MOD( nfi - 1, iprint) == 0 ) ) THEN  
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "From vofrho:"
            WRITE( stdout,*) "cell parameters h"
            WRITE( stdout,5555) (a1(i),a2(i),a3(i),i=1,3)
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(tot)"
            WRITE( stdout,5555) ((detot(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( detot, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*)
            WRITE( stdout,*) "derivative of e(kin)"
            WRITE( stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dekin, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(h)"
            WRITE( stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dh, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
             !
            WRITE( stdout,*) "derivative of e(sr)"
            WRITE( stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dsr, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(ps)"
            WRITE( stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dps, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(nl)"
            WRITE( stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( denl, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
            !
            WRITE( stdout,*) "derivative of e(xc)"
            WRITE( stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "kbar"
            detmp = -1.0d0 * MATMUL( dxc, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
         ENDIF
      ENDIF

      DEALLOCATE( rhotmp )
      DEALLOCATE( vtemp )
      DEALLOCATE( v )
      DEALLOCATE( vs )

      RETURN

5555  FORMAT(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!

      END SUBROUTINE vofrho

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
      subroutine ldaU_init
!-----------------------------------------------------------------------
!
      USE constants,        ONLY: autoev
      use ldaU,             ONLY: n_atomic_wfc, atomwfc,lda_plus_u, Hubbard_U
      use ldaU,             ONLY: Hubbard_lmax, Hubbard_l, ns, vupsi
      use input_parameters, ONLY: atom_label, lda_plus_u_ => lda_plus_u
      use input_parameters, ONLY: Hubbard_U_ => Hubbard_U
      use ions_base,        only: na, nsp, nat
      use gvecw,            only: ngw
      use electrons_base,   only: nspin, nx => nbspx
      USE uspp_param,       ONLY: upf
      !
      implicit none
      integer is, nb, l
      integer, external :: set_Hubbard_l

! allocate vupsi

      lda_plus_u = lda_plus_u_

      allocate(vupsi(ngw,nx))

      vupsi=(0.0d0,0.0d0)
      ! allocate(vpsi_con(ngw,nx)) ! step_constraint 
      n_atomic_wfc=0

      do is=1,nsp
         !
         Hubbard_U( is ) = Hubbard_U_( is )/autoev
         !
         do nb = 1,upf(is)%nwfc
            l = upf(is)%lchi(nb)
            n_atomic_wfc = n_atomic_wfc + (2*l+1)*na(is)
         end do
         !
      end do
!
      allocate(atomwfc(ngw,n_atomic_wfc))

      if (lda_plus_u) then
         Hubbard_lmax = -1
         do is=1,nsp
            if (Hubbard_U(is).ne.0.d0) then 
!                Hubbard_l(is)=2
               Hubbard_l(is) = set_Hubbard_l( atom_label(is) )
                Hubbard_lmax = max(Hubbard_lmax,Hubbard_l(is))
               write (6,*) ' HUBBARD L FOR TYPE ',atom_label(is),' IS ',&
     &                       Hubbard_l(is)
            end if
         end do
         write (6,*) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
         if (Hubbard_lmax.eq.-1) call errore                            &
     &        ('setup','lda_plus_u calculation but Hubbard_l not set',1)
      end if
      l = 2 * Hubbard_lmax + 1
      allocate(ns(nat,nspin,l,l))
      return
      end subroutine ldaU_init
!
!-----------------------------------------------------------------------
integer function set_Hubbard_l(psd) result (hubbard_l)
!-----------------------------------------------------------------------
!
implicit none
character*3 :: psd
!
! TRANSITION METALS
!
if (psd.eq.'V'  .or. psd.eq.'Cr' .or. psd .eq.'Mn' .or. psd.eq.'Fe' .or. &
    psd.eq.'Co' .or. psd.eq.'Ni' .or. psd .eq.'Cu'.or. psd .eq.'Fe1'.or. &
    psd .eq.'Fe2' ) then
    hubbard_l = 2
!
! RARE EARTHS
!
elseif (psd .eq.'Ce') then
   hubbard_l =  3
!
! OTHER ELEMENTS
!
elseif (psd .eq.'H') then
   hubbard_l =  0
elseif (psd .eq.'O') then
   hubbard_l = 1
else
   hubbard_l = -1
   call errore ('set_Hubbard_l','pseudopotential not yet inserted', 1)
endif
return
end function set_Hubbard_l
!
!-----------------------------------------------------------------------
      subroutine new_ns(c,eigr,betae,hpsi,hpsi_con,forceh)
!-----------------------------------------------------------------------
!
! This routine computes the on site occupation numbers of the Hubbard ions.
! It also calculates the contribution of the Hubbard Hamiltonian to the
! electronic potential and to the forces acting on ions.
!
      use control_flags,      ONLY: tfor, tprnfor
      use kinds,              ONLY: DP        
      use parameters,         ONLY: nsx
      use ions_base,          only: na, nat, nsp
      use gvecw,              only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      USE uspp,               ONLY: nhsa=>nkb
      USE uspp_param,         ONLY: upf
      use electrons_base,     only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE ldaU,               ONLY: lda_plus_u, Hubbard_U, Hubbard_l
      USE ldaU,               ONLY: n_atomic_wfc, ns, e_hubbard
!
      implicit none
#ifdef __PARA
      include 'mpif.h'
#endif
      integer, parameter :: ldmx = 7
      complex(DP), intent(in) :: c(ngw,nx), eigr(ngw,nat),      &
     &                               betae(ngw,nhsa)
      complex(DP), intent(out) :: hpsi(ngw,nx), hpsi_con(1,1)
      real(DP) forceh(3,nat)
!
      complex(DP), allocatable:: wfc(:,:), swfc(:,:),dphi(:,:,:),   &
     &                               spsi(:,:)
      real(DP), allocatable   :: becwfc(:,:), bp(:,:),              &
     &                               dbp(:,:,:), wdb(:,:,:)
      real(DP), allocatable   :: dns(:,:,:,:)
      real(DP), allocatable   :: e(:), z(:,:),                      &
     &                               proj(:,:), temp(:)
      real(DP), allocatable   :: ftemp1(:), ftemp2(:)
      real(DP)                :: lambda(ldmx), somma, SSUM, ntot,   &
     &                               nsum, nsuma, x_value, g_value,     &
     &                               step_value
      real(DP) :: f1 (ldmx, ldmx), vet (ldmx, ldmx)
      integer is, ia, iat, nb, isp, l, m, m1, m2, k, i, counter, err, ig
      integer iv, jv, inl, jnl,alpha,alpha_a,alpha_s,ipol
      integer, allocatable ::  offset (:,:)
      complex(DP) :: tempsi
!
!
      allocate(wfc(ngw,n_atomic_wfc))
      allocate(ftemp1(ldmx))
      allocate(ftemp2(ldmx))
!
! calculate wfc = atomic states
!
!!!      call ewfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
      allocate(becwfc(nhsa,n_atomic_wfc))
!!!      call nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)
!
      allocate(swfc(ngw,n_atomic_wfc))
!!!      call s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate proj = <c|S|wfc>
!
      allocate(proj(n,n_atomic_wfc))
      CALL projwfc_hub( c, nx, eigr, betae, n, n_atomic_wfc,            &
     & wfc, becwfc, swfc, proj ) !@@
!
      allocate(offset(nsp,nat))
      counter = 0
      do is = 1, nsp
         do ia = 1, na(is)
            do i = 1, upf(is)%nwfc
               l = upf(is)%lchi(i)
               if (l.eq.Hubbard_l(is)) offset (is,ia) = counter
               counter = counter + 2 * l + 1
            end do
         end do
      end do
      if (counter.ne.n_atomic_wfc)                                      &
     &                 call errore ('new_ns','nstart<>counter',1)
      ns(:,:,:,:) = 0.d0
      iat = 0
      do is = 1,nsp
         do ia = 1,na(is)
            iat = iat + 1
            if (Hubbard_U(is).ne.0.d0) then 
               k = offset(is,ia)
               do m1 = 1, 2*Hubbard_l(is) + 1
                  do m2 = m1, 2*Hubbard_l(is) + 1
                     do i = 1,n
!                      write(6,*) i,ispin(i),f(i)
                      ns(iat,ispin(i),m1,m2) = ns(iat,ispin(i),m1,m2) + &
     &                               f(i) * proj(i,k+m2) * proj(i,k+m1)
                     end do
!                     ns(iat,:,m2,m1) = ns(iat,:,m1,m2)
                     ns(iat,1,m2,m1) = ns(iat,1,m1,m2)
                     ns(iat,2,m2,m1) = ns(iat,2,m1,m2)
                  end do
               end do
            end if
         end do
      end do
      if (nspin.eq.1) ns = 0.5d0 * ns
! Contributions to total energy
      e_hubbard = 0.d0
      iat = 0
      do is = 1,nsp
         do ia = 1,na(is)
            iat=iat + 1
            if (Hubbard_U(is).ne.0.d0) then
                k = offset(is,ia)
                do isp = 1,nspin
                   do m1 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard + 0.5d0 * Hubbard_U(is) *    &
     &                           ns(iat,isp,m1,m1)
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        e_hubbard = e_hubbard - 0.5d0 * Hubbard_U(is) * &
     &                              ns(iat,isp,m1,m2) * ns(iat,isp,m2,m1)
                     end do
                   end do
                end do
             end if
         end do
       end do
       if (nspin.eq.1) e_hubbard = 2.d0*e_hubbard
!       if (nspin.eq.1) e_lambda = 2.d0*e_lambda
!
!      Calculate the potential and forces on wavefunctions due to U
!
      hpsi(:,:)=(0.d0,0.d0)
      iat=0
      do is = 1, nsp
         do ia=1, na(is)
            iat = iat + 1
            if (Hubbard_U(is).ne.0.d0) then
               do i=1, n
                  do m1 = 1, 2 * Hubbard_l(is) + 1
                     tempsi = proj (i,offset(is,ia)+m1)
                     do m2 = 1, 2 * Hubbard_l(is) + 1
                        tempsi = tempsi - 2.d0 * ns(iat,ispin(i),m1,m2)*&
     &                                proj (i,offset(is,ia)+m2)
                     enddo
                     tempsi = tempsi * Hubbard_U(is)/2.d0*f(i)
                     call ZAXPY (ngw,tempsi,swfc(1,offset(is,ia)+m1),1, &
     &                           hpsi(1,i),1)
                  enddo
               enddo
            endif
         enddo
      enddo
!
!      Calculate the potential and energy due to constraint
!
      hpsi_con(:,:)=0.d0
!
! Calculate the contribution to forces on ions due to U and constraint
!
      forceh=0.d0
      if ((tfor).or.(tprnfor)) then
        allocate (bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3))
        allocate(dns(nat,nspin,ldmx,ldmx))
        allocate (spsi(ngw,n))
!
        call nlsm1 (n,1,nsp,eigr,c,bp)
        call s_wfc(n,bp,betae,c,spsi)
        call nlsm2(ngw,nhsa,n,eigr,c,dbp,.true.)
        call nlsm2(ngw,nhsa,n_atomic_wfc,eigr,wfc,wdb,.true.)
!
        alpha=0
        do alpha_s = 1, nsp
         do alpha_a = 1, na(alpha_s)
            alpha=alpha+1
            do ipol = 1,3
               call dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,      &
     &                    offset,c,wfc,eigr,betae,proj,ipol,dns)
               iat=0
               do is = 1, nsp
                  do ia=1, na(is)
                     iat = iat + 1
                     if (Hubbard_U(is).ne.0.d0) then
                        do isp = 1,nspin
                           do m2 = 1,2*Hubbard_l(is) + 1
                              forceh(ipol,alpha) = forceh(ipol,alpha) -            &
     &                        Hubbard_U(is) * 0.5d0 * dns(iat,isp,m2,m2)
                              do m1 = 1,2*Hubbard_l(is) + 1
                                 forceh(ipol,alpha) = forceh(ipol,alpha) +         &
     &                           Hubbard_U(is)*ns(iat,isp,m2,m1)*       &
     &                           dns(iat,isp,m1,m2)
                              end do
                           end do
                        end do
                     end if
! Occupation constraint add here
                  end do
               end do
            end do
         end do
        end do
        if (nspin.eq.1) then
           forceh = 2.d0 * forceh
        end if
!
        deallocate ( wfc, becwfc, spsi, proj, offset, swfc, dns, bp, dbp, wdb)
      end if
      return
      end subroutine new_ns
!
!
!
!-----------------------------------------------------------------------
      subroutine write_ns
!-----------------------------------------------------------------------
!
! This routine computes the occupation numbers on atomic orbitals.
! It also write the occupation number in the output file.
!
      USE kinds,            only: DP
      USE constants,        ONLY: autoev
      use electrons_base,   only: nspin
      use electrons_base,   only: n => nbsp 
      use ions_base,        only: na, nat, nsp
      use gvecw,            only: ngw
      USE ldaU,             ONLY: lda_plus_u, Hubbard_U, Hubbard_l
      USE ldaU,             ONLY: n_atomic_wfc, ns, e_hubbard
      USE ldaU,             ONLY: Hubbard_lmax
      use dspev_module,     only : dspev_drv

      implicit none

  integer :: is, isp, ia, m1, m2, ldim, iat, err, k
! cpunter on atoms type
! counter on spin component
! counter on atoms
! counter on wavefn
! counters on d components
  integer, parameter :: ldmx = 7
  real(DP), allocatable   :: ftemp1(:), ftemp2(:)
  real(DP) :: f1 (ldmx * ldmx), vet (ldmx, ldmx)
  real(DP) :: lambda (ldmx), nsum, nsuma
  write (*,*) 'enter write_ns'

  if ( 2 * Hubbard_lmax + 1 .gt. ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

!  if (step_con) then
!     do isp=1,nspin
!        write (6,'(6(a,i2,a,i2,a,f8.4,6x))') &
!        ('A_con(',is,',',isp,') =', A_con(is,isp),is=1,nsp)
!     enddo
!     write (6,'(6(a,i2,a,f8.4,6x))') &
!           ('sigma_con(',is,') =', sigma_con(is), is=1,nsp)
!     write (6,'(6(a,i2,a,f8.4,6x))') &
!        ('alpha_con(',is,') =', alpha_con(is), is=1,nsp)
!  endif
  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('U(',is,') =', Hubbard_U(is) * autoev, is=1,nsp)
!  write (6,'(6(a,i2,a,f8.4,6x))') &
!        ('alpha(',is,') =', Hubbard_alpha(is) * autoev, is=1,nsp)
      nsum = 0.d0
      allocate(ftemp1(ldmx))
      allocate(ftemp2(ldmx))
      iat = 0
      write(6,*) 'nsp',nsp
      do is = 1,nsp
         do ia = 1, na(is)
            nsuma = 0.d0
            iat = iat + 1
!        if (iat.eq.1) then
            if (Hubbard_U(is).ne.0.d0) then
               do isp = 1, nspin
                   do m1 = 1, 2 * Hubbard_l(is) + 1
                      nsuma = nsuma + ns (iat, isp, m1, m1)
                   end do
               end do
               if (nspin.eq.1) nsuma = 2.d0 * nsuma
               write(6,'(a,x,i2,2x,a,f11.7)') 'atom', iat,              &
     &                                      ' Tr[ns(na)]= ',nsuma
               nsum = nsum + nsuma
!
               do isp = 1, nspin

                  k = 0
                  do m1 = 1, 2 * Hubbard_l(is) + 1
                     do m2 = m1, 2 * Hubbard_l(is) + 1
                        k = k + 1
                        f1 ( k ) = ns (iat, isp, m2, m1)
                     enddo
                  enddo

                  CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1, f1, lambda, vet, ldmx  )

                  write(6,'(a,x,i2,2x,a,x,i2)') 'atom', iat, 'spin', isp
                  write(6,'(a,7f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,&
     &                                2 * Hubbard_l(is) + 1)
                  write(6,*) 'eigenvectors'
                  do m2 = 1, 2*Hubbard_l(is)+1
                     write(6,'(i2,2x,7(f10.7,x))') m2,(real(vet(m1,m2)),&
     &                            m1=1,2 * Hubbard_l(is) + 1)
                  end do
                  write(6,*) 'occupations'
                  do m1 = 1, 2*Hubbard_l(is)+1
                     write (6,'(7(f6.3,x))') (ns(iat,isp,m1,m2),m2=1,    &
     &                     2*Hubbard_l(is)+1)
                  end do
               end do
            end if
!        end if
         end do
      end do
      deallocate ( ftemp1, ftemp2)
      return
      end subroutine write_ns
!-----------------------------------------------------------------------
      subroutine genatwfc(n_atomic_wfc,atwfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space, in the same order as used in new_ns
!
      use ions_base,          only: na, nsp
      use gvecw,              only: ngw
      use reciprocal_vectors, only: g, gx, ng0 => gstart
      use cell_base,          only: omega, tpiba
      use constants,          only: fpi
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
      USE kinds,              ONLY: DP
!
      implicit none
      integer, intent(in) :: n_atomic_wfc
      complex(DP), intent(out):: atwfc(ngw,n_atomic_wfc)
!
      integer natwfc, is, ia, ir, nb, l, m, lm, i, lmax_wfc, ig
      real(DP), allocatable::  ylm(:,:), q(:), jl(:), vchi(:),        &
     &     chiq(:), gxn(:,:)
!
!
      allocate(q(ngw))
      allocate(gxn(3,ngw))
      allocate(chiq(ngw))
!
      do ig=1,ngw
         q(ig) = sqrt(g(ig))*tpiba
      end do
      if (ng0.eq.2) gxn(1,:)=0.0d0
      do ig=ng0,ngw
         gxn(:,ig) = gx(:,ig)/sqrt(g(ig)) !ik<=>ig
      end do
!
      natwfc=0
!@@@@@
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX (lmax_wfc, MAXVAL ( upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, gx, g, ylm)
!@@@@@

      do is = 1, nsp
         ALLOCATE  ( jl(rgrid(is)%mesh), vchi(rgrid(is)%mesh) )
         do ia=1,na(is)
!
!   radial fourier transform of the chi functions
!   NOTA BENE: chi is r times the radial part of the atomic wavefunction
!              bess requires l+1, not l, on input
!
            do nb = 1,upf(is)%nwfc
               l = upf(is)%lchi(nb)
               do i=1,ngw
                  call sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
                  do ir=1,rgrid(is)%mesh
                     vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
                  enddo
                  call simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
               enddo
!
!   multiply by angular part and structure factor
!   NOTA BENE: the factor i^l MUST be present!!!
!
               do m = 1,2*l+1
                  lm = l**2 + m
!                  call ylmr2b(lm,ngw,ngw,gxn,ylm)
                  natwfc = natwfc + 1
                  atwfc(:,natwfc) = (0.d0,1.d0)**l * ylm(:,lm)*chiq(:)
               enddo
            enddo
         end do
         DEALLOCATE  ( vchi, jl )
      end do
!
      do i = 1,natwfc
        call DSCAL(2*ngw,fpi/sqrt(omega),atwfc(1,i),1)
      end do
!
      if (natwfc.ne.n_atomic_wfc)                                       &
     &     call errore('atomic_wfc','unexpected error',natwfc)
!
      deallocate(ylm)
      deallocate(chiq)
      deallocate(gxn)
      deallocate(q)
!
      return
      end subroutine genatwfc
!
!-------------------------------------------------------------------------
      subroutine dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,         &
     &                  offset,c,wfc,                                   &
     &                  eigr,betae,                                     &
     &                  proj,ipol,dns)
!-----------------------------------------------------------------------
!
! This routine computes the derivative of the ns with respect to the ionic
! displacement tau(alpha,ipol) used to obtain the Hubbard contribution to the
! atomic forces.
!
      use ions_base, only: na, nat, nsp
      use gvecw, only: ngw
      use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE uspp,           ONLY: nhsa=>nkb
      USE ldaU,           ONLY: Hubbard_U, Hubbard_l
      USE ldaU,           ONLY: n_atomic_wfc, ns
      USE kinds,          ONLY: DP
!
      implicit none
      integer, parameter :: ldmx = 7
      integer ibnd,is,i,ia,counter, m1,m2, l, iat, alpha, ldim
! input
      integer,      intent(in) :: offset(nsp,nat)
      integer,      intent(in) :: alpha_a,alpha_s,ipol
      real(DP),     intent(in) :: wfc(ngw,n_atomic_wfc),  c(2,ngw,nx),  &
     &                            eigr(2,ngw,nat),betae(2,ngw,nhsa),    &
     &                            becwfc(nhsa,n_atomic_wfc),            &
     &                            bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3)
      real(DP),     intent(in) :: proj(n,n_atomic_wfc)
      complex (DP), intent(in) :: spsi(ngw,n)
! output
      real (DP),   intent(out) :: dns(nat,nspin,ldmx,ldmx)
!
!     dns !derivative of ns(:,:,:,:) w.r.t. tau
!
      real (DP),   allocatable :: dproj(:,:)
!
!     dproj(n,n_atomic_wfc) ! derivative of proj(:,:) w.r.t. tau 
!
      allocate (dproj(n,n_atomic_wfc) )
!
      dns(:,:,:,:) = 0.d0
!
          call dprojdtau(c,wfc,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,     &
     &                   alpha_s,ipol,offset(alpha_s,alpha_a),dproj)
!
! compute the derivative of occupation numbers (the quantities dn(m1,m2))
! of the atomic orbitals. They are real quantities as well as n(m1,m2)
!
      iat=0
      do is=1,nsp
         do ia = 1,na(is)
            iat=iat+1
            if (Hubbard_U(is).ne.0.d0) then
               ldim = 2*Hubbard_l(is) + 1
               do m1 = 1, ldim
                  do m2 = m1, ldim
                     do ibnd = 1,n
                        dns(iat,ispin(ibnd),m1,m2) =                    &
     &                  dns(iat,ispin(ibnd),m1,m2) +                    &
     &                   f(ibnd)*REAL(  proj(ibnd,offset(is,ia)+m1) *   &
     &                   (dproj(ibnd,offset(is,ia)+m2))  +              &
     &                         dproj(ibnd,offset(is,ia)+m1)  *          &
     &                         (proj(ibnd,offset(is,ia)+m2)) )
                     end do
                     dns(iat,:,m2,m1) = dns(iat,:,m1,m2)
                  end do
               end do
            end if
         end do
      end do
!
      deallocate (dproj)
      return
      end subroutine dndtau
!
!
!-----------------------------------------------------------------------
      subroutine dprojdtau(c,wfc,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,    &
     &                     alpha_s,ipol,offset,dproj)
!-----------------------------------------------------------------------
!
! This routine computes the first derivative of the projection
! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
!
      use ions_base, only: na, nat
      use gvecw, only: ngw
      use reciprocal_vectors, only: g, gx, ng0 => gstart
      use electrons_base, only: n => nbsp, nx => nbspx
!      use gvec
!      use constants
      USE uspp,           ONLY: nhsa=>nkb, qq
      use cvan,           ONLY: ish
      USE ldaU,           ONLY: Hubbard_U, Hubbard_l
      USE ldaU,           ONLY: n_atomic_wfc
      use cell_base,      ONLY: tpiba
      USE uspp_param,     only: nh !@@@@
      use mp_global,      only: intra_image_comm
      use mp,             only: mp_sum
      USE kinds,          ONLY: DP
!
       implicit none
       integer, parameter :: ldmx = 7
       integer alpha_a, alpha_s,ipol, offset
! input: the displaced atom
! input: the component of displacement
! input: the offset of the wfcs of the atom "alpha_a,alpha_s"
       complex (DP), intent(in) :: spsi(ngw,n),                     &
     &                  c(ngw,nx), eigr(ngw,nat)
! input: the atomic wfc
! input: S|evc>
       real(DP), intent(in) ::becwfc(nhsa,n_atomic_wfc),            &
     &                            wfc(2,ngw,n_atomic_wfc),              &
     &            bp(nhsa,n), dbp(nhsa,n,3), wdb(nhsa,n_atomic_wfc,3)
       real(DP), intent(out) :: dproj(n,n_atomic_wfc)
! output: the derivative of the projection
!
      integer i,ig,m1,ibnd,iwf,ia,is,iv,jv,ldim,alpha,l,m,k,inl
!
      real (DP)  a1, a2
      real(kind=8), allocatable :: gk(:)
!
      complex (DP), allocatable :: dwfc(:,:)
      real (DP), allocatable :: betapsi(:,:),                       &
     &                              dbetapsi(:,:),                      &
     &                              wfcbeta(:,:),wfcdbeta(:,:),temp(:)
!      dwfc(ngw,ldmx),             ! the derivative of the atomic d wfc
!      betapsi(nh,n),              ! <beta|evc>
!      dbetapsi(nh,n),             ! <dbeta|evc>
!      wfcbeta(n_atomic_wfc,nh),   ! <wfc|beta>
!      wfcdbeta(n_atomic_wfc,nh),  ! <wfc|dbeta>
      ldim = 2 * Hubbard_l(alpha_s) + 1
      allocate ( dwfc(ngw,ldmx),betapsi(nh(alpha_s),n))
      allocate ( dbetapsi(nh(alpha_s),n),                               &
     &           wfcbeta(n_atomic_wfc,nh(alpha_s)))
      allocate (wfcdbeta(n_atomic_wfc,nh(alpha_s)) )
      dproj(:,:)=0.d0
!
! At first the derivative of the atomic wfc is computed
!
!
      allocate(gk(ngw))
      allocate(temp(ngw))
!
      if (Hubbard_U(alpha_s).ne.0.d0) then
!
         do ig=1,ngw
            gk(ig)=gx(ipol,ig)*tpiba 
!
            do m1=1,ldim
                  dwfc(ig,m1) = cmplx (gk(ig)*wfc(2,ig,offset+m1),      &
     &                  -1*gk(ig)*wfc(1,ig,offset+m1) )
            end do
         end do
!
         do ibnd=1,n
            do m1=1,ldim
               temp(:)=real(conjg(dwfc(:,m1))*spsi(:,ibnd))
               dproj(ibnd,offset+m1)=2.d0*SUM(temp) 
               if (ng0.eq.2) dproj(ibnd,offset+m1)=dproj(ibnd,offset+m1)-temp(1)
            end do
         end do
         call mp_sum( dproj, intra_image_comm )
      end if
      do iv=1,nh(alpha_s)
         inl=ish(alpha_s)+(iv-1)*na(alpha_s)+alpha_a
         do i=1,n
            betapsi(iv,i)=bp(inl,i)
            dbetapsi(iv,i)=dbp(inl,i,ipol)
         end do
         do m=1,n_atomic_wfc
!                 do m1=1,2**Hubbard_l(is) + 1
            wfcbeta(m,iv)=becwfc(inl,m)
            wfcdbeta(m,iv)=wdb(inl,m,ipol)
         end do
      end do
      do ibnd=1,n
         do iv=1,nh(alpha_s)
            do jv=1,nh(alpha_s)
               do m=1,n_atomic_wfc
!                       do m1=1,2**Hubbard_l(is) + 1
                  dproj(ibnd,m) =                                       &
     &                        dproj(ibnd,m) + qq(iv,jv,alpha_s) *       &
     &                         ( wfcdbeta(m,iv)*betapsi(jv,ibnd) +      &
     &                           wfcbeta(m,iv)*dbetapsi(jv,ibnd) )
               end do
            end do
         end do
      end do
      deallocate(temp, gk)
      deallocate (betapsi)
      deallocate (dwfc)
      deallocate (dbetapsi)
      deallocate (wfcbeta)
      deallocate (wfcdbeta)
      return
      end subroutine dprojdtau
!
!
!-----------------------------------------------------------------------
      subroutine stepfn(A,sigma,x_value,g_value,step_value)
!-----------------------------------------------------------------------
!     This subroutine calculates the value of the gaussian and step
!     functions with a given x_value. A and sigma are given in the
!     input file. ... to be used in occupation_constraint...
!
      USE constants, ONLY : pi
      implicit none
      real(kind=8) A, sigma, x_value, g_value, step_value
      real(kind=8) x
      integer i
      step_value=0.0d0
      g_value=0.0d0
!
      do i=1,100000
         x=x_value + (i-100000)/100000.0d0*(x_value + 5.d0*sigma)
!
! Integrate from 5 sigma before the x_value
!
         g_value=A*dexp(-x*x/(2*sigma*sigma))/(sigma*dsqrt(2*pi))
!         write(6,*) 'step', step_value,'g',g_value
!         if (g_value.le.0.0) g_value=0.0
         if ((x_value+5*sigma).ge.0.0d0) then
         step_value=step_value+g_value/100000.0d0*(x_value+5.d0*sigma)
         end if
      end do
      return
      end subroutine stepfn
!
!-----------------------------------------------------------------------
      SUBROUTINE projwfc_hub( c, nx, eigr, betae, n, n_atomic_wfc,  &
     & wfc, becwfc, swfc, proj )
!-----------------------------------------------------------------------
      !
      ! Projection on atomic wavefunctions
      ! Atomic wavefunctions are not orthogonized
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: autoev
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE ions_base,          ONLY: nsp, na, nat
      USE uspp,               ONLY: nhsa => nkb
!
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: nx, n
      COMPLEX(DP), INTENT(IN) :: c( ngw, nx ), eigr(ngw,nat), betae(ngw,nhsa)
!
      COMPLEX(DP), INTENT(OUT):: wfc(ngw,n_atomic_wfc),    &
     & swfc( ngw, n_atomic_wfc )
      real(DP), intent(out):: becwfc(nhsa,n_atomic_wfc) !DEBUG
      REAL(DP),    ALLOCATABLE :: overlap(:,:), e(:), z(:,:)
      REAL(DP),    ALLOCATABLE :: temp(:)
      REAL(DP)                 :: somma, proj(n,n_atomic_wfc)
      INTEGER :: n_atomic_wfc
      INTEGER :: is, ia, nb, l, m, k, i
      !
      ! calculate number of atomic states
      !
      !
      IF ( n_atomic_wfc .EQ. 0 ) RETURN
      !
      !
      ! calculate wfc = atomic states
      !
      CALL atomic_wfc_northo( eigr, n_atomic_wfc, wfc )
      !
      ! calculate bec = <beta|wfc>
      !
      CALL nlsm1( n_atomic_wfc, 1, nsp, eigr, wfc, becwfc )
      !
      ! calculate swfc = S|wfc>
      !
      CALL s_wfc( n_atomic_wfc, becwfc, betae, wfc, swfc )
      !
      ! calculate proj = <c|S|wfc>
      !
      ALLOCATE(temp(ngw))
      DO m=1,n
         DO l=1,n_atomic_wfc
            temp(:)=DBLE(CONJG(c(:,m))*swfc(:,l)) !@@@@
            proj(m,l)=2.d0*SUM(temp)
            IF (gstart == 2) proj(m,l)=proj(m,l)-temp(1)
         END DO
      END DO
      DEALLOCATE(temp)
      CALL mp_sum( proj, intra_image_comm )
!
      RETURN
      END SUBROUTINE projwfc_hub
!
!-----------------------------------------------------------------------
      SUBROUTINE atomic_wfc_northo( eigr, n_atomic_wfc, wfc )
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
! Atomic wavefunctions not orthogonalized
!
      USE kinds,              ONLY: DP
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart, g, gx
      USE ions_base,          ONLY: nsp, na, nat
      USE cell_base,          ONLY: tpiba, omega !@@@@
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
!@@@@@
      USE constants,          ONLY: fpi
!@@@@@
!
      IMPLICIT NONE
      INTEGER,     INTENT(in) :: n_atomic_wfc
      COMPLEX(DP), INTENT(in) :: eigr( ngw, nat )
      COMPLEX(DP), INTENT(out):: wfc( ngw, n_atomic_wfc )
!
      INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
      REAL(DP), ALLOCATABLE ::  ylm(:,:), q(:), jl(:), vchi(:), chiq(:)
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, gx, g, ylm)
      ndm = MAXVAL(rgrid(1:nsp)%mesh)
      !
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
         DO ia = 1 + isa, na(is) + isa
            DO nb = 1,upf(is)%nwfc
               l = upf(is)%lchi(nb)
               DO i=1,ngw
                  CALL sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
                  DO ir=1,rgrid(is)%mesh
                     vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
                  ENDDO
                  CALL simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
               ENDDO
               !
               !   multiply by angular part and structure factor
               !   NOTA BENE: the factor i^l MUST be present!!!
               !
               DO m = 1,2*l+1
                  lm = l**2 + m
                  !DO ia = 1 + isa, na(is) + isa
                  natwfc = natwfc + 1
                  wfc(:,natwfc) = (0.d0,1.d0)**l * eigr(:,ia)* ylm(:,lm)*chiq(:)
                  !ENDDO
               ENDDO
            ENDDO
         ENDDO
         isa = isa + na(is)
      ENDDO
!
      IF (natwfc.NE.n_atomic_wfc)                                       &
     &     CALL errore('atomic_wfc','unexpected error',natwfc)
!
!@@@@@
      do i = 1,n_atomic_wfc
        call DSCAL(2*ngw,fpi/sqrt(omega),wfc(1,i),1)
      end do
!@@@@@
      DEALLOCATE(q, chiq, vchi, jl, ylm)
!
      RETURN
      END SUBROUTINE atomic_wfc_northo
