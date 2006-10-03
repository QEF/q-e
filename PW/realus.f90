!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE realus
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  ! ... module originally written by Antonio Suriano and Stefano de Gironcoli
  ! ... modified by Carlo Sbraccia
  !
  INTEGER,  ALLOCATABLE :: box(:,:), maxbox(:)
  REAL(DP), ALLOCATABLE :: qsave(:,:)
  REAL(DP), ALLOCATABLE :: boxradius(:), boxdistance(:,:), xyz(:,:,:)     
  REAL(DP), ALLOCATABLE :: spher(:,:,:)
  LOGICAL               :: tqr
  !
  CONTAINS
    !
    SUBROUTINE deallocatenewdreal()
      !
      IF ( ALLOCATED( box ) )         DEALLOCATE( box )
      IF ( ALLOCATED( boxdistance ) ) DEALLOCATE( boxdistance )
      IF ( ALLOCATED( maxbox ) )      DEALLOCATE( maxbox )
      IF ( ALLOCATED( qsave ) )       DEALLOCATE( qsave )
      IF ( ALLOCATED( boxradius ) )   DEALLOCATE( boxradius )
      IF ( ALLOCATED( xyz ) )         DEALLOCATE( xyz )
      IF ( ALLOCATED( spher ) )       DEALLOCATE( spher )
      !
    END SUBROUTINE deallocatenewdreal
    !
    !------------------------------------------------------------------------
    SUBROUTINE qpointlist()
      !------------------------------------------------------------------------
      !
      ! ... This subroutine is the driver routine of the box system in this 
      ! ... implementation of US in real space.
      ! ... All the variables common in the module are computed and stored for 
      ! ... reusing. 
      ! ... This routine has to be called every time the atoms are moved and of
      ! ... course at the beginning.
      ! ... A set of spherical boxes are computed for each atom. 
      ! ... In boxradius there are the radii of the boxes.
      ! ... In maxbox the upper limit of leading index, namely the number of
      ! ... points of the fine mesh contained in each box.
      ! ... In xyz there are the coordinates of the points with origin in the
      ! ... centre of atom.
      ! ... In boxdistance the distance from the centre.
      ! ... In spher the spherical harmonics computed for each box
      ! ... In qsave the q value interpolated in these boxes.
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi, fpi, eps8, eps16
      USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
      USE uspp,       ONLY : okvan
      USE uspp,       ONLY : indv, nhtol, nhtolm, ap, nhtoj, lpx, lpl
      USE uspp_param, ONLY : lmaxq, nh, nhm, tvanp, kkbeta, nbeta, &
                             qfunc, dion, lmaxkb, qfcoef, nqf, nqlc, &
                             lll, rinner
      USE atom,       ONLY : r
      USE pfft,       ONLY : npp
      USE mp_global,  ONLY : me_pool
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      INTEGER               :: qsdim, ia, ib, ibia, it, mbia
      INTEGER               :: indm, inbrx1, inbrx2, idimension, &
                               ilm, ih, jh, iih, lllnbnt, lllmbnt
      INTEGER               :: roughestimate, goodestimate, lamx2, l, nt
      INTEGER,  ALLOCATABLE :: bufferpoints(:,:)
      REAL(DP), ALLOCATABLE :: bufferdistance(:,:), bufferxyz(:,:,:)
      REAL(DP)              :: distance, interpolatedqtot, first, second
      INTEGER               :: index0, index, indproc, ir
      INTEGER               :: i, j, k, i0, j0, k0, ipol, lm, nb, mb, ilast
      REAL(DP)              :: posi(3)
      REAL(DP), ALLOCATABLE :: tau0(:,:), rl(:,:), rl2(:)
      REAL(DP), ALLOCATABLE :: tempspher(:,:), qtot(:,:,:), &
                               xsp(:), ysp(:), wsp(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
      !
      !
      IF ( .NOT. okvan ) RETURN
      !
      CALL start_clock( 'qpointlist' )
      !
      ! ... qsave is deallocated here to free the memory for the buffers
      !
      IF ( ALLOCATED( qsave ) ) DEALLOCATE( qsave )
      !
      IF ( .NOT. ALLOCATED( boxradius ) ) THEN
         !
         ! ... here we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         !
         ALLOCATE( boxradius( ntyp ) )
         !
         boxradius(:) = 0.D0
         !
         DO it = 1, ntyp
            DO inbrx1 = 1, nbeta(it)
               DO inbrx2 = 1, nbeta(it)
                  DO indm = kkbeta(it), 1, -1
                     !
                     IF ( ABS( qfunc(indm,inbrx1,inbrx2,it) ) > eps16 ) THEN
                        !
                        boxradius(it) = MAX( r(indm,it), boxradius(it) )
                        !
                     END IF
                     !
                  END DO
               END DO
            END DO
         END DO
         !
      END IF
      !
      ! ... bring all the atomic positions on the first unit cell
      !
      ALLOCATE( tau0( 3, nat ) )
      !
      tau0(:,:) = tau(:,:)
      !
      CALL cryst_to_cart( nat, tau0, bg, -1 )
      !
      tau0(:,:) = tau0(:,:) - ANINT( tau0(:,:) )
      !
      CALL cryst_to_cart( nat, tau0, at, 1 )
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is provided here
      !
      mbr = MAXVAL( boxradius(:) )
      !
      mbx = mbr*&
            SQRT( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 ) / alat
      mby = mbr*&
            SQRT( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 ) / alat
      mbz = mbr*&
            SQRT( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 ) / alat
      !
      dmbx = 2*ANINT( mbx*nrx1 ) + 1
      dmby = 2*ANINT( mby*nrx2 ) + 1
      dmbz = 2*ANINT( mbz*nrx3 ) + 1
      !
      roughestimate = ANINT( DBLE( dmbx*dmby*dmbz ) * pi / 6.D0 )
      !
     ! PRINT *, roughestimate, 10*INT( ( (1/omega)*nrxx*( mbr**3 ) ) )
      !
      ALLOCATE( bufferpoints(   roughestimate, nat ) )
      ALLOCATE( bufferdistance( roughestimate, nat ) )
      ALLOCATE( bufferxyz( 3,   roughestimate, nat ) )
      !
      bufferpoints(:,:) = 0
      bufferdistance(:,:) = 0.D0
      !
      ! ... from soubroutine make_pointlist we have copied these lines and
      ! ... adapted to our needs; the mesh is multiplied by 27 to cover 
      ! ... also border atoms on properly.
      !
      IF ( .NOT.ALLOCATED( maxbox ) ) ALLOCATE( maxbox( nat ) )
      !
      maxbox(:) = 0
      !
      ! ... now we find the points
      !
      index0 = 0
      !
#if defined (__PARA)
      !
      DO i = 1, me_pool
         index0 = index0 + nrx1*nrx2*npp(i)
      END DO
      !
#endif
      !
      DO ia = 1, nat
         !
         DO ir = 1, nrxx
            !
            index = index0 + ir - 1
            k0    = index / (nrx1*nrx2)
            index = index - (nrx1*nrx2)*k0
            j0    = index / nrx1
            index = index - nrx1*j0
            i0    = index
            !
            DO i = i0-nr1, i0+nr1, nr1
               DO j = j0-nr2, j0+nr2, nr2
                  DO k = k0-nr3, k0+nr3, nr3
                     !
                     DO ipol = 1, 3
                        posi(ipol) = DBLE( i ) / DBLE( nr1 ) * at(ipol,1) + &
                                     DBLE( j ) / DBLE( nr2 ) * at(ipol,2) + &
                                     DBLE( k ) / DBLE( nr3 ) * at(ipol,3)
                     END DO
                     !
                     distance = SQRT( ( posi(1) - tau0(1,ia) )**2 + &
                                      ( posi(2) - tau0(2,ia) )**2 + &
                                      ( posi(3) - tau0(3,ia) )**2 )*alat
                     !
                     IF ( distance < boxradius(ityp(ia)) ) THEN
                        !
                        mbia = maxbox(ia) + 1
                        !
                        maxbox(ia) = mbia
                        bufferpoints(mbia,ia) = ir
                        bufferdistance(mbia,ia) = distance
                        bufferxyz(:,mbia,ia) = ( posi(:) - tau0(:,ia) )*alat
                        !
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      DEALLOCATE( tau0 )
      !
      goodestimate = MAXVAL( maxbox )
     ! PRINT *, "GOOD ESTIMATE = ", goodestimate
      !
      IF ( goodestimate > roughestimate ) &
         CALL errore( 'qpointlist', 'rough-estimate is too rough', 2 )
      !
      ! ... now store them in a more convenient place
      !
      IF ( ALLOCATED( box ) )         DEALLOCATE( box )
      IF ( ALLOCATED( boxdistance ) ) DEALLOCATE( boxdistance )
      IF ( ALLOCATED( xyz ) )         DEALLOCATE( xyz )
      !
      ALLOCATE( box(         goodestimate, nat ) )
      ALLOCATE( boxdistance( goodestimate, nat ) )
      ALLOCATE( xyz( 3,      goodestimate, nat ) )
      !
      box(:,:)         = bufferpoints(1:goodestimate,:)
      boxdistance(:,:) = bufferdistance(1:goodestimate,:)
      xyz(:,:,:)       = bufferxyz(:,1:goodestimate,:)
      !
      DEALLOCATE( bufferpoints )
      DEALLOCATE( bufferdistance )
      DEALLOCATE( bufferxyz )
      !
      ! ... now it computes the spherical harmonics
      !
      lamx2 = lmaxq*lmaxq
      !
      IF ( ALLOCATED( spher ) ) DEALLOCATE( spher )
      !
      ALLOCATE( spher( goodestimate, lamx2, nat ) )
      !
      spher(:,:,:) = 0.D0
      !
      DO ia = 1, nat
         !
         idimension = maxbox(ia)
         !
         ALLOCATE( rl( 3, idimension ), rl2( idimension ) )
         !
         DO ir = 1, idimension
            !
            rl(:,ir) = xyz(:,ir,ia)
            !
            rl2(ir) = rl(1,ir)**2 + rl(2,ir)**2 + rl(3,ir)**2
            !
         END DO
         !
         ALLOCATE( tempspher( idimension, lamx2 ) )
         !
         CALL ylmr2( lamx2, idimension, rl, rl2, tempspher )
         !
         spher(1:idimension,:,ia) = tempspher(:,:)
         !
         DEALLOCATE( rl, rl2, tempspher )
         !
      END DO
      !
      DEALLOCATE( xyz )
      !
      ! ... let's do the main work
      !
      qsdim = 0
      DO ia = 1, nat
         IF ( maxbox(ia) == 0 ) CYCLE
         nt = ityp(ia)
         IF ( .NOT. tvanp(nt) ) CYCLE
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
               qsdim = qsdim + 1
            END DO
         END DO
      END DO
      !      
      PRINT *, "QSAVE SIZE : ", goodestimate*qsdim
      !
      ALLOCATE( qsave( goodestimate, qsdim ) )
      !
      qsave(:,:) = 0.D0
      !
      ! ... the source is inspired by init_us_1
      !
      ! ... we perform two steps: first we compute for each l the qtot
      ! ... (radial q), then we interpolate it in our mesh, and then we
      ! ... add it to qsave with the correct spherical harmonics
      !
      ! ... Q is read from pseudo and it is divided into two parts:
      ! ... in the inner radius a polinomial representation is known and so
      ! ... strictly speaking we do not use interpolation but just compute
      ! ... the correct value
      !
      ib   = 0
      ibia = 0
      !
      DO ia = 1, nat
         !
         IF ( maxbox(ia) == 0 ) CYCLE
         !
         nt = ityp(ia)
         !
         IF ( .NOT. tvanp(nt) ) CYCLE
         !
         ALLOCATE( qtot( kkbeta(nt), nbeta(nt), nbeta(nt) ) )
         !
         ! ... variables used for spline interpolation
         !
         ALLOCATE( xsp( kkbeta(nt) ), ysp( kkbeta(nt) ), wsp( kkbeta(nt ) ) )
         !
         ! ... the radii in x
         !
         xsp(:) = r(1:kkbeta(nt),nt)
         !
         DO l = 0, nqlc(nt) - 1
            !
            ! ... first we build for each nb,mb,l the total Q(|r|) function
            ! ... note that l is the true (combined) angular momentum
            ! ... and that the arrays have dimensions 1..l+1
            !
            DO nb = 1, nbeta(nt)
               DO mb = nb, nbeta(nt)
                  !
                  lllnbnt = lll(nb,nt)
                  lllmbnt = lll(mb,nt)
                  !
                  IF ( .NOT. ( l >= ABS( lllnbnt - lllmbnt ) .AND. &
                               l <= lllnbnt + lllmbnt        .AND. &
                               MOD( l + lllnbnt + lllmbnt, 2 ) == 0 ) ) CYCLE
                  !
                  DO ir = 1, kkbeta(nt)
                     IF ( r(ir,nt) >= rinner(l+1,nt) ) THEN
                        qtot(ir,nb,mb) = qfunc(ir,nb,mb,nt) / r(ir,nt)**2
                     ELSE
                        ilast = ir
                     END IF
                  END DO
                  !
                  IF ( rinner(l+1,nt) > 0.D0 ) &
                     CALL setqfcorr( qfcoef(1,l+1,nb,mb,nt), &
                                     qtot(1,nb,mb), r(1,nt), nqf(nt), l, ilast )
                  !
                  ! ... we save the values in y
                  !
                  ysp(:) = qtot(1:kkbeta(nt),nb,mb)
                  !
                  ! ... compute the first derivative in first point
                  !
                  CALL setqfcorrpointfirst( qfcoef(1,l+1,nb,mb,nt), &
                                            first, r(1,nt), nqf(nt), l )
                  !
                  ! ... compute the second derivative in second point
                  !
                  CALL setqfcorrpointsecond( qfcoef(1,l+1,nb,mb,nt), &
                                             second, r(1,nt), nqf(nt), l )
                  !
                  ! ... call spline
                  !
                  CALL spline( xsp, ysp, first, second, wsp )
                  !
                  DO ir = 1, maxbox(ia)
                     !
                     IF ( boxdistance(ir,ia) < rinner(l+1,nt) ) THEN
                        !
                        ! ... if in the inner radius just compute the
                        ! ... polynomial
                        !
                        CALL setqfcorrpoint( qfcoef(1,l+1,nb,mb,nt), &
                                             interpolatedqtot,       &
                                             boxdistance(ir,ia), nqf(nt), l )
                        !
                     ELSE   
                        !
                        ! ... spline interpolation
                        !
                        interpolatedqtot = splint( xsp, ysp, &
                                                   wsp, boxdistance(ir,ia) )
                        !
                     END IF
                     !
                     ib = ibia
                     !
                     DO ih = 1, nh(nt)
                        DO jh = ih, nh(nt)
                           !
                           ib = ib + 1
                           !
                           IF ( .NOT.( nb == indv(ih,nt) .AND. &
                                       mb == indv(jh,nt) ) ) CYCLE
                           !
                           DO lm = l*l+1, (l+1)*(l+1)
                              !
                              qsave(ir,ib) = qsave(ir,ib) + &
                                             interpolatedqtot*spher(ir,lm,ia)*&
                                             ap(lm,nhtolm(ih,nt),nhtolm(jh,nt))
                              !
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         !
         ibia = ib
         !
         DEALLOCATE( qtot )
         DEALLOCATE( xsp )
         DEALLOCATE( ysp )
         DEALLOCATE( wsp )
         !
      END DO
      !
      DEALLOCATE( boxdistance )
      DEALLOCATE( spher )
      !
      CALL stop_clock( 'qpointlist' )
      !
    END SUBROUTINE qpointlist
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_r()
      !------------------------------------------------------------------------
      !
      ! ... this subroutine is the version of newd in real space
      !
      USE constants,        ONLY : pi, fpi
      USE ions_base,        ONLY : nat, ityp
      USE cell_base,        ONLY : omega
      USE gvect,            ONLY : nr1, nr2, nr3, nrxx
      USE lsda_mod,         ONLY : nspin
      USE scf,              ONLY : vr, vltot
      USE uspp,             ONLY : okvan
      USE uspp,             ONLY : deeq, deeq_nc, dvan, dvan_so
      USE uspp_param,       ONLY : nh, nhm, tvanp
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : so, domag, lspinorb
      !
      IMPLICIT NONE
      !
      REAL(DP), ALLOCATABLE :: aux(:,:)
      INTEGER               :: ia, ih, jh, ijh, is, ir, nt, nspin0
      INTEGER               :: ib, mbia, nht, nhnt
      !
      IF ( .NOT. okvan ) THEN
         !
         ! ... no ultrasoft potentials: use bare coefficients for projectors
         !
         DO ia = 1, nat
            !
            nt  = ityp(ia)
            nht = nh(nt)
            !
            IF ( lspinorb ) THEN
               !
               deeq_nc(1:nht,1:nht,ia,1:nspin) = dvan_so(1:nht,1:nht,1:nspin,nt)
               !
            ELSE IF ( noncolin ) THEN
               !
               deeq_nc(1:nht,1:nht,ia,1) = dvan(1:nht,1:nht,nt)
               deeq_nc(1:nht,1:nht,ia,2) = ( 0.D0, 0.D0 )
               deeq_nc(1:nht,1:nht,ia,3) = ( 0.D0, 0.D0 )
               deeq_nc(1:nht,1:nht,ia,4) = dvan(1:nht,1:nht,nt)
               !
            ELSE
               !
               DO is = 1, nspin
                  !
                  deeq(1:nht,1:nht,ia,is) = dvan(1:nht,1:nht,nt)
                  !
               END DO
               !
            END IF
            !
         END DO
         !
         ! ... early return
         !
         RETURN
         !
      END IF
      !
      CALL start_clock( 'newd' )
      !
      nspin0 = nspin
      !
      IF ( noncolin .AND..NOT. domag ) nspin0 = 1
      !
      deeq(:,:,:,:) = 0.D0
      !
      ALLOCATE( aux( nrxx, nspin0 ) )
      !
      DO is = 1, nspin0
         !
         IF ( nspin0 == 4 .AND. is /= 1 ) THEN
            aux(:,is) = vr(:,is)
         ELSE
            aux(:,is) = vltot(:) + vr(:,is)
         END IF
         !
      END DO
      !
      ib = 0
      !
      DO ia = 1, nat
         !
         mbia = maxbox(ia)
         !
         IF ( mbia == 0 ) CYCLE
         !
         nt = ityp(ia)
         !
         IF ( .NOT. tvanp(nt) ) CYCLE
         !
         nhnt = nh(nt)
         !
         DO ih = 1, nhnt
            DO jh = ih, nhnt
               ib = ib + 1
               DO is = 1, nspin0
                  DO ir = 1, mbia
                     deeq(ih,jh,ia,is)= deeq(ih,jh,ia,is) + &
                                        qsave(ir,ib)*aux(box(ir,ia),is)
                  END DO
                  deeq(jh,ih,ia,is) = deeq(ih,jh,ia,is)
               END DO
            END DO
         END DO
      END DO
      !
      deeq(:,:,:,:) = deeq(:,:,:,:)*omega/(nr1*nr2*nr3)
      !
      DEALLOCATE( aux )
      !
      CALL reduce( nhm*nhm*nat*nspin0, deeq )
      !
      DO ia = 1, nat
         !
         nt = ityp(ia)
         !
         IF ( noncolin ) THEN
            !
            IF ( so(nt) ) THEN
               CALL newd_so( ia )
            ELSE
               CALL newd_nc( ia )
            END IF
            !
         ELSE
            !
            nhnt = nh(nt)
            !
            DO is = 1, nspin0
               DO ih = 1, nhnt
                  DO jh = ih, nhnt
                     deeq(ih,jh,ia,is) = deeq(ih,jh,ia,is) + dvan(ih,jh,nt)
                     deeq(jh,ih,ia,is) = deeq(ih,jh,ia,is)
                  END DO
               END DO
            END DO
            !
         END IF
      ENDDO
      !
      CALL stop_clock( 'newd' )
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE newd_so( ia )
          !--------------------------------------------------------------------
          !
          USE spin_orb, ONLY : fcoef, so, domag, lspinorb
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: ia
          INTEGER             :: ijs, is1, is2, kh, lh
          !
          !
          nt = ityp(ia)
          ijs = 0
          !
          DO is1 = 1, 2
             DO is2 = 1, 2
                !
                ijs = ijs + 1
                !
                IF ( domag ) THEN
                   !
                   DO ih = 1, nh(nt)
                      DO jh = 1, nh(nt)
                         !
                         deeq_nc(ih,jh,ia,ijs) = dvan_so(ih,jh,ijs,nt)
                         !
                         DO kh = 1, nh(nt)
                            DO lh = 1, nh(nt)
                               !
                               deeq_nc(ih,jh,ia,ijs) = deeq_nc(ih,jh,ia,ijs) + &
                                deeq (kh,lh,ia,1)*                             &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) + &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                                deeq (kh,lh,ia,2)*                             &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt) + &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                (0.D0,-1.D0)*deeq (kh,lh,ia,3)*                &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt) - &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                deeq (kh,lh,ia,4)*                             &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) - &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))   
                               !
                            END DO
                         END DO
                      END DO
                   END DO
                   !
                ELSE
                   !
                   DO ih = 1, nh(nt)
                      DO jh = 1, nh(nt)
                         !
                         deeq_nc(ih,jh,ia,ijs) = dvan_so(ih,jh,ijs,nt)
                         !
                         DO kh = 1, nh(nt)
                            DO lh = 1, nh(nt)
                               !
                               deeq_nc(ih,jh,ia,ijs) = deeq_nc(ih,jh,ia,ijs) + &
                                deeq (kh,lh,ia,1)*                             &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) + &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt) ) 
                               !
                            END DO
                         END DO
                      END DO
                   END DO
                   !
                END IF
                !
             END DO
          END DO
          !
          RETURN
          !
        END SUBROUTINE newd_so
        !
        !--------------------------------------------------------------------
        SUBROUTINE newd_nc( ia )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: ia
          !
          nt = ityp(ia)
          !
          DO ih = 1, nh(nt)
             DO jh = 1, nh(nt)
                !
                IF ( lspinorb ) THEN
                   !
                   deeq_nc(ih,jh,ia,1) = dvan_so(ih,jh,1,nt) + &
                                         deeq(ih,jh,ia,1) + deeq(ih,jh,ia,4)
                   deeq_nc(ih,jh,ia,4) = dvan_so(ih,jh,4,nt) + &
                                         deeq(ih,jh,ia,1) - deeq(ih,jh,ia,4)
                   !
                ELSE
                   !
                   deeq_nc(ih,jh,ia,1) = dvan(ih,jh,nt) + &
                                         deeq(ih,jh,ia,1) + deeq(ih,jh,ia,4)
                   deeq_nc(ih,jh,ia,4) = dvan(ih,jh,nt) + &
                                         deeq(ih,jh,ia,1) - deeq(ih,jh,ia,4)
                   !
                END IF
                !
                deeq_nc(ih,jh,ia,2) = deeq(ih,jh,ia,2) - &
                                      ( 0.D0, 1.D0 ) * deeq(ih,jh,ia,3)
                !                      
                deeq_nc(ih,jh,ia,3) = deeq(ih,jh,ia,2) + &
                                      ( 0.D0, 1.D0 ) * deeq(ih,jh,ia,3)
                !                      
             END DO
          END DO
          !
          RETURN
          !
        END SUBROUTINE newd_nc
        !
    END SUBROUTINE newd_r
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorr( qfcoef, rho, r, nqf, ltot, mesh )
      !-----------------------------------------------------------------------
      !
      ! ... This routine compute the first part of the Q function up to rinner.
      ! ... On output it contains  Q
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN):: nqf, ltot, mesh
        ! input: the number of coefficients
        ! input: the angular momentum
        ! input: the number of mesh point
      REAL(DP), INTENT(IN) :: r(mesh), qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(OUT) :: rho(mesh)
        ! output: the function to be computed
      !
      INTEGER  :: ir, i
      REAL(DP) :: rr
      !
      DO ir = 1, mesh
         !
         rr = r(ir)**2
         !
         rho(ir) = qfcoef(1)
         !
         DO i = 2, nqf
            rho(ir) = rho(ir) + qfcoef(i)*rr**(i-1)
         END DO
         !
         rho(ir) = rho(ir)*r(ir)**ltot
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE setqfcorr
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrpoint( qfcoef, rho, r, nqf, ltot )
      !------------------------------------------------------------------------
      !
      ! ... This routine compute the first part of the Q function at the
      ! ... point r. On output it contains  Q
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN):: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(IN) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(OUT) :: rho 
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      !
      rho = qfcoef(1)
      !
      DO i = 2, nqf
         rho = rho + qfcoef(i)*rr**(i-1)
      END DO
      !
      rho = rho*r**ltot
      !
      RETURN
      !
    END SUBROUTINE setqfcorrpoint
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrpointfirst( qfcoef, rho, r, nqf, ltot )
      !------------------------------------------------------------------------
      !
      ! ... On output it contains  Q'  (probably wrong)
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(IN) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(OUT) :: rho 
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      !
      DO i = MAX( 1, 2-ltot ), nqf
         rho = rho + qfcoef(i)*rr**(i-2+ltot)*(i-1+ltot)
      END DO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrpointfirst
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrpointsecond( qfcoef, rho, r, nqf, ltot )
      !------------------------------------------------------------------------
      !
      ! ... On output it contains  Q''
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(IN) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(OUT) :: rho
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      !
      DO i = MAX( 3-ltot, 1 ), nqf
         rho = rho + qfcoef(i)*rr**(i-3+ltot)*(i-1+ltot)*(i-2+ltot)
      END DO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrpointsecond
    !
    !------------------------------------------------------------------------
    SUBROUTINE addusdens_r()
      !------------------------------------------------------------------------
      !
      ! ... This routine adds to the charge density the part which is due to
      ! ... the US augmentation.
      !
      USE ions_base,        ONLY : nat, ityp
      USE lsda_mod,         ONLY : nspin
      USE scf,              ONLY : rho
      USE uspp,             ONLY : okvan
      USE uspp,             ONLY : becsum
      USE uspp_param,       ONLY : tvanp, nh
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      !
      IMPLICIT NONE
      !
      INTEGER :: ia, ib, nt, ir, irb, ih, jh, ijh, is, nspin0, mbia, nhnt
      !
      !
      IF ( .NOT. okvan ) RETURN
      !
      CALL start_clock( 'addusdens' )
      !
      nspin0 = nspin
      !
      IF ( noncolin .AND..NOT. domag ) nspin0 = 1
      !
      DO is = 1, nspin0
         !
         ib = 0
         !
         DO ia = 1, nat
            !
            mbia = maxbox(ia)
            !
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            !
            IF ( .NOT. tvanp(nt) ) CYCLE
            !
            nhnt = nh(nt)
            !
            ijh = 0
            !
            DO ih = 1, nhnt
               DO jh = ih, nhnt
                  !
                  ijh = ijh + 1
                  ib  = ib  + 1
                  !
                  DO ir = 1, mbia
                     irb = box(ir,ia)
                     rho(irb,is) = rho(irb,is) + &
                                   qsave(ir,ib)*becsum(ijh,ia,is)
                  END DO
               END DO
            END DO
         END DO
         !
      END DO
      !
      CALL stop_clock( 'addusdens' )
      !
      RETURN
      !
    END SUBROUTINE addusdens_r
    !
END MODULE realus
