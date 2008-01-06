!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
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
  REAL(DP), ALLOCATABLE :: qsave(:)
  REAL(DP), ALLOCATABLE :: boxrad(:)
  REAL(DP), ALLOCATABLE :: boxdist(:,:), xyz(:,:,:)     
  REAL(DP), ALLOCATABLE :: spher(:,:,:)
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocatenewdreal()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( box ) )     DEALLOCATE( box )
      IF ( ALLOCATED( maxbox ) )  DEALLOCATE( maxbox )
      IF ( ALLOCATED( qsave ) )   DEALLOCATE( qsave )
      IF ( ALLOCATED( boxrad ) )  DEALLOCATE( boxrad )
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
      ! ... In boxdist the distance from the centre.
      ! ... In spher the spherical harmonics computed for each box
      ! ... In qsave the q value interpolated in these boxes.
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi, fpi, eps8, eps16
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
      USE uspp,       ONLY : okvan, indv, nhtol, nhtolm, ap, nhtoj, lpx, lpl
      USE uspp_param, ONLY : upf, lmaxq, nh, nhm
      USE atom,       ONLY : rgrid
      USE fft_base,   ONLY : dfftp
      USE mp_global,  ONLY : me_pool
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      INTEGER               :: qsdim, ia, mbia, iqs, iqsia
      INTEGER               :: indm, idimension, &
                               ih, jh, ijh, lllnbnt, lllmbnt
      INTEGER               :: roughestimate, goodestimate, lamx2, l, nt
      INTEGER,  ALLOCATABLE :: buffpoints(:,:)
      REAL(DP), ALLOCATABLE :: buffdist(:,:)
      REAL(DP)              :: distsq, qtot_int, first, second
      INTEGER               :: idx0, idx, ir
      INTEGER               :: i, j, k, ipol, lm, nb, mb, ijv, ilast
      REAL(DP)              :: posi(3)
      REAL(DP), ALLOCATABLE :: rl(:,:), rl2(:)
      REAL(DP), ALLOCATABLE :: tempspher(:,:), qtot(:,:,:), &
                               xsp(:), ysp(:), wsp(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
      REAL(DP)              :: inv_nr1, inv_nr2, inv_nr3, tau_ia(3), boxradsq_ia
      !
      !
      IF ( .NOT. okvan ) RETURN
      !
      CALL start_clock( 'realus' )
      !
      ! ... qsave is deallocated here to free the memory for the buffers
      !
      IF ( ALLOCATED( qsave ) ) DEALLOCATE( qsave )
      !
      IF ( .NOT. ALLOCATED( boxrad ) ) THEN
         !
         ! ... here we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         !
         ALLOCATE( boxrad( nsp ) )
         !
         boxrad(:) = 0.D0
         !
         DO nt = 1, nsp
            DO ijv = 1, upf(nt)%nbeta*(upf(nt)%nbeta+1)/2 
               DO indm = upf(nt)%kkbeta, 1, -1
                  !
                  IF ( ABS( upf(nt)%qfunc(indm,ijv) ) > eps16 ) THEN
                     !
                     boxrad(nt) = MAX( rgrid(nt)%r(indm), boxrad(nt) )
                     !
                     CYCLE
                     !
                  END IF
                  !
               END DO
            END DO
         END DO
         !
         boxrad(:) = boxrad(:) / alat
         !
      END IF
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is provided here
      !
      mbr = MAXVAL( boxrad(:) )
      !
      mbx = mbr*SQRT( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
      mby = mbr*SQRT( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
      mbz = mbr*SQRT( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
      !
      dmbx = 2*ANINT( mbx*nrx1 ) + 2
      dmby = 2*ANINT( mby*nrx2 ) + 2
      dmbz = 2*ANINT( mbz*nrx3 ) + 2
      !
      roughestimate = ANINT( DBLE( dmbx*dmby*dmbz ) * pi / 6.D0 )
      !
      CALL start_clock( 'realus:boxes' )
      !
      ALLOCATE( buffpoints( roughestimate, nat ) )
      ALLOCATE( buffdist(   roughestimate, nat ) )
      !
      ALLOCATE( xyz( 3, roughestimate, nat ) )
      !
      buffpoints(:,:) = 0
      buffdist(:,:) = 0.D0
      !
      IF ( .NOT.ALLOCATED( maxbox ) ) ALLOCATE( maxbox( nat ) )
      !
      maxbox(:) = 0
      !
      ! ... now we find the points
      !
      idx0 = 0
      !
#if defined (__PARA)
      !
      DO i = 1, me_pool
         idx0 = idx0 + nrx1*nrx2*dfftp%npp(i)
      END DO
      !
#endif
      !
      inv_nr1 = 1.D0 / DBLE( nr1 )
      inv_nr2 = 1.D0 / DBLE( nr2 )
      inv_nr3 = 1.D0 / DBLE( nr3 )
      !
      DO ia = 1, nat
         !
         nt = ityp(ia)
         !
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         !
         boxradsq_ia = boxrad(nt)**2
         !
         tau_ia(1) = tau(1,ia)
         tau_ia(2) = tau(2,ia)
         tau_ia(3) = tau(3,ia)
         !
         DO ir = 1, nrxx
            !
            ! ... three dimensional indexes
            !
            idx   = idx0 + ir - 1
            k     = idx / (nrx1*nrx2)
            idx   = idx - (nrx1*nrx2)*k
            j     = idx / nrx1
            idx   = idx - nrx1*j
            i     = idx
            !
            DO ipol = 1, 3
               posi(ipol) = DBLE( i )*inv_nr1*at(ipol,1) + &
                            DBLE( j )*inv_nr2*at(ipol,2) + &
                            DBLE( k )*inv_nr3*at(ipol,3)
            END DO
            !
            posi(:) = posi(:) - tau_ia(:)
            !
            ! ... minimum image convenction
            !
            CALL cryst_to_cart( 1, posi, bg, -1 )
            !
            posi(:) = posi(:) - ANINT( posi(:) )
            !
            CALL cryst_to_cart( 1, posi, at, 1 )
            !
            distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
            !
            IF ( distsq < boxradsq_ia ) THEN
               !
               mbia = maxbox(ia) + 1
               !
               maxbox(ia)          = mbia
               buffpoints(mbia,ia) = ir
               buffdist(mbia,ia)   = SQRT( distsq )*alat
               xyz(:,mbia,ia)      = posi(:)*alat
               !
            END IF
         END DO
      END DO
      !
      goodestimate = MAXVAL( maxbox )
      !
      IF ( goodestimate > roughestimate ) &
         CALL errore( 'qpointlist', 'rough-estimate is too rough', 2 )
      !
      ! ... now store them in a more convenient place
      !
      IF ( ALLOCATED( box ) )     DEALLOCATE( box )
      IF ( ALLOCATED( boxdist ) ) DEALLOCATE( boxdist )
      !
      ALLOCATE( box(     goodestimate, nat ) )
      ALLOCATE( boxdist( goodestimate, nat ) )
      !
      box(:,:)     = buffpoints(1:goodestimate,:)
      boxdist(:,:) = buffdist(1:goodestimate,:)
      !
      DEALLOCATE( buffpoints )
      DEALLOCATE( buffdist )
      !
      CALL stop_clock( 'realus:boxes' )
      CALL start_clock( 'realus:spher' )
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
         nt = ityp(ia)
         !
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
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
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:qsave' )
      !
      ! ... let's do the main work
      !
      qsdim = 0
      DO ia = 1, nat
         mbia = maxbox(ia)
         IF ( mbia == 0 ) CYCLE
         nt = ityp(ia)
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
               qsdim = qsdim + mbia
            END DO
         END DO
      END DO
      !      
!      PRINT *, "QSAVE SIZE : ", qsdim
      !
      ALLOCATE( qsave( qsdim ) )
      !
      qsave(:) = 0.D0
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
      iqs   = 0
      iqsia = 0
      !
      DO ia = 1, nat
         !
         mbia = maxbox(ia)
         !
         IF ( mbia == 0 ) CYCLE
         !
         nt = ityp(ia)
         !
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         !
         ALLOCATE( qtot( upf(nt)%kkbeta, upf(nt)%nbeta, upf(nt)%nbeta ) )
         !
         ! ... variables used for spline interpolation
         !
         ALLOCATE( xsp( upf(nt)%kkbeta ), ysp( upf(nt)%kkbeta ), &
                   wsp( upf(nt)%kkbeta ) )
         !
         ! ... the radii in x
         !
         xsp(:) = rgrid(nt)%r(1:upf(nt)%kkbeta)
         !
         DO l = 0, upf(nt)%nqlc - 1
            !
            ! ... first we build for each nb,mb,l the total Q(|r|) function
            ! ... note that l is the true (combined) angular momentum
            ! ... and that the arrays have dimensions 1..l+1
            !
            DO nb = 1, upf(nt)%nbeta
               DO mb = nb, upf(nt)%nbeta
                  ijv = mb * (mb-1) /2 + nb
                  !
                  lllnbnt = upf(nt)%lll(nb)
                  lllmbnt = upf(nt)%lll(mb)
                  !
                  IF ( .NOT. ( l >= ABS( lllnbnt - lllmbnt ) .AND. &
                               l <= lllnbnt + lllmbnt        .AND. &
                               MOD( l + lllnbnt + lllmbnt, 2 ) == 0 ) ) CYCLE
                  !
                  DO ir = 1, upf(nt)%kkbeta
                     IF ( rgrid(nt)%r(ir) >= upf(nt)%rinner(l+1) ) THEN
                        qtot(ir,nb,mb) = upf(nt)%qfunc(ir,ijv) / &
                                         rgrid(nt)%r(ir)**2
                     ELSE
                        ilast = ir
                     END IF
                  END DO
                  !
                  IF ( upf(nt)%rinner(l+1) > 0.D0 ) &
                     CALL setqfcorr( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                        qtot(1,nb,mb), rgrid(nt)%r(1), upf(nt)%nqf, l, ilast )
                  !
                  ! ... we save the values in y
                  !
                  ysp(:) = qtot(1:upf(nt)%kkbeta,nb,mb)
                  !
                  ! ... compute the first derivative in first point
                  !
                  CALL setqfcorrptfirst( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                                   first, rgrid(nt)%r(1), upf(nt)%nqf, l )
                  !
                  ! ... compute the second derivative in second point
                  !
                  CALL setqfcorrptsecond( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                                   second, rgrid(nt)%r(1), upf(nt)%nqf, l )
                  !
                  ! ... call spline
                  !
                  CALL spline( xsp, ysp, first, second, wsp )
                  !
                  DO ir = 1, maxbox(ia)
                     !
                     IF ( boxdist(ir,ia) < upf(nt)%rinner(l+1) ) THEN
                        !
                        ! ... if in the inner radius just compute the
                        ! ... polynomial
                        !
                        CALL setqfcorrpt( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                                   qtot_int, boxdist(ir,ia), upf(nt)%nqf, l )
                        !
                     ELSE   
                        !
                        ! ... spline interpolation
                        !
                        qtot_int = splint( xsp, ysp, wsp, boxdist(ir,ia) )
                        !
                     END IF
                     !
                     ijh = 0
                     !
                     DO ih = 1, nh(nt)
                        DO jh = ih, nh(nt)
                           !
                           iqs = iqsia + ijh*mbia + ir
                           ijh = ijh + 1
                           !
                           IF ( .NOT.( nb == indv(ih,nt) .AND. &
                                       mb == indv(jh,nt) ) ) CYCLE
                           !
                           DO lm = l*l+1, (l+1)*(l+1)
                              !
                              qsave(iqs) = qsave(iqs) + &
                                           qtot_int*spher(ir,lm,ia)*&
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
         iqsia = iqs
         !
         DEALLOCATE( qtot )
         DEALLOCATE( xsp )
         DEALLOCATE( ysp )
         DEALLOCATE( wsp )
         !
      END DO
      !
      DEALLOCATE( boxdist )
      DEALLOCATE( spher )
      !
      CALL stop_clock( 'realus:qsave' )
      CALL stop_clock( 'realus' )
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
      USE scf,              ONLY : v, vltot
      USE uspp,             ONLY : okvan, deeq, deeq_nc, dvan, dvan_so
      USE uspp_param,       ONLY : upf, nh, nhm
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : so, domag, lspinorb
      !
      IMPLICIT NONE
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER               :: ia, ih, jh, is, ir, nt, nspin0
      INTEGER               :: mbia, nht, nhnt, iqs
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
      ALLOCATE( aux( nrxx ) )
      !
      DO is = 1, nspin0
         !
         IF ( nspin0 == 4 .AND. is /= 1 ) THEN
            aux(:) = v%of_r(:,is)
         ELSE
            aux(:) = vltot(:) + v%of_r(:,is)
         END IF
         !
         iqs = 0
         !
         DO ia = 1, nat
            !
            mbia = maxbox(ia)
            !
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            !
            IF ( .NOT. upf(nt)%tvanp ) CYCLE
            !
            nhnt = nh(nt)
            !
            DO ih = 1, nhnt
               DO jh = ih, nhnt
                  DO ir = 1, mbia
                     iqs = iqs + 1
                     deeq(ih,jh,ia,is)= deeq(ih,jh,ia,is) + &
                                        qsave(iqs)*aux(box(ir,ia))
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
    SUBROUTINE setqfcorrpt( qfcoef, rho, r, nqf, ltot )
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
    END SUBROUTINE setqfcorrpt
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrptfirst( qfcoef, rho, r, nqf, ltot )
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
      rho = 0.D0
      !
      DO i = MAX( 1, 2-ltot ), nqf
         rho = rho + qfcoef(i)*rr**(i-2+ltot)*(i-1+ltot)
      END DO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrptfirst
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrptsecond( qfcoef, rho, r, nqf, ltot )
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
      rho = 0.D0
      !
      DO i = MAX( 3-ltot, 1 ), nqf
         rho = rho + qfcoef(i)*rr**(i-3+ltot)*(i-1+ltot)*(i-2+ltot)
      END DO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrptsecond
    !
    !------------------------------------------------------------------------
    SUBROUTINE addusdens_r()
      !------------------------------------------------------------------------
      !
      ! ... This routine adds to the charge density the part which is due to
      ! ... the US augmentation.
      !
      USE ions_base,        ONLY : nat, ityp
      USE cell_base,        ONLY : omega
      USE lsda_mod,         ONLY : nspin
      USE scf,              ONLY : rho
      USE klist,            ONLY : nelec
      USE gvect,            ONLY : nr1, nr2, nr3
      USE uspp,             ONLY : okvan, becsum
      USE uspp_param,       ONLY : upf, nh
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      !
      IMPLICIT NONE
      !
      INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, nspin0, mbia, nhnt, iqs
      REAL(DP) :: charge
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
         iqs = 0
         !
         DO ia = 1, nat
            !
            mbia = maxbox(ia)
            !
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            !
            IF ( .NOT. upf(nt)%tvanp ) CYCLE
            !
            nhnt = nh(nt)
            !
            ijh = 0
            !
            DO ih = 1, nhnt
               DO jh = ih, nhnt
                  !
                  ijh = ijh + 1
                  !
                  DO ir = 1, mbia
                     !
                     irb = box(ir,ia)
                     iqs = iqs + 1
                     !
                     rho%of_r(irb,is) = rho%of_r(irb,is) + qsave(iqs)*becsum(ijh,ia,is)
                  END DO
               END DO
            END DO
         END DO
         !
      END DO
      !
      ! ... check the integral of the total charge
      !
      charge = SUM( rho%of_r(:,1:nspin0) )*omega / ( nr1*nr2*nr3 )
      !
      CALL reduce( 1, charge )
      !
      IF ( ABS( charge - nelec ) / charge > 1.D-4 ) THEN
         !
         ! ... the error on the charge is too large
         !
         CALL errore( 'addusdens_r', 'charge is wrong: increase ecutrho', 1 )
         !
      ELSE
         !
         ! ... rescale the density to impose the correct number of electrons
         !
         rho%of_r(:,:) = rho%of_r(:,:) / charge * nelec
         !
      END IF
      !
      CALL stop_clock( 'addusdens' )
      !
      RETURN
      !
    END SUBROUTINE addusdens_r
    !
END MODULE realus
