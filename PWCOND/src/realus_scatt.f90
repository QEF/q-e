!
! Copyright (C) 2009 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE realus_scatt
!
! Some extra subroutines to the module realus
! needed for the scattering problem
!
 INTEGER,  ALLOCATABLE :: orig_or_copy(:,:)

 CONTAINS

 SUBROUTINE realus_scatt_0()
!
! Calculates orig_or_copy array
!
   USE kinds,            ONLY : dp
   USE constants,        ONLY : pi
   USE ions_base,        ONLY : nat, tau, ityp
   USE cell_base,        ONLY : at, bg
   USE realus,           ONLY : qpointlist, tabp, boxrad
   USE uspp,             ONLY : okvan
   USE uspp_param,       ONLY : upf
   USE mp_global,        ONLY : me_pool
   USE fft_base,         ONLY : dfftp

   IMPLICIT NONE

   INTEGER  :: ia, ir, mbia, roughestimate, idx0, idx, i, j, k, i_lr, ipol
   REAL(DP) :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz, distsq
   REAL(DP) :: inv_nr1, inv_nr2, inv_nr3, boxradsq_ia, posi(3)

   IF ( .NOT. okvan ) RETURN

   CALL qpointlist(dfftp, tabp)

!--   Finds roughestimate
   mbr = MAXVAL( boxrad(:) )
   mbx = mbr*SQRT( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
   mby = mbr*SQRT( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
   mbz = mbr*SQRT( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
   dmbx = 2*ANINT( mbx*dfftp%nr1x ) + 2
   dmby = 2*ANINT( mby*dfftp%nr2x ) + 2
   dmbz = 2*ANINT( mbz*dfftp%nr3x ) + 2
   roughestimate = ANINT( DBLE( dmbx*dmby*dmbz ) * pi / 6.D0 )
!--

   IF (ALLOCATED(orig_or_copy)) DEALLOCATE( orig_or_copy )
   ALLOCATE( orig_or_copy( roughestimate, nat ) )

   ! idx0 = starting index of real-space FFT arrays for this processor
   idx0 = dfftp%nr1x*dfftp%nr2x * dfftp%ipp(me_pool+1)

   inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
   inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
   inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )

   DO ia = 1, nat
       IF ( .NOT. upf(ityp(ia))%tvanp ) CYCLE
       mbia = 0
       boxradsq_ia = boxrad(ityp(ia))**2
       DO ir = 1, dfftp%nr1x*dfftp%nr2x * dfftp%npl
         idx   = idx0 + ir - 1
         k     = idx / (dfftp%nr1x*dfftp%nr2x)
         idx   = idx - (dfftp%nr1x*dfftp%nr2x)*k
         j     = idx / dfftp%nr1x
         idx   = idx - dfftp%nr1x*j
         i     = idx
         DO ipol = 1, 3
           posi(ipol) = DBLE( i )*inv_nr1*at(ipol,1) + &
                        DBLE( j )*inv_nr2*at(ipol,2) + &
                        DBLE( k )*inv_nr3*at(ipol,3)
         END DO
         posi(:) = posi(:) - tau(:,ia)
         CALL cryst_to_cart( 1, posi, bg, -1 )
         IF ( abs(ANINT(posi(3))).gt.1.d-6 ) THEN
           i_lr = 0
         ELSE
           i_lr = 1
         END IF
         posi(:) = posi(:) - ANINT( posi(:) )
         CALL cryst_to_cart( 1, posi, at, 1 )
         distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
         IF ( distsq < boxradsq_ia ) THEN
            mbia = mbia + 1
            orig_or_copy(mbia,ia) = i_lr
         END IF
       END DO
   END DO

   RETURN
 END SUBROUTINE realus_scatt_0

 SUBROUTINE realus_scatt_1(becsum_orig)
!
! Augments the charge and spin densities.
!
   USE kinds,            ONLY : dp
   USE ions_base,        ONLY : nat, ityp
   USE lsda_mod,         ONLY : nspin
   USE scf,              ONLY : rho
   USE realus,           ONLY : tabp
   USE uspp,             ONLY : okvan, becsum, ijtoh
   USE uspp_param,       ONLY : upf, nhm, nh
   USE noncollin_module, ONLY : noncolin
   USE spin_orb,         ONLY : domag

   IMPLICIT NONE

   INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, nspin0, mbia, nhnt, iqs
   REAL(DP) :: becsum_orig(nhm*(nhm+1)/2,nat,nspin)

   IF (.NOT.okvan) RETURN

   nspin0 = nspin
   IF (noncolin.AND..NOT.domag) nspin0 = 1
   DO is = 1, nspin0
     iqs = 0
     DO ia = 1, nat
        mbia = tabp(ia)%maxbox
        IF ( mbia == 0 ) CYCLE
        nt = ityp(ia)
        IF ( .NOT. upf(nt)%tvanp ) CYCLE
        nhnt = nh(nt)
        ijh = 0
        DO ih = 1, nhnt
           DO jh = ih, nhnt
              ijh = ijh + 1
              DO ir = 1, mbia
                 irb = tabp(ia)%box(ir)
                 iqs = iqs + 1
                 if(orig_or_copy(ir,ia).eq.1) then
                  rho%of_r(irb,is) = rho%of_r(irb,is) + tabp(ia)%qr(ir,ijtoh(ih,jh,nt))*becsum_orig(ijh,ia,is)
                 else
                  rho%of_r(irb,is) = rho%of_r(irb,is) + tabp(ia)%qr(ir,ijtoh(ih,jh,nt))*becsum(ijh,ia,is)
                 endif
              ENDDO
           ENDDO
        ENDDO
     ENDDO
   ENDDO
   RETURN
 END SUBROUTINE realus_scatt_1

END MODULE realus_scatt

