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
  ! ... modified by O. Baris Malcioglu (2008)
  ! ... TODO : Write the k points part
  INTEGER,  ALLOCATABLE :: box(:,:), maxbox(:)
  REAL(DP), ALLOCATABLE :: qsave(:)
  REAL(DP), ALLOCATABLE :: boxrad(:)
  REAL(DP), ALLOCATABLE :: boxdist(:,:), xyz(:,:,:)     
  REAL(DP), ALLOCATABLE :: spher(:,:,:)
  ! Beta function in real space 
  INTEGER,  ALLOCATABLE :: box_beta(:,:), maxbox_beta(:) 
  REAL(DP), ALLOCATABLE :: betasave(:,:,:)
  REAL(DP), ALLOCATABLE :: boxrad_beta(:)
  REAL(DP), ALLOCATABLE :: boxdist_beta(:,:), xyz_beta(:,:,:)     
  REAL(DP), ALLOCATABLE :: spher_beta(:,:,:)
  !General
  !LOGICAL               :: tnlr        ! old hidden variable, should be removed soon
  LOGICAL               :: real_space  ! When this flag is true, real space versions of the corresponding 
                                       ! calculations are performed
  INTEGER               :: real_space_debug ! remove this, for debugging purposes
  INTEGER               :: initialisation_level ! init_realspace_vars sets this to 3 qpointlist adds 5
                                                ! betapointlist adds 7 so the value should be 15 if the
                                                ! real space routine is initalised properly
                                            
  integer, allocatable :: &
       igk_k(:,:),&                   ! The g<->k correspondance for each k point 
       npw_k(:)                       ! number of plane waves at each k point
                           ! They are (used many times, it is much better to hold them in memory 
  
  !REAL(DP), ALLOCATABLE :: psic_rs(:) !In order to prevent mixup, a redundant copy of psic`
  !
  complex(DP), allocatable :: tg_psic(:)
  complex(DP), allocatable :: psic_temp(:),tg_psic_temp(:) !Copies of psic and tg_psic
  complex(DP), allocatable :: tg_vrs(:) !task groups linear V memory
  complex(DP), allocatable :: psic_box_temp(:),tg_psic_box_temp(:)
  !
  CONTAINS
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE init_realspace_vars()
     !---------------------------------------------------------------------------
     !This subroutine should be called to allocate/reset real space related variables. 
     !---------------------------------------------------------------------------
     use wvfct,                only : npwx,npw, igk, g2kin
     use klist,                only : nks,xk
     use gvect,                only : ngm, g, ecutwfc
     use cell_base,            only : tpiba2
     USE gvect,                ONLY : nrxx
     USE control_flags,        ONLY : tqr, use_task_groups
     USE fft_base,             ONLY : dffts
     USE mp_global,            ONLY : nogrp
     use wavefunctions_module, only : psic
     USE io_global,            ONLY : stdout
     
     
     implicit none
     
     integer :: ik
     
     print *, "<<<<<init_realspace_vars>>>>>>>"
     
     IF ( ALLOCATED( igk_k ) )     DEALLOCATE( igk_k )
     IF ( ALLOCATED( npw_k ) )     DEALLOCATE( npw_k )
     
     allocate(igk_k(npwx,nks))
     allocate(npw_k(nks))
     !allocate (psic_temp(size(psic)))
     !real space, allocation for task group fft work arrays
     IF( use_task_groups ) THEN
        !
        If (allocated( tg_psic ) ) deallocate( tg_psic ) 
        !
        ALLOCATE( tg_psic( dffts%nnrx * nogrp ) )
        !ALLOCATE( tg_psic_temp( dffts%nnrx * nogrp ) )
        ALLOCATE( tg_vrs( dffts%nnrx * nogrp ) )

        !
     ENDIF
     !allocate (psic_rs( nrxx))
     !at this point I can not decide if I should preserve a redundant copy of the real space psi, or transform it whenever required, 
     do ik=1,nks
      !
      CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
      !
      npw_k(ik) = npw
      !
      igk_k(:,ik) = igk(:)
      !
     !
     enddo
     
     !tqr = .true.
     initialisation_level = initialisation_level + 7
     if (real_space_debug > 20 .and. real_space_debug < 30) then 
       real_space=.false.
       if (tqr) then
         tqr = .false.
         write(stdout,'("Debug level forced tqr to be set false")')
       else
         write(stdout,'("tqr was already set false")')
       endif  
       real_space_debug=real_space_debug-20
     endif
     !print *, "Real space = ", real_space
     !print *, "Real space debug ", real_space_debug 

     
    END SUBROUTINE init_realspace_vars
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
      REAL(DP), ALLOCATABLE :: rl(:,:), rl2(:), d1y(:), d2y(:)
      REAL(DP), ALLOCATABLE :: tempspher(:,:), qtot(:,:,:), &
                               xsp(:), ysp(:), wsp(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz, aux
      REAL(DP)              :: inv_nr1, inv_nr2, inv_nr3, tau_ia(3), boxradsq_ia
      !
      !
      initialisation_level = 3
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
            IF ( .NOT. upf(nt)%tvanp ) CYCLE
            DO ijv = 1, upf(nt)%nbeta*(upf(nt)%nbeta+1)/2 
               DO indm = upf(nt)%mesh,1,-1
                  !
                  IF( upf(nt)%q_with_l ) THEN
                     aux = SUM(ABS( upf(nt)%qfuncl(indm,ijv,:) ))
                  ELSE
                     aux = ABS( upf(nt)%qfunc(indm,ijv) )
                  ENDIF
                  IF ( aux > eps16 ) THEN
                     !
                     boxrad(nt) = MAX( rgrid(nt)%r(indm), boxrad(nt) )
                     !
                     EXIT
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
#if defined (__PARA)
      idx0 = nrx1*nrx2 * SUM ( dfftp%npp(1:me_pool) )
#else
      idx0 = 0
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
            ! ... three dimensional indices (i,j,k)
            !
            idx   = idx0 + ir - 1
            k     = idx / (nrx1*nrx2)
            idx   = idx - (nrx1*nrx2)*k
            j     = idx / nrx1
            idx   = idx - nrx1*j
            i     = idx
            !
            ! ... do not include points outside the physical range!
            !
            IF ( i >= nr1 .or. j >= nr2 .or. k >= nr3 ) CYCLE
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
                  IF( upf(nt)%q_with_l ) THEN
                      qtot(1:upf(nt)%kkbeta,nb,mb) = &
                          upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l) &
                           / rgrid(nt)%r(1:upf(nt)%kkbeta)**2
                  ELSE
                      DO ir = 1, upf(nt)%kkbeta
                        IF ( rgrid(nt)%r(ir) >= upf(nt)%rinner(l+1) ) THEN
                            qtot(ir,nb,mb) = upf(nt)%qfunc(ir,ijv) / &
                                            rgrid(nt)%r(ir)**2
                        ELSE
                            ilast = ir
                        END IF
                      END DO
                  ENDIF
                  !
                  IF ( upf(nt)%rinner(l+1) > 0.D0 ) &
                     CALL setqfcorr( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                        qtot(1,nb,mb), rgrid(nt)%r, upf(nt)%nqf, l, ilast )
                  !
                  ! ... we save the values in y
                  !
                  ysp(:) = qtot(1:upf(nt)%kkbeta,nb,mb)
                  !
                  IF ( upf(nt)%nqf > 0 ) THEN
                      !
                      ! ... compute the first derivative in first point
                      !
                      CALL setqfcorrptfirst( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                                      first, rgrid(nt)%r(1), upf(nt)%nqf, l )
                      !
                      ! ... compute the second derivative in first point
                      !
                      CALL setqfcorrptsecond( upf(nt)%qfcoef(1:,l+1,nb,mb), &
                                      second, rgrid(nt)%r(1), upf(nt)%nqf, l )
                  ELSE
                      !
                      ! ... if we don't have the analitical coefficients, try to do
                      ! ... the same numerically (note that setting first=0.d0 and 
                      ! ... second=0.d0 makes almost no difference)
                      !
                      ALLOCATE( d1y(upf(nt)%kkbeta), d2y(upf(nt)%kkbeta) )
                      CALL radial_gradient(ysp(1:upf(nt)%kkbeta), d1y, &
                                           rgrid(nt)%r, upf(nt)%kkbeta, 1)
                      CALL radial_gradient(d1y, d2y, rgrid(nt)%r, upf(nt)%kkbeta, 1)
                      !
                      first = d1y(1) ! first derivative in first point
                      second =d2y(1) ! second derivative in first point
                      DEALLOCATE( d1y, d2y )
                  ENDIF
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
    SUBROUTINE betapointlist()
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
      USE constants,  ONLY : pi, eps8, eps16
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE gsmooth,    ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
      USE uspp,       ONLY : okvan, indv, nhtol, nhtolm, ap
      USE uspp_param, ONLY : upf, lmaxq, nh, nhm
      USE atom,       ONLY : rgrid
      !USE pffts,      ONLY : npps
      USE fft_base,   ONLY : dffts
      USE mp_global,  ONLY : me_pool
      USE splinelib,  ONLY : spline, splint
      USE ions_base,             ONLY : ntyp => nsp
      !
      IMPLICIT NONE
      !
      INTEGER               :: betasdim, ia, it, mbia, iqs
      INTEGER               :: indm, inbrx, idimension, &
                               ilm, ih, jh, iih, ijh
      INTEGER               :: roughestimate, goodestimate, lamx2, l, nt
      INTEGER,  ALLOCATABLE :: buffpoints(:,:)
      REAL(DP), ALLOCATABLE :: buffdist(:,:)
      REAL(DP)              :: distsq, qtot_int, first, second
      INTEGER               :: index0, index, indproc, ir
      INTEGER               :: i, j, k, i0, j0, k0, ipol, lm, nb, mb, ijv, ilast
      REAL(DP)              :: posi(3)
      REAL(DP), ALLOCATABLE :: rl(:,:), rl2(:)
      REAL(DP), ALLOCATABLE :: tempspher(:,:), qtot(:,:,:), &
                               xsp(:), ysp(:), wsp(:), d1y(:), d2y(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
      REAL(DP)              :: inv_nr1s, inv_nr2s, inv_nr3s, tau_ia(3), boxradsq_ia
      !Delete Delete
      character(len=256) :: filename
      character(len=256) :: tmp
      !Delete Delete
      !
      !
      initialisation_level = initialisation_level + 5
      IF ( .NOT. okvan ) RETURN
      !
      print *, "<<<betapointlist>>>"
      !
      CALL start_clock( 'betapointlist' )
      !
      ! ... qsave is deallocated here to free the memory for the buffers
      !
      IF ( ALLOCATED( betasave ) ) DEALLOCATE( betasave )
      !
      IF ( .NOT. ALLOCATED( boxrad_beta ) ) THEN
         !
         ! ... here we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         !
         ALLOCATE( boxrad_beta( nsp ) )
         !
         boxrad_beta(:) = 0.D0
         !
         DO it = 1, nsp
            DO inbrx = 1, upf(it)%nbeta
               DO indm = upf(it)%kkbeta, 1, -1
                  !
                  IF ( ABS( upf(it)%beta(indm,inbrx) ) > 0.d0 ) THEN
                     !
                     boxrad_beta(it) = MAX( rgrid(it)%r(indm), boxrad_beta(it) )
                     !
                     CYCLE
                     !
                  END IF
                  !
               END DO
            END DO
         END DO
         !
         boxrad_beta(:) = boxrad_beta(:) / alat
         !
      END IF
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is provided here
      !
      mbr = MAXVAL( boxrad_beta(:) )
      !
      mbx = mbr*SQRT( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
      mby = mbr*SQRT( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
      mbz = mbr*SQRT( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
      !
      dmbx = 2*ANINT( mbx*nrx1s ) + 2
      dmby = 2*ANINT( mby*nrx2s ) + 2
      dmbz = 2*ANINT( mbz*nrx3s ) + 2
      !
      roughestimate = ANINT( DBLE( dmbx*dmby*dmbz ) * pi / 6.D0 )
      !
      CALL start_clock( 'realus:boxes' )
      !
      ALLOCATE( buffpoints( roughestimate, nat ) )
      ALLOCATE( buffdist(   roughestimate, nat ) )
      !
      ALLOCATE( xyz_beta( 3, roughestimate, nat ) )
      !
      buffpoints(:,:) = 0
      buffdist(:,:) = 0.D0
      !
      IF ( .NOT.ALLOCATED( maxbox_beta ) ) ALLOCATE( maxbox_beta( nat ) )
      !
      maxbox_beta(:) = 0
      !
      ! ... now we find the points
      !
      index0 = 0
      !
#if defined (__PARA)
      !
      DO i = 1, me_pool
         index0 = index0 + nrx1s*nrx2s*dffts%npp(i)
      END DO
      !
#endif
      !
      inv_nr1s = 1.D0 / DBLE( nr1s )
      inv_nr2s = 1.D0 / DBLE( nr2s )
      inv_nr3s = 1.D0 / DBLE( nr3s )
      !
      DO ia = 1, nat
         !
         IF ( .NOT. upf(ityp(ia))%tvanp ) CYCLE
         !
         boxradsq_ia = boxrad_beta(ityp(ia))**2
         !
         tau_ia(1) = tau(1,ia)
         tau_ia(2) = tau(2,ia)
         tau_ia(3) = tau(3,ia)
         !
         DO ir = 1, nrxxs
            !
            ! ... three dimensional indexes
            !
            index = index0 + ir - 1
            k     = index / (nrx1s*nrx2s)
            index = index - (nrx1s*nrx2s)*k
            j     = index / nrx1s
            index = index - nrx1s*j
            i     = index
            !
            DO ipol = 1, 3
               posi(ipol) = DBLE( i )*inv_nr1s*at(ipol,1) + &
                            DBLE( j )*inv_nr2s*at(ipol,2) + &
                            DBLE( k )*inv_nr3s*at(ipol,3)
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
               mbia = maxbox_beta(ia) + 1
               !
               maxbox_beta(ia)     = mbia
               buffpoints(mbia,ia) = ir
               buffdist(mbia,ia)   = SQRT( distsq )*alat
               xyz_beta(:,mbia,ia) = posi(:)*alat
               !
            END IF
         END DO
      END DO
      !
      goodestimate = MAXVAL( maxbox_beta )
      !
      IF ( goodestimate > roughestimate ) &
         CALL errore( 'betapointlist', 'rough-estimate is too rough', 2 )
      !
      ! ... now store them in a more convenient place
      !
      IF ( ALLOCATED( box_beta ) )     DEALLOCATE( box_beta )
      IF ( ALLOCATED( boxdist_beta ) ) DEALLOCATE( boxdist_beta )
      !
      ALLOCATE( box_beta    ( goodestimate, nat ) )
      ALLOCATE( boxdist_beta( goodestimate, nat ) )
      !
      box_beta(:,:)     = buffpoints(1:goodestimate,:)
      boxdist_beta(:,:) = buffdist(1:goodestimate,:)
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
      IF ( ALLOCATED( spher_beta ) ) DEALLOCATE( spher_beta )
      !
      ALLOCATE( spher_beta( goodestimate, lamx2, nat ) )
      !
      spher_beta(:,:,:) = 0.D0
      !
      DO ia = 1, nat
         !
         IF ( .NOT. upf(ityp(ia))%tvanp ) CYCLE
         !
         idimension = maxbox_beta(ia)
         !
         ALLOCATE( rl( 3, idimension ), rl2( idimension ) )
         !
         DO ir = 1, idimension
            !
            rl(:,ir) = xyz_beta(:,ir,ia)
            !
            rl2(ir) = rl(1,ir)**2 + rl(2,ir)**2 + rl(3,ir)**2
            !
         END DO
         !
         ALLOCATE( tempspher( idimension, lamx2 ) )
         !
         CALL ylmr2( lamx2, idimension, rl, rl2, tempspher )
         !
         spher_beta(1:idimension,:,ia) = tempspher(:,:)
         !
         DEALLOCATE( rl, rl2, tempspher )
         !
      END DO
      !
      DEALLOCATE( xyz_beta )
      !
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:qsave' )
      !
      ! ... let's do the main work
      !
      betasdim = 0
      DO ia = 1, nat
         mbia = maxbox_beta(ia)
         IF ( mbia == 0 ) CYCLE
         nt = ityp(ia)
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         DO ih = 1, nh(nt)
               betasdim = betasdim + mbia
         END DO
      END DO
      !      
      !PRINT *, "BETASAVE SIZE : ", betasdim
      !
      ALLOCATE( betasave( nat, nhm, goodestimate )  )
      !
      betasave = 0.D0
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
      !
      DO ia = 1, nat
         !
         mbia = maxbox_beta(ia)
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
         ALLOCATE( xsp( upf(nt)%kkbeta ), ysp( upf(nt)%kkbeta ), wsp( upf(nt)%kkbeta ) )
         !
         ! ... the radii in x
         !
         xsp(:) = rgrid(nt)%r(1:upf(nt)%kkbeta)
         !
         DO ih = 1, nh (nt)
            !
            lm = nhtolm(ih, nt)
            nb = indv(ih, nt)
            !
            ! The following is the proper way of interpolating to grid, 
            ! suggested by Lorenzo Paulatto, however does not work for 
            ! me at the moment. I will debug as soon as possible, please
            ! do not remove. 
            !ysp(:) = upf(nt)%beta(1:upf(nt)%kkbeta,nb)
            !ALLOCATE( d1y(upf(nt)%kkbeta), d2y(upf(nt)%kkbeta) )
            !CALL radial_gradient(ysp(1:upf(nt)%kkbeta), d1y, &
            !                     rgrid(nt)%r, upf(nt)%kkbeta, 1)
            !CALL radial_gradient(d1y, d2y, rgrid(nt)%r, upf(nt)%kkbeta, 1)
            !
            !first = d1y(1) ! first derivative in first point
            !second =d2y(1) ! second derivative in first point
            !DEALLOCATE( d1y, d2y )
            !
            !
            !CALL spline( xsp, ysp, first, second, wsp )
            !
            !end of lorenzo interpolation
            !  
            ! This is a homegrown solution to interpolation problem
            ! should be replaced soon
            !
            !OBM rgrid(nt)%r(1) == 0 ????? attempting correction
            if (rgrid(nt)%r(1)==0) then 
             ysp(2:) = upf(nt)%beta(2:upf(nt)%kkbeta,nb) / rgrid(nt)%r(2:upf(nt)%kkbeta)
             ysp(1)=0.d0
            else
             ysp(:) = upf(nt)%beta(1:upf(nt)%kkbeta,nb) / rgrid(nt)%r(1:upf(nt)%kkbeta)
            endif
            !print *, "ysp1",ysp(1),"=",upf(nt)%beta(1,nb),"//",rgrid(nt)%r(1)
            
            !
            first = (ysp(2) - ysp(1)) / (xsp(2) - xsp(1))
            !print *,"first",first
!            first = 0.d0
            !
            CALL spline( xsp, ysp, first, 0.d0, wsp )
            ! end of OBM workaround
 
            DO ir = 1, mbia
               !
               ! ... spline interpolation
               !
               qtot_int = splint( xsp, ysp, wsp, boxdist_beta(ir,ia) )
               !
               iqs = iqs + 1
               !
               betasave(ia,ih,ir) = qtot_int*spher_beta(ir,lm,ia)
               !print *, "qtot check=",qtot_int
               !
            END DO
         END DO
         !
         DEALLOCATE( qtot )
         DEALLOCATE( xsp )
         DEALLOCATE( ysp )
         DEALLOCATE( wsp )
         !
      END DO
      !
      DEALLOCATE( boxdist_beta )
      DEALLOCATE( spher_beta )
      !
      CALL stop_clock( 'realus:qsave' )
      CALL stop_clock( 'betapointlist' )
! Delete Delete Delete Delete
!          print *,"Dumping Betasave v1"
!          DO ia = 1, nat           
!           nt = ityp(ia)
!           print *, "atom ",ia," has ", nh (nt), " beta functions"
!           DO ih = 1, nh (nt)
!           write(tmp,fmt='("at",i1,"beta",i1)'),ia,ih
!           filename="rsv1-" // trim(tmp) // ".dump"
!           OPEN(UNIT=47,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
!             DO ir = 1, nrxxs
!                !
!                ! ... three dimensional indexes
!                !
!               index = ir - 1
!               k     = index / (nrx1s*nrx2s)
!               index = index - (nrx1s*nrx2s)*k
!               j     = index / nrx1s
!               index = index - nrx1s*j
!               i     = index
!               !
!               
!               DO ipol = 1, 3
!                 posi(ipol) = DBLE( i )*inv_nr1s*at(ipol,1) + &
!                 DBLE( j )*inv_nr2s*at(ipol,2) + &
!                 DBLE( k )*inv_nr3s*at(ipol,3)
!               END DO
!               posi(:) = posi(:)*alat*0.52917
!               !if (posi(3) /= 0 ) cycle
!               DO jh = 1, mbia                    
!                  IF(box_beta(jh,ia) .eq. ir) then
!                    write(unit=47,FMT='(3(F8.4,"  "),E14.5)'),posi(1),posi(2),posi(3),betasave(ia,ih,jh)
!                  endif
!               ENDDO
!             ENDDO
!           close(47)
!          ENDDO
!         ENDDO
!
!         
!         
!
! Delete Delete Delete Delete
      !
    END SUBROUTINE betapointlist

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
      USE noncollin_module, ONLY : noncolin, nspin_mag
      USE spin_orb,         ONLY : domag, lspinorb
      USE mp_global,        ONLY : intra_pool_comm
      USE mp,               ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER               :: ia, ih, jh, is, ir, nt
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
      deeq(:,:,:,:) = 0.D0
      !
      ALLOCATE( aux( nrxx ) )
      !
      DO is = 1, nspin_mag
         !
         IF ( nspin_mag == 4 .AND. is /= 1 ) THEN
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
      CALL mp_sum(  deeq(:,:,:,1:nspin_mag) , intra_pool_comm )
      !
      DO ia = 1, nat
         !
         nt = ityp(ia)
         !
         IF ( noncolin ) THEN
            !
            IF ( upf(nt)%has_so ) THEN
               CALL newd_so( ia )
            ELSE
               CALL newd_nc( ia )
            END IF
            !
         ELSE
            !
            nhnt = nh(nt)
            !
            DO is = 1, nspin_mag
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
          USE spin_orb, ONLY : fcoef, domag, lspinorb
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
      USE noncollin_module, ONLY : noncolin, nspin_mag
      USE spin_orb,         ONLY : domag
      USE mp_global,        ONLY : intra_pool_comm
      USE mp,               ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, mbia, nhnt, iqs
      CHARACTER(LEN=80) :: msg
      REAL(DP) :: charge
      !
      !
      IF ( .NOT. okvan ) RETURN
      !
      CALL start_clock( 'addusdens' )
      !
      DO is = 1, nspin_mag
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
      charge = SUM( rho%of_r(:,1:nspin_mag) )*omega / ( nr1*nr2*nr3 )
      !
      CALL mp_sum(  charge , intra_pool_comm )
      !
      IF ( ABS( charge - nelec ) / charge > 1.D-4 ) THEN
         !
         ! ... the error on the charge is too large
         !
         WRITE (msg,'("expected ",f13.8,", found ",f13.8)') &
            nelec, charge
         CALL errore( 'addusdens_r', &
            TRIM(msg)//': wrong charge, increase ecutrho', 1 )
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
    !--------------------------------------------------------------------------
    SUBROUTINE calbec_rs_gamma ( ibnd, m, rbecp )
    
  !--------------------------------------------------------------------------
  ! 
  ! Subroutine written by Dario Rocca Stefano de Gironcoli, modified by O. Baris Malcioglu 
  !
  ! Calculates rbecp in real space
  ! Requires BETASAVE (the beta functions at real space) calculated by betapointlist() (added to realus) 
  ! ibnd is an index that runs over the number of bands, which is given by m
  ! So you have to call this subroutine inside a cycle with index ibnd
  ! In this cycle you have to perform a Fourier transform of the orbital
  ! corresponding to ibnd, namely you have to transform the orbital to
  ! real space and store it in the global variable psic. 
  ! Remember that in the gamma_only case you
  ! perform two fast Fourier transform at the same time, and so you have
  ! that the real part correspond to ibnd, and the imaginary part to ibnd+1
  ! 
  ! WARNING: For the sake of speed, there are no checks performed in this routine, check beforehand!
    USE kinds,                 ONLY : DP
    USE cell_base,             ONLY : omega
    USE wavefunctions_module,  ONLY : psic
    USE ions_base,             ONLY : nat, ntyp => nsp, ityp
    USE gsmooth,               ONLY : nr1s, nr2s, nr3s
    USE uspp_param,            ONLY : nh, nhm
    !USE becmod,                ONLY : rbecp
    USE control_flags,         ONLY : use_task_groups
    USE task_groups,           ONLY : tg_gather
    USE mp_global,             ONLY : nogrp, ogrp_comm, me_pool, nolist
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in) :: ibnd, m
    INTEGER :: iqs, iqsp, ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), allocatable, dimension(:) :: wr, wi
    REAL(DP) :: bcr, bci
    REAL(DP), dimension(:,:), intent(out) :: rbecp
    !COMPLEX(DP), allocatable, dimension(:) :: bt
    !integer :: ir, k
    !
    REAL(DP), external :: DDOT
    !
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
  
     CALL errore( 'calbec_rs_gamma', 'task_groups not implemented', 1 )  
  
    ELSE !non task groups part starts here

    fac = SQRT(omega) / (nr1s*nr2s*nr3s)
    !
    rbecp(:,ibnd)=0.d0
    IF ( ibnd+1 .le. m ) rbecp(:,ibnd+1)=0.d0
    ! Clearly for an odd number of bands for ibnd=nbnd=m you don't have
    ! anymore bands, and so the imaginary part equal zero
    !
       !
       iqs = 1
       ikb = 0
       !
       DO nt = 1, ntyp
          !
           DO ia = 1, nat
             !
             IF ( ityp(ia) == nt ) then
                !
                mbia = maxbox_beta(ia)
                
                ! maxbox_beta contains the maximum number of real space points necessary
                ! to describe the beta function corresponding to the atom ia
                ! Namely this is the number of grid points for which beta is
                ! different from zero
                !
                ALLOCATE( wr(mbia), wi(mbia) )
                ! just working arrays
                !
                DO ih = 1, nh(nt)
                   ! nh is the number of beta functions, or something similar
                   !
                   ikb = ikb + 1
                   iqsp = iqs+mbia-1
                   wr(:) = DBLE ( psic( box_beta(1:mbia,ia) ) )
                   wi(:) = AIMAG( psic( box_beta(1:mbia,ia) ) )
                   !print *, "betasave check", betasave(ia,ih,:)
                   ! box_beta contains explictly the points of the real space grid in
                   ! which the beta functions are differet from zero. Remember
                   ! that dble(psic) corresponds to ibnd, and aimag(psic) to ibnd+1:
                   ! this is the standard way to perform fourier transform in pwscf
                   ! in the gamma_only case
                   bcr  = DDOT( mbia, betasave(ia,ih,:), 1, wr(:) , 1 )
                   bci  = DDOT( mbia, betasave(ia,ih,:), 1, wi(:) , 1 )
                   ! in the previous two lines the real space integral is performed, using
                   ! few points of the real space mesh only
                   rbecp(ikb,ibnd)   = fac * bcr
                   IF ( ibnd+1 .le. m ) rbecp(ikb,ibnd+1) = fac * bci
                   ! It is necessary to multiply by fac which to obtain the integral in real
                   ! space
                   !print *, rbecp(ikb,ibnd)
                   iqs = iqsp + 1
                   !
                END DO
                !
                DEALLOCATE( wr, wi )
                !
             END IF
             !
          END DO
          !
       END DO
       !
       ! 
    ENDIF
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN 

  END SUBROUTINE calbec_rs_gamma
    !
    SUBROUTINE calbec_rs_k ( ibnd, m )
  !--------------------------------------------------------------------------
  ! The k_point generalised version of calbec_rs_gamma. Basically same as above, but becp is used instead
  ! of rbecp, skipping the gamma point reduction
  ! derived from above by OBM 051108
    USE kinds,                 ONLY : DP
    USE cell_base,             ONLY : omega
    USE wavefunctions_module,  ONLY : psic
    USE ions_base,             ONLY : nat, ntyp => nsp, ityp
    USE gsmooth,               ONLY : nr1s, nr2s, nr3s
    USE uspp_param,            ONLY : nh, nhm
    USE becmod,                ONLY : becp
    USE control_flags,         ONLY : use_task_groups
    USE task_groups,           ONLY : tg_gather
    USE mp_global,             ONLY : nogrp, ogrp_comm, me_pool, nolist
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in) :: ibnd, m
    INTEGER :: iqs, iqsp, ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), allocatable, dimension(:) :: wr, wi
    REAL(DP) :: bcr, bci
    !REAL(DP), dimension(:,:), intent(out) :: rbecp
    !COMPLEX(DP), allocatable, dimension(:) :: bt
    !integer :: ir, k
    !
    REAL(DP), external :: DDOT
    !
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
  
     CALL errore( 'calbec_rs_k', 'task_groups not implemented', 1 )  
  
    ELSE !non task groups part starts here

    fac = SQRT(omega) / (nr1s*nr2s*nr3s)
    !
    becp(:,ibnd)=0.d0
       iqs = 1
       ikb = 0
       !
       DO nt = 1, ntyp
          !
           DO ia = 1, nat
             !
             IF ( ityp(ia) == nt ) then
                !
                mbia = maxbox_beta(ia)
                ALLOCATE( wr(mbia), wi(mbia) )
                DO ih = 1, nh(nt)
                   ! nh is the number of beta functions, or something similar
                   !
                   ikb = ikb + 1
                   iqsp = iqs+mbia-1
                   wr(:) = DBLE ( psic( box_beta(1:mbia,ia) ) )
                   wi(:) = AIMAG( psic( box_beta(1:mbia,ia) ) )
                   bcr  = DDOT( mbia, betasave(ia,ih,:), 1, wr(:) , 1 )
                   bci  = DDOT( mbia, betasave(ia,ih,:), 1, wi(:) , 1 )
                   becp(ikb,ibnd)   = fac * CMPLX( bcr, bci)
                   iqs = iqsp + 1
                   !
                END DO
                !
                DEALLOCATE( wr, wi )
                !
             END IF
             !
          END DO
          !
       END DO
       !
       ! 
    ENDIF
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN 

  END SUBROUTINE calbec_rs_k
    !--------------------------------------------------------------------------
    SUBROUTINE s_psir_gamma ( ibnd, m )
  !--------------------------------------------------------------------------
  !  
  ! ... This routine applies the S matrix to m wavefunctions psi in real space (in psic), 
  ! ... and puts the results again in psic for backtransforming.
  ! ... Requires rbecp (calbecr in REAL SPACE) and betasave (from betapointlist in realus)
  ! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu
  ! WARNING ! for the sake of speed, no checks performed in this subroutine

      USE kinds,                  ONLY : DP
      USE cell_base,              ONLY : omega
      USE wavefunctions_module,   ONLY : psic
      USE ions_base,              ONLY : nat, ntyp => nsp, ityp
      USE uspp_param,             ONLY : nh
      USE lsda_mod,               ONLY : current_spin
      USE uspp,                   ONLY : qq
      USE becmod,                 ONLY : rbecp
      USE control_flags,          ONLY : use_task_groups
      USE task_groups,            ONLY : tg_gather
      USE mp_global,              ONLY : nogrp, ogrp_comm, me_pool, nolist
      !
      IMPLICIT NONE
      !
      INTEGER, intent(in) :: ibnd, m
      !
      INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      REAL(DP), allocatable, dimension(:) :: w1, w2, bcr, bci
      !
      real(DP), external :: DDOT
      !
       
      
      
      CALL start_clock( 's_psir' )
      IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
   
        CALL errore( 's_psir_gamma', 'task_groups not implemented', 1 )  
  
      ELSE !non task groups part starts here

      !
      fac = sqrt(omega)
      !
      ikb = 0
      iqs = 0
      jqs = 0
      !
      DO nt = 1, ntyp
         !
         DO ia = 1, nat
            !
            IF ( ityp(ia) == nt ) THEN
               !
               mbia = maxbox_beta(ia)
               !print *, "mbia=",mbia
               ALLOCATE( w1(nh(nt)),  w2(nh(nt)) )
               w1 = 0.D0
               w2 = 0.D0
               !
               DO ih = 1, nh(nt)
                  !
                  DO jh = 1, nh(nt)
                     !
                     jkb = ikb + jh
                     w1(ih) = w1(ih) + qq(ih,jh,nt) * rbecp(jkb, ibnd) 
                     IF ( ibnd+1 .le. m ) w2(ih) = w2(ih) + qq(ih,jh,nt) * rbecp(jkb, ibnd+1) 
                     !
                  END DO
                  !
               END DO
               !
               w1 = w1 * fac
               w2 = w2 * fac
               ikb = ikb + nh(nt)
               !
               DO ih = 1, nh(nt)
                  !
                  DO ir = 1, mbia
                     !
                     iqs = jqs + ir
                     psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*CMPLX( w1(ih), w2(ih) )
                     !
                  END DO
                  !
                  jqs = iqs
                  !
               END DO
               !
               DEALLOCATE( w1, w2 )
               !
            END IF
            !
         END DO
         !
      END DO
      !
      ENDIF
      CALL stop_clock( 's_psir' )
      !
      RETURN
      !
  END SUBROUTINE s_psir_gamma
  !
  SUBROUTINE s_psir_k ( ibnd, m )
  !--------------------------------------------------------------------------
  ! Same as s_psir_gamma but for generalised k point scheme i.e.:
  ! 1) Only one band is considered at a time 
  ! 2) Becp is a complex entity now
  ! Derived from s_psir_gamma by OBM 061108
      USE kinds,                  ONLY : DP
      USE cell_base,              ONLY : omega
      USE wavefunctions_module,   ONLY : psic
      USE ions_base,              ONLY : nat, ntyp => nsp, ityp
      USE uspp_param,             ONLY : nh
      USE lsda_mod,               ONLY : current_spin
      USE uspp,                   ONLY : qq
      USE becmod,                 ONLY : becp
      USE control_flags,          ONLY : use_task_groups
      USE task_groups,            ONLY : tg_gather
      USE mp_global,              ONLY : nogrp, ogrp_comm, me_pool, nolist
      !
      IMPLICIT NONE
      !
      INTEGER, intent(in) :: ibnd, m
      !
      INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      REAL(DP), allocatable, dimension(:) :: bcr, bci
      COMPLEX(DP) , allocatable, dimension(:) :: w1
      !
      real(DP), external :: DDOT
      !
       
      
      
      CALL start_clock( 's_psir' )
      IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
   
        CALL errore( 's_psir_k', 'task_groups not implemented', 1 )  
  
      ELSE !non task groups part starts here

      !
      fac = sqrt(omega)
      !
      ikb = 0
      iqs = 0
      jqs = 0
      !
      DO nt = 1, ntyp
         !
         DO ia = 1, nat
            !
            IF ( ityp(ia) == nt ) THEN
               !
               mbia = maxbox_beta(ia)
               ALLOCATE( w1(nh(nt)) )
               w1 = 0.D0
               !
               DO ih = 1, nh(nt)
                  !
                  DO jh = 1, nh(nt)
                     !
                     jkb = ikb + jh
                     w1(ih) = w1(ih) + qq(ih,jh,nt) * becp(jkb, ibnd) 
                     !
                  END DO
                  !
               END DO
               !
               w1 = w1 * fac
               ikb = ikb + nh(nt)
               !
               DO ih = 1, nh(nt)
                  !
                  DO ir = 1, mbia
                     !
                     iqs = jqs + ir
                     psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*w1(ih)
                     !
                  END DO
                  !
                  jqs = iqs
                  !
               END DO
               !
               DEALLOCATE( w1 )
               !
            END IF
            !
         END DO
         !
      END DO
      !
      ENDIF
      CALL stop_clock( 's_psir' )
      !
      RETURN
      !
  END SUBROUTINE s_psir_k
  ! 
  SUBROUTINE add_vuspsir_gamma ( ibnd, m )
  !--------------------------------------------------------------------------
  !
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector transformed in real space contained in psic. 
  !    ibnd is an index that runs over the number of bands, which is given by m
  !    Requires the products of psi with all beta functions
  !    in array rbecp(nkb,m) (calculated by calbecr in REAL SPACE)
  ! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu
  ! WARNING ! for the sake of speed, no checks performed in this subroutine

  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE wavefunctions_module,   ONLY : psic
  USE ions_base,              ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,             ONLY : nh
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq
  USE becmod,                 ONLY : rbecp, becp
  USE control_flags,          ONLY : use_task_groups
  USE task_groups,            ONLY : tg_gather
  USE mp_global,              ONLY : nogrp, ogrp_comm, me_pool, nolist
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: ibnd, m
  !
  INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  REAL(DP), allocatable, dimension(:) :: w1, w2, bcr, bci
  ! 
  real(DP), external :: DDOT
  !
  CALL start_clock( 'add_vuspsir' )
  
  IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
  
    CALL errore( 'add_vuspsir_gamma', 'task_groups not implemented', 1 )  
  
  ELSE !non task groups part starts here

   !
   fac = sqrt(omega)
   !
   ikb = 0
   iqs = 0
   jqs = 0
   !
   DO nt = 1, ntyp
      !
      DO ia = 1, nat
         !
         IF ( ityp(ia) == nt ) THEN
            !
            mbia = maxbox_beta(ia)
            ALLOCATE( w1(nh(nt)),  w2(nh(nt)) )
            w1 = 0.D0
            w2 = 0.D0
            !
            DO ih = 1, nh(nt)
               !
               DO jh = 1, nh(nt)
                  !
                  jkb = ikb + jh
                  !
                  w1(ih) = w1(ih) + deeq(ih,jh,ia,current_spin) * rbecp(jkb,ibnd) 
                  IF ( ibnd+1 .le. m )  w2(ih) = w2(ih) + deeq(ih,jh,ia,current_spin) * rbecp(jkb,ibnd+1) 
                  !
               END DO
               !
            END DO
            !
            w1 = w1 * fac
            w2 = w2 * fac
            ikb = ikb + nh(nt)
            !
            DO ih = 1, nh(nt)
               !
               DO ir = 1, mbia
                  ! 
                  iqs = jqs + ir
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*CMPLX( w1(ih), w2(ih) )
                  ! 
               END DO
                  !
               jqs = iqs
               !
            END DO
            !
            DEALLOCATE( w1, w2 )
            !
         END IF
         !
      END DO
      !
   END DO
   !
  END IF 
  CALL stop_clock( 'add_vuspsir' )
  !
  RETURN
  !
  END SUBROUTINE add_vuspsir_gamma
  !
  SUBROUTINE add_vuspsir_k ( ibnd, m )
  !--------------------------------------------------------------------------
  !
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector transformed in real space contained in psic. 
  !    ibnd is an index that runs over the number of bands, which is given by m
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated by calbecr in REAL SPACE)
  ! Subroutine written by Stefano de Gironcoli, modified by O. Baris Malcioglu
  ! WARNING ! for the sake of speed, no checks performed in this subroutine
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE wavefunctions_module,   ONLY : psic
  USE ions_base,              ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,             ONLY : nh
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq
  USE becmod,                 ONLY : becp
  USE control_flags,          ONLY : use_task_groups
  USE task_groups,            ONLY : tg_gather
  USE mp_global,              ONLY : nogrp, ogrp_comm, me_pool, nolist
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: ibnd, m
  !
  INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  REAL(DP), allocatable, dimension(:) ::  bcr, bci
  ! 
  COMPLEX(DP), allocatable, dimension(:) :: w1
  !
  real(DP), external :: DDOT
  !
  CALL start_clock( 'add_vuspsir' )
  
  IF( ( use_task_groups ) .AND. ( m >= nogrp ) ) THEN
  
    CALL errore( 'add_vuspsir_k', 'task_groups not implemented', 1 )  
  
  ELSE !non task groups part starts here

   !
   fac = sqrt(omega)
   !
   ikb = 0
   iqs = 0
   jqs = 0
   !
   DO nt = 1, ntyp
      !
      DO ia = 1, nat
         !
         IF ( ityp(ia) == nt ) THEN
            !
            mbia = maxbox_beta(ia)
            ALLOCATE( w1(nh(nt)) )
            w1 = (0.d0, 0d0)
            !
            DO ih = 1, nh(nt)
               !
               DO jh = 1, nh(nt)
                  !
                  jkb = ikb + jh
                  !
                  w1(ih) = w1(ih) + deeq(ih,jh,ia,current_spin) * becp(jkb,ibnd) 
                  !
               END DO
               !
            END DO
            !
            w1 = w1 * fac
            ikb = ikb + nh(nt)
            !
            DO ih = 1, nh(nt)
               !
               DO ir = 1, mbia
                  ! 
                  iqs = jqs + ir
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*w1(ih)
                  ! 
               END DO
                  !
               jqs = iqs
               !
            END DO
            !
            DEALLOCATE( w1)
            !
         END IF
         !
      END DO
      !
   END DO
  ENDIF
  CALL stop_clock( 'add_vuspsir' )
  RETURN
  !
  END SUBROUTINE add_vuspsir_k
  
  !--------------------------------------------------------------------------
  SUBROUTINE fft_orbital_gamma (orbital, ibnd, nbnd, conserved)
   !--------------------------------------------------------------------------
    !
    ! OBM 241008
    ! This driver subroutine transforms the given orbital using fft and puts the result in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    ! orbital: the orbital to be transformed
    ! ibnd: band index
    ! nbnd: total number of bands
    use wavefunctions_module,     only : psic
    use gsmooth,                  only : nr1s,nr2s,nr3s,nrx1s,nrx2s,&
       nrx3s,nrxxs,nls,nlsm,doublegrid  
    USE kinds,         ONLY : DP
    USE fft_parallel,  ONLY : tg_cft3s
    USE fft_base,      ONLY : dffts
    USE control_flags, ONLY : use_task_groups
    USE task_groups,   ONLY : tg_gather
    USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist

    implicit none
     
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
    complex(DP),intent(in) :: orbital(:,:)
    logical, optional :: conserved !if this flag is true, the orbital is stored in temporary memory
   
    !integer :: ig
    
    !Internal temporary variables
    complex(DP) :: fp, fm,alpha

    integer :: i, j, incr, ierr, idx, ioff, nsiz
    logical :: use_tg
    
    complex(DP), allocatable :: psic_temp2(:)

    !Task groups
    !COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
    INTEGER :: recv_cnt( nogrp ), recv_displ( nogrp )
    INTEGER :: v_siz

    
   !The new task group version based on vloc_psi
   !print *, "->Real space"
   CALL start_clock( 'fft_orbital' )
        use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp )
        
        IF( use_tg ) THEN
        !

        tg_psic = (0.d0, 0.d0)
        ioff   = 0  !print *, "cft3s sign",sign
        !
        DO idx = 1, 2*nogrp, 2

           IF( idx + ibnd - 1 < nbnd ) THEN
              DO j = 1, npw_k(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =        orbital(j,idx+ibnd-1) + (0.0d0,1.d0) * orbital(j,idx+ibnd) 
                 tg_psic(nlsm(igk_k(j,1))+ioff) = CONJG( orbital(j,idx+ibnd-1) - (0.0d0,1.d0) * orbital(j,idx+ibnd) )
              END DO
           ELSE IF( idx + ibnd - 1 == nbnd ) THEN
              DO j = 1, npw_k(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =        orbital(j,idx+ibnd-1)
                 tg_psic(nlsm(igk_k(j,1))+ioff) = CONJG( orbital(j,idx+ibnd-1) )
              END DO
           END IF

           ioff = ioff + dffts%nnrx

        END DO
        !
        !
        call tg_cft3s ( tg_psic, dffts, 2, use_tg )
        !
        !
        if (present(conserved)) then
         if (conserved .eqv. .true.) then
          if (.not. allocated(tg_psic_temp)) ALLOCATE( tg_psic_temp( dffts%nnrx * nogrp ) )
          tg_psic_temp=tg_psic
         endif
        endif

     ELSE !Task groups not used
        !
        psic(:) = (0.d0, 0.d0)
           
!           alpha=(0.d0,1.d0)
!           if (ibnd .eq. nbnd) alpha=(0.d0,0.d0)
!           
!           allocate (psic_temp2(npw_k(1)))
!           call ZCOPY(npw_k(1),orbital(:, ibnd),1,psic_temp2,1)
!           call ZAXPY(npw_k(1),alpha,orbital(:, ibnd+1),1,psic_temp2,1)
!           psic (nls (igk_k(:,1)))=psic_temp2(:)
!           call ZAXPY(npw_k(1),(-2.d0,0.d0)*alpha,orbital(:, ibnd+1),1,psic_temp2,1)
!           psic (nlsm (igk_k(:,1)))=conjg(psic_temp2(:))
!           deallocate(psic_temp2)
           

        if (ibnd < nbnd) then
           ! two ffts at the same time
           !print *,"alpha=",alpha
           do j = 1, npw_k(1)
              psic (nls (igk_k(j,1))) =       orbital(j, ibnd) + (0.0d0,1.d0)*orbital(j, ibnd+1)
              psic (nlsm(igk_k(j,1))) = CONJG(orbital(j, ibnd) - (0.0d0,1.d0)*orbital(j, ibnd+1))
              !print *, nls (igk_k(j,1))
           enddo
           !CALL errore( 'fft_orbital_gamma', 'bye bye', 1 )  
        else
           do j = 1, npw_k(1)
              psic (nls (igk_k(j,1))) =       orbital(j, ibnd)
              psic (nlsm(igk_k(j,1))) = CONJG(orbital(j, ibnd))
           enddo
        end if
        !
        !
       call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
        !
        if (present(conserved)) then
         if (conserved .eqv. .true.) then
           if (.not. allocated(psic_temp) ) allocate (psic_temp(size(psic))) 
           call ZCOPY(size(psic),psic,1,psic_temp,1)
         endif
        endif

     END IF
    
    !if (.not. allocated(psic)) CALL errore( 'fft_orbital_gamma', 'psic not allocated', 2 )
    ! OLD VERSION ! Based on an algorithm found somewhere in the TDDFT codes, generalised to k points
    !
    !   psic(:) =(0.0d0,0.0d0)
    !   if(ibnd<nbnd) then
    !      do ig=1,npw_k(1)
    !         !
    !         psic(nls(igk_k(ig,1)))=orbital(ig,ibnd)+&
    !              (0.0d0,1.0d0)*orbital(ig,ibnd+1)
    !         psic(nlsm(igk_k(ig,1)))=conjg(orbital(ig,ibnd)-&
    !              (0.0d0,1.0d0)*orbital(ig,ibnd+1))
    !         !
    !      enddo
    !   else
    !      do ig=1,npw_k(1)
    !         !
    !         psic(nls(igk_k(ig,1)))=orbital(ig,ibnd)
    !         psic(nlsm(igk_k(ig,1)))=conjg(orbital(ig,ibnd))
    !         !
    !      enddo
    !   endif
    !   !
    !   call cft3s(psic,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,2)
    CALL stop_clock( 'fft_orbital' )

    END SUBROUTINE fft_orbital_gamma
    ! 
    !
    !--------------------------------------------------------------------------
    SUBROUTINE bfft_orbital_gamma (orbital, ibnd, nbnd,conserved)
    !--------------------------------------------------------------------------
    !
    ! OBM 241008
    ! This driver subroutine -back- transforms the given orbital using fft using the already existent data
    ! in psic. Warning! This subroutine does not reset the orbital, use carefully! 
    ! Warning 2! In order to be fast, no checks on the supplied data are performed!
    ! Variables:
    ! orbital: the orbital to be transformed
    ! ibnd: band index
    ! nbnd: total number of bands
    use wavefunctions_module,     only : psic
    use gsmooth,                  only : nr1s,nr2s,nr3s,nrx1s,nrx2s,&
       nrx3s,nrxxs,nls,nlsm,doublegrid  
    USE kinds,         ONLY : DP
    USE fft_parallel,  ONLY : tg_cft3s
    USE fft_base,      ONLY : dffts
    USE control_flags, ONLY : use_task_groups
    USE task_groups,   ONLY : tg_gather
    USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist

    implicit none
     
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
    complex(DP),intent(out) :: orbital(:,:)
   
    !integer :: ig
    
    logical, optional :: conserved !if this flag is true, the orbital is stored in temporary memory
    
    !Internal temporary variables
    complex(DP) :: fp, fm
    integer :: i,  j, incr, ierr, idx, ioff, nsiz
    logical :: use_tg

    !Task groups
    INTEGER :: recv_cnt( nogrp ), recv_displ( nogrp )
    INTEGER :: v_siz
    !print *, "->fourier space"
    CALL start_clock( 'bfft_orbital' )
    !New task_groups versions
    use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp )
    IF( use_tg ) THEN
      call tg_cft3s ( tg_psic, dffts, -2, use_tg )
              !
        ioff   = 0
        !
        DO idx = 1, 2*nogrp, 2
           !
           IF( idx + ibnd - 1 < nbnd ) THEN
              DO j = 1, npw_k(1)
                 fp= ( tg_psic( nls(igk_k(j,1)) + ioff ) +  tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk_k(j,1)) + ioff ) -  tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 orbital (j, ibnd+idx-1) =  CMPLX( DBLE(fp), AIMAG(fm))
                 orbital (j, ibnd+idx  ) =  CMPLX(AIMAG(fp),- DBLE(fm))
              END DO
           ELSE IF( idx + ibnd - 1 == nbnd ) THEN
              DO j = 1, npw_k(1)
                 orbital (j, ibnd+idx-1) =  tg_psic( nls(igk_k(j,1)) + ioff )
              END DO
           END IF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
           !
        END DO
        !
        if (present(conserved)) then
         if (conserved .eqv. .true.) then
          if (allocated(tg_psic_temp)) deALLOCATE( tg_psic_temp )
         endif
        endif

    ELSE !Non task_groups version
          !larger memory slightly faster
          call cft3s(psic,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
           

          if (ibnd < nbnd) then

           ! two ffts at the same time
           do j = 1, npw_k(1)
              fp = (psic (nls(igk_k(j,1))) + psic (nlsm(igk_k(j,1))))*0.5d0
              fm = (psic (nls(igk_k(j,1))) - psic (nlsm(igk_k(j,1))))*0.5d0
              orbital( j, ibnd)   = CMPLX( DBLE(fp), AIMAG(fm))
              orbital( j, ibnd+1) = CMPLX(AIMAG(fp),- DBLE(fm))
           enddo
        else
           do j = 1, npw_k(1)
              orbital(j, ibnd)   =  psic (nls(igk_k(j,1)))
           enddo
        end if
        if (present(conserved)) then
         if (conserved .eqv. .true.) then
           if (allocated(psic_temp) ) deallocate(psic_temp) 
         endif
        endif
    ENDIF
    !! OLD VERSION Based on the algorithm found in lr_apply_liovillian
    !!print * ,"a"
    !call cft3s(psic,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
    !!
    !!print *, "b"
    !if (ibnd<nbnd) then
    !   !
    !   do ig=1,npw_k(1)
    !      !
    !      fp=(psic(nls(igk_k(ig,1)))&
    !           +psic(nlsm(igk_k(ig,1))))*(0.50d0,0.0d0)
    !      !
    !      fm=(psic(nls(igk_k(ig,1)))&
    !           -psic(nlsm(igk_k(ig,1))))*(0.50d0,0.0d0)
    !      !
    !      orbital(ig,ibnd)=cmplx(dble(fp),aimag(fm),dp)
    !      !
    !      orbital(ig,ibnd+1)=cmplx(aimag(fp),-dble(fm),dp)
    !      !
    !   enddo
    !   !
    !else
    !   !
    !   do ig=1,npw_k(1)
    !      !
    !      orbital(ig,ibnd)=psic(nls(igk_k(ig,1)))
    !      !
    !   enddo
    !   !
    !endif
    !print * , "c"
    !
    !
    CALL stop_clock( 'bfft_orbital' )
    
    END SUBROUTINE bfft_orbital_gamma
    !
  !--------------------------------------------------------------------------
  SUBROUTINE fft_orbital_k (orbital, ibnd, nbnd,conserved)
   !--------------------------------------------------------------------------
    !
    ! OBM 110908
    ! This subroutine transforms the given orbital using fft and puts the result in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    ! orbital: the orbital to be transformed
    ! ibnd: band index
    ! nbnd: total number of bands
    use wavefunctions_module,     only : psic
    use gsmooth,                  only : nr1s,nr2s,nr3s,nrx1s,nrx2s,&
       nrx3s,nrxxs,nls,nlsm,doublegrid  
    USE kinds,         ONLY : DP
    USE fft_parallel,  ONLY : tg_cft3s
    USE fft_base,      ONLY : dffts
    USE control_flags, ONLY : use_task_groups
    USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
    USE wvfct,         ONLY : igk
    
    implicit none
    
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
    complex(DP),intent(in) :: orbital(:,:)
    logical, optional :: conserved !if this flag is true, the orbital is stored in temporary memory          
    
    ! Internal variables
    integer :: j, ioff, idx
    logical :: use_tg
          
    CALL start_clock( 'fft_orbital' )    
    use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp ) 
    
    IF( use_tg ) THEN
             !
             tg_psic = ( 0.D0, 0.D0 )
             ioff   = 0
             !
             DO idx = 1, nogrp
                !
                IF( idx + ibnd - 1 <= nbnd ) THEN
                   !DO j = 1, size(orbital,1)
                      tg_psic( nls( igk(:) ) + ioff ) = orbital(:,idx+ibnd-1)
                   !END DO
                END IF

                ioff = ioff + dffts%nnrx

             END DO
             !
             call tg_cft3s ( tg_psic, dffts, 2, use_tg )
             if (present(conserved)) then
              if (conserved .eqv. .true.) then
               if (.not. allocated(tg_psic_temp)) ALLOCATE( tg_psic_temp( dffts%nnrx * nogrp ) )
               tg_psic_temp=tg_psic
              endif
             endif
             !
    ELSE  !non task_groups version
             !
             psic(1:nrxxs) = ( 0.D0, 0.D0 )
             !
             psic(nls(igk(:))) = orbital(:,ibnd)
             !
             CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
             if (present(conserved)) then
              if (conserved .eqv. .true.) then
               if (.not. allocated(psic_temp) ) allocate (psic_temp(size(psic))) 
               psic_temp=psic
              endif
             endif
             !
    ENDIF
    CALL stop_clock( 'fft_orbital' )
    END SUBROUTINE fft_orbital_k
   !--------------------------------------------------------------------------
   SUBROUTINE bfft_orbital_k (orbital, ibnd, nbnd,conserved)
   !--------------------------------------------------------------------------
    !
    ! OBM 110908
    ! This subroutine transforms the given orbital using fft and puts the result in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    ! orbital: the orbital to be transformed
    ! ibnd: band index
    ! nbnd: total number of bands
    use wavefunctions_module,     only : psic
    use gsmooth,                  only : nr1s,nr2s,nr3s,nrx1s,nrx2s,&
       nrx3s,nrxxs,nls,nlsm,doublegrid  
    USE kinds,         ONLY : DP
    USE fft_parallel,  ONLY : tg_cft3s
    USE fft_base,      ONLY : dffts
    USE control_flags, ONLY : use_task_groups
    USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
    USE wvfct,         ONLY : igk
    
    implicit none
    
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
    complex(DP),intent(out) :: orbital(:,:)
    logical, optional :: conserved !if this flag is true, the orbital is stored in temporary memory          
    
    ! Internal variables
    integer :: j, ioff, idx
    logical :: use_tg
          
   CALL start_clock( 'bfft_orbital' )    
    use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp ) 
    
    IF( use_tg ) THEN
             !
             call tg_cft3s ( tg_psic, dffts, -2, use_tg )
             !
             ioff   = 0 
             !  
             DO idx = 1, nogrp
                !  
                IF( idx + ibnd - 1 <= nbnd ) THEN
                   orbital (:, ibnd+idx-1) = tg_psic( nls(igk(:)) + ioff )
                END IF 
                !
                ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
                !
             END DO 
            if (present(conserved)) then
             if (conserved .eqv. .true.) then
              if (allocated(tg_psic_temp)) deALLOCATE( tg_psic_temp )
             endif
            endif
         !
    ELSE !non task groups version
             !
             CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
             !
             orbital(:,ibnd) = psic(nls(igk(:)))
             !
        if (present(conserved)) then
         if (conserved .eqv. .true.) then
           if (allocated(psic_temp) ) deallocate(psic_temp) 
         endif
        endif
    ENDIF
    CALL stop_clock( 'bfft_orbital' )
    END SUBROUTINE bfft_orbital_k
  !--------------------------------------------------------------------------
  SUBROUTINE v_loc_psir (ibnd, nbnd)
    !--------------------------------------------------------------------------
    ! Basically the same thing as v_loc but without implicit fft
    ! modified for real space implementation
    ! OBM 241008
    !
    use wavefunctions_module,     only : psic
    use gsmooth,                  only : nr1s,nr2s,nr3s,nrx1s,nrx2s,&
       nrx3s,nrxxs,nls,nlsm,doublegrid  
    USE kinds,         ONLY : DP
    USE fft_parallel,  ONLY : tg_cft3s
    USE fft_base,      ONLY : dffts
    USE control_flags, ONLY : use_task_groups
    USE task_groups,   ONLY : tg_gather
    USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
    USE scf,           ONLY : vrs 
    USE lsda_mod,      ONLY : current_spin
    

    implicit none
     
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
        
    !Internal temporary variables
    complex(DP) :: fp, fm
    integer :: i,  j, incr, ierr, idx, ioff, nsiz
    logical :: use_tg

    !Task groups
    REAL(DP),    ALLOCATABLE :: tg_v(:)
    INTEGER :: recv_cnt( nogrp ), recv_displ( nogrp )
    INTEGER :: v_siz
    CALL start_clock( 'v_loc_psir' )
        use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp )
     IF( use_tg ) THEN
        if (ibnd == 1 ) then 
          CALL tg_gather( dffts, vrs(:,current_spin), tg_v ) !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        endif
        !
        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) + tg_psic_temp (j) * tg_v(j)
        enddo
        !
        DEALLOCATE( tg_v )
     ELSE
        !   product with the potential v on the smooth grid
        !
        do j = 1, nrxxs
           psic (j) = psic (j) + psic_temp (j) * vrs(j,current_spin)
        enddo
     END IF
  CALL stop_clock( 'v_loc_psir' )
  END SUBROUTINE v_loc_psir
    !--------------------------------------------------------------------------
    
    
    
    
    
    
    
    !ERASE THE SUBROUTINES AFTER THIS LINE 
    SUBROUTINE check_fft_orbital_gamma(orbital, ibnd, nbnd)
     USE gvect, only : nrxx
     use gsmooth, only : nrxxs
     use wavefunctions_module,     only : psic
    implicit none
    integer, intent(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
                          
    complex(DP) :: orbital(:,:)
    
    complex(DP),allocatable :: compare(:)
    
    
    integer :: cnt
    if (size(orbital,1) .eq. nrxx) then 
     print *,"Fine grid"
    else if (size(orbital,1) .eq. nrxxs) then
     print *,"Smooth grid"
    else 
     print *,"UNKNOWN GRID" 
    endif
    !print *, orbital
    allocate(compare(size(orbital,1)))
    print *, "Check started"
    compare(1:size(orbital,1))=orbital(:,ibnd) 
    print *, "Step 1 complete (array copying)"
    call fft_orbital_gamma(orbital,ibnd,nbnd,.true.)
    print *, "Step 2 complete (fft)"
    do cnt=1, size(psic_temp,1)
     if( ABS(psic(cnt) - psic_temp(cnt)) > 1e-6 ) print *, "fft check failed"
    enddo
    orbital=0
    call bfft_orbital_gamma(orbital,ibnd,nbnd)
    print *, "Step 3 complete (bfft)"
    do cnt=1, size(orbital,1)
     if( ABS(orbital(cnt,ibnd) - compare(cnt)) > 1e-6 ) print *, "fft check failed"
    enddo
    END SUBROUTINE check_fft_orbital_gamma


    !
END MODULE realus
