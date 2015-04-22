!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  ! ... modified by P. Umari and G. Stenuit (2009)
  ! ... GWW-specific stuff moved out by P. Giannozzi (2015)
  !
  REAL(DP), ALLOCATABLE :: boxrad(:) ! radius of boxes, does not depend on the grid
  ! Beta function in real space
  INTEGER,  ALLOCATABLE :: box_beta(:,:), maxbox_beta(:)
  REAL(DP), ALLOCATABLE :: betasave(:,:,:)
  REAL(DP), ALLOCATABLE :: boxrad_beta(:)
  REAL(DP), ALLOCATABLE :: boxdist_beta(:,:), xyz_beta(:,:,:)
  REAL(DP), ALLOCATABLE :: spher_beta(:,:,:)
  !General
  LOGICAL               :: real_space
  INTEGER :: real_space_debug = 0 ! FIXME: must disappear
  ! if true perform calculations in real spave
  INTEGER               :: initialisation_level
  ! init_realspace_vars sets this to 3; qpointlist adds 5; betapointlist adds 7
  ! so the value should be 15 if the real space routine is initalised properly

  INTEGER, ALLOCATABLE :: &
       igk_k(:,:),&       ! The g<->k correspondance for each k point
       npw_k(:)           ! number of plane waves at each k point
  ! They are (used many times, it is much better to hold them in memory
  ! FIXME: npw_k is redundant. it is already there
  !
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  COMPLEX(DP), ALLOCATABLE :: psic_temp(:),tg_psic_temp(:) !Copies of psic and tg_psic
  COMPLEX(DP), ALLOCATABLE :: tg_vrs(:) !task groups linear V memory
  COMPLEX(DP), ALLOCATABLE :: psic_box_temp(:),tg_psic_box_temp(:)
  !
  ! Contains the augmentation functions and related quantities, in realspace for a single atom
  TYPE realsp_augmentation
     REAL(DP),ALLOCATABLE :: qr(:,:) 
        ! (maxbox,number_of_q_funcs) the Q functions sampled over R points
     INTEGER,ALLOCATABLE  :: box(:) 
        ! (maxbox) Index of point R in the global order of the R-space grid
     INTEGER :: maxbox = 0
        !  number of R points in the augmentation sphere of this atom
  END TYPE realsp_augmentation
  ! Augmentation functions on the RHO (HARD) grid for all atoms
  TYPE(realsp_augmentation),POINTER :: tabp(:) => null()
  ! Augmentation functions on the SMOOTH grid for all atoms
  TYPE(realsp_augmentation),POINTER :: tabs(:) => null()
  !
  CONTAINS
    !------------------------------------------------------------------------
    SUBROUTINE generate_qpointlist
      !------------------------------------------------------------------------
      USE fft_base,     ONLY : dfftp, dffts
      USE funct,        ONLY : dft_is_hybrid
      USE gvecs,        ONLY : doublegrid
      USE io_global,    ONLY : stdout
      IMPLICIT NONE
      !
      ! 1. initialize hard grid
      WRITE(stdout, '(/,5x,a)') "Initializing real-space augmentation for DENSE grid"
      CALL qpointlist(dfftp, tabp) 
      !
      ! 2. initialize smooth grid (only for EXX at this moment)
      IF ( dft_is_hybrid() ) THEN
        IF(doublegrid)THEN
          WRITE(stdout, '(5x,a)') "Initializing real-space augmentation for SMOOTH grid"
          CALL qpointlist(dffts, tabs)
        ELSE
          ! smooth and rho grid are the same if not double grid
          WRITE(stdout, '(7x,a)') " SMOOTH grid -> DENSE grid"
          tabs => tabp
        ENDIF
      ENDIF
      !
      RETURN
      !------------------------------------------------------------------------
    END SUBROUTINE generate_qpointlist
    !------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE init_realspace_vars()
    !---------------------------------------------------------------------------
    !This subroutine should be called to allocate/reset real space related variables.
    !---------------------------------------------------------------------------
     USE wvfct,                ONLY : npwx,npw, igk, g2kin, ecutwfc
     USE klist,                ONLY : nks, xk
     USE gvect,                ONLY : ngm, g
     USE cell_base,            ONLY : tpiba2
     USE control_flags,        ONLY : tqr
     USE fft_base,             ONLY : dffts
     USE io_global,            ONLY : stdout


     IMPLICIT NONE

     INTEGER :: ik

     !print *, "<<<<<init_realspace_vars>>>>>>>"

     IF ( allocated( igk_k ) )     DEALLOCATE( igk_k )
     IF ( allocated( npw_k ) )     DEALLOCATE( npw_k )

     ALLOCATE(igk_k(npwx,nks))
     ALLOCATE(npw_k(nks))
     !real space, allocation for task group fft work arrays
     IF( dffts%have_task_groups ) THEN
        !
        IF (allocated( tg_psic ) ) DEALLOCATE( tg_psic )
        !
        ALLOCATE( tg_psic( dffts%tg_nnr * dffts%nogrp ) )
        ALLOCATE( tg_vrs( dffts%tg_nnr * dffts%nogrp ) )
        !
     ENDIF
     !
     DO ik=1,nks
      !
      CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
      npw_k(ik) = npw
      igk_k(:,ik) = igk(:)
      !
     ENDDO

     initialisation_level = initialisation_level + 7
     IF (real_space_debug > 20 .and. real_space_debug < 30) THEN
       real_space=.false.
       IF (tqr) THEN
         tqr = .false.
         WRITE(stdout,'("Debug level forced tqr to be set false")')
       ELSE
         WRITE(stdout,'("tqr was already set false")')
       ENDIF
       real_space_debug=real_space_debug-20
     ENDIF

    END SUBROUTINE init_realspace_vars
    !------------------------------------------------------------------------
    SUBROUTINE deallocatenewdreal()
      !------------------------------------------------------------------------
      !
      USE ions_base,  ONLY : nat
      USE gvecs,      ONLY : doublegrid
      IMPLICIT NONE
      INTEGER :: ia
      IF ( allocated( boxrad ) )  DEALLOCATE( boxrad )
      ! ------------ tabp ------------ 
      IF ( associated( tabp) ) THEN
        DO ia=1,nat
          IF(allocated(tabp(ia)%qr))  DEALLOCATE(tabp(ia)%qr)
          IF(allocated(tabp(ia)%box)) DEALLOCATE(tabp(ia)%box)
          tabp(ia)%maxbox = 0
        ENDDO
        DEALLOCATE(tabp)
      ENDIF
      ! ------------ tabs ------------ 
      IF(doublegrid)THEN
        IF ( associated( tabs) ) THEN
          DO ia=1,nat
            IF(allocated(tabs(ia)%qr))  DEALLOCATE(tabs(ia)%qr)
            IF(allocated(tabs(ia)%box)) DEALLOCATE(tabs(ia)%box)
            tabs(ia)%maxbox = 0
          ENDDO
          DEALLOCATE(tabs)
        ENDIF
      ELSE
        NULLIFY(tabs)
      ENDIF

      !
    END SUBROUTINE deallocatenewdreal
    !
    !------------------------------------------------------------------------
    SUBROUTINE qpointlist(dfft, tabp)
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
      ! ... In tabp the q value interpolated in these boxes.
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi, fpi, eps8, eps16
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE uspp,       ONLY : okvan, indv, nhtol, nhtolm, ap, nhtoj, lpx, lpl
      USE uspp_param, ONLY : upf, lmaxq, nh, nhm
      USE atom,       ONLY : rgrid
      USE fft_types,  ONLY : fft_dlay_descriptor
      USE mp_bands,   ONLY : me_bgrp
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      TYPE(fft_dlay_descriptor),INTENT(in)    :: dfft
      TYPE(realsp_augmentation),POINTER,INTENT(inout) :: tabp(:)
      !
      INTEGER               :: qsdim, ia, mbia, iqs, iqsia, nfuncs
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
      REAL(DP), ALLOCATABLE :: boxdist(:,:), xyz(:,:,:)
      REAL(DP), ALLOCATABLE :: spher(:,:,:)
      !
      initialisation_level = 3
      IF ( .not. okvan ) RETURN
      !
      CALL start_clock( 'realus' )
      !
      ! ... tabp is deallocated here to free the memory for the buffers
      !
      IF( associated(tabp) ) THEN
        DO ia=1,nat
          IF(allocated(tabp(ia)%qr))  DEALLOCATE(tabp(ia)%qr)
          IF(allocated(tabp(ia)%box)) DEALLOCATE(tabp(ia)%box)
        ENDDO
        DEALLOCATE(tabp)
      ENDIF
      !
      ALLOCATE(tabp(nat))

      !
      IF ( .not. allocated( boxrad ) ) THEN
         !
         ! ... here we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         !
         ALLOCATE( boxrad( nsp ) )
         !
         boxrad(:) = 0.D0
         !
         DO nt = 1, nsp
            IF ( .not. upf(nt)%tvanp ) CYCLE
            DO ijv = 1, upf(nt)%nbeta*(upf(nt)%nbeta+1)/2
               DO indm = upf(nt)%mesh,1,-1
                  !
                  IF( upf(nt)%q_with_l ) THEN
                     aux = sum(abs( upf(nt)%qfuncl(indm,ijv,:) ))
                  ELSE
                     aux = abs( upf(nt)%qfunc(indm,ijv) )
                  ENDIF
                  IF ( aux > eps16 ) THEN
                     boxrad(nt) = max( rgrid(nt)%r(indm), boxrad(nt) )
                     exit
                  ENDIF
                  !
               ENDDO
            ENDDO
         ENDDO
         !
         boxrad(:) = boxrad(:) / alat
         !
      ENDIF
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is provided here
      !
      mbr = maxval( boxrad(:) )
      !
      mbx = mbr*sqrt( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
      mby = mbr*sqrt( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
      mbz = mbr*sqrt( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
      !
      dmbx = 2*anint( mbx*dfft%nr1x ) + 2
      dmby = 2*anint( mby*dfft%nr2x ) + 2
      dmbz = 2*anint( mbz*dfft%nr3x ) + 2
      !
      roughestimate = anint( dble( dmbx*dmby*dmbz ) * pi / 6.D0 )
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
      ! ... now we find the points
      !
#if defined (__MPI)
      idx0 = dfft%nr1x*dfft%nr2x * sum ( dfft%npp(1:me_bgrp) )
#else
      idx0 = 0
#endif
      !
      inv_nr1 = 1.D0 / dble( dfft%nr1 )
      inv_nr2 = 1.D0 / dble( dfft%nr2 )
      inv_nr3 = 1.D0 / dble( dfft%nr3 )
      !
      DO ia = 1, nat
         !
         nt = ityp(ia)
         !
         IF ( .not. upf(nt)%tvanp ) CYCLE
         !
         boxradsq_ia = boxrad(nt)**2
         !
         tau_ia(1) = tau(1,ia)
         tau_ia(2) = tau(2,ia)
         tau_ia(3) = tau(3,ia)
         !
         DO ir = 1, dfft%nnr
            !
            ! ... three dimensional indices (i,j,k)
            !
            idx   = idx0 + ir - 1
            k     = idx / (dfft%nr1x*dfft%nr2x)
            idx   = idx - (dfft%nr1x*dfft%nr2x)*k
            j     = idx / dfft%nr1x
            idx   = idx - dfft%nr1x*j
            i     = idx
            !
            ! ... do not include points outside the physical range!
            !
            IF ( i >= dfft%nr1 .or. j >= dfft%nr2 .or. k >= dfft%nr3 ) CYCLE
            !
            DO ipol = 1, 3
               posi(ipol) = dble( i )*inv_nr1*at(ipol,1) + &
                            dble( j )*inv_nr2*at(ipol,2) + &
                            dble( k )*inv_nr3*at(ipol,3)
            ENDDO
            !
            posi(:) = posi(:) - tau_ia(:)
            !
            ! ... minimum image convention
            !
            CALL cryst_to_cart( 1, posi, bg, -1 )
            posi(:) = posi(:) - anint( posi(:) )
            CALL cryst_to_cart( 1, posi, at, 1 )
            distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
            !
            IF ( distsq < boxradsq_ia ) THEN
               !
               mbia = tabp(ia)%maxbox + 1
               IF(mbia>roughestimate) THEN
                  CALL errore( 'qpointlist', 'rough-estimate is too rough', 3 )
               ENDIF
               tabp(ia)%maxbox     = mbia
               buffpoints(mbia,ia) = ir
               buffdist(mbia,ia)   = sqrt( distsq )*alat
               xyz(:,mbia,ia)      = posi(:)*alat
               !
            ENDIF
         ENDDO
      ENDDO
      !
      ! ... count the points and store them in a more convenient place
      !
      goodestimate = 0
      DO ia = 1,nat
        goodestimate = max( goodestimate, tabp(ia)%maxbox )
        !
        !IF ( allocated( tabp(ia)%box ) ) DEALLOCATE( tabp(ia)%box )
        ALLOCATE( tabp(ia)%box(tabp(ia)%maxbox) )
        tabp(ia)%box(:) = buffpoints(1:tabp(ia)%maxbox,ia)
      ENDDO
      !
      IF ( goodestimate > roughestimate ) &
         CALL errore( 'qpointlist', 'rough-estimate is too rough', 2 )
      !
      ALLOCATE( boxdist( goodestimate, nat ) )
      !
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
      IF ( allocated( spher ) ) DEALLOCATE( spher )
      !
      ALLOCATE( spher( goodestimate, lamx2, nat ) )
      !
      spher(:,:,:) = 0.D0
      !
      DO ia = 1, nat
         !
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
         idimension = tabp(ia)%maxbox
         ALLOCATE( rl( 3, idimension ), rl2( idimension ) )

         DO ir = 1, idimension
            rl(:,ir) = xyz(:,ir,ia)
            rl2(ir) = rl(1,ir)**2 + rl(2,ir)**2 + rl(3,ir)**2
         ENDDO
         !
         ALLOCATE( tempspher( idimension, lamx2 ) )
         CALL ylmr2( lamx2, idimension, rl, rl2, tempspher )
         spher(1:idimension,:,ia) = tempspher(:,:)
         DEALLOCATE( rl, rl2, tempspher )
         !
      ENDDO
      !
      DEALLOCATE( xyz )
      !
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:tabp' )
      !
      ! ... let's do the main work
      !
      qsdim = 0
      DO ia = 1, nat
         mbia = tabp(ia)%maxbox
         IF ( mbia == 0 ) CYCLE
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
               qsdim = qsdim + mbia
            ENDDO
         ENDDO
      ENDDO
      !
      ! ... the source is inspired by init_us_1
      !
      ! ... we perform two steps: first we compute for each l the qtot
      ! ... (radial q), then we interpolate it in our mesh, and then we
      ! ... add it to tabp with the correct spherical harmonics
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
         mbia = tabp(ia)%maxbox
         IF ( mbia == 0 ) CYCLE
         !
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
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
         nfuncs = (  nh(nt)**2+nh(nt)  )/2
         ALLOCATE(tabp(ia)%qr(mbia, nfuncs))
         tabp(ia)%qr=0._dp
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
                  IF ( .not. ( l >= abs( lllnbnt - lllmbnt ) .and. &
                               l <= lllnbnt + lllmbnt        .and. &
                               mod( l + lllnbnt + lllmbnt, 2 ) == 0 ) ) CYCLE
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
                        ENDIF
                      ENDDO
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
                  DO ir = 1, tabp(ia)%maxbox
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
                     ENDIF
                     !
                     ijh = 0
                     !
                     DO ih = 1, nh(nt)
                        DO jh = ih, nh(nt)
                           !
                           iqs = iqsia + ijh*mbia + ir
                           ijh = ijh + 1
                           !
                           IF ( .not.( nb == indv(ih,nt) .and. &
                                       mb == indv(jh,nt) ) ) CYCLE
                           !
                           DO lm = l**2+1, (l+1)**2
                              tabp(ia)%qr(ir,ijh) = tabp(ia)%qr(ir,ijh) + &
                                           qtot_int*spher(ir,lm,ia)*&
                                           ap(lm,nhtolm(ih,nt),nhtolm(jh,nt))
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !
         iqsia = iqs
         !
         DEALLOCATE( qtot )
         DEALLOCATE( xsp )
         DEALLOCATE( ysp )
         DEALLOCATE( wsp )
         !
      ENDDO
      !
      DEALLOCATE( boxdist )
      DEALLOCATE( spher )
      !
      CALL stop_clock( 'realus:tabp' )
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
      ! ... In tabp the q value interpolated in these boxes.
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      ! The source inspired by qsave
      !
      USE constants,  ONLY : pi, eps8, eps16
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE uspp,       ONLY : okvan, indv, nhtol, nhtolm, ap
      USE uspp_param, ONLY : upf, lmaxq, nh, nhm
      USE atom,       ONLY : rgrid
      USE fft_base,   ONLY : dffts
      USE mp_bands,   ONLY : me_bgrp
      USE splinelib,  ONLY : spline, splint
      USE ions_base,  ONLY : ntyp => nsp
      !
      IMPLICIT NONE
      !
      INTEGER               :: betasdim, ia, it, mbia
      INTEGER               :: indm, inbrx, idimension, ih
      INTEGER               :: roughestimate, goodestimate, lamx2, nt
      INTEGER,  ALLOCATABLE :: buffpoints(:,:)
      REAL(DP), ALLOCATABLE :: buffdist(:,:)
      REAL(DP)              :: distsq, qtot_int, first, second
      INTEGER               :: idx0, idx, ir
      INTEGER               :: i, j, k, ipol, lm, nb
      REAL(DP)              :: posi(3)
      REAL(DP), ALLOCATABLE :: rl(:,:), rl2(:)
      REAL(DP), ALLOCATABLE :: tempspher(:,:), qtot(:,:,:), &
                               xsp(:), ysp(:), wsp(:), d1y(:), d2y(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
      REAL(DP)              :: inv_nr1s, inv_nr2s, inv_nr3s, tau_ia(3), boxradsq_ia
      !
      initialisation_level = initialisation_level + 5
      IF ( .not. okvan ) CALL errore &
                        ('betapointlist','real space routines for USPP only',1)
      !
      !print *, "<<<betapointlist>>>"
      !
      CALL start_clock( 'betapointlist' )
      !
      ! ... betasave is deallocated here to free the memory for the buffers
      !
      IF ( allocated( betasave ) ) DEALLOCATE( betasave )
      !
      IF ( .not. allocated( boxrad_beta ) ) THEN
         !
         ! ... here we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         !
         ALLOCATE( boxrad_beta( nsp ) )
         boxrad_beta(:) = 0.D0
         !
         DO it = 1, nsp
            DO inbrx = 1, upf(it)%nbeta
               DO indm = upf(it)%kkbeta, 1, -1
                  IF ( abs( upf(it)%beta(indm,inbrx) ) > 0.d0 ) THEN
                     boxrad_beta(it) = max( rgrid(it)%r(indm), boxrad_beta(it) )
                     CYCLE
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         !
         boxrad_beta(:) = boxrad_beta(:) / alat
         !
      ENDIF
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is provided here
      !
      mbr = maxval( boxrad_beta(:) )
      !
      mbx = mbr*sqrt( bg(1,1)**2 + bg(1,2)**2 + bg(1,3)**2 )
      mby = mbr*sqrt( bg(2,1)**2 + bg(2,2)**2 + bg(2,3)**2 )
      mbz = mbr*sqrt( bg(3,1)**2 + bg(3,2)**2 + bg(3,3)**2 )
      !
      dmbx = 2*anint( mbx*dffts%nr1x ) + 2
      dmby = 2*anint( mby*dffts%nr2x ) + 2
      dmbz = 2*anint( mbz*dffts%nr3x ) + 2
      !
      roughestimate = anint( dble( dmbx*dmby*dmbz ) * pi / 6.D0 )
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
      IF ( .not.allocated( maxbox_beta ) ) ALLOCATE( maxbox_beta( nat ) )
      !
      maxbox_beta(:) = 0
      !
      ! ... now we find the points
      !
      ! The beta functions are treated on smooth grid
#if defined (__MPI)
      idx0 = dffts%nr1x*dffts%nr2x * dffts%ipp(me_bgrp+1)
#else
      idx0 = 0
#endif
      !
      inv_nr1s = 1.D0 / dble( dffts%nr1 )
      inv_nr2s = 1.D0 / dble( dffts%nr2 )
      inv_nr3s = 1.D0 / dble( dffts%nr3 )
      !
      DO ia = 1, nat
         !
         IF ( .not. upf(ityp(ia))%tvanp ) CYCLE
         !
         boxradsq_ia = boxrad_beta(ityp(ia))**2
         !
         tau_ia(1) = tau(1,ia)
         tau_ia(2) = tau(2,ia)
         tau_ia(3) = tau(3,ia)
         !
         DO ir = 1, dffts%nr1x*dffts%nr2x * dffts%npl
            !
            ! ... three dimensional indexes
            !
            idx = idx0 + ir - 1
            k   = idx / (dffts%nr1x*dffts%nr2x)
            idx = idx - (dffts%nr1x*dffts%nr2x)*k
            j   = idx / dffts%nr1x
            idx = idx - dffts%nr1x*j
            i   = idx
            !
            DO ipol = 1, 3
               posi(ipol) = dble( i )*inv_nr1s*at(ipol,1) + &
                            dble( j )*inv_nr2s*at(ipol,2) + &
                            dble( k )*inv_nr3s*at(ipol,3)
            ENDDO
            !
            posi(:) = posi(:) - tau_ia(:)
            !
            ! ... minimum image convenction
            !
            CALL cryst_to_cart( 1, posi, bg, -1 )
            !
            posi(:) = posi(:) - anint( posi(:) )
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
               buffdist(mbia,ia)   = sqrt( distsq )*alat
               xyz_beta(:,mbia,ia) = posi(:)*alat
               !
            ENDIF
         ENDDO
      ENDDO
      !
      goodestimate = maxval( maxbox_beta )
      !
      IF ( goodestimate > roughestimate ) &
         CALL errore( 'betapointlist', 'rough-estimate is too rough', 2 )
      !
      ! ... now store them in a more convenient place
      !
      IF ( allocated( box_beta ) )     DEALLOCATE( box_beta )
      IF ( allocated( boxdist_beta ) ) DEALLOCATE( boxdist_beta )
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
      IF ( allocated( spher_beta ) ) DEALLOCATE( spher_beta )
      !
      ALLOCATE( spher_beta( goodestimate, lamx2, nat ) )
      !
      spher_beta(:,:,:) = 0.D0
      !
      DO ia = 1, nat
         !
         IF ( .not. upf(ityp(ia))%tvanp ) CYCLE
         !
         idimension = maxbox_beta(ia)
         ALLOCATE( rl( 3, idimension ), rl2( idimension ) )
         !
         DO ir = 1, idimension
            rl(:,ir) = xyz_beta(:,ir,ia)
            rl2(ir) = rl(1,ir)**2 + rl(2,ir)**2 + rl(3,ir)**2
         ENDDO
         !
         ALLOCATE( tempspher( idimension, lamx2 ) )
         CALL ylmr2( lamx2, idimension, rl, rl2, tempspher )
         spher_beta(1:idimension,:,ia) = tempspher(:,:)
         DEALLOCATE( rl, rl2, tempspher )
         !
      ENDDO
      !
      DEALLOCATE( xyz_beta )
      !
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:tabp' )
      !
      ! ... let's do the main work
      !
      betasdim = 0
      DO ia = 1, nat
         mbia = maxbox_beta(ia)
         IF ( mbia == 0 ) CYCLE
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
         DO ih = 1, nh(nt)
            betasdim = betasdim + mbia
         ENDDO
      ENDDO
      !
      ALLOCATE( betasave( nat, nhm, goodestimate )  )
      !
      betasave = 0.D0
      ! Box is set, Y_lm is known in the box, now the calculation can commence
      ! Reminder: In real space
      ! |Beta_lm(r)>=f_l(r).Y_lm(r)
      ! In q space (calculated in init_us_1 and then init_us_2 )
      ! |Beta_lm(q)>= (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
      ! Where
      ! f_l(q)=\int_0 ^\infty dr r^2 f_l (r) j_l(q.r)
      !
      ! We know f_l(r) and Y_lm(r) for certain points,
      ! basically we interpolate the known values to new mesh using splint
      !      iqs   = 0
      !
      DO ia = 1, nat
         !
         mbia = maxbox_beta(ia)
         IF ( mbia == 0 ) CYCLE
         !
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
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
            !OBM rgrid(nt)%r(1) == 0, attempting correction
            ! In the UPF file format, beta field is r*|beta>
            IF (rgrid(nt)%r(1)==0) THEN
             ysp(2:) = upf(nt)%beta(2:upf(nt)%kkbeta,nb) / rgrid(nt)%r(2:upf(nt)%kkbeta)
             ysp(1)=0.d0
            ELSE
             ysp(:) = upf(nt)%beta(1:upf(nt)%kkbeta,nb) / rgrid(nt)%r(1:upf(nt)%kkbeta)
            ENDIF

            ALLOCATE( d1y(upf(nt)%kkbeta), d2y(upf(nt)%kkbeta) )
            CALL radial_gradient(ysp(1:upf(nt)%kkbeta), d1y, &
                                 rgrid(nt)%r, upf(nt)%kkbeta, 1)
            CALL radial_gradient(d1y, d2y, rgrid(nt)%r, upf(nt)%kkbeta, 1)

            first = d1y(1) ! first derivative in first point
            second =d2y(1) ! second derivative in first point
            DEALLOCATE( d1y, d2y )


            CALL spline( xsp, ysp, first, second, wsp )


            DO ir = 1, mbia
               !
               ! ... spline interpolation
               !
               qtot_int = splint( xsp, ysp, wsp, boxdist_beta(ir,ia) ) !the value of f_l(r) in point ir in atom ia
               !
               !iqs = iqs + 1
               !
               betasave(ia,ih,ir) = qtot_int*spher_beta(ir,lm,ia) !spher_beta is the Y_lm in point ir for atom ia
               !
            ENDDO
         ENDDO
         !
         DEALLOCATE( qtot )
         DEALLOCATE( xsp )
         DEALLOCATE( ysp )
         DEALLOCATE( wsp )
         !
      ENDDO
      !
      DEALLOCATE( boxdist_beta )
      DEALLOCATE( spher_beta )
      !
      CALL stop_clock( 'realus:tabp' )
      CALL stop_clock( 'betapointlist' )
      !
    END SUBROUTINE betapointlist
    !------------------------------------------------------------------------
    SUBROUTINE newq_r(vr,deeq,skip_vltot)
    !
    !   This routine computes the integral of the perturbed potential with
    !   the Q function in real space
    !
      USE cell_base,        ONLY : omega
      USE fft_base,         ONLY : dfftp
      USE lsda_mod,         ONLY : nspin
      USE ions_base,        ONLY : nat, ityp
      USE uspp_param,       ONLY : upf, nh, nhm
      USE uspp,             ONLY : ijtoh
      USE control_flags,    ONLY : tqr
      USE noncollin_module, ONLY : nspin_mag
      USE scf,              ONLY : vltot
      USE mp_bands,         ONLY : intra_bgrp_comm
      USE mp,               ONLY : mp_sum

          IMPLICIT NONE
      !
      ! Input: potential , output: contribution to integral
      REAL(kind=dp), INTENT(in)  :: vr(dfftp%nnr,nspin)
      REAL(kind=dp), INTENT(out) :: deeq( nhm, nhm, nat, nspin )
      LOGICAL, INTENT(in) :: skip_vltot !If .false. vltot is added to vr when necessary
      !Internal
      REAL(DP), ALLOCATABLE :: aux(:)
      !
      INTEGER               :: ia, ih, jh, is, ir, nt
      INTEGER               :: mbia, nhnt, iqs
      !
      IF (tqr .and. .not. associated(tabp)) THEN
         CALL generate_qpointlist()
      ENDIF

      deeq(:,:,:,:) = 0.D0
      !
      ALLOCATE( aux( dfftp%nnr ) )
      !
      DO is = 1, nspin_mag
         !
         IF ( (nspin_mag == 4 .and. is /= 1) .or. skip_vltot ) THEN
            aux(:) = vr(:,is)
         ELSE
            aux(:) = vltot(:) + vr(:,is)
         ENDIF
         !
         iqs = 0
         !
         DO ia = 1, nat
            !
            mbia = tabp(ia)%maxbox
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            IF ( .not. upf(nt)%tvanp ) CYCLE
            !
            nhnt = nh(nt)
            !
            DO ih = 1, nhnt
               DO jh = ih, nhnt
                  DO ir = 1, mbia
                     iqs = iqs + 1
                     deeq(ih,jh,ia,is)= deeq(ih,jh,ia,is) + &
                        tabp(ia)%qr(ir,ijtoh(ih,jh,nt))*aux(tabp(ia)%box(ir))
                  ENDDO
                  deeq(jh,ih,ia,is) = deeq(ih,jh,ia,is)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      deeq(:,:,:,:) = deeq(:,:,:,:)*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      DEALLOCATE( aux )
      CALL mp_sum(  deeq(:,:,:,1:nspin_mag) , intra_bgrp_comm )
      !
    END SUBROUTINE newq_r
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
      INTEGER,  INTENT(in):: nqf, ltot, mesh
        ! input: the number of coefficients
        ! input: the angular momentum
        ! input: the number of mesh point
      REAL(DP), INTENT(in) :: r(mesh), qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(out) :: rho(mesh)
        ! output: the function to be computed
      !
      INTEGER  :: ir, i
      REAL(DP) :: rr
      !
      DO ir = 1, mesh
         rr = r(ir)**2
         rho(ir) = qfcoef(1)
         DO i = 2, nqf
            rho(ir) = rho(ir) + qfcoef(i)*rr**(i-1)
         ENDDO
         rho(ir) = rho(ir)*r(ir)**ltot
      ENDDO
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
      INTEGER,  INTENT(in):: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(in) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(out) :: rho
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      !
      rho = qfcoef(1)
      DO i = 2, nqf
         rho = rho + qfcoef(i)*rr**(i-1)
      ENDDO
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
      INTEGER,  INTENT(in) :: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(in) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(out) :: rho
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      rho = 0.D0
      DO i = max( 1, 2-ltot ), nqf
         rho = rho + qfcoef(i)*rr**(i-2+ltot)*(i-1+ltot)
      ENDDO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrptfirst
    !
    !------------------------------------------------------------------------
    SUBROUTINE setqfcorrptsecond( qfcoef, rho, r, nqf, ltot )
      !------------------------------------------------------------------------
      !
      ! ... On output it contains  Q
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(in) :: nqf, ltot
        ! input: the number of coefficients
        ! input: the angular momentum
      REAL(DP), INTENT(in) :: r, qfcoef(nqf)
        ! input: the radial mesh
        ! input: the coefficients of Q
      REAL(DP), INTENT(out) :: rho
        ! output: the function to be computed
      !
      INTEGER  :: i
      REAL(DP) :: rr
      !
      rr = r*r
      rho = 0.D0
      DO i = max( 3-ltot, 1 ), nqf
         rho = rho + qfcoef(i)*rr**(i-3+ltot)*(i-1+ltot)*(i-2+ltot)
      ENDDO
      !
      RETURN
      !
    END SUBROUTINE setqfcorrptsecond
    !
    !------------------------------------------------------------------------
    SUBROUTINE addusdens_r(rho_1,rescale)
      !------------------------------------------------------------------------
      !
      ! ... This routine adds to the charge density the part which is due to
      ! ... the US augmentation.
      !
      USE ions_base,        ONLY : nat, ityp
      USE cell_base,        ONLY : omega
      USE lsda_mod,         ONLY : nspin
      !USE scf,              ONLY : rho
      USE klist,            ONLY : nelec
      USE fft_base,         ONLY : dfftp
      USE uspp,             ONLY : okvan, becsum
      USE uspp_param,       ONLY : upf, nh
      USE noncollin_module, ONLY : noncolin, nspin_mag, nspin_lsda
      USE spin_orb,         ONLY : domag
      USE mp_pools,         ONLY : inter_pool_comm
      USE mp_bands,         ONLY : intra_bgrp_comm
      USE mp,               ONLY : mp_sum

      !
      IMPLICIT NONE
      ! The charge density to be augmented
      REAL(kind=dp), INTENT(inout) :: rho_1(dfftp%nnr,nspin_mag) 
      ! If this is the ground charge density, enable rescaling
      LOGICAL, INTENT(in) :: rescale
      !
      INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, mbia, nhnt, iqs
      CHARACTER(len=80) :: msg
      REAL(DP) :: charge
      REAL(DP) :: tolerance
      !
      !
      IF ( .not. okvan ) RETURN
      tolerance = 1.d-3
      IF ( real_space ) tolerance = 1.d-2 !Charge loss in real_space case is even worse.
      !Final verdict: Mixing of Real Space paradigm and      !Q space paradigm results in fast but not so
      ! accurate code. Not giving up though, I think
      ! I can still increase the accuracy a bit...
      !
      CALL start_clock( 'addusdens' )
      !
      DO is = 1, nspin_mag
         !
         iqs = 0
         DO ia = 1, nat
            !
            mbia = tabp(ia)%maxbox
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            IF ( .not. upf(nt)%tvanp ) CYCLE
            !
            nhnt = nh(nt)
            ijh = 0
            DO ih = 1, nhnt
               DO jh = ih, nhnt
                  ijh = ijh + 1
                  DO ir = 1, mbia
                     irb = tabp(ia)%box(ir)
                     iqs = iqs + 1
                     rho_1(irb,is) = rho_1(irb,is) + tabp(ia)%qr(ir,ijh)*becsum(ijh,ia,is)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !
      ENDDO
      !
      ! ... check the integral of the total charge
      !
      IF (rescale) THEN
      ! RHO IS NOT NECESSARILY GROUND STATE CHARGE DENSITY, thus rescaling is optional
       charge = sum( rho_1(:,1:nspin_lsda) )*omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       CALL mp_sum(  charge , intra_bgrp_comm )
       CALL mp_sum(  charge , inter_pool_comm )

       IF ( abs( charge - nelec ) / nelec > tolerance ) THEN
          !
          ! ... the error on the charge is too large
          !
          WRITE (msg,'("expected ",f10.6,", found ",f10.6)') &
             nelec, charge
          CALL errore( 'addusdens_r', 'WRONG CHARGE '//trim(msg)//&
                       ': ions may be overlapping or increase ecutrho', 1 )
          !
       ELSE
          !
          ! ... rescale the density to impose the correct number of electrons
          !
          rho_1(:,:) = rho_1(:,:) / charge * nelec
          !
       ENDIF
      ENDIF
      !
      CALL stop_clock( 'addusdens' )
      !
      RETURN
      !
    END SUBROUTINE addusdens_r
    !--------------------------------------------------------------------------
    SUBROUTINE calbec_rs_gamma ( ibnd, m, becp_r )

  !--------------------------------------------------------------------------
  !
  ! Subroutine written by Dario Rocca Stefano de Gironcoli, modified by O. Baris Malcioglu
  !
  ! Calculates becp_r in real space
  ! Requires BETASAVE (the beta functions at real space) calculated by betapointlist()
  ! (added to realus)
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
    USE uspp_param,            ONLY : nh, nhm
    USE fft_base,              ONLY : tg_gather, dffts
    USE mp_bands,              ONLY : me_bgrp, intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, m
    INTEGER :: iqs, iqsp, ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr, wi
    REAL(DP) :: bcr, bci
    REAL(DP), DIMENSION(:,:), INTENT(out) :: becp_r
    !
    REAL(DP), EXTERNAL :: ddot
    !
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN

     CALL errore( 'calbec_rs_gamma', 'task_groups not implemented', 1 )

    ELSE !non task groups part starts here

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    becp_r(:,ibnd)=0.d0
    IF ( ibnd+1 <= m ) becp_r(:,ibnd+1)=0.d0
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
             IF ( ityp(ia) == nt ) THEN
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
                   wr(:) = dble ( psic( box_beta(1:mbia,ia) ) )
                   wi(:) = aimag( psic( box_beta(1:mbia,ia) ) )
                   !print *, "betasave check", betasave(ia,ih,:)
                   ! box_beta contains explictly the points of the real space grid in
                   ! which the beta functions are differet from zero. Remember
                   ! that dble(psic) corresponds to ibnd, and aimag(psic) to ibnd+1:
                   ! this is the standard way to perform fourier transform in pwscf
                   ! in the gamma_only case
                   bcr  = ddot( mbia, betasave(ia,ih,:), 1, wr(:) , 1 )
                   bci  = ddot( mbia, betasave(ia,ih,:), 1, wi(:) , 1 )
                   ! in the previous two lines the real space integral is performed, using
                   ! few points of the real space mesh only
                   becp_r(ikb,ibnd)   = fac * bcr
                   IF ( ibnd+1 <= m ) becp_r(ikb,ibnd+1) = fac * bci
                   ! It is necessary to multiply by fac which to obtain the integral
                   ! in real space
                   !print *, becp_r(ikb,ibnd)
                   iqs = iqsp + 1
                   !
                ENDDO
                !
                DEALLOCATE( wr, wi )
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       !
    ENDIF
    CALL mp_sum( becp_r( :, ibnd ), intra_bgrp_comm )
    IF ( ibnd+1 <= m ) CALL mp_sum( becp_r( :, ibnd+1 ), intra_bgrp_comm )
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN

  END SUBROUTINE calbec_rs_gamma
    !
    SUBROUTINE calbec_rs_k ( ibnd, m )
    !--------------------------------------------------------------------------
    ! The k_point generalised version of calbec_rs_gamma. Basically same as above,
    ! but becp is used instead of becp_r, skipping the gamma point reduction
    ! derived from above by OBM 051108
    USE kinds,                 ONLY : DP
    USE cell_base,             ONLY : omega
    USE wavefunctions_module,  ONLY : psic
    USE ions_base,             ONLY : nat, ntyp => nsp, ityp
    USE uspp_param,            ONLY : nh, nhm
    USE becmod,                ONLY : bec_type, becp
    USE fft_base,              ONLY : tg_gather, dffts
    USE mp_bands,              ONLY : me_bgrp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, m
    INTEGER :: iqs, iqsp, ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr, wi
    REAL(DP) :: bcr, bci
    !COMPLEX(DP), allocatable, dimension(:) :: bt
    !integer :: ir, k
    !
    REAL(DP), EXTERNAL :: ddot
    !
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN

     CALL errore( 'calbec_rs_k', 'task_groups not implemented', 1 )

    ELSE !non task groups part starts here

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    becp%k(:,ibnd)=0.d0
       iqs = 1
       ikb = 0
       !
       DO nt = 1, ntyp
          !
           DO ia = 1, nat
             !
             IF ( ityp(ia) == nt ) THEN
                !
                mbia = maxbox_beta(ia)
                ALLOCATE( wr(mbia), wi(mbia) )
                DO ih = 1, nh(nt)
                   ! nh is the number of beta functions, or something similar
                   !
                   ikb = ikb + 1
                   iqsp = iqs+mbia-1
                   wr(:) = dble ( psic( box_beta(1:mbia,ia) ) )
                   wi(:) = aimag( psic( box_beta(1:mbia,ia) ) )
                   bcr  = ddot( mbia, betasave(ia,ih,:), 1, wr(:) , 1 )
                   bci  = ddot( mbia, betasave(ia,ih,:), 1, wi(:) , 1 )
                   becp%k(ikb,ibnd)   = fac * cmplx( bcr, bci,kind=DP)
                   iqs = iqsp + 1
                   !
                ENDDO
                DEALLOCATE( wr, wi )
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
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
    ! ... This routine applies the S matrix to m wavefunctions psi in real space 
    ! ... (in psic), and puts the results again in psic for backtransforming.
    ! ... Requires becp%r (calbecr in REAL SPACE) and betasave (from betapointlist
    ! ... in realus)
    ! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu
    ! WARNING ! for the sake of speed, no checks performed in this subroutine

      USE kinds,                  ONLY : DP
      USE cell_base,              ONLY : omega
      USE wavefunctions_module,   ONLY : psic
      USE ions_base,              ONLY : nat, ntyp => nsp, ityp
      USE uspp_param,             ONLY : nh
      USE lsda_mod,               ONLY : current_spin
      USE uspp,                   ONLY : qq
      USE becmod,                 ONLY : bec_type, becp
      USE fft_base,               ONLY : tg_gather, dffts
      USE mp_bands,               ONLY : me_bgrp
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, m
      !
      INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
      !
      REAL(DP), EXTERNAL :: ddot
      !
      CALL start_clock( 's_psir' )
      IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN
         CALL errore( 's_psir_gamma', 'task_groups not implemented', 1 )
      ELSE
      ! non task groups part starts here
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
                  DO jh = 1, nh(nt)
                     jkb = ikb + jh
                     w1(ih) = w1(ih) + qq(ih,jh,nt) * becp%r(jkb, ibnd)
                     IF ( ibnd+1 <= m ) w2(ih) = w2(ih) + qq(ih,jh,nt) * becp%r(jkb, ibnd+1)
                  ENDDO
               ENDDO
               !
               w1 = w1 * fac
               w2 = w2 * fac
               ikb = ikb + nh(nt)
               !
               DO ih = 1, nh(nt)
                  !
                  DO ir = 1, mbia
                     iqs = jqs + ir
                     psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*cmplx( w1(ih), w2(ih) ,kind=DP)
                  ENDDO
                  !
                  jqs = iqs
                  !
               ENDDO
               !
               DEALLOCATE( w1, w2 )
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
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
      USE becmod,                 ONLY : bec_type, becp
      USE fft_base,               ONLY : tg_gather, dffts
      USE mp_bands,               ONLY : me_bgrp
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, m
      !
      INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      COMPLEX(DP) , ALLOCATABLE, DIMENSION(:) :: w1
      !
      REAL(DP), EXTERNAL :: ddot
      !

      CALL start_clock( 's_psir' )
      IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN
        CALL errore( 's_psir_k', 'task_groups not implemented', 1 )
      ELSE
      !non task groups part starts here
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
                  DO jh = 1, nh(nt)
                     jkb = ikb + jh
                     w1(ih) = w1(ih) + qq(ih,jh,nt) * becp%k(jkb, ibnd)
                  ENDDO
               ENDDO
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
                  ENDDO
                  !
                  jqs = iqs
                  !
               ENDDO
               !
               DEALLOCATE( w1 )
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
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
  !    in array becp%r(nkb,m) (calculated by calbecr in REAL SPACE)
  ! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu
  ! WARNING ! for the sake of speed, no checks performed in this subroutine

  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE wavefunctions_module,   ONLY : psic
  USE ions_base,              ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,             ONLY : nh
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq
  USE becmod,                 ONLY : bec_type, becp
  USE fft_base,               ONLY : tg_gather, dffts
  USE mp_bands,               ONLY : me_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, m
  !
  INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
  !
  REAL(DP), EXTERNAL :: ddot
  !
  CALL start_clock( 'add_vuspsir' )

  IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN

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
                  w1(ih) = w1(ih) + deeq(ih,jh,ia,current_spin) * becp%r(jkb,ibnd)
                  IF ( ibnd+1 <= m )  w2(ih) = w2(ih) + deeq(ih,jh,ia,current_spin)* &
                       becp%r(jkb,ibnd+1)
                  !
               ENDDO
               !
            ENDDO
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
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + &
                       betasave(ia,ih,ir)*cmplx( w1(ih), w2(ih) ,kind=DP)
                  !
               ENDDO
                  !
               jqs = iqs
               !
            ENDDO
            !
            DEALLOCATE( w1, w2 )
            !
         ENDIF
         !
      ENDDO
      !
   ENDDO
   !
  ENDIF
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
  USE becmod,                 ONLY : bec_type, becp
  USE fft_base,               ONLY : tg_gather, dffts
  USE mp_bands,               ONLY : me_bgrp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, m
  !
  INTEGER :: ih, jh, iqs, jqs, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  !
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: w1
  !
  REAL(DP), EXTERNAL :: ddot
  !
  CALL start_clock( 'add_vuspsir' )

  IF( ( dffts%have_task_groups ) .and. ( m >= dffts%nogrp ) ) THEN
    CALL errore( 'add_vuspsir_k', 'task_groups not implemented', 1 )
  ELSE
   ! non task groups part starts here
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
                  w1(ih) = w1(ih) + deeq(ih,jh,ia,current_spin) * becp%k(jkb,ibnd)
                  !
               ENDDO
               !
            ENDDO
            !
            w1 = w1 * fac
            ikb = ikb + nh(nt)
            !
            DO ih = 1, nh(nt)
               !
               DO ir = 1, mbia
                  !
                  iqs = jqs + ir
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + &
                       betasave(ia,ih,ir)*w1(ih)
                  !
               ENDDO
                  !
               jqs = iqs
               !
            ENDDO
            !
            DEALLOCATE( w1)
            !
         ENDIF
         !
      ENDDO
      !
   ENDDO
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
  ! This driver subroutine transforms the given orbital using fft and puts the
  ! result in psic
  ! Warning! In order to be fast, no checks on the supplied data are performed!
  ! orbital: the orbital to be transformed
  ! ibnd: band index
  ! nbnd: total number of bands
    USE wavefunctions_module,     ONLY : psic
    USE gvecs,                  ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, tg_gather
    USE fft_interfaces,ONLY : invfft
    USE mp_bands,      ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform

    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    LOGICAL, OPTIONAL :: conserved
    ! if this flag is true, the orbital is stored in temporary memory

    !Internal temporary variables
    INTEGER :: j, idx, ioff
    LOGICAL :: use_tg

    !Task groups
    INTEGER :: v_siz

    !The new task group version based on vloc_psi
    !print *, "->Real space"
    CALL start_clock( 'fft_orbital' )
    !
    ! The following is dirty trick to prevent usage of task groups if
    ! the number of bands is smaller than the number of task groups
    !
    use_tg = dffts%have_task_groups
    dffts%have_task_groups = ( dffts%have_task_groups ) .and. ( nbnd >= dffts%nogrp )

    IF( dffts%have_task_groups ) THEN
        !

        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*dffts%nogrp, 2

           IF( idx + ibnd - 1 < nbnd ) THEN
              DO j = 1, npw_k(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =      orbital(j,idx+ibnd-1) +&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd)
                 tg_psic(nlsm(igk_k(j,1))+ioff) =conjg(orbital(j,idx+ibnd-1) -&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == nbnd ) THEN
              DO j = 1, npw_k(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =        orbital(j,idx+ibnd-1)
                 tg_psic(nlsm(igk_k(j,1))+ioff) = conjg( orbital(j,idx+ibnd-1))
              ENDDO
           ENDIF

           ioff = ioff + dffts%tg_nnr

        ENDDO
        !
        !
        CALL invfft ('Wave', tg_psic, dffts)
        !
        !
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
          IF (.not. allocated(tg_psic_temp)) ALLOCATE( tg_psic_temp( dffts%tg_nnr * dffts%nogrp ) )
          tg_psic_temp=tg_psic
         ENDIF
        ENDIF

    ELSE !Task groups not used
        !
        psic(:) = (0.d0, 0.d0)

        IF (ibnd < nbnd) THEN
           ! two ffts at the same time
           DO j = 1, npw_k(1)
              psic (nls (igk_k(j,1))) =       orbital(j, ibnd) + (0.0d0,1.d0)*orbital(j, ibnd+1)
              psic (nlsm(igk_k(j,1))) = conjg(orbital(j, ibnd) - (0.0d0,1.d0)*orbital(j, ibnd+1))
           ENDDO
        ELSE
           DO j = 1, npw_k(1)
              psic (nls (igk_k(j,1))) =       orbital(j, ibnd)
              psic (nlsm(igk_k(j,1))) = conjg(orbital(j, ibnd))
           ENDDO
        ENDIF
        !
       CALL invfft ('Wave', psic, dffts)
        !
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
           IF (.not. allocated(psic_temp) ) ALLOCATE (psic_temp(size(psic)))
           CALL zcopy(size(psic),psic,1,psic_temp,1)
         ENDIF
        ENDIF

    ENDIF

    dffts%have_task_groups = use_tg

    CALL stop_clock( 'fft_orbital' )

  END SUBROUTINE fft_orbital_gamma
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE bfft_orbital_gamma (orbital, ibnd, nbnd,conserved)
  !--------------------------------------------------------------------------
  !
  ! OBM 241008
  ! This driver subroutine -back- transforms the given orbital using fft using
  ! the already existent data in psic. 
  ! Warning! This subroutine does not reset the orbital, use carefully!
  ! Warning 2! In order to be fast, no checks on the supplied data are performed!
  ! Variables:
  ! orbital: the orbital to be transformed
  ! ibnd: band index
  ! nbnd: total number of bands
    USE wavefunctions_module,     ONLY : psic
    USE gvecs,                  ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, tg_gather
    USE fft_interfaces,ONLY : fwfft
    USE mp_bands,      ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform

    COMPLEX(DP),INTENT(out) :: orbital(:,:)

    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    !Internal temporary variables
    COMPLEX(DP) :: fp, fm
    INTEGER :: j, idx, ioff
    LOGICAL :: use_tg

    !Task groups
    INTEGER :: v_siz
    !print *, "->fourier space"
    CALL start_clock( 'bfft_orbital' )
    !New task_groups versions
    use_tg = dffts%have_task_groups
    dffts%have_task_groups = ( dffts%have_task_groups ) .and. ( nbnd >= dffts%nogrp )
    IF( dffts%have_task_groups ) THEN
       !
        CALL fwfft ('Wave', tg_psic, dffts )
        !
        ioff   = 0
        !
        DO idx = 1, 2*dffts%nogrp, 2
           !
           IF( idx + ibnd - 1 < nbnd ) THEN
              DO j = 1, npw_k(1)
                 fp= ( tg_psic( nls(igk_k(j,1)) + ioff ) +  &
                      tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk_k(j,1)) + ioff ) -  &
                      tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 orbital (j, ibnd+idx-1) =  cmplx( dble(fp), aimag(fm),kind=DP)
                 orbital (j, ibnd+idx  ) =  cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == nbnd ) THEN
              DO j = 1, npw_k(1)
                 orbital (j, ibnd+idx-1) =  tg_psic( nls(igk_k(j,1)) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
           !
        ENDDO
        !
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
          IF (allocated(tg_psic_temp)) DEALLOCATE( tg_psic_temp )
         ENDIF
        ENDIF

    ELSE !Non task_groups version
         !larger memory slightly faster
        CALL fwfft ('Wave', psic, dffts)


        IF (ibnd < nbnd) THEN

           ! two ffts at the same time
           DO j = 1, npw_k(1)
              fp = (psic (nls(igk_k(j,1))) + psic (nlsm(igk_k(j,1))))*0.5d0
              fm = (psic (nls(igk_k(j,1))) - psic (nlsm(igk_k(j,1))))*0.5d0
              orbital( j, ibnd)   = cmplx( dble(fp), aimag(fm),kind=DP)
              orbital( j, ibnd+1) = cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           DO j = 1, npw_k(1)
              orbital(j, ibnd)   =  psic (nls(igk_k(j,1)))
           ENDDO
        ENDIF
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
           IF (allocated(psic_temp) ) DEALLOCATE(psic_temp)
         ENDIF
        ENDIF
    ENDIF
    dffts%have_task_groups = use_tg
    !
    CALL stop_clock( 'bfft_orbital' )

  END SUBROUTINE bfft_orbital_gamma
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fft_orbital_k (orbital, ibnd, nbnd, ik, conserved)
  !--------------------------------------------------------------------------
  !
  ! OBM 110908
  ! This subroutine transforms the given orbital using fft and puts the result
  ! in psic
  ! Warning! In order to be fast, no checks on the supplied data are performed!
  ! orbital: the orbital to be transformed
  ! ibnd: band index
  ! nbnd: total number of bands
    USE kinds,                    ONLY : DP
    USE wavefunctions_module,     ONLY : psic
    USE gvecs,                    ONLY : nls, nlsm, doublegrid
    USE fft_base,                 ONLY : dffts
    USE fft_interfaces,           ONLY : invfft
    USE mp_bands,                 ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! Index of the band currently being transformed
                           nbnd,& ! Total number of bands you want to transform
                           ik     ! kpoint index of the bands

    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    ! Internal variables
    INTEGER :: ioff, idx
    LOGICAL :: use_tg

    CALL start_clock( 'fft_orbital' )
    use_tg = dffts%have_task_groups
    dffts%have_task_groups = ( dffts%have_task_groups ) .and. ( nbnd >= dffts%nogrp )

    IF( dffts%have_task_groups ) THEN
       !
       tg_psic = ( 0.D0, 0.D0 )
       ioff   = 0
       !
       DO idx = 1, dffts%nogrp
          !
          IF( idx + ibnd - 1 <= nbnd ) THEN
             !DO j = 1, size(orbital,1)
             tg_psic( nls( igk_k(:, ik) ) + ioff ) = orbital(:,idx+ibnd-1)
             !END DO
          ENDIF

          ioff = ioff + dffts%tg_nnr

       ENDDO
       !
       CALL invfft ('Wave', tg_psic, dffts)
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (.not. allocated(tg_psic_temp)) &
                  &ALLOCATE( tg_psic_temp( dffts%tg_nnr * dffts%nogrp ) )
             tg_psic_temp=tg_psic
          ENDIF
       ENDIF
       !
    ELSE  !non task_groups version
       !
       psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
       !
       psic(nls(igk_k(1:npw_k(ik), ik))) = orbital(1:npw_k(ik),ibnd)
       !
       CALL invfft ('Wave', psic, dffts)
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (.not. allocated(psic_temp) ) ALLOCATE (psic_temp(size(psic)))
             psic_temp=psic
          ENDIF
       ENDIF
       !
    ENDIF
    dffts%have_task_groups = use_tg
    CALL stop_clock( 'fft_orbital' )
  END SUBROUTINE fft_orbital_k
  !--------------------------------------------------------------------------
  SUBROUTINE bfft_orbital_k (orbital, ibnd, nbnd, ik, conserved)
    !--------------------------------------------------------------------------
    !
    ! OBM 110908
    ! This subroutine transforms the given orbital using fft and puts the result
    ! in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    ! orbital: the orbital to be transformed
    ! ibnd: band index
    ! nbnd: total number of bands
    USE wavefunctions_module,     ONLY : psic
    USE gvecs,                    ONLY : nls, nlsm, doublegrid
    USE kinds,                    ONLY : DP
    USE fft_base,                 ONLY : dffts
    USE fft_interfaces,           ONLY : fwfft
    USE mp_bands,                 ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! Index of the band currently being transformed
                           nbnd,& ! Total number of bands you want to transform
                           ik     ! kpoint index of the bands
    COMPLEX(DP),INTENT(out) :: orbital(:,:)
    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    ! Internal variables
    INTEGER :: ioff, idx
    LOGICAL :: use_tg

   CALL start_clock( 'bfft_orbital' )
   use_tg = dffts%have_task_groups
   dffts%have_task_groups = ( dffts%have_task_groups ) .and. ( nbnd >= dffts%nogrp )

    IF( dffts%have_task_groups ) THEN
       !
       CALL fwfft ('Wave', tg_psic, dffts)
       !
       ioff   = 0
       !
       DO idx = 1, dffts%nogrp
          !
          IF( idx + ibnd - 1 <= nbnd ) THEN
             orbital (:, ibnd+idx-1) = tg_psic( nls(igk_k(:,ik)) + ioff )
          ENDIF
          !
          ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
          !
       ENDDO
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (allocated(tg_psic_temp)) DEALLOCATE( tg_psic_temp )
          ENDIF
       ENDIF
       !
    ELSE !non task groups version
       !
       CALL fwfft ('Wave', psic, dffts)
       !
       orbital(1:npw_k(ik),ibnd) = psic(nls(igk_k(1:npw_k(ik),ik)))
       !
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (allocated(psic_temp) ) DEALLOCATE(psic_temp)
          ENDIF
       ENDIF
    ENDIF
    dffts%have_task_groups = use_tg
    CALL stop_clock( 'bfft_orbital' )

  END SUBROUTINE bfft_orbital_k

  !--------------------------------------------------------------------------
  SUBROUTINE v_loc_psir (ibnd, nbnd)
    !--------------------------------------------------------------------------
    ! Basically the same thing as v_loc but without implicit fft
    ! modified for real space implementation
    ! OBM 241008
    !
    USE wavefunctions_module,     ONLY : psic
    USE gvecs,       ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, tg_gather
    USE mp_bands,      ONLY : me_bgrp
    USE scf,           ONLY : vrs
    USE lsda_mod,      ONLY : current_spin


    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! Current index of the band currently being transformed
                           nbnd ! Total number of bands you want to transform
    !Internal temporary variables
    INTEGER :: j
    !Task groups
    REAL(DP),    ALLOCATABLE :: tg_v(:)
    INTEGER :: v_siz
    CALL start_clock( 'v_loc_psir' )

    IF( dffts%have_task_groups .and. nbnd >= dffts%nogrp  ) THEN
        IF (ibnd == 1 ) THEN
          CALL tg_gather( dffts, vrs(:,current_spin), tg_v )
          !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        ENDIF
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_psic (j) + tg_psic_temp (j) * tg_v(j)
        ENDDO
        !
        DEALLOCATE( tg_v )
     ELSE
        !   product with the potential v on the smooth grid
        !
        DO j = 1, dffts%nnr
           psic (j) = psic (j) + psic_temp (j) * vrs(j,current_spin)
        ENDDO
     ENDIF
  CALL stop_clock( 'v_loc_psir' )
  END SUBROUTINE v_loc_psir
    !--------------------------------------------------------------------------
  !
END MODULE realus
