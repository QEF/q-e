!
! Copyright (C) 2004-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
  ! ... Cleanup, GWW-specific stuff moved out by P. Giannozzi (2015)
  ! ... Computation of dQ/dtau_i needed for forces added by P. Giannozzi (2015)
  ! ... Some comments about the way some routines act added by S. de Gironcoli  (2015)
  ! ... extended to generic k by S. de Gironcoli (2016)
  !
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: boxrad(:) ! radius of boxes, does not depend on the grid
  ! Beta function in real space
  INTEGER,  ALLOCATABLE :: box_beta(:,:), maxbox_beta(:)
  REAL(DP), ALLOCATABLE :: betasave(:,:,:)
  REAL(DP), ALLOCATABLE :: boxrad_beta(:)
  REAL(DP), ALLOCATABLE :: boxdist_beta(:,:), xyz_beta(:,:,:)
  REAL(DP), ALLOCATABLE :: spher_beta(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: xkphase(:,:)  ! kpoint-related phase factor around each atom
  INTEGER               :: current_phase_kpoint=-1 ! the kpoint index for which the xkphase is currently set
                                                    ! negative initial value  means not set
  !General
  LOGICAL               :: real_space
  ! if true perform calculations in real spave
  LOGICAL               :: do_not_use_spline_inside_rinner = .false.
  INTEGER               :: real_space_debug = 0 ! FIXME: must disappear
  INTEGER               :: initialisation_level
  ! init_realspace_vars sets this to 3; qpointlist adds 5; betapointlist adds 7
  ! so the value should be 15 if the real space routine is initialised properly

  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  COMPLEX(DP), ALLOCATABLE :: psic_temp(:),tg_psic_temp(:) !Copies of psic and tg_psic
  COMPLEX(DP), ALLOCATABLE :: tg_vrs(:) !task groups linear V memory
  COMPLEX(DP), ALLOCATABLE :: psic_box_temp(:),tg_psic_box_temp(:)
  !
  ! Contains the augmentation functions and related quantities, in realspace for a single atom
  TYPE realsp_augmentation
     INTEGER :: maxbox = 0
        !  number of R points in the augmentation sphere of this atom
     INTEGER,ALLOCATABLE  :: box(:) 
        ! (maxbox) Index of point R in the global order of the R-space grid
     REAL(DP),ALLOCATABLE :: dist(:) 
        ! (maxbox) distances |R-tau_i| (box is centered on atom tau_i)
     REAL(DP),ALLOCATABLE :: xyz(:,:) 
        ! (3,maxbox) coordinates of R-tau_i
     REAL(DP),ALLOCATABLE :: qr(:,:) 
        ! (maxbox,number_of_q_funcs) the Q functions sampled over R points
  END TYPE realsp_augmentation
  ! Augmentation functions on the RHO (HARD) grid for all atoms
  TYPE(realsp_augmentation),POINTER :: tabp(:) => null()
  ! Augmentation functions on the SMOOTH grid for all atoms
  !TYPE(realsp_augmentation),POINTER :: tabs(:) => null()
  ! Augmentation functions on the EXX grid for all atoms
  TYPE(realsp_augmentation),POINTER :: tabxx(:) => null()
  !
  PRIVATE
  ! variables for real-space Q, followed by routines
  PUBLIC :: tabp, tabxx, boxrad, realsp_augmentation
  PUBLIC :: generate_qpointlist, qpointlist, addusdens_r, newq_r, &
       addusforce_r, real_space_dq, deallocate_realsp
  ! variables for real-space beta, followed by routines
  PUBLIC :: real_space, initialisation_level, real_space_debug, &
       tg_psic, betasave, maxbox_beta, box_beta
  PUBLIC :: betapointlist, init_realspace_vars, v_loc_psir, v_loc_psir_inplace
  PUBLIC :: invfft_orbital_gamma, fwfft_orbital_gamma, s_psir_gamma, &
            calbec_rs_gamma, add_vuspsir_gamma, invfft_orbital_k,    &
            fwfft_orbital_k, s_psir_k, calbec_rs_k, add_vuspsir_k
  !
  CONTAINS
  
    !------------------------------------------------------------------------
    SUBROUTINE generate_qpointlist
      !------------------------------------------------------------------------
      USE fft_base,     ONLY : dfftp, dffts
      USE funct,        ONLY : dft_is_hybrid
      USE gvecs,        ONLY : doublegrid
      USE io_global,    ONLY : stdout
      !USE exx,  ONLY : exx_fft
      IMPLICIT NONE
      !
      ! 1. initialize hard grid
      WRITE(stdout, '(/,5x,a)') "Initializing real-space augmentation for DENSE grid"
      CALL qpointlist(dfftp, tabp) 
      !
      ! NOTE: nothing is using the real space smooth grid at the moment, as EXX
      !       now uses its own custom grid. In case this is re-introduced, please
      !       also modify exx_fft_create to recycle this grid when ecutfock = ecutwfc
!       ! 2. initialize smooth grid (only for EXX at this moment)
!       IF ( dft_is_hybrid() ) THEN
!         IF(doublegrid)THEN
!           WRITE(stdout, '(5x,a)') "Initializing real-space augmentation for SMOOTH grid"
!           CALL qpointlist(dffts, tabs)
!         ELSE
!           ! smooth and rho grid are the same if not double grid
!           WRITE(stdout, '(7x,a)') " SMOOTH grid -> DENSE grid"
!           tabs => tabp
!         ENDIF
!       ENDIF
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
     USE control_flags,        ONLY : tqr
     USE fft_base,             ONLY : dtgs
     USE io_global,            ONLY : stdout


     IMPLICIT NONE

     INTEGER :: ik

     !print *, "<<<<<init_realspace_vars>>>>>>>"

     !real space, allocation for task group fft work arrays

     IF( dtgs%have_task_groups ) THEN
        !
        IF (allocated( tg_psic ) ) DEALLOCATE( tg_psic )
        !
        ALLOCATE( tg_psic( dtgs%tg_nnr * dtgs%nogrp ) )
        ALLOCATE( tg_vrs( dtgs%tg_nnr * dtgs%nogrp ) )
        !
     ENDIF
     !
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
    SUBROUTINE deallocate_realsp()
      !------------------------------------------------------------------------
      !
      USE gvecs,      ONLY : doublegrid
      IMPLICIT NONE
      !
      IF ( allocated( boxrad ) )  DEALLOCATE( boxrad )
      !      
      CALL deallocate_realsp_aug ( tabp )
!       IF ( doublegrid ) THEN
!          CALL deallocate_realsp_aug ( tabs )
!       ELSE
!          NULLIFY(tabs)
!       END IF
      !
    END SUBROUTINE deallocate_realsp

    !------------------------------------------------------------------------
    SUBROUTINE deallocate_realsp_aug ( tab )
      !------------------------------------------------------------------------
      TYPE(realsp_augmentation), POINTER :: tab(:)
      INTEGER :: ia
      !
      IF ( associated( tab ) ) THEN
         DO ia = 1, SIZE(tab)
            IF(allocated(tab(ia)%qr ))  DEALLOCATE(tab(ia)%qr)
            IF(allocated(tab(ia)%box))  DEALLOCATE(tab(ia)%box)
            IF(allocated(tab(ia)%dist)) DEALLOCATE(tab(ia)%dist)
            IF(allocated(tab(ia)%xyz))  DEALLOCATE(tab(ia)%xyz)
            tab(ia)%maxbox = 0
         ENDDO
         DEALLOCATE(tab)
      ENDIF
      !
    END SUBROUTINE deallocate_realsp_aug
    !
    !------------------------------------------------------------------------
    SUBROUTINE qpointlist(dfft, tabp)
      !------------------------------------------------------------------------
      !
      ! ... This subroutine computes and stores into "tabp" the values of Q(r)
      ! ... on atom-centered boxes in the real-space grid of descriptor "dfft".
      ! ... Must be called at startup and every time atoms (tau_i) move.
      ! ... A set of spherical boxes are computed for each atom. The Q functions
      ! ... Q(r-tau_i) are then interpolated on the real-space grid points.
      ! ... On output:
      ! ... - boxrad    contains the radii of the boxes for each atom type
      ! ... - tabp(nat) contains pointers to structure tabp, each contaning
      ! ...   % maxbox = number of real-space grid points in the atom-centered box
      ! ...   % qr(maxbox) = interpolated values of Q(r) on the points of the box
      ! ...   % box(maxbox)= index of grid points in the real-space grid "dfft"
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi, fpi, eps16
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, omega, alat
      USE uspp,       ONLY : okvan
      USE uspp_param, ONLY : upf, nh
      USE atom,       ONLY : rgrid
      USE fft_types,  ONLY : fft_type_descriptor
      USE mp_bands,   ONLY : me_bgrp
      !
      IMPLICIT NONE
      !
      TYPE(fft_type_descriptor),INTENT(in)    :: dfft
      TYPE(realsp_augmentation),POINTER,INTENT(inout) :: tabp(:)
      !
      INTEGER               :: ia, mbia
      INTEGER               :: indm, ih, jh, ijh
      INTEGER               :: roughestimate, nt, l
      INTEGER,  ALLOCATABLE :: buffpoints(:)
      INTEGER               :: idx0, idx, ir
      INTEGER               :: i, j, k, ipol, ijv
      REAL(DP)              :: distsq, posi(3)
      REAL(DP), ALLOCATABLE :: boxdist(:), xyz(:,:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz, aux
      REAL(DP)              :: inv_nr1, inv_nr2, inv_nr3, boxradsq_ia
      !
      initialisation_level = 3
      IF ( .not. okvan ) RETURN
      !
      CALL start_clock( 'realus' )
      !
      ! ... tabp is deallocated here to free the memory for the buffers
      !
      CALL deallocate_realsp_aug ( tabp )
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
                  aux = maxval(abs( upf(nt)%qfuncl(indm,ijv,:) ))
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
      ALLOCATE( buffpoints( roughestimate ) )
      ALLOCATE( boxdist (   roughestimate ) )
      ALLOCATE( xyz( 3, roughestimate ) )
      !
      ! idx0 = starting index of real-space FFT arrays for this processor
      idx0 = dfft%nr1x*dfft%nr2x * dfft%ipp(me_bgrp+1)
      !
      inv_nr1 = 1.D0 / dble( dfft%nr1 )
      inv_nr2 = 1.D0 / dble( dfft%nr2 )
      inv_nr3 = 1.D0 / dble( dfft%nr3 )
      !
      DO ia = 1, nat
         !
         CALL start_clock( 'realus:boxes' )
         !
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
         !
         ! ... here we find the points in the box surrounding atom "ia"
         !
         buffpoints(:) = 0
         boxdist(:) = 0.D0
         !
         boxradsq_ia = boxrad(nt)**2
         !
         mbia = 0
         !
         ! ... The do loop includes only planes that belong to this processor
         !
         DO ir = 1, dfft%nr1x*dfft%nr2x * dfft%npl
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
            ! ... do not include points outside the physical range 
            !
            IF ( i >= dfft%nr1 .OR. j >= dfft%nr2 .OR. k >= dfft%nr3 ) CYCLE
            !
            DO ipol = 1, 3
               posi(ipol) = dble( i )*inv_nr1*at(ipol,1) + &
                            dble( j )*inv_nr2*at(ipol,2) + &
                            dble( k )*inv_nr3*at(ipol,3)
            ENDDO
            !
            posi(:) = posi(:) - tau(:,ia)
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
               mbia = mbia + 1
               IF( mbia > roughestimate ) &
                    CALL errore('qpointlist', 'rough-estimate is too rough', 3)
               buffpoints(mbia)    = ir
               boxdist(mbia)       = sqrt( distsq )*alat
               xyz(:,mbia)         = posi(:)*alat
               !
            ENDIF
            !
         ENDDO
         !
         IF ( mbia == 0 ) CYCLE
         !
         ! ... store points in the appropriate place
         !
         tabp(ia)%maxbox = mbia
         ALLOCATE( tabp(ia)%box(mbia) )
         tabp(ia)%box(:) = buffpoints(1:mbia)
         ALLOCATE( tabp(ia)%dist(mbia) )
         tabp(ia)%dist(:) = boxdist(1:mbia)
         ALLOCATE( tabp(ia)%xyz(3,mbia) )
         tabp(ia)%xyz(:,:) = xyz(:,1:mbia)
         CALL stop_clock( 'realus:boxes' )
         !
         ! ... compute the Q functions for this box
         !
         CALL real_space_q( nt, ia, mbia, tabp )
         !
      ENDDO
      !
      DEALLOCATE( xyz )
      DEALLOCATE( boxdist )
      DEALLOCATE( buffpoints )
      CALL stop_clock( 'realus' )
      !
    END SUBROUTINE qpointlist
    !
    !------------------------------------------------------------------------
    SUBROUTINE real_space_q( nt, ia, mbia, tab )
      !------------------------------------------------------------------------
      !
      ! ... Compute Q_lm(r-tau_ia) centered around  atom ia of type nt. 
      ! ... First we compute the radial Q(r) (qtot) for each l; then we
      ! ... we interpolate it in our mesh (contained in tabp); finally
      ! ... we multiply by the spherical harmonics and store it into tabp
      !
      ! ... Q is read from pseudo and it is divided into two parts:
      ! ... in the inner radius a polinomial representation is known and so
      ! ... strictly speaking we do not use interpolation but just compute
      ! ... the correct value
      !
      USE constants,  ONLY : eps16
      USE uspp,       ONLY : indv, nhtol, nhtolm, ap, nhtoj
      USE uspp_param, ONLY : upf, lmaxq, nh
      USE atom,       ONLY : rgrid
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ia, nt, mbia
      TYPE(realsp_augmentation), INTENT(INOUT), POINTER :: tab(:)
      !
      INTEGER  :: l, nb, mb, ijv, lllnbnt, lllmbnt, ir, nfuncs, lm, &
           i, ijh, ih, jh, ipol
      REAL(dp) :: first, second, qtot_int
      REAL(dp), ALLOCATABLE :: qtot(:), dqtot(:), xsp(:), wsp(:), &
           rl2(:), spher(:,:)

      CALL start_clock( 'realus:tabp' )
      !
      nfuncs = nh(nt)*(nh(nt)+1)/2
      ALLOCATE(tab(ia)%qr(mbia, nfuncs))
      !
      ! ... compute the spherical harmonics
      !
      ALLOCATE( spher(mbia, lmaxq**2) )
      spher (:,:) = 0.0_dp
      !
      ALLOCATE( rl2( mbia ) )
      DO ir = 1, mbia
         rl2(ir) = tab(ia)%xyz(1,ir)**2 + &
                   tab(ia)%xyz(2,ir)**2 + &
                   tab(ia)%xyz(3,ir)**2
      ENDDO
      CALL ylmr2 ( lmaxq**2, mbia, tab(ia)%xyz, rl2, spher )
      DEALLOCATE( rl2 )

      tab(ia)%qr(:,:)=0.0_dp
      ALLOCATE( qtot(upf(nt)%kkbeta), dqtot(upf(nt)%kkbeta) )
      !
      ! ... variables used for spline interpolation
      !
      ALLOCATE( xsp( upf(nt)%kkbeta ), wsp( upf(nt)%kkbeta ) )
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
               IF ( .not. ( l >= abs( lllnbnt - lllmbnt ) .and. &
                            l <= lllnbnt + lllmbnt        .and. &
                            mod( l + lllnbnt + lllmbnt, 2 ) == 0 ) ) CYCLE
               !
               IF( upf(nt)%tvanp ) THEN
                  qtot(2:upf(nt)%kkbeta) = &
                       upf(nt)%qfuncl(2:upf(nt)%kkbeta,ijv,l) &
                       / rgrid(nt)%r(2:upf(nt)%kkbeta)**2
                  ! prevents division by zero if r(1)=0
                  IF (rgrid(nt)%r(1)< eps16) THEN
                     qtot(1) = qtot(2)
                  ELSE
                     qtot(1) = upf(nt)%qfuncl(1,ijv,l)/rgrid(nt)%r(1)**2
                  END IF
               ENDIF
               !
               ! ... compute the first derivative
               !
               ! ... analytical derivative up to r = rinner, numerical beyond
               !
               CALL radial_gradient(qtot, dqtot, rgrid(nt)%r, upf(nt)%kkbeta, 1)
               !
               ! ... we need the first and second derivatives in the first point
               !
               first = dqtot(1)
               !
               ! ... if we don't have the analitical coefficients, try the same
               ! ... numerically (note that setting first=0.0 and second=0.0
               ! ... makes almost no difference) - wsp is used as work space
               !
               CALL radial_gradient(dqtot, wsp, rgrid(nt)%r, upf(nt)%kkbeta, 1)
               second = wsp(1) ! second derivative in first point
               !
               ! ... call spline for interpolation of Q(r)
               !
               CALL spline( xsp, qtot, first, second, wsp )
               !
               DO ir = 1, tab(ia)%maxbox
                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IF .not. do_not_use_spline_inside_rinner IT DIFFERS ! it changes the energy by 1.4d-5
                  IF ( tab(ia)%dist(ir) < upf(nt)%rinner(l+1) .and. do_not_use_spline_inside_rinner) THEN
                     !
                     ! ... if in the inner radius just compute the
                     ! ... polynomial
                     !
                     CALL setqfnew( upf(nt)%nqf, upf(nt)%qfcoef(1,l+1,nb,mb),&
                          1, tab(ia)%dist(ir), l, 0, qtot_int )
                     !
                  ELSE
                     !
                     ! ... spline interpolation
                     !
                     qtot_int = splint( xsp, qtot, wsp, tab(ia)%dist(ir) )
                     !
                  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                  !
                  ijh = 0
                  DO ih = 1, nh(nt)
                     DO jh = ih, nh(nt)
                        !
                        ijh = ijh + 1
                        !
                        IF ( .not.( nb == indv(ih,nt) .and. &
                                    mb == indv(jh,nt) ) ) CYCLE
                        !
                        DO lm = l**2+1, (l+1)**2
                           tab(ia)%qr(ir,ijh) = tab(ia)%qr(ir,ijh) + &
                                        qtot_int*spher(ir,lm)*&
                                        ap(lm,nhtolm(ih,nt),nhtolm(jh,nt))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE( qtot, dqtot )
      DEALLOCATE( wsp, xsp )
      DEALLOCATE( spher )
      !
      CALL stop_clock( 'realus:tabp' )
      !
    END SUBROUTINE real_space_q
    !
    !------------------------------------------------------------------------
    SUBROUTINE real_space_dq( nt, ia, mbia, nfuncs, dqr )
      !------------------------------------------------------------------------
      !
      ! ... compute dQ(r-tau_i)/dtau_{i,a} using the formula
      !     dQ(r-tau_i)/dtau_{i,a}= - dq/dr(|r-tau_i|)(r_a-tau_{i,a})/|r-tau_i| 
      !                               * Y_lm(r-tau_i)
      !                             - q(|r-tau_i|) * (dY_lm(r-tau_i)/dtau_{i,a})
      ! ... This routine follows the same logic of real_space_q 
      !
      USE constants,  ONLY : eps8, eps16
      USE uspp,       ONLY : indv, nhtol, nhtolm, ap, nhtoj
      USE uspp_param, ONLY : upf, lmaxq, nh
      USE atom,       ONLY : rgrid
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ia, nt, mbia, nfuncs
      REAL(dp),INTENT(OUT):: dqr(mbia,nfuncs,3)
      !
      INTEGER  :: l, nb, mb, ijv, lllnbnt, lllmbnt, ir, lm, &
           i, ijh, ih, jh, ipol
      REAL(dp) :: first, second, qtot_int, dqtot_int, dqxyz(3)
      REAL(dp), ALLOCATABLE :: qtot(:), dqtot(:), xsp(:), wsp(:,:), &
           rl2(:), spher(:,:), dspher(:,:,:)
      !
      CALL start_clock( 'realus:tabp' )
      !
      ! ... compute the spherical harmonics
      !
      ALLOCATE( spher(mbia, lmaxq**2), dspher(mbia, lmaxq**2, 3) )
      spher (:,:) = 0.0_dp
      dspher(:,:,:) = 0.0_dp
      !
      ALLOCATE( rl2( mbia ) )
      DO ir = 1, mbia
         rl2(ir) = tabp(ia)%xyz(1,ir)**2 + &
                   tabp(ia)%xyz(2,ir)**2 + &
                   tabp(ia)%xyz(3,ir)**2
      ENDDO
      CALL ylmr2 ( lmaxq**2, mbia, tabp(ia)%xyz, rl2, spher )
      do ipol = 1, 3
         CALL dylmr2( lmaxq**2, mbia, tabp(ia)%xyz, rl2, dspher(1,1,ipol),ipol )
      END DO
      DEALLOCATE( rl2 )

      dqr(:,:,:)=0.0_dp
      ALLOCATE( qtot(upf(nt)%kkbeta), dqtot(upf(nt)%kkbeta) )
      !
      ! ... variables used for spline interpolation
      !
      ALLOCATE( xsp( upf(nt)%kkbeta ), wsp( upf(nt)%kkbeta, 2 ) )
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
               IF ( .not. ( l >= abs( lllnbnt - lllmbnt ) .and. &
                            l <= lllnbnt + lllmbnt        .and. &
                            mod( l + lllnbnt + lllmbnt, 2 ) == 0 ) ) CYCLE
               !
               IF( upf(nt)%tvanp ) THEN
                  qtot(1:upf(nt)%kkbeta) = &
                       upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l) &
                       / rgrid(nt)%r(1:upf(nt)%kkbeta)**2
                  if (rgrid(nt)%r(1)< eps16) qtot(1) = qtot(2)
               ENDIF
               !
               ! ... compute the first derivative
               !
               ! ... analytical derivative up to r = rinner, numerical beyond
               !
               CALL radial_gradient(qtot, dqtot, rgrid(nt)%r, upf(nt)%kkbeta, 1)
               !
               ! ... we need the first and second derivatives in the first point
               !
               first = dqtot(1)
               !
               ! ... if we don't have the analitical coefficients, try the same
               ! ... numerically (note that setting first=0.0 and second=0.0
               ! ... makes almost no difference) - wsp is used as work space
               !
               CALL radial_gradient(dqtot, wsp, rgrid(nt)%r, upf(nt)%kkbeta, 1)
               second = wsp(1,1) ! second derivative in first point
               !
               ! ... call spline for interpolation of Q(r)
               !
               CALL spline( xsp, qtot, first, second, wsp(:,1) )
               !
               ! ... same for dQ/dr (third derivative at first point is set to zero)
               !
               CALL spline( xsp, dqtot, second, 0.0_dp, wsp(:,2) )
               !
               DO ir = 1, mbia
                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IF .not. do_not_use_spline_inside_rinner IT DIFFERS ! it changes the force on Oxygen by 0.004
                  IF ( tabp(ia)%dist(ir) < upf(nt)%rinner(l+1) .and. do_not_use_spline_inside_rinner ) THEN
                     !
                     ! ... if in the inner radius just compute the
                     ! ... polynomial
                     !
                     CALL setqfnew( upf(nt)%nqf, upf(nt)%qfcoef(1,l+1,nb,mb),&
                          1, tabp(ia)%dist(ir), l, 0, qtot_int )
                     CALL setdqf( upf(nt)%nqf, upf(nt)%qfcoef(1,l+1,nb,mb), &
                          1, tabp(ia)%dist(ir), l, dqtot_int)
                  ELSE
                     !
                     ! ... spline interpolation
                     !
                     qtot_int = splint( xsp, qtot, wsp(:,1), tabp(ia)%dist(ir) )
                     dqtot_int= splint( xsp,dqtot, wsp(:,2), tabp(ia)%dist(ir) )
                     !
                  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                  !
                  ! ... prevent floating-point error if dist = 0
                  !
                  IF ( tabp(ia)%dist(ir) > eps8 ) THEN
                     dqxyz(:) = dqtot_int*tabp(ia)%xyz(:,ir) / &
                                     tabp(ia)%dist(ir)
                  ELSE
                     dqxyz(:) = 0.0_dp
                  ENDIF
                  !
                  ijh = 0
                  DO ih = 1, nh(nt)
                     DO jh = ih, nh(nt)
                        !
                        ijh = ijh + 1
                        !
                        IF ( .not.( nb == indv(ih,nt) .and. &
                                    mb == indv(jh,nt) ) ) CYCLE
                        !
                        DO lm = l**2+1, (l+1)**2
                           DO ipol = 1,3
                              dqr(ir,ijh,ipol) = &
                                   dqr(ir,ijh,ipol) - &
                                   ( qtot_int*dspher(ir,lm,ipol) +  &
                                     dqxyz(ipol)*spher(ir,lm) ) * &
                                     ap(lm,nhtolm(ih,nt),nhtolm(jh,nt))
                           END DO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               !
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE( qtot, dqtot )
      DEALLOCATE( wsp, xsp )
      DEALLOCATE( dspher, spher )
      !
      CALL stop_clock( 'realus:tabp' )
      !
    END SUBROUTINE real_space_dq
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
      USE constants,  ONLY : pi
      USE control_flags, ONLY : gamma_only
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
      INTEGER               :: ia, it, mbia
      INTEGER               :: indm, inbrx, idimension, ih
      INTEGER               :: roughestimate, goodestimate, lamx2, nt
      INTEGER,  ALLOCATABLE :: buffpoints(:,:)
      REAL(DP), ALLOCATABLE :: buffdist(:,:)
      REAL(DP), ALLOCATABLE :: buff_xyz_beta(:,:,:)
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
      IF ( allocated( xyz_beta ) ) DEALLOCATE( xyz_beta )
      ALLOCATE( buff_xyz_beta( 3, roughestimate, nat ) )
      !
      buffpoints(:,:) = 0
      buffdist(:,:) = 0.D0
      !
      IF ( .not.allocated( maxbox_beta ) ) ALLOCATE( maxbox_beta( nat ) )
      !
      maxbox_beta(:) = 0
      !
      ! idx0 = starting index of real-space FFT arrays for this processor
      ! The beta functions are treated on smooth grid
      !
      idx0 = dffts%nr1x*dffts%nr2x * dffts%ipp(me_bgrp+1)
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
               buff_xyz_beta(:,mbia,ia) = posi(:)*alat
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
      ALLOCATE( xyz_beta ( 3, goodestimate, nat ) )
      ALLOCATE( box_beta    ( goodestimate, nat ) )
      ALLOCATE( boxdist_beta( goodestimate, nat ) )
      ALLOCATE( xkphase     ( goodestimate, nat ) )
      !
      xyz_beta(:,:,:)   = buff_xyz_beta(:,1:goodestimate,:)
      box_beta(:,:)     = buffpoints(1:goodestimate,:)
      boxdist_beta(:,:) = buffdist(1:goodestimate,:)
      call set_xkphase(1)
      !
      DEALLOCATE( buff_xyz_beta )
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
      if (gamma_only) DEALLOCATE( xyz_beta )
      !
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:tabp' )
      !
      ! ... let's do the main work
      !
      ALLOCATE( betasave( goodestimate, nhm, nat )  )
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
               betasave(ir,ih,ia) = qtot_int*spher_beta(ir,lm,ia) !spher_beta is the Y_lm in point ir for atom ia
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
      INTEGER               :: mbia
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
         DO ia = 1, nat
            !
            mbia = tabp(ia)%maxbox
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            IF ( .not. upf(nt)%tvanp ) CYCLE
            !
            DO ih = 1, nh(nt)
               DO jh = ih, nh(nt)
                  DO ir = 1, mbia
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
      INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, mbia
      CHARACTER(len=80) :: msg
      REAL(DP) :: charge
      REAL(DP) :: tolerance
      !
      !
      IF ( .not. okvan ) RETURN
      tolerance = 1.d-3
      ! Charge loss with real-space betas is worse
      IF ( real_space ) tolerance = 1.d-2
      !
      CALL start_clock( 'addusdens' )
      !
      DO is = 1, nspin_mag
         !
         DO ia = 1, nat
            !
            mbia = tabp(ia)%maxbox
            IF ( mbia == 0 ) CYCLE
            !
            nt = ityp(ia)
            IF ( .not. upf(nt)%tvanp ) CYCLE
            !
            ijh = 0
            DO ih = 1, nh(nt)
               DO jh = ih, nh(nt)
                  ijh = ijh + 1
                  DO ir = 1, mbia
                     irb = tabp(ia)%box(ir)
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
    !----------------------------------------------------------------------
    SUBROUTINE addusforce_r (forcenl)
      !----------------------------------------------------------------------
      !
      !   This routine computes the contribution to atomic forces due
      !   to the dependence of the Q function on the atomic position.
      !   F_j,at=-\int sum_lm V^*(r) dQ_lm(r)/dtau_at dr becsum(lm,at)
      !   where becsum(lm,at) = sum_i <psi_i|beta_l>w_i<beta_m|psi_i>
      !   On output: the contribution is added to forcenl
      !
      USE kinds,      ONLY : DP
      USE cell_base,  ONLY : omega
      USE ions_base,  ONLY : nat, ntyp => nsp, ityp
      USE fft_base,   ONLY : dfftp
      USE noncollin_module,   ONLY : nspin_mag
      USE scf,        ONLY : v, vltot
      USE uspp,       ONLY : becsum, okvan
      USE uspp_param, ONLY : upf,  nh
      USE mp_bands,   ONLY : intra_bgrp_comm
      USE mp,         ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT) :: forcenl (3, nat)
      !
      INTEGER :: na, nt, ih, jh, ijh, nfuncs, mbia, ir, is, irb
      REAL(dp), ALLOCATABLE:: dqr(:,:,:), forceq(:,:)
      REAL(dp) :: dqrforce(3)
      !
      IF (.not.okvan) RETURN
      !
      ALLOCATE ( forceq(3,nat) )
      forceq(:,:) = 0.0_dp
      !
      DO na = 1, nat
         nt = ityp(na)
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         !
         mbia = tabp(na)%maxbox
         IF ( mbia == 0 ) CYCLE
         nfuncs = nh(nt)*(nh(nt)+1)/2
         ALLOCATE ( dqr(mbia,nfuncs,3) )
         !
         CALL real_space_dq( nt, na, mbia, nfuncs, dqr )
         !
         ijh = 0
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
               ijh = ijh + 1
               DO is = 1, nspin_mag
                  dqrforce(:) = 0.0_dp
                  IF (nspin_mag==4.and.is/=1) THEN
                     DO ir = 1, mbia
                        irb = tabp(na)%box(ir)
                        dqrforce(:) = dqrforce(:) + dqr(ir,ijh,:) * &
                             v%of_r(irb,is)
                     ENDDO
                  ELSE
                     DO ir = 1, mbia
                        irb = tabp(na)%box(ir)
                        dqrforce(:) = dqrforce(:) + dqr(ir,ijh,:) * &
                             (vltot(irb) + v%of_r(irb,is) )
                     ENDDO
                  END IF
                  forceq(:,na) = forceq(:,na) - dqrforce(:)*becsum(ijh,na,is)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( dqr )
      ENDDO
      !
      CALL mp_sum ( forceq, intra_bgrp_comm )
      !
      forcenl(:,:) = forcenl(:,:) + forceq(:,:) * &
           omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      DEALLOCATE (forceq)
      !
      RETURN
    END SUBROUTINE addusforce_r
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_xkphase(ik)
    !--------------------------------------------------------------------------
    ! in the calculation of becp or when performing add_vuspsir the wavefunction
    ! psi_k and  not its periodic part (which is what we get from the FFT) should be
    ! used. A k-dependent phase exp(-xk(current_k*(r-tau(ia))) ) must be added
    !
    USE kinds,      ONLY : DP
    USE klist,      ONLY : xk
    USE cell_base,  ONLY : tpiba
    USE ions_base,  ONLY : nat

    IMPLICIT NONE
  
    INTEGER, INTENT (IN) :: ik

    INTEGER :: ia, mbia, ir
    REAL(DP) :: arg

    if (.not.allocated ( xkphase ) ) call errore ('set_xkphase',' array not allocated yes',1)
    if (ik .eq. current_phase_kpoint ) return
    !
    DO ia = 1, nat
       mbia = maxbox_beta(ia)
       do ir =1, mbia
          arg = ( xk(1,ik) * xyz_beta(1,ir,ia) + &
                  xk(2,ik) * xyz_beta(2,ir,ia) + &
                  xk(3,ik) * xyz_beta(3,ir,ia) ) * tpiba
          xkphase( ir, ia ) = CMPLX(COS(arg),-SIN(arg),KIND=dp)
       end do
    end do
    !
    current_phase_kpoint = ik
    !
    return
    END SUBROUTINE set_xkphase

    !--------------------------------------------------------------------------
    SUBROUTINE calbec_rs_gamma ( ibnd, last, becp_r )

  !--------------------------------------------------------------------------
  !
  ! Subroutine written by Dario Rocca, Stefano de Gironcoli, modified by O. Baris Malcioglu
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
    USE fft_base,              ONLY : dffts, dtgs
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, last
    INTEGER :: ikb, nt, ia, ih, mbia
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
    IF( dtgs%have_task_groups ) THEN

     CALL errore( 'calbec_rs_gamma', 'task_groups not implemented', 1 )

    ELSE !non task groups part starts here

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    becp_r(:,ibnd)=0.d0
    IF ( ibnd+1 <= last ) becp_r(:,ibnd+1)=0.d0
    ! Clearly for an odd number of bands for ibnd=nbnd=last you don't have
    ! anymore bands, and so the imaginary part equal zero
    !
       !
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
                ! just working arrays to order the points in the clever way
                wr(:) = dble ( psic( box_beta(1:mbia,ia) ) )
                wi(:) = aimag( psic( box_beta(1:mbia,ia) ) )
                !
                !
                DO ih = 1, nh(nt)
                   ! nh is the number of beta functions, or something similar
                   !
                   ikb = ikb + 1
                   !print *, "betasave check", betasave(ia,ih,:)
                   ! box_beta contains explictly the points of the real space grid in
                   ! which the beta functions are differet from zero. Remember
                   ! that dble(psic) corresponds to ibnd, and aimag(psic) to ibnd+1:
                   ! this is the standard way to perform fourier transform in pwscf
                   ! in the gamma_only case
                   bcr  = ddot( mbia, betasave(:,ih,ia), 1, wr(:) , 1 )
                   bci  = ddot( mbia, betasave(:,ih,ia), 1, wi(:) , 1 )
                   ! in the previous two lines the real space integral is performed, using
                   ! few points of the real space mesh only
                   becp_r(ikb,ibnd)   = fac * bcr
                   IF ( ibnd+1 <= last ) becp_r(ikb,ibnd+1) = fac * bci
                   ! It is necessary to multiply by fac which to obtain the integral
                   ! in real space
                   !print *, becp_r(ikb,ibnd)
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
    IF ( ibnd+1 <= last ) CALL mp_sum( becp_r( :, ibnd+1 ), intra_bgrp_comm )
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN

  END SUBROUTINE calbec_rs_gamma
    !
    SUBROUTINE calbec_rs_k ( ibnd, last )
    !--------------------------------------------------------------------------
    ! The k_point generalised version of calbec_rs_gamma. Basically same as above,
    ! but becp is used instead of becp_r, skipping the gamma point reduction
    ! derived from above by OBM 051108
    ! k-point phase factor fixed by SdG 030816
    !
    USE kinds,                 ONLY : DP
    USE wvfct,                 ONLY : current_k
    USE cell_base,             ONLY : omega
    USE wavefunctions_module,  ONLY : psic
    USE ions_base,             ONLY : nat, ntyp => nsp, ityp
    USE uspp_param,            ONLY : nh, nhm
    USE becmod,                ONLY : bec_type, becp
    USE fft_base,              ONLY : dffts, dtgs
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, last
    INTEGER :: ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr, wi
    REAL(DP) :: bcr, bci
    !COMPLEX(DP), allocatable, dimension(:) :: bt
    integer :: ir
    !
    REAL(DP), EXTERNAL :: ddot
    !
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( dtgs%have_task_groups ) CALL errore( 'calbec_rs_k', 'task_groups not implemented', 1 )

    call set_xkphase(current_k)

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    becp%k(:,ibnd)=0.d0
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
                   wr(:) = dble ( psic( box_beta(1:mbia,ia) ) * CONJG(xkphase(1:mbia,ia)))
                   wi(:) = aimag( psic( box_beta(1:mbia,ia) ) * CONJG(xkphase(1:mbia,ia)))
                   bcr  = ddot( mbia, betasave(:,ih,ia), 1, wr(:) , 1 )
                   bci  = ddot( mbia, betasave(:,ih,ia), 1, wi(:) , 1 )
                   becp%k(ikb,ibnd)   = fac * cmplx( bcr, bci,kind=DP)
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
    CALL mp_sum( becp%k( :, ibnd ), intra_bgrp_comm )
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN

  END SUBROUTINE calbec_rs_k
    !--------------------------------------------------------------------------
    SUBROUTINE s_psir_gamma ( ibnd, last )
    !--------------------------------------------------------------------------
    !
    ! ... This routine applies the S matrix to wfc ibnd (and wfc ibnd+1 if <= last)
    ! ... stored in real space in psic, and puts the results again in psic for 
    ! ... backtransforming.
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
      USE fft_base,               ONLY : dffts, dtgs
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, last
      !
      INTEGER :: ih, jh, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
      !
      REAL(DP), EXTERNAL :: ddot
      !
      CALL start_clock( 's_psir' )

      IF( dtgs%have_task_groups ) CALL errore( 's_psir_gamma', 'task_groups not implemented', 1 )

      !
      fac = sqrt(omega)
      !
      ikb = 0
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
                     IF ( ibnd+1 <= last ) w2(ih) = w2(ih) + qq(ih,jh,nt) * becp%r(jkb, ibnd+1)
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
                     psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ir,ih,ia)*cmplx( w1(ih), w2(ih) ,kind=DP)
                  ENDDO
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
      CALL stop_clock( 's_psir' )
      !
      RETURN
      !
  END SUBROUTINE s_psir_gamma
  !
  SUBROUTINE s_psir_k ( ibnd, last )
  !--------------------------------------------------------------------------
  ! Same as s_psir_gamma but for generalised k point scheme i.e.:
  ! 1) Only one band is considered at a time
  ! 2) Becp is a complex entity now
  ! Derived from s_psir_gamma by OBM 061108
  ! k-point phase factor fixed by SdG 030816
      USE kinds,                  ONLY : DP
      USE wvfct,                  ONLY : current_k
      USE cell_base,              ONLY : omega
      USE wavefunctions_module,   ONLY : psic
      USE ions_base,              ONLY : nat, ntyp => nsp, ityp
      USE uspp_param,             ONLY : nh
      USE lsda_mod,               ONLY : current_spin
      USE uspp,                   ONLY : qq
      USE becmod,                 ONLY : bec_type, becp
      USE fft_base,               ONLY : dffts, dtgs
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, last
      !
      INTEGER :: ih, jh, ikb, jkb, nt, ia, ir, mbia
      REAL(DP) :: fac
      COMPLEX(DP) , ALLOCATABLE :: w1(:)
      !
      REAL(DP), EXTERNAL :: ddot
      !

      CALL start_clock( 's_psir' )
   
      IF( dtgs%have_task_groups ) CALL errore( 's_psir_k', 'task_groups not implemented', 1 )

      call set_xkphase(current_k)

      !
      fac = sqrt(omega)
      !
      ikb = 0
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
                     psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + xkphase(ir,ia)*betasave(ir,ih,ia)*w1(ih)
                     !
                  ENDDO
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
      CALL stop_clock( 's_psir' )
      !
      RETURN
      !
  END SUBROUTINE s_psir_k
  !
  SUBROUTINE add_vuspsir_gamma ( ibnd, last )
  !--------------------------------------------------------------------------
  !
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector transformed in real space contained in psic.
  !    ibnd is an index that runs up to band last. 
  !    Requires the products of psi with all beta functions
  !    in array becp%r(nkb,last) (calculated by calbecr in REAL SPACE)
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
  USE fft_base,               ONLY : dffts, dtgs
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, last
  !
  INTEGER :: ih, jh, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
  !
  REAL(DP), EXTERNAL :: ddot
  !
  CALL start_clock( 'add_vuspsir' )

  IF( dtgs%have_task_groups ) THEN

    CALL errore( 'add_vuspsir_gamma', 'task_groups not implemented', 1 )

  ELSE !non task groups part starts here

   !
   fac = sqrt(omega)
   !
   ikb = 0
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
                  IF ( ibnd+1 <= last )  w2(ih) = w2(ih) + deeq(ih,jh,ia,current_spin)* &
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
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + &
                       betasave(ir,ih,ia)*cmplx( w1(ih), w2(ih) ,kind=DP)
                  !
               ENDDO
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
  SUBROUTINE add_vuspsir_k ( ibnd, last )
  !--------------------------------------------------------------------------
  !
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector transformed in real space contained in psic.
  !    ibnd is an index that runs up to band las.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,last) (calculated by calbecr in REAL SPACE)
  ! Subroutine written by Stefano de Gironcoli, modified by O. Baris Malcioglu
  ! WARNING ! for the sake of speed, no checks performed in this subroutine
  !
  ! k-point phase factor fixed by SdG 030816
  !
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE cell_base,              ONLY : omega
  USE wavefunctions_module,   ONLY : psic
  USE ions_base,              ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,             ONLY : nh
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq
  USE becmod,                 ONLY : bec_type, becp
  USE fft_base,               ONLY : dffts, dtgs
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, last
  !
  INTEGER :: ih, jh, ikb, jkb, nt, ia, ir, mbia
  REAL(DP) :: fac
  !
  COMPLEX(DP), ALLOCATABLE :: w1(:)
  !
  REAL(DP), EXTERNAL :: ddot
  !
  CALL start_clock( 'add_vuspsir' )

  IF( dtgs%have_task_groups ) CALL errore( 'add_vuspsir_k', 'task_groups not implemented', 1 )

  call set_xkphase(current_k)
   !
   fac = sqrt(omega)
   !
   ikb = 0
   !
   DO nt = 1, ntyp
      !
      DO ia = 1, nat
         !
         IF ( ityp(ia) == nt ) THEN
            !
            mbia = maxbox_beta(ia)

            ALLOCATE( w1(nh(nt)))
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
                  psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + xkphase(ir,ia)*betasave(ir,ih,ia)*w1(ih)
                  !
               ENDDO
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
   CALL stop_clock( 'add_vuspsir' )
   RETURN
  !
  END SUBROUTINE add_vuspsir_k

  !--------------------------------------------------------------------------
  SUBROUTINE invfft_orbital_gamma (orbital, ibnd, last, conserved)
  !--------------------------------------------------------------------------
  !
  ! OBM 241008
  ! This driver subroutine transforms the given orbital using fft and puts the
  ! result in psic
  ! Warning! In order to be fast, no checks on the supplied data are performed!
  !
  ! orbital: the array of orbitals to be transformed
  ! ibnd: band index of the band currently being transformed
  ! last: index of the last band you want to transform (usually the total number 
  !       of bands but can be different in band parallelization)
  !
    USE wavefunctions_module, &
                       ONLY : psic
    USE gvecs,         ONLY : nls,nlsm,doublegrid
    USE klist,         ONLY : ngk, igk_k
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, dtgs
    USE fft_interfaces,ONLY : invfft

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being transformed
                           last   ! index of the last band that you want to transform

    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    LOGICAL, OPTIONAL :: conserved
    ! if this flag is true, the orbital is stored in temporary memory

    !Internal temporary variables
    INTEGER :: j, idx, ioff

    !Task groups

    !The new task group version based on vloc_psi
    !print *, "->Real space"
    CALL start_clock( 'invfft_orbital' )
    !
    IF( dtgs%have_task_groups ) THEN
        !

        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*dtgs%nogrp, 2

           IF( idx + ibnd - 1 < last ) THEN
              DO j = 1, ngk(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =      orbital(j,idx+ibnd-1) +&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd)
                 tg_psic(nlsm(igk_k(j,1))+ioff) =conjg(orbital(j,idx+ibnd-1) -&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == last ) THEN
              DO j = 1, ngk(1)
                 tg_psic(nls (igk_k(j,1))+ioff) =        orbital(j,idx+ibnd-1)
                 tg_psic(nlsm(igk_k(j,1))+ioff) = conjg( orbital(j,idx+ibnd-1))
              ENDDO
           ENDIF

           ioff = ioff + dtgs%tg_nnr

        ENDDO
        !
        !
        CALL invfft ('Wave', tg_psic, dffts, dtgs)
        !
        !
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
          IF (.not. allocated(tg_psic_temp)) ALLOCATE( tg_psic_temp( dtgs%tg_nnr * dtgs%nogrp ) )
          tg_psic_temp=tg_psic
         ENDIF
        ENDIF

    ELSE !Task groups not used
        !
        psic(:) = (0.d0, 0.d0)

        IF (ibnd < last) THEN
           ! two ffts at the same time
           DO j = 1, ngk(1)
              psic (nls (igk_k(j,1))) =       orbital(j, ibnd) + (0.0d0,1.d0)*orbital(j, ibnd+1)
              psic (nlsm(igk_k(j,1))) = conjg(orbital(j, ibnd) - (0.0d0,1.d0)*orbital(j, ibnd+1))
           ENDDO
        ELSE
           DO j = 1, ngk(1)
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

    CALL stop_clock( 'invfft_orbital' )

  END SUBROUTINE invfft_orbital_gamma
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fwfft_orbital_gamma (orbital, ibnd, last,conserved)
  !--------------------------------------------------------------------------
  !
  ! OBM 241008
  ! This driver subroutine -back- transforms the given orbital using fft using
  ! the already existent data in psic. 
  ! Warning! This subroutine does not reset the orbital, use carefully!
  ! Warning 2! In order to be fast, no checks on the supplied data are performed!
  ! Variables:
  !
  ! orbital: the array of orbitals to be transformed
  ! ibnd: band index of the band currently being transformed
  ! last: index of the last band you want to transform (usually the total number 
  !       of bands but can be different in band parallelization)
  !
    USE wavefunctions_module, &
                       ONLY : psic
    USE klist,         ONLY : ngk, igk_k
    USE gvecs,         ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, dtgs
    USE fft_interfaces,ONLY : fwfft
    USE mp_bands,      ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being transformed
                           last   ! index of the last band that you want to transform

    COMPLEX(DP),INTENT(out) :: orbital(:,:)

    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    !Internal temporary variables
    COMPLEX(DP) :: fp, fm
    INTEGER :: j, idx, ioff

    !Task groups
    !print *, "->fourier space"
    CALL start_clock( 'fwfft_orbital' )
    !New task_groups versions
    IF( dtgs%have_task_groups ) THEN
       !
        CALL fwfft ('Wave', tg_psic, dffts, dtgs )
        !
        ioff   = 0
        !
        DO idx = 1, 2*dtgs%nogrp, 2
           !
           IF( idx + ibnd - 1 < last ) THEN
              DO j = 1, ngk(1)
                 fp= ( tg_psic( nls(igk_k(j,1)) + ioff ) +  &
                      tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk_k(j,1)) + ioff ) -  &
                      tg_psic( nlsm(igk_k(j,1)) + ioff ) ) * 0.5d0
                 orbital (j, ibnd+idx-1) =  cmplx( dble(fp), aimag(fm),kind=DP)
                 orbital (j, ibnd+idx  ) =  cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == last ) THEN
              DO j = 1, ngk(1)
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


        IF (ibnd < last) THEN

           ! two ffts at the same time
           DO j = 1, ngk(1)
              fp = (psic (nls(igk_k(j,1))) + psic (nlsm(igk_k(j,1))))*0.5d0
              fm = (psic (nls(igk_k(j,1))) - psic (nlsm(igk_k(j,1))))*0.5d0
              orbital( j, ibnd)   = cmplx( dble(fp), aimag(fm),kind=DP)
              orbital( j, ibnd+1) = cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           DO j = 1, ngk(1)
              orbital(j, ibnd)   =  psic (nls(igk_k(j,1)))
           ENDDO
        ENDIF
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
           IF (allocated(psic_temp) ) DEALLOCATE(psic_temp)
         ENDIF
        ENDIF
    ENDIF
    !
    CALL stop_clock( 'fwfft_orbital' )

  END SUBROUTINE fwfft_orbital_gamma
  !
  !--------------------------------------------------------------------------
  SUBROUTINE invfft_orbital_k (orbital, ibnd, last, ik, conserved)
    !--------------------------------------------------------------------------
    !
    ! OBM 110908
    ! This subroutine transforms the given orbital using fft and puts the result
    ! in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    !
    ! orbital: the array of orbitals to be transformed
    ! ibnd: band index of the band currently being transformed
    ! last: index of the last band you want to transform (usually the total number 
    !       of bands but can be different in band parallelization)
    !
    !  current_k  variable  must contain the index of the desired kpoint
    !
    USE kinds,                    ONLY : DP
    USE wavefunctions_module,     ONLY : psic
    USE klist,                    ONLY : ngk, igk_k
    USE wvfct,                    ONLY : current_k
    USE gvecs,                    ONLY : nls, nlsm, doublegrid
    USE fft_base,                 ONLY : dffts, dtgs
    USE fft_interfaces,           ONLY : invfft

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being transformed
                           last   ! index of the last band that you want to transform

    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    INTEGER, OPTIONAL :: ik
    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    ! Internal variables
    INTEGER :: ioff, idx, ik_

    CALL start_clock( 'invfft_orbital' )

    ik_ = current_k ; if (present(ik)) ik_ = ik

    IF( dtgs%have_task_groups ) THEN
       !
       tg_psic = ( 0.D0, 0.D0 )
       ioff   = 0
       !
       DO idx = 1, dtgs%nogrp
          !
          IF( idx + ibnd - 1 <= last ) THEN
             !DO j = 1, size(orbital,1)
             tg_psic( nls( igk_k(:, ik_) ) + ioff ) = orbital(:,idx+ibnd-1)
             !END DO
          ENDIF

          ioff = ioff + dtgs%tg_nnr

       ENDDO
       !
       CALL invfft ('Wave', tg_psic, dffts, dtgs)
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (.not. allocated(tg_psic_temp)) &
                  &ALLOCATE( tg_psic_temp( dtgs%tg_nnr * dtgs%nogrp ) )
             tg_psic_temp=tg_psic
          ENDIF
       ENDIF
       !
    ELSE  !non task_groups version
       !
       psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
       !
       psic(nls(igk_k(1:ngk(ik_), ik_))) = orbital(1:ngk(ik_),ibnd)
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
    CALL stop_clock( 'invfft_orbital' )
  END SUBROUTINE invfft_orbital_k
  !--------------------------------------------------------------------------
  SUBROUTINE fwfft_orbital_k (orbital, ibnd, last, ik, conserved)
    !--------------------------------------------------------------------------
    !
    ! OBM 110908
    ! This subroutine transforms the given orbital using fft and puts the result
    ! in psic
    ! Warning! In order to be fast, no checks on the supplied data are performed!
    !
    ! orbital: the array of orbitals to be transformed
    ! ibnd: band index of the band currently being transformed
    ! last: index of the last band you want to transform (usually the total number 
    !       of bands but can be different in band parallelization)
    !
    !  current_k  variable  must contain the index of the desired kpoint
    !
    USE wavefunctions_module,     ONLY : psic
    USE klist,                    ONLY : ngk, igk_k
    USE wvfct,                    ONLY : current_k
    USE gvecs,                    ONLY : nls, nlsm, doublegrid
    USE kinds,                    ONLY : DP
    USE fft_base,                 ONLY : dffts, dtgs
    USE fft_interfaces,           ONLY : fwfft
    USE mp_bands,                 ONLY : me_bgrp

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being transformed
                           last   ! index of the last band that you want to transform

    COMPLEX(DP),INTENT(out) :: orbital(:,:)
    INTEGER, OPTIONAL :: ik
    LOGICAL, OPTIONAL :: conserved
    !if this flag is true, the orbital is stored in temporary memory

    ! Internal variables
    INTEGER :: ioff, idx, ik_

    CALL start_clock( 'fwfft_orbital' )

    ik_ = current_k ; if (present(ik)) ik_ = ik

    IF( dtgs%have_task_groups ) THEN
       !
       CALL fwfft ('Wave', tg_psic, dffts, dtgs)
       !
       ioff   = 0
       !
       DO idx = 1, dtgs%nogrp
          !
          IF( idx + ibnd - 1 <= last ) THEN
             orbital (:, ibnd+idx-1) = tg_psic( nls(igk_k(:,ik_)) + ioff )

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
       orbital(1:ngk(ik_),ibnd) = psic(nls(igk_k(1:ngk(ik_),ik_)))
       !
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (allocated(psic_temp) ) DEALLOCATE(psic_temp)
          ENDIF
       ENDIF
    ENDIF
    CALL stop_clock( 'fwfft_orbital' )

  END SUBROUTINE fwfft_orbital_k

  !--------------------------------------------------------------------------
  SUBROUTINE v_loc_psir (ibnd, last)
    !--------------------------------------------------------------------------
    ! Basically the same thing as v_loc but without implicit fft
    ! modified for real space implementation
    ! OBM 241008
    !
    USE wavefunctions_module, &
                       ONLY : psic
    USE gvecs,         ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, dtgs
    USE fft_parallel,  ONLY : tg_gather
    USE mp_bands,      ONLY : me_bgrp
    USE scf,           ONLY : vrs
    USE lsda_mod,      ONLY : current_spin

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being operated on
                           last   ! index of the last band that you want to operate on
    !Internal temporary variables
    INTEGER :: j
    !Task groups
    REAL(DP),    ALLOCATABLE :: tg_v(:)
    CALL start_clock( 'v_loc_psir' )

    IF( dtgs%have_task_groups ) THEN
        IF (ibnd == 1 ) THEN
          CALL tg_gather( dffts, dtgs, vrs(:,current_spin), tg_v )
          !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        ENDIF
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
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
  SUBROUTINE v_loc_psir_inplace (ibnd, last)
    !--------------------------------------------------------------------------
    ! The same thing as v_loc_psir but 
    ! - on input  psic contains the wavefunction
    ! - on output psic overwritten to contain v_loc_psir 
    ! Therefore must be the first term to be considered whn building hpsi
    ! SdG 290716
    !
    USE wavefunctions_module, &
                       ONLY : psic
    USE gvecs,         ONLY : nls,nlsm,doublegrid
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts, dtgs
    USE fft_parallel,  ONLY : tg_gather
    USE mp_bands,      ONLY : me_bgrp
    USE scf,           ONLY : vrs
    USE lsda_mod,      ONLY : current_spin

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd,& ! index of the band currently being operated on
                           last   ! index of the last band that you want to operate on
    !Internal temporary variables
    INTEGER :: j
    !Task groups
    REAL(DP),    ALLOCATABLE :: tg_v(:)
    CALL start_clock( 'v_loc_psir' )

    IF( dtgs%have_task_groups ) THEN
        IF (ibnd == 1 ) THEN
          CALL tg_gather( dffts, dtgs, vrs(:,current_spin), tg_v )
          !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        ENDIF
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_v(j) * tg_psic(j)
        ENDDO
        !
        DEALLOCATE( tg_v )
     ELSE
        !   product with the potential v on the smooth grid
        !
        DO j = 1, dffts%nnr
           psic (j) = vrs(j,current_spin) * psic(j)
        ENDDO
     ENDIF
  CALL stop_clock( 'v_loc_psir' )
  END SUBROUTINE v_loc_psir_inplace
    !--------------------------------------------------------------------------
  !
END MODULE realus
