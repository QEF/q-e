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
  !! Module Real Ultra-Soft.
  !
  !! Originally written by Antonio Suriano and Stefano de Gironcoli;
  !! modified by Carlo Sbraccia;
  !! modified by O. Baris Malcioglu (2008);
  !! modified by P. Umari and G. Stenuit (2009);
  !! cleanup, GWW-specific stuff moved out by P. Giannozzi (2015);
  !! computation of dQ/dtau_i needed for forces added by P. Giannozzi (2015);
  !! some comments about the way some routines act added by S. de Gironcoli (2015);
  !! extended to generic k by S. de Gironcoli (2016);
  !! addusstress_r added by S. de Gironcoli (2016);
  !! task group reorganization by S. de Gironcoli (2016);
  !! additional parallelization work (OpenMP and MPI) by S. de Gironcoli (2020).
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: boxrad(:)
  !! radius of Q boxes, does not depend on the grid
  !
  ! memory space for Beta function in real space
  REAL(DP),ALLOCATABLE :: boxrad_beta(:)     ! (ntypx) radius of the boxes
  INTEGER              :: boxtot             ! total dimension of the boxes
  INTEGER              :: bboxtot            ! total dimension of the broken boxes
  ! nat dependent array
  INTEGER, ALLOCATABLE :: box0(:),  box_s(:),  box_e(:) ! offset, start, and end index of a given box section
  INTEGER, ALLOCATABLE :: bbox0(:), bbox_s(:), bbox_e(:)! offset, start, and end index for broken boxes
  INTEGER, ALLOCATABLE :: maxbbox_beta(:)    ! number of points in a given broken-box atomic section
  ! grid dependent arrays
  INTEGER, ALLOCATABLE :: bbox_ind(:)        ! index for internal bbox redistribution
  COMPLEX(DP),ALLOCATABLE :: box_psic(:)     ! a box-friendly version of psic
  COMPLEX(DP),ALLOCATABLE :: bbox_aux(:)     ! an auxiliary space to temporarily store bbox stuff
  INTEGER, ALLOCATABLE :: box_beta(:)        ! (combined) list of points in the atomic boxes (in box)
  REAL(DP),ALLOCATABLE :: xyz_beta(:,:)      ! coordinates of the points in box_beta (in bbox)
  REAL(DP),ALLOCATABLE :: betasave(:,:)      ! beta funtions on the (combined) list of points (in bbox)
  !! beta function
  !!
  COMPLEX(DP), ALLOCATABLE :: xkphase(:)     ! (combined) phases for all atomic boxes (in bbox)
  !! kpoint-related phase factor around each atom
  INTEGER :: current_phase_kpoint=-1
  !! the kpoint index for which the xkphase is currently set.
  !! Negative initial value  means not set
  !
  !! MPI parallelization stuff
  INTEGER :: intra_bbox_comm, me_bbox, nproc_bbox
  !! MPI communicator used to break boxes of the beta functions
  INTEGER, ALLOCATABLE :: bbox_rank(:)
  !! collection of the ranks of the members of the bbox group in order of appearance
  INTEGER, ALLOCATABLE :: send_id(:), recv_id(:), datasize(:), break_at(:), break_ir0(:), offset(:)
  !! MPI send/recv info
  INTEGER, ALLOCATABLE :: send_requests(:)
  !! MPI list of active non-blocking send requests in box_to_bbox/bbox_to_box conversions;
  !! used to delay completion in order to hide communication time
  LOGICAL :: tbbox
  !! .true. if atomic boxes are broken and redistributed
  !
  ! ... General
  !
  LOGICAL :: real_space = .FALSE.
  !! if true perform calculations in real space
  INTEGER :: initialisation_level
  !! init_realspace_vars sets this to 3; qpointlist adds 5; betapointlist adds 7
  !! so the value should be 15 if the real space routine is initialised properly
  !
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  COMPLEX(DP), ALLOCATABLE :: psic_temp(:)
  !! Copy of psic
  COMPLEX(DP), ALLOCATABLE :: tg_psic_temp(:)
  !! Copy of tg_psic
  COMPLEX(DP), ALLOCATABLE :: tg_vrs(:)
  !! task groups linear V memory
  COMPLEX(DP), ALLOCATABLE :: psic_box_temp(:),tg_psic_box_temp(:)
  !
  TYPE realsp_augmentation
  !! Contains the augmentation functions and related quantities, in real space for
  !! a single atom.
     INTEGER :: maxbox = 0
     !!  number of R points in the augmentation sphere of this atom
     INTEGER,ALLOCATABLE  :: box(:)
     !! (maxbox) Index of point R in the global order of the R-space grid
     REAL(DP),ALLOCATABLE :: dist(:)
     !! (maxbox) distances |R-tau_i| (box is centered on atom tau_i)
     REAL(DP),ALLOCATABLE :: xyz(:,:)
     !! (3,maxbox) coordinates of R-tau_i
     REAL(DP),ALLOCATABLE :: qr(:,:)
     !! (maxbox,number_of_q_funcs) the Q functions sampled over R points
  END TYPE realsp_augmentation
  ! 
  TYPE(realsp_augmentation),POINTER :: tabp(:) => null()
  !! Augmentation functions on the RHO (HARD) grid for all atoms
  !
  !TYPE(realsp_augmentation),POINTER :: tabs(:) => null()
  ! !Augmentation functions on the SMOOTH grid for all atoms
  !
  TYPE(realsp_augmentation),POINTER :: tabxx(:) => null()
  !! Augmentation functions on the EXX grid for all atoms
  !
  PRIVATE
  ! variables for real-space Q, followed by routines
  PUBLIC :: tabp, tabxx, boxrad, realsp_augmentation
  PUBLIC :: generate_qpointlist, qpointlist, addusdens_r, newq_r, &
       addusforce_r, addusstress_r, real_space_dq, deallocate_realsp
  ! variables for real-space beta, followed by routines
  PUBLIC :: real_space, initialisation_level, tg_psic, betasave, box_beta
  PUBLIC :: betapointlist, init_realspace_vars, v_loc_psir, v_loc_psir_inplace
  PUBLIC :: bbox0, bbox_s, bbox_e, maxbbox_beta, box_psic
  PUBLIC :: invfft_orbital_gamma, fwfft_orbital_gamma, s_psir_gamma, &
            calbec_rs_gamma, add_vuspsir_gamma, invfft_orbital_k,    &
            fwfft_orbital_k, s_psir_k, calbec_rs_k, add_vuspsir_k, add_box_to_psic
  !
  CONTAINS
  
    !------------------------------------------------------------------------
    SUBROUTINE generate_qpointlist
      !------------------------------------------------------------------------
      USE fft_base,     ONLY : dfftp  !, dffts
      USE funct,        ONLY : dft_is_hybrid
      !USE gvecs,        ONLY : doublegrid
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
    !
    !----------------------------------------------------------------------------
    SUBROUTINE init_realspace_vars()
     !---------------------------------------------------------------------------
     !! This subroutine should be called to allocate/reset real space related
     !! variables.
     !
     USE fft_base,             ONLY : dffts

     IMPLICIT NONE

     !print *, "<<<<<init_realspace_vars>>>>>>>"

     !real space, allocation for task group fft work arrays

     IF( dffts%has_task_groups ) THEN
        !
        IF (allocated( tg_psic ) ) DEALLOCATE( tg_psic )
        !
        ALLOCATE( tg_psic( dffts%nnr_tg ) )
        ALLOCATE( tg_vrs ( dffts%nnr_tg ) )
        !
     ENDIF
     !
     initialisation_level = initialisation_level + 7

    END SUBROUTINE init_realspace_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_realsp()
      !------------------------------------------------------------------------
      !! Deallocate real space related variables.
      !
!      USE gvecs,      ONLY : doublegrid
      !
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
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_realsp_aug( tab )
      !------------------------------------------------------------------------
      !! Deallocate \(\text{realsp_augmentation}\) type variable.
      !
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
    SUBROUTINE qpointlist( dfft, tabp )
      !------------------------------------------------------------------------
      !! This subroutine computes and stores into \(\text{tabp}\) the values of \(Q(r)\)
      !! on atom-centered boxes in the real-space grid of descriptor "dfft".
      !! Must be called at startup and every time atoms (tau_i) move.
      !! A set of spherical boxes are computed for each atom. The Q functions
      !! \(Q(r-\tau_i)\) are then interpolated on the real-space grid points.
      !! On output:
      !
      !! * \(\text{boxrad}\) contains the radii of the boxes for each atom type;
      !! * \(\text{tabp(nat)}\) contains pointers to structure tabp, each contaning:
      !!    * \(\text{maxbox}\) = number of real-space grid points in the atom-centered box;
      !!    * \(\text{qr(maxbox)}\) = interpolated values of Q(r) on the points of the box;
      !!    * \(\text{box(maxbox)}\) = index of grid points in the real-space grid "dfft".
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi, fpi, eps16, eps6
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, alat
      USE uspp,       ONLY : okvan, qq_at, qq_nt, nhtol
      USE uspp_param, ONLY : upf, nh
      USE atom,       ONLY : rgrid
      USE fft_types,  ONLY : fft_type_descriptor
      USE mp_bands,   ONLY : intra_bgrp_comm
      USE mp,         ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      TYPE(fft_type_descriptor),INTENT(in)    :: dfft
      TYPE(realsp_augmentation),POINTER,INTENT(inout) :: tabp(:)
      !
      INTEGER               :: ia, mbia
      INTEGER               :: indm, ih, jh, ijh, il, jl
      INTEGER               :: roughestimate, nt
      INTEGER,  ALLOCATABLE :: buffpoints(:)
      INTEGER               :: ir, i, j, k, ijv
      INTEGER               :: imin, imax, ii, jmin, jmax, jj, kmin, kmax, kk
      REAL(DP)              :: distsq, posi(3)
      REAL(DP), ALLOCATABLE :: boxdist(:), xyz(:,:), diff(:)
      REAL(DP)              :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz, aux
      REAL(DP)              :: inv_nr1, inv_nr2, inv_nr3, boxradsq_ia, boxrad_ia
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
#if defined(__DEBUG)
      write (stdout,*) 'BOXRAD = ', boxrad(:) * alat
#endif
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
      inv_nr1 = 1.D0 / dble( dfft%nr1 )
      inv_nr2 = 1.D0 / dble( dfft%nr2 )
      inv_nr3 = 1.D0 / dble( dfft%nr3 )
      !
      ! the qq_at matrices are recalculated  here. initialize them
      qq_at(:,:,:) =0.d0

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
         boxrad_ia   = boxrad(nt)
         boxradsq_ia = boxrad_ia**2
         !
         mbia = 0
         !
         ! ... compute the needed ranges for i,j,k  indices around the atom position in crystal coordinates
         !
         posi(:) = tau(:,ia) ; CALL cryst_to_cart( 1, posi, bg, -1 )
         imin = NINT((posi(1) - boxrad_ia * sqrt(bg(1,1)*bg(1,1)+bg(2,1)*bg(2,1)+bg(3,1)*bg(3,1)))*dfft%nr1)
         imax = NINT((posi(1) + boxrad_ia * sqrt(bg(1,1)*bg(1,1)+bg(2,1)*bg(2,1)+bg(3,1)*bg(3,1)))*dfft%nr1)
         jmin = NINT((posi(2) - boxrad_ia * sqrt(bg(1,2)*bg(1,2)+bg(2,2)*bg(2,2)+bg(3,2)*bg(3,2)))*dfft%nr2)
         jmax = NINT((posi(2) + boxrad_ia * sqrt(bg(1,2)*bg(1,2)+bg(2,2)*bg(2,2)+bg(3,2)*bg(3,2)))*dfft%nr2)
         kmin = NINT((posi(3) - boxrad_ia * sqrt(bg(1,3)*bg(1,3)+bg(2,3)*bg(2,3)+bg(3,3)*bg(3,3)))*dfft%nr3)
         kmax = NINT((posi(3) + boxrad_ia * sqrt(bg(1,3)*bg(1,3)+bg(2,3)*bg(2,3)+bg(3,3)*bg(3,3)))*dfft%nr3)
#if defined(__DEBUG)
         write (*,*) tau(:,ia)
         write (*,*) posi
         write (*,*) imin,imax,jmin,jmax,kmin,kmax
#endif
         DO k = kmin, kmax
            kk = modulo(k,dfft%nr3) - dfft%my_i0r3p
            if (kk .LT. 0 .OR. kk .ge. dfft%my_nr3p ) cycle
            DO j = jmin, jmax
               jj = modulo(j,dfft%nr2) - dfft%my_i0r2p
               if (jj .LT. 0 .OR. jj .ge. dfft%my_nr2p ) cycle
               DO i = imin, imax
                  ii = modulo(i,dfft%nr1)
                  !
                  posi(:) = i * inv_nr1*at(:,1) + j * inv_nr2*at(:,2) + k * inv_nr3*at(:,3) - tau(:,ia)
                  !
                  distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
                  IF ( distsq < boxradsq_ia ) THEN
                     ! compute fft index ir from ii,jj,kk
                     ir = 1 + ii + jj * dfft%nr1x + kk * dfft%nr1x * dfft%my_nr2p
                     !
                     mbia = mbia + 1
                     IF( mbia > roughestimate ) CALL errore('qpointlist', 'rough-estimate is too rough', 3)
                     !
                     buffpoints(mbia)    = ir
                     boxdist(mbia)       = sqrt( distsq )*alat
                     xyz(:,mbia)         = posi(:)*alat
                     !
                  ENDIF
               END DO
            END DO
         END DO
#if defined(__DEBUG)
         WRITE (stdout,*) 'QPOINTLIST  ATOM ', ia, ' MBIA =', mbia
#endif
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
      ! collect the result of the qq_at matrix across processors
      CALL mp_sum( qq_at, intra_bgrp_comm )
      ! and test that they don't differ too much from the result computed on the atomic grid
      ALLOCATE (diff(nsp))
      diff(:)=0.0_dp
      do ia=1,nat
         nt = ityp(ia)
         ijh = 0
         do ih = 1, nh(nt)
            il = nhtol (ih, nt)
            do jh = ih, nh(nt)
               jl = nhtol (jh, nt)
               ijh=ijh+1
               if (abs(qq_at(ih,jh,ia)-qq_nt(ih,jh,nt)) .gt. eps6*10 ) then
                  diff(nt) = MAX(diff(nt), abs(qq_at(ih,jh,ia)-qq_nt(ih,jh,nt)))
               end if
            end do
         end do
      end do
      IF ( ANY( diff(:) > 0.0_dp ) ) THEN
         CALL infomsg('real_space_q', &
                      'integrated real space q too different from target q')
         WRITE(stdout, &
         '(5X,"Largest difference: ",f10.6," for atom type #",i2)') &
            MAXVAL(diff), MAXLOC(diff(:))
      END IF
      DEALLOCATE( diff )
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
      !! Compute \(Q_{lm}(r-\tau_{ia})\) centered around  atom \(\text{ia}\) of
      !! type \(\text{nt}\).
      !! First we compute the radial \(Q(r)\) (\(\text{qtot}\)) for each \(l\),
      !! then we interpolate it in our mesh (contained in \(\text{tabp}\)), finally
      !! we multiply by the spherical harmonics and store it into \(\text{tab}\).
      !
      USE constants,  ONLY : eps16, eps6
      USE uspp,       ONLY : indv, nhtolm, ap, qq_at
      USE uspp_param, ONLY : upf, lmaxq, nh
      USE atom,       ONLY : rgrid
      USE splinelib,  ONLY : spline, splint
      USE cell_base,  ONLY : omega
      USE fft_base,   ONLY : dfftp
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ia, nt, mbia
      TYPE(realsp_augmentation), INTENT(INOUT), POINTER :: tab(:)
      !
      INTEGER  :: l, nb, mb, ijv, lllnbnt, lllmbnt, ir, nfuncs, lm, i0, ijh, ih, jh
      REAL(dp) :: first, second, qtot_int
      REAL(dp), ALLOCATABLE :: qtot(:), dqtot(:), xsp(:), wsp(:), rl2(:), spher(:,:)

      CALL start_clock( 'realus:tabp' )

      if (mbia.ne.tab(ia)%maxbox) then
         write (stdout,*) 'real space q: ia,nt,mbia,tab(ia)%maxbox ', ia, nt, mbia, tab(ia)%maxbox
         call errore('real_space_q','inconsistent tab(ia)%maxbox dimension',ia)
      end if
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
      ! ... prevent FP errors if r(1) = 0  (may happen for some grids)
      !
      IF ( rgrid(nt)%r(1) < eps16 ) THEN
         i0 = 2
      ELSE
         i0 = 1
      END IF
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
                  qtot(i0:upf(nt)%kkbeta) = &
                       upf(nt)%qfuncl(i0:upf(nt)%kkbeta,ijv,l) &
                       / rgrid(nt)%r(i0:upf(nt)%kkbeta)**2
                  IF ( i0 == 2 ) qtot(1) = qtot(2)
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
               DO ir = 1, mbia
                  !
                  ! ... spline interpolation
                  !
                  qtot_int = splint( xsp, qtot, wsp, tab(ia)%dist(ir) )
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
      ! compute qq_at interating qr in the sphere
      ijh = 0
      DO ih = 1, nh(nt)
         DO jh = ih, nh(nt)
            ijh = ijh + 1
            qq_at(ih,jh,ia) = sum (tab(ia)%qr(1:mbia,ijh) )
            qq_at(jh,ih,ia) = qq_at(ih,jh,ia)
         END DO
      END DO
      qq_at(:,:,ia) = qq_at(:,:,ia) * omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

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
      !! Compute \(dQ(r-\tau_i)/d\tau_{ia}\) using the formula:
      !! \begin{equation}\notag
      !!   \frac{dQ(r-\tau_i)}{d\tau_{ia}} = - \frac{dq}{dr}(|r-\tau_i|)\frac{
      !!                                       (r_a-\tau_{ia})}{|r-\tau_i|}
      !!                                       Y_{lm}(r-\tau_i)
      !!                                     - q(|r-\tau_i|)\cdot
      !!                                       \frac{dY_{lm}(r-\tau_i)}{d\tau_{ia}}
      !! \end{equation}
      !! This routine follows the same logic of \(\texttt{real_space_q}\).
      !
      USE constants,  ONLY : eps8, eps16
      USE uspp,       ONLY : indv, nhtolm, ap
      USE uspp_param, ONLY : upf, lmaxq, nh
      USE atom,       ONLY : rgrid
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ia, nt, mbia, nfuncs
      REAL(dp),INTENT(OUT):: dqr(mbia,nfuncs,3)
      !
      INTEGER  :: l, nb, mb, ijv, lllnbnt, lllmbnt, ir, lm, ijh, ih, jh, ipol
      REAL(dp) :: first, second, qtot_int, dqtot_int, dqxyz(3)
      REAL(dp), ALLOCATABLE :: qtot(:), dqtot(:), xsp(:), wsp(:,:), &
           rl2(:), spher(:,:), dspher(:,:,:)
      !
      CALL start_clock( 'realus:tabp' )
      if (mbia.ne.tabp(ia)%maxbox) then
         write (stdout,*) 'real space dq: ia,nt,mbia,tabp(ia)%maxbox ', ia, nt, mbia, tabp(ia)%maxbox
         call errore('real_space_dq','inconsistent tab(ia)%maxbox dimension',ia)
      end if
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
                  if (rgrid(nt)%r(1)< eps16) then
!                      if (l==0) then
                         qtot(1) = qtot(2)
                         qtot(1) = qtot(3)
#if defined(__DEBUG)
                         write (stdout,*) qtot(2),qtot(3)
#endif
!                      else
!                         qtot(1) = 0.d0
!                      end if
                  end if
               ENDIF
               !
               ! ... compute the first derivative
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
                  ! ... spline interpolation
                  !
                  qtot_int = splint( xsp, qtot, wsp(:,1), tabp(ia)%dist(ir) )
                  dqtot_int= splint( xsp,dqtot, wsp(:,2), tabp(ia)%dist(ir) )
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
                                   ( qtot_int*dspher(ir,lm,ipol) + &
                                     dqxyz(ipol)*spher(ir,lm) ) * &
                                     ap(lm,nhtolm(ih,nt),nhtolm(jh,nt))
                           ENDDO
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
      !! This subroutine is the driver routine of the box system in this
      !! implementation of beta functions terms in real space.
      !! US case only for now but should not be.
      !
      !! * all the variables common in the module are computed and stored for
      !!   reusing;
      !! * this routine has to be called every time the atoms are moved and of
      !!   course at the beginning;
      !! * a set of spherical boxes are computed for each atom;
      !! * in \(\text{boxrad_beta}\) there are the radii of the boxes;
      !! * in \(\text{maxbox_beta}\) the upper limit of leading index, namely the
      !!   number of points of the fine mesh contained in each box;
      !! * in \(\text{xyz_beta}\) there are the coordinates of the points with
      !!   origin in the centre of atom;
      !
      ! ... Most of time is spent here; the calling routines are faster.
      !
      USE constants,  ONLY : pi
      USE control_flags, ONLY : gamma_only
      USE ions_base,  ONLY : nat, nsp, ityp, tau
      USE cell_base,  ONLY : at, bg, alat
      USE uspp,       ONLY : okvan, indv,nhtolm
      USE uspp_param, ONLY : upf, lmaxq, nh, nhm
      USE atom,       ONLY : rgrid
      USE fft_base,   ONLY : dffts
      USE splinelib,  ONLY : spline, splint
      !
      IMPLICIT NONE
      !
      LOGICAL, PARAMETER   :: tprint = .false.  ! whether to print info on betapointlist
      !
      INTEGER              :: ia, it, mbia
      INTEGER              :: indm, inbrx, ih, kkbeta, kkbeta_max
      INTEGER              :: roughestimate, lamx2, nt
      INTEGER, ALLOCATABLE :: maxbox_beta(:)  ! number of points in a given box atomic section
      INTEGER, ALLOCATABLE :: tmp_box_beta(:,:)
      REAL(DP),ALLOCATABLE :: tmp_xyz_beta(:,:,:), tmp_boxdist_beta(:,:)
      REAL(DP),ALLOCATABLE :: spher_beta(:,:), boxdist_beta(:)
      REAL(DP)             :: distsq, qtot_int, first, second
      INTEGER              :: ir, i, j, k, lm, nb, box_ir
      INTEGER              :: imin, imax, ii, jmin, jmax, jj, kmin, kmax, kk
      REAL(DP)             :: posi(3)
      REAL(DP),ALLOCATABLE :: xsp(:), ysp(:), wsp(:), d1y(:), d2y(:)
      REAL(DP)             :: mbr, mbx, mby, mbz, dmbx, dmby, dmbz
      REAL(DP)             :: inv_nr1s, inv_nr2s, inv_nr3s, tau_ia(3), boxradsq_ia, boxrad_ia

      integer              :: nsend ! number of non-blocking sends left behind by box_to_bbox
      !
      initialisation_level = initialisation_level + 5
      IF ( .not. okvan ) CALL errore &
                        ('betapointlist','real space routines for USPP only',1)
      !
      !print *, "<<<betapointlist>>>"
      !
      CALL start_clock( 'betapointlist' )
      !
      ! ... betasave and other arrays from previous step are deallocated here for a fresh start
      !
      IF ( allocated( betasave ) ) DEALLOCATE( betasave )
      IF ( allocated( xyz_beta ) ) DEALLOCATE( xyz_beta )
      IF ( allocated( box_beta ) ) DEALLOCATE( box_beta )
      IF ( allocated( bbox_aux ) ) DEALLOCATE( bbox_aux )
      IF ( allocated( box_psic ) ) DEALLOCATE( box_psic )
      IF ( allocated( bbox_ind ) ) DEALLOCATE( bbox_ind )
      IF ( allocated( xkphase ) )  DEALLOCATE( xkphase )
      !
      !  ... only the first time the routine is being called
      IF ( .not. allocated( boxrad_beta ) ) THEN
         ! ... we calculate the radius of each spherical box ( one
         ! ... for each non-local projector )
         ALLOCATE( boxrad_beta( nsp ) ) ; boxrad_beta(:) = 0.D0
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
         boxrad_beta(:) = boxrad_beta(:) / alat
      ENDIF
      !  ... and allocate the arrays that only depend on nat
      IF (.not. ALLOCATED( box0  ) ) ALLOCATE ( box0  ( nat ) )
      IF (.not. ALLOCATED( box_s ) ) ALLOCATE ( box_s ( nat ) )
      IF (.not. ALLOCATED( box_e ) ) ALLOCATE ( box_e ( nat ) )
      IF (.not. ALLOCATED( bbox0  ) ) ALLOCATE ( bbox0  ( nat ) )
      IF (.not. ALLOCATED( bbox_s ) ) ALLOCATE ( bbox_s ( nat ) )
      IF (.not. ALLOCATED( bbox_e ) ) ALLOCATE ( bbox_e ( nat ) )
      IF (.not. ALLOCATED( maxbbox_beta) ) ALLOCATE( maxbbox_beta( nat ) )
      ALLOCATE( maxbox_beta( nat ) ) ; maxbox_beta(:) = 0 ! this will be deallocated at the end
      !
      ! ... a rough estimate for the number of grid-points per box
      ! ... is computed here
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
      ALLOCATE( tmp_box_beta( roughestimate, nat ) )     ; tmp_box_beta(:,:) = 0
      ALLOCATE( tmp_xyz_beta( 3, roughestimate, nat ) )  ; tmp_xyz_beta(:,:,:) = 0.D0
      ALLOCATE( tmp_boxdist_beta( roughestimate, nat ) ) ; tmp_boxdist_beta(:,:) = 0.D0
      !
      ! The beta functions are treated on smooth grid
      !
      inv_nr1s = 1.D0 / dble( dffts%nr1 )
      inv_nr2s = 1.D0 / dble( dffts%nr2 )
      inv_nr3s = 1.D0 / dble( dffts%nr3 )
      !
      DO ia = 1, nat
         !
         IF ( .not. upf(ityp(ia))%tvanp ) CYCLE
         !
         boxrad_ia = boxrad_beta(ityp(ia))
         boxradsq_ia = boxrad_ia**2
         !
         tau_ia(1) = tau(1,ia)
         tau_ia(2) = tau(2,ia)
         tau_ia(3) = tau(3,ia)
         !
         ! ... compute the needed ranges for i,j,k  indices around the atom position in crystal coordinates
         !
         posi(:) = tau_ia(:) ; CALL cryst_to_cart( 1, posi, bg, -1 )
         imin = NINT((posi(1) - boxrad_ia * sqrt(bg(1,1)*bg(1,1)+bg(2,1)*bg(2,1)+bg(3,1)*bg(3,1)))*dffts%nr1)
         imax = NINT((posi(1) + boxrad_ia * sqrt(bg(1,1)*bg(1,1)+bg(2,1)*bg(2,1)+bg(3,1)*bg(3,1)))*dffts%nr1)
         jmin = NINT((posi(2) - boxrad_ia * sqrt(bg(1,2)*bg(1,2)+bg(2,2)*bg(2,2)+bg(3,2)*bg(3,2)))*dffts%nr2)
         jmax = NINT((posi(2) + boxrad_ia * sqrt(bg(1,2)*bg(1,2)+bg(2,2)*bg(2,2)+bg(3,2)*bg(3,2)))*dffts%nr2)
         kmin = NINT((posi(3) - boxrad_ia * sqrt(bg(1,3)*bg(1,3)+bg(2,3)*bg(2,3)+bg(3,3)*bg(3,3)))*dffts%nr3)
         kmax = NINT((posi(3) + boxrad_ia * sqrt(bg(1,3)*bg(1,3)+bg(2,3)*bg(2,3)+bg(3,3)*bg(3,3)))*dffts%nr3)
         if (tprint) then
            write (*,*) tau_ia
            write (*,*) posi
            write (*,*) imin,imax,jmin,jmax,kmin,kmax
         end if
         DO k = kmin, kmax
            kk = modulo(k,dffts%nr3) - dffts%my_i0r3p
            if (kk .LT. 0 .OR. kk .ge. dffts%my_nr3p ) cycle
            DO j = jmin, jmax
               jj = modulo(j,dffts%nr2) - dffts%my_i0r2p
               if (jj .LT. 0 .OR. jj .ge. dffts%my_nr2p ) cycle
               DO i = imin, imax
                  ii = modulo(i,dffts%nr1)
                  !
                  posi(:) = i * inv_nr1s*at(:,1) + j * inv_nr2s*at(:,2) + k * inv_nr3s*at(:,3) - tau_ia(:)
                  !
                  distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
                  IF ( distsq < boxradsq_ia ) THEN
                     ! compute fft index ir from ii,jj,kk
                     ir = 1 + ii + jj * dffts%nr1x + kk * dffts%nr1x * dffts%my_nr2p
                     !
                     mbia = maxbox_beta(ia) + 1
                     maxbox_beta(ia) = mbia
                     !
                     tmp_box_beta(mbia,ia)     = ir
                     tmp_boxdist_beta(mbia,ia) = distsq * alat**2 ! the square is first needed in ylm2
                     tmp_xyz_beta(:,mbia,ia)   = posi(:)*alat
                     !
                  ENDIF
               END DO
            END DO
         END DO
         if (tprint) WRITE (*,*) 'BETAPOINTLIST: ATOM ',ia, ' MAXBOX_BETA =', maxbox_beta(ia)
      ENDDO
      !
      boxtot = sum(maxbox_beta(1:nat)) ; WRITE( stdout,* )  'BOXTOT', boxtot
      box0(1) = 0 ; box_s(1) = 1 ; box_e(1) = maxbox_beta(1)
      do ia =2,nat
         box0(ia) = box_e(ia-1) ; box_s(ia) = box0(ia) + 1 ;  box_e(ia) = box0(ia) + maxbox_beta(ia)
      end do
      !
      ! ... deal with the unbalance of box dimensions by breaking them
      !
      CALL beta_box_breaking (nat, maxbox_beta, maxbbox_beta, tbbox)
      !
      bboxtot = sum(maxbbox_beta(1:nat)) ; WRITE( stdout,* )  'BBOXTOT', bboxtot
      bbox0(1) = 0 ; bbox_s(1) = 1 ; bbox_e(1) = maxbbox_beta(1)
      do ia =2,nat
         bbox0(ia) = bbox_e(ia-1) ; bbox_s(ia) = bbox0(ia)+1 ;  bbox_e(ia) = bbox0(ia) + maxbbox_beta(ia)
      end do
      !
      ! ... now store grid point in a more convenient place (sometime in box sometime in bbox) ...
      !
      ALLOCATE( box_beta     ( boxtot ) )            ! the box_psic -> psic index (in box)
      ALLOCATE( xyz_beta  ( 3, bboxtot ) )           ! the relative position (in bbox)
      ALLOCATE( boxdist_beta ( bboxtot ) )           ! and its length squared (in bbox)
      ALLOCATE( xkphase      ( bboxtot ) )           ! the phase (in bbox)
      ALLOCATE( box_psic( max( bboxtot, boxtot ) ) ) ! used to communicate between box and bbox
      !
      ! ... fix the bbox turn-around needed for box_to_bbox and bbox_to_box operations
      !
      call beta_box_turning(nat, maxbox_beta, maxbbox_beta, tbbox)
      !
      do ia =1,nat
         ! store  box_psic -> psic index (used in add_box_to_psic)
         box_beta ( box_s(ia):box_e(ia) ) =  tmp_box_beta(  1:maxbox_beta(ia),ia)
      end do
      ! convert positions and distances from box to bbox and store them (uses box_psic as vector)
      !-xy components
      do ia =1,nat
         box_psic ( box_s(ia):box_e(ia) ) = CMPLX( tmp_xyz_beta( 1, 1:maxbox_beta(ia), ia), &
                                                   tmp_xyz_beta( 2, 1:maxbox_beta(ia), ia)  )
      end do
      if (tbbox) call box_to_bbox ( nsend )
      do ia =1,nat
         xyz_beta ( 1, bbox_s(ia):bbox_e(ia) ) = dble ( box_psic ( bbox_s(ia):bbox_e(ia) ) )
         xyz_beta ( 2, bbox_s(ia):bbox_e(ia) ) = aimag( box_psic ( bbox_s(ia):bbox_e(ia) ) )
      end do
      if (tbbox) call wait_bbox_send_to_complete ( nsend )
      !-z component & dist^2
      do ia =1,nat
         box_psic ( box_s(ia):box_e(ia) ) = CMPLX( tmp_xyz_beta( 3, 1:maxbox_beta(ia), ia), &
                                                   tmp_boxdist_beta( 1:maxbox_beta(ia), ia) )
      end do
      if (tbbox) call box_to_bbox ( nsend )
      do ia =1,nat
         xyz_beta ( 3, bbox_s(ia):bbox_e(ia) ) = dble ( box_psic ( bbox_s(ia):bbox_e(ia) ) )
         boxdist_beta( bbox_s(ia):bbox_e(ia) ) = aimag( box_psic ( bbox_s(ia):bbox_e(ia) ) )
      end do
      if (tbbox) call wait_bbox_send_to_complete ( nsend )
      !
      call set_xkphase(1)
      !
      DEALLOCATE( tmp_xyz_beta,  tmp_box_beta, tmp_boxdist_beta )
      DEALLOCATE( maxbox_beta )
      !
      CALL stop_clock( 'realus:boxes' )
      CALL start_clock( 'realus:spher' )
      !
      ! ... now it computes the spherical harmonics (in bbox)
      !
      lamx2 = lmaxq*lmaxq
      !
      ALLOCATE( spher_beta( bboxtot, lamx2 ) ) ! ( boxtot,lmax2 )
      !
      CALL ylmr2( lamx2, bboxtot, xyz_beta, boxdist_beta, spher_beta )
      !
      if (gamma_only) DEALLOCATE( xyz_beta )  ! not needed anymore (no further set_xkphase call)
      !
      CALL stop_clock( 'realus:spher' )
      CALL start_clock( 'realus:tabp' )
      !
      ! ... let's do the main work
      !
      ! ... set boxdist_beta to be the distance (needed in the spline interpolation)
      !$omp parallel do
      DO box_ir = 1, bboxtot
         boxdist_beta( box_ir ) = sqrt ( boxdist_beta (box_ir ) )
      ENDDO
      !$omp end parallel do

      ALLOCATE( betasave( bboxtot, nhm )  ) ! (in bbox)
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
      ! ... allocate variables used for spline interpolation
      !
      kkbeta_max = maxval ( upf(1:nsp)%kkbeta )
      ALLOCATE( xsp( kkbeta_max ), ysp( kkbeta_max ), wsp( kkbeta_max ) )
      ALLOCATE( d1y( kkbeta_max ), d2y( kkbeta_max ) )

      DO ia = 1, nat
         !
         mbia = maxbbox_beta(ia); IF ( mbia == 0 ) CYCLE
         !
         nt = ityp(ia)
         IF ( .not. upf(nt)%tvanp ) CYCLE
         kkbeta = upf(nt)%kkbeta
         !
         ! ... the radii in x
         !
         xsp(1:kkbeta) = rgrid(nt)%r(1:kkbeta)
         !
         DO ih = 1, nh (nt)
            !
            lm = nhtolm(ih, nt)
            nb = indv(ih, nt)
            !
            !OBM rgrid(nt)%r(1) == 0, attempting correction
            ! In the UPF file format, beta field is r*|beta>
            IF (rgrid(nt)%r(1)==0) THEN
             ysp(2:kkbeta) = upf(nt)%beta(2:kkbeta,nb) / rgrid(nt)%r(2:kkbeta)
             ysp(1)=0.d0
            ELSE
             ysp(1:kkbeta) = upf(nt)%beta(1:kkbeta,nb) / rgrid(nt)%r(1:kkbeta)
            ENDIF

            CALL radial_gradient(ysp(1:kkbeta), d1y, rgrid(nt)%r, kkbeta, 1)
            CALL radial_gradient(d1y, d2y, rgrid(nt)%r, kkbeta, 1)

            first = d1y(1) ! first derivative in first point
            second =d2y(1) ! second derivative in first point

            CALL spline( xsp, ysp, first, second, wsp, kkbeta ) ! NB: optional kkbeta is needed

            DO box_ir = bbox_s(ia), bbox_e(ia)
               ! ... spline interpolate the value of f_l(r) in point ir in atom ia
               qtot_int = splint( xsp,ysp,wsp, boxdist_beta(box_ir), kkbeta ) ! NB: optional kkbeta is needed
               ! ... spher_beta is the Y_lm in point ir for atom ia
               betasave(box_ir,ih) = qtot_int * spher_beta(box_ir,lm) 
               !
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE( d1y, d2y )
      DEALLOCATE( xsp, ysp, wsp )
      !
      DEALLOCATE( spher_beta, boxdist_beta )
      !
      CALL stop_clock( 'realus:tabp' )
      CALL stop_clock( 'betapointlist' )
      !
    END SUBROUTINE betapointlist
    !
    SUBROUTINE beta_box_breaking (nat, maxbox_beta, maxbbox_beta, tbbox)
      !! This routine determines whether it is needed to re-shuffle the
      !! beta point grids in real space.
      USE mp_bands,   ONLY : me_bgrp, intra_bgrp_comm, nproc_bgrp
      USE mp,         ONLY : mp_sum, mp_comm_split, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      ! I/O variables
      !
      INTEGER, INTENT(IN)  :: nat                 ! number of atom
      INTEGER, INTENT(IN)  :: maxbox_beta(nat)    ! box dimension of each atom
      INTEGER, INTENT(OUT) :: maxbbox_beta(nat)   ! (broken) box dimension of each atom
      LOGICAL, INTENT(OUT) :: tbbox               ! .TRUE. if boxes are broken
      !
      ! local variables
      !
      REAL(DP)             :: boxtot_avg, unbalance, opt_unbalance, maximum_nn
      INTEGER              :: ia, ib, in, nn, ip, np, id, my_id, boxtot_max, opt_nn, my_color
      INTEGER, ALLOCATABLE :: ind(:), ind0(:), color(:), data(:), boxtot_beta(:)
      INTEGER              :: recv, send, delta, delta_tot
      INTEGER              :: target_ndata, recv_target, send_target
      !
      INTEGER, ALLOCATABLE :: mbox(:,:)            ! working array for box_to_bbox and back
      !
      REAL(DP), PARAMETER  :: tollerable_unbalance = 1.25D0
      INTEGER, PARAMETER   :: maximum_nn_param = 6
      !
      ALLOCATE ( boxtot_beta ( nproc_bgrp ) )
      boxtot_beta (:) = 0 ; boxtot_beta ( me_bgrp+1 ) = SUM(maxbox_beta(1:nat))
      CALL mp_sum(  boxtot_beta, intra_bgrp_comm )

      boxtot_avg = SUM ( boxtot_beta ( 1:nproc_bgrp ) ) / dble ( nproc_bgrp )
      boxtot_max = MAXVAL( boxtot_beta ( 1:nproc_bgrp ) )
      unbalance  =  boxtot_max / boxtot_avg
!      WRITE ( stdout ,* ) ' me_bgrp ', me_bgrp ; FLUSH(6)
      WRITE ( stdout, * ) 'BETAPOINTLIST SUMMARY: ', boxtot_max , boxtot_avg, unbalance ; FLUSH(6)
      WRITE ( stdout, * ) 'BETATOT LIST:', boxtot_beta(:) ; FLUSH(6)

      nn = 1; opt_nn = 1 ; opt_unbalance = unbalance

      ALLOCATE ( ind(nproc_bgrp), color(nproc_bgrp),  data(nproc_bgrp) )
      data(:) = - boxtot_beta(:) ; ind(1) = 0
      call ihpsort ( nproc_bgrp, data, ind )

      maximum_nn = min ( maximum_nn_param, nproc_bgrp )
      DO WHILE ( unbalance > tollerable_unbalance .and. nn < maximum_nn )

         nn = nn + 1

         np = nproc_bgrp / nn ; IF ( nproc_bgrp .ne. np * nn  ) CYCLE

         color(ind(1:np)) = ind(1:np)
         DO in = nn, 2, -1
            DO ip = 1, np
               data(ip) = data(ip) + data (in*np+1-ip)
               color(ind(in*np+1-ip)) = ind(ip)
            END DO
            call ihpsort ( np, data, ind )
         END DO
         unbalance =  -data(1)/(boxtot_avg*nn)
!         WRITE ( stdout, * ) nn, ' is a factor ',-data(1), -data(1)/boxtot_avg, unbalance ; FLUSH(6)
         IF ( unbalance < opt_unbalance ) THEN
            opt_unbalance = unbalance ; opt_nn = nn ; my_color = color(me_bgrp+1)
         END IF
         DO ip = 1, np
            data(ip) = -boxtot_beta(ind(ip))
         END DO

      END DO
      nn = opt_nn
      WRITE ( stdout, * ) ' opt_nn is', opt_nn, 'with opt_unbalance', opt_unbalance ; FLUSH(6)

      IF ( nn > 1 ) THEN

! collect optimal colors for reporting and redistribution
         color(:) = 0 ; color(me_bgrp+1) = my_color
         CALL mp_sum( color, intra_bgrp_comm )

         np = nproc_bgrp / nn

         !- report
         DO ip = 1, np
            WRITE ( stdout, '(A,i4,A,$)' ) ' color ', color(ind(ip)), ':'  ; FLUSH(6)
            target_ndata = 0 ; id = 0
            DO in = 1, nproc_bgrp
               IF ( color(in) == color(ind(ip)) ) THEN
                  WRITE (stdout, '(i4,$)' ) in ; FLUSH(6)
                  id = id + 1
                  target_ndata = target_ndata + boxtot_beta(in)
               END IF
            END DO
            WRITE ( stdout, * ) ' ', target_ndata ; FLUSH(6)
            IF ( id .ne. nn ) CALL errore( 'beta_box_breaking', '# of terms .ne. nn', 1)
         END DO
         !- create new communicators (labeled by my_color) inside intra_bgrp_comm
         CALL mp_comm_split( intra_bgrp_comm, my_color, me_bgrp, intra_bbox_comm )
         me_bbox =  mp_rank( intra_bbox_comm ); WRITE ( stdout, * ) 'me_bbox', me_bbox ; FLUSH(6)
         nproc_bbox = mp_size( intra_bbox_comm ); WRITE ( stdout, * ) ' nproc_bbox', nproc_bbox ; FLUSH(6)
         if ( nproc_bbox .ne. nn ) call errore( 'beta_box_breaking', 'nproc_bbox .ne. nn', nn)

         DEALLOCATE ( data, ind ) ; ALLOCATE ( data ( nn ), ind( nn ), ind0 ( nn ) )

         !- compute target dimensions of the broken boxes of my_color
         !- define an ordered identifier (id) in the group and set my_id
         id = 0 ; target_ndata = 0 ; data = 0 ; ind = 0
         DO in = 1, nproc_bgrp
            IF ( color(in) == my_color ) THEN
               id       = id + 1
               ind(id)  = in
               data(id) = boxtot_beta ( in )
               target_ndata = target_ndata + data(id)
               IF ( in == me_bgrp + 1 ) my_id = id
            END IF
         END DO
         IF ( id .ne. nn ) CALL errore( 'beta_box_breaking', '# of terms .ne. nn', 2)

         !- link the order in the group to the mp_rank in the communicator
         ALLOCATE ( bbox_rank( nproc_bbox ) ) ; bbox_rank=0; bbox_rank(my_id) = me_bbox
         call mp_sum( bbox_rank, intra_bbox_comm )

         !- define the target dimension of the broken boxes
         recv_target = 1 + ( target_ndata - 1 )/ nproc_bbox
         send_target =       target_ndata / nproc_bbox

         WRITE ( stdout, * ) 'targets ', target_ndata, recv_target, send_target ; FLUSH(6)
         !
         ALLOCATE ( send_id(nproc_bbox-1), recv_id(nproc_bbox-1), datasize(nproc_bbox-1) )
         ALLOCATE ( break_at(nproc_bbox-1), break_ir0(nproc_bbox-1), offset(nproc_bbox-1) )
         ALLOCATE ( send_requests(300) ) ! hopefully a safe allocation. should be fixed

         ! collected atomic box dimensions to plan grid point redistribution
         ALLOCATE ( mbox(nat,nproc_bbox) )
         mbox = 0 ; mbox(:,my_id) = maxbox_beta(:); CALL mp_sum( mbox, intra_bbox_comm )

         write ( stdout, * ) ind(1:nproc_bbox) ; FLUSH(6)
         do ia = 1, nat
            write ( stdout, * ) (mbox(ia,id),id=1,nproc_bbox), '|',sum(mbox(ia,1:nproc_bbox)) ; FLUSH(6)
         end do
         write ( stdout, * ) "--------------------------------------------------------------"
         write ( stdout, * ) (sum(mbox(1:nat,id)),id=1,nproc_bbox) ; FLUSH(6)
         write ( stdout, * ) ind(1:nproc_bbox) ; FLUSH(6)
         write ( stdout, * ) (sum(mbox(1:nat,id)),id=1,nproc_bbox) ; FLUSH(6)
         write ( stdout, * ) (data(id),id=1,nproc_bbox) ; FLUSH(6)

         !- compute dimensions of the broken boxes of my_color and the needed exchanges
         delta_tot = 0 ; ind0 = 0
         do ip = 1, nn - 1   ! NB: nn = nproc_bbox
            CALL ihpsort ( nn, data, ind0 ) ! the first time ind0 is initialized, after it is reshuffled
!            do id =1, nn     ! write the bbox group members (bgrp and bbox ordering) and their data dim
!               WRITE ( stdout, * ) ind(ind0(id)), bbox_rank(ind0(id)),':', data(id) ; FLUSH(6)
!            end do
            send = data(nn) - send_target
            recv = recv_target - data(1)
            delta = min ( recv, send )
            if ( data(nn) == recv_target ) delta = 0 ! avoid having a last point going around (see NB later)
!            write( stdout, * )'recv, send, delta', recv, send, delta  ; FLUSH(6)
            delta_tot = delta_tot + delta
            data(1)  = data(1)  + delta
            data(nn) = data(nn) - delta
            write ( stdout, * ) ind(ind0(nn)),bbox_rank(ind0(nn)),' -> ', ind(ind0(1)),bbox_rank(ind0(1)), delta
            ! set the communication plan
            send_id (ip) = ind0(nn) ! who sends
            recv_id (ip) = ind0(1)  ! who receives
            datasize(ip) = delta    ! transmitted data size
            offset(ip)   = data(nn) ! the global offset of the breaking point
            ! determine the the effect on the atomic data distribution (up to a permutation)
            ! NB: as who gains never loses later (delta=0 fix!) permutations are immaterial in the count
            do ia =nat,1,-1
               delta = delta - mbox(ia,ind0(nn))
               mbox(ia,ind0(1))  = mbox(ia,ind0(1)) + mbox(ia,ind0(nn))
               mbox(ia,ind0(nn)) = 0
               if ( delta > 0 ) cycle
               break_at (ip) = ia         ! the atomic index of the breaking point
               break_ir0(ip) = -delta     ! the offset of the breaking point
               mbox(ia,ind0(1))  = mbox(ia,ind0(1)) - break_ir0(ip)
               mbox(ia,ind0(nn)) = break_ir0(ip)
               write ( stdout, * ) '     break_at, break_ir0:', break_at(ip), break_ir0(ip); FLUSH(6)
               write ( stdout, * ) ind(1:nproc_bbox); FLUSH(6)
               do ib = 1, nat
                  write ( stdout, * ) (mbox(ib,id),id=1,nproc_bbox), '|',sum(mbox(ib,1:nproc_bbox)); FLUSH(6)
               end do
               write ( stdout, * ) "--------------------------------------------------------------"
               write ( stdout, * ) (sum(mbox(1:nat,id)),id=1,nproc_bbox); FLUSH(6)
               write ( stdout, * ) ind(ind0(1:nproc_bbox)); FLUSH(6)
               write ( stdout, * ) (sum(mbox(1:nat,ind0(id))),id=1,nproc_bbox); FLUSH(6)
               write ( stdout, * ) (data(id),id=1,nproc_bbox); FLUSH(6)
               exit
            end do
         end do
         do id =1, nproc_bbox
            do ib =1, nproc_bbox
               if (ind0(ib) == id ) WRITE ( stdout, * ) ind(id),':', data(ib); FLUSH(6)
            end do
         end do
         WRITE ( stdout, * ) 'total number of points moved ', delta_tot; FLUSH(6)
         maxbbox_beta(:) = mbox(:, my_id )  ! these are the atomic dimensions after redistribution

         DEALLOCATE ( mbox )
         DEALLOCATE ( ind0 )

         tbbox = .TRUE.

      ELSE

         intra_bbox_comm = 0
         maxbbox_beta(:) = maxbox_beta(:)
         tbbox = .FALSE.

      END IF

      DEALLOCATE ( data, color, ind, boxtot_beta )

      RETURN
      !
    END SUBROUTINE beta_box_breaking
    !
    SUBROUTINE beta_box_turning(nat, maxbox_beta, maxbbox_beta, tbbox )
      !! This routine complete the box/bbox redistribution by findig the permutation needed to
      !! pack the grid points in bbox in contiguous lists for each atom
      !
      IMPLICIT NONE
      !
      ! I/O variables
      !
      INTEGER, INTENT(IN) :: nat                 ! number of atom
      INTEGER, INTENT(IN) :: maxbox_beta(nat)    ! box dimension of each atom
      INTEGER, INTENT(IN) :: maxbbox_beta(nat)   ! (broken) box dimension of each atom
      LOGICAL, INTENT(IN) :: tbbox               ! .TRUE. if boxes are broken
      !
      ! local variables
      !
      INTEGER :: ia, box_ir, capel
      INTEGER, ALLOCATABLE :: grid_check(:)

      !$omp parallel
      CALL threaded_barrier_memset(box_psic, 0.D0, max(boxtot,bboxtot)*2)
      DO ia = 1, nat
         !$omp do
         DO box_ir = box_s(ia), box_e(ia)
            box_psic ( box_ir ) =  ia
         END DO
         !$omp end do
      END DO
      !$omp end parallel

      ALLOCATE ( grid_check(nat) ); grid_check(1:nat) = 0
      do box_ir =  1, boxtot
         ia = nint( dble (box_psic ( box_ir ) ) )
         grid_check(ia) = grid_check(ia) + 1
      end do
      capel = MAXVAL (ABS ( grid_check(1:nat)-maxbox_beta(1:nat) ) )
      IF ( capel > 0 ) THEN
         write (stdout,*) 'box atomic grids'
         write (stdout,*) (grid_check(ia),ia=1,nat)
         write (stdout,*) (maxbox_beta(ia),ia=1,nat)
         call errore ('beta_box_turning',' wrong box grid number ', capel )
      END IF

      if ( ALLOCATED ( bbox_ind ) ) DEALLOCATE( bbox_ind )
      if ( ALLOCATED ( bbox_aux ) ) DEALLOCATE( bbox_aux )

      if (tbbox) call box_to_bbox ()

      ALLOCATE ( bbox_ind ( bboxtot ), bbox_aux ( bboxtot ) )
      grid_check(1:nat) = 0 ; bbox_is_ordered = .true.
      WRITE(stdout,*) 'BBOXTOT', bboxtot
      do box_ir =  1, bboxtot
         ia = nint( dble (box_psic (box_ir ) ) )
         grid_check(ia) = grid_check(ia) + 1
         bbox_ind ( box_ir ) = bbox0(ia)+grid_check(ia)
      end do
      capel = MAXVAL (ABS ( grid_check(1:nat)-maxbbox_beta(1:nat) ) )
      IF ( capel > 0 ) THEN
         write (stdout,*) 'bbox atomic grids'
         write (stdout,*) (grid_check(ia),ia=1,nat)
         write (stdout,*) (maxbbox_beta(ia),ia=1,nat)
         call errore ('beta_box_turning',' wrong bbox grid number ', capel )
      END IF

      DEALLOCATE ( grid_check )

      RETURN

    END SUBROUTINE beta_box_turning
    !
    !------------------------------------------------------------------------
    SUBROUTINE box_to_bbox ( nsend_ )
      !-----------------------------------------------------------------------
      !! This routine redistributes box_psic inplace from the box to the bbox ordering
      !! a final index turn-around is neeed to bring same atom points in contiguous positions
      !
      USE mp,         ONLY : mp_send, mp_recv, mp_barrier
      !
      IMPLICIT NONE
      !
      INTEGER, OPTIONAL, INTENT(OUT) :: nsend_
      INTEGER :: ip, id_send, id_recv, ndata, ir0, msg_tag, my_size, box_ir, nsend
      LOGICAL :: sending, recving
      !
      if ( present(nsend_) ) nsend_ = 0 ! initialize nsend_ if present
      if ( .not. tbbox ) return

!      write ( stdout, * ) 'enter  box_to_bbox '; FLUSH(6)
!      write ( stdout, * ) 'me_bbox = ', me_bbox; FLUSH(6)
      call mp_barrier ( intra_bbox_comm )

      my_size = boxtot

      nsend = 0
      do ip =1, nproc_bbox -1

         msg_tag = ip

         id_send = send_id(ip) ; sending = ( me_bbox == bbox_rank(send_id(ip)) )
         id_recv = recv_id(ip) ; recving = ( me_bbox == bbox_rank(recv_id(ip)) )
         ndata = datasize(ip)
         ir0   = offset(ip)

         !NB this is the intended direct operation so variable names mean what they say
!         write (stdout,*) 'msg_tag is',msg_tag; FLUSH(6)
!         write (stdout,*) 'sending is',id_send, bbox_rank(id_send); FLUSH(6)
!         write (stdout,*) 'recving is',id_recv, bbox_rank(id_recv); FLUSH(6)
!         write (stdout,*) 'ndata is', datasize(ip); FLUSH(6)
!         write (stdout,*) 'offset is', ir0; FLUSH(6)
!         write (stdout,*) 'sending is', sending,' recving is', recving ; FLUSH(6)

         if ( ndata == 0 ) cycle

!recv should be blocking as box_psic is accesssed immediately after
         if (recving) then ! id_recv is receiving ndata from id_send
!            write ( stdout, * )  id_recv, ' receives ',ndata,' from ', id_send ; FLUSH(6)
            call mp_recv ( box_psic(my_size+1:my_size+ndata), bbox_rank(id_send), msg_tag, intra_bbox_comm )
            my_size = my_size + ndata
         end if
!send could be non-blocking so that computation can start with the unmodified part of box_psic
         if (sending) then ! id_send is sending ndata to id_recv
!            write ( stdout, * )  id_send, ' sends ',ndata,' to ', id_recv ; FLUSH(6)
            my_size = my_size - ndata
            if (ir0 .ne. my_size) call errore('box_to_bbox',' wrong size ', ip)
            call mp_send ( box_psic(my_size+1:my_size+ndata), bbox_rank(id_recv), msg_tag, intra_bbox_comm, &
                                                                                       send_requests, nsend )
         end if

!         write(stdout,*) 'done', my_size ; FLUSH(6)

      end do

      if (present(nsend_)) then ! return nsend_ for a delayed completion of the non-blocking sends
         nsend_ = nsend
      else ! wait here for the completion of the non-blocking sends (if any)
         call wait_bbox_send_to_complete ( nsend )
      end if

!      call mp_barrier ( intra_bbox_comm )

      if (my_size .ne. bboxtot) call errore('box_to_bbox',' wrong size ', nproc_bbox )

!      write ( stdout, * ) 'exit  box_to_bbox '; FLUSH(6)

      if (.not. ALLOCATED( bbox_ind ) )  RETURN

      CALL threaded_memcpy( bbox_aux, box_psic, 2 * bboxtot )
      !$omp parallel do
      do box_ir =1, bboxtot
         box_psic ( bbox_ind( box_ir ) ) = bbox_aux ( box_ir )
      end do
      !$omp end parallel do
      !
      RETURN
    END SUBROUTINE box_to_bbox
    !
    !------------------------------------------------------------------------
    SUBROUTINE wait_bbox_send_to_complete ( nsend )
      !-----------------------------------------------------------------------
      !! wait for the non blocking send in box_to_bbox to complete
      !! so that box_psic can be freely modified again
      !
      USE mp,         ONLY : mp_waitall
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nsend

      if (nsend > 0) call mp_waitall( send_requests(1:nsend) )

      RETURN
    END SUBROUTINE wait_bbox_send_to_complete
    !
    !------------------------------------------------------------------------
    SUBROUTINE bbox_to_box ( nsend )
      !-----------------------------------------------------------------------
      !! This routine redistributes box_psic inplace from the box to the bbox ordering
      !! an initial index turn-around brings atom points in the needed send/recv locations
      !
      USE mp,         ONLY : mp_send, mp_recv, mp_barrier
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(OUT) :: nsend
      !
      INTEGER :: ip, id_send, id_recv, ndata, ir0, msg_tag, my_size, box_ir
      LOGICAL :: sending, recving
      !
      nsend = 0 ! initialize nsend
      if ( .not. tbbox ) return

!      write ( stdout, * ) 'enter  bbox_to_box '; FLUSH(6)
!      write ( stdout, * ) 'me_bbox = ', me_bbox; FLUSH(6)

      if (.not. ALLOCATED( bbox_ind ) ) CALL errore ( 'bbox_to_box',' bbox_ind not yet initialized',1)
!
      !$omp parallel do
      do box_ir =1, bboxtot
         bbox_aux ( box_ir ) = box_psic ( bbox_ind( box_ir ) )
      end do
      !$omp end parallel do
      CALL threaded_memcpy( box_psic, bbox_aux, 2 * bboxtot )

      call mp_barrier ( intra_bbox_comm )

      my_size = bboxtot

      do ip =nproc_bbox -1, 1, -1

         msg_tag = ip

         id_send = send_id(ip) ; sending = ( me_bbox == bbox_rank(send_id(ip)) )
         id_recv = recv_id(ip) ; recving = ( me_bbox == bbox_rank(recv_id(ip)) )
         ndata = datasize(ip)
         ir0   = offset(ip)

         !NB this is the inverse of box_to_bbox so variable names mean their opposite
!         write (stdout,*) 'msg_tag is',msg_tag; FLUSH(6)
!         write (stdout,*) 'un-sending is',id_send, bbox_rank(id_send); FLUSH(6)
!         write (stdout,*) 'un-recving is',id_recv, bbox_rank(id_recv); FLUSH(6)
!         write (stdout,*) 'ndata is', datasize(ip); FLUSH(6)
!         write (stdout,*) 'offset is', ir0; FLUSH(6)
!         write (stdout,*) 'un-sending is', sending,' un-recving is', recving ; FLUSH(6)

         if ( ndata == 0 ) cycle

!send could be non-blocking so it can start before receive proc is ready and data are immediately ready.
         if (recving) then ! id_recv is actually sending ndata to id_send
!            write ( stdout, * )  id_recv, ' sends ',ndata,' to ', id_send ; FLUSH(6)
            my_size = my_size - ndata
            call mp_send ( box_psic(my_size+1:my_size+ndata), bbox_rank(id_send), msg_tag, intra_bbox_comm, &
                                                                                       send_requests, nsend )
         end if
!receive must be blocking as data are needed immediately after and it should be the bottleneck
         if (sending) then ! id_send is actually receiving ndata from id_recv
!            write ( stdout, * )  id_send, ' receives ',ndata,' from ', id_recv ; FLUSH(6)
            if (ir0 .ne. my_size) call errore('bbox_to_box',' wrong size ', ip)
            call mp_recv ( box_psic(my_size+1:my_size+ndata), bbox_rank(id_recv), msg_tag, intra_bbox_comm )
            my_size = my_size + ndata
         end if

!         write(stdout,*) 'done',my_size ; FLUSH(6)

      end do

      call mp_barrier ( intra_bbox_comm )

      if (my_size .ne. boxtot) call errore('bbox_to_box',' wrong size ', nproc_bbox )

!      write ( stdout, * ) 'exit  bbox_to_box '; FLUSH(6)
      !
      RETURN
    END SUBROUTINE bbox_to_box
    !
    !------------------------------------------------------------------------
    SUBROUTINE newq_r( vr,deeq, skip_vltot )
      !-----------------------------------------------------------------------
      !! This routine computes the integral of the perturbed potential with
      !! the \(Q\) function in real space.
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
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: vr(dfftp%nnr,nspin)
      !! The input potential
      REAL(DP), INTENT(OUT) :: deeq(nhm,nhm,nat,nspin)
      !! Contribution to the integral
      LOGICAL, INTENT(in) :: skip_vltot
      !! If .FALSE. \(\text{vltot}\) is added to \(\text{vr}\) when necessary
      !
      ! ... local variables
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      !
      INTEGER :: ia, ih, jh, is, ir, nt
      INTEGER :: mbia
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
    SUBROUTINE addusdens_r( rho )
      !------------------------------------------------------------------------
      !! This routine adds to the charge density the part which is due to
      !! the US augmentation, in real space.
      !
      USE ions_base,        ONLY : nat, ityp
      USE uspp,             ONLY : okvan, becsum
      USE uspp_param,       ONLY : upf, nh
      USE noncollin_module, ONLY : nspin_mag
      USE fft_interfaces,   ONLY : fwfft
      USE fft_base,         ONLY : dfftp
      USE wavefunctions,    ONLY : psic
#if defined (__DEBUG)
      USE noncollin_module, ONLY : nspin_lsda
      USE constants,        ONLY : eps6
      USE klist,            ONLY : nelec
      USE cell_base,        ONLY : omega
      USE mp_bands,         ONLY : intra_bgrp_comm
      USE mp,               ONLY : mp_sum
      USE gvect,            ONLY : gstart
#endif
      !
      IMPLICIT NONE
      ! The charge density to be augmented (in G-space)
      COMPLEX(kind=dp), INTENT(inout) :: rho(dfftp%ngm,nspin_mag)
      !
      INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, mbia
      REAL(kind=dp), ALLOCATABLE :: rhor(:,:)
#if defined (__DEBUG)
      CHARACTER(len=80) :: msg
      REAL(kind=dp) :: charge
#endif
      !
      !
      IF ( .not. okvan ) RETURN
      !
      CALL start_clock( 'addusdens' )
      !
      ALLOCATE ( rhor(dfftp%nnr,nspin_mag) )
      rhor(:,:) = 0.0_dp
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
                     rhor(irb,is) = rhor(irb,is) + tabp(ia)%qr(ir,ijh)*becsum(ijh,ia,is)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !
      ENDDO
      !
      !
      DO is = 1, nspin_mag
         psic(:) = rhor(:,is)
         CALL fwfft ('Rho', psic, dfftp)
         rho(:,is) = rho(:,is) + psic(dfftp%nl(:))
      END DO
      !
      DEALLOCATE ( rhor )
#if defined (__DEBUG)
      !
      ! ... check the total charge (must not be summed on k-points)
      !
      IF ( gstart == 2) THEN
         charge = SUM(rho(1,1:nspin_lsda) )*omega
      ELSE
         charge = 0.0_dp
      ENDIF
      CALL mp_sum(  charge , intra_bgrp_comm )
      write (stdout,*) 'charge before rescaling ', charge
      IF ( abs(charge - nelec) > MAX(eps6,eps6*ABS(nelec)) ) THEN
         !
         ! ... the error on the charge is too large, stop and complain
         !
         WRITE (msg,'("expected ",f10.6,", found ",f10.6)') nelec, charge
         CALL errore( 'addusdens_r', 'WRONG CHARGE '//trim(msg)//&
              ': ions may be overlapping or increase ecutrho', 1 )
      ENDIF
      !
#endif
      CALL stop_clock( 'addusdens' )
      !
      RETURN
      !
    END SUBROUTINE addusdens_r
    !
    !----------------------------------------------------------------------
    SUBROUTINE addusforce_r( forcenl )
      !----------------------------------------------------------------------
      !! This routine computes the contribution to atomic forces due
      !! to the dependence of the \(Q\) function on the atomic position.
      !! \begin{equation}\notag
      !! \begin{split}
      !!    F_{j,at} = &-\int \sum_{lm} \frac{dQ_{lm}(r)}{d\tau_{at}}
      !!                 \text{becsum}(lm,at)\ V(r)\ dr \\
      !!               &+\int \sum_{lm} \frac{dQ_{lm}(r)}{d\tau_{at}}
      !!                 \text{ebecsum}(lm,at)\ dr 
      !! \end{split}
      !! \end{equation}
      !! where:
      !! \begin{equation}\notag
      !! \begin{split}
      !!    &\text{becsum}(lm,at) = \sum_i \langle\psi_i|\beta_l\rangle w_i
      !!                                   \langle\beta_m|\psi_i\rangle, \\
      !!    &\text{ebecsum}(lm,at) = \sum_i\langle\psi_i|\beta_l\rangle w_i
      !!                                   \langle\beta_m|\psi_i\rangle et_i
      !! \end{split}
      !! \end{equation}
      !! On output the contribution is added to \(\text{forcenl}\).
      !
      USE kinds,      ONLY : DP
      USE cell_base,  ONLY : omega
      USE ions_base,  ONLY : nat, ityp
      USE fft_base,   ONLY : dfftp
      USE noncollin_module,   ONLY : nspin_mag
      USE scf,        ONLY : v, vltot
      USE uspp,       ONLY : becsum, ebecsum, okvan
      USE uspp_param, ONLY : upf,  nh
      USE mp_bands,   ONLY : intra_bgrp_comm
      USE mp,         ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT) :: forcenl (3, nat)
      !
      INTEGER :: na, nt, ijh, nfuncs, mbia, ir, is, irb
      REAL(dp), ALLOCATABLE:: dqr(:,:,:), forceq(:,:)
      REAL(dp) :: dqrforce(3), dqb(3), dqeb(3), v_eff
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
         dqrforce(:) = 0.0_dp
         !
         mbia = tabp(na)%maxbox
         IF ( mbia == 0 ) CYCLE
!         write (stdout,*) ' inside addusforce na, mbia ', na,mbia
         nfuncs = nh(nt)*(nh(nt)+1)/2
         ALLOCATE ( dqr(mbia,nfuncs,3) )
         CALL real_space_dq( nt, na, mbia, nfuncs, dqr )
         !
         DO ir = 1, mbia
            irb = tabp(na)%box(ir)

            DO is = 1, nspin_mag
               dqb = 0.d0; dqeb = 0.d0
               do ijh =1, nfuncs
                  dqb (:) = dqb (:) +  dqr(ir,ijh,:) *  becsum(ijh,na,is)
                  dqeb(:) = dqeb(:) +  dqr(ir,ijh,:) * ebecsum(ijh,na,is) 
               end do
               v_eff = vltot(irb) + v%of_r(irb,is) ;  IF (nspin_mag==4.and.is/=1) v_eff = v%of_r(irb,is)
               dqrforce(:) = dqrforce(:) + dqb(:) * v_eff - dqeb(:)
            ENDDO
 
         ENDDO
         DEALLOCATE ( dqr )
         !
         forceq(:,na) = - dqrforce(:) * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      ENDDO
      CALL mp_sum ( forceq, intra_bgrp_comm )
      !
      forcenl(:,:) = forcenl(:,:) + forceq(:,:)
      !
      DEALLOCATE (forceq )
      !
      RETURN
    END SUBROUTINE addusforce_r
    !
    !----------------------------------------------------------------------
    SUBROUTINE addusstress_r( sigmanl )
      !---------------------------------------------------------------------
      !! This routine computes the contribution to atomic stress due
      !! to the dependence of the Q function on the atomic position.
      !! \begin{equation}\notag
      !! \begin{split}
      !!    S_{i,j} = &+\int \sum_{lm} \left[r_i \frac{dQ_{lm}(r)}{d\tau_{at_j}} +
      !!                 Q_{lm}(r)\right] \text{becsum}(lm,at)\ V(ir)\ dr \\
      !!              &-\int \sum_{lm} \left[r_i \frac{dQ_{lm}(r)}{d\tau_{at_j}} +
      !!                 Q_{lm}(r)\right] \text{ebecsum}(lm,at)\ dr
      !! \end{split}
      !! \end{equation}
      !! where:
      !! \begin{equation}\notag
      !! \begin{split}
      !!    &\text{becsum}(lm,at) = \sum_i \langle\psi_i|\beta_l\rangle w_i
      !!                                   \langle\beta_m|\psi_i\rangle, \\
      !!    &\text{ebecsum}(lm,at) = \sum_i\langle\psi_i|\beta_l\rangle w_i
      !!                                   \langle\beta_m|\psi_i\rangle et_i
      !! \end{split}
      !! \end{equation}
      !! On output the contribution is added to \(\text{stressnl}\).
      !! NB: a factor -1.0/omega is later applied to \(\text{stressnl}\) as a whole.
      !
      ! FOR REASONS I DON'T UNDERSTAND THE SECOND TERM APPEARS TO NEED THE + SIGN TOO
      !
      USE kinds,      ONLY : DP
      USE cell_base,  ONLY : omega
      USE ions_base,  ONLY : nat, ityp
      USE fft_base,   ONLY : dfftp
      USE noncollin_module,  ONLY : nspin_mag
      USE scf,        ONLY : v, vltot
      USE uspp,       ONLY : becsum, ebecsum, okvan
      USE uspp_param, ONLY : upf,  nh
!      USE mp_bands,   ONLY : intra_bgrp_comm
!      USE mp,         ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT) :: sigmanl (3,3)
      !
      INTEGER :: ipol, na, nt, nfuncs, mbia, ir, is, irb , ijh
      REAL(dp), ALLOCATABLE:: dqr(:,:,:)
      REAL(dp) :: sus(3,3), sus_at(3,3), qb, qeb, dqb(3), dqeb(3), v_eff
      !
      IF (.not.okvan) RETURN
      !
      sus(:,:) = 0.0_dp
      !
      DO na = 1, nat
         nt = ityp(na)
         IF ( .NOT. upf(nt)%tvanp ) CYCLE
         !
         mbia = tabp(na)%maxbox
!         write (stdout,*) ' inside addusstress na, mbia ', na,mbia
         nfuncs = nh(nt)*(nh(nt)+1)/2
         ALLOCATE ( dqr(mbia,nfuncs,3) )
         CALL real_space_dq( nt, na, mbia, nfuncs, dqr )
         !
         sus_at(:,:) = 0.0_dp
         DO ir = 1, mbia
            irb = tabp(na)%box(ir)
            DO is = 1, nspin_mag
               qb = 0.d0; qeb =0.d0; dqb = 0.d0; dqeb = 0.d0
               do ijh =1, nfuncs
                  qb      = qb      +  tabp(na)%qr(ir,ijh) *  becsum(ijh,na,is)
                  qeb     = qeb     +  tabp(na)%qr(ir,ijh) * ebecsum(ijh,na,is)
                  dqb (:) = dqb (:) +  dqr(ir,ijh,:) *  becsum(ijh,na,is)
                  dqeb(:) = dqeb(:) +  dqr(ir,ijh,:) * ebecsum(ijh,na,is)
               end do
               v_eff = vltot(irb) + v%of_r(irb,is) ;  IF (nspin_mag==4.and.is/=1) v_eff = v%of_r(irb,is)
               do ipol=1, 3
                  sus_at(:, ipol) = sus_at(:,ipol) - tabp(na)%xyz(:,ir) * ( dqb(ipol) * v_eff + dqeb(ipol) )
                  sus_at(ipol,ipol) = sus_at(ipol,ipol) +                 ( qb * v_eff + qeb )
               end do
            END DO
         ENDDO
         DEALLOCATE ( dqr )

         sus (:,:) = sus(:,:) + sus_at(:,:)

      ENDDO
      !
!      CALL mp_sum ( sus, intra_bgrp_comm )
      sus(:,:) = sus(:,:) * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      sigmanl(:,:) = sigmanl(:,:) + sus(:,:)
      !
      RETURN
    END SUBROUTINE addusstress_r
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_xkphase( ik )
    !--------------------------------------------------------------------------
    !! In the calculation of \(\text{becp}\) or when performing \(\texttt{add_vuspsir}\)
    !! the wavefunction \(\text{psi_k}\) and not its periodic part (which is what
    !! we get from the FFT) should be used.
    !! A k-dependent phase \(\text{exp}[-xk(\text{current_k}\cdot(r-\tau_\text{ia}))]\) must
    !! be added.
    !
    USE kinds,      ONLY : DP
    USE klist,      ONLY : xk
    USE cell_base,  ONLY : tpiba

    IMPLICIT NONE
  
    INTEGER, INTENT (IN) :: ik

    INTEGER :: box_ir
    REAL(DP) :: arg

    if (.not.allocated ( xkphase ) ) call errore ('set_xkphase',' array not allocated yes',1)
    if (ik .eq. current_phase_kpoint ) return
    !
    !$omp parallel do
    do box_ir =1, bboxtot
       arg = ( xk(1,ik) * xyz_beta(1,box_ir) + &
               xk(2,ik) * xyz_beta(2,box_ir) + &
               xk(3,ik) * xyz_beta(3,box_ir) ) * tpiba
       xkphase( box_ir ) = CMPLX(COS(arg),-SIN(arg),KIND=dp)
    end do
    !$omp end parallel do
    !
    current_phase_kpoint = ik
    !
    return
    END SUBROUTINE set_xkphase

    !--------------------------------------------------------------------------
    SUBROUTINE calbec_rs_gamma( ibnd, last, becp_r )
      !--------------------------------------------------------------------------
      !! Calculates \(\text{becp_r}\) in real space.
      !
      !! * Requires \(\text{betasave}\), the beta functions at real space calculated
      !!   by \(\texttt{betapointlist()}\);
      !! * \(\text{ibnd}\) is an index that runs over the number of bands, which
      !!   is given by m, so you have to call this subroutine inside a cycle with
      !!   index \(\text{ibnd}\);
      !! * In this cycle you have to perform a Fourier transform of the orbital
      !!   corresponding to \(\text{ibnd}\), namely you have to transform the orbital
      !!   to real space and store it in the global variable \(\text{psic}\);
      !! * Remember that in the gamma_only case you perform two FFT at the same time,
      !!   so you have that the real part corresponds to \(\text{ibnd}\), and the
      !!   imaginary part to \(\text{ibnd}+1\).
      !
      !! WARNING: For the sake of speed, there are no checks performed in this
      !!          routine, check beforehand!
      !
      !! Subroutine written by Dario Rocca, Stefano de Gironcoli, modified by O.
      !! Baris Malcioglu.
      !! Some speedup and restructuring by SdG April 2020.
      !!
      !
    USE kinds,                 ONLY : DP
    USE cell_base,             ONLY : omega
    USE wavefunctions,         ONLY : psic
    USE ions_base,             ONLY : nat, nsp, ityp
    USE uspp_param,            ONLY : nh
    USE fft_base,              ONLY : dffts
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    USE uspp,                  ONLY : indv_ijkb0
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, last
    INTEGER :: ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr, wi
    REAL(DP) :: bcr, bci
    REAL(DP), DIMENSION(:,:), INTENT(out) :: becp_r
    integer :: ir, box_ir, maxbox, ijkb0, nh_nt, nsend
    !
    REAL(DP), EXTERNAL :: ddot
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( dffts%has_task_groups ) CALL errore( 'calbec_rs_gamma', 'task_groups not implemented', 1 )

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    maxbox = MAXVAL(maxbbox_beta(1:nat))
    !
    becp_r(:,ibnd)=0.d0
    IF ( ibnd+1 <= last ) becp_r(:,ibnd+1)=0.d0
    ! Clearly for an odd number of bands for ibnd=nbnd=last you don't have
    ! anymore bands, and so the imaginary part equal zero
    !
!
! copy psic into a box-friendly array (and scatter it across intra_bbox_comm if needed)
!
    !$omp parallel do
    DO box_ir =1, boxtot
       box_psic ( box_ir ) = psic( box_beta( box_ir ) )
    END DO
    !$omp end parallel do

    if (tbbox) call box_to_bbox ( nsend )

    ALLOCATE( wr(maxbox), wi(maxbox) )
    ! working arrays to order the points in the clever way
    DO nt = 1, nsp
       !
       nh_nt = nh(nt)
       !
       DO ia = 1, nat
          !
          IF ( ityp(ia) == nt ) THEN
             !
             mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE
             !
             ijkb0 = indv_ijkb0(ia)
             !$omp parallel default(shared) private(ih,ikb,ir,bcr,bci)
             !$omp do 
             DO ir =1, mbia
                wr(ir) = dble ( box_psic( bbox0(ia)+ir) )
             END DO
             !$omp end do
             !$omp do
             DO ih = 1, nh_nt
                !
                ikb = ijkb0 + ih
                bcr = ddot( mbia, betasave(bbox_s(ia):bbox_e(ia),ih), 1, wr(:) , 1 )
                becp_r(ikb,ibnd)   = fac * bcr
                !
             ENDDO
             !$omp end do nowait
             IF ( ibnd+1 <= last ) THEN
                !$omp do
                DO ir =1, mbia
                   wi(ir) = aimag( box_psic( bbox0(ia)+ir) )
                END DO
                !$omp end do
                !$omp do
                DO ih = 1, nh_nt
                   !
                   ikb = ijkb0 + ih
                   bci = ddot( mbia, betasave(bbox_s(ia):bbox_e(ia),ih), 1, wi(:) , 1 )
                   becp_r(ikb,ibnd+1) = fac * bci
                   !
                ENDDO
                !$omp end do
             END IF
             !$omp end parallel
             !
          ENDIF
          !
       ENDDO
       !
    ENDDO
    DEALLOCATE( wr, wi )

    if (tbbox) call wait_bbox_send_to_complete ( nsend )
    !
    CALL mp_sum( becp_r( :, ibnd ), intra_bgrp_comm )
    IF ( ibnd+1 <= last ) CALL mp_sum( becp_r( :, ibnd+1 ), intra_bgrp_comm )
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN

  END SUBROUTINE calbec_rs_gamma
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE calbec_rs_k( ibnd, last )
    !-------------------------------------------------------------------------------
    !! The k-point generalised version of \(\texttt{calbec_rs_gamma}\). Basically
    !! same as above, but \(\text{becp}\) is used instead of \(\text{becp_r}\),
    !! skipping the gamma point reduction derived from above by OBM 051108.
    !
    !! k-point phase factor fixed by SdG 030816.
    !! Some speedup and restrucuring by SdG April 2020.
    !
    USE kinds,                 ONLY : DP
    USE wvfct,                 ONLY : current_k
    USE cell_base,             ONLY : omega
    USE wavefunctions,         ONLY : psic
    USE ions_base,             ONLY : nat, nsp, ityp
    USE uspp_param,            ONLY : nh
    USE becmod,                ONLY : bec_type, becp
    USE fft_base,              ONLY : dffts
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    USE uspp,                  ONLY : indv_ijkb0
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd, last
    INTEGER :: ikb, nt, ia, ih, mbia
    REAL(DP) :: fac
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr, wi
    REAL(DP) :: bcr, bci
    integer :: ir, box_ir, maxbox, ijkb0, nh_nt, nsend
    !
    REAL(DP), EXTERNAL :: ddot
    !
    CALL start_clock( 'calbec_rs' )
    !
    IF( dffts%has_task_groups ) CALL errore( 'calbec_rs_k', 'task_groups not implemented', 1 )

    call set_xkphase(current_k)

    fac = sqrt(omega) / (dffts%nr1*dffts%nr2*dffts%nr3)
    !
    maxbox = MAXVAL(maxbbox_beta(1:nat))
    !
    becp%k(:,ibnd)=0.d0
!
! copy psic into a box-friendly array (and scatter it across intra_bbox_comm if needed)
!
    !$omp parallel do
    DO box_ir =1, boxtot
       box_psic ( box_ir ) = psic( box_beta( box_ir ) )
    END DO
    !$omp end parallel do

    if (tbbox) call box_to_bbox ( nsend )

    ALLOCATE( wr(maxbox), wi(maxbox) )
    ! working arrays to order the points in the clever way
    DO nt = 1, nsp
       !
       nh_nt = nh(nt)
       !
       DO ia = 1, nat
          !
          IF ( ityp(ia) == nt ) THEN
             !
             mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE
             !
             ijkb0 = indv_ijkb0(ia)
             !
             !$omp parallel default(shared) private(ih,ikb,ir,bcr,bci)
             !$omp do
             DO ir =1, mbia
                wr(ir) = dble ( box_psic( bbox0(ia)+ir ) * CONJG(xkphase(bbox0(ia)+ir)))
                wi(ir) = aimag( box_psic( bbox0(ia)+ir ) * CONJG(xkphase(bbox0(ia)+ir)))
             END DO
             !$omp end do
             !$omp do
             DO ih = 1, nh_nt
                ikb = ijkb0 + ih
                bcr = ddot( mbia, betasave(bbox_s(ia):bbox_e(ia),ih), 1, wr(:) , 1 )
                bci = ddot( mbia, betasave(bbox_s(ia):bbox_e(ia),ih), 1, wi(:) , 1 )
                becp%k(ikb,ibnd)   = fac * cmplx( bcr, bci, kind=DP )
             ENDDO
             !$omp end do
             !$omp end parallel
             !
          ENDIF
          !
       ENDDO
       !
    ENDDO
    DEALLOCATE( wr, wi )

    if (tbbox) call wait_bbox_send_to_complete ( nsend )
    !
    CALL mp_sum( becp%k( :, ibnd ), intra_bgrp_comm )
    CALL stop_clock( 'calbec_rs' )
    !
    RETURN

  END SUBROUTINE calbec_rs_k
    !
    !--------------------------------------------------------------------------
    SUBROUTINE s_psir_gamma( ibnd, last )
    !--------------------------------------------------------------------------
    !! This routine applies the \(S\) matrix to wfc ibnd (and wfc ibnd+1 if
    !! \(\le\text{last}\) ) stored in real space in psic, and puts the results
    !! again in psic for backtransforming.
    !! Requires \(\text{becp%r}\) (\(\texttt{calbecr}\) in REAL SPACE) and
    !! \(\text{betasave}\) (from \(\texttt{betapointlist}\) in \(\texttt{realus}\)).
    !
    !! WARNING: for the sake of speed, no checks performed in this subroutine.
    !
    !! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu,
    !! and S. de Gironcoli(2020).
    !
      USE kinds,                  ONLY : DP
      USE cell_base,              ONLY : omega
      USE ions_base,              ONLY : nat, nsp, ityp
      USE uspp_param,             ONLY : nh, nhm
      USE uspp,                   ONLY : qq_at, indv_ijkb0
      USE becmod,                 ONLY : bec_type, becp
      USE fft_base,               ONLY : dffts
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, last
      !
      INTEGER :: ih, nt, ia, mbia, ijkb0, box_ir
      REAL(DP) :: fac
      REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
      !
      CALL start_clock( 's_psir' )

      IF( dffts%has_task_groups ) CALL errore( 's_psir_gamma', 'task_groups not implemented', 1 )

      ALLOCATE( w1(nhm), w2(nhm) ) ; if ( ibnd+1 > last) w2 = 0.D0
      !
      fac = sqrt(omega)
      !
      DO nt = 1, nsp
         !
         DO ia = 1, nat
            !
            IF ( ityp(ia) == nt ) THEN
               !
               mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE
               !print *, "mbia=",mbia
               !
               ijkb0 = indv_ijkb0(ia)
               !
               !$omp parallel
               !$omp do
               DO ih = 1, nh(nt)
                  w1(ih) = fac * SUM(qq_at(ih,1:nh(nt),ia) * becp%r(ijkb0+1:ijkb0+nh(nt), ibnd))
                  IF ( ibnd+1 <= last ) &
                     w2(ih) = fac * SUM(qq_at(ih,1:nh(nt),ia) * becp%r(ijkb0+1:ijkb0+nh(nt), ibnd+1))
               ENDDO
               !$omp end do
               !$omp do
               DO box_ir = bbox_s(ia), bbox_e(ia)
                  box_psic( box_ir ) = SUM(betasave(box_ir,1:nh(nt))*cmplx(w1(1:nh(nt)),w2(1:nh(nt)),kind=DP))
               ENDDO
               !$omp end do
               !$omp end parallel
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      DEALLOCATE( w1, w2 )

      call add_box_to_psic ( )

      CALL stop_clock( 's_psir' )
      RETURN
      !
  END SUBROUTINE s_psir_gamma
  !
  !---------------------------------------------------------------------------
  SUBROUTINE s_psir_k( ibnd, last )
      !--------------------------------------------------------------------------
      !! Same as \(\texttt{s_psir_gamma}\) but for generalised k-point scheme i.e.:
      !! 1) Only one band is considered at a time;
      !! 2) \(\text{becp}\) is a complex entity now.
      !
      !! Derived from \(\texttt{s_psir_gamma}\) by OBM 061108.
      !! k-point phase factor fixed by SdG 030816.
      !
      USE kinds,                  ONLY : DP
      USE wvfct,                  ONLY : current_k
      USE cell_base,              ONLY : omega
      USE ions_base,              ONLY : nat, nsp, ityp
      USE uspp_param,             ONLY : nh, nhm
      USE uspp,                   ONLY : qq_at, indv_ijkb0
      USE becmod,                 ONLY : bec_type, becp
      USE fft_base,               ONLY : dffts
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ibnd, last
      !
      INTEGER :: ih, nt, ia, mbia, ijkb0, box_ir
      REAL(DP) :: fac
      COMPLEX(DP) , ALLOCATABLE :: w1(:)
      !
      REAL(DP), EXTERNAL :: ddot
      !

      CALL start_clock( 's_psir' )

      IF( dffts%has_task_groups ) CALL errore( 's_psir_k', 'task_groups not implemented', 1 )

      call set_xkphase(current_k)
      !
      fac = sqrt(omega)
      !
      ALLOCATE( w1(nhm) )
      DO nt = 1, nsp
         !
         DO ia = 1, nat
            !
            IF ( ityp(ia) == nt ) THEN
               !
               mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE
               !print *, "mbia=",mbia
               !
               ijkb0 = indv_ijkb0(ia)
               !
               !$omp parallel
               !$omp do
               DO ih = 1, nh(nt)
                  w1(ih) = fac * SUM(qq_at(ih,1:nh(nt),ia) * becp%k(ijkb0+1:ijkb0+nh(nt),ibnd))
               ENDDO
               !$omp end do
               !$omp do
               DO box_ir = bbox_s(ia), bbox_e(ia)
                  box_psic( box_ir ) = SUM(xkphase(box_ir)*betasave(box_ir,1:nh(nt))*w1(1:nh(nt)))
               ENDDO
               !$omp end do
               !$omp end parallel
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      DEALLOCATE( w1 )

      call add_box_to_psic ( )

      CALL stop_clock( 's_psir' )
      RETURN
      !
  END SUBROUTINE s_psir_k
  !
  !---------------------------------------------------------------------------
  SUBROUTINE add_vuspsir_gamma( ibnd, last )
  !--------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a vector transformed
  !! in real space contained in \(\text{psic}\).
  !! The index \(\text{ibnd}\) runs up to band \(\text{last}\).
  !! Requires the products of psi with all beta functions in array
  !! \(\text{becp%r(nkb,last)}\) (calculated by \(\texttt{calbecr}\) in
  !! real space).
  !
  !! WARNING: for the sake of speed, no checks performed in this subroutine.
  !
  !! Subroutine written by Dario Rocca, modified by O. Baris Malcioglu.
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE ions_base,              ONLY : nat, nsp, ityp
  USE uspp_param,             ONLY : nh, nhm
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq, indv_ijkb0
  USE becmod,                 ONLY : bec_type, becp
  USE fft_base,               ONLY : dffts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, last
  !
  INTEGER :: ih, nt, ia, mbia, ijkb0, box_ir
  REAL(DP) :: fac
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: w1, w2
  !
  CALL start_clock( 'add_vuspsir' )

  IF( dffts%has_task_groups ) CALL errore( 'add_vuspsir_gamma', 'task_groups not implemented', 1 )
  !
  fac = sqrt(omega)
  !
  ALLOCATE( w1(nhm), w2(nhm) ) ; IF ( ibnd+1 > last) w2 = 0.D0
  DO nt = 1, nsp
     !
     DO ia = 1, nat
        !
        IF ( ityp(ia) == nt ) THEN
           !
           mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE
           !
           ijkb0 = indv_ijkb0(ia)
           !
           !$omp parallel
           !$omp do
           DO ih = 1, nh(nt)
              w1(ih) = fac * SUM(deeq(ih,1:nh(nt),ia,current_spin) * becp%r(ijkb0+1:ijkb0+nh(nt),ibnd))
              IF ( ibnd+1 <= last ) &
                 w2(ih) = fac * SUM(deeq(ih,1:nh(nt),ia,current_spin) * becp%r(ijkb0+1:ijkb0+nh(nt),ibnd+1))
           ENDDO
           !$omp end do
           !$omp do
           DO box_ir = bbox_s(ia), bbox_e(ia)
              box_psic( box_ir ) = SUM(betasave(box_ir,1:nh(nt))*cmplx( w1(1:nh(nt)), w2(1:nh(nt)) ,kind=DP))
           ENDDO
           !$omp end do
           !$omp end parallel
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  DEALLOCATE( w1, w2 )

  call add_box_to_psic ( )

  CALL stop_clock( 'add_vuspsir' )
  RETURN
  !
  END SUBROUTINE add_vuspsir_gamma
  !
  !----------------------------------------------------------------------------
  SUBROUTINE add_vuspsir_k( ibnd, last )
  !--------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a vector transformed
  !! in real space contained in \(\text{psic}\).
  !! \(\text{ibnd}\) is an index that runs up to band \(\text{last}\).
  !! Requires the products of psi with all beta functions in array
  !! \(\text{becp(nkb,last)}\) (calculated by \(\texttt{calbecr}\) in real space).
  !
  !! WARNING: for the sake of speed, no checks performed in this subroutine.
  !
  !! Subroutine written by Stefano de Gironcoli, modified by O. Baris Malcioglu.
  !! k-point phase factor fixed by SdG 030816.
  !
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE cell_base,              ONLY : omega
  USE ions_base,              ONLY : nat, nsp, ityp
  USE uspp_param,             ONLY : nh, nhm
  USE lsda_mod,               ONLY : current_spin
  USE uspp,                   ONLY : deeq, indv_ijkb0
  USE becmod,                 ONLY : bec_type, becp
  USE fft_base,               ONLY : dffts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ibnd, last
  !
  INTEGER :: ih, nt, ia, mbia, ijkb0, box_ir
  REAL(DP) :: fac
  !
  COMPLEX(DP), ALLOCATABLE :: w1(:)
  !
  CALL start_clock( 'add_vuspsir' )

  IF( dffts%has_task_groups ) CALL errore( 'add_vuspsir_k', 'task_groups not implemented', 1 )
  !
  call set_xkphase(current_k)
  !
  fac = sqrt(omega)
  !
  ALLOCATE( w1(nhm))
  DO nt = 1, nsp
     !
     DO ia = 1, nat
        !
        IF ( ityp(ia) == nt ) THEN
           !
           mbia = maxbbox_beta(ia) ; IF ( mbia == 0 ) CYCLE

           ijkb0 = indv_ijkb0(ia)

           !$omp parallel
           !$omp do
           DO ih = 1, nh(nt)
              w1(ih) = fac * SUM( deeq(ih,1:nh(nt),ia,current_spin) * becp%k(ijkb0+1:ijkb0+nh(nt),ibnd))
           ENDDO
           !$omp end do
           !$omp do
           DO box_ir = bbox_s(ia), bbox_e(ia)
              box_psic( box_ir ) = SUM(xkphase(box_ir)*betasave(box_ir,1:nh(nt))*w1(1:nh(nt)))
           ENDDO
           !$omp end do
           !$omp end parallel
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  DEALLOCATE( w1 )

  call add_box_to_psic ( )

  CALL stop_clock( 'add_vuspsir' )
  RETURN
  !
  END SUBROUTINE add_vuspsir_k

  !--------------------------------------------------------------------------
  SUBROUTINE add_box_to_psic ( )
    !--------------------------------------------------------------------------
    !! the box-friendly array box_psic is (collected from across intra_bbox_comm
    !! and) used to update psic.
    !
    USE wavefunctions, ONLY : psic
    USE ions_base,     ONLY : nat
    !
    IMPLICIT NONE
    !
    INTEGER :: ia, box_ir, nsend

    if (tbbox) call bbox_to_box ( nsend )

    !$omp parallel
    DO ia = 1, nat
       !$omp do
       DO box_ir = box_s(ia), box_e(ia)
          psic( box_beta( box_ir ) ) = psic( box_beta( box_ir ) ) + box_psic ( box_ir )
       END DO
       !$omp end do
    END DO
    !$omp end parallel

    if (tbbox) call wait_bbox_send_to_complete ( nsend )

    RETURN

  END SUBROUTINE add_box_to_psic
  !
  !--------------------------------------------------------------------------
  SUBROUTINE invfft_orbital_gamma( orbital, ibnd, last, conserved )
    !--------------------------------------------------------------------------
    !! This driver subroutine transforms the given orbital using FFT and puts the
    !! result in \(\text{psic}\).
    !
    !! WARNING: in order to be fast, no checks on the supplied data are performed!
    !
    !! OBM 241008.
    !
    USE wavefunctions, ONLY : psic
    USE klist,         ONLY : ngk, igk_k
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts
    USE fft_interfaces,ONLY : invfft
    USE fft_helper_subroutines,   ONLY : fftx_ntgrp, tg_get_recip_inc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being transformed
    INTEGER, INTENT(in) :: last
    !! index of the last band that you want to transform (usually the
    !! total number of bands but can be different in band parallelization)
    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    !! the array of orbitals to be transformed
    LOGICAL, OPTIONAL :: conserved
    !! if this flag is true, the orbital is stored in temporary memory
    !
    ! Internal temporary variables
    INTEGER :: j, idx, ioff, right_inc, ntgrp

    !Task groups

    !The new task group version based on vloc_psi
    !print *, "->Real space"
    CALL start_clock( 'invfft_orbital' )
    !
    IF( dffts%has_task_groups ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        CALL tg_get_recip_inc( dffts, right_inc )
        ntgrp = fftx_ntgrp(dffts)
        !
        DO idx = 1, 2*ntgrp, 2

           IF( idx + ibnd - 1 < last ) THEN
              DO j = 1, ngk(1)
                 tg_psic(dffts%nl (igk_k(j,1))+ioff) =      orbital(j,idx+ibnd-1) +&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd)
                 tg_psic(dffts%nlm(igk_k(j,1))+ioff) =conjg(orbital(j,idx+ibnd-1) -&
                      (0.0d0,1.d0) * orbital(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == last ) THEN
              DO j = 1, ngk(1)
                 tg_psic(dffts%nl (igk_k(j,1))+ioff) =        orbital(j,idx+ibnd-1)
                 tg_psic(dffts%nlm(igk_k(j,1))+ioff) = conjg( orbital(j,idx+ibnd-1) )
              ENDDO
           ENDIF

           ioff = ioff + right_inc

        ENDDO
        !
        CALL invfft ('tgWave', tg_psic, dffts)
        !
        !
        IF (present(conserved)) THEN
         IF (conserved .eqv. .true.) THEN
          IF (.not. allocated(tg_psic_temp)) ALLOCATE( tg_psic_temp( dffts%nnr_tg ) )
          tg_psic_temp=tg_psic
         ENDIF
        ENDIF

    ELSE !Task groups not used
        !
        !$omp parallel default(shared) private(j)
        CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)

        IF (ibnd < last) THEN
           ! two ffts at the same time
           !$omp do
           DO j = 1, ngk(1)
              psic (dffts%nl (igk_k(j,1))) =       orbital(j, ibnd) + (0.0d0,1.d0)*orbital(j, ibnd+1)
              psic (dffts%nlm(igk_k(j,1))) = conjg(orbital(j, ibnd) - (0.0d0,1.d0)*orbital(j, ibnd+1))
           ENDDO
           !$omp end do
        ELSE
           !$omp do
           DO j = 1, ngk(1)
              psic (dffts%nl (igk_k(j,1))) =       orbital(j, ibnd)
              psic (dffts%nlm(igk_k(j,1))) = conjg(orbital(j, ibnd))
           ENDDO
           !$omp end do
        ENDIF
        !$omp end parallel
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
  !--------------------------------------------------------------------------
  SUBROUTINE fwfft_orbital_gamma( orbital, ibnd, last, conserved, add_to_orbital )
    !--------------------------------------------------------------------------
    !! This driver subroutine -back- transforms the given contribution using FFT from
    !! the already existent data in \(\text{psic}\) and return it in (or optionally
    !! add it to) orbital.
    !
    !! WARNING 1: this subroutine does not reset the orbital, use carefully!
    !! WARNING 2: in order to be fast, no checks on the supplied data are performed!
    !
    !! OBM 241008,
    !! SdG 130420.
    !
    USE wavefunctions, ONLY : psic
    USE klist,         ONLY : ngk, igk_k
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts
    USE fft_interfaces,ONLY : fwfft
    USE fft_helper_subroutines,ONLY : fftx_ntgrp, tg_get_recip_inc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being transformed
    INTEGER, INTENT(IN) :: last
    !! index of the last band that you want to transform (usually the
    !! total number of bands but can be different in band parallelization)
    COMPLEX(DP),INTENT(inout) :: orbital(:,:)
    !! the array of orbitals to be returned (or updated)
    LOGICAL, OPTIONAL :: conserved
    !! if this flag is true, the orbital is stored in temporary memory
    LOGICAL, OPTIONAL :: add_to_orbital
    !! if this flag is true, the result is added to (rather than stored into) orbital
    !
    ! Internal temporary variables
    COMPLEX(DP) :: fp, fm
    INTEGER :: j, idx, ioff, right_inc, ntgrp
    LOGICAL :: add_to_orbital_

    !Task groups
    !print *, "->fourier space"
    CALL start_clock( 'fwfft_orbital' )
    !
    add_to_orbital_=.FALSE. ; IF( present(add_to_orbital)) add_to_orbital_ = add_to_orbital
    !
    !New task_groups versions
    IF ( dffts%has_task_groups ) THEN
       !
       CALL fwfft ('tgWave', tg_psic, dffts )
       !
       ioff   = 0
       CALL tg_get_recip_inc( dffts, right_inc )
       ntgrp = fftx_ntgrp(dffts)
       !
       DO idx = 1, 2*ntgrp, 2
          !
          IF ( idx + ibnd - 1 < last ) THEN
             DO j = 1, ngk(1)
                fp= ( tg_psic( dffts%nl(igk_k(j,1)) + ioff ) +  &
                      tg_psic( dffts%nlm(igk_k(j,1)) + ioff ) ) * 0.5d0
                fm= ( tg_psic( dffts%nl(igk_k(j,1)) + ioff ) -  &
                      tg_psic( dffts%nlm(igk_k(j,1)) + ioff ) ) * 0.5d0
                IF ( add_to_orbital_ ) THEN
                   orbital (j, ibnd+idx-1) = orbital (j, ibnd+idx-1) + cmplx( dble(fp), aimag(fm),kind=DP)
                   orbital (j, ibnd+idx  ) = orbital (j, ibnd+idx  ) + cmplx(aimag(fp),- dble(fm),kind=DP)
                ELSE
                   orbital (j, ibnd+idx-1) = cmplx( dble(fp), aimag(fm),kind=DP)
                   orbital (j, ibnd+idx  ) = cmplx(aimag(fp),- dble(fm),kind=DP)
                ENDIF
             ENDDO
          ELSEIF( idx + ibnd - 1 == last ) THEN
             DO j = 1, ngk(1)
                IF ( add_to_orbital_ ) THEN
                   orbital (j, ibnd+idx-1) = orbital (j, ibnd+idx-1) + tg_psic( dffts%nl(igk_k(j,1)) + ioff )
                ELSE
                   orbital (j, ibnd+idx-1) = tg_psic( dffts%nl(igk_k(j,1)) + ioff )
                ENDIF
             ENDDO
          ENDIF
          !
          ioff = ioff + right_inc
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

          IF( add_to_orbital_ ) THEN
             !$omp parallel do
             DO j = 1, ngk(1)
                fp = (psic (dffts%nl(igk_k(j,1))) + psic (dffts%nlm(igk_k(j,1))))*0.5d0
                fm = (psic (dffts%nl(igk_k(j,1))) - psic (dffts%nlm(igk_k(j,1))))*0.5d0
                orbital( j, ibnd)   = orbital( j, ibnd)   + cmplx( dble(fp), aimag(fm),kind=DP)
                orbital( j, ibnd+1) = orbital( j, ibnd+1) + cmplx(aimag(fp),- dble(fm),kind=DP)
             ENDDO
             !$omp end parallel do
          ELSE
             !$omp parallel do
             DO j = 1, ngk(1)
                fp = (psic (dffts%nl(igk_k(j,1))) + psic (dffts%nlm(igk_k(j,1))))*0.5d0
                fm = (psic (dffts%nl(igk_k(j,1))) - psic (dffts%nlm(igk_k(j,1))))*0.5d0
                orbital( j, ibnd)   = cmplx( dble(fp), aimag(fm),kind=DP)
                orbital( j, ibnd+1) = cmplx(aimag(fp),- dble(fm),kind=DP)
             ENDDO
             !$omp end parallel do
          ENDIF
       ELSE
          IF( add_to_orbital_ ) THEN
             !$omp parallel do
             DO j = 1, ngk(1)
                orbital(j, ibnd)   =  orbital(j, ibnd) +  psic (dffts%nl(igk_k(j,1)))
             ENDDO
             !$omp end parallel do
          ELSE
             !$omp parallel do
             DO j = 1, ngk(1)
                orbital(j, ibnd)   =  psic (dffts%nl(igk_k(j,1)))
             ENDDO
             !$omp end parallel do
          ENDIF
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
  SUBROUTINE invfft_orbital_k( orbital, ibnd, last, ik, conserved )
    !--------------------------------------------------------------------------
    !! This subroutine transforms the given orbital using FFT and puts the result
    !! in \(\text{psic}\).  
    !
    !! WARNING: in order to be fast, no checks on the supplied data are performed!
    !
    !! OBM 110908
    !
    USE kinds,                    ONLY : DP
    USE wavefunctions,            ONLY : psic
    USE klist,                    ONLY : ngk, igk_k
    USE wvfct,                    ONLY : current_k
    USE fft_base,                 ONLY : dffts
    USE fft_interfaces,           ONLY : invfft
    USE fft_helper_subroutines,   ONLY : fftx_ntgrp, tg_get_recip_inc

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being transformed
    INTEGER, INTENT(in) :: last
    !! index of the last band that you want to transform (usually the total
    !! number of bands but can be different in band parallelization)
    COMPLEX(DP),INTENT(in) :: orbital(:,:)
    !! the array of orbitals to be transformed
    INTEGER, OPTIONAL :: ik
    !! index of the desired k-point
    LOGICAL, OPTIONAL :: conserved
    !! if this flag is true, the orbital is stored in temporary memory
    !
    ! Internal variables
    INTEGER :: ioff, idx, ik_ , right_inc, ntgrp, ig

    CALL start_clock( 'invfft_orbital' )

    ! current_k  variable  must contain the index of the desired kpoint
    ik_ = current_k ; if (present(ik)) ik_ = ik

    IF( dffts%has_task_groups ) THEN
       !
       tg_psic = ( 0.D0, 0.D0 )
       ioff   = 0
       CALL tg_get_recip_inc( dffts, right_inc )
       ntgrp = fftx_ntgrp(dffts)
       !
       DO idx = 1, ntgrp
          !
          IF( idx + ibnd - 1 <= last ) THEN
             !DO j = 1, size(orbital,1)
             tg_psic( dffts%nl( igk_k(:, ik_) ) + ioff ) = orbital(:,idx+ibnd-1)
             !END DO
          ENDIF

          ioff = ioff + right_inc

       ENDDO
       !
       CALL invfft ('tgWave', tg_psic, dffts)
       IF (present(conserved)) THEN
          IF (conserved .eqv. .true.) THEN
             IF (.not. allocated(tg_psic_temp)) &
                  &ALLOCATE( tg_psic_temp( dffts%nnr_tg ) )
             tg_psic_temp=tg_psic
          ENDIF
       ENDIF
       !
    ELSE  !non task_groups version
       !
       !$omp parallel default(shared) private(ig)
       CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
       !$omp do
       do ig = 1, ngk(ik_)
          psic(dffts%nl(igk_k(ig, ik_))) = orbital(ig,ibnd)
       end do
       !$omp end do
       !$omp end parallel
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
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fwfft_orbital_k( orbital, ibnd, last, ik, conserved, add_to_orbital )
    !-------------------------------------------------------------------------
    !! This driver subroutine -back- transforms the given contribution using FFT from
    !! the already existent data in \(\text{psic}\) and return it in (or optionally
    !! add it to) orbital.
    !
    !! WARNING 1: this subroutine does not reset the orbital, use carefully!
    !! WARNING 2: in order to be fast, no checks on the supplied data are performed!
    !
    !! OBM 241008,
    !! SdG 130420.
    !
    USE wavefunctions,            ONLY : psic
    USE klist,                    ONLY : ngk, igk_k
    USE wvfct,                    ONLY : current_k
    USE kinds,                    ONLY : DP
    USE fft_base,                 ONLY : dffts
    USE fft_interfaces,           ONLY : fwfft
    USE fft_helper_subroutines,   ONLY : fftx_ntgrp, tg_get_recip_inc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being transformed
    INTEGER, INTENT(in) :: last
    !! index of the last band that you want to transform (usually the
    !! total number of bands but can be different in band parallelization)
    COMPLEX(DP),INTENT(inout) :: orbital(:,:)
    !! the array of orbitals to be returned (or updated)
    INTEGER, OPTIONAL :: ik
    !! the index of the desired kpoint
    LOGICAL, OPTIONAL :: conserved
    !! if this flag is true, the orbital is stored in temporary memory
    LOGICAL, OPTIONAL :: add_to_orbital
    !! if this flag is true, the result is added to (rather than stored into) orbital
    !
    ! Internal variables
    INTEGER :: ioff, idx, ik_ , right_inc, ntgrp, ig
    LOGICAL :: add_to_orbital_
    !
    CALL start_clock( 'fwfft_orbital' )
    !
    add_to_orbital_=.FALSE. ; IF( present(add_to_orbital)) add_to_orbital_ = add_to_orbital
    !
    ! current_k  variable  must contain the index of the desired kpoint
    ik_ = current_k ; if (present(ik)) ik_ = ik

    IF( dffts%has_task_groups ) THEN
       !
       CALL fwfft ('tgWave', tg_psic, dffts)
       !
       ioff   = 0
       CALL tg_get_recip_inc( dffts, right_inc )
       ntgrp = fftx_ntgrp(dffts)
       !
       DO idx = 1, ntgrp
          !
          IF( idx + ibnd - 1 <= last ) THEN
             IF( add_to_orbital_ ) THEN
                orbital (:, ibnd+idx-1) = orbital (:, ibnd+idx-1) + tg_psic( dffts%nl(igk_k(:,ik_)) + ioff )
             ELSE
                orbital (:, ibnd+idx-1) = tg_psic( dffts%nl(igk_k(:,ik_)) + ioff )
             END IF
          ENDIF
          !
          ioff = ioff + right_inc
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
       IF( add_to_orbital_ ) THEN
          !$omp parallel do default(shared) private(ig)
          do ig=1,ngk(ik_)
             orbital(ig,ibnd) = orbital(ig,ibnd) + psic(dffts%nl(igk_k(ig,ik_)))
          end do
          !$omp end parallel do
       ELSE
          !$omp parallel do default(shared) private(ig)
          do ig=1,ngk(ik_)
             orbital(ig,ibnd) = psic(dffts%nl(igk_k(ig,ik_)))
          end do
          !$omp end parallel do
       END IF
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
  SUBROUTINE v_loc_psir( ibnd, last )
    !--------------------------------------------------------------------------
    !! Basically the same thing as \(\texttt{v_loc}\) but without implicit FFT
    !! modified for real space implementation.
    !
    !! OBM 241008
    !
    USE wavefunctions, ONLY : psic
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts
    USE scf,           ONLY : vrs
    USE lsda_mod,      ONLY : current_spin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being operated on
    INTEGER, INTENT(in) :: last
    !! index of the last band that you want to operate on
    !
    ! Internal temporary variables
    INTEGER :: j
    !Task groups
    REAL(DP), ALLOCATABLE :: tg_v(:)
    !
    CALL start_clock( 'v_loc_psir' )

    IF( dffts%has_task_groups ) THEN
        IF (ibnd == 1 ) THEN
          CALL tg_gather( dffts, vrs(:,current_spin), tg_v )
          !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        ENDIF
        !
        !$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
           tg_psic (j) = tg_psic (j) + tg_psic_temp (j) * tg_v(j)
        ENDDO
        !$omp end parallel do
        !
        DEALLOCATE( tg_v )
     ELSE
        !   product with the potential v on the smooth grid
        !
        !$omp parallel do
        DO j = 1, dffts%nnr
           psic (j) = psic (j) + psic_temp (j) * vrs(j,current_spin)
        ENDDO
        !$omp end parallel do
     ENDIF
  CALL stop_clock( 'v_loc_psir' )
  END SUBROUTINE v_loc_psir
  !
  !--------------------------------------------------------------------------
  SUBROUTINE v_loc_psir_inplace( ibnd, last )
    !--------------------------------------------------------------------------
    !! The same thing as \(\texttt{v_loc_psir}\), but:
    !! - on input \(\text{psic}\) contains the wavefunction;
    !! - on output \(\text{psic}\) is overwritten to contain \(\texttt{v_loc_psir}\).
    !! Therefore must be the first term to be considered when building \(\text{hpsi}\).
    !
    !! SdG 290716
    !
    USE wavefunctions, ONLY : psic
    USE kinds,         ONLY : DP
    USE fft_base,      ONLY : dffts
    USE scf,           ONLY : vrs
    USE lsda_mod,      ONLY : current_spin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ibnd
    !! index of the band currently being operated on
    INTEGER, INTENT(in) :: last
    !! index of the last band that you want to operate on
    !
    ! ... local variables
    !
    INTEGER :: j  !Task groups
    REAL(DP), ALLOCATABLE :: tg_v(:)
    !
    CALL start_clock( 'v_loc_psir' )

    IF( dffts%has_task_groups ) THEN
        IF (ibnd == 1 ) THEN
          CALL tg_gather( dffts, vrs(:,current_spin), tg_v )
          !if ibnd==1 this is a new calculation, and tg_v should be distributed.
        ENDIF
        !
        !$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
           tg_psic (j) = tg_v(j) * tg_psic(j)
        ENDDO
        !$omp end parallel do
        !
        DEALLOCATE( tg_v )
    ELSE
       !   product with the potential v on the smooth grid
       !
       !$omp parallel do
       DO j = 1, dffts%nnr
          psic (j) = vrs(j,current_spin) * psic(j)
       ENDDO
       !$omp end parallel do
    ENDIF
  CALL stop_clock( 'v_loc_psir' )
  END SUBROUTINE v_loc_psir_inplace
    !--------------------------------------------------------------------------
  !
END MODULE realus
