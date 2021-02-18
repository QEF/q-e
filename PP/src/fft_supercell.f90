!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!-----------------------------------------------------------------------
MODULE fft_supercell
  !---------------------------------------------------------------------
  !
  ! ...  This module contains all the quantities and routines 
  ! ...  related to the supercell FFT.
  !
  USE kinds,               ONLY : DP
  USE fft_types,           ONLY : fft_type_descriptor
  USE stick_base,          ONLY : sticks_map, sticks_map_deallocate
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PUBLIC
  !
  PRIVATE :: fftcp_base_info, gveccp_init, gshells_cp, compare_dfft
  !
  TYPE( fft_type_descriptor ) :: dfftcp
  TYPE( sticks_map ) :: smap_cp
  !
  REAL(DP) :: at_cp(3,3)
  REAL(DP) :: bg_cp(3,3)
  REAL(DP) :: omega_cp
  LOGICAL :: gamma_only_x     ! CHECK ALSO npwxcp WHEN USING GAMMA TRICK
  !
  INTEGER :: nat_cp
  INTEGER, ALLOCATABLE :: ityp_cp(:)
  REAL(DP), ALLOCATABLE :: tau_cp(:,:) 
  !
  ! in the following all the G-vectors related quantities
  ! connected to the dfftcp are defined
  !
  INTEGER :: npwxcp = 0                      ! eq to npwx (in wvfct)
  INTEGER :: ngmcp = 0                       ! eq to ngm (in gvect)
  INTEGER :: ngmcpx = 0                      ! eq to ngmx (in gvect)
  INTEGER :: ngmcp_g = 0                     ! eq to ngm_g (in gvect)
  INTEGER :: nglcp = 0                       ! eq to ngl (in gvect)
  INTEGER :: gstart_cp = 2                   ! eq to gstart (in gvect)
  REAL(DP), ALLOCATABLE, TARGET :: gg_cp(:)                    ! eq to gg (in gvect)
  REAL(DP), ALLOCATABLE, TARGET :: g_cp(:,:)                   ! eq to g (in gvect)
  REAL(DP), POINTER, PROTECTED :: gl_cp(:)                     ! eq to gl (in gvect)
  INTEGER, ALLOCATABLE, TARGET :: mill_cp(:,:)                 ! eq to mill (in gvect)
  INTEGER, ALLOCATABLE, TARGET :: ig_l2g_cp(:)                 ! eq to ig_l2g (in gvect)
  INTEGER, ALLOCATABLE, TARGET, PROTECTED :: igtongl_cp(:)     ! eq to igtongl (in gvect)
  !
  INTEGER :: iunwann=224           ! unit for supercell Wannier functions
  INTEGER :: nwordwann             ! record length for Wannier functions
  !
  LOGICAL :: check_fft=.false.     ! if .true. dfftcp is built with the same inputs
                                   ! of dffts. Useful for comparison and to check errors
  !
  ! ...  end of module-scope declarations
  !
CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE setup_scell_fft
    !-------------------------------------------------------------------
    !
    ! ...  Here we set up the fft descriptor (dfftcp) for the supercell
    ! ...  commensurate to the Brillouin zone sampling.
    !
    USE io_global,           ONLY : stdout, ionode
    USE fft_base,            ONLY : smap, dffts
    USE fft_types,           ONLY : fft_type_init
    USE mp_bands,            ONLY : nproc_bgrp, intra_bgrp_comm, nyfft, &
                                    ntask_groups
    USE mp_pools,            ONLY : inter_pool_comm
    USE mp,                  ONLY : mp_max
    USE gvect,               ONLY : gcutm
    USE gvecs,               ONLY : gcutms
    USE gvecw,               ONLY : gcutw, gkcut
    USE ions_base,           ONLY : nat, tau, ityp
    USE parameters,          ONLY : ntypx
    USE recvec_subs,         ONLY : ggen, ggens
    USE klist,               ONLY : nks, xk
    USE control_flags,       ONLY : gamma_only
    USE cell_base,           ONLY : at, bg, omega
    USE cellmd,              ONLY : lmovecell
    USE realus,              ONLY : real_space
    USE symm_base,           ONLY : fft_fact
    USE read_wannier,        ONLY : num_kpts, kgrid
    USE command_line_options, ONLY : nmany_
    !
    !
    IMPLICIT NONE
    !
    INTEGER, EXTERNAL :: n_plane_waves
    INTEGER :: i, j, k, ir, ik, iat
    INTEGER :: ngmcp_
    INTEGER :: nkscp
    REAL(DP) :: rvec(3)
    REAL(DP), ALLOCATABLE :: xkcp(:,:)
    LOGICAL :: lpara
    !
    !
    CALL dealloc_sc_fft
    !
    ! ...  we find the supercell lattice vectors and volume
    !
    DO i = 1, 3
      at_cp(:,i) = at(:,i) * kgrid(i)
      bg_cp(:,i) = bg(:,i) / kgrid(i)
    ENDDO
    omega_cp = omega * num_kpts
    !
    lpara =  ( nproc_bgrp > 1 )
    gkcut = gcutw
    !
    ! ...  determine atomic positions and types in the supercell
    !
    nat_cp = nat * num_kpts
    IF ( .not. ALLOCATED(tau_cp) ) ALLOCATE( tau_cp(3,nat_cp) )
    IF ( .not. ALLOCATED(ityp_cp) ) ALLOCATE( ityp_cp(nat_cp) )
    !
    ir = 0
    !
    DO i = 1, kgrid(1)
      DO j = 1, kgrid(2)
        DO k = 1, kgrid(3)
          !
          rvec(:) = (/ i-1, j-1, k-1 /)
          CALL cryst_to_cart( 1, rvec, at, 1 )
          !
          DO iat = 1, nat
            !
            ityp_cp( ir*nat + iat ) = ityp( iat )
            !
            tau_cp(:, ir*nat + iat ) = tau(:,iat) + rvec(:)
            !
          ENDDO
          !
          ir = ir + 1
          !
        ENDDO
      ENDDO
    ENDDO
    !
    !
    ! ... uncomment the following line in order to realize dfftcp in the
    ! ... same way of dffts and check possible errors or incosistencies
    !check_fft = .true.
    !
    IF ( check_fft ) THEN
      ! 
      gamma_only_x = gamma_only
      at_cp(:,:) = at(:,:)
      bg_cp(:,:) = bg(:,:)
      omega_cp = omega
      !
      nkscp = nks
      IF ( .not. ALLOCATED(xkcp) ) ALLOCATE( xkcp(3,nkscp) )
      xkcp(:,:) = xk(:,:)
      !
      nat_cp = nat
      DEALLOCATE( tau_cp, ityp_cp )
      ALLOCATE( tau_cp(3,nat_cp) )
      ALLOCATE( ityp_cp(nat_cp) )
      ityp_cp(:) = ityp(:)
      tau_cp(:,:) = tau(:,:)
      !
      ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
      !
      IF (nks == 0) THEN
        !
        ! k-point list not available:
        ! use max(bg)/2 as an estimate of the largest k-point
        !
        gkcut = 0.5d0 * max ( &
            sqrt (sum(bg_cp (1:3, 1)**2) ), &
            sqrt (sum(bg_cp (1:3, 2)**2) ), &
            sqrt (sum(bg_cp (1:3, 3)**2) ) )
      ELSE
        gkcut = 0.0d0
        DO ik = 1, nks
          gkcut = max (gkcut, sqrt ( sum(xk (1:3, ik)**2) ) )
        ENDDO
      ENDIF
      gkcut = (sqrt (gcutw) + gkcut)**2
      !
      ! ... find maximum value among all the processors
      !
      CALL mp_max (gkcut, inter_pool_comm )
      !
    ELSE
      !
      nkscp = 1
      IF ( .not. ALLOCATED(xkcp) ) ALLOCATE( xkcp(3,nkscp) )
      xkcp(:,1) = (/ 0.D0, 0.D0, 0.D0 /)
      !
      ! force the supercell FFT grid to be a multiple of the 
      ! primitive cell FFT grid
      !
      dfftcp%nr1 = dffts%nr1 * kgrid(1)
      dfftcp%nr2 = dffts%nr2 * kgrid(2)
      dfftcp%nr3 = dffts%nr3 * kgrid(3)
      !
    ENDIF
    !
    ! ... set up the supercell fft descriptor, including parallel 
    ! ... stuff: sticks, planes, etc.
    !
    ! task group are disabled if real_space calculation of calbec is used
    dfftcp%has_task_groups = (ntask_groups >1) .and. .not. real_space
    CALL fft_type_init( dfftcp, smap_cp, "wave", gamma_only_x, lpara, intra_bgrp_comm,&
         at_cp, bg_cp, gkcut, gcutms/gkcut, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_ )
    dfftcp%rho_clock_label='ffts' ; dfftcp%wave_clock_label='fftw'
    !
    ngmcp_ = dfftcp%ngl( dfftcp%mype + 1 )
    IF ( gamma_only_x ) ngmcp_ = ( ngmcp_ + 1 ) / 2
    CALL gveccp_init( ngmcp_, intra_bgrp_comm )
    !
    !
    ! Some checks (done normally in allocate_fft)
    !
    IF (dfftcp%nnr < ngmcp) THEN
      WRITE( stdout, '(/,4x," nr1cp=",i4," nr2cp= ", i4, " nr3cp=",i4, &
          &" nnrcp= ",i8," ngmcp=",i8)') dfftcp%nr1, dfftcp%nr2, dfftcp%nr3, dfftcp%nnr, ngmcp
      CALL errore( 'setup_scell_fft', 'the nrs"s are too small!', 1 )
    ENDIF 
    !
    IF (ngmcp  <= 0)      CALL errore( 'setup_scell_fft', 'wrong ngmcp' , 1 )
    IF (dfftcp%nnr <= 0) CALL errore( 'setup_scell_fft', 'wrong nnr',  1 )
    !
    ! NB: ggen normally would have dfftp and ggens dffts... here I put dfftcp in both !!!
    CALL ggen ( dfftcp, gamma_only_x, at_cp, bg_cp, gcutm, ngmcp_g, ngmcp, &
         g_cp, gg_cp, mill_cp, ig_l2g_cp, gstart_cp )
    CALL ggens( dfftcp, gamma_only_x, at_cp, g_cp, gg_cp, mill_cp, gcutms, ngmcp )
    CALL gshells_cp( lmovecell )
    CALL fftcp_base_info( ionode, stdout )
    !
    !
    ! find the number of PWs for the supercell wfc
    npwxcp = n_plane_waves( gcutw, nkscp, xkcp, g_cp, ngmcp )
    !
    IF ( check_fft ) CALL compare_dfft( dffts, dfftcp )
    !
    !
  END SUBROUTINE setup_scell_fft
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE fftcp_base_info( ionode, stdout )
    !-------------------------------------------------------------------
    !
    LOGICAL, INTENT(IN) :: ionode
    INTEGER, INTENT(IN) :: stdout
    !
    !  Display fftcp basic information
    !
    IF (ionode) THEN
       WRITE( stdout,*)
       IF ( check_fft ) THEN
         WRITE( stdout, '(5X,"Info about the pcell FFT")')
         WRITE( stdout, '(5X,"------------------------")')
       ELSE
         WRITE( stdout, '(5X,"Info about the supercell FFT")')
         WRITE( stdout, '(5X,"----------------------------")')
       ENDIF
       WRITE( stdout, '(5X,"sticks:   smooth     PW", &
                      & 5X,"G-vecs:    smooth      PW")')
       IF ( dfftcp%nproc > 1 ) THEN
          WRITE( stdout,'(5X,"Min",5X,I8,I7,13X,I9,I8)') &
             minval(dfftcp%nsp), minval(dfftcp%nsw), &
             minval(dfftcp%ngl), minval(dfftcp%nwl)
          WRITE( stdout,'(5X,"Max",5X,I8,I7,13X,I9,I8)') &
             maxval(dfftcp%nsp), maxval(dfftcp%nsw), &
             maxval(dfftcp%ngl), maxval(dfftcp%nwl)
       END IF
       WRITE( stdout,'(5X,"Sum",5X,I8,I7,13X,I9,I8)') &
          sum(dfftcp%nsp), sum(dfftcp%nsw), &
          sum(dfftcp%ngl), sum(dfftcp%nwl)
       WRITE( stdout, '(/5x,"grid: ",i10," G-vectors", 5x, &
       &               "FFT dimensions: (",i4,",",i4,",",i4,")")') &
       &         ngmcp_g, dfftcp%nr1, dfftcp%nr2, dfftcp%nr3
       IF ( .not. check_fft ) THEN
         IF ( gamma_only_x ) THEN
           WRITE( stdout, '(/5X,"Gamma-only algorithm is used &
                          --> real wavefunctions")')
         ELSE
           WRITE( stdout, '(/5X,"Gamma-only algorithm is not used &
                          --> complex wavefunctions")')
         ENDIF
       ENDIF
    ENDIF
    !
    IF(ionode) WRITE( stdout,*)
    !
    RETURN
    !
    !
  END SUBROUTINE fftcp_base_info
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE gveccp_init( ngmcp_ , comm )
    !-------------------------------------------------------------------
    !
    ! Set local and global dimensions, allocate arrays
    !
    USE mp,                  ONLY : mp_max, mp_sum
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ngmcp_
    INTEGER, INTENT(IN) :: comm  ! communicator of the group on which g-vecs are distributed
    !
    !
    ngmcp = ngmcp_
    !
    !  calculate maximum over all processors
    !
    ngmcpx = ngmcp
    CALL mp_max( ngmcpx, comm )
    !
    !  calculate sum over all processors
    !
    ngmcp_g = ngmcp
    CALL mp_sum( ngmcp_g, comm )
    !
    !  allocate arrays - only those that are always kept until the end
    !
    ALLOCATE( gg_cp(ngmcp) )
    ALLOCATE( g_cp(3, ngmcp) )
    ALLOCATE( mill_cp(3, ngmcp) )
    ALLOCATE( ig_l2g_cp(ngmcp) )
    ALLOCATE( igtongl_cp(ngmcp) )
    !
    RETURN
    !
    !
  END SUBROUTINE gveccp_init
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE gshells_cp ( vc )
    !----------------------------------------------------------------------
    !
    ! calculate number of G shells for the supercell: nglcp, and the 
    ! index ng = igtongl_cp(ig) that gives the shell index ng for 
    ! (local) G-vector of index ig
    !
    USE constants,          ONLY : eps8
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: vc
    !
    INTEGER :: ng, igl
    !
    !
    IF ( vc ) THEN
      !
      ! in case of a variable cell run each G vector has its shell
      !
      nglcp = ngmcp
      gl_cp => gg_cp
      !
      DO ng = 1, ngmcp
        igtongl_cp (ng) = ng
      ENDDO
      !
    ELSE
      !
      ! G vectors are grouped in shells with the same norm
      !
      nglcp = 1
      igtongl_cp (1) = 1
      !
      DO ng = 2, ngmcp
        IF (gg_cp (ng) > gg_cp (ng - 1) + eps8) THEN
          nglcp = nglcp + 1
        ENDIF
        igtongl_cp (ng) = nglcp
      ENDDO
      !
      ALLOCATE ( gl_cp(nglcp) )
      gl_cp (1) = gg_cp (1)
      igl = 1
      !
      DO ng = 2, ngmcp
        IF (gg_cp (ng) > gg_cp (ng - 1) + eps8) THEN
          igl = igl + 1
          gl_cp (igl) = gg_cp (ng)
        ENDIF
      ENDDO
      !
    IF (igl /= nglcp) CALL errore ('gshells_cp', 'igl <> ngl', nglcp)
    !
    ENDIF
    !
    !
  END SUBROUTINE gshells_cp 
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE dealloc_sc_fft
    !----------------------------------------------------------------------
    !
    !
    IMPLICIT NONE
    !
    !
    IF ( ALLOCATED(gg_cp) ) DEALLOCATE( gg_cp )
    IF ( ALLOCATED(g_cp) ) DEALLOCATE( g_cp )
    IF ( ALLOCATED(mill_cp) ) DEALLOCATE( mill_cp )
    IF ( ALLOCATED(ig_l2g_cp) ) DEALLOCATE( ig_l2g_cp )
    IF ( ALLOCATED(igtongl_cp) ) DEALLOCATE( igtongl_cp )
    !
    !
  END SUBROUTINE dealloc_sc_fft 
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE compare_dfft( dfft1, dfft2 )
    !-------------------------------------------------------------------
    !
    USE io_global,           ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    TYPE( fft_type_descriptor ), INTENT(IN) :: dfft1, dfft2
    !
    !
    IF ( dfft1%nr1 .ne. dfft2%nr1 ) THEN
      WRITE(stdout,*) "Mismatch in nr1", dfft1%nr1, dfft2%nr1
    ENDIF
    !
    IF ( dfft1%nr2 .ne. dfft2%nr2 ) THEN
      WRITE(stdout,*) "Mismatch in nr2", dfft1%nr2, dfft2%nr2
    ENDIF
    !
    IF ( dfft1%nr3 .ne. dfft2%nr3 ) THEN
      WRITE(stdout,*) "Mismatch in nr3", dfft1%nr3, dfft2%nr3
    ENDIF
    !
    IF ( dfft1%nr1x .ne. dfft2%nr1x ) THEN
      WRITE(stdout,*) "Mismatch in nr1x", dfft1%nr1x, dfft2%nr1x
    ENDIF
    !
    IF ( dfft1%nr2x .ne. dfft2%nr2x ) THEN
      WRITE(stdout,*) "Mismatch in nr2x", dfft1%nr2x, dfft2%nr2x
    ENDIF
    !
    IF ( dfft1%nr3x .ne. dfft2%nr3x ) THEN
      WRITE(stdout,*) "Mismatch in nr3x", dfft1%nr3x, dfft2%nr3x
    ENDIF
    !
    IF ( dfft1%lpara .ne. dfft2%lpara ) THEN
      WRITE(stdout,*) "Mismatch in lpara", dfft1%lpara, dfft2%lpara 
    ENDIF
    !
    IF ( dfft1%lgamma .ne. dfft2%lgamma ) THEN
      WRITE(stdout,*) "Mismatch in lgamma", dfft1%lgamma, dfft2%lgamma
    ENDIF
    !
    IF ( dfft1%root .ne. dfft2%root ) THEN
      WRITE(stdout,*) "Mismatch in root", dfft1%root, dfft2%root
    ENDIF
    !
    IF ( dfft1%comm .ne. dfft2%comm ) THEN
      WRITE(stdout,*) "Mismatch in comm", dfft1%comm, dfft2%comm
    ENDIF
    !
    IF ( dfft1%comm2 .ne. dfft2%comm2 ) THEN
      WRITE(stdout,*) "Mismatch in comm2", dfft1%comm2, dfft2%comm2
    ENDIF
    !
    IF ( dfft1%comm3 .ne. dfft2%comm3 ) THEN
      WRITE(stdout,*) "Mismatch in comm3", dfft1%comm3, dfft2%comm3
    ENDIF
    !
    IF ( dfft1%nproc .ne. dfft2%nproc ) THEN
      WRITE(stdout,*) "Mismatch in nproc", dfft1%nproc, dfft2%nproc
    ENDIF
    !
    IF ( dfft1%nproc2 .ne. dfft2%nproc2 ) THEN
      WRITE(stdout,*) "Mismatch in nproc2", dfft1%nproc2, dfft2%nproc2
    ENDIF
    !
    IF ( dfft1%nproc3 .ne. dfft2%nproc3 ) THEN
      WRITE(stdout,*) "Mismatch in nproc3", dfft1%nproc3, dfft2%nproc3
    ENDIF
    !
    IF ( dfft1%mype .ne. dfft2%mype ) THEN
      WRITE(stdout,*) "Mismatch in mype", dfft1%mype, dfft2%mype
    ENDIF
    !
    IF ( dfft1%mype2 .ne. dfft2%mype2 ) THEN
      WRITE(stdout,*) "Mismatch in mype2", dfft1%mype2, dfft2%mype2
    ENDIF
    !
    IF ( dfft1%mype3 .ne. dfft2%mype3 ) THEN
      WRITE(stdout,*) "Mismatch in mype3", dfft1%mype3, dfft2%mype3
    ENDIF
    !
    IF ( ANY(dfft1%iproc .ne. dfft2%iproc) ) THEN
      WRITE(stdout,*) "Mismatch in iproc", dfft1%iproc, dfft2%iproc
    ENDIF
    !
    IF ( ANY(dfft1%iproc2 .ne. dfft2%iproc2) ) THEN
      WRITE(stdout,*) "Mismatch in iproc2", dfft1%iproc2, dfft2%iproc2
    ENDIF
    !
    IF ( ANY(dfft1%iproc3 .ne. dfft2%iproc3) ) THEN
      WRITE(stdout,*) "Mismatch in iproc3", dfft1%iproc3, dfft2%iproc3
    ENDIF
    !
    IF ( dfft1%my_nr3p .ne. dfft2%my_nr3p ) THEN
      WRITE(stdout,*) "Mismatch in my_nr3p", dfft1%my_nr3p, dfft2%my_nr3p
    ENDIF
    !
    IF ( dfft1%my_nr2p .ne. dfft2%my_nr2p ) THEN
      WRITE(stdout,*) "Mismatch in my_nr2p", dfft1%my_nr2p, dfft2%my_nr2p
    ENDIF
    !
    IF ( dfft1%my_i0r3p .ne. dfft2%my_i0r3p ) THEN
      WRITE(stdout,*) "Mismatch in my_i0r3p", dfft1%my_i0r3p, dfft2%my_i0r3p
    ENDIF
    !
    IF ( dfft1%my_i0r2p .ne. dfft2%my_i0r2p ) THEN
      WRITE(stdout,*) "Mismatch in my_i0r2p", dfft1%my_i0r2p, dfft2%my_i0r2p
    ENDIF
    !
    IF ( ANY(dfft1%nr3p .ne. dfft2%nr3p) ) THEN
      WRITE(stdout,*) "Mismatch in nr3p", dfft1%nr3p, dfft2%nr3p
    ENDIF
    !
    IF ( ANY(dfft1%nr3p_offset .ne. dfft2%nr3p_offset) ) THEN
      WRITE(stdout,*) "Mismatch in nr3p_offset", dfft1%nr3p_offset, dfft2%nr3p_offset
    ENDIF
    !
    IF ( ANY(dfft1%nr2p .ne. dfft2%nr2p) ) THEN
      WRITE(stdout,*) "Mismatch in nr2p", dfft1%nr2p, dfft2%nr2p
    ENDIF
    !
    IF ( ANY(dfft1%nr2p_offset .ne. dfft2%nr2p_offset) ) THEN
      WRITE(stdout,*) "Mismatch in nr2p_offset", dfft1%nr2p_offset, dfft2%nr2p_offset
    ENDIF
    !
    IF ( ANY(dfft1%nr1p .ne. dfft2%nr1p) ) THEN
      WRITE(stdout,*) "Mismatch in nr1p", dfft1%nr1p, dfft2%nr1p
    ENDIF
    !
    IF ( ANY(dfft1%nr1w .ne. dfft2%nr1w) ) THEN
      WRITE(stdout,*) "Mismatch in nr1w", dfft1%nr1w, dfft2%nr1w
    ENDIF
    !
    IF ( dfft1%nr1w_tg .ne. dfft2%nr1w_tg ) THEN
      WRITE(stdout,*) "Mismatch in nr1w_tg", dfft1%nr1w_tg, dfft2%nr1w_tg
    ENDIF
    !
    IF ( ANY(dfft1%i0r3p .ne. dfft2%i0r3p) ) THEN
      WRITE(stdout,*) "Mismatch in i0r3p", dfft1%i0r3p, dfft2%i0r3p
    ENDIF
    !
    IF ( ANY(dfft1%i0r2p .ne. dfft2%i0r2p) ) THEN
      WRITE(stdout,*) "Mismatch in i0r2p", dfft1%i0r2p, dfft2%i0r2p
    ENDIF
    !
    IF ( ANY(dfft1%ir1p .ne. dfft2%ir1p) ) THEN
      WRITE(stdout,*) "Mismatch in ir1p", dfft1%ir1p, dfft2%ir1p
    ENDIF
    !
    IF ( ANY(dfft1%indp .ne. dfft2%indp) ) THEN
      WRITE(stdout,*) "Mismatch in indp", dfft1%indp, dfft2%indp
    ENDIF
    !
    IF ( ANY(dfft1%ir1w .ne. dfft2%ir1w) ) THEN
      WRITE(stdout,*) "Mismatch in ir1w", dfft1%ir1w, dfft2%ir1w
    ENDIF
    !
    IF ( ANY(dfft1%indw .ne. dfft2%indw) ) THEN
      WRITE(stdout,*) "Mismatch in indw", dfft1%indw, dfft2%indw
    ENDIF
    !
    IF ( ANY(dfft1%ir1w_tg .ne. dfft2%ir1w_tg) ) THEN
      WRITE(stdout,*) "Mismatch in ir1w_tg", dfft1%ir1w_tg, dfft2%ir1w_tg
    ENDIF
    !
    IF ( ANY(dfft1%indw_tg .ne. dfft2%indw_tg) ) THEN
      WRITE(stdout,*) "Mismatch in indw_tg", dfft1%indw_tg, dfft2%indw_tg
    ENDIF
    !
    IF ( dfft1%nst .ne. dfft2%nst ) THEN
      WRITE(stdout,*) "Mismatch in nst", dfft1%nst, dfft2%nst
    ENDIF
    !
    IF ( ANY(dfft1%nsp .ne. dfft2%nsp) ) THEN
      WRITE(stdout,*) "Mismatch in nsp", dfft1%nsp, dfft2%nsp
    ENDIF
    !
    IF ( ANY(dfft1%nsp_offset .ne. dfft2%nsp_offset) ) THEN
      WRITE(stdout,*) "Mismatch in nsp_offset", dfft1%nsp_offset, dfft2%nsp_offset
    ENDIF
    !
    IF ( ANY(dfft1%nsw .ne. dfft2%nsw) ) THEN
      WRITE(stdout,*) "Mismatch in nsw", dfft1%nsw, dfft2%nsw
    ENDIF
    !
    IF ( ANY(dfft1%nsw_offset .ne. dfft2%nsw_offset) ) THEN
      WRITE(stdout,*) "Mismatch in nsw_offset", dfft1%nsw_offset, dfft2%nsw_offset
    ENDIF
    !
    IF ( ANY(dfft1%nsw_tg .ne. dfft2%nsw_tg) ) THEN
      WRITE(stdout,*) "Mismatch in nsw_tg", dfft1%nsw_tg, dfft2%nsw_tg
    ENDIF
    !
    IF ( ANY(dfft1%ngl .ne. dfft2%ngl) ) THEN
      WRITE(stdout,*) "Mismatch in ngl", dfft1%ngl, dfft2%ngl
    ENDIF
    !
    IF ( ANY(dfft1%nwl .ne. dfft2%nwl) ) THEN
      WRITE(stdout,*) "Mismatch in nwl", dfft1%nwl, dfft2%nwl
    ENDIF
    !
    IF ( dfft1%ngm .ne. dfft2%ngm ) THEN
      WRITE(stdout,*) "Mismatch in ngm", dfft1%ngm, dfft2%ngm
    ENDIF
    !
    IF ( dfft1%ngw .ne. dfft2%ngw ) THEN
      WRITE(stdout,*) "Mismatch in ngw", dfft1%ngw, dfft2%ngw
    ENDIF
    !
    IF ( ANY(dfft1%iplp .ne. dfft2%iplp) ) THEN
      WRITE(stdout,*) "Mismatch in iplp", dfft1%iplp, dfft2%iplp
    ENDIF
    !
    IF ( ANY(dfft1%iplw .ne. dfft2%iplw) ) THEN
      WRITE(stdout,*) "Mismatch in iplw", dfft1%iplw, dfft2%iplw
    ENDIF
    !
    IF ( dfft1%nnp .ne. dfft2%nnp ) THEN
      WRITE(stdout,*) "Mismatch in nnp", dfft1%nnp, dfft2%nnp
    ENDIF
    !
    IF ( dfft1%nnr .ne. dfft2%nnr ) THEN
      WRITE(stdout,*) "Mismatch in nnr", dfft1%nnr, dfft2%nnr
    ENDIF
    !
    IF ( dfft1%nnr_tg .ne. dfft2%nnr_tg ) THEN
      WRITE(stdout,*) "Mismatch in nnr_tg", dfft1%nnr_tg, dfft2%nnr_tg
    ENDIF
    !
    IF ( ANY(dfft1%iss .ne. dfft2%iss) ) THEN
      WRITE(stdout,*) "Mismatch in iss", dfft1%iss, dfft2%iss
    ENDIF
    !
    IF ( ANY(dfft1%isind .ne. dfft2%isind) ) THEN
      WRITE(stdout,*) "Mismatch in isind", dfft1%isind, dfft2%isind
    ENDIF
    !
    IF ( ANY(dfft1%ismap .ne. dfft2%ismap) ) THEN
      WRITE(stdout,*) "Mismatch in ismap", dfft1%ismap, dfft2%ismap
    ENDIF
    !
    IF ( ANY(dfft1%nl .ne. dfft2%nl) ) THEN
      WRITE(stdout,*) "Mismatch in nl", dfft1%nl, dfft2%nl
    ENDIF
    !
    IF ( ANY(dfft1%nlm .ne. dfft2%nlm) ) THEN
      WRITE(stdout,*) "Mismatch in nlm", dfft1%nlm, dfft2%nlm
    ENDIF
    !
    IF ( ANY(dfft1%tg_snd .ne. dfft2%tg_snd) ) THEN
      WRITE(stdout,*) "Mismatch in tg_snd", dfft1%tg_snd, dfft2%tg_snd
    ENDIF
    !
    IF ( ANY(dfft1%tg_rcv .ne. dfft2%tg_rcv) ) THEN
      WRITE(stdout,*) "Mismatch in tg_rcv", dfft1%tg_rcv, dfft2%tg_rcv
    ENDIF
    !
    IF ( ANY(dfft1%tg_sdsp .ne. dfft2%tg_sdsp) ) THEN
      WRITE(stdout,*) "Mismatch in tg_sdsp", dfft1%tg_sdsp, dfft2%tg_sdsp
    ENDIF
    !
    IF ( ANY(dfft1%tg_rdsp .ne. dfft2%tg_rdsp) ) THEN
      WRITE(stdout,*) "Mismatch in tg_rdsp", dfft1%tg_rdsp, dfft2%tg_rdsp
    ENDIF
    !
    IF ( dfft1%has_task_groups .ne. dfft2%has_task_groups ) THEN
      WRITE(stdout,*) "Mismatch in has_task_groups", dfft1%has_task_groups, dfft2%has_task_groups
    ENDIF
    !
    IF ( dfft1%rho_clock_label .ne. dfft2%rho_clock_label ) THEN
      WRITE(stdout,*) "Mismatch in rho_clock_label", dfft1%rho_clock_label, dfft2%rho_clock_label
    ENDIF
    !
    IF ( dfft1%wave_clock_label .ne. dfft2%wave_clock_label ) THEN
      WRITE(stdout,*) "Mismatch in wave_clock_label", dfft1%wave_clock_label, dfft2%wave_clock_label
    ENDIF
    !
    IF ( dfft1%grid_id .ne. dfft2%grid_id ) THEN
      WRITE(stdout,*) "Mismatch in grid_id", dfft1%grid_id, dfft2%grid_id
    ENDIF
    !
    !
  END SUBROUTINE compare_dfft
  !
  !
END MODULE fft_supercell
