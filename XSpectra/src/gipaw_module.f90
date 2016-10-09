!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE gipaw_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for GIPAW calculations
  !
  USE kinds, ONLY : DP
  USE parameters, ONLY : npk, ntypx, lmaxx

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER:: natx=20000 ! max number of atoms
  INTEGER, PARAMETER:: nbrx=14   ! max number of beta functions
  ! alpha
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp

  ! speed of light in atomic units: c = 1/alpha
  !REAL(DP), PARAMETER :: c = 137.03599911d0

  ! avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142e23_dp

  ! number of occupied bands at each k-point
  INTEGER :: nbnd_occ(npk)

  ! alpha shift of the projector on the valence wfcs
  REAL(DP) :: alpha_pv

  ! eigenvalues and eigenfunctions at k+q
  COMPLEX(DP), ALLOCATABLE :: evq(:,:)

  ! induced current (bare term) and induced magnetic field
  REAL(DP), ALLOCATABLE :: j_bare(:,:,:,:), b_ind_r(:,:,:)

  ! induced magnetic field in reciprocal space
  COMPLEX(DP), ALLOCATABLE :: b_ind(:,:,:)

  ! convergence threshold for diagonalizationa and greenfunction
  REAL(DP) :: conv_threshold

  ! q for the perturbation (in bohrradius^{-1})
  REAL(DP) :: q_gipaw

  ! q for the EFG
  REAL(DP) :: q_efg ( ntypx )

  ! verbosity
  INTEGER :: iverbosity

  ! diagonalization method
  INTEGER :: isolve

  ! job: nmr, g_tensor, efg, hyperfine
  CHARACTER(80) :: job

  ! format for a rank-2 tensor
  CHARACTER(*), PARAMETER :: tens_fmt = '(3(5X,3(F14.4,2X)/))'

  ! for plotting the induced current and induced field
  CHARACTER(80) :: filcurr, filfield

  ! macroscopic shape for the NMR
  LOGICAL :: use_nmr_macroscopic_shape
  REAL(DP) :: nmr_macroscopic_shape ( 3, 3 )

  ! parametres for hyper-fine interaction
  CHARACTER ( LEN = 10 ) :: hfi_input_unit
  CHARACTER ( LEN = 10 ) :: hfi_output_unit
  INTEGER :: hfi_isotope(natx)
  REAL ( dp ) :: hfi_nuclear_g_factor(natx)
  LOGICAL :: radial_integral_splines
  LOGICAL :: hfi_via_reconstruction_only
  INTEGER :: hfi_extrapolation_npoints

  !<apsi>
  CHARACTER(256) :: file_reconstruction ( ntypx )
  REAL(dp) :: rc(ntypx,0:lmaxx)
  COMPLEX(dp), ALLOCATABLE :: paw_becp2 ( :, : )
  REAL(dp), ALLOCATABLE, DIMENSION ( :, : ) :: lx, ly, lz
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_rmc(:,:,:)
  !<apsi>

  ! contribution to NMR chemical shift due to core contribution
  REAL ( dp ) :: nmr_shift_core(ntypx)

CONTAINS

  !-----------------------------------------------------------------------
  ! Read in the gipaw input file.
  ! Format: &inputgipaw
  !              prefix = '...'     prefix of SCF calculation
  !              tmp_dir = '...'    scratch directory
  !              job = 'nmr' or 'g_tensor'
  !              conv_threshold = 1e-14
  !              q_gipaw = 0.01
  !              filcurr = '...'
  !              filfield = '...'
  !              iverbosity = 0
  !         /
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_readin()
    USE io_files,      ONLY : nd_nmbr, prefix, tmp_dir
    USE io_global,     ONLY : ionode
    USE us,            ONLY : spline_ps
    IMPLICIT NONE
    INTEGER :: ios
    NAMELIST /inputgipaw/ job, prefix, tmp_dir, conv_threshold, &
                         q_gipaw, iverbosity, filcurr, filfield, &
                         file_reconstruction, use_nmr_macroscopic_shape, &
                         nmr_macroscopic_shape, spline_ps, isolve, &
                         q_efg, hfi_output_unit, hfi_isotope, &
                         hfi_nuclear_g_factor, radial_integral_splines, &
                         hfi_via_reconstruction_only, &
                         hfi_extrapolation_npoints

    if ( .not. ionode ) goto 400

    CALL input_from_file()

    job = ' '
    prefix = 'pwscf'
    CALL get_environment_variable( 'ESPRESSO_TMPDIR', tmp_dir )
    IF ( TRIM( tmp_dir ) == ' ' ) tmp_dir = './scratch/'
    conv_threshold = 1e-14_dp
    q_gipaw = 0.01_dp
    iverbosity = 0
    filcurr = ' '
    filfield = ' '
    file_reconstruction ( : ) = " "
    use_nmr_macroscopic_shape = .TRUE.
    nmr_macroscopic_shape = 2.0_dp / 3.0_dp
    spline_ps = .true.    ! TRUE in this case!!!!!
    isolve = 0

    hfi_output_unit = 'MHz'
    hfi_isotope ( : ) = 0
    hfi_nuclear_g_factor ( : ) = 1.0
    radial_integral_splines = .TRUE.
    hfi_via_reconstruction_only = .TRUE.
    hfi_extrapolation_npoints = 10000
    q_efg = 1.0

    read( 5, inputgipaw, err = 200, iostat = ios )

200 call errore( 'gipaw_readin', 'reading inputgipaw namelist', abs( ios ) )

400 continue

    call gipaw_bcast_input

  END SUBROUTINE gipaw_readin

  !-----------------------------------------------------------------------
  ! Broadcast input data to all processors
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_bcast_input
#if defined(__MPI)
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : world_comm
    USE io_files,      ONLY : prefix, tmp_dir
    USE us,            ONLY : spline_ps
    implicit none
    integer :: root = 0
    call mp_bcast(job, root, world_comm)
    call mp_bcast(prefix, root, world_comm)
    call mp_bcast(tmp_dir, root, world_comm)
    call mp_bcast(conv_threshold, root, world_comm)
    call mp_bcast(q_gipaw, root, world_comm)
    call mp_bcast(iverbosity, root, world_comm)
    call mp_bcast(filcurr, root, world_comm)
    call mp_bcast(filfield, root, world_comm)
    call mp_bcast(file_reconstruction, root, world_comm)
    call mp_bcast(use_nmr_macroscopic_shape, root, world_comm)
    call mp_bcast(nmr_macroscopic_shape, root, world_comm)
    call mp_bcast(spline_ps, root, world_comm)
    call mp_bcast(isolve, root, world_comm)
    call mp_bcast(hfi_output_unit, root, world_comm)
    call mp_bcast(hfi_isotope, root, world_comm)
    call mp_bcast(hfi_nuclear_g_factor, root, world_comm)
    call mp_bcast(radial_integral_splines, root, world_comm)
    CALL mp_bcast(hfi_via_reconstruction_only, root, world_comm)
    CALL mp_bcast(hfi_extrapolation_npoints, root, world_comm)
#endif
  END SUBROUTINE gipaw_bcast_input


  !-----------------------------------------------------------------------
  ! Allocate memory for GIPAW
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_allocate
    USE lsda_mod,      ONLY : nspin, lsda
    USE gvect,         ONLY : ngm
    USE wvfct,         ONLY : nbnd, npwx
    USE ions_base,     ONLY : ntyp => nsp
    USE paw_gipaw,     ONLY : paw_recon
    USE fft_base,      ONLY : dffts

    IMPLICIT NONE

    allocate(evq(npwx,nbnd))
    allocate(j_bare(dffts%nnr,3,3,nspin), b_ind_r(dffts%nnr,3,3), b_ind(ngm,3,3))
    if (.not. allocated(paw_recon) ) allocate ( paw_recon(ntyp) )

  END SUBROUTINE gipaw_allocate


  !-----------------------------------------------------------------------
  ! Print a short summary of the calculation
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_summary
    USE io_global,     ONLY : stdout
    IMPLICIT NONE

    FLUSH( stdout )
  END SUBROUTINE gipaw_summary


  !-----------------------------------------------------------------------
  ! Open files needed for GIPAW
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_openfil
    USE io_global,        ONLY : stdout
    USE wvfct,            ONLY : nbnd, npwx
    USE ldaU,             ONLY : lda_plus_U, nwfcU
    USE io_files,         ONLY : prefix, iunhub, iunwfc, &
                                 nwordwfcU, nwordwfc, seqopn
    USE noncollin_module, ONLY : npol
    USE buffers,          ONLY : open_buffer
    USE control_flags,    ONLY : io_level
    IMPLICIT NONE
    LOGICAL            :: exst
    !
    ! ... nwordwfc is the record length (IN REAL WORDS)
    ! ... for the direct-access file containing wavefunctions
    ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
    !
    nwordwfc = nbnd*npwx*npol
    CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )

    ! ... Needed for LDA+U
    ! ... iunhub contains the (orthogonalized) atomic wfcs * S
    !
    nwordwfcU = npwx*nwfcU*npol
    IF ( lda_plus_u ) &
       CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )
    !
    RETURN
  END SUBROUTINE gipaw_openfil
  !-----------------------------------------------------------------------
  ! Print timings
  !-----------------------------------------------------------------------
  SUBROUTINE print_clock_gipaw
    USE io_global,  ONLY : stdout
    IMPLICIT NONE

    WRITE( stdout, * )
    call print_clock ('GIPAW')
    WRITE( stdout,  * ) '    INITIALIZATION: '
    call print_clock ('gipaw_setup')
    WRITE( stdout, * )
    call print_clock ('greenf')
    call print_clock ('cgsolve')
    call print_clock ('ch_psi')
    call print_clock ('h_psi')
    WRITE( stdout, * )
    call print_clock ('u_kq')
    WRITE( stdout, * )
    call print_clock ('apply_p')
    call print_clock ('apply_vel')
    WRITE( stdout, * )
    call print_clock ('j_para')
    call print_clock ('biot_savart')
    call print_clock ('c_sigma')
    WRITE( stdout, * )
    WRITE( stdout, * ) '     General routines'
    call print_clock ('calbec')
    call print_clock ('fft')
    call print_clock ('ffts')
    call print_clock ('fftw')
    call print_clock ('cinterpolate')
    call print_clock ('davcio')
    call print_clock ('write_rec')
    WRITE( stdout, * )
#if defined(__MPI)
    WRITE( stdout,  * ) '     Parallel routines'
    call print_clock ('reduce')
#endif
  END SUBROUTINE print_clock_gipaw



  !-----------------------------------------------------------------------
  ! GIPAW setup
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_setup
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode
    USE ions_base,     ONLY : tau, nat, ntyp => nsp, atm
    USE atom,          ONLY : rgrid
    USE wvfct,         ONLY : nbnd, et, wg, npwx
    USE lsda_mod,      ONLY : nspin, lsda
    USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
    USE gvect,         ONLY : ngm
    USE fft_base,      ONLY : dfftp
    USE gvecs,         ONLY : doublegrid
    USE klist,         ONLY : xk, degauss, ngauss, nks, nelec
    USE constants,     ONLY : degspin, pi
    USE paw_gipaw,     ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp, &
                              read_recon, set_paw_upf
    USE symm_base,     ONLY : nsym, s
    USE uspp_param,    ONLY : upf
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_max, mp_min
    USE dfunct,                 only : newd

    IMPLICIT none
    integer :: ik, nt, ibnd
    ! logical :: nlcc_any
    real(dp) :: emin, emax

    !<apsi>
    integer :: il, lm, l, m, lm1, lm2, m1, m2, abs_m1, abs_m2
    integer :: sign_m1, sign_m2, il1, il2, l1, l2, j, kkpsi, nrc
    real(dp) :: alpha_lm, beta_lm
    integer, allocatable :: lm2l ( : ), lm2m ( : )
    integer :: kkpsi_max
    real(dp), allocatable :: work(:), kinetic_aephi(:), kinetic_psphi(:)
    real(dp), allocatable :: aephi_dvloc_dr(:), psphi_dvloc_dr(:)

    real(dp) :: mysum1 ( 3, lmaxx ) !TMPTMPTMP
    real(dp) :: mysum2 ( 3, 1:lmaxx ) !TMPTMPTMP
    logical :: vloc_set

    INTEGER :: core_orb
    REAL ( dp ) :: integral, occupation
    !</apsi>

    call start_clock ('gipaw_setup')

    ! Test whether the symmetry operations map the Cartesian axis to each
    ! other - if not, remove them (how to check the k point mesh then? - oops!

!*apsi*    CALL test_symmetries ( s, nsym )

    ! initialize pseudopotentials
    call init_us_1

    ! initialise data, also for the case that no GIPAW is present
    IF ( .NOT. ALLOCATED ( paw_recon ) ) ALLOCATE ( paw_recon(ntyp) )
    paw_recon(:)%gipaw_data_in_upf_file = .FALSE.
    paw_recon(:)%paw_nbeta = 0
    paw_recon(:)%paw_nh = 0
    paw_recon(:)%paw_kkbeta = 0
    paw_recon(:)%gipaw_ncore_orbital = 0
    paw_recon(:)%vloc_present = .FALSE.
    DO nt = 1, ntyp
       paw_recon(nt)%paw_lll(:) = 0
    END DO

    ! Read in qe format
    DO nt = 1, ntyp
       CALL set_paw_upf (nt, upf(nt))
       CALL read_recon ( file_reconstruction(nt), nt, paw_recon(nt) )
       !!!paw_recon(nt)%paw_nbeta = 0
    END DO

    if ( iverbosity > 20 ) then
       ! Write the wave functions and local potentials for debugging/testing
       do nt = 1, ntyp
          do il = 1, paw_recon(nt)%paw_nbeta
             do j = 1, MIN(SIZE(rgrid(nt)%r),SIZE(paw_recon(nt)%aephi(il)%psi))
                write(1000+nt,*) rgrid(nt)%r(j), paw_recon(nt)%aephi(il)%psi(j)
                write(1100+nt,*) rgrid(nt)%r(j), paw_recon(nt)%psphi(il)%psi(j)
             end do
             write(1000+nt,*) " "
             write(1100+nt,*) " "
          end do
       end do

       do nt = 1, ntyp
          do j = 1, MIN(SIZE(rgrid(nt)%r),SIZE(paw_recon(nt)%gipaw_ae_vloc))
             write(1200+nt,*) rgrid(nt)%r(j), paw_recon(nt)%gipaw_ae_vloc(j)
             write(1300+nt,*) rgrid(nt)%r(j), paw_recon(nt)%gipaw_ps_vloc(j)
          end do
          write(1200+nt,*) " "
          write(1300+nt,*) " "
       end do
    end if

    ! initialize paw
    do nt = 1, ntyp
       do il = 1, paw_recon(nt)%paw_nbeta
          IF ( paw_recon(nt)%psphi(il)%label%rc < -0.99_dp ) THEN
             rc(nt,paw_recon(nt)%psphi(il)%label%l) = 1.6_dp
             rc(nt,paw_recon(nt)%aephi(il)%label%l) = 1.6_dp
             paw_recon(nt)%psphi(il)%label%rc &
                  = rc(nt,paw_recon(nt)%psphi(il)%label%l)
             paw_recon(nt)%aephi(il)%label%rc &
                  = rc(nt,paw_recon(nt)%aephi(il)%label%l)
          ELSE
             rc(nt,paw_recon(nt)%psphi(il)%label%l) &
                  = paw_recon(nt)%psphi(il)%label%rc
             rc(nt,paw_recon(nt)%aephi(il)%label%l) &
                  = paw_recon(nt)%aephi(il)%label%rc
          END IF
       enddo
    enddo

    call init_gipaw_1()

    allocate ( paw_vkb(npwx,paw_nkb) )
    allocate ( paw_becp(paw_nkb,nbnd) )
    allocate ( paw_becp2(paw_nkb,nbnd) )

    allocate ( radial_integral_diamagnetic(paw_nkb,paw_nkb,ntypx) )
    allocate ( radial_integral_paramagnetic(paw_nkb,paw_nkb,ntypx) )
    allocate ( radial_integral_diamagnetic_so(paw_nkb,paw_nkb,ntypx) )
    allocate ( radial_integral_paramagnetic_so(paw_nkb,paw_nkb,ntypx) )
    allocate ( radial_integral_rmc(paw_nkb,paw_nkb,ntypx) )
    radial_integral_diamagnetic = 0.0_dp
    radial_integral_paramagnetic = 0.0_dp
    radial_integral_diamagnetic_so = 0.0_dp
    radial_integral_paramagnetic_so = 0.0_dp
    radial_integral_rmc = 0.0_dp

    do nt=1, ntyp

       do il1=1, paw_recon(nt)%paw_nbeta
          l1 = paw_recon(nt)%psphi(il1)%label%l
          kkpsi = paw_recon(nt)%aephi(il1)%kkpsi

          nrc = paw_recon(nt)%psphi(il1)%label%nrc

          allocate ( work(kkpsi) )

          do il2 = 1, paw_recon(nt)%paw_nbeta
             l2 = paw_recon(nt)%psphi(il2)%label%l

             IF ( l1 /= l2 ) CYCLE

             !
             ! NMR shielding, diamagnetic
             !
             do j = 1, nrc
                work(j) = (paw_recon(nt)%aephi(il1)%psi(j)*paw_recon(nt)%aephi(il2)%psi(j)-&
                     paw_recon(nt)%psphi(il1)%psi(j)*paw_recon(nt)%psphi(il2)%psi(j))/rgrid(nt)%r(j)
             enddo

             CALL simpson( nrc, work, rgrid(nt)%rab(:nrc), &
                  radial_integral_diamagnetic(il1,il2,nt) )
             if (iverbosity > 10) then
                write(stdout,*) "DIA (NMR) :", nt, l1, l2, &
                     radial_integral_diamagnetic(il1,il2,nt) &
                     * alpha ** 2 * 1e6 * 4
             end if

             !
             ! NMR shielding, paramagnetic
             !
             ! calculate radial integration on atom site
             ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
             !
             do j = 1, nrc
                work(j) = &
                     ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j) ** 3
             end do

             call simpson( nrc, work, rgrid(nt)%rab(:nrc), &
                  radial_integral_paramagnetic(il1,il2,nt) )
             if (iverbosity > 10) then
                write(stdout,*) "PARA (NMR):", nt, l1, l2, &
                     radial_integral_paramagnetic(il1,il2,nt) &
                     * alpha ** 2 * 1e6 * 4
             end if

             ! Calculate the radial integral only if the radial potential
             !    is present
             IF ( .NOT. paw_recon(nt)%vloc_present ) CYCLE

             !
             ! g tensor, relativistic mass correction
             !
             ALLOCATE ( kinetic_aephi ( kkpsi ), kinetic_psphi ( kkpsi ) )
             CALL radial_kinetic_energy ( l2, rgrid(nt)%r(:nrc), &
                  paw_recon(nt)%aephi(il2)%psi(:nrc), kinetic_aephi(:nrc) )
             CALL radial_kinetic_energy ( l2, rgrid(nt)%r(:nrc), &
                  paw_recon(nt)%psphi(il2)%psi(:nrc), kinetic_psphi(:nrc) )

             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * kinetic_aephi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * kinetic_psphi(j) )
             end do
             DEALLOCATE ( kinetic_aephi, kinetic_psphi )

             CALL simpson ( nrc, work, rgrid(nt)%rab(:nrc), &
                  radial_integral_rmc(il1,il2,nt) )
             if (iverbosity > 10) then
                write(stdout,*) "RMC (SO)  :", nt, l1, l2, &
                     radial_integral_rmc(il1,il2,nt)
             end if

             ALLOCATE ( aephi_dvloc_dr ( nrc ), psphi_dvloc_dr ( nrc ) )

             CALL radial_derivative ( rgrid(nt)%r(:nrc), &
                  paw_recon(nt)%gipaw_ae_vloc(:nrc), &
                  aephi_dvloc_dr(:nrc) )
             CALL radial_derivative ( rgrid(nt)%r(:nrc), &
                  paw_recon(nt)%gipaw_ps_vloc(:nrc), &
                  psphi_dvloc_dr ( :nrc ) )

             !
             ! g tensor, diamagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     * rgrid(nt)%r(j)
             end do

             call simpson( nrc, work, rgrid(nt)%rab(:nrc), &
                  radial_integral_diamagnetic_so(il1,il2,nt) )
             if (iverbosity > 10) then
                write(stdout,*) "DIA (SO)  :", nt, l1, l2, &
                     radial_integral_diamagnetic_so(il1,il2,nt) * alpha
             end if

             !
             ! g tensor, paramagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j)
             end do

             if ( iverbosity > 20 ) then
                if ( l1 == 0 ) then
                   do j = 1, nrc
                      write(90,*) rgrid(nt)%r(j), work(j)*rgrid(nt)%r(j)**2
                   end do
                   write(90,*) ""
                end if
             end if

             call simpson( nrc,work,rgrid(nt)%rab(:nrc), &
                  radial_integral_paramagnetic_so(il1,il2,nt) )
             if ( iverbosity > 10 ) then
                write(stdout,*) "PARA (SO) :", nt, l1, l2, &
                     radial_integral_paramagnetic_so(il1,il2,nt) * alpha
             end if

             DEALLOCATE ( aephi_dvloc_dr, psphi_dvloc_dr )

          enddo

          DEALLOCATE ( work )

       enddo
    enddo

    ! terms for paramagnetic

    allocate ( lx ( lmaxx**2, lmaxx**2 ) )
    allocate ( ly ( lmaxx**2, lmaxx**2 ) )
    allocate ( lz ( lmaxx**2, lmaxx**2 ) )

    allocate ( lm2l ( lmaxx**2 ), lm2m ( lmaxx**2 ) )

    lm = 0
    do l = 0, lmaxx - 1
       do m = 0, l
          lm = lm + 1
          lm2l ( lm ) = l
          lm2m ( lm ) = m
          if ( m /= 0 ) then
             lm = lm + 1
             lm2l ( lm ) = l
             lm2m ( lm ) = - m
          end if
       end do
    end do

    lx = 0.0_dp
    ly = 0.0_dp
    lz = 0.0_dp
    do lm2 = 1, lmaxx**2
       do lm1 = 1, lmaxx**2
          if ( lm2l ( lm1 ) /= lm2l ( lm2 ) ) cycle

          l = lm2l ( lm1 )

          m1 = lm2m ( lm1 )
          m2 = lm2m ( lm2 )

          ! L_x, L_y
          if ( m2 == 0 ) then
             if ( m1 == -1 ) then
                lx ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             else if ( m1 == +1 ) then
                ly ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             end if
          else if ( m1 == 0 ) then
             if ( m2 == -1 ) then
                lx ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             else if ( m2 == +1 ) then
                ly ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
             end if
          else
             abs_m1 = abs ( m1 )
             abs_m2 = abs ( m2 )
             sign_m1 = sign ( 1, m1 )
             sign_m2 = sign ( 1, m2 )
             alpha_lm = sqrt(real(l*(l+1)-abs_m2*(abs_m2+1),dp))
             beta_lm  = sqrt(real(l*(l+1)-abs_m2*(abs_m2-1),dp))
             if ( abs_m1 == abs_m2 + 1 ) then
                lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * alpha_lm / 4.0_dp
                ly ( lm1, lm2 ) = ( sign_m2 + sign_m1 ) * alpha_lm / 4.0_dp &
                     / sign_m2
             else if ( abs_m1 == abs_m2 - 1 ) then
                lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * beta_lm / 4.0_dp
                ly ( lm1, lm2 ) =-( sign_m2 + sign_m1 ) * beta_lm / 4.0_dp &
                     / sign_m2
             end if
          end if

          ! L_z
          if ( m1 == - m2 ) then
             lz ( lm1, lm2 ) = - m2
          end if

       end do
    end do

    if (iverbosity > 20) then
       write(stdout,'(A)') "lx:"
       write(stdout,'(9F8.5)') lx

       write(stdout,'(A)') "ly:"
       write(stdout,'(9F8.5)') ly

       write(stdout,'(A)') "lz:"
       write(stdout,'(9F8.5)') lz

       ! Checks
       mysum1 = 0
       mysum2 = 0
       do lm2 = 1, lmaxx**2
          do lm1 = 1, lmaxx**2
             if ( lm2l ( lm1 ) /= lm2l ( lm2 ) ) cycle
             l = lm2l ( lm2 )

             mysum1(1,l+1) = mysum1(1,l+1) + lx(lm1,lm2)
             mysum2(1,l+1) = mysum2(1,l+1) + lx(lm1,lm2)**2
             mysum1(2,l+1) = mysum1(2,l+1) + ly(lm1,lm2)
             mysum2(2,l+1) = mysum2(2,l+1) + ly(lm1,lm2)**2
             mysum1(3,l+1) = mysum1(3,l+1) + lz(lm1,lm2)
             mysum2(3,l+1) = mysum2(3,l+1) + lz(lm1,lm2)**2
          end do
       end do
       write(stdout,'(A,9F8.4)') "Debug, sum1: x = ", mysum1(1,:)
       write(stdout,'(A,9F8.4)') "Debug, sum1: y = ", mysum1(2,:)
       write(stdout,'(A,9F8.4)') "Debug, sum1: z = ", mysum1(3,:)
       write(stdout,'(A,9F8.4)') "Debug, sum2: x = ", mysum2(1,:)
       write(stdout,'(A,9F8.4)') "Debug, sum2: y = ", mysum2(2,:)
       write(stdout,'(A,9F8.4)') "Debug, sum2: z = ", mysum2(3,:)
    end if

    deallocate ( lm2l, lm2m )

    ! Compute the shift due to core orbitals
    DO nt = 1, ntyp

       IF ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) CYCLE

       ALLOCATE ( work(rgrid(nt)%mesh) )

       nmr_shift_core(nt) = 0.0

       DO core_orb = 1, paw_recon(nt)%gipaw_ncore_orbital

          DO j = 1, SIZE(work)
             work(j) = paw_recon(nt)%gipaw_core_orbital(j,core_orb) ** 2 &
                  / rgrid(nt)%r(j)
          END DO

          CALL simpson( SIZE(work), work, rgrid(nt)%rab(:), &
               integral )

          occupation = 2 * ( &
               2 * paw_recon(nt)%gipaw_core_orbital_l(core_orb) + 1 )
          nmr_shift_core(nt) = nmr_shift_core(nt) + occupation * integral
       END DO

       DEALLOCATE ( work )

       nmr_shift_core(nt) = nmr_shift_core(nt) * 17.75045395 * 1e-6

    END DO

    WRITE ( stdout, '()' )
    DO nt = 1, ntyp
       IF ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) THEN
          WRITE ( stdout, '( 3A, F8.3)' ) "NMR: species ", &
               TRIM ( atm(nt) ), &
               ", no information on the core"
       ELSE
          WRITE ( stdout, '( 3A, F8.3)' ) "NMR: species ", &
               TRIM ( atm(nt) ), &
               ", contribution to shift due to core = ", &
               nmr_shift_core(nt) * 1e6
       END IF
    END DO
    WRITE ( stdout, '()' )

    ! computes the total local potential (external+scf) on the smooth grid
    call setlocal
    call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)

    ! compute the D for the pseudopotentials
    call newd

    ! set non linear core correction stuff
    !! nlcc_any = ANY ( upf(1:ntyp)%nlcc )
    !!if (nlcc_any) allocate (drc( ngm, ntyp))

    !! setup all gradient correction stuff
    !!call setup_dgc

    ! computes the number of occupied bands for each k point
    nbnd_occ (:) = 0
    do ik = 1, nks
       do ibnd = 1, nbnd
          if ( wg(ibnd,ik) > 1e-6 ) then
             nbnd_occ(ik) = ibnd
          end if
       end do
    end do

    ! computes alpha_pv
    emin = et (1, 1)
    do ik = 1, nks
      do ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
      enddo
    enddo
#if defined(__MPI)
    ! find the minimum across pools
    call mp_min( emin, inter_pool_comm )
#endif
    if (degauss.ne.0.0_dp) then
      call errore('gipaw_setup', 'implemented only for insulators', -1)
    else
      emax = et (1, 1)
      do ik = 1, nks
        do ibnd = 1, nbnd
          emax = max (emax, et (ibnd, ik) )
        enddo
      enddo
#if defined(__MPI)
      ! find the maximum across pools
      call mp_max( emax, inter_pool_comm )
#endif
      alpha_pv = 2.0_dp * (emax - emin)
    endif
    ! avoid zero value for alpha_pv
    alpha_pv = max (alpha_pv, 1.0d-2)

    call stop_clock('gipaw_setup')

  CONTAINS

    SUBROUTINE radial_kinetic_energy ( l, rdata, ydata, kin_ydata )
      USE splinelib

      INTEGER, INTENT ( IN ) :: l
      REAL(dp), INTENT ( IN ) :: rdata ( : ), ydata ( : )
      REAL(dp), INTENT ( OUT ) :: kin_ydata ( : )

      REAL(dp) :: d1

      d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
      CALL spline ( rdata, ydata, 0.0_dp, d1, kin_ydata )

      kin_ydata = - kin_ydata + l*(l+1) * ydata / rdata ** 2

    END SUBROUTINE radial_kinetic_energy

    SUBROUTINE radial_derivative ( rdata, ydata, dydata_dr )
      USE splinelib

      ! ydata passed as y * r
      REAL(dp), INTENT ( IN ) :: rdata ( : ), ydata ( : )
      REAL(dp), INTENT ( OUT ) :: dydata_dr ( : )

      INTEGER :: j
      REAL(dp) :: d1, tab_d2y ( SIZE ( ydata ) )

      d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
      CALL spline ( rdata, ydata, 0.0_dp, d1, tab_d2y )

      DO j = 1, SIZE ( ydata )
         dydata_dr ( j ) = &
              ( splint_deriv ( rdata, ydata, tab_d2y, rdata ( j ) ) &
              - ydata ( j ) / rdata ( j ) ) / rdata ( j )
      END DO

    END SUBROUTINE radial_derivative

  END SUBROUTINE gipaw_setup

  !****************************************************************************

  FUNCTION spline_integration ( xdata, ydata )

    USE splinelib, ONLY : spline

    IMPLICIT NONE

    ! Return type
    REAL ( dp ) :: spline_integration

    ! Arguments
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:)

    ! Local
    INTEGER  :: n
    REAL ( dp ) :: startd, startu
    REAL ( dp ) :: d2y(SIZE(xdata))

    !--------------------------------------------------------------------------

    IF ( SIZE ( xdata ) /= SIZE ( ydata ) ) &
         CALL errore ( "spline_interpolation", &
         "sizes of arguments do not match", 1 )

    startu = ( ydata(2) - ydata(1) ) / ( xdata(2) - xdata(1) )
    startu = 0.0
    startd = 0.0

    CALL spline ( xdata, ydata, startu, startd, d2y )

    spline_integration = 0.0
    DO n = 1, SIZE ( xdata )
       spline_integration = spline_integration &
            + splint_integr ( xdata, ydata, d2y, xdata(n) )
    END DO

  END FUNCTION spline_integration

  !****************************************************************************

  FUNCTION spline_integration_mirror ( xdata, ydata )

    !
    ! Like 'spline_integration' but assumes the function to be symmetric
    !    on x axis [i.e. f(-x) = f(x) ]

    USE splinelib, ONLY : spline

    IMPLICIT NONE

    ! Return type
    REAL ( dp ) :: spline_integration_mirror

    ! Arguments
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:)

    ! Local
    INTEGER  :: n
    REAL ( dp ) :: startd, startu
    REAL ( dp ) :: xdata2(2*SIZE(xdata)), ydata2(2*SIZE(ydata))
    REAL ( dp ) :: d2y(2*SIZE(xdata))

    !--------------------------------------------------------------------------

    IF ( SIZE ( xdata ) /= SIZE ( ydata ) ) &
         CALL errore ( "spline_interpolation", &
         "sizes of arguments do not match", 1 )

    xdata2(SIZE ( xdata ):1:-1) = - xdata
    xdata2(SIZE ( xdata )+1:) = xdata

    ydata2(SIZE ( xdata ):1:-1) = ydata
    ydata2(SIZE ( xdata )+1:) = ydata

    startu = ( ydata(2) - ydata(1) ) / ( xdata(2) - xdata(1) )
    startu = 0.0
    startd = 0.0

    CALL spline ( xdata2, ydata2, startu, startd, d2y )

    spline_integration_mirror = 0.0
    DO n = 1, SIZE ( xdata2 )
       spline_integration_mirror = spline_integration_mirror &
            + splint_integr ( xdata2, ydata2, d2y, xdata2(n) )
    END DO

    spline_integration_mirror = spline_integration_mirror / 2

  END FUNCTION spline_integration_mirror

  !****************************************************************************

  FUNCTION splint_integr ( xdata, ydata, d2y, x )

    USE splinelib, ONLY : spline

    IMPLICIT NONE

    ! Return type
    REAL ( dp ) :: splint_integr

    ! Arguments
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:), d2y(:), x

    ! Local
    INTEGER  :: khi, klo, xdim
    REAL ( dp ) :: a, b, da, db, h

    !--------------------------------------------------------------------------

    xdim = SIZE( xdata )

    klo = 1
    khi = xdim

    klo = MAX( MIN( locate( xdata, x ), ( xdim - 1 ) ), 1 )

    khi = klo + 1

    h = xdata(khi) - xdata(klo)

    a = ( xdata(khi) - x ) / h
    b = ( x - xdata(klo) ) / h
    da = -1 / h
    db =  1 / h

    splint_integr = -0.5 * ydata(klo) / da + 0.5 * ydata(khi) / db &
         + (  0.25 / da * d2y(klo) &
         + ( -0.25 / db * d2y(khi) ) ) &
         * ( h**2 ) / 6.D0

  END FUNCTION splint_integr

  !****************************************************************************

         !-------------------------------------------------------------------
         FUNCTION locate( xx, x )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           REAL(DP), INTENT(IN) :: xx(:)
           REAL(DP), INTENT(IN) :: x
           !
           INTEGER :: locate
           INTEGER :: n, jl, jm, ju
           LOGICAL :: ascnd
           !
           !
           n     = SIZE( xx )
           ascnd = ( xx(n) >= xx(1) )
           jl    = 0
           ju    = n + 1
           !
           main_loop: DO
              !
              IF ( ( ju - jl ) <= 1 ) EXIT main_loop
              !
              jm = ( ju + jl ) / 2
              !
              IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
                 !
                 jl = jm
                 !
              ELSE
                 !
                 ju = jm
                 !
              END IF
              !
           END DO main_loop
           !
           IF ( x == xx(1) ) THEN
              !
              locate = 1
              !
           ELSE IF ( x == xx(n) ) THEN
              !
              locate = n - 1
              !
           ELSE
              !
              locate = jl
              !
           END IF
           !
         END FUNCTION locate

  !****************************************************************************

  FUNCTION spline_mirror_extrapolate ( xdata, ydata, x )

    USE splinelib, ONLY : spline, splint

    IMPLICIT NONE

    ! Return type
    REAL ( dp ) :: spline_mirror_extrapolate

    ! Arguments
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:), x

    ! Local
    INTEGER  :: n
    REAL ( dp ) :: startd, startu
    REAL ( dp ) :: xdata2(2*SIZE(xdata)), ydata2(2*SIZE(ydata))
    REAL ( dp ) :: d2y(2*SIZE(xdata))

    !--------------------------------------------------------------------------

    IF ( SIZE ( xdata ) /= SIZE ( ydata ) ) &
         CALL errore ( "spline_mirror_extrapolate", &
         "sizes of arguments do not match", 1 )

    n = SIZE ( xdata )

    xdata2(n:1:-1) = - xdata
    xdata2(n+1:) = xdata

    ydata2(n:1:-1) = ydata
    ydata2(n+1:) = ydata

    startu = 0.0
    startd = 0.0

    CALL spline ( xdata2, ydata2, startu, startd, d2y )

    spline_mirror_extrapolate = splint ( xdata2, ydata2, d2y, x )

  END FUNCTION spline_mirror_extrapolate

  !****************************************************************************

  SUBROUTINE radial_extrapolation ( x, y, x_extrapolate, y_extrapolate, &
       norders )

    IMPLICIT NONE

    ! Arguments
    REAL ( dp ), INTENT ( IN ) :: x(:), y(:)
    REAL ( dp ), INTENT ( IN ) :: x_extrapolate(:)
    REAL ( dp ), INTENT ( OUT ) :: y_extrapolate(:)
    INTEGER, INTENT ( IN ) :: norders

    ! Local
    INTEGER :: n, i, j

    REAL ( dp ) :: a(0:norders,0:norders), b(0:norders), c(0:norders)

    !--------------------------------------------------------------------------

    ! Dirac delta function
    do j = 0, norders
       a(0:norders,j) = x(1:norders+1) ** j
       b(0:norders) = y(1:norders+1)
    END DO

    CALL invert ( a, norders + 1 )

    c = MATMUL ( a, b )

    y_extrapolate = 0.0
    DO i = 1, SIZE ( x_extrapolate )
       DO j = 0, norders
          y_extrapolate(i) = y_extrapolate(i) + c(j) * x_extrapolate(i) ** j
       END DO
    END DO

  END SUBROUTINE radial_extrapolation

  !****************************************************************************

  SUBROUTINE invert ( a, n )

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT ( IN ) :: n
    REAL ( dp ), INTENT ( INOUT ) :: a(n,n)

    ! Local
    INTEGER :: ipvt(n)
    INTEGER :: info
    REAL ( dp ) :: cwork(n)

    !--------------------------------------------------------------------------

    CALL dgetrf ( n, n, a, n, ipvt, info )
    IF ( info /= 0 ) &
         CALL errore ( "invert_in_hypefile", "dgetrf failed", ABS ( info ) )

    CALL dgetri ( n, a, n, ipvt, cwork, n, info )
    IF ( info /= 0 ) &
         CALL errore ( "invert_in_hypefile", "dgetri failed", ABS ( info ) )

  END SUBROUTINE invert


  !-----------------------------------------------------------------------
  ! Test whether symmetry axis map into each other
  !-----------------------------------------------------------------------
  SUBROUTINE test_symmetries ( s, nsym )

    USE cell_base, ONLY : at, bg

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT ( INOUT ) :: s(3,3,48)
    INTEGER, INTENT ( INOUT ) :: nsym

    ! Local
    INTEGER :: i, isym, j
    INTEGER :: s_new(3,3,48)
    INTEGER :: isym_now
    REAL ( dp ) :: vect(3,3), vec2(3), vec3(3), vec4(3)
    REAL ( dp ) :: vect_len, vec4_len

    !--------------------------------------------------------------------------

    vect = 0.0
    vect(1,1) = 1.0
    vect(2,2) = 1.0
    vect(3,3) = 1.0

    isym_now = 0

    DO isym = 1, nsym
       DO i = 1, 3

          write(6,*) ">>>>>>>>>> ", i, isym
          write(6,'(A,3F12.6)') "AAA ", vect(:,i)

          ! Transfer into reciprocal lattice vectors
          vec2 = MATMUL ( TRANSPOSE ( bg ), vect(:,i) )
          write(6,'(A,3F12.6)') "BBB ", vec2

          ! Multiply with symmetry operation
          vec3 = MATMUL ( s(:,:,isym), vec2 )
          write(6,'(A,3F12.6)') "CCC ", vec3

          ! Transfer into Cartesian coordinates
          vec4 = MATMUL ( at, vec3 )
          write(6,'(A,3F12.6)') "DDD ", vec4

          ! Check if the length is the same...
          vect_len = SQRT ( SUM ( vect(:,i) ** 2 ) )
          vec4_len = SQRT ( SUM ( vec4(:) ** 2 ) )

          IF ( ABS ( vect_len - vec4_len ) > 1e-8 ) THEN
             write(6,*) "DEBUG: ", vect_len, vec4_len
             CALL errore ( "test_symmetries", &
                  "length of vectors cannot change", -1 )
          END IF

          vec4 = vec4 / vec4_len
          DO j = 1, 3
             vec4(j) = ABS ( vec4(j) )
          END DO
          IF ( ABS ( vec4(1) - 1 ) + ABS ( vec4(2) ) + ABS ( vec4(3) ) > 1e-8 &
               .AND. ABS ( vec4(1) ) + ABS ( vec4(2) - 1 ) + ABS ( vec4(3) ) > 1e-8 &
               .AND. ABS ( vec4(1) ) + ABS ( vec4(2) ) + ABS ( vec4(3) - 1 ) > 1e-8 ) THEN
             ! Symmetry operation does not map one axis to another - remove it
          ELSE
             isym_now = isym_now + 1
             s_new(:,:,isym_now) = s(:,:,isym)
          END IF

          write(6,*) "<<<<<<<<<<<< ", i, isym

       END DO
    END DO

!    nsym = isym_now
!    s = s_new

  END SUBROUTINE test_symmetries


!-----------------------------------------------------------------------
END MODULE gipaw_module
!-----------------------------------------------------------------------
