!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
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
  USE constants, ONLY : a0_to_cm => bohr_radius_cm
  USE parameters, ONLY : npk, ntypx, lmaxx
  
  IMPLICIT NONE
  SAVE
  
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
  
  ! job: nmr, g_tensor, efg, hyperfine
  CHARACTER(80) :: job
  
  ! format for a rank-2 tensor
  CHARACTER(*), PARAMETER :: tens_fmt = '(3(5X,3(F14.4,2X)/))'
  
  ! for plotting the induced current and induced field
  CHARACTER(80) :: filcurr, filfield
  
  ! macroscopic shape for the NMR
  LOGICAL :: use_nmr_macroscopic_shape
  REAL(DP) :: nmr_macroscopic_shape ( 3, 3 )
  
  !<apsi>
  CHARACTER(256) :: file_reconstruction ( ntypx )
  LOGICAL :: read_recon_in_paratec_fmt
  REAL(dp) :: rc(ntypx,0:lmaxx)
  COMPLEX(dp), ALLOCATABLE :: paw_becp2 ( :, : )
  REAL(dp), ALLOCATABLE, DIMENSION ( :, : ) :: lx, ly, lz
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_rmc(:,:,:)
  !<apsi>
  
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
                         read_recon_in_paratec_fmt, &
                         file_reconstruction, use_nmr_macroscopic_shape, &
                         nmr_macroscopic_shape, spline_ps
    
    if ( .not. ionode ) goto 400
    
    job = ''
    prefix = 'pwscf'
    tmp_dir = './scratch/'
    conv_threshold = 1e-14_dp
    q_gipaw = 0.01_dp
    iverbosity = 0
    filcurr = ''
    filfield = ''
    read_recon_in_paratec_fmt = .FALSE.
    file_reconstruction ( : ) = " "
    nmr_macroscopic_shape = 2.0_dp / 3.0_dp
    spline_ps = .true.    ! TRUE in this case!!!!!
    
    read( 5, inputgipaw, err = 200, iostat = ios )
    
200 call errore( 'gipaw_readin', 'reading inputgipaw namelist', abs( ios ) )
    
400 continue
    
    call gipaw_bcast_input
    
  END SUBROUTINE gipaw_readin

  
  !-----------------------------------------------------------------------
  ! Broadcast input data to all processors 
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_bcast_input
#ifdef __PARA
#include "f_defs.h"
    USE mp,            ONLY : mp_bcast
    USE io_files,      ONLY : prefix, tmp_dir
    USE us,            ONLY : spline_ps
    implicit none
    integer :: root = 0
    call mp_bcast(job, root)
    call mp_bcast(prefix, root)
    call mp_bcast(tmp_dir, root)
    call mp_bcast(conv_threshold, root)
    call mp_bcast(q_gipaw, root)
    call mp_bcast(iverbosity, root)
    call mp_bcast(filcurr, root)
    call mp_bcast(filfield, root)
    call mp_bcast(read_recon_in_paratec_fmt, root)
    call mp_bcast(file_reconstruction, root)
    call mp_bcast(use_nmr_macroscopic_shape, root)
    call mp_bcast(nmr_macroscopic_shape, root)
    call mp_bcast(spline_ps, root)
#endif
  END SUBROUTINE gipaw_bcast_input


  !-----------------------------------------------------------------------
  ! Allocate memory for GIPAW
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_allocate
    USE becmod,        ONLY : becp
    USE lsda_mod,      ONLY : nspin, lsda
    USE pwcom
    USE ions_base,     ONLY : ntyp => nsp
    USE paw,           ONLY : paw_recon
    
    IMPLICIT NONE
    
    allocate(evq(npwx,nbnd))
    ! DO NOT not allocate becp here! it's allocated and deallocated
    ! by a lot of routines (crazy!)
    ! allocate(becp(nkb,nbnd))
    allocate(j_bare(nrxxs,3,3,nspin), b_ind_r(nrxxs,3,3), b_ind(ngm,3,3))
    allocate ( paw_recon(ntyp) )
    
  END SUBROUTINE gipaw_allocate
  

  !-----------------------------------------------------------------------
  ! Print a short summary of the calculation
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_summary
    USE io_global,     ONLY : stdout
    IMPLICIT NONE

    CALL flush_unit( stdout )
  END SUBROUTINE gipaw_summary
  

  !-----------------------------------------------------------------------
  ! Open files needed for GIPAW
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_openfil
    USE wvfct,          ONLY : npwx
    USE uspp,           ONLY : nkb
    IMPLICIT NONE

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
    call print_clock ('h_psiq')
    WRITE( stdout, * )
    call print_clock ('u_kq')
    call print_clock ('h_psi')
    WRITE( stdout, * )
    call print_clock ('apply_p')
    call print_clock ('apply_vel')
    WRITE( stdout, * )
    call print_clock ('j_para')
    call print_clock ('biot_savart')
    call print_clock ('c_sigma')
    WRITE( stdout, * )
    WRITE( stdout, * ) '     General routines'
    call print_clock ('ccalbec')
    call print_clock ('cft3')
    call print_clock ('cft3s')
    call print_clock ('cinterpolate')
    call print_clock ('davcio')
    call print_clock ('write_rec')
    WRITE( stdout, * )
#ifdef __PARA
    WRITE( stdout,  * ) '     Parallel routines'
    call print_clock ('reduce')
    call print_clock ('poolreduce')
#endif
  END SUBROUTINE print_clock_gipaw



  !-----------------------------------------------------------------------
  ! GIPAW setup
  !-----------------------------------------------------------------------
  SUBROUTINE gipaw_setup
    USE kinds,         ONLY : DP
    USE io_global,  ONLY : stdout
    USE ions_base,     ONLY : tau, nat, ntyp => nsp
    USE atom,          ONLY : nlcc
    USE wvfct,         ONLY : nbnd, et, wg, npwx
    USE lsda_mod,      ONLY : nspin, lsda
    USE scf,           ONLY : vr, vrs, vltot, rho, rho_core
    USE gvect,         ONLY : nrxx, ngm
    USE gsmooth,       ONLY : doublegrid
    USE klist,         ONLY : xk, degauss, ngauss, nks, nelec
    USE constants,     ONLY : degspin, pi
    !<apsi>
!    USE paw,           ONLY : paw_nbeta, aephi, psphi, paw_vkb, &
!                              paw_becp, paw_nkb
    USE paw,           ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp, &
                              read_recon, read_recon_paratec
    USE atom,          ONLY : r, rab
    
    !</apsi>
    
    IMPLICIT none
    integer :: ik, nt, ibnd
    logical :: nlcc_any
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
    !</apsi>
    
    call start_clock ('gipaw_setup')
    
    ! initialize pseudopotentials
    call init_us_1
    
    !<apsi>
    
    ! Read in qe format
    DO nt = 1, ntyp
       IF ( read_recon_in_paratec_fmt ) THEN
          
          ! Read in paratec format
          CALL read_recon_paratec ( file_reconstruction(nt), nt, &
               paw_recon(nt), vloc_set )
          IF ( .NOT. vloc_set .AND. job == "g-tensor" ) THEN
             CALL errore ( "gipaw_setup", &
                  "no local potential set in read_recon_paratec", 1 )
             stop
          END IF
          
       ELSE
          CALL read_recon ( file_reconstruction(nt), nt, paw_recon(nt) )
       END IF
    END DO
    
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
    
    call init_paw_1()
    
    allocate ( paw_vkb(npwx,paw_nkb) )
    allocate ( paw_becp(paw_nkb,nbnd) )
    allocate ( paw_becp2(paw_nkb,nbnd) )
    
    !<apsi>
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
                     paw_recon(nt)%psphi(il1)%psi(j)*paw_recon(nt)%psphi(il2)%psi(j))/r(j,nt)
             enddo
             
             CALL simpson( nrc, work, rab(:nrc,nt), &
                  radial_integral_diamagnetic(il1,il2,nt) )
             
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
                     / r(j,nt) ** 3
             end do
             
             call simpson( nrc, work, rab(:,nt), &
                  radial_integral_paramagnetic(il1,il2,nt) )
             if (iverbosity > 10) then
                write(stdout,*) "WWW1: ", l2, il1, il2, &
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
             CALL radial_kinetic_energy ( l2, r(:nrc,nt), &
                  paw_recon(nt)%aephi(il2)%psi(:nrc), kinetic_aephi(:nrc) )
             CALL radial_kinetic_energy ( l2, r(:nrc,nt), &
                  paw_recon(nt)%psphi(il2)%psi(:nrc), kinetic_psphi(:nrc) )
             
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * kinetic_aephi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * kinetic_psphi(j) )
             end do
             DEALLOCATE ( kinetic_aephi, kinetic_psphi )
             
             CALL simpson ( nrc, work, rab(:,nt), &
                  radial_integral_rmc(il1,il2,nt) )
             if (iverbosity > 0) then
                write(stdout,*) "WWW2: ", l2, il1, il2, &
                     radial_integral_rmc(il1,il2,nt)
             end if
             
             ALLOCATE ( aephi_dvloc_dr ( nrc ), psphi_dvloc_dr ( nrc ) )
             
             CALL radial_derivative ( r(:nrc,nt), &
                  paw_recon(nt)%gipaw_ae_vloc(:nrc), &
                  aephi_dvloc_dr(:nrc) )
             CALL radial_derivative ( r(:nrc,nt), &
                  paw_recon(nt)%gipaw_ps_vloc(:nrc), &
                  psphi_dvloc_dr ( :nrc ) )
             
             !
             ! g tensor, diamagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     * r(j,nt)
             end do
             
             call simpson( nrc, work, rab(:,nt), &
                  radial_integral_diamagnetic_so(il1,il2,nt) )
             if (iverbosity > 0) then
                write(stdout,*) "WWW3: ", l2, il1, il2, &
                     radial_integral_diamagnetic_so(il1,il2,nt) * alpha
             end if
             
             !
             ! g tensor, paramagnetic
             !
             do j = 1, nrc
                work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) &
                     * paw_recon(nt)%aephi(il2)%psi(j) - paw_recon(nt)%psphi(il1)%psi(j) &
                     * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / r(j,nt)
             end do
             if (iverbosity > 10) then
                if ( l1 == 0 ) then
                   do j = 1, nrc
                      write(90,*) r(j,nt), work(j)*r(j,nt)**2
                   end do
                   write(90,*) ""
                end if
             end if
             
             call simpson( nrc,work,rab(:,nt), &
                  radial_integral_paramagnetic_so(il1,il2,nt) )
             if (iverbosity > 0) then
                write(stdout,*) "WWW4: ", l2, il1, il2, &
                     radial_integral_paramagnetic_so(il1,il2,nt) * alpha
             end if
             
             DEALLOCATE ( aephi_dvloc_dr, psphi_dvloc_dr )
             
          enddo
          
          DEALLOCATE ( work )
          
       enddo
    enddo
    
    !</apsi>
    
    ! terms for paramagnetic
    
    !<apsi>
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
    
    if (iverbosity > 10) then
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
       write(stdout,'(A,9F8.4)') "DDD1x: ", mysum1(1,:)
       write(stdout,'(A,9F8.4)') "DDD1y: ", mysum1(2,:)
       write(stdout,'(A,9F8.4)') "DDD1z: ", mysum1(3,:)
       write(stdout,'(A,9F8.4)') "DDD2x: ", mysum2(1,:)
       write(stdout,'(A,9F8.4)') "DDD2y: ", mysum2(2,:)
       write(stdout,'(A,9F8.4)') "DDD2z: ", mysum2(3,:)
    end if
    
    deallocate ( lm2l, lm2m )
    
    ! check whether the symmetry operations map the coordinate axis to each
    ! other - if not, remove them (how to check the k point mesh then? - oops!
    
    
    ! computes the total local potential (external+scf) on the smooth grid
    call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)

    ! compute the D for the pseudopotentials
    call newd

    ! set non linear core correction stuff
    nlcc_any = .false.
    do nt = 1, ntyp
       nlcc_any = nlcc_any.or.nlcc (nt)
    enddo
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
#ifdef __PARA
    ! find the minimum across pools
    call poolextreme (emin, -1)
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
#ifdef __PARA
      ! find the maximum across pools
      call poolextreme (emax, + 1)
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
  
 
!-----------------------------------------------------------------------
END MODULE gipaw_module
!-----------------------------------------------------------------------
