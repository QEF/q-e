!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE nmr_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for NMR calculations
  !
  USE kinds, ONLY : DP
  USE parameters, ONLY: npk
  IMPLICIT NONE
  SAVE
  
  ! speed of light in atomic units: c = 1/alpha
  REAL(DP), PARAMETER :: c = 137.03599911d0

  ! avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142d23

  ! bohrradius to cm
  REAL(DP), PARAMETER :: a0_to_cm = 0.529177d-8

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
  REAL(DP) :: q_nmr

  ! verbosity
  INTEGER :: iverbosity

  ! job: nmr, g_tensor
  CHARACTER(80) :: job

  CONTAINS
  !-----------------------------------------------------------------------
  ! Read in the nmr input file.
  ! Format: &inputmagn
  !              prefix = '...'     prefix of SCF calculation
  !              tmp_dir = '...'    scratch directory
  !              job = 'nmr' or 'g_tensor'
  !              conv_threshold = 1d-14
  !              q_nmr = 0.01d0
  !              iverbosity = 0
  !         /
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_readin()
    USE io_files,      ONLY : nd_nmbr, prefix, tmp_dir  
    USE io_global,     ONLY : ionode
    IMPLICIT NONE
    INTEGER :: ios
    NAMELIST /inputmagn/ job, prefix, tmp_dir, conv_threshold, &
                         q_nmr, iverbosity
   
    if ( .not. ionode ) goto 400

    job = ''
    prefix = 'pwscf'
    tmp_dir = './scratch/'
    conv_threshold = 1d-14
    q_nmr = 0.01d0
    iverbosity = 0
    read( 5, inputmagn, err = 200, iostat = ios )
200 call errore( 'nmr_readin', 'reading inputmagn namelist', abs( ios ) )
400 continue
    call nmr_bcast_input
  END SUBROUTINE nmr_readin



  !-----------------------------------------------------------------------
  ! Broadcast input data to all processors 
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_bcast_input
#ifdef __PARA
#include "f_defs.h"
    USE mp,            ONLY : mp_bcast
    USE io_files,      ONLY : prefix, tmp_dir
    implicit none
    integer :: root = 0
    call mp_bcast(job, root)
    call mp_bcast(prefix, root)
    call mp_bcast(tmp_dir, root)
    call mp_bcast(conv_threshold, root)
    call mp_bcast(q_nmr, root)
    call mp_bcast(iverbosity, root)
#endif
  END SUBROUTINE nmr_bcast_input



  !-----------------------------------------------------------------------
  ! Allocate memory for NMR
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_allocate
    USE becmod, ONLY : becp
    USE lsda_mod,      ONLY : nspin, lsda
    USE pwcom
    IMPLICIT NONE
  
    allocate(evq(npwx,nbnd))
    allocate(becp(nkb,nbnd))
    allocate(j_bare(nrxxs,3,3,nspin), b_ind_r(nrxxs,3,3), b_ind(ngm,3,3))
  END SUBROUTINE nmr_allocate
 


  !-----------------------------------------------------------------------
  ! Print a short summary of the calculation
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_summary
    USE io_global,     ONLY : stdout
    IMPLICIT NONE

    CALL flush_unit( stdout )
  END SUBROUTINE nmr_summary



  !-----------------------------------------------------------------------
  ! Open files needed for NMR
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_openfil
    USE wvfct,          ONLY : npwx
    USE uspp,           ONLY : nkb
    IMPLICIT NONE

  END SUBROUTINE nmr_openfil



  !-----------------------------------------------------------------------
  ! Print timings
  !-----------------------------------------------------------------------
  SUBROUTINE print_clock_nmr
    USE io_global,  ONLY : stdout
    IMPLICIT NONE

    WRITE( stdout, * )
    call print_clock ('NMR') 
    WRITE( stdout,  * ) '    INITIALIZATION: '
    call print_clock ('nmr_setup')
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
  END SUBROUTINE print_clock_nmr



  !-----------------------------------------------------------------------
  ! NMR setup
  !-----------------------------------------------------------------------
  SUBROUTINE nmr_setup
    USE kinds,         ONLY : DP
    USE ions_base,     ONLY : tau, nat, ntyp => nsp
    USE atom,          ONLY : nlcc
    USE wvfct,         ONLY : nbnd, et, wg
    USE lsda_mod,      ONLY : nspin, lsda
    USE scf,           ONLY : vr, vrs, vltot, rho, rho_core
    USE gvect,         ONLY : nrxx, ngm
    USE gsmooth,       ONLY : doublegrid
    USE klist,         ONLY : xk, degauss, ngauss, nks, nelec
    USE constants,     ONLY : degspin, pi
    
    IMPLICIT none
    integer :: ik, nt, ibnd
    logical :: nlcc_any
    real(dp) :: emin, emax

    call start_clock ('nmr_setup')

    ! initialize pseudopotentials
    call init_us_1

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
    if (degauss.ne.0.d0) then
      call errore('nmr_setup', 'implemented only for insulators', -1)
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
      alpha_pv = 2.d0 * (emax - emin)
    endif
    ! avoid zero value for alpha_pv
    alpha_pv = max (alpha_pv, 1.0d-2)

    call stop_clock('nmr_setup')
  END SUBROUTINE nmr_setup
 
!-----------------------------------------------------------------------
END MODULE nmr_module
!-----------------------------------------------------------------------
