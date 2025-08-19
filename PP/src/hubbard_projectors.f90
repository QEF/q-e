!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
!
!-----------------------------------------------------------------------
SUBROUTINE hubbard_projectors (filplot, plot_num, nc, n0, nplots)
  !-----------------------------------------------------------------------
  !
  ! This routine reads the Bloch sums of WFs at k points in the G-space 
  ! and then builds the localized WFs in r-space.
  ! The visualization of the localized Wfs is done in the supercell of size 
  ! nc(1) x nc(2) x nc(3). These localized WFs are shifted by n0(1), n0(2), n0(3) 
  ! in units of the real-space lattice vectors.
  !
  ! Inspired by PP/src/wannier_plot.f90 and Modules/compute_dipole.f90
  ! However, here the algorithm for working with the r-grids is redesigned 
  ! compared to PP/src/wannier_plot.f90, because the latter does not provide 
  ! the correct results.
  !
  ! Currently, this routine does not parallelize the calculations over 
  ! r-grid points, hence the calculation can be slow for dense r and k grids.
  !
  ! Note: custom_data_structure is not really needed since we do not 
  ! distribute the r-grid calculations over multiple processors
  ! (that would require the call of gather_grid to collect the results).
  ! If the distribution is done in the future, then custom_data_structure
  ! would be needed for the custom FFT.  
  !
  ! Written by I. Timrov (November 2024)
  ! 
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : at, bg, celldm, ibrav
  USE io_global,        ONLY : stdout, ionode
  USE ions_base,        ONLY : nat, ntyp => nsp, atm, ityp, zv, tau
  USE constants,        ONLY : tpi
  USE gvect,            ONLY : gcutm
  USE gvecs,            ONLY : dual
  USE gvecw,            ONLY : ecutwfc
  USE io_files,         ONLY : diropn, iunhub_noS, nwordwfcU, prefix
  USE run_info,         ONLY : title
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_wave,         ONLY : invfft_wave
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : inter_pool_comm
  USE klist,            ONLY : nks, nkstot, ngk, igk_k, xk, wk
  USE ldaU,             ONLY : nwfcU, wfcU
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : npol, noncolin
  USE scatter_mod,      ONLY : gather_grid
  USE lsda_mod,         ONLY : isk, nspin
  USE symm_base,        ONLY : nsym
  USE control_flags,    ONLY : gamma_only
  USE wannier,          ONLY : dfftp_c, dffts_c
  !
  IMPLICIT NONE
  CHARACTER(len=256), INTENT(IN) :: filplot
  INTEGER, INTENT(IN) :: plot_num
  INTEGER, INTENT(OUT) :: nplots
  INTEGER, INTENT(IN) :: nc(3), & ! dimension of the supercell 
                         n0(3)    ! position of the center
  COMPLEX(DP), ALLOCATABLE :: psic(:),            &
                              psic_all(:),        &
                              psic3(:,:,:),       &
                              psic3_sum(:,:,:,:,:)
  REAL(DP), ALLOCATABLE ::    rho(:)
  REAL(DP) :: arg
  REAL(DP) :: r(3)
  CHARACTER(len=256) :: filplot2, spin_label(2)
  COMPLEX(DP) :: phase
  LOGICAL :: exst
  INTEGER :: ik, i, j, k, m, ir, i1, j1, k1, is, nkstot_
  !
  nwordwfcU = nwfcU * npwx * npol
  !
  nplots = nwfcU
  !
  ! The results with symmetry differ a bit from the one without symmetry.
  ! Maybe the symmetrization of rho is needed?
  IF (nsym > 1) CALL errore ('hubbard_projectors', &
             & 'Symmetry is currently not supported',1)
  !
  IF (noncolin) CALL errore ('hubbard_projectors', &
             & 'Noncollinear case is currently not supported',1)
  !
  IF (nspin==1) THEN
     nkstot_ = nkstot
  ELSEIF (nspin==2) THEN
     nkstot_ = nkstot/2
  ENDIF
  IF (nc(1)*nc(2)*nc(3) < nkstot_) THEN
     WRITE(stdout,'(/5x,A)') 'Warning!!! nc(1)*nc(2)*nc(3) is smaller than the k-grid size nk1*nk2*nk3'
     WRITE(stdout,'(5x,A)') 'This will lead to the wrong visualization of the WFs.'
     ! Explanation: From tests it seems that the R-grid and the k-grid must match, otherwise the 
     ! visualized WFs look wrong and weird.
     ! 
  ENDIF
  !
  WRITE(stdout,'(/5x,A)') 'Plotting the squared modulus of the Hubbard projector functions'
  !
  spin_label(:) = ''
  IF (nspin==2) THEN
     spin_label(1) = '_up'
     spin_label(2) = '_down'
  ENDIF
  !
  ! Define the data structure for the custom FFT
  CALL custom_data_structure (gamma_only, nc)   
  !
  ALLOCATE (wfcU(npwx*npol,nwfcU))
  ALLOCATE (psic(dffts%nnr))
  ALLOCATE (psic_all(dffts%nr1x * dffts%nr2x * dffts%nr3x))
  ALLOCATE (psic3(dffts%nr1x, dffts%nr2x, dffts%nr3x))
  ALLOCATE (psic3_sum(dffts_c%nr1x, dffts_c%nr2x, dffts_c%nr3x, nwfcU, nspin))
  ALLOCATE (rho(dffts_c%nr1x * dffts_c%nr2x * dffts_c%nr3x))
  !
  psic3 = ZERO
  psic3_sum = ZERO
  !
  ! Open files to store the new Hubbard projectors (WFs).
  CALL diropn( iunhub_noS, 'hubnoS', 2*nwordwfcU, exst )
  IF (.NOT.exst ) CALL errore ('hubbard_projectors', &
             & 'File '//trim( prefix )//'.hubnoS'//' not found',1)
  !
  WRITE(stdout,'(/5x,a,I5)') 'Number of local k points = ', nks
  DO ik = 1, nks
     !
     WRITE (stdout,'(i8)',advance='no') ik
     IF ( MOD(ik,10) == 0 ) WRITE (stdout,*)
     FLUSH(stdout)
     !
     is  = isk(ik)
     !
     ! Read Bloch sums of WFs in G-space at k point ik
     CALL davcio (wfcU, 2*nwordwfcU, iunhub_noS, ik, -1)
     !
     DO m = 1, nwfcU
        !
        ! Inverse FFT: G-space -> r-space
        CALL invfft_wave (npwx, ngk(ik), igk_k(:,ik), wfcU(:,m), psic)
        !
#if defined(__MPI)
        CALL gather_grid (dffts, psic, psic_all)
#else
        psic_all = psic
#endif
        ! 
        DO k = 1, dffts%nr3x
           DO j = 1, dffts%nr2x
              DO i = 1, dffts%nr1x
                 ir = i + (j-1)*dffts%nr1x + &
                          (k-1)*dffts%nr2x*dffts%nr1x
                 psic3(i,j,k) = psic_all(ir)
              ENDDO
           ENDDO
        ENDDO
        !
        DO k1 = 1, dffts_c%nr3x
           DO j1 = 1, dffts_c%nr2x
              DO i1 = 1, dffts_c%nr1x
                 !
                 r = n0(1)*at(:,1) + &
                     n0(2)*at(:,2) + &
                     n0(3)*at(:,3)
                 !
                 ! Compute the real-space position vector "r" for the larger grid
                 r = r + DBLE(i1-1)*at(:,1)*nc(1) / DBLE(dffts_c%nr1x) + &
                         DBLE(j1-1)*at(:,2)*nc(2) / DBLE(dffts_c%nr2x) + &
                         DBLE(k1-1)*at(:,3)*nc(3) / DBLE(dffts_c%nr3x)
                 !
                 ! Compute the phase factor exp(ik*r)
                 arg = ( xk(1,ik) * r(1) + &
                         xk(2,ik) * r(2) + &
                         xk(3,ik) * r(3) ) * tpi
                 phase = CMPLX(cos(arg), sin(arg), kind=DP)
                 !
                 ! Map the larger grid indices (i1,j1,k1) to the smaller grid indices (i,j,k) 
                 i = MOD(i1-1, dffts%nr1x) + 1
                 j = MOD(j1-1, dffts%nr2x) + 1
                 k = MOD(k1-1, dffts%nr3x) + 1
                 !
                 psic3_sum(i1,j1,k1,m,is) = psic3_sum(i1,j1,k1,m,is) + &
                                            CMPLX(wk(ik),0.d0,kind=DP) &
                                            * psic3(i,j,k) * phase
                 !
              ENDDO
           ENDDO
        ENDDO
        !
     ENDDO ! m
     !
  ENDDO ! ik
  !
  ! Collect the data among all k-points pools
  CALL mp_sum (psic3_sum, inter_pool_comm)
  !
  WRITE(stdout,*) ""
  !
  IF ( ionode ) THEN
     DO is = 1, nspin
        DO m = 1, nwfcU
           !
           DO i1 = 1, dffts_c%nr1x
              DO j1 = 1, dffts_c%nr2x
                 DO k1 = 1, dffts_c%nr3x
                    !
                    ir = i1 + (j1-1)*(dffts_c%nr1x) + &
                              (k1-1)*(dffts_c%nr2x)*(dffts_c%nr1x)
                    ! 
                    rho(ir) = DBLE(psic3_sum(i1,j1,k1,m,is))**2 + &
                              AIMAG(psic3_sum(i1,j1,k1,m,is))**2
                 ENDDO
              ENDDO
           ENDDO
           !
           WRITE(filplot2,'(A, I0.3, A)') TRIM(filplot), m, TRIM(spin_label(is))
           !
           CALL plot_io (TRIM(filplot2), title, dffts_c%nr1x, dffts_c%nr2x, dffts_c%nr3x, &
                dffts_c%nr1, dffts_c%nr2, dffts_c%nr3, nat, ntyp, ibrav, celldm, at, &
                gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, rho, + 1)
           !
        ENDDO
     ENDDO
  ENDIF
  !
  CLOSE ( UNIT=iunhub_noS, STATUS='KEEP' )
  DEALLOCATE (wfcU)
  DEALLOCATE (psic)
  DEALLOCATE (psic_all)
  DEALLOCATE (psic3)
  DEALLOCATE (psic3_sum)
  DEALLOCATE (rho)
  !
  RETURN
  !
CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE custom_data_structure (gamma_only_, nc_)
  !-----------------------------------------------------------------------
  !! This routine sets the data structure for the custom FFT arrays (both the
  !! smooth and the dense grid).  
  !! In the parallel case, it distributes columns to processes too.  
  !! BEWARE: to compute \(\text{gkcut}\), \(\text{nks}\) and the list of
  !! k-points or \(\text{nks}=0\) and the primitive lattice vectors, 
  !! \(\text{bg}\) are needed.
  !! Adapted from PW/src/data_structure.f90
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : nproc_bgrp, intra_bgrp_comm, nyfft, ntask_groups
  USE mp_pools,   ONLY : inter_pool_comm
  USE fft_base,   ONLY : dfftp, dffts, fft_base_info, smap
  USE fft_types,  ONLY : fft_type_init
  USE cell_base,  ONLY : at, bg, tpiba
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : gcutm, gvect_init, deallocate_gvect
  USE gvecs,      ONLY : gcutms, gvecs_init, doublegrid
  USE gvecw,      ONLY : gcutw, gkcut
  USE io_global,  ONLY : stdout, ionode
  ! FIXME: find a better way to transmit these three variables, or remove them
  USE realus,     ONLY : real_space
  USE symm_base,  ONLY : fft_fact
  USE command_line_options, ONLY: pencil_decomposition_
  USE command_line_options, ONLY : nmany_
  USE wannier,              ONLY : dfftp_c, dffts_c
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: gamma_only_
  INTEGER, INTENT(IN) :: nc_(3)
  !
  ! ... local variables
  !
  INTEGER :: ik, i, ngm_, ngs_
  LOGICAL :: lpara
  REAL(DP) :: at_(3,3), bg_(3,3)
  !
  lpara = ( nproc_bgrp > 1 )
  !
  DO i = 1, 3
     at_(:,i) = at(:,i) * nc_(i)
     bg_(:,i) = bg(:,i) / nc_(i)
  ENDDO
  !
  ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
  !
  IF (nks == 0) THEN
     !
     ! k-point list not available:
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * MAX( SQRT( SUM(bg(1:3,1)**2) ), &
                          SQRT( SUM(bg(1:3,2)**2) ), &
                          SQRT( SUM(bg(1:3,3)**2) )  )
  ELSE
     gkcut = 0.0d0
     DO ik = 1, nks
        gkcut = MAX(gkcut, SQRT( SUM(xk(1:3,ik)**2) ) )
     ENDDO
  ENDIF
  gkcut = (SQRT(gcutw) + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max( gkcut, inter_pool_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  ! task group are disabled if real_space calculation of calbec is used
  dffts_c%has_task_groups = (ntask_groups >1) .AND. .NOT. real_space
  !
  ! The custom FFT dimensions must exactly match the original
  ! FFT dimensions with an integer multiplicative factors nc(:).
  !
  dffts_c%nr1x = nc(1) * dffts%nr1x
  dffts_c%nr2x = nc(2) * dffts%nr2x
  dffts_c%nr3x = nc(3) * dffts%nr3x
  dffts_c%nr1  = nc(1) * dffts%nr1
  dffts_c%nr2  = nc(2) * dffts%nr2
  dffts_c%nr3  = nc(3) * dffts%nr3
  dffts_c%nnr  = nc(1) * nc(2) * nc(3) * dffts%nnr
  !
  dfftp_c%nr1x = nc(1) * dfftp%nr1x
  dfftp_c%nr2x = nc(2) * dfftp%nr2x
  dfftp_c%nr3x = nc(3) * dfftp%nr3x
  dfftp_c%nr1  = nc(1) * dfftp%nr1
  dfftp_c%nr2  = nc(2) * dfftp%nr2
  dfftp_c%nr3  = nc(3) * dfftp%nr3
  dfftp_c%nnr  = nc(1) * nc(2) * nc(3) * dfftp%nnr
  !
  CALL fft_type_init( dffts_c, smap, "wave", gamma_only_, lpara, intra_bgrp_comm, &
       at_, bg_, gkcut, gcutms/gkcut, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_, use_pd=pencil_decomposition_  )
  !
  CALL fft_type_init( dfftp_c, smap, "rho" , gamma_only_, lpara, intra_bgrp_comm, &
       at_, bg_, gcutm , 4.d0, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_, use_pd=pencil_decomposition_ )
  !
  ! define the clock labels ( this enables the corresponding fft too ! )
  dffts_c%rho_clock_label = 'ffts' ; dffts_c%wave_clock_label = 'fftw'
  dfftp_c%rho_clock_label = 'fft'
  ! this makes so that interpolation is just a copy.
  IF (.NOT.doublegrid) dfftp_c%grid_id = dffts_c%grid_id
  !
  !CALL fft_base_info( ionode, stdout )
  ngs_ = dffts_c%ngl( dffts_c%mype + 1 )
  ngm_ = dfftp_c%ngl( dfftp_c%mype + 1 )
  !
  IF( gamma_only_ ) THEN
     ngs_ = (ngs_ + 1)/2
     ngm_ = (ngm_ + 1)/2
  ENDIF
  !
  ! ... on output, ngm_ and ngs_ contain the local number of G-vectors
  ! for the two grids. Re-initialize local and global number of G-vectors
  !
  CALL deallocate_gvect()
  !
  CALL gvect_init( ngm_, intra_bgrp_comm )
  CALL gvecs_init( ngs_, intra_bgrp_comm )
  !
  RETURN
  !
END SUBROUTINE custom_data_structure

END SUBROUTINE hubbard_projectors
