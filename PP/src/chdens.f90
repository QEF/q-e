!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE chdens (filplot,plot_num)
  !-----------------------------------------------------------------------
  !      Writes the charge density (or potential, or polarisation)
  !      into a file format suitable for plotting
  !-----------------------------------------------------------------------
  !
  !      DESCRIPTION of the INPUT: see file INPUT_PP in Doc/
  !
  USE kinds,      ONLY : dp
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : nd_nmbr
  USE mp_pools,   ONLY : nproc_pool
  USE mp_world,   ONLY : world_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_bcast
  USE parameters, ONLY : ntypx
  USE constants,  ONLY : pi, fpi
  USE cell_base,  ONLY : at, bg, celldm, ibrav, alat, omega, tpiba, tpiba2
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE lsda_mod,   ONLY : nspin
  USE fft_base,   ONLY : dfftp, dffts
  USE scatter_mod,   ONLY : scatter_grid
  USE fft_interfaces,  ONLY : fwfft
  USE fft_types,  ONLY : fft_type_allocate
  USE gvect,      ONLY : ngm, nl, g, gcutm
  USE gvecs,      ONLY : gcutms, doublegrid, dual, ecuts 
  USE recvec_subs,ONLY: ggen 
  USE gvecw,      ONLY: ecutwfc
  USE run_info,   ONLY: title
  USE control_flags, ONLY: gamma_only
  USE wavefunctions_module,  ONLY: psic

  IMPLICIT NONE
  CHARACTER (len=256), INTENT(in) :: filplot
  !
  ! If plot_num=-1 the dimensions and structural data are read from the charge
  ! or potential file, otherwise it uses the data already read from
  ! the files in outdir.
  !
  INTEGER, INTENT(in) :: plot_num
  !
  INTEGER, PARAMETER :: nfilemax = 7
  ! maximum number of files with charge

  INTEGER :: ounit, iflag, ios, ipol, nfile, ifile, nx, ny, nz, &
       na, i, output_format, idum, direction

  real(DP) :: e1(3), e2(3), e3(3), x0 (3), radius, m1, m2, m3, &
       weight (nfilemax), isovalue,heightmin,heightmax

  real(DP), ALLOCATABLE :: aux(:)

  CHARACTER (len=256) :: fileout
  CHARACTER (len=13), DIMENSION(0:7) :: formatname = &
       (/ 'gnuplot      ', &
          'contour.x    ', &
          'plotrho.x    ', &
          'XCrySDen     ', &
          'gOpenMol     ', &
          'XCrySDen     ', &
          'Gaussian cube', & 
          'gnuplot x,y,f' /)
  CHARACTER (len=20), DIMENSION(0:4) :: plotname = &
       (/ '1D spherical average', &
          '1D along a line     ', &
          '2D contour          ', &
          '3D                  ', &
          '2D polar on a sphere'/)

  real(DP) :: celldms (6), gcutmsa, duals, zvs(ntypx), ats(3,3)
  real(DP), ALLOCATABLE :: taus (:,:), rhor(:), rhos(:)
  INTEGER :: ibravs, nr1sxa, nr2sxa, nr3sxa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  INTEGER, ALLOCATABLE :: ityps (:)
  CHARACTER (len=3) :: atms(ntypx)
  CHARACTER (len=256) :: filepp(nfilemax)
  CHARACTER (len=20) :: interpolation
  real(DP) :: rhotot
  COMPLEX(DP), ALLOCATABLE:: rhog (:)
  ! rho or polarization in G space
  LOGICAL :: fast3d, isostm_flag

  NAMELIST /plot/  &
       nfile, filepp, weight, iflag, e1, e2, e3, nx, ny, nz, x0, &
       radius, output_format, fileout, interpolation, &
       isostm_flag, isovalue, heightmin, heightmax, direction

  !
  !   set the DEFAULT values
  !
  nfile         = 1
  filepp(1)     = filplot
  weight(1)     = 1.0d0
  iflag         = 0
  radius        = 1.0d0
  output_format = -1
  fileout       = ' '
  e1(:)         = 0.d0
  e2(:)         = 0.d0
  e3(:)         = 0.d0
  x0(:)         = 0.d0
  nx            = 0
  ny            = 0
  nz            = 0
  interpolation = 'fourier'
  isostm_flag   = .false.
  isovalue      = 0.d0
  heightmin     = 0.0d0
  heightmax     = 1.0d0
  direction     = 1
  !
  !    read and check input data
  !
  ! reading the namelist 'plot'
  !
  IF (ionode) READ (5, plot, iostat = ios)
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  CALL mp_bcast( nfile, ionode_id, world_comm )

  IF (ios /= 0) THEN
     IF (nfile > nfilemax) THEN
        ! if this happens the reading of the namelist will fail
        ! tell to user why
        CALL infomsg('chdens ', 'nfile is too large, exiting')
     ELSE
        CALL infomsg ('chdens', 'namelist plot not found or invalid, exiting')
     ENDIF
     RETURN
  ENDIF

  CALL mp_bcast( filepp, ionode_id, world_comm )
  CALL mp_bcast( weight, ionode_id, world_comm )
  CALL mp_bcast( iflag, ionode_id, world_comm )
  CALL mp_bcast( radius, ionode_id, world_comm )
  CALL mp_bcast( output_format, ionode_id, world_comm )
  CALL mp_bcast( fileout, ionode_id, world_comm )
  CALL mp_bcast( e1, ionode_id, world_comm )
  CALL mp_bcast( e2, ionode_id, world_comm )
  CALL mp_bcast( e3, ionode_id, world_comm )
  CALL mp_bcast( x0, ionode_id, world_comm )
  CALL mp_bcast( nx, ionode_id, world_comm )
  CALL mp_bcast( ny, ionode_id, world_comm )
  CALL mp_bcast( nz, ionode_id, world_comm )
  CALL mp_bcast( interpolation, ionode_id, world_comm )
  CALL mp_bcast( isostm_flag, ionode_id, world_comm )
  CALL mp_bcast( isovalue, ionode_id, world_comm )
  CALL mp_bcast( heightmin, ionode_id, world_comm )
  CALL mp_bcast( heightmax, ionode_id, world_comm )
  CALL mp_bcast( direction, ionode_id, world_comm )  

  IF (output_format == -1 .or. iflag == -1) THEN
     CALL infomsg ('chdens', 'output format not set, exiting' )
     RETURN
  ENDIF
  !
  ! check for number of files
  !
  IF (nfile < 1 .or. nfile > nfilemax) &
       CALL errore ('chdens ', 'nfile is wrong ', 1)

  ! check for iflag

  IF (iflag <= 1) THEN

     ! 1D plot : check variables

     IF (e1(1)**2 + e1(2)**2 + e1(3)**2 < 1d-6) &
         CALL errore ('chdens', 'missing e1 vector', 1)
     IF (nx <= 0 )   CALL errore ('chdens', 'wrong nx', 1)

  ELSEIF (iflag == 2) THEN

     ! 2D plot : check variables

     IF (e1(1)**2 + e1(2)**2 + e1(3)**2 <  1d-6 .or. &
         e2(1)**2 + e2(2)**2 + e2(3)**2 <  1d-6)     &
         CALL errore ('chdens', 'missing e1/e2 vectors', 1)
     IF (abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6) &
         CALL errore ('chdens', 'e1 and e2 are not orthogonal', 1)
     IF (nx <= 0 .or. ny <= 0 )   CALL errore ('chdens', 'wrong nx/ny', 2)

  ELSEIF (iflag == 3) THEN

     ! 3D plot : check variables

     IF ( abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6 .or. &
          abs(e1(1)*e3(1) + e1(2)*e3(2) + e1(3)*e3(3)) > 1d-6 .or. &
          abs(e2(1)*e3(1) + e2(2)*e3(2) + e2(3)*e3(3)) > 1d-6 )    &
         CALL errore ('chdens', 'e1, e2, e3 are not orthogonal', 1)

     IF ((iflag==3) .and.(output_format < 3 .or. output_format > 6)) &
        CALL errore ('chdens', 'incompatible iflag/output_format', 1)
     IF ((iflag/=3) .and. ((output_format == 5) .or. (output_format == 6))) &
        CALL errore ('chdens', 'output_format=5/6, iflag<>3', 1)

  ELSEIF (iflag  == 4) THEN

     IF (nx <= 0 .or. ny <= 0 )   CALL errore ('chdens', 'wrong nx/ny', 4)

  ELSE

     CALL errore ('chdens', 'iflag not implemented', 1)

  ENDIF

  ! check interpolation
  if (trim(interpolation) /= 'fourier' .and. trim(interpolation) /= 'bspline') &
     call errore('chdens', 'wrong interpolation: ' // trim(interpolation), 1)

  ! if isostm_flag checks whether the input variables are set
  IF (isostm_flag) THEN
     IF (heightmax > 1.0 .or. heightmin > 1.0 .or. heightmin < 0.0 &
               .or. heightmax < 0.0 ) THEN
         CALL errore('isostm','problem with heightmax/min',1)
     ENDIF
      
     IF (direction /= 1 .and. direction /= -1) THEN
         CALL errore('isostm','direction not equal to +- 1',1)
     ENDIF
  END IF

  !
  ! Read the header and allocate objects
  !
  IF (plot_num==-1) THEN
     IF (ionode) &
        CALL read_io_header(filepp (1), title, dfftp%nr1x, dfftp%nr2x, &
                dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp,&
                ibrav, celldm, at, gcutm, dual, ecutwfc, idum )
     CALL mp_bcast( title, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr1x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr2x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr3x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr1, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr2, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr3, ionode_id, world_comm )
     CALL mp_bcast( nat, ionode_id, world_comm )
     CALL mp_bcast( ntyp, ionode_id, world_comm )
     CALL mp_bcast( ibrav, ionode_id, world_comm )
     CALL mp_bcast( celldm, ionode_id, world_comm )
     CALL mp_bcast( at, ionode_id, world_comm )
     CALL mp_bcast( gcutm, ionode_id, world_comm )
     CALL mp_bcast( dual, ionode_id, world_comm )
     CALL mp_bcast( ecutwfc, ionode_id, world_comm )
     !
     ! ... see comment above
     !
     ALLOCATE(tau (3, nat))
     ALLOCATE(ityp(nat))
     !
     CALL latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm (1) ! define alat
     at = at / alat    ! bring at in units of alat

     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2
     doublegrid = dual>4.0d0
     IF (doublegrid) THEN
        gcutms = 4.d0 * ecutwfc / tpiba2
     ELSE
        gcutms = gcutm
     ENDIF

     nspin = 1

     CALL recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     CALL volume (alat, at(1,1), at(1,2), at(1,3), omega)
     CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm )
     CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm)
  ENDIF

  ALLOCATE  (rhor(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  ALLOCATE  (rhos(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  ALLOCATE  (taus( 3 , nat))
  ALLOCATE  (ityps( nat))
  !
  rhor (:) = 0.0_DP
  !
  ! Read files, verify consistency
  ! Note that only rho is read; all other quantities are discarded
  !
  DO ifile = 1, nfile
     !
     CALL plot_io (filepp (ifile), title, nr1sxa, nr2sxa, nr3sxa, &
          nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
          duals, ecuts, idum, atms, ityps, zvs, taus, rhos, - 1)

     IF (ifile==1.and.plot_num==-1) THEN
        atm=atms
        ityp=ityps
        zv=zvs
        tau=taus
     ENDIF
     !
     IF (nats>nat) CALL errore ('chdens', 'wrong file order? ', 1)
     IF (dfftp%nr1x/=nr1sxa.or.dfftp%nr2x/=nr2sxa) CALL &
          errore ('chdens', 'incompatible nr1x or nr2x', 1)
     IF (dfftp%nr1/=nr1sa.or.dfftp%nr2/=nr2sa.or.dfftp%nr3/=nr3sa) CALL &
          errore ('chdens', 'incompatible nr1 or nr2 or nr3', 1)
     IF (ibravs/=ibrav) CALL errore ('chdens', 'incompatible ibrav', 1)
     IF (abs(gcutmsa-gcutm)>1.d-8.or.abs(duals-dual)>1.d-8.or.&
         abs(ecuts-ecutwfc)>1.d-8) &
          CALL errore ('chdens', 'incompatible gcutm or dual or ecut', 1)
     IF (ibravs /= 0 ) THEN
        DO i = 1, 6
           IF (abs( celldm (i)-celldms (i) ) > 1.0d-7 ) &
              CALL errore ('chdens', 'incompatible celldm', 1)
        ENDDO
     ENDIF
     !
     rhor (:) = rhor (:) + weight (ifile) * rhos (:)
  ENDDO
  DEALLOCATE (ityps)
  DEALLOCATE (taus)
  DEALLOCATE (rhos)
  !
  ! open output file, i.e., "fileout"
  !
  IF (ionode) THEN
     IF (fileout /= ' ') THEN
        ounit = 1
        OPEN (unit=ounit, file=fileout, form='formatted', status='unknown')
        WRITE( stdout, '(/5x,"Writing data to be plotted to file ",a)') &
             trim(fileout)
     ELSE
        ounit = 6
     ENDIF
  ENDIF

  ! the isostm subroutine is called only when isostm_flag is true and the
  ! charge density is related to an STM image (5) or is read from a file 
  IF ( (isostm_flag) .AND. ( (plot_num == -1) .OR. (plot_num == 5) ) ) THEN
     IF ( .NOT. (iflag == 2))&
        CALL errore ('chdens', 'isostm should have iflag = 2', 1)
     CALL isostm_plot(rhor, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
           isovalue, heightmin, heightmax, direction)     
  END IF

  
  !
  !    At this point we start the calculations, first we normalize the
  !    vectors defining the plotting region.
  !    If these vectors have 0 length, replace them with crystal axis
  !

  m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  IF (abs(m1) < 1.d-6) THEN
     e1 (:) = at(:,1)
     m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  ENDIF
  e1 (:) = e1 (:) / m1
  !
  m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  IF (abs(m2) < 1.d-6) THEN
     e2 (:) = at(:,2)
     m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  ENDIF
  e2 (:) = e2 (:) / m2
  !
  m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  IF (abs(m3) < 1.d-6) THEN
     e3 (:) = at(:,3)
     m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  ENDIF
  e3 (:) = e3 (:) / m3
  !
  ! are vectors defining the plotting region aligned along xyz ?
  !
  fast3d = ( e1(2) == 0.d0  .and.  e1(3) == 0.d0) .and. &
           ( e2(1) == 0.d0  .and.  e2(3) == 0.d0) .and. &
           ( e3(1) == 0.d0  .and.  e3(2) == 0.d0)
  !
  ! are crystal axis aligned along xyz ?
  !
  fast3d = fast3d .and. &
       ( at(2,1) == 0.d0  .and.  at(3,1) == 0.d0) .and. &
       ( at(1,2) == 0.d0  .and.  at(3,2) == 0.d0) .and. &
       ( at(1,3) == 0.d0  .and.  at(2,3) == 0.d0)

  fast3d = fast3d .and. (trim(interpolation) == 'fourier')
  !
  !    Initialise FFT for rho(r) => rho(G) conversion if needed
  !
  IF (.not. ( iflag == 3 .and. ( output_format == 5 .or. &
                                 output_format == 6 .or. &
                                 fast3d ) ) ) THEN
     IF (plot_num==-1) THEN
        !
        gamma_only=.false.
!       nproc_pool=1
        !
        CALL data_structure ( gamma_only )
        CALL allocate_fft()
        !
        !    and rebuild G-vectors in reciprocal space
        !
        CALL ggen ( gamma_only, at, bg )
        !
        !    here we compute the fourier components of the quantity to plot
        !
     ELSE
        !
        IF (gamma_only .and. (trim(interpolation) == 'fourier')) THEN
             WRITE(stdout,'(/"BEWARE: plot requiring G-space interpolation",&
                            &" not implemented for Gamma only!",/, &
                            &"SOLUTION: restart this calculation with", &
                            &" emtpy namelist &inputpp")')
             CALL errore ('chdens','Not implemented, please read above',1)
        ENDIF
        !
     ENDIF
#if defined(__MPI)
     ALLOCATE(aux(dfftp%nnr))
     CALL scatter_grid(dfftp, rhor, aux)
     psic(:) = cmplx(aux(:), 0.d0,kind=DP)
     DEALLOCATE(aux)
#else
     psic(:) = cmplx(rhor(:), 0.d0,kind=DP)
#endif
     CALL fwfft ('Dense', psic, dfftp)
     !
     !    we store the fourier components in the array rhog
     !
     ALLOCATE (rhog( ngm))
     rhog (:) = psic (nl (:) )
     !
  ENDIF
  !
  !     And now the plot (rhog in G-space, rhor in real space)
  !
  IF (iflag <= 1) THEN

     IF (TRIM(interpolation) == 'fourier') THEN
        CALL plot_1d (nx, m1, x0, e1, ngm, g, rhog, alat, iflag, ounit)
     ELSE
        CALL plot_1d_bspline (nx, m1, x0, e1, rhor, alat, iflag, ounit)
     ENDIF

  ELSEIF (iflag == 2) THEN

     IF (TRIM(interpolation) == 'fourier') THEN
       CALL plot_2d (nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
            at, nat, tau, atm, ityp, output_format, ounit)
     ELSE
       CALL plot_2d_bspline (nx, ny, m1, m2, x0, e1, e2, rhor, alat, &
            at, nat, tau, atm, ityp, output_format, ounit)
     ENDIF
     IF (output_format == 2.and.ionode) THEN
        WRITE (ounit, '(i4)') nat
        WRITE (ounit, '(3f8.4,i3)') ( (tau(ipol,na), ipol=1,3), 1, na=1,nat)
        WRITE (ounit, '(f10.6)') celldm (1)
        WRITE (ounit, '(3(3f12.6/))') at
     ENDIF

  ELSEIF (iflag == 3) THEN

     IF (output_format == 4.and.ionode) THEN

        ! gopenmol wants the coordinates in a separate file

        IF (fileout /= ' ') THEN
           OPEN (unit = ounit+1, file = trim(fileout)//'.xyz', &
                form = 'formatted', status = 'unknown')
           WRITE( stdout, '(5x,"Writing coordinates to file ",a)') &
                trim(fileout)//'.xyz'
        ELSE
           OPEN (unit = ounit+1, file = 'coord.xyz', &
                form = 'formatted', status = 'unknown')
           WRITE( stdout, '("Writing coordinates to file coord.xyz")')
        ENDIF
     ENDIF


     IF (output_format == 5.and.ionode) THEN
        !
        ! XCRYSDEN FORMAT
        !
        CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
        CALL xsf_fast_datagrid_3d &
             (rhor, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, ounit)

     ELSEIF (output_format == 6.and.ionode ) THEN
        !
        ! GAUSSIAN CUBE FORMAT
        !
        IF (TRIM(interpolation) == 'fourier') THEN
           CALL write_cubefile (alat, at, bg, nat, tau, atm, ityp, rhor, &
             dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, ounit)
        ELSE
           CALL plot_3d_bspline(celldm(1), at, nat, tau, atm, ityp, rhor,&
                nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, &
                ounit, rhotot)
        END IF

     ELSEIF (ionode) THEN
        !
        ! GOPENMOL OR XCRYSDEN FORMAT
        !
        IF (fast3d) THEN

           CALL plot_fast (celldm (1), at, nat, tau, atm, ityp, &
               dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, rhor, &
               bg, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, &
               rhotot)
        ELSE
           IF (nx<=0 .or. ny <=0 .or. nz <=0) &
               CALL errore("chdens","nx,ny,nz, required",1)

           IF (TRIM(interpolation) == 'fourier') THEN
              CALL plot_3d (celldm (1), at, nat, tau, atm, ityp, ngm, g, rhog,&
                   nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, &
                   ounit, rhotot)
           ELSE
              CALL plot_3d_bspline(celldm(1), at, nat, tau, atm, ityp, rhor,&
                   nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, &
                   ounit, rhotot)
           ENDIF
           !
        ENDIF
     ENDIF

  ELSEIF (iflag == 4) THEN
     radius = radius / alat
     CALL plot_2ds (nx, ny, radius, ngm, g, rhog, output_format, ounit)
  ELSE

     CALL errore ('chdens', 'wrong iflag', 1)

  ENDIF
  !
  WRITE(stdout, '(5x,"Plot Type: ",a,"   Output format: ",a)') &
       plotname(iflag), formatname(output_format)
  !
  IF (allocated(rhog)) DEALLOCATE(rhog)
  DEALLOCATE(rhor)
  DEALLOCATE(tau)
  DEALLOCATE(ityp)

END SUBROUTINE chdens
!
!-----------------------------------------------------------------------
SUBROUTINE plot_1d (nx, m1, x0, e, ngm, g, rhog, alat, iflag, ounit)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY:  pi
  USE io_global, ONLY : stdout, ionode
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE
  INTEGER :: nx, ngm, iflag, ounit
  ! number of points along the line
  ! number of G vectors
  ! type of plot
  ! output unit

  real(DP) :: e (3), x0 (3), m1, alat, g (3, ngm)
  ! vector defining the line
  ! origin of the line
  ! modulus of e
  ! lattice parameter
  ! G-vectors

  COMPLEX(DP) :: rhog (ngm)
  ! rho or polarization in G space
  INTEGER :: i, ig
  real(DP) :: rhomin, rhomax, rhoint, rhoim, xi, yi, zi, deltax, arg, gr, gg
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated charge
  ! integrated imaginary charge
  ! coordinates of a 3D point
  ! steps along the line
  ! the argument of the exponential
  ! |G|*|r|

  COMPLEX(DP) :: rho0g, carica (nx)

  deltax = m1 / (nx - 1)
  carica(:) = (0.d0,0.d0)
  IF (iflag == 1) THEN
     DO i = 1, nx
        xi = x0 (1) + (i - 1) * deltax * e (1)
        yi = x0 (2) + (i - 1) * deltax * e (2)
        zi = x0 (3) + (i - 1) * deltax * e (3)
        !
        !     for each point we compute the charge from the Fourier components
        !
        DO ig = 1, ngm
           !
           !     NB: G are in 2pi/alat units, r are in alat units
           !
           arg = 2.d0 * pi * ( xi*g(1,ig) + yi*g(2,ig) + zi*g(3,ig) )
           carica(i) = carica(i) + rhog (ig) * cmplx(cos(arg),sin(arg),kind=DP)
        ENDDO
     ENDDO
  ELSEIF (iflag == 0) THEN
     !
     !     spherically averaged charge: rho0(|r|) = int rho(r) dOmega
     !     rho0(r) = 4pi \sum_G rho(G) j_0(|G||r|)
     !
     !     G =0 term
     gg=sqrt(g(1,1)**2+g(2,1)**2+g(3,1)**2)
     IF (gg<1.d-10) THEN
        DO i = 1, nx
           carica (i) = 4.d0 * pi * rhog (1)
        ENDDO
     ENDIF
     !     G!=0 terms
     DO ig = 2, ngm
        arg = 2.d0 * pi * ( x0(1)*g(1,ig) + x0(2)*g(2,ig) + x0(3)*g(3,ig) )
        !     This displaces the origin into x0
        rho0g = rhog (ig) * cmplx(cos(arg),sin(arg),kind=DP)
        !     r =0 term
        carica (1) = carica (1) + 4.d0 * pi * rho0g
        !     r!=0 terms
        DO i = 2, nx
           gr = 2.d0 * pi * sqrt(g(1,ig)**2 + g(2,ig)**2 + g(3,ig)**2) * &
                       (i-1) * deltax
           carica (i) = carica (i) + 4.d0 * pi * rho0g * sin (gr) / gr
        ENDDO

     ENDDO
  ELSE
     CALL errore ('plot_1d', ' bad type of plot', 1)
  ENDIF
  CALL mp_sum( carica, intra_bgrp_comm )
  !
  !    Here we check the value of the resulting charge
  !
  rhomin =  1.0d10
  rhomax = -1.0d10

  rhoim = 0.d0
  DO i = 1, nx
     rhomin = min (rhomin,  dble (carica (i) ) )
     rhomax = max (rhomax,  dble (carica (i) ) )
     rhoim = rhoim + abs (aimag (carica (i) ) )
  ENDDO

  rhoim = rhoim / nx
  WRITE(stdout, '(5x,"Min, Max, imaginary charge: ",3f12.6)') &
                                          rhomin, rhomax, rhoim
  !
  !       we print the charge on output
  !
  IF (ionode) THEN
     IF (iflag == 1) THEN
        DO i = 1, nx
           WRITE (ounit, '(2f20.10)') deltax*dble(i-1), dble(carica(i))
        ENDDO
     ELSE
        rhoint = 0.d0
        DO i = 1, nx
           !
           !       simple trapezoidal rule: rhoint=int carica(i) r^2(i) dr
           !
           rhoint = rhoint + dble(carica(i)) * (i-1)**2 * (deltax*alat)**3
           WRITE (ounit, '(3f20.10)') deltax*dble(i-1), dble(carica(i)), rhoint
        ENDDO
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE plot_1d
!
!-----------------------------------------------------------------------
SUBROUTINE plot_2d (nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
     at, nat, tau, atm, ityp, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout, ionode
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  INTEGER :: nx, ny, ngm, nat, ityp (nat), output_format, ounit
  ! number of points along x
  ! number of points along y
  ! number of G vectors
  ! number of atoms
  ! types of atoms
  ! output unit
  ! output format
  CHARACTER(len=3) :: atm(*) ! atomic symbols
  real(DP) :: e1(3), e2(3), x0(3), m1, m2, g(3,ngm), alat, &
       tau(3,nat), at(3,3)
  ! vectors e1, e2 defining the plane
  ! origin
  ! modulus of e1
  ! modulus of e2
  ! G-vectors

  COMPLEX(DP) :: rhog (ngm)
  ! rho or polarization in G space
  INTEGER :: i, j, ig

  real(DP) :: rhomin, rhomax, rhoim, deltax, deltay
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
  COMPLEX(DP), ALLOCATABLE :: eigx (:), eigy (:), carica(:,:)

  ALLOCATE (eigx(  nx))
  ALLOCATE (eigy(  ny))
  ALLOCATE (carica( nx , ny))

  deltax = m1 / (nx - 1)
  deltay = m2 / (ny - 1)

  carica(:,:) = (0.d0,0.d0)
  DO ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
     ! These factors are calculated and stored in order to save CPU time
     !
     DO i = 1, nx
        eigx (i) = exp ( (0.d0, 1.d0) * 2.d0 * pi * ( (i - 1) * deltax * &
             (e1(1) * g(1,ig) + e1(2) * g(2,ig) + e1(3) * g(3,ig) ) + &
             (x0 (1) * g(1,ig) + x0 (2) * g(2,ig) + x0 (3) * g(3,ig) ) ) )
     ENDDO
     DO j = 1, ny
        eigy (j) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (j - 1) * deltay * &
             (e2(1) * g(1,ig) + e2(2) * g(2,ig) + e2(3) * g(3,ig) ) )
     ENDDO
     DO j = 1, ny
        DO i = 1, nx
           carica (i, j) = carica (i, j) + rhog (ig) * eigx (i) * eigy (j)
        ENDDO
     ENDDO
  ENDDO
  CALL mp_sum( carica, intra_bgrp_comm )
  !
  !    Here we check the value of the resulting charge
  !
  rhomin =  1.0d10
  rhomax = -1.0d10

  rhoim = 0.d0
  DO i = 1, nx
     DO j = 1, ny
        rhomin = min (rhomin,  dble (carica (i, j) ) )
        rhomax = max (rhomax,  dble (carica (i, j) ) )
        rhoim = rhoim + abs (aimag (carica (i, j) ) )
     ENDDO

  ENDDO

  rhoim = rhoim / nx / ny
  WRITE(stdout, '(5x,"Min, Max, imaginary charge: ",3f12.6)') &
                 rhomin, rhomax, rhoim

  !
  !     and we print the charge on output
  !
  IF (ionode) THEN
     IF (output_format == 0) THEN
        !
        !     gnuplot format
        !
        !         write(ounit,'(2i6)') nx,ny
        DO i = 1, nx
           WRITE (ounit, '(e25.14)') (  dble(carica(i,j)), j = 1, ny )
           WRITE (ounit, * )
        ENDDO
     ELSEIF (output_format == 1) THEN
        !
        !     contour.x format
        !
        WRITE (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
        WRITE (ounit, '(4e25.14)') ( (  dble(carica(i,j)), j = 1, ny ), i = 1, nx )
     ELSEIF (output_format == 2) THEN
        !
        !     plotrho format
        !
        WRITE (ounit, '(2i4)') nx - 1, ny - 1
        WRITE (ounit, '(8f8.4)') (deltax * (i - 1) , i = 1, nx)
        WRITE (ounit, '(8f8.4)') (deltay * (j - 1) , j = 1, ny)
        WRITE (ounit, '(6e12.4)') ( (  dble(carica(i,j)), i = 1, nx ), j = 1, ny )
        WRITE (ounit, '(3f8.4)') x0
        WRITE (ounit, '(3f8.4)') (m1 * e1 (i) , i = 1, 3)
        WRITE (ounit, '(3f8.4)') (m2 * e2 (i) , i = 1, 3)

     ELSEIF (output_format == 3) THEN
        !
        ! XCRYSDEN's XSF format
        !
        CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
        CALL xsf_datagrid_2d (carica, nx, ny, m1, m2, x0, e1, e2, alat, ounit)
     ELSEIF (output_format == 7) THEN
        !
        !     gnuplot format : x, y, f(x,y)
        !
        DO i=1, nx
           DO j=1, ny 
              WRITE (ounit, '(3e20.8)')  alat*deltax * (i - 1), &
                      alat*deltay * (j - 1), dble(carica(i,j))
           ENDDO
           WRITE(ounit, *)
        ENDDO
     ELSE
        CALL errore('plot_2d', 'wrong output_format', 1)
     ENDIF
  ENDIF

  DEALLOCATE (carica)
  DEALLOCATE (eigy)
  DEALLOCATE (eigx)
  RETURN
END SUBROUTINE plot_2d
!
!-----------------------------------------------------------------------
SUBROUTINE plot_2ds (nx, ny, x0, ngm, g, rhog, output_format, ounit)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY:  pi
  USE io_global, ONLY : stdout, ionode
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  INTEGER :: nx, ny, ngm, ounit, output_format
  ! number of points along x
  ! number of points along y
  ! number of G vectors
  ! output unit

  real(DP) :: x0, g (3, ngm)
  ! radius of the sphere
  ! G-vectors

  COMPLEX(DP) :: rhog (ngm)
  ! rho or polarization in G space
  INTEGER :: i, j, ig

  real(DP), ALLOCATABLE :: r (:,:,:)
  real(DP) :: theta, phi, rhomin, rhomax, rhoim, deltax, deltay
  ! the point in space
  ! the position on the sphere
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
  COMPLEX(DP), ALLOCATABLE :: carica (:,:)
  COMPLEX(DP) :: eig

  ALLOCATE (carica( nx , ny))
  ALLOCATE (r (3, nx , ny))

  deltax = 2.d0 * pi / (nx - 1)

  deltay = pi / (ny - 1)

  carica(:,:) = (0.d0,0.d0)
  DO j = 1, ny
     DO i = 1, nx
        phi = (i - 1) * deltax
        theta = (j - 1) * deltay
        r (1, i, j) = x0 * sin (theta) * cos (phi)
        r (2, i, j) = x0 * sin (theta) * sin (phi)
        r (3, i, j) = x0 * cos (theta)
     ENDDO
  ENDDO
  DO ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
     ! These factors are calculated and stored in order to save CPU time
     !
     DO j = 1, ny
        DO i = 1, nx
           eig = exp ( (0.d0,1.d0) * 2.d0 * pi * &
               ( r(1,i,j)*g(1,ig) + r(2,i,j)*g(2,ig) + r(3,i,j)*g(3,ig) ) )
           carica (i, j) = carica (i, j) + rhog (ig) * eig
        ENDDO
     ENDDO
  ENDDO
  CALL mp_sum( carica, intra_bgrp_comm )
  !
  !    Here we check the value of the resulting charge
  !
  rhomin =  1.0d10
  rhomax = -1.0d10

  rhoim = 0.d0
  DO i = 1, nx
     DO j = 1, ny
        rhomin = min (rhomin,  dble (carica (i, j) ) )
        rhomax = max (rhomax,  dble (carica (i, j) ) )
        rhoim = rhoim + abs (aimag (carica (i, j) ) )
     ENDDO

  ENDDO

  rhoim = rhoim / nx / ny
  WRITE(stdout, '(5x,"Min, Max, imaginary charge: ",3f12.6)') &
                  rhomin, rhomax, rhoim
  !
  !     and we print the charge on output
  !
  IF (ionode) THEN
     IF (output_format==0) THEN
        !
        !     gnuplot format
        !
        WRITE (ounit, '(2i8)') nx, ny
        DO i = 1, nx
           WRITE (ounit, '(e25.14)') (  dble(carica(i,j)), j = 1, ny )
        ENDDO
     ELSEIF (output_format==1) THEN
        !
        !     contour.x format
        !
        WRITE (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
        WRITE (ounit, '(4e25.14)') ( (  dble(carica(i,j)), j = 1, ny ), i = 1, nx )
     ELSE
        CALL errore ('plot_2ds', 'not implemented plot', 1)

     ENDIF
  ENDIF
  DEALLOCATE (carica)
  DEALLOCATE (r)
  RETURN

END SUBROUTINE plot_2ds
!
!-----------------------------------------------------------------------
SUBROUTINE plot_3d (alat, at, nat, tau, atm, ityp, ngm, g, rhog, &
     nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, &
     rhotot)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY:  pi
  USE io_global, ONLY : stdout, ionode
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  INTEGER :: nat, ityp (nat), ngm, nx, ny, nz, output_format, ounit
  ! number of atoms
  ! type of atoms
  ! number of G vectors
  ! number of points along x, y, z
  ! output format
  ! output unit
  CHARACTER(len=3) :: atm(*)

  real(DP) :: alat, tau(3,nat), at(3,3), g(3,ngm), x0(3), &
                   e1(3), e2(3), e3(3), m1, m2, m3
  ! lattice parameter
  ! atomic positions
  ! lattice vectors
  ! G-vectors
  ! origin
  ! vectors e1,e2,e3 defining the parallelepiped
  ! moduli of e1,e2,e3

  COMPLEX(DP) :: rhog (ngm)
  ! rho or polarization in G space
  INTEGER :: i, j, k, ig

  real(DP) :: rhomin, rhomax, rhotot, rhoabs, deltax, deltay, deltaz
  ! min, max value of the charge, total charge, total absolute charge
  ! steps along e1, e2, e3
  COMPLEX(DP), ALLOCATABLE :: eigx (:), eigy (:), eigz (:)
  real(DP), ALLOCATABLE :: carica (:,:,:)
  real(DP) :: omega

  ALLOCATE (eigx(  nx))
  ALLOCATE (eigy(  ny))
  ALLOCATE (eigz(  nz))
  ALLOCATE (carica( nx , ny , nz))

  deltax = m1 / nx
  deltay = m2 / ny
  deltaz = m3 / nz

  carica = 0.d0
  DO ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=exp(iG*e2), eigz=exp(iG*e3)
     ! These factors are calculated and stored in order to save CPU time
     !
     DO i = 1, nx
        eigx (i) = exp( (0.d0,1.d0) * 2.d0 * pi * ( (i-1) * deltax * &
             (e1(1)*g(1,ig)+e1(2)*g(2,ig)+e1(3)*g(3,ig)) + &
             ( x0(1)*g(1,ig)+ x0(2)*g(2,ig)+ x0(3)*g(3,ig)) ) )
     ENDDO
     DO j = 1, ny
        eigy (j) = exp( (0.d0,1.d0) * 2.d0 * pi * (j-1) * deltay * &
             (e2(1)*g(1,ig)+e2(2)*g(2,ig)+e2(3)*g(3,ig)) )
     ENDDO
     DO k = 1, nz
        eigz (k) = exp( (0.d0,1.d0) * 2.d0 * pi * (k-1) * deltaz * &
             (e3(1)*g(1,ig)+e3(2)*g(2,ig)+e3(3)*g(3,ig)) )
     ENDDO
     DO k = 1, nz
        DO j = 1, ny
           DO i = 1, nx
              carica (i, j, k) = carica (i, j, k) + &
                    dble (rhog (ig) * eigz (k) * eigy (j) * eigx (i) )
           ENDDO
        ENDDO
     ENDDO

  ENDDO
  CALL mp_sum( carica, intra_bgrp_comm )
  !
  !    Here we check the value of the resulting charge
  !

  CALL volume(alat,e1(1),e2(1),e3(1),omega)

  rhomin = max ( minval (carica), 1.d-10 )
  rhomax = maxval (carica)
  rhotot = sum (carica(:,:,:)) * omega * deltax * deltay * deltaz
  rhoabs = sum (abs(carica(:,:,:))) * omega * deltax * deltay * deltaz

  WRITE(stdout, '(/5x,"Min, Max, Total, Abs charge: ",2f10.6,2x, 2f10.4)')&
     rhomin, rhomax, rhotot, rhoabs

  IF (ionode) THEN
     IF (output_format == 4) THEN
        !
        ! "gOpenMol" file
        !

        CALL write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
             m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
     ELSE
        ! user has calculated for very long, be nice and write some output even
        ! if the output_format is wrong; use XSF format as default

        !
        ! XCRYSDEN's XSF format
        !
        CALL xsf_struct      (alat, at, nat, tau, atm, ityp, ounit)
        CALL xsf_datagrid_3d &
             (carica, nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, alat, ounit)
     ENDIF
  ENDIF

  DEALLOCATE (carica)
  DEALLOCATE (eigz)
  DEALLOCATE (eigy)
  DEALLOCATE (eigx)
  RETURN
END SUBROUTINE plot_3d
!
!-----------------------------------------------------------------------
SUBROUTINE plot_fast (alat, at, nat, tau, atm, ityp,&
     nr1x, nr2x, nr3x, nr1, nr2, nr3, rho, bg, m1, m2, m3, &
     x0, e1, e2, e3, output_format, ounit, rhotot)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nat, ityp(nat), nr1x, nr2x, nr3x, nr1, nr2, nr3, &
       output_format, ounit
  CHARACTER(len=3) :: atm(*)

  real(DP) :: alat, tau (3, nat), at (3, 3), rho(nr1x,nr2x,nr3x), &
       bg (3, 3), e1(3), e2(3), e3(3), x0 (3), m1, m2, m3

  INTEGER :: nx, ny, nz, nx0, ny0, nz0, nx1, ny1, nz1, nix, niy, niz, i, j, k, i1, j1, k1
  real(DP) :: rhomin, rhomax, rhotot, rhoabs
  real(DP), ALLOCATABLE :: carica (:,:,:)
  real(DP) :: deltax, deltay, deltaz
  real(DP) :: omega

  ! find FFT grid point closer to X0 (origin of the parallelepiped)
  ! (add 1 because r=0 correspond to n=1)

  nx0 = nint ( (x0(1)*bg(1,1) + x0(2)*bg(2,1) + x0(3)*bg(3,1) )*nr1) + 1
  ny0 = nint ( (x0(1)*bg(1,2) + x0(2)*bg(2,2) + x0(3)*bg(3,2) )*nr2) + 1
  nz0 = nint ( (x0(1)*bg(1,3) + x0(2)*bg(2,3) + x0(3)*bg(3,3) )*nr3) + 1
  !
  IF ( e1(2) /= 0.d0  .or.  e1(3) /= 0.d0 .or. &
       e2(1) /= 0.d0  .or.  e2(3) /= 0.d0 .or. &
       e3(1) /= 0.d0  .or.  e3(2) /= 0.d0 )   &
       CALL errore ('plot_fast','need vectors along x,y,z',1)

  ! find FFT grid points closer to X0 + e1, X0 + e2, X0 + e3
  ! (the opposite vertex of the parallelepiped)

  nx1 = nint ( ((x0(1)+m1)*bg(1,1)+x0(2)*bg(2,1)+x0(3)*bg(3,1) )*nr1)
  ny1 = nint ( (x0(1)*bg(1,2)+(x0(2)+m2)*bg(2,2)+x0(3)*bg(3,2) )*nr2)
  nz1 = nint ( (x0(1)*bg(1,3)+x0(2)*bg(2,3)+(x0(3)+m3)*bg(3,3) )*nr3)

  ! find number of intervals between points
  nix = nx1 - nx0 + 1
  niy = ny1 - ny0 + 1
  niz = nz1 - nz0 + 1

  IF ( output_format == 3 ) THEN
     ! XSF grids require one more point at the end of the parallelepiped sides
     nx1 = nx1 + 1
     ny1 = ny1 + 1
     nz1 = nz1 + 1
     nx = nix + 1
     ny = niy + 1
     nz = niz + 1
  ELSE
     nx = nix
     ny = niy
     nz = niz
  END IF

  ALLOCATE ( carica(nx, ny, nz) )

  carica = 0.d0
  DO k = nz0, nz1
     k1 = mod(k, nr3)
     IF (k1<=0) k1 = k1 + nr3
     DO j = ny0, ny1
        j1 = mod(j, nr2)
        IF (j1<=0) j1 = j1 + nr2
        DO i = nx0, nx1
           i1 = mod(i, nr1)
           IF (i1<=0) i1 = i1 + nr1
           carica (i-nx0+1, j-ny0+1, k-nz0+1) = rho(i1, j1, k1)
        ENDDO
     ENDDO
  ENDDO
  !
  ! recalculate m1, m2, m3 (the sides of the parallelepiped divided by alat)
  ! consistent with the FFT grid
  !
  WRITE( stdout,'(5x,"Requested parallelepiped sides : ",3f8.4)') m1, m2,m3
  m1 = nix * sqrt (at(1, 1) **2 + at(2, 1) **2 + at(3, 1) **2) / nr1
  m2 = niy * sqrt (at(1, 2) **2 + at(2, 2) **2 + at(3, 2) **2) / nr2
  m3 = niz * sqrt (at(1, 3) **2 + at(2, 3) **2 + at(3, 3) **2) / nr3
  WRITE( stdout,'(5x,"Redefined parallelepiped sides : ",3f8.4)') m1, m2,m3
  !
  ! recalculate x0 (the origin of the parallelepiped)
  ! consistent with the FFT grid
  !
  WRITE( stdout,'(5x,"Requested parallelepiped origin: ",3f8.4)') x0
  x0(1)=(nx0-1)*at(1,1)/ nr1 +(ny0-1)*at(1,2)/ nr2 +(nz0-1)*at(1,3)/ nr3
  x0(2)=(nx0-1)*at(2,1)/ nr1 +(ny0-1)*at(2,2)/ nr2 +(nz0-1)*at(2,3)/ nr3
  x0(3)=(nx0-1)*at(3,1)/ nr1 +(ny0-1)*at(3,2)/ nr2 +(nz0-1)*at(3,3)/ nr3
  WRITE( stdout,'(5x,"Redefined parallelepiped origin: ",3f8.4)') x0

  deltax = m1/nix
  deltay = m2/niy
  deltaz = m3/niz
  !
  !    Here we check the value of the resulting charge
  !
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)

  rhomin = max ( minval (carica), 1.d-10 )
  rhomax = maxval (carica)
  rhotot = sum (carica(1:nix,1:niy,1:niz)) * omega * deltax * deltay * deltaz
  rhoabs = sum (abs(carica(1:nix,1:niy,1:niz))) * omega * deltax * deltay * deltaz

  WRITE(stdout, '(/5x,"Min, Max, Total, Abs charge: ",4f10.6)') rhomin, &
       rhomax, rhotot, rhoabs

  IF (output_format == 4) THEN
     !
     !     "gopenmol" file
     !
     CALL write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
          m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
  ELSE
     !
     ! write XSF format
     !
     CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
     CALL xsf_datagrid_3d (carica, nx, ny, nz, m1, m2, m3, x0, &
          e1, e2, e3, alat, ounit)
  ENDIF
  !
  DEALLOCATE (carica)
  RETURN

END SUBROUTINE plot_fast
!
!-----------------------------------------------------------------------
SUBROUTINE write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
     m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
  !-----------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : bohr => BOHR_RADIUS_ANGS, eps4
  IMPLICIT NONE
  INTEGER :: nat, ityp (nat), nx, ny, nz, ounit
  real(DP) :: alat, tau (3, nat), at (3, 3), rhomax, x0 (3), &
       m1, m2, m3, carica (nx, ny, nz)
  CHARACTER(len=3) :: atm(*)
  !
  INTEGER, PARAMETER :: MAXATOMS = 999
  INTEGER :: natoms
  CHARACTER(len=2) TYPE (MAXATOMS)
  INTEGER :: n1, n2, n3, na, i
  real(DP) :: atoms (3, MAXATOMS), r (3), x, y, z
  real(DP) :: sidex, sidey, sidez
  !
  !   sides of the parallelepiped in A
  !
  sidex = m1 * alat * bohr
  sidey = m2 * alat * bohr
  sidez = m3 * alat * bohr

  ! really bad algorithm to  generate (hopefully) all atoms
  ! that are inside the visualization box

  natoms = 0
  DO n1 = - 3, + 3
     DO n2 = - 3, + 3
        DO n3 = - 3, + 3
           DO i = 1, 3
              r (i) = n1 * at (i, 1) + n2 * at (i, 2) + n3 * at (i, 3)
           ENDDO
           DO na = 1, nat
              ! x,y,z are in A
              x = (tau (1, na) + r (1) - x0 (1) ) * alat * bohr
              y = (tau (2, na) + r (2) - x0 (2) ) * alat * bohr
              z = (tau (3, na) + r (3) - x0 (3) ) * alat * bohr
              IF ( x > -eps4 .and. x < sidex+eps4 .and. &
                   y > -eps4 .and. y < sidey+eps4 .and. &
                   z > -eps4 .and. z < sidez+eps4 ) THEN
                 natoms = natoms + 1
                 IF (natoms>MAXATOMS) THEN
                    WRITE(stdout, '(" MAXATOMS (",i4,") Exceeded, " &
                         &       ,"Truncating " )') MAXATOMS
                    natoms = MAXATOMS
                    GOTO 10
                 ENDIF
                 !
                 atoms (1, natoms) = x
                 atoms (2, natoms) = y
                 atoms (3, natoms) = z
                 !
                 TYPE(natoms)=atm(ityp(na))
              ENDIF
           ENDDO
        ENDDO
     ENDDO

  ENDDO

10 WRITE( stdout,'(5x,"Found ",i4," atoms in the box")') natoms
  WRITE(ounit,'("  3 2")')
  WRITE(ounit,'(3i5)') nz,ny,nx
  WRITE(ounit,'(6f10.4)') 0.0d0,sidez,0.0d0,sidey,0.0d0,sidex
  DO n3=1,nz
     DO n2 = 1, ny
        DO n1 = 1, nx
           WRITE (ounit, '(f20.10)') carica (n1, n2, n3)
        ENDDO
     ENDDO
  ENDDO
  !
  ! gopenmol needs atomic positions in a separate file
  !
  WRITE(ounit+1,'(i4,/)') natoms
  WRITE(ounit+1,'(2x,a2,3f9.4)') (TYPE(na),( atoms(i,na), i=1,3 ), na=1,natoms )
  !
  RETURN
END SUBROUTINE write_openmol_file
!
SUBROUTINE isostm_plot(rhor, nr1x, nr2x, nr3x, & 
                      isovalue, heightmin, heightmax, direction)
  !-----------------------------------------------------------------------
  !
  !   Written by Andrea Cepellotti (2011), modified by Marco Pividori (2014) 
  !   to better interface with the postprocessing suite of QE
  !
  !      This subroutine calculates 2D images of STM as isosurface of 
  !      integrated ldos.
  !      It receives as input the STM charge density (that will be
  !      overwritten!) and the dimension of the grid in the real space.
  !      Works only for surfaces perpendicular to idir=3, searching for the
  !      highest isovalue found from heightmax to heightmin or viceversa
  !      according to the variable direction.
  !      
  !
  !      DESCRIPTION of the INPUT CARD  ISOSTM :
  !
  !      isovalue     ! (real) value of the charge of the isosurface
  !                   ! default value -> 0.0d0
  !      heightmin    ! (real) minimum value of the plane in which searching for the isosurface
  !                   ! default value -> 0.0d0
  !      heightmax    ! (real) maximum value of the plane in which searching for the isosurface
  !                   ! default value -> 1.0d0
  !                   ! the two parameters above are in percentage with respect to the 
  !                   ! height of the cell, i.e. between 0.0 and 1.0.
  !                   ! If heightmax < heightmin, it treats it as if it's in the 
  !                   !  upper periodically repeated slab.
  !                   ! Put heightmin somewhere in the bulk and heightmax in the vacuum
  !      direction    ! (integer) direction along z of the scan for the stm image:
  !                   ! if direction = 1 generates the isosurface as seen from heightmax to heightmin
  !                   ! if direction =-1 generates the isosurface as seen from heightmin to heightmax
  !                   ! default value -> 1

  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp,         ONLY : mp_bcast

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nr1x, nr2x, nr3x
  !dimension of the grid in the REAL SPACE

  INTEGER :: ios

  real(DP) :: rhor(nr1x*nr2x*nr3x)
  ! charge in R space

  REAL(DP), ALLOCATABLE :: image (:,:)  
  ! array for storing z coordinates

  REAL(DP), ALLOCATABLE :: reorder (:)
  ! temporary array used to reorder z coord if  heightmax < heightmin
  INTEGER :: kmin,kmax,deltakz, ir, ir2,direction
  ! min fft z value
  ! max fft z value
  ! difference between kmin and kmax
  ! counters on grid
  ! direction of scan
   
  REAL(DP) :: isovalue,heightmin,heightmax
  ! input parameters
  REAL(DP) :: maximum,minimum
  ! max and min value of iLDOS
  LOGICAL :: saturation			! check on the image

  INTEGER :: i, j, k

  !
  ! algorithm to find the isovalue
  !

  kmin=NINT(heightmin*nr3x)
  kmax=NINT(heightmax*nr3x)   
  deltakz=0

  ! if heightmin > heightmax, translate the z coordinates so that heightmin < heightmax

  IF ( heightmin > heightmax ) THEN
    ALLOCATE (reorder(nr1x*nr2x*nr3x))
    kmin=NINT(heightmin*nr3x)
    kmax=NINT(heightmax*nr3x)       
    deltakz=nr3x-kmin+1
   
    DO k = 1,nr3x
      DO j = 1,nr2x
        DO i = 1,nr1x
            ir  = i + (j - 1) * nr1x + (k - 1) * nr1x * nr2x
            ir2 = i + (j - 1) * nr1x + ( mod((k + deltakz),nr3x) &
                     - 1) * nr1x * nr2x
            reorder(ir2) = rhor(ir)               
        ENDDO    
      ENDDO
    ENDDO
    rhor=reorder
    DEALLOCATE (reorder)
   
    kmin= mod( kmin + deltakz, nr3x)
    kmax= mod( kmax + deltakz, nr3x)
 
  ENDIF

  IF (kmax > nr3x .or. kmin > nr3x .or. kmax <0 .or. kmin <0) THEN
    CALL errore('isostm','problem with heightmax/min',1)
  ENDIF

  !
  ! now search for the isosurface
  !

  ! if heightmin is 0.0d0, the lower limit is set to alat/nr3x
  IF (kmin == 0) kmin = 1  

  ALLOCATE (image(nr1x,nr2x))
  image=0.d0
  minimum=10.d0
  maximum=0.d0
  saturation=.false.

  DO k = kmin, kmax, direction
     DO j = 1, nr2x
        DO i = 1, nr1x
           ir = i + (j - 1) * nr1x + (k - 1) * nr1x * nr2x
           IF ( dble (rhor (ir) ) >= isovalue ) THEN
              image (i,j) = k
              IF (k==kmax) THEN
                 saturation=.true.
              ENDIF
              IF (k < NINT(heightmin*nr3x) ) THEN
                 image (i,j) = image (i,j) + mod((k + deltakz),nr3x)
              ENDIF
           ENDIF
           IF (dble (rhor (ir) ) < minimum ) THEN
              minimum = dble (rhor(ir))
           ELSE IF (dble (rhor (ir) ) > maximum ) THEN
              maximum = dble (rhor (ir))
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  WRITE( stdout, * )
  WRITE( stdout, * )  '    image of z coordinates of the constant isovalue'
  WRITE( stdout, * )  '    -----------------------------------------------'  
  WRITE( stdout, * )  '    max density found: ',maximum       
  WRITE( stdout, * )  '    min density found: ',minimum
  WRITE( stdout, * )  '    isovalue:          ', isovalue
 
  IF (minimum > isovalue) CALL errore('isostm','too low isovalue',1)
  IF (maximum < isovalue) CALL errore('isostm','too high isovalue',1)     

  IF (saturation) THEN
      WRITE( stdout, * )  '!! WARNING: possibly saturated image, change heights or isovalue'
  ENDIF

  !--------  
  !WARNING!  We overwrite image(x,y) in the 3D real grid to use the FFT3D algorithm
  !--------

  !overwriting charge with image(x,y) (z is a dummy variable)
  DO k = 1, nr3x
     DO j = 1, nr2x
        DO i = 1, nr1x
           ir = i + (j - 1) * nr1x + (k - 1) * nr1x * nr2x
           rhor(ir) = image(i,j)
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(image)

END SUBROUTINE isostm_plot
