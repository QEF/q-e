
! Copyright (C) 2001-2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM average
  !-----------------------------------------------------------------------
  !
  !      Compute planar and macroscopic averages of a quantity (e.g. charge)
  !      in real space on a 3D FFT mesh. The quantity is read from a file
  !      produced by "pp.x", or from multiple files as follows:
  !          Q(i,j,k) = \sum_n w_n q_n(i,j,k)
  !      where q_n is the quantity for file n, w_n is a user-supplied weight
  !      The planar average is defined as
  !         p(k) = \sum_{i=1}^{N_1} \sum_{j=1}^{N_2} Q(i,j,k) / (N_1 N_2)
  !      along direction 3, and the like for directions 1 and 2;
  !      N_1, N_2, N_3 are the three dimensions of the 3D FFT.
  !      Note that if Q is a charge density whose integral is Z_v:
  !         Z_v = \int p(z) dV = \sum_k p(k) \Omega/N_3
  !      where \Omega is the size of the unit cell (or supercell)
  !      The planar average is then interpolated on the specified number
  !      of points supplied in input and written to file "avg.dat"
  !      The macroscopic average is defined as
  !         m(z) = \int_z^{z+a} p(z) dz
  !      where a is the size of the window (supplied in input)
  !
  !      Input variables
  !
  !      nfile        the number of files contaning the desired quantities
  !                   All files must refer to the same physical system!
  ! for each file:
  !      filename     the name of the n-th file
  !      weight       the weight w_n of the quantity read from n-th file
  !      .
  !      .
  ! end
  !      npt          the number of points for the final interpolation of
  !                   the planar and macroscopic averages, as written to file
  !                   If npt <= N_idir (see below) no interpolation is done,
  !                   the N_idir FFT points in direction idir are printed.
  !      idir         1,2 or 3. Planar average is done in the plane orthogonal
  !                   to direction "idir", as defined for the crystal cell
  !      awin         the size of the window for macroscopic average (a.u.)
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE parameters,           ONLY : ntypx
  USE constants,            ONLY : pi
  USE run_info,             ONLY : title
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : ibrav, alat, omega, celldm, tpiba, &
                                   tpiba2, at, bg
  USE gvect,                ONLY : gcutm
  USE gvecs,                ONLY : doublegrid, gcutms, dual
  USE gvecw,                ONLY : ecutwfc
  USE fft_base,             ONLY : dfftp
  USE fft_types,            ONLY : fft_type_allocate
  USE fft_base,             ONLY : dffts
  USE ions_base,            ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE lsda_mod,             ONLY : nspin
  USE wavefunctions_module, ONLY : psic
  USE io_files,             ONLY : iunpun
  USE scf,                  ONLY : rho
  USE mp_global,            ONLY : mp_startup
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE environment,          ONLY : environment_start, environment_end
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: ibravs, nr1sxa, nr2sxa, nr3sxa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  INTEGER :: npt, inunit, plot_num, ios, nfile, ifile, nmacro,  &
       ir, i, j, k
  ! number of points
  ! number of input unit
  ! number of the plot
  ! integer unit for I/O control
  ! the number of files
  ! counter on the files
  ! points in the window
  ! counter on atoms
  ! counter on mesh points
  ! counters on directions

  REAL(DP) :: awin, deltaz
  ! length of the window
  ! the delta on the thick mesh
  REAL (dp), ALLOCATABLE :: weight (:)
  ! the weight of each file
  REAL(dp), ALLOCATABLE :: gre(:), gim(:), macros(:)
  ! the function to average in thick mesh (real part)
  ! the function to average in thick mesh (im. part)
  ! the macroscopic average
  REAL(DP), ALLOCATABLE :: funcr (:), funci (:)
  ! the function to average (real part)
  ! the function to average (im. part)

  REAL(DP) :: celldms (6), gcutmsa, duals, ecuts, zvs (ntypx), ats(3,3)
  REAL(DP) :: leng
  REAL(DP), ALLOCATABLE :: taus (:,:)
  INTEGER, ALLOCATABLE :: ityps (:)
  CHARACTER (len=3) :: atms(ntypx)

  INTEGER :: nfft, nfftx, idir

  CHARACTER (len=256), ALLOCATABLE :: filename (:)
  ! names of the files with the charge
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'AVERAGE' )
  !
  ! Works for parallel machines but only for one processor !!!
  !
  IF ( ionode ) THEN
     !
     inunit = 5
     READ (inunit, *, err = 1100, iostat = ios) nfile
     IF ( nfile <=0 ) CALL errore ('average ', 'nfile is wrong ', 1)
     ALLOCATE ( filename (nfile) )
     ALLOCATE ( weight (nfile) )
     DO ifile = 1, nfile
        READ (inunit, '(a)', err = 1100, iostat = ios) filename (ifile)
        READ (inunit, *, err = 1100, iostat = ios) weight (ifile)
     ENDDO
     READ (inunit, *, err = 1100, iostat = ios) npt

     IF ( npt < 0 ) CALL errore ('average', ' wrong npt', 1)
     READ (inunit, *, err = 1100, iostat = ios) idir
     READ (inunit, *, err = 1100, iostat = ios) awin

1100 CALL errore ('average', 'readin input', abs (ios) )

     CALL read_io_header(filename (1), title, dfftp%nr1x, dfftp%nr2x, &
                         dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
          nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num)
     nspin = 1
     CALL latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm(1)  ! define alat
     at = at / alat    ! bring at in units of alat

     CALL recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2

     IF (idir==1) THEN
        nfft=dfftp%nr1
        nfftx=dfftp%nr1x
        leng=alat*sqrt(at(1,1)**2+at(2,1)**2+at(3,1)**2)
     ELSEIF (idir==2) THEN
        nfft=dfftp%nr2
        nfftx=dfftp%nr2x
        leng=alat*sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2)
     ELSEIF (idir==3) THEN
        nfft=dfftp%nr3
        nfftx=dfftp%nr3x
        leng=alat*sqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2)
     ELSE
        CALL errore('average','idir is wrong',1)
     ENDIF
     IF ( npt < nfft ) THEN
        WRITE (stdout, '("Notice: npt was too small and is set to ",i4, &
             & " (number of points in 1D FFT)")' ) nfft
        npt = nfft
     END IF
     ALLOCATE(tau (3, nat))
     ALLOCATE(ityp(nat))
     doublegrid = dual>4.d0
     IF (doublegrid) THEN
        gcutms = 4.d0 * ecutwfc / tpiba2
     ELSE
        gcutms = gcutm
     ENDIF
     ! not sure whether this is the correct thing to do in presence
     ! of a double grid, but the info on nrXs is not read from file!
     dffts%nr1 = dfftp%nr1 ; dffts%nr2 = dfftp%nr2 ; dffts%nr3 = dfftp%nr3
     ! as above: this can be used in allocate_fft
     nks = 0

     CALL volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)

     CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm )
     CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm)
     CALL data_structure ( gamma_only )
     CALL allocate_fft ( )
     !
     rho%of_r = 0.d0
     !
     ! Read first file
     !
     CALL plot_io (filename (1), title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, &
          dfftp%nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
          plot_num, atm, ityp, zv, tau, rho%of_r, -1)
     !
     DO ir = 1, dfftp%nnr
        psic (ir) = weight (1) * cmplx(rho%of_r(ir, 1),0.d0,kind=DP)
     ENDDO
     !
     !       Now we open all the other files
     !
     iunpun = 4
     !
     ! Read following files (if any), verify consistency
     ! Note that only rho is read; all other quantities are discarded
     !
     DO ifile = 2, nfile
        ALLOCATE  (taus( 3 , nat))
        ALLOCATE  (ityps( nat))
        !
        CALL plot_io (filename (ifile), title, nr1sxa, nr2sxa, nr3sxa, &
             nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
             duals, ecuts, plot_num, atms, ityps, zvs, taus, rho%of_r, - 1)
        !
        DEALLOCATE (ityps)
        DEALLOCATE (taus)
        !
        IF (nats>nat) CALL errore ('average', 'wrong file order? ', 1)
        IF (dfftp%nr1x/=nr1sxa.or.dfftp%nr2x/=nr2sxa) &
             CALL errore ('average', 'incompatible nr1x or nr2x', 1)
        IF (dfftp%nr1/=nr1sa.or.dfftp%nr2/=nr2sa.or.dfftp%nr3/=nr3sa) &
             CALL errore ('average', 'incompatible nr1 or nr2 or nr3', 1)
        IF (ibravs/=ibrav) CALL errore ('average', 'incompatible ibrav', 1)
        IF (gcutmsa/=gcutm.or.duals/=dual.or.ecuts/=ecutwfc ) &
             CALL errore ('average', 'incompatible gcutm or dual or ecut', 1)
        DO i = 1, 6
           IF (abs( celldm (i)-celldms (i) ) > 1.0d-7 ) &
                CALL errore ('average', 'incompatible celldm', 1)
        ENDDO
        DO ir = 1, dfftp%nnr
           psic (ir) = psic (ir) + weight(ifile) * cmplx(rho%of_r(ir, 1),0.d0,kind=DP)
        ENDDO
     ENDDO
     DEALLOCATE ( filename )
     DEALLOCATE ( weight )
     !
     !   compute the direct and reciprocal lattices
     !
     ALLOCATE (funcr(nfftx))
     ALLOCATE (funci(nfftx))
     !
     !     At this point we start the calculations, first we compute the
     !     planar averages
     !
     IF (idir==1) THEN
        DO i = 1, dfftp%nr1
           funcr (i) = 0.d0
           funci (i) = 0.d0
           DO j = 1, dfftp%nr2
              DO k = 1, dfftp%nr3
                 ir = i + (j - 1) * dfftp%nr1x + (k - 1) * dfftp%nr1x * dfftp%nr2x
                 funcr (i) = funcr (i) + dble (psic(ir))
              ENDDO
           ENDDO
           funcr (i) = funcr (i) / (dble (dfftp%nr2 * dfftp%nr3))
        ENDDO
     ELSEIF (idir==2) THEN
        DO j = 1, dfftp%nr2
           funcr (j) = 0.d0
           funci (j) = 0.d0
           DO i = 1, dfftp%nr1
              DO k = 1, dfftp%nr3
                 ir = i + (j - 1) * dfftp%nr1x + (k - 1) * dfftp%nr1x * dfftp%nr2x
                 funcr (j) = funcr (j) + dble (psic (ir) )
              ENDDO
           ENDDO
           funcr (j) = funcr (j) / (dble (dfftp%nr1 * dfftp%nr3) )
        ENDDO
     ELSEIF (idir==3) THEN
        DO k = 1, dfftp%nr3
           funcr (k) = 0.d0
           funci (k) = 0.d0
           DO j = 1, dfftp%nr2
              DO i = 1, dfftp%nr1
                 ir = i + (j - 1) * dfftp%nr1x + (k - 1) * dfftp%nr1x * dfftp%nr2x
                 funcr (k) = funcr (k) + dble (psic (ir) )
              ENDDO
           ENDDO
           funcr (k) = funcr (k) / (dble (dfftp%nr1 * dfftp%nr2) )
        ENDDO
     ELSE
        CALL errore('average','wrong idir',1)
     ENDIF
     !
     ALLOCATE ( gre(npt) )
     ALLOCATE ( gim(npt) )
     IF ( npt == nfft) THEN
        gre(:) = funcr
        gim(:) = funci
     ELSE
        !
        !     add more points to compute the macroscopic average
        !
        CALL cft (funcr, funci, nfft, nfft, nfft, - 1)
        CALL dscal (nfft, 1.d0 / nfft, funcr, 1)
        CALL dscal (nfft, 1.d0 / nfft, funci, 1)
        !
        DO k = 1, npt
           IF (k<=nfft / 2) THEN
              gre (k) = funcr (k)
              gim (k) = funci (k)
           ELSEIF (k>npt - nfft / 2) THEN
              gre (k) = funcr (k - npt + nfft)
              gim (k) = funci (k - npt + nfft)
           ELSE
              gre (k) = 0.d0
              gim (k) = 0.d0
           ENDIF
        ENDDO
        IF (mod (nfft, 2) ==0) THEN
           gre (nfft / 2 + 1) = 0.5d0 * funcr (nfft / 2 + 1)
           gim (nfft / 2 + 1) = 0.5d0 * funci (nfft / 2 + 1)
           gre (npt - nfft / 2 + 1) = gre (nfft / 2 + 1)
           gim (npt - nfft / 2 + 1) = - gim (nfft / 2 + 1)
        ELSE
           gre (nfft / 2 + 1) = funcr (nfft / 2 + 1)
           gim (nfft / 2 + 1) = funci (nfft / 2 + 1)
        ENDIF
        !
        CALL cft (gre, gim, npt, npt, npt, 1)
        !
     END IF
     !
     !     compute the macroscopic average
     !
     ALLOCATE ( macros(npt) )
     nmacro = npt * (awin / leng )
     IF (nmacro<=0) CALL errore ('average ', 'nmacro is too small ', 1)
     DO i = 1, npt
        macros (i) = 0.d0
        DO j = - nmacro / 2, nmacro / 2
           k = i + j
           IF (k <= 0) k = k + npt
           IF (k >npt) k = k - npt

           IF ( (2*j==nmacro) .or. (2*j==-nmacro) ) THEN
              macros (i) = macros (i) + 0.5d0 * gre(k)
           ELSE
              macros (i) = macros (i) + gre (k)
           ENDIF
        ENDDO
        macros (i) = macros (i) / dble (nmacro)
     ENDDO
     !
     !     print the results on output
     !
     deltaz = leng / dble (npt)
     WRITE( stdout, '(5x,"Output written to file avg.dat")')
     OPEN (unit=4, file='avg.dat', form='formatted', status='unknown')
     WRITE(4, '(3f15.9)') (deltaz*(i-1), gre(i), macros(i), i = 1, npt)
     CLOSE (unit=4, status='keep')
     !
     DEALLOCATE(macros)
     DEALLOCATE(gre)
     DEALLOCATE(gim)
     DEALLOCATE(funci)
     DEALLOCATE(funcr)
     !
  ENDIF
  !
  CALL environment_end ( 'AVERAGE' )
  !
  CALL stop_pp
  !
END PROGRAM average
