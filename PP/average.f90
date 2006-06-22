
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
PROGRAM average
  !-----------------------------------------------------------------------
  !
  !      This program calculates planar and macroscopic averages
  !      of a quantity defined on a 3D-FFT mesh. 
  !      The planar average is done on FFT mesh planes. 
  !      It reads the quantity to average, or several quantities, from
  !      one or several files and adds them with the given weights.
  !      It computes the planar average of the resulting quantity
  !      averaging on planes defined by the FFT mesh points and by one
  !      direction perpendicular to the planes.
  !      The planar average can be interpolated on a  
  !      1D-mesh with an arbitrary number of points.
  !      Finally, it computes the macroscopic average. The size
  !      of the averaging window is given as input.
  !
  !      It receive as input the following variables:
  !
  !      nfile        ! the number of 3D-FFT files
  ! for each file:
  !      filename     ! the name of the 3D-FFT file
  !      weight       ! the weight of the quantity in this file
  !      .
  !      .
  ! end
  !      npt          ! the number of points of the thick mesh
  !      idir         ! 1,2 or 3. It is the fixed index which defines
  !                   ! the planes of the planar average
  !      awin         ! the size of the window for macroscopic averages.
  !
  USE kinds,                ONLY : DP
  USE parameters,           ONLY : ntypx
  USE constants,            ONLY : pi
  USE char,                 ONLY : title
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : ibrav, alat, omega, celldm, tpiba, &
                                   tpiba2, at, bg
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   gcutm, ecutwfc, dual
  USE gsmooth,              ONLY : doublegrid, gcutms
  USE ions_base,            ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE lsda_mod,             ONLY : nspin
  USE wavefunctions_module, ONLY : psic
  USE io_files,             ONLY : nd_nmbr, iunpun
  USE scf,                  ONLY : rho
  USE mp_global,            ONLY : mpime, root
  !
  IMPLICIT NONE
  !
  INTEGER :: npixmax, nfilemax
  ! maximum number of pixel
  ! maximum number of files with charge
  !
  PARAMETER (npixmax = 5000, nfilemax = 7)
  !
  INTEGER :: ibravs, nrx1sa, nrx2sa, nrx3sa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  INTEGER :: npt, inunit, plot_num, ios, nfile, ifile, nmacro, na, &
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

  REAL(DP) :: rhodum, awin, deltaz, weight (nfilemax), gre(npixmax), &
       gim(npixmax), macros(npixmax)
  ! length of the window
  ! the delta on the thick mesh
  ! the weight of each file
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

  CHARACTER (len=256) :: filename (nfilemax)
  ! names of the files with the charge
  !
  CALL start_postproc (nd_nmbr)
  !
  ! Works for parallel machines but only for one processor !!!
  !
  IF ( mpime == root ) THEN
     !
     inunit = 5
     READ (inunit, *, err = 1100, iostat = ios) nfile
     IF (nfile.LE.0.OR.nfile.GT.nfilemax) CALL errore ('average ', &
          'nfile is wrong ', 1)
     DO ifile = 1, nfile
        READ (inunit, '(a)', err = 1100, iostat = ios) filename (ifile)
        READ (inunit, *, err = 1100, iostat = ios) weight (ifile)
     ENDDO
     READ (inunit, *, err = 1100, iostat = ios) npt

     IF (npt.LT.0.OR.npt.GT.npixmax) CALL errore ('average', ' wrong npt', 1)
     READ (inunit, *, err = 1100, iostat = ios) idir
     READ (inunit, *, err = 1100, iostat = ios) awin

1100 CALL errore ('average', 'readin input', ABS (ios) )

     CALL read_io_header(filename (1), title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
          nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num)
     nspin = 1
     CALL latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm(1)  ! define alat
     at = at / alat    ! bring at in units of alat

     CALL recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2

     IF (idir.EQ.1) THEN
        nfft=nr1
        nfftx=nrx1
        leng=alat*SQRT(at(1,1)**2+at(2,1)**2+at(3,1)**2)
     ELSEIF (idir.EQ.2) THEN
        nfft=nr2
        nfftx=nrx2
        leng=alat*SQRT(at(1,2)**2+at(2,2)**2+at(3,2)**2)
     ELSEIF (idir.EQ.3) THEN
        nfft=nr3
        nfftx=nrx3
        leng=alat*SQRT(at(1,3)**2+at(2,3)**2+at(3,3)**2)
     ELSE
        CALL errore('average','idir is wrong',1)
     ENDIF
     IF (npt.LT.nfft) CALL errore ('average', 'npt smaller than nfft', 1)

     ALLOCATE(tau (3, nat))
     ALLOCATE(ityp(nat))
     doublegrid = dual.GT.4.d0
     IF (doublegrid) THEN
        gcutms = 4.d0 * ecutwfc / tpiba2
     ELSE
        gcutms = gcutm
     ENDIF


     CALL volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)

     CALL set_fft_dim

     CALL allocate_fft
     !
     rho = 0.d0
     !
     ! Read first file
     !
     CALL plot_io (filename (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
          nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
          plot_num, atm, ityp, zv, tau, rho, -1)
     !
     DO ir = 1, nrxx
        psic (ir) = weight (1) * CMPLX (rho (ir, 1),0.d0)
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
        CALL plot_io (filename (ifile), title, nrx1sa, nrx2sa, nrx3sa, &
             nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
             duals, ecuts, plot_num, atms, ityps, zvs, taus, rho, - 1)
        !
        DEALLOCATE (ityps)
        DEALLOCATE (taus)
        !
        IF (nats.GT.nat) CALL errore ('chdens', 'wrong file order? ', 1)
        IF (nrx1.NE.nrx1sa.OR.nrx2.NE.nrx2sa) &
             CALL errore ('average', 'incompatible nrx1 or nrx2', 1)
        IF (nr1.NE.nr1sa.OR.nr2.NE.nr2sa.OR.nr3.NE.nr3sa) &
             CALL errore ('average', 'incompatible nr1 or nr2 or nr3', 1)
        IF (ibravs.NE.ibrav) CALL errore ('average', 'incompatible ibrav', 1)
        IF (gcutmsa.NE.gcutm.OR.duals.NE.dual.OR.ecuts.NE.ecutwfc ) &
             CALL errore ('average', 'incompatible gcutm or dual or ecut', 1)
        DO i = 1, 6
           IF (ABS( celldm (i)-celldms (i) ) .GT. 1.0e-7 ) &
                CALL errore ('chdens', 'incompatible celldm', 1)
        ENDDO
        DO ir = 1, nrxx
           psic (ir) = psic (ir) + weight(ifile) * CMPLX(rho(ir, 1),0.d0)
        ENDDO
     ENDDO
     !
     !   compute the direct and reciprocal lattices
     !
     ALLOCATE (funcr(nfftx))    
     ALLOCATE (funci(nfftx))    
     !
     !     At this point we start the calculations, first we compute the
     !     planar averages
     !
     IF (idir.EQ.1) THEN
        DO i = 1, nr1
           funcr (i) = 0.d0
           funci (i) = 0.d0
           DO j = 1, nr2
              DO k = 1, nr3
                 ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
                 funcr (i) = funcr (i) + DBLE (psic(ir))
              ENDDO
           ENDDO
           funcr (i) = funcr (i) / (DBLE (nr2 * nr3))
        ENDDO
     ELSEIF (idir.EQ.2) THEN
        DO j = 1, nr2
           funcr (j) = 0.d0
           funci (j) = 0.d0
           DO i = 1, nr1
              DO k = 1, nr3
                 ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
                 funcr (j) = funcr (j) + DBLE (psic (ir) )
              ENDDO
           ENDDO
           funcr (j) = funcr (j) / (DBLE (nr1 * nr3) )
        ENDDO
     ELSEIF (idir.EQ.3) THEN
        DO k = 1, nr3
           funcr (k) = 0.d0
           funci (k) = 0.d0
           DO j = 1, nr2
              DO i = 1, nr1
                 ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
                 funcr (k) = funcr (k) + DBLE (psic (ir) )
              ENDDO
           ENDDO
           funcr (k) = funcr (k) / (DBLE (nr1 * nr2) )
        ENDDO
     ELSE
        CALL errore('average','wrong idir',1)
     ENDIF
     !
     !     add more points to compute the macroscopic average
     !
     CALL cft (funcr, funci, nfft, nfft, nfft, - 1)
     CALL DSCAL (nfft, 1.d0 / nfft, funcr, 1)
     CALL DSCAL (nfft, 1.d0 / nfft, funci, 1)
     DO k = 1, npt
        IF (k.LE.nfft / 2) THEN
           gre (k) = funcr (k)
           gim (k) = funci (k)
        ELSEIF (k.GT.npt - nfft / 2) THEN
           gre (k) = funcr (k - npt + nfft)
           gim (k) = funci (k - npt + nfft)
        ELSE
           gre (k) = 0.d0
           gim (k) = 0.d0
        ENDIF
     ENDDO
     IF (MOD (nfft, 2) .EQ.0) THEN
        gre (nfft / 2 + 1) = 0.5d0 * funcr (nfft / 2 + 1)
        gim (nfft / 2 + 1) = 0.5d0 * funci (nfft / 2 + 1)
        gre (npt - nfft / 2 + 1) = gre (nfft / 2 + 1)
        gim (npt - nfft / 2 + 1) = - gim (nfft / 2 + 1)
     ELSE
        gre (nfft / 2 + 1) = funcr (nfft / 2 + 1)
        gim (nfft / 2 + 1) = funci (nfft / 2 + 1)
     ENDIF


     CALL cft (gre, gim, npt, npt, npt, 1)
     !
     !     compute the macroscopic average
     !
     nmacro = npt * (awin / leng )
     IF (nmacro.LE.0) CALL errore ('average ', 'nmacro is too small ', 1)
     DO i = 1, npt
        macros (i) = 0.d0
        DO j = - nmacro / 2, nmacro / 2
           k = i + j
           IF (k.LE.0) k = k + npt
           IF (k.GT.npt) k = k - npt

           if ( (2*j==nmacro) .or. (2*j==-nmacro) ) then
              macros (i) = macros (i) + 0.5d0 * gre(k)
           else
              macros (i) = macros (i) + gre (k)
           end if
        ENDDO
        macros (i) = macros (i) / DBLE (nmacro)
     ENDDO
     !
     !     print the results on output
     !
     deltaz = leng / DBLE (npt)


     WRITE( stdout, '(3f15.9)') (deltaz * (i - 1) , gre (i) , macros (i) , &
          i = 1, npt)
     DEALLOCATE(funci)
     DEALLOCATE(funcr)
     !
  END IF
  !
  CALL stop_pp
  !
END PROGRAM average
