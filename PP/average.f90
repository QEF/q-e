
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program average
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
#include "machine.h"
  use parameters, only: DP
  use pwcom
  USE wavefunctions,  ONLY: psic
  use io_files, only: nd_nmbr
#ifdef __PARA
  use para, only: me
#endif
  implicit none
  integer :: npixmax, nfilemax
  ! maximum number of pixel
  ! maximum number of files with charge

  parameter (npixmax = 5000, nfilemax = 7)

  integer :: ibravs, nrx1sa, nrx2sa, nrx3sa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  integer :: npt, inunit, plot_num, ios, nfile, ifile, nmacro, na, &
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

  real(kind=DP) :: rhodum, awin, deltaz, weight (nfilemax), gre(npixmax), &
       gim(npixmax), macros(npixmax)
  ! length of the window
  ! the delta on the thick mesh
  ! the weight of each file
  ! the function to average in thick mesh (real part)
  ! the function to average in thick mesh (im. part)
  ! the macroscopic average
  real(kind=DP), allocatable :: funcr (:), funci (:)
  ! the function to average (real part)
  ! the function to average (im. part)

  real(kind=DP) :: celldms (6), gcutmsa, duals, ecuts, zvs (ntypx), ats(3,3)
  real(kind=DP) :: leng
  real(kind=DP), allocatable :: taus (:,:)
  integer, allocatable :: ityps (:)
  character (len=3) :: atms(ntypx)

  integer :: nfft, nfftx, idir

  character (len=80) :: filename (nfilemax)
  ! names of the files with the charge
  !
  call start_postproc (nd_nmbr)
#ifdef __PARA
  !
  ! Works for parallel machines but only for one processor !!!
  !
  if (me == 1) then
#endif
  inunit = 5
  read (inunit, *, err = 1100, iostat = ios) nfile
  if (nfile.le.0.or.nfile.gt.nfilemax) call errore ('average ', &
       'nfile is wrong ', 1)
  do ifile = 1, nfile
     read (inunit, '(a)', err = 1100, iostat = ios) filename (ifile)
     read (inunit, *, err = 1100, iostat = ios) weight (ifile)
  enddo
  read (inunit, *, err = 1100, iostat = ios) npt

  if (npt.lt.0.or.npt.gt.npixmax) call errore ('average', ' wrong npt', 1)
  read (inunit, *, err = 1100, iostat = ios) idir
  read (inunit, *, err = 1100, iostat = ios) awin

1100 call errore ('average', 'readin input', abs (ios) )

  call plot_io (filename (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
       nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
       plot_num, atm, ityp, zv, tau, rhodum, 0)

  nspin = 1
  if (ibrav.gt.0) call latgen (ibrav, celldm, at (1, 1), &
                                       at (1, 2), at (1, 3), omega )
  alat = celldm(1)
  at=at/alat
  call recips (at (1, 1), at (1, 2), at (1, 3), bg (1, 1), bg (1, 2) &
       , bg (1, 3) )

  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2

  if (idir.eq.1) then
     nfft=nr1
     nfftx=nrx1
     leng=alat*sqrt(at(1,1)**2+at(2,1)**2+at(3,1)**2)
  elseif (idir.eq.2) then
     nfft=nr2
     nfftx=nrx2
     leng=alat*sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2)
  elseif (idir.eq.3) then
     nfft=nr3
     nfftx=nrx3
     leng=alat*sqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2)
  else
     call errore('average','idir is wrong',1)
  endif
  if (npt.lt.nfft) call errore ('average', 'npt smaller than nfft', 1)

  allocate(tau (3, nat))
  allocate(ityp(nat))
  doublegrid = dual.gt.4.d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif


  call volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)

  call set_fft_dim

  call allocate_fft
  !
  rho = 0.d0
  !
  ! Read first file
  !
  call plot_io (filename (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
       nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
       plot_num, atm, ityp, zv, tau, rho, -1)
  !
  do ir = 1, nrxx
     psic (ir) = weight (1) * cmplx (rho (ir, 1),0.d0)
  enddo
  !
  !       Now we open all the other files
  !
  iunpun = 4
  !
  ! Read following files (if any), verify consistency
  ! Note that only rho is read; all other quantities are discarded
  !
  do ifile = 2, nfile
     allocate  (taus( 3 , nat))    
     allocate  (ityps( nat))    
     !
     call plot_io (filename (ifile), title, nrx1sa, nrx2sa, nrx3sa, &
          nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
          duals, ecuts, plot_num, atms, ityps, zvs, taus, rho, - 1)
     !
     deallocate (ityps)
     deallocate (taus)
     !
     if (nats.gt.nat) call errore ('chdens', 'wrong file order? ', 1)
     if (nrx1.ne.nrx1sa.or.nrx2.ne.nrx2sa) &
          call errore ('average', 'incompatible nrx1 or nrx2', 1)
     if (nr1.ne.nr1sa.or.nr2.ne.nr2sa.or.nr3.ne.nr3sa) &
          call errore ('average', 'incompatible nr1 or nr2 or nr3', 1)
     if (ibravs.ne.ibrav) call errore ('average', 'incompatible ibrav', 1)
     if (gcutmsa.ne.gcutm.or.duals.ne.dual.or.ecuts.ne.ecutwfc ) &
          call errore ('average', 'incompatible gcutm or dual or ecut', 1)
     do i = 1, 6
        if (abs( celldm (i)-celldms (i) ) .gt. 1.0e-7 ) &
             call errore ('chdens', 'incompatible celldm', 1)
     enddo
     do ir = 1, nrxx
        psic (ir) = psic (ir) + weight(ifile) * cmplx(rho(ir, 1),0.d0)
     enddo
  enddo
  !
  !   compute the direct and reciprocal lattices
  !
  allocate (funcr(nfftx))    
  allocate (funci(nfftx))    
  !
  !     At this point we start the calculations, first we compute the
  !     planar averages
  !
  if (idir.eq.1) then
     do i = 1, nr1
        funcr (i) = 0.d0
        funci (i) = 0.d0
        do j = 1, nr2
           do k = 1, nr3
              ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
              funcr (i) = funcr (i) + real (psic(ir))
           enddo
        enddo
        funcr (i) = funcr (i) / (float (nr2 * nr3))
     enddo
  elseif (idir.eq.2) then
     do j = 1, nr2
        funcr (j) = 0.d0
        funci (j) = 0.d0
        do i = 1, nr1
           do k = 1, nr3
              ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
              funcr (j) = funcr (j) + real (psic (ir) )
           enddo
        enddo
        funcr (j) = funcr (j) / (float (nr1 * nr3) )
     enddo
  elseif (idir.eq.3) then
     do k = 1, nr3
        funcr (k) = 0.d0
        funci (k) = 0.d0
        do j = 1, nr2
           do i = 1, nr1
              ir = i + (j - 1) * nrx1 + (k - 1) * nrx1 * nrx2
              funcr (k) = funcr (k) + real (psic (ir) )
           enddo
        enddo
        funcr (k) = funcr (k) / (float (nr1 * nr2) )
     enddo
  else
     call errore('average','wrong idir',1)
  endif
!
!     add more points to compute the macroscopic average
!
  call cft (funcr, funci, nfft, nfft, nfft, - 1)
  call DSCAL (nfft, 1.d0 / nfft, funcr, 1)
  call DSCAL (nfft, 1.d0 / nfft, funci, 1)
  do k = 1, npt
     if (k.le.nfft / 2) then
        gre (k) = funcr (k)
        gim (k) = funci (k)
     elseif (k.gt.npt - nfft / 2) then
        gre (k) = funcr (k - npt + nfft)
        gim (k) = funci (k - npt + nfft)
     else
        gre (k) = 0.d0
        gim (k) = 0.d0
     endif
  enddo
  if (mod (nfft, 2) .eq.0) then
     gre (nfft / 2 + 1) = 0.5d0 * funcr (nfft / 2 + 1)
     gim (nfft / 2 + 1) = 0.5d0 * funci (nfft / 2 + 1)
     gre (npt - nfft / 2 + 1) = gre (nfft / 2 + 1)
     gim (npt - nfft / 2 + 1) = - gim (nfft / 2 + 1)
  else
     gre (nfft / 2 + 1) = funcr (nfft / 2 + 1)
     gim (nfft / 2 + 1) = funci (nfft / 2 + 1)
  endif


  call cft (gre, gim, npt, npt, npt, 1)
  !
  !     compute the macroscopic average
  !
  nmacro = npt * (awin / leng )
  if (nmacro.le.0) call errore ('average ', 'nmacro is too small ', 1)
  do i = 1, npt
     macros (i) = 0.d0
     do j = - nmacro / 2, nmacro / 2
        k = i + j
        if (k.le.0) k = k + npt
        if (k.gt.npt) k = k - npt
        macros (i) = macros (i) + gre (k)
     enddo
     if (mod (nmacro, 2) .eq.0) then
        macros (i) = macros (i) / float (nmacro + 1)
     else
        macros (i) = macros (i) / float (nmacro)
     endif
  enddo
  !
  !     print the results on output
  !
  deltaz = leng / float (npt)


  write (6, '(3f15.9)') (deltaz * (i - 1) , gre (i) , macros (i) , &
       i = 1, npt)
  deallocate(funci)
  deallocate(funcr)
#ifdef __PARA
  end if
#endif
  call stop_pp
end program average
