
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine average
  !-----------------------------------------------------------------------
  !
  !      This program is used to compute the macroscopic averages of the
  !      charge density along a given direction.
  !      It read the charge density, or more than one charge from
  !      input and adds them with appropriate weights.
  !      Then it computes the planar averages of this charge
  !      along the xy planes. As a third step it interpolates
  !      the averages on a mesh with a given number of points.
  !      At the end it computes a macroscopic average with a
  !      window of fixed dimension.
  !
  !      It receive as input the following file:
  !
  !      nfile        ! the number of charge files
  !      filename     ! the name of the charge file
  !      weight       ! the weight of this charge
  !      .
  !      .
  !      npt          ! the number of points of the thick mesh
  !      awin         ! the lenght of window which computes the averages.
  !      ionflag      ! if true the ionic charge is added
  !
#include "machine.h"
  use parameters, only: DP
  use pwcom
  use io

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
  ! counter on directions

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
  real(kind=DP), allocatable :: taus (:,:)
  integer, allocatable :: ityps (:)
  character (len=3) :: atms(ntypx)

  logical :: ionflag
  ! if true the ionic charge is added

  character (len=80) :: filename (nfilemax)
  ! names of the files with the charge
  !
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
  read (inunit, *, err = 1100, iostat = ios) awin
  read (inunit, *, err = 1100, iostat = ios) ionflag

1100 call errore ('average', 'readin input', abs (ios) )

  call plot_io (filename (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
       nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
       plot_num, atm, ityp, zv, tau, rhodum, 0)
  if (npt.le.nr3) call errore ('average', 'npt smaller than nr3', 1)

  allocate(tau (3, nat) )
  allocate(ityp(nat) )
  alat = celldm (1)
  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  doublegrid = dual.gt.4.d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif

  nspin = 1
  if (ibrav.gt.0) call latgen (ibrav, celldm, at (1, 1), &
                                              at (1, 2), at (1, 3) )
  call recips (at (1, 1), at (1, 2), at (1, 3), bg (1, 1), bg (1, 2) &
       , bg (1, 3) )

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
  !       Now we open the input file and we read what has been written
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
     if (nrx1.ne.nrx1sa.or.nrx2.ne.nrx2sa) call &
          errore ('average', 'incompatible nrx1 or nrx2', 1)
     if (nr1.ne.nr1sa.or.nr2.ne.nr2sa.or.nr3.ne.nr3sa) call errore ( &
          'average', 'incompatible nr1 or nr2 or nr3', 1)
     if (ibravs.ne.ibrav) call errore ('average', 'incompatible ibrav', 1)
     if (gcutmsa.ne.gcutm.or.duals.ne.dual.or.ecuts.ne.ecutwfc ) &
          call errore ('average', 'incompatible gcutm or dual or ecut', 1)
     do i = 1, 6
        if (abs( celldm (i)-celldms (i) ) .gt. 1.0e-7 ) call errore &
             ('chdens', 'incompatible celldm', 1)
     enddo
     do ir = 1, nrxx
        psic (ir) = psic (ir) + weight (ifile) * cmplx (rho (ir, 1), &
             0.d0)
     enddo
  enddo
  !
  !   compute the direct and reciprocal lattices
  !
  allocate (funcr(nrx3))    
  allocate (funci(nrx3))    
  !
  !     At this point we start the calculations, first we compute the
  !     planar averages
  !
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
  do k = 1, nr3
     write (6, * ) k, funcr (k)
  enddo
  !
  !     add more points to compute the macroscopic average
  !
  call cft (funcr, funci, nr3, nr3, nr3, - 1)
  call DSCAL (nr3, 1.d0 / nr3, funcr, 1)
  call DSCAL (nr3, 1.d0 / nr3, funci, 1)
  do k = 1, npt
     if (k.le.nr3 / 2) then
        gre (k) = funcr (k)
        gim (k) = funci (k)
     elseif (k.gt.npt - nr3 / 2) then
        gre (k) = funcr (k - npt + nr3)
        gim (k) = funci (k - npt + nr3)
     else
        gre (k) = 0.d0
        gim (k) = 0.d0
     endif
  enddo
  if (mod (nr3, 2) .eq.0) then
     gre (nr3 / 2 + 1) = 0.5d0 * funcr (nr3 / 2 + 1)
     gim (nr3 / 2 + 1) = 0.5d0 * funci (nr3 / 2 + 1)
     gre (npt - nr3 / 2 + 1) = gre (nr3 / 2 + 1)
     gim (npt - nr3 / 2 + 1) = - gim (nr3 / 2 + 1)
  else
     gre (nr3 / 2 + 1) = funcr (nr3 / 2 + 1)
     gim (nr3 / 2 + 1) = funci (nr3 / 2 + 1)
  endif


  call cft (gre, gim, npt, npt, npt, 1)
  !
  !     compute the macroscopic average
  !
  nmacro = npt * (awin / (alat * celldm (3) ) )
  if (nmacro.le.0) call errore ('average ', 'nmacro is too small ', &
       1)
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
  deltaz = alat * celldm (3) / float (npt - 1)


  write (6, '(3f15.9)') (deltaz * (i - 1) , gre (i) , macros (i) , &
       i = 1, npt)
  deallocate(funci)
  deallocate(funcr)
  return
end subroutine average
