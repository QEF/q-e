!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine do_chdens
  !-----------------------------------------------------------------------
  !      Charge density/polarization plotting program
  !-----------------------------------------------------------------------
  !      The input data of this program are read from standard input
  !      or from a file and have the following format.
  !
  !      DESCRIPTION of the INPUT:
  !
  !-&input           Namelist &input;
  !                  BELOW IS THE DESCRIPTION OF THE VARIABLES OF THE &INPUT NAMELIST
  !
  !      nfile       the number of data files
  !
  !----FOR i = 1, nfile:
  !
  !      filepp(i)   file containing the 3D charge (produced by pp.x)
  !      weight(i)   weight - The quantity to be plotted will be
  !                  weight(1)*rho(1) + weight(2)*rho(2) + weight(3)*rho(3) + ...
  !
  !                  BEWARE: atomic coordinates are read from the first file;
  !                  if their number is different for different files,
  !                  the first file must have the largest number of atoms
  !
  !----END_FOR
  !	
  !	
  !      iflag     1 if a 1D plot is required
  !                2 if a 2D plot is required
  !                3 if a 3D plot is required
  !               30 if a "fast" 3D plot is required: the 3D plot points
  !                    are not recalculated from the Fourier components but
  !                    directly extracted from the FFT grid of the system
  !               31 if a "fast" 3D plot is required: the same as iflag=30,
  !                    but the whole unit-cell is taken. Does not require the
  !                    input of the e1,e2,e3,x0 and nx,ny,nz
  !                    NOTE: works only for XCRYSDEN format (output_format=3)
  !                4 if a 2D polar plot on a sphere is required
  !
  !      plot_out  0   plot the spherical average of the charge density
  !                1   plot the charge density
  !                2   plot the induced polarization along x
  !                3   plot the induced polarization along y
  !                4   plot the induced polarization along z
  !
  !      output_format  (ignored on 1D plot)
  !                0  format suitable for gnuplot
  !                1  format suitable for contour.x
  !                2  format suitable for plotrho
  !                3  format suitable for XCRYSDEN
  !                4  format suitable for gOpenMol
  !
  !      fileout   name of the file to which the plot is written
  !
  !----IF plot_out=2,3,4:
  !
  !      epsilon   the dielectric constant for polarization computation
  !
  !      filepol   name of an output file to which the induced polarization
  !                is written (in postproc format) for further processing
  !                (macroscopic average)
  !
  !----END_IF
  !
  !-/              END of namelist &input
  !
  !--------------- DESCRIPTION OF FURTHER STANDARD INPUT
  !
  !-IF iflag = 1      the following cards are
  !
  !         e1        3D vector which determines the plotting line
  !         x0        3D vector, origin of the line
  !         nx        number of points in the line:
  !                   rho(i) = rho( x0 + e1 * (i-1)/(nx-1) ), i=1, nx
  !
  !-ELSEIF iflag = 2  the following cards are
  !
  !         e1        3D vectors which determine the plotting plane
  !         e2           NB: the two axes must be orthogonal
  !         x0        3D vector, origin of the plane
  !         nx, ny    number of points in the plane:
  !                   rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)
  !                                      + e2 * (j-1)/(ny-1) ), i = 1, nx ; j = 1, ny
  !
  !-ELSEIF iflag = 3 or iflag = 30 the following cards are
  !         e1
  !         e2        3D vectors which determine the plotting parallelepiped
  !         e3        REQUIREMENT: the three axes must be orthogonal !!!
  !         x0        3D vector, origin of the parallelepiped
  !         nx,ny,nz  number of points in the parallelepiped:
  !                   rho(i,j,k) = rho( x0 + e1 * (i-1)/(nx-1)
  !                                        + e2 * (j-1)/(ny-1)
  !                                        + e3 * (k-1)/(nz-1) ),
  !                                i = 1, nx ; j = 1, ny ; k = 1, nz
  !
  !         IMPORTANT: all 3D vectors in a0 (alat) units
  !         BEWARE: iflag = 30 (fast plot)
  !            - works only if crystal axes are orthogonal
  !            - works only if e1 is along x, e2 along y, e3 along z
  !            - the plotting grid is fixed and determined by the FFT grid
  !            - nx, ny, nz  are ignored
  !            - the parallelepiped defined by e1, e2, e3, x0 is displaced
  !              and stretched so as to fit the FFT grid
  !         BEWARE: iflag = 30 (fast plot, whole cell)
  !            - works for all crystal axes
  !            - works only for XCRYSDEN XSF output format (i.e. output_format = 3)
  !
  !-ELSEIF iflag = 4  the following cards are
  !
  !         radius    Radius of the sphere (alat units), centered at (0,0,0)
  !         nx, ny    number of points in the polar plane:
  !                   phi(i)  = 2 pi * (i - 1)/(nx-1), i=1, nx
  !                   theta(j)=   pi * (j - 1)/(ny-1), j=1, ny
  !
  !-ENDIF
  !
  !EOF

#include "machine.h"
  use pwcom
  use io

  implicit none
  integer, parameter :: nfilemax = 7
  ! maximum number of files with charge

  integer :: inunit, ounit, iflag, ios, ipol, nfile, ifile, nx, ny, nz, &
       na, ir, i, j, ig, plot_out, output_format, plot_num

  real(kind=DP) :: e (3, 3), x0 (3), radius, m1, m2, m3, &
       weight (nfilemax), epsilon

  character (len=80) :: fileout, filepol, filename (nfilemax)

  real(kind=DP) :: celldms (6), gcutmsa, duals, ecuts, zvs(ntypx), ats(3,3)
  real(kind=DP), allocatable :: taus (:,:)
  integer :: ibravs, nrx1sa, nrx2sa, nrx3sa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  integer, allocatable :: ityps (:)
  character (len=3) :: atms(ntypx)
  character (len=80) :: filepp(nfilemax)
  real(kind=DP) :: rhodum
  complex(kind=DP), allocatable:: vgc (:)
  ! rho or polarization in G space
  logical :: fast3d

  namelist /input/  &
       nfile, filepp, weight, iflag, &
       plot_out, output_format, fileout, epsilon, filepol

  !
  !   set the DEFAULT values
  !
  nfile         = 1
  filepp(1)   = ' '
  weight(1)     = 1.0d0
  iflag         = 1
  radius        = 1.0d0
  plot_out      = 1
  output_format = 0
  fileout       = ' '
  epsilon       = 1.0d0
  filepol       = ' '

  !
  !    read and check input data
  !
  inunit = 5

  !
  ! reading the namelist input
  !
  read (5, input, err = 200, iostat = ios)
200 call error ('chdens', 'reading input namelist', abs (ios) )

  ! check for number of files
  if (nfile.le.0.or.nfile.gt.nfilemax) &
       call error ('chdens ', 'nfile is wrong ', 1)

  ! check for iflag
  if (iflag.eq.30) then
     iflag = 3
     fast3d = .true.
  end if
  if ((iflag.lt.1.or.iflag.gt.4) .and. iflag.ne.31) &
       call error ('chdens', 'iflag not implemented', 1)


  ! reading the rest of input (spanning vectors, origin, number-of points)
  if (iflag.lt.4) then
     read (inunit, *, err = 1100, iostat = ios) (e (ipol, 1), &
          ipol = 1, 3)
     if (e (1, 1) **2 + e (2, 1) **2 + e (3, 1) **2.lt.1d-3) call &
          error ('chdens', 'zero vector', 1)
  endif
  if (iflag.eq.1) then
     !
     !    reading for the 1D plot
     !
     read (inunit, *, err = 1100, iostat = ios) (x0 (ipol), ipol = 1, 3)
     read (inunit, *, err = 1100, iostat = ios) nx
  endif
  !
  if (iflag.ge.2.and.iflag.lt.4) then
     !
     !    reading for the 2D and 3D plots
     !
     read (inunit, *, err = 1100, iostat = ios) (e (ipol, 2), &
          ipol = 1, 3)
     !
     !    here we control that the vectors are not on the same line
     !
     if ( (abs (e (1, 1)  * e (2, 2)  - e (2, 1)  * e (1, 2) ) .lt.1e-7 ) &
          .and. (abs (e (3, 1)  * e (1, 2)  - e (1, 1)  * e (3, 2) ) .lt.1e-7) &
          .and. (abs (e (3, 1)  * e (2, 2)  - e (2, 1)  * e (3, 2) ) .lt.1e-7) ) &
          call error ('chdens', 'vectors on the same line', 1)
     !
     !    and here that they are orthogonal
     !
     if (abs (e (1, 1)  * e (1, 2)  + e (2, 1)  * e (2, 2)  + e (3, 1) &
          & * e (3, 2) ) .gt.1e-4) call error ('chdens', &
          'vectors are not orthogonal', 1)
     !
     if (iflag.eq.3) then
        !
        !    reading for the 3D plot
        !
        read (inunit, *, err = 1100, iostat = ios) (e (ipol, 3), &
             ipol = 1, 3)

        !
        !    here we control that the vectors are not on the same line
        !
        if ( (abs (e (1, 1) * e (2, 3) - e (2, 1) * e (1, 3) ) &
             .lt.1e-7) .and. (abs (e (3, 1) * e (1, 3) - e (1, 1) &
             * e (3, 3) ) .lt.1e-7) .and. (abs (e (3, 1) * e (2, 3) &
             - e (2, 1) * e (3, 3) ) .lt.1e-7) ) call error ('chdens', &
             'vectors on the same line', 2)
        !
        !    and here that they are orthogonal
        !
        if (abs (e (1, 1) * e (1, 3) + e (2, 1) * e (2, 3) + e (3, &
             1) * e (3, 3) ) .gt.1e-4.or.abs (e (1, 2) * e (1, 3) &
             + e (2, 2) * e (2, 3) + e (3, 2) * e (3, 3) ) .gt.1e-4) &
             call error ('chdens', 'vectors are not orthogonal', 2)
     endif
     read (inunit, *, err = 1100, iostat = ios) (x0 (ipol), ipol = &
          1, 3)
     if (iflag.eq.2) then
        read (inunit, *, err = 1100, iostat = ios) nx, ny
     elseif (iflag.eq.3) then
        read (inunit, *, err = 1100, iostat = ios) nx, ny, nz
     endif
  endif

  if (iflag.eq.4) then
     read (inunit, *, err = 1100, iostat = ios) radius
     read (inunit, *, err = 1100, iostat = ios) nx, ny
  endif

  ! check for plot_out
  if (plot_out.lt.0.or.plot_out.gt.4) call error ('chdens', &
       'plot_out wrong', 1)

  !
  ! Read the header and allocate objects
  !
  call plot_io (filepp (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
       nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
       plot_num, atm, ityp, zv, tau, rhodum, 0)
  !
  allocate(tau (3, nat))
  allocate(ityp(nat))
  !
  alat = celldm (1)
  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  doublegrid = dual.gt.4.0d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif

  nspin = 1
  if (ibrav.gt.0) call latgen (ibrav, celldm, at (1, 1), &
                                            & at (1, 2), at (1, 3) )
  call recips (at (1, 1), at (1, 2), at (1, 3), bg (1, 1), bg (1, 2) &
       , bg (1, 3) )
  call volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)

  call set_fft_dim

  call allocate_fft

  rho = 0.d0
  !
  ! Read first file
  !
  call plot_io (filepp (1), title, nrx1, nrx2, nrx3, nr1, nr2, &
       nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
       plot_num, atm, ityp, zv, tau, rho, -1)
  !
  do ir = 1, nrxx
     psic (ir) = weight (1) * cmplx (rho (ir, 1),0.d0)
  enddo
  !
  ! Read following files (if any), verify consistency
  ! Note that only rho is read; all other quantities are discarded
  !
  do ifile = 2, nfile
     allocate  (taus( 3 , nat))    
     allocate  (ityps( nat))    
     !
     call plot_io (filepp (ifile), title, nrx1sa, nrx2sa, nrx3sa, &
          nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
          duals, ecuts, plot_num, atms, ityps, zvs, taus, rho, - 1)
     !
     deallocate (ityps)
     deallocate (taus)
     !
     if (nats.gt.nat) call error ('chdens', 'wrong file order? ', 1)
     if (nrx1.ne.nrx1sa.or.nrx2.ne.nrx2sa) call &
          error ('chdens', 'incompatible nrx1 or nrx2', 1)
     if (nr1.ne.nr1sa.or.nr2.ne.nr2sa.or.nr3.ne.nr3sa) call &
          error ('chdens', 'incompatible nr1 or nr2 or nr3', 1)
     if (ibravs.ne.ibrav) call error ('chdens', 'incompatible ibrav', 1)
     if (gcutmsa.ne.gcutm.or.duals.ne.dual.or.ecuts.ne.ecutwfc ) &
          call error ('chdens', 'incompatible gcutm or dual or ecut', 1)
     do i = 1, 6
        if (abs( celldm (i)-celldms (i) ) .gt. 1.0e-7 ) call error &
             ('chdens', 'incompatible celldm', 1)
     enddo
     !
     do ir = 1, nrxx
        psic (ir) = psic (ir) + weight (ifile) * cmplx (rho (ir, 1), &
             0.d0)
     enddo
  enddo

  !
  ! open output file, i.e., "fileout"
  !
  if (fileout.ne.' ') then
     ounit = 1
     open (unit = ounit, file = fileout, form = 'formatted', status &
          = 'unknown')
     write (6, '(5x,"Writing data on file ",a)') fileout
  else
     ounit = 6
  endif


  ! ----------------------------------------------------------------
  ! iflag=31
  ! ----------------------------------------------------------------
  !   at this point we are ready to print the whole FFT mesh (density only)
  !   TODO: check the consistency of iflag=31 with plot_out
  if (iflag.eq.31 .and. plot_out.eq.1) then

     if (output_format .ne. 3) then
        ! so far iflag.eq.31 works only for XCRYSDEN's XSF format
        call error ('chdens', 'wrong output_format for iflag.eq.31; works only if output_format.eq.3', 1)
     endif

     call plot_whole_cell (alat, at, nat, tau, atm, ityp, &
          nr1, nr2, nr3, nrx1, nrx2, nrx3, psic, output_format, ounit)
     return
  endif

  !
  !    At this point we start the calculations, first we normalize the vec
  !
  if (iflag.lt.4) then
     m1 = sqrt (e (1, 1) **2 + e (2, 1) **2 + e (3, 1) **2)
     if (iflag.ge.2) m2 = sqrt (e (1, 2) **2 + e (2, 2) **2 + e (3, &
          2) **2)

     if (iflag.eq.3) m3 = sqrt (e (1, 3) **2 + e (2, 3) **2 + e (3, &
          3) **2)

     do ipol = 1, 3
        e (ipol, 1) = e (ipol, 1) / m1
        if (iflag.ge.2) e (ipol, 2) = e (ipol, 2) / m2
        if (iflag.eq.3) e (ipol, 3) = e (ipol, 3) / m3
     enddo
  endif
  !
  !    and rebuild G-vectors in reciprocal space
  !
  call ggen
  !
  !    here we compute the fourier component of the quantity to plot
  !
  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  !    we store the fourier components in the array vgc
  !

  allocate (vgc( ngm))    
  if (plot_out.le.1) then
     do ig = 1, ngm
        vgc (ig) = psic (nl (ig) )
     enddo
  else
     ipol = plot_out - 1
     do ig = 2, ngm
        vgc (ig) = psic (nl (ig) ) * g (ipol, ig) / gg (ig)
     enddo

     vgc (1) = (epsilon - 1.d0) / fpi
     if (filepol.ne.' ') then
        !
        !    bring the quantity in real space and write the output file
        !
        call setv (2 * nrxx, 0.d0, psic, 1)
        do ig = 1, ngm
           psic (nl (ig) ) = vgc (ig)
        enddo
        call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        !
        do ir = 1, nrxx
           rho (ir, 1) = real (psic (ir) )
        enddo
        call plot_io (filepol, title, nrx1, nrx2, nrx3, nr1, nr2, &
             nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
             plot_num, atm, ityp, zv, tau, rho, + 1)
     endif
  endif

  !
  !     And now the plot
  !
  if (iflag.eq.1) then

     call plot_1d (nx, m1, x0, e, ngm, g, vgc, alat, plot_out, &
          ounit)

  elseif (iflag.eq.2) then

     call plot_2d (nx, ny, m1, m2, x0, e, ngm, g, vgc, alat, &
          at, nat, tau, atm, ityp, output_format, ounit)
     if (output_format.eq.2) then
        write (ounit, '(i4)') nat
        write (ounit, '(3f8.4,i3)') ( (tau (ipol, na) , ipol = 1, 3) &
             , 1, na = 1, nat)
        write (ounit, '(f10.6)') celldm (1)
        write (ounit, '(3(3f12.6/))') at
     endif

  elseif (iflag.eq.3) then
     if (output_format.eq.4) then

        ! gopenmol wants the coordinates in a separate file

        if (fileout .ne. ' ') then
           open (unit = ounit+1, file = trim(fileout)//'.xyz', &
                form = 'formatted', status = 'unknown')
           write (6, '(5x,"Writing coordinates on file ",a)') &
                trim(fileout)//'.xyz'
        else
           open (unit = ounit+1, file = 'coord.xyz', &
                form = 'formatted', status = 'unknown')
           write (6, '("Writing coordinates on file coord.xyz")')
        end if
     endif

     if (fast3d) then
        call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)
        call plot_fast (celldm (1), at, nat, tau, atm, ityp, &
             nrx1, nrx2, nrx3, nr1, nr2, nr3, psic, &
             bg, m1, m2, m3, x0, e, output_format, ounit)
     else
        call plot_3d (celldm (1), at, nat, tau, atm, ityp, ngm, g, vgc, &
             nx, ny, nz, m1, m2, m3, x0, e, output_format, ounit)
     end if

  elseif (iflag.eq.4) then
     radius = radius / alat

     call plot_2ds (nx, ny, radius, ngm, g, vgc, output_format, &
          ounit)
  else
     call error ('chdens', 'wrong iflag', 1)

  endif

  deallocate(vgc)
  return
1100 call error ('chdens', 'reading input data', abs (ios) )
end subroutine do_chdens
!
!-----------------------------------------------------------------------
subroutine plot_1d (nx, m1, x0, e, ngm, g, vgc, alat, plot_out, &
     ounit)
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  integer :: nx, ngm, plot_out, ounit
  ! number of points along the line
  ! number of G vectors
  ! type of plot
  ! output unit

  real(kind=DP) :: e (3), x0 (3), m1, alat, g (3, ngm)
  ! vector defining the line
  ! origin of the line
  ! modulus of e
  ! lattice parameter
  ! G-vectors

  complex(kind=DP) :: vgc (ngm)
  ! rho or polarization in G space
  integer :: i, ig
  real(kind=DP) :: rhomin, rhomax, rhoint, rhoim, xi, yi, zi, deltax, arg, gr
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated charge
  ! integrated imaginary charge
  ! coordinates of a 3D point
  ! steps along the line
  ! the argument of the exponential
  ! |G|*|r|

  real, parameter :: pi = 3.14159265358979d0
  complex(kind=DP) :: rho0g, carica (nx)

  deltax = m1 / (nx - 1)
  call setv (2 * nx, 0.d0, carica, 1)
  if (plot_out.eq.1) then
     do i = 1, nx
        xi = x0 (1) + (i - 1) * deltax * e (1)
        yi = x0 (2) + (i - 1) * deltax * e (2)
        zi = x0 (3) + (i - 1) * deltax * e (3)
        !
        !     for each point we compute the charge from the Fourier components
        !
        do ig = 1, ngm
           !
           !     NB: G are in 2pi/alat units, r are in alat units
           !
           arg = 2.d0 * pi * (xi * g (1, ig) + yi * g (2, ig) + zi * g (3, &
                ig) )
           carica (i) = carica (i) + vgc (ig) * cmplx (cos (arg), sin ( &
                arg) )
        enddo
     enddo
  else
     !
     !     spherically averaged charge: rho0(|r|) = int rho(r) dOmega
     !     rho0(r) = 4pi \sum_G rho(G) j_0(|G||r|)
     !
     !     G =0 term
     do i = 1, nx
        carica (i) = 4.d0 * pi * vgc (1)
     enddo
     !     G!=0 terms
     do ig = 2, ngm
        arg = 2.d0 * pi * (x0 (1) * g (1, ig) + x0 (2) * g (2, ig) &
             + x0 (3) * g (3, ig) )
        !     This displaces the origin into x0
        rho0g = vgc (ig) * cmplx (cos (arg), sin (arg) )
        !     r =0 term
        carica (1) = carica (1) + 4.d0 * pi * rho0g
        !     r!=0 terms
        do i = 2, nx
           gr = 2.d0 * pi * sqrt (g (1, ig) **2 + g (2, ig) **2 + g (3, &
                ig) **2) * (i - 1) * deltax
           carica (i) = carica (i) + 4.d0 * pi * rho0g * sin (gr) / gr
        enddo

     enddo
  endif
  !
  !    Here we check the value of the resulting charge
  !
  rhomin = 1.0d10
  rhomax = - 1.0d10

  rhoim = 0.d0
  do i = 1, nx
     rhomin = min (rhomin, dreal (carica (i) ) )
     rhomax = max (rhomax, dreal (carica (i) ) )
     rhoim = rhoim + abs (DIMAG (carica (i) ) )

  enddo

  rhoim = rhoim / nx
  print '(5x,"Min, Max, imaginary charge: ",3f12.6)', rhomin, &
       rhomax, rhoim
  !
  !       we print the charge on output
  !
  if (plot_out.eq.1) then
     do i = 1, nx
        write (ounit, '(2f20.10)') deltax * float (i - 1) , real ( &
             carica (i) )
     enddo
  else
     rhoint = 0.d0
     do i = 1, nx
        !
        !       simple trapezoidal rule: rhoint=int carica(i) r^2(i) dr
        !
        rhoint = rhoint + real (carica (i) ) * ( (i - 1) * deltax) **2 &
             * deltax * alat**3
        write (ounit, '(3f20.10)') deltax * float (i - 1) , real ( &
             carica (i) ) , rhoint
     enddo

  endif

  return

end subroutine plot_1d
!
!-----------------------------------------------------------------------
subroutine plot_2d (nx, ny, m1, m2, x0, e, ngm, g, vgc, alat, &
     at, nat, tau, atm, ityp, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  integer :: nx, ny, ngm, nat, ityp (nat), output_format, ounit
  ! number of points along x
  ! number of points along y
  ! number of G vectors
  ! number of atoms
  ! types of atoms
  ! output unit
  ! output format
  character(len=3) :: atm(*) ! atomic symbols
  real(kind=DP)    :: e (3, 2), x0 (3), m1, m2, g (3, ngm), &
       alat, tau (3, nat), at (3, 3)
  ! vectors e1, e2 defining the plane
  ! origin
  ! modulus of e1
  ! modulus of e2
  ! G-vectors

  complex(kind=DP) :: vgc (ngm)
  ! rho or polarization in G space
  integer :: i, j, ig

  real(kind=DP) :: rhomin, rhomax, rhoim, deltax, deltay
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
  real(kind=DP) :: pi

  parameter (pi = 3.14159265358979d0)
  complex(kind=DP), allocatable :: eigx (:), eigy (:), carica(:,:)

  allocate (eigx(  nx))    
  allocate (eigy(  ny))    
  allocate (carica( nx , ny))    

  deltax = m1 / (nx - 1)
  deltay = m2 / (ny - 1)

  call setv (2 * nx * ny, 0.d0, carica, 1)
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
     ! These factors are calculated and stored in order to save CPU time
     !
     do i = 1, nx
        eigx (i) = exp ( (0.d0, 1.d0) * 2.d0 * pi * ( (i - 1) * deltax * &
             (e(1,1) * g(1,ig) + e(2,1) * g(2,ig) + e(3,1) * g(3,ig) ) + &
             (x0 (1) * g(1,ig) + x0 (2) * g(2,ig) + x0 (3) * g(3,ig) ) ) )
     enddo
     do j = 1, ny
        eigy (j) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (j - 1) * deltay * &
             (e(1,2) * g(1,ig) + e(2,2) * g(2,ig) + e(3,2) * g(3,ig) ) )
     enddo
     do j = 1, ny
        do i = 1, nx
           carica (i, j) = carica (i, j) + vgc (ig) * eigx (i) * eigy (j)
        enddo
     enddo
  enddo
  !
  !    Here we check the value of the resulting charge
  !
  rhomin = 1.0d10
  rhomax = - 1.0d10

  rhoim = 0.d0
  do i = 1, nx
     do j = 1, ny
        rhomin = min (rhomin, dreal (carica (i, j) ) )
        rhomax = max (rhomax, dreal (carica (i, j) ) )
        rhoim = rhoim + abs (dimag (carica (i, j) ) )
     enddo

  enddo

  rhoim = rhoim / nx / ny
  print '(5x,"Min, Max, imaginary charge: ",3f12.6)', rhomin, &
       rhomax, rhoim
  print '(5x,"Output format: ",i3)', output_format

  !
  !     and we print the charge on output
  !
  if (output_format.eq.0) then
     !
     !     gnuplot format
     !
     !         write(ounit,'(2i6)') nx,ny
     do i = 1, nx
        write (ounit, '(e25.14)') (dreal (carica (i, j) ) , j = 1, ny)
        write (ounit, * )
     enddo
  elseif (output_format.eq.1) then
     !
     !     contour.x format
     !
     write (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
     write (ounit, '(4e25.14)') ( (dreal (carica (i, j) ) , j = 1, &
          ny) , i = 1, nx)
  elseif (output_format.eq.2) then
     !
     !     plotrho format
     !
     write (ounit, '(2i4)') nx - 1, ny - 1
     write (ounit, '(8f8.4)') (deltax * (i - 1) , i = 1, nx)
     write (ounit, '(8f8.4)') (deltay * (j - 1) , j = 1, ny)
     write (ounit, '(6e12.4)') ( (dreal (carica (i, j) ) , i = 1, &
          nx) , j = 1, ny)
     write (ounit, '(3f8.4)') x0
     write (ounit, '(3f8.4)') (m1 * e (i, 1) , i = 1, 3)
     write (ounit, '(3f8.4)') (m2 * e (i, 2) , i = 1, 3)

  elseif (output_format.eq.3) then
     !
     ! XCRYSDEN's XSF format
     !
     call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
     call xsf_datagrid_2d (carica, nx, ny, m1, m2, x0, e, alat, ounit)
  else
     call error('plot_2d', 'wrong output_format', 1)
  endif

  deallocate (carica)
  deallocate (eigy)
  deallocate (eigx)
  return
end subroutine plot_2d
!
!-----------------------------------------------------------------------
subroutine plot_2ds (nx, ny, x0, ngm, g, vgc, output_format, &
     ounit)
  !-----------------------------------------------------------------------
  use parameters, only : DP
  !
  implicit none
  integer :: nx, ny, ngm, ounit, output_format
  ! number of points along x
  ! number of points along y
  ! number of G vectors
  ! output unit

  real(kind=DP) :: x0, g (3, ngm)
  ! radius of the sphere
  ! G-vectors

  complex(kind=DP) :: vgc (ngm)
  ! rho or polarization in G space
  integer :: i, j, ig

  real(kind=DP), allocatable :: r (:,:,:)
  real(kind=DP) :: theta, phi, rhomin, rhomax, rhoim, &
       deltax, deltay
  ! the point in space
  ! the position on the sphere
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
  real(kind=DP) :: pi

  parameter (pi = 3.14159265358979d0)
  complex(kind=DP), allocatable :: carica (:,:)
  complex(kind=DP) :: eig

  allocate (carica( nx , ny))    
  allocate (r (3, nx , ny))    

  deltax = 2.d0 * pi / (nx - 1)

  deltay = pi / (ny - 1)

  call setv (2 * nx * ny, 0.d0, carica, 1)
  do j = 1, ny
     do i = 1, nx
        phi = (i - 1) * deltax
        theta = (j - 1) * deltay
        r (1, i, j) = x0 * sin (theta) * cos (phi)
        r (2, i, j) = x0 * sin (theta) * sin (phi)
        r (3, i, j) = x0 * cos (theta)
     enddo
  enddo
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
     ! These factors are calculated and stored in order to save CPU time
     !
     do j = 1, ny
        do i = 1, nx
           eig = exp ( (0.d0,1.d0) * 2.d0 * pi * (r (1, i, j) * g (1, ig) &
                + r (2, i, j) * g (2, ig) + r (3, i, j) * g (3, ig) ) )
           carica (i, j) = carica (i, j) + vgc (ig) * eig
        enddo
     enddo
  enddo
  !
  !    Here we check the value of the resulting charge
  !
  rhomin = 1.0d10
  rhomax = - 1.0d10

  rhoim = 0.d0
  do i = 1, nx
     do j = 1, ny
        rhomin = min (rhomin, dreal (carica (i, j) ) )
        rhomax = max (rhomax, dreal (carica (i, j) ) )
        rhoim = rhoim + abs (dimag (carica (i, j) ) )
     enddo

  enddo

  rhoim = rhoim / nx / ny
  print '(5x,"Min, Max, imaginary charge: ",3f12.6)', rhomin, &
       rhomax, rhoim
  !
  !     and we print the charge on output
  !
  if (output_format.eq.0) then
     !
     !     gnuplot format
     !
     write (ounit, '(2i8)') nx, ny
     do i = 1, nx
        write (ounit, '(e25.14)') (dreal (carica (i, j) ) , j = 1, ny)
     enddo
  elseif (output_format.eq.1) then
     !
     !     contour.x format
     !
     write (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
     write (ounit, '(4e25.14)') ( (dreal (carica (i, j) ) , j = 1, &
          ny) , i = 1, nx)
  else
     call error ('plot_2ds', 'not implemented plot', 1)

  endif
  deallocate (carica)
  deallocate (r)
  return

end subroutine plot_2ds
!
!-----------------------------------------------------------------------
subroutine plot_3d (alat, at, nat, tau, atm, ityp, ngm, g, vgc, &
     nx, ny, nz, m1, m2, m3, x0, e, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  integer :: nat, ityp (nat), ngm, nx, ny, nz, output_format, ounit
  ! number of atoms
  ! type of atoms
  ! number of G vectors
  ! number of points along x, y, z
  ! output format
  ! output unit
  character(len=3) :: atm(*)

  real(kind=DP) :: alat, tau (3, nat), at (3, 3), g (3, ngm), e (3, 3), &
       x0 (3), m1, m2, m3
  ! lattice parameter
  ! atomic positions
  ! lattice vectors
  ! G-vectors
  ! vectors e1,e2,e3 defining the parallelepiped
  ! origin
  ! moduli of e1,e2,e3

  complex(kind=DP) :: vgc (ngm)
  ! rho or polarization in G space
  integer :: i, j, k, ig

  real(kind=DP) :: rhomin, rhomax, rhotot, rhoabs, deltax, deltay, deltaz
  ! min, max value of the charge, total charge, total absolute charge
  ! steps along e1, e2, e3
  real(kind=DP), parameter :: pi = 3.14159265358979d0
  complex(kind=DP), allocatable :: eigx (:), eigy (:), eigz (:)
  real(kind=DP), allocatable :: carica (:,:,:)

  allocate (eigx(  nx))    
  allocate (eigy(  ny))    
  allocate (eigz(  nz))    
  allocate (carica( nx , ny , nz))    

  deltax = m1 / (nx - 1)
  deltay = m2 / (ny - 1)
  deltaz = m3 / (nz - 1)

  carica = 0.d0
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=exp(iG*e2), eigz=exp(iG*e3)
     ! These factors are calculated and stored in order to save CPU time
     !
     do i = 1, nx
        eigx (i) = exp ( (0.d0, 1.d0) * 2.d0 * pi * ( (i - 1) * deltax * &
             (e (1, 1) * g (1, ig) + e (2, 1) * g (2, ig) + e (3, 1) * g (3, ig)) &
             + (x0 (1) * g (1, ig) + x0 (2) * g (2, ig) + x0 (3) * g (3, ig) ) ) )
     enddo
     do j = 1, ny
        eigy (j) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (j - 1) * deltay * &
             (e (1, 2) * g (1, ig) + e (2, 2) * g (2, ig) + e (3, 2) * g (3, ig) ) )
     enddo
     do k = 1, nz
        eigz (k) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (k - 1) * deltaz * &
             (e (1, 3) * g (1, ig) + e (2, 3) * g (2, ig) + e (3, 3) * g (3, ig) ) )
     enddo
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              carica (i, j, k) = carica (i, j, k) + &
                   dreal (vgc (ig) * eigz (k) * eigy (j) * eigx (i) )
           enddo
        enddo
     enddo

  enddo
  !
  !    Here we check the value of the resulting charge
  !
  rhomin = 1.0d10
  rhomax =-1.0d10
  rhotot = 0.d0
  rhoabs = 0.d0
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           rhomin = min (rhomin, carica (i, j, k) )
           rhomax = max (rhomax, carica (i, j, k) )
           rhotot = rhotot + carica (i, j, k)
           rhoabs = rhoabs + abs (carica (i, j, k) )
        enddo
     enddo
  enddo

  rhotot = rhotot / nx / ny / nz * m1 * m2 * m3 * alat**3
  rhoabs = rhoabs / nx / ny / nz * m1 * m2 * m3 * alat**3
  print '(/5x,"Min, Max, Total, Abs charge: ",4f10.6)', rhomin, &
       rhomax, rhotot, rhoabs

  if (output_format.eq.4) then
     !
     ! "gOpenMol" file
     !

     call write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
          m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
  else
     ! user has calculated for very long, be nice and write some output even
     ! if the output_format is wrong; use XSF format as default

     !
     ! XCRYSDEN's XSF format
     !
     call xsf_struct      (alat, at, nat, tau, atm, ityp, ounit)
     call xsf_datagrid_3d (carica, nx, ny, nz, m1, m2, m3, x0, e, alat, ounit)
  endif

  deallocate (carica)
  deallocate (eigz)
  deallocate (eigy)
  deallocate (eigx)
  return
end subroutine plot_3d
!
!-----------------------------------------------------------------------
subroutine plot_fast (alat, at, nat, tau, atm, ityp, &
     nrx1, nrx2, nrx3, nr1, nr2, nr3, rho, &
     bg, m1, m2, m3, x0, e, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  integer :: nat, ityp (nat), nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       output_format, ounit
  character(len=3) :: atm(*)

  real(kind=DP) :: alat, tau (3, nat), at (3, 3), rho(2, nrx1,nrx2,nrx3), &
       bg (3, 3), e (3, 3), x0 (3), m1, m2, m3

  integer :: nx, ny, nz, nx0, ny0, nz0, nx1, ny1, nz1, i, j, k, i1, j1, k1
  real(kind=DP) :: rhomin, rhomax, rhotot, rhoabs
  real(kind=DP), allocatable :: carica (:,:,:)

  ! find FFT grid point closer to X0 (origin of the parallelepiped)
  ! (add 1 because r=0 correspond to n=1)

  nx0 = nint ( x0(1)*bg(1,1)*nr1 + x0(2)*bg(2,1)*nr1 + x0(3)*bg(3,1)*nr1 ) + 1
  ny0 = nint ( x0(1)*bg(1,2)*nr2 + x0(2)*bg(2,2)*nr2 + x0(3)*bg(3,2)*nr2 ) + 1
  nz0 = nint ( x0(1)*bg(1,3)*nr3 + x0(2)*bg(2,3)*nr3 + x0(3)*bg(3,3)*nr3 ) + 1
  !
  if ( e(2,1) .ne. 0.d0  .or.  e(3,1) .ne. 0.d0 .or. &
       e(1,2) .ne. 0.d0  .or.  e(3,2) .ne. 0.d0 .or. &
       e(1,3) .ne. 0.d0  .or.  e(2,3) .ne. 0.d0 )   &
       call error ('plot_fast','need vectors along x,y,z',1)

  ! find FFT grid points closer to X0 + e1, X0 + e2, X(0 + e3
  ! (the opposite vertex of the parallelepiped)

  nx1 = nint ((x0(1)+m1)*bg(1,1)*nr1 + x0(2)*bg(2,1)*nr1 + x0(3)*bg(3,1)*nr1 ) + 1
  ny1 = nint ( x0(1)*bg(1,2)*nr2 +(x0(2)+m2)*bg(2,2)*nr2 + x0(3)*bg(3,2)*nr2 ) + 1
  nz1 = nint ( x0(1)*bg(1,3)*nr3 + x0(2)*bg(2,3)*nr3 +(x0(3)+m3)*bg(3,3)*nr3 ) + 1

  nx = nx1 - nx0 + 1
  ny = ny1 - ny0 + 1
  nz = nz1 - nz0 + 1

  allocate (carica( nx, ny, nz))    

  carica = 0.d0
  do k = nz0, nz1
     k1 = mod(k-1, nr3) + 1
     if (k1.le.0) k1 = k1 + nr3
     do j = ny0, ny1
        j1 = mod(j-1, nr2) + 1
        if (j1.le.0) j1 = j1 + nr2
        do i = nx0, nx1
           i1 = mod(i-1,nr1)+1
           if (i1.le.0) i1 = i1 + nr1
           carica (i-nx0+1, j-ny0+1, k-nz0+1) = rho(1, i1, j1, k1)
        enddo
     enddo
  enddo
  !
  ! recalculate m1, m2, m3 (the sides of the parallelepiped divided by alat)
  ! consistent with the FFT grid
  !
  write(6,'(5x,"Requested parallelepiped sides : ",3f8.4)') m1, m2,m3
  m1 = (nx-1) * sqrt (at(1, 1) **2 + at(2, 1) **2 + at(3, 1) **2) / nr1
  m2 = (ny-1) * sqrt (at(1, 2) **2 + at(2, 2) **2 + at(3, 2) **2) / nr2
  m3 = (nz-1) * sqrt (at(1, 3) **2 + at(2, 3) **2 + at(3, 3) **2) / nr3
  write(6,'(5x,"Redefined parallelepiped sides : ",3f8.4)') m1, m2,m3
  !
  ! recalculate x0 (the origin of the parallelepiped)
  ! consistent with the FFT grid
  !
  write(6,'(5x,"Requested parallelepiped origin: ",3f8.4)') x0
  x0(1) = (nx0-1) * at(1,1)/nr1 + (ny0-1) * at(1,2)/nr2 + (nz0-1) * at(1,3)/nr3
  x0(2) = (nx0-1) * at(2,1)/nr1 + (ny0-1) * at(2,2)/nr2 + (nz0-1) * at(2,3)/nr3
  x0(3) = (nx0-1) * at(3,1)/nr1 + (ny0-1) * at(3,2)/nr2 + (nz0-1) * at(3,3)/nr3
  write(6,'(5x,"Redefined parallelepiped origin: ",3f8.4)') x0
  !
  !    Here we check the value of the resulting charge
  !
  rhomin = 1.0d10
  rhomax =-1.0d10
  rhotot = 0.d0
  rhoabs = 0.d0
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           rhomin = min (rhomin, carica (i, j, k) )
           rhomax = max (rhomax, carica (i, j, k) )
           rhotot = rhotot + carica (i, j, k)
           rhoabs = rhoabs + abs (carica (i, j, k) )
        enddo
     enddo
  enddo

  rhotot = rhotot / nx / ny / nz * m1 * m2 * m3 * alat**3
  rhoabs = rhoabs / nx / ny / nz * m1 * m2 * m3 * alat**3
  print '(/5x,"Min, Max, Total, Abs charge: ",4f10.6)', rhomin, &
       rhomax, rhotot, rhoabs

  if (output_format.eq.4) then
     !
     !     "gopenmol" file
     !

     call write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
          m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
  else
     !
     ! write XSF format
     !
     call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
     call xsf_datagrid_3d (carica, nx, ny, nz, m1, m2, m3, x0, e, alat, ounit)
  endif

  !
  deallocate (carica)
  return

end subroutine plot_fast
!-----------------------------------------------------------------------

subroutine write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
     m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
  !-----------------------------------------------------------------------

  use parameters, only : DP
  implicit none
  integer :: nat, ityp (nat), nx, ny, nz, ounit
  real(kind=DP) :: alat, tau (3, nat), at (3, 3), rhomax, x0 (3), &
       m1, m2, m3, carica (nx, ny, nz)
  character(len=3) :: atm(*)
  !
  integer, parameter :: MAXATOMS = 999
  real, parameter :: bohr = 0.529177
  integer :: natoms
  character(len=2) type (MAXATOMS)
  integer :: n1, n2, n3, na, i
  real(kind=DP) :: atoms (3, MAXATOMS), r (3), x, y, z
  real(kind=DP) :: sidex, sidey, sidez
  !
  !   sides of the parallelepiped in A
  !
  sidex = m1 * alat * bohr
  sidey = m2 * alat * bohr
  sidez = m3 * alat * bohr

  ! really bad algorithm to  generate (hopefully) all atoms
  ! that are inside the visualization box

  natoms = 0
  do n1 = - 3, + 3
     do n2 = - 3, + 3
        do n3 = - 3, + 3
           do i = 1, 3
              r (i) = n1 * at (i, 1) + n2 * at (i, 2) + n3 * at (i, 3)
           enddo
           do na = 1, nat
              ! x,y,z are in A
              x = (tau (1, na) + r (1) - x0 (1) ) * alat * bohr
              y = (tau (2, na) + r (2) - x0 (2) ) * alat * bohr
              z = (tau (3, na) + r (3) - x0 (3) ) * alat * bohr
              if ( x.gt.0d0 .and. x.lt.sidex .and. &
                   y.gt.0d0 .and. y.lt.sidey .and. &
                   z.gt.0d0 .and. z.lt.sidez) then
                 natoms = natoms + 1
                 if (natoms.gt.MAXATOMS) then
                    print '(" MAXATOMS (",i4,") Exceeded, " &
                         &       ,"Truncating " )', MAXATOMS
                    natoms = MAXATOMS
                    goto 10
                 endif
                 !
                 atoms (1, natoms) = x
                 atoms (2, natoms) = y
                 atoms (3, natoms) = z
                 !
                 type(natoms)=atm(ityp(na))
              endif
           enddo
        enddo
     enddo

  enddo

10 write(6,'(5x,"Found ",i4," atoms in the box")') natoms
  write(ounit,'("  3 2")')
  write(ounit,'(3i5)') nz,ny,nx
  write(ounit,'(6f10.4)') 0.0,sidez,0.0,sidey,0.0,sidex
  do n3=1,nz
     do n2 = 1, ny
        do n1 = 1, nx
           write (ounit, '(f20.10)') carica (n1, n2, n3)
        enddo
     enddo
  enddo
  !
  ! gopenmol needs atomic positions in a separate file
  !
  write(ounit+1,'(i4)') natoms
  write(ounit+1,'(2x,a2,3f9.4)') (type(na),( atoms(i,na), i=1,3 ), na=1,natoms )
  !
  return
end subroutine write_openmol_file
