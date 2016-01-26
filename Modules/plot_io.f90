!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine plot_io (filplot, title, nr1x, nr2x, nr3x, nr1, nr2, &
     nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecut, plot_num, atm, &
     ityp, zv, tau, plot, iflag)
  !-----------------------------------------------------------------------
  !
  !     iflag >0 : write header and the quantity to be plotted ("plot")
  !                to file "filplot"
  !     iflag< 0 : read everything (requires that all variables that are
  !                read are allocated with the correct dimensions!)
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  character (len=*) :: filplot
  character (len=75) :: title
!   integer :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav, &
!        plot_num, ityp (nat), iflag, i
  integer :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav, &
       plot_num, ityp (*), iflag, i
  character (len=3) :: atm(*)
!   real(DP) :: celldm (6), gcutm, dual, ecut, zv (ntyp), tau (3, nat) &
!        , plot (nr1x * nr2x * nr3x), at(3,3)
  real(DP) :: celldm (6), gcutm, dual, ecut, zv (*), tau (3, *) &
       , plot (*), at(3,3)
  !
  integer :: iunplot, ios, ipol, na, nt, ir, ndum
  !
  if (filplot == ' ') call errore ('plot_io', 'filename missing', 1)
  !
  iunplot = 4
  if (iflag == 0 ) call errore('plot_io',&
                        ' iflag==0 not allowed, use read_io_header ',1)
  if (iflag > 0) then
     WRITE( stdout, '(5x,"Writing data to file  ",a)') TRIM(filplot)
     open (unit = iunplot, file = filplot, form = 'formatted', &
          status = 'unknown', err = 100, iostat = ios)
  else
     WRITE( stdout, '(5x,"Reading data from file  ",a)') TRIM(filplot)
     open (unit = iunplot, file = filplot, form = 'formatted', &
          status = 'old', err = 100, iostat = ios)
  endif

100 call errore ('plot_io', 'opening file '//TRIM(filplot), abs (ios) )

  rewind (iunplot)
  if (iflag > 0) then
     write (iunplot, '(a)') title
     write (iunplot, '(8i8)') nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
     write (iunplot, '(i6,2x,6f16.8)') ibrav, celldm
     if (ibrav == 0) then
        do i = 1,3
           write ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
        enddo
     endif
     write (iunplot, '(3f20.10,i6)') gcutm, dual, ecut, plot_num
     write (iunplot, '(i4,3x,a2,3x,f5.2)') &
          (nt, atm (nt), zv (nt), nt=1, ntyp)
     write (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
          (tau (ipol, na), ipol = 1, 3), ityp (na), na = 1, nat)
     write (iunplot, '(5(1pe17.9))') (plot (ir) , ir = 1, nr1x * nr2x * nr3)
  else
     read (iunplot, '(a)') title
     read (iunplot, * ) nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
     read (iunplot, * ) ibrav, celldm
     if (ibrav == 0) then
        do i = 1,3
           read ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
        enddo
     endif
     read (iunplot, * ) gcutm, dual, ecut, plot_num
     read (iunplot, '(i4,3x,a2,3x,f5.2)') &
             (ndum, atm(nt), zv(nt), nt=1, ntyp)
     read (iunplot, *) (ndum,  (tau (ipol, na), ipol = 1, 3), &
             ityp(na), na = 1, nat)
     read (iunplot, * ) (plot (ir), ir = 1, nr1x * nr2x * nr3)
  endif

  close (unit = iunplot)

  return
end subroutine plot_io
!-----------------------------------------------------------------------
subroutine read_io_header(filplot, title, nr1x, nr2x, nr3x, nr1, nr2, nr3, &
                  nat, ntyp, ibrav, celldm, at, gcutm, dual, ecut, plot_num)
     
  !-----------------------------------------------------------------------
  !
  !     read header of file "filplot" 
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  character (len=*) :: filplot
  character (len=75) :: title
  integer :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav, plot_num, i
  real(DP) :: celldm (6), gcutm, dual, ecut, at(3,3)
  !
  integer :: iunplot, ios, ipol
  !
  if (filplot == ' ') call errore ('read_io_h', 'filename missing', 1)
  !
  iunplot = 4
  WRITE( stdout, '(5x,"Reading header from file  ",a)') TRIM(filplot)
  open (unit = iunplot, file = filplot, form = 'formatted', &
          status = 'old', err = 100, iostat = ios)
100 call errore ('plot_io', 'opening file '//TRIM(filplot), abs (ios) )

  rewind (iunplot)
  read (iunplot, '(a)') title
  read (iunplot, * ) nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
  read (iunplot, * ) ibrav, celldm
  if (ibrav == 0) then
     do i = 1,3
        read ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
     enddo
  endif
  read (iunplot, * ) gcutm, dual, ecut, plot_num
  close (unit = iunplot)

  return
end subroutine read_io_header
