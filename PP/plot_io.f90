!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine plot_io (filplot, title, nrx1, nrx2, nrx3, nr1, nr2, &
     nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecut, plot_num, atm, &
     ityp, zv, tau, plot, iflag)
  !-----------------------------------------------------------------------
  !
  !     iflag >0 : write header and the quantity to be plotted ("plot")
  !                on file "filplot"
  !     iflag= 0 : read only the header and not "plot"
  !     iflag< 0 : read everything (requires that all variables that are
  !                read are allocated with the correct dimensions!)
  !
  use parameters, only : DP
  implicit none
  character (len=14) :: filplot
  character (len=75) :: title
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, ntyp, ibrav, &
       plot_num, ityp (ntyp), iflag, i
  character (len=3) :: atm(ntyp)
  real(kind=DP) :: celldm (6), gcutm, dual, ecut, zv (ntyp), tau (3, nat) &
       , plot (nrx1 * nrx2 * nrx3), at(3,3)
  !
  integer :: iunplot, ios, ipol, na, nt, ir, ndum

  if (filplot.eq.' ') call error ('plot_io', 'filename missing', 1)
  !
  iunplot = 4
  if (iflag.gt.0) then
     write (6, '(5x,"Writing data on file  ",a)') filplot
     open (unit = iunplot, file = filplot, form = 'formatted', &
          status = 'unknown', err = 100, iostat = ios)
  else
     if (iflag.lt.0) then
        write (6, '(5x,"Reading data from file  ",a)') filplot
     else
        write (6, '(5x,"Reading header from file  ",a)') filplot
     endif
     open (unit = iunplot, file = filplot, form = 'formatted', &
          status = 'old', err = 100, iostat = ios)
  endif

100 call error ('plot_io', 'opening file '//filplot, abs (ios) )

  rewind (iunplot)
  if (iflag.gt.0) then
     write (iunplot, '(a)') title
     write (iunplot, '(8i8)') nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, &
          ntyp
     write (iunplot, '(i6,6f12.8)') ibrav, celldm
     if (ibrav.eq.0) then
        do i = 1,3
           write ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
        enddo
     endif
     write (iunplot, '(3f20.10,i6)') gcutm, dual, ecut, plot_num
     write (iunplot, '(i4,3x,a2,3x,f5.2)') &
          (nt, atm (nt), zv (nt), nt=1, ntyp)
     write (iunplot, '(i4,3x,3f14.10,3x,i2)') (na, &
          (tau (ipol, na), ipol = 1, 3), ityp (na), na = 1, nat)
     if (plot_num.ne.9) write (iunplot, '(5(1pe17.9))') (plot (ir) , &
          ir = 1, nrx1 * nrx2 * nr3)
  else
     read (iunplot, '(a)') title
     read (iunplot, * ) nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, ntyp
     read (iunplot, * ) ibrav, celldm
     if (ibrav.eq.0) then
        do i = 1,3
           read ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
        enddo
     endif
     read (iunplot, * ) gcutm, dual, ecut, plot_num
     if (iflag.lt.0) then
        read (iunplot, '(i4,3x,a2,3x,f5.2)') &
             (ndum, atm(nt), zv(nt), nt=1, ntyp)
        read (iunplot, *) (ndum,  (tau (ipol, na), ipol = 1, 3), &
             ityp(na), na = 1, nat)
        read (iunplot, * ) (plot (ir), ir = 1, nrx1 * nrx2 * nr3)
     endif
  endif

  if (plot_num.ne.9) close (unit = iunplot)
  return
end subroutine plot_io
