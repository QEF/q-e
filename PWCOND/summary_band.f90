!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine summary_band(ik,ien)
!
! It gives a PWCOND summary on complex band structures of leads.
!
  USE io_global,  ONLY :  stdout
  USE noncollin_module, ONLY : npol
  use cond
  implicit none

  character(len=14) :: extension
  integer ::  i, k, irun, ik, ien, nstl, nstr
  real(DP) :: dim, dre, eev

  eev = earr(ien)

  nstl = n2d+npol*nocrosl
  nstr = n2d+npol*nocrosr
!
!  Output of complex bands in a separate file
!
  if(band_file.ne.' '.and.ikind.eq.0) then
      if(ik*ien.eq.1) then
        extension = '.re'
        open (3,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        extension = '.im'
        open (4,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        extension = '.co_re'
        open (11,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        extension = '.co_im'
        open (12,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        extension = '.3d'
        open (13,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        write(3,'("# Re(k), E-Ef")')
        write(4,'("# Im(k), E-Ef")')
        write(11,'("# Re(k), E-Ef")')
        write(12,'("# Im(k), E-Ef")')
        write(13,'("# Re(k), Im(k), E-Ef")')
      endif
      if(ien.eq.1) then
        write(3,'("# k-point", i5)') ik
        write(4,'("# k-point", i5)') ik
        write(11,'("# k-point", i5)') ik
        write(12,'("# k-point", i5)') ik
        write(13,'("# k-point", i5)') ik
      endif

!---
!     Propagating states

      do i = 1, nchanl
        write(3,'(2f10.4)') DBLE(kvall(i)), eev
        write(13,'(3f10.4)') DBLE(kvall(i)), AIMAG(kvall(i)), eev
        k = nstl + i
        write(3,'(2f10.4)') DBLE(kvall(k)), eev
        write(13,'(3f10.4)') DBLE(kvall(k)), AIMAG(kvall(k)), eev
      enddo
!---

!---
!     Evanescent states

      do k = nchanl+1, nstl
        do irun = 0, 1
          i = k + irun*nstl
          dre = DBLE(kvall(i))
          dim = abs(AIMAG(kvall(i)))
          if(dim.le.cutplot) then
            if(abs(dre).le.1.d-3)  then
              write(4,'(2f10.4)') -dim, eev
            elseif(abs(dre-0.5d0).le.1.d-3.or.abs(dre+0.5d0).le.1.d-3) then
              write(4,'(2f10.4)')  0.5d0+dim, eev
            else
              write(11,'(2f10.4)')  dre, eev
              write(12,'(2f10.4)') -0.5d0-dim, eev
            endif
            write(13,'(3f10.4)') DBLE(kvall(i)), AIMAG(kvall(i)), eev
          endif
        enddo
      enddo
!---

      if(ik*ien.eq.nkpts*nenergy) then
        close(unit=3)
        close(unit=4)
        close(unit=11)
        close(unit=12)
      endif
  endif
!
! Output of complex k onto common file
!
  WRITE( stdout,*) 'Nchannels of the left tip = ', nchanl
  WRITE( stdout,*) 'Right moving states:'
  WRITE( stdout,'(2x, a10, 2x, a10, 2x, a10)') 'k1(2pi/a)',  &
                               'k2(2pi/a)', 'E-Ef (eV)'
  do i = 1, nchanl
    WRITE( stdout,'(3f12.7)')  DBLE(kvall(i)), AIMAG(kvall(i)), eev
  enddo
  WRITE( stdout,*) 'Left moving states:'
  WRITE( stdout,'(2x, a10, 2x, a10, 2x, a10)') 'k1(2pi/a)',  &
                               'k2(2pi/a)', 'E-Ef (eV)'
  do i = nstl+1, nstl+nchanl
    WRITE( stdout,'(3f12.7)')  DBLE(kvall(i)), AIMAG(kvall(i)), eev
  enddo
  WRITE(stdout,*)
  if(ikind.eq.2) then
    WRITE( stdout,*) 'Nchannels of the right tip = ', nchanr
    WRITE( stdout,*) 'Right moving states:'
    WRITE( stdout,'(2x, a10, 2x, a10, 2x, a10)') 'k1(2pi/a)',  &
                                 'k2(2pi/a)', 'E-Ef (eV)'
    do i = 1, nchanr
      WRITE( stdout,'(3f12.7)')  DBLE(kvalr(i)), AIMAG(kvalr(i)), eev
    enddo
    WRITE( stdout,*) 'Left moving states:'
    WRITE( stdout,'(2x, a10, 2x, a10, 2x, a10)') 'k1(2pi/a)',  &
                               'k2(2pi/a)', 'E-Ef (eV)'
    do i = nstr+1, nstr+nchanr
      WRITE( stdout,'(3f12.7)')  DBLE(kvalr(i)), AIMAG(kvalr(i)), eev
    enddo
    WRITE(stdout,*)
  endif

  return
end subroutine summary_band
