!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine summary_band(ik,ien)
!
! It gives a PWCOND summary and does some final transmittance
! calculations.
! 
#include "f_defs.h"
  USE io_global,  ONLY :  stdout
  USE io_files, ONLY: band_file
  use pwcom
  use cond
  implicit none

  character(len=14) :: extension
  integer ::  i, ik, ien, nstl, nstr 
  real(kind=DP) :: dim, dre, eev

  eev = earr(ien)

  nstl = n2d+nocrosl
  nstr = n2d+nocrosr
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
        extension = '.co'
        open (11,file=trim(band_file)//extension,form='formatted', &
                      status='unknown')
        write(3,'("# Re(k), E-Ef")')
        write(4,'("# Im(k), E-Ef")')
        write(11,'("# Re (Im(k)), E-Ef")')
      endif
      if(ien.eq.1) then
        write(3,'("# k-point", i5)') ik
        write(4,'("# k-point", i5)') ik
        write(11,'("# k-point", i5)') ik 
      endif

      do i=1, nchanl
       write(3,'(2f10.4)') abs(DREAL(kvall(i))), eev
      enddo

      do i=nchanl+1, nstl
       dre = abs(DREAL(kvall(i)))
       dim = abs(DIMAG(kvall(i)))
       if(dim.le.cutplot) then
        if(dre.gt.1.d-4.and.dre.le.0.4999d0) then                    
          write(11,'(2f10.4)')  dre, eev
          write(11,'(2f10.4)') -dim, eev
        else
          if(dre.le.0.1d0) write(4,'(2f10.4)') -dim, eev  
          if(dre.gt.0.4d0) write(4,'(2f10.4)')  0.5d0+dim, eev
        endif
       endif
      enddo

      if(ik*ien.eq.nkpts*nenergy) then
        close(unit=3)
        close(unit=4)
        close(unit=11)
      endif
  endif
!
! Output of complex k onto common file
!
  WRITE( stdout,*) 'Nchannels of the left tip = ', nchanl
  WRITE( stdout,'(7x, a10, 3x, a10, 6x, a10)') 'k1(2pi/a)',  &
                               'k2(2pi/a)', 'E-Ef (eV)'
  WRITE( stdout,*)
  do i=1, nchanl
    WRITE( stdout,'(3f12.7)') DREAL(kvall(i)), DIMAG(kvall(i)), eev
  enddo
  if(ikind.eq.2) then
    WRITE( stdout,*) 'Nchannels of the right tip = ', nchanr
    WRITE( stdout,'(7x, a10, 3x, a10, 6x, a10)') 'k1(2pi/a)',&
                               'k2(2pi/a)', 'E-Ef (eV)'          
    WRITE( stdout,*)
    do i=1, nchanr
      WRITE( stdout,'(3f12.7)') DREAL(kvalr(i)), DIMAG(kvalr(i)), eev
    enddo         
  endif

  return
end subroutine summary_band            

