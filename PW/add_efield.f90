! Copyright (C) 2003 J. Tobik 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
  subroutine add_efield
!--------------------------------------------------------------------------
!
!   This routine adds an electric field to the ionic potential. The
!   field is made artificially periodic by introducing a saw-tooth
!   potential. The field is parallel to a reciprocal lattice vector bg, 
!   according to the index edir.
!
!
#include "machine.h"
    use pwcom
#ifdef __PARA
    use para
#endif

  implicit none

  integer :: npoints, nmax, ndesc
  integer :: ii, ij, ik, itmp, ir, izlb, izub
  real(kind=dp) :: length, vamp, value
  real(kind=dp), parameter :: eps=1.d-8

#ifndef __PARA
  integer me, npp(1)
  me=1
  npp(1)=nr3
#endif
  if(edir.eq.1) then
     npoints=nr1
     length=alat/sqrt(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
  else if (edir.eq.2) then
     npoints=nr2
     length=alat/sqrt(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)
  elseif (edir.eq.3) then
     npoints=nr3
     length=alat/sqrt(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
  else
     call errore('setlocal',' wrong edir',1)
  endif

  nmax =int(real(npoints,dp)*(emaxpos-eps))+1
  if (nmax.lt.1.or.nmax.gt.npoints) &
     call errore('setlocal','nmax out of range',1)

  ndesc=int(real(npoints,dp)*(eopreg-eps))+1
  if (ndesc.lt.1.or.ndesc.gt.npoints) &
     call errore('setlocal','ndesc out of range',1)
!
!    The electric field is assumed in a.u.( 1 a.u. of field change the
!    potential energy of an electron of 1 Hartree in a distance of 1 Bohr. 
!    The factor 2 converts potential energy to Ry. 
!    
  vamp=2.0d0*eamp*length*real(npoints-ndesc,dp)/real(npoints,dp)

  write(6,*)
  write(6,*) 'Adding an external electric field'
  write(6,*) 'Intensity [a.u.]: ', eamp
  write(6,*) 'Amplitude vamp [Ry]', vamp
  write(6,*) 'Total length [points] ', npoints
  write(6,*) 'Total length [bohr rad] ', length
  write(6,*) 'Field is reversed between points ', nmax, nmax+ndesc
!
! in this case in x direction
!
  if(edir.eq.1) then
    do ij=1,nr2
      do ik=1,npp(me)
        do ii=nmax,nmax+ndesc-1
           value=vamp*(real(nmax+ndesc-ii,dp)/real(ndesc,dp)-0.5d0)
           itmp=ii
           if (itmp.gt.nr1) itmp=itmp-nr1
           ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
           vltot(ir)=vltot(ir)+value
        end do
        do ii=nmax+ndesc,nmax+nr1-1
           value=vamp*(real(ii-nmax-ndesc,dp)/real(nr1-ndesc,dp)-0.5d0)
           itmp=ii
           if (itmp.gt.nr1) itmp=itmp-nr1
           ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
           vltot(ir)=vltot(ir)+value
        end do
      end do
    end do
!
! in this case in y direction
!
  else if (edir.eq.2) then
    do ii=1,nr1
      do ik=1,npp(me)
        do ij=nmax,nmax+ndesc-1
           value=vamp*(real(nmax+ndesc-ij,dp)/real(ndesc,dp)-0.5d0)
           itmp=ij
           if (itmp.gt.nr2) itmp=itmp-nr2
           ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
           vltot(ir)=vltot(ir)+value
        end do
        do ij=nmax+ndesc,nmax+nr2-1
           value=vamp*(real(ij-nmax-ndesc,dp)/real(nr2-ndesc,dp)-0.5d0)
           itmp=ij
           if (itmp.gt.nr2) itmp=itmp-nr2
           ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
           vltot(ir)=vltot(ir)+value
        end do
      end do
    end do
!
! and in other cases in z direction
!
  elseif (edir.eq.3) then
#ifdef __PARA
    izub=0
    do itmp=1,me
       izlb=izub+1
       izub=izub+npp(itmp)
    end do
#else
    izlb=1
    izub=nr3
#endif
!
!  now we have set up boundaries - let's calculate potential
!
    do ii=1,nr1
      do ij=1,nr2
        do ik=nmax,nmax+ndesc-1
          value=vamp*(real(nmax+ndesc-ik,dp)/real(ndesc,dp)-0.5d0)
          itmp=ik
          if (itmp.gt.nr3) itmp=itmp-nr3
          if((itmp.ge.izlb).and.(itmp.le.izub)) then
!
! Yes - this point belongs to me
!
            itmp=itmp-izlb+1
            ir=ii+(ij-1)*nrx1+(itmp-1)*nrx1*nrx2
            vltot(ir)=vltot(ir)+value
          end if
        end do
        do ik=nmax+ndesc,nmax+nr3-1
           value=vamp*(real(ik-nmax-ndesc,dp)/real(nr3-ndesc,dp)-0.5d0)
           itmp=ik
           if (itmp.gt.nr3) itmp=itmp-nr3
           if((itmp.ge.izlb).and.(itmp.le.izub)) then
              itmp=itmp-izlb+1
              ir=ii+(ij-1)*nrx1+(itmp-1)*nrx1*nrx2
              vltot(ir)=vltot(ir)+value
           end if
        end do
      end do
    end do
  else
     call errore('setlocal', 'wrong edir', 1)
  endif
  return
end subroutine add_efield

