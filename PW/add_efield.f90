! Copyright (C) 2003 J. Tobik 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
!
!--------------------------------------------------------------------------
  subroutine add_efield(vpoten)
!--------------------------------------------------------------------------
!
!   This routine adds an electric field to the local potential. The
!   field is made artificially periodic by introducing a saw-tooth
!   potential. The field is parallel to a reciprocal lattice vector bg, 
!   according to the index edir.
!
!   if dipfield is false the electric field correction is added to the
!   potential given as input (the bare local potential) only
!   at the first call to this routine. In the following calls
!   the routine exit.
!
!   if dipfield is true the dipole moment per unit surface is calculated
!   and used to cancel the electric field due to periodic boundary
!   conditions. This potential is added to the Hartree and xc potential
!   in v_of_rho. NB: in this case the electric field contribution to the 
!   band energy is subtracted by deband.
!
!
#include "machine.h"
    USE kinds, ONLY : DP
    USE constants, ONLY: fpi
    USE basis, ONLY : nat, ityp, zv
    USE cell_base, ONLY : alat, bg, omega
    USE extfield, ONLY: tefield, dipfield, edir, eamp, emaxpos, eopreg, &
         etotefield, forcefield
    USE force_mod, ONLY: lforce
    USE gvect, ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
    USE io_global,  ONLY : stdout
    USE control_flags, ONLY: mixing_beta
#ifdef __PARA
    USE mp_global,     ONLY : intra_image_comm
    use para
    use mp
#endif

  implicit none

  integer :: npoints, nmax, ndesc
  integer :: ii, ij, ik, itmp, ir, izlb, izub, na, ipol, n3
  real(kind=dp) :: length, vamp, value
  real(kind=dp), parameter :: eps=1.d-8
  real(kind=dp) :: vpoten(nrxx) ! the ef is added to this potential
  real(kind=dp) :: dip, dipion, bmod, z0, dipold
  real(kind=dp) :: deltal

  logical :: first=.true.
  save first, dipold

#ifndef __PARA
  integer me, npp(1)
  me=1
  npp(1)=nr3
#endif
  

  if (.not.tefield) return
  if ((.not.dipfield).and. (.not.first)) return

  bmod=sqrt(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)
  if(edir.eq.1) then
     npoints=nr1
  else if (edir.eq.2) then
     npoints=nr2
  elseif (edir.eq.3) then
     npoints=nr3
  else
     call errore('add_efield',' wrong edir',1)
  endif
  length=alat/bmod
  deltal=length/npoints

  nmax =int(real(npoints,dp)*(emaxpos-eps))+1
  if (nmax.lt.1.or.nmax.gt.npoints) &
     call errore('add_efield','nmax out of range',1)

  ndesc=int(real(npoints,dp)*(eopreg-eps))+1
  if (ndesc.lt.1.or.ndesc.gt.npoints) &
     call errore('add_efield','ndesc out of range',1)

  dip=0.d0
  dipion=0.d0
  n3=nmax+ndesc+(nr3-ndesc)/2
  if (n3.gt.nr3) n3=n3-nr3
  z0=(n3-1)*deltal
  if (mod(nr3-ndesc,2).ne.0) z0=z0+deltal*0.5d0
  z0=z0/alat

  if (first.and.dipfield) z0=0.d0
  call compute_dip(dip,dipion,z0)
!
!  This is used to reach self-consistency. Mixing the dipole field improves
!  convergence. 
!
  if (first) then
     dipold=dip
  else
     dip=dip*mixing_beta+dipold*(1.d0-mixing_beta)
     dipold=dip
  endif
#ifdef __PARA
  call mp_bcast(dip,0,intra_image_comm)
#endif
  if (.not.dipfield) then
     etotefield=-2.d0*dipion*eamp*omega/fpi 
     dip=0.d0
  else
     etotefield=-2.d0*(eamp-dip/2.d0)*dip*omega/fpi 
  endif

  if (lforce) then
     do na=1,nat
        do ipol=1,3
           forcefield(ipol,na)=2.d0*(eamp-dip) &
                                    *zv(ityp(na))*bg(ipol,edir)/bmod
        enddo
     enddo
  endif
     
!
!    The electric field is assumed in a.u.( 1 a.u. of field changes the
!    potential energy of an electron of 1 Hartree in a distance of 1 Bohr. 
!    The factor 2 converts potential energy to Ry. 
!    NB: dip is the dipole moment per unit area divided by length and
!        multiplied by four pi
!    
  vamp=2.0d0*(eamp-dip)*length*real(npoints-ndesc,dp)&
                                               /real(npoints,dp)
  if (first) then
     WRITE( stdout,*)
     WRITE( stdout,'(5x,"Adding an external electric field")')
     WRITE( stdout,'(5x,"Intensity [a.u.]: ",f15.8)') eamp
  endif
  if (dipfield) WRITE( stdout,'(5x,"Dipole field [a.u.]: ", f15.8)') dip
  if (first) then
     WRITE( stdout,'(5x,"Potential amplitude [Ry]: ", f15.8)') vamp
     WRITE( stdout,'(5x,"Total length [points]: ", i5)') npoints
     WRITE( stdout,'(5x,"Total length [bohr rad]: ", f15.8)') length
     WRITE( stdout,'(5x,"Field is reversed between points: ",2i6)')nmax, nmax+ndesc
  endif
!
! in this case x direction
!
  if(edir.eq.1) then
    do ij=1,nr2
      do ik=1,npp(me)
        do ii=nmax,nmax+ndesc-1
           value=vamp*(real(nmax+ndesc-ii,dp)/real(ndesc,dp)-0.5d0)
           itmp=ii
           if (itmp.gt.nr1) itmp=itmp-nr1
           ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
           vpoten(ir)=vpoten(ir)+value
        end do
        do ii=nmax+ndesc,nmax+nr1-1
           value=vamp*(real(ii-nmax-ndesc,dp)/real(nr1-ndesc,dp)-0.5d0)
           itmp=ii
           if (itmp.gt.nr1) itmp=itmp-nr1
           ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
           vpoten(ir)=vpoten(ir)+value
        end do
      end do
    end do
!
! in this case y direction
!
  else if (edir.eq.2) then
    do ii=1,nr1
      do ik=1,npp(me)
        do ij=nmax,nmax+ndesc-1
           value=vamp*(real(nmax+ndesc-ij,dp)/real(ndesc,dp)-0.5d0)
           itmp=ij
           if (itmp.gt.nr2) itmp=itmp-nr2
           ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
           vpoten(ir)=vpoten(ir)+value
        end do
        do ij=nmax+ndesc,nmax+nr2-1
           value=vamp*(real(ij-nmax-ndesc,dp)/real(nr2-ndesc,dp)-0.5d0)
           itmp=ij
           if (itmp.gt.nr2) itmp=itmp-nr2
           ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
           vpoten(ir)=vpoten(ir)+value
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
            vpoten(ir)=vpoten(ir)+value
          end if
        end do
        do ik=nmax+ndesc,nmax+nr3-1
           value=vamp*(real(ik-nmax-ndesc,dp)/real(nr3-ndesc,dp)-0.5d0)
           itmp=ik
           if (itmp.gt.nr3) itmp=itmp-nr3
           if((itmp.ge.izlb).and.(itmp.le.izub)) then
              itmp=itmp-izlb+1
              ir=ii+(ij-1)*nrx1+(itmp-1)*nrx1*nrx2
              vpoten(ir)=vpoten(ir)+value
           end if
        end do
      end do
    end do
  else
     call errore('add_efield', 'wrong edir', 1)
  endif
  first=.false.
  return
end subroutine add_efield

