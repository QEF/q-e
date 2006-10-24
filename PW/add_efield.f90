!
! Copyright (C) 2003-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by J. Tobik
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
#include "f_defs.h"
!
!--------------------------------------------------------------------------
SUBROUTINE add_efield(rho,vpoten,etotefield,iflag)
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
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, bg, omega
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
  USE force_mod,     ONLY : lforce
  USE gvect,         ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_global,     ONLY : intra_image_comm, me_pool
  USE pfft,          ONLY : npp
  USE mp,            ONLY : mp_bcast

  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP) :: rho(nrxx,nspin) ! the density whose dipole is computed
  REAL(DP) :: vpoten(nrxx) ! the ef is added to this potential
  REAL(DP) :: etotefield   ! the contribution to etot due to ef

  INTEGER :: iflag
  !
  ! local variables
  !
  INTEGER :: npoints, nmax, ndesc
  INTEGER :: ii, ij, ik, itmp, ir, izlb, izub, na, ipol, n3
  REAL(DP) :: length, vamp, value
  REAL(DP) :: dip, dipion, bmod, z0
  REAL(DP) :: deltal

  LOGICAL :: first=.TRUE.
  SAVE first

#ifndef __PARA
  npp(1) = nr3
#endif


  IF (.NOT.tefield) RETURN
  IF ((.NOT.dipfield).AND. (.NOT.first).and. (iflag.EQ.0)) RETURN

  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)
  IF(edir.EQ.1) THEN
     npoints=nr1
  ELSE IF (edir.EQ.2) THEN
     npoints=nr2
  ELSEIF (edir.EQ.3) THEN
     npoints=nr3
  ELSE
     CALL errore('add_efield',' wrong edir',1)
  ENDIF
  length=alat/bmod
  deltal=length/npoints

  nmax =INT(REAL(npoints,dp)*(emaxpos-eps8))+1
  IF (nmax.LT.1.OR.nmax.GT.npoints) &
       CALL errore('add_efield','nmax out of range',1)

  ndesc=INT(REAL(npoints,dp)*(eopreg-eps8))+1
  IF (ndesc.LT.1.OR.ndesc.GT.npoints) &
       CALL errore('add_efield','ndesc out of range',1)

  dip=0.d0
  dipion=0.d0
  n3=nmax+ndesc+(npoints-ndesc)/2
  IF (n3.GT.nmax+ndesc) n3=n3-npoints
  z0=(n3-1)*deltal
  IF (MOD(npoints-ndesc,2).NE.0) z0=z0+deltal*0.5d0
  z0=z0/alat

  IF (first.AND.dipfield) z0=0.d0
  CALL compute_dip(rho,dip,dipion,z0)
  !
#ifdef __PARA
  CALL mp_bcast(dip,0,intra_image_comm)
#endif
  IF (.NOT.dipfield) THEN
     etotefield=-2.d0*dipion*eamp*omega/fpi 
     dip=0.d0
  ELSE
     etotefield=-2.d0*(eamp-dip/2.d0)*dip*omega/fpi 
  ENDIF

  IF (lforce) THEN
     DO na=1,nat
        DO ipol=1,3
           forcefield(ipol,na)=2.d0*(eamp-dip) &
                *zv(ityp(na))*bg(ipol,edir)/bmod
        ENDDO
     ENDDO
  ENDIF

  !
  !    The electric field is assumed in a.u.( 1 a.u. of field changes the
  !    potential energy of an electron of 1 Hartree in a distance of 1 Bohr. 
  !    The factor 2 converts potential energy to Ry. 
  !    NB: dip is the dipole moment per unit area divided by length and
  !        multiplied by four pi
  !    
  vamp=2.0d0*(eamp-dip)*length*REAL(npoints-ndesc,dp)&
       /REAL(npoints,dp)
  IF (first) THEN
     WRITE( stdout,*)
     WRITE( stdout,'(5x,"Adding an external electric field")')
     WRITE( stdout,'(5x,"Intensity [a.u.]: ",f15.8)') eamp
  ENDIF
  IF (dipfield) WRITE( stdout,'(5x,"Dipole field [a.u.]: ", f15.8)') dip
  IF (first) THEN
     WRITE( stdout,'(5x,"Potential amplitude [Ry]: ", f15.8)') vamp
     WRITE( stdout,'(5x,"Total length [points]: ", i5)') npoints
     WRITE( stdout,'(5x,"Total length [bohr rad]: ", f15.8)') length
     WRITE( stdout,'(5x,"Field is reversed between points: ",2i6)')nmax, nmax+ndesc
  ENDIF
  !
  ! in this case x direction
  !
  IF(edir.EQ.1) THEN
     DO ij=1,nr2
        DO ik=1,npp(me_pool + 1)
           DO ii=nmax,nmax+ndesc-1
              value=vamp*(REAL(nmax+ndesc-ii,dp)/REAL(ndesc,dp)-0.5d0)
              itmp=ii
              IF (itmp.GT.nr1) itmp=itmp-nr1
              ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
              vpoten(ir)=vpoten(ir)+value
           END DO
           DO ii=nmax+ndesc,nmax+nr1-1
              value=vamp*(REAL(ii-nmax-ndesc,dp)/REAL(nr1-ndesc,dp)-0.5d0)
              itmp=ii
              IF (itmp.GT.nr1) itmp=itmp-nr1
              ir=itmp+(ij-1)*nrx1+(ik-1)*nrx1*nrx2
              vpoten(ir)=vpoten(ir)+value
           END DO
        END DO
     END DO
     !
     ! in this case y direction
     !
  ELSE IF (edir.EQ.2) THEN
     DO ii=1,nr1
        DO ik=1,npp(me_pool + 1)
           DO ij=nmax,nmax+ndesc-1
              value=vamp*(REAL(nmax+ndesc-ij,dp)/REAL(ndesc,dp)-0.5d0)
              itmp=ij
              IF (itmp.GT.nr2) itmp=itmp-nr2
              ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
              vpoten(ir)=vpoten(ir)+value
           END DO
           DO ij=nmax+ndesc,nmax+nr2-1
              value=vamp*(REAL(ij-nmax-ndesc,dp)/REAL(nr2-ndesc,dp)-0.5d0)
              itmp=ij
              IF (itmp.GT.nr2) itmp=itmp-nr2
              ir=ii+(itmp-1)*nrx1+(ik-1)*nrx1*nrx2
              vpoten(ir)=vpoten(ir)+value
           END DO
        END DO
     END DO
     !
     ! and in other cases in z direction
     !
  ELSEIF (edir.EQ.3) THEN
#ifdef __PARA
     izub=0
     DO itmp=1,me_pool + 1
        izlb=izub+1
        izub=izub+npp(itmp)
     END DO
#else
     izlb=1
     izub=nr3
#endif
     !
     !  now we have set up boundaries - let's calculate potential
     !
     DO ii=1,nr1
        DO ij=1,nr2
           DO ik=nmax,nmax+ndesc-1
              value=vamp*(REAL(nmax+ndesc-ik,dp)/REAL(ndesc,dp)-0.5d0)
              itmp=ik
              IF (itmp.GT.nr3) itmp=itmp-nr3
              IF((itmp.GE.izlb).AND.(itmp.LE.izub)) THEN
                 !
                 ! Yes - this point belongs to me
                 !
                 itmp=itmp-izlb+1
                 ir=ii+(ij-1)*nrx1+(itmp-1)*nrx1*nrx2
                 vpoten(ir)=vpoten(ir)+value
              END IF
           END DO
           DO ik=nmax+ndesc,nmax+nr3-1
              value=vamp*(REAL(ik-nmax-ndesc,dp)/REAL(nr3-ndesc,dp)-0.5d0)
              itmp=ik
              IF (itmp.GT.nr3) itmp=itmp-nr3
              IF((itmp.GE.izlb).AND.(itmp.LE.izub)) THEN
                 itmp=itmp-izlb+1
                 ir=ii+(ij-1)*nrx1+(itmp-1)*nrx1*nrx2
                 vpoten(ir)=vpoten(ir)+value
              END IF
           END DO
        END DO
     END DO
  ELSE
     CALL errore('add_efield', 'wrong edir', 1)
  ENDIF
  first=.FALSE.
  RETURN
END SUBROUTINE add_efield

