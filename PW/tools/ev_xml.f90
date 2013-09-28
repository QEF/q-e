!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ev_xml
!
!   This module contains routines to write the information obtained by the
!   ev.x program in an xml file.
!
USE iotk_module
!
USE kinds,     ONLY : DP

IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: write_evdata_xml

  INTEGER :: iunout
  !

  CONTAINS
!-----------------------------------------------------------------------
    SUBROUTINE write_evdata_xml &
        (npt,fac,v0,etot,efit,istat,par,npar,emin,chisq,filout, ierr)
!-----------------------------------------------------------------------
!
  USE constants, ONLY : ry_kbar, bohr_radius_angs
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npt, istat, npar
  REAL(DP), INTENT(in):: v0(npt), etot(npt), efit(npt), emin, chisq, fac
  REAL(DP), INTENT(in):: par(npar)
  CHARACTER(len=256), INTENT(IN) :: filout
  INTEGER, INTENT(out) :: ierr
  !
  REAL(DP) :: p(npt), volume(2), a0(2), alldata(6,npt)
  INTEGER :: i
  CHARACTER(len=256) :: filename
  REAL(DP), EXTERNAL :: birch, keane

  IF (filout/=' ') THEN
     filename = TRIM(filout) // '.xml'
  ELSE
     filename = 'ev.xml'
  ENDIF

  ierr=0
  CALL iotk_free_unit( iunout, ierr )
  IF ( ierr /= 0 ) THEN
     ierr = 11 
     RETURN
  END IF
  !
  ! ... open XML descriptor
  !
  CALL iotk_open_write( iunout, FILE = TRIM( filename ), &
                          BINARY = .FALSE., IERR = ierr )
  IF ( ierr /= 0 ) THEN
     ierr = 12 
     RETURN
  END IF

  CALL iotk_write_begin(iunout, "EQUATIONS_OF_STATE" )
  IF (istat==1) THEN
     CALL iotk_write_dat(iunout, "EQUATION_TYPE", "Birch 1st order")
  ELSEIF (istat==2) THEN
     CALL iotk_write_dat(iunout, "EQUATION_TYPE", "Birch 2nd order")
  ELSEIF (istat==3) THEN
     CALL iotk_write_dat(iunout, "EQUATION_TYPE", "Keane")
  ELSEIF (istat==4) THEN
     CALL iotk_write_dat(iunout, "EQUATION_TYPE", "Murnaghan")
  ENDIF
  CALL iotk_write_dat(iunout, "CHI_SQUARE", chisq)
  CALL iotk_write_end(iunout, "EQUATIONS_OF_STATE" )

  IF (istat==1 .or. istat==2) THEN
     DO i=1,npt
        p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ELSE
     DO i=1,npt
        p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ENDIF

  DO i=1,npt
     alldata (1,i) = v0(i)
     alldata (2,i) = etot(i) 
     alldata (3,i) = efit(i)
     alldata (4,i) = etot(i) - efit(i)
     alldata (5,i) = p(i) 
     alldata (6,i) = etot(i) + p(i) * v0(i) / ry_kbar
  ENDDO

  CALL iotk_write_begin(iunout, "EQUATIONS_PARAMETERS" )

  volume(1)=par(1)
  volume(2)=par(1)*bohr_radius_angs**3
  CALL iotk_write_dat(iunout, "EQUILIBRIUM_VOLUME_AU_A", volume(:), COLUMNS=2 )
  CALL iotk_write_dat(iunout, "BULK_MODULUS_KBAR", par(2))
  CALL iotk_write_dat(iunout, "DERIVATIVE_BULK_MODULUS", par(3))
  CALL iotk_write_dat(iunout, "SECOND_DERIVATIVE_BULK_MODULUS", par(4))
  CALL iotk_write_dat(iunout, "MINIMUM_ENERGY_RY", emin)
  CALL iotk_write_dat(iunout, "CELL_FACTOR", fac)
  IF (fac /= 0.0_DP) THEN
     a0(1) = (par(1)/fac)**(1d0/3d0)
     a0(2) = (par(1)/fac)**(1d0/3d0) * bohr_radius_angs
     CALL iotk_write_dat(iunout, "CELL_PARAMETER_AU_A", a0, COLUMNS=2 )
  ENDIF
  CALL iotk_write_end(iunout, "EQUATIONS_PARAMETERS" )

  CALL iotk_write_begin(iunout, "FIT_CHECK" )
  CALL iotk_write_dat(iunout, "NUMBER_OF_DATA", npt )
  CALL iotk_write_dat(iunout, "VOL_ENE_EFIT_DELTA_P_GIBBS", &
                alldata(:,:), COLUMNS=6 )

  CALL iotk_write_end(iunout, "FIT_CHECK" )
  CALL iotk_close_write( iunout )

  RETURN
END SUBROUTINE write_evdata_xml
!
END MODULE ev_xml

