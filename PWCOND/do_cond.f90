!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
SUBROUTINE do_cond(nodenumber)
!  
!   This is the main driver of the pwcond.x program.
!   It calculates the complex band structure, solves the
!   scattering problem and calculates the transmission coefficient. 
!
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau
  USE pwcom
  USE cond 
  USE io_files 
  USE io_global, ONLY : ionode, ionode_id

  USE mp

  IMPLICIT NONE

  CHARACTER(len=3) nodenumber
  REAL(kind=DP), PARAMETER :: eps=1.d-8
  REAL(kind=DP) :: wtot
  INTEGER :: ik, ien, ios, orbin, orbfin 
  LOGICAL :: write0

  NAMELIST /inputcond/ outdir, prefix, band_file, tran_file, fil_loc,  & 
                       lwrite_loc, lread_loc, ikind, iofspin, llocal,  & 
                       bdl1, bdl2, bdr1, bdr2, nz1, dnslab, energy0,   &
                       denergy, ecut2d, ewind, epsproj, delgep,cutplot,&
                       llapack
                                                                               
  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr, ILEN
  INTEGER, EXTERNAL   :: iargc

  nd_nmbr=nodenumber
  CALL init_clocks(.TRUE.)
  CALL start_clock('PWCOND')
  CALL start_clock('init')
!
!   set default values for variables in namelist
!                                             
  outdir = './'
  prefix = ' '
  band_file = ' '
  tran_file = ' '
  fil_loc = ' '
  lwrite_loc = .FALSE.
  lread_loc = .FALSE.
  iofspin = 1
  ikind = 0
  bdl1 = 0.d0
  bdl2 = 0.d0
  bdr1 = 0.d0
  bdr2 = 0.d0
  nz1 = 11
  dnslab = 0
  energy0 = 0.d0
  denergy = 0.d0
  ecut2d = 0.d0
  ewind = 0.d0
  llocal = .FALSE.
  llapack = .FALSE.
  epsproj = 1.d-3
  delgep = 0.d0
  cutplot = 2.d0

  IF ( ionode )  THEN
     !
     ! ... Input from file ?
     !
     nargs = iargc()
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, input_file )
        IF ( TRIM( input_file ) == '-input' .OR. &
             TRIM( input_file ) == '-inp'   .OR. &
             TRIM( input_file ) == '-in' ) THEN
           !
           CALL getarg( ( iiarg + 1 ) , input_file )
           OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
                STATUS = 'OLD', IOSTAT = ierr )
           CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
           !
        END IF
        !
     END DO
     !
     !     reading the namelist inputpp
     !
     READ (5, inputcond, err=200, iostat=ios )
200  CALL errore ('do_cond','reading inputcond namelist',ABS(ios))
     tmp_dir=TRIM(outdir)
     !
     !     Reading 2D k-point
     READ(5, *, err=300, iostat=ios ) nkpts
     ALLOCATE( xyk(2,nkpts) )
     ALLOCATE( wkpt(nkpts) )
     wtot = 0.d0
     DO ik = 1, nkpts
        READ(5, *, err=300, iostat=ios) xyk(1,ik), xyk(2,ik), wkpt(ik)
        wtot = wtot + wkpt(ik)
     ENDDO
     DO ik = 1, nkpts
        wkpt(ik) = wkpt(ik)/wtot
     ENDDO
300  CALL errore ('do_cond','2D k-point',ABS(ios))

     !
     !     To form the array of energies for calculation
     !
     READ(5, *, err=400, iostat=ios ) nenergy
     ALLOCATE( earr(nenergy) )
     ALLOCATE( tran_tot(nenergy) )
     IF(ABS(denergy).LE.1.d-8) THEN
        !   the list of energies is read
        DO ien = 1, nenergy
           READ(5, *, err=400, iostat=ios) earr(ien)
        ENDDO
     ELSE
        !   the array of energies is automatically formed
        DO ien = 1, nenergy
           earr(ien) = energy0 + (ien-1)*denergy
           tran_tot(ien) = 0.d0 
        ENDDO
     ENDIF
400  CALL errore ('do_cond','reading energy list',ABS(ios))
     !
  END IF

!
! ... Broadcast variables
!
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( band_file, ionode_id )
  CALL mp_bcast( tran_file, ionode_id )
  CALL mp_bcast( fil_loc, ionode_id )
  CALL mp_bcast( lwrite_loc, ionode_id )
  CALL mp_bcast( lread_loc, ionode_id )
  CALL mp_bcast( ikind, ionode_id )
  CALL mp_bcast( iofspin, ionode_id )
  CALL mp_bcast( llocal, ionode_id )
  CALL mp_bcast( bdl1, ionode_id )
  CALL mp_bcast( bdl2, ionode_id )
  CALL mp_bcast( bdr1, ionode_id )
  CALL mp_bcast( bdr2, ionode_id )
  CALL mp_bcast( nz1, ionode_id )
  CALL mp_bcast( dnslab, ionode_id )
  CALL mp_bcast( energy0, ionode_id )
  CALL mp_bcast( denergy, ionode_id )
  CALL mp_bcast( ecut2d, ionode_id )
  CALL mp_bcast( ewind, ionode_id )
  CALL mp_bcast( epsproj, ionode_id )
  CALL mp_bcast( delgep, ionode_id )
  CALL mp_bcast( cutplot, ionode_id )
  CALL mp_bcast( llapack, ionode_id )
  CALL mp_bcast( nkpts, ionode_id )
  CALL mp_bcast( nenergy, ionode_id )

  IF ( .NOT. ionode ) THEN
     ALLOCATE( xyk(2,nkpts) )
     ALLOCATE( wkpt(nkpts) )
     ALLOCATE( earr(nenergy) )
     ALLOCATE( tran_tot(nenergy) )
  ENDIF
  CALL mp_bcast( xyk, ionode_id )
  CALL mp_bcast( wkpt, ionode_id )
  CALL mp_bcast( earr, ionode_id )
  CALL mp_bcast( tran_tot, ionode_id )

!
! Now allocate space for pwscf variables, read and check them.
!
  
  CALL read_file
  CALL openfil
  CALL struc_fact (nat,tau,ntyp,ityp,ngm,g,bg,                &
                   nr1,nr2,nr3,strf,eigts1,eigts2,eigts3)
  CALL init_us_1
  CALL newd

!
! Allocation for pwcond variables
!
  CALL allocate_cond
  CALL init_cond
  CALL stop_clock('init')

  IF (llocal) &
    CALL local_set(nocrosl,noinsl,noinss,nocrosr,noinsr,norb,norbs)
  CALL poten

  DO ik=1, nkpts

    CALL init_gper(ik)
!
! The main loop
!
    eryd = earr(1)/rytoev + ef
    CALL local
    DO ien=1, nenergy
      eryd = earr(ien)/rytoev + ef
      CALL form_zk(n2d, nrzp, zkr, zk, eryd, tpiba)
      orbin=1
      orbfin=orbin-1+2*nocrosl+noinsl
      CALL scatter_forw(bdl1, bdl2, orbin, orbfin)
      orbin=1
      orbfin=orbin-1+2*nocrosl+noinsl   
      
      CALL compbs(0, bdl1, bdl2, nocrosl,   & 
              orbfin-orbin+1, orbin, nchanl, kvall,    &
              kfunl, kfundl, kintl, kcoefl)
      IF (ikind.EQ.2) THEN
       orbin=2*nocrosl+noinsl+noinss+1
       orbfin=orbin-1+2*nocrosr+noinsr
       CALL scatter_forw(bdr1, bdr2, orbin, orbfin)
       orbin=2*nocrosl+noinsl+noinss+1
       orbfin=orbin-1+2*nocrosr+noinsr    
       CALL compbs(1, bdr1, bdr2, nocrosr,  &
            orbfin-orbin+1, orbin, nchanr, kvalr,        &
            kfunr, kfundr, kintr, kcoefr)
      ENDIF
      CALL summary_band(ik,ien)
      IF (ikind.NE.0) THEN
       orbin=nocrosl+noinsl+1
       orbfin=orbin-1+nocrosl+noinss+nocrosr
       CALL scatter_forw(bdl2, bdr1, orbin, orbfin)
       CALL transmit(ik,ien) 
      ENDIF                                   
    ENDDO                    
   CALL free_mem

  ENDDO

  IF(ikind.GT.0.AND.tran_file.NE.' ')  &
   CALL summary_tran(tran_file,nenergy,earr,tran_tot)

  CALL print_clock_pwcond()
  CALL stop_clock('PWCOND')
  RETURN
END SUBROUTINE do_cond                                    

