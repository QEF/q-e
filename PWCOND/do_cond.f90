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
  USE noncollin_module, ONLY : noncolin
  USE io_global, ONLY : stdout, ionode, ionode_id

  USE mp

  IMPLICIT NONE

  CHARACTER(len=3) nodenumber
  REAL(DP) :: wtot
  INTEGER :: ik, ipol, ien, ios 

  NAMELIST /inputcond/ outdir, prefixt, prefixl, prefixs, prefixr,     &
                       band_file, tran_file, save_file, fil_loc,       &
                       lwrite_loc, lread_loc, lwrite_cond, lread_cond, & 
                       orbj_in,orbj_fin,ikind,iofspin,llocal,          & 
                       bdl, bds, bdr, nz1, energy0, denergy, ecut2d,   &
                       ewind, epsproj, delgep, cutplot, lorb
                                                                               
  nd_nmbr=nodenumber
  CALL init_clocks(.TRUE.)
  CALL start_clock('PWCOND')
  CALL start_clock('init')
!
!   set default values for variables in namelist
!                                             
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefixt = ' '
  prefixl = ' '
  prefixs = ' '
  prefixr = ' '
  band_file = ' '
  tran_file = ' '
  save_file = ' ' 
  fil_loc = ' '
  lwrite_loc = .FALSE.
  lread_loc = .FALSE.
  lwrite_cond = .FALSE.
  lread_cond  = .FALSE.
  orbj_in = 0
  orbj_fin = 0
  iofspin = 1
  ikind = 0
  bdl = 0.d0
  bds = 0.d0
  bdr = 0.d0
  nz1 = 11
  energy0 = 0.d0
  denergy = 0.d0
  ecut2d = 0.d0
  ewind = 1.d0
  llocal = .FALSE.
  epsproj = 1.d-3
  delgep = 5.d-10
  cutplot = 2.d0
  lorb=.FALSE.

  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     !     reading the namelist inputcond
     !
     READ (5, inputcond, err=200, iostat=ios )
200  CALL errore ('do_cond','reading inputcond namelist',ABS(ios))
     tmp_dir=trimcheck (outdir)
     !
     !     Reading 2D k-point
     READ(5, *, err=250, iostat=ios ) nkpts
250     CALL errore ('do_cond','reading number of kpoints',ABS(ios))
     IF (nkpts>0) THEN
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
300     CALL errore ('do_cond','2D k-point',ABS(ios))
     ELSE
        ALLOCATE( xyk(2,npk) )
        ALLOCATE( wkpt(npk) )
        READ(5, *, err=350, iostat=ios) nk1ts, nk2ts, k1ts, k2ts
350     CALL errore ('do_cond','2D k-point size or shift',ABS(ios))
     ENDIF

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
  CALL mp_bcast( prefixt, ionode_id )
  CALL mp_bcast( prefixl, ionode_id )
  CALL mp_bcast( prefixs, ionode_id )
  CALL mp_bcast( prefixr, ionode_id )
  CALL mp_bcast( band_file, ionode_id )
  CALL mp_bcast( tran_file, ionode_id )
  CALL mp_bcast( fil_loc, ionode_id )
  CALL mp_bcast( save_file, ionode_id )
  CALL mp_bcast( lwrite_loc, ionode_id )
  CALL mp_bcast( lread_loc, ionode_id )
  CALL mp_bcast( lwrite_cond, ionode_id )
  CALL mp_bcast( lread_cond, ionode_id )
  CALL mp_bcast( ikind, ionode_id )
  CALL mp_bcast( iofspin, ionode_id )
  CALL mp_bcast( orbj_in, ionode_id )
  CALL mp_bcast( orbj_fin, ionode_id )
  CALL mp_bcast( llocal, ionode_id )
  CALL mp_bcast( bdl, ionode_id )
  CALL mp_bcast( bds, ionode_id )
  CALL mp_bcast( bdr, ionode_id )
  CALL mp_bcast( nz1, ionode_id )
  CALL mp_bcast( energy0, ionode_id )
  CALL mp_bcast( denergy, ionode_id )
  CALL mp_bcast( ecut2d, ionode_id )
  CALL mp_bcast( ewind, ionode_id )
  CALL mp_bcast( epsproj, ionode_id )
  CALL mp_bcast( delgep, ionode_id )
  CALL mp_bcast( cutplot, ionode_id )
  CALL mp_bcast( nkpts, ionode_id )
  CALL mp_bcast( nenergy, ionode_id )
  CALL mp_bcast( nk1ts, ionode_id )
  CALL mp_bcast( nk2ts, ionode_id )
  CALL mp_bcast( k1ts, ionode_id )
  CALL mp_bcast( k2ts, ionode_id )

  IF ( .NOT. ionode ) THEN
     IF (nkpts>0) THEN
        ALLOCATE( xyk(2,nkpts) )
        ALLOCATE( wkpt(nkpts) )
     ELSE
        ALLOCATE( xyk(2,npk) )
        ALLOCATE( wkpt(npk) )
     ENDIF
     ALLOCATE( earr(nenergy) )
     ALLOCATE( tran_tot(nenergy) )
  ENDIF
  IF (nkpts>0) THEN
     CALL mp_bcast( xyk, ionode_id )
     CALL mp_bcast( wkpt, ionode_id )
  ENDIF
  CALL mp_bcast( earr, ionode_id )
  CALL mp_bcast( tran_tot, ionode_id )


!
! Now allocate space for pwscf variables, read and check them.
!
  
IF (lread_cond) THEN
  call save_cond (.false.,1,efl,nrzl,nocrosl,noinsl,   &
                  norbl,rl,rabl,betarl)
  if(ikind.eq.1) then
    call save_cond (.false.,2,efs,nrzs,-1,      &
                             noinss,norbs,rs,rabs,betars)
    norbr = norbl
    nocrosr = nocrosl
    noinsr = noinsl
  endif
  if(ikind.eq.2) then
    call save_cond (.false.,2,efs,nrzs,-1,      &
                             noinss,norbs,rs,rabs,betars)

    call save_cond (.false.,3,efr,nrzr,nocrosr,&
                             noinsr,norbr,rr,rabr,betarr)
  endif
ELSE
  IF (prefixt.ne.' ') then
    prefix = prefixt
    call read_file
    call init_us_1
    call newd
    IF (ikind.eq.0) then
      CALL init_cond(1,'t')
    ELSEIF (ikind.eq.1) then
      CALL init_cond(2,'t')
    ELSEIF (ikind.eq.2) then
      CALL init_cond(3,'t')
    ENDIF
    CALL clean_pw(.true.)
  ENDIF
  IF (prefixl.ne.' ') then
    prefix = prefixl
    call read_file
    call init_us_1
    call newd
    CALL init_cond(1,'l')
    CALL clean_pw(.true.)
  ENDIF
  IF (prefixs.ne.' ') then
    prefix = prefixs
    call read_file
    call init_us_1
    call newd
    CALL init_cond(1,'s')
    CALL clean_pw(.true.)
  ENDIF
  IF (prefixr.ne.' ') then
    prefix = prefixr
    call read_file
    call init_us_1
    call newd
    CALL init_cond(1,'r')
    CALL clean_pw(.true.)
  ENDIF
ENDIF

IF (lorb.and.okvan) call errore('do_cond','lorb not working with US-PP',1)
#ifdef __PARA
   IF (lorb) call errore('do_cond','lorb not working in parallel',1)
#endif

IF (nkpts==0) THEN
   IF (ionode) THEN
      CALL kpoint_grid( nsym, s, bg, npk, k1ts, k2ts, 0, &
                        nk1ts, nk2ts, 1, nkpts, xk, wkpt )
      DO ik=1,nkpts
         xyk(1,ik)=xk(1,ik)
         xyk(2,ik)=xk(2,ik)
      ENDDO
   ENDIF
   CALL mp_bcast( xyk, ionode_id )
   CALL mp_bcast( wkpt, ionode_id )
ENDIF


IF (lwrite_cond) then
  call save_cond (.true.,1,efl,nrzl,nocrosl,noinsl,         &
                  norbl,rl,rabl,betarl)
  if(ikind.gt.0) call save_cond (.true.,2,efs,nrzs,-1,      &
                             noinss,norbs,rs,rabs,betars)
  if(ikind.gt.1) call save_cond (.true.,3,efr,nrzr,nocrosr,&
                             noinsr,norbr,rr,rabr,betarr)
endif

  CALL cond_out

  CALL stop_clock('init')

  IF (llocal) &
  CALL local_set(nocrosl,noinsl,norbl,noinss,norbs,nocrosr,noinsr,norbr)

  do ik=1, nkpts

    WRITE( stdout, '(/,8x,"k(",i5,") = (",2f12.7,"), wk =",f12.7,/)') ik, &
         (xyk (ipol, ik) , ipol = 1, 2) , wkpt (ik)

    CALL init_gper(ik)

    CALL local 

    DO ien=1, nenergy
      eryd = earr(ien)/rytoev + efl
      CALL form_zk(n2d, nrzpl, zkrl, zkl, eryd, tpiba)
      CALL scatter_forw(nrzl, nrzpl, zl, psiperl, zkl, norbl,     &
                        tblml, crosl, taunewl, rl, rabl, betarl) 
      CALL compbs(1, nocrosl, norbl, nchanl, kvall, kfunl, kfundl, &
                  kintl, kcoefl) 

      IF (ikind.EQ.2) THEN
        eryd = earr(ien)/rytoev + efr
        CALL form_zk(n2d, nrzpr, zkrr, zkr, eryd, tpiba)
        CALL scatter_forw(nrzr, nrzpr, zr, psiperr, zkr, norbr,    &
                          tblmr, crosr, taunewr, rr, rabr, betarr)
        CALL compbs(0, nocrosr, norbr, nchanr, kvalr, kfunr, kfundr,&
                     kintr, kcoefr) 
      ENDIF

      CALL summary_band(ik,ien)

      IF (ikind.NE.0) THEN
         eryd = earr(ien)/rytoev + efs
         CALL form_zk(n2d, nrzps, zkrs, zks, eryd, tpiba)
         CALL scatter_forw(nrzs, nrzps, zs, psipers, zks, norbs,    &
                          tblms, cross, taunews, rs, rabs, betars)

         WRITE(stdout,*) 'to transmit'

         CALL transmit(ik, ien)
      ENDIF


   ENDDO
   CALL free_mem
  ENDDO

  IF(ikind.GT.0.AND.tran_file.NE.' ') CALL summary_tran()

  CALL print_clock_pwcond()
  CALL stop_clock('PWCOND')

  RETURN

END SUBROUTINE do_cond


