!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine do_cond(nodenumber)
!  
!   This is the main driver of the pwcond.x program.
!   It calculates the complex band structure, solves the
!   scattering problem and calculates the transmission coefficient. 
!
#include "f_defs.h"
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau
  use pwcom
  use cond 
  use io_files 
  use io_global, only : ionode_id
#ifdef __PARA
  use para, only: me
  use mp
#endif
implicit none
  character(len=3) nodenumber
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: wtot
  integer :: ik, ien, ios, orbin, orbfin 
  logical :: write0

  namelist /inputcond/ outdir, prefix, band_file, tran_file, fil_loc,  & 
                       lwrite_loc, lread_loc, ikind, iofspin, llocal,  & 
                       bdl1, bdl2, bdr1, bdr2, nz1, dnslab, energy0,   &
                       denergy, ecut2d, ewind, epsproj, delgep,cutplot,&
                       llapack
                                                                                
  CHARACTER (LEN=256) :: input_file
  INTEGER             :: nargs, iiarg, ierr, ilen
  INTEGER, EXTERNAL   :: iargc

  nd_nmbr=nodenumber
!
!   set default values for variables in namelist
!                                             
   outdir = './'
   prefix = ' '
   band_file = ' '
   tran_file = ' '
   fil_loc = ' '
   lwrite_loc = .false.
   lread_loc = .false.
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
   llocal = .false.
   llapack = .false.
   epsproj = 1.d-3
   delgep = 0.d0
   cutplot = 2.d0

   call start_clock('cond')


#ifdef __PARA
  if (me == 1)  then
#endif
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
   read (5, inputcond, err=200, iostat=ios )
200   call errore ('do_cond','reading inputcond namelist',abs(ios))
   tmp_dir=trim(outdir)
!
!     Reading 2D k-point
  read(5, *, err=300, iostat=ios ) nkpts
  allocate( xyk(2,nkpts) )
  allocate( wkpt(nkpts) )
  wtot = 0.d0
  do ik = 1, nkpts
   read(5, *, err=300, iostat=ios) xyk(1,ik), xyk(2,ik), wkpt(ik)
   wtot = wtot + wkpt(ik)
  enddo
  do ik = 1, nkpts
   wkpt(ik) = wkpt(ik)/wtot
  enddo
300 call errore ('do_cond','2D k-point',abs(ios))

!
!     To form the array of energies for calculation
!
  read(5, *, err=400, iostat=ios ) nenergy
  allocate( earr(nenergy) )
  allocate( tran_tot(nenergy) )
  if(abs(denergy).le.1.d-8) then
!   the list of energies is read
    do ien = 1, nenergy
      read(5, *, err=400, iostat=ios) earr(ien)
    enddo
  else
!   the array of energies is automatically formed
    do ien = 1, nenergy
      earr(ien) = energy0 + (ien-1)*denergy
      tran_tot(ien) = 0.d0 
    enddo
  endif
400 call errore ('do_cond','reading energy list',abs(ios))
#ifdef __PARA
  end if

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
  if (me .ne. 1) then
     allocate( xyk(2,nkpts) )
     allocate( wkpt(nkpts) )
     allocate( earr(nenergy) )
     allocate( tran_tot(nenergy) )
  endif
  CALL mp_bcast( xyk, ionode_id )
  CALL mp_bcast( wkpt, ionode_id )
  CALL mp_bcast( earr, ionode_id )
  CALL mp_bcast( tran_tot, ionode_id )
#endif


!
! Now allocate space for pwscf variables, read and check them.
!
  
  call read_file
  call openfil
  call struc_fact (nat,tau,ntyp,ityp,ngm,g,bg,                &
                   nr1,nr2,nr3,strf,eigts1,eigts2,eigts3)
  call init_us_1
  call newd

!
! Allocation for pwcond variables
!
  call allocate_cond
  call init_cond

  if (llocal) &
    call local_set(nocrosl,noinsl,noinss,nocrosr,noinsr,norb,norbs)
  call poten

  do ik=1, nkpts

    call init_gper(ik)
!
! The main loop
!
    eryd = earr(1)/rytoev + ef
    call local
    do ien=1, nenergy
      eryd = earr(ien)/rytoev + ef
      call form_zk(n2d, nrzp, zkr, zk, eryd, tpiba)
      orbin=1
      orbfin=orbin-1+2*nocrosl+noinsl
      call scatter_forw(bdl1, bdl2, orbin, orbfin)
      orbin=1
      orbfin=orbin-1+2*nocrosl+noinsl   
      call compbs(0, bdl1, bdl2, nocrosl,   & 
              orbfin-orbin+1, orbin, nchanl, kvall,    &
              kfunl, kfundl, kintl, kcoefl)
      if (ikind.eq.2) then
       orbin=2*nocrosl+noinsl+noinss+1
       orbfin=orbin-1+2*nocrosr+noinsr
       call scatter_forw(bdr1, bdr2, orbin, orbfin)
       orbin=2*nocrosl+noinsl+noinss+1
       orbfin=orbin-1+2*nocrosr+noinsr    
       call compbs(1, bdr1, bdr2, nocrosr,  &
            orbfin-orbin+1, orbin, nchanr, kvalr,        &
            kfunr, kfundr, kintr, kcoefr)
      endif
      call summary_band(ik,ien)
      if (ikind.ne.0) then
       orbin=nocrosl+noinsl+1
       orbfin=orbin-1+nocrosl+noinss+nocrosr
       call scatter_forw(bdl2, bdr1, orbin, orbfin)
       call transmit(ik,ien) 
      endif                                   
    enddo                    
   call free_mem

  enddo

  if(ikind.gt.0)  &
   call summary_tran(tran_file,nenergy,earr,tran_tot)

  call print_clock('cond')          
  call stop_clock('cond')
  return
end subroutine do_cond                                    

