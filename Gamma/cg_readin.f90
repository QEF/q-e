!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cg_readin
  !-----------------------------------------------------------------------
  !
  USE ions_base, ONLY : nat
  use pwcom
  use cgcom
  use io_files, only: tmp_dir, prefix
#ifdef __PARA
  use para, only: me
  use mp, only: mp_bcast
#endif
  implicit none
  integer :: iunit =5, ionode_id = 0
  character(len=256) :: outdir
  namelist /inputph/ prefix, fildyn, trans, epsil, raman, nmodes,     &
            tr2_ph, niter_ph, amass, outdir, asr, deltatau, nderiv, &
            first, last
                                                                                
  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr
  INTEGER, EXTERNAL   :: iargc

  !
  call start_clock('cg_readin')
  !
  outdir = './'
  prefix = 'pwscf'
  fildyn = 'matdyn'
  epsil  = .true.
  trans  = .true.
  raman  = .false.
  asr    = .false.
  tr2_ph = 1.0e-12
  niter_ph= 50
  nmodes =  0
  deltatau= 0.0
  nderiv = 2
  first  = 1
  last   = 0
#ifdef __PARA
  if (me == 1) then
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

  read(iunit,'(a)') title_ph
  read(iunit,inputph)
  !
  tmp_dir = trim(outdir)
#ifdef __PARA
  endif
  !
  call mp_bcast(prefix,ionode_id)
  call mp_bcast(fildyn,ionode_id)
  call mp_bcast(trans,ionode_id)
  call mp_bcast(epsil,ionode_id)
  call mp_bcast(raman,ionode_id)
  call mp_bcast(nmodes,ionode_id)
  call mp_bcast(tr2_ph,ionode_id)
  call mp_bcast(niter_ph,ionode_id)
  call mp_bcast(amass,ionode_id)
  call mp_bcast(tr2_ph,ionode_id)
  call mp_bcast(tmp_dir,ionode_id)
  call mp_bcast(asr,ionode_id)
  call mp_bcast(deltatau,ionode_id)
  call mp_bcast(nderiv,ionode_id)
  call mp_bcast(first,ionode_id)
  call mp_bcast(last,ionode_id)
  !
#endif
  !
  !  read the input file produced by 'punch' subroutine in pwscf program
  !  allocate memory and recalculate what is needed
  !
  call read_file
  !
  !  various checks
  !
  if (.not.trans .and. .not.epsil)                                  &
       &     call errore('data','nothing to do',1)
  if (nks.ne.1) call errore('data','too many k-points',1)
  !      if (xk(1,1).ne.0.0 .or. xk(2,1).ne.0.0 .or. xk(3,1).ne.0.0)
  !     &    call errore('data','only k=0 allowed',1)
  if (nmodes.gt.3*nat .or. nmodes.lt.0)                             &
       &     call errore('data','wrong number of normal modes',1)
  if (epsil .and. nmodes.ne.0) call errore('data','not allowed',1)
  if (raman .and. deltatau.le.0.d0)                                 &
       &     call errore('data','deltatau > 0 needed for raman CS',1)
  if (nderiv.ne.2 .and. nderiv.ne.4) &
       call errore('data','nderiv not allowed',1)
  !
  if (last.eq.0) last=3*nat
  !
  call cg_readmodes(iunit)
  !
  call stop_clock('cg_readin')
  !
  return
end subroutine cg_readin
!
!-----------------------------------------------------------------------
subroutine cg_readmodes(iunit)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE ions_base, ONLY : nat
  USE kinds, only: DP
  use pwcom
  use cgcom
#ifdef __PARA
  use para, only: me
  use mp, only: mp_bcast
#endif
  !
  implicit none
  integer :: iunit
  !
  integer :: ionode_id = 0, na, nu, mu
  real(kind=DP) utest, unorm, DDOT
  !
  ! allocate space for modes, dynamical matrix, auxiliary stuff
  !
  allocate  (u(  3*nat, 3*nat))    
  allocate  (dyn(3*nat, 3*nat))    
  allocate  (equiv_atoms( nat, nat))    
  allocate  (n_equiv_atoms( nat))    
  allocate  (has_equivalent(nat))    
  !
  ! nmodes not given: use defaults (all modes) as normal modes ...
  !
  if (nmodes.eq.0) then
     call find_equiv_sites (nat,nat,nsym,irt,has_equivalent,        &
          &      n_diff_sites,n_equiv_atoms,equiv_atoms)
     if (n_diff_sites .le. 0 .or. n_diff_sites .gt. nat)            &
          &      call errore('equiv.sites','boh!',1)
     !
     ! these are all modes, but only independent modes are calculated
     !
     nmodes = 3*nat
     u(:,:) = 0.d0
     do nu = 1,nmodes
        u(nu,nu) = 1.0
     end do
     ! look if ASR can be exploited to reduce the number of calculations
     ! we need to locate an independent atom with no equivalent atoms
     nasr=0
     if (asr.and.n_diff_sites.gt.1) then
        do na = 1, n_diff_sites
           if (n_equiv_atoms(na).eq.1 ) then
              nasr = equiv_atoms(na, 1)
              go to 1
           end if
        end do
 1      continue
     end if
  else
     if (asr) call errore('readin','warning: asr disabled',-1)
     nasr=0
     !
     ! ... otherwise read normal modes from input
     !
     do na = 1,nat
        has_equivalent(na) = 0
     end do
#ifdef __PARA
     if (me == 0) then
#endif
     do nu = 1,nmodes
        read (iunit,*,end=10,err=10) (u(mu,nu), mu=1,3*nat)
     end do
#ifdef __PARA
     endif
     call mp_bcast(u,ionode_id)
#endif
     do nu = 1,nmodes
        do mu = 1, nu-1
           utest = DDOT(3*nat,u(1,nu),1,u(1,mu),1)
           if (abs(utest).gt.1.0e-10) then
              print *, ' warning: input modes are not orthogonal'
              call DAXPY(3*nat,-utest,u(1,mu),1,u(1,nu),1)
           end if
        end do
        unorm = sqrt(DDOT(3*nat,u(1,nu),1,u(1,nu),1))
        if (abs(unorm).lt.1.0e-10) go to 10
        call DSCAL(3*nat,1.0/unorm,u(1,nu),1)
     end do
     go to 20
10   call errore('phonon','wrong data read',1)
  endif
20 continue
  !
  return
end subroutine cg_readmodes
