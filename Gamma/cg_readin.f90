!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE cg_readin()
  !-----------------------------------------------------------------------
  !
  USE ions_base, ONLY : nat, amass
  USE pwcom
  USE cgcom
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE noncollin_module, ONLY : noncolin
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  INTEGER :: iunit =5
  CHARACTER(len=256) :: outdir
  NAMELIST /inputph/ prefix, fildyn, trans, epsil, raman, nmodes,     &
            tr2_ph, niter_ph, amass, outdir, asr, deltatau, nderiv, &
            first, last, recover
  !
  CALL start_clock('cg_readin')
  !
  outdir = './'
  prefix = 'pwscf'
  fildyn = 'matdyn'
  epsil  = .TRUE.
  trans  = .TRUE.
  raman  = .FALSE.
  asr    = .FALSE.
  tr2_ph = 1.0e-12
  niter_ph= 50
  nmodes =  0
  deltatau= 0.0
  nderiv = 2
  first  = 1
  last   = 0
  recover=.FALSE.
  !
  IF ( ionode ) THEN
     !
     CALL input_from_file ( )
     !
     READ(iunit,'(a)') title_ph
     READ(iunit,inputph)
     !
     tmp_dir = TRIM(outdir)
     !
  END IF
  !
  CALL mp_bcast(prefix,ionode_id)
  CALL mp_bcast(fildyn,ionode_id)
  CALL mp_bcast(trans,ionode_id)
  CALL mp_bcast(epsil,ionode_id)
  CALL mp_bcast(raman,ionode_id)
  CALL mp_bcast(nmodes,ionode_id)
  CALL mp_bcast(tr2_ph,ionode_id)
  CALL mp_bcast(niter_ph,ionode_id)
  CALL mp_bcast(amass,ionode_id)
  CALL mp_bcast(tr2_ph,ionode_id)
  CALL mp_bcast(tmp_dir,ionode_id)
  CALL mp_bcast(asr,ionode_id)
  CALL mp_bcast(deltatau,ionode_id)
  CALL mp_bcast(nderiv,ionode_id)
  CALL mp_bcast(first,ionode_id)
  CALL mp_bcast(last,ionode_id)
  CALL mp_bcast(recover,ionode_id)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  CALL read_file

  if (noncolin) call errore('cg_readin','noncolinear version not available',1)
  !
  !  various checks
  !
  IF (.NOT. gamma_only) CALL errore('cg_readin', &
      'need pw.x data file produced using Gamma tricks',1)
  IF (okvan) CALL errore('cg_readin', &
      'ultrasoft pseudopotential not implemented',1)
  IF (doublegrid) &
      CALL errore('cg_readin', 'double grid not implemented',1)
  IF (.NOT.trans .AND. .NOT.epsil)                                  &
       &     CALL errore('cg_readin','nothing to do',1)
  IF (nks.NE.1) CALL errore('cg_readin','too many k-points',1)
  !      if (xk(1,1).ne.0.0 .or. xk(2,1).ne.0.0 .or. xk(3,1).ne.0.0)
  !     &    call errore('data','only k=0 allowed',1)
  IF (nmodes.GT.3*nat .OR. nmodes.LT.0)                             &
       &     CALL errore('cg_readin','wrong number of normal modes',1)
  IF (epsil .AND. nmodes.NE.0) CALL errore('cg_readin','not allowed',1)
  !
  IF (raman .AND. deltatau.LE.0.d0)                                 &
       &     CALL errore('cg_readin','deltatau > 0 needed for raman CS',1)
  IF (nderiv.NE.2 .AND. nderiv.NE.4) &
       CALL errore('cg_readin','nderiv not allowed',1)
  !
  IF (last.EQ.0) last=3*nat
  !
  CALL cg_readmodes(iunit)
  !
  CALL stop_clock('cg_readin')
  !
  RETURN
END SUBROUTINE cg_readin
!
!-----------------------------------------------------------------------
SUBROUTINE cg_readmodes(iunit)
  !-----------------------------------------------------------------------
  !
  USE ions_base, ONLY : nat
  USE kinds,     ONLY : DP
  USE pwcom
  USE cgcom
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: iunit
  !
  INTEGER :: na, nu, mu
  REAL(DP) utest, unorm, DDOT
  !
  ! allocate space for modes, dynamical matrix, auxiliary stuff
  !
  ALLOCATE  (u(  3*nat, 3*nat))    
  ALLOCATE  (dyn(3*nat, 3*nat))    
  ALLOCATE  (equiv_atoms( nat, nat))    
  ALLOCATE  (n_equiv_atoms( nat))    
  ALLOCATE  (has_equivalent(nat))    
  !
  ! nmodes not given: use defaults (all modes) as normal modes ...
  !
  IF (nmodes.EQ.0) THEN
     CALL find_equiv_sites (nat,nat,nsym,irt,has_equivalent,        &
          &      n_diff_sites,n_equiv_atoms,equiv_atoms)
     IF (n_diff_sites .LE. 0 .OR. n_diff_sites .GT. nat)            &
          &      CALL errore('equiv.sites','boh!',1)
     !
     ! these are all modes, but only independent modes are calculated
     !
     nmodes = 3*nat
     u(:,:) = 0.d0
     DO nu = 1,nmodes
        u(nu,nu) = 1.0
     END DO
     ! look if ASR can be exploited to reduce the number of calculations
     ! we need to locate an independent atom with no equivalent atoms
     nasr=0
     IF (asr.AND.n_diff_sites.GT.1) THEN
        DO na = 1, n_diff_sites
           IF (n_equiv_atoms(na).EQ.1 ) THEN
              nasr = equiv_atoms(na, 1)
              go to 1
           END IF
        END DO
 1      CONTINUE
     END IF
  ELSE
     IF (asr) CALL infomsg ('readin','warning: asr disabled', -1)
     nasr=0
     !
     ! ... otherwise read normal modes from input
     !
     DO na = 1,nat
        has_equivalent(na) = 0
     END DO

     IF ( ionode ) THEN
        !
        DO nu = 1,nmodes
           READ (iunit,*,END=10,err=10) (u(mu,nu), mu=1,3*nat)
        END DO
        !
     END IF
     CALL mp_bcast(u,ionode_id)
     DO nu = 1,nmodes
        DO mu = 1, nu-1
           utest = DDOT(3*nat,u(1,nu),1,u(1,mu),1)
           IF (ABS(utest).GT.1.0e-10) THEN
              PRINT *, ' warning: input modes are not orthogonal'
              CALL DAXPY(3*nat,-utest,u(1,mu),1,u(1,nu),1)
           END IF
        END DO
        unorm = SQRT(DDOT(3*nat,u(1,nu),1,u(1,nu),1))
        IF (ABS(unorm).LT.1.0e-10) go to 10
        CALL DSCAL(3*nat,1.0/unorm,u(1,nu),1)
     END DO
     go to 20
10   CALL errore('phonon','wrong data read',1)
  ENDIF
20 CONTINUE
  !
  RETURN
END SUBROUTINE cg_readmodes
