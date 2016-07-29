!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE cg_readin()
  !-----------------------------------------------------------------------
  !
  USE ions_base, ONLY: nat, amass
  USE fft_base,  ONLY: dffts
  USE gvecs,     ONLY: doublegrid
  USE klist,     ONLY: nks
  USE control_flags, ONLY: gamma_only, llondon
  USE uspp,      ONLY: okvan
  USE io_files,  ONLY: tmp_dir, prefix
  USE io_global, ONLY: ionode, ionode_id
  USE noncollin_module, ONLY: noncolin
  USE mp_bands,  ONLY: nbgrp, ntask_groups
  USE mp,        ONLY: mp_bcast
  USE mp_world,  ONLY: world_comm
  USE cgcom
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: iunit =5
  CHARACTER(len=256) :: outdir
  NAMELIST /inputph/ prefix, fildyn, trans, epsil, raman, nmodes,     &
            tr2_ph, niter_ph, amass, outdir, asr, deltatau, nderiv, &
            first, last, recover
  !
  CALL start_clock('cg_readin')
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  prefix = 'pwscf'
  fildyn = 'matdyn'
  epsil  = .true.
  trans  = .true.
  raman  = .false.
  asr    = .false.
  tr2_ph = 1.0d-12
  niter_ph= 50
  nmodes =  0
  deltatau= 0.0d0
  nderiv = 2
  first  = 1
  last   = 0
  recover=.false.
  !
  IF ( ionode ) THEN
     !
     CALL input_from_file ( )
     !
     READ(iunit,'(a)') title_ph
     READ(iunit,inputph)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast(prefix,ionode_id,world_comm)
  CALL mp_bcast(fildyn,ionode_id,world_comm)
  CALL mp_bcast(trans,ionode_id,world_comm)
  CALL mp_bcast(epsil,ionode_id,world_comm)
  CALL mp_bcast(raman,ionode_id,world_comm)
  CALL mp_bcast(nmodes,ionode_id,world_comm)
  CALL mp_bcast(tr2_ph,ionode_id,world_comm)
  CALL mp_bcast(niter_ph,ionode_id,world_comm)
  CALL mp_bcast(amass,ionode_id,world_comm)
  CALL mp_bcast(tr2_ph,ionode_id,world_comm)
  CALL mp_bcast(tmp_dir,ionode_id,world_comm)
  CALL mp_bcast(asr,ionode_id,world_comm)
  CALL mp_bcast(deltatau,ionode_id,world_comm)
  CALL mp_bcast(nderiv,ionode_id,world_comm)
  CALL mp_bcast(first,ionode_id,world_comm)
  CALL mp_bcast(last,ionode_id,world_comm)
  CALL mp_bcast(recover,ionode_id,world_comm)
  !
  !  read the input file produced by the pwscf program
  !  allocate memory and recalculate what is needed
  !
  CALL read_file

  IF (noncolin) CALL errore('cg_readin','noncolinear version not available',1)
  !
  !  various checks
  !
  IF (.not. gamma_only) CALL errore('cg_readin', &
      'need pw.x data file produced using Gamma tricks',1)
  IF ( llondon ) CALL errore('cg_readin', &
      'phonons with DFT-D not implemented',1)
  !
  !   band group not available
  !
  IF (nbgrp /=1 ) &
     CALL errore('cg_readin','band parallelization not available',1)

  IF (okvan) CALL errore('cg_readin', &
      'ultrasoft pseudopotential not implemented',1)

  IF (doublegrid) &
      CALL errore('cg_readin', 'double grid not implemented',1)
  IF (.not.trans .and. .not.epsil)                                  &
       &     CALL errore('cg_readin','nothing to do',1)
  IF (nks/=1) CALL errore('cg_readin','too many k-points',1)
  !      if (xk(1,1).ne.0.0 .or. xk(2,1).ne.0.0 .or. xk(3,1).ne.0.0)
  !     &    call errore('data','only k=0 allowed',1)
  IF (nmodes>3*nat .or. nmodes<0)                             &
       &     CALL errore('cg_readin','wrong number of normal modes',1)
  IF (epsil .and. nmodes/=0) CALL errore('cg_readin','not allowed',1)
  !
  IF (raman .and. deltatau<=0.d0)                                 &
       &     CALL errore('cg_readin','deltatau > 0 needed for raman CS',1)
  IF (nderiv/=2 .and. nderiv/=4) &
       CALL errore('cg_readin','nderiv not allowed',1)
  !
  IF (last==0) last=3*nat
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
  USE symm_base, ONLY : nsym, s, irt
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE cgcom
  !
  IMPLICIT NONE
  !
  INTEGER :: iunit
  !
  INTEGER :: na, nu, mu
  REAL(DP) utest, unorm, ddot
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
  IF (nmodes==0) THEN
     CALL find_equiv_sites (nat,nsym,irt,has_equivalent,        &
          &      n_diff_sites,n_equiv_atoms,equiv_atoms)
     IF (n_diff_sites <= 0 .or. n_diff_sites > nat)            &
          &      CALL errore('equiv.sites','boh!',1)
     !
     ! these are all modes, but only independent modes are calculated
     !
     nmodes = 3*nat
     u(:,:) = 0.d0
     DO nu = 1,nmodes
        u(nu,nu) = 1.0d0
     ENDDO
     ! look if ASR can be exploited to reduce the number of calculations
     ! we need to locate an independent atom with no equivalent atoms
     nasr=0
     IF (asr.and.n_diff_sites>1) THEN
        DO na = 1, n_diff_sites
           IF (n_equiv_atoms(na)==1 ) THEN
              nasr = equiv_atoms(na, 1)
              GOTO 1
           ENDIF
        ENDDO
 1      CONTINUE
     ENDIF
  ELSE
     IF (asr) CALL infomsg ('readin','warning: asr disabled')
     nasr=0
     !
     ! ... otherwise read normal modes from input
     !
     DO na = 1,nat
        has_equivalent(na) = 0
     ENDDO

     IF ( ionode ) THEN
        !
        DO nu = 1,nmodes
           READ (iunit,*,END=10,err=10) (u(mu,nu), mu=1,3*nat)
        ENDDO
        !
     ENDIF
     CALL mp_bcast(u,ionode_id,world_comm)
     DO nu = 1,nmodes
        DO mu = 1, nu-1
           utest = ddot(3*nat,u(1,nu),1,u(1,mu),1)
           IF (abs(utest)>1.0d-10) THEN
              PRINT *, ' warning: input modes are not orthogonal'
              CALL daxpy(3*nat,-utest,u(1,mu),1,u(1,nu),1)
           ENDIF
        ENDDO
        unorm = sqrt(ddot(3*nat,u(1,nu),1,u(1,nu),1))
        IF (abs(unorm)<1.0d-10) GOTO 10
        CALL dscal(3*nat,1.0d0/unorm,u(1,nu),1)
     ENDDO
     GOTO 20
10   CALL errore('phonon','wrong data read',1)
  ENDIF
20 CONTINUE
  !
  RETURN
END SUBROUTINE cg_readmodes
