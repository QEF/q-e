!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE summary_3drism()
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output the 3D-RISM's information obtained
  !    from the input file and from the setup routine, before starting the
  !    iterative calculation.
  !
  !    if iverbosity < 1 only a partial summary is done.
  !
  USE cell_base,     ONLY : at, alat
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE rism,          ONLY : ITYPE_LAUERISM, CLOSURE_HNC, CLOSURE_KH
  USE rism3d_facade, ONLY : rism3t, niter, epsv, conv_level, &
                          & mdiis_size, mdiis_step, ecutsolv, laue_nfit, &
                          & ireference, IREFERENCE_AVERAGE, IREFERENCE_RIGHT, IREFERENCE_LEFT
  USE solute,        ONLY : iwall, wall_tau, IWALL_RIGHT, IWALL_LEFT
  USE solvmol,       ONLY : solVs, get_nuniq_in_solVs, iuniq_to_isite, &
                          & isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  LOGICAL                       :: laue
  INTEGER                       :: nsite
  INTEGER                       :: isite
  INTEGER                       :: iisite
  INTEGER                       :: isolV
  INTEGER                       :: iatom
  REAL(DP)                      :: zstart
  REAL(DP)                      :: zend
  REAL(DP)                      :: zstep
  REAL(DP)                      :: zedge1
  REAL(DP)                      :: zedge2
  CHARACTER(LEN=8), ALLOCATABLE :: asite(:)
  CHARACTER(LEN=3)              :: sclosure
  CHARACTER(LEN=32)             :: sbound
  CHARACTER(LEN=32)             :: sreference
  !
  CALL print_solu_info(iverbosity)
  !
  IF (rism3t%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  ELSE
    laue = .FALSE.
  END IF
  !
  IF (rism3t%closure == CLOSURE_HNC) THEN
    sclosure = 'HNC'
  ELSE IF (rism3t%closure == CLOSURE_KH) THEN
    sclosure = 'KH '
  ELSE
    sclosure = '???'
  END IF
  !
  IF (laue) THEN
    IF (rism3t%lfft%xleft .AND. rism3t%lfft%xright) THEN
      sbound = 'Solvent-Slab-Solvent'
    ELSE IF (rism3t%lfft%xleft) THEN
      sbound = 'Solvent-Slab-Vacuum'
    ELSE IF (rism3t%lfft%xright) THEN
      sbound = 'Vacuum-Slab-Solvent'
    ELSE
      sbound = 'Vacuum-Slab-Vacuum'
    END IF
    !
    IF (ireference == IREFERENCE_AVERAGE) THEN
      sreference = 'average of both sides'
    ELSE IF (ireference == IREFERENCE_RIGHT) THEN
      sreference = 'right-hand side'
    ELSE IF (ireference == IREFERENCE_LEFT) THEN
      sreference = 'left-hand side'
    ELSE
      sreference = 'not defined'
    END IF
    !
    zstep  = rism3t%lfft%zstep
    zedge1 = rism3t%lfft%zleft
    zedge2 = rism3t%lfft%zleft + rism3t%lfft%zstep
  END IF
  !
  nsite = get_nuniq_in_solVs()
  !
  ALLOCATE(asite(nsite))
  DO isite = 1, nsite
    iisite = iuniq_to_isite(1, isite)
    isolV  = isite_to_isolV(iisite)
    iatom  = isite_to_iatom(iisite)
    asite(isite) = TRIM(ADJUSTL(solVs(isolV)%aname(iatom)))
    asite(isite) = ADJUSTR(asite(isite))
  END DO
  !
  WRITE(stdout, '()')
  IF (.NOT. laue) THEN
  WRITE(stdout, '(5X,"3D-RISM info")')
  WRITE(stdout, '(5X,"------------")')
  ELSE
  WRITE(stdout, '(5X,"Laue-RISM info")')
  WRITE(stdout, '(5X,"--------------")')
  END IF
  WRITE(stdout, '(5X,"closure equation        = ",A12)')               TRIM(sclosure)
  WRITE(stdout, '(5X,"temperature             = ",F12.4,"  kelvin")')  rism3t%temp
  WRITE(stdout, '(5X,"coulomb smearing radius = ",F12.4,"  bohr")')    rism3t%tau
  WRITE(stdout, '(5X,"number of solvent sites = ",I12)')               nsite
  WRITE(stdout, '(5X,"solvent sites           = ",4A8)')               asite
  IF (iverbosity > 0) THEN
  WRITE(stdout, '(5X,"#R-grids in local       = ",I12)')               rism3t%nr
  IF (laue) THEN
  WRITE(stdout, '(5X,"#Rz-grids (unit-cell)   = ",I12)')               rism3t%nrzs
  WRITE(stdout, '(5X,"#Rz-grids (expand-cell) = ",I12)')               rism3t%nrzl
  END IF
  WRITE(stdout, '(5X,"#G-grids in local       = ",I12)')               rism3t%ng
  IF (.NOT. laue) THEN
  WRITE(stdout, '(5X,"#G-shells in local      = ",I12)')               rism3t%ngs
  ELSE
  WRITE(stdout, '(5X,"#Gxy-grids in local     = ",I12)')               rism3t%ngxy
  WRITE(stdout, '(5X,"#Gxy-shells in local    = ",I12)')               rism3t%ngs
  END IF
  END IF
  WRITE(stdout, '(5X,"number of iterations    = ",I12)')               niter
  WRITE(stdout, '(5X,"convergence threshold   = ",1PE12.1)')           epsv
  WRITE(stdout, '(5X,"convergence level       = ",0PF12.4)')           conv_level
  WRITE(stdout, '(5X,"size of MDIIS           = ",I12)')               mdiis_size
  WRITE(stdout, '(5X,"step of MDIIS           = ",0PF12.4)')           mdiis_step
  WRITE(stdout, '(5X,"solvent cutoff          = ",F12.4,"  Ry")')      ecutsolv
  WRITE(stdout, '(5X,"Custom grid: ",I8," G-vectors",5X,"FFT dimensions: (",I4,",",I4,",",I4,")")') &
    &                                                                  rism3t%gvec%ngm_g,    &
    &                                                                  rism3t%dfft%nr1, &
    &                                                                  rism3t%dfft%nr2, &
    &                                                                  rism3t%dfft%nr3
  IF (laue) THEN
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Boundary Conditions        : ",A)')              TRIM(sbound)
  WRITE(stdout, '(5X,"unit-cell (in bohr)        : [",2F11.6,"]")')    (-0.5_DP * at(3, 3) * alat), &
    &                                                                  (+0.5_DP * at(3, 3) * alat)
  WRITE(stdout, '(5X,"expand-cell (in bohr)      : [",2F11.6,"]")')    (rism3t%lfft%zleft  * alat), &
    &                                                                  (rism3t%lfft%zright * alat)
  IF (rism3t%lfft%xleft) THEN
  zstart = zedge1
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%izleft_gedge - 1)
  WRITE(stdout, '(5X,"solvent of left (in bohr)  : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
  IF (rism3t%lfft%izleft_gedge /= rism3t%lfft%izleft_end) THEN
  zstart = zedge1 + zstep * DBLE(rism3t%lfft%izleft_gedge    )
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%izleft_end   - 1)
  WRITE(stdout, '(5X,"buffer  of left (in bohr)  : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
  END IF
#if !defined (__DEBUG_RISM)
  IF (iverbosity > 0) THEN
#endif
  zstart = zedge1 + zstep * DBLE(rism3t%lfft%izleft_start0 - 1)
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%izleft_end0   - 1)
  WRITE(stdout, '(5X,"Gxy = 0 of left (in bohr)  : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
#if !defined (__DEBUG_RISM)
  END IF
#endif
  END IF
  IF (rism3t%lfft%xright) THEN
  zstart = zedge1 + zstep * DBLE(rism3t%lfft%izright_gedge - 1)
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%nrz           - 1)
  WRITE(stdout, '(5X,"solvent of right (in bohr) : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
  IF (rism3t%lfft%izright_start /= rism3t%lfft%izright_gedge) THEN
  zstart = zedge1 + zstep * DBLE(rism3t%lfft%izright_start - 1)
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%izright_gedge - 2)
  WRITE(stdout, '(5X,"buffer  of right (in bohr) : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
  END IF
#if !defined (__DEBUG_RISM)
  IF (iverbosity > 0) THEN
#endif
  zstart = zedge1 + zstep * DBLE(rism3t%lfft%izright_start0 - 1)
  zend   = zedge2 + zstep * DBLE(rism3t%lfft%izright_end0   - 1)
  WRITE(stdout, '(5X,"Gxy = 0 of right (in bohr) : [",2F11.6,"]")')    (zstart * alat), &
    &                                                                  (zend   * alat)
#if !defined (__DEBUG_RISM)
  END IF
#endif
  END IF
  IF (iwall == IWALL_RIGHT) THEN
  WRITE(stdout, '(5X,"repulsive wall (in bohr)   : [",F11.6,A11,"]")') (wall_tau * alat), "+Infinity"
  ELSE IF (iwall == IWALL_LEFT) THEN
  WRITE(stdout, '(5X,"repulsive wall (in bohr)   : [",A11,F11.6,"]")') "-Infinity", (wall_tau * alat)
  END IF
  WRITE(stdout, '(5X,"reference of potential     : ",A)')              TRIM(sreference)
  WRITE(stdout, '(5X,"#grids to fit at edges of unit-cell =",I5)')     laue_nfit
  WRITE(stdout, '(5X,"Laue grid: ",I8," Gxy-vectors",5X,"FFT dimensions: (",I4,",",I4,",",I4,")")') &
    &                                                                  rism3t%lfft%ngxy_g,    &
    &                                                                  rism3t%dfft%nr1, &
    &                                                                  rism3t%dfft%nr2, &
    &                                                                  rism3t%lfft%nrz
  END IF
  WRITE(stdout, '()')
  !
  DEALLOCATE(asite)
  !
  CALL print_mp_rism3d_info(iverbosity)
  !
  FLUSH(stdout)
  !
END SUBROUTINE summary_3drism
!
!-----------------------------------------------------------------------
SUBROUTINE print_solu_info(iverbosity)
  !-----------------------------------------------------------------------
  !
  USE ions_base,      ONLY : nat, atm, ityp
  USE constants,      ONLY : BOHR_RADIUS_ANGS
  USE io_global,      ONLY : stdout
  USE molecule_const, ONLY : RY_TO_KCALMOLm1
  USE solute,         ONLY : solU_ljeps, solU_ljsig, solU_ljname, rmax_lj, &
                             iwall, wall_rho, wall_ljeps, wall_ljsig, wall_lj6, IWALL_NULL
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  !
  INTEGER :: ia
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Solute:")')
  WRITE(stdout, '(5X,A)') '  # atom     E (kcal/mol)  S (angs)'
  !
  DO ia = 1, nat
    !
    WRITE(stdout, '(5X,I3,2X,A4,2F14.8,4X,A)') ia, &
      & ADJUSTL(atm(ityp(ia))) // '    ', &
      & solU_ljeps(ia) * RY_TO_KCALMOLm1, &
      & solU_ljsig(ia) * BOHR_RADIUS_ANGS, &
      & TRIM(ADJUSTL(solU_ljname(ia)))
    !
  END DO
  !
  IF (iwall /= IWALL_NULL) THEN
    !
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Wall:")')
    WRITE(stdout, '(5X,"  Density =",F14.8,"  bohr^-3")')  wall_rho
    WRITE(stdout, '(5X,"  E       =",F14.8,"  kcal/mol")') wall_ljeps * RY_TO_KCALMOLm1
    WRITE(stdout, '(5X,"  S       =",F14.8,"  angs")')     wall_ljsig * BOHR_RADIUS_ANGS
    IF (wall_lj6) THEN
    WRITE(stdout, '(5X,"  the attractive term of -(1/r)^6 is used")')
    ELSE
    END IF
    !
  END IF
  !
  WRITE(stdout, '()')
  !
  IF (iverbosity > 0) THEN
    WRITE(stdout, '(5X,"L.J. cutoff radius = ",F12.4,"  S")') rmax_lj
    WRITE(stdout, '()')
  END IF
  !
END SUBROUTINE print_solu_info
!
!-----------------------------------------------------------------------
SUBROUTINE print_mp_rism3d_info(iverbosity)
  !-----------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout
  USE rism3d_facade, ONLY : rism3t
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  !
  IF (iverbosity < 1) THEN
    RETURN
  END IF
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"MPI for site:")')
  WRITE(stdout, '(5X,"number of site groups     = ",I12)') rism3t%mp_site%nsitg
  WRITE(stdout, '(5X,"#procs in a site group    = ",I12)') rism3t%mp_site%nproc_sitg
  WRITE(stdout, '(5X,"this proc in a site group = ",I12)') rism3t%mp_site%me_sitg
  WRITE(stdout, '(5X,"the root in a site group  = ",I12)') rism3t%mp_site%root_sitg
  WRITE(stdout, '(5X,"this site group           = ",I12)') rism3t%mp_site%my_sitg_id
  WRITE(stdout, '(5X,"inter-site group comm.    = ",I12)') rism3t%mp_site%inter_sitg_comm
  WRITE(stdout, '(5X,"intra-site group comm.    = ",I12)') rism3t%mp_site%intra_sitg_comm
  WRITE(stdout, '(5X,"total number of sites     = ",I12)') rism3t%mp_site%nsite
  WRITE(stdout, '(5X,"starting index of sites   = ",I12)') rism3t%mp_site%isite_start
  WRITE(stdout, '(5X,"ending index of sites     = ",I12)') rism3t%mp_site%isite_end
  WRITE(stdout, '()')
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"MPI for task:")')
  WRITE(stdout, '(5X,"#procs in a task group    = ",I12)') rism3t%mp_task%nproc_task
  WRITE(stdout, '(5X,"this proc in a task group = ",I12)') rism3t%mp_task%me_task
  WRITE(stdout, '(5X,"the root in a task group  = ",I12)') rism3t%mp_task%root_task
  WRITE(stdout, '(5X,"task group comm.          = ",I12)') rism3t%mp_task%itask_comm
  WRITE(stdout, '()')
  !
END SUBROUTINE print_mp_rism3d_info
