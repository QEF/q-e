!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE summary_1drism()
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output the 1D-RISM's information obtained
  !    from the input file and from the setup routine, before starting the
  !    iterative calculation.
  !
  !    if iverbosity < 1 only a partial summary is done.
  !
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_size
  USE rism,          ONLY : CLOSURE_HNC, CLOSURE_KH
  USE rism1d_facade, ONLY : rism1t, niter, epsv, bond_width, &
                          & dielectric, molesize, mdiis_size, mdiis_step
  USE solvmol,       ONLY : get_nsite_in_solVs
  !
  IMPLICIT NONE
  !
  INTEGER          :: nsite
  INTEGER          :: ngrid
  CHARACTER(LEN=3) :: sclosure
  !
  IF (.NOT. rism1t%is_intra) THEN
    RETURN
  END IF
  !
  CALL print_solv_info(iverbosity)
  !
  IF (rism1t%closure == CLOSURE_HNC) THEN
    sclosure = 'HNC'
  ELSE IF (rism1t%closure == CLOSURE_KH) THEN
    sclosure = 'KH '
  ELSE
    sclosure = '???'
  END IF
  !
  nsite = get_nsite_in_solVs()
  !
  ngrid = rism1t%rfft%ngrid
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"1D-RISM info")')
  WRITE(stdout, '(5X,"------------")')
  WRITE(stdout, '(5X,"closure equation        = ",A12)')               TRIM(sclosure)
  WRITE(stdout, '(5X,"temperature             = ",F12.4,"  kelvin")')  rism1t%temp
  WRITE(stdout, '(5X,"coulomb smearing radius = ",F12.4,"  bohr")')    rism1t%tau
  WRITE(stdout, '(5X,"number of solvent sites = ",I12)')               nsite
  IF (iverbosity > 0) THEN
  WRITE(stdout, '(5X,"nv * (nv + 1) / 2       = ",I12)')               rism1t%nsite
  END IF
  WRITE(stdout, '(5X,"number of grids         = ",I12)')               ngrid
  WRITE(stdout, '(5X,"maximum of R-space      = ",F12.4,"  bohr")')    rism1t%rfft%rgrid(ngrid)
  WRITE(stdout, '(5X,"maximum of G-space      = ",F12.4,"  bohr^-1")') rism1t%rfft%ggrid(ngrid)
  IF (iverbosity > 0) THEN
  WRITE(stdout, '(5X,"#R-grids in local       = ",I12)')               rism1t%nr
  WRITE(stdout, '(5X,"#G-grids in local       = ",I12)')               rism1t%ng
  END IF
  WRITE(stdout, '(5X,"number of iterations    = ",I12)')               niter
  WRITE(stdout, '(5X,"convergence threshold   = ",1PE12.1)')           epsv
  WRITE(stdout, '(5X,"Gaussian width of bonds = ",F12.4,"  bohr")')    bond_width
  WRITE(stdout, '(5X,"size of MDIIS           = ",I12)')               mdiis_size
  WRITE(stdout, '(5X,"step of MDIIS           = ",0PF12.4)')           mdiis_step
  WRITE(stdout, '(5X,"number of processes     = ",I12)')               mp_size(rism1t%intra_comm)
  IF (dielectric > 0.0_DP) THEN
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"--- Dielectrically Consistent RISM is used. ---")')
  WRITE(stdout, '(5X,"dielectric constant     = ",F12.4)')             dielectric
  WRITE(stdout, '(5X,"size of molecule        = ",F12.4,"  bohr")')    molesize
  END IF
  WRITE(stdout, '()')
  !
  CALL print_radfft_info(iverbosity)
  !
  CALL print_mp_rism1d_info(iverbosity)
  !
  FLUSH(stdout)
  !
END SUBROUTINE summary_1drism
!
!-----------------------------------------------------------------------
SUBROUTINE print_solv_info(iverbosity)
  !-----------------------------------------------------------------------
  !
  USE cell_base,      ONLY : omega
  USE constants,      ONLY : BOHR_RADIUS_ANGS, BOHR_RADIUS_SI, ELECTRON_SI, AU_DEBYE, eps32
  USE io_files,       ONLY : pseudo_dir, molfile
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : BOHRm3_TO_MOLLm1, BOHRm3_TO_MOLCMm3, RY_TO_KCALMOLm1
  USE solvmol,        ONLY : nsolV, solVs, get_nsite_in_solVs, get_nuniq_in_solVs, &
                           & isite_to_isolV, isite_to_iatom, iuniq_to_nsite, iuniq_to_isite
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  !
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: nsite, msite
  INTEGER  :: nuniq, iuniq
  INTEGER  :: irho, nrho
  REAL(DP) :: rhov
  REAL(DP) :: dipv
  !
  DO isolV = 1, nsolV
    !
    WRITE(stdout, '()')
    !
    WRITE(stdout, '(5X,"Molecule #",I2," for ",A," read from file:")') isolV, TRIM(solVs(isolV)%name)
    WRITE(stdout, '(5X,A)') TRIM(pseudo_dir) // TRIM (molfile(isolV))
    !
    IF (ABS(solVs(isolV)%density - solVs(isolV)%subdensity) < eps32) THEN
      nrho = 1
    ELSE
      nrho = 2
    END IF
    !
    DO irho = 1, nrho
      IF (irho == 1) THEN
        rhov = solVs(isolV)%density
      ELSE
        rhov = solVs(isolV)%subdensity
      END IF
      !
      IF (nrho == 1) THEN
        WRITE(stdout, '(5X,"Density:")')
      ELSE IF (irho == 1) THEN
        WRITE(stdout, '(5X,"Density (the right-hand side):")')
      ELSE
        WRITE(stdout, '(5X,"Density (the left-hand side):")')
      END IF
      WRITE(stdout, '(5X,2X,E16.8," cell^-1")') rhov * omega
      WRITE(stdout, '(5X,2X,E16.8," bohr^-3")') rhov
      WRITE(stdout, '(5X,2X,E16.8," mol/L")')   rhov * BOHRm3_TO_MOLLm1
      WRITE(stdout, '(5X,2X,E16.8," g/cm^3")')  rhov * solVs(isolV)%mass * BOHRm3_TO_MOLCMm3
    END DO
    !
    IF (solVs(isolV)%permittivity > 0.0_DP) THEN
      WRITE(stdout, '(5X,"Permittivity:")')
      WRITE(stdout, '(5X,F12.6)') solVs(isolV)%permittivity
    END IF
    !
    IF (solVs(isolV)%is_polar) THEN
      dipv = solVs(isolV)%dipole
      WRITE(stdout, '(5X,"Dipole-moment:")')
      WRITE(stdout, '(5X,2X,E16.8," e*bohr")') dipv
      WRITE(stdout, '(5X,2X,E16.8," debye")')  dipv * AU_DEBYE
      WRITE(stdout, '(5X,2X,E16.8," C*m")')    dipv * ELECTRON_SI * BOHR_RADIUS_SI
    END IF
    !
    WRITE(stdout, '(5X,"Number of atoms: ",I3)') solVs(isolV)%natom
    !
    WRITE(stdout, '(5X,"Atoms:")')
    WRITE(stdout, '(5X,A)') &
      & '  #  atom  ' // &
      & '    X (angs)      Y (angs)      Z (angs)      Q (e)     ' // &
      & '    E (kcal/mol)  S (angs)'
    !
    DO iatom = 1, solVs(isolV)%natom
      WRITE(stdout, '(5X,I3,2X,A6,6F14.8)') iatom, &
        & solVs(isolV)%aname(iatom) // '    ', &
        & solVs(isolV)%coord(1:3, iatom) * BOHR_RADIUS_ANGS, &
        & solVs(isolV)%charge(iatom), &
        & solVs(isolV)%ljeps(iatom) * RY_TO_KCALMOLm1, &
        & solVs(isolV)%ljsig(iatom) * BOHR_RADIUS_ANGS
    END DO
    !
    WRITE(stdout, '()')
    !
  END DO
  !
  IF (iverbosity > 0) THEN
    !
    nsite = get_nsite_in_solVs()
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Number of total sites: ",I3)') nsite
    WRITE(stdout, '(5X,"Index site -> solvent: ")')
    WRITE(stdout, '(5X,20I3)') isite_to_isolV
    WRITE(stdout, '(5X,"Index site -> atom (in a solvent): ")')
    WRITE(stdout, '(5X,20I3)') isite_to_iatom
    WRITE(stdout, '()')
    !
    nuniq = get_nuniq_in_solVs()
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Number of unique sites: ",I3)') nuniq
    WRITE(stdout, '(5X,"Multiplicity of unique site: ")')
    WRITE(stdout, '(5X,20I3)') iuniq_to_nsite
    WRITE(stdout, '(5X,"Index unique site -> site: ")')
    DO iuniq = 1, nuniq
      msite = iuniq_to_nsite(iuniq)
      WRITE(stdout, '(5X,I3,":",20I3)') iuniq, iuniq_to_isite(1:msite, iuniq)
    END DO
    WRITE(stdout, '()')
    !
  END IF
  !
END SUBROUTINE print_solv_info
!
!-----------------------------------------------------------------------
SUBROUTINE print_radfft_info(iverbosity)
  !-----------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout
  USE rism1d_facade, ONLY : rism1t
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  !
  INTEGER :: nt
  INTEGER :: nt1
  INTEGER :: nt2
  !
  INTEGER, PARAMETER :: MT = 10
  !
  IF (iverbosity < 1) THEN
    RETURN
  END IF
  !
  WRITE(stdout, '()')
  !
  WRITE(stdout, '(5X,"Radial FFT:")')
  WRITE(stdout, '(5X,"number of radial grids   = ",I12)') rism1t%rfft%ngrid
  WRITE(stdout, '(5X,"number of FFT grids      = ",I12)') rism1t%rfft%mgrid
  WRITE(stdout, '(5X,"good number of FFT grids = ",I12)') rism1t%rfft%lgrid
  !
  nt  = rism1t%rfft%ngrid
  nt1 = MIN(MT, nt)
  nt2 = MAX(nt - MT + 1, nt1 + 1)
  !
  WRITE(stdout, '(5X,"R-space grids:")')
  WRITE(stdout, '(5X,5E16.8)') rism1t%rfft%rgrid(1:nt1)
  IF (nt2 <= nt) THEN
    WRITE(stdout, '(5X,5("  .............."))')
    WRITE(stdout, '(5X,5E16.8)') rism1t%rfft%rgrid(nt2:nt)
  END IF
  !
  WRITE(stdout, '(5X,"G-space grids:")')
  WRITE(stdout, '(5X,5E16.8)') rism1t%rfft%ggrid(1:nt1)
  IF (nt2 <= nt) THEN
    WRITE(stdout, '(5X,5("  .............."))')
    WRITE(stdout, '(5X,5E16.8)') rism1t%rfft%ggrid(nt2:nt)
  END IF
  !
  WRITE(stdout, '()')
  !
END SUBROUTINE print_radfft_info
!
!-----------------------------------------------------------------------
SUBROUTINE print_mp_rism1d_info(iverbosity)
  !-----------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout
  USE rism1d_facade, ONLY : rism1t
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
  WRITE(stdout, '(5X,"number of site groups     = ",I12)') rism1t%mp_site%nsitg
  WRITE(stdout, '(5X,"#procs in a site group    = ",I12)') rism1t%mp_site%nproc_sitg
  WRITE(stdout, '(5X,"this proc in a site group = ",I12)') rism1t%mp_site%me_sitg
  WRITE(stdout, '(5X,"the root in a site group  = ",I12)') rism1t%mp_site%root_sitg
  WRITE(stdout, '(5X,"this site group           = ",I12)') rism1t%mp_site%my_sitg_id
  WRITE(stdout, '(5X,"inter-site group comm.    = ",I12)') rism1t%mp_site%inter_sitg_comm
  WRITE(stdout, '(5X,"intra-site group comm.    = ",I12)') rism1t%mp_site%intra_sitg_comm
  WRITE(stdout, '(5X,"total number of sites     = ",I12)') rism1t%mp_site%nsite
  WRITE(stdout, '(5X,"starting index of sites   = ",I12)') rism1t%mp_site%isite_start
  WRITE(stdout, '(5X,"ending index of sites     = ",I12)') rism1t%mp_site%isite_end
  WRITE(stdout, '()')
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"MPI for task:")')
  WRITE(stdout, '(5X,"#procs in a task group    = ",I12)') rism1t%mp_task%nproc_task
  WRITE(stdout, '(5X,"this proc in a task group = ",I12)') rism1t%mp_task%me_task
  WRITE(stdout, '(5X,"the root in a task group  = ",I12)') rism1t%mp_task%root_task
  WRITE(stdout, '(5X,"task group comm.          = ",I12)') rism1t%mp_task%itask_comm
  WRITE(stdout, '(5X,"total number of vectors   = ",I12)') rism1t%mp_task%nvec
  WRITE(stdout, '(5X,"starting index of vectors = ",I12)') rism1t%mp_task%ivec_start
  WRITE(stdout, '(5X,"ending index of vectors   = ",I12)') rism1t%mp_task%ivec_end
  WRITE(stdout, '(5X,"lengths of vectors        = ")')
  WRITE(stdout, '(5X,10I5)')  rism1t%mp_task%ilen_vecs
  WRITE(stdout, '(5X,"displacement of vectors   = ")')
  WRITE(stdout, '(5X,10I5)')  rism1t%mp_task%idis_vecs
  WRITE(stdout, '()')
  !
END SUBROUTINE print_mp_rism1d_info
