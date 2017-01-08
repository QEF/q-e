!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
PROGRAM do_dos
  !--------------------------------------------------------------------
  !
  ! Calculates the Density of States (DOS),
  ! separated into up and down components for LSDA
  !
  ! See files INPUT_DOS.* in Doc/ directory for usage
  ! IMPORTANT: since v.5 namelist name is &dos and no longer &inputpp
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : prefix, tmp_dir
  USE constants,  ONLY : rytoev
  USE ener,       ONLY : ef, ef_up, ef_dw 
  USE kinds,      ONLY : DP
  USE klist,      ONLY : xk, wk, degauss, ngauss, lgauss, ltetra, nks, nkstot,&
                         two_fermi_energies
  USE wvfct,      ONLY : nbnd, et
  USE lsda_mod,   ONLY : lsda, nspin
  USE noncollin_module, ONLY: noncolin
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE mp_global,     ONLY : mp_startup
  USE environment,   ONLY : environment_start, environment_end
  ! following modules needed for generation of tetrahedra
  USE ktetra,     ONLY : tetra, tetra_type, tetra_init, tetra_dos_t, &
       opt_tetra_init, opt_tetra_dos_t
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev
  USE cell_base,  ONLY : at, bg
  USE start_k,    ONLY : k1, k2, k3, nk1, nk2, nk3
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: fildos, outdir
  CHARACTER(LEN=33) :: fermi_str
  REAL(DP) :: E, DOSofE (2), DOSint, DeltaE, Emin, Emax, &
              degauss1, E_unset=1000000.d0
  INTEGER :: nks2, n, ndos, ngauss1, ios

  NAMELIST /dos/ outdir, prefix, fildos, degauss, ngauss, &
       Emin, Emax, DeltaE
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'DOS' )
  !
  ios = 0
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     fildos =' '
     Emin   =-E_unset
     Emax   = E_unset
     DeltaE = 0.01d0
     ngauss = 0
     degauss= 0.d0
     !
     CALL input_from_file ( )
     !
     READ (5, dos, iostat=ios )
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1 = degauss
     ngauss1  = ngauss
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios /= 0 ) CALL errore('dos','reading dos namelist',abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  !
  CALL read_xml_file( )
  !
  IF ( ionode ) THEN
     !
     IF (nks /= nkstot) &
        CALL errore ('dos', 'pools not implemented, or incorrect file read', 1)
     !
     IF (degauss1/=0.d0) THEN
        degauss=degauss1
        ngauss =ngauss1
        WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
        ltetra=.false.
        lgauss=.true.
     ELSEIF (ltetra) THEN
        !
        ! info on tetrahedra is no longer saved to file and must be rebuilt
        !
        ! workaround for old xml file, to be removed
        IF ( ALLOCATED ( tetra ) ) DEALLOCATE (tetra)
        ! in the lsda case, only the first half of the k points
        ! are needed in the input of "tetrahedra"
        !
        IF ( lsda ) THEN
           nks2 = nks/2
        ELSE
           nks2 = nks
        END IF
        !
        IF(tetra_type == 0) THEN
           WRITE( stdout,'(/5x,"Tetrahedra used"/)')
           CALL tetra_init ( nsym, s, time_reversal, t_rev, at, bg, nks, &
                k1,k2,k3, nk1,nk2,nk3, nks2, xk )
        ELSE
           IF(tetra_type == 1) THEN 
              WRITE( stdout,'(/5x,"Linear tetrahedron method is used"/)')
           ELSE
              WRITE( stdout,'(/5x,"Optimized tetrahedron method used"/)')
           END IF
           CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, nks, &
                &                k1, k2, k3, nk1, nk2, nk3, nks2, xk, 1)
           !
        END IF
        !
     ELSEIF (lgauss) THEN
        WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ELSE
        degauss=DeltaE/rytoev
        ngauss =0
        WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
        ltetra=.false.
        lgauss=.true.
     ENDIF
     !
     ! find min and max energy for plot (band extrema if not set)
     !
     IF ( Emin == -E_unset ) THEN
        Emin = MINVAL ( et(1, 1:nks) )
        IF ( degauss > 0.0_dp ) Emin = Emin - 3.0_dp * degauss
     ELSE
        Emin = Emin/rytoev
     END IF
     IF ( Emax  == E_unset ) THEN
        Emax = MINVAL ( et(nbnd, 1:nks) )
        IF ( degauss > 0.0_dp ) Emax = Emax + 3.0_dp * degauss
     ELSE 
        Emax = Emax/rytoev
     END IF
     !
     DeltaE = DeltaE / rytoev
     ndos = nint ( (Emax - Emin) / DeltaE+0.500001d0)
     DOSint = 0.d0
     !
     IF ( fildos == ' ' ) fildos = trim(prefix)//'.dos'
     OPEN (unit = 4, file = fildos, status = 'unknown', form = 'formatted')
     IF ( two_fermi_energies ) THEN
        WRITE(fermi_str,'(" EFermi = ",2f7.3," eV")') ef_up*rytoev, ef_dw*rytoev
     ELSE
        WRITE(fermi_str,'(" EFermi = ",f7.3," eV")') ef*rytoev
     ENDIF

     IF (nspin==1.or.nspin==4) THEN
        WRITE(4,'("#  E (eV)   dos(E)     Int dos(E)",A)') TRIM(fermi_str)
     ELSE
        WRITE(4,'("#  E (eV)   dosup(E)     dosdw(E)   Int dos(E)",A)') &
        &          TRIM(fermi_str)
     ENDIF
     !
     DO n= 1, ndos
        E = Emin + (n - 1) * DeltaE
        IF (ltetra) THEN
           IF (tetra_type == 0) THEN
              CALL tetra_dos_t( et, nspin, nbnd, nks, E, DOSofE)
           ELSE
              CALL opt_tetra_dos_t( et, nspin, nbnd, nks, E, DOSofE)
           END IF
        ELSE
           CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, E, DOSofE)
        ENDIF
        IF (nspin==1.or.nspin==4) THEN
           DOSint = DOSint + DOSofE (1) * DeltaE
           WRITE (4, '(f8.3,2e12.4)') E * rytoev, DOSofE(1)/rytoev, DOSint
        ELSE
           DOSint = DOSint + (DOSofE (1) + DOSofE (2) ) * DeltaE
           WRITE (4, '(f8.3,3e12.4)') E * rytoev, DOSofE/rytoev, DOSint
        ENDIF
     ENDDO

     CLOSE (unit = 4)
     !
  ENDIF
  !
  CALL environment_end ( 'DOS' )
  !
  CALL stop_pp
  !
END PROGRAM do_dos

