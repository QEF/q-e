!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE read_pseudo_mod
!=----------------------------------------------------------------------------=!
  !
  !! read pseudopotential files and store the data on internal variables of the 
  !! program. Note that all processors read the same file!
  !
  USE io_files,     ONLY: pseudo_dir, pseudo_dir_cur, psfile, tmp_dir
  USE ions_base,    ONLY: ntyp => nsp
  !! global variables  required on input 
  !
  USE atom,         ONLY: msh, rgrid
  USE ions_base,    ONLY: zv
  USE uspp_param,   ONLY: upf, nvb
  USE uspp,         ONLY: okvan, nlcc_any
  !! global variables modified on output 
  ! 
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  PUBLIC :: readpp, check_order
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
SUBROUTINE readpp ( input_dft, printout, ecutwfc_pp, ecutrho_pp )
  !-----------------------------------------------------------------------
  !
  !! Reads PP files and puts the result into the "upf" structure of module uspp_param
  !! Sets  DFT to input_dft if present, to the value read in PP files otherwise
  !! Sets  number of valence electrons Zv, control variables okvan and nlcc_any,
  !! compatibility variable nvb
  !! Optionally returns cutoffs read from PP files into ecutwfc_pp, ecutrho_pp
  !
  USE kinds,        ONLY: DP
  USE mp,           ONLY: mp_bcast, mp_sum
  USE mp_images,    ONLY: intra_image_comm
  USE io_global,    ONLY: stdout, ionode, ionode_id
  USE pseudo_types, ONLY: pseudo_upf, deallocate_pseudo_upf
  USE funct,        ONLY: enforce_input_dft, set_dft_from_name, &
       get_iexch, get_icorr, get_igcx, get_igcc, get_inlc
  use radial_grids, ONLY: deallocate_radial_grid, nullify_radial_grid
  USE wrappers,     ONLY: md5_from_file, f_remove
  USE read_upf_v1_module,   ONLY: read_upf_v1
  USE read_upf_new_module,  ONLY: read_upf_new
  USE upf_auxtools, ONLY: upf_get_pp_format, upf_check_atwfc_norm
  USE upf_to_internal,  ONLY: add_upf_grid, set_upf_q
  USE read_uspp_module, ONLY: readvan, readrrkj
  USE m_gth,            ONLY: readgth
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(INOUT) :: input_dft
  LOGICAL, OPTIONAL, INTENT(IN) :: printout
  REAL(DP), OPTIONAL, INTENT(OUT) :: ecutwfc_pp, ecutrho_pp  
  !
  REAL(DP), parameter :: rcut = 10.d0 
  ! 2D Coulomb cutoff: modify this (at your own risks) if problems with cutoff 
  ! being smaller than pseudo rcut. original value=10.0
  CHARACTER(len=512) :: file_pseudo ! file name complete with path
  CHARACTER(len=512) :: file_fixed, msg
  LOGICAL :: printout_ = .FALSE., exst, is_xml
  INTEGER :: iunps, isupf, nt, nb, ir, ios
  INTEGER :: iexch_, icorr_, igcx_, igcc_, inlc_
  !
  ! ... initializations, allocations, etc
  !
  iunps = 4
  !
  IF( ALLOCATED( upf ) ) THEN
     DO nt = 1, SIZE( upf )
        CALL deallocate_pseudo_upf( upf( nt ) )
     END DO
     DEALLOCATE( upf )
  END IF
  !
  ALLOCATE ( upf( ntyp ) )
  !
  IF ( PRESENT(printout) ) THEN
     printout_ = printout .AND. ionode
  END IF
  IF ( printout_) THEN
     WRITE( stdout,"(//,3X,'Atomic Pseudopotentials Parameters',/, &
                   &    3X,'----------------------------------' )" )
  END IF
  !
  DO nt = 1, ntyp
     !
     ! try first pseudo_dir_cur if set: in case of restart from file,
     ! this is where PP files should be located
     !
     ios = 1
     IF ( pseudo_dir_cur /= ' ' ) THEN
        file_pseudo  = TRIM (pseudo_dir_cur) // TRIM (psfile(nt))
        INQUIRE(file = file_pseudo, EXIST = exst) 
        IF (exst) ios = 0
        CALL mp_sum (ios,intra_image_comm)
        IF ( ios /= 0 ) CALL infomsg &
                     ('readpp', 'file '//TRIM(file_pseudo)//' not found')
     END IF
     !
     ! file not found? no panic (yet): try the original location pseudo_dir
     ! as set in input (it should already contain a slash at the end)
     !
     IF ( ios /= 0 ) THEN
        file_pseudo = TRIM (pseudo_dir) // TRIM (psfile(nt))
        INQUIRE ( file = file_pseudo, EXIST = exst) 
        IF (exst) ios = 0
        CALL mp_sum (ios,intra_image_comm)
        CALL errore('readpp', 'file '//TRIM(file_pseudo)//' not found',ABS(ios))
     END IF
     !
     IF( printout_ ) THEN
        WRITE( stdout, "(/,3X,'Reading pseudopotential for specie # ',I2, &
                       & ' from file :',/,3X,A)") nt, TRIM(file_pseudo)
     END IF
     !
     IF ( ionode ) THEN
        isupf = 0
        CALL  read_upf_new( file_pseudo, upf(nt), isupf )
        !
        !! start reading - check  first if files are readable as xml files,
        !! then as UPF v.2, then as UPF v.1
        !
        IF (isupf ==-81 ) THEN
           !! error code -81 means that the file is not xml or UPF v.2 
           !! (the funny code value is for compatibility with FoX)
           CALL  read_upf_v1 (file_pseudo, upf(nt), isupf )
           !! try to read UPF v.1 file
           IF ( isupf == 0 ) isupf = -1
        END IF
        !
     END IF
     !
     CALL mp_bcast (isupf,ionode_id,intra_image_comm)
     !
     IF (isupf == -2 .OR. isupf == -1 .OR. isupf == 0) THEN
        !
        CALL upf_bcast(upf(nt), ionode, ionode_id, intra_image_comm)
        !! broadcast the pseudopotential to all processors
        !
        IF( printout_) THEN
           IF ( isupf == 0 ) THEN
              WRITE( stdout, "(3X,'file type is xml')") 
           ELSE
              WRITE( stdout, "(3X,'file type is UPF v.',I1)") ABS(isupf) 
           END IF
        END IF
        !
     ELSE
        !
        OPEN ( UNIT = iunps, FILE = file_pseudo, STATUS = 'old', FORM = 'formatted' )
        !
        !     The type of the pseudopotential is determined by the file name:
        !    *.xml or *.XML  UPF format with schema              pp_format=0
        !    *.upf or *.UPF  UPF format                          pp_format=1
        !    *.vdb or *.van  Vanderbilt US pseudopotential code  pp_format=2
        !    *.gth           Goedecker-Teter-Hutter NC pseudo    pp_format=3
        !    *.RRKJ3         Andrea's   US new code              pp_format=4
        !    none of the above: PWSCF norm-conserving format     pp_format=5
        !
        IF ( upf_get_pp_format( psfile(nt) ) == 2  ) THEN
           !
           IF( printout_ ) &
              WRITE( stdout, "(3X,'file type is Vanderbilt US PP')")
           CALL readvan (iunps, nt, upf(nt))
           !
        ELSE IF ( upf_get_pp_format( psfile(nt) ) == 3 ) THEN
           !
           IF( printout_ ) &
              WRITE( stdout, "(3X,'file type is GTH (analytical)')")
           CALL readgth (iunps, nt, upf(nt))
           !
        ELSE IF ( upf_get_pp_format( psfile(nt) ) == 4 ) THEN
           !
           IF( printout_ ) &
              WRITE( stdout, "(3X,'file type is RRKJ3')")
           CALL readrrkj (iunps, nt, upf(nt))
           !
        ELSE IF ( upf_get_pp_format( psfile(nt) ) == 5 ) THEN
           !
           IF( printout_ ) &
              WRITE( stdout, "(3X,'file type is old PWscf NC format')")
           CALL read_ncpp (iunps, nt, upf(nt))
           !
        ELSE
           !
           CALL errore('readpp', 'file '//TRIM(file_pseudo)//' not readable',1)
           !
        ENDIF
        !
        ! end of reading
        !
        CLOSE (iunps)
        !
     ENDIF
     !
     ! reconstruct Q(r) if needed
     !
     CALL set_upf_q (upf(nt))
     !
     ! Calculate MD5 checksum for this pseudopotential
     !
     CALL md5_from_file(file_pseudo, upf(nt)%md5_cksum)
     !
  END DO
  !
  ! end of PP reading - now set up more variables
  !
  ! radial grids - 
  !
  IF( ALLOCATED( rgrid ) ) THEN
     DO nt = 1, SIZE( rgrid )
        CALL deallocate_radial_grid( rgrid( nt ) )
        CALL nullify_radial_grid( rgrid( nt ) )
     END DO
     DEALLOCATE( rgrid )
     if(allocated(msh)) DEALLOCATE( msh )
  END IF
  ALLOCATE( rgrid( ntyp ), msh( ntyp ) )
  !
  nvb = 0
  DO nt = 1, ntyp
     !
     CALL nullify_radial_grid( rgrid( nt ) )
     CALL add_upf_grid (upf(nt), rgrid(nt))
     !
     ! the radial grid is defined up to r(mesh) but we introduce 
     ! an auxiliary variable msh to limit the grid up to rcut=10 a.u. 
     ! This is used to cut off the numerical noise arising from the
     ! large-r tail in cases like the integration of V_loc-Z/r
     !
     DO ir = 1, rgrid(nt)%mesh
        IF (rgrid(nt)%r(ir) > rcut) THEN
           msh (nt) = ir
           GOTO 5
        END IF
     END DO
     msh (nt) = rgrid(nt)%mesh 
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1
     !
     ! msh is forced to be odd for simpson integration (maybe obsolete?)
     !
     ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
     ! ... is set equal to Zp = pseudo-charge of the pseudopotential
     !
     zv(nt) = upf(nt)%zp
     !
     ! ... count US species (obsolete?)
     !
     IF (upf(nt)%tvanp) nvb=nvb+1
     !
     ! check for zero atomic wfc, 
     ! check that (occupied) atomic wfc are properly normalized
     !
     CALL upf_check_atwfc_norm(upf(nt),psfile(nt))
     !
  END DO
  !
  ! ... set DFT value
  !
  IF (input_dft /='none') CALL enforce_input_dft (input_dft)
  !
  DO nt = 1, ntyp
     !
     CALL set_dft_from_name( upf(nt)%dft )
     !
     ! ... Check for DFT consistency - ignored if dft enforced from input
     !
     IF (nt == 1) THEN
        iexch_ = get_iexch()
        icorr_ = get_icorr()
        igcx_  = get_igcx()
        igcc_  = get_igcc()
        inlc_  = get_inlc()
     ELSE
        IF ( iexch_ /= get_iexch() .OR. icorr_ /= get_icorr() .OR. &
             igcx_  /= get_igcx()  .OR. igcc_  /= get_igcc()  .OR.  &
             inlc_  /= get_inlc() ) THEN
           CALL errore( 'readpp','inconsistent DFT read from PP files', nt)
        END IF
     END IF
     !
  END DO
  !
  ! more initializations
  !
  okvan = ( nvb > 0 )
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  ! return cutoff read from PP file, if required
  !
  IF ( PRESENT(ecutwfc_pp) ) THEN
     ecutwfc_pp = MAXVAL ( upf(1:ntyp)%ecutwfc )
  END IF
  IF ( PRESENT(ecutrho_pp) ) THEN
     ecutrho_pp = MAXVAL ( upf(1:ntyp)%ecutrho )
  END IF
  !
END SUBROUTINE readpp
!
SUBROUTINE check_order
   ! CP-specific check
   IF ( ANY(upf(1:ntyp)%tpawp) ) CALL errore ('readpp','PAW not implemented',1) 
END SUBROUTINE check_order
!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------+
SUBROUTINE upf_bcast(upf, ionode, ionode_id, comm)
  !---------------------------------------------+
  !
  !! Broadcast the "upf" structure, read on processor "ionode_id",
  !! to all other processors in the communicator "comm".
  !
  USE kinds,        ONLY: DP
  USE pseudo_types, ONLY: pseudo_upf
  USE mp,           ONLY: mp_bcast
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf),INTENT(INOUT) :: upf
  !! pseudo_upf type structure storing the pseudo data
  LOGICAL, INTENT(in) :: ionode
  !! true if we are on the processor that broadcasts
  !! upf is allocated if (ionode), must be allocated otherwise
  INTEGER, INTENT(in) :: ionode_id
  !! ID of the processor that broadcasts
  INTEGER, INTENT(in) :: comm
  !! MPI communicator
  !
  CALL mp_bcast (upf%nv, ionode_id, comm )
  CALL mp_bcast (upf%generated, ionode_id, comm )
  CALL mp_bcast (upf%author, ionode_id, comm )
  CALL mp_bcast (upf%date, ionode_id, comm )
  CALL mp_bcast (upf%comment, ionode_id, comm )
  CALL mp_bcast (upf%psd, ionode_id, comm )
  CALL mp_bcast (upf%typ, ionode_id, comm )
  CALL mp_bcast (upf%rel, ionode_id, comm )
  CALL mp_bcast (upf%tvanp, ionode_id, comm )
  CALL mp_bcast (upf%tpawp, ionode_id, comm )
  CALL mp_bcast (upf%tcoulombp, ionode_id, comm )
  CALL mp_bcast (upf%is_gth, ionode_id, comm )
  CALL mp_bcast (upf%is_multiproj, ionode_id, comm )
  CALL mp_bcast (upf%has_so, ionode_id, comm )
  CALL mp_bcast (upf%has_wfc, ionode_id, comm )
  CALL mp_bcast (upf%has_gipaw, ionode_id, comm )
  CALL mp_bcast (upf%paw_as_gipaw, ionode_id, comm )
  CALL mp_bcast (upf%nlcc, ionode_id, comm )
  CALL mp_bcast (upf%dft, ionode_id, comm )
  CALL mp_bcast (upf%zp, ionode_id, comm )
  CALL mp_bcast (upf%etotps, ionode_id, comm )
  CALL mp_bcast (upf%ecutwfc, ionode_id, comm )
  CALL mp_bcast (upf%ecutrho, ionode_id, comm )
  CALL mp_bcast (upf%lmax, ionode_id, comm )
  CALL mp_bcast (upf%lmax_rho, ionode_id, comm )
  CALL mp_bcast (upf%lloc, ionode_id, comm )
  CALL mp_bcast (upf%mesh, ionode_id, comm )
  CALL mp_bcast (upf%nwfc, ionode_id, comm )
  CALL mp_bcast (upf%nbeta, ionode_id, comm )
  CALL mp_bcast (upf%dx, ionode_id, comm )
  CALL mp_bcast (upf%xmin, ionode_id, comm )
  CALL mp_bcast (upf%rmax, ionode_id, comm )
  CALL mp_bcast (upf%zmesh, ionode_id, comm )
  !
  IF ( .NOT. ionode) ALLOCATE( upf%r( upf%mesh ), upf%rab( upf%mesh ) )
  CALL mp_bcast (upf%r,   ionode_id, comm )
  CALL mp_bcast (upf%rab, ionode_id, comm )
  !
  IF ( .NOT. ionode) ALLOCATE( upf%rho_atc(upf%mesh) )
  CALL mp_bcast (upf%rho_atc, ionode_id, comm )
  !
  IF(.not. upf%tcoulombp) THEN
     IF ( .NOT. ionode) ALLOCATE( upf%vloc(upf%mesh) )
     CALL mp_bcast (upf%vloc, ionode_id, comm )
  ENDIF
  !
  IF ( .not. ionode) THEN
     IF ( upf%nbeta == 0) THEN
        upf%nqf = 0
        upf%nqlc= 0
        upf%qqq_eps= -1._dp
        upf%kkbeta = 0  
        ALLOCATE( upf%kbeta(1),     &
             upf%lll(1),           &
             upf%beta(upf%mesh,1), &
             upf%dion(1,1),        &
             upf%rcut(1),          &
             upf%rcutus(1),        &
             upf%els_beta(1) )
     ELSE
        ALLOCATE( upf%kbeta(upf%nbeta),     &
             upf%lll(upf%nbeta),            &
             upf%beta(upf%mesh, upf%nbeta), &
             upf%dion(upf%nbeta, upf%nbeta),&
             upf%rcut(upf%nbeta),           &
             upf%rcutus(upf%nbeta),         &
             upf%els_beta(upf%nbeta) )
     END IF
  END IF
  !
  CALL mp_bcast (upf%beta, ionode_id, comm )
  CALL mp_bcast (upf%kbeta, ionode_id, comm )
  CALL mp_bcast (upf%els_beta, ionode_id, comm )
  CALL mp_bcast (upf%lll, ionode_id, comm )
  CALL mp_bcast (upf%rcut, ionode_id, comm )
  CALL mp_bcast (upf%rcutus, ionode_id, comm )
  CALL mp_bcast (upf%dion, ionode_id, comm )
  
  IF(upf%tvanp .or. upf%tpawp) THEN
     CALL mp_bcast (upf%q_with_l, ionode_id, comm )
     CALL mp_bcast (upf%nqf, ionode_id, comm )
     CALL mp_bcast (upf%nqlc, ionode_id, comm )
     IF (upf%tpawp) THEN
        IF ( .not. ionode) ALLOCATE &
             ( upf%paw%augmom(upf%nbeta,upf%nbeta, 0:2*upf%lmax) )
        CALL mp_bcast (upf%paw%augshape, ionode_id, comm )
        CALL mp_bcast (upf%paw%raug, ionode_id, comm )
        CALL mp_bcast (upf%paw%iraug, ionode_id, comm )
        CALL mp_bcast (upf%paw%lmax_aug, ionode_id, comm )
        CALL mp_bcast (upf%paw%augmom, ionode_id, comm )
     END IF
     CALL mp_bcast (upf%qqq_eps, ionode_id, comm )
     IF ( .not. ionode) THEN
        IF ( upf%nbeta == 0 ) THEN
           ALLOCATE(upf%rinner(1), &
             upf%qqq(1,1),         &
             upf%qfunc(upf%mesh,1),&
             upf%qfcoef(1,1,1,1) )
             IF ( upf%q_with_l ) &
                ALLOCATE( upf%qfuncl ( upf%mesh, 1, 1 ) )
        ELSE
           ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
           IF ( upf%q_with_l ) THEN
              ALLOCATE( upf%qfuncl ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:2*upf%lmax ) )
           ELSE
              ALLOCATE( upf%qfunc (upf%mesh, upf%nbeta*(upf%nbeta+1)/2) )
           ENDIF
           ALLOCATE( upf%rinner( upf%nqlc ) )
           IF(upf%nqf <= 0) THEN
              ALLOCATE( upf%qfcoef(1,1,1,1) )
           ELSE
              ALLOCATE( upf%qfcoef( upf%nqf, upf%nqlc, &
                   upf%nbeta, upf%nbeta ) )
           END IF
        END IF
     ENDIF
     CALL mp_bcast (upf%qqq   , ionode_id, comm )
     CALL mp_bcast (upf%rinner, ionode_id, comm )
     CALL mp_bcast (upf%qfcoef, ionode_id, comm )
     IF (upf%q_with_l) THEN 
        CALL mp_bcast (upf%qfuncl, ionode_id, comm )
     ELSE
        CALL mp_bcast (upf%qfunc , ionode_id, comm )
     END IF
     !
  END IF
  upf%kkbeta = MAXVAL(upf%kbeta(1:upf%nbeta))
  IF(upf%tpawp) upf%kkbeta = MAX(upf%kkbeta, upf%paw%iraug)
  
  IF ( .not. ionode ) THEN
     ALLOCATE( upf%chi(upf%mesh,upf%nwfc) )
     ALLOCATE( upf%els(upf%nwfc), &
          upf%oc(upf%nwfc), &
          upf%lchi(upf%nwfc), &
          upf%nchi(upf%nwfc), &
          upf%rcut_chi(upf%nwfc), &
          upf%rcutus_chi(upf%nwfc), &
          upf%epseu(upf%nwfc) )
  END IF
  CALL mp_bcast (upf%chi,ionode_id, comm )
  CALL mp_bcast (upf%els, ionode_id, comm )
  CALL mp_bcast (upf%oc,ionode_id, comm )
  CALL mp_bcast (upf%lchi,ionode_id, comm )
  CALL mp_bcast (upf%nchi,ionode_id, comm )
  CALL mp_bcast (upf%rcut_chi,ionode_id, comm )
  CALL mp_bcast (upf%rcutus_chi,ionode_id, comm )
  CALL mp_bcast (upf%epseu,ionode_id, comm )
  !
  IF(upf%has_wfc) THEN
     IF ( .not. ionode) THEN
        ALLOCATE( upf%aewfc(upf%mesh, upf%nbeta) )
        ALLOCATE( upf%pswfc(upf%mesh, upf%nbeta) )
        IF (upf%has_so .and. upf%tpawp) ALLOCATE &
             ( upf%paw%aewfc_rel(upf%mesh, upf%nbeta) )
     END IF
     IF (upf%has_so .and. upf%tpawp) CALL mp_bcast &
          (upf%paw%aewfc_rel,ionode_id,comm )
     CALL mp_bcast &
          (upf%aewfc,ionode_id,comm )
     CALL mp_bcast &
          (upf%pswfc,ionode_id,comm )
  END IF
  !
  IF ( .not. ionode) ALLOCATE( upf%rho_at(upf%mesh) )
  CALL mp_bcast (upf%rho_at,ionode_id,comm )
  
  IF (upf%has_so) THEN
     IF ( .NOT. ionode) THEN
        ALLOCATE (upf%nn(upf%nwfc))
        ALLOCATE (upf%jchi(upf%nwfc))
        ALLOCATE(upf%jjj(upf%nbeta))
     END IF
     CALL mp_bcast (upf%nn,ionode_id,comm )
     CALL mp_bcast (upf%jchi,ionode_id,comm )
     CALL mp_bcast (upf%jjj,ionode_id,comm )
  END IF
  
  IF (upf%tpawp) THEN
     CALL mp_bcast (upf%paw_data_format,ionode_id,comm )
     CALL mp_bcast (upf%paw%core_energy,ionode_id,comm )
     IF ( .not. ionode ) THEN
        ALLOCATE( upf%paw%oc(upf%nbeta) )
        ALLOCATE( upf%paw%ae_rho_atc(upf%mesh) )
        ALLOCATE( upf%paw%ae_vloc(upf%mesh) )
        ALLOCATE( upf%paw%pfunc(upf%mesh, upf%nbeta,upf%nbeta) )
        ALLOCATE(upf%paw%ptfunc(upf%mesh, upf%nbeta,upf%nbeta) )
        IF (upf%has_so) &
             ALLOCATE(upf%paw%pfunc_rel(upf%mesh, upf%nbeta,upf%nbeta) )
     END IF
     CALL mp_bcast (upf%paw%oc,ionode_id,comm )
     CALL mp_bcast (upf%paw%ae_rho_atc,ionode_id,comm )
     CALL mp_bcast (upf%paw%ae_vloc,ionode_id,comm )
     CALL mp_bcast (upf%paw%pfunc,ionode_id,comm )
     CALL mp_bcast (upf%paw%ptfunc,ionode_id,comm )
     IF (upf%has_so) &
          CALL mp_bcast (upf%paw%pfunc_rel,ionode_id,comm )
  END IF
  
  IF (upf%has_gipaw) THEN
     CALL mp_bcast (upf%gipaw_data_format,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_ncore_orbitals,ionode_id,comm )
     IF ( .not. ionode) THEN
        ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
        ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
        ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
        ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
     END IF
     CALL mp_bcast (upf%gipaw_core_orbital_n ,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_core_orbital_el,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_core_orbital_l ,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_core_orbital   ,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_nchannels  ,ionode_id,comm )
     IF ( .not. ionode) THEN
        ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
        ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
        ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
     END IF
     CALL mp_bcast (upf%gipaw_wfs_el, ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_ll, ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_rcut  , ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_rcutus, ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_ae,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_wfs_ps,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_vlocal_ae,ionode_id,comm )
     CALL mp_bcast (upf%gipaw_vlocal_ps,ionode_id,comm )
     !
  END IF
  !
END SUBROUTINE upf_bcast

END MODULE read_pseudo_mod

