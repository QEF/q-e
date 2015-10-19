!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE read_pseudo_mod
!=----------------------------------------------------------------------------=!
  !
  ! read pseudopotential files. Note that all processors read the same file!
  !
  ! Required on input:
  USE io_files,     ONLY: pseudo_dir, pseudo_dir_cur, psfile
  USE ions_base,    ONLY: ntyp => nsp
  ! Modified on output:
  USE atom,         ONLY: msh, rgrid
  USE ions_base,    ONLY: zv
  USE uspp_param,   ONLY: upf, newpseudo, oldvan, nvb
  USE uspp,         ONLY: okvan, nlcc_any

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
  ! Reads PP files and puts the result into the "upf" structure
  ! Sets  DFT to input_dft if present, to the value read in PP files otherwise
  ! Sets  number of valence electrons Zv, control variables okvan and nlcc_any,
  !       compatibility variables newpseudo, oldvan, nvb
  ! Optionally returns cutoffs read from PP files into ecutwfc_pp, ecutrho_pp
  !
  USE kinds,        ONLY: DP
  USE mp,           ONLY: mp_bcast, mp_sum
  USE mp_images,    ONLY: intra_image_comm
  USE io_global,    ONLY: stdout, ionode
  USE pseudo_types, ONLY: pseudo_upf, nullify_pseudo_upf, deallocate_pseudo_upf
  USE funct,        ONLY: enforce_input_dft, &
                          get_iexch, get_icorr, get_igcx, get_igcc, get_inlc
  use radial_grids, ONLY: deallocate_radial_grid, nullify_radial_grid
  USE wrappers,     ONLY: md5_from_file
  USE upf_module,   ONLY: read_upf
  USE upf_to_internal,  ONLY: set_pseudo_upf
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
  CHARACTER(len=256) :: file_pseudo ! file name complete with path
  LOGICAL :: printout_ = .FALSE.
  INTEGER :: iunps, isupf, nt, nb, ir, ios
  INTEGER :: iexch_, icorr_, igcx_, igcc_, inlc_
  !
  ! ... initialization: allocate radial grids etc
  !
  iunps = 4
  IF( ALLOCATED( rgrid ) ) THEN
     DO nt = 1, SIZE( rgrid )
        CALL deallocate_radial_grid( rgrid( nt ) )
        CALL nullify_radial_grid( rgrid( nt ) )
     END DO
     DEALLOCATE( rgrid )
     if(allocated(msh)) DEALLOCATE( msh )
  END IF

  ALLOCATE( rgrid( ntyp ), msh( ntyp ) )

  DO nt = 1, ntyp
     CALL nullify_radial_grid( rgrid( nt ) )
  END DO

  IF( ALLOCATED( upf ) ) THEN
     DO nt = 1, SIZE( upf )
        CALL deallocate_pseudo_upf( upf( nt ) )
        CALL nullify_pseudo_upf( upf( nt ) )
     END DO
     DEALLOCATE( upf )
  END IF
  !
  ALLOCATE ( upf( ntyp ) )
  !
  !  nullify upf objects as soon as they are instantiated
  !
  do nt = 1, ntyp 
     CALL nullify_pseudo_upf( upf( nt ) )
  end do
  !
  IF (input_dft /='none') CALL enforce_input_dft (input_dft)
  !
  IF ( PRESENT(printout) ) THEN
     printout_ = printout .AND. ionode
  END IF
  IF ( printout_) THEN
     WRITE( stdout,"(//,3X,'Atomic Pseudopotentials Parameters',/, &
                   &    3X,'----------------------------------' )" )
  END IF
  !
  nvb = 0
  do nt = 1, ntyp
     !
     ! variables not necessary for USPP, but necessary for PAW;
     ! will be read from file if it is a PAW dataset.
     !
     rgrid(nt)%xmin = 0.d0
     rgrid(nt)%dx = 0.d0
     !
     ! try first pseudo_dir_cur if set: in case of restart from file,
     ! this is where PP files should be located
     !
     ios = 1
     IF ( pseudo_dir_cur /= ' ' ) THEN
        file_pseudo  = TRIM (pseudo_dir_cur) // TRIM (psfile(nt))
        OPEN  (unit = iunps, file = file_pseudo, status = 'old', &
               form = 'formatted', action='read', iostat = ios)
        CALL mp_sum (ios,intra_image_comm)
        IF ( ios /= 0 ) CALL infomsg &
                     ('readpp', 'file '//TRIM(file_pseudo)//' not found')
        !
        ! file not found? no panic (yet): if the restart file is not visible
        ! to all processors, this may happen. Try the original location
     END IF
     !
     ! try the original location pseudo_dir, as set in input
     ! (it should already contain a slash at the end)
     !
     IF ( ios /= 0 ) THEN
        file_pseudo = TRIM (pseudo_dir) // TRIM (psfile(nt))
        OPEN  (unit = iunps, file = file_pseudo, status = 'old', &
               form = 'formatted', action='read', iostat = ios)
        CALL mp_sum (ios,intra_image_comm)
        CALL errore('readpp', 'file '//TRIM(file_pseudo)//' not found',ABS(ios))
     END IF
     !
     upf(nt)%grid => rgrid(nt)
     !
     ! start reading - UPF first: the UPF format is detected via the
     ! presence of the keyword '<PP_HEADER>' at the beginning of the file
     !
     IF( printout_ ) THEN
        WRITE( stdout, "(/,3X,'Reading pseudopotential for specie # ',I2, &
                       & ' from file :',/,3X,A)") nt, TRIM(file_pseudo)
     END IF
     !
     call read_upf(upf(nt), rgrid(nt), isupf, unit=iunps)
     !
     upf(nt)%is_gth=.false.
     if (isupf ==-1 .OR. isupf== 0) then
        !
        IF( printout_ ) &
           WRITE( stdout, "(3X,'file type is UPF v.',i1)") isupf+2
        call set_pseudo_upf (nt, upf(nt))
        ! 
        ! UPF is assumed to be multi-projector
        !
        newpseudo (nt) = .true.
        !
     else
        !
        rewind (unit = iunps)
        !
        !     The type of the pseudopotential is determined by the file name:
        !    *.vdb or *.van  Vanderbilt US pseudopotential code  pseudo_type=1
        !    *.RRKJ3         Andrea's   US new code              pseudo_type=2
        !    *.gth           Goedecker-Teter-Hutter NC pseudo    pseudo_type=3
        !    none of the above: PWSCF norm-conserving format     pseudo_type=0
        !
        if ( pseudo_type (psfile (nt) ) == 1 .or. &
             pseudo_type (psfile (nt) ) == 2 ) then
           !
           ! PPs produced by Andrea Dal Corso's atomic code are assumed to
           ! be multiprojector; NCPP produced by Vanderbilt's core are not
           !    
           newpseudo (nt) = ( pseudo_type (psfile (nt) ) == 2 )
           !
           IF ( newpseudo (nt) ) THEN
              IF( printout_ ) &
                 WRITE( stdout, "(3X,'file type is RRKJ3')")
              call readrrkj (iunps, nt, upf(nt))
           ELSE
              IF( printout_ ) &
                 WRITE( stdout, "(3X,'file type is Vanderbilt US PP')")
              CALL readvan (iunps, nt, upf(nt))
           ENDIF
           CALL set_pseudo_upf (nt, upf(nt), rgrid(nt))
           !
        elseif ( pseudo_type (psfile (nt) ) == 3 ) then
           newpseudo (nt) = .true.
           !
           CALL readgth (iunps, nt, upf(nt))
           !
           CALL set_pseudo_upf (nt, upf(nt), rgrid(nt))
           !
        else
           newpseudo (nt) = .false.
           IF( printout_ ) &
              WRITE( stdout, "(3X,'file type is old PWscf NC format')")
           ! 
           call read_ncpp (iunps, nt, upf(nt))
           !
           CALL set_pseudo_upf (nt, upf(nt), rgrid(nt)) 
           !
        endif
        !
     endif
     !
     ! end of reading
     !
     close (iunps)
     !
     ! Calculate MD5 checksum for this pseudopotential
     !
     CALL md5_from_file(file_pseudo, upf(nt)%md5_cksum)
     !
     ! ... Zv = valence charge of the (pseudo-)atom, read from PP files,
     ! ... is set equal to Zp = pseudo-charge of the pseudopotential
     !
     zv(nt) = upf(nt)%zp
     !
     ! ... count US species
     !
     IF (upf(nt)%tvanp) nvb=nvb+1
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
     ! the radial grid is defined up to r(mesh) but we introduce 
     ! an auxiliary variable msh to limit the grid up to rcut=10 a.u. 
     ! This is used to cut off the numerical noise arising from the
     ! large-r tail in cases like the integration of V_loc-Z/r
     !
     do ir = 1, rgrid(nt)%mesh
        if (rgrid(nt)%r(ir) > rcut) then
           msh (nt) = ir
           goto 5
        endif
     enddo
     msh (nt) = rgrid(nt)%mesh 
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1
     !
     ! msh is forced to be odd for simpson integration (maybe obsolete?)
     !
     ! check for zero atomic wfc, 
     ! check that (occupied) atomic wfc are properly normalized
     !
     call check_atwfc_norm(nt)
     !
  enddo
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
  return
end subroutine readpp
!-----------------------------------------------------------------------
integer function pseudo_type (psfile)
  !-----------------------------------------------------------------------
  implicit none
  character (len=*) :: psfile
  integer :: l
  !
  l = len_trim (psfile)
  pseudo_type = 0
  if (psfile (l - 3:l) .eq.'.vdb'.or.psfile (l - 3:l) .eq.'.van') &
       pseudo_type = 1
  if (psfile (l - 3:l) .eq.'.gth') pseudo_type = 3
  if (l > 5) then
     if (psfile (l - 5:l) .eq.'.RRKJ3') pseudo_type = 2
  end if
  !
  return

end function pseudo_type

!---------------------------------------------------------------
SUBROUTINE check_atwfc_norm(nt)
  !---------------------------------------------------------------
  !  check for the presence of zero wavefunctions first
  !  check the normalization of the atomic wfc (only those with non-negative
  !  occupations) and renormalize them if the calculated norm is incorrect 
  !  by more than eps6 (10^{-6})
  !
  USE kinds,        ONLY : dp
  USE constants,    ONLY : eps6, eps8
  USE io_global,    ONLY : stdout

  implicit none

  integer,intent(in) :: nt ! index of the pseudopotential to be checked
  !
  integer ::             &
     mesh, kkbeta,       & ! auxiliary indices of integration limits
     l,                  & ! orbital angular momentum 
     iwfc, ir,           & ! counter on atomic wfcs and on radial mesh
     ibeta, ibeta1, ibeta2 ! counters on betas
  logical :: &
     match                 ! a logical variable 
  real(DP) :: &
     norm,               & ! the norm
     j                     ! total (spin+orbital) angular momentum
  real(DP), allocatable :: &
     work(:), gi(:)        ! auxiliary variable for becp
  character (len=80) :: renorm
  !
  allocate (work(upf(nt)%nbeta), gi(upf(nt)%grid%mesh) )

  ! define indices for integration limits
  mesh = upf(nt)%grid%mesh
  kkbeta = upf(nt)%kkbeta
  !
  renorm = ' '
  DO iwfc = 1, upf(nt)%nwfc
     l = upf(nt)%lchi(iwfc)
     if ( upf(nt)%has_so ) j = upf(nt)%jchi(iwfc)
     !
     ! the smooth part first ..
     gi(1:mesh) = upf(nt)%chi(1:mesh,iwfc) * upf(nt)%chi(1:mesh,iwfc)
     call simpson (mesh, gi, upf(nt)%grid%rab, norm)
     !
     IF ( norm < eps8 ) then
        WRITE( stdout,'(5X,"WARNING: atomic wfc # ",i2, &
             & " for atom type",i2," has zero norm")') iwfc, nt
       !
       ! set occupancy to a small negative number so that this wfc
       ! is not going to be used for starting wavefunctions
       !
       upf(nt)%oc (iwfc) = -eps8
     END IF
     !
     IF ( upf(nt)%oc(iwfc) < 0.d0) CYCLE ! only occupied states are normalized
     !
     if (  upf(nt)%tvanp ) then
        !
        ! the US part if needed
        do ibeta = 1, upf(nt)%nbeta
           match = l.eq.upf(nt)%lll(ibeta)
           if (upf(nt)%has_so) match=match.and.abs(j-upf(nt)%jjj(ibeta)) < eps6
           if (match) then
              gi(1:kkbeta)= upf(nt)%beta(1:kkbeta,ibeta) * &
                            upf(nt)%chi (1:kkbeta,iwfc) 
              call simpson (kkbeta, gi, upf(nt)%grid%rab, work(ibeta))
           else
              work(ibeta)=0.0_dp
           endif
        enddo
        do ibeta1=1,upf(nt)%nbeta
           do ibeta2=1,upf(nt)%nbeta
              norm=norm+upf(nt)%qqq(ibeta1,ibeta2)*work(ibeta1)*work(ibeta2)  
           enddo
        enddo
     end if
     norm=sqrt(norm)
     if (abs(norm-1.0_dp) > eps6 ) then
        renorm = TRIM(renorm) // ' ' // upf(nt)%els(iwfc)
        upf(nt)%chi(1:mesh,iwfc)=upf(nt)%chi(1:mesh,iwfc)/norm
     end if
  end do
  deallocate (work, gi )
  IF ( LEN_TRIM(renorm) > 0 ) WRITE( stdout, &
     '(15x,"file ",a,": wavefunction(s) ",a," renormalized")') &
     TRIM(psfile(nt)),TRIM(renorm)
  RETURN
  !
END SUBROUTINE check_atwfc_norm

SUBROUTINE check_order
   ! CP-specific check
   INTEGER :: nt
   IF ( ANY(upf(1:ntyp)%tpawp) ) CALL errore ('readpp','PAW not implemented',1) 
   DO nt =2, ntyp
      IF ( (.NOT. upf(nt-1)%tvanp) .AND. upf(nt)%tvanp ) THEN
        CALL errore ('readpp', 'ultrasoft PPs must precede norm-conserving',nt)
      ENDIF
   END DO
END SUBROUTINE check_order
END MODULE read_pseudo_mod
