!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
      MODULE write_upf_v2_module
!=----------------------------------------------------------------------------=!
!  this module handles the writing of pseudopotential data

! ...   declare modules
#ifdef __STANDALONE
        USE kinds, ONLY: DP, errore
#else
        USE kinds,        ONLY: DP
#endif
        USE pseudo_types, ONLY: pseudo_upf
        USE radial_grids, ONLY: radial_grid_type
        USE iotk_module
        !
        IMPLICIT NONE
        !
        PRIVATE
        PUBLIC :: write_upf_v2
        PUBLIC :: pseudo_upf_meta_info
        PUBLIC :: default_meta_info
        PUBLIC :: allocate_meta_info

      TYPE :: pseudo_upf_meta_info
         ! basic info:
         integer           :: relativistic ! 0=no, 1=scalar, 2=full
         integer           :: lloc ! l of local channel, pseudized AE potential if < 0
         real(DP)          :: rcloc! cutoff radius for local channel
         integer           :: nvo  ! number of valence orbitals
         ! Valence orbitals configuration:
         character(len=2),pointer  :: &
             els(:)           ! label e.g. 3S, 2P ...
         integer,pointer :: &
             nns(:), lls(:)   ! pseudo quantum numbers n and l
         real(DP),pointer :: &
            ocs(:), enls(:), &! occupation and reference energy
            rcut(:), rcutus(:)! norm-conserving and US cutoff radiuses
      END TYPE pseudo_upf_meta_info

 CONTAINS

!------------------------------------------+
SUBROUTINE write_upf_v2(u, upf, meta_info) !
   !---------------------------------------+
   ! Write pseudopotential in UPF format version 2, uses iotk
   !
   IMPLICIT NONE
   INTEGER,INTENT(IN)                      :: u   ! i/o unit
   TYPE(pseudo_upf),INTENT(IN)             :: upf ! the pseudo data
   TYPE(pseudo_upf_meta_info),OPTIONAL,INTENT(IN) :: meta_info ! additional non mandatory info
   !
   CHARACTER(len=iotk_attlenx) :: attr
   INTEGER :: ierr      ! /= 0 if something went wrong
   !
   ! Initialize the file
   CALL iotk_write_attr(attr, 'version', 2, first=.true.)
   CALL iotk_open_write(u, attr=attr, root='UPF', skip_head=.true.)
   !
   ! Write human-readable header
   CALL write_info(u, upf, meta_info)
   !
   ! Write machine-readable header
   CALL write_header(u, upf)
   ! Write radial grid mesh
   CALL write_mesh(u, upf)
   ! Write non-linear core correction charge
   IF(upf%nlcc) CALL iotk_write_dat(u, 'PP_NLCC', upf%rho_atc, columns=4)
   ! Write local potential
   IF(.not. upf%tcoulombp) THEN
      CALL iotk_write_dat(u, 'PP_LOCAL', upf%vloc, columns=4)
   ELSE
         CALL iotk_write_attr(attr, 'type', '1/r', first=.true.)
         CALL iotk_write_attr(attr, 'comment', 'Coulomb 1/r potential')
      CALL iotk_write_empty(u, 'PP_NLCC', attr=attr)
   ENDIF
   ! Write nonlocal components: projectors, augmentation, hamiltonian elements
   CALL write_nonlocal(u, upf)
   ! Write initial pseudo wavefunctions
   ! (usually only wfcs with occupancy > 0)
   CALL write_pswfc(u, upf)
   ! Write valence atomic density (used for initial density)
   CALL iotk_write_dat(u, 'PP_RHOATOM', upf%rho_at, columns=4)
   ! Write additional info for full-relativistic calculation
   CALL write_spin_orb(u, upf)
   ! Write additional data for PAW (All-electron charge, wavefunctions, vloc..)
   CALL write_paw(u, upf)
   ! Write additional data for GIPAW reconstruction
   CALL write_gipaw(u, upf)
   !
   ! Close the file (not the unit!)
   CALL iotk_close_write(u)

   CONTAINS
   !
   SUBROUTINE write_info(u, upf, meta_info)
      ! Write human-readable header
      ! The header is written directly, not via iotk
      IMPLICIT NONE
      INTEGER,INTENT(IN)          :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN) :: upf  ! the pseudo data
      TYPE(pseudo_upf_meta_info),OPTIONAL,TARGET,INTENT(IN) ::&
                                     meta_info ! additional non mandatory info
      !
      INTEGER :: nb ! aux counter
      INTEGER :: ierr ! /= 0 if something went wrong
      TYPE(pseudo_upf_meta_info),POINTER :: mi
      !
      CALL iotk_write_begin(u,'PP_INFO')
      !
      IF(.not. present(meta_info)) THEN
         CALL iotk_write_comment(u,' Meta info was not specified ')
      ELSE
         mi => meta_info
      ENDIF

      WRITE(u, '(4x,a)', err=100) TRIM(upf%generated)
      WRITE(u, '(4x,a)', err=100) &
         'Author: '//TRIM(upf%author)
      WRITE(u, '(4x,a)', err=100) &
         'Generation date: '//TRIM(upf%date)
      WRITE(u, '(4x,a)', err=100) &
         'Pseudopotential type: '//TRIM(upf%typ)
      WRITE(u, '(4x,a)', err=100) &
         'Element: '//TRIM(upf%psd)
      WRITE(u,'()')
      !
      ! Cutoff Information
      WRITE(u, '(4x,a,f5.0,a)') &
         'Suggested minimum cutoff for wavefunctions:',upf%ecutwfc,' Ry'
      WRITE(u, '(4x,a,f5.0,a)') &
         'Suggested minimum cutoff for charge density:',&
         upf%ecutrho,' Ry'
      IF(TRIM(upf%comment) /= ' ') THEN
         WRITE(u, '(4x,a)', err=100) upf%comment
      END IF
      !
      ! Write relativistic information
      IF (mi%relativistic==2) THEN
         WRITE(u, '(4x,a)', err=100) &
               "The Pseudo was generated with a Fully-Relativistic Calculation"
      ELSE IF (mi%relativistic==1) THEN
         WRITE(u, '(4x,a)', err=100) &
               "The Pseudo was generated with a Scalar-Relativistic Calculation"
      ELSE
         WRITE(u, '(4x,a)', err=100) &
               "The Pseudo was generated with a Non-Relativistic Calculation"
      ENDIF
      !
      ! Write local potential information
      IF (mi%lloc >= 0 ) THEN
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "L component and cutoff radius for Local Potential:", mi%lloc, mi%rcloc
      ELSE IF (mi%lloc == -1 ) THEN
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "Local Potential by smoothing AE potential with Bessel fncs, cutoff radius:", mi%rcloc
      ELSE IF (mi%lloc == -2 ) THEN
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "Local Potential according to Troullier-Martins recipe, cutoff radius:", mi%rcloc
      ELSE
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "Local Potential: unknown format, L component and cutoff radius:",mi%lloc, mi%rcloc
      ENDIF
      !
      IF (upf%has_so) &
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "Pseudopotential contains additional information for spin-orbit calculations."
      IF (upf%has_gipaw) &
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
               "Pseudopotential contains additional information for GIPAW reconstruction."

      !
      ! Write valence orbitals information
      WRITE(u, '(4x,a2,2a3,a6,3a19)', err=100) &
            "nl"," pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
      DO nb = 1, mi%nvo
            WRITE(u, '(4x,a2,2i3,f6.2,3f19.11)') mi%els(nb), mi%nns(nb), &
               mi%lls(nb), mi%ocs(nb), mi%rcut(nb), mi%rcutus(nb), mi%enls(nb)
      END DO
      !
      CALL iotk_write_end(u,'PP_INFO')
      CALL iotk_write_comment(u,'                               ')
      CALL iotk_write_comment(u,' END OF HUMAN READABLE SECTION ')
      CALL iotk_write_comment(u,'                               ')
      !
      RETURN
100   CALL errore('write_upf_v2::write_info', 'Writing pseudo file', 1)
   END SUBROUTINE write_info
   !
   SUBROUTINE write_header(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nw
      !
      ! Write HEADER section with some initialization data
         CALL iotk_write_attr(attr, 'version',       upf%nv, first=.true.)
         CALL iotk_write_attr(attr, 'element',       upf%psd)
         CALL iotk_write_attr(attr, 'pseudo_type',   TRIM(upf%typ))
         !
         CALL iotk_write_attr(attr, 'is_ultrasoft',  upf%tvanp)
         CALL iotk_write_attr(attr, 'is_paw',        upf%tpawp)
         CALL iotk_write_attr(attr, 'is_coulomb',    upf%tcoulombp)
         !
         CALL iotk_write_attr(attr, 'has_so',        upf%has_so)
         CALL iotk_write_attr(attr, 'has_gipaw',     upf%has_gipaw)
         !
         CALL iotk_write_attr(attr, 'nlcc',          upf%nlcc)
         CALL iotk_write_attr(attr, 'functional',    upf%dft)
         CALL iotk_write_attr(attr, 'z_valence',     upf%zp)
         CALL iotk_write_attr(attr, 'total_psenergy',upf%etotps)
         CALL iotk_write_attr(attr, 'wfc_cutoff',    upf%ecutwfc)
         CALL iotk_write_attr(attr, 'rho_cutoff',    upf%ecutrho)
         CALL iotk_write_attr(attr, 'l_max',         upf%lmax)
         CALL iotk_write_attr(attr, 'l_max_rho',     upf%lmax_rho)
         CALL iotk_write_attr(attr, 'mesh_size',     upf%mesh)
         CALL iotk_write_attr(attr, 'number_of_wfc', upf%nwfc)
         CALL iotk_write_attr(attr, 'number_of_proj',upf%nbeta)
      CALL iotk_write_begin(u, 'PP_HEADER', attr=attr)
      !
         CALL iotk_write_attr(attr, 'value', TRIM(upf%generated), first=.true.)
         CALL iotk_write_empty(u, 'PP_GENERATED',  attr=attr)
         !
         CALL iotk_write_attr(attr, 'value', TRIM(upf%author), first=.true.)
         CALL iotk_write_empty(u, 'PP_AUTHOR',     attr=attr)
         !
         CALL iotk_write_attr(attr, 'value', TRIM(upf%date), first=.true.)
         CALL iotk_write_empty(u, 'PP_DATE',       attr=attr)
         !
         CALL iotk_write_attr(attr, 'value', TRIM(upf%comment), first=.true.)
         CALL iotk_write_empty(u, 'PP_COMMENT',    attr=attr)
      !
      CALL iotk_write_end(u, 'PP_HEADER')
      !
      RETURN
   END SUBROUTINE write_header
   !
   SUBROUTINE write_mesh(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
         CALL iotk_write_attr(attr, 'dx',   upf%dx, first=.true.)
         CALL iotk_write_attr(attr, 'mesh', upf%mesh)
         CALL iotk_write_attr(attr, 'xmin', upf%xmin)
         CALL iotk_write_attr(attr, 'rmax', upf%rmax)
         CALL iotk_write_attr(attr, 'zmesh',upf%zmesh)
      CALL iotk_write_begin(u, 'PP_MESH', attr=attr)
      !
      CALL iotk_write_dat(u, 'PP_R',   upf%r,   columns=4)
      CALL iotk_write_dat(u, 'PP_RAB', upf%rab, columns=4)
      !
      CALL iotk_write_end(u, 'PP_MESH')
      !
      RETURN
   END SUBROUTINE write_mesh
   !
   SUBROUTINE write_nonlocal(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb,mb,ln,lm,l,nmb
      LOGICAL :: isnull
      !
      IF (upf%tcoulombp) RETURN
      !
      CALL iotk_write_begin(u, 'PP_NONLOCAL')
      !
      ! Write the projectors:
      DO nb = 1,upf%nbeta
            CALL iotk_write_attr(attr, 'index',                  nb, first=.true.)
            CALL iotk_write_attr(attr, 'label',                  upf%els_beta(nb))
            CALL iotk_write_attr(attr, 'angular_momentum',       upf%lll(nb))
            CALL iotk_write_attr(attr, 'cutoff_radius_index',    upf%kbeta(nb))
            CALL iotk_write_attr(attr, 'cutoff_radius',          upf%rcut(nb))
            CALL iotk_write_attr(attr, 'norm_conserving_radius', upf%rcutus(nb))
         CALL iotk_write_dat(u, 'PP_BETA'//iotk_index( nb ), &
                             upf%beta(:,nb), attr=attr, columns=4)
      ENDDO
      !
      ! Write the hamiltonian terms D_ij
         CALL iotk_write_attr(attr, 'non_zero_elements', upf%nd, first=.true.)
      CALL iotk_write_dat(u, 'PP_DIJ', upf%dion, attr=attr, columns=4)
      !
      ! Write the augmentation charge section
      augmentation : &
      IF(upf%tvanp .or. upf%tpawp) THEN
         CALL iotk_write_attr(attr, 'q_with_l', upf%q_with_l, first=.true.)
         CALL iotk_write_attr(attr, 'nqf',      upf%nqf)
         CALL iotk_write_attr(attr, 'nqlc',     upf%nqlc)
         IF (upf%tpawp) THEN
            CALL iotk_write_attr(attr,'shape',          TRIM(upf%paw%augshape))
            CALL iotk_write_attr(attr,'cutoff_r',       upf%paw%raug)
            CALL iotk_write_attr(attr,'cutoff_r_index', upf%paw%iraug)
            CALL iotk_write_attr(attr,'augmentation_epsilon',upf%qqq_eps)
            CALL iotk_write_attr(attr,'lmax_aug',       upf%paw%lmax_aug)
         ENDIF
         !
      CALL iotk_write_begin(u, 'PP_AUGMENTATION', attr=attr)
      !
      ! Write the integrals of the Q functions
      CALL iotk_write_dat(u, 'PP_Q',upf%qqq, columns=4)
      !
      ! Write charge multipoles (only if PAW)
      IF ( upf%tpawp ) THEN
         CALL iotk_write_comment(u, ' augmentation charge multipoles (only for PAW) ')
         CALL iotk_write_dat(u, 'PP_MULTIPOLES', &
                             upf%paw%augmom, attr=attr, columns=4)
      ENDIF
      !
      ! Write polinomial coefficients for Q_ij expansion at small radius
      IF ( upf%nqf > 0) THEN
         CALL iotk_write_comment(u, ' polinomial expansion of Q_ij at small radius ')
         CALL iotk_write_dat(u, 'PP_QFCOEF',upf%qfcoef, attr=attr, columns=4)
         CALL iotk_write_dat(u, 'PP_RINNER',upf%rinner, attr=attr, columns=4)
      ENDIF
      !
      ! Write augmentation charge Q_ij
      DO nb = 1,upf%nbeta
      ln = upf%lll(nb)
      DO mb = nb,upf%nbeta
      lm = upf%lll(mb)
      nmb = mb * (mb-1) /2 + nb
         IF( upf%q_with_l ) THEN
            DO l = abs(ln-lm),ln+lm,2 ! only even terms
               CALL iotk_write_attr(attr, 'first_index',  nb, first=.true.)
               CALL iotk_write_attr(attr, 'second_index', mb)
               CALL iotk_write_attr(attr, 'composite_index', nmb)
               CALL iotk_write_attr(attr, 'angular_momentum', l)
               !
               isnull = .false. ! omit functions that are multiplied by zero
               IF( upf%tpawp ) isnull = (abs(upf%paw%augmom(nb,mb,l)) < upf%qqq_eps)
               !
               IF ( isnull ) THEN
                  CALL iotk_write_attr(attr, 'is_null', isnull)
                  CALL iotk_write_empty(u, 'PP_QIJL'//iotk_index((/nb,mb,l/)),&
                                       attr=attr)
               ELSE
                  CALL iotk_write_dat(u, 'PP_QIJL'//iotk_index((/nb,mb,l/)),&
                                    upf%qfuncl(:,nmb,l),attr=attr, columns=4)
               ENDIF
            ENDDO
         ELSE
            CALL iotk_write_attr(attr, 'first_index',  nb, first=.true.)
            CALL iotk_write_attr(attr, 'second_index', mb)
            CALL iotk_write_attr(attr, 'composite_index', nmb)
            !
            isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
            IF ( isnull ) THEN
               CALL iotk_write_attr(attr, 'is_null', isnull)
               CALL iotk_write_empty(u, 'PP_QIJ'//iotk_index((/nb,mb/)),&
                                 attr=attr)
            ELSE
               CALL iotk_write_dat(u, 'PP_QIJ'//iotk_index((/nb,mb/)),&
                                 upf%qfunc(:,nmb),attr=attr, columns=4)
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      !
      CALL iotk_write_end(u, 'PP_AUGMENTATION')
      !
      ENDIF augmentation
      !
      CALL iotk_write_end(u, 'PP_NONLOCAL')
      !
      RETURN
   END SUBROUTINE write_nonlocal
   !
   SUBROUTINE write_pswfc(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nw
      !
      CALL iotk_write_begin(u, 'PP_PSWFC')
      !
      DO nw = 1,upf%nwfc
         CALL iotk_write_attr(attr, 'index',      nw, first=.true.)
         CALL iotk_write_attr(attr, 'label',      upf%els(nw))
         CALL iotk_write_attr(attr, 'l',          upf%lchi(nw))
         CALL iotk_write_attr(attr, 'occupation', upf%oc(nw))
         CALL iotk_write_dat(u, 'PP_CHI'//iotk_index(nw), &
                              upf%chi(:,nw), columns=4, attr=attr)
      ENDDO
      !
      CALL iotk_write_end(u, 'PP_PSWFC')
      !
      RETURN
   END SUBROUTINE write_pswfc
   !
   SUBROUTINE write_spin_orb(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nw, nb
      !
      IF (.not. upf%has_so) RETURN
      !
      CALL iotk_write_begin(u, 'PP_SPIN_ORB')
      !
      DO nw = 1,upf%nwfc
            CALL iotk_write_attr(attr, 'index', nw, first=.true.)
            CALL iotk_write_attr(attr, 'els',   upf%els(nw))
            CALL iotk_write_attr(attr, 'nn',    upf%nn(nw))
            CALL iotk_write_attr(attr, 'lchi',  upf%lchi(nw))
            CALL iotk_write_attr(attr, 'jchi',  upf%jchi(nw))
            CALL iotk_write_attr(attr, 'oc',    upf%oc(nw))
         CALL iotk_write_empty(u, 'PP_RELWFC'//iotk_index(nw),&
                               attr=attr)
      ENDDO
      !
      DO nb = 1,upf%nbeta
            CALL iotk_write_attr(attr, 'index', nb, first=.true.)
            CALL iotk_write_attr(attr, 'lll',   upf%lll(nb))
            CALL iotk_write_attr(attr, 'jjj',   upf%jjj(nb))
         CALL iotk_write_empty(u, 'PP_RELBETA'//iotk_index(nb),&
                               attr=attr)
      ENDDO
      !
      CALL iotk_write_end(u, 'PP_SPIN_ORB')
      !
      RETURN
   END SUBROUTINE write_spin_orb
   !
   SUBROUTINE write_paw(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb

      IF (.not. upf%tpawp ) RETURN

         CALL iotk_write_attr(attr, 'paw_data_format', upf%paw_data_format, first=.true.)
         CALL iotk_write_attr(attr, 'core_energy',     upf%paw%core_energy)
      CALL iotk_write_begin(u, 'PP_PAW', attr=attr)
      ! Full occupation (not only > 0 ones)
      CALL iotk_write_dat(u, 'PP_OCCUPATIONS',upf%paw%oc, columns=4)
      ! All-electron core charge
      CALL iotk_write_dat(u, 'PP_AE_NLCC', upf%paw%ae_rho_atc, columns=4)
      ! All-electron local potential
      CALL iotk_write_dat(u, 'PP_AE_VLOC', upf%paw%ae_vloc,columns=4)
      ! All-electron wavefunctions
      DO nb = 1,upf%nbeta
         CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
         CALL iotk_write_attr(attr, 'label',      upf%els_beta(nb))
         CALL iotk_write_attr(attr, 'l',          upf%lll(nb))
         CALL iotk_write_attr(attr, 'occupation', upf%paw%oc(nb))
         CALL iotk_write_dat(u, 'PP_AEWFC'//iotk_index(nb), &
                              upf%aewfc(:,nb), columns=4, attr=attr)
      ENDDO
      ! Pseudo wavefunctions (not only the ones for oc > 0)
      DO nb = 1,upf%nbeta
         CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
         CALL iotk_write_attr(attr, 'label',      upf%els_beta(nb))
         CALL iotk_write_attr(attr, 'l',          upf%lll(nb))
         CALL iotk_write_attr(attr, 'occupation', upf%paw%oc(nb))
         CALL iotk_write_dat(u, 'PP_PSWFC'//iotk_index(nb), &
                              upf%pswfc(:,nb), columns=4, attr=attr)
      ENDDO
      ! Finalize
      CALL iotk_write_end(u, 'PP_PAW')

      RETURN
   END SUBROUTINE write_paw
!
   SUBROUTINE write_gipaw(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb
      IF (.not. upf%has_gipaw ) RETURN

         CALL iotk_write_attr(attr, 'gipaw_data_format', upf%gipaw_data_format, first=.true.)
      CALL iotk_write_begin(u, 'PP_GIPAW', attr=attr)

         CALL iotk_write_attr(attr, 'number_of_core_orbitals', upf%gipaw_ncore_orbitals, first=.true.)
      CALL iotk_write_begin(u, 'PP_GIPAW_CORE_ORBITALS', attr=attr)
      DO nb = 1,upf%gipaw_ncore_orbitals
            CALL iotk_write_attr(attr, 'index', nb, first=.true.)
            CALL iotk_write_attr(attr, 'label', upf%gipaw_core_orbital_el(nb))
            CALL iotk_write_attr(attr, 'n',     upf%gipaw_core_orbital_n(nb))
            CALL iotk_write_attr(attr, 'l',     upf%gipaw_core_orbital_l(nb))
         CALL iotk_write_dat(u, 'PP_GIPAW_CORE_ORBITAL'//iotk_index(nb), &
                              upf%gipaw_core_orbital(:,nb), columns=4, attr=attr)
      ENDDO
      CALL iotk_write_end(u, 'PP_GIPAW_CORE_ORBITALS')
      !
      ! Write valence all-electron and pseudo orbitals
         CALL iotk_write_attr(attr, 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels, first=.true.)
      CALL iotk_write_begin(u, 'PP_GIPAW_ORBITALS', attr=attr)
      DO nb = 1,upf%gipaw_wfs_nchannels
            CALL iotk_write_attr(attr, 'index', nb, first=.true.)
            CALL iotk_write_attr(attr, 'label', upf%gipaw_wfs_el(nb))
            CALL iotk_write_attr(attr, 'l',     upf%gipaw_wfs_ll(nb))
            CALL iotk_write_attr(attr, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb))
            CALL iotk_write_attr(attr, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb))
         CALL iotk_write_begin(u, 'PP_GIPAW_CORE_ORBITAL'//iotk_index(nb), attr=attr)
         !
         CALL iotk_write_dat(u, 'PP_GIPAW_WFS_AE', upf%gipaw_wfs_ae(:,nb), columns=4)
         CALL iotk_write_dat(u, 'PP_GIPAW_WFS_PS', upf%gipaw_wfs_ps(:,nb), columns=4)
         !
         CALL iotk_write_end(u, 'PP_GIPAW_ORBITAL'//iotk_index(nb))
      ENDDO
      CALL iotk_write_end(u, 'PP_GIPAW_ORBITALS')
      !
      ! Write all-electron and pseudo local potentials
      CALL iotk_write_begin(u, 'PP_GIPAW_VLOCAL')
      CALL iotk_write_dat(u, 'PP_GIPAW_VLOCAL_AE'//iotk_index(nb), &
                           upf%gipaw_vlocal_ae(:), columns=4)
      CALL iotk_write_dat(u, 'PP_GIPAW_VLOCAL_PS'//iotk_index(nb), &
                           upf%gipaw_vlocal_ae(:), columns=4)
      CALL iotk_write_end(u, 'PP_GIPAW_VLOCAL')
      !
      CALL iotk_write_end(u, 'PP_GIPAW')

      RETURN
   END SUBROUTINE write_gipaw
!
END SUBROUTINE write_upf_v2
!##############################################################################
!##############################################################################
SUBROUTINE default_meta_info(mi, upf)
   TYPE(pseudo_upf_meta_info),INTENT(INOUT) :: mi
   TYPE(pseudo_upf),OPTIONAL,INTENT(IN) :: upf
   INTEGER :: nw
         mi%relativistic=0
         mi%lloc= -9
         mi%rcloc=0.0_dp
         IF ( .not. present(upf)) THEN
            mi%nvo = 0
            CALL allocate_meta_info(mi,0)
         ELSE
            mi%nvo = size(upf%oc)
            CALL allocate_meta_info(mi,mi%nvo)
            mi%ocs(:) = upf%oc(:)
            !
            IF(associated(upf%els)) THEN
               mi%els(:) = upf%els(:)
            ELSE
               mi%els(:) = 'nX'
            ENDIF
            IF(associated(upf%lchi)) THEN
               mi%lls(:) = upf%lchi(:)
               mi%nns(:) = upf%lchi(:)+1
            ELSE
               mi%lls(:) = -1
               mi%nns(:) = -1
            ENDIF
            IF(associated(upf%epseu) THEN
               mi%enls(:)= upf%epseu(:)
            ELSE
               mi%enls(:)= 0.0_dp
            ENDIF
            IF(associated(upf%rcut)) THEN
               mi%rcut(:)= upf%rcut(:)
            ELSE
               mi%rcut(:)= 0.0_dp
            ENDIF
            IF(associated(upf%rcutus)) THEN
               mi%rcutus(:)=upf%rcutus(:)
            ELSE
               mi%rcutus(:)=0.0_dp
            ENDIF
         ENDIF
END SUBROUTINE default_meta_info
SUBROUTINE allocate_meta_info(mi,nvo)
   INTEGER,INTENT(IN) :: nvo
   TYPE(pseudo_upf_meta_info),INTENT(INOUT) :: mi
         allocate(              &
            mi%els(mi%nvo),   &
            mi%nns(mi%nvo),   &
            mi%lls(mi%nvo),   &
            mi%ocs(mi%nvo),   &
            mi%enls(mi%nvo),  &
            mi%rcut(mi%nvo),  &
            mi%rcutus(mi%nvo) &
         )
END SUBROUTINE allocate_meta_info

END MODULE write_upf_v2_module
