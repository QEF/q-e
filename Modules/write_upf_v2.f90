!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
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
  USE kinds,        ONLY: DP
  USE pseudo_types, ONLY: pseudo_upf
  USE radial_grids, ONLY: radial_grid_type
  USE iotk_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_upf_v2, pseudo_config, deallocate_pseudo_config

  TYPE pseudo_config
     INTEGER :: nwfs
     CHARACTER(len=32) :: pseud
     CHARACTER(len=2),POINTER :: els(:) !=> null()    ! label
     INTEGER,POINTER          :: nns(:) !=> null()    ! n
     INTEGER,POINTER          :: lls(:) !=> null()    ! l
     REAL(DP),POINTER         :: ocs(:) !=> null()    ! occupation
     REAL(DP),POINTER         :: rcut(:) !=> null()   ! NC cutoff radius
     REAL(DP),POINTER         :: rcutus(:) !=> null() ! US cutoff radius
     REAL(DP),POINTER         :: enls(:) !=> null()   ! energy
  END TYPE pseudo_config

CONTAINS

  !-------------------------------+
  SUBROUTINE write_upf_v2(u, upf, conf, u_input)
    !----------------------------+
    ! Write pseudopotential in UPF format version 2, uses iotk
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: u   ! unit for writing
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    ! optional: configuration used to generate the pseudopotential
    TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
    ! optional: unit pointing to input file containing generation data
    INTEGER, OPTIONAL, INTENT(IN) :: u_input
    !
    CHARACTER(len=iotk_attlenx) :: attr
    !
    ! Initialize the file
    CALL iotk_write_attr(attr, 'version', TRIM(upf%nv), first=.true.)
    CALL iotk_open_write(u, attr=attr, root='UPF', skip_head=.true.)
    !
    ! Write human-readable header
    CALL write_info(u, upf, conf, u_input)
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
    ! Write potentials in semilocal form (if existing)
    IF ( upf%typ == "SL" ) CALL write_semilocal(u, upf)
    ! Write nonlocal components: projectors, augmentation, hamiltonian elements
    CALL write_nonlocal(u, upf)
    ! Write initial pseudo wavefunctions
    ! (usually only wfcs with occupancy > 0)
    CALL write_pswfc(u, upf)
    ! If included, write all-electron and pseudo wavefunctions
    CALL write_full_wfc(u, upf)
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
    SUBROUTINE write_info(u, upf, conf, u_input)
      ! Write human-readable header
      ! The header is written directly, not via iotk
      IMPLICIT NONE
      INTEGER,INTENT(IN)          :: u    ! i/o unit: write to unit u
      TYPE(pseudo_upf),INTENT(IN) :: upf  ! the pseudo data
      ! optional: configuration used to generate the pseudopotential
      TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
      INTEGER, OPTIONAL, INTENT(IN) :: u_input ! read input data from u_input
      !
      INTEGER :: nb ! aux counter
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(len=256) :: line
      LOGICAL :: opnd
      !
      CALL iotk_write_begin(u,'PP_INFO')
      !
      WRITE(u, '(4x,a)', err=100) TRIM(CHECK(upf%generated))
      WRITE(u, '(4x,a)', err=100) &
           'Author: '//TRIM(CHECK(upf%author))
      WRITE(u, '(4x,a)', err=100) &
           'Generation date: '//TRIM(CHECK(upf%date))
      WRITE(u, '(4x,a)', err=100) &
           'Pseudopotential type: '//TRIM(CHECK(upf%typ))
      WRITE(u, '(4x,a)', err=100) &
           'Element: '//TRIM(CHECK(upf%psd))
      WRITE(u, '(4x,a)', err=100) &
           'Functional: '//TRIM(CHECK(upf%dft))
      WRITE(u,'()')
      !
      ! Cutoff Information
      WRITE(u, '(4x,a,f5.0,a)') &
           'Suggested minimum cutoff for wavefunctions:',upf%ecutwfc,' Ry'
      WRITE(u, '(4x,a,f5.0,a)') &
           'Suggested minimum cutoff for charge density:',&
           upf%ecutrho,' Ry'
      !
      ! Write relativistic information
      IF (TRIM(upf%rel)=='full') THEN
         WRITE(u, '(4x,a)', err=100) &
              "The Pseudo was generated with a Fully-Relativistic Calculation"
      ELSE IF (TRIM(upf%rel)=='scalar') THEN
         WRITE(u, '(4x,a)', err=100) &
              "The Pseudo was generated with a Scalar-Relativistic Calculation"
      ELSE
         WRITE(u, '(4x,a)', err=100) &
              "The Pseudo was generated with a Non-Relativistic Calculation"
      ENDIF
      !
      ! Write local potential information
      IF (upf%lloc >= 0 ) THEN
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
              "L component and cutoff radius for Local Potential:", upf%lloc, upf%rcloc
      ELSE IF (upf%lloc == -1 ) THEN
         WRITE(u, '(4x,a,f9.4)', err=100) &
              "Local Potential by smoothing AE potential with Bessel fncs, cutoff radius:", upf%rcloc
      ELSE IF (upf%lloc == -2 ) THEN
         WRITE(u, '(4x,a,f9.4)', err=100) &
              "Local Potential according to Troullier-Martins recipe, cutoff radius:", upf%rcloc
      ELSE
         WRITE(u, '(4x,a,i3,f9.4)', err=100) &
              "Local Potential: unknown format, L component and cutoff radius:",upf%lloc, upf%rcloc
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
      WRITE(u, '(/,4x,a)') 'Valence configuration: '
      WRITE(u, '(4x,a2,2a3,a6,2a11,1a13)', err=100) &
           "nl"," pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
      DO nb = 1, upf%nwfc
         IF(upf%oc(nb) >= 0._dp) THEN
            WRITE(u, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
                 CHECK(upf%els(nb)), upf%nchi(nb), &
                 upf%lchi(nb), upf%oc(nb), upf%rcut_chi(nb), &
                 upf%rcutus_chi(nb), upf%epseu(nb)
         ENDIF
      END DO
      IF( present(conf) ) THEN
         WRITE(u, '(4x,a)') 'Generation configuration:'
         DO nb = 1,conf%nwfs
            WRITE(u, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
                 CHECK(conf%els(nb)), conf%nns(nb), &
                 conf%lls(nb), conf%ocs(nb), conf%rcut(nb), &
                 conf%rcutus(nb), conf%enls(nb)
         ENDDO
         WRITE(u,'(/,4x,2a)') 'Pseudization used: ',TRIM(CHECK(conf%pseud))
      ELSE
         WRITE(u, '(/,4x,a)') 'Generation configuration: not available.'
      ENDIF

      IF(TRIM(upf%comment) /= ' ') THEN
         WRITE(u, '(4x,"Comment:",/,4x,a)', err=100) TRIM(CHECK(upf%comment))
      END IF
      !
      IF ( PRESENT(u_input) ) THEN
         !
         ! copy content of input file used in pseudopotential generation
         !
         INQUIRE (unit=u_input, opened=opnd)
         IF (opnd) THEN
            WRITE (u,'("<PP_INPUTFILE>")')
            REWIND (unit=u_input)
10          READ (u_input, '(A)',end=20,err=25) line
            WRITE (u, '(A)') TRIM(CHECK(line))
            GO TO 10
25          CALL infomsg('write_upf_v2::write_inputfile', 'problem writing input data')
20          WRITE (u,'("</PP_INPUTFILE>")')
         ELSE
            CALL infomsg('write_upf_v2::write_inputfile', 'input file not open')
         END IF
         !
      END IF
      !
      CALL iotk_write_end(u,'PP_INFO')
      CALL iotk_write_comment(u,'                               ')
      CALL iotk_write_comment(u,' END OF HUMAN READABLE SECTION ')
      CALL iotk_write_comment(u,'                               ')
      !
      RETURN
100   CALL errore('write_upf_v2::write_info', 'Writing pseudo file', 1)
      !
    END SUBROUTINE write_info
    !
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
      !CALL iotk_write_attr(attr, 'version',       upf%nv, first=.true., newline=.true.)
      CALL iotk_write_attr(attr, 'generated',      TRIM(upf%generated),first=.true.)
      CALL iotk_write_attr(attr, 'author',         TRIM(upf%author),   newline=.true.)
      CALL iotk_write_attr(attr, 'date',           TRIM(upf%date),     newline=.true.)
      CALL iotk_write_attr(attr, 'comment',        TRIM(upf%comment),  newline=.true.)
      !
      CALL iotk_write_attr(attr, 'element',        upf%psd,        newline=.true.)
      CALL iotk_write_attr(attr, 'pseudo_type',    TRIM(upf%typ),  newline=.true.)
      CALL iotk_write_attr(attr, 'relativistic',   TRIM(upf%rel),  newline=.true.)
      !
      CALL iotk_write_attr(attr, 'is_ultrasoft',   upf%tvanp,      newline=.true.)
      CALL iotk_write_attr(attr, 'is_paw',         upf%tpawp,      newline=.true.)
      CALL iotk_write_attr(attr, 'is_coulomb',     upf%tcoulombp,  newline=.true.)
      !
      CALL iotk_write_attr(attr, 'has_so',         upf%has_so,     newline=.true.)
      CALL iotk_write_attr(attr, 'has_wfc',        upf%has_wfc,    newline=.true.)
      !EMINE
      CALL iotk_write_attr(attr, 'has_gipaw',      upf%has_gipaw,  newline=.true.)
      CALL iotk_write_attr(attr, 'paw_as_gipaw',      upf%paw_as_gipaw,  newline=.true.)
      !
      CALL iotk_write_attr(attr, 'core_correction',upf%nlcc,       newline=.true.)
      CALL iotk_write_attr(attr, 'functional',     TRIM(upf%dft),  newline=.true.)
      CALL iotk_write_attr(attr, 'z_valence',      upf%zp,         newline=.true.)
      CALL iotk_write_attr(attr, 'total_psenergy', upf%etotps,     newline=.true.)
      CALL iotk_write_attr(attr, 'wfc_cutoff',     upf%ecutwfc,    newline=.true.)
      CALL iotk_write_attr(attr, 'rho_cutoff',     upf%ecutrho,    newline=.true.)
      CALL iotk_write_attr(attr, 'l_max',          upf%lmax,       newline=.true.)
      CALL iotk_write_attr(attr, 'l_max_rho',      upf%lmax_rho,   newline=.true.)
      CALL iotk_write_attr(attr, 'l_local',        upf%lloc,       newline=.true.)
      CALL iotk_write_attr(attr, 'mesh_size',      upf%mesh,       newline=.true.)
      CALL iotk_write_attr(attr, 'number_of_wfc',  upf%nwfc,       newline=.true.)
      CALL iotk_write_attr(attr, 'number_of_proj', upf%nbeta,      newline=.true.)
      CALL iotk_write_empty(u, 'PP_HEADER', attr=attr)
      !
      !CALL iotk_write_end(u, 'PP_HEADER')
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
    SUBROUTINE write_semilocal(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      INTEGER :: nb, l, ind
      !
      CALL iotk_write_begin(u, 'PP_SEMILOCAL')
      !
      ! Write V_l(r)
      DO nb = 1,upf%nbeta
         l = upf%lll(nb)
         ind = 1
         CALL iotk_write_attr(attr, 'L',l, first=.true.)
         IF ( upf%has_so ) THEN
            CALL iotk_write_attr(attr, 'J', upf%jjj(nb))
            IF ( l > 0 .AND. ABS (upf%jjj(nb)-l-0.5_dp) < 0.001_dp) ind = 2
         ENDIF
         CALL iotk_write_dat(u, 'PP_VNL'//iotk_index(l), &
              upf%vnl(:,l,ind), attr=attr, columns=4)
      END DO
      !
      CALL iotk_write_end(u, 'PP_SEMILOCAL')
      !
    END SUBROUTINE write_semilocal
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
         CALL iotk_write_attr(attr, 'ultrasoft_cutoff_radius',upf%rcutus(nb))
         CALL iotk_write_dat(u, 'PP_BETA'//iotk_index( nb ), &
              upf%beta(:,nb), attr=attr, columns=4)
      ENDDO
      !
      ! Write the hamiltonian terms D_ij
      CALL iotk_write_dat(u, 'PP_DIJ', upf%dion, columns=4)
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
         CALL iotk_write_attr(attr,'l_max_aug',      upf%paw%lmax_aug)
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
         CALL iotk_write_dat(u, 'PP_MULTIPOLES', upf%paw%augmom, columns=4)
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
               isnull = .false. ! omit functions that are multiplied by zero
               IF( upf%tpawp ) isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
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
      CALL iotk_write_attr(attr, 'index',         nw, first=.true.)
      CALL iotk_write_attr(attr, 'label',         upf%els(nw))
      CALL iotk_write_attr(attr, 'l',             upf%lchi(nw))
      CALL iotk_write_attr(attr, 'occupation',    upf%oc(nw))
      CALL iotk_write_attr(attr, 'n',             upf%nchi(nw))
      CALL iotk_write_attr(attr, 'pseudo_energy', upf%epseu(nw))
      CALL iotk_write_attr(attr, 'cutoff_radius', upf%rcut_chi(nw))
      CALL iotk_write_attr(attr, 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw))
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

 SUBROUTINE write_full_wfc(u, upf)
   IMPLICIT NONE
   INTEGER,INTENT(IN)           :: u    ! i/o unit
   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
   INTEGER :: ierr ! /= 0 if something went wrong
   !
   CHARACTER(len=iotk_attlenx) :: attr
   !
   INTEGER :: nb

   IF(.not. upf%has_wfc) RETURN

   CALL iotk_write_attr(attr, 'number_of_wfc', upf%nbeta, first=.true.)
   CALL iotk_write_begin(u, 'PP_FULL_WFC', attr=attr)
   ! All-electron wavefunctions corresponding to beta functions
   DO nb = 1,upf%nbeta
      CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
      CALL iotk_write_attr(attr, 'label',      upf%els_beta(nb))
      CALL iotk_write_attr(attr, 'l',          upf%lll(nb))
      CALL iotk_write_dat(u, 'PP_AEWFC'//iotk_index(nb), &
           upf%aewfc(:,nb), columns=4, attr=attr)
   ENDDO
   IF (upf%has_so.and.upf%tpawp) THEN
      DO nb = 1,upf%nbeta
         CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
         CALL iotk_write_attr(attr, 'label',      upf%els_beta(nb))
         CALL iotk_write_attr(attr, 'l',          upf%lll(nb))
         CALL iotk_write_attr(attr, 'j',          upf%jjj(nb))
         CALL iotk_write_dat(u, 'PP_AEWFC_REL'//iotk_index(nb), &
              upf%paw%aewfc_rel(:,nb), columns=4, attr=attr)
      ENDDO
   ENDIF
   ! Pseudo wavefunctions 
   DO nb = 1,upf%nbeta
      CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
      CALL iotk_write_attr(attr, 'label',      upf%els_beta(nb))
      CALL iotk_write_attr(attr, 'l',          upf%lll(nb))
      CALL iotk_write_dat(u, 'PP_PSWFC'//iotk_index(nb), &
           upf%pswfc(:,nb), columns=4, attr=attr)
   ENDDO
   ! Finalize
   CALL iotk_write_end(u, 'PP_FULL_WFC')

 END SUBROUTINE write_full_wfc

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
   !
   CALL iotk_write_end(u, 'PP_PAW')
   !
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
   IF (.not. upf%has_gipaw) RETURN

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
   ! Only write core orbitals in the PAW as GIPAW case
   IF (upf%paw_as_gipaw) THEN
      CALL iotk_write_end(u, 'PP_GIPAW')
      RETURN
   ENDIF
   !
   ! Write valence all-electron and pseudo orbitals
   CALL iotk_write_attr(attr, 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels, first=.true.)
   CALL iotk_write_begin(u, 'PP_GIPAW_ORBITALS', attr=attr)
   !
   DO nb = 1,upf%gipaw_wfs_nchannels
      CALL iotk_write_attr(attr, 'index', nb, first=.true.)
      CALL iotk_write_attr(attr, 'label', upf%gipaw_wfs_el(nb))
      CALL iotk_write_attr(attr, 'l',     upf%gipaw_wfs_ll(nb))
      CALL iotk_write_attr(attr, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb))
      CALL iotk_write_attr(attr, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb))
      CALL iotk_write_begin(u, 'PP_GIPAW_ORBITAL'//iotk_index(nb), attr=attr)
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
   CALL iotk_write_dat(u, 'PP_GIPAW_VLOCAL_AE', &
        upf%gipaw_vlocal_ae(:), columns=4)
   CALL iotk_write_dat(u, 'PP_GIPAW_VLOCAL_PS', &
        upf%gipaw_vlocal_ps(:), columns=4)
   CALL iotk_write_end(u, 'PP_GIPAW_VLOCAL')
   !
   CALL iotk_write_end(u, 'PP_GIPAW')

   RETURN
 END SUBROUTINE write_gipaw
 !
 ! Remove '<' and '>' from string, replacing them with '/', necessary
 ! or iotk will complain while read-skipping PP_INFO section.
 FUNCTION CHECK(in) RESULT (out)
   CHARACTER(len=*) :: in
#if defined(__PGI)
   INTEGER, PARAMETER :: length = 255
   CHARACTER(len=length) :: out
#else
   CHARACTER(len=len(in)) :: out
#endif   
   INTEGER :: i
   DO i = 1,len(in)
      IF ( in(i:i) == '<' .or. in(i:i) == '>' ) THEN
         out(i:i) = '/'
      ELSE
         out(i:i) = in(i:i)
      ENDIF
   ENDDO
 END FUNCTION CHECK
END SUBROUTINE write_upf_v2

SUBROUTINE deallocate_pseudo_config(conf)
 TYPE(pseudo_config),INTENT(INOUT) :: conf
 if (associated(conf%els)   ) deallocate(conf%els)
 if (associated(conf%nns)   ) deallocate(conf%nns)
 if (associated(conf%lls)   ) deallocate(conf%lls)
 if (associated(conf%ocs)   ) deallocate(conf%ocs)
 if (associated(conf%rcut)  ) deallocate(conf%rcut)
 if (associated(conf%rcutus)) deallocate(conf%rcutus)
 if (associated(conf%enls)  ) deallocate(conf%enls)

END SUBROUTINE deallocate_pseudo_config

END MODULE write_upf_v2_module
