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
  USE pseudo_types, ONLY: pseudo_upf, pseudo_config, deallocate_pseudo_config
  USE radial_grids, ONLY: radial_grid_type
  USE FoX_wxml
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_upf_v2 


  INTERFACE write_data
   MODULE PROCEDURE  write_columns, write_columns_2, write_columns_3, write_columns_4         
  END INTERFACE 
                             

CONTAINS



  !--------------------------------------------------------+
  SUBROUTINE write_upf_v2(xf, upf, conf, u_input)
    !------------------------------------------------------+
    !! Write pseudopotential in UPF format version 2, uses FoX_wxml
    !
    IMPLICIT NONE
    TYPE(xmlf_t),INTENT(IN OUT)     :: xf   
    !! FoX_xml file descriptor  
    !! fortran unit identifier, used for direct writin
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    ! optional: configuration used to generate the pseudopotential
    TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
    ! optional: unit pointing to input file containing generation data
    INTEGER, OPTIONAL, INTENT(IN) :: u_input
    INTEGER                       :: i
    !
    !
    ! Initialize the file
    CALL xml_newElement(xf, "UPF")
      CALL xml_addAttribute(xf, 'version', TRIM(upf%nv))
      !
      ! Write human-readable header
      CALL write_info(xf, upf, conf, u_input)
      ! Write machine-readable header
      CALL write_header(xf, upf)
      ! Write radial grid mesh
      CALL write_mesh(xf, upf)
      ! Write non-linear core correction charge
      IF(upf%nlcc) THEN
         CALL xml_newElement(xf, 'PP_NLCC')
            CALL xml_addAttribute(xf,'columns',4)
            CALL write_data(xf, upf%rho_atc) 
         CALL xml_EndElement(xf, 'PP_NLCC')
      ENDIF 
      ! Write local potential
      IF(.not. upf%tcoulombp) THEN
         CALL xml_newElement(xf, 'PP_LOCAL')
            CALL xml_addAttribute(xf, 'columns', 4)
            CALL write_data(xf, upf%vloc, tag = 'PP_LOCAL')
      ELSE
         CALL xml_newElement(xf,'PP_LOCAL')
            CALL xml_addAttribute(xf,'type','1/r')
            CALL xml_addComment(xf,'Coulomb 1/r potential')
         CALL xml_EndElement(xf,'PP_LOCAL')
      ENDIF
      ! Write potentials in semilocal form (if existing)
      IF ( upf%typ == "SL" ) CALL write_semilocal(xf, upf)
      ! Write nonlocal components: projectors, augmentation, hamiltonian elements
      CALL write_nonlocal(xf, upf)
      ! Write initial pseudo wavefunctions
      ! (usually only wfcs with occupancy > 0)
      CALL write_pswfc(xf, upf)
      ! If included, write all-electron and pseudo wavefunctions
      CALL write_full_wfc(xf, upf)
      ! Write valence atomic density (used for initial density)
      CALL xml_newElement(xf, 'PP_RHOATOM')
         CALL write_data(xf, upf%rho_at, tag = 'PP_RHOATOM')
      ! Write additional info for full-relativistic calculation
      CALL write_spin_orb(xf, upf)
      ! Write additional data for PAW (All-electron charge, wavefunctions, vloc..)
      CALL write_paw(xf, upf)
      ! Write additional data for GIPAW reconstruction
      CALL write_gipaw(xf, upf)
      !
      ! Close the file (not the unit!)
      CALL xml_EndElement(xf, 'UPF')
      CALL xml_Close(xf)
      !
  CONTAINS
    !
    SUBROUTINE write_info(u, upf, conf, u_input)
      ! Write human-readable header
      !
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)          :: u   ! i/o unit: write to unit u
      TYPE(pseudo_upf),INTENT(IN) :: upf  ! the pseudo data
      ! optional: configuration used to generate the pseudopotential
      TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
      INTEGER, OPTIONAL, INTENT(IN) :: u_input ! read input data from u_input
      !
      INTEGER :: nb ! aux counter
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(len=4096) :: char_buff
      CHARACTER(LEN=256)  :: line
      LOGICAL :: opnd
      !
      CALL xml_newElement(u, 'PP_INFO')
      !
      char_buff = ""
      char_buff = trim(char_buff)// new_line('a') // TRIM(upf%generated)
      char_buff = trim(char_buff)// new_line('a') // 'Author: '//TRIM(upf%author)
      char_buff = trim(char_buff)// new_line('a') // 'Generation date: '//TRIM(upf%date) 
      char_buff = trim(char_buff)// new_line('a') // 'Pseudopotential type: '//TRIM(upf%typ) 
      char_buff = trim(char_buff)// new_line('a') // 'Element: '//TRIM(upf%psd) 
      char_buff = trim(char_buff)// new_line('a') //'Functional: '//TRIM(upf%dft) 
      !
      ! Cutoff Information
      WRITE(line, '(4x,a,f5.0,a)') &
           'Suggested minimum cutoff for wavefunctions:',upf%ecutwfc,' Ry'
      char_buff = trim(char_buff) // new_line('a') // trim(line)
      !
      WRITE(line, '(4x,a,f5.0,a)') &
           'Suggested minimum cutoff for charge density:',&
           upf%ecutrho,' Ry'
      char_buff = trim(char_buff) // new_line('a') // trim(line)
      !
      ! Write relativistic information
      IF (TRIM(upf%rel)=='full') THEN
         WRITE(line, '(4x,a)', err=100) &
              "The Pseudo was generated with a Fully-Relativistic Calculation"
      ELSE IF (TRIM(upf%rel)=='scalar') THEN
         WRITE(line, '(4x,a)', err=100) &
              "The Pseudo was generated with a Scalar-Relativistic Calculation"
      ELSE
         WRITE(line, '(4x,a)', err=100) &
              "The Pseudo was generated with a Non-Relativistic Calculation"
      ENDIF
      char_buff = trim(char_buff) // new_line('a') // trim(line)
      !
      ! Write local potential information
      IF (upf%lloc >= 0 ) THEN
         WRITE(line, '(4x,a,i3,f9.4)', err=100) &
              "L component and cutoff radius for Local Potential:", upf%lloc, upf%rcloc
      ELSE IF (upf%lloc == -1 ) THEN
         WRITE(line, '(4x,a,f9.4)', err=100) &
              "Local Potential by smoothing AE potential with Bessel fncs, cutoff radius:", upf%rcloc
      ELSE IF (upf%lloc == -2 ) THEN
         WRITE(line, '(4x,a,f9.4)', err=100) &
              "Local Potential according to Troullier-Martins recipe, cutoff radius:", upf%rcloc
      ELSE
         WRITE(line, '(4x,a,i3,f9.4)', err=100) &
              "Local Potential: unknown format, L component and cutoff radius:",upf%lloc, upf%rcloc
      ENDIF
      char_buff = trim(char_buff) // new_line('a') // trim(line)
      !
      IF (upf%has_so) THEN
         WRITE(line, '(4x,a,i3,f9.4)', err=100) &
           "Pseudopotential contains additional information for spin-orbit calculations."
         char_buff = trim(char_buff) // new_line('a') // trim(line)
      ENDIF
      !
      IF (upf%has_gipaw) THEN
         WRITE(line, '(4x,a,i3,f9.4)', err=100) &
           "Pseudopotential contains additional information for GIPAW reconstruction."
         char_buff = trim(char_buff) // new_line('a') // trim(line)
      END IF
      !
      ! Write valence orbitals information
      WRITE(line, '(4x,a)') 'Valence configuration: '
      char_buff = trim(char_buff) // new_line('a') // trim(line)
      !
      WRITE(line, '(4x,a2,2a3,a6,2a11,1a13)', err=100) &
           "nl"," pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
      char_buff = trim(char_buff) // new_line('a') // trim(line) 
      DO nb = 1, upf%nwfc
         IF(upf%oc(nb) >= 0._dp) THEN
            WRITE(line, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
                 upf%els(nb), upf%nchi(nb), &
                 upf%lchi(nb), upf%oc(nb), upf%rcut_chi(nb), &
                 upf%rcutus_chi(nb), upf%epseu(nb)
            char_buff = TRIM(char_buff) // new_line('a') // TRIM(line)
         ENDIF
      END DO
      IF( present(conf) ) THEN
         WRITE(line, '(4x,a)') 'Generation configuration:'
         char_buff = TRIM(char_buff) //new_line('a') //  TRIM(line)
         DO nb = 1,conf%nwfs
            WRITE(line, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
                 conf%els(nb), conf%nns(nb), &
                 conf%lls(nb), conf%ocs(nb), conf%rcut(nb), &
                 conf%rcutus(nb), conf%enls(nb)
            char_buff = TRIM(char_buff) // new_line('a') // TRIM(line)
         ENDDO
         WRITE(line,'(4x,2a)') 'Pseudization used: ',TRIM(conf%pseud)
         char_buff = TRIM(char_buff) // new_line('a') // TRIM(line)
      ELSE
         WRITE(line, '(4x,a)') 'Generation configuration: not available.'
         char_buff = TRIM(char_buff) // new_line('a') // TRIM(line)
      ENDIF

      IF(TRIM(upf%comment) /= ' ') THEN
         WRITE(line, '(4x,"Comment:",/,4x,a)', err=100) TRIM(upf%comment)
         char_buff = TRIM(char_buff) // new_line('a') // TRIM(line)
      END IF
      char_buff = TRIM(char_buff) // new_line('a')
      CALL xml_AddCharacters(u, TRIM(char_buff), parsed = .FALSE.) 
      !
      IF ( PRESENT(u_input) ) THEN
         !
         ! copy content of input file used in pseudopotential generation
         !
         INQUIRE (unit=u_input, opened=opnd)
         IF (opnd) THEN
            char_buff="" 
            CALL xml_newElement(u, 'PP_INPUTFILE') 
            REWIND (unit=u_input)
10          READ (u_input, '(A)',end=20,err=25) line
            char_buff = TRIM(char_buff) // new_line('a') // TRIM(line) 
            GO TO 10
25          CALL infomsg('write_upf_v2::write_inputfile', 'problem writing input data')
            !
20          char_buff = TRIM(char_buff) // new_line('a') 
            CALL xml_AddCharacters(u, TRIM(char_buff), parsed = .FALSE.)
            CALL xml_EndElement (u, 'PP_INPUTFILE')
         ELSE
            CALL infomsg('write_upf_v2::write_inputfile', 'input file not open')
         END IF
         !
      END IF
      !
      CALL xml_EndElement(u, 'PP_INFO')
      CALL xml_AddComment(u, new_line('a')//' END OF HUMAN READABLE SECTION '//new_line('a'))
      !
      RETURN
100   CALL errore('write_upf_v2::write_info', 'Writing pseudo file', 1)
      !
    END SUBROUTINE write_info
    !
    !
    SUBROUTINE write_header(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      INTEGER :: nw
      !
      ! Write HEADER section with some initialization data
      CALL xml_newElement  (u, 'PP_HEADER')
      CALL xml_addAttribute(u, 'generated',  TRIM(upf%generated) )
      CALL xml_addAttribute(u, 'author',     TRIM(upf%author)    )
      CALL xml_addAttribute(u, 'date',       TRIM(upf%date)   )
      CALL xml_addAttribute(u, 'comment',    TRIM(upf%comment) )
      !
      CALL xml_addAttribute(u, 'element',  upf%psd )
      CALL xml_addAttribute(u, 'pseudo_type',    TRIM(upf%typ) )
      CALL xml_addAttribute (u, 'relativistic',   TRIM(upf%rel) )
      !
      CALL xml_addAttribute(u, 'is_ultrasoft',   upf%tvanp )
      CALL xml_addAttribute(u, 'is_paw',         upf%tpawp )
      CALL xml_addAttribute(u, 'is_coulomb',     upf%tcoulombp )
      !
      CALL xml_addAttribute(u, 'has_so',         upf%has_so )
      CALL xml_addAttribute(u, 'has_wfc',        upf%has_wfc )
      !EMINE
      CALL xml_addAttribute(u, 'has_gipaw',      upf%has_gipaw )
      CALL xml_addAttribute(u, 'paw_as_gipaw',      upf%paw_as_gipaw )
      !
      CALL xml_addAttribute(u, 'core_correction',upf%nlcc        )
      CALL xml_addAttribute(u, 'functional',     TRIM(upf%dft)   )
      CALL xml_addAttribute(u, 'z_valence',      upf%zp          )
      CALL xml_addAttribute(u, 'total_psenergy', upf%etotps      )
      CALL xml_addAttribute(u, 'wfc_cutoff',     upf%ecutwfc     )
      CALL xml_addAttribute(u, 'rho_cutoff',     upf%ecutrho     )
      CALL xml_addAttribute(u, 'l_max',          upf%lmax        )
      CALL xml_addAttribute(u, 'l_max_rho',      upf%lmax_rho    )
      CALL xml_addAttribute(u, 'l_local',        upf%lloc        )
      CALL xml_addAttribute(u, 'mesh_size',      upf%mesh        )
      CALL xml_addAttribute(u, 'number_of_wfc',  upf%nwfc        )
      CALL xml_addAttribute(u, 'number_of_proj', upf%nbeta       )
      CALL xml_EndElement(u, 'PP_HEADER' )
      RETURN
    END SUBROUTINE write_header
    !
    SUBROUTINE write_mesh(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      CALL xml_newElement(u, 'PP_MESH')
      IF (upf%dx    .GT. 0.d0)  CALL xml_addAttribute(u, 'dx',   upf%dx )
      IF (upf%mesh  .GT. 0 )    CALL xml_addAttribute(u, 'mesh', upf%mesh)
      IF (upf%dx    .GT. 0.d0)  CALL xml_addAttribute(u, 'xmin', upf%xmin)
      IF (upf%rmax  .GT. 0.d0)  CALL xml_addAttribute(u, 'rmax', upf%rmax)
      IF (upf%zmesh .GT. 0.d0)  CALL xml_addAttribute(u, 'zmesh',upf%zmesh)
      !
      CALL xml_newElement(u, 'PP_R') 
         CALL write_data(u, upf%r, tag = 'PP_R')
      !
      CALL xml_newElement(u, 'PP_RAB')
         CALL write_data(u, upf%rab, tag = 'PP_RAB') 
      !
      CALL xml_EndElement(u, 'PP_MESH' )
      !
      RETURN
    END SUBROUTINE write_mesh
    !
    SUBROUTINE write_semilocal(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    
      TYPE(pseudo_upf),INTENT(IN)  :: upf  
      INTEGER :: ierr 
      CHARACTER(LEN=10)  :: tag_vnl 
      INTEGER :: nb, l, ind
      !
      CALL xml_newElement(u, 'PP_SEMILOCAL')
      !
      ! Write V_l(r)
      DO nb = 1,upf%nbeta
         l = upf%lll(nb)
         ind = 1
         IF (upf%has_so) THEN
            IF ( l > 0 .AND. ABS (upf%jjj(nb)-l-0.5_dp) < 0.001_dp) ind = 2
         END IF 
         WRITE(tag_vnl,'("PP_VNL.",I0)') ind 
         CALL xml_newElement(u, TRIM(tag_vnl))
            CALL xml_addAttribute(u, 'L',l )
            IF (upf%has_so) CALL xml_addAttribute(u, 'J', upf%jjj(nb))
            CALL write_data(u, upf%vnl(:,l,ind), tag = trim(tag_vnl)) 
      END DO
      !
      CALL xml_EndElement(u, 'PP_SEMILOCAL')
      !
    END SUBROUTINE write_semilocal
    !
    SUBROUTINE write_nonlocal(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u  
      TYPE(pseudo_upf),INTENT(IN)  :: upf
      INTEGER :: ierr 
      CHARACTER(LEN=10) :: tag_beta
      CHARACTER(LEN=15) :: tag_qijl
      !
      INTEGER :: nb,mb,ln,lm,l,nmb
      LOGICAL :: isnull
      REAL(DP),ALLOCATABLE :: aux(:)
      INTEGER              :: auxlen
      !
      IF (upf%tcoulombp) RETURN
      !
      CALL xml_newElement(u, 'PP_NONLOCAL')
      !
      ! Write the projectors:
      DO nb = 1,upf%nbeta
         WRITE (tag_beta, '("PP_BETA.",I0)') nb
         !tag_beta="PP_BETA."
         CALL xml_newElement(u, trim(tag_beta)) 
            CALL xml_addAttribute(u, 'index',                  nb )
            CALL xml_addAttribute(u, 'label',                  upf%els_beta(nb))
            CALL xml_addAttribute(u, 'angular_momentum',       upf%lll(nb))
            CALL xml_addAttribute(u, 'cutoff_radius_index',    upf%kbeta(nb))
            CALL xml_addAttribute(u, 'cutoff_radius',          upf%rcut(nb))
            CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius',upf%rcutus(nb))
            CALL write_data(u, upf%beta(:,nb), tag = TRIM(tag_beta)) 
      ENDDO
      !
      ! Write the hamiltonian terms D_ij
      CALL xml_newElement(u, 'PP_DIJ')
         CALL write_data(u, upf%dion, tag = 'PP_DIJ') 
      !
      ! Write the augmentation charge section
      augmentation : IF(upf%tvanp .or. upf%tpawp) THEN
         CALL xml_newElement(u, 'PP_AUGMENTATION')
            CALL xml_addAttribute(u, 'q_with_l', upf%q_with_l )
            CALL xml_addAttribute(u, 'nqf',      upf%nqf)
            CALL xml_addAttribute(u, 'nqlc',     upf%nqlc)
            IF (upf%tpawp) THEN
               CALL xml_addAttribute(u,'shape',          TRIM(upf%paw%augshape))
               CALL xml_addAttribute(u,'cutoff_r',       upf%paw%raug, fmt ='s16')
               CALL xml_addAttribute(u,'cutoff_r_index', upf%paw%iraug)
               CALL xml_addAttribute(u,'augmentation_epsilon',upf%qqq_eps, fmt='s16')
               CALL xml_addAttribute(u,'l_max_aug',      upf%paw%lmax_aug)
            ENDIF
            !
            !
            ! Write the integrals of the Q functions
            CALL xml_newElement(u, 'PP_Q')
            CALL write_data(u, upf%qqq, tag = 'PP_Q')
            !
            ! Write charge multipoles (only if PAW)
            IF ( upf%tpawp ) THEN
               CALL xml_addComment(u, ' augmentation charge multipoles (only for PAW) ')
               CALL xml_newElement(u, 'PP_MULTIPOLES') 
                  CALL write_data(u, upf%paw%augmom, tag = 'PP_MULTIPOLES')
            ENDIF
            !
            ! Write polinomial coefficients for Q_ij expansion at small radius
            IF ( upf%nqf > 0) THEN
               CALL xml_addComment(u, ' polinomial expansion of Q_ij at small radius ')
               CALL xml_newElement(u, 'PP_QFCOEF')
                  CALL write_data(u, upf%qfcoef, tag = 'PP_QFCOEF')
               !
               CALL xml_newElement(u, 'PP_RINNER')
                  CALL write_data(u, upf%rinner, tag = 'PP_RINNER')
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
                        WRITE(tag_qijl, '("PP_QIJL.",I0,".",I0,".",I0)') nb,mb,l
                        CALL xml_newElement(u, TRIM(tag_qijl)) 
                           CALL xml_addAttribute(u, 'first_index',  nb )
                           CALL xml_addAttribute(u, 'second_index', mb)
                           CALL xml_addAttribute(u, 'composite_index', nmb)
                           CALL xml_addAttribute(u, 'angular_momentum', l)
                           !
                           isnull = .false. ! omit functions that are multiplied by zero
                           IF( upf%tpawp ) isnull = (abs(upf%paw%augmom(nb,mb,l)) < upf%qqq_eps)
                           !
                           IF ( isnull ) THEN
                              CALL xml_addAttribute(u, 'is_null', isnull)
                           ELSE
                              CALL write_data(u, upf%qfuncl(:,nmb,l))
                           ENDIF
                        CALL xml_EndElement(u, TRIM(tag_qijl)) 
                     ENDDO
                  ELSE
                     WRITE(tag_qijl, '("PP_QIJ.",I0,".",I0)') nb, mb
                     CALL xml_newElement(u, TRIM(tag_qijl)) 
                        CALL xml_addAttribute(u, 'first_index',  nb )
                        CALL xml_addAttribute(u, 'second_index', mb)
                        CALL xml_addAttribute(u, 'composite_index', nmb)
                        !
                        isnull = .false. ! omit functions that are multiplied by zero
                        IF( upf%tpawp ) isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
                        IF ( isnull ) THEN
                           CALL xml_addAttribute(u, 'is_null', isnull)
                        ELSE
                           CALL write_data(u, upf%qfunc(:,nmb))
                        ENDIF
                     CALL xml_EndElement(u, TRIM(tag_qijl)) 
                  ENDIF
               ENDDO
            ENDDO
            !
            CALL xml_EndElement(u, 'PP_AUGMENTATION')
            !
         ENDIF augmentation
         !
      CALL xml_EndElement(u, 'PP_NONLOCAL')
      !
      RETURN
    END SUBROUTINE write_nonlocal
   !
   SUBROUTINE write_pswfc(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t), INTENT(INOUT)  :: u    
      TYPE(pseudo_upf),INTENT(IN)  :: upf  
      INTEGER :: ierr
      CHARACTER(LEN=12) ::   tag_chi
      !
      INTEGER :: nw
      !
      CALL xml_newElement(u, 'PP_PSWFC')
      !
      DO nw = 1,upf%nwfc
         WRITE(tag_chi, '("PP_CHI.",I0)') nw 
         CALL xml_newElement(u, TRIM(tag_chi)) 
            CALL xml_addAttribute(u, 'index',         nw )
            CALL xml_addAttribute(u, 'label',         upf%els(nw))
            CALL xml_addAttribute(u, 'l',             upf%lchi(nw))
            CALL xml_addAttribute(u, 'occupation',    upf%oc(nw))
            IF ( upf%nchi(nw) .GT. upf%lchi(nw) )  CALL xml_addAttribute(u, 'n',             upf%nchi(nw))
            IF ( upf%epseu(nw) .GT. 0.0_DP)        CALL xml_addAttribute(u, 'pseudo_energy', upf%epseu(nw))
            IF ( upf%rcut_chi(nw) .GT. 0.0_DP )    CALL xml_addAttribute(u, 'cutoff_radius', upf%rcut_chi(nw))
            IF ( upf%rcutus_chi(nw) .GT. 0.0_DP )   CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw))
            CALL write_data(u, upf%chi(:,nw), tag = TRIM(tag_chi))
      ENDDO
      !
      CALL xml_EndElement(u, 'PP_PSWFC')
      !
      RETURN
   END SUBROUTINE write_pswfc
   !
   SUBROUTINE write_spin_orb(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u   
      TYPE(pseudo_upf),INTENT(IN)  :: upf 
      INTEGER :: ierr 
      CHARACTER(LEN=15) :: tag_relwfc, tag_relbeta
      !
      INTEGER :: nw, nb
      !
      IF (.not. upf%has_so) RETURN
      !
      CALL xml_newElement(u, 'PP_SPIN_ORB')
      !
         DO nw = 1,upf%nwfc
            WRITE(tag_relwfc, '("PP_RELWFC.",I0)') nw
            CALL xml_newElement(u, TRIM(tag_relwfc) ) 
               CALL xml_addAttribute(u, 'index' , nw)
               CALL xml_addAttribute(u, 'els',   upf%els(nw))
               CALL xml_addAttribute(u, 'nn',    upf%nn(nw))
               CALL xml_addAttribute(u, 'lchi',  upf%lchi(nw))
               CALL xml_addAttribute(u, 'jchi',  upf%jchi(nw))
               CALL xml_addAttribute(u, 'oc',    upf%oc(nw))
            CALL xml_EndElement(u, TRIM(tag_relwfc) )
         ENDDO
         !
         DO nb = 1,upf%nbeta
            WRITE(tag_relbeta, '("PP_RELBETA.",I0)') nb
            CALL xml_newElement(u, TRIM(tag_relbeta) ) 
               CALL xml_addAttribute(u, 'index', nb )
               CALL xml_addAttribute(u, 'lll',   upf%lll(nb))
               CALL xml_addAttribute(u, 'jjj',   upf%jjj(nb))
            CALL xml_EndElement(u, TRIM(tag_relbeta))
         ENDDO
         !
      CALL xml_EndElement(u, 'PP_SPIN_ORB')
      !
      RETURN
   END SUBROUTINE write_spin_orb
   !
   SUBROUTINE write_full_wfc(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    
      TYPE(pseudo_upf),INTENT(IN)  :: upf  
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      !
      INTEGER :: nb
      CHARACTER(LEN=20)  :: tag_aewfc 
      IF(.not. upf%has_wfc) RETURN
      CALL xml_newElement(u, 'PP_FULL_WFC')
         CALL xml_addAttribute(u, 'number_of_wfc', upf%nbeta )
         ! All-electron wavefunctions corresponding to beta functions
         DO nb = 1,upf%nbeta
            WRITE(tag_aewfc, '("PP_AEWFC.",I0)') nb  
            CALL xml_newElement(u, TRIM(tag_aewfc) )
               CALL xml_addAttribute(u, 'index',      nb )
               CALL xml_addAttribute(u, 'label',      upf%els_beta(nb))
               CALL xml_addAttribute(u, 'l',          upf%lll(nb))
               CALL write_data(u, upf%aewfc(:,nb), tag = TRIM(tag_aewfc))
         ENDDO
         IF (upf%has_so.and.upf%tpawp) THEN
            DO nb = 1,upf%nbeta
               WRITE(tag_aewfc, '("PP_AEWFC_REL.",I0)' ) nb
               CALL xml_newElement(u, TRIM(tag_aewfc)) 
                  CALL xml_addAttribute(u, 'index',      nb )
                  CALL xml_addAttribute(u, 'label',      upf%els_beta(nb))
                  CALL xml_addAttribute(u, 'l',          upf%lll(nb))
                  CALL xml_addAttribute(u, 'j',          upf%jjj(nb))
                  CALL write_data(u, upf%paw%aewfc_rel(:,nb), tag = TRIM(tag_aewfc))  
            ENDDO
         ENDIF
         ! Pseudo wavefunctions 
         DO nb = 1,upf%nbeta
            WRITE(tag_aewfc, '("PP_PSWFC.",I0)') nb 
            CALL xml_newElement(u, TRIM(tag_aewfc)) 
               CALL xml_addAttribute(u, 'index',      nb )
               CALL xml_addAttribute(u, 'label',      upf%els_beta(nb))
               CALL xml_addAttribute(u, 'l',          upf%lll(nb))
               CALL write_data(u, upf%pswfc(:,nb), tag = TRIM(tag_aewfc))
         ENDDO
         ! Finalize
      CALL xml_EndElement(u, 'PP_FULL_WFC')
   END SUBROUTINE write_full_wfc

   SUBROUTINE write_paw(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u   
      TYPE(pseudo_upf),INTENT(IN)  :: upf 
      INTEGER :: ierr 
      !
      !
      INTEGER :: nb
      IF (.not. upf%tpawp ) RETURN
      CALL xml_newElement(u, 'PP_PAW') 
         CALL xml_addAttribute(u, 'paw_data_format', upf%paw_data_format )
         CALL xml_addAttribute(u, 'core_energy',     upf%paw%core_energy)
         ! Full occupation (not only > 0 ones)
         CALL xml_newElement(u, 'PP_OCCUPATIONS') 
            CALL write_data(u, upf%paw%oc, tag = 'PP_OCCUPATIONS')
         ! All-electron core charge
         CALL xml_newElement(u,'PP_AE_NLCC')
            CALL write_data(u, upf%paw%ae_rho_atc, tag = 'PP_AE_NLCC')
         ! All-electron local potential
         CALL xml_newElement(u, 'PP_AE_VLOC')
            CALL write_data(u, upf%paw%ae_vloc, tag = 'PP_AE_VLOC')
         !
      CALL xml_EndElement(u, 'PP_PAW')
      !
      RETURN
   END SUBROUTINE write_paw
   !
   SUBROUTINE write_gipaw(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(LEN=30)  :: tag_aux
      INTEGER            :: auxlen
      !
      !
      INTEGER :: nb
      IF (.not. upf%has_gipaw) RETURN
      CALL xml_newElement(u, 'PP_GIPAW')
         CALL xml_addAttribute(u, 'gipaw_data_format', upf%gipaw_data_format )
         CALL xml_newElement(u, 'PP_GIPAW_CORE_ORBITALS') 
            CALL xml_addAttribute(u, 'number_of_core_orbitals', upf%gipaw_ncore_orbitals )
            DO nb = 1,upf%gipaw_ncore_orbitals
               WRITE(tag_aux, '("PP_GIPAW_CORE_ORBITAL.",I0)') nb
               CALL xml_newElement(u, TRIM(tag_aux)) 
                  CALL xml_addAttribute(u, 'index', nb )
                  CALL xml_addAttribute(u, 'label', upf%gipaw_core_orbital_el(nb))
                  CALL xml_addAttribute(u, 'n',     upf%gipaw_core_orbital_n(nb))
                  CALL xml_addAttribute(u, 'l',     upf%gipaw_core_orbital_l(nb))
                  CALL write_data(u, upf%gipaw_core_orbital(:,nb), tag = TRIM(tag_aux))
            ENDDO
         CALL xml_EndElement(u, 'PP_GIPAW_CORE_ORBITALS')
         !
         ! Only write core orbitals in the PAW as GIPAW case
      IF (upf%paw_as_gipaw) THEN
         CALL xml_EndElement(u, 'PP_GIPAW')
         RETURN
      ENDIF
      !
      ! Write valence all-electron and pseudo orbitals
         CALL xml_newElement(u, 'PP_GIPAW_ORBITALS') 
            CALL xml_addAttribute(u, 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels )
            !
            DO nb = 1,upf%gipaw_wfs_nchannels
               WRITE(tag_aux, '("PP_GIPAW_ORBITAL.",I0)') nb 
               CALL xml_newElement(u, TRIM(tag_aux)) 
                  CALL xml_addAttribute(u, 'index', nb )
                  CALL xml_addAttribute(u, 'label', upf%gipaw_wfs_el(nb))
                  CALL xml_addAttribute(u, 'l',     upf%gipaw_wfs_ll(nb))
                  CALL xml_addAttribute(u, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb))
                  CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb))
                  CALL xml_newElement(u, 'PP_GIPAW_WFS_AE')
                     CALL write_data(u, upf%gipaw_wfs_ae(:,nb), tag = 'PP_GIPAW_WFS_AE')
                  CALL xml_newElement(u, 'PP_GIPAW_WFS_PS')
                     CALL write_data(u, upf%gipaw_wfs_ps(:,nb), tag = 'PP_GIPAW_WFS_PS')
               !
               CALL xml_EndElement(u, trim(tag_aux) )
            ENDDO
         CALL xml_EndElement(u, 'PP_GIPAW_ORBITALS')
   !
   ! Write all-electron and pseudo local potentials
   CALL xml_newElement(u, 'PP_GIPAW_VLOCAL')
      CALL xml_newElement(u, 'PP_GIPAW_VLOCAL_AE')
         CALL write_data(u, upf%gipaw_vlocal_ae(:), tag = 'PP_GIPAW_VLOCAL_AE')
      CALL xml_newElement(u, 'PP_GIPAW_VLOCAL_PS')
         CALL write_data(u, upf%gipaw_vlocal_ps(:), tag = 'PP_GIPAW_VLOCAL_PS')
   CALL xml_EndElement(u, 'PP_GIPAW_VLOCAL')
   !
   CALL xml_EndElement(u, 'PP_GIPAW')

   RETURN
 END SUBROUTINE write_gipaw
 !
END SUBROUTINE write_upf_v2

 
SUBROUTINE write_columns(xf, data, tag, columns )
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)       :: xf
   CHARACTER(LEN=*),OPTIONAL        :: tag
   REAL(DP), INTENT(IN)             :: data(:)
   INTEGER,OPTIONAL,INTENT(IN)      :: columns
   ! 
   INTEGER                          :: inc = 4, length, i  
   IF (PRESENT(columns) ) inc = columns
   length = SIZE(data) 
   CALL xml_addNewLine(xf) 
   DO i = 1, length, inc 
      CALL xml_AddCharacters(xf, data(i:MIN(length,i+inc-1) ), fmt = 's16')
      CALL xml_addNewLine(xf) 
   END DO
   IF (PRESENT(tag)) CALL xml_EndElement(xf, TRIM(tag)) 
END SUBROUTINE write_columns


SUBROUTINE write_columns_2(xf, data, tag, columns) 
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)       :: xf
   CHARACTER(LEN=*),OPTIONAL        :: tag
   REAL(DP),INTENT(IN)              :: data(:,:)
   INTEGER,OPTIONAL,INTENT(IN)      :: columns
   !
   REAL(DP),ALLOCATABLE             :: aux(:) 
   ALLOCATE (aux(SIZE(data)))
   aux = RESHAPE(data, [SIZE(data)]) 
   CALL write_columns(xf, aux, tag, columns) 
   DEALLOCATE(aux)
END SUBROUTINE write_columns_2

SUBROUTINE write_columns_3(xf, data, tag, columns) 
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)       :: xf
   CHARACTER(LEN=*)                 :: tag
   REAL(DP),INTENT(IN)              :: data(:,:,:)
   INTEGER,OPTIONAL,INTENT(IN)      :: columns
   !
   REAL(DP),ALLOCATABLE             :: aux(:) 
   ALLOCATE (aux(SIZE(data))) 
   aux = RESHAPE(data, [SIZE(data)]) 
   CALL write_columns(xf, aux, tag, columns) 
   DEALLOCATE(aux)
END SUBROUTINE write_columns_3

SUBROUTINE write_columns_4(xf, data, tag, columns) 
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)       :: xf
   CHARACTER(LEN=*)                 :: tag
   REAL(DP),INTENT(IN)               :: data(:,:,:,:)
   INTEGER,OPTIONAL,INTENT(IN)      :: columns
   !
   REAL(DP),ALLOCATABLE             :: aux(:) 
   ALLOCATE (aux(SIZE(data))) 
   aux = RESHAPE(data, [SIZE(data)]) 
   CALL write_columns(xf, aux, tag, columns) 
   DEALLOCATE(aux)
END SUBROUTINE write_columns_4


END MODULE write_upf_v2_module
