!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
! TB
! included monopole potential similar to the saw function of the sawtooth
! electric field, search 'TB'
!------------------------------------------------------------------------

!------------------------------------------------------------------------------!
  MODULE cell_base
!------------------------------------------------------------------------------!

    USE kinds, ONLY : DP
    USE constants, ONLY : pi, bohr_radius_angs
    USE io_global, ONLY : stdout
!
    IMPLICIT NONE
    SAVE
    !
    !  ibrav: index of the bravais lattice (see latgen.f90)
    INTEGER          :: ibrav
    !  celldm: old-style parameters of the simulation cell (se latgen.f90)
    REAL(DP) :: celldm(6) = (/ 0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP /)
    !  traditional crystallographic cell parameters (alpha=cosbc and so on)
    
    REAL(DP) :: a, b, c, cosab, cosac, cosbc
    ! format of input cell parameters:
    ! 'alat','bohr','angstrom'
    CHARACTER(len=80) :: cell_units
    !  alat: lattice parameter - often used to scale quantities, or
    !  in combination to other parameters/constants to define new units
    REAL(DP) :: alat = 0.0_DP
    ! omega: volume of the simulation cell
    REAl(DP) :: omega = 0.0_DP
    ! tpiba: 2 PI/alat, tpiba2=tpiba^2
    REAL(DP) :: tpiba  = 0.0_DP, tpiba2 = 0.0_DP
    !  direct and reciprocal lattice primitive vectors
    !  at(:,i) are the lattice vectors of the simulation cell, a_i,
    !          in alat units: a_i(:) = at(:,i)/alat
    !  bg(:,i) are the reciprocal lattice vectors, b_i,
    !          in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
    REAL(DP) :: at(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    REAL(DP) :: bg(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    !
    ! parameters for reference cell 
    REAL(DP) :: ref_tpiba2 = 0.0_DP
    REAL(DP) :: ref_at(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    REAL(DP) :: ref_bg(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    !
    ! parameter to store tpiba2 calculated from the input cell parameter 
    ! used in emass_preconditioning, required for restarting variable cell calculation correctly in CP
    REAL(DP) :: init_tpiba2 = 0.0_DP
    !
! -------------------------------------------------------------------------
! ...  periodicity box
! ...  In the matrix "a" every row is the vector of each side of 
! ...  the cell in the real space

        TYPE boxdimensions
          REAL(DP) :: a(3,3)    ! direct lattice generators
          REAL(DP) :: m1(3,3)   ! reciprocal lattice generators
          REAL(DP) :: omega     ! cell volume = determinant of a
          REAL(DP) :: g(3,3)    ! metric tensor
          REAL(DP) :: gvel(3,3) ! metric velocity
          REAL(DP) :: pail(3,3) ! stress tensor ( scaled coor. )
          REAL(DP) :: paiu(3,3) ! stress tensor ( cartesian coor. )
          REAL(DP) :: hmat(3,3) ! cell parameters ( transpose of "a" )
          REAL(DP) :: hvel(3,3) ! cell velocity
          REAL(DP) :: hinv(3,3)
          REAL(DP) :: deth
          INTEGER :: perd(3)
        END TYPE boxdimensions

        !  The following relations should always be kept valid:
        !     h = at*alat; ainv = h^(-1); ht=transpose(h)
        REAL(DP) :: h(3,3)    = 0.0_DP ! simulation cell at time t 
        REAL(DP) :: ainv(3,3) = 0.0_DP
        REAL(DP) :: hold(3,3) = 0.0_DP ! simulation cell at time t-delt
        REAL(DP) :: hnew(3,3) = 0.0_DP ! simulation cell at time t+delt
        REAL(DP) :: velh(3,3) = 0.0_DP ! simulation cell velocity
        REAL(DP) :: deth      = 0.0_DP ! determinant of h ( cell volume )

        INTEGER   :: iforceh(3,3) = 1  ! if iforceh( i, j ) = 0 then h( i, j ) 
                                       ! is not allowed to move
        LOGICAL   :: fix_volume = .FALSE.! True if cell volume is kept fixed
        LOGICAL   :: fix_area = .FALSE.  ! True if area in xy plane is kept constant
        LOGICAL   :: isotropic = .FALSE. ! True if volume option is chosen for cell_dofree 
        REAL(DP) :: wmass = 0.0_DP     ! cell fictitious mass
        REAL(DP) :: press = 0.0_DP     ! external pressure 

        REAL(DP) :: frich  = 0.0_DP    ! friction parameter for cell damped dynamics
        REAL(DP) :: greash = 1.0_DP    ! greas parameter for damped dynamics

        LOGICAL :: tcell_base_init = .FALSE.

        INTERFACE cell_init
          MODULE PROCEDURE cell_init_ht, cell_init_a
        END INTERFACE

        INTERFACE pbcs
          MODULE PROCEDURE pbcs_components, pbcs_vectors
        END INTERFACE

        INTERFACE s_to_r
          MODULE PROCEDURE s_to_r1, s_to_r1b, s_to_r3
        END INTERFACE

        INTERFACE r_to_s
          MODULE PROCEDURE r_to_s1, r_to_s1b, r_to_s3
        END INTERFACE
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!
!
  SUBROUTINE cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
               cosbc_, trd_ht, rd_ht, cell_units_ )
    !
    ! ... initialize cell_base module variables, set up crystal lattice
    !

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ibrav_
    REAL(DP), INTENT(IN) :: celldm_ (6)
    LOGICAL, INTENT(IN) :: trd_ht
    REAL(DP), INTENT(IN) :: rd_ht (3,3)
    CHARACTER(LEN=*), INTENT(IN) :: cell_units_
    REAL(DP), INTENT(IN) :: a_ , b_ , c_ , cosab_, cosac_, cosbc_

    REAL(DP) :: units
    !
    IF ( ibrav_ == 0 .and. .not. trd_ht ) THEN
       CALL errore('cell_base_init', 'ibrav=0: must read cell parameters', 1)
    ELSE IF ( ibrav_ /= 0 .and. trd_ht ) THEN
       CALL errore('cell_base_init', 'redundant data for cell parameters', 2)
    END IF
    !
    ibrav  = ibrav_
    celldm = celldm_
    a = a_ ; b = b_ ; c = c_ ; cosab = cosab_ ; cosac = cosac_ ; cosbc = cosbc_
    cell_units = cell_units_
    units = 0.0_DP
    !
    IF ( trd_ht ) THEN
      !
      ! ... crystal lattice vectors read from input: find units
      !
      SELECT CASE ( TRIM( cell_units ) )
        CASE ( 'bohr' )
          IF( celldm( 1 ) /= 0.0_DP .OR. a /= 0.0_dp ) CALL errore &
              ('cell_base_init','lattice parameter specified twice',1)
          units = 1.0_DP
        CASE ( 'angstrom' )
          IF( celldm( 1 ) /= 0.0_DP .OR. a /= 0.0_dp ) CALL errore &
              ('cell_base_init','lattice parameter specified twice',2)
          units = 1.0_DP / bohr_radius_angs
        CASE ( 'alat' ) 
          IF( celldm( 1 ) /= 0.0_DP ) THEN
             units = celldm( 1 )
          ELSE IF ( a /= 0.0_dp ) THEN
             units = a / bohr_radius_angs
          ELSE
            CALL errore ('cell_base_init', &
                         'lattice parameter not specified',1) 
          END IF
          ! following case is deprecated and should be removed
        CASE ( 'none' ) 
          ! cell_units is 'none' if nothing was specified
          IF( celldm( 1 ) /= 0.0_DP ) THEN
             units = celldm( 1 )
             cell_units = 'alat'
          ELSE IF ( a /= 0.0_dp ) THEN
             units = a / bohr_radius_angs
             cell_units = 'alat'
          ELSE
             units = 1.0_DP
             cell_units = 'bohr'
          END IF
          !
        CASE DEFAULT
          CALL errore ('cell_base_init', &
                       'unexpected cell_units '//TRIM(cell_units),1) 
     END SELECT
     !
     ! ... Beware the transpose operation between matrices ht and at!
     !
     at = TRANSPOSE( rd_ht ) * units
     !
     ! ... at is in atomic units: find alat, bring at to alat units, find omega
     !
     IF( celldm( 1 ) /= 0.0_DP ) THEN
        alat = celldm( 1 )
     ELSE IF ( a /= 0.0_dp ) THEN
        alat = a / bohr_radius_angs
     ELSE
        alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 )
     END IF
     ! for compatibility: celldm still used in phonon etc 
     celldm(1) = alat
     !
     at(:,:) = at(:,:) / alat
     CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
     !
  ELSE
  ! ... crystal lattice via celldm or crystallographica parameters
  !
     IF ( celldm(1) == 0.D0 .and. a /= 0.D0 ) THEN
        !
        celldm(1) = a / bohr_radius_angs
        celldm(2) = b / a
        celldm(3) = c / a
        IF ( (ABS(cosab) > 1.0_dp) .OR. (ABS(cosac) > 1.0_dp) .OR. &
             (ABS(cosbc) > 1.0_dp) ) CALL errore ('cell_base_init',&
                         'incorrect values for cosab, cosac, cosbc',1)
        !
        IF ( ibrav == 14 ) THEN
           !
           ! ... triclinic lattice
           !
           celldm(4) = cosbc
           celldm(5) = cosac
           celldm(6) = cosab
           !
        ELSE IF ( ibrav ==-12 ) THEN
           !
           ! ... monoclinic P lattice, unique axis b
           !
           celldm(5) = cosac
           !
        ELSE
           !
           ! ... trigonal and monoclinic lattices, unique axis c
           !
           celldm(4) = cosab
           !
        ENDIF
        !
     ELSE IF ( celldm(1) /= 0.D0 .and. a /= 0.D0 ) THEN
        !
        CALL errore( 'input', 'do not specify both celldm and a,b,c!', 1 )
        !
     END IF
     !
     ! ... generate at (in atomic units) from ibrav and celldm
     !
     CALL latgen( ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
     !
     ! ... define lattice constants alat, divide at by alat
     !
     alat = celldm(1)
     at(:,:) = at(:,:) / alat
     !
  END IF
  !
  ! ... Generate the reciprocal lattice vectors
  !
  CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
  tpiba  = 2.0_DP * pi / alat
  tpiba2 = tpiba * tpiba
  init_tpiba2 = tpiba2 ! BS : this is used in CPV/src/init_run.f90 
  RETURN
  !
  END SUBROUTINE cell_base_init
  !
  SUBROUTINE ref_cell_base_init( ref_alat, rd_ref_ht, ref_cell_units )
      !
      ! ... initialize cell_base module variables, set up crystal lattice
      !

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rd_ref_ht (3,3)
      REAL(DP), INTENT(INOUT) :: ref_alat
      CHARACTER(LEN=*), INTENT(IN) :: ref_cell_units

      REAL(DP) :: units, ref_omega
      !
      ! ... reference cell lattice vectors read from REF_CELL_PARAMETERS Card: find units
      !
      SELECT CASE ( TRIM( ref_cell_units ) )
      !
      CASE ( 'bohr' )
        units = 1.0_DP
      CASE ( 'angstrom' )
        units = 1.0_DP / bohr_radius_angs
      CASE DEFAULT
        IF( ref_alat .GT. 0.0_DP ) THEN
          units = ref_alat 
        ELSE
          CALL errore('ref_cell_base_init', 'ref_alat must be set to a positive value (in A.U.) in SYSTEM namelist', 1)
        END IF
        !
      END SELECT
      !
      ! ... Beware the transpose operation between matrices ht and at!
      !
      ref_at = TRANSPOSE( rd_ref_ht ) * units
      !
      ! ... ref_at is in atomic units: find ref_alat, bring ref_at to ref_alat units
      !
      ref_alat = SQRT ( ref_at(1,1)**2+ref_at(2,1)**2+ref_at(3,1)**2 )
      !
      ref_at(:,:) = ref_at(:,:) / ref_alat
      !
      ! ... Generate the reciprocal lattice vectors from the reference cell
      !
      CALL recips( ref_at(1,1), ref_at(1,2), ref_at(1,3), ref_bg(1,1), ref_bg(1,2), ref_bg(1,3) )
      !
      ref_tpiba2  = (2.0_DP * pi / ref_alat)**2
      !
      CALL volume( ref_alat, ref_at(1,1), ref_at(1,2), ref_at(1,3), ref_omega )
      !
      WRITE( stdout, * )
      WRITE( stdout, '(3X,"Reference Cell read from REF_CELL_PARAMETERS Card")' )
      WRITE( stdout, '(3X,"Reference Cell alat  =",F14.8,1X,"A.U.")' ) ref_alat
      WRITE( stdout, '(3X,"ref_cell_a1 = ",1X,3f14.8)' ) ref_at(:,1)*ref_alat
      WRITE( stdout, '(3X,"ref_cell_a2 = ",1X,3f14.8)' ) ref_at(:,2)*ref_alat
      WRITE( stdout, '(3X,"ref_cell_a3 = ",1X,3f14.8)' ) ref_at(:,3)*ref_alat
      WRITE( stdout, * )
      WRITE( stdout, '(3X,"ref_cell_b1 = ",1X,3f14.8)' ) ref_bg(:,1)/ref_alat
      WRITE( stdout, '(3X,"ref_cell_b2 = ",1X,3f14.8)' ) ref_bg(:,2)/ref_alat
      WRITE( stdout, '(3X,"ref_cell_b3 = ",1X,3f14.8)' ) ref_bg(:,3)/ref_alat
      WRITE( stdout, '(3X,"Reference Cell Volume",F16.8,1X,"A.U.")' ) ref_omega 
      !
      RETURN
      !
  END SUBROUTINE ref_cell_base_init
!------------------------------------------------------------------------------!
! ...     set box
! ...     box%m1(i,1) == b1(i)   COLUMN are B vectors
! ...     box%a(1,i)  == a1(i)   ROW are A vector
! ...     box%omega   == volume
! ...     box%g(i,j)  == metric tensor G
!------------------------------------------------------------------------------!

        SUBROUTINE cell_init_ht( what, box, hval )
          TYPE (boxdimensions) :: box
          REAL(DP),  INTENT(IN) :: hval(3,3)
          CHARACTER, INTENT(IN) :: what
            IF( what == 't' .OR. what == 'T' ) THEN
               !  hval == ht
               box%a = hval
               box%hmat = TRANSPOSE( hval )
            ELSE
               !  hval == hmat
               box%hmat = hval
               box%a = TRANSPOSE( hval )
            END IF
            CALL gethinv( box )
            box%g = MATMUL( box%a(:,:), box%hmat(:,:) )
            box%gvel = 0.0_DP
            box%hvel = 0.0_DP
            box%pail = 0.0_DP
            box%paiu = 0.0_DP
          RETURN
        END SUBROUTINE cell_init_ht
          
!------------------------------------------------------------------------------!

        SUBROUTINE cell_init_a( alat, at, box )
          TYPE (boxdimensions) :: box
          REAL(DP), INTENT(IN) :: alat, at(3,3)
          INTEGER :: i
            DO i=1,3
             ! this is HT: the rows are the lattice vectors
              box%a(1,i) = at(i,1)*alat
              box%a(2,i) = at(i,2)*alat
              box%a(3,i) = at(i,3)*alat
              ! this is H : the column are the lattice vectors
              box%hmat(i,1) = at(i,1)*alat
              box%hmat(i,2) = at(i,2)*alat
              box%hmat(i,3) = at(i,3)*alat
            END DO
            box%pail = 0.0_DP
            box%paiu = 0.0_DP
            box%hvel = 0.0_DP
            CALL gethinv(box)
            box%g    = MATMUL( box%a(:,:), box%hmat(:,:) )
            box%gvel = 0.0_DP
          RETURN
        END SUBROUTINE cell_init_a

!------------------------------------------------------------------------------!

        SUBROUTINE r_to_s1 (r,s,box)
          REAL(DP), intent(out) ::  S(3)
          REAL(DP), intent(in) :: R(3)
          type (boxdimensions), intent(in) :: box
          integer i,j
          DO I=1,3
            S(I) = 0.0_DP
            DO J=1,3
              S(I) = S(I) + R(J)*box%m1(J,I)
            END DO
          END DO
          RETURN
        END SUBROUTINE r_to_s1

!------------------------------------------------------------------------------!

        SUBROUTINE r_to_s3 ( r, s, na, nsp, hinv )
          REAL(DP), intent(out) ::  S(:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(DP), intent(in) :: R(:,:)
          REAL(DP), intent(in) :: hinv(:,:)    ! hinv = TRANSPOSE( box%m1 )
          integer :: i, j, ia, is, isa
          isa = 0
          DO is = 1, nsp
            DO ia = 1, na(is)
              isa = isa + 1
              DO I=1,3
                S(I,isa) = 0.0_DP
                DO J=1,3
                  S(I,isa) = S(I,isa) + R(J,isa)*hinv(i,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE r_to_s3

!------------------------------------------------------------------------------!

        SUBROUTINE r_to_s1b ( r, s, hinv )
          REAL(DP), intent(out) ::  S(:)
          REAL(DP), intent(in) :: R(:)
          REAL(DP), intent(in) :: hinv(:,:)    ! hinv = TRANSPOSE( box%m1 )
          integer :: i, j
          DO I=1,3
            S(I) = 0.0_DP
            DO J=1,3
              S(I) = S(I) + R(J)*hinv(i,j)
            END DO
          END DO
          RETURN
        END SUBROUTINE r_to_s1b


!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r1 (S,R,box)
          REAL(DP), intent(in) ::  S(3)
          REAL(DP), intent(out) :: R(3)
          type (boxdimensions), intent(in) :: box
          integer i,j
          DO I=1,3
            R(I) = 0.0_DP
            DO J=1,3
              R(I) = R(I) + S(J)*box%a(J,I)
            END DO
          END DO
          RETURN
        END SUBROUTINE s_to_r1

!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r1b (S,R,h)
          REAL(DP), intent(in) ::  S(3)
          REAL(DP), intent(out) :: R(3)
          REAL(DP), intent(in) :: h(:,:)    ! h = TRANSPOSE( box%a )
          integer i,j
          DO I=1,3
            R(I) = 0.0_DP
            DO J=1,3
              R(I) = R(I) + S(J)*h(I,j)
            END DO
          END DO
          RETURN
        END SUBROUTINE s_to_r1b

!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r3 ( S, R, na, nsp, h )
          REAL(DP), intent(in) ::  S(:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(DP), intent(out) :: R(:,:)
          REAL(DP), intent(in) :: h(:,:)    ! h = TRANSPOSE( box%a )
          integer :: i, j, ia, is, isa
          isa = 0
          DO is = 1, nsp
            DO ia = 1, na(is)
              isa = isa + 1
              DO I = 1, 3
                R(I,isa) = 0.0_DP
                DO J = 1, 3
                  R(I,isa) = R(I,isa) + S(J,isa) * h(I,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE s_to_r3


!
!------------------------------------------------------------------------------!
!

      SUBROUTINE gethinv(box)
        USE matrix_inversion
        IMPLICIT NONE
        TYPE (boxdimensions), INTENT (INOUT) :: box
        !
        CALL invmat( 3, box%a, box%m1, box%omega )
        box%deth = box%omega
        box%hinv = TRANSPOSE( box%m1 )
        !
        RETURN
      END SUBROUTINE gethinv


      FUNCTION get_volume( hmat )
         IMPLICIT NONE
         REAL(DP) :: get_volume
         REAL(DP) :: hmat( 3, 3 )
         get_volume = hmat(1,1)*(hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2)) + &
                      hmat(1,2)*(hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3)) + &
                      hmat(1,3)*(hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1))
         RETURN
      END FUNCTION get_volume
!
!------------------------------------------------------------------------------!
!
      FUNCTION pbc(rin,box,nl) RESULT (rout)
        IMPLICIT NONE
        TYPE (boxdimensions) :: box
        REAL (DP) :: rin(3)
        REAL (DP) :: rout(3), s(3)
        INTEGER, OPTIONAL :: nl(3)

        s = matmul(box%hinv(:,:),rin)
        s = s - box%perd*nint(s)
        rout = matmul(box%hmat(:,:),s)
        IF (present(nl)) THEN
          s = REAL( nl, DP )
          rout = rout + matmul(box%hmat(:,:),s)
        END IF
      END FUNCTION pbc

!
!------------------------------------------------------------------------------!
!
      FUNCTION saw(emaxpos,eopreg,x) RESULT (sawout)
        IMPLICIT NONE
        REAL(DP) :: emaxpos,eopreg,x
        REAL(DP) :: y, sawout, z
        
        z = x - emaxpos 
        y = z - floor(z)
        
        if (y.le.eopreg) then
        
            sawout = (0.5_DP - y/eopreg) * (1._DP-eopreg)
        
        else 
!
! I would use:   sawout = y - 0.5_DP * ( 1.0_DP + eopreg )
!
            sawout = (-0.5_DP + (y-eopreg)/(1._DP-eopreg)) * (1._DP-eopreg)
        
        end if
        
      END FUNCTION saw

!TB - start
!------------------------------------------------------------------------------!
!mopopla - add a potential of a monopole plane (kflag = .true.)
!          or the compensating background charge (kflag = .false.)
!          I split those in order to plot both independently
! cite PRB 89, 245406 (2014)
!
      FUNCTION mopopla(zmon,x,kflag) RESULT (mopoplaout)
        IMPLICIT NONE
        REAL(DP) :: zmon,x
        REAL(DP) :: mopoplaout, z
        LOGICAL  :: kflag

        DO ! is x within the cell?
          IF (x>1.0) x=x-1.0
          IF (x<0.0) x=x+1.0
          IF (x<=1.0.and.x>=0.0) EXIT
        ENDDO

        z = (x - zmon)

        !Monopole-plane
        ! if z < 0, we are below the plane
        !    z > 0, above
        !    z < -0.5, the potential is again the same as for z > 0
        !              in order to make it periodic
        !    z > 0.5, the same as for z < 0
        !
        IF (z<=-0.5) z=z+1
        IF (z>=0.5) z=z-1
        IF (z.LE.0) THEN
           IF (kflag) THEN
              mopoplaout = ( 1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ELSE
           IF (kflag) THEN
              mopoplaout = ( -1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ENDIF

      END FUNCTION mopopla

!TB - end

!
!------------------------------------------------------------------------------!
!
          SUBROUTINE get_cell_param(box,cell,ang)
          IMPLICIT NONE
          TYPE(boxdimensions), INTENT(in) :: box
          REAL(DP), INTENT(out), DIMENSION(3) :: cell
          REAL(DP), INTENT(out), DIMENSION(3), OPTIONAL :: ang
! This code gets the cell parameters given the h-matrix:
! a
          cell(1)=sqrt(box%hmat(1,1)*box%hmat(1,1)+box%hmat(2,1)*box%hmat(2,1) &
                      +box%hmat(3,1)*box%hmat(3,1))
! b
          cell(2)=sqrt(box%hmat(1,2)*box%hmat(1,2)+box%hmat(2,2)*box%hmat(2,2) &
                      +box%hmat(3,2)*box%hmat(3,2))
! c
          cell(3)=sqrt(box%hmat(1,3)*box%hmat(1,3)+box%hmat(2,3)*box%hmat(2,3) &
                      +box%hmat(3,3)*box%hmat(3,3))
          IF (PRESENT(ang)) THEN
! gamma
             ang(1)=acos((box%hmat(1,1)*box%hmat(1,2)+ &
                          box%hmat(2,1)*box%hmat(2,2) &
                      +box%hmat(3,1)*box%hmat(3,2))/(cell(1)*cell(2)))
! beta
             ang(2)=acos((box%hmat(1,1)*box%hmat(1,3)+ &
                          box%hmat(2,1)*box%hmat(2,3) &
                      +box%hmat(3,1)*box%hmat(3,3))/(cell(1)*cell(3)))
! alpha
           ang(3)=acos((box%hmat(1,2)*box%hmat(1,3)+ &
                        box%hmat(2,2)*box%hmat(2,3) &
                      +box%hmat(3,2)*box%hmat(3,3))/(cell(2)*cell(3)))
!           ang=ang*180.0_DP/pi

          ENDIF
          END SUBROUTINE get_cell_param

!------------------------------------------------------------------------------!

      SUBROUTINE pbcs_components(x1, y1, z1, x2, y2, z2, m)
! ... This subroutine compute the periodic boundary conditions in the scaled
! ... variables system
        USE kinds
        INTEGER, INTENT(IN)  :: M
        REAL(DP),  INTENT(IN)  :: X1,Y1,Z1
        REAL(DP),  INTENT(OUT) :: X2,Y2,Z2
        REAL(DP) MIC
        MIC = REAL( M, DP )
        X2 = X1 - DNINT(X1/MIC)*MIC
        Y2 = Y1 - DNINT(Y1/MIC)*MIC
        Z2 = Z1 - DNINT(Z1/MIC)*MIC
        RETURN
      END SUBROUTINE pbcs_components

!------------------------------------------------------------------------------!

      SUBROUTINE pbcs_vectors(v, w, m)
! ... This subroutine compute the periodic boundary conditions in the scaled
! ... variables system
        USE kinds
        INTEGER, INTENT(IN)  :: m
        REAL(DP),  INTENT(IN)  :: v(3)
        REAL(DP),  INTENT(OUT) :: w(3)
        REAL(DP) :: MIC
        MIC = REAL( M, DP )
        w(1) = v(1) - DNINT(v(1)/MIC)*MIC
        w(2) = v(2) - DNINT(v(2)/MIC)*MIC
        w(3) = v(3) - DNINT(v(3)/MIC)*MIC
        RETURN
      END SUBROUTINE pbcs_vectors

!------------------------------------------------------------------------------!

  SUBROUTINE set_h_ainv()
    !
    ! CP-PW compatibility: align CP arrays H and ainv to at and bg
    !
    IMPLICIT NONE
    !
    !write(stdout,*) 'alat=',alat
    !write(stdout,*) 'at=',at
    !write(stdout,*) 'bg=',bg
    !
    h(:,:) = at(:,:)*alat
    !
    ainv(1,:) = bg(:,1)/alat
    ainv(2,:) = bg(:,2)/alat
    ainv(3,:) = bg(:,3)/alat
    !
  END SUBROUTINE set_h_ainv

!------------------------------------------------------------------------------!

  SUBROUTINE cell_dyn_init( trd_ht, rd_ht, wc_ , total_ions_mass , press_ , &
               frich_ , greash_ , cell_dofree )

    USE constants, ONLY: au_gpa, amu_au
    USE io_global, ONLY: stdout

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: cell_dofree
    LOGICAL, INTENT(IN) :: trd_ht
    REAL(DP), INTENT(IN) :: rd_ht (3,3)
    REAL(DP),  INTENT(IN) :: wc_ , frich_ , greash_ , total_ions_mass
    REAL(DP),  INTENT(IN) :: press_ ! external pressure from input 
                                    ! ( in KBar = 0.1 GPa )
    INTEGER   :: j
    !
    press  = press_ / 10.0_DP ! convert press in KBar to GPa
    press  = press  / au_gpa  ! convert to AU
    !  frich  = frich_   ! for the time being this is set elsewhere
    greash = greash_

    WRITE( stdout, 105 )
    WRITE( stdout, 110 ) press_
105 format(/,3X,'Simulation Cell Parameters (from input)')
110 format(  3X,'external pressure       = ',f15.2,' [KBar]')

    wmass  = wc_
    IF( wmass == 0.0_DP ) THEN
      wmass = 3.0_DP / (4.0_DP * pi**2 ) * total_ions_mass
      wmass = wmass * AMU_AU
      WRITE( stdout,130) wmass
    ELSE
      WRITE( stdout,120) wmass
    END IF
120 format(3X,'wmass (read from input) = ',f15.2,' [AU]')
130 format(3X,'wmass (calculated)      = ',f15.2,' [AU]')

    IF( wmass <= 0.0_DP ) &
      CALL errore(' cell_dyn_init',' wmass out of range ',0)

    IF ( trd_ht ) THEN
      !
      WRITE( stdout, 210 )
      WRITE( stdout, 220 ) ( rd_ht( 1, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( rd_ht( 2, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( rd_ht( 3, j ), j = 1, 3 )
      !
210   format(3X,'initial cell from CELL_PARAMETERS card')
220   format(3X,3F14.8)
      !
    END IF
    !
    ainv(1,:) = bg(:,1)/alat
    ainv(2,:) = bg(:,2)/alat
    ainv(3,:) = bg(:,3)/alat
    !
    CALL init_dofree ( cell_dofree ) 
    !
    tcell_base_init = .TRUE.

    WRITE( stdout, 300 ) ibrav
    WRITE( stdout, 305 ) alat
    WRITE( stdout, 310 ) at(:,1)*alat
    WRITE( stdout, 320 ) at(:,2)*alat
    WRITE( stdout, 330 ) at(:,3)*alat
    WRITE( stdout, *   )
    WRITE( stdout, 350 ) bg(:,1)/alat
    WRITE( stdout, 360 ) bg(:,2)/alat
    WRITE( stdout, 370 ) bg(:,3)/alat
    WRITE( stdout, 340 ) omega
300 FORMAT( 3X, 'ibrav = ',I4)
305 FORMAT( 3X, 'alat  = ',F14.8)
310 FORMAT( 3X, 'a1    = ',3F14.8)
320 FORMAT( 3X, 'a2    = ',3F14.8)
330 FORMAT( 3X, 'a3    = ',3F14.8)
350 FORMAT( 3X, 'b1    = ',3F14.8)
360 FORMAT( 3X, 'b2    = ',3F14.8)
370 FORMAT( 3X, 'b3    = ',3F14.8)
340 FORMAT( 3X, 'omega = ',F16.8)


    RETURN
  END SUBROUTINE cell_dyn_init

!------------------------------------------------------------------------------!
  SUBROUTINE init_dofree ( cell_dofree ) 

     ! set constraints on cell dynamics/optimization

     CHARACTER(LEN=*), INTENT(IN) :: cell_dofree

     SELECT CASE ( TRIM( cell_dofree ) )

            CASE ( 'all', 'default' )
              iforceh = 1
            CASE ( 'shape' )
              iforceh = 1
              fix_volume = .true.
! 2DSHAPE: CASE FOR SHAPE CHANGE IN xy PLANE WITH CONST AREA
! contribution from Richard Charles Andrew
! Physics Department, University of Pretoria
! South Africa, august 2012.
            CASE ( '2Dshape' )
              iforceh = 1
              iforceh(3,3) = 0
              iforceh(1,3) = 0
              iforceh(3,1) = 0
              iforceh(2,3) = 0
              iforceh(3,2) = 0
              fix_area = .true.
! 2DSHAPE
            CASE ( 'volume' )
              !CALL errore(' init_dofree ', &
              !   ' cell_dofree = '//TRIM(cell_dofree)//' not yet implemented ', 1 )
              IF ( ibrav /= 1 ) THEN
                CALL errore('cell_dofree', 'Isotropic expansion is only allowed for ibrav=1; i.e. for simple cubic', 1)
              END IF
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
              iforceh(3,3) = 1
              isotropic    = .TRUE.
            CASE ('x')
              iforceh      = 0
              iforceh(1,1) = 1
            CASE ('y')
              iforceh      = 0
              iforceh(2,2) = 1
            CASE ('z')
              iforceh      = 0
              iforceh(3,3) = 1
            CASE ('xy')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
              ! ... if you want the entire xy plane to be free, uncomment:
              ! iforceh(1,2) = 1
              ! iforceh(2,1) = 1
! 2DSHAPE THE ENTIRE xy PLANE IS FREE
            CASE ('2Dxy')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
              iforceh(1,2) = 1
              iforceh(2,1) = 1
! 2DSHAPE
            CASE ('xz')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(3,3) = 1
            CASE ('yz')
              iforceh      = 0
              iforceh(2,2) = 1
              iforceh(3,3) = 1
            CASE ('xyz')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
              iforceh(3,3) = 1
            CASE DEFAULT
              CALL errore(' init_dofree ',' unknown cell_dofree '//TRIM(cell_dofree), 1 )

     END SELECT
  END SUBROUTINE init_dofree 
      
!------------------------------------------------------------------------------!

  SUBROUTINE cell_base_reinit( ht )

    USE control_flags, ONLY: iverbosity

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ht (3,3)

    INTEGER   :: j

    alat   =  sqrt( ht(1,1)*ht(1,1) + ht(1,2)*ht(1,2) + ht(1,3)*ht(1,3) )
    tpiba  = 2.0_DP * pi / alat
    tpiba2 = tpiba * tpiba
    !
    IF( iverbosity > 2 ) THEN
      WRITE( stdout, 210 )
      WRITE( stdout, 220 ) ( ht( 1, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( ht( 2, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( ht( 3, j ), j = 1, 3 )
    END IF
210 format(3X,'Simulation cell parameters with the new cell:')
220 format(3X,3F14.8)

    !    matrix "ht" used in CP is the transpose of matrix "at"
    !    times the lattice parameter "alat"; matrix "ainv" is "bg" divided alat
    !
    at     = TRANSPOSE( ht ) / alat
    !
    CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
    CALL volume( alat, at(1,1), at(1,2), at(1,3), deth )
    omega = deth
    !
    ainv(1,:) = bg(:,1)/alat
    ainv(2,:) = bg(:,2)/alat
    ainv(3,:) = bg(:,3)/alat
    !
    IF( iverbosity > 2 ) THEN
      WRITE( stdout, 305 ) alat
      WRITE( stdout, 310 ) at(:,1)*alat
      WRITE( stdout, 320 ) at(:,2)*alat
      WRITE( stdout, 330 ) at(:,3)*alat
      WRITE( stdout, *   )
      WRITE( stdout, 350 ) bg(:,1)/alat
      WRITE( stdout, 360 ) bg(:,2)/alat
      WRITE( stdout, 370 ) bg(:,3)/alat
      WRITE( stdout, 340 ) omega
    END IF

305 FORMAT( 3X, 'alat  = ',F14.8)
310 FORMAT( 3X, 'a1    = ',3F14.8)
320 FORMAT( 3X, 'a2    = ',3F14.8)
330 FORMAT( 3X, 'a3    = ',3F14.8)
350 FORMAT( 3X, 'b1    = ',3F14.8)
360 FORMAT( 3X, 'b2    = ',3F14.8)
370 FORMAT( 3X, 'b3    = ',3F14.8)
340 FORMAT( 3X, 'omega = ',F14.8)


    RETURN
  END SUBROUTINE cell_base_reinit



!------------------------------------------------------------------------------!

  SUBROUTINE cell_steepest( hnew, h, delt, iforceh, fcell )
    REAL(DP), INTENT(OUT) :: hnew(3,3)
    REAL(DP), INTENT(IN) :: h(3,3), fcell(3,3)
    INTEGER,      INTENT(IN) :: iforceh(3,3)
    REAL(DP), INTENT(IN) :: delt
    INTEGER      :: i, j
    REAL(DP) :: dt2,fiso
    dt2 = delt * delt
    !
    IF( isotropic ) THEN
      !
      ! Isotropic force on the cell
      !
      fiso = (fcell(1,1)+fcell(2,2)+fcell(3,3))/3.0_DP
      !
      DO j=1,3
        DO i=1,3
          hnew(i,j) = h(i,j) + dt2 * fiso * REAL( iforceh(i,j), DP )
        ENDDO
      ENDDO
      !
    ELSE
      !
      DO j=1,3
        DO i=1,3
          hnew(i,j) = h(i,j) + dt2 * fcell(i,j) * REAL( iforceh(i,j), DP )
        ENDDO
      ENDDO
      !
    END IF
    !
    RETURN
  END SUBROUTINE cell_steepest


!------------------------------------------------------------------------------!

  SUBROUTINE cell_verlet( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, hnos )
    REAL(DP), INTENT(OUT) :: hnew(3,3)
    REAL(DP), INTENT(IN) :: h(3,3), hold(3,3), hnos(3,3), fcell(3,3)
    INTEGER,      INTENT(IN) :: iforceh(3,3)
    REAL(DP), INTENT(IN) :: frich, delt
    LOGICAL,      INTENT(IN) :: tnoseh

    REAL(DP) :: htmp(3,3)
    REAL(DP) :: verl1, verl2, verl3, dt2, ftmp, v1, v2, v3, fiso
    INTEGER      :: i, j
  
    dt2 = delt * delt

    IF( tnoseh ) THEN
      ftmp = 0.0_DP
      htmp = hnos
    ELSE
      ftmp = frich
      htmp = 0.0_DP
    END IF

    verl1 = 2.0_DP / ( 1.0_DP + ftmp ) 
    verl2 = 1.0_DP - verl1
    verl3 = dt2 / ( 1.0_DP + ftmp )
    verl1 = verl1 - 1.0_DP

    IF( isotropic ) THEN
      !
      fiso = (fcell(1,1)+fcell(2,2)+fcell(3,3))/3.0_DP
      !
      DO j=1,3
        DO i=1,3
          v1 = verl1 * h(i,j)
          v2 = verl2 * hold(i,j)
          v3 = verl3 * ( fiso - htmp(i,j) )
          hnew(i,j) = h(i,j) + ( v1 + v2 + v3 ) * REAL( iforceh(i,j), DP )
        ENDDO
      ENDDO
      !
    ELSE
      !
      DO j=1,3
        DO i=1,3
          v1 = verl1 * h(i,j)
          v2 = verl2 * hold(i,j)
          v3 = verl3 * ( fcell(i,j) - htmp(i,j) )
          hnew(i,j) = h(i,j) + ( v1 + v2 + v3 ) * REAL( iforceh(i,j), DP )
        ENDDO
      ENDDO
      !
    END IF
  
    RETURN
  END SUBROUTINE cell_verlet

!------------------------------------------------------------------------------!

  subroutine cell_hmove( h, hold, delt, iforceh, fcell )
    REAL(DP), intent(out) :: h(3,3)
    REAL(DP), intent(in) :: hold(3,3), fcell(3,3)
    REAL(DP), intent(in) :: delt
    integer, intent(in) :: iforceh(3,3)
    REAL(DP) :: dt2by2, fac
    integer :: i, j
    dt2by2 = 0.5_DP * delt * delt
    fac = dt2by2
    do i=1,3
      do j=1,3
        h(i,j) = hold(i,j) + fac * iforceh(i,j) * fcell(i,j)
      end do
    end do
    return
  end subroutine cell_hmove

!------------------------------------------------------------------------------!

  subroutine cell_force( fcell, ainv, stress, omega, press, wmassIN )
    USE constants, ONLY : eps8
    REAL(DP), intent(out) :: fcell(3,3)
    REAL(DP), intent(in) :: stress(3,3), ainv(3,3)
    REAL(DP), intent(in) :: omega, press
    REAL(DP), intent(in), optional :: wmassIN
    integer        :: i, j
    REAL(DP) :: wmass
    IF (.not. present(wmassIN)) THEN
      wmass = 1.0
    ELSE
      wmass = wmassIN
    END IF
    do j=1,3
      do i=1,3
        fcell(i,j) = ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) + ainv(j,3)*stress(i,3)
      end do
    end do
    do j=1,3
      do i=1,3
        fcell(i,j) = fcell(i,j) - ainv(j,i) * press
      end do
    end do
    IF( wmass < eps8 ) &
       CALL errore( ' movecell ',' cell mass is less than 0 ! ', 1 )
    fcell = omega * fcell / wmass
    return
  end subroutine cell_force

!------------------------------------------------------------------------------!

  subroutine cell_move( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )
    REAL(DP), intent(out) :: hnew(3,3)
    REAL(DP), intent(in) :: h(3,3), hold(3,3), fcell(3,3)
    REAL(DP), intent(in) :: vnhh(3,3), velh(3,3)
    integer,  intent(in) :: iforceh(3,3)
    REAL(DP), intent(in) :: frich, delt
    logical,  intent(in) :: tnoseh, tsdc

    REAL(DP) :: hnos(3,3)

    hnew = 0.0

    if( tnoseh ) then
      hnos = vnhh * velh
    else
      hnos = 0.0_DP
    end if
!
    IF( tsdc ) THEN
      call cell_steepest( hnew, h, delt, iforceh, fcell )
    ELSE
      call cell_verlet( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, hnos )
    END IF

    return
  end subroutine cell_move

!------------------------------------------------------------------------------!

  SUBROUTINE cell_gamma( hgamma, ainv, h, velh )
    !
    ! Compute hgamma = g^-1 * dg/dt
    ! that enters in the ions equation of motion
    !
    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: hgamma(3,3)
    REAL(DP), INTENT(IN)  :: ainv(3,3), h(3,3), velh(3,3)
    REAL(DP) :: gm1(3,3), gdot(3,3)
    ! 
    !  g^-1 inverse of metric tensor = (ht*h)^-1 = ht^-1 * h^-1
    !
    gm1    = MATMUL( ainv, TRANSPOSE( ainv ) )  
    ! 
    !  dg/dt = d(ht*h)/dt = dht/dt*h + ht*dh/dt ! derivative of metrix tensor
    !
    gdot   = MATMUL( TRANSPOSE( velh ), h ) + MATMUL( TRANSPOSE( h ), velh )
    !
    hgamma = MATMUL( gm1, gdot )
    !
    RETURN
  END SUBROUTINE cell_gamma

!------------------------------------------------------------------------------!

  SUBROUTINE cell_update_vel( htp, ht0, htm, delt, velh )
    !
    IMPLICIT NONE
    TYPE (boxdimensions)  :: htp, ht0, htm
    REAL(DP), INTENT(IN)  :: delt
    REAL(DP), INTENT(OUT) :: velh( 3, 3 )

    velh(:,:) = ( htp%hmat(:,:) - htm%hmat(:,:) ) / ( 2.0d0 * delt )
    htp%gvel = ( htp%g(:,:) - htm%g(:,:) ) / ( 2.0d0 * delt )
    ht0%hvel = velh

    RETURN
  END SUBROUTINE cell_update_vel


!------------------------------------------------------------------------------!

  subroutine cell_kinene( ekinh, temphh, velh )
    use constants, only: k_boltzmann_au
    implicit none
    REAL(DP), intent(out) :: ekinh, temphh(3,3)
    REAL(DP), intent(in)  :: velh(3,3)
    integer :: i,j
    ekinh = 0.0_DP
    do j=1,3
       do i=1,3
          ekinh = ekinh + 0.5_DP*wmass*velh(i,j)*velh(i,j)
          temphh(i,j) = wmass*velh(i,j)*velh(i,j)/k_boltzmann_au
       end do
    end do
    return
  end subroutine cell_kinene

!------------------------------------------------------------------------------!

  function cell_alat( )
    real(DP) :: cell_alat
    if( .NOT. tcell_base_init ) &
      call errore( ' cell_alat ', ' alat has not been set ', 1 )
    cell_alat = alat
    return 
  end function cell_alat
!
!------------------------------------------------------------------------------!
   END MODULE cell_base
!------------------------------------------------------------------------------!
