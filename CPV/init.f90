!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Oct 31 12:18:54 MET 1999
!  ----------------------------------------------

     MODULE init_fpmd

       USE kinds
       USE io_global, ONLY: ionode, stdout

       IMPLICIT NONE
       PRIVATE       
       SAVE

       REAL(dbl) :: cclock
       EXTERNAL  :: cclock

       PUBLIC :: init1s, init0s

     CONTAINS

!  BEGIN manual

      SUBROUTINE init0s(gv, kp, ps, atoms_m, atoms_0, atoms_p, wfill, &
        wempt, ht_m2, ht_m, ht, fnl, eigr, nspin)

!  this routine handles data initialization
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE parameters, ONLY: nspinx
      USE phase_factors_module, ONLY: strucf
      USE nose_ions, ONLY: nosepinit
      USE nose_electrons, ONLY: noseeinit
      USE cp_types
      USE atoms_type_module, ONLY: atoms_type
      USE time_step, ONLY: delt
      USE cell_module, ONLY: cell_init, get_lattice_vectors
      USE cell_base, ONLY: omega, alat
      USE electrons_module, ONLY: electron_mass_init, band_init, n_emp
      USE electrons_base, ONLY: nupdwn
      use small_box, only: a1b, a2b, a3b, omegab, ainvb, tpibab, small_box_set
      use small_box, only: alatb, b1b, b2b, b3b
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrbx
      USE reciprocal_space_mesh, ONLY:  newg
      USE reciprocal_vectors, ONLY:  ngwt, gstart, gzero, ngm, ngmt, ngw
      USE pseudopotential, ONLY: formf, nsanl, ngh, pseudopotential_init
      USE ions_module, ONLY: atoms_init
      USE ions_base, ONLY: nsp, na, nat
      USE pseudo_projector, ONLY: allocate_projector, projector
      USE diis, ONLY: allocate_diis, delt_diis
      USE cell_module, only: boxdimensions
      USE mp_global, ONLY: mpime, root
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor, wave_descriptor_init, &
        wave_descriptor_info
      USE descriptors_module, ONLY: get_local_dims, get_global_dims
      USE control_flags, ONLY: nbeg, prn, tbeg, timing, t_diis
      USE input_parameters, ONLY: rd_ht
      USE turbo, ONLY: tturbo, init_turbo
      USE ions_base, ONLY: tau_srt, tau_units, ind_srt, if_pos, atm
      USE gvecp, ONLY: ecutp
      USE gvecb, ONLY: gcutb
      USE stick, ONLY: dfftp
      USE grid_dimensions, ONLY: nr1, nr2, nr3

      IMPLICIT NONE

! ... declare subroutine arguments

      TYPE (atoms_type)   :: atoms_0, atoms_p, atoms_m
      TYPE (wave_descriptor) :: wfill, wempt
      TYPE (pseudo) :: ps
      TYPE (projector) :: fnl(:,:)
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv
      TYPE (kpoints) :: kp
      TYPE (boxdimensions) :: ht_m2, ht_m, ht
      INTEGER :: nspin

! ... declare other variables
      REAL(dbl) :: s1, s2, s3, s4, s5
      REAL(dbl) :: annee
      REAL(dbl) :: a1(3), a2(3), a3(3)
      INTEGER :: i, ispin, isym
      LOGICAL :: tk
      COMPLEX(dbl), ALLOCATABLE :: sfac(:,:)
      REAL(dbl) :: rat1, rat2, rat3
      INTEGER :: neupdwn( nspinx )


!  end of declarations
!  ----------------------------------------------

      s1 = cclock()

! ... report information
      call get_lattice_vectors( a1, a2, a3 )

      rat1 = DBLE( nr1b ) / DBLE( nr1 )
      rat2 = DBLE( nr2b ) / DBLE( nr2 )
      rat3 = DBLE( nr3b ) / DBLE( nr3 )
      CALL small_box_set( alat, omega, a1, a2, a3, rat1, rat2, rat3 )

      gcutb = ecutp / tpibab / tpibab
      CALL ggenb (b1b, b2b, b3b, nr1b ,nr2b, nr3b, nr1bx ,nr2bx, nr3bx, gcutb )  


! ... arrange for reciprocal lattice vectors
      tk = .NOT. ( kp%scheme == 'gamma' )
      CALL allocate_recvecs(gv, ngm, ngmt, ngw, ngwt, tk, kp%nkpt)

! ... count reciprocal lattice vectors (generated from celldm)
      CALL ggenew(gv, kp, prn)

! ... Allocate + Initialize pseudopotentials
      CALL pseudopotential_init(ps, na, nsp, gv, kp)

! ... Arrange for phase factors exp(i G dot r)
      CALL allocate_phfac(eigr, nr1, nr2, nr3, nsp, nat, ngw, ngm)
      ALLOCATE( sfac( nsp, ngm ) )

      s2 = cclock()

! ... Initialize simulation cell
! ... if tbeg = .TRUE. cell parameters have been specified in the input file
! ... (see card 'TBEG')
! ... But the g-space grid is genereted according to celldm

      IF( tbeg ) THEN
        CALL cell_init( ht, rd_ht )
        CALL cell_init( ht_m, rd_ht )
        CALL cell_init( ht_m2, rd_ht )
      ELSE
        CALL cell_init( ht, a1, a2, a3 )
        CALL cell_init( ht_m, a1, a2, a3 )
        CALL cell_init( ht_m2, a1, a2, a3 )
      END IF

! ... initialize atomic configuration (should be called after metric_init)
      CALL atoms_init( atoms_m, atoms_0, atoms_p, tau_srt, ind_srt, if_pos, atm, tau_units, alat, ht )

      CALL nosepinit( REAL(atoms_0%doft, dbl) )

! ... compute reciprocal lattice vectors
      CALL newg(gv, kp, ht%m1)

      s3 = cclock()

! ... compute structure factors
      CALL strucf(sfac, atoms_0, eigr, gv)

      s4 = cclock()

! ... compute local form factors
      CALL formf(ht, gv, kp, ps)

      s5 = cclock()

      IF(ionode) THEN
        WRITE( stdout,'(/,"   ggen   (sec) : ",F8.3)') (s2-s1)
        WRITE( stdout,'(  "   newg   (sec) : ",F8.3)') (s3-s2)
        WRITE( stdout,'(  "   strucf (sec) : ",F8.3)') (s4-s3)
        WRITE( stdout,'(  "   formf  (sec) : ",F8.3)') (s5-s4)
      END IF

      isym = 0
      IF( tk ) isym = 1

      !  empty states, always same number of spin up and down states
      neupdwn( 1:nspin ) = n_emp

      CALL wave_descriptor_init( wfill, gv%ngw_l, gv%ngw_g, nupdwn,  nupdwn, &
             kp%nkpt, kp%nkpt, nspin, isym, gv%gzero )
      CALL wave_descriptor_init( wempt, gv%ngw_l, gv%ngw_g, neupdwn, neupdwn, &
             kp%nkpt, kp%nkpt, nspin, isym, gv%gzero )

      IF( prn ) THEN
        CALL wave_descriptor_info( wfill, 'wfill', stdout )
        CALL wave_descriptor_info( wempt, 'wempt', stdout )
      END IF

! ... if tturbo=.TRUE. some data is stored in memory instead of being
! ... recalculated (see card 'TURBO')
      IF( tturbo ) THEN
        CALL init_turbo( dfftp%nr1x, dfftp%nr2x, dfftp%npl )
      ENDIF


      CALL cpflush  ! flush output streams
      DEALLOCATE( sfac )

      RETURN
      END SUBROUTINE


!  BEGIN manual

    SUBROUTINE init1s(gv, kp, ps, atoms_m, atoms_0, atoms_p, cm, c0, wfill, &
      ce, wempt, ht_m2, ht_m, ht, fnl, eigr, occ)

!  this routine handles data initialization
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE phase_factors_module, ONLY: strucf
      USE wave_init, ONLY: pw_atomic_init
      USE nose_ions, ONLY: nosepinit
      USE nose_electrons, ONLY: noseeinit
      USE cp_types
      USE atoms_type_module, ONLY: atoms_type
      USE time_step, ONLY: delt
      USE cell_module, ONLY: cell_init, get_lattice_vectors, alat
      USE electrons_module, ONLY: electron_mass_init, band_init, nbnd
      USE reciprocal_space_mesh, ONLY:  newg
      USE reciprocal_vectors, ONLY:  ngwt, gstart, gzero, ngm, ngmt, ngw
      USE pseudopotential, ONLY: formf, nsanl, ngh, pseudopotential_init
      USE ions_module, ONLY: atoms_init
      USE ions_base, ONLY: nsp, na, nat
      USE pseudo_projector, ONLY: allocate_projector, projector
      USE diis, ONLY: allocate_diis, delt_diis
      USE cell_module, only: boxdimensions
      USE mp_global, ONLY: mpime, root
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE descriptors_module, ONLY: get_local_dims, get_global_dims
      USE control_flags, ONLY: nbeg, prn, tbeg, timing, t_diis

      IMPLICIT NONE

! ... declare subroutine arguments

      TYPE (atoms_type)   :: atoms_0, atoms_p, atoms_m
      COMPLEX(dbl) :: cm(:,:,:,:), c0(:,:,:,:), ce(:,:,:,:)
      TYPE (wave_descriptor) :: wfill, wempt
      TYPE (pseudo) :: ps
      REAL(dbl) :: occ(:,:,:)
      TYPE (projector) :: fnl(:,:)
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv
      TYPE (kpoints) :: kp
      TYPE (boxdimensions) :: ht_m2, ht_m, ht

! ... declare other variables
      REAL(dbl) :: s1, s2, s3, s4, s5
      REAL(dbl) :: annee
      REAL(dbl) :: a1(3), a2(3), a3(3)
      INTEGER :: i
      LOGICAL :: tk
      REAL(dbl),  ALLOCATABLE :: hg_g(:)   ! squared length
      INTEGER,    ALLOCATABLE :: mill(:,:) ! Miller index, axis x, y, z
      COMPLEX(dbl), ALLOCATABLE :: sfac(:,:)


!  end of declarations
!  ----------------------------------------------

      s1 = cclock()

! ... initialize bands
      CALL band_init( occ )

! ...   initialize wave functions
      CALL pw_atomic_init(nbeg, cm, c0, wfill, ce, wempt, gv, kp, eigr%xyz)

! ... initialize the electronic fictitious mass and time step for 
! ... electron dynamics. Note that for 'diis' the electronic time step is
! ... different from that of ions (delt), this is because in the diis
! ... the time step for electrons is simply a convergence parameter.
      IF (t_diis) THEN
        CALL electron_mass_init(alat, gv%hg_l, gv%ngw_l)
      ELSE
        CALL electron_mass_init(alat, gv%hg_l, gv%ngw_l)
      END IF

! ... initialize nonlocal pseudopotentials coefficients
      CALL allocate_projector(fnl, nsanl, nbnd, ngh, kp%gamma_only)

      annee = 0.0d0
      CALL noseeinit( annee )

      IF(t_diis) THEN
! ...   arrange for DIIS minimization
        CALL allocate_diis(gv%ngw_l, nbnd, kp%nkpt)
      END IF

      CALL cpflush  ! flush output streams

      RETURN
    END SUBROUTINE


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Oct  8 14:58:38 MDT; 1999
!  ----------------------------------------------
!  SUBROUTINE ggenew(ht_0,gv,nr1,nr2,nr3,prn)
!  SUBROUTINE compute_hg(i,j,k,b1,b2,b3,nhg,ngw, &
!                        nr1,nr2,nr3)
!  ----------------------------------------------
!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE ggenew(gv, kp, prn)

!  this routine generates reciprocal lattice vectors with length squared
!  less than gcut, and sorts them in order of increasing length
!  their Fast Fourier Transform indices n1,n2,n3 are stored
!  when using only the Gamma point, only the positive half space is
!  spanned and indices for -G are provided
!  N.B: THE G'S ARE IN UNITS OF 2PI/A.
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE fft,       ONLY: fft_setup
      USE reciprocal_space_mesh, ONLY : gindexset
      USE reciprocal_vectors, ONLY : mill_g, g2_g
      USE cp_types
      USE cell_module, only: boxdimensions
      USE cell_module, ONLY : get_lattice_vectors, alat
      USE brillouin, ONLY: kpoints
      USE gvecp, ONLY:  gcutp
      USE gvecw, ONLY:  gkcut
      USE gvecs, ONLY:  gcuts
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE smooth_grid_dimensions, ONLY: nr1s,  nr2s,  nr3s
      USE control_flags, ONLY: tdipole
      USE berry_phase, ONLY: berry_setup

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: prn

      TYPE (recvecs), INTENT(OUT) :: gv
      TYPE (kpoints), INTENT(IN) :: kp

! ... declare other variables
      REAL(dbl) :: b1(3), b2(3), b3(3)
      REAL(dbl) :: a1(3), a2(3), a3(3)

      INTEGER   :: i

      REAL(dbl) :: s1, s2, s3, s4, s5, s6


! ... end of declarations
!  ----------------------------------------------

! ... reciprocal lattice generators
      call get_lattice_vectors( a1, a2, a3 )
      call recips(a1, a2, a3, b1, b2, b3 )

      !  Normally if a1, a2 and a3 are in cartesian coordinates
      !  and in a.u. units the corresponding bs are in cartesian
      !  coordinate too and in unit of 2 PI / a.u.
      !  now bring b1, b2 and b3 in units of 2 PI / alat

      b1 = b1 * alat 
      b2 = b2 * alat
      b3 = b3 * alat

! ... b1(i) = alat * ht_0%m1(i,1)

      s1 = cclock()

      CALL ggencp ( b1, b2, b3, nr1, nr2, nr3, nr1s, nr2s, nr3s,               &
     &      gcutp, gcuts, gkcut, kp%gamma_only )

!
! ... diagnostics
      IF( ionode ) THEN 
        WRITE( stdout,100)
        WRITE( stdout,150)
        WRITE( stdout,250) 1,(a1(i)/alat, i=1,3)
        WRITE( stdout,250) 2,(a2(i)/alat, i=1,3)
        WRITE( stdout,250) 3,(a3(i)/alat, i=1,3)
        WRITE( stdout,110)
        WRITE( stdout,200) 1,(b1(i), i=1,3)
        WRITE( stdout,200) 2,(b2(i), i=1,3)
        WRITE( stdout,200) 3,(b3(i), i=1,3)
      END IF

 100  FORMAT(/' * Reciprocal space initialization ( from ggen ) :',/)
 110  FORMAT(/'   Primitive vectors (units: 2*PI/ALAT):')
 150  FORMAT(/'   Real Space vectors (cart. coord., units: ALAT):')
 200  FORMAT(10X,'B',I1,' = ',3(2X,F6.2))
 250  FORMAT(10X,'A',I1,' = ',3(2X,F6.2))

      s3 = cclock()

! ... setup reciprocal space mesh
      b1 = b1 / alat
      b2 = b2 / alat
      b3 = b3 / alat
      CALL gindexset(gv, mill_g, b1, b2, b3)
      s4 = cclock()

      IF( tdipole ) THEN
        IF( ionode ) WRITE( stdout,800)
        CALL berry_setup( gv%ngw_l, gv%ngw_g, nr1, nr2, nr3, mill_g )
  800 FORMAT(   3X,'Polarizability using berry phase')
      END IF

      IF( ALLOCATED( g2_g ) ) deallocate(g2_g)
      IF( ALLOCATED( mill_g ) ) deallocate(mill_g)

      CALL fft_setup(gv, kp)
      s5 = cclock()

      RETURN
      END SUBROUTINE

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE compute_hg(hg_g, mill, i,j,k,b1,b2,b3,gcut,gkcut,ig,igw)

!  this routine computes squared lengths and indices for Fourier
!  transforms of G points
!  maximum values for nxh indices are also returned, for checking purpose
!  ----------------------------------------------

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(dbl), INTENT(OUT) :: hg_g(:)
      INTEGER, INTENT(OUT) :: mill(:,:)
      INTEGER, INTENT(IN)  :: i, j, k
      INTEGER, INTENT(INOUT) :: ig, igw
      REAL(dbl)  :: b1(3), b2(3), b3(3)
      REAL(dbl)  :: gcut, gkcut

! ... declare functions
      REAL(dbl)  miller2gsqr

! ... declare other variables
      REAL(dbl) :: gsq

! ... end of declarations
!  ----------------------------------------------

      gsq = ( REAL(i) * b1(1) + REAL(j) * b2(1) + REAL(k) * b3(1) )**2 + &
          & ( REAL(i) * b1(2) + REAL(j) * b2(2) + REAL(k) * b3(2) )**2 + &
          & ( REAL(i) * b1(3) + REAL(j) * b2(3) + REAL(k) * b3(3) )**2 
      IF(gsq .LE. gcut) THEN
        IF(gsq .LE. gkcut) igw = igw + 1
        ig = ig + 1
        IF(ig .GT. SIZE(hg_g)) THEN
          CALL errore(' compute_hg ',' too many G counts ',ig)
        END IF
        hg_g(ig) = gsq
        mill(1,ig) = i;  mill(2,ig) = j;  mill(3,ig) = k
      END IF

      RETURN
      END SUBROUTINE compute_hg


     END MODULE init_fpmd
