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

       REAL(dbl), EXTERNAL :: cclock

       PUBLIC :: init1s, init0s

     CONTAINS

!  BEGIN manual

      SUBROUTINE init0s(gv, kp, ps, atoms_m, atoms_0, atoms_p, wfill, &
        wempt, ht_m, ht, fnl, eigr, nspin)

!  this routine handles data initialization
!  ----------------------------------------------
!  END manual

! ... declare modules
      use mp_global, only: nproc
      USE parameters, ONLY: nspinx
      USE phase_factors_module, ONLY: strucf
      USE cp_types
      USE atoms_type_module, ONLY: atoms_type
      USE time_step, ONLY: delt
      USE cell_module, ONLY: cell_init, get_lattice_vectors
      USE cell_base, ONLY: omega, alat
      USE electrons_module, ONLY: electron_mass_init, band_init, n_emp, bmeshset
      USE electrons_base, ONLY: nupdwn
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
      USE control_flags, ONLY: nbeg, tbeg, timing, t_diis, iprsta
      USE input_parameters, ONLY: rd_ht
      USE turbo, ONLY: tturbo, allocate_turbo
      USE ions_base, ONLY: tau_srt, tau_units, ind_srt, if_pos, atm
      USE stick, ONLY: dfftp
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE problem_size, ONLY: cpsizes
      USE reciprocal_space_mesh, ONLY : gindexset

      IMPLICIT NONE

! ... declare subroutine arguments

      TYPE (atoms_type)   :: atoms_0, atoms_p, atoms_m
      TYPE (wave_descriptor) :: wfill, wempt
      TYPE (pseudo) :: ps
      TYPE (projector) :: fnl(:,:)
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv
      TYPE (kpoints) :: kp
      TYPE (boxdimensions) :: ht_m, ht
      INTEGER :: nspin

! ... declare other variables
      REAL(dbl) :: s1, s2, s3, s4, s5
      REAL(dbl) :: a1(3), a2(3), a3(3)
      real(dbl) :: b1(3), b2(3), b3(3)
      INTEGER :: i, ispin, isym
      LOGICAL :: tk
      COMPLEX(dbl), ALLOCATABLE :: sfac(:,:)
      REAL(dbl) :: rat1, rat2, rat3
      INTEGER :: neupdwn( nspinx )


!  end of declarations
!  ----------------------------------------------

      s1 = cclock()

      call get_lattice_vectors( a1, a2, a3 )
      call recips( a1, a2, a3, b1, b2, b3 )

! ... arrange for reciprocal lattice vectors
      tk = .NOT. ( kp%scheme == 'gamma' )
      CALL allocate_recvecs(gv, ngm, ngmt, ngw, ngwt, tk, kp%nkpt)

      CALL gindexset( gv, b1, b2, b3 )

      ! ... set the bands mesh
      !
      CALL bmeshset( )

      CALL cpsizes( nproc )
      !
      CALL cpflush( )


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
      ELSE
        CALL cell_init( ht, a1, a2, a3 )
        CALL cell_init( ht_m, a1, a2, a3 )
      END IF

! ... initialize atomic configuration (should be called after metric_init)
      CALL atoms_init( atoms_m, atoms_0, atoms_p, tau_srt, ind_srt, if_pos, atm, tau_units, alat, ht )

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

      IF( iprsta > 2 ) THEN
        CALL wave_descriptor_info( wfill, 'wfill', stdout )
        CALL wave_descriptor_info( wempt, 'wempt', stdout )
      END IF

! ... if tturbo=.TRUE. some data is stored in memory instead of being
! ... recalculated (see card 'TURBO')
      IF( tturbo ) THEN
        CALL allocate_turbo( dfftp%nr1x, dfftp%nr2x, dfftp%npl )
      ENDIF


      CALL cpflush  ! flush output streams
      DEALLOCATE( sfac )

      RETURN
      END SUBROUTINE


!  BEGIN manual

    SUBROUTINE init1s(gv, kp, ps, atoms_m, atoms_0, atoms_p, cm, c0, wfill, &
      ce, wempt, ht_m, ht, fnl, eigr, occ)

!  this routine handles data initialization
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE phase_factors_module, ONLY: strucf
      USE wave_init, ONLY: pw_atomic_init
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
      USE control_flags, ONLY: nbeg, tbeg, timing, t_diis

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
      TYPE (boxdimensions) :: ht_m, ht

! ... declare other variables
      REAL(dbl) :: s1, s2, s3, s4, s5
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

      IF(t_diis) THEN
! ...   arrange for DIIS minimization
        CALL allocate_diis(gv%ngw_l, nbnd, kp%nkpt)
      END IF

      CALL cpflush  ! flush output streams

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------=!
   END MODULE init_fpmd
!=----------------------------------------------------------------------=!




!=----------------------------------------------------------------------=!
!
!   CP90 / FPMD common init subroutine 
!
!=----------------------------------------------------------------------=!


  subroutine init_dimensions(  )

      !
      !     initialize G-vectors and related quantities
      !

      use io_global, only: stdout, ionode
      use control_flags, only: program_name, gamma_only
      use grid_dimensions, only: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr => nnrx
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use small_box, only: a1b, a2b, a3b, omegab, ainvb, tpibab, small_box_set
      use small_box, only: alatb, b1b, b2b, b3b
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, nnrsx
      USE grid_subroutines, ONLY: realspace_grids_init, realspace_grids_para
      USE reciprocal_space_mesh, ONLY: gmeshinfo
      USE reciprocal_vectors, ONLY : mill_g, g2_g, bi1, bi2, bi3
      USE recvecs_subroutines, ONLY: recvecs_init
      use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix, gcutw, gkcut
      use gvecp, only: ecut => ecutp, gcut => gcutp
      use gvecs, only: gcuts
      use gvecb, only: gcutb
      use fft_scalar, only: good_fft_dimension, good_fft_order
      USE fft_base, ONLY: dfftp, dffts, fft_dlay_descriptor
      USE fft,       ONLY: fft_setup
      USE stick_base, ONLY: pstickset
      USE control_flags, ONLY: tdipole
      USE berry_phase, ONLY: berry_setup
      USE real_space_mesh, ONLY: realspace_procgrid_init
      USE bands_mesh, ONLY: bands_procgrid_init

      implicit none
! 
      integer  :: i
      real(kind=8) :: rat1, rat2, rat3
      real(kind=8) :: b1(3), b2(3), b3(3)
      integer :: ng_ , ngs_ , ngm_ , ngw_

      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( 3X,'Simulation dimensions initialization',/, &
                3X,'------------------------------------' )
      END IF
      !
      ! ... Initialize processor grid for parallel linear algebra 
      !     used electronic states lagrange multiplier matrixes
      !

      CALL bands_procgrid_init( )

      !
      ! ... Initialize (global) real and compute global reciprocal dimensions
      !

      CALL realspace_grids_init( alat, a1, a2, a3, gcut, gcuts, ng_ , ngs_ )

      !
      ! ... Initialize real space processor grid
      !

      CALL realspace_procgrid_init( )

      !
      ! ... cell dimensions and lattice vectors
      !

      call recips( a1, a2, a3, b1, b2, b3 )

      !     Store the base vectors used to generate the reciprocal space

      bi1 = b1
      bi2 = b2
      bi3 = b3

      !     Change units:  b1, b2, b3  are the 3 basis vectors generating 
      !     the reciprocal lattice in 2pi/alat units
      !
      !      Normally if a1, a2 and a3 are in cartesian coordinates
      !      and in a.u. units the corresponding bs are in cartesian
      !      coordinate too and in unit of 2 PI / a.u.
      !      now bring b1, b2 and b3 in units of 2 PI / alat

      b1 = b1 * alat
      b2 = b2 * alat
      b3 = b3 * alat

      IF( ionode ) THEN

        WRITE( stdout,210) 
210     format(/,3X,'unit vectors of full simulation cell',&
              &/,3X,'in real space:',25x,'in reciprocal space (units 2pi/alat):')
        WRITE( stdout,'(3X,I1,1X,3f10.4,10x,3f10.4)') 1,a1,b1
        WRITE( stdout,'(3X,I1,1X,3f10.4,10x,3f10.4)') 2,a2,b2
        WRITE( stdout,'(3X,I1,1X,3f10.4,10x,3f10.4)') 3,a3,b3

      END IF

      !
      do i=1,3
         ainv(1,i)=b1(i)/alat
         ainv(2,i)=b2(i)/alat
         ainv(3,i)=b3(i)/alat
      end do

      !
      ! ainv  is transformation matrix from cartesian to crystal coordinates
      !       if r=x1*a1+x2*a2+x3*a3 => x(i)=sum_j ainv(i,j)r(j)
      !       Note that ainv is really the inverse of a=(a1,a2,a3)
      !       (but only if the axis triplet is right-handed, otherwise
      !        for a left-handed triplet, ainv is minus the inverse of a)
      !

      ! ... set the sticks mesh and distribute g vectors among processors
      !

      CALL pstickset( dfftp, dffts, alat, a1, a2, a3, gcut, gkcut, gcuts, &
        nr1, nr2, nr3, nr1x, nr2x, nr3x, nr1s, nr2s, nr3s, nr1sx, nr2sx,   &
        nr3sx, ngw_ , ngm_ , ngs_ )
      !
      !
      ! ... Initialize reciprocal space local and global dimensions
      !     NOTE in a parallel run ngm_ , ngw_ , ngs_ here are the 
      !     local number of reciprocal vectors
      !
      CALL recvecs_init( ngm_ , ngw_ , ngs_ )
      !
      !
      ! ... Initialize (local) real space dimensions
      !
      CALL realspace_grids_para( dfftp, dffts )
      !
      !
      ! ... Initialize FFT module
      !
      CALL fft_setup( gamma_only , ngm_ , ngs_ , ngw_ )
      !
      ! ... printout g vector distribution summary
      !
      CALL gmeshinfo()
      !
      !
      ! ... generate g-space
      !
      call ggencp( b1, b2, b3, nr1, nr2, nr3, nr1s, nr2s, nr3s, gcut, gcuts, gkcut, gamma_only )

      ! 
      !  Allocate index required to compute polarizability
      !
      IF( tdipole ) THEN
        CALL berry_setup( ngw_ , mill_g )
      END IF

      !
      !     global arrays are no more needed
      !
      if( allocated( g2_g ) )   deallocate( g2_g )
      if( allocated( mill_g ) ) deallocate( mill_g )
      
      !
      !     generation of little box g-vectors
      !
      !  sets the small box parameters

      rat1 = DBLE( nr1b ) / DBLE( nr1 )
      rat2 = DBLE( nr2b ) / DBLE( nr2 )
      rat3 = DBLE( nr3b ) / DBLE( nr3 )
      CALL small_box_set( alat, omega, a1, a2, a3, rat1, rat2, rat3 )

      !  now set gcutb
      !
      gcutb = ecut / tpibab / tpibab
      !
      CALL ggenb ( b1b, b2b, b3b, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, gcutb )
      
      IF( ionode ) THEN
        WRITE( stdout, 110 )
 110    FORMAT( /,3X,'Simulation dimensions initialization completed',// )
      END IF
      !
      !   Flush stdout
      !
      CALL cpflush()
      !

      return
      end subroutine



!-----------------------------------------------------------------------
      subroutine init ( ibrav, celldm, ecut, ecutw, ndr, nbeg,  &
                        tfirst, tau0, taus, delt, tps, iforce )
!-----------------------------------------------------------------------
!
!     initialize G-vectors and related quantities
!     use ibrav=0 for generic cell vectors given by the matrix h(3,3)
!
      use control_flags, only: iprint, thdyn
      use io_global, only: stdout
      use gvecw, only: ngw
      use ions_base, only: na, pmass, nsp
      use cell_base, only: ainv, a1, a2, a3, r_to_s, s_to_r
      use constants, only: pi, fpi
      use cell_base, only: hold, h
      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix
      use betax, only: mmx, refg
      use cp_restart, only: cp_read_cell
      use parameters, only: nacx, nsx, natx, nhclm
      use electrons_base, only: f

      implicit none
! input/output
      integer ibrav, ndr, nbeg
      logical tfirst
      real(kind=8) tau0(3,natx), taus(3,natx)
      integer iforce(3,natx)
      real(kind=8) celldm(6), ecut, ecutw
      real(kind=8) delt, tps
! local
      integer i, j, ia, is, nfi, isa, isat
      real(kind=8) gvel(3,3)
      real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
!
!
! taus = scaled, tau0 = alat units
!
      CALL r_to_s( tau0, taus, na, nsp, ainv )
!
      refg = 1.0d0 * ecut / ( mmx - 1 )
      WRITE( stdout,*) '   NOTA BENE: refg, mmx = ',refg,mmx
!
      if( nbeg >= 0 ) then

        !
        ! read only h and hold from file ndr
        !

        CALL cp_read_cell( ndr, ' ', .TRUE., h, hold, velh, gvel, xnhh0, xnhhm, vnhh )

        h     = TRANSPOSE( h    )
        hold  = TRANSPOSE( hold )
        velh  = TRANSPOSE( velh )

        WRITE( stdout,344) ibrav
        do i=1,3
          WRITE( stdout,345) (h(i,j),j=1,3)
        enddo
        WRITE( stdout,*)

      else

        !
        ! with variable-cell we use h to describe the cell
        !

        do i = 1, 3
            h(i,1) = a1(i)
            h(i,2) = a2(i)
            h(i,3) = a3(i)
        enddo

        hold = h

      end if
!
!     ==============================================================
!     ==== generate true g-space                                ====
!     ==============================================================
!
      call newinit( ibrav )
!
      !
 344  format(' ibrav = ',i4,'       cell parameters ',/)
 345  format(3(4x,f10.5))
      return
      end



!-----------------------------------------------------------------------
      subroutine newinit(ibrav)
!-----------------------------------------------------------------------
!     re-initialization of lattice parameters and g-space vectors.
!     Note that direct and reciprocal lattice primitive vectors
!     a1,a2,a3, ainv, and corresponding quantities for small boxes
!     are recalculated according to the value of cell parameter h
!
      use control_flags, only: iprint, iprsta
      use io_global, only: stdout
      use gvec
      use grid_dimensions, only: nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use constants, only: pi, fpi
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
      use small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
      use cell_base, only: h, deth
      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix
!
      implicit none
      integer ibrav
!
! local
      integer i, j
      real(kind=8) alatb, gmax, b1(3),b2(3),b3(3), b1b(3),b2b(3),b3b(3)
      real(kind=8) ddum
!
!
      alat = sqrt( h(1,1)*h(1,1) + h(2,1)*h(2,1) + h(3,1)*h(3,1) )

!     ==============================================================
      tpiba  = 2.d0 * pi / alat
      tpiba2 = tpiba * tpiba

!     ==============================================================
!     ==== generate g-space                                     ====
!     ==============================================================
      call invmat (3, h, ainv, deth)
      omega = deth
!
      do i = 1, 3
         a1(i) = h(i,1)
         a2(i) = h(i,2)
         a3(i) = h(i,3)
      enddo
!
      call recips( a1, a2, a3, b1, b2, b3 )
      b1 = b1 * alat
      b2 = b2 * alat
      b3 = b3 * alat
      call gcal( b1, b2, b3, gmax )
!
!     ==============================================================
!     generation of little box g-vectors
!     ==============================================================
!
      call newgb( a1, a2, a3, omega, alat )

!     ==============================================================
      if(iprsta.ge.4)then
         WRITE( stdout,34) ibrav,alat,omega
         if(ibrav.eq.0) then
            WRITE( stdout,344)
            do i=1,3
               WRITE( stdout,345) (h(i,j),j=1,3)
            enddo
            WRITE( stdout,*)
         endif
      endif
!
 34   format(' initialization ',//,                                     &
     &       ' ibrav=',i3,' alat=',f7.3,' omega=',f10.4,//)
 344  format(' cell parameters ',/)
 345  format(3(4x,f10.5))
!
      return
      end

