!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

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
      USE reciprocal_vectors, ONLY : mill_g, g2_g, bi1, bi2, bi3
      USE recvecs_subroutines, ONLY: recvecs_init
      use gvecw, only: gcutw, gkcut
      use gvecp, only: ecut => ecutp, gcut => gcutp
      use gvecs, only: gcuts
      use gvecb, only: gcutb
      use fft_scalar, only: good_fft_dimension, good_fft_order
      USE fft_base, ONLY: dfftp, dffts, fft_dlay_descriptor
      USE stick_base,       ONLY: pstickset
      USE control_flags,    ONLY: tdipole
      USE berry_phase,      ONLY: berry_setup
      USE electrons_module, ONLY: bmeshset
      USE problem_size,     ONLY: cpsizes
      USE mp_global,        ONLY: nproc_image

      implicit none
! 
      integer  :: i
      real(8) :: rat1, rat2, rat3
      real(8) :: b1(3), b2(3), b3(3)
      integer :: ng_ , ngs_ , ngm_ , ngw_

      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( //, &
                3X,'Simulation dimensions initialization',/, &
                3X,'------------------------------------' )
      END IF
      !
      ! ... Initialize bands indexes for parallel linear algebra 
      ! ... (distribute bands to processors)
      !
      CALL bmeshset( )

      !
      ! ... Initialize (global) real and compute global reciprocal dimensions
      !

      CALL realspace_grids_init( alat, a1, a2, a3, gcut, gcuts, ng_ , ngs_ )

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

      ! ... printout g vector distribution summary
      !
      CALL gmeshinfo()
      !
      IF( program_name == 'FPMD' ) THEN
        !
        CALL cpsizes( nproc_image )
        !
      END IF
      !
      !   Flush stdout
      !
      CALL flush_unit( stdout )
      !
      return
      end subroutine init_dimensions




!-----------------------------------------------------------------------
      subroutine init_geometry ( )
!-----------------------------------------------------------------------
!
      USE kinds,            ONLY: DP
      use control_flags,    only: iprint, thdyn, ndr, nbeg, program_name, tbeg
      use io_global,        only: stdout, ionode
      use mp_global,        only: nproc_image
      USE io_files,         ONLY: scradir     
      use ions_base,        only: na, nsp, nat, tau_srt, ind_srt, if_pos, atm
      use cell_base,        only: a1, a2, a3, r_to_s, cell_init

      use cell_base,        only: ibrav, ainv, h, hold, tcell_base_init
      USE ions_positions,   ONLY: allocate_ions_positions, &
                                  atoms0, atomsm, atomsp
      use cp_restart,       only: cp_read_cell
      USE fft_base,         ONLY: dfftb, fft_dlay_descriptor
      USE fft_types,        ONLY: fft_box_allocate
      USE cp_main_variables,     ONLY: ht0, htm, taub
      USE brillouin,             ONLY: kp
      USE ions_module,           ONLY: atoms_init
      USE atoms_type_module,     ONLY: atoms_type

      implicit none
      !
      ! local
      !
      integer :: i, j
      real(DP) :: gvel(3,3), ht(3,3)
      real(DP) :: xnhh0(3,3), xnhhm(3,3), vnhh(3,3), velh(3,3)
      REAL(DP), ALLOCATABLE :: taus_srt( :, : )

      IF( .NOT. tcell_base_init ) &
         CALL errore( ' init_geometry ', ' cell_base_init has not been call yet! ', 1 )

      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( //, &
                3X,'System geometry initialization',/, &
                3X,'------------------------------' )
      END IF

      CALL cell_init( ht0, a1, a2, a3 )
      CALL cell_init( htm, a1, a2, a3 )

      CALL allocate_ions_positions( nsp, nat )
      ! 
      ! Scale positions that have been read from standard input 
      ! according to the cell given in the standard input too
      ! taus_srt = scaled, tau_srt = atomic units
      !
      ALLOCATE( taus_srt( 3, nat ) )
      
      CALL r_to_s( tau_srt, taus_srt, na, nsp, ainv )

      CALL atoms_init( atomsm, atoms0, atomsp, taus_srt, ind_srt, if_pos, atm, ht0%hmat )
      !
      DEALLOCATE( taus_srt )
      !
      !  Allocate box descriptor
      !
      ALLOCATE( taub( 3, nat ) )
      !
      CALL fft_box_allocate( dfftb, nproc_image, nat )
      !
      !  if tbeg = .true.  the geometry is given in the standard input even if
      !  we are restarting a previous run
      !
      if( ( nbeg > -1 ) .and. ( .not. tbeg ) ) then
        !
        ! read only h and hold from restart file "ndr"
        !
        CALL cp_read_cell( ndr, scradir, .TRUE., ht, hold, velh, gvel, xnhh0, xnhhm, vnhh )

        CALL cell_init( ht0, ht   )
        CALL cell_init( htm, hold )
        ht0%hvel = velh  !  set cell velocity
        ht0%gvel = gvel 

        h     = TRANSPOSE( ht   )
        ht    = TRANSPOSE( hold )
        hold  = ht
        ht    = TRANSPOSE( velh )
        velh  = ht

        WRITE( stdout,344) ibrav
        do i=1,3
          WRITE( stdout,345) (h(i,j),j=1,3)
        enddo
        WRITE( stdout,*)


      else
        !
        ! geometry is set to the cell parameters read from stdin ( a1, a2, a3 )
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
      call newinit( h )
!
      !
 344  format(3X,'ibrav = ',i4,'       cell parameters ',/)
 345  format(3(4x,f10.5))
      return
      end subroutine init_geometry



!-----------------------------------------------------------------------

    subroutine newinit( h )
      !
      !     re-initialization of lattice parameters and g-space vectors.
      !     Note that direct and reciprocal lattice primitive vectors
      !     a1,a2,a3, ainv, and corresponding quantities for small boxes
      !     are recalculated according to the value of cell parameter h
      !
      USE cell_base,             ONLY : a1, a2, a3, omega, alat, cell_base_reinit
      USE control_flags,         ONLY : program_name
      USE reciprocal_space_mesh, ONLY : newgk
      !
      implicit none
      !
      real(8) :: h(3,3)

      ! local
      !
      real(8) :: gmax, b1(3), b2(3), b3(3)
      !
      !  re-initialize the cell base module with the new geometry
      !
      CALL cell_base_reinit( TRANSPOSE( h ) )
      !
      call recips( a1, a2, a3, b1, b2, b3 )
      !
      call gcal( alat, b1, b2, b3, gmax )
      !
      !   generation of little box g-vectors
      !
      call newgb( a1, a2, a3, omega, alat )
      !
      IF( program_name == 'FPMD' ) THEN
        !
        CALL newgk( )
        !
      END IF
      !
      return
    end subroutine newinit
