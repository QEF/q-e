!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE read_file_xspectra(xread_wf)
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only, io_level
  USE ions_base,            ONLY : nat, nsp, ityp, tau, if_pos, extfor
  USE basis,                ONLY : natomwfc
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk, npk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg, npwx
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym
  USE ktetra,               ONLY : tetra, ntetra 
  USE extfield,             ONLY : forcefield, tefield
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : gg, ngm, g, gcutm,&
                                   eigts1, eigts2, eigts3, nl, gstart
  USE gvecs,                ONLY : ngms, nls, gcutms
  USE grid_subroutines,     ONLY : realspace_grids_init
  USE recvec_subs,          ONLY : ggen
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions_module, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE buffers,              ONLY : open_buffer, close_buffer
  USE uspp_param,           ONLY : upf
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE mp_pools,             ONLY : kunit
  USE pw_restart,           ONLY : pw_readfile
  USE xml_io_base,          ONLY : pp_check_file
  USE uspp,                 ONLY : okvan, becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : lda_plus_u, eth, oatwfc
  USE dfunct,               ONLY : newd
  USE realus,               ONLY : qpointlist,betapointlist,init_realspace_vars,real_space
!<CG>
  USE paw_gipaw,            ONLY : set_paw_upf
!</CG>
  USE read_pseudo_mod,      ONLY : readpp       

  !
  IMPLICIT NONE
  !
  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, epaw, charge
  REAL(DP) :: sr(3,3,48)
  LOGICAL  :: exst
!<MCB>
  LOGICAL  :: xread_wf
!</MCB>
  !
  CHARACTER(len=80) :: input_dft = 'none'
  !
  ! ... first we get the version of the qexml file
  !     if not already read
  !
  CALL pw_readfile( 'header', ierr )
  CALL errore( 'read_file ', 'unable to determine qexml version', ABS(ierr) )
  !
  !
  ! ... first we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) &
     CALL infomsg( 'read_file', 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing' )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  ! ... a reset of the internal flags is necessary because some codes call
  ! ... read_file() more than once
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat < 0 ) &
     CALL errore( 'read_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  !
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  !
  ALLOCATE( irt( 48, nat ) )
  ALLOCATE( tetra( 4, MAX( ntetra, 1 ) ) )
  !
  ! ... here we read all the variables defining the system
  ! ... in parallel execution, only root proc read the file
  ! ... and then broadcast the values to all other procs
  !
  !-------------------------------------------------------------------------------
  ! ... XML punch-file
  !-------------------------------------------------------------------------------
  !
  CALL set_dimensions()
  CALL realspace_grids_init (dfftp,dffts,at, bg, gcutm, gcutms )  
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  if (cell_factor == 0.d0) cell_factor = 1.D0
!  lmovecell = .FALSE.
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  CALL pw_readfile( 'nowave', ierr )
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )
  !
  CALL readpp( input_dft )
  !
  !<CG>
  DO nt = 1, nsp
    call set_paw_upf(nt, upf(nt))
  ENDDO  
  !</CG>
  !
  okvan = ANY ( upf(:)%tvanp )
  okpaw = ANY ( upf(1:nsp)%tpawp )

  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... check for spin-orbit pseudopotentials
  !

!<CG>  this seems to be obsolete
!  DO nt = 1, nsp
!     !
!     so(nt) = upf(nt)%has_so  
!     !
!  END DO
!</CG>
  !

  !
  !<MCB>
  IF ( .NOT. xread_wf ) THEN
     CALL read_k_points()
     !<CG>
     nkstot=nks
     !</CG>
     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
     DEALLOCATE ( et, wg )
     ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  END IF
  !</MCB>


  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL allocate_fft()
  CALL ggen ( gamma_only, at, bg )
  CALL gshells ( lmovecell )
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF

  !<MCB> new in 4.0
  IF ( lda_plus_u ) THEN
     ALLOCATE ( oatwfc(nat) )
     CALL offset_atom_wfc ( nat, oatwfc )
  ENDIF
  !
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL pw_readfile( 'rho', ierr )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  CALL init_vloc()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  !
  CALL setlocal()
  !
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     !
     CALL fwfft ('Dense', psic, dfftp)
     !
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )

  !
  ! ... reads the wavefunctions and writes them in 'distributed' form 
  ! ... to unit iunwfc (for compatibility)
  !
  !<MCB>
  if(xread_wf) then

     nwordwfc = nbnd*npwx*npol
  !
     CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
     CALL pw_readfile( 'wave', ierr )
  !
     
     nks=nkstot
     if(lsda) then
        nks=nkstot/2
        CALL set_kup_and_kdw( xk, wk, isk, nks, npk )
     endif

     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )

  endif
!</MCB>

  CALL init_us_1()
  IF (okpaw) then
!     CALL compute_becsum(1)
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF
  if ( real_space ) THEN !initialisation of real space related stuff
    !OBM - correct parellism issues
    !call qpointlist()
    call betapointlist()
    call init_realspace_vars()
    !call betapointlist_v2()
    write(stdout,'(5X,"Real space initialisation completed")')
  endif
  
  CALL newd()

  !<MCB>
  if(xread_wf) CALL close_buffer  ( iunwfc, 'KEEP' )
  !</MCB>
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE wvfct,     ONLY : ecutwfc
      USE gvect,     ONLY : gcutm
      USE gvecs,   ONLY : gcutms, dual, doublegrid
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutm = dual * ecutwfc / tpiba2
      !
      doublegrid = ( dual > 4.D0 )
      !
      IF ( doublegrid ) THEN
         !
         gcutms = 4.D0 * ecutwfc / tpiba2
         !
      ELSE
         !
         gcutms = gcutm
         !
      END IF
      !
    END SUBROUTINE set_dimensions
    !

    !<MCB>


    subroutine read_k_points
      USE start_k,            ONLY : nk1, nk2, nk3, k1, k2, k3
      USE io_global,          ONLY : ionode_id 
      USE klist,              ONLY : npk
      USE constants,          ONLY : degspin
      USE parser,             ONLY : read_line
      USE mp,                 ONLY : mp_bcast
      USE mp_world,           ONLY : world_comm ! not sure about this
      implicit none
      INTEGER               :: npool, nkl, nkr, nkbl, iks, ike
      CHARACTER(LEN=256)         :: input_line
      INTEGER i,j,k,n
      !   Define new k-point mesh
      
      CALL read_line( input_line )
      READ(input_line, *) nk1, nk2, nk3, k1, k2 ,k3
      IF ( k1 < 0 .OR. k1 > 1 .OR. &
           k2 < 0 .OR. k2 > 1 .OR. &
           k3 < 0 .OR. k3 > 1 ) CALL errore &
           ('card_kpoints', 'invalid offsets: must be 0 or 1', 1)
      IF ( nk1 <= 0 .OR. nk2 <= 0 .OR. nk3 <= 0 ) CALL errore &
           ('card_kpoints', 'invalid values for nk1, nk2, nk3', 1)
      
      CALL mp_bcast( k1,  ionode_id, world_comm )
      CALL mp_bcast( k2,  ionode_id, world_comm )
      CALL mp_bcast( k3,  ionode_id, world_comm )
      CALL mp_bcast( nk1,  ionode_id,world_comm )
      CALL mp_bcast( nk2,  ionode_id,world_comm )
      CALL mp_bcast( nk3,  ionode_id,world_comm )
      
      !     if(nsym.eq.1) then
      nks=nk1*nk2*nk3
      do i=1,nk1
         do j=1,nk2
            do k=1,nk3
               !  this is nothing but consecutive ordering
               n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
               !  xkg are the components of the complete grid in crystal axis
               xk(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
               xk(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
               xk(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
            end do
         end do
      end do
      wk(1:nks)=1.0/DBLE(nks)
      call cryst_to_cart(nks,xk,bg,1)
      !     else
      !       CALL kpoint_grid( nsym, s, bg, npk, k1, k2, k3, &
      !          nk1, nk2, nk3, nks, xk, wk )
      !     endif
      
      if(lsda) then
         CALL set_kup_and_kdw( xk, wk, isk, nks, npk )
      ELSE IF ( noncolin ) THEN
         call errore('define_and_distribute_k_points', &
              'noncolinear not implemented', 1 )
      else 
         isk(1:nks)=1
         wk(1:nks)=2.d0/dfloat(nks)
      endif
      
      
    end subroutine read_k_points

    !</MCB>



  END SUBROUTINE read_file_xspectra
