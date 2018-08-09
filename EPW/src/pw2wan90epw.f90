  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PP/pw2wannier - Quantum-ESPRESSO group               
  !------------------------------------------------------------------------
  SUBROUTINE pw2wan90epw 
  !------------------------------------------------------------------------
  !! This is the interface to the Wannier90 code: see http://www.wannier.org
  !!
  !!
  !! 10/2008  Parellel computation of Amn and Mmn 
  !! 12/2008  Added phase setting of overlap matrix elements
  !! 02/2009  works with standard nk1*nk2*nk3 grids
  !! 12/2009  works with USPP 
  !! 12/2014  RM: Imported the noncolinear case implemented by xlzhang
  !! 06/2016  SP: Debug of SOC + print/reading of nnkp file
  !!
  !------------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : prefix
  USE epwcom,     ONLY : write_wfn
  USE noncollin_module, ONLY : noncolin
  USE wannierEPW, ONLY : seedname2, wvfn_formatted, reduce_unk, ispinw, &
                         ikstart, ikstop, iknum
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=4)   :: spin_component
  !! Determine the spin case
  CHARACTER(len=256) :: outdir
  !! Name of the output directory
  !
  !
  outdir = './'
  seedname2 = prefix
  spin_component = 'none'
  wvfn_formatted = .false.
  reduce_unk= .false.
  !
  !
  WRITE(stdout,*)
  SELECT CASE ( TRIM( spin_component ) )
  CASE ( 'up' )
     WRITE(stdout,*) '    Spin CASE ( up )'
     ispinw  = 1
     ikstart = 1
     ikstop  = nkstot/2
     iknum   = nkstot/2
  CASE ( 'down' )
     WRITE(stdout,*) '    Spin CASE ( down )'
     ispinw = 2
     ikstart = nkstot/2 + 1
     ikstop  = nkstot
     iknum   = nkstot/2
  CASE DEFAULT
     IF (noncolin) THEN
        WRITE(stdout,*) '    Spin CASE ( non-collinear )'
     ELSE
        WRITE(stdout,*) '    Spin CASE ( default = unpolarized )'
     ENDIF
     ispinw = 0
     ikstart = 1
     ikstop  = nkstot
     iknum   = nkstot
  END SELECT
  !
  WRITE(stdout,*) 
  WRITE(stdout,*) '    Initializing Wannier90'
  WRITE(stdout,*) 
  CALL setup_nnkp
  CALL ylm_expansion
  CALL compute_amn_para
  CALL compute_mmn_para
  !
  CALL phases_a_m
  !
  CALL write_band
  !
  IF (write_wfn) CALL write_plot
  !
  WRITE(stdout,*)
  WRITE(stdout,*) '    Running Wannier90'
  CALL run_wannier
  !
  CALL lib_dealloc
  !
  END SUBROUTINE pw2wan90epw
!
!-------------------------------------------------------------------------
SUBROUTINE lib_dealloc
  !-----------------------------------------------------------------------
  !! Routine to de-allocate Wannier related matrices. 
  !
  USE wannierEPW
  !
  IMPLICIT NONE
  ! 
  IF (ALLOCATED(m_mat) )     DEALLOCATE(m_mat)
  IF (ALLOCATED(u_mat) )     DEALLOCATE(u_mat)
  IF (ALLOCATED(u_mat_opt) ) DEALLOCATE(u_mat_opt)
  IF (ALLOCATED(a_mat) )     DEALLOCATE(a_mat)
  IF (ALLOCATED(eigval) )    DEALLOCATE(eigval)
  IF (ALLOCATED(lwindow) )   DEALLOCATE(lwindow)
  !
END SUBROUTINE lib_dealloc
!
!-------------------------------------------------------------------------
SUBROUTINE setup_nnkp (  )
  !-----------------------------------------------------------------------
  !! 
  !! This routine write and read the .nnkp file. 
  !! The file specifies
  !! 1) The initial projections functions in the format
  !! num_proj
  !! proj_site(1,i),proj_site(2,i),proj_site(3,i) proj_l(i),proj_m(i),proj_radial(i)
  !! proj_z(1,i),proj_z(2,i),proj_z(3,i),proj_x(1,i),proj_x(2,i),proj_x(3,i), proj_zona(i)
  !! proj_s(i), proj_s_qaxis(1,i),proj_s_qaxis(2,i),proj_s_qaxis(3,i) 
  !!
  !! 2) begin nnkpts: the nearest neighbours of each         
  !! k-point, and therefore provides the information required to      
  !! calculate the M_mn(k,b) matrix elements 
  !! 
  ! ---------------------------------------------------------------------- 
  USE io_global, ONLY : meta_ionode, stdout, meta_ionode_id
  USE mp_world,  ONLY : world_comm
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps6, tpi
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE mp,        ONLY : mp_bcast
  USE wvfct,     ONLY : nbnd, npwx
  USE wannierEPW, ONLY : num_nnmax, mp_grid, atcart, atsym, kpb, g_kpb, &
                         center_w, alpha_w, l_w, mr_w, r_w, zaxis,      &
                         xaxis, excluded_band, rlatt, glatt, gf,        &
                         csph, ig_, iknum, seedname2, kpt_latt, nnb,    &
                         num_bands, n_wannier, nexband, nnbx, n_proj,   &
                         spin_eig, spin_qaxis
  USE noncollin_module, ONLY : noncolin
  USE constants_epw,    ONLY : bohr
  USE mp_global,        ONLY : intra_pool_comm, mp_sum
  USE epwcom,           ONLY : nbndskip
!  USE w90_PARAMETERs,   ONLY : postproc_setup
  USE w90_io,           ONLY : post_proc_flag
  ! 
  IMPLICIT NONE
  REAL(DP) :: g_(3)
  !! Temporary vector G_k+b, g_(:) = g_kpb(:,ik,ib) 
  REAL(DP) :: gg_
  !! Square of g_(3)
  INTEGER :: ik
  !! Counter on k-points
  INTEGER :: ib
  !! Counter on b-vectors
  INTEGER :: ig
  !! Index on G_k+b vectors
  INTEGER :: iw
  !! Counter on number of projections 
  INTEGER :: ia
  !! Counter on number of atoms
  INTEGER  :: type
  !! Type of atom
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER  :: indexb
  !! Index of exclude_bands 
  INTEGER  :: ipol
  !! Counter on polarizations
  INTEGER  :: idum
  !! Dummy index for reading nnkp file
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  !! Temporary index on G_k+b vectors
  REAL(DP) :: xnorm
  !! Norm of xaxis
  REAL(DP) :: znorm
  !! Norm of zaxis
  REAL(DP) :: coseno
  !! Cosine between xaxis and zaxis
  INTEGER  :: exclude_bands(nbnd)
  !! Bands excluded from the calcultion of WFs
  INTEGER  :: i
  !! Counter on polarizations
  INTEGER  :: iun_nnkp
  !! Unit for .nnkp file
  LOGICAL  :: have_nnkp
  !! Check if the .nnkp file exists.
  LOGICAL  :: found
  !! Check if the section in the .nnkp file is found. 
  INTEGER, EXTERNAL :: find_free_unit
  !! Look for a free unit for the .nnkp file.

! SP: An interface is required because the Wannier routine has optional
!     arguments
  Interface
    Subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc,&
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)

     USE kinds,       ONLY : dp
     USE wannierEPW,  ONLY : num_nnmax

     IMPLICIT NONE
     !
     CHARACTER(len=*), INTENT(in) :: seed__name
     INTEGER, DIMENSION(3), INTENT(in) :: mp_grid_loc
     INTEGER, INTENT(in) :: num_kpts_loc
     REAL(kind=DP), DIMENSION(3,3), INTENT(in) :: REAL_lattice_loc
     REAL(kind=DP), DIMENSION(3,3), INTENT(in) :: recip_lattice_loc
     REAL(kind=DP), DIMENSION(3,num_kpts_loc), INTENT(in) :: kpt_latt_loc
     INTEGER, INTENT(in) :: num_bands_tot
     INTEGER, INTENT(in) :: num_atoms_loc
     CHARACTER(len=*), DIMENSION(num_atoms_loc), INTENT(in) :: atom_symbols_loc
     REAL(kind=DP), DIMENSION(3,num_atoms_loc), INTENT(in) :: atoms_cart_loc
     LOGICAL, INTENT(in) :: gamma_only_loc
     LOGICAL, INTENT(in) :: spinors_loc
     INTEGER, INTENT(out) :: nntot_loc
     INTEGER, DIMENSION(num_kpts_loc,num_nnmax), INTENT(out) :: nnlist_loc
     INTEGER, DIMENSION(3,num_kpts_loc,num_nnmax), INTENT(out) :: nncell_loc
     INTEGER, INTENT(out) :: num_bands_loc
     INTEGER, INTENT(out) :: num_wann_loc
     REAL(kind=DP), DIMENSION(3,num_bands_tot), INTENT(out) :: proj_site_loc
     INTEGER, DIMENSION(num_bands_tot), INTENT(out) :: proj_l_loc
     INTEGER, DIMENSION(num_bands_tot), INTENT(out) :: proj_m_loc
     INTEGER, DIMENSION(num_bands_tot), INTENT(out) :: proj_radial_loc
     REAL(kind=DP), DIMENSION(3,num_bands_tot), INTENT(out) :: proj_z_loc
     REAL(kind=DP), DIMENSION(3,num_bands_tot), INTENT(out) :: proj_x_loc
     REAL(kind=DP), DIMENSION(num_bands_tot), INTENT(out) :: proj_zona_loc
     INTEGER, DIMENSION(num_bands_tot), INTENT(out) :: exclude_bands_loc
     INTEGER, DIMENSION(num_bands_tot), OPTIONAL, INTENT(out) :: proj_s_loc
     REAL(kind=DP), DIMENSION(3,num_bands_tot), OPTIONAL, INTENT(out) :: proj_s_qaxis_loc

    End Subroutine wannier_setup
  End Interface

  num_nnmax = 32

  ! aam: translations between PW2Wannier90 and Wannier90
  ! pw2wannier90   <==>   Wannier90
  !    nbnd                num_bands_tot
  !    n_wannier           num_wann
  !    num_bands           num_bands
  !    nat                 num_atoms
  !    iknum               num_kpts
  !    rlatt               transpose(real_lattice)
  !    glatt               transpose(recip_lattice)
  !    kpt_latt            kpt_latt
  !    nnb                 nntot
  !    kpb                 nnlist
  !    g_kpb               nncell
  !    mp_grid             mp_grid
  !    center_w            proj_site
  !    l_w,mr_w,r_w        proj_l,proj_m,proj_radial
  !    xaxis,zaxis         proj_x,proj_z
  !    alpha_w             proj_zona
  !    exclude_bands       exclude_bands
  !    atcart              atoms_cart
  !    atsym               atom_symbols

  ALLOCATE( atcart(3,nat), atsym(nat) )
  ALLOCATE( kpb(iknum,num_nnmax), g_kpb(3,iknum,num_nnmax) )
  ALLOCATE( center_w(3,nbnd), alpha_w(nbnd), l_w(nbnd), &
       mr_w(nbnd), r_w(nbnd), zaxis(3,nbnd), xaxis(3,nbnd) )
  ALLOCATE( excluded_band(nbnd) )

  ! real lattice (Cartesians, Angstrom)
  rlatt(:,:) = transpose(at(:,:))*alat*bohr
  ! reciprocal lattice (Cartesians, Angstrom)
  glatt(:,:) = transpose(bg(:,:))*tpi/(alat*bohr)
  ! atom co-ordinates in Cartesian co-ords and Angstrom units
  atcart(:,:) = tau(:,:)*bohr*alat
  ! atom symbols
  DO ia = 1, nat
     type = ityp(ia)
     atsym(ia) = atm(type)
  ENDDO

  IF (meta_ionode) THEN
    !postproc_setup = .true.
    post_proc_flag = .true.
    !print*,'postproc_setup ',postproc_setup
    CALL wannier_setup(seedname2, mp_grid, iknum,        &  ! input
           rlatt, glatt, kpt_latt, nbnd,                 &  ! input
           nat, atsym, atcart, .false., noncolin,        &  ! input
           nnb, kpb, g_kpb, num_bands, n_wannier,        &  ! output
           center_w, l_w, mr_w, r_w, zaxis,              &  ! output
           xaxis, alpha_w, exclude_bands)                   ! output
   ! SP: In wannier_setup, the .nnkp file is produced.
  ENDIF
   
  CALL mp_bcast(nnb,          meta_ionode_id, world_comm )
  CALL mp_bcast(kpb,          meta_ionode_id, world_comm )
  CALL mp_bcast(g_kpb,        meta_ionode_id, world_comm )
  CALL mp_bcast(num_bands,    meta_ionode_id, world_comm )
  CALL mp_bcast(n_wannier,    meta_ionode_id, world_comm )
  CALL mp_bcast(center_w,     meta_ionode_id, world_comm )
  CALL mp_bcast(l_w,          meta_ionode_id, world_comm )
  CALL mp_bcast(mr_w,         meta_ionode_id, world_comm )
  CALL mp_bcast(r_w,          meta_ionode_id, world_comm )
  CALL mp_bcast(zaxis,        meta_ionode_id, world_comm )
  CALL mp_bcast(xaxis,        meta_ionode_id, world_comm )
  CALL mp_bcast(alpha_w,      meta_ionode_id, world_comm )
  CALL mp_bcast(exclude_bands,meta_ionode_id, world_comm )
  CALL mp_bcast(noncolin,     meta_ionode_id, world_comm )
  !
  ! SP: Commented because we now write on file the .nnkp file and read from it.
  ! 
  ! n_proj = nr. of projections (=#WF unless spinors then =#WF/2) 
  !IF (noncolin) THEN
  !   n_proj = n_wannier/2
  !ELSE
     n_proj = n_wannier
  !ENDIF
  !
  WRITE (stdout,*)
  WRITE (stdout,*) '    Initial Wannier projections'
  WRITE (stdout,*)
  !
  DO iw = 1, n_proj
     WRITE (stdout, '(5x,"(",3f10.5,") :  l = ",i3, " mr = ", i3)') center_w(:,iw), l_w(iw), mr_w(iw)
  ENDDO
  !
  excluded_band(1:nbnd) = .false.
  nexband = 0
  band_loop: DO ibnd = 1, nbnd
     indexb = exclude_bands(ibnd)
     IF (indexb>nbnd .or. indexb<0) THEN
        CALL errore('setup_nnkp',' wrong excluded band index ', 1)
     ELSEIF (indexb.eq.0) THEN
        EXIT band_loop
     ELSE
        nexband = nexband + 1
        excluded_band(indexb) = .true.
     ENDIF
  ENDDO band_loop

  WRITE(stdout,'(/,"      - Number of bands is (",i3,")")') num_bands
  WRITE(stdout,'("      - Number of total bands is (",i3,")")') nbnd
  WRITE(stdout,'("      - Number of excluded bands is (",i3,")")') nexband
  WRITE(stdout,'("      - Number of wannier functions is (",i3,")")') n_wannier
  IF ((nexband .gt. 0) .and. (nbndskip .ne. nexband)) THEN
    WRITE(stdout,'(/5x,"Warning: check if nbndskip = ", i3 " makes sense since ", i3, &
                  " bands are excluded from wannier projection")') nbndskip, nexband
  ENDIF
  !
  IF ( (nbnd-nexband).ne.num_bands ) &
      CALL errore('setup_nnkp',' something wrong with num_bands',1)
  ! 
  ! Now we read the .nnkp file 
  !
  IF (meta_ionode) THEN  ! Read nnkp file on ionode only
    INQUIRE(file=trim(seedname2)//".nnkp",exist=have_nnkp)
    IF(.not. have_nnkp) THEN
       CALL errore( 'pw2wannier90', 'Could not find the file '&
          &//trim(seedname2)//'.nnkp', 1 )
    ENDIF
    iun_nnkp = find_free_unit()
    OPEN(unit=iun_nnkp, file=trim(seedname2)//".nnkp",form='formatted')
  ENDIF
  !
  IF (meta_ionode) THEN   ! read from ionode only
    IF (noncolin) THEN
       CALL scan_file_to(iun_nnkp,'spinor_projections',found)
       IF(.not. found) THEN
         !try old style projections
         CALL scan_file_to(iun_nnkp,'projections',found)
         IF(.not. found) THEN
           CALL errore( 'pw2wannier90', 'Could not find projections block in '&
              &//trim(seedname2)//'.nnkp', 1 )
         ENDIF
       ENDIF
    ELSE
      CALL scan_file_to(iun_nnkp,'projections',found)
      IF (.not. found) THEN
        CALL errore( 'pw2wannier90', 'Could not find projections block in '&
           &//trim(seedname2)//'.nnkp', 1 )
      ENDIF
    ENDIF
    READ(iun_nnkp,*) n_proj
  ENDIF
  CALL mp_bcast(n_proj, meta_ionode_id, world_comm)
  ! 
  ALLOCATE( gf(npwx,n_proj), csph(16,n_proj) )
  IF(noncolin) ALLOCATE( spin_eig(n_proj), spin_qaxis(3,n_proj) )
  ! 
  IF (meta_ionode) THEN   ! read from ionode only
    DO iw = 1, n_proj
      READ(iun_nnkp,*) (center_w(i,iw), i=1,3), l_w(iw), mr_w(iw), r_w(iw)
      READ(iun_nnkp,*) (zaxis(i,iw),i=1,3),(xaxis(i,iw),i=1,3),alpha_w(iw)
      xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
           xaxis(3,iw)*xaxis(3,iw))
      IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
      znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
           zaxis(3,iw)*zaxis(3,iw))
      IF (znorm < eps6) CALL errore ('read_nnkp',' |zaxis| < eps ',1)
      coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
           xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
      IF (abs(coseno) > eps6) &
           CALL errore('read_nnkp',' xaxis and zaxis are not orthogonal !',1)
      IF (alpha_w(iw) < eps6) &
           CALL errore('read_nnkp',' zona value must be positive', 1)
      ! convert wannier center in cartesian coordinates (in unit of alat)
      CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
      IF (noncolin) THEN
         READ(iun_nnkp,*) spin_eig(iw),(spin_qaxis(i,iw),i=1,3)
         xnorm = sqrt(spin_qaxis(1,iw)*spin_qaxis(1,iw) + spin_qaxis(2,iw)*spin_qaxis(2,iw) + &
           spin_qaxis(3,iw)*spin_qaxis(3,iw))
         IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
         spin_qaxis(:,iw)=spin_qaxis(:,iw)/xnorm
      ENDIF
    ENDDO
  ENDIF
  ! 
  WRITE(stdout,*) '     - All guiding functions are given '
  ! 
  ! Broadcast
  CALL mp_bcast(center_w,meta_ionode_id, world_comm)
  CALL mp_bcast(l_w,meta_ionode_id, world_comm)
  CALL mp_bcast(mr_w,meta_ionode_id, world_comm)
  CALL mp_bcast(r_w,meta_ionode_id, world_comm)
  CALL mp_bcast(zaxis,meta_ionode_id, world_comm)
  CALL mp_bcast(xaxis,meta_ionode_id, world_comm)
  CALL mp_bcast(alpha_w,meta_ionode_id, world_comm)
  IF(noncolin) THEN
     CALL mp_bcast(spin_eig,meta_ionode_id, world_comm)
     CALL mp_bcast(spin_qaxis,meta_ionode_id, world_comm)
  ENDIF
  !
  IF (meta_ionode) THEN   ! read from ionode only
    CALL scan_file_to(iun_nnkp,'nnkpts',found)
    IF(.not.found) THEN
       CALL errore( 'pw2wannier90epw', 'Could not find nnkpts block in '&
          &//trim(seedname2)//'.nnkp', 1 )
    ENDIF
    READ (iun_nnkp,*) nnb
  ENDIF
  ! Broadcast
  CALL mp_bcast(nnb,meta_ionode_id, world_comm)
  !
  !
  nnbx = 0
  nnbx = max (nnbx, nnb )
  ALLOCATE ( ig_(iknum,nnbx), ig_check(iknum,nnbx) )
  !
  ! Read data about neighbours
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading data about k-point neighbours '
  WRITE(stdout,*)
  IF (meta_ionode) THEN
    DO ik = 1, iknum
      DO ib = 1, nnb
        READ(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
      ENDDO
    ENDDO
  ENDIF
  ! Broadcast
  CALL mp_bcast(kpb, meta_ionode_id, world_comm)
  CALL mp_bcast(g_kpb, meta_ionode_id, world_comm)
  ! 
  DO ik =1, iknum
    DO ib = 1, nnb
      g_(:) = REAL( g_kpb(:,ik,ib) )
      CALL cryst_to_cart (1, g_, bg, 1)
      gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
      ig_(ik,ib) = 0
      ig = 1
      DO WHILE  (gg(ig) <= gg_ + eps6)
        IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
             (abs(g(2,ig)-g_(2)) < eps6) .and.  &
             (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
        ig = ig + 1
      ENDDO
    ENDDO
  ENDDO
  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik = 1, iknum
    DO ib = 1, nnb
      IF (ig_check(ik,ib) ==0) &
        CALL errore('read_nnkp', &
                    ' g_kpb vector is not in the list of Gs', 100*ik+ib )
    ENDDO
  ENDDO
  DEALLOCATE (ig_check)
  !
  WRITE(stdout,*) '     - All neighbours are found '
  WRITE(stdout,*)
  !
  IF (meta_ionode) THEN
    CLOSE (iun_nnkp)
  ENDIF
  RETURN
  !
END SUBROUTINE setup_nnkp
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
SUBROUTINE scan_file_to (iun_nnkp,keyword,found)
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: keyword
  !! Keyword searched for
  INTEGER, INTENT(in)          :: iun_nnkp
  !! Unit for .nnkp file
  LOGICAL, INTENT(out)         :: found
  !! Check if the section in the .nnkp file is found.
  CHARACTER(len=80)            :: line1, line2

! by uncommenting the following line the file scan restarts every time
! from the beginning thus making the reading independent on the order
! of data-blocks
!  rewind (iun_nnkp)
!
10 CONTINUE
  READ(iun_nnkp,*,end=20) line1, line2
  IF(line1/='begin')  GOTO 10
  IF(line2/=keyword) GOTO 10
  found=.true.
  RETURN
20 found=.false.
  rewind (iun_nnkp)

END SUBROUTINE scan_file_to
! ------------------------------------------------------------------------
!
! ------------------------------------------------------------------------
SUBROUTINE run_wannier
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, meta_ionode_id, meta_ionode
  USE ions_base, ONLY : nat
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE cell_base, ONLY : celldm
  USE io_files,  ONLY : prefix
  USE io_epw,    ONLY : QPeig_read
  USE pwcom,     ONLY : nkstot
  USE wannierEPW,ONLY : u_mat, lwindow, wann_centers, wann_spreads, eigval,  &
                        n_wannier, spreads, nnb, rlatt, glatt, kpt_latt,     &
                        iknum, seedname2, num_bands, u_mat_opt, atsym, a_mat,&
                        atcart, m_mat, mp_grid
  USE epwcom,    ONLY : eig_read
  USE wvfct,     ONLY : nbnd
  USE constants_epw, ONLY : czero, bohr
  !
  IMPLICIT NONE
  ! 
  INTEGER :: iw 
  !! Counter on wannier functions
  INTEGER :: ik
  !! Counter of k-point index
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: ios
  !! Integer variable for I/O control
  CHARACTER (len=256) :: tempfile
  !! Temporary file
  CHARACTER (len=80) :: line
  !! Temporary character
  !
  ALLOCATE(u_mat(n_wannier,n_wannier,iknum))
  ALLOCATE(u_mat_opt(num_bands,n_wannier,iknum))
  ALLOCATE(lwindow(num_bands,iknum))
  ALLOCATE(wann_centers(3,n_wannier))
  ALLOCATE(wann_spreads(n_wannier))
  !
  u_mat = czero
  u_mat_opt = czero
  !
  IF (meta_ionode) THEN
    ! read in external eigenvalues, e.g.  GW
    IF (eig_read) then
      WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
           nbnd, ",", nkstot,")"
      tempfile=trim(prefix)//'.eig'
      OPEN(QPeig_read, file=tempfile, form='formatted', action='read', iostat=ios)
      IF (ios /= 0) CALL errore ('run_wannier','error opening' // tempfile, 1)
      READ (QPeig_read,'(a)') line
      DO ik = 1, nkstot
        ! We do not save the k-point for the moment ==> should be read and
        ! tested against the current one  
        READ (QPeig_read,'(a)') line
        READ (QPeig_read,*) eigval (:,ik)
      ENDDO
      CLOSE(QPeig_read)
    ENDIF

! SP : This file is not used for now. Only required to build the UNK file
!      tempfile=trim(prefix)//'.mmn'
!      OPEN(iummn, file=tempfile, iostat=ios, form='unformatted')
!      WRITE(iummn) m_mat
!      CLOSE(iummn)

    CALL wannier_run(seedname2, mp_grid, iknum,   &                 ! input
         rlatt, glatt, kpt_latt, num_bands,       &                 ! input
         n_wannier, nnb, nat, atsym,              &                 ! input
         atcart, .false., m_mat, a_mat, eigval,   &                 ! input
         u_mat, u_mat_opt, lwindow, wann_centers, &                 ! output
         wann_spreads, spreads)                                     ! output
  ENDIF
  !
  CALL mp_bcast(u_mat,       meta_ionode_id, world_comm )
  CALL mp_bcast(u_mat_opt,   meta_ionode_id, world_comm )
  CALL mp_bcast(lwindow,     meta_ionode_id, world_comm )
  CALL mp_bcast(wann_centers,meta_ionode_id, world_comm )
  CALL mp_bcast(wann_spreads,meta_ionode_id, world_comm )
  CALL mp_bcast(spreads,     meta_ionode_id, world_comm )
  !
  !
  ! output the results of the wannierization
  !
  WRITE (stdout,*)
  WRITE (stdout,*) '    Wannier Function centers (cartesian, alat) and spreads (ang):'
  WRITE (stdout,*)
  ! RM - loop is up to n_wannier according to W90/wannierise.F90
  DO iw = 1, n_wannier
     WRITE (stdout, '(5x,"(",3f10.5,") :  ",f8.5)') wann_centers(:,iw)/celldm(1)/bohr, wann_spreads(iw)
  ENDDO
  WRITE (stdout,*)
  !
  ! store the final minimisation matrix on disk for later use
  !
  CALL write_filukk
  !
  RETURN
  !
END SUBROUTINE run_wannier
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE compute_amn_para
!-----------------------------------------------------------------------
!!  adapted from compute_amn in pw2wannier90.f90
!!  parallelization on k-points has been added
!!  10/2008 Jesse Noffsinger UC Berkeley
!!  01/2018 Roxana Margine updated 
!!
  USE io_global,       ONLY : stdout 
  USE kinds,           ONLY : DP
  USE klist,           ONLY : xk, nks, igk_k
  USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
  USE gvecw,           ONLY : ecutwfc
  USE wavefunctions,  ONLY : evc
  USE gvect,           ONLY : g, ngm
  USE cell_base,       ONLY : tpiba2
  USE uspp,            ONLY : nkb, vkb
  USE becmod,          ONLY : becp, calbec, deallocate_bec_type, allocate_bec_type
  USE wannierEPW,      ONLY : csph, excluded_band, gf, num_bands, &
                              n_wannier, iknum, n_proj, a_mat,  spin_qaxis, &
                              spin_eig, gf_spinor, sgf_spinor
  USE uspp_param,      ONLY : upf
  USE noncollin_module,ONLY : noncolin, npol
  USE constants_epw,   ONLY : czero, zero, eps6
#if defined(__NAG)
  USE f90_unix_io,     ONLY : flush
#endif
  USE mp_global,       ONLY : my_pool_id, npool, intra_pool_comm, inter_pool_comm
  USE mp,              ONLY : mp_sum
  ! 
  IMPLICIT NONE
  ! 
  INTEGER :: iw
  !! Counter on number of projections
  INTEGER :: ik
  !! Counter of k-point index
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: ibnd1
  !! Band index
  INTEGER :: ipool
  !! Index of current pool
  INTEGER ::  nkq
  !! Index of k-point in the pool 
  INTEGER ::  nkq_abs
  !! Absolute index of k-point 
  INTEGER ::  ik_g
  !! Temporary index of k-point, ik_g = nkq_abs
  INTEGER  :: ipol
  !! Index of spin-up/spin-down polarizations
  INTEGER ::  istart
  !! Index on plane waves
  !
  REAL(kind=DP) :: zero_vect(3)
  !! Temporary zero vector
  !
  COMPLEX(DP) :: amn
  !! Element of A_mn matrix
  COMPLEX(DP) :: zdotc
  !! <psi_mk|g_n>
  COMPLEX(DP) :: amn_tmp
  !! Element of A_mn matrix
  COMPLEX(DP) :: fac(2)
  !! Factor for spin quantization axis
  COMPLEX(DP), ALLOCATABLE :: sgf(:,:)
  !! Guiding functions
  !
  LOGICAL :: any_uspp
  !! Check if uspp are used
  LOGICAL :: spin_z_pos
  !! Detect if spin quantisation axis is along +z
  LOGICAL :: spin_z_neg
  !! Detect if spin quantisation axis is along -z
  !
  !nocolin: we have half as many projections g(r) defined as wannier
  !         functions. We project onto (1,0) (ie up spin) and then onto
  !         (0,1) to obtain num_wann projections. jry
  !
  any_uspp = ANY( upf(:)%tvanp )
  !
! RM - to remove this once USP is implemented
  IF (any_uspp .and. noncolin) CALL errore('pw2wan90epw',&
      'noncolin calculation not implemented with USP',1)
  !
  ALLOCATE(a_mat(num_bands,n_wannier,iknum))
  ALLOCATE(sgf(npwx,n_proj) )
  ALLOCATE(gf_spinor(2*npwx,n_proj))
  ALLOCATE(sgf_spinor(2*npwx,n_proj))
  !
  ! initialize
  a_mat = czero
  zero_vect = zero
  !
  WRITE (stdout,'(5x,a)') 'AMN'
  !
  IF (any_uspp) THEN
     CALL deallocate_bec_type ( becp )
     CALL allocate_bec_type ( nkb, n_wannier, becp )
     CALL init_us_1
  ENDIF
  !
#if defined(__MPI)
  WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
  ! 
  DO ik = 1, nks
    ! returns in-pool index nkq and absolute index nkq_abs of xk
    CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, nkq_abs )
    ik_g = nkq_abs
    !
    WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
    CALL flush(6)
    ! SP: Replaced by our wrapper to deal with parallel
    CALL readwfc(my_pool_id+1, ik, evc) 
    !
    ! sorts k+G vectors in order of increasing magnitude, up to ecut
    CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), g2kin)
    !
    CALL generate_guiding_functions(ik)   ! they are called gf(npw,n_proj)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (noncolin) THEN
      sgf_spinor = czero
      CALL orient_gf_spinor(npw)
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  USPP
    !
    IF (any_uspp) THEN
       CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
       ! below we compute the product of beta functions with trial func.
       IF (noncolin) THEN
          CALL calbec ( npw, vkb, gf_spinor, becp, n_proj )
       ELSE
          CALL calbec ( npw, vkb, gf, becp, n_proj )
       ENDIF
       ! and we use it for the product S|trial_func>
       IF (noncolin) THEN
          CALL s_psi (npwx, npw, n_proj, gf_spinor, sgf_spinor)
       ELSE
          CALL s_psi (npwx, npw, n_proj, gf, sgf)
       ENDIF
    ELSE
       sgf(:,:) = gf(:,:)
    ENDIF
    !
    IF (noncolin) THEN
      ! we do the projection as g(r)*a(r) and g(r)*b(r)
      DO iw = 1, n_proj
         spin_z_pos=.false.;spin_z_neg=.false.
         ! detect if spin quantisation axis is along z
         IF ((abs(spin_qaxis(1,iw)-0.0d0)<eps6).and.(abs(spin_qaxis(2,iw)-0.0d0)<eps6) &
              .and.(abs(spin_qaxis(3,iw)-1.0d0)<eps6)) THEN
            spin_z_pos=.true.
         ELSEIF (abs(spin_qaxis(1,iw)-0.0d0)<eps6.and.abs(spin_qaxis(2,iw)-0.0d0)<eps6 &
              .and.abs(spin_qaxis(3,iw)+1.0d0)<eps6) THEN
            spin_z_neg=.true.
         ENDIF
         IF (spin_z_pos .or. spin_z_neg) THEN
            ibnd1 = 0
            DO ibnd = 1, nbnd
               IF (excluded_band(ibnd)) CYCLE
               IF (spin_z_pos) THEN
                  ipol = (3-spin_eig(iw))/2
               ELSE
                  ipol = (3+spin_eig(iw))/2
               ENDIF
               istart = (ipol-1)*npwx + 1
               amn = czero
               IF (any_uspp) THEN
                  amn = zdotc(npw, evc(1,ibnd), 1, sgf_spinor(1,iw), 1)
                  amn = amn + zdotc(npw, evc(npwx+1,ibnd), 1, sgf_spinor(npwx+1,iw), 1)
               ELSE
                  amn = zdotc(npw, evc(istart,ibnd), 1, sgf(1,iw), 1)
               ENDIF
               CALL mp_sum(amn, intra_pool_comm)
               ibnd1 = ibnd1 + 1
               !
               ! changed since n_proj is now read from .nnkp file (n_proj=n_wannier) 
               ! a_mat(ibnd1,iw+n_proj*(ipol-1),ik_g) = amn
               a_mat(ibnd1,iw,ik_g) = amn
            ENDDO
         ELSE
           ! general routine
           ! for quantisation axis (a,b,c) 
           ! 'up'    eigenvector is 1/sqrt(1+c) [c+1,a+ib]
           ! 'down'  eigenvector is 1/sqrt(1-c) [c-1,a+ib]
           IF (spin_eig(iw)==1) THEN
              fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw)+1)*cmplx(1.0d0,0.0d0,dp)
              fac(2)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
           ELSE
              fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw))*cmplx(1.0d0,0.0d0,dp)
              fac(2)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
           ENDIF
           ibnd1 = 0
           DO ibnd = 1, nbnd
              IF (excluded_band(ibnd)) CYCLE
              amn = czero
              DO ipol = 1, npol
                 istart = (ipol-1)*npwx + 1
                 amn_tmp = czero
                 IF (any_uspp) THEN 
                    amn_tmp = zdotc(npw, evc(istart,ibnd), 1, sgf_spinor(istart,iw), 1)
                    CALL mp_sum(amn_tmp, intra_pool_comm)
                    amn = amn + amn_tmp
                 ELSE
                    amn_tmp = zdotc(npw, evc(istart,ibnd), 1, sgf(1,iw), 1)
                    CALL mp_sum(amn_tmp, intra_pool_comm)
                    amn = amn + fac(ipol)*amn_tmp
                 ENDIF
              ENDDO
              ibnd1 = ibnd1 + 1
              ! changed since n_proj is now read from .nnkp file (n_proj=n_wannier)
              !a_mat(ibnd1,iw+n_proj*(ipol-1),ik_g) = amn
              a_mat(ibnd1,iw,ik_g) = amn
           ENDDO
         ENDIF ! spin_z_pos
      ENDDO
! -----------
    ELSE ! scalar wavefunction
      DO iw = 1, n_proj
         ibnd1 = 0
         DO ibnd = 1, nbnd
            IF (excluded_band(ibnd)) CYCLE
            amn = czero
            amn = zdotc(npw, evc(1,ibnd), 1, sgf(1,iw), 1)
            CALL mp_sum(amn, intra_pool_comm)
            ibnd1 = ibnd1 + 1
            a_mat(ibnd1,iw,ik_g) = amn
         ENDDO !bands
      ENDDO !wannier fns
    ENDIF
  ENDDO  ! k-points
  !
  DEALLOCATE (sgf)
  DEALLOCATE (csph)
  DEALLOCATE (gf_spinor) 
  DEALLOCATE (sgf_spinor)
  !
  IF (any_uspp) CALL deallocate_bec_type ( becp )
  !
  CALL mp_sum(a_mat, inter_pool_comm)
  !
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a)') 'AMN calculated'
  !
  RETURN
!-----------------------------------------------------------------------
END SUBROUTINE compute_amn_para
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE orient_gf_spinor(npw)
!-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE constants_epw,    ONLY : czero, eps6
  USE noncollin_module, ONLY : npol
  USE wvfct,            ONLY : npwx
  USE wannierEPW,       ONLY : gf, n_proj, spin_qaxis, spin_eig, gf_spinor
  ! 
  IMPLICIT NONE
  ! 
  INTEGER :: npw
  !! Counter on plane waves 
  INTEGER :: iw
  !! Counter on number of projections
  INTEGER  :: ipol
  !! Index of spin-up/spin-down polarizations
  INTEGER ::  istart
  !! Index on plane waves
  LOGICAL :: spin_z_pos
  !! Detect if spin quantisation axis is along +z
  LOGICAL :: spin_z_neg
  !! Detect if spin quantisation axis is along -z
  COMPLEX(DP) :: fac(2)

  gf_spinor = czero
  DO iw = 1, n_proj
     spin_z_pos=.false.;spin_z_neg=.false.
     ! detect if spin quantisation axis is along z
     IF ((abs(spin_qaxis(1,iw)-0.0d0)<eps6).and.(abs(spin_qaxis(2,iw)-0.0d0)<eps6) &
          .and.(abs(spin_qaxis(3,iw)-1.0d0)<eps6)) THEN 
          spin_z_pos=.true.
     ELSEIF (abs(spin_qaxis(1,iw)-0.0d0)<eps6.and.abs(spin_qaxis(2,iw)-0.0d0)<eps6 &
          .and.abs(spin_qaxis(3,iw)+1.0d0)<eps6) then
          spin_z_neg=.true.
     ENDIF
     IF (spin_z_pos .or. spin_z_neg) THEN
        IF (spin_z_pos) THEN
           ipol = (3-spin_eig(iw))/2
        ELSE
           ipol = (3+spin_eig(iw))/2
        ENDIF
        istart = (ipol-1)*npwx + 1
        gf_spinor(istart:istart+npw-1, iw) = gf(1:npw, iw)
     ELSE
        IF(spin_eig(iw)==1) THEN
          fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw)+1)*cmplx(1.0d0,0.0d0,dp)
          fac(2)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
        ELSE
          fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw))*cmplx(1.0d0,0.0d0,dp)
          fac(2)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
        ENDIF
        gf_spinor(1:npw, iw) = gf(1:npw, iw) * fac(1)
        gf_spinor(npwx + 1:npwx + npw, iw) = gf(1:npw, iw) * fac(2)
     ENDIF
  ENDDO
  !
  RETURN
!-----------------------------------------------------------------------
END SUBROUTINE orient_gf_spinor
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE compute_mmn_para
!-----------------------------------------------------------------------
!
!  adapted from compute_mmn in pw2wannier90.f90
!  parallelization on k-points has been added
!  10/2008 Jesse Noffsinger UC Berkeley
!
  USE io_global,       ONLY : stdout, meta_ionode
  USE io_files,        ONLY : diropn, prefix
  USE kinds,           ONLY : DP
  USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
  USE gvecw,           ONLY : ecutwfc
  USE wavefunctions, ONLY : evc, psic, psic_nc
  USE units_lr,        ONLY : lrwfc, iuwfc
  USE fft_base,        ONLY : dffts
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE klist,           ONLY : nkstot, xk, nks, igk_k
  USE gvect,           ONLY : g, ngm, gstart
  USE cell_base,       ONLY : tpiba2, omega, tpiba, bg
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
  USE uspp,            ONLY : nkb, vkb
  USE uspp_param,      ONLY : upf, lmaxq, nh, nhm
  USE becmod,          ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
  USE noncollin_module,ONLY : noncolin, npol
  USE spin_orb,        ONLY : lspinorb
  USE wannierEPW,      ONLY : m_mat, num_bands, nnb, iknum, g_kpb, kpb, ig_, &
                              excluded_band, write_mmn
  USE constants_epw,   ONLY : czero, cone, twopi, zero
  USE io_epw,          ONLY : iummn
#if defined(__NAG)
  USE f90_unix_io,     ONLY : flush
#endif
  USE mp,              ONLY : mp_sum
  USE mp_global,       ONLY : my_pool_id, npool, intra_pool_comm, inter_pool_comm
  ! 
  IMPLICIT NONE
  !
  INTEGER :: mmn_tot 
  !! mmn_tot=iknum*nnb*nbnd*nbnd
  INTEGER :: ik
  !! Counter on k-points
  INTEGER :: ikp
  !! Index of k-point nearest neighbour
  INTEGER :: ib
  !! Counter on b-vectors
  INTEGER :: npwq
  !! Number of planes waves at k+b
  INTEGER :: m
  !! Counter over bands
  INTEGER :: n
  !! Counter over bands
  INTEGER :: ibnd_m
  !! Band index
  INTEGER :: ibnd_n
  !! Band index
  INTEGER :: ih, jh
  !! Counters over the beta functions per atomic type
  INTEGER :: ikb, jkb, ijkb0
  !! Indexes over beta functions
  INTEGER :: na
  !! Counter over atoms
  INTEGER :: nt
  !! Number of type of atoms
  INTEGER :: ind
  !! Counter over nearest neighbours of all k-points
  INTEGER :: nbt
  !!
  INTEGER :: ipool
  !! Index of current pool
  INTEGER ::  nkq
  !! Index of k-point in the pool
  INTEGER ::  nkq_abs
  !! Absolute index of k-point
  INTEGER  :: ipol
  !! Index of spin-up/spin-down polarizations
  INTEGER ::  istart
  !! Index on plane waves ?
  INTEGER :: iend
  !!
  INTEGER ::  ik_g
  !! Temporary index of k-point, ik_g = nkq_abs
  INTEGER :: ikp_g
  !! Temporary index of k+b, ikp_g = ikp
  INTEGER :: ind0 
  !! Starting index for k-point nearest neighbours in each pool
  INTEGER, ALLOCATABLE  :: igkq(:)
  !!
  !
  REAL(DP), ALLOCATABLE :: dxk(:,:)
  !! Temporary vector k+b - k 
  REAL(DP), ALLOCATABLE :: qg(:)
  !!  Square of dxk
  REAL(DP), ALLOCATABLE :: ylm(:,:)
  !!
  REAL(DP) :: xktot(3,nkstot)
  !! Coordinates of k-points
  REAL(DP) :: zero_vect(3)
  !! Temporary zero vector
  REAL(DP) :: arg
  !! 2*pi*(k+b - k)*R
  REAL(DP) :: g_(3)
  !! Temporary vector G_k+b, g_(:) = g_kpb(:,ik,ib)
  !
  COMPLEX(DP), ALLOCATABLE :: phase(:), aux(:), aux_nc(:,:)
  !!
  COMPLEX(DP), ALLOCATABLE :: evcq(:,:)
  !! Wave functions psi_k+b
  COMPLEX(DP), ALLOCATABLE :: Mkb(:,:)
  !! Element of M_mn matrix
  COMPLEX(DP), ALLOCATABLE :: becp2(:,:)
  !! Beta functions
  COMPLEX(DP), ALLOCATABLE :: becp2_nc(:,:,:)
  !! Beta functions, noncolin
  COMPLEX(DP), ALLOCATABLE :: qb(:,:,:,:), qgm(:), qq_so(:,:,:,:)
  !! Local variables for uspp 
  COMPLEX(DP) :: mmn
  !! Element of M_mn matrix
  COMPLEX(DP), ALLOCATABLE :: m_mat_tmp(:,:,:,:)
  !! M_mn matrix
  COMPLEX(DP) :: zdotc
  !! 
  COMPLEX(DP) :: phase1
  !! e^{i*2*pi*(k+b - k)*R}
  !
  LOGICAL :: any_uspp
  !! Check if uspp are used
  LOGICAL :: exst
  !!
  CHARACTER (LEN=80) :: filmmn
  !! file containing M_mn(k,b) matrix
  ! 
  any_uspp = ANY( upf(:)%tvanp )
  !
! RM - to remove this once USP is implemented
  IF (any_uspp .and. noncolin) CALL errore('pw2wan90epw',&
      'noncolin calculation not implimented with USP',1)
  !
  ALLOCATE( phase(dffts%nnr) ) 
  ALLOCATE( igkq(npwx) )
  ALLOCATE( evcq(npol*npwx,nbnd) )
  !
  IF (noncolin) THEN
     ALLOCATE( aux_nc(npwx,npol) )
  ELSE
     ALLOCATE( aux(npwx) )
  ENDIF
  !
  ALLOCATE(m_mat(num_bands,num_bands,nnb,iknum))
  !
  ! close all the wfc files to allow access for each pool to all wfs
  CLOSE (unit = iuwfc,  status = 'keep')
  !
  WRITE (stdout,*)
  WRITE (stdout,'(5x,a)') 'MMN'
  !
  ! Get all the k-vector coords to each pool via xktot
  xktot = 0.d0
  IF (meta_ionode) then
     DO ik = 1, nkstot
        xktot(:,ik) = xk(:,ik)
     ENDDO
  ENDIF
  CALL mp_sum(xktot, inter_pool_comm)
  !
  zero_vect = zero
  m_mat = czero
  !
  mmn_tot = iknum * nnb * nbnd * nbnd
  !
  !   USPP
  !
  IF (any_uspp) THEN
     CALL init_us_1
     CALL allocate_bec_type ( nkb, nbnd, becp )
     IF (noncolin) THEN
        ALLOCATE ( becp2_nc(nkb,2,nbnd) )
     ELSE
        ALLOCATE ( becp2(nkb,nbnd) )
     ENDIF
  ENDIF
  !
  !     qb is  FT of Q(r)
  !
  nbt = nnb * iknum
  !
  ALLOCATE( qg(nbt) )
  ALLOCATE( dxk(3,nbt) )
  !
  ind = 0
  DO ik = 1, iknum ! loop over k-points
     DO ib = 1, nnb ! loop over nearest neighbours for each k-point
        ind = ind + 1
        ikp = kpb(ik,ib) 
        !
        g_(:) = REAL( g_kpb(:,ik,ib) )
        ! bring g_ to cartesian
        CALL cryst_to_cart (1, g_, bg, 1)
        dxk(:,ind) = xktot(:,ikp) + g_(:) - xktot(:,ik) 
        qg(ind) = dxk(1,ind)*dxk(1,ind)+dxk(2,ind)*dxk(2,ind)+dxk(3,ind)*dxk(3,ind)
     ENDDO
  ENDDO
  !
  !  USPP
  !
  IF (any_uspp) THEN
     !
     ALLOCATE( ylm(nbt,lmaxq*lmaxq), qgm(nbt) )
     ALLOCATE( qb(nkb, nkb, ntyp, nbt) )
     ALLOCATE( qq_so(nhm, nhm, 4, ntyp) )
     !
     CALL ylmr2 (lmaxq*lmaxq, nbt, dxk, qg, ylm)
     qg(:) = sqrt(qg(:)) * tpiba
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO ih = 1, nh (nt)
              DO jh = 1, nh (nt)
                 CALL qvan2 (nbt, ih, jh, nt, qg, qgm, ylm)
                 qb (ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     DEALLOCATE (qg, qgm, ylm )
     !
  ENDIF
  !
  ALLOCATE( Mkb(nbnd,nbnd) )
  !
#if defined(__MPI)
  WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
  ! returns in-pool index nkq and absolute index nkq_abs of first k-point in this pool 
  CALL ktokpmq( xk(:, 1), zero_vect, +1, ipool, nkq, nkq_abs )
  ind0 = (nkq_abs-1) * nnb
  !
  ind = ind0
  DO ik = 1, nks 
     ! returns in-pool index nkq and absolute index nkq_abs of xk
     CALL ktokpmq( xk(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
     ik_g = nkq_abs
     !
     WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
     CALL flush(6)
     !
     ! read wfc at k
     CALL readwfc(my_pool_id+1, ik, evc)
     !
     ! sorts k+G vectors in order of increasing magnitude, up to ecut
     CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), g2kin)
     !
     !  USPP
     !
     IF (any_uspp) THEN
        CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
        ! below we compute the product of beta functions with |psi>
        CALL calbec(npw, vkb, evc, becp) 
     ENDIF
     !
     !
     DO ib = 1, nnb  ! loop on finite diff b-vectors
        ind = ind + 1
        !
        ikp = kpb(ik_g,ib)
        ikp_g = ikp 
        !
        CALL ktokpmq( xk(:, ik), xktot(:,ikp_g)-xk(:,ik), +1, ipool, nkq, nkq_abs )
        !
        ! read wfc at k+b
        CALL readwfc( ipool, nkq, evcq )
        !
        CALL gk_sort (xktot(1,ikp_g), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
        !
        ! compute the phase
        phase(:) = czero
        IF ( ig_(ik_g,ib)>0) phase(dffts%nl(ig_(ik_g,ib)) ) = cone
        CALL invfft ('Wave', phase, dffts)
        !
        !  USPP
        !
        IF (any_uspp) THEN
           CALL init_us_2 (npwq, igkq, xk(1,ikp), vkb)
           ! below we compute the product of beta functions with |psi>
           IF (noncolin) THEN
              CALL calbec ( npwq, vkb, evcq, becp2_nc )

              IF (lspinorb) THEN
                 qq_so = czero
                 CALL transform_qq_so(qb(:,:,:,ind), qq_so)
              ENDIF

           ELSE
              CALL calbec ( npwq, vkb, evcq, becp2 )
           ENDIF
        ENDIF
        !
        Mkb(:,:) = czero
        !
        IF (any_uspp) THEN
           ijkb0 = 0
           DO nt = 1, ntyp
              IF ( upf(nt)%tvanp ) THEN
                 DO na = 1, nat
                    !
                    arg = dot_product( dxk(:,ind), tau(:,na) ) * twopi
                    phase1 = cmplx ( COS(arg), -SIN(arg), kind=DP )
                    !
                    IF ( ityp(na) == nt ) THEN
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             !
                             DO m = 1,nbnd
                                IF (excluded_band(m)) CYCLE
                                IF (noncolin) THEN
                                   DO n = 1, nbnd
                                      IF (excluded_band(n)) CYCLE
                                      IF (lspinorb) THEN
                                         Mkb(m,n) = Mkb(m,n) + phase1 * &
                                          ( qq_so(ih,jh,1,nt) * conjg( becp%nc(ikb, 1, m) ) * becp2_nc(jkb, 1, n) &
                                          + qq_so(ih,jh,2,nt) * conjg( becp%nc(ikb, 1, m) ) * becp2_nc(jkb, 2, n) &
                                          + qq_so(ih,jh,3,nt) * conjg( becp%nc(ikb, 2, m) ) * becp2_nc(jkb, 1, n) &
                                          + qq_so(ih,jh,4,nt) * conjg( becp%nc(ikb, 2, m) ) * becp2_nc(jkb, 2, n) )
                                      ELSE
                                         Mkb(m,n) = Mkb(m,n) + phase1 * qb(ih,jh,nt,ind) * &
                                          ( conjg( becp%nc(ikb, 1, m) ) * becp2_nc(jkb, 1, n) &
                                          + conjg( becp%nc(ikb, 2, m) ) * becp2_nc(jkb, 2, n) )
                                      ENDIF
                                   ENDDO
                                ELSE
                                   DO n = 1, nbnd
                                      IF (excluded_band(n)) CYCLE
                                      Mkb(m,n) = Mkb(m,n) + phase1 * qb(ih,jh,nt,ind) * &
                                       conjg( becp%k(ikb,m) ) *becp2(jkb,n)

                                   ENDDO
                                ENDIF
                             ENDDO ! m
                          ENDDO !ih
                       ENDDO !jh
                       ijkb0 = ijkb0 + nh(nt)
                    ENDIF  !ityp
                 ENDDO  !nat 
              ELSE  !tvanp
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                 ENDDO
              ENDIF !tvanp
           ENDDO !ntyp
        ENDIF ! any_uspp
        !
        !
        DO m = 1, nbnd  ! loop on band m
           IF (excluded_band(m)) CYCLE
           !
           IF (noncolin) THEN
              psic_nc(:,:) = czero
              DO ipol = 1, 2!npol
                 istart = (ipol-1)*npwx + 1
                 iend = istart + npw - 1
                 psic_nc(dffts%nl(igk_k (1:npw,ik) ),ipol ) = evc(istart:iend, m)
                 CALL invfft ('Wave', psic_nc(:,ipol), dffts)
                 psic_nc(1:dffts%nnr,ipol) = psic_nc(1:dffts%nnr,ipol) * &
                                               phase(1:dffts%nnr)
                 CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
                 aux_nc(1:npwq,ipol) = psic_nc( dffts%nl(igkq(1:npwq)), ipol )
              ENDDO
           ELSE
              psic(:) = czero
              psic(dffts%nl(igk_k (1:npw,ik) ) ) = evc (1:npw, m)
              CALL invfft ('Wave', psic, dffts)
              psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
              CALL fwfft ('Wave', psic, dffts)
              aux(1:npwq) = psic(dffts%nl(igkq(1:npwq) ) )
           ENDIF
         !  aa = 0.d0
           !
           !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
           !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 > 
           !
           IF (noncolin) THEN
              DO n = 1, nbnd
                 IF (excluded_band(n)) CYCLE
                 mmn = czero
!                DO ipol = 1, 2
!                    mmn = mmn + zdotc(npwq, aux_nc(1,ipol), 1, evcq_nc(1,ipol,n), 1)
                    mmn = mmn + zdotc(npwq, aux_nc(1,1), 1, evcq(1,n),      1) &
                              + zdotc(npwq, aux_nc(1,2), 1, evcq(npwx+1,n), 1)
!                ENDDO
                 CALL mp_sum(mmn, intra_pool_comm)
                 Mkb(m,n) = mmn + Mkb(m,n)
               !  aa = aa + abs(mmn)**2
              ENDDO
           ELSE
              DO n = 1, nbnd
                 IF (excluded_band(n)) CYCLE
                 mmn = zdotc(npwq, aux, 1, evcq(1,n), 1)
                 CALL mp_sum(mmn, intra_pool_comm)
                 Mkb(m,n) = mmn + Mkb(m,n)
               !  aa = aa + abs(mmn)**2
              ENDDO
           ENDIF
        ENDDO   ! m
        !
        ibnd_n = 0
        DO n = 1, nbnd
           IF (excluded_band(n)) CYCLE
           ibnd_n = ibnd_n + 1
           ibnd_m = 0
           DO m = 1, nbnd
              IF (excluded_band(m)) CYCLE
              ibnd_m = ibnd_m + 1
              m_mat(ibnd_m,ibnd_n,ib,ik_g) = Mkb(m,n)
           ENDDO
        ENDDO
        !
     ENDDO !ib
  ENDDO !ik
  !
  CALL mp_sum(m_mat, inter_pool_comm)
  !
  ! RM - write mmn to file
  IF (meta_ionode) THEN
    write_mmn = .true.
    IF (write_mmn) THEN
      ALLOCATE ( m_mat_tmp(nbnd,nbnd,nnb,nkstot) )
      m_mat_tmp = czero
      !
      filmmn=trim(prefix)//'.mmn'
      OPEN (unit = iummn, file = filmmn, form = 'formatted')
      DO ik = 1, nkstot
        DO ib = 1, nnb
          !
          ibnd_n = 0
          DO n = 1, nbnd
            IF (excluded_band(n)) CYCLE
            ibnd_n = ibnd_n + 1
            !
            ibnd_m = 0
            DO m = 1, nbnd
              IF (excluded_band(m)) CYCLE
              ibnd_m = ibnd_m + 1
              !   
              m_mat_tmp(m,n,ib,ik) =  m_mat(ibnd_m,ibnd_n,ib,ik)
              !
            ENDDO
          ENDDO
          !
          DO n = 1, nbnd
            DO m = 1, nbnd
              WRITE (iummn,*) m_mat_tmp(m,n,ib,ik)
            ENDDO
          ENDDO
          !
        ENDDO
      ENDDO
      !
      CLOSE(iummn)
      DEALLOCATE (m_mat_tmp)
    ENDIF
  ENDIF
  !
  DEALLOCATE (Mkb, dxk, phase, evcq, igkq)
  IF (noncolin) THEN
     DEALLOCATE(aux_nc)
  ELSE
     DEALLOCATE(aux)
  ENDIF
  IF (any_uspp) THEN
     DEALLOCATE(qb)
     DEALLOCATE(qq_so)
     CALL deallocate_bec_type (becp)
     IF (noncolin) THEN
        DEALLOCATE(becp2_nc)
     ELSE
        DEALLOCATE(becp2)
     ENDIF
  ENDIF
  !
  WRITE(stdout,'(5x,a)') 'MMN calculated'
  !
  ! reopen wfc here, leaving unit=20 in the same state
  iuwfc = 20
  CALL diropn(iuwfc,'wfc',lrwfc,exst)  
  !
  RETURN
END SUBROUTINE compute_mmn_para
!------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE compute_pmn_para
!-----------------------------------------------------------------------
!
!!
!!  Computes dipole matrix elements.
!!  This can be used to compute the velocity in the local approximation.
!!  The commutator with the non-local psp is neglected.
!!
!!  06/2010  Jesse Noffsinger: adapted from compute_amn_para
!!  08/2016  Samuel Ponce: adapted to work with SOC
!! 
!
  USE io_global,       ONLY : stdout
  USE mp,              ONLY : mp_sum
  USE kinds,           ONLY : DP
  USE klist,           ONLY : xk, nks, igk_k
  USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
  USE gvecw,           ONLY : ecutwfc
  USE wavefunctions,  ONLY : evc
  USE units_lr,        ONLY : lrwfc, iuwfc
  USE gvect,           ONLY : g, ngm
  USE cell_base,       ONLY : tpiba2, tpiba
  USE noncollin_module,ONLY : noncolin
  USE elph2,           ONLY : dmec
  USE constants_epw,   ONLY : czero
  USE uspp_param,      ONLY : upf
  USE becmod,          ONLY : becp, deallocate_bec_type, allocate_bec_type
  USE uspp,            ONLY : nkb
  USE wannierEPW,      ONLY : n_wannier
  USE io_global,    ONLY : meta_ionode
  !
  IMPLICIT NONE
  !  
  LOGICAL :: any_uspp
  !! Check if USPP is present
  !
  INTEGER :: ik
  !! Counter on k-point
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: ig
  !! Counter on G vector
  ! 
  COMPLEX(kind=DP) :: dipole_aux(3,nbnd,nbnd)
  !! Auxilary dipole
  COMPLEX(kind=DP) :: caux 
  !! Wavefunction squared
  !
  any_uspp = ANY( upf(:)%tvanp )
  !
  IF (any_uspp .and. noncolin) CALL errore('pw2wan90epw',&
             'noncolin calculation not implimented with USP',1)
  !
  ALLOCATE( dmec(3,nbnd,nbnd,nks) )
  !
  ! initialize
  dmec = czero
  dipole_aux(:,:,:) = czero
  !
  IF (any_uspp) THEN
    CALL deallocate_bec_type ( becp )
    CALL allocate_bec_type ( nkb, n_wannier, becp )
    CALL init_us_1
  ENDIF
  !
  DO ik = 1, nks
    !
    ! read wfc for the given kpt
    CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
    !
    ! setup k+G grids for each kpt
    CALL gk_sort (xk(:,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(:,ik), g2kin)
    !
    dipole_aux = czero
    DO jbnd = 1,nbnd
      DO ibnd = 1,nbnd
        !
        IF ( ibnd .eq. jbnd ) CYCLE
        !
        ! taken from PP/epsilon.f90 subroutine dipole_calc
        DO ig = 1, npw
          IF (igk_k(ig,ik) > SIZE(g,2) .or. igk_k(ig,ik) < 1) CYCLE
          !
          caux = conjg(evc(ig,ibnd))*evc(ig,jbnd) 
          !
          IF (noncolin) THEN
             !
             caux = caux + conjg(evc(ig+npwx,ibnd))*evc(ig+npwx,jbnd)
             !
          ENDIF
          !
          dipole_aux(:,ibnd,jbnd) = dipole_aux(:,ibnd,jbnd) + &
                                  ( g(:,igk_k(ig,ik)) ) * caux
          !
        ENDDO
        !
      ENDDO !bands i
    ENDDO ! bands j
    ! metal diagonal part
    DO ibnd = 1, nbnd
       DO ig = 1, npw
         IF (igk_k(ig,ik) > SIZE(g,2) .or. igk_k(ig,ik) < 1) CYCLE
         !
         caux = conjg(evc(ig,ibnd))*evc(ig,ibnd) 
         !
         IF (noncolin) THEN
            !
            caux = caux + conjg(evc(ig+npwx,ibnd))*evc(ig+npwx,ibnd)
            !
         ENDIF
         !
         dipole_aux(:,ibnd,ibnd) = dipole_aux(:,ibnd,ibnd) + &
                                         ( g(:,igk_k(ig,ik)) + xk(:,ik) ) * caux
         !
       ENDDO
    ENDDO
    ! need to divide by 2pi/a to fix the units
    dmec(:,:,:,ik) = dipole_aux(:,:,:) * tpiba
    !
  ENDDO  ! k-points
  !
  !
  WRITE(stdout,'(/5x,a)') 'Dipole matrix elements calculated'
  WRITE(stdout,*)
  !DBSP 
  !WRITE(stdout,*) 'dmec ',sum(dmec)
  !IF (meta_ionode) THEN
  !   WRITE(stdout,*) 'dmec(:,:,:,1) ',sum(dmec(:,:,:,1))
  !   WRITE(stdout,*) 'dmec(:,:,:,2) ',sum(dmec(:,:,:,2))
  !ENDIF
  !
  RETURN
END SUBROUTINE compute_pmn_para
!-----------------------------------------------------------------------
!!
!-----------------------------------------------------------------------
SUBROUTINE write_filukk
!-----------------------------------------------------------------------
!
! Here we compute and write out the final ukk matrix which is used by
! epw.x to localize the electron wavefuctions (and therefore the ep-matrix 
! elements)
! 10/2008 Jesse Noffsinger UC Berkeley
! 07/2010 Fixed the rotation for ndimwin when lower bands are not included
!
  USE kinds,        ONLY : DP
  USE io_epw,       ONLY : iuukk
  USE wvfct,        ONLY : nbnd
  USE wannierEPW,   ONLY : n_wannier, iknum, u_mat, u_mat_opt, lwindow, &
                           excluded_band, num_bands
  USE epwcom,       ONLY : filukk
  USE constants_epw,ONLY : czero
  USE io_global,    ONLY : meta_ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! Counter of k-point index
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: ibnd1
  !! Band index
  INTEGER :: iw
  !! Counter on number of Wannier functions
  INTEGER :: ndimwin(iknum)
  !! Number of bands within outer window at each k-point
  !
  COMPLEX(kind=DP), ALLOCATABLE :: u_kc_tmp(:,:,:)
  !! Temporaty rotation matrix
  COMPLEX(kind=DP), ALLOCATABLE :: u_kc(:,:,:)
  !! Rotation matrix
  !
  LOGICAL, ALLOCATABLE :: lwindow_tmp(:,:)
  !! .true. if the band ibnd lies within the outer window at k-point ik 
  !
  ! RM: Band-dimension of u_mat_opt and lwindow is num_bands while
  !     that of u_kc is nbnd to be compatible when reading umat in loadumat. 
  !     Care needs to be taken if exclude_bands is used because
  !     num_bands = nbnd - nexband
  !
  !     u_mat_opt(num_bands,n_wannier,nkstot) in wannier_run
  !     u_mat(n_wannier,n_wannier,nkstot) in wannier_run
  !     lwindow(num_bands,nkstot) in wannier_run
  !     u_kc_tmp(num_bands,n_wannier,nkstot) in write_filukk
  !     u_kc(nbnd,n_wannier,nkstot) in write_filukk
  !
  IF (meta_ionode) THEN
    !
    ndimwin(:) = 0
    DO ik = 1, iknum
       DO ibnd = 1, num_bands
          IF (lwindow(ibnd,ik)) ndimwin(ik) = ndimwin(ik) + 1
       ENDDO
    ENDDO
    !
    ! get the final rotation matrix, which is the product of the optimal
    ! subspace and the rotation among the n_wannier wavefunctions
    !
    ALLOCATE( u_kc_tmp(num_bands, n_wannier, iknum) )
    u_kc_tmp(:,:,:) = czero
    !
    DO ik = 1, iknum
       !
       u_kc_tmp(1:ndimwin(ik),1:n_wannier,ik) = &
            matmul(u_mat_opt(1:ndimwin(ik),:,ik), u_mat(:,1:n_wannier,ik))
       !
    ENDDO
    !
    ALLOCATE( u_kc(nbnd, n_wannier, iknum) )
    u_kc(:,:,:) = czero
    !
    OPEN (unit = iuukk, file = filukk, form = 'formatted')
    ! u_kc(1:num_bands,:,:) = u_kc_tmp(1:num_bands,:,:)
    ! u_kc(num_bands+1:nbnd,:,:) = czero
    DO ik = 1, iknum
       DO iw = 1, n_wannier
          ibnd1 = 0
          DO ibnd = 1, nbnd
             IF (excluded_band(ibnd)) CYCLE
             ibnd1 = ibnd1 + 1
             u_kc(ibnd,iw,ik) = u_kc_tmp(ibnd1,iw,ik)
          ENDDO
       ENDDO
       !
       DO ibnd = 1, nbnd
          DO iw = 1, n_wannier
             WRITE (iuukk,*) u_kc(ibnd,iw,ik)
          ENDDO
       ENDDO
    ENDDO
    !
    ! needs also lwindow when disentanglement is used
    ALLOCATE( lwindow_tmp(nbnd, iknum) )
    lwindow_tmp(:,:) = .false.
    !
    DO ik = 1, iknum
       ibnd1 = 0
       DO ibnd = 1, nbnd
          IF (excluded_band(ibnd)) CYCLE
          ibnd1 = ibnd1 + 1
          lwindow_tmp(ibnd,ik) = lwindow(ibnd1,ik) 
       ENDDO
       !
       DO ibnd = 1, nbnd
          WRITE (iuukk,*) lwindow_tmp(ibnd,ik)
       ENDDO
    ENDDO
    !
    DO ibnd = 1, nbnd
       WRITE(iuukk,*) excluded_band(ibnd)
    ENDDO
    !
    CLOSE (iuukk)
    DEALLOCATE(u_kc_tmp)
    DEALLOCATE(u_kc)
    DEALLOCATE(lwindow_tmp)
    !
  ENDIF
  !
  RETURN
END SUBROUTINE write_filukk
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE phases_a_m
!-----------------------------------------------------------------------
! will set phases here on the matrices.  Should not affect the spreads and
! centers found in w90, but will leave the u_mat_opt and u_mat to reflect the
! known phases
!
!
  USE mp_global,       ONLY : inter_pool_comm
  USE mp,              ONLY : mp_sum
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : meta_ionode
  USE klist,           ONLY : nkstot, xk, nks
  USE wvfct,           ONLY : nbnd
  USE wannierEPW,      ONLY : a_mat, m_mat, num_bands, n_wannier, n_proj, & 
                              nnb, kpb, iknum, excluded_band
  USE elph2,           ONLY : umat, umat_all
  USE constants_epw,   ONLY : czero, cone, zero

  IMPLICIT NONE

  INTEGER :: ik
  !! Counter on k-points
  INTEGER :: ikb
  !! Index of k-point nearest neighbour
  INTEGER :: ib
  !! Counter on b-vectors
  INTEGER :: ipool
  !! Index of current pool
  INTEGER ::  nkq
  !! Index of k-point in the pool
  INTEGER ::  nkq_abs
  !! Absolute index of k-point
  INTEGER ::  ik_g
  !! Temporary index of k-point, ik_g = nkq_abs
  INTEGER :: m
  !! Counter over bands
  INTEGER :: n
  !! Counter over bands
  INTEGER :: ibnd_m
  !! Band index
  INTEGER :: ibnd_n
  !! Band index
  INTEGER :: iw
  !! Counter on number of projections
  !
  REAL(DP) :: xktot(3,nkstot)
  !! Coordinates of k-points
  REAL(DP) :: zero_vect(3)
  !! Temporary zero vector
  !
  COMPLEX(kind=DP), ALLOCATABLE :: a_mat_tmp(:,:,:), a_mat_tmp1(:,:,:)
  !! Temporary a_mat matrices
  COMPLEX(kind=DP), ALLOCATABLE :: m_mn_tmp1(:,:), m_mn_tmp2(:,:), & 
                                   m_mn_tmp3(:,:,:,:), m_mat_tmp(:,:,:,:)
  !! Temporary m_mat matrices
  !
  xktot = 0.d0
  IF (meta_ionode) then
     DO ik = 1, nkstot
        xktot(:,ik) = xk(:,ik)
     ENDDO
  ENDIF
  CALL mp_sum(xktot, inter_pool_comm)
  !
  ! RM: Band-dimension of a_mat and m_mat is num_bands while that of
  !     umat and umat_all is nbnd. This causes a problem if exclude_bands 
  !     is used because num_bands = nbnd - nexband 
  !     a_mat(num_bands,n_wannier,nkstot) in compute_amn_para    
  !     m_mat(num_bands,n_wannier,nnb,nkstot) in compute_mmn_para
  !     umat(nbnd,nbnd,nks) and umat_mat(nbnd,nbnd,nkstot) in setphases_wrap
  !
  ALLOCATE(a_mat_tmp(nbnd,n_wannier,iknum)) 
  ALLOCATE(a_mat_tmp1(nbnd,n_wannier,iknum))
  ALLOCATE(m_mat_tmp(nbnd,nbnd,nnb,iknum))
  ALLOCATE(m_mn_tmp1(nbnd,nbnd))
  ALLOCATE(m_mn_tmp2(nbnd,nbnd))
  ALLOCATE(m_mn_tmp3(nbnd,nbnd,nnb,iknum))
  !
  ! zero all temporary/work quantities
  !
  zero_vect = zero
  a_mat_tmp = czero
  a_mat_tmp1 = czero
  m_mn_tmp1 = czero
  m_mn_tmp2 = czero
  m_mn_tmp3 = czero
  m_mat_tmp = czero
  !
  ! full size a_mat_tmp1 and m_mat_tmp matrices to nbnd bands
  !
  DO ik = 1, nkstot
     DO iw = 1, n_proj
        ibnd_m = 0
        DO m = 1, nbnd
           IF (excluded_band(m)) CYCLE
           ibnd_m = ibnd_m + 1
           a_mat_tmp1(m,iw,ik) = a_mat(ibnd_m,iw,ik)
        ENDDO
     ENDDO
     !
     DO ib = 1, nnb
        ibnd_n = 0
        DO n = 1, nbnd
           IF (excluded_band(n)) CYCLE
           ibnd_n = ibnd_n + 1
           ibnd_m = 0
           DO m = 1, nbnd
              IF (excluded_band(m)) CYCLE
              ibnd_m = ibnd_m + 1
              m_mat_tmp(m,n,ib,ik) = m_mat(ibnd_m,ibnd_n,ib,ik)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DO ik = 1, nks
     ! returns in-pool index nkq and absolute index nkq_abs of xk
     CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, nkq_abs )
     ik_g = nkq_abs
     !
     !  GF_n are the guiding functions which are our initial guesses 
     !  Amn(k) = <psi_k,m|GF_n>.  
     !  We want U(k)^\dagger<psi_k,m|GF_m>
     !
     ! CALL zgemm ('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:,:,ik), & 
     !      nbnd, a_mat(:,:,ik_g), nbnd, czero, a_mat_tmp(:,:,ik_g), nbnd)
     CALL zgemm ('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:,:,ik), &
          nbnd, a_mat_tmp1(:,:,ik_g), nbnd, czero, a_mat_tmp(:,:,ik_g), nbnd)
     !
     DO ib = 1,nnb
        ikb = kpb(ik_g,ib)
        !
        ! Mmn(k,k+b)  = <psi_k_m| psi_(k+b)_n> so we need
        !  (U(k)^\dagger <psi_k_m| ) * (|psi_k+b_n> U(k+b)
        ! = U(k)^\dagger (M_mn) = m_mat_tmp, 
        ! Mmn(k,k+b)' = m_mat_tmp*U(k+b) 
        !
        ! CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umat(:,:,ik), & 
        !      nbnd, m_mat(:,:,ib,ik_g), nbnd, czero, m_mn_tmp1(:,:), nbnd)
        CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umat(:,:,ik), &
             nbnd, m_mat_tmp(:,:,ib,ik_g), nbnd, czero, m_mn_tmp1(:,:), nbnd)
        CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, m_mn_tmp1(:,:), & 
             nbnd, umat_all(:,:,ikb), nbnd, czero, m_mn_tmp2(:,:), nbnd)
        ! 
        ! m_mn_tmp1 = matmul ( conjg( transpose (umat(:,:,ik) )), m_mat(:,:,ib,ik_g ) )
        ! m_mn_tmp2 = matmul ( m_mn_tmp1, umat_g(:,:,ikb) )
        !
        m_mn_tmp3(:,:,ib,ik_g) = m_mn_tmp2(:,:)
     ENDDO
  ENDDO
  CALL mp_sum(a_mat_tmp, inter_pool_comm)
  CALL mp_sum(m_mn_tmp3, inter_pool_comm)
  !
  ! a_mat(:,:,:) = a_mat_tmp(:,:,:)
  ! m_mat(:,:,:,:) = m_mn_tmp3(:,:,:,:)
  !
  ! slim down a_mat and m_mat matrices to num_bands=nbnd-nexband bands
  ! 
  DO ik = 1, nkstot
     DO iw = 1, n_proj
        ibnd_m = 0
        DO m = 1, nbnd
           IF (excluded_band(m)) CYCLE
           ibnd_m = ibnd_m + 1
           a_mat(ibnd_m,iw,ik) = a_mat_tmp(m,iw,ik)
        ENDDO
     ENDDO
     !
     DO ib = 1, nnb
        ibnd_n = 0
        DO n = 1, nbnd
           IF (excluded_band(n)) CYCLE
           ibnd_n = ibnd_n + 1
           ibnd_m = 0
           DO m = 1, nbnd
              IF (excluded_band(m)) CYCLE
              ibnd_m = ibnd_m + 1
              m_mat(ibnd_m,ibnd_n,ib,ik) = m_mn_tmp3(m,n,ib,ik) 
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(a_mat_tmp)
  DEALLOCATE(a_mat_tmp1)
  DEALLOCATE(m_mat_tmp)
  DEALLOCATE(m_mn_tmp1)
  DEALLOCATE(m_mn_tmp2)
  DEALLOCATE(m_mn_tmp3)
  !
  RETURN
END SUBROUTINE phases_a_m
!-----------------------------------------------------------------------
!---------------------------------------
!
! SUBROUTINES below are largely unchanged
! from pw2wannier90.x of esp-4.0.1
!
!---------------------------------------
!
SUBROUTINE generate_guiding_functions(ik)
  !
  USE kinds,          ONLY : DP
  USE mp_global,      ONLY : intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE constants,      ONLY : tpi, eps8
  USE wvfct,          ONLY : npw
  USE gvect,          ONLY : g
  USE cell_base,      ONLY : tpiba
  USE wannierEPW,     ONLY : n_proj, gf, center_w, csph, alpha_w, &
                             r_w
  USE klist,          ONLY : xk, igk_k 
  USE constants_epw,  ONLY : czero

  IMPLICIT NONE

  INTEGER, PARAMETER :: lmax=3, lmax2=(lmax+1)**2
  INTEGER :: iw, ig, ik, l
  INTEGER :: lm, iig
  REAL(DP) :: arg, anorm
  COMPLEX(DP) :: zdotc, lphase
  REAL(DP), ALLOCATABLE :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:) 
  !
  ALLOCATE( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax) )
  !
  DO ig = 1, npw
     gk (1,ig) = xk(1, ik) + g(1, igk_k(ig,ik) )
     gk (2,ig) = xk(2, ik) + g(2, igk_k(ig,ik) )
     gk (3,ig) = xk(3, ik) + g(3, igk_k(ig,ik) )
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  ENDDO

  CALL ylmr2 (lmax2, npw, gk, qg, ylm)
  ! define qg as the norm of (k+g) in a.u.
  qg(:) = sqrt(qg(:)) * tpiba
  ! RM changed according to QE4.0.3/PP/pw2wannier90 
  DO iw = 1, n_proj
  !DO iw = 1, n_wannier
     !
     gf(:,iw) = czero
     !
     CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial) 
     !
     DO lm = 1, lmax2
        IF ( abs(csph(lm,iw)) < eps8 ) CYCLE
        l = int (sqrt( lm-1.d0))
        lphase = (0.d0,-1.d0)**l
        !
        DO ig=1,npw
           gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
        ENDDO !ig
     ENDDO ! lm
     !
     DO ig=1,npw
        iig = igk_k(ig,ik)
        arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                          gk(3,ig)*center_w(3,iw) ) * tpi
        ! center_w are cartesian coordinates in units of alat 
        sk(ig) = CMPLX(cos(arg), -sin(arg) , kind=DP)
        gf(ig,iw) = gf(ig,iw) * sk(ig) 
     ENDDO
     anorm = REAL(ZDOTC(npw,gf(1,iw),1,gf(1,iw),1))
     CALL mp_sum(anorm, intra_pool_comm)
     gf(:,iw) = gf(:,iw) / dsqrt(anorm)
  ENDDO
  !
  DEALLOCATE ( gk, qg, ylm, sk, radial)
  !
  RETURN
END SUBROUTINE generate_guiding_functions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_band
!-----------------------------------------------------------------------
  USE wvfct,         ONLY : nbnd, et
  USE constants_epw, ONLY : ryd2ev
  USE constants,     ONLY : rytoev
  USE wannierEPW,    ONLY : ikstart, ikstop, iknum, num_bands, eigval, &
                            excluded_band
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! Counter on k-point
  INTEGER ikevc
  !! k-point index
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ibnd1
  !! Band index
  !
  ALLOCATE(eigval(num_bands,iknum))
  !
  DO ik = ikstart, ikstop
     ikevc = ik - ikstart + 1
     ibnd1 = 0
     DO ibnd = 1, nbnd
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
! RM - same value for rytoev as in wannier90
!        eigval(ibnd1,ikevc) = et(ibnd,ik)*ryd2ev
        eigval(ibnd1,ikevc) = et(ibnd,ik)*rytoev
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE write_band
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_plot
!-----------------------------------------------------------------------
!
! JN 06/2009: 
! added a couple of calls -- now works with multiple
! pools/procs (but one proc per pool)
!
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, meta_ionode
  USE io_epw,          ONLY : iun_plot
  USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
  USE gvecw,           ONLY : ecutwfc
  USE wavefunctions, ONLY : evc, psic, psic_nc
  USE wannierEPW,      ONLY : reduce_unk, wvfn_formatted, ispinw, nexband, &
                              excluded_band 
  USE klist,           ONLY : xk, nks, igk_k
  USE gvect,           ONLY : g, ngm 
  USE cell_base,       ONLY : tpiba2
  USE fft_base,        ONLY : dffts
  USE fft_interfaces,  ONLY : invfft
  USE noncollin_module,ONLY : noncolin, npol
  USE scatter_mod,     ONLY : gather_grid
  USE constants_epw,   ONLY : czero, zero
  USE mp_global,       ONLY : my_pool_id
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! Counter on k-point
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ibnd1
  !! Band index
  INTEGER :: ipool
  !! Index of current pool
  INTEGER ::  nkq
  !! Index of k-point in the pool
  INTEGER ::  nkq_abs
  !! Absolute index of k-point
  INTEGER ::  ik_g
  !! Temporary index of k-point, ik_g = nkq_abs
  INTEGER :: spin
  !! Spin case
  INTEGER :: ipol
  !! Counter on polarizations
  !
  REAL(kind=DP) :: zero_vect(3)
  !
  CHARACTER(len=20) wfnname
  !
  ! aam: 1/5/06: for writing smaller unk files 
  INTEGER :: n1by2, n2by2, n3by2, i, j, k, idx, pos
  COMPLEX(DP), ALLOCATABLE :: psic_small(:), psic_nc_small(:,:)   
  !-------------------------------------------!

#if defined(__MPI)
  INTEGER nxxs
  COMPLEX(DP),ALLOCATABLE :: psic_all(:), psic_nc_all(:,:)
  nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
  IF (.NOT.noncolin) THEN
     ALLOCATE(psic_all(nxxs) )
  ELSE
     ALLOCATE(psic_nc_all(nxxs,npol) )
  ENDIF
#endif
  !
  CALL start_clock( 'write_unk' )
  ! 
  zero_vect = zero
  WRITE (stdout,*) '    Writing out UNK plot files'
  !
  IF (reduce_unk) then
     WRITE(stdout,'(3(a,i5))') 'nr1s =',dffts%nr1,'nr2s=',dffts%nr2,'nr3s=',dffts%nr3
     n1by2=(dffts%nr1+1)/2
     n2by2=(dffts%nr2+1)/2
     n3by2=(dffts%nr3+1)/2
     WRITE(stdout,'(3(a,i5))') 'n1by2=',n1by2,'n2by2=',n2by2,'n3by2=',n3by2
     IF (.NOT.noncolin) THEN                                                               
        ALLOCATE(psic_small(n1by2*n2by2*n3by2))                                            
       psic_small = czero                                                      
     ELSE                                                                                  
        ALLOCATE(psic_nc_small(n1by2*n2by2*n3by2,npol))                                    
        psic_nc_small = czero                                                   
     ENDIF 
  ENDIF
  !
  DO ik = 1, nks
     ! returns in-pool index nkq and absolute index nkq_abs of xk
     CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, nkq_abs )
     ik_g = nkq_abs
     !
     spin = ispinw
     IF (ispinw.eq.0) spin = 1
     IF (.NOT.noncolin) THEN                                                               
        WRITE(wfnname,200) ik_g, spin                                                     
     ELSE                                                                                  
        WRITE(wfnname,201) ik_g                                                           
     ENDIF                                                                                 
201  FORMAT ('UNK',i5.5,'.','NC')                                                          
200  FORMAT ('UNK',i5.5,'.',i1)
     !
     IF (meta_ionode) THEN 
        IF (wvfn_formatted) THEN 
           OPEN (unit=iun_plot, file=wfnname,form='formatted')
           IF (reduce_unk) THEN
              WRITE(iun_plot,*)  n1by2, n2by2, n3by2, ik_g, nbnd-nexband
           ELSE
              WRITE(iun_plot,*)  dffts%nr1, dffts%nr2, dffts%nr3, ik_g, nbnd-nexband
           ENDIF
        ELSE
           OPEN (unit=iun_plot, file=wfnname,form='unformatted')
           IF (reduce_unk) THEN
              WRITE(iun_plot)  n1by2, n2by2, n3by2, ik_g, nbnd-nexband
           ELSE
              WRITE(iun_plot)  dffts%nr1, dffts%nr2, dffts%nr3, ik_g, nbnd-nexband
           ENDIF
        ENDIF
     ENDIF
     !
     CALL readwfc(my_pool_id+1, ik, evc)
     CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), g2kin)
     !
     ibnd1 = 0
     DO ibnd = 1, nbnd
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
        IF (.NOT.noncolin) THEN
           psic(:) = czero
           psic(dffts%nl (igk_k (1:npw,ik) ) ) = evc (1:npw, ibnd)
           CALL invfft ('Wave', psic, dffts)
        ELSE
           psic_nc(:,:) = czero
           DO ipol = 1, npol                                                               
              psic_nc(dffts%nl (igk_k (1:npw,ik) ), ipol) = evc (1+npwx*(ipol-1):npw+npwx*(ipol-1), ibnd) 
              CALL invfft ('Wave', psic_nc(:,ipol), dffts)                                 
           ENDDO 
        ENDIF
        IF (reduce_unk) pos = 0
#if defined(__MPI)
        IF (.NOT.noncolin) THEN
           CALL gather_grid(dffts,psic,psic_all)
        ELSE
           DO ipol = 1, npol
              CALL gather_grid(dffts,psic_nc(:,ipol),psic_nc_all(:,ipol))
           ENDDO
        ENDIF
        !
        IF (reduce_unk) THEN
           DO k = 1, dffts%nr3, 2
              DO j = 1, dffts%nr2, 2
                 DO i = 1, dffts%nr1, 2
                    idx = (k-1)*dffts%nr3*dffts%nr2 + (j-1)*dffts%nr2 + i
                    pos = pos + 1
                    IF (.NOT.noncolin) THEN
                       psic_small(pos) = psic_all(idx)
                    ELSE
                       DO ipol = 1, npol
                          psic_nc_small(pos,ipol) = psic_nc_all(idx,ipol)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
        !
        IF (meta_ionode) THEN
           IF (wvfn_formatted) THEN
              IF (reduce_unk) THEN
                 IF (.NOT.noncolin) THEN
                    WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
                 ELSE
                    DO ipol = 1, npol
                       WRITE (iun_plot,'(2ES20.10)') (psic_nc_small(j,ipol),j=1,n1by2*n2by2*n3by2)
                    ENDDO
                 ENDIF
              ELSE
                 IF (.NOT.noncolin) THEN
                    WRITE (iun_plot,'(2ES20.10)') (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                 ELSE
                    DO ipol = 1, npol
                       WRITE (iun_plot,'(2ES20.10)') (psic_nc_all(j,ipol),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                    ENDDO
                 ENDIF
              ENDIF
           ELSE
              IF (reduce_unk) THEN
                 IF (.NOT.noncolin) THEN
                    WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
                 ELSE
                    DO ipol = 1, npol
                       WRITE (iun_plot) (psic_nc_small(j,ipol),j=1,n1by2*n2by2*n3by2)
                    ENDDO
                 ENDIF
              ELSE
                 IF (.NOT.noncolin) THEN
                    WRITE (iun_plot) (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                 ELSE
                    DO ipol = 1, npol
                       WRITE (iun_plot) (psic_nc_all(j,ipol),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
#else
        IF (reduce_unk) THEN
           DO k = 1, dffts%nr3, 2
              DO j = 1, dffts%nr2, 2
                 DO i = 1, dffts%nr1, 2
                    idx = (k-1)*dffts%nr3*dffts%nr2 + (j-1)*dffts%nr2 + i
                    pos = pos + 1
                    IF (.NOT.noncolin) THEN
                       psic_small(pos) = psic(idx)
                    ELSE
                       DO ipol = 1, npol
                          psic_nc_small(pos,ipol) = psic_nc(idx,ipol)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
        !
        IF (meta_ionode) THEN
           IF (wvfn_formatted) THEN
              IF (.NOT.noncolin) THEN
                 IF (reduce_unk) THEN
                    WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
                 ELSE
                    WRITE (iun_plot,'(2ES20.10)') (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                 ENDIF
              ELSE
                 DO ipol = 1, npol
                    IF (reduce_unk) THEN
                       WRITE (iun_plot,'(2ES20.10)') (psic_nc_small(j,ipol),j=1,n1by2*n2by2*n3by2)
                    ELSE
                       WRITE (iun_plot,'(2ES20.10)') (psic_nc(j,ipol),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                    ENDIF
                 ENDDO
              ENDIF
           ELSE
              IF (.NOT.noncolin) THEN
                 IF (reduce_unk) THEN
                    WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
                 ELSE
                    WRITE (iun_plot) (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                 ENDIF
              ELSE
                 DO ipol = 1, npol
                    IF (reduce_unk) THEN
                       WRITE (iun_plot) (psic_nc_small(j,ipol),j=1,n1by2*n2by2*n3by2)
                    ELSE
                       WRITE (iun_plot) (psic_nc(j,ipol),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
                    ENDIF
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
#endif
     ENDDO !ibnd
     !
     IF (meta_ionode) CLOSE (iun_plot)
     !
  ENDDO  !ik
  !
  IF (reduce_unk) THEN
     IF (.NOT.noncolin) THEN
        DEALLOCATE(psic_small)
     ELSE
        DEALLOCATE(psic_nc_small)
     ENDIF
  ENDIF
  !
#if defined(__MPI)
  IF (.NOT.noncolin) THEN
     DEALLOCATE( psic_all )
  ELSE
     DEALLOCATE( psic_nc_all )
  ENDIF
#endif
  !
  WRITE(stdout,'(/)')
  WRITE(stdout,*) ' UNK written'
  !
  CALL stop_clock( 'write_unk' )
  !
  RETURN
END SUBROUTINE write_plot
!---------------------------------------
SUBROUTINE ylm_expansion 
  USE kinds,         ONLY : DP
  USE random_numbers,ONLY : randy
  USE wannierEPW,       ONLY : n_proj , xaxis, zaxis, csph, l_w, mr_w
  USE matrix_inversion, ONLY: invmat
  !
  IMPLICIT NONE
  !
  ! local variables
  INTEGER, PARAMETER :: lmax2=16
  INTEGER ::  i, ir, iw
  REAL(DP), ALLOCATABLE :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
  REAL(DP) :: u(3,3)
  !
  ALLOCATE (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2))
  ALLOCATE (ylm(lmax2,lmax2), mly(lmax2,lmax2))

  ! generate a set of nr=lmax2 random vectors
  DO ir = 1, lmax2
     DO i = 1,3
        r(i,ir) = randy() -0.5d0
     ENDDO
  ENDDO
  rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
  !- compute ylm(ir,lm)
  CALL ylmr2(lmax2, lmax2, r, rr, ylm)
  !- store the inverse of ylm(ir,lm) in mly(lm,ir)
  CALL invmat(lmax2, ylm, mly)
  !- check that r points are independent
  CALL check_inverse(lmax2, ylm, mly)
  !
  DO iw = 1, n_proj 
     !
     !- define the u matrix that rotate the reference frame
     CALL set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
     !- find rotated r-vectors 
     rp(:,:) = matmul ( u(:,:) , r(:,:) )
     !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
     CALL ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2) 
     !
     csph(:,iw) = matmul (mly(:,:), ylm_w(:))
     !
!     WRITE (stdout,*)
!     WRITE (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
!     WRITE (stdout,'(16i6)')   (lm, lm=1,lmax2)
!     WRITE (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)
     !
  ENDDO
  DEALLOCATE (r, rp, rr, ylm_w, ylm, mly )
  !
  RETURN
END SUBROUTINE ylm_expansion
!-----------------------------------------------
SUBROUTINE check_inverse(lmax2, ylm, mly)
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER :: lmax2
  REAL(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
  ! local variables
  REAL(DP), ALLOCATABLE :: uno(:,:)
  REAL(DP) :: capel
  INTEGER :: lm
  !
  ALLOCATE (uno(lmax2,lmax2) )
  uno = matmul(mly, ylm)
  capel = 0.d0
  DO lm = 1, lmax2
     uno(lm,lm) = uno(lm,lm) - 1.d0
  ENDDO
  capel = capel + SUM ( abs(uno(1:lmax2,1:lmax2) ) )
!  WRITE (stdout,*) "capel = ", capel
  IF (capel > eps8) CALL errore('ylm_expansion', &
                   ' inversion failed: r(*,1:nr) are not all independent !!',1)
  DEALLOCATE (uno)
  !
  RETURN
END SUBROUTINE check_inverse
!--------------------------------------------   
SUBROUTINE set_u_matrix(x,z,u)
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  !
  ! I/O variables
  REAL(DP) :: x(3),z(3),u(3,3)
  ! local variables
  REAL(DP) :: xx, zz, y(3), coseno
  !
  xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
  IF (xx < eps8) CALL errore ('set_u_matrix',' |xaxis| < eps ',1)
!  x(:) = x(:)/xx
  zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
  IF (zz < eps8) CALL errore ('set_u_matrix',' |zaxis| < eps ',1)
!  z(:) = z(:)/zz
  !
  coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
  IF (abs(coseno) > eps8) CALL errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)
  !
  y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
  y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
  y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz
  !
  u(1,:) = x(:)/xx
  u(2,:) = y(:)
  u(3,:) = z(:)/zz
  !
!  WRITE (stdout,'(3f10.7)') u(:,:)
  !
  RETURN
   !
END SUBROUTINE set_u_matrix
!--------------------------------------------
SUBROUTINE ylm_wannier(ylm,l,mr,r,nr) 
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr) 
! of the spherical harmonic identified  by indices (l,mr) 
! in table 3.1 of the wannierf90 specification.
! 
! No reference to the particular ylm ordering internal to quantum-espresso
! is assumed. 
!
! If ordering in wannier90 code is changed or extended this should be the 
! only place to be modified accordingly
!
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi, eps8
  !
  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER :: l, mr, nr
  REAL(DP) :: ylm(nr), r(3,nr)
  !
  ! local variables
  REAL(DP), EXTERNAL :: s, p_z, px, py, dz2, dxz, dyz, dx2my2, dxy, &
                       fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
  REAL(DP) :: rr, cost, phi
  INTEGER :: ir
  REAL(DP) :: bs2, bs3, bs6, bs12
  bs2 = 1.d0/sqrt(2.d0)
  bs3 = 1.d0/sqrt(3.d0)
  bs6 = 1.d0/sqrt(6.d0)
  bs12 = 1.d0/sqrt(12.d0)
  !
  IF (l > 3 .OR. l < -5 ) CALL errore('ylm_wannier',' l out of range ', 1)
  IF (l>=0) THEN
     IF (mr < 1 .OR. mr > 2*l+1) CALL errore('ylm_wannier','mr out of range' ,1)
  ELSE
     IF (mr < 1 .OR. mr > abs(l)+1 ) CALL errore('ylm_wannier','mr out of range',1)
  ENDIF
  !
  DO ir =1, nr
     rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
     IF (rr < eps8) CALL errore('ylm_wannier',' rr too small ',1)
     !
     cost =  r(3,ir) / rr
     !
     !  beware the arc tan, it is defined modulo pi
     !
     IF (r(1,ir) > eps8) THEN
        phi = atan( r(2,ir)/r(1,ir) )
     ELSEIF (r(1,ir) < -eps8 ) THEN
        phi = atan( r(2,ir)/r(1,ir) ) + pi
     ELSE
        phi = sign( pi/2.d0,r(2,ir) )
     ENDIF
     !
     IF (l==0) THEN   ! s orbital
                   ylm(ir) = s()  
     ENDIF
     IF (l==1) THEN   ! p orbitals
        IF (mr==1) ylm(ir) = p_z(cost) 
        IF (mr==2) ylm(ir) = px(cost,phi)
        IF (mr==3) ylm(ir) = py(cost,phi)
     ENDIF
     IF (l==2) THEN   ! d orbitals
        IF (mr==1) ylm(ir) = dz2(cost)
        IF (mr==2) ylm(ir) = dxz(cost,phi)
        IF (mr==3) ylm(ir) = dyz(cost,phi)
        IF (mr==4) ylm(ir) = dx2my2(cost,phi)
        IF (mr==5) ylm(ir) = dxy(cost,phi)
     ENDIF
     IF (l==3) THEN   ! f orbitals
        IF (mr==1) ylm(ir) = fz3(cost)
        IF (mr==2) ylm(ir) = fxz2(cost,phi)
        IF (mr==3) ylm(ir) = fyz2(cost,phi)
        IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
        IF (mr==5) ylm(ir) = fxyz(cost,phi)
        IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
        IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
     ENDIF
     IF (l==-1) THEN  !  sp hybrids
        IF (mr==1) ylm(ir) = bs2 * ( s() + px(cost,phi) ) 
        IF (mr==2) ylm(ir) = bs2 * ( s() - px(cost,phi) ) 
     ENDIF
     IF (l==-2) THEN  !  sp2 hybrids 
        IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
        IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
        IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
     ENDIF
     IF (l==-3) THEN  !  sp3 hybrids
        IF (mr==1) ylm(ir) = 0.5d0*(s()+px(cost,phi)+py(cost,phi)+p_z(cost))
        IF (mr==2) ylm(ir) = 0.5d0*(s()+px(cost,phi)-py(cost,phi)-p_z(cost))
        IF (mr==3) ylm(ir) = 0.5d0*(s()-px(cost,phi)+py(cost,phi)-p_z(cost))
        IF (mr==4) ylm(ir) = 0.5d0*(s()-px(cost,phi)-py(cost,phi)+p_z(cost))
     ENDIF
     IF (l==-4) THEN  !  sp3d hybrids
        IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
        IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
        IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
        IF (mr==4) ylm(ir) = bs2*p_z(cost)+bs2*dz2(cost)
        IF (mr==5) ylm(ir) =-bs2*p_z(cost)+bs2*dz2(cost)
     ENDIF
     IF (l==-5) THEN  ! sp3d2 hybrids
        IF (mr==1) ylm(ir) = bs6*s()-bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
        IF (mr==2) ylm(ir) = bs6*s()+bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
        IF (mr==3) ylm(ir) = bs6*s()-bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
        IF (mr==4) ylm(ir) = bs6*s()+bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
        IF (mr==5) ylm(ir) = bs6*s()-bs2*p_z(cost)+bs3*dz2(cost)
        IF (mr==6) ylm(ir) = bs6*s()+bs2*p_z(cost)+bs3*dz2(cost)
     ENDIF
     !
  ENDDO
  !
  RETURN
END SUBROUTINE ylm_wannier

!======== l = 0 =====================================================================
FUNCTION s()
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) :: s
  s = 1.d0/ sqrt(fpi)
  RETURN
END FUNCTION s
!======== l = 1 =====================================================================
FUNCTION p_z(cost)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::p_z, cost
  p_z =  sqrt(3.d0/fpi) * cost
  RETURN
END FUNCTION p_z
FUNCTION px(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::px, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  px =  sqrt(3.d0/fpi) * sint * cos(phi)
  RETURN
END FUNCTION px
FUNCTION py(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::py, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  py =  sqrt(3.d0/fpi) * sint * sin(phi)
  RETURN
END FUNCTION py
!======== l = 2 =====================================================================
FUNCTION dz2(cost)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::dz2, cost
  dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
  RETURN
END FUNCTION dz2
FUNCTION dxz(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::dxz, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
  RETURN
END FUNCTION dxz
FUNCTION dyz(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::dyz, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
  RETURN
END FUNCTION dyz
FUNCTION dx2my2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::dx2my2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
  RETURN
END FUNCTION dx2my2
FUNCTION dxy(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : fpi
  IMPLICIT NONE
  REAL(DP) ::dxy, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
  RETURN
END FUNCTION dxy
!======== l = 3 =====================================================================
FUNCTION fz3(cost)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fz3, cost
  fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
  RETURN
END FUNCTION fz3
FUNCTION fxz2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fxz2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
  RETURN
END FUNCTION fxz2
FUNCTION fyz2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fyz2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
  RETURN
END FUNCTION fyz2
FUNCTION fzx2my2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fzx2my2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
  RETURN
END FUNCTION fzx2my2
FUNCTION fxyz(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fxyz, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
  RETURN
END FUNCTION fxyz
FUNCTION fxx2m3y2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fxx2m3y2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
  RETURN
END FUNCTION fxx2m3y2
FUNCTION fy3x2my2(cost,phi)
  USE kinds, ONLY :  DP
  USE constants_epw, ONLY : pi
  IMPLICIT NONE
  REAL(DP) ::fy3x2my2, cost, phi, sint
  sint = sqrt(abs(1.d0 - cost*cost))
  fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
  RETURN
END FUNCTION fy3x2my2
!
!-----------------------------------------------------------------------
SUBROUTINE radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the radial functions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  !
  IMPLICIT NONE
  ! I/O
  INTEGER :: ng, rvalue, lmax
  REAL(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  REAL(DP), PARAMETER :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  REAL(DP) :: rad_int, pref, x
  INTEGER :: l, ir, ig, mesh_r
  REAL(DP), ALLOCATABLE :: bes(:), func_r(:), r(:), rij(:), aux(:)
  ! 
  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  ALLOCATE ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
  ALLOCATE ( aux(mesh_r))
  !
  !    compute the radial mesh
  !
  DO ir = 1, mesh_r
     x = xmin  + DBLE (ir - 1) * dx 
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  ENDDO
  !
  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * & 
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(2.0d0/3.0d0) * &
                     (1.d0 - 1.5d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(:)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  DO l = 0, lmax
     DO ig=1,ng
       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:) * r(:)
       CALL simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     ENDDO
  ENDDO
  !
  DEALLOCATE (bes, func_r, r, rij, aux )
  !
  RETURN
END SUBROUTINE radialpart

