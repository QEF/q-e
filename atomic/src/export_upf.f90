!
! Copyright (C) 2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE export_upf(iunps, unit_loc)
  !---------------------------------------------------------------------
  !
  use constants, only : fpi
  use kinds, only : dp
  use radial_grids, only : radial_grid_COPY, nullify_radial_grid, &
                           deallocate_radial_grid
  use ld1inc, only : author, nlcc, zval, lpaw, write_coulomb, &
                     etots, rel, ecutwfc, ecutrho, iswitch, &
                     nwfts, nbeta, lmax, which_augfun, elts, octs, llts, &
                     nnts, rcutusts, rcutts, rcut, rcutus, els, ikk, nwfs, &
                     lls, nns, ocs, beta, bmat, qq, qvan, qvanl, rcloc, lloc, &
                     betas, grid, rhos, phits, psipaw, vpsloc, phis, &
                     rmatch_augfun, etot, etots, jjs, pawsetup, nn, &
                     core_state, ll, el, nwf, psi, vpot, nconf, zed, &
                     jjts, vpstot, lltsc, rcuttsc, rcutustsc, eltsc, &
                     lsave_wfc, wfc_ae_recon, wfc_ps_recon, tm, enlts, &
                     nstoaets, pseudotype, enls, rhoc, vnl, vpsloc, &
                     lgipaw_reconstruction, use_paw_as_gipaw
  use funct, only: get_dft_name
  use global_version, only: version_number, svn_revision
  !
  use pseudo_types
  use write_upf_v2_module, only: write_upf_v2, &
                                 pseudo_config, deallocate_pseudo_config
  !
  implicit none
  !
  !CHARACTER(len=*),INTENT(IN) :: filename
  INTEGER,INTENT(IN)::iunps, unit_loc
  !
  integer :: ibeta, jbeta, kbeta, l, ind, l1, l2
  !
  !     Local variables
  !
  integer :: nb, mesh
  TYPE (pseudo_upf)              :: upf
  TYPE (pseudo_config)           :: at_conf
  TYPE (radial_grid_type),TARGET :: internal_grid
  CHARACTER(len=2), external :: atom_name
  CHARACTER(len=9) :: day, hour

  call date_and_tim(day,hour)
  !
  CALL nullify_pseudo_upf( upf )
  CALL nullify_radial_grid( internal_grid )
  !
  IF (iswitch < 4 ) THEN
     upf%generated='Generated using "atomic" code by A. Dal Corso &
                  & v.' // TRIM (version_number)
    IF ( TRIM (svn_revision) /= 'unknown' ) upf%generated = &
         TRIM (upf%generated) // ' svn rev. ' // TRIM (svn_revision)
 
  ELSE IF (iswitch==4) THEN
     upf%generated='Generated using LDA-1/2 implemented by Leonardo&
                  & Matheus Marion Jorge'
  ENDIF
  upf%author=trim(author)
  upf%date=trim(day)
  upf%nv = "2.0.1" ! format version
  !
  upf%zp   = zval
  upf%nlcc = nlcc
  upf%dft  = get_dft_name()
  upf%psd  = atom_name(nint(zed))

  if( pseudotype == 3) then
     upf%tvanp = .true.
     upf%typ='USPP'
  else
     upf%tvanp = .false.
     upf%typ='NC'
  endif
  if(lpaw)          upf%typ='PAW'
  if(write_coulomb) upf%typ='1/r'

  upf%tpawp = lpaw
  upf%tcoulombp = write_coulomb
  upf%has_gipaw = lgipaw_reconstruction
  upf%paw_as_gipaw = use_paw_as_gipaw
  upf%etotps = etots
  upf%has_so = (rel == 2)
  IF (rel == 2) THEN
      upf%rel='full'
  ELSE IF (rel == 1) THEN
      upf%rel='scalar'
  ELSE IF (rel < 1) THEN
      upf%rel='no'
  ELSE
      call errore('export_upf', 'Unknown relativistic',1)
  ENDIF
  !
  upf%ecutwfc = ecutwfc
  upf%ecutrho = max(ecutrho, ecutwfc*4._dp)
  !
  upf%nwfc = nwfts 
  upf%nbeta = nbeta
  !
  if (.not. lpaw) then
   upf%lmax = lmax
   upf%q_with_l = (which_augfun == 'PSQ')
  else
   upf%lmax = pawsetup%lmax
   upf%q_with_l = .true.
  endif
  upf%lmax_rho = 2*upf%lmax
  upf%nqlc = 2* upf%lmax+1

  call radial_grid_COPY(grid, internal_grid)
  !
  upf%grid => internal_grid
  upf%mesh  = upf%grid%mesh
  upf%dx    = upf%grid%dx
  upf%xmin  = upf%grid%xmin
  upf%zmesh = upf%grid%zmesh
  upf%rmax  = upf%grid%rmax
  !
  upf%r   => upf%grid%r
  upf%rab => upf%grid%rab
  !
  ! when possible, write semilocal PP's in the UPF file - may be
  ! useful if one wants to use PPs in the UPF format in other codes
  !
  if( pseudotype == 1 ) then
      if ( rel == 2 ) then
        allocate(upf%vnl(1:grid%mesh, 0:upf%lmax,2))
     else
        allocate(upf%vnl(1:grid%mesh, 0:upf%lmax,1))
     end if
     do nb=1, nbeta
        l=lls(nb)
        if ( rel < 2 .or. l == 0 .or. &
             abs(jjs(nb)-l+0.5_dp) < 0.001_dp) then
           ind = 1
        else if ( rel == 2 .and. l > 0 .and. &
                  abs(jjs(nb)-l-0.5_dp) < 0.001_dp) then
           ind = 2
        endif
        upf%vnl(1:grid%mesh,l,ind) = vnl(1:grid%mesh,l,ind) + &
                                     vpsloc(1:grid%mesh)
     end do
  end if
  !
  allocate(upf%lll(nbeta))
  upf%lll(1:nbeta) = lls(1:nbeta)
  !
  ! *initial* wavefunctions indexes and parameters
  allocate(upf%els(upf%nwfc), upf%oc(upf%nwfc), &
           upf%nchi(upf%nwfc), upf%lchi(upf%nwfc), &
           upf%epseu(upf%nwfc), upf%rcut_chi(upf%nwfc), &
           upf%rcutus_chi(upf%nwfc) )
  upf%els(1:upf%nwfc)   = elts(1:upf%nwfc)
  upf%oc(1:upf%nwfc)    = octs(1:upf%nwfc)
  upf%lchi(1:upf%nwfc)  = llts(1:upf%nwfc)
  upf%nchi(1:upf%nwfc)  = nnts(1:upf%nwfc)
  upf%epseu(1:upf%nwfc) = enlts(1:upf%nwfc)
  upf%rcut_chi(1:upf%nwfc)   = rcutts(1:upf%nwfc)
  upf%rcutus_chi(1:upf%nwfc) = rcutusts(1:upf%nwfc)
  !
  ! projectors indexes and parameters
  !
  allocate(upf%kbeta(nbeta), upf%els_beta(nbeta),&
           upf%rcut(nbeta), upf%rcutus(nbeta))
  do nb=1,nbeta
     upf%kbeta(nb)   = ikk(nb)
     upf%els_beta(nb)= els(nb)
     upf%rcut(nb)    = rcut(nb)
     upf%rcutus(nb)  = rcutus(nb)
  end do
  upf%kkbeta = maxval(upf%kbeta(1:nbeta))
  !
  ! Save GENERATION configuration: not needed to use the pseudopotential, 
  ! but must be saved for reference and for re-generating the pseudo
  !
   at_conf%nwfs  = nwfs
   if (tm) then
      at_conf%pseud = 'troullier-martins'
   else
      at_conf%pseud = 'rrkj'
   endif

   allocate(at_conf%els   (nwfs),&
            at_conf%nns   (nwfs),&
            at_conf%lls   (nwfs),&
            at_conf%ocs   (nwfs),&
            at_conf%rcut  (nwfs),&
            at_conf%rcutus(nwfs),&
            at_conf%enls  (nwfs))
   at_conf%els   (1:nwfs) = els   (1:nwfs) ! label (char*2)
   at_conf%nns   (1:nwfs) = nns   (1:nwfs) ! n
   at_conf%lls   (1:nwfs) = lls   (1:nwfs) ! l
   at_conf%ocs   (1:nwfs) = ocs   (1:nwfs) ! occupation
   at_conf%rcut  (1:nwfs) = rcut  (1:nwfs) ! inner cutoff radius
   at_conf%rcutus(1:nwfs) = rcutus(1:nwfs) ! outer cutoff radius
   at_conf%enls  (1:nwfs) = enls  (1:nwfs) ! one-particle energy


  ! projectors
  allocate(upf%beta(grid%mesh, upf%nbeta))
  upf%beta(1:grid%mesh, 1:upf%nbeta) = betas(1:grid%mesh, 1:nbeta)
  !
  ! hamiltonian terms
  allocate(upf%dion(upf%nbeta, upf%nbeta))
  upf%dion(1:upf%nbeta, 1:upf%nbeta) = bmat(1:nbeta, 1:nbeta)
  !
  if (pseudotype.eq.3) then
     allocate(upf%qqq(upf%nbeta, upf%nbeta))
     upf%qqq(1:upf%nbeta,1:upf%nbeta) = qq(1:nbeta,1:nbeta)
     !
     upf%qqq_eps = 1.e-12_dp ! (hardcoded)
     upf%nqf = 0             ! polinomial expansion of aug.charge is not supported by atomic
     !
     if (upf%q_with_l .or. lpaw) then
        allocate(upf%qfuncl(upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:2*upf%lmax))
     else
        allocate(upf%qfunc(upf%mesh, upf%nbeta*(upf%nbeta+1)/2))
     endif
     !
     if(lpaw) qvanl(1:grid%mesh,:,:,:) = pawsetup%augfun(1:grid%mesh,:,:,:)
     do ibeta=1,nbeta
        do jbeta=ibeta,nbeta
           kbeta = jbeta * (jbeta-1) / 2 + ibeta
           if (upf%q_with_l .or. lpaw) then
              l1=upf%lll(ibeta)
              l2=upf%lll(jbeta)
              do l=abs(l1-l2), l1+l2
                 upf%qfuncl(1:grid%mesh,kbeta,l) = qvanl(1:grid%mesh,ibeta,jbeta,l)
              enddo
           else
              upf%qfunc(1:grid%mesh,kbeta) = qvan (1:grid%mesh, ibeta, jbeta)
           endif
        enddo
     enddo
     !
  endif
  !
  allocate(upf%rho_atc(upf%mesh))
  if (upf%nlcc) then
     upf%rho_atc(1:grid%mesh) = rhoc(1:grid%mesh)/fpi/grid%r2(1:grid%mesh)
  else
     upf%rho_atc(:) = 0.0_dp
  end if

  allocate(upf%rho_at(upf%mesh))
  upf%rho_at (1:grid%mesh) = rhos (1:grid%mesh,1)
  !
  allocate(upf%chi(upf%mesh,upf%nwfc))
  upf%chi(1:grid%mesh,1:upf%nwfc) = phits(1:grid%mesh,1:upf%nwfc)
  !
  allocate(upf%vloc(upf%mesh))
  upf%vloc(1:grid%mesh) = vpsloc(1:grid%mesh)
  upf%lloc = lloc
  upf%rcloc = rcloc
  !
  !
  if (upf%has_so)    CALL export_upf_so()
  if (upf%tpawp)     CALL export_upf_paw()
  if (upf%has_gipaw) CALL export_upf_gipaw()
  upf%has_wfc = lsave_wfc
  if (upf%has_wfc)   CALL export_upf_wfc()
  !
  CALL write_upf_v2( iunps, upf, at_conf, unit_loc )
  !
  CALL deallocate_pseudo_upf( upf )
  CALL deallocate_radial_grid( internal_grid )
  CALL deallocate_pseudo_config( at_conf )

   RETURN

 CONTAINS

   SUBROUTINE export_upf_wfc
      ALLOCATE( upf%aewfc(upf%mesh, upf%nbeta), upf%pswfc(upf%mesh, upf%nbeta) )
      upf%aewfc(1:upf%mesh,1:upf%nbeta) = psipaw(1:upf%mesh,1:upf%nbeta)
      upf%pswfc(1:upf%mesh,1:upf%nbeta) = phis(1:upf%mesh,1:upf%nbeta)
   END SUBROUTINE export_upf_wfc

   SUBROUTINE export_upf_so
      ALLOCATE( upf%nn(upf%nwfc), upf%jchi(upf%nwfc), upf%jjj(upf%nbeta) )

      upf%els(1:upf%nwfc)  = elts(1:upf%nwfc)
      upf%nn(1:upf%nwfc)   = nnts(1:upf%nwfc)
      upf%lchi(1:upf%nwfc) = llts(1:upf%nwfc)
      upf%jchi(1:upf%nwfc) = jjts(1:upf%nwfc)
      !
      upf%lll(1:upf%nbeta) = lls(1:upf%nbeta)
      upf%jjj(1:upf%nbeta) = jjs(1:upf%nbeta)

   END SUBROUTINE export_upf_so
   !
   SUBROUTINE export_upf_paw
      INTEGER :: co,n   !EMINE
      upf%paw_data_format = 2
      !
      upf%paw%core_energy = etot -etots
      upf%paw%lmax_aug = 2*upf%lmax
      upf%paw%augshape = which_augfun
      upf%paw%raug     = rmatch_augfun
      upf%paw%iraug    = pawsetup%irc

      allocate(upf%paw%ae_rho_atc(upf%mesh))
      upf%paw%ae_rho_atc(1:upf%mesh) = pawsetup%aeccharge(1:upf%mesh)/fpi/grid%r2(1:grid%mesh)
      !
      allocate(upf%paw%ae_vloc(upf%mesh))
      upf%paw%ae_vloc(1:upf%mesh)    = pawsetup%aeloc(1:upf%mesh)
      !
      allocate(upf%paw%oc(upf%nbeta))
      do nb = 1,upf%nbeta
         upf%paw%oc(nb)  = max(pawsetup%oc(nb),0._dp)
      enddo
      !
      allocate(upf%paw%augmom(upf%nbeta, upf%nbeta, 0:2*upf%lmax))
      upf%paw%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax) &
            = pawsetup%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax)
      !
      upf%kkbeta = max(upf%kkbeta, upf%paw%iraug)

      IF (upf%has_so) THEN
         ALLOCATE( upf%paw%aewfc_rel(upf%mesh, upf%nbeta) )
         upf%paw%aewfc_rel(1:upf%mesh,1:upf%nbeta) = &
                              pawsetup%aewfc_rel(1:upf%mesh,1:upf%nbeta)
      ENDIF
      !
      !upf%paw%pfunc(:)  = not used when writing, reconstructed from upf%aewfc
      !upf%paw%ptfunc(:) = not used when writing, reconstructed from upf%pswfc
      !===============================================================
      !For PAW pseudopotentials, now we also include core information:
      !even when lgipaw_reconstruction = .false.
      !EMINE
      upf%gipaw_ncore_orbitals = COUNT(core_state(1:nwf))
      co = upf%gipaw_ncore_orbitals
      ALLOCATE ( &
         upf%gipaw_core_orbital_n(co), &
         upf%gipaw_core_orbital_l(co), &
         upf%gipaw_core_orbital_el(co), &
         upf%gipaw_core_orbital(upf%mesh,co))
      upf%gipaw_core_orbital_n(1:co)  = nn(1:co)
      upf%gipaw_core_orbital_l(1:co)  = ll(1:co)
      upf%gipaw_core_orbital_el(1:co) = el(1:co)
      DO n = 1,co
         upf%gipaw_core_orbital(1:upf%mesh,n) &
            = psi(1:upf%mesh,1,n)
      ENDDO
      !================================================================
      RETURN
   END SUBROUTINE export_upf_paw
   SUBROUTINE export_upf_gipaw
      INTEGER :: co,nw,n,nb

      IF ( nconf /= 1 ) CALL errore ( "write_gipaw_orbitals", &
            "The (GI)PAW reconstruction requires one test configuration", abs(nconf) )

      upf%gipaw_data_format = 2  ! The version of the format

      upf%gipaw_ncore_orbitals = COUNT(core_state(1:nwf))
      co = upf%gipaw_ncore_orbitals
      upf%gipaw_wfs_nchannels  = nwfts
      nw = upf%gipaw_wfs_nchannels

      ALLOCATE ( &
         upf%gipaw_core_orbital_n(co), &
         upf%gipaw_core_orbital_l(co), &
         upf%gipaw_core_orbital_el(co), &
         upf%gipaw_core_orbital(upf%mesh,co), &
         upf%gipaw_wfs_el(nw), &
         upf%gipaw_wfs_ll(nw), &
         upf%gipaw_wfs_rcut(nw), &
         upf%gipaw_wfs_rcutus(nw), &
         upf%gipaw_wfs_ae(upf%mesh,nw), &
         upf%gipaw_wfs_ps(upf%mesh,nw), &
         upf%gipaw_vlocal_ae(upf%mesh), &
         upf%gipaw_vlocal_ps(upf%mesh) &
         )

      upf%gipaw_core_orbital_n(1:co)  = nn(1:co)
      upf%gipaw_core_orbital_l(1:co)  = ll(1:co)
      upf%gipaw_core_orbital_el(1:co) = el(1:co)

      DO n = 1,co
         upf%gipaw_core_orbital(1:upf%mesh,n) &
            = psi(1:upf%mesh,1,n)
      ENDDO

      upf%gipaw_vlocal_ae(1:upf%mesh) &
            = grid%r(1:upf%mesh) * vpot(1:mesh,1)
      upf%gipaw_vlocal_ps(1:upf%mesh) &
            = grid%r(1:upf%mesh) * vpstot(1:mesh,1)

      upf%gipaw_wfs_el(1:nw)     = elts(1:nw)
      upf%gipaw_wfs_ll(1:nw)     = lltsc(1:nw,1)

     ! Find the suitable core radii to be written out
     !*apsi WARNING: DOES NOT WORK WITH VANDERBILT PP YET
      DO nb = 1,nw
         upf%gipaw_wfs_rcut(nb) = -1.0_dp
         DO n = 1, nwfs
            IF ( els(n) == eltsc(nb,1) ) THEN
               upf%gipaw_wfs_rcut(nb)   = rcuttsc(nb,1)
               upf%gipaw_wfs_rcutus(nb) = rcutustsc(nb,1)
            END IF
         END DO
         !
         IF ( upf%gipaw_wfs_rcut(nb) < 0.0_dp ) THEN
            DO n = 1, nwfs
               ! If there is one with the same l...
               IF ( lltsc(nb,1) == lls(n) ) THEN
                  upf%gipaw_wfs_rcut(nb)   = rcuttsc(nb,1)
                  upf%gipaw_wfs_rcutus(nb) = rcutustsc(nb,1)
               END IF
            END DO
         END IF
      ENDDO

      DO n = 1,nw
         !
         upf%gipaw_wfs_ae(1:upf%mesh,n) = wfc_ae_recon(1:upf%mesh,nstoaets(n))
         !
         upf%gipaw_wfs_ps(1:upf%mesh,n) = wfc_ps_recon(1:upf%mesh,n)
      ENDDO


      RETURN
   END SUBROUTINE export_upf_gipaw
   !
END SUBROUTINE export_upf
