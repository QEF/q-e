!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM projwfc
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! calculates Lowdin charges, spilling parameter, projected DOS
  ! or computes the LDOS in a volume given in input as function of energy
  !
  !
  ! Input (namelist &inputpp ... / ):                         Default value
  !
  !    prefix        prefix of input file produced by pw.x    'pwscf'
  !                    (wavefunctions are needed)
  !    outdir        directory containing the input file       ./
  !    ngauss        type of gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            =-99  Fermi-Dirac function
  !    degauss       gaussian broadening, Ry (not eV!)          0.0
  !    Emin, Emax    min, max energy (eV) for DOS plot          band extrema
  !    DeltaE        energy grid step (eV)                      none
  !    lsym          if true the projections are symmetrized    .true.
  !    filproj       file containing the projections            none
  !    filpdos       prefix for output files containing PDOS(E) prefix
  !    lgww          if .true. take energies from previous GWW calculation
  !                  (file bands.dat)
  !    kresolveddos  if .true. the DOS is written as function   .false.
  !                  of the k-point (not summed over all of them)
  !                  all k-points but results
  !    tdosinboxes   if .true., the local DOS in specified      .false.
  !                  volumes (boxes) is computed
  !    n_proj_boxes  number of volumes for the local DOS        0
  !    irmin, irmax  first and last point of the FFT grid       1, 0
  !                  included in the volume
  !    plotboxes     if .true., the volumes are given in output .false.
  !
  !
  ! Output:
  !
  !   Projections are written to standard output,
  !   and also to file filproj if given as input.
  !   The total DOS and the sum of projected DOS are written to file
  !   "filpdos".pdos_tot.
  !   The format for the collinear, spin-unpolarized case and the
  !   non-collinear, spin-orbit case is
  !        E DOS(E) PDOS(E)
  !   The format for the collinear, spin-polarized case is
  !        E DOSup(E) DOSdw(E)  PDOSup(E) PDOSdw(E)
  !   The format for the non-collinear, non spin-orbit case is
  !        E DOS(E) PDOSup(E) PDOSdw(E)
  !
  !   In the collinear case and the non-collinear, non spin-orbit case
  !   projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   (one file per atomic wavefunction found in the pseudopotential file)
  !   - The format for the collinear, spin-unpolarized case is
  !        E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
  !     where LDOS = \sum m=1,2l+1 PDOS_m(E)
  !     and PDOS_m(E) = projected DOS on atomic wfc with component m
  !   - The format for the collinear, spin-polarized case and the
  !     non-collinear, non spin-orbit case is as above with
  !     two components for both  LDOS(E) and PDOS_m(E)
  !
  !   In the non-collinear, spin-orbit case (i.e. if there is at least one
  !   fully relativistic pseudopotential) wavefunctions are projected
  !   onto eigenstates of the total angular-momentum.
  !   Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l_j),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   and j is the value of the total angular momentum.
  !   In this case the format is
  !      E LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
  !
  !   All DOS(E) are in states/eV plotted vs E in eV
  !
  !   If the kresolveddos option is used, the k-point index is prepended
  !   to the formats above, e.g. (collinear, spin-unpolarized case)
  !        ik E DOS(E) PDOS(E)
  !
  !   If the local DOS(E) is computed (tdosinboxes),
  !   projections are written to file "filproj" if given as input.
  !   Volumes are written as xsf files with 3D datagrids, valued 1.0
  !   inside the box volume and 0 outside, named filpdos.box#n.xsf
  !   The local DOS(E) is written to file "filpdos".ldos_boxes, with format
  !      E totDOS(E) SumLDOS(E) LDOS_1(E) ... LDOS_n(E)
  !
  !  Order of m-components for each l in the output:
  !
  !  1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)
  !
  !  where phi is the polar angle:x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
  !  This is determined in file flib/ylmr2.f90 that calculates spherical harm.
  !      L=1 :
  !  1 pz     (m=0)
  !  2 px     (real combination of m=+/-1 with cosine)
  !  3 py     (real combination of m=+/-1 with sine)
  !      L=2 :
  !  1 dz2    (m=0)
  !  2 dzx    (real combination of m=+/-1 with cosine)
  !  3 dzy    (real combination of m=+/-1 with sine)
  !  4 dx2-y2 (real combination of m=+/-2 with cosine)
  !  5 dxy    (real combination of m=+/-1 with sine)
  !
  ! Important notice:
  !
  !    The tetrahedron method is presently not implemented.
  !    Gaussian broadening is used in all cases:
  !    - if degauss is set to some value in namelist &inputpp, that value
  !      (and the optional value for ngauss) is used
  !    - if degauss is NOT set to any value in namelist &inputpp, the
  !      value of degauss and of ngauss are read from the input data
  !      file (they will be the same used in the pw.x calculations)
  !    - if degauss is NOT set to any value in namelist &inputpp, AND
  !      there is no value of degauss and of ngauss in the input data
  !      file, degauss=DeltaE (in Ry) and ngauss=0 will be used
  ! Obsolete variables, ignored:
  !   io_choice
  !   smoothing
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : degauss, ngauss, lgauss
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup, nproc_ortho
  USE environment,      ONLY : environment_start
  !
  ! for GWW
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filpdos, filproj, io_choice, outdir
  REAL (DP)      :: Emin, Emax, DeltaE, degauss1, smoothing
  INTEGER :: ngauss1, ios
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)

  !
  ! for GWW
  INTEGER :: iun, idum
  REAL(DP) :: rdum1,rdum2,rdum3
  LOGICAL :: lex, lgww
  !
  !
  NAMELIST / inputpp / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, DeltaE, io_choice, smoothing, filpdos, filproj, &
             lgww, & !if .true. use GW QP energies from file bands.dat
             kresolveddos, tdosinboxes, n_proj_boxes, irmin, irmax, plotboxes
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PROJWFC' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.01d0
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  lgww   = .false.
  kresolveddos = .false.
  tdosinboxes = .false.
  plotboxes   = .false.
  n_proj_boxes= 1
  irmin(:,:)  = 1
  irmax(:,:)  = 0
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id )
  !
  IF (ios /= 0) CALL errore ('projwfc', 'reading inputpp namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix,  ionode_id )
  CALL mp_bcast( filproj,  ionode_id )
  CALL mp_bcast( ngauss1, ionode_id )
  CALL mp_bcast( degauss1,ionode_id )
  CALL mp_bcast( DeltaE,  ionode_id )
  CALL mp_bcast( lsym,  ionode_id )
  CALL mp_bcast( Emin, ionode_id )
  CALL mp_bcast( Emax, ionode_id )
  ! for GWW
  CALL mp_bcast( lgww, ionode_id )
  ! for projection on boxes
  CALL mp_bcast( tdosinboxes, ionode_id )
  CALL mp_bcast( n_proj_boxes, ionode_id )
  CALL mp_bcast( irmin, ionode_id )
  CALL mp_bcast( irmax, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  CALL openfil_pp ( )
  !
  !   decide Gaussian broadening
  !
  IF (degauss1/=0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ELSEIF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  ELSE
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !
  IF ( tdosinboxes ) THEN
     IF( nproc_ortho > 1 ) THEN
        CALL errore ('projwfc', 'nproc_ortho > 1 not yet implemented', 1)
     ELSE
        CALL projwave_boxes (filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes)
     ENDIF
  ELSE
     IF (noncolin) THEN
        CALL projwave_nc(filproj, lsym )
     ELSE
        IF( nproc_ortho > 1 ) THEN
           CALL pprojwave (filproj, lsym)
        ELSE
           CALL projwave (filproj, lsym, lgww)
        ENDIF
     ENDIF
  ENDIF
  !
  IF ( ionode ) THEN
     IF ( tdosinboxes ) THEN
        CALL partialdos_boxes (Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
     ELSE
        IF ( lsym ) THEN
           !
           IF (noncolin) THEN
              CALL partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
           ELSE
              CALL partialdos (Emin, Emax, DeltaE, kresolveddos, filpdos)
           ENDIF
           !
        ENDIF
     ENDIF
  ENDIF
  !
  CALL stop_pp
  !
END PROGRAM projwfc
!
MODULE projections
  USE kinds, ONLY : DP

  TYPE wfc_label
     INTEGER na, n, l, m
  END TYPE wfc_label
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)

END MODULE projections
!
MODULE projections_nc
  USE kinds, ONLY : DP

  TYPE wfc_label_nc
     INTEGER na, n, l, m, ind
     REAL (DP) jj
  END TYPE wfc_label_nc
  TYPE(wfc_label_nc), ALLOCATABLE :: nlmchi(:)

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)

END MODULE projections_nc
!
MODULE projections_ldos
  USE kinds, ONLY : DP
  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
END MODULE projections_ldos
!
!-----------------------------------------------------------------------
SUBROUTINE projwave( filproj, lsym, lgww )
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE fft_base, ONLY : dfftp
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symm_base, ONLY: nsym, irt, d1, d2, d3
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE spin_orb, ONLY: lspinorb
  USE wavefunctions_module, ONLY: evc
  !
  USE projections
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  CHARACTER (len=*) :: filproj
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, nwfc,&
       nwfc1, lmax_wfc, is, ios, iunproj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE ::roverlap(:,:), rwork1(:),rproj0(:,:)
  ! ... or for gamma-point.
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2)
  INTEGER  :: nksinit, nkslast
  CHARACTER(len=256) :: filename
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  INTEGER, ALLOCATABLE :: idx(:)
  LOGICAL :: lsym
  !
  !
  ! for GWW
  INTEGER :: iun, idum
  REAL(DP) :: rdum1,rdum2,rdum3
  LOGICAL :: lex, lgww
  !
  !
  WRITE( stdout, '(/5x,"Calling projwave .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  !
  ! fill structure nlmchi
  !
  ALLOCATE (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           lmax_wfc = max (lmax_wfc, l )
           DO m = 1, 2 * l + 1
              nwfc=nwfc+1
              nlmchi(nwfc)%na = na
              nlmchi(nwfc)%n  =  n
              nlmchi(nwfc)%l  =  l
              nlmchi(nwfc)%m  =  m
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  IF (lmax_wfc > 3) CALL errore ('projwave', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projwave', 'wrong # of atomic wfcs?', 1)
  !
  ALLOCATE( proj (natomwfc, nbnd, nkstot) )
  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj      = 0.d0
  proj_aux  = (0.d0, 0.d0)
  !
  IF (.not. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) )
  ALLOCATE(wfcatom (npwx, natomwfc) )
  ALLOCATE(overlap (natomwfc, natomwfc) )
  overlap= (0.d0,0.d0)
  !
  IF ( gamma_only ) THEN
     ALLOCATE(roverlap (natomwfc, natomwfc) )
     roverlap= 0.d0
  ENDIF
  CALL allocate_bec_type (nkb, natomwfc, becp )
  ALLOCATE(e (natomwfc) )
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  DO ik = 1, nks
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)

     CALL atomic_wfc (ik, wfcatom)

     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL calbec ( npw, vkb, wfcatom, becp)

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     IF ( gamma_only ) THEN
        CALL calbec ( npw, wfcatom, swfcatom, roverlap )
        overlap(:,:)=cmplx(roverlap(:,:),0.0_dp, kind=dp)
        ! TEMP: diagonalization routine for real matrix should be used instead
     ELSE
        CALL calbec ( npw, wfcatom, swfcatom, overlap )
     ENDIF

     !
     ! calculate O^{-1/2}
     !
     ALLOCATE(work (natomwfc, natomwfc) )
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO
     DO i = 1, natomwfc
        DO j = i, natomwfc
           overlap (i, j) = (0.d0, 0.d0)
           DO k = 1, natomwfc
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * conjg (work (i, k) )
           ENDDO
           IF (j /= i) overlap (j, i) = conjg (overlap (i, j))
        ENDDO
     ENDDO
     DEALLOCATE (work)
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     IF ( gamma_only ) THEN
        roverlap(:,:)=REAL(overlap(:,:),DP)
        ! TEMP: diagonalization routine for real matrix should be used instead
        CALL DGEMM ('n', 't', 2*npw, natomwfc, natomwfc, 1.d0 , &
             swfcatom, 2*npwx,  roverlap, natomwfc, 0.d0, wfcatom, 2*npwx)
     ELSE
        CALL ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , &
             swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx)
     ENDIF

     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>
     !
     IF ( gamma_only ) THEN
        !
        ALLOCATE(rproj0(natomwfc,nbnd), rwork1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, rproj0)
        !
        proj_aux(:,:,ik) = cmplx( rproj0(:,:), 0.0_dp, kind=dp )
        !
     ELSE
        !
        ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, proj0)
        !
        proj_aux(:,:,ik) = proj0(:,:)
        !
     ENDIF
     !
     ! symmetrize the projections
     !
     IF (lsym) THEN
        DO nwfc = 1, natomwfc
           !
           !  atomic wavefunction nwfc is on atom na
           !
           na= nlmchi(nwfc)%na
           n = nlmchi(nwfc)%n
           l = nlmchi(nwfc)%l
           m = nlmchi(nwfc)%m
           !
           DO isym = 1, nsym
              nb = irt (isym, na)
              DO nwfc1 =1, natomwfc
                 IF (nlmchi(nwfc1)%na == nb             .and. &
                      nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                      nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                      nlmchi(nwfc1)%m == 1 ) GOTO 10
              ENDDO
              CALL errore('projwave','cannot symmetrize',1)
10            nwfc1=nwfc1-1
              !
              !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
              !
              IF ( gamma_only ) THEN
                 IF (l == 0) THEN
                    rwork1(:) = rproj0 (nwfc1 + 1,:)
                 ELSEIF (l == 1) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 3
                       rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 2) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 5
                       rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 3) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 7
                       rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         rwork1(ibnd) * rwork1(ibnd) / nsym
                 ENDDO
              ELSE
                 IF (l == 0) THEN
                    work1(:) = proj0 (nwfc1 + 1,:)
                 ELSEIF (l == 1) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 3
                       work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 2) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 5
                       work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 3) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 7
                       work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         work1(ibnd) * conjg (work1(ibnd)) / nsym
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ELSE
        IF ( gamma_only ) THEN
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(rproj0(nwfc,ibnd))**2
              ENDDO
           ENDDO
        ELSE
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(proj0(nwfc,ibnd))**2
              ENDDO
           ENDDO
        ENDIF
     ENDIF
     IF ( gamma_only ) THEN
        DEALLOCATE (rwork1)
        DEALLOCATE (rproj0)
     ELSE
        DEALLOCATE (work1)
        DEALLOCATE (proj0)
     ENDIF
     ! on k-points
  ENDDO
  !
  DEALLOCATE (e)
  IF ( gamma_only ) THEN
     DEALLOCATE (roverlap)
  ENDIF
  CALL deallocate_bec_type (becp)
  DEALLOCATE (overlap)
  DEALLOCATE (wfcatom)
  IF (.not. lda_plus_u) DEALLOCATE (swfcatom)
  !
  !   vectors et and proj are distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (et,       nbnd, nkstot, nks)
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks)
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  !
  !!!! for GWW
  IF(lgww) THEN
    INQUIRE ( file='bands.dat', EXIST=lex )
    WRITE(stdout,*) 'lex=', lex
    CALL flush_unit(stdout)
    !
    IF(lex) THEN
       WRITE(stdout,*) 'Read the file bands.dat => GWA Eigenvalues used.'
       CALL flush_unit(stdout)
       iun = find_free_unit()
       OPEN(unit=iun, file='bands.dat', status='unknown', form='formatted', IOSTAT=ios)
       READ(iun,*) idum
       DO i=1, nbnd
         READ(iun,*) idum,rdum1,rdum2,et(i,1),rdum3
       ENDDO
       et(:,1)=et(:,1)/rytoev !! because in bands.dat file, the QP energies are in eV
    ELSE
       WRITE(stdout,*) 'The file bands.dat does not exist.'
       WRITE(stdout,*) 'Eigenergies are not modified'
       CALL flush_unit(stdout)
    ENDIF
    !!!! end GWW
    !
  ENDIF
  IF ( ionode ) THEN
     !
     ! write on the file filproj
     !
     IF (filproj/=' ') THEN
        DO is=1,nspin
           IF (nspin==2) THEN
              IF (is==1) filename=trim(filproj)//'.up'
              IF (is==2) filename=trim(filproj)//'.down'
              nksinit=(nkstot/2)*(is-1)+1
              nkslast=(nkstot/2)*is
           ELSE
              filename=trim(filproj)
              nksinit=1
              nkslast=nkstot
           ENDIF
           iunproj=33
           CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual,   &
                ecutwfc, nkstot/nspin, nbnd, natomwfc)
           DO nwfc = 1, natomwfc
              WRITE(iunproj,'(2i5,a3,3i5)') &
                  nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                  nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
              DO ik=nksinit,nkslast
                 DO ibnd=1,nbnd
                   WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, &
                                                 abs(proj(nwfc,ibnd,ik))
                 ENDDO
              ENDDO
           ENDDO
           CLOSE(iunproj)
        ENDDO
     ENDIF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     !
     ! write on the standard output file
     !
     WRITE( stdout,'(/5x,"Atomic states used for projection")')
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
     DO nwfc = 1, natomwfc
        WRITE(stdout,1000) &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
     ENDDO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")")
     !
     ALLOCATE(idx(natomwfc), proj1 (natomwfc) )
     DO ik = 1, nkstot
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
        DO ibnd = 1, nbnd
           WRITE( stdout, '("==== e(",i4,") = ",f11.5," eV ==== ")') &
              ibnd, et (ibnd, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           DO nwfc = 1, natomwfc
              idx (nwfc) = 0
              proj1 (nwfc) = - proj (nwfc, ibnd, ik)
           ENDDO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (natomwfc, proj1, idx, eps4)
           !
           !  only projections that are larger than 0.001 are written
           !
           DO nwfc = 1, natomwfc
              proj1 (nwfc) = - proj1(nwfc)
              IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
           ENDDO
           nwfc = natomwfc + 1
20         nwfc = nwfc -1
           !
           ! fancy (?!?) formatting
           !
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') &
                (proj1 (i), idx(i), i = 1, min(5,nwfc))
           DO j = 1, (nwfc-1)/5
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') &
                   (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
           ENDDO
           psum = 0.d0
           DO nwfc = 1, natomwfc
              psum = psum + proj (nwfc, ibnd, ik)
           ENDDO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        ENDDO
     ENDDO
     DEALLOCATE (idx, proj1)
     !
     ! estimate partial charges (Loewdin) on each atom
     !
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin ) )
     charges = 0.0d0
     DO ik = 1, nkstot
        IF ( nspin == 1 ) THEN
           current_spin = 1
        ELSEIF ( nspin == 2 ) THEN
           current_spin = isk ( ik )
        ELSE
           CALL errore ('projwfc_nc',' called in the wrong case ',1)
        ENDIF
        DO ibnd = 1, nbnd
           DO nwfc = 1, natomwfc
              na= nlmchi(nwfc)%na
              l = nlmchi(nwfc)%l
              charges(na,l,current_spin) = charges(na,l,current_spin) + &
                   wg (ibnd,ik) * proj (nwfc, ibnd, ik)
           ENDDO
        ENDDO
     ENDDO
     !
     WRITE( stdout, '(/"Lowdin Charges: "/)')
     !
     DO na = 1, nat
        DO is = 1, nspin
           totcharge(is) = sum(charges(na,0:lmax_wfc,is))
        ENDDO
        IF ( nspin == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
        ELSEIF ( nspin == 2) THEN
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
           WRITE( stdout, 2002) totcharge(2), &
                ( l_label(l), charges(na,l,2), l= 0,lmax_wfc)
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        ENDIF
     ENDDO
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4,4(", ",a1," =",f8.4))
2001 FORMAT (15x,"  spin up      = ",f8.4,4(", ",a1," =",f8.4))
2002 FORMAT (15x,"  spin down    = ",f8.4,4(", ",a1," =",f8.4))
2003 FORMAT (15x,"  polarization = ",f8.4,4(", ",a1," =",f8.4))
     !
     psum = sum(charges(:,:,:)) / nelec
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     DEALLOCATE (charges)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE projwave
!
!-----------------------------------------------------------------------
SUBROUTINE projwave_nc(filproj, lsym )
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE printout_base, ONLY: title
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE fft_base, ONLY : dfftp
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU
  USE lsda_mod, ONLY: nspin
  USE noncollin_module, ONLY: noncolin, npol, angle1, angle2
  USE symm_base, ONLY: nsym, irt, t_rev
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  USE spin_orb,   ONLY: lspinorb, domag
  USE projections_nc
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: filproj
  LOGICAL :: lsym
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, ind, n, m, m1, n1, &
             n2, l, nwfc, nwfc1, lmax_wfc, is, nspin0, iunproj,    &
             ind0
  REAL(DP) :: jj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2), fact(2), spinor, compute_mj
  INTEGER, ALLOCATABLE :: idx(:)
  !
  COMPLEX(DP) :: d12(2, 2, 48), d32(4, 4, 48), d52(6, 6, 48), &
                      d72(8, 8, 48)
  COMPLEX(DP) :: d012(2, 2, 48), d112(6, 6, 48), d212(10, 10, 48), &
                      d312(14, 14, 48)
  !
  !
  !
  IF (.not.noncolin) CALL errore('projwave_nc','called in the wrong case',1)
  IF (gamma_only) CALL errore('projwave_nc','gamma_only not yet implemented',1)
  WRITE( stdout, '(/5x,"Calling projwave_nc .... ")')
  !
  ! fill structure nlmchi
  !
  ALLOCATE (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     n2 = 0
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           lmax_wfc = max (lmax_wfc, l )
           IF (lspinorb) THEN
              IF (upf(nt)%has_so) THEN
                jj = upf(nt)%jchi (n)
                ind = 0
                DO m = -l-1, l
                   fact(1) = spinor(l,jj,m,1)
                   fact(2) = spinor(l,jj,m,2)
                   IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                      nwfc = nwfc + 1
                      ind = ind + 1
                      nlmchi(nwfc)%na = na
                      nlmchi(nwfc)%n  =  n
                      nlmchi(nwfc)%l  =  l
                      nlmchi(nwfc)%m  =  m
                      nlmchi(nwfc)%ind  =  ind
                      nlmchi(nwfc)%jj  =  jj
                   ENDIF
                ENDDO
              ELSE
                DO n1 = l, l+1
                  jj= dble(n1) - 0.5d0
                  ind = 0
                  IF (jj>0.d0)  THEN
                    n2 = n2 + 1
                    DO m = -l-1, l
                      fact(1) = spinor(l,jj,m,1)
                      fact(2) = spinor(l,jj,m,2)
                      IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                        nwfc = nwfc + 1
                        ind = ind + 1
                        nlmchi(nwfc)%na = na
                        nlmchi(nwfc)%n  =  n2
                        nlmchi(nwfc)%l  =  l
                        nlmchi(nwfc)%m  =  m
                        nlmchi(nwfc)%ind  =  ind
                        nlmchi(nwfc)%jj  =  jj
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
           ELSE
              DO m = 1, 2 * l + 1
                 nwfc=nwfc+1
                 nlmchi(nwfc)%na = na
                 nlmchi(nwfc)%n  =  n
                 nlmchi(nwfc)%l  =  l
                 nlmchi(nwfc)%m  =  m
                 nlmchi(nwfc)%ind  =  m
                 nlmchi(nwfc)%jj  =  0.d0
                 nlmchi(nwfc+2*l+1)%na = na
                 nlmchi(nwfc+2*l+1)%n  =  n
                 nlmchi(nwfc+2*l+1)%l  =  l
                 nlmchi(nwfc+2*l+1)%m  =  m
                 nlmchi(nwfc+2*l+1)%ind  =  m+2*l+1
                 nlmchi(nwfc+2*l+1)%jj  =  0.d0
              ENDDO
              nwfc=nwfc+2*l+1
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !
  IF (lmax_wfc > 3) CALL errore ('projwave_nc', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projwave_nc','wrong # of atomic wfcs?',1)
  !
  ALLOCATE(wfcatom (npwx*npol,natomwfc) )
  IF (.not. lda_plus_u) ALLOCATE(swfcatom (npwx*npol, natomwfc ) )
  CALL allocate_bec_type (nkb, natomwfc, becp )
  ALLOCATE(e (natomwfc) )
  ALLOCATE(work (natomwfc, natomwfc) )
  !
  ALLOCATE(overlap (natomwfc, natomwfc) )
  ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) )
  ALLOCATE(proj (natomwfc, nbnd, nkstot) )
  ALLOCATE(proj_aux (natomwfc, nbnd, nkstot) )
  overlap  = (0.d0,0.d0)
  proj0    = (0.d0,0.d0)
  proj     = 0.d0
  proj_aux = (0.d0,0.d0)
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  IF (lspinorb) THEN
     !
     ! initialize D_Sj for j=1/2, j=3/2, j=5/2 and j=7/2
     !
     CALL d_matrix_so (d12, d32, d52, d72)
     !
  ELSE
     !
     ! initialize D_Sl for l=0, l=1, l=2 and l=3
     !
     CALL d_matrix_nc (d012, d112, d212, d312)
     !
  ENDIF
  !
  DO ik = 1, nks
     wfcatom = (0.d0,0.d0)
     swfcatom= (0.d0,0.d0)
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     CALL atomic_wfc_nc_proj (ik, wfcatom)
     !
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL calbec ( npw, vkb, wfcatom, becp )

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     CALL ZGEMM ('C', 'N', natomwfc, natomwfc, npwx*npol, (1.d0, 0.d0), wfcatom, &
       npwx*npol, swfcatom, npwx*npol, (0.d0, 0.d0), overlap, natomwfc)
     CALL mp_sum ( overlap, intra_pool_comm )
     !
     ! calculate O^{-1/2}
     !
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO
     DO i = 1, natomwfc
        DO j = i, natomwfc
          overlap (i, j) = (0.d0, 0.d0)
          DO k = 1, natomwfc
            overlap(i, j) = overlap(i, j) + e(k) * work(j, k) * conjg(work (i, k) )
          ENDDO
          IF (j /= i) overlap (j, i) = conjg (overlap (i, j))
        ENDDO
     ENDDO
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     CALL ZGEMM ('n', 't', npwx*npol, natomwfc, natomwfc, (1.d0, 0.d0) , &
     swfcatom, npwx*npol,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx*npol)
     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>
     !
     CALL ZGEMM ('C','N',natomwfc, nbnd, npwx*npol, (1.d0, 0.d0), wfcatom, &
                 npwx*npol, evc, npwx*npol, (0.d0, 0.d0), proj0, natomwfc)
     CALL mp_sum ( proj0( :, 1:nbnd ), intra_pool_comm )
     !
     proj_aux(:,:,ik) = proj0(:,:)

     !
     IF (lsym) THEN
        DO nwfc = 1, natomwfc
           !
           !  atomic wavefunction nwfc is on atom na
           !
           IF (lspinorb) THEN
              na= nlmchi(nwfc)%na
              n = nlmchi(nwfc)%n
              l = nlmchi(nwfc)%l
              m = nlmchi(nwfc)%m
              ind0 = nlmchi(nwfc)%ind
              jj = nlmchi(nwfc)%jj
              !
              DO isym = 1, nsym
!-- check for the time reversal
                 IF (t_rev(isym) == 1) THEN
                   ind = 2*jj + 2 - ind0
                 ELSE
                   ind = ind0
                 ENDIF
!--
                 nb = irt (isym, na)
                 DO nwfc1 =1, natomwfc
                    IF (nlmchi(nwfc1)%na == nb             .and. &
                         nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                         nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                         nlmchi(nwfc1)%jj == nlmchi(nwfc)%jj .and. &
                         nlmchi(nwfc1)%ind == 1 ) GOTO 10
                 ENDDO
                 CALL errore('projwave_nc','cannot symmetrize',1)
10               nwfc1=nwfc1-1
                 !
                 !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
                 !
                 IF (abs(jj-0.5d0)<1.d-8) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 2
                       work1(:)=work1(:)+d12(m1,ind,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (abs(jj-1.5d0)<1.d-8) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 4
                       work1(:)=work1(:)+d32(m1,ind,isym)*proj0(nwfc1 + m1,:)
                    ENDDO
                 ELSEIF (abs(jj-2.5d0)<1.d-8) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 6
                       work1(:)=work1(:)+d52(m1,ind,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (abs(jj-3.5d0)<1.d-8) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 8
                       work1(:)=work1(:)+d72(m1,ind,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         work1(ibnd) * conjg (work1(ibnd)) / nsym
                 ENDDO
                 ! on symmetries
!--  in a nonmagnetic case - another loop with the time reversal
                 IF (.not.domag.and.ind==ind0) THEN
                   ind = 2*jj + 2 - ind0
                   nwfc1 = nwfc1 + 1
                   GOTO 10
                 ENDIF
!--
              ENDDO
!--  in a nonmagnetic case - rescale
              IF (.not.domag) THEN
                DO ibnd = 1, nbnd
                  proj(nwfc,ibnd,ik) = 0.5d0*proj(nwfc,ibnd,ik)
                ENDDO
              ENDIF
!--
           ELSE
              na= nlmchi(nwfc)%na
              n = nlmchi(nwfc)%n
              l = nlmchi(nwfc)%l
              m = nlmchi(nwfc)%m
              ind0 = nlmchi(nwfc)%ind
              !
              DO isym = 1, nsym
!-- check for the time reversal
                 IF (t_rev(isym) == 1) THEN
                   ind = 2*m - ind0 + 2*l + 1
                 ELSE
                   ind = ind0
                 ENDIF
!--
                 nb = irt (isym, na)
                 DO nwfc1 =1, natomwfc
                    IF (nlmchi(nwfc1)%na == nb             .and. &
                        nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                        nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                        nlmchi(nwfc1)%m == 1 .and. &
                        nlmchi(nwfc1)%ind == 1) GOTO 15
                 ENDDO
                 CALL errore('projwave_nc','cannot symmetrize',1)
15               nwfc1=nwfc1-1
                 IF (l == 0) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 2
                       work1(:) = work1(:) + d012 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    ENDDO
                 ELSEIF (l == 1) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 6
                       work1(:) = work1(:) + d112 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    ENDDO
                 ELSEIF (l == 2) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 10
                       work1(:) = work1(:) + d212 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    ENDDO
                 ELSEIF (l == 3) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 14
                       work1(:) = work1(:) + d312 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         work1(ibnd) * conjg (work1(ibnd)) / nsym
                 ENDDO
                 ! on symmetries
              ENDDO
           ENDIF
           ! on atomic wavefunctions
        ENDDO
     ELSE
        DO nwfc=1,natomwfc
           DO ibnd=1,nbnd
              proj(nwfc,ibnd,ik)=abs(proj0(nwfc,ibnd))**2
           ENDDO
        ENDDO
     ENDIF
     ! on k-points
  ENDDO
  !
  DEALLOCATE (work)
  DEALLOCATE (work1)
  DEALLOCATE (proj0)
  DEALLOCATE (e)
  CALL deallocate_bec_type (becp)
  DEALLOCATE (overlap)
  DEALLOCATE (wfcatom)
  IF (.not. lda_plus_u) DEALLOCATE (swfcatom)
  !
  !   vectors et and proj are distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (et,       nbnd, nkstot, nks)
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks)
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  !
  IF ( ionode ) THEN
     !
     ! write on the file filproj
     !
     IF (filproj/=' ') THEN
        iunproj=33
        CALL write_io_header(filproj, iunproj, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
           dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
           nkstot,nbnd,natomwfc)
        DO nwfc = 1, natomwfc
           IF (lspinorb) THEN
              WRITE(iunproj,1000) &
                   nwfc, nlmchi(nwfc)%na,atm(ityp(nlmchi(nwfc)%na)), &
                   nlmchi(nwfc)%n,nlmchi(nwfc)%jj,nlmchi(nwfc)%l, &
                   compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
           ELSE
              WRITE(iunproj,1500) &
                   nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                   nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m, &
                   0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
           ENDIF
           DO ik=1,nkstot
              DO ibnd=1,nbnd
                 WRITE(iunproj,'(2i8,f20.10)') ik,ibnd,abs(proj(nwfc,ibnd,ik))
              ENDDO
           ENDDO
        ENDDO
        CLOSE(iunproj)
     ENDIF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     !
     ! write on the standard output file
     !
     WRITE( stdout,'(/5x,"Atomic states used for projection")')
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
     IF (lspinorb) THEN
        DO nwfc = 1, natomwfc
           WRITE(stdout,1000) &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%jj, nlmchi(nwfc)%l,   &
             compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
        ENDDO
1000    FORMAT (5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                   " (j=",f3.1," l=",i1," m_j=",f4.1,")")
     ELSE
        DO nwfc = 1, natomwfc
           WRITE(stdout,1500) &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m, &
             0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
        ENDDO
1500    FORMAT (5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                   " (l=",i1," m=",i2," s_z=",f4.1,")")
     ENDIF
     !
     ALLOCATE(idx (natomwfc), proj1 (natomwfc) )
     DO ik = 1, nkstot
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
        DO ibnd = 1, nbnd
           WRITE( stdout, '("==== e(",i4,") = ",f11.5," eV ==== ")') &
              ibnd, et (ibnd, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           DO nwfc = 1, natomwfc
              idx (nwfc) = 0
              proj1 (nwfc) = - proj (nwfc, ibnd, ik)
           ENDDO
           CALL hpsort_eps (natomwfc, proj1, idx, eps4)
           !
           !  only projections that are larger than 0.001 are written
           !
           DO nwfc = 1, natomwfc
              proj1 (nwfc) = - proj1(nwfc)
              IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
           ENDDO
           nwfc = natomwfc + 1
20         nwfc = nwfc -1
           !
           ! fancy (?!?) formatting
           !
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i3,"]+"))') &
                (proj1 (i), idx(i), i = 1, min(5,nwfc))
           DO j = 1, (nwfc-1)/5
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i3,"]+"))') &
                   (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
           ENDDO
           psum = 0.d0
           DO nwfc = 1, natomwfc
              psum = psum + proj (nwfc, ibnd, ik)
           ENDDO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        ENDDO
     ENDDO
     DEALLOCATE (idx, proj1)
     !
     ! estimate partial charges (Loewdin) on each atom
     !
     IF (lspinorb) THEN
        nspin0 = 1
     ELSE
        nspin0 = 2
     ENDIF
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin0 ) )
     charges = 0.0d0
     IF (lspinorb) THEN
        DO ik = 1, nkstot
           DO ibnd = 1, nbnd
              DO nwfc = 1, natomwfc
                 na= nlmchi(nwfc)%na
                 l = nlmchi(nwfc)%l
                 charges(na,l,1) = charges(na,l,1) + &
                      wg (ibnd,ik) * proj (nwfc, ibnd, ik)
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO ik = 1, nkstot
           DO ibnd = 1, nbnd
              DO nwfc = 1, natomwfc
                 na= nlmchi(nwfc)%na
                 l = nlmchi(nwfc)%l
                 IF ( nlmchi(nwfc)%ind<=(2*l+1)) THEN
                    charges(na,l,1) = charges(na,l,1) + &
                         wg (ibnd,ik) * proj (nwfc, ibnd, ik)
                 ELSE
                    charges(na,l,2) = charges(na,l,2) + &
                         wg (ibnd,ik) * proj (nwfc, ibnd, ik)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
     WRITE( stdout, '(/"Lowdin Charges: "/)')
     !
     DO na = 1, nat
        DO is = 1, nspin0
           totcharge(is) = sum(charges(na,0:lmax_wfc,is))
        ENDDO
        IF ( nspin0 == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( charges(na,l,1), l= 0,lmax_wfc)
        ELSEIF ( nspin0 == 2) THEN
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( charges(na,l,1), l= 0,lmax_wfc)
           WRITE( stdout, 2002) totcharge(2), &
                ( charges(na,l,2), l= 0,lmax_wfc)
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                 ( charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        ENDIF
     ENDDO
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4 ,&
          & ", s, p, d, f = ",4f8.4)
2001 FORMAT (15x,"  spin up      = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4)
2002 FORMAT (15x,"  spin down    = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4)
2003 FORMAT (15x,"  polarization = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4)
     !
     psum = sum(charges(:,:,:)) / nelec
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     DEALLOCATE (charges)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE projwave_nc
!
!-----------------------------------------------------------------------
SUBROUTINE  partialdos (Emin, Emax, DeltaE, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev
  !
  USE projections
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  INTEGER :: ik, ibnd,  m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff
  REAL(DP) :: etev, delta, Elw, Eup, wkeff
  REAL(DP), ALLOCATABLE :: dostot(:,:,:), pdos(:,:,:,:), pdostot(:,:,:), &
       ldos(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  DO ik = 2, nkstot
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  ENDDO
  IF (degauss/=0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = max (Emin/rytoev, Elw)
  Emax = min (Emax/rytoev, Eup)
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (kresolveddos) THEN
     IF ( nspin==2 ) THEN
        nkseff=nkstot/2
     ELSE
        nkseff=nkstot
     ENDIF
  ELSE
     nkseff=1
  ENDIF
  !
  ALLOCATE (pdos(0:ne,natomwfc,nspin,nkseff))
  ALLOCATE (dostot(0:ne,nspin,nkseff), pdostot(0:ne,nspin,nkseff), ldos(0:ne,nspin,nkseff) )
  pdos(:,:,:,:) = 0.d0
  dostot(:,:,:) = 0.d0
  pdostot(:,:,:)= 0.d0
  !
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1

  DO ik = 1,nkstot
     !
     IF (kresolveddos) THEN
        ! set equal weight to all k-points
        wkeff=1.D0
        !
        IF (( nspin==2 ).AND.( isk(ik)==2 )) THEN
           ikeff=ik-nkstot/2
        ELSE
           ikeff=ik
        ENDIF
     ELSE
        ! use true weights
        wkeff=wk(ik)
        ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
        ikeff=1
     ENDIF
     !
     IF ( nspin == 2 ) current_spin = isk ( ik )
     DO ibnd = 1, nbnd
        etev = et(ibnd,ik)
        ie_mid = nint( (etev-Emin)/DeltaE )
        DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           ! pdos(:,nwfc,ns,ik) = DOS (states/eV) for spin "ns"
           !                      projected over atomic wfc "nwfc"
           !                      for k-point "ik" (or summed over all kp)
           !
           DO nwfc = 1, natomwfc
              pdos(ie,nwfc,current_spin,ikeff) = pdos(ie,nwfc,current_spin,ikeff) + &
                   wkeff * delta * proj (nwfc, ibnd, ik)
           ENDDO
           !
           ! dostot(:,ns,ik) = total DOS (states/eV) for spin "ns"
           !                   for k-point "ik" (or summed over all kp)
           !
           dostot(ie,current_spin,ikeff) = dostot(ie,current_spin,ikeff) + &
                wkeff * delta
        ENDDO
     ENDDO
  ENDDO
  !
  ! pdostot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin
        DO ie=0,ne
           pdostot(ie,is,ik) = sum(pdos(ie,:,is,ik))
        ENDDO
     ENDDO
  ENDDO

  DO nwfc = 1, natomwfc
     IF (nlmchi(nwfc)%m == 1) THEN
        filextension='.pdos_atm#'
        !             12345678901
        c_tab = 11
        IF (nlmchi(nwfc)%na < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na
           c_tab = c_tab + 1
        ELSEIF (nlmchi(nwfc)%na < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na
           c_tab = c_tab + 2
        ELSEIF (nlmchi(nwfc)%na < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na
           c_tab = c_tab + 3
        ELSEIF (nlmchi(nwfc)%na < 10000) THEN
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na
           c_tab = c_tab + 4
        ELSE
           CALL errore('partialdos',&
                'file extension not supporting so many atoms', nwfc)
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') &
             '(',trim(atm(ityp(nlmchi(nwfc)%na)))
        c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1
        IF (nlmchi(nwfc)%n >= 10) &
             CALL errore('partialdos',&
             'file extension not supporting so many atomic wfc', nwfc)
        IF (nlmchi(nwfc)%l > 3) &
             CALL errore('partialdos',&
             'file extension not supporting so many l', nwfc)
        WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
        fileout = trim(filpdos)//trim(filextension)
        OPEN (4,file=fileout,form='formatted', &
             status='unknown')

        IF (kresolveddos) THEN
           WRITE (4,'("# ik   ",$)')
        ELSE
           WRITE (4,'("#",$)')
        ENDIF
        IF (nspin == 1) THEN
           WRITE (4,'(" E (eV)   ldos(E)  ",$)')
        ELSE
           WRITE (4,'(" E (eV)  ldosup(E)  ldosdw(E)",$)')
        ENDIF
        DO m=1,2 * nlmchi(nwfc)%l + 1
           IF (nspin == 1) THEN
              WRITE(4,'(" pdos(E)   ",$)')
           ELSE
              WRITE(4,'(" pdosup(E) ",$)')
              WRITE(4,'(" pdosdw(E) ",$)')
           ENDIF
        ENDDO
        WRITE(4,*)
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:,:) = 0.d0
        DO ik=1,nkseff
           DO ie= 0, ne
              DO is=1, nspin
                 DO m=1,2 * nlmchi(nwfc)%l + 1
                    ldos  (ie, is, ik) = ldos  (ie, is, ik) + pdos(ie,nwfc+m-1,is,ik)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        DO ik=1,nkseff
           DO ie= 0, ne
              IF (kresolveddos) THEN
                 WRITE (4,'(i5," ",$)') ik
              ENDIF
              etev = Emin + ie * DeltaE
              WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  &
                   (ldos(ie,is,ik), is=1,nspin), &
                   ((pdos(ie,nwfc+m-1,is,ik), is=1,nspin), &
                   m=1,2*nlmchi(nwfc)%l+1)
           ENDDO
           IF (kresolveddos) WRITE (4,*)
        ENDDO
        CLOSE (4)
     ENDIF
  ENDDO
  fileout = trim(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown')
  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ",$)')
  ELSE
     WRITE (4,'("#",$)')
  ENDIF
  IF (nspin == 1) THEN
     WRITE (4,'(" E (eV)  dos(E)    pdos(E)")')
  ELSE
     WRITE (4,'(" E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)")')
  ENDIF
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ",$)') ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4e11.3)') etev*rytoev, (dostot(ie,is,ik), is=1,nspin), &
             (pdostot(ie,is,ik), is=1,nspin)
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi)
  DEALLOCATE (proj)
  DEALLOCATE (proj_aux)
  !
  RETURN
END SUBROUTINE partialdos
!
!-----------------------------------------------------------------------
SUBROUTINE  partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev
  !
  USE spin_orb,   ONLY: lspinorb
  USE projections_nc
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  INTEGER :: ik, ibnd, ind, m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff, nspin0
  REAL(DP) :: etev, delta, Elw, Eup, wkeff, fact(2), spinor
  REAL(DP), ALLOCATABLE :: dostot(:,:), pdos(:,:,:,:), pdostot(:,:,:), &
       ldos(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  DO ik = 2, nkstot
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  ENDDO
  IF (degauss/=0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = max (Emin/rytoev, Elw)
  Emax = min (Emax/rytoev, Eup)
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (lspinorb) THEN
     nspin0 = 1
  ELSE
     nspin0 = 2
  ENDIF
  !
  IF (kresolveddos) THEN
     nkseff=nkstot
  ELSE
     nkseff=1
  ENDIF
  !
  ALLOCATE (pdos(0:ne,natomwfc,nspin0,nkseff))
  ALLOCATE (dostot(0:ne,nkseff), pdostot(0:ne,nspin0,nkseff), ldos(0:ne,nspin0,nkseff) )
  pdos(:,:,:,:) = 0.d0
  dostot(:,:) = 0.d0
  pdostot(:,:,:)= 0.d0
  ie_delta = 5 * degauss / DeltaE + 1

  DO ik = 1,nkstot
     !
     IF (kresolveddos) THEN
        ! set equal weight to all k-points
        wkeff=1.D0
        ikeff=ik
     ELSE
        wkeff=wk(ik)
        ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
        ikeff=1
     ENDIF
     !
     DO ibnd = 1, nbnd
        etev = et(ibnd,ik)
        ie_mid = nint( (etev-Emin)/DeltaE )
        DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           ! pdos(:,nwfc,ns,ik) = DOS (states/eV) for spin "ns"
           !                      projected over atomic wfc "nwfc"
           !                      for k-point "ik" (or summed over all kp)
           !
           !
           ! dostot(:,ik) = total DOS (states/eV)
           !                for k-point "ik" (or summed over all kp)
           !
           IF (lspinorb) THEN
              DO nwfc = 1, natomwfc
                 pdos(ie,nwfc,1,ikeff) = pdos(ie,nwfc,1,ikeff) + &
                      wkeff * delta * proj (nwfc, ibnd, ik)
              ENDDO
              dostot(ie,ikeff) = dostot(ie,ikeff) + wkeff * delta
           ELSE
              DO nwfc = 1, natomwfc
                 IF ( nlmchi(nwfc)%ind<=(2* nlmchi(nwfc)%l+1)) THEN
                    pdos(ie,nwfc,1,ikeff) = pdos(ie,nwfc,1,ikeff) + &
                        wkeff * delta * proj (nwfc, ibnd, ik)
                    pdos(ie,nwfc,2,ikeff) = 0.d0
                 ELSE
                    pdos(ie,nwfc,1,ikeff) = 0.d0
                    pdos(ie,nwfc,2,ikeff) = pdos(ie,nwfc,2,ikeff) + &
                        wkeff * delta * proj (nwfc, ibnd, ik)
                 ENDIF
              ENDDO
              dostot(ie,ikeff) = dostot(ie,ikeff) + wkeff * delta
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  ! pdostot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin0
        DO ie=0,ne
           pdostot(ie,is,ik) = sum(pdos(ie,:,is,ik))
        ENDDO
     ENDDO
  ENDDO

  DO nwfc = 1, natomwfc
     IF (nlmchi(nwfc)%ind == 1) THEN
        filextension='.pdos_atm#'
        !             12345678901
        c_tab = 11
        IF (nlmchi(nwfc)%na < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na
           c_tab = c_tab + 1
        ELSEIF (nlmchi(nwfc)%na < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na
           c_tab = c_tab + 2
        ELSEIF (nlmchi(nwfc)%na < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na
           c_tab = c_tab + 3
        ELSEIF (nlmchi(nwfc)%na < 10000) THEN
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na
           c_tab = c_tab + 4
        ELSE
           CALL errore('partialdos_nc',&
                'file extension not supporting so many atoms', nwfc)
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') &
             '(',trim(atm(ityp(nlmchi(nwfc)%na)))
        c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1
        IF (nlmchi(nwfc)%n >= 10) &
             CALL errore('partialdos_nc',&
             'file extension not supporting so many atomic wfc', nwfc)
        IF (nlmchi(nwfc)%l > 3) &
             CALL errore('partialdos_nc',&
             'file extension not supporting so many l', nwfc)
        IF (lspinorb) THEN
           WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,"_j",f3.1,")")') &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l),nlmchi(nwfc)%jj
        ELSE
           WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
        ENDIF
        fileout = trim(filpdos)//trim(filextension)
        OPEN (4,file=fileout,form='formatted', &
             status='unknown')

        IF (kresolveddos) THEN
           WRITE (4,'("# ik   ",$)')
        ELSE
           WRITE (4,'("#",$)')
        ENDIF
        IF (nspin0 == 1) THEN
           WRITE (4,'(" E(eV)   ldos(E)   ",$)')
        ELSE
           WRITE (4,'(" E(eV)  ldosup(E)  ldosdw(E)",$)')
        ENDIF
        IF (lspinorb) THEN
           ind = 0
           DO m = -nlmchi(nwfc)%l-1, nlmchi(nwfc)%l
              fact(1) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,1)
              fact(2) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,2)
              IF (abs(fact(1))>1.d-8.or.abs(fact(2))>1.d-8) THEN
                 ind = ind + 1
                 WRITE(4,'("pdos(E)_",i1,"   ",$)') ind
              ENDIF
           ENDDO
        ELSE
           DO ind=1,2 * nlmchi(nwfc)%l + 1
              WRITE(4,'(" pdosup(E) ",$)')
              WRITE(4,'(" pdosdw(E) ",$)')
           ENDDO
        ENDIF
        WRITE(4,*)
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:,:) = 0.d0
        IF (lspinorb) THEN
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0)<1.d-8) THEN
                    DO ind = 1, 2 * nlmchi(nwfc)%l + 2
                       ldos  (ie, 1, ik) = ldos  (ie, 1, ik) + pdos(ie,nwfc+ind-1,1,ik)
                    ENDDO
                 ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0)<1.d-8) THEN
                    DO ind = 1, 2 * nlmchi(nwfc)%l
                       ldos  (ie, 1, ik) = ldos  (ie, 1, ik) + pdos(ie,nwfc+ind-1,1,ik)
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (kresolveddos) THEN
                    WRITE (4,'(i5," ",$)') ik
                 ENDIF
                 etev = Emin + ie * DeltaE
                 IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0)<1.d-8) THEN
                    WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1,ik), &
                         (pdos(ie,nwfc+ind-1,1,ik), ind=1,2*nlmchi(nwfc)%l+2)
                 ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0)<1.d-8) THEN
                    WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1,ik), &
                         (pdos(ie,nwfc+ind-1,1,ik), ind=1,2*nlmchi(nwfc)%l)
                 ENDIF
              ENDDO
              IF (kresolveddos) WRITE (4,*)
           ENDDO
        ELSE
           DO ik=1,nkseff
              DO ie= 0, ne
                 DO is=1, nspin0
                    DO ind=1,4 * nlmchi(nwfc)%l + 2
                       ldos  (ie, is, ik) = ldos  (ie, is, ik) + pdos(ie,nwfc+ind-1,is, ik)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (kresolveddos) THEN
                    WRITE (4,'(i5," ",$)') ik
                 ENDIF
                 etev = Emin + ie * DeltaE
                 WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  &
                      (ldos(ie,is,ik), is=1,nspin0), &
                      ((pdos(ie,nwfc+ind-1+(is-1)*(2*nlmchi(nwfc)%l+1),is,ik), is=1,nspin0), &
                      ind=1,2*nlmchi(nwfc)%l+1)
              ENDDO
              IF (kresolveddos) WRITE (4,*)
           ENDDO
        ENDIF
        CLOSE (4)
     ENDIF
  ENDDO
  fileout = trim(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown')
  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ",$)')
  ELSE
     WRITE (4,'("#",$)')
  ENDIF
  IF (nspin0 == 1) THEN
     WRITE (4,'(" E (eV)  dos(E)    pdos(E)")')
  ELSE
     WRITE (4,'(" E (eV)  dos(E)   pdosup(E)  pdosdw(E)")')
  ENDIF
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ",$)') ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4e11.3)') etev*rytoev, dostot(ie,ik), &
             (pdostot(ie,is,ik), is=1,nspin0)
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi)
  DEALLOCATE (proj)
  DEALLOCATE (proj_aux)
  !
  RETURN
END SUBROUTINE partialdos_nc
!
!-----------------------------------------------------------------------
SUBROUTINE write_io_header(filplot, iunplot, title, nr1x, nr2x, nr3x, &
           nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
           nkstot,nbnd,natomwfc)
   !-----------------------------------------------------------------------

   USE kinds, ONLY: DP
   USE ions_base, ONLY : zv, atm, tau, ityp
   USE noncollin_module, ONLY: noncolin
   USE spin_orb, ONLY: lspinorb

   IMPLICIT NONE
   CHARACTER (len=*) :: filplot
   CHARACTER (len=*) :: title
   INTEGER :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav
   REAL(DP) :: celldm (6), gcutm, dual, ecutwfc, at(3,3)
   INTEGER :: iunplot, ios, na, nt, i
   INTEGER :: nkstot,nbnd,natomwfc
   !
   IF (filplot == ' ') CALL errore ('write_io_h', 'filename missing', 1)

   OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
         STATUS = 'unknown', ERR = 101, IOSTAT = ios)
   101     CALL errore ('write_io_h', 'opening file '//trim(filplot), abs (ios) )
   WRITE (iunplot, '(a)') title
   WRITE (iunplot, '(8i8)') nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
   WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
   IF  (ibrav == 0) THEN
       WRITE ( iunplot, * ) at(:,1)
       WRITE ( iunplot, * ) at(:,2)
       WRITE ( iunplot, * ) at(:,3)
   ENDIF
   WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, 9
   WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
         (nt, atm (nt), zv (nt), nt=1, ntyp)
   WRITE (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
         (tau (i, na), i = 1, 3), ityp (na), na = 1, nat)
   WRITE (iunplot, '(3i8)') natomwfc, nkstot, nbnd
   WRITE (iunplot, '(2l5)') noncolin, lspinorb

   RETURN
END SUBROUTINE write_io_header
!
!-----------------------------------------------------------------------
FUNCTION compute_mj(j,l,m)
   !-----------------------------------------------------------------------
   USE kinds, ONLY: DP
   IMPLICIT NONE
   !
   REAL(DP) :: compute_mj, j
   INTEGER  :: l, m

   IF (abs(j-l-0.5d0)<1.d-4) THEN
       compute_mj=m+0.5d0
   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
      compute_mj=m-0.5d0
   ELSE
      CALL errore('compute_mj','l and j not compatible',1)
   ENDIF

   RETURN
END FUNCTION compute_mj
!
!-----------------------------------------------------------------------
SUBROUTINE  write_proj (filename, projs)
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir
  USE basis,            ONLY : natomwfc
  USE cell_base
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : nspin, isk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE iotk_module
  IMPLICIT NONE

  CHARACTER(*),  INTENT(in) :: filename
  COMPLEX(DP),   INTENT(in) :: projs(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256)          :: tmp
  CHARACTER(iotk_attlenx) :: attr
  INTEGER :: ik, ik_eff, ia, ierr, num_k_points

!
! subroutine body
!

  tmp = trim( tmp_dir ) // trim( prefix ) // '.save/' //trim(filename)
  !
  CALL iotk_open_write(iun, FILE=trim(tmp), ROOT="ATOMIC_PROJECTIONS", IERR=ierr )
  IF ( ierr /= 0 ) RETURN
  !
  !
  num_k_points = nkstot
  IF ( nspin == 2 ) num_k_points = nkstot / 2
  !
  CALL iotk_write_begin(iun, "HEADER")
  !
  CALL iotk_write_dat(iun, "NUMBER_OF_BANDS", nbnd)
  CALL iotk_write_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
  CALL iotk_write_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL iotk_write_dat(iun, "NON-COLINEAR_CALCULATION",noncolin)
  CALL iotk_write_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL iotk_write_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_ENERGY", ATTR=attr)
  CALL iotk_write_dat(iun, "FERMI_ENERGY", ef )
  !
  CALL iotk_write_end(iun, "HEADER")
  !
  !
  CALL iotk_write_dat(iun, "K-POINTS", xk(:,1:num_k_points) , COLUMNS=3 )
  CALL iotk_write_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points), COLUMNS=8 )
  !
  CALL iotk_write_begin(iun, "EIGENVALUES")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_dat( iun, "EIG.1", et(:,ik) )
        CALL iotk_write_dat( iun, "EIG.2", et(:,ik_eff) )
        !
     ELSE
        !
        CALL iotk_write_dat( iun, "EIG", et(:,ik) )
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "EIGENVALUES")

  !
  ! main loop atomic wfc
  !
  CALL iotk_write_begin(iun, "PROJECTIONS")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        CALL iotk_write_begin ( iun, "SPIN.1" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.1" )
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_begin ( iun, "SPIN.2" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.2" )
        !
     ELSE
        !
        DO ia = 1,natomwfc
            CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
        ENDDO
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "PROJECTIONS")

  !
  ! closing the file
  !
  CALL iotk_close_write(iun)

END SUBROUTINE write_proj


!
!  projwave with distributed matrixes
!
!-----------------------------------------------------------------------
SUBROUTINE pprojwave( filproj, lsym )
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE fft_base, ONLY : dfftp
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symm_base, ONLY: nsym, irt, d1, d2, d3
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE spin_orb, ONLY: lspinorb
  USE mp,       ONLY: mp_bcast
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, root_pool, &
                               intra_pool_comm, me_image, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id, &
                               leg_ortho, mpime
  USE wavefunctions_module, ONLY: evc
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst, zsqmcll, dsqmsym
  USE zhpev_module,     ONLY : pzhpev_drv, zhpev_drv
  USE descriptors,      ONLY : la_descriptor, descla_init
  USE projections
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  COMPLEX(DP), PARAMETER :: zero = ( 0.0d0, 0.0d0 )
  COMPLEX(DP), PARAMETER :: one  = ( 1.0d0, 0.0d0 )

  CHARACTER (len=*) :: filproj
  INTEGER :: ik, ibnd, i, j, na, nb, nt, isym, n,  m, m1, l, nwfc,&
       nwfc1, lmax_wfc, is, iunproj, iunaux
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: work1(:), proj0(:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap_d(:,:), work_d(:,:), diag(:,:), vv(:,:)
  COMPLEX(DP), ALLOCATABLE :: e_work_d(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE ::rwork1(:),rproj0(:,:)
  REAL   (DP), ALLOCATABLE ::roverlap_d(:,:)
  ! ... or for gamma-point.
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2)
  INTEGER  :: nksinit, nkslast
  CHARACTER(len=256) :: filename
  CHARACTER(len=256) :: auxname
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  INTEGER, ALLOCATABLE :: idx(:)
  LOGICAL :: lsym
  TYPE(la_descriptor) :: desc
  TYPE(la_descriptor), ALLOCATABLE :: desc_ip( :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx, nrl, nrlx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  !
  WRITE( stdout, '(/5x,"Calling pprojwave .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  ! Open file as temporary storage
  !
  iunaux = find_free_unit()
  WRITE( auxname, fmt='(I6.1)' ) mpime
  auxname = TRIM(tmp_dir) // TRIM(ADJUSTL(prefix)) // '.AUX' // TRIM(ADJUSTL(auxname))
  OPEN( unit=iunaux, file=trim(auxname), status='unknown', form='unformatted')
  !
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ) )
  ALLOCATE( notcnv_ip( np_ortho(2) ) )
  ALLOCATE( desc_ip( np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( natomwfc, desc, desc_ip )
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  !
  ! fill structure nlmchi
  !
  ALLOCATE (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           lmax_wfc = max (lmax_wfc, l )
           DO m = 1, 2 * l + 1
              nwfc=nwfc+1
              nlmchi(nwfc)%na = na
              nlmchi(nwfc)%n  =  n
              nlmchi(nwfc)%l  =  l
              nlmchi(nwfc)%m  =  m
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  IF (lmax_wfc > 3) CALL errore ('projwave', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projwave', 'wrong # of atomic wfcs?', 1)
  !
  !
  IF( ionode ) THEN
     WRITE( stdout, * )
     WRITE( stdout, * ) ' Problem Sizes '
     WRITE( stdout, * ) ' natomwfc = ', natomwfc
     WRITE( stdout, * ) ' nbnd     = ', nbnd
     WRITE( stdout, * ) ' nkstot   = ', nkstot
     WRITE( stdout, * ) ' npwx     = ', npwx
     WRITE( stdout, * ) ' nkb      = ', nkb
     WRITE( stdout, * )
  ENDIF
  !
  ALLOCATE( proj (natomwfc, nbnd, nkstot) )
  proj      = 0.d0
  !
  IF (.not. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) )
  ALLOCATE(wfcatom (npwx, natomwfc) )
  !
  ALLOCATE(e (natomwfc) )
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  DO ik = 1, nks
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)

     CALL atomic_wfc (ik, wfcatom)

     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL allocate_bec_type ( nkb, natomwfc, becp )
     CALL calbec ( npw, vkb, wfcatom, becp)

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     CALL deallocate_bec_type (becp)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     IF( la_proc ) THEN
        ALLOCATE(overlap_d (nx, nx) )
     ELSE
        ALLOCATE(overlap_d (1, 1) )
     ENDIF
     overlap_d = (0.d0,0.d0)
     IF ( gamma_only ) THEN
        IF( la_proc ) THEN
           ALLOCATE(roverlap_d (nx, nx) )
        ELSE
           ALLOCATE(roverlap_d (1, 1) )
        ENDIF
        roverlap_d = 0.d0
        CALL calbec_ddistmat( npw, wfcatom, swfcatom, natomwfc, nx, roverlap_d )
        overlap_d(:,:)=cmplx(roverlap_d(:,:),0.0_dp, kind=dp)
        ! TEMP: diagonalization routine for real matrix should be used instead
     ELSE
        CALL calbec_zdistmat( npw, wfcatom, swfcatom, natomwfc, nx, overlap_d )
     ENDIF
     !
     ! calculate O^{-1/2}
     !
     IF ( desc%active_node > 0 ) THEN
        !
        !  Compute local dimension of the cyclically distributed matrix
        !
        ALLOCATE(work_d (nx, nx) )

        nrl  = desc%nrl
        nrlx = desc%nrlx

        ALLOCATE( diag( nrlx, natomwfc ) )
        ALLOCATE( vv( nrlx, natomwfc ) )
        !
        CALL blk2cyc_zredist( natomwfc, diag, nrlx, natomwfc, overlap_d, nx, nx, desc )
        !
        CALL pzhpev_drv( 'V', diag, nrlx, e, vv, nrlx, nrl, natomwfc, &
           desc%npc * desc%npr, desc%mype, desc%comm )
        !
        CALL cyc2blk_zredist( natomwfc, vv, nrlx, natomwfc, work_d, nx, nx, desc )
        !
        DEALLOCATE( vv )
        DEALLOCATE( diag )
        !
     ELSE
        ALLOCATE(work_d (1, 1) )
     ENDIF

     CALL mp_bcast( e, root_pool, intra_pool_comm )

     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO

     IF ( desc%active_node > 0 ) THEN
        ALLOCATE(e_work_d (nx, nx) )
        DO j = 1, desc%nc
           DO i = 1, desc%nr
              e_work_d( i, j ) = e( j + desc%ic - 1 ) * work_d( i, j )
           ENDDO
        ENDDO
        CALL sqr_zmm_cannon( 'N', 'C', natomwfc, ONE, e_work_d, nx, work_d, nx, ZERO, overlap_d, nx, desc )
        CALL zsqmher( natomwfc, overlap_d, nx, desc )
        DEALLOCATE( e_work_d )
     ENDIF
     !
     DEALLOCATE( work_d )
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     IF ( gamma_only ) THEN
        ! TEMP: diagonalization routine for real matrix should be used instead
        roverlap_d(:,:)=REAL(overlap_d(:,:),DP)
        CALL wf_times_roverlap( swfcatom, roverlap_d, wfcatom )
        DEALLOCATE( roverlap_d )
     ELSE
        CALL wf_times_overlap( swfcatom, overlap_d, wfcatom )
        DEALLOCATE( overlap_d )
     ENDIF

     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>
     !
     IF ( gamma_only ) THEN
        !
        ALLOCATE( rproj0(natomwfc,nbnd), rwork1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, rproj0)
        !
        WRITE( iunaux ) rproj0
        !
     ELSE
        !
        ALLOCATE( proj0(natomwfc,nbnd), work1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, proj0)
        !
        WRITE( iunaux ) proj0
        !
     ENDIF
     !
     ! symmetrize the projections
     !
     IF (lsym) THEN
        !
        DO nwfc = 1, natomwfc
           !
           !  atomic wavefunction nwfc is on atom na
           !
           na= nlmchi(nwfc)%na
           n = nlmchi(nwfc)%n
           l = nlmchi(nwfc)%l
           m = nlmchi(nwfc)%m
           !
           DO isym = 1, nsym
              !
              nb = irt (isym, na)
              DO nwfc1 =1, natomwfc
                 IF (nlmchi(nwfc1)%na == nb             .and. &
                      nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                      nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                      nlmchi(nwfc1)%m == 1 ) GOTO 10
              ENDDO
              CALL errore('projwave','cannot symmetrize',1)
10            nwfc1=nwfc1-1
              !
              !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
              !
              IF ( gamma_only ) THEN
                 IF (l == 0) THEN
                    rwork1(:) = rproj0 (nwfc1 + 1,:)
                 ELSEIF (l == 1) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 3
                       rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 2) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 5
                       rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 3) THEN
                    rwork1(:) = 0.d0
                    DO m1 = 1, 7
                       rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         rwork1(ibnd) * rwork1(ibnd) / nsym
                 ENDDO
              ELSE
                 IF (l == 0) THEN
                    work1(:) = proj0 (nwfc1 + 1,:)
                 ELSEIF (l == 1) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 3
                       work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 2) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 5
                       work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ELSEIF (l == 3) THEN
                    work1(:) = 0.d0
                    DO m1 = 1, 7
                       work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:)
                    ENDDO
                 ENDIF
                 DO ibnd = 1, nbnd
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                         work1(ibnd) * conjg (work1(ibnd)) / nsym
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
        !
     ELSE
        !
        IF ( gamma_only ) THEN
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(rproj0(nwfc,ibnd))**2
              ENDDO
           ENDDO
        ELSE
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(proj0(nwfc,ibnd))**2
              ENDDO
           ENDDO
        ENDIF
        !
     ENDIF
     !
     IF ( gamma_only ) THEN
        DEALLOCATE (rwork1)
        DEALLOCATE (rproj0)
     ELSE
        DEALLOCATE (work1)
        DEALLOCATE (proj0)
     ENDIF
     !
  ENDDO
  !
  !
  DEALLOCATE (e)
  !
  DEALLOCATE (wfcatom)
  IF (.not. lda_plus_u) DEALLOCATE (swfcatom)
  !
  CLOSE( unit=iunaux )
  !
  !
  !   vectors et and proj are distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (et,       nbnd, nkstot, nks)
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks)

  !
  !  Recover proj_aux
  !
  OPEN( unit=iunaux, file=trim(auxname), status='old', form='unformatted')

  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj_aux  = (0.d0, 0.d0)
  DO ik = 1, nks
     !
     IF( gamma_only ) THEN
        ALLOCATE( rproj0( natomwfc, nbnd ) )
        READ( iunaux ) rproj0(:,:)
        proj_aux(:,:,ik) = cmplx( rproj0(:,:), 0.00_dp, kind=dp )
        DEALLOCATE ( rproj0 )
     ELSE
        READ( iunaux ) proj_aux(:,:,ik)
     ENDIF
     !
  ENDDO
  !
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  !
  CLOSE( unit=iunaux, status='delete' )
  !
  IF ( ionode ) THEN
     !
     ! write on the file filproj
     !
     IF (filproj/=' ') THEN
        DO is=1,nspin
           IF (nspin==2) THEN
              IF (is==1) filename=trim(filproj)//'.up'
              IF (is==2) filename=trim(filproj)//'.down'
              nksinit=(nkstot/2)*(is-1)+1
              nkslast=(nkstot/2)*is
           ELSE
              filename=trim(filproj)
              nksinit=1
              nkslast=nkstot
           ENDIF
           iunproj=33
           CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, &
                ecutwfc, nkstot/nspin,nbnd,natomwfc)
           DO nwfc = 1, natomwfc
              WRITE(iunproj,'(2i5,a3,3i5)') &
                  nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                  nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
              DO ik=nksinit,nkslast
                 DO ibnd=1,nbnd
                   WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, &
                                                 abs(proj(nwfc,ibnd,ik))
                 ENDDO
              ENDDO
           ENDDO
           CLOSE(iunproj)
        ENDDO
     ENDIF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     !
     ! write on the standard output file
     !
     WRITE( stdout,'(/5x,"Atomic states used for projection")')
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
     DO nwfc = 1, natomwfc
        WRITE(stdout,1000) &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
     ENDDO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")")
     !
     ALLOCATE(idx(natomwfc), proj1 (natomwfc) )
     !
     DO ik = 1, nkstot
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
        DO ibnd = 1, nbnd
           WRITE( stdout, '(5x,"e = ",f11.5," eV")') et (ibnd, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           DO nwfc = 1, natomwfc
              idx (nwfc) = 0
              proj1 (nwfc) = - proj (nwfc, ibnd, ik)
           ENDDO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (natomwfc, proj1, idx, eps4)
           !
           !  only projections that are larger than 0.001 are written
           !
           DO nwfc = 1, natomwfc
              proj1 (nwfc) = - proj1(nwfc)
              IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
           ENDDO
           nwfc = natomwfc + 1
20         nwfc = nwfc -1
           !
           ! fancy (?!?) formatting
           !
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') &
                (proj1 (i), idx(i), i = 1, min(5,nwfc))
           DO j = 1, (nwfc-1)/5
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') &
                   (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
           ENDDO
           psum = 0.d0
           DO nwfc = 1, natomwfc
              psum = psum + proj (nwfc, ibnd, ik)
           ENDDO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        ENDDO
     ENDDO
     !
     DEALLOCATE (idx, proj1)
     !
     ! estimate partial charges (Loewdin) on each atom
     !
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin ) )
     charges = 0.0d0
     DO ik = 1, nkstot
        IF ( nspin == 1 ) THEN
           current_spin = 1
        ELSEIF ( nspin == 2 ) THEN
           current_spin = isk ( ik )
        ELSE
           CALL errore ('projwfc_nc',' called in the wrong case ',1)
        ENDIF
        DO ibnd = 1, nbnd
           DO nwfc = 1, natomwfc
              na= nlmchi(nwfc)%na
              l = nlmchi(nwfc)%l
              charges(na,l,current_spin) = charges(na,l,current_spin) + &
                   wg (ibnd,ik) * proj (nwfc, ibnd, ik)
           ENDDO
        ENDDO
     ENDDO
     !
     WRITE( stdout, '(/"Lowdin Charges: "/)')
     !
     DO na = 1, nat
        DO is = 1, nspin
           totcharge(is) = sum(charges(na,0:lmax_wfc,is))
        ENDDO
        IF ( nspin == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
        ELSEIF ( nspin == 2) THEN
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
           WRITE( stdout, 2002) totcharge(2), &
                ( l_label(l), charges(na,l,2), l= 0,lmax_wfc)
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        ENDIF
     ENDDO
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4,4(", ",a1," =",f8.4))
2001 FORMAT (15x,"  spin up      = ",f8.4,4(", ",a1," =",f8.4))
2002 FORMAT (15x,"  spin down    = ",f8.4,4(", ",a1," =",f8.4))
2003 FORMAT (15x,"  polarization = ",f8.4,4(", ",a1," =",f8.4))
     !
     psum = sum(charges(:,:,:)) / nelec
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     DEALLOCATE (charges)
     !
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(in)  :: nsiz
     TYPE(la_descriptor), INTENT(out) :: desc
     TYPE(la_descriptor), INTENT(out) :: desc_ip(:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
     !
     nx = desc%nrcx
     !
     DO j = 0, desc%npc - 1
        DO i = 0, desc%npr - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(i+1,j+1), desc%n, desc%nx, np_ortho, coor_ip, ortho_comm, 1 )
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        ENDDO
     ENDDO
     !
     la_proc = .false.
     IF( desc%active_node > 0 ) la_proc = .true.
     !
     RETURN
  END SUBROUTINE desc_init
  !

  SUBROUTINE calbec_zdistmat( npw, v, w, n, nx, dm )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     USE mp, ONLY : mp_root_sum
     !
     IMPLICIT NONE
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, ldv, ldw
     INTEGER, INTENT(in) :: npw ! local number of plane wave
     INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
     COMPLEX(DP), INTENT(out) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = zero
     !
     ldv = size( v, 1 )
     ldw = size( w, 1 )
     !
     DO ipc = 1, desc%npc !  loop on column procs
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, npw, ONE , &
                       v(1,ir), ldv, w(1,ic), ldw, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        ENDDO
        !
     ENDDO
     !
     CALL zsqmher( n, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE calbec_zdistmat
  !

  SUBROUTINE calbec_ddistmat( npw, v, w, n, nx, dm )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     USE mp, ONLY : mp_root_sum
     USE gvect, ONLY : gstart
     !
     IMPLICIT NONE
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, ldv, ldw, npw2, npwx2
     INTEGER, INTENT(in) :: npw ! local number of plane wave
     INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
     REAL(DP), INTENT(out) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     npw2  = 2*npw
     npwx2 = 2*npwx
     !
     work = zero
     !
     ldv = size( v, 1 )
     ldw = size( w, 1 )
     !
     DO ipc = 1, desc%npc !  loop on column procs
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0 , &
                       v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        ENDDO
        !
     ENDDO
     !
     CALL dsqmsym( n, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE calbec_ddistmat
  !
  !
  !
  SUBROUTINE wf_times_overlap( swfc, ovr, wfc )

     COMPLEX(DP) :: swfc( :, : ), ovr( :, : ), wfc( :, : )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        beta = ZERO

        DO ipr = 1, desc%npr
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc%myr .and. ipc-1 == desc%myc .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, ONE, &
                          swfc(1,ir), npwx, ovr, nx, beta, wfc(1,ic), npwx )
           ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, ONE, &
                       swfc(1,ir), npwx, vtmp, nx, beta, wfc(1,ic), npwx )
           ENDIF
           !
           beta = ONE

        ENDDO
        !
     ENDDO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_overlap

  !
  SUBROUTINE wf_times_roverlap( swfc, ovr, wfc )

     USE gvect, ONLY : gstart

     COMPLEX(DP) :: swfc( :, : ), wfc( :, : )
     REAL(DP)    :: ovr( :, : )
     !
     INTEGER :: ipc, ipr, npw2, npwx2
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     npw2  = 2*npw
     npwx2 = 2*npwx

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        beta = 0.0d0

        DO ipr = 1, desc%npr
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc%myr .and. ipc-1 == desc%myc .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, ovr, nx, beta, wfc(1,ic), npwx2 )
              !
           ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, vtmp, nx, beta, wfc(1,ic), npwx2 )
              !
           ENDIF
           !
           beta = 1.0d0

        ENDDO
        !
     ENDDO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_roverlap
  !
END SUBROUTINE pprojwave
!
!-----------------------------------------------------------------------
SUBROUTINE projwave_boxes( filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes )
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE atom
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev
  USE gvect
  USE gvecs,   ONLY: dual
  USE klist, ONLY: xk, nks, nkstot
  USE lsda_mod, ONLY: nspin, isk, current_spin, lsda
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: okvan
  USE noncollin_module, ONLY: noncolin, npol
  USE wavefunctions_module, ONLY: evc,    psic
  USE wavefunctions_module, ONLY: psic_nc
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE scf,                  ONLY : rho
  USE projections_ldos
  USE fft_base,             ONLY : grid_scatter, dfftp
  USE fft_interfaces,       ONLY : invfft
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
!
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  CHARACTER (len=256) :: filpdos
  CHARACTER (len=*) :: filproj
  INTEGER :: n_proj_boxes, irmin(3,*), irmax(3,*)
  LOGICAL :: plotboxes
  !
  INTEGER :: ik, ibnd, i, ir, ig, ipol, ibox, ir1, ir2, ir3, c_tab, is, iunproj
  INTEGER :: nri(3)
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  COMPLEX(DP), ALLOCATABLE :: caux(:)
  REAL(DP), ALLOCATABLE :: thetabox(:), raux(:), thetathisproc(:,:), union(:), intersection(:)
  LOGICAL, ALLOCATABLE :: isInside(:,:)
  REAL(DP), EXTERNAL :: DDOT
  REAL(DP), ALLOCATABLE :: boxvolume(:), boxcharge(:)
  !
  WRITE( stdout, '(/5x,"Calling projwave_boxes .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  IF (noncolin) THEN
     WRITE( stdout, '(/5x,"Non spin-resolved DOS will be computed")')
  ENDIF
  !
  IF (okvan) THEN
     CALL errore( 'projwave_boxes', 'Augmentation contributions are currently not included to the DOS in boxes',-1)
  ENDIF
  !
  IF ( ( n_proj_boxes > N_MAX_BOXES ) .or. ( n_proj_boxes < 1 ) ) &
     CALL errore ('projwave_boxes', 'n_proj_boxes not correct', abs (n_proj_boxes) )
  !
  ! ... Define functions with values 1.0
  ! ... on the specified boxes and 0.0 elsewhere.
  !
  ALLOCATE( thetabox (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  !
  ALLOCATE( thetathisproc(dfftp%nnr,1:n_proj_boxes) )
  !
  ALLOCATE ( isInside ( max(dfftp%nr1,dfftp%nr2,dfftp%nr3), 3 ) )
  !
  DO ibox = 1, n_proj_boxes
     !
     ! A. Do the three directions independently:
     nri(1)=dfftp%nr1
     nri(2)=dfftp%nr2
     nri(3)=dfftp%nr3
     DO i = 1, 3
        ! boxes include the points in [irmin,irmax] if irmin<=irmax
        ! and the points in [1,irmax] and [irmin,nr] if irmin > irmax
        irmin(i,ibox)=mod(irmin(i,ibox),nri(i))
        IF (irmin(i,ibox)<=0) irmin(i,ibox)=irmin(i,ibox)+nri(i)
        irmax(i,ibox)=mod(irmax(i,ibox),nri(i))
        IF (irmax(i,ibox)<=0) irmax(i,ibox)=irmax(i,ibox)+nri(i)
        DO ir = 1, nri(i)
           IF (irmin(i,ibox)<=irmax(i,ibox)) THEN
              isInside(ir,i)=(ir>=irmin(i,ibox)).and.(ir<=irmax(i,ibox))
           ELSE
              isInside(ir,i)=(ir>=irmin(i,ibox)).or. (ir<=irmax(i,ibox))
           ENDIF
        ENDDO
     ENDDO
     !
     ! B. Combine the conditions for the three directions to form a box
     ir=0
     DO ir3 = 1, dfftp%nr3
        DO ir2 = 1, dfftp%nr2
           DO ir1 = 1, dfftp%nr1
              ir=ir+1
              IF ( isInside(ir1,1) .and. &
                   isInside(ir2,2) .and. &
                   isInside(ir3,3)         ) THEN
                 thetabox(ir)=1._DP
              ELSE
                 thetabox(ir)=0._DP
              ENDIF
           ENDDO
        ENDDO
        !
     ENDDO
     !
     ! C. Output the functions thetabox in the XCrySDen format,
     ! so that the projection boxes can be visualised.
     IF ( ionode .and. plotboxes ) THEN
        filextension='.box#'
        !             123456
        c_tab = 6
        IF (ibox < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') ibox
           c_tab = c_tab + 1
        ELSEIF (ibox < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') ibox
           c_tab = c_tab + 2
        ELSEIF (ibox < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') ibox
           c_tab = c_tab + 3
        ELSE
           CALL errore('projwave_boxes',&
                'file extension not supporting so many boxes', n_proj_boxes)
        ENDIF
        !
        fileout = trim(filpdos)//trim(filextension)//'.xsf'
        OPEN (4,file=fileout,form='formatted', status='unknown')
        CALL xsf_struct (alat, at, nat, tau, atm, ityp, 4)
        CALL xsf_fast_datagrid_3d(thetabox(1:dfftp%nr1x*dfftp%nr2x*dfftp%nr3x),&
                 dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, 4)
        CLOSE (4)
        !
     ENDIF
     !
     CALL grid_scatter ( thetabox(:), thetathisproc(:,ibox) )
     !
  ENDDO
  !
  DEALLOCATE ( isInside )
  DEALLOCATE ( thetabox )
  !
  !
  ! ... For each box output the volume and the electronic charge contained
  !
  ALLOCATE ( boxvolume (1:n_proj_boxes) )
  ALLOCATE ( boxcharge (1:n_proj_boxes) )
  ALLOCATE ( raux (dfftp%nnr) )
  !
  ! A. Integrate the volume
  DO ibox = 1, n_proj_boxes
     boxvolume(ibox) = sum(thetathisproc(1:dfftp%nnr,ibox))
     CALL mp_sum ( boxvolume(ibox) , intra_pool_comm )
  ENDDO
  !
  ! B1. Copy the total charge density to raux
  IF (noncolin) THEN
     CALL DCOPY (dfftp%nnr, rho%of_r, 1, raux, 1)
  ELSE
     CALL DCOPY (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
     DO is = 2, nspin
        CALL DAXPY (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
     ENDDO
  ENDIF
  !
  ! B2. Integrate the charge
  !     the correct integral has dv = omega/(nr1*nr2*nr3)
  !     not  omega/(nr1x*nr2x*nr3x) . PG 24 Oct 2010
  DO ibox = 1, n_proj_boxes
     boxcharge(ibox) = DDOT(dfftp%nnr,raux(:),1,thetathisproc(:,ibox),1) &
          &   * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
     CALL mp_sum ( boxcharge(ibox) , intra_pool_comm )
  ENDDO
  !
  ! C. Write the result
  IF (ionode) THEN
     WRITE (stdout,*)
     DO ibox = 1, n_proj_boxes
        WRITE (stdout, &
             '(5x,"Box #",i3," : vol ",f10.6," % = ",f14.6," (a.u.)^3; ",e13.6," elec")') &
             ibox, 100* boxvolume(ibox) /(dfftp%nr1*dfftp%nr2*dfftp%nr3), &
             omega* boxvolume(ibox)/(dfftp%nr1*dfftp%nr2*dfftp%nr3), boxcharge(ibox)
     ENDDO
  ENDIF
  !
  DEALLOCATE ( boxvolume , boxcharge )
  !
  ! ... Here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge in the specified box
  !
  ALLOCATE( proj(1:n_proj_boxes,nbnd,nkstot) )
  proj(:,:,:)=0._DP
  !
  ALLOCATE( caux(dfftp%nnr) )
  !
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     bnd_loop: DO ibnd = 1, nbnd
        !
        IF (noncolin) THEN
           !
           psic_nc = (0.d0,0.d0)
           DO ig = 1, npw
              psic_nc(nl(igk(ig)),1)=evc(ig     ,ibnd)
              psic_nc(nl(igk(ig)),2)=evc(ig+npwx,ibnd)
           ENDDO
           raux=0._DP
           DO ipol=1,npol
              CALL invfft ('Dense', psic_nc(:,ipol), dfftp)
              raux(:) = raux(:)+dble( psic_nc(:,ipol) )**2 &
                             + aimag( psic_nc(:,ipol) )**2
           ENDDO
           !
        ELSE
           !
           caux(1:dfftp%nnr) = (0._DP,0._DP)
           DO ig = 1, npw
              caux (nl (igk (ig) ) ) = evc (ig, ibnd)
           ENDDO
           IF (gamma_only) THEN
              DO ig = 1, npw
                 caux (nlm(igk (ig) ) ) = conjg(evc (ig, ibnd))
              ENDDO
           ENDIF
           CALL invfft ('Dense', caux, dfftp)
           !
           raux(:) = dble( caux(:) )**2 + aimag( caux(:) )**2
           !
        ENDIF
        !
        ! The contribution of this wavefunction to the LDOS
        ! integrated in the volume is the projection of the
        ! squared wfc on a function =1 in the volume itself:
        !
        DO ibox = 1, n_proj_boxes
           proj(ibox,ibnd,ik) = DDOT(dfftp%nnr,raux(:),1,thetathisproc(:,ibox),1) &
                &               / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
        ENDDO
        !
     ENDDO bnd_loop
     !
     CALL mp_sum ( proj(:,:,ik) , intra_pool_comm )
     !
  ENDDO k_loop
  !
  DEALLOCATE ( caux )
  DEALLOCATE ( raux )
  DEALLOCATE ( thetathisproc )
  !
  !   vector proj is distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (proj, n_proj_boxes*nbnd, nkstot, nks)
  !
  ! Output the projections
  IF ( ionode ) THEN
     IF (filproj/=' ') THEN
        iunproj=33
        CALL write_io_header(filproj, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
           dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
           celldm, at, gcutm, dual, ecutwfc, nkstot,nbnd,natomwfc)
        DO ibox = 1, n_proj_boxes
           WRITE (iunproj,'(3i6)') ibox, n_proj_boxes
           WRITE (iunproj,'(i6,i6,f9.4,e13.6)') &
                ((ik,ibnd,et(ibnd,ik)*rytoev,proj(ibox,ibnd,ik),ibnd=1,nbnd),ik=1,nkstot)
        ENDDO
        CLOSE (iunproj)
     ENDIF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE projwave_boxes
!
!-----------------------------------------------------------------------
SUBROUTINE partialdos_boxes(Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev
  USE projections_ldos
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  INTEGER :: n_proj_boxes
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  !
  INTEGER :: ik, ibnd, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff, ibox, nspin0
  REAL(DP) :: etev, delta, Elw, Eup, wkeff
  REAL(DP), ALLOCATABLE :: dostot(:,:,:), dosbox(:,:,:,:), dosboxtot(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  DO ik = 2, nkstot
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  ENDDO
  IF (degauss/=0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = max (Emin/rytoev, Elw)
  Emax = min (Emax/rytoev, Eup)
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (nspin==2) THEN
     nspin0 = 2
  ELSE
     nspin0 = 1
  ENDIF
  !
  IF (kresolveddos) THEN
     IF ( nspin==2 ) THEN
        nkseff=nkstot/2
     ELSE
        nkseff=nkstot
     ENDIF
  ELSE
     nkseff=1
  ENDIF
  !
  ALLOCATE (dosbox(0:ne,1:n_proj_boxes,nspin0,nkseff))
  ALLOCATE (dostot(0:ne,nspin0,nkseff), dosboxtot(0:ne,nspin0,nkseff) )
  dosbox(:,:,:,:) = 0.d0
  dostot(:,:,:) = 0.d0
  dosboxtot(:,:,:)= 0.d0
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1
  !
  DO ik = 1,nkstot
     !
     IF (kresolveddos) THEN
        ! set equal weight to all k-points
        wkeff=1.D0
        !
        IF (( nspin==2 ).AND.( isk(ik)==2 )) THEN
           ikeff=ik-nkstot/2
        ELSE
           ikeff=ik
        ENDIF
     ELSE
        ! use true weights
        wkeff=wk(ik)
        ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
        ikeff=1
     ENDIF
     !
     IF ( nspin == 2 ) current_spin = isk ( ik )
     DO ibnd = 1, nbnd
        etev = et(ibnd,ik)
        ie_mid = nint( (etev-Emin)/DeltaE )
        DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           DO ibox = 1, n_proj_boxes
              dosbox(ie,ibox,current_spin,ikeff)               = &
                   dosbox(ie,ibox,current_spin,ikeff)          + &
                   wkeff * delta * proj (ibox, ibnd, ik)
           ENDDO
           !
           ! dostot(:,ns,ik) = total DOS (states/eV) for spin "ns"
           !                   for k-point "ik" (or summed over all kp)
           !
           dostot(ie,current_spin,ikeff) = dostot(ie,current_spin,ikeff) + &
                wkeff * delta
        ENDDO
     ENDDO
  ENDDO
  !
  ! dosboxtot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin0
        DO ie=0,ne
           dosboxtot(ie,is,ik) = sum(dosbox(ie,1:n_proj_boxes,is,ik))
        ENDDO
     ENDDO
  ENDDO
  !
  fileout = trim(filpdos)//'.ldos_boxes'
  !
  OPEN (4,file=fileout,form='formatted', &
       status='unknown')

  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ",$)')
  ELSE
     WRITE (4,'("#",$)')
  ENDIF
  IF (nspin0 == 2) THEN
     WRITE (4,'(" E (eV)  tot_up(E)  tot_dw(E)  totldos_up totldos_dw ",$)')
  ELSE
     WRITE (4,'(" E (eV)  tot(E)     totldos    ",$)')
  ENDIF
  DO ibox=1, n_proj_boxes
     IF (nspin0 == 2) THEN
        WRITE(4,'("#",i3," up(E) ",$)')  ibox
        WRITE(4,'("#",i3," dw(E) ",$)')  ibox
     ELSE
        WRITE(4,'("#",i3," (E)   ",$)')  ibox
     ENDIF
  ENDDO
  WRITE (4,*)
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ",$)') ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4(2e11.3),999(2e11.3))') etev*rytoev,  &
             dostot(ie,1:nspin0,ik), dosboxtot(ie,1:nspin0,ik), &
             ( dosbox(ie,ibox,1:nspin0,ik), ibox = 1, n_proj_boxes )
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (dostot, dosboxtot)
  DEALLOCATE (dosbox)
  !
  DEALLOCATE (proj)
  !
  RETURN
END SUBROUTINE partialdos_boxes
