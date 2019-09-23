  !-----------------------------------------------------------------------
  SUBROUTINE dynmat_asr(iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq, &
                              invs, s, irt, rtau, sumr) 
  !-----------------------------------------------------------------------
  !!
  !! read dynamical matrix for the q points, either in plain text or xml.
  !! iq_first, iq_first+1, ... iq_first+nq-1
  !!
  !-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  use io_files,         ONLY : prefix 
  USE cell_base,        ONLY : ibrav, celldm, omega, at, bg
  USE ions_base,        ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE elph2,            ONLY : dynq, zstar, epsi
  USE symm_base,        ONLY : nsym
  USE epwcom,           ONLY : dvscf_dir, lpolar, lifc, nqc1, nqc2, nqc3
  USE modes,            ONLY : nmodes
  USE control_flags,    ONLY : iverbosity
  USE noncollin_module, ONLY : nspin_mag
  USE io_epw,           ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                               read_dyn_mat, check_is_xml_file
  USE io_global,        ONLY : ionode, stdout
  USE constants_epw,    ONLY : cone, czero, twopi, rydcm1, eps6, zero
  USE io_var,           ONLY : iudyn
  USE wan2bloch,        ONLY : dynifc2blochc
  USE low_lvl,          ONLY : set_ndnmbr, eqvect_strict
  !
  IMPLICIT NONE
  !
  ! Input
  INTEGER, INTENT(in) :: iq_irr
  !! The index of the irreducible q point
  INTEGER, INTENT(in) :: nqc_irr
  !! The number of irreducible qpoints 
  INTEGER, INTENT(in) :: nq
  !! The number of q points in the star of q
  INTEGER, INTENT(in) :: iq_first
  !! The index of the first qpoint to be read in the uniform q-grid
  INTEGER, INTENT(in) :: imq
  !! Flag which tells whether we have to consider the -q vectors
  INTEGER, INTENT(in) :: isq(48)
  !! Index of q in the star for a given symmetry.
  INTEGER, INTENT(in) :: invs(48)
  !! Matrix of inversed symmetry 
  INTEGER, INTENT(in) :: s(3, 3, 48)
  !! Set of symmetry operations
  INTEGER, INTENT(in) :: irt(48, nat)
  !! For each atoms give the rotated atoms
  REAL(KIND = DP), INTENT(inout) :: sxq(3, 48)
  !! Symmetry matrix 
  REAL(KIND = DP), INTENT(inout) :: rtau(3, 48, nat)
  !! the relative position of the rotated atom to the original one
  REAL(KIND = DP), INTENT(inout) :: sumr(2, 3, nat, 3)
  !! Sum to impose the ASR
  !
  ! Local variables 
  LOGICAL :: found
  !! Is it found 
  LOGICAL :: lrigid
  !! ifc 
  LOGICAL :: lraman
  !! is raman present
  LOGICAL :: nog
  !! 
  LOGICAL :: is_xml_file
  !! Is the file XML
  CHARACTER(LEN = 3) :: atm
  !! 
  CHARACTER(LEN = 3) :: filelab
  !! 
  CHARACTER(LEN = 80) :: line
  !! 
  CHARACTER(LEN = 256) :: tempfile
  !! 
  INTEGER :: neig
  !! The total number of eigenvalues found
  INTEGER :: info
  !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
  INTEGER :: ifail(nmodes)
  !! Contains the indices of the eigenvectors that failed to converge
  INTEGER :: iwork(5 * nmodes)
  !! Integer work array
  INTEGER :: isym
  !! Symmetry number
  INTEGER :: m1, m2, m3
  !! Supercell index
  INTEGER :: sna
  !! Rotation + translation
  INTEGER :: jsym 
  !! Symmetry index
  INTEGER :: nsq
  !! nsq is the degeneracy of the small group for this iq in the star
  INTEGER :: sym_sgq(48)
  !! Degenerate sym 
  INTEGER :: current_iq
  !! the q index in the dynq matrix
  INTEGER :: mq 
  !! Size of q
  INTEGER :: ism1
  !! Inversion symmetry
  INTEGER :: imode
  !! i-mode
  INTEGER :: jmode
  !! j-mode
  INTEGER :: ntyp_
  !! 
  INTEGER :: nat_
  !! 
  INTEGER :: ibrav_
  !! 
  INTEGER :: ityp_
  !! 
  INTEGER :: ios
  !! 
  INTEGER :: iq, jq
  !! 
  INTEGER :: nt
  !! Number of type of atoms
  INTEGER :: na, nb
  !! 
  INTEGER :: naa, nbb
  !! 
  INTEGER :: nu, mu
  !! 
  INTEGER :: i, j 
  !! 
  INTEGER :: ipol, jpol
  !! 
  INTEGER :: nqs 
  !! 
  INTEGER :: axis
  !! 
  INTEGER :: nrws
  !! 
  INTEGER :: ierr
  !! Error status
  INTEGER, parameter :: ntypx = 10
  !! 
  INTEGER, PARAMETER :: nrwsx = 200
  !! 
  REAL(KIND = DP) :: rwork(7 * nmodes)
  !! Real work array
  REAL(KIND = DP) :: arg
  !! 
  REAL(KIND = DP) :: aq(3)
  !! 
  REAL(KIND = DP) :: saq(3)
  !! 
  REAL(KIND = DP) :: raq(3)
  !! 
  REAL(KIND = DP) :: xq(3)
  !! The rotated q-vector
  REAL(KIND = DP) :: scart(3, 3)
  !! 
  REAL(KIND = DP) :: massfac
  !! 
  REAL(KIND = DP) :: w1(nmodes)
  !! 
  REAL(KIND = DP) :: wtmp(nmodes)
  !! 
  REAL(KIND = DP) :: celldm_(6)
  !! 
  REAL(KIND = DP) :: amass_
  !! 
  REAL(KIND = DP) :: tau_(3)
  !! 
  REAL(KIND = DP) :: q(3, nqc1 * nqc2 * nqc3)
  !! 
  REAL(KIND = DP) :: dynr(2, 3, nat, 3, nat)
  !! 
  REAL(KIND = DP) :: sumz
  !! 
  REAL(KIND = DP) :: qout(3)
  !! 
  REAL(KIND = DP) :: amass2(ntypx)
  !! 
  REAL(KIND = DP) :: rws(0:3, nrwsx)
  !! 
  REAL(KIND = DP) :: atws(3, 3)
  !! WS cell dist. 
  REAL(KIND = DP) :: zstar_(3, 3, nat) 
  !!  Temporary array to store Z*
  REAL(KIND = DP) :: epsi_(3, 3)
  !! Temporary array to store dielectric function
  REAL(KIND = DP), ALLOCATABLE :: m_loc(:, :)
  !! 
  REAL(KIND = DP), ALLOCATABLE :: dchi_dtau(:, :, :, :)
  !! 
  COMPLEX(KIND = DP) :: cfac
  !!  
  COMPLEX(KIND = DP) :: cwork(2 * nmodes)
  !! Complex work array
  COMPLEX(KIND = DP) :: gamma(nmodes, nmodes) 
  !! the Gamma matrix for the symmetry operation on the dyn mat 
  COMPLEX(KIND = DP) :: dynq_tmp(nmodes, nmodes)
  !!  
  COMPLEX(KIND = DP) :: cz1(nmodes, nmodes)
  !!  
  COMPLEX(KIND = DP) :: dynp(nmodes * (nmodes + 1) / 2)
  !!  
  COMPLEX(KIND = DP), ALLOCATABLE :: dyn(:, :, :, :) ! 3,3,nat,nat
  !! Dynamical matrix
  !
  axis = 3 
  ! 
  ! the call to set_ndnmbr is just a trick to get quickly
  ! a file label by exploiting an existing subroutine
  ! (if you look at the sub you will find that the original
  ! purpose was for pools and nodes)
  !
  CALL set_ndnmbr(0, iq_irr, 1, nqc_irr, filelab)
  tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dyn_q' // TRIM(filelab)
  ! The following function will check either or not the file is formatted in
  ! xml. If no file is found, an error is raised
  CALL check_is_xml_file(tempfile, is_xml_file)
  ! 
  IF (is_xml_file) THEN
    CALL read_dyn_mat_param(tempfile, ntyp, nat)
    ALLOCATE(m_loc(3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('dynmat', 'Error allocating m_loc', 1)
    ALLOCATE(dchi_dtau(3, 3, 3,nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('dynmat', 'Error allocating dchi_dtau', 1)
    CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
            celldm, at, bg, omega, atm, amass2, tau, ityp, &
            m_loc, nqs, lrigid, epsi_, zstar_, lraman, dchi_dtau)
    ALLOCATE(dyn(3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('dynmat', 'Error allocating dyn', 1)
    IF (ionode) THEN
      DO nt = 1, ntyp
        IF (amass(nt) <= 0.0d0) amass(nt) = amass2(nt)
      ENDDO
    ENDIF
    !
    IF (lrigid) THEN
      WRITE(stdout, '(5x,a)') 'Read dielectric tensor and effective charges'
      zstar = zstar_
      epsi = epsi_
      !ASR on effective charges
      DO i = 1, 3
        DO j = 1, 3
          sumz = 0.0d0
          DO na = 1,nat
            sumz = sumz + zstar(i, j, na)
          ENDDO
          DO na = 1,nat
            zstar(i, j, na) = zstar(i, j, na) - sumz / nat
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! If time-reversal is not included in the star of q, then double the nq to
    ! search from.
    IF (imq == 0) THEN
      mq = 2 * nq
    ELSE
      mq = nq
    ENDIF
    !
    ! First read the data and then sort them. This is because the order of 
    ! q in the star has changed (SP).
    ! 
    DO iq = 1, mq
      q(:, iq) = 0.0d0
      CALL read_dyn_mat(nat, iq, qout, dyn(:, :, :, :))
      q(:, iq) = qout(:)
      !
      DO na = 1,nat
        DO ipol = 1, 3
          DO jpol = 1, 3
            dynr(1, ipol, na, jpol, :) = REAL(dyn(ipol, jpol, na, :))
            dynr(2, ipol, na, jpol, :) = REAL(AIMAG(dyn(ipol, jpol, na, :)))
          ENDDO
        ENDDO
      ENDDO
      !
      ! Impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
      ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
      !
      IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
        WRITE(stdout, '(5x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
        IF (lpolar .AND. .NOT. lrigid) CALL errore('dynmat', &
          &'You set lpolar = .TRUE. but did not put epsil = true in the PH calculation at Gamma. ',1)
      ENDIF
      DO na = 1, nat
        DO ipol = 1, 3
          DO jpol = ipol, 3
            !
            IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6 ) THEN
              sumr(1, ipol, na, jpol) = SUM(dynr(1, ipol, na, jpol, :))
              sumr(2, ipol, na, jpol) = SUM(dynr(2, ipol, na, jpol, :))
            ENDIF
            !
            dynr(:, ipol, na, jpol, na) = dynr(:, ipol, na, jpol, na) - sumr(:, ipol, na, jpol)
            !
          ENDDO
        ENDDO
      ENDDO
      !
      ! Fill the two-indices dynamical matrix in cartesian coordinates
      ! the proper index in the complete list is iq_first+iq-1
      !
      DO na = 1, nat
        DO nb = 1, nat
          DO ipol = 1, 3
            DO jpol = 1, 3
              !
              mu = (na - 1) * 3 + ipol
              nu = (nb - 1) * 3 + jpol
              ! Only store the dyn of the first q in the star.
              IF (iq == 1) THEN
                dynq(mu, nu, iq_first) = DCMPLX(dynr(1, ipol, na, jpol, nb), dynr(2, ipol, na, jpol, nb))
              ENDIF
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO !  iq = 1, mq
    !  
  ELSE ! not a xml file
    OPEN(UNIT = iudyn, FILE = tempfile, STATUS = 'old', IOSTAT = ios)
    IF (ios /= 0) CALL errore('dynmat', 'opening file' // tempfile, ABS(ios))
    !
    !  read header and run some checks
    !
    READ(iudyn, '(a)') line
    READ(iudyn, '(a)') line
    READ(iudyn, *) ntyp_, nat_, ibrav_, celldm_
    ! 
    ! We stop testing celldm as it can be different between scf and nscf
    !IF (ntyp/=ntyp_.OR.nat/=nat_.OR.ibrav_/=ibrav.OR.abs ( &
    !   celldm_ (1) - celldm (1) )  > 1.0d-5) call errore ('dynmat', &
    !   'inconsistent data', 1)
    IF (ntyp /= ntyp_ .OR. nat /= nat_ .OR. ibrav_ /= ibrav) CALL errore ('dynmat', &
       'inconsistent data', 1)
    ! 
    !  skip reading of cell parameters here
    ! 
    IF (ibrav_ == 0) THEN
      DO i = 1, 4
        READ(iudyn, *) line
      ENDDO
    ENDIF
    DO nt = 1, ntyp
      READ(iudyn, *) i, atm, amass_
      IF (nt /=i .OR. ABS(amass_ - amass(nt))  > 1.0d-2) THEN
        WRITE(stdout,*) amass_, amass(nt)
        CALL errore('dynmat', 'inconsistent data', 1)
      ENDIF
    ENDDO
    DO na = 1, nat
      READ(iudyn, * ) i, ityp_, tau_
      IF (na /= i .OR. ityp_ /= ityp(na)) CALL errore('dynmat', &
          'inconsistent data (names)', 10 + na)
    ENDDO
    !
    ! Read dyn mat only for the first q in the star and then reconstruct the
    ! other using symmetry.
    ! 
    ! If time-reversal is not included in the star of q, then double the nq to
    ! search from.
    IF (imq == 0) THEN
      mq = 2 * nq
    ELSE
      mq = nq
    ENDIF
    !
    ! Read the dyn for the first q in the star. 
    ! For other q in the star, only read the value of the q for checking purposes.
    ! 
    DO iq = 1, mq
      !
      READ(iudyn, '(///a)') line
      READ(line(11:80), *) (q(i, iq), i = 1, 3)
      READ(iudyn, '(a)') line
      !
      DO na = 1, nat
        DO nb = 1, nat
          READ(iudyn, *) naa, nbb
          IF (na /= naa .OR. nb /= nbb) CALL errore('readmat', 'error reading file', nb)
          READ(iudyn, *)((dynr(1, i, na, j, nb), dynr(2, i, na, j, nb), j = 1, 3), i = 1, 3)
        ENDDO
      ENDDO
      !
      ! impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
      ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
      !
      IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
        WRITE(stdout, '(5x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
      ENDIF
      DO na = 1, nat
        DO ipol = 1,3
          DO jpol = ipol, 3
            !
            IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
              sumr(1, ipol, na, jpol) = SUM(dynr(1, ipol, na, jpol, :))
              sumr(2, ipol, na, jpol) = SUM(dynr(2, ipol, na, jpol, :))
            ENDIF
            !
            dynr(:, ipol, na, jpol, na) = dynr(:, ipol, na, jpol, na) - sumr(:, ipol, na, jpol)
            !
          ENDDO
        ENDDO
      ENDDO
      !
      !  fill the two-indices dynamical matrix in cartesian coordinates
      !  the proper index in the complete list is iq_first+iq-1
      !
      DO na = 1, nat
        DO nb = 1, nat
          DO ipol = 1, 3
            DO jpol = 1, 3
              !
              mu = (na - 1) * 3 + ipol
              nu = (nb - 1) * 3 + jpol
              ! Only store the dyn of the first q in the star.
              IF (iq == 1) THEN 
                dynq(mu, nu, iq_first) = DCMPLX(dynr(1, ipol, na, jpol, nb), dynr(2, ipol, na, jpol, nb))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    IF (ABS(q(1, 1)) < eps6 .AND. ABS(q(2, 1)) < eps6 .AND. ABS(q(3, 1)) < eps6) THEN
      !  read dielectric tensor and effective charges if present
      !  SP: Warning zstar is not properly bcast at the moment
      READ(iudyn, '(a)') line
      READ(iudyn, '(a)') line
      !
      IF (line(6:15) .EQ. 'Dielectric') THEN
        READ(iudyn, '(a)') line
        READ(iudyn, *) ((epsi(i, j), j = 1, 3), i = 1, 3)
        READ(iudyn, '(a)') line
        READ(iudyn, '(a)') line
        READ(iudyn, '(a)') line
        DO na = 1, nat
          READ(iudyn, '(a)') line
          READ(iudyn, *) ((zstar(i, j, na), j = 1, 3), i = 1, 3)
        ENDDO
        WRITE(stdout,'(5x,a)') 'Read dielectric tensor and effective charges'
        !
        !ASR on effective charges
        DO i = 1, 3
          DO j = 1, 3
            sumz = 0.0d0
            DO na = 1, nat
              sumz = sumz + zstar(i, j, na)
            ENDDO
            DO na = 1, nat
              zstar(i, j, na) = zstar(i, j, na) - sumz / nat
            ENDDO
          ENDDO
        ENDDO
        !
      ELSE 
        IF (lpolar) CALL errore('dynmat', &
           'You set lpolar = .TRUE. but did not put epsil = true in the PH calculation at Gamma. ', 1)
      ENDIF
    ENDIF
    CLOSE(iudyn)
  ENDIF ! col
  !
  ! Now check that the dyn file is consistent with the current EPW run (SP)
  ! SP: Be careful, here time-reversal is not actual time reversal but is due to 
  !     change in order and values of the q in the star between QE 4 and 5.
  !
  CALL cryst_to_cart(nqc1 * nqc2 * nqc3, q, at, -1)
  CALL cryst_to_cart(48, sxq, at, -1)
  !
  current_iq = iq_first
  !
  DO iq = 1, nq
    ! 
    found = .FALSE.
    DO jq = 1, mq
      DO m1 = -2, 2
        DO m2 = -2, 2
          DO m3 = -2, 2
            IF ((ABS(q(1, jq) - (sxq(1, iq) + m1)) < eps6 .AND. &
                 ABS(q(2, jq) - (sxq(2, iq) + m2)) < eps6 .AND. &
                 ABS(q(3, jq) - (sxq(3, iq) + m3)) < eps6 )) THEN
              found = .TRUE.
              EXIT ! exit loop
            ENDIF
          ENDDO
          IF (found) EXIT
        ENDDO
        IF (found) EXIT
      ENDDO
      IF (found) EXIT
    ENDDO
    current_iq = current_iq + 1
    IF (found .EQV. .FALSE.) THEN
      CALL errore('dynmat', 'wrong qpoint', 1)
    ENDIF
  ENDDO
  ! Transform back the sxq in Cartesian 
  CALL cryst_to_cart(48, sxq, bg, 1)
  !
  ! In case of reading of the IFC to impose the asr in real space
  ! We still call the above just to make the checks. The content of dynq 
  ! will be re-written just below and NOT read from the dyn from the /save folder
  IF (lifc) THEN
    !
    ! build the WS cell corresponding to the force constant grid
    atws(:, 1) = at(:, 1) * DBLE(nqc1)
    atws(:, 2) = at(:, 2) * DBLE(nqc2)
    atws(:, 3) = at(:, 3) * DBLE(nqc3)
    ! initialize WS r-vectors
    CALL wsinit(rws, nrwsx, nrws, atws)
    ! dynifc2blochc requires ifc
    CALL dynifc2blochc(nmodes, rws, nrws, q(:, 1), dynq_tmp)
    dynq(:, :, iq_first) = dynq_tmp
    WRITE(stdout, '(5x,a)') "Dyn mat calculated from ifcs"
    !
  ENDIF
  !
  ! Now construct the other dyn matrix for the q in the star using sym. 
  ! For this we use the gamma matrix. 
  ! 
  current_iq = iq_first 
  ! 
  DO iq = 1, nq 
    !
    xq = sxq(:, iq)
    nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
    sym_sgq(:) = 0
    DO jsym = 1, nsym
      IF (isq(jsym) == iq ) then
        nsq = nsq + 1
        sym_sgq(nsq) = jsym
      ENDIF
    ENDDO
    !
    ! SP: We now need to select one symmetry among the small group of q (i.e. that respect 
    !     Sq0+G=q ) that has G=0. There should always be such symmetry. 
    !     We enforce this for later easiness. 
    ! 
    aq = sxq(:, 1) ! This is xq0
    saq = xq
    call cryst_to_cart(1, aq, at, - 1)
    CALL cryst_to_cart(1, saq, at, -1)
    ! Initialize isym
    isym = 1
    DO jsym = 1, nsq
      ism1 = invs(sym_sgq(jsym))
      raq = 0.d0
      DO ipol = 1, 3
        DO jpol = 1, 3
          raq(ipol) = raq(ipol) + s(ipol, jpol, ism1) * aq(jpol)
        ENDDO
      ENDDO
      nog = eqvect_strict(raq, saq, eps6)
      IF (nog) THEN ! This is the symmetry such that Sq=q
        isym = sym_sgq(jsym)
        EXIT
      ENDIF
      ! If we enter into that loop it means that we have not found 
      ! such symmetry within the small group of Q. 
      IF (jsym == nsq) THEN
        CALL errore('dynmat ', 'No sym. such that Sxq0=iq was found in the sgq !', 1)
      ENDIF
    ENDDO
    !
    !  -----------------------------------------------------------------------
    !  the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)   
    !  -----------------------------------------------------------------------
    ! 
    ism1 = invs(isym)
    !
    !  the symmetry matrix in cartesian coordinates 
    !  (so that we avoid going back and forth with the dynmat)  
    !  note the presence of both at and bg in the transform!
    !
    scart = DBLE(s(:, :, ism1))
    scart = MATMUL(MATMUL(bg, scart), TRANSPOSE(at))
    ! 
    gamma = czero
    DO na = 1, nat
      !
      ! the corresponding of na in {S|v}
      sna = irt(isym, na)
      !
      ! cfac = exp[iSq*(tau_K - {S|v} tau_k)]   (Maradudin&Vosko RMP Eq. 2.33)
      ! [v can be ignored since it cancels out, see endnotes. xq is really Sq]
      ! rtau(:,isym,na) = s(:,:,invs(isym)) * tau(:, na) - tau(:,irt(isym,na))) (cartesian)
      !
      arg = twopi * DOT_PRODUCT(xq, rtau (:, isym, na))
      cfac = DCMPLX(COS(arg),-SIN(arg))
      !
      !  the submatrix (sna,na) contains the rotation scart 
      !
      gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart
      !
    ENDDO  
    !
    !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5) 
    ! 
    CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
           nmodes, dynq(:, :, iq_first), nmodes, czero, dynq_tmp, nmodes)
    CALL zgemm ('n', 'c', nmodes, nmodes, nmodes, cone, dynq_tmp, &
           nmodes, gamma, nmodes, czero, dynq(:, :,current_iq), nmodes)
    !
    DO nu = 1, nmodes
      DO mu = 1, nmodes
        IF (mu /= nu .AND. ABS(dynq(mu, nu, current_iq)) > eps6 ) CALL errore &
          ('rotate_eigenm', 'problem with rotated eigenmodes', 0)
        ENDDO
    ENDDO
    !
    ! DBSP-----------------------------------------------
    !  a simple check on the frequencies
    !
    IF (iverbosity == 1) THEN
      DO na = 1, nat
        DO nb = 1, nat
          massfac = 1.d0 / SQRT(amass(ityp(na)) * amass(ityp(nb)))
          dynq_tmp(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
          dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, current_iq) * massfac
        END DO
      END DO
      !
      DO jmode = 1, nmodes
        DO imode = 1, jmode
           dynp(imode + (jmode - 1) * jmode / 2) = &
               (dynq_tmp(imode, jmode) + CONJG(dynq_tmp(jmode, imode))) / 2.d0
        ENDDO
      ENDDO
      !
      CALL ZHPEVX('V', 'A', 'U', nmodes, dynp, 0.0, 0.0, 0, 0, -1.0, neig, w1, cz1, nmodes, cwork, &
                 rwork, iwork, ifail, info)
      ! 
      DO nu = 1, nmodes
        IF (w1(nu) > 0.d0) THEN
          wtmp(nu) =  SQRT(ABS(w1(nu)))
        ELSE
          wtmp(nu) = -SQRT(ABS(w1(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the current q in the star (cm^-1)")')
      WRITE(stdout, '(6(2x,f10.5))' ) (wtmp(nu) * rydcm1, nu = 1, nmodes)
    ENDIF
    !END --------------------------------------------------
    current_iq = current_iq + 1
    ! 
    ! SP Repeat the same but for minus_q one
    IF (imq == 0) then
      !       
      xq = -sxq(:, iq)
      ism1 = invs(isym)
      scart = DBLE(s(:, :, ism1))
      scart = MATMUL(MATMUL(bg, scart), TRANSPOSE(at))
      ! 
      gamma = czero
      DO na = 1, nat
        ! 
        sna = irt(isym, na)
        arg = twopi * DOT_PRODUCT(xq, rtau(:, isym, na))
        cfac = dcmplx(COS(arg), -SIN(arg))
        gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart
        !
      ENDDO
      !
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma  , &
             nmodes, CONJG(dynq(:, :, iq_first)) , nmodes, czero , dynq_tmp, nmodes)
      CALL ZGEMM('n', 'c', nmodes, nmodes, nmodes, cone, dynq_tmp, &
             nmodes, gamma, nmodes, czero , dynq(:, :, current_iq), nmodes)
      !
      DO nu = 1, nmodes
        DO mu = 1, nmodes
          IF (mu /= nu .AND. ABS(dynq(mu, nu, current_iq)) > eps6 ) CALL errore &
            ('rotate_eigenm', 'problem with rotated eigenmodes', 0)
        ENDDO
      ENDDO
      current_iq = current_iq + 1
    ENDIF
  ENDDO ! iq
  ! 
  !---------------------------------------------------------------------------------
  END SUBROUTINE dynmat_asr
  !---------------------------------------------------------------------------------

