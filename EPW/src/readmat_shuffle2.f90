  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine readmat_shuffle2 ( iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq,&
                                invs, s, irt, rtau)
  !-----------------------------------------------------------------------
  !!
  !! read dynamical matrix for the q points  
  !! iq_first, iq_first+1, ... iq_first+nq-1
  !!
  !-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  use io_files,         ONLY : prefix 
  USE cell_base,        ONLY : ibrav, celldm, omega, at, bg
  USE ions_base,        ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE elph2,            ONLY : dynq, sumr, zstar, epsi
  USE symm_base,        ONLY : nsym
  USE epwcom,           ONLY : dvscf_dir, lpolar, lifc
  USE modes,            ONLY : nmodes
  USE control_flags,    ONLY : iverbosity
  USE phcom,            ONLY : nq1, nq2, nq3
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE io_dyn_mat2,      ONLY : read_dyn_mat_param, read_dyn_mat_header,&
                               read_dyn_mat
  USE io_global,        ONLY : ionode, stdout
  USE constants_epw,    ONLY : cone, czero, twopi, rydcm1
  USE io_epw,           ONLY : iudyn
  !
  implicit none
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
  INTEGER, INTENT(in) :: isq (48)
  !! Index of q in the star for a given symmetry.
  INTEGER, INTENT(in) :: invs (48)
  !! Matrix of inversed symmetry 
  INTEGER, INTENT(in) :: s(3, 3, 48)
  !! Set of symmetry operations
  INTEGER, INTENT(in) :: irt(48,nat)
  !! For each atoms give the rotated atoms
  !
  REAL(kind=DP), INTENT(inout) :: sxq(3,48)
  !! Symmetry matrix 
  REAL(kind=DP), INTENT(inout) :: rtau(3, 48, nat)
  !! the relative position of the rotated atom to the original one
  !
  ! Work variables 
  INTEGER :: isym, m1,m2,m3, sna, jsym, &
             neig, iwork( 5*nmodes ),  info, ifail( nmodes), nsq, sym_sgq(48)
  INTEGER ::  current_iq, mq
  !  the q index in the dynq matrix
  INTEGER :: ism1, imode, jmode
  ! 
  LOGICAL :: eqvect_strict
  !! True if we are using time reversal
  !
  REAL(kind=DP) ::  arg, rwork( 7*nmodes ), aq(3), saq(3), raq(3)
  REAL(kind=DP) :: xq(3)
  !!  the rotated q vector
  REAL(kind=DP) :: scart(3,3), massfac, w1( nmodes ),wtmp( nmodes )
  !  the q vectors in the star
  COMPLEX(kind=DP) :: cfac, cwork( 2*nmodes )
  COMPLEX(kind=DP) :: gamma(nmodes, nmodes), dynq_tmp(nmodes, nmodes), &
                      cz1(nmodes, nmodes), dynp( nmodes*(nmodes+1)/2 )
  !  the Gamma matrix for the symmetry operation on the dyn mat
  !
  ! output
  !
  ! dynq (nmode,nmodes,nqc) (in elph2.mod)
  !
  ! local variable
  !
  CHARACTER(len=3)    :: atm, filelab
  CHARACTER(len=80)   :: line
  CHARACTER(len=256)  :: tempfile
  LOGICAL             :: found, lrigid, lraman, nog
  INTEGER             :: ntyp_, nat_, ibrav_, ityp_, ios, iq, jq,  &
                         nt, na, nb, naa, nbb, nu, mu, i, j, ipol,jpol
  INTEGER, parameter  :: ntypx = 10
  INTEGER             :: nqs, axis
  INTEGER, PARAMETER  :: nrwsx=200
  INTEGER             :: nrws
  REAL(kind=DP)       :: eps, celldm_ (6), amass_, tau_ (3), q(3,nq1*nq2*nq3), &
                         dynr (2, 3, nat, 3, nat), sumz
  REAL(DP), allocatable :: m_loc(:,:), dchi_dtau(:,:,:,:)
  REAL(DP)              :: qout(3), amass2(ntypx)
  REAL(kind=DP)         :: rws(0:3,nrwsx), atws(3,3)
  COMPLEX(KIND=DP),allocatable   :: dyn(:,:,:,:) ! 3,3,nat,nat
  !
  axis = 3 
  !
  eps = 1.d-5
  ! 
!DBSP
  ! SP: If noncolin, the dynamical matrix are printed in xml format by QE
  IF (noncolin) THEN
    CALL set_ndnmbr ( 0, iq_irr, 1, nqc_irr, filelab)
    tempfile = trim(dvscf_dir) // trim(prefix) // '.dyn_q' // filelab
    CALL read_dyn_mat_param(tempfile,ntyp,nat)
    ALLOCATE (m_loc(3,nat))
    ALLOCATE (dchi_dtau(3,3,3,nat) )
    CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
            celldm, at, bg, omega, atm, amass2, tau, ityp, &
            m_loc, nqs, lrigid, epsi, zstar, lraman, dchi_dtau)
    ALLOCATE (dyn(3,3,nat,nat) )
    IF (ionode) THEN
       DO nt=1, ntyp
          IF (amass(nt) <= 0.0d0) amass(nt)=amass2(nt)
       END DO
    END IF
    !
    IF (lrigid) THEN
       WRITE (6,'(8x,a)') 'Read dielectric tensor and effective charges'
       !ASR on effective charges
        DO i=1,3
           DO j=1,3
              sumz=0.0d0
              DO na=1,nat
                 sumz = sumz + zstar(i,j,na)
              ENDDO
              DO na=1,nat
                 zstar(i,j,na) = zstar(i,j,na) - sumz/nat
              ENDDO
           ENDDO
        ENDDO
    ENDIF
    !
    ! If time-reversal is not included in the star of q, then double the nq to
    ! search from.
    IF (imq.eq.0) then
      mq = 2*nq
    ELSE
      mq = nq
    ENDIF
    !
    ! First read the data and then sort them. This is because the order of 
    ! q in the star has changed (SP).
    ! 
    DO iq = 1, mq
      q(:,iq) = 0.0d0
      CALL read_dyn_mat(nat,iq,qout,dyn(:,:,:,:))
      q(:,iq) = qout(:)
      !
      DO na = 1,nat
       DO ipol = 1,3
        DO jpol = 1,3
          dynr(1,ipol,na,jpol,:) = REAL(dyn(ipol,jpol,na,:))
          dynr(2,ipol,na,jpol,:) = REAL(AIMAG(dyn(ipol,jpol,na,:)))
        ENDDO
       ENDDO
      ENDDO
      !
      ! Impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
      ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
      !
      IF ( abs(q(1,iq)).lt.eps .and. abs(q(2,iq)).lt.eps .and. abs(q(3,iq)).lt.eps ) THEN
        WRITE(6,'(8x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
      ENDIF
      DO na = 1,nat
       DO ipol = 1,3
        DO jpol = ipol,3
         !
          IF ( abs(q(1,iq)).lt.eps .and. abs(q(2,iq)).lt.eps .and. abs(q(3,iq)).lt.eps ) then
             IF ( .not. allocated(sumr) ) allocate ( sumr(2,3,nat,3) )
             sumr(1,ipol,na,jpol) = sum ( dynr (1,ipol,na,jpol,:) )
             sumr(2,ipol,na,jpol) = sum ( dynr (2,ipol,na,jpol,:) )
          ENDIF
          !
          dynr (:,ipol,na,jpol,na) = dynr (:,ipol,na,jpol,na) - sumr(:,ipol,na,jpol)
          !
        ENDDO
       ENDDO
      ENDDO
      !
      !  fill the two-indices dynamical matrix in cartesian coordinates
      !  the proper index in the complete list is iq_first+iq-1
      !
      DO na = 1,nat
        DO nb = 1,nat
          DO ipol = 1,3
           DO jpol = 1,3
             !
             mu = (na-1)*3+ipol
             nu = (nb-1)*3+jpol
             ! Only store the dyn of the first q in the star.
             if (iq == 1) then
               dynq ( mu, nu, iq_first) = &
                 dcmplx ( dynr (1,ipol,na,jpol,nb), dynr (2,ipol,na,jpol,nb) )
             endif
             !
           ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO !  iq = 1, mq
    !  
  ELSE ! noncolin
!END
    !
    !  the call to set_ndnmbr is just a trick to get quickly
    !  a file label by exploiting an existing subroutine
    !  (if you look at the sub you will find that the original
    !  purpose was for pools and nodes)
    !
    CALL set_ndnmbr ( 0, iq_irr, 1, nqc_irr, filelab)
    tempfile = trim(dvscf_dir) // trim(prefix) // '.dyn_q' // filelab
    !
    open (unit = iudyn, file = tempfile, status = 'old', iostat = ios)
    IF (ios /=0)  call errore ('readmat_shuffle2', 'opening file'//tempfile, abs (ios) )
    !
    !  read header and run some checks
    !
    read (iudyn, '(a)') line
    read (iudyn, '(a)') line
    read (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
    IF (ntyp.ne.ntyp_.or.nat.ne.nat_.or.ibrav_.ne.ibrav.or.abs ( &
       celldm_ (1) - celldm (1) ) .gt.1.0d-5) call errore ('readmat2', &
       'inconsistent data', 1)
    ! 
    !  skip reading of cell parameters here
    ! 
    IF (ibrav_ .eq. 0) then
       DO i = 1,4
          read (iudyn, * ) line
       ENDDO
    ENDIF
    DO nt = 1, ntyp
       read (iudyn, * ) i, atm, amass_
       IF (nt.ne.i.or.abs (amass_ - amass (nt) ) .gt.1.0d-2) then
       write (6,*) amass_, amass(nt)
       call errore ('readmat', 'inconsistent data', 0)
    endif
    ENDDO
    DO na = 1, nat
       read (iudyn, * ) i, ityp_, tau_
       IF (na.ne.i.or.ityp_.ne.ityp (na) ) call errore ('readmat2', &
          'inconsistent data (names)', 10 + na)
    ENDDO
    !
    ! Read dyn mat only for the first q in the star and then reconstruct the
    ! other using symmetry.
    ! 
    ! If time-reversal is not included in the star of q, then double the nq to
    ! search from.
    IF (imq.eq.0) then
      mq = 2*nq
    ELSE
      mq = nq
    ENDIF
    !
    ! Read the dyn for the first q in the star. 
    ! For other q in the star, only read the value of the q for checking purposes.
    ! 
    DO iq = 1, mq
      !
      !
      read (iudyn, '(///a)') line
      read (line (11:80), * ) (q (i,iq), i = 1, 3)
      !
      read (iudyn, '(a)') line
      !
      !
      DO na = 1, nat
         DO nb = 1, nat
            read (iudyn, * ) naa, nbb
            IF (na.ne.naa.or.nb.ne.nbb) call errore &
                 ('readmat', 'error reading file', nb)
            read (iudyn, * ) &
                ((dynr (1,i,na,j,nb), dynr (2,i,na,j,nb), j = 1, 3), i = 1, 3)
         ENDDO
      ENDDO
      !
      ! impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
      ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
      !
      IF ( abs(q(1,iq)).lt.eps .and. abs(q(2,iq)).lt.eps .and. abs(q(3,iq)).lt.eps ) then
        WRITE(6,'(8x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
      ENDIF
      DO na = 1,nat
       DO ipol = 1,3
        DO jpol = ipol,3
          !
          IF ( abs(q(1,iq)).lt.eps .and. abs(q(2,iq)).lt.eps .and. abs(q(3,iq)).lt.eps ) then
             IF ( .not. allocated(sumr) ) allocate ( sumr(2,3,nat,3) )
             sumr(1,ipol,na,jpol) = sum ( dynr (1,ipol,na,jpol,:) )
             sumr(2,ipol,na,jpol) = sum ( dynr (2,ipol,na,jpol,:) )
          ENDIF
          !
          dynr (:,ipol,na,jpol,na) = dynr (:,ipol,na,jpol,na) - sumr(:,ipol,na,jpol)
          !
        ENDDO
       ENDDO
      ENDDO
      !
      !  fill the two-indices dynamical matrix in cartesian coordinates
      !  the proper index in the complete list is iq_first+iq-1
      !
      DO na = 1,nat
       DO nb = 1,nat
          DO ipol = 1,3
           DO jpol = 1,3
             !
             mu = (na-1)*3+ipol
             nu = (nb-1)*3+jpol
             ! Only store the dyn of the first q in the star.
             if (iq == 1) then 
               dynq ( mu, nu, iq_first) = &
                 dcmplx ( dynr (1,ipol,na,jpol,nb), dynr (2,ipol,na,jpol,nb) )
             endif
             !
           ENDDO
          ENDDO
       ENDDO
      ENDDO
      !
    ENDDO
    !
    IF ( abs(q(1,1)).lt.eps .and. abs(q(2,1)).lt.eps .and. abs(q(3,1)).lt.eps ) THEN
      !  read dielectric tensor and effective charges if present
      !  SP: Warning zstar is not properly bcast at the moment
      read (iudyn,'(a)') line
      read (iudyn,'(a)') line
      !
      IF (lpolar) THEN
        IF (line(6:15).EQ.'Dielectric') THEN
          read (iudyn,'(a)') line
          read (iudyn,*) ((epsi(i,j), j=1,3), i=1,3)
          read (iudyn,'(a)') line
          read (iudyn,'(a)') line
          read (iudyn,'(a)') line
          DO na = 1,nat
             read (iudyn,'(a)') line
             read (iudyn,*) ((zstar(i,j,na), j=1,3), i=1,3)
          ENDDO
          WRITE (6,'(8x,a)') 'Read dielectric tensor and effective charges'
          !
          !ASR on effective charges
          DO i=1,3
             DO j=1,3
                sumz=0.0d0
                DO na=1,nat
                   sumz = sumz + zstar(i,j,na)
                ENDDO
                DO na=1,nat
                   zstar(i,j,na) = zstar(i,j,na) - sumz/nat
                ENDDO
             ENDDO
          ENDDO
          !
        ELSE 
          call errore('readmat_shuffle2','You set lpolar = .true. but did not put epsil = true in the PH calculation at Gamma. ',1)
        ENDIF
        !
      ENDIF
      !   
    ENDIF
    !
    CLOSE(iudyn)
  ENDIF ! col
  !
  ! Now check that the dyn file is consistent with the current EPW run (SP)
  ! SP: Be careful, here time-reversal is not actual time reversal but is due to 
  !     change in order and values of the q in the star between QE 4 and 5.
  !
  CALL cryst_to_cart (nq1*nq2*nq3, q, at, -1)
  CALL cryst_to_cart (48, sxq, at, -1)
  !
  current_iq = iq_first
  !
  DO iq = 1, nq
    ! 
    found = .false.
    DO jq = 1, mq
      DO m1=-2,2
        DO m2=-2,2
          DO m3=-2,2
            IF ((abs(q(1,jq)-(sxq(1,iq)+m1)).lt.eps .and. &
                 abs(q(2,jq)-(sxq(2,iq)+m2)).lt.eps .and. &
                 abs(q(3,jq)-(sxq(3,iq)+m3)).lt.eps )) THEN
                 found = .true.
                 exit ! exit loop
            ENDIF
          ENDDO
          if (found) exit
        ENDDO
        if (found) exit
      ENDDO
      if (found) exit
    ENDDO
    current_iq = current_iq + 1
    IF (found .eqv. .false.) THEN
      call errore ('readmat_shuffle2', 'wrong qpoint', 1)
    ENDIF
  ENDDO
  ! Transform back the sxq in Cartesian 
  CALL cryst_to_cart (48, sxq, bg, 1)
  !
  ! In case of reading of the IFC to impose the asr in real space
  ! We still call the above just to make the checks. The content of dynq 
  ! will be re-written just below and NOT read from the dyn from the /save folder
  IF (lifc) THEN
     !
     ! build the WS cell corresponding to the force constant grid
     atws(:,1) = at(:,1)*DBLE(nq1)
     atws(:,2) = at(:,2)*DBLE(nq2)
     atws(:,3) = at(:,3)*DBLE(nq3)
     ! initialize WS r-vectors
     CALL wsinit(rws,nrwsx,nrws,atws)
     CALL dynifc2blochc (nmodes, rws, nrws, q(:,1), dynq_tmp)
     dynq(:,:,iq_first)=dynq_tmp
     write(stdout,*) " Dyn mat calculated from ifcs"
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
    xq = sxq(:,iq)
    !
    nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
    !
    sym_sgq(:) = 0
    DO jsym = 1, nsym
      IF ( isq(jsym) .eq. iq ) then
        nsq = nsq + 1
        sym_sgq(nsq) = jsym
      ENDIF
    ENDDO
    !
    ! SP: We now need to select one symmetry among the small group of q (i.e. that respect 
    !     Sq0+G=q ) that has G=0. There should always be such symmetry. 
    !     We enforce this for later easiness. 
    ! 
    aq = sxq(:,1) ! This is xq0
    saq = xq
    call cryst_to_cart (1, aq, at, - 1)
    CALL cryst_to_cart (1, saq, at, -1)
    ! Initialize isym
    isym = 1
    DO jsym=1, nsq
      ism1 = invs (sym_sgq(jsym))
      raq = 0.d0
      DO ipol = 1, 3
         DO jpol = 1, 3
            raq (ipol) = raq (ipol) + s (ipol, jpol, ism1) * aq (jpol)
         ENDDO
      ENDDO
      nog = eqvect_strict( raq, saq)
      IF (nog) THEN ! This is the symmetry such that Sq=q
        isym = sym_sgq(jsym)
        EXIT
      ENDIF
      ! If we enter into that loop it means that we have not found 
      ! such symmetry within the small group of Q. 
      IF (jsym == nsq) THEN
        CALL errore( 'readmat_shuffle2 ', 'No sym. such that Sxq0=iq was found in the sgq !', 1 )
      ENDIF
    ENDDO
    !
    !  -----------------------------------------------------------------------
    !  the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)   
    !  -----------------------------------------------------------------------
    ! 
    ism1 = invs (isym)
    !
    !  the symmetry matrix in cartesian coordinates 
    !  (so that we avoid going back and forth with the dynmat)  
    !  note the presence of both at and bg in the transform!
    !
    scart = dble ( s ( :, :, ism1) )
    scart = matmul ( matmul ( bg, scart), transpose (at) )
    ! 
    gamma = czero
    DO na = 1, nat
      !
      ! the corresponding of na in {S|v}
      sna = irt (isym, na)
      !
      ! cfac = exp[iSq*(tau_K - {S|v} tau_k)]   (Maradudin&Vosko RMP Eq. 2.33)
      ! [v can be ignored since it cancels out, see endnotes. xq is really Sq]
      ! rtau(:,isym,na) = s(:,:,invs(isym)) * tau(:, na) - tau(:,irt(isym,na))) (cartesian)
      !
      arg = twopi * dot_product (xq, rtau (:, isym, na))
      cfac = dcmplx (cos(arg),-sin(arg))
      !
      !  the submatrix (sna,na) contains the rotation scart 
      !
      gamma ( 3*(sna-1)+1:3*sna, 3*(na-1)+1:3*na ) = cfac * scart
      !
    ENDDO  
    !
    !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5) 
    ! 
    CALL zgemm ('n', 'n', nmodes, nmodes, nmodes , cone  , gamma  , &
           nmodes, dynq ( :, :, iq_first) , nmodes, czero , dynq_tmp, nmodes)
    CALL zgemm ('n', 'c', nmodes, nmodes, nmodes , cone  , dynq_tmp, &
           nmodes, gamma, nmodes, czero , dynq (:, :,current_iq), nmodes)
    !
    DO nu = 1, nmodes
      DO mu = 1, nmodes
      IF ( mu.ne.nu .and. abs(dynq(mu,nu,current_iq)).gt.eps ) call errore &
        ('rotate_eigenm','problem with rotated eigenmodes',0)
      ENDDO
    ENDDO
    !
    ! DBSP-----------------------------------------------
    !  a simple check on the frequencies
    !
    IF (iverbosity.eq.1) THEN
      DO na = 1, nat
       DO nb = 1, nat
         massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
         dynq_tmp(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
         dynq(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, current_iq) * massfac
       END DO
      END DO
      !
      DO jmode = 1, nmodes
       DO imode = 1, jmode
          dynp (imode + (jmode - 1) * jmode/2 ) = &
          ( dynq_tmp ( imode, jmode) + conjg ( dynq_tmp ( jmode, imode) ) ) / 2.d0
       ENDDO
      ENDDO
      !
      CALL zhpevx ('V', 'A', 'U', nmodes, dynp , 0.0, 0.0, &
                 0, 0,-1.0, neig, w1, cz1, nmodes, cwork, &
                 rwork, iwork, ifail, info)


      DO nu = 1, nmodes
        IF ( w1 (nu) .gt. 0.d0 ) then
           wtmp(nu) =  sqrt(abs( w1 (nu) ))
        ELSE
           wtmp(nu) = -sqrt(abs( w1 (nu) ))
        ENDIF
      ENDDO
      WRITE ( stdout, '(5x,"Frequencies of the matrix for the current q in the star")')
      WRITE ( stdout, '(6(2x,f10.5))' ) (wtmp(nu)*rydcm1, nu = 1, nmodes)
    ENDIF
    !END --------------------------------------------------
    current_iq = current_iq + 1
    ! 
    ! SP Repeat the same but for minus_q one
    IF (imq.eq.0) then
      !       
      xq = -sxq(:,iq)
      !saq = xq
      !CALL cryst_to_cart (1, saq, at, -1)
      !
      ism1 = invs (isym)
      !
      scart = dble ( s ( :, :, ism1) )
      scart = matmul ( matmul ( bg, scart), transpose (at) )
      ! 
      gamma = czero
      DO na = 1, nat
        !
        sna = irt (isym, na)
        !
        arg = twopi * dot_product (xq, rtau (:, isym, na))
        cfac = dcmplx (cos(arg),-sin(arg))
        !
        gamma ( 3*(sna-1)+1:3*sna, 3*(na-1)+1:3*na ) = cfac * scart
        !
      ENDDO
      !
      CALL zgemm ('n', 'n', nmodes, nmodes, nmodes , cone  , gamma  , &
             nmodes, CONJG(dynq ( :, :, iq_first)) , nmodes, czero , dynq_tmp, nmodes)
             !nmodes, CONJG(TRANSPOSE(dynq ( :, :, iq_first))) , nmodes, czero , dynq_tmp, nmodes)
      CALL zgemm ('n', 'c', nmodes, nmodes, nmodes , cone  , dynq_tmp, &
             nmodes, gamma, nmodes, czero , dynq (:, :,current_iq), nmodes)
      !
      DO nu = 1, nmodes
        DO mu = 1, nmodes
        IF ( mu.ne.nu .and. abs(dynq(mu,nu,current_iq)).gt.eps ) call errore &
          ('rotate_eigenm','problem with rotated eigenmodes',0)
        ENDDO
      ENDDO
      !
      current_iq = current_iq + 1
    ENDIF
    !
  ENDDO ! iq
  !
  end subroutine readmat_shuffle2
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE read_ifc
  !---------------------------------------------------------------------------------
  !!
  !! Read IFC in real space from the file generated by q2r. 
  !! Adapted from PH/matdyn.x by C. Verdi and S. Ponce
  !! 
  !
  USE kinds,     ONLY : DP
  USE elph2,     ONLY : ifc, zstar, epsi
  USE epwcom,    ONLY : asr_typ, dvscf_dir
  USE ions_base, ONLY : nat
  USE cell_base, ONLY : ibrav, omega, at, bg, celldm
  USE phcom,     ONLY : nq1, nq2, nq3
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iunifc
  USE constants_epw, ONLY :  czero
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE io_dyn_mat2,      ONLY : read_dyn_mat_param, read_dyn_mat_header,&
                               read_dyn_mat, read_ifc_xml, read_ifc_param
#if defined(__NAG)
  USE f90_unix_io,ONLY : flush
#endif
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
  !
  implicit none
  !
  LOGICAL             :: exst, lpolar_, has_zstar
  CHARACTER (len=80)  :: line
  CHARACTER(len=256)  :: tempfile
  INTEGER             :: ipol, ios, i, j, m1,m2,m3, na,nb, &
                         idum, ibid, jbid, nabid, nbbid, m1bid, m2bid, m3bid, &
                         ntyp_, nat_, ibrav_, ityp_(nat), nqs
  INTEGER, parameter  :: ntypx = 10
  REAL(kind=DP)       :: tau_(3,nat), alat, amass2(ntypx)
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)
  REAL(DP), ALLOCATABLE :: m_loc(:,:)
  !! he magnetic moments of each atom
  COMPLEX(kind=DP)    :: sumasr
  !
  WRITE(stdout,'(/5x,"Reading interatomic force constants"/)')
  CALL FLUSH(6)
  ! 
  ! This is important in restart mode as zstar etc has not been allocated
  IF (.NOT. ALLOCATED (zstar) ) ALLOCATE( zstar(3,3,nat) )
  IF (.NOT. ALLOCATED (epsi) ) ALLOCATE( epsi(3,3) )
  IF (.not. ALLOCATED (ifc)) ALLOCATE ( ifc ( nq1, nq2, nq3, 3, 3, nat, nat ) )
  zstar=0.d0
  epsi=0.d0

  IF (mpime.eq.ionode_id) THEN
    IF (noncolin) THEN
      !
      tempfile = trim(dvscf_dir) // 'ifc.q2r'
      CALL read_dyn_mat_param(tempfile,ntyp_,nat_)
      ALLOCATE (m_loc(3,nat_))
      ALLOCATE (atm(ntyp_))
      CALL read_dyn_mat_header(ntyp_, nat_, ibrav, nspin_mag, &
               celldm, at, bg, omega, atm, amass2, &
               tau_, ityp_,  m_loc, nqs, has_zstar, epsi, zstar )
      alat=celldm(1)
      call volume(alat,at(1,1),at(1,2),at(1,3),omega)
      CALL read_ifc_param(nq1,nq2,nq3)
      CALL read_ifc_xml(nq1,nq2,nq3,nat_,ifc)
      ! 
    ELSE
      !
      tempfile = trim(dvscf_dir) // 'ifc.q2r' 
      OPEN(unit=iunifc,file=tempfile,status='old',iostat=ios)
      IF (ios /= 0) call errore ('read_ifc', 'error opening ifc.q2r',iunifc)
      !
      !  read real-space interatomic force constants
      !
      READ(iunifc,'(3i4)') ntyp_ , nat_ , ibrav_
      IF (ibrav_ .eq. 0) then
         DO i = 1,4
            read (iunifc, * ) line
         ENDDO
      ENDIF
      DO i=1,ntyp_
         READ(iunifc,'(a)') line
      ENDDO
      DO na=1,nat
         READ(iunifc,*) idum, idum, (tau_(j,na),j=1,3)
              ENDDO
      READ(iunifc,*) lpolar_
      !
      IF (lpolar_) THEN
         READ (iunifc,*) ((epsi(i,j), j=1,3), i=1,3)
         DO na = 1, nat
            READ (iunifc,*) idum
            READ (iunifc,*) ((zstar(i,j,na), j=1,3), i=1,3)
         ENDDO
         WRITE(stdout,*) " Read Z*, epsilon"
      ENDIF
      !
      READ (iunifc,*) idum
      !
      ifc = 0.d0
      DO i=1,3
         DO j=1,3
            DO na=1,nat
               DO nb=1,nat
                  READ (iunifc,*) ibid, jbid, nabid, nbbid
                  IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                     na.NE.nabid .OR. nb.NE.nbbid)                      &
                     CALL errore  ('read_epw','error in reading ifc',1)
                  READ (iunifc,*) (((m1bid, m2bid, m3bid,        &
                              ifc(m1,m2,m3,i,j,na,nb),                  &
                               m1=1,nq1),m2=1,nq2),m3=1,nq3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
    ENDIF ! noncol
  ENDIF
  !
  ! It has to be casted like this because mpi cannot cast 7 indices
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           DO nb=1,nat
              CALL mp_bcast (ifc(:,:,:,i,j,na,nb), ionode_id, inter_pool_comm)
              CALL mp_bcast (ifc(:,:,:,i,j,na,nb), root_pool, intra_pool_comm)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  CALL mp_bcast (zstar, ionode_id, inter_pool_comm)
  CALL mp_bcast (zstar, root_pool, intra_pool_comm)
  CALL mp_bcast (epsi, ionode_id, inter_pool_comm)
  CALL mp_bcast (epsi, root_pool, intra_pool_comm)
  CALL mp_bcast (tau_, ionode_id, inter_pool_comm)
  CALL mp_bcast (tau_, root_pool, intra_pool_comm)
  CALL mp_bcast (ibrav_, ionode_id, inter_pool_comm)
  CALL mp_bcast (ibrav_, root_pool, intra_pool_comm)

  !
  write(stdout,*) ' IFC last ', ifc(nq1,nq2,nq3,3,3,nat,nat)
  !
  CALL set_asr2 (asr_typ, nq1, nq2, nq3, ifc, zstar, &
             nat, ibrav_, tau_)
  !
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
    CLOSE(iunifc)
  ENDIF
  !
  WRITE(stdout,'(/5x,"Finished reading ifcs"/)')
 write(stdout,*) " IFC ", ifc(1,1,1,1,1,1,1), ifc(nq1,nq2,nq3,3,3,nat,nat)
  !
  END SUBROUTINE read_ifc
!-------------------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE set_asr2 (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !!
  !! Set the acoustic sum rule. 
  !! Taken directly from PHonon/PH/q2trans.f90
  !! It would be better to take the set_asr for /Modules/. 
  !! However they are different (frc) and to be consitent with q2r.x we take this one.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), intent(in) :: tau(3,nat)
  REAL(DP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(DP) :: zeu_new(3,3,nat)
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     real(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements
  ! (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(DP) :: scal,norm2, sum
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if
  ! asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if
  ! axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
               do nb=1,nat
                  do n1=1,nr1
                     do n2=1,nr2
                        do n3=1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        end do
                     end do
                  end do
               end do
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
            end do
         end do
      end do
    write(stdout,*) " Imposed simple ASR"
      !
      return
      !
   end if

  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'asr: rotational axis may be wrong'
     endif
     write(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)')
  !  k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry
                       ! constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec(:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)')
  !  l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)')
  !  k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
   write(stdout,*) " Imposed crystal ASR"
  !
  return
end subroutine set_asr2
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered
  ! as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual
  ! way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered
  ! as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual
  ! way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3
!

