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
  USE epwcom,           ONLY : dvscf_dir, lpolar
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
  integer :: isym, invs (48), m1,m2,m3, s(3, 3, 48), sna, irt(48,nat), jsym, &
             neig, iwork( 5*nmodes ),  info, ifail( nmodes), isq (48), nsq, sym_sgq(48)
  ! 
  logical :: eqvect_strict
  ! True if we are using time reversal
  !
  real(kind=DP) ::  arg, rwork( 7*nmodes ), aq(3), saq(3), raq(3)
  complex(kind=DP) :: cfac, cwork( 2*nmodes )
  real(kind=DP) :: rtau (3, 48, nat)
  !  the relative position of the rotated atom to the original one
  real(kind=DP) :: xq(3)
  !  the rotated q vector
  !
  integer :: iq_irr, nqc_irr, nq, iq_first, imq, current_iq, mq
  !  the index of the irreducible q point 
  !  the number of irreducible qpoints 
  !  the number of q points in the star of q
  !  the index of the first qpoint to be read in the uniform q-grid
  !  flag which tells whether we have to consider the -q vectors
  !  the q index in the dynq matrix
  integer :: ism1, imode, jmode
  ! 
  real(kind=DP) :: sxq(3,48), scart(3,3), massfac, w1( nmodes ),wtmp( nmodes )
  !  the q vectors in the star
  complex(kind=DP) :: gamma(nmodes, nmodes), dynq_tmp(nmodes, nmodes), &
                      cz1(nmodes, nmodes), dynp( nmodes*(nmodes+1)/2 )
  !  the Gamma matrix for the symmetry operation on the dyn mat
  !
  ! output
  !
  ! dynq (nmode,nmodes,nqc) (in elph2.mod)
  !
  ! local 
  !
  character(len=3)    :: atm, filelab
  character(len=80)   :: line
  character (len=256) :: tempfile
  logical             :: found, lrigid, lraman, nog
  integer             :: ntyp_, nat_, ibrav_, ityp_, ios, iq, jq,  &
                         nt, na, nb, naa, nbb, nu, mu, i, j, ipol,jpol
  real(kind=DP)       :: eps, celldm_ (6), amass_, tau_ (3), q(3,nq1*nq2*nq3), &
                         dynr (2, 3, nat, 3, nat), sumz
  COMPLEX(KIND=DP),allocatable   :: dyn(:,:,:,:) ! 3,3,nat,nat
  integer, parameter    :: ntypx = 10
  real(DP), allocatable :: m_loc(:,:), dchi_dtau(:,:,:,:)
  integer               :: nqs, axis
  real(DP)              :: qout(3), amass2(ntypx)
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
