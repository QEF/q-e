MODULE restart

  IMPLICIT NONE
  SAVE

CONTAINS

!-----------------------------------------------------------------------
      subroutine writefile_new                                         &
     &     ( ndw,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,  &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm, &
     &       fion, tps)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      use ions_base, only: nsp, na
      use cell_base, only: boxdimensions, s_to_r, cell_init
      use elct, only: n, nx, nspin, nel
      use gvecw, only: ngw, ngwt
!      use reciprocal_vectors, only: gstart
      use parameters, only: natx
      use grid_dimensions, ONLY: nr1, nr2, nr3
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s
      use gvecp, only: &
        ng => ngm, &
        ngl => ngml, &
        ng_g => ngmt
      use gvec, ONLY: mill_g, mill_l, bi1, bi2, bi3, ig_l2g
      use io_base, only: write_restart_header, write_restart_ions, &
          write_restart_cell, write_restart_electrons, &
          write_restart_gvec, write_restart_gkvec, write_restart_charge, &
          write_restart_wfc, write_restart_symmetry, &
          write_restart_xdim, write_restart_pseudo, write_restart_tetra
      use mp, only: mp_sum
      use mp_global
      use io_global
      USE atom, ONLY: r, rab
      use control_flags, only: twfcollect
      USE parser, ONLY: int_to_char
      use input_parameters, only: scradir
!
      implicit none
      integer :: ndw, nfi
      real(kind=8) :: h(3,3), hold(3,3)
      complex(kind=8) :: c0(:,:), cm(:,:)
      real(kind=8) :: tausm(:,:),taus(:,:), fion(:,:)
      real(kind=8) :: vels(:,:), velsm(:,:)
      real(kind=8) :: acc(:),lambda(:,:), lambdam(:,:)
      real(kind=8) :: xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8), INTENT(in) :: ecut, ecutw, delt
      real(kind=8), INTENT(in) :: pmass(:)
      real(kind=8), INTENT(in) :: celldm(:)
      real(kind=8), INTENT(in) :: tps
      integer, INTENT(in) :: ibrav
      integer :: nk = 1
      integer :: ngwkg(1), nbnd, nelt, nelu, neld, ntyp, nb_g
      integer :: nat = 0
      integer :: nacx = 10
      real(kind=8) :: ecutwfc, ecutrho

      REAL(kind=8), ALLOCATABLE :: stau0(:,:), staum(:,:), svel0(:,:), svelm(:,:), tautmp(:,:)
      REAL(kind=8), ALLOCATABLE :: fiontmp(:,:)
      type (boxdimensions) :: box
      real(kind=8) :: ht(3,3), htvel(3,3), ht0(3,3), htm(3,3), htm2(3,3)
      real(kind=8) :: xdum
      real(kind=8) :: hdum(3,3)
      real(kind=8) :: cdmi(3)
      real(kind=8) :: mass(nsp)
      real(kind=8), allocatable :: occ(:), occm(:), lamtmp(:,:), lamtmpm(:,:), eigtmp(:)
      LOGICAL :: tocc, tlam, trho, tv, tw0, twm
      integer, allocatable :: mill(:,:), igk(:)
      real(kind=8) :: xk(3), wk
      complex(kind=8), allocatable :: rhog(:), vg(:)

      LOGICAL :: lstres, lforce
      character(len=80) :: title, crystal, tmp_dir
!
      integer i, ia, is, j, ispin, ik
      integer :: strlen
      character(len=80) :: filename

      INTEGER ::  k1, k2, k3, nk1, nk2, nk3
      REAL(kind=8) :: dgauss
      INTEGER :: ngauss
      INTEGER :: ntetra
      INTEGER :: natomwfc
      LOGICAL :: doublegrid, tupf
      REAL(kind=8) :: gcutm, gcuts, dual
      INTEGER :: modenum
      REAL(kind=8) :: alat
      REAL(kind=8) :: ef, rnel, wfc_scal_cp90
      character(len=4) :: atom_label(nsp)

      LOGICAL :: tscal
      LOGICAL :: teig
      LOGICAL :: tmill
      LOGICAL :: lgauss
      LOGICAL :: ltetra
      LOGICAL :: lgamma
      LOGICAL :: lda_plus_u
      LOGICAL :: noncolin, lspinorb
      INTEGER, ALLOCATABLE :: ityp(:)
      INTEGER :: isk
      REAL(kind=8) :: zmesh_, xmin_, dx_
      REAL(kind=8) :: ainv(3,3), deth
      INTEGER :: iswitch = 0
      LOGICAL :: tfixed_occ_, tefield_, dipfield_
      INTEGER :: edir_
      REAL(kind=8) :: emaxpos_, eopreg_, eamp_

!
! Do not write restart file if the unit number 
! is negative, this is used mainly for benchmarks
! and tests
!

      if ( ndw < 1 ) then
        return
      end if

!
! Only the first node writes
!
      if (ionode) then
       !  open (unit=ndw,status='unknown',form='unformatted')
         filename = 'fort.'//int_to_char( ndw )
         strlen  = index(scradir,' ') - 1 
         if( strlen >= 1 ) then
           filename = scradir(1:strlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1 
         OPEN(unit=ndw, file=filename(1:strlen), form='unformatted', status='unknown')
         REWIND ndw
      end if

      ht = TRANSPOSE(h)
      call cell_init(box,ht)

!       ==--------------------------------------------------------------==
!       ==  WRITE HEADER INFORMATION                                    ==
!       ==--------------------------------------------------------------==

      ngwkg(1) = ngwt
      nbnd = n
      IF( nspin > 1 ) THEN
        nelu = nel(1)
        neld = nel(2)
      ELSE
        nelu = 0
        neld = 0
      END IF
      nelt = nel(1)+nel(2)
      rnel = REAL(nelt)
      ntyp = nsp
      nat = SUM( na(1:ntyp) )
      ecutwfc = ecutw
      ecutrho = ecut
      title   = ' '
      crystal = ' '
      tmp_dir = ' '
      kunit = 1
      lgauss = .FALSE.
      ltetra = .FALSE.
      ntetra = 0
      tupf   = .TRUE.
      lgamma = .TRUE.
      noncolin = .FALSE.
      lspinorb = .FALSE.
      lda_plus_u = .FALSE.
      tfixed_occ_ = .FALSE.
      tefield_ = .FALSE.
      dipfield_ = .FALSE.
      edir_ = 0
      emaxpos_ = 0.0d0
      eopreg_ = 0.0d0
      eamp_ = 0.0d0
      CALL write_restart_header(ndw, nfi, iswitch, tps, nr1, nr2, nr3, &
        nr1s, nr2s, nr3s, ng_g, nk, ngwkg, nspin, nbnd, rnel, nelu, &
        neld, nat, ntyp, na, acc, nacx, ecutwfc, ecutrho, alat, ekincm, &
        kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
        natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, &
        title, crystal, tmp_dir, tupf, lgamma, noncolin, lspinorb, &
        lda_plus_u, &
        tfixed_occ_, tefield_, dipfield_, edir_, emaxpos_, eopreg_, eamp_, twfcollect )

!     ==--------------------------------------------------------------==
!     ==  MAX DIMENSIONS                                              ==
!     ==--------------------------------------------------------------==

      CALL write_restart_xdim( ndw )

!     ==--------------------------------------------------------------==
!     ==  CELL & METRIC                                               ==
!     ==--------------------------------------------------------------==

      hdum = 0.0d0
      ht = TRANSPOSE(h) 
      htm = TRANSPOSE(hold) 
      htvel = TRANSPOSE(velh) 
      htm2 = 0.0d0
      CALL write_restart_cell( ndw, ibrav, celldm, ht, htm, &
          htm2, htvel, vnhh, xnhh0, xnhhm, hdum)

!     ==--------------------------------------------------------------==
!     ==  IONS                                                        ==
!     ==--------------------------------------------------------------==

      ALLOCATE( stau0(3, nat) )
      ALLOCATE( staum(3, nat) )
      ALLOCATE( svel0(3, nat) )
      ALLOCATE( svelm(3, nat) )
      ALLOCATE( tautmp(3, nat) )
      ALLOCATE( fiontmp(3, nat) )
      ALLOCATE( ityp(nat) )

      ia = 0
      DO i = 1, nsp
        DO j = 1, na(i)
          ia = ia + 1
          stau0(:,ia) = taus(:,ia)
          staum(:,ia) = tausm(:,ia)
          svel0(:,ia) = vels(:,ia)
          svelm(:,ia) = velsm(:,ia)
          CALL s_to_r( taus(:,ia), tautmp(:,ia), box )
          fiontmp(:,ia) = fion(:,ia)
          ityp(ia) = i
        END DO
        mass(i) = pmass(i)
      END DO

      xdum = 0.0d0
      cdmi = 0.0d0
      tscal = .TRUE.
      CALL write_restart_ions(ndw, atom_label, tscal, stau0, svel0, &
        staum, svelm, tautmp, fiontmp, cdmi, nat, ntyp, ityp, na, mass,     &
        vnhp, xnhp0, xnhpm, xdum)

      DEALLOCATE( stau0, staum, svel0, svelm, tautmp, fiontmp, ityp )

!     ==--------------------------------------------------------------==
!     ==  SYMMETRIES                                                  ==
!     ==--------------------------------------------------------------==

      CALL write_restart_symmetry( ndw )

!     ==--------------------------------------------------------------==
!     ==  PSEUDOPOTENTIALS                                            ==
!     ==--------------------------------------------------------------==

      DO i = 1, nsp
!       CALL write_restart_pseudo( ndw, &
!         zmesh_, xmin_, dx_, r(:,i), rab(:,i), vnl(:,:,i), chi(:,:,i), oc(:,i), &
!         rho_at(:,i), rho_atc(:,i), mesh(i), msh(i), nchi(i), lchi(:,i), &
!         numeric(i), cc(:,i), alpc(:,i), zp(i), aps(:,:,i), alps(:,:,i), &
!         zv(i), nlc(i), nnl(i), lmax(i), lloc(i), dion(:,:,i), &
!         betar(:,:,i), qqq(:,:,i), qfunc(:,:,:,i), qfcoef(:,:,:,:,i), &
!         rinner(:,i), nh(i), nbeta(i), kkbeta(i), nqf(i), nqlc(i), ifqopt(i), &
!         lll(:,i), iver(:,i), tvanp(i), okvan, newpseudo(i), iexch, icorr, &
!         igcx, igcc, lsda, a_nlcc(i), b_nlcc(i), alpha_nlcc(i), nlcc(i), psd(i) )
       CALL write_restart_pseudo( ndw )
      END DO


!     ==--------------------------------------------------------------==
!     ==  OCCUPATION NUMBER                                           ==
!     ==--------------------------------------------------------------==

      ALLOCATE( occ(nbnd), eigtmp(nbnd) )

      occ = 0.0d0
      ik = 1
      ispin = 1
      rnel = REAL(nelt)
      tocc = .FALSE.
      tlam = .TRUE.
      teig = .FALSE.
      CALL write_restart_electrons( ndw, occ, occ, tocc, lambda, lambdam, &
        nx, tlam, nbnd, ispin, nspin, ik, nk, rnel, nelu, neld, vnhe, xnhe0, xnhem, xdum, &
        ef, teig, eigtmp, eigtmp)

      DEALLOCATE( occ, eigtmp )


!     ==--------------------------------------------------------------==
!     ==  G-Vectors                                                   ==
!     ==--------------------------------------------------------------==

      ALLOCATE( mill(3,ng_g) )
      mill = 0
      DO i = 1, ng
        mill(:,ig_l2g(i)) = mill_l(:,i)
      END DO
      CALL mp_sum( mill )
      tmill = .TRUE.
      CALL write_restart_gvec( ndw, ng_g, bi1, bi2, bi3, &
        bi1, bi2, bi3, tmill, mill )
      DEALLOCATE( mill )


!     ==--------------------------------------------------------------==
!     ==  (G+k)-Vectors                                               ==
!     ==--------------------------------------------------------------==

      DO i = 1, nk
        xk(1) = 0.0d0
        xk(2) = 0.0d0
        xk(3) = 0.0d0
        wk = 1.0d0
        isk = 1
        CALL write_restart_gkvec(ndw, i, nk, ngwkg(i), xk, wk, isk)
      END DO

!     ==--------------------------------------------------------------==
!     ==  Tetrahedra                                                  ==
!     ==--------------------------------------------------------------==

      CALL write_restart_tetra(ndw) 

!     ==--------------------------------------------------------------==
!     ==  CHARGE DENSITY AND POTENTIALS                               ==
!     ==--------------------------------------------------------------==

      trho = .FALSE.
      tv   = .FALSE.
      DO j = 1, nspin
        ALLOCATE( rhog(ng), vg(ng) )
!        CALL fft_initialize
!        CALL pfwfft( rhog(:,i), rho(i)%r(:,:,:) )
!        CALL pfwfft( vg(:,i)  , vpot(:,:,:,i) )
        CALL write_restart_charge(ndw, rhog, trho, vg, tv, ng_g, &
          j, nspin, ig_l2g, ng )
        DEALLOCATE( rhog, vg )
      END DO


!     ==--------------------------------------------------------------==
!     ==  WAVEFUNCTIONS                                               ==
!     ==--------------------------------------------------------------==

      tw0 = .TRUE.
      twm = .TRUE.
        
      call invmat (3, h, ainv, deth)
      wfc_scal_cp90 = 1.0d0 / SQRT(ABS(deth))
      DO j = 1, nspin
        DO i = 1, nk
          nb_g = nx
          CALL write_restart_wfc(ndw, i, nk, kunit, j, nspin, &
            wfc_scal_cp90, c0, tw0, cm, twm, ngwt, nb_g, ig_l2g, ngw )
        END DO
      END DO

      if(ionode) then
         close (unit=ndw)
      end if


      return
      end subroutine

!-----------------------------------------------------------------------
      subroutine readfile_new                                         &
     &     ( flag, ndr,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm, &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,&
     &       fion, tps)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!


      use elct, only: n, nx, nspin, nel
      use gvecw, only: ngw, ngwt
!      use reciprocal_vectors, only: gstart
      use ions_base, only: nsp, na
      use parameters, only: natx, nsx
      use grid_dimensions, ONLY: nr1, nr2, nr3
      use gvecp, only: &
        ng => ngm, &
        ngl => ngml, &
        ng_g => ngmt
      use gvec, ONLY: mill_g, mill_l, bi1, bi2, bi3, ig_l2g
      use io_base, only: read_restart_header, read_restart_ions, &
          read_restart_cell, read_restart_electrons, &
          read_restart_gvec, read_restart_gkvec, read_restart_charge, &
          read_restart_wfc, read_restart_xdim, read_restart_tetra, &
          read_restart_symmetry, read_restart_pseudo
      use mp, only: mp_sum
      use mp_global
      use io_global
      use cell_base, only: boxdimensions, s_to_r, cell_init, r_to_s
      use control_flags, only: twfcollect
      USE parser, ONLY: int_to_char
      use input_parameters, only: scradir

!
      implicit none
      integer :: ndr, nfi, flag
      real(kind=8) :: h(3,3), hold(3,3)
      complex(kind=8) :: c0(:,:), cm(:,:)
      real(kind=8) :: tausm(:,:),taus(:,:), fion(:,:)
      real(kind=8) :: vels(:,:), velsm(:,:)
      real(kind=8) :: acc(:),lambda(:,:), lambdam(:,:)
      real(kind=8) :: xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8), INTENT(in) :: ecut, ecutw, delt
      real(kind=8), INTENT(in) :: pmass(:)
      real(kind=8), INTENT(in) :: celldm(6)
      integer, INTENT(in) :: ibrav
      real(kind=8) :: tps

      integer :: nx_, nbnd_, nelt, nb_g
      integer :: nacx = 10
      real(kind=8) :: ecutwfc, ecutrho

      REAL(kind=8), ALLOCATABLE :: stau0(:,:), staum(:,:), svel0(:,:), svelm(:,:), tautmp(:,:)
      REAL(kind=8), ALLOCATABLE :: fiontmp(:,:)
      type (boxdimensions) :: box, boxm
      real(kind=8) :: ht(3,3), htvel(3,3), ht0(3,3), htm(3,3), htm2(3,3)
      real(kind=8) :: xdum
      real(kind=8) :: hdum(3,3)
      real(kind=8) :: cdmi(3)
      real(kind=8) :: mass(nsp)
      real(kind=8), allocatable :: occ(:), eigtmp(:)
      LOGICAL :: tocc, tlam, trho, tv, tw0, twm
      integer, allocatable :: mill(:,:), igk(:)
      real(kind=8) :: xk(3), wk
      complex(kind=8), allocatable :: rhog(:), vg(:)

      INTEGER ::  k1_, k2_, k3_, nk1_, nk2_, nk3_
      REAL(kind=8) :: dgauss_
      INTEGER :: ngauss_
      LOGICAL :: lgauss_
      INTEGER :: ntetra_
      LOGICAL :: ltetra_
      INTEGER :: natomwfc_
      LOGICAL :: doublegrid_
      REAL(kind=8) :: gcutm_, gcuts_, dual_
      INTEGER :: modenum_
      REAL(kind=8) :: ef, rnel
      LOGICAL :: lstres_, lforce_
      character(len=80) :: title_, crystal_, tmp_dir_
      character(len=4) :: atom_label(nsp)
      real(kind=8) :: bi1_(3), bi2_(3), bi3_(3)
!
      integer :: i, ia, is, j
      integer :: nfi_, ik_, nk, nk_, ispin_, nspin_, isk_
      real(kind=8) :: acc_(10), celldm_(6)
      real(kind=8) :: vnhp_, xnhp0_, xnhpm_
      integer :: strlen, ibrav_, kunit_
      character(len=80) :: filename
      LOGICAL :: tscal
      LOGICAL :: tmill, tigl, lgamma_, lda_plus_u_
      LOGICAL :: noncolin_, lspinorb_
      LOGICAL :: teig, tupf_
      INTEGER, ALLOCATABLE :: ityp(:)
      REAL(kind=8) :: wfc_scal, wfc_scal_cp90

      LOGICAL :: tfixed_occ_, tefield_, dipfield_
      INTEGER :: edir_
      REAL(kind=8) :: emaxpos_, eopreg_, eamp_
      INTEGER :: nr1_, nr2_, nr3_, nr1s_, nr2s_, nr3s_, ng_g_, ngwkg_(1)
      INTEGER :: ngwt_
      REAL(kind=8) :: rnel_
      INTEGER :: nelu_, neld_
      INTEGER :: nat_, ntyp_, na_( nsx )
      INTEGER :: nacx_
      REAL(kind=8) :: ecutwfc_, ecutrho_
      REAL(kind=8) :: alat_, ekincm_ 
      LOGICAL :: twfcollect_

      INTEGER :: iswitch_ 


!
! Only the first node read
!
      if (ionode) then
         ! open (unit=ndr, status='old', form='unformatted')
         filename = 'fort.'//int_to_char( ndr ) 
         strlen  = index(scradir,' ') - 1 
         if( strlen >= 1 ) then
           filename = scradir(1:strlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1 
         OPEN(unit=ndr, file=filename(1:strlen), form='unformatted', status='old')
         REWIND (ndr)
         WRITE( stdout,10)
 10      FORMAT(/,3X,'READING FROM RESTART FILE ...')
      end if

      if ( flag == -1 ) then
         WRITE( stdout,'((a,i3,a))') ' ### reading from file ',ndr,' only h  ##'
      else if ( flag == 0 ) then
         WRITE( stdout,'((a,i3,a))') ' ## reading from file ',ndr,' only c0  ##'
      else
         WRITE( stdout,'((a,i3))') ' ## reading from file ',ndr
      end if

!     ==--------------------------------------------------------------==
!     ==  READ HEADER INFORMATION                                     ==
!     ==--------------------------------------------------------------==

      CALL read_restart_header( ndr, nfi_, iswitch_, tps, &
        nr1_, nr2_, nr3_, nr1s_, nr2s_, nr3s_, ng_g_, nk_, ngwkg_, nspin_, nbnd_, &
        rnel_, nelu_, neld_, nat_, ntyp_, na_, acc_, nacx_, ecutwfc_, ecutrho_, &
        alat_, ekincm_, kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_, dgauss_, ngauss_, &
        lgauss_, ntetra_, ltetra_, natomwfc_, gcutm_, gcuts_, dual_, doublegrid_, &
        modenum_, lstres_, lforce_, title_, crystal_, tmp_dir_, tupf_, &
        lgamma_, noncolin_, lspinorb_, lda_plus_u_, &
        tfixed_occ_, tefield_, dipfield_, edir_, emaxpos_, eopreg_, eamp_, twfcollect_ )

      if( .not. lgamma_ ) &
        call errore(' readfile ',' restart contains a system not at gamma ',1)

      if( nbnd_ /= n ) &
        call errore(' readfile ',' inconsistent number of bands in restart ',1)
      if( nk_ /= 1 ) &
        call errore(' readfile ',' inconsistent number of kpoints in restart ',1)

      ekincm = ekincm_
      if ( flag > -1 ) then
        nfi = nfi_
        acc = acc_
      end if
      nk = nk_

!     ==--------------------------------------------------------------==
!     ==  MAX DIMENSIONS                                              ==
!     ==--------------------------------------------------------------==

      CALL read_restart_xdim( ndr )
      
!     ==--------------------------------------------------------------==
!     ==  CELL & METRIC                                               ==
!     ==--------------------------------------------------------------==

      hdum = 0.0d0
      htm2 = 0.0d0
        
      CALL read_restart_cell( ndr, ibrav_, celldm_, ht, htm, &
        htm2, htvel, vnhh, xnhh0, xnhhm, hdum)

      h    = TRANSPOSE( ht    ) 
      hold = TRANSPOSE( htm   ) 
      velh = TRANSPOSE( htvel ) 

      CALL cell_init(box, ht)
      CALL cell_init(boxm, htm)

!     ==--------------------------------------------------------------==
!     ==  IONS                                                        ==
!     ==--------------------------------------------------------------==

        ALLOCATE( stau0(3, nat_) )
        ALLOCATE( staum(3, nat_) )
        ALLOCATE( svel0(3, nat_) )
        ALLOCATE( svelm(3, nat_) )
        ALLOCATE( tautmp(3, nat_) )
        ALLOCATE( fiontmp(3, nat_) )
        ALLOCATE( ityp(nat_) )

        xdum = 0.0d0
        cdmi = 0.0d0
        CALL read_restart_ions(ndr, atom_label, tscal, stau0, svel0, &
          staum, svelm, tautmp, fiontmp, cdmi, nat_, ntyp_, ityp, na_, mass, vnhp_,   &
          xnhp0_, xnhpm_, xdum)

        IF( flag > 0 ) THEN
          ia = 0
          DO i = 1, nsp
            DO j = 1, na(i)
              ia = ia + 1
              IF( tscal ) THEN
                taus(:,ia) = stau0(:,ia)
                tausm(:,ia) = staum(:,ia)
                vels(:,ia) = svel0(:,ia)
                velsm(:,ia) = svelm(:,ia)
              ELSE
                CALL r_to_s( stau0(:,ia), taus(:,ia), box )
                CALL r_to_s( staum(:,ia), tausm(:,ia), boxm )
                CALL r_to_s( svel0(:,ia), vels(:,ia), box )
                CALL r_to_s( svelm(:,ia), velsm(:,ia), boxm )
              END IF
            END DO
          END DO
          xnhp0 = xnhp0_
          xnhpm = xnhpm_
          vnhp  = vnhp_
        END IF

        DEALLOCATE( stau0, staum, svel0, svelm, tautmp, ityp, fiontmp )


!       ==--------------------------------------------------------------==
!       ==  SYMMETRIES                                                  ==
!       ==--------------------------------------------------------------==

        CALL read_restart_symmetry( ndr )

!       ==--------------------------------------------------------------==
!       ==  PSEUDOPOTENTIALS                                            ==
!       ==--------------------------------------------------------------==

        DO i = 1, nsp
          CALL read_restart_pseudo( ndr )
        END DO

!       ==--------------------------------------------------------------==
!       ==  OCCUPATION NUMBER                                           ==
!       ==--------------------------------------------------------------==

        ALLOCATE( occ(n), eigtmp(n) )

        occ = 0.0d0
        IF( flag > 0 ) THEN
          tocc = .FALSE.
          tlam = .TRUE.
        ELSE
          tocc = .FALSE.
          tlam = .FALSE.
        END IF
        teig = .FALSE.
        CALL read_restart_electrons( ndr, occ, occ, tocc, lambda, &
          lambdam, nx_, tlam, nbnd_, ispin_, nspin_, ik_, nk_, rnel_, nelu_, neld_,    &
          vnhe, xnhe0, xnhem, xdum, ef, teig, eigtmp, eigtmp)

        DEALLOCATE( occ, eigtmp )

!       ==--------------------------------------------------------------==
!       ==  G-Vectors                                                   ==
!       ==--------------------------------------------------------------==

        ALLOCATE( mill(3,ng_g) )
        mill = 0
        DO i = 1, ng
          mill(:,ig_l2g(i)) = mill_l(:,i)
        END DO
        CALL mp_sum( mill )
        tmill = .FALSE.
        CALL read_restart_gvec( ndr, ng_g_, bi1_, bi2_, bi3_,  &
          bi1_, bi2_, bi3_, tmill, mill )
        DEALLOCATE( mill )

!       ==--------------------------------------------------------------==
!       ==  (G+k)-Vectors                                               ==
!       ==--------------------------------------------------------------==

        DO i = 1, nk
          xk(1) = 0.0d0
          xk(2) = 0.0d0
          xk(3) = 0.0d0
          wk = 1.0d0
          CALL read_restart_gkvec(ndr, ik_, nk_, ngwkg_(i), &
            xk, wk, isk_)
        END DO

!       ==--------------------------------------------------------------==
!       ==  Tetrahedra                                                  ==
!       ==--------------------------------------------------------------==

        CALL read_restart_tetra( ndr )

!       ==--------------------------------------------------------------==
!       ==  CHARGE DENSITY AND POTENTIALS                               ==
!       ==--------------------------------------------------------------==

        trho = .FALSE.
        tv   = .FALSE.
        DO j = 1, nspin
          ALLOCATE( rhog(ng), vg(ng) )
          CALL read_restart_charge(ndr, rhog, trho, vg, tv, ng_g_, &
            ispin_, nspin_, ig_l2g, ng )
!          CALL fft_initialize
!          CALL pinvfft( vpot(:,:,:,i), rhog(:,i) )
!          CALL pinvfft( rho(i)%r(:,:,:), vg(:,i) )
          DEALLOCATE( rhog, vg )
        END DO

!       ==--------------------------------------------------------------==
!       ==  WAVEFUNCTIONS                                               ==
!       ==--------------------------------------------------------------==

        DO j = 1, nspin
          tigl = .FALSE.
          IF( flag == -1 ) THEN
            tw0 = .FALSE.
            twm = .FALSE.
          ELSE IF( flag == 0 ) THEN
            tw0 = .TRUE.
            twm = .FALSE.
          ELSE
            tw0 = .TRUE.
            twm = .TRUE.
          END IF
          DO i = 1, nk
            CALL read_restart_wfc(ndr, ik_, nk_, kunit_, ispin_, nspin_, &
              wfc_scal, c0, tw0, cm, twm, ngwt_, nbnd_, ig_l2g, ngw )
          END DO
        END DO


      if(ionode) then
         close (unit=ndr)
      end if

      return
      end subroutine



!=----------------------------------------------------------------------------=!

        LOGICAL FUNCTION check_restartfile( scradir, ndr )

          USE io_global, ONLY: ionode, ionode_id
          USE mp, ONLY: mp_bcast
          USE parser, ONLY: int_to_char

          IMPLICIT NONE

          INTEGER, INTENT(IN) :: ndr
          CHARACTER(LEN=*) :: scradir
          CHARACTER(LEN=256) :: filename
          LOGICAL :: lval
          INTEGER :: strlen

          lval = .FALSE.
          IF ( ionode ) THEN
            filename = 'fort.' // int_to_char( ndr )
            strlen  = INDEX( scradir, ' ' ) - 1
            filename = scradir( 1 : strlen ) // '/' // filename
            INQUIRE( FILE = TRIM( filename ), EXIST = lval )
          END IF
          CALL mp_bcast( lval, ionode_id )
          check_restartfile = lval
          RETURN
        END FUNCTION


END MODULE restart
