!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
  MODULE restart_file
!=----------------------------------------------------------------------------=!

  USE kinds, ONLY: dbl

  IMPLICIT NONE

  PRIVATE

  SAVE

  REAL(dbl) :: cclock
  EXTERNAL  :: cclock

  PUBLIC :: writefile, readfile, check_restartfile

  INTERFACE readfile
    MODULE PROCEDURE readfile_cp, readfile_fpmd
  END INTERFACE

  INTERFACE writefile
    MODULE PROCEDURE writefile_cp, writefile_fpmd
  END INTERFACE
  
!=----------------------------------------------------------------------------=!
     CONTAINS
!=----------------------------------------------------------------------------=!

#ifdef __OLD_RESTART

!-----------------------------------------------------------------------
      subroutine writefile_cp_new                                         &
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
      character(len=80) :: title, crystal
      character(len=256) :: tmp_dir !
      integer i, ia, is, j, ispin, ik
      integer :: strlen
      character(len=256) :: filename

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
      LOGICAL :: tfixed_occ_, tefield_, dipfield_
      INTEGER :: edir_
      REAL(kind=8) :: emaxpos_, eopreg_, eamp_

!
! Do not write restart file if the unit number 
! is negative, this is used mainly for benchmarks
! and tests
!

        IF( ndw < 1 ) RETURN
        !
        !   this is used for benchmarking and debug
        !   if ndw < 1 Do not save wave functions and other system
        !   properties on the writefile subroutine

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
      CALL write_restart_header(ndw, nfi, tps, nr1, nr2, nr3, &
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
      subroutine readfile_cp_new                                         &
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
      character(len=80) :: title_, crystal_
      character(len=256) :: tmp_dir_
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

      CALL read_restart_header( ndr, nfi_, tps, &
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


#endif

!=----------------------------------------------------------------------------=!

!-----------------------------------------------------------------------
      subroutine writefile_cp                                         &
     &     ( ndw,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,  &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm, &
     &       fion, tps, mat_z, occ_f )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      use ions_base, only: nsp, na
      use cell_base, only: s_to_r
      use cp_restart, only: cp_writefile, cp_write_wfc
      USE electrons_base, ONLY: nspin, nbnd, iupdwn, nupdwn
!
      implicit none
      integer, INTENT(IN) :: ndw, nfi
      real(kind=8), INTENT(IN) :: h(3,3), hold(3,3)
      complex(kind=8), INTENT(IN) :: c0(:,:), cm(:,:)
      real(kind=8), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
      real(kind=8), INTENT(IN) :: vels(:,:), velsm(:,:)
      real(kind=8), INTENT(IN) :: acc(:), lambda(:,:), lambdam(:,:)
      real(kind=8), INTENT(IN) :: xnhe0, xnhem, vnhe, xnhp0, xnhpm, vnhp, ekincm
      real(kind=8), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8), INTENT(in) :: ecut, ecutw, delt
      real(kind=8), INTENT(in) :: pmass(:)
      real(kind=8), INTENT(in) :: celldm(:)
      real(kind=8), INTENT(in) :: tps
      integer, INTENT(in) :: ibrav
      real(kind=8), INTENT(in) :: mat_z(:,:,:), occ_f(:)

      real(kind=8) :: ht(3,3), htm(3,3), htvel(3,3), htm2_(3,3) , xnhhp_(3,3)
      integer :: nk = 1, ispin, i
      real(kind=8) :: xk(3,1) = 0.0d0, wk(1) = 1.0d0
      real(kind=8) :: cdmi_ (3) = 0.0d0
      real(kind=8), ALLOCATABLE :: taui_ (:,:) 
      real(kind=8) :: xnhpp_ , xnhpm2_
      real(kind=8), ALLOCATABLE :: occ_ ( :, :, : )
      real(kind=8) :: xnhep_ , xnhem2_
      real(kind=8) :: htm1(3,3), b1(3) , b2(3), b3(3), omega

!
! Do not write restart file if the unit number 
! is negative, this is used mainly for benchmarks
! and tests
!

      if ( ndw < 1 ) then
        return
      end if

      ht     = TRANSPOSE( h ) 
      htm    = TRANSPOSE( hold ) 
      htvel  = TRANSPOSE( velh ) 
      htm2_  = htm
      xnhhp_ = 0.0d0

      ALLOCATE( taui_ ( 3, SIZE( taus, 2 ) ) )
      CALL s_to_r( taus, taui_ , na, nsp, h )

      cdmi_ = 0.0d0
      xnhpp_  = xnhp0
      xnhpm2_ = xnhpm

      ALLOCATE( occ_ ( nbnd, 1, nspin ) )
      occ_ = 0.0d0
      do ispin = 1, nspin
        do i = iupdwn ( ispin ), iupdwn ( ispin ) - 1 + nupdwn ( ispin )
          occ_ ( i - iupdwn ( ispin ) + 1, 1, ispin ) = occ_f( i ) 
        end do
      end do


      CALL invmat( 3, ht, htm1, omega )
      b1 = htm1(:,1)
      b2 = htm1(:,2)
      b3 = htm1(:,3)

      xnhep_ = 0.0d0
      xnhem2_ = 0.0d0

      CALL cp_writefile( ndw, ' ', .TRUE., nfi, tps, acc, nk, xk, wk, &
        ht, htm, htm2_ , htvel, xnhh0, xnhhm, xnhhp_ , vnhh, taui_ , cdmi_ , taus, &
        vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, xnhpp_ , xnhpm2_ , occ_ , &
        occ_ , lambda, lambdam, b1, b2, b3, &
        xnhe0, xnhem, xnhep_ , xnhem2_ , vnhe, ekincm, mat_z  )

      DEALLOCATE( taui_ )
      DEALLOCATE( occ_ )

      IF( nspin /= 1 ) CALL errore( ' writefile ', ' nspin > 1 ', 1 )

      DO ispin = 1, nspin
        CALL cp_write_wfc( ndw, ' ', 1, 1, ispin, nspin, c0, '0' )
        CALL cp_write_wfc( ndw, ' ', 1, 1, ispin, nspin, cm, 'm' )
      END DO

      return
      end subroutine

!-----------------------------------------------------------------------
      subroutine readfile_cp                                        &
     &     ( flag, ndr,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm, &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,&
     &       fion, tps, mat_z, occ_f )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!

      use electrons_base, only: nbnd, nspin, nel, nupdwn, iupdwn
      use gvecw, only: ngw, ngwt
      use ions_base, only: nsp, na
      use cp_restart, only: cp_readfile, cp_read_cell, cp_read_wfc
      use ensemble_dft, only: tens
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
      real(kind=8), INTENT(OUT) :: tps
      real(kind=8), INTENT(INOUT) :: mat_z(:,:,:), occ_f(:)
      !
      real(kind=8) :: ht(3,3), htm(3,3), htvel(3,3), htm2_(3,3) , xnhhp_(3,3)
      integer :: nk = 1, ispin, i
      real(kind=8) :: xk(3,1) = 0.0d0, wk(1) = 1.0d0
      real(kind=8) :: cdmi_ (3) = 0.0d0
      real(kind=8), ALLOCATABLE :: taui_ (:,:)
      real(kind=8) :: xnhpp_ , xnhpm2_
      real(kind=8), ALLOCATABLE :: occ_ ( :, :, : )
      real(kind=8) :: xnhep_ , xnhem2_
      real(kind=8) :: htm1(3,3), b1(3) , b2(3), b3(3), omega


      IF( nspin /= 1 ) CALL errore( ' writefile ', ' nspin > 1 ', 1 )

      IF( flag == -1 ) THEN
        CALL cp_read_cell( ndr, ' ', .TRUE., ht, htm, htvel, xnhh0, xnhhm, xnhhp_ , vnhh )
        h     = TRANSPOSE( ht )
        hold  = TRANSPOSE( htm )
        velh  = TRANSPOSE( htvel )
        RETURN
      ELSE IF ( flag == 0 ) THEN
        DO ispin = 1, nspin
          CALL cp_read_wfc( ndr, ' ', 1, 1, ispin, nspin, cm, 'm' )
        END DO
        RETURN
      END IF

      ALLOCATE( taui_ ( 3, SIZE( taus, 2 ) ) )
      ALLOCATE( occ_ ( nbnd, 1, nspin ) )

      CALL cp_readfile( ndr, ' ', .TRUE., nfi, tps, acc, nk, xk, wk, &
        ht, htm, htm2_ , htvel, xnhh0, xnhhm, xnhhp_ , vnhh, taui_ , cdmi_ , taus, &
        vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, xnhpp_ , xnhpm2_ , occ_ , &
        occ_ , lambda, lambdam, b1, b2, b3, &
        xnhe0, xnhem, xnhep_ , xnhem2_ , vnhe, ekincm, mat_z, tens  )

      DEALLOCATE( taui_ )

      do ispin = 1, nspin
        do i = iupdwn ( ispin ), iupdwn ( ispin ) - 1 + nupdwn ( ispin )
          occ_f( i ) = occ_ ( i - iupdwn ( ispin ) + 1, 1, ispin )
        end do
      end do
      DEALLOCATE( occ_ )

      h     = TRANSPOSE( ht )
      hold  = TRANSPOSE( htm )
      velh  = TRANSPOSE( htvel )

      DO ispin = 1, nspin
        CALL cp_read_wfc( ndr, ' ', 1, 1, ispin, nspin, c0, '0' )
        CALL cp_read_wfc( ndr, ' ', 1, 1, ispin, nspin, cm, 'm' )
      END DO

      return
      end subroutine


!=----------------------------------------------------------------------------=!

#ifdef __OLD_RESTART

        SUBROUTINE writefile_fpmd_new( nfi, trutime, &
          c0, cm, cdesc, occ, atoms_0, atoms_m, acc, taui, cdmi, ibrav, celldm, &
          ht_m2, ht_m, ht_0, rho, desc, vpot, gv, kp)
                                                                        
        use io_files, only: scradir
        use electrons_module, only: nel, nelt
        use nose_ions, only: get_nose_ions 
        use nose_electrons, only: get_nose_electrons 
        use cp_types, only: recvecs
        USE cell_module, only: boxdimensions, r_to_s
        USE brillouin, only: kpoints
        use parameters, ONLY: nacx, nsx, npkx
        USE mp, ONLY: mp_max, mp_bcast, mp_sum
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp_wave, ONLY: mergewf
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndw, gamma_only
        USE atoms_type_module, ONLY: atoms_type
        USE io_base, only: write_restart_header, write_restart_ions, &
            write_restart_cell, write_restart_electrons, &
            write_restart_gvec, write_restart_gkvec, write_restart_charge, &
            write_restart_wfc, write_restart_symmetry, &
            write_restart_xdim, write_restart_pseudo, write_restart_tetra
        USE io_global, ONLY: ionode, ionode_id
        USE gvecw, ONLY: ecutwfc => ecutw
        USE gvecp, ONLY: ecutrho => ecutp
        USE fft, ONLY : pfwfft, pinvfft
        USE charge_types, ONLY: charge_descriptor
        USE control_flags, ONLY: twfcollect, force_pairing
        USE io_global, ONLY: stdout
        USE parser, ONLY: int_to_char
        USE cell_nose, ONLY: vnhh, xnhh0, xnhhm
        USE grid_dimensions, ONLY: nr1, nr2, nr3
        USE cp_restart, ONLY: cp_writefile
                                                                        
        IMPLICIT NONE 
 
        INTEGER, INTENT(IN) :: nfi, ibrav
        COMPLEX(dbl), INTENT(IN) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(dbl), INTENT(IN) :: occ(:,:,:)
        TYPE (kpoints), INTENT(IN) :: kp 
        TYPE (boxdimensions), INTENT(IN) :: ht_m2, ht_m, ht_0
        TYPE (recvecs), INTENT(IN) :: gv
        TYPE (atoms_type), INTENT(IN) :: atoms_0, atoms_m
        REAL(dbl), INTENT(IN) :: rho(:,:,:,:)
        TYPE (charge_descriptor), INTENT(IN) :: desc
        TYPE (wave_descriptor) :: cdesc
        REAL(dbl), INTENT(INOUT) :: vpot(:,:,:,:)
                                                                        
        REAL(dbl), INTENT(IN) :: taui(:,:)
        REAL(dbl), INTENT(IN) :: celldm(:)
        REAL(dbl), INTENT(IN) :: acc(:), cdmi(:) 
        REAL(dbl), INTENT(IN) :: trutime

        INTEGER  :: ngw_g, nb_g, nk, nkt, nspin, nx, ns
        INTEGER  :: nxu, nxd, ispin, ng, ng_g, ng_l, nel_, nelu_, neld_, nbnd
        INTEGER  :: ngwk_g( npkx ), nat, ntyp, na( nsx )
        INTEGER  :: i, ik, ikt, ig, ip, ios, ierr, ia, j
        REAL(dbl)  :: rnel
        REAL(dbl) :: xnosm2, xnosm, xnos0, xnosp 
        REAL(dbl) :: xenosm2, xenosm, xenos0, xenosp 
        REAL(dbl) :: xhdum(3,3)
        REAL(dbl) :: mass( nsx ), xk(3), wk
        REAL(dbl), ALLOCATABLE :: stau0(:,:), staum(:,:), svel0(:,:), svelm(:,:), tautmp(:,:)
        REAL(dbl), ALLOCATABLE :: occtmp(:), lambda(:,:)
        COMPLEX(dbl), ALLOCATABLE :: rhog(:), vg(:)
        INTEGER, ALLOCATABLE :: mill(:,:), igk(:)
        REAL(dbl)  :: s0, s1
        LOGICAL :: tocc, tlam, tw0, twm, trho, tv, teig

        INTEGER ::  kunit, k1, k2, k3, nk1, nk2, nk3
        REAL(dbl) :: dgauss
        INTEGER :: ngauss
        LOGICAL :: lgauss
        INTEGER :: ntetra
        LOGICAL :: ltetra
        INTEGER :: natomwfc
        LOGICAL :: doublegrid
        REAL(dbl) :: gcutm, gcuts, dual
        INTEGER :: modenum
        REAL(dbl) :: alat
        REAL(dbl) :: ef
        REAL(dbl) :: ekincm
        REAL(dbl) :: wfc_scal_fpmd
        LOGICAL :: tscal


        INTEGER :: strlen, ldim
        CHARACTER(LEN=256) :: filename 
        LOGICAL :: lstres = .FALSE.
        LOGICAL :: lforce = .FALSE.
        CHARACTER(LEN=80) :: title
        CHARACTER(LEN=80) :: crystal
        CHARACTER(LEN=256) :: tmp_dir
        CHARACTER(LEN=4) :: atom_label(nsx)
        INTEGER, ALLOCATABLE :: ityp(:)
        INTEGER :: isk
        LOGICAL :: tmill, tupf, lgamma, lda_plus_u
        LOGICAL :: noncolin, lspinorb
        LOGICAL :: tfixed_occ_, tefield_, dipfield_
        INTEGER :: edir_
        REAL(dbl) :: emaxpos_, eopreg_, eamp_
        INTEGER :: nspin_wfc

                                                                        
!       ==--------------------------------------------------------------==
!       OPEN THE APPROPRIATE FILE (FORT.NDW.ME) IN THE DIRECTORY CDIR     
!       ==--------------------------------------------------------------==

        s0 = cclock()

        IF( ndw < 1 ) RETURN
        !
        !   this is used for benchmarking and debug
        !   if ndw < 1 Do not save wave functions and other system
        !   properties on the writefile subroutine

        ngw_g = cdesc%ngwt
        ng_l  = gv%ng_l
        ng_g  = gv%ng_g

!       ==--------------------------------------------------------------==
!       ==  ROOT PROCESSOR WRITE TO THE SCRATCH GENERAL VARIABLES                                  
!       ==--------------------------------------------------------------==


        IF ( ionode ) THEN 
          filename = 'fort.' // int_to_char( ndw ) 
          strlen  = index(scradir,' ') - 1 
          IF ( strlen > 1 ) THEN 
            filename = scradir(1:strlen) // '/' // filename 
          END IF 
          strlen  = index(filename,' ') - 1 
          OPEN(unit=ndw, file=filename(1:strlen), form='unformatted', &
            status='unknown', iostat = ierr)
        END IF
        CALL mp_bcast( ierr, ionode_id )
        IF( ierr /= 0 ) THEN
          CALL errore(' writefile ', ' unable to open file '//filename , ierr ) 
        END IF 

!       ==--------------------------------------------------------------==
!       ==  WRITE HEADER INFORMATION                                    ==
!       ==--------------------------------------------------------------==

        nk   = cdesc%nkl
        ngwk_g = 0
        ng   = gv%ng_g
        DO i = 1, nk
          ngwk_g(i) = cdesc%ngwt
        END DO
        nspin = cdesc%nspin
        nbnd    = MAXVAL( cdesc%nbt )
        nel_    = nelt
        IF( nspin > 1 ) THEN
          nelu_   = nel(1)
          neld_   = nel(2)
          nkt     = 2 * nk
        ELSE
          nelu_   = 0
          neld_   = 0
          nkt     = nk
        END IF
        nat  = atoms_0%nat
        ntyp = atoms_0%nsp
        na = 0
        DO i = 1, ntyp
          na(i) = atoms_0%na(i)
        END DO
        alat = celldm(1)
        kunit = 1    
        k1 = 0
        k2 = 0
        k3 = 0
        nk1 = 1
        nk2 = 1
        nk3 = 1
        dgauss = 0.0d0
        ngauss = 0
        lgauss = .FALSE.
        ntetra = 0
        ltetra = .FALSE.
        natomwfc = 0
        gcutm = 0.0d0
        gcuts = 0.0d0
        dual = 0.0d0
        doublegrid = .FALSE.
        modenum = 0
        ekincm = 0.0d0
        lstres = .FALSE.
        lforce = .FALSE.
        title   = ' '
        crystal = ' '
        tmp_dir = ' '
        rnel = REAL(nelt)
        tupf = .TRUE.
        lgamma = gamma_only
        lda_plus_u = .FALSE.
        noncolin =.FALSE.
        lspinorb =.FALSE.
        tfixed_occ_ = .FALSE.
        tefield_ = .FALSE.
        dipfield_ = .FALSE.
        edir_ = 0
        emaxpos_ = 0.0d0
        eopreg_ = 0.0d0
        eamp_ = 0.0d0

        CALL write_restart_header( ndw, nfi, trutime, nr1, nr2, nr3, &
          nr1, nr2, nr3, ng_g, nkt, ngwk_g, nspin, nbnd, rnel, nelu_, &
          neld_, nat, ntyp, na, acc, nacx, ecutwfc, ecutrho, alat, ekincm, &
          kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
          natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, &
          title, crystal, tmp_dir, tupf, lgamma, noncolin, lspinorb, &
          lda_plus_u, &
          tfixed_occ_, tefield_, dipfield_, edir_, emaxpos_, eopreg_, eamp_, twfcollect  )

!       ==--------------------------------------------------------------==
!       ==  MAX DIMENSIONS                                              ==
!       ==--------------------------------------------------------------==

        CALL write_restart_xdim( ndw )

!       ==--------------------------------------------------------------==
!       ==  CELL & METRIC                                               ==
!       ==--------------------------------------------------------------==

        xhdum = 0.0d0
        CALL write_restart_cell( ndw, ibrav, celldm, ht_0%a, ht_m%a, &
          ht_m2%a, ht_0%hvel, vnhh, xnhh0, xnhhm, xhdum )

!       ==--------------------------------------------------------------==
!       ==  IONS                                                        ==
!       ==--------------------------------------------------------------==

        ALLOCATE( stau0(3, nat) )
        ALLOCATE( staum(3, nat) )
        ALLOCATE( svel0(3, nat) )
        ALLOCATE( svelm(3, nat) )
        ALLOCATE( tautmp(3, nat) )
        ALLOCATE( ityp(nat) )

        ia = 0
        DO i = 1, atoms_0%nsp
          DO j = 1, atoms_0%na(i)
            ia = ia + 1
            stau0(:,ia) = atoms_0%taus(:,ia)
            staum(:,ia) = atoms_m%taus(:,ia)
            svel0(:,ia) = atoms_0%vels(:,ia) 
            svelm(:,ia) = atoms_m%vels(:,ia)
            tautmp(:,ia) = taui(:,ia)
            ityp(ia) = atoms_0%ityp(ia)
          END DO
          mass(i) = atoms_0%m(i)
          ! ... WRITE( stdout,*) ' *** DEBUG write ions ', mass(i) ! debug
        END DO

        CALL get_nose_ions(XNOS0,XNOSM,XNOSM2,XNOSP) 

        tscal = .TRUE.
        CALL write_restart_ions( ndw, atom_label(1:ntyp), tscal, stau0, svel0, &
          staum, svelm, tautmp, atoms_0%for(1:3, 1:nat), cdmi, nat, ntyp, ityp, na, &
          mass, xnosp, xnos0, xnosm, xnosm2 )

        DEALLOCATE( stau0, staum, svel0, svelm, tautmp, ityp )

!       ==--------------------------------------------------------------==
!       ==  SYMMETRIES                                                  ==
!       ==--------------------------------------------------------------==

        CALL write_restart_symmetry( ndw )

!       ==--------------------------------------------------------------==
!       ==  PSEUDOPOTENTIALS                                            ==
!       ==--------------------------------------------------------------==

        DO i = 1, ntyp
          CALL write_restart_pseudo( ndw )
        END DO

!       ==--------------------------------------------------------------==
!       ==  OCCUPATION NUMBER                                           ==
!       ==--------------------------------------------------------------==

        CALL get_nose_electrons(XENOS0, XENOSM, XENOSM2, XENOSP) 

        tocc = .FALSE.
        tlam = .TRUE.
        teig = .FALSE.
        ldim = nbnd

        ikt = 0
        DO ispin = 1, cdesc%nspin
          DO ik = 1, cdesc%nkl
            ikt = ikt + 1
            ALLOCATE( occtmp(nbnd) )
            ALLOCATE( lambda(nbnd,nbnd) )
            lambda = 0.0d0
            DO i = 1, MIN( cdesc%nbt( ispin), nbnd )
              occtmp(i) = occ( i, ik, ispin )
              lambda(i,i) = 1.0d0
            END DO
            CALL write_restart_electrons( ndw, occtmp, occtmp, tocc, lambda, lambda, &
              ldim, tlam, nbnd, ispin, nspin, ikt, nkt, rnel, nelu_, neld_, xenosp, xenos0, &
              xenosm, xenosm2, ef, teig, occtmp, occtmp)
            DEALLOCATE( occtmp )
            DEALLOCATE( lambda )
          END DO
        END DO

!       ==--------------------------------------------------------------==
!       ==  G-Vectors                                                   ==
!       ==--------------------------------------------------------------==

        ALLOCATE( mill(3,ng) )
        mill = 0
        DO i = 1, gv%ng_l
          mill(:,gv%ig(i)) = gv%mill(:,i) 
        END DO
        CALL mp_sum( mill )
        tmill  = .TRUE.
        CALL write_restart_gvec( ndw, ng, gv%bi1, gv%bi2, gv%bi3, &
          gv%b1, gv%b2, gv%b3, tmill, mill )
        DEALLOCATE( mill )

!       ==--------------------------------------------------------------==
!       ==  (G+k)-Vectors                                               ==
!       ==--------------------------------------------------------------==

        ikt = 0
        DO ispin = 1, nspin
          DO i = 1, nk
            ikt = ikt + 1
            xk = kp%xk(:,i)
            wk = kp%weight(i)
            isk = ispin
            CALL write_restart_gkvec(ndw, ikt, nkt, ngwk_g(i), xk, wk, isk)
          END DO
        END DO

!       ==--------------------------------------------------------------==
!       ==  Tetrahedra                                                  ==
!       ==--------------------------------------------------------------==

        CALL write_restart_tetra( ndw )

!       ==--------------------------------------------------------------==
!       ==  CHARGE DENSITY AND POTENTIALS                               ==
!       ==--------------------------------------------------------------==

        trho = .FALSE.
        tv   = .FALSE.
        DO j = 1, nspin
          ALLOCATE( rhog(ng_l) )
          ALLOCATE( vg(ng_l) )
          CALL pfwfft( rhog(:), rho(:,:,:,j) )
          CALL pfwfft( vg(:)  , vpot(:,:,:,j) )
          CALL write_restart_charge(ndw, rhog, trho, vg, tv, ng_g, &
            j, nspin, gv%ig, ng_l )
          DEALLOCATE( rhog )
          DEALLOCATE( vg )
        END DO

!       ==--------------------------------------------------------------==
!       ==  WAVEFUNCTIONS                                               ==
!       ==--------------------------------------------------------------==

        tw0 = .TRUE.
        twm = .TRUE.
        wfc_scal_fpmd = SQRT(ht_0%omega)
        ikt = 0 
        IF( force_pairing ) THEN
          nspin_wfc = 1
        ELSE
          nspin_wfc = nspin
        END IF
        DO j = 1, nspin_wfc
          DO i = 1, nk
            ikt   = ikt + 1
            nb_g  = cdesc%nbt( j )
            ngw_g = cdesc%ngwt
            CALL write_restart_wfc( ndw, ikt, nkt, nkt, j, nspin_wfc, &
              wfc_scal_fpmd, c0(:,:,i,j), tw0, cm(:,:,i,j), twm, ngw_g, nb_g, &
              gv%ig, gv%ngw_l )
          END DO
        END DO

!       ==--------------------------------------------------------------==
!       ==  CLOSE RESTART FILE                                          ==
!       ==--------------------------------------------------------------==

        IF( ionode ) THEN 
          CLOSE ( ndw ) 
        END IF

        s1 = cclock()

!       ==--------------------------------------------------------------==
        IF( ionode ) THEN 
          WRITE( stdout,10) filename(1:INDEX(filename,' ')), (s1-s0)
   10     FORMAT(/,3X,'RESTART FILE WRITTEN ON UNIT ',A,' COMPLETED IN ',F8.3,' SEC.',/) 
        END IF 
!       ==--------------------------------------------------------------==

        RETURN 
        END SUBROUTINE writefile_fpmd_new

#endif

!=------------------------------------------------------------------


   SUBROUTINE writefile_fpmd( nfi, trutime, c0, cm, cdesc, occ, &
     atoms_0, atoms_m, acc, taui, cdmi, ibrav, celldm, &
     ht_m2, ht_m, ht_0, rho, desc, vpot, gv, kp)
                                                                        
        use electrons_module, only: nspin
        use cp_types, only: recvecs
        USE cell_module, only: boxdimensions, r_to_s
        USE brillouin, only: kpoints
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndw, gamma_only
        USE control_flags, ONLY: twfcollect, force_pairing
        USE atoms_type_module, ONLY: atoms_type
        USE io_global, ONLY: ionode, ionode_id
        USE io_global, ONLY: stdout
        USE charge_types, ONLY: charge_descriptor
        USE electrons_nose, ONLY: xnhe0, xnhem, xnhep, xnhem2, vnhe
        USE cell_nose, ONLY: xnhh0, xnhhm, xnhhp, vnhh
        USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, xnhpp, xnhpm2
        USE cp_restart, ONLY: cp_writefile, cp_write_wfc
                                                                        
        IMPLICIT NONE 
 
        INTEGER, INTENT(IN) :: nfi, ibrav
        COMPLEX(dbl), INTENT(IN) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(dbl), INTENT(IN) :: occ(:,:,:)
        TYPE (kpoints), INTENT(IN) :: kp 
        TYPE (boxdimensions), INTENT(IN) :: ht_m2, ht_m, ht_0
        TYPE (recvecs), INTENT(IN) :: gv
        TYPE (atoms_type), INTENT(IN) :: atoms_0, atoms_m
        REAL(dbl), INTENT(IN) :: rho(:,:,:,:)
        TYPE (charge_descriptor), INTENT(IN) :: desc
        TYPE (wave_descriptor) :: cdesc
        REAL(dbl), INTENT(INOUT) :: vpot(:,:,:,:)
                                                                        
        REAL(dbl), INTENT(IN) :: taui(:,:)
        REAL(dbl), INTENT(IN) :: celldm(:)
        REAL(dbl), INTENT(IN) :: acc(:), cdmi(:) 
        REAL(dbl), INTENT(IN) :: trutime

        INTEGER  :: nbnd
        REAL(dbl), ALLOCATABLE :: lambda(:,:)
        REAL(dbl) S0, S1
        REAL(dbl) :: ekincm
        REAL(dbl) :: mat_z(1,1,nspin)
             
        s0 = cclock()

        IF( ndw < 1 ) RETURN
        !
        !   this is used for benchmarking and debug
        !   if ndw < 1 Do not save wave functions and other system
        !   properties on the writefile subroutine

        nbnd    = MAXVAL( cdesc%nbt )
        ALLOCATE( lambda(nbnd,nbnd) )
        lambda  = 0.0d0
        ekincm = 0.0d0
        mat_z = 0.0d0

        CALL cp_writefile( ndw, ' ', .TRUE., nfi, trutime, acc, kp%nkpt, kp%xk, kp%weight, &
          ht_0%a, ht_m%a, ht_m2%a, ht_0%hvel, xnhh0, xnhhm, xnhhp, vnhh, taui, cdmi, &
          atoms_0%taus, atoms_0%vels, atoms_m%taus, atoms_m%vels, atoms_0%for, vnhp, &
          xnhp0, xnhpm, xnhpp, xnhpm2, occ, occ, lambda, lambda, gv%b1, gv%b2, gv%b3, &
          xnhe0, xnhem, xnhep, xnhem2, vnhe, ekincm, mat_z )

        DEALLOCATE( lambda )

        CALL cp_write_wfc( ndw, ' ', kp%nkpt, nspin, c0, cm )

        s1 = cclock()

!       ==--------------------------------------------------------------==
        IF( ionode ) THEN 
          WRITE( stdout,10) (s1-s0)
   10     FORMAT(/,3X,'RESTART FILE WRITTEN COMPLETED IN ',F8.3,' SEC.',/) 
        END IF 
!       ==--------------------------------------------------------------==

     RETURN 
   END SUBROUTINE writefile_fpmd



!=----------------------------------------------------------------------------=!

#ifdef __OLD_RESTART

        SUBROUTINE readfile_fpmd_new( nfi, trutime, &
          c0, cm, cdesc, occ, atoms_0, atoms_m, acc, taui, cdmi, ibrav, celldm, &
          ht_m2, ht_m, ht_0, rho, desc, vpot, gv, kp)
                                                                        
        use io_files, only: scradir
        use nose_ions, only: set_nose_ions 
        use nose_electrons, only: set_nose_electrons 
        use electrons_module, only: nel, nelt
        use cp_types, ONLY: recvecs
        USE cell_module, only: boxdimensions, cell_init, r_to_s, s_to_r
        USE brillouin, only: kpoints
        use parameters, only: npkx, nsx
        USE mp, ONLY: mp_sum, mp_barrier
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp_wave, ONLY: mergewf
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndr, tbeg, taurdr, gamma_only
        USE atoms_type_module, ONLY: atoms_type
        USE io_base, ONLY:  read_restart_header, read_restart_ions, &
          read_restart_cell, read_restart_electrons, &
          read_restart_gvec, read_restart_gkvec, read_restart_charge, &
          read_restart_wfc, read_restart_xdim, &
          read_restart_symmetry, read_restart_pseudo, read_restart_tetra
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE gvecw, ONLY: ecutwfc => ecutw
        USE gvecp, ONLY: ecutrho => ecutp
        USE fft, ONLY : pfwfft, pinvfft
        USE charge_types, ONLY: charge_descriptor
        USE ions_base, ONLY: nat, nsp, na
        USE electrons_module, ONLY: nspin
        USE control_flags, ONLY: twfcollect, force_pairing
        USE parser, ONLY: int_to_char
        USE wave_functions, ONLY: gram
        USE cell_nose, ONLY: vnhh, xnhh0, xnhhm
        USE grid_dimensions, ONLY: nr1, nr2, nr3
 
        IMPLICIT NONE 
 
        INTEGER, INTENT(INOUT) :: ibrav
        INTEGER, INTENT(OUT) :: nfi
        COMPLEX(dbl), INTENT(INOUT) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(dbl), INTENT(INOUT) :: occ(:,:,:)
        TYPE (kpoints), INTENT(INOUT) :: kp 
        TYPE (boxdimensions), INTENT(INOUT) :: ht_m2, ht_m, ht_0
        TYPE (recvecs), INTENT(INOUT) :: gv
        TYPE (atoms_type), INTENT(INOUT) :: atoms_0, atoms_m
        REAL(dbl), INTENT(INOUT) :: rho(:,:,:,:)
        TYPE (charge_descriptor), INTENT(IN) :: desc
        TYPE (wave_descriptor) :: cdesc
        REAL(dbl), INTENT(INOUT) :: vpot(:,:,:,:)
                                                                        
        REAL(dbl), INTENT(OUT) :: taui(:,:)
        REAL(dbl), INTENT(INOUT) :: celldm(:)
        REAL(dbl), INTENT(OUT) :: acc(:), cdmi(:) 
        REAL(dbl), INTENT(OUT) :: trutime

        COMPLEX(dbl), ALLOCATABLE :: rhog(:), vg(:)
        REAL(dbl) :: s0, s1

        !  variables, being read from file

        REAL(dbl) :: xk_(3), wk_
        INTEGER,   ALLOCATABLE :: igk_(:)
        REAL(dbl) :: b1_(3), b2_(3), b3_(3), bi1_(3), bi2_(3), bi3_(3)
        INTEGER,   ALLOCATABLE :: mill_(:,:)
        REAL(dbl), ALLOCATABLE :: stau0_(:,:), staum_(:,:), svel0_(:,:), svelm_(:,:), tautmp_(:,:)
        REAL(dbl) :: xnosm2_, xnosm_, xnos0_, xnosp_
        REAL(dbl) :: mass_( nsx )
        REAL(dbl) :: cdmi_( 3 )
        INTEGER  :: ngwk_g_( npkx ), ntyp_, na_( nsx ), nat_
        INTEGER ::  kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_, nkt_
        INTEGER  :: nbnd_, neld_, nelu_
        REAL(dbl) :: dgauss_, trutime_
        INTEGER :: ngauss_
        LOGICAL :: lgauss_
        INTEGER :: ntetra_
        LOGICAL :: ltetra_
        INTEGER :: natomwfc_
        LOGICAL :: doublegrid_
        REAL(dbl) :: gcutm_, gcuts_, dual_
        INTEGER :: modenum_
        REAL(dbl) :: alat_
        REAL(dbl) :: dum_
        LOGICAL :: lstres_ = .FALSE.
        LOGICAL :: lforce_ = .FALSE.
        CHARACTER(LEN=80) :: title_
        CHARACTER(LEN=80) :: crystal_
        CHARACTER(LEN=256) :: tmp_dir_
        REAL(dbl) :: celldm_(6)
        INTEGER :: ik_, nk_, ispin_, nspin_, isk_
        INTEGER :: nr1_, nr2_, nr3_, ibrav_, ngg_
        INTEGER :: nr1s_, nr2s_, nr3s_, ng_, nacc_
        INTEGER :: nfi_, ldim_
        REAL(dbl) :: rnel_
        REAL(dbl) :: ecutwfc_, ecutrho_, acc_( SIZE( acc ) )
        LOGICAL :: tupf_, lgamma_, lda_plus_u_
        LOGICAL :: noncolin_, lspinorb_
        REAL(dbl) :: xhdum(3,3)
        REAL(dbl) :: htvel_(3,3)
        REAL(dbl) :: hp0_(3,3), hm1_(3,3), hm2_(3,3)
        REAL(dbl), ALLOCATABLE :: occtmp_(:), lambda_(:,:)
        REAL(dbl) :: ef_
        REAL(dbl) :: xenosm2_, xenosm_, xenos0_, xenosp_ 
        INTEGER, ALLOCATABLE :: ityp_(:)
        CHARACTER(LEN=4) :: atom_label_( nsx )
        LOGICAL :: tscal_

        !  variables  that are not read from files

        INTEGER :: strlen
        CHARACTER(LEN=256) :: filename 
        LOGICAL :: tread

        REAL(dbl) :: wfc_scal, wfc_scal_fpmd
        LOGICAL :: tmill, tigl
        LOGICAL :: tocc, tlam, tw0, twm, trho, tv, teig


        INTEGER  :: ngw_g, nb_g, nk, nx, ns, nb_g_
        INTEGER  :: nxu, nxd, ispin, ng, ng_l, ng_g, maxig
        INTEGER  :: i, ik, ig, ip, ios, ierr, ia, j, ib
        REAL(dbl) :: rr1, rr2, rranf 
        EXTERNAL rranf

        LOGICAL :: tfixed_occ_, tefield_, dipfield_
        INTEGER :: edir_
        REAL(dbl) :: emaxpos_, eopreg_, eamp_
        INTEGER :: nspin_wfc

                                                                        
!       ==--------------------------------------------------------------==
!       OPEN THE APPROPRIATE FILE (FORT.NDW.ME) IN THE DIRECTORY CDIR     
!       ==--------------------------------------------------------------==

        CALL mp_barrier()
        s0 = cclock()

        ngw_g = cdesc%ngwt
        ng_l  = gv%ng_l
        ng_g  = gv%ng_g

!       ==--------------------------------------------------------------==
!       ==  ROOT PROCESSOR WRITE TO THE SCRATCH GENERAL VARIABLES                                  
!       ==--------------------------------------------------------------==

        IF ( ionode ) THEN 
          filename = 'fort.' // int_to_char( ndr ) 
          strlen  = index( scradir, ' ' ) - 1 
          IF ( strlen > 1 ) THEN 
            filename = scradir( 1 : strlen ) // '/' // filename 
          END IF 
          strlen  = index(filename,' ') - 1 
          OPEN( unit = ndr, file = filename(1:strlen), form = 'unformatted', status = 'old')
          REWIND (ndr)
          WRITE( stdout,10)
 10       FORMAT(/,3X,'READING FROM RESTART FILE ...')
        END IF

!       ==--------------------------------------------------------------==
!       ==  READ HEADER INFORMATION                                    ==
!       ==--------------------------------------------------------------==

        CALL read_restart_header( ndr, nfi_, trutime_, nr1_, nr2_, nr3_, &
          nr1s_, nr2s_, nr3s_, ngg_, nk_, ngwk_g_, nspin_, nbnd_, rnel_, &
          nelu_, neld_, nat_, ntyp_, na_, acc_, nacc_, ecutwfc_, ecutrho_, alat_, dum_, &
          kunit_, k1_, k2_, k3_, nk1_, nk2_, nk3_, dgauss_, ngauss_, lgauss_, ntetra_, ltetra_, & 
          natomwfc_, gcutm_, gcuts_, dual_, doublegrid_, modenum_, lstres_, lforce_, &
          title_, crystal_, tmp_dir_, tupf_, lgamma_, noncolin_, lspinorb_, &
          lda_plus_u_, & 
          tfixed_occ_, tefield_, dipfield_, edir_, emaxpos_, eopreg_, eamp_, twfcollect  )

        IF( lgamma_ .AND. .NOT. gamma_only )  &
          CALL errore(' readfile ', ' restart at gamma while using kpoints ', 1)
        IF( .NOT. lgamma_ .AND. gamma_only )  &
          CALL errore(' readfile ', ' restart with kpoints while using gamma ', 1)

        nfi     = nfi_
        trutime = trutime_

        IF( nspin /= nspin_ ) &
          CALL errore(' readfile ', ' restart spins differs ', 1)
        IF( nsp /= ntyp_ ) &
          CALL errore(' readfile ', ' atomic species differs ', 1)
        IF( nat /= nat_ ) &
          CALL errore(' readfile ', ' number of atoms differs ', 1)
        IF( nacc_ > SIZE( acc ) ) &
          CALL errore(' readfile ', ' invalid number of accumulators ', 1)

        acc( 1 : nacc_ ) = acc_ ( 1 : nacc_ )


!       ==--------------------------------------------------------------==
!       ==  MAX DIMENSIONS                                              ==
!       ==--------------------------------------------------------------==

        CALL read_restart_xdim( ndr )

!       ==--------------------------------------------------------------==
!       ==  CELL & METRIC                                               ==
!       ==--------------------------------------------------------------==

        CALL read_restart_cell( ndr, ibrav_, celldm_, hp0_, hm1_, &
          hm2_, htvel_, vnhh, xnhh0, xnhhm, xhdum )

        IF( ibrav_ /= ibrav ) THEN
          WRITE( stdout,fmt="(3X,'W: read_restart_cell, ibrav changed' )")
          WRITE( stdout,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ibrav_, ibrav
        END IF

        IF( ANY( celldm_ /= celldm ) ) THEN
          WRITE( stdout,fmt="(3X,'W: read_restart_cell, celldm changed' )")
          WRITE( stdout,fmt="(3X,'W: old = ',6F14.8 )") celldm_
          WRITE( stdout,fmt="(3X,'W: new = ',6F14.8 )") celldm
        END IF

        IF( .NOT. tbeg ) THEN
          CALL cell_init( ht_0, hp0_ )
          IF( ALL( hm2_ == 0.0d0 ) ) THEN
            CALL cell_init( ht_m, hp0_ )
          ELSE
            CALL cell_init( ht_m, hm1_ )
          END IF
          IF( ALL( hm2_ == 0.0d0 ) ) THEN
            CALL cell_init( ht_m2, hp0_ )
          ELSE
            CALL cell_init( ht_m2, hm2_ )
          END IF
          ht_0%hvel = htvel_  !  set cell velocity
        END IF

!       ==--------------------------------------------------------------==
!       ==  IONS                                                        ==
!       ==--------------------------------------------------------------==

        ALLOCATE( stau0_(3, nat) )
        ALLOCATE( staum_(3, nat) )
        ALLOCATE( svel0_(3, nat) )
        ALLOCATE( svelm_(3, nat) )
        ALLOCATE( tautmp_(3, nat) )
        ALLOCATE( ityp_( nat ) )

        CALL read_restart_ions(ndr, atom_label_, tscal_, stau0_, svel0_, staum_, &
          svelm_, tautmp_, atoms_0%for(1:3, 1:nat), cdmi_, nat_, ntyp_, ityp_, na_, mass_, &
          xnosp_, xnos0_, xnosm_, xnosm2_ )

        DO i = 1, ntyp_
          IF( na_( i ) /= atoms_0%na( i ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_ions, Number of atoms chenged' )")
            WRITE( stdout,fmt="(3X,'W: is = ',I3,' old = ',I10,' new = ',I10 )") i, na_(i), atoms_0%na(i)
          END IF
        END DO

        cdmi( 1:3 ) = cdmi_ ( 1:3 )

        ia = 0
        DO i = 1, atoms_0%nsp
          DO j = 1, atoms_0%na(i)
            ia = ia + 1
            atoms_0%ityp(ia) = ityp_(ia)
            atoms_m%ityp(ia) = ityp_(ia)
            IF( tscal_ ) THEN
              IF( .NOT. taurdr ) THEN
                atoms_0%taus(:,ia) = stau0_(:,ia)
                atoms_m%taus(:,ia) = staum_(:,ia)
                CALL s_to_r( stau0_(:,ia), atoms_0%taur(:,ia), ht_0 )
                CALL s_to_r( staum_(:,ia), atoms_m%taur(:,ia), ht_m )
              END IF
              atoms_0%vels(:,ia) = svel0_(:,ia)
              atoms_m%vels(:,ia) = svelm_(:,ia)
            ELSE
              IF( .NOT. taurdr ) THEN
                CALL r_to_s( stau0_(:,ia), atoms_0%taus(:,ia), ht_0 )
                CALL r_to_s( staum_(:,ia), atoms_m%taus(:,ia), ht_m )
                atoms_0%taur(:,ia) = stau0_(:,ia)
                atoms_m%taur(:,ia) = staum_(:,ia)
              END IF
              CALL r_to_s( svel0_(:,ia), atoms_0%vels(:,ia), ht_0 )
              CALL r_to_s( svelm_(:,ia), atoms_m%vels(:,ia), ht_m )
            END IF
            taui(:,ia) = tautmp_(:,ia)
          END DO
        END DO

        CALL set_nose_ions( xnos0_ , xnosm_ , xnosm2_ , xnosp_ ) 

        DEALLOCATE( stau0_ , staum_ , svel0_ , svelm_ , tautmp_ , ityp_ )

!       ==--------------------------------------------------------------==
!       ==  SYMMETRIES                                                  ==
!       ==--------------------------------------------------------------==

        CALL read_restart_symmetry( ndr )

!       ==--------------------------------------------------------------==
!       ==  PSEUDOPOTENTIALS                                            ==
!       ==--------------------------------------------------------------==

        DO i = 1, ntyp_
          CALL read_restart_pseudo( ndr )
        END DO

!       ==--------------------------------------------------------------==
!       ==  OCCUPATION NUMBER                                           ==
!       ==--------------------------------------------------------------==


        tocc = .FALSE.
        tlam = .FALSE.
        teig = .FALSE.

        DO ispin = 1, cdesc%nspin

          DO ik = 1, cdesc%nkl

            ALLOCATE( occtmp_ ( MAXVAL( cdesc%nbt ) ) )
            ALLOCATE( lambda_ ( 1, 1 ) )

            CALL read_restart_electrons( ndr, occtmp_, occtmp_, tocc, lambda_, lambda_, &
              ldim_, tlam, nbnd_, ispin_, nspin_, ik_, nkt_, rnel_, nelu_, neld_, xenosp_, xenos0_, &
              xenosm_, xenosm2_, ef_, teig, occtmp_, occtmp_ )

            IF( tocc ) THEN
              DO i = 1, MIN( cdesc%nbt( ispin ), nbnd_ )
                occ( i, ik, ispin ) = occtmp_( i )
              END DO
            END IF

            DEALLOCATE( occtmp_ )
            DEALLOCATE( lambda_ )

          END DO

        END DO

        CALL set_nose_electrons( xenos0_, xenosm_, xenosm2_, xenosp_ ) 


!       ==--------------------------------------------------------------==
!       ==  G-Vectors                                                   ==
!       ==--------------------------------------------------------------==

        ALLOCATE( mill_(3, gv%ng_g) )

        ! maxig = MAXVAL( gv%ig )
        ! IF( ( maxig > ng ) .OR. ( gv%ng_l > ng ) .OR. ( gv%ng_l > SIZE( gv%mill, 2 ) ) ) &
        !   CALL errore(' readfile ',' wrong G vector index array ', 1)

        tmill = .FALSE.

        CALL read_restart_gvec( ndr, ng_, bi1_, bi2_, bi3_, &
          b1_, b2_, b3_, tmill, mill_ )

          IF( ng_ /= gv%ng_g ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, ng changed' )")
            WRITE( stdout,fmt="(3X,'W: old = ',I10,' new = ',I10 )") ng_, gv%ng_g
          END IF
          IF( ANY( bi1_ /= gv%bi1 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, bi1 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi1_, gv%bi1
          END IF
          IF( ANY( bi2_ /= gv%bi2 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, bi2 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi2_, gv%bi2
          END IF
          IF( ANY( bi3_ /= gv%bi3 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, bi3 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") bi3_, gv%bi3
          END IF
          IF( ANY( b1_ /= gv%b1 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, b1 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b1_, gv%b1
          END IF
          IF( ANY( b2_ /= gv%b2 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, b2 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b2_, gv%b2
          END IF
          IF( ANY( b3_ /= gv%b3 ) ) THEN
            WRITE( stdout,fmt="(3X,'W: read_restart_gvec, b3 changed' )")
            WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") b3_, gv%b3
          END IF

        DEALLOCATE( mill_ )

!       ==--------------------------------------------------------------==
!       ==  (G+k)-Vectors                                               ==
!       ==--------------------------------------------------------------==

        DO ispin = 1, nspin

          DO i = 1, kp%nkpt

            CALL read_restart_gkvec(ndr, ik_, nk_, ngwk_g_( i ), &
              xk_, wk_, isk_)

            IF( ANY( xk_ /= kp%xk(:,i) ) ) THEN
              WRITE( stdout,fmt="(3X,'W: read_restart_gkvec, xk changed' )")
              WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") xk_, kp%xk(:,i)
            END IF

            IF( wk_ /= kp%weight(i) ) THEN
              WRITE( stdout,fmt="(3X,'W: read_restart_gkvec, wk changed' )")
              WRITE( stdout,fmt="(3X,'W: old, new = ', 3F10.4, 2X, 3F10.4 )") wk_, kp%weight(i)
            END IF

          END DO
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
        tread = .FALSE.
        DO j = 1, nspin
          IF( tread ) THEN
            ALLOCATE( rhog(ng_l), vg(ng_l) )
            CALL read_restart_charge(ndr, rhog, trho, vg, tv, ng_g, &
              ispin_, nspin_, gv%ig, ng_l )
            CALL pinvfft( vpot(:,:,:,j), rhog(:) )
            CALL pinvfft( rho(:,:,:,j), vg(:) )
            DEALLOCATE( rhog, vg )
          ELSE
            CALL read_restart_charge( ndr )
          END IF
        END DO

!       ==--------------------------------------------------------------==
!       ==  WAVEFUNCTIONS                                               ==
!       ==--------------------------------------------------------------==

        wfc_scal_fpmd = 1.0d0 / SQRT(ht_0%omega)
        IF( force_pairing ) THEN
          nspin_wfc = 1
        ELSE
          nspin_wfc = nspin
        END IF
        DO j = 1, nspin_wfc
          DO i = 1, cdesc%nkt
            tw0 = .TRUE.
            twm = .TRUE.
            c0(:,:,i,j) = 0.0d0
            ik_  = i
            nk_  = cdesc%nkt
            nb_g = cdesc%nbt( j )
            CALL read_restart_wfc(ndr, ik_, nk_, nk_, ispin_, nspin_, &
              wfc_scal, c0(:,:,i,j), tw0, cm(:,:,i,j), twm, ngw_g, nb_g_, gv%ig, gv%ngw_l )
            IF( tw0 .AND. .NOT. twm ) THEN
              cm(:,:,i,j) = c0(:,:,i,j)
            END IF
            IF( wfc_scal == 1.0d0 ) THEN
              c0(:,:,i,j) = c0(:,:,i,j) * wfc_scal_fpmd
              cm(:,:,i,j) = cm(:,:,i,j) * wfc_scal_fpmd
            END IF
            IF( nb_g > nb_g_ ) THEN
              DO ib = nb_g_, nb_g
                DO ig = 1, SIZE( c0, 1 )
                  rr1 = 0.5d0 - rranf()
                  rr2 = rranf()
                  c0( ig, ib, i, j ) =  DCMPLX(rr1, rr2)
                END DO
              END DO
              CALL gram( j, c0( :, :, i, j ), cdesc )
              cm( :, nb_g_ : nb_g, i, j ) = c0( :, nb_g_ : nb_g, i, j )
            END IF
            ! WRITE( *, * ) 'DEBUG: ', SUM(  c0(:,:,i,j) ), SUM( cm(:,:,i,j) ) ! DEBUG
          END DO
        END DO

!       ==--------------------------------------------------------------==
!       ==  CLOSE RESTART FILE                                          ==
!       ==--------------------------------------------------------------==

        IF( ionode ) THEN 
          CLOSE ( ndr ) 
        END IF

        CALL mp_barrier()
        s1 = cclock()

!       ==--------------------------------------------------------------==
        IF( ionode ) THEN 
          WRITE( stdout,20) filename(1:INDEX(filename,' ')), (s1-s0)
   20     FORMAT(3X,'DISK READ FROM UNIT ',A,' COMPLETED IN ',F8.3,' SEC.',/) 
        END IF 
!       ==--------------------------------------------------------------==

        RETURN 
        END SUBROUTINE readfile_fpmd_new

#endif

!=----------------------------------------------------------------------------=!

        SUBROUTINE readfile_fpmd( nfi, trutime, &
          c0, cm, cdesc, occ, atoms_0, atoms_m, acc, taui, cdmi, ibrav, celldm, &
          ht_m2, ht_m, ht_0, rho, desc, vpot, gv, kp)
                                                                        
        use electrons_module, only: nel, nelt
        use cp_types, ONLY: recvecs
        USE cell_module, only: boxdimensions, cell_init, r_to_s, s_to_r
        USE brillouin, only: kpoints
        use parameters, only: npkx, nsx
        USE mp, ONLY: mp_sum, mp_barrier
        USE mp_global, ONLY: mpime, nproc, group, root
        USE mp_wave, ONLY: mergewf
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndr, tbeg, taurdr, gamma_only
        USE atoms_type_module, ONLY: atoms_type
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE gvecw, ONLY: ecutwfc => ecutw
        USE gvecp, ONLY: ecutrho => ecutp
        USE fft, ONLY : pfwfft, pinvfft
        USE charge_types, ONLY: charge_descriptor
        USE ions_base, ONLY: nat, nsp, na
        USE electrons_module, ONLY: nspin
        USE control_flags, ONLY: twfcollect, force_pairing
        USE wave_functions, ONLY: gram
        USE grid_dimensions, ONLY: nr1, nr2, nr3
        USE electrons_nose, ONLY: xnhe0, xnhem, xnhep, xnhem2, vnhe
        USE cell_nose, ONLY: xnhh0, xnhhm, xnhhp, vnhh
        USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, xnhpp, xnhpm2
        USE cp_restart, ONLY: cp_readfile, cp_read_wfc
 
        IMPLICIT NONE 
 
        INTEGER, INTENT(INOUT) :: ibrav
        INTEGER, INTENT(OUT) :: nfi
        COMPLEX(dbl), INTENT(INOUT) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(dbl), INTENT(INOUT) :: occ(:,:,:)
        TYPE (kpoints), INTENT(INOUT) :: kp 
        TYPE (boxdimensions), INTENT(INOUT) :: ht_m2, ht_m, ht_0
        TYPE (recvecs), INTENT(INOUT) :: gv
        TYPE (atoms_type), INTENT(INOUT) :: atoms_0, atoms_m
        REAL(dbl), INTENT(INOUT) :: rho(:,:,:,:)
        TYPE (charge_descriptor), INTENT(IN) :: desc
        TYPE (wave_descriptor) :: cdesc
        REAL(dbl), INTENT(INOUT) :: vpot(:,:,:,:)
                                                                        
        REAL(dbl), INTENT(OUT) :: taui(:,:)
        REAL(dbl), INTENT(INOUT) :: celldm(:)
        REAL(dbl), INTENT(OUT) :: acc(:), cdmi(:) 
        REAL(dbl), INTENT(OUT) :: trutime


        REAL(dbl) :: s0, s1
        REAL(dbl), ALLOCATABLE :: lambda_ ( : , : )
        INTEGER :: nbnd_
        REAL(dbl) :: ekincm
        REAL(dbl) :: mat_z_(1,1,nspin)
        LOGICAL :: tens = .FALSE.

        CALL mp_barrier()
        s0 = cclock()

        nbnd_   = MAXVAL( cdesc%nbt )
        ALLOCATE( lambda_( nbnd_ , nbnd_ ) )
        lambda_  = 0.0d0

        CALL cp_readfile( ndr, ' ', .TRUE., nfi, trutime, acc, kp%nkpt, kp%xk, kp%weight, &
          ht_0%a, ht_m%a, ht_m2%a, ht_0%hvel, xnhh0, xnhhm, xnhhp, vnhh, taui, cdmi, &
          atoms_0%taus, atoms_0%vels, atoms_m%taus, atoms_m%vels, atoms_0%for, vnhp, &
          xnhp0, xnhpm, xnhpp, xnhpm2, occ, occ, lambda_ , lambda_ , gv%b1, gv%b2,   &
          gv%b3, xnhe0, xnhem, xnhep, xnhem2, vnhe, ekincm, mat_z_ , tens )

        DEALLOCATE( lambda_ )

        CALL cp_read_wfc( ndr, ' ', kp%nkpt, nspin, c0, cm )

        CALL mp_barrier()
        s1 = cclock()

!       ==--------------------------------------------------------------==
        IF( ionode ) THEN 
          WRITE( stdout,20)  (s1-s0)
   20     FORMAT(3X,'DISK READ COMPLETED IN ',F8.3,' SEC.',/) 
        END IF 
!       ==--------------------------------------------------------------==

        RETURN 
        END SUBROUTINE readfile_fpmd


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

          IF ( ionode ) THEN
            filename = 'fort.' // int_to_char( ndr )
            strlen  = INDEX( scradir, ' ' ) - 1
            filename = scradir( 1 : strlen ) // '/' // filename
            INQUIRE( FILE = TRIM( filename ), EXIST = lval )
            ! WRITE(6,*) '  checking file ', lval, ' ', TRIM( filename )
          END IF
          CALL mp_bcast( lval, ionode_id )
          check_restartfile = lval
          RETURN 
        END FUNCTION

!=----------------------------------------------------------------------------=!
     END MODULE restart_file
!=----------------------------------------------------------------------------=!
