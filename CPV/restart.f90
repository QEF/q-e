MODULE restart
  IMPLICIT NONE
  SAVE

CONTAINS

!-----------------------------------------------------------------------
      subroutine writefile_new                                         &
     &     ( ndw,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      use elct, only: n, nx, ngw, ng0, nspin, nel, ngw_g
      use ions_module, only: nsp, na, nax
      use parm, ONLY: nr1, nr2, nr3
      use gvec, ONLY: ng, ngl, mill_g, ng_g, mill_l, bi1, bi2, bi3, ig_l2g
      use io_base, only: write_restart_header, write_restart_ions, &
          write_restart_cell, write_restart_electrons, &
          write_restart_gvec, write_restart_gkvec, write_restart_charge, &
          write_restart_wfc, write_restart_symmetry, &
          write_restart_xdim, write_restart_pseudo
      use mp, only: mp_sum
      use mp_global
      use io_global
      use cell_module, only: boxdimensions, s_to_r, cell_init
      USE ncprm, ONLY: r, rab
!
      implicit none
      integer :: ndw, nfi
      real(kind=8) :: h(3,3), hold(3,3)
      complex(kind=8) :: c0(:,:), cm(:,:)
      real(kind=8) :: tausm(:,:,:),taus(:,:,:), fion(:,:,:)
      real(kind=8) :: vels(:,:,:), velsm(:,:,:)
      real(kind=8) :: acc(:),lambda(:,:), lambdam(:,:)
      real(kind=8) :: xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8), INTENT(in) :: ecut, ecutw, delt
      real(kind=8), INTENT(in) :: pmass(:)
      real(kind=8), INTENT(in) :: celldm(:)
      integer, INTENT(in) :: ibrav
      integer :: nbeg = 0
      integer :: nk = 1
      integer :: ngwkl(1), ngwkg(1), nbnd, nelt, nelu, neld, ntyp, nb_g
      integer :: nat = 0
      integer :: nacx = 10
      real(kind=8) :: trutime = 0.0d0
      real(kind=8) :: ecutwfc, ecutrho

      REAL(dbl), ALLOCATABLE :: stau0(:,:), staum(:,:), svel0(:,:), svelm(:,:), tautmp(:,:)
      REAL(dbl), ALLOCATABLE :: fiontmp(:,:)
      type (boxdimensions) :: box
      real(dbl) :: ht(3,3), htvel(3,3), ht0(3,3), htm(3,3), htm2(3,3)
      real(dbl) :: xdum
      real(dbl) :: hdum(3,3)
      real(dbl) :: cdmi(3)
      real(dbl) :: mass(nsp)
      real(dbl), allocatable :: occ(:), occm(:), lamtmp(:,:), lamtmpm(:,:), eigtmp(:)
      LOGICAL :: tocc, tlam, trho, tv, tw0, twm
      integer, allocatable :: mill(:,:), igk(:)
      real(dbl) :: xk(3), wk
      complex(kind=8), allocatable :: rhog(:), vg(:)

      LOGICAL :: lstres, lforce
      character(len=80) :: title, crystal, tmp_dir
!
      integer i, ia, is, j, ispin, ik
      integer :: strlen
      character(len=80) :: filename

      LOGICAL :: tovrw = .FALSE.
      INTEGER ::  k1, k2, k3, nk1, nk2, nk3
      REAL(dbl) :: dgauss
      INTEGER :: ngauss
      INTEGER :: ntetra
      INTEGER :: natomwfc
      LOGICAL :: doublegrid, tupf
      REAL(dbl) :: gcutm, gcuts, dual
      INTEGER :: modenum, kunit
      REAL(dbl) :: alat
      REAL(dbl) :: ef, rnel, wfc_scal_cp90
      character(len=4) :: atom_label(nsp)

      LOGICAL :: twrite
      LOGICAL :: tscal
      LOGICAL :: teig
      LOGICAL :: tmill
      LOGICAL :: lgauss
      LOGICAL :: ltetra
      LOGICAL :: lgamma
      INTEGER, ALLOCATABLE :: ityp(:)
      INTEGER :: isk, tetra(4)
      REAL(dbl) :: zmesh_, xmin_, dx_
      REAL(dbl) :: ainv(3,3), deth

!
! Only the first node writes
!
      if (ionode) then
       !  open (unit=ndw,status='unknown',form='unformatted')
         CALL cpitoa(ndw, filename)
         filename = 'fort.'//filename
         strlen  = index(filename,' ') - 1 
         OPEN(unit=ndw, file=filename(1:strlen), form='unformatted', status='unknown')
         REWIND ndw
      end if

      ht = TRANSPOSE(h)
      call cell_init(box,ht)

!       ==--------------------------------------------------------------==
!       ==  WRITE HEADER INFORMATIONS                                   ==
!       ==--------------------------------------------------------------==

      ngwkg(1) = ngw_g
      ngwkl(1) = ngw 
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
      title = ''
      crystal = ''
      tmp_dir = ''
      kunit = 1
      lgauss = .FALSE.
      ltetra = .FALSE.
      twrite = .TRUE.
      tupf   = .TRUE.
      lgamma = .TRUE.
      CALL write_restart_header(ndw, twrite, nfi, trutime, nbeg, nr1, nr2, nr3, &
        nr1, nr2, nr3, ng, ng_g, nk, nk, ngwkl, ngwkg, nspin, nbnd, rnel, nelu, &
        neld, nat, ntyp, na, acc, nacx, ecutwfc, ecutrho, alat, ekincm, &
        kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
        natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, &
        title, crystal, tmp_dir, tupf, lgamma)

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
      twrite = .TRUE. 
      CALL write_restart_cell( ndw, twrite, ibrav, celldm, ht, htm, &
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
          stau0(:,ia) = taus(:,j,i)
          staum(:,ia) = tausm(:,j,i)
          svel0(:,ia) = vels(:,j,i)
          svelm(:,ia) = velsm(:,j,i)
          CALL s_to_r( taus(:,j,i), tautmp(:,ia), box )
          fiontmp(:,ia) = fion(:,j,i)
          ityp(ia) = i
        END DO
        mass(i) = pmass(i)
      END DO

      xdum = 0.0d0
      cdmi = 0.0d0
      tscal = .TRUE.
      twrite = .TRUE. 
      CALL write_restart_ions(ndw, twrite, atom_label, tscal, stau0, svel0, &
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
!       CALL write_restart_pseudo( ndw, twrite, &
!         zmesh_, xmin_, dx_, r(:,i), rab(:,i), vnl(:,:,i), chi(:,:,i), oc(:,i), &
!         rho_at(:,i), rho_atc(:,i), mesh(i), msh(i), nchi(i), lchi(:,i), &
!         numeric(i), cc(:,i), alpc(:,i), zp(i), aps(:,:,i), alps(:,:,i), &
!         zv(i), nlc(i), nnl(i), lmax(i), lloc(i), bhstype(i), dion(:,:,i), &
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
      twrite = .TRUE.
      CALL write_restart_electrons( ndw, twrite, occ, occ, tocc, lambda, lambdam, &
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
      twrite = .TRUE.
      CALL write_restart_gvec( ndw, twrite, ng_g, bi1, bi2, bi3, &
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
        tetra = 0.0d0
        isk = 1
        twrite = .TRUE.
        CALL write_restart_gkvec(ndw, twrite, i, nk, ngwkg(i), xk, wk, tetra, isk)
      END DO

!     ==--------------------------------------------------------------==
!     ==  CHARGE DENSITY AND POTENTIALS                               ==
!     ==--------------------------------------------------------------==

      trho = .FALSE.
      tv   = .FALSE.
      twrite = .TRUE.
      DO j = 1, nspin
        ALLOCATE( rhog(ng), vg(ng) )
!        CALL fft_initialize
!        CALL pfwfft( rhog(:,i), rho(i)%r(:,:,:) )
!        CALL pfwfft( vg(:,i)  , vpot(:,:,:,i) )
        CALL write_restart_charge(ndw, twrite, rhog, trho, vg, tv, ng_g, &
          j, nspin, ig_l2g, ng )
        DEALLOCATE( rhog, vg )
      END DO


!     ==--------------------------------------------------------------==
!     ==  WAVEFUNCTIONS                                               ==
!     ==--------------------------------------------------------------==

      tw0 = .TRUE.
      twm = .TRUE.
      twrite = .TRUE.
        
      call invmat3(h,ainv,deth)
      wfc_scal_cp90 = 1.0d0 / SQRT(ABS(deth))
      DO j = 1, nspin
        DO i = 1, nk
          nb_g = nx
          CALL write_restart_wfc(ndw, twrite, i, nk, kunit, j, nspin, &
            wfc_scal_cp90, c0, tw0, cm, twm, ngw_g, nb_g, ig_l2g, ngw )
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
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      use elct, only: n, nx, ngw, ng0, nspin, nel, ngw_g
      use ions_module, only: nsp, na, nax
      use parm, ONLY: nr1, nr2, nr3
      use gvec, ONLY: ng, ngl, mill_g, ng_g, mill_l, bi1, bi2, bi3, ig_l2g
      use io_base, only: read_restart_header, read_restart_ions, &
          read_restart_cell, read_restart_electrons, &
          read_restart_gvec, read_restart_gkvec, read_restart_charge, &
          read_restart_wfc, read_restart_xdim, &
          read_restart_symmetry, read_restart_pseudo
      use mp, only: mp_sum
      use mp_global
      use io_global
      use cell_module, only: boxdimensions, s_to_r, cell_init, r_to_s
!
      implicit none
      integer :: ndr, nfi, flag
      real(kind=8) :: h(3,3), hold(3,3)
      complex(kind=8) :: c0(:,:), cm(:,:)
      real(kind=8) :: tausm(:,:,:),taus(:,:,:), fion(:,:,:)
      real(kind=8) :: vels(:,:,:), velsm(:,:,:)
      real(kind=8) :: acc(:),lambda(:,:), lambdam(:,:)
      real(kind=8) :: xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm
      real(kind=8) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8), INTENT(in) :: ecut, ecutw, delt
      real(kind=8), INTENT(in) :: pmass(:)
      real(kind=8), INTENT(in) :: celldm(6)
      integer, INTENT(in) :: ibrav
      integer :: nbeg = 0
      integer :: nk = 1
      integer :: ngwkg(1), ngwkl(1), nbnd, nelt, nelu, neld, ntyp, nb_g
      integer :: nat = 0
      integer :: nacx = 10
      real(kind=8) :: trutime = 0.0d0
      real(kind=8) :: ecutwfc, ecutrho

      REAL(dbl), ALLOCATABLE :: stau0(:,:), staum(:,:), svel0(:,:), svelm(:,:), tautmp(:,:)
      REAL(dbl), ALLOCATABLE :: fiontmp(:,:)
      type (boxdimensions) :: box, boxm
      real(dbl) :: ht(3,3), htvel(3,3), ht0(3,3), htm(3,3), htm2(3,3)
      real(dbl) :: xdum
      real(dbl) :: hdum(3,3)
      real(dbl) :: cdmi(3)
      real(dbl) :: mass(nsp)
      real(dbl), allocatable :: occ(:), eigtmp(:)
      LOGICAL :: tocc, tlam, trho, tv, tw0, twm
      integer, allocatable :: mill(:,:), igk(:)
      real(dbl) :: xk(3), wk
      complex(kind=8), allocatable :: rhog(:), vg(:)

      LOGICAL :: tovrw = .FALSE.
      INTEGER ::  k1, k2, k3, nk1, nk2, nk3
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
      REAL(dbl) :: ef, rnel
      LOGICAL :: lstres, lforce
      character(len=80) :: title, crystal, tmp_dir
      character(len=4) :: atom_label(nsp)
!
      integer :: i, ia, is, j
      integer :: nfi_, ik_, nk_, ispin_, nspin_, isk_, tetra_(4)
      real(kind=8) :: acc_(10), celldm_(6)
      real(kind=8) :: vnhp_, xnhp0_, xnhpm_
      integer :: strlen, ibrav_, kunit
      character(len=80) :: filename
      LOGICAL :: tread
      LOGICAL :: tscal
      LOGICAL :: tmill, tigl, lgamma
      LOGICAL :: teig, tupf
      INTEGER, ALLOCATABLE :: ityp(:)
      REAL(dbl) :: wfc_scal, wfc_scal_cp90
!
! Only the first node read
!
      if (ionode) then
         ! open (unit=ndr, status='old', form='unformatted')
         CALL cpitoa(ndr, filename)
         filename = 'fort.'//filename 
         strlen  = index(filename,' ') - 1 
         OPEN(unit=ndr, file=filename(1:strlen), form='unformatted', status='old')
         REWIND (ndr)
         WRITE(6,10)
 10      FORMAT(/,3X,'READING FROM RESTART FILE ...')
      end if

      if (flag.eq.-1) then
         write(6,'((a,i3,a))') ' ### reading from file ',ndr,' only h  ##'
      else if (flag.eq.0) then
         write(6,'((a,i3,a))') ' ## reading from file ',ndr,' only c0  ##'
      else
         write(6,'((a,i3))') ' ## reading from file ',ndr
      end if

!     ==--------------------------------------------------------------==
!     ==  READ HEADER INFORMATIONS                                   ==
!     ==--------------------------------------------------------------==

      ngwkl(1) = ngw 
      ngwkg(1) = ngw_g
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
      kunit = nk
      ecutwfc = ecutw
      ecutrho = ecut
      tread = .TRUE.
      tovrw = .FALSE.
      CALL read_restart_header(ndr, tovrw, tread, nfi_, trutime, nbeg, nr1, nr2, nr3, &
        nr1, nr2, nr3, ng, ng_g, nk, nk, ngwkl, ngwkg, nspin, nbnd, rnel, nelu, neld, &
        nat, ntyp, na, acc_, nacx, ecutwfc, ecutrho, alat, ekincm, &
        kunit, k1, k2, k3, nk1, nk2, nk3, dgauss, ngauss, lgauss, ntetra, ltetra, &
        natomwfc, gcutm, gcuts, dual, doublegrid, modenum, lstres, lforce, &
        title, crystal, tmp_dir, tupf, lgamma)

!      if( .not. lgamma ) &
!        call errore(' readfile ',' restart contains a system not at gamma ',1)

      if (flag > -1) then
        nfi = nfi_
        acc = acc_
      end if

!     ==--------------------------------------------------------------==
!     ==  MAX DIMENSIONS                                              ==
!     ==--------------------------------------------------------------==

      CALL read_restart_xdim( ndr )
      
!     ==--------------------------------------------------------------==
!     ==  CELL & METRIC                                               ==
!     ==--------------------------------------------------------------==

      hdum = 0.0d0
      htm2 = 0.0d0
        
      tovrw = .FALSE.
      tread = .TRUE.
      celldm_ = celldm
      CALL read_restart_cell( ndr, tovrw, tread, ibrav_, celldm_, ht, htm, &
        htm2, htvel, vnhh, xnhh0, xnhhm, hdum)

      h = TRANSPOSE(ht) 
      hold = TRANSPOSE(htm) 
      velh = TRANSPOSE(htvel) 

      CALL cell_init(box, ht)
      CALL cell_init(boxm, htm)

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

        DO i = 1, nsp
          mass(i) = pmass(i)
        END DO
        xdum = 0.0d0
        cdmi = 0.0d0
        tovrw = .FALSE.
        tread = .TRUE.
        CALL read_restart_ions(ndr, tovrw, tread, atom_label, tscal, stau0, svel0, &
          staum, svelm, tautmp, fiontmp, cdmi, nat, ntyp, ityp, na, mass, vnhp_,   &
          xnhp0_, xnhpm_, xdum)

        IF( flag > 0 ) THEN
          ia = 0
          DO i = 1, nsp
            DO j = 1, na(i)
              ia = ia + 1
              IF( tscal ) THEN
                taus(:,j,i) = stau0(:,ia)
                tausm(:,j,i) = staum(:,ia)
                vels(:,j,i) = svel0(:,ia)
                velsm(:,j,i) = svelm(:,ia)
              ELSE
                CALL r_to_s( stau0(:,ia), taus(:,j,i), box )
                CALL r_to_s( staum(:,ia), tausm(:,j,i), boxm )
                CALL r_to_s( svel0(:,ia), vels(:,j,i), box )
                CALL r_to_s( svelm(:,ia), velsm(:,j,i), boxm )
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

        ALLOCATE( occ(nbnd), eigtmp(nbnd) )

        occ = 0.0d0
        IF( flag > 0 ) THEN
          tocc = .FALSE.
          tlam = .TRUE.
        ELSE
          tocc = .FALSE.
          tlam = .FALSE.
        END IF
        tovrw = .FALSE.
        teig = .FALSE.
        tread = .TRUE.
        CALL read_restart_electrons( ndr, tovrw, tread, occ, occ, tocc, lambda, &
          lambdam, nx, tlam, nbnd, ispin_, nspin, ik_, nk, rnel, nelu, neld,    &
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
        tread = .TRUE.
        tovrw = .FALSE.
        tmill = .FALSE.
        CALL read_restart_gvec( ndr, tovrw, tread, ng_g, bi1, bi2, bi3,  &
          bi1, bi2, bi3, tmill, mill )
        DEALLOCATE( mill )

!       ==--------------------------------------------------------------==
!       ==  (G+k)-Vectors                                               ==
!       ==--------------------------------------------------------------==

        DO i = 1, nk
          xk(1) = 0.0d0
          xk(2) = 0.0d0
          xk(3) = 0.0d0
          wk = 1.0d0
          tread = .TRUE.
          tovrw = .FALSE.
          CALL read_restart_gkvec(ndr, tovrw, tread, ik_, nk_, ngwkg(i), &
            xk, wk, tetra_, isk_)
        END DO

!       ==--------------------------------------------------------------==
!       ==  CHARGE DENSITY AND POTENTIALS                               ==
!       ==--------------------------------------------------------------==

        trho = .FALSE.
        tv   = .FALSE.
        DO j = 1, nspin
          ALLOCATE( rhog(ng), vg(ng) )
          tread = .TRUE.
          tovrw = .FALSE.
          CALL read_restart_charge(ndr, tovrw, tread, rhog, trho, vg, tv, ng_g, &
            ispin_, nspin, ig_l2g, ng )
!          CALL fft_initialize
!          CALL pinvfft( vpot(:,:,:,i), rhog(:,i) )
!          CALL pinvfft( rho(i)%r(:,:,:), vg(:,i) )
          DEALLOCATE( rhog, vg )
        END DO

!       ==--------------------------------------------------------------==
!       ==  WAVEFUNCTIONS                                               ==
!       ==--------------------------------------------------------------==

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
        DO j = 1, nspin
          DO i = 1, nk
            nb_g = nx
            tread = .TRUE.
            tovrw = .FALSE.
            CALL read_restart_wfc(ndr, tovrw, tread, ik_, nk_, kunit, ispin_, nspin_, &
              wfc_scal, c0, tw0, cm, twm, ngw_g, nb_g, ig_l2g, tigl, ngw )
          END DO
        END DO


      if(ionode) then
         close (unit=ndr)
      end if

      return
      end subroutine

END MODULE restart
