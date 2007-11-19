
! Copyright (C) 2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  contains two modules:
!      MODULE read_paw_module
!      MODULE paw_to_internal

!=----------------------------------------------------------------------------=!
      MODULE read_paw_module
!=----------------------------------------------------------------------------=!

!  this module handles the reading of pseudopotential data

! ...   declare modules
        USE kinds, ONLY: DP
        IMPLICIT NONE
        SAVE
        PRIVATE
        PUBLIC :: paw_io, nullify_pseudo_paw, allocate_pseudo_paw, deallocate_pseudo_paw
      CONTAINS

  !============================================================================
  !
  ! Read/write the PAW dataset
  !
  SUBROUTINE paw_io(pawset_,un,what,a_mesh,a_nwfc,a_lmax)
    use pseudo_types, only : paw_t
    use radial_grids, only : write_grid_on_file, read_grid_from_file
    USE io_global,  ONLY : stdout, ionode
    IMPLICIT NONE
    TYPE(paw_t), INTENT(INOUT) :: pawset_
    INTEGER, INTENT(IN) :: un
    CHARACTER(LEN=3), INTENT(IN) :: what
    INTEGER,OPTIONAL,INTENT(IN) :: a_mesh,a_nwfc,a_lmax
    INTEGER :: n, ns, ns1, l
    INTEGER :: mesh, nwfc, lmax
    CHARACTER(len=12) :: dummy
    SELECT CASE (what)
    CASE ("OUT")
       ! minimal data needed to alloc:
        write(un,'(a)') "init:"
       WRITE(un,'(3i8)')    pawset_%grid%mesh, pawset_%nwfc, pawset_%lmax
       ! the pseudopot:
       WRITE(un,'(A)') pawset_%symbol
        write(un,'(a)') "scalars:"
       WRITE(un,'(e20.10)') pawset_%zval
       WRITE(un,'(e20.10)') pawset_%z
       WRITE(un,'(L)')      pawset_%nlcc
       WRITE(un,'(i8)')     pawset_%nwfc
       WRITE(un,'(i8)')     pawset_%irc
       WRITE(un,'(i8)')     pawset_%lmax
       WRITE(un,'(e20.10)') pawset_%rmatch_augfun
        write(un,'(a)') "grid:"
       ! write the radial grid data
       call write_grid_on_file(un,pawset_%grid)
        write(un,'(a)') "l:"
       WRITE(un,'(i8)')     (pawset_%l(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "ikk:"
       WRITE(un,'(i8)')     (pawset_%ikk(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "oc:"
       WRITE(un,'(e20.10)') (pawset_%oc(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "els:"
       WRITE(un,'(a2)') (pawset_%els(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "jj:"
       WRITE(un,'(e20.10)') (pawset_%jj(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "rcutus:"
       WRITE(un,'(e20.10)') (pawset_%rcutus(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "enl:"
       WRITE(un,'(e20.10)') (pawset_%enl(ns), ns=1,pawset_%nwfc)
        write(un,'(a)') "aewfc:"
       WRITE(un,'(e20.10)') ((pawset_%aewfc(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        write(un,'(a)') "pswfc:"
       WRITE(un,'(e20.10)') ((pawset_%pswfc(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        write(un,'(a)') "proj:"
       WRITE(un,'(e20.10)') ((pawset_%proj(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        write(un,'(a)') "augfun:"
       WRITE(un,'(e20.10)') ((((pawset_%augfun(n,ns,ns1,l), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc),&
                                ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
        write(un,'(a)') "augmom:"
       WRITE(un,'(e20.10)') (((pawset_%augmom(ns,ns1,l), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
        write(un,'(a)') "aeccharge:"
       WRITE(un,'(e20.10)') (pawset_%aeccharge(n), n=1,pawset_%grid%mesh)
       IF (pawset_%nlcc) THEN
           write(un,'(a)') "psccharge:"
          WRITE(un,'(e20.10)') (pawset_%psccharge(n), n=1,pawset_%grid%mesh)
       ENDIF
        write(un,'(a)') "pscharge:"
       WRITE(un,'(e20.10)') (pawset_%pscharge(n), n=1,pawset_%grid%mesh)
        write(un,'(a)') "aeloc:"
       WRITE(un,'(e20.10)') (pawset_%aeloc(n), n=1,pawset_%grid%mesh)
        write(un,'(a)') "psloc:"
       WRITE(un,'(e20.10)') (pawset_%psloc(n), n=1,pawset_%grid%mesh)
        write(un,'(a)') "kdiff:"
       WRITE(un,'(e20.10)') ((pawset_%kdiff(ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
        write(un,'(a)') "dion:"
       WRITE(un,'(e20.10)') ((pawset_%dion (ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
        write(un,'(a)') "dft:"
       WRITE(un,'(A)') TRIM(pawset_%dft)
    CASE ("INP")
        read(un, '(a)') dummy
       READ(un,'(3i8)')    mesh, nwfc, lmax
       call deallocate_pseudo_paw( pawset_ )
       IF (present(a_mesh) .and. present(a_nwfc) .and. present(a_lmax) ) THEN
!             write(0,*) "Allocating from INPUT parameters"
            call allocate_pseudo_paw( pawset_, a_mesh, a_nwfc, a_lmax)
       ELSE
!             write(0,*) "Allocating from FILE parameters"
            call allocate_pseudo_paw( pawset_, mesh, nwfc, lmax)
       ENDIF
        pawset_%l(:) = 0._dp
        pawset_%ikk(:) = 0._dp
        pawset_%oc(:) = 0._dp
        pawset_%enl(:) = 0._dp
        pawset_%aewfc(:,:) = 0._dp
        pawset_%pswfc(:,:) = 0._dp
        pawset_%proj(:,:) = 0._dp
        pawset_%augfun(:,:,:,:) = 0._dp
        pawset_%augmom(:,:,:) = 0._dp
        pawset_%aeccharge(:) = 0._dp
        pawset_%psccharge(:) = 0._dp
        pawset_%pscharge(:) = 0._dp
        pawset_%aeloc(:) = 0._dp
        pawset_%psloc(:) = 0._dp
        pawset_%kdiff(:,:) = 0._dp
        pawset_%dion (:,:) = 0._dp
       READ (un,'(A)')      pawset_%symbol
        if (ionode) then
            write(stdout,"(5x,a)") &
                "Reading PAW setup for: "//TRIM(pawset_%symbol)
            write(stdout,"(7x,a,3i5)") &
                "dimensions (mesh, nwfc, lmax): ",mesh, nwfc, lmax
        endif
        read(un, '(a)') dummy
       READ (un,'(e20.10)') pawset_%zval
       READ (un,'(e20.10)') pawset_%z
       READ (un,'(L)')      pawset_%nlcc
       READ (un,'(i8)')     pawset_%nwfc
       READ (un,'(i8)')     pawset_%irc
       READ (un,'(i8)')     pawset_%lmax
       READ (un,'(e20.10)') pawset_%rmatch_augfun
       ! write the radial grid data
        read(un, '(a)') dummy
       call read_grid_from_file(un,pawset_%grid)
        read(un, '(a)') dummy
       READ (un,'(i8)')     (pawset_%l(ns), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(i8)')     (pawset_%ikk(ns), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%oc(ns), ns=1,pawset_%nwfc)
        read(un,'(a)') dummy
       READ(un,'(a2)') (pawset_%els(ns), ns=1,pawset_%nwfc)
        read(un,'(a)') dummy
       READ(un,'(e20.10)') (pawset_%jj(ns), ns=1,pawset_%nwfc)
        read(un,'(a)') dummy
       READ(un,'(e20.10)') (pawset_%rcutus(ns), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%enl(ns), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((pawset_%aewfc(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((pawset_%pswfc(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((pawset_%proj(n,ns), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((((pawset_%augfun(n,ns,ns1,l), n=1,pawset_%grid%mesh), ns=1,pawset_%nwfc),&
                                ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (((pawset_%augmom(ns,ns1,l), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%aeccharge(n), n=1,pawset_%grid%mesh)
       IF (pawset_%nlcc) THEN
           read(un, '(a)') dummy
          READ (un,'(e20.10)') (pawset_%psccharge(n), n=1,pawset_%grid%mesh)
       ENDIF
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%pscharge(n), n=1,pawset_%grid%mesh)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%aeloc(n), n=1,pawset_%grid%mesh)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') (pawset_%psloc(n), n=1,pawset_%grid%mesh)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((pawset_%kdiff(ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
        read(un, '(a)') dummy
       READ (un,'(e20.10)') ((pawset_%dion (ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
       pawset_%dft="                                                                                "
        read(un, '(a)') dummy
       READ (un,'(A)') pawset_%dft
       pawset_%dft = TRIM(pawset_%dft)
        if(ionode) write(*,"(7x,a,3i5)") "functional used: "//TRIM(pawset_%dft)
    CASE DEFAULT
       CALL errore ('paw_io','specify (INP)ut or (OUT)put',1)
    END SELECT

  END SUBROUTINE paw_io

SUBROUTINE nullify_pseudo_paw( paw )
  USE pseudo_types, ONLY: paw_t

  TYPE( paw_t ), INTENT(INOUT) :: paw
  NULLIFY( paw%l, paw%ikk )
  NULLIFY( paw%oc, paw%enl, paw%aewfc, paw%pswfc, paw%proj )
  NULLIFY( paw%augfun, paw%augmom, paw%aeccharge, paw%psccharge, paw%pscharge )
  NULLIFY( paw%aeloc, paw%psloc, paw%kdiff, paw%dion )
  RETURN
END SUBROUTINE nullify_pseudo_paw

SUBROUTINE allocate_pseudo_paw( paw, size_mesh, size_nwfc, size_lmax )
  USE pseudo_types, ONLY: paw_t

  TYPE( paw_t ), INTENT(INOUT) :: paw
  INTEGER, INTENT(IN) :: size_mesh, size_nwfc, size_lmax
  !WRITE(0,"(a,3i5)") "Allocating PAW setup: ",size_mesh, size_nwfc, size_lmax
  ALLOCATE ( paw%l(size_nwfc) )
  ALLOCATE ( paw%jj(size_nwfc) )
  ALLOCATE ( paw%ikk(size_nwfc) )
  ALLOCATE ( paw%oc(size_nwfc) )
  ALLOCATE ( paw%rcutus(size_nwfc) )
  ALLOCATE ( paw%els(size_nwfc) )
  ALLOCATE ( paw%enl(size_nwfc) )
  ALLOCATE ( paw%aewfc(size_mesh,size_nwfc) )
  ALLOCATE ( paw%pswfc(size_mesh,size_nwfc) )
  ALLOCATE ( paw%proj (size_mesh,size_nwfc) )
  ALLOCATE ( paw%augfun(size_mesh,size_nwfc,size_nwfc,0:2*size_lmax+2) )
  ALLOCATE ( paw%augmom(size_nwfc,size_nwfc,0:2*size_lmax+1) )
  ALLOCATE ( paw%aeccharge(size_mesh) )
  ALLOCATE ( paw%psccharge(size_mesh) )
  ALLOCATE ( paw%pscharge(size_mesh) )
  ALLOCATE ( paw%aeloc(size_mesh) )
  ALLOCATE ( paw%psloc(size_mesh) )
  ALLOCATE ( paw%kdiff(size_nwfc,size_nwfc) )
  ALLOCATE ( paw%dion (size_nwfc,size_nwfc) )
END SUBROUTINE allocate_pseudo_paw

SUBROUTINE deallocate_pseudo_paw( paw )
  USE pseudo_types, ONLY: paw_t

  TYPE( paw_t ), INTENT(INOUT) :: paw
  IF( ASSOCIATED( paw%l ) ) DEALLOCATE( paw%l )
  IF( ASSOCIATED( paw%jj ) ) DEALLOCATE( paw%jj )
  IF( ASSOCIATED( paw%ikk ) ) DEALLOCATE( paw%ikk )
  IF( ASSOCIATED( paw%oc ) ) DEALLOCATE( paw%oc )
  IF( ASSOCIATED( paw%els ) ) DEALLOCATE( paw%els )
  IF( ASSOCIATED( paw%rcutus ) ) DEALLOCATE( paw%rcutus )
  IF( ASSOCIATED( paw%enl ) ) DEALLOCATE( paw%enl )
  IF( ASSOCIATED( paw%aewfc ) ) DEALLOCATE( paw%aewfc )
  IF( ASSOCIATED( paw%pswfc ) ) DEALLOCATE( paw%pswfc )
  IF( ASSOCIATED( paw%proj ) ) DEALLOCATE( paw%proj )
  IF( ASSOCIATED( paw%augfun ) ) DEALLOCATE( paw%augfun )
  IF( ASSOCIATED( paw%augmom ) ) DEALLOCATE( paw%augmom )
  IF( ASSOCIATED( paw%aeccharge ) ) DEALLOCATE( paw%aeccharge )
  IF( ASSOCIATED( paw%psccharge ) ) DEALLOCATE( paw%psccharge )
  IF( ASSOCIATED( paw%pscharge ) ) DEALLOCATE( paw%pscharge )
  IF( ASSOCIATED( paw%aeloc ) ) DEALLOCATE( paw%aeloc )
  IF( ASSOCIATED( paw%psloc ) ) DEALLOCATE( paw%psloc )
  IF( ASSOCIATED( paw%kdiff ) ) DEALLOCATE( paw%kdiff )
  IF( ASSOCIATED( paw%dion ) ) DEALLOCATE( paw%dion )
  RETURN
END SUBROUTINE deallocate_pseudo_paw

!---------------------------------------------------------------------

!=----------------------------------------------------------------------------=!
      END MODULE read_paw_module
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
      MODULE paw_to_internal
!=----------------------------------------------------------------------------=!
! ...   declare modules
        USE kinds, ONLY: DP
        IMPLICIT NONE
        SAVE
        PRIVATE
        PUBLIC :: set_pseudo_paw
      CONTAINS

!---------------------------------------------------------------------
subroutine set_pseudo_paw (is, pawset)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE atom,  ONLY: rgrid, msh
  USE uspp_param, ONLY: upf, tvanp
  USE funct, ONLY: set_dft_from_name, dft_is_meta, dft_is_hybrid
  USE io_global,  ONLY : stdout, ionode
  !
  USE pseudo_types
  USE constants, ONLY: FPI
  !
  implicit none
  !
  real(DP), parameter :: rcut = 10.d0
  integer :: is, ir
  !
  !     Local variables
  !
  integer :: nb, mb, ijv
  TYPE (paw_t) :: pawset
  integer :: i,j, l, nrc, nrs
  integer :: nwfc, mesh
  real (DP) :: pow
  !
  ! Cutoffing: WARNING: arbitrary right now, for grid calculation
  pow = 1.d0
  nrs =  Count(pawset%grid%r(1:pawset%grid%mesh).le. (pawset%grid%r(pawset%irc)*1.2d0))
  nrc =  Count(pawset%grid%r(1:pawset%grid%mesh).le. (pawset%grid%r(pawset%irc)*1.8d0))
  !
  upf(is)%zp = pawset%zval
  upf(is)%psd = pawset%symbol
  upf(is)%tvanp=.true.
  upf(is)%nlcc = pawset%nlcc
  upf(is)%tpawp = .true.
  call set_dft_from_name( pawset%dft )
  !
  IF ( dft_is_meta() ) &
    CALL errore( 'set_pseudo_paw', 'META-GGA not implemented for PAW', 1 )
#if defined (EXX)
#else
  IF ( dft_is_hybrid() ) &
    CALL errore( 'set_pseudo_paw', 'HYBRID XC not implemented for PAW', 1 )
#endif

  ! to make this file barely readable:
  mesh = pawset%grid%mesh 
  nwfc = pawset%nwfc

  ! set radial grid data
  rgrid(is)%mesh = mesh
  rgrid(is)%xmin = pawset%grid%xmin
  rgrid(is)%rmax = pawset%grid%rmax
  rgrid(is)%zmesh= pawset%grid%zmesh
  rgrid(is)%dx   = pawset%grid%dx

!   upf(is)%paw%rmax  = pawset%grid%rmax
!   upf(is)%paw%zmesh = pawset%grid%zmesh
!   upf(is)%paw%xmin  = pawset%grid%xmin
!   upf(is)%paw%dx    = pawset%grid%dx
  !
  rgrid(is)%r  (1:mesh) = pawset%grid%r(1:mesh)
  rgrid(is)%rab(1:mesh) = pawset%grid%r(1:mesh)*pawset%grid%dx
  !rgrid(is)%rab(1:mesh) = pawset%grid%rab(1:mesh)
  !
  rgrid(is)%r2(1:mesh)  = pawset%grid%r2(1:mesh)
  rgrid(is)%sqr(1:mesh) = sqrt(pawset%grid%r(1:mesh))
  ! these speed up a lot a few calculations (paw XC and GCXC):
  rgrid(is)%rm1(1:mesh) = pawset%grid%rm1(1:mesh)
  rgrid(is)%rm2(1:mesh) = pawset%grid%rm2(1:mesh)
  rgrid(is)%rm3(1:mesh) = pawset%grid%rm3(1:mesh)

  !
  ! ... Copy wavefunctions used for PAW construction.
  ! ... Copy also the unoccupied ones, e.g.
  ! ... corresponding to second energy for the same channel
  ! ... (necessary to set starting occupations correctly)
  !
  upf(is)%nwfc = nwfc
  ALLOCATE ( upf(is)%lchi(nwfc), upf(is)%oc(nwfc),upf(is)%paw%oc(nwfc) )
  ALLOCATE ( upf(is)%chi( mesh, nwfc) )
  do i=1, nwfc
     upf(is)%lchi(i)=pawset%l(i)
     upf(is)%oc(i)=MAX(pawset%oc(i),0._DP)
     upf(is)%paw%oc(i)=MAX(pawset%oc(i),0._DP)
     upf(is)%chi(1:mesh, i) = pawset%pswfc(1:mesh, i)
  end do
  !
  ! Augmentation charge cutoff:
  upf(is)%paw%iraug = pawset%irc
  !
  upf(is)%nbeta= nwfc
  ! 
  allocate ( upf(is)%kbeta(nwfc) )
  do nb=1,nwfc
     upf(is)%kbeta(nb)=pawset%ikk(nb)
  end do
  ! kkbeta is the maximum distance from nucleus where AE 
  ! and PS wavefunction may not match:
  upf(is)%kkbeta=MAXVAL (upf(is)%kbeta(:))
  ! WARNING! In paw kkbeta may be smaller than the cutoff radius of augmentation function
  ! we have to keep this into account!!!
  upf(is)%paw%irmax=MAX(upf(is)%kkbeta, upf(is)%paw%iraug)
  !
  if (ionode) then
    WRITE(stdout,"(7x,a)") 'PAW functions cut-off radii:'
    WRITE(stdout,"(9x,a,i5,f12.6)") "max pfunc radius: ", upf(is)%kkbeta, pawset%grid%r(upf(is)%kkbeta)
    WRITE(stdout,"(9x,a,i5,f12.6)") "aug sphere radius:", pawset%irc, pawset%grid%r(pawset%irc)
    WRITE(stdout,"(9x,a,i5,f12.6)") "max radius:       ", upf(is)%paw%irmax, pawset%grid%r(upf(is)%paw%irmax)
  ! Currently unused:
!   WRITE(*,"(9x,a,i5,f12.6)") "inner stepping radius: ", nrs, pawset%grid%r(nrs)
!   WRITE(*,"(9x,a,i5,f12.6)") "outer stepping radius: ", nrc, pawset%grid%r(nrc)
  endif
  


  allocate (upf(is)%beta(1:mesh, 1:nwfc))
  upf(is)%beta(1:mesh, 1:nwfc) = &
                    pawset%proj(1:mesh, 1:nwfc)

  allocate(upf(is)%dion(1:nwfc, 1:nwfc))
  upf(is)%dion(1:nwfc, 1:nwfc) = pawset%dion(1:nwfc, 1:nwfc)

  allocate(upf(is)%paw%kdiff(1:nwfc, 1:nwfc))
  upf(is)%paw%kdiff(1:nwfc, 1:nwfc) = pawset%kdiff(1:nwfc, 1:nwfc)

  upf(is)%nqlc = 2*pawset%lmax+1
  upf(is)%nqf = 0                   !! no rinner, all numeric

  allocate (upf(is)%lll(nwfc) )
  upf(is)%lll(1:nwfc) = pawset%l(1:nwfc)
  allocate (upf(is)%rinner(upf(is)%nqlc))
  upf(is)%rinner(1:upf(is)%nqlc) = 0._dp  !! no rinner, all numeric

  ! integral of augmentation charges vanishes for different values of l
  allocate ( upf(is)%qqq(nwfc,nwfc))
  do i = 1, nwfc
     do j = 1, nwfc
        if (pawset%l(i)==pawset%l(j)) then
           upf(is)%qqq(i,j) = pawset%augmom(i,j,0) !!gf spherical approximation
        else
           upf(is)%qqq(i,j) = 0._dp
        end if
     end do
  end do


  ! Import total multipoles of AE-pseudo density (actually unused)
  allocate( upf(is)%paw%augmom(nwfc,nwfc, 0:2*pawset%lmax) )
  do l = 0, 2*pawset%lmax
     do i = 1, nwfc
        do j = 1, nwfc
           upf(is)%paw%augmom(i,j,l) = pawset%augmom(i,j,l)
        end do 
     end do
  end do

  ! Triangularize matrix of qfunc's
  ! FIXME!!! Probably this loop is unnecessary as the qfunc are reconstructed
  ! directly from the full augfuncs in init_us_1!!
!   allocate ( upf(is)%qfunc(1:mesh,nwfc*(nwfc+1)/2) )
!   do nb = 1, nwfc
!       do mb = nb, nwfc
!           ijv = mb * (mb-1) / 2 + nb
!           upf(is)%qfunc (1:mesh,ijv) = &
!              pawset%augfun(1:mesh,nb,mb,0)
!       enddo
!   enddo

  ! Augmentation functions for PAW depend on angular momentum!
  ! FIXME: actually they are allocated up to mesh, but they only
  ! need to be allocated up to irc!
  ! FIXME(2): augfun form a triangular matrix, so they could use a
  ! composite index ijh = nwfc*(ih-1) - ih*(ih-1)/2 + jh as for
  ! becsum, doing it would complicate paw_onecenter quite a bit...
  allocate(upf(is)%paw%aug(1:mesh,1:nwfc,1:nwfc,0:2*pawset%lmax))
  upf(is)%paw%aug(1:mesh,1:nwfc,1:nwfc,0:2*pawset%lmax) &
           = pawset%augfun(1:mesh,1:nwfc,1:nwfc,0:2*pawset%lmax)

  ! pfunc and ptfunc are couple-wise products of atomic wavefunctions
  ! FIXME: same as augfun!
  allocate( upf(is)%paw%pfunc (1:mesh, 1:nwfc,1:nwfc),&
            upf(is)%paw%ptfunc(1:mesh, 1:nwfc,1:nwfc) )
  do i=1,nwfc
     do j=1,nwfc
        upf(is)%paw%pfunc(:,i,j)  = 0._dp
        upf(is)%paw%ptfunc(:,i,j) = 0._dp
        !
        upf(is)%paw%pfunc (1:pawset%irc, i, j) = &
             pawset%aewfc(1:pawset%irc, i) * pawset%aewfc(1:pawset%irc, j)
        upf(is)%paw%ptfunc (1:pawset%irc, i, j) = &
             pawset%pswfc(1:pawset%irc, i) * pawset%pswfc(1:pawset%irc, j)
        !write(20000+100*i+10*j,'(f15.7)') upf(is)%paw%pfunc(:,i,j)
     enddo
  enddo
  !
  ! ... Add augmentation charge to ptfunc already here.
  ! ... One should not need \tilde{n}^1 alone in any case.
  !
  ! nqf is always 0 for this PAW format
  ! qfcoef(1:pawset%nqf, 1:pawset%nqlc, 1:nwfc, 1:nwfc, is ) = 0._dp
  !

  ! NO spin orbit PAW implemented right now (oct 2007)
!!$  if (lspinorb.and..not.pawset%has_so) &
!!$     call infomsg ('pawset_to_internal','At least one non s.o. pseudo')
  allocate (upf(is)%jchi(1:nwfc))
  allocate (upf(is)%jjj(1:nwfc))

!!$  lspinorb=lspinorb.and.pawset%has_so
!!$  if (pawset%has_so) then
!!$     jchi(1:nwfc, is) = pawset%jchi(1:nwfc)
!!$     jjj(1:pawset%nbeta, is) = pawset%jjj(1:pawset%nbeta)
!!$  else

  upf(is)%jchi(1:nwfc) = 0._dp
  upf(is)%jjj(1:nwfc) = 0._dp

!!$  endif

  ! Core charge: AE core charge is always present, while pseudo core can be
  ! omitted. Actually PAW uses zero core when it's not used, but this will be
  ! fixed sooner or later.
  allocate ( upf(is)%rho_atc(mesh) )
  allocate ( upf(is)%paw%ae_rho_atc(mesh) )
  !
  if ( pawset%nlcc) then
     upf(is)%rho_atc(1:mesh) = &
         pawset%psccharge(1:mesh) / FPI / pawset%grid%r2(1:mesh)
  else
     upf(is)%rho_atc(1:mesh) = 0._dp
  end if
  !
  upf(is)%paw%ae_rho_atc(1:mesh) = &
          pawset%aeccharge(1:mesh) / FPI / pawset%grid%r2(1:mesh)

  ! Atomic charge:
  allocate ( upf(is)%rho_at(mesh) )
  upf(is)%rho_at (1:mesh) = pawset%pscharge(1:mesh)


  ! Pseudo and All-electron local potentials, on radial grid
  allocate (upf(is)%vloc(1:mesh))
  allocate (upf(is)%paw%ae_vloc(mesh))
  upf(is)%vloc(1:mesh) = pawset%psloc(1:mesh)
  upf(is)%paw%ae_vloc(1:mesh) = pawset%aeloc(1:mesh)

  ! Minimum mesh point larger than cutoff radius
  msh (is) = rgrid(is)%mesh
  ir_loop : do ir = 1, rgrid(is)%mesh
    if (rgrid(is)%r (ir) > rcut) then
        msh (is) = ir
        exit ir_loop
    endif
  enddo ir_loop
  !
  ! force msh to be odd for simpson integration
  msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

end subroutine set_pseudo_paw

!=----------------------------------------------------------------------------=!
      END MODULE paw_to_internal
!=----------------------------------------------------------------------------=!

