
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
        write(6,"(a)") " reading pseudopotential for: "//TRIM(pawset_%symbol)
        write(6,"(a,3i5)") " dimensions (mesh, nwfc, lmax): ",mesh, nwfc, lmax
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
    CASE DEFAULT
       CALL errore ('paw_io','specify (INP)ut or (OUT)put',1)
    END SELECT
#undef __DEBUG_PAW_IO
#ifdef __DEBUG_PAW_IO
    call check_paw_io(pawset_)

    !CLOSE(un) ! it should be parent's business to close the file
    CONTAINS
    SUBROUTINE check_paw_io(pawset_)
        TYPE(paw_t), INTENT(IN) :: pawset_
        INTEGER :: un = 0 ! = stderr
       WRITE(un,'(a,3i8)')    "init :", pawset_%grid%mesh, pawset_%nwfc, pawset_%lmax
       ! the pseudopot:
       WRITE(un,'(a,A)')      "symbol           :",pawset_%symbol
       WRITE(un,'(a,e20.10)') "zval             :",pawset_%zval
       WRITE(un,'(a,e20.10)') "z                :",pawset_%z
       WRITE(un,'(a,L)')      "nlcc             :",pawset_%nlcc
       WRITE(un,'(a,i8)')     "nwfc             :",pawset_%nwfc
       WRITE(un,'(a,i8)')     "irc              :",pawset_%irc
       WRITE(un,'(a,i8)')     "lmax             :",pawset_%lmax
       WRITE(un,'(a,e20.10)') "rmatch_augfun    :",pawset_%rmatch_augfun
       ! write the radial grid data
        WRITE(un,'(a,i8)')     "grid%mesh  :",pawset_%grid%mesh
        WRITE(un,'(a,e20.10)') "grid%dx    :",pawset_%grid%dx
        WRITE(un,'(a,e20.10)') "grid%xmin  :",pawset_%grid%xmin
        WRITE(un,'(a,e20.10)') "grid%zmesh :",pawset_%grid%zmesh
        WRITE(un,'(a,e20.10)') "grid%r     :",MAXVAL(pawset_%grid%r(:))
        WRITE(un,'(a,e20.10)') "grid%r     :",MAXVAL(pawset_%grid%r2(:))
        WRITE(un,'(a,e20.10)') "grid%sqr   :",MAXVAL(pawset_%grid%sqr(:))
       WRITE(un,'(a,i8)')     "l            :",MAXVAL(pawset_%l(:))
       WRITE(un,'(a,i8)')     "ikk          :",MAXVAL(pawset_%ikk(:))
       WRITE(un,'(a,e20.10)') "oc           :",MAXVAL(pawset_%oc(:))
       WRITE(un,'(a,e20.10)') "enl          :",MAXVAL(pawset_%enl(:))
       WRITE(un,'(a,e20.10)') "aewfc        :",MAXVAL(pawset_%aewfc(:,:))
       WRITE(un,'(a,e20.10)') "pswfc        :",MAXVAL(pawset_%pswfc(:,:))
       WRITE(un,'(a,e20.10)') "proj         :",MAXVAL(pawset_%proj(:,:))
       WRITE(un,'(a,e20.10)') "augfun       :",MAXVAL(pawset_%augfun(:,:,:,:))
       WRITE(un,'(a,e20.10)') "augmom       :",MAXVAL(pawset_%augmom(:,:,:))
       WRITE(un,'(a,e20.10)') "aeccharge    :",MAXVAL(pawset_%aeccharge(:))
       IF (pawset_%nlcc) WRITE(un,'(a,e20.10)') "psccharge :",MAXVAL(pawset_%psccharge(:))
       WRITE(un,'(a,e20.10)') "pscharge :",MAXVAL(pawset_%pscharge(:))
       WRITE(un,'(a,e20.10)') "aeloc    :",MAXVAL(pawset_%aeloc(:))
       WRITE(un,'(a,e20.10)') "psloc    :",MAXVAL(pawset_%psloc(:))
       WRITE(un,'(a,e20.10)') "kdiff    :",MAXVAL(pawset_%kdiff(:,:))
       WRITE(un,'(a,e20.10)') "dion     :",MAXVAL(pawset_%dion (:,:))
    END SUBROUTINE check_paw_io
#endif
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
!   WRITE(0,"(a,3i5)") "Allocating PAW setup: ",size_mesh, size_nwfc, size_lmax
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
!#define __DO_NOT_CUTOFF_PAW_FUNC
#define __SHARP_CUTOFF_PAW_FUNC
subroutine set_pseudo_paw (is, pawset)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE atom,  ONLY: rgrid, msh
  USE uspp_param, ONLY: upf
  USE funct, ONLY: set_dft_from_name, dft_is_meta, dft_is_hybrid
  !
  USE pseudo_types
  USE constants, ONLY: FPI
  !
  USE grid_paw_variables, ONLY : tpawp, pfunc, ptfunc, aevloc_at, psvloc_at, &
                                 aerho_atc, psrho_atc, kdiff, &
                                 augmom, nraug, step_f, aug 
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
  real (DP) :: pow
  real (DP), ALLOCATABLE :: aux (:)
  !
#if defined __DO_NOT_CUTOFF_PAW_FUNC
  PRINT '(A)', 'WARNING __DO_NOT_CUTOFF_PAW_FUNC'
#else
#ifdef __SHARP_CUTOFF_PAW_FUNC
  PRINT '(A)', 'WARNING __SHARP_CUTOFF_PAW_FUNC'
#endif
#endif
  !
  ! Cutoffing: WARNING: arbitrary right now, for grid calculation
  pow = 1.d0
  nrs =  Count(pawset%grid%r(1:pawset%grid%mesh).le. (pawset%grid%r(pawset%irc)*1.2d0))
  nrc =  Count(pawset%grid%r(1:pawset%grid%mesh).le. (pawset%grid%r(pawset%irc)*1.8d0))
  PRINT *, 'PAW CUTOFF parameters'
  PRINT *, pawset%irc, pawset%grid%r(pawset%irc)
  PRINT *, nrs, pawset%grid%r(nrs)
  PRINT *, nrc, pawset%grid%r(nrc)
  !
  upf(is)%zp = pawset%zval
  upf(is)%psd = pawset%symbol
  upf(is)%tvanp=.true.
  upf(is)%nlcc = pawset%nlcc
  tpawp(is)=.true.
  call set_dft_from_name( pawset%dft )
  !
  IF ( dft_is_meta() ) &
    CALL errore( 'upf_to_internal', 'META-GGA not implemented in PWscf', 1 )
#if defined (EXX)
#else
  IF ( dft_is_hybrid() ) &
    CALL errore( 'upf_to_internal', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  rgrid(is)%mesh = pawset%grid%mesh
  !
  ! ... Copy wavefunctions used for PAW construction.
  ! ... Copy also the unoccupied ones, e.g.
  ! ... corresponding to second energy for the same channel
  ! ... (necessary to set starting occupations correctly)
  !
  upf(is)%nwfc = pawset%nwfc
  ALLOCATE ( upf(is)%lchi(pawset%nwfc), upf(is)%oc(pawset%nwfc) )
  ALLOCATE ( upf(is)%chi( pawset%grid%mesh, pawset%nwfc) )
  do i=1, pawset%nwfc
     upf(is)%lchi(i)=pawset%l(i)
     upf(is)%oc(i)=MAX(pawset%oc(i),0._DP)
     upf(is)%chi(1:pawset%grid%mesh, i) = pawset%pswfc(1:pawset%grid%mesh, i)
  end do
  !
  upf(is)%nbeta= pawset%nwfc
  allocate ( upf(is)%kbeta(pawset%nwfc) )
  do nb=1,pawset%nwfc
     upf(is)%kbeta(nb)=pawset%ikk(nb)
  end do
  ! kkbeta is the maximum distance from nucleus where AE 
  ! and PS wavefunction may not match:
  upf(is)%kkbeta=MAXVAL (upf(is)%kbeta(:))

  allocate (upf(is)%beta(1:pawset%grid%mesh, 1:pawset%nwfc))
  upf(is)%beta(1:pawset%grid%mesh, 1:pawset%nwfc) = &
                    pawset%proj(1:pawset%grid%mesh, 1:pawset%nwfc)

  allocate(upf(is)%dion(1:pawset%nwfc, 1:pawset%nwfc))
  upf(is)%dion(1:pawset%nwfc, 1:pawset%nwfc) = pawset%dion(1:pawset%nwfc, 1:pawset%nwfc)
  kdiff(1:pawset%nwfc, 1:pawset%nwfc, is) = pawset%kdiff(1:pawset%nwfc, 1:pawset%nwfc)

  upf(is)%nqlc = 2*pawset%lmax+1
  upf(is)%nqf = 0                   !! no rinner, all numeric

  allocate (upf(is)%lll(pawset%nwfc) )
  upf(is)%lll(1:pawset%nwfc) = pawset%l(1:pawset%nwfc)
  allocate (upf(is)%rinner(upf(is)%nqlc))
  upf(is)%rinner(1:upf(is)%nqlc) = 0._dp  !! no rinner, all numeric

  ! integral of augmentation charges vanishes for different values of l
  allocate ( upf(is)%qqq(pawset%nwfc,pawset%nwfc))
  do i = 1, pawset%nwfc
     do j = 1, pawset%nwfc
        if (pawset%l(i)==pawset%l(j)) then
           upf(is)%qqq(i,j) = pawset%augmom(i,j,0) !!gf spherical approximation
        else
           upf(is)%qqq(i,j) = 0._dp
        end if
     end do
  end do

  ! import augmentation charge:
  nraug(is) = pawset%irc
  rgrid(is)%dx = pawset%grid%dx
  rgrid(is)%r2 (1:pawset%grid%mesh) = pawset%grid%r2  (1:pawset%grid%mesh)
  do l = 0, 2*pawset%lmax
     do i = 1, pawset%nwfc
        do j = 1, pawset%nwfc
           augmom(i,j,l,is) = pawset%augmom(i,j,l)
        end do 
     end do
  end do

  ! triangularize matrix of qfunc's
  allocate ( upf(is)%qfunc(1:pawset%grid%mesh,pawset%nwfc*(pawset%nwfc+1)/2) )
  do nb = 1, pawset%nwfc
      do mb = nb, pawset%nwfc
          ijv = mb * (mb-1) / 2 + nb
          upf(is)%qfunc (1:pawset%grid%mesh,ijv) = &
             pawset%augfun(1:pawset%grid%mesh,nb,mb,0)
      enddo
  enddo

  ! new sparse allocation for augmentation functions
  if (allocated(aug(is)%fun)) deallocate(aug(is)%fun)
  allocate(aug(is)%fun(1:pawset%grid%mesh,1:pawset%nwfc,1:pawset%nwfc,0:2*pawset%lmax))
  !
    aug(is)%fun(1:pawset%grid%mesh,1:pawset%nwfc,1:pawset%nwfc,0:2*pawset%lmax) &
       = &
  pawset%augfun(1:pawset%grid%mesh,1:pawset%nwfc,1:pawset%nwfc,0:2*pawset%lmax)


  !
#if ! defined __DO_NOT_CUTOFF_PAW_FUNC
  ALLOCATE (aux(1:pawset%grid%mesh))
#endif
  do i=1,pawset%nwfc
     do j=1,pawset%nwfc
#if defined __DO_NOT_CUTOFF_PAW_FUNC
        pfunc (1:pawset%grid%mesh, i, j, is) = &
             pawset%aewfc(1:pawset%grid%mesh, i) * pawset%aewfc(1:pawset%grid%mesh, j)
        ptfunc (1:pawset%grid%mesh, i, j, is) = &
             pawset%pswfc(1:pawset%grid%mesh, i) * pawset%pswfc(1:pawset%grid%mesh, j)
#else
#if defined __SHARP_CUTOFF_PAW_FUNC
        pfunc(:,i,j,is) = 0._dp
        ptfunc(:,i,j,is) = 0._dp
        pfunc (1:pawset%irc, i, j, is) = &
             pawset%aewfc(1:pawset%irc, i) * pawset%aewfc(1:pawset%irc, j)
        ptfunc (1:pawset%irc, i, j, is) = &
             pawset%pswfc(1:pawset%irc, i) * pawset%pswfc(1:pawset%irc, j)
#else
        aux(1:pawset%grid%mesh) = pawset%aewfc(1:pawset%grid%mesh, i) * &
             pawset%aewfc(1:pawset%grid%mesh, j)
        CALL step_f( pfunc(1:pawset%grid%mesh,i,j,is), aux(1:pawset%grid%mesh), &
             pawset%grid%r(1:pawset%grid%mesh), nrs, nrc, pow, pawset%grid%mesh)
        aux(1:pawset%grid%mesh) = pawset%pswfc(1:pawset%grid%mesh, i) * &
             pawset%pswfc(1:pawset%grid%mesh, j)
        CALL step_f( ptfunc(1:pawset%grid%mesh,i,j,is), aux(1:pawset%grid%mesh), &
             pawset%grid%r(1:pawset%grid%mesh), nrs, nrc, pow, pawset%grid%mesh)
#endif
#endif
     end do
  end do
  !
  ! ... Add augmentation charge to ptfunc already here.
  ! ... One should not need \tilde{n}^1 alone in any case.
  !
  ! nqf is always 0 for this PAW format
  ! qfcoef(1:pawset%nqf, 1:pawset%nqlc, 1:pawset%nwfc, 1:pawset%nwfc, is ) = 0._dp
  !
  rgrid(is)%r  (1:pawset%grid%mesh) = pawset%grid%r  (1:pawset%grid%mesh)
  rgrid(is)%rab(1:pawset%grid%mesh) = pawset%grid%r  (1:pawset%grid%mesh)*pawset%grid%dx

!
! set radial grid data
!
  rgrid(is)%mesh = pawset%grid%mesh
  rgrid(is)%xmin = pawset%grid%xmin
  rgrid(is)%rmax = pawset%grid%rmax
  rgrid(is)%zmesh= pawset%grid%zmesh
  rgrid(is)%dx   = pawset%grid%dx
  rgrid(is)%r(1:pawset%grid%mesh)   = pawset%grid%r(1:pawset%grid%mesh)
  rgrid(is)%r2(1:pawset%grid%mesh)  = pawset%grid%r2(1:pawset%grid%mesh)
  rgrid(is)%rab(1:pawset%grid%mesh) = pawset%grid%rab(1:pawset%grid%mesh)
  rgrid(is)%sqr(1:pawset%grid%mesh) = sqrt(pawset%grid%r(1:pawset%grid%mesh))
! these speed up a lot a few calculations:
  rgrid(is)%rm1(1:pawset%grid%mesh)   = pawset%grid%rm1(1:pawset%grid%mesh)
  rgrid(is)%rm2(1:pawset%grid%mesh)   = pawset%grid%rm2(1:pawset%grid%mesh)
  rgrid(is)%rm3(1:pawset%grid%mesh)   = pawset%grid%rm3(1:pawset%grid%mesh)

  ! NO spin orbit PAW implemented right now (oct 2007)
!!$  if (lspinorb.and..not.pawset%has_so) &
!!$     call infomsg ('pawset_to_internal','At least one non s.o. pseudo')
  allocate (upf(is)%jchi(1:pawset%nwfc))
  allocate (upf(is)%jjj(1:pawset%nwfc))

!!$  lspinorb=lspinorb.and.pawset%has_so
!!$  if (pawset%has_so) then
!!$     jchi(1:pawset%nwfc, is) = pawset%jchi(1:pawset%nwfc)
!!$     jjj(1:pawset%nbeta, is) = pawset%jjj(1:pawset%nbeta)
!!$  else

  upf(is)%jchi(1:pawset%nwfc) = 0._dp
  upf(is)%jjj(1:pawset%nwfc) = 0._dp

!!$  endif
  !
  if ( pawset%nlcc) then
     allocate ( upf(is)%rho_atc(pawset%grid%mesh) )
     upf(is)%rho_atc(1:pawset%grid%mesh) = &
	pawset%psccharge(1:pawset%grid%mesh)   &
        / FPI / pawset%grid%r2(1:pawset%grid%mesh)
  end if

  aerho_atc(1:pawset%grid%mesh, is) = pawset%aeccharge(1:pawset%grid%mesh) &
       &                         / FPI / pawset%grid%r2(1:pawset%grid%mesh)
  if ( pawset%nlcc) then
     psrho_atc(1:pawset%grid%mesh, is) = pawset%psccharge(1:pawset%grid%mesh) &
          &                         / FPI / pawset%grid%r2(1:pawset%grid%mesh)
  else
     psrho_atc(:,is) = 0._dp
  end if
  !
  allocate ( upf(is)%rho_at(pawset%grid%mesh) )
  upf(is)%rho_at (1:pawset%grid%mesh) = pawset%pscharge(1:pawset%grid%mesh)

  allocate (upf(is)%vloc(1:pawset%grid%mesh))
  upf(is)%vloc(1:pawset%grid%mesh) = pawset%psloc(1:pawset%grid%mesh)
#if defined __DO_NOT_CUTOFF_PAW_FUNC
  aevloc_at(1:pawset%grid%mesh,is) = pawset%aeloc(1:pawset%grid%mesh)
  psvloc_at(1:pawset%grid%mesh,is) = pawset%psloc(1:pawset%grid%mesh)
#else
#if defined __SHARP_CUTOFF_PAW_FUNC
  aevloc_at(1:pawset%grid%mesh,is) = pawset%aeloc(1:pawset%grid%mesh)
  psvloc_at(1:pawset%grid%mesh,is) = pawset%psloc(1:pawset%grid%mesh)
 aux(1:pawset%grid%mesh) = pawset%aeloc(1:pawset%grid%mesh)
#else
 CALL step_f( aevloc_at(1:pawset%grid%mesh,is), aux(1:pawset%grid%mesh), &
      pawset%grid%r(1:pawset%grid%mesh), nrs, nrc, pow, pawset%grid%mesh)
 aux(1:pawset%grid%mesh) = pawset%psloc(1:pawset%grid%mesh)
 CALL step_f( psvloc_at(1:pawset%grid%mesh,is), aux(1:pawset%grid%mesh), &
      pawset%grid%r(1:pawset%grid%mesh), nrs, nrc, pow, pawset%grid%mesh)
 DEALLOCATE (aux)
#endif
#endif

  do ir = 1, rgrid(is)%mesh
    if (rgrid(is)%r (ir) .gt.rcut) then
        msh (is) = ir
        goto 5
    endif
  enddo
  msh (is) = rgrid(is)%mesh
  !
  ! force msh to be odd for simpson integration
  !
5 msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

end subroutine set_pseudo_paw

!=----------------------------------------------------------------------------=!
      END MODULE paw_to_internal
!=----------------------------------------------------------------------------=!

