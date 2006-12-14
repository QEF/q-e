!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ensemble_dft

  USE kinds, ONLY: DP

  IMPLICIT NONE
  SAVE

      logical      :: tens       = .false. ! whether to do ensemble calculations.
      logical      :: tgrand     = .false. ! whether to do grand canonical
                                           ! ensemble calculations.
      integer      :: ninner     = 0       ! number of inner loops per CP step.
      integer      :: ismear     = 2       ! type of smearing:
                                           !  1 => gaussian
                                           !  2 => fermi-dirac
                                           !  3 => hermite-delta_function
                                           !  4 => gaussian splines
                                           !  5 => cold smearing i
                                           !  6 => cold smearing ii
                                           ! (only 2 works).
      real(DP) :: etemp      = 0       ! smearing temperature.
      real(DP) :: ef         = 0       ! Fermi energy (relevant if tgrand=.true.).
      logical      :: tdynz      = .false. ! whether to do dynamics for the
                                           ! rotational degrees of freedom.
      logical      :: tdynf      = .false. ! whether to do dynamics for the
                                           ! unitary degrees of freedom.
      real(DP) :: zmass      = 0       ! mass for the rotational degrees of freedom
                                           ! in CP Lagrangian.
      real(DP) :: fmass      = 0       ! mass for the occupational degrees of freedom
                                           ! in CP Lagrangian.
      real(DP) :: fricz      = 0       ! unitary degrees of freedom damping.
      real(DP) :: fricf      = 0       ! occupational degrees of freedom damping.
      logical :: l_blockocc!if true consider lowest n_block states fixed
      integer :: n_blockocc(2) !number of fixed states for spin channel

!***ensemble-DFT
      real(DP), allocatable::                 bec0(:,:)
      real(DP), allocatable::                 becm(:,:)
      real(DP), allocatable::           becdrdiag(:,:,:)
      real(DP), allocatable::                  z0(:,:,:)
      real(DP), allocatable::                  id(:,:,:)
      real(DP), allocatable::               fion2(:,:)
      complex(DP), allocatable::             c0diag(:,:)
      real(DP), allocatable::               becdiag(:,:)
      real(DP), allocatable::               c0hc0(:,:,:)
      real(DP), allocatable::              c0h0c0(:,:,:)
      real(DP), allocatable::             c0hxcc0(:,:,:)
      complex(DP), allocatable::               h0c0(:,:)
      complex(DP), allocatable::              hxcc0(:,:)
      real(DP), allocatable::                  z1(:,:,:)
      real(DP), allocatable::                  zx(:,:,:)
      real(DP), allocatable::                 zxt(:,:,:)
      real(DP), allocatable::                zaux(:,:,:)
      real(DP), allocatable::                    dval(:)
      real(DP), allocatable::                      e0(:)
      real(DP), allocatable::                      e1(:)
      real(DP), allocatable::                      ex(:)
      real(DP), allocatable::                      dx(:)
      real(DP), allocatable::                      f0(:)
      real(DP), allocatable::                      f1(:)
      real(DP), allocatable::                      fx(:)
      real(DP), allocatable::                    faux(:)
      real(DP), allocatable::               fmat0(:,:,:)
      real(DP), allocatable::               fmat1(:,:,:)
      real(DP), allocatable::               fmatx(:,:,:)
      real(DP), allocatable::               dfmat(:,:,:)
      real(DP), allocatable::                     v0s(:)
      real(DP), allocatable::                 vhxcs(:,:)
      real(DP), allocatable::               epsi0(:,:,:)
      real(DP) :: atot0,atot1,atotmin,etot0,etot1,etotmin
      real(DP) :: ef1,enocc
      real(DP) :: dadx1,dedx1,dentdx1,eqa,eqb,eqc
      real(DP) :: etot2,entropy2
      real(DP) :: f2,x,xx,xmin
      complex(DP) :: c0doti,c0dotk
      integer ::  niter,nss,istart,il
      real(DP) :: gibbsfe


CONTAINS


  SUBROUTINE compute_entropy( entropy, f, nspin )
    implicit none
    real(DP), intent(out) :: entropy
    real(DP), intent(in) :: f
    integer, intent(in) :: nspin
    real(DP) :: f2
    entropy=0.0d0
    if ((f.gt.1.0d-20).and.(f.lt.(2.0/DBLE(nspin)-1.0d-20))) then
       f2=DBLE(nspin)*f/2.0d0
       entropy=-f2*log(f2)-(1.d0-f2)*log(1.d0-f2)
    end if
    entropy=-etemp*2.0d0*entropy/DBLE(nspin)
  END SUBROUTINE compute_entropy


  SUBROUTINE compute_entropy2( entropy, f, n, nspin )
    implicit none
    real(DP), intent(out) :: entropy
    real(DP), intent(in) :: f(:)
    integer, intent(in) :: n, nspin
    real(DP) :: f2
    integer :: i
    entropy=0.0d0
    do i=1,n
      if ((f(i).gt.1.0d-20).and.(f(i).lt.(2.0/DBLE(nspin)-1.0d-20))) then
        f2=DBLE(nspin)*f(i)/2.0d0
        entropy=entropy-f2*log(f2)-(1.d0-f2)*log(1.d0-f2)
      end if
    end do
    entropy=-etemp*2.0d0*entropy/DBLE(nspin)
    return
  END SUBROUTINE compute_entropy2


  SUBROUTINE compute_entropy_der( ex, fx, n, nspin )
    implicit none
    real(DP), intent(out) :: ex(:)
    real(DP), intent(in) :: fx(:)
    integer, intent(in) :: n, nspin
    real(DP) :: f2,xx
    integer :: i
    !     calculation of the entropy derivative at x
    do i=1,n
    if ((fx(i).gt.1.0d-200).and.(fx(i).lt.(2.0/DBLE(nspin)-1.0d-200))) then
      ex(i)=(log((2.0d0/DBLE(nspin)-fx(i))/fx(i)))
    else if (fx(i).le.1.0d-200) then
      xx=1.0d-200
      ex(i)=log(2.0d0/DBLE(nspin)/xx-1)
    else
      !      the calculation of ex_i is done using ex_i=-log(mf/(1-f_i)-1)
      !                                      instead of ex_i=log(mf/f_i-1)
      !      to avoid numerical errors
      xx=1.0d-200
      ex(i)=-log(2.0d0/DBLE(nspin)/xx-1)
    end if
    end do

    return
  END SUBROUTINE compute_entropy_der



  SUBROUTINE id_matrix_init( nupdwn, nspin )
    ! initialization of the matrix identity
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nupdwn(2), nspin
    INTEGER :: is, nss, i
      id(:,:,:)=0.0d0
      do  is=1,nspin
        nss=nupdwn(is)
        do i=1,nss
          id(i,i,is)=1.d0
        end do
      end do
      z0 = id  ! initialize rotation matrix to a default value
    RETURN
  END SUBROUTINE id_matrix_init


  SUBROUTINE ensemble_initval &
    ( occupations_ , n_inner_ , fermi_energy_ , rotmass_ , occmass_ , rotation_damping_ , &
      occupation_damping_ , occupation_dynamics_ , rotation_dynamics_ ,  degauss_ , smearing_ ,&
       l_blockocc_, n_blockocc_)     
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: occupations_
    CHARACTER(LEN=*), INTENT(IN) :: rotation_dynamics_
    CHARACTER(LEN=*), INTENT(IN) :: occupation_dynamics_
    CHARACTER(LEN=*), INTENT(IN) :: smearing_
    INTEGER, INTENT(IN) :: n_inner_
    REAL(DP), INTENT(IN) :: fermi_energy_ , rotmass_ , occmass_ , rotation_damping_
    REAL(DP), INTENT(IN) :: occupation_damping_ , degauss_
    LOGICAL :: l_blockocc_
    INTEGER :: n_blockocc_(2)

      SELECT CASE ( TRIM( occupations_ ) )
          !
      CASE ('bogus')
          !
      CASE ('from_input')
          !
      CASE ('fixed')
          !
      CASE ('grand-canonical','g-c','gc')
          !
          tens    =.true.
          tgrand  =.true.
          CALL errore(' ensemble_initval ','grand-canonical not yet implemented ', 1 )
          !
      CASE ('ensemble','ensemble-dft','edft')
          !
          tens    =.true.
          ninner  = n_inner_
          etemp   = degauss_
          ef      = fermi_energy_
          fricz   = rotation_damping_
          fricf   = occupation_damping_
          zmass   = rotmass_
          fmass   = occmass_
          l_blockocc = l_blockocc_
          n_blockocc(1:2) = n_blockocc_(1:2)

          SELECT CASE ( TRIM( rotation_dynamics_ ) )
            CASE ( 'line-minimization','l-m','lm' )
              tdynz = .FALSE.
              fricz = 0.0d0
              zmass = 0.0d0
            CASE DEFAULT
              CALL errore(' ensemble_initval ',' rotation_dynamics not implemented ', 1 )
          END SELECT

          SELECT CASE ( TRIM( occupation_dynamics_ ) )
            CASE ( 'line-minimization','l-m','lm' )
              tdynf = .FALSE.
              fricf = 0.0d0
              fmass = 0.0d0
            CASE DEFAULT
              CALL errore(' ensemble_initval ',' occupation_dynamics not implemented ', 1 )
          END SELECT

          SELECT CASE ( TRIM( smearing_ ) )
            CASE ( 'gaussian','g' )
              ismear = 1
            CASE ( 'fermi-dirac','f-d', 'fd' )
              ismear = 2
            CASE ( 'hermite-delta','h-d','hd' )
              ismear = 3
            CASE ( 'gaussian-splines','g-s','gs' )
              ismear = 4
            CASE ( 'cold-smearing','c-s','cs','cs1' )
              ismear = 5
            CASE ( 'marzari-vanderbilt','m-v','mv','cs2' )
              ismear = 6
            CASE ( '0')
              ismear = 0
            CASE ( '-1')
              ismear = -1
            CASE DEFAULT
              CALL errore(' ensemble_initval ',' smearing not implemented', 1 )
          END SELECT
          !
      CASE DEFAULT
          !
          CALL errore(' ensemble_initval ',' occupation method not implemented', 1 )
          !
      END SELECT
     
      IF(tens) CALL ensemble_dft_info()

    RETURN
  END SUBROUTINE ensemble_initval


  SUBROUTINE ensemble_dft_info()
    USE io_global, ONLY: stdout
      write(stdout,250) tens
      write(stdout,252) tgrand
!     write(stdout,253) tdynz
!     write(stdout,254) tdynf
250   format (4x,'  ensemble-DFT calculation     =',l5)
252   format (4x,'  grand-canonical calculation  =',l5)
253   format (4x,'  CP rotational evolution      =',l5)
254   format (4x,'  CP occupational evolution    =',l5)

      if(tens) then
         write (stdout,251) ninner,etemp,ismear,ef                         
      endif
251   format (/4x,'====================================='                          &
     &        /4x,'|      ensemble-DFT parameters      |'                          &
     &        /4x,'====================================='                          &
     &        /4x,'| ninner       =',i10,'          |'                             &
     &        /4x,'| etemp        =',f10.5,' a.u.     |'                           &
     &        /4x,'| ismear       =',i10,'          |'                             &
     &        /4x,'| fermi energy =',f10.5,' a.u.     |'                           &
!    &        /4x,'| zmass        =',f10.5,' a.u.     |'                           &
!    &        /4x,'| fmass        =',f10.5,' a.u.     |'                           &
!    &        /4x,'| fricz        =',f10.5,'          |'                           &
!    &        /4x,'| fricf        =',f10.5,'          |'                           &
     &        /4x,'=====================================')

    RETURN
  END SUBROUTINE ensemble_dft_info


  SUBROUTINE allocate_ensemble_dft( nhsa, n, ngw, nudx, nspin, nx, nnrsx, nat )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nhsa, n, ngw, nudx, nspin, nx, nnrsx, nat
      allocate( bec0(nhsa,n))
      allocate( becm(nhsa,n))
      allocate(c0diag(ngw,n))
      allocate(becdrdiag(nhsa,n,3))
      allocate(id(nudx,nudx,nspin))
      allocate(z0(nudx,nudx,nspin))
      allocate(fion2(3,nat))
      allocate(becdiag(nhsa,n))
      allocate(c0hc0(nudx,nudx,nspin))
      allocate(c0h0c0(nudx,nudx,nspin))
      allocate(c0hxcc0(nudx,nudx,nspin))
      allocate(h0c0(ngw,nx))
      allocate(hxcc0(ngw,nx))
      allocate(z1(nudx,nudx,nspin))
      allocate(zx(nudx,nudx,nspin))
      allocate(zxt(nudx,nudx,nspin))
      allocate(zaux(nudx,nudx,nspin))
      allocate(dval(nx))
      allocate(e0(nx))
      allocate(e1(nx))
      allocate(ex(nx))
      allocate(dx(nx))
      allocate(f0(nx))
      allocate(f1(nx))
      allocate(fx(nx))
      allocate(faux(nx))
      allocate(fmat0(nudx,nudx,nspin))
      allocate(fmat1(nudx,nudx,nspin))
      allocate(fmatx(nudx,nudx,nspin))
      allocate(dfmat(nudx,nudx,nspin))
      allocate(v0s(nnrsx))
      allocate(vhxcs(nnrsx,nspin))
      allocate(epsi0(nudx,nudx,nspin))
    RETURN
  END SUBROUTINE allocate_ensemble_dft



  SUBROUTINE deallocate_ensemble_dft( )
    IMPLICIT NONE
    IF( ALLOCATED( bec0 ) )  deallocate( bec0)
    IF( ALLOCATED( becm ) )  deallocate( becm )
    IF( ALLOCATED( c0diag ) )  deallocate(c0diag )
    IF( ALLOCATED( becdrdiag ) )  deallocate(becdrdiag )
    IF( ALLOCATED( id ) )  deallocate(id )
    IF( ALLOCATED( z0 ) )  deallocate(z0 )
    IF( ALLOCATED( fion2 ) )  deallocate(fion2 )
    IF( ALLOCATED( becdiag ) )  deallocate(becdiag )
    IF( ALLOCATED( c0hc0 ) )  deallocate(c0hc0 )
    IF( ALLOCATED( c0h0c0 ) )  deallocate(c0h0c0 )
    IF( ALLOCATED( c0hxcc0 ) )  deallocate(c0hxcc0 )
    IF( ALLOCATED( h0c0 ) )  deallocate(h0c0 )
    IF( ALLOCATED( hxcc0 ) )  deallocate(hxcc0 )
    IF( ALLOCATED( z1 ) )  deallocate(z1 )
    IF( ALLOCATED( zx ) )  deallocate(zx )
    IF( ALLOCATED( zxt ) )  deallocate(zxt )
    IF( ALLOCATED( zaux ) )  deallocate(zaux )
    IF( ALLOCATED( dval ) )  deallocate(dval )
    IF( ALLOCATED( e0 ) )  deallocate(e0 )
    IF( ALLOCATED( e1 ) )  deallocate(e1 )
    IF( ALLOCATED( ex ) )  deallocate(ex )
    IF( ALLOCATED( dx ) )  deallocate(dx )
    IF( ALLOCATED( f0 ) )  deallocate(f0 )
    IF( ALLOCATED( f1 ) )  deallocate(f1 )
    IF( ALLOCATED( fx ) )  deallocate(fx )
    IF( ALLOCATED( faux ) )  deallocate(faux )
    IF( ALLOCATED( fmat0 ) )  deallocate(fmat0 )
    IF( ALLOCATED( fmat1 ) )  deallocate(fmat1 )
    IF( ALLOCATED( fmatx ) )  deallocate(fmatx )
    IF( ALLOCATED( dfmat ) )  deallocate(dfmat )
    IF( ALLOCATED( v0s ) )  deallocate(v0s )
    IF( ALLOCATED( vhxcs ) )  deallocate(vhxcs )
    IF( ALLOCATED( epsi0 ) )  deallocate(epsi0 )
    RETURN
  END SUBROUTINE deallocate_ensemble_dft
  

END MODULE ensemble_dft
