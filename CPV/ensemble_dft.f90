MODULE ensemble_dft

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
      real(kind=8) :: etemp      = 0       ! smearing temperature.
      real(kind=8) :: ef         = 0       ! Fermi energy (relevant if tgrand=.true.).
      logical      :: tdynz      = .false. ! whether to do dynamics for the
                                           ! rotational degrees of freedom.
      logical      :: tdynf      = .false. ! whether to do dynamics for the
                                           ! unitary degrees of freedom.
      real(kind=8) :: zmass      = 0       ! mass for the rotational degrees of freedom
                                           ! in CP Lagrangian.
      real(kind=8) :: fmass      = 0       ! mass for the occupational degrees of freedom
                                           ! in CP Lagrangian.
      real(kind=8) :: fricz      = 0       ! unitary degrees of freedom damping.
      real(kind=8) :: fricf      = 0       ! occupational degrees of freedom damping.

!***ensemble-DFT
      real(kind=8), allocatable::                 bec0(:,:)
      real(kind=8), allocatable::                 becm(:,:)
      real(kind=8), allocatable::           becdrdiag(:,:,:)
      real(kind=8), allocatable::                  z0(:,:,:)
      real(kind=8), allocatable::                  id(:,:,:)
      real(kind=8), allocatable::               fion2(:,:)
      complex(kind=8), allocatable::             c0diag(:,:)
      real(kind=8), allocatable::               becdiag(:,:)
      real(kind=8), allocatable::               c0hc0(:,:,:)
      real(kind=8), allocatable::              c0h0c0(:,:,:)
      real(kind=8), allocatable::             c0hxcc0(:,:,:)
      complex(kind=8), allocatable::               h0c0(:,:)
      complex(kind=8), allocatable::              hxcc0(:,:)
      real(kind=8), allocatable::                  z1(:,:,:)
      real(kind=8), allocatable::                  zx(:,:,:)
      real(kind=8), allocatable::                 zxt(:,:,:)
      real(kind=8), allocatable::                zaux(:,:,:)
      real(kind=8), allocatable::                    dval(:)
      real(kind=8), allocatable::                      e0(:)
      real(kind=8), allocatable::                      e1(:)
      real(kind=8), allocatable::                      ex(:)
      real(kind=8), allocatable::                      dx(:)
      real(kind=8), allocatable::                      f0(:)
      real(kind=8), allocatable::                      f1(:)
      real(kind=8), allocatable::                      fx(:)
      real(kind=8), allocatable::                    faux(:)
      real(kind=8), allocatable::               fmat0(:,:,:)
      real(kind=8), allocatable::               fmat1(:,:,:)
      real(kind=8), allocatable::               fmatx(:,:,:)
      real(kind=8), allocatable::               dfmat(:,:,:)
      real(kind=8), allocatable::                     v0s(:)
      real(kind=8), allocatable::                 vhxcs(:,:)
      real(kind=8), allocatable::               epsi0(:,:,:)
      real(kind=8) :: atot0,atot1,atotmin,etot0,etot1,etotmin
      real(kind=8) :: ef1,enocc
      real(kind=8) :: dadx1,dedx1,dentdx1,eqa,eqb,eqc
      real(kind=8) :: etot2,entropy2
      real(kind=8) :: f2,x,xx,xmin
      complex(kind=8) :: c0doti,c0dotk
      integer ::  niter,nss,istart,il
      real(kind=8) :: gibbsfe


CONTAINS


  SUBROUTINE compute_entropy( entropy, f, nspin )
    implicit none
    real(kind=8), intent(out) :: entropy
    real(kind=8), intent(in) :: f
    integer, intent(in) :: nspin
    real(kind=8) :: f2
    entropy=0.0
    if ((f.gt.1.0d-20).and.(f.lt.(2.0/float(nspin)-1.0d-20))) then
       f2=float(nspin)*f/2.0
       entropy=-f2*log(f2)-(1.-f2)*log(1.-f2)
    end if
    entropy=-etemp*2.0*entropy/float(nspin)
  END SUBROUTINE


  SUBROUTINE compute_entropy2( entropy, f, n, nspin )
    implicit none
    real(kind=8), intent(out) :: entropy
    real(kind=8), intent(in) :: f(:)
    integer, intent(in) :: n, nspin
    real(kind=8) :: f2
    integer :: i
    entropy=0.0
    do i=1,n
      if ((f(i).gt.1.0d-20).and.(f(i).lt.(2.0/float(nspin)-1.0d-20))) then
        f2=float(nspin)*f(i)/2.0
        entropy=entropy-f2*log(f2)-(1.-f2)*log(1.-f2)
      end if
    end do
    entropy=-etemp*2.0*entropy/float(nspin)
    return
  END SUBROUTINE


  SUBROUTINE compute_entropy_der( ex, fx, n, nspin )
    implicit none
    real(kind=8), intent(out) :: ex(:)
    real(kind=8), intent(in) :: fx(:)
    integer, intent(in) :: n, nspin
    real(kind=8) :: f2,xx
    integer :: i
    !     calculation of the entropy derivative at x
    do i=1,n
    if ((fx(i).gt.1.0d-200).and.(fx(i).lt.(2.0/float(nspin)-1.0d-200))) then
      ex(i)=(log((2.0/float(nspin)-fx(i))/fx(i)))
    else if (fx(i).le.1.0d-200) then
      xx=1.0d-200
      ex(i)=log(2.0/float(nspin)/xx-1)
    else
      !      the calculation of ex_i is done using ex_i=-log(mf/(1-f_i)-1)
      !                                      instead of ex_i=log(mf/f_i-1)
      !      to avoid numerical errors
      xx=1.0d-200
      ex(i)=-log(2.0/float(nspin)/xx-1)
    end if
    end do

    return
  END SUBROUTINE



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
  END SUBROUTINE



  SUBROUTINE enemble_dft_info()
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
  END SUBROUTINE


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
  END SUBROUTINE



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
  END SUBROUTINE
  

END MODULE
