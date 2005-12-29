!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE cg_module

 USE kinds, ONLY: DP

  IMPLICIT NONE
  SAVE

      logical      :: tcg        = .false.   ! if true do conjugate gradient minimization for electrons
      integer      :: maxiter    = 100      ! maximum number of iterations
      real(DP) :: etresh    = 1.d-5   !energy treshold 
      real(DP) :: passop    =0.3d0    !small step for conjugate gradient

!***
!***  Conjugate Gradient
!***
      real(DP)  esse,essenew !factors in c.g.
      COMPLEX(DP), ALLOCATABLE :: gi(:,:)!conjugates 
      COMPLEX(DP), ALLOCATABLE :: hi(:,:)!gradients 
      COMPLEX(DP), ALLOCATABLE :: c0old(:,:)!old wfcs for extrapolation
      COMPLEX(DP), ALLOCATABLE :: hpsi(:,:) !terms H|Psi_i>
      real(DP), allocatable::               s_minus1(:,:)!factors for inverting US S matrix
      real(DP), allocatable::               k_minus1(:,:)!factors for inverting US preconditioning matrix 
      real(DP)  ene0,ene1,dene0,enever,enesti !energy terms for linear minimization along hi
      real(DP)  passof,passov !step to minimum: effective, estimated
      integer itercg !iteration number
      logical ltresh!flag for convergence on energy
      real(DP) passo!step to minimum
      real(DP) etotnew,etotold!energies
      real(DP) spasso!sign of small step
      logical tcutoff!
      real(DP), ALLOCATABLE :: emme(:,:)!matrix used for cal_emme style of projection
      logical restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer numok!counter on converged iterations
      real(DP) pcnum,pcden
      integer iter3
      real(DP) ebanda
      logical ene_ok!if .true. do not recalculate energy
      integer ninner_ef


CONTAINS


  SUBROUTINE cg_init( tcg_ , maxiter_ , etresh_ , passop_ )
    USE kinds, ONLY: DP
    LOGICAL, INTENT(IN) :: tcg_
    INTEGER, INTENT(IN) :: maxiter_
    REAL(DP), INTENT(IN) :: etresh_ , passop_
    tcg=tcg_
    maxiter=maxiter_
    etresh=etresh_
    passop=passop_
    IF (tcg) CALL cg_info()
    RETURN
  END SUBROUTINE cg_init

  SUBROUTINE cg_info()
    USE io_global, ONLY: stdout 
    if(tcg) then
       write (stdout,400) maxiter,etresh,passop                         
    endif
400 format (/4x,'====================================='                          &
   &        /4x,'|  CONJUGATE GRADIENT               |'                          &
   &        /4x,'====================================='                          &
   &        /4x,'| iterations   =',i10,'            |'                             &
   &        /4x,'| etresh       =',f10.5,' a.u.     |'                           &
   &        /4x,'| passop       =',f10.5,' a.u.     |'                           &
   &        /4x,'=====================================')
    RETURN
  END SUBROUTINE cg_info


  SUBROUTINE allocate_cg( ngw, nx, nhsavb )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx, nhsavb
    allocate(hi(ngw,nx))
    allocate(gi(ngw,nx))
    allocate(c0old(ngw,nx))
    allocate( emme(nx,nx))
    allocate( hpsi(ngw,nx))
    allocate( s_minus1(nhsavb,nhsavb))
    allocate( k_minus1(nhsavb,nhsavb))
    RETURN
  END SUBROUTINE allocate_cg

  SUBROUTINE deallocate_cg( )
    IMPLICIT NONE
    IF( ALLOCATED( hi ) ) deallocate(hi )
    IF( ALLOCATED( gi ) ) deallocate(gi )
    IF( ALLOCATED( c0old ) ) deallocate(c0old )
    IF( ALLOCATED( emme ) ) deallocate( emme )
    IF( ALLOCATED( hpsi ) ) deallocate( hpsi )
    IF( ALLOCATED( s_minus1) ) deallocate( s_minus1)
    IF( ALLOCATED( k_minus1) ) deallocate( k_minus1)
    RETURN
  END SUBROUTINE deallocate_cg

  SUBROUTINE cg_update( tfirst, nfi, c0 )
    use gvecw, only: ngw
    use electrons_base, only: n => nbsp
    IMPLICIT NONE
    COMPLEX(DP) :: c0( :, :, :, : )
    INTEGER :: nfi
    LOGICAL :: tfirst
    INTEGER :: i, ig
    if(.not. tfirst.and.(mod(nfi,10).ne.1)) then
      call DSWAP(2*ngw*n,c0,1,c0old,1)
      do i=1,n
        do ig=1,ngw
          c0(ig,i,1,1)=-c0(ig,i,1,1)+2.d0*c0old(ig,i)
        enddo
      enddo
    else
      do i=1,n
        do ig=1,ngw
          c0old(ig,i)=c0(ig,i,1,1)
        enddo
      enddo
    endif
    RETURN
  END SUBROUTINE cg_update


END MODULE cg_module
