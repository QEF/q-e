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
      real(DP) :: conv_thr    = 1.d-5   !energy treshold 
      real(DP) :: passop    =0.3d0    !small step for conjugate gradient
      integer :: niter_cg_restart = 20!frequency (in iterations) for restarting the cg algorith

!***
!***  Conjugate Gradient
!***
     
      COMPLEX(DP), ALLOCATABLE :: c0old(:,:)!old wfcs for extrapolation
      logical ene_ok!if .true. do not recalculate energy
      REAL(DP) :: enever!used to pass data to/from inner_loop
      INTEGER :: itercg!number of cg iterations

   !   real(DP)  ene0,ene1,dene0,enever,enesti !energy terms for linear minimization along hi
   !   real(DP)  passof,passov !step to minimum: effective, estimated
   !   integer itercg !iteration number
   !   logical ltresh!flag for convergence on energy
   !   real(DP) passo!step to minimum
   !   real(DP) etotnew,etotold!energies
   !   real(DP) spasso!sign of small step
   !   logical tcutoff!
   !   logical restartcg!if .true. restart again the CG algorithm, performing a SD step
   !   integer numok!counter on converged iterations
   !   real(DP) pcnum,pcden
   !   integer iter3
   !   real(DP) ebanda
     
      integer ninner_ef


CONTAINS


  SUBROUTINE cg_init( tcg_ , maxiter_ , conv_thr_ , passop_ ,niter_cg_restart_)
    USE kinds, ONLY: DP
    LOGICAL, INTENT(IN) :: tcg_
    INTEGER, INTENT(IN) :: maxiter_
    REAL(DP), INTENT(IN) :: conv_thr_ , passop_
    INTEGER :: niter_cg_restart_
    tcg=tcg_
    maxiter=maxiter_
    conv_thr=conv_thr_
    passop=passop_
    niter_cg_restart=niter_cg_restart_
    IF (tcg) CALL cg_info()
    RETURN
  END SUBROUTINE cg_init

  SUBROUTINE cg_info()
    USE io_global, ONLY: stdout 
    if(tcg) then
       write (stdout,400) maxiter,conv_thr,passop,niter_cg_restart                        
    endif
400 format (/4x,'========================================'                          &
   &        /4x,'|  CONJUGATE GRADIENT                  |'                          &
   &        /4x,'========================================'                          &
   &        /4x,'| iterations   =',i14,'         |'                             &
   &        /4x,'| conv_thr     =',f14.11,' a.u.    |'                           &
   &        /4x,'| passop       =',f14.5,' a.u.    |'                           &
   &        /4x,'| niter_cg_restart =',i4,'      |'                           &
   &        /4x,'========================================')
    RETURN
  END SUBROUTINE cg_info


  SUBROUTINE allocate_cg( ngw, nx, nhsavb )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ngw, nx, nhsavb
    allocate(c0old(ngw,nx))
    RETURN
  END SUBROUTINE allocate_cg

  SUBROUTINE deallocate_cg( )
    IMPLICIT NONE
    IF( ALLOCATED( c0old ) ) deallocate(c0old )
    RETURN
  END SUBROUTINE deallocate_cg

  SUBROUTINE cg_update( tfirst, nfi, c0 )
    use gvecw, only: ngw
    use electrons_base, only: n => nbsp
    IMPLICIT NONE
    COMPLEX(DP) :: c0( :, : )
    INTEGER :: nfi
    LOGICAL :: tfirst
    INTEGER :: i, ig
    if(.not. tfirst.and.(mod(nfi,10).ne.1)) then
      call DSWAP(2*ngw*n,c0,1,c0old,1)
      do i=1,n
        do ig=1,ngw
          c0(ig,i)=-c0(ig,i)+2.d0*c0old(ig,i)
        enddo
      enddo
    else
      do i=1,n
        do ig=1,ngw
          c0old(ig,i)=c0(ig,i)
        enddo
      enddo
    endif
    RETURN
  END SUBROUTINE cg_update

  SUBROUTINE print_clock_tcg()
  CALL print_clock( 'runcg_uspp')
  CALL print_clock( 'inner_loop')
  CALL print_clock( 'rotate' )
  CALL print_clock( 'calcmt' )
  CALL print_clock( 'calcm' )
  CALL print_clock( 'pc2' )
  CALL print_clock( 'pcdaga2' )
  CALL print_clock( 'set_x_minus1' )
  CALL print_clock( 'xminus1' )
  CALL print_clock( 'emass_p_tpa' )
  return
  END SUBROUTINE print_clock_tcg

END MODULE cg_module
