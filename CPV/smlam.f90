!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************************************
!
! subroutine SMLAMBDA for O-sesame/CP
!
!       042304   Y. Kanai
!
!       iterative algorithme to calculate the lambda for the
!       string method constraints.
!
!       Version : As in the paper, but the initial guess for constant C 
!                 is calculated from unparametrized string at t+dt.
!
!
!*********************************************************************


SUBROUTINE SMLAMBDA(statep,state,tan,con_ite,err_const)

  use ions_base, ONLY: na, nsp, nat
  use parameters, only: nsx
  use path_variables, ONLY: &
        sm_p => smd_p, &
        ptr  => smd_ptr, &
        maxlm => smd_maxlm


  IMPLICIT NONE

  integer :: i,j,is,ia,a,ite, isa
  integer :: sm_k,smpm

  integer :: n_const,exit_sign
  integer, intent(out) :: con_ite

  real(8), intent(out) :: err_const(sm_p)
  real(8) :: cons_c 

  type(ptr) :: statep(0:sm_p)
  type(ptr) :: state(0:sm_p)
  type(ptr) :: tan(0:sm_p)

  real(8) :: mov(3,nat,0:sm_p) 
  real(8) :: lambda(0:sm_p), dotp1, dotp2
  real(8) :: dalpha(0:sm_p),t_alpha
  real(8), ALLOCATABLE :: dotp1a(:), dotp2a(:), lam_f(:), lam_b(:)

  !_______________________________________
  !***************************************

  call start_clock( 'smlam')
  ALLOCATE( dotp1a(0:sm_p))
  ALLOCATE( dotp2a(0:sm_p))
  ALLOCATE( lam_f(0:sm_p))
  ALLOCATE( lam_b(0:sm_p))
  smpm = sm_p -1


  ! ... Number of constraints ... 

  n_const = smpm 



  ! ... Initialization ... 

  exit_sign = 0

  lambda(0:sm_p) = 0.d0
  lam_f(0:sm_p) = 0.0d0
  lam_b(0:sm_p) = 0.0d0


  ! ... Copy ... 

  DO sm_k=0,sm_p
     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              mov(i,isa,sm_k) = statep(sm_k)%d3(i,isa)
           ENDDO
        ENDDO
     ENDDO
  ENDDO



  ! ... Iteration loop ... 

  ite = 0

  ITERATION_LOOP : DO  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>!

     ite = ite+1

     IF(ite > maxlm) THEN
        ite = ite -1
        GOTO 9090 
     ENDIF

     exit_sign = 0

     call CALC(mov,n_const,exit_sign,err_const,cons_c)

     IF(exit_sign == 1) THEN
        ite = ite-1
        GOTO 9090
     ENDIF
     
     
     ! ... calculate dalpha(i) = phi(i) - phi(i-1) ...
     
     call ARC(statep,dalpha,t_alpha,1)

     ! ...  Calculate the const C ...

     call CALC(mov,n_const,exit_sign,err_const,cons_c)



     CONSTRAINT_LOOP : DO a=0,n_const  ! >>>>>>>>>>>>>>>>>>>>>>> !



        ! ... tan(l)*dphi(l) ...

        dotp1 = 0.d0

        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3

                 dotp1 = dotp1 & 
                      &    + (statep(a+1)%d3(i,isa)-statep(a)%d3(i,isa)) &
                      &    * tan(a)%d3(i,isa) 

              ENDDO
           ENDDO
        ENDDO


        ! ... tan(l+1)*dphi(l) ...

        dotp2 = 0.d0

        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3

                 dotp2 = dotp2 &
                      &    + (statep(a+1)%d3(i,isa)-statep(a)%d3(i,isa)) &
                      &    * tan(a+1)%d3(i,isa)

              ENDDO
           ENDDO
        ENDDO


        dotp1a(a) = dotp1
        dotp2a(a) = dotp2

     ENDDO CONSTRAINT_LOOP  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        
     ! do the forward loop, i.e. get a+1 from a
     CONSTRAINT_LOOP2f : DO a=0,(n_const-1)  ! >>>>>>>>>>>>>>>>>>>>>>> !
        ! ... Lagrange multiplier ...

        lam_f(a+1) = &
             & ( cons_c - (t_alpha*dalpha(a+1))**2.d0 + 2.d0 *lam_f(a) * dotp1a(a)) &
             & / (2.d0 *dotp2a(a))  
        ! ... Update ...
     ENDDO CONSTRAINT_LOOP2f  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!     write(6,*) 'lam_f after the update ',lam_f

     ! do the backward loop, i.e. get a from a+1
     CONSTRAINT_LOOP2b : DO a=n_const,1,-1  ! >>>>>>>>>>>>>>>>>>>>>>> !
        ! ... Lagrange multiplier ...

        lam_b(a) = &
             & -( cons_c - (t_alpha*dalpha(a+1))**2.d0 - 2.d0 *lam_b(a+1) * dotp2a(a)) &
             & / (2.d0 *dotp1a(a))  
        ! ... Update ...
     ENDDO CONSTRAINT_LOOP2b  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!     write(6,*) 'lam_b after the update ',lam_b

     lambda = 0.5d0*(lam_f+lam_b)
     CONSTRAINT_LOOP3 : DO a=0,n_const-1  ! >>>>>>>>>>>>>>>>>>>>>>> !
         isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3
                 mov(i,isa,a+1) = statep(a+1)%d3(i,isa) + lambda(a+1)*tan(a+1)%d3(i,isa)
              ENDDO
           ENDDO
        ENDDO
    ENDDO CONSTRAINT_LOOP3  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<!



!     write(6,*) 'lambda after the update ',lambda

  ENDDO ITERATION_LOOP  ! <<<<<<<<<<<<<<<<<<<<<<<< !



  ! ... FINAL UPDATE ...  


9090 IF(ite < maxlm) THEN
  ELSE
     call errore(' SUB. smlam ',': Maxlm exceeded. ',ite)
  ENDIF

  DO sm_k=1,smpm
     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              statep(sm_k)%d3(i,isa) = mov(i,isa,sm_k) 
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  con_ite = ite

  DEALLOCATE( dotp1a,dotp2a,lam_f,lam_b)
  call stop_clock( 'smlam')
  RETURN


END SUBROUTINE SMLAMBDA


!======================================================================================

SUBROUTINE CALC(state,n_const,exit_sign,err_const,cons) 


  use ions_base, ONLY: na,nsp, nat
  use parameters, only: nsx
  use path_variables, ONLY: &
        sm_p => smd_p, &
        tol => smd_tol


  IMPLICIT NONE

  integer :: i,is,ia,sm_k,sm_kk,smpm,ace, isa
  integer, intent(out) :: exit_sign
  integer, intent(in) :: n_const

  real(8), intent(out) :: err_const(sm_p)
  real(8), intent(out) :: cons

  real(8) :: state(3,nat,0:sm_p),temp(3,nat)
  real(8) :: dalpha(0:sm_p),t_alpha,alpha(0:sm_p)
  real(8) :: diff, total


  ! ... ARC C ... 

  ! -- seg.

  call start_clock( 'calc')
  dalpha(0) = 0.d0

  DO sm_k=1,sm_p

     dalpha(sm_k) = 0.d0

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              temp(i,isa) = state(i,isa,sm_k) - state(i,isa,sm_k-1)
           ENDDO
        ENDDO
     ENDDO

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              dalpha(sm_k) = dalpha(sm_k) + temp(i,isa)*temp(i,isa)
           ENDDO
        ENDDO
     ENDDO

     dalpha(sm_k) = DSQRT(dalpha(sm_k))
  ENDDO


  ! -- total.

  t_alpha = 0.d0
  alpha = 0.d0

  DO sm_k=0,sm_p
     DO sm_kk=0,sm_k
        alpha(sm_k) = alpha(sm_k) + dalpha(sm_kk)
     ENDDO
  ENDDO

  t_alpha = alpha(sm_p) 


  ! -- Norm.

  DO sm_k=1,sm_p
     alpha(sm_k) = alpha(sm_k)/t_alpha
     dalpha(sm_k) = dalpha(sm_k)/t_alpha
  ENDDO


  ! ** Check if the constraint is ok

  ace = 0
  exit_sign = 0

  DO sm_k=1,sm_p
     diff = DABS(dalpha(sm_k) - 1.d0/DBLE(sm_p))
     err_const(sm_k) = diff
     IF(diff <= tol) ace = ace+1
  ENDDO

  IF(ace == sm_p) THEN
     exit_sign = 1
     call stop_clock( 'calc')
     RETURN
  ENDIF


  ! ... Calc const.C

  total = 0.d0

  DO sm_k = 1,sm_p
     total = total + (dalpha(sm_k)*t_alpha)**2.d0
  ENDDO

  cons = total/DBLE(sm_p)

  call stop_clock( 'calc')
  RETURN

END SUBROUTINE CALC



