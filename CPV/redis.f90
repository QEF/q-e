!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!******************************************************************************
!
! subroutine REDIS
!
!
!   For SMGD.
!
!       Library dep.    :       LAPACK & BLAS
!
!
!   NOTE :::
!
!   There is a minor int. parameter associated with polynomial interpolation
!   in SUB. INTPOL. The parameter 'mi' controls the number of replicas
!   used to interpolate a particular geom. mi is set equal to 4 by default.
!
!   This number seems to work well for the string consisting more than 7 total
!   replicas, by my experience. Yosuke
!
!*******************************************************************************


SUBROUTINE REDIS(state)

  use ions_base, ONLY: na,nsp
  use path_variables, ONLY: &
        sm_p => smd_p, &
        ptr  => smd_ptr

  IMPLICIT NONE

  integer :: a,i,is,ia,sm_k,smpm
  integer :: info,md,isa

  type(ptr) :: state(0:sm_p)


  real(8), allocatable :: t_state(:,:)
  real(8), allocatable :: state_out(:,:)

  real(8) :: dis(0:sm_p),dis_out(0:sm_p)
  real(8) :: monitor(sm_p)
  real(8) :: tmp



  ! A.    Calculate the arclength of each replicas along the string ========= 


  ! Calc. arclength between each adjucent replicas -------------------


  dis(0) = 0.d0

  DO sm_k=1,sm_p

     dis(sm_k) = 0.d0

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              dis(sm_k) = dis(sm_k) + (state(sm_k)%d3(i,isa)-state(sm_k-1)%d3(i,isa))**2.d0
           ENDDO
        ENDDO
     ENDDO

     dis(sm_k) = SQRT(dis(sm_k))
  ENDDO


  ! Calc. total arclength ----------------------------------------

  DO sm_k=1,sm_p
     dis(sm_k) = dis(sm_k-1) + dis(sm_k)
  ENDDO


  ! Normalize the arclength between replicas ------------------------

  DO sm_k=1,sm_p
     dis(sm_k) = dis(sm_k)/dis(sm_p)
  ENDDO


  ! B.    Redistribute the mesh ===========================================   

  !* replicas are always distributed equally on the mesh. Mesh actually
  !* determines how replicas lie on String.  monitor(i) is the monitor
  !* function at each half point. The mesh is equally 
  !* distributed along the string
  !* if monitor(i) = 1.0

  DO sm_k=1,sm_p

     monitor(sm_k) = 1.d0

     !IF(ci_op) THEN
     !   tmp = (ene(k-1)+ene(k))/2.d0
     !   monitor(k) = 1.d0 + tmp**5.d0 
     !ENDIF   

  ENDDO

  call MESH(monitor,dis,dis_out)


  ! C.   Interpolation to get new replicas on arc using new mesh ===============



  ! calculate the dimension MD -------- 

  md = 0

  DO is=1,nsp 
     md = md + 3*na(is)
  ENDDO

  ! Allocate arrays------------------

  allocate(t_state(md,0:sm_p))
  allocate(state_out(md,0:sm_p))


  ! Map the real coord. into this MD space --------


  DO sm_k=0,sm_p

     a = 1

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              t_state(a,sm_k) = state(sm_k)%d3(i,isa)

              a = a+1

           ENDDO
        ENDDO
     ENDDO

  ENDDO



  ! Interpolate ---------------------------

  call INTPOL(md,sm_p,dis,t_state,sm_p,dis_out,state_out)


  ! Just translate new state to original state array -----------


  DO sm_k=0,sm_p

     a = 1

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              state(sm_k)%d3(i,isa) = state_out(a,sm_k)

              a = a+1

           ENDDO
        ENDDO
     ENDDO

  ENDDO

  ! Deallocation ----------------

  if(allocated(t_state)) deallocate(t_state)
  if(allocated(state_out)) deallocate(state_out)

  RETURN

END SUBROUTINE REDIS



!*********************************************************************************
!* dis          : old mesh
!* dis_out      : new mesh such that [dis_out(i)-dis_out(i+1)]monitor(i) = constant
!* They hold the distances from the beginning of the string,
!* not distance between each
!******************************************************************************

SUBROUTINE MESH(monitor,dis,dis_out)

  use ions_base, ONLY: na,nsp
  use path_variables, ONLY: &
        sm_p => smd_p

  IMPLICIT NONE

  integer :: i,is,ia,sm_k,j,info 
  real(8) :: monitor(sm_p),dis(0:sm_p),dis_out(0:sm_p)
  real(8) :: d(sm_p-1),e(sm_p-2),b(0:sm_p)
  real(8) :: tmp(sm_p),tmp2(0:sm_p)


  ! A.    PREPARATION ==============================================================

  DO sm_k=1,sm_p
     tmp(sm_k) = 1.d0/(monitor(sm_k)*(dis(sm_k)-dis(sm_k-1)))
  ENDDO


  ! Sub-diagonal line ------------------------------------------------

  DO sm_k=1,sm_p-1
     d(sm_k) = tmp(sm_k)+tmp(sm_k+1)
  ENDDO


  ! Diagonal Line -----------------------------------------------------

  DO sm_k=1,sm_p-2
     e(sm_k) = -tmp(sm_k+1) 
  ENDDO


  ! B vector -----------------------------------------------------------

  DO sm_k=1,sm_p-1
     b(sm_k) = 0.d0 
  ENDDO

  b(sm_p-1) = tmp(sm_p)



  ! B.   SOLVE AX = B for X   ===================================================

  !* Solving AX = B matrix for X
  !* where A is tridiagonal matrix :
  !*
  !*       d e 0 0 0      0
  !*      e d e 0 0      0
  !*      0 e d e 0      0
  !*   A =   0 0 e d e    B =    0 
  !*      0 0 0 e d             b(n-1)   
  !*
  !* This is actually solving 
  !*
  !* [x(i)-x(i-1)]/[dis(i)-dis(i-1)] - [x(i+1)-x(i)]/[dis(i+1)-dis(i)] = 0
  !*
  !* with x(0) = 0, and x(n) = 1 


  ! using LAPACK subroutine

  call DPTSV(sm_p-1,1,d,e,b(1),sm_p-1,info)

  ! Now, b holds the x's (see LAPACK manual for detail)


  IF(info /= 0) THEN
     ! write(stdout,*) "SM SUB. REDIS>MESH : 001 >  building mesh failed [", info, "]"
     ! write(stdout,*) "SM SUB. REDIS>MESH : 002 >  Abnormal termination"
     call errore('REDIS - ERROR ', ': MESH failed ', info )
  ENDIF

  b(0) = 0.d0
  b(sm_p) = 1.d0

  !* If monitor(i) = 1.0, dis(i) and b(i) are the same


  ! make the index normalized 

  DO i = 0,sm_p
     tmp2(i) = DBLE(i)/DBLE(sm_p)
  ENDDO



  ! C.   interpolate to make new mesh ================================================

  call INTPOL(1,sm_p,b,dis,sm_p,tmp2,dis_out)


  !* In case of monitor(i) = 1.0, 
  !* dis() usually is no longer equally spaced line after few iterations,
  !* The above process fix that and dis_out is then equally space line.


  RETURN

END SUBROUTINE MESH



!******************************************************************************

SUBROUTINE INTPOL(md,n0,xi0,x0,n1,xi1,x1)

  ! Mi is a parameter which decides
  ! how many replicas are used to interpolate.


  use path_variables, ONLY: mi => smmi


  IMPLICIT NONE
  integer :: jlo,kk,k,i,j,n1,n0
  integer, intent(in) :: md
  real(8) :: xi0(0:n0),x0(md,0:n0),xi1(0:n1),x1(md,0:n1)
  real(8) :: xx(mi),dy

  DO i=1,n1-1

     call INTPOL_HUNT(xi0,n0+1,xi1(i),jlo)
     kk=min(max(jlo-(mi-1)/2,1),n0+2-mi)-1


     DO k=1,md
        DO j=1,mi
           xx(j)=x0(k,kk+j-1)
        ENDDO

        call INTPOL_POLINT(xi0(kk),xx(1),mi,xi1(i),x1(k,i),dy)

     ENDDO
  ENDDO

  DO k=1,md
     x1(k,0)=x0(k,0)
     x1(k,n1)=x0(k,n0)
  ENDDO

  RETURN

END SUBROUTINE INTPOL


!******************************************************************************

!
! It is haunting for the closest replicas to interplolating points.
!

SUBROUTINE INTPOL_HUNT(xx,n,x,jlo)

  IMPLICIT NONE
  integer :: inc,jhi,jlo,n,jm
  real(8) :: x,xx(n)
  integer :: skip
  logical :: ascnd

  skip = 0

  ascnd=xx(n)>=xx(1)   ! stands for ascendant


  IF(jlo <= 0 .or. jlo > n) THEN
     jlo=0
     jhi=n+1
     skip = 1   
  ENDIF


  ! Skip this part if skip=1 

  IF(skip == 0) THEN

     inc = 1

     IF(x >= xx(jlo) .eqv. ascnd) THEN

        DO 
           jhi=jlo+inc

           IF(jhi > n) THEN
              jhi = n+1
              EXIT
           ELSE IF(x >= xx(jhi) .eqv. ascnd) THEN    
              jlo = jhi
              inc = inc + inc
           ELSE
              EXIT
           ENDIF
        ENDDO
     ELSE   

        jhi=jlo

        DO 
           jlo=jhi-inc

           IF(jlo < 1) THEN
              jlo = 0
              EXIT
           ELSE IF(x < xx(jlo) .eqv. ascnd) THEN
              jhi = jlo
              inc = inc + inc
           ELSE
              EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDIF



  DO 
     IF((jhi-jlo) == 1) THEN

        IF(x == xx(n)) jlo = n-1
        IF(x == xx(1)) jlo = 1

        RETURN 
     ENDIF

     jm=(jhi+jlo)/2

     IF(x >= xx(jm) .eqv. ascnd) THEN
        jlo=jm
     ELSE   
        jhi=jm
     ENDIF
  ENDDO

END SUBROUTINE INTPOL_HUNT


!******************************************************************************

! See "numerical recipes" for details under polint().

SUBROUTINE INTPOL_POLINT(xa,ya,n,x,y,dy)


  IMPLICIT NONE
  integer, parameter :: nmax=10
  integer, intent(in) :: n
  integer :: i,ns,m
  real(8) :: dy,x,y,xa(n),ya(n)
  real(8) :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)

  ns = 1
  dif = abs(x-xa(1))

  DO i=1,n
     dift = abs(x-xa(i))

     IF(dift < dif) THEN
        ns = i 
        dif = dift
     ENDIF

     c(i) = ya(i)
     d(i) = ya(i) 
  ENDDO

  y = ya(ns)
  ns = ns-1

  DO m=1,n-1
     DO i=1,n-m
        ho = xa(i)-x
        hp = xa(i+m)-x
        w = c(i+1)-d(i)
        den = ho-hp


        IF(den == 0) &
             & call errore('SUB. INTPOL_POLINT',': 001 > interpolation failed',1) 

        den = w/den
        d(i) = hp*den
        c(i) = ho*den
     ENDDO

     IF((2*ns) < (n-m)) THEN

        dy = c(ns+1)
     ELSE 
        dy = d(ns)
        ns = ns-1
     ENDIF

     y = y+dy
  ENDDO

  RETURN

END SUBROUTINE INTPOL_POLINT

!****************************************************************************
