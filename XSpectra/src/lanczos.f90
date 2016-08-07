
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE lanczos (a,b,npw,psi,ncalcv,terminator)
  !----------------------------------------------------------------------------
  ! subroutine written by CG
  USE kinds,     ONLY: DP
  USE constants, ONLY: rytoev
  USE wvfct,     ONLY: npwx, nbnd
  USE becmod,    ONLY: becp
  USE uspp,      ONLY: vkb, nkb
  !USE cell_base, ONLY: omega
  USE xspectra,  ONLY: xniter, &
                       xnepoint, &
                       xcheck_conv,&
                       xnitermax,&
                       xemin,&
                       xemax,&
                       xgamma,&
                       xerror
  USE mp_global, ONLY: intra_pool_comm  
  USE mp,        ONLY: mp_sum
  USE io_global, ONLY: stdout

  IMPLICIT NONE
  ! 
  REAL(dp), DIMENSION (xnitermax), INTENT(INOUT) :: a, b
  COMPLEX(dp), DIMENSION (npwx),   INTENT(INOUT) :: psi
  INTEGER, INTENT(INOUT) :: ncalcv
  INTEGER, INTENT(IN)    :: npw
  LOGICAL, INTENT(IN)    :: terminator
  LOGICAL                :: ldummy

  !... Local variables
  LOGICAL  :: converge
  LOGICAL  :: iconv
  INTEGER  :: ibnd, j, i, m
  REAL(dp) :: norm, error, xemin_ry, xemax_ry, xgamma_ry
  COMPLEX (dp) :: ac,bc
  REAL(dp), ALLOCATABLE :: comp(:)
  COMPLEX (dp), ALLOCATABLE :: hpsi(:), u(:)

  REAL (dp) :: ddot
  COMPLEX (DP) :: zdotc
  EXTERNAL :: zdotc,ddot
  EXTERNAL :: h_psi


  ALLOCATE(hpsi(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))

  hpsi(:)=(0.d0,0.d0)
  u(:)=(0.d0,0.d0)
  a(:)=0.d0
  b(:)=0.d0

  xemax_ry=xemax/rytoev
  xemin_ry=xemin/rytoev
  xgamma_ry=xgamma/rytoev

  iconv=.false.

  ! ------------------------  1st iteration --------------------------
  !
  ! -- computes H*Psi  (|u0>=|Psi>)

  CALL h_psi( npwx, npw,1, psi, hpsi )


  ! -- computes a_(1)=<Psi|HPsi>

  a(1)=dble(zdotc(npw,psi,1,hpsi,1))

  CALL mp_sum(a(1), intra_pool_comm)

  ac=-a(1)*(1.d0,0.d0)
  !
  ! -- computes t3=H*Psi-a_1*Psi

  CALL zaxpy(npw,ac,psi,1,hpsi,1)

  !
  ! -- computes the norm of t3

  b(1) = dble(zdotc(npw,hpsi,1,hpsi,1))
  CALL mp_sum( b(1), intra_pool_comm )
  b(1) = SQRT( b(1) )
  !
  ! -- computes the vector |u1>

  CALL zdscal(npw,1.d0/b(1),hpsi,1)

  !
  ! -- saves |u1>


  u(1:npw)=hpsi(1:npw)

  hpsi(:)=(0.d0,0.d0)


  !
  ! -------------------------- Next iterations -----------------------
  !


  comp(:)=0.d0
  comp(1)=1.d0



  DO i=2,xniter


     !From here below we have:
     !   u=psi_j
     !   hpsi=empty
     !   psi=psi_j-1

     CALL h_psi( npwx, npw, 1, u, hpsi )


     !
     ! I compute hpsi=hpsi-b_{j-1}*psi_{j-1} (j is actually i-1)

     bc=-b(i-1)*(1.d0,0.d0)
     CALL zaxpy(npw,bc,psi,1,hpsi,1)

     !
     ! computes a(i)=<t2|t3>=<t2|H|t2>


     a(i)=REAL(zdotc(npw,hpsi,1,u,1),dp)
     CALL mp_sum( a(i), intra_pool_comm )
     !
     ! computes t3=t3-a(i)*t2
     !

     ac=-a(i)*(1.d0,0.d0)
     CALL zaxpy(npw,ac,u,1,hpsi,1)


     !
     ! computes b(i) the norm of t3
     !

     b(i)=REAL(zdotc(npw,hpsi,1,hpsi,1),dp)
     CALL mp_sum( b(i), intra_pool_comm )
     b(i) = SQRT( b(i) )


     !
     ! saves initial vector in t1 (t1 = t2)

     psi(1:npw)=u(1:npw)
     !
     ! computes vector t3/norm(t3) and saves it in t2

     CALL zdscal(npw,1.d0/b(i),hpsi,1)
     u(1:npw)=hpsi(1:npw)
     hpsi(1:npw)=(0.d0,0.d0)


     !
     !   I should gather all the ncalcv,error at the end and write only
     !   then.
     !
     IF(mod(i,xcheck_conv).EQ.0) THEN
        IF(converge(a,b,i,comp,error,xemin_ry,xemax_ry,&
                    xgamma_ry,xnepoint,xerror,terminator)) THEN
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   => CONVERGED at iter ',i,&
                                             ' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '|   Estimated error at iter ',i,&
                                             ': ', error
        ENDIF
     ENDIF

  ENDDO

  IF(.NOT.iconv) THEN
     ldummy=converge(a,b,i-1,comp,error,xemin_ry,xemax_ry,&
                   xgamma_ry,xnepoint,xerror,terminator) 
     WRITE(stdout,'(8x,a,i6,a)') '!   XANES not converged after', i-1,&
                                 ' iterations'
     WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   Estimated final error after ',&
          i-1,'iterations: ', error
     ncalcv=i-1
  ENDIF

  DEALLOCATE(hpsi)
  DEALLOCATE(u)
  DEALLOCATE(comp)

END SUBROUTINE lanczos


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE lanczos_uspp (a,b,npw,psi,ncalcv,terminator)
  !----------------------------------------------------------------------------
  ! subroutine written by CG
  USE kinds,     ONLY: DP
  USE wvfct,     ONLY: npwx, nbnd
  USE becmod,    ONLY: becp, calbec
  USE constants, ONLY: rytoev
  USE uspp,      ONLY: vkb, nkb
  !USE cell_base, ONLY: omega
  USE xspectra,  ONLY: xniter,&
                       xnepoint,&
                       xcheck_conv,&
                       xnitermax,&
                       xemin,&
                       xemax,&
                       xgamma,&
                       xerror
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  USE io_global, ONLY: stdout

  IMPLICIT NONE
  ! 
  REAL (dp), DIMENSION (xnitermax), INTENT (inout) :: a, b
  COMPLEX (dp), DIMENSION (npwx),   INTENT (inout) :: psi
  INTEGER, INTENT (inout) :: ncalcv
  INTEGER, INTENT (in)    :: npw
  LOGICAL, INTENT (in)    :: terminator

  !... local variables
  LOGICAL :: converge
  LOGICAL :: iconv
  LOGICAL :: recalc
  INTEGER :: ibnd,j,i,m
  REAL(dp) :: norm,error,xemin_ry,xemax_ry,xgamma_ry
  COMPLEX(dp):: ac,bc
  REAL(dp), ALLOCATABLE :: comp(:)
  COMPLEX(dp), ALLOCATABLE :: u(:), v1(:), v2(:), v3(:)
  COMPLEX(dp) :: vecteuraux1(npwx,1), vecteuraux2(npwx,1)

  REAL (dp) :: ddot
  COMPLEX(dp) :: zdotc
  EXTERNAL :: zdotc,ddot
  EXTERNAL :: h_psi

  ALLOCATE(v1(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))
  ALLOCATE(v2(npwx))
  ALLOCATE(v3(npwx))

  v1(:) = (0.d0,0.d0)
  v2(:) = (0.d0,0.d0)
  u(:)  = (0.d0,0.d0)
  a(:)  = 0.d0
  b(:)  = 0.d0

  xemax_ry = xemax/rytoev
  xemin_ry = xemin/rytoev
  xgamma_ry= xgamma/rytoev

  iconv=.false.

  ! ------------------------  1st iteration --------------------------
  !
  recalc=.true.
  CALL sm1_psi(recalc,npwx, npw, 1, psi, v1) ! -- computes v1= S^-1 psi
  CALL h_psi( npwx, npw,1, v1, u )           ! -- computes u = H v1

  ! -- computes a_(1)=<v1|u>

  a(1)=dble(zdotc(npw,v1,1,u,1))
  CALL mp_sum(a(1), intra_pool_comm)
  ac=-a(1)*(1.d0,0.d0)
  !
  ! -- computes u= S^-1 (u-a_1*psi)

  CALL zaxpy(npw,ac,psi,1,u,1)               ! -- computes u=u-a_1*psi
  recalc=.false.
  CALL sm1_psi(recalc,npwx, npw, 1, u,v1 )  ! -- computes u=S^-1 u
  !
  ! -- computes the norm

  b(1) = zdotc(npw,u,1,v1,1)                 ! -- computes b_1 =sqrt(<u|v1>)
  CALL mp_sum( b(1), intra_pool_comm )
  b(1) = SQRT( b(1) )
  !
  ! -- computes the vector |u1>
  CALL zdscal(npw,1.d0/b(1),u,1)             ! -- computes u=u/b1
  CALL zdscal(npw,1.d0/b(1),v1,1)            ! -- computes v1=v1/b1
  !
  !
  ! -------------------------- Next iterations -----------------------
  !
  comp(:)=0.d0
  comp(1)=1.d0

  DO i=2,xniter

     CALL h_psi( npwx, npw,1, v1,v2)         ! -- computes v2= H v1

     a(i)=REAL(zdotc(npw,v1,1,v2,1),dp)      ! -- computes a_i=<v1|v2>
     CALL mp_sum( a(i), intra_pool_comm )
     !
     ! I compute hpsi=hpsi-b_{j-1}*psi_{j-1} (j is actually i-1)

     bc=-b(i-1)*(1.d0,0.d0)
     CALL zaxpy(npw,bc,psi,1,v2,1)           ! -- computes v2=v2-b_{i-1}*psi

     ac=-a(i)*(1.d0,0.d0)
     CALL zaxpy(npw,ac,u,1,v2,1)             ! -- computes v2=v2-a_i*u

     v1(:)=(0.d0,0.d0)
     recalc=.false.
     CALL sm1_psi(recalc,npwx, npw, 1,v2 ,v1 ) ! computes v1= S^-1 v2

     b(i)=REAL(zdotc(npw,v2,1,v1,1),dp)      ! -- computes b_i=sqrt(<v2|v1>)
     CALL mp_sum( b(i), intra_pool_comm )
     b(i) = SQRT( b(i) )

     !
     ! saves initial vector in t1 (t1 = t2)

     psi(1:npwx)=u(1:npwx)                   ! Psi -> u
     !
     CALL zdscal(npw,1.d0/b(i),v2,1)         ! -- computes v2=v2/b_i
     CALL zdscal(npw,1.d0/b(i),v1,1)         ! -- computes v1=v1/b_i
     u(1:npwx)=v2(1:npwx)                    ! v2 -> u
     !
     !   I should gather all the ncalcv,error at the end and write only
     !   then.
     !
     IF(mod(i,xcheck_conv).EQ.0) THEN
        IF(converge(a, b, i, comp,error, xemin_ry, xemax_ry,&
                    xgamma_ry, xnepoint, xerror, terminator)) THEN
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   => CONVERGED at iter ',i,&
                                             ' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '|   Estimated error at iter ',i,&
                                             ': ', error
        ENDIF
     ENDIF

  ENDDO

  IF(.NOT.iconv) THEN
     WRITE(stdout,'(8x,a,i6,a)') '!   XANES not converged after', i-1,&
                                 ' iterations'
     WRITE(stdout,'(8x,a,i6,a,l1)') '!   Estimated final error after ',&
          i-1,'iterations: ', &
          converge(a,b,i-1,comp,error,xemin_ry,xemax_ry,&
                   xgamma_ry,xnepoint,xerror,terminator)
     ncalcv=i-1
  ENDIF

  DEALLOCATE(v1,v2,u)
  DEALLOCATE(comp)

END SUBROUTINE lanczos_uspp

!</CG>

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
LOGICAL FUNCTION converge(a,b,m,comp,estimated_error,xemin,xemax,xgamma,&
                          xnepoint,xerror,use_term)
  !----------------------------------------------------------------------------
  USE kinds, ONLY : dp
  USE xspectra, ONLY : xnitermax
  IMPLICIT NONE

  !
  ! INPUT:
  ! -----
  !
  INTEGER :: xnepoint
  REAL(dp) :: xemin,xemax,xgamma,xerror

  !
  ! INPUT / OUTPUT:
  ! ----------------
  !
  REAL(dp) :: a(xnitermax),b(xnitermax)
  INTEGER :: m        ! number of calculated vectors
  REAL(dp) :: comp(xnepoint)
  REAL(dp) :: estimated_error  
  !
  ! ---------------------- local variables ------------------------------
  !
  REAL(dp) :: deltae,tmp,area,e,err
  REAL(dp)  :: continued_fraction
  INTEGER  :: n
  LOGICAL :: use_term
  !
  deltae = (xemax-xemin) / xnepoint
  err = 0.d0
  area = 0.d0
  n = 0
  e = xemin
  DO n = 1, xnepoint
     e = e + deltae
     tmp = continued_fraction(a,b,e,xgamma,m,use_term)
     err = err + abs(comp(n)-tmp)
     area = area + abs(tmp)
     comp(n) = tmp
  ENDDO
  err = err / area
  estimated_error = err
  IF(err < xerror) THEN
     converge = .true.
  ELSE
     converge = .false.
  ENDIF
END FUNCTION converge
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION continued_fraction(a,b,e,gamma,m, term)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY: dp
  USE xspectra, ONLY: xnitermax, xcheck_conv
  USE io_global,ONLY: stdout
  IMPLICIT NONE

  ! computes the continued fraction.
  !
  ! INPUT:
  ! -----
  !
  REAL(dp) :: continued_fraction
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e
  LOGICAL :: term
  !
  ! ---------------------- local variables ------------------------------
  !
  INTEGER :: i, p,q
  COMPLEX(dp) :: res ,lastterm ! intermediate variable
  REAL(dp) :: aa, bb

  q=xcheck_conv/2
  IF (term) THEN
     aa=0.0
     bb=0.0
     DO p=1, q
        aa=aa+a(m-p)
        bb=bb+b(m-p)
     ENDDO
     aa=aa/q
     bb=bb/q
     res=lastterm(aa-e,bb*bb,gamma)
  ELSE
     res = CMPLX(a(m)-e,gamma,kind=DP)
  ENDIF
  DO i = 1, m -1
     res = CMPLX(a(m-i)-e, -gamma,kind=DP) -b(m-i)*b(m-i)/res
  ENDDO
  continued_fraction = AIMAG(1/res)
END FUNCTION continued_fraction
