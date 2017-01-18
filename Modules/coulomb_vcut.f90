!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Giovanni Bussi
! Adapted to QE by Andrea Ferretti & Layla Martin Samos
!
!----------------------------------
  MODULE coulomb_vcut_module
  !----------------------------------
  !
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose parameters
  !
  INTEGER,  PARAMETER :: DP=KIND(1.0d0)
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: TPI    = 2.0_DP * pi
  REAL(DP), PARAMETER :: FPI    = 4.0_DP * pi
  REAL(DP), PARAMETER :: e2     = 2.0_DP
  REAL(DP), PARAMETER :: eps6   = 1.0E-6_DP
 
  !
  ! definitions
  ! 
  TYPE vcut_type
      REAL(DP)          :: a(3,3)
      REAL(DP)          :: b(3,3)
      REAL(DP)          :: a_omega
      REAL(DP)          :: b_omega
      REAL(DP), POINTER :: corrected(:,:,:)
      REAL(DP)          :: cutoff
      LOGICAL           :: orthorombic
  END TYPE vcut_type

  !
  PUBLIC :: vcut_type
  PUBLIC :: vcut_init
  PUBLIC :: vcut_get
  PUBLIC :: vcut_spheric_get
  PUBLIC :: vcut_destroy
  PUBLIC :: vcut_info

CONTAINS

!------------------------------------------
  SUBROUTINE vcut_init(vcut,a,cutoff)
  !------------------------------------------
  !
  TYPE(vcut_type),   INTENT(OUT) :: vcut
  REAL(DP),          INTENT(IN)  :: a(3,3)
  REAL(DP),          INTENT(IN)  :: cutoff

  INTEGER      :: n1,n2,n3
  INTEGER      :: i1,i2,i3
  INTEGER      :: ierr
  REAL(DP)     :: q(3)
  CHARACTER(9) :: subname='vcut_init'     
  REAL(DP)     :: mod2a(3)

  vcut%cutoff=cutoff

  vcut%a=a
  vcut%b= TPI * transpose(num_inverse(vcut%a))
  vcut%b_omega=num_determinant(vcut%b)
  vcut%a_omega=num_determinant(vcut%a)

  ! automatically finds whether the cell is orthorombic or not
  vcut%orthorombic=.false.
  !
  mod2a=sqrt(sum(vcut%a**2,1))
  if(abs(sum(vcut%a(:,1)*vcut%a(:,2)))/(mod2a(1)*mod2a(2))<eps6 .and. &
     abs(sum(vcut%a(:,2)*vcut%a(:,3)))/(mod2a(2)*mod2a(3))<eps6 .and. &
     abs(sum(vcut%a(:,3)*vcut%a(:,1)))/(mod2a(3)*mod2a(1))<eps6) vcut%orthorombic=.true.
  !
  if (.not.vcut%orthorombic) call errore(subname,"non-orthorombic case untested",1)

  n1=ceiling(vcut%cutoff*sqrt(sum(vcut%a(1,:)**2))/(2.0*pi))
  n2=ceiling(vcut%cutoff*sqrt(sum(vcut%a(2,:)**2))/(2.0*pi))
  n3=ceiling(vcut%cutoff*sqrt(sum(vcut%a(3,:)**2))/(2.0*pi))

  ALLOCATE(vcut%corrected(-n1:n1,-n2:n2,-n3:n3), STAT=ierr)
  IF ( ierr/=0 ) CALL errore(subname,'allocating cvut%corrected',ABS(ierr))
  !
  vcut%corrected=0.0

  !
  ! define the Fourier component of the modified Coulomb potential
  !
  DO i1=-n1,n1
    DO i2=-n2,n2
      DO i3=-n3,n3
        !
        q = MATMUL(vcut%b,(/i1,i2,i3/)) 
        !
        IF( SUM(q**2) > vcut%cutoff**2 ) CYCLE
        !
        vcut%corrected(i1,i2,i3) = &
             vcut_formula(q,vcut%a,vcut%b,vcut%a_omega,vcut%orthorombic)
        !
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE vcut_init

!------------------------------------------
  SUBROUTINE vcut_info(iun, vcut)
  !------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN) :: iun
  TYPE(vcut_type), INTENT(IN) :: vcut
  !
  INTEGER :: i, n(3)
  !
  IF ( ASSOCIATED( vcut%corrected ) ) THEN
     !
     DO i = 1, 3
        n(i) = ( SIZE( vcut%corrected, i) -1 ) / 2
     ENDDO
     !
     WRITE(iun, "(  2x,'Cutoff: ',f6.2,4x,'  n grid: ',3i4,/)") vcut%cutoff, n(:)
     !
  ENDIF
  !
END SUBROUTINE vcut_info

!------------------------------------------
  SUBROUTINE vcut_destroy(vcut)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(INOUT) :: vcut
  INTEGER :: ierr
  !
  DEALLOCATE(vcut%corrected, STAT=ierr)
  IF ( ierr/=0 ) CALL errore('vcut_destroy','deallocating vcut',ABS(ierr))
  !
END SUBROUTINE vcut_destroy

!------------------------------------------
  FUNCTION vcut_get(vcut,q) RESULT(res)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(IN) :: vcut
  REAL(DP),        INTENT(IN) :: q(3)
  REAL(DP)                    :: res
  !
  REAL(DP)     :: i_real(3)
  INTEGER      :: i(3)
  CHARACTER(8) :: subname='vcut_get'
  !
  i_real=(MATMUL(TRANSPOSE(vcut%a),q))/ TPI
  i=NINT(i_real)
  !
  ! internal check
  IF( SUM( (i-i_real)**2 ) > eps6 ) &
     CALL errore(subname,'q vector out of the grid',10)
  !
  IF( SUM(q**2) > vcut%cutoff**2 ) THEN
     !
     ! usual form of Coulomb potential
     res = FPI * e2 / SUM(q**2) 
     !
  ELSE
     !
     IF( i(1)>ubound(vcut%corrected,1) .OR. i(1)<lbound(vcut%corrected,1) .OR. &
         i(2)>ubound(vcut%corrected,2) .OR. i(2)<lbound(vcut%corrected,2) .OR. &
         i(3)>ubound(vcut%corrected,3) .OR. i(3)<lbound(vcut%corrected,3)) THEN
         CALL errore(subname,'index out of bound', 10) 
     ENDIF
     !
     res=vcut%corrected(i(1),i(2),i(3))
     !
  ENDIF
  !
END FUNCTION vcut_get

!------------------------------------------
  FUNCTION vcut_spheric_get(vcut,q) RESULT(res)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(IN) :: vcut
  REAL(DP),        INTENT(IN) :: q(3)
  REAL(DP)                    :: res 
  !
  REAl(DP) :: a(3,3), Rcut, kg2 
  LOGICAL  :: limit
  !
  !
  a = vcut%a
  !
  Rcut=0.5*minval(sqrt(sum(a**2,1)))
  Rcut=Rcut-Rcut/50.0
  limit=.false.
  kg2=sum(q**2)
  if(kg2<eps6) then
    limit=.true.
  endif
  if(.not.limit) then
    res=FPI*e2/kg2*(1.0-cos(Rcut*sqrt(kg2)))
  else
    res=FPI*e2*Rcut**2/2.0
  endif
  !
END FUNCTION vcut_spheric_get

!---------------------------------------------------------
  FUNCTION vcut_formula(q,a,b,a_omega,orthorombic) result(res)
  !---------------------------------------------------------
  !
  ! Define the FT of the Coulomb potential according to the
  ! current lattice.
  !
  REAL(DP), INTENT(IN) :: q(3)
  REAL(DP), INTENT(IN) :: a(3,3)
  REAL(DP), INTENT(IN) :: b(3,3)
  REAL(DP), INTENT(IN) :: a_omega
  LOGICAL,  INTENT(IN) :: orthorombic
  REAL(DP)             :: res
  !
  real(dp) :: rwigner
  real(dp) :: sigma

  rwigner=0.5*sqrt(1.0/maxval(sum(b**2,1)))*2*pi

  !
  ! 3.0 is set to give a short range contribution inside the WS cell
  !
  sigma=3.0/rwigner

  ! compute longrange and shortrange contributions
  res=vcut_formula_longrange(q,a,b,a_omega,sigma,6.0D0,orthorombic) &
     +vcut_formula_shortrange(q,sigma)

END FUNCTION vcut_formula

!---------------------------------------------------------
  FUNCTION vcut_formula_longrange(q,a,b,a_omega,sigma,security,orthorombic) result(res)
  !---------------------------------------------------------
  ! compute the longrange contribution
  real(dp), intent(in) :: q(3)
  real(dp), intent(in) :: a(3,3)
  real(dp), intent(in) :: b(3,3)
  real(dp), intent(in) :: a_omega
  real(dp), intent(in) :: sigma
  real(dp), intent(in) :: security ! it determines the grid for the real-space sum; a reasonable value is 4.0
  logical,  intent(in) :: orthorombic
  real(dp)             :: res
  integer :: n1,n2,n3
  integer :: i1,i2,i3
  real(dp) :: d1,d2,d3,weight,factor
  real(dp) :: r(3),r2,modr
  logical :: n1_is_even,n1_is_odd
  real(dp) :: tmp, rtmp(3)
  logical, parameter :: shifted=.false.
  integer :: n1max
  real(dp) :: i1_real,i2_real,i3_real
  real(DP), external :: qe_erf

  n1=security*sqrt(sum(a(:,1)**2))*sigma
  n2=security*sqrt(sum(a(:,2)**2))*sigma
  n3=security*sqrt(sum(a(:,3)**2))*sigma

  n1_is_even=(n1/2)*2==n1
  n1_is_odd=.not.n1_is_even

  d1=1.0/real(n1,dp)
  d2=1.0/real(n2,dp)
  d3=1.0/real(n3,dp)
  res=0.0
  weight=a_omega*d1*d2*d3
! the only symmetry which can be used for any value of q is inversion
! NON-SHIFTED:
!  if n1 is even: loop between 0 and n1/2, with weight=2.0 for all points except 0 and n1max
!  if n2 is odd:  loop between 0 and (n1+1)/2, with weight=2.0 for all points except 0
! SHIFTED:
!  if n1 is even: loop between 0 and n1/2-1, with weight=2.0 for all points
!  if n2 is odd:  loop between 0 and (n1+1)/2, with weight=2.0 for all points except n1max

  if(shifted)then
    if(n1_is_even) n1max=n1/2-1
    if(n1_is_odd)  n1max=(n1+1)/2-1
  else
    if(n1_is_even) n1max=n1/2
    if(n1_is_odd)  n1max=(n1+1)/2
  end if
  do i1=0,n1max
    factor=2.0
    if(shifted) then
      if(n1_is_odd .and. i1==n1max) factor=1.0
    else
      if(i1==0) factor=1.0
      if(n1_is_even .and. i1==n1max) factor=1.0
    end if
    i1_real=i1
    if(shifted) i1_real=i1_real+0.5
    do i2=0,n2-1
      i2_real=i2
      if(shifted) i2_real=i2_real+0.5
      do i3=0,n3-1
        i3_real=i3
        if(shifted) i3_real=i3_real+0.5
        rtmp=matmul(a,(/i1_real*d1,i2_real*d2,i3_real*d3/))
        r=vcut_minimal_image(a,b,rtmp,orthorombic)
        r2=sum(r**2)
        modr=sqrt(r2)
        if(modr*sigma<eps6) then
          tmp=e2*sqrt(2.0/pi)*sigma
        else
          tmp=e2*qe_erf(sigma*sqrt(0.5)*modr)/modr
        end if
        res=res+weight*factor*tmp*cos(sum(r*q))
      end do
    end do
  end do
 END FUNCTION vcut_formula_longrange

!---------------------------------------------------------
FUNCTION vcut_formula_shortrange(q,sigma) result(res)
  !---------------------------------------------------------
  real(dp), intent(in) :: q(3)
  real(dp), intent(in) :: sigma
  real(dp)             :: res
  if(sum(q**2/(sigma*sigma))<eps6) then
! analytic limit for small q
    res=e2*tpi/(sigma*sigma)
  else
    res=e2*fpi/sum(q**2)*(1-exp(-0.5*sum(q**2)/(sigma*sigma)))
  end if
END FUNCTION vcut_formula_shortrange

!---------------------------------------------------------
FUNCTION vcut_minimal_image(a,b,r,orthorombic) result(res)
  !---------------------------------------------------------
  real(dp), intent(in) :: a(3,3)
  real(dp), intent(in) :: b(3,3)
  real(dp), intent(in) :: r(3)
  logical,  intent(in) :: orthorombic
  real(dp)             :: res(3)
  real(dp) :: r_minimal(3)
  real(dp) :: r2_minimal
  real(dp) :: r_try(3)
  real(dp) :: r2_try
  real(dp) :: r_components(3)
  integer :: i1,i2,i3
  integer, parameter :: max_displacement=1
  if(orthorombic) then
! NINT ALGORITHM FOR ORTHOROMBIC CELL
    r_components=(matmul(transpose(b),r))/(2.0*pi)
    r_components=r_components-nint(r_components)
    r_minimal=matmul(a,r_components)
  else
! POOR MAN ALGORITHM FOR GENERIC CELL
    r_minimal=r
    r2_minimal=sum(r_minimal**2)
! loop over the possible neighbours
    do i1=-max_displacement,max_displacement
      do i2=-max_displacement,max_displacement
        do i3=-max_displacement,max_displacement
          if(i1==0 .and. i2==0 .and. i3==0) cycle
            r_try=r+matmul(a,(/i1,i2,i3/))
            r2_try=sum(r_try**2)
            if(r2_try<r2_minimal) then
              r2_minimal=r2_try
              r_minimal=r_try
            endif
        end do
      end do
    end do
  end if
  res=r_minimal
END FUNCTION vcut_minimal_image



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! tools from sax

function num_inverse(a) result(inv)
  real(dp)              :: inv(0:2,0:2)
  real(dp), intent(in)  :: a(0:2,0:2)
  real(dp) :: tmp(0:2,0:2)
  real(dp) :: det
  real(dp),parameter :: eye3(3,3) = reshape((/ 1,0,0,0,1,0,0,0,1/),(/3,3/))
  integer i,j
  do i=0,2
    do j=0,2
      tmp(i,j) = a(modulo(i+1,3),modulo(j+1,3)) * a(modulo(i+2,3),modulo(j+2,3)) &
  &            - a(modulo(i+1,3),modulo(j+2,3)) * a(modulo(i+2,3),modulo(j+1,3))
    end do
  end do
  det = num_determinant(a)
  inv = transpose(tmp) / det
  if(sum((matmul(inv,a))**2-eye3) > 1d-5) then
    write(0,*) "AHIA",sum((matmul(inv,a)-eye3)**2)
    write(0,*) "A",a
    write(0,*) "inv",inv
    write(0,*)">>", matmul(inv,a)
    stop
  end if
end function num_inverse

function num_determinant(a) result(det)
  real(dp), intent(in) :: a(3,3)
  real(dp)             :: det
  det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) &
    - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3) - a(1,3)*a(2,2)*a(3,1)
end function num_determinant

!!! end tools from sax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE coulomb_vcut_module

