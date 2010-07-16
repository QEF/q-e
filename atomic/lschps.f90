!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE lschps (mode, z, grid, nin, mch, n, l, e, v, vtau, u, nstop)
  !
  ! integrates radial pauli-type scalar-relativistic equation
  ! on a logarithmic grid
  ! modified routine to be used in finding norm-conserving
  ! pseudopotential
  !
  ! on input:
  !   mode = 1 find energy and wavefunction of bound states,
  !            scalar-relativistic (all-electron)
  !   mode = 2 find energy and wavefunction of bound state,
  !            nonrelativistic (pseudopotentials)
  !   mode = 3 fixed-energy calculation, for logarithmic derivatives
  !   mode = 4 find energy which produces a specified logarithmic
  !            derivative (nonrelativistic, pseudopotentials)
  !   mode = 5 is for pseudopotential to produce wavefunction beyond
  !            radius used for pseudopotential construction
  !   grid =   structure containing radial grid information
  !   z    = atomic number
  !   l, n = main and angular quantum numbers
  !   v(i) = self-consistent potential
  !   vtau(i) = metaGGA potential  D E_XC / D tau
  !   e    = starting estimate of the energy (mode=1,2)
  !          fixed energy at which the wavefctn is calculated (mode=3,4)
  !   mch  = matching at index mch
  !
  ! on output:
  !   e    = final energy (mode=1,2)
  !   u(i) = radial wavefunction (defined as the radial part of the wavefct
  !          multiplied by r)
  !   nstop= 0 if regular termination, 1 otherwise
  !   nin  = last grid point for which the wavefct is calculated?
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: radial_grid_type
  USE ld1inc, ONLY : cau_fact
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT (in) :: mode, n, l
  real(DP), INTENT(in) :: z
  TYPE (radial_grid_type), INTENT(in) :: grid
  real(DP), INTENT(inout) :: v(grid%mesh), vtau(grid%mesh)
  INTEGER, INTENT(inout) :: nin, mch
  real(DP), INTENT(inout) :: e
  INTEGER, INTENT(out) :: nstop
  real (DP), INTENT(out) :: u(grid%mesh)
  !
  ! local variables
  !
  real(DP),PARAMETER:: e2=2.0_dp
  real(DP), EXTERNAL:: aei, aeo, aii, aio
  real(DP):: al, als, ammax,  cn
  real(DP):: de, emax, emin
  real(DP):: eps, fss, gamma, ro, sc
  real(DP):: sls, sn, uld, uout,  upin, upout
  real(DP):: xkap
  INTEGER:: i, it, mmax, nit, node, ierr

  ! these arrays are used as work space
  real(DP),ALLOCATABLE :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)
  !
  !
  nstop=0
  ammax= exp(grid%dx)
  al   = grid%dx
  mmax = grid%mesh

  ALLOCATE(up(mmax), stat=ierr)
  ALLOCATE(upp(mmax), stat=ierr)
  ALLOCATE(cf(mmax), stat=ierr)
  ALLOCATE(dv(mmax), stat=ierr)
  ALLOCATE(fr(mmax), stat=ierr)
  ALLOCATE(frp(mmax), stat=ierr)

  uld=0.0_dp
  v=v/e2
  e=e/e2
  !
  ! convergence factor for solution of schroedinger eq.  if calculated
  ! correction to eigenvalue is smaller in magnitude than eps times
  ! the magnitude of the current guess, the current guess is not changed.
  eps=1.0e-8_dp
  !
  ! relativistic - non-relativistic switch
  !
  IF(mode == 1 .or. mode == 3) THEN
!     fss=(1.0_dp/137.036_dp)**2
     fss=(1.0_dp/cau_fact)**2
     IF(l == 0) gamma=sqrt(1.0_dp-fss*z**2)
     IF(l > 0) gamma=(l*sqrt(l**2-fss*z**2) + &
          (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
  ELSE
     fss=1.0e-20_dp
     gamma=l+1
  ENDIF
  !
  sls=l*(l+1)
  !
  IF(mode == 1 .or. mode == 2) THEN
     emax=v(mmax)+0.5_dp*sls/grid%r(mmax)**2
     emin=0.0_dp
     DO i=1,mmax
        emin=min(emin,v(i)+0.5_dp*sls/grid%r(i)**2)
        !           if (l.eq.0)  write(6,*) grid%r(i),v(i)
     ENDDO
     !         if (l.eq.0) stop
     IF(e > emax) e=1.25_dp*emax
     IF(e < emin) e=0.75_dp*emin
     IF(e > emax) e=0.5_dp*(emax+emin)
  ELSEIF(mode == 4) THEN
     emax=e + 10.0_dp
     emin=e - 10.0_dp
  ENDIF
  !
  DO i=1,4
     u(i)=0.0_dp
     up(i)=0.0_dp
     upp(i)=0.0_dp
  ENDDO
  nit=0
  als=al**2
  !
  ! return point for bound state convergence
10 nit=nit+1
  IF(nit > 60) THEN
     PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
     u=0.0_dp
     nstop=1
     GOTO 999
  ENDIF
  !
  ! coefficient array for u in differential eq.
  DO i=1,mmax
     cf(i)=als*sls + 2.0_dp*als*(v(i)-e)*grid%r(i)**2
  ENDDO
  !
  ! calculate dv/dr for darwin correction
  dv(1)=(-50.0_dp*v(1)+96.0_dp*v(2)-72.0_dp*v(3)+32.0_dp*v(4) &
       -6.0_dp*v(5))/(24.0_dp*al*grid%r(1))
  dv(2)=(-6.0_dp*v(1)-20.0_dp*v(2)+36.0_dp*v(3)-12.0_dp*v(4) &
       +2.0_dp*v(5))/(24.0_dp*al*grid%r(2))
  !
  DO i=3,mmax-2
     dv(i)=(2.0_dp*v(i-2)-16.0_dp*v(i-1)+16.0_dp*v(i+1) &
          -2.0_dp*v(i+2))/(24.0_dp*al*grid%r(i))
  ENDDO
  dv(mmax-1)=( 3.0_dp*v(mmax)+10.0_dp*v(mmax-1)-18.0_dp*v(mmax-2)+ &
       6.0_dp*v(mmax-3)-v(mmax-4))/(12.0_dp*al*grid%r(mmax-1))
  dv(mmax)=( 25.0_dp*v(mmax)-48.0_dp*v(mmax-1)+36.0_dp*v(mmax-2)-&
       16.0_dp*v(mmax-3)+3.0_dp*v(mmax-4))/(12.0_dp*al*grid%r(mmax))
  !
  !  relativistic coefficient arrays for u (fr) and up (frp).
  DO i=1,mmax
     fr(i)=als*(grid%r(i)**2)*(-fss*(v(i)-e)**2 + 0.5_dp*fss*dv(i)/ &
          (grid%r(i)*(1.0_dp+0.5_dp*fss*(e-v(i)))))
     frp(i)=-al*grid%r(i)*0.5_dp*fss*dv(i)/(1.0_dp+0.5_dp*fss*(e-v(i)))
  ENDDO
  !
  ! find classical turning point for matching
  IF(mode == 1 .or. mode == 2) THEN
     DO i=mmax,2,-1
        IF(cf(i-1) <= 0.0_dp .and. cf(i) > 0.0_dp) THEN
           mch=i
           GOTO 40
        ENDIF
     ENDDO
     PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
     e=0.0_dp
     DO i=1,mmax
        u (i)=0.0_dp
     ENDDO
     nstop=1
     GOTO 999
  ELSE
     nin=mch
  ENDIF
40 CONTINUE
  !
  ! start wavefunction with series
  !
  DO i=1,4
     u(i)=grid%r(i)**gamma
     up(i)=al*gamma*grid%r(i)**gamma
     upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
  ENDDO
  !
  ! outward integration using predictor once, corrector
  ! twice
  node=0
  !
  DO i=4,mch-1
     u(i+1)=u(i)+aeo(up,i)
     up(i+1)=up(i)+aeo(upp,i)
     DO it=1,2
        upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
        up(i+1)=up(i)+aio(upp,i)
        u(i+1)=u(i)+aio(up,i)
     ENDDO
     IF(u(i+1)*u(i) <= 0.0_dp) node=node+1
  ENDDO
  !
  uout=u(mch)
  upout=up(mch)
  !
  !
  IF(node-n+l+1 == 0 .or. mode == 3 .or. mode == 5) THEN
     !
     IF(mode == 1 .or. mode == 2) THEN
        !
        ! start inward integration at 10*classical turning
        ! point with simple exponential
        nin=mch+2.3_dp/al
        IF(nin+4 > mmax) nin=mmax-4
        xkap=sqrt(sls/grid%r(nin)**2 + 2.0_dp*(v(nin)-e))
        !
        DO i=nin,nin+4
           u(i)=exp(-xkap*(grid%r(i)-grid%r(nin)))
           up(i)=-grid%r(i)*al*xkap*u(i)
           upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
        ENDDO
        !
        ! integrate inward
        !
        DO i=nin,mch+1,-1
           u(i-1)=u(i)+aei(up,i)
           up(i-1)=up(i)+aei(upp,i)
           DO it=1,2
              upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
              up(i-1)=up(i)+aii(upp,i)
              u(i-1)=u(i)+aii(up,i)
           ENDDO
        ENDDO
        !
        ! scale outside wf for continuity
        sc=uout/u(mch)
        !
        DO i=mch,nin
           up(i)=sc*up(i)
           u (i)=sc*u (i)
        ENDDO
        !
        upin=up(mch)
        !
     ELSE
        !
        upin=uld*uout
        !
     ENDIF
     !
     ! perform normalization sum
     !
     ro=grid%r(1)/sqrt(ammax)
     sn=ro**(2.0_dp*gamma+1.0_dp)/(2.0_dp*gamma+1.0_dp)
     !
     DO i=1,nin-3
        sn=sn+al*grid%r(i)*u(i)**2
     ENDDO
     !
     sn=sn + al*(23.0_dp*grid%r(nin-2)*u(nin-2)**2 &
          + 28.0_dp*grid%r(nin-1)*u(nin-1)**2 &
          +  9.0_dp*grid%r(nin  )*u(nin  )**2)/24.0_dp
     !
     ! normalize u
     cn=1.0_dp/sqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin
     !
     DO i=1,nin
        up(i)=cn*up(i)
        u(i)=cn*u(i)
     ENDDO
     DO i=nin+1,mmax
        u(i)=0.0_dp
     ENDDO
     !
     ! exit for fixed-energy calculation
     !
     IF(mode == 3 .or. mode == 5) GOTO 999

     ! perturbation theory for energy shift
     de=0.5_dp*uout*(upout-upin)/(al*grid%r(mch))
     !
     ! convergence test and possible exit
     !
     IF ( abs(de) < max(abs(e),0.2_dp)*eps) GOTO 999
     !
     IF(de > 0.0_dp) THEN
        emin=e
     ELSE
        emax=e
     ENDIF
     e=e+de
     IF(e > emax .or. e < emin) e=0.5_dp*(emax+emin)
     !
     ! loop back to converge e
     !
     GOTO 10
     !
  ELSEIF(node-n+l+1 < 0) THEN
     ! too few nodes
     emin=e
     e=0.5_dp*(emin+emax)
     GOTO 10
     !
  ELSE
     ! too many nodes
     emax=e
     e=0.5_dp*(emin+emax)
     GOTO 10
  ENDIF
  !
  ! deallocate arrays and exit
  !
999 CONTINUE
  DEALLOCATE(frp)
  DEALLOCATE(fr)
  DEALLOCATE(dv)
  DEALLOCATE(cf)
  DEALLOCATE(upp)
  DEALLOCATE(up)
  e=e*e2
  v=v*e2
  RETURN

END SUBROUTINE lschps
!
FUNCTION aei(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER j
  real(DP):: y(j+3), aei
  !
  aei=-(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j+1) &
       +37.0_dp*y(j+2)-9.0_dp*y(j+3))
  RETURN
END FUNCTION aei
!
! adams extrapolation and interpolation formulas for
! outward and inward integration, abramowitz and
! stegun, p. 896
FUNCTION aeo(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP):: y(j), aeo
  !
  aeo=(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j-1) &
       +37.0_dp*y(j-2)-9.0_dp*y(j-3))
  RETURN
END FUNCTION aeo
!
FUNCTION aii(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP) :: y(j+2), aii
  !
  aii=-(4.16666666667e-2_dp)*(9.0_dp*y(j-1)+19.0_dp*y(j) &
       -5.0_dp*y(j+1)+y(j+2))
  RETURN
END FUNCTION aii
!
FUNCTION aio(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: j
  real(DP):: y(j+1), aio
  !
  aio=(4.16666666667e-2_dp)*(9.0_dp*y(j+1)+19.0_dp*y(j) &
       -5.0_dp*y(j-1)+y(j-2))
  RETURN
END FUNCTION aio
