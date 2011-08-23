!
! Copyright (C) 2002 Vanderbilt group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine has been modified in order to be compatible with the
! ld1 code. The numerical algorithm is unchanged.
! ADC Nov 2003
!
!-------------------------------------------------------------------------
!
subroutine dirsol(idim1,mesh,ncur,lcur,jcur,it,e0,thresh,grid,snl,ruae,nstop)
!
!     subroutine to compute solutions to the full dirac equation
!
!     the dirac equation in rydberg units reads:
!
!     df(r)     k                     alpha
!     ----- = + - f(r) - ( e - v(r) ) ----- g(r)
!      dr       r                       2
!
!     dg(r)     k                          4      alpha
!     ----- = - - g(r) + ( e - v(r) +  -------- ) ----- f(r)
!      dr       r                      alpha**2     2
!
!     where 
!            alpha is the fine structure constant
!            f(r) is r*minor component
!            g(r) is r*major component
!            k is quantum number 
!               k = - (l+1)    if  j = l+0.5
!               k = + l        if  j = l-0.5
!     IMPORTANT: on output, snl(:,1) contains the MAJOR component
!                           snl(:,2) contains the MINOR component
!
!----------------------------------------------------------------------------
!
use io_global, only : stdout
use kinds, only : DP
use radial_grids, only: radial_grid_type
use ld1inc, only : cau_fact, zed
implicit none
integer :: idim1 
type(radial_grid_type),intent(in)::grid
real(DP) :: ruae(idim1),  &   ! the all electron potential
            snl(idim1,2)       ! the wavefunction

real(DP) :: e0,       &     ! the starting energy eigenvalue
                 jcur,     &     ! the j of the state
                 thresh
                  
integer ::  mesh,  &          ! the dimension of the mesh 
            it,    &          ! the iteration
            ncur,  &          ! the n of the state
            lcur              ! the l of the state

real(DP) :: tbya, abyt,        &  
                 emin, emax,        &
                 zz(idim1,2,2),     &
                 tolinf,alpha2,alpha,  &
                 yy(idim1,2),       &
                 vzero,             &
                 f0,f1,f2,g0,g1,g2, &
                 gout, gpout,       &
                 gin, gpin,         &
                 factor,gamma0,     &
                 ecur,              &
                 decur, decurp          
real(DP) :: f(idim1), int_0_inf_dr

integer :: itmax, &     ! maximum number of iterations
           iter,  &     ! current iteration
           ir,    &     ! counter
           ig,    &     ! auxiliary
           kcur,  &     ! current k
           nctp,  &     ! index of the classical turning point
           nodes, &     ! the number of nodes
           ninf,  &     ! practical infinite
           nstop        ! 0 if all ok, 1 otherwise
!
!               r o u t i n e  i n i t i a l i s a t i o n
if (mesh.ne.grid%mesh) call errore('dirsol','mesh dimension is not as expected',1)
nstop=0
!
!     set the maximum number of iterations for improving wavefunctions
!
itmax = 100
!
!     set ( 2 / fine structure constant )
!tbya = 2.0_DP * 137.04_DP
tbya = 2.0_DP * cau_fact

!     set ( fine structure constant / 2 )
abyt = 1.0_DP / tbya

if (jcur.eq.lcur+0.5_DP) then
    kcur = - ( lcur + 1 )
else
    kcur = lcur
endif
!
!       set initial upper and lower bounds for the eigen value
emin = - 1.0e10_DP
emax = 1.0_DP
ecur=e0
!
do iter = 1,itmax
   yy = 0.0_DP
!
!         define the zz array
!         ===================
!
  if ( iter .eq. 1 ) then
    do ir = 1,mesh
       zz(ir,1,1) = grid%rab(ir) * DBLE(kcur) / grid%r(ir)
       zz(ir,2,2) = - zz(ir,1,1)
    enddo
  endif
  do ir = 1,mesh
     zz(ir,1,2) = - grid%rab(ir) * ( ecur - ruae(ir) ) * abyt
     zz(ir,2,1) = - zz(ir,1,2) + grid%rab(ir) * tbya
  enddo
!
!   ==============================================
!   classical turning point and practical infinity
!   ==============================================
!
  do nctp = mesh,10,-1
     if ( zz(nctp,1,2) .lt. 0.0_DP ) goto 240
  enddo
  call errore('dirsol', 'no classical turning point found', 1)
!
!         jump point out of classical turning point loop
240   continue

  if ( nctp .gt. mesh - 10 ) then 
!     write(stdout,*) 'State nlk=', ncur, lcur, kcur, nctp, mesh
!     write(stdout,*) 'ecur, ecurmax=', ecur, ruae(mesh-10)
     write(stdout,*) 'classical turning point too close to mesh',ncur,lcur,kcur
     snl=0.0_DP
     e0=e0-0.2_DP
     nstop=1
     goto 700
  endif
!
  tolinf = log(thresh) ** 2
  do ninf = nctp+10,mesh
     alpha2 = (ruae(ninf)-ecur) * (grid%r(ninf) - grid%r(nctp))**2
     if ( alpha2 .gt. tolinf ) goto 260
  enddo
!
!         jump point out of practical infinity loop
260     continue
!
  if (ninf.gt.mesh) ninf=mesh
!
!         ===========================================================
!         analytic start up of minor and major components from origin
!         ===========================================================
!
!         with finite nucleus so potential constant at origin we have
!
!         f(r) = sum_n f_n r ** ( ig + n )
!         g(r) = sum_n g_n r ** ( ig + n )
!
!         with
!
!         f_n+1 = - (ecur-v(0)) * abyt * g_n / ( ig - kcur + n + 1 )
!         g_n+1 = (ecur-v(0)+tbya**2 ) * abyt * f_n / ( ig + kcur + n + 1)
!
!         if kcur > 0  ig = + kcur , f_0 = 1 , g_0 = 0
!         if kcur < 0  ig = - kcur , f_0 = 0 , g_1 = 1
!
!    Uncomment the following instruction if you need this boundary conditions
!
!  vzero = ruae(1) 
!
!         set f0 and g0
!  if ( kcur .lt. 0 ) then
!     ig = - kcur
!     f0 = 0
!     g0 = 1
!  else
!     ig = kcur
!     f0 = 1
!     g0 = 0
!  endif
! 
!  f1 = - (ecur-vzero) * abyt * g0 / DBLE( ig - kcur + 1 )
!  g1 = (ecur-vzero+tbya**2) * abyt * f0 / DBLE( ig + kcur + 1 )
!  f2 = - (ecur-vzero) * abyt * g1 / DBLE( ig - kcur + 2 )
!  g2 = (ecur-vzero+tbya**2) * abyt * f1 / DBLE( ig + kcur + 2 )
!
!
!  do ir = 1,5
!     yy(ir,1) = grid%r(ir)**ig * ( f0 + grid%r(ir) * ( f1 + grid%r(ir) * f2 ) )
!     yy(ir,2) = grid%r(ir)**ig * ( g0 + grid%r(ir) * ( g1 + grid%r(ir) * g2 ) )
!  enddo

!
!  The following boundary conditions are for the Coulomb nuclear potential
!  ADC 05/10/2007
!
   vzero = ruae(1)+zed*2.0_DP/grid%r(1) 
!
  gamma0=sqrt(kcur**2-4.0_DP*(abyt*zed)**2)
  if ( kcur .lt. 0 ) then
     ig = - kcur
     f0 = (kcur+gamma0)/(2.0_DP*abyt*zed)
     g0 = 1.0_DP
  else
     ig = kcur
     f0 = 1.0_DP
     g0 = (kcur-gamma0)/(2.0_DP*abyt*zed)
  endif
! 
  f1 = - (ecur-vzero) * abyt * g0 / ( gamma0 - kcur + 1.0_DP )
  g1 = (ecur-vzero+tbya**2) * abyt * f0 / ( gamma0 + kcur + 1.0_DP )
  f2 = - (ecur-vzero) * abyt * g1 / ( gamma0 - kcur + 2.0_DP )
  g2 = (ecur-vzero+tbya**2) * abyt * f1 / ( gamma0 + kcur + 2.0_DP )
!
!
  do ir = 1,5
     yy(ir,1) = grid%r(ir)**gamma0*(f0+grid%r(ir)*(f1+grid%r(ir)*f2))
     yy(ir,2) = grid%r(ir)**gamma0*(g0+grid%r(ir)*(g1+grid%r(ir)*g2))
  enddo

!         ===========================
!         outward integration to nctp
!         ===========================
!
!         fifth order predictor corrector integration routine
  call cfdsol(zz,yy,6,nctp,idim1)
!
!         save major component and its gradient at nctp
  gout = yy(nctp,2)
  gpout = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
  gpout = gpout / grid%rab(nctp)
!
!   ==============================================
!   start up of wavefunction at practical infinity
!   ==============================================
!
  do ir = ninf,ninf-4,-1
     alpha = sqrt( ruae(ir) - ecur )
     yy(ir,2) = exp ( - alpha * ( grid%r(ir) - grid%r(nctp) ) )
     yy(ir,1) = ( DBLE(kcur)/grid%r(ir) - alpha ) * yy(ir,2)*tbya / &
  &               ( ecur - ruae(ir) + tbya ** 2 )
  enddo
!
!         ==========================
!         inward integration to nctp
!         ==========================
!
!         fifth order predictor corrector integration routine
  call cfdsol(zz,yy,ninf-5,nctp,idim1)
!
!         save major component and its gradient at nctp
  gin = yy(nctp,2)
  gpin = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
  gpin = gpin / grid%rab(nctp)
!
!
!         ===============================================
!         rescale tail to make major component continuous
!         ===============================================
!
  factor = gout / gin
  do ir = nctp,ninf
     yy(ir,1) = factor * yy(ir,1)
     yy(ir,2) = factor * yy(ir,2)
  enddo
!
  gpin = gpin * factor
!
!         =================================
!         check that the number of nodes ok
!         =================================
!
!         count the number of nodes in major component
  call nodeno(yy(1,2),1,ninf,nodes,idim1)
 
  if ( nodes .lt. ncur - lcur - 1 ) then
!           energy is too low
     emin = ecur
!         write(stdout,*) 'energy too low'
!         write(stdout,'(i5,3f12.5,2i5)') &
!    &         iter,emin,ecur,emax,nodes,ncur-lcur-1
     if ( ecur * 0.9_DP .gt. emax ) then
         ecur = 0.5_DP * ecur + 0.5_DP * emax 
     else
         ecur = 0.9_DP * ecur
     endif
     goto 370
  endif
!
  if ( nodes .gt. ncur - lcur - 1 ) then
!           energy is too high
     emax = ecur
!         
!         write(stdout,*) 'energy too high'
!         write(stdout,'(i5,3f12.5,2i5)') &
!    &         iter,emin,ecur,emax,nodes,ncur-lcur-1
     if ( ecur * 1.1_DP .lt. emin ) then
        ecur = 0.5_DP * ecur + 0.5_DP * emin
     else
        ecur = 1.1_DP * ecur
     endif
     goto 370
  endif
!
!
!         =======================================================
!         find normalisation of wavefunction 
!         =======================================================
!
   do ir = 1,ninf
      f(ir) = (yy(ir,1)**2 + yy(ir,2)**2)
   enddo
   factor=int_0_inf_dr(f,grid,ninf,2*ig)
!
!
!         =========================================
!         variational improvement of the eigenvalue
!         =========================================
!
   decur = gout * ( gpout - gpin ) / factor
!
!         to prevent convergence problems:
!         do not allow decur to exceed 20% of | ecur |
!         do not allow decur to exceed 70% of distance to emin or emax
   if (decur.gt.0.0_DP) then
      emin=ecur
      decurp=min(decur,-0.2_DP*ecur,0.7_DP*(emax-ecur))
   else
      emax=ecur
      decurp=-min(-decur,-0.2_DP*ecur,0.7_DP*(ecur-emin))
   endif
!
!         write(stdout,'(i5,3f12.5,1p2e12.4)') &
!    &         iter,emin,ecur,emax,decur,decurp
!
!         test to see whether eigenvalue converged
   if ( abs(decur) .lt. thresh ) goto 400
 
   ecur = ecur + decurp
!
!         jump point from node check
370  continue
!
!         =======================================================
!         check that the iterative loop is not about to terminate
!         =======================================================
!
   if ( iter .eq. itmax ) then
!           eigenfunction has not converged in allowed number of iterations
            
!      write(stdout,999) it,ncur,lcur,jcur,e0,ecur
!999   format('iter',i4,' state',i4,i4,f4.1,' could not be converged.',/,   &
!    &      ' starting energy for calculation was',f10.5, &
!    &      ' and end value =',f10.5)
       write(stdout,*) 'state nlj',ncur,lcur,jcur, ' not converged'
       snl=0.0_DP
       nstop=1
       goto 700
   endif
!
!       close iterative loop
enddo
!
!       jump point on successful convergence of eigenvalue
400   continue
!
!   normalize the wavefunction and exit. Note also that in this routine
!   yy(.,1) is the small component
!   yy(.,2) is the large component
!   
!      
snl=0.0_DP
do ir=1,mesh
   snl(ir,1)=yy(ir,2)/sqrt(factor)
   snl(ir,2)=yy(ir,1)/sqrt(factor)
enddo
e0=ecur
700 continue
return
end subroutine dirsol
