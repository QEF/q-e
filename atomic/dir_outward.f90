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
subroutine dir_outward(idim1,mesh,lcur,jcur,e0,dx,snl,r,rab,ruae)
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
use kinds, only : DP
use ld1inc, only : cau_fact
implicit none
integer :: idim1 
real(DP) :: r(idim1),     &   ! the radial mesh
                 rab(idim1),   &   ! derivative of the radial mesh
                 ruae(idim1),  &   ! the all electron potential
                 snl(idim1,2)       ! the wavefunction

real(DP) :: e0,       &     ! the starting energy eigenvalue
                 dx,       &     ! dx mesh value
                 jcur            ! the j of the state
                  
integer ::  mesh,  &          ! the dimension of the mesh 
            lcur              ! the l of the state

real(DP) :: tbya, abyt,        &  
                 zz(idim1,2,2),     &
                 tolinf,alpha2,alpha,  &
                 yy(idim1,2),       &
                 vzero,             &
                 f0,f1,f2,g0,g1,g2, &
                 ecur                
real(DP) :: r2(idim1), f(idim1), int_0_inf_dr

integer :: ir,    &     ! counter
           ig,    &     ! auxiliary
           kcur         ! current k
!
!               r o u t i n e  i n i t i a l i s a t i o n
!
do ir=1,mesh
   ruae(ir)=ruae(ir)*r(ir)
enddo
!
!     set ( 2 / fine structure constant )
!tbya = 2.0_DP * 137.04_DP
tbya = 2.0_DP * cau_fact
!     set ( fine structure constant / 2 )
abyt = 1.0_DP / tbya

r2=r**2

if (jcur.eq.lcur+0.5_DP) then
    kcur = - ( lcur + 1 )
else
    kcur = lcur
endif
!
!       set initial upper and lower bounds for the eigen value
ecur=e0
!
yy = 0.0_DP
!
!         define the zz array
!         ===================
!
do ir = 1,mesh
   zz(ir,1,1) = rab(ir) * DBLE(kcur) / r(ir)
   zz(ir,2,2) = - zz(ir,1,1)
   zz(ir,1,2) = - rab(ir) * ( ecur - ruae(ir) / r(ir) ) * abyt
   zz(ir,2,1) = - zz(ir,1,2) + rab(ir) * tbya
enddo
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
vzero = ruae(1) / r(1)
!
!         set f0 and g0
if ( kcur .lt. 0 ) then
   ig = - kcur
   f0 = 0
   g0 = 1
else
   ig = kcur
   f0 = 1
   g0 = 0
endif
 
f1 = - (ecur-vzero) * abyt * g0 / DBLE( ig - kcur + 1 )
g1 = (ecur-vzero+tbya**2) * abyt * f0 / DBLE( ig + kcur + 1 )
f2 = - (ecur-vzero) * abyt * g1 / DBLE( ig - kcur + 2 )
g2 = (ecur-vzero+tbya**2) * abyt * f1 / DBLE( ig + kcur + 2 )
!
!
do ir = 1,5
   yy(ir,1) = r(ir)**ig * ( f0 + r(ir) * ( f1 + r(ir) * f2 ) )
   yy(ir,2) = r(ir)**ig * ( g0 + r(ir) * ( g1 + r(ir) * g2 ) )
enddo

!         ===========================
!         outward integration to mesh
!         ===========================
!
!         fifth order predictor corrector integration routine
call cfdsol(zz,yy,6,mesh,idim1)
!
!         =======================================================
!         copy the wavefunction 
!         =======================================================
!
snl=0.0_DP
do ir=1,mesh
   snl(ir,2)=yy(ir,1)
   snl(ir,1)=yy(ir,2)
enddo
do ir=1,mesh
   ruae(ir)=ruae(ir)/r(ir)
enddo
return
end subroutine dir_outward
