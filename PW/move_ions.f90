!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine move_ions  
  !-----------------------------------------------------------------------
  !
  !     This routine moves the ions according to the requested scheme:
  !
  !     iswitch = 1      bfgs minimizations or conjugate gradient
  !     iswitch = 2      constrained bfgs minimization:
  !                      the user must supply the routine 'constrain' which
  !                      defines the constraint equation and the gradient
  !                      the constraint function gv(tau), dgv(i,tau) such
  !                      that:
  !
  !                                gv({tau}) - target = 0,
  !
  !                      and
  !
  !                                             D gv( {tau} )
  !                                dgv(i,na) = ---------------.
  !                                             D tau(i,na)
  !
  !     iswitch = 3      molecular dynamics, (verlet of vcsmd)
  !     iswitch = 4      molecular dynamics with one constraint,
  !                      the same conventions as iswitch = 2
  !
#include "machine.h"
  use pwcom  
  implicit none
  real(kind=DP) :: dummy, gv
  real(kind=DP) , allocatable :: dgv (:,:)
  real(kind=DP) :: dgv2, theta0  
  ! auxiliar variable
  ! gv=0 defines the constrain
  ! the gradient of gv
  ! its square modulus
  ! the value of the one-dimensional const

  conv_ions = .false.  
  if (iswitch.eq.2.or.iswitch.eq.4) then  
     allocate ( dgv(3,nat) ) 
     !
     !     gv is the function which define the constrain, now first of all we
     !     find the constrain theta0 such that gv=0, and find the gradient of
     !     gv, dgv
     !
     dummy = 0.d0  
     call constrain (theta0, gv, dgv, dgv2, dummy, nat, tau, alat)  
     !
     !     find the constrained forces
     !
     call new_force (dgv, dgv2)  
     deallocate (dgv)  
  endif
  !
  !     do the minimization / dynamics step
  !

  if (lmovecell.and. (iswitch.eq.2.or.iswitch.eq.4) ) call error ( &
       'move_ions', 'variable cell and constrain not implemented', 1)

  if (iswitch.eq.1.or.iswitch.eq.2) call bfgs  
  if (iswitch.eq.3.or.iswitch.eq.4) then  
     if (calc.eq.' ') call dynamics  ! verlet dynamics
     if (calc.ne.' ') call vcsmd     ! variable cell shape md
  endif
  if (iswitch.gt.4 .or. iswitch.le.0) then  
     call error ('move_ions', 'iswitch value not implemented or wrong', 1)  
  endif
  !
  !     check if the new positions satisfy the constrain equation, in
  !     the CP case this is done inside the routine "cp"
  !

  if (iswitch.eq.2.or.iswitch.eq.4) call check_constrain (alat, tau, atm, &
       ityp, theta0, nat)
  !
  !     before leaving check that the new positions still transform
  !     according to the symmetry of the system.
  !

  call checkallsym (nsym, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, irt, ftau)
  return  

end subroutine move_ions
!-------------------------------------------------------------------
subroutine new_force (dg, dg2)  
!-------------------------------------------------------------------
!
!     find the lagrange multiplier lambda for the problem with one const
!
!                force*dg
!     lambda = - --------,
!                 |dg|^2
!
!     and redefine the forces:
!
!     force = force + lambda*dg
!
!     where dg is the gradient of the constraint function
!
use pwcom  
integer :: na, i, ipol  

real(kind=DP) :: dg (3, nat), lambda, dg2, DDOT, sum  

lambda = 0.d0  
if (dg2.ne.0.d0) then  
   lambda = - DDOT (3 * nat, force, 1, dg, 1) / dg2  
   call DAXPY (3 * nat, lambda, dg, 1, force, 1)  
   if (DDOT (3 * nat, force, 1, dg, 1) **2.gt.1.d-30) then  
call error ('new_force', 'force is not orthogonal to constrain', - 1)
      print *, DDOT (3 * nat, force, 1, dg, 1) **2  
   endif  
   do ipol = 1, 3  
   sum = 0.d0  
   do na = 1, nat  
   sum = sum + force (ipol, na)  
   enddo  
!
!     impose total force = 0
!
   do na = 1, nat  
   force (ipol, na) = force (ipol, na) - sum / nat  
   enddo  
   enddo  
!
! resymmetrize (should not be needed, but...)
!
   if (nsym.gt.1) then  
      do na = 1, nat  
      call trnvect (force (1, na), at, bg, - 1)  
      enddo  
      call symvect (nat, force, nsym, s, irt)  
      do na = 1, nat  
      call trnvect (force (1, na), at, bg, 1)  
      enddo  
   endif  
   write (6, '(/5x,"Constrained forces")')  
   do na = 1, nat  
   write (6, '(3f14.8)') (force (i, na) , i = 1, 3)  
   enddo  

endif  
return  

end subroutine new_force
!---------------------------------------------------------------------

subroutine check_constrain (alat, tau, atm, ityp, theta0, nat)  
  !---------------------------------------------------------------------
  !
  !     update tau so that the constraint equation g=0 is satisfied,
  !     use the recursion formula:
  !
  !                      g(tau)
  !     tau' = tau -  ------------ * dg(tau)
  !                    |dg(tau)|^2
  !
  !     in normal cases the constraint equation should be always satisfied
  !     the very first iteration.
  !
  use parameters
  implicit none  
  integer :: ityp ( * ), nat, na, i, maxiter  
  character(len=3 ) ::  atm(*) 

  real(kind=DP) :: tau (3, nat)
  real(kind=DP), allocatable :: dg (:,:)
  real(kind=DP) :: alat, dg2, g, theta0, dummy, eps

  parameter (eps = 1.d-15, maxiter = 250)  
  allocate ( dg(3,nat) ) 
  call constrain (dummy, g, dg, dg2, theta0, nat, tau, alat)  
  write (6, '(5x,"G = ",1pe9.2," iteration # ",i3)') g, 0  
  do i = 1, maxiter  
     !
     ! check if g=0
     !
     if (abs (g) .lt.eps) goto 14  
     !
     ! if g<>0 find new tau = tau - g*dg/dg2 and check again
     !
     call DAXPY (3 * nat, - g / dg2, dg, 1, tau, 1)  
     call constrain (dummy, g, dg, dg2, theta0, nat, tau, alat)  
     write (6, '(5x,"G = ",1pe9.2," iteration # ",i3)') g, i  
  enddo
  call error ('new_dtau', 'g=0 is not satisfied g=', - 1)  
14 continue  
  !     write(6,'(5x,"G = ",1pe9.2)')g
  write (6, '(5x,"Number of step(s): ",i3)') i - 1  
  !
  !     if the atomic positions have been corrected write them on output
  !
  if (i.gt.1) then  
     write (6, '(/5x,"Corrected atomic positions:",/)')  
     do na = 1, nat  
        write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau(i,na), i=1,3)
     enddo

  endif
  deallocate (dg)  
  return  
end subroutine check_constrain

