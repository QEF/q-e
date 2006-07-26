!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine find_qi(logderae,xc,ik,lam,ncn,flag,iok)
  !--------------------------------------------------------------------------
  !
  !      This routine finds three values of q such that the
  !      functions f_l have a logarithmic derivative equal to
  !      logderae at the point ik
  !
  !      if flag=0 f_l=j_l(r)
  !      if flag=1 f_l=r*j_l(r)
  !

  use ld1inc
  implicit none

  integer,parameter ::  ncmax=10   ! maximum allowed nc

  integer ::      &
       ik,    & ! input: the point corresponding to rcut
       ncn,   & ! input: the number of qi to compute
       flag,  & ! input: the type of function
       iok,   & ! output: if 0 the calculation in this routine is ok
       lam      ! input: the angular momentum

  real(DP) :: &
       xc(ncn),& ! output: the values of qi
       logderae  ! input: the logarithmic derivative


  real(DP) ::   &
       j1(ncmax),& ! the bessel function in three points
       qmax,qmin,& ! the limits of the q search
       logdermax,logdermin,& ! the maximum and minimum logder
       logder, & ! the actual logder
       compute_log, &! function for log derivative
       dq      ! the step to braket the q

  integer ::    &
       nc,  &    ! counter on the q found
       imax,&   ! maximum number of iteration to braket
       iq      ! counter on iteration

  iok=0
  if (ncn.gt.ncmax) &
       call errore('find_qi','ncn is too large',1)

  if (flag.eq.0.and.lam.ne.0) &
       call errore('find_qi','lam too large for this iflag',1)

  if (lam.gt.3) &
       call errore('find_qi','l not programmed',1)
  !
  !    fix deltaq and the maximum step number
  !
  dq=0.05_dp
  imax=600
  !
  !    prepare for the first iteration  
  !
  qmax=0.1_dp
  call sph_bes(7,r(ik-3),qmax,lam,j1)
  if (flag.ne.0) then 
     j1(1:7) = j1(1:7)*r(ik-3:ik+3)
  endif
  logdermax=compute_log(j1,r(ik),dx)-logderae

  do nc=1,ncn
     !
     !    bracket the zero
     !
200  qmin=qmax
     logdermin=logdermax
     do iq=1,imax
        xc(nc)=qmin+dq*iq
        call sph_bes(7,r(ik-3),xc(nc),lam,j1)
        if (flag.ne.0) then
           j1(1:7) = j1(1:7)*r(ik-3:ik+3)
        endif
        logdermax=compute_log(j1,r(ik),dx)-logderae
        !
        !    the zero has been bracketed?
        !
        if (logdermax*logdermin.lt.0.0_dp) then
           qmax=xc(nc)
           goto 100
        endif
     enddo
     call infomsg ('find_qi','qmax not found ', -1)
     iok=1
     return
100  continue
     !
     !      start bisection loop
     !
     xc(nc)=(qmax+qmin)/2.0_dp
     call sph_bes(7,r(ik-3),xc(nc),lam,j1)
     if (flag.ne.0) then
        j1(1:7) = j1(1:7)*r(ik-3:ik+3)
     endif
     logder=compute_log(j1,r(ik),dx)-logderae
     if (logder*logdermin.lt.0.0_dp) then
        qmax=xc(nc)
        logdermax=logder
     else
        qmin=xc(nc)
        logdermin=logder
     endif
     !
     !    avoid the asintotes
     !
     if (abs(logdermin-logdermax).gt.1.e3_dp) then 
        qmax=xc(nc)  
        logdermax=logder
        goto 200
     endif
     !
     !    check for convergence
     !
     if (abs(logdermax-logdermin).gt.1.e-8_dp) goto 100
  enddo

  return
end subroutine find_qi

function compute_log(j1,rj,dx)
  use kinds, only : DP
  implicit none
  real(DP) ::   &
       compute_log, &
       deriv_7pts,  &
       dx,          &
       j1(7),       &
       rj(2)    

  !      compute_log=(j1(2)-j1(1))*2.0_dp/( (rj(2)-rj(1))*(j1(2)+j1(1)) )
  compute_log=deriv_7pts(j1,4,rj(1),dx)/j1(4)

  return
end function compute_log
