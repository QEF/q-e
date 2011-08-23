
!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine normalize(phi,l,j,ns)
  !---------------------------------------------------------------
  !
  !     normalize the US wavefunction so that <phis|S|phis>=1
  !
  use kinds, only : dp
  use ld1_parameters, only : nwfsx
  use radial_grids, only : ndmx
  use io_global,    only : stdout
  use ld1inc, only: grid, qq, betas, ikk, lls, jjs, nbeta, pseudotype

  implicit none

  real(DP) ::    &
       phi(ndmx), &  ! function to normalize
       j            ! total angular momentum
  integer ::    &
       l,       &   ! orbital angular momentum
       ns           ! state index
  integer ::    &
       n,n1,n2, &   ! counters on beta and mesh function
       nst,ikl   ! counter on wavefunctions

  real(DP) :: &
       work(nwfsx), & ! auxiliary variable for becp
       work1,       & ! the norm
       int_0_inf_dr,& ! integration function
       gi(ndmx)        ! used to compute the integrals

  if (pseudotype.ne.3) return 
  !
  !    if US pseudopotential compute the augmentation part
  !
  nst=(l+1)*2
  do n1=1,nbeta
     if (l.eq.lls(n1).and.abs(j-jjs(n1)).lt.1.e-7_dp) then
        ikl=ikk(n1)
        do n=1,ikl
           gi(n)=betas(n,n1)*phi(n)
        enddo
        work(n1)=int_0_inf_dr(gi,grid,ikl,nst)
     else
        work(n1)=0.0_dp
     endif
  enddo
  do n=1,grid%mesh
     gi(n)=phi(n)*phi(n)
  enddo
  work1=int_0_inf_dr(gi,grid,grid%mesh,nst)
  !
  !   and adding to the charge density
  !
  do n1=1,nbeta
     do n2=1,nbeta
        work1=work1+qq(n1,n2)*work(n1)*work(n2)  
     enddo
  enddo
  if (abs(work1) < 1e-10_dp) then
     !call infomsg('normalize','zero norm: not a true US PP ?')
     write(stdout,'(7x,a,i3,a,i3,a,f3.1,a)') &
     'Zero norm: self consistency problem; state:',ns,' (l=' ,l,', j=',j,')'
     work1=1.0_dp
  else if (work1 <= -1e-10_dp) then
     call errore('normalize','negative norm?',ns)   
  end if
  work1=sqrt(work1)
  do n=1,grid%mesh
     phi(n)=phi(n)/work1
  enddo

  return
end subroutine normalize
