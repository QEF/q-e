!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine normalize
  !---------------------------------------------------------------
  !
  !     normalize the US wavefunction so that <phis|S|phis>=1
  !
  !
  !
  use ld1inc
  implicit none

  integer :: &
       n,n1,n2, &   ! counters on beta and mesh function
       ns,nst,ikl  ! counter on wavefunctions

  real(DP) :: &
       work(nwfsx), & ! auxiliary variable for becp
       work1,       & ! the norm
       int_0_inf_dr,& ! integration function
       gi(ndm)       ! used to compute the integrals


  if (pseudotype.ne.3) return 
  !
  !    if US pseudopotential compute the augmentation part
  !
  do ns=1,nwfts
     if (octs(ns).gt.0.0_dp) then
        nst=(llts(ns)+1)*2
        do n1=1,nbeta
           if (llts(ns).eq.lls(n1).and.abs(jjts(ns)-jjs(n1)).lt.1.e-7_dp) then
              ikl=ikk(n1)
              do n=1,ikl
                 gi(n)=betas(n,n1)*phis(n,ns)
              enddo
              work(n1)=int_0_inf_dr(gi,r,r2,dx,ikl,nst)
           else
              work(n1)=0.0_dp
           endif
        enddo
        do n=1,mesh
           gi(n)=phis(n,ns)*phis(n,ns)
        enddo
        work1=int_0_inf_dr(gi,r,r2,dx,mesh,nst)
        !
        !   and adding to the charge density
        !
        do n1=1,nbeta
           do n2=1,nbeta
              work1=work1+qq(n1,n2)*work(n1)*work(n2)  
           enddo
        enddo
        if (abs(work1) < 1e-10_dp) then
             call infomsg('normalize','zero norm: not a true US PP ?',-1)   
             work1=1.0_dp
        else if (work1 < -1e-10_dp) then
             call errore('normalize','negative norm?',1)   
        end if
        work1=sqrt(work1)
        do n=1,mesh
           phis(n,ns)=phis(n,ns)/work1
        enddo
     endif
  enddo

  return
end subroutine normalize
