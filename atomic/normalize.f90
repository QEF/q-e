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

  real(kind=dp) :: &
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
        if (work1.lt.1e-10_dp)  &
             call errore('normalize','negative or zero norm?',1)   
        work1=sqrt(work1)
        do n=1,mesh
           phis(n,ns)=phis(n,ns)/work1
        enddo
     endif
  enddo

  return
end subroutine normalize
