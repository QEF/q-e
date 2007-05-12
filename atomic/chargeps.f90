!
! Copyright (C) 2007 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine chargeps(rho_i,phi_i,nwf_i,ll_i,jj_i,oc_i,iswf_i)
  !---------------------------------------------------------------
  !
  !   calculate the (spherical) pseudo charge density 
  !
  use ld1inc

  integer :: &
       nwf_i,        & ! input: the number of wavefunctions
       iswf_i(nwfsx),& ! input: their spin
       ll_i(nwfsx)     ! input: their angular momentum

  real(DP) ::  &
       jj_i(nwfsx), & ! input: their total angular momentum
       oc_i(nwfsx), & ! input: the occupation
       phi_i(ndm,nwfsx), & ! input: the functions to add
       rho_i(ndm,2)   ! output: the (nspin) components of the charge

  integer ::     &
       is,     &   ! counter on spin
       n,n1,n2,&   ! counters on beta and mesh function
       ns,nst,ikl  ! counter on wavefunctions

  real(DP) ::    &
       work(nwfsx), & ! auxiliary variable for becp
       int_0_inf_dr,& ! integration function
       gi(ndm)        ! used to compute the integrals


  rho_i=0.0_dp
  !
  !    compute the square modulus of the eigenfunctions
  !
  do ns=1,nwf_i
     if (oc_i(ns).gt.0.0_dp) then
        is=iswf_i(ns)
        do n=1,mesh
           rho_i(n,is)=rho_i(n,is)+oc_i(ns)*phi_i(n,ns)**2
        end do
     endif
  enddo
  !
  !    if US pseudopotential compute the augmentation part
  !
  if (pseudotype.eq.3) then
     do ns=1,nwf_i
        if (oc_i(ns).gt.0.0_dp) then
           is=iswf_i(ns)
           do n1=1,nbeta
              if (ll_i(ns).eq.lls(n1).and. &
                   abs(jj_i(ns)-jjs(n1)).lt.1.e-7_dp) then
                 nst=(ll_i(ns)+1)*2
                 ikl=ikk(n1)
                 do n=1,ikl
                    gi(n)=betas(n,n1)*phi_i(n,ns)
                 enddo
                 work(n1)=int_0_inf_dr(gi,r,r2,dx,ikl,nst)
              else
                 work(n1)=0.0_dp
              endif
           enddo
           !
           !   and adding to the charge density
           !
           do n1=1,nbeta
              do n2=1,nbeta
                 do n=1,mesh
                    rho_i(n,is)=rho_i(n,is)+qvan(n,n1,n2)*oc_i(ns)* &
                         work(n1)*work(n2)
                 enddo
              enddo
           enddo
        endif
     enddo
  endif

  return
end subroutine chargeps
