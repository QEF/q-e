!
!---------------------------------------------------------------
subroutine chargeps(nwff,lli,jji,oci,iswfi)
  !---------------------------------------------------------------
  !
  !   calculate the (spherical) pseudo charge density and
  !   spin polarization
  !
  use ld1inc

  integer :: &
       nwff,       & ! input: the number of wavefunctions
       iswfi(nwff),& ! input: their spin
       lli(nwff)     ! input: their angular momentum

  real(kind=dp) ::  &
       jji(nwff), & ! input: their total angular momentum
       oci(nwff)    ! input: the occupation

  integer ::     &
       is,     &   ! counter on spin
       n,n1,n2,&   ! counters on beta and mesh function
       ns,nst,ikl  ! counter on wavefunctions

  real(kind=dp) ::    &
       work(nwfsx), & ! auxiliary variable for becp
       int_0_inf_dr,& ! integration function
       gi(ndm)        ! used to compute the integrals


  rhos=0.0_dp
  !
  !    compute the square modulus of the eigenfunctions
  !
  do ns=1,nwff
     if (oci(ns).gt.0.0_dp) then
        is=iswfi(ns)
        do n=1,mesh
           rhos(n,is)=rhos(n,is)+oci(ns)*phis(n,ns)**2
        end do
     endif
  enddo
  !
  !    if US pseudopotential compute the augmentation part
  !
  if (pseudotype.eq.3) then
     do ns=1,nwff
        if (oci(ns).gt.0.0_dp) then
           is=iswfi(ns)
           do n1=1,nbeta
              if (lli(ns).eq.lls(n1).and. &
                   abs(jji(ns)-jjs(n1)).lt.1.e-7_dp) then
                 nst=(lli(ns)+1)*2
                 ikl=ikk(n1)
                 do n=1,ikl
                    gi(n)=betas(n,n1)*phis(n,ns)
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
                    rhos(n,is)=rhos(n,is)+qvan(n,n1,n2)*oci(ns)* &
                         work(n1)*work(n2)
                 enddo
              enddo
           enddo
        endif
     enddo
  endif

  return
end subroutine chargeps
