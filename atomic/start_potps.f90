!
!---------------------------------------------------------------
subroutine start_potps
  !---------------------------------------------------------------
  !
  !     This routine computes an initial estimate of the screening
  !     potential
  !
  !
  !
  use ld1inc
  implicit none

  integer :: &
       ns, &       ! counter on pseudowavefunctions
       is, &       ! counter on spin
       n,  &       ! counter on mesh
       ib,jb,nst,  &   ! counter on lambda
       nnode,  &     ! the number of nodes in lambda
       ik,ikus,lam,nwf0 ! initial phi

  real(kind=dp) ::    &
       xc(8),       & ! coefficients of bessel
       gi(ndm),     & ! auxiliary
       int_0_inf_dr,& ! integral function
       vnew(ndm,2)    ! the potential
  !
  !    compute an initial estimate of the potential
  !
  !
  do ns=1,nwfts
     if (octs(ns).gt.0.0_dp) then
        lam=llts(ns)
        nwf0=nstoae(ns)
        !
        !        compute the ik closer to r_cut
        !
        ik=0
        ikus=0
        do n=1,mesh
           if (r(n).lt.rcutts(ns)) ik=n
           if (r(n).lt.rcutusts(ns)) ikus=n
        enddo
        if (mod(ik,2).eq.0) ik=ik+1
        if (mod(ikus,2).eq.0) ikus=ikus+1
        if (ikus.gt.mesh) &
             call errore('starting potential','ik is wrong ',1)
        !
        !    compute the phi functions
        !
        call compute_phi(lam,ik,nwf0,ns,xc,0,nnode,octs(ns))
        if (pseudotype.eq.3) then
           !
           !   US only on the components where ikus <> ik
           !
           do n=1,mesh
              psipsus(n,ns)=phis(n,ns)
           enddo
           if (ikus.ne.ik) call compute_phius(lam,ikus,ns,xc,0)
        endif
     endif
  enddo

  call normalize
  call chargeps(nwfts,llts,jjts,octs,iswts)
  call new_potential(ndm,mesh,r,r2,sqr,dx,0.0_dp,vxt,lsd,nlcc, &
       latt,enne,rhoc,rhos,vh,vnew)

  do is=1,nspin
     do n=1,mesh
        vpstot(n,is)=vpsloc(n)+vnew(n,is)
        !      if (is.eq.1.and.nspin.eq.1.and.n.lt.420.and.n.gt.410) &
        !                  write(6,'(3f25.16)') r(n), rhos(n,1),vpstot(n,1)
        !      if (is.eq.2.and.nspin.eq.2)  &
        !            write(6,'(3f25.16)') 2.0_dp*rhos(n,1),vpstot(n,1),vpstot(n,2)
     enddo
  enddo
  !
  !    screening the D coefficients
  !
  if (pseudotype.eq.3) then
     do ib=1,nbeta
        do jb=1,ib
           if (lls(ib).eq.lls(jb).and.abs(jjs(ib)-jjs(jb)).lt.1.e-7_dp) then
              lam=lls(ib)
              nst=(lam+1)*2
              do is=1,nspin
                 do n=1,ikk(ib)
                    gi(n)=qvan(n,ib,jb)*vpstot(n,is)
                 enddo
                 ddd(ib,jb,is)= bmat(ib,jb) &
                      + int_0_inf_dr(gi,r,r2,dx,ikk(ib),nst)
                 ddd(jb,ib,is)=ddd(ib,jb,is)
              enddo
           endif
        enddo
     enddo
  endif

  return
end subroutine start_potps
