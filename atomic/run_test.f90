subroutine run_test
!
!   This routine is a driver to the tests of the pseudopotential
!

use ld1inc
implicit none

integer  &
      n, &  ! counter on wavefunctions
      n1,&  ! counter on mesh points
      ir,&  ! counter on mesh points
      im,&  ! position of the maximum
      nc    ! counter on configurations

real(kind=dp) :: dum

do nc=1,nconf
   nwfts=nwftsc(nc)
   do n=1,nwf
      oc(n)=oc_old(n)
   enddo
   do n=1,nwfts
      nnts(n)=nntsc(n,nc)
      llts(n)=lltsc(n,nc)
      elts(n)=eltsc(n,nc)
      jjts(n)=jjtsc(n,nc)
      iswts(n)=iswtsc(n,nc)
      octs(n)=octsc(n,nc)
      if (rel==2) jjts(n)=jjtsc(n,nc)
      nstoae(n)=nstoaec(n,nc)
      oc(nstoae(n))=octs(n)
   enddo
   call all_electron(.true.)
   if (nc.eq.1) etot0=etot
!
!   choose the cut-off radius for the initial estimate of the wavefunctions
!   find the maximum of the all electron wavefunction
!
   do n=1,nwfts
      do n1=1,nbeta
         if (els(n1).eq.elts(n).and.rcut(n1).gt.1.d-3) then
            rcutts(n)=rcut(n1)
            rcutusts(n)=rcutus(n1)
            goto 20
         endif
      enddo
      dum=0.d0
      do ir=1,mesh-1
         dum=abs(psi(ir+1,nstoae(n)))
         if(dum.gt.abs(psi(ir,nstoae(n)))) im=ir+1
      enddo
      if (pseudotype.lt.3) then
        rcutts(n)=r(im)*1.1d0
        rcutusts(n)=r(im)*1.1d0
      else
        if (ll(nstoae(n)).eq.0) then
           rcutts(n)=r(im)*1.6d0
           rcutusts(n)=r(im)*1.7d0
        elseif (ll(nstoae(n)).eq.1) then
           rcutts(n)=r(im)*1.6d0
           rcutusts(n)=r(im)*1.7d0
        elseif (ll(nstoae(n)).eq.2) then
           rcutts(n)=r(im)*2.0d0
           rcutusts(n)=r(im)*2.2d0
           if (el(nstoae(n)).eq.'3D') then
              rcutts(n)=r(im)*2.5d0
              rcutusts(n)=r(im)*3.0d0
           endif
        endif
      endif
20 enddo
!
!   and run the pseudopotential test
!
   call run_pseudo
   if (nc.eq.1) etots0=etots
enddo

return
end
