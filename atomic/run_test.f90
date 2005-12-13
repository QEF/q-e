!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine run_test
  !
  !   This routine is a driver to the tests of the pseudopotential
  !---------------------------------------------------------------
  !

  use ld1inc
  implicit none

  real(DP) :: oc_old(nwfx)

  integer  &
       n, &  ! counter on wavefunctions
       n1,&  ! counter on mesh points
       ir,&  ! counter on mesh points
       im,&  ! position of the maximum
       nc    ! counter on configurations
  integer :: ios
  character(len=1) :: nch
  real(DP) :: dum

  do n=1,nwf
     oc_old(n)=oc(n)
  enddo

  file_tests = trim(prefix)//'.test'
  open(unit=13, file=file_tests, iostat=ios, err=1111, status='unknown')
1111 call errore('ld1_setup','opening file_tests',abs(ios))

  do nc=1,nconf
     if (nconf == 1) then
        file_wavefunctions  = trim(prefix)//'.wfc'
        file_wavefunctionsps= trim(prefix)//'ps.wfc'
        file_logder   = trim(prefix)//'.dlog'
        file_logderps = trim(prefix)//'ps.dlog'
     else
        if (nc < 10) then
           write (nch, '(i1)') nc 
        else
           nch='0'
           call errore ('run_test', &
                'results for some configs not written to file',-1)
        endif
        file_wavefunctions  = trim(prefix)//nch//'.wfc'
        file_wavefunctionsps= trim(prefix)//nch//'ps.wfc'
        file_logder   = trim(prefix)//nch//'.dlog'
        file_logderps = trim(prefix)//nch//'ps.dlog'
     endif
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
           if (els(n1).eq.elts(n).and.rcut(n1).gt.1.e-3_dp) then
              rcutts(n)=rcut(n1)
              rcutusts(n)=rcutus(n1)
              goto 20
           endif
        enddo
        dum=0.0_dp
        im=2
        do ir=1,mesh-1
           dum=abs(psi(ir+1,1,nstoae(n)))
           if(dum.gt.abs(psi(ir,1,nstoae(n)))) im=ir+1
        enddo
        if (pseudotype.lt.3) then
           rcutts(n)=r(im)*1.1_dp
           rcutusts(n)=r(im)*1.1_dp
        else
           if (ll(nstoae(n)).eq.0) then
              rcutts(n)=r(im)*1.6_dp
              rcutusts(n)=r(im)*1.7_dp
           elseif (ll(nstoae(n)).eq.1) then
              rcutts(n)=r(im)*1.6_dp
              rcutusts(n)=r(im)*1.7_dp
              if (el(nstoae(n)).eq.'2P') then
                 rcutts(n)=r(im)*1.7_dp
                 rcutusts(n)=r(im)*1.8_dp
              endif
           elseif (ll(nstoae(n)).eq.2) then
              rcutts(n)=r(im)*2.0_dp
              rcutusts(n)=r(im)*2.2_dp
              if (el(nstoae(n)).eq.'3D') then
                 rcutts(n)=r(im)*2.5_dp
                 rcutusts(n)=r(im)*3.0_dp
              endif
           endif
        endif
20   enddo
     !
     !   and run the pseudopotential test
     !
     call run_pseudo
     !
     if (nc.eq.1) etots0=etots
     !
     !   print results
     !
     call write_resultsps 
     !
  enddo
  close (unit = 13)  

  return
end subroutine run_test
