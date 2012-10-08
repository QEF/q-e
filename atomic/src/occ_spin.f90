!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine occ_spin(nwf,nwfx,el,nn,ll,oc,isw)
  !---------------------------------------------------------------
  !
  !  This routine splits the occupations of the states between spin-up
  !  and spin down. If the occupations are lower than 2*l+1 it does
  !  nothing, otherwise 2*l+1 states are assumed with spin up and
  !  the difference with spin down. 
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx)
  character(len=2) :: el(nwfx)

  integer :: nwf0, n, n1
  logical :: ok

  nwf0=nwf
  do n=1,nwf0
     if (oc(n) > (2*ll(n)+1)) then
        !
        !    check that the new state is not already available
        !
        do n1=n+1,nwf0
           if (el(n1)==el(n)) call errore('occ_spin','wrong occupations',1)
        enddo
        !
        !    and add it
        !
        nwf=nwf+1
        if (nwf > nwfx) call errore('occ_spin','too many wavefunctions',1)
        el(nwf)=el(n)
        nn(nwf)=nn(n)
        ll(nwf)=ll(n)
        oc(nwf)=oc(n)-2*ll(n)-1
        oc(n)=2*ll(n)+1
        if (isw(n) == 1) isw(nwf)=2 
        if (isw(n) == 2) isw(nwf)=1 
     else
        ok=.true.
        do n1=1,nwf0
           if (n1 /= n) ok=ok.and.(el(n1) /= el(n))  
        enddo
        if (ok) then
           nwf=nwf+1
           if (nwf > nwfx) &
                & call errore('occ_spin','too many wavefunctions',1)
           el(nwf)=el(n)
           nn(nwf)=nn(n)
           ll(nwf)=ll(n)
           oc(nwf)=0.0_dp
           if (oc(n)<0.0_dp) oc(nwf)=oc(n)
           if (isw(n) == 1) isw(nwf)=2 
           if (isw(n) == 2) isw(nwf)=1 
        endif
     endif
  enddo
  return
end subroutine occ_spin
!
subroutine occ_spinorb(nwf,nwfx,el,nn,ll,jj,oc,isw,rel_dist)
  !
  ! This subroutine splits the states according to their j as needed to
  ! make a spin orbit calculation.
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx)
  character(len=2) :: el(nwfx)
  character(len=20) :: rel_dist

  integer :: nwftot, nwf0, n, m
  logical :: ok

  nwftot=0
  do n=1,nwf
     nwftot=nwftot+1
     if (ll(n).gt.0 .and. abs(jj(n)).lt.1.e-2_dp) nwftot=nwftot+1
  enddo
  ok=.true.
  nwf0=nwf
  do n=1,nwftot
     if (ok) then
        if (ll(n).gt.0.and.abs(jj(n)).lt.1.e-2_dp) then
           ok=.false.
           jj(n)=ll(n)-0.5_dp
           nwf0=nwf0+1
           if (nwf0.gt.nwfx) call errore('occ_spinorb','too many wfc',1)
           do m=nwf0-1,n+1,-1
              nn(m+1)=nn(m)
              ll(m+1)=ll(m)
              jj(m+1)=jj(m)
              el(m+1)=el(m)
              isw(m+1)=isw(m)
              oc(m+1)=oc(m)
           enddo
           nn(n+1)=nn(n)
           ll(n+1)=ll(n)
           jj(n+1)=ll(n)+0.5_dp  
           el(n+1)=el(n)
           isw(n+1)=isw(n)
           if (oc(n).gt.-1.e-3_dp) then
              if (rel_dist=='average'.or.rel_dist=='AVERAGE' &
                          .or.rel_dist=='Average') then
                 oc(n+1)=oc(n)*(2.0_DP*(ll(n)+1))/(2.0_DP*(2*ll(n)+1)) 
                 oc(n)=oc(n)*(2.0_DP*ll(n))/(2.0_DP*(2*ll(n)+1)) 
              else
                 oc(n+1)=max(0.0_dp,oc(n)-min(2.0_dp*ll(n),oc(n)))
                 oc(n)=min(2.0_dp*ll(n),oc(n))
              endif
           else
              oc(n+1)=oc(n)
           endif
        else
           if (abs( jj(n) ).lt. 1.e-2_dp) jj(n)=0.5_dp
        endif
     else
        ok=.true.
     endif
  enddo
  nwf=nwf0
  return
end subroutine occ_spinorb
!---------------------------------------------------------------
subroutine occ_spinorbps(nwf,nwfx,el,nn,ll,jj,oc,rcut,rcutus,enls,isw,rel_dist)
!---------------------------------------------------------------
  !
  ! This subroutine splits the states according to their j as needed to
  ! make a spin orbit calculation.
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx), rcut(nwfx), rcutus(nwfx), enls(nwfx)
  character(len=2) :: el(nwfx)
  character(len=20) :: rel_dist

  integer :: nwftot, nwf0, n, m
  logical :: ok

  nwftot=0
  do n=1,nwf
     nwftot=nwftot+1
     if (ll(n).gt.0.and.abs(jj(n)).lt.1.e-2_dp) nwftot=nwftot+1
  enddo
  ok=.true.
  nwf0=nwf
  do n=1,nwftot
     if (ok) then
        if (ll(n).gt.0.and.abs(jj(n)).lt.1.e-2_dp) then
           ok=.false.
           jj(n)=ll(n)-0.5_dp
           nwf0=nwf0+1
           if (nwf0.gt.nwfx) call errore('occ_spinorbps','too many wfc',1)
           do m=nwf0-1,n+1,-1
              nn(m+1)=nn(m)
              ll(m+1)=ll(m)
              jj(m+1)=jj(m)
              el(m+1)=el(m)
              isw(m+1)=isw(m)
              oc(m+1)=oc(m)
              rcut(m+1)=rcut(m)
              rcutus(m+1)=rcutus(m)
              enls(m+1)=enls(m)
           enddo
           nn(n+1)=nn(n)
           ll(n+1)=ll(n)
           jj(n+1)=ll(n)+0.5_dp
           el(n+1)=el(n)
           rcut(n+1)=rcut(n)
           rcutus(n+1)=rcutus(n)
           enls(n+1)=enls(n)
           isw(n+1)=isw(n)
           if (rel_dist=='average'.or.rel_dist=='AVERAGE' &
                          .or.rel_dist=='Average') then
              oc(n+1)=oc(n)*(2.0_DP*(ll(n)+1))/(2.0_DP*(2*ll(n)+1)) 
              oc(n)=oc(n)*(2.0_DP*ll(n))/(2.0_DP*(2*ll(n)+1)) 
           else
              oc(n+1)=max(0.0_dp,oc(n)-min(2.0_dp*ll(n),oc(n)))
              oc(n)=min(2.0_dp*ll(n),oc(n))
           endif
        else
           if (abs(jj(n)).lt.1.e-2_dp) jj(n)=0.5_dp
        endif
     else
        ok=.true.
     endif
  enddo
  nwf=nwf0
  return
end subroutine occ_spinorbps

!---------------------------------------------------------------
subroutine occ_spin_tot(nwf,nwfx,el,nn,ll,oc,isw,enl,psi)
  !---------------------------------------------------------------
  !
  !  This routine is similar to occ_spin, but updates also energies
  !  and wavefunctions for frozen_core calculations
  !
  use radial_grids, only : ndmx
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), enl(nwfx), psi(ndmx,2,nwfx)
  character(len=2) :: el(nwfx)

  integer :: nwf0, n, n1
  logical :: ok

  nwf0=nwf
  do n=1,nwf0
     if (oc(n) > (2*ll(n)+1)) then
        !
        !    check that the new state is not already available
        !
        do n1=n+1,nwf0
           if (el(n1)==el(n)) call errore('occ_spin_tot','wrong occupations',1)
        enddo
        !
        !    and add it
        !
        nwf=nwf+1
        if (nwf > nwfx) call errore('occ_spin_tot','too many wavefunctions',1)
        el(nwf)=el(n)
        nn(nwf)=nn(n)
        ll(nwf)=ll(n)
        oc(nwf)=oc(n)-2*ll(n)-1
        oc(n)=2*ll(n)+1
        if (isw(n) == 1) isw(nwf)=2 
        if (isw(n) == 2) isw(nwf)=1 
        enl(nwf)=enl(n)
        psi(:,1,nwf)=psi(:,1,n)
     else
        ok=.true.
        do n1=1,nwf0
           if (n1 /= n) ok=ok.and.(el(n1) /= el(n))  
        enddo
        if (ok) then
           nwf=nwf+1
           if (nwf > nwfx) &
                & call errore('occ_spin_tot','too many wavefunctions',1)
           el(nwf)=el(n)
           nn(nwf)=nn(n)
           ll(nwf)=ll(n)
           oc(nwf)=0.0_dp
           if (oc(n)<0.0_DP) oc(nwf)=oc(n)
           if (isw(n) == 1) isw(nwf)=2 
           if (isw(n) == 2) isw(nwf)=1 
           enl(nwf)=enl(n)
           psi(:,1,nwf)=psi(:,1,n)
        endif
     endif
  enddo
  return
end subroutine occ_spin_tot
!
