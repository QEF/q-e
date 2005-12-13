!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ld1_setup
  !---------------------------------------------------------------
  !
  !     this routine computes the general variables needed in
  !     the atomic calculation
  !
  !
  use ld1inc
  use funct, only : get_iexch, dft_is_meta !, set_dft_from_name
  implicit none

  integer n, n1, nc, m, nwftot, ios
  logical ok, hf, oep
  !
  !     transform dft in a series of codes for the exchange and
  !     correlation routine
  !
!  if (iswitch /= 2 .or. pseudotype == 1) call set_dft_from_name(dft)
  !
  if (dft_is_meta()) call errore('setyp','meta-GGA not implemented yet', 1)
  hf  = get_iexch().eq.5
  if (hf)     call errore('setup','HF not implemented yet',1)
  oep = get_iexch().eq.4
  if (oep.and.iswitch>1) &
     call errore('setup','OEP is implemented only for all-electron calc.',1)
  if (oep.and.rel>0) &
     call errore('setup','OEP is implemented only for non-relativistic calc.',1)
  !
  call set_sl3(sl3,lmx)
  !
  !  make the correspondence all-electron pseudopotential
  !
  if (iswitch >= 2) then
     do nc=1, nconf
        do n=1,nwftsc(nc)
           nstoaec(n,nc)=0
           do n1=1,nwf
              if (lsd.eq.1) then
                 if (eltsc(n,nc).eq.el(n1) &
                      .and.iswtsc(n,nc).eq.isw(n1)) then
                    nstoaec(n,nc)=n1
                 endif
              else
                 if (rel==2) then
                    if (eltsc(n,nc).eq.el(n1).and.jjtsc(n,nc).eq.jj(n1)) &
                         & nstoaec(n,nc)=n1
                 else
                    if (eltsc(n,nc).eq.el(n1))  nstoaec(n,nc)=n1
                 endif
              endif
           enddo
           if (nstoaec(n,nc).eq.0) call errore('ld1_setup', &
                'all electron wfc corresponding to pseudo-state ' &
          &     //eltsc(n,nc)//' not found',nc)
        enddo
     enddo
  endif
  !
  new(:)=.false.
  !
  !  divide the core and valence states
  !
  if (iswitch == 3) then 
     isws=1
     do n=1,nwf
        core_state(n)=.true.
     enddo
     do n=1,nwfs
        nstoae(n)=0
        do n1=1,nwf
           if (rel==2) then
              if (els(n).eq.el(n1).and.jjs(n).eq.jj(n1)) then
                 nstoae(n)=n1
                 core_state(n1)=.false.
              endif
           else
              if (els(n).eq.el(n1)) then
                 nstoae(n)=n1
                 core_state(n1)=.false.
              endif
           endif
        enddo
        if (nstoae(n).eq.0) call errore('ld1_setup', &
             'no all electron for this ps',n)
        if (enls(n).ne.0.0_dp) then
           new(n)=.true.
        else
           new(n)=.false.
        endif
     enddo
     if (lloc > -1) then
        nsloc=nwfs
        nbeta=nwfs-1
        if (rel==2.and.lloc.ne.0) then
           nsloc=nwfs-1
           nbeta=nwfs-2
           if (lls(nsloc+1).ne.lloc) &
             call errore('ld1_setup','mismatch between lloc and l of ' // &
           &           'spin-orbit split wfc chosen for local potential',nsloc)
        endif
        if (lls(nsloc).ne.lloc) then
           if (rel==2) then
              call errore('ld1_setup','mismatch between lloc and l of ' // &
           &            'spin-orbit split wfc chosen for local potential',nsloc)
           else
              call errore('ld1_setup','mismatch between lloc and l of ' // &
           &            'the wavefunction chosen for local potential',nsloc)
           end if
        end if
     else
        nsloc=-1   
        nbeta=nwfs
     endif
     !
     !     test the occupations: for pseudopotential generation
     !     all-electron and pseudopotential occupations must match
     !
     do n=1,nwfs
        if (.not.new(n)) then
           if (abs(oc(nstoae(n))-ocs(n)) > 1.0d-8 ) call errore &
	   ('ld1_setup','mismatched all-electron/pseudo occupations',n)
        endif
     enddo
  endif
  !
  !     zero the external potential and total energies
  !
  vxt=0.0_dp
  etot0=0.0_dp
  etots0=0.0_dp
  !
  !    initialize file names (used in all_electron)
  !
  file_wavefunctions  = trim(prefix)//'.wfc'
  file_wavefunctionsps= trim(prefix)//'ps.wfc'
  file_logder   = trim(prefix)//'.dlog'
  file_logderps = trim(prefix)//'ps.dlog'

  return
end subroutine ld1_setup

subroutine occ_spinorb(nwf,nwfx,el,nn,ll,jj,oc,isw)
  !
  ! This subroutine splits the states according to their j as needed to
  ! make a spin orbit calculation.
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx)
  character(len=2) :: el(nwfx)

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
           if (nwf0.gt.nwfx) call errore('ld1_setup','too many wfc',1)
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
              oc(n+1)=max(0.0_dp,oc(n)-min(2.0_dp*ll(n),oc(n)))
              oc(n)=min(2.0_dp*ll(n),oc(n))
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

subroutine occ_spinorbps(nwf,nwfx,el,nn,ll,jj,oc,rcut,rcutus,enls,isw)
  !
  ! This subroutine splits the states according to their j as needed to
  ! make a spin orbit calculation.
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx), rcut(nwfx), rcutus(nwfx), enls(nwfx)
  character(len=2) :: el(nwfx)

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
           if (nwf0.gt.nwfx) call errore('ld1_setup','too many wfc',1)
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
           oc(n+1)=max(0.0_dp,oc(n)-min(2.0_dp*ll(n),oc(n)))
           oc(n)=min(2.0_dp*ll(n),oc(n))
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
