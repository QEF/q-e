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
  use kinds, only : dp
  use ld1inc, only : iswitch, lsd, rel,sl3, lmx, nconf, &
                nbeta, lloc, nsloc, etot, etot0, etots0, vxt,   &
                file_wavefunctions, file_logder,  &
                file_wavefunctionsps, file_logderps, &
                prefix, new, core_state, enls, enl,  &
                nwf,                       jj,   el,   isw,   oc, nstoae, &
                nwfs,          lls,        jjs,  els,  isws,  ocs, &
                nwfts,  nnts,  llts, nnts, jjts, elts, iswts, octs, nstoaets, &
                nwftsc, nntsc, lltsc, jjtsc, eltsc, iswtsc, octsc, nstoaec, &
                lpaw
  use funct, only : get_iexch, dft_is_meta, start_exx !, set_dft_from_name
  implicit none

  integer n, n1, nc, m, nwftot, ios
  logical ok, hf, oep
  !
  !     transform dft in a series of codes for the exchange and
  !     correlation routine
  !
!  if (iswitch /= 2 .or. pseudotype == 1) call set_dft_from_name(dft)
  !
  if (dft_is_meta()) call errore('setup','meta-GGA not implemented yet', 1)
  hf  = get_iexch().eq.5
  if (hf)     call errore('setup','HF not implemented yet',1)
  oep = get_iexch().eq.4
  if (oep) call start_exx
  if (oep.and.iswitch>1) &
     call errore('setup','OEP is implemented only for all-electron calc.',1)
  if (oep.and.rel>0) &
     call errore('setup','OEP is implemented only for non-relativistic calc.',1)
  !
  call set_sl3(sl3,lmx)
  !
  !  make the correspondence all-electron pseudopotential
  !
  if (iswitch >= 3) then
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
                 if (eltsc(n,nc).eq.el(n1).and.ABS(jjtsc(n,nc)-jj(n1))<1.d-7) &
                         & nstoaec(n,nc)=n1
              endif
           enddo
           if (nstoaec(n,nc).eq.0) call errore('ld1_setup', &
                'all electron wfc corresponding to pseudo-state ' &
          &     //eltsc(n,nc)//' not found',nc)
        enddo
     enddo
! 
!  set the test configuration for descreening
!
     if (iswitch == 3) then 
        nwfts=nwftsc(1)
        do n=1,nwfts
           nnts(n)=nntsc(n,1)
           llts(n)=lltsc(n,1)
           elts(n)=eltsc(n,1)
           jjts(n) = jjtsc(n,1)
           iswts(n)=iswtsc(n,1)
           octs(n)=octsc(n,1)
           nstoaets(n)=nstoaec(n,1)
        enddo
     endif
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
        if (lpaw .and. ocs(nsloc)>0.0_DP) &
           call errore('ld1_setup','Paw generation with electrons' // &
           &    'in the local channel is not available',1)
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
  enl=0.0_dp
  !
  !    initialize file names (used in all_electron)
  !
  file_wavefunctions  = trim(prefix)//'.wfc'
  file_wavefunctionsps= trim(prefix)//'ps.wfc'
  file_logder   = trim(prefix)//'.dlog'
  file_logderps = trim(prefix)//'ps.dlog'

  return
end subroutine ld1_setup

