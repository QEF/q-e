!
! Copyright (C) 2004-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE run_test
  !
  !   This routine is a driver to the tests of the pseudopotential
  !---------------------------------------------------------------
  !
  use kinds, only : dp
  use io_global, only : ionode, ionode_id
  use mp,        only : mp_bcast
  use mp_world,  only : world_comm
  use radial_grids, only : ndmx
  use ld1_parameters, only : nwfx
  use ld1inc,    only : file_tests, prefix, nconf, rel, etot0, &
                    nbeta, grid, psi, pseudotype, els, zed, &
                    rcut, rcutus, rcutts, rcutusts,  etot, etots0, etots, &
                    nwf,         ll,                oc, el, &
                    nwfts, nnts, llts, jjts, iswts, octs, elts, nstoaets, &
                    nwftsc, nntsc, lltsc, jjtsc, iswtsc, octsc, eltsc,nstoaec, &
                    file_wavefunctions, file_logder, &
                    file_wavefunctionsps, file_logderps, &
                    core_state, use_paw_as_gipaw !EMINE
  implicit none

  integer  &
       n, &  ! counter on wavefunctions
       n1,&  ! counter on mesh points
       ir,&  ! counter on mesh points
       im,&  ! position of the maximum
       nc    ! counter on configurations
  integer   ::      &
       nn_old(nwfx), ll_old(nwfx), nwf_old, isw_old(nwfx), lsd_old
  real(DP) ::              &
       jj_old(nwfx), oc_old(nwfx), enl_old(nwfx), psi_old(ndmx,2,nwfx)
  logical ::  &
       core_state_old(nwfx)
  integer :: ios
  character(len=1) :: nch
  real(DP) :: dum

  file_tests = trim(prefix)//'.test'
  if (ionode) &
     open(unit=13, file=file_tests, iostat=ios, err=1111, status='unknown')
1111 call mp_bcast(ios, ionode_id, world_comm)
     call errore('run_test','opening file_tests',abs(ios))

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
     if (nc>1) call save_ae(nwf_old,nn_old,ll_old,jj_old,enl_old,   &
                            oc_old,isw_old, core_state_old,psi_old,lsd_old,-1)
     !EMINE
     core_state_old = core_state
     call set_conf(nc)
     call all_electron(.true.,nc)
     if (nc.eq.1) then
         etot0=etot
         call save_ae(nwf_old,nn_old,ll_old,jj_old,enl_old,oc_old,isw_old, &
                core_state_old,psi_old,lsd_old,1)
     endif
     !
     !   choose the cut-off radius for the initial estimate of the wavefunctions
     !   find the maximum of the all electron wavefunction
     !
     do n=1,nwfts
        do n1=1,nbeta
           if (els(n1).eq.elts(n).and.(rcut(n1).gt.1.e-3_dp.or.&
                                      rcutus(n1)>1.e-3_DP)) then
              rcutts(n)=rcut(n1)
              rcutusts(n)=rcutus(n1)
              goto 20
           endif
        enddo
        dum=0.0_dp
        im=2
        do ir=1,grid%mesh-1
           dum=abs(psi(ir+1,1,nstoaets(n)))
           if(dum.gt.abs(psi(ir,1,nstoaets(n)))) im=ir+1
        enddo
        if (pseudotype.lt.3) then
           rcutts(n)=grid%r(im)*1.1_dp
           rcutusts(n)=grid%r(im)*1.1_dp
           if (el(nstoaets(n))=='6S'.or.el(nstoaets(n))=='5S') then
              rcutts(n)=grid%r(im)*1.2_dp
              rcutusts(n)=grid%r(im)*1.2_dp
           endif
        else
           if (ll(nstoaets(n)).eq.0) then
              rcutts(n)=grid%r(im)*1.6_dp
              rcutusts(n)=grid%r(im)*1.7_dp
           elseif (ll(nstoaets(n)).eq.1) then
              rcutts(n)=grid%r(im)*1.6_dp
              rcutusts(n)=grid%r(im)*1.7_dp
              if (el(nstoaets(n)).eq.'2P') then
                 rcutts(n)=grid%r(im)*1.7_dp
                 rcutusts(n)=grid%r(im)*1.8_dp
              endif
           elseif (ll(nstoaets(n)).eq.2) then
              rcutts(n)=grid%r(im)*2.0_dp
              rcutusts(n)=grid%r(im)*2.2_dp
              if (el(nstoaets(n)).eq.'3D') then
                 rcutts(n)=grid%r(im)*2.5_dp
                 if (zed>28) then
                    rcutusts(n)=grid%r(im)*3.4_dp
                 else
                    rcutusts(n)=grid%r(im)*3.0_dp
                 endif
              endif
           endif
        endif
20   continue
!     write(6,*) n, rcutts(n), rcutusts(n)
     enddo
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
     call test_bessel ( )
     !
  enddo
  if (ionode) close (unit = 13)  

  !EMINE
  if(use_paw_as_gipaw)core_state = core_state_old

END SUBROUTINE run_test


SUBROUTINE save_ae(nwf_old,nn_old,ll_old,jj_old,enl_old,oc_old,isw_old, &
               core_state_old,psi_old,lsd_old,iflag)
!
! This routine saves the all-electron configuration, or copy it in the
! all-electron variables 

use kinds, only : dp
use ld1_parameters, only : nwfx
use radial_grids, only : ndmx
use ld1inc, only : nwf, nn, ll, jj, enl, oc, isw, core_state, psi, lsd
implicit none
integer, intent(in) :: iflag
integer   ::      &
       nn_old(nwfx),   &   ! the main quantum number
       ll_old(nwfx),   &   ! the orbital angular momentum
       nwf_old,        &   ! the number of wavefunctions
       lsd_old,        &   ! if 1 the calculation has spin
       isw_old(nwfx)       ! spin of the wfc. if(.not.lsd) all 1 (default)

real(DP) ::              &
       jj_old(nwfx),     & ! the total angular momentum
       oc_old(nwfx),     & ! the occupations of the all-electron atom
       enl_old(nwfx),        & ! the energies of the all-electron atom
       psi_old(ndmx,2,nwfx)    ! the all-electron (dirac) wavefunctions

logical ::   &
       core_state_old(nwfx) 

if (iflag==1) then
   nwf_old=nwf
   lsd_old=lsd
   nn_old=nn
   ll_old=ll
   jj_old=jj
   oc_old=oc
   isw_old=isw
   enl_old=enl
   psi_old=psi
else
   nwf=nwf_old
   lsd=lsd_old
   nn=nn_old
   ll=ll_old
   jj=jj_old
   oc=oc_old
   isw=isw_old
   enl=enl_old
   psi=psi_old
endif

END SUBROUTINE save_ae


SUBROUTINE set_conf(nc)
!
!  This routine copy the variables of the current configuration in the
!  test variables. If we pass from a non-spin polarized to a spin
!  polarized calculation it splits the occupations of the all-electron states.
!
use ld1_parameters, only : nwfx
use ld1inc, only : nwf, nn, ll, oc, isw, el, enl, psi, nstoaets, nwftsc,  &
                   core_state, lsdts, eltsc, iswtsc, nnts, llts, jjts,   &
                   octs, elts, iswts, nntsc, lltsc, jjtsc, octsc, &
                   jj, frozen_core, lsd, nwfts, nspin
implicit none
integer, intent(in) :: nc
integer :: n, n1

if (lsdts(nc)==1) then
   if (frozen_core.and.nc>1) then
      call occ_spin_tot(nwf,nwfx,el,nn,ll,oc,isw,enl,psi)
   else
      call occ_spin(nwf,nwfx,el,nn,ll,oc,isw)
   endif 
   lsd=1
   nspin=2
else
   lsd=0
   nspin=1
endif

do n=1,nwf
   core_state(n)=.true.
enddo
do n=1,nwftsc(nc)
   nstoaets(n)=0
   do n1=1,nwf
      if (lsdts(nc).eq.1) then
         if (eltsc(n,nc).eq.el(n1) &
                    .and.iswtsc(n,nc).eq.isw(n1)) then
            nstoaets(n)=n1
            core_state(n1)=.false.
         endif
      else
         if (eltsc(n,nc).eq.el(n1).and.jjtsc(n,nc).eq.jj(n1)) then
            nstoaets(n)=n1
            core_state(n1)=.false.
         endif
      endif
   enddo
   if (nstoaets(n).eq.0) call errore('set_conf', &
                'all electron wfc corresponding to pseudo-state ' &
          &     //eltsc(n,nc)//' not found',nc)
enddo

do n=1,nwfts
   nnts(n)=nntsc(n,nc)
   llts(n)=lltsc(n,nc)
   elts(n)=eltsc(n,nc)
   jjts(n)=jjtsc(n,nc)
   iswts(n)=iswtsc(n,nc)
   octs(n)=octsc(n,nc)
   oc(nstoaets(n))=octs(n)
enddo

END SUBROUTINE set_conf

