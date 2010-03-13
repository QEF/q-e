!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine pseudovloc
  !--------------------------------------------------------------------------
  !
  !     This routine generate a local pseudopotential 
  !     The output of the routine are:
  !     vpsloc: the local pseudopotential
  !      
  use kinds, only : DP
  use radial_grids, only : ndmx
  use io_global, only : stdout
  use ld1inc, only : lloc, rcloc, grid, vpot, vpsloc, rel, nsloc, &
                     phis, els, chis, psipsus, &
                     jjs, nstoae, enls, new, psi, enl, rcut, psipaw, &
                     psipaw_rel
  implicit none

  integer :: &
       nwf0, &  ! used to specify the all electron function
       nst,  &  ! auxiliary
       ik       ! the point corresponding to rc

  real(DP) ::             &
       xc(8),              &  ! the coefficients of the fit
       vaux(ndmx,2),        &  ! keeps the potential
       psi_in(ndmx)            ! auxiliary

  integer ::         &
       n,        &  ! counter on mesh points
       ns,       &  ! auxiliary
       indi,rep, &  ! auxiliary
       indns(0:1)    ! auxiliary

  if (lloc < 0) then
     !
     !   Compute the potential by smoothing the AE potential
     !
     !   Compute the ik which correspond to this cutoff radius
     !
     write(stdout, &
          "(/,5x,' Generating local potential from pseudized AE potential:',&
            &  /,5x,' Matching radius rcloc = ',f8.4)") rcloc
     ik=0
     do n=1,grid%mesh
        if (grid%r(n) < rcloc) ik=n
     enddo
     if (mod(ik,2) == 0) ik=ik+1
     if (ik <= 1 .or. ik > grid%mesh) &
          call errore('pseudovloc','wrong matching point',1)
!
!  smooth the potential before ik.
!
!  ... with the original recipe
     if (lloc==-1) call compute_potps(ik,vpot,vpsloc,xc)
!  ... or with a modified recipe that enforce V''(0)=0 as suggested by TM
     if (lloc==-2) write(stdout,"(5x,' Enforcing V''''(0)=0 (lloc=-2)')")
     if (lloc==-2) call compute_potps_new(ik,vpot,vpsloc,xc)
     write(stdout, 110) grid%r(ik),xc(5)**2 
110  format (/5x, ' Local pseudo, rcloc=',f6.3, &
          ' Estimated cut-off energy= ', f8.2,' Ry')
  else
     !
     !    if a given angular momentum gives the local component this is done 
     !    here
     !
     nst=(lloc+1)*2
     if (rel==2 .and. lloc > 0) then
        rep=1
        indns(0)=nsloc
        indns(1)=nsloc+1
        if (jjs(nsloc) > jjs(nsloc+1) ) then
           indns(0)=nsloc+1
           indns(1)=nsloc
        endif
     else
        rep=0
        indns(0)=nsloc
     endif
     vpsloc=0.0_dp
     vaux=0.0_dp
     do indi=0,rep
        nwf0=nstoae(nsloc+indi)
        if (enls(nsloc+indi) == 0.0_dp)  enls(nsloc+indi)=enl(nwf0)
        !
        !    compute the ik closer to r_cut
        !
        ik=0
        do n=1,grid%mesh
           if (grid%r(n) < rcut(nsloc+indi)) ik=n
        enddo
        if (mod(ik,2).eq.0) ik=ik+1
        if (ik <= 1 .or. ik > grid%mesh) &
           call errore('pseudovloc','wrong matching point',1)
        rcloc=rcut(nsloc+indi)
        if (rep == 0) then
           write(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
                  & ', matching radius rcloc = ',f8.4)") lloc, rcloc
        else
           if (rel==2) then
              write(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
                  &', j=',f5.2,', matching radius rcloc = ',f8.4)") &
                  lloc, lloc-0.5d0+indi, rcloc
           else
              write(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
                  &', spin=',i1,', matching radius rcloc = ',f8.4)") &
                  lloc, indi+1, rcloc
           endif
        endif
        !
        !   compute the phi functions
        !
        ns=indns(indi)
        if (new(ns)) then
           call set_psi_in(ik,lloc,jjs(ns),enls(ns),psi_in,psipaw_rel)
        else
           psi_in(:)=psi(:,1,nwf0)
        endif
        psipaw(:,ns)=psi_in(:)
        !
        !  compute the phi and chi functions
        !
        call compute_phi_tm(lloc,ik,psi_in,phis(1,ns),0,xc,enls(ns),els(ns))
        call compute_chi_tm(lloc,ik,ik+10,phis(1,ns),chis(1,ns),xc,enls(ns))
        !
        !     set the local potential equal to the all-electron one at large r
        !
        do n=1,grid%mesh
           if (grid%r(n) > rcloc) then
              vaux(n,indi+1)=vpot(n,1)
           else
              vaux(n,indi+1)=chis(n,ns)/phis(n,ns)
           endif
        enddo
        psipsus(:,ns)=phis(:,ns)
     enddo
     if (rep==0) then
        do n=1,grid%mesh
           vpsloc(n)=vaux(n,1)
        enddo
     else
        do n=1,grid%mesh
           vpsloc(n)=(lloc*vaux(n,1)+(lloc+1.0_dp)*vaux(n,2))/ &
                (2.0_dp*lloc+1.0_dp)
        enddo
     endif
  endif

  return
end subroutine pseudovloc
