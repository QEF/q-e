!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rho_core
  !-----------------------------------------------------------------------
  !
  !      input : all-electron wavefunctions + valence states
  !      output: smoothed core charge for r > rcore
  !
  use kinds, only : dp
  use constants, only : pi
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use mp_world,  only : world_comm
  use ld1inc, only : nlcc, grid, rhoc, aeccharge, psccharge, rcore, &
                     nwf, oc, rel, core_state, psi, file_core, new_core_ps,&
                     lpaw, lnc2paw
  implicit none

  real(DP) :: drho, const, br1, br2, &
       eps1, eps2, br12, xc(8), a, b, eps12, totrho
  real(DP), allocatable:: rhov(:)
  real(DP), external :: int_0_inf_dr
  integer :: i, ik, n, ns, ios

  if (nlcc) then
     write(stdout,'(/,5x,'' Computing core charge for nlcc: '')')
  else
     if (lpaw) write(stdout,'(/,5x,'' Computing core charge for PAW: '')')
  end if
  allocate (rhov(grid%mesh))
  !
  !      calculates core charge density
  !
  do n=1,grid%mesh
     rhov(n) = 0.0_dp
     rhoc(n) = 0.0_dp
     do ns=1,nwf
        if (rel==2) then
           if (core_state(ns)) then
              rhoc(n)=rhoc(n)+oc(ns)*(psi(n,1,ns)**2+psi(n,2,ns)**2)
           else
              rhov(n)=rhov(n)+oc(ns)*(psi(n,1,ns)**2+psi(n,2,ns)**2)
           endif
        else
           if (core_state(ns)) then
              rhoc(n) = rhoc(n) + oc(ns)*psi(n,1,ns)**2
           else
              rhov(n) = rhov(n) + oc(ns)*psi(n,1,ns)**2
           endif
        endif
     enddo
  enddo
  totrho = int_0_inf_dr(rhoc,grid,grid%mesh,2)
  if (totrho<1.d-6.and.lpaw) then
!
!  All valence charge for this atom (mainly for H)
!
     aeccharge(1:grid%mesh) = 0.0_DP
     psccharge(1:grid%mesh) = 0.0_DP
     goto 1100
  endif

!  write(stdout,'("Integrated core charge",f15.10)') totrho
  aeccharge(1:grid%mesh) = rhoc(1:grid%mesh)
  !
  if (rcore > 0.0_dp) then
     !      rcore read on input
     do ik=1,grid%mesh
        if (grid%r(ik) > rcore) go to 100
     enddo
     call infomsg('set_rho_core','rcore too big')
     return
  else
     !      rcore determined by the condition  rhoc(rcore) = 2*rhov(rcore)
     do ik=1,grid%mesh
        if (rhoc(ik) < 2.0 * rhov(ik)) go to 100
     enddo
  end if
100 rcore=grid%r(ik)
  drho = ( rhoc(ik+1)/grid%r2(ik+1) - rhoc(ik)/grid%r2(ik) ) / grid%dx / grid%r(ik)
  !
  !   true_rho = rhoc(r)/r**2/4 pi
  !      (factor 1/r from logarithmic mesh)
  !   smoothened core charge density for nonlinear core correction:
  !      rhoc(r) = core charge        for r > rcore
  !      rhoc(r) = r^2 a sin(b r)/r   for r < rcore
  !
  if (new_core_ps) then
     call compute_phius(1,ik,aeccharge,rhoc,xc,0,'  ')
  else
     if (drho > 0.0_dp) then
        call infomsg('set_rho_core','d rho/ d r > 0')
        return
     endif
     const= grid%r(ik)*drho / ( rhoc(ik)/grid%r2(ik) ) + 1.0_dp
     if (const > 0.0_dp) then
        br1 = 0.00001_dp
        br2 = pi/2.0_dp-0.00001_dp
     else
        br1 = pi/2.0_dp+0.00001_dp
        br2 = pi
     end if
     do n=1, 15
        eps1 = br1 - const*tan(br1)
        eps2 = br2 - const*tan(br2)
        br12 = (br1+br2)/2.0_dp
        eps12 = br12 - const*tan(br12)
        if(eps1*eps12 < 0.0_dp) then
           br2 = br12
        else if(eps12*eps2 < 0.0_dp) then
           br1 = br12
        else
           call errore('set_rho_core','error in bisection',n)
        end if
     end do
     b = br12/grid%r(ik)
     a = ( rhoc(ik)/grid%r2(ik) ) * grid%r(ik)/sin(br12)
     do n=1,ik
        rhoc(n) = a*sin(b*grid%r(n))/grid%r(n) * grid%r2(n)
     end do
  endif
  if (lpaw) then
     if (lnc2paw) then
        ! Mimic NC calculation. If NLCC, the pseudized core charge.
        if (nlcc) then
           aeccharge(1:grid%mesh) = rhoc(1:grid%mesh)
           ! Here one could set another pseudized ccharge, for
           ! example with a larger matching radius. Right now,
           ! just take the same as the AE (ie NC) one:
           psccharge(1:grid%mesh) = rhoc(1:grid%mesh)
        else
           ! Reference NC calculation does not have core charge.
           aeccharge(1:grid%mesh) = 0._dp
           psccharge(1:grid%mesh) = 0._dp
        end if
     else
!        aeccharge(1:grid%mesh) = rhoc(1:grid%mesh)
        if (nlcc) then
           psccharge(1:grid%mesh) = rhoc(1:grid%mesh)
        else
           psccharge(1:grid%mesh) = 0._dp
        end if
     end if
  end if
  write(stdout,'(/,5x,''  r > '',f4.2,'' : true rho core'')') grid%r(ik)
  if (new_core_ps) then
     write(stdout,'(6x,"Core charge pseudized with two Bessel functions")')
  else
     write(stdout,110) grid%r(ik), a, b
  endif
110 format (5x, '  r < ',f4.2,' : rho core = a sin(br)/r', &
       '    a=',f7.2,'  b=',f7.2/)
1100 continue
  if (file_core .ne. ' ') then
     write(stdout,'(6x, "***Writing file ",a, " ***")') trim(file_core)
     if (ionode) &
        open(unit=26,file=file_core, status='unknown', iostat=ios, err=300 )
300  call mp_bcast(ios, ionode_id, world_comm)
     call errore('set_rho_core','opening file '//file_core,abs(ios))
     if (ionode) then
        if (totrho>1.d-6) then
           do n=1,grid%mesh
              write(26,'(4f20.10)') grid%r(n),rhoc(n),rhov(n),aeccharge(n)
           enddo
        else
           do n=1,grid%mesh
              write(26,'(2f20.10)') grid%r(n),rhov(n)
           enddo
        endif
        close(26)
     endif
  endif
  totrho = int_0_inf_dr(rhoc,grid,grid%mesh,2)
  write(stdout,'(6x,''Integrated core pseudo-charge : '',f6.2)')  totrho
  if (.not.nlcc) rhoc(1:grid%mesh) = 0.0_dp
  deallocate (rhov)
  return
end subroutine set_rho_core
