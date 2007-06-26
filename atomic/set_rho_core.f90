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
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use ld1inc
  implicit none

  real(DP) :: drho, const, br1, br2, pi, &
       eps1, eps2, br12, a, b, eps12, totrho
  real(DP), allocatable:: rhov(:), rhoco(:)
  real(DP), external :: int_0_inf_dr
  integer :: i, ik, n, ns, ios

  if (nlcc) then
     write(stdout,'(/,5x,'' Computing core charge for nlcc: '')')
  else
     if (lpaw) write(stdout,'(/,5x,'' Computing core charge for PAW: '')')
  end if
  pi = 4.0_dp*atan(1.0_dp)
  allocate (rhov(mesh), rhoco(mesh))
  !
  !      calculates core charge density
  !
  do n=1,mesh
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
!  totrho = int_0_inf_dr(rhoc,r,r2,dx,mesh,2)
!  write(stdout,'("Integrated core charge",f15.10)') totrho
  rhoco(:) = rhoc(1:mesh)
  if (lpaw) aeccharge(1:mesh) = rhoc(1:mesh)
  !
  if (rcore > 0.0_dp) then
     !      rcore read on input
     do ik=1,mesh
        if (r(ik) > rcore) go to 100
     enddo
     call infomsg('set_rho_core','rcore too big', -1)
     return
  else
     !      rcore determined by the condition  rhoc(rcore) = 2*rhov(rcore)
     do ik=1,mesh
        if (rhoc(ik) < 2.0 * rhov(ik)) go to 100
     enddo
  end if
100 rcore=r(ik)
  drho = ( rhoc(ik+1)/r2(ik+1) - rhoc(ik)/r2(ik) ) / dx / r(ik)
  !
  !   true_rho = rhoc(r)/r**2/4 pi
  !      (factor 1/r from logarithmic mesh)
  !   smoothened core charge density for nonlinear core correction:
  !      rhoc(r) = core charge        for r > rcore
  !      rhoc(r) = r^2 a sin(b r)/r   for r < rcore
  !
  if (drho > 0.0_dp) then
     call infomsg('set_rho_core','d rho/ d r > 0', -1)
     return
  endif
  const= r(ik)*drho / ( rhoc(ik)/r2(ik) ) + 1.0_dp
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
  b = br12/r(ik)
  a = ( rhoc(ik)/r2(ik) ) * r(ik)/sin(br12)
  do n=1,ik
     rhoc(n) = a*sin(b*r(n))/r(n) * r2(n)
  end do
  if (lpaw) psccharge(1:mesh) = rhoc(1:mesh)
  write(stdout,'(/,5x,''  r > '',f4.2,'' : true rho core'')') r(ik)
  write(stdout,110) r(ik), a, b
110 format (5x, '  r < ',f4.2,' : rho core = a sin(br)/r', &
       '    a=',f7.2,'  b=',f7.2/)
  if (file_core .ne. ' ') then
     write(stdout,*) '***Writing file ',trim(file_core),' ***'
     if (ionode) &
        open(unit=26,file=file_core, status='unknown', iostat=ios, err=300 )
300  call mp_bcast(ios, ionode_id)
     call errore('set_rho_core','opening file '//file_core,abs(ios))
     if (ionode) then
        do n=1,mesh
           write(26,'(4f20.10)') r(n),rhoc(n),rhov(n),rhoco(n)
        enddo
        close(26)
     endif
  endif
  deallocate (rhoco, rhov)
  totrho = int_0_inf_dr(rhoc,r,r2,dx,mesh,2)
  write(stdout,'(13x,''integrated core pseudo-charge : '',f6.2)')  totrho
  if (.not.nlcc) rhoc(1:mesh) = 0.0_dp
  return
end subroutine set_rho_core
