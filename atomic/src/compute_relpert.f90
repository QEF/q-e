! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------
subroutine compute_relpert( evel, edar, eso )
  !------------------------------------------
  !
  ! Starting from the all-electron wave-functions and eigenvalues
  ! (calculated with rel=0 and lsd=0), computes the relativistic
  ! corrections to the eigenvalues in perturbation theory (in v/c).
  !
  ! For each wavefunction iwf:
  !   evel(iwf): velocity (p^4) correction
  !   edar(iwf): Darwin term
  !   eso(iwf):  spin-orbit term 
  !
  ! Formulas and notation follow the Herman-Skillman book:
  !
  !  Herman, Frank and Skillman, Sherwood, "Atomic Structure Calculations", 
  !  Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1963
  ! 

  use kinds, only : DP
  use radial_grids, only : series
  use ld1inc, only : rel, isic, zed, cau_fact, &
                     grid, vpot, vsic, psi, &
                     nwf, nn, ll, enl 
  implicit none

  ! output variables
  real(DP), intent(out) :: evel(nwf), edar(nwf), eso(nwf)
  ! local variables
  integer :: nrx
  integer :: iwf, l
  real(DP) :: pre, relct1, relct2, b(0:3) 
  real(DP) :: gvpot(grid%mesh), vtmp(grid%mesh), gtmp(grid%mesh)

  real(DP), external :: int_0_inf_dr

  nrx = grid%mesh
  pre =  0.25_DP / cau_fact**2

  if (isic.eq.0) then
    ! calculate d/dr[ V(r) ], where V = V_ion + V_H + V_xc + V_ext
    call grad_log(vpot(1:nrx,1),gvpot,grid%rm1(1:nrx),grid%dx,nrx,4)
  endif

  !!! loop on the wave functions
  do iwf=1,nwf

    l = ll(iwf) ! angular quantum number
    
    ! if SIC allowed, every wfc has a different potential
    if (isic /= 0) then
      ! calculate d/dr[ V(r) ], where V = V_ion + V_H + V_xc + V_ext
      call grad_log(vpot(1:nrx,1)-vsic(1:nrx,iwf),gvpot,grid%rm1(1:nrx),grid%dx,nrx,2)
    endif

    ! if calculation has no relativistic terms (rel==0)
    if (rel.eq.0) then
      !
      ! calculate velocity (p^4) correction
      if (isic.eq.0) then
        vtmp = (enl(iwf)-vpot(1:nrx,1))**2 * psi(1:nrx,1,iwf)**2
      else 
        vtmp = (enl(iwf)-vpot(1:nrx,1)-vsic(1:nrx,iwf))**2 * psi(1:nrx,1,iwf)**2
      endif
      evel(iwf) = pre * int_0_inf_dr(vtmp,grid,nrx,2*l)
      !
      ! compute Darwin term:
      ! A) "hydrogrenoid" term (vanishes for ang. mom. GT 0)
      if (l.eq.0) then
        call series(psi(1:4,1,iwf)*grid%rm1(1:4),grid%r(1:4),grid%r2(1:4),b)
        ! lim_{r\to 0}[ r^2 * dV/dr ] = 2*zed
        relct1 = -zed * pre * b(0)*b(0)
        !!!
        ! exact value for hydrogenoid atom (one electron and isic=1)
        !relct2 = -4._DP * pre * (zed)**4 / nn(iwf)**3
        !!!
      else
        relct1 = 0._DP
        !!!
        !relct2 = 0._DP
        !!!
      endif
      ! B)
      ! calculate d/dr[ r^2 * dV/dr ]
      vtmp = gvpot(1:nrx)*grid%r2(1:nrx)
      ! third order polynomial extrapolation for the first point
      call series(vtmp(2:5),grid%r(2:5),grid%r2(2:5),b)
      vtmp(1) = b(0) + b(1)*grid%r(1) + b(2)*grid%r2(1) +&
                b(3)*grid%r(1)*grid%r2(1)
      call radial_gradient(vtmp,gtmp,grid%r(1:nrx),nrx,1)
      vtmp = gtmp(1:nrx) * (psi(1:nrx,1,iwf)*grid%rm1(1:nrx))**2
      ! in the hydrogenoid atom (one electron and isic=1), 
      !  this term should be vanishing (xmin-> -\inf)
      edar(iwf) = -0.5_DP * pre * int_0_inf_dr(vtmp,grid,nrx,2*l+1)
      edar(iwf) = edar(iwf) + relct1
      !
    endif

    ! compute spin-orbit term (only for ang. mom. greater than zero)
    if (l.eq.0) then
      eso(iwf) = 0._DP
    else
      if (rel.eq.0) then
        vtmp = grid%rm1(1:nrx)*gvpot(1:nrx)*psi(1:nrx,1,iwf)**2
        eso(iwf) = -1._DP * pre * int_0_inf_dr(vtmp,grid,nrx,2*l-1)
      else
         call errore('compute_relpert','not programmed for rel>0!!!',1)
      endif
    endif

  enddo
  !!! end loop on the wave functions

end subroutine compute_relpert
