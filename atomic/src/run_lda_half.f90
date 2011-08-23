!
! Copyright (C) 2010 Quantm-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine run_lda_half
  !
  !   This routine is a driver to correct pseudopotentials with LDA-1/2
  !   Courtesy of Leonardo Matheus Marion jorge, University of Sao Paolo (Brazil)
  !   L. G. Ferreira, M. Marques and L. K. Teles, Phys. Rev. B 78 125116 (2008)
  !---------------------------------------------------------------
  !
  use kinds, only : dp
  use io_global, only : ionode, ionode_id, stdout
  use mp,        only : mp_bcast
  use radial_grids
  use ld1_parameters, only : nwfx
  use ld1inc,    only : file_tests, prefix, nconf, rel, etot0, &
                    nbeta, grid, psi, pseudotype, els, zed, bmat, &
                    rcut, rcutus, rcutts, rcutusts,  etot, etots0, etots, &
                    nwf, lls, ikk, betas, ll, file_potscf, oc, el, &
                    nwfts, nnts, llts, jjts, iswts, octs, elts, nstoaets, &
                    nwftsc, nntsc, lltsc, jjtsc, iswtsc, octsc, eltsc,nstoaec, &
                    file_wavefunctions, file_logder, file_pseudopw, &
                    file_wavefunctionsps, file_logderps, vpot, vpsloc, rcutv
  implicit none

  integer  &
       n, &  ! counter on wavefunctions
       n1,&  ! counter on mesh points
       ir,&  ! counter on mesh points
       im,&  ! position of the maximum
       nc,&  ! counter on configurations
       nb    ! counter on betas
  integer   ::      &
       nn_old(nwfx), ll_old(nwfx), nwf_old, isw_old(nwfx), lsd_old
  real(DP) ::              &
       jj_old(nwfx), oc_old(nwfx), enl_old(nwfx), psi_old(ndmx,2,nwfx), beta2, f
  logical ::  &
       core_state_old(nwfx)
  integer :: ios, ncut
  character(len=1) :: nch
  real(DP) :: dum, wrcutv
  real(DP) :: dvpot(ndmx,2)
!  file_tests = trim(prefix)//'.test'
!  if (ionode) &
!     open(unit=13, file=file_tests, iostat=ios, err=1111, status='unknown')
!1111 call mp_bcast(ios, ionode_id)
!     call errore('ld1_setup','opening file_tests',abs(ios))

  do nc=1,nconf
     write (nch, '(i1)') nc 
     nwfts=nwftsc(nc)
     call set_conf(nc)
     call all_electron(.true.,nc)
     !
     if (nc.eq.1) then
        dvpot = vpot
     elseif (nc .eq. 2) then 
        dvpot = dvpot - vpot
     endif
  enddo
  ncut = 8
  do ir=1, grid%mesh
     if (grid%r(ir).le.rcutv) then
        wrcutv = (1.0_dp - (grid%r(ir)/rcutv)**ncut)**3
        dvpot(ir,1) = dvpot(ir,1)*wrcutv
        vpsloc(ir) = vpsloc(ir) - dvpot(ir,1)
     endif
  enddo
!
! re-call set conf to write in the pseudo file the Valence
! config without LDA-1/2
!
  call set_conf(1)
  return
end subroutine run_lda_half
