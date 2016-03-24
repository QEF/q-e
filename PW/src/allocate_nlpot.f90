!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_nlpot
  !-----------------------------------------------------------------------
  !
  ! This routine computes the dimension of the Hamiltonian matrix and
  ! allocates arrays containing the non-local part of the pseudopotential
  !
  ! It computes the following global quantities:
  !
  !     ngk           !  number of plane waves (for each k point)
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE cellmd,           ONLY : cell_factor
  USE gvect,            ONLY : ngm, gcutm, g
  USE klist,            ONLY : xk, wk, ngk, nks, qnorm, igk_k
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : Hubbard_lmax
  USE scf,              ONLY : rho
  USE noncollin_module, ONLY : noncolin
  USE wvfct,            ONLY : npwx, npw, igk, g2kin
  USE gvecw,            ONLY : gcutw, ecutwfc
  USE us,               ONLY : qrad, tab, tab_d2y, tab_at, dq, nqx, &
                               nqxq, spline_ps
  USE uspp,             ONLY : indv, nhtol, nhtolm, ijtoh, qq, dvan, deeq, &
                               vkb, indv_ijkb0, okvan, nkb, nkbus, nhtoj, &
                               becsum, qq_so,dvan_so, deeq_nc
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nbetam
  USE spin_orb,         ONLY : lspinorb, fcoef
  USE control_flags,    ONLY : program_name
  USE io_global,        ONLY : stdout
  !
  implicit none
  !
  integer :: nwfcm
  !
  !   calculate number of PWs for all kpoints
  !
  allocate (ngk( nks ))
  !
  call n_plane_waves (gcutw, nks, xk, g, ngm, npwx, ngk)
  !
  !   igk relates the index of PW k+G to index in the list of G vector
  !
  allocate ( igk_k( npwx,nks ) )
  allocate (igk( npwx ), g2kin ( npwx ) )    
  !
  ! Note: computation of the number of beta functions for
  ! each atomic type and the maximum number of beta functions
  ! and the number of beta functions of the solid has been
  ! moved to init_run.f90 : pre_init()
  !
  allocate (indv( nhm, nsp))    
  allocate (nhtol(nhm, nsp))    
  allocate (nhtolm(nhm, nsp))    
  allocate (nhtoj(nhm, nsp))    
  allocate (ijtoh(nhm, nhm, nsp))
  allocate (indv_ijkb0(nat))
  allocate (deeq( nhm, nhm, nat, nspin))    
  if (noncolin) then
     allocate (deeq_nc( nhm, nhm, nat, nspin))    
  endif
  allocate (qq(   nhm, nhm, nsp))    
  if (lspinorb) then
    allocate (qq_so(nhm, nhm, 4, nsp))    
    allocate (dvan_so( nhm, nhm, nspin, nsp))    
    allocate (fcoef(nhm,nhm,2,2,nsp))
  else
    allocate (dvan( nhm, nhm, nsp))    
  endif
  ! GIPAW needs a slighly larger q-space interpolation for quantities calculated
  ! at k+q_gipaw
  if (trim(program_name) == 'GIPAW') then
    if (cell_factor == 1.d0) cell_factor = 1.1d0
    write(stdout,"(5X,'q-space interpolation up to ',F8.2,' Rydberg')") ecutwfc*cell_factor
  endif
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = INT( ( (sqrt(gcutm) + qnorm) / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  !
  if (lmaxq > 0) allocate (qrad( nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp))    
  allocate (vkb( npwx,  nkb))    
  allocate (becsum( nhm * (nhm + 1)/2, nat, nspin))    
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = INT( (sqrt (ecutwfc) / dq + 4) * cell_factor )

  allocate (tab( nqx , nbetam , nsp))

  ! d2y is for the cubic splines
  if (spline_ps) allocate (tab_d2y( nqx , nbetam , nsp))

  nwfcm = MAXVAL ( upf(1:nsp)%nwfc )
  allocate (tab_at( nqx , nwfcm , nsp))

  return
end subroutine allocate_nlpot

