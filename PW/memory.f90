!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program pwmemory
  !-----------------------------------------------------------------------
  !
  USE wvfct,    ONLY : gamma_only
  USE klist,    ONLY : nkstot, xk
  USE lsda_mod, ONLY : nspin
  USE io_files
  USE mp,       ONLY : mp_end
  USE global_version
  !
  implicit none
  !
  character(len=9) :: code = 'memory'
  !
  call startup (nd_nmbr, code, version_number)
  !
  call iosys
  call setup
  !
  gamma_only = ( (nkstot == 1) .and. (nspin == 1) .and. &
                 ( xk(1,1)**2 + xk(2,1)**2 + xk(3,1)**2 < 1.0e-6 ) ) &
          .OR. ( (nkstot == 2) .and. (nspin == 2) .and. &
                 ( xk(1,1)**2 + xk(2,1)**2 + xk(3,1)**2 < 1.0e-6 ) .AND. &
                 ( xk(1,2)**2 + xk(2,2)**2 + xk(3,2)**2 < 1.0e-6 ) )
  !
  call data_structure(gamma_only)
  !
  call setup2()
  !
  call memory_estimate()
  !
  call mp_end ()
  stop
end program pwmemory
!
!-----------------------------------------------------------------------
subroutine setup2()
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega
  USE cellmd,     ONLY : cell_factor
  USE constants,  ONLY :tpi, fpi
  USE gvect,      ONLY : ecutwfc, gcutm, ngl, ngm
  USE klist,      ONLY : xqq
  USE uspp_param, ONLY : lmaxkb, nbeta, lll, lmaxq, nhm, tvanp, nh
  USE us,         ONLY : nqx, nqxq, dq
  USE uspp,       ONLY : nkb
  USE wvfct,      ONLY : npwx
  !
  implicit none
  !
  real(DP) :: omegaBZ
  integer :: nb, nt
  !
  !     calculate the number of beta functions for each atomic type
  !
  lmaxkb = - 1
  do nt = 1, ntyp
     nh (nt) = 0
     do nb = 1, nbeta (nt)
        nh (nt) = nh (nt) + 2 * lll (nb, nt) + 1
        lmaxkb = max (lmaxkb, lll (nb, nt) )
     enddo
  enddo
  lmaxq = 2*lmaxkb+1
  !
  ! maximum number of projectors (beta) per atom
  !
  nhm = MAXVAL (nh(1:ntyp))
  !
  ! total number of projectors (beta)
  !
  nkb = SUM (nh(ityp(1:nat)))
  !
  ! number of points for interpolation tables: qrad and table
  ! (including a possible factor coming from cell contraction
  !  during variable cell relaxation/MD)
  !
  nqxq = ( (sqrt(gcutm) + sqrt(xqq(1)**2 + xqq(2)**2 + xqq(3)**2) ) &
          / dq + 4) * cell_factor
  nqx = (sqrt (ecutwfc) / dq + 4) * cell_factor
  !
  ! estimate npwx as (volume in g-space/volume of the BZ)
  !
  omegaBZ = tpi**3/omega
  npwx = fpi /3.d0 * ecutwfc**(1.5d0) / omegaBZ
  !
  ! ngl is assumed to have its maximum value (ngm)
  !
  ngl = ngm
  !
  return
end subroutine setup2

!
!-----------------------------------------------------------------------
subroutine memory_estimate ( )
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE parameters,    ONLY : ndmx, npsx, lmaxx, nchix, nbrx, nqfx, lqmax
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE basis,         ONLY : natomwfc
  USE cellmd,        ONLY : lmovecell
  USE klist,         ONLY : npk, nks, nkstot
  USE ktetra,        ONLY : ntetra
  USE gvect,         ONLY : nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                            nrxx, ngl, ngm, ngm_l
  USE gsmooth,       ONLY : ngms, doublegrid
  USE ldaU,          ONLY : Hubbard_lmax, lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE us,            ONLY : nqx, nqxq
  USE uspp,          ONLY : nkb, okvan
  USE uspp_param,    ONLY : lmaxq, nhm
  USE control_flags, ONLY : nmix, isolve, diis_ndim
  USE wvfct,         ONLY : gamma_only, npwx, nbnd, nbndx
  USE mp_global,     ONLY : nproc_pool, npool
  USE pfft,          ONLY : nct, ncplane
  USE pffts,         ONLY : ncts, ncplanes
  !
  implicit none
  !
  ! The following variables are real in order to prevent
  ! integer overflow for very large jobs
  !
  real(DP), parameter :: real_size = 8, int_size = 8
  real(DP), parameter :: comp_size = 2*real_size
  real(DP) :: total_mem, scalable_mem, nonscalable_mem
  real(DP) :: scalable_wspace, nonscalable_wspace
  real(DP) :: wspace_diago, wspace_mix, diis_steps
  !
  ! fixed memory, or memory allocated in iosys, setup
  !
  nonscalable_mem = &
       real_size * 2 * 3 * nat +  int_size * nat +   & ! tau force ityp
       int_size * 48 * nat + int_size * 4 * ntetra + & ! irt tetra
       real_size * 4 * npk + int_size * 2 * npk +    & ! xk wk ngk isk
       real_size * ndmx * npsx * (4+(lmaxx+1)+nchix)+ & ! atomic PP
       real_size * ndmx * npsx * (nbrx + nbrx*nbrx) + & ! atomic USPP
       real_size * nqfx * lqmax * nbrx * nbrx * npsx       ! qfcoef
  !
  ! dynamically allocated memory that does not scale with N CPUs
  !
  nonscalable_mem = nonscalable_mem + int_size * ngm_l + & ! ig_l2g
       real_size * nqx * nbrx * ntyp +                   & ! tab
       real_size * nhm * (nhm + 1)/2 * nat * nspin +     & ! becsum
       real_size * nqxq* nbrx * (nbrx+1)/2 * lmaxq * ntyp +& ! qrad
       comp_size * ( 2*nr1 + 2*nr2 + 2*nr3 + 3) * nat +  & ! eigts
       real_size * ( nhm * nhm * ( nat*nspin + 2*ntyp) )+& ! qq dvan deeq
       real_size * 2 * nbnd * nkstot  +                  & ! et wg
       comp_size * nkb * nbndx                             ! becp
  if (lda_plus_u) & 
     nonscalable_mem = nonscalable_mem + &  ! ns, nsnew
          real_size * 2 * nat * nspin * (2 * Hubbard_lmax + 1)**2 
  !
  ! dynamically allocated memory that scales with N CPUs
  !
  scalable_mem = real_size * 4 * ngm + int_size * 5 * ngm 
  ! g, gg, nl, igtongl, ig1, ig2, ig3
  scalable_mem =  scalable_mem + &
       real_size * nrxx * (5 * nspin + 3) + & ! rho and v in real space
       comp_size * (ngl + ngm) * ntyp +     & ! vloc strf
       int_size  * npwx * (nks + 1 ) +      & ! igk igk_l2g
       real_size * npwx +                   & ! g2kin 
       comp_size * ngm  +                   & ! qgm
       comp_size * npwx * ( nkb + nbnd)       ! vkb, evc
  if (lda_plus_u) &
       scalable_mem = scalable_mem + comp_size * npwx * natomwfc ! swfcatom
  if (doublegrid) scalable_mem = scalable_mem + int_size * ngms ! nls
  if (lmovecell)  scalable_mem = scalable_mem + real_size * ngl !  gl

#ifdef __PARA
  nonscalable_mem=nonscalable_mem + int_size * (ncplane + ncplanes) ! ipc, ipcs
  scalable_mem = scalable_mem + int_size * (nct + ncts ) ! icpl. icpls
#endif
  !
  ! workspace : diagonalization
  !
  if (isolve == 0) then
     if (okvan) then 
        wspace_diago = comp_size * 3 * npwx * nbndx  ! cegter: psi hpsi spsi
     else
        wspace_diago = comp_size * 2 * npwx * nbndx  ! cegter: psi hpsi
     end if
  else if (isolve == 1 ) then
     wspace_diago = comp_size * 7 * npwx
  else if (isolve == 2) then
     wspace_diago = comp_size * npwx * 4*diis_ndim
  end if
  nonscalable_wspace=comp_size * 3 * nbndx * nbndx  ! hc sc vc
#ifdef __PARA
  nonscalable_wspace = nonscalable_wspace + real_size * nrx1*nrx2*nrx3
  ! psymrho, io_pot
#endif
  !
  ! workspace : mixing (mix_rho, save on file)
  !
  wspace_mix = comp_size * 2 * ngm * nspin * (nmix + 1)
  !
  scalable_wspace = max (wspace_mix, wspace_diago)
  !
  total_mem = scalable_mem + nonscalable_mem + &
              scalable_wspace + nonscalable_wspace
#ifdef __PARA
  print '(5x,"Number of processors/pools:",2i4)', nproc_pool, npool
#endif
  if (gamma_only) then
     print '(5x,"Estimated Max memory (Gamma-only code): ",f8.2,"Mb")', &
          total_mem/1024/1024
  else
     print '(5x,"Estimated Max memory (k-point code): ",f8.2,"Mb")', &
          total_mem/1024/1024
  end if
  !
  print '(5x,"nonscalable memory =",f8.2,"Mb")', &
       nonscalable_mem/1024/1024
  print '(5x,"   scalable memory =",f8.2,"Mb")', &
       scalable_mem/1024/1024
  print '(5x,"nonscalable wspace =",f8.2,"Mb")', &
       nonscalable_wspace/1024/1024
  print '(5x,"   scalable wspace =",f8.2,"Mb",   &
      & "   (diag:",f8.2,"Mb, mix:",f8.2,"Mb)")', &
       scalable_wspace/1024/1024, &
       wspace_diago/1024/1024,    &
       wspace_mix/1024/1024
  print '(5x,"Estimates based on ",i2,"-byte complex, ",i2,"-byte reals, ", &
       & i2,"-byte integers")', nint(comp_size), nint(real_size), nint(int_size)
  !
  return
end subroutine memory_estimate
