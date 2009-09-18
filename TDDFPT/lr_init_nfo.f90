!--------------------------------------------------------------
!OBM This subroutine initialises stuff related to open shell 
! calculations (kpoint > 1 degauss/=0 or nspin/=1)
!-------------------------------------------------------------
#include "f_defs.h"
subroutine lr_init_nfo()
  !
  USE kinds, ONLY : DP
  use klist,                only : nks,degauss,lgauss,ngauss,xk, nelec
  USE wvfct,                ONLY : nbnd, et, igk, npw, g2kin 
  use realus,               only : npw_k, igk_k
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : pi, degspin
  USE noncollin_module,     ONLY : noncolin
  USE mp,                   ONLY : mp_max, mp_min
  USE mp_global,            ONLY : inter_pool_comm
  USE gvect,                ONLY : ngm,g, ecutwfc
  USE cell_base,            ONLY : bg, tpiba, tpiba2, omega
  USE ener,                 ONLY : Ef
  USE ktetra,               ONLY : ltetra
  USE lsda_mod,             ONLY : lsda
  USE realus,               ONLY : real_space
  USE control_ph,            ONLY : alpha_pv, nbnd_occ
 !
  implicit none
  !
  ! local variables
  real(kind=DP) :: small, emin, emax, xmax, fac, targ
  integer       :: ik,ibnd, ipol
  !
  if (.not. real_space) then
  do ik=1,nks
      !
      CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
      !
      npw_k(ik) = npw
      !
      igk_k(:,ik) = igk(:)
      !
     !
  enddo
  endif
  !OBM!! The following part is derived from phonon phq_setup
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  !if (.not. allocated (nbnd_occ) allocate( nbnd_occ (nks) )
  if (lgauss) then
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     if (ngauss.eq. - 99) then
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     endif
     targ = ef + xmax * degauss
     do ik = 1, nks
        do ibnd = 1, nbnd
           if (et (ibnd, ik) .lt.targ) nbnd_occ (ik) = ibnd
        enddo
        if (nbnd_occ (ik) .eq.nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     enddo
  else if (ltetra) then
     call errore('lr_init_nfo','phonon + tetrahedra not implemented', 1)
  else
     if (lsda) call infomsg('lr_init_nfo','occupation numbers probably wrong')
     if (noncolin) then
        nbnd_occ = nint (nelec) 
     else
        do ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / degspin
        enddo
     endif
  endif
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  do ik = 1, nks
     do ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     enddo
  enddo
#ifdef __PARA
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif
  if (lgauss) then
     emax = targ
     alpha_pv = emax - emin
  else
     emax = et (1, 1)
     do ik = 1, nks
        do ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        enddo
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  endif
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
return
end subroutine lr_init_nfo
