
subroutine init_us_1all

   implicit none
   call init_us_1()
   call init_us_1a()

end subroutine

! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_1a
   !----------------------------------------------------------------------
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: fpi, sqrt2
   USE atom, ONLY: rgrid
   USE ions_base, ONLY: ntyp => nsp, zv
   USE cell_base, ONLY: omega, tpiba
   USE gvect, ONLY: g, gg
   USE lsda_mod, ONLY: nspin
   USE us, ONLY: nqxq, dq, nqx, tab, tab_d2y, qrad, spline_ps
   USE splinelib
   USE uspp_param, ONLY: upf, lmaxq, nbetam, nh, nhm, lmaxkb
   USE mp_global, ONLY: intra_bgrp_comm
   USE mp, ONLY: mp_sum
   use zero_mod, ONLY: l_tab, tabr, tabr_d2y, &
                       tablocal_hg, tablocal_d2y_hg

   USE ions_base, ONLY: nsp
   !
   implicit none
   !
   !     here a few local variables
   !
   integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, iq, is, startq, &
              lastq, ilast, ndm, rm
   ! various counters
   real(DP), allocatable :: aux(:), besr(:)
   ! various work space
   real(DP) ::  pref, q, qi
   ! the prefactor of the q functions
   ! the prefactor of the beta functions
   ! the modulus of g for each shell
   ! q-point grid for interpolation
   real(DP) ::  vqint
   ! interpolated value
   integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
              lk, mk, vk, kh, lh
   !
   real(DP), allocatable :: xdata(:)
   real(DP) :: d1
   !
   call start_clock('init_us_1a')
   !    Initialization of the variables
   !

!!!!!!!aggiunta di codice per la zero_current
!tabr(nqxq,nbetam,nsp,-1:1) e tablocal(nqxq,nsp,0:2)
!contiene f_(la,lb)(q)=\int _0 ^\infty dr r^3 f_la(r) j_lb(q.r)
!dove la è il momento angolare della beta function (ce ne possono
!essere più di una con il medesimo l) e
!lb è il momento angolare della funzione di bessel
!si può avere solo lb=la+1 lb=la-1 o lb=la
!e questi tre casi sono rappresentati dall'ultimo indice
!3.71
   ndm = MAXVAL(upf(:)%kkbeta)
   allocate (aux(ndm))
   allocate (besr(ndm))

!
!inizializzazione di tabr
! questo può essere estratto. Il resto della routine è init_us_1.f90 di PW
!
   pref = fpi/sqrt(omega)
   call divide(intra_bgrp_comm, nqxq, startq, lastq)
   tabr(:, :, :, :) = 0.d0
   do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
         l = upf(nt)%lll(nb)
         do iq = startq, lastq
            qi = (iq - 1)*dq
!inizializzazione tabella con l_bessel=l
            !                  kkbeta dice dove bessel va a zero
            call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
            do ir = 1, upf(nt)%kkbeta
               aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir) ! nel file upf c'è x per il proiettorie
               !quindi qui c'è solo r^2)
            enddo
            call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
            tabr(iq, nb, nt, 0) = vqint*pref
!l_bessel=l+1
            call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l + 1, besr)
            do ir = 1, upf(nt)%kkbeta
               aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
            enddo
            call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
            tabr(iq, nb, nt, 1) = vqint*pref
!l_bessel=l-1, solo se l>0
            if (l > 0) then
               call sph_bes(upf(nt)%kkbeta, rgrid(nt)%r, qi, l - 1, besr)
               do ir = 1, upf(nt)%kkbeta
                  aux(ir) = upf(nt)%beta(ir, nb)*besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
               enddo
               call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
               tabr(iq, nb, nt, -1) = vqint*pref
            end if
         enddo
      enddo
   enddo
#ifdef __MPI
   call mp_sum(tabr, intra_bgrp_comm)
#endif
   ! initialize spline interpolation
   if (spline_ps) then
      CALL errore('init_us_1a', 'splines for tabr not implemented', 1)
   endif
!!!!!!!!!!!!!!
!inizializzazione di tablocal_hg 3.17 3.20
   !attenzione alla nuova griglia che deve essere compatibile con quella
!dello pseudo locale
   rm = MAXVAL(rgrid(:)%mesh)
   deallocate (besr)
   allocate (besr(rm))
   deallocate (aux)
   allocate (aux(rm))
!
   pref = fpi/omega
   call divide(intra_bgrp_comm, nqxq, startq, lastq)
   tablocal_hg(:, :, :) = 0.d0
   do nt = 1, ntyp
      do iq = startq, lastq
         qi = (iq - 1)*dq
!inizializzazione con l=0
         call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 0, besr)
         do ir = 1, rgrid(nt)%mesh
            aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                      besr(ir)*rgrid(nt)%r(ir)
         end do
         call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
         tablocal_hg(iq, nt, 0) = vqint*pref
!inizializzazione con l=1
         call sph_bes(rgrid(nt)%mesh, rgrid(nt)%r, qi, 1, besr)
         do ir = 1, rgrid(nt)%mesh
            aux(ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir) + 2.d0*zv(nt))* &
                      besr(ir)*rgrid(nt)%r(ir)*rgrid(nt)%r(ir)
         end do
         call simpson(rgrid(nt)%mesh, aux, rgrid(nt)%rab, vqint)
         tablocal_hg(iq, nt, 1) = vqint*pref
      end do
   end do
#ifdef __MPI
   call mp_sum(tablocal_hg, intra_bgrp_comm)
#endif
   ! initialize spline interpolation
   if (spline_ps) then
      CALL errore('init_us_1a', 'splines for tablocal not implemented', 1)
   endif
   deallocate (besr)
   deallocate (aux)

   call stop_clock('init_us_1a')
   return
end subroutine init_us_1a

