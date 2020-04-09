subroutine init_hartree()
   use qpoint, only : ikqs
   implicit none

   if ( .not. allocated(ikqs)) allocate(ikqs(1))
   ikqs(1)=1

end subroutine


subroutine routine_hartree()
   use kinds, only: DP
   use hartree_mod
   use wvfct, only: nbnd, npw, npwx
   use io_files, only: nwordwfc, diropn, iunwfc, prefix, tmp_dir
   use fft_base, only: dffts
   use wavefunctions, only: psic, evc
   use mp, only: mp_sum, mp_barrier
   use io_global, only: stdout, ionode
   use gvect, only: g, ngm, gstart
   !USE gvecs, ONLY: dffts%nl, dffts%nlm   !moved to fft type
   USE cell_base, ONLY: tpiba, omega, tpiba2, alat, at
   USE constants, ONLY: e2, fpi, pi
   USE fft_interfaces, ONLY: fwfft, invfft
   use mp_global, only: mp_global_end
   use ions_base, only: tau
   use gvecw, only: ecutwfc
!
   use uspp, ONLY: vkb, nkb
   use klist, only: xk, igk_k
   use wvfct, ONLY: g2kin, et
   use becmod
   use constants, only: rytoev
   USE eqv, ONLY: dpsi, dvpsi
   USE mp_pools, ONLY: intra_pool_comm
   
   implicit none
   
   real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)
   !complex(kind=DP), allocatable :: evc_uno(:, :), evc_due(:, :), evp(:, :), tmp(:, :)
   complex(kind=DP), allocatable ::  evp(:, :), tmp(:, :)
   integer :: iun, iv, igm, ibnd, i
   logical :: exst
   real(kind=DP), allocatable ::charge(:), fac(:)
   real(kind=DP), allocatable ::excharge_r(:), exgradcharge_r(:, :), exdotcharge_r(:)
   real(kind=DP) :: qq_fact, q_tot
   real(kind=DP) :: update(1:3), update_a(1:3), update_b(1:3)
   real(kind=DP) :: amodulus
!
   complex(kind=DP), allocatable ::charge_g(:), v_uno(:), v_due(:), v_point(:), v_mean(:), charge_g_due(:)
!
   complex(kind=DP), allocatable ::exgradcharge_g(:, :)
!
   integer ::icoord, enne, ir
!
   logical  :: l_test
   real(DP) :: M(nbnd, nbnd), emme(nbnd, nbnd), kcurrent(3), s(nbnd, nbnd), ss(nbnd, nbnd), &
               kcurrent_a(3), kcurrent_b(3), ecurrent(3)
   integer  :: inbd, jbnd, ig, ipol
   integer, external :: find_free_unit

   write (stdout, *) 'INIZIO ROUTINE HARTREE & KOHN'
   npw=npwx ! only gamma
   call start_clock('routine_hartree')
   call start_clock('hartree_current')
   call init_hartree()
!-------------------------------HARTREE---------------------------------------------------------
!call flush_unit(stdout)
!Il calcolo della corrente di Hartree e quella di Kohn-Sham
!richiedono che in memoria esistano contemporaneamente le funzioni d'onda a due tempi diversi t e t'= t-Dt,
!dove Dt e' l'intervallo usato per la derivata numerica.
!Quantita' con indice uno o prive di indici si riferiscono a quantita' calcolate al tempo t (i.e. prefix,charge_g,evc_uno)
!e l'indice due si riferisce a quantita' calcolate al tempo t'=t-Dt (i.e. prefix_due,charge_due,evc_due)

!-------STEP1: lettura e allocazione (npwx = number of plane waves (npwx>npw per comodita'), nbnd = number of bands (n_electrons/2 per insulators))
   allocate (tmp(npwx, nbnd))
   allocate (evp(npwx, nbnd))
   !allocate (evc_uno(npwx, nbnd))
   !allocate (evc_due(npwx, nbnd))
   allocate (charge_g(ngm))
   allocate (charge_g_due(ngm))
   allocate (charge(dffts%nnr))
   allocate (fac(ngm))
   allocate (v_uno(ngm))
   allocate (v_due(ngm))
   allocate (v_point(ngm))
   allocate (v_mean(ngm))
!
   allocate (excharge_r(dffts%nnr))
   allocate (exgradcharge_r(3, dffts%nnr))
   allocate (exgradcharge_g(3, ngm))
   allocate (exdotcharge_r(dffts%nnr))
!
! legge le funzioni d'onda uno e due (evc_uno, evc_due)
!   iun = find_free_unit()
!   call diropn_due(prefix_due, iun, 'wfc', 2*nwordwfc, exst, tmp_dir)
!   call davcio(evc_due, 2*nwordwfc, iun, 1, -1)
!   close (iun)
!
!   close (iunwfc)
!   call diropn(iunwfc, 'wfc', 2*nwordwfc, exst, tmp_dir)
!   call davcio(evc_uno, 2*nwordwfc, iunwfc, 1, -1)
!   close (iunwfc)

!-------STEP2.1: inizializzazione di charge_g, la carica al tempo t.
! charge_g = densità di carica nello spazio reciproco. Si ottiene dalla FFT di |evc(r)|^2, dove evc(r)=IFFT(evc)
! per ottimizzare: fa 2 bande alla volta con una sola IFFT


   charge = 0.d0
   do iv = 1, nbnd, 2
      psic = 0.d0
      if (iv == nbnd) then
         psic(dffts%nl(1:npw)) = evc_uno(1:npw, iv)
         psic(dffts%nlm(1:npw)) = CONJG(evc_uno(1:npw, iv))
      else
         psic(dffts%nl(1:npw)) = evc_uno(1:npw, iv) + (0.D0, 1.D0)*evc_uno(1:npw, iv + 1)
         psic(dffts%nlm(1:npw)) = CONJG(evc_uno(1:npw, iv) - (0.D0, 1.D0)*evc_uno(1:npw, iv + 1))
      end if
      call invfft('Wave', psic, dffts)
      charge(1:dffts%nnr) = charge(1:dffts%nnr) + dble(psic(1:dffts%nnr))**2.0
      if (iv /= nbnd) then
         charge(1:dffts%nnr) = charge(1:dffts%nnr) + dimag(psic(1:dffts%nnr))**2.0
      end if
   end do
   q_tot = 0.
   do i = 1, dffts%nnr
      q_tot = q_tot + charge(i)
   end do
   q_tot = q_tot/(dffts%nr1*dffts%nr2*dffts%nr3)
   call mp_sum(q_tot, intra_pool_comm)
   IF (ionode) THEN
      print *, 'check_charge', q_tot
   ENDIF
!moltiplico per due causa degenerazione di spin
   charge(1:dffts%nnr) = charge(1:dffts%nnr)*2.d0
!
!calcolo carica in spazio reciproco (FFT di psic)
   psic = 0.d0
   psic(1:dffts%nnr) = dcmplx(charge(1:dffts%nnr), 0.d0)
   call fwfft('Rho', psic, dffts) ! TODO: smooth does not exist anymore
   charge_g(1:ngm) = psic(dffts%nl(1:ngm))

!!!!!!!!!!!!------------primo exchange-corr_intermezzo 1/2  -------------!!!!!!!!!!!!!!!!!!
!copio la densità di carica in spazio reale
   excharge_r(1:dffts%nnr) = charge(1:dffts%nnr)
   do icoord = 1, 3
!
      do igm = 1, ngm
         exgradcharge_g(icoord, igm) = charge_g(igm)*(0.d0, 1.d0)*g(icoord, igm)*tpiba
      end do
!
      psic = 0.d0
      psic(dffts%nl(1:ngm)) = exgradcharge_g(icoord, 1:ngm)
      psic(dffts%nlm(1:ngm)) = CONJG(exgradcharge_g(icoord, 1:ngm))
!
      call invfft('Rho', psic, dffts)
!
      exgradcharge_r(icoord, 1:dffts%nnr) = dble(psic(1:dffts%nnr))
   end do
!
!!!!!!!!!!!!------------- fine intermezzo----------------- !!!!!!!!!!!!!!!!!!!!

!-------STEP2.2-------inizializzazione di chargeg_due, la carica al tempo t-Dt.
   charge = 0.d0
   do iv = 1, nbnd, 2
      psic = 0.d0
      if (iv == nbnd) then
         psic(dffts%nl(1:npw)) = evc_due(1:npw, iv)
         psic(dffts%nlm(1:npw)) = CONJG(evc_due(1:npw, iv))
      else
         psic(dffts%nl(1:npw)) = evc_due(1:npw, iv) + (0.D0, 1.D0)*evc_due(1:npw, iv + 1)
         psic(dffts%nlm(1:npw)) = CONJG(evc_due(1:npw, iv) - (0.D0, 1.D0)*evc_due(1:npw, iv + 1))
      end if
      call invfft('Wave', psic, dffts)
      charge(1:dffts%nnr) = charge(1:dffts%nnr) + dble(psic(1:dffts%nnr))**2.0
      if (iv /= nbnd) then
         charge(1:dffts%nnr) = charge(1:dffts%nnr) + dimag(psic(1:dffts%nnr))**2.0
      end if
   end do
   q_tot = 0.
   do i = 1, dffts%nnr
      q_tot = q_tot + charge(i)
   end do
   q_tot = q_tot/(dffts%nr1*dffts%nr2*dffts%nr3)
   call mp_sum(q_tot, intra_pool_comm)
   if (ionode) print *, 'check_charge', q_tot
!moltiplico per due causa degenerazione di spin
   charge(1:dffts%nnr) = charge(1:dffts%nnr)*2.d0

!!!!!!!!!!!!------------secondo exchange-corr_intermezzo 2/2 -------------!!!!!!!!!!!!!!!!!!
!
   exdotcharge_r(1:dffts%nnr) = (excharge_r(1:dffts%nnr) - charge(1:dffts%nnr))/delta_t
!
!!!!!!!!!!!!------------- fine intermezzo----------------- !!!!!!!!!!!!!!!!!!!!

!
!calcolo carica due in spazio reciproco (FFT di psic)
   psic = 0.d0
   psic(1:dffts%nnr) = dcmplx(charge(1:dffts%nnr), 0.d0)
   call fwfft('Rho', psic, dffts)
   charge_g_due(1:ngm) = psic(dffts%nl(1:ngm))

!-------STEP3----- calcolo dei potenziali di Hartree a partire dalle cariche appena trovate.

!calcolo del potenziale v_uno e di fac
! fac(G) = e2*fpi/(tpiba2*G^2*omega)
! v(G) = charge(G)*fac
   do igm = 1, ngm
      qq_fact = g(1, igm)**2.d0 + g(2, igm)**2.d0 + g(3, igm)**2.d0
      if (qq_fact > 1.d-8) then
         fac(igm) = (e2*fpi/(tpiba2*qq_fact))
      else
         fac(igm) = 0.d0
      end if
   end do
   fac(:) = fac(:)/omega

   do igm = 1, ngm
      v_uno(igm) = charge_g(igm)*fac(igm)
   end do

!medesimo codice per il calcolo di v_due
   do igm = 1, ngm
      v_due(igm) = charge_g_due(igm)*fac(igm)
   end do

!-------STEP4----- derivate numeriche dei potenziali

!A questo punto si esegue il calcolo di v_point e v_mean
   do igm = 1, ngm
      v_point(igm) = (v_uno(igm) - v_due(igm))/delta_t
   end do
   do igm = 1, ngm
      v_mean(igm) = (v_uno(igm) + v_due(igm))/2.d0
   end do

!-------STEP 5----- applicazione formula finale, formula 1.17 del manuale.

   J_hartree = 0.d0
   do igm = gstart, ngm
      J_hartree(1:3) = J_hartree(1:3) - omega/(4.d0*pi)*dimag(CONJG(v_point(igm))*v_mean(igm))*g(1:3, igm)*tpiba
   end do
   call mp_sum(J_hartree, intra_pool_comm)

   call stop_clock('hartree_current')
   call print_clock('hartree_current')
   if (ionode) print *, 'CORRENTE HARTREE CALCOLATA'

!-----------------EXCHANGE-CORRELATION-----------------------------------------------
   call start_clock('xc_current')
!
!!!!!questo passaggio e' necessario? discuss with Paolo!!!!!
!
   exdotcharge_r(1:dffts%nnr) = exdotcharge_r(1:dffts%nnr)/omega
   excharge_r(1:dffts%nnr) = excharge_r(1:dffts%nnr)/omega
   exgradcharge_r(1:3, 1:dffts%nnr) = exgradcharge_r(1:3, 1:dffts%nnr)/omega
!
   J_xc(1:3) = 0.d0
!
   do ir = 1, dffts%nnr
      if (excharge_r(ir) > 1.0E-10_DP) then
         call pbex_current(abs(excharge_r(ir)), exgradcharge_r(1:3, ir), 1, update_a(1:3))
         call pbec_current(abs(excharge_r(ir)), exgradcharge_r(1:3, ir), 1, update_b(1:3))
         update(1:3) = e2*(update_a(1:3) + update_b(1:3))
         update(1:3) = -update(1:3)*exdotcharge_r(ir)*excharge_r(ir)
         J_xc(1:3) = J_xc(1:3) + update(1:3)
      end if
   end do
!
   call mp_sum(J_xc, intra_pool_comm)
!
!elemento di volume
!
   J_xc(1:3) = J_xc(1:3)*omega/(dffts%nr1*dffts%nr2*dffts%nr3)
!

   call stop_clock('xc_current')
   call print_clock('xc_current')
   if (ionode) print *, 'CORRENTE X-C CALCOLATA'

!---------------------------------KOHN------------------------------------------------
   call start_clock('kohn_current')
   allocate (dpsi(npwx, nbnd))
   allocate (dvpsi(npwx, nbnd))
! per il preconditioning:
   do ig = 1, npw
      igk_k(ig, 1) = ig
      g2kin(ig) = ((xk(1, 1) + g(1, igk_k(ig, 1)))**2 + &
                   (xk(2, 1) + g(2, igk_k(ig, 1)))**2 + &
                   (xk(3, 1) + g(3, igk_k(ig, 1)))**2)*tpiba2
   end do
! inizializzazioni potenziali che servono per fare  H|psi>
   call init_us_1()
   call init_us_2(npw, igk_k(1, 1), xk(1, 1), vkb)
   call allocate_bec_type(nkb, nbnd, becp)
   call calbec(npw, vkb, evc, becp)

   l_test = .true.
   if (l_test) then
      s = 0.d0
! calcola s = < evc_due, evc_uno >, toglie un contributo a G=0
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evc_due, 2*npwx, evc_uno, 2*npwx, 0.d0, s, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               s(ibnd, jbnd) = s(ibnd, jbnd) - dble(conjg(evc_due(1, ibnd))*evc_uno(1, jbnd))
            end do
         end do
      end if
      call mp_sum(s, intra_pool_comm)
! calcola phi_v^c_punto usando <s,s>
      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, s, nbnd, s, nbnd, 0.d0, ss, nbnd)
      evp = 0.d0
      do ibnd = 1, nbnd
         do jbnd = 1, nbnd
            do ig = 1, npw
               evp(ig, ibnd) = evp(ig, ibnd) - evc_due(ig, jbnd)*s(jbnd, ibnd)
               evp(ig, ibnd) = evp(ig, ibnd) + evc_uno(ig, jbnd)*ss(ibnd, jbnd)
            end do
         end do
      end do
      evp(:, :) = evp(:, :)/delta_t
   end if



   tmp(:, :) = (0.d0, 0.d0)
!fa H|evp>  cioe' H|phi_v^c_punto> e lo mette in tmp:
   call h_psi(npwx, npw, nbnd, evp, tmp)
   kcurrent = 0.d0
   kcurrent_a = 0.d0
   kcurrent_b = 0.d0
   ecurrent = 0.d0
   polariz: do ipol = 1, 3
      call start_clock('project')
      ! calcola la proiezione sulla banda di conduzione
      call project(ipol)
      call stop_clock('project')
      call print_clock('project')

      call start_clock('kohn-first')
      emme = 0.d0
      ! prodotto scalare
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, dvpsi, 2*npwx, tmp, 2*npwx, 0.d0, emme, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(dvpsi(1, ibnd))*tmp(1, jbnd))
            end do
         end do
      end if
      call mp_sum(emme, intra_pool_comm)

      do ibnd = 1, nbnd
         kcurrent(ipol) = kcurrent(ipol) + emme(ibnd, ibnd)
         kcurrent_a(ipol) = kcurrent_a(ipol) + emme(ibnd, ibnd)
      end do
      call stop_clock('kohn-first')
      call print_clock('kohn-first')

      call start_clock('kohn-second')
      emme = 0.d0
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evp, 2*npwx, dvpsi, 2*npwx, 0.d0, emme, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(evp(1, ibnd))*dvpsi(1, jbnd))
            end do
         end do
      end if
      call mp_sum(emme, intra_pool_comm)

      ! et contiene gli autovalori
      do ibnd = 1, nbnd
         kcurrent(ipol) = kcurrent(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
         kcurrent_b(ipol) = kcurrent_b(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
         ecurrent(ipol) = ecurrent(ipol) + emme(ibnd, ibnd)
      end do

      call stop_clock('kohn-second')
      call print_clock('kohn-second')

      if (ionode) write (*, "('  KOHN POLARIZZAZIONE ',I3,' COMPLETATA')") ipol



   end do polariz
   J_kohn = 0.d0
   J_kohn_a = 0.d0
   J_kohn_b = 0.d0
   J_electron = 0.d0
   do ipol = 1, 3
!     at(:, ipol) / amodulus is the versor along direction ipol
!      amodulus = sqrt(at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2)
      amodulus = at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2
      J_kohn(:) = J_kohn(:) + 2.d0*at(:, ipol)*real(kcurrent(ipol)) / amodulus
      J_kohn_a(:) = J_kohn_a(:) + 2.d0*at(:, ipol)*real(kcurrent_a(ipol)) / amodulus
      J_kohn_b(:) = J_kohn_b(:) + 2.d0*at(:, ipol)*real(kcurrent_b(ipol)) / amodulus
      J_electron(:) = J_electron(:) + 2.d0*e2*at(:, ipol)*real(ecurrent(ipol)) / amodulus
!    J_kohn(ipol)=J_kohn(ipol)+2.d0*real(kcurrent(ipol))
   end do
   call stop_clock('kohn_current')
   call print_clock('kohn_current')
   if (ionode) print *, 'CORRENTE KOHN CALCOLATA'
!---------------------------------------------------------------------------
!!!!!!!!!! temporaneamente tolto per test

   if (ionode) then
      iun = find_free_unit()
      open (iun, file=trim(file_output), position='append')
      write (iun, *) 'Passo: '
      write (iun, '(A,10E20.12)') 'h&K-XC', J_xc(:)
      write (iun, '(A,10E20.12)') 'h&K-H', J_hartree(:)
      write (iun, '(A,1F15.7,9E20.12)') 'h&K-K', delta_t, J_kohn(1:3), J_kohn_a(1:3), J_kohn_b(1:3)
      write (iun, '(A,3E20.12)') 'h&K-ELE', J_electron(1:3)
      close (iun)
   end if


   !deallocate (evc_uno)
   !deallocate (evc_due)
   deallocate (charge)
   deallocate (charge_g)
   deallocate (fac)
   deallocate (v_uno)
   deallocate (v_due)
   deallocate (v_point)
   deallocate (v_mean)
   call deallocate_bec_type(becp)
   deallocate (dpsi)
   deallocate (dvpsi)
   deallocate (excharge_r)
   deallocate (exgradcharge_r)
   deallocate (exgradcharge_g)
   deallocate (exdotcharge_r)

   call stop_clock('routine_hartree')
   call print_clock('routine_hartree')
end subroutine routine_hartree

