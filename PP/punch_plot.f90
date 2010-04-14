!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE punch_plot (filplot, plot_num, sample_bias, z, dz, &
     emin, emax, kpoint, kband, spin_component, lsign, epsilon)
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output several quantities
  !     in a real space 3D mesh for subsequent processing or plotting
  !     The integer variable plot_num is used to choose the output quantity
  !     See file Doc/INPUT_PP.* for a description of plotted quantities
  !
  !     The output quantity is written (formatted) on file filplot.
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev
  USE cell_base,        ONLY : at, bg, omega, alat, celldm, ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE printout_base,    ONLY : title 
  USE extfield,         ONLY : tefield, dipfield
  USE gvect
  USE klist,            ONLY : nks, nkstot, xk
  USE lsda_mod,         ONLY : nspin, current_spin
  USE ener,             ONLY : ehart
  USE io_global,        ONLY : stdout, ionode
  USE scf,              ONLY : rho, vltot, v
  USE wvfct,            ONLY : npw, nbnd, wg, igk
  USE noncollin_module, ONLY : noncolin
  USE fft_base,         ONLY : grid_gather
  USE paw_postproc,     ONLY : PAW_make_ae_charge

  IMPLICIT NONE
  CHARACTER(len=*) :: filplot
  INTEGER :: kpoint, kband, spin_component, plot_num
  LOGICAL :: lsign
  REAL(DP) :: sample_bias, z, dz, dummy
  REAL(DP) :: emin, emax, wf, charge, epsilon

  INTEGER :: is, ipol
#ifdef __PARA
  ! auxiliary vector (parallel case)
  REAL(DP), ALLOCATABLE :: raux1 (:)

#endif
  ! auxiliary vector
  REAL(DP), ALLOCATABLE :: raux (:)


  IF (filplot == ' ') RETURN
#ifdef __PARA
  ALLOCATE (raux1( nrx1 * nrx2 * nrx3))    
#endif

  WRITE( stdout, '(/5x,"Calling punch_plot, plot_num = ",i3)') plot_num
  IF (plot_num == 7 ) &
     WRITE( stdout, '(/5x,"Plotting k_point = ",i3,"  band =", i3  )') &
                                                   kpoint, kband
  IF (plot_num == 7 .AND. noncolin .AND. spin_component .NE. 0 ) &
     WRITE( stdout, '(/5x,"Plotting spin magnetization ipol = ",i3)') &
                                                          spin_component
  !
  ALLOCATE (raux( nrxx))    
  !
  !     Here we decide which quantity to plot
  !
  IF (plot_num == 0) THEN
     !
     !      plot of the charge density
     !
     IF (noncolin) THEN
        call dcopy (nrxx, rho%of_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (nrxx, rho%of_r (1, 1), 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (nrxx, 1.d0, rho%of_r (1, is), 1, raux, 1)
           ENDDO
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (nrxx, rho%of_r (1, current_spin), 1, raux, 1)
           CALL dscal (nrxx, 0.5d0 * nspin, raux, 1)
        ENDIF
     ENDIF

  ELSEIF (plot_num == 1) THEN
     !
     !       The total self-consistent potential V_H+V_xc on output
     !
     IF (noncolin) THEN
        call dcopy (nrxx, v%of_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (nrxx, v%of_r, 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (nrxx, 1.0d0, v%of_r (1, is), 1, raux, 1)
           ENDDO
           CALL dscal (nrxx, 1.d0 / nspin, raux, 1)
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (nrxx, v%of_r (1, current_spin), 1, raux, 1)
        ENDIF
     ENDIF
     CALL daxpy (nrxx, 1.0d0, vltot, 1, raux, 1)

  ELSEIF (plot_num == 2) THEN
     !
     !       The local pseudopotential on output
     !
     CALL dcopy (nrxx, vltot, 1, raux, 1)

  ELSEIF (plot_num == 3) THEN
     !
     !       The local density of states at e_fermi on output
     !
     if (noncolin) call errore('punch_plot','not implemented yet',2)
     CALL local_dos (6, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 4) THEN
     !
     !       The local density of electronic entropy on output
     !
     if (noncolin) call errore('punch_plot','not implemented yet',1)
     CALL local_dos (2, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 5) THEN

     if (noncolin) call errore('punch_plot','not implemented yet',1)
     CALL work_function (wf)
#ifdef __PARA
     CALL stm (wf, sample_bias, z, dz, raux1)
#else
     CALL stm (wf, sample_bias, z, dz, raux)
#endif
     WRITE (title, '(" Bias in eV = ",f10.4," # states",i4)') &
             sample_bias * rytoev, NINT (wf)

  ELSEIF (plot_num == 6) THEN
     !
     !      plot of the spin polarisation
     !
     IF (nspin == 2) THEN
        CALL dcopy (nrxx, rho%of_r (1, 1), 1, raux, 1)
        CALL daxpy (nrxx, - 1.d0, rho%of_r (1, 2), 1, raux, 1)
     ELSE
        raux(:) = 0.d0
     ENDIF

  ELSEIF (plot_num == 7) THEN

     IF (noncolin) THEN
        IF (spin_component==0) THEN
           CALL local_dos (0, lsign, kpoint, kband, spin_component, emin, emax, raux)
        ELSE
           CALL local_dos_mag (spin_component, kpoint, kband, raux)
        ENDIF
     ELSE
        !!!CALL local_dos (0, lsign, kpoint, kband, spin_component, emin, emax, raux)
     CALL write_all_states (0, lsign, kpoint, kband, spin_component, emin, emax, raux)
     END IF
  ELSEIF (plot_num == 8) THEN

     if (noncolin) &
        call errore('punch_plot','elf+noncolin not yet implemented',1)
     CALL do_elf (raux)

  ELSEIF (plot_num == 9) THEN

     call errore('punch_plot','no longer implemented, see PP/plan_avg.f90',1)

  ELSEIF (plot_num == 10) THEN

     CALL local_dos (3, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 11) THEN

     raux(:) = vltot(:) 
     IF (nspin == 2) THEN
        rho%of_g(:,1) =  rho%of_g(:,1) +  rho%of_g(:,2)
        rho%of_r (:,1) =  rho%of_r (:,1) +  rho%of_r (:,2)
        nspin = 1
     END IF
     CALL v_h (rho%of_g, ehart, charge, raux)
     IF (tefield.AND.dipfield) CALL add_efield(raux,dummy,rho%of_r,.true.)

  ELSEIF (plot_num == 12) THEN

     raux=0.d0
     IF (tefield) THEN
         CALL add_efield(raux,dummy,rho%of_r,.true.)
     ELSE
         CALL infomsg ('punch_plot','e_field is not calculated')
     ENDIF

  ELSEIF (plot_num == 13) THEN

     IF (noncolin) THEN
        IF (spin_component==0) THEN
           raux(:) = SQRT(rho%of_r(:,2)**2 + rho%of_r(:,3)**2 + rho%of_r(:,4)**2 )
        ELSEIF (spin_component >= 1 .OR. spin_component <=3) THEN
           raux(:) = rho%of_r(:,spin_component+1)
        ELSE
           CALL errore('punch_plot','spin_component not allowed',1)
        ENDIF
     ELSE
        CALL errore('punch_plot','noncollinear spin required',1)
     ENDIF

  ELSEIF (plot_num == 14 .OR. plot_num == 15 .OR. plot_num == 16 ) THEN

     ipol = plot_num - 13
     call polarization ( spin_component, ipol, epsilon, raux )

  ELSEIF (plot_num == 17) THEN
     write(stdout, '(7x,a)') "Reconstructing all-electron valence charge."
     ! code partially duplicate from plot_num=0, should be unified
     CALL PAW_make_ae_charge(rho)
     !
     IF (spin_component == 0) THEN
         CALL dcopy (nrxx, rho%of_r (1, 1), 1, raux, 1)
         DO is = 2, nspin
            CALL daxpy (nrxx, 1.d0, rho%of_r (1, is), 1, raux, 1)
         ENDDO
      ELSE
         IF (nspin == 2) current_spin = spin_component
         CALL dcopy (nrxx, rho%of_r (1, current_spin), 1, raux, 1)
         CALL dscal (nrxx, 0.5d0 * nspin, raux, 1)
      ENDIF
  ELSEIF (plot_num == 18) THEN

     IF (noncolin) THEN
        IF (spin_component==0) THEN
           raux(:) = SQRT(v%of_r(:,2)**2 + v%of_r(:,3)**2 + v%of_r(:,4)**2 )
        ELSEIF (spin_component >= 1 .OR. spin_component <=3) THEN
           raux(:) = v%of_r(:,spin_component+1)
        ELSE
           CALL errore('punch_plot','spin_component not allowed',1)
        ENDIF
     ELSE
        CALL errore('punch_plot','B_xc available only when noncolin=.true.',1)
     ENDIF
  ELSE

     CALL infomsg ('punch_plot', 'plot_num not implemented')

  ENDIF

#ifdef __PARA
  IF (.NOT. (plot_num == 5 ) ) CALL grid_gather (raux, raux1)
  IF ( ionode ) &
     CALL plot_io (filplot, title, nrx1, &
         nrx2, nrx3, nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, &
         dual, ecutwfc, plot_num, atm, ityp, zv, tau, raux1, + 1)
  DEALLOCATE (raux1)
#else

  CALL plot_io (filplot, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num, &
       atm, ityp, zv, tau, raux, + 1)

#endif

  DEALLOCATE (raux)
  RETURN
END SUBROUTINE punch_plot

SUBROUTINE polarization ( spin_component, ipol, epsilon, raux )
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : fpi
  USE gvect, ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3, nl, nlm, &
       ngm, nrxx, gstart, g, gg
  USE lsda_mod,  ONLY : nspin
  USE scf, ONLY: rho
  USE control_flags,    ONLY : gamma_only
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  INTEGER :: spin_component, ipol, ig
  REAL(DP) :: epsilon, raux (nrxx)
  !
  IF (ipol < 1 .OR. ipol > 3) CALL errore('polarization', &
       'wrong component',1)
  !
  IF (spin_component == 0) THEN
     IF (nspin == 1 .OR. nspin == 4 ) THEN
        psic(:) = CMPLX(rho%of_r(:,1), 0.d0,kind=DP)
     ELSE IF (nspin == 2) THEN
        psic(:) = CMPLX(rho%of_r(:,1) + rho%of_r(:,2), 0.d0,kind=DP) 
     END IF
  ELSE 
     IF (spin_component > nspin .OR. spin_component < 1) &
          CALL errore('polarization', 'wrong spin component',1)
     psic(:) = CMPLX(rho%of_r(:,spin_component), 0.d0,kind=DP)
  END IF
  !
  !   transform to G space
  !
  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  IF (gstart == 2) psic (1) = (epsilon - 1.d0) / fpi
  DO ig = gstart, ngm
     psic (nl (ig) ) = psic (nl (ig) ) * g (ipol, ig) / gg (ig) &
       / (0.d0, 1.d0)
     if (gamma_only) psic (nlm(ig) ) = CONJG ( psic (nl (ig) ) )
  END DO
  !
  CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  raux (:) =  DBLE (psic (:) )
  !
  RETURN
  !
END SUBROUTINE polarization


! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!--------------------------------------------------------------------
subroutine write_all_states(iflag, lsign, kpoint, kband, spin_component, &
                      emin, emax, dos)
  !--------------------------------------------------------------------
  !
  !     iflag=0: calculates |psi|^2 for band "kband" at point "kpoint"
  !     iflag=1: calculates the local density of state at e_fermi
  !              (only for metals)
  !     iflag=2: calculates the local density of  electronic entropy
  !              (only for metals with fermi spreading)
  !     iflag=3: calculates the integral of local dos from "emin" to "emax"
  !              (emin, emax in Ry)
  !
  !     lsign:   if true and k=gamma and iflag=0, write |psi|^2 * sign(psi)
  !     spin_component: for iflag=3 and LSDA calculations only
  !                     0 for up+down dos,  1 for up dos, 2 for down dos
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega, tpiba2
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : ef
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, g, ecutwfc
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk, nkstot
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symm_base,            ONLY : nsym, s, ftau
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, g2kin
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE io_files,             ONLY : iunwfc, nwordwfc, outdir
  USE mp_global,            ONLY : me_pool, nproc_pool, my_pool_id, npool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE becmod,               ONLY : calbec, becp, allocate_bec_type, &
                                   deallocate_bec_type
  implicit none
  !
  ! input variables
  !
  integer, intent(in) :: iflag, kpoint, kband, spin_component
  logical, intent(in) :: lsign
  real(DP), intent(in) :: emin, emax
  !
  real(DP), intent(out) :: dos (nrxx)
  !
  !    local variables
  !
  integer :: ikb, jkb, ijkb0, ih, jh, kh, na, ijh, np
  ! counters for US PPs
  integer :: ir, is, ig, ibnd, ik, irm, isup, isdw, ipol, kkb, is1, is2
  ! counters
!NANNI
  real(DP) :: sum
  ! temporary sum
  integer :: jbnd,iunit
  ! counters for output files with states
  character(LEN=80)  :: filename,filenameS
  character(LEN=5)   :: fileindex 
  ! names for output files with states
!END NANNI
  real(DP) :: w, w1, modulus
  real(DP), allocatable :: segno(:), maxmod(:)
!NANNI
  complex(DP), allocatable :: spsi(:,:)
!  complex(DP), allocatable :: spsic(:)
!END NANNI
  integer :: who_calculate, iproc
  complex(DP) :: phase 
  real(DP), external :: w0gauss, w1gauss
  logical :: i_am_the_pool
  integer :: which_pool, kpoint_pool
  !
  ! input checks
  !
  if (noncolin.and. lsign) call errore('write_all_states','not available',1)
  if (noncolin.and. gamma_only) call errore('write_all_states','not available',1)
  !
  if ( (iflag == 0) .and. ( kband < 1 .or. kband > nbnd ) ) &
       call errore ('write_all_states', 'wrong band specified', 1)
  if ( (iflag == 0) .and. ( kpoint < 1 .or. kpoint > nkstot ) ) &
       call errore ('write_all_states', 'wrong kpoint specified', 1)
  if (lsign) then
     if (iflag /= 0) call errore ('write_all_states', 'inconsistent flags', 1)
     if (sqrt(xk(1,kpoint)**2+xk(2,kpoint)**2+xk(3,kpoint)**2) > 1d-9 )  &
        call errore ('write_all_states', 'k must be zero', 1)
  end if
  !
  CALL allocate_bec_type ( nkb, nbnd, becp )
  rho%of_r(:,:) = 0.d0
  dos(:) = 0.d0
  becsum(:,:,:) = 0.d0
  if (lsign) allocate(segno(nrxx))
!NANNI
  allocate(spsi(npwx,nbnd)) 
!  allocate(spsic(nrxx)) 
!END NANNI

  !
  !   calculate the correct weights
  !
  if (iflag /= 0 .and. .not.lgauss) call errore ('write_all_states', &
       'gaussian broadening needed', 1)
  if (iflag == 2 .and. ngauss /= -99) call errore ('write_all_states', &
       ' beware: not using Fermi-Dirac function ',  - ngauss)
  do ik = 1, nks
     do ibnd = 1, nbnd
        if (iflag == 0) then
           wg (ibnd, ik) = 0.d0
        elseif (iflag == 1) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss) / degauss
        elseif (iflag == 2) then
           wg (ibnd, ik) = - wk (ik) * w1gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss)
        elseif (iflag == 3) then
           if (et (ibnd, ik) <=  emax .and. et (ibnd, ik) >= emin) then
              wg (ibnd, ik) = wk (ik)
           else
              wg (ibnd, ik) = 0.d0
           endif
        else
           call errore ('write_all_states', ' iflag not allowed', abs (iflag) )
        endif
     enddo
  enddo

  IF (npool>1) THEN
     CALL xk_pool( kpoint, nkstot, kpoint_pool,  which_pool )
     if (kpoint_pool<1 .or. kpoint_pool> nks) &
        CALL errore('write_all_states','problems with xk_pool',1)
     i_am_the_pool=(my_pool_id==which_pool)
  ELSE
     i_am_the_pool=.true.
     kpoint_pool=kpoint
  ENDIF

  if (iflag == 0.and.i_am_the_pool) wg (kband, kpoint_pool) = 1.d0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the density of states
  !
  do ik = 1, nks
     if (ik == kpoint_pool .and.i_am_the_pool.or. iflag /= 0) then
        if (lsda) current_spin = isk (ik)
        call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call init_us_2 (npw, igk, xk (1, ik), vkb)
!
        call calbec ( npw, vkb, evc, becp )
     !
     !     here we computed the density of states
! Ho toccato solo da qui fino a FINE_NANNI
     !     now for each state ibnd, the wfc is printed in state_XXXXX.bin
     !     and S(psi_i) is printed in Sstate_XXXXX.bin
     ! noncolin case only, gamma_only case only
     !
        iunit=70
        jbnd=nbnd
!PAOLO: se entra qui si pianta (segmentation fault). Come mai?
 print *,'npwx = ',npwx
 print *,'npw  = ',npw
 print *,'nbnd = ',nbnd
        CALL s_psi( npwx, npw, nbnd, evc, spsi )
        do ibnd = 1, nbnd
              jbnd=jbnd-1
              write(fileindex, fmt='(i5.5)') jbnd
              filename = 'state_' // fileindex // '.bin'
              open (iunit, file=filename, status='unknown', form='unformatted')
!              open (iunit, file=filename, status='unknown', form='formatted')
              psic(1:nrxxs) = (0.d0,0.d0)
              do ig = 1, npw
                 psic (nls (igk (ig) ) ) = evc (ig, ibnd)
              enddo
              do ig = 1, npw
                 psic (nlsm(igk (ig) ) ) = CONJG(evc (ig, ibnd))
              enddo
              call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
              w1 = wg (ibnd, ik) / omega
!
!  save the wavefunction at the gamma point
                 !  psi(r) is real by construction
!Se stampo qui, la vedo reale, ok...
!interpolate funziona solo con reali
!quindi prendo la parte reale di psic usando segno
                 segno(1:nrxxs) = DBLE(psic(1:nrxxs))
                 if (doublegrid) call interpolate (segno, segno, 1)
!a questo punto e' reale per forza...
!stampa psi_i
                 write(iunit) (segno(ir),ir=1,nrxx)
!                 print *,'w1= ',w1
!                 write(iunit,'(f)') (segno(ir),ir=1,nrxx)
!                 close(iunit)
!inizio S(psi)
!              spsic(1:nrxxs) = (0.d0,0.d0)
!              do ig = 1, npw
!                 spsic (nls (igk (ig) ) ) = spsi (ig, ibnd)
!              enddo
!              do ig = 1, npw
!                 spsic (nlsm(igk (ig) ) ) = CONJG(spsi (ig, ibnd))
!              enddo
!              call cft3s (spsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!              if (doublegrid) call interpolate (spsic, spsic, 1)
!calcolo l'integrale di psi^* S(psi) su tutto lo spazio
!per adesso provo con psi^* psi
              sum=0.d0
              do ir=1,nrxx
               sum = sum + segno(ir)*segno(ir) 
              enddo
              print *,'ibnd, sum= ',ibnd,sum/nr1/nr2/nr3
              !
        enddo !enddo ibnd
     endif !endif ik
  enddo !enddo ik
!FINE_NANNI

  CALL deallocate_bec_type ( becp ) 

  if (doublegrid) then
     if (noncolin) then
       call interpolate(rho%of_r, rho%of_r, 1)
     else
       do is = 1, nspin
         call interpolate(rho%of_r(1, is), rho%of_r(1, is), 1)
       enddo
     endif
  endif

  return

end subroutine write_all_states
