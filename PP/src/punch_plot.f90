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
     emin, emax, kpoint, kband, spin_component, lsign)
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
  USE run_info,         ONLY : title
  USE extfield,         ONLY : tefield, dipfield
  USE fft_base,         ONLY : dfftp
  USE scatter_mod,      ONLY : gather_grid
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : gcutm
  USE gvecs,            ONLY : dual
  USE klist,            ONLY : nks, nkstot, xk
  USE lsda_mod,         ONLY : nspin, current_spin
  USE ener,             ONLY : ehart
  USE io_global,        ONLY : stdout, ionode
  USE scf,              ONLY : rho, vltot, v
  USE wvfct,            ONLY : nbnd, wg
  USE gvecw,            ONLY : ecutwfc
  USE noncollin_module, ONLY : noncolin
  USE paw_postproc,     ONLY : PAW_make_ae_charge

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: filplot
  INTEGER, INTENT(IN) :: plot_num, kpoint, kband, spin_component
  LOGICAL, INTENT(IN) :: lsign
  REAL(DP), INTENT(IN) :: sample_bias, z, dz, &
      emin, emax
  REAL(DP) :: dummy, charge
  INTEGER :: is, ipol, istates
#if defined(__MPI)
  ! auxiliary vector (parallel case)
  REAL(DP), ALLOCATABLE :: raux1 (:)
#endif
  ! auxiliary vector
  REAL(DP), ALLOCATABLE :: raux (:), raux2(:,:)


  IF (filplot == ' ') RETURN
#if defined(__MPI)
  ALLOCATE (raux1(  dfftp%nr1x *  dfftp%nr2x *  dfftp%nr3x))
#endif

  WRITE( stdout, '(/5x,"Calling punch_plot, plot_num = ",i3)') plot_num
  IF (plot_num == 3 ) &
     WRITE(stdout, '(/5x,"Energy =", f10.5, " eV, broadening =", f10.5, "eV" )') &
                       emin * rytoev, emax * rytoev
  IF (plot_num == 7) &
      WRITE( stdout, '(/5x,"Plotting k_point = ",i3,"  band =", i3  )') &
                                                   kpoint, kband
  IF ((plot_num == 7) .and. noncolin .and. spin_component /= 0 ) &
     WRITE( stdout, '(/5x,"Plotting spin magnetization ipol = ",i3)') &
                                                          spin_component
  !
  ALLOCATE (raux(dfftp%nnr))
  !
  !     Here we decide which quantity to plot
  !
  IF (plot_num == 0) THEN
     !
     !      plot of the charge density
     !
     IF (noncolin) THEN
        CALL dcopy (dfftp%nnr, rho%of_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
           ENDDO
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (dfftp%nnr, rho%of_r (1, current_spin), 1, raux, 1)
           CALL dscal (dfftp%nnr, 0.5d0 * nspin, raux, 1)
        ENDIF
     ENDIF

  ELSEIF (plot_num == 1) THEN
     !
     !       The total self-consistent potential V_H+V_xc on output
     !
     IF (noncolin) THEN
        CALL dcopy (dfftp%nnr, v%of_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (dfftp%nnr, v%of_r, 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (dfftp%nnr, 1.0d0, v%of_r (1, is), 1, raux, 1)
           ENDDO
           CALL dscal (dfftp%nnr, 1.d0 / nspin, raux, 1)
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (dfftp%nnr, v%of_r (1, current_spin), 1, raux, 1)
        ENDIF
     ENDIF
     CALL daxpy (dfftp%nnr, 1.0d0, vltot, 1, raux, 1)

  ELSEIF (plot_num == 2) THEN
     !
     !       The local pseudopotential on output
     !
     CALL dcopy (dfftp%nnr, vltot, 1, raux, 1)

  ELSEIF (plot_num == 3) THEN
     !
     !       The local density of states at emin, with broadening emax
     !
     WRITE (title, '(" Energy = ",f8.4," eV, ", "broadening = ",f8.4," eV")') &
                   emin * rytoev, emax * rytoev
     IF (noncolin) CALL errore('punch_plot','not implemented yet',1)
     CALL local_dos(1, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 4) THEN
     !
     !       The local density of electronic entropy on output
     !
     IF (noncolin) CALL errore('punch_plot','not implemented yet',1)
     CALL local_dos (2, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 5) THEN

     IF (noncolin) CALL errore('punch_plot','not implemented yet',1)
#if defined(__MPI)
     CALL stm (sample_bias, raux1, istates)
#else
     CALL stm (sample_bias, raux,  istates)
#endif
     WRITE (title, '(" Bias in eV = ",f10.4," # states",i4)') &
             sample_bias * rytoev, istates

  ELSEIF (plot_num == 6) THEN
     !
     !      plot of the spin polarisation
     !
     IF (nspin == 2) THEN
        CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
        CALL daxpy (dfftp%nnr, - 1.d0, rho%of_r (1, 2), 1, raux, 1)
     ELSE
        raux(:) = 0.d0
     ENDIF

  ELSEIF (plot_num == 7) THEN
     WRITE (title, '("k_point ",i4,", band ",i4)') kpoint ,kband

     IF (noncolin) THEN
        IF (spin_component==0) THEN
           CALL local_dos (0, lsign, kpoint, kband, spin_component, emin, emax, raux)
        ELSE
           CALL local_dos_mag (spin_component, kpoint, kband, raux)
        ENDIF
     ELSE
        CALL local_dos (0, lsign, kpoint, kband, spin_component, emin, emax, raux)
     ENDIF
  ELSEIF (plot_num == 8) THEN

     IF (noncolin) &
        CALL errore('punch_plot','elf+noncolin not yet implemented',1)
     CALL do_elf (raux)

  ELSEIF (plot_num == 9) THEN
     !
     !      plot of the charge density minus the atomic rho
     !
     allocate (raux2(dfftp%nnr,nspin))
     raux2 = 0.d0
     call atomic_rho(raux2, nspin)
     rho%of_r(:,:) = rho%of_r(:,:) - raux2
     deallocate (raux2)

     IF (noncolin) THEN
        CALL dcopy (dfftp%nnr, rho%of_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
           ENDDO
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (dfftp%nnr, rho%of_r (1, current_spin), 1, raux, 1)
           CALL dscal (dfftp%nnr, 0.5d0 * nspin, raux, 1)
        ENDIF
     ENDIF

  ELSEIF (plot_num == 10) THEN

     CALL local_dos (3, lsign, kpoint, kband, spin_component, emin, emax, raux)

  ELSEIF (plot_num == 11) THEN

     raux(:) = vltot(:)
     IF (nspin == 2) THEN
        rho%of_g(:,1) =  rho%of_g(:,1) +  rho%of_g(:,2)
        rho%of_r (:,1) =  rho%of_r (:,1) +  rho%of_r (:,2)
        nspin = 1
     ENDIF
     CALL v_h (rho%of_g, ehart, charge, raux)
     IF (tefield.and.dipfield) CALL add_efield(raux,dummy,rho%of_r,.true.)

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
           raux(:) = sqrt(rho%of_r(:,2)**2 + rho%of_r(:,3)**2 + rho%of_r(:,4)**2 )
        ELSEIF (spin_component >= 1 .or. spin_component <=3) THEN
           raux(:) = rho%of_r(:,spin_component+1)
        ELSE
           CALL errore('punch_plot','spin_component not allowed',1)
        ENDIF
     ELSE
        CALL errore('punch_plot','noncollinear spin required',1)
     ENDIF

  ELSEIF (plot_num == 14 .or. plot_num == 15 .or. plot_num == 16 ) THEN

     CALL errore('punch_plot','polarization no longer implemented',1)
     ! ipol = plot_num - 13
     ! CALL polarization ( spin_component, ipol, epsilon, raux )

  ELSEIF (plot_num == 17 .or. plot_num == 21) THEN
     WRITE(stdout, '(7x,a)') "Reconstructing all-electron valence charge."
     ! code partially duplicate from plot_num=0, should be unified
     CALL init_us_1()
     CALL PAW_make_ae_charge(rho,(plot_num==21))
     !
     IF (spin_component == 0) THEN
         CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
         DO is = 2, nspin
            CALL daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
         ENDDO
      ELSE
         IF (nspin == 2) current_spin = spin_component
         CALL dcopy (dfftp%nnr, rho%of_r (1, current_spin), 1, raux, 1)
         CALL dscal (dfftp%nnr, 0.5d0 * nspin, raux, 1)
      ENDIF
  ELSEIF (plot_num == 18) THEN

     IF (noncolin) THEN
        IF (spin_component==0) THEN
           raux(:) = sqrt(v%of_r(:,2)**2 + v%of_r(:,3)**2 + v%of_r(:,4)**2 )
        ELSEIF (spin_component >= 1 .or. spin_component <=3) THEN
           raux(:) = v%of_r(:,spin_component+1)
        ELSE
           CALL errore('punch_plot','spin_component not allowed',1)
        ENDIF
     ELSE
        CALL errore('punch_plot','B_xc available only when noncolin=.true.',1)
     ENDIF

  ELSEIF (plot_num == 19) THEN
     !
     ! Reduced density gradient
     !
     IF (noncolin) CALL errore('punch_plot','rdg+noncolin not yet implemented',1)
     CALL do_rdg (raux)           ! in elf.f90

  ELSEIF (plot_num == 20) THEN
     !
     ! Density * second eigenvalue of Hessian of density (for coloring RDG plots)
     !
     IF (noncolin) CALL errore('punch_plot','rdg+noncolin not yet implemented',1)
     CALL do_sl2rho (raux)        ! in elf.f90

  ELSEIF (plot_num == 22) THEN
     !
     !      plot of the kinetic charge density
     !
     IF (noncolin) THEN
        CALL dcopy (dfftp%nnr, rho%kin_r, 1, raux, 1)
     ELSE
        IF (spin_component == 0) THEN
           CALL dcopy (dfftp%nnr, rho%kin_r (1, 1), 1, raux, 1)
           DO is = 2, nspin
              CALL daxpy (dfftp%nnr, 1.d0, rho%kin_r (1, is), 1, raux, 1)
           ENDDO
        ELSE
           IF (nspin == 2) current_spin = spin_component
           CALL dcopy (dfftp%nnr, rho%kin_r (1, current_spin), 1, raux, 1)
           CALL dscal (dfftp%nnr, 0.5d0 * nspin, raux, 1)
        ENDIF
     ENDIF

  ELSE

     CALL infomsg ('punch_plot', 'plot_num not implemented')

  ENDIF

#if defined(__MPI)
  IF (.not. (plot_num == 5 ) ) CALL gather_grid (dfftp, raux, raux1)
  IF ( ionode ) &
     CALL plot_io (filplot, title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x, &
         dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
         gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, raux1, + 1)
  DEALLOCATE (raux1)
#else

  CALL plot_io (filplot, title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
        dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
        gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, raux, + 1)

#endif

  DEALLOCATE (raux)
  RETURN
END SUBROUTINE punch_plot
