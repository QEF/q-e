!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by J. Tobik
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
!          25/06/2009 (Riccardo Sabatini)
!               reformulation using a unique saw(x) function (included in 
!               cell_base) in all e-field related routines and inclusion of 
!               a macroscopic electronic dipole contribution in the mixing 
!               scheme. 
!

!
!--------------------------------------------------------------------------
SUBROUTINE add_efield( vpoten, etotefield, rho, iflag )
  !--------------------------------------------------------------------------
  !! This routine adds an electric field to the local potential. The
  !! field is made artificially periodic by introducing a saw-tooth
  !! potential. The field is parallel to a reciprocal lattice vector bg, 
  !! according to the index edir.
  !
  !! * If \(\textit{dipfield}\) is false the electric field correction is 
  !!   added to the potential given as input (the bare local potential) 
  !!   only at the first call to this routine. In the following calls the
  !!   routine exit.
  !! * If \(\textit{dipfield}\) is true the dipole moment per unit surface
  !!   is calculated and used to cancel the electric field due to periodic
  !!   boundary conditions. This potential is added to the Hartree and 
  !!   xc potential in v_of_rho. NB: in this case the electric field 
  !!   contribution to the band energy is subtracted by deband.
  !
  USE kinds,         ONLY: DP
  USE constants,     ONLY: fpi, eps8, e2, au_debye
  USE ions_base,     ONLY: nat, ityp, zv
  USE cell_base,     ONLY: alat, at, omega, bg
  USE extfield,      ONLY: tefield, dipfield, edir, eamp, emaxpos, saw, &
                           eopreg, forcefield, el_dipole, ion_dipole, tot_dipole
  USE io_global,     ONLY: stdout,ionode
  USE control_flags, ONLY: mixing_beta, lforce => tprnfor
  USE lsda_mod,      ONLY: nspin
  USE mp_images,     ONLY: intra_image_comm
  USE mp_bands,      ONLY: me_bgrp
  USE fft_base,      ONLY: dfftp
  USE fft_types,     ONLY: fft_index_to_3d
  USE mp,            ONLY: mp_bcast, mp_sum
  USE control_flags, ONLY: iverbosity
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: vpoten(dfftp%nnr)
  !! ef is added to this potential
  REAL(DP), INTENT(INOUT) :: etotefield
  !! contribution to etot due to ef
  REAL(DP), INTENT(IN) :: rho(dfftp%nnr)
  !! the density whose dipole is computed
  LOGICAL,INTENT(IN) :: iflag
  !! set to true to force recalculation of field
  !
  ! ... local variables
  !
  INTEGER :: i, j, k
  INTEGER :: ir, na, ipol
  REAL(DP) :: length, vamp, value, sawarg, bmod
  !
  LOGICAL :: offrange, first=.TRUE.
  SAVE first
  !
  !---------------------
  !  Execution control
  !---------------------
  !
  IF (.NOT.tefield) RETURN
  ! efield only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.
  IF ((.NOT.dipfield).AND.(.NOT.first) .AND..NOT. iflag) RETURN
  first=.FALSE.
  !
  IF ((edir<1).OR.(edir>3)) THEN
     CALL errore( 'add_efield', ' wrong edir', 1 )
  ENDIF
  !
  !---------------------
  !  Variable initialization
  !---------------------
  !
  bmod = SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)
  !
  tot_dipole = 0._dp
  el_dipole  = 0._dp
  ion_dipole = 0._dp
  !
  !---------------------
  !  Calculate dipole
  !---------------------
  ! 
  IF (dipfield) THEN
     !
     ! dipole correction is active 
     !
     CALL compute_el_dip( emaxpos, eopreg, edir, rho, el_dipole )
     CALL compute_ion_dip( emaxpos, eopreg, edir, ion_dipole )
     !
     tot_dipole = -el_dipole + ion_dipole
     CALL mp_bcast( tot_dipole, 0, intra_image_comm )
     !  
     !  E_{TOT} = -e^{2} \left( eamp - dip \right) dip \frac{\Omega}{4\pi} 
     !
     etotefield = -e2 * (eamp-tot_dipole/2.d0) * tot_dipole * omega/fpi 
     !
     !---------------------
     !  Define forcefield
     !  
     !  F_{s} = e^{2} \left( eamp - dip \right) z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
     !---------------------
     !
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *(eamp - tot_dipole) &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF
     !
  ELSE
     !
     ! dipole correction is not active
     !
     CALL compute_ion_dip( emaxpos, eopreg, edir, ion_dipole )
     !  
     !  E_{TOT} = -e^{2} eamp * iondip \frac{\Omega}{4\pi} 
     !
     etotefield=-e2*eamp*ion_dipole*omega/fpi 
     !---------------------
     !  Define forcefield
     !  
     !  F_{s} = e^{2}  eamp z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
     !---------------------
     !
     IF (lforce) THEN
        DO na = 1, nat
           DO ipol = 1, 3
              forcefield(ipol,na) = e2 * eamp * &
                                    zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF
     !
  ENDIF
  !
  !  Calculate potential and print values 
  !   
  length = (1._dp-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
  !
  vamp = e2*(eamp-tot_dipole)*length
  !
  IF (ionode) THEN
       !
       ! Output data
       !
       WRITE( stdout,*)
       WRITE( stdout,'(5x,"Adding external electric field":)')
       !
       IF (dipfield) THEN
          WRITE( stdout,'(/5x,"Computed dipole along edir(",i1,") : ")' ) edir
          !
          !  If verbose prints also the different components
          !
          IF ( iverbosity > 0 ) THEN
              WRITE( stdout, '(8X,"Elec. dipole ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            el_dipole, (el_dipole*au_debye)
              WRITE( stdout, '(8X,"Ion. dipole  ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                          ion_dipole, (ion_dipole*au_debye)
          ENDIF
          !
          WRITE( stdout, '(8X,"Dipole       ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            (tot_dipole* (omega/fpi)),   &
                                            ((tot_dipole* (omega/fpi))*au_debye)  
          !
          WRITE( stdout, '(8x,"Dipole field ", 1F15.4," Ry au, ")') &
                                             tot_dipole
          WRITE( stdout,*)
          !
       ENDIF

       IF (ABS(eamp)>0._dp) WRITE( stdout, &
          '(8x,"E field amplitude [Ha a.u.]: ", es11.4)') eamp 
       !
       WRITE( stdout,'(8x,"Potential amp.   ", f11.4," Ry")') vamp 
       WRITE( stdout,'(8x,"Total length     ", f11.4," bohr")') length
       WRITE( stdout,*)     
  ENDIF
  !
  !------------------------------
  !  Add potential
  !  
  !  V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v}
  !          Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod} 
  !          
  !---------------------
  !
  ! Loop in the charge array
  !
  !
  DO ir = 1, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
     !
     ! ... three dimensional indexes
     !
     CALL fft_index_to_3d (ir, dfftp, i,j,k, offrange)
     IF ( offrange ) CYCLE
     !
     IF (edir==1) sawarg = DBLE(i)/DBLE(dfftp%nr1)
     IF (edir==2) sawarg = DBLE(j)/DBLE(dfftp%nr2)
     IF (edir==3) sawarg = DBLE(k)/DBLE(dfftp%nr3)
     !
     value = e2*(eamp - tot_dipole)*saw(emaxpos,eopreg,sawarg) * (alat/bmod)
     !
     vpoten(ir) = vpoten(ir) + value
     !
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE add_efield
! 

