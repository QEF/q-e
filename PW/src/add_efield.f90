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
SUBROUTINE add_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds an electric field to the local potential. The
  !   field is made artificially periodic by introducing a saw-tooth
  !   potential. The field is parallel to a reciprocal lattice vector bg, 
  !   according to the index edir.
  !
  !   if dipfield is false the electric field correction is added to the
  !   potential given as input (the bare local potential) only
  !   at the first call to this routine. In the following calls
  !   the routine exit.
  !
  !   if dipfield is true the dipole moment per unit surface is calculated
  !   and used to cancel the electric field due to periodic boundary
  !   conditions. This potential is added to the Hartree and xc potential
  !   in v_of_rho. NB: in this case the electric field contribution to the 
  !   band energy is subtracted by deband.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity
  
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr)! ef is added to this potential
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  REAL(DP),INTENT(IN)    :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  LOGICAL,INTENT(IN)     :: iflag ! set to true to force recalculation of field
  !
  ! local variables
  !
  INTEGER :: idx0, idx,  i, j, k
  INTEGER :: ir, na, ipol
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
  REAL(DP) :: tot_dipole, bmod

  LOGICAL :: first=.TRUE.
  SAVE first
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT.tefield) RETURN
  ! efield only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.
  IF ((.NOT.dipfield).AND.(.NOT.first) .AND..NOT. iflag) RETURN
  first=.FALSE.

  IF ((edir.lt.1).or.(edir.gt.3)) THEN
     CALL errore('add_efield',' wrong edir',1)
  ENDIF

  !---------------------
  !  Variable initialization
  !---------------------

  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)

  tot_dipole=0._dp
  e_dipole  =0._dp
  ion_dipole=0._dp
  
  !---------------------
  !  Calculate dipole
  !---------------------
  
  if (dipfield) then
  !
  ! dipole correction is active 
  !
     CALL compute_el_dip(emaxpos, eopreg, edir, rho, e_dipole)
     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)
    
     tot_dipole  = -e_dipole + ion_dipole
     CALL mp_bcast(tot_dipole, 0, intra_image_comm)
  !  
  !  E_{TOT} = -e^{2} \left( eamp - dip \right) dip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*(eamp-tot_dipole/2.d0)*tot_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2} \left( eamp - dip \right) z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *(eamp - tot_dipole) &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  else
  !
  ! dipole correction is not active
  !

     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)

  !  
  !  E_{TOT} = -e^{2} eamp * iondip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*eamp*ion_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2}  eamp z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *eamp &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  end if

  !
  !  Calculate potential and print values 
  !   
  
  length=(1._dp-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
  
  vamp=e2*(eamp-tot_dipole)*length

  IF (ionode) THEN
       !
       ! Output data
       !
       WRITE( stdout,*)
       WRITE( stdout,'(5x,"Adding external electric field":)')

       IF (dipfield) then
          WRITE( stdout,'(/5x,"Computed dipole along edir(",i1,") : ")' ) edir

          !
          !  If verbose prints also the different components
          !
          IF ( iverbosity > 0 ) THEN
              WRITE( stdout, '(8X,"Elec. dipole ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            e_dipole, (e_dipole*au_debye)
              WRITE( stdout, '(8X,"Ion. dipole  ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                          ion_dipole, (ion_dipole*au_debye)
          ENDIF

          WRITE( stdout, '(8X,"Dipole       ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
                                            (tot_dipole* (omega/fpi)),   &
                                            ((tot_dipole* (omega/fpi))*au_debye)  

          WRITE( stdout, '(8x,"Dipole field ", 1F15.4," Ry au, ")') &
                                             tot_dipole
          WRITE( stdout,*)

       ENDIF

       IF (abs(eamp)>0._dp) WRITE( stdout, &
          '(8x,"E field amplitude [Ha a.u.]: ", es11.4)') eamp 
        
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
  ! idx0 = starting index of real-space FFT arrays for this processor
  !
  idx0 = dfftp%nr1x*dfftp%nr2x*dfftp%ipp(me_bgrp+1)
  !
  DO ir = 1, dfftp%nr1x*dfftp%nr2x*dfftp%npl
     !
     ! ... three dimensional indexes
     !
     idx = idx0 + ir - 1
     k   = idx / (dfftp%nr1x*dfftp%nr2x)
     idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
     j   = idx / dfftp%nr1x
     idx = idx - dfftp%nr1x*j
     i   = idx

     ! ... do not include points outside the physical range

     IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
 
     if (edir.eq.1) sawarg = DBLE(i)/DBLE(dfftp%nr1)
     if (edir.eq.2) sawarg = DBLE(j)/DBLE(dfftp%nr2)
     if (edir.eq.3) sawarg = DBLE(k)/DBLE(dfftp%nr3)
     
     value = e2*(eamp - tot_dipole)*saw(emaxpos,eopreg,sawarg) * (alat/bmod)

     vpoten(ir) = vpoten(ir) + value

  END DO
  
  
  RETURN

END SUBROUTINE add_efield
! 
!------------------------------------------------------------------------------------------------
SUBROUTINE init_dipole_info (dipole_info, rho) 
!------------------------------------------------------------------------------------------------
! 
  USE kinds,           ONLY : DP
  USE constants,       ONLY : e2, fpi 
  USE qes_types_module,ONLY : dipoleOutput_type, scalarQuantity_type
  USE qes_libs_module, ONLY : qes_init_scalarQuantity, qes_reset_scalarQuantity
  USE extfield,        ONLY : edir, eamp, emaxpos, eopreg   
  USE fft_base,        ONLY : dfftp
  USE lsda_mod,        ONLY : nspin
  USE cell_base,       ONLY : alat, at, omega
  ! 
  IMPLICIT NONE  
  ! 
  TYPE ( dipoleOutput_type ), INTENT(OUT)  :: dipole_info
  REAL(DP),INTENT(IN)                      :: rho(dfftp%nnr,nspin)
  ! 
  REAL(DP)                                 :: ion_dipole, el_dipole, tot_dipole, length, vamp, fac
  TYPE ( scalarQuantity_type)              :: temp_qobj
  ! 
  CALL compute_ion_dip (emaxpos, eopreg, edir, ion_dipole)
  CALL compute_el_dip ( emaxpos, eopreg, edir, rho, el_dipole ) 
  tot_dipole = -el_dipole+ion_dipole
  ! 
  dipole_info%idir = edir  
  fac=omega/fpi
  CALL qes_init_scalarQuantity(dipole_info%ion_dipole,"ion_dipole" , units="Atomic Units", scalarQuantity= ion_dipole*fac)
  CALL qes_init_scalarQuantity(dipole_info%elec_dipole,"elec_dipole" , units="Atomic Units", scalarQuantity= el_dipole*fac)
  CALL qes_init_scalarQuantity(dipole_info%dipole,"dipole" , units="Atomic Units", scalarQuantity= tot_dipole*fac)
  CALL qes_init_scalarQuantity(dipole_info%dipoleField,"dipoleField" , units="Atomic Units", scalarQuantity= tot_dipole)
  ! 
  length=(1._DP-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
  vamp=e2*(eamp-tot_dipole)*length
  !
  CALL qes_init_scalarQuantity(dipole_info%potentialAmp,"potentialAmp" , units="Atomic Units", scalarQuantity= vamp)
  CALL qes_init_scalarQuantity(dipole_info%totalLength, "totalLength", units = "Bohr", scalarQuantity = length ) 
  
  END SUBROUTINE init_dipole_info
