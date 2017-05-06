!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by T. Brumme with add_efield.f90 as template
!
!
!--------------------------------------------------------------------------
SUBROUTINE add_monofield(vpoten,etotmonofield,linear,quadratic)
  !--------------------------------------------------------------------------
  !
  ! This routine adds an electric field of a monopole-plane to the local potential.
  ! (in the system setup - setlocal.f90)
  ! Furthermore, a harmonic potential is added as background charge.
  ! USE ONLY WITH ELECTRIC FIELD IF DIPOLE CORRECTION IS ACITVE
  !
  ! see PRB 89, 245406 (2014)
  !
  ! compare also with add_efield as this file here is more or less a copy of it
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye, tpi
  USE ions_base,     ONLY : nat, ityp, zv, tau
  USE cell_base,     ONLY : alat, at, omega, bg
  USE extfield,      ONLY : zmon, monopole, dipfield, forcemono, mopopla, &
                            relaxz, block, block_height, block_1, block_2, eopreg
  USE klist,         ONLY : nelec
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
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr)! field is added to this potential
  REAL(DP),INTENT(INOUT) :: etotmonofield    ! contribution to etot due to field
  LOGICAL,INTENT(IN)     :: linear    ! set to true to calculate the linear part
  LOGICAL,INTENT(IN)     :: quadratic ! set to true to calculate the quadratic part
  !
  ! local variables
  !
  INTEGER :: idx0, idx,  i, j, k
  INTEGER :: ir, na, ipol
  REAL(DP) :: value, monoarg, ion_dipole, monoamp, block_size
  REAL(DP) :: bmod, ionic_charge, area, sgn1, zvia, tvectb
  REAL(DP), ALLOCATABLE :: xau(:,:)

  LOGICAL :: first=.TRUE.

  SAVE first
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT.monopole) RETURN

  ! if we are not in the first step, and we do not want to calculate the linear or quadratic part
  IF ( (.NOT.first) .AND. (.NOT.linear) .AND. (.NOT.quadratic) ) RETURN ! why are we here?
  ! btw check if 'first' is needed at all

  ! check if we have a proper unit cell with the third vector being orthogonal to the first two
  ! and being along z, i.e.
  !  CELL_PARAMETERS
  !  x1 y1 0
  !  x2 y2 0
  !  0  0  z
  ! there might be an easier way of checking this instead of checking all components
  ! but I'm too lazy atm and it is only checked once in the beginning of an scf
  IF (.not.(at(3,1)==0.0 .and. at(3,2)==0 .and. at(1,3)==0 .and. at(2,3)==0)) &
     CALL errore( 'add_monofield', '3. lattice vector has to be orthogonal to others and along z' , 1 )

  !---------------------
  !  Variable initialization
  !---------------------

  bmod=SQRT(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
  ion_dipole=0.0
  ionic_charge = SUM( zv(ityp(1:nat)) )
  area = ABS((at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2)
  monoamp = (-(nelec-ionic_charge)/area*tpi)

  IF (block) THEN
     block_size=block_2-block_1
  ENDIF
  
  !-----------------------------
  !  Calculate energy and forces
  !-----------------------------
  
  !energy of the ions in the field of the monopole
  DO na = 1, nat
     zvia = zv(ityp(na))
     tvectb = tau(1,na)*bg(1,3)+tau(2,na)*bg(2,3)+tau(3,na)*bg(3,3)
     ion_dipole = ion_dipole+zvia*(mopopla(zmon,tvectb,.true.)+1/6)*(alat/bmod)*(fpi/omega) !the linear part
     ion_dipole = ion_dipole+zvia*mopopla(zmon,tvectb,.false.)*(alat/bmod)*(fpi/omega)      !the quadratic part
  END DO
  etotmonofield = -e2 * monoamp * ion_dipole * omega/fpi

  !energy of the monopole in its own potential
  etotmonofield = etotmonofield - e2 * (nelec-ionic_charge) * monoamp * (alat/bmod) / 6.0

  !---------------------
  !  Define forcefield
  !---------------------
    
  IF (lforce) THEN

     forcemono=0
     ALLOCATE (xau(3,nat))
     !     Compute the coordinates of each atom in the basis of the
     !     direct lattice vectors
     DO na = 1, nat
        DO ipol = 1, 3
           xau(ipol,na) = bg(1,ipol)*tau(1,na) + bg(2,ipol)*tau(2,na) + bg(3,ipol)*tau(3,na)
        ENDDO
     ENDDO
     DO na=1,nat
        sgn1=0
        sgn1=xau(3,na)
        ! atom position within cell?
        DO
          IF (sgn1>1) sgn1=sgn1-1
          IF (sgn1<0) sgn1=sgn1+1
          IF (sgn1<=1.and.sgn1>=0) EXIT
        ENDDO
        sgn1=sgn1-zmon

        ! periodicity of the potential, thus the force
        IF (sgn1<=-0.5) sgn1=sgn1+1
        IF (sgn1>0.5) sgn1=sgn1-1

        ! two parts here, first for the linear part of the potential, ie the monopole
        ! and second part is related to the harmonic potential, ie the background
        DO ipol=1,3
          forcemono(ipol,na) = -zv(ityp(na)) * monoamp * e2 * (SIGN(1._dp,sgn1)-2.0*sgn1) * bg(ipol,3)/bmod
        ENDDO

     ENDDO
     !
     !   deallocate work space
     !
     DEALLOCATE(xau)

  ENDIF


  IF (ionode) THEN
       !
       ! Output data
       !
       WRITE( stdout,*)
       WRITE( stdout,'(5x,"Adding monopole plane to compensate the charge of the system - add_monofield")')
       WRITE( stdout,'(5x,"see PRB 89, 245406 (2014), works only for 2D systems perpendicular to z")')
       WRITE( stdout,'(5x,"                           i.e. 3rd lattice vector is (0,0,c)")')
       WRITE( stdout,*)

       WRITE( stdout, '(8x,"prefactor of the potential in [Ha a.u.]: ", f12.6)') monoamp
       WRITE( stdout, '(8x,"   position of the monopole within cell:    ", f8.5)') zmon
       WRITE( stdout, '(8x," ion-monopole + monopole-monopole contribution to the total energy: ", f12.6)') etotmonofield
       WRITE( stdout, '(8x,"         monopole-monopole contribution: ", f12.6)') &
              (- (nelec-ionic_charge) * monoamp * (alat/bmod) / 6.0)

       IF ((lforce).AND.(relaxz)) THEN
          WRITE( stdout,*)
          WRITE( stdout,'(8x,"Allow relaxation in z-direction (i.e. disabled control for total force = 0) ")')
       ENDIF

       IF (block) THEN
          WRITE( stdout,*)
          IF (dipfield) THEN
             WRITE( stdout,'(8x,"Adding potential to prevent charge spilling into region of the monopole")')
             WRITE( stdout,'(8x,"Potential is linearly increased and decreased within the first/last eopreg of block_size")')
          ELSE
             WRITE( stdout,'(8x,"Adding potential to prevent interaction between lower and upper part of unit cell")')
             WRITE( stdout,'(8x,"Potential is linearly increased and decreased within the first/last 10% of block_size")')
          ENDIF
          WRITE( stdout,'(8x,"block_size = ", f8.5," in units of unit cell length")') block_size
          WRITE( stdout,'(8x,"block_height = ", f8.5," in Ry")') block_height
       ENDIF

       WRITE( stdout,*)
  ENDIF


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

     ! ... do not include points outside the physical range

     IF ( k >= dfftp%nr3 ) CYCLE
 
     monoarg = DBLE(k)/DBLE(dfftp%nr3)

     ! PRB 89, 245406 (2014)
     !
     ! V(z') = - n e^2 / (2 e0) (-z' + z'^2 / L)
     !           |          ||  \------||------/
     !           |        1/4pi      mopopla
     !           |                  L = alat/bmod [in bohr]
     !           ||                (z' is in units of L => (-z' + z'^2 / L)=(-z' + z'^2)*L
     ! total charge = n * e * area
     ! n > 0 electron doping, n < 0 hole doping
     ! => total charge = (nelec - ionic_charge) * e
     !
     ! V(z') = -(nelec-ionic_charge)/area * e2 / (2 * 1/fpi) * mopopla * L
     !       = -(nelec-ionic_charge)/area * e2 * tpi * mopopla * L
     !       = monoamp * e2 * mopopla * (alat/bmod)

     ! for the linear part we have different case
     !
     ! 1. no dipole correction, i.e. dipfield .ne. true
     ! 1.1 with potential barrier to block tunneling in z-direction or avoid spilling
     ! 1.2 without, which is the same as with but in the part where Vb = 0
     ! 2. dipole correction active
     ! 2.1 with potential barrier to avoid the relaxation of the ions towards the monopole
     ! 2.2 without (no relaxation), which is the same as with but in the part where Vb = 0
     !
     IF (linear) THEN

        ! case 1.1
        IF ((block).AND.(monoarg>=block_1).AND.(monoarg<=block_2).AND.(.NOT.(dipfield))) THEN
           IF (monoarg-zmon<=-block_size/2.0*0.9) THEN ! smooth increase within the first 10%
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) &
                      + block_height * ((monoarg-zmon)+block_size/2.0)/(0.1*block_size/2.0)
           ELSEIF (monoarg-zmon>=block_size/2.0*0.9) THEN ! smooth decrease
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) &
                      + block_height * ((monoarg-zmon)-block_size/2.0)/(-0.1*block_size/2.0)
           ELSE ! block
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) + block_height
           ENDIF

        ! case 2.1 - blocking used to model dielectric?
        ELSEIF ((block).AND.(monoarg>=block_1).AND.(monoarg<=block_2).AND.(dipfield)) THEN

           IF (monoarg<=(block_1+eopreg)) THEN ! smooth increase within the first eopreg
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) &
                      + (monoarg-block_1)/(eopreg) * block_height

           ELSEIF (monoarg>=(block_2-eopreg)) THEN ! smooth decrease within the last eopreg
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) &
                      + (block_2-monoarg )/(eopreg) * block_height
           ELSE
              value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod) + block_height
           ENDIF

        ! linear part of the potential (case 1.2 and 2.2)
        ELSE
           value = monoamp * e2 * ( mopopla(zmon,monoarg,.true.) + 1.0/6.0 ) * (alat/bmod)

        ENDIF
        vpoten(ir) = vpoten(ir) + value
     ENDIF

     IF (quadratic) THEN   ! quadratic part of the potential, i.e. the background charge
        value = monoamp * e2 * mopopla(zmon,monoarg,.false.) * (alat/bmod)
        vpoten(ir) = vpoten(ir) + value
     ENDIF

  END DO

  first=.FALSE.

  RETURN

END SUBROUTINE add_monofield
