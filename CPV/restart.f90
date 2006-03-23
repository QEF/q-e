!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE restart_file
!=----------------------------------------------------------------------------=!

  USE kinds, ONLY: DP

  IMPLICIT NONE

  PRIVATE

  SAVE

  REAL(DP) :: cclock
  EXTERNAL  :: cclock

  PUBLIC :: writefile, readfile

  INTERFACE readfile
    MODULE PROCEDURE readfile_cp, readfile_fpmd
  END INTERFACE

  INTERFACE writefile
    MODULE PROCEDURE writefile_cp, writefile_fpmd
  END INTERFACE
  
!=----------------------------------------------------------------------------=!
     CONTAINS
!=----------------------------------------------------------------------------=!

!-----------------------------------------------------------------------
      subroutine writefile_cp                                         &
     &     ( ndw,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm, &
     &       fion, tps, mat_z, occ_f, rho )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE ions_base,        ONLY: nsp, na, cdmi, taui
      USE cell_base,        ONLY: s_to_r
      USE cp_restart,       ONLY: cp_writefile
      USE electrons_base,   ONLY: nspin, nbnd, nbsp, iupdwn, nupdwn
      USE electrons_module, ONLY: ei
      USE io_files,         ONLY: scradir
      USE ensemble_dft,     ONLY: tens
!
      implicit none
      integer, INTENT(IN) :: ndw, nfi
      REAL(DP), INTENT(IN) :: h(3,3), hold(3,3)
      complex(8), INTENT(IN) :: c0(:,:), cm(:,:)
      REAL(DP), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
      REAL(DP), INTENT(IN) :: vels(:,:), velsm(:,:)
      REAL(DP), INTENT(IN) :: acc(:), lambda(:,:,:), lambdam(:,:,:)
      REAL(DP), INTENT(IN) :: xnhe0, xnhem, vnhe, ekincm
      REAL(DP), INTENT(IN) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer,      INTENT(in) :: nhpcl, nhpdim
      REAL(DP), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(in) :: ecut, ecutw, delt
      REAL(DP), INTENT(in) :: pmass(:)
      REAL(DP), INTENT(in) :: celldm(:)
      REAL(DP), INTENT(in) :: tps
      REAL(DP), INTENT(in) :: rho(:,:)
      integer, INTENT(in) :: ibrav
      REAL(DP), INTENT(in) :: occ_f(:)
      REAL(DP), INTENT(in) :: mat_z(:,:,:)

      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      integer :: nk = 1, ispin, i, ib
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 1.0d0
      REAL(DP), ALLOCATABLE :: occ_ ( :, :, : )
      REAL(DP) :: htm1(3,3), omega

!
! Do not write restart file if the unit number 
! is negative, this is used mainly for benchmarks
! and tests
!

      if ( ndw < 1 ) then
        return
      end if

      ht     = TRANSPOSE( h ) 
      htm    = TRANSPOSE( hold ) 
      htvel  = TRANSPOSE( velh ) 
      gvel   = 0.0d0
      

      ALLOCATE( occ_ ( nbnd, 1, nspin ) )
      occ_ = 0.0d0
      do ispin = 1, nspin
        do i = iupdwn ( ispin ), iupdwn ( ispin ) - 1 + nupdwn ( ispin )
          occ_ ( i - iupdwn ( ispin ) + 1, 1, ispin ) = occ_f( i ) 
        end do
      end do

      IF( tens ) THEN
        CALL cp_writefile( ndw, scradir, .TRUE., nfi, tps, acc, nk, xk, wk, &
          ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus, &
          vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_ , &
          occ_ , lambda, lambdam, xnhe0, xnhem, vnhe, ekincm, ei, &
          rho, c02 = c0, cm2 = cm, mat_z = mat_z  )
      ELSE
        CALL cp_writefile( ndw, scradir, .TRUE., nfi, tps, acc, nk, xk, wk, &
          ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus, &
          vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_ , &
          occ_ , lambda, lambdam, xnhe0, xnhem, vnhe, ekincm, ei, &
          rho, c02 = c0, cm2 = cm  )
      END IF

      DEALLOCATE( occ_ )

      return
      end subroutine writefile_cp

!-----------------------------------------------------------------------
      subroutine readfile_cp                                        &
     &     ( flag, ndr,h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,&
     &       fion, tps, mat_z, occ_f )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE io_files,       ONLY : scradir
      USE electrons_base, ONLY : nbnd, nbsp, nspin, nupdwn, iupdwn
      USE gvecw,          ONLY : ngw, ngwt
      USE ions_base,      ONLY : nsp, na, cdmi, taui
      USE cp_restart,     ONLY : cp_readfile, cp_read_cell, cp_read_wfc
      USE ensemble_dft,   ONLY : tens
      USE autopilot,      ONLY : event_step, event_index, max_event_step
      USE autopilot,      ONLY : employ_rules
!
      implicit none
      integer :: ndr, nfi, flag
      REAL(DP) :: h(3,3), hold(3,3)
      complex(8) :: c0(:,:), cm(:,:)
      REAL(DP) :: tausm(:,:),taus(:,:), fion(:,:)
      REAL(DP) :: vels(:,:), velsm(:,:)
      REAL(DP) :: acc(:),lambda(:,:,:), lambdam(:,:,:)
      REAL(DP) :: xnhe0,xnhem,vnhe
      REAL(DP) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer, INTENT(inout) :: nhpcl,nhpdim
      REAL(DP) :: ekincm
      REAL(DP) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(in) :: ecut, ecutw, delt
      REAL(DP), INTENT(in) :: pmass(:)
      REAL(DP), INTENT(in) :: celldm(6)
      integer, INTENT(in) :: ibrav
      REAL(DP), INTENT(OUT) :: tps
      REAL(DP), INTENT(INOUT) :: mat_z(:,:,:), occ_f(:)
      !
      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      integer :: nk = 1, ispin, i, ib
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 1.0d0
      REAL(DP), ALLOCATABLE :: occ_ ( :, :, : )
      REAL(DP) :: htm1(3,3), b1(3) , b2(3), b3(3), omega
        
      LOGICAL::lopen

      IF( flag == -1 ) THEN
        CALL cp_read_cell( ndr, scradir, .TRUE., ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh )
        h     = TRANSPOSE( ht )
        hold  = TRANSPOSE( htm )
        velh  = TRANSPOSE( htvel )
        RETURN
      ELSE IF ( flag == 0 ) THEN
        DO ispin = 1, nspin
          CALL cp_read_wfc( ndr, scradir, 1, 1, ispin, nspin, c2 = cm(:,:), tag = 'm' )
        END DO
        RETURN
      END IF

      ALLOCATE( occ_ ( nbnd, 1, nspin ) )

      IF( tens ) THEN
         CALL cp_readfile( ndr, scradir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c02 = c0, cm2 = cm, mat_z = mat_z )
      ELSE
         CALL cp_readfile( ndr, scradir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c02 = c0, cm2 = cm )
      END IF

      ! AutoPilot (Dynamic Rules) Implementation
      event_index = 1

      do while (event_step(event_index) <= nfi)
         ! Assuming that the remaining dynamic parm values are set correctly by reading 
         ! the the restart file.
         ! if this is not true, employ rules as events are updated right here as:
         call employ_rules()
         event_index = event_index + 1
         if(  event_index > max_event_step ) then
            call errore(' readfile ','maximum events exceeded for dynamic rules', 1)
         endif
      enddo
      
      do ispin = 1, nspin
        do i = iupdwn ( ispin ), iupdwn ( ispin ) - 1 + nupdwn ( ispin )
          occ_f( i ) = occ_ ( i - iupdwn ( ispin ) + 1, 1, ispin )
        end do
      end do
      DEALLOCATE( occ_ )

      h     = TRANSPOSE( ht )
      hold  = TRANSPOSE( htm )
      velh  = TRANSPOSE( htvel )

      return
      end subroutine readfile_cp


!=----------------------------------------------------------------------------=!


   SUBROUTINE writefile_fpmd( nfi, trutime, c0, cm, cdesc, occ, &
     atoms_0, atoms_m, acc, taui, cdmi, &
     ht_m, ht_0, rho, vpot)
                                                                        
        USE cell_module, only: boxdimensions, r_to_s
        USE brillouin, only: kpoints, kp
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndw, gamma_only
        USE control_flags, ONLY: twfcollect, force_pairing
        USE atoms_type_module, ONLY: atoms_type
        USE io_global, ONLY: ionode, ionode_id
        USE io_global, ONLY: stdout
        USE electrons_nose, ONLY: xnhe0, xnhem, vnhe
        USE electrons_base, ONLY: nbsp, nspin, nudx
        USE cell_nose, ONLY: xnhh0, xnhhm, vnhh
        USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, nhpcl, nhpdim
        USE cp_restart, ONLY: cp_writefile
        USE electrons_module, ONLY: ei
        USE io_files, ONLY: scradir
        USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x

        IMPLICIT NONE 
 
        INTEGER, INTENT(IN) :: nfi
        COMPLEX(DP), INTENT(IN) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(DP), INTENT(IN) :: occ(:,:,:)
        TYPE (boxdimensions), INTENT(IN) :: ht_m, ht_0
        TYPE (atoms_type), INTENT(IN) :: atoms_0, atoms_m
        REAL(DP), INTENT(IN) :: rho(:,:)
        TYPE (wave_descriptor) :: cdesc
        REAL(DP), INTENT(INOUT) :: vpot(:,:)
                                                                        
        REAL(DP), INTENT(IN) :: taui(:,:)
        REAL(DP), INTENT(IN) :: acc(:), cdmi(:) 
        REAL(DP), INTENT(IN) :: trutime

        REAL(DP), ALLOCATABLE :: lambda(:,:,:)
        REAL(DP) :: ekincm
        INTEGER   :: i, j, k, iss, ir
             
        IF( ndw < 1 ) RETURN
        !
        !   this is used for benchmarking and debug
        !   if ndw < 1 Do not save wave functions and other system
        !   properties on the writefile subroutine

        ALLOCATE( lambda(nudx,nudx,nspin) )
        lambda  = 0.0d0
        ekincm = 0.0d0
        !
        !
        CALL cp_writefile( ndw, scradir, .TRUE., nfi, trutime, acc, kp%nkpt, kp%xk, kp%weight, &
          ht_0%a, ht_m%a, ht_0%hvel, ht_0%gvel, xnhh0, xnhhm, vnhh, taui, cdmi, &
          atoms_0%taus, atoms_0%vels, atoms_m%taus, atoms_m%vels, atoms_0%for, vnhp, &
          xnhp0, xnhpm, nhpcl, nhpdim, occ, occ, lambda, lambda,  &
          xnhe0, xnhem, vnhe, ekincm, ei, rho, c04 = c0, cm4 = cm )

        DEALLOCATE( lambda )

     RETURN 
   END SUBROUTINE writefile_fpmd



!=----------------------------------------------------------------------------=!

        SUBROUTINE readfile_fpmd( nfi, trutime, &
          c0, cm, cdesc, occ, atoms_0, atoms_m, acc, taui, cdmi, &
          ht_m, ht_0, rho, vpot )
                                                                        
        use electrons_base, only: nbsp, nspin, nudx
        USE cell_module, only: boxdimensions, cell_init, r_to_s, s_to_r
        USE brillouin, only: kpoints, kp
        use parameters, only: npkx, nsx
        USE mp_global, ONLY: mpime, nproc, intra_image_comm, root
        USE mp_wave, ONLY: mergewf
        USE wave_types, ONLY: wave_descriptor
        USE control_flags, ONLY: ndr, tbeg, gamma_only
        USE atoms_type_module, ONLY: atoms_type
        USE io_global, ONLY: ionode
        USE io_global, ONLY: stdout
        USE gvecw, ONLY: ecutwfc => ecutw
        USE gvecp, ONLY: ecutrho => ecutp
        USE ions_base, ONLY: nat, nsp, na
        USE control_flags, ONLY: twfcollect, force_pairing
        USE grid_dimensions, ONLY: nr1, nr2, nr3
        USE electrons_nose, ONLY: xnhe0, xnhem, vnhe
        USE cell_nose, ONLY: xnhh0, xnhhm, vnhh
        USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, nhpcl, nhpdim
        USE cp_restart, ONLY: cp_readfile
        USE io_files, ONLY: scradir
 
        IMPLICIT NONE 
 
        INTEGER, INTENT(OUT) :: nfi
        COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:,:), cm(:,:,:,:) 
        REAL(DP), INTENT(INOUT) :: occ(:,:,:)
        TYPE (boxdimensions), INTENT(INOUT) :: ht_m, ht_0
        TYPE (atoms_type), INTENT(INOUT) :: atoms_0, atoms_m
        REAL(DP), INTENT(INOUT) :: rho(:,:)
        TYPE (wave_descriptor) :: cdesc
        REAL(DP), INTENT(INOUT) :: vpot(:,:)
                                                                        
        REAL(DP), INTENT(OUT) :: taui(:,:)
        REAL(DP), INTENT(OUT) :: acc(:), cdmi(:) 
        REAL(DP), INTENT(OUT) :: trutime

        REAL(DP), ALLOCATABLE :: lambda_ ( : , : , : )
        REAL(DP) :: ekincm
        REAL(DP) :: hp0_ (3,3)
        REAL(DP) :: hm1_ (3,3)
        REAL(DP) :: gvel_ (3,3)
        REAL(DP) :: hvel_ (3,3)
        REAL(DP) :: b1(3), b2(3), b3(3)
        LOGICAL :: tens = .FALSE.

        ALLOCATE( lambda_( nudx , nudx, nspin ) )
        lambda_  = 0.0d0

        CALL cp_readfile( ndr, scradir, .TRUE., nfi, trutime, acc, kp%nkpt, kp%xk, kp%weight, &
          hp0_ , hm1_ , hvel_ , gvel_ , xnhh0, xnhhm, vnhh, taui, cdmi, &
          atoms_0%taus, atoms_0%vels, atoms_m%taus, atoms_m%vels, atoms_0%for, vnhp, &
          xnhp0, xnhpm, nhpcl, nhpdim, occ, occ, lambda_ , lambda_ , b1, b2,   &
          b3, xnhe0, xnhem, vnhe, ekincm, c04 = c0, cm4 = cm )

        DEALLOCATE( lambda_ )

        IF( .NOT. tbeg ) THEN
          CALL cell_init( ht_0, hp0_ )
          CALL cell_init( ht_m, hm1_ )
          ht_0%hvel = hvel_  !  set cell velocity
          ht_0%gvel = gvel_  !  set cell velocity
        END IF

        RETURN 
        END SUBROUTINE readfile_fpmd

!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
     END MODULE restart_file
!=----------------------------------------------------------------------------=!
