!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_fft
  !-----------------------------------------------------------------------
  !! This routine allocates memory for FFT-related arrays.  
  !! IMPORTANT:  routine "data_structure" must be called before it in
  !! order to set the proper dimensions and grid distribution across
  !! processors these dimensions.
  !
  USE io_global,        ONLY : stdout
  USE gvect,            ONLY : ngm, g, gg, mill, igtongl
  USE gvecs,            ONLY : ngms
  USE fft_base,         ONLY : dfftp, dffts
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : rho, v, vnew, vltot, vrs, rho_core, rhog_core, &
                               kedtau, create_scf_type
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : pointlist, factlist, report, noncolin, npol
  USE wavefunctions,    ONLY : psic, psic_nc
  USE xc_lib,           ONLY : xclib_dft_is
  !
  USE scf_gpum,  ONLY : using_vrs
  !
  USE wavefunctions_gpum, ONLY : using_psic, using_psic_nc, &
                                        using_psic_d, using_psic_nc_d
  !
  IMPLICIT NONE
  !
  ! ... First a bunch of checks
  !
  IF (dfftp%nnr < ngm) THEN
     WRITE( stdout, '(/,4x," nr1=",i4," nr2= ", i4, " nr3=",i4, &
          &" nrxx = ",i8," ngm=",i8)') dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, ngm
     CALL errore( 'allocate_fft', 'the nr"s are too small!', 1 )
  !
  ENDIF
  !
  IF (dffts%nnr < ngms) THEN
     WRITE( stdout, '(/,4x," nr1s=",i4," nr2s= ", i4, " nr3s=",i4, &
          &" nrxxs = ",i8," ngms=",i8)') dffts%nr1, dffts%nr2, dffts%nr3, dffts%nnr, ngms
     CALL errore( 'allocate_fft', 'the nrs"s are too small!', 1 )
  ENDIF
  !
  IF (ngm  <= 0)      CALL errore( 'allocate_fft', 'wrong ngm' , 1 )
  IF (ngms <= 0)      CALL errore( 'allocate_fft', 'wrong ngms', 1 )
  IF (dfftp%nnr <= 0) CALL errore( 'allocate_fft', 'wrong nnr',  1 )
  IF (dffts%nnr<= 0)  CALL errore( 'allocate_fft', 'wrong smooth nnr', 1 )
  IF (nspin<= 0)      CALL errore( 'allocate_fft', 'wrong nspin', 1 )
  !
  ! ... Allocate memory for all kind of stuff.
  !
  CALL create_scf_type( rho )
  CALL create_scf_type( v,    do_not_allocate_becsum = .TRUE. )
  CALL create_scf_type( vnew, do_not_allocate_becsum = .TRUE. )
  !
  ALLOCATE( vltot(dfftp%nnr) )
  ALLOCATE( rho_core(dfftp%nnr) )
  IF ( xclib_dft_is('meta') ) THEN
     ALLOCATE( kedtau(dffts%nnr,nspin) )
  ELSE
     ALLOCATE( kedtau(1,nspin) )
  ENDIF
  ALLOCATE( rhog_core(ngm)  )
  ALLOCATE( psic(dfftp%nnr) )
  ALLOCATE( vrs(dfftp%nnr,nspin) )
#if defined(__CUDA)
  CALL using_vrs(2)
  CALL using_psic(2); CALL using_psic_d(0)
#endif
  !
  IF (noncolin) ALLOCATE( psic_nc(dfftp%nnr,npol) )
#if defined(__CUDA)
  IF (noncolin) THEN
     CALL using_psic_nc(2)
     CALL using_psic_nc_d(0)
  END IF
#endif
  !
  IF ( report /= 0 ) THEN
     !
     ! ... In order to print out local quantities, integrated around the atoms,
     ! we need the following variables
     !
     ALLOCATE( pointlist(dfftp%nnr) )
     ALLOCATE( factlist(dfftp%nnr)  )
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE allocate_fft
