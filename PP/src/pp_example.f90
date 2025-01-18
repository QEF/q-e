!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM post_processing_example
  !-----------------------------------------------------------------------
  !
  ! Sample code, showing how to read QE data and re-use QE variables:
  ! 1. reads the data directory of QE, then
  ! 2. fills the hamiltonian matrix and diagonalizes it for each k-point
  !    (conventional, not iterative diagonalization)
  !
  ! Input: namelist &inputpp [outdir=...] [prefix=...] / as in QE input
  ! (default values as in QE).
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : tmp_dir, prefix
  USE mp_global,  ONLY : mp_startup
  USE mp_images,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: needwf = .true.
  INTEGER :: ios
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / inputpp / outdir, prefix
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL environment_start ( 'PPEXAMPLE' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm)
  IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  !
  !   Read xml file, allocate and initialize general variables
  !
  CALL read_file_new ( needwf )
  !
  CALL re_diagonalize ( )
  !
  CALL environment_end ( 'PPEXAMPLE' )
  !
  CALL stop_pp()
  !
END PROGRAM post_processing_example
!
!-----------------------------------------------------------------------
SUBROUTINE re_diagonalize ( )
  !-----------------------------------------------------------------------
  !
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE becmod,     ONLY : becp, calbec, allocate_bec_type
  USE fft_base,   ONLY : dfftp
  USE klist,      ONLY : xk, nks, nkstot, igk_k, ngk
  USE lsda_mod,   ONLY : nspin, isk, current_spin
  USE io_files,   ONLY : restart_dir
  USE scf,        ONLY : vrs, vltot, v, kedtau
  USE gvecs,      ONLY : doublegrid
  USE uspp,       ONLY : nkb, vkb
  USE uspp_init,  ONLY : init_us_2
  USE wvfct,      ONLY : npwx, nbnd, current_k
  USE mp_bands,   ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
  USE wavefunctions,  ONLY : evc
  USE pw_restart_new, ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.h'
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  INTEGER :: ik, npw
  !
  ALLOCATE( aux(npwx, nbnd ) )
  ALLOCATE( hc( nbnd, nbnd) )    
  ALLOCATE( sc( nbnd, nbnd) )    
  ALLOCATE( vc( nbnd, nbnd) )    
  ALLOCATE( en( nbnd ) )
  CALL allocate_bec_type (nkb, nbnd, becp )
  CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)

  DO ik = 1, nkstot
     !
     CALL read_collected_wfc ( restart_dir() , ik, evc )
     !
     npw = ngk(ik)
     current_k=ik
     current_spin  = isk(ik)
     !
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     CALL calbec ( npw, vkb, evc, becp)
     CALL g2_kin (ik)
     !
     CALL h_psi( npwx, npw, nbnd, evc, aux )
     CALL calbec ( npw, evc, aux, hc )
     CALL s_psi( npwx, npw, nbnd, evc, aux )
     CALL calbec ( npw, evc, aux, sc )
     !
     CALL diaghg( nbnd, nbnd, hc, sc, nbnd, en, vc, me_bgrp, &
          root_bgrp, intra_bgrp_comm )
     !
     print '(/,12x,"k =",3f7.4," (",i6," PWs)  bands (eV):",/)', xk(:,ik),npw
     print '(8f9.4)', en(:)*rytoev
     !
  END DO
  !
END SUBROUTINE re_diagonalize
