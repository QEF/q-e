!
! Copyright (C) 2003-2024 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
PROGRAM wannier2pw
  !---------------------------------------------------------------------------
  !
  ! This is the post-processing program that prepares the Wannier functions 
  ! from Wannier90 for calculations using PW and other codes.
  ! This program is written using the same logic and style as 
  ! PP/src/pw2wannier90.f90.
  !
  ! Written by I. Timrov (November 2024)
  !
  USE io_global,             ONLY : stdout, ionode, ionode_id
  USE io_files,              ONLY : prefix, tmp_dir
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : mp_startup
  USE mp_world,              ONLY : world_comm
  USE mp_bands,              ONLY : nbgrp
  USE mp_pools,              ONLY : npool
  USE control_flags,         ONLY : gamma_only
  USE noncollin_module,      ONLY : noncolin
  USE klist,                 ONLY : nkstot
  USE wannier,               ONLY : seedname, hubbard, exclude_ks_bands, &
                                    wan2hub, ispinw, iknum
  USE environment,           ONLY : environment_start, environment_end
  USE read_namelists_module, ONLY : check_namelist_read
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ios
  CHARACTER(len=4) :: spin_component
  CHARACTER(len=256) :: outdir
  LOGICAL :: needwf = .TRUE.
  !
  NAMELIST / inputpp / outdir, prefix, spin_component, seedname, &
                       hubbard, exclude_ks_bands, wan2hub
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !
  CALL environment_start ( 'WANNIER2PW' )
  !
  CALL start_clock( 'init_wan2pw' )
  !
  ! Read input on i/o node and broadcast to the rest
  !
  IF (ionode) THEN
     !
     ! Check to see if we are reading from a file
     !
     CALL input_from_file()
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     !
     prefix = ' '
     seedname = 'wannier'
     spin_component = 'none'
     hubbard = .false.
     exclude_ks_bands = 0
     wan2hub(:) = .false.
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat=ios)
     !
     !     Check of namelist variables
     !
     tmp_dir = trimcheck(outdir)
     ! back to all nodes
  ENDIF
  !
  CALL mp_bcast(ios,ionode_id, world_comm)
  CALL check_namelist_read(ios, 5, "inputpp")
  !
  ! broadcast input variable to all nodes
  !
  CALL mp_bcast(outdir,ionode_id, world_comm)
  CALL mp_bcast(tmp_dir,ionode_id, world_comm)
  CALL mp_bcast(prefix,ionode_id, world_comm)
  CALL mp_bcast(seedname,ionode_id, world_comm)
  CALL mp_bcast(spin_component,ionode_id, world_comm)
  CALL mp_bcast(hubbard,ionode_id, world_comm)
  CALL mp_bcast(exclude_ks_bands,ionode_id, world_comm)
  CALL mp_bcast(wan2hub,ionode_id, world_comm)
  !
  ! Check: bands distribution not implemented
  IF (nbgrp > 1) CALL errore('wannier2pw', 'bands (-nb) not implemented', nbgrp)
  !
  ! Read the data produced by pw.x
  CALL read_file_new ( needwf )
  !
  IF (npool > 1 .and. gamma_only) CALL errore('wannier2pw', &
      'pools not compatible with gamma_only', 1)
  !
  IF (noncolin.and.gamma_only) CALL errore('wannier2pw',&
       'Non-collinear and gamma_only not implemented',1)
  !
  SELECT CASE ( trim( spin_component ) )
  CASE ( 'up' )
     WRITE(stdout, '(5x,A)') ' Spin CASE ( up )'
     ispinw  = 1
     iknum   = nkstot/2
  CASE ( 'down' )
     WRITE(stdout, '(5x,A)') ' Spin CASE ( down )'
     ispinw = 2
     iknum   = nkstot/2
  CASE DEFAULT
     IF(noncolin) THEN
        WRITE(stdout, '(5x,A)') ' Spin CASE ( non-collinear )'
     ELSE
        WRITE(stdout, '(5x,A)') ' Spin CASE ( default = unpolarized )'
     ENDIF
     ispinw = 0
     iknum   = nkstot
  END SELECT
  !
  IF (hubbard) THEN
     ! DFT+U-related part
     CALL wannier2hubbard ()
  ELSE
     ! This program can be used for other purposes in the future. 
     ! For now just stop if hubbard = .false.
     CALL errore('wannier2pw', &
      'wannier2pw currently works only when hubbard = .true.', 1)
  ENDIF
  !  
  CALL environment_end ( 'WANNIER2PW' )
  !   
  CALL stop_pp
  !
  STOP
  !
END PROGRAM wannier2pw
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE wannier2hubbard ()
   !-------------------------------------------------------------------------
   !
   ! Initially, the Hubbard projectors of the DFT+U scheme are defined
   ! using the (ortho-)atomic orbitals (read from PP) as a basis set.
   ! This subroutine replaces the (ortho-)atomic orbitals by the
   ! Wannier functions (WFs) which are generated using the WANNIER90 code.
   ! The matrices U are stored in u_matrix: U(i,j) = <KS_i|wann_j>
   ! The WFs are written in the files $PREFIX.hub which are
   ! located in the temporary directory.
   ! Inspired by PP/src/poormanwannier.f90
   ! For more detail see Ref. [1] arXiv:2411.03937
   !
   ! Written by Iurii Timrov (November 2024)
   ! Extended to the disentanglement case by Alberto Carta (November 2024)
   !
   USE kinds,                      ONLY : DP
   USE io_global,                  ONLY : stdout, ionode, ionode_id
   USE mp_world,                   ONLY : world_comm
   USE mp,                         ONLY : mp_bcast
   USE ldaU,                       ONLY : nwfcU, wfcU
   USE cell_base,                  ONLY : bg
   USE klist,                      ONLY : nks, nkstot, xk, ngk
   USE io_files,                   ONLY : nwordwfcU, iunhub, &
                                          diropn, prefix, restart_dir
   USE noncollin_module,           ONLY : npol, noncolin
   USE wvfct,                      ONLY : npwx, nbnd
   USE wannier,                    ONLY : seedname, exclude_ks_bands, wan2hub, &
                                          ispinw, iknum
   USE wavefunctions,              ONLY : evc
   USE lsda_mod,                   ONLY : lsda, isk
   USE pw_restart_new,             ONLY : read_collected_wfc
   !
   IMPLICIT NONE
   !
   ! ... Local variables ...
   !
   COMPLEX(DP), ALLOCATABLE :: u_matrix(:,:,:), u_dis(:,:,:), pp(:,:,:)
   ! U gauge matrix, disentanglement matrix and projector matrix (product of the two)
   REAL(DP), ALLOCATABLE :: kpt_latt(:,:)
   ! k points coordinates
   REAL(DP) :: U_re, U_im, diff
   !string for the filenames for umat and udis
   CHARACTER(LEN=500) :: filename_umat, filename_udis
   INTEGER :: i, j, k, ig, matunit, matunit2, ierr, ierr2, &
      ik,          & ! dummy index to number k points
      ik2,         & ! dummy index to number k points
      iknum_,      & ! total number of k points
      n_wannier,   & ! number of Wannier functions
      n_dis_bands, & ! number of bands to be disentangled
      npw,         & ! number of plane waves
      counter        ! counter of Wannier functions which are selected
                     ! for the Hubbard manifold
   INTEGER :: igk_(npwx) ! dummy array (not used)
   LOGICAL :: exst, udis_exists, umat_exists
   ! 
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   ! Translations between Wannier2PW and Wannier90
   !    Wannier2PW   <==>   Wannier90
   !    nbnd                num_bands_tot
   !    n_wannier           num_wann
   !    num_bands           num_bands
   !    nat                 num_atoms
   !    iknum_              num_kpts
   !
   WRITE( stdout, '(/5x,79("="))')
   WRITE( stdout, '(/5x,"Please cite this paper, which describes the Wannier90 to DFT+Hubbard interface,")' )
   WRITE( stdout, '(/5x,"in publications and presentations arising from this work:")' )
   WRITE( stdout, '(/5x,"A. Carta, I. Timrov, P. Mlkvik, A. Hampel, C. Ederer, arXiv:2411.03937.")')
   WRITE( stdout, '(/5x,79("="))')
   !
   WRITE( stdout, '(/5x,"Setting up the Hubbard projectors using Wannier functions ... ")')
   !
   nwordwfcU = npwx * nwfcU * npol
   !
   IF (noncolin) CALL errore('wannier2hubbard', &
                'Noncollinear case is currently not supported', 1)
   !
   ! Open files to store the new Hubbard projectors (WFs).
   ! Note: If these files already exist, they contain ortho-atomic projectors.
   ! In this case they will be overwritten by the Wannier projectors.
   CALL diropn( iunhub, 'hub', 2*nwordwfcU, exst )
   !
   ! Read the matrices generated by Wannier90
   !
   IF (ionode) THEN
      !
      filename_umat = TRIM(seedname)//'_u.mat'
      umat_exists = file_exists(filename_umat)
      WRITE( stdout, '(/5x,"Checking the existence of the Umat file...")')
      WRITE( stdout, '(5x,A)' ) TRIM(filename_umat)
      IF (umat_exists) THEN
         WRITE( stdout, '(5x,"Umat file exists")')
      ELSE
         WRITE( stdout, '(5x,"Umat file does not exist")')
         CALL errore('wannier2hubbard', 'Umat file does not exist', 1)
      ENDIF
      !
      filename_udis = TRIM(seedname)//'_u_dis.mat'
      udis_exists = file_exists(filename_udis)
      WRITE( stdout, '(/5x,"Checking the existence of the Udis file...")')
      WRITE( stdout, '(5x,A)' ) TRIM(filename_udis)
      IF (udis_exists) THEN
         WRITE( stdout, '(5x,"Udis file exists")')
      ELSE
         WRITE( stdout, '(/5x,"Udis file does not exist")')
      ENDIF
      !
      OPEN(NEWUNIT=matunit, FILE=TRIM(seedname)//'_u.mat', &
          FORM='formatted', STATUS = 'old', IOSTAT=ierr)
      IF (ierr /= 0 ) CALL errore('wannier2hubbard', &
         'Error opening'//TRIM(seedname)//'_u.mat', ABS(ierr) )
      READ(matunit,*) ! skip one line which contains the header
      READ(matunit,*) iknum_, n_wannier, n_wannier
      !
      IF (udis_exists) THEN
         ! The disentanglement case
         OPEN(NEWUNIT=matunit2, FILE=TRIM(seedname)//'_u_dis.mat', &
            FORM='formatted', STATUS = 'old', IOSTAT=ierr2)
         IF (ierr2 /= 0 ) CALL errore('wannier2hubbard', &
            'Error opening'//TRIM(seedname)//'_u_dis.mat', ABS(ierr2) )
         READ(matunit2,*) ! skip one line which contains the header
         READ(matunit2,*) iknum_, n_wannier, n_dis_bands
      ELSE
         ! No disentanglement
         ! n_wannier is the same as the number of bands without disenatnglement
         n_dis_bands = n_wannier
      ENDIF
      !
      ! Write a summary
      WRITE( stdout, '(/,5(5x,A,5x,I8,/))')                         &
         'total number of k points                    ', nkstot,    &
         'total number of atomic wfcs with U          ', nwfcU,     &
         'total number of WFs                         ', n_wannier, &
         'total number of Kohn-Sham bands             ', nbnd,      &
         'number of Kohn-Sham bands used to build WFs ', n_dis_bands
      !
      ! Check on the number of k points
      IF (iknum_ /= iknum) CALL errore('wannier2hubbard', &
         'Mismatch in the number of k points: Input vs U matrix', 1)
      ! Check on the number of Wannier functions
      IF (n_wannier < nwfcU) CALL errore('wannier2hubbard', &
         'Too few Wannier functions: n_wannier < nwfcU', 1)
      !
   ENDIF
   !
   CALL mp_bcast (n_wannier, ionode_id, world_comm)
   CALL mp_bcast (n_dis_bands, ionode_id, world_comm)
   ! 
   ! Checks on exclude_ks_bands
   IF ( exclude_ks_bands < 0 ) CALL errore ('wannier2hubbard', &
            & 'exclude_ks_bands cannot be negative',1)
   IF ( exclude_ks_bands > nbnd ) CALL errore ('wannier2hubbard', &
            & 'exclude_ks_bands cannot be larger than nbnd',1)
   IF ( exclude_ks_bands > (nbnd-n_dis_bands) ) THEN
      WRITE(stdout,'(5x, "exclude_ks_bands =", I4)') exclude_ks_bands
      WRITE(stdout,'(5x,a)') 'WARNING: exclude_ks_bands is larger than the number' 
      WRITE(stdout,'(5x,a)') 'of available KS bands outside the Wannierization window' 
      WRITE(stdout,'(5x,a)') 'Possible reason: It is likely that you are re-running this calculation'
      WRITE(stdout,'(5x,a)') 'after DFT+U with WFs with a different number of Kohn-Sham bands (nbnd)'
      WRITE(stdout,'(5x,a)') 'Solution: Start everything from scratch'
      CALL errore ('wannier2hubbard', &
            & 'Mismatch in the number of bands',1)
   ENDIF
   !
   ! Read the U matrices on the IO rank
   ALLOCATE (kpt_latt(3,iknum))
   ALLOCATE (u_matrix(n_wannier,n_wannier,iknum))
   ALLOCATE (pp(n_dis_bands,n_wannier,iknum))
   IF (udis_exists) ALLOCATE (u_dis(n_dis_bands,n_wannier,iknum))
   !
   IF ( ionode ) THEN
      DO ik = 1, iknum
         READ(matunit,*)                ! skip one empty line
         READ(matunit,*) kpt_latt(:,ik) ! read the k point coordinates (crystal)
         DO j = 1, n_wannier
            DO i = 1, n_wannier
               READ(matunit,*) U_re, U_im
               u_matrix(i,j,ik) = CMPLX(U_re, U_im, kind=DP)
            ENDDO
         ENDDO
         IF (udis_exists) THEN
            ! Disentanglement case
            READ(matunit2,*)                ! skip one empty line
            READ(matunit2,*) kpt_latt(:,ik) ! read the k point coordinates (crystal)
            DO j = 1, n_wannier
               DO i = 1, n_dis_bands
                  READ(matunit2,*) U_re, U_im
                  u_dis(i,j,ik) = CMPLX(U_re, U_im, kind=DP)
               ENDDO
            ENDDO
            pp(:,:,ik) = MATMUL(u_dis(:,:,ik),u_matrix(:,:,ik))
         ELSE
            ! No disentanglement
            pp(:,:,ik) = u_matrix(:,:,ik)
         ENDIF
      ENDDO
      !
      ! Transform k points from crystal to cartesian coordinates
      CALL cryst_to_cart(iknum, kpt_latt, bg, 1)
      diff = 0.0d0
      DO ik = 1, iknum
         diff = diff + (SUM(kpt_latt(:,ik)-xk(:,ik)))
      ENDDO
      IF (diff.GT.1.D-06) CALL errore('wannier2hubbard', &
         'Mismatch between k points',1)
      !
      CLOSE( UNIT=matunit, STATUS='keep' )
      IF (udis_exists) CLOSE( UNIT=matunit2, STATUS='keep' )
      !
   ENDIF
   !
   DEALLOCATE (kpt_latt)
   DEALLOCATE (u_matrix)
   IF (udis_exists) DEALLOCATE (u_dis)
   !
   CALL mp_bcast (pp, ionode_id, world_comm)
   !
   ! Various checks
   !
   IF (SIZE(wan2hub) < n_wannier) &
      CALL errore('wannier2hubbard', &
      'Increase the size of the array wan2hub and recompile the code',1)
   counter = 0
   DO i = 1, SIZE(wan2hub)
      IF (wan2hub(i)) counter = counter + 1
   ENDDO
   WRITE( stdout, '(2(/5x,A,5x,I4,/))')           &
      'number of basis functions needed for the Hubbard projectors          ', nwfcU, &
      'number of Wannier functions selected to build the Hubbard projectors ', counter
   IF (counter.NE.nwfcU) CALL errore('wannier2hubbard', &
      'Mismatch between the number of selected WFs and the size of the Hubbard projectors basis',1)
   !
   ALLOCATE(wfcU(npwx*npol,nwfcU))
   !
   ! Rotate the KS wavefunctions using the pp matrix
   !
   WRITE(stdout,'(5x,a)') 'Building up the Wannier functions...'
   WRITE(stdout,'(5x,a,I5)') 'Number of local k points = ', nks
   DO ik = 1, nks
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      WRITE (stdout,'(i8)',advance='no') ik
      IF ( MOD(ik,10) == 0 ) WRITE (stdout,*)
      FLUSH(stdout)
      !
      IF (ispinw==0) THEN
         ! spin-unpolarized case
         ik2 = global_kpoint_index(nkstot, ik) 
      ELSEIF (ispinw==1) THEN
         ! spin up
         ik2 = global_kpoint_index(nkstot, ik)
      ELSEIF (ispinw==2) THEN
         ! spin down
         ik2 = global_kpoint_index(nkstot, ik) - iknum
      ELSE
        CALL errore('wannier2hubbard', 'unknown case for ispinw',1)
      ENDIF
      !
      npw = ngk(ik)
      !
      ! Read Kohn-Sham (KS) wavefunctions at k point ik
      CALL read_collected_wfc ( restart_dir(), ik, evc )
      !
      ! Apply the U matrix to KS wavefunctions and then sum up
      ! |phi_i> = \sum_j |psi_j> * U_ji
      ! The result is stored in the array wfcU which is
      ! used in DFT+U to store the basis functions of the Hubbard projectors.
      ! See Eq. (4) in Ref. [1]. Note, the phase factor of the Kohn-Sham
      ! wavefunction exp(ik*r) is not included; i.e. we work here only with the
      ! lattice-periodic part of the Kohn-Sham wavefunction.
      !
      k = 0
      DO i = 1, n_wannier
         IF (wan2hub(i)) THEN
            k = k + 1
            wfcU(:,k) = (0.0d0, 0.0d0)
            DO j = 1, n_dis_bands
               DO ig = 1, npw
                  wfcU(ig,k) = wfcU(ig,k) + evc(ig,exclude_ks_bands+j) * pp(j,i,ik2)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      !
      ! Write Bloch sums of Wannier functions to files $prefix.hub
      ! |WF_{i,k}> (no ultrasoft S)
      CALL davcio (wfcU, 2*nwordwfcU, iunhub, ik, 1)
      !
   ENDDO
   !
   DEALLOCATE (pp)
   !
   WRITE( stdout, *)
   WRITE( stdout, '(/5x,"Wannier functions are successfully written to files .hub in outdir!")')
   WRITE( stdout, '(5x,"Now this is the new basis for the Hubbard projectors.")')
   WRITE( stdout, '(5x,"Set HUBBARD{wf} in the PW input file to perform DFT+U calculations with these projectors...")')
   !
   CLOSE ( UNIT=iunhub, STATUS='KEEP' )
   !
   ! Since we called read_file, many arrays were allocated,
   ! and now we need to deallocate all of them.
   CALL clean_pw(.true.)
   !
   RETURN
   !
CONTAINS
   !
   FUNCTION file_exists(filename) RESULT(res)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: filename
      LOGICAL                     :: res
      ! Check if the file exists
      INQUIRE ( FILE=TRIM(filename), EXIST=res)
   END FUNCTION
   !
END SUBROUTINE wannier2hubbard












