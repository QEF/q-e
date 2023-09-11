!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! --------------------------------------------------------------
MODULE add_dmft_occ
  ! Written by Sophie D. Beck (2020-2021)
  ! arXiv:2111.10289 (2021)
  !
  ! This module contains a subroutine to read DMFT occupation uddaptes
  ! and to change the fermi weights and the wavefunctions accordingly.
  !
  !
  USE wvfct,                ONLY : wg, nbnd, et
  USE klist,                ONLY : wk, nkstot, nks
  USE kinds,                ONLY : dp
  USE io_global,			ONLY : ionode, ionode_id, stdout
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_images,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE control_flags,        ONLY : restart
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: dmft
  !! if .TRUE.: update occupations/eigenvectors once
  LOGICAL :: dmft_updated = .FALSE.
  !! if .TRUE.: the update was done
  CHARACTER(len=256) :: dmft_prefix
  !! name of hdf5 input archive containing the occupation update and bandwindow
  COMPLEX(DP), ALLOCATABLE :: v_dmft(:,:,:)
  !! transformation matrix from previous to new eigenvectors
  !
  PUBLIC :: dmft, dmft_updated, v_dmft
  PUBLIC :: dmft_update
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dmft_update()
#if defined(__HDF5)
    !----------------------------------------------------------------------------
    !! Reads delta_N and band_window from file, diagonalizes non-diagonal occupation
    !! matrix, writes new weights wg
    !! 
    !
    USE HDF5
	USE ISO_C_BINDING

    IMPLICIT NONE

    INTEGER :: ibnd,    &! band counter
               ik,      &! k-point counter
               nbnd_c,  &! number of bands in correlated subset from DMFT
               bnd_low, &! lower band index from where nbnd_c occupations are added
               ierr      ! hdf5 read error
    CHARACTER(LEN=32)               :: hdfg_bnd_low
    CHARACTER(LEN=32)               :: hdfd_bnd_low
    CHARACTER(LEN=32)               :: hdfg_delta_n
    CHARACTER(LEN=32)               :: hdfd_delta_n
    COMPLEX(DP), ALLOCATABLE	    :: n_dmft_root(:,:,:)
    COMPLEX(DP), ALLOCATABLE	    :: n_dmft(:,:,:)
    REAL(DP), ALLOCATABLE, TARGET   :: delta_n(:,:,:,:)
    COMPLEX(DP), ALLOCATABLE        :: n_evec_c(:,:,:)
    REAL(DP), ALLOCATABLE           :: n_eval_c(:,:)
    REAL(DP), TARGET                :: band_window(2,nkstot,1)
    !
    INTEGER(HID_T)                  :: f_id, g_id, d_id ! hdf5-related
    TYPE(C_PTR)                     :: f_ptr ! hdf5-related
    !
    ! ... allocate space to store HDF5 data, read dataset 'band_window in group 'dft_misc_input'
    ! ... and 'delta_N' in group 'dft_update' from h5, save it as complex 'n_dmft'
    !
    IF (.NOT. restart ) CALL errore('dmft_update', &
        'dmft_update not implemented for "restart_mode = from_scratch"', 1)
    !
    IF (.NOT. ALLOCATED(v_dmft)) ALLOCATE(v_dmft(nbnd, nbnd, nks))
    !
    ! ... HDF5 structure of band window with dataset dimension ( 1 x nkstot x 2) 
    ! ... modify names here if your input has a different structure
    !
    hdfg_bnd_low = 'dft_misc_input' ! group
    hdfd_bnd_low = 'band_window'    ! dataset
    !
    ! ... HDF5 structure of delta_N with dataset dimension ( nkstot x nbnd_c x nbnd_c x 2 )
    ! ... modify names here if your input has a different structure
    !
    hdfg_delta_n = 'dft_update'     ! group
    hdfd_delta_n = 'delta_N'        ! dataset
    !
    IF ( ionode ) THEN
        !
        CALL h5open_f(ierr)
        CALL h5fopen_f(TRIM(dmft_prefix)//'.h5', h5f_acc_rdwr_f, f_id, ierr)
        !
        CALL errore('dmft_update', 'hdf5 archive '//TRIM(dmft_prefix)//'.h5 not &
                    &found', abs(ierr))
        !
        CALL h5gopen_f(f_id, TRIM(hdfg_bnd_low), g_id, ierr)
        !
        CALL errore('dmft_update', 'hdf5 group '//TRIM(dmft_prefix)//'.h5/'//TRIM(hdfg_bnd_low)//' &
                    &not found', abs(ierr))
        !
        CALL h5dopen_f(g_id, TRIM(hdfd_bnd_low), d_id, ierr)
        !
        CALL errore('dmft_update', 'hdf5 dataset '//TRIM(dmft_prefix)//'.h5/'//TRIM(hdfg_bnd_low)//'/&
                    &'//TRIM(hdfd_bnd_low)//' not found', abs(ierr))
        !
        f_ptr = c_loc(band_window(1,1,1))
        CALL h5dread_f(d_id, h5t_native_double, f_ptr, ierr)
        !
        CALL errore('dmft_update', 'number of bands "nbnd" must be larger or equal than upper limit &
                     in hdf5 dataset '//TRIM(dmft_prefix)//'.h5/'//TRIM(hdfg_bnd_low)//'/&
                     &'//TRIM(hdfd_bnd_low)//'' , INT(band_window(2,1,1)) - nbnd)
        !
        ! ... for now assume that all bands have the same band_window, so take data from first k-point
        !
        bnd_low = INT(band_window(1,1,1)) - 1
        nbnd_c = INT(band_window(2,1,1)) - bnd_low
        !
        CALL h5dclose_f(d_id, ierr)
        CALL h5gclose_f(g_id, ierr)
        !
        IF (.NOT. ALLOCATED(n_dmft_root)) ALLOCATE(n_dmft_root(nbnd_c, nbnd_c, nkstot))
        IF (.NOT. ALLOCATED(delta_n)) ALLOCATE(delta_n(2,nbnd_c,nbnd_c,nkstot))
        !
        CALL h5gopen_f(f_id, TRIM(hdfg_delta_n), g_id, ierr)
        !
        CALL errore('dmft_update', 'hdf5 group '//TRIM(dmft_prefix)//'.h5/'//TRIM(hdfg_delta_n)//' &
                    &not found', abs(ierr))
        !
        CALL h5dopen_f(g_id, TRIM(hdfd_delta_n), d_id, ierr)
        !
        CALL errore('dmft_update', 'hdf5 dataset '//TRIM(dmft_prefix)//'.h5/'//TRIM(hdfg_delta_n)//'/&
                    &'//TRIM(hdfd_delta_n)//' not found', abs(ierr))
        !
        f_ptr = C_LOC(delta_n(1,1,1,1))
        CALL h5dread_f(d_id, h5t_native_double, f_ptr, ierr)
        !
        ! ... set corrections to the DFT occupations from DMFT for each band in nbnd_c at each ik
        !
        n_dmft_root(:,:,:) = CMPLX(delta_n(1,:,:,:), delta_n(2,:,:,:), KIND=DP)
        !
        CALL h5dclose_f(d_id, ierr)
        CALL h5gclose_f(g_id, ierr)
        CALL h5fclose_f(f_id, ierr)
        !
        IF (ALLOCATED(delta_n)) DEALLOCATE(delta_n)
        !
    ENDIF
    !
    CALL mp_bcast(ierr, ionode_id, intra_image_comm)
    CALL errore('dmft_update', 'reading of hdf5 archive failed', abs(ierr))
    !
    CALL mp_bcast(bnd_low, ionode_id, intra_image_comm)
    CALL mp_bcast(nbnd_c, ionode_id, intra_image_comm)
    !
    ! ... make transition from ionode n_dmft_root for nkstot to n_dmft with nks
    !
    IF ( .NOT. ionode ) ALLOCATE(n_dmft_root(nbnd_c, nbnd_c, nkstot))
    CALL mp_bcast(n_dmft_root, ionode_id, intra_image_comm)
    !
    IF (.NOT. ALLOCATED(n_dmft)) ALLOCATE(n_dmft(nbnd_c, nbnd_c, nks))
    CALL poolscatter_matrix(nbnd_c, nkstot, n_dmft_root, nks, n_dmft)
    IF (ALLOCATED(n_dmft_root)) DEALLOCATE(n_dmft_root)
    !
    IF (.NOT. ALLOCATED(n_evec_c)) ALLOCATE(n_evec_c(nbnd_c, nbnd_c, nks))
    IF (.NOT. ALLOCATED(n_eval_c)) ALLOCATE(n_eval_c(nbnd_c, nks))
    !
    ! ... for correlated subspace add DFT Fermi weights, diagonalize, sort according to
    ! ... eigenvalues and write adjusted eigenvalues and weights for these nbnd_c bands
    !
    DO ik = 1, nks
       !
       n_dmft(:,:,ik) = TRANSPOSE(n_dmft(:,:,ik))
       !
       DO ibnd = 1, nbnd_c
           n_dmft(ibnd,ibnd,ik) = n_dmft(ibnd,ibnd,ik) + CMPLX(wg(bnd_low+ibnd,ik)/wk(ik),0d0,kind=DP)
       ENDDO
       !
	   CALL cdiagh(nbnd_c, n_dmft(:,:,ik), nbnd_c, n_eval_c(:,ik), n_evec_c(:,:,ik))
       !
       DO ibnd = 1, nbnd_c
           wg(bnd_low+ibnd, ik) = n_eval_c(nbnd_c+1-ibnd, ik) * wk(ik)
       ENDDO
       !
    ENDDO
    !
    ! ... integrate diagonalization matrix n_evec_c for correlated subspace into scope of
    ! ... the total number of bands, using diagonal matrix
    !
    v_dmft(:,:,:) = (0.0_dp,0.0_dp)
    DO ik = 1, nks
        DO ibnd = 1, nbnd
            IF (ibnd > bnd_low .AND. ibnd < bnd_low + nbnd_c + 1) THEN
                v_dmft(bnd_low+1:bnd_low+nbnd_c,ibnd,ik) = n_evec_c(:,nbnd_c+1-(ibnd-bnd_low),ik)
            ELSE
                v_dmft(ibnd,ibnd,ik) = CMPLX(1d0, 0d0, kind=DP)
            ENDIF
        ENDDO
    ENDDO
	!
	IF (ALLOCATED(n_dmft)) DEALLOCATE(n_dmft)
	IF (ALLOCATED(n_evec_c)) DEALLOCATE(n_evec_c)
    IF (ALLOCATED(n_eval_c)) DEALLOCATE(n_eval_c)
    !
    CALL poolrecover( wg, nbnd, nkstot, nks )
    CALL poolrecover( et, nbnd, nkstot, nks )
    !
#endif
  END SUBROUTINE dmft_update
  !
  !------------------------------------------------------------------------
END MODULE add_dmft_occ
