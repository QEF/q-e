!
! Copyright (C) 2016-2017 Quantum ESPRESSO Foundation 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_base
  !----------------------------------------------------------------------------
  !
  ! ... subroutines used to read and write binary data produced by QE
  ! ... Author: Paolo Giannozzi, based on previous work by Carlo Cavazzoni
  !
  USE kinds,     ONLY : dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_wfc, read_wfc, write_rhog, read_rhog
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, filename, ik, xk, ispin, nspin, wfc, ngw,   &
                          gamma_only, nbnd, igl, ngwl, b1,b2,b3, mill_k,    &
                          scalef, ionode_in_group, root_in_group, intra_group_comm)
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf, mergekg
      USE mp,         ONLY : mp_size, mp_rank, mp_max
      !
#if defined(__HDF5)
      USE qeh5_base_module,  ONLY  : qeh5_file, qeh5_dataset, qeh5_openfile, qeh5_open_dataset, &
                             qeh5_add_attribute, qeh5_write_dataset, qeh5_close, qeh5_set_space, &
                             qeh5_set_file_hyperslab              
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      CHARACTER(LEN=*),   INTENT(IN) :: filename
      INTEGER,            INTENT(IN) :: ik, ispin, nspin
      REAL(DP),           INTENT(IN) :: xk(:)
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      INTEGER,            INTENT(IN) :: mill_k(:,:)
      REAL(DP),           INTENT(IN) :: b1(3), b2(3), b3(3)    
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode_in_group
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      INTEGER                  :: j, ierr
      INTEGER                  :: igwx, npwx, npol
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      INTEGER, ALLOCATABLE     :: itmp(:,:)
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
#if defined(__HDF5) 
      TYPE (qeh5_file)         :: h5file
      TYPE (qeh5_dataset)      :: evc_dset, igw_dset
      !  
#endif 
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      npol = 1
      IF ( nspin == 4 ) npol = 2
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN
#if defined  __HDF5
         CALL qeh5_openfile(h5file, TRIM(filename)//'.hdf5',action = 'write') 
         CALL qeh5_add_attribute( h5file%id, "ik", ik, 1, [1])
         CALL qeh5_add_attribute( h5file%id, "xk", xk, 1, [3]) 
         CALL qeh5_add_attribute( h5file%id, "ispin", ispin ) 
         IF (gamma_only) THEN 
            CALL qeh5_add_attribute(h5file%id, "gamma_only", ".TRUE.") 
         ELSE 
            CALL qeh5_add_attribute( h5file%id, "gamma_only", ".FALSE." )
         END IF 
         CALL qeh5_add_attribute( h5file%id, "scale_factor", scalef, 1, [1]) 
         CALL qeh5_add_attribute( h5file%id, "ngw", ngw, 1, [1])
         CALL qeh5_add_attribute( h5file%id, "igwx", igwx, 1, [1])
         CALL qeh5_add_attribute( h5file%id, "npol", npol, 1, [1])
         CALL qeh5_add_attribute( h5file%id, "nbnd", nbnd, 1, [1])
#else
         OPEN ( UNIT = iuni, FILE = TRIM(filename)//'.dat', &
              FORM='unformatted', STATUS = 'unknown' )
         WRITE(iuni) ik, xk, ispin, gamma_only, scalef
         WRITE(iuni) ngw, igwx, npol, nbnd
#endif
         !
      END IF
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( itmp( 3, MAX (igwx,1) ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      itmp (:,:) = 0
      CALL mergekg( mill_k, itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         igw_dset%name = "MillerIndices"
         CALL qeh5_set_space( igw_dset, itmp(1,1), RANK = 2, DIMENSIONS = [3,igwx])
         CALL qeh5_open_dataset (h5file, igw_dset, ACTION = 'write') 
         CALL qeh5_add_attribute( igw_dset%id, "bg1", b1, RANK =1, DIMS = [3]) 
         CALL qeh5_add_attribute( igw_dset%id, "bg2", b2, RANK =1, DIMS = [3])
         CALL qeh5_add_attribute( igw_dset%id, "bg3", b3, RANK =1, DIMS = [3]) 
         CALL qeh5_add_attribute( igw_dset%id, "doc","Miller Indices of the wave-vectors, &
                                   same ordering as wave-function components") 
         CALL qeh5_write_dataset(itmp, igw_dset)
         CALL qeh5_close(igw_dset)   
       
#else
         WRITE(iuni) b1, b2, b3
         WRITE(iuni) itmp(1:3,1:igwx)
#endif
      END IF
      DEALLOCATE( itmp )
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( wtmp( MAX( npol*igwx, 1 ) ) )
      ELSE
         ALLOCATE( wtmp( 1 ) )
      ENDIF
      wtmp = 0.0_DP
      !
#if defined(__HDF5)
      IF ( ionode_in_group) THEN 
         CALL qeh5_set_space ( evc_dset, wtmp(1), 2, [npol*igwx, nbnd], MODE = 'f')
         CALL qeh5_set_space ( evc_dset, wtmp(1), 1, [npol*igwx], MODE = 'm')
         CALL qeh5_open_dataset (h5file, evc_dset, ACTION = 'write', NAME = 'evc' )       
         CALL qeh5_add_attribute( evc_dset%id, "doc:","Wave Functions, (npwx,nbnd), &
                               each contiguous line represents a wave function,  &
                               each complex coefficient is given by a couple of contiguous floats") 
      END IF
#endif 
      DO j = 1, nbnd
         !
         IF ( npol == 2 ) THEN
            !
            ! Quick-and-dirty noncolinear case - mergewf should be modified
            !
            CALL mergewf( wfc(1:npwx,       j), wtmp(1:igwx),       ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            CALL mergewf( wfc(npwx+1:2*npwx,j), wtmp(igwx+1:2*igwx), ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
         ELSE
            !
            CALL mergewf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                 nproc_in_group, root_in_group, intra_group_comm )
            !
         END IF
         !
         IF ( ionode_in_group ) THEN
#if defined(__HDF5)
            CALL qeh5_set_file_hyperslab ( evc_dset,  OFFSET = [0,j-1], COUNT = [2*npol*igwx,1] ) 
            CALL qeh5_write_dataset ( wtmp, evc_dset)   
#else
            WRITE(iuni) wtmp(1:npol*igwx)
#endif
         END IF
         !
      END DO
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL qeh5_close ( evc_dset) 
         CALL qeh5_close (h5file)   
#else 
         CLOSE (UNIT = iuni, STATUS = 'keep' )
#endif
      END IF
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !



    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, filename, ik, xk, ispin, npol, wfc, ngw, &
                         gamma_only, nbnd, igl, ngwl, b1, b2, b3, mill_k,&
                         scalef, ionode_in_group, root_in_group, &
                         intra_group_comm, ierr )
      ! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp_wave,     ONLY : splitwf, splitkg
      USE mp,          ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
#if defined (__HDF5)
      USE  qeh5_base_module
#endif

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: ik
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, npol
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      REAL(DP),           INTENT(OUT)   :: xk(3)
      REAL(DP),           INTENT(OUT)   :: b1(3), b2(3), b3(3)
      INTEGER,            INTENT(OUT)   :: mill_k(:,:)
      LOGICAL,            INTENT(OUT)   :: gamma_only
      LOGICAL,            INTENT(IN)    :: ionode_in_group
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      COMPLEX(DP), ALLOCATABLE          :: wtmp(:)
      INTEGER, ALLOCATABLE              :: itmp(:,:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, ik_
      INTEGER                           :: me_in_group, nproc_in_group
#if defined(__HDF5)
      TYPE (qeh5_file)    ::   h5file
      TYPE (qeh5_dataset) ::   h5dset_wfc, h5dset_mill
      INTEGER             ::   h5_err
      CHARACTER(LEN=8)    ::   char_buf 
#endif  
      !
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      !
#if !defined __HDF5
      IF ( ionode_in_group ) OPEN ( UNIT = iuni, FILE=TRIM(filename)//'.dat', &
           FORM='unformatted', STATUS = 'old', IOSTAT = ierr_)
      CALL mp_bcast( ierr_, root_in_group, intra_group_comm )
      IF ( PRESENT(ierr) ) THEN
         ierr = ierr_
         IF ( ierr /= 0 ) RETURN
      ELSE
         CALL errore( 'read_wfc ', &
              'cannot open restart file for reading', ierr_ )
      END IF
#endif
      IF ( ionode_in_group ) THEN
          !
#if defined  __HDF5
         CALL qeh5_openfile( h5file, TRIM(filename)//'.hdf5', ACTION = 'read', ERROR = h5_err)
         IF (h5_err /= 0 ) THEN 
            IF ( PRESENT (ierr) ) THEN 
               ierr = h5_err 
               RETURN 
            END IF
         ELSE
           CALL errore ( 'read_wfc', 'unable to open '//TRIM(filename)//'.hdf5',h5_err)  
         END IF 
         CALL qeh5_read_attribute (h5file%id, "ik", ik_)
         CALL qeh5_read_attribute (h5file%id, "xk",xk, RANK =1, DIMENSIONS = [3])
         CALL qeh5_read_attribute (h5file%id, "ispin", ispin)
         CALL qeh5_read_attribute (h5file%id, "gamma_only", char_buf, MAXLEN = len(char_buf) )
         IF (TRIM(char_buf) =='.TRUE.' .OR. TRIM(char_buf)=='.true.') THEN 
            gamma_only = .TRUE. 
         ELSE 
            gamma_only = .FALSE.
         END IF
         CALL qeh5_read_attribute (h5file%id, "scale_factor",scalef)
         CALL qeh5_read_attribute (h5file%id, "ngw", ngw)
         CALL qeh5_read_attribute (h5file%id, "nbnd", nbnd)
         CALL qeh5_read_attribute (h5file%id, "npol",npol)
         CALL qeh5_read_attribute (h5file%id, "igwx",igwx_)
#else
         READ (iuni) ik_, xk, ispin, gamma_only, scalef
         READ (iuni) ngw, igwx_, npol, nbnd
#endif
      END IF
      !
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( xk,     root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only, root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( npol,   root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd,   root_in_group, intra_group_comm )
      !
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( itmp( 3,MAX( igwx_, igwx ) ) )
#if defined(__HDF5)
       CALL qeh5_open_dataset(h5file, h5dset_mill, ACTION = 'read', NAME = 'MillerIndices')
       IF ( h5dset_mill%filespace%dims(2) .GT. MAX(igwx_, igwx)  ) &
          CALL errore ( 'read_wfc', 'real dimensions of Miller Indices dataset do not  match with igwx attribute', 8) 
       ! no reading of b1, b2, and b3 from file. They should be already set. 
       CALL qeh5_read_dataset ( itmp(:,1), h5dset_mill) 
       CALL qeh5_close ( h5dset_mill) 
#else
         READ (iuni) b1, b2, b3
         READ (iuni) itmp(1:3,1:igwx_)
#endif
         IF ( igwx > igwx_ ) itmp(1:3,igwx_+1:igwx) = 0
      ELSE
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      CALL splitkg( mill_k(:,:), itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      DEALLOCATE (itmp)
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( wtmp( npol*MAX( igwx_, igwx ) ) )
      ELSE
         ALLOCATE( wtmp(1) )
      ENDIF
#if defined (__HDF5) 
      CALL qeh5_open_dataset( h5file, h5dset_wfc, ACTION = 'read', NAME = 'evc')
      CALL qeh5_set_space ( h5dset_wfc, wtmp(1), RANK = 1, DIMENSIONS = [npol*MAX(igwx_, igwx)], MODE = 'm') 
#endif
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
#if defined __HDF5
               CALL qeh5_set_file_hyperslab (h5dset_wfc, OFFSET = [0,j-1], COUNT = [npol*igwx,1] )
               CALL qeh5_read_dataset (wtmp, h5dset_wfc )  
#else
               READ (iuni) wtmp(1:npol*igwx_) 
#endif
               IF ( igwx > igwx_ ) wtmp((npol*igwx_+1):npol*igwx) = 0.0_DP
               !
            END IF
            !
            IF ( npol == 2 ) THEN
               CALL splitwf( wfc(1:npwx,       j), wtmp(1:igwx_       ),   &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
               CALL splitwf( wfc(npwx+1:2*npwx,j), wtmp(igwx_+1:2*igwx_),  &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
            ELSE
               CALL splitwf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                    nproc_in_group, root_in_group, intra_group_comm )
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode_in_group ) THEN
#if defined (__HDF5)
         CALL qeh5_close(h5dset_wfc) 
         CALL qeh5_close(h5file)
#else
         CLOSE ( UNIT = iuni, STATUS = 'keep' )
#endif
      END IF
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rhog ( dirname, b1, b2, b3, gamma_only, mill, ig_l2g, &
         rho )
      !------------------------------------------------------------------------
      !! Write rho(G) in reciprocal space and related information to file
      !! 'charge-density.*' (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Quick-and-dirty version, allocates a large array on all mpi processes
      !! Processor "ionode" collects data from band group, writes to file
      !
      USE mp,                   ONLY : mp_sum, mp_bcast
      USE mp_bands,             ONLY : me_bgrp, root_bgrp, nproc_bgrp, &
                                       intra_bgrp_comm
      USE mp_wave,              ONLY : mergewf, mergekg
      USE io_global,            ONLY : ionode, ionode_id
#if defined (__HDF5)
      USE qeh5_base_module
#endif
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !! directory name where file is written - must end by '/'
      REAL(dp),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !!  b1, b2, b3 are the three primitive vectors in a.u.
      INTEGER,          INTENT(IN) :: mill(:,:)
      !! Miller indices for local G-vectors
      !! G = mill(1)*b1 + mill(2)*b2 + mill(3)*b3
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      LOGICAL,          INTENT(IN) :: gamma_only
      !! if true, only the upper half of G-vectors (z >=0) is present
      COMPLEX(dp),      INTENT(IN) :: rho(:,:)
      !! rho(G) on this processor
      !
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      !! Local rho(G), with LSDA workaround
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      !! Global rho(G) collected on root proc
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !! Global Miller indices collected on root proc
      INTEGER                  :: ngm, nspin, ngm_g, igwx
      INTEGER                  :: iun, ns, ig, ierr
      CHARACTER(LEN=320)       :: filename
      !
#if defined __HDF5
      TYPE (qeh5_file)          ::  h5file
      TYPE (qeh5_dataset)       ::  h5dset_mill, h5dset_rho_g
      CHARACTER(LEN=10)          :: bool_char = ".FALSE.", datasets(2) = ['rhotot_g  ','rhodiff_g ']
#endif
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (mill, 2) .OR. ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('write_rhog', 'inconsistent input dimensions', 1)
      nspin= SIZE (rho, 2)
      iun  = 4
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      filename = TRIM( dirname ) // 'charge-density.dat'
      ierr = 0
#if defined (__HDF5)
      IF ( ionode ) CALL qeh5_openfile(h5file, FILE = TRIM(dirname) //'charge-density.hdf5', &
                                                ACTION = 'write', &
                                                ERROR = ierr) 
#else
      IF ( ionode ) OPEN ( UNIT = iun, FILE = TRIM( filename ), &
                FORM = 'unformatted', STATUS = 'unknown', iostat = ierr )
#endif
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error opening file ' &
           & // TRIM( filename ), 1 )
      IF ( ionode ) THEN
#if defined(__HDF5)
          IF ( gamma_only) bool_char = '.TRUE.'
          CALL qeh5_add_attribute (h5file%id, NAME = "gamma_only", TEXT = TRIM(bool_char) )
          CALL qeh5_add_attribute (h5file%id, "ngm_g",  ngm_g ) 
          CALL qeh5_add_attribute (h5file%id, "nspin",  nspin )       
#else
          WRITE (iun, iostat=ierr) gamma_only, ngm_g, nspin
          WRITE (iun, iostat=ierr) b1, b2, b3
#endif
      END IF
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 1 )
      !
      ! ... collect all G-vectors across processors within the band group
      !
      IF ( me_bgrp == root_bgrp ) THEN
         ALLOCATE( mill_g( 3, ngm_g ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( mill_g( 3, 1 ) )
      END IF
      !
      ! ... mergekg collects distributed array mill(1:3,ig) where ig is the
      ! ... local index, into array mill_g(1:3,ig_g), where ig_g=ig_l2g(ig)
      ! ... is the global index. mill_g is collected on root_bgrp only
      !
      CALL mergekg( mill, mill_g, ngm, ig_l2g, me_bgrp, &
           nproc_bgrp, root_bgrp, intra_bgrp_comm )
      !
      ! ... write G-vectors
      !
      IF ( ionode ) THEN
#if defined(__HDF5)
         CALL qeh5_set_space ( h5dset_mill, mill_g(1,1), RANK = 2, DIMENSIONS = [3,ngm_g] ) 
         CALL qeh5_open_dataset ( h5file, h5dset_mill, NAME = "MillerIndices" , ACTION = 'write')
         !
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg1', VALUE = b1(1), RANK = 1, DIMS = [3]) 
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg2', VALUE = b2(1), RANK = 1, DIMS = [3])
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg3', VALUE = b3(1), RANK = 1, DIMS = [3]) 
         !
         CALL qeh5_write_dataset( mill_g, h5dset_mill )
         ! 
         CALL qeh5_close( h5dset_mill)    
#else
         WRITE (iun, iostat=ierr) mill_g(1:3,1:ngm_g)
#endif
      END IF
      CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... deallocate to save memory
      !
      DEALLOCATE( mill_g )
      !
      ! ... now collect all G-vector components of the charge density
      ! ... (one spin at the time to save memory) using the same logic
      !
      IF ( me_bgrp == root_bgrp ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         ! Workaround for LSDA, while waiting for much-needed harmonization:
         ! we have rhoup and rhodw, we write rhotot=up+dw and rhodif=up-dw
         ! 
         IF ( ns == 1 .AND. nspin == 2 ) THEN
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns) + rho(ig,ns+1)
            END DO
         ELSE IF ( ns == 2 .AND. nspin == 2 ) THEN
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns-1) - rho(ig,ns)
            END DO
        ELSE
            DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns)
            END DO
         END IF
         !
         rho_g = 0
         CALL mergewf( rhoaux, rho_g, ngm, ig_l2g, me_bgrp, &
              nproc_bgrp, root_bgrp, intra_bgrp_comm )
         !
         IF ( ionode ) THEN
#if defined(__HDF5)
         CALL qeh5_set_space ( h5dset_rho_g, rho_g(1), RANK = 1 , DIMENSIONS = [ngm_g] ) 
         CALL qeh5_open_dataset( h5file, h5dset_rho_g, NAME = datasets(ns) , ACTION = 'write', ERROR = ierr )
         CALL qeh5_write_dataset(rho_g, h5dset_rho_g) 
         CALL qeh5_close( h5dset_rho_g)     
#else
            WRITE (iun, iostat=ierr) rho_g(1:ngm_g)
#endif
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_bgrp_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
      END DO
      !
#if defined(__HDF5) 
      IF (ionode) CALL qeh5_close(h5file) 
#else
      IF (ionode) CLOSE (UNIT = iun, status ='keep' )
#endif
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE write_rhog
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rhog ( dirname, ig_l2g, nspin, rho )
      !------------------------------------------------------------------------
      !! Read rho(G) in reciprocal space from file  'charge-density.*' 
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Quick-and-dirty version, allocates a large array on all mpi processes
      !! All pools read the file: 1 proc reads, broadcasts to other procs
      !! Works only if there is a single band group per pool
      !
      USE mp,                   ONLY : mp_bcast
      USE mp_wave,              ONLY : splitwf
      USE mp_pools,             ONLY : me_pool, root_pool, nproc_pool, &
                                       intra_pool_comm
      !
#if defined (__HDF5) 
      USE qeh5_base_module
#endif
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !! directory name where file is read - must end by '/'
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      INTEGER,          INTENT(IN) :: nspin
      !! read up to nspin components
      COMPLEX(dp),  INTENT(INOUT) :: rho(:,:)
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      COMPLEX(dp)              :: rhoup, rhodw
      REAL(dp)                 :: b1(3), b2(3), b3(3)
      INTEGER                  :: ngm, nspin_, ngm_g, isup, isdw
      INTEGER                  :: iun, mill_dum, ns, ig, ierr
      LOGICAL                  :: ionode_k, gamma_only
      CHARACTER(LEN=320)       :: filename
      !
#if defined __HDF5
      TYPE ( qeh5_file)       :: h5file
      TYPE ( qeh5_dataset)    :: h5dset_mill, h5dset_rho_g
      CHARACTER(LEN=10)       :: tempchar, datasets(2) = ['rhotot_g  ', 'rhodiff_g ']
      !
      filename = TRIM( dirname ) // 'charge-density.hdf5'
#else 
      filename = TRIM( dirname ) // 'charge-density.dat'
#endif 
      !
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('read_rhog', 'inconsistent input dimensions', 1)
      !
      iun  = 4
      ierr = 0
      !
      ! ... the root processor of each pool reads
      !
      ionode_k = (me_pool == root_pool)
      !
      IF ( ionode_k ) THEN
#if defined (__HDF5) 
         CALL qeh5_openfile(h5file, TRIM(filename), ACTION = 'read', error = ierr)
         CALL qeh5_read_attribute (h5file%id, "gamma_only", tempchar, MAXLEN = len(tempchar)  )
         CALL qeh5_read_attribute (h5file%id, "ngm_g", ngm_g ) 
         CALL qeh5_read_attribute (h5file%id, "nspin", nspin_)  
         SELECT CASE (TRIM(tempchar) )  
            CASE ('.true.', '.TRUE.' ) 
                gamma_only = .TRUE.
            CASE DEFAULT
                gamma_only = .FALSE.
         END SELECT    
#else
         OPEN ( UNIT = iun, FILE = TRIM( filename ), &
              FORM = 'unformatted', STATUS = 'old', iostat = ierr )
         IF ( ierr /= 0 ) THEN
            ierr = 1
            GO TO 10
         END IF
         READ (iun, iostat=ierr) gamma_only, ngm_g, nspin_
         IF ( ierr /= 0 ) THEN
            ierr = 2
            GO TO 10
         END IF
         READ (iun, iostat=ierr) b1, b2, b3
         IF ( ierr /= 0 ) ierr = 3
#endif
10       CONTINUE 
      END IF
      !
      CALL mp_bcast( ierr, root_pool, intra_pool_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), ierr )
      CALL mp_bcast( ngm_g, root_pool, intra_pool_comm )
      CALL mp_bcast( nspin_, root_pool, intra_pool_comm )
      !
      IF ( nspin > nspin_ ) &
         CALL infomsg('read_rhog', 'some spin components not found')
      IF ( ngm_g < MAXVAL (ig_l2g(:)) ) &
           CALL errore('read_rhog', 'some G-vectors are missing', 1)
      !
      ! ... skip record containing G-vector indices
      !
      IF ( ionode_k ) THEN
#if !defined(__HDF5)
         READ (iun, iostat=ierr) mill_dum
#endif
      END IF
      CALL mp_bcast( ierr, root_pool, intra_pool_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... now read, broadcast and re-order G-vector components
      ! ... of the charge density (one spin at the time to save memory)
      !
      IF ( ionode_k ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         IF ( ionode_k ) THEN
#if defined(__HDF5)
            CALL qeh5_open_dataset( h5file, h5dset_rho_g, NAME = datasets(ns), ACTION = 'read', ERROR = ierr) 
            CALL qeh5_read_dataset ( rho_g , h5dset_rho_g )
            CALL qeh5_close ( h5dset_rho_g )  
#else 
            READ (iun, iostat=ierr) rho_g(1:ngm_g)
#endif
         END IF
         CALL mp_bcast( ierr, root_pool, intra_pool_comm )
         IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
              & // TRIM( filename ), 2+ns )
         !
         CALL splitwf( rhoaux, rho_g, ngm, ig_l2g, me_pool, &
              nproc_pool, root_pool, intra_pool_comm )
         DO ig = 1, ngm
            rho(ig,ns) = rhoaux(ig)
         END DO
         !
         ! Workaround for LSDA, while waiting for much-needed harmonization:
         ! if file contains rhotot=up+dw and rhodif=up-dw (nspin_=2), and
         ! if we want rhoup and rho down (nspin=2), convert 
         ! 
         IF ( nspin_ == 2 .AND. nspin == 2 .AND. ns == 2 ) THEN
            DO ig = 1, ngm
               rhoup = (rho(ig,ns-1) + rhoaux(ig)) / 2.0_dp
               rhodw = (rho(ig,ns-1) - rhoaux(ig)) / 2.0_dp
               rho(ig,ns-1)= rhoup
               rho(ig,ns  )= rhodw
            END DO
         END IF
      END DO
      !
      IF ( ionode_k ) CLOSE (UNIT = iun, status ='keep' )
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE read_rhog
    !
  END MODULE io_base
