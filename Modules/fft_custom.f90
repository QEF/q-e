!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
! Module containing routines for fft with a custom energy cutoff
!--------------------------------------------------------------------
!
MODULE fft_custom

  USE kinds, ONLY: DP
  USE parallel_include
  
  USE fft_types, ONLY: fft_type_descriptor
  
  IMPLICIT NONE

  TYPE fft_cus
  
     ! ... data structure containing all information
     ! ... about fft data distribution for a given
     ! ... potential grid, and its wave functions sub-grid.

     TYPE ( fft_type_descriptor ) :: dfftt 
     ! descriptor for the custom grid

     REAL(kind=DP) :: ecutt
     ! Custom cutoff (rydberg)
     REAL(kind=DP) :: dual_t
     ! Dual factor
     REAL(kind=DP) :: gcutmt
     INTEGER :: nr1t,nr2t,nr3t
     INTEGER :: nrx1t,nrx2t,nrx3t
     INTEGER :: nrxxt
     INTEGER :: ngmt,ngmt_l,ngmt_g
     INTEGER, DIMENSION(:), POINTER :: nlt,nltm
     REAL(kind=DP), DIMENSION(:), POINTER :: ggt
     REAL(kind=DP), DIMENSION(:,:),POINTER :: gt
     INTEGER, DIMENSION(:), POINTER :: ig_l2gt 
     INTEGER :: gstart_t
     INTEGER,  DIMENSION(:), POINTER :: ig1t,ig2t,ig3t
     INTEGER :: nlgt
     INTEGER :: npwt,npwxt
     LOGICAL :: initialized = .FALSE.
     
  END TYPE fft_cus


!--------------------------------------------------------------------
CONTAINS
!=----------------------------------------------------------------------------=!

     SUBROUTINE gvec_init( fc, ngm_, comm )
       !
       ! Set local and global dimensions, allocate arrays
       !
       USE mp, ONLY: mp_max, mp_sum
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngm_
       INTEGER, INTENT(IN) :: comm  ! communicator of the group on which g-vecs are distributed
       TYPE(fft_cus), INTENT(INOUT) :: fc
       !
       fc%ngmt = ngm_
       !
       !  calculate maximum over all processors
       !
       fc%ngmt_l = ngm_
       CALL mp_max( fc%ngmt_l, comm )
       !
       !  calculate sum over all processors
       !
       fc%ngmt_g = ngm_
       CALL mp_sum( fc%ngmt_g, comm )
       !
       !  allocate arrays - only those that are always kept until the end
       !
       ALLOCATE( fc%ggt(fc%ngmt) )
       ALLOCATE( fc%gt (3, fc%ngmt) )
!       ALLOCATE( mill(3, fc%ngmt) )
       ALLOCATE( fc%nlt (fc%ngmt) )
       ALLOCATE( fc%nltm(fc%ngmt) )
       ALLOCATE( fc%ig_l2gt(fc%ngmt) )
!       ALLOCATE( igtongl(fc%ngmt) )
       !
       RETURN 
       !
     END SUBROUTINE gvec_init



  !
  !--------------------------------------------------------------------
  SUBROUTINE set_custom_grid(fc)
    !-----------------------------------------------------------------------
    !     This routine computes the dimensions of the minimum FFT grid
    !     compatible with the input cut-off
    !
    !     NB: The values of nr1, nr2, nr3 are computed only if they are not
    !     given as input parameters. Input values are kept otherwise.
    !
    USE cell_base,   ONLY : at, tpiba2
    USE fft_support, ONLY : allowed
    
    IMPLICIT NONE
    
    TYPE(fft_cus) :: fc
    
    INTEGER, PARAMETER :: nmax = 5000
    ! an unreasonably big number for a FFT grid
    !
    ! the values of nr1, nr2, nr3 are computed only if they are not given
    ! as input parameters
    !

    fc%nr1t=0
    fc%nr2t=0
    fc%nr3t=0
    
    IF (fc%nr1t == 0) THEN
       !
       ! estimate nr1 and check if it is an allowed value for FFT
       !
       fc%nr1t = INT(2 * SQRT(fc%gcutmt) * SQRT(at(1, 1)**2 + &
            &at(2, 1)**2 + at(3, 1)**2) ) + 1  
10     CONTINUE
       IF (fc%nr1t > nmax) &
            CALL errore ('set_custom_grid', 'nr1 is unreasonably large', fc%nr1t)
       IF (allowed (fc%nr1t) ) GOTO 15
       fc%nr1t = fc%nr1t + 1
       GOTO 10
    ELSE
       IF (.NOT.allowed (fc%nr1t) ) CALL errore ('set_custom_grid', &
            'input nr1t value not allowed', 1)
    ENDIF
15  CONTINUE
    !
    IF (fc%nr2t == 0) THEN
       !
       ! estimate nr1 and check if it is an allowed value for FFT
       !
       fc%nr2t = INT(2 * SQRT(fc%gcutmt) * SQRT(at(1, 2)**2 + &
            &at(2, 2)**2 + at(3, 2)**2) ) + 1 
20     CONTINUE
       IF (fc%nr2t > nmax) &
            CALL errore ('set_custom_grid', 'nr2t is unreasonably large', fc%nr2t)
       IF (allowed (fc%nr2t) ) GOTO 25
       fc%nr2t = fc%nr2t + 1
       GOTO 20
    ELSE
       IF (.NOT.allowed (fc%nr2t) ) CALL errore ('set_fft_dim', &
            'input nr2t value not allowed', 2)
    ENDIF
25  CONTINUE
    !
    IF (fc%nr3t == 0) THEN
       !
       ! estimate nr3 and check if it is an allowed value for FFT
       !
       fc%nr3t = INT(2 * SQRT(fc%gcutmt) * SQRT(at(1, 3) **2 + &
            &at(2, 3)**2 + at(3, 3) **2) ) + 1
30     CONTINUE
       IF (fc%nr3t > nmax) &
            CALL errore ('set_custom_grid', 'nr3 is unreasonably large', fc%nr3t)
       IF (allowed (fc%nr3t) ) GOTO 35
       fc%nr3t = fc%nr3t + 1
       GOTO 30
    ELSE
       IF (.NOT.allowed (fc%nr3t) ) CALL errore ('set_custom_grid', &
            'input nr3t value not allowed', 3)
    ENDIF
35  CONTINUE
    !
    !    here we compute nr3s if it is not in input
    !
    RETURN
  END SUBROUTINE set_custom_grid
  !
  !--------------------------------------------------------------------
  SUBROUTINE ggent(fc)
    !--------------------------------------------------------------------
    !
    USE kinds,              ONLY : DP
    USE cell_base,          ONLY : at, bg, tpiba2
    USE control_flags,      ONLY : gamma_only
    USE constants,          ONLY : eps8
    
    IMPLICIT NONE
    
    TYPE(fft_cus) :: fc
    
    !
    REAL(DP) ::  t (3), tt, swap
    !
    INTEGER :: ngmx, n1, n2, n3, n1s, n2s, n3s
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing all g vectors, on all processors: replicated data
    INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors:
    !     replicated data
    INTEGER, ALLOCATABLE :: igsrt(:)
    !

#if defined(__MPI)
    INTEGER :: m1, m2, mc
    !
#endif
    INTEGER :: i, j, k, ipol, ng, igl, iswap, indsw, ni, nj, nk
    
    
!    ALLOCATE( fc%gt(3,fc%ngmt), fc%ggt(fc%ngmt) )
!    ALLOCATE( fc%ig_l2gt( fc%ngmt_l ) )
    ALLOCATE( mill_g( 3, fc%ngmt_g ), mill_unsorted( 3, fc%ngmt_g ) )
    ALLOCATE( igsrt( fc%ngmt_g ) )
    ALLOCATE( g2sort_g( fc%ngmt_g ) )
    ALLOCATE( fc%ig1t(fc%ngmt), fc%ig2t(fc%ngmt), fc%ig3t(fc%ngmt) )
   
    g2sort_g(:) = 1.0d20
    !
    ! save present value of ngm in ngmx variable
    !
    ngmx = fc%ngmt
    !
    fc%ngmt = 0
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (fc%dfftt%nr1-1)/2
    nj = (fc%dfftt%nr2-1)/2
    nk = (fc%dfftt%nr3-1)/2
    !
    iloop: DO i = -ni, ni
       !
       ! gamma-only: exclude space with x < 0
       !
       IF ( gamma_only .AND. i < 0) CYCLE iloop
       jloop: DO j = -nj, nj
          !
          ! gamma-only: exclude plane with x = 0, y < 0
          !
          IF ( gamma_only .AND. i == 0 .AND. j < 0) CYCLE jloop
          kloop: DO k = -nk, nk
             !
             ! gamma-only: exclude line with x = 0, y = 0, z < 0
             !
             IF ( gamma_only .AND. i == 0 .AND. j == 0 .AND. k < 0) CYCLE kloop
             t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
             tt = SUM(t(:)**2)
             IF (tt <= fc%gcutmt) THEN
                fc%ngmt = fc%ngmt + 1
                IF (fc%ngmt > fc%ngmt_g) CALL errore ('ggent', 'too many g-vectors', fc%ngmt)
                mill_unsorted( :, fc%ngmt ) = (/ i,j,k /)
                IF ( tt > eps8 ) THEN
                   g2sort_g(fc%ngmt) = tt
                ELSE
                   g2sort_g(fc%ngmt) = 0.d0
                ENDIF
             ENDIF
          ENDDO kloop
       ENDDO jloop
    ENDDO iloop
    
    IF (fc%ngmt  /= fc%ngmt_g ) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(fc%ngmt - fc%ngmt_g))

    igsrt(1) = 0
    CALL hpsort_eps( fc%ngmt_g, g2sort_g, igsrt, eps8 )
    mill_g(1,:) = mill_unsorted(1,igsrt(:))
    mill_g(2,:) = mill_unsorted(2,igsrt(:))
    mill_g(3,:) = mill_unsorted(3,igsrt(:))
    DEALLOCATE( g2sort_g, igsrt, mill_unsorted )
    fc%ngmt = 0
    
    ngloop: DO ng = 1, fc%ngmt_g

       i = mill_g(1, ng)
       j = mill_g(2, ng)
       k = mill_g(3, ng)
       
#if defined(__MPI)
       m1 = MOD (i, fc%dfftt%nr1) + 1
       IF (m1 < 1) m1 = m1 + fc%dfftt%nr1
       m2 = MOD (j, fc%dfftt%nr2) + 1
       IF (m2 < 1) m2 = m2 + fc%dfftt%nr2
       mc = m1 + (m2 - 1) * fc%dfftt%nr1x
       IF ( fc%dfftt%isind ( mc ) == 0) CYCLE ngloop
#endif
       
       fc%ngmt = fc%ngmt + 1
       
       !  Here map local and global g index !!!
       !  N.B. the global G vectors arrangement depends on the number of processors
       !
       fc%ig_l2gt( fc%ngmt ) = ng
       
       fc%gt (1:3, fc%ngmt) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
       fc%ggt (fc%ngmt) = SUM(fc%gt (1:3, fc%ngmt)**2)
       
       IF (fc%ngmt > ngmx) CALL errore ('ggent', 'too many g-vectors', fc%ngmt)
    ENDDO ngloop

    IF (fc%ngmt /= ngmx) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(fc%ngmt - ngmx))
    !
    !     determine first nonzero g vector
    !
    IF (fc%ggt(1).LE.eps8) THEN
       fc%gstart_t=2
    ELSE
       fc%gstart_t=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    DO ng = 1, fc%ngmt
       n1 = NINT (SUM(fc%gt (:, ng) * at (:, 1))) + 1
       fc%ig1t (ng) = n1 - 1
       IF (n1<1) n1 = n1 + fc%dfftt%nr1
       
       
       n2 = NINT (SUM(fc%gt (:, ng) * at (:, 2))) + 1
       fc%ig2t (ng) = n2 - 1
       IF (n2<1) n2 = n2 + fc%dfftt%nr2
       
       
       n3 = NINT (SUM(fc%gt (:, ng) * at (:, 3))) + 1
       fc%ig3t (ng) = n3 - 1
       IF (n3<1) n3 = n3 + fc%dfftt%nr3
       
       
       IF (n1>fc%dfftt%nr1 .OR. n2>fc%dfftt%nr2 .OR. n3>fc%dfftt%nr3) &
            CALL errore('ggent','Mesh too small?',ng)
       
#if defined (__MPI) && !defined (__USE_3D_FFT)
       fc%nlt (ng) = n3 + ( fc%dfftt%isind (n1 + (n2 - 1) * fc%dfftt%nr1x)&
            & - 1) * fc%dfftt%nr3x
#else
       fc%nlt (ng) = n1 + (n2 - 1) * fc%dfftt%nr1x + (n3 - 1) * &
            & fc%dfftt%nr1x * fc%dfftt%nr2x 
       
#endif
    ENDDO
    !
    DEALLOCATE( mill_g )
    !
    ! calculate number of G shells: ngl
    
    IF ( gamma_only) CALL index_minusg_custom(fc)
       
    !set npwt,npwxt
    !This should eventually be calculated somewhere else with 
    !n_plane_waves() but it is good enough for gamma_only

    IF(gamma_only) THEN
       fc%npwt=0
       fc%npwxt=0
       DO ng = 1, fc%ngmt
          tt = (fc%gt (1, ng) ) **2 + (fc%gt (2, ng) ) **2 + (fc%gt&
               & (3, ng) ) **2
          IF (tt <= fc%ecutt / tpiba2) THEN
             !
             ! here if |k+G|^2 <= Ecut increase the number of G
             !  inside the sphere
             !
             fc%npwt = fc%npwt + 1
          ENDIF
       ENDDO
       fc%npwxt=fc%npwt
    ENDIF

!    IF( ALLOCATED( ngmpe ) ) DEALLOCATE( ngmpe )

    RETURN
    !    
  END SUBROUTINE ggent

  !-----------------------------------------------------------------------
  SUBROUTINE index_minusg_custom(fc)
    !----------------------------------------------------------------------
    !
    !     compute indices nlm and nlms giving the correspondence
    !     between the fft mesh points and -G (for gamma-only calculations)
    !
    !
    IMPLICIT NONE
    !
    TYPE(fft_cus), INTENT(INOUT) :: fc
    !
    INTEGER :: n1, n2, n3, n1s, n2s, n3s, ng
    !
    DO ng = 1, fc%ngmt
       n1 = -fc%ig1t (ng) + 1
       IF (n1 < 1) n1 = n1 + fc%dfftt%nr1
       
       n2 = -fc%ig2t (ng) + 1
       IF (n2 < 1) n2 = n2 + fc%dfftt%nr2
       
       n3 = -fc%ig3t (ng) + 1
       IF (n3 < 1) n3 = n3 + fc%dfftt%nr3
       
       IF (n1>fc%dfftt%nr1 .OR. n2>fc%dfftt%nr2 .OR. n3>fc%dfftt%nr3) THEN
          CALL errore('index_minusg_custom','Mesh too small?',ng)
       ENDIF
       
#if defined (__MPI) && !defined (__USE_3D_FFT)
       fc%nltm(ng) = n3 + (fc%dfftt%isind (n1 + (n2 - 1) * fc&
            &%dfftt%nr1x) - 1) * fc%dfftt%nr3x
       
#else
       fc%nltm(ng) = n1 + (n2 - 1) * fc%dfftt%nr1x + (n3 - 1) * fc&
            &%dfftt%nr1x * fc%dfftt%nr1x
       
#endif
    ENDDO
    
  END SUBROUTINE index_minusg_custom
  
  SUBROUTINE deallocate_fft_custom(fc)
    !this subroutine deallocates all the fft custom stuff
    USE fft_types, ONLY : fft_type_deallocate
    
    IMPLICIT NONE

    TYPE(fft_cus) :: fc

    IF(.NOT. fc%initialized) RETURN

    DEALLOCATE(fc%nlt,fc%nltm)
    CALL fft_type_deallocate(fc%dfftt)
    DEALLOCATE(fc%ig_l2gt,fc%ggt,fc%gt)
    DEALLOCATE(fc%ig1t,fc%ig2t,fc%ig3t)
    fc%initialized=.FALSE.

    RETURN

  END SUBROUTINE deallocate_fft_custom
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reorderwfp_col ( nbands, npw1, npw2, pw1, pw2, ngwl1, ngwl2,&
       & ig_l2g1, ig_l2g2, n_g, mpime, nproc, comm )  
    !--------------------------------------------------------------------------
    !
    ! A routine using collective mpi calls that reorders the
    ! wavefunction in pw1 on a grid specified by ig_l2g1 and puts it
    ! in pw2 in the order required by ig_l2g2.
    !
    ! Can transform multiple bands at once, as specifed by the nbands  
    ! option.
    !
    ! This operation could previously be performed by calls to
    ! mergewf and splitwf however that scales very badly with number
    ! of procs.
    !
    ! Written by P. Umari, documentationa added by S. Binnie
    !
    
    USE kinds
    USE parallel_include
    USE io_global, ONLY : stdout

    IMPLICIT NONE

    INTEGER, INTENT(in)         :: npw1, npw2
    INTEGER, INTENT(IN)         :: nbands ! Number of bands to be transformed

    COMPLEX(DP), INTENT(IN)     :: pw1(npw1,nbands) ! Input wavefunction
    COMPLEX(DP), INTENT(INOUT)  :: pw2(npw2,nbands) ! Output

    INTEGER, INTENT(IN) :: mpime ! index of calling proc (starts at 0)
    INTEGER, INTENT(IN) :: nproc ! number of procs in the communicator
    INTEGER, INTENT(IN) :: comm  ! communicator

    INTEGER, INTENT(IN) :: ngwl1,ngwl2
    INTEGER, INTENT(IN) :: ig_l2g1(ngwl1),ig_l2g2(ngwl2)
    ! Global maximum number of G vectors for both grids
    INTEGER, INTENT(in) :: n_g

    
    ! Local variables
    INTEGER :: ngwl1_max, ngwl2_max, npw1_max, npw2_max, ngwl_min
    INTEGER :: gid,ierr
    INTEGER, ALLOCATABLE :: npw1_loc(:),npw2_loc(:)
    INTEGER, ALLOCATABLE :: ig_l2g1_tot(:,:),ig_l2g2_tot(:,:), itmp(:)

    INTEGER :: ii,ip,ilast,iband
    COMPLEX(kind=DP), ALLOCATABLE :: pw1_tot(:,:),pw2_tot(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: pw1_tmp(:),pw2_tmp(:), pw_global(:)


#if defined(__MPI)

    gid=comm

    ALLOCATE(npw1_loc(nproc),npw2_loc(nproc))
    !
    ! Calculate the size of the global correspondance arrays
    !
    CALL MPI_ALLREDUCE( ngwl1, ngwl1_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
    CALL MPI_ALLREDUCE( ngwl2, ngwl2_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
    CALL MPI_ALLREDUCE( npw1, npw1_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
    CALL MPI_ALLREDUCE( npw2, npw2_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
    CALL MPI_ALLGATHER( npw1, 1, MPI_INTEGER, npw1_loc, 1,&
         & MPI_INTEGER, gid, IERR )
    CALL MPI_ALLGATHER( npw2, 1, MPI_INTEGER, npw2_loc, 1,&
         & MPI_INTEGER, gid, IERR )
    !
    ALLOCATE(ig_l2g1_tot(ngwl1_max,nproc),ig_l2g2_tot(ngwl2_max&
         &,nproc))
    !
    ! All procs gather correspondance arrays
    !
    ALLOCATE(itmp(ngwl1_max))
    itmp(1:ngwl1)=ig_l2g1(1:ngwl1)
    CALL MPI_ALLGATHER( itmp, ngwl1_max, MPI_INTEGER, ig_l2g1_tot,&
         & ngwl1_max, MPI_INTEGER, gid, IERR )
    DEALLOCATE(itmp)
    !
    ALLOCATE(itmp(ngwl2_max))
    itmp(1:ngwl2)=ig_l2g2(1:ngwl2)
    CALL MPI_ALLGATHER( itmp, ngwl2_max, MPI_INTEGER, ig_l2g2_tot,&
         & ngwl2_max, MPI_INTEGER, gid, IERR)
    DEALLOCATE(itmp)
    !
    !
    ALLOCATE( pw1_tot(npw1_max,nproc), pw2_tot(npw2_max,nproc) )
    ALLOCATE( pw1_tmp(npw1_max), pw2_tmp(npw2_max) )
    ALLOCATE( pw_global(n_g) )
    !
    DO ii=1, nbands, nproc
       !
       ilast=MIN(nbands,ii+nproc-1)
       !
       ! Gather the input wavefunction.
       !
       DO iband=ii, ilast
          !
          ip = MOD(iband,nproc)      ! ip starts from 1 to nproc-1
          pw1_tmp(1:npw1)=pw1(1:npw1,iband)
          CALL MPI_GATHER( pw1_tmp, npw1_max, MPI_DOUBLE_COMPLEX,&
               & pw1_tot, npw1_max, MPI_DOUBLE_COMPLEX, ip, gid, ierr )
          !
       ENDDO
       !
       pw_global = ( 0.d0, 0.d0 )
       !
       ! Put the gathered wavefunction into the standard order.
       !
       DO ip=1,nproc
          !
          pw_global( ig_l2g1_tot(1:npw1_loc(ip), ip) ) = &
               & pw1_tot( 1:npw1_loc(ip), ip )
          !
       ENDDO
       !
       ! Now put this into the correct order for output.
       !
       DO ip=1,nproc
          !
          pw2_tot( 1:npw2_loc(ip), ip ) = &
               & pw_global ( ig_l2g2_tot(1:npw2_loc(ip),ip) )
          !
       ENDDO
       !
       ! Scatter the output wavefunction across the processors.
       !
       DO iband=ii,ilast
          !
          ip=MOD(iband,nproc)
          CALL MPI_SCATTER( pw2_tot, npw2_max, MPI_DOUBLE_COMPLEX,&
               & pw2_tmp, npw2_max, MPI_DOUBLE_COMPLEX, ip, gid, ierr )
          pw2(1:npw2,iband)=pw2_tmp(1:npw2)
          !
       ENDDO
       !
    ENDDO
    !    
    DEALLOCATE(npw1_loc,npw2_loc)
    DEALLOCATE(ig_l2g1_tot,ig_l2g2_tot)
    DEALLOCATE(pw1_tot,pw2_tot)
    DEALLOCATE(pw1_tmp,pw2_tmp)
    DEALLOCATE(pw_global)
    !
#else
    !
    ngwl_min = MIN( ngwl1, ngwl2 )
    !
    pw2(:, 1:nbands) = ( 0.0d0, 0.0d0 )    
    pw2( ig_l2g2(1:ngwl_min), 1:nbands ) = pw1( ig_l2g1(1:ngwl_min), 1:nbands ) 
    !
#endif
    !
    RETURN
    !
  END SUBROUTINE reorderwfp_col
  !----------------------------------------------------------------------------  
END MODULE fft_custom
