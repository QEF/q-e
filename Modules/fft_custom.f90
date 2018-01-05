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
  USE fft_ggen, ONLY: fft_set_nl, fft_set_nlm
  
  IMPLICIT NONE

!--------------------------------------------------------------------
CONTAINS
!=----------------------------------------------------------------------------=!
  !
  !--------------------------------------------------------------------
  SUBROUTINE ggent(dfftt, gcutmt, gcutt, ngmt_g, gt, ggt, gstart_t, npwt )
    !--------------------------------------------------------------------
    !
    ! Initialize g-vectors for custom grid
    !
    ! FIXME: Should be merged with ggen
    !
    USE kinds,              ONLY : DP
    USE cell_base,          ONLY : at, bg
    USE control_flags,      ONLY : gamma_only
    USE constants,          ONLY : eps8
    USE mp,                 ONLY: mp_max, mp_sum
    
    IMPLICIT NONE
    
    TYPE (fft_type_descriptor), INTENT(INOUT) :: dfftt 
    REAL(DP), INTENT(IN):: gcutmt, gcutt
    ! Total number of G-vectors in custom grid
    INTEGER, INTENT(OUT):: ngmt_g
    REAL(kind=DP), INTENT(OUT), DIMENSION(:), POINTER :: ggt
    REAL(kind=DP), INTENT(OUT), DIMENSION(:,:),POINTER :: gt
    INTEGER, INTENT(OUT):: gstart_t, npwt
    !
    INTEGER,  DIMENSION(:), ALLOCATABLE :: mill(:,:)
    INTEGER :: ngmt, ngmx, n1, n2, n3, n1s, n2s, n3s
    REAL(DP):: t (3), tt, swap
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing all g vectors, on all processors: replicated data
    INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors:
    !     replicated data
    INTEGER, ALLOCATABLE :: igsrt(:)
    !
    INTEGER :: m1, m2, mc
    INTEGER :: i, j, k, ipol, ng, igl, iswap, indsw, ni, nj, nk
    !    
    ngmt = dfftt%ngm
    !
    !  calculate sum over all processors
    !
    ngmt_g = ngmt
    CALL mp_sum( ngmt_g, dfftt%comm )
    !
    !  allocate arrays - only those that are always kept until the end
    !
    ALLOCATE( ggt(ngmt) )
    ALLOCATE( gt (3, ngmt) )
    !
    ALLOCATE( mill_g( 3, ngmt_g ) )
    ALLOCATE( mill_unsorted( 3, ngmt_g ) )
    ALLOCATE( igsrt( ngmt_g ) )
    ALLOCATE( g2sort_g( ngmt_g ) )
   
    g2sort_g(:) = 1.0d20
    !
    ! save present value of ngm in ngmx variable
    !
    ngmx = ngmt
    !
    ngmt = 0
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftt%nr1-1)/2
    nj = (dfftt%nr2-1)/2
    nk = (dfftt%nr3-1)/2
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
             IF (tt <= gcutmt) THEN
                ngmt = ngmt + 1
                IF (ngmt > ngmt_g) CALL errore ('ggent', 'too many g-vectors', ngmt)
                mill_unsorted( :, ngmt ) = (/ i,j,k /)
                IF ( tt > eps8 ) THEN
                   g2sort_g(ngmt) = tt
                ELSE
                   g2sort_g(ngmt) = 0.d0
                ENDIF
             ENDIF
          ENDDO kloop
       ENDDO jloop
    ENDDO iloop
    
    IF (ngmt  /= ngmt_g ) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(ngmt - ngmt_g))

    igsrt(1) = 0
    CALL hpsort_eps( ngmt_g, g2sort_g, igsrt, eps8 )
    mill_g(1,:) = mill_unsorted(1,igsrt(:))
    mill_g(2,:) = mill_unsorted(2,igsrt(:))
    mill_g(3,:) = mill_unsorted(3,igsrt(:))
    DEALLOCATE( g2sort_g, igsrt, mill_unsorted )
    ngmt = 0
    
    ngloop: DO ng = 1, ngmt_g

       i = mill_g(1, ng)
       j = mill_g(2, ng)
       k = mill_g(3, ng)
       
       IF ( dfftt%lpara ) THEN
          m1 = MOD (i, dfftt%nr1) + 1
          IF (m1 < 1) m1 = m1 + dfftt%nr1
          m2 = MOD (j, dfftt%nr2) + 1
          IF (m2 < 1) m2 = m2 + dfftt%nr2
          mc = m1 + (m2 - 1) * dfftt%nr1x
          IF ( dfftt%isind ( mc ) == 0) CYCLE ngloop
       END IF
       
       ngmt = ngmt + 1
       
       !  To map local (ngmt) and global (ng) G-vector index: 
       !     ig_l2gt( ngmt ) = ng
       !  The global G-vector arrangement depends on the number of processors
       !
       
       gt (1:3, ngmt) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
       ggt (ngmt) = SUM( gt(1:3, ngmt)**2 )
       
       IF (ngmt > ngmx) CALL errore ('ggent', 'too many g-vectors', ngmt)
    ENDDO ngloop

    DEALLOCATE( mill_g )

    IF (ngmt /= ngmx) &
         CALL errore ('ggent', 'g-vectors missing !', ABS(ngmt - ngmx))
    !
    !     determine first nonzero g vector
    !
    IF (ggt(1).LE.eps8) THEN
       gstart_t=2
    ELSE
       gstart_t=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    ALLOCATE( mill(3,ngmt) )
    !
    CALL fft_set_nl ( dfftt, at, gt, mill  )
    IF ( gamma_only) CALL fft_set_nlm (dfftt, mill)
    !
    !     Miller indices no longer needed
    !
    DEALLOCATE ( mill )
    !
    ! set npwt - it should eventually be calculated somewhere else with 
    ! n_plane_waves() but it is good enough for gamma_only case

    IF(gamma_only) THEN
       npwt=0
       DO ng = 1, ngmt
          tt = gt(1,ng)**2 + gt(2,ng)**2 + gt(3, ng)**2
          IF (tt <= gcutt) THEN
             npwt = npwt + 1
          ENDIF
       ENDDO
    ENDIF

    RETURN
    !    
  END SUBROUTINE ggent
  
  !----------------------------------------------------------------------------  
END MODULE fft_custom
