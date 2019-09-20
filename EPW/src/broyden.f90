  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .         
  !                                                                            
  ! Adapted from QE  
  !
  !
  !----------------------------------------------------------------------
  MODULE broyden
  !----------------------------------------------------------------------
  !! 
  !! This module contains the routines associated with Broyden's method 
  !! for potential/charge density mixing
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    ! 
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nsiter
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reache
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output Delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! Delta at previous iteration
    !
    ! Local variables
    INTEGER, PARAMETER :: maxter = 8
    !! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER :: n
    !! 
    INTEGER :: i
    !! 
    INTEGER :: j
    !! 
    INTEGER :: iwork(maxter)
    !! 
    INTEGER :: info
    !! 
    INTEGER :: iter_used
    !! 
    INTEGER :: ipos
    !! 
    INTEGER :: inext
    !! 
    !
    ! work space containing info from previous iterations:
    ! must be kept in memory and saved between calls
    REAL(KIND = DP), ALLOCATABLE, SAVE :: df(:, :), dv(:, :)
    !
    REAL(KIND = DP), ALLOCATABLE :: deltainsave(:)
    REAL(KIND = DP) :: beta(maxter,maxter), gammamix, work(maxter), norm
    REAL(KIND = DP), EXTERNAL :: DDOT, DNRM2
    ! adjustable PARAMETERs as suggested in the original paper
    REAL(KIND = DP) wg(maxter), wg0
    DATA wg0 / 0.01d0 /, wg / maxter * 1.d0 /
    !
    IF (iter < 1 ) CALL errore('mix_broyden','n_iter is smaller than 1',1)
    IF (n_iter > maxter ) CALL errore('mix_broyden','n_iter is too big',1)
    IF (ndim <= 0 ) CALL errore('mix_broyden','ndim <= 0',1)
    !
    IF (iter == 1) THEN
      IF (.NOT. ALLOCATED(df) ) ALLOCATE(df(ndim,n_iter) )    
      IF (.NOT. ALLOCATED(dv) ) ALLOCATE(dv(ndim,n_iter) )    
    ENDIF
    IF (conv .OR. iter == nsiter) THEN
      IF (ALLOCATED(df) ) DEALLOCATE(df)
      IF (ALLOCATED(dv) ) DEALLOCATE(dv)
      RETURN
    ENDIF
    IF (.NOT. ALLOCATED(deltainsave) ) ALLOCATE(deltainsave(ndim) )    
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = MIN(iter - 1, n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    !
    ipos = iter - 1 - (( iter - 2) / n_iter) * n_iter
    !
    DO n = 1, ndim
      deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
       DO n = 1, ndim
          df(n,ipos) = deltaout(n) - df(n,ipos)
          dv(n,ipos) = deltain(n)  - dv(n,ipos)
       ENDDO
       norm = ( DNRM2( ndim, df(1,ipos), 1 ) )**2.d0
       norm = SQRT(norm)
       CALL DSCAL( ndim, 1.d0/norm, df(1,ipos), 1 )
       CALL DSCAL( ndim, 1.d0/norm, dv(1,ipos), 1 )
    ENDIF
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(i,j) = wg(i) * wg(j) * DDOT( ndim, df(1,j), 1, df(1,i), 1 )
       ENDDO
       beta(i,i) = wg0**2.d0 + wg(i)**2.d0
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix 
    !
    CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
    CALL errore('broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(j, i) = beta(i, j)
       ENDDO
    ENDDO
    !
    DO i = 1, iter_used
       work(i) = DDOT( ndim, df(1,i), 1, deltaout, 1 )
    ENDDO
    !
    DO n = 1, ndim
       deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
       gammamix = 0.d0
       DO j = 1, iter_used
          gammamix = gammamix + beta(j,i) * wg(j) * work(j)
       ENDDO
       !
       DO n = 1, ndim
          deltain(n) = deltain(n) - wg(i) * gammamix * ( alphamix * df(n,i) + dv(n,i) )
       ENDDO
    ENDDO
    !
    inext = iter - ( ( iter - 1 ) / n_iter) * n_iter
    df(:,inext) = deltaout(:)
    dv(:,inext) = deltainsave(:)
    !
    IF (ALLOCATED(deltainsave) ) DEALLOCATE(deltainsave)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_broyden
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden2(ndim, deltaout, deltain, alphamix, iter, n_iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nsiter
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reache
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output Delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! Delta at previous iteration  
    !
    !   Here the local variables
    !
    ! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER, PARAMETER :: maxter = 8
    !
    INTEGER ::  n, i, j, iwork(maxter), info, iter_used, ipos, inext 
    ! work space containing info from previous iterations:
    ! must be kept in memory and saved between calls
    REAL(KIND = DP), ALLOCATABLE, SAVE :: df2(:, :), dv2(:, :)
    !
    REAL(KIND = DP), ALLOCATABLE :: deltainsave(:)
    REAL(KIND = DP) :: beta(maxter,maxter), gammamix, work(maxter), norm
    REAL(KIND = DP), EXTERNAL :: DDOT, DNRM2
    ! adjustable PARAMETERs as suggested in the original paper
    REAL(KIND = DP) wg(maxter), wg0
    DATA wg0 / 0.01d0 /, wg / maxter * 1.d0 /
    !
    IF (iter < 1 ) CALL errore('mix_broyden2','n_iter is smaller than 1',1)
    IF (n_iter > maxter ) CALL errore('mix_broyden2','n_iter is too big',1)
    IF (ndim <= 0 ) CALL errore('mix_broyden2','ndim <= 0',1)
    !
    IF (iter == 1) THEN
       IF (.NOT. ALLOCATED(df2) ) ALLOCATE(df2(ndim,n_iter) )    
       IF (.NOT. ALLOCATED(dv2) ) ALLOCATE(dv2(ndim,n_iter) )    
    ENDIF
    IF (conv .OR. iter == nsiter) THEN
       IF (ALLOCATED(df2) ) DEALLOCATE(df2)
       IF (ALLOCATED(dv2) ) DEALLOCATE(dv2)
       RETURN
    ENDIF
    IF (.NOT. ALLOCATED(deltainsave) ) ALLOCATE(deltainsave(ndim) )    
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = min(iter-1,n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    !
    ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
    !
    DO n = 1, ndim
       deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
       DO n = 1, ndim
          df2(n,ipos) = deltaout(n) - df2(n,ipos)
          dv2(n,ipos) = deltain(n)  - dv2(n,ipos)
       ENDDO
       norm = ( DNRM2( ndim, df2(1,ipos), 1 ) )**2.d0
       norm = SQRT(norm)
       CALL DSCAL( ndim, 1.d0/norm, df2(1,ipos), 1 )
       CALL DSCAL( ndim, 1.d0/norm, dv2(1,ipos), 1 )
    ENDIF
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(i,j) = wg(i) * wg(j) * DDOT( ndim, df2(1,j), 1, df2(1,i), 1 )
       ENDDO
       beta(i,i) = wg0**2.d0 + wg(i)**2.d0
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix 
    !
    CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
    CALL errore('broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(j, i) = beta(i, j)
       ENDDO
    ENDDO
    !
    DO i = 1, iter_used
       work(i) = DDOT( ndim, df2(1,i), 1, deltaout, 1 )
    ENDDO
    !
    DO n = 1, ndim
       deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
       gammamix = 0.d0
       DO j = 1, iter_used
          gammamix = gammamix + beta(j,i) * wg(j) * work(j)
       ENDDO
       !
       DO n = 1, ndim
          deltain(n) = deltain(n) - wg(i) * gammamix * ( alphamix * df2(n,i) + dv2(n,i) )
       ENDDO
    ENDDO
    !
    inext = iter - ( ( iter - 1 ) / n_iter) * n_iter
    df2(:,inext) = deltaout(:)
    dv2(:,inext) = deltainsave(:)
    !
    IF (ALLOCATED(deltainsave) ) DEALLOCATE(deltainsave)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_broyden2
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden_aniso(ik, ibnd, ndim, deltaout, deltain, alphamix, iter, n_iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nsiter
    USE eliashbergcom, ONLY : nkfs, nbndfs
    !
    IMPLICIT NONE
    ! 
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reache
    !
    INTEGER, INTENT(in) :: ik
    !! K-point index
    INTEGER, INTENT(in) :: ibnd
    !! Band index
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output Delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! Delta at previous iteration   
    ! 
    !   Here the local variables
    !
    ! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER, PARAMETER :: maxter = 8
    !
    INTEGER ::  n, i, j, iwork(maxter), info, iter_used, ipos, inext 
    ! work space containing info from previous iterations:
    ! must be kept in memory and saved between calls
    REAL(KIND = DP), ALLOCATABLE, SAVE :: df(:, :, :, :), dv(:, :, :, :)
    !
    REAL(KIND = DP), ALLOCATABLE :: deltainsave(:)
    REAL(KIND = DP) :: beta(maxter,maxter), gammamix, work(maxter), norm
    REAL(KIND = DP), EXTERNAL :: DDOT, DNRM2
    ! adjustable PARAMETERs as suggested in the original paper
    REAL(KIND = DP) wg(maxter), wg0
    DATA wg0 / 0.01d0 /, wg / maxter * 1.d0 /
    REAL(KIND = DP) :: df_(ndim,n_iter), dv_(ndim,n_iter)
    !
    IF (iter < 1 ) CALL errore('mix_broyden','n_iter is smaller than 1',1)
    IF (n_iter > maxter ) CALL errore('mix_broyden','n_iter is too big',1)
    IF (ndim <= 0 ) CALL errore('mix_broyden','ndim <= 0',1)
    !
    IF (iter == 1) THEN
       IF (.NOT. ALLOCATED(df) ) ALLOCATE(df(nbndfs,nkfs,ndim,n_iter) )    
       IF (.NOT. ALLOCATED(dv) ) ALLOCATE(dv(nbndfs,nkfs,ndim,n_iter) )    
    ENDIF
    IF (conv .OR. iter == nsiter) THEN
       IF (ALLOCATED(df)) DEALLOCATE(df)
       IF (ALLOCATED(dv)) DEALLOCATE(dv)
       RETURN
    ENDIF
    IF (.NOT. ALLOCATED(deltainsave) ) ALLOCATE(deltainsave(ndim) )    
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = min(iter-1,n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    !
    ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
    !
    DO n = 1, ndim
       deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
       DO n = 1, ndim
          df(ibnd,ik,n,ipos) = deltaout(n) - df(ibnd,ik,n,ipos)
          dv(ibnd,ik,n,ipos) = deltain(n)  - dv(ibnd,ik,n,ipos)
       ENDDO
       df_(:, :) = df(ibnd,ik,:,:)
       dv_(:, :) = dv(ibnd,ik,:,:)
       norm = ( DNRM2( ndim, df_(1,ipos), 1 ) )**2.d0
       norm = SQRT(norm)
       CALL DSCAL( ndim, 1.d0/norm, df_(1,ipos), 1 )
       CALL DSCAL( ndim, 1.d0/norm, dv_(1,ipos), 1 )
    ENDIF
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(i,j) = wg(i) * wg(j) * DDOT( ndim, df_(1,j), 1, df_(1,i), 1 )
       ENDDO
       beta(i,i) = wg0**2.d0 + wg(i)**2.d0
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix 
    !
    CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
    CALL errore('broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(j, i) = beta(i, j)
       ENDDO
    ENDDO
    !
    DO i = 1, iter_used
       work(i) = DDOT( ndim, df_(1,i), 1, deltaout, 1 )
    ENDDO
    !
    DO n = 1, ndim
       deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
       gammamix = 0.d0
       DO j = 1, iter_used
          gammamix = gammamix + beta(j,i) * wg(j) * work(j)
       ENDDO
       !
       DO n = 1, ndim
          deltain(n) = deltain(n) - wg(i) * gammamix * ( alphamix * df_(n,i) + dv_(n,i) )
       ENDDO
    ENDDO
    !
    inext = iter - ( ( iter - 1 ) / n_iter) * n_iter
    df(ibnd,ik,:,inext) = deltaout(:)
    dv(ibnd,ik,:,inext) = deltainsave(:)
    !
    IF (ALLOCATED(deltainsave) ) DEALLOCATE(deltainsave)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_broyden_aniso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden2_aniso(ik, ibnd, ndim, deltaout, deltain, alphamix, iter, n_iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nsiter
    USE eliashbergcom, ONLY : nkfs, nbndfs
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reache
    !
    INTEGER, INTENT(in) :: ik
    !! K-point index
    INTEGER, INTENT(in) :: ibnd
    !! Band index
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output Delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! Delta at previous iteration
    !
    !   Here the local variables
    !
    ! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER, PARAMETER :: maxter = 8
    !
    INTEGER ::  n, i, j, iwork(maxter), info, iter_used, ipos, inext 
    ! work space containing info from previous iterations:
    ! must be kept in memory and saved between calls
    REAL(KIND = DP), ALLOCATABLE, SAVE :: df2(:, :, :, :), dv2(:, :, :, :)
    !
    REAL(KIND = DP), ALLOCATABLE :: deltainsave(:)
    REAL(KIND = DP) :: beta(maxter,maxter), gammamix, work(maxter), norm
    REAL(KIND = DP), EXTERNAL :: DDOT, DNRM2
    ! adjustable PARAMETERs as suggested in the original paper
    REAL(KIND = DP) wg(maxter), wg0
    DATA wg0 / 0.01d0 /, wg / maxter * 1.d0 /
    REAL(KIND = DP) :: df_(ndim,n_iter), dv_(ndim,n_iter)
    !
    IF (iter < 1 ) CALL errore('mix_broyden','n_iter is smaller than 1',1)
    IF (n_iter > maxter ) CALL errore('mix_broyden','n_iter is too big',1)
    IF (ndim <= 0 ) CALL errore('mix_broyden','ndim <= 0',1)
    !
    IF (iter == 1) THEN
       IF (.NOT. ALLOCATED(df2) ) ALLOCATE(df2(nbndfs,nkfs,ndim,n_iter) )
       IF (.NOT. ALLOCATED(dv2) ) ALLOCATE(dv2(nbndfs,nkfs,ndim,n_iter) )
    ENDIF
    IF (conv .OR. iter == nsiter) THEN
       IF (ALLOCATED(df2)) DEALLOCATE(df2)
       IF (ALLOCATED(dv2)) DEALLOCATE(dv2)
       RETURN
    ENDIF
    IF (.NOT. ALLOCATED(deltainsave) ) ALLOCATE(deltainsave(ndim) )
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = min(iter-1,n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    !
    ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
    !
    DO n = 1, ndim
       deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
       DO n = 1, ndim
          df2(ibnd,ik,n,ipos) = deltaout(n) - df2(ibnd,ik,n,ipos)
          dv2(ibnd,ik,n,ipos) = deltain(n)  - dv2(ibnd,ik,n,ipos)
       ENDDO
       df_(:, :) = df2(ibnd,ik,:,:)
       dv_(:, :) = dv2(ibnd,ik,:,:)
       norm = ( DNRM2( ndim, df_(1,ipos), 1 ) )**2.d0
       norm = SQRT(norm)
       CALL DSCAL( ndim, 1.d0/norm, df_(1,ipos), 1 )
       CALL DSCAL( ndim, 1.d0/norm, dv_(1,ipos), 1 )
    ENDIF
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(i,j) = wg(i) * wg(j) * DDOT( ndim, df_(1,j), 1, df_(1,i), 1 )
       ENDDO
       beta(i,i) = wg0**2.d0 + wg(i)**2.d0
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix 
    !
    CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
    CALL errore('broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(j, i) = beta(i, j)
       ENDDO
    ENDDO
    !
    DO i = 1, iter_used
       work(i) = DDOT( ndim, df_(1,i), 1, deltaout, 1 )
    ENDDO
    !
    DO n = 1, ndim
       deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
       gammamix = 0.d0
       DO j = 1, iter_used
          gammamix = gammamix + beta(j,i) * wg(j) * work(j)
       ENDDO
       !
       DO n = 1, ndim
          deltain(n) = deltain(n) - wg(i) * gammamix * ( alphamix * df_(n,i) + dv_(n,i) )
       ENDDO
    ENDDO
    !
    inext = iter - ( ( iter - 1 ) / n_iter) * n_iter
    df2(ibnd,ik,:,inext) = deltaout(:)
    dv2(ibnd,ik,:,inext) = deltainsave(:)
    !
    IF (ALLOCATED(deltainsave) ) DEALLOCATE(deltainsave)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_broyden2_aniso
    !-----------------------------------------------------------------------
    ! 
  !-----------------------------------------------------------------------------
  END MODULE broyden
  !-----------------------------------------------------------------------------

