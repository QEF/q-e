!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! Copyright (C) 2017 Andrea Dal Corso (for removing the history reset
! and the removal of maxter)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE mix_potential_eels (ndim, vout, vin, alphamix, dr2, tr2, iter, &
           file_extension, conv)
  !-----------------------------------------------------------------------
  !
  ! Modified Broyden's method for potential/charge density mixing
  !             D.D.Johnson, PRB 38, 12807 (1988)
  !
  ! This version differs from the QE routine because it does not
  ! reset the history at any iteration. The history is always saved on
  ! disk. The weights are not used. (ADC)
  !
  ! On input :
  !    ndim      dimension of arrays vout, vin
  !    vout      output potential/rho at current iteration
  !    vin       potential/rho at previous iteration
  !    alphamix  mixing factor (0 < alphamix <= 1)
  !    tr2       threshold for selfconsistency
  !    iter      current iteration number
  !    file_extension  if present save previous iterations on 
  !                    file 'prefix'.'file_extension'
  !                    otherwise keep everything in memory
  ! On output:
  !    dr2       [(vout-vin)/ndim]^2
  !    vin       mixed potential
  !    vout      vout-vin
  !    conv      true if dr2.le.tr2

  USE kinds,           ONLY : DP
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp_images,       ONLY : intra_image_comm
  USE io_global,       ONLY : ionode, ionode_id
  USE io_files,        ONLY : tmp_dir
  USE mp,              ONLY : mp_sum, mp_bcast
  USE io_files,        ONLY : diropn
  IMPLICIT NONE
  !
  !   First the dummy variables
  !
  CHARACTER (LEN=256) :: file_extension
  INTEGER  :: ndim, iter
  REAL(DP) :: vout (ndim), vin (ndim), alphamix, dr2, tr2
  LOGICAL  :: conv
  !
  !   Here the local variables
  !
  INTEGER, ALLOCATABLE :: iwork(:)
  INTEGER :: iunit, iunmix, iunmix_beta, n, i, j, info, iter_used, ipos, &
             ndimtot
  INTEGER :: find_free_unit
 
  REAL(DP), ALLOCATABLE :: df (:,:), dv (:,:), beta(:,:)
  !
  REAL(DP), ALLOCATABLE :: vinsave (:)
  REAL(DP), ALLOCATABLE :: work(:), w(:)
  REAL(DP) :: gamma, norm
  LOGICAL  :: opnd, exst
  REAL(DP), EXTERNAL :: ddot, dnrm2
  ! adjustable parameters as suggested in the original paper
  REAL(DP) :: w0
  DATA w0 / 0.01d0 /

  CHARACTER (LEN=256) :: filebeta

  !
  CALL start_clock ('mix_pot')

  IF (iter.LT.1) CALL errore ('mix_potential', 'iter is wrong', 1)
  IF (ndim.LE.0) CALL errore ('mix_potential', 'ndim .le. 0', 3)
  IF (file_extension==' ') CALL errore('mix_potential','a filename is needed',1)

  filebeta=TRIM(tmp_dir)//'/'//TRIM(file_extension)//'_beta'
  !
  !
  DO n = 1, ndim
     vout (n) = vout (n) - vin (n)
  ENDDO
  dr2 = dnrm2 (ndim, vout, 1)**2

  ndimtot = ndim
  !
  CALL mp_sum (dr2, intra_bgrp_comm)
  CALL mp_sum (ndimtot, intra_bgrp_comm)
  !
  dr2 = (SQRT(dr2) / ndimtot) **2

  conv = dr2 < tr2

  iunmix=find_free_unit()

  IF (conv) THEN
     ! remove temporary file (open and close it)
     CALL diropn (iunmix, file_extension, ndim, exst)
     CLOSE (UNIT=iunmix, STATUS='delete')
     IF (ionode) THEN
        OPEN (UNIT=iunmix, FILE=TRIM(filebeta), FORM='unformatted')
        CLOSE (UNIT=iunmix, STATUS='delete')
     ENDIF
     RETURN
  ENDIF
  !
  iter_used = iter - 1 
  CALL diropn (iunmix, file_extension, ndim, exst)

  IF (iter>1.AND..NOT.exst) THEN
     CALL infomsg ('mix_potential', 'file not found, restarting')
     iter = 1
  ENDIF
  iunmix_beta=find_free_unit()
!
!  after the first iteration we save space for df, dv, beta, w and the
!  working space for beta inversion
!


  IF (iter_used>0) THEN
     INQUIRE (UNIT = iunmix_beta, EXIST = exst)  

     IF (.NOT.exst) THEN
        CALL infomsg ('mix_potential', 'beta file not found, restarting')
        iter = 1
     ENDIF

     OPEN(UNIT=iunmix_beta, FILE=TRIM(filebeta), FORM='unformatted')


     ALLOCATE (df( ndim , iter_used))    
     ALLOCATE (dv( ndim , iter_used))    
     ALLOCATE (beta(iter_used, iter_used))
     ALLOCATE (work(iter_used))
     ALLOCATE (w(iter_used))
     ALLOCATE (iwork(iter_used))
     df = 0.0d0
     dv = 0.0d0
     w=1.0_DP
  ENDIF
  !
  ! ipos is the position in which results from the present iteraction
  ! are stored. In this routine it coincide with iter_used because we
  ! never reset the history.
  !
  ipos = iter_used
  !

  !  compute df and dv of the last iteration 
  !
  IF (iter_used>0) THEN
     CALL davcio (df (1, ipos), ndim, iunmix, 1, - 1)
     CALL davcio (dv (1, ipos), ndim, iunmix, 2, - 1)
     DO n = 1, ndim
        df (n, ipos) = vout (n) - df (n, ipos)
        dv (n, ipos) = vin (n) - dv (n, ipos)
     ENDDO
     norm = dnrm2 (ndim, df (1, ipos), 1) **2
     CALL mp_sum (norm, intra_bgrp_comm)
     norm=SQRT(norm)
     CALL dscal (ndim, 1.d0 / norm, df (1, ipos), 1)
     CALL dscal (ndim, 1.d0 / norm, dv (1, ipos), 1)
  ENDIF


  !
  ! read df and dv for all previous iterations
  !
  DO i = 1, iter_used-1
     CALL davcio (df (1, i), ndim, iunmix, 2 * i + 1, - 1)
     CALL davcio (dv (1, i), ndim, iunmix, 2 * i + 2, - 1)
  ENDDO
!
!  In positions 1 and 2 save vin and vout of the present iteration
!
  CALL davcio (vout, ndim, iunmix, 1, 1)
  CALL davcio (vin, ndim, iunmix, 2, 1)
!
!  and save df and dv of the present iteration if available
!
  IF (iter_used>0) THEN
     CALL davcio (df (1, ipos), ndim, iunmix, 2 * ipos + 1, 1)
     CALL davcio (dv (1, ipos), ndim, iunmix, 2 * ipos + 2, 1)
  ENDIF
!
!  now read the part of beta that was computed in previous iteration
!
  IF (iter_used>0) beta=0.0_DP
  IF (iter_used>2) THEN
     IF (ionode) THEN
        DO i=1, iter_used-1
           READ(iunmix_beta) (beta(i,j), j=i+1,iter_used-1)
        ENDDO
     ENDIF
     CALL mp_bcast(beta, ionode_id, intra_image_comm)
  ENDIF
!
!  add the last column of beta and recompute all the diagonal part
!
  DO i = 1, iter_used-1
     j=iter_used
     beta (i, j) = w (i) * w (j) * ddot (ndim, df (1, j), 1, df (1, i), 1)
     CALL mp_sum ( beta (i, j), intra_bgrp_comm )
     beta (i, i) = w0**2 + w (i) **2
  ENDDO
  IF (iter_used>0) beta(iter_used,iter_used)= w0**2 + w (iter_used) **2
!
!  save the beta computed so far for the next iterations
!
  IF (iter_used>1) THEN
     IF (ionode) THEN
        REWIND(iunmix_beta)
        DO i=1, iter_used
           WRITE(iunmix_beta) (beta(i,j), j=i+1,iter_used)
        ENDDO
     ENDIF
  ENDIF
!
! and invert beta (only if of dimension larger than 1)
!
  IF (iter_used>1) THEN
     CALL DSYTRF ('U', iter_used, beta, iter_used, iwork, work, iter_used, info)
     CALL errore ('broyden', 'factorization', info)
     CALL DSYTRI ('U', iter_used, beta, iter_used, iwork, work, info)
     CALL errore ('broyden', 'DSYTRI', info)
  ELSEIF (iter_used==1) THEN
     beta(1,1)=1.0_DP/beta(1,1)
  ENDIF
  !
  DO i = 1, iter_used
     DO j = i + 1, iter_used
        beta (j, i) = beta (i, j)
     ENDDO
  ENDDO
!
!  Here update vin for the next iteration. First linear mixing and then
!  subtract the terms with the previous u 
!
  DO i = 1, iter_used
     work (i) = ddot (ndim, df (1, i), 1, vout, 1)
  ENDDO
  IF (iter_used>0) CALL mp_sum ( work(1:iter_used), intra_bgrp_comm )
  !
  DO n = 1, ndim
     vin (n) = vin (n) + alphamix * vout (n)
  ENDDO
  !
  DO i = 1, iter_used
     gamma = 0.d0
     DO j = 1, iter_used
        gamma = gamma + beta (j, i) * w (j) * work (j)
     ENDDO
     !
     DO n = 1, ndim
        vin (n) = vin (n) - w (i) * gamma * (alphamix * df (n, i) + dv (n, i) )
     ENDDO
  ENDDO
  !
  CLOSE (UNIT = iunmix, STATUS='keep')
  IF (iter_used>0) THEN
     CLOSE (UNIT = iunmix_beta, STATUS='keep')
     DEALLOCATE(dv)
     DEALLOCATE(df)
     DEALLOCATE(beta)
     DEALLOCATE(iwork)
     DEALLOCATE(w)
     DEALLOCATE(work)
  ENDIF
  CALL stop_clock ('mix_pot')
  RETURN
END SUBROUTINE mix_potential_eels
