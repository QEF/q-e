!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! -----------------------------------------------------------------
! This program reads the prefix.wfc in G-space written by QE and 
! writes it in real space prefix.wfc_r.
! Warning: The wfc is written out in real space on the smooth
! grid, as such it occupies much more disk space then that in G-space.
!
! Program written by Matteo Calandra.
! 
!-----------------------------------------------------------------------
PROGRAM wfck2r
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE mp_global, ONLY : npool, mp_startup,  intra_image_comm
  USE wvfct,     ONLY : nbnd, npwx
  USE klist,     ONLY : xk, nks, ngk, igk_k
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_world,  ONLY : world_comm
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect, ONLY : ngm, g 
  USE gvecs, ONLY : nls
  USE noncollin_module, ONLY : npol, nspin_mag, noncolin
  USE environment,ONLY : environment_start, environment_end
  USE fft_base,  only : dffts
  USE scatter_mod,  only : gather_grid
  USE fft_interfaces, ONLY : invfft

  !
  IMPLICIT NONE
  CHARACTER (len=256) :: outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER            :: npw, iunitout,ios,ik,i,iuwfcr,lrwfcr,ibnd, ig, is
  LOGICAL            :: exst
  COMPLEX(DP), ALLOCATABLE :: evc_r(:,:), dist_evc_r(:,:)

  NAMELIST / inputpp / outdir, prefix


  !
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'WFCK2R' )

  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'

  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     ! 
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('WFCK2R', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir = trimcheck (outdir)
     ! 
  END IF

  !
  ! ... Broadcast variables
  !

  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  call openfil_pp

  exst=.false.

  filename='wfc_r'
  write(6,*) 'filename=',filename
  iuwfcr=877
  lrwfcr = 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x * npol
  ! lrwfc = 2 * nbnd * npwx * npol
  write(6,*) dffts%nnr, npwx
  
 

  write(6,*) 'length of wfc in real space/per band', nks*lrwfcr*8
  write(6,*) 'length of wfc in k space', 2*nbnd*npwx*nks*8
  CALL init_us_1

!
!define lrwfcr
!
  IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst)

  ALLOCATE ( evc_r(dffts%nnr,npol) )
  ALLOCATE ( dist_evc_r(dffts%nr1x*dffts%nr2x*dffts%nr3x,nspin_mag) )
  

  DO ik = 1,nks
     
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     do ibnd=1,nbnd 
        !
        ! perform the fourier transform
        !
        evc_r = (0.d0, 0.d0)     
        do ig = 1, npw
           evc_r (nls (igk_k(ig,ik) ),1 ) = evc (ig,ibnd)
        enddo
        CALL invfft ('Wave', evc_r(:,1), dffts)
        IF (noncolin) THEN
           DO ig = 1, npw
              evc_r (nls(igk_k(ig,ik)),2) = evc (ig+npwx, ibnd)
           ENDDO
           CALL invfft ('Wave', evc_r(:,2), dffts)
        ENDIF

        dist_evc_r=(0.d0,0.d0)

#if defined (__MPI)
        DO is = 1, nspin_mag
           !
           CALL gather_grid( dffts, evc_r(:,is), dist_evc_r(:,is) )
           !
        END DO
#else
        dist_evc_r(1:dffts%nnr,:)=evc_r(1:dffts%nnr,:)
#endif

        if(ionode)     call davcio (dist_evc_r, lrwfcr, iuwfcr, (ik-1)*nbnd+ibnd, +1)
     enddo
        
     !
     ! ... First task is the only task allowed to write the file
     !

  enddo

  if(ionode) close(iuwfcr)
  DEALLOCATE (evc_r)

  CALL environment_end ( 'WFCK2R' )

  CALL stop_pp
  STOP

end PROGRAM wfck2r
