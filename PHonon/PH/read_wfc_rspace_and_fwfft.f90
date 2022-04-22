!
!
subroutine read_wfc_rspace_and_fwfft( evc , ik , lrec ,  iunit , npw , igmap )
  !! This routine reads a wavefunction in real space and transform it in
  !! Fourier space.  
  !! Not tested for the non-collinear case.
  !
  ! Matteo Calandra
  !
  use kinds,           ONLY : DP
  use wvfct,       ONLY : npwx, nbnd
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE fft_base,            ONLY : dffts
  USE scatter_mod,            ONLY : scatter_grid
  USE fft_interfaces,      ONLY : fwfft
  USE io_global,             ONLY : ionode_id, ionode
  USE mp_pools,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_bcast
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik
  !! k-point to read
  INTEGER, INTENT(IN) :: lrec
  !! length of the record
  INTEGER, INTENT(IN) :: npw
  !! number of plane waves
  INTEGER, INTENT(IN) :: iunit
  !! input iunit from where to read
  INTEGER, INTENT(IN) :: igmap(npwx)
  !! index for the mapping of the g
  COMPLEX(DP), INTENT(OUT) :: evc(npol*npwx,nbnd)
  !! wavefunction in g space
  !
  ! ... local variables
  !
  INTEGER                   :: ibnd, ig, is
  COMPLEX(DP), ALLOCATABLE  :: evc_r(:,:), dist_evc_r(:,:)

  allocate( evc_r( dffts%nnr, npol ) )
  allocate( dist_evc_r( dffts%nr1x*dffts%nr2x*dffts%nr3x , npol) )
  
  !
  ! Fourier transform it in reciprocal space
  !

  do ibnd=1,nbnd
     !
     ! read wfc in real space
     !
#if defined(__MPI)
     !
     ! ... First task reads and broadcasts ddrho to all pools
     !

     IF ( ionode ) &
          CALL davcio (dist_evc_r, lrec, iunit, (ik-1)*nbnd+ibnd, - 1)
     
     CALL mp_bcast( dist_evc_r, ionode_id, inter_pool_comm )
          !
     ! ... distributes ddrho between between the tasks of the pool
     !
     DO is = 1, npol
        !
        CALL scatter_grid ( dffts, dist_evc_r(:,is), evc_r(:,is) )
        !
     END DO
     !
     !     call mp_bcast( evc_r, ionode_id, inter_pool_comm )
#else    
     CALL davcio (evc_r, lrec, iunit, (ik-1)*nbnd+ibnd, - 1)
#endif
     
     call fwfft('Wave',evc_r(:,1),dffts)
     do ig = 1, npw
        evc (ig,ibnd) = evc_r (dffts%nl (igmap (ig) ), 1 )
     enddo
     
     IF (noncolin) THEN
        CALL fwfft ('Wave', evc_r(:,2), dffts)
        DO ig = 1, npw
           evc (ig+npwx,ibnd) = evc_r (dffts%nl(igmap(ig)),2)
        ENDDO
     ENDIF
     
  enddo

  deallocate( evc_r )
  deallocate( dist_evc_r )

end subroutine read_wfc_rspace_and_fwfft
!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! -----------------------------------------------------------------
! This program reads wavefunctions in G-space written by QE,
! re-writes then in real space
! Warning: The wfc is written out in real space on the smooth
! grid, as such it occupies much more disk space than in G-space.
!
! input: a namelist like 
! &inputpp
!   prefix='MgB2',
!   outdir='./tmp',
! /
! with "prefix" and "outdir" as in the scf/nscf/band calculation. 
! A file "prefix".wfc_r1 will be created in "outdir" with wfcs in real space
! The code prints on screen the dimension of the grid and of the wavefunctions

! Other namelist variables
! To select a subset of k-points and bands (by default, everything is written):
!        * first_k
!        * last_k
!        * first_band
!        * last_band
!
! Program written by Matteo Calandra.
! Modified by D. Ceresoli (2017)
! 
!-----------------------------------------------------------------------
subroutine wfck2r_ep()
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE save_ph,              ONLY : tmp_dir_save
  USE io_files,  ONLY : prefix, diropn
  USE wvfct,     ONLY : nbnd, npwx, et, wg
  USE klist,     ONLY : xk, nks, ngk, igk_k, wk
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_global, ONLY : mp_startup
  USE mp_images, ONLY : intra_image_comm
  USE mp_pools,  ONLY : npool
  USE wavefunctions, ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect, ONLY : ngm, g 
  USE noncollin_module, ONLY : npol, noncolin
  USE environment,ONLY : environment_start, environment_end
  USE fft_base,  only : dffts
  USE scatter_mod,  only : gather_grid
  USE fft_interfaces, ONLY : invfft
  USE ener, ONLY: efermi => ef
  USE pw_restart_new,ONLY : read_collected_wfc
USE el_phon,    ONLY : elph_nbnd_min,elph_nbnd_max
  USE dfile_star,    ONLY : dvscf_star
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: dirname,outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER            :: npw, iunitout,ios,ik,i,iuwfcr,lrwfcr,ibnd, ig, is
  LOGICAL            :: needwf= .TRUE., exst
  COMPLEX(DP), ALLOCATABLE :: evc_r(:,:), dist_evc_r(:,:)
  INTEGER :: first_k, last_k, first_band, last_band
  INTEGER, EXTERNAL :: find_free_unit


  !
  !


  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
!  IF ( ionode )  THEN
     !
     ! set defaults:
     first_k = 1
     last_k = nks
     write(*,*) 'nks:', nks
     first_band = elph_nbnd_min
     last_band = elph_nbnd_max
     !
     ! 
!  END IF

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
 ! CALL read_file_new ( needwf ) !GIO CAREFUL

  filename='wfc_r'
  write(6,*) 'filename              = ', trim(filename)
  iuwfcr=find_free_unit()
  lrwfcr = 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x * npol
  ! lrwfc = 2 * nbnd * npwx * npol
  write(6,*) 'dffts%nnr, npwx       =', dffts%nnr, npwx
 
  write(6,*) 'first_k, last_k       =', first_k, last_k
  write(6,*) 'first_band, last_band =', first_band, last_band
  write(6,*)

  write(6,*) 'length of wfc in real space/per band', (last_k-first_k+1)*lrwfcr*8
  write(6,*) 'length of wfc in k space', 2*(last_band-first_band+1)*npwx*nks*8

!
!define lrwfcr
!
  exst=.false.
      !outdir=trimcheck ('Rotated_DVSCF')
      !print*,dvscf_star%dir
 
  IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst,dvscf_star%dir)
 !CALL diropn (iuwfcr, filename, lrwfcr, exst,outdir)
  ALLOCATE ( evc_r(dffts%nnr,npol) )
  ALLOCATE ( dist_evc_r(dffts%nr1x*dffts%nr2x*dffts%nr3x,npol) )
  

  DO ik = first_k, last_k
     
     npw = ngk(ik)
      dirname=trimcheck ( TRIM(tmp_dir_save) // TRIM(prefix) // &
                                & '.save')
     CALL read_collected_wfc ( dirname, ik, evc )

     do ibnd = first_band, last_band
        !
        ! perform the fourier transform
        !
        evc_r = (0.d0, 0.d0)     
        do ig = 1, npw
           evc_r (dffts%nl (igk_k(ig,ik) ),1 ) = evc (ig,ibnd)
        enddo
        CALL invfft ('Wave', evc_r(:,1), dffts)
        IF (noncolin) THEN
           DO ig = 1, npw
              evc_r (dffts%nl(igk_k(ig,ik)),2) = evc (ig+npwx, ibnd)
           ENDDO
           CALL invfft ('Wave', evc_r(:,2), dffts)
        ENDIF

        dist_evc_r=(0.d0,0.d0)

#if defined (__MPI)
        DO is = 1, npol
           !
           CALL gather_grid( dffts, evc_r(:,is), dist_evc_r(:,is) )
           !
        END DO
#else
        dist_evc_r(1:dffts%nnr,:)=evc_r(1:dffts%nnr,:)
#endif

         !call davcio (dist_evc_r, lrwfcr, iuwfcr, (ik-1)*nbnd+ibnd, +1)
        if (ionode) call davcio (dist_evc_r, lrwfcr, iuwfcr, (ik-1)*nbnd+ibnd, +1)
     enddo
        
     !
     ! ... First task is the only task allowed to write the file
     !

  enddo

  if (ionode) close(iuwfcr)
  DEALLOCATE (evc_r)
  deallocate( dist_evc_r )
return
end subroutine wfck2r_ep

subroutine wfck2r_clean_files()
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE dfile_star,    ONLY : dvscf_star
  USE io_global, ONLY : ionode, stdout
  USE el_phon,         ONLY :  elph_mat
  USE io_files,      ONLY :  prefix
  !
  IMPLICIT NONE
  LOGICAL :: exst
  INTEGER :: un, ios
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER, EXTERNAL :: find_free_unit

  IF(elph_mat.and.ionode)then

       ! ...  search for file wfc_r and delete it
       filename=trim ( TRIM(dvscf_star%dir) // trim(prefix) //'.wfc_r1')
       WRITE(stdout,'(5x,"Deleting: ",a)')filename
       INQUIRE( FILE=TRIM(filename), EXIST=exst )
       IF( exst ) THEN
          un = find_free_unit()
          OPEN( UNIT=un, FILE=TRIM(filename), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=un, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: ",a," file could not be deleted")')filename
          END IF
        ELSE
             WRITE(stdout,'(5x,"Remark: ",a," file not exist?!")')filename
       END IF
  endif
end subroutine wfck2r_clean_files
