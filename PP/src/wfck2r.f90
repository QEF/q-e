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
! To create a file that is readable by octave (false by default):
!        * loctave=.true.
!
! Program written by Matteo Calandra.
! Modified by D. Ceresoli (2017)
! 
!-----------------------------------------------------------------------
PROGRAM wfck2r
  !-----------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : prefix, tmp_dir, diropn, restart_dir
  USE wvfct,            ONLY : nbnd, npwx, et, wg
  USE klist,            ONLY : xk, nks, ngk, igk_k, wk
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast, mp_barrier
  USE mp_global,        ONLY : mp_startup
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : npool
  USE wavefunctions,    ONLY : evc
  USE io_files,         ONLY : nwordwfc, iunwfc
  USE gvect,            ONLY : ngm, g
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : npol, noncolin
  USE environment,      ONLY : environment_start, environment_end
  USE fft_base,         ONLY : dffts
  USE scatter_mod,      ONLY : gather_grid
  USE fft_interfaces,   ONLY : invfft
  USE ener,             ONLY : efermi => ef
  USE pw_restart_new,   ONLY : read_collected_wfc
  USE constants,        ONLY : rytoev
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER            :: npw, iunitout,ios,ik,i,iuwfcr,lrwfcr,ibnd, ig, is
  LOGICAL            :: needwf= .TRUE., exst
  COMPLEX(DP), ALLOCATABLE :: evc_r(:,:), dist_evc_r(:,:)
  COMPLEX(DP) :: val
  INTEGER :: first_k, last_k, first_band, last_band
  INTEGER :: nevery(3), i1, i2, i3
  LOGICAL :: loctave

  NAMELIST / inputpp / outdir, prefix, first_k, last_k, first_band, last_band, loctave, nevery


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
     ! set defaults:
     first_k = 0
     last_k = 0
     first_band = 0
     last_band = 0
     nevery(1:3) = 1
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

  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  CALL mp_bcast( first_k, ionode_id, intra_image_comm )
  CALL mp_bcast( last_k, ionode_id, intra_image_comm )
  CALL mp_bcast( first_band, ionode_id, intra_image_comm )
  CALL mp_bcast( last_band, ionode_id, intra_image_comm )
  CALL mp_bcast( loctave, ionode_id, intra_image_comm )
  CALL mp_bcast( nevery, ionode_id, intra_image_comm )

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file_new ( needwf )

  filename='wfc_r'
  write(6,*) 'filename              = ', trim(filename)
  iuwfcr=877
  lrwfcr = 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x * npol
  ! lrwfc = 2 * nbnd * npwx * npol
  write(6,*) 'dffts%nnr, npwx       =', dffts%nnr, npwx
 
  if (first_k <= 0) first_k = 1 
  if (last_k <= 0) last_k = nks
  if (first_band <= 0) first_band = 1 
  if (last_band <= 0) last_band = nbnd
  write(6,*) 'first_k, last_k       =', first_k, last_k
  write(6,*) 'first_band, last_band =', first_band, last_band
  if (mod(dffts%nr1x,nevery(1)) /= 0) call errore("wfck2r", "nr1x not multiple of nevery(1)", 1)
  if (mod(dffts%nr2x,nevery(2)) /= 0) call errore("wfck2r", "nr2x not multiple of nevery(2)", 1)
  if (mod(dffts%nr3x,nevery(3)) /= 0) call errore("wfck2r", "nr3x not multiple of nevery(3)", 1)
  write(6,*)

  write(6,*) 'length of wfc in real space/per band', (last_k-first_k+1)*lrwfcr*8
  write(6,*) 'length of wfc in k space', 2*(last_band-first_band+1)*npwx*nks*8

!
!define lrwfcr
!
  exst=.false.
  IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst)
  IF (loctave .and. ionode) then
     open(unit=iuwfcr+1, file='wfck2r.oct', status='unknown', form='formatted')
     write(iuwfcr+1,'(A)') '# created by wfck2r.x of Quantum-Espresso'
     ! Fermi energy
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,E20.10,//)') 'efermi', efermi*rytoev
     ! k-points
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,I4,//)') 'nkpoints', (last_k-first_k+1)
     write(iuwfcr+1,'("# name: ",A,/,"# type: matrix")') 'xk'
     write(iuwfcr+1,'("# rows: ",I5)') last_k-first_k+1
     write(iuwfcr+1,'("# columns: ",I5)') 3
     do ik = first_k, last_k
        write(iuwfcr+1,'(E20.12)') (xk(i,ik), i=1,3)
     enddo
     write(iuwfcr+1,*)
     write(iuwfcr+1,'("# name: ",A,/,"# type: matrix")') 'wk'
     write(iuwfcr+1,'("# rows: ",I5)') last_k-first_k+1
     write(iuwfcr+1,'("# columns: ",I5)') 1
     do ik = first_k, last_k
        write(iuwfcr+1,'(E20.12)') wk(ik)
     enddo
     write(iuwfcr+1,*)
     ! bands
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,I4,//)') 'nbands', (last_band-first_band+1)
     write(iuwfcr+1,'("# name: ",A,/,"# type: matrix")') 'eigs'
     write(iuwfcr+1,'("# rows: ",I5)') last_k-first_k+1
     write(iuwfcr+1,'("# columns: ",I5)') last_band-first_band+1
     do ik = first_k, last_k
        write(iuwfcr+1,'(E20.12)') (et(i,ik)*rytoev, i=first_band,last_band)
     enddo
     write(iuwfcr+1,*)
     write(iuwfcr+1,'("# name: ",A,/,"# type: matrix")') 'occup'
     write(iuwfcr+1,'("# rows: ",I5)') last_k-first_k+1
     write(iuwfcr+1,'("# columns: ",I5)') last_band-first_band+1
     do ik = first_k, last_k
        write(iuwfcr+1,'(F20.12)') (wg(i,ik)/wk(ik), i=first_band,last_band)
     enddo
     write(iuwfcr+1,*)
     ! FFT mesh
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,I3,//)') 'nr1x', dffts%nr1x / nevery(1)
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,I3,//)') 'nr2x', dffts%nr2x / nevery(2)
     write(iuwfcr+1,'("# name: ",A,/,"# type: scalar",/,I3,//)') 'nr3x', dffts%nr3x / nevery(3)
     ! wavefunctions
     write(iuwfcr+1,'("# name: ",A,/,"# type: complex matrix")') 'unkr'
     write(iuwfcr+1,'("# ndims: 5")')
     write(iuwfcr+1,'(5I10)') dffts%nr1x/nevery(1), dffts%nr2x/nevery(2), dffts%nr3x/nevery(3), &
                              last_band-first_band+1, last_k-first_k+1
  endif

  ALLOCATE ( evc_r(dffts%nnr,npol) )
  ALLOCATE ( dist_evc_r(dffts%nr1x*dffts%nr2x*dffts%nr3x,npol) )
  

  DO ik = first_k, last_k
     
     npw = ngk(ik)
     CALL read_collected_wfc ( restart_dir(), ik, evc )

     do ibnd = first_band, last_band
        !
        ! perform the fourier transform
        !
        evc_r = (0.d0, 0.d0)     
        do ig = 1, npw
           evc_r(dffts%nl(igk_k(ig,ik)),1) = evc(ig,ibnd)
        enddo
        if (gamma_only) then
           do ig = 1, npw
              evc_r(dffts%nlm(igk_k(ig,ik)),1) = conjg(evc(ig,ibnd))
           enddo
        endif
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

        if (ionode .and. .not. loctave) call davcio (dist_evc_r, lrwfcr, iuwfcr, (ik-1)*nbnd+ibnd, +1)
        if (ionode .and. loctave) then
           do i3 = 1, dffts%nr3x, nevery(3)
           do i2 = 1, dffts%nr2x, nevery(2)
           do i1 = 1, dffts%nr1x, nevery(1)
               !val = local_average()
               val = dist_evc_r(i1 + (i2-1)*dffts%nr1x + (i3-1)*dffts%nr1x*dffts%nr2x,1)
               write(iuwfcr+1,'("(",E20.12,",",E20.12,")")') val
           enddo
           enddo
           enddo
        endif
     enddo
        
     !
     ! ... First task is the only task allowed to write the file
     !

  enddo

  if (ionode .and. .not. loctave) close(iuwfcr)
  if (loctave .and. ionode) close(iuwfcr+1)
  DEALLOCATE (evc_r)

  CALL environment_end ( 'WFCK2R' )

  CALL stop_pp
  STOP

end PROGRAM wfck2r
