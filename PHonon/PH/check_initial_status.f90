!
! Copyright (C) 2012-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE check_initial_status(auxdyn)
  !-----------------------------------------------------------------------
  !
  ! This routine checks the initial status of the phonon run and sets
  ! the variables that control the run, dealing with the image
  ! and GRID parallelization features of the phonon code.
  ! 
  ! The size of the grid is determined by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  !
  ! nfs : the number of imaginary frequencies
  ! fiu : which frequencies 
  !
  ! The flags that control which tensors to calculate
  !
  ! In a recover calculation the q grid variables are already known, 
  ! read from file in phq_readin. In a calculation starting from
  ! scratch this routine sets them. The frequencies variables and the
  ! tensors flags are read from input. 
  ! The amount of work to do for each representation of each q
  ! point depends on the size of the representation and the 
  ! order of the small group of q. In a recover calculation
  ! these information are on file, when recover=.false. this
  ! routine writes the modes and their degeneration on files 
  ! and calculates the order of the small group of q. The following
  ! variables are set
  !
  ! irr_iq : for each q point how many irreducible representations
  ! npert_irr_iq : how many perturbation per representation and per q
  ! nsymq_iq : the order of the small group of q for each q
  !
  ! The following variables are set by this routine on the basis of
  ! start_irr, last_irr, start_iq, last_iq, OR of modenum, OR of ifat and 
  ! atomo:
  !
  ! comp_iq : =.TRUE. if the q point is calculated in this run
  ! comp_irr_iq : =.TRUE. if the representation is calculated in this run
  ! comp_iu : =.TRUE. if this frequency is calculated in this run
  !                   NB: start_iu, last_iu is not yet programmed
  ! 
  ! After knowing this info the routine divides the total work among
  ! the images (when nimage > 1) INDEPENDENTLY of what has been already
  ! calculated and is available on file.
  !
  ! Then, when recover=.true., the routine looks on files for pieces
  ! already calculated and sets the array
  !
  ! done_irr_iq : =.TRUE. if the representation has been already calculated
  ! done_iq : =.TRUE. if the q point has been already calculated
  ! done_iu : =.TRUE. if already calculated
  ! done_bands_iq : .TRUE. if the bands for the q point are on file.
  !
  ! If recover=.false. all these array are initialized to .false.
  !
  ! Finally this routine creates a file fildyn0 and writes the q mesh, if
  ! this file is not present in the current directory, or if recover=.false..
  ! It also creates a directory for each q inside outdir/_ph# 
  ! if this directory does not exist and lqdir=.true.
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : modenum
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE scf,             ONLY : rho
  USE disp,            ONLY : nqs, x_q, comp_iq, nq1, nq2, nq3, &
                              done_iq, lgamma_iq
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma
  USE output,          ONLY : fildyn
  USE control_ph,      ONLY : ldisp, recover, where_rec, rec_code, &
                              start_q, last_q, current_iq, tmp_dir_ph, &
                              ext_recover, ext_restart, tmp_dir_phq, lqdir, &
                              start_irr, last_irr, newgrid, qplot, &
                              done_zeu, done_start_zstar, done_epsil, &
                              done_zue, with_ext_images, always_run
  USE save_ph,         ONLY : tmp_dir_save
  USE units_ph,        ONLY : iudyn
  USE ph_restart,      ONLY : check_directory_phsave, check_available_bands,&
                              allocate_grid_variables, ph_writefile
  USE freq_ph,         ONLY : current_iu
  USE io_rho_xml,      ONLY : write_scf
  USE mp_images,       ONLY : nimage, intra_image_comm
  USE io_global,       ONLY : ionode, ionode_id
  USE io_files,        ONLY : prefix
  USE mp,              ONLY : mp_bcast
  USE xml_io_base,     ONLY : create_directory
  USE mp_global,       ONLY : mp_global_end
  USE el_phon,         ONLY : elph_mat
  ! YAMBO >
  USE YAMBO,           ONLY : elph_yambo,dvscf_yambo
  ! YAMBO <
  !
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn, filename
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER :: iq, iq_start, ierr
  !
  tmp_dir=tmp_dir_ph
  !
  ! If this not a recover run, we generate the q mesh. Otherwise at this
  ! point the code has read the q mesh from the files contained in 
  ! prefix.phsave
  !
  IF (.NOT.recover) THEN
     !
     ! recover file not found or not looked for
     !
     current_iu=1
     current_iq=1
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        IF(elph_mat) then
           CALL q_points_wannier()
        ELSE
           IF (.NOT. qplot) CALL q_points()
        END IF
        !
        ! YAMBO >
     ELSE IF (.NOT.elph_yambo .AND. .NOT. dvscf_yambo) then
        ! YAMBO <
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        ALLOCATE(lgamma_iq(1))
        x_q(:,1)=xq(:)
        lgamma_iq(1)=lgamma
        !
     END IF
     !
     !   Save the mesh of q and the control flags on file
     !
     CALL ph_writefile('init',0,0,ierr)
     !
     !   Initialize the representations and write them on file.
     !
     CALL init_representations()
     !
     IF ((start_irr==0).AND.(last_irr==0)) THEN
        where_rec='init_rep..'
        rec_code=-50
        current_iq=1
        current_iu=1
        CALL ph_writefile('status_ph',current_iq,0,ierr)
        CALL clean_pw(.FALSE.)
        CALL close_files(.FALSE.)
        CALL mp_global_end()
        STOP
     ENDIF
  ENDIF

  IF (last_q<1.or.last_q>nqs) last_q=nqs
  IF (start_q<1.or.start_q>last_q) call errore('check_initial_status', &
     'wrong start_q',1)
!
!  now we allocate the variables needed to describe the grid
!
  CALL allocate_grid_variables()
!
!  This routine assumes that the modes are on file, either written by 
!  init_representation or written by a previous run. It takes care
!  of dealing with start_irr, last_irr flags and ifat or modenum
!  restricted  computation, moreover it sets the size of each 
!  representation and the size of the small group of q for each point.
!
  CALL initialize_grid_variables()
!
! If there are more than one image, divide the work among the images
!
  IF (nimage > 1 .AND. .NOT. with_ext_images) CALL image_q_irr()
!
  IF (recover) THEN
!
! ... Checking the status of the calculation
!
!  sets which q point and representations have been already calculated
!
     CALL check_directory_phsave()
!
!  If a recover or a restart file exists the first q point is the current one.
!
     IF ((.NOT.lgamma_iq(current_iq).OR. newgrid).AND.lqdir) THEN
        tmp_dir_phq= TRIM (tmp_dir_ph) //TRIM(prefix)//&
                          & '.q_' // TRIM(int_to_char(current_iq))//'/'
        tmp_dir=tmp_dir_phq
        CALL check_restart_recover(ext_recover, ext_restart)
        tmp_dir=tmp_dir_ph
     ELSE
        CALL check_restart_recover(ext_recover, ext_restart)
     ENDIF
     IF (.NOT.ext_recover.AND..NOT.ext_restart) THEN
        current_iq=start_q
     ELSE
!
!   Check that the representations from start_q to current_iq have been done
!
        DO iq=start_q, current_iq-1
           IF (comp_iq(iq) .AND. .NOT.done_iq(iq)) &
              CALL errore('check_initial_status',&
                      & 'recover file found, change in start_q not allowed',1)
           comp_iq(iq)=.FALSE.
        ENDDO
     ENDIF
     iq_start=current_iq
!
!  check which bands are available and set the array done_bands
!
    CALL check_available_bands()
!
! write the information on output
    IF (iq_start<=last_q.AND.iq_start>0) THEN
        WRITE(stdout, &
            '(5x,i4," /",i4," q-points for this run, from", i3,&
               & " to", i3,":")') last_q-iq_start+1, nqs, iq_start, last_q
        WRITE(stdout, '(5x,"  N       xq(1)         xq(2)         xq(3) " )')
        DO iq = 1, nqs
           WRITE(stdout, '(5x,i3, 3f14.9,l6)') iq, x_q(1,iq), x_q(2,iq), &
                             x_q(3,iq)
        END DO
        WRITE(stdout, *)
     ELSEIF (iq_start>last_q) THEN
        WRITE(stdout, &
            '(5x,"Starting q",i4," larger than total number of q points", i4, &
               & " or of last q ", i3)') iq_start, nqs, last_q
     ELSEIF (iq_start<0) THEN
        CALL errore('check_initial_status','wrong iq_start',1)
     ENDIF
  ELSE
     done_zeu=.FALSE.
     done_start_zstar=.FALSE.
     done_epsil=.FALSE.
     done_zue=.FALSE.
  ENDIF
  !
  !  Create a new directory where the ph variables are saved and copy
  !  the charge density there.
!!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
  IF (acfdt_is_active) THEN
     ! ACFDT -test always write rho on file
     IF (acfdt_num_der) THEN
        CALL write_scf( rho, nspin )
     ELSE 
        IF ((ldisp.OR..NOT.lgamma.OR.modenum/=0).AND.(.NOT.lqdir)) &
                                           CALL write_scf( rho, nspin )
     ENDIF   
  ELSE  
     ! this is the standard treatment
     IF ( ( ( ldisp.OR..NOT.lgamma .OR. modenum/=0 ) .AND. (.NOT.lqdir) ) &
          .OR. newgrid .OR. always_run ) CALL write_scf( rho, nspin )
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! END OF ACFDT TEST !!!!!!!!!!!!!!!!
  !
  !  Write the file fildyn0 with the mesh of q points. This file is used
  !  by postprocessing programs such as q2r.x and we write it again if
  !  it is not found in the running directory.
  !
  filename=TRIM(fildyn)//'0'
  ierr=0 
  IF (ionode.and..NOT.elph_mat) THEN
     INQUIRE (FILE = TRIM(filename), EXIST = exst)
     ierr=0
     IF ((.NOT. exst .OR. .NOT. recover).AND.ldisp) THEN
        iudyn=26
        OPEN (unit=iudyn, file=TRIM(filename), status='unknown', iostat=ierr)
        IF ( ierr == 0 ) THEN
           WRITE (iudyn, '(3i4)' ) nq1, nq2, nq3
           WRITE (iudyn, '( i4)' ) nqs
           DO  iq = 1, nqs
               WRITE (iudyn, '(3e24.15)') x_q(1,iq), x_q(2,iq), x_q(3,iq)
           END DO
           CLOSE (unit=iudyn)
        ENDIF
     ENDIF
  END IF
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  CALL errore ('check_initial_status','cannot open file ' // TRIM(filename),&
                abs(ierr))
  !
  !  The following commands deal with the flag lqdir=.true. In this case
  !  each q point works on a different directory. We create the directories
  !  if they do not exist and copy the self consistent charge density
  !  there.
  !
  DO iq = 1,nqs
     IF (.NOT.comp_iq(iq)) CYCLE
     lgamma = lgamma_iq(iq) 
     !
     ! ... each q /= gamma works on a different directory. We create them
     ! here and copy the charge density inside
     !
     IF ((.NOT.lgamma.OR. newgrid).AND.lqdir) THEN
        tmp_dir_phq= TRIM (tmp_dir_ph) //TRIM(prefix)//&
                          & '.q_' // TRIM(int_to_char(iq))//'/'
        filename=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/charge-density.dat'
        IF (ionode) inquire (file =TRIM(filename), exist = exst)
        !
        CALL mp_bcast( exst, ionode_id, intra_image_comm )
        !
        IF (.NOT. exst) THEN
           CALL create_directory( tmp_dir_phq )
           tmp_dir=tmp_dir_phq
           CALL write_scf( rho, nspin )
           tmp_dir=tmp_dir_save
        ENDIF
     ENDIF
  ENDDO
  !
  auxdyn = fildyn
  RETURN
  END SUBROUTINE check_initial_status

  SUBROUTINE image_q_irr()
!
!  This routine is an example of the load balancing among images.
!  It decides which image makes which q and which irreducible representation
!  The algorithm at the moment is straightforward. Possibly better
!  methods could be found.
!  It receives as input:
!  nsym  : the dimension of the point group
!  nsymq_iq  : the dimension of the small group of q for each q
!  irr_iq : the number of irreps for each q
!  npert_irr_iq : for each q and each irrep its dimension
!  It provides as output the two arrays
!  comp_iq : if this q has to be calculated by the present image
!  comp_irr_iq : for each q the array to be copied into comp_irr

   USE ions_base, ONLY : nat
   USE disp, ONLY : comp_iq, nqs, nq1, nq2, nq3
   USE grid_irr_iq, ONLY : irr_iq, npert_irr_iq, comp_irr_iq, nsymq_iq
   USE control_ph, ONLY : start_q, last_q
   USE io_global,  ONLY : stdout
   USE mp_images,  ONLY : nimage, my_image_id
   USE symm_base,  ONLY : nsym

   IMPLICIT NONE
   INTEGER :: total_work,  &  ! total amount of work to do
              total_nrapp, &  ! total number of representations
              work_per_image  ! approximate minimum work per image

   INTEGER, ALLOCATABLE :: image_iq(:,:), work(:)
   INTEGER :: iq, irr, image, work_so_far, actual_diff, diff_for_next
   CHARACTER(LEN=256) :: string
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   ALLOCATE (image_iq(0:3*nat,nqs))
   ALLOCATE (work(0:nimage-1))

   total_work=0
   total_nrapp=0
   DO iq = start_q, last_q
      DO irr = 1, irr_iq(iq)
         IF (comp_irr_iq(irr,iq)) THEN
            total_work = total_work + npert_irr_iq(irr, iq) * nsym / nsymq_iq(iq)
            IF (irr==1) total_work = total_work + nsym / nsymq_iq(iq)
            total_nrapp = total_nrapp + 1
         ENDIF
      END DO
   END DO
   IF (nimage > total_nrapp) &
      CALL errore('image_q_irr','some images have no rapp', 1)

   work_per_image = total_work / nimage
!
!  If nimage=total_nrapp we put one representation per image
!  No load balancing is possible. Otherwise we try to minimize the number of
!  different q per image doing all representations of a given q until
!  the work becomes too large.
!  The initialization is done by the image with the first representation of
!  each q point.
!
   image=0
   work=0
   work_so_far=0
   DO iq = start_q, last_q
      DO irr = 1, irr_iq(iq)
       IF (comp_irr_iq(irr,iq)) THEN
         image_iq(irr,iq) = image
         work(image)=work(image) + npert_irr_iq(irr, iq) * nsym / nsymq_iq(iq)
         work_so_far=work_so_far + npert_irr_iq(irr, iq) * nsym / nsymq_iq(iq)
         IF (irr==1) THEN
            image_iq(0,iq)=image
            work(image)=work(image) + nsym / nsymq_iq(iq)
            work_so_far=work_so_far + nsym / nsymq_iq(iq)
         ENDIF

!
!  The logic is the following. We know how much work the current image
!  has already accumulated and we calculate how far it is from the target.
!  Note that actual_diff is a positive number in the usual case in which
!  we are below the target. Then we calculate the work that the current
!  image would do if we would give it the next representation. If the work is
!  still below the target, diff_for_next is negative and we give the
!  representation to the current image. If the work is above the target,
!  we give it to the current image only if its distance from the target
!  is less than actual_diff.
!
         actual_diff=-work(image)+work_per_image
         IF (irr<irr_iq(iq)) THEN
            diff_for_next= work(image)+npert_irr_iq(irr+1, iq)*nsym/nsymq_iq(iq) &
                           - work_per_image
         ELSEIF (irr==irr_iq(iq).and.iq<last_q) THEN
            diff_for_next= work(image)+npert_irr_iq(1, iq+1)* &
                       nsym/nsymq_iq(iq+1) + nsym/nsymq_iq(iq+1)-work_per_image
         ELSE
            diff_for_next=0
         ENDIF

         IF ((nimage==total_nrapp.OR.diff_for_next>actual_diff).AND. &
                              (image < nimage-1)) THEN
            work_per_image= (total_work-work_so_far) / (nimage-image-1)
            image=image+1
         ENDIF
       ENDIF
      ENDDO
   ENDDO
!
!  Here we actually distribute the work. This image makes only
!  the representations calculated before.
!
   DO iq = start_q, last_q
      DO irr = 0, irr_iq(iq)
         IF (image_iq(irr,iq)/=my_image_id ) THEN
            comp_irr_iq(irr,iq)=.FALSE.
         ENDIF
      ENDDO
   ENDDO

   comp_iq = .FALSE.
   DO iq = start_q, last_q
      DO irr = 0, irr_iq(iq)
         IF (comp_irr_iq(irr,iq).AND..NOT.comp_iq(iq)) THEN
            comp_iq(iq)=.TRUE.
         ENDIF
      ENDDO
   ENDDO

   WRITE(stdout, &
            '(/,5x," Image parallelization. There are", i3,&
               & " images", " and ", i5, " representations")') nimage, &
               total_nrapp
   WRITE(stdout, &
            '(5x," The estimated total work is ", i5,&
               & " self-consistent (scf) runs")') total_work

   WRITE(stdout, '(5x," I am image number ",i5," and my work is about",i5, &
                      &  " scf runs. I calculate: ")') &
                        my_image_id, work(my_image_id)

   DO iq = 1, nqs
      IF (comp_iq(iq)) THEN
         WRITE(stdout, '(5x," q point number ", i5, ", representations:")') iq
         string=' '
         DO irr=0, irr_iq(iq)
            IF (comp_irr_iq(irr, iq)) &
                string=TRIM(string) // " " // TRIM(int_to_char(irr))
         ENDDO
         WRITE(stdout,'(6x,A)') TRIM(string)
      ENDIF
   ENDDO

   DEALLOCATE(image_iq)
   DEALLOCATE(work)
   RETURN
   END SUBROUTINE image_q_irr

   SUBROUTINE collect_grid_files()
   !
   !  This subroutine collects all the xml files contained in different
   !  directories and created by the diffent images in the phsave directory
   !  of the image 0
   !
   USE io_files,  ONLY : tmp_dir, prefix
   USE control_ph, ONLY : tmp_dir_ph
   USE save_ph,   ONLY : tmp_dir_save
   USE disp,      ONLY : nqs
   USE grid_irr_iq,  ONLY : comp_irr_iq, irr_iq
   USE el_phon,     ONLY : elph
   USE wrappers,  ONLY : f_copy
   USE mp,        ONLY : mp_barrier
   USE mp_images, ONLY : my_image_id, nimage, intra_image_comm
   USE io_global, ONLY : stdout, ionode

   IMPLICIT NONE

   INTEGER :: iq, irr, ios
   LOGICAL :: exst
   CHARACTER(LEN=256) :: file_input, file_output
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   CALL mp_barrier(intra_image_comm)
   IF (nimage == 1) RETURN
   IF (my_image_id==0) RETURN

   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         IF (comp_irr_iq(irr,iq).and.ionode) THEN
            file_input=TRIM( tmp_dir_ph ) // &
                    & TRIM( prefix ) // '.phsave/dynmat.'  &
                    &  // TRIM(int_to_char(iq))&
                    &  // '.' // TRIM(int_to_char(irr)) // '.xml'

            file_output=TRIM( tmp_dir_save ) // '/_ph0/' &
                    &    // TRIM( prefix ) // '.phsave/dynmat.' &
                    &    // TRIM(int_to_char(iq))  &
                    &    // '.' // TRIM(int_to_char(irr)) // '.xml'

            INQUIRE (FILE = TRIM(file_input), EXIST = exst)
            IF (exst) ios = f_copy(file_input, file_output)
            IF ( elph .AND. irr>0 ) THEN

               file_input=TRIM( tmp_dir_ph ) // &
                    & TRIM( prefix ) // '.phsave/elph.'  &
                    &  // TRIM(int_to_char(iq))&
                    &  // '.' // TRIM(int_to_char(irr)) // '.xml'

               file_output=TRIM( tmp_dir_save ) // '/_ph0/' // &
                    &   TRIM( prefix ) // '.phsave/elph.' &
                    &    // TRIM(int_to_char(iq))  &
                    &    // '.' // TRIM(int_to_char(irr)) // '.xml'

               INQUIRE (FILE = TRIM(file_input), EXIST = exst)
               IF (exst) ios = f_copy(file_input, file_output)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
   RETURN
   END SUBROUTINE collect_grid_files
