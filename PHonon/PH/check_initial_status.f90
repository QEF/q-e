!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE check_initial_status(auxdyn)
  !-----------------------------------------------------------------------
  !
  ! This routine checks the initial status of the phonon run and prepares
  ! the control of the dispersion calculation. The grid is determined
  ! by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  ! comp_iq : =1 if this q point is calculated in this run, 0 otherwise
  ! done_iq : =1 if already calculated, 0 otherwise
  ! rep_iq  : for each q point how many irreducible representations
  ! done_rep_iq : =1 if the representation has been already calculated
  ! nfs : the number of imaginary frequencies
  ! The last four variables are calculated only at gamma if fpol is .true.
  ! If recover is true this routine checks also that the control parameters
  ! read in input match the values given in the recover file.
  ! Finally it sets the variable:
  ! done_bands : if true the bands have been already calculated for this q
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : modenum
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE scf,             ONLY : rho
  USE disp,            ONLY : nqs, x_q, comp_iq, comp_irr_iq
  USE qpoint,          ONLY : xq
  USE output,          ONLY : fildyn
  USE control_ph,      ONLY : ldisp, recover, done_bands,  &
                              start_q, last_q, current_iq, tmp_dir_ph, lgamma, &
                              ext_recover, ext_restart, tmp_dir_phq, lqdir, &
                              start_irr, last_irr, newgrid
  USE save_ph,         ONLY : tmp_dir_save
  USE ph_restart,      ONLY : ph_readfile, check_status_run, init_status_run, &
                              ph_writefile
  USE save_ph,         ONLY : save_ph_input_variables
  USE io_rho_xml,      ONLY : write_rho
  USE mp_global,       ONLY : nimage, my_image_id, intra_image_comm
  USE io_global,       ONLY : ionode, ionode_id
  USE io_files,        ONLY : prefix
  USE mp,              ONLY : mp_bcast
  USE xml_io_base,     ONLY : create_directory
  USE mp_global,       ONLY : mp_global_end
  USE el_phon,         ONLY : elph_mat
  !
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn, filename
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER :: iq, iq_start, ierr
  INTEGER :: iu
  !
  ! Initialize local variables
  !
  tmp_dir=tmp_dir_ph
  !
  ! ... Checking the status of the calculation
  !
  IF (recover) THEN
!
!  check if a recover file exists. In this case the first q point is
!  the current one.
!
     IF (.NOT.ext_recover.AND..NOT.ext_restart) THEN
        iq_start=start_q
        done_bands=.FALSE.
     ELSE
        iq_start=current_iq
     ENDIF
!
!  check which representation files are available on the disk and
!  sets which q points and representations have been already calculated
!
     CALL check_status_run()
!
! write the information on output
!
     IF (last_q<1.OR.last_q>nqs) last_q=nqs
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
  ENDIF
  !
  !  Create a new directory where the ph variables are saved and copy
  !  the charge density there.
!!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
  IF (acfdt_is_active) THEN
     ! ACFDT -test always write rho on file
     IF (acfdt_num_der) THEN
        CALL write_rho( rho, nspin )
     ELSE 
        IF ((ldisp.OR..NOT.lgamma.OR.modenum/=0).AND.(.NOT.lqdir)) &
                                           CALL write_rho( rho, nspin )
     ENDIF   
  ELSE  
     ! this is the standard treatment
     IF ( ( ( ldisp .OR. .NOT.lgamma .OR. modenum/=0 ) .AND. (.NOT.lqdir) ) &
          .OR. newgrid ) CALL write_rho( rho, nspin )
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! END OF ACFDT TEST !!!!!!!!!!!!!!!!
  !
  CALL save_ph_input_variables()
  !
  IF (.NOT.recover) THEN
     !
     ! recover file not found or not looked for
     !
     done_bands=.FALSE.
     iq_start=start_q
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        IF(elph_mat) then
           CALL q_points_wannier()
        ELSE
           CALL q_points()
        END IF

        IF (last_q<1.or.last_q>nqs) last_q=nqs
        !
     ELSE
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
        !
     END IF
     !
     ! This routine initialize the grid control in order to
     ! calculate all q and all representations. The representations are
     ! written on file and read again by phq_setup.
     !
     CALL init_status_run()
     CALL init_representations()
     IF ((start_irr==0).AND.(last_irr==0)) THEN
        CALL ph_writefile('init',0)
        CALL clean_pw(.FALSE.)
        CALL close_files(.FALSE.)
        CALL mp_global_end()
        STOP
     ENDIF
     !
  END IF
  !
  !  Set the q points to calculate. If there is the recover file, start from
  !  the q point of the recover file.
  !
  IF (nimage==1) THEN
     comp_iq=0
     DO iq=iq_start,last_q
        comp_iq(iq)=1
     ENDDO
  ELSE
     CALL image_q_irr(iq_start)
  ENDIF
  !
  DO iq=1,nqs
     IF (comp_iq(iq).ne.1) CYCLE
     lgamma = ( x_q(1,iq) == 0.D0 .AND. x_q(2,iq) == 0.D0 .AND. &
                x_q(3,iq) == 0.D0 )
     !
     ! ... each q /= gamma works on a different directory. We create them
     ! here and copy the charge density inside
     !
     IF ((.NOT.lgamma.OR. newgrid).AND.lqdir) THEN
        tmp_dir_phq= TRIM (tmp_dir_ph) //TRIM(prefix)//&
                          & '_q' // TRIM(int_to_char(iq))//'/'
        filename=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/charge-density.dat'
        IF (ionode) inquire (file =TRIM(filename), exist = exst)
        !
        CALL mp_bcast( exst, ionode_id, intra_image_comm )
        !
        IF (.NOT. exst) THEN
           CALL create_directory( tmp_dir_phq )
           tmp_dir=tmp_dir_phq
           CALL write_rho( rho, nspin )
           tmp_dir=tmp_dir_save
        ENDIF
     ENDIF
  ENDDO
  !
  auxdyn = fildyn
  !
  RETURN
  END SUBROUTINE check_initial_status

  SUBROUTINE image_q_irr(iq_start)
!
!  This routine is an example of the load balancing among images.
!  It decides which image makes which q and which irreducible representations
!  The algorithm at the moment is straightforward. Possibly better
!  methods could be found.
!  It receives as input:
!  nsym  : the dimension of the point group
!  nsymq_iq  : the dimension of the small group of q for each q
!  rep_iq : the number of representation for each q
!  npert_iq : for each q and each irrep its dimension
!  It provides as output the two arrays
!  comp_iq : if this q has to be calculated by the present image
!  comp_irr_iq : for each q the array to be copied into comp_irr

   USE ions_base, ONLY : nat
   USE disp, ONLY : rep_iq, npert_iq, comp_iq, comp_irr_iq, nqs, &
                    nq1, nq2, nq3, nsymq_iq
   USE control_ph, ONLY : start_q, last_q
   USE io_global,  ONLY : stdout
   USE mp_global,  ONLY : nimage, my_image_id
   USE symm_base,  ONLY : nsym

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: iq_start ! the calculation start from this q.
                               ! It can be different from start_q in a
                               ! recovered run. The division of the work
                               ! is done from start_q to last_q.
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
      DO irr = 1, rep_iq(iq)
         total_work = total_work + npert_iq(irr, iq) * nsym / nsymq_iq(iq)
         IF (irr==1) total_work = total_work + nsym / nsymq_iq(iq)
         total_nrapp = total_nrapp + 1
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
      DO irr = 1, rep_iq(iq)
         image_iq(irr,iq) = image
         work(image)=work(image) + npert_iq(irr, iq) * nsym / nsymq_iq(iq)
         work_so_far=work_so_far + npert_iq(irr, iq) * nsym / nsymq_iq(iq)
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
         IF (irr<rep_iq(iq)) THEN
            diff_for_next= work(image)+npert_iq(irr+1, iq)*nsym/nsymq_iq(iq) &
                           - work_per_image
         ELSEIF (irr==rep_iq(iq).and.iq<last_q) THEN
            diff_for_next= work(image)+npert_iq(1, iq+1)* &
                       nsym/nsymq_iq(iq+1) + nsym/nsymq_iq(iq+1)-work_per_image
         ELSE
            diff_for_next=0
         ENDIF

         IF ((nimage==total_nrapp.OR.diff_for_next>actual_diff).AND. &
                              (image < nimage-1)) THEN
            work_per_image= (total_work-work_so_far) / (nimage-image-1)
            image=image+1
         ENDIF
      ENDDO
   ENDDO
!
!  Here we actually distribute the work. This image makes only
!  the representations calculated before.
!
   comp_iq = 0
   comp_irr_iq=0
   DO iq = iq_start, last_q
      DO irr = 0, rep_iq(iq)
         IF (image_iq(irr,iq)==my_image_id ) THEN
            comp_iq(iq)=1
            comp_irr_iq(irr,iq)=1
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

   WRITE(stdout, '(5x," I am image number ", i5 " and my work is about",i5, &
                      &  " scf runs. I calculate: ")') &
                        my_image_id, work(my_image_id)

   DO iq = 1, nqs
      IF (comp_iq(iq)==1) THEN
         WRITE(stdout, '(5x," q point number ", i5, ", representations:")') iq
         string=' '
         DO irr=0, rep_iq(iq)
            IF (comp_irr_iq(irr, iq)==1) &
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
   USE io_files,  ONLY : tmp_dir, xmlpun, prefix
   USE control_ph, ONLY : tmp_dir_ph
   USE save_ph,   ONLY : tmp_dir_save
   USE disp,      ONLY : nqs, comp_irr_iq, rep_iq
   USE xml_io_base, ONLY : copy_file
   USE mp,        ONLY : mp_barrier
   USE mp_global, ONLY : my_image_id, nimage, intra_image_comm
   USE io_global, ONLY : stdout, ionode

   IMPLICIT NONE

   INTEGER :: iq, irr
   LOGICAL :: exst
   CHARACTER(LEN=256) :: file_input, file_output
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   CALL mp_barrier(intra_image_comm)
   IF (nimage == 1) RETURN
   IF (my_image_id==0) RETURN

   DO iq=1,nqs
      DO irr=0, rep_iq(iq)
         IF (comp_irr_iq(irr,iq)==1.and.ionode) THEN
            file_input=TRIM( tmp_dir_ph ) // &
                    & TRIM( prefix ) // '.phsave' // '/' // TRIM( xmlpun ) &
                    &  // '.' // TRIM(int_to_char(iq))&
                    &  // '.' // TRIM(int_to_char(irr))

            file_output=TRIM( tmp_dir_save ) // '/' // '_ph0' // &
                    &   TRIM( prefix ) // '.phsave' // '/' // TRIM( xmlpun ) &
                    &    // '.' // TRIM(int_to_char(iq))&
                    &    // '.' // TRIM(int_to_char(irr))

            INQUIRE (FILE = TRIM(file_input), EXIST = exst)
            IF (exst) CALL copy_file(file_input, file_output)
         ENDIF
      ENDDO
   ENDDO
   RETURN
   END SUBROUTINE
