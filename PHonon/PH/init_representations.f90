!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine init_representations()
  !-----------------------------------------------------------------------
  !
  !  This subroutine initializes the modes of all irreducible representations
  !  for all q points. It writes the files patterns.#q.xml in the outdir 
  !  directory. It is used by unrecovered  phonon runs. The small group of 
  !  q must be calculated for each q. Note that all images receives the 
  !  same modes calculated by the root processor and save them on file. 
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat
  USE cell_base,     ONLY : at, bg
  USE io_global,     ONLY : stdout
  USE symm_base,     ONLY : nsym, sr, ftau, irt, time_reversal, t_rev, s
  USE control_ph,    ONLY : search_sym, current_iq, u_from_file, &
                            search_sym_save
  USE modes,         ONLY : u, npert, nirr, nmodes, name_rap_mode, &
                            num_rap_mode
  USE disp,          ONLY : x_q, nqs, lgamma_iq
  USE cryst_ph,      ONLY : magnetic_sym
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : modenum, noinv
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : root, world_comm

  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_lr,    ONLY : lgamma

  implicit none

  integer ::  isym, irr, iq
  ! counters
  LOGICAL, EXTERNAL :: symmorphic_or_nzb
  integer :: ierr

  call start_clock ('init_rep')

  allocate (rtau ( 3, 48, nat))
  allocate (u ( 3 * nat, 3 * nat))
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat))
  allocate (npert ( 3 * nat))

  u_from_file=.FALSE.
  !
  ! allocate and calculate rtau, the rotated position of each atom
  !
  nmodes = 3 * nat
  minus_q = (modenum .eq. 0)
  IF ( .not. time_reversal ) minus_q = .false.
  ! if minus_q=.t. set_irr will search for Sq=-q+G symmetry.
  ! On output minus_q=.t. if such a symmetry has been found
  DO iq=1, nqs
     xq(1:3)  = x_q(1:3,iq)
     lgamma = lgamma_iq(iq)
!
!    search for the small group of q
!
     CALL set_small_group_of_q(nsymq,invsymq,minus_q)
!
!    calculate rtau with the new symmetry order
!
     CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
!
!    and calculate the vectors G associated to the symmetry Sq = q + G
!    if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
!
     CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!
!    Check if we can search symmetry for this q point
!
     search_sym = search_sym_save .AND. symmorphic_or_nzb()
     num_rap_mode=-1
     name_rap_mode=' '
     IF (search_sym) CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)

     CALL find_irrep()
!
!  Only the modes calculated by node zero are sent to all images
!
     CALL mp_bcast (u, root, world_comm)
     CALL mp_bcast (nsymq, root, world_comm)
     CALL mp_bcast (npert, root, world_comm)
     CALL mp_bcast (nirr, root, world_comm)
     CALL mp_bcast (name_rap_mode, root, world_comm)
     CALL mp_bcast (num_rap_mode, root, world_comm)

     CALL ph_writefile('data_u',iq,0,ierr)
  ENDDO
  u_from_file=.TRUE.
  search_sym=search_sym_save

  DEALLOCATE (rtau)
  DEALLOCATE (u)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (npert)

  CALL stop_clock ('init_rep')
  RETURN
END SUBROUTINE init_representations

!-----------------------------------------------------------------------
subroutine initialize_grid_variables()
  !-----------------------------------------------------------------------
  !
  !  This subroutine initializes the grid variables by reading the
  !  modes from file. It uses the routine check_if_partial_dyn to 
  !  set the modes to compute according to start_irr, last_irr or
  !  modenum and ifat flags.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat
  USE modes,         ONLY : u, npert, nirr, name_rap_mode, num_rap_mode
  USE disp,          ONLY : nqs, comp_iq
  USE partial,       ONLY : comp_irr
  USE grid_irr_iq,   ONLY : nsymq_iq, irr_iq, npert_irr_iq, comp_irr_iq
  USE ph_restart,    ONLY : ph_readfile
  USE control_ph,    ONLY : start_q, last_q, epsil, zeu
  USE el_phon,       ONLY : elph
  USE io_global,     ONLY : stdout
  USE mp_global,     ONLY : mp_global_end
  USE environment,   ONLY : environment_end

  USE lr_symm_base, ONLY : nsymq

  implicit none

  INTEGER ::  irr, iq
  ! counters
  INTEGER :: ierr
  LOGICAL :: something_to_do

  allocate (u ( 3 * nat, 3 * nat))
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat))
  allocate (npert ( 3 * nat))

  DO iq=1, nqs
!
!  Read from file the modes and the representations
!
     CALL ph_readfile('data_u', iq, 0, ierr)
     IF (ierr /= 0) call errore('initialize_grid_variables',&
                                'problems reading u',1)

     nsymq_iq(iq) = nsymq
     irr_iq(iq) = nirr
     DO irr=1, nirr
        npert_irr_iq(irr,iq)=npert(irr)
     ENDDO
!
!    here we deal with start_irr, last_irr, OR of modenum OR of ifat atomo
!
     CALL check_if_partial_dyn(u, nirr, npert, comp_irr)
     comp_irr_iq(:,iq)=comp_irr(:)
  ENDDO
!
!  here deal with the start_q, last_q flags
!
  comp_iq=.FALSE.
  something_to_do=.FALSE.
  DO iq=1,nqs
     IF (iq>=start_q.AND.iq<=last_q) THEN
        DO irr=0,irr_iq(iq)
           IF (comp_irr_iq(irr,iq)) THEN
              comp_iq(iq)=.TRUE.
              something_to_do=.TRUE.
           ENDIF
        ENDDO
     ELSE
        comp_irr_iq(:,iq)=.FALSE.
     ENDIF
  ENDDO

  DEALLOCATE (u)
  DEALLOCATE (npert)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  IF (.NOT.(something_to_do.OR.epsil.OR.zeu.OR.elph)) THEN
     write(stdout,'(/,5x, "The code stops because there is nothing to do")') 
     CALL clean_pw(.FALSE.)
     CALL close_files(.FALSE.)
     CALL environment_end('PHONON')
     CALL mp_global_end()
     STOP
  ENDIF

  RETURN
END SUBROUTINE initialize_grid_variables

