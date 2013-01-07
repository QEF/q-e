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
  !  for all q points. It writes the file data-file.#q.x. It is used by 
  !  unrecovered phonon runs.
  !
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp, amass
  USE cell_base,     ONLY : at, bg
  USE symm_base,     ONLY : nrot, nsym, sr, ftau, irt, t_rev, time_reversal, &
                            sname, invs, s  
  USE control_ph,    ONLY : rec_code, search_sym, search_sym_save, lgamma, &
                            where_rec, current_iq, u_from_file
  USE modes,         ONLY : u, npertx, npert, gi, gimq, nirr, &
                            t, tmq, irotmq, minus_q, invsymq, &
                            nsymq, nmodes, rtau, name_rap_mode, num_rap_mode
  USE qpoint,        ONLY : xq
  USE disp,          ONLY : x_q, nqs, nsymq_iq, rep_iq, npert_iq
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : domag
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : modenum, noinv
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : root, world_comm

  implicit none

  real(DP) :: sr_is(3,3,48)

  integer :: ir,  isym, jsym, &
       mu, nu, irr, na, it, nt, is, js, nsym_is, iq
  ! counters

  logical :: sym (48), magnetic_sym, is_symmorphic
  ! the symmetry operations
  integer :: ierr

  call start_clock ('init_rep')

  allocate (rtau ( 3, 48, nat))
  allocate (u ( 3 * nat, 3 * nat))
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat))
  allocate (npert ( 3 * nat))

  name_rap_mode=' '
  u_from_file=.FALSE.

  magnetic_sym = noncolin .AND. domag
  !
  ! allocate and calculate rtau, the rotated position of each atom
  !
  nmodes = 3 * nat
  minus_q = (modenum .eq. 0)
  ! if minus_q=.t. set_irr will search for Sq=-q+G symmetry.
  ! On output minus_q=.t. if such a symmetry has been found
  ! TEMP: set_irr_* should not find again the small group of q
  !
  DO iq=1, nqs
     xq(1:3)  = x_q(1:3,iq)
     lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
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
!    if minus_q is true calculate also irotmq and the G associated to Sq=-q+G
!
     CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)


     is_symmorphic=.NOT.(ANY(ftau(:,1:nsymq) /= 0))
     search_sym=search_sym_save
     IF (.NOT.is_symmorphic) THEN
        DO isym=1,nsymq
           search_sym=( search_sym.and.(abs(gi(1,isym))<1.d-8).and.  &
                                       (abs(gi(2,isym))<1.d-8).and.  &
                                       (abs(gi(3,isym))<1.d-8) )
        END DO
     END IF
     num_rap_mode=-1
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

     nsymq_iq(iq) = nsymq
     rep_iq(iq) = nirr
     DO irr=1, nirr
        npert_iq(irr,iq)=npert(irr)
     ENDDO

     current_iq=iq
     where_rec='init_rep..'
     rec_code=-50
     CALL ph_writefile('data',0)
     CALL deallocate_pert()
  ENDDO
  search_sym=search_sym_save
  u_from_file=.TRUE.

  DEALLOCATE (rtau)
  DEALLOCATE (u)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (npert)

  CALL stop_clock ('init_rep')
  RETURN
END SUBROUTINE init_representations
