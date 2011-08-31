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
  !  for all q points. It writes the file xml.#q. It is used by unrecovered
  !  phonon runs.
  !
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp, amass
  USE cell_base,     ONLY : at, bg
  USE io_global,     ONLY : stdout
  USE symm_base,     ONLY : nrot, nsym, sr, ftau, irt, t_rev, time_reversal, &
                            sname, invs, s, inverse_s, copy_sym, &
                            s_axis_to_cart
  USE control_ph,    ONLY : rec_code, lgamma_gamma, search_sym, lgamma, &
                            where_rec, current_iq, u_from_file
  USE modes,         ONLY : u, npertx, npert, gi, gimq, nirr, &
                            t, tmq, irotmq, irgq, minus_q, &
                            nsymq, nmodes, rtau, name_rap_mode, num_rap_mode
  USE qpoint,        ONLY : xq
  USE disp,          ONLY : x_q, nqs, nsymq_iq, rep_iq, npert_iq
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE spin_orb,      ONLY : domag
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : iverbosity, modenum, noinv
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : root, world_comm

  implicit none

  real(DP) :: sr_is(3,3,48)

  integer :: ir,  isym, jsym, &
       mu, nu, irr, na, it, nt, is, js, nsym_is, iq
  ! counters

  real(DP), allocatable :: w2(:)

  logical :: sym (48), magnetic_sym
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
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! allocate and calculate rtau, the rotated position of each atom
  !
  !
  nmodes = 3 * nat
  minus_q = (modenum .eq. 0)
  ! if minus_q=.t. set_irr will search for Sq=-q+G symmetry.
  ! On output minus_q=.t. if such a symmetry has been found
  ! TEMP: set_irr_* should not find again the small group of q
  allocate( w2(3*nat) )
  DO iq=1, nqs
     xq(1:3)  = x_q(1:3,iq)
     lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
     sym(1:nsym)=.true.
     call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)
     IF ( .not. time_reversal ) minus_q = .false.
     nsymq = copy_sym ( nsym, sym )
     call inverse_s ( )
     call s_axis_to_cart ( )
     sym (1:nsym) = .true.
     call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)

     if (nsym > 1.and..not.lgamma_gamma) then
        call set_irr (nat, at, bg, xq, s, sr, tau, ntyp, ityp, ftau, invs, &
                    nsym, rtau, irt, irgq, nsymq, minus_q, irotmq, u, npert,  &
                    nirr, gi, gimq, iverbosity, u_from_file, w2, search_sym,  &
                    nspin_mag, t_rev, amass, num_rap_mode, name_rap_mode)
        npertx = 0
        DO irr = 1, nirr
           npertx = max (npertx, npert (irr) )
        ENDDO
        CALL allocate_pert()
        CALL set_irr_sym (nat, at, bg, xq, s, rtau, irt, irgq, nsymq,  &
                          minus_q, irotmq, t, tmq, u, npert, nirr, npertx )
     else
        npertx=1
        CALL allocate_pert()
        call set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, npert, &
             nirr, gi, gimq, iverbosity)
     endif
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
  u_from_file=.TRUE.

  DEALLOCATE(w2)
  DEALLOCATE (rtau)
  DEALLOCATE (u)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (npert)

  CALL stop_clock ('init_rep')
  RETURN
END SUBROUTINE init_representations
