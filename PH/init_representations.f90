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
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp, pmass
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout
  USE symm_base,     ONLY : nrot, nsym, s, ftau, irt, t_rev, time_reversal, &
                            sname, invs, inverse_s, copy_sym
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr,&
                                  char_mat, name_rap, gname, name_class, ir_ram
  USE rap_point_group_is,   ONLY : code_group_is, gname_is
  USE control_ph,    ONLY : rec_code, lgamma_gamma, search_sym, lgamma, &
                            where_rec, current_iq
  USE modes,         ONLY : u, npertx, npert, gi, gimq, nirr, &
                            t, tmq, irotmq, irgq, minus_q, &
                            nsymq, nmodes, rtau, name_rap_mode
  USE qpoint,        ONLY : xq
  USE disp,          ONLY : x_q, nqs, nsymq_iq, rep_iq, npert_iq
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : domag
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : iverbosity, modenum, noinv
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : root, world_comm

  implicit none

  real(DP) :: sr(3,3,48), sr_is(3,3,48)

  integer :: ir,  isym, jsym, &
       mu, nu, irr, na, it, nt, is, js, nsym_is, iq
  ! counters
  
real(DP), allocatable :: w2(:)

  logical :: sym (48), is_symmorphic, magnetic_sym, u_from_file
  ! the symmetry operations
  integer :: ierr

  call start_clock ('init_rep')

  allocate (rtau ( 3, 48, nat))    
  allocate (u ( 3 * nat, 3 * nat))    
  allocate (name_rap_mode( 3 * nat))    
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
     sym (1:nsym) = .true.
     call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)

     if (nsym > 1.and..not.lgamma_gamma) then
        call set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, u, npert, &
             nirr, gi, gimq, iverbosity,u_from_file,w2)
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
     is_symmorphic=.true.
     DO isym=1,nsymq
        is_symmorphic=( is_symmorphic.and.(ftau(1,irgq(isym))==0).and.  &
                                          (ftau(2,irgq(isym))==0).and.  &
                                          (ftau(3,irgq(isym))==0) )
  
     END DO
     search_sym=.true.
     IF (.not.is_symmorphic) THEN
        DO isym=1,nsymq
           search_sym=( search_sym.and.(abs(gi(1,irgq(isym)))<1.d-8).and.  &
                                       (abs(gi(2,irgq(isym)))<1.d-8).and.  &
                                       (abs(gi(3,irgq(isym)))<1.d-8) )
        END DO
     END IF
     IF (search_sym) THEN
        DO isym=1,nsymq
           CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
        END DO
        CALL find_group(nsymq,sr,gname,code_group)
        CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
        CALL divide_class(code_group,nsymq,sr,nclass,nelem,elem,which_irr)
        IF (noncolin .and. domag) THEN
           nsym_is=0.d0
           DO isym=1,nsymq
              IF (t_rev(isym)==0) THEN
                 nsym_is=nsym_is+1
                 CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
              ENDIF
           END DO
           CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
        ENDIF
        IF (.not.lgamma_gamma.and.modenum==0) &
              CALL find_mode_sym (u, w2, at, bg, nat, nsymq, &
                        s, irt, xq, rtau, pmass, ntyp, ityp, 0)
     ENDIF
!
!  Only the modes calculated by node zero are sent to all images
!
     CALL mp_bcast (u, root, world_comm)
     CALL mp_bcast (nsymq, root, world_comm)
     CALL mp_bcast (npert, root, world_comm)
     CALL mp_bcast (nirr, root, world_comm)
     CALL mp_bcast (name_rap_mode, root, world_comm)
  
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

  DEALLOCATE(w2)
  DEALLOCATE (rtau)    
  DEALLOCATE (u)    
  DEALLOCATE (name_rap_mode)    
  DEALLOCATE (npert)    

  CALL stop_clock ('init_rep')
  RETURN
END SUBROUTINE init_representations
