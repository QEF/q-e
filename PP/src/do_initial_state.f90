!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE do_initial_state (excite)
  !----------------------------------------------------------------------
  !
  !    This routine is a driver routine which computes the initial state
  !    contribution to the core level shift.
  !
  !    contains five parts which are computed by different routines:
  !    a)   add_shift_lc,   contribution due to the local potential
  !    b)   add_shift_cc,   contribution due to NLCC
  !    c)   add_shift_us ,  contribution due to the non-local potential
  !    d)   add_shift_ew,   contribution due to the electrostatic ewald term
  !
  !
  USE kinds, ONLY : DP
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : at, bg, alat, omega
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE gvect,      ONLY : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm, eigts1, eigts2, eigts3
  USE fft_base,   ONLY : dfftp
  USE lsda_mod,   ONLY : nspin
  USE symme,      ONLY : symscalar
  USE vlocal,     ONLY : strf, vloc
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : lda_plus_u
  USE extfield,   ONLY : tefield, forcefield
  USE uspp,       ONLY : nkb, vkb
  USE uspp_param, ONLY : nh
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE wvfct,      ONLY : npwx
  USE ener,       ONLY : ef
  USE parameters, ONLY : ntypx
  USE control_flags, ONLY: gamma_only
  USE DFUNCT,     ONLY : newd
  USE constants,  ONLY : rytoev
  !
  IMPLICIT NONE
  !
  INTEGER :: excite(ntypx)
  INTEGER, ALLOCATABLE :: ityp_gs(:), ityp_excited(:)
  REAL(DP), ALLOCATABLE :: shift(:), &
                                shift_ef (:), &
                                shift_nl (:), &
                                shift_lc (:), &
                                shift_cc (:), &
                                shift_ion (:), &
                                shift_hub(:), &
                                delta_zv(:)
  !
  ! nonlocal, local, core-correction, ewald, and scf correction terms
  !
  INTEGER :: ipol, na, nt, ik
  ! counter on polarization
  ! counter on atoms
  LOGICAL :: first
  !
  CALL start_clock( 'do_shift' )
  !
  ALLOCATE( shift(nat), shift_ef(nat), shift_nl(nat), shift_lc(nat), &
            shift_cc(nat), shift_hub(nat), shift_ion(nat), delta_zv(ntyp) )
  ALLOCATE ( ityp_gs(nat), ityp_excited(nat) )

  ityp_gs(:) = ityp(:)
  DO nt =1,ntyp
     IF (excite(nt)<0 .or. excite(nt)>ntyp) &
        CALL errore ('do_initial_state', ' wrong excite value ', nt )
  ENDDO
  DO nt=ntyp+1, ntypx
     IF (excite(nt)/=0 ) &
        CALL errore ('do_initial_state', ' cannot exicte nt>ntyp ', nt )
  ENDDO

  ityp_gs(:) = ityp(:)
  ityp_excited(:) = ityp(:)
  DO na=1,nat
    IF (excite(ityp(na))/=0) ityp_excited (na) = excite(ityp(na))
  ENDDO

  delta_zv(:) = 0.d0
  DO nt=1,ntyp
     IF (excite(nt)/=0) delta_zv(nt) = zv(excite(nt)) - zv(nt)
  ENDDO
  !
  shift_ef(:)   = 0.D0
  shift_nl(:)   = 0.D0
  shift_lc(:)   = 0.D0
  shift_cc(:)   = 0.D0
  shift_hub(:)   = 0.D0
  shift_ion(:)   = 0.D0
  !
  WRITE( stdout, '(/,5x,"INITIAL STATE CONTRIBUTION TO", &
                 & /,5x,"CORE LEVEL SHIFT ON ATOMS:", / )')

  DO na=1,nat
     shift_ef(na) = ef * delta_zv(ityp(na))
  ENDDO

    first = .true.
 10 CONTINUE
  !
  ! ... The  nonlocal contribution is computed here
  !
  CALL add_shift_us( shift_nl )
  !
  ! ... The local contribution
  !
  CALL add_shift_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 dfftp%nnr, g, rho%of_r, nl, nspin, &
                 gstart, gamma_only, vloc, shift_lc )
  !
  ! ... The NLCC contribution
  !
  CALL add_shift_cc( shift_cc )
  !
  ! ... The Hubbard contribution
  !
  IF ( lda_plus_u ) CALL errore('do_initial_state','LDA+U not implemented',1)

  !
  ! change atomic type and recompute needed quantities
  !
  IF ( first ) THEN
     ityp(:) = ityp_excited(:)
     CALL newd()
     nkb = 0
     DO na = 1, nat
        nkb = nkb + nh (ityp(na))
     ENDDO
     DEALLOCATE(vkb)
     IF(nkb>0) ALLOCATE(vkb(npwx,nkb))

     IF ( nks == 1 ) THEN
        ik = 1
        IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     ENDIF
     shift_nl = - shift_nl
     shift_lc = - shift_lc
     shift_cc = - shift_cc
     shift_hub= - shift_hub
     first = .false.
     GOTO 10
  ELSE
     ityp(:) = ityp_gs(:)
     CALL newd()
     nkb = 0
     DO na = 1, nat
        nkb = nkb + nh (ityp(na))
     ENDDO
     DEALLOCATE(vkb)
     IF(nkb>0) ALLOCATE(vkb(npwx,nkb))
     IF ( nks == 1 ) THEN
        ik = 1
        IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     ENDIF
  ENDIF

  !
  ! ... The ionic contribution is computed here
  !
!  call infomsg ('do_initial_state',' EWALD term is still missing')
  CALL do_shift_ew (alat, nat, ntyp, ityp, zv, delta_zv, at, bg, tau, &
     omega, g, gg, ngm, gcutm, gstart, gamma_only, shift_ion)
  !
  ! ... here we sum all the contributions and compute the total force acting
  ! ... on the crstal
  !
  DO na = 1, nat
     shift(na) = shift_ef(na)  + &
                 shift_nl(na)  + &
                 shift_ion(na) + &
                 shift_lc(na)  + &
                 shift_cc(na)  + &
                 shift_hub(na)
  ENDDO
  !
  ! ... resymmetrize (should not be needed, but ...)
  !
  CALL symscalar( nat, shift )
  !
  ! ... write on output the initial state core level shifts
  !
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na),  shift(na), shift(na)*rytoev
  ENDDO
  WRITE (stdout,*)
#define DEBUG
#ifdef DEBUG
  WRITE( stdout, '(5x,"The FERMI ENERGY contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_ef(na), shift_ef(na)*rytoev
  ENDDO
  WRITE( stdout, '(5x,"The NON LOCAL contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_nl(na), shift_nl(na)*rytoev
  ENDDO
  WRITE( stdout, '(5x,"The LOCAL contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_lc(na), shift_lc(na)*rytoev
  ENDDO
  WRITE( stdout, '(5x,"The IONIC contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_ion(na), shift_ion(na)*rytoev
  ENDDO
  WRITE( stdout, '(5x,"The CC contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_cc(na), shift_cc(na)*rytoev
  ENDDO
  WRITE( stdout, '(5x,"The Hubbard contribution to shift")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), shift_hub(na), shift_hub(na)*rytoev
  ENDDO
#endif
  !
  DEALLOCATE( shift_ef, shift_nl, shift_lc, shift_cc, shift_hub, &
              shift_ion, delta_zv )
  DEALLOCATE (ityp_gs, ityp_excited)
  !
  CALL stop_clock( 'do_shift' )
  !
  RETURN
  !
9035 FORMAT(5X,'atom ',I3,' type ',I2,'   shift =',F13.6,' Ry, =',F13.5,' eV')
!
END SUBROUTINE do_initial_state

