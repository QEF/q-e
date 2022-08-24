!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_if_partial_dyn(u, nirr, npert, comp_irr)
!
!! This routine decides which irreducible representation have to
!! be computed for each q point on the basis of start_irr and last_irr
!! or on the basis of nat_todo. It sets the array comp_irr.  
!! If this is part of a dispersion calculation the routine has to be
!! called separately for each q and in that case different displacement
!! patterns are given as input. Note that this routine is called before
!! distributing the irrep among the images and only the irrep selected
!! here are distributed.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE symm_base, ONLY : irt
USE partial, ONLY : nat_todo, atomo
USE control_ph, ONLY : start_irr, last_irr, ldiag
USE control_flags, ONLY : modenum

USE lr_symm_base, ONLY : nsymq

IMPLICIT NONE

COMPLEX(DP), INTENT(IN) :: u( 3*nat, 3*nat )
INTEGER, INTENT(IN) :: nirr, npert( 3*nat )
LOGICAL, INTENT(OUT) :: comp_irr(0:3*nat)

INTEGER, ALLOCATABLE :: ifat(:)

INTEGER :: na, isym, mu, nu, ipol, imode0, irr, ipert
INTEGER :: last_irr_eff
!
!  If modenum is specified, only this mode is calculated and nat_todo is
!  ignored
!
  comp_irr = .FALSE.
  comp_irr(0) = .TRUE.
  IF (modenum /= 0) THEN
     comp_irr(modenum)=.TRUE.
     RETURN
  ENDIF

  ALLOCATE ( ifat(nat) )
  IF (nat_todo > 0) THEN
     !
     !   Sets the atoms which must be computed: the requested
     !   atoms and all the symmetry related atoms
     !
     ifat = 0
     DO na = 1, nat_todo
        IF ( atomo(na)>nat .OR. atomo(na)<1 ) &
           CALL errore('phq_setup', &
                       'one of atoms to do (nat_todo) is < 0 or > nat', 1)
        ifat (atomo (na) ) = 1
        DO isym = 1, nsymq
           ifat (irt (isym, atomo (na) ) ) = 1
        ENDDO
     ENDDO
     !
     !  Find the irreducible representations where the required atoms moves
     !
     imode0 = 0
     do irr = 1, nirr
        do ipert = 1, npert (irr)
           mu = imode0 + ipert
           do na = 1, nat
              if (ifat (na) == 1 .and. .NOT.comp_irr (irr) ) then
                 do ipol = 1, 3
                    nu = 3 * (na - 1) + ipol
                    if (abs (u (nu, mu) ) > 1.d-6)  comp_irr (irr) = .TRUE.
                 enddo
              endif
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  ELSE
     !
     !  nat_todo is not given. In principle all representation have
     !  to be calculated if not found somewhere else or limited below
     !
     comp_irr=.TRUE.
  ENDIF
!
! The representations that are smaller than start_irr or larger than
! last_irr_eff must be removed
!
  last_irr_eff=last_irr
  IF (last_irr > nirr.or.last_irr<0) last_irr_eff=nirr
  IF (start_irr > 1) comp_irr(0:start_irr-1) = .FALSE.
  IF (last_irr_eff < nirr ) comp_irr(last_irr_eff+1:nirr) = .FALSE.

!
! When ldiag=.true. the partial dynamical matrix is diagonalized.
! It must contain also the part calculated by dynmat0
!  
  IF (ldiag) comp_irr(0)=.TRUE.

  DEALLOCATE(ifat)

  RETURN
END SUBROUTINE check_if_partial_dyn

SUBROUTINE set_ifat(nat,   nat_todo, atomo, nsym, irt, ifat)
  !! sets to 1 the value in ifat for each atom in atomo list and for their 
  !! symmetry equivalents. 
  IMPLICIT NONE
  INTEGER,INTENT(IN)   :: nat 
  !! total of number of atoms, size of ifat 
  INTEGER,INTENT(IN)   :: nat_todo
  !! number of displacicng atoms, size of atomo 
  INTEGER,INTENT(IN)   :: atomo(nat_todo) 
  !! contains the atomic index of the displacing atoms 
  INTEGER,INTENT(IN)   :: nsym 
  !! number of symmetries 
  INTEGER,INTENT(IN)   :: irt(48,nat)
  !! list of equivalent atoms per each operation symmetries,
  !! the routine uses only the irt(1:nsym,1:nat) slice   
  INTEGER,INTENT(OUT)  :: ifat(nat) 
  !! on output contains 0 for fixed atoms and 1 for displacing ones 
  ! 
  INTEGER :: na, isym  
  !  
  IF (nat_todo == 0) THEN 
   ifat = 1 
   RETURN  
  END IF 

  IF (MAXVAL(atomo) > nat .OR. MINVAL(atomo) < 1) & 
   CALL errore("set_ifat:", "internal error: atomo list is inconsistent",1)
  !
  ifat = 0  
  DO na = 1, nat_todo
    DO isym = 1, nsym 
      ifat( irt( isym, atomo(na)) ) = 1
    END DO 
  END DO 
END SUBROUTINE set_ifat  
   
SUBROUTINE set_local_atomo(nat, nat_todo, atomo, nsym, irt,  nat_l, atomo_l) 
   IMPLICIT NONE 
   INTEGER,INTENT(IN)               :: nat, nat_todo, nsym, atomo(nat_todo), irt(48,nat)
   !! :nat: total number of atoms
   !! :nat_todo: number of atoms effectively displaced 
   !! :nsym: number of symmetries in the system
   !! :atomo: list of atoms to be displaced before symmetrization 
   !! :irt: atoms corresponding atom for each sym operation and atom
   INTEGER,INTENT(OUT)              :: nat_l 
   !! actual number of atoms to be displaced considering symmetries
   INTEGER,ALLOCATABLE,INTENT(OUT)  :: atomo_l(:)
   !! list with the indeces of all the atoms to be displaced 
   !
   INTEGER,ALLOCATABLE  :: ifat(:)
   INTEGER              :: na
   ALLOCATE(ifat(nat))
   CALL set_ifat(nat, nat_todo, atomo, nsym, irt, ifat)
   nat_l = COUNT(ifat == 1)
   ALLOCATE (atomo_l(nat_l))
   atomo_l = PACK([(na,na=1,nat)], MASK = ifat == 1)
   DEALLOCATE(ifat)
END SUBROUTINE set_local_atomo 