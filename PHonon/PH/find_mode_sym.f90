!
! Copyright (C) 2006-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym_new (u, w2, tau, nat, nsym, s, sr, irt, xq,    &
     rtau, amass, ntyp, ityp, flag, lmolecule, lstop, num_rap_mode, ierr)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of eigenvectors of the dynamical
  !   matrix. It does NOT work at zone border in non symmorphic space groups.
  !   if flag=1 the true displacements are given in input, otherwise the
  !   eigenvalues of the dynamical matrix are given.
  !   The output of this routine is only num_rap_mode, the number of
  !   the irreducible representation for each mode.
  !   error conditions:
  !   num_rap_mode(i)= 0   ! the routine could not determine mode symmetry
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, RY_TO_CMM1
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_is, ONLY : gname_is
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::             &
       nat,         &     ! number of atoms
       nsym,        &     ! number of symmetries
       flag,        &     ! if 1 u are displacements, if 0 u are eigenvectors
       ntyp,        &     ! number of atomic types
       ityp(nat),   &     ! the type of each atom
       irt(48,nat)        ! the rotated of each atom
  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )

  INTEGER, INTENT(OUT) :: ierr ! 0 if the routine determined mode symmetry

  REAL(DP), INTENT(IN) ::   &
       xq(3),          &  ! the q vector of the modes
       tau(3,nat),     &  ! the atomic coordinates
       rtau(3,48,nat), &  ! the R vector for each rotated atom
       amass(ntyp),    &  ! the mass of the atoms
       w2(3*nat),      &  ! the square of the frequencies
       sr(3,3,48)         ! the rotation matrices in real space.

  COMPLEX(DP), INTENT(IN) ::  &
       u(3*nat, 3*nat)       ! The eigenvectors or the displacement pattern

  LOGICAL, INTENT(IN) :: lmolecule, & ! if .true. these are eigenvalues of an
                                   ! isolated system and do not find the
                                   ! symmetry of the first six eigenvectors,
                                   ! or five for a linear molecule.
                         lstop     ! if .true. the routine stops if it
                                   ! does not understand the symmetry of a 
                                   ! mode

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ngroup,    &   ! number of different frequencies groups
       s(3,3,48), &   ! rotation matrices
       nmodes,    &   ! number of modes
       imode,     &   ! counter on modes
       igroup,    &   ! counter on groups
       nu_i, mu,  &   ! counters on modes
       irot,      &   ! select a rotation
       irap,      &   ! counter on representations
       iclass,    &   ! counter on classes
       na,        &   ! counter on atoms
       i              ! generic counter

  INTEGER, ALLOCATABLE :: istart(:), dim_rap(:)

  COMPLEX(DP) :: times              ! safe dimension
  ! in case of accidental degeneracy
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), ALLOCATABLE :: w1(:)
  COMPLEX(DP), ALLOCATABLE ::  rmode(:,:), trace(:,:), z(:,:)
  LOGICAL :: is_linear
  INTEGER :: counter, counter_s
  LOGICAL :: found
  INTEGER :: invs(48), ss(3,3), isym, jsym
  !
  !    Divide the modes on the basis of the mode degeneracy.
  !
  ierr=0
  num_rap_mode=0
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(dim_rap(nmodes))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(w1(nmodes))
  ALLOCATE(rmode(nmodes,nmodes))
  ALLOCATE(trace(48,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvalues of the dynmaical matrix
     !  Note that amass is in amu; amu_ry converts it to Ry au
     !
     DO nu_i = 1, nmodes
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           z (mu, nu_i) = u (mu, nu_i) * SQRT (amu_ry*amass (ityp (na) ) )
        END DO
     END DO
  ELSE
     z=u
  ENDIF

  DO isym = 1, nsym
     found = .false.
     DO jsym = 1, nsym
        !
        ss = matmul (s(:,:,jsym),s(:,:,isym))
        ! s(:,:,1) is the identity
        IF ( all ( s(:,:,1) == ss(:,:) ) ) THEN
           invs (isym) = jsym
           found = .true.
        ENDIF
     ENDDO
     IF ( .NOT.found) CALL errore ('inverse_s', ' Not a group', 1)
  ENDDO

!
!  Compute the mode frequency in cm-1. Two modes are considered degenerate
!  if their frequency is lower 0.05 cm-1
! 
  w1(:)=SIGN(SQRT(ABS(w2(:)))*RY_TO_CMM1,w2(:))

  ngroup=1
  istart(ngroup)=1
!
!  The symmetry of these modes is not computed
!
  IF (lmolecule) THEN
     istart(1)=7
     IF(is_linear(nat,tau)) istart(1)=6
  ENDIF
!
! The other modes are divided into groups of degenerate modes
!
  DO imode=istart(1)+1,nmodes
     IF (ABS(w1(imode)-w1(imode-1)) > 5.0d-2) THEN
        ngroup=ngroup+1
        istart(ngroup)=imode
     END IF
  END DO
  istart(ngroup+1)=nmodes+1
  !
  !  Find the character of one symmetry operation per class
  !
  DO igroup=1,ngroup
     dim_rap(igroup)=istart(igroup+1)-istart(igroup)
     DO iclass=1,nclass
        irot=elem(1,iclass)
!
!   rotate all modes together
!
        CALL rotate_mod(z,rmode,sr(1,1,irot),irt,rtau,xq,nat,invs(irot))
        trace(iclass,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap(igroup)
           nu_i=istart(igroup)+i-1
           trace(iclass,igroup)=trace(iclass,igroup) + &
                zdotc(3*nat,z(1,nu_i),1,rmode(1,nu_i),1)
        END DO
!              write(6,*) 'group,class',igroup, iclass, trace(iclass,igroup)
     END DO
  END DO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  DO igroup=1,ngroup
     counter=istart(igroup)
!
!   If the frequency is so small probably it has not been calculated.
!   This value 
!
     IF (ABS(w1(counter))<1.d-3) CYCLE
     DO irap=1,nclass
        times=(0.d0,0.d0)
        DO iclass=1,nclass
           times=times+trace(iclass,igroup)*CONJG(char_mat(irap, &
                which_irr(iclass)))*nelem(iclass)
           !         write(6,*) igroup, irap, iclass, which_irr(iclass)
        ENDDO
        times=times/nsym
!
!   times must be a positive integer or zero, otherwise some error occured
!   somewhere
!
        IF ((ABS(NINT(ABS(DBLE(times)))-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN 
           IF (lstop) THEN
              CALL errore('find_mode_sym','unknown mode symmetry',1)
           ELSE
              counter=counter + dim_rap(igroup)-1
              ierr=1
           ENDIF
        ELSE
!
!    If the code arrives here, no error occured and we can set the mode
!    symmetry for all the modes of the group
!
           IF (ABS(times) > eps) THEN
              IF (ABS(NINT(DBLE(times))-DBLE(times)) < 1.d-4) THEN
                 counter_s=counter
                 DO imode=counter_s, counter_s+NINT(DBLE(times))*&
                                              NINT(DBLE(char_mat(irap,1)))-1
                    num_rap_mode(imode) = irap
                    counter=counter+1
                 ENDDO
              END IF
           END IF
        END IF
     END DO
  END DO

100 CONTINUE

  DEALLOCATE(trace)
  DEALLOCATE(z)
  DEALLOCATE(w1)
  DEALLOCATE(rmode)
  DEALLOCATE(dim_rap)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym_new

SUBROUTINE rotate_mod(mode,rmode,sr,irt,rtau,xq,nat,irot)

  USE kinds, ONLY : DP
  USE constants, ONLY: tpi
  USE cell_base, ONLY : bg
  !
  !  irot must be the number of the rotation S^-1
  !
  IMPLICIT NONE

  INTEGER :: nat,        & !  number of atoms
             irot,       & !  the index of the inverse of the rotation sr
             irt(48,nat)   !  the rotated of each atom for all rotations

  COMPLEX(DP) :: mode(3*nat,3*nat), & ! the mode to rotate
                 rmode(3*nat,3*nat)  ! the rotated mode

  REAL(DP)  :: sr(3,3),  &  ! the symmetry matrices in cartesian coordinates
               rtau(3,48,nat), & ! the vector R = S tau - tau for all rotations
               xq(3)       ! the q vector

  COMPLEX(DP) :: phase   ! auxiliary phase
  REAL(DP)    :: arg     ! an auxiliary argument

  INTEGER :: na, nb,     & ! counters on atoms
             ipol, jpol, & ! counters on coordinates
             mu_i, mu_j    ! counters on modes coordinates

  rmode=(0.d0,0.d0)
  DO na=1,nat
     nb=irt(irot,na)
     arg = ( xq(1)*rtau(1,irot,na) + xq(2)*rtau(2,irot,na) +  &
             xq(3)*rtau(3,irot,na) ) * tpi
     phase = CMPLX(cos(arg), sin(arg), kind=DP)
     DO ipol=1,3
        mu_i=3*(na-1)+ipol
        DO jpol=1,3
           mu_j=3*(nb-1)+jpol
           rmode(mu_i,:)=rmode(mu_i,:) + sr(ipol,jpol)*mode(mu_j,:)*phase
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE rotate_mod

FUNCTION is_linear(nat,tau)
  !
  !  This function is true if the nat atoms are all on the same line
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  LOGICAL :: is_linear
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: tau(3,nat)
  REAL(DP) :: u(3), v(3), umod, vmod
  INTEGER :: na

  is_linear=.TRUE.
  IF (nat<=2) RETURN

  u(:)=tau(:,2)-tau(:,1)
  umod=sqrt(u(1)**2+u(2)**2+u(3)**2)
  DO na=3,nat
     v(:)=tau(:,na)-tau(:,1)
     vmod=sqrt(v(1)**2+v(2)**2+v(3)**2)
     is_linear=is_linear.AND.(abs(1.0_DP- &
          abs(u(1)*v(1)+u(2)*v(2)+u(3)*v(3))/umod/vmod)<1.d-4)
  ENDDO

  RETURN
END FUNCTION is_linear

SUBROUTINE print_mode_sym(w2, num_rap_mode, lir)
!
!  This routine prints the eigenvalues of the dynamical matrix and the 
!  symmetry of their eigenvectors. If lir is true it writes also 
!  which modes are infrared and/or raman active.
!
USE kinds, ONLY : DP
USE constants, ONLY : ry_to_cmm1
USE noncollin_module, ONLY : nspin_mag
USE ions_base, ONLY : nat
USE io_global, ONLY : stdout
USE rap_point_group, ONLY : char_mat, name_rap, gname, ir_ram
USE rap_point_group_is, ONLY : gname_is

IMPLICIT NONE
REAL(DP), INTENT(IN) :: w2( 3*nat )
INTEGER, INTENT(IN) :: num_rap_mode( 3*nat )
LOGICAL, INTENT(IN) :: lir

REAL(DP) :: w1( 3*nat )
INTEGER :: next, irap, imode
CHARACTER(LEN=3) :: cdum
!
!  Transform the frequencies to cm^-1
!
w1(:)=SIGN(SQRT(ABS(w2(:)))*ry_to_cmm1,w2(:))
!
!  prints the name of the point group 
!
IF ( nspin_mag == 4 ) THEN
   WRITE(stdout,  &
          '(/,5x,"Mode symmetry, ",a11," [",a11,"] magnetic point group:",/)') &
          gname, gname_is
ELSE
      WRITE(stdout,'(/,5x,"Mode symmetry, ",a11," point group:",/)') gname
END IF
!
! for each mode, or group of degenerate modes, writes the name of the
! irreducible representation
!
next=0
DO imode = 1, 3 * nat
   IF ( imode < next .OR. ABS(w1(imode)) < 1.d-3 ) CYCLE
   IF (num_rap_mode(imode) == 0)  THEN
      WRITE(stdout,'(5x,"freq (",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') imode, imode, w1(imode)
   ELSE
      irap=num_rap_mode(imode)
      next = imode + NINT(DBLE(char_mat(irap,1)))
      cdum="   "
      IF (lir) cdum=TRIM(ir_ram(irap))
      WRITE(stdout,'(5x,"freq (",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",a19)') &
           imode, next-1, w1(imode), name_rap(irap)//" "//cdum
   ENDIF
ENDDO

RETURN
END SUBROUTINE print_mode_sym
