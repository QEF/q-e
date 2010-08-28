!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym (u, w2, at, bg, tau, nat, nsym, sr, irt, xq, &
     rtau, amass, ntyp, ityp, flag, lri, lmolecule, nspin_mag,        &
     name_rap_mode, num_rap_mode)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of eigenvectors of the dynamical
  !   matrix. It does NOT work at zone border in non symmorphic space groups.
  !   if flag=1 the true displacements are given in input, otherwise the
  !   eigenvalues of the dynamical matrix are given.
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : RY_TO_CMM1
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_is, ONLY : gname_is
  IMPLICIT NONE

  CHARACTER(15), INTENT(OUT) :: name_rap_mode( 3 * nat )
  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )
  INTEGER, INTENT(IN) :: nspin_mag

  INTEGER, INTENT(IN) ::             &
       nat,         &
       nsym,        &
       flag,        &
       ntyp,        &
       ityp(nat),   &
       irt(48,nat)

  REAL(DP), INTENT(IN) ::   &
       at(3,3),        &
       bg(3,3),        &
       xq(3),          &
       tau(3,nat),     &
       rtau(3,48,nat), &
       amass(ntyp),    &
       w2(3*nat),      &
       sr(3,3,48)

  COMPLEX(DP), INTENT(IN) ::  &
       u(3*nat, 3*nat)       ! The eigenvectors or the displacement pattern
  LOGICAL, INTENT(IN) :: lri      ! if .true. print the Infrared/Raman flag
  LOGICAL, INTENT(IN) :: lmolecule ! if .true. these are eigenvalues of an
                                   ! isolated system

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ngroup, &   ! number of different frequencies groups
       nmodes, &   ! number of modes
       imode, imode1, igroup, dim_rap, nu_i,  &
       irot, irap, iclass, mu, na, i

  INTEGER, ALLOCATABLE :: istart(:)

  COMPLEX(DP) :: times              ! safe dimension
  ! in case of accidental degeneracy
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), ALLOCATABLE :: w1(:)
  COMPLEX(DP), ALLOCATABLE ::  rmode(:), trace(:,:), z(:,:)
  LOGICAL :: is_linear
  CHARACTER(3) :: cdum
  INTEGER :: counter, counter_s
  !
  !    Divide the modes on the basis of the mode degeneracy.
  !
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(w1(nmodes))
  ALLOCATE(rmode(nmodes))
  ALLOCATE(trace(48,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvalues of the dynmaical matrix
     !
     DO nu_i = 1, nmodes
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           z (mu, nu_i) = u (mu, nu_i) * SQRT (amass (ityp (na) ) )
        END DO
     END DO
  ELSE
     z=u
  ENDIF

  DO imode=1,nmodes
     w1(imode)=SIGN(SQRT(ABS(w2(imode)))*RY_TO_CMM1,w2(imode))
  ENDDO

  ngroup=1
  istart(ngroup)=1
  imode1=1
  IF (lmolecule) THEN
     istart(ngroup)=7
     imode1=6
     IF(is_linear(nat,tau)) istart(ngroup)=6
  ENDIF
  DO imode=imode1+1,nmodes
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
     dim_rap=istart(igroup+1)-istart(igroup)
     DO iclass=1,nclass
        irot=elem(1,iclass)
        trace(iclass,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap
           nu_i=istart(igroup)+i-1
           CALL rotate_mod(z(1,nu_i),rmode,sr(1,1,irot),irt,rtau,xq,nat,irot)
           trace(iclass,igroup)=trace(iclass,igroup) + &
                zdotc(3*nat,z(1,nu_i),1,rmode,1)
        END DO
!              write(6,*) igroup,iclass, trace(iclass,igroup)
     END DO
  END DO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  IF (nspin_mag==4) THEN
     IF (flag==1) WRITE(stdout,  &
          '(/,5x,"Mode symmetry, ",a11," [",a11,"] magnetic point group:",/)') &
          gname, gname_is
  ELSE
     IF (flag==1) WRITE(stdout,'(/,5x,"Mode symmetry, ",a11," point group:",/)') gname
  END IF
  num_rap_mode=-1
  counter=1
  DO igroup=1,ngroup
     IF (ABS(w1(istart(igroup)))<1.d-3) CYCLE
     DO irap=1,nclass
        times=(0.d0,0.d0)
        DO iclass=1,nclass
           times=times+CONJG(trace(iclass,igroup))*char_mat(irap, &
                which_irr(iclass))*nelem(iclass)
           !         write(6,*) igroup, irap, iclass, which_irr(iclass)
        ENDDO
        times=times/nsym
        cdum="   "
        IF (lri) cdum=ir_ram(irap)
        IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN
           IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') &
                istart(igroup), istart(igroup+1)-1, w1(istart(igroup))
        ENDIF

        IF (ABS(times) > eps) THEN
           IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
              IF (flag==1) WRITE(stdout,'(5x, "omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",a19)') &
                   istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3*nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           ELSE
              IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",i3,a19)') &
                   istart(igroup), istart(igroup+1)-1, &
                   w1(istart(igroup)), NINT(DBLE(times)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(times))*&
                                            NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3 * nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           END IF
        END IF
     END DO
  END DO
  IF (flag==1) WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(z)
  DEALLOCATE(w1)
  DEALLOCATE(rmode)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym

SUBROUTINE rotate_mod(mode,rmode,sr,irt,rtau,xq,nat,irot)
  USE kinds, ONLY : DP
  USE constants, ONLY: tpi
  IMPLICIT NONE

  INTEGER :: nat, irot, irt(48,nat)
  COMPLEX(DP) :: mode(3*nat), rmode(3*nat), phase
  REAL(DP)  :: sr(3,3), rtau(3,48,nat), xq(3), arg
  INTEGER :: na, nb, ipol, kpol, mu_i, mu_k

  rmode=(0.d0,0.d0)
  DO na=1,nat
     nb=irt(irot,na)
     arg = ( xq(1)*rtau(1,irot,na) + xq(2)*rtau(2,irot,na)+  &
          xq(3)*rtau(3,irot,na) ) * tpi
     phase = CMPLX(cos(arg), sin(arg), kind=DP)
     DO ipol=1,3
        mu_i=3*(na-1)+ipol
        DO kpol=1,3
           mu_k=3*(nb-1)+kpol
           rmode(mu_i)=rmode(mu_i) + sr(kpol,ipol)*mode(mu_k)*phase
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
