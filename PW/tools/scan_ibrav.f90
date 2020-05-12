!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM ibrav2cell
!----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, ANGSTROM_AU
  USE powell, ONLY : POWELL_MIN
  !
  IMPLICIT NONE
  CHARACTER(len=1024) :: line
  INTEGER,PARAMETER :: npar = 6+3 ! celldm(6) and 3 angles
  INTEGER,PARAMETER :: nibrav = 20
  INTEGER,PARAMETER :: ibrav_list(nibrav) =  (/1,2,3,-3,4,5,-5,6,7,8,9,-9,91,10,11,12,-12,13,-13,14/)
  INTEGER :: ibrav, ios, ii, i
  REAL(DP) :: celldm(6), angle(3), alat, chisq, chisq_aux
  !
  REAL(DP) :: at(3,3), omega, R(3,3), par(npar), par_aux(npar)
  REAL(DP),PARAMETER :: grad_to_rad = pi/180
  REAL(DP) :: xi(npar,npar)
  !
  LOGICAL,EXTERNAL :: matches

  WRITE(*,*) "Enter the unit of measur (angstrom, bohr) or alat in bohr units, or alat in Angstrom units, followed by ' A'"
  READ(*,"(a1024)") line
  IF(matches("angstrom",line)) THEN
    alat=ANGSTROM_AU
  ELSE IF(matches("bohr",line)) THEN
    alat=1._dp
  ELSE
    READ(line,*,iostat=ios) alat
    IF(ios/=0) CALL errore("scan_ibrav","Could not understand alat '"//TRIM(line)//"'",1)
    IF(matches("A",line)) alat=alat*ANGSTROM_AU
  ENDIF
  WRITE(*,'("alat (bohr)",f12.6)') alat
  WRITE(*,*) "Enter the cell basis vectors (one per line)"
  DO i = 1,3
    READ(*,*) at(:,i)
  ENDDO
  at = at*alat
  WRITE(*,'("Requested axes in Bohr units:")') 
  WRITE(*, '("at1", 6f14.6)') at(:,1)
  WRITE(*, '("at2", 6f14.6)') at(:,2)
  WRITE(*, '("at3", 6f14.6)') at(:,3)

  DO ii = 1, nibrav
    ibrav = ibrav_list(ii)
    xi = 0._dp
    FORALL(i=1:npar) xi(i,i) = 1._dp
    par(1) = SQRT(SUM(at**2))/3
    par(2) = 1._dp
    par(3) = 1._dp
    par(4) = 0.1_dp
    par(5) = 0.1_dp
    par(6) = 0.1_dp
    par(7) = 0._dp
    par(8) = 0._dp
    par(9) = 0._dp
    CALL POWELL_MIN(optimize_this,par,xi,npar,npar,1.d-36,i,chisq)
    IF(chisq<1.d-6)THEN
      WRITE(*,'("____________ MATCH (chisq=",g7.1,") ____________")') chisq
      WRITE(*, '("  ibrav = ",i3)') ibrav
      DO i = 1,6
        par_aux = par
        par_aux(i) = par_aux(i) * .5_dp
        chisq_aux = optimize_this(par_aux)
        IF(chisq_aux/=chisq)THEN 
          WRITE(*, '("    celldm(",i2,") = ", f14.9)') i,par(i)
        ENDIF
      ENDDO
      IF(ANY(ABS(par(7:9))>1.d-12))THEN
        WRITE(*, '("angles", 6f14.3)') par(7:9)
      ENDIF
      WRITE(*, '("at1", 6f14.6)') at(:,1)
      WRITE(*, '("at2", 6f14.6)') at(:,2)
      WRITE(*, '("at3", 6f14.6)') at(:,3)
    ENDIF
  ENDDO


 !
 CONTAINS

 REAL(DP) FUNCTION optimize_this(parameters_) RESULT(chi2)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: parameters_(npar)
    REAL(DP) :: celldm_(6), angle_(3), at_(3,3), R(3,3), omega_, pars_(npar), penality
    ! Global variables from main function: ibrav, at

    pars_ = parameters_ ! parameters_ is read only!
    !WRITE(*, '("A:", 6f10.4,3x,3f10.4)') pars_
    CALL check_bounds(pars_, penality)
    !WRITE(*, '("B:", 6f10.4,3x,3f10.4)') pars_
    
    celldm_ = pars_(1:6)
    angle_  = pars_(7:9)*grad_to_rad   
    
    !WRITE(*, '(6(f14.4,2x),3x,3f10.4)') celldm_, angle_
    CALL latgen_internal( ibrav, celldm_, at_(:,1), at_(:,2), at_(:,3), omega_)
    !
    IF (ANY(angle_/=0._dp)) THEN
      R = rot(angle_(1), angle_(2), angle_(3)) 
      at_ = matmul(R,at_)
    ENDIF

    chi2 = SUM( (at-at_)**2 )*penality
    !WRITE(*, '(3(3f10.4,3x))') at

 END FUNCTION

 SUBROUTINE check_bounds(pars_, penalty)
   IMPLICIT NONE
   REAL(DP),INTENT(inout) :: pars_(npar), penalty
   REAL(DP),PARAMETER :: infty = 1.d+100, eps=1.d-6
   REAL(DP),PARAMETER :: par_min(npar) = (/ eps, eps, eps, -.5_dp+eps, -.5_dp+eps, -1._dp, -180._dp, -180._dp, -180._dp /)
   REAL(DP),PARAMETER :: par_max(npar) = (/ infty, infty, infty,  1._dp-eps,  1._dp-eps,  1._dp-eps,  180._dp,  180._dp,  180._dp /)
   INTEGER :: i
   penalty = 1._dp
   DO i = 1, npar
     IF(pars_(i) < par_min(i)) THEN
        pars_(i) = par_min(i)
        penalty=penalty*10
     ENDIF
     IF(pars_(i) > par_max(i)) THEN
        pars_(i) = par_max(i)
        penalty=penalty*10
     ENDIF
   ENDDO
 END SUBROUTINE

 function rotx (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/ 1._dp,          0._dp,            0._dp/)
   R(:,2) = (/ 0._dp, cos(theta), -sin(theta) /)
   R(:,3) = (/ 0._dp, sin(theta), cos(theta)  /)
 endfunction
 function roty (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/ cos(theta), 0._dp, sin(theta) /)
   R(:,2) = (/ 0._dp,1._dp,0._dp/)
   R(:,3) = (/ -sin(theta), 0._dp, cos(theta)/)
 endfunction
 function rotz (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/cos(theta),-sin(theta), 0._dp /)
   R(:,2) = (/sin(theta), cos(theta), 0._dp /)
   R(:,3) = (/0._dp,0._dp,1._dp /)
 endfunction
 function rot(alpha,beta,gamma) RESULT(R)
   IMPLICIT NONE
   REAL(DP),INTENT(in) :: alpha,beta,gamma
   REAL(DP) :: R(3,3)
   R = matmul(matmul(rotx(alpha),roty(beta)), rotz(gamma))
 endfunction

END PROGRAM ibrav2cell
!----------------------------------------------------------------------




!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE latgen_internal(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  real(DP), INTENT(inout) :: celldm(6)
  real(DP), INTENT(inout) :: a1(3), a2(3), a3(3)
  real(DP), INTENT(out) :: omega
  !
  real(DP), PARAMETER:: sr2 = 1.414213562373d0, &
                        sr3 = 1.732050807569d0
  INTEGER :: i,j,k,l,iperm,ir
  real(DP) :: term, cbya, s, term1, term2, singam, sen
  !
  !  user-supplied lattice vectors
  !
  IF (ibrav == 0) THEN
     IF (sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  &
         THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  &
         THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  &
         THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF

     IF ( celldm(1) /= 0.D0 ) THEN
     !
     ! ... input at are in units of alat => convert them to a.u.
     !
         a1(:) = a1(:) * celldm(1)
         a2(:) = a2(:) * celldm(1)
         a3(:) = a3(:) * celldm(1)
     ELSE
     !
     ! ... input at are in atomic units: define celldm(1) from a1
     !
         celldm(1) = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
     ENDIF
     !
  ELSE
     a1(:) = 0.d0
     a2(:) = 0.d0
     a3(:) = 0.d0
  ENDIF
  !
  IF (celldm (1) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
  !
  !  index of bravais lattice supplied
  !
  IF (ibrav == 1) THEN
     !
     !     simple cubic lattice
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  ELSEIF (ibrav == 2) THEN
     !
     !     fcc lattice
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  ELSEIF (abs(ibrav) == 3) THEN
     !
     !     bcc lattice
     !
     term=celldm(1)/2.d0
     DO ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     ENDDO
     IF ( ibrav < 0 ) THEN
        a1(1)=-a1(1)
        a2(2)=-a2(2)
        a3(3)=-a3(3)
     ELSE
        a2(1)=-a2(1)
        a3(1)=-a3(1)
        a3(2)=-a3(2)
     ENDIF
     !
  ELSEIF (ibrav == 4) THEN
     !
     !     hexagonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (abs(ibrav) == 5) THEN
     !
     !     trigonal lattice
     !
     IF (celldm (4) <= -0.5_dp .or. celldm (4) >= 1.0_dp) &
          THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     term1=sqrt(1.0_dp + 2.0_dp*celldm(4))
     term2=sqrt(1.0_dp - celldm(4))
     !
     IF ( ibrav == 5) THEN
        !     threefold axis along c (001)
        a2(2)=sr2*celldm(1)*term2/sr3
        a2(3)=celldm(1)*term1/sr3
        a1(1)=celldm(1)*term2/sr2
        a1(2)=-a1(1)/sr3
        a1(3)= a2(3)
        a3(1)=-a1(1)
        a3(2)= a1(2)
        a3(3)= a2(3)
     ELSEIF ( ibrav == -5) THEN
        !     threefold axis along (111)
        ! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
        ! does not yield the x,y,z axis, but an equivalent rotated triplet:
        !    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
        ! If you prefer the x,y,z axis as cubic limit, you should modify the
        ! definitions of a1(1) and a1(2) as follows:'
        !    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
        !    a1(2) = celldm(1)*(term1-term2)/3.0_dp
        ! (info by G. Pizzi and A. Cepellotti)
        !
        a1(1) = celldm(1)*(term1-2.0_dp*term2)/3.0_dp
        a1(2) = celldm(1)*(term1+term2)/3.0_dp
        a1(3) = a1(2)
        a2(1) = a1(3)
        a2(2) = a1(1)
        a2(3) = a1(2)
        a3(1) = a1(2)
        a3(2) = a1(3)
        a3(3) = a1(1)
     ENDIF
  ELSEIF (ibrav == 6) THEN
     !
     !     tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (ibrav == 7) THEN
     !
     !     body centered tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  ELSEIF (ibrav == 8) THEN
     !
     !     Simple orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF ( abs(ibrav) == 9) THEN
     !
     !     One face (base) centered orthorhombic lattice  (C type)
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     IF ( ibrav == 9 ) THEN
        !   old PWscf description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) = a1(1) * celldm(2)
        a2(1) = - a1(1)
        a2(2) = a1(2)
     ELSE
        !   alternate description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) =-a1(1) * celldm(2)
        a2(1) = a1(1)
        a2(2) =-a1(2)
     ENDIF
     a3(3) = celldm(1) * celldm(3)
     !
  ELSEIF ( ibrav == 91 ) THEN
     !
     !     One face (base) centered orthorhombic lattice  (A type)
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     a1(1) = celldm(1)
     a2(2) = celldm(1) * celldm(2) * 0.5_DP
     a2(3) = - celldm(1) * celldm(3) * 0.5_DP
     a3(2) = a2(2)
     a3(3) = - a2(3)
     !
  ELSEIF (ibrav == 10) THEN
     !
     !     All face centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     a2(1) = 0.5d0 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 11) THEN
     !
     !     Body centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 12) THEN
     !
     !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF (ibrav ==-12) THEN
     !
     !     Simple monoclinic lattice, unique axis: b (more common)
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     sen=sqrt(1.d0-celldm(5)**2)
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(3)=celldm(1)*celldm(3)*sen
     !
  ELSEIF (ibrav == 13) THEN
     !
     !     One face centered monoclinic lattice unique axis c
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(3) =-a1(1) * celldm(3)
     a2(1) = celldm(1) * celldm(2) * celldm(4)
     a2(2) = celldm(1) * celldm(2) * sen
     a3(1) = a1(1)
     a3(3) =-a1(3)
  ELSEIF (ibrav == -13) THEN
!     CALL infomsg('latgen','BEWARE: axis for ibrav=-13 changed, see documentation!')
     !
     !     One face centered monoclinic lattice unique axis b
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     sen = sqrt( 1.d0 - celldm(5) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a2(1) =-a1(1)
     a2(2) = a1(2)
     a3(1) = celldm(1) * celldm(3) * celldm(5)
     a3(3) = celldm(1) * celldm(3) * sen
     !
  ELSEIF (ibrav == 14) THEN
     !
     !     Triclinic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     IF (abs(celldm(6))>=1.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)
     IF (term < 0.d0) THEN; a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN; ENDIF
     term= sqrt(term/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  ELSE
     !
     a1=0._dp; a2=0._dp; a3=0._dp; omega=0._dp; RETURN
     !
  ENDIF
  !
  !  calculate unit-cell volume omega
  !
  CALL volume (1.0_dp, a1, a2, a3, omega)
  !
  RETURN
  !
END SUBROUTINE latgen_internal
!
