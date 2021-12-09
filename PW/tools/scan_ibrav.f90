!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM scan_ibrav
!----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, ANGSTROM_AU
  !USE powell, ONLY : POWELL_MIN
  USE lmdif_module, ONLY : lmdif0
  !
  IMPLICIT NONE
  CHARACTER(len=1024) :: line
  INTEGER,PARAMETER :: npar = 6+3 ! celldm(6) and 3 angles
  INTEGER,PARAMETER :: nfnc = 9
  INTEGER,PARAMETER :: nibrav = 20
  INTEGER,PARAMETER :: ibrav_list(nibrav) =  (/1,2,3,-3,4,5,-5,6,7,8,9,-9,91,10,11,12,-12,13,-13,14/)
  INTEGER :: ibrav, ios, ii, i, info
  REAL(DP) :: celldm(6), angle(3), alat, chisq, chisq_aux
  !
  REAL(DP) :: at(3,3), at_new(3,3), omega, R(3,3), par(npar), par_aux(npar), celldiff(nfnc)
  REAL(DP),PARAMETER :: grad_to_rad = pi/180
!  REAL(DP) :: xi(npar,npar)
  INTEGER :: lwa, iwa(npar)
  REAL(DP),ALLOCATABLE :: wa(:)  
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

  ALLOCATE(wa(lwa))
      
  DO ii = 1, nibrav
    ibrav = ibrav_list(ii)
    !xi = 0._dp
    !FORALL(i=1:npar) xi(i,i) = 1._dp
!    par(1) = SQRT(SUM(at(:,1)**2))
    CALL  at2celldm (ibrav,alat,at(:,1), at(:,2), at(:,3),par)
!    par(1) = SQRT(SUM(at**2))/3
!    par(2) = 1._dp
!    par(3) = 1._dp
!    par(4) = 0.1_dp
!    par(5) = 0.1_dp
!    par(6) = 0.1_dp
    par(7) = 0._dp
    par(8) = 0._dp
    par(9) = 0._dp
    WRITE(*,'(2/,"Scanning ibrav ",i3)') ibrav
    CALL lmdif0(optimize_this_s, nfnc, npar, par, celldiff, 1.d-8, info)
    
    IF(info>0 .and. info<5) THEN
       !PRINT*, "Minimization succeeded"
    ELSEIF(info>=5) THEN
      WRITE(*,'(a)')  "Minimization stopped before convergence"
    ELSEIF(info<=0) THEN 
      WRITE(*,'(a,i6)') "Minimization error", info
      !STOP
    ENDIF
    chisq = SUM(celldiff**2)
      
    IF(chisq<1.d-3) THEN
      WRITE(*,'("Minimization succeeded  (chisq=",g7.1,")")') chisq
      WRITE(*, '("  ibrav = ",i3)') ibrav
      DO i = 1,6
        par_aux = par
        par_aux(i) = par_aux(i) * .5_dp
        !chisq_aux = optimize_this(par_aux)
        CALL optimize_this_s(nfnc, npar, par_aux, celldiff, info)
        chisq_aux = SUM(celldiff**2)
        IF(chisq_aux/=chisq)THEN 
          WRITE(*, '("    celldm(",i2,") = ", f14.9)') i,par(i)
        ENDIF
      ENDDO
      IF(ANY(ABS(par(7:9))>1.d-12))THEN
        WRITE(*,'(a,/,a)') "WARNING! Cell is rotated. Atomic positions will also need",&
                           " to be rotated if they are not in crystal coords!"
        WRITE(*, '("angles (around x,y,z)", 6f14.3)') par(7:9)
      ENDIF
      WRITE(*, '("at1", 6f14.6)') at_new(:,1)
      WRITE(*, '("at2", 6f14.6)') at_new(:,2)
      WRITE(*, '("at3", 6f14.6)') at_new(:,3)
      WRITE(*,*)
    ELSE
      WRITE(*,'("The best cell with this ibrav is not good enough (chisq=",g7.1,")")') chisq
    ENDIF
  ENDDO


 !
 CONTAINS

  SUBROUTINE optimize_this_s(m_,n_,p_,f_,i_)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: m_, n_
    INTEGER,INTENT(inout) :: i_
    REAL(DP),INTENT(out) :: f_(m_)
    REAL(DP),INTENT(in)  :: p_(n_)
    
    REAL(DP) :: celldm_(6), angle_(3), at_(3,3), R(3,3), omega_, pars_(n_), penality
    INTEGER :: ierr, i,j,k
    CHARACTER(len=32) :: errormsg
    ! Global variables from main function: ibrav, at

    pars_ = p_ ! p_ is read only!
    !WRITE(*, '("A:", 6f10.4,3x,3f10.4)') pars_
    !CALL check_bounds(pars_, penality)
    !WRITE(*, '("B:", 6f10.4,3x,3f10.4)') pars_
    
    celldm_ = pars_(1:6)
    angle_  = pars_(7:9)*grad_to_rad   
    
    !WRITE(*, '(6(f14.4,2x),3x,3f10.4)') celldm_, angle_
    CALL latgen_lib( ibrav, celldm_, at_new(:,1), at_new(:,2), at_new(:,3), omega_, ierr, errormsg)
    !IF(ierr/=0) penalty=penalty*10
    !
    IF (ANY(angle_/=0._dp)) THEN
      R = rot(angle_(1), angle_(2), angle_(3)) 
      at_new = matmul(R,at_new)
    ENDIF

    !f_ = 0._dp
    !f_(1) = SUM( (at-at_)**2 )*penality
    k=0
    DO i = 1,3
    DO j = 1,3
      k=k+1
      f_(k) = at(i,j)-at_new(i,j)
    ENDDO
    ENDDO

 END SUBROUTINE
 
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

END PROGRAM scan_ibrav
!----------------------------------------------------------------------

