!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dyndiar (dyn,nat3,nmodes,u,nat,ityp,amass,w2,dynout)
  !-----------------------------------------------------------------------
  !
  !   diagonalizes the dynamical matrix "dyn", returns energies in "w2"
  !   and mode displacements in "dynout". dyn is unchanged on output.
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, ry_to_thz, ry_to_cmm1
  USE io_global,  ONLY : stdout
  IMPLICIT NONE
  INTEGER :: nmodes, nat3, nat,ityp(nat), iudyn
  real(DP):: dyn(nat3,nmodes), u(nat3,nmodes), amass(*)
  real(DP):: dynout(nat3,nmodes), w2(nat3)
  !
  INTEGER:: nu_i, nu_j, mu, na, nb, nt, i, j
  real(DP), ALLOCATABLE :: m(:,:), z(:,:)
  real(DP) :: w1, unorm, sum, dif
  !
  ALLOCATE  ( m  ( nat3, nat3))
  ALLOCATE  ( z  ( nat3, nat3))
  !
  CALL dcopy(nat3*nmodes,dyn,1,dynout,1)
  !
  !  Impose symmetry to the matrix
  !
  dif=0.d0
  DO nu_i=1,nmodes
     DO nu_j=1,nu_i-1
        dif = dif + abs(dynout(nu_i,nu_j)-dynout(nu_j,nu_i))
        dynout(nu_j,nu_i) = 0.5d0*(dynout(nu_i,nu_j)+dynout(nu_j,nu_i))
        dynout(nu_i,nu_j) = dynout(nu_j,nu_i)
     ENDDO
  ENDDO
  WRITE( stdout,9000) dif
  !
  !  Impose Acoustic Sum Rule
  !
  dif=0.d0
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           sum=0.d0
           DO nb=1,nat
              IF (na/=nb) sum=sum+dynout((na-1)*3+i,(nb-1)*3+j)
           ENDDO
           dif = dif + abs(dynout((na-1)*3+i,(na-1)*3+j) + sum)
           dynout((na-1)*3+i,(na-1)*3+j) = -sum
        ENDDO
     ENDDO
  ENDDO
  WRITE( stdout,9005) dif
  !
  !  fill the mass matrix (masses are in amu, amu_ry converts to a.u.)
  !
  DO nu_i = 1,nmodes
     DO nu_j = 1,nmodes
        m(nu_i,nu_j) = 0.0d0
        DO mu = 1,3*nat
           na = (mu-1)/3+1
           nt = ityp(na)
           m(nu_i,nu_j) = m(nu_i,nu_j) + amu_ry*amass(nt)*u(mu,nu_i)*u(mu,nu_j)
        ENDDO
     ENDDO
  ENDDO
  !
  !  solve the generalized eigenvalue problem w2*(M*z) = (Cz)
  !  Note that z are eigendisplacements in the base of input
  !  modes u and that they are normalized as <z|M|z>=I
  !
  CALL rdiaghg (nat3, nmodes, dynout, m, nat3, w2, z)
  !
  !  write frequencies
  !
  WRITE( stdout,'(5x,"diagonalizing the dynamical matrix ..."//)')
  WRITE( stdout,'(1x,74("*"))')
  !
  dynout (:,:) = 0.0d0
  DO nu_i = 1,nmodes
     w1 = sqrt(abs(w2(nu_i)))
     IF (w2(nu_i)<0.0) w1 = -w1
     WRITE( stdout,9010) nu_i, w1*ry_to_thz, w1*ry_to_cmm1
     !  bring eigendisplacements in cartesian axis
     DO mu = 1,3*nat
        DO i = 1,nmodes
           dynout(mu,nu_i) = dynout(mu,nu_i) + z(i,nu_i)*u(mu,i)
        ENDDO
     ENDDO
  ENDDO
  WRITE( stdout,'(1x,74("*"))')
  !
  DEALLOCATE(z)
  DEALLOCATE(m)
  RETURN
  !
9000 FORMAT ('  Symmetry violation  sum_ij |D_ij-D_ji| :',f15.6)
9005 FORMAT ('  ASR violation  sum_i |D_ij| :',f15.6)
9010 FORMAT(5x,'omega(',i3,') =',f10.6,' [THz] =',f11.6,' [cm-1]')
  !
END SUBROUTINE dyndiar
