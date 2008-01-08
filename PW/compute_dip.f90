!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_dip(rho, dip, dipion, z0)
  !
  ! This routine computes the integral 1/Omega \int \rho(r) r d^3r 
  ! and gives as output the projection of the dipole in the direction of
  ! the electric field. (This routine is called only if tefield is true)
  ! The direction is the reciprocal lattice vector bg(.,edir)
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : fpi
  USE ions_base, ONLY : nat, ityp, tau, zv
  USE cell_base, ONLY : alat, at, bg, omega
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,  ONLY : nspin
  USE extfield,  ONLY : edir
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : grid_gather
  
  IMPLICIT NONE
!
! I/O variables
!
  REAL(DP) :: rho(nrxx,nspin)
  REAL(DP) :: dip, dipion,z0
!
! local variables
!
  REAL(DP), ALLOCATABLE :: rrho (:)
#ifdef __PARA
  REAL(DP), ALLOCATABLE :: aux(:), rws(:,:)
  INTEGER  :: nrws, nrwsx
!  REAL(DP) :: wsweight
#endif
  REAL(DP) :: dipol_ion(3), dipol(3)

  INTEGER:: ipol, i, j, k, i1, j1, k11, na, is, ir
  REAL(DP) :: deltax, deltay, deltaz, rijk(3), bmod, proj, x0(3)
  REAL(DP) :: weight
  !
  !  calculate ionic dipole
  !
  x0=0.d0
  x0(3)=-z0
  dipol_ion=0.d0
  DO na=1,nat
     DO ipol=1,3
        dipol_ion(ipol)=dipol_ion(ipol)+zv(ityp(na))*(tau(ipol,na)+x0(ipol))*alat
     ENDDO
  ENDDO
  !
  !  collect the charge density: sum over spin and collect in parallel case
  !
  ALLOCATE (rrho(nrx1*nrx2*nrx3))
  rrho(:) = 0.d0
#ifdef __PARA
  ALLOCATE(aux(nrxx))
  aux(:) =0.d0
  DO is=1,nspin
     aux(:) = aux(:) + rho(:,is)
  ENDDO
  CALL grid_gather (aux, rrho)
  DEALLOCATE(aux)
  IF ((me_pool+1).EQ.1) THEN
#else
     DO is=1,nspin
        rrho=rrho+rho(:,is)
     ENDDO
#endif

!     nrwsx=125
!     allocate(rws(0:3,nrwsx))
!     call wsinit(rws,nrwsx,nrws,at)

     deltax=1.d0/REAL(nr1,dp)
     deltay=1.d0/REAL(nr2,dp)
     deltaz=1.d0/REAL(nr3,dp)
     dipol=0.d0
     DO i1 = -nr1/2, nr1/2
        i=i1+1
        IF (i.LT.1) i=i+nr1
        DO j1 =  -nr2/2,nr2/2
           j=j1+1
           IF (j.LT.1) j=j+nr2
           DO k11 = -nr3/2, nr3/2
              k=k11+1
              IF (k.LT.1) k=k+nr3
              ir=i + (j-1)*nrx1 + (k-1)*nrx1*nrx2
              DO ipol=1,3
                 rijk(ipol) = REAL(i1,dp)*at(ipol,1)*deltax + &
                      REAL(j1,dp)*at(ipol,2)*deltay + &
                      REAL(k11,dp)*at(ipol,3)*deltaz
              ENDDO
              weight = 1.d0
              IF(i1.EQ.-REAL(nr1,dp)/2.d0.OR.i1.EQ.REAL(nr1,dp)/2.d0) &
                 weight = weight*0.5d0
              IF(j1.EQ.-REAL(nr2,dp)/2.d0.OR.j1.EQ.REAL(nr2,dp)/2.d0) &
                 weight = weight*0.5d0
              IF(k11.EQ.-REAL(nr3,dp)/2.d0.OR.k11.EQ.REAL(nr3,dp)/2.d0) &
                 weight = weight*0.5d0
!              weight=wsweight(rijk,rws,nrws)
              DO ipol=1,3
                 dipol(ipol)=dipol(ipol)+weight*(rijk(ipol)+x0(ipol))*rrho(ir)
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     dipol=dipol*alat*omega/nr1/nr2/nr3
     WRITE( stdout,'(5x,"electron", 3f15.5)') dipol(1), dipol(2), dipol(3)
     WRITE( stdout,'(5x,"ion     ", 3f15.5)') dipol_ion(1), dipol_ion(2), dipol_ion(3)
     WRITE( stdout,'(5x,"total   ", 3f15.5)') dipol_ion(1)-dipol(1), &
          dipol_ion(2)-dipol(2), &
          dipol_ion(3)-dipol(3)

     bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)
     proj=0.d0
     DO ipol=1,3
        proj=proj+(dipol_ion(ipol)-dipol(ipol))*bg(ipol,edir)
     ENDDO
     proj=proj/bmod
     dip= fpi*proj/omega

     proj=0.d0
     DO ipol=1,3
        proj=proj+dipol_ion(ipol)*bg(ipol,edir)
     ENDDO
     proj=proj/bmod
     dipion= fpi*proj/omega
!     deallocate(rws)

#ifdef __PARA
  ENDIF
#endif

  DEALLOCATE (rrho)

  RETURN
END SUBROUTINE compute_dip
