!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... original code written by Giovanni Cantele and Paolo Cazzato
! ... adapted to work in the parallel case by Carlo Sbraccia
! ... code for the calculation of the vacuum level written by Carlo Sbraccia
! ... code ported from PW to CP by Federico Zipoli
!
SUBROUTINE makov_payne(etot)
!
! CP Modules
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, zv, ityp, ind_srt
  USE ions_positions,    ONLY : tau0
  USE io_global,         ONLY : stdout, ionode, ionode_id
  USE constants,         ONLY : pi, autoev, au_debye
  USE cp_main_variables, ONLY : rhor
  USE electrons_base,    ONLY : nspin
  USE cell_base,         ONLY : at, bg, omega, alat, ibrav
  USE parallel_include
  USE gvecw ,            ONLY : ngw
  USE fft_base,          ONLY : dfftp
#if defined __MPI
  USE mp_global,         ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
  USE mp,                ONLY : mp_barrier
  USE mp_world,          ONLY : world_comm
#endif
!
IMPLICIT NONE
INTEGER :: nspecie
INTEGER :: i,j,k,l,m,n,ip
INTEGER, ALLOCATABLE, DIMENSION(:) :: zvv
REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: rhof
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: r
REAL(DP) :: h(3,3),volumetto,a(3,3)
REAL(DP) :: usunx,usuny,usunz,R0(3),qq,aa,bb
REAL(DP) :: charge, charge_ion, charge_el
REAL(DP) :: dipole(3), dipole_ion(3), dipole_el(3)
REAL(DP) :: quadrupole, quadrupole_ion, quadrupole_el
REAL(DP) :: corr1, corr2, etot
REAL(DP), ALLOCATABLE, DIMENSION(:) :: rgx,rgy,rgz
INTEGER :: ir, is
INTEGER :: ierr
#if defined __MPI
INTEGER :: proc
INTEGER, ALLOCATABLE:: displs(:), recvcount(:)
#endif
REAL(KIND=DP), ALLOCATABLE:: rhodist1(:)
REAL(KIND=DP), ALLOCATABLE:: rhodist2(:)
!
 IF(ibrav.NE.1)THEN
  WRITE(*,*)" WARNING Makov-Payne implemented in CP only when ibrav=1 "
  RETURN
 ENDIF
!
 usunx=1.0D0/DBLE(dfftp%nr1x)
 usuny=1.0D0/DBLE(dfftp%nr2x)
 usunz=1.0D0/DBLE(dfftp%nr3x)
 ALLOCATE ( r(nat,3),rhof(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x),&
  & rgx(dfftp%nr1x),rgy(dfftp%nr2x),rgz(dfftp%nr3x),zvv(nat) )
 !
 DO i=1,nat
  zvv(i)=zv(ityp(ind_srt(i))) 
  DO j=1,3
   r(i,j)=tau0(j,i)
  ENDDO
 ENDDO
 !
 ip=0
 rhof=0.0D0
!
!--------------------------------------------------------------------
ALLOCATE(rhodist1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
IF (nspin.EQ.2) ALLOCATE(rhodist2(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
#if defined __MPI
ALLOCATE( displs( nproc_bgrp ), recvcount( nproc_bgrp ) )
!
      do proc=1,nproc_bgrp
         recvcount(proc) =  dfftp%nnp  * ( dfftp%npp(proc) )
         if (proc.eq.1) then
            displs(proc)=0
         else
            displs(proc)=displs(proc-1) + recvcount(proc-1)
         end if
      end do
!
! gather the charge density on the first node
!
   call mp_barrier( world_comm )
   call mpi_gatherv( rhor(1,1), recvcount(me_bgrp+1), MPI_DOUBLE_PRECISION,&
 &                rhodist1,recvcount, displs, MPI_DOUBLE_PRECISION,&
 &                ionode_id, intra_bgrp_comm, ierr)
   call errore('mpi_gatherv','ierr<>0',ierr)
!
IF(nspin .eq. 2)THEN
         call mp_barrier( world_comm )
         call mpi_gatherv( rhor(1,2), recvcount(me_bgrp+1), MPI_DOUBLE_PRECISION,        &
     &                  rhodist2,recvcount, displs, MPI_DOUBLE_PRECISION,    &
     &                  ionode_id, intra_bgrp_comm, ierr)
         call errore('mpi_gatherv','ierr<>0',ierr)
ENDIF
#else
  rhodist1=rhor(:,1)
 IF(nspin .eq. 2) rhodist2=rhor(:,2)
#endif
!
#if defined __MPI
IF ( ionode ) THEN
#endif
 DO k = 1, dfftp%nr3x
  DO j = 1, dfftp%nr2x
   DO i = 1, dfftp%nr1x
    ip=ip+1
    IF (nspin == 1 )rhof(i,j,k)=rhodist1(ip)
    IF (nspin == 2 )rhof(i,j,k)=rhodist1(ip)+rhodist2(ip)
   ENDDO
  ENDDO
 ENDDO
 ip=0
 DO i=1,dfftp%nr1x
  rgx(i)=DBLE(i-1)*usunx*alat
 ENDDO
 DO i=1,dfftp%nr2x
  rgy(i)=DBLE(i-1)*usuny*alat
 ENDDO
 DO i=1,dfftp%nr3x
  rgz(i)=DBLE(i-1)*usunz*alat
 ENDDO
 !
 !----------------------------------------------------------
 !
 ! center of charge of the ions
 !
 R0=0.0D0
 DO i=1,nat
  R0(1)=R0(1)+zvv(i)*r(i,1)
  R0(2)=R0(2)+zvv(i)*r(i,2)
  R0(3)=R0(3)+zvv(i)*r(i,3)
 ENDDO
  R0=R0/SUM(zvv(1:nat))
 !
 ! shift of the ions (no PBC)
 !
 DO i=1,nat
  r(i,1)=(r(i,1)-R0(1))
  r(i,2)=(r(i,2)-R0(2))
  r(i,3)=(r(i,3)-R0(3))
 ENDDO
 !
 ! shift of the electon density
 !
 DO i=1,dfftp%nr1x
  rgx(i)=(rgx(i)-R0(1))-alat*anint( (rgx(i)-R0(1))/alat )
 ENDDO
 DO i=1,dfftp%nr2x
  rgy(i)=(rgy(i)-R0(2))-alat*anint( (rgy(i)-R0(2))/alat )
 ENDDO
 DO i=1,dfftp%nr3x
  rgz(i)=(rgz(i)-R0(3))-alat*anint( (rgz(i)-R0(3))/alat )
 ENDDO
 
 !
 ! ions
 !
   charge_ion     = SUM(zvv(1:nat))
   dipole_ion     = 0.D0
   quadrupole_ion = 0.D0
 
  DO i = 1, nat
   DO j = 1, 3
    dipole_ion(j) = dipole_ion(j) + zvv(i)*r(i,j)
    quadrupole_ion = quadrupole_ion + zvv(i)*(r(i,j))**2
   ENDDO
  ENDDO
 
 !
 ! electrons
 !
   charge_el = 0.0D0
   dipole_el = 0.0D0
   quadrupole_el = 0.0D0
 
  DO i = 1, dfftp%nr1x
   DO j = 1, dfftp%nr2x
    DO k = 1, dfftp%nr3x
 
     charge_el = charge_el + rhof(i,j,k)
     dipole_el(1) = dipole_el(1) + rgx(i)*rhof(i,j,k)
     dipole_el(2) = dipole_el(2) + rgy(j)*rhof(i,j,k)
     dipole_el(3) = dipole_el(3) + rgz(k)*rhof(i,j,k)
 
     quadrupole_el = quadrupole_el + rhof(i,j,k) * &
   & ( (rgx(i))**2 + (rgy(j))**2 + (rgz(k))**2 )
 
    ENDDO
   ENDDO
  ENDDO
   charge_el=charge_el*alat**3/DBLE(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
   dipole_el=dipole_el*alat**3/DBLE(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
   quadrupole_el=quadrupole_el*alat**3/DBLE(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
 
   ! ... compute ionic+electronic total charge, dipole and quadrupole moments
   !
   charge = -charge_el + charge_ion
   dipole  = -dipole_el + dipole_ion
   quadrupole = -quadrupole_el  + quadrupole_ion
   !
   !
   WRITE( stdout, * )"total charge of the system ",charge
   WRITE( stdout, '(/5X,"charge density inside the ", &
        &               "Wigner-Seitz cell:",3F14.8," el.")' ) charge_el
   !
   WRITE( stdout, &
          '(/5X,"reference position (R0):",5X,3F14.8," bohr")' ) R0(:)
   !
   ! ... A positive dipole goes from the - charge to the + charge.
   !
   WRITE( stdout, '(/5X,"Dipole moments (with respect to x0):")' )
   WRITE( stdout, '( 5X,"Elect",3F10.4," au,  ", 3F10.4," Debye")' ) &
       (-dipole_el(ip), ip = 1, 3), (-dipole_el(ip)*au_debye, ip = 1, 3 )
   WRITE( stdout, '( 5X,"Ionic",3F10.4," au,  ", 3F10.4," Debye")' ) &
       ( dipole_ion(ip),ip = 1, 3), ( dipole_ion(ip)*au_debye,ip = 1, 3 )
   WRITE( stdout, '( 5X,"Total",3F10.4," au,  ", 3F10.4," Debye")' ) &
       ( dipole(ip),    ip = 1, 3), ( dipole(ip)*au_debye,    ip = 1, 3 )
   !
   ! ... print the electronic, ionic and total quadrupole moments
   !
   WRITE( stdout, '(/5X,"Electrons quadrupole moment",F20.8," a.u.")' )  &
       -quadrupole_el
   WRITE( stdout, '( 5X,"     Ions quadrupole moment",F20.8," a.u.")' ) &
       quadrupole_ion
   WRITE( stdout, '( 5X,"    Total quadrupole moment",F20.8," a.u.")' ) &
       quadrupole
   !
   ! ... Makov-Payne correction, PRB 51, 43014 (1995)
   ! ... Note that Eq. 15 has the wrong sign for the quadrupole term
   !
   ! 1 / 2 Ry -> a.u.
   corr1 = - 2.8373D0 / alat * charge**2 / 2.0D0
   !
   aa = quadrupole
   bb = dipole(1)**2 + dipole(2)**2 + dipole(3)**2
   !
   corr2 = ( 2.D0 / 3.D0 * pi )*( charge*aa - bb ) / alat**3
   !
   ! ... print the Makov-Payne correction
   !
   WRITE( stdout, '(/,5X,"*********    MAKOV-PAYNE CORRECTION    *********")' )
   !
   WRITE( stdout,'(/5X,"Makov-Payne correction ",F14.8," a.u. = ",F6.3, &
        &              " eV (1st order, 1/a0)")'   ) -corr1, -corr1*autoev
   WRITE( stdout,'( 5X,"                       ",F14.8," a.u. = ",F6.3, &
        &              " eV (2nd order, 1/a0^3)")' ) -corr2, -corr2*autoev
   WRITE( stdout,'( 5X,"                       ",F14.8," a.u. = ",F6.3, &
        &              " eV (total)")' ) -corr1-corr2, (-corr1-corr2)*autoev
   !
   WRITE( stdout,'(/5X,"corrected Total energy = ",F14.8," a.u.")' ) &
       etot - corr1 - corr2
!
#if defined __MPI
ENDIF ! ionode
#endif
!
IF ( ALLOCATED( rhodist1 ) )   DEALLOCATE( rhodist1 )
IF ( ALLOCATED( rhodist2 ) )   DEALLOCATE( rhodist2 )
#if defined __MPI
IF ( ALLOCATED( displs ) )     DEALLOCATE( displs )
IF ( ALLOCATED( recvcount ) )  DEALLOCATE( recvcount )
#endif
IF ( ALLOCATED( r ) )          DEALLOCATE( r )
IF ( ALLOCATED( rgx ) )        DEALLOCATE( rgx )
IF ( ALLOCATED( rgy ) )        DEALLOCATE( rgy )
IF ( ALLOCATED( rgz ) )        DEALLOCATE( rgz )
IF ( ALLOCATED( zvv ) )        DEALLOCATE( zvv )
IF ( ALLOCATED( rhof ) )       DEALLOCATE( rhof )
!
RETURN
END
