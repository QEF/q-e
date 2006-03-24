!
! Copyright (C) 2001-2005 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
    SUBROUTINE poles ( dipole_moment, dipole_vec, quadrupole, rhortot, &
         ion_flag, tau, coc_flag)
      !------------------------------------------------------------------------
      !
      ! The subroutine computes dipole and quadrupole moments of a charge
      ! distribution 'rhortot'. It was designed originally to be used in 
      ! the context of the Makov-Payne correction. 
      !
      ! If 'ion_flag' is .true., the 
      ! contribution to the dipole due to the nuclei core charge is added.
      !
      ! If coc_flag==.true. the dipole is calculated around the Center Of
      !                     Charge (hence COC)
      !
      ! Note: rhortot id the TOTAL charge density, i.e. spin up + down
      !
      USE cell_base,        ONLY : a1, a2, a3, ainv, omega, alat, at, bg
      USE ions_base,        ONLY : nat, zv, ityp
      USE kinds,            ONLY : DP
      USE electrons_base,   ONLY : qbac
      USE fft_base,         ONLY : dfftp 
      USE mp_global,        ONLY : mpime
      USE mp,               ONLY : mp_sum
#ifdef DFT_CP
      USE grid_dimensions,  ONLY : nr1, nr2, nr3, nr1x, nr2x, nnr=> nnrx
#endif
#ifdef DFT_PW
      USE gvect,            ONLY : nr1, nr2, nr3, nr1x=>nrx1, nr2x=>nrx2, nnr=> nrxx
#endif
      !
      IMPLICIT NONE
      !REAL(KIND=DP), PARAMETER  :: debye=1./0.39344, angs=1./0.52917726
      !
      ! ... input variales
      !
      REAL(KIND=DP), intent(IN)  :: rhortot(nnr), tau(3,nat)
      logical,       intent(IN)  :: ion_flag, coc_flag
      !
      ! ... output variales
      !
      REAL(KIND=DP), intent(OUT) :: dipole_moment,quadrupole,dipole_vec(3)
      !
      ! ... local variales
      !
      REAL(KIND=DP)              :: quad(6), coc(3), tot_charge, ionic_dipole(3)
      REAL(KIND=DP)              :: ax,ay,az,XG0,YG0,ZG0,X,Y,Z,D,rzero,x0,y0,z0
      REAL(KIND=DP)              :: pass1, pass2, pass3, rin(3),rout(3)
      REAL(KIND=DP), ALLOCATABLE :: dip(:)
      INTEGER                    :: ix,ir, i, j, k, me
      !
      me = mpime + 1
      !
      ALLOCATE(dip(nnr))
      tot_charge = 0.d0
      coc = 0.d0
#ifdef DFT_PW
      a1=at(:,1)*alat
      a2=at(:,2)*alat
      a3=at(:,3)*alat
      ainv=transpose(bg)/alat
#endif
      !
      ! ... compute the cores center of charge
      !
      if (coc_flag) then
         do i=1,nat
            coc(:)     = coc(:)+zv(ityp(i))*tau(:,i)
            tot_charge = tot_charge + zv(ityp(i))
         end do
      end if
      coc = coc / tot_charge
      !
      ! ... (1) compute the dipole moment
      !
      ax=a1(1)
      ay=a2(2)
      az=a3(3)
      !
      if (.not.coc_flag) then
         XG0 = -ax/2.
         YG0 = -ay/2.
         ZG0 = -az/2.
      end if
      pass1=ax/nr1
      pass2=ax/nr2
      pass3=ax/nr3
      !
      DO ix=1,3
         ir=1
         !
         print *, dfftp%ipp(me)+1, dfftp%ipp(me)+ dfftp%npp(me)
         DO k = dfftp%ipp(me)+1, dfftp%ipp(me)+ dfftp%npp(me)
            DO j=1,nr2x
               DO i=1,nr1x
                  X=XG0+(i-1)*pass1
                  Y=YG0+(j-1)*pass2
                  Z=ZG0+(k-1)*pass3
                  IF (ix.EQ.1) D=X
                  IF (ix.EQ.2) D=Y
                  IF (ix.EQ.3) D=Z
                  rin(ix) = D-coc(ix) 
                  call vib_pbc(rin,a1,a2,a3,ainv,rout)
                  dip(ir)=rout(ix)*rhortot(ir)
                  ir=ir+1
               END DO
            END DO
         END DO
         !
         dipole_vec(ix)=SUM(dip(1:nnr))
         !
      END DO !!!!!!! ix
      !
      CALL mp_sum(dipole_vec)
      !
      DO ix=1,3
         dipole_vec(ix)=dipole_vec(ix)*omega/DBLE(nr1*nr2*nr3)
      END DO
      !
      IF (ion_flag .and. (.not.coc_flag)) THEN
         !
         ! ... when coc_flag=.true. the ionic contribution
         !     vanishes by definition
         !
         ionic_dipole = 0.d0
         DO ix = 1,nat
            ionic_dipole(:) = ionic_dipole(:) - & 
                 zv(ityp(ix)) * (tau(:,ix)-coc(:))
         END DO
         dipole_vec = - dipole_vec + ionic_dipole
         ! note: electron charge is taken as positive in the code, and thus 
         !       the core is negative. When ion_flag is .true. we reverse the 
         !       sign in the output such that the dipole will have the correct
         !       physical sign.
      END IF
      !
      dipole_moment=SQRT(dipole_vec(1)**2+dipole_vec(2)**2+dipole_vec(3)**2)
      !
      !
      !       compute the coordinates which put the dipole moment to zero
      !
      IF (ABS(qbac).GT.1.d-05) THEN
         x0=dipole_vec(1)/ABS(qbac)
         y0=dipole_vec(2)/ABS(qbac)
         z0=dipole_vec(3)/ABS(qbac)
         rzero=x0**2+y0**2+z0**2
      ELSE
         rzero=0.
      END IF
      !
      ! ... (2) compute the quadrupole moment
      !
      DO ix=1,6
         !
         ir=1
         DO k=dfftp%ipp(me)+1, dfftp%ipp(me) + dfftp%npp(me)
            DO j=1,nr2x
               DO i=1,nr1x
                  !
                  X=XG0+(i-1)*pass1
                  Y=YG0+(j-1)*pass2
                  Z=ZG0+(k-1)*pass3
                  !
                  IF (ix.EQ.1) D=X*X
                  IF (ix.EQ.2) D=Y*Y
                  IF (ix.EQ.3) D=Z*Z
                  IF (ix.EQ.4) D=X*Y
                  IF (ix.EQ.5) D=X*Z
                  IF (ix.EQ.6) D=Y*Z
                  !
                  dip(ir)=D*rhortot(ir)
                  !
                  ir=ir+1
               END DO
            END DO
         END DO
         !
         quad(ix)=SUM(dip(1:nnr))
      END DO
      !
      CALL mp_sum(quad)

      DO ix=1,6
         quad(ix)=quad(ix)*omega/DBLE(nr1*nr2*nr3)
      END DO
      !
      quadrupole=quad(1)+quad(2)+quad(3)-rzero*qbac
      !
      !  only the diagonal elements contribute to the inetaction energy
      !  the term rzero*qbac is subtracted to zero the dipole moment
      !
      WRITE (*,1001)(dipole_vec(ix),ix=1,3)
      WRITE (*,1002) dipole_moment
      WRITE (*,*) ' '
      WRITE (*,1003)(quad(ix),ix=1,3)
      WRITE (*,1004)(quad(ix),ix=4,6)
      WRITE (*,1005) quadrupole,rzero*qbac
      !
1001  FORMAT('DIPOLE XYZ-COMPONENTS (HARTREE A.U.)',f10.4,2x,f10.4,2x,f10.4)
1002  FORMAT('DIPOLE MOMENT         (HARTREE A.U.)',f10.4)
1003  FORMAT('QUADRUPOLE XX-YY-ZZ COMPONENTS (HARTREE A.U.)',             &
           &f9.4,2x,f9.4,2x,f9.4)
1004  FORMAT('QUADRUPOLE XY-XZ-YZ COMPONENTS (HARTREE A.U.)',             &
           &f9.4,2x,f9.4,2x,f9.4)
1005  FORMAT('QUADRUPOLE MOMENT              (HARTREE A.U.)',2f9.4)
      !
      DEALLOCATE(dip)
      !
      RETURN
    END SUBROUTINE poles

