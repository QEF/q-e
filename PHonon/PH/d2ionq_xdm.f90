!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine calculates the XDM contribution to the dynamical matrix.
! It uses the XDM dispersion coefficients and van der Waals radii in
! the prefix.xdm file, which is written by pw.x.
! This code is based on the d2ionq_mm.f90 file by Fabrizio Masullo and
! Paolo Giannozzi.
!
SUBROUTINE d2ionq_xdm(alat,nat,at,tau,q,der2xdm)
  USE kinds, ONLY : DP
  USE io_global, ONLY: ionode, ionode_id, stdout
  USE io_files, ONLY : seqopn, postfix
  USE constants, ONLY : tpi, eps8
  USE mp_images, ONLY : me_image , nproc_image , intra_image_comm
  USE mp, ONLY : mp_sum, mp_bcast
  USE save_ph, ONLY: tmp_dir_save
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nat ! number of atoms in the unit cell
  REAL(DP), INTENT(IN) :: alat ! cell parameter (celldm(1))
  REAL(DP), INTENT(IN) :: tau(3,nat) ! atomic positions in alat units
  REAL(DP), INTENT(IN) :: at(3,3) ! at(:,i) is lattice vector i in alat units
  REAL(DP), INTENT(IN) :: q(3) ! wavevector in 2pi/alat units
  COMPLEX(DP), INTENT(OUT) :: der2xdm(3,nat,3,nat) ! XDM contribution to the dynamical matrix

  INTEGER :: k, l ! cell atom indices (1 -> nat)
  INTEGER :: aa, bb ! coordinate indices (1 -> 3)
  INTEGER :: nr ! lattice vector index (1 -> nvec)
  INTEGER :: first, last, resto, divid ! for parallelization over atoms
  REAL(DP) :: dd, d2 ! atom-atom distance and square distance
  REAL(DP) :: dtau(3) ! cell atom-cell atom vector
  REAL(DP) :: r(3) ! cell atom-env atom vector
  COMPLEX(DP) :: eiqr ! phase factor for the Fourier transform of the dyn. mat.
  ! accumulation auxiliary variables
  REAL(DP) :: auxr
  REAL(DP) :: aux(3,3,nat)
  COMPLEX(DP) :: aux2(3,3,nat)
  REAL(DP) :: g ! pairwise energy contribution g(d), Exdm = -1/2 sum_ij g_ij(d_ij)
  REAL(DP) :: gp ! derivative of g(d) wrt distance
  REAL(DP) :: h ! gp(d) / d
  REAL(DP) :: hp ! derivative of h(d) wrt distance
  ! for reading the xdm.dat file
  INTEGER :: iunxdm, ierr, iver 
  LOGICAL :: lexist
  ! environment info from the XDM module (via .xdm file)
  INTEGER :: nvec ! number of lattice vectors in the environment
  INTEGER :: nat0 ! number of atoms (equal to nat)
  INTEGER, ALLOCATABLE :: lvec(:,:) ! lvec(:,i) is the ith environment lattice vector (cryst. coords.)
  REAL(DP), ALLOCATABLE :: cx(:,:,:) ! cx(i,j,2:4) is nth dispersion coefficient between cell atoms i and j (2=c6,3=c8,4=c10)
  REAL(DP), ALLOCATABLE :: rvdw(:,:) ! rvdw(i,j) is the sum of vdw radii of cell atoms i and j
  REAL(DP) :: rmax2 ! max. distance for energy sum - to be consistent with pw.x
  REAL(DP) :: ene ! total energy (for debug only)

  INTEGER, EXTERNAL :: find_free_unit


  ! write activity message to output
  IF (ionode) THEN
     WRITE (stdout,'(/,5X,"Calculating the XDM contribution to the dynamical matrix.")')
  END IF

  ! read the XDM environment, coefficients, and Rvdw
  IF (ionode) THEN
     iunxdm = find_free_unit ()
     CALL seqopn(iunxdm,postfix(2:6)//'xdm.dat','UNFORMATTED',lexist,tmp_dir_save)
     IF (.NOT.lexist) CALL errore('d2ionq_xdm','could not open xdm data file',1)
     READ (iunxdm,iostat=ierr) iver
     IF (ierr /= 0) CALL errore('d2ionq_xdm','reading xdm.dat 1',1)
     READ (iunxdm,iostat=ierr) nvec, nat0, rmax2
     IF (ierr /= 0) CALL errore('d2ionq_xdm','reading xdm.dat 2',1)
  END IF
  CALL mp_bcast(nvec, ionode_id, intra_image_comm)
  CALL mp_bcast(nat0, ionode_id, intra_image_comm)
  CALL mp_bcast(rmax2, ionode_id, intra_image_comm)
  if (nat /= nat0) CALL errore('d2ionq_xdm','inconsistent number of atoms in .xdm file, nat /= nat0',1)

  ALLOCATE(lvec(3,nvec),cx(nat0,nat0,2:4),rvdw(nat0,nat0))
  IF (ionode) THEN
     READ (iunxdm,iostat=ierr) lvec, cx, rvdw
     IF (ierr /= 0) CALL errore('d2ionq_xdm','reading xdm.dat 3',1)
     CLOSE (UNIT=iunxdm, STATUS='KEEP')
  ENDIF
  CALL mp_bcast(lvec, ionode_id, intra_image_comm)
  CALL mp_bcast(cx, ionode_id, intra_image_comm)
  CALL mp_bcast(rvdw, ionode_id, intra_image_comm)

  ene = 0d0
  der2xdm = 0d0
  ! parallelize atoms over processors in this image
#if defined __MPI
  resto = MOD(nat,nproc_image)
  divid = nat/nproc_image
  IF (me_image+1 <= resto) THEN
     first = (divid+1) * me_image + 1
     last = (divid+1) * (me_image+1)
  ELSE
     first = ((divid+1) * resto) + divid * (me_image-resto) + 1
     last  = (divid+1) * resto + divid * (me_image-resto+1)
  ENDIF
#else
  first = 1
  last = nat
#endif

  DO k = first, last
     aux = 0.d0
     aux2 = 0.d0
     DO l = 1, nat
        dtau = tau(:,k) - tau(:,l)

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1600)
!$omp parallel do private(r,d2,dd,g,gp,h,hp,eiqr,auxr) default(shared), reduction(+:aux), reduction(+:aux2), reduction(+:ene)
#endif
        DO nr = 1, nvec
           r = lvec(1,nr) * at(:,1) + lvec(2,nr) * at(:,2) + lvec(3,nr) * at(:,3) - dtau
           d2 = (r(1)*r(1) + r(2)*r(2) + r(3)*r(3)) * alat * alat
           dd = sqrt(d2)

           IF (dd <= eps8 .OR. d2 > rmax2) CYCLE
           call calcgh_xdm(k,l,dd,g,gp,h,hp)
           ene = ene - 0.5d0 * g

           eiqr = EXP(tpi * (0_dp,1_dp) * (q(1)*(r(1)+dtau(1))+q(2)*(r(2)+dtau(2))+q(3)*(r(3)+dtau(3))))
           DO aa = 1 , 3
              DO bb = 1 , 3
                 IF (aa /= bb) THEN
                    auxr = hp * r(aa) * alat * r(bb) * alat / dd
                 ELSE
                    auxr = hp * r(aa) * alat * r(bb) * alat / dd + h
                 ENDIF
                 aux(aa,bb,l) = aux(aa,bb,l) + auxr
                 aux2(aa,bb,l) = aux2(aa,bb,l) + auxr*eiqr
              ENDDO ! bb
           ENDDO ! aa
        ENDDO ! nr
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1600)
!$omp end parallel do
#endif
        DO aa =1, 3
           DO bb = 1, 3
              der2xdm(aa,k,bb,l) = aux2(aa,bb,l)
           ENDDO ! bb
        ENDDO ! aa
     ENDDO ! l

     DO l = 1, nat
        DO aa = 1, 3
           DO bb = 1, 3
              der2xdm(aa,k,bb,k) = der2xdm(aa,k,bb,k) - aux(aa,bb,l)
           ENDDO ! bb
        ENDDO ! aa
     ENDDO ! l
  ENDDO ! k
  
  CALL mp_sum(ene, intra_image_comm)
  CALL mp_sum(der2xdm, intra_image_comm)
  
  IF (ionode) THEN
     WRITE (stdout,'(5X,"XDM energy = ",F17.8," Ry")') ene
     WRITE (stdout,'(5X,"Done."/)')
  END IF

CONTAINS

  SUBROUTINE calcgh_xdm(i,j,d,g,gp,h,hp)
    IMPLICIT NONE
    integer, intent(in) :: i, j
    real*8, intent(in) :: d
    real*8, intent(out) :: g, gp, h, hp

    real*8 :: c6, c8, c10, rr
    real*8 :: d2, d4, d6, d8, d10, dpr6, dpr8, dpr10, r2, r6, r8, r10
    real*8 :: d5, d7, d9, dpr6sq, dpr8sq, dpr10sq, d17, d13, d3, dpr6cub
    real*8 :: dpr8cub, dpr10cub

    c6 = cx(i,j,2)
    c8 = cx(i,j,3)
    c10 = cx(i,j,4)
    rr = rvdw(i,j)

    r2 = rr * rr
    r6 = r2 * r2 * r2
    r8 = r6 * r2
    r10 = r8 * r2

    d2 = d * d
    d3 = d2 * d
    d4 = d2 * d2
    d5 = d4 * d
    d6 = d4 * d2
    d7 = d6 * d
    d8 = d6 * d2
    d9 = d8 * d
    d10 = d8 * d2
    d13 = d6 * d7
    d17 = d10 * d7

    dpr6 = r6 + d6
    dpr8 = r8 + d8
    dpr10 = r10 + d10
    dpr6sq = dpr6 * dpr6
    dpr8sq = dpr8 * dpr8
    dpr10sq = dpr10 * dpr10
    dpr6cub = dpr6sq * dpr6
    dpr8cub = dpr8sq * dpr8
    dpr10cub = dpr10sq * dpr10

    g = c6 / dpr6 + c8 / dpr8 + c10 / dpr10
    gp = -10d0 * c10 * d9 / dpr10sq - 8d0 * c8 * d7 / dpr8sq - 6d0 * c6 * d5 / dpr6sq
    h = gp / d
    hp = -80d0 * c10 * d7 / dpr10sq + 200d0 * c10 * d17 / dpr10cub - 48d0 * c8 * d5 / dpr8sq &
       + 128d0 * c8 * d13 / dpr8cub - 24 * c6 * d3 / dpr6sq + 72d0 * c6 * d9 / dpr6cub

  END SUBROUTINE calcgh_xdm

END SUBROUTINE d2ionq_xdm

