!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE d2ionq_dispd3( alat, nat, ityp, at, bg, tau, q, der2disp )
  !------------------------------------------------------------------------
  !! This routine calculates the Grimme-D3 contribution to the dynamical matrix.
  !
  USE kinds,            ONLY: DP
  USE io_global,        ONLY: stdout
  USE symme,            ONLY: symvector
  USE mp,               ONLY: mp_stop
  USE input_parameters, ONLY: dftd3_version, dftd3_threebody
  USE funct,            ONLY: get_dft_short
  USE dftd3_api,        ONLY: dftd3_init, dftd3_set_functional, dftd3_pbc_dispersion, get_atomic_number
  USE dftd3_qe,         ONLY: dftd3_xc, dftd3, dftd3_in
  USE ions_base,        ONLY: atm
  USE ener,             ONLY: edftd3
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the unit cell
  REAL(DP), INTENT(IN) :: alat
  !! cell parameter (celldm(1))
  INTEGER, INTENT(IN) :: ityp(nat)
  !! atomic types for atoms in the unit cell
  REAL(DP), INTENT(IN) :: at(3,3)
  !! at(:,i) is lattice vector i in alat units
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! bg(:,i) is reciprocal lattice vector i in 2pi/alat units
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions in alat units
  REAL(DP), INTENT(IN) :: q(3)
  !! wavevector in 2pi/alat units
  COMPLEX(DP), INTENT(OUT) :: der2disp(3,nat,3,nat)
  !! dispersion contribution to the (massless) dynamical matrix
  ! 
  !  Local variables  
  CHARACTER(LEN=256):: dft_ , formt
  REAL(DP) :: latvecs(3,3)
  REAL(DP) :: step, eerr, eerl, eelr, eell
  INTEGER:: atnum(1:nat), iat, ixyz, jat, jxyz
  REAL(DP) :: xyz(3,nat)
  REAL(DP) :: stress_d3(3,3)
  REAL(DP), ALLOCATABLE :: force_d3(:,:), force_num(:,:)
  REAL(DP), ALLOCATABLE :: der2disp_ene(:,:,:,:), der2disp_frc(:,:,:,:) 
  
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)

  write(stdout,'(A)') 'Adding Grimme-D3 contribution to the dynamical matrix'
  WRITE ( stdout , 9078 ) edftd3 
  !
  ! Setting DFT-D3 functional dependent parameters
  !
  if (dftd3_version==2) dftd3_threebody=.false.
  dftd3_in%threebody = dftd3_threebody
  CALL dftd3_init(dftd3, dftd3_in)
  dft_ = get_dft_short( )
  dft_ = dftd3_xc ( dft_ )
  CALL dftd3_set_functional(dftd3, func=dft_, version=dftd3_version,tz=.false.)

  CALL start_clock('force_dftd3')
  ALLOCATE( force_d3(3, nat), force_num(3, nat) )
  force_num = 0._dp
  force_d3(:,:) = 0.0_DP
  latvecs(:,:)=at(:,:)*alat
  xyz(:,:)=tau(:,:)*alat
  DO iat = 1, nat
     atnum(iat) = get_atomic_number(TRIM(atm(ityp(iat))))
  ENDDO
  call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, edftd3, force_d3, stress_d3)
  edftd3=edftd3*2.d0
  force_d3 = -2.d0*force_d3
  
  step=1.d-6
  !step=2.d-5

  do iat = 1, nat
    do ixyz = 1, 3 
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerr)
      eerr = eerr*2.0d0
      xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eell)
      eell = eell*2.0d0
      force_num(ixyz,iat)=-0.5*(eerr-eell)/step
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
    end do 
  end do 
  !CALL symvector( nat, force_num )
  CALL stop_clock('force_dftd3')

  WRITE ( stdout , 9078 ) edftd3  
  WRITE( stdout, '(/,5x,"DFT-D3 dispersion contribution to forces (analytical):")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_d3(ixyz,iat), ixyz = 1, 3)
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 dispersion contribution to forces (numerical):")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_num(ixyz,iat), ixyz = 1, 3)
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 analytical vs numerical err:")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_d3(ixyz,iat)-force_num(ixyz,iat), ixyz = 1, 3)
  ENDDO

  step=2.d-5
  CALL start_clock('hessian_dftd3')
  ALLOCATE( der2disp_ene(3,nat,3,nat), der2disp_frc(3,nat,3,nat) )
  der2disp = 0._dp
  der2disp_ene = 0._dp
  der2disp_frc = 0._dp
  do iat = 1, nat
    do ixyz = 1, 3 
      write(*,*) 'displacing forces: ', iat, ixyz

      xyz(ixyz,iat)=xyz(ixyz,iat)+step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerr, force_d3, stress_d3)
      force_d3 = -2.d0*force_d3
      xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eell, force_num, stress_d3)
      force_num = -2.d0*force_num

      der2disp_frc(ixyz,iat,1:3,1:nat) = -0.5 * (force_d3(1:3,1:nat) - force_num(1:3,1:nat) ) / step
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
    end do 
  end do 

  do iat = 1, nat
    do ixyz = 1, 3 
      do jat = 1, nat
        do jxyz = 1, 3 

          write(*,*) 'displacing energy: ', iat, ixyz, jat, jxyz

          xyz(ixyz,iat)=xyz(ixyz,iat)+step
          xyz(jxyz,jat)=xyz(jxyz,jat)+step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerr)  ! i+s  j+s
          eerr = eerr*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eelr)  ! i-s  j+s
          eelr = eelr*2.0d0

          xyz(jxyz,jat)=xyz(jxyz,jat)-2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eell)  ! i-s  j-s
          eell = eell*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)+2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerl)  ! i+s  j-s
          eerl = eerl*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)-step
          xyz(jxyz,jat)=xyz(jxyz,jat)+step

          der2disp_ene(ixyz,iat,jxyz,jat) = (eerr - eerl - eelr + eell) / 4.0d0 / step / step

        end do 
      end do 
    end do 
  end do 
  CALL stop_clock('hessian_dftd3')

  write(formt,'(A,I5,A)') '(',3*nat,'F12.8)'
  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Forces):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_frc(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Energies):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_ene(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO


  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Forces vs from Energies):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_frc(ixyz,iat,1:3,1:nat)-der2disp_ene(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO


  DEALLOCATE( force_d3, force_num, der2disp_ene, der2disp_frc)

  Call mp_stop(555)

  RETURN

END SUBROUTINE d2ionq_dispd3
!---------------------------------------------------------------------------
SUBROUTINE d2ionq_disp( alat, nat, ityp, at, bg, tau, q, der2disp )
  !------------------------------------------------------------------------
  !! This routine calculates the XDM contribution to the dynamical matrix.
  !! It uses the XDM dispersion coefficients and Van der Waals radii in
  !! the \(\texttt{prefix.xdm}\) file, which is written by \(\texttt{pw.x}\).
  !
  !! This code is based on the \(\texttt{d2ionq_mm.f90}\) file by Fabrizio
  !! Masullo and Paolo Giannozzi.
  !
  USE london_module, ONLY: init_london, dealloca_london, mxr, dist2, r_cut, r
  USE kinds,         ONLY: DP
  USE io_global,     ONLY: ionode, ionode_id, stdout
  USE io_files,      ONLY: seqopn, postfix
  USE control_flags, ONLY: llondon, lxdm
  USE constants,     ONLY: tpi, eps8
  USE mp_images,     ONLY: me_image , nproc_image , intra_image_comm
  USE mp,            ONLY: mp_sum, mp_bcast
  USE save_ph,       ONLY: tmp_dir_save
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the unit cell
  REAL(DP), INTENT(IN) :: alat
  !! cell parameter (celldm(1))
  INTEGER, INTENT(IN) :: ityp(nat)
  !! atomic types for atoms in the unit cell
  REAL(DP), INTENT(IN) :: at(3,3)
  !! at(:,i) is lattice vector i in alat units
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! bg(:,i) is reciprocal lattice vector i in 2pi/alat units
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions in alat units
  REAL(DP), INTENT(IN) :: q(3)
  !! wavevector in 2pi/alat units
  COMPLEX(DP), INTENT(OUT) :: der2disp(3,nat,3,nat)
  !! dispersion contribution to the (massless) dynamical matrix
  
  ! ... local variables
  
  INTEGER :: ii, jj, kk ! some indices
  INTEGER :: k, l ! cell atom indices (1 -> nat)
  INTEGER :: aa, bb ! coordinate indices (1 -> 3)
  INTEGER :: nr ! lattice vector index (1 -> nvec)
  INTEGER :: first, last, resto, divid ! for parallelization over atoms
  REAL(DP) :: dd, d2 ! atom-atom distance and square distance
  REAL(DP) :: dtau(3) ! cell atom-cell atom vector
  REAL(DP) :: rr(3) ! cell atom-env atom vector
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
  INTEGER :: lmax(3) ! max. lattice vector in the environment for each crystallographic axis
  INTEGER, ALLOCATABLE :: lvec(:,:) ! lvec(:,i) is the ith environment lattice vector (cryst. coords.)
  REAL(DP), ALLOCATABLE :: cx(:,:,:) ! cx(i,j,2:4) is nth dispersion coefficient between cell atoms i and j (2=c6,3=c8,4=c10)
  REAL(DP), ALLOCATABLE :: rvdw(:,:) ! rvdw(i,j) is the sum of vdw radii of cell atoms i and j
  REAL(DP) :: rmax2 ! max. distance for energy sum - to be consistent with pw.x
  REAL(DP) :: ene ! total energy (for debug only)
  CHARACTER*10 :: namestr ! name of the dispersion correction

  INTEGER, EXTERNAL :: find_free_unit


  ! initialization
  IF (llondon) THEN
     ! D2 does not require any saved info; just init the module
     CALL init_london()
     namestr = "D2"

  ELSE IF (lxdm) THEN
     ! read the XDM environment, coefficients, and Rvdw
     ALLOCATE(cx(nat,nat,2:4),rvdw(nat,nat))
     IF (ionode) THEN
        iunxdm = find_free_unit ()
        CALL seqopn(iunxdm,postfix(2:6)//'xdm.dat','UNFORMATTED',lexist,tmp_dir_save)
        IF (.NOT.lexist) CALL errore('d2ionq_disp','could not open xdm data file',1)
        READ (iunxdm,iostat=ierr) iver
        IF (ierr /= 0) CALL errore('d2ionq_disp','reading xdm.dat 1',1)
        READ (iunxdm,iostat=ierr) lmax, rmax2
        IF (ierr /= 0) CALL errore('d2ionq_disp','reading xdm.dat 2',1)
        READ (iunxdm,iostat=ierr) cx, rvdw
        IF (ierr /= 0) CALL errore('d2ionq_disp','reading xdm.dat 3',1)
        CLOSE (UNIT=iunxdm, STATUS='KEEP')
     END IF
     CALL mp_bcast(iver, ionode_id, intra_image_comm)
     CALL mp_bcast(lmax, ionode_id, intra_image_comm)
     CALL mp_bcast(rmax2, ionode_id, intra_image_comm)
     CALL mp_bcast(cx, ionode_id, intra_image_comm)
     CALL mp_bcast(rvdw, ionode_id, intra_image_comm)
     namestr = "XDM"

     ! pre-calculate the list of lattice vectors
     ALLOCATE(lvec(3,PRODUCT(2*lmax + 1)))
     nvec = 0
     DO ii = -lmax(1), lmax(1)
        DO jj = -lmax(2), lmax(2)
           DO kk = -lmax(3), lmax(3)
              nvec = nvec + 1
              lvec(:,nvec) = (/ii,jj,kk/)
           END DO
        END DO
     END DO
  ELSE
     CALL errore('d2ionq_disp','Dispersion correction not one of D2 or XDM',1)
  ENDIF

  IF (ionode) THEN
     WRITE (stdout,'(/,5X,"Calculating the ",A," contribution to the dynamical matrix.")') TRIM(namestr)
  END IF

  ene = 0._dp
  der2disp = 0._dp
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
     aux = 0._dp
     aux2 = 0._dp
     DO l = 1, nat
        dtau = tau(:,k) - tau(:,l)

        IF (llondon) THEN
           CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nvec )
        END IF
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1600)
!$omp parallel do private(rr,d2,dd,g,gp,h,hp,eiqr,auxr) default(shared), reduction(+:aux), reduction(+:aux2), reduction(+:ene)
#endif
        DO nr = 1, nvec
           IF (llondon) THEN
              rr = r(:,nr)
              d2 = dist2(nr) * alat * alat
           ELSE
              rr = lvec(1,nr) * at(:,1) + lvec(2,nr) * at(:,2) + lvec(3,nr) * at(:,3) - dtau
              d2 = (rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3)) * alat * alat
           END IF
           dd = SQRT(d2)

           IF (dd <= eps8 .OR. (lxdm .AND. d2 > rmax2)) CYCLE
           IF (lxdm) THEN
              CALL calcgh_xdm(k,l,dd,g,gp,h,hp)
           ELSE
              CALL calcgh_d2(k,l,dd,g,gp,h,hp)
           END IF
           ene = ene - 0.5_dp * g

           eiqr = EXP(tpi * (0_dp,1_dp) * (q(1)*(rr(1)+dtau(1))+q(2)*(rr(2)+dtau(2))+q(3)*(rr(3)+dtau(3))))
           DO aa = 1 , 3
              DO bb = 1 , 3
                 IF (aa /= bb) THEN
                    auxr = hp * rr(aa) * alat * rr(bb) * alat / dd
                 ELSE
                    auxr = hp * rr(aa) * alat * rr(bb) * alat / dd + h
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
              der2disp(aa,k,bb,l) = aux2(aa,bb,l)
           ENDDO ! bb
        ENDDO ! aa
     ENDDO ! l

     DO l = 1, nat
        DO aa = 1, 3
           DO bb = 1, 3
              der2disp(aa,k,bb,k) = der2disp(aa,k,bb,k) - aux(aa,bb,l)
           ENDDO ! bb
        ENDDO ! aa
     ENDDO ! l
  ENDDO ! k
  
  CALL mp_sum(ene, intra_image_comm)
  CALL mp_sum(der2disp, intra_image_comm)
  
  IF (ionode) THEN
     WRITE (stdout,'(5X,A," energy = ",F17.8," Ry")') TRIM(namestr), ene
     WRITE (stdout,'(5X,"Done."/)')
  END IF

  ! cleanup
  IF (llondon) THEN
     CALL dealloca_london ()
  ENDIF

CONTAINS

  SUBROUTINE calcgh_xdm(i,j,d,g,gp,h,hp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    REAL(DP), INTENT(IN) :: d
    REAL(DP), INTENT(OUT) :: g, gp, h, hp

    REAL(DP) :: c6, c8, c10, rr
    REAL(DP) :: d2, d4, d6, d8, d10, dpr6, dpr8, dpr10, r2, r6, r8, r10
    REAL(DP) :: d5, d7, d9, dpr6sq, dpr8sq, dpr10sq, d17, d13, d3, dpr6cub
    REAL(DP) :: dpr8cub, dpr10cub

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
    gp = -10._dp * c10 * d9 / dpr10sq - 8._dp * c8 * d7 / dpr8sq - 6._dp * c6 * d5 / dpr6sq
    h = gp / d
    hp = -80._dp * c10 * d7 / dpr10sq + 200._dp * c10 * d17 / dpr10cub - 48._dp * c8 * d5 / dpr8sq &
       + 128._dp * c8 * d13 / dpr8cub - 24 * c6 * d3 / dpr6sq + 72._dp * c6 * d9 / dpr6cub

  END SUBROUTINE calcgh_xdm
  
  SUBROUTINE calcgh_d2(ii,jj,d,g,gp,h,hp)
    USE london_module, ONLY: beta, R_sum, C6_ij, scal6
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii, jj
    REAL(DP), INTENT(IN) :: d
    REAL(DP), INTENT(OUT) :: g, gp, h, hp

    INTEGER :: i, j
    REAL(DP) :: ed, fij, d6, d7, d2

    i = ityp(ii)
    j = ityp(jj)
    d2 = d * d
    d6 = d**6
    d7 = d6 * d
    ed = EXP(-beta * (d / R_sum(i,j) - 1._dp))
    fij = 1._dp / (1._dp + ed)
    g = C6_ij(i,j) * scal6 / d6 * fij
    gp = C6_ij(i,j) * scal6 / d6 / (1._dp + ed) * (beta * ed / R_sum(i,j) / (1._dp + ed) - 6._dp / d)
    h = gp / d
    hp = C6_ij(i,j) * scal6 / d7 / (1._dp + ed) * (48._dp / d2 - &
       13._dp * beta * ed / R_sum(i,j) / d / (1._dp + ed) - &
       beta**2 * ed / R_sum(i,j)**2 / (1._dp + ed)**2 * (1._dp - ed))

  END SUBROUTINE calcgh_d2

END SUBROUTINE d2ionq_disp

