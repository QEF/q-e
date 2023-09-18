!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE d2ionq_dispd3( alat, nat, at, q, der2disp )
  !------------------------------------------------------------------------
  USE kinds,         ONLY: DP
  USE io_global,     ONLY: ionode, ionode_id, stdout
  USE control_ph,    ONLY: dftd3_hess
  USE constants,     ONLY: tpi
  USE control_lr,    ONLY: lgamma
  USE dftd3_qe,      ONLY: print_dftd3_hessian
  USE mp_images,     ONLY: intra_image_comm
  USE mp,            ONLY: mp_bcast

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: alat
  !! cell parameter (celldm(1))
  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the unit cell
  REAL(DP), INTENT(IN) :: at(3,3)
  !! at(:,i) is lattice vector i in alat units
  REAL(DP), INTENT(IN) :: q(3)
  !! wavevector in 2pi/alat units
  COMPLEX(DP), INTENT(INOUT) :: der2disp(3,nat,3,nat)
  !! dispersion contribution to the (massless) dynamical matrix
  ! 
  INTEGER :: n, nn, nnn, rep(3), nrep, nhess, irep, jrep, krep, irp, jrp, krp
  INTEGER :: i, j, iat, jat, ixyz, jxyz
  INTEGER :: iprint
  CHARACTER(LEN=100) :: string
  REAL(DP), ALLOCATABLE :: d3hess(:,:,:,:,:,:,:), buffer(:)
  COMPLEX(DP), ALLOCATABLE :: mmat(:,:,:,:)
  COMPLEX(DP) :: eiqr, tt(3)
  LOGICAL :: q_gamma ! whether the Hessian stored in the file has been computed for q=0,0,0 only 
  ! 
  if( ionode ) then 
    !
    der2disp = (0._dp, 0._dp)
    !
    write(stdout,'(/,5x,2A)') 'Reading Grimme-D3 Hessian from file: ', TRIM(dftd3_hess)
    !
    ! Reading Hessian from file
    !
    OPEN (unit = 1, file = dftd3_hess, status = 'unknown')
    READ(1, * ) 
    READ(1, * ) string, rep(1:3), n, nn, q_gamma
    !
    ! some consistency checks before allocation
    IF((q_gamma.eqv..true.) .and. (lgamma.eqv..false.)) THEN
      Call errore('d2ionq_dispd3', 'The Hessian in the file is only good for q=0,0,0. Recompute it with q_gamma=.false.', 1)
    ELSE 
      WRITE( stdout, '(/,5x,A,3I4)') 'Number of cells replicated along each semiaxis: ', rep(1), rep(2), rep(3)
    END IF
    !
    iprint = 1
    if(q_gamma) iprint = 0
    !
    nrep = (2*rep(1)+1) * (2*rep(2)+1) * (2*rep(3)+1)  ! number of unit cells in the supercell 
    nnn = nat * nrep                                   ! number of atoms in the supercell 
    nhess = (3 * nat)**2 * nrep                        ! Hessian dimensions
    IF((n.ne.nat).or.(nn.ne.nnn)) THEN
      WRITE(stdout, '(/,5x,4I9)' ) nat, n, nn, nnn
      Call errore('d2ionq_dispd3', 'Wrong cell or supercell size', 1)
    END IF

    WRITE( stdout, '(5x,A,I9)') 'Number of cells in the supercell: ', nrep 
    WRITE( stdout, '(5x,A,I9)') 'Number of atoms in the supercell: ', nn  
    WRITE( stdout, '(5x,A,I9)') 'Hessian allocation dimension: ', nhess 

    ALLOCATE( d3hess(-rep(3):rep(3),-rep(2):rep(2),-rep(1):rep(1), 3,nat,3,nat), buffer(3*nat), mmat(3,nat,3,nat) )
    IF( size(d3hess) .ne. nhess ) Call errore('d2ionq_dispd3', "Wrong Hessian dimensions", 1)
    d3hess(:,:,:,:,:,:,:)=0.0_dp
    !
    DO irep = -rep(1), rep(1)
      DO jrep = -rep(2), rep(2)
        DO krep = -rep(3), rep(3)
          !
          READ(1, * ) string,string, string,irp, string,jrp, string,krp
          IF(irep.ne.irp .or. jrep.ne.jrp .or. krep.ne.krp ) Call errore('d2ionq_dispd3', "Wrong Hessian I/O", 1)
          !
          DO i = 1, 3*nat         ! 1 2 3 4 5 6 7 8 9 ... 3*nat
            iat  = (i+2)/3        ! 1 1 1 2 2 2 3 3 3 ... nat 
            ixyz = i - 3* (iat-1) ! 1 2 3 1 2 3 1 2 3 ... 3 
            READ(1, * ) buffer(1:3*nat) 
            !
            DO j = 1, 3*nat
              jat  = (j+2)/3 
              jxyz = j - 3* (jat-1) 
              d3hess(krep,jrep,irep,ixyz,iat,jxyz,jat) = buffer(j)
            END DO 
            !  
          END DO 
          !
        END DO 
      END DO 
    END DO 
    !
    CLOSE (1)
    !
    DEALLOCATE( buffer) 
    !
    write(stdout,'(/,5x,A)') 'Grimme-D3 Hessian read '
    !
    ! Computing dynamical matrix 
    !
    WRITE(stdout,'(/,5x,A,3f12.6)') 'Computing dynamical matrix for q: ', q(1:3) 
    !
    mmat = (0.0_dp, 0.0_dp)
    DO krep = -rep(3), rep(3)
      DO jrep = -rep(2), rep(2)
        DO irep = -rep(1), rep(1)
          !
          tt(:) = alat * ( irep * at(:,1) + jrep * at(:,2) + krep * at(:,3) )
          eiqr = EXP(- tpi * (0_dp,1_dp) * ( q(1)*tt(1)+q(2)*tt(2)+q(3)*tt(3) ) )
          DO ixyz = 1, 3
            DO iat = 1, nat
              DO jxyz = 1, 3
                DO jat = 1, nat
                  der2disp(ixyz,iat,jxyz,jat) = der2disp(ixyz,iat,jxyz,jat) &
                            + d3hess(krep,jrep,irep,ixyz,iat,jxyz,jat) * eiqr
 
                  mmat(ixyz,iat,jxyz,jat) = mmat(ixyz,iat,jxyz,jat) &
                            + d3hess(krep,jrep,irep,ixyz,iat,jxyz,jat) * eiqr
                END DO 
              END DO 
            END DO 
          END DO 
          !
        END DO 
      END DO 
    END DO 
    !
    if(q_gamma) then 
      write(string,'(A)' ) 'true'
    else
      write(string,'(A)' ) 'false'
    end if 
    CALL print_dftd3_hessian( mmat, n, string )
    !
    DEALLOCATE( d3hess, mmat ) 
    !
    WRITE(stdout,'(/,5x,A)') 'Dynamical matrix computed'
    !
  end if 
  CALL mp_bcast(der2disp, ionode_id, intra_image_comm)
  !
  RETURN
  !
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

