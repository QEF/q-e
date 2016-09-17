!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM molecularpdos
  !-----------------------------------------------------------------------
  !
  ! Takes the projections onto orthogonalized atomic wavefunctions
  ! as computed by projwfc.x (see outdir/prefix.save/atomic_proj.xml)
  ! to build an LCAO-like representation of the eigenvalues of a system
  ! "full" and "part" of it (each should provide its own atomic_proj.xml file).
  ! Then the eigenvectors of the full system are projected onto the ones of the
  ! part.
  !
  ! An explanation of the keywords and the implementation is provided in
  ! Scientific Reports | 6:24603 | DOI: 10.1038/srep24603 (2016) (Supp. Info)
  !
  ! Typical application: decompose the PDOS of an adsorbed molecule into
  ! its molecular orbital, as determined by a gas-phase calculation.
  !
  ! The user has to specify which atomic functions (range beg:end) to use in
  ! both the full system and the part (the same atomic set should be used).
  !
  ! MOPDOS(E,ibnd_part) = \sum_k w_k [ \sum_{ibnd_full}
  !                                    <psi_{ibnd_part,k}|psi_{ibnd_full,k}>
  !                                    * \delta (E-\epsilon_{ibnd_full,k}) *
  !                                    <psi_{ibnd_full,k}|psi_{ibnd_part,k}> ]
  !
  ! where <psi_{ibnd_part,k}|psi_{ibnd_full,k}> are computed by using the LCAO
  ! representations:
  !
  ! |psi_{ibnd_full,k}> = 
  !        \sum_iatmwfc projs_full(iatmwfc,ibnd_full,k) |phi_{iatmwfc}>
  ! |psi_{ibnd_part,k}> = 
  !        \sum_iatmwfc projs_part(iatmwfc,ibnd_part,k) |phi_{iatmwfc}>
  !
  ! <psi_{ibnd_part,k}|psi_{ibnd_full,k}> =: projs_mo(ibnd_part,ibnd_full,k)
  !      = \sum_iatmwfc CONJG(projs_part(iatmwfc,ibnd_part,k))
  !                         * projs_full(iatmwfc,ibnd_full,k)
  !
  ! If kresolveddos=.true. from input, the summation over k is not performed
  ! and individual k-resolved contributions are given in output.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : PI, RYTOEV, eps4
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  USE mp_global,   ONLY : mp_startup
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  ! inputmopdos namelist
  CHARACTER(LEN=256) :: xmlfile_full, xmlfile_part, fileout
  INTEGER  :: i_atmwfc_beg_full, i_atmwfc_end_full, i_bnd_beg_full, i_bnd_end_full
  INTEGER  :: i_atmwfc_beg_part, i_atmwfc_end_part, i_bnd_beg_part, i_bnd_end_part 
  REAL(DP) :: degauss, Emax, Emin, DeltaE
  INTEGER  :: ngauss
  LOGICAL  :: kresolveddos
  !
  NAMELIST / inputmopdos / &
       xmlfile_full, i_atmwfc_beg_full, i_atmwfc_end_full, i_bnd_beg_full, i_bnd_end_full, &
       xmlfile_part, i_atmwfc_beg_part, i_atmwfc_end_part, i_bnd_beg_part, i_bnd_end_part, &
       fileout, Emin, Emax, DeltaE, ngauss, degauss, kresolveddos
  !
  INTEGER  :: nbnd_full, nkstot_full, num_k_points_full, nspin_full, natomwfc_full
  INTEGER  :: nbnd_part, nkstot_part, num_k_points_part, nspin_part, natomwfc_part
  REAL(DP), ALLOCATABLE :: xk_full(:,:), wk_full(:), et_full(:,:), &
       &                   xk_part(:,:), wk_part(:), et_part(:,:)
  !
  ! The read-from-file and the computed projections
  COMPLEX(DP), ALLOCATABLE :: projs_full(:,:,:), projs_part(:,:,:), projs_mo(:,:,:)  
  REAL(DP), ALLOCATABLE :: projs_mo_sq(:,:,:), mopdos(:,:,:,:), mopdostot(:,:,:)
  !
  ! For sorting projections
  INTEGER,  ALLOCATABLE :: idx(:)
  REAL(DP), ALLOCATABLE :: proj1 (:)
  !
  INTEGER  :: nkstot, natmwfc,num_k_points, nspin, ns, nwfc
  INTEGER  :: ibnd_full, ibnd_part, iatmwfc, ik, ik_eff, ik0, j, is, i, ios
  INTEGER  :: nksum, iksum
  REAL(DP) :: wksum
  !
  REAL(DP) :: Elw, Eup, delta, etev, psum
  INTEGER  :: ne, ie_delta, ie_mid, ie
  !
  REAL(DP), EXTERNAL :: w0gauss
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'MOLECULARPDOS' )
  !
  !    set default values for variables in namelist
  !
  i_atmwfc_beg_full=1
  i_atmwfc_end_full=0
  i_bnd_beg_full=1
  i_bnd_end_full=0
  i_atmwfc_beg_part=1
  i_atmwfc_end_part=0
  i_bnd_beg_part=1
  i_bnd_end_part=0
  fileout='molecularpdos'
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.01d0
  ngauss = 0
  degauss= 0.d0
  kresolveddos = .false.
  !
  ios = 0
  IF ( ionode ) THEN
     CALL input_from_file ( )
     READ (5, inputmopdos, iostat = ios)
  END IF
  !
  CALL mp_bcast (ios, ionode_id, world_comm )
  IF (ios /= 0) CALL errore ('molecularpdos', 'reading inputmopdos namelist', abs (ios) )
  !
  IF ( ionode ) THEN
     !
     ! This program is indeed serial
     !
     !
     ! Read projection files
     CALL readprojfile(xmlfile_full, nbnd_full, nkstot_full, num_k_points_full,&
              nspin_full, natomwfc_full, xk_full, wk_full, et_full, projs_full)
     CALL readprojfile(xmlfile_part, nbnd_part, nkstot_part, num_k_points_part,&
              nspin_part, natomwfc_part, xk_part, wk_part, et_part, projs_part)
     !
     ! Defaults ranges are maximum ones
     IF (i_atmwfc_end_full<1) i_atmwfc_end_full=natomwfc_full
     IF (i_atmwfc_end_part<1) i_atmwfc_end_part=natomwfc_part
     IF (i_bnd_end_full<1) i_bnd_end_full=nbnd_full
     IF (i_bnd_end_part<1) i_bnd_end_part=nbnd_part
     !
     ! Perform some of the several possible consistency checks.
     ! It is the user's responsibility to use the same setup for the
     ! SCF calculation of the two systems.
     IF ( nspin_full .GT. 2 ) &
          CALL errore ('molecularpdos','nspin not allowed',abs(nspin_full))
     IF ( nspin_full /= nspin_part ) &
          CALL errore ('molecularpdos','nspin does not match',abs(nspin_full-nspin_part))
     IF ( nkstot_full /= nkstot_part ) &
          CALL errore ('molecularpdos','nkstot does not match',abs(nkstot_full-nkstot_part))
     IF ( num_k_points_full /= num_k_points_part ) &
          CALL errore ('molecularpdos','num_k_points does not match',abs(num_k_points_full-num_k_points_part))
     IF ( i_atmwfc_end_part - i_atmwfc_beg_part /= i_atmwfc_end_full - i_atmwfc_beg_full ) &
          CALL errore ('molecularpdos','number of atomic wavefunctions does not match', &
          &            abs((i_atmwfc_end_part-i_atmwfc_beg_part)-(i_atmwfc_end_full-i_atmwfc_beg_full)) )
     !
     nkstot = nkstot_full
     num_k_points = num_k_points_full
     nspin = nspin_full 
     natmwfc=i_atmwfc_end_part-i_atmwfc_beg_part+1
     nbnd_part=i_bnd_end_part-i_bnd_beg_part+1
     ALLOCATE (projs_mo(i_bnd_beg_part:i_bnd_end_part, i_bnd_beg_full:i_bnd_end_full, nkstot))
     ALLOCATE (projs_mo_sq(i_bnd_beg_part:i_bnd_end_part, i_bnd_beg_full:i_bnd_end_full, nkstot))
     !
     WRITE( stdout,'(/5x,"Molecular orbitals used for projection")' )
     WRITE( stdout,'( 5x,"(data for the full system from file ",A,")")' ) TRIM(xmlfile_full)
     WRITE( stdout,'( 5x,"Atomic wavefunctions used: ",i5," - ",i5/)') i_atmwfc_beg_full, i_atmwfc_end_full
     !
     WRITE( stdout,'( 5x,"Projecting onto eigenvectors number: ",i5," - ",i5)') i_bnd_beg_part, i_bnd_end_part
     WRITE( stdout,'( 5x,"(of the subsytem described in file ",A,")")' ) TRIM(xmlfile_part)
     WRITE( stdout,'( 5x,"Atomic wavefunctions used: ",i5," - ",i5/)') i_atmwfc_beg_part, i_atmwfc_end_part
     !
     ! Compute the projection of full-system eigenvectors over part-system ones
     DO ik=1,nkstot
        DO ibnd_part=i_bnd_beg_part,i_bnd_end_part
           DO ibnd_full=i_bnd_beg_full,i_bnd_end_full
              !
              projs_mo(ibnd_part,ibnd_full,ik)=0
              DO iatmwfc=1, natmwfc
                 !
                 projs_mo(ibnd_part,ibnd_full,ik) = projs_mo(ibnd_part,ibnd_full,ik) &
                      + DCONJG(projs_part(iatmwfc+i_atmwfc_beg_part-1,ibnd_part,ik)) &
                      &      * projs_full(iatmwfc+i_atmwfc_beg_full-1,ibnd_full,ik)
                 !
              END DO
              !
              projs_mo_sq(ibnd_part,ibnd_full,ik) = &
                   DCONJG(projs_mo(ibnd_part,ibnd_full,ik)) &
                   &    * projs_mo(ibnd_part,ibnd_full,ik)
              !
           END DO
        END DO
     END DO
     !
     ! Write projections on standard output
     !
     ALLOCATE(idx(nbnd_part), proj1 (nbnd_part) )
     DO ik=1,nkstot
        !
        ! ik0 is the k-point index in [1:num_k_points]
        ik0 = MOD(ik,num_k_points)
        IF (ik0==0) ik0=num_k_points
        !
        WRITE( stdout, '(/" k = ",3f14.10)') (xk_full (i, ik0) , i = 1, 3)
        DO ibnd_full = i_bnd_beg_full, i_bnd_end_full
           WRITE( stdout, '("==== e(",i4,") = ",f11.5," eV ==== ")') &
                ibnd_full, et_full (ibnd_full, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           DO nwfc = 1, nbnd_part
              idx (nwfc) = 0
              proj1 (nwfc) = - projs_mo_sq (nwfc+i_bnd_beg_part-1, ibnd_full, ik)
           ENDDO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (nbnd_part, proj1, idx, eps4)
           !
           !  only projections that are larger than 0.001 are written
           !
           DO nwfc = 1, nbnd_part
              proj1 (nwfc) = - proj1(nwfc)
              IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
           ENDDO
           nwfc = nbnd_part + 1
20         nwfc = nwfc -1
           !
           idx(:) = idx(:) + i_bnd_beg_part - 1
           !
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') &
                (proj1 (i), idx(i), i = 1, min(5,nwfc))
           DO j = 1, (nwfc-1)/5
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') &
                   (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
           ENDDO
           psum = SUM ( projs_mo_sq(i_bnd_beg_part:i_bnd_end_part, ibnd_full, ik) )
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        END DO
     END DO
     DEALLOCATE(idx, proj1)
     !
     ! Prepare to plot the molecular orbital projected DOS
     !
     ! Find band extrema
     Elw = et_full (1, 1)
     Eup = et_full (nbnd_full, 1)
     DO ik = 2, nkstot_full
        Elw = min (Elw, et_full (1, ik) )
        Eup = max (Eup, et_full (nbnd_full, ik) )
     ENDDO
     IF (degauss/=0.d0) THEN
        Eup = Eup + 3d0 * degauss
        Elw = Elw - 3d0 * degauss
     ENDIF
     Emin = max (Emin/rytoev, Elw)
     Emax = min (Emax/rytoev, Eup)
     DeltaE = DeltaE/rytoev
     ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
     ie_delta = 5 * degauss / DeltaE + 1
     !
     IF (kresolveddos) THEN
        nksum=num_k_points
     ELSE
        nksum=1
     ENDIF
     !
     ALLOCATE (mopdos(i_bnd_beg_part:i_bnd_end_part, 0:ne,nspin,nksum))
     mopdos(:,:,:,:)=0d0
     !
     ! Compute mopdos(E)
     DO ik=1,num_k_points
        IF (kresolveddos) THEN
           ! set equal weight to all k-points
           wksum=1.D0
           ! do not sum over k-points
           iksum=ik
        ELSE
           ! use true weights
           wksum=wk_full(ik)
           ! contributions from all k-points are summed in mopdos(:,:,:,iksum)
           iksum=1
        ENDIF
        !
        DO ibnd_full=i_bnd_beg_full,i_bnd_end_full
           !
           etev = et_full(ibnd_full,ik)
           ie_mid = nint( (etev-Emin)/DeltaE )
           !
           DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
              delta=w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                   / degauss / rytoev
              DO ibnd_part=i_bnd_beg_part,i_bnd_end_part
                 mopdos(ibnd_part,ie,1,iksum) = mopdos(ibnd_part,ie,1,iksum)  &
                      +  wksum * delta * projs_mo_sq(ibnd_part,ibnd_full,ik)
              END DO
           END DO
           !
           IF ( nspin == 2 ) THEN
              !
              ik_eff = ik + num_k_points
              etev = et_full(ibnd_full,ik_eff)
              ie_mid = nint( (etev-Emin)/DeltaE )
              !
              DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
                 delta=w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                      / degauss / rytoev
                 DO ibnd_part=i_bnd_beg_part,i_bnd_end_part
                    mopdos(ibnd_part,ie,2,iksum) = mopdos(ibnd_part,ie,2,iksum)  &
                         +  wksum * delta * projs_mo_sq(ibnd_part,ibnd_full,ik_eff)
                 END DO
              END DO
              !
           ENDIF
           !
        END DO
     END DO
     !
     ! Ouput mopdos(E)
     OPEN(UNIT=12, FILE=TRIM(fileout)//".mopdos", ACTION="write", STATUS="replace")
     !
     IF (kresolveddos) THEN
        WRITE (12,'("# ik   ")', advance="NO")
     ENDIF
     !
     IF ( nspin == 2 ) THEN
        WRITE (12,'("# ibnd_part  E (eV)  tot_up(E)  tot_dw(E)")')
     ELSE
        WRITE (12,'("# ibnd_part  E (eV)  tot(E)")')
     END IF
     !
     DO iksum=1, nksum
        !
        DO ibnd_part=i_bnd_beg_part,i_bnd_end_part
           !
           DO ie=0,ne
              etev = Emin + ie * DeltaE
              IF (kresolveddos) THEN
                 WRITE (12,'(i5," ")', advance="NO") iksum
              ENDIF
              WRITE(12,'(i11," ")', advance="NO") ibnd_part
              WRITE(12,'(f7.3)', advance="NO") etev*rytoev
              DO is=1,nspin
                 WRITE(12,'(e11.3)', advance="NO") mopdos(ibnd_part,ie,is,iksum)
              END DO
              WRITE(12,'()')
           END DO
           WRITE (12,'()')
        END DO
        IF (kresolveddos) WRITE (12,'()')
     END DO
     CLOSE(12)
     !
     ! Compute the total mopdos(E)
     ALLOCATE (mopdostot(0:ne,nspin,nksum))
     mopdostot(:,:,:)=0d0
     DO iksum=1, nksum
        DO ns=1, nspin
           DO ie=0,ne
              DO ibnd_part=i_bnd_beg_part,i_bnd_end_part
                 mopdostot(ie,ns,iksum)=mopdostot(ie,ns,iksum) + mopdos(ibnd_part,ie,ns,iksum)
              END DO
           END DO
        END DO
     END DO
     !
     ! Ouput the total mopdos(E)
     OPEN(UNIT=13, FILE=TRIM(fileout)//".mopdos_tot", ACTION="write", STATUS="replace")
     !
     IF (kresolveddos) THEN
        WRITE (13,'("# ik   ")', advance="NO")
     ENDIF
     !
     IF ( nspin == 2 ) THEN
        WRITE (13,'("# E (eV)  tot_up(E)  tot_dw(E) ")')
     ELSE
        WRITE (13,'("# E (eV)  tot(E) ")')
     END IF
     DO iksum=1, nksum
        DO ie=0,ne
           etev = Emin + ie * DeltaE
           IF (kresolveddos) THEN
              WRITE (13,'(i5," ")', advance="NO") iksum
           ENDIF
           WRITE(13,'(f7.3)', advance="NO") etev*rytoev
           DO is=1,nspin
              WRITE(13,'(e11.3)', advance="NO") mopdostot(ie,is,iksum)
           END DO
           WRITE(13,'()')
        END DO
        IF (kresolveddos) WRITE(13,'()')
     END DO
     CLOSE(13)
     !
     DEALLOCATE (xk_full, wk_full, et_full, projs_full)
     DEALLOCATE (xk_part, wk_part, et_part, projs_part)
     DEALLOCATE (projs_mo, projs_mo_sq, mopdos, mopdostot)
     !
  END IF
  !
  CALL environment_end ( 'MOLECULARPDOS' ) 
  !
  CALL stop_pp
  !
CONTAINS
  !-----------------------------------------------------------------------  
  SUBROUTINE readprojfile(filexml,nbnd,nkstot,num_k_points,nspin,natomwfc,xk,wk,et,projs)
    !-----------------------------------------------------------------------
    USE iotk_module
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: filexml
    INTEGER  :: nbnd, nkstot, num_k_points, nspin, natomwfc
    REAL(DP),    ALLOCATABLE :: xk(:,:), wk(:), et(:,:)
    COMPLEX(DP), ALLOCATABLE :: projs(:,:,:)
    !
    REAL(DP)  :: ef, nelec
    INTEGER :: ik, ik_eff, ia, ierr, iun
    LOGICAL :: noncolin
    !
    CALL iotk_open_read(iun, FILE=TRIM(filexml), &
         BINARY=.FALSE., IERR=ierr )
    IF ( ierr /= 0 ) STOP 'error reading file'
    CALL iotk_scan_begin(iun, "HEADER")
    CALL iotk_scan_dat(iun, "NUMBER_OF_BANDS", nbnd)
    CALL iotk_scan_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
    CALL iotk_scan_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
    CALL iotk_scan_dat(iun, "NON-COLINEAR_CALCULATION", noncolin )
    IF (noncolin) STOP 'noncolin not implemented in molecularpdos'
    CALL iotk_scan_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
    CALL iotk_scan_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
    CALL iotk_scan_dat(iun, "FERMI_ENERGY", ef )
    CALL iotk_scan_end(iun, "HEADER")
    !
    nkstot = num_k_points
    IF ( nspin == 2 ) nkstot = nkstot * 2
    !
    ALLOCATE (xk(3,num_k_points))
    ALLOCATE (wk(num_k_points))
    ALLOCATE (et(nbnd,nspin*nkstot))
    ALLOCATE (projs(natomwfc,nbnd,nkstot))
    !
    CALL iotk_scan_dat(iun, "K-POINTS", xk(:,1:num_k_points) )
    CALL iotk_scan_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points) )
    !
    CALL iotk_scan_begin(iun, "EIGENVALUES")
    !
    DO ik=1, num_k_points
       CALL iotk_scan_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
       IF ( nspin == 2 ) THEN
          !
          ik_eff = ik + num_k_points
          !
          CALL iotk_scan_dat( iun, "EIG.1", et(:,ik) )
          CALL iotk_scan_dat( iun, "EIG.2", et(:,ik_eff) )
          !
       ELSE
          !
          CALL iotk_scan_dat( iun, "EIG", et(:,ik) )
          !
       ENDIF
       !
       CALL iotk_scan_end( iun, "K-POINT"//trim(iotk_index(ik)) )
       !
    ENDDO
    !
    CALL iotk_scan_end(iun, "EIGENVALUES")
    !
    CALL iotk_scan_begin(iun, "PROJECTIONS")
    !
    DO ik=1,num_k_points
       !
       CALL iotk_scan_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
       !
       IF ( nspin == 2 ) THEN
          !
          CALL iotk_scan_begin ( iun, "SPIN.1" )
          !
          DO ia = 1, natomwfc
             CALL iotk_scan_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
          ENDDO
          !
          CALL iotk_scan_end ( iun, "SPIN.1" )
          !
          ik_eff = ik + num_k_points
          !
          CALL iotk_scan_begin ( iun, "SPIN.2" )
          !
          DO ia = 1, natomwfc
             CALL iotk_scan_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
          ENDDO
          !
          CALL iotk_scan_end ( iun, "SPIN.2" )
          !
       ELSE
          !
          DO ia = 1,natomwfc
             CALL iotk_scan_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
          ENDDO
          !
       ENDIF
       !
       CALL iotk_scan_end( iun, "K-POINT"//trim(iotk_index(ik)) )
       !
    ENDDO
    !
    CALL iotk_scan_end(iun, "PROJECTIONS")
    !
    CALL iotk_close_read(iun)
    !
  END SUBROUTINE readprojfile

END PROGRAM molecularpdos
