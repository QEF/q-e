!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Program "fermisurface", at the end of module fs, writes the Fermi Surface
! in a format suitable for plotting by XCrySDen, from the xml data file of a
! non-scf calculation done with a uniform grid of k-points. The complete grid,
! including corners, is reconstructed from the symmetry-reduced PWscf one.
!
! The pw.x run should be a 'scf' or 'nscf' calculation, with smearing or
! tetrahedra, and a sufficiently dense automatic unshifted grid as follows:
! K_POINTS automatic
! nk1 nk2 nk3 0 0 0
!
! Written by Paolo Giannozzi, based on the "bands_FS" utility by Eyvaz Isaev
!
MODULE fs
  !
  ! Input variables. On input, namelist &fermi ... / contains:
  ! outdir : as given in the preceding pw.x run (default as in pw.x)
  ! prefix : as given in the preceding pw.x run (default as in pw.x)
  ! file_fs: filename for XCrySDen output (default: same as "prefix")
  ! deltaE : include all bands crossing Ef within +/- deltaE (eV, default: 1 eV)
  !
  USE kinds, ONLY: dp
  CHARACTER(LEN=80) :: file_fs
  REAL(dp) :: deltaE
  !
  ! Internal variables
  !
  ! Number of k-points in the complete grid, including corners
  ! nkfs = (nk1+1)*(nk2+1)*(nk3+1), with nk1,nk2,nk3 as in the pw.x input
  INTEGER :: nkfs
  ! Index, in the k-point list, of the k-points of the complete grid 
  INTEGER, ALLOCATABLE:: equivalent_kpoint(:)
  !
  PRIVATE
  !
  ! Subroutines
  !
  PUBLIC :: read_input_fs, fill_fs_grid, write_xcrysden_fs
  !
CONTAINS
  !--------------------------------------------------------------------
  SUBROUTINE read_input_fs
    !--------------------------------------------------------------------
    !
    USE io_global,  ONLY : ionode, ionode_id
    USE io_files,   ONLY : prefix, tmp_dir
    USE mp,         ONLY : mp_bcast
    USE mp_world,   ONLY : world_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), EXTERNAL :: trimcheck
    CHARACTER(len=256) :: outdir
    INTEGER :: ios
    !
    NAMELIST /fermi/ outdir, prefix, file_fs, deltaE
    !
    ios = 0
    IF ( ionode ) THEN
       !
       !   set default values for variables in namelist
       !
       CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
       IF ( trim( outdir ) == ' ' ) outdir = './'
       prefix ='pwscf'
       file_fs=' '
       deltaE = 1.0_dp
       !
       CALL input_from_file ( )
       !
       READ (5, fermi, iostat=ios )
       !
       tmp_dir = trimcheck (outdir)
       IF ( file_fs == ' ' ) file_fs = trim(prefix)//'_fs'
       !
    ENDIF
    !
    CALL mp_bcast( ios, ionode_id, world_comm )
    IF ( ios /= 0 ) CALL errore('fermi','reading fermi namelist',abs(ios))
    !
    ! ... Broadcast variables
    !
    CALL mp_bcast( file_fs,ionode_id, world_comm )
    CALL mp_bcast( tmp_dir,ionode_id, world_comm )
    CALL mp_bcast( prefix, ionode_id, world_comm )
    CALL mp_bcast( deltaE, ionode_id, world_comm )
    !
  !
  END SUBROUTINE read_input_fs
  !
  !-----------------------------------------------------------------------
  SUBROUTINE fill_fs_grid
    !-----------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : nsym, s, time_reversal, t_rev
    USE cell_base, ONLY : at
    USE start_k,   ONLY : k1,k2,k3, nk1,nk2,nk3
    USE klist,     ONLY : xk, nks, nkstot
    USE lsda_mod,  ONLY : nspin
    !
    IMPLICIT NONE
    !
    REAL(DP) :: xkr(3), deltap(3), deltam(3)
    REAL(DP), PARAMETER:: eps=1.0d-5
    REAL(DP), ALLOCATABLE :: xkg(:,:)
    INTEGER :: nks_
    INTEGER :: i,j,k, ns, n, nk, ierr
    !
    ! nks_ is the true number of k-points (they are doubled in LSDA)
    !
    IF ( nspin == 2 ) THEN
       nks_ = nkstot/2
    ELSE
       nks_ = nkstot
    END IF
    !
    ! Generate a uniform grid of nkfs k-points, xkg, including corners
    !
    IF ( nk1 == 0 .OR. nk2 == 0 .OR. nk3 == 0 .OR. &
          k1 == 1 .OR.  k2 == 1 .OR.  k3 == 1 )    &
       CALL errore ('fill_fs_grid', 'uniform unshifted k-point grid expected',1)
    !
    nkfs = (nk1+1)*(nk2+1)*(nk3+1)
    ALLOCATE (equivalent_kpoint( nkfs))
    ALLOCATE (xkg( 3,nkfs))
    !
    DO i=1,nk1+1
       DO j=1,nk2+1
          DO k=1,nk3+1
             !  this is nothing but consecutive ordering
             n = (k-1) + (j-1)*(nk3+1) + (i-1)*(nk2+1)*(nk3+1) + 1
             !  xkg are the components of the complete grid in crystal axis
             xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
             xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
             xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
          ENDDO
       ENDDO
    ENDDO
    !
    !  locate k-points of the uniform grid in the list of irreducible k-points
    !
    !  bring irreducible k-points to crystal axis (beware: overwrites xk)
    ! (but there is no need to bring k-points back to cartesian axis)
    !
    CALL cryst_to_cart (nks,xk,at,-1)
    !
    ierr = 0
    DO nk=1,nkfs
       DO n=1,nks_
          DO ns=1,nsym
             DO i=1,3
                xkr(i) = s(i,1,ns) * xk(1,n) + &
                     s(i,2,ns) * xk(2,n) + &
                     s(i,3,ns) * xk(3,n)
             ENDDO
             IF(t_rev(ns)==1) xkr = -xkr
             !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
             DO i=1,3
                deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
                deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
             ENDDO
             !  deltap is the difference vector, brought back in the first BZ
             !  deltam is the same but with k => -k (for time reversal)
             IF ( sqrt ( deltap(1)**2 + &
                         deltap(2)**2 + &
                         deltap(3)**2 ) < eps .or. ( time_reversal .and. &
                  sqrt ( deltam(1)**2 +  &
                         deltam(2)**2 +  &
                         deltam(3)**2 ) < eps ) ) THEN
                !  equivalent irreducible k-point found
                equivalent_kpoint(nk) = n
                GOTO 15
             ENDIF
          ENDDO
       ENDDO
       !  equivalent irreducible k-point not found - something wrong
       CALL errore('fill_fs_grid','cannot locate  k point',nk)
15     CONTINUE
    ENDDO
    
    DEALLOCATE(xkg)
 
    DO n=1,nks_
       DO nk=1,nkfs
          IF (equivalent_kpoint(nk)==n) GOTO 20
       ENDDO
       !  this failure of the algorithm may indicate that the displaced grid
       !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
       CALL errore('fill_fs_grid','cannot remap grid on k-point list',n)
20     CONTINUE
    ENDDO
    
  END SUBROUTINE fill_fs_grid
  !
  !----------------------------------------------------------------------- 
  SUBROUTINE write_xcrysden_fs
    !----------------------------------------------------------------------- 
    !
    USE constants, ONLY : rytoev
    USE cell_base, ONLY : bg
    USE klist,     ONLY : nkstot
    USE lsda_mod,  ONLY : nspin
    USE ener,      ONLY : Ef 
    USE wvfct,     ONLY : nbnd, et
    USE start_k,   ONLY : nk1, nk2, nk3
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE :: e_fs(:,:)
    REAL(dp) :: emin, emax, efermi
    INTEGER :: first_band, last_band, i, ik, nkb2
    CHARACTER(LEN=256) :: filename
    !
    efermi = Ef*rytoev ! Ef is in a.u., convert to eV
    !
    ! unpolarized or noncolinear calculations
    !
    IF ( nspin /= 2 ) THEN
       ! find bands that cross the Fermi surface
       first_band = 0
       last_band  = 0
       DO i=1,nbnd
          emin = MINVAL(et(i,1:nkstot))*rytoev
          emax = MAXVAL(et(i,1:nkstot))*rytoev
          IF ( emin-deltaE < efermi .AND. emax+deltaE > efermi ) THEN
             ! include all bands crossing the Fermi energy 
             ! within a tolerance +/- deltaE
             IF ( first_band == 0 ) first_band = i
             last_band = i
          END IF
       END DO
       WRITE(stdout,'(5X,i4," bands found crossing Ef =",f12.6)') &
             last_band-first_band+1, efermi
       ALLOCATE (e_fs(last_band-first_band+1,nkfs))
       !
       DO ik=1,nkfs
          DO i=1,last_band-first_band+1
             e_fs(i,ik) = rytoev * et(i+first_band-1,equivalent_kpoint(ik))
          END DO
       END DO
       filename = TRIM(file_fs)//'.bxsf'
       CALL xsf_fs ( filename, Ef*rytoev, first_band, last_band, &
                     nk1,nk2,nk3, bg(:,1),bg(:,2),bg(:,3), e_fs )
    ELSE
       !
       ! spin up
       !
       nkb2=nkstot/2
       ! find bands that cross the Fermi surface
       first_band = 0
       last_band  = 0
       DO i=1,nbnd
          emin = MINVAL(et(i,1:nkb2))*rytoev
          emax = MAXVAL(et(i,1:nkb2))*rytoev
          IF ( emin-deltaE < efermi .AND. emax+deltaE > efermi ) THEN
             ! include all bands crossing the Fermi energy 
             ! within a tolerance +/- deltaE
             IF ( first_band == 0 ) first_band = i
             last_band = i
          END IF
       END DO
       WRITE(stdout,'(5X,i4," spin-up bands found crossing Ef =",f12.6)') &
             last_band-first_band+1, efermi
       ALLOCATE (e_fs(last_band-first_band+1,nkfs))
       DO ik=1,nkfs
          DO i=1,last_band-first_band+1
             e_fs(i,ik) = rytoev * et(i+first_band-1,equivalent_kpoint(ik))
          END DO
       END DO
       filename = TRIM(file_fs)//'up.bxsf'
       CALL xsf_fs ( filename, Ef*rytoev, first_band, last_band, &
                     nk1,nk2,nk3, bg(:,1),bg(:,2),bg(:,3), e_fs )
       DEALLOCATE (e_fs)
       !
       ! spin down (second half of the k-points)
       !
       ! find bands that cross the Fermi surface
       first_band = 0
       last_band  = 0
       DO i=1,nbnd
          emin = MINVAL(et(i,nkb2+1:nkstot))*rytoev
          emax = MAXVAL(et(i,nkb2+1:nkstot))*rytoev
          IF ( emin-deltaE < efermi .AND. emax+deltaE > efermi ) THEN
             ! include all bands crossing the Fermi energy 
             ! within a tolerance +/- deltaE
             IF ( first_band == 0 ) first_band = i
             last_band = i
          END IF
       END DO
       WRITE(stdout,'(5X,i4," spin-down bands found crossing Ef =",f12.6)') &
             last_band-first_band+1, efermi
       ALLOCATE (e_fs(last_band-first_band+1,nkfs))
       DO ik=1,nkfs
          DO i=1,last_band-first_band+1
             e_fs(i,ik) = rytoev * et(i+first_band-1,nkb2+equivalent_kpoint(ik))
          END DO
       END DO
       filename = TRIM(file_fs)//'dw.bxsf'
       CALL xsf_fs ( filename, Ef*rytoev, first_band, last_band, &
                     nk1,nk2,nk3, bg(:,1),bg(:,2),bg(:,3), e_fs )
    END IF
    DEALLOCATE (e_fs)
    DEALLOCATE (equivalent_kpoint)
      
  END SUBROUTINE write_xcrysden_fs
  !
  !----------------------------------------------------------------------- 
  SUBROUTINE xsf_fs (filename, Ef, n_start, n_last, na,nb,nc, &
       x, y, z, e_fs )
    !----------------------------------------------------------------------- 
    !     Write  header file here
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n_start, n_last, na, nb, nc
    REAL(dp), INTENT(IN) :: Ef, x(:), y(:), z(:), e_fs(:,:)
    REAL(dp) :: x0=0.0_dp, y0=0.0_dp, z0=0.0_dp
    INTEGER :: i, ik, fsunit = 4
    !
    OPEN (unit = fsunit, file = filename, status='unknown', form='formatted')
    WRITE(fsunit, '(" BEGIN_INFO")')
    WRITE(fsunit, '("   #")') 
    WRITE(fsunit, '("   # this is a Band-XCRYSDEN-Structure-File")')
    WRITE(fsunit, '("   # aimed at Visualization of Fermi Surface")')
    WRITE(fsunit, '("   #")') 
    WRITE(fsunit, '("   # Case:   ",A)')     TRIM(filename)
    WRITE(fsunit, '("   #")') 
    WRITE(fsunit, '(" Fermi Energy:    ", f12.4)') Ef
    WRITE(fsunit, '(" END_INFO")') 
     
    WRITE(fsunit, '(" BEGIN_BLOCK_BANDGRID_3D")')
    WRITE(fsunit, '(" band_energies")')
    WRITE(fsunit, '(" BANDGRID_3D_BANDS")')
    WRITE(fsunit, '(I5)')  n_last-n_start+1
    WRITE(fsunit, '(3I5)') na+1, nb+1, nc+1
    WRITE(fsunit, '(3f10.6)') x0,  y0,  z0
    WRITE(fsunit, '(3f10.6)') x(1), x(2), x(3)
    WRITE(fsunit, '(3f10.6)') y(1), y(2), y(3)
    WRITE(fsunit, '(3f10.6)') z(1), z(2), z(3)
    DO i=n_start, n_last
       WRITE(fsunit, '("BAND:", i4)') i
       WRITE(fsunit, '(6f10.4)') (e_fs(i-n_start+1,ik),ik=1,nkfs)
    END DO
    WRITE(fsunit, '(" END_BANDGRID_3D")')
    WRITE(fsunit, '(" END_BLOCK_BANDGRID_3D")')
    !
    CLOSE  (unit = fsunit, status='keep' )
    !
  END SUBROUTINE xsf_fs

END MODULE fs
!
!--------------------------------------------------------------------
PROGRAM fermisurface
  !--------------------------------------------------------------------
  !
  USE fs
  USE io_global,  ONLY : ionode
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start, environment_end
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL environment_start ( 'FERMI' )
  !
  CALL read_input_fs ( )
  CALL read_xml_file ( )
  CALL fill_fs_grid ( )
  IF ( ionode ) CALL write_xcrysden_fs ( )
  !
  CALL environment_end ( 'FERMI' )
  !
  CALL stop_pp
  !
END PROGRAM fermisurface
!
