!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM pw2gt
  !-----------------------------------------------------------------------
  ! Convert QE data files to the file format used by the GreenT code.
  ! The output is a single, simple, formatted and self-contained file.
  ! Can be useful for anybody wanting to access and process QE data
  ! with little hassle. 
  !
  ! IMPORTANT NOTICE 1: must be run on a single processor ONLY
  ! IMPORTANT NOTICE 2: the file may easily become VERY large!
  !
  ! Input: namelist &inputpp ... /
  !        outdir   as in QE input (default values as in QE).
  !        prefix     "     "       "        "       "
  !        check    .false. (default) / .true.
  !                 if .true. do not write the file but read it,
  !                 fill and diagonalize the hamiltonian matrix
  !                 BEWARE! do not attempt to do that for systems 
  !                         exceeding a few thousands plane waves
  !
  ! Output file written in "outdir"/"prefix".save/output.dat
  ! NOTA BENE: complex numbers are written as "a b", not "(a,b)"
  !
  ! Written by Paolo Giannozzi, with help from Marco Pala
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : tmp_dir, prefix, restart_dir
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE mp_images,  ONLY : intra_image_comm, nproc_image
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: needwf = .true., check
  INTEGER :: ios
  CHARACTER(LEN=256) :: outdir, fileout
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / inputpp / outdir, prefix, check
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL environment_start ( 'PW2GT' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  check = .false.
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm)
  IF ( ios /= 0) CALL errore ('pw2gt', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  !
  IF ( nproc_image > 1 ) CALL errore ('pw2gt', &
       'must be run on a single processor only', nproc_image )
  !
  fileout = TRIM ( restart_dir() ) //  'output.dat'
  IF ( check ) THEN
     !
     CALL simple_diag ( fileout )
     !
  ELSE
     !
     !   Read xml file, allocate and initialize general variables
     !
     CALL read_file_new ( needwf )
     CALL simple_output ( fileout )
     !
  END IF
  !  
  CALL environment_end ( 'PW2GT' )
  !
  CALL mp_global_end ()
  !
END PROGRAM pw2gt
!----------------------------------------------------------------------------
SUBROUTINE simple_output ( fileout  )
  !----------------------------------------------------------------------------
  !
  ! Not-so-smart but easy-to-read output file for simple cases (e.g. Si)
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : restart_dir
  USE cell_base,    ONLY : at, bg, alat, tpiba
  USE ions_base,    ONLy : nat, nsp, ityp, atm, tau
  USE gvect,        ONLY : ngm, mill
  USE gvecs,        ONLY : doublegrid
  USE scf,          ONLY : vrs, vltot, v, kedtau
  USE fft_rho,      ONLY : rho_r2g
  USE lsda_mod,     ONLY : nspin, isk
  USE klist,        ONLY : nks, xk, ngk, igk_k
  USE uspp,         ONLY : okvan, nkb, vkb, dvan, dvan_so
  USE uspp_param,   ONLY : nh
  USE uspp_init,    ONLY : init_us_2
  USE fft_base,     ONLY : dfftp, dffts
  USE wvfct,        ONLY : nbnd, et, npwx
  USE wavefunctions,     ONLY : evc
  USE noncollin_module,  ONLY : lspinorb, domag
  USE pw_restart_new,    ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), intent(in) :: fileout
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:)
  INTEGER :: iun=4, ig, is, ik, ikb, ibnd, na, nt, npw
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing simple output data file ",A)' ) &
       trim(fileout)
  OPEN ( UNIT = iun, FORM = 'formatted', STATUS = 'unknown', &
       FILE = fileout )
  WRITE(iun,'("# Primitive lattice vectors a_1, a_2, a_3 (a.u.)")')
  WRITE(iun,*) alat*at(:,1), alat*at(:,2), alat*at(:,3)
  WRITE(iun,'("# Reciprocal lattice vectors b_1, b_2, b_3 (a.u.)")')
  WRITE(iun,*) tpiba*bg(:,1), tpiba*bg(:,2), tpiba*bg(:,3)
  WRITE(iun,'("# Number of types of atom")')
  WRITE(iun,*) nsp
  WRITE(iun,'("# Number of atoms")')
  WRITE(iun,*) nat
  WRITE(iun,'("# Atomic species and positions (x, y, z, in a.u.)")')
  DO na =1, nat
     nt = ityp(na)
     WRITE(iun,*) nt
     WRITE(iun,'(a,3e25.15)') atm(nt), alat*tau(:,na)
  END DO
  WRITE(iun,'("# number of G-vectors")')
  WRITE(iun,*) ngm
  WRITE(iun,'("# Miller indices: G=i_1*b_1+i_2*b_2+i_3*b_3")')
  WRITE(iun,'(3i8)') (mill(:,ig), ig=1,ngm)
  WRITE(iun,'("# number of spin components")')
  WRITE(iun,*) nspin
  IF ( nspin == 4 ) THEN
     WRITE(iun,'("# magnetic, spin-orbit?")')
     WRITE(iun,*) domag, lspinorb
  END IF
  CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)
  ALLOCATE (vaux(ngm,nspin))
  CALL rho_r2g (dffts, vrs, vaux )
  WRITE(iun,'("# Local potential V(G) (one column per spin component)")')
  DO is=1,nspin
     WRITE(iun,'("# spin component n.",i4)') is
     ! NOTE: free format is not used to write complex number here and below
     !       complex numbers are written as real and imaginary parts instead
     WRITE(iun,'(2e25.15)') (vaux(ig,is), ig=1,ngm)
  END DO
  DEALLOCATE (vaux)
  WRITE(iun,'("# US PP?")')
  WRITE(iun,*) okvan
  WRITE(iun,'("# Number of projectors")')
  WRITE(iun,*) nh(1:nsp)
  WRITE(iun,'("# Nonlocal PP coefficients Dlm")')
  DO nt = 1, nsp
     WRITE(iun,'("# Atom  type, number of projectors")')
     WRITE(iun,*) nt, nh(nt)
     IF ( lspinorb ) THEN
        DO is = 1, nspin
           WRITE(iun,'("# spin component n.",i4)') is
           WRITE(iun,'(2e25.15)') dvan_so(1:nh(nt),1:nh(nt),is,nt)
        END DO
     ELSE
        WRITE(iun,*) dvan(1:nh(nt),1:nh(nt),nt)
     END IF
  END DO
  WRITE(iun,'("# number of beta functions")')
  WRITE(iun,*) nkb
  !
  WRITE(iun,'("# number of k-points")')
  WRITE(iun,*) nks
  WRITE(iun,'("# number of plane waves for all k-points")')
  WRITE(iun,*) ngk(1:nks)
  DO ik=1,nks
     IF( nspin == 2 ) THEN
        WRITE(iun,'("# k-point n.",i4," spin ",i1)') ik,isk(ik)
     ELSE
        WRITE(iun,'("# k-point n.",i4)') ik
     ENDIF
     WRITE(iun,*) tpiba*xk(:,ik)
     WRITE(iun,'("# number of plane waves")')
     npw = ngk(ik)
     WRITE(iun,*) npw
     WRITE(iun,'("# index of k+G: (k+G)_i = k + G_index(i)")')
     WRITE(iun,'(i8)') (igk_k(ig,ik), ig=1,npw)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     DO ikb=1,nkb
        WRITE(iun,'("# beta function n.",i4)') ikb
        WRITE(iun,'(2e25.15)') vkb(1:npw,ikb)
     END DO
     CALL read_collected_wfc ( restart_dir() , ik, evc )
     WRITE(iun,'("# number of bands")')
     WRITE(iun,*) nbnd
     DO ibnd=1,nbnd
        WRITE(iun,'("# band n.",i4)') ibnd
        WRITE(iun,'("# eigenvalue (eV):",e25.15)') et(ibnd,ik)*13.6058
        WRITE(iun,'(2e25.15)') evc(1:npw,ibnd)
        IF ( nspin == 4 ) &
           WRITE(iun,'(2e25.15)') evc(1+npwx:npw+npwx,ibnd)
     END DO
  END DO
  WRITE(iun,'("# end of file")')
  !
  CLOSE (unit=iun, STATUS='keep')
  !
END SUBROUTINE simple_output
!
!----------------------------------------------------------------------------
SUBROUTINE simple_diag ( fileout )
  !----------------------------------------------------------------------------
  !
  ! Check: read the simple data file, re-diagonalize the matrix 
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), intent(in) :: fileout
  !
  ! These variables are read from file
  INTEGER :: npw, nbnd, nspin, current_spin, ngm, nat, nkb, nsp, npol, nks
  INTEGER, ALLOCATABLE :: ngk(:), ityp(:), nh(:), mill(:,:), igk(:)
  LOGICAL :: lspinorb, domag, okvan
  REAL(dp) :: at(3,3), bg(3,3), xk(3)
  REAL(dp), ALLOCATABLE :: et(:), dvan(:,:,:)
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:), evc(:,:), vkb(:,:), dvan_so(:,:,:,:)
  ! In order to avoid trouble with the format of complex numbers,
  ! these are written as "a b", not "(a,b)"
  REAL(dp), ALLOCATABLE :: raux(:,:)
  REAL(dp) :: dvan_re, dvan_im
  !
  CHARACTER(LEN=80) :: line
  INTEGER :: iun=4, ig, is, ik, ikb, ibnd, na, nt, nt_, i, j, ii, jj, ij, &
       nhm, ih, jh, i1, i2, i3, ipol
  INTEGER, ALLOCATABLE :: limm(:,:,:)
  LOGICAL :: debug = .false., skip_diag = .false.
  REAL(dp) :: g(3)
  REAL(dp), ALLOCATABLE :: e(:)
  COMPLEX(dp), ALLOCATABLE :: h(:,:), v(:,:) 
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Checking simple output data file ",A)' ) &
       trim(fileout)
  OPEN ( UNIT = iun, FORM = 'formatted', STATUS = 'old', FILE = fileout )
  READ(iun,'(a)')  line
  READ(iun,*) at(:,1), at(:,2), at(:,3)
  WRITE(stdout,*) trim(line), at(:,1), at(:,2), at(:,3)
  READ(iun,'(a)')  line
  READ(iun,*) bg(:,1), bg(:,2), bg(:,3)
  WRITE(stdout,*) trim(line), bg(:,1), bg(:,2), bg(:,3)
  READ(iun,'(a)')  line
  READ(iun,*) nsp
  WRITE(stdout,*) trim(line), nsp
  READ(iun,'(a)')  line
  READ(iun,*) nat
  WRITE(stdout,*) trim(line), nat
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  ALLOCATE (ityp(nat))
  DO na =1, nat
     READ(iun,*) nt
     ityp(na) = nt
     READ(iun,'(a)')  line
     ! READ(iun,'(a,3e25.15)') atm(nt), tau(:,na)
     IF (debug) WRITE(stdout,'(a)') line
  END DO
  READ(iun,'(a)')  line
  READ(iun,*) ngm
  WRITE(stdout,*) trim(line), ngm
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  ALLOCATE ( mill(3,ngm) ) 
  READ(iun,'(3i8)') (mill(:,ig), ig=1,ngm)
  ! if i1=mill(1,ig), i2=mill(2,ig), i3=mill(3,ig), then:
  ! limm(i1,i2,i3) = ig
  i1 = MAXVAL( ABS(mill(1,:)) )
  i2 = MAXVAL( ABS(mill(2,:)) )
  i3 = MAXVAL( ABS(mill(3,:)) )
  ALLOCATE (limm(-i1:i1,-i2:i2,-i3:i3))
  limm = 0
  DO ig=1,ngm
     limm( mill(1,ig), mill(2,ig), mill(3,ig) ) = ig
  ENDDO
  !
  READ(iun,'(a)')  line
  READ(iun,*) nspin
  WRITE(stdout,*) trim(line), nspin
  IF ( nspin <= 2 ) THEN
     npol = 1
     domag=.FALSE.
     lspinorb=.FALSE.
  ELSE
     READ(iun,'(a)')  line
     READ(iun,*) domag, lspinorb
     IF ( .NOT.lspinorb ) skip_diag = .true.
     npol = 2
  ENDIF
  ALLOCATE (vaux(ngm,nspin))
  ALLOCATE (raux(2,ngm))
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  DO is=1,nspin
     READ(iun,'(a)')  line
     IF (debug) WRITE(stdout,'(a)') line
     READ(iun,*) (raux(1,ig),raux(2,ig), ig=1,ngm)
     DO ig=1,ngm
        vaux(ig,is) = CMPLX ( raux(1,ig), raux(2,ig) )
     END DO
  END DO
  DEALLOCATE (raux)
  READ(iun,'(a)')  line
  READ(iun,*)  okvan
  IF ( okvan ) skip_diag = .true.
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)',advance='no') trim(line)
  ALLOCATE (nh(nsp))
  READ(iun,*) nh(:)
  WRITE(stdout,*) nh(:)
  nhm = maxval(nh(:))
  IF ( lspinorb ) THEN
     ALLOCATE (dvan_so(nhm,nhm,nspin,nsp) )
  ELSE
     ALLOCATE (dvan(nhm,nhm,nsp) )
  END IF
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  DO nt = 1, nsp
     READ(iun,'(a)')  line
     IF (debug) WRITE(stdout,'(a)') line
     READ(iun,*) nt_, nh(nt)
     IF ( lspinorb ) THEN
        DO is=1,nspin
           READ(iun,'(a)')  line
           IF (debug) WRITE(stdout,'(a)') line
           !READ(iun,*) dvan_so(1:nh(nt),1:nh(nt),is,nt)
           DO jh=1,nh(nt)
              DO ih=1,nh(nt)
                 READ(iun,*) dvan_re, dvan_im
                 dvan_so(ih,jh,is,nt) = CMPLX(dvan_re,dvan_im)
              END DO
           END DO
        END DO
     ELSE
        READ(iun,*) dvan(1:nh(nt),1:nh(nt),nt)
     END IF
  END DO
  READ(iun,'(a)')  line
  READ(iun,*)  nkb
  WRITE(stdout,*) trim(line), nkb
  !
  READ(iun,'(a)')  line
  READ(iun,*) nks
  WRITE(stdout,*) trim(line), nks
  READ(iun,'(a)')  line
  ALLOCATE ( ngk(nks) )
  READ(iun,*) ngk(1:nks)
  DO ik=1,nks
     READ(iun,'(a)')  line
     ! For LSDA, read spin index for this k-point
     IF ( nspin == 2 ) THEN
        READ(line(21:),*) current_spin
        IF (current_spin < 1 .or. current_spin > 2 ) &
            CALL errore ('simple_diag','mismatch in spin index',ik)
     ELSE
        current_spin = 1
     END IF
     READ(iun,*) xk(:)
     WRITE(stdout,*) trim(line), xk(:)
     READ(iun,'(a)')  line
     READ(iun,*) npw
     WRITE(stdout,*) trim(line), npw
     IF ( npw /= ngk(ik) ) CALL errore ('simple_diag','mismatch in Npw',ik)
     READ(iun,'(a)')  line
     WRITE(stdout,'(a)') line
     ALLOCATE ( igk(npw) )
     READ(iun,'(i8)') (igk(ig), ig=1,npw)
     ALLOCATE ( vkb(npw,nkb) )
     ALLOCATE ( raux(2,npw) )
     DO ikb=1,nkb
        READ(iun,'(a)')  line
        IF (debug) WRITE(stdout,'(a)') line
        READ(iun,*) (raux(1,ig), raux(2,ig), ig=1,npw)
        DO ig=1,npw
          vkb(ig,ikb) = CMPLX ( raux(1,ig), raux(2,ig) )
        END DO
     END DO
     DEALLOCATE (raux)
     READ(iun,'(a)')  line
     READ(iun,*)  nbnd
     WRITE(stdout,*) trim(line), nbnd
     ALLOCATE ( et(nbnd), evc(npol*npw,nbnd) )
     ALLOCATE ( raux(2,npol*npw) )
     DO ibnd=1,nbnd
        READ(iun,'(a)')  line
        IF (debug) WRITE(stdout,'(a)') line
        READ(iun,'(a)')  line
        READ(line(19:),*) et(ibnd)
        READ(iun,*) (raux(1,ig), raux(2,ig), ig=1,npol*npw)
        DO ig=1,npol*npw
          evc(ig,ibnd) = CMPLX ( raux(1,ig), raux(2,ig) )
        END DO
     END DO
     DEALLOCATE (raux)
     WRITE(stdout,'("Data read for k-point #",i4,", eigenvalues :")') ik
     WRITE(stdout,'(6f12.6)') et(:)
     ALLOCATE ( h(npol*npw,npol*npw), v(npol*npw,npol*npw), e(npol*npw) )
     h(:,:) = (0.0_dp,0.0_dp)
     DO j=1,npw
        !
        ! kinetic energy
        !
        g(:) = mill(1,igk(j))*bg(:,1) + &
               mill(2,igk(j))*bg(:,2) + &
               mill(3,igk(j))*bg(:,3)
        h(j,j)= ( xk(1)+g(1) )**2 + &
                ( xk(2)+g(2) )**2 + &
                ( xk(3)+g(3) )**2
        IF ( npol == 2 ) h(npw+j,npw+j) = h(j,j)
        !
        ! nonlocal potential
        !
        ikb=0
        DO nt=1,nsp
           DO na=1,nat
              IF ( nt == ityp(na) ) THEN
                 DO ih=1,nh(nt)
                    DO jh=1,nh(nt)
                       IF ( lspinorb ) THEN
                          !
                          ! noncollinear spin-orbit
                          !
                          DO i=1,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,1,nt) *DCONJG(vkb(j,ikb+jh))
                             h(i+npw,j+npw) = h(i+npw,j+npw) + vkb(i,ikb+ih)* &
                                  dvan_so(ih,jh,4,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                          DO i=1,npw
                             h(i,j+npw) = h(i,j+npw) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,2,nt) *DCONJG(vkb(j,ikb+jh))
                             h(i+npw,j) = h(i+npw,j) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,3,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       ELSE IF ( npol == 2 ) THEN
                          !
                          ! noncollinear but no spin-orbit
                          !
                          DO i=1,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan(ih,jh,nt) *DCONJG(vkb(j,ikb+jh))
                             h(i+npw,j+npw) = h(i+npw,j+npw) + vkb(i,ikb+ih)* &
                                  dvan(ih,jh,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       ELSE
                          !
                          ! collinear (LSDA) or unpolarized
                          !
                          DO i=1,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan(ih,jh,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       END IF
                    END DO
                 END DO
                 ikb = ikb + nh(nt)
              END IF
           END DO
        END DO
        !
        ! local potential
        !
        DO i=1,npw
           i1 = mill(1,igk(i)) - mill(1,igk(j))
           i2 = mill(2,igk(i)) - mill(2,igk(j))
           i3 = mill(3,igk(i)) - mill(3,igk(j))
           IF ( ABS(i1) > SIZE(limm,1) .OR. &
                ABS(i2) > SIZE(limm,2) .OR. &
                ABS(i3) > SIZE(limm,3) ) &
                CALL errore ('simple_diag','internal error (1)',i)
           ij = limm ( i1,i2,i3 )
           IF ( ij <= 0 .OR. ij > ngm ) &
                CALL errore ('simple_diag','internal error (2)',i)
           IF (npol == 2) THEN
              ii = npw + i
              jj = npw + j
              IF (domag) THEN
                 !
                 ! noncollinear magnetic
                 !
                 h( i, j) = h( i, j) + vaux(ij,1) + vaux(ij,4)
                 h( i,jj) = h( i,jj) + vaux(ij,2) - (0.0_dp,1.0_dp)*vaux(ij,3)
                 h(ii, j) = h( i,jj) + vaux(ij,2) + (0.0_dp,1.0_dp)*vaux(ij,3)
                 h(ii,jj) = h(ii,jj) + vaux(ij,1) - vaux(ij,4)
              ELSE
                 !
                 ! noncollinear (spin-orbit) not magnetic
                 !
                 h( i, j) = h( i, j) + vaux(ij,1)
                 h(ii,jj) = h(ii,jj) + vaux(ij,1)
              END IF
           ELSE
              !
              ! collinear (LSDA) or unpolarized
              !
              h( i, j) = h( i, j) + vaux(ij,current_spin)
           END IF
        END DO
     END DO
     !
     IF ( .not. skip_diag ) THEN
        CALL cdiagh ( npol*npw, h, npol*npw, e, v)
        WRITE(stdout,'("Recomputed eigenvalues:")')
        WRITE(stdout,'(6f12.6)') e(1:nbnd)*13.6058
     ELSE
        CALL infomsg ('simple_diag','not implemented: diagonalization skipped')
     END IF
     DEALLOCATE ( e, v, h, evc, et, vkb, igk )
  END DO
  DEALLOCATE ( ngk, nh, vaux, limm, mill, ityp )
  IF ( lspinorb ) THEN
     DEALLOCATE ( dvan_so ) 
  ELSE
     DEALLOCATE ( dvan ) 
  END IF
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  !
  CLOSE (unit=iun, STATUS='keep')
  !
END SUBROUTINE simple_diag
