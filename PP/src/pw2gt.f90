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
  ! The output file is a single, formatted and self-contained file
  ! (BEWARE! may be quite large!) that can be useful for anybody else
  ! not willing to figure out how read QE data files.
  ! IMPORTANT NOTICE: must be run on a single processor ONLY.
  ! Uncomment next line to test the file: the code will re-read the file
  ! just written, fill and diagonalize the hamiltonian matrix
#define __DEBUG
  !
  ! Input: namelist &inputpp [outdir=...] [prefix=...] / as in QE input
  ! (default values as in QE).
  ! Output file written in "outdir"/"prefix".save/output.dat
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : tmp_dir, prefix
  USE mp_global,  ONLY : mp_startup
  USE mp_images,  ONLY : intra_image_comm, nproc_image
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: needwf = .true.
  INTEGER :: ios
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / inputpp / outdir, prefix
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
  !   Read xml file, allocate and initialize general variables
  !
  CALL read_file_new ( needwf )
  !
  CALL simple_output ()
  !
#if defined (__DEBUG)
  CALL simple_diag ()
#endif
  !  
  CALL environment_end ( 'PW2GT' )
  !
  CALL stop_pp()
  !
END PROGRAM pw2gt
!----------------------------------------------------------------------------
SUBROUTINE simple_output (  )
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
  USE lsda_mod,     ONLY : nspin
  USE klist,        ONLY : nks, xk, ngk, igk_k
  USE uspp,         ONLY : okvan, nkb, vkb, dvan, dvan_so
  USE uspp_param,   ONLY : nh
  USE uspp_init,    ONLY : init_us_2
  USE fft_base,     ONLY : dfftp, dffts
  USE wvfct,        ONLY : nbnd, et
  USE wavefunctions,     ONLY : evc
  USE noncollin_module,  ONLY : lspinorb, domag
  USE pw_restart_new,    ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:)
  CHARACTER(LEN=256) :: fileout
  INTEGER :: iun, ig, is, ik, ikb, ibnd, na, nt, npw
  !
  fileout = TRIM ( restart_dir() ) //  'output.dat'
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing simple output data file ",A)' ) &
       fileout
  OPEN ( NEWUNIT = iun, FORM = 'formatted', STATUS = 'unknown', &
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
     IF ( nspin /= 4 ) THEN
        WRITE(iun,*) dvan(1:nh(nt),1:nh(nt),nt)
     ELSE
        DO is = 1, nspin
           WRITE(iun,'("# spin component n.",i4)') is
           WRITE(iun,*) dvan_so(1:nh(nt),1:nh(nt),is,nt)
        END DO
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
     WRITE(iun,'("# k-point n.",i4)') ik
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
     END DO
  END DO
  WRITE(iun,'("# end of file")')
  !
  CLOSE (unit=iun, STATUS='keep')
  !
END SUBROUTINE simple_output
!
!----------------------------------------------------------------------------
SUBROUTINE simple_diag (  )
  !----------------------------------------------------------------------------
  !
  ! Check: read the simple data file, re-diagonalize the matrix 
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : restart_dir
  USE cell_base,    ONLY : at, bg, alat, tpiba
  USE ions_base,    ONLy : ityp, atm, tau
  USE gvect,        ONLY : mill
  USE klist,        ONLY : xk, igk_k
  USE uspp,         ONLY : vkb, dvan, dvan_so
  USE uspp_param,   ONLY : nh
  USE wavefunctions,ONLY : evc
  !
  IMPLICIT NONE
  !
  ! These variables are read from file
  INTEGER :: npw, nbnd, nspin, ngm, nat, nkb, nsp, npol, nks
  LOGICAL :: lspinorb, domag, okvan
  REAL(dp), ALLOCATABLE :: et(:)
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:)
  COMPLEX(dp) :: gfortran_merda
  !
  CHARACTER(LEN=256) :: fileout
  CHARACTER(LEN=80) :: line
  COMPLEX(dp), ALLOCATABLE :: h(:,:), v(:,:) 
  INTEGER :: iun, ig, is, ik, ikb, ibnd, na, nt, nt_, i, j, ii, jj, ij, &
       ih, jh, i1, i2, i3, ipol
  INTEGER, ALLOCATABLE :: limm(:,:,:)
  REAL(dp) :: g(3)
  LOGICAL :: debug = .true.
  !
  fileout = TRIM ( restart_dir() ) //  'output.dat'
  WRITE( UNIT = stdout, FMT = '(/,5X,"Checking simple output data file ",A)' ) &
       fileout
  OPEN ( NEWUNIT = iun, FORM = 'formatted', STATUS = 'old', FILE = fileout )
  READ(iun,'(a)')  line
  READ(iun,*) at(:,1), at(:,2), at(:,3)
  IF ( debug) WRITE(stdout,*) trim(line), at(:,1), at(:,2), at(:,3)
  READ(iun,'(a)')  line
  READ(iun,*) bg(:,1), bg(:,2), bg(:,3)
  IF ( debug) WRITE(stdout,*) trim(line), bg(:,1), bg(:,2), bg(:,3)
  READ(iun,'(a)')  line
  READ(iun,*) nsp
  IF ( debug) WRITE(stdout,*) trim(line), nsp
  READ(iun,'(a)')  line
  READ(iun,*) nat
  IF ( debug) WRITE(stdout,*) trim(line), nat
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO na =1, nat
     READ(iun,'(a)')  line
     READ(iun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
  END DO
  READ(iun,'(a)')  line
  READ(iun,*) ngm
  IF ( debug) WRITE(stdout,*) trim(line), ngm
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
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
  IF ( debug) WRITE(stdout,*) trim(line), nspin
  IF ( nspin <= 2 ) THEN
     npol = 1
  ELSE
     READ(iun,'(a)')  line
     READ(iun,*) domag, lspinorb
     IF ( domag .OR. .NOT.lspinorb ) &
          CALL errore ('simple_diag','spin-orbit with no magnetization only',1)
     npol = 2
  ENDIF
  ALLOCATE (vaux(ngm,nspin))
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO is=1,nspin
     READ(iun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iun,'(2e25.15)') (vaux(ig,is), ig=1,ngm)
     ! should be READ(iun,*) (vaux(ig,is), ig=1,ngm)
  END DO
  READ(iun,'(a)')  line
  READ(iun,*)  okvan
  IF ( okvan ) CALL errore ('simple_diag','US PP not implemented',1)
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  READ(iun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO nt = 1, nsp
     READ(iun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iun,*) nt_, nh(nt)
     IF ( nspin /= 4 ) THEN
        READ(iun,*) dvan(1:nh(nt),1:nh(nt),nt)
     ELSE
        DO is=1,nspin
           READ(iun,'(a)')  line
           IF ( debug) WRITE(stdout,'(a)') line
           READ(iun,*) dvan_so(1:nh(nt),1:nh(nt),is,nt)
        END DO
     END IF
  END DO
  READ(iun,'(a)')  line
  READ(iun,*)  nkb
  IF ( debug) WRITE(stdout,*) trim(line), nkb
  !
  READ(iun,'(a)')  line
  READ(iun,*) nks
  IF ( debug) WRITE(stdout,*) trim(line), nks
  READ(iun,'(a)')  line
  READ(iun,'(a)')  line
  DO ik=1,nks
     READ(iun,'(a)')  line
     READ(iun,*) xk(:,ik)
     IF ( debug) WRITE(stdout,*) trim(line), xk(:,ik)
     READ(iun,'(a)')  line
     READ(iun,*) npw
     IF ( debug) WRITE(stdout,*) trim(line), npw
     READ(iun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iun,'(i8)') (igk_k(ig,ik), ig=1,npw)
     DO ikb=1,nkb
        READ(iun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iun,'(2e25.15)') vkb(1:npw,ikb)
     END DO
     READ(iun,'(a)')  line
     READ(iun,*)  nbnd
     IF ( debug) WRITE(stdout,*) trim(line), nbnd
     DO ibnd=1,nbnd
        READ(iun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iun,'(2e25.15)') evc(1:npw,ibnd)
     END DO
     WRITE(stdout,'("# Data read for k-point #",i4,", diagonalizing...")') ik
     ALLOCATE ( h(npol*npw,npol*npw), v(npol*npw,npol*npw), et(npol*npw) )
     h(:,:) = (0.0_dp,0.0_dp)
     DO j=1,npw
        !
        ! kinetic energy
        !
        g(:) = mill(1,igk_k(j,ik))*bg(:,1) + &
               mill(2,igk_k(j,ik))*bg(:,2) + &
               mill(3,igk_k(j,ik))*bg(:,3)
        h(j,j)= ( xk(1,ik)+g(1) )**2 + &
                ( xk(2,ik)+g(2) )**2 + &
                ( xk(3,ik)+g(3) )**2
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
                       IF ( nspin /= 4 ) THEN
                          DO i=j,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan(ih,jh,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       ELSE
                          DO i=j,npw
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
        DO i=j,npw
           i1 = mill(1,igk_k(i,ik)) - mill(1,igk_k(j,ik))
           i2 = mill(2,igk_k(i,ik)) - mill(2,igk_k(j,ik))
           i3 = mill(3,igk_k(i,ik)) - mill(3,igk_k(j,ik))
           IF ( ABS(i1) > SIZE(limm,1) .OR. &
                ABS(i2) > SIZE(limm,2) .OR. &
                ABS(i3) > SIZE(limm,3) ) &
                CALL errore ('simple_diag','internal error (1)',i)
           ij = limm ( i1,i2,i3 )
           IF ( ij <= 0 .OR. ij > ngm ) &
                CALL errore ('simple_diag','internal error (2)',i)
           DO ipol = 1, npol
              ii = (ipol-1)*npw + i
              jj = (ipol-1)*npw + j
              h(ii,jj) = h(ii,jj) + vaux(ij,1)
              IF ( i > j ) h(jj,ii) =DCONJG(h(ii,jj))
           END DO
        END DO
     END DO
     !
     CALL cdiagh ( npol*npw, h, npol*npw, et, v)
     WRITE(stdout,'(4f12.6)') (et(ibnd)*13.6058, ibnd=1,nbnd)
     DEALLOCATE ( et, v, h )
  END DO
  DEALLOCATE (limm, vaux)
  READ(iun,'(a)')  line
  WRITE(stdout,'(a)') line
  !
  CLOSE (unit=iun, STATUS='keep')
  !
END SUBROUTINE simple_diag
