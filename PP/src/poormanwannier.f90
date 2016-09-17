!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ONE  (1.D0,0.D0)
#define ZERO (0.D0,0.D0)
!
!-----------------------------------------------------------------------
PROGRAM pmw
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto atomic wavefunctions,
  !
  ! input: namelist "&inputpp", with variables
  !   prefix      prefix of input files saved by program pwscf
  !   outdir      temporary directory where files resides
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : prefix, tmp_dir
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, nproc_pool, nproc_pool_file
  USE environment, ONLY : environment_start, environment_end
  USE wvfct,      ONLY : nbnd
  USE ldaU,       ONLY : lda_plus_U
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: outdir
  INTEGER :: ios
  INTEGER :: first_band, last_band
  REAL(DP) :: min_energy, max_energy, sigma
  LOGICAL :: writepp
  NAMELIST / inputpp / outdir, prefix, first_band, last_band, writepp, &
                      min_energy, max_energy, sigma
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PMW' )
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  first_band=-1
  last_band =-1
  min_energy = -9.d99
  max_energy =  9.d99
  sigma      = -1.d0
  writepp = .FALSE.
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios/=0 ) CALL errore ('pmwannier', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( first_band, ionode_id, world_comm )
  CALL mp_bcast( last_band, ionode_id, world_comm )
  CALL mp_bcast( min_energy, ionode_id, world_comm )
  CALL mp_bcast( max_energy, ionode_id, world_comm )
  CALL mp_bcast( sigma, ionode_id, world_comm )
  CALL mp_bcast( writepp, ionode_id, world_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  ! Here we trap restarts from a different number of nodes.
  !
  IF (nproc_pool /= nproc_pool_file)  &
     CALL errore('pmw', &
     'pw.x run on a different number of procs/pools', nproc_pool_file)
  !
  ! Check on correctness and consistency of the input
  !
  IF ( first_band == -1 )  first_band = 1
  IF ( last_band  == -1 )  last_band  = nbnd
  IF ( first_band > last_band ) CALL errore ('pmw',' first_band > last_band',1)
  IF ( first_band < 0 ) CALL errore ('pmw',' first_band < 0 ', first_band)
  IF ( last_band > nbnd ) CALL errore ('pmw',' last_band > nbnd ', nbnd)
  IF ( sigma > 0.d0 ) THEN
     IF ( min_energy > max_energy ) CALL errore ('pmw',' min_energy > max_energy',1)
  END IF
  ! Check on compatibilities
  IF ( noncolin ) CALL errore('pmw','non-colinear not implemented / not tested', 1)
  ! Currently, WF projectors are built for Hubbard species only
  IF ( .NOT.lda_plus_U ) CALL errore('pmw','Hubbard U calculation required', 1)
  !
  CALL openfil_pp ( )
  !
  CALL projection( first_band, last_band, min_energy, max_energy, sigma, writepp)
  !
  CALL environment_end ( 'PMW' )
  !
  CALL stop_pp
  !
END PROGRAM pmw

!-----------------------------------------------------------------------
SUBROUTINE projection (first_band, last_band, min_energy, max_energy, sigma, iopp)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout, ionode
  USE uspp_param, ONLY : upf
  USE ions_base,  ONLY : nat, ityp
  USE basis,      ONLY : natomwfc, swfcatom
  USE cell_base
  USE constants,  ONLY: rytoev
  USE gvect
  USE klist
  USE wvfct,      ONLY : nbnd, npwx, et
  USE ldaU,       ONLY : is_Hubbard, Hubbard_lmax, Hubbard_l, &
                         oatwfc, offsetU, nwfcU, wfcU, copy_U_wfc
  USE lsda_mod
  USE symm_base
  USE mp_pools,   ONLY : me_pool, root_pool, my_pool_id, kunit, npool
  USE control_flags, ONLY: gamma_only
  USE uspp,       ONLY: nkb, vkb
  USE becmod,     ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files,   ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, &
                        iunhub, nwordwfcU, nwordatwfc, diropn
  USE wavefunctions_module, ONLY: evc

  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  REAL(DP), EXTERNAL :: wgauss
  !
  ! I/O variables
  !
  INTEGER, INTENT(IN) :: first_band, last_band
  REAL(DP), INTENT(IN) :: min_energy, max_energy, sigma
  LOGICAL, INTENT(IN) :: iopp
  !
  ! local variables
  !
  INTEGER :: npw, ibnd, ik, na, nt, n, m, l, nwfc, lmax_wfc, &
             ldim1, ldim2, lwork, i, j, info, counter, counter_ldau
  LOGICAL :: exst
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  !
  COMPLEX(DP), ALLOCATABLE ::  proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE :: rproj0(:,:)
  ! ... or for gamma-point.
  COMPLEX(DP), ALLOCATABLE :: pp(:,:), u_m(:,:), w_m(:,:), work(:)
  ! the overlap matrix pp
  ! left unitary matrix in the SVD of sp_m
  ! right unitary matrix in the SVD of sp_m
  ! workspace for ZGESVD
  REAL(DP), ALLOCATABLE :: ew(:), rwork(:), gk(:)
  ! the eigenvalues of pp
  ! workspace for ZGESVD
  REAL (DP) :: capel
  REAL (DP) :: e
  !!
  INTEGER :: iun_pp, nks1, nks1tot, nks2, nks2tot, nbase, rest
  CHARACTER(len=9) :: kptstr
  INTEGER, ALLOCATABLE :: kstatus(:)
  !!
  !
  WRITE( stdout, '(/5x,"Calling projection to compute Wannier projectors ... ")')
  !
  IF ( gamma_only ) WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  !
  ! Set dimensions of the problem:
  ! ldim1 = number of atomic-projectors
  ! ldim2 = number of KS-bands
  !! N.B.: ldim1 = nwfcU = size(wfcU,2), in current implementation
  ldim1 = nwfcU
  ldim2 = last_band + 1 - first_band
  IF (ldim1 > ldim2 ) CALL errore( 'projection','too few bands',ldim1-ldim2)
  !
  lmax_wfc = Hubbard_lmax
  IF (lmax_wfc > 3) CALL errore ('projection', 'l > 3 not yet implemented', 1)
  !
  ! write initial summary
  IF ( ionode ) THEN
     WRITE( stdout, '(/,6(5x,A,5x,I8,/))') &
        'number of k-points            ', nkstot, &
        'number of atomic wfcs with U  ', nwfcU, &
        'total number of atomic wfcs   ', natomwfc, &
        'number of Wannier functions   ', ldim1, &
        'number of bands for projecting', ldim2, &
        'total number bands            ', nbnd
     IF (sigma > 0.d0) THEN
        WRITE( stdout, '(5x,A,/,3(5x,A,f10.4,A,/))') &
           'projection limited to an energy window',&
           'smoothing simgma                     ', sigma, ' eV', &
           'minimum energy                       ', min_energy,' eV', &
           'maximum energy                       ', max_energy,' eV'
     END IF
  ENDIF
  !
  ! initializations needed to get correct output when npool>1
  ALLOCATE( kstatus(nkstot) )
  kstatus(:) = -1
  nks1 = kunit * ( nkstot / kunit / npool )
  rest = ( nkstot - nks1 * npool ) / kunit
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  IF ( nks1 .NE. nks ) call errore('projection','problems with nks1',1)
  nbase = nks * my_pool_id
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! Delete .hub files (if present) to prevent size mismatch or other problems
  INQUIRE( UNIT=iunhub, OPENED=exst )
  IF ( .NOT. exst ) CALL diropn( iunhub, 'hub', 2*nwordwfcU, exst )
  IF ( ionode .AND. exst ) WRITE( stdout, '(5x,A)' ) '.hub files will be overwritten'
  CLOSE( UNIT=iunhub, STATUS='delete' )
  ! Create and open output files to store the projectors
  CALL diropn( iunhub, 'hub', 2*nwordwfcU, exst )
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  ! Allocate temporary arrays and workspace
  !
  lwork = 5 * max(ldim1,ldim2)
  ALLOCATE (pp(ldim1,ldim2), u_m(ldim1,ldim1), w_m(ldim2,ldim2), &
            work(lwork), ew(ldim1), rwork(lwork))
  pp   = 0.d0
  ALLOCATE(wfcatom (npwx, natomwfc) )
  ALLOCATE(swfcatom (npwx , ldim1 ) )
  ! Allocate the array containing <beta|wfcatom>
  CALL allocate_bec_type ( nkb, ldim1, becp)
  !
  ! Main loop (on k-points)
  !
  ALLOCATE (gk(npwx))
  DO ik = 1, nks
     !
     npw = ngk(ik)

     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     CALL atomic_wfc (ik, wfcatom)

     ! select Hubbard wavefunctions and copy then into wfcU (ldaU module)
     ! (in order to avoid computing S|phi> for non-Hubbard atomic wfcs)
     CALL copy_U_wfc (wfcatom)
     !CALL copy_U_wfc (swfcatom, noncolin) ! not yet implemented/tested

     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

     CALL calbec ( npw, vkb, wfcU, becp )

     CALL s_psi (npwx, npw, ldim1, wfcU, swfcatom)
     !
     ! wfcU = |phi_i> , swfcatom = \hat S |phi_i>
     !
     ! make the projection <psi_i| \hat S |phi_j>
     !
     IF ( gamma_only ) THEN
        ALLOCATE(rproj0(ldim1,nbnd) )
        CALL calbec ( npw, swfcatom, evc, rproj0 )
        pp(:,:) = cmplx(rproj0(:,first_band:last_band),0.d0,kind=DP)
        DEALLOCATE (rproj0)
     ELSE
        ALLOCATE(proj0(ldim1,nbnd) )
        CALL calbec ( npw, swfcatom, evc, proj0 )
        pp(:,:) = proj0(:,first_band:last_band)
        DEALLOCATE (proj0)
     ENDIF
!
! add a damping factor if we want to select an energy window
!
     if (sigma > 0.d0) then 
        do i=1,ldim2
           ibnd = i + first_band -1
           e = et(ibnd,ik) * rytoev
           pp(:,i) = pp(:,i) * wgauss((e-min_energy)/sigma,0) * &
                               wgauss((max_energy-e)/sigma,0)
        end do
     end if

     !
     ! Use S.V.D. to make the orthonormalization of projectors
     !
     CALL ZGESVD( 'A', 'A', ldim1, ldim2, pp, ldim1, ew, u_m, ldim1, &
                  w_m, ldim2, work, lwork, rwork, info )
     IF ( abs(info) .ne. 0 ) kstatus(ik) = -2
     !CALL errore ('projection','Singular Value Decomposition failed', abs(info))
     !!DEBUG
!      WRITE ( * , * ) (ew(i),i=1,ldim1)
     !DO i = 1, ldim1
     !   WRITE ( * , * ) ew(i)
     !   WRITE ( * , '(8(2f5.2,2x))') u_m(:,i)
     !   WRITE ( * , '(8(2f5.2,2x))') w_m(i,:)
     !ENDDO
     !
     ! ... use sp_m to store u_m * w_m
     !
     CALL zgemm( 'N', 'N', ldim1, ldim2, ldim1, ONE, u_m, ldim1, w_m, &
                    ldim2, ZERO, pp, ldim1 )
     ! ... check orthogonality
     CALL zgemm( 'N', 'C', ldim1, ldim1, ldim2, ONE, pp, ldim1, pp, &
                    ldim1, ZERO, u_m, ldim1 )
     capel = 0.d0
     DO i=1,ldim1
        u_m(i,i) = u_m(i,i) -1.d0
        DO j=1,ldim1
           capel = capel + abs( u_m(i,j) )
        ENDDO
        u_m(i,i) = u_m(i,i) +1.d0
     ENDDO

     IF ( kstatus(ik) == -1 ) THEN 
        IF ( capel > 1.d-10 ) THEN
           kstatus(ik) = -3
           !!DEBUG
           !WRITE (*,*) " ORTHOGONALITY CHECK FAILED"
           !WRITE (*,*) " CAPEL = ", capel
           !DO i=1,ldim1
           !   WRITE (*, '(8(2f5.2,2x))') u_m(:,i)
           !ENDDO
        ELSE
           kstatus(ik) = 0
        ENDIF
     ENDIF

     ! ... compute wave functions |ophi_j> for the orthonormalized projectors
     CALL zgemm( 'N', 'C', npw, ldim1, ldim2, ONE, evc(1,first_band), npwx, &
                  pp(1,1), ldim1, ZERO, wfcU, npwx )

     CALL calbec ( npw, vkb, wfcU, becp )

     CALL s_psi (npwx, npw, ldim1, wfcU, swfcatom)

     ! write \hat S |ophi_j> into iunhub unit
     CALL davcio (swfcatom, 2*nwordwfcU, iunhub, ik, 1)

     ! ... write U matrices to disk, if required
     ! (each k-point into a separate file named <prefix.pp?>)
     IF ( me_pool == root_pool .AND. iopp ) THEN
        WRITE(kptstr,'(I9)') nbase + ik
        iun_pp = find_free_unit()
        OPEN (unit=iun_pp, file=trim(prefix)//'.pp'//trim(adjustl(kptstr)), &
           form='formatted')

        DO i = 1,ldim1
           WRITE(iun_pp,'(2f22.15)') ( pp(i,j), j=1,ldim2 )
           WRITE(iun_pp,'(/)')
        ENDDO

        CLOSE(iun_pp)
     ENDIF
     !
     ! on k-points
  ENDDO
  DEALLOCATE (gk)
  !
  ! Check if everything went OK (on all k-points)
  !
  CALL ipoolrecover( kstatus, 1, nkstot, nks)
  IF ( ionode ) WRITE(stdout, '(/,5x,A)') 'Orthogonality check'
  DO ik = 1,nkstot
     IF ( ionode ) WRITE(stdout, '(7x,A,I8,A)', advance='no') &
        'k-point', ik, ' : '
     SELECT CASE ( kstatus(ik) )
        CASE (  0 )
           IF ( ionode ) WRITE( stdout, '(A)' ) 'ok'
        CASE ( -1 )
           CALL errore ('projection','kstatus not set?', -1)
        CASE ( -2 )
           !IF ( ionode ) WRITE( stdout, '(A)' ) 'Singular Value Dec. *FAILED*'
           CALL errore ('projection','Singular Value Decomposition failed', -2)
        CASE ( -3 )
           !IF ( ionode ) WRITE( stdout, '(A)' ) '*FAILED*'
           CALL errore ('projection','orthogonality check failed', ik)
        CASE DEFAULT
           CALL errore ('projection','invalid kstatus code', kstatus(ik))
     END SELECT
  ENDDO
  !
  ! if required, write symmetries to disk (into file prefix.pp0)
  IF ( ionode .AND. iopp ) THEN
     iun_pp = find_free_unit()
     OPEN (unit=iun_pp, file=trim(prefix)//".pp0", form='formatted')
     WRITE(iun_pp,'("# nrot,nsym,nsym_ns,nsym_na",4I4)') nrot, nsym, nsym_ns, nsym_na
     DO i = 1,nsym
        WRITE(iun_pp,'("#symm",I3," : ",A)') i, trim(sname(i))
        WRITE(iun_pp,'(3I3,I7,5x,3F7.2,F9.2)') ( s(j,:,i), ftau(j,i), sr(j,:,i), ft(j,i), j=1,3 )
        WRITE(iun_pp,'(99I3)') irt(i,1:nat)
     ENDDO
     CLOSE(iun_pp)
  ENDIF

  !
  CALL deallocate_bec_type (becp)
  !
  DEALLOCATE (pp, u_m, w_m, work, ew, rwork)
  DEALLOCATE (swfcatom)
  DEALLOCATE (wfcatom)

  RETURN
END SUBROUTINE projection
