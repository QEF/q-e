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
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: outdir
  INTEGER :: ios
  INTEGER :: first_band, last_band
  NAMELIST / inputpp / outdir, prefix, first_band, last_band
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PMW' )
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  first_band=-1
  last_band=-1
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
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  CALL openfil_pp ( )
  !
  CALL projection( first_band, last_band)
  !
  CALL environment_end ( 'PMW' )
  !
  CALL stop_pp
  !
END PROGRAM pmw

!-----------------------------------------------------------------------
SUBROUTINE projection (first_band, last_band)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE uspp_param, ONLY : upf
  USE ions_base,  ONLY : nat, ityp
  USE basis,      ONLY : natomwfc, swfcatom
  USE cell_base
  USE constants,  ONLY: rytoev
  USE gvect
  USE klist
  USE ldaU,       ONLY : lda_plus_u, &
                         Hubbard_lmax, Hubbard_l, Hubbard_alpha, Hubbard_U
  USE lsda_mod
  USE symm_base,  ONLY: nsym, irt, d1, d2, d3
  USE wvfct
  USE control_flags, ONLY: gamma_only
  USE uspp,       ONLY: nkb, vkb
  USE becmod,     ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files,   ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, &
                        iunsat, nwordatwfc, diropn
  USE wavefunctions_module, ONLY: evc

  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER :: first_band, last_band
  !
  ! local variables
  !
  INTEGER :: ik, na, nt, n, m, l, nwfc, lmax_wfc, &
             ldim1, ldim2, lwork, i, j, info, counter, counter_ldau
  LOGICAL :: exst
  COMPLEX(DP), ALLOCATABLE :: proj (:,:,:)
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
  REAL(DP), ALLOCATABLE :: ew(:), rwork(:)
  ! the eigenvalues of pp
  ! workspace for ZGESVD
  REAL (DP) :: capel
  !
  WRITE( stdout, '(/5x,"Calling projection .... ")')
  IF ( gamma_only ) WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  !
  nwordatwfc = npwx * natomwfc
  CALL diropn( iunsat, 'hub', 2*nwordatwfc, exst )
  !
  ALLOCATE(proj (natomwfc, nbnd, nkstot) )
  ALLOCATE(wfcatom (npwx, natomwfc) )
  ALLOCATE(swfcatom (npwx , natomwfc ) )
  ! Allocate the array containing <beta|wfcatom>
  CALL allocate_bec_type ( nkb, natomwfc, becp)

  IF (first_band == -1)  first_band = 1
  IF (last_band  == -1)  last_band  = nbnd
  IF (first_band > last_band ) CALL errore ('pmw',' first_band > last_band',1)
  IF (first_band < 0         ) CALL errore ('pmw',' first_band < 0 ',       1)
  IF (last_band > nbnd       ) CALL errore ('pmw',' last_band > nbnd ',     1)


  counter = 0
  counter_ldaU = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           IF ( (Hubbard_U(nt)/=0.d0 .or. Hubbard_alpha(nt)/=0.d0) .and. &
                                            l==Hubbard_l(nt) )THEN
               counter_ldaU = counter_ldaU + 2 * l + 1
           ENDIF
           counter = counter + 2 * l + 1
        ENDIF
     ENDDO
  ENDDO

  WRITE( stdout, *) "    NBND = ", nbnd
  WRITE( stdout, *) "    NATOMWFC =", natomwfc
  WRITE( stdout, *) "    NKSTOT =", nkstot

  ldim1 = counter_ldaU
  ldim2 = last_band + 1 - first_band
  WRITE( stdout, *) ldim1, ldim2

  IF (ldim1 > ldim2 ) CALL errore( 'projection','too few bands',ldim1-ldim2)
  lwork = 5 * max(ldim1,ldim2)
  ALLOCATE (pp(ldim1,ldim2), u_m(ldim1,ldim1), w_m(ldim2,ldim2), &
            work(lwork), ew(ldim1), rwork(lwork))
  proj   = 0.d0
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  WRITE (stdout,*) " Hubbard_lmax = ", Hubbard_lmax, lda_plus_u
  nwfc=0
  lmax_wfc = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN
           l = upf(nt)%lchi (n)
           lmax_wfc = max (lmax_wfc, l )
           DO m = 1, 2 * l + 1
              nwfc=nwfc+1
              WRITE(stdout,*) " ATOMIC WFC #", nwfc,":", na,n,l,m
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  IF (lmax_wfc > 3) CALL errore ('projection', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projection', 'wrong # of atomic wfcs?', 1)
  !
  !    loop on k points
  !
  CALL init_us_1
  CALL init_at_1
  !
  DO ik = 1, nks
     WRITE ( stdout, * ) "KPOINT =", ik
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     CALL atomic_wfc (ik, wfcatom)

     CALL init_us_2 (npw, igk, xk (1, ik), vkb)

     CALL calbec ( npw, vkb, wfcatom, becp )

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>

     !
     ! make the projection <psi_i| \hat S | phi_j>
     !
     IF ( gamma_only ) THEN
        ALLOCATE(rproj0(natomwfc,nbnd) )
        CALL calbec ( npw, swfcatom, evc, rproj0 )
        proj(:,:,ik) = cmplx(rproj0(:,:),0.d0,kind=DP)
        DEALLOCATE (rproj0)
     ELSE
        ALLOCATE(proj0(natomwfc,nbnd) )
        CALL calbec ( npw, swfcatom, evc, proj0 )
        proj(:,:,ik) = proj0(:,:)
        DEALLOCATE (proj0)
     ENDIF

     counter = 0
     counter_ldaU = 0
     DO na = 1, nat
        nt = ityp (na)
        DO n = 1, upf(nt)%nwfc
           IF (upf(nt)%oc (n) >= 0.d0) THEN
              l = upf(nt)%lchi (n)
              IF ( (Hubbard_U(nt)/=0.d0.or.Hubbard_alpha(nt)/=0.d0) .and. &
                                            l==Hubbard_l(nt) )THEN
                  pp(counter_ldaU+1:counter_ldaU+2*l+1, 1:ldim2) = &
                      proj(counter+1:counter+2*l+1,first_band:last_band,ik)
                  counter_ldaU = counter_ldaU + 2 * l + 1
              ENDIF
              counter = counter + 2 * l + 1
           ENDIF
        ENDDO
     ENDDO
     IF (counter_ldaU /= ldim1) CALL errore ('projection','wrong counter',1)

     CALL ZGESVD( 'A', 'A', ldim1, ldim2, pp, ldim1, ew, u_m, ldim1, &
                  w_m, ldim2, work, lwork, rwork, info )
     CALL errore ('projection','Singular Value Decomposition failed', abs(info))
     DO i = 1, ldim1
        WRITE ( stdout, * ) ew(i)
        WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        WRITE ( stdout, '(8(2f5.2,2x))') w_m(i,:)
     ENDDO
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

     IF (capel < 1.d-10) THEN
        WRITE ( stdout, *) " ORTHOGONALITY CHECK PASSED "
     ELSE
        WRITE ( stdout, *) " ORTHOGONALITY CHECK FAILED"
        WRITE ( stdout, *) " CAPEL = ", capel
        DO i=1,ldim1
           WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        ENDDO
     ENDIF
     counter = 0
     counter_ldaU = 0
     DO na = 1, nat
        nt = ityp (na)
        DO n = 1, upf(nt)%nwfc
           IF (upf(nt)%oc (n) >= 0.d0) THEN
              l = upf(nt)%lchi (n)
              IF ( (Hubbard_U(nt)/=0.d0.or.Hubbard_alpha(nt)/=0.d0) .and. &
                                            l==Hubbard_l(nt) )THEN
                  CALL zgemm( 'N', 'C', npw, 2*l+1, ldim2, ONE, &
                              evc(1,first_band), npwx, &
                              pp(counter_ldaU+1,1), ldim1, ZERO, &
                              wfcatom(1,counter+1), npwx )
                  counter_ldaU = counter_ldaU + 2 * l + 1
              ENDIF
              counter = counter + 2 * l + 1
           ENDIF
        ENDDO
     ENDDO

     CALL calbec ( npw, vkb, wfcatom, becp )

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     CALL davcio (swfcatom, 2*nwordatwfc, iunsat, ik, 1)

     ! on k-points
  ENDDO
  !
  CALL deallocate_bec_type (becp)
  !
  DEALLOCATE (pp, u_m, w_m, work, ew, rwork)
  DEALLOCATE (swfcatom)
  DEALLOCATE (wfcatom)
  DEALLOCATE (proj)

  RETURN
END SUBROUTINE projection
