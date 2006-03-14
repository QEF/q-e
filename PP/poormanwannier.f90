!
! Copyright (C) 2001-2003 PWSCF group 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#include "f_defs.h" 
#define ONE  (1.D0,0.D0)
#define ZERO (0.D0,0.D0)  
!
!----------------------------------------------------------------------- 
PROGRAM poormanwannier 
  !----------------------------------------------------------------------- 
  ! 
  ! projects wavefunctions onto atomic wavefunctions, 
  ! 
  ! input: namelist "&inputpp", with variables 
  !   prefix      prefix of input files saved by program pwscf 
  !   outdir      temporary directory where files resides 
  ! 
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE kinds,      ONLY : DP 
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir, trimcheck
  USE mp,         ONLY : mp_bcast
  !
  IMPLICIT NONE 
  CHARACTER(len=256) :: outdir
  INTEGER :: ios
  INTEGER :: first_band, last_band
  NAMELIST / inputpp / outdir, prefix, first_band, last_band
  ! 
  CALL start_postproc (nd_nmbr) 
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf' 
  outdir = './' 
  first_band=-1
  last_band=-1
  ! 
  IF ( ionode )  THEN 
     ! 
     READ (5, inputpp, err = 200, iostat = ios) 
200  CALL errore ('pmwannier', 'reading inputpp namelist', ABS (ios) ) 
     ! 
     tmp_dir = trimcheck (outdir) 
     ! 
  END IF 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( first_band, ionode_id ) 
  CALL mp_bcast( last_band, ionode_id )  
  ! 
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  CALL read_file ( )
  !
  CALL openfil_pp ( ) 
  !
  CALL projection( first_band, last_band)
  ! 
  CALL stop_pp 
  ! 
END PROGRAM poormanwannier
 
!----------------------------------------------------------------------- 
SUBROUTINE projection (first_band, last_band)
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout 
  USE atom 
  USE ions_base, ONLY : nat, ityp
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev 
  USE gvect 
  USE klist 
  USE ldaU 
  USE lsda_mod 
  USE symme, ONLY: nsym, irt 
  USE wvfct 
  USE uspp, ONLY: nkb, vkb
  USE becmod,   ONLY: becp, rbecp
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, &
                      iunat, nwordatwfc
  USE wavefunctions_module, ONLY: evc 

  IMPLICIT NONE 
  !
  ! I/O variables 
  !
  INTEGER :: first_band, last_band
  !
  ! local variables
  !
  INTEGER :: ik, ia, ib, na, nt, n, m, l, nwfc, lmax_wfc, &
             ldim1, ldim2, lwork, i, j, info, counter, counter_ldau
  LOGICAL :: exst 
  COMPLEX(DP), ALLOCATABLE :: proj (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:) 
  INTEGER, ALLOCATABLE :: INDEX(:) 
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
  iunat = 13
  nwordatwfc = 2 * npwx * natomwfc
  CALL diropn( iunat, 'atwfc', nwordatwfc, exst )
  !
  ALLOCATE(proj (natomwfc, nbnd, nkstot) ) 
  ALLOCATE(wfcatom (npwx, natomwfc) ) 
  ! Allocate the array containing <beta|wfcatom>
  IF ( gamma_only ) THEN 
     ALLOCATE (rbecp (nkb,natomwfc)) 
  ELSE
     ALLOCATE ( becp (nkb,natomwfc)) 
  END IF 

  IF (first_band == -1)  first_band = 1
  IF (last_band  == -1)  last_band  = nbnd
  IF (first_band > last_band ) CALL errore ('pmw',' first_band > last_band',1)
  IF (first_band < 0         ) CALL errore ('pmw',' first_band < 0 ',       1)
  IF (last_band > nbnd       ) CALL errore ('pmw',' last_band > nbnd ',     1)
 

  counter = 0
  counter_ldaU = 0
  DO na = 1, nat
     nt = ityp (na)
     DO n = 1, nchi (nt)
        IF (oc (n, nt) >= 0.d0) THEN
           l = lchi (n, nt)
           IF ( (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) .AND. &
                                            l.EQ.Hubbard_l(nt) )THEN
               counter_ldaU = counter_ldaU + 2 * l + 1
           END IF
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
  lwork = 5 * MAX(ldim1,ldim2)
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
     DO n = 1, nchi (nt) 
        IF (oc (n, nt) >= 0.d0) THEN 
           l = lchi (n, nt) 
           lmax_wfc = MAX (lmax_wfc, l ) 
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
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
 
     CALL atomic_wfc (ik, wfcatom) 
 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb) 
 
     IF ( gamma_only ) THEN 
        CALL pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, wfcatom, npwx, rbecp, nkb)   
     ELSE 
        CALL ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom) 
     END IF 
 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
 
     ! 
     ! make the projection <psi_i| \hat S | phi_j> 
     ! 
     IF ( gamma_only ) THEN 
        ALLOCATE(rproj0(natomwfc,nbnd) ) 
        CALL pw_gemm ('Y', natomwfc, nbnd, npw, swfcatom, npwx, evc, npwx, &
                       rproj0, natomwfc) 
        proj(:,:,ik) = CMPLX(rproj0(:,:),0.d0)
        DEALLOCATE (rproj0) 
     ELSE 
        ALLOCATE(proj0(natomwfc,nbnd) ) 
        CALL ccalbec (natomwfc, npwx, npw, nbnd, proj0, swfcatom, evc) 
        proj(:,:,ik) = proj0(:,:)
        DEALLOCATE (proj0) 
     END IF 

     counter = 0
     counter_ldaU = 0
     DO na = 1, nat
        nt = ityp (na)
        DO n = 1, nchi (nt)
           IF (oc (n, nt) >= 0.d0) THEN
              l = lchi (n, nt)
              IF ( (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) .AND. &
                                            l.EQ.Hubbard_l(nt) )THEN
                  pp(counter_ldaU+1:counter_ldaU+2*l+1, 1:ldim2) = &
                      proj(counter+1:counter+2*l+1,first_band:last_band,ik)
                  counter_ldaU = counter_ldaU + 2 * l + 1
              END IF
              counter = counter + 2 * l + 1
           ENDIF
        ENDDO
     ENDDO
     IF (counter_ldaU .NE. ldim1) CALL errore ('projection','wrong counter',1)

     CALL ZGESVD( 'A', 'A', ldim1, ldim2, pp, ldim1, ew, u_m, ldim1, &
                  w_m, ldim2, work, lwork, rwork, info )
     CALL errore ('projection','Singular Value Deconposition failed', ABS(info))
     DO i = 1, ldim1
        WRITE ( stdout, * ) ew(i)
        WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        WRITE ( stdout, '(8(2f5.2,2x))') w_m(i,:)
     END DO
     !
     ! ... use sp_m to store u_m * w_m
     !
     CALL ZGEMM( 'N', 'N', ldim1, ldim2, ldim1, ONE, u_m, ldim1, w_m, &
                    ldim2, ZERO, pp, ldim1 )
     ! ... check orthogonality
     CALL ZGEMM( 'N', 'C', ldim1, ldim1, ldim2, ONE, pp, ldim1, pp, &
                    ldim1, ZERO, u_m, ldim1 )
     capel = 0.d0
     DO i=1,ldim1
        u_m(i,i) = u_m(i,i) -1.d0
        DO j=1,ldim1
           capel = capel + ABS( u_m(i,j) )
        END DO
        u_m(i,i) = u_m(i,i) +1.d0
     END DO

     IF (capel < 1.d-10) THEN
        WRITE ( stdout, *) " ORTHOGONALITY CHECK PASSED "
     ELSE
        WRITE ( stdout, *) " ORTHOGONALITY CHECK FAILED"
        WRITE ( stdout, *) " CAPEL = ", capel
        DO i=1,ldim1
           WRITE ( stdout, '(8(2f5.2,2x))') u_m(:,i)
        END DO
     END IF
     counter = 0
     counter_ldaU = 0
     DO na = 1, nat
        nt = ityp (na)
        DO n = 1, nchi (nt)
           IF (oc (n, nt) >= 0.d0) THEN
              l = lchi (n, nt)
              IF ( (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) .AND. &
                                            l.EQ.Hubbard_l(nt) )THEN
                  CALL ZGEMM( 'N', 'C', npw, 2*l+1, ldim2, ONE, &
                              evc(1,first_band), npwx, &
                              pp(counter_ldaU+1,1), ldim1, ZERO, &
                              wfcatom(1,counter+1), npwx )
                  counter_ldaU = counter_ldaU + 2 * l + 1
              END IF
              counter = counter + 2 * l + 1
           ENDIF
        ENDDO
     ENDDO
     IF ( gamma_only ) THEN 
        CALL pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, wfcatom, npwx, rbecp, nkb)   
     ELSE 
        CALL ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom) 
     END IF 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 

     CALL davcio (swfcatom, nwordatwfc, iunat, ik, 1)


     ! on k-points 
  ENDDO 
  ! 
  IF ( gamma_only ) THEN 
     DEALLOCATE (rbecp) 
  ELSE
     DEALLOCATE (becp)
   END IF
  ! 

  DEALLOCATE (wfcatom) 
 
  DEALLOCATE (proj) 
  RETURN 
END SUBROUTINE projection
