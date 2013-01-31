!--------------------------------------------------------------
!OBM This subroutine initialises stuff related to open shell
! calculations (kpoint > 1 degauss/=0 or nspin/=1)
!-------------------------------------------------------------
SUBROUTINE lr_init_nfo()
!
!Created by Osman Baris Malcioglu (2009)
!
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks,degauss,lgauss,ngauss,xk, nelec
  USE wvfct,                ONLY : nbnd, et, igk, npw, g2kin
  USE realus,               ONLY : npw_k, igk_k
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : pi, degspin
  USE noncollin_module,     ONLY : noncolin
  USE mp,                   ONLY : mp_max, mp_min
  USE mp_global,            ONLY : inter_pool_comm
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : bg, tpiba, tpiba2, omega
  USE ener,                 ONLY : Ef
  USE ktetra,               ONLY : ltetra
  USE lsda_mod,             ONLY : lsda
  USE realus,               ONLY : real_space
  USE control_ph,           ONLY : alpha_pv, nbnd_occ
  USE wvfct,                ONLY : npwx, ecutwfc
  USE klist,                ONLY : nks
  USE io_files,             ONLY : iunigk, seqopn
  !
  IMPLICIT NONE
  !
  ! local variables
  REAL(kind=DP) :: small, emin, emax, xmax, fac, targ
  INTEGER       :: ik,ibnd, ipol
  LOGICAL       :: exst
  !
  ! Open shell related
  IF ( .not. allocated( igk_k ) )    ALLOCATE(igk_k(npwx,nks))
  IF ( .not. allocated( npw_k ) )    ALLOCATE(npw_k(nks))
  !IF ( .not. ALLOCATED( nbnd_occ ) ) allocate (nbnd_occ (nks))

  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )

  IF (.not. real_space) THEN
  DO ik=1,nks
      !
      CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
      !
      npw_k(ik) = npw
      !
      igk_k(:,ik) = igk(:)
      !
      ! For systems with more than one kpoint, we also write 
      ! igk to iunigk. This is required by exx_init().
      IF ( nks > 1 ) WRITE( iunigk ) igk
      !
  ENDDO
  ENDIF
  !OBM!! The following part is derived from phonon phq_setup
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  !if (.not. allocated (nbnd_occ) allocate( nbnd_occ (nks) )
  IF (lgauss) THEN
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss== - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     targ = ef + xmax * degauss
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) <targ) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) ==nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
  ELSEIF (ltetra) THEN
     CALL errore('lr_init_nfo','phonon + tetrahedra not implemented', 1)
  ELSE
     IF (lsda) CALL infomsg('lr_init_nfo','occupation numbers probably wrong')
     IF (noncolin) THEN
        nbnd_occ = nint (nelec)
     ELSE
        DO ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / degspin
        ENDDO
     ENDIF
  ENDIF
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
#ifdef __MPI
  ! find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
#endif
  IF (lgauss) THEN
     emax = targ
     alpha_pv = emax - emin
  ELSE
     emax = et (1, 1)
     DO ik = 1, nks
        DO ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        ENDDO
     ENDDO
#ifdef __MPI
     ! find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  ENDIF
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
RETURN
END SUBROUTINE lr_init_nfo
