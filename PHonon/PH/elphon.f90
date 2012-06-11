!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE elphon()
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : celldm, omega, ibrav
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, amass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE noncollin_module, ONLY : nspin_mag
  USE dynmat, ONLY : dyn, w2
  USE qpoint, ONLY : xq
  USE modes,  ONLY : npert, nirr, u
  USE control_ph, ONLY : trans
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  USE dfile_star,    ONLY : dvscf_star
  USE io_global, ONLY: stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: irr, imode0, ipert, is
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)


  CALL start_clock ('elphon')

  if(dvscf_star%basis.eq.'cartesian') then
     write(stdout,*) 'Setting patterns to identity'
     u=CMPLX(0.d0,0.d0)
     do irr=1,3*nat
        u(irr,irr)=CMPLX(1.d0,0.d0)
     enddo
  endif
  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin (dfftp%nnr, nspin_mag , npert(irr)) )
     DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                           imode0 + ipert,  -1 )
     END DO
     IF (doublegrid) THEN
        ALLOCATE (dvscfins (dffts%nnr, nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     CALL newdq (dvscfin, npert(irr))
     CALL elphel (npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  IF (.NOT.trans) CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
       ityp, omega, amass, tau, xq, w2, dyn)
  !
  CALL stop_clock ('elphon')
  RETURN
END SUBROUTINE elphon
!
!-----------------------------------------------------------------------
SUBROUTINE readmat (iudyn, ibrav, celldm, nat, ntyp, ityp, omega, &
     amass, tau, q, w2, dyn)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : amconv
  IMPLICIT NONE
  ! Input
  INTEGER :: iudyn, ibrav, nat, ntyp, ityp (nat)
  REAL(DP) :: celldm (6), amass (ntyp), tau (3, nat), q (3), &
       omega
  ! output
  REAL(DP) :: w2 (3 * nat)
  COMPLEX(DP) :: dyn (3 * nat, 3 * nat)
  ! local (control variables)
  INTEGER :: ntyp_, nat_, ibrav_, ityp_
  REAL(DP) :: celldm_ (6), amass_, tau_ (3), q_ (3)
  ! local
  REAL(DP) :: dynr (2, 3, nat, 3, nat)
  CHARACTER(len=80) :: line
  CHARACTER(len=3)  :: atm
  INTEGER :: nt, na, nb, naa, nbb, nu, mu, i, j
  !
  !
  REWIND (iudyn)
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  IF ( ntyp.NE.ntyp_ .OR. nat.NE.nat_ .OR.ibrav_.NE.ibrav .OR. &
       ABS ( celldm_ (1) - celldm (1) ) > 1.0d-5) &
          CALL errore ('readmat', 'inconsistent data', 1)
  DO nt = 1, ntyp
     READ (iudyn, * ) i, atm, amass_
     IF ( nt.NE.i .OR. ABS (amass_ - amconv*amass (nt) ) > 1.0d-5) &
        CALL errore ( 'readmat', 'inconsistent data', 1 + nt)
  ENDDO
  DO na = 1, nat
     READ (iudyn, * ) i, ityp_, tau_
     IF (na.NE.i.OR.ityp_.NE.ityp (na) ) CALL errore ('readmat', &
          'inconsistent data', 10 + na)
  ENDDO
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (line (11:80), * ) (q_ (i), i = 1, 3)
  READ (iudyn, '(a)') line
  DO na = 1, nat
     DO nb = 1, nat
        READ (iudyn, * ) naa, nbb
        IF (na.NE.naa.OR.nb.NE.nbb) CALL errore ('readmat', 'error reading &
             &file', nb)
        READ (iudyn, * ) ( (dynr (1, i, na, j, nb), dynr (2, i, na, j, nb) &
             , j = 1, 3), i = 1, 3)
     ENDDO
  ENDDO
  !
  ! divide the dynamical matrix by the (input) masses (in amu)
  !
  DO nb = 1, nat
     DO j = 1, 3
        DO na = 1, nat
           DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) ) / amconv
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) ) / amconv
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! solve the eigenvalue problem.
  ! NOTA BENE: eigenvectors are overwritten on dyn
  !
  CALL cdiagh (3 * nat, dynr, 3 * nat, w2, dyn)
  !
  ! divide by sqrt(mass) to get displacements
  !
  DO nu = 1, 3 * nat
     DO mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu) = dyn (mu, nu) / SQRT ( amconv * amass (ityp (na) ) )
     ENDDO
  ENDDO
  !
  !
  RETURN
END SUBROUTINE readmat
!
!-----------------------------------------------------------------------
SUBROUTINE elphel (npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !         <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !      Original routine written by Francesco Mauri
  !
  USE kinds, ONLY : DP
  USE fft_base, ONLY : dffts
  USE wavefunctions_module,  ONLY: evc
  USE io_files, ONLY: iunigk
  USE klist, ONLY: xk
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct, ONLY: nbnd, npw, npwx, igk
  USE uspp, ONLY : vkb
  USE el_phon, ONLY : el_ph_mat
  USE modes, ONLY : u
  USE units_ph, ONLY : iubar, lrbar, lrwfc, iuwfc
  USE eqv,      ONLY : dvpsi, evq
  USE qpoint,   ONLY : igkq, npwq, nksq, ikks, ikqs
  USE control_ph, ONLY : trans, lgamma
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  IMPLICIT NONE
  !
  INTEGER :: npe, imode0
  COMPLEX(DP) :: dvscfins (dffts%nnr, nspin_mag, npe)
  ! LOCAL variables
  INTEGER :: nrec, ik, ikk, ikq, ipert, mode, ibnd, jbnd, ir, ig, &
       ios
  COMPLEX(DP) , ALLOCATABLE :: aux1 (:,:), elphmat (:,:,:)
  COMPLEX(DP), EXTERNAL :: zdotc
  !
  ALLOCATE (aux1    (dffts%nnr, npol))
  ALLOCATE (elphmat ( nbnd , nbnd , npe))
  !
  !  Start the loops over the k-points
  !
  IF (nksq.GT.1) REWIND (unit = iunigk)
  DO ik = 1, nksq
     IF (nksq.GT.1) THEN
        READ (iunigk, err = 100, iostat = ios) npw, igk
100     CALL errore ('elphel', 'reading igk', ABS (ios) )
     ENDIF
     !
     !  ik = counter of k-points with vector k
     !  ikk= index of k-point with vector k
     !  ikq= index of k-point with vector k+q
     !       k and k+q are alternated if q!=0, are the same if q=0
     !
     IF (lgamma) npwq = npw
     ikk = ikks(ik)
     ikq = ikqs(ik)
     IF (lsda) current_spin = isk (ikk)
     IF (.NOT.lgamma.AND.nksq.GT.1) THEN
        READ (iunigk, err = 200, iostat = ios) npwq, igkq
200     CALL errore ('elphel', 'reading igkq', ABS (ios) )
     ENDIF
     !
     CALL init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     !
     ! read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     IF (nksq.GT.1) THEN
        IF (lgamma) THEN
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
        ELSE
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
           CALL davcio (evq, lrwfc, iuwfc, ikq, - 1)
        ENDIF
     ENDIF
     !
     DO ipert = 1, npe
        nrec = (ipert - 1) * nksq + ik
        !
        !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
        !
        IF (trans) THEN
           CALL davcio (dvpsi, lrbar, iubar, nrec, - 1)
        ELSE
           mode = imode0 + ipert
           ! TODO : .false. or .true. ???
           CALL dvqpsi_us (ik, u (1, mode), .FALSE. )
        ENDIF
        !
        ! calculate dvscf_q*psi_k
        !
        DO ibnd = 1, nbnd
           CALL cft_wave (evc(1, ibnd), aux1, +1)
           CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipert), current_spin)
           CALL cft_wave (dvpsi(1, ibnd), aux1, -1)
        END DO
        CALL adddvscf (ipert, ik)

        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        DO ibnd =1, nbnd
           DO jbnd = 1, nbnd
              elphmat (jbnd, ibnd, ipert) = zdotc (npwq, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
              IF (noncolin) &
                 elphmat (jbnd, ibnd, ipert) = elphmat (jbnd, ibnd, ipert)+ &
                   zdotc (npwq, evq(npwx+1,jbnd),1,dvpsi(npwx+1,ibnd), 1)
           ENDDO
        ENDDO
     ENDDO
     !
     CALL mp_sum (elphmat, intra_pool_comm)
     !
     !  save all e-ph matrix elements into el_ph_mat
     !
     DO ipert = 1, npe
        DO jbnd = 1, nbnd
           DO ibnd = 1, nbnd
              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE (elphmat)
  DEALLOCATE (aux1)
  !
  RETURN
END SUBROUTINE elphel
!
!-----------------------------------------------------------------------
SUBROUTINE elphsum ( )
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri, modified by PG
  !      New version by  Malgorzata Wierzbowska
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : pi, rytoev, degspin
  USE ions_base,   ONLY : nat, ityp, tau
  USE cell_base,   ONLY : at, bg
  USE lsda_mod,    ONLY: isk, nspin
  USE klist,       ONLY: nks, nkstot, xk, wk, nelec
  USE start_k,     ONLY: nk1, nk2, nk3
  USE symm_base,   ONLY: s, irt, nsym, invs
  USE noncollin_module, ONLY: nspin_lsda, nspin_mag
  USE wvfct,       ONLY: nbnd, et
  USE parameters,  ONLY : npk
  USE el_phon,     ONLY : el_ph_mat
  USE qpoint,      ONLY : xq, nksq
  USE modes,       ONLY : u, minus_q, nsymq, rtau
  USE dynmat,      ONLY : dyn, w2
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp_global,   ONLY : my_pool_id, npool, kunit, intra_image_comm
  USE mp,          ONLY : mp_bcast
  USE control_ph,  ONLY : lgamma, tmp_dir_phq, xmldyn
  USE save_ph,     ONLY : tmp_dir_save
  USE io_files,    ONLY : prefix, tmp_dir, seqopn
  !
  IMPLICIT NONE
  ! epsw = 20 cm^-1, in Ry
  REAL(DP), PARAMETER :: Rytocm1 = 109737.57990d0, RytoGHz = 3.289828D6, &
       RytoTHz = RytoGHz/1000.d0, epsw = 20.d0 / Rytocm1, eps = 1.0d-6
  !
  INTEGER :: iuna2Fsave  = 40
  !
  REAL(DP), allocatable :: gam(:,:), lamb(:,:)
  !
  ! Quantities ending with "fit" are relative to the "dense" grid
  !
  REAL(DP), allocatable :: xkfit(:,:)
  REAL(DP), allocatable, target :: etfit(:,:), wkfit(:)
  INTEGER :: nksfit, nk1fit, nk2fit, nk3fit, nkfit, nksfit_real
  INTEGER, allocatable :: eqkfit(:), eqqfit(:), sfit(:)
  !
  integer :: nq, isq (48), imq
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)
  real(DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  ! workspace used for symmetrisation
  !
  COMPLEX(DP), allocatable :: g1(:,:,:), g2(:,:,:), g0(:,:), gf(:,:,:)
  COMPLEX(DP), allocatable :: point(:), noint(:), ctemp(:)
  COMPLEX(DP) :: dyn22(3*nat,3*nat)
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, nsig, iuelph, ios, i,k,j, ii, jj
  INTEGER :: nkBZ, nti, ntj, ntk, nkr, itemp1, itemp2, nn, &
       qx,qy,qz,iq,jq,kq
  INTEGER, ALLOCATABLE :: eqBZ(:), sBZ(:)
  REAL(DP) :: weight, wqa, w0g1, w0g2, degauss1, dosef, &
       ef1, lambda, gamma
  REAL(DP) :: deg(10), effit(10), dosfit(10), etk, etq
  REAL(DP), EXTERNAL :: dos_ef, efermig, w0gauss
  character(len=80) :: name
  LOGICAL  :: exst, xmldyn_save
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat)

  COMPLEX(DP), POINTER :: el_ph_mat_collect(:,:,:,:)
  REAL(DP), ALLOCATABLE :: xk_collect(:,:), wk_collect(:)
  REAL(DP), POINTER :: wkfit_dist(:), etfit_dist(:,:)
  INTEGER :: nksfit_dist, rest, kunit_save
  INTEGER :: nks_real, ispin, nksqtot
  !
  !
  WRITE (6, '(5x,"electron-phonon interaction  ..."/)')
  ngauss1 = 0
  nsig = 10

  ALLOCATE(xk_collect(3,nkstot))
  ALLOCATE(wk_collect(nkstot))
  IF (npool==1) THEN
!
!  no pool, just copy old variable on the new ones
!
     nksqtot=nksq
     xk_collect(:,1:nks) = xk(:,1:nks)
     wk_collect(1:nks) = wk(1:nks)
     el_ph_mat_collect => el_ph_mat
  ELSE  
!
!  pools, allocate new variables and collect the results. All the rest
!  remain unchanged.
!
     IF (lgamma) THEN
        nksqtot=nkstot
     ELSE
        nksqtot=nkstot/2
     ENDIF
     ALLOCATE(el_ph_mat_collect(nbnd,nbnd,nksqtot,3*nat))
     CALL xk_wk_collect(xk_collect,wk_collect,xk,wk,nkstot,nks)
     CALL el_ph_collect(el_ph_mat,el_ph_mat_collect,nksqtot,nksq)
  ENDIF
  !
  ! read eigenvalues for the dense grid
  ! FIXME: this might be done from the xml file, not from a specialized file
  ! parallel case: only first node reads
  !
  IF ( ionode ) THEN
     tmp_dir=tmp_dir_save
     CALL seqopn( iuna2Fsave, 'a2Fsave', 'FORMATTED', exst )
     tmp_dir=tmp_dir_phq
     READ(iuna2Fsave,*) ibnd, nksfit
  END IF
  !
  CALL mp_bcast (ibnd, ionode_id, intra_image_comm)
  CALL mp_bcast (nksfit, ionode_id, intra_image_comm)
  if ( ibnd /= nbnd ) call errore('elphsum','wrong file read',iuna2Fsave)
  allocate (etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
  !
  IF ( ionode ) THEN
     READ(iuna2Fsave,*) etfit
     READ(iuna2Fsave,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
     READ(iuna2Fsave,*) wkfit
     READ(iuna2Fsave,*) nk1fit, nk2fit, nk3fit
     CLOSE( UNIT = iuna2Fsave, STATUS = 'KEEP' )
  END IF
  !
  ! broadcast all variables read
  !
  CALL mp_bcast (etfit, ionode_id, intra_image_comm)
  CALL mp_bcast (xkfit, ionode_id, intra_image_comm)
  CALL mp_bcast (wkfit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk1fit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk2fit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk3fit, ionode_id, intra_image_comm)
  !
  nkfit=nk1fit*nk2fit*nk3fit
  !
  ! efermig and dos_ef require scattered points and eigenvalues
  ! isk is neither read nor used. phonon with two Fermi energies is
  ! not yet implemented.
  !
  nksfit_dist  = ( nksfit / npool )
  rest = ( nksfit - nksfit_dist * npool ) 
  IF ( ( my_pool_id + 1 ) <= rest ) nksfit_dist = nksfit_dist + 1
  kunit_save=kunit
  kunit=1

#ifdef __MPI
  ALLOCATE(etfit_dist(nbnd,nksfit_dist))
  ALLOCATE(wkfit_dist(nksfit_dist))
  CALL poolscatter( 1, nksfit, wkfit, nksfit_dist, wkfit_dist )
  CALL poolscatter( nbnd, nksfit, etfit, nksfit_dist, etfit_dist )
#else
   wkfit_dist => wkfit
   etfit_dist => etfit
#endif
  !
  do isig=1,nsig
     !
     ! recalculate Ef = effit and DOS at Ef N(Ef) = dosfit using dense grid
     ! for value "deg" of gaussian broadening
     !
     deg(isig) = isig * 0.005d0
     !
     effit(isig) = efermig &
          ( etfit_dist, nbnd, nksfit_dist, nelec, wkfit_dist, &
              deg(isig), ngauss1, 0, isk)
     dosfit(isig) = dos_ef ( ngauss1, deg(isig), effit(isig), etfit_dist, &
          wkfit_dist, nksfit_dist, nbnd) / 2.0d0
  enddo
#ifdef __MPI
  DEALLOCATE(etfit_dist)
  DEALLOCATE(wkfit_dist)
#endif
  kunit=kunit_save
  allocate (eqkfit(nkfit), eqqfit(nkfit), sfit(nkfit))
  !
  ! map k-points in the IBZ to k-points in the complete uniform grid
  !
  nksfit_real=nksfit/nspin_lsda
  call lint ( nsym, s, .true., at, bg, npk, 0,0,0, &
       nk1fit,nk2fit,nk3fit, nksfit_real, xkfit, 1, nkfit, eqkfit, sfit)
  deallocate (sfit, xkfit, wkfit)
  !
  ! find epsilon(k+q) in the dense grid
  !
  call cryst_to_cart (1, xq, at, -1)
  qx = nint(nk1fit*xq(1))
  if (abs(qx-nk1fit*xq(1)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',1)
  if (qx < 0) qx = qx + nk1fit
  if (qx > nk1fit) qx = qx - nk1fit
  qy = nint(nk2fit*xq(2))
  if (abs(qy-nk2fit*xq(2)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',2)
  if (qy < 0) qy = qy + nk2fit
  if (qy > nk2fit) qy = qy - nk2fit
  qz = nint(nk3fit*xq(3))
  if (abs(qz-nk3fit*xq(3)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',3)
  if (qz < 0) qz = qz + nk3fit
  if (qz > nk3fit) qz = qz - nk3fit
  call cryst_to_cart (1, xq, bg, 1)
  !
  eqqfit(:) = 0
  do i=1,nk1fit
     do j=1,nk2fit
        do k=1,nk3fit
           ik = k-1 + (j-1)*nk3fit + (i-1)*nk2fit*nk3fit + 1
           iq = i+qx
           if (iq > nk1fit) iq = iq - nk1fit
           jq = j+qy
           if (jq > nk2fit) jq = jq - nk2fit
           kq = k+qz
           if (kq > nk3fit) kq = kq - nk3fit
           nn = (kq-1)+(jq-1)*nk3fit+(iq-1)*nk2fit*nk3fit + 1
           eqqfit(ik) = eqkfit(nn)
        enddo
     enddo
  enddo
  !
  ! calculate the electron-phonon coefficient using the dense grid
  !
  nti  = nk1fit/nk1
  ntj  = nk2fit/nk2
  ntk  = nk3fit/nk3
  nkBZ  = nk1*nk2*nk3
  allocate (eqBZ(nkBZ), sBZ(nkBZ))
  !
  nks_real=nkstot/nspin_lsda
  IF ( lgamma ) THEN
     call lint ( nsymq, s, minus_q, at, bg, npk, 0,0,0, &
          nk1,nk2,nk3, nks_real, xk_collect, 1, nkBZ, eqBZ, sBZ)
  ELSE
     call lint ( nsymq, s, minus_q, at, bg, npk, 0,0,0, &
          nk1,nk2,nk3, nks_real, xk_collect, 2, nkBZ, eqBZ, sBZ)
  END IF
  !
  allocate (gf(3*nat,3*nat,nsig))
  gf = (0.0d0,0.0d0)
  !
  wqa  = 1.0d0/nkfit
  IF (nspin==1) wqa=degspin*wqa
  !
  do ibnd = 1, nbnd
     do jbnd = 1, nbnd
        allocate (g2(nkBZ*nspin_lsda,3*nat,3*nat))
        allocate (g1(nksqtot,3*nat,3*nat))
        do ik = 1, nksqtot
           do ii = 1, 3*nat
              do jj = 1, 3*nat
                 g1(ik,ii,jj)=CONJG(el_ph_mat_collect(jbnd,ibnd,ik,ii))* &
                      el_ph_mat_collect(jbnd,ibnd,ik,jj)
              enddo    ! ipert
           enddo    !jpert
        enddo   ! ik
        !
        allocate (g0(3*nat,3*nat))
        do i=1,nk1
           do j=1,nk2
              do k=1,nk3
                 do ispin=1,nspin_lsda
                    nn = k-1 + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                    itemp1 = eqBZ(nn)
                    if (ispin==2) itemp1=itemp1+nksqtot/2
                    g0(:,:) = g1(itemp1,:,:)
                    itemp2 = sBZ(nn)
                    call symm ( g0, u, xq, s, itemp2, rtau, irt, &
                         at, bg, nat)
                    if (ispin==2) nn=nn+nkBZ 
                    g2(nn,:,:) = g0(:,:)
                 enddo
              enddo ! k
           enddo !j
        enddo !i
        deallocate (g0)
        deallocate (g1)
        !
        allocate ( point(nkBZ), noint(nkfit), ctemp(nkfit) )
        do jpert = 1, 3 * nat
           do ipert = 1, 3 * nat
              do ispin=1,nspin_lsda
              !
              point(1:nkBZ) = &
                  g2(1+nkBZ*(ispin-1):nkBZ+nkBZ*(ispin-1),ipert,jpert)
              !
              CALL clinear(nk1,nk2,nk3,nti,ntj,ntk,point,noint)
              !
              do isig = 1, nsig
                 degauss1 = deg(isig)
                 do ik=1,nkfit
                    etk = etfit(ibnd,eqkfit(ik)+nksfit*(ispin-1)/2)
                    etq = etfit(jbnd,eqqfit(ik)+nksfit*(ispin-1)/2)
                    w0g1 = w0gauss( (effit(isig)-etk) &
                                   / degauss1,ngauss1) / degauss1
                    w0g2 = w0gauss( (effit(isig)-etq) &
                                   / degauss1,ngauss1) / degauss1
                    ctemp(ik) = noint(ik)* wqa * w0g1 * w0g2
                 enddo
                 gf(ipert,jpert,isig) = gf(ipert,jpert,isig) + &
                      SUM (ctemp)
              enddo ! isig
              enddo ! ispin
           enddo    ! ipert
        enddo    !jpert
        deallocate (point, noint, ctemp)
        deallocate (g2)
        !
     enddo    ! ibnd
  enddo    ! jbnd

  deallocate (eqqfit, eqkfit)
  deallocate (etfit)
  deallocate (eqBZ, sBZ)
!
  allocate (gam(3*nat,nsig), lamb(3*nat,nsig))
  lamb(:,:) = 0.0d0
  gam (:,:) = 0.0d0
  do isig= 1,nsig
     do nu = 1,3*nat
        gam(nu,isig) = 0.0d0
        do mu = 1, 3 * nat
           do vu = 1, 3 * nat
              gam(nu,isig) = gam(nu,isig) + DBLE(conjg(dyn(mu,nu)) * &
                   gf(mu,vu,isig) * dyn(vu,nu))
           enddo
        enddo
        gam(nu,isig) = gam(nu,isig) *  pi/2.0d0
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        if (sqrt(abs(w2(nu))) > epsw) then
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lamb(nu,isig) = gam(nu,isig)/pi/w2(nu)/dosfit(isig)
        else
           lamb(nu,isig) = 0.0d0
        endif
        gam(nu,isig) = gam(nu,isig)*RytoGHz
     enddo  !nu
  enddo  ! isig
  !
  do isig= 1,nsig
     WRITE (6, 9000) deg(isig), ngauss1
     WRITE (6, 9005) dosfit(isig), effit(isig) * rytoev
     do nu=1,3*nat
        WRITE (6, 9010) nu, lamb(nu,isig), gam(nu,isig)
     enddo
  enddo
  ! Isaev: save files in suitable format for processing by lambda.x
  write(name,'(A5,f9.6,A1,f9.6,A1,f9.6)') 'elph.',xq(1),'.',xq(2),'.',xq(3)
  open(12,file=name, form='formatted', status='unknown')

  write(12, "(5x,3f14.6,2i6)") xq(1),xq(2),xq(3), nsig, 3*nat
  write(12, "(6e14.6)") (w2(i), i=1,3*nat)

  do isig= 1,nsig
     WRITE (12, 9000) deg(isig), ngauss1
     WRITE (12, 9005) dosfit(isig), effit(isig) * rytoev
     do nu=1,3*nat
        WRITE (12, 9010) nu, lamb(nu,isig), gam(nu,isig)
     enddo
  enddo
  close (unit=12,status='keep')
  ! Isaev end
  deallocate (gam)
  deallocate (lamb)
  write(stdout,*)
  !
  !    Prepare interface to q2r and matdyn
  !
  call star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE. )
  !
  do isig=1,nsig
     write(name,"(A7,I2)") 'a2Fq2r.',50 + isig
     if (ionode) then
        iuelph = 4
        open(iuelph, file=name, STATUS = 'unknown', FORM = 'formatted', &
             POSITION='append')
     else
        !
        ! this node doesn't write: unit 6 is redirected to /dev/null
        !
        iuelph =6
     end if
     dyn22(:,:) = gf(:,:,isig)
     write(iuelph,*) deg(isig), effit(isig), dosfit(isig)
     IF ( imq == 0 ) THEN
        write(iuelph,*) 2*nq
     ELSE
        write(iuelph,*) nq
     ENDIF
     xmldyn_save=xmldyn
     xmldyn=.FALSE.
     call q2qstar_ph (dyn22, at, bg, nat, nsym, s, invs, &
          irt, rtau, nq, sxq, isq, imq, iuelph)
     xmldyn=xmldyn_save
     if (ionode) CLOSE( UNIT = iuelph, STATUS = 'KEEP' )
  enddo
  deallocate (gf)
  DEALLOCATE(xk_collect)
  DEALLOCATE(wk_collect)
  IF (npool /= 1) DEALLOCATE(el_ph_mat_collect)

  !
9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i2,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  RETURN
END SUBROUTINE elphsum

!-----------------------------------------------------------------------
SUBROUTINE elphsum_simple
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !      Rewritten by Matteo Calandra
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, ry_to_cmm1, rytoev
  USE ions_base, ONLY : nat, ityp, tau,amass,tau, ntyp => nsp, atm
  USE cell_base, ONLY : at, bg, ibrav, celldm 
  USE fft_base,  ONLY: dfftp
  USE symm_base, ONLY : s, sr, irt, nsym, time_reversal, invs
  USE klist, ONLY : xk, nelec, nks, wk
  USE wvfct, ONLY : nbnd, et
  USE el_phon
  USE mp_global, ONLY : me_pool, root_pool, inter_pool_comm, npool, intra_pool_comm
  USE io_global, ONLY : stdout
  USE klist, only : degauss,ngauss
  USE control_flags, ONLY : modenum, noinv
  USE units_ph,       ONLY :iudyn
  USE io_files,  ONLY : prefix
  USE qpoint, ONLY : xq, nksq
  USE dynmat, ONLY : dyn, w2
  USE modes, ONLY : u,rtau, irgq, nsymq,irotmq, minus_q
  USE control_ph, only : lgamma
  USE lsda_mod, only : isk,nspin, current_spin,lsda
  USE mp,        ONLY: mp_sum
  !
  IMPLICIT NONE
  REAL(DP), PARAMETER :: eps = 20_dp/ry_to_cmm1 ! eps = 20 cm^-1, in Ry
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, nsig, iuelph, ios, iuelphmat,icnt,i,j,rrho,nt,k
  INTEGER :: na,nb,icar,jcar,iu_sym,nmodes
  INTEGER :: iu_Delta_dyn,iu_analdyn,iu_nonanaldyn
  INTEGER :: io_file_unit
  !   for star_q
  INTEGER :: nsymloc, sloc(3,3,48), invsloc(48), irtloc(48,nat), &
             nqloc, isqloc(48), imqloc
  REAL(DP) :: rtauloc(3,48,nat), sxqloc(3,48)
  !   end of star_q definitions
  REAL(DP) :: weight, w0g1, w0g2, w0gauss, wgauss,degauss1, dosef, &
       ef1, phase_space, lambda, gamma, wg1, w0g,wgp,deltae
  REAL(DP), EXTERNAL :: dos_ef, efermig
  REAL(DP) xk_dummy(3)
  COMPLEX(DP), allocatable :: phi(:,:,:,:),phi_nonanal(:,:,:,:)
  COMPLEX(DP), allocatable :: dyn_mat_r(:,:),zz(:,:)
  CHARACTER(len=20) :: char_deg
  CHARACTER(len=1) :: char_ng
  character(len=80) :: filelph
  CHARACTER(len=256) ::  file_elphmat
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit

  nmodes=3*nat


  write(filelph,'(A5,f9.6,A1,f9.6,A1,f9.6)') 'elph.',xq(1),'.',xq(2),'.',xq(3)

  ! parallel case: only first node writes
  IF ( me_pool /= root_pool ) THEN
     iuelph = 0
  ELSE
     !
     iuelph = find_free_unit()
     OPEN (unit = iuelph, file = filelph, status = 'unknown', err = &
          100, iostat = ios)
100  CALL errore ('elphon', 'opening file '//filelph, ABS (ios) )
     REWIND (iuelph)
     !
  END IF

  WRITE (iuelph, '(3f15.8,2i8)') xq, nsig, 3 * nat
  WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
  

  ngauss1=0
  DO isig = 1, el_ph_nsigma
     !     degauss1 = 0.01 * isig
     degauss1 = el_ph_sigma * isig
     write(stdout,*) degauss1
     el_ph_sum(:,:) = (0.d0, 0.d0)
     phase_space = 0.d0
     !
     ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
     ! for this gaussian broadening
     !
     ! Note that the weights of k+q points must be set to zero for the
     ! following call to yield correct results
     !
      
     ef1 = efermig (et, nbnd, nks, nelec, wk, degauss1, el_ph_ngauss, 0, isk)
     dosef = dos_ef (el_ph_ngauss, degauss1, ef1, et, wk, nks, nbnd)
     ! N(Ef) is the DOS per spin, not summed over spin
     dosef = dosef / 2.d0
     !
     ! Sum over bands with gaussian weights
     !
     
     DO ik = 1, nksq
        
        !
        ! see subroutine elphel for the logic of indices
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        DO ibnd = 1, nbnd
           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
                / degauss1
           xk_dummy(:)=xk(:,ikk)
           call cryst_to_cart(1,xk_dummy,at,-1)
           DO jbnd = 1, nbnd
              w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                   / degauss1
              ! note that wk(ikq)=wk(ikk)
              weight = wk (ikk) * w0g1 * w0g2
              DO jpert = 1, 3 * nat
                 DO ipert = 1, 3 * nat
                    el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert)  +  weight * &
                         CONJG (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                         el_ph_mat (jbnd, ibnd, ik, jpert)
                 ENDDO
              ENDDO
              phase_space = phase_space+weight
           ENDDO
        ENDDO
        
     ENDDO
     
     ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
     !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
     !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
     !
     ! collect contributions from all pools (sum over k-points)
     !
     
!     CALL poolreduce (2 * 3 * nat * 3 * nat, el_ph_sum)
!     CALL poolreduce (1, phase_space)
     call mp_sum ( el_ph_sum , inter_pool_comm )
     call mp_sum ( phase_space , inter_pool_comm )

     !
     ! symmetrize el_ph_sum(mu,nu) : it transforms as the dynamical matrix
     !
     
     CALL symdyn_munu (el_ph_sum, u, xq, s, invs, rtau, irt, irgq, at, &
          bg, nsymq, nat, irotmq, minus_q)
     !
     WRITE (6, 9000) degauss1, ngauss1
     WRITE (6, 9005) dosef, ef1 * rytoev
     WRITE (6, 9006) phase_space
     IF (iuelph.NE.0) THEN
        WRITE (iuelph, 9000) degauss1, ngauss1
        WRITE (iuelph, 9005) dosef, ef1 * rytoev
     ENDIF
     
     DO nu = 1, nmodes
        gamma = 0.0
        DO mu = 1, 3 * nat
           DO vu = 1, 3 * nat
              gamma = gamma + DBLE (CONJG (dyn (mu, nu) ) * el_ph_sum (mu, vu)&
                   * dyn (vu, nu) )
           ENDDO
        ENDDO
        write(819+mu,*) gamma
        gamma = pi * gamma / 2.d0
        write(6,*) 'gamma*pi/2=',gamma
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        IF (SQRT (ABS (w2 (nu) ) ) > eps) THEN
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lambda = gamma / pi / w2 (nu) / dosef
        ELSE
           lambda = 0.0
        ENDIF
        ! 3.289828x10^6 is the conversion factor from Ry to GHz
        WRITE (6, 9010) nu, lambda, gamma * 3.289828d6
        IF (iuelph.NE.0) WRITE (iuelph, 9010) nu, lambda, gamma * &
             3.289828d6
     ENDDO
  ENDDO
  

9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i2,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  !
  IF (iuelph.NE.0) CLOSE (unit = iuelph)
  RETURN
  

     
     !          call star_q(x_q(1,iq), at, bg, nsym , s , invs , nq, sxq, &
     !               isq, imq, .FALSE. )
     

END SUBROUTINE elphsum_simple
   



!-----------------------------------------------------------------------
FUNCTION dos_ef (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  IMPLICIT NONE
  REAL(DP) :: dos_ef
  INTEGER :: ngauss, nbnd, nks
  REAL(DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  INTEGER :: ik, ibnd
  REAL(DP), EXTERNAL :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef = 0.0d0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef = dos_ef + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  !    Collects partial sums on k-points from all pools
  !
  CALL mp_sum ( dos_ef, inter_pool_comm )
  !
  RETURN
END FUNCTION dos_ef

!a2F
subroutine lint ( nsym, s, minus_q, at, bg, npk, k1,k2,k3, &
     nk1,nk2,nk3, nks, xk, kunit, nkBZ, eqBZ, sBZ)
  !-----------------------------------------------------------------------
  !
  ! Find which k-points of a uniform grid are in the IBZ
  !
  use kinds, only : DP
  implicit none
  integer, intent (IN) :: nks, nsym, s(3,3,48), npk, k1, k2, k3, &
       nk1, nk2, nk3, kunit, nkBZ
  logical, intent (IN) :: minus_q
  real(kind=DP), intent(IN):: at(3,3), bg(3,3), xk(3,npk)
  integer, INTENT(OUT) :: eqBZ(nkBZ), sBZ(nkBZ)
  !
  real(kind=DP) :: xkr(3), deltap(3), deltam(3)
  real(kind=DP), parameter:: eps=1.0d-5
  real(kind=DP), allocatable :: xkg(:,:), xp(:,:)
  integer ::  i,j,k, ns, n, nk
  integer :: nkh
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  allocate (xkg( 3,nkBZ))
  !
  if(kunit < 1 .or. kunit > 2) call errore('lint','bad kunit value',kunit)
  !
  ! kunit=2: get only "true" k points, not k+q points, from the list
  !
  nkh = nks/kunit
  allocate (xp(3,nkh))
  if (kunit == 1) then
     xp(:,1:nkh) = xk(:,1:nkh)
  else
     do j=1,nkh
        xp(:,j) = xk(:,2*j-1)
     enddo
  end if
  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        end do
     end do
  end do

  call cryst_to_cart (nkh,xp,at,-1)

  do nk=1,nkBZ
     do n=1,nkh
        do ns=1,nsym
           do i=1,3
              xkr(i) = s(i,1,ns) * xp(1,n) + &
                       s(i,2,ns) * xp(2,n) + &
                       s(i,3,ns) * xp(3,n)
           end do
           do i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           end do
           if ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( minus_q .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) then
              eqBZ(nk) = n
              sBZ(nk) = ns
              go to 15
           end if
        end do
     end do
     call errore('lint','cannot locate  k point  xk',nk)
15   continue
  end do

  do n=1,nkh
     do nk=1,nkBZ
        if (eqBZ(nk) == n) go to 20
     end do
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     call errore('lint','cannot remap grid on k-point list',n)
20   continue
  end do

  deallocate(xkg)
  deallocate(xp)

  return
end subroutine lint
