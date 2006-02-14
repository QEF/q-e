!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE gvect, ONLY: nrxx
  USE gsmooth, ONLY: nrxxs, doublegrid
  USE lsda_mod, ONLY: nspin
  USE phcom
  USE el_phon
  !
  IMPLICIT NONE
  !
  INTEGER :: irr, imode0, ipert, is
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)


  CALL start_clock ('elphon')

  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  !!! rewind (iudvscf)
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin ( nrxx , nspin , npert(irr)) )
     !!! read (iudvscf) dvscfin
     DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                           imode0 + ipert, -1 )
     END DO
     IF (doublegrid) THEN
        ALLOCATE (dvscfins ( nrxxs , nspin , npert(irr)) )
        DO is = 1, nspin
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
  IF (ntyp.NE.ntyp_.OR.nat.NE.nat_.OR.ibrav_.NE.ibrav.OR.ABS ( &
       celldm_ (1) - celldm (1) ) .GT.1.0d-5) CALL errore ('readmat', &
       'inconsistent data', 1)
  DO nt = 1, ntyp
     READ (iudyn, * ) i, atm, amass_
     IF (nt.NE.i.OR.ABS (amass_ - amass (nt) ) .GT.1.0d-5) CALL errore ( &
          'readmat', 'inconsistent data', 1 + nt)
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
  ! divide the dynamical matrix by the masses
  !
  DO nb = 1, nat
     DO j = 1, 3
        DO na = 1, nat
           DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) )
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
        dyn (mu, nu) = dyn (mu, nu) / SQRT (amass (ityp (na) ) )
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
  USE gsmooth, ONLY: nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls
  USE wavefunctions_module,  ONLY: evc
  USE io_files, ONLY: iunigk
  USE klist, ONLY: xk
  USE lsda_mod, ONLY: nspin, lsda, current_spin, isk
  USE wvfct, ONLY: nbnd, npw, igk
  USE uspp, ONLY : vkb
  USE phcom
  USE el_phon
  IMPLICIT NONE
  !
  INTEGER :: npe, imode0
  COMPLEX(DP) :: dvscfins (nrxxs, nspin, npe)
  ! LOCAL variables
  INTEGER :: nrec, ik, ikk, ikq, ipert, mode, ibnd, jbnd, ir, ig, &
       ios
  COMPLEX(DP) , ALLOCATABLE :: aux1 (:), elphmat (:,:,:)
  COMPLEX(DP) :: ZDOTC
  !
  ALLOCATE (aux1    ( nrxxs))    
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
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
        npwq = npw
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
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
           CALL dvqpsi_us (ik, mode, u (1, mode), .FALSE. )
        ENDIF
        !
        ! calculate dvscf_q*psi_k
        !
        DO ibnd = 1, nbnd
           aux1(:) = (0.d0, 0.d0)
           DO ig = 1, npw
              aux1 (nls (igk (ig) ) ) = evc (ig, ibnd)
           ENDDO
           CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
           DO ir = 1, nrxxs
              aux1 (ir) = aux1 (ir) * dvscfins (ir, current_spin, ipert)
           ENDDO
           CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
           DO ig = 1, npwq
              dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux1 (nls (igkq (ig) ) )
           ENDDO
        END DO
        CALL adddvscf (ipert, ik)

        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        DO ibnd =1, nbnd
           DO jbnd = 1, nbnd
              elphmat (jbnd, ibnd, ipert) = ZDOTC (npwq, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
           ENDDO
           !
        ENDDO
     ENDDO
     !
     CALL reduce (2 * nbnd * nbnd * npe, elphmat)
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
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi, rytoev
  USE ions_base,     ONLY : nat, ityp, tau
  USE cell_base,     ONLY : at, bg, ibrav, symm_type
  USE gvect, ONLY: nr1, nr2, nr3
  USE lsda_mod, ONLY: isk
  USE klist, ONLY: nks, xk, wk, nelec
  USE ktetra, ONLY: nk1, nk2, nk3
  USE symme, ONLY: s, irt, nsym
  USE wvfct, ONLY: nbnd, et
  USE phcom
  USE el_phon
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp_global, ONLY : npool
  USE mp, ONLY : mp_bcast
  USE control_flags, ONLY : modenum, noinv
  USE control_ph, ONLY : lgamma
  USE io_files,  ONLY : prefix
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
  REAL(DP), allocatable :: etfit(:,:), xkfit(:,:), wkfit(:)
  INTEGER :: nksfit, nk1fit, nk2fit, nk3fit, nkfit
  INTEGER, allocatable :: eqkfit(:), eqqfit(:), sfit(:)
  !
  ! Quantities ending with "gam" are symmetries for Gamma (q=0) 
  !
  INTEGER :: nsymgam, sgam(3,3,48), invsgam(48), tablegam(48,48), &
       irtgam(48,nat)
  REAL(DP) :: rtaugam(3,48,nat)
  LOGICAL  :: symgam(48)
  !
  ! workspace used for symmetrisation
  !
  COMPLEX(DP), allocatable :: g0(:,:,:), g3(:,:,:), g1(:,:), g2(:,:), gf(:,:,:)
  COMPLEX(DP), allocatable :: point(:), noint(:), ctemp(:)
  COMPLEX(DP) :: dyn22(3*nat,3*nat)
  !
  ! Quantities ending with "loc" are ... what???
  !
  integer :: nsymloc, sloc(3,3,48), invsloc(48), irtloc(48,nat), &
             nqloc, isqloc(48), imqloc
  real(DP) :: rtauloc(3,48,nat), sxqloc(3,48)
  !
  ! Misc auxiliary variables
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, ngauss2, nsig, iuelph, ios, i,k,j, ii, jj
  INTEGER :: nkBZ, nti, ntj, ntk, nkr, itemp1, itemp2, nn, &
       qx,qy,qz,iq,jq,kq
  INTEGER, ALLOCATABLE :: eqBZ(:), sBZ(:)
  REAL(DP) :: weight, wqa, w0g1, w0g2, degauss1, degauss2, dosef, &
       ef1, phase_space, lambda, gamma
  REAL(DP) :: deg(10), effit(10), dosfit(10), etk, etq
  REAL(DP), EXTERNAL :: dos_ef, efermig, w0gauss
  character(len=9) :: name
  LOGICAL  :: exst
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat)
  !
  !
  WRITE (6, '(5x,"electron-phonon interaction  ..."/)')
  ngauss1 = 1
  ngauss2 = 1
  nsig = 10
  !
  IF(la2F) THEN
     !
     IF (npool > 1) CALL errore ('elphsum', 'pools and a2F not implemented', 1)
     !
     ! read eigenvalues for the dense grid
     ! parallel case: only first node reads
     !
     IF ( ionode ) THEN
        CALL seqopn( iuna2Fsave, 'a2Fsave', 'FORMATTED', exst )
        READ(iuna2Fsave,*) ibnd, nksfit
     END IF
     !
     CALL mp_bcast (ibnd, ionode_id)
     CALL mp_bcast (nksfit, ionode_id)
     if ( ibnd /= nbnd ) call errore('elphsum','wrong file read',iuna2Fsave)
     allocate (etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
     !
     IF ( ionode ) THEN
        READ(iuna2Fsave,*) etfit
        READ(iuna2Fsave,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
        READ(iuna2Fsave,*) wkfit
        READ(iuna2Fsave,*) nk1fit, nk2fit, nk3fit
        ! 
        !    Read symmetries for q=0 from file (needed for symmetrization)
        !
        READ( iuna2Fsave, * )  nsymgam
        do k=1,nsymgam
           READ( iuna2Fsave, * )  ((sgam(i,j,k),j=1,3),i=1,3)
        enddo
        READ( iuna2Fsave, * )  ((irtgam(i,j),i=1,nsymgam),j=1,nat)
        !
        CLOSE( UNIT = iuna2Fsave, STATUS = 'KEEP' )
     END IF
     !
     ! broadcast all variables read
     !
     CALL mp_bcast (etfit, ionode_id)
     CALL mp_bcast (xkfit, ionode_id)
     CALL mp_bcast (wkfit, ionode_id)
     CALL mp_bcast (nk1fit, ionode_id)
     CALL mp_bcast (nk2fit, ionode_id)
     CALL mp_bcast (nk3fit, ionode_id)
     CALL mp_bcast (nsymgam, ionode_id)
     CALL mp_bcast (sgam, ionode_id)
     CALL mp_bcast (irtgam, ionode_id)
     !
     nkfit=nk1fit*nk2fit*nk3fit
     !
     ! find S^{-1} for q=0 
     !
     call multable (nsymgam, sgam, tablegam)
     do k = 1, nsymgam
        do nn = 1, nsymgam
           if (tablegam (k, nn) ==  1 ) invsgam (k) = nn
        enddo
     enddo
     !
     ! find rtau = S tau for q=0 
     !
     symgam(1:nsymgam) = .true.
     CALL sgam_ph (at, bg, nsymgam, sgam, irtgam, tau, rtaugam, nat, symgam)
     !
     ! calculate Ef and DOS(Ef) using dense grid
     !
     do isig=1,nsig
        deg(isig) = 0.005d0*isig
        effit(isig) = efermig &
             ( etfit, nbnd, nksfit, nelec, wkfit, deg(isig), ngauss1, 0, isk)
        dosfit(isig) = dos_ef ( ngauss1, deg(isig), effit(isig), etfit, &
             wkfit, nksfit, nbnd) / 2.0d0
     enddo
     allocate (eqkfit(nkfit), eqqfit(nkfit), sfit(nkfit))
     !
     ! map k-points in the IBZ to k-points in the complete uniform grid
     !
     call lint ( nsymgam, sgam, .true., at, bg, npk, 0,0,0, &
          nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
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
     IF ( lgamma ) THEN
        call lint ( nsym, s, minus_q, at, bg, npk, 0,0,0, &
             nk1,nk2,nk3, nks, xk, 1, nkBZ, eqBZ, sBZ)
     ELSE
        call lint ( nsym, s, minus_q, at, bg, npk, 0,0,0, &
             nk1,nk2,nk3, nks, xk, 2, nkBZ, eqBZ, sBZ)
     END IF
     !
     allocate (gf(3*nat,3*nat,nsig))
     gf = (0.0d0,0.0d0)
     !
     wqa  = 2.0d0/nkfit
     !
     do ibnd = 1, nbnd
        do jbnd = 1, nbnd
           allocate (g3(nkBZ,3*nat,3*nat))
           allocate (g0(nksq,3*nat,3*nat))
           do ik = 1, nksq
              do ii = 1, 3*nat
                 do jj = 1, 3*nat
                    g0(ik,ii,jj)=conjg(el_ph_mat(jbnd,ibnd,ik,ii))* &
                         el_ph_mat(jbnd,ibnd,ik,jj)
                 enddo    ! ipert
              enddo    !jpert
           enddo   ! ik
           !
           allocate (g1(3*nat,3*nat))
           allocate (g2(3*nat,3*nat))
           do i=1,nk1
              do j=1,nk2
                 do k=1,nk3
                    nn = k-1 + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                    itemp1 = eqBZ(nn)
                    g1(:,:) = g0(itemp1,:,:)
                    itemp2 = sBZ(nn)
                    call symm( g1, g2, u, xq, s, itemp2, rtau, irt, &
                         at, bg, nat)
                    g3(nn,:,:) = g2(:,:)
                 enddo ! k
              enddo !j
           enddo !i
           deallocate (g2)
           deallocate (g1)
           deallocate (g0)
           !
           allocate ( point(nkBZ), noint(nkfit), ctemp(nkfit) )
           do jpert = 1, 3 * nat
              do ipert = 1, 3 * nat
                 !
                 point(:) = g3(:,ipert,jpert)
                 !
                 CALL clinear(nk1,nk2,nk3,nti,ntj,ntk,point,noint)
                 !
                 do isig = 1, nsig
                    degauss2 = deg(isig)
                    do ik=1,nkfit
                       etk = etfit(ibnd,eqkfit(ik))
                       etq = etfit(jbnd,eqqfit(ik))
                       w0g1 = w0gauss( (effit(isig)-etk) &
                                      / degauss2,ngauss2) / degauss2 
                       w0g2 = w0gauss( (effit(isig)-etq) &
                                      / degauss2,ngauss2) / degauss2
                       ctemp(ik) = noint(ik)* wqa * w0g1 * w0g2
                    enddo
                    gf(ipert,jpert,isig) = gf(ipert,jpert,isig) + &
                         SUM (ctemp) 
                 enddo ! isig 
              enddo    ! ipert
           enddo    !jpert
           deallocate (point, noint, ctemp)
           deallocate (g3)
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
           if (sqrt(abs(w2(nu))) > epsw) then
              ! lambda is the adimensional el-ph coupling for mode nu:
              ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
              lamb(nu,isig) = gam(nu,isig)/pi/w2(nu)/dosfit(isig)
           else
              lamb(nu,isig) = 0.0d0
           endif
        enddo  !nu
     enddo  ! isig
     !
     do isig= 1,nsig
        degauss2 = deg(isig) 
        write(stdout,'(A13,F6.4)') ' Broadening ',degauss2
        do nu=1,3*nat
           write(stdout,'(A6,I4,3F20.10)') ' mode ',nu, &
                sqrt(w2(nu))*RytoTHz, lamb(nu,isig), gam(nu,isig)*RytoGHz
           gam(nu,isig) = gam(nu,isig)*RytoGHz
        enddo
     enddo
     write(stdout,*)
     do isig= 1,nsig
        degauss2 = deg(isig)
        write(stdout,'(F6.4,3F20.10)') degauss2, (gam(nu,isig),nu=1,3*nat)
     enddo
     deallocate (gam)
     deallocate (lamb)
     write(stdout,*)
     write(stdout,*)
     write(stdout,*)
     !
     !    Prepare interface to q2r and matdyn
     !
     call star_q (xq, at, bg, ibrav, symm_type, nat, tau, ityp, nr1, &
          nr2, nr3, nsymloc, sloc, invsloc, irtloc, rtauloc, nqloc, sxqloc, &
          isqloc, imqloc, noinv, modenum)
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
        write(iuelph,*) nqloc
        call q2qstar_ph (dyn22, at, bg, nat, nsymloc, sgam, invsgam, &
             irtgam, rtaugam, nqloc, sxqloc, isqloc, imqloc, iuelph)
        if (ionode) CLOSE( UNIT = iuelph, STATUS = 'KEEP' )
     enddo
     deallocate (gf)
     !    
  ELSE
     !
     ! ======================================================================  
     !
     IF (filelph.NE.' ') THEN
        ! parallel case: only first node writes
        IF ( ionode ) THEN
           !
           iuelph = 4
           OPEN (unit = iuelph, file = filelph, status = 'unknown', err = &
                100, iostat = ios)
100        CALL errore ('elphon', 'opening file'//filelph, ABS (ios) )
           REWIND (iuelph)
           WRITE (iuelph, '(3f15.8,2i8)') xq, nsig, 3*nat
           WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
           !
        ELSE
           iuelph = 0
        END IF
        !
     ELSE
        iuelph = 0
     ENDIF
     !
     DO isig = 1, nsig
        degauss1 = 0.01 * isig
        el_ph_sum(:,:) = (0.d0, 0.d0)
        phase_space = 0.d0
        !
        ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
        ! for this gaussian broadening
        !
        ! Note that the weights of k+q points must be set to zero for the
        ! following call to yield correct results
        !
        ef1 = efermig (et, nbnd, nks, nelec, wk, degauss1, ngauss1, 0, isk)
        dosef = dos_ef (ngauss1, degauss1, ef1, et, wk, nks, nbnd)
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
              DO jbnd = 1, nbnd
                 w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                      / degauss1
                 ! note that wk(ikq)=wk(ikk)
                 weight = wk (ikk) * w0g1 * w0g2
                 DO jpert = 1, 3 * nat
                    DO ipert = 1, 3 * nat
                       el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert) + weight * &
                            CONJG (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                            el_ph_mat (jbnd, ibnd, ik, jpert)
                    ENDDO
                 ENDDO
                 phase_space = phase_space + weight
              ENDDO
           ENDDO
           
        ENDDO  ! nksq
        !
        ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
        !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
        !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
        !
        ! collect contributions from all pools (sum over k-points)
        !
        CALL poolreduce (2 * 3 * nat * 3 * nat, el_ph_sum)
        CALL poolreduce (1, phase_space)
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
        !
        DO nu = 1, nmodes
           gamma = 0.0
           DO mu = 1, 3 * nat
              DO vu = 1, 3 * nat
                 gamma = gamma + DBLE (CONJG (dyn (mu, nu) ) * el_ph_sum (mu, vu)&
                      * dyn (vu, nu) )
              ENDDO
           ENDDO
           gamma = pi * gamma / 2.d0
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
           IF (SQRT (ABS (w2 (nu) ) ) > epsw) THEN
              ! lambda is the adimensional el-ph coupling for mode nu:
              ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
              lambda = gamma / pi / w2 (nu) / dosef
           ELSE
              lambda = 0.0
           ENDIF
           WRITE (6, 9010) nu, lambda, gamma * RytoGHz
           IF (iuelph.NE.0) WRITE (iuelph, 9010) nu, lambda, gamma*RytoGHz
        ENDDO
     ENDDO
9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i2,')=',f8.4,'   gamma=',f8.2,' GHz')
     !
     !
  ENDIF
  !
  IF (iuelph.NE.0) CLOSE (unit = iuelph)
  !
  RETURN
END SUBROUTINE elphsum
!-----------------------------------------------------------------------
FUNCTION dos_ef (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
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
  dos_ef = 0.0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef = dos_ef + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  !    Collects partial sums on k-points from all pools
  !
  CALL poolreduce (1, dos_ef)
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
  integer ::  i,j,k, ns, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
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
