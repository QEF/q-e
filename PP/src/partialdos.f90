!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE  partialdos (Emin, Emax, DeltaE, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss, ltetra
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev
  USE ktetra, ONLY: tetra_type, opt_tetra_partialdos
  !
  USE projections
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  INTEGER :: nproj, ik, ibnd,  m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff
  REAL(DP) :: etev, delta, Elw, Eup, wkeff
  REAL(DP), ALLOCATABLE :: dostot(:,:,:), pdos(:,:,:,:), pdostot(:,:,:), &
       ldos(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !
  ! this can be either natomwfc or nkb, depending upon the projection
  !
  nproj = SIZE(proj,1)
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  DO ik = 2, nkstot
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  ENDDO
  IF (degauss/=0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = max (Emin/rytoev, Elw)
  Emax = min (Emax/rytoev, Eup)
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (kresolveddos) THEN
     IF ( nspin==2 ) THEN
        nkseff=nkstot/2
     ELSE
        nkseff=nkstot
     ENDIF
  ELSE
     nkseff=1
  ENDIF
  !
  ALLOCATE (pdos(0:ne,nproj,nspin,nkseff))
  ALLOCATE (dostot(0:ne,nspin,nkseff), pdostot(0:ne,nspin,nkseff), ldos(0:ne,nspin,nkseff) )
  pdos(:,:,:,:) = 0.d0
  dostot(:,:,:) = 0.d0
  pdostot(:,:,:)= 0.d0
  !
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1

  IF ( ltetra .AND. tetra_type == 1 .OR. tetra_type == 2 ) THEN
     !
     ! Tetrahedron method
     !
     CALL opt_tetra_partialdos(nspin,kresolveddos,ne,nproj,nkseff, &
          &                    Emin,DeltaE, proj, pdos, dostot, nspin)
     !
  ELSE
     !
     DO ik = 1,nkstot
        !
        IF (kresolveddos) THEN
           ! set equal weight to all k-points
           wkeff=1.D0
           !
           IF (( nspin==2 ).AND.( isk(ik)==2 )) THEN
              ikeff=ik-nkstot/2
           ELSE
              ikeff=ik
           ENDIF
        ELSE
           ! use true weights
           wkeff=wk(ik)
           ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
           ikeff=1
        ENDIF
        !
        IF ( nspin == 2 ) current_spin = isk ( ik )
        DO ibnd = 1, nbnd
           etev = et(ibnd,ik)
           ie_mid = nint( (etev-Emin)/DeltaE )
           DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
              delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                   / degauss / rytoev
           !
           ! pdos(:,nwfc,ns,ik) = DOS (states/eV) for spin "ns"
           !                      projected over atomic wfc "nwfc"
           !                      for k-point "ik" (or summed over all kp)
           !
              DO nwfc = 1, nproj
                 pdos(ie,nwfc,current_spin,ikeff) = pdos(ie,nwfc,current_spin,ikeff) + &
                      wkeff * delta * proj (nwfc, ibnd, ik)
              ENDDO
           !
           ! dostot(:,ns,ik) = total DOS (states/eV) for spin "ns"
           !                   for k-point "ik" (or summed over all kp)
           !
              dostot(ie,current_spin,ikeff) = dostot(ie,current_spin,ikeff) + &
                   wkeff * delta
           ENDDO
        ENDDO
     ENDDO
     !
  END IF
  !
  ! pdostot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin
        DO ie=0,ne
           pdostot(ie,is,ik) = sum(pdos(ie,:,is,ik))
        ENDDO
     ENDDO
  ENDDO

  DO nwfc = 1, nproj
     IF (nlmchi(nwfc)%m == 1) THEN
        filextension='.pdos_atm#'
        !             12345678901
        c_tab = 11
        IF (nlmchi(nwfc)%na < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na
           c_tab = c_tab + 1
        ELSEIF (nlmchi(nwfc)%na < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na
           c_tab = c_tab + 2
        ELSEIF (nlmchi(nwfc)%na < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na
           c_tab = c_tab + 3
        ELSEIF (nlmchi(nwfc)%na < 10000) THEN
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na
           c_tab = c_tab + 4
        ELSE
           CALL errore('partialdos',&
                'file extension not supporting so many atoms', nwfc)
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') &
             '(',trim(atm(ityp(nlmchi(nwfc)%na)))
        c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1
        IF (nlmchi(nwfc)%l > 3) &
             CALL errore('partialdos',&
             'file extension not supporting so many l', nwfc)
        IF (nlmchi(nwfc)%n < 10) THEN
           WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  &
                nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
        ELSE IF (nlmchi(nwfc)%n < 100) THEN
           WRITE (filextension(c_tab:),'(")_wfc#",i2,"(",a1,")")')  &
                nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
        ELSE
           CALL errore('partialdos',&
             'file extension not supporting so many atomic wfc', nwfc)
        END IF
        fileout = trim(filpdos)//trim(filextension)
        OPEN (4,file=fileout,form='formatted', status='unknown')

        IF (kresolveddos) THEN
           WRITE (4,'("# ik   ")', advance="NO")
        ELSE
           WRITE (4,'("#")', advance="NO")
        ENDIF
        IF (nspin == 1) THEN
           WRITE (4,'(" E (eV)   ldos(E)  ")', advance="NO")
        ELSE
           WRITE (4,'(" E (eV)  ldosup(E)  ldosdw(E)")', advance="NO")
        ENDIF
        DO m=1,2 * nlmchi(nwfc)%l + 1
           IF (nspin == 1) THEN
              WRITE(4,'(" pdos(E)   ")', advance="NO")
           ELSE
              WRITE(4,'(" pdosup(E) ")', advance="NO")
              WRITE(4,'(" pdosdw(E) ")', advance="NO")
           ENDIF
        ENDDO
        WRITE(4,*)
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:,:) = 0.d0
        DO ik=1,nkseff
           DO ie= 0, ne
              DO is=1, nspin
                 DO m=1,2 * nlmchi(nwfc)%l + 1
                    ldos  (ie, is, ik) = ldos  (ie, is, ik) + pdos(ie,nwfc+m-1,is,ik)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        DO ik=1,nkseff
           DO ie= 0, ne
              IF (kresolveddos) THEN
                 WRITE (4,'(i5," ")', advance="NO") ik
              ENDIF
              etev = Emin + ie * DeltaE
              WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  &
                   (ldos(ie,is,ik), is=1,nspin), &
                   ((pdos(ie,nwfc+m-1,is,ik), is=1,nspin), &
                   m=1,2*nlmchi(nwfc)%l+1)
           ENDDO
           IF (kresolveddos) WRITE (4,*)
        ENDDO
        CLOSE (4)
     ENDIF
  ENDDO
  fileout = trim(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown')
  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ")', advance="NO")
  ELSE
     WRITE (4,'("#")', advance="NO")
  ENDIF
  IF (nspin == 1) THEN
     WRITE (4,'(" E (eV)  dos(E)    pdos(E)")')
  ELSE
     WRITE (4,'(" E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)")')
  ENDIF
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ")', advance="NO") ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4e11.3)') etev*rytoev, (dostot(ie,is,ik), is=1,nspin), &
             (pdostot(ie,is,ik), is=1,nspin)
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi)
  DEALLOCATE (proj)
  !
  RETURN
END SUBROUTINE partialdos
!
!-----------------------------------------------------------------------
SUBROUTINE  partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss, ltetra
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev
  USE ktetra, ONLY: opt_tetra_partialdos
  !
  USE spin_orb,   ONLY: lspinorb
  USE projections
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  INTEGER :: ik, ibnd, ind, m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff, nspin0
  REAL(DP) :: etev, delta, Elw, Eup, wkeff, fact(2), spinor
  REAL(DP), ALLOCATABLE :: dostot(:,:), pdos(:,:,:,:), pdostot(:,:,:), &
       ldos(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
  !
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  DO ik = 2, nkstot
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  ENDDO
  IF (degauss/=0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = max (Emin/rytoev, Elw)
  Emax = min (Emax/rytoev, Eup)
  DeltaE = DeltaE/rytoev
  ne = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (lspinorb) THEN
     nspin0 = 1
  ELSE
     nspin0 = 2
  ENDIF
  !
  IF (kresolveddos) THEN
     nkseff=nkstot
  ELSE
     nkseff=1
  ENDIF
  !
  ALLOCATE (pdos(0:ne,natomwfc,nspin0,nkseff))
  ALLOCATE (dostot(0:ne,nkseff), pdostot(0:ne,nspin0,nkseff), ldos(0:ne,nspin0,nkseff) )
  pdos(:,:,:,:) = 0.d0
  dostot(:,:) = 0.d0
  pdostot(:,:,:)= 0.d0
  ie_delta = 5 * degauss / DeltaE + 1

  IF(ltetra) THEN
     !
     CALL opt_tetra_partialdos(nspin0,kresolveddos,ne,natomwfc,nkseff, &
     &                         Emin,DeltaE, proj, pdos, dostot, 1)
     !
     IF(.NOT. lspinorb) THEN
        DO nwfc = 1, natomwfc
           IF ( nlmchi(nwfc)%ind > (2 * nlmchi(nwfc)%l+1)) THEN
              pdos(0:ne,nwfc,2,1:nkseff) = pdos(0:ne,nwfc,1,1:nkseff)
              pdos(0:ne,nwfc,1,1:nkseff) = 0.0_dp
           ENDIF
        ENDDO
     END IF
     !
  ELSE
     DO ik = 1,nkstot
        !
        IF (kresolveddos) THEN
           ! set equal weight to all k-points
           wkeff=1.D0
           ikeff=ik
        ELSE
           wkeff=wk(ik)
           ! contributions from all k-points are summed in pdos(:,:,:,ikeff)
           ikeff=1
        ENDIF
        !
        DO ibnd = 1, nbnd
           etev = et(ibnd,ik)
           ie_mid = nint( (etev-Emin)/DeltaE )
           DO ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
              delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                   / degauss / rytoev
           !
           ! pdos(:,nwfc,ns,ik) = DOS (states/eV) for spin "ns"
           !                      projected over atomic wfc "nwfc"
           !                      for k-point "ik" (or summed over all kp)
           !
           !
           ! dostot(:,ik) = total DOS (states/eV)
           !                for k-point "ik" (or summed over all kp)
           !
              IF (lspinorb) THEN
                 DO nwfc = 1, natomwfc
                    pdos(ie,nwfc,1,ikeff) = pdos(ie,nwfc,1,ikeff) + &
                         wkeff * delta * proj (nwfc, ibnd, ik)
                 ENDDO
                 dostot(ie,ikeff) = dostot(ie,ikeff) + wkeff * delta
              ELSE
                 DO nwfc = 1, natomwfc
                    IF ( nlmchi(nwfc)%ind<=(2* nlmchi(nwfc)%l+1)) THEN
                       pdos(ie,nwfc,1,ikeff) = pdos(ie,nwfc,1,ikeff) + &
                            wkeff * delta * proj (nwfc, ibnd, ik)
                       pdos(ie,nwfc,2,ikeff) = 0.d0
                    ELSE
                       pdos(ie,nwfc,1,ikeff) = 0.d0
                       pdos(ie,nwfc,2,ikeff) = pdos(ie,nwfc,2,ikeff) + &
                            wkeff * delta * proj (nwfc, ibnd, ik)
                    ENDIF
                 ENDDO
                 dostot(ie,ikeff) = dostot(ie,ikeff) + wkeff * delta
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     !
  END IF
  !
  ! pdostot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin0
        DO ie=0,ne
           pdostot(ie,is,ik) = sum(pdos(ie,:,is,ik))
        ENDDO
     ENDDO
  ENDDO

  DO nwfc = 1, natomwfc
     IF (nlmchi(nwfc)%ind == 1) THEN
        filextension='.pdos_atm#'
        !             12345678901
        c_tab = 11
        IF (nlmchi(nwfc)%na < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na
           c_tab = c_tab + 1
        ELSEIF (nlmchi(nwfc)%na < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na
           c_tab = c_tab + 2
        ELSEIF (nlmchi(nwfc)%na < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na
           c_tab = c_tab + 3
        ELSEIF (nlmchi(nwfc)%na < 10000) THEN
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na
           c_tab = c_tab + 4
        ELSE
           CALL errore('partialdos_nc',&
                'file extension not supporting so many atoms', nwfc)
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') &
             '(',trim(atm(ityp(nlmchi(nwfc)%na)))
        c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1
        IF (nlmchi(nwfc)%l > 3) &
             CALL errore('partialdos_nc',&
             'file extension not supporting so many l', nwfc)
        IF (nlmchi(nwfc)%n < 10) THEN
           IF (lspinorb) THEN
             WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,"_j",f3.1,")")') &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l),nlmchi(nwfc)%jj
           ELSE
             WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
           ENDIF
        ELSE IF (nlmchi(nwfc)%n < 100) THEN
           IF (lspinorb) THEN
             WRITE (filextension(c_tab:),'(")_wfc#",i2,"(",a1,"_j",f3.1,")")') &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l),nlmchi(nwfc)%jj
           ELSE
             WRITE (filextension(c_tab:),'(")_wfc#",i2,"(",a1,")")')  &
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
           ENDIF
        ELSE
           CALL errore('partialdos_nc',&
             'file extension not supporting so many atomic wfc', nwfc)
        ENDIF
        fileout = trim(filpdos)//trim(filextension)
        OPEN (4,file=fileout,form='formatted', status='unknown')

        IF (kresolveddos) THEN
           WRITE (4,'("# ik   ")', advance="NO")
        ELSE
           WRITE (4,'("#")', advance="NO")
        ENDIF
        IF (nspin0 == 1) THEN
           WRITE (4,'(" E(eV)   ldos(E)   ")', advance="NO")
        ELSE
           WRITE (4,'(" E(eV)  ldosup(E)  ldosdw(E)")', advance="NO")
        ENDIF
        IF (lspinorb) THEN
           ind = 0
           DO m = -nlmchi(nwfc)%l-1, nlmchi(nwfc)%l
              fact(1) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,1)
              fact(2) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,2)
              IF (abs(fact(1))>1.d-8.or.abs(fact(2))>1.d-8) THEN
                 ind = ind + 1
                 WRITE(4,'("pdos(E)_",i1,"   ")', advance="NO") ind
              ENDIF
           ENDDO
        ELSE
           DO ind=1,2 * nlmchi(nwfc)%l + 1
              WRITE(4,'(" pdosup(E) ")', advance="NO")
              WRITE(4,'(" pdosdw(E) ")', advance="NO")
           ENDDO
        ENDIF
        WRITE(4,*)
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:,:) = 0.d0
        IF (lspinorb) THEN
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0)<1.d-8) THEN
                    DO ind = 1, 2 * nlmchi(nwfc)%l + 2
                       ldos  (ie, 1, ik) = ldos  (ie, 1, ik) + pdos(ie,nwfc+ind-1,1,ik)
                    ENDDO
                 ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0)<1.d-8) THEN
                    DO ind = 1, 2 * nlmchi(nwfc)%l
                       ldos  (ie, 1, ik) = ldos  (ie, 1, ik) + pdos(ie,nwfc+ind-1,1,ik)
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (kresolveddos) THEN
                    WRITE (4,'(i5," ")', advance="NO") ik
                 ENDIF
                 etev = Emin + ie * DeltaE
                 IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0)<1.d-8) THEN
                    WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1,ik), &
                         (pdos(ie,nwfc+ind-1,1,ik), ind=1,2*nlmchi(nwfc)%l+2)
                 ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0)<1.d-8) THEN
                    WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1,ik), &
                         (pdos(ie,nwfc+ind-1,1,ik), ind=1,2*nlmchi(nwfc)%l)
                 ENDIF
              ENDDO
              IF (kresolveddos) WRITE (4,*)
           ENDDO
        ELSE
           DO ik=1,nkseff
              DO ie= 0, ne
                 DO is=1, nspin0
                    DO ind=1,4 * nlmchi(nwfc)%l + 2
                       ldos  (ie, is, ik) = ldos  (ie, is, ik) + pdos(ie,nwfc+ind-1,is, ik)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           DO ik=1,nkseff
              DO ie= 0, ne
                 IF (kresolveddos) THEN
                    WRITE (4,'(i5," ")', advance="NO") ik
                 ENDIF
                 etev = Emin + ie * DeltaE
                 WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  &
                      (ldos(ie,is,ik), is=1,nspin0), &
                      ((pdos(ie,nwfc+ind-1+(is-1)*(2*nlmchi(nwfc)%l+1),is,ik), is=1,nspin0), &
                      ind=1,2*nlmchi(nwfc)%l+1)
              ENDDO
              IF (kresolveddos) WRITE (4,*)
           ENDDO
        ENDIF
        CLOSE (4)
     ENDIF
  ENDDO
  fileout = trim(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown')
  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ")', advance="NO")
  ELSE
     WRITE (4,'("#")', advance="NO")
  ENDIF
  IF (nspin0 == 1) THEN
     WRITE (4,'(" E (eV)  dos(E)    pdos(E)")')
  ELSE
     WRITE (4,'(" E (eV)  dos(E)   pdosup(E)  pdosdw(E)")')
  ENDIF
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ")', advance="NO") ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4e11.3)') etev*rytoev, dostot(ie,ik), &
             (pdostot(ie,is,ik), is=1,nspin0)
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi)
  DEALLOCATE (proj)
  !
  RETURN
END SUBROUTINE partialdos_nc
