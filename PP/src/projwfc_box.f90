!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! real-space local DOS code courtesy of Guido Fratesi
!
MODULE projections_ldos
  USE kinds, ONLY : DP
  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
END MODULE projections_ldos
!
!-----------------------------------------------------------------------
SUBROUTINE projwave_boxes( filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE run_info,  ONLY: title
  USE atom
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev
  USE gvect
  USE gvecs,     ONLY: dual
  USE gvecw,     ONLY: ecutwfc
  USE klist,     ONLY: xk, nks, nkstot, ngk, igk_k
  USE lsda_mod,  ONLY: nspin, isk, current_spin, lsda
  USE wvfct, ONLY: npwx, nbnd, et, wg
  USE control_flags, ONLY: gamma_only
  USE uspp,      ONLY: okvan
  USE noncollin_module, ONLY: noncolin, npol
  USE wavefunctions_module, ONLY: evc,    psic
  USE wavefunctions_module, ONLY: psic_nc
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE scf,                  ONLY : rho
  USE projections_ldos,     ONLY : proj
  USE fft_base,             ONLY : dfftp
  USE scatter_mod,          ONLY : scatter_grid
  USE fft_interfaces,       ONLY : invfft
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
!
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  CHARACTER (len=256) :: filpdos
  CHARACTER (len=*) :: filproj
  INTEGER :: n_proj_boxes, irmin(3,*), irmax(3,*), nri(3)
  LOGICAL :: plotboxes
  !
  INTEGER :: ik, ibnd, i, ir, ig, ipol, npw
  INTEGER :: ibox, ir1, ir2, ir3, c_tab, is, iunproj
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  COMPLEX(DP), ALLOCATABLE :: caux(:)
  REAL(DP), ALLOCATABLE :: thetabox(:), raux(:), thetathisproc(:,:), union(:), intersection(:)
  LOGICAL, ALLOCATABLE :: isInside(:,:)
  REAL(DP), EXTERNAL :: DDOT
  REAL(DP), ALLOCATABLE :: boxvolume(:), boxcharge(:)
  !
  WRITE( stdout, '(/5x,"Calling projwave_boxes .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  IF (noncolin) THEN
     WRITE( stdout, '(/5x,"Non spin-resolved DOS will be computed")')
  ENDIF
  !
  IF (okvan) THEN
     CALL errore( 'projwave_boxes', 'Augmentation contributions are currently not included to the DOS in boxes',-1)
  ENDIF
  !
  IF ( ( n_proj_boxes > N_MAX_BOXES ) .or. ( n_proj_boxes < 1 ) ) &
     CALL errore ('projwave_boxes', 'n_proj_boxes not correct', abs (n_proj_boxes) )
  !
  ! ... Define functions with values 1.0
  ! ... on the specified boxes and 0.0 elsewhere.
  !
  ALLOCATE( thetabox (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  !
  ALLOCATE( thetathisproc(dfftp%nnr,1:n_proj_boxes) )
  !
  ALLOCATE ( isInside ( max(dfftp%nr1,dfftp%nr2,dfftp%nr3), 3 ) )
  !
  DO ibox = 1, n_proj_boxes
     !
     ! A. Do the three directions independently:
     nri(1)=dfftp%nr1
     nri(2)=dfftp%nr2
     nri(3)=dfftp%nr3
     DO i = 1, 3
        ! boxes include the points in [irmin,irmax] if irmin<=irmax
        ! and the points in [1,irmax] and [irmin,nr] if irmin > irmax
        irmin(i,ibox)=mod(irmin(i,ibox),nri(i))
        IF (irmin(i,ibox)<=0) irmin(i,ibox)=irmin(i,ibox)+nri(i)
        irmax(i,ibox)=mod(irmax(i,ibox),nri(i))
        IF (irmax(i,ibox)<=0) irmax(i,ibox)=irmax(i,ibox)+nri(i)
        DO ir = 1, nri(i)
           IF (irmin(i,ibox)<=irmax(i,ibox)) THEN
              isInside(ir,i)=(ir>=irmin(i,ibox)).and.(ir<=irmax(i,ibox))
           ELSE
              isInside(ir,i)=(ir>=irmin(i,ibox)).or. (ir<=irmax(i,ibox))
           ENDIF
        ENDDO
     ENDDO
     !
     ! B. Combine the conditions for the three directions to form a box
     ir=0
     DO ir3 = 1, dfftp%nr3
        DO ir2 = 1, dfftp%nr2
           DO ir1 = 1, dfftp%nr1
              ir=ir+1
              IF ( isInside(ir1,1) .and. &
                   isInside(ir2,2) .and. &
                   isInside(ir3,3)         ) THEN
                 thetabox(ir)=1._DP
              ELSE
                 thetabox(ir)=0._DP
              ENDIF
           ENDDO
        ENDDO
        !
     ENDDO
     !
     ! C. Output the functions thetabox in the XCrySDen format,
     ! so that the projection boxes can be visualised.
     IF ( ionode .and. plotboxes ) THEN
        filextension='.box#'
        !             123456
        c_tab = 6
        IF (ibox < 10) THEN
           WRITE (filextension( c_tab : c_tab ),'(i1)') ibox
           c_tab = c_tab + 1
        ELSEIF (ibox < 100) THEN
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') ibox
           c_tab = c_tab + 2
        ELSEIF (ibox < 1000) THEN
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') ibox
           c_tab = c_tab + 3
        ELSE
           CALL errore('projwave_boxes',&
                'file extension not supporting so many boxes', n_proj_boxes)
        ENDIF
        !
        fileout = trim(filpdos)//trim(filextension)//'.xsf'
        OPEN (4,file=fileout,form='formatted', status='unknown')
        CALL xsf_struct (alat, at, nat, tau, atm, ityp, 4)
        CALL xsf_fast_datagrid_3d(thetabox(1:dfftp%nr1x*dfftp%nr2x*dfftp%nr3x),&
                 dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, 4)
        CLOSE (4)
        !
     ENDIF
     !
#if defined(__MPI)
     CALL scatter_grid ( dfftp, thetabox(:), thetathisproc(:,ibox) )
#else
     thetathisproc(:,ibox) = thetabox(1:dfftp%nnr)      
#endif
     !
  ENDDO
  !
  DEALLOCATE ( isInside )
  DEALLOCATE ( thetabox )
  !
  !
  ! ... For each box output the volume and the electronic charge contained
  !
  ALLOCATE ( boxvolume (1:n_proj_boxes) )
  ALLOCATE ( boxcharge (1:n_proj_boxes) )
  ALLOCATE ( raux (dfftp%nnr) )
  !
  ! A. Integrate the volume
  DO ibox = 1, n_proj_boxes
     boxvolume(ibox) = sum(thetathisproc(1:dfftp%nnr,ibox))
     CALL mp_sum ( boxvolume(ibox) , intra_pool_comm )
  ENDDO
  !
  ! B1. Copy the total charge density to raux
  IF (noncolin) THEN
     CALL DCOPY (dfftp%nnr, rho%of_r, 1, raux, 1)
  ELSE
     CALL DCOPY (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
     DO is = 2, nspin
        CALL DAXPY (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
     ENDDO
  ENDIF
  !
  ! B2. Integrate the charge
  !     the correct integral has dv = omega/(nr1*nr2*nr3)
  !     not  omega/(nr1x*nr2x*nr3x) . PG 24 Oct 2010
  DO ibox = 1, n_proj_boxes
     boxcharge(ibox) = DDOT(dfftp%nnr,raux(:),1,thetathisproc(:,ibox),1) &
          &   * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
     CALL mp_sum ( boxcharge(ibox) , intra_pool_comm )
  ENDDO
  !
  ! C. Write the result
  IF (ionode) THEN
     WRITE (stdout,*)
     DO ibox = 1, n_proj_boxes
        WRITE (stdout, &
             '(5x,"Box #",i3," : vol ",f10.6," % = ",f14.6," (a.u.)^3; ",e13.6," elec")') &
             ibox, 100* boxvolume(ibox) /(dfftp%nr1*dfftp%nr2*dfftp%nr3), &
             omega* boxvolume(ibox)/(dfftp%nr1*dfftp%nr2*dfftp%nr3), boxcharge(ibox)
     ENDDO
  ENDIF
  !
  DEALLOCATE ( boxvolume , boxcharge )
  !
  ! ... Here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge in the specified box
  !
  ALLOCATE( proj(1:n_proj_boxes,nbnd,nkstot) )
  proj(:,:,:)=0._DP
  !
  ALLOCATE( caux(dfftp%nnr) )
  !
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     bnd_loop: DO ibnd = 1, nbnd
        !
        IF (noncolin) THEN
           !
           psic_nc = (0.d0,0.d0)
           DO ig = 1, npw
              psic_nc(nl(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
              psic_nc(nl(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
           ENDDO
           raux=0._DP
           DO ipol=1,npol
              CALL invfft ('Dense', psic_nc(:,ipol), dfftp)
              raux(:) = raux(:)+dble( psic_nc(:,ipol) )**2 &
                             + aimag( psic_nc(:,ipol) )**2
           ENDDO
           !
        ELSE
           !
           caux(1:dfftp%nnr) = (0._DP,0._DP)
           DO ig = 1, npw
              caux (nl (igk_k (ig,ik) ) ) = evc (ig, ibnd)
           ENDDO
           IF (gamma_only) THEN
              DO ig = 1, npw
                 caux (nlm(igk_k (ig,ik) ) ) = conjg(evc (ig, ibnd))
              ENDDO
           ENDIF
           CALL invfft ('Dense', caux, dfftp)
           !
           raux(:) = dble( caux(:) )**2 + aimag( caux(:) )**2
           !
        ENDIF
        !
        ! The contribution of this wavefunction to the LDOS
        ! integrated in the volume is the projection of the
        ! squared wfc on a function =1 in the volume itself:
        !
        DO ibox = 1, n_proj_boxes
           proj(ibox,ibnd,ik) = DDOT(dfftp%nnr,raux(:),1,thetathisproc(:,ibox),1) &
                &               / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
        ENDDO
        !
     ENDDO bnd_loop
     !
     CALL mp_sum ( proj(:,:,ik) , intra_pool_comm )
     !
  ENDDO k_loop
  !
  DEALLOCATE ( caux )
  DEALLOCATE ( raux )
  DEALLOCATE ( thetathisproc )
  !
  !   vector proj is distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (proj, n_proj_boxes*nbnd, nkstot, nks)
  !
  ! Output the projections
  IF ( ionode ) THEN
     IF (filproj/=' ') THEN
        iunproj=33
        CALL write_io_header(filproj, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
           dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
           celldm, at, gcutm, dual, ecutwfc, nkstot,nbnd,natomwfc)
        DO ibox = 1, n_proj_boxes
           WRITE (iunproj,'(3i6)') ibox, n_proj_boxes
           WRITE (iunproj,'(i6,i6,f9.4,e13.6)') &
                ((ik,ibnd,et(ibnd,ik)*rytoev,proj(ibox,ibnd,ik),ibnd=1,nbnd),ik=1,nkstot)
        ENDDO
        CLOSE (iunproj)
     ENDIF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE projwave_boxes
!
!-----------------------------------------------------------------------
SUBROUTINE partialdos_boxes(Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE klist,      ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod,   ONLY: nspin, isk, current_spin
  USE wvfct,      ONLY: et, nbnd
  USE constants,  ONLY: rytoev
  USE projections_ldos, ONLY: proj
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  LOGICAL :: kresolveddos
  INTEGER :: n_proj_boxes
  !
  CHARACTER (len=33) :: filextension
  CHARACTER (len=256):: fileout
  !
  INTEGER :: ik, ibnd, ne, ie_mid, ie_delta, ie, is, nkseff, ikeff, ibox, nspin0
  REAL(DP) :: etev, delta, Elw, Eup, wkeff
  REAL(DP), ALLOCATABLE :: dostot(:,:,:), dosbox(:,:,:,:), dosboxtot(:,:,:)
  REAL(DP), EXTERNAL :: w0gauss
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
  IF (nspin==2) THEN
     nspin0 = 2
  ELSE
     nspin0 = 1
  ENDIF
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
  ALLOCATE (dosbox(0:ne,1:n_proj_boxes,nspin0,nkseff))
  ALLOCATE (dostot(0:ne,nspin0,nkseff), dosboxtot(0:ne,nspin0,nkseff) )
  dosbox(:,:,:,:) = 0.d0
  dostot(:,:,:) = 0.d0
  dosboxtot(:,:,:)= 0.d0
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1
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
           DO ibox = 1, n_proj_boxes
              dosbox(ie,ibox,current_spin,ikeff)               = &
                   dosbox(ie,ibox,current_spin,ikeff)          + &
                   wkeff * delta * proj (ibox, ibnd, ik)
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
  ! dosboxtot(:,ns,ik) = sum of all projected DOS
  !
  DO ik=1,nkseff
     DO is=1,nspin0
        DO ie=0,ne
           dosboxtot(ie,is,ik) = sum(dosbox(ie,1:n_proj_boxes,is,ik))
        ENDDO
     ENDDO
  ENDDO
  !
  fileout = trim(filpdos)//'.ldos_boxes'
  !
  OPEN (4,file=fileout,form='formatted', &
       status='unknown')

  IF (kresolveddos) THEN
     WRITE (4,'("# ik   ")', advance="NO")
  ELSE
     WRITE (4,'("#")', advance="NO")
  ENDIF
  IF (nspin0 == 2) THEN
     WRITE (4,'(" E (eV)  tot_up(E)  tot_dw(E)  totldos_up totldos_dw ")', advance="NO")
  ELSE
     WRITE (4,'(" E (eV)  tot(E)     totldos    ")', advance="NO")
  ENDIF
  DO ibox=1, n_proj_boxes
     IF (nspin0 == 2) THEN
        WRITE(4,'("#",i3," up(E) ")', advance="NO")  ibox
        WRITE(4,'("#",i3," dw(E) ")', advance="NO")  ibox
     ELSE
        WRITE(4,'("#",i3," (E)   ")', advance="NO")  ibox
     ENDIF
  ENDDO
  WRITE (4,*)
  DO ik=1,nkseff
     DO ie= 0, ne
        IF (kresolveddos) THEN
           WRITE (4,'(i5," ")', advance="NO") ik
        ENDIF
        etev = Emin + ie * DeltaE
        WRITE (4,'(f7.3,4(2e11.3),999(2e11.3))') etev*rytoev,  &
             dostot(ie,1:nspin0,ik), dosboxtot(ie,1:nspin0,ik), &
             ( dosbox(ie,ibox,1:nspin0,ik), ibox = 1, n_proj_boxes )
     ENDDO
     IF (kresolveddos) WRITE (4,*)
  ENDDO
  CLOSE (4)
  DEALLOCATE (dostot, dosboxtot)
  DEALLOCATE (dosbox)
  !
  DEALLOCATE (proj)
  !
  RETURN
END SUBROUTINE partialdos_boxes
