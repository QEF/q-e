!
! Copyright (C) 2009 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine scat_states_plot(ik,ien,norb,nocros,nchan,vec,veceig,left_to_right)
!
! Writes the the XY integrated and the 3D charge and spin densities of
! right-moving scattering states (or Bloch states if ikind = 0).
!
  use kinds,            ONLY : DP
  USE constants,        ONLY : tpi, rytoev
  use io_global,        ONLY : stdout, ionode
  USE ions_base,        ONLY : ityp, tau, nat, atm
  use noncollin_module, ONLY : noncolin, npol
  USE spin_orb,         ONLY : domag
  use lsda_mod,         ONLY : nspin
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : alat, at
  USE cond,             ONLY : ikind, n2d, nrzpl, nrzps, taunewl, taunews,  &
                               lorb3d, funz0, lcharge, denergy, rho_scatt,  &
                               sarea, nenergy, nkpts, wkpt

  implicit none

  INTEGER  :: ik, ien, nspin0, ichan, nchan, ipol, ix, iy, iz, ij, ounit, ios, norb, nocros, &
              c_tab
  real(DP) :: raux1, raux2
  real(DP), allocatable    :: spin_mag(:,:,:), zdata(:), mag(:), aux_plot(:)
  COMPLEX(DP), PARAMETER   :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  COMPLEX(DP) :: vec(4*n2d+npol*(norb+2*nocros),nchan), veceig(nchan,nchan)
  CHARACTER(LEN=50) :: filename
  LOGICAL :: norm_flag, left_to_right

  ALLOCATE( spin_mag(nspin,nchan,dfftp%nr1*dfftp%nr2*dfftp%nr3) )

!---
! Construct the states

  if (ikind.eq.0) then
    if (lcharge) then
     norm_flag = .FALSE.
    else
     norm_flag = .TRUE.
    endif
    CALL scat_states_comp(nchan, nrzpl, norb, nocros, &
               taunewl, vec, veceig, spin_mag, norm_flag, left_to_right)
  else
    CALL scat_states_comp(nchan, nrzps, norb, nocros, &
               taunews, vec, veceig, spin_mag, .FALSE., left_to_right)
  endif
!---

  allocate( zdata(dfftp%nr3) )
  allocate( aux_plot(dfftp%nr1*dfftp%nr2*dfftp%nr3) )

!-- z-mezh (in \AA)
  raux1 = at(3,3)*alat*0.5291772108d0 / dfftp%nr3
  zdata(1) = 0.d0
  do iz = 2, dfftp%nr3
   zdata(iz) = zdata(iz-1) + raux1
  enddo
!---

!-- the densities integrated in the XY plane
  write(stdout,*)
  if (ikind.eq.0) then
    if (left_to_right) then
     write(stdout,*) 'RIGHT MOVING Bloch states (integrated in XY) as a function of z:'
    else
     write(stdout,*) 'LEFT MOVING Bloch states (integrated in XY) as a function of z:'
    endif
  else
    if (left_to_right) then
     write(stdout,*) 'RIGHT MOVING scatt. states (integrated in XY) as a function of z:'
    else
     write(stdout,*) 'LEFT MOVING scatt. states (integrated in XY) as a function of z:'
    endif
  endif
  if (noncolin) then
   nspin0 = 4
   write(stdout,'(2a11,a13,a12,a16)') '--- z, Ang','n(z)','m_x(z)','m_y(z)','m_z(z) ---'
  else
   nspin0 = 1
   write(stdout,'(2a11,a13)') '--- z, Ang','n(z)'
  endif
  allocate( mag(nspin0) )
!-- Put correct multiplier in case of integrating over energies and k-points
  raux2 = wkpt(ik)
  if (nenergy.gt.1) raux2 = raux2 * ABS(denergy) / (tpi*rytoev)
!--
  do ichan = 1, nchan
    write(stdout,*) 'Channel ', ichan
    do iz = 1, dfftp%nr3
      do ipol = 1, nspin0
        mag(ipol) = 0.d0
        do ix = 1, dfftp%nr1
           do iy = 1, dfftp%nr2
             ij = ix + (iy - 1) * dfftp%nr1 + (iz - 1) * dfftp%nr1 * dfftp%nr2
             mag(ipol) = mag(ipol) + spin_mag(ipol,ichan,ij)
           enddo
        enddo
        mag(ipol) = mag(ipol) * sarea / (dfftp%nr1*dfftp%nr2)
        if (lcharge) rho_scatt(iz,ipol) = rho_scatt(iz,ipol) + mag(ipol)*raux2
      enddo
      write(stdout,'(5f12.6)') zdata(iz), (mag(ipol), ipol = 1, nspin0)
    enddo
  enddo
!--

!-- Total charge and magnetization densities
  if (lcharge) then
    write(stdout,*)
    mag(:) = 0.d0
    if (noncolin) then
     write(stdout,*) '-- Total charge and magnetization --'
     write(stdout,'(2a11,a13,a12,a16)') '--- z, Ang','n(z)','m_x(z)','m_y(z)','m_z(z) ---'
    else
     write(stdout,*) '-- Total charge --'
     write(stdout,'(2a11,a13)') '--- z, Ang','n(z)'
    endif
    do iz = 1, dfftp%nr3
      write(stdout,'(5f12.6)') zdata(iz), (rho_scatt(iz,ipol), ipol = 1, nspin0)
      do ipol = 1, nspin0
        mag(ipol) = mag(ipol) + rho_scatt(iz,ipol)
      enddo
    enddo
    mag(:) = mag(:) / dfftp%nr3 * at(3,3) * alat
    if (noncolin) then
     write(stdout,'(''Nelec and Magn. Moment '',4f12.6)') (mag(ipol), ipol = 1, nspin0)
    else
     write(stdout,'(''Nelec '',f12.6)') (mag(ipol), ipol = 1, nspin0)
    endif
  endif
!--

!--- 3D output for XCRYSDENS plot
  IF (lorb3d.and.ionode) THEN
    do ichan = 1, nchan
      ounit = 34
!-- Filename
!
      if (left_to_right) then
        filename='wfc_lr_k'
      else
        filename='wfc_rl_k'
      endif
      c_tab = 9
      IF (ik>99) THEN
        write(filename( c_tab : c_tab+2 ),'(i3)') ik
        c_tab = c_tab + 3
      ELSEIF (ik>9) THEN
        write(filename( c_tab : c_tab+1 ),'(i2)') ik
        c_tab = c_tab + 2
      ELSE
        write(filename( c_tab : c_tab ),'(i1)') ik
        c_tab = c_tab + 1
      ENDIF
      WRITE (filename(c_tab:c_tab),'(a1)') 'e'
      c_tab = c_tab + 1

      IF (ien>99) THEN
        write(filename( c_tab : c_tab+2 ),'(i3)') ien
        c_tab = c_tab + 3
      ELSEIF (ien>9) THEN
        write(filename( c_tab : c_tab+1 ),'(i2)') ien
        c_tab = c_tab + 2
      ELSE
        write(filename( c_tab : c_tab ),'(i1)') ien
        c_tab = c_tab + 1
      ENDIF
      WRITE (filename(c_tab:c_tab),'(a1)') 'n'
      c_tab = c_tab + 1

      IF (ichan>99) THEN
        write(filename( c_tab : c_tab+2 ),'(i3)') ichan
        c_tab = c_tab + 3
      ELSEIF (ichan>9) THEN
        write(filename( c_tab : c_tab+1 ),'(i2)') ichan
        c_tab = c_tab + 2
      ELSE
        write(filename( c_tab : c_tab ),'(i1)') ichan
        c_tab = c_tab + 1
      ENDIF
!--
      filename=TRIM(filename)
      OPEN (UNIT=ounit, FILE=filename, FORM='formatted', &
               STATUS='unknown', ERR=100, IOSTAT=ios)
100   CALL errore('write_states','opening file'//filename,ABS(ios))
      call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
      ix = 1
      IF (noncolin.AND.domag) ix = nspin
      do ipol = 1, ix
        do ij = 1, dfftp%nr1*dfftp%nr2*dfftp%nr3
           aux_plot(ij) = spin_mag(ipol,ichan,ij)
        enddo
        call xsf_fast_datagrid_3d &
          (aux_plot, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1, dfftp%nr2, dfftp%nr3, at, alat, ounit)
      enddo
      CLOSE(ounit)
    enddo
  ENDIF
!---

  DEALLOCATE( spin_mag )
  DEALLOCATE( zdata )
  DEALLOCATE( aux_plot )

  return
end subroutine scat_states_plot

SUBROUTINE scat_states_comp(nchan, nrzp, norb, nocros, taunew, vec, &
                            veceig, spin_mag_tot, norm_flag, left_to_right)
!
! Calculates the charge and spin densities of scattering states.
!

 USE kinds,     ONLY : DP
 USE constants, ONLY : tpi
 USE noncollin_module, ONLY : noncolin, npol
 use lsda_mod,  only : nspin
 USE mp_global, ONLY : nproc_pool, me_pool, intra_pool_comm
 USE mp,        ONLY : mp_sum
 USE fft_base,  ONLY : dffts, dfftp
 USE scatter_mod, ONLY : gather_grid
 USE cond,      ONLY : ngper, newbg, intw1, intw2, &
                       nl_2ds, nl_2d, korbl, korbr, funz0, kfunl, xyk, ikind, &
                       n2d, kvall
 USE realus_scatt
 USE cell_base, ONLY : omega
 USE scf,       ONLY : rho
 USE uspp_param,ONLY : upf, nhm, nh
 USE uspp,      ONLY : nkb, vkb, becsum
 USE ions_base, ONLY : ityp, zv, nat, ntyp => nsp, tau, atm
 USE fft_scalar,ONLY : cft_2xy
 USE splinelib, ONLY : spline, splint
! USE becmod,    ONLY : bec_type, becp

 IMPLICIT NONE
 INTEGER :: nocros, nchan, nrzp, norb, irun, nrun
 COMPLEX(DP) :: x1, vec(4*n2d+npol*(norb+2*nocros),nchan), veceig(nchan,nchan)
 LOGICAL :: norm_flag, left_to_right
 INTEGER :: ik, ig, ir, mu, ig1, ichan, ichan1, ipol, iat, ih, jh, ijh, np
 INTEGER :: iorb, iorb1
 REAL(DP) :: r_aux1, r_aux2, r_aux3, taunew(4,norb)
 COMPLEX(DP), PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
 REAL(DP) :: spin_mag_tot(nspin,nchan,dfftp%nr1*dfftp%nr2*dfftp%nr3)

 INTEGER :: is, js, ix, jx
 COMPLEX(DP), ALLOCATABLE :: fung(:), amat(:,:), aux_proc(:,:), &
                             becsum_nc(:,:,:,:), vec1(:), kfunz(:,:,:)
 COMPLEX(DP), ALLOCATABLE :: funr(:)

 REAL(DP), ALLOCATABLE :: spin_mag(:,:), becsum_orig(:,:,:)
 REAL(DP), ALLOCATABLE :: xdata(:), xdatax(:), ydata(:), ydatax(:), y2d(:)
 INTEGER, ALLOCATABLE :: ipiv(:), ind(:)


 ALLOCATE( ind(nat) )
 ALLOCATE( fung(ngper) )
 ALLOCATE( funr(dfftp%nr1*dfftp%nr2) )
 ALLOCATE( aux_proc(dfftp%nr1*dfftp%nr2*nrzp,npol) )
 ALLOCATE( amat(norb*npol,norb*npol) )
 ALLOCATE( vec1(norb*npol) )
 ALLOCATE( ipiv(norb*npol) )
 ALLOCATE( spin_mag(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nspin) )
 ALLOCATE( xdata(dffts%nr3+1) )
 ALLOCATE( ydata(dffts%nr3+1) )
 ALLOCATE( xdatax(dfftp%nr3) )
 ALLOCATE( ydatax(dfftp%nr3) )
 ALLOCATE( y2d(dffts%nr3+1) )
 IF (noncolin) ALLOCATE(becsum_nc(nhm*(nhm+1)/2,nat,npol,npol))
 ALLOCATE( becsum_orig(nhm*(nhm+1)/2,nat,nspin) )

 ALLOCATE( kfunz(n2d,nchan,nrzp) )

!---  PS functions for all channels kfunz(n2d, nchan, nrzp)
 do ik = 1, nrzp
    CALL zgemm('n', 'n', n2d, nchan, 2*n2d+npol*norb, one, funz0(1,1,ik), &
            n2d, vec, 4*n2d+npol*(norb+2*nocros), zero, kfunz(1,1,ik), n2d)
 enddo
!---

!-- Relation between PW atom and PWCOND (first) orbital, iat --> ind(iat)
 ind = 0
 do iat = 1, nat
  np = ityp(iat)
  if (upf(np)%tvanp) then
    r_aux1 = 1.d0
    iorb = 0
    do while (r_aux1.gt.1.d-6)
      iorb = iorb + 1
      r_aux1 = (taunew(1,iorb)-tau(1,iat))**2+(taunew(2,iorb)-tau(2,iat))**2+&
               (taunew(3,iorb)-tau(3,iat))**2
    enddo
    ind(iat) = iorb
!   write(6,*) '...', iat, ind(iat)
  endif
 enddo
!--

 call realus_scatt_0()

!--------
! Constructs n(r) and m(r) for all propagating channels

 DO ichan = 1, nchan

!---------
!   constructs nonlocal coefficients
!
    if(ikind.eq.0) then
      if (left_to_right) then
        ichan1 = ichan
      else
        ichan1 = n2d + npol*nocros + ichan
      endif
    endif

    vec1 = 0.d0
    !--- inside orbitals
    do iorb = npol*nocros+1, npol*(norb-nocros)
      do ir=1, 2*n2d
        vec1(iorb)=vec1(iorb)+       &
                      intw1(iorb,ir)*vec(ir,ichan)
      enddo
      do ir=1, norb*npol
        vec1(iorb)=vec1(iorb)+       &
           intw2(iorb,ir)*vec(2*n2d+ir,ichan)
      enddo
    enddo
    !--- right-crossing orbitals
    if(ikind.eq.0) then
      do iorb = 1, npol*nocros
        vec1(npol*(norb-nocros)+iorb) = korbr(iorb,ichan1)
      enddo
    else
     do iorb=1, npol*nocros
        iorb1 = npol*(norb-nocros)+iorb
        do ig=1, n2d+npol*nocros
          vec1(iorb1)=vec1(iorb1) + korbr(iorb,ig)*&
             vec(3*n2d+npol*(norb+nocros)+ig,ichan)
        enddo
     enddo
    endif
    !--- left-crossing orbitals
    if(ikind.eq.0) then
      do iorb=1, npol*nocros
        vec1(iorb) = korbl(iorb,ichan1)
      enddo
    else
     !-- reflected part
     do ig=1, n2d+npol*nocros
       do iorb=1, npol*nocros
         vec1(iorb)=vec1(iorb) + korbl(iorb,n2d+npol*nocros+ig)*&
            vec(2*n2d+npol*norb+ig,ichan)
       enddo
     enddo
     !-- incident part
     if(left_to_right) then
       do iorb=1, npol*nocros
         do ig=1, nchan
           vec1(iorb)=vec1(iorb) + korbl(iorb,ig)*&
              veceig(ig,ichan)
         enddo
       enddo
     else
       do iorb=1, npol*nocros
         iorb1 = npol*(norb-nocros)+iorb
         do ig=1, nchan
           vec1(iorb1)=vec1(iorb1) + korbr(iorb,n2d+npol*nocros+ig)*&
              veceig(ig,ichan)
         enddo
       enddo
     endif
     !---
    endif
!----------

!----------
!   Construct becsum_orig and becsum for original atom and for its copy
!   in the z direction if the atom crosses the cell boundary

 becsum = 0.d0
 becsum_orig = 0.d0

 do iat = 1, nat
   np = ityp(iat)
   if (upf(np)%tvanp) then
     iorb = (ind(iat)-1)*npol + 1
     nrun = 1
!--  crossing orbs need 2 runs (for the original and copy atomic positions)
     if (iorb.le.nocros*npol.or.iorb.gt.(norb-nocros)*npol) nrun = 2
!--
     do irun = 1, nrun
        ijh = 1
        do ih = 1, nh(np)
         if(noncolin) then
           DO is=1,npol
              DO js=1,npol
                 becsum_nc(ijh,iat,is,js) =         &
                    CONJG(vec1(iorb-1+npol*(ih-1)+is)) * &
                          vec1(iorb-1+npol*(ih-1)+js)
              END DO
           END DO
         else
           becsum(ijh,iat,1) = DBLE(CONJG(vec1(iorb-1+ih))*vec1(iorb-1+ih))
         endif
         ijh = ijh + 1
         do jh = ih+1, nh(np)
           if (noncolin) THEN
              DO is=1,npol
                 DO js=1,npol
                   becsum_nc(ijh,iat,is,js) =         &
                       CONJG(vec1(iorb-1+npol*(ih-1)+is)) *  &
                             vec1(iorb-1+npol*(jh-1)+js)
                 END DO
              END DO
           else
             becsum(ijh,iat,1) = 2.d0*DBLE(CONJG(vec1(iorb-1+ih))*vec1(iorb-1+jh))
           endif
           ijh = ijh + 1
         enddo
        enddo
        IF (noncolin) THEN
          IF (upf(np)%has_so) THEN
            CALL transform_becsum_so(becsum_nc,becsum,iat)
          ELSE
            CALL transform_becsum_nc(becsum_nc,becsum,iat)
          ENDIF
        ENDIF
!--  fill in the original atomic position (becsum_orig)
        if (irun.eq.1) then
         do is = 1, nhm*(nhm+1)/2
          do js = 1, nspin
            becsum_orig(is,iat,js) = becsum(is,iat,js)
            becsum(is,iat,js) = 0.d0
          enddo
         enddo
        endif
!--  another run for the copy (becsum)
        if (nrun.eq.2.and.irun.eq.1) then
         IF (iorb.le.nocros*npol) THEN
          iorb = iorb + (norb-nocros)*npol
         ELSE
          iorb = iorb - (norb-nocros)*npol
         ENDIF
        endif
     enddo
   endif
 enddo
!-------------

!----------
!   compute and collect the densities from CPUs
!
    spin_mag = 0.d0
    rho%of_r(:,:) = 0.d0
    call realus_scatt_1(becsum_orig)
    do ipol = 1, nspin
#if defined(__MPI)
     CALL gather_grid (dfftp, rho%of_r(:,ipol),spin_mag(:,ipol))
#else
     do ig = 1, dfftp%nnr
       spin_mag(ig,ipol) = rho%of_r(ig,ipol)
     enddo
#endif
    enddo
!   noaug
!   spin_mag = 0.d0
!---

!-- Pseudo wave function in the real space
    DO ik=1,nrzp
       DO ipol=1,npol
          fung=(0.d0,0.d0)
          DO ig=1,ngper
             ig1=ig+(ipol-1)*ngper
             DO mu=1,n2d
                fung(ig)=fung(ig) + kfunz(mu,ichan,ik) * newbg(ig1,mu)
             END DO
          END DO
          funr=(0.d0,0.d0)
          DO ig=1,ngper
             funr(nl_2d(ig))=fung(ig)
          END DO
          call cft_2xy(funr, 1, dfftp%nr1, dfftp%nr2, dfftp%nr1, dfftp%nr2, 1)
          DO ix=1,dfftp%nr1
             DO jx=1,dfftp%nr2
                ir=ix + (jx - 1) * dfftp%nr1 + (ik - 1) * dfftp%nr1 * dfftp%nr2
                ig=ix+(jx-1)*dfftp%nr1
                aux_proc(ir,ipol) = funr(ig)
             END DO
          END DO
       END DO
    END DO
!--


!----------
! calculates PS density and spin magnetization
  rho%of_r(:,:) = 0.d0
  do ik = 1, nrzp
    DO ix=1,dfftp%nr1
      DO jx=1,dfftp%nr2
       ig=ix + (jx - 1) * dfftp%nr1 + (ik - 1) * dfftp%nr1 * dfftp%nr2
       ig1=ix + (jx - 1) * dfftp%nr1x + (ik - 1) * dfftp%nr1x * dfftp%nr2x
       rho%of_r(ig1,1) = DBLE(aux_proc(ig,1))**2+AIMAG(aux_proc(ig,1))**2
       IF (noncolin) THEN
        rho%of_r(ig1,2) = 2.D0*(DBLE(aux_proc(ig,1))*DBLE(aux_proc(ig,2)) + &
                           AIMAG(aux_proc(ig,1))*AIMAG(aux_proc(ig,2)))
        rho%of_r(ig1,3) = 2.D0*(DBLE(aux_proc(ig,1))*AIMAG(aux_proc(ig,2)) - &
                           DBLE(aux_proc(ig,2))*AIMAG(aux_proc(ig,1)))
        rho%of_r(ig1,4) = DBLE(aux_proc(ig,1))**2+AIMAG(aux_proc(ig,1))**2 - &
                     DBLE(aux_proc(ig,2))**2-AIMAG(aux_proc(ig,2))**2
       ENDIF
      enddo
    enddo
  enddo
!----------

!-------
! collecting the density and magnetizations on fine mesh

  r_aux1 = 1.d0/dffts%nr3
  do ig = 1, dffts%nr3+1
   xdata(ig) = r_aux1*(ig-1)
  enddo
  r_aux2 = 1.d0/dfftp%nr3
  do ig = 1, dfftp%nr3
   xdatax(ig) = r_aux2*(ig-1)
  enddo

  if(me_pool.eq.0) then
   ik = 0
  else
   ik = SUM(dffts%npp(1:me_pool))
  endif

  do ipol = 1, nspin
   DO ix = 1, dfftp%nr1
    DO jx = 1, dfftp%nr2

     ydata = 0.d0
     do ig = 1, nrzp
      ig1=ix + (jx - 1) * dfftp%nr1x + (ig - 1) * dfftp%nr1x * dfftp%nr2x
      ydata(ik+ig) = rho%of_r(ig1,ipol)
     enddo
     CALL mp_sum(ydata, intra_pool_comm )
     if(ikind.eq.0) then
       ydata(dffts%nr3+1) = ydata(1)
     else
       ydata(dffts%nr3+1) = 2.d0*ydata(dffts%nr3)-ydata(dffts%nr3-1)
     endif
     r_aux3 = 0.d0
     CALL spline( xdata, ydata, 0.d0, r_aux3, y2d )
     do ig = 1, dfftp%nr3
       ydatax(ig) = splint( xdata, ydata, y2d, xdatax(ig) )
     enddo
     do ig = 1, dfftp%nr3
       ig1=ix + (jx - 1) * dfftp%nr1x + (ig - 1) * dfftp%nr1x * dfftp%nr2x
       spin_mag(ig1,ipol) = spin_mag(ig1,ipol) + ydatax(ig)
     enddo

    END DO
   END DO
  enddo
!-------

!-------
! Normalization
 IF (norm_flag) THEN
    r_aux1 = SUM(spin_mag(:,1))
    r_aux1 = r_aux1*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
    r_aux1 = 1.d0/r_aux1
    CALL dscal(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x*nspin,r_aux1,spin_mag,1)
 END IF
!-------

!----
! save in the spin_mag_tot array
  DO ix = 1, dfftp%nr1
    DO jx = 1, dfftp%nr2
      do ig = 1, dfftp%nr3
        ir = ix + (jx - 1) * dfftp%nr1 + (ig - 1) * dfftp%nr1 * dfftp%nr2
        ig1= ix + (jx - 1) * dfftp%nr1x + (ig - 1) * dfftp%nr1x * dfftp%nr2x
        do ipol = 1, nspin
          spin_mag_tot(ipol,ichan,ir) = spin_mag(ig1,ipol)
        enddo
      enddo
    ENDDO
  ENDDO
!----

 END DO
!------------ ichan loop

 DEALLOCATE(ind)
 DEALLOCATE(funr)
 DEALLOCATE(fung)
 DEALLOCATE(kfunz)
 DEALLOCATE(amat)
 DEALLOCATE(vec1)
 DEALLOCATE(ipiv)
 DEALLOCATE(aux_proc)
 DEALLOCATE(spin_mag)
 IF (noncolin) DEALLOCATE(becsum_nc)
 DEALLOCATE( becsum_orig )
 DEALLOCATE( xdata, ydata, y2d )
 DEALLOCATE( xdatax, ydatax )

 RETURN
END SUBROUTINE scat_states_comp
