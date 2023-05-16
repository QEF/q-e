!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ep_matrix_element_wannier()
  !-----------------------------------------------------------------------
  !! Electron-phonon calculation from data saved in \(\texttt{fildvscf}\).
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : celldm, omega, ibrav,at
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, amass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fft_interpolate
  USE noncollin_module, ONLY : nspin_mag, noncolin
  USE dynmat, ONLY : dyn, w2
  USE modes,  ONLY : npert, nirr, u
  USE control_ph, ONLY : trans
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  USE io_global, ONLY : stdout
  USE mp_pools,  ONLY : me_pool, root_pool
  USE klist, ONLY : xk
  USE el_phon, ONLY: elph_mat, kpq, g_kpq, igqg, xk_gamma
  USE uspp,                 ONLY: okvan
  USE paw_variables, ONLY : okpaw
  USE uspp_param, ONLY : nhm
  USE lsda_mod, ONLY : nspin
  USE lrus,   ONLY : int3, int3_nc, int3_paw
  USE qpoint, ONLY : xq, nksq, ikks

  USE apply_dpot_mod,   ONLY : apply_dpot_allocate, apply_dpot_deallocate, apply_dpot_bands
  USE becmod,          ONLY : bec_type, becp, calbec, &
    allocate_bec_type, deallocate_bec_type
  USE wavefunctions, ONLY: evc
  USE eqv,             ONLY :  dvpsi
  USE el_phon, ONLY : iunwfcwann,lrwfcr
  USE wvfct, ONLY: nbnd, npwx
  USE io_global,       ONLY : ionode
  USE uspp,            ONLY : okvan, nkb, vkb
  USE klist,           ONLY : xk, ngk, igk_k
  USE lrus,            ONLY : becp1
  USE uspp_init,        ONLY : init_us_2
  USE io_files,              ONLY : prefix, diropn
  USE control_lr, ONLY : lgamma
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE dfile_star,    ONLY : drho_star, dvscf_star
  USE units_ph,        ONLY : iudrho
  USE dv_of_drho_lr,         ONLY : dv_of_drho

  USE qexsd_module,         ONLY : qexsd_readschema
  USE qexsd_copy,         ONLY : qexsd_copy_efermi, qexsd_copy_kpoints
  USE qes_types_module,   ONLY : output_type
  USE io_global,          ONLY : stdout, ionode, ionode_id
  USE mp_images,          ONLY : intra_image_comm
  USE qes_libs_module,    ONLY : qes_reset
  USE mp,                 ONLY : mp_bcast
  !
  IMPLICIT NONE

  !
  LOGICAL :: read_dvscf_cart, ascii_dvscf
  INTEGER :: irr, imode0, ipert, is, ik
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)

  type(bec_type) :: becp2 
  integer :: ipol,iundvscf_e
  INTEGER ::  ikk, ikq, ikqg
  LOGICAL :: exst,lsda_
  integer :: npw,npwq, ierr
  INTEGER, EXTERNAL :: find_free_unit

  TYPE(output_type)  :: output_obj 
  LOGICAL              :: two_fermi_energies_
  REAL(dp)             :: nelec_, ef_, ef_up_, ef_dw_

  INTEGER               :: nks_start_,nk1_, nk2_, nk3_, k1_, k2_, k3_,nkstot_
  REAL(dp), ALLOCATABLE  :: xk_start_(:,:), wk_start_(:)
  REAL(dp)              :: degauss_,homo_,lumo_
  CHARACTER(LEN=250) :: occupations_, smearing_, fname
  REAL(dp) , ALLOCATABLE  :: xk_(:,:), wk_(:)
  !REAL(dp)                 ALLOCATABLE ::  wg(:,:), et(:,:)
       !
  CALL start_clock ('elphon')

  fname = trim(dvscf_star%dir)//trim(prefix)//'.xml'
 


  INQUIRE ( file=fname, exist=exst )
  IF (exst ) THEN
      write(stdout,*)
      write(stdout,*) 'Reading xml file:  '
      write(stdout,*) fname
      write(stdout,*)
     !
     ! ... in these cases, we need (or it is useful) to read the Fermi energy
     !
     if(ionode)CALL qexsd_readschema ( fname, ierr, output_obj )
     CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF ( ierr > 0 ) CALL errore ( 'ep_matrix_element_wannier', 'problem reading  from file '//fname, ierr )
     if(ionode)CALL qexsd_copy_efermi ( output_obj%band_structure, &
       nelec_, ef_, two_fermi_energies_, ef_up_, ef_dw_ )
     if(ionode)call qexsd_copy_kpoints ( output_obj%band_structure, nks_start_, xk_start_,&
       wk_start_, nk1_, nk2_, nk3_, k1_, k2_, k3_, occupations_, smearing_, degauss_ )

     if(ionode)then
       homo_ = output_obj%band_structure%highestOccupiedLevel
       lumo_ = output_obj%band_structure%lowestUnoccupiedLevel
     end if



     CALL qes_reset  ( output_obj )
     if (occupations_.eq.'fixed')then
       smearing_ = 'none'
       degauss_ = -100
     endif

     if (ionode) call print_ph_input_param('scf_dfpt.2epiq.in',nks_start_,occupations_,&
       smearing_,degauss_,homo_,lumo_,xk_start_(:,:),wk_start_, nk1_, nk2_, nk3_, k1_, k2_, k3_,nelec_,ef_,lsda_)
   endif


  if (iudvscf.ne.0)then
    ascii_dvscf=.false.
    if(elph_mat) read_dvscf_cart=.true.
    if(read_dvscf_cart) then
      write(stdout,*)
      write(stdout,*) 'Reading dvscf in cartesian coordinates !'
      write(stdout,*)

      u=(0.d0,0.d0)
      do irr=1,3*nat
        u(irr,irr)=(1.d0,0.d0)
      enddo



    endif



    !
    ! read Delta Vscf and calculate electron-phonon coefficients
    !
    imode0 = 0
    DO irr = 1, nirr
      ALLOCATE (dvscfin (dfftp%nnr, nspin_mag , npert(irr)) )
      IF (okvan) THEN
        ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, npert(irr)))
        IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, npert(irr)))
        IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, npert(irr)))
      ENDIF

      DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
          imode0 + ipert,  -1 )
      END DO
      !     endif
      IF (doublegrid) THEN
        ALLOCATE (dvscfins (dffts%nnr, nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
          DO ipert = 1, npert(irr)
            CALL fft_interpolate (dfftp, dvscfin(:,is,ipert), dffts, dvscfins(:,is,ipert))
          ENDDO
        ENDDO
      ELSE
        dvscfins => dvscfin
      ENDIF
      CALL newdq (dvscfin, npert(irr))
      CALL elphel_refolded (npert (irr), imode0, dvscfins)
      !
      imode0 = imode0 + npert (irr)
      IF (doublegrid) DEALLOCATE (dvscfins)
      DEALLOCATE (dvscfin)

      IF (okvan) THEN
        DEALLOCATE (int3)
        IF (okpaw) DEALLOCATE (int3_paw)
        IF (noncolin) DEALLOCATE(int3_nc)
      ENDIF

    ENDDO
    !
    ! now read the eigenvalues and eigenvectors of the dynamical matrix
    ! calculated in a previous run
    !
    IF (.NOT.trans) CALL readmat_findq (iudyn, ibrav, celldm, nat, ntyp, &
      ityp, omega, amass, tau, xq, w2, dyn)
    !
  endif











  IF(lgamma) then
    write(stdout,*)
    write(stdout,*)  '     Writing momentum operator'
    write(stdout,*)  '     (estimated size:',3*nksq*nbnd*nbnd*0.059571,'kB)'


    IF((iudrho .ne. 0)) then
      allocate (dvscfins ( dfftp%nnr , nspin_mag , 3))
      dvscfins = (0.0_DP,0.0_DP)

      do ipol = 1, 3
        CALL davcio_drho ( dvscfins(:,:,ipol),  lrdrho, iudrho,  ipol,  -1 )
        call dv_of_drho ( dvscfins(:,:,ipol), .false.)
      end do

      CALL apply_dpot_allocate()
    end if

    do ik = 1, nksq
      npw = ngk(ik)
      npwq= npw     ! q=0 always in this routine
      !
      ikk = ikks(ik)
      ikqg = kpq(ik)
      npw = ngk(ikk)
      IF (lsda) current_spin = isk (ikk)
      !
      CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)

      evc=(0.d0,0.d0)
      call read_wfc_rspace_and_fwfft( evc , ik , lrwfcr , iunwfcwann , npw , igk_k(1,ikk) )


      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      do ipol = 1, 3

        dvpsi=(0.d0, 0.d0)
        call allocate_bec_type ( nkb, nbnd, becp2)
        call commutator_Hx_psi (ikk, nbnd, at(:, ipol), becp1(ik), becp2, dvpsi(:, 1:nbnd) )

        IF (nkb > 0) call deallocate_bec_type (becp2)

        call oper2epiq("comm_op",ik,ipol, dvpsi,cmplx(1,0,DP), evc) ! comm

        IF((iudrho .ne. 0)) then
          CALL apply_dpot_bands(ik, nbnd, dvscfins(:, :, ipol), evc, dvpsi)

          call adddvscf (ipol, ik)
          call oper2epiq("dvhxc_e",ik,ipol, dvpsi,cmplx(1,0,DP),evc) ! dvhxc_e
        end if

      end do!pol
    end do !ik



    IF((iudrho .ne. 0)) then
      deallocate (dvscfins)
      CALL apply_dpot_deallocate()
    end if

    call symm_dump()
  end if

  deallocate(xk_gamma)
  deallocate(kpq,g_kpq,igqg)

  ! CALL stop_clock ('elphon')
  RETURN
END SUBROUTINE ep_matrix_element_wannier

subroutine print_ph_input_param(fname,nkstot_,occupations_,&
    smearing_,degauss_,homo_,lumo_,xk_,wk_, nk1_, nk2_, nk3_, k1_, k2_, k3_,nelec_,ef_,lsda_)
  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev
  implicit none
  LOGICAL :: lsda_
  INTEGER               :: io_un, ik 
  INTEGER               :: nkstot_,ngauss_
  INTEGER               :: nk1_, nk2_, nk3_, k1_, k2_, k3_
  REAL(dp)             :: nelec_, ef_, ef_up_, ef_dw_
  REAL(dp)              :: degauss_,homo_,lumo_
  REAL(dp)                :: xk_(3,nkstot_), wk_(nkstot_)
  CHARACTER(LEN=*) :: occupations_, smearing_, fname
  INTEGER, EXTERNAL :: find_free_unit

  io_un = find_free_unit()
  open(unit=io_un,file=fname,status='unknown')
  write(io_un,'(A)') "! parameter of the SCF DFPT calculation useful for EPIq"
  write(io_un,'(A)') "&Diff_Start_Param"
  
  write(io_un,'(3x,A,f12.6,A)') "efermi=", ef_*2*rytoev,", ! in (eV)"
  write(io_un,'(3x,A,f12.6,A)') "nel_r=", nelec_,","
  if(trim(occupations_).eq.'fixed')then
    write(io_un,'(3x,A,f12.6,A,A)') "homo=", homo_*2*rytoev,", ! in (eV)"
    write(io_un,'(3x,A,f12.6,A,A)') "lumo=", lumo_*2*rytoev,", ! in (eV)"
  else
    write(io_un,'(3x,A,f12.6,A)') "sigma_ph=", degauss_*2,", ! in (Rydberg)"
    if(trim(smearing_).eq.'gauss')then
      ngauss_ = 0
    else if(trim(smearing_).eq.'mp')then
      ngauss_ = 1
    else if (trim(smearing_).eq.'fd')then
      ngauss_ = -99
    else
      ngauss_ = 66
    endif
    write(io_un,'(3x,A,I3,A)') "ngauss_ph=", ngauss_,", ! "//trim(smearing_)
  end if
  write(io_un,'(A)') "/"
  write(io_un,'(/,A)') "KPOINTS"
  if(nkstot_ .eq. 0)then 
    write(io_un,'(A)') "automatic"
    write(io_un,'(3(2x,i6),3x,3(2x,i3))')  nk1_, nk2_, nk3_, k1_, k2_, k3_
  else
    write(io_un,'(A)') "tpiba"
    write(io_un,'(6x,i9)') nkstot_
    DO ik =1,nkstot_
      IF ( lsda_) THEN
        write(io_un,'(3x,4(f12.6,2x))') xk_(:,ik),wk_(ik)
      ELSE 
        write(io_un,'(3x,4(f12.6,2x))') xk_(:,ik),wk_(ik)
      END IF
      !
    END DO
  end if
  close(io_un)

  

  end subroutine
!-----------------------------------------------------------------------
SUBROUTINE elphsum_wannier(q_index)
  !-----------------------------------------------------------------------
  !! Sum over BZ of the electron-phonon matrix elements \(\text{el_ph_mat}\).
  !
  !! Original routine written by Francesco Mauri.  
  !! Adapted to wannier functions by Matteo Calandra.  
  !! Dev. Comment: missing calc_sigma_yet
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, ityp, tau,amass,tau, ntyp => nsp, atm
  USE cell_base, ONLY : at, bg, ibrav, celldm 
  USE symm_base, ONLY : s, sr, irt, nsym, time_reversal, invs, copy_sym, inverse_s
  USE klist, ONLY : xk, nelec
  USE wvfct, ONLY : nbnd, et
  USE el_phon
  USE mp_pools,  ONLY : me_pool, root_pool, inter_pool_comm, npool
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE io_global, ONLY : stdout,ionode
  USE io_files,  ONLY : prefix
  USE dynmat, ONLY : dyn, w2
  USE modes, ONLY : u
  USE lsda_mod, only : isk,nspin, current_spin,lsda
  USE mp,        ONLY: mp_sum

  USE lr_symm_base, ONLY : irotmq, irgq, gimq, gi
  USE qpoint,     ONLY : xq, nksq, ikks, ikqs
  USE control_lr, ONLY : lgamma
  USE noncollin_module, ONLY : noncolin
  USE units_ph, ONLY : iudvscf
  !
  IMPLICIT NONE
  !
  LOGICAL :: lborn,write_ascii
  INTEGER :: q_index,iobabby
  !
  !
  logical :: minus_qloc,sym (48)
  integer :: nq, imq, isq(48)
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, ipert, jpert, nu, mu, &
    ios, iuelphmat,i,j,nt,k
  INTEGER :: iu_sym,nmodes,nsymq
  INTEGER :: io_file_unit
  !   for star_q
  REAL(DP) :: rtauloc(3,48,nat)
  !   end of star_q definitions
  real(DP) :: sxq (3, 48)
  REAL(DP) xk_dummy(3)
  character(len=80) :: filelph
  CHARACTER(len=256) ::  file_elphmat
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER (LEN=6), EXTERNAL :: int_to_char

  if(iudvscf.eq.0)return

  write_ascii=.false.
  write_ascii=.true.

  nmodes=3*nat

  write(filelph,'(A5,f9.6,A1,f9.6,A1,f9.6)') 'elph.',xq(1),'.',xq(2),'.',xq(3)
  file_elphmat=trim(adjustl(prefix))//'_elph.mat.q_'// TRIM( int_to_char( q_index ) )

  lborn=.false.
  ! parallel case: only first node writes
  IF ( .not.ionode ) THEN
    iuelphmat = 0
    iobabby = 0
  ELSE
    !
    ! First I dump information for the electron-phonon interaction
    !

    if(write_ascii)then
      iobabby = find_free_unit()
      write(filelph,'(A,A)') trim(file_elphmat),'.asci'

      OPEN (unit = iobabby, file = filelph, status = 'unknown', err = &
        1111, iostat = ios, form='formatted')
      1111  CALL errore ('elphsum_wannier', 'opening file'//"dump_elph", ABS (ios) )
    end if

    iuelphmat = find_free_unit()
    OPEN (unit = iuelphmat, file = file_elphmat, status = 'unknown', err = &
      111, iostat = ios, form='unformatted')
    111  CALL errore ('elphsum_wannier', 'opening file'//file_elphmat, ABS (ios) )
    REWIND (iuelphmat)
    xk_dummy(:)=xq(:)
    call cryst_to_cart(1,xk_dummy,at,-1)
    WRITE (iuelphmat) xk_dummy
    WRITE (iuelphmat) noncolin, nspin, lborn
    WRITE (iuelphmat) nelec
    WRITE (iuelphmat) elph_nbnd_min,elph_nbnd_max,nbnd
    WRITE (iuelphmat) nmodes, nksq, nat, ntyp
    WRITE (iuelphmat) ibrav,(celldm(j),j=1,6)
    WRITE(iuelphmat)  (atm(j),j=1,ntyp),(amass(j),j=1,ntyp), &
      (ityp(j),j=1,nat),((tau(j,i),j=1,3),i=1,nat)
    WRITE (iuelphmat) (w2 (nu) , nu = 1, nmodes)
    WRITE (iuelphmat) ((u(ipert,jpert),ipert=1,nmodes),jpert=1,nmodes)
    WRITE (iuelphmat) ((dyn(ipert,jpert),ipert=1,3*nat),jpert=1,3*nat)
    do ik=1,nksq
      ikk=ikks(ik)
      ikq=ikqs(ik) 
      xk_dummy(:)=xk(:,ikk)
      call cryst_to_cart(1,xk_dummy,at,-1)
      WRITE (iuelphmat) (xk_dummy(ipert),ipert=1,3)
      WRITE (iuelphmat) (et(ibnd,ikk),ibnd=elph_nbnd_min,elph_nbnd_max)

      do nu=1,nmodes
        WRITE (iuelphmat) &
          ((el_ph_mat (jbnd, ibnd, ik, nu),jbnd=elph_nbnd_min,elph_nbnd_max),&
          ibnd=elph_nbnd_min,elph_nbnd_max)

        if(write_ascii)then
          do ibnd=elph_nbnd_min,elph_nbnd_max
            do jbnd=elph_nbnd_min,elph_nbnd_max
              WRITE (iobabby, FMT = "(4I5,2ES30.15)") ik, nu, ibnd, jbnd, el_ph_mat (ibnd, jbnd, ik, nu)
            enddo
          enddo
        endif

      enddo
    enddo

    if(write_ascii)    close(iobabby)



    !
    ! Then I dump symmetry operations
    !
    minus_qloc = .true.
    sym = .false.
    sym(1:nsym) = .true.

    call smallg_q (xq, 0, at, bg, nsym, s, sym, minus_qloc)
    nsymq = copy_sym(nsym, sym)
    ! recompute the inverses as the order of sym.ops. has changed
    CALL inverse_s ( )

    ! part 2: this redoes most of the above, plus it computes irgq, gi, gimq
    CALL smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, &
      minus_qloc, gi, gimq)

    sym(1:nsym)=.true.
    call sgam_ph (at, bg, nsym, s, irt, tau, rtauloc, nat, sym)
    call star_q(xq, at, bg, nsym , s , invs , nq, sxq, &
      isq, imq, .FALSE. )


    do j=1,3
      write(iuelphmat) (at(i,j),i=1,3)
    enddo
    do j=1,3
      write(iuelphmat) (bg(i,j),i=1,3)
    enddo
    write(iuelphmat) nsym,nq,imq
    do i=1,nsym
      write(iuelphmat) i,invs(i),isq(i)
      do j=1,3
        do k=1,3
          write(iuelphmat) k,j, s(k,j,i)
        enddo
      enddo
      do j=1,nat
        write(iuelphmat) j, irt(i,j) 
      enddo
      do j=1,3
        do k=1,nat
          write(iuelphmat) j,i, rtauloc(j,i,k)  
        enddo
      enddo
      do j=1,3
        write(iuelphmat) j, sxq(j,i)
      enddo
    enddo

    close(iuelphmat)
  endif

  !
  !


  RETURN


END SUBROUTINE ELPHSUM_WANNIER

!
!-----------------------------------------------------------------------
SUBROUTINE elphel_refolded (npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !! Calculation of the electron-phonon matrix elements \(\text{el_ph_mat}\):
  !! $$ \langle \psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)\rangle $$
  !
  !! Original routine written by Francesco Mauri.
  !
  USE kinds, ONLY : DP
  USE fft_base, ONLY : dffts
  USE wavefunctions,  ONLY: evc
  USE io_files, ONLY: prefix, diropn
  USE klist, ONLY: xk, ngk, igk_k
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE buffers, ONLY : get_buffer
  USE wvfct, ONLY: nbnd, npwx
  USE uspp, ONLY : vkb
  USE el_phon, ONLY : el_ph_mat, iunwfcwann, igqg, kpq, g_kpq, &
    xk_gamma, npwq_refolded, lrwfcr
  USE modes, ONLY : u
  USE units_ph, ONLY : iubar, lrbar
  USE control_ph, ONLY : trans
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp_pools,  ONLY: me_pool, root_pool
  USE mp,        ONLY: mp_sum
  USE ions_base, ONLY : nat
  USE io_global, ONLY : stdout

  USE eqv,        ONLY : dvpsi!, evq
  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE control_lr, ONLY : lgamma
  USE lrus,       ONLY : becp1
  USE phus,       ONLY : alphap
  USE apply_dpot_mod,   ONLY : apply_dpot_allocate, apply_dpot_deallocate, apply_dpot_bands
  USE uspp_init,        ONLY : init_us_2

  IMPLICIT NONE
  !
  INTEGER :: npe, imode0
  COMPLEX(DP) :: dvscfins (dffts%nnr, nspin_mag, npe)
  COMPLEX(DP), allocatable :: evq(:,:)


  ! LOCAL variables
  logical :: exst
  INTEGER :: npw, npwq
  INTEGER :: nrec, ik, ikk, ikq, ikqg,ipert, mode, ibnd, jbnd, ir, ig, &
    ios

  COMPLEX(DP) , ALLOCATABLE :: elphmat (:,:,:),aux_psi(:,:)
  COMPLEX(DP), EXTERNAL :: zdotc
  INTEGER, EXTERNAL :: find_free_unit
  !
  allocate (evq(npol*npwx,nbnd))
  ALLOCATE (elphmat ( nbnd , nbnd , 3*nat))
  allocate (aux_psi(npol*npwx,nbnd))
  CALL apply_dpot_allocate()
  !
  !$acc enter data create(dvscfins(1:dffts%nnr, 1:nspin_mag, 1:npe), dvpsi(1:npwx*npol, 1:nbnd))
  !

  ! iunwfcwann=find_free_unit()
  ! CALL diropn (iunwfcwann, 'wfc', lrwfc, exst, dvscf_dir)
  ! IF (.NOT.exst) THEN
  !    CALL errore ('elphel_refolded', 'file '//trim(prefix)//'.wfc not found in Rotated_DVSCF', 1)
  ! END IF
  !
  !  Start the loops over the k-points
  !

  DO ik = 1, nksq
    !
    !  ik = counter of k-points with vector k
    !  ikk= index of k-point with vector k
    !  ikq= index of k-point with vector k+q
    !       k and k+q are alternated if q!=0, are the same if q=0
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    ikqg = kpq(ik)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    IF (lsda) current_spin = isk (ikk)
    !
    CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
    !
    ! read unperturbed wavefuctions psi(k) and psi(k+q)
    !
    evc=(0.d0,0.d0)

    !    Warning error in reading wfc, this could explain.
    !    We read here the wfc at the Gamma point, that is
    !    that saved by Wannier.

    !     CALL davcio (evc, lrwfc, iunwfcwann, ik, - 1)
    !     CALL davcio (evq, lrwfc, iunwfcwann, ikqg, - 1)



    !     IF (nksq.GT.1) THEN
    !        IF (lgamma) THEN
    !           CALL davcio (evc, lrwfc, iunwfcwann, ikk, - 1)
    !        ELSE
    !           CALL davcio (evc, lrwfc, iunwfcwann, ik, - 1)
    !           CALL davcio (evq, lrwfc, iunwfcwann, ikqg, - 1)
    !        ENDIF
    !     ENDIF
    !

    call read_wfc_rspace_and_fwfft( evc , ik , lrwfcr , iunwfcwann , npw , igk_k(1,ikk) )


    call calculate_and_apply_phase(ik, ikqg, igqg, npwq_refolded, g_kpq,xk_gamma, evq, .true.)


    DO ipert = 1, npe
      nrec = (ipert - 1) * nksq + ik
      !
      !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
      !
      IF (trans) THEN
        CALL get_buffer (dvpsi, lrbar, iubar, nrec)
      ELSE
        mode = imode0 + ipert
        ! FIXME : .false. or .true. ???
        CALL dvqpsi_us (ik, u (1, mode), .FALSE., becp1, alphap)
      ENDIF
      !
      ! calculate dvscf_q*psi_k
      !
      CALL apply_dpot_bands(ik, nbnd, dvscfins(:, :, ipert), evc, aux_psi)
      dvpsi(:,:)=dvpsi(:,:)+aux_psi(:,:)
      !
      CALL adddvscf (ipert, ik)
      !
      ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
      !
      DO ibnd =1, nbnd
        DO jbnd = 1, nbnd
          elphmat (jbnd, ibnd, ipert) = zdotc (npwq_refolded, evq (1, jbnd), 1, &
            dvpsi (1, ibnd), 1)
          IF (noncolin) &
            elphmat (jbnd, ibnd, ipert) = elphmat (jbnd, ibnd, ipert)+ &
            zdotc (npwq_refolded, evq(npwx+1,jbnd),1,dvpsi(npwx+1,ibnd), 1)
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum (elphmat, intra_bgrp_comm)
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
  !$acc exit data delete(dvscfins, dvpsi)

  !  CLOSE( UNIT = iunwfcwann, STATUS = 'KEEP' )
  !
  CALL apply_dpot_deallocate()
  DEALLOCATE (elphmat)
  DEALLOCATE(evq)
  !
  RETURN
END SUBROUTINE elphel_refolded
!
subroutine get_equivalent_kpq(xk,xq,kpq,g_kpq, igqg)
  !==================================================================
  !! Set up the \(k+q\) shell for electron-phonon coupling
  !
  !! This routine finds the G vectors such that \(k+q+G=k'\)  with 
  !! \(k\) and \(k'\) belonging to \(\text{nksq}\) for each \(k\),
  !! the G vector is stored in \(\text{g_kpq}\) and \(k'=\text{kpq}(ik)\)
  !! and finally \(\text{igqg}(ik)\) is the index that allows to find
  !! the g-vector \(\text{g_kpq}\) in the list of all the G-vectors.
  !
  !! Matteo Calandra
  !===================================================================  
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE cell_base,   ONLY : at, bg
  USE qpoint, ONLY : nksq, ikks
  USE gvect, ONLY: g, gg
  USE qpoint, ONLY : nksq
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp, ONLY : mp_sum
  ! WARNING g_kpq mesh is an integer
  implicit none

  ! Variables that are private

  integer :: iqx,iqy,iqz,i,j,k,n,nn,iq,ik, ig
  integer :: kpq(nksq),g_kpq(3,nksq),igqg(nksq)
  integer, allocatable :: ig_check(:)
  real(kind=dp) :: gg_
  real(kind=dp) :: xq(3), xk(3,*)
  real(kind=dp) :: xkpq(3),Gvec(3),xq_crys(3)
  real(kind=dp), allocatable :: xk_crys(:,:)

  !
  ! nksq = number of k point per pool withour k+q
  !

  ! 
  ! The xk_point entering here must be the k and not
  !    the k+q
  !

  xq_crys=xq

  call cryst_to_cart (1, xq_crys, at, -1)


  allocate(xk_crys(3,nksq))

  do ik=1,nksq
    xk_crys(1:3,ik)=xk(1:3,ik)
  enddo
  call cryst_to_cart (nksq, xk_crys, at, -1)

  !
  ! kpt_latt are the BZ vectors in crystalline coordinates
  ! xq is the q vector in crystalline coordinates
  !

  do iq=1,nksq
    xkpq(:)=xk_crys(:,iq)+xq_crys(:)
    do i=1,nksq
      do iqx=-4,4
        do iqy=-4,4
          do iqz=-4,4
            Gvec(1)=real(iqx,dp)+xkpq(1)
            Gvec(2)=real(iqy,dp)+xkpq(2)
            Gvec(3)=real(iqz,dp)+xkpq(3)

            if(dabs(xk_crys(1,i)-Gvec(1)).lt.1.d-6.and. &
              dabs(xk_crys(2,i)-Gvec(2)).lt.1.d-6.and. &
              dabs(xk_crys(3,i)-Gvec(3)).lt.1.d-6) then
              kpq(iq)=i
              g_kpq(1,iq)=-iqx
              g_kpq(2,iq)=-iqy
              g_kpq(3,iq)=-iqz
              goto 99
            endif
          enddo
        enddo
      enddo
    enddo
    CALL errore ('get_equivalent_kpq', 'cannot find index k+q ', 2 )
    stop
    99   continue
  enddo

  !
  ! here between all the g-vectors I find the index of that
  ! related to the translation in the Brillouin zone.
  ! Warning: if G does not belong to the processor igqg is zero.
  !

  igqg=0
  do ik=1,nksq
    Gvec(:) = REAL( g_kpq(:,ik),dp )
    call cryst_to_cart (1, Gvec, bg, 1)
    gg_ = Gvec(1)*Gvec(1) + Gvec(2)*Gvec(2) + Gvec(3)*Gvec(3)
    igqg(ik)=0
    ig=1
    do while  (gg(ig) <= gg_ + 1.d-6) 
      if ( (abs(g(1,ig)-Gvec(1)) < 1.d-6) .and.  &
        (abs(g(2,ig)-Gvec(2)) < 1.d-6) .and.  &
        (abs(g(3,ig)-Gvec(3)) < 1.d-6)  ) then
        igqg(ik) = ig

      endif
      ig= ig +1
    end do
  end do

  allocate(ig_check(nksq))
  ig_check=igqg
  CALL mp_sum( ig_check, intra_bgrp_comm )
  do ik=1,nksq
    if(ig_check(ik).eq.0) &
      CALL errore('get_equivalent_kpq', &
      ' g_kpq vector is not in the list of Gs', 100*ik )

  enddo
  deallocate(xk_crys)


end subroutine get_equivalent_kpq

subroutine calculate_and_apply_phase(ik, ikqg, igqg, npwq_refolded, g_kpq, xk_gamma, evq, lread)
  USE kinds, ONLY : DP
  USE fft_base, ONLY : dffts
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE wvfct, ONLY: nbnd, npwx
  USE gvect, ONLY : ngm, g
  USE gvecw, ONLY : gcutw
  USE cell_base, ONLY : bg
  USE qpoint, ONLY : nksq
  USE wavefunctions, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE el_phon, ONLY:iunwfcwann, lrwfcr

  IMPLICIT NONE

  LOGICAL :: lread
  INTEGER :: ik, ikqg,   npwq_refolded
  INTEGER :: igqg(nksq)
  INTEGER :: g_kpq(3,nksq)
  REAL (DP) :: xk_gamma(3,nksq)
  complex(dp) :: evq(npwx*npol,nbnd)
  !  internal
  INTEGER :: npw_, m,i
  INTEGER, allocatable :: igk_(:), igkq_(:)
  REAL(DP) :: xkqg(3), g_(3), g_scra(3,ngm)
  REAL(DP), ALLOCATABLE :: gk(:)
  COMPLEX (DP), allocatable :: psi_scratch(:)
  complex(DP), allocatable :: phase(:)

  allocate(igk_(npwx), igkq_(npwx), gk(npwx) )
  allocate (psi_scratch ( dffts%nnr) )
  allocate (phase(dffts%nnr))
  FLUSH (6)

  g_scra=g

  g_(:)=real( g_kpq(:,ik), dp )
  call cryst_to_cart (1, g_, bg, 1)
  xkqg(:)=xk_gamma(:,ikqg)+g_(:)

  npw_=0
  npwq_refolded=0
  igk_=0
  igkq_=0


  call gk_sort (xk_gamma(1,ikqg), ngm, g_scra, gcutw, npw_, igk_, gk)

  if(lread) then
    call read_wfc_rspace_and_fwfft( evq , ikqg , lrwfcr , iunwfcwann , npw_ , igk_ )
  endif

  call gk_sort (xkqg, ngm, g_scra, gcutw, npwq_refolded, igkq_, gk)

  phase(:) = (0.d0,0.d0)

  if ( igqg(ik)>0) then
    phase( dffts%nl(igqg(ik)) ) = (1.d0,0.d0)
  endif


  CALL invfft ('Wave', phase, dffts)
  !  call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
  phase(:)=conjg(phase(:))


  if(npwq_refolded.ne.npw_) call errore('calculate_and_apply_phase', 'Warning : npwq_refolded \= npw_',-1)

  do m=1,nbnd
    psi_scratch = (0.d0, 0.d0)
    psi_scratch(dffts%nl (igk_ (1:npw_) ) ) = evq (1:npw_, m)
    !     psi_scratch(dffts%nl (igk_ (1:npw) ) ) = evq (1:npw, m)
    CALL invfft ('Wave', psi_scratch, dffts)
    !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
    psi_scratch(1:dffts%nnr) = psi_scratch(1:dffts%nnr) * phase(1:dffts%nnr)
    !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
    CALL fwfft ('Wave', psi_scratch, dffts)
    evq(1:npwq_refolded,m) = psi_scratch(dffts%nl (igkq_(1:npwq_refolded) ) )
  enddo

  if(noncolin) then
    do m=1,nbnd
      psi_scratch = (0.d0, 0.d0)
      psi_scratch(dffts%nl (igk_ (1:npw_) ) ) = evq (npwx+1:npwx+npw_, m)
      !       psi_scratch(dffts%nl (igk_ (1:npw) ) ) = evq (1:npw, m)
      CALL invfft ('Wave', psi_scratch, dffts)
      !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
      psi_scratch(1:dffts%nnr) = psi_scratch(1:dffts%nnr) * phase(1:dffts%nnr)
      !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
      CALL fwfft ('Wave', psi_scratch, dffts)
      !       evq(npwx+1:npwx+npwq_refolded,m) = psi_scratch(dffts%nl (igkq_(1:npwq_refolded) ) )
      evq((npwx+1):(npwx+npwq_refolded),m) = psi_scratch(dffts%nl (igkq_(1:npwq_refolded) ) )
    enddo
  endif

  deallocate(psi_scratch)
  DEALLOCATE(phase)
  deallocate(gk, igk_, igkq_)

  return
end subroutine calculate_and_apply_phase


!-----------------------------------------------------------------------
SUBROUTINE readmat_findq (iudyn, ibrav, celldm, nat, ntyp, ityp, omega, &
    amass, tau, q, w2, dyn)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry
  USE control_ph, ONLY : xmldyn
  USE output, ONLY : fildyn
  USE io_dyn_mat,  ONLY : read_dyn_mat_param, read_dyn_mat_header, &
    read_dyn_mat, read_dyn_mat_tail
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
  INTEGER :: nspin_mag, nqs
  REAL(DP) :: at(3,3)
  REAL(DP) :: bg(3,3)
  REAL(DP) :: m_loc(3,nat)
  INTEGER :: ityp__ (nat)
  REAL(DP) :: amass__ (ntyp)
  INTEGER :: iq
  REAL(DP) :: xq(3)
  COMPLEX(DP) :: u(3*nat,3*nat)
  REAL(DP) :: dynr (2, 3, nat, 3, nat), err_q(3)
  COMPLEX(DP) :: dynr_c(3,3,nat,nat)

  CHARACTER(len=80) :: line
  CHARACTER(len=3)  :: atm
  INTEGER :: nt, na, nb, naa, nbb, nu, mu, i, j
  LOGICAL :: lfound

  !
  !
  IF(xmldyn) THEN
    CALL read_dyn_mat_param(fildyn, ntyp_, nat_ )
    CALL read_dyn_mat_header(ntyp_, nat_, ibrav_, nspin_mag,  &
      celldm_, at, bg, omega, atm, amass__, tau_, ityp__, m_loc, &
      nqs )
    IF ( ntyp.NE.ntyp_ .OR. nat.NE.nat_ .OR.ibrav_.NE.ibrav .OR. &
      ABS ( celldm_ (1) - celldm (1) ) > 1.0d-5) &
      CALL errore ('readmat', 'inconsistent data a', 1)
    DO nt = 1, ntyp
      IF ( ABS (amass__ (nt) - amass (nt) ) > 1.0d-5) &
        CALL errore ( 'readmat', 'inconsistent data  b', 1 + nt)
    ENDDO
    DO na = 1, nat
      IF (ityp__ (na).NE.ityp (na) ) CALL errore ('readmat', &
        'inconsistent data c',  na)
    ENDDO

  ELSE
    REWIND (iudyn)
    READ (iudyn, '(a)') line
    READ (iudyn, '(a)') line
    READ (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
    IF ( ntyp.NE.ntyp_ .OR. nat.NE.nat_ .OR.ibrav_.NE.ibrav .OR. &
      ABS ( celldm_ (1) - celldm (1) ) > 1.0d-5) &
      CALL errore ('readmat', 'inconsistent data', 1)
  IF ( ibrav_ == 0 ) THEN
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
  END IF
    DO nt = 1, ntyp
      READ (iudyn, * ) i, atm, amass_
      IF ( nt.NE.i .OR. ABS (amass_ - amu_ry*amass (nt) ) > 1.0d-5) then
        !CALL errore ( 'readmat', 'inconsistent data mass', 1 + nt)
        write (6, *)'readmat inconsistent data mass using read ones'
        amass (nt) = amass_/amu_ry
        write (6, *)'amass',nt,amass (nt)*amu_ry
      endif
    ENDDO
    DO na = 1, nat
      READ (iudyn, * ) i, ityp_, tau_
      IF (na.NE.i.OR.ityp_.NE.ityp (na) ) CALL errore ('readmat', &
        'inconsistent data', 10 + na)
    ENDDO

  ENDIF

  lfound=.false.
  iq=0

  do while(.not.lfound)

    IF(xmldyn) THEN

      iq = iq+1
      CALL read_dyn_mat(nat,iq,xq,dynr_c)
      !     CALL read_dyn_mat_tail(nat,omega,u)
      err_q(1:3)=dabs(xq(1:3)-q(1:3))

      if(err_q(1).lt.1.d-7.and.err_q(2).lt.1.d-7.and.err_q(3).lt.1.d-7) lfound=.true.

      DO nb = 1, nat
        DO j = 1, 3
          DO na = 1, nat
            DO i = 1, 3
              dynr (1, i, na, j, nb) = REAL(dynr_c(i, j, na, nb))
              dynr (2, i, na, j, nb) = AIMAG(dynr_c(i, j, na, nb))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ELSE
      READ (iudyn, '(a)') line
      READ (iudyn, '(a)') line
      READ (iudyn, '(a)') line
      READ (iudyn, '(a)') line

      READ (line (11:80), * ) (q_ (i), i = 1, 3)

      err_q(1:3)=dabs(q_(1:3)-q(1:3))

      if(err_q(1).lt.1.d-7.and.err_q(2).lt.1.d-7.and.err_q(3).lt.1.d-7) lfound=.true.



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

    ENDIF

    if(lfound) then
      !
      ! divide the dynamical matrix by the (input) masses (in amu)
      !
      DO nb = 1, nat
        DO j = 1, 3
          DO na = 1, nat
            DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                ityp (na) ) * amass (ityp (nb) ) ) / amu_ry
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                ityp (na) ) * amass (ityp (nb) ) ) / amu_ry
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
          dyn (mu, nu) = dyn (mu, nu) / SQRT ( amu_ry * amass (ityp (na) ) )
        ENDDO
      ENDDO
      !
      !
    endif
  enddo


  RETURN
END SUBROUTINE readmat_findq





subroutine symm_dump()
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, ityp, tau,amass,tau, ntyp => nsp, atm
  USE cell_base, ONLY : at, bg, ibrav, celldm 
  USE symm_base, ONLY : s, sr, irt, nsym, time_reversal, invs, copy_sym, inverse_s
  USE klist, ONLY : xk, nelec
  USE wvfct, ONLY : nbnd, et
  !USE el_phon
  USE mp_pools,  ONLY : me_pool, root_pool, inter_pool_comm, npool
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE io_global, ONLY : stdout,ionode
  USE io_files,  ONLY : prefix
  !USE dynmat, ONLY : dyn, w2
  !USE modes, ONLY : u
  USE lsda_mod, only : isk,nspin, current_spin,lsda
  USE mp,        ONLY: mp_sum

  USE lr_symm_base, ONLY : irotmq, irgq, gimq, gi
  USE qpoint,     ONLY : xq, nksq, ikks, ikqs
  USE control_lr, ONLY : lgamma
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  LOGICAL :: lborn
  INTEGER :: q_index
  !
  !
  logical :: minus_qloc,sym (48)
  integer :: nq, imq, isq(48)
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, ipert, jpert, nu, mu, &
    ios, iuelphmat,i,j,nt,k
  INTEGER :: iu_sym,nmodes,nsymq
  INTEGER :: io_file_unit
  !   for star_q
  REAL(DP) :: rtauloc(3,48,nat)
  !   end of star_q definitions
  real(DP) :: sxq (3, 48)
  !REAL(DP) xk_dummy(3)
  !character(len=80) :: filelph
  CHARACTER(len=256) ::  file_elphmat
  !
  !COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER (LEN=6), EXTERNAL :: int_to_char


  !
  ! Then I dump symmetry operations
  !
  minus_qloc = .true.
  sym = .false.
  sym(1:nsym) = .true.

  !write(iuelphmat, *) 
  !write(iuelphmat, *) xq, 0,  nsym, s, sym, minus_qloc

  call smallg_q (xq, 0, at, bg, nsym, s, sym, minus_qloc)
  nsymq = copy_sym(nsym, sym)
  ! recompute the inverses as the order of sym.ops. has changed
  CALL inverse_s ( )

  ! part 2: this redoes most of the above, plus it computes irgq, gi, gimq
  CALL smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, &
    minus_qloc, gi, gimq)

  sym(1:nsym)=.true.
  call sgam_ph (at, bg, nsym, s, irt, tau, rtauloc, nat, sym)
  call star_q(xq, at, bg, nsym , s , invs , nq, sxq, &
    isq, imq, .FALSE. )

  file_elphmat= 'symm_dump.dat'

  iuelphmat = find_free_unit()

  OPEN (unit = iuelphmat, file = file_elphmat, status = 'unknown', err = &
    111, iostat = ios, FORM='FORMATTED')
  !111, iostat = ios, form='unformatted')


  111  CALL errore ('symm_dump', 'opening file'//file_elphmat, ABS (ios) )
  do j=1,3
    !write(iuelphmat) (at(i,j),i=1,3)
    write(iuelphmat, '(3E16.8)') (at(i,j),i=1,3)
  enddo
  do j=1,3
    !write(iuelphmat) (bg(i,j),i=1,3)
    write(iuelphmat, '(3E16.8)') (bg(i,j),i=1,3)
  enddo


  write(iuelphmat,*) nsym,nq,imq
  do i=1,nsym
    write(iuelphmat,*) i,invs(i),isq(i)
    do j=1,3
      do k=1,3
        write(iuelphmat,*) k,j, s(k,j,i)
      enddo
    enddo
    do j=1,nat
      write(iuelphmat,*) j, irt(i,j) 
    enddo
    do j=1,3
      do k=1,nat
        write(iuelphmat,*) j,i, rtauloc(j,i,k)  
      enddo
    enddo
    do j=1,3
      write(iuelphmat,*) j, sxq(j,i)
    enddo
  enddo

  write(iuelphmat,*) 'irotmq, minus_qloc'
  write(iuelphmat,*) irotmq, minus_qloc

  close(iuelphmat)



  file_elphmat= 'symm_dump.bin'

  iuelphmat = find_free_unit()

  OPEN (unit = iuelphmat, file = file_elphmat, status = 'unknown', err = &
    111, iostat = ios, form='unformatted')

  do j=1,3
    write(iuelphmat) (at(i,j),i=1,3)
  enddo
  do j=1,3
    write(iuelphmat) (bg(i,j),i=1,3)
  enddo
  write(iuelphmat) nsym,nq,imq
  do i=1,nsym
    write(iuelphmat) i,invs(i),isq(i)
    do j=1,3
      do k=1,3
        write(iuelphmat) k,j, s(k,j,i)
      enddo
    enddo
    do j=1,nat
      write(iuelphmat) j, irt(i,j) 
    enddo
    do j=1,3
      do k=1,nat
        write(iuelphmat) j,i, rtauloc(j,i,k)  
      enddo
    enddo
    do j=1,3
      write(iuelphmat) j, sxq(j,i)
    enddo
  enddo

  close(iuelphmat)

end subroutine symm_dump




!-----------------------------------------------------------------------
SUBROUTINE oper2epiq(fname,ik,ipol, dpsi,coeff,evc_in)
  !-----------------------------------------------------------------------
  !
  !      extract and write the operator from dpsi projecting in onto the PW basis (evc)
  !

  USE kinds,          ONLY  : DP
  !USE wavefunctions,  ONLY  : evc
  USE wvfct,          ONLY  : npwx, nbnd
  USE noncollin_module,ONLY :  npol
  USE io_global,       ONLY : stdout, ionode
  USE mp_world,        ONLY : mpime, world_comm
  USE mp_pools,        ONLY : me_pool, my_pool_id, npool
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,             ONLY  : mp_sum,mp_rank, mp_size
  USE io_global,       ONLY : stdout
  USE qpoint,          ONLY : nksq
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE buffers,          ONLY : get_buffer


  USE klist, ONLY: xk, ngk, igk_k
  USE qpoint,     ONLY : ikks, ikqs
  USE io_files,              ONLY : prefix, diropn
  USE fft_base,              ONLY : dfftp, dffts
  USE noncollin_module,      ONLY :  npol 
  USE dfile_star,    ONLY : drho_star, dvscf_star

  implicit none
  CHARACTER (LEN=*),INTENT(IN)      :: fname 
  CHARACTER(LEN=256)                  :: str , file_name_peid
  INTEGER     ,INTENT(IN)             :: ipol, ik !, comm
  !INTEGER     ,INTENT(IN),OPTIONAL   :: iunwfcr
  COMPLEX(DP)                         ::     evc_in(npol*npwx,nbnd)
  COMPLEX(DP) ,INTENT(IN)         :: coeff  
  COMPLEX(DP) ,INTENT(IN)    :: dpsi(npwx*npol,nbnd)

  INTEGER                             :: i, j, ios, pe_id, file_unit
  COMPLEX(DP), ALLOCATABLE            :: oper(:,:)  
  COMPLEX(DP), POINTER           :: evc(:,:)

  !COMPLEX(DP), ALLOCATABLE            :: gathered_oper(:,:,:)  

  INTEGER                             :: iunwfcwann, ikk, npw,lrwfcr
  logical                         :: exst

  allocate(oper(nbnd,nbnd))

  !call get_buffer (evc, lrwfc, iuwfc, ik)



  !if (.not.present(evc_in))then
  !ALLOCATE(evc(npol * npwx, nbnd))
  !if (present(iunwfcr))then
  !iunwfcwann=iunwfcr
  !else
  !iunwfcwann=188
  !if (ionode)then
  !lrwfcr= 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x *npol
  !CALL diropn (iunwfcwann, 'wfc_r', lrwfcr, exst, dvscf_star%dir)
  !IF (.NOT.exst) THEN
  !CALL errore ('oper2epiq', 'file .wfc_r not found in '//dvscf_star%dir, 1)
  !END IF
  !END IF
  !end if

  !!
  !! read unperturbed wavefuctions psi(k) and psi(k+q)
  !!
  !evc=(0.d0,0.d0)

  !ikk = ikks(ik)
  !npw = ngk(ikk)

  !call read_wfc_rspace_and_fwfft( evc , ik , lrwfcr , iunwfcwann , npw , igk_k(1,ikk) )!lrwfcr
  !end if

  !call get_buffer (evc, lrwfc, iuwfc, ik)


  !oper = (0_dp,0_dp)
  !if (present(evc_in))then
  CALL zgemm( 'C', 'N', nbnd, nbnd , npwx, &
    (1.d0,0.d0), evc_in, npwx, dpsi, npwx, &
    (0.d0,0.d0), oper, nbnd )
  !else
  !CALL zgemm( 'C', 'N', nbnd, nbnd , npwx, &
  !(1.d0,0.d0), evc, npwx, dpsi, npwx, &
  !(0.d0,0.d0), oper, nbnd )
  !end if
  !subroutine zgemm 	( 	
  !   character  	TRANSA,
  !		character  	TRANSB,
  !		integer  	M,
  !		integer  	N,
  !		integer  	K,
  !		complex*16  	ALPHA,
  !		complex*16, dimension(lda,*)  	A,
  !		integer  	LDA,
  !		complex*16, dimension(ldb,*)  	B,
  !		integer  	LDB,
  !		complex*16  	BETA,
  !		complex*16, dimension(ldc,*)  	C,
  !		integer  	LDC 
  !	) 		
  !C := alpha*op( A )*op( B ) + beta*C,

  ! TRANSA = 'N' or 'n',  op( A ) = A.
  ! TRANSA = 'T' or 't',  op( A ) = A**T.
  ! TRANSA = 'C' or 'c',  op( A ) = A**H.
  ! BETA is COMPLEX*16
  !  On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !  supplied as zero then C need not be set on input.

  CALL mp_sum(oper,intra_bgrp_comm)

  pe_id = my_pool_id + 1

  if (npool .eq. 1) then 
    if (ionode) then 
      write(str , '(i10)') pe_id
      file_name_peid = adjustl(trim(fname))//trim(adjustl(trim(str)))//adjustl(trim('.dat'))
      file_unit = 88 + pe_id !io_file_unit()
      ios=0
      open (file_unit, FILE=file_name_peid, position="append",  IOSTAT=ios)
      if(ios.ne.0)  write(stdout,*) 'problem opening file for printing operator ', file_name_peid

      DO i=1,nbnd !size(oper,1)
        DO j =1 , nbnd !size (oper,2)
          WRITE(UNIT=file_unit, FMT = "(4I5,2ES20.10)") ik  , ipol , i , j , coeff*oper(i,j)
        ENDDO
      ENDDO

      close (file_unit)
    end if
  else
    write(str , '(i10)') pe_id
    file_name_peid = adjustl(trim(fname))//trim(adjustl(trim(str)))//adjustl(trim('.dat'))
    file_unit = 88 + pe_id !io_file_unit()
    ios=0
    open (file_unit, FILE=file_name_peid, position="append",  IOSTAT=ios)
    if(ios.ne.0)  write(stdout,*) 'problem opening file for printing operator ', file_name_peid

    DO i=1,nbnd !size(oper,1)
      DO j =1 , nbnd !size (oper,2)
        !WRITE(UNIT=file_unit, FMT = "(4I5,2ES40.26)") ik + (pe_id-1)*nksq , ipol , i , j , coeff*oper(i,j)
        WRITE(UNIT=file_unit, FMT = "(4I5,2ES20.10)") ik + (pe_id-1)*nksq , ipol , i , j , coeff*oper(i,j)
      ENDDO
    ENDDO
    close (file_unit)
  end if

  deallocate(oper)


  !if (.not.present(iunwfcr).and..not.present(evc_in))&
  !CLOSE( UNIT = iunwfcwann, STATUS = 'KEEP' )

  !------ if all k are already computed-------------------- 

  !rank = MPI_Comm_rank(MPI_COMM_WORLD);
  !mpi_size = MPI_Comm_size(MPI_COMM_WORLD);

  !allocate(oper(nks_pp[rank],nbnd,nbnd))
  !if (ionode) then 
  !  allocate(gathered_oper(nks,nbnd,nbnd))
  !endif

  !i_tmp = 0 
  !do i=0,mpi_size
  !  i_tmp = i_tmp + nks_pp[rank]
  !  receive_counts[i] = i_tmp 
  !  receive_displacements[i] = i_tmp 
  !enddo

  !MPI_Gatherv(sendarray, nks_pp[rank], MPI_DOUBLE_COMPLEX, &
  !  oper, receive_counts, receive_displacements, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ios)

  !ios=0
  !call mpi_barrier(MPI_COMM_WORLD,ios)
  !if(ios.ne.0)  write(stdout,*) 'problem in mpi_barrier '     

  !------ if all k are already computed-------------------- 
END SUBROUTINE oper2epiq


SUBROUTINE oper2epiq_header(fname)
  USE cell_base,       ONLY : tpiba, at, bg
  USE io_global,       ONLY : stdout, ionode
  implicit none
  CHARACTER (LEN=*),INTENT(IN)      :: fname 
  CHARACTER(LEN=256)                  :: str , file_name_peid
  INTEGER                             :: i, j, ios, pe_id, file_unit



  if (ionode) then 
    pe_id = 0
    write(str , '(i10)') pe_id
    file_name_peid = adjustl(trim(fname))//trim(adjustl(trim(str)))//adjustl(trim('.dat'))
    file_unit = 88 + pe_id !io_file_unit()
    ios=0
    !open (file_unit, FILE=file_name_peid, position="append",  IOSTAT=ios)
    open (file_unit, FILE=file_name_peid,   IOSTAT=ios)
    if(ios.ne.0)  write(stdout,*) 'problem opening file for printing operator ', file_name_peid


    WRITE(UNIT=file_unit, FMT = "(a,ES20.10)") 'tpiba', tpiba

    WRITE(UNIT=file_unit, FMT = "(a)") 'at'
    do i=1,3
      WRITE(UNIT=file_unit, FMT = "(3ES20.10)") at(i,1:3)
    end do 

    WRITE(UNIT=file_unit, FMT = "(a)") 'bg'
    do i=1,3
      WRITE(UNIT=file_unit, FMT = "(3ES20.10)") bg(i,1:3)
    end do 
    close (file_unit)
  end if

END SUBROUTINE oper2epiq_header



SUBROUTINE write_zeu2epiq(zstarue)
  !export zeu for epiq

  USE kinds,      ONLY : DP
  USE io_files,   ONLY : prefix
  USE klist,      ONLY : degauss, ngauss, nelec, lgauss
  USE cell_base,        ONLY : celldm
  USE ions_base,  ONLY : nat
  USE ener,       ONLY : ef         
  USE io_global,  ONLY : stdout, ionode
  implicit none

  integer :: ipol, jpol, icart, jcart, na, nu, mu, ierr,i,j, ios, iunit
  INTEGER, EXTERNAL :: find_free_unit
  REAL (DP)              :: zstarue(3, nat, 3) 
  real(DP) :: ehomo, elumo

  if(ionode)then
    iunit=find_free_unit()
    open(unit=iunit,file=trim(adjustl(prefix))//'.zeu.2epik',status='unknown',form='formatted',IOSTAT=ios)
    if(ios.ne.0) then
      write(stdout,*) 'ERROR reading ', trim(adjustl(prefix))//'.zeu.2epik'
    end if

    CALL get_homo_lumo (ehomo, elumo)

    write(iunit,*) '# nat celldm(1) efermi (Ryd) nelec'
    if(lgauss)then
      write(iunit,*) nat, celldm(1), ef, nelec
    else
      write(iunit,*) nat, celldm(1), ehomo, nelec
    end if
    write(iunit,*) '# sigma ngauss  omega eta'
    write(iunit,*) degauss, ngauss,0._dp, 0._dp
    !write(iunit,*) degauss, ngauss,omegaLuca, etaLuca
    write(iunit,*) '# Born effective charges'
    do na = 1, nat
      do ipol = 1, 3
        write (iunit, '(3(e24.12,"  0.0 "))')  (zstarue (ipol, na, jpol) , jpol = 1, 3)
        !write (iunit, '(6e24.12)')  (zstarue_cmplx (ipol, na, jpol) , jpol = 1, 3)
      enddo
    enddo

    close(iunit)
  end if


END SUBROUTINE write_zeu2epiq


!   CALL elph_scdft_gather_r(gg,3 * nat * nbnd_fs * nbnd_fs * nksq2,gg_col,nrcv, &
!   &                me_image, nproc_image, intra_image_comm)
! 
! SUBROUTINE elph_scdft_gather_r(snd,nsnd,rcv,nrcv,mype,npe,comm)
!   !----------------------------------------------------------------------
!   !
!   ! This routine gathers a real matrix to PE 0.
!   !
!   USE kinds, ONLY : dp
!   USE mp, ONLY : mp_sum, mp_gather
!   !
!   INTEGER,INTENT(IN) :: nsnd, nrcv, mype, npe, comm
!   REAL(dp),INTENT(IN) :: snd(nsnd)
!   REAL(dp),INTENT(OUT) :: rcv(nrcv)
!   !
!   INTEGER :: cnt(0:npe - 1), dsp(0:npe - 1), ipe
!   !
!   cnt(0:npe - 1) = 0
!   cnt(mype) = nsnd
!   !
!   CALL mp_sum(cnt, comm)
!   !
!   dsp(0) = 0
!   DO ipe = 1, npe - 1
!      dsp(ipe) = dsp(ipe - 1) + cnt(ipe - 1)
!   END DO
!   !
!   CALL mp_gather(snd(1:nsnd), rcv(1:nrcv), cnt, dsp, 0, comm)
!   !
! END SUBROUTINE elph_scdft_gather_r
! 



