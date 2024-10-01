!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kcw_setup
  !-----------------------------------------------------------------------
  !
  !!  This subroutine prepares several variables which are needed in the
  !!  KC calculation:
  !!  * computes the total local potential (external+scf) on the smooth
  !!    grid to be used in h_psi and similia
  !!  * set the inverse of every matrix invs
  !!  * for metals sets the occupied bands
  !!  * Open buffer for the KS and eventualy the minimizing WFCs
  !!  * Open buffer for the KS states in the WANNIER gauge
  !!  * Read the U matrix from Wannier
  !!  * Rotate the KS state to the localized gauge
  !!  * Compute the periodic part of the wannier orbital and store in the buffer 
  !!    a separate directory is generated for each q 
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization, lsda, isk
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE gvecs,             ONLY : doublegrid, ngms
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
  USE gvect,             ONLY : ig_l2g, mill
  USE wvfct,             ONLY : nbnd
  !USE funct,             ONLY : dft_is_gradient
  USE xc_lib,            ONLY : xclib_dft_is
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx
  USE control_flags,     ONLY : io_level, gamma_only
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer, get_buffer
  USE control_kcw,       ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, kcw_iverbosity, lgamma_iq, &
                                spin_component, isq, read_unitary_matrix, x_q, tmp_dir_save, & 
                                num_wann, num_wann_occ, occ_mat, tmp_dir_kcw, tmp_dir_kcwq, &
                                io_sp, io_real_space, nrho, nkstot_eff!, wq, nqstot
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : nkstot, xk, nelec, nelup, neldw
  USE cell_base,         ONLY : at, omega, bg, tpiba
  USE fft_base,          ONLY : dffts
  !
  USE control_lr,       ONLY : nbnd_occ, lgamma
  USE scf,              ONLY : rho
  USE io_global,        ONLY : ionode, ionode_id
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_files,         ONLY : create_directory
  USE io_rho_xml,       ONLY : write_scf
  !
  USE io_kcw,           ONLY : write_rhowann, write_rhowann_g
  !
  USE mp,               ONLY : mp_sum
  USE control_lr,       ONLY : lrpa
  USE coulomb,          ONLY : setup_coulomb
  !
  USE mp_pools,         ONLY : my_pool_id
  USE mp_bands,         ONLY : my_bgrp_id, root_bgrp_id, &
                               root_bgrp, intra_bgrp_comm, &
                               inter_bgrp_comm
  USE symm_base,        ONLY : s, time_reversal, sr, ft
  USE control_kcw,      ONLY : s_w, nqstot_ibz, xq_ibz, wq_ibz, nsym_w
  USE start_k,          ONLY : nk1, nk2, nk3, k1, k2, k3
  USE parameters,       ONLY : npk

  !
  implicit none

  integer :: na, i, ik, ip, ipp
  ! counters
  !
  INTEGER   :: lrwfc, iun_qlist!, nkstot_eff
  LOGICAL   :: exst, exst_mem
  INTEGER :: iq, nqs
  REAL(DP) :: xq(3)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:), rhowann_aux(:)
  CHARACTER (LEN=256) :: filename, file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER:: mp1, mp2, mp3
  INTEGER :: nsymwi
  INTEGER :: s_wi(3,3,48)
  REAL(DP) :: Gvector(3)
  INTEGER :: iq_
  INTEGER, ALLOCATABLE :: ibz2fbz(:)

  !
  INTEGER :: &
       igk_k_all(npwx,nkstot),&    ! index of G corresponding to a given index of k+G
       ngk_all(nkstot)             ! number of plane waves for each k point
  !
  ! Auxiliary variables for SH calculation
  COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:), sh(:),rhor(:,:), rhog_aux(:)
  ! The periodic part of the wannier orbital density
  COMPLEX(DP), ALLOCATABLE  :: rhowann_g(:,:,:)
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
  !
  ! The weight of each q point
  REAL(DP), ALLOCATABLE :: weight(:)
  LOGICAL :: lrpa_save
  REAL(DP) :: xq_(3)
  INTEGER  :: lrrho
  COMPLEX(DP) :: struct_fact, int_wann, int_rho
  COMPLEX(DP), ALLOCATABLE :: rho_c(:,:,:),wann_c(:,:,:)
  COMPLEX(DP) :: phase(dffts%nnr)
  INTEGER :: iwann
  !
  ALLOCATE ( delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
  !WRITE(*,*), 'NGMS', ngms
  !
  call start_clock ('kcw_setup')
  !
  IF (nspin == 4) THEN
    nkstot_eff = nkstot
    nrho = 4
  ELSE
    nkstot_eff = nkstot/nspin
    nrho = 1
  ENDIF 
  ALLOCATE(rhor(dffts%nnr,nrho), rhog(ngms,nrho), rhog_aux(ngms))
  !rho_c
  ! ... Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin_mag, doublegrid)
  !
  !  ... If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     ux=0.0_DP
     IF ( xclib_dft_is('gradient') ) CALL compute_ux(m_loc,ux,nat)
     !if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! ... Computes the number of occupied bands for each k point
  !
  !call setup_nbnd_occ ( ) 
  WRITE(stdout, '(/,5X, "REPORT # of electrons")')
  WRITE(stdout, '(  5X, "nelec= ", F15.8)') nelec
  IF (nspin == 2) THEN 
    WRITE(stdout, '(  5X, "nelup= ", F15.8)') nelup
    WRITE(stdout, '(  5X, "neldw= ", F15.8)') neldw
  ELSE IF( nspin ==1 )  THEN
    WRITE(stdout, '(  5X, "nelup= ", F15.8)') nelec/2
    WRITE(stdout, '(  5X, "neldw= ", F15.8)') nelec/2
  ENDIF
  WRITE(stdout, '(  5X, "nkstot= ", I5)') nkstot
  WRITE(stdout, '(  5X, "nspin=  ", I5)') nspin
  !WRITE(stdout, '(  5X, "nbnd_occ(up)= ", I5)') nbnd_occ(1)
  !WRITE(stdout, '(  5X, "nbnd_occ(dw)= ", I5)') nbnd_occ(nkstot/nspin+1)

  !
  ! ... Open buffers for the KS and eventualy the minimizing WFCs 
  ! 
  iuwfc = 20; lrwfc = nbnd * npwx * npol; io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  !
  IF (.NOT.exst.AND..NOT.exst_mem) &
       CALL errore ('kcw_setup', 'file '//trim(prefix)//'.wfc not found', 1)
  !
  if (kcw_iverbosity .gt. 1) &
       WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! ... READ the U matrix from Wannier and set the toal number of WFs
  !
  IF (read_unitary_matrix) THEN 
    !
    CALL read_wannier ( )
    if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
    !
  ELSE
    !
    num_wann = nbnd
    num_wann_occ = nint (nelec)/2 ! spin upolarized
    IF (nspin_mag == 4) THEN 
      num_wann_occ = nint (nelec) 
    ELSE IF (nspin==2) THEN 
      num_wann_occ = nint (nelup) ! Assuming insulating
      IF (spin_component == 2) num_wann_occ = nint (neldw) ! Assuming insulating
    END IF
    !
  ENDIF
  ! 
  !  ... Allocate relevant quantities ...
  ! 
  ALLOCATE (rhowann ( dffts%nnr, num_wann, nrho), rhowann_aux(dffts%nnr) ) ! was * nrho
  ALLOCATE ( evc0(npwx*npol, num_wann) )
  ALLOCATE ( rhowann_g (ngms, num_wann, nrho) )
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )
  ALLOCATE ( sh(num_wann) )
  ALLOCATE ( rho_c(dffts%nnr,num_wann,nrho) )
  ALLOCATE ( wann_c(dffts%nnr,num_wann,nrho) )
  sh(:) = CMPLX(0.D0,0.D0,kind=DP)
  occ_mat = 0.D0
  !
  ! ... Open a new buffer to store the KS states in the WANNIER gauge
  !
  iuwfc_wann = 21
  io_level = 1
  lrwfc = num_wann * npwx * npol
  CALL open_buffer ( iuwfc_wann, 'wfc_wann', lrwfc, io_level, exst )
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs, OPENED")')
  !
  ! ... Open an other buffer for the KS states in the WANNIER gauge which contains
  !     all the k points (not just the one in this pool). This is needed for each k-point 
  !     to have access to all the other k-points. MEMORY INTENSE
  !
  iuwfc_wann_allk = 210
  io_level = 1
  lrwfc = num_wann * npwx * npol
  CALL open_buffer ( iuwfc_wann_allk, 'wfc_wann_allk', lrwfc, io_level, exst )
  !
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs ALL-k, OPENED")')
  !
  ! ... Rotate the KS state to the localized gauge nd save on a buffer
  !
  CALL rotate_ks () 
  !
  ! ... pass all the WFs to all the pool (needed to have pool parallelization)
  !
  CALL bcast_wfc ( igk_k_all, ngk_all )
  !
  !DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
  call setup_coulomb()
  !
  ! 8)
  !CALL compute_map_ikq ()
  !
  ! ... Compute the periodic part of the wannier orbital and store it on file 
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: PREPARING THE KCW CALCULATION ...")')
  !
  !     Write the list of q points on a file and store the q coordinates
  !
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir_kcw)//'qlist.txt')
  !
  ! 
  nqs = nkstot_eff
  ALLOCATE (x_q (3, nqs), isq(nqs) )
  ALLOCATE (weight(nqs) )
  iq=1
  IF (ionode) THEN 
     WRITE(iun_qlist,'(i5)') nkstot_eff
     DO ik = 1, nkstot
       IF (lsda .AND. isk(ik) /= spin_component) CYCLE
       WRITE(iun_qlist, '(3f12.8)') xk(:,ik)
       x_q(:,iq) = xk(:,ik)  
       isq(iq) = isk(ik) 
       iq = iq + 1
     ENDDO 
  ENDIF
  CALL mp_bcast (x_q, ionode_id, intra_image_comm)
  CALL mp_bcast (isq, ionode_id, intra_image_comm)
  !
  ALLOCATE ( lgamma_iq(nqs) )
  lgamma_iq(:) = .FALSE.
  !
  WRITE( stdout, '(/, 5X,"INFO: Compute Wannier-orbital Densities ...")')
  !
  rho_c = ZERO
  wann_c = ZERO
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    xq = x_q(:,iq)
    !
    ! IF (ionode) WRITE(iun_qlist,'(3f12.8)') xq
    !
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    WRITE( stdout, '(/,8X, 78("="))')
    WRITE( stdout, '(  8X, "iq = ", i5)') iq
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cart ]")') x_q(:,iq)
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cryst]")') xq(:)
    WRITE( stdout, '(  8X, 78("="),/)')
    !
    CALL compute_map_ikq_single (iq)
    ! The map to identify which k point in the 1BZ corresponds to k+q and the G vector that produce the mapping
    ! The results are stored in the global variable map_ikq and shift_1bz (used inside rho_of_q) 
    ! can (should) be moved inside rho_of_q ( )
    !
    rhowann(:,:,:)=ZERO
    !! Initialize the periodic part of the wannier orbtal density at this q
    !
    CALL rho_of_q (rhowann, ngk_all, igk_k_all)
    ! Compute the peridic part rho_q(r) of the wannier density rho(r)
    ! rho(r)   = \sum_q exp[iqr]rho_q(r)
    !
    WRITE( stdout, '(8X,"INFO: rho_q(r) DONE ",/)')
    !
    ! Compute the Self Hartree
    weight(iq) = 1.D0/nqs ! No SYMM 
    lrpa_save=lrpa
    lrpa = .true.

    !!! Start code to check sum over q for wannier objects
    !tmp_dir = tmp_dir_save  ! the periodic part are written on the original outdir 
    !lrrho=num_wann*dffts%nnr
    !CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    !tmp_dir = tmp_dir_kcwq   ! go back to the q-specific directory
    !
    xq_(:) = -x_q(:,iq)
    ! calculate_phase has a - sign inside
    phase=ZERO
    CALL calculate_phase(xq_, phase)
    CALL structure_factor(iq, struct_fact)
    CALL cryst_to_cart(1, xq_, at, -1)
    IF (kcw_iverbosity .gt. 1) WRITE(stdout,'(8X, "INFO: iq = ", i5, 3x, "Structure Factor S(q) [Re, Im] = ", 2f12.8,/)') & 
                                                                iq, struct_fact
    !
     DO iwann = 1, num_wann
       !
       DO ip=1,nrho
         wann_c(:,iwann,ip) = wann_c(:,iwann,ip) + phase(:)*rhowann(:,iwann,ip)*weight(iq)
         rho_c(:,iwann,ip)  = rho_c(:,iwann,ip)  + phase(:)*rhowann(:,iwann,ip)*struct_fact*weight(iq)
         !WRITE(*,*), phase(:)*rhowann(:,iwann,ip)*weight(iq)
       ENDDO 
       !
    ENDDO    
    !!! End code to check sum over q for wannier objects
    !
    DO i = 1, num_wann
      !
      rhog(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
      delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
      vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
      rhor(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
      !
      rhor(:,:) = rhowann(:,i,:) 
      !! The periodic part of the orbital density in real space
      !
      CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
      rhowann_g(:,i,:) = rhog
      !! The periodic part of the perturbation DeltaV_q(G)
      ! 
      sh(i) = sh(i) + 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:) )*weight(iq)*omega
      !
    ENDDO
    !
    lrpa=lrpa_save
    !
    ! ... each q /= gamma is saved on a different directory
    lgamma = lgamma_iq(iq)
    !
    tmp_dir_kcwq= TRIM (tmp_dir_kcw) //'q' &
                  & // TRIM(int_to_char(iq))//'/'
    filename=TRIM(tmp_dir_kcwq)//'charge-density.dat'
    IF (ionode) inquire (file =TRIM(filename), exist = exst)
    !
    CALL mp_bcast( exst, ionode_id, intra_image_comm )
    !
    CALL create_directory( tmp_dir_kcwq )
    tmp_dir=tmp_dir_kcwq
    CALL write_scf( rho, nspin )
    ! write the periodic part of the wannier orbital density on file
    DO i = 1, num_wann
      !
      IF ( .NOT. io_real_space) THEN 
        !
        DO ip = 1, nrho
          rhog_aux = rhowann_g(:,i, ip)
          file_base=TRIM(tmp_dir_kcwq)//'rhowann_g_iwann_'//TRIM(int_to_char((i-1)*nrho+ip))
          IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
               CALL write_rhowann_g( file_base, &
               root_bgrp, intra_bgrp_comm, &
               bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
               gamma_only, mill, ig_l2g, rhog_aux(:), .FALSE. )
        ENDDO
        !
      ELSE  
        !
        DO ip = 1, nrho
          rhowann_aux (:) = rhowann(:,i,ip)
          file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char((i-1)*nrho+ip))
          CALL write_rhowann( file_base, rhowann_aux, dffts, ionode, inter_bgrp_comm )
        ENDDO
        !
      ENDIF
      !
    ENDDO
    tmp_dir=tmp_dir_save
    !
  ENDDO
  !
  IF (kcw_iverbosity .gt. 1) THEN
  !  write(*,'(/,"DEBUG")')
    WRITE(stdout,'(5X, "INFO: Wannier density number")')
    DO iwann = 1, num_wann
       DO ip = 1, nrho
         int_rho = (0.D0,0.D0)
         int_rho = SUM (rho_c(:,iwann,ip))/(dffts%nr1*dffts%nr2*dffts%nr3)
         int_wann = (0.D0,0.D0)
         int_wann = SUM (wann_c(:,iwann,ip))/(dffts%nr1*dffts%nr2*dffts%nr3)
         CALL mp_sum( int_rho, intra_bgrp_comm )
         CALL mp_sum( int_wann, intra_bgrp_comm )
         IF (ip == 1) THEN
            WRITE(stdout,'(5X, "iwann= ", i3, 3x, "ipol= ", i3, 3x, "int rho_wann[ipol](r) [Re, Im] =", 2f18.8)') &
                   iwann, ip, int_rho
            !WRITE(stdout,'(5X, "iwann= ", i3, 3x, "ipol= ", i3, 3x, "int rho_wann[ipol](r) [Im]     =", 2f18.8)') &
            !      iwann, ip, AIMAG(int_wann)
         ELSE
            WRITE(stdout,'(5X, 13X,                "ipol= ", i3, 3x, "int rho_wann[ipol](r) [Re, Im] =", 2f18.8)') &
                   ip, int_rho
            !WRITE(stdout,'(5X,13X,                "ipol= ", i3, 3x, "int rho_wann[ipol](r) [Im]     =", 2f18.8)') &
            !      ip, AIMAG(int_wann)
         ENDIF
       ENDDO
    ENDDO
  ENDIF
  !
  ! Print on output the self-Hatree
  CALL mp_sum (sh, intra_bgrp_comm)
  !
  OPEN (128, file = TRIM(tmp_dir_kcw)//'sh.txt')
  WRITE(stdout,'(/, 5X, "INFO: Orbital Self-Hartree (SH)")') 
  DO i = 1, num_wann
    WRITE(stdout,'(5X, "orb ", 1i5, 5X, "SH ", 1F10.6)') i, REAL(sh(i))
    WRITE(128,*) sh(i)
  ENDDO
  CLOSE (128)
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  WRITE(stdout, '(/,5X,"INFO: PREPARING THE KCW CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  !
  ! For symmetry
  !
  CALL symmetries_of_wannier_function()

  mp1=2
  mp2=2
  mp3=2
  i=6
  ALLOCATE(xq_ibz(3,mp1*mp2*mp3))
  ALLOCATE(wq_ibz(mp1*mp2*mp3))
  WRITE(*,*) "nsym_w(i)", nsym_w(i), "time_reversal", time_reversal
  nsymwi = nsym_w(i)
  s_wi(:,:,:) = s_w(:,:,:,i)
  CALL kpoint_grid ( nsymwi, time_reversal, .false., s_wi, 0, bg, &
                     mp1*mp2*mp3, 0,0,0, mp1,mp2,mp3, nqstot_ibz, xq_ibz, wq_ibz )


                     WRITE(*,*)
                     WRITE(*,*)
                     WRITE(*,*)
                     WRITE(*,*)
                     WRITE(*,*)
                     WRITE(*,*)
                     WRITE(*,*)
                 
WRITE(*,*) "nqstot_ibz", nqstot_ibz
!DO iq=1, nsymwi
!  WRITE(*,*) "isym", iq, s_wi(:,:,iq)
!END DO

WRITE(*,*) "xq_ibz"
DO iq=1, nqstot_ibz
  WRITE(*,*) xq_ibz(:, iq)
END DO                   

WRITE(*,*) "wq_ibz"
DO iq=1, nqstot_ibz
  WRITE(*,*) wq_ibz(iq)
END DO         

WRITE(*,*) "xk"
DO ik=1, mp1*mp2*mp3
  WRITE(*,*) xk(:,ik)
END DO

DO ik=nkstot+1, nkstot+nqstot_ibz
  xk(:, ik) = xq_ibz(:, ik-nkstot)
END DO                   

ALLOCATE(ibz2fbz(nqstot_ibz))
WRITE(*,*) "         iq        ik         Gvector"
DO iq=nkstot+1, nkstot+nqstot_ibz
  CALL rotate_xk(iq, 1, i, Gvector)
  ibz2fbz(iq-nkstot) = i
  WRITE(*,*) iq-nkstot, i, Gvector
END DO

sh=0
DO iq_=1, nqstot_ibz
  DO i = 6,6
    iq = ibz2fbz(iq_)
    xq = x_q(:,iq)
    !
    ! IF (ionode) WRITE(iun_qlist,'(3f12.8)') xq
    !
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    !
    CALL compute_map_ikq_single (iq)
    ! The map to identify which k point in the 1BZ corresponds to k+q and the G vector that produce the mapping
    ! The results are stored in the global variable map_ikq and shift_1bz (used inside rho_of_q) 
    ! can (should) be moved inside rho_of_q ( )
    !
    rhowann(:,:)=ZERO
    CALL rho_of_q (rhowann, ngk_all, igk_k_all)
    ! Compute the peridic part rho_q(r) of the wannier density rho(r)
    ! rho(r)   = \sum_q exp[iqr]rho_q(r)
    !
    WRITE( stdout, '(8X,"INFO: rho_q(r) DONE ",/)')
    !
    ! Compute the Self Hartree
!    weight(iq) = 1.D0/nqs ! No SYMM 
    lrpa_save=lrpa
    lrpa = .true.
    !
    rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
    delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
    vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
    rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
    !
    rhor(:) = rhowann(:,i) 
    !! The periodic part of the orbital desity in real space
    !
    CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
    rhowann_g(:,i) = rhog
    !! The periodic part of the perturbation DeltaV_q(G)
    ! 
    sh(i) = sh(i) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:) )*wq_ibz(iq_)*omega
    !
  ENDDO
END DO

WRITE(stdout,'(5X, "INFO: Orbital Self-Hartree (SH)")') 
DO i = 1, num_wann
  WRITE(stdout,'(5X, "orb ", 1i5, 5X, "SH ", 1F10.6)') i, REAL(sh(i))
END DO



  !CALL compute_dmn()
  !
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kcw_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux )
  DEALLOCATE (rhowann_g)
  DEALLOCATE (rhog , delta_vg, vh_rhog, delta_vg_, rhog_aux )
  DEALLOCATE (sh, rho_c, wann_c) 
  !
  RETURN
  !
END SUBROUTINE kcw_setup
