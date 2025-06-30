!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kcw_setup_screen
  !-----------------------------------------------------------------------
  !
  !! As kcw_setup.f90 plus screening specific setups
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
  USE fft_interfaces,    ONLY : invfft
  USE gvecs,             ONLY : doublegrid, ngms
  USE gvect,             ONLY : ig_l2g
  USE wvfct,             ONLY : nbnd
  USE xc_lib,            ONLY : xclib_dft_is
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx
  USE control_flags,     ONLY : io_level, gamma_only
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer
  USE control_kcw,       ONLY : alpha_final, iurho_wann, kcw_iverbosity, io_real_space, &
                                read_unitary_matrix, num_wann, num_wann_occ, i_orb, iorb_start, &
                                iorb_end, nqstot, occ_mat, l_do_alpha, group_alpha, &
                                tmp_dir_kcw, tmp_dir_kcwq, x_q, lgamma_iq, io_sp, irr_bz, nrho, spin_component
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : xk, nkstot, nelec, nelup, neldw
  USE cell_base,         ONLY : at, omega !, bg
  USE fft_base,          ONLY : dffts
  !
  USE control_lr,        ONLY : nbnd_occ
  USE mp,                ONLY : mp_bcast
  USE io_kcw,            ONLY : read_rhowann, read_rhowann_g
  !
  USE coulomb,           ONLY : setup_coulomb
  !
  USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: na, i, ip
  ! counters
  !
  INTEGER   :: lrwfc, lrrho, iun_qlist
  ! Lenght record sor wfc and density, iunit for the q point file (written by wann2kcw.x)
  !
  LOGICAL   :: exst, exst_mem
  ! Check on the existence of the buffers
  !
  INTEGER :: iq, nqs
  ! Counters on the q points, total number of q points
  !
  REAL(DP) :: xq(3)
  ! the q-point coordinatew
  !
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:), rhowann_aux(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  ! the periodic part of the wannier orbital density
  !
  CHARACTER (LEN=256) :: file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: iorb_from_input
  !
  call start_clock ('kcw_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin_mag, doublegrid)
  !
  !  3) If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  !
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
  ENDIF
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  !call setup_nbnd_occ ( ) 
  !
  ALLOCATE ( alpha_final(nbnd) )
  alpha_final(:) = 1.D0
  !
  ! ... Open buffers for the KS and eventualy the minimizing WFCs 
  !
  iuwfc = 20
  lrwfc = nbnd * npwx * npol
  io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('kcw_setup_screen', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! 8) READ the U matrix from Wannier and set the toal number of WFs
  !
  IF (read_unitary_matrix) THEN 
    CALL read_wannier ( )
    if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
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
  ! ... Open a buffer for the wannier orbital densities. Those have been written by wann2kcw
  ! ... and must be in the outdir. If not STOP
  !
  iurho_wann = 22
  io_level = 1
  lrrho=num_wann*dffts%nnr*nrho
  CALL open_buffer ( iurho_wann, 'rho_wann', lrrho, io_level, exst )
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WF rho, OPENED")')
  !
  ALLOCATE (rhowann ( dffts%nnr, num_wann, nrho), rhowann_aux(dffts%nnr) )
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )  !<--- nkstot_eff??
  ALLOCATE ( rhog (ngms) )
  ALLOCATE (l_do_alpha(num_wann), group_alpha(num_wann) ) 
  !
  l_do_alpha = .TRUE.
  DO i = 1, num_wann; group_alpha(i)=i; ENDDO
  occ_mat = 0.D0
  !
  ! ... Set up the coulomb kernel. If l_vcut=.true. the Gygi Balderschi scheme is used.
  ! ... Otherwise the g=0 component is set to zero 
  !
  call setup_coulomb()
  ! OLD implementation (only isotropic eps) for reference
  !call setup_coulomb_exx()
  !
  !DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
  ! ... Read the q-point grid written by wann2kcw
  !
  WRITE( stdout, '(/, 5X,"INFO: READING Wannier-orbital Densities ...")') 
  !
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir_kcw)//'qlist.txt')
  !
  READ(iun_qlist,'(i5)') nqs
  nqstot = nqs
  !
  ALLOCATE (x_q (3, nqs) ) 
  ALLOCATE ( lgamma_iq(nqs) )
  lgamma_iq(:) = .FALSE.
  !
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    READ (iun_qlist,'(3f12.8)') xq
    x_q(:,iq) = xq ! Store the q vectors
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    !
    WRITE( stdout,'(/,8X, 78("="))')
    WRITE( stdout, '( 8X, "iq = ", i5)') iq
    WRITE( stdout, '( 8X,"The  Wannier density at  q = ",3F12.7, "  [Cart ]")')  xk(:,iq)
    WRITE( stdout, '( 8X,"The  Wannier density at  q = ",3F12.7, "  [Cryst]")')  xq(:)
    WRITE( stdout, '( 8X, 78("="),/)')
    ! 
    tmp_dir_kcwq= TRIM (tmp_dir_kcw) // 'q' &
                & // TRIM(int_to_char(iq))//'/'
    !
    DO i = 1, num_wann
      !
      IF ( .NOT. io_real_space ) THEN 
        !
        DO ip = 1, nrho
          file_base=TRIM(tmp_dir_kcwq)//'rhowann_g_iwann_'//TRIM(int_to_char((i-1)*nrho+ip))
          CALL read_rhowann_g( file_base, &
               root_bgrp, intra_bgrp_comm, &
               ig_l2g, 1, rhog(:), .FALSE., gamma_only )
          rhowann_aux=(0.d0,0.d0)
          rhowann_aux(dffts%nl(:)) = rhog(:)
          CALL invfft ('Rho', rhowann_aux, dffts)
          rhowann(:,i,ip) = rhowann_aux(:)*omega
        ENDDO
        !
      ELSE 
        !
        DO ip = 1, nrho
          file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char((i-1)*nrho+ip))
          CALL read_rhowann( file_base, dffts, rhowann_aux )
          rhowann(:,i,ip) = rhowann_aux(:)
        ENDDO
        !
      ENDIF
      !
    ENDDO
    !
    ! ... Save the rho_q on a direct access file
    !
    lrrho=num_wann*dffts%nnr*nrho
    CALL save_buffer (rhowann, lrrho, iurho_wann, iq)
    !
  ENDDO
  !
  !read qlist_ibz
  ! 
  IF(irr_bz) CALL read_qlist_ibz()
  !
  ! ... Which orbital to compute
  !
  iorb_from_input = .FALSE. 
  !
  IF (i_orb == -1 ) THEN 
    !
    iorb_start = 1
    iorb_end = num_wann
    !
  ELSE IF  (i_orb .gt. num_wann ) THEN 
    !
    CALL errore('kcw_setup','i_orb > num_wann' ,ABS(num_wann))
    !
  ELSE
    !
    iorb_start = i_orb
    iorb_end = i_orb
    iorb_from_input = .TRUE. 
    !
  ENDIF
  !
  IF (iorb_from_input) WRITE( stdout, '(5X,"INFO: Orbital to compute specified from input", i5)') i_orb
  WRITE( stdout, '(5X,"INFO: total number of wannier to compute", i5, " from ", i5, " to ", i5)')  &
  iorb_end-iorb_start+1, iorb_start, iorb_end
  !
  ! ... Group the orbitals as a function of the Self-Hartree 
  !
  IF( .NOT. iorb_from_input) CALL group_orbitals ( )
  !
  WRITE( stdout, '(/,5X,"INFO: PREPARING THE KCW CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  ! 
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kcw_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux)
  DEALLOCATE (rhog) 
  !
  RETURN
  !
END SUBROUTINE kcw_setup_screen
