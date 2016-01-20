
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE dfile_star
!----------------------------------------------------------------------------
!
  USE kinds, ONLY : DP

  TYPE open_star_descriptor
     LOGICAL            :: open
     LOGICAL            :: pat
     CHARACTER(len=256) :: dir
     CHARACTER(len=10)  :: basis
     CHARACTER(len=256) :: ext
  END TYPE open_star_descriptor

  ! NOTE: default values are set in phq_readin
  TYPE(open_star_descriptor) :: &
      drho_star, & ! 
      dvscf_star   !

  TYPE rotated_pattern_repr
    !
    INTEGER, ALLOCATABLE ::  npert (:), irgq (:)
    ! output: the dimension of each represe
    ! output: the small group of q
    INTEGER :: nsymq, irotmq, nirr, npertx
    ! output: the order of the small group
    ! output: the symmetry sending q -> -q+
    ! output: the number of irr. representa
    ! output: the max number of irreps
    !
    REAL(DP), ALLOCATABLE :: gi (:,:), gimq (:), eigen(:)
    ! output: [S(irotq)*q - q]
    ! output: [S(irotmq)*q + q]
    ! output: eigenvalues of the dynmat
    !
    COMPLEX(DP), ALLOCATABLE :: u(:,:), t(:,:,:,:), tmq (:,:,:)
    ! output: the pattern vectors
    ! output: the symmetry matrices
    ! output: the matrice sending q -> -q+G
    LOGICAL :: minus_q
    ! output: if true one symmetry send q -> -q+G
    INTEGER, ALLOCATABLE :: num_rap_mode(:)
    CHARACTER(len=15), ALLOCATABLE :: name_rap_mode(:)
    ! output: the number of the representation of each mode
    ! output: the name of the representation for each group of modes
  END TYPE rotated_pattern_repr

CONTAINS
SUBROUTINE allocate_rotated_pattern_repr(rpat, nat, npertx)
  TYPE(rotated_pattern_repr),INTENT(inout) :: rpat
  INTEGER,INTENT(in) :: nat, npertx
  !
  ALLOCATE(rpat%npert(3*nat))
  ALLOCATE(rpat%irgq(48))
  ALLOCATE(rpat%gi(3,48))
  ALLOCATE(rpat%gimq(3))
  ALLOCATE(rpat%eigen(3*nat))
  ALLOCATE(rpat%u(3*nat, 3*nat))
  ALLOCATE(rpat%t(npertx, npertx, 48, 3*nat))
  ALLOCATE(rpat%tmq(npertx, npertx, 3*nat))
  ALLOCATE(rpat%num_rap_mode(3*nat))
  ALLOCATE(rpat%name_rap_mode(3*nat))
  !
END SUBROUTINE allocate_rotated_pattern_repr
SUBROUTINE deallocate_rotated_pattern_repr(rpat)
  TYPE(rotated_pattern_repr),INTENT(inout) :: rpat
  !
  DEALLOCATE(rpat%npert)
  DEALLOCATE(rpat%irgq)
  DEALLOCATE(rpat%gi)
  DEALLOCATE(rpat%gimq)
  DEALLOCATE(rpat%eigen)
  DEALLOCATE(rpat%u)
  DEALLOCATE(rpat%t)
  DEALLOCATE(rpat%tmq)
  DEALLOCATE(rpat%num_rap_mode)
  DEALLOCATE(rpat%name_rap_mode)
  !
END SUBROUTINE deallocate_rotated_pattern_repr

!-----------------------------------------------------------------------
SUBROUTINE write_dfile_star(descr, source, nsym, xq, u, nq, sxq, isq, s, &
            sr, invs, irt, ntyp, ityp, dfile_minus_q, iq_ )
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in dfile_rot
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : at, bg
  USE ions_base,        ONLY : nat, tau, amass
  USE symm_base,        ONLY : ftau, t_rev
  USE lsda_mod,         ONLY : nspin
  USE modes,            ONLY : nirr, npert, npertx
  USE units_ph,         ONLY : lrdrho
  USE io_global,        ONLY : stdout , ionode, ionode_id
  use io_files,         ONLY : diropn, prefix
  USE constants,        ONLY : tpi
  USE dfile_autoname,   ONLY : dfile_name
  USE control_ph,       ONLY : search_sym
  USE noncollin_module, ONLY : nspin_mag
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE wrappers,         ONLY : f_mkdir_safe

  USE lr_symm_base, ONLY : rtau
  !
  IMPLICIT NONE
  ! input variables:
  TYPE(open_star_descriptor),INTENT(inout) :: descr
  ! what to do, and where to put it
  CHARACTER(len=*),INTENT(in) :: source
  ! the frile where the quantity to rotate is currently stored
  INTEGER,INTENT(in) :: nsym, nq, isq(48)
  ! number of symmetry operations
  ! number of q in the star
  ! symmetry op. giving the rotated q
  INTEGER          :: iq_
  INTEGER,INTENT(in) :: s(3, 3, 48), invs(48), irt(48, nat)
  ! symmetry matrices and their inverses
  REAL(DP) :: sr(3,3,48)
  ! symmetry matrices in cartesian coordinates
  REAL(DP),INTENT(in) :: xq(3), sxq (3, 48)
  ! corrent q-point at which drho has been caclulated
  ! list of the q in the star
  COMPLEX(DP),INTENT(in) :: u(3*nat, 3*nat)
  ! the modes of the starting drho
  !
  INTEGER,INTENT(in) :: ntyp, ityp(nat)
  LOGICAL,INTENT(in) :: dfile_minus_q
  ! if .true. also use time reversal to save drho(-q) = conjg(drho(q))

  ! local variables
  INTEGER :: na, i, j
  INTEGER :: isym, nsymrot, iudfile_rot, iudfile
  INTEGER, EXTERNAL :: find_free_unit
  ! auxiliary xq\cdot\tau and \xq_s\cdot\tau
  REAL(DP) :: xq_tau,sxq_tau
  !
  INTEGER :: irr, imode0, ipert, is,k,n,nn,ri,rj,rk,isym_inv
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), ALLOCATABLE :: dfile_at(:,:,:), dfile_rot(:,:,:), dfile_rot_scr(:,:,:)
  LOGICAL :: exst
  CHARACTER(LEN=256) :: dfile_rot_name
  COMPLEX(DP) :: phase_xq
  INTEGER     :: ipol,iq,index0,nar
  INTEGER     :: ichosen_sym(48)
  COMPLEX(DP), ALLOCATABLE :: phase_sxq(:)
  ! fake vars for cartesian "patterns"
  TYPE(rotated_pattern_repr) :: rpat
  ! functions:
  CHARACTER(len=256),EXTERNAL :: trimcheck
  !
  IF ( .not. descr%open ) RETURN
  IF (descr%ext(1:5) /= 'auto:') descr%ext = 'auto:'//descr%ext
  !
  IF(nsym==1) &
    CALL errore('write_dfile_star', 'this subroutine produces random garbage without symmetry!', 1)
  !
  !
  ! create a directory to store the files
  IF (TRIM(descr%dir) ==' ') CALL errore('dfile_star', 'directory not specified', 1)
  ! the next line is not needed in phonon, but may be needed if this code is reused
  descr%dir = trimcheck(descr%dir)
  !
  IF (ionode) INQUIRE(file=trimcheck(descr%dir)//'.', exist = exst)
  CALL mp_bcast( exst, ionode_id, intra_image_comm )
  !if(.not.exst) CALL create_directory(descr%dir)
  if(.not.exst) is = f_mkdir_safe(descr%dir)
  !
  ! ionode does all the work from here on, the other nodes are aly required for
  ! calling set_irr which includes a mp broadcast
  ONLY_IONODE_1 : IF (ionode) THEN
  !
  !  Between all the possible symmetries I chose the first one
  !        (all of them lead to the same rotated dwhatever)
  !
  DO iq=1,nq
     nsymrot=0
     DO isym=1,nsym
        IF (isq(isym) == iq) THEN
           nsymrot=nsymrot+1
           IF (nsymrot == 1) ichosen_sym(iq)=isym
        ENDIF
     ENDDO
     if(nsymrot == 0) THEN
        call errore('dfile_star','no symmetry relates q at star(q)',iq)
     ENDIF
     !
  ENDDO
  !
  ALLOCATE(     dfile_at(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin, 3*nat))
  ALLOCATE(    dfile_rot(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin, 3*nat))
  ALLOCATE(dfile_rot_scr(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin, 3*nat))
  !
  dfile_at = (0._dp,0._dp)
  !
  ! Open the drho file for reading
  iudfile = 90334 !find_free_unit()
  CALL diropn(iudfile, source, lrdrho, exst)
  !
  imode0 = 0
  DO irr = 1, nirr
    ! read in drho for all the irreps
    DO is = 1, nspin
      DO ipert = 1, npert(irr)
        CALL davcio( dfile_at(:,:,imode0+ipert), lrdrho, iudfile, imode0 + ipert, -1 )
      END DO
    ENDDO
    !
    imode0 = imode0 + npert(irr)
    !
  ENDDO
  CLOSE(iudfile) 
  !
  ! Transform from the basis of the patterns to cartesian basis
  dfile_rot = (0._dp,0._dp)
  DO i=1,3*nat
     DO j=1,3*nat
        dfile_rot(:,:,i) = dfile_rot(:,:,i) + CONJG(u(i,j))*dfile_at(:,:,j)
     ENDDO
  ENDDO
  !
  ! Transform to crystalline coordinates (necessary in order to apply s)
  dfile_at =(0._dp,0._dp)
  DO i = 1,nat
     na=(i-1)*3
     DO j=1,3
        dfile_at(:,:,na+j)=dfile_rot(:,:,na+1)*at(1,j) + &
                           dfile_rot(:,:,na+2)*at(2,j) + &
                           dfile_rot(:,:,na+3)*at(3,j)
     ENDDO
  ENDDO
  !
  ! take away the phase due to the q-point
  dfile_rot = (0._dp,0._dp)
  DO i = 1,nat
    !
    xq_tau=tpi*SUM(xq*tau(:,i))
    phase_xq= CMPLX (cos(xq_tau),sin(xq_tau))
    !
    DO ipol=1,3
       imode0 = (i-1)*3 + ipol
       dfile_rot(:,:,imode0) = phase_xq*dfile_at(:,:,imode0)
    ENDDO
  ENDDO
  !
  dfile_at=dfile_rot
  !
  ! Now I rotate the dvscf
  !
  ALLOCATE(phase_sxq(nat))
  !
  ENDIF ONLY_IONODE_1 
  !
  ! This part has to be done by all cpus because some parts of rpat% are used by set_irr which 
  ! calls mp_bcast.
  CALL allocate_rotated_pattern_repr(rpat, nat, npertx)
  ! 
  Q_IN_THE_STAR : &
  DO iq=1,nq
    ONLY_IONODE_2 : IF (ionode) THEN
    dfile_rot = (0._dp,0._dp)
    !
    ! note that below isym is S and isym_inv refers to S^-1
    isym=ichosen_sym(iq)
    isym_inv=invs(ichosen_sym(iq))
    !
    DO k=1,nat  
      sxq_tau=(sxq(1,iq)*tau(1,k)+ &
            sxq(2,iq)*tau(2,k)+ &
            sxq(3,iq)*tau(3,k))*tpi
      phase_sxq(k)=1._dp/CMPLX(cos(sxq_tau),sin(sxq_tau))
    ENDDO
    !
    DO is=1,nspin
      KLOOP : DO k = 1, dfftp%nr3
        JLOOP : DO j = 1, dfftp%nr2
          ILOOP : DO i = 1, dfftp%nr1
            !
            ! Here I rotate r
            !
            CALL ruotaijk(s(1,1,isym_inv), ftau(1,isym_inv), i, j, k, &
                          dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
            !
            n  = (i-1)  + (j-1)*dfftp%nr1  + (k-1)*dfftp%nr2*dfftp%nr1  + 1
            nn = (ri-1) + (rj-1)*dfftp%nr1 + (rk-1)*dfftp%nr2*dfftp%nr1 + 1
            !
            DO na=1,nat
              nar=irt(isym_inv,na)
              index0=(nar-1)*3
              !
              DO ipol=1,3
                  imode0=(na-1)*3+ipol
                  !
                  dfile_rot(n,is,imode0) = dfile_rot(n,is,imode0) + &
                      ( s(ipol, 1, isym_inv) * dfile_at(nn,is,index0+1) + &
                        s(ipol, 2, isym_inv) * dfile_at(nn,is,index0+2) + &
                        s(ipol, 3, isym_inv) * dfile_at(nn,is,index0+3) )
                  !
              ENDDO
            ENDDO
            !
          ENDDO ILOOP
        ENDDO JLOOP
      ENDDO KLOOP
      !
    ENDDO
    !
    ! Add back the phase factor for the new q-point
    !
    DO na=1,nat
      DO ipol=1,3
        imode0=(na-1)*3+ipol
        dfile_rot_scr(:,:,imode0 )=dfile_rot(:,:,imode0)*phase_sxq(na)
      ENDDO
    ENDDO
    !
    ! Back to cartesian coordinates
    !
    dfile_rot=CMPLX(0._dp,0._dp)
    DO i=1,nat
      imode0=(i-1)*3
      DO j=1,3
        dfile_rot(:,:,imode0+j)=dfile_rot_scr(:,:,imode0+1)*bg(j,1) +&
              dfile_rot_scr(:,:,imode0+2)*bg(j,2) + dfile_rot_scr(:,:,imode0+3)*bg(j,3)
      ENDDO
    ENDDO
    !
    !
    ENDIF ONLY_IONODE_2
    !
    ! This part has to be done on all nodes because set_irr calls mp_bcast!
    ! NOTE: the new set_irr_new subroutine woudl not work here as it uses global variables!!
    IF (descr%basis=='modes') THEN
      !
      ! Transform to the basis of the patterns at the new q...
      !
      CALL set_irr (nat, at, bg, xq, s, sr, tau, ntyp, ityp, ftau, invs, nsym, &
                    rtau, irt, rpat%irgq, rpat%nsymq, rpat%minus_q, rpat%irotmq, rpat%u, rpat%npert,   &
                    rpat%nirr, rpat%gi, rpat%gimq, 0, .false., rpat%eigen, search_sym,&
                    nspin_mag, t_rev, amass, rpat%num_rap_mode, rpat%name_rap_mode)
      !
      ONLY_IONODE_2b : IF (ionode) THEN
      dfile_rot_scr = (0._dp, 0._dp)
      DO i=1,3*nat
        DO j=1,3*nat
            dfile_rot_scr(:,:,i) = dfile_rot_scr(:,:,i) + rpat%u(j,i)*dfile_rot(:,:,j)
        ENDDO
      ENDDO
      dfile_rot = dfile_rot_scr
      ENDIF ONLY_IONODE_2b
    !
    ELSE IF (descr%basis=='cartesian') THEN
      !
      ! ...or leave in the basis of cartesian displacements
      !
      rpat%nirr = 3*nat
      rpat%u = (0._dp,0._dp)
      DO i = 1,3*nat
        rpat%u(i,i) = (1._dp, 0._dp)
      ENDDO
      rpat%npert = 0
      rpat%npert(1:nirr) = 1
    ELSE
      CALL errore('dfile_star', 'basis can only be "modes" or "cartesian"', 3)
    ENDIF
    !
    !
    !  Opening files and writing
    !
    ONLY_IONODE_3 : IF (ionode) THEN
    !
    dfile_rot_name = dfile_name(sxq(:,iq), at, TRIM(descr%ext), &
         TRIM(descr%dir)//prefix, generate=.true., index_q=iq_ )
       
    iudfile_rot = find_free_unit()
    CALL diropn (iudfile_rot, TRIM(dfile_rot_name), lrdrho, exst, descr%dir)
    WRITE(stdout, '(7x,a,3f10.6,3a)') "Writing drho for q = (",sxq(:,iq),') on file "',&
                          TRIM(dfile_rot_name),'"'
    !
    DO na=1,nat
      DO ipol=1,3
        imode0=(na-1)*3+ipol
            CALL davcio( dfile_rot(:,:,imode0), lrdrho, iudfile_rot, imode0, + 1 )
      ENDDO
    ENDDO
    !
    IF(descr%pat) CALL io_pattern(nat, dfile_rot_name, rpat%nirr, rpat%npert, &
                                  rpat%u, sxq(:,iq), descr%dir, +1)
    !
    CLOSE(iudfile_rot)
    !
    ! Also store drho(-q) if necessary
    MINUS_Q : &
    IF (dfile_minus_q .and. xq(1)**2+xq(2)**2+xq(3)**2 > 1.d-5 )  THEN
      !
      dfile_rot_name = dfile_name(-sxq(:,iq), at, TRIM(descr%ext), &
                        TRIM(descr%dir)//prefix, generate=.true., index_q=iq_)
      !
      iudfile_rot = find_free_unit()
      CALL diropn (iudfile_rot, TRIM(dfile_rot_name), lrdrho, exst, descr%dir)
      WRITE(stdout, '(7x,a,3f10.6,3a)') "Writing drho for q = (",-sxq(:,iq),') on file "',&
                            TRIM(dfile_rot_name),'"'
      !
      DO na=1,nat
        DO ipol=1,3
          imode0=(na-1)*3+ipol
              CALL davcio( CONJG(dfile_rot(:,:,imode0)), lrdrho, iudfile_rot, imode0, + 1 )
        ENDDO
      ENDDO
      !
      IF(descr%pat) CALL io_pattern(nat, dfile_rot_name, rpat%nirr, rpat%npert, &
                                    CONJG(rpat%u), -sxq(:,iq), descr%dir, +1)
      !
      CLOSE(iudfile_rot)
    ENDIF &
    MINUS_Q
    !
    ENDIF ONLY_IONODE_3
    !
  ENDDO &
  Q_IN_THE_STAR
  !
  IF (ionode) THEN
    DEALLOCATE(dfile_rot, dfile_rot_scr, dfile_at)
    DEALLOCATE(phase_sxq)
  ENDIF
  CALL deallocate_rotated_pattern_repr(rpat)
  !
  RETURN
  !----------------------------------------------------------------------------
END SUBROUTINE write_dfile_star
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
END MODULE dfile_star
!----------------------------------------------------------------------------
