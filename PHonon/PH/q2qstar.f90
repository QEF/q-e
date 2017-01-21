!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! A small utility that read the first q from a dynamical matrix file (either xml or plain text),
! recomputes the system symmetry (starting from the lattice) and generates the star of q.
!
! Useful for debugging and for producing the star of the wannier-phonon code output.
!
! Syntax:
!   q2qstar.x filein [fileout]
!
! fileout default: filein.rot (old format) or filein.rot.xml (new format)
!
!----------------------------------------------------------------------------
PROGRAM Q2QSTAR
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : amu_ry
  USE parameters,         ONLY : ntypx
  USE mp,                 ONLY : mp_bcast
  USE mp_global,          ONLY : mp_startup, mp_global_end
  USE mp_world,           ONLY : world_comm
  USE io_global,          ONLY : ionode_id, ionode, stdout
  USE environment,        ONLY : environment_start, environment_end
  ! symmetry
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, irt, ftau, copy_sym, nrot, inverse_s
  ! for reading the dyn.mat.
  USE cell_base,          ONLY : at, bg, celldm, ibrav, omega
  USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
  ! as above, unused here
  USE control_ph,         ONLY : xmldyn
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  !
  USE dynmat,             ONLY : w2
  !
  ! for non-xml file only:
  USE dynamicalq,         ONLY : dq_phiq => phiq, dq_tau => tau, dq_ityp => ityp, zeu 
  ! fox xml files only
  USE io_dyn_mat,         ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                                 read_dyn_mat, read_dyn_mat_tail, &
                                 write_dyn_mat_header
  ! small group symmetry
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq
  !
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: CODE="Q2QSTAR"
  CHARACTER(len=256) :: fildyn, filout
  INTEGER :: ierr, nargs
  !
  INTEGER       :: nqs, isq (48), imq, nqq
  REAL(DP)      :: sxq(3, 48), xq(3), xqs(3,48), epsil(3,3)
  !
  LOGICAL :: sym(48), lrigid
  LOGICAL, EXTERNAL :: has_xml
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:,:,:,:), d2(:,:)
  INTEGER :: i,j, icar,jcar, na,nb
  !
  NAMELIST / input / fildyn
  !
  CALL mp_startup()
  CALL environment_start(CODE)
  !
  nargs = command_argument_count()
  IF(nargs < 1) CALL errore(CODE, 'Argument is missing! Syntax: "q2qstar dynfile [outfile]"', 1)
  !
  CALL get_command_argument(1, fildyn)
  CALL mp_bcast(fildyn, ionode_id,world_comm)
  !
  ! check input
  IF (fildyn == ' ')  CALL errore (CODE,' bad fildyn',1)
  xmldyn=has_xml(fildyn)
  !
  ! set up output
  IF (nargs > 1) THEN
    CALL get_command_argument(2, filout)
  ELSE
      filout = TRIM(fildyn)//".rot"
  ENDIF
  CALL mp_bcast(filout, ionode_id,world_comm)
  !
  ! ######################### reading ######################### 
  XML_FORMAT_READ : &
  IF (xmldyn) THEN
    ! read params
    CALL read_dyn_mat_param(fildyn,ntyp,nat)
    ALLOCATE(m_loc(3,nat))
    ALLOCATE(tau(3,nat))
    ALLOCATE(ityp(nat))
    ALLOCATE(zeu(3,3,nat))
    ALLOCATE(phi(3,3,nat,nat))
    ! read system information
    CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                             celldm, at, bg, omega, atm, amass, tau, ityp, &
                             m_loc, nqs, lrigid, epsil, zeu )
    ! read dyn.mat.
    CALL read_dyn_mat(nat,1,xq,phi)
    ! close file
    CALL read_dyn_mat_tail(nat)
    !
  ELSE XML_FORMAT_READ
    ! open file
    IF (ionode)OPEN (unit=1, file=fildyn,status='old',form='formatted',iostat=ierr)
    CALL mp_bcast(ierr, ionode_id,world_comm)
    IF (ierr /= 0) CALL errore(CODE,'file '//TRIM(fildyn)//' missing!',1)
    ! read everything, this use global variables
    ntyp = ntypx
    CALL read_dyn_from_file (nqs, xqs, epsil, lrigid,  &
        ntyp, nat, ibrav, celldm, at, atm, amass)
    !
    IF (ionode) CLOSE(unit=1)
    !
    xq = xqs(:,1)
    ALLOCATE(phi(3,3,nat,nat))
    ALLOCATE(tau(3,nat))
    ALLOCATE(ityp(nat))
    phi  = dq_phiq(:,:,:,:,1)
    tau =  dq_tau
    ityp = dq_ityp
    !zeu =  dq_zeu ! note: zeu from dynamicalq is a real(dp) array, zeu from control_ph is a flag (logical)
    amass = amass/amu_ry
    !
  ENDIF XML_FORMAT_READ
  !
  ! regenerate the lattice
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  at = at / celldm(1)  !  bring at in units of alat
  CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
  CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
!   IF( nqs > 1) CALL errore(CODE, 'This code can in principle read dyn.mat. with the star of q, but it makes no sense', 1)
  WRITE(stdout,'(//,5x,a,3f14.9/)') "Dynamical matrix at q =", xq
  !
  ! ######################### symmetry setup #########################
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
  CALL set_sym_bl ( )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  CALL find_sym ( nat, tau, ityp, .false., m_loc )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !
  ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
  ! part 1: call smallg_q and the copy_sym, 
  minus_q = .true.
  sym = .false.
  sym(1:nsym) = .true.
  CALL smallg_q(xq, 0, at, bg, nsym, s, ftau, sym, minus_q)
  nsymq = copy_sym(nsym, sym)
  ! recompute the inverses as the order of sym.ops. has changed
  CALL inverse_s ( ) 
  ! part 2: this computes gi, gimq
  call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
  IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G:"
  !
  ! finally this does some of the above again and also computes rtau...
  ALLOCATE(rtau( 3, 48, nat))
  CALL sgam_ph_new(at, bg, nsym, s, irt, tau, rtau, nat)
  !
  ! ######################### star of q #########################
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, - 1)
     enddo
  enddo
  CALL symdynph_gq_new (xq, phi, s, invs, rtau, irt, nsymq, nat, &
       irotmq, minus_q)
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, + 1)
     enddo
  enddo
  !
  CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, .true. )
  !
  XML_FORMAT_WRITE : &
  IF (xmldyn) THEN
     nqq=nqs
     IF (imq==0) nqq=2*nqs
!      IF (lgamma.AND.done_epsil.AND.done_zeu) THEN
!         CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
!              celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
!              nqq, epsilon, zstareu, lraman, ramtns)
!      ELSE
        CALL write_dyn_mat_header( filout, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau,ityp,m_loc,nqq)
!      ENDIF
  ELSE XML_FORMAT_WRITE
      OPEN (unit=1, file=filout,status='unknown',form='formatted',iostat=ierr)
      IF (ierr /= 0) CALL errore(CODE,'opening output file',1)
      CALL write_old_dyn_mat_head(1)
  ENDIF XML_FORMAT_WRITE
  !
  ! repack phi to 3*nat,3*nat so that it can be repacked and then rerepacked again in q2qstar_ph
  ALLOCATE(d2(3*nat, 3*nat))
  DO i = 1, 3 * nat
    na = (i - 1) / 3 + 1
    icar = i - 3 * (na - 1)
    DO j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcar = j - 3 * (nb - 1)
        d2 (i, j) = phi(icar, jcar, na, nb)
    ENDDO
  ENDDO
  !
  CALL q2qstar_ph (d2, at, bg, nat, nsym, s, invs, irt, rtau, &
                   nqs, sxq, isq, imq, 1)
  ALLOCATE(w2(3*nat))
  CALL dyndia (xq, 3*nat, nat, ntyp, ityp, amass, 1, d2, w2)

  IF (.not.xmldyn) THEN
    WRITE(1, '(/,3a,/)') "File generated with q2qstar.x from '", TRIM(fildyn), "'" ! <-- to prevent crash with old versions of q2r.x
    CLOSE(1)
  ENDIF
  !
  DEALLOCATE(phi, d2, w2)
  DEALLOCATE(rtau, tau, ityp)
  IF( .not.xmldyn ) DEALLOCATE(dq_phiq, dq_tau, dq_ityp, zeu) ! from read_dyn_from_file
  IF( xmldyn) DEALLOCATE(zeu, m_loc)
  DEALLOCATE(irt) ! from symm_base
  !----------------------------------------------------------------------------
 END PROGRAM Q2QSTAR
!----------------------------------------------------------------------------
!
