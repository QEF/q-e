!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
program dynmat
!--------------------------------------------------------------------
!
!  This program
!  - reads a dynamical matrix file produced by the phonon code
!  - adds the nonanalytical part (if Z* and epsilon are read from file),
!    applies the chosen Acoustic Sum Rule (if q=0)
!  - diagonalise the dynamical matrix
!  - calculates IR and Raman cross sections (if Z* and Raman tensors
!    are read from file, respectively)
!  - writes the results to files, both for inspection and for plotting
!
!  Input data (namelist "input")
!
!  fildyn  character input file containing the dynamical matrix
!                    (default: fildyn='matdyn')
!  q(3)      real    calculate LO modes (add nonanalytic terms) along
!                    the direction q (cartesian axis, default: q=(0,0,0) )
!  amass(nt) real    mass for atom type nt, amu
!                    (default: amass is read from file fildyn)
!  asr   character   indicates the type of Acoustic Sum Rule imposed
!                     - 'no': no Acoustic Sum Rules imposed (default)
!                     - 'simple':  previous implementation of the asr used
!                     (3 translational asr imposed by correction of
!                     the diagonal elements of the dynamical matrix)
!                     - 'crystal': 3 translational asr imposed by optimized
!                     correction of the dyn. matrix (projection).
!                     - 'one-dim': 3 translational asr + 1 rotational asr
!                     imposed by optimized correction of the dyn. mat. (the
!                     rotation axis is the direction of periodicity; it
!                     will work only if this axis considered is one of
!                     the cartesian axis).
!                     - 'zero-dim': 3 translational asr + 3 rotational asr
!                     imposed by optimized correction of the dyn. mat.
!                     Note that in certain cases, not all the rotational asr
!                     can be applied (e.g. if there are only 2 atoms in a
!                     molecule or if all the atoms are aligned, etc.).
!                     In these cases the supplementary asr are cancelled
!                     during the orthonormalization procedure (see below).
!                     Finally, in all cases except 'no' a simple correction
!                     on the effective charges is performed (same as in the
!                     previous implementation).
!  axis    integer    indicates the rotation axis for a 1D system
!                     (1=Ox, 2=Oy, 3=Oz ; default =3)
!  lperm   logical    .true. to calculate Gamma-point mode contributions to
!                     dielectric permittivity tensor
!                     (default: lperm=.false.)
!  lplasma logical    .true. to calculate Gamma-point mode effective plasma 
!                     frequencies, automatically triggers lperm = .true. 
!                     (default: lplasma=.false.)
!  filout  character output file containing phonon frequencies and normalized
!                    phonon displacements (i.e. eigenvectors divided by the
!                    square root of the mass and then normalized; they are
!                    not orthogonal)
!                    (default: filout='dynmat.out')
!  fileig  character output file containing phonon frequencies and eigenvectors
!                    of the dynamical matrix (they are orthogonal)
!                    (default: fileig=' ')
!  filmol  character as above, in a format suitable for 'molden'
!                    (default: filmol='dynmat.mold')
!  filxsf  character as above, in axsf format suitable for xcrysden
!                    (default: filxsf='dynmat.axsf')
!
  USE kinds, ONLY: DP
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE environment, ONLY : environment_start, environment_end
  USE io_dyn_mat,  ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail
  USE constants,   ONLY : amu_ry
  USE dynamical
  !
  implicit none
  integer, parameter :: ntypx = 10
  character(len=256):: fildyn, filout, filmol, filxsf, fileig
  character(len=3) :: atm(ntypx)
  character(len=10) :: asr
  logical :: lread, gamma
  complex(DP), allocatable :: z(:,:)
  real(DP) :: amass(ntypx), amass_(ntypx), eps0(3,3), a0, omega, &
       at(3,3), bg(3,3), q(3), q_(3)
  real(DP), allocatable :: w2(:)
  integer :: nat, na, nt, ntyp, iout, axis, nspin_mag, ios
  real(DP) :: celldm(6)
  logical :: xmldyn, lrigid, lraman, lperm, lplasma
  logical, external :: has_xml
  integer :: ibrav, nqs
  integer, allocatable :: itau(:)
  namelist /input/ amass, asr, axis, fildyn, filout, filmol, filxsf, &
                   fileig, lperm, lplasma, q
  !
  ! code is parallel-compatible but not parallel
  !
  CALL mp_startup()
  CALL environment_start('DYNMAT')
  !
  IF (ionode) CALL input_from_file ( )
  !
  asr  = 'no'
  axis = 3
  fildyn='matdyn'
  filout='dynmat.out'
  filmol='dynmat.mold'
  filxsf='dynmat.axsf'
  fileig=' '
  amass(:)=0.0d0
  q(:)=0.0d0
  lperm=.false.
  lplasma=.false.
  !
  IF (ionode) read (5,input, iostat=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('dynmat', 'reading input namelist', ABS(ios))
  !
  CALL mp_bcast(asr,ionode_id, world_comm)
  CALL mp_bcast(axis,ionode_id, world_comm)
  CALL mp_bcast(amass,ionode_id, world_comm)
  CALL mp_bcast(fildyn,ionode_id, world_comm)
  CALL mp_bcast(filout,ionode_id, world_comm)
  CALL mp_bcast(filmol,ionode_id, world_comm)
  CALL mp_bcast(fileig,ionode_id, world_comm)
  CALL mp_bcast(filxsf,ionode_id, world_comm)
  CALL mp_bcast(q,ionode_id, world_comm)
  !
  IF (ionode) inquire(file=fildyn,exist=lread)
  CALL mp_bcast(lread, ionode_id, world_comm)
  IF (lread) THEN
     IF (ionode) WRITE(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file '&
                                     , TRIM(fildyn)
  ELSE
     CALL errore('dynmat', 'File '//TRIM(fildyn)//' not found', 1)
  END IF
  !
  ntyp = ntypx ! avoids spurious out-of-bound errors
  xmldyn=has_xml(fildyn)
  IF (xmldyn) THEN
     CALL read_dyn_mat_param(fildyn,ntyp,nat)
     ALLOCATE (m_loc(3,nat))
     ALLOCATE (tau(3,nat))
     ALLOCATE (ityp(nat))
     ALLOCATE (zstar(3,3,nat))
     ALLOCATE (dchi_dtau(3,3,3,nat) )
     CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass_, tau, ityp, &
             m_loc, nqs, lrigid, eps0, zstar, lraman, dchi_dtau)
     IF (nqs /= 1) CALL errore('dynmat','only q=0 matrix allowed',1)
     a0=celldm(1) ! define alat
     ALLOCATE (dyn(3,3,nat,nat) )
     CALL read_dyn_mat(nat,1,q_,dyn(:,:,:,:))
     CALL read_dyn_mat_tail(nat)
     IF(asr.ne.'no') THEN
         CALL set_asr ( asr, axis, nat, tau, dyn, zstar )
     END IF
     IF (ionode) THEN
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)
        END DO
     END IF
  ELSE
     IF (ionode) THEN
        CALL readmat2 ( fildyn, asr, axis, nat, ntyp, atm, a0, &
                        at, omega, amass_, eps0, q_ )
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)/amu_ry
        END DO
     END IF
  ENDIF
  !
  IF (ionode) THEN
     !
     ! from now on, execute on a single processor
     !
     gamma = ( abs( q_(1)**2+q_(2)**2+q_(3)**2 ) < 1.0d-8 )
     !
     IF (gamma) THEN
        ALLOCATE (itau(nat))
        DO na=1,nat
           itau(na)=na
        END DO
        CALL nonanal ( nat, nat, itau, eps0, q, zstar, omega, dyn )
        DEALLOCATE (itau)
     END IF
     !
     ALLOCATE ( z(3*nat,3*nat), w2(3*nat) )
     CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2,z)
     !
     IF (filout.eq.' ') then
        iout=6
     ELSE
        iout=4
        OPEN (unit=iout,file=filout,status='unknown',form='formatted')
     END IF
     CALL writemodes(nat,q_,w2,z,iout)
     IF(iout .ne. 6) close(unit=iout)
     IF (fileig .ne. ' ') THEN
       OPEN (unit=15,file=TRIM(fileig),status='unknown',form='formatted')
       CALL write_eigenvectors (nat,ntyp,amass,ityp,q_,w2,z,15)
       CLOSE (unit=15)
     ENDIF
     CALL writemolden (filmol, gamma, nat, atm, a0, tau, ityp, w2, z)
     CALL writexsf (filxsf, gamma, nat, atm, a0, at, tau, ityp, z)
     IF (gamma) THEN 
        CALL RamanIR (nat, omega, w2, z, zstar, eps0, dchi_dtau)
        IF (lperm .OR. lplasma) THEN
            CALL polar_mode_permittivity(nat,eps0,z,zstar,w2,omega, &
                                         lplasma)
            IF ( ABS( q(1)**2+q(2)**2+q(3)**2 ) > 1.0d-8 ) &
               WRITE(6,'(5x,a)') 'BEWARE: phonon contribution to &
               & permittivity computed with TO-LO splitting'
        ENDIF
     ENDIF
  ENDIF
  !
  IF (xmldyn) THEN
     DEALLOCATE (m_loc)
     DEALLOCATE (tau)
     DEALLOCATE (ityp)
     DEALLOCATE (zstar)
     DEALLOCATE (dchi_dtau)
     DEALLOCATE (dyn)
  ENDIF
  !
  CALL environment_end('DYNMAT')
  !
  CALL mp_global_end()
  !
end program dynmat
