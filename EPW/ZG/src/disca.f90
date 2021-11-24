! Copyright (C) 2020-2021 Marios Zacharias, Feliciano Giustino 
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!
Module ifconstants
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  REAL(DP), ALLOCATABLE :: frc(:, :, :, :, :, :, :), tau_blk(:, :),  zeu(:, :, :), &
               m_loc(:, :)
  ! frc : interatomic force constants in REAL space
  ! tau_blk : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  INTEGER, ALLOCATABLE  :: ityp_blk(:)
  ! ityp_blk : atomic types for each atom of the original cell
  !
  CHARACTER(LEN = 3), ALLOCATABLE :: atm(:)
end Module ifconstants
!
!---------------------------------------------------------------------
PROGRAM diff_sca
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !! authors: Marios Zacharias and Feliciano Giustino 
  !! acknowledgement: Hyungjun Lee for help packaging this release
  !! version: v0.2
  !! license: GNU
  !
  !  This program generates diffuse scattering maps for a list of
  !  scattering vectors that are comensurate with the q-grid used to sample each
  !  Brillouin zone. The first part of the code is obtained by 
  !  modifying the matdyn.f90 subroutine of QE. The program starts from the 
  !  interatomic force constants generated from the DFPT phonon code through
  !  the companion program q2r.
  !
  !  disca generates OUTPUTS ... 
  !  Data required for disca are read from the force constant 
  !
  !  Input cards: namelist &input
  !     flfrc     file produced by q2r containing force constants (needed)
  !               It is the same as in the input of q2r.x (+ the .xml extension
  !               if the dynamical matrices produced by ph.x were in xml
  !               format). No default value: must be specified.
  !      asr      (character) indicates the type of Acoustic Sum Rule imposed
  !               - 'no': no Acoustic Sum Rules imposed (default)
  !               - 'simple':  previous implementation of the asr used
  !                  (3 translational asr imposed by correction of
  !                  the diagonal elements of the force constants matrix)
  !               - 'crystal': 3 translational asr imposed by optimized
  !                  correction of the force constants (projection).
  !               - 'one-dim': 3 translational asr + 1 rotational asr
  !                  imposed by optimized correction of the force constants
  !                  (the rotation axis is the direction of periodicity;
  !                   it will work only if this axis considered is one of
  !                   the cartesian axis).
  !               - 'zero-dim': 3 translational asr + 3 rotational asr
  !                  imposed by optimized correction of the force constants
  !               Note that in certain cases, not all the rotational asr
  !               can be applied (e.g. if there are only 2 atoms in a
  !               molecule or if all the atoms are aligned, etc.).
  !               In these cases the supplementary asr are cancelled
  !               during the orthonormalization procedure (see below).
  !     nk1,nk2,nk3  uniform q-point grid for DOS calculation (includes q=0)
  !     at        supercell lattice vectors - must form a superlattice of the
  !               original lattice (default: use original cell)
  !     l1,l2,l3  supercell lattice vectors are original cell vectors times
  !               l1, l2, l3 respectively (default: 1, ignored if at specified)
  !     ntyp      number of atom types in the supercell (default: ntyp of the
  !               original cell)
  !     amass     masses of atoms in the supercell (a.m.u.), one per atom type
  !               (default: use masses read from file flfrc)
  !     readtau   read  atomic positions of the supercell from input
  !               (used to specify different masses) (default: .false.)
  !     q_in_band_form if .true. the q points are given in band form:
  !               Only the first and last point of one or more lines 
  !               are given. See below. (default: .false.).
  !     q_in_cryst_coord if .true. input q points are in crystalline 
  !              coordinates (default: .false.)
  !     eigen_similarity: use similarity of the displacements to order 
  !                       frequencies  (default: .false.)
  !                NB: You cannot use this option with the symmetry
  !                analysis of the modes.
  !     fd         (logical) if .t. the ifc come from the finite displacement calculation
  !     na_ifc     (logical) add non analitic contributions to the interatomic force 
  !                constants if finite displacement method is used (as in Wang et al.
  !                Phys. Rev. B 85, 224303 (2012)) [to be used in conjunction with fd.x]
  !     nosym      if .true., no symmetry and no time reversal are imposed
  !     loto_2d    set to .true. to activate two-dimensional treatment of LO-TO splitting.
  !     loto_disable (logical) if .true. do not apply LO-TO splitting for q=0
  !                  (default: .false.)
  !
  !  if (readtau) atom types and positions in the supercell follow:
  !     (tau(i,na),i= 1,3), ityp(na)
  !  IF (q_in_band_form.and..not.dos) THEN
  !     nq     ! number of q points
  !     (q(i,n),i= 1,3), nptq   nptq is the number of points between this point
  !                            and the next. These points are automatiCALLy
  !                            generated. the q points are given in Cartesian
  !                            coordinates, 2pi/a units (a=lattice parameters)
  !  ELSE, if (.not.dos) :
  !     nq         number of q-points
  !     (q(i,n), i= 1,3)    nq q-points in cartesian coordinates, 2pi/a units
  !  If q = 0, the direction qhat (q=>0) for the non-analytic part
  !  is extracted from the sequence of q-points as follows:
  !     qhat = q(n) - q(n-1)  or   qhat = q(n) - q(n+1)
  !  depending on which one is available and nonzero.
  !  For low-symmetry crystals, specify twice q = 0 in the list
  !  if you want to have q = 0 results for two different directions
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, mp_global_end, inter_pool_comm
  USE environment, ONLY : environment_start, environment_end
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_ifc_param, read_ifc
  USE cell_base,  ONLY : at, bg, celldm
  USE constants,  ONLY : RY_TO_THZ, RY_TO_CMM1, amu_ry
  USE symm_base,  ONLY : set_sym
  USE rap_point_group,  ONLY : code_group
  USE bz_form,    ONLY : transform_label_coord
  USE rigid,      ONLY : dyndiag, nonanal, nonanal_ifc
  USE parser,     ONLY : read_line

  USE ifconstants, ONLY : frc, atm, zeu, tau_blk, ityp_blk, m_loc
  !
  IMPLICIT NONE
  !
  INTEGER :: gid
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  INTEGER:: nax, nax_blk
  INTEGER, PARAMETER:: ntypx = 10, nrwsx = 200
  REAL(DP), PARAMETER :: eps = 1.0d-6
  INTEGER :: nr1, nr2, nr3, nsc, nk1, nk2, nk3, ibrav, qstart, qfinal
  CHARACTER(LEN = 256) :: flfrc, filename
  CHARACTER(LEN = 10)  :: asr
  LOGICAL :: has_zstar, q_in_cryst_coord, eigen_similarity, loto_disable
  COMPLEX(DP), ALLOCATABLE :: dyn(:, :, :, :), dyn_blk(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE :: z(:,:), frc_ifc(:, :, :, :)
  REAL(DP), ALLOCATABLE:: tau(:, :), q(:, :), w2(:, :), freq(:,:), wq(:)
  REAL(DP), ALLOCATABLE:: q_super(:, :), q_strft(:, :), q_strft_part(:, :)
  INTEGER, ALLOCATABLE:: ityp(:), itau_blk(:)
  REAL(DP) ::     omega,alat, &! cell parameters and volume
                  at_blk(3, 3), bg_blk(3, 3),  &! original cell
                  omega_blk,                 &! original cell volume
                  epsil(3, 3),                &! dielectric tensor
                  amass(ntypx),              &! atomic masses
                  amass_blk(ntypx),          &! original atomic masses
                  atws(3, 3),      &! lattice vector for WS initialization
                  rws(0 : 3, nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  !
  INTEGER :: nat, nat_blk, ntyp, ntyp_blk, &
             l1, l2, l3,                   &! supercell dimensions
             nrws,                         &! number of nearest neighbor
             code_group_old

  INTEGER :: nspin_mag, nqs, ios
  !
  LOGICAL :: readtau, xmlifc, lo_to_split, na_ifc, fd, nosym, loto_2d
  !
  REAL(DP) :: qhat(3), qh, E, qq
  REAL(DP) :: delta 
  REAL(DP), ALLOCATABLE :: xqaux(:, :)
  INTEGER, ALLOCATABLE :: nqb(:)
  INTEGER :: n, i, j, it, nq, nq_super, nq_strft, nqx, na, nb, nqtot
  LOGICAL, EXTERNAL :: has_xml
  INTEGER, ALLOCATABLE :: num_rap_mode(:, :)
  LOGICAL, ALLOCATABLE :: high_sym(:)
  LOGICAL :: q_in_band_form
  ! .... variables for band plotting based on similarity of eigenvalues
  COMPLEX(DP), ALLOCATABLE :: tmp_z(:, :)
  REAL(DP), ALLOCATABLE :: ABS_similarity(:, :), tmp_w2(:)
  COMPLEX(DP), ALLOCATABLE :: f_of_q(:, :, :, :)
  INTEGER :: location(1), isig
  CHARACTER(LEN=6) :: int_to_char
  LOGICAL, ALLOCATABLE :: mask(:)
  INTEGER            :: npk_label, nch
  CHARACTER(LEN=3), ALLOCATABLE :: letter(:)
  INTEGER, ALLOCATABLE :: label_list(:)
  LOGICAL :: tend, terr
  CHARACTER(LEN=256) :: input_line, buffer
  CHARACTER(LEN= 10) :: point_label_type
  CHARACTER(len=80) :: k_points = 'tpiba'
  ! mz_b
  COMPLEX(DP), ALLOCATABLE :: z_nq_zg(:, :, :) ! nomdes,nmodes,nq
  REAL(DP),    ALLOCATABLE :: q_nq_zg(:, :), qlist_strf(:, :), qlist_strf_cart(:, :) ! 3,nq
  REAL(DP)                 :: atmsf_a(ntypx, 5), atmsf_b(ntypx, 5)
  LOGICAL                  :: disca, q_external, print_raw_data
  LOGICAL                  :: zero_one_phonon, atom_resolved, mode_resolved, full_phonon
  INTEGER                  :: dim1, dim2, dim3, nks1, nks2, nks3, nksf1, nksf2, nksf3
  INTEGER                  :: lower_bnd, upper_bnd
  INTEGER                  :: T, plane_dir
  REAL(DP)                 :: plane_val, eps2
  CHARACTER(LEN=3)         :: atm_zg(ntypx)
  ! inputs for pp_disca, rotate
  INTEGER                  :: nrots, kres1, kres2, col1, col2, Np
  REAL(DP)                 :: kmin, kmax
  ! mz_e
  !
  NAMELIST /input/ flfrc, amass, asr, at, &
       &           nk1, nk2, nk3, l1, l2, l3, ntyp, readtau, &
       &           q_in_band_form, q_in_cryst_coord, &
       &           eigen_similarity, na_ifc, fd, point_label_type, nosym, &
       &           loto_2d, loto_disable, &
  ! mz_b we add the inputs for diffuse scattering
       &           disca, dim1, dim2, dim3, atm_zg, mode_resolved, & 
       &           T, q_external, zero_one_phonon, print_raw_data, & 
       &           nks1, nks2, nks3, plane_val, plane_dir, qstart, qfinal, atom_resolved, & 
       &           full_phonon, atmsf_a, atmsf_b, nksf1, nksf2, nksf3, eps2
  ! disca --> if true compute the diffuse_scattering 
  !we add a new list for the inputs og rotate and pp_disca
  NAMELIST /pp_disca/ nrots, kres1, kres2, kmin, kmax, col1, col2, Np
  ! mz_e 
  !
  CALL mp_startup()
  CALL environment_start('DISCA')
  !
  IF (ionode) CALL input_from_file ( )
     !
     ! ... all calculations are done by the first cpu
     !
     ! set namelist default
     !
     nk1 = 0
     nk2 = 0
     nk3 = 0
     asr = 'no'
     readtau=.FALSE.
     flfrc=' '
     amass(:) =0.d0
     amass_blk(:) =0.d0
     at(:,:) = 0.d0
     ntyp = 0
     l1 = 1
     l2 = 1
     l3 = 1
     q_in_band_form = .FALSE.
     eigen_similarity = .FALSE.
     q_in_cryst_coord = .FALSE.
     na_ifc = .FALSE.
     fd = .FALSE.
     point_label_type = 'SC'
     nosym = .FALSE.
     loto_2d = .FALSE.
     loto_disable = .FALSE.
     ! mz_b
     disca = .TRUE.
     T = 0
     q_external = .FALSE.
     zero_one_phonon = .TRUE.
     full_phonon = .FALSE.
     atom_resolved = .FALSE.
     mode_resolved = .FALSE.
     print_raw_data = .FALSE.
     dim1 = 0
     dim2 = 0
     dim3 = 0
     atm_zg = "Element"
     atmsf_a = 0.d0
     atmsf_b = 0.d0
     plane_val = 0.d0
     plane_dir = 3
     nks1 = 0
     nks2 = 0
     nks3 = 0
     nksf1 = 6
     nksf2 = 6
     nksf3 = 6
     qstart = 1
     qfinal = 2
     eps2 = 1.0d-15
     ! 
     nrots = 1
     kres1 = 250
     kres2 = 250
     kmin = -5
     kmax = 10
     col1 = 1
     col2 = 2
     Np = 100
     ! mz_e 
     !
     !
     IF (ionode) READ (5,input,IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, world_comm) 
     CALL errore('disca', 'reading input namelist', ABS(ios))
     IF (ionode) READ (5,pp_disca,IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, world_comm)
     CALL errore('pp_disca', 'reading pp_disca namelist', ABS(ios))
     CALL mp_bcast(nk1,ionode_id, world_comm)
     CALL mp_bcast(nk2,ionode_id, world_comm)
     CALL mp_bcast(nk3,ionode_id, world_comm)
     CALL mp_bcast(asr,ionode_id, world_comm)
     CALL mp_bcast(readtau,ionode_id, world_comm)
     CALL mp_bcast(flfrc,ionode_id, world_comm)
     CALL mp_bcast(amass,ionode_id, world_comm)
     CALL mp_bcast(amass_blk,ionode_id, world_comm)
     CALL mp_bcast(at,ionode_id, world_comm)
     CALL mp_bcast(ntyp,ionode_id, world_comm)
     CALL mp_bcast(l1,ionode_id, world_comm)
     CALL mp_bcast(l2,ionode_id, world_comm)
     CALL mp_bcast(l3,ionode_id, world_comm)
     CALL mp_bcast(na_ifc,ionode_id, world_comm)
     CALL mp_bcast(fd,ionode_id, world_comm)
     CALL mp_bcast(q_in_band_form,ionode_id, world_comm)
     CALL mp_bcast(eigen_similarity,ionode_id, world_comm)
     CALL mp_bcast(q_in_cryst_coord,ionode_id, world_comm)
     CALL mp_bcast(point_label_type,ionode_id, world_comm)
     CALL mp_bcast(loto_2d,ionode_id, world_comm)
     CALL mp_bcast(loto_disable,ionode_id, world_comm)
     ! mz_b
     CALL mp_bcast(disca,ionode_id, world_comm)
     CALL mp_bcast(T,ionode_id, world_comm)
     CALL mp_bcast(q_external, ionode_id, world_comm)
     CALL mp_bcast(zero_one_phonon, ionode_id, world_comm)
     CALL mp_bcast(atom_resolved, ionode_id, world_comm)
     CALL mp_bcast(mode_resolved, ionode_id, world_comm)
     CALL mp_bcast(print_raw_data, ionode_id, world_comm)
     CALL mp_bcast(full_phonon, ionode_id, world_comm)
     CALL mp_bcast(atmsf_a, ionode_id, world_comm)
     CALL mp_bcast(atmsf_b, ionode_id, world_comm)
     CALL mp_bcast(dim1, ionode_id, world_comm)
     CALL mp_bcast(dim2, ionode_id, world_comm)
     CALL mp_bcast(dim3, ionode_id, world_comm)
     CALL mp_bcast(nks1, ionode_id, world_comm)
     CALL mp_bcast(nks2, ionode_id, world_comm)
     CALL mp_bcast(nks3, ionode_id, world_comm)
     CALL mp_bcast(nksf1, ionode_id, world_comm)
     CALL mp_bcast(nksf2, ionode_id, world_comm)
     CALL mp_bcast(nksf3, ionode_id, world_comm)
     CALL mp_bcast(plane_val, ionode_id, world_comm)
     CALL mp_bcast(plane_dir, ionode_id, world_comm)
     CALL mp_bcast(qstart, ionode_id, world_comm)
     CALL mp_bcast(qfinal, ionode_id, world_comm)
     CALL mp_bcast(eps2, ionode_id, world_comm)
     ! 
     CALL mp_bcast(nrots, ionode_id, world_comm)
     CALL mp_bcast(kres1, ionode_id, world_comm)
     CALL mp_bcast(kres2, ionode_id, world_comm)
     CALL mp_bcast(kmin, ionode_id, world_comm)
     CALL mp_bcast(kmax, ionode_id, world_comm)
     CALL mp_bcast(col1, ionode_id, world_comm)
     CALL mp_bcast(col2, ionode_id, world_comm)
     CALL mp_bcast(Np, ionode_id, world_comm)
     ! 
     IF (loto_2d .AND. loto_disable) CALL errore('disca', &
         'loto_2d and loto_disable cannot be both true', 1)
     ! To check that use specify supercell dimensions
     IF (disca) THEN 
        IF ((dim1 .LT. 1)  .OR. (dim2 .LT. 1) .OR. (dim3 .LT. 1)) CALL errore('disca', 'reading supercell size', dim1)
        IF ((qfinal .LT. qstart)) CALL errore('disca', 'qfinal smaller than qstart', qfinal)
     ENDIF
     !
     ! read force constants
     !
     ntyp_blk = ntypx ! avoids fake out-of-bound error
     xmlifc=has_xml(flfrc)
     IF (xmlifc) THEN
        CALL read_dyn_mat_param(flfrc,ntyp_blk,nat_blk)
        ALLOCATE (m_loc(3,nat_blk))
        ALLOCATE (tau_blk(3,nat_blk))
        ALLOCATE (ityp_blk(nat_blk))
        ALLOCATE (atm(ntyp_blk))
        ALLOCATE (zeu(3,3,nat_blk))
        CALL read_dyn_mat_header(ntyp_blk, nat_blk, ibrav, nspin_mag, &
                 celldm, at_blk, bg_blk, omega_blk, atm, amass_blk, &
                 tau_blk, ityp_blk,  m_loc, nqs, has_zstar, epsil, zeu )
        alat=celldm(1)
        CALL volume(alat,at_blk(1,1),at_blk(1,2),at_blk(1,3),omega_blk)
        CALL read_ifc_param(nr1,nr2,nr3)
        ALLOCATE(frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk))
        CALL read_ifc(nr1,nr2,nr3,nat_blk,frc)
     ELSE
        CALL readfc ( flfrc, nr1, nr2, nr3, epsil, nat_blk, &
            ibrav, alat, at_blk, ntyp_blk, &
            amass_blk, omega_blk, has_zstar)
     ENDIF
     !
     CALL recips ( at_blk(1,1),at_blk(1,2),at_blk(1,3),  &
          bg_blk(1,1),bg_blk(1,2),bg_blk(1,3) )
     !
     ! set up (super)cell
     !
     if (ntyp < 0) then
        CALL errore ('disca','wrong ntyp ', ABS(ntyp))
     else if (ntyp == 0) then
        ntyp=ntyp_blk
     end if
     !
     ! masses (for mass approximation)
     !
     DO it= 1,ntyp
        IF (amass(it) < 0.d0) THEN
           CALL errore ('disca','wrong mass in the namelist',it)
        ELSE IF (amass(it) == 0.d0) THEN
           IF (it.LE.ntyp_blk) THEN
              WRITE (stdout,'(a,i3,a,a)') ' mass for atomic type ',it,      &
                   &                     ' not given; uses mass from file ',flfrc
              amass(it) = amass_blk(it)
           ELSE
              CALL errore ('disca','missing mass in the namelist',it)
           ENDIF
        ENDIF
     ENDDO
     !
     ! lattice vectors
     !
     IF (SUM(ABS(at(:,:))) == 0.d0) THEN
        IF (l1.LE.0 .OR. l2.LE.0 .OR. l3.LE.0) CALL                    &
             &             errore ('disca',' wrong l1,l2 or l3',1)
        at(:,1) = at_blk(:,1)*DBLE(l1)
        at(:,2) = at_blk(:,2)*DBLE(l2)
        at(:,3) = at_blk(:,3)*DBLE(l3)
     ENDIF
     !
     CALL check_at(at,bg_blk,alat,omega)
     !
     ! the supercell contains "nsc" times the original unit cell
     !
     nsc = NINT(omega/omega_blk)
     IF (ABS(omega/omega_blk-nsc) > eps) &
          CALL errore ('disca', 'volume ratio not INTEGER', 1)
     !
     ! read/generate atomic positions of the (super)cell
     !
     nat = nat_blk * nsc
     nax = nat
     nax_blk = nat_blk
     !
     ALLOCATE ( tau (3, nat), ityp(nat), itau_blk(nat) )
     !
     IF (readtau) THEN
        CALL read_tau &
             (nat, nat_blk, ntyp, bg_blk, tau, tau_blk, ityp, itau_blk)
     ELSE
        CALL set_tau  &
             (nat, nat_blk, at, at_blk, tau, tau_blk, ityp, ityp_blk, itau_blk)
     ENDIF
     !
     !
     ! reciprocal lattice vectors
     !
     CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
     !
     ! build the WS cell corresponding to the force constant grid
     !
     atws(:,1) = at_blk(:,1)*DBLE(nr1)
     atws(:,2) = at_blk(:,2)*DBLE(nr2)
     atws(:,3) = at_blk(:,3)*DBLE(nr3)
     ! initialize WS r-vectors
     CALL wsinit(rws,nrwsx,nrws,atws)
     !
     ! end of (super)cell setup
     !
     IF (.NOT. q_external) THEN 
       CALL qpoint_gen1(dim1, dim2, dim3, nq_super) ! nq_super are the q-points
                                                    ! in sets AB commensurate to the
                                                    ! supercell size
       ! nq_super = ctrAB
       CALL mp_bcast(nq_super, ionode_id, world_comm)
       !
       IF (ionode) WRITE (*,*) "Brillouin zone q-points", nq_super
       !
       ALLOCATE ( q_super(3, nq_super) )        
       CALL mp_bcast(q_super, ionode_id, world_comm)
       !
       CALL qpoint_gen2(dim1, dim2, dim3, nq_super, q_super) 
       !
       CALL mp_bcast(q_super, ionode_id, world_comm)
       !
       CALL cryst_to_cart(nq_super, q_super, bg, +1) ! convert them to Cartesian
       ! 
       ! we first read nq for structure factor q-points - S
       ALLOCATE(qlist_strf(dim1 * dim2 * dim3 * (nksf1 - nks1) * (nksf2 - nks2) * (nksf3 - nks3), 3))
       ALLOCATE(qlist_strf_cart(dim1 * dim2 * dim3 * (nksf1 - nks1) * (nksf2 - nks2) * (nksf3 - nks3), 3)) 
       IF (ionode) CALL qpoint_gen_map_1(dim1, dim2, dim3, nks1, nks2, nks3, nksf1, nksf2, nksf3, & 
                         plane_val, plane_dir, nq_strft, qlist_strf, qlist_strf_cart) 
       ! nq_strft is the # of q-points
       CALL mp_bcast(nq_strft, ionode_id, world_comm)
       !
       IF (ionode) WRITE (*,*) 
       IF (ionode) WRITE (*,*) "Generating", nq_strft, "q-points for the full map ..." 
       !
       IF (qfinal .GT. nq_strft) CALL errore('disca', 'qfinal larger than full # of q-points', qfinal)
       !
       ALLOCATE (q_strft(3, nq_strft), q_strft_part(3, qfinal - qstart + 1 ))
       !
       IF (ionode) CALL qpoint_gen_map_2(dim1, dim2, dim3, nks1, nks2, nks3, nksf1, nksf2, nksf3, & 
                                         plane_val, plane_dir, nq_strft, qlist_strf, qlist_strf_cart, q_strft) ! q_strft array with q-points
       DEALLOCATE(qlist_strf, qlist_strf_cart)
       !
       CALL mp_bcast(q_strft, ionode_id, world_comm)
       IF (ionode) then
         OPEN (unit = 30, file = 'qpts_strf.dat', status = 'unknown', form = 'formatted')
         DO n = 1, nq_strft
            WRITE(30, '(3F10.6)') q_strft(:, n)
         ENDDO
         CLOSE(30)
       ENDIF
       q_strft_part(:,:) = q_strft(:, qstart : qfinal)
       nq_strft = qfinal - qstart + 1
       CALL mp_bcast(nq_strft, ionode_id, world_comm)
       !
       IF (ionode) WRITE (*,*) "The calculation is for", qfinal - qstart + 1, "q-points ..." 
       !
       ! convert them to Cartesian
       CALL cryst_to_cart(nq_strft, q_strft_part, bg, +1)
       !
       ! Now ALLOCATE and read full q-list containing 
       ! supercell q-points and structure factor q-points 
       nq = nq_super + nq_strft
       ALLOCATE (q(3, nq))
       DO n = 1, nq_super 
          q(:, n) = q_super(:, n)
       ENDDO
       !
       DO n = nq_super + 1, nq
          q(:, n) = q_strft_part(:, n - nq_super) 
       ENDDO
       ! 
       DEALLOCATE(q_strft_part, q_super)
     ELSE
     ! mz_ends
        IF (ionode) READ (5,*) nq
        CALL mp_bcast(nq, ionode_id, world_comm)
        ALLOCATE ( q(3,nq) )
        IF (.NOT.q_in_band_form) THEN
           DO n = 1,nq
     ! mz_edits
              IF (ionode) READ (5,*) (q(i,n),i= 1,3)
     ! IF (ionode) READ (5,'(3F10.6)') q(:,n) 
     ! mz_done
           ENDDO
           CALL mp_bcast(q, ionode_id, world_comm)
           !
           IF (q_in_cryst_coord)  CALL cryst_to_cart(nq,q,bg,+1)
        ELSE
           ALLOCATE( nqb(nq) )
           ALLOCATE( xqaux(3,nq) )
           ALLOCATE( letter(nq) )
           ALLOCATE( label_list(nq) )
           npk_label=0
           DO n = 1, nq
              CALL read_line( input_line, end_of_file = tend, error = terr )
              IF (tend) CALL errore('disca','Missing lines',1)
              IF (terr) CALL errore('disca','Error reading q points',1)
              DO j= 1,256   ! loop over all characters of input_line
                 IF ( (ICHAR(input_line(j:j)) < 58 .AND. &   ! a digit
                       ICHAR(input_line(j:j)) > 47)      &
                   .OR.ICHAR(input_line(j:j)) == 43 .OR. &   ! the + sign
                       ICHAR(input_line(j:j)) == 45 .OR. &   ! the - sign
                       ICHAR(input_line(j:j)) == 46 ) THEN   ! a dot .
!
!   This is a digit, therefore this line contains the coordinates of the
!   k point. We read it and exit from the loop on characters
!
                     READ(input_line,*) xqaux(1,n), xqaux(2,n), xqaux(3,n), &
                                                    nqb(n)
                     EXIT
                 ELSEIF ((ICHAR(input_line(j:j)) < 123 .AND. &
                          ICHAR(input_line(j:j)) > 64))  THEN
!
!   This is a letter, not a space character. We read the next three 
!   characters and save them in the letter array, save also which k point
!   it is
!
                    npk_label=npk_label+1
                    READ(input_line(j:),'(a3)') letter(npk_label)
                    label_list(npk_label)=n
!
!  now we remove the letters from input_line and read the number of points
!  of the line. The next two line should account for the case in which
!  there is only one space between the letter and the number of points.
!
                    nch=3
                    IF ( ICHAR(input_line(j+1:j+1))==32 .OR. &
                         ICHAR(input_line(j+2:j+2))==32 ) nch=2
                    buffer=input_line(j+nch:)
                    READ(buffer,*,err=20,iostat=ios) nqb(n)
20                  IF (ios /=0) CALL errore('disca',&
                                      'problem reading number of points',1)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
           IF (q_in_cryst_coord) k_points='crystal'
           IF ( npk_label > 0 ) &
              CALL transform_label_coord(ibrav, celldm, xqaux, letter, &
                   label_list, npk_label, nq, k_points, point_label_type )

           DEALLOCATE(letter)
           DEALLOCATE(label_list)

           CALL mp_bcast(xqaux, ionode_id, world_comm)
           CALL mp_bcast(nqb, ionode_id, world_comm)
           IF (q_in_cryst_coord)  CALL cryst_to_cart(nq,xqaux,bg,+1)
           nqtot=SUM(nqb(1:nq-1))+1
           DO i= 1,nq-1
              IF (nqb(i)==0) nqtot=nqtot+1
           ENDDO
           DEALLOCATE(q)
           ALLOCATE(q(3,nqtot))
           ALLOCATE(wq(nqtot))
           CALL generate_k_along_lines(nq, xqaux, nqb, q, wq, nqtot)
           nq=nqtot
           DEALLOCATE(xqaux)
           DEALLOCATE(nqb)
        ENDIF
        ! 
     ENDIF ! q_external
     !
     IF (asr /= 'no') THEN
        CALL set_asr (asr, nr1, nr2, nr3, frc, zeu, &
             nat_blk, ibrav, tau_blk)
     ENDIF
     !
     ALLOCATE ( dyn(3,3,nat,nat), dyn_blk(3,3,nat_blk,nat_blk) )
     ALLOCATE ( z(3*nat,3*nat), w2(3*nat,nq), f_of_q(3,3,nat,nat) )
     ALLOCATE ( tmp_w2(3*nat), ABS_similarity(3*nat,3*nat), mask(3*nat) )
     ! mz_b
     w2 = 0.d0
     IF (disca) THEN
      ALLOCATE ( z_nq_zg(3 * nat, 3 * nat, nq),q_nq_zg(3, nq))
      z_nq_zg(:,:,:)= (0.d0, 0.d0)
      q_nq_zg(:,:) = 0.d0
     ENDIF
     ! mz_e

     IF (xmlifc) CALL set_sym(nat, tau, ityp, nspin_mag, m_loc )

     ALLOCATE(num_rap_mode(3*nat,nq))
     ALLOCATE(high_sym(nq))
     num_rap_mode=-1
     high_sym=.TRUE.

     CALL fkbounds( nq, lower_bnd, upper_bnd )
     !
     DO n= lower_bnd, upper_bnd ! 1, nq
        dyn(:,:,:,:) = (0.d0, 0.d0)

        lo_to_split=.FALSE.
        f_of_q(:,:,:,:) = (0.d0,0.d0)

        IF(na_ifc) THEN

           qq=SQRT(q(1,n)**2+q(2,n)**2+q(3,n)**2)
           if(ABS(qq) < 1d-8) qq= 1.0
           qhat(1)=q(1,n)/qq
           qhat(2)=q(2,n)/qq
           qhat(3)=q(3,n)/qq

           CALL nonanal_ifc (nat,nat_blk,itau_blk,epsil,qhat,zeu,omega,dyn, &
                           nr1, nr2, nr3,f_of_q)
        ENDIF

        CALL setupmat (q(1,n), dyn, nat, at, bg, tau, itau_blk, nsc, alat, &
             dyn_blk, nat_blk, at_blk, bg_blk, tau_blk, omega_blk,  &
             loto_2d, &
             epsil, zeu, frc, nr1,nr2,nr3, has_zstar, rws, nrws, na_ifc,f_of_q,fd)

        IF (.not.loto_2d) THEN 
          qhat(1) = q(1,n)*at(1,1)+q(2,n)*at(2,1)+q(3,n)*at(3,1)
          qhat(2) = q(1,n)*at(1,2)+q(2,n)*at(2,2)+q(3,n)*at(3,2)
          qhat(3) = q(1,n)*at(1,3)+q(2,n)*at(2,3)+q(3,n)*at(3,3)
          IF ( ABS( qhat(1) - NINT (qhat(1) ) ) <= eps .AND. &
               ABS( qhat(2) - NINT (qhat(2) ) ) <= eps .AND. &
               ABS( qhat(3) - NINT (qhat(3) ) ) <= eps ) THEN
             !
             ! q = 0 : we need the direction q => 0 for the non-analytic part
             !
             IF ( n == 1 ) THEN
                ! if q is the first point in the list
                IF ( nq > 1 ) THEN
                   ! one more point
                   qhat(:) = q(:,n) - q(:,n+1)
                ELSE
                   ! no more points
                   qhat(:) = 0.d0
                ENDIF
             ELSE IF ( n > 1 ) THEN
                ! if q is not the first point in the list
                IF ( q(1,n-1)==0.d0 .AND. &
                     q(2,n-1)==0.d0 .AND. &
                     q(3,n-1)==0.d0 .AND. n < nq ) THEN
                   ! if the preceding q is also 0 :
                   qhat(:) = q(:,n) - q(:,n+1)
                ELSE
                   ! if the preceding q is npt 0 :
                   qhat(:) = q(:,n) - q(:,n-1)
                ENDIF
             ENDIF
             qh = SQRT(qhat(1)**2+qhat(2)**2+qhat(3)**2)
             ! WRITE(*,*) ' qh,  has_zstar ',qh,  has_zstar
             IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
             IF (qh /= 0.d0 .AND. .NOT. has_zstar) THEN
                  IF (ionode) WRITE(*,*)
                  CALL infomsg  &
                  ('disca','Z* not found in file '//TRIM(flfrc)// &
                            ', TO-LO splitting at q=0 will be absent!')
             ELSEIF (loto_disable) THEN
                CALL infomsg('disca', &
                    'loto_disable is true. Disable LO-TO splitting at q=0.')
             ELSE
                lo_to_split=.TRUE.
             ENDIF
             !
             IF (lo_to_split) CALL nonanal (nat, nat_blk, itau_blk, epsil, qhat, zeu, omega, dyn)
             !
          ENDIF
        !
        ENDIF
        !
        CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2(1,n),z)
        !
        ! Fill a 3D matrix with all eigenvectors
        !
        IF (disca) THEN
           z_nq_zg(:,:,n) = z(:,:)               
           q_nq_zg(:,n) = q(:,n)
        ENDIF
        ! Cannot use the small group of \Gamma to analize the symmetry
        ! of the mode if there is an electric field.
        !
        IF (xmlifc.AND..NOT.lo_to_split) THEN
             WRITE(stdout,'(10x,"xq=",3F8.4)') q(:,n)
             CALL find_representations_mode_q(nat,ntyp,q(:,n), &
                       w2(:,n),z,tau,ityp,amass, num_rap_mode(:,n), nspin_mag)
            IF (code_group==code_group_old.OR.high_sym(n-1)) high_sym(n)=.FALSE.
            code_group_old=code_group
        ENDIF

        IF (eigen_similarity) THEN
           ! ... order phonon dispersions using similarity of eigenvalues
           ! ... Courtesy of Takeshi Nishimatsu, IMR, Tohoku University 
           IF (.NOT.ALLOCATED(tmp_z)) THEN
              ALLOCATE(tmp_z(3*nat,3*nat))
           ELSE
              ABS_similarity = ABS ( MATMUL ( CONJG( TRANSPOSE(z)), tmp_z ) )
              mask(:) = .true.
              DO na= 1,3*nat
                 location = maxloc( ABS_similarity(:,na), mask(:) )
                 mask(location(1)) = .false.
                 tmp_w2(na) = w2(location(1),n)
                 tmp_z(:,na) = z(:,location(1))
              ENDDO
              w2(:,n) = tmp_w2(:)
              z(:,:) = tmp_z(:,:)
           ENDIF
           tmp_z(:,:) = z(:,:)
        ENDIF
        !
        !
        !
     ENDDO  !nq
     ! 
     CALL mp_sum(z_nq_zg, inter_pool_comm)
     CALL mp_sum(q_nq_zg, inter_pool_comm)
     CALL mp_sum(w2, inter_pool_comm)
     CALL mp_barrier(inter_pool_comm)
     DEALLOCATE (tmp_w2, ABS_similarity, mask)
     IF (eigen_similarity) DEALLOCATE(tmp_z)
     !
     !mz_b
     IF (disca) CALL diffuse_scattering(nq, nq_super, nq_strft, nat, ntyp, amass, ityp, q_nq_zg, & 
                         w2, z_nq_zg, q_external, zero_one_phonon, full_phonon, atmsf_a, atmsf_b, & 
                         dim1, dim2, dim3, tau, alat, atm_zg, atom_resolved, mode_resolved, & 
                         ntypx, at, q_in_cryst_coord, T, eps2, print_raw_data, & 
                         nrots, kres1, kres2, kmin, kmax, col1, col2, Np)
     !mz_e
     ! 
     !
     DEALLOCATE (z, w2, dyn, dyn_blk)
     ! mz_b
     IF (disca) DEALLOCATE (z_nq_zg, q_nq_zg) 
     ! mz_e
     !
     !
     !DEALLOCATE ( freq)
     DEALLOCATE(num_rap_mode)
     DEALLOCATE(high_sym)
  !

  CALL environment_end('DISCA')
  !
  CALL mp_global_end()
  !
  STOP
  !
END PROGRAM diff_sca
!
!-----------------------------------------------------------------------
SUBROUTINE readfc ( flfrc, nr1, nr2, nr3, epsil, nat,    &
                    ibrav, alat, at, ntyp, amass, omega, has_zstar )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ifconstants,ONLY : tau => tau_blk, ityp => ityp_blk, frc, zeu
  USE cell_base,  ONLY : celldm
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp,         ONLY : mp_bcast 
  USE mp_world,   ONLY : world_comm 
  USE constants,  ONLY : amu_ry
  !
  IMPLICIT NONE
  ! I/O variable
  CHARACTER(LEN=256) :: flfrc
  INTEGER :: ibrav, nr1,nr2,nr3,nat, ntyp
  REAL(DP) :: alat, at(3,3), epsil(3,3)
  LOGICAL :: has_zstar
  ! local variables
  INTEGER :: i, j, na, nb, m1,m2,m3
  INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  REAL(DP) :: amass(ntyp), amass_from_file, omega
  INTEGER :: nt
  CHARACTER(LEN=3) :: atm
  !
  !
  IF (ionode) OPEN (unit= 1,file=flfrc,status='old',form='formatted')
  !
  !  read cell data
  !
  IF (ionode)THEN
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i= 1,6)
     if (ibrav==0) then
        read(1,*) ((at(i,j),i= 1,3),j=1,3)
     end if
  ENDIF
  CALL mp_bcast(ntyp, ionode_id, world_comm)
  CALL mp_bcast(nat, ionode_id, world_comm)
  CALL mp_bcast(ibrav, ionode_id, world_comm)
  CALL mp_bcast(celldm, ionode_id, world_comm)
  IF (ibrav==0) THEN
     CALL mp_bcast(at, ionode_id, world_comm)
  ENDIF
  !
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  !  read atomic types, positions and masses
  !
  DO nt = 1,ntyp
     IF (ionode) READ(1,*) i,atm,amass_from_file
     CALL mp_bcast(i,ionode_id, world_comm)
     CALL mp_bcast(atm,ionode_id, world_comm)
     CALL mp_bcast(amass_from_file,ionode_id, world_comm)
     IF (i.NE.nt) CALL errore ('readfc','wrong data read',nt)
     IF (amass(nt).EQ.0.d0) THEN
        amass(nt) = amass_from_file/amu_ry
     ELSE
        WRITE(stdout,*) 'for atomic type',nt,' mass from file not used'
     ENDIF
  ENDDO
  !
  ALLOCATE (tau(3,nat), ityp(nat), zeu(3,3,nat))
  !
  DO na= 1,nat
     IF (ionode) READ(1,*) i,ityp(na),(tau(j,na),j= 1,3)
     CALL mp_bcast(i,ionode_id, world_comm)
     IF (i.NE.na) CALL errore ('readfc','wrong data read',na)
  ENDDO
  CALL mp_bcast(ityp,ionode_id, world_comm)
  CALL mp_bcast(tau,ionode_id, world_comm)
  !
  !  read macroscopic variable
  !
  IF (ionode) READ (1,*) has_zstar
  CALL mp_bcast(has_zstar,ionode_id, world_comm)
  IF (has_zstar) THEN
     IF (ionode) READ(1,*) ((epsil(i,j),j= 1,3),i=1,3)
     CALL mp_bcast(epsil,ionode_id, world_comm)
     IF (ionode) THEN
        DO na= 1,nat
           READ(1,*)
           READ(1,*) ((zeu(i,j,na),j= 1,3),i=1,3)
        ENDDO
     ENDIF
     CALL mp_bcast(zeu,ionode_id, world_comm)
  ELSE
     zeu  (:,:,:) = 0.d0
     epsil(:,:) = 0.d0
  ENDIF
  !
  IF (ionode) READ (1,*) nr1,nr2,nr3
  CALL mp_bcast(nr1,ionode_id, world_comm)
  CALL mp_bcast(nr2,ionode_id, world_comm)
  CALL mp_bcast(nr3,ionode_id, world_comm)
  !
  !  read REAL-space interatomic force constants
  !
  ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
  frc(:,:,:,:,:,:,:) = 0.d0
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           DO nb= 1,nat
              IF (ionode) READ (1,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid,ionode_id, world_comm)
              CALL mp_bcast(jbid,ionode_id, world_comm)
              CALL mp_bcast(nabid,ionode_id, world_comm)
              CALL mp_bcast(nbbid,ionode_id, world_comm)
              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                 na.NE.nabid .OR. nb.NE.nbbid)                      &
                 CALL errore  ('readfc','error in reading',1)
              IF (ionode) READ (1,*) (((m1bid, m2bid, m3bid,        &
                          frc(m1,m2,m3,i,j,na,nb),                  &
                           m1= 1,nr1),m2=1,nr2),m3=1,nr3)
               
              CALL mp_bcast(frc(:,:,:,i,j,na,nb),ionode_id, world_comm)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  IF (ionode) CLOSE(unit= 1)
  !
  RETURN
END SUBROUTINE readfc
!
!-----------------------------------------------------------------------
SUBROUTINE setupmat (q,dyn,nat,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &         loto_2d, & 
     &         epsil,zeu,frc,nr1,nr2,nr3,has_zstar,rws,nrws,na_ifc,f_of_q,fd)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : celldm
  USE rigid,      ONLY : rgd_blk
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER:: nr1, nr2, nr3, nat, nat_blk, nsc, nrws, itau_blk(nat)
  REAL(DP) :: q(3), tau(3,nat), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nat_blk), rws(0:3,nrws),   &
                  frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk)
  REAL(DP) :: tau_blk(3,nat_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  COMPLEX(DP) dyn_blk(3,3,nat_blk,nat_blk), f_of_q(3,3,nat,nat)
  COMPLEX(DP) ::  dyn(3,3,nat,nat)
  LOGICAL :: has_zstar, na_ifc, fd, loto_2d
  !
  ! local variables
  !
  REAL(DP) :: arg
  COMPLEX(DP) :: cfac(nat)
  INTEGER :: i,j,k, na,nb, na_blk, nb_blk, iq
  REAL(DP) :: qp(3), qbid(3,nsc) ! automatic array
  !
  !
  CALL q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  DO iq= 1,nsc
     !
     DO k= 1,3
        qp(k)= q(k) + qbid(k,iq)
     ENDDO
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     CALL frc_blk (dyn_blk,qp,tau_blk,nat_blk,              &
          &              nr1,nr2,nr3,frc,at_blk,bg_blk,rws,nrws,f_of_q,fd)
      IF (has_zstar .and. .not. na_ifc) &
           CALL rgd_blk(nr1,nr2,nr3,nat_blk,dyn_blk,qp,tau_blk,   &
                         epsil,zeu,bg_blk,omega_blk,celldm(1), loto_2d,+1.d0)
           ! LOTO 2D added celldm(1)=alat to passed arguments
     !
     DO na= 1,nat
        na_blk = itau_blk(na)
        DO nb= 1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = CMPLX(COS(arg),SIN(arg),kind=DP)/nsc
           !
        ENDDO ! nb
        !
        DO i= 1,3
           DO j= 1,3
              !
              DO nb= 1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              ENDDO ! nb
              !
           ENDDO ! j
        ENDDO ! i
     ENDDO ! na
     !
  ENDDO ! iq
  !
  RETURN
END SUBROUTINE setupmat
!
!
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN= 10), INTENT(in) :: asr
  INTEGER, INTENT(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), INTENT(in) :: tau(3,nat)
  REAL(DP), INTENT(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(DP) :: zeu_new(3,3,nat)
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     REAL(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  INTEGER :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  INTEGER, allocatable :: ind_v(:,:,:)
  REAL(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We DO so in order
  ! to limit the amount of memory used.
  !
  REAL(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  REAL(DP) :: scal,norm2, sum
  !
  REAL(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  INTEGER :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  REAL(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     CALL errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  ENDIF
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     DO i= 1,3
        DO j= 1,3
           sum=0.0d0
           DO na= 1,nat
              sum = sum + zeu(i,j,na)
           ENDDO
           DO na= 1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           ENDDO
        ENDDO
     ENDDO
     !
     ! Simple Acoustic Sum Rule on force constants in REAL space
     !
     DO i= 1,3
        DO j= 1,3
           DO na= 1,nat
              sum=0.0d0
               DO nb= 1,nat
                  DO n1= 1,nr1
                     DO n2= 1,nr2
                        DO n3= 1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               WRITE(6,*) ' na, i, j, sum = ',na,i,j,sum
            ENDDO
         ENDDO
      ENDDO
      !
      RETURN
      !
   end if
  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typiCALLy, ibrav= 1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis= 1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        CALL errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     ENDIF
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        WRITE(stdout,*) 'asr: rotational axis may be wrong'
     ENDIF
     WRITE(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  ENDIF
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        ENDDO
     ENDDO
  ENDDO
  !
  p=0
  DO i= 1,3
     DO j= 1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)= 1.0d0
        !
     ENDDO
  ENDDO
  !
  if (n.eq.4) then
     DO i= 1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        DO na= 1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        ENDDO
        !
     ENDDO
  ENDIF
  !
  if (n.eq.6) then
     DO i= 1,3
        DO j= 1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           DO na= 1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           ENDDO
           !
        ENDDO
     ENDDO
  ENDIF
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  DO k= 1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     DO q= 1,k-1
        r= 1
        DO izeu_less= 1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        ENDDO
        if (r.ne.0) then
           CALL sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        ENDIF
     ENDDO
     CALL sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     ENDIF
  ENDDO
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  DO k= 1,p
     r= 1
     DO izeu_less= 1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     ENDDO
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        CALL sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     ENDIF
  ENDDO
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  CALL sp_zeu(zeu_w,zeu_w,nat,norm2)
  IF (ionode) WRITE(*,*)
  WRITE(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !WRITE(6,'("Check projection of zeu")')
  !DO k= 1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  CALL sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) WRITE(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !ENDDO
  !
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  DO k= 1,18*nat
     ALLOCATE(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  ENDDO
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           DO nb= 1,nat
              DO n1= 1,nr1
                 DO n2= 1,nr2
                    DO n3= 1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  p=0
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)= 1.0d0
           !
        ENDDO
     ENDDO
  ENDDO
  !
  if (n.eq.4) then
     DO i= 1,3
        DO na= 1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           DO nb= 1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           ENDDO
           !
        ENDDO
     ENDDO
  ENDIF
  !
  if (n.eq.6) then
     DO i= 1,3
        DO j= 1,3
           DO na= 1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              DO nb= 1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              ENDDO
              !
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  ALLOCATE (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           DO nb= 1,nat
              DO n1= 1,nr1
                 DO n2= 1,nr2
                    DO n3= 1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q= 1
                       l= 1
                       DO while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       ENDDO
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)= 1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  ALLOCATE (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  DO k= 1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     DO l= 1,m
        !
        CALL sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        DO r= 1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        ENDDO
     ENDDO
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        ENDIF
     ENDIF
     DO q= 1,k-1
        r= 1
        DO i_less= 1,n_less
           if (u_less(i_less).eq.q) r=0
        ENDDO
        if (r.ne.0) then
           CALL sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec (:,:,:,:,:,:,:)
        ENDIF
     ENDDO
     CALL sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     ENDIF
  ENDDO
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  DO l= 1,m
     CALL sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     DO r= 1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     ENDDO
  ENDDO
  DO k= 1,p
     r= 1
     DO i_less= 1,n_less
        if (u_less(i_less).eq.k) r=0
     ENDDO
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        CALL sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     ENDIF
     DEALLOCATE(u(k) % vec)
  ENDDO
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  CALL sp1(w,w,nr1,nr2,nr3,nat,norm2)
  WRITE(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !WRITE(6,'("Check projection IFC")')
  !DO l= 1,m
  !  CALL sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) WRITE(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !ENDDO
  !DO k= 1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  CALL sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) WRITE(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  DEALLOCATE(u(k) % vec)
  !ENDDO
  !
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           DO nb= 1,nat
              DO n1= 1,nr1
                 DO n2= 1,nr2
                    DO n3= 1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE (x, w)
  DEALLOCATE (v, ind_v)
  DEALLOCATE (frc_new)
  !
  RETURN
end subroutine set_asr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  INTEGER i,j,na,nat
  REAL(DP) zeu_u(3,3,nat)
  REAL(DP) zeu_v(3,3,nat)
  REAL(DP) scal
  !
  !
  scal=0.0d0
  DO i= 1,3
    DO j= 1,3
      DO na= 1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  INTEGER nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  REAL(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  REAL(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  REAL(DP) scal
  !
  !
  scal=0.0d0
  DO i= 1,3
    DO j= 1,3
      DO na= 1,nat
        DO nb= 1,nat
          DO n1= 1,nr1
            DO n2= 1,nr2
              DO n3= 1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as EXPlained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  INTEGER nr1,nr2,nr3,i,nat
  REAL(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  INTEGER ind_v(2,7)
  REAL(DP) v(2)
  REAL(DP) scal
  !
  !
  scal=0.0d0
  DO i= 1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  ENDDO
  !
  RETURN
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  INTEGER nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  REAL(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  REAL(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  REAL(DP) scal
  !
  !
  scal=0.0d0
  DO j= 1,3
    DO nb= 1,nat
      DO n1= 1,nr1
        DO n2= 1,nr2
          DO n3= 1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !
end subroutine sp3
!
!-----------------------------------------------------------------------
SUBROUTINE q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: nsc
  REAL(DP) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  INTEGER, PARAMETER:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  REAL(DP), PARAMETER:: eps= 1.0d-7
  INTEGER :: i, j, k,i1, i2, i3, idum(nrm), iq
  REAL(DP) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  LOGICAL lbho
  !
  i = 0
  DO i1=-nr1,nr1
     DO i2=-nr2,nr2
        DO i3=-nr3,nr3
           i = i + 1
           DO j= 1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           ENDDO ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           DO j= 1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           ENDDO ! j
           !
           idum(i) = 1
           !
        ENDDO ! i3
     ENDDO ! i2
  ENDDO ! i1
  !
  DO i= 1,nrm-1
     IF (idum(i).EQ.1) THEN
        DO j=i+1,nrm
           IF (idum(j).EQ.1) THEN
              lbho=.TRUE.
              DO k= 1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.AND. (ABS(NINT(delta)-delta).LT.eps)
              ENDDO ! k
              IF (lbho) THEN
                 IF(qnorm(i).GT.qnorm(j)) THEN
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 ENDIF
                 idum(j) = 0
              ENDIF
           ENDIF
        ENDDO ! j
     ENDIF
  ENDDO ! i
  !
  iq = 0
  DO i= 1,nrm
     IF (idum(i).EQ.1) THEN
        iq=iq+1
        qbid(1,iq)= bg_blk(1,1)*qbd(1,i) +  &
                    bg_blk(1,2)*qbd(2,i) +  &
                    bg_blk(1,3)*qbd(3,i)
        qbid(2,iq)= bg_blk(2,1)*qbd(1,i) +  &
                    bg_blk(2,2)*qbd(2,i) +  &
                    bg_blk(2,3)*qbd(3,i)
        qbid(3,iq)= bg_blk(3,1)*qbd(1,i) +  &
                    bg_blk(3,2)*qbd(2,i) +  &
                    bg_blk(3,3)*qbd(3,i)
     ENDIF
  ENDDO ! i
  !
  IF (iq.NE.nsc) CALL errore('q_gen',' probably nr1,nr2,nr3 too small ', iq)
  RETURN
END SUBROUTINE q_gen
!
!-----------------------------------------------------------------------
SUBROUTINE check_at(at,bg_blk,alat,omega)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at(3,3), bg_blk(3,3), alat, omega
  REAL(DP) :: work(3,3)
  INTEGER :: i,j
  REAL(DP), PARAMETER :: small= 1.d-6
  !
  work(:,:) = at(:,:)
  CALL cryst_to_cart(3,work,bg_blk,-1)
  !
  DO j= 1,3
     DO i = 1,3
        IF ( ABS(work(i,j)-NINT(work(i,j))) > small) THEN
           WRITE (stdout,'(3f9.4)') work(:,:)
           CALL errore ('check_at','at not multiple of at_blk',1)
        ENDIF
     ENDDO
  ENDDO
  !
  omega =alat**3 * ABS(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
                       at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
                       at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
  !
  RETURN
END SUBROUTINE check_at
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau (nat, nat_blk, at, at_blk, tau, tau_blk, &
     ityp, ityp_blk, itau_blk)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nat, nat_blk,ityp(nat),ityp_blk(nat_blk), itau_blk(nat)
  REAL(DP) at(3,3),at_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(DP) bg(3,3), r(3) ! work vectors
  INTEGER i,i1,i2,i3,na,na_blk
  REAL(DP) small
  INTEGER NN1,NN2,NN3
  PARAMETER (NN1=8, NN2=8, NN3=8, small= 1.d-8)
  !
  CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  na = 0
  !
  DO i1 = -NN1,NN1
     DO i2 = -NN2,NN2
        DO i3 = -NN3,NN3
           r(1) = i1*at_blk(1,1) + i2*at_blk(1,2) + i3*at_blk(1,3)
           r(2) = i1*at_blk(2,1) + i2*at_blk(2,2) + i3*at_blk(2,3)
           r(3) = i1*at_blk(3,1) + i2*at_blk(3,2) + i3*at_blk(3,3)
           CALL cryst_to_cart(1,r,bg,-1)
           !
           IF ( r(1).GT.-small .AND. r(1).LT.1.d0-small .AND.          &
                r(2).GT.-small .AND. r(2).LT.1.d0-small .AND.          &
                r(3).GT.-small .AND. r(3).LT.1.d0-small ) THEN
              CALL cryst_to_cart(1,r,at,+1)
              !
              DO na_blk= 1, nat_blk
                 na = na + 1
                 IF (na.GT.nat) CALL errore('set_tau','too many atoms',na)
                 tau(1,na)    = tau_blk(1,na_blk) + r(1)
                 tau(2,na)    = tau_blk(2,na_blk) + r(2)
                 tau(3,na)    = tau_blk(3,na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              ENDDO
              !
           ENDIF
           !
        ENDDO
     ENDDO
  ENDDO
  !
  IF (na.NE.nat) CALL errore('set_tau','too few atoms: increase NNs',na)
  !
  RETURN
END SUBROUTINE set_tau
!
!-----------------------------------------------------------------------
SUBROUTINE read_tau &
     (nat, nat_blk, ntyp, bg_blk, tau, tau_blk, ityp, itau_blk)
  !---------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode_id, ionode
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER nat, nat_blk, ntyp, ityp(nat),itau_blk(nat)
  REAL(DP) bg_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(DP) r(3) ! work vectors
  INTEGER i,na,na_blk
  !
  REAL(DP) small
  PARAMETER ( small = 1.d-6 )
  !
  DO na= 1,nat
     IF (ionode) READ(5,*) (tau(i,na),i= 1,3), ityp(na)
     CALL mp_bcast(tau(:,na),ionode_id, world_comm)
     CALL mp_bcast(ityp(na),ionode_id, world_comm)
     IF (ityp(na).LE.0 .OR. ityp(na) .GT. ntyp) &
          CALL errore('read_tau',' wrong atomic type', na)
     DO na_blk= 1,nat_blk
        r(1) = tau(1,na) - tau_blk(1,na_blk)
        r(2) = tau(2,na) - tau_blk(2,na_blk)
        r(3) = tau(3,na) - tau_blk(3,na_blk)
        CALL cryst_to_cart(1,r,bg_blk,-1)
        IF (ABS( r(1)-NINT(r(1)) ) .LT. small .AND.                 &
            ABS( r(2)-NINT(r(2)) ) .LT. small .AND.                 &
            ABS( r(3)-NINT(r(3)) ) .LT. small ) THEN
           itau_blk(na) = na_blk
           go to 999
        ENDIF
     ENDDO
     CALL errore ('read_tau',' wrong atomic position ', na)
999  CONTINUE
  ENDDO
  !
  RETURN
END SUBROUTINE read_tau
!
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
subroutine setgam (q, gam, nat, at,bg,tau,itau_blk,nsc,alat, &
     &             gam_blk, nat_blk, at_blk,bg_blk,tau_blk,omega_blk, &
     &             frcg, nr1,nr2,nr3, rws,nrws, fd)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi
  implicit none
  !
  ! I/O variables
  !
  INTEGER        :: nr1, nr2, nr3, nat, nat_blk,  &
                    nsc, nrws, itau_blk(nat)
  REAL(DP)       :: q(3), tau(3,nat), at(3,3), bg(3,3), alat, rws(0:3,nrws)
  REAL(DP)       :: tau_blk(3,nat_blk), at_blk(3,3), bg_blk(3,3), omega_blk, &
                    frcg(nr1,nr2,nr3,3,3,nat_blk,nat_blk)
  COMPLEX(DP)    :: gam_blk(3,3,nat_blk,nat_blk),f_of_q(3,3,nat,nat)
  COMPLEX(DP)    :: gam(3,3,nat,nat)
  LOGICAL        :: fd
  !
  ! local variables
  !
  REAL(DP)        :: arg
  COMPLEX(DP)     :: cfac(nat)
  INTEGER         :: i,j,k, na,nb, na_blk, nb_blk, iq
  REAL(DP)        :: qp(3), qbid(3,nsc) ! automatic array
  !
  !
  CALL q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  f_of_q=(0.0_DP,0.0_DP)
  DO iq= 1,nsc
     !
     DO k= 1,3
        qp(k)= q(k) + qbid(k,iq)
     ENDDO
     !
     gam_blk(:,:,:,:) = (0.d0,0.d0)
     CALL frc_blk (gam_blk,qp,tau_blk,nat_blk,              &
                   nr1,nr2,nr3,frcg,at_blk,bg_blk,rws,nrws,f_of_q,fd)
     !
     DO na= 1,nat
        na_blk = itau_blk(na)
        DO nb= 1,nat
           nb_blk = itau_blk(nb)
           !
           arg = tpi * ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = CMPLX(cos(arg),sin(arg), kind=dp)/nsc
           !
        end DO ! nb
        DO nb= 1,nat
           DO i= 1,3
              DO j= 1,3
                 nb_blk = itau_blk(nb)
                 gam(i,j,na,nb) = gam(i,j,na,nb) + cfac(nb) * &
                     gam_blk(i,j,na_blk,nb_blk)
              end DO ! j
              end DO ! i
        end DO ! nb
     end DO ! na
     !
  end DO ! iq
  !
  RETURN
end subroutine setgam
!
!--------------------------------------------------------------------
!
subroutine readfg ( ifn, nr1, nr2, nr3, nat, frcg )
  !-----------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  implicit none
  ! I/O variable
  INTEGER, INTENT(in) ::  nr1,nr2,nr3, nat
  REAL(DP), INTENT(out) :: frcg(nr1,nr2,nr3,3,3,nat,nat)
  ! local variables
  INTEGER i, j, na, nb, m1,m2,m3, ifn
  INTEGER ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  !
  !
  IF (ionode) READ (ifn,*) m1, m2, m3
  CALL mp_bcast(m1, ionode_id, world_comm)
  CALL mp_bcast(m2, ionode_id, world_comm)
  CALL mp_bcast(m3, ionode_id, world_comm)
  if ( m1 /= nr1 .or. m2 /= nr2 .or. m3 /= nr3) &
       CALL errore('readfg','inconsistent nr1, nr2, nr3 read',1)
  DO i= 1,3
     DO j= 1,3
        DO na= 1,nat
           DO nb= 1,nat
              IF (ionode) read (ifn,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid, ionode_id, world_comm)
              CALL mp_bcast(jbid, ionode_id, world_comm)
              CALL mp_bcast(nabid, ionode_id, world_comm)
              CALL mp_bcast(nbbid, ionode_id, world_comm)
              
              if(i.ne.ibid.or.j.ne.jbid.or.na.ne.nabid.or.nb.ne.nbbid)  then
                  WRITE(stdout,*) i,j,na,nb,'  <>  ', ibid, jbid, nabid, nbbid
                  CALL errore  ('readfG','error in reading',1)
              else
                  IF (ionode) read (ifn,*) (((m1bid, m2bid, m3bid,     &
                                 frcg(m1,m2,m3,i,j,na,nb), &
                                 m1= 1,nr1),m2=1,nr2),m3=1,nr3)
              ENDIF
              CALL mp_bcast(frcg(:,:,:,i,j,na,nb), ionode_id, world_comm)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  IF (ionode) CLOSE(ifn)
  !
  RETURN
end subroutine readfg
!
!
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, num_rap_mode, nspin_mag )

  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, sr, ft, irt, nsym, nrot, t_rev, time_reversal,&
                         sname, copy_sym, s_axis_to_cart

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, ntyp, nspin_mag
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3,nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN) :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)
  REAL(DP) :: gi (3, 48), gimq (3), sr_is(3,3,48), rtau(3,48,nat)
  INTEGER :: irotmq, nsymq, nsym_is, isym, i, ierr
  LOGICAL :: minus_q, search_sym, sym(48), magnetic_sym
!
!  find the small group of q
!
  time_reversal=.TRUE.
  IF (.NOT.time_reversal) minus_q=.FALSE.

  sym(1:nsym)=.true.
  CALL smallg_q (xq, 0, at, bg, nsym, s, sym, minus_q)
  nsymq=copy_sym(nsym,sym )
  CALL s_axis_to_cart ()
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!
!  if the small group of q is non symmorphic,
!  search the symmetries only if there are no G such that Sq -> q+G
!
  search_sym=.TRUE.
  IF ( ANY ( ABS(ft(:,1:nsymq)) > 1.0d-8 ) ) THEN
     DO isym= 1,nsymq
        search_sym=( search_sym.and.(ABS(gi(1,isym))<1.d-8).and.  &
                                    (ABS(gi(2,isym))<1.d-8).and.  &
                                    (ABS(gi(3,isym))<1.d-8) )
     ENDDO
  ENDIF
!
!  Set the representations tables of the small group of q and
!  find the mode symmetry
!
  IF (search_sym) THEN
     magnetic_sym=(nspin_mag==4)
     CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
     sym (1:nsym) = .TRUE.
     CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
     CALL find_mode_sym_new (u, w2, tau, nat, nsymq, s, sr, irt, xq,    &
             rtau, amass, ntyp, ityp, 1, .FALSE., .FALSE., num_rap_mode, ierr)

  ENDIF
  RETURN
  END SUBROUTINE find_representations_mode_q

!mz adds this routine
SUBROUTINE diffuse_scattering(nq, nq_super, nq_strft, nat, ntyp, amass, ityp, q, w2, &
               z_nq_zg, q_external, zero_one_phonon, full_phonon, atmsf_a, atmsf_b, & 
               dim1, dim2, dim3, tau, alat, atm, atom_resolved, mode_resolved, &
               ntypx, at, q_in_cryst_coord, T, eps2, print_raw_data, &
               nrots, kres1, kres2, kmin, kmax, col1, col2, Np) 
  ! we start here with the WRITE_eigenvectors.f90 routine
  USE kinds,      ONLY : dp
  USE constants,  ONLY : amu_ry, ry_to_thz, ry_to_cmm1, H_PLANCK_SI, &  
                         K_BOLTZMANN_SI, AMU_SI, pi, BOHR_RADIUS_ANGS
  USE cell_base,  ONLY : bg
  USE io_global,  ONLY : ionode, ionode_id
  USE mp_global,  ONLY : inter_pool_comm
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,   ONLY : world_comm
  ! 
  IMPLICIT NONE
  ! input
  LOGICAL, INTENT(in)          :: q_in_cryst_coord, q_external, print_raw_data
  LOGICAL, INTENT(in)          :: zero_one_phonon, full_phonon
  LOGICAL, INTENT(in)          :: atom_resolved, mode_resolved
  INTEGER, INTENT(in)          :: dim1, dim2, dim3, T 
  INTEGER, INTENT(in)          :: nq, nat, ntyp, ntypx, nq_super, nq_strft
  INTEGER, INTENT(in)          :: nrots, kres1, kres2, col1, col2, Np
  REAL(DP), INTENT(in)         :: kmin, kmax
  REAL(DP), INTENT(in)         :: alat, eps2
  REAL(DP), INTENT(in)         :: at(3, 3)
  REAL(DP), intent(in)         :: atmsf_a(ntypx, 5), atmsf_b(ntypx, 5)
  REAL(DP), INTENT(in)         :: q(3, nq), w2(3 * nat, nq), amass(ntyp), tau(3, nat)
  ! nq is the number of qpoints in sets A and B
  CHARACTER(LEN=3), INTENT(in) :: atm(ntypx)
  INTEGER ityp(nat)
  COMPLEX(DP), INTENT(in)      :: z_nq_zg(3 * nat, 3 * nat, nq)
  !
  ! local
  !
  INTEGER                  :: nat3, na, nta, ipol, i, j, ii 
  INTEGER                  :: nq_tot, qp, qp2, p, p1, ntap, k, kp 
  INTEGER                  :: ctr, ctr2, lower_bnd, upper_bnd, qlistA(nq)
  ! nq_tot total number of q-points (including sets A,B,C)
  REAL(DP)                 :: freq(3 * nat, nq), ph_w(3 * nat, nq),l_q(3 * nat, nq) 
  REAL(DP)                 :: q_A(3), q_B(3), atm_ffct(nat, nq), p_q(3 * nat, nq), e_nq(3 * nat, nq)
  REAL(DP)                 :: hbar, ang, dotp, dotp2, dotp3, dotp4, ph_w_mean, units
  REAL(DP), PARAMETER      :: eps = 1.0d-4
  ! l_q is the amplitude \sigma at temperature T, e_nq --> to calculate total vibrational energy
  ! p_q is the momentum on the nuclei \hbar\2\l_\bq\nu \SQRT(n_{q\nu,T}+1/2)
  !
  COMPLEX(DP)              :: imagi
  COMPLEX(DP)              :: z_zg(3 * nat,3 * nat, nq)
  CHARACTER(len=256)       :: filename, pt_mz, pt_mz2
  ! ALLOCATE TABLES
  REAL(DP),    ALLOCATABLE :: DW_T_fact_kkp(:, :), DW_T_fact_q0_kkp(:, :), DW_T_fact_q0(:)
  REAL(DP),    ALLOCATABLE :: sigma_DW(:, :, :), Q_ppkkaa(:, :, :)
  REAL(DP),    ALLOCATABLE :: strf_kkp(:, :, :), strf(:), strf_j(:, :)
  REAL(DP),    ALLOCATABLE :: strf_q0(:), strf_q0_kkp(:, :, :), strf_rot(:,:)
  COMPLEX(DP), ALLOCATABLE :: strf_full(:), strf_full_kk(:, :, :)
  ! for displacements
  REAL(DP), ALLOCATABLE    :: Rlist(:, :)
  !
  !  
  hbar   = 0.5 * H_PLANCK_SI / pi ! reduce Plnack constant
  units  = DBLE(2.d0 * pi / alat / BOHR_RADIUS_ANGS)
  ang    = 1.0E-10            ! angstrom units
  imagi  = (0.0d0, 1.0d0) !imaginary unit
  ! Inititialize eigenvectors matrix
  z_zg   = (0.d0, 0.d0)
  ! Set intitial values
  nq_tot = dim1 * dim2 * dim3
  nat3   = 3 * nat
  !
  !
  freq = 0.d0
  CALL fkbounds( nq, lower_bnd, upper_bnd )
  DO qp =  lower_bnd, upper_bnd ! 1, nq
  ! convert eigenvectors to mass-unsCALLed
    DO i = 1, nat3
      DO na = 1, nat
         nta = ityp(na)
        DO ipol = 1, 3
           z_zg((na - 1) * 3 + ipol, i, qp) = z_nq_zg((na - 1) * 3 + ipol, i, qp) * SQRT(amu_ry * amass(nta))
        ENDDO
     ENDDO
    ENDDO
  !
    DO i = 1, nat3 
  !    if (w2(i, qp) < 0.0) freq(i,qp) = -freq(i,qp)
      freq(i, qp)= SQRT(ABS(w2(i, qp)))
    ENDDO
  ENDDO ! qp
  CALL mp_sum(freq, inter_pool_comm)
  CALL mp_sum(z_zg, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  IF (ionode) THEN
   WRITE(*,*) nq, nat3
   DO i = 1, nat3
    DO qp =  1, nq
      IF (w2(i, qp) < 0.0 .AND. ionode) WRITE(*,*) "WARNING: Negative frequencies", w2(i, qp)
    ENDDO
   ENDDO
  ENDIF
  ! Einstein model set polarization vectors to unity
  !
  ! z_zg = (1.d0, 1.d0)/DBLE(SQRT(real(2.d0*nat3))) ! For Einstein model
  !
  ph_w = freq * ry_to_thz * (1.0E12) * 2 * pi ! correct frequency for phonons in SI
  ph_w_mean = 0.d0
  DO qp = 1, nq
    DO i = 1, nat3
      ph_w_mean = ph_w_mean + ph_w(i,qp) 
    ENDDO
  ENDDO
  ! Einstein model mean frequency
  ph_w_mean = DBLE(ph_w_mean / nq / nat3) 
  !ph_w(:,:) = ph_w_mean ! For Einstein model
  IF (ionode) WRITE(*,*)
  IF (ionode) WRITE(*,'(A, 1F10.6)') "Mean frequency (THz):", ph_w_mean / 1.0E12 / 2 / pi
  !
  !
  ! set amplitudes of displacements l_q = \sigma_\bq\nu and momenta 
  IF (q_external) THEN
      IF (q_in_cryst_coord .EQV. .FALSE.) THEN
      ! in both cases convert them to crystal 
          CALL cryst_to_cart(nq, q, at, -1)
      ELSE
          CALL cryst_to_cart(nq, q, at, -1)
      ENDIF
    ELSE
    CALL cryst_to_cart(nq, q, at, -1)
  ENDIF
  !
  l_q = 0.d0
  p_q = 0.d0
  CALL fkbounds( nq, lower_bnd, upper_bnd )
  DO qp =  lower_bnd, upper_bnd ! 1, nq
    DO i = 1, nat3 
      IF (w2(i, qp) .lt. 0.0d0 + eps2) THEN
          l_q(i, qp) = 0.d0
          p_q(i, qp) = 0.d0
      ELSE
          l_q(i, qp) = SQRT(hbar / ph_w(i, qp) / 2.0d0 / AMU_SI / ang**2.0d0) & 
                     * SQRT(1.0d0 + 2.0d0 / (EXP(hbar * ph_w(i,qp)/ (K_BOLTZMANN_SI * T)) - 1)) 
      ENDIF
       ! l_qn SQRT(2n_qn +1) = sigma_qn
    ENDDO
   ! First set all displ. amplitudes to zero for accoustic modes for every G 
    IF ( (( ABS(MODULO(ABS(q(1, qp)),1.0d0)) .LT. eps) .OR. & 
          ( ABS(MODULO(ABS(q(1, qp)),1.0d0) - 1.0d0) .LT. eps)) .AND. &
         (( ABS(MODULO(ABS(q(2, qp)),1.0d0)) .LT. eps) .OR. & 
          ( ABS(MODULO(ABS(q(2, qp)),1.0d0) - 1.0d0) .LT. eps)) .AND. &
         (( ABS(MODULO(ABS(q(3, qp)),1.0d0)) .LT. eps) .OR. & 
          ( ABS(MODULO(ABS(q(3, qp)),1.0d0)-1.0d0) .LT. eps)) ) THEN
            l_q(1, qp) = 0.0d0
            l_q(2, qp) = 0.0d0
            l_q(3, qp) = 0.0d0
            p_q(1, qp) = 0.0d0
            p_q(2, qp) = 0.0d0
            p_q(3, qp) = 0.0d0
    ENDIF
  ENDDO ! qp
  CALL mp_sum(l_q, inter_pool_comm)
  CALL mp_sum(p_q, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  ! 
  ! To determine which points belong in set A
  qlistA = 0
  CALL fkbounds( nq, lower_bnd, upper_bnd )
    DO qp =  lower_bnd, upper_bnd ! 1, nq
      q_A(1) = q(1, qp) +  q(1, qp)
      q_A(2) = q(2, qp) +  q(2, qp)
      q_A(3) = q(3, qp) +  q(3, qp)
    ! 
    !  q_A = q(:, qp) +  q(:, qp)
      IF ( (( ABS(MODULO(ABS(q_A(1)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(1)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(2)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(2)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(3)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(3)), 1.0d0) - 1.0d0) .LT. eps)) ) THEN
            qlistA(qp) = 1
            !WRITE(*,*) "qAAA", q(:, qp)
      ENDIF
    ENDDO
  CALL mp_sum(qlistA, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !  
  !
  ! Exponent of DW factors 
  ALLOCATE(sigma_DW(nat, 3, 3))
  ! sigma_DW: thermal displacement tensor
  CALL fkbounds( nq - nq_strft, lower_bnd, upper_bnd )
  !
  sigma_DW = 0.0d0
  DO qp = lower_bnd, upper_bnd  ! for q-points in supercell 
  ! DO qp = 1, nq_super
    DO k= 1, nat
      nta = ityp(k)
      DO i= 1, 3  ! i is for cart directions
        DO p = 1, 3 ! p is also for cartesian
          DO j = 1, nat3 ! modes \nu
            IF (qlistA(qp) .EQ. 1) THEN
              sigma_DW(k, i, p) = sigma_DW(k, i, p) + 1.0 / DBLE(nq_tot) / DBLE(2.0d0 * amass(nta)) * & 
                                  DBLE(z_zg((k - 1) * 3 + i, j, qp) * & 
                                  CONJG(z_zg((k - 1) * 3 + p, j, qp))) * l_q(j, qp)**2
            ELSE 
              sigma_DW(k, i, p) = sigma_DW(k, i, p) + 1.0/DBLE(nq_tot) / DBLE(amass(nta)) * & 
                                  DBLE(z_zg((k - 1) * 3 + i, j, qp) * & 
                                  CONJG(z_zg((k - 1) * 3 + p, j, qp))) * l_q(j, qp)**2
               ! an extra factor of 2 is need because we have only points in set B
            ENDIF
          ENDDO ! j
        ENDDO ! p
      ENDDO ! i
    ENDDO ! k
  ENDDO ! qp
  !
  CALL mp_sum(sigma_DW, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  IF (ionode) THEN
    WRITE(*,*) "Mean-squared displacement (Ang**2)"
    DO k= 1, nat
      WRITE(*,*) "atom:", k
      DO i= 1,3  ! i is for cart directions
         WRITE(*,'(3F14.8)') (2.0d0 * sigma_DW(k, i, p), p = 1, 3)
      ENDDO ! i
    ENDDO ! k
  ENDIF
  !!!!!!!!!!!!!
  IF (ionode) THEN
    atm_ffct = 0.d0
    DO k = 1, nat
      nta = ityp(k)
      DO qp = nq_super + 1, nq !1, qpts_strf
         DO ii = 1, 5
           q_A = q(:, qp) 
           CALL cryst_to_cart(1, q_A, bg, +1)
           q_A = q_A * units 
           atm_ffct(k, qp) = atm_ffct(k, qp) + (atmsf_a(nta, ii) * & 
                             EXP(-atmsf_b(nta, ii) * (NORM2(q_A) / 4.d0 / pi)**2))
         ENDDO
      ENDDO
    ENDDO
  ENDIF !(ionode)
  CALL mp_bcast(atm_ffct, ionode_id, world_comm)
  !
  IF (zero_one_phonon) THEN
    ALLOCATE(DW_T_fact_q0_kkp(nq, nat),DW_T_fact_q0(nq))
    ALLOCATE(strf_q0(nq))
    IF (atom_resolved) ALLOCATE(strf_q0_kkp(nq, nat, nat))
     DW_T_fact_q0_kkp = 0.0d0
     DW_T_fact_q0 = 0.0d0
     strf_q0 = 0.0d0
    !  
    CALL fkbounds( nq - nq_super, lower_bnd, upper_bnd )
    DO qp = lower_bnd + nq_super, upper_bnd + nq_super
      q_A(1) = q(1, qp) 
      q_A(2) = q(2, qp) 
      q_A(3) = q(3, qp) 
      !
      IF ( (( ABS(MODULO(ABS(q_A(1)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(1)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(2)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(2)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(3)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(3)), 1.0d0) - 1.0d0) .LT. eps)) ) THEN
          ! 
          CALL cryst_to_cart(1, q(:, qp), bg, +1)
          DO k= 1 , nat ! k represents the atom
             nta = ityp(k)       
                   dotp = 0.0d0
                   DO i= 1, 3  ! i is for cart directions
                    DO p = 1, 3 ! p is also for cartesian
                      dotp = dotp + q(i, qp) * q(p, qp) * sigma_DW(k, i, p) * units**2
                     ! this will be COMPLEX otherwise
                    ENDDO ! p
                   ENDDO ! i
                   DW_T_fact_q0_kkp(qp, k) = DW_T_fact_q0_kkp(qp, k) + dotp ! DW factor for each q
          ENDDO ! k
          CALL cryst_to_cart(1, q(:,qp), at, -1)
      ENDIF
    ENDDO ! qp
  CALL mp_sum(DW_T_fact_q0_kkp, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
    !
    ! evaluate structure factor for zero-phonon contribution
  DO qp = lower_bnd + nq_super, upper_bnd + nq_super
      q_A(1) = q(1, qp) 
      q_A(2) = q(2, qp) 
      q_A(3) = q(3, qp) 
      !
      IF ( (( ABS(MODULO(ABS(q_A(1)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(1)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(2)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(2)), 1.0d0) - 1.0d0) .LT. eps)) .AND. &
           (( ABS(MODULO(ABS(q_A(3)), 1.0d0)) .LT. eps) .OR. &
            ( ABS(MODULO(ABS(q_A(3)), 1.0d0) - 1.0d0) .LT. eps)) ) THEN
      CALL cryst_to_cart(1, q(:, qp), bg, +1)
      DO k = 1, nat ! k for the atom
         DO kp = 1, nat ! kp for the atom
               dotp = 0.0d0
                DO ii = 1, 3
                   dotp = dotp + q(ii, qp) * (tau(ii, k) - tau(ii, kp))! 
                ENDDO
               IF (atom_resolved) THEN
                 strf_q0_kkp(qp, k, kp) = EXP(-DW_T_fact_q0_kkp(qp, k))*EXP(-DW_T_fact_q0_kkp(qp, kp))& 
                                        * nq_tot * nq_tot * COS(2.d0 * pi * dotp) * atm_ffct(k, qp) * atm_ffct(kp, qp)
                                         !*nq_tot*nq_tot*EXP(imagi*2*pi*dotp) * atm_ffct(k, qp) * atm_ffct(kp, qp)
               ENDIF
               strf_q0(qp) = strf_q0(qp) + EXP(-DW_T_fact_q0_kkp(qp,k))*EXP(-DW_T_fact_q0_kkp(qp,kp))& 
                           * nq_tot * nq_tot * COS(2.d0 * pi * dotp) * atm_ffct(k, qp) * atm_ffct(kp, qp)  
          ENDDO ! kp
      ENDDO ! k
      CALL cryst_to_cart(1, q(:,qp), at, -1)
      ENDIF
  ENDDO ! qp
  !
  IF (atom_resolved) CALL mp_sum(strf_q0_kkp, inter_pool_comm)
  !
  CALL mp_sum(strf_q0, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  ENDIF ! if zero_one_phonon
  ! 
  ! convert to Cartesian
  CALL cryst_to_cart(nq, q, bg, +1)
  ! now for structure factor one and full phonon contribution
  ALLOCATE(DW_T_fact_kkp(nq, nat)) 
  DW_T_fact_kkp = 0.0d0!(0.0d0, 0.0d0)
  CALL fkbounds( nq - nq_super, lower_bnd, upper_bnd )
  DO qp = lower_bnd + nq_super, upper_bnd + nq_super
    DO k = 1, nat ! k represents the atom
      dotp = 0.0d0
      DO i= 1, 3  ! i is for Cart. directions
        DO p = 1, 3 ! p is also for Cart. directions
          dotp = dotp + q(i, qp) * q(p, qp) * sigma_DW(k, i, p) * units**2
          ctr2 = ctr2 + 1
        ENDDO ! p
      ENDDO ! i
      DW_T_fact_kkp(qp, k) = DW_T_fact_kkp(qp, k) + dotp
    ENDDO ! k
  ENDDO ! qp
  !
  CALL mp_sum(DW_T_fact_kkp, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  DEALLOCATE(sigma_DW)
  ! 
  ! compute one_phonon contribution
  IF (zero_one_phonon) THEN
  ! Here q_pts should be in cartessian
    ALLOCATE(strf_j(nq,nat3))
    IF (atom_resolved)  ALLOCATE(strf_kkp(nq, nat, nat))
    !WRITE(*,*) "mz6"
    IF (atom_resolved) strf_kkp = 0.0d0
    strf_j =  0.0d0
  
    CALL fkbounds( nq - nq_super, lower_bnd, upper_bnd )
    DO qp = lower_bnd + nq_super, upper_bnd + nq_super
      DO k = 1, nat ! k represents the atom
        nta = ityp(k)       
        DO kp = 1, nat ! k represents the atom
          ntap = ityp(kp)       
          dotp2 = 0.0d0
          DO ii = 1, 3
            dotp2 = dotp2 + q(ii, qp) * (tau(ii, k) - tau(ii, kp))! 
          ENDDO
          DO j = 1, nat3 ! modes \nu
            dotp = 0.0d0
            DO i = 1, 3  ! i is for cart directions
              DO p = 1, 3 ! p is also for cartesian
                dotp = dotp + q(i, qp) * q(p, qp) * units**2 * & ! 
                     DBLE(z_zg((k - 1) * 3 + i, j, qp) * CONJG(z_zg((kp - 1) * 3 + p, j, qp)) * & 
                     EXP(-imagi * 2.0d0 * pi * dotp2)) 
                       
              ENDDO ! p
            ENDDO ! i
            IF (atom_resolved) THEN 
               strf_kkp(qp, k, kp) = strf_kkp(qp, k, kp) + & 
                                     EXP(-DW_T_fact_kkp(qp, k)) * EXP(-DW_T_fact_kkp(qp, kp)) * & 
                                     nq_tot / DBLE(SQRT(amass(nta) * amass(ntap))) * dotp * &
                                     l_q(j, qp)**2  * atm_ffct(k, qp) * atm_ffct(kp, qp) 
            ENDIF
            !
               strf_j(qp, j) = strf_j(qp, j) + EXP(-DW_T_fact_kkp(qp, k)) * EXP(-DW_T_fact_kkp(qp, kp)) * & 
                               nq_tot / DBLE(SQRT(amass(nta) * amass(ntap))) * dotp * & 
                               l_q(j, qp)**2 * atm_ffct(k, qp) * atm_ffct(kp, qp) 
          ENDDO ! j
        ENDDO ! kp
      ENDDO ! k
    ENDDO ! qp
    IF (atom_resolved) CALL mp_sum(strf_kkp, inter_pool_comm)
    CALL mp_sum(strf_j, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
  ENDIF ! zero_one_phonon
  !
  !
  IF (zero_one_phonon) THEN
    ALLOCATE(strf(nq))
    strf = 0.0d0
    DO qp = nq_super + 1, nq
       DO j = 1,nat3 
          strf(qp) = strf(qp) + strf_j(qp, j)
          !
          strf_j(qp, j) = strf_j(qp, j) + strf_q0(qp)
       ENDDO
       strf(qp) = strf(qp) + strf_q0(qp)
       IF (atom_resolved) strf_kkp(qp, :, :) = strf_kkp(qp, :, :) + strf_q0_kkp(qp, :, :)
    ENDDO
    ! APPLY BROADENING and print outputs
    ctr = 0
    DO p = nq_super + 1, nq
      IF ((q(col1, p) .GT. 0.d0 - eps) .AND. (q(col2, p) .GT. 0.d0 - eps)) THEN
        ctr = ctr + 1
      ENDIF
    ENDDO
    ALLOCATE(strf_rot(ctr * nrots,4))
    !
    strf_rot = 0.d0
    IF (ionode) CALL rotate(strf, q, nq, nq_super, nrots, &
                            ctr, strf_rot, col1, col2)
    CALL mp_bcast(strf_rot, ionode_id, world_comm)
    CALL disca_broadening(strf_rot, ctr * nrots, kres1, kres2, alat, &
                          kmin, kmax, col1, col2, Np, 'strf_one-ph_broad.dat')
    ! for mode_resolved
    IF (mode_resolved) THEN
      DO j = 1, nat3
        strf_rot = 0.d0
        WRITE( pt_mz,'(i2.2)') j
        filename = 'strf_mode_'//TRIM(pt_mz)//'_broad.dat'
        IF (ionode) CALL rotate(strf_j(:, j), q, nq, nq_super, nrots, & 
                                ctr, strf_rot, col1, col2)
        !
        CALL mp_bcast(strf_rot, ionode_id, world_comm)
        CALL disca_broadening(strf_rot, ctr * nrots, kres1, kres2, alat, &
                              kmin, kmax, col1, col2, Np, filename)
      ENDDO
    ENDIF ! mode_resolved
    !
    IF (atom_resolved) THEN
      DO k = 1, nat
        DO kp = 1, k
        strf_rot = 0.d0
        WRITE( pt_mz,'(i1.1)') k
        WRITE( pt_mz2,'(i1.1)') kp
        filename = 'strf_one_k_' // TRIM(pt_mz) // '_' & 
                                      // TRIM(pt_mz2) // '_broad.dat'
        ! 
        IF (ionode) CALL rotate(strf_kkp(:, k, kp), q, nq, nq_super, nrots, & 
                                ctr, strf_rot, col1, col2)
        !
        CALL mp_bcast(strf_rot, ionode_id, world_comm)
        CALL disca_broadening(strf_rot, ctr * nrots, kres1, kres2, alat, &
                              kmin, kmax, col1, col2, Np, filename)
                             !'zero_one_k_'//TRIM()//'_'//TRIM(pt_mz2)//'_broad.dat')
        ENDDO ! kp
      ENDDO
    ENDIF
    ! 
    DEALLOCATE(strf_rot)
    ! 
    IF (print_raw_data) THEN
      IF (ionode) THEN 
        filename = 'Bragg_scattering.dat'
        OPEN (unit = 90, file = filename, status = 'unknown', form = 'formatted')
        filename = 'strf_q_nu_one-phonon.dat'
        OPEN (unit = 95, file = filename, status = 'unknown', form = 'formatted')
        filename = 'strf_one-phonon.dat'
        OPEN (unit = 97, file = filename, status = 'unknown', form = 'formatted')
        !
        DO qp = nq_super + 1, nq
           WRITE (90, '(40f26.6)') q(:,qp) * units, strf_q0(qp)
           WRITE (95, '(40f26.6)') q(:,qp) * units, ( strf_j(qp, j),  j = 1, nat3)
           WRITE (97, '(40f26.6)') q(:,qp) * units,  strf(qp) 
      ! 
        ENDDO
      !
        CLOSE(90) 
        CLOSE(95)
        CLOSE(97)
        !
        IF (atom_resolved) THEN
          DO k = 1, nat
            WRITE( pt_mz,'(i1.1)') k
            filename = 'strf_one_k_resolved_' // TRIM( pt_mz ) //'.dat' !'.fp'
            OPEN (unit = 66, file = filename, status = 'unknown', form = 'formatted')
            DO qp = nq_super + 1, nq
             WRITE (66, '(40f26.6)') q(:, qp) * units, &
                            ( DBLE((strf_kkp(qp, k, kp) + strf_kkp(qp, kp, k)) / 2.0d0), kp = 1, k)
             ! NOTE: We should take the real part if and only if we sum
             ! [strf_full_kk(qp, 2, 1) + strf_full_kk(qp, 1, 2)]/2.
             ! strf_full_kk(qp, 2, 1) is not necessarilly real
            ENDDO ! qp
            CLOSE(66)
          ENDDO
        ENDIF ! k-resolved
        !
      ENDIF
    ENDIF ! print_raw_data
  ENDIF ! zero_one_phonon 
  !
  !
  IF (zero_one_phonon) THEN
    DEALLOCATE(DW_T_fact_q0_kkp, DW_T_fact_q0)
    DEALLOCATE(strf_j, strf,strf_q0)
    IF (atom_resolved)  DEALLOCATE(strf_kkp, strf_q0_kkp)
  ENDIF 
  ! Compute Full strf
  !
  IF (full_phonon) THEN
    ! Generate lattice vectors in crystal coordinates   
    ALLOCATE(Rlist(3, nq_tot))
    !
    ctr2 = 1
    DO i = 0, dim3 - 1
     DO j = 0, dim2 - 1
      DO  k= 0, dim1 - 1
       Rlist(1, ctr2) = k
       Rlist(2, ctr2) = j
       Rlist(3, ctr2) = i 
       !
       ctr2= ctr2 + 1
      ENDDO
     ENDDO
    ENDDO
    !
    ! Convert to Cartesian
    CALL cryst_to_cart(nq_tot, Rlist, at, +1)
    ! 
    ! Convert to crystal
    ALLOCATE(strf_full(nq))
    IF (atom_resolved) ALLOCATE(strf_full_kk(nq, nat, nat))
    !
    CALL fkbounds( nq_tot, lower_bnd, upper_bnd )
    !
    ALLOCATE(Q_ppkkaa(upper_bnd - lower_bnd + 1, nat3, nat3))
    ! allocate array in parallel / minimize memory burden
    Q_ppkkaa = 0.0d0
    ! DO qp2 = lower_bnd, upper_bnd  ! for q-points in supercell 
    ctr = 1
    DO p1 = lower_bnd, upper_bnd !1, nq_tot
        DO qp2 = 1, nq_super !lower_bnd, upper_bnd  ! for q-points in supercell 
        dotp4 = 0.0d0
        DO ii = 1, 3
          dotp4 = dotp4 + q(ii, qp2) * Rlist(ii, p1) ! to get EXP (iq.(Rp-Rp')) 
        ENDDO
        !
        DO j = 1, nat3 ! modes \nu
        !DO j = 7, 7 ! modes \nu
          DO i = 1, 3  ! i is for cart directions
            DO p = 1, 3 ! p is also for cartesian
              DO k = 1, nat ! k represents the atom
                nta = ityp(k)       
                DO kp = 1, nat ! k represents the atom
                  ntap = ityp(kp)       
                  IF (qlistA(qp2) .EQ. 1) THEN
                    Q_ppkkaa(ctr, (k - 1) * 3 + i, (kp - 1) * 3 + p) = &
                             Q_ppkkaa(ctr, (k - 1) * 3 + i, (kp - 1) * 3 + p) + & 
                             1.0d0 / DBLE(SQRT(amass(nta) * amass(ntap))) / nq_tot * l_q(j, qp2)**2.0d0 * & 
                             DBLE(z_zg((k - 1) * 3 + i, j, qp2) * CONJG(z_zg((kp - 1) * 3 + p, j, qp2)) * & 
                             EXP(imagi * 2.0d0 * pi * dotp4))   
                  ELSE 
                    Q_ppkkaa(ctr, (k - 1) * 3 + i, (kp - 1) * 3 + p) = & 
                             Q_ppkkaa(ctr, (k - 1) * 3 + i, (kp - 1) * 3 + p) + & 
                             2.0d0 / DBLE(SQRT(amass(nta) * amass(ntap))) / nq_tot * l_q(j, qp2)**2.0d0 * & 
                             DBLE(z_zg((k - 1) * 3 + i, j, qp2) * CONJG(z_zg((kp - 1) * 3 + p, j, qp2)) * & 
                             EXP(imagi * 2.0d0 * pi * dotp4))   
                  ENDIF
                ENDDO ! kp 
              ENDDO ! k 
            ENDDO ! p
          ENDDO ! i
        ENDDO ! j
      ENDDO ! qp2
    ctr = ctr + 1
    ENDDO ! p1
    ! 
    strf_full = 0.0d0
    !
    IF (atom_resolved) strf_full_kk = 0.0d0
    !
    ctr = 1
    DO p1 = lower_bnd, upper_bnd !1, nq_tot
       DO qp = nq_super + 1, nq !lower_bnd + nq_super, upper_bnd + nq_super
          dotp3 = 0.0d0
          DO ii = 1, 3
            dotp3 = dotp3 + q(ii, qp) * Rlist(ii, p1)! to get EXP (iS.(Rp-Rp'))
          ENDDO
          DO k = 1, nat ! k represents the atom
           DO kp = 1, nat ! kp represents the atom
           dotp2 = 0.0d0
           DO ii = 1, 3
              dotp2 = dotp2 + q(ii, qp) * (tau(ii, k) - tau(ii, kp))! 
           ENDDO
              dotp = 0.0d0
              DO i = 1, 3  ! i is for cart directions
                DO p = 1, 3 ! p is also for cartesian
                   dotp = dotp +  q(i, qp) * q(p, qp) * units**2.0d0 &  ! dotp is P_p,kk' 
                               * Q_ppkkaa(ctr, (k - 1) * 3 + i, (kp - 1) * 3 + p)
                ENDDO
              ENDDO
              !
              strf_full(qp) = strf_full(qp) + EXP(-DW_T_fact_kkp(qp, k)) * EXP(-DW_T_fact_kkp(qp, kp)) & 
                                   * EXP(dotp) * EXP(imagi * 2.0d0 * pi * dotp2) * EXP(imagi * 2.0d0 * pi * dotp3) & 
                                   * atm_ffct(k, qp) * atm_ffct(kp, qp) * nq_tot 
                                   ! here we can replace EXP(dotp) with 1 + dotp
              IF (atom_resolved) THEN
                strf_full_kk(qp, k, kp) = strf_full_kk(qp, k, kp) & 
                                   + EXP(-DW_T_fact_kkp(qp, k)) * EXP(-DW_T_fact_kkp(qp, kp)) &
                                   * EXP(dotp) * EXP(imagi * 2.0d0 * pi * dotp2) * EXP(imagi * 2.0d0 * pi * dotp3) &
                                   * atm_ffct(k, qp) * atm_ffct(kp, qp) * nq_tot 
              ENDIF
            ENDDO ! kp
          ENDDO ! k
        ENDDO ! qp
      ctr = ctr + 1
    ENDDO ! p1
    CALL mp_sum( strf_full, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    IF (atom_resolved) CALL mp_sum( strf_full_kk, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    ! To allocate a new matrix without negative values for q coordinates 
    ctr2 = 0
    DO p = nq_super + 1, nq
      IF ((q(col1, p) .GT. 0.d0 - eps) .AND. (q(col2, p) .GT. 0.d0 - eps)) THEN
        ctr2 = ctr2 + 1
      ENDIF
    ENDDO
    ! To rotate, apply broadening, and print the map
    ALLOCATE(strf_rot(ctr2 * nrots, 4))
    !
    strf_rot = 0.d0
    !
    ! strf_full = strf_full - strf_q0
    IF (ionode) CALL rotate(DBLE(strf_full), q, nq, nq_super, nrots, ctr2, &
                            strf_rot, col1, col2)
    CALL mp_bcast(strf_rot, ionode_id, world_comm)
    CALL disca_broadening(strf_rot, ctr2 * nrots, kres1, kres2, alat, &
                          kmin, kmax, col1, col2, Np, 'strf_all-ph_broad.dat')
    ! 
    IF (atom_resolved) THEN
      DO k = 1, nat
        DO kp = 1, k
        strf_rot = 0.d0
        WRITE( pt_mz,'(i1.1)') k
        WRITE( pt_mz2,'(i1.1)') kp
        filename = 'strf_all_k_'//TRIM(pt_mz)//'_'//TRIM(pt_mz2)//'_broad.dat'
        ! 
        IF (ionode) CALL rotate(DBLE(strf_full_kk(:, k, kp)), &
                                q, nq, nq_super, nrots, ctr2, strf_rot, col1, col2)
        CALL mp_bcast(strf_rot, ionode_id, world_comm)
        CALL disca_broadening(strf_rot, ctr2 * nrots, kres1, kres2, alat, &
                              kmin, kmax, col1, col2, Np, filename)
                             !'zero_one_k_'//TRIM()//'_'//TRIM(pt_mz2)//'_broad.dat')
        ENDDO ! kp
      ENDDO
    ENDIF
    !
    DEALLOCATE(strf_rot)
    !
    !
    IF (print_raw_data) THEN
      IF (ionode) THEN 
        filename = 'strf_all-phonon.dat'
        OPEN (unit = 66, file = filename, status = 'unknown', form = 'formatted')
        DO qp = nq_super + 1, nq 
           WRITE (66, '(40f26.6)') q(:,qp) * units, DBLE(strf_full(qp)) !, aimag(strf_full(qp))
        ENDDO
        CLOSE(66)
        IF (atom_resolved) THEN
          DO k = 1, nat
            WRITE( pt_mz,'(i1.1)') k
            filename = 'strf_q_full_k_resolved_' // TRIM( pt_mz ) // '.dat' !'.fp'
            OPEN (unit = 66, file = filename, status = 'unknown', form = 'formatted')
            DO qp = nq_super + 1, nq 
             WRITE (66, '(40f26.6)') q(:, qp) * units, & 
                            ( DBLE((strf_full_kk(qp, k, kp) + strf_full_kk(qp, kp, k)) / 2.0d0), kp = 1, k)  
             ! We take the real part iff we sum [strf_full_kk(qp, 2, 1) + strf_full_kk(qp, 1, 2)]/2.
             ! strf_full_kk(qp, 2, 1) is not necessarilly real
            ENDDO ! qp
            CLOSE(66)
          ENDDO
        ENDIF ! atom_resolved
      ENDIF
      ENDIF ! print_raw_data
      !
      DEALLOCATE(Rlist, strf_full, Q_ppkkaa)
    IF (atom_resolved) DEALLOCATE(strf_full_kk)
  !
  ENDIF ! full_phonon
  !
  DEALLOCATE(DW_T_fact_kkp)
  RETURN ! if we DO not want to go further
    !
    !
END SUBROUTINE diffuse_scattering
!
!SUBROUTINE  qpoint_gen1_serial(dim1, dim2, dim3, ctrAB) 
!!
!  use kinds, only: dp
!  
!  IMPLICIT NONE
!  ! input
!  INTEGER, INTENT(in)             :: dim1, dim2, dim3
!  INTEGER, INTENT(out)            :: ctrAB
!  ! local
!  INTEGER                         :: i, j, k, n, ctr, nqs
!  REAL(DP), ALLOCATABLE           :: q_all(:,:)
!  REAL(DP)                        :: q_B(3), q_A(3), eps
!  !
!  nqs = dim1 * dim2 * dim3  
!  eps = 1.0E-06
!  !
!  ALLOCATE(q_all(3, nqs))
!  !
!  DO i = 1, dim1
!      DO j = 1, dim2
!         DO k = 1, dim3
!            !  this is nothing but consecutive ordering
!            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
!            !  q_all are the components of the complete grid in crystal axis
!            q_all(1, n) = DBLE(i - 1) / dim1 ! + DBLE(k1)/2/dim1
!            q_all(2, n) = DBLE(j - 1) / dim2 ! + DBLE(k2)/2/dim2
!            q_all(3, n) = DBLE(k - 1) / dim3 ! + DBLE(k3)/2/dim3 ! k1 , k2 , k3 is for the shift
!         ENDDO
!      ENDDO
!  ENDDO
!  !
!  ctr = 0
!  ctrAB = 0
!  DO i = 1, nqs
!    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
!    IF (((ABS(q_A(1)) .LT. eps) .OR. (ABS(abs(q_A(1)) - 1) .LT. eps)) .AND. &
!        ((ABS(q_A(2)) .LT. eps) .OR. (ABS(abs(q_A(2)) - 1) .LT. eps)) .AND. &
!        ((ABS(q_A(3)) .LT. eps) .OR. (ABS(abs(q_A(3)) - 1) .LT. eps))) THEN
!        ctrAB = ctrAB + 1
!    ELSE
!     DO j = i + 1, nqs
!        q_B = q_all(:, i) + q_all(:, j)       
!       IF (((ABS(q_B(1)) .LT. eps) .OR. (ABS(abs(q_B(1)) - 1) .LT. eps)) .AND. &
!           ((ABS(q_B(2)) .LT. eps) .OR. (ABS(abs(q_B(2)) - 1) .LT. eps)) .AND. &
!           ((ABS(q_B(3)) .LT. eps) .OR. (ABS(abs(q_B(3)) - 1) .LT. eps))) THEN
!           ctr = ctr + 1
!           ctrAB = ctrAB + 1
!       ENDIF 
!      ENDDO
!    ENDIF
!  ENDDO
!  !
!  DEALLOCATE(q_all)
!  ! 
!  !
!END SUBROUTINE qpoint_gen1_serial
!
!
!SUBROUTINE  qpoint_gen2_serial(dim1, dim2, dim3, ctrAB, q_AB) 
!!
!  use kinds, only: dp
!  
!  IMPLICIT NONE
!  ! input
!  INTEGER, INTENT(in)             :: dim1, dim2, dim3, ctrAB
!  REAL(DP), INTENT(out)           :: q_AB(3, ctrAB)
!  ! local
!  INTEGER                         :: i, j, k, n, ctr, nqs
!  REAL(DP), ALLOCATABLE           :: q_all(:, :)
!  REAL(DP)                        :: q_B(3), q_A(3), eps
!  !
!  nqs = dim1 * dim2 * dim3  
!  eps = 1.0E-06
!  !
!  ALLOCATE(q_all(3, nqs))
!  DO i = 1, dim1
!      DO j = 1, dim2
!         DO k = 1, dim3
!            !  this is nothing but consecutive ordering
!            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
!            !  q_all are the components of the complete grid in crystal axis
!            q_all(1, n) = DBLE(i - 1) / dim1 ! + DBLE(k1)/2/dim1
!            q_all(2, n) = DBLE(j - 1) / dim2 ! + DBLE(k2)/2/dim2
!            q_all(3, n) = DBLE(k - 1) / dim3 ! + DBLE(k3)/2/dim3 ! k1 , k2 , k3 is for the shift
!         ENDDO
!      ENDDO
!  ENDDO
!  !
!  ctr = 0
!  DO i = 1, nqs
!    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
!    IF (((ABS(q_A(1)) .LT. eps) .OR. (ABS(abs(q_A(1)) - 1) .LT. eps)) .AND. &
!        ((ABS(q_A(2)) .LT. eps) .OR. (ABS(abs(q_A(2)) - 1) .LT. eps)) .AND. &
!        ((ABS(q_A(3)) .LT. eps) .OR. (ABS(abs(q_A(3)) - 1) .LT. eps))) THEN
!        ctr = ctr + 1
!        q_AB(:, ctr) = q_all(:, i)
!        !WRITE(*,*) "A", q_AB(:, ctr)
!!        WRITE(*,'(A,3F10.6)') "A", q_all(:,i)
!    ELSE
!     DO j = i + 1, nqs
!        q_B = q_all(:, i) + q_all(:, j)       
!       IF (((ABS(q_B(1)) .LT. eps) .OR. (ABS(abs(q_B(1)) - 1) .LT. eps)) .AND. &
!           ((ABS(q_B(2)) .LT. eps) .OR. (ABS(abs(q_B(2)) - 1) .LT. eps)) .AND. &
!           ((ABS(q_B(3)) .LT. eps) .OR. (ABS(abs(q_B(3)) - 1) .LT. eps))) THEN
!           ctr = ctr + 1
!           q_AB(:, ctr) = q_all(:, i)
!         !  WRITE(*,*) q_AB(:, ctr)
!!           WRITE(*,'(A,3F10.6)') "B", q_all(:,i)
!       ENDIF 
!      ENDDO
!    ENDIF
!  ENDDO
!  !
!  DEALLOCATE(q_all)
!
!END SUBROUTINE qpoint_gen2_serial

  SUBROUTINE fkbounds( nktot, lower_bnd, upper_bnd )
  !-----------------------------------------------------------------------
  !!
  !!   Subroutine from EPW finds the lower and upper bounds a k-grid in parallel
  !!
  !! @ Note: 
  !!    If you have 19 kpts and 2 pool, this routine will RETURN
  !!    lower_bnd= 1 and upper_bnd= 10 for the first pool
  !!    lower_bnd= 1 and upper_bnd=9 for the second pool
  !-----------------------------------------------------------------------
  !
  USE mp_global,    ONLY: my_pool_id, npool
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: nktot
  !! nktot k-points splited over pools
  INTEGER, INTENT (out) :: lower_bnd
  !! Lower kpt bounds for that image pool 
  INTEGER, INTENT (out) :: upper_bnd
  !! Upper kpt for that image pool
  !
#if defined(__MPI)
  !
  INTEGER :: nkl, nkr
  !
  ! find the bounds of k-dependent arrays in the parallel case
  ! number of kpoint blocks, kpoints per pool and reminder
  !
  nkl = nktot / npool
  nkr = nktot - nkl * npool
  !
  ! the reminder goes to the first nkr pools (0...nkr-1)
  !
  IF (my_pool_id < nkr ) nkl = nkl + 1
  !
  ! the index of the first k point in this pool
  !
  lower_bnd = my_pool_id * nkl + 1
  IF ( my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
  !
  ! the index of the last k point in this pool
  !
  upper_bnd = lower_bnd + nkl - 1
  !
#else  
  !     
  ! In serial the definitions are much easier 
  !     
  lower_bnd = 1
  upper_bnd = nktot
  !     
#endif 
  !
  RETURN
  !
  END SUBROUTINE fkbounds
!

SUBROUTINE qpoint_gen_map_1(nk1, nk2, nk3, k1, k2, k3, kf1, kf2, kf3, & 
                plane_val, plane_dir, nq_strft, qlist_strf, qlist_strf_cart)
!
USE kinds, ONLY: DP
USE cell_base,  ONLY : bg
IMPLICIT NONE
 
 INTEGER, intent(in)     :: nk1, nk2, nk3, k1, k2, k3, kf1, kf2, kf3, plane_dir
 REAL(DP), intent(in)    :: plane_val
 INTEGER, intent(out)    :: nq_strft
 REAL(DP), intent(out)   :: qlist_strf(nk1 * nk2 * nk3 * (kf1 - k1) * (kf2 - k2) * (kf3 - k3), 3)
 REAL(DP), intent(out)   :: qlist_strf_cart(nk1 * nk2 * nk3 * (kf1 - k1) * (kf2 - k2) * (kf3 - k3), 3)
 INTEGER                 :: i, j, k, n, ctr, p, nk_tot
 REAL(DP), ALLOCATABLE   :: xkg(:, :), xkg_or(:, :)
 REAL(DP)                :: q_B(3), q_A(3), eps, B(3, 3)
!
! Reciprocal axis, as printed form espresso but transpose
! 
B = (bg)
!
! 
 eps = 1.0E-05
 !
 !
 nk_tot= nk1 * nk2 * nk3 ! nk1, nk2 ,nk3 define the resolution of the q-grid
 !
 ALLOCATE(xkg(3, nk_tot), xkg_or(3, nk_tot))
 DO i = 1, nk1
     DO j = 1, nk2
        DO k = 1, nk3
           !  this is nothing but consecutive ordering
           n = (k - 1) + (j - 1) * nk3 + (i - 1) * nk2 * nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1, n) = DBLE(i - 1) / nk1 
           xkg(2, n) = DBLE(j - 1) / nk2
           xkg(3, n) = DBLE(k - 1) / nk3 
        ENDDO
     ENDDO
 ENDDO
 xkg_or = xkg 
 !
 ctr=1
 !
 DO p = k1, kf1 - 1
   DO j = k2, kf2 - 1
      DO k = k3, kf3 - 1
         xkg(1, :) =  xkg_or(1, :) + p
         xkg(2, :) =  xkg_or(2, :) + j
         xkg(3, :) =  xkg_or(3, :) + k
        DO i = 1, nk_tot 
           qlist_strf(ctr, :) = xkg(:, i)
           qlist_strf_cart(ctr, :) = B(:, 1) * xkg(1, i) + B(:, 2) * xkg(2, i) + B(:, 3) * xkg(3, i)
           ctr= ctr+1
        END DO
      END DO
   END DO
 END DO
 ! 
 ! To allocate matrix size
 nq_strft = 0
 DO i = 1, nk_tot * (kf1 - k1) * (kf2 - k2) * (kf3 - k3)
   IF ( (qlist_strf_cart(i, plane_dir) .LT. plane_val + eps) .AND. & 
        (qlist_strf_cart(i, plane_dir) .GT. plane_val - eps) ) THEN
      nq_strft = nq_strft + 1
   ENDIF
 ENDDO
 !
 DEALLOCATE(xkg, xkg_or)
 !
END SUBROUTINE qpoint_gen_map_1
!
!
SUBROUTINE qpoint_gen_map_2(nk1, nk2, nk3, k1, k2, k3, kf1, kf2, kf3, &
               plane_val, plane_dir, nq_strft, qlist_strf, qlist_strf_cart, qlist_strf_cryst_z)
!
USE kinds, ONLY: DP
USE cell_base,  ONLY : at
IMPLICIT NONE
 
 INTEGER,  intent(in)    :: nk1, nk2, nk3, k1, k2, k3, nq_strft, plane_dir, kf1, kf2, kf3
 REAL(DP), intent(in)    :: plane_val, qlist_strf(nk1 * nk2 * nk3 * (kf1 - k1) * (kf2 - k2) * (kf3 - k3), 3)
 REAL(DP), intent(in)    :: qlist_strf_cart(nk1 * nk2 * nk3 * (kf1 - k1) * (kf2 - k2) * (kf3 - k3), 3)
 REAL(DP), intent(out)   :: qlist_strf_cryst_z(3, nq_strft)
 INTEGER                 :: i, ctr, nk_tot
 REAL(DP)                :: eps, invB(3, 3)
! 
invB = TRANSPOSE(at)
! Crystal axis, as printed form espresso but NOT transpose
! 
! 
 eps = 1.0E-05
 !
 nk_tot= nk1 * nk2 * nk3 
 ! nk1, nk2 ,nk3 define the resolution of the q-grid
 !
 ctr = 0 
 DO i = 1, nk_tot * (kf1 - k1) * (kf2 - k2) * (kf3 - k3)
    IF ( (qlist_strf_cart(i, plane_dir) .LT. plane_val + eps) .AND. & 
         (qlist_strf_cart(i, plane_dir) .GT. plane_val - eps)  ) THEN
        ctr = ctr + 1
!       qlist_strf_cart_z(:, ctr)  = qlist_strf_cart(i,:)
        qlist_strf_cryst_z(:, ctr) = invB(:, 1) * qlist_strf_cart(i, 1) + & 
                                     invB(:, 2) * qlist_strf_cart(i, 2) + &
                                     invB(:, 3) * qlist_strf_cart(i, 3)
     ENDIF
 ENDDO
 !
 !
 !
END SUBROUTINE qpoint_gen_map_2

SUBROUTINE  qpoint_gen1(dim1, dim2, dim3, ctrAB) 
!
  use kinds, only: dp
  USE mp,   ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm
  
  IMPLICIT NONE
  ! input
  INTEGER, intent(in)             :: dim1, dim2, dim3
  INTEGER, intent(out)            :: ctrAB
  ! local
  INTEGER                         :: i, j, k, n, nqs
  INTEGER                         :: lower_bnd, upper_bnd
  REAL(DP), ALLOCATABLE           :: q_all(:, :)
  REAL(DP)                        :: q_B(3), q_A(3), eps
  !
  nqs = dim1 * dim2 * dim3  
  eps = 1.0E-06
  !
  ALLOCATE(q_all(3, nqs))
  q_all = 0.d0
  !
  CALL fkbounds( dim1, lower_bnd, upper_bnd )
!  DO i = 1, dim1
  DO i = lower_bnd, upper_bnd
      DO j = 1, dim2
         DO k = 1, dim3
            !  this is nothing but consecutive ordering
            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
            !  q_all are the components of the complete grid in crystal axis
            q_all(1, n) = DBLE(i - 1) / dim1 ! + DBLE(k1)/2/dim1
            q_all(2, n) = DBLE(j - 1) / dim2 ! + DBLE(k2)/2/dim2
            q_all(3, n) = DBLE(k - 1) / dim3 ! + DBLE(k3)/2/dim3 ! k1 , k2 , k3 is for the shift
         ENDDO
      ENDDO
  ENDDO
  CALL mp_sum(q_all, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  ctrAB = 0
  CALL fkbounds( nqs, lower_bnd, upper_bnd )
     !
  DO i = lower_bnd, upper_bnd !1, nqs
    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
    IF (((ABS(q_A(1)) .LT. eps) .OR. (ABS(ABS(q_A(1)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(2)) .LT. eps) .OR. (ABS(ABS(q_A(2)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(3)) .LT. eps) .OR. (ABS(ABS(q_A(3)) - 1) .LT. eps))) THEN
        ctrAB = ctrAB + 1
    ELSE
     DO j = i + 1, nqs
        q_B = q_all(:, i) + q_all(:, j)       
       IF (((ABS(q_B(1)) .LT. eps) .OR. (ABS(ABS(q_B(1)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(2)) .LT. eps) .OR. (ABS(ABS(q_B(2)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(3)) .LT. eps) .OR. (ABS(ABS(q_B(3)) - 1) .LT. eps))) THEN
           ctrAB = ctrAB + 1
       END IF 
      END DO
    END IF
  END DO
  CALL mp_sum(ctrAB, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  DEALLOCATE(q_all)
  ! 
  !
END SUBROUTINE qpoint_gen1

SUBROUTINE  qpoint_gen2(dim1, dim2, dim3, ctrAB, q_AB) 
!
  use kinds,      ONLY : dp
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm
  
  IMPLICIT NONE
  ! input
  INTEGER, intent(in)             :: dim1, dim2, dim3, ctrAB
  REAL(DP), intent(out)           :: q_AB(3, ctrAB)
  ! local
  INTEGER                         :: i, j, k, n, ctr, nqs
  INTEGER                         :: lower_bnd, upper_bnd
  REAL(DP), ALLOCATABLE           :: q_all(:, :), q_AB_TMP(:, :)
  REAL(DP)                        :: q_B(3), q_A(3), eps
  !
  nqs = dim1 * dim2 * dim3  
  eps = 1.0E-06
  !
  ALLOCATE(q_all(3, nqs), q_AB_TMP(3, nqs))
  ! 
  q_all = 0.d0
  ! 
  CALL fkbounds( dim1, lower_bnd, upper_bnd )
  DO i = lower_bnd, upper_bnd
  !DO i = 1, dim1
      DO j = 1, dim2
         DO k = 1, dim3
            !  this is nothing but consecutive ordering
            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
            !  q_all are the components of the complete grid in crystal axis
            q_all(1, n) = DBLE(i - 1) / dim1 
            q_all(2, n) = DBLE(j - 1) / dim2 
            q_all(3, n) = DBLE(k - 1) / dim3 
         ENDDO
      ENDDO
  ENDDO
  CALL mp_sum(q_all, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  CALL fkbounds( nqs, lower_bnd, upper_bnd )
     !
  q_AB_TMP = 0.d0
  DO i = lower_bnd, upper_bnd !1, nqs
    !DO i = 1, nqs
    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
    IF (((ABS(q_A(1)) .LT. eps) .OR. (ABS(ABS(q_A(1)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(2)) .LT. eps) .OR. (ABS(ABS(q_A(2)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(3)) .LT. eps) .OR. (ABS(ABS(q_A(3)) - 1) .LT. eps))) THEN
        q_AB_TMP(:, i) = q_all(:, i)
    ELSE
     DO j = i + 1, nqs
        q_B = q_all(:, i) + q_all(:, j)       
       IF (((ABS(q_B(1)) .LT. eps) .OR. (ABS(ABS(q_B(1)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(2)) .LT. eps) .OR. (ABS(ABS(q_B(2)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(3)) .LT. eps) .OR. (ABS(ABS(q_B(3)) - 1) .LT. eps))) THEN
           q_AB_TMP(:, i) = q_all(:, i)
       END IF 
      END DO
    END IF
  ! 
  END DO
  CALL mp_sum(q_AB_TMP, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  ctr = 1 ! so that Gamma is the first entry
  q_AB = 0.d0
  DO i = 1, nqs
    IF ((SUM(ABS(q_AB_TMP(:, i)))) .GT. eps ) THEN
        ctr = ctr + 1
        q_AB(:, ctr) = q_AB_TMP(:, i)
    ENDIF
  ENDDO
  !
  DEALLOCATE(q_all, q_AB_TMP)
!
END SUBROUTINE qpoint_gen2

!
!-----------------------------------------------------------------------
SUBROUTINE frc_blk(dyn,q,tau,nat,nr1,nr2,nr3,frc,at,bg,rws,nrws,f_of_q,fd)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nat, n1, n2, n3, nr1_, nr2_, nr3_, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws, nax
  COMPLEX(DP) dyn(3,3,nat,nat), f_of_q(3,3,nat,nat)
  REAL(DP) frc(nr1,nr2,nr3,3,3,nat,nat), tau(3,nat), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws), alat
  REAL(DP), EXTERNAL :: wsweight
  REAL(DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: ttt(:,:,:,:,:), tttx(:,:)
  LOGICAL,SAVE :: first=.true.
  LOGICAL :: fd
  !
  nr1_=2*nr1
  nr2_=2*nr2
  nr3_=2*nr3
  FIRST_TIME : IF (first) THEN
    first=.false.
    ALLOCATE( wscache(-nr3_:nr3_, -nr2_:nr2_, -nr1_:nr1_, nat,nat) )
    DO na=1, nat
       DO nb=1, nat
          total_weight=0.0d0
          !
          ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY VERY SAFE RANGE!
          !
          DO n1=-nr1_,nr1_
             DO n2=-nr2_,nr2_
                DO n3=-nr3_,nr3_
                   DO i=1, 3
                      r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                      if (fd) r_ws(i) = r(i) + tau(i,nb)-tau(i,na)
                   END DO
                   wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
                   total_weight=total_weight + wscache(n3,n2,n1,nb,na) 
                ENDDO
             ENDDO
          ENDDO
          IF (ABS(total_weight-nr1*nr2*nr3).GT.1.0d-8) THEN
             WRITE(stdout,*) na,nb,total_weight
             CALL errore ('frc_blk','wrong total_weight',1)
          END IF
      ENDDO
    ENDDO
  ENDIF FIRST_TIME
  !
  ALLOCATE(ttt(3,nat,nr1,nr2,nr3))
  ALLOCATE(tttx(3,nat*nr1*nr2*nr3))
  ttt(:,:,:,:,:)=0.d0

  DO na=1, nat
     DO nb=1, nat
        DO n1=-nr1_,nr1_
           DO n2=-nr2_,nr2_
              DO n3=-nr3_,nr3_
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                 END DO

                 weight = wscache(n3,n2,n1,nb,na) 
                 IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nr1)
                    IF(m1.LE.0) m1=m1+nr1
                    m2 = MOD(n2+1,nr2)
                    IF(m2.LE.0) m2=m2+nr2
                    m3 = MOD(n3+1,nr3)
                    IF(m3.LE.0) m3=m3+nr3
                 !   write(*,'(6i4)') n1,n2,n3,m1,m2,m3
                    !
                    ! FOURIER TRANSFORM
                    !
                    DO i=1,3
                       ttt(i,na,m1,m2,m3)=tau(i,na)+m1*at(i,1)+m2*at(i,2)+m3*at(i,3)
                       ttt(i,nb,m1,m2,m3)=tau(i,nb)+m1*at(i,1)+m2*at(i,2)+m3*at(i,3)
                    END DO

                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) = dyn(ipol,jpol,na,nb) +                &
                               (frc(m1,m2,m3,ipol,jpol,na,nb)+f_of_q(ipol,jpol,na,nb)) &
                               *CMPLX(COS(arg),-SIN(arg),kind=DP)*weight
                       END DO
                    END DO

                 END IF
              END DO
           END DO
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE frc_blk
!
SUBROUTINE rotate(strf, q, nq, nq_super, nrots, & 
                   ctr, strf_rot, col1, col2)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: nq, nrots, col1, col2, ctr, nq_super
  REAL(DP), INTENT(in)   :: strf(nq), q(3, nq)
  REAL(DP), INTENT(out)  :: strf_rot(ctr * nrots, 4)
  INTEGER                :: i, p, ctr2
  REAL(DP)               :: str_f(ctr,4), str_fp(ctr, 4)
  !
  REAL(DP)               :: Rmat(2,2), theta, eps
  ! 
  eps        = 1d-5
  !
  !WRITE(*,*) "Number of rotations provided:", nrots
  theta      = tpi / FLOAT(nrots)
  !
  !
  str_f = 0.d0
  ctr2 = 0
  DO p = nq_super + 1, nq
    IF ((q(col1, p) .GT. 0.0 - eps) .AND. (q(col2, p) .GT. 0.0 - eps)) THEN
      ctr2 = ctr2 + 1
      str_f(ctr2, 1:3) = q(:, p)
      str_f(ctr2, 4) = strf(p)
    ENDIF
  ENDDO
  !
  ! To remove double contribution upon rotation
  !
  str_fp = 0.d0
  str_fp(1, :) = str_f(1, :)
  DO p = 2, ctr
    IF (ATAN(str_f(p, 1) / str_f(p, 2)) .LT. ( tpi / float(nrots) - eps) ) THEN
    str_fp(p, :) = str_f(p, :)
    ENDIF
  ENDDO
  !
  !
  strf_rot = 0.d0
  ctr2 = 1
  DO i = 0, nrots - 1
    Rmat(1, :) = (/ COS(i * theta), -SIN(i * theta) /)
    Rmat(2, :) = (/ SIN(i * theta),  COS(i * theta) /)
    DO p = 1, ctr
      strf_rot(ctr2, 1 : 2) = MATMUL(Rmat, str_fp(p, 1 : 2))
      strf_rot(ctr2, 3) = str_fp(p, 3)
      strf_rot(ctr2, 4) = str_fp(p, 4)
      ctr2 = ctr2 + 1
    END DO
  ENDDO
END SUBROUTINE
!
SUBROUTINE disca_broadening(strf_rot, steps, kres1, kres2, alat, &
                            kmin, kmax, col1, col2, Np, flstrfout)
!-------------------------------------------------------------------------
!! authors: Marios Zacharias, Feliciano Giustino 
!! acknowledgement: Hyungjun Lee for help packaging this release
!! version: v0.1
!! license: GNU
!
USE kinds,       ONLY : dp
USE mp_global,   ONLY : inter_pool_comm
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_bcast, mp_barrier, mp_sum
USE io_global,   ONLY : ionode, ionode_id, stdout
USE constants,   ONLY : pi, BOHR_RADIUS_ANGS
!
IMPLICIT NONE
!
 CHARACTER(LEN=256), INTENT(IN)  :: flstrfout
 INTEGER, INTENT(IN)             :: steps, kres1, kres2, Np, col1, col2
 REAL(DP), INTENT(IN)            :: kmin, kmax, alat
 REAL(DP), INTENT(IN)            :: strf_rot(steps, 4)
 INTEGER                         :: ii, lower_bnd, upper_bnd, ik, iky
 REAL(DP)                        :: jump, sf_smearingx, sf_smearingy, maxv, units !, pi 
 REAL(DP), ALLOCATABLE           :: kgridx(:), kgridy(:), strf_out(:, :)
!
!
!
units = DBLE(2.d0*pi/alat/BOHR_RADIUS_ANGS)
ALLOCATE( kgridx(kres1), kgridy(kres2))
ALLOCATE(strf_out(kres1,kres2))
!kmin = -10.0
!kmax = 10.0

jump = (kmax - kmin) / DBLE(kres1 - 1)
DO ik = 1, kres1
  kgridx(ik) = kmin + (ik - 1) * jump
ENDDO
sf_smearingx = (kmax - kmin) / DBLE(kres1)
!
jump = (kmax - kmin) / DBLE(kres2 - 1)
DO ik = 1, kres2
  kgridy(ik) = kmin + (ik- 1)*jump
ENDDO
sf_smearingy = (kmax - kmin) / DBLE(kres2)
!! 
!
strf_out = 0.d0
!
CALL fkbounds( steps, lower_bnd, upper_bnd )
!
DO ii = lower_bnd, upper_bnd
  DO ik = 1, kres1 !
    DO iky = 1, kres2
    !
    strf_out(ik, iky) =  strf_out(ik, iky) +  &
                          strf_rot(ii, 4) / sf_smearingx / SQRT(2.0d0 * pi) / sf_smearingy / SQRT(2.0d0 * pi)* &
                          (EXP(-(strf_rot(ii, col1) * units - kgridx(ik))**2.d0 / sf_smearingx**2.d0 / 2.d0))*&
                          (EXP(-(strf_rot(ii, col2) * units - kgridy(iky))**2.d0 / sf_smearingy**2.d0 / 2.d0))
    !
    ENDDO
  ENDDO
ENDDO
!
CALL mp_sum(strf_out, inter_pool_comm)
CALL mp_barrier(inter_pool_comm)
!
IF (ionode) THEN
  OPEN(46,FILE=flstrfout)
  maxv =  maxval(strf_out)
  WRITE(46,*) "#", maxv, maxval(strf_out)
  DO ik = 1, kres1
    DO iky = 1, kres2
      WRITE(46,'(3F28.12)') kgridx(ik), kgridy(iky), strf_out(ik,iky) * Np**(-2.0d0) ! 
                                         !Np**(-2.0d0) ! / maxv
    ENDDO
    WRITE(46,*)
  ENDDO
  CLOSE(46)
ENDIF
!
DEALLOCATE(strf_out, kgridx, kgridy)
!
!
END SUBROUTINE
